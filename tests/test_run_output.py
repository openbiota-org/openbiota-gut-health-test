"""What `openbiota run` puts on the terminal: little, and the right things.

A run used to print every step with a timestamp, then a table of every
stage, then several thousand lines of text summary, with the outcome
somewhere above all of that. The default is now the short form: one line
per stage, a failure once and in words a person can act on, an outcome in a
line or two, and the report's path last. Everything else is still written -
to ``run.log`` and ``summary.txt`` - and ``--verbose`` / ``--summary`` put
it back on the screen.
"""

from __future__ import annotations

import argparse
import io
import os
import re
from pathlib import Path

import pytest

from openbiota import cli, logging_util
from openbiota.completeness import Completeness
from openbiota.logging_util import Reporter, Style, plain_reason

REPO = Path(__file__).resolve().parents[1]
HMMSEARCH = ("singlem pipe failed (exit 1): hmm - returned non-zero exit status 127. STDERR was: "
             "b'bash: line 1: hmmsearch: command not found\\n'")


class _Terminal(io.StringIO):
    def isatty(self) -> bool:
        return True


def _short(monkeypatch, *, terminal: bool = False) -> tuple[Reporter, io.StringIO]:
    stream = _Terminal() if terminal else io.StringIO()
    monkeypatch.setattr(logging_util.sys, "stderr", stream)
    monkeypatch.setenv("TERM", "xterm")
    reporter = Reporter(detail=False)
    assert reporter._redraw is terminal  # noqa: SLF001
    return reporter, stream


# --------------------------------------------------------------------------- #
# the short form
# --------------------------------------------------------------------------- #


def test_short_form_is_one_line_per_stage_without_the_detail(monkeypatch):
    reporter, stream = _short(monkeypatch)
    monkeypatch.setattr(logging_util, "STAGE_SHOWN_AFTER_S", 0.0)
    reporter.info("OpenBiota Gut Health Test v0")
    reporter.headline("S: lite report from 2 FASTQ files (1.0GB).")
    with reporter.stage("input validation"):
        reporter.info("  R1: profiling S_1.fastq.gz (9.8MB)")
        reporter.ok("fastp: 120,000 pairs")
        with reporter.stage("a step inside it"):
            reporter.step("deeper still")
    shown = stream.getvalue().splitlines()
    assert shown[0] == "S: lite report from 2 FASTQ files (1.0GB)."
    assert re.fullmatch(r"  ✓ input validation +\d+\.\ds", shown[1])
    assert len(shown) == 2, shown
    assert "[" not in stream.getvalue(), "no timestamps in the short form"


def test_a_stage_too_quick_to_matter_leaves_no_line(monkeypatch):
    reporter, stream = _short(monkeypatch)
    with reporter.stage("tally and normalisation"):
        pass
    assert stream.getvalue() == ""
    assert any("tally and normalisation done in" in line for line in reporter._lines)  # noqa: SLF001


def test_into_a_log_a_long_stage_is_announced_once_then_finished(monkeypatch):
    """Not a terminal, so no redrawn line: the last line of the log still says what is running."""
    reporter, stream = _short(monkeypatch)
    monkeypatch.setattr(logging_util, "STAGE_ANNOUNCED_AFTER_S", 0.0)
    with reporter.stage("competitive confirmation"):
        reporter._refresh()  # noqa: SLF001 - the once-a-second tick
        reporter._refresh()  # noqa: SLF001
        assert stream.getvalue() == "  · competitive confirmation …\n"
    shown = stream.getvalue().splitlines()
    assert len(shown) == 2 and shown[1].startswith("  ✓ competitive confirmation ")
    assert "\r" not in stream.getvalue() and "\033" not in stream.getvalue()


def test_the_log_keeps_every_line_whatever_was_shown(monkeypatch):
    reporter, _ = _short(monkeypatch)
    with reporter.stage("input validation"):
        reporter.info("  R1: profiling")
        with reporter.stage("a step inside it"):
            reporter.ok("deep detail")
    log = "\n".join(reporter._lines)  # noqa: SLF001
    for text in ("· input validation …", "R1: profiling", "· a step inside it …", "✓ deep detail",
                 "✓ a step inside it done in", "✓ input validation done in"):
        assert text in log, text
    assert all(line.startswith("[") for line in reporter._lines)  # noqa: SLF001


def test_the_log_is_plain_text_even_when_the_terminal_is_coloured(monkeypatch):
    monkeypatch.setattr(logging_util.sys, "stderr", io.StringIO())
    reporter = Reporter(_style=Style(True))
    reporter.ok("coloured on screen")
    reporter.warn("and this")
    assert "\033" not in "\n".join(reporter._lines)  # noqa: SLF001


def test_on_a_terminal_the_running_stage_is_one_line_redrawn_in_place(monkeypatch):
    reporter, stream = _short(monkeypatch, terminal=True)
    monkeypatch.setattr(logging_util, "STAGE_SHOWN_AFTER_S", 0.0)
    with reporter.stage("competitive confirmation"):
        reporter.info("  mapping reads to 212 genomes")
        reporter._refresh()  # noqa: SLF001 - the once-a-second tick
    text = stream.getvalue()
    assert "\r\033[2K  · competitive confirmation  0s" in text
    assert "mapping reads to 212 genomes" in text, "the latest detail rides beside the stage name"
    permanent = [line for line in text.split("\n") if "✓" in line]
    assert len(permanent) == 1 and permanent[0].endswith("s")
    assert text.count("\n") == 1, "only the finished stage takes a line"


def test_a_step_inside_a_stage_shows_beside_it_not_under_it(monkeypatch):
    reporter, stream = _short(monkeypatch, terminal=True)
    with reporter.stage("taxonomic engine (MetaPhlAn)"), reporter.stage("host filtering"):
        reporter._draw_live()  # noqa: SLF001
    assert "· taxonomic engine (MetaPhlAn)  0s  host filtering" in re.sub(r"\033\[[0-9;]*m", "", stream.getvalue())
    assert "\n" not in stream.getvalue()


def test_on_a_terminal_an_instant_stage_leaves_no_line(monkeypatch):
    reporter, stream = _short(monkeypatch, terminal=True)
    with reporter.stage("tally and normalisation"):
        pass
    assert "\n" not in stream.getvalue()
    assert stream.getvalue().endswith("\r\033[2K"), "and the status line is wiped, not left behind"


def test_work_between_stages_shows_as_working_on_a_terminal(monkeypatch):
    """Half a minute of unstaged work with nothing on the screen reads as a hang."""
    reporter, stream = _short(monkeypatch, terminal=True)
    monkeypatch.setattr(logging_util, "IDLE_SHOWN_AFTER_S", 0.0)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue() == "", "nothing before the first stage"
    with reporter.stage("nucleotide locus mapping"):
        pass
    reporter._refresh()  # noqa: SLF001
    assert re.sub(r"\033\[[0-9;]*m", "", stream.getvalue()).endswith("\r\033[2K  · working  0s"[9:])
    reporter.warn("a permanent line wipes it first")
    assert stream.getvalue().endswith("\r\033[2K  ! a permanent line wipes it first\n")


def test_close_wipes_the_status_line_for_good(monkeypatch):
    reporter, stream = _short(monkeypatch, terminal=True)
    monkeypatch.setattr(logging_util, "IDLE_SHOWN_AFTER_S", 0.0)
    with reporter.stage("PDF report"):
        pass
    reporter._refresh()  # noqa: SLF001
    reporter.close()
    assert stream.getvalue().endswith("\r\033[2K")
    before = stream.getvalue()
    reporter._refresh()  # noqa: SLF001
    reporter.headline("✓ Done in 1m00s")
    assert stream.getvalue() == before + "✓ Done in 1m00s\n", "nothing is drawn after the run's last word"


def test_the_refresh_thread_ends_with_its_reporter(monkeypatch):
    import gc
    import threading
    import time

    reporter, _ = _short(monkeypatch, terminal=True)
    real_sleep = time.sleep
    monkeypatch.setattr(logging_util.time, "sleep", lambda _s: real_sleep(0.01))
    with reporter.stage("a stage"):
        ticker = reporter._ticker  # noqa: SLF001
        assert isinstance(ticker, threading.Thread) and ticker.is_alive()
    del reporter
    gc.collect()
    ticker.join(timeout=5)
    assert not ticker.is_alive()


def test_the_status_line_never_wraps(monkeypatch):
    reporter, stream = _short(monkeypatch, terminal=True)
    monkeypatch.setattr(logging_util.shutil, "get_terminal_size", lambda _fallback=None: os.terminal_size((40, 24)))
    with reporter.stage("a stage with a very long name indeed, longer than the window is wide"):
        reporter.info("x" * 500)
        reporter._draw_live()  # noqa: SLF001
    for drawn in stream.getvalue().split("\r\033[2K"):
        plain = re.sub(r"\033\[[0-9;]*m", "", drawn.split("\n")[0])
        if "·" in plain:
            assert len(plain) <= 39, (len(plain), plain)


def test_a_failed_stage_says_why_once_in_words(monkeypatch):
    reporter, stream = _short(monkeypatch)
    label = "unrepresented-lineage rescue (SingleM / GlobDB)"
    try:
        with reporter.stage(label):
            raise RuntimeError(HMMSEARCH)
    except RuntimeError as exc:
        # what the lane runner does next: it explains the failure it caught
        reporter.warn(f"{label} failed: RuntimeError: {str(exc)[:200]}")
    shown = stream.getvalue().splitlines()
    assert shown[-1] == f"  ✗ {label}: the program 'hmmsearch' is not installed - run ./install.sh to add it"
    assert sum("hmmsearch" in line for line in shown) == 1, shown
    assert "Traceback" not in stream.getvalue() and "exit status" not in stream.getvalue()


def test_a_warning_is_one_line_and_said_once(monkeypatch):
    reporter, stream = _short(monkeypatch)
    for _ in range(3):
        reporter.warn("host index not built: human reads are not removed\nmore\nlines")
    reporter.warn("x" * 400)
    shown = stream.getvalue().splitlines()
    assert shown[0] == "  ! host index not built: human reads are not removed"
    assert len(shown) == 2
    assert shown[1].endswith("...") and len(shown[1]) < 240


def test_an_unrelated_warning_after_a_failure_is_still_shown(monkeypatch):
    reporter, stream = _short(monkeypatch)
    with pytest.raises(ValueError, match="broken"), reporter.stage("stage one"):
        raise ValueError("broken")
    reporter.warn("something else entirely")
    assert "  ! something else entirely" in stream.getvalue().splitlines()


def test_an_error_shows_its_first_line_and_logs_the_rest(monkeypatch):
    reporter, stream = _short(monkeypatch)
    reporter.error("PDF report failed: ValueError: x\nTraceback (most recent call last):\n  File ...")
    assert stream.getvalue().splitlines() == ["  ✗ PDF report failed: ValueError: x"]
    assert "Traceback" in "\n".join(reporter._lines)  # noqa: SLF001


def test_ctrl_c_in_a_stage_adds_no_failure_line(monkeypatch):
    reporter, stream = _short(monkeypatch)
    with pytest.raises(KeyboardInterrupt), reporter.stage("a long stage"):
        raise KeyboardInterrupt
    assert "✗" not in stream.getvalue()


# --------------------------------------------------------------------------- #
# the other two sizes
# --------------------------------------------------------------------------- #


def test_the_detailed_form_is_every_line_with_its_time(monkeypatch):
    stream = io.StringIO()
    monkeypatch.setattr(logging_util.sys, "stderr", stream)
    reporter = Reporter()
    with reporter.stage("input validation"):
        reporter.info("  R1: profiling")
    shown = stream.getvalue().splitlines()
    assert len(shown) == 3
    assert re.match(r"\[ +\d+\.\ds\] · input validation …", shown[0])
    assert "R1: profiling" in shown[1]
    assert re.search(r"✓ input validation done in \d", shown[2])


def test_quiet_prints_nothing_in_either_form(monkeypatch):
    for detail in (True, False):
        stream = io.StringIO()
        monkeypatch.setattr(logging_util.sys, "stderr", stream)
        reporter = Reporter(verbose=False, detail=detail)
        reporter.headline("h")
        with reporter.stage("s"):
            reporter.info("i")
        reporter.warn("w")
        reporter.error("e")
        assert stream.getvalue() == ""
        assert len(reporter._lines) >= 5  # noqa: SLF001


# --------------------------------------------------------------------------- #
# reasons a person can act on
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("detail", "expected"), [
    (f"failed: RuntimeError: {HMMSEARCH}", "the program 'hmmsearch' is not installed - run ./install.sh to add it"),
    ("/bin/sh: orfm: not found", "the program 'orfm' is not installed - run ./install.sh to add it"),
    ("FileNotFoundError: [Errno 2] No such file or directory: 'bracken'",
     "missing file or program: bracken - run ./install.sh"),
    ("ModuleNotFoundError: No module named 'singlem'", "a Python package is missing (singlem) - run ./install.sh"),
    ("failed: RuntimeError: kraken2 read no sequences from S_1.fastq", "kraken2 read no sequences from S_1.fastq"),
    ("Traceback (most recent call last):\n  File \"x.py\", line 1\nValueError: the real message", "the real message"),
    ("status stage_error", "status stage_error"),
    ("failed: RuntimeError: MetaPhlAn Jan26 failed (exit 1): Killed", logging_util.OUT_OF_MEMORY),
    ("bowtie2-align died with signal 9", logging_util.OUT_OF_MEMORY),
    ("sh: line 1: 4411 Killed: 9   bowtie2 --seed 1992", logging_util.OUT_OF_MEMORY),
    ("Command '['kraken2']' returned non-zero exit status 137.", logging_util.OUT_OF_MEMORY),
    ("RuntimeError: 12 reads were killed by the filter", "12 reads were killed by the filter"),
    ("", "failed without a message; see run.log"),
])
def test_plain_reason(detail: str, expected: str):
    assert plain_reason(detail) == expected


def test_plain_reason_is_one_bounded_line():
    out = plain_reason("RuntimeError: " + "y" * 1000 + "\nsecond line")
    assert "\n" not in out and len(out) <= 163 and out.endswith("...")


# --------------------------------------------------------------------------- #
# the run's first and last words
# --------------------------------------------------------------------------- #


def _stages(*missing: tuple[str, str], n: int = 19) -> Completeness:
    rows = [{"key": f"s{i}", "label": f"stage {i}", "ran": True, "detail": "ok"} for i in range(n - len(missing))]
    rows += [{"key": label, "label": label, "ran": False, "detail": detail} for label, detail in missing]
    return Completeness(rows, mode="lite")


def test_a_complete_run_ends_in_one_line(monkeypatch):
    reporter, stream = _short(monkeypatch)
    cli._short_outcome(reporter, completeness=_stages(), total=3057.0,  # noqa: SLF001
                       results={"organism_inventory": {"n_organisms": 332}})
    assert stream.getvalue().splitlines() == [
        "", "✓ Done in 50m57s: lite report, all 19 stages ran, 332 organisms named."]


def test_an_incomplete_run_names_only_what_is_missing_and_why(monkeypatch):
    reporter, stream = _short(monkeypatch)
    label = "unrepresented-lineage rescue (SingleM / GlobDB)"
    cli._short_outcome(reporter, completeness=_stages((label, f"failed: RuntimeError: {HMMSEARCH}")),  # noqa: SLF001
                       total=3057.0, results={})
    shown = stream.getvalue().splitlines()
    assert shown[1] == "✗ Finished in 50m57s, but 1 of 19 stages produced no result:"
    assert shown[2] == f"    {label}: the program 'hmmsearch' is not installed - run ./install.sh to add it"
    assert shown[3].startswith("  The report was written all the same")
    assert len(shown) == 4, "not the table of every stage"
    assert "full install" not in stream.getvalue(), "a failure is not something --full adds"


def test_the_first_line_says_what_is_being_made_and_how_long_it_takes():
    class Mate:
        pass

    class Sample:
        sample = "EM1"
        mates = (Mate(), Mate())

    args = argparse.Namespace(report_mode="lite", subsample=None)
    line = cli._run_headline(args, Sample(), 10 * 1024**3)  # noqa: SLF001
    assert line == ("EM1: lite report from 2 FASTQ files (10.0GB). Usually one to two hours for a new "
                    "sample, about five minutes for one run before.")
    args.mode_reason = "chosen at the prompt"
    assert cli._run_headline(args, Sample(), 10 * 1024**3) == "EM1: lite report from 2 FASTQ files (10.0GB).", \
        "the prompt has just given the times"  # noqa: SLF001
    args.report_mode = None
    assert cli._run_headline(args, Sample(), 1024) == "EM1: report from 2 FASTQ files (1.0KB)."  # noqa: SLF001


def test_the_long_summary_and_the_step_by_step_are_opt_in():
    parser = cli.build_parser()
    plain = parser.parse_args(["run"])
    assert plain.summary is False and plain.verbose is False
    asked = parser.parse_args(["run", "--summary", "-v"])
    assert asked.summary is True and asked.verbose is True


def test_run_prints_the_text_summary_only_on_request():
    """The stdout block of cmd_run: --json, --compact, --summary, else nothing."""
    source = (REPO / "openbiota" / "cli.py").read_text()
    block = source[source.index("    # ---- stdout"):source.index("    if pdf_error is not None:\n        return 1")]
    assert 'elif getattr(args, "summary", False):' in block
    assert block.index('getattr(args, "summary", False)') < block.index("render_summary(")
    assert "\n    else:\n" not in block, "no branch prints unasked"
    assert "Reporter(verbose=not args.quiet, detail=bool(getattr(args, \"verbose\", False)))" in source
