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
import contextlib
import gc
import io
import re
import threading
import time
import weakref
from collections.abc import Iterator
from pathlib import Path

import pytest

from openbiota import cli, logging_util
from openbiota.completeness import Completeness
from openbiota.logging_util import Reporter, Style, plain_reason

REPO = Path(__file__).resolve().parents[1]
HMMSEARCH = ("singlem pipe failed (exit 1): hmm - returned non-zero exit status 127. STDERR was: "
             "b'bash: line 1: hmmsearch: command not found\\n'")


#: Colour and weight. Not the line-clearing sequence, which is part of what is being tested.
SGR = re.compile(r"\033\[[0-9;]*m")
#: Back to the start of the line and wipe it: how the running stage's line is redrawn.
CLEAR = "\r\033[2K"


class _Terminal(io.StringIO):
    def isatty(self) -> bool:
        return True


class _Clock:
    """The reporter's clock, moved by hand: nothing here depends on how fast the machine is."""

    def __init__(self) -> None:
        self.now = 1000.0

    def monotonic(self) -> float:
        return self.now

    def sleep(self, _seconds: float) -> None:
        time.sleep(0.01)   # the refresh thread's one-second wait, when a test lets that thread run


class _Short:
    """Reporters in the short form, in an environment the test states rather than inherits.

    Whether the program colours its output, redraws a line in place, and how
    wide it takes the window to be are read from the environment (NO_COLOR,
    TERM, COLUMNS) and from the stream. A developer's shell, an editor's
    agent and a CI runner set those differently - two of these tests passed
    where NO_COLOR happened to be set and failed on GitHub, where it is not -
    so every reporter made here is told all of them. Time is the test's too:
    the clock moves only when `advance` is called, and no refresh thread
    runs unless the test asks for it, so "0s" cannot become "1s" on a busy
    machine.
    """

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self._monkeypatch = monkeypatch
        self.clock = _Clock()
        self._made: list[weakref.ReferenceType[Reporter]] = []
        monkeypatch.setattr(logging_util, "time", self.clock)

    def __call__(self, *, terminal: bool = False, colour: bool = False, columns: int = 100,
                 ticking: bool = False) -> tuple[Reporter, io.StringIO]:
        patch = self._monkeypatch
        stream = _Terminal() if terminal else io.StringIO()
        patch.setattr(logging_util.sys, "stderr", stream)
        patch.setenv("TERM", "xterm-256color")
        patch.setenv("COLUMNS", str(columns))
        if colour:
            patch.delenv("NO_COLOR", raising=False)
        else:
            patch.setenv("NO_COLOR", "1")
        # Unless the test is about the refresh thread, its place is taken and
        # none starts: the test calls _refresh() where a second would pass.
        reporter = Reporter(detail=False, _t0=self.clock.now,
                            _ticker=None if ticking else threading.current_thread())
        assert reporter._redraw is terminal  # noqa: SLF001
        assert reporter.style.enabled is (terminal and colour)
        self._made.append(weakref.ref(reporter))
        return reporter, stream

    def advance(self, seconds: float) -> None:
        self.clock.now += seconds

    def close(self) -> None:
        for ref in self._made:
            reporter = ref()
            if reporter is not None:
                reporter.close()


@pytest.fixture
def short(monkeypatch: pytest.MonkeyPatch) -> Iterator[_Short]:
    session = _Short(monkeypatch)
    yield session
    session.close()


def _done(name: str, took: str) -> str:
    """A finished stage's line: the name, then the duration in a column of its own."""
    return f"  ✓ {name:<58} {took:>8}"


# --------------------------------------------------------------------------- #
# the short form
# --------------------------------------------------------------------------- #


def test_short_form_is_one_line_per_stage_without_the_detail(short):
    reporter, stream = short()
    reporter.info("OpenBiota Gut Health Test v0")
    reporter.headline("S: lite report from 2 FASTQ files (1.0GB).")
    with reporter.stage("input validation"):
        reporter.info("  R1: profiling S_1.fastq.gz (9.8MB)")
        reporter.ok("fastp: 120,000 pairs")
        with reporter.stage("a step inside it"):
            reporter.step("deeper still")
        short.advance(72)
    assert stream.getvalue().splitlines() == [
        "S: lite report from 2 FASTQ files (1.0GB).",
        "  ✓ input validation" + " " * 46 + "1m12s",
    ]
    assert _done("input validation", "1m12s") == stream.getvalue().splitlines()[1]
    assert "[" not in stream.getvalue(), "no timestamps in the short form"


def test_a_stage_too_quick_to_matter_leaves_no_line(short):
    reporter, stream = short()
    with reporter.stage("tally and normalisation"):
        short.advance(0.9)
    assert stream.getvalue() == ""
    assert any("tally and normalisation done in 0.9s" in line for line in reporter._lines)  # noqa: SLF001
    with reporter.stage("per-panel hit views"):
        short.advance(1.0)
    assert stream.getvalue() == _done("per-panel hit views", "1.0s") + "\n", "a second is long enough to say"


def test_into_a_log_a_long_stage_is_announced_once_then_finished(short):
    """Not a terminal, so no redrawn line: the last line of the log still says what is running."""
    reporter, stream = short()
    with reporter.stage("competitive confirmation"):
        short.advance(1)
        reporter._refresh()  # noqa: SLF001 - what the refresh thread does once a second
        assert stream.getvalue() == "", "not yet: most stages are over in a moment"
        short.advance(1)
        reporter._refresh()  # noqa: SLF001
        assert stream.getvalue() == "  · competitive confirmation …\n"
        short.advance(60)
        reporter._refresh()  # noqa: SLF001
        assert stream.getvalue() == "  · competitive confirmation …\n", "said once, however long it runs"
    assert stream.getvalue().splitlines() == [
        "  · competitive confirmation …", _done("competitive confirmation", "1m02s")]
    assert "\r" not in stream.getvalue() and "\033" not in stream.getvalue()


def test_the_log_keeps_every_line_whatever_was_shown(short):
    reporter, _ = short()
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


def test_on_a_terminal_the_running_stage_is_one_line_redrawn_in_place(short):
    reporter, stream = short(terminal=True)
    with reporter.stage("competitive confirmation"):
        assert stream.getvalue() == CLEAR + "  · competitive confirmation  0s"
        reporter.info("  mapping reads to 212 genomes")
        short.advance(65)
        reporter._refresh()  # noqa: SLF001 - what the refresh thread does once a second
        assert stream.getvalue().endswith(
            CLEAR + "  · competitive confirmation  1m05s  mapping reads to 212 genomes"
        ), "the clock, and the latest detail beside the stage name"
        assert "\n" not in stream.getvalue(), "redrawn in place: nothing has scrolled"
    assert stream.getvalue().endswith(CLEAR + _done("competitive confirmation", "1m05s") + "\n")
    assert stream.getvalue().count("\n") == 1, "only the finished stage takes a line"


def test_a_step_inside_a_stage_shows_beside_it_not_under_it(short):
    reporter, stream = short(terminal=True)
    with reporter.stage("taxonomic engine (MetaPhlAn)"), reporter.stage("host filtering"):
        reporter._refresh()  # noqa: SLF001
    assert CLEAR + "  · taxonomic engine (MetaPhlAn)  0s  host filtering" in stream.getvalue()
    assert "\n" not in stream.getvalue()


def test_on_a_terminal_an_instant_stage_leaves_no_line(short):
    reporter, stream = short(terminal=True)
    with reporter.stage("tally and normalisation"):
        pass
    assert stream.getvalue() == CLEAR + "  · tally and normalisation  0s" + CLEAR, \
        "drawn, then wiped: nothing is left behind and nothing has scrolled"


def test_work_between_stages_shows_as_working_on_a_terminal(short):
    """Half a minute of unstaged work with nothing on the screen reads as a hang."""
    reporter, stream = short(terminal=True)
    short.advance(30)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue() == "", "nothing before the first stage: the run has not begun"
    with reporter.stage("nucleotide locus mapping"):
        pass
    quiet = stream.getvalue()
    short.advance(1)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue() == quiet, "a pause of a second is not worth a line"
    short.advance(33)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue() == quiet + CLEAR + "  · working  34s"
    reporter.warn("a permanent line wipes it first")
    assert stream.getvalue().endswith("  · working  34s" + CLEAR + "  ! a permanent line wipes it first\n")
    short.advance(1)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue().endswith("\n" + CLEAR + "  · working  35s"), "and it comes back under the warning"


def test_close_wipes_the_status_line_for_good(short):
    reporter, stream = short(terminal=True)
    with reporter.stage("PDF report"):
        pass
    short.advance(5)
    reporter._refresh()  # noqa: SLF001
    assert stream.getvalue().endswith(CLEAR + "  · working  5s")
    reporter.close()
    assert stream.getvalue().endswith("  · working  5s" + CLEAR)
    before = stream.getvalue()
    short.advance(5)
    reporter._refresh()  # noqa: SLF001
    reporter.headline("✓ Done in 1m00s")
    assert stream.getvalue() == before + "✓ Done in 1m00s\n", "nothing is drawn after the run's last word"


def test_the_refresh_thread_redraws_the_running_stage(short):
    """The one place the real thread runs: it is what calls _refresh() in a real run."""
    reporter, stream = short(terminal=True, ticking=True)
    with reporter.stage("a long stage"):
        deadline = time.monotonic() + 30
        while stream.getvalue().count(CLEAR) < 3 and time.monotonic() < deadline:
            time.sleep(0.01)
        assert stream.getvalue().count(CLEAR) >= 3, "the stage's line is redrawn while it runs"
    assert "\n" not in stream.getvalue()


def test_the_refresh_thread_ends_with_its_reporter(short):
    reporter, _ = short(terminal=True, ticking=True)
    with reporter.stage("a stage"):
        ticker = reporter._ticker  # noqa: SLF001
        assert isinstance(ticker, threading.Thread) and ticker is not threading.current_thread()
        assert ticker.is_alive()
        assert ticker.daemon, "and it can never keep the program from exiting"
    del reporter
    gc.collect()
    ticker.join(timeout=30)
    assert not ticker.is_alive()


def test_the_status_line_never_wraps(short):
    reporter, stream = short(terminal=True, columns=40)
    with reporter.stage("a stage with a very long name indeed, longer than the window is wide"):
        reporter.info("x" * 500)
        reporter._refresh()  # noqa: SLF001
        short.advance(3700)
        reporter._refresh()  # noqa: SLF001
        drawn = stream.getvalue().split(CLEAR)[1:]
    assert len(drawn) == 3, "at the start, and at each of the two refreshes"
    assert drawn[0].startswith("  · a stage with a very long name") and drawn[2].endswith("  1h01m40s")
    assert max(len(line) for line in drawn) <= 39, [len(line) for line in drawn]


def test_colour_changes_nothing_but_the_colour(short):
    """The same run on a terminal with colour and without: the same text, once the colour is removed."""

    def a_run(reporter: Reporter) -> None:
        reporter.headline("S: lite report from 2 FASTQ files (1.0GB).")
        with reporter.stage("translated search (diamond blastx)"):
            reporter.info("  R1: 7,358 hit lines")
            short.advance(400)
            reporter._refresh()  # noqa: SLF001
        short.advance(3)
        reporter._refresh()  # noqa: SLF001
        reporter.warn("host index not built: human reads are not removed")
        with contextlib.suppress(RuntimeError), reporter.stage("unrepresented-lineage rescue (SingleM / GlobDB)"):
            raise RuntimeError(HMMSEARCH)
        reporter.error("PDF report failed: ValueError: x")
        reporter.close()

    plain_reporter, plain = short(terminal=True)
    a_run(plain_reporter)
    coloured_reporter, coloured = short(terminal=True, colour=True)
    a_run(coloured_reporter)
    assert not SGR.search(plain.getvalue())
    for code, what in (("36", "the dot of a running stage"), ("2", "its clock"), ("32", "the tick"),
                       ("33", "a warning's mark"), ("31", "a failure's cross")):
        assert f"\033[{code}m" in coloured.getvalue(), what
    assert SGR.sub("", coloured.getvalue()) == plain.getvalue()
    assert plain.getvalue().count("\n") == 5, "headline, stage, warning, failed stage, error"


@pytest.mark.parametrize(("terminal", "no_color", "expected"), [
    (True, None, True),      # an ordinary terminal
    (True, "1", False),      # one whose user has asked for no colour
    (False, None, False),    # a pipe or a log file: GitHub's runner, `openbiota run 2> log`
    (False, "1", False),
])
def test_colour_is_for_terminals_that_have_not_declined_it(monkeypatch, terminal: bool, no_color: str | None,
                                                           expected: bool):
    monkeypatch.setattr(logging_util.sys, "stderr", _Terminal() if terminal else io.StringIO())
    if no_color is None:
        monkeypatch.delenv("NO_COLOR", raising=False)
    else:
        monkeypatch.setenv("NO_COLOR", no_color)
    assert Reporter().style.enabled is expected


@pytest.mark.parametrize(("terminal", "term", "expected"), [
    (True, "xterm-256color", True),
    (True, "dumb", False),             # a terminal that cannot move the cursor
    (False, "xterm-256color", False),  # not a terminal at all
])
def test_a_line_is_redrawn_only_where_it_can_be(monkeypatch, terminal: bool, term: str, expected: bool):
    monkeypatch.setattr(logging_util.sys, "stderr", _Terminal() if terminal else io.StringIO())
    monkeypatch.setenv("TERM", term)
    assert Reporter(detail=False)._redraw is expected  # noqa: SLF001


def test_a_failed_stage_says_why_once_in_words(short):
    reporter, stream = short()
    label = "unrepresented-lineage rescue (SingleM / GlobDB)"
    try:
        with reporter.stage(label):
            raise RuntimeError(HMMSEARCH)
    except RuntimeError as exc:
        # what the lane runner does next: it explains the failure it caught
        reporter.warn(f"{label} failed: RuntimeError: {str(exc)[:200]}")
    assert stream.getvalue().splitlines() == [
        f"  ✗ {label}: the program 'hmmsearch' is not installed - run ./install.sh to add it"]


def test_a_warning_is_one_line_and_said_once(short):
    reporter, stream = short()
    for _ in range(3):
        reporter.warn("host index not built: human reads are not removed\nmore\nlines")
    reporter.warn("x" * 400)
    shown = stream.getvalue().splitlines()
    assert shown == ["  ! host index not built: human reads are not removed", "  ! " + "x" * 220 + "..."]


def test_an_unrelated_warning_after_a_failure_is_still_shown(short):
    reporter, stream = short()
    with pytest.raises(ValueError, match="broken"), reporter.stage("stage one"):
        raise ValueError("broken")
    reporter.warn("something else entirely")
    assert stream.getvalue().splitlines() == ["  ✗ stage one: broken", "  ! something else entirely"]


def test_an_error_shows_its_first_line_and_logs_the_rest(short):
    reporter, stream = short()
    reporter.error("PDF report failed: ValueError: x\nTraceback (most recent call last):\n  File ...")
    assert stream.getvalue().splitlines() == ["  ✗ PDF report failed: ValueError: x"]
    assert "Traceback" in "\n".join(reporter._lines)  # noqa: SLF001


def test_ctrl_c_in_a_stage_adds_no_failure_line(short):
    reporter, stream = short()
    with pytest.raises(KeyboardInterrupt), reporter.stage("a long stage"):
        raise KeyboardInterrupt
    assert stream.getvalue() == ""


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


def test_a_complete_run_ends_in_one_line(short):
    reporter, stream = short()
    cli._short_outcome(reporter, completeness=_stages(), total=3057.0,  # noqa: SLF001
                       results={"organism_inventory": {"n_organisms": 332}})
    assert stream.getvalue().splitlines() == [
        "", "✓ Done in 50m57s: lite report, all 19 stages ran, 332 organisms named."]


def test_an_incomplete_run_names_only_what_is_missing_and_why(short):
    reporter, stream = short()
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
