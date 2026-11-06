"""`openbiota run` from the shell: inputs, choices, modes, the last line.

What a person at a terminal meets, exercised without a terminal: the
positional FASTQ paths and how they resolve, the refusals for inputs that
are not FASTQ, the several-samples message when nobody can be asked, the
full/lite decision in every non-interactive case, the census per mode, and
the report path printed last. The two prompts are driven through a real
pty in ``test_interactive_prompts_on_a_pty``.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import pty
import select
import sys
import time
from pathlib import Path

import pytest

from openbiota import cli, completeness, reportmode
from openbiota.errors import InputError
from openbiota.search import validate_fastq_inputs

REPO = Path(__file__).resolve().parents[1]
READ = b"@LH00469:636:23VMTGLT4:8:1101:1003:20056 %d:N:0:ACGT\nACGTACGTAC\n+\nJJJJJJJJJJ\n"


def _fastq(path: Path, mate: int = 1, n: int = 3, *, gz: bool = True) -> Path:
    body = READ.replace(b"%d", str(mate).encode()) * n
    path.write_bytes(gzip.compress(body) if gz else body)
    return path


def _args(**over) -> argparse.Namespace:
    ns = cli.build_parser().parse_args(["run", *over.pop("argv", [])])
    for k, v in over.items():
        setattr(ns, k, v)
    return ns


# --------------------------------------------------------------------------- #
# positional inputs
# --------------------------------------------------------------------------- #


def test_two_files_become_the_pair_mate_one_first(tmp_path: Path):
    r2 = _fastq(tmp_path / "S_2.fastq.gz", 2)
    r1 = _fastq(tmp_path / "S_1.fastq.gz", 1)
    args = _args(argv=[str(r2), str(r1)])          # given the wrong way round
    cli._resolve_inputs(args)
    assert (args.r1, args.r2) == (r1, r2)


def test_one_file_finds_its_mate_beside_it(tmp_path: Path):
    r1 = _fastq(tmp_path / "S_R1.fastq.gz", 1)
    r2 = _fastq(tmp_path / "S_R2.fastq.gz", 2)
    args = _args(argv=[str(r1)])
    cli._resolve_inputs(args)
    assert (args.r1, args.r2) == (r1, r2)


def test_one_file_without_a_mate_is_single_end(tmp_path: Path):
    r1 = _fastq(tmp_path / "lonely_1.fastq.gz", 1)
    args = _args(argv=[str(r1)])
    cli._resolve_inputs(args)
    assert (args.r1, args.r2) == (r1, None)


def test_a_folder_becomes_the_fastq_dir(tmp_path: Path):
    args = _args(argv=[str(tmp_path)])
    cli._resolve_inputs(args)
    assert args.fastq_dir == tmp_path and args.r1 is None


@pytest.mark.parametrize("n", [3, 4])
def test_more_than_two_files_is_refused_with_the_way_out(tmp_path: Path, n: int):
    files = [_fastq(tmp_path / f"s{i}_1.fastq.gz") for i in range(n)]
    with pytest.raises(InputError, match="run-all"):
        cli._resolve_inputs(_args(argv=[str(f) for f in files]))


def test_positional_and_flag_inputs_cannot_be_mixed(tmp_path: Path):
    r1 = _fastq(tmp_path / "S_1.fastq.gz")
    with pytest.raises(InputError, match="not both"):
        cli._resolve_inputs(_args(argv=[str(r1), "--r1", str(r1)]))


# --------------------------------------------------------------------------- #
# what is refused before an hour is spent on it
# --------------------------------------------------------------------------- #


def test_fasta_is_named_as_such(tmp_path: Path):
    p = tmp_path / "x.fa"
    p.write_text(">seq\nACGT\n")
    with pytest.raises(InputError, match="FASTA"):
        validate_fastq_inputs([p])


@pytest.mark.parametrize("content,word", [
    (b"<!DOCTYPE html><html>", "web page"),
    (b"%PDF-1.4 ...", "PDF"),
    (b"hello world\n", "first line is 'hello world'"),
])
def test_other_files_are_named_for_what_they_are(tmp_path: Path, content: bytes, word: str):
    p = tmp_path / "x.fastq"
    p.write_bytes(content)
    with pytest.raises(InputError, match=word):
        validate_fastq_inputs([p])


def test_empty_and_missing_files_are_refused(tmp_path: Path):
    empty = tmp_path / "e.fastq.gz"
    empty.write_bytes(b"")
    with pytest.raises(InputError, match="empty"):
        validate_fastq_inputs([empty])
    with pytest.raises(InputError, match="not found"):
        validate_fastq_inputs([tmp_path / "nope.fastq.gz"])


def test_a_bad_gzip_is_explained(tmp_path: Path):
    p = tmp_path / "x.fastq.gz"
    p.write_bytes(b"\x1f\x8b\x08\x00garbage-not-a-stream")
    with pytest.raises(InputError, match="cannot be read"):
        validate_fastq_inputs([p])


def test_two_different_samples_are_not_a_pair(tmp_path: Path):
    a = _fastq(tmp_path / "A_1.fastq.gz", 1)
    b = tmp_path / "B_2.fastq.gz"
    b.write_bytes(gzip.compress(b"@OTHER:1:1:1:1 2:N:0:A\nACGT\n+\nJJJJ\n"))
    with pytest.raises(InputError, match="two mates of one sample"):
        validate_fastq_inputs([a, b])


def test_a_real_pair_passes(tmp_path: Path):
    validate_fastq_inputs([_fastq(tmp_path / "S_1.fastq.gz", 1), _fastq(tmp_path / "S_2.fastq.gz", 2)])
    validate_fastq_inputs([_fastq(tmp_path / "plain_1.fastq", 1, gz=False)])


# --------------------------------------------------------------------------- #
# which sample, when nobody can be asked
# --------------------------------------------------------------------------- #


def _two_samples(tmp_path: Path) -> Path:
    for name in ("alice", "bob"):
        d = tmp_path / name
        d.mkdir()
        _fastq(d / f"{name}_1.fastq.gz", 1)
        _fastq(d / f"{name}_2.fastq.gz", 2)
    return tmp_path


def test_several_samples_without_a_terminal_says_how_to_choose(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    args = _args(fastq_dir=_two_samples(tmp_path))
    with pytest.raises(InputError) as excinfo:
        cli._select_samples(args)
    text = str(excinfo.value)
    assert "alice, bob" in text and "--sample NAME" in text and "run-all" in text


def test_sample_flag_picks_one_and_names_the_rest_when_wrong(tmp_path: Path):
    args = _args(fastq_dir=_two_samples(tmp_path), sample="bob")
    assert [s.sample for s in cli._select_samples(args)] == ["bob"]
    with pytest.raises(InputError, match="alice, bob"):
        cli._select_samples(_args(fastq_dir=tmp_path, sample="carol"))


def test_one_sample_needs_no_question(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: pytest.fail("should not be consulted"))
    d = tmp_path / "only"
    d.mkdir()
    _fastq(d / "only_1.fastq.gz", 1)
    _fastq(d / "only_2.fastq.gz", 2)
    assert [s.sample for s in cli._select_samples(_args(fastq_dir=tmp_path))] == ["only"]


# --------------------------------------------------------------------------- #
# full or lite
# --------------------------------------------------------------------------- #


@pytest.fixture(autouse=True)
def _standard_install_present(monkeypatch):
    """Unless a test says otherwise, the machine has the standard install."""
    monkeypatch.setattr(reportmode, "standard_install_missing", lambda _repo=None: [])


def test_lite_switches_the_full_only_lanes_off_and_keeps_the_catalogues(monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    args = _args(argv=["--mode", "lite"])
    assert cli._apply_mode(args, allow_prompt=True) == "lite"
    assert args.no_extended_catalogue and args.no_pathogens and args.no_mycobiome
    assert not getattr(args, "no_expanded_lanes", False), "the catalogue lanes run in lite"
    # Settled once: a second call changes nothing.
    args.mode = "full"
    assert cli._apply_mode(args, allow_prompt=True) == "lite"


def test_auto_is_lite_without_the_full_install(monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    monkeypatch.setattr(reportmode, "full_install_missing", lambda _repo=None: ["pathogens", "mycobiome_refs"])
    args = _args(argv=["--mode", "auto"])
    assert cli._apply_mode(args, allow_prompt=True) == "lite"
    assert "not complete" in args.mode_reason and "pathogens" in args.mode_reason


def test_auto_is_full_with_the_full_install(monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    monkeypatch.setattr(reportmode, "full_install_missing", lambda _repo=None: [])
    args = _args()   # no --mode, no terminal: auto
    assert cli._apply_mode(args, allow_prompt=True) == "full"
    assert not args.no_pathogens


def test_full_without_the_install_is_refused_with_the_command(monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    monkeypatch.setattr(reportmode, "full_install_missing", lambda _repo=None: ["pathogens"])
    with pytest.raises(InputError, match="openbiota install --full"):
        cli._apply_mode(_args(argv=["--mode", "full"]), allow_prompt=True)


def test_a_run_refuses_to_start_without_the_standard_install(monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    monkeypatch.setattr(reportmode, "standard_install_missing", lambda _repo=None: ["catalogue_refs"])
    with pytest.raises(InputError, match="catalogue_refs"):
        cli._apply_mode(_args(argv=["--mode", "lite"]), allow_prompt=True)


def test_mode_of_reads_the_record_and_defaults_to_full():
    assert reportmode.mode_of({"run": {"report_mode": "lite"}}) == "lite"
    assert reportmode.mode_of({"run": {"report_mode": "FULL"}}) == "full"
    assert reportmode.mode_of({}) == "full"
    assert reportmode.mode_of({"run": {"report_mode": "banana"}}) == "full"


def test_the_split_is_consistent():
    """Every lite-excluded stage is one the full-only sections depend on,
    and every full-only section has a stage behind it."""
    keys = {s.key for s in completeness.REQUIRED}
    assert keys >= reportmode.FULL_ONLY_STAGES
    assert reportmode.excluded_stages("full") == frozenset()
    assert reportmode.omitted_sections("full") == frozenset()
    assert {"pathogens", "mycobiome", "strain_typing"} <= reportmode.FULL_ONLY_STAGES
    assert {"pathogens", "mycobiome", "strains"} <= reportmode.FULL_ONLY_SECTIONS
    # The five 2025 catalogue lanes, the whole-genome profile, confirmation and
    # assembly belong to every report: the standard install carries them.
    for key in ("lane_jan26", "lane_globdb", "lane_motus", "lane_kraken", "lane_singlem",
                "genome_profile", "confirmation", "assembly"):
        assert key not in reportmode.FULL_ONLY_STAGES, key


def test_the_census_asks_a_lite_run_only_about_lite_stages():
    lite = completeness.audit({"run": {"report_mode": "lite"}})
    full = completeness.audit({})
    assert lite.mode == "lite" and full.mode == "full"
    assert len(lite.stages) + len(lite.not_in_mode) == len(full.stages) == len(completeness.REQUIRED)
    assert "pathogen screen" in lite.not_in_mode
    assert len(lite.not_in_mode) == 5 and len(lite.stages) == 19
    assert lite.to_json()["not_in_mode"] == lite.not_in_mode


def test_zero_triggers_is_a_result_not_a_missing_stage():
    r = {"findings_and_evidence": {"triggers": [], "findings": {"x": 1}}}
    ran, _ = completeness._findings(r)
    assert ran
    ran, _ = completeness._findings({"findings_and_evidence": {"triggers": [], "findings": {}}})
    assert not ran


# --------------------------------------------------------------------------- #
# the last line, and the viewer
# --------------------------------------------------------------------------- #


def test_the_report_path_is_the_last_thing_printed(tmp_path: Path, monkeypatch, capsys):
    pdf = tmp_path / "results" / "S" / "S_report.pdf"
    pdf.parent.mkdir(parents=True)
    pdf.write_bytes(b"%PDF-1.4")
    r1 = _fastq(tmp_path / "S_1.fastq.gz", 1)
    r2 = _fastq(tmp_path / "S_2.fastq.gz", 2)

    def fake_run(args, **_kw):
        print("... many lines of summary ...")
        args._last_report = pdf
        return 0

    monkeypatch.setattr(cli, "cmd_run", fake_run)
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    opened: list[Path] = []
    monkeypatch.setattr(cli, "_open_report", lambda p: opened.append(p) or True)
    args = _args(argv=[str(r1), str(r2), "--mode", "lite"])
    assert cli.cmd_run_entry(args) == 0
    lines = [line for line in capsys.readouterr().out.splitlines() if line.strip()]
    assert lines[-1] == f"Report: {pdf.resolve()}"
    assert opened == [], "no terminal, so nothing is opened"


def test_the_viewer_opens_only_at_a_terminal_and_only_once(tmp_path: Path, monkeypatch, capsys):
    pdf = tmp_path / "S_report.pdf"
    pdf.write_bytes(b"%PDF-1.4")
    r1 = _fastq(tmp_path / "S_1.fastq.gz", 1)
    monkeypatch.setattr(cli, "cmd_run", lambda args, **_kw: setattr(args, "_last_report", pdf) or 0)
    monkeypatch.setattr(cli, "_interactive", lambda: True)
    opened: list[Path] = []
    monkeypatch.setattr(cli, "_open_report", lambda p: opened.append(p) or True)
    assert cli.cmd_run_entry(_args(argv=[str(r1), "--mode", "lite"])) == 0
    assert opened == [pdf]
    assert "Opening it in your PDF viewer." in capsys.readouterr().out
    opened.clear()
    assert cli.cmd_run_entry(_args(argv=[str(r1), "--mode", "lite", "--no-open"])) == 0
    assert opened == []


def test_open_report_respects_the_environment(monkeypatch, tmp_path: Path):
    launched: list[list[str]] = []
    monkeypatch.setattr(cli.subprocess, "Popen", lambda cmd, **_kw: launched.append(cmd))
    pdf = tmp_path / "x.pdf"
    monkeypatch.setenv("OPENBIOTA_NO_OPEN", "1")
    assert cli._open_report(pdf) is False and launched == []
    monkeypatch.delenv("OPENBIOTA_NO_OPEN")
    monkeypatch.setenv("CI", "true")
    assert cli._open_report(pdf) is False
    monkeypatch.delenv("CI")
    monkeypatch.setattr(cli.sys, "platform", "linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    assert cli._open_report(pdf) is False, "no desktop, nothing to open with"
    monkeypatch.setenv("DISPLAY", ":0")
    monkeypatch.setattr(cli.shutil, "which", lambda _name: "/usr/bin/xdg-open")
    assert cli._open_report(pdf) is True and launched[-1][0] == "/usr/bin/xdg-open"
    monkeypatch.setattr(cli.sys, "platform", "darwin")
    assert cli._open_report(pdf) is True and launched[-1][0] == "open"


def test_a_failed_run_says_so_instead_of_a_path(tmp_path: Path, monkeypatch, capsys):
    r1 = _fastq(tmp_path / "S_1.fastq.gz", 1)
    monkeypatch.setattr(cli, "cmd_run", lambda *_a, **_kw: 1)
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    assert cli.cmd_run_entry(_args(argv=[str(r1), "--mode", "lite"])) == 1
    assert "No report for S" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# the two prompts, on a real pty
# --------------------------------------------------------------------------- #


def _drive(argv: list[str], answers: list[tuple[str, str]], cwd: Path, timeout: float = 30.0) -> str:
    """Run argv on a pty, answer each prompt when it appears, return everything printed."""
    pid, fd = pty.fork()
    if pid == 0:  # pragma: no cover - the child
        os.chdir(cwd)
        os.environ["OPENBIOTA_NO_OPEN"] = "1"
        os.execv(argv[0], argv)
    out = b""
    i = 0
    deadline = time.time() + timeout
    while time.time() < deadline:
        ready, _, _ = select.select([fd], [], [], 0.3)
        if not ready:
            continue
        try:
            chunk = os.read(fd, 65536)
        except OSError:
            break
        if not chunk:
            break
        out += chunk
        if i < len(answers) and out.rstrip().endswith(answers[i][0].encode()):
            os.write(fd, (answers[i][1] + "\n").encode())
            i += 1
    try:
        os.waitpid(pid, os.WNOHANG)
        os.kill(pid, 9)
    except (ChildProcessError, ProcessLookupError):
        pass
    return out.decode("utf-8", "replace")


@pytest.mark.skipif(not hasattr(os, "fork"), reason="needs a pty")
def test_interactive_prompts_on_a_pty(tmp_path: Path):
    """With two samples and no flags, a person is asked which sample and
    which report; 'q' at either prompt stops cleanly."""
    _two_samples(tmp_path)
    exe = str(REPO / ".venv" / "bin" / "openbiota")
    if not Path(exe).exists():
        pytest.skip("no .venv in this checkout")
    text = _drive([exe, "run", "--fastq-dir", str(tmp_path)], [("sample:", "2"), ("[1]", "q"), ("[2]", "q")], cwd=REPO)
    assert "Several samples are in the FASTQ folder" in text
    assert " 1)  alice" in text and " 2)  bob" in text
    assert "Which report?" in text and "1) full" in text and "2) lite" in text
    assert "the first time" in text, "each mode states how long it takes"
    assert "no report mode chosen" in text


def test_the_published_sample_reports_are_made_in_full_mode():
    script = (REPO / "scripts" / "generate_sample_reports.py").read_text()
    assert '"--mode", "full"' in script and '"--no-open"' in script


def test_results_record_the_mode_for_the_census(tmp_path: Path):
    payload = {"run": {"report_mode": "lite"}}
    (tmp_path / "results.json").write_text(json.dumps(payload))
    assert completeness.audit(json.loads((tmp_path / "results.json").read_text())).mode == "lite"


def test_bare_launcher_prints_help_not_a_run():
    launcher = (REPO / "bin" / "openbiota").read_text()
    assert "set -- run" not in launcher


def test_python_version_marker():
    assert sys.version_info >= (3, 11)
