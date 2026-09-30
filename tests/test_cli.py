"""CLI argument handling, command dispatch and error presentation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.cli import (
    _flags_for_display,
    _split_list,
    build_parser,
    main,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- #
# implicit `run`
# --------------------------------------------------------------------------- #


def test_bare_flags_imply_the_run_subcommand():
    parser = build_parser()
    args = parser.parse_args(["run", "--subsample", "1000"])
    assert args.command == "run"
    assert args.subsample == 1000


@pytest.mark.parametrize("argv", [[], ["--threads", "4"], ["--compact"]])
def test_main_inserts_run_for_bare_invocations(argv, monkeypatch):
    seen: dict[str, object] = {}

    def fake_run(args):
        seen["command"] = args.command
        return 0

    monkeypatch.setattr("openbiota.cli.cmd_run", fake_run)
    assert main(argv) == 0
    assert seen["command"] is None or seen["command"] == "run"


def test_explicit_subcommands_are_not_rewritten(capsys):
    assert main(["panels"]) == 0
    out = capsys.readouterr().out
    assert "normaliser:" in out
    assert "urda" in out


# --------------------------------------------------------------------------- #
# list parsing
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, None),
        ("", None),
        ("all", None),
        ("ALL", None),
        ("urda", ["urda"]),
        ("urda,cutc", ["urda", "cutc"]),
        (" urda , cutc ", ["urda", "cutc"]),
        ("urda,,cutc", ["urda", "cutc"]),
    ],
)
def test_split_list(value, expected):
    assert _split_list(value) == expected


# --------------------------------------------------------------------------- #
# flag display
# --------------------------------------------------------------------------- #


def test_flags_for_display_elides_per_run_paths():
    command = (
        "/usr/local/bin/diamond", "blastx",
        "--db", "/abs/path/combined.dmnd",
        "--query", "/abs/path/A02_1.fastq.gz",
        "--out", "/abs/path/hits.tsv",
        "--max-target-seqs", "1",
        "--threads", "32",
        "--quiet",
    )
    shown = _flags_for_display(command)

    assert shown.startswith("blastx")
    assert "--max-target-seqs 1" in shown
    assert "--threads 32" in shown
    assert "--quiet" in shown
    # paths and their flags are gone, and the executable path is not shown
    for absent in ("--db", "--query", "--out", "/abs/path", "/usr/local/bin"):
        assert absent not in shown


def test_flags_for_display_handles_empty():
    assert _flags_for_display(()) == "n/a"


# --------------------------------------------------------------------------- #
# panels command
# --------------------------------------------------------------------------- #


def test_panels_lists_excluded_genes(capsys):
    assert main(["panels"]) == 0
    out = capsys.readouterr().out
    assert "Deliberately excluded" in out
    for gene in ("ldh", "gadB", "speE", "murI"):
        assert gene in out
    assert "urolithin a" in out.lower()


def test_panels_show_prints_queries(capsys):
    assert main(["panels", "--show", "urda"]) == 0
    out = capsys.readouterr().out
    assert "urocanate reductase" in out
    assert "Q8CVD0" in out
    assert "373" in out
    assert "upper bound" in out


def test_panels_show_json_is_valid(capsys):
    assert main(["panels", "--show", "cutc", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "cutc"
    assert {t["id"] for t in payload["targets"]} == {"CUTC", "CNTA"}
    assert payload["decoys"]


def test_panels_json_lists_everything(capsys):
    assert main(["panels", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["normalizer"]["id"] == "RPOB"
    assert len(payload["panels"]) >= 8


def test_unknown_panel_exits_cleanly(capsys):
    assert main(["panels", "--show", "nosuchpanel"]) == 2
    err = capsys.readouterr().err
    assert "unknown panel" in err
    assert "available:" in err
    assert "Traceback" not in err


# --------------------------------------------------------------------------- #
# error presentation
# --------------------------------------------------------------------------- #


def test_missing_fastq_directory_exits_cleanly(capsys, tmp_path: Path):
    code = main(["run", "--fastq-dir", str(tmp_path / "absent"), "--panels-dir", str(REPO_ROOT / "panels")])
    assert code == 2
    err = capsys.readouterr().err
    assert "openbiota:" in err
    assert "not found" in err
    assert "Traceback" not in err


def test_malformed_panels_directory_exits_cleanly(capsys, tmp_path: Path):
    (tmp_path / "_normalizer.yaml").write_text("normalizer: {}\n", encoding="utf-8")
    code = main(["run", "--panels-dir", str(tmp_path)])
    assert code == 2
    err = capsys.readouterr().err
    assert "openbiota:" in err
    assert "Traceback" not in err


# --------------------------------------------------------------------------- #
# defaults
# --------------------------------------------------------------------------- #


def test_run_defaults_match_the_documented_diamond_flags():
    args = build_parser().parse_args(["run"])
    assert args.block_size == 8.0
    assert args.index_chunks == 1
    assert args.sensitivity == "default"
    assert args.evalue == 1e-5
    assert args.threads >= 1
    assert args.panels == "all"
    assert args.db_panels == "all"
    assert args.subsample is None
    assert args.fastq_dir.name == "fastq"


def test_doctor_runs_and_reports(capsys):
    code = main(["doctor"])
    out = capsys.readouterr().out
    assert "python" in out
    assert "diamond" in out
    assert "panels" in out
    assert code in (0, 1)
