"""Every stage of the report ran, for every report on disk.

One report went out saying "strain typing did not run" on page one, from
a run that exited 0. Another's simulation was never solved. The
completeness gate exists so that cannot happen quietly again; these tests
hold the gate itself to account, and then hold every stored report to it.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openbiota.completeness import REQUIRED, Completeness, audit, describe
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
RESULTS = sorted(
    p for p in (REPO / "results").glob("*/results.json")
    if not p.parent.name.startswith(("_", "fmt_"))
)
SAMPLES = [p.parent.name for p in RESULTS]


def _load(sample: str) -> dict:
    return json.loads((results_file(sample)).read_text())


# --------------------------------------------------------------------------- #
# the gate itself
# --------------------------------------------------------------------------- #


def test_every_required_stage_has_a_distinct_key_and_label() -> None:
    keys = [s.key for s in REQUIRED]
    labels = [s.label for s in REQUIRED]
    assert len(set(keys)) == len(keys)
    assert len(set(labels)) == len(labels)


def test_an_empty_results_object_fails_every_stage() -> None:
    """Absence is never mistaken for success."""
    report = audit({})
    assert not report.complete
    assert len(report.missing) == len(REQUIRED)


def test_a_broken_check_is_a_failed_stage_not_a_crash() -> None:
    report = audit({"panels": "not a list"})
    assert any(not s["ran"] for s in report.stages)


@pytest.mark.skipif(not SAMPLES, reason="no results on disk")
def test_strain_typing_not_run_is_caught() -> None:
    """The exact failure that motivated the gate."""
    results = copy.deepcopy(_load(SAMPLES[0]))
    results["strain_resolution"] = {"status": "not_run"}
    report = audit(results)
    assert "strain typing (marker consensus)" in report.missing


@pytest.mark.skipif(not SAMPLES, reason="no results on disk")
def test_an_unsolved_simulation_is_caught() -> None:
    """The other failure that motivated the gate: a cached-mode dry run."""
    results = copy.deepcopy(_load(SAMPLES[0]))
    results["extension"]["views"]["simulation"]["scenarios"] = {"execution_state": "not_run"}
    report = audit(results)
    assert "model-assisted simulation" in report.missing


@pytest.mark.skipif(not SAMPLES, reason="no results on disk")
def test_an_abstained_age_still_counts_as_having_run() -> None:
    """Abstaining on an out-of-distribution sample is the model working."""
    results = copy.deepcopy(_load(SAMPLES[0]))
    results["microbiome_age"]["status"] = "abstained_ood"
    assert "estimated age of biota" not in audit(results).missing
    results["microbiome_age"]["status"] = "not_computable"
    assert "estimated age of biota" in audit(results).missing


def test_the_json_record_carries_what_the_reader_needs() -> None:
    report = Completeness([
        {"key": "a", "label": "stage a", "ran": True, "detail": "ok"},
        {"key": "b", "label": "stage b", "ran": False, "detail": "absent"},
    ])
    record = report.to_json()
    assert record["complete"] is False
    assert record["n_ran"] == 1 and record["n_required"] == 2
    assert record["missing"] == ["stage b"]
    text = describe(report)
    assert "MISSING" in text and "stage b" in text


def test_the_run_exits_non_zero_when_incomplete() -> None:
    """The gate has to reach the exit code, or it is a log line."""
    src = (REPO / "openbiota" / "cli.py").read_text()
    assert "if not completeness.complete:\n        return 2" in src
    assert '["completeness"] = completeness.to_json()' in src


def test_strain_typing_is_a_stage_of_the_run_not_a_second_pass() -> None:
    """The reason one report said it did not run: nobody ran it."""
    src = (REPO / "openbiota" / "cli.py").read_text()
    assert 'reporter.stage("strain typing (marker consensus)")' in src
    assert "run_sample2markers(" in src
    assert "sam_out=marker_sam" in src, "the SAM must come off the MetaPhlAn 4 pass"


# --------------------------------------------------------------------------- #
# every report on disk
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_stage_ran_for(sample: str) -> None:
    report = audit(_load(sample))
    assert report.complete, (
        f"{sample}: {len(report.missing)} stage(s) did not produce a result\n"
        + describe(report)
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_run_recorded_its_own_completeness(sample: str) -> None:
    """A report built after the gate carries the gate's verdict in the file."""
    record = (_load(sample).get("run") or {}).get("completeness")
    assert record, f"{sample} was built before the completeness gate; regenerate it"
    assert record["complete"] is True, f"{sample}: {record['missing']}"
