"""The legacy result objects are a contract, held by construction.

Spec 0.8.3 protects every pre-existing result object - scores, profiles,
organism calls, pathogen results, functional panels - from the extension.
These tests state that contract on synthetic data, so they hold on any
machine without a single sample: nothing here is a copy of a real
person's results.

AT001 / AT002 / AT003: with the extension attached, every protected leaf is
byte-identical to what it was before.
AT004: nothing new can enter a legacy object; the extension lives under its
own key and refuses to be anything else.
AT006: the command-line entrypoints and the legacy result keys consumers read
are still there.
"""

from __future__ import annotations

import copy
import json

import pytest

from openbiota.extension import engine, preservation


def _legacy_results() -> dict:
    """A results object with every protected key present and populated."""
    out: dict = {"sample": "SYNTHETIC_0001", "version": "0.8.4"}
    for i, key in enumerate(preservation.PROTECTED_KEYS):
        out[key] = {"score": 0.5 + i, "items": [{"id": f"{key}-{j}", "value": j * 1.5} for j in range(3)],
                    "label": f"{key} label", "nested": {"percentile": 42.0 + i, "band": "mid"}}
    return out


def test_at001_at002_at003_attaching_the_extension_leaves_every_protected_leaf_identical() -> None:
    results = _legacy_results()
    before = copy.deepcopy(results)
    engine.attach(results, engine.ExtensionResult(sample="SYNTHETIC_0001", mode="cached"))
    assert engine.RESULT_KEY in results
    report = preservation.compare(before, results)
    assert report.preserved, report.summary()
    assert report.compared_leaves > 0
    assert set(report.protected_keys_present) == set(preservation.PROTECTED_KEYS)
    for key in preservation.PROTECTED_KEYS:
        assert json.dumps(before[key], sort_keys=True) == json.dumps(results[key], sort_keys=True), key


def test_at004_the_extension_cannot_be_a_protected_object() -> None:
    assert engine.RESULT_KEY not in preservation.PROTECTED_KEYS
    # a protected object that the candidate has lost is a failure, never silently accepted
    before = _legacy_results()
    after = copy.deepcopy(before)
    del after["organism_inventory"]
    report = preservation.compare(before, after)
    assert not report.preserved and "organism_inventory" in report.protected_keys_missing
    # and a changed calibrated leaf inside a protected object is a difference the compare names
    after = copy.deepcopy(before)
    after["pathogens"]["nested"]["percentile"] = 99.0
    report = preservation.compare(before, after)
    assert not report.preserved
    assert any("pathogens" in d.path and "percentile" in d.path for d in report.differences)


def test_at006_the_entrypoints_and_legacy_keys_consumers_read_are_still_there() -> None:
    from openbiota.cli import build_parser

    parser = build_parser()
    sub = next(a for a in parser._actions if getattr(a, "choices", None) and "run" in a.choices)
    for command in ("run", "panels", "build-db", "doctor", "cohort", "validate", "depth-check", "run-all"):
        assert command in sub.choices, command
    # `run` still accepts the flags the documentation and the regeneration scripts use
    run_flags = {opt for action in sub.choices["run"]._actions for opt in action.option_strings}
    for flag in ("--r1", "--r2", "--sample", "--threads", "--simulate", "--no-pdf"):
        assert flag in run_flags, flag
    # the keys the report and the web interface read are the protected set
    for key in ("organism_inventory", "organism_verdicts", "pathogens", "profile_similarity",
                "findings_and_evidence", "biofilm", "mycobiome", "microbiome_age", "sequencing_quality"):
        assert key in preservation.PROTECTED_KEYS, key


@pytest.mark.parametrize("key", sorted(preservation.PROTECTED_KEYS))
def test_every_protected_key_is_a_plain_identifier(key: str) -> None:
    assert key.replace("_", "").isalnum() and key == key.lower()
