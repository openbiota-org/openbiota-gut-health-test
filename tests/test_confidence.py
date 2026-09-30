"""Detection confidence tiers and the thresholds behind them."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.confidence import (
    CONFIRMED,
    INDETERMINATE,
    MIN_FRAGMENTS,
    MIN_MEAN_IDENTITY,
    NOT_DETECTED,
    PROVISIONAL,
    classify_confidence,
    reasons_below_confirmed,
)

VALIDATION = Path(__file__).resolve().parent.parent / "docs" / "validation_results.json"


# --------------------------------------------------------------------------- #
# tier assignment
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("fragments", "identity", "expected"),
    [
        (0, None, NOT_DETECTED),
        (0, 99.0, NOT_DETECTED),
        (1, 99.0, PROVISIONAL),          # a single read is not evidence
        (2, 99.0, CONFIRMED),
        (2, MIN_MEAN_IDENTITY, CONFIRMED),
        (2, MIN_MEAN_IDENTITY - 0.1, PROVISIONAL),
        (1000, 74.9, PROVISIONAL),        # many reads, but distant homologs
        (1366, 82.9, CONFIRMED),
        (52, 68.3, PROVISIONAL),          # the real urdA false positive
        (144, 66.2, PROVISIONAL),         # the real lcdA false positive
        (64, 70.6, PROVISIONAL),          # the true positive given up at this tier
        (2, None, PROVISIONAL),
    ],
)
def test_tier_assignment(fragments, identity, expected):
    assert classify_confidence(fragments=fragments, mean_identity=identity) == expected


def test_indeterminate_overrides_everything():
    assert (
        classify_confidence(fragments=5000, mean_identity=99.0, indeterminate=True)
        == INDETERMINATE
    )


def test_thresholds_are_the_documented_values():
    """These are the operating point, so a silent change would invalidate the docs."""
    assert MIN_FRAGMENTS == 2
    assert MIN_MEAN_IDENTITY == 75.0


# --------------------------------------------------------------------------- #
# explanations
# --------------------------------------------------------------------------- #


def test_reasons_name_the_failing_criterion():
    assert "1 matching fragment" in " ".join(reasons_below_confirmed(1, 99.0))
    assert "identity" in " ".join(reasons_below_confirmed(500, 60.0))

    both = reasons_below_confirmed(1, 60.0)
    assert len(both) == 2

    assert reasons_below_confirmed(500, 99.0) == []


def test_reasons_handle_missing_identity():
    assert reasons_below_confirmed(500, None)


# --------------------------------------------------------------------------- #
# the thresholds against the actual validation data
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not VALIDATION.is_file(), reason="run `openbiota validate` first")
def test_confirmed_tier_has_no_false_positives_in_validation():
    """The whole justification for these thresholds, checked against the data."""
    payload = json.loads(VALIDATION.read_text(encoding="utf-8"))
    confirmed = payload["operating_points"]["confirmed"]

    assert confirmed["confusion"]["FP"] == 0, (
        "the confirmed tier is defined as the point with zero false positives; "
        f"validation now reports {confirmed['confusion']['FP']}"
    )
    assert confirmed["metrics"]["specificity"] == 1.0
    assert confirmed["metrics"]["precision"] == 1.0
    assert confirmed["metrics"]["sensitivity"] >= 0.90


@pytest.mark.skipif(not VALIDATION.is_file(), reason="run `openbiota validate` first")
def test_confirmed_tier_beats_the_permissive_one_on_specificity():
    payload = json.loads(VALIDATION.read_text(encoding="utf-8"))
    permissive = payload["operating_points"]["permissive"]
    confirmed = payload["operating_points"]["confirmed"]

    assert confirmed["confusion"]["FP"] < permissive["confusion"]["FP"]
    assert confirmed["metrics"]["specificity"] > permissive["metrics"]["specificity"]


@pytest.mark.skipif(not VALIDATION.is_file(), reason="run `openbiota validate` first")
def test_every_validation_false_positive_is_below_the_identity_floor():
    """The empirical basis for the identity criterion."""
    payload = json.loads(VALIDATION.read_text(encoding="utf-8"))
    permissive_fps = [
        outcome
        for experiment in payload["experiments"]
        for outcome in experiment["genes"]
        if outcome["classification"] == "FP"
    ]
    assert permissive_fps, "expected the permissive tier to have some false positives"
    for outcome in permissive_fps:
        assert classify_confidence(
            fragments=outcome["reported_fragments"],
            mean_identity=outcome["mean_identity"],
        ) == PROVISIONAL, f"{outcome['entry_key']} would be reported as confirmed"


@pytest.mark.skipif(not VALIDATION.is_file(), reason="run `openbiota validate` first")
def test_residue_filter_improves_accuracy():
    """The reference-side residue check must measurably reduce overcounting."""
    payload = json.loads(VALIDATION.read_text(encoding="utf-8"))
    effect = payload["residue_filter_effect"]
    headline = effect["headline_vs_residue_filtered_truth"]["slope"]
    strict = effect["residue_consistent_vs_residue_filtered_truth"]["slope"]

    assert headline is not None and strict is not None
    # both overcount, but the filtered figure must be closer to 1.0
    assert abs(strict - 1.0) < abs(headline - 1.0)


@pytest.mark.skipif(not VALIDATION.is_file(), reason="run `openbiota validate` first")
def test_detection_limit_reaches_low_abundance():
    payload = json.loads(VALIDATION.read_text(encoding="utf-8"))
    limits = payload["detection_limit"]
    assert limits, "expected detection-limit results"
    assert all(row["detected"] for row in limits), "a spiked carrier was missed entirely"
    assert min(row["carrier_cell_fraction"] for row in limits) <= 0.002
