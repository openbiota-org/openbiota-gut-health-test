"""ROC machinery and the cross-profile specificity verdicts."""

from __future__ import annotations

import pytest

from openbiota.profilevalidate import (
    AucResult,
    LabelledScores,
    SpecificityMatrix,
    auc_confidence_interval,
    evaluate_auc,
    roc_auc,
    spread_verdict,
)

# --------------------------------------------------------------------------- #
# AUC
# --------------------------------------------------------------------------- #


def test_perfect_separation_is_one() -> None:
    assert roc_auc([1.0, 2.0, 3.0, 4.0], [0, 0, 1, 1]) == pytest.approx(1.0)


def test_reversed_separation_is_zero() -> None:
    assert roc_auc([4.0, 3.0, 2.0, 1.0], [0, 0, 1, 1]) == pytest.approx(0.0)


def test_all_ties_is_a_coin_flip() -> None:
    assert roc_auc([1.0] * 6, [1, 1, 1, 0, 0, 0]) == pytest.approx(0.5)


def test_partial_ties_are_corrected() -> None:
    # cases {2,2}, controls {1,2}: one clear win, two ties
    assert roc_auc([2.0, 2.0, 1.0, 2.0], [1, 1, 0, 0]) == pytest.approx(0.75)


def test_single_class_has_no_auc() -> None:
    assert roc_auc([1.0, 2.0, 3.0], [1, 1, 1]) is None
    assert roc_auc([1.0, 2.0, 3.0], [0, 0, 0]) is None


def test_auc_matches_the_mann_whitney_definition() -> None:
    """AUC is the probability a random case outranks a random control."""
    cases = [3.0, 5.0, 7.0]
    controls = [1.0, 4.0, 6.0, 8.0]
    wins = sum(
        1.0 if c > k else 0.5 if c == k else 0.0 for c in cases for k in controls
    )
    expected = wins / (len(cases) * len(controls))
    got = roc_auc(cases + controls, [1] * len(cases) + [0] * len(controls))
    assert got == pytest.approx(expected)


def test_confidence_interval_brackets_the_estimate() -> None:
    import random

    rng = random.Random(11)
    cases = [rng.gauss(1.0, 1.0) for _ in range(120)]
    controls = [rng.gauss(0.0, 1.0) for _ in range(120)]
    scores = cases + controls
    labels = [1] * 120 + [0] * 120
    auc = roc_auc(scores, labels)
    low, high = auc_confidence_interval(scores, labels, replicates=400)
    assert low is not None and high is not None
    assert low < auc < high
    assert high - low < 0.30


def test_confidence_interval_needs_enough_data() -> None:
    assert auc_confidence_interval([1.0, 2.0], [0, 1]) == (None, None)


# --------------------------------------------------------------------------- #
# per-study evaluation
# --------------------------------------------------------------------------- #


def _scores() -> LabelledScores:
    out = LabelledScores(profile="testprofile", condition="TEST")
    # study A separates perfectly, study B not at all
    for value, label in ((3.0, 1), (4.0, 1), (1.0, 0), (2.0, 0)):
        out.scores.append(value)
        out.labels.append(label)
        out.studies.append("StudyA")
        out.sample_ids.append(f"A{len(out.scores)}")
    for value, label in ((1.0, 1), (2.0, 1), (1.0, 0), (2.0, 0)):
        out.scores.append(value)
        out.labels.append(label)
        out.studies.append("StudyB")
        out.sample_ids.append(f"B{len(out.scores)}")
    return out


def test_case_and_control_counts() -> None:
    scores = _scores()
    assert scores.n_cases == 4
    assert scores.n_controls == 4


def test_per_study_auc_is_reported_separately() -> None:
    result = evaluate_auc(_scores())
    assert result.per_study["StudyA"]["auc"] == pytest.approx(1.0)
    assert result.per_study["StudyB"]["auc"] == pytest.approx(0.5)


def test_study_holdout_mean_averages_the_studies() -> None:
    """Pooling and splitting randomly would leak recruitment-batch signal."""
    result = evaluate_auc(_scores())
    assert result.holdout_auc == pytest.approx(0.75)


def test_single_class_study_is_flagged_not_scored() -> None:
    scores = LabelledScores(profile="p", condition="c")
    for value in (1.0, 2.0, 3.0):
        scores.scores.append(value)
        scores.labels.append(1)
        scores.studies.append("OnlyCases")
        scores.sample_ids.append("x")
    for value in (0.5, 0.7):
        scores.scores.append(value)
        scores.labels.append(0)
        scores.studies.append("OnlyControls")
        scores.sample_ids.append("y")
    result = evaluate_auc(scores)
    assert result.per_study["OnlyCases"]["auc"] is None
    assert "single-class" in result.per_study["OnlyCases"]["note"]


def test_auc_result_serialises() -> None:
    payload = evaluate_auc(_scores(), note="a note").to_json()
    assert payload["profile"] == "testprofile"
    assert payload["auc"] is not None
    assert payload["note"] == "a note"


# --------------------------------------------------------------------------- #
# specificity verdicts
# --------------------------------------------------------------------------- #


def _matrix(row: dict[str, float], profile: str = "crc") -> SpecificityMatrix:
    matrix = SpecificityMatrix()
    for condition, auc in row.items():
        matrix.add(
            profile,
            condition,
            AucResult(
                profile=profile, condition=condition, n_cases=100, n_controls=100,
                auc=auc, ci_low=auc - 0.05, ci_high=auc + 0.05,
            ),
        )
    return matrix


def test_condition_specific_profile_is_recognised() -> None:
    verdict = spread_verdict(_matrix({"CRC": 0.85, "IBD": 0.55, "T2D": 0.52}), "crc")
    assert verdict.startswith("condition-specific")
    assert "0.85" in verdict


def test_profile_scoring_higher_elsewhere_is_called_out() -> None:
    """The real CRC result: it separates IBD better than CRC."""
    verdict = spread_verdict(_matrix({"CRC": 0.62, "IBD": 0.76, "T2D": 0.52}), "crc")
    assert verdict.startswith("NOT condition-specific")
    assert "IBD" in verdict
    assert "general-disturbance" in verdict


def test_indistinguishable_scores_mean_general_dysbiosis() -> None:
    verdict = spread_verdict(_matrix({"CRC": 0.70, "IBD": 0.68, "T2D": 0.52}), "crc")
    assert verdict.startswith("NOT condition-specific")
    assert "general dysbiosis" in verdict


def test_condition_matching_is_case_insensitive() -> None:
    """cMD labels are upper case; profile names are lower."""
    verdict = spread_verdict(_matrix({"CRC": 0.85, "IBD": 0.55}), "crc")
    assert "UNMEASURED" not in verdict


def test_missing_own_condition_is_reported_unmeasured() -> None:
    verdict = spread_verdict(_matrix({"IBD": 0.55, "T2D": 0.52}, profile="mecfs"), "mecfs")
    assert "UNMEASURED" in verdict


def test_unmeasured_profile_scoring_high_elsewhere_is_flagged() -> None:
    verdict = spread_verdict(
        _matrix({"IBD": 0.80, "CRC": 0.75}, profile="longcovid"), "longcovid"
    )
    assert "UNMEASURED" in verdict
    assert "not specific" in verdict


def test_one_condition_is_not_enough_to_judge() -> None:
    verdict = spread_verdict(_matrix({"CRC": 0.85}), "crc")
    assert "not enough labelled conditions" in verdict


def test_matrix_serialises_with_its_interpretation() -> None:
    payload = _matrix({"CRC": 0.62, "IBD": 0.76}).to_json()
    assert "matrix" in payload
    assert "general gut disturbance" in payload["interpretation"]
