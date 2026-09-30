"""Shape analysis, competition annotation and the confounder ledger (spec 7.8, 12).

The analysis only reads a handful of attributes from each ``ProfileResult``,
so the tests build light stand-ins rather than scoring real profiles: what is
under test is the verdict logic, not the scoring.
"""

from __future__ import annotations

from types import SimpleNamespace

from openbiota.shape import (
    HIGH_AT,
    VERY_HIGH_AT,
    annotate_competition,
    confounder_ledger,
    empirical_specificity,
    shape_analysis,
)


def _result(name, pct, *, family=None, abstained=False, spec=None, empirical=None,
            triggered=()):
    vector = None
    if spec is not None or empirical is not None:
        vector = SimpleNamespace(
            empirical_specificity=empirical, literature_uniqueness_prior=spec
        )
    return SimpleNamespace(
        profile=SimpleNamespace(name=name, label=name.upper(), family=family or name),
        combined_percentile=None if abstained else pct,
        reportable=not abstained and pct is not None,
        abstention=SimpleNamespace(abstained=abstained, triggered=list(triggered)),
        vector=vector,
        competing=(),
        top_two_margin=None,
    )


# --------------------------------------------------------------------------- #
# shape verdicts
# --------------------------------------------------------------------------- #


def test_thresholds_match_the_report_bands():
    # 75 is the top of the typical band, 95 where "notably high" begins.
    assert HIGH_AT == 75.0
    assert VERY_HIGH_AT == 95.0


def test_quiet_when_nothing_is_above_the_typical_band():
    results = [_result(f"p{i}", 40 + i) for i in range(10)]
    shape = shape_analysis(results, anchor=None)
    assert shape.verdict == "quiet"
    assert shape.n_above_typical == 0
    assert "typical band" in shape.explanation


def test_specific_when_one_distinctive_profile_stands_alone():
    results = [_result(f"p{i}", 30 + i) for i in range(10)]
    results.append(_result("crc", 96, spec=0.9))
    shape = shape_analysis(results, anchor=None)
    assert shape.verdict == "specific"
    assert shape.n_above_typical == 1
    assert shape.n_notably_high == 1
    assert shape.high_scorers == ("crc",)


def test_diffuse_when_most_profiles_score_high_on_shared_features():
    results = [_result(f"p{i}", 80 + i, spec=0.2) for i in range(6)]
    results += [_result(f"q{i}", 40, spec=0.2) for i in range(4)]
    anchor = SimpleNamespace(score=-1.2, band="dysbiotic")
    shape = shape_analysis(results, anchor=anchor)
    assert shape.verdict == "diffuse"
    assert shape.n_above_typical == 6
    assert "dysbiosis index" in shape.explanation
    assert shape.dysbiosis_band == "dysbiotic"


def test_mixed_when_several_score_high_with_partial_specificity():
    results = [_result(f"p{i}", 80 + i, spec=0.6) for i in range(3)]
    results += [_result(f"q{i}", 40, spec=0.6) for i in range(10)]
    shape = shape_analysis(results, anchor=None)
    assert shape.verdict == "mixed"
    assert "3 of 13" in shape.headline


def test_unscored_when_every_profile_abstained():
    results = [_result("a", None, abstained=True), _result("b", None, abstained=True)]
    shape = shape_analysis(results, anchor=None)
    assert shape.verdict == "unscored"
    assert shape.n_abstained == 2
    assert shape.n_scored == 0


def test_specificity_basis_prefers_empirical_over_prior():
    results = [_result("crc", 90, spec=0.3, empirical=0.8), _result("ibd", 88, spec=0.3)]
    shape = shape_analysis(results, anchor=None)
    assert "1 of 2 from challenge-cohort AUROCs" in shape.specificity_basis
    assert shape.mean_specificity_of_high == (0.8 + 0.3) / 2


def test_to_json_uses_band_vocabulary_not_raw_thresholds():
    shape = shape_analysis([_result("a", 50)], anchor=None)
    payload = shape.to_json()
    assert "n_above_typical_band" in payload
    assert "n_notably_high" in payload
    assert "n_above_60th" not in payload


# --------------------------------------------------------------------------- #
# competition
# --------------------------------------------------------------------------- #


def test_competition_lists_other_families_and_top_two_margin():
    crohns = _result("crohns", 90, family="ibd")
    ibd = _result("ibd", 88, family="ibd")
    crc = _result("crc", 70)
    quiet = _result("t2d", 20)
    annotate_competition([crohns, ibd, crc, quiet])
    # Same-family siblings are not "competitors".
    assert all(name != "ibd" for name, _ in crohns.competing)
    assert crohns.competing[0] == ("crc", 70)
    assert crohns.top_two_margin == 2
    assert quiet.top_two_margin == 2


def test_competition_with_a_single_scored_profile_has_no_margin():
    only = _result("crc", 80)
    annotate_competition([only, _result("x", None, abstained=True)])
    assert only.competing == ()
    assert only.top_two_margin is None


# --------------------------------------------------------------------------- #
# empirical specificity from the validation matrix
# --------------------------------------------------------------------------- #


def test_empirical_specificity_penalises_off_target_discrimination():
    validation = {
        "specificity_matrix": {
            "matrix": {
                "crc": {"CRC": {"auc": 0.82}, "IBD": {"auc": 0.75}, "T2D": {"auc": 0.55}},
            }
        }
    }
    value, basis = empirical_specificity("crc", validation)
    # worst off-target 0.75 -> excess 0.5 -> specificity 0.5
    assert value == 0.5
    assert "0.75" in basis and "own condition 0.82" in basis


def test_empirical_specificity_absent_when_profile_not_in_matrix():
    value, basis = empirical_specificity("ms", {"specificity_matrix": {"matrix": {"crc": {}}}})
    assert value is None
    assert "not yet scored" in basis
    assert empirical_specificity("crc", None) == (None, "no challenge-cohort matrix available")


# --------------------------------------------------------------------------- #
# confounder ledger
# --------------------------------------------------------------------------- #


def test_ledger_separates_recorded_missing_and_abstaining():
    metadata = {"country": "USA", "age": 42, "sex": "male", "metformin": None, "ppi": ""}
    results = [
        _result("t2d", None, abstained=True, triggered=["metformin status unknown"]),
        _result("crc", 60),
    ]
    ledger = confounder_ledger(metadata, results)
    recorded = {k for k, _, _ in ledger.recorded}
    missing = {k for k, _, _ in ledger.missing}
    assert {"country", "age", "sex"} <= recorded
    assert {"metformin", "ppi", "antibiotics"} <= missing
    assert ledger.abstentions == (("T2D", "metformin status unknown"),)
    payload = ledger.to_json()
    assert payload["profiles_abstained"][0]["profile"] == "T2D"
    assert any(m["field"] == "metformin" for m in payload["missing"])
