"""The community-health adjustment on top of the species age model.

The species model reads composition only and returns close to its adult
training mean for nearly every sample — including four children, whom it
put in their thirties. These tests pin the arithmetic that moves it, the
bound that stops it running away, and the honesty requirements: the
unadjusted estimate survives, and nothing here claims to be calibrated.
"""

from __future__ import annotations

from typing import Any

import pytest

from openbiota.age import (
    GMWI2_YEARS_PER_POINT,
    MIN_REPORTABLE_AGE,
    AgeResult,
    adjust_for_community_health,
)


def _scored(point: float, half80: float = 19.6) -> AgeResult:
    """A scored result with a symmetric 80% interval and nothing else set."""
    return AgeResult(
        status="scored", reason=None,
        predicted_chronological_age_years=point,
        prediction_interval_80=(point - half80, point + half80),
        prediction_interval_95=(point - 30.0, point + 30.0),
        actual_chronological_age_years=None, age_residual_years=None,
        matched_stratum_mae_years=None, matched_stratum_n=None, matched_stratum=None,
        calibration_population="", ood_status="in_distribution", ood_detail={},
        model_bundle_id="test",
        contributions_increasing_this_prediction=[],
        contributions_decreasing_this_prediction=[],
        attribution_method_and_background_digest="",
        species_matched=0, species_in_sample=0, nearest_reference={},
        interval_spans_multiple_life_stages=False,
    )


def test_disease_patterns_move_a_year_each_with_a_dead_band() -> None:
    """Above the 90th adds a year, below the 79th takes one, between: nothing."""
    out = adjust_for_community_health(
        _scored(40.0), gmwi2=None,
        profile_percentiles=[95.0, 92.0, 90.0, 85.0, 80.0, 79.0, 50.0, 10.0],
    )
    d = out.health_adjustment
    assert d is not None
    assert d["profiles_adverse_at_or_above_90th"] == 3
    assert d["profiles_equivocal_79th_to_89th"] == 3
    assert d["profiles_typical_below_79th"] == 2
    assert d["disease_years"] == pytest.approx(1.0)
    assert out.health_adjusted_age_years == pytest.approx(41.0)


def test_boundaries_land_on_the_stated_side() -> None:
    """Exactly 90 counts against; exactly 79 is equivocal, not favourable."""
    at_90 = adjust_for_community_health(
        _scored(40.0), gmwi2=None, profile_percentiles=[90.0]
    ).health_adjustment
    at_79 = adjust_for_community_health(
        _scored(40.0), gmwi2=None, profile_percentiles=[79.0]
    ).health_adjustment
    just_under_79 = adjust_for_community_health(
        _scored(40.0), gmwi2=None, profile_percentiles=[78.9]
    ).health_adjustment
    assert at_90 is not None and at_79 is not None and just_under_79 is not None
    assert at_90["disease_years"] == pytest.approx(1.0)
    assert at_79["disease_years"] == pytest.approx(0.0)
    assert just_under_79["disease_years"] == pytest.approx(-1.0)


def test_a_healthier_community_reads_younger() -> None:
    """The GMWI2 term's sign is the whole point of including it."""
    healthy = adjust_for_community_health(_scored(40.0), gmwi2=2.0)
    dysbiotic = adjust_for_community_health(_scored(40.0), gmwi2=-2.0)
    assert healthy.health_adjusted_age_years is not None
    assert dysbiotic.health_adjusted_age_years is not None
    assert healthy.health_adjusted_age_years < 40.0 < dysbiotic.health_adjusted_age_years
    # Symmetric and at the stated rate.
    assert healthy.health_adjusted_age_years == pytest.approx(
        40.0 - 2.0 * GMWI2_YEARS_PER_POINT
    )
    assert dysbiotic.health_adjusted_age_years == pytest.approx(
        40.0 + 2.0 * GMWI2_YEARS_PER_POINT
    )


def test_the_answer_cannot_leave_the_models_own_90_percent_interval() -> None:
    """This is an adjustment, not a second model.

    A large pile of favourable readings must not walk the estimate to an
    age the model never considered plausible.
    """
    out = adjust_for_community_health(
        _scored(40.0), gmwi2=5.0, profile_percentiles=[5.0] * 40
    )
    d = out.health_adjustment
    assert d is not None and out.prediction_interval_90 is not None
    lo90, hi90 = out.prediction_interval_90
    assert d["unclamped_years"] < lo90, "the test must actually exercise the floor"
    assert out.health_adjusted_age_years == pytest.approx(lo90)
    assert d["clamped"] and d["clamped_to"] == "floor"
    # And the band really is the 90% one: tighter than 95, wider than 80.
    pi80, pi95 = out.prediction_interval_80, out.prediction_interval_95
    assert pi80 is not None and pi95 is not None
    assert pi80[0] > lo90 > pi95[0]
    assert pi80[1] < hi90 < pi95[1]


def test_the_floor_never_goes_below_a_liveable_age() -> None:
    out = adjust_for_community_health(
        _scored(12.0), gmwi2=5.0, profile_percentiles=[5.0] * 40
    )
    assert out.health_adjusted_age_years is not None
    assert out.health_adjusted_age_years >= MIN_REPORTABLE_AGE


def test_the_unadjusted_estimate_survives_and_the_working_is_reported() -> None:
    """A reader must be able to see what the model said on its own."""
    out = adjust_for_community_health(
        _scored(38.7), gmwi2=-1.354, profile_percentiles=[95.0] * 18 + [50.0] * 11
    )
    assert out.predicted_chronological_age_years == pytest.approx(38.7)
    assert out.reported_age_years == out.health_adjusted_age_years
    d = out.health_adjustment
    assert d is not None
    # Every term of the sum is recoverable from the record.
    assert d["model_point_years"] + d["disease_years"] + d["gmwi2_years"] == pytest.approx(
        d["unclamped_years"], abs=0.05
    )
    assert "not a validated calibration" in d["tuning_basis"]


def test_an_unscored_result_is_left_exactly_alone() -> None:
    """No model estimate, nothing to adjust — and no invented one."""
    bad = AgeResult(
        status="not_computable", reason="too few species",
        predicted_chronological_age_years=None,
        prediction_interval_80=None, prediction_interval_95=None,
        actual_chronological_age_years=None, age_residual_years=None,
        matched_stratum_mae_years=None, matched_stratum_n=None, matched_stratum=None,
        calibration_population="", ood_status="in_distribution", ood_detail={},
        model_bundle_id="test",
        contributions_increasing_this_prediction=[],
        contributions_decreasing_this_prediction=[],
        attribution_method_and_background_digest="",
        species_matched=0, species_in_sample=0, nearest_reference={},
        interval_spans_multiple_life_stages=False,
    )
    out = adjust_for_community_health(bad, gmwi2=2.0, profile_percentiles=[95.0])
    assert out is bad
    assert out.health_adjusted_age_years is None
    assert out.reported_age_years is None


def test_the_cli_reads_the_percentile_field_that_actually_exists() -> None:
    """Guard against the silent-zero failure mode.

    The first wiring of this asked each scored pattern for `.percentile`,
    which `ProfileResult` does not have. `getattr(r, "percentile", None)`
    returned None for every one, the disease term came out as exactly zero,
    and the reports looked plausible while carrying none of the signal that
    makes the estimate work. A wrong attribute name must break a test, not a
    report.
    """
    import ast
    import inspect

    from openbiota import cli
    from openbiota.scoring import ProfileResult

    assert "combined_percentile" in ProfileResult.__dataclass_fields__
    assert isinstance(
        inspect.getattr_static(ProfileResult, "reportable"), property
    )

    tree = ast.parse(inspect.getsource(cli))
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "adjust_for_community_health"
    ]
    assert len(calls) == 1, "one wiring site, so one place to get this wrong"
    kwargs = {k.arg: k for k in calls[0].keywords}
    assert "profile_percentiles" in kwargs and "gmwi2" in kwargs
    used = {
        n.attr
        for n in ast.walk(kwargs["profile_percentiles"].value)
        if isinstance(n, ast.Attribute)
    }
    assert "combined_percentile" in used, (
        f"the CLI reads {used}, none of which is the field that holds the "
        "percentile"
    )
    # And no bare getattr default, which is what hid the mistake.
    assert not [
        n for n in ast.walk(kwargs["profile_percentiles"].value)
        if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "getattr"
    ], "a getattr default here turns a typo into a silent zero"


def test_the_adjustment_is_reported_in_the_results_json() -> None:
    out = adjust_for_community_health(_scored(40.0), gmwi2=1.0, profile_percentiles=[95.0])
    doc: dict[str, Any] = out.to_json()
    assert doc["health_adjusted_age_years"] == pytest.approx(39.0)
    assert doc["prediction_interval_90"] is not None
    assert doc["health_adjustment"]["gmwi2_years"] == pytest.approx(-GMWI2_YEARS_PER_POINT)
    # The unadjusted figure stays in the record.
    assert doc["predicted_chronological_age_years"] == pytest.approx(40.0)
