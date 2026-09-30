"""The scoring engine: robust values, aggregation, percentiles, abstention.

Arithmetic is checked against values computed by hand, and the two bugs this
engine actually shipped with are pinned as regression tests:

* a per-sample CLR pseudocount, which made the transform depth-dependent
* treating an undetected taxon as a number rather than a category, which
  turned reference-frame noise into a similarity signal
"""

from __future__ import annotations

import math
import statistics

import pytest

from openbiota.profiles import parse_profile
from openbiota.refcohort import CLR_FLOOR_PERCENT, clr, prevalent_frame
from openbiota.scoring import (
    Z_CLIP,
    FeatureReference,
    ReferenceBundle,
    SampleMeasurements,
    _probit,
    aggregate_module,
    combined_null,
    ecological_metrics,
    percentile_against,
    score_feature,
    score_profile,
    shannon,
    summarise_reference,
)

# --------------------------------------------------------------------------- #
# the compositional transform
# --------------------------------------------------------------------------- #


def test_clr_sums_to_zero() -> None:
    """A centred log-ratio is centred: its components sum to zero."""
    out = clr([10.0, 20.0, 30.0, 40.0])
    assert sum(out) == pytest.approx(0.0, abs=1e-9)


def test_clr_is_scale_invariant() -> None:
    """Doubling every part leaves the CLR unchanged — that is the point."""
    a = clr([5.0, 10.0, 85.0])
    b = clr([10.0, 20.0, 170.0])
    assert a == pytest.approx(b)


def test_clr_zero_floor_does_not_depend_on_the_samples_minimum() -> None:
    """Regression: the floor must be a constant, not derived per sample.

    The textbook "half the smallest observed value" rule makes the transform
    depth-dependent: a deeper sample detects rarer taxa, gets a lower floor,
    and maps its zeros further down than a shallow one. For a taxon absent in
    most samples, that difference *is* the reported signal.

    The invariant is that a zero component's pre-centring value equals
    ``log(CLR_FLOOR_PERCENT)`` whatever else the sample contains. Both
    compositions below have the same components; only the smallest detected
    value differs, by three orders of magnitude.
    """
    for parts in ([60.0, 39.0, 1.0, 0.0], [60.0, 39.999, 0.001, 0.0]):
        out = clr(parts)
        # Undo the centring: the mean of the logs, floor included.
        logs = [math.log(v if v > 0 else CLR_FLOOR_PERCENT) for v in parts]
        recovered = out[3] + statistics.fmean(logs)
        assert recovered == pytest.approx(math.log(CLR_FLOOR_PERCENT))


def test_clr_pseudocount_default_is_not_derived_from_the_data() -> None:
    """Passing the default explicitly must change nothing."""
    parts = [70.0, 20.0, 10.0, 0.0]
    assert clr(parts) == pytest.approx(clr(parts, pseudocount=CLR_FLOOR_PERCENT))


def test_clr_floor_matches_the_documented_constant() -> None:
    out = clr([100.0, 0.0])
    # log(floor) - mean(log) for a two-part composition
    expected = math.log(CLR_FLOOR_PERCENT) - statistics.fmean(
        [math.log(100.0), math.log(CLR_FLOOR_PERCENT)]
    )
    assert out[1] == pytest.approx(expected)


def test_prevalent_frame_needs_enough_members() -> None:
    """A frame of three taxa is noise; fall back to the whole composition."""
    abundance = [[1.0, 1.0], [0.0, 0.0], [2.0, 2.0]]
    assert prevalent_frame(abundance) == []


def test_prevalent_frame_selects_prevalent_taxa() -> None:
    # 30 taxa, first 25 present in both samples, last 5 in neither
    abundance = [[1.0, 1.0] for _ in range(25)] + [[0.0, 0.0] for _ in range(5)]
    frame = prevalent_frame(abundance)
    assert frame == list(range(25))


def test_clr_frame_restricts_the_denominator() -> None:
    values = [10.0, 10.0, 0.0, 0.0]
    unframed = clr(values)
    framed = clr(values, frame=[0, 1])
    assert framed[0] == pytest.approx(0.0)  # equal to the frame's geometric mean
    assert framed != pytest.approx(unframed)


# --------------------------------------------------------------------------- #
# probit
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("p", "expected"),
    [(0.5, 0.0), (0.975, 1.959964), (0.025, -1.959964), (0.8413447, 1.0), (0.99, 2.326348)],
)
def test_probit_inverts_the_normal_cdf(p: float, expected: float) -> None:
    assert _probit(p) == pytest.approx(expected, abs=1e-4)


def test_probit_is_zero_outside_the_unit_interval() -> None:
    assert _probit(0.0) == 0.0
    assert _probit(1.0) == 0.0


# --------------------------------------------------------------------------- #
# per-feature scoring
# --------------------------------------------------------------------------- #


def _feature(direction: str = "decreased", **factors):
    data = {
        "name": "testprofile",
        "version": "0.1.0",
        "label": "T",
        "status": "research_only",
        "summary": "s",
        "citation": "c",
        "modules": {
            "taxonomic": {
                "weight_in_combined": 1.0,
                "features": [
                    {
                        "name": "Test_species",
                        "level": "species",
                        "direction": direction,
                        **{"w": 1.0, "q": 1.0, "s": 1.0, "c": 1.0, **factors},
                    }
                ],
            }
        },
    }
    return parse_profile(data).features(module="taxonomic")[0]


def test_mad_z_matches_the_spec_formula() -> None:
    """z = (value - median) / (1.4826 * MAD)."""
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference("Test_species", values, "test", prevalence=1.0)
    median = statistics.median(values)
    mad = statistics.median([abs(v - median) for v in values])
    assert reference.robust_z(80.0) == pytest.approx((80.0 - median) / (1.4826 * mad))


def test_v_is_direction_signed_and_clipped() -> None:
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference("Test_species", values, "test", prevalence=1.0)

    # A sample far below the median, for a feature reported as decreased,
    # is strongly concordant: v -> +1.
    low = score_feature(
        _feature("decreased"), transformed_value=-500.0, raw_value=0.0, reference=reference
    )
    assert low.v == pytest.approx(1.0)

    # The same value for a feature reported as increased is discordant.
    high = score_feature(
        _feature("increased"), transformed_value=-500.0, raw_value=0.0, reference=reference
    )
    assert high.v == pytest.approx(-1.0)


def test_v_never_exceeds_one_in_magnitude() -> None:
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference("Test_species", values, "test", prevalence=1.0)
    for value in (-1e9, -5.0, 50.0, 1e9):
        score = score_feature(
            _feature(), transformed_value=value, raw_value=value, reference=reference
        )
        assert score.v is not None
        assert -1.0 <= score.v <= 1.0


def test_unmeasured_feature_is_not_scored() -> None:
    score = score_feature(_feature(), transformed_value=None, raw_value=None, reference=None)
    assert not score.measured
    assert score.v is None
    assert score.unmeasured_reason


def test_missing_reference_is_reported_not_guessed() -> None:
    score = score_feature(_feature(), transformed_value=1.0, raw_value=1.0, reference=None)
    assert not score.measured
    assert "no reference distribution" in score.unmeasured_reason


# --------------------------------------------------------------------------- #
# zero inflation: absence is a category, not a number
# --------------------------------------------------------------------------- #


def _zero_inflated_reference(prevalence: float = 0.05, n: int = 1000) -> FeatureReference:
    """A rare taxon: present in 5% of the reference, absent in the rest.

    Absent samples get slightly different CLR values, because the CLR carries
    the geometric-mean denominator and each sample detected a different number
    of other taxa. That spread is the noise the prevalence-aware path removes.
    """
    import random

    rng = random.Random(7)
    n_present = int(n * prevalence)
    values: list[float] = []
    detected: list[bool] = []
    for index in range(n):
        if index < n_present:
            values.append(rng.uniform(2.0, 6.0))
            detected.append(True)
        else:
            values.append(-9.0 + rng.uniform(-1.5, 1.5))  # floor plus frame noise
            detected.append(False)
    return summarise_reference(
        "Test_species", values, "test", prevalence=prevalence, detected=detected
    )


def test_zero_inflated_reference_is_recognised() -> None:
    reference = _zero_inflated_reference()
    assert reference.zero_inflated
    assert reference.n_absent == 950


def test_absent_sample_lands_mid_tie_block() -> None:
    """Regression: an absent taxon must not score as if it were elevated.

    With the taxon absent in 95% of the reference, a sample that also lacks it
    sits at the mid-rank of that block — about the 47th percentile — not at
    the 75th, which is what comparing raw CLR values produced.
    """
    reference = _zero_inflated_reference(prevalence=0.05)
    percentile = reference.prevalence_aware_percentile(-9.3, detected=False)
    assert percentile == pytest.approx(47.5, abs=0.1)


def test_absent_sample_contributes_almost_nothing() -> None:
    reference = _zero_inflated_reference(prevalence=0.05)
    score = score_feature(
        _feature("increased"), transformed_value=-9.3, raw_value=0.0,
        reference=reference, detected=False,
    )
    assert score.v is not None
    assert abs(score.v) < 0.05, "absence must not read as evidence of enrichment"


def test_detected_sample_ranks_above_every_absent_one() -> None:
    reference = _zero_inflated_reference(prevalence=0.05)
    percentile = reference.prevalence_aware_percentile(4.0, detected=True)
    assert percentile > 95.0


def test_detecting_a_rare_taxon_is_strong_evidence() -> None:
    reference = _zero_inflated_reference(prevalence=0.05)
    score = score_feature(
        _feature("increased"), transformed_value=5.5, raw_value=0.4,
        reference=reference, detected=True,
    )
    assert score.v is not None
    assert score.v > 0.4


def test_prevalent_feature_still_uses_the_mad_path() -> None:
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference(
        "Test_species", values, "test", prevalence=1.0, detected=[True] * 100
    )
    assert not reference.zero_inflated
    median = statistics.median(values)
    mad = statistics.median([abs(v - median) for v in values])
    assert reference.robust_z(80.0, detected=True) == pytest.approx(
        (80.0 - median) / (1.4826 * mad)
    )


# --------------------------------------------------------------------------- #
# aggregation
# --------------------------------------------------------------------------- #


def test_module_aggregation_is_the_weighted_mean() -> None:
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference("Test_species", values, "test", prevalence=1.0)
    a = score_feature(_feature("increased", w=1.0, q=1.0, s=1.0, c=1.0),
                      transformed_value=1e9, raw_value=1.0, reference=reference)
    b = score_feature(_feature("increased", w=0.5, q=1.0, s=1.0, c=1.0),
                      transformed_value=-1e9, raw_value=0.0, reference=reference)
    # (1.0*(+1) + 0.5*(-1)) / (1.0 + 0.5)
    assert aggregate_module([a, b]) == pytest.approx((1.0 - 0.5) / 1.5)


def test_aggregation_ignores_unmeasured_features() -> None:
    values = [float(v) for v in range(1, 101)]
    reference = summarise_reference("Test_species", values, "test", prevalence=1.0)
    measured = score_feature(_feature("increased"), transformed_value=1e9,
                             raw_value=1.0, reference=reference)
    missing = score_feature(_feature("increased"), transformed_value=None,
                            raw_value=None, reference=reference)
    assert aggregate_module([measured, missing]) == pytest.approx(1.0)


def test_aggregation_of_nothing_is_none() -> None:
    assert aggregate_module([]) is None


# --------------------------------------------------------------------------- #
# percentiles
# --------------------------------------------------------------------------- #


def test_percentile_against_null_uses_ranking() -> None:
    null = [0.0] * 50 + [1.0] * 50
    assert percentile_against(null, -1.0) == pytest.approx(0.0)
    assert percentile_against(null, 2.0) == pytest.approx(100.0)
    assert percentile_against(null, 0.0) == pytest.approx(25.0)  # mid-rank of the tie


def test_percentile_of_empty_null_is_none() -> None:
    assert percentile_against([], 0.5) is None
    assert percentile_against([1.0], None) is None


def test_combined_null_respects_module_weights() -> None:
    nulls = {"taxonomic": [1.0] * 100, "functional": [0.0] * 100}
    weights = {"taxonomic": 0.5, "functional": 0.35}
    null = combined_null(nulls, weights, draws=200)
    expected = (0.5 * 1.0 + 0.35 * 0.0) / 0.85
    assert all(v == pytest.approx(expected) for v in null)


def test_combined_null_is_empty_without_usable_modules() -> None:
    assert combined_null({}, {}) == []


# --------------------------------------------------------------------------- #
# ecological metrics
# --------------------------------------------------------------------------- #


def test_shannon_of_a_single_taxon_is_zero() -> None:
    assert shannon([100.0]) == pytest.approx(0.0)


def test_shannon_is_maximal_when_even() -> None:
    even = shannon([25.0] * 4)
    skewed = shannon([97.0, 1.0, 1.0, 1.0])
    assert even == pytest.approx(math.log(4))
    assert skewed < even


def test_ecological_metrics_reports_three_values() -> None:
    metrics = ecological_metrics({"a": 50.0, "b": 30.0, "c": 20.0})
    assert set(metrics) == {"shannon_diversity", "richness", "evenness"}
    assert metrics["richness"] == 3.0
    assert 0.0 < metrics["evenness"] <= 1.0


# --------------------------------------------------------------------------- #
# abstention
# --------------------------------------------------------------------------- #


def _profile_with_abstention():
    return parse_profile(
        {
            "name": "testprofile",
            "version": "0.1.0",
            "label": "T",
            "status": "research_only",
            "summary": "s",
            "citation": "c",
            "modules": {
                "taxonomic": {
                    "weight_in_combined": 1.0,
                    "features": [
                        {"name": "Test_species", "level": "species", "direction": "decreased"}
                    ],
                }
            },
            "abstain_if": [
                {"recent_antibiotics_days_lt": 90},
                {"usable_nonhost_reads_lt": 500000},
                {"missing_metadata": ["onset_date"]},
            ],
        }
    )


def _bundle(n: int = 500) -> ReferenceBundle:
    values = [float(v) for v in range(n)]
    bundle = ReferenceBundle(
        taxonomic={"Test_species": summarise_reference("Test_species", values, "t", prevalence=1.0)},
        taxonomic_source="test",
        taxonomic_n=n,
    )
    bundle.match_info = {"matched": True, "reason": "matched"}
    return bundle


def test_abstains_on_recent_antibiotics() -> None:
    result = score_profile(
        profile=_profile_with_abstention(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            antibiotics_days_ago=14,
            usable_nonhost_reads=10_000_000,
            metadata={"onset_date": "2026-01-01"},
        ),
        references=_bundle(),
    )
    assert result.abstention.abstained
    assert any("antibiotic" in r.lower() for r in result.abstention.triggered)
    assert result.abstention.remedies


def test_abstains_on_insufficient_depth() -> None:
    result = score_profile(
        profile=_profile_with_abstention(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=100_000,
            metadata={"onset_date": "2026-01-01"},
        ),
        references=_bundle(),
    )
    assert result.abstention.abstained
    assert any("500,000" in r for r in result.abstention.triggered)


def test_abstains_on_missing_required_metadata() -> None:
    result = score_profile(
        profile=_profile_with_abstention(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=10_000_000,
            metadata={},
        ),
        references=_bundle(),
    )
    assert result.abstention.abstained
    assert any("onset_date" in r for r in result.abstention.triggered)


def test_does_not_abstain_when_everything_is_supplied() -> None:
    result = score_profile(
        profile=_profile_with_abstention(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            antibiotics_days_ago=400,
            usable_nonhost_reads=10_000_000,
            metadata={"onset_date": "2026-01-01"},
        ),
        references=_bundle(),
    )
    assert not result.abstention.abstained
    assert result.combined_percentile is not None


def test_abstains_when_the_reference_cohort_is_too_small() -> None:
    result = score_profile(
        profile=_profile_with_abstention(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            antibiotics_days_ago=400,
            usable_nonhost_reads=10_000_000,
            metadata={"onset_date": "2026-01-01"},
        ),
        references=_bundle(n=40),
    )
    assert result.abstention.abstained
    assert any("reference cohort" in r for r in result.abstention.triggered)


# --------------------------------------------------------------------------- #
# confidence
# --------------------------------------------------------------------------- #


def test_missing_metadata_downgrades_confidence() -> None:
    profile = parse_profile(
        {
            "name": "testprofile", "version": "0.1.0", "label": "T",
            "status": "research_only", "summary": "s", "citation": "c",
            "modules": {
                "taxonomic": {
                    "weight_in_combined": 1.0,
                    "features": [
                        {"name": "Test_species", "level": "species", "direction": "decreased"}
                    ],
                }
            },
        }
    )
    bare = score_profile(
        profile=profile,
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=10_000_000,
        ),
        references=_bundle(),
    )
    assert bare.confidence.grade != "HIGH"
    assert any("missing metadata" in r for r in bare.confidence.reasons)


def test_full_metadata_and_matched_reference_gives_high_confidence() -> None:
    profile = parse_profile(
        {
            "name": "testprofile", "version": "0.1.0", "label": "T",
            "status": "research_only", "summary": "s", "citation": "c",
            "references": {"primary": {"source": "PRJNA1", "note": "exists"}},
            "modules": {
                "taxonomic": {
                    "weight_in_combined": 1.0,
                    "features": [
                        {"name": "Test_species", "level": "species", "direction": "decreased"}
                    ],
                }
            },
        }
    )
    result = score_profile(
        profile=profile,
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=10_000_000,
            metadata={
                "antibiotics": True, "stool_form": "4", "onset_date": "2026-01-01",
                "sampling_date": "2026-02-01", "medications": "none",
            },
        ),
        references=_bundle(),
    )
    assert result.confidence.grade == "HIGH", result.confidence.reasons


def test_z_clip_constant_is_three() -> None:
    assert Z_CLIP == 3.0
