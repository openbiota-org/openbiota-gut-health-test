"""Chronic urticaria (chronic hives) increment — build spec v4.1 acceptance tests.

Chronic urticaria is chronic hives, and the report is required to say so
wherever it names the condition. Beyond the wording, each test here pins one
promise:

* the two registered tasks are added by union and adding them twice is
  idempotent — no unrelated profile arrives to make a count look right;
* the feature manifest is exactly what the source table says: fifteen unique
  species, five increases and ten decreases, with exactly eight surviving
  multiple-testing correction;
* a published association coefficient is provenance, never a scoring weight;
* no genus figure may fill a species slot, and no slot overlaps another;
* nothing that was not measured is scored as zero — it lowers coverage and
  widens the bounds instead;
* the numbers themselves match the specification's fixtures to 1e-10;
* no card promises a cure, an eradication or a transplant, and no
  microbiome result selects a treatment.

Test IDs from the specification's acceptance table appear in the docstrings.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

from openbiota import cuindex as cu
from openbiota.interventions import load_registry
from openbiota.profiles import load_profile_set

REPO = Path(__file__).resolve().parents[1]
PROFILES = REPO / "profiles"
INTERVENTIONS = REPO / "interventions"

URTICARIA_PROFILES = ("csu", "symptomatic_dermographism")


@pytest.fixture(scope="module")
def profile_set():
    return load_profile_set(PROFILES)


@pytest.fixture(scope="module")
def registry():
    return load_registry(INTERVENTIONS)


@pytest.fixture(scope="module")
def controls() -> list[dict[str, float]]:
    """Twenty-five synthetic control rows, each a distinct participant.

    Synthetic numerical fixtures only. Nothing here is reference data, and
    production enforces distinct participant identity before fitting.
    """
    return [
        {"a": 0.001 * (1.0 + 0.05 * i), "b": 0.01 * (1.0 + 0.05 * i)}
        for i in range(25)
    ]


# --------------------------------------------------------------------------- #
# registration (CU01, CU02, CU56)
# --------------------------------------------------------------------------- #


def test_both_urticaria_tasks_are_registered(profile_set) -> None:
    """CU01/CU02: two new tasks, and nothing unrelated added beside them."""
    names = {p.name for p in profile_set.profiles}
    assert set(URTICARIA_PROFILES) <= names
    # Both are hives profiles, and nothing else in the library claims to be.
    hives = {p.name for p in profile_set.profiles if "urticaria" in p.label.lower() or "hives" in p.label.lower()}
    assert hives == set(URTICARIA_PROFILES)


def test_loading_twice_is_idempotent() -> None:
    """CU01: migration runs twice without duplicating a task."""
    first = load_profile_set(PROFILES)
    second = load_profile_set(PROFILES)
    names_first = sorted(p.name for p in first.profiles)
    names_second = sorted(p.name for p in second.profiles)
    assert names_first == names_second
    assert len(names_first) == len(set(names_first))


def test_every_urticaria_mention_also_says_hives(profile_set) -> None:
    """The reader is never left to know that urticaria means hives.

    Any profile that names urticaria in its label, summary or caveats must
    say hives in the same place.
    """
    for name in URTICARIA_PROFILES:
        profile = next(p for p in profile_set.profiles if p.name == name)
        assert "hives" in profile.label.lower(), f"{name} label does not say hives"
        assert "urticaria" in profile.label.lower() or name == "symptomatic_dermographism"
        blocks = [profile.summary, *profile.caveats, profile.effect_size_note]
        for block in blocks:
            if "urticaria" in block.lower():
                assert "hives" in block.lower(), f"{name}: a block names urticaria without hives"
        # The summary explains what the condition actually is.
        assert "hives" in profile.summary.lower()


def test_unsupported_inducible_subtypes_are_named_and_not_scored(profile_set) -> None:
    """CU56: an unsupported subtype gets no score and no substitute.

    Cold, delayed-pressure, solar, heat, vibratory, cholinergic, contact and
    aquagenic urticaria have no stool-signature cohort. They must be named as
    not established, and neither hives profile may stand in for them.
    """
    names = {p.name for p in profile_set.profiles}
    for subtype in (
        "cold_urticaria", "delayed_pressure_urticaria", "solar_urticaria", "heat_urticaria",
        "vibratory_angioedema", "cholinergic_urticaria", "contact_urticaria", "aquagenic_urticaria",
        "acute_urticaria", "urticarial_vasculitis", "mastocytosis", "mcas",
    ):
        assert subtype not in names, f"{subtype} must not be a scored profile"

    text = " ".join(
        block.lower()
        for name in URTICARIA_PROFILES
        for profile in [next(p for p in profile_set.profiles if p.name == name)]
        for block in (profile.summary, *profile.caveats)
    )
    for subtype in ("cold", "solar", "cholinergic", "aquagenic", "vibratory", "contact"):
        assert subtype in text, f"{subtype} urticaria is not declared unsupported anywhere"
    assert "evidence not established" in text or "no stool signature" in text or "no stool-signature" in text


# --------------------------------------------------------------------------- #
# feature manifest (CU04, CU05, CU06, CU07, CU08, CU19, CU20, CU22)
# --------------------------------------------------------------------------- #


def test_csu_panel_is_fifteen_species_five_up_ten_down() -> None:
    """CU04: fifteen unique species, five positive and ten negative signs."""
    ids = cu.CSU_PANEL.ids
    assert len(ids) == 15 and len(set(ids)) == 15
    assert sum(1 for f in cu.CSU_PANEL.features if f.direction == 1) == 5
    assert sum(1 for f in cu.CSU_PANEL.features if f.direction == -1) == 10
    assert all(f.rank == "species" for f in cu.CSU_PANEL.features)


def test_strict_subset_is_exactly_the_eight_corrected_species() -> None:
    """CU05: filtering at q < 0.05 yields exactly the eight listed species."""
    strict = set(cu.CSU_STRICT_PANEL.ids)
    assert strict == {
        "Alistipes_onderdonkii", "Alistipes_putredinis", "Alistipes_shahii",
        "Bacteroides_intestinalis", "Coprococcus_catus", "Ruminococcus_obeum",
        "Veillonella_parvula", "Ruminococcus_gnavus",
    }
    assert len(strict) == 8
    # The mechanistically emphasised Enterobacteriaceae are not in it.
    assert "Escherichia_coli" not in strict and "Klebsiella_pneumoniae" not in strict


def test_source_statistics_reconcile_and_are_never_weights() -> None:
    """CU08: coefficients and q-values match the source table, as provenance."""
    by_id = {f.source_id: f for f in cu.CSU_PANEL.features}
    expected = {
        "Alistipes_onderdonkii": (-2.469669346, 0.015384911),
        "Coprococcus_catus": (-2.167381929, 0.002208388),
        "Ruminococcus_gnavus": (3.057857054, 0.002208388),
        "Bacteroides_stercoris": (1.827477087, 0.119496077),
        "Bilophila_wadsworthia": (-1.083777733, 0.119496077),
    }
    for name, (coef, q) in expected.items():
        assert abs(by_id[name].coef - coef) < 1e-9
        assert abs(by_id[name].q_value - q) < 1e-9
    # Provenance, not a weight: every feature's declared rule is equal weight,
    # and the record says so in machine-readable form.
    assert cu.CSU_PANEL.weight_rule == "equal"
    assert all(f.to_json()["coef_is_a_weight"] is False for f in cu.CSU_PANEL.features)


def test_profile_features_all_carry_equal_scoring_weight(profile_set) -> None:
    """The equal-weight rule reaches the scored profile, not just the manifest."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    panel = profile.modules["zhu_2024_full"]
    assert len(panel.features) == 15
    assert {f.w for f in panel.features} == {1.0}
    assert {f.factor for f in panel.features} == {1.0}


def test_bacteroides_stercoris_is_higher_and_bilophila_lower() -> None:
    """CU06/CU07: the corrected directions, with no invented conflict flag.

    B. stercoris rises in cases, so a higher abundance raises its
    contribution. B. wadsworthia falls in cases, so a higher abundance lowers
    it. The earlier draft had these the wrong way round.
    """
    by_id = {f.source_id: f for f in cu.CSU_PANEL.features}
    assert by_id["Bacteroides_stercoris"].direction == 1
    assert by_id["Bilophila_wadsworthia"].direction == -1

    ref = {"mu": -3.0, "scale": 1.0}
    low = cu.support(cu.z_value(1e-4, ref), 1)
    high = cu.support(cu.z_value(1e-2, ref), 1)
    assert high > low  # more B. stercoris, more alignment
    low_w = cu.support(cu.z_value(1e-4, ref), -1)
    high_w = cu.support(cu.z_value(1e-2, ref), -1)
    assert high_w < low_w  # more B. wadsworthia, less alignment


def test_no_genus_may_fill_a_ruminococcus_species_slot(profile_set) -> None:
    """CU19: three Ruminococcus species, three separate measurements."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    tax = [f for f in profile.features() if f.engine == "metaphlan"]
    assert tax, "the hives profile must have taxonomic features"
    assert all(f.level != "genus" for f in tax), "no genus feature belongs in a shotgun panel"
    # R. gnavus rises and R. obeum falls: a genus total would cancel them.
    by_name = {f.name: f for f in tax}
    assert by_name["Ruminococcus_gnavus"].direction == "increased"
    assert by_name["Ruminococcus_obeum"].direction == "decreased"


def test_reclassified_species_keeps_both_names(profile_set) -> None:
    """CU20/§4.2: the source name identifies the slot, the current name measures it."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    feature = next(f for f in profile.features() if f.name == "Ruminococcus_obeum")
    assert feature.measured_as == "Blautia_obeum"
    assert feature.key == "Blautia_obeum"
    assert feature.mapping_reason and "Blautia" in feature.mapping_reason
    assert "Ruminococcus obeum" in feature.display_name.replace("_", " ")
    assert "Blautia obeum" in feature.display_name.replace("_", " ")
    # No parent total is ever summed with a child.
    assert not any(f.name == "Ruminococcus" for f in profile.features())


def test_sd_panel_slots_do_not_overlap(profile_set) -> None:
    """CU22: the genus aggregate is non-overlapping and excludes R. bromii."""
    assert len(cu.SD_PANEL.features) == 2
    subdo, bromii = cu.SD_PANEL.features
    assert subdo.rank == "genus" and bromii.rank == "species"
    assert not bromii.source_id.startswith("Subdoligranulum")

    profile = next(p for p in profile_set.profiles if p.name == "symptomatic_dermographism")
    tax = [f for f in profile.features() if f.engine == "metaphlan"]
    assert {f.name for f in tax} == {"Subdoligranulum", "Ruminococcus_bromii"}
    # A genus feature is only tolerated because the module declares the transport.
    genus = next(f for f in tax if f.level == "genus")
    module = profile.modules[genus.module]
    assert module.assay_transport == "amplicon_taxon_to_shotgun"


def test_no_diversity_term_in_either_hives_profile(profile_set) -> None:
    """CU49: Shannon and evenness carry no weight in either index."""
    for name in URTICARIA_PROFILES:
        profile = next(p for p in profile_set.profiles if p.name == name)
        features = {f.name for f in profile.features()}
        assert not features & {"shannon", "shannon_diversity", "evenness", "richness"}


def test_sd_transport_warning_is_permanent(profile_set) -> None:
    """The 16S-to-shotgun label is on the module, so no figure escapes it."""
    profile = next(p for p in profile_set.profiles if p.name == "symptomatic_dermographism")
    module = profile.modules["liu_sd_2021_translated"]
    assert module.assay_transport == "amplicon_taxon_to_shotgun"
    assert profile.is_transported
    assert profile.evidence_maturity == "P1"


def test_predicted_histidine_pathway_is_not_a_histamine_module(profile_set) -> None:
    """CU76: an inferred pathway is not relabelled as measured production."""
    profile = next(p for p in profile_set.profiles if p.name == "symptomatic_dermographism")
    assert not any(f.name in ("histamine", "hdc") for f in profile.features())
    source = next(s for s in profile.unbound_sources if "histidine" in s.source.lower())
    assert source.evidence_type == "PATH_ABUND"
    note = source.note.lower()
    assert "prediction" in note or "predicted" in note
    assert "not a measurement" in note or "rather than a measurement" in note
    for claim in ("diamine-oxidase deficiency", "histamine intolerance", "mast-cell activation"):
        assert claim in note


def test_overlapping_lps_pathways_are_not_independent_terms(profile_set) -> None:
    """CU52: the six envelope pathways count once, and never in the index."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    source = next(s for s in profile.unbound_sources if "LPSSYN-PWY" in s.source)
    assert source.evidence_type == "PATH_ABUND"
    assert "not six independent mechanisms" in source.note
    # The carrier proxy that stands in for them is a mechanism lane at zero
    # weight, and its two features share one cluster.
    lane = profile.modules["gram_negative_envelope_lane"]
    assert lane.weight_in_combined == 0.0
    assert {f.cluster_key for f in lane.features} == {"gram_negative_envelope"}


def test_rgnavus_carries_no_canonical_lps_claim(profile_set) -> None:
    """CU53: R. gnavus is Gram-positive; its polysaccharide is not LPS."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    feature = next(f for f in profile.features() if f.name == "Ruminococcus_gnavus")
    note = feature.note.lower()
    assert "gram-positive" in note
    assert "not classical lipopolysaccharide" in note or "not classical lps" in note


def test_enterobacteriaceae_are_not_called_infection(profile_set) -> None:
    """CU57: an abundance is not an infection or an eradication target."""
    profile = next(p for p in profile_set.profiles if p.name == "csu")
    for name in ("Escherichia_coli", "Klebsiella_pneumoniae"):
        note = next(f for f in profile.features() if f.name == name).note.lower()
        assert "not an infection" in note or "is carriage, not infection" in note
        assert "eradication" in note or "antimicrobials" in note
    joined = " ".join(profile.caveats).lower()
    assert "not an infection" in joined and "eradication" in joined


# --------------------------------------------------------------------------- #
# the numerical contract (CU23–CU40, CU72, CU82)
# --------------------------------------------------------------------------- #


def test_specification_self_test_passes() -> None:
    """CU72: the specification's own numerical assertions."""
    cu.self_test()


def test_percent_and_fraction_agree_under_declared_units() -> None:
    """CU23: 1% and 0.01 are the same measurement once units are declared."""
    ref = {"mu": -2.0, "scale": 1.0}
    assert cu.z_value(1.0 / 100.0, ref) == cu.z_value(0.01, ref)


@pytest.mark.parametrize("bad", [True, False, -0.1, 1.1, float("nan"), float("inf"), "0.5", None])
def test_invalid_abundance_is_refused(bad) -> None:
    """CU24: nothing invalid becomes a finite patient score."""
    with pytest.raises(cu.CuIndexError):
        cu.fraction(bad)


def test_nineteen_controls_cannot_fit_a_reference(controls) -> None:
    """CU25: twenty distinct eligible controls is the engineering minimum."""
    with pytest.raises(cu.CuIndexError):
        cu.fit_reference(controls[:19], ["a"])
    assert cu.fit_reference(controls[:20], ["a"])


def test_constant_control_feature_uses_the_scale_floor() -> None:
    """CU27: a constant reference gives a finite transform, not a divide by zero."""
    ref = cu.fit_reference([{"a": 0.01}] * 20, ["a"])
    assert ref["a"]["scale"] == cu.SCALE_FLOOR
    assert math.isfinite(cu.z_value(0.01, ref["a"]))


def test_query_at_the_reference_median_contributes_fifty() -> None:
    """CU28: neutral signed deviation is fifty, which is not a health grade."""
    ref = cu.fit_reference([{"a": 0.01}] * 20, ["a"])
    assert cu.support(cu.z_value(0.01, ref["a"]), 1) == 50.0
    assert cu.support(cu.z_value(0.01, ref["a"]), -1) == 50.0


def test_two_standard_deviations_gives_the_specified_contributions() -> None:
    """CU29: the exact fixtures, both directions, to 1e-10."""
    assert abs(cu.support(2.0, 1) - 88.07970779778825) < 1e-10
    assert abs(cu.support(2.0, -1) - 11.92029220221176) < 1e-10


def test_partial_coverage_bounds_match_the_specified_fixture() -> None:
    """CU34: one of two features, point 80 over 70–90, gives 80 and 35–95."""
    result = cu.summarize({"a": (80.0, 70.0, 90.0)}, ["a", "b"])
    assert result["index"] == 80.0
    assert result["coverage"] == 0.5
    assert result["full_panel_bounds"] == [35.0, 95.0]


def test_no_usable_features_is_null_and_not_zero() -> None:
    """CU35: a null score with bounds 0–100, never a healthy-looking zero."""
    result = cu.summarize({}, ["a", "b"])
    assert result["index"] is None
    assert result["full_panel_bounds"] == [0.0, 100.0]

    empty = cu.score_panel(cu.CSU_PANEL, {}, None)
    assert empty.status == "reference_unavailable"
    assert empty.score is None
    assert empty.to_json()["score_0_100"] is None
    assert empty.to_json()["full_panel_bounds"] == [0.0, 100.0]


def test_duplicate_feature_ids_and_malformed_intervals_are_refused() -> None:
    """CU82: a compile failure, not silent acceptance."""
    ref = cu.fit_reference([{"a": 0.01, "b": 0.02}] * 20, ["a", "b"])
    with pytest.raises(cu.CuIndexError):
        cu.evaluate(
            [{"id": "a", "direction": 1}, {"id": "a", "direction": -1}],
            {"a": {"point": 0.01, "lower": 0.01, "upper": 0.01}},
            ref,
        )
    with pytest.raises(cu.CuIndexError):
        cu.evaluate(
            [{"id": "a", "direction": 1}],
            {"a": {"point": 0.05, "lower": 0.01, "upper": 0.02}},
            ref,
        )
    with pytest.raises(cu.CuIndexError):
        cu.summarize({}, [])


def test_midrank_percentile_fixture() -> None:
    """CU40: equal weights, ties at the query value, mid-rank 50."""
    assert cu.weighted_midrank(50, [20, 50, 50, 80], [1, 1, 1, 1]) == 50.0


def test_percentiles_are_null_without_held_out_scores() -> None:
    """CU39: no calibration set, no percentile, and never a percent sign."""
    ref = cu.fit_reference([{"Ruminococcus_bromii": 0.01}] * 20, ["Ruminococcus_bromii"])
    frozen = cu.FrozenReference(
        bundle_id="fixture", scope="discovery_controls_only", stats=ref,
        prevalence={"Ruminococcus_bromii": 1.0}, n_participants=20,
    )
    result = cu.score_panel(cu.SD_PANEL, {"Ruminococcus_bromii": 0.01}, frozen)
    payload = result.to_json()
    assert payload["control_percentile"] is None and payload["case_percentile"] is None
    assert payload["score_is_not_a_percentage"] is True
    assert "reconstruction" in payload["percentile_note"]


# --------------------------------------------------------------------------- #
# missingness and censoring against the real catalogue (CU31–CU33, CU36, CU38)
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def frozen_csu():
    stats = {
        f.source_id: {"mu": -3.0, "scale": 1.0, "n_controls": 30.0}
        for f in cu.CSU_PANEL.features
    }
    return cu.FrozenReference(
        bundle_id="fixture:n=30",
        scope="compatible_independent_controls",
        stats=stats,
        median_abundance={f.source_id: 1e-3 for f in cu.CSU_PANEL.features},
        prevalence={f.source_id: 0.9 for f in cu.CSU_PANEL.features},
        n_participants=30,
    )


def _abundances(**overrides: float) -> dict[str, float]:
    """Every panel key detected at the reference median, then overridden."""
    out = {f.key: 1e-3 for f in cu.CSU_PANEL.features}
    out.update(overrides)
    return out


def test_a_taxon_missing_from_the_sample_table_is_not_zero_imputed(frozen_csu) -> None:
    """CU31: a species the profiler did not list is not scored as zero.

    The profiler reports only what it detected. The reference catalogue does
    contain this organism, so its absence here is a nondetection with no
    validated limit: the point is dropped and the bounds are kept.
    """
    abundances = _abundances()
    del abundances["Coprococcus_catus"]
    result = cu.score_panel(cu.CSU_PANEL, abundances, frozen_csu)
    row = next(r for r in result.rows if r.feature.source_id == "Coprococcus_catus")
    assert row.state == "below_reporting_limit"
    assert row.abundance is None and row.contribution is None
    assert "not a measured absence" in row.note
    assert result.usable_count == 14 and result.coverage == pytest.approx(14 / 15)
    lo, hi = result.full_panel_bounds
    assert lo < result.score < hi  # the missing slot widened the bounds


def test_a_taxon_with_no_catalogue_slot_is_not_zero_imputed() -> None:
    """CU31/CU21: no slot in the catalogue means the absence was never measured."""
    stats = {"Ruminococcus_bromii": {"mu": -3.0, "scale": 1.0, "n_controls": 30.0}}
    frozen = cu.FrozenReference(
        bundle_id="fixture", scope="compatible_independent_controls", stats=stats,
        median_abundance={"Ruminococcus_bromii": 1e-3},
        prevalence={"Ruminococcus_bromii": 0.9}, n_participants=30,
    )
    result = cu.score_panel(cu.SD_PANEL, {"Ruminococcus_bromii": 1e-3}, frozen)
    subdo = next(r for r in result.rows if r.feature.source_id == "Subdoligranulum")
    assert subdo.state == "reference_ineligible"
    assert subdo.contribution is None
    assert result.usable_count == 1


def test_nondetection_with_a_known_limit_is_censored_in_the_right_direction(frozen_csu) -> None:
    """CU32: the interval endpoints are transformed then reordered by direction.

    Since v4.2 (V42-24) a censored nondetection is interval-only: it tightens
    the panel bounds but no point is invented for it.
    """
    limit = 1e-5
    result = cu.score_panel(
        cu.CSU_PANEL, _abundances(Ruminococcus_gnavus=0.0), frozen_csu, reporting_limit=limit
    )
    row = next(r for r in result.rows if r.feature.source_id == "Ruminococcus_gnavus")
    assert row.state == "below_reporting_limit"
    assert row.reporting_limit == limit
    assert row.contribution is None
    low, high = row.contribution_bounds
    assert low <= high
    # R. gnavus is an increase feature, so a nondetection bounds low alignment.
    assert high < 50.0
    assert result.usable_count == 14
    assert result.interval_coverage > result.coverage
    assert "not biological absence" in row.note


def test_nondetection_without_a_limit_keeps_only_its_bounds(frozen_csu) -> None:
    """CU33: no validated limit means the point is omitted, not guessed."""
    result = cu.score_panel(cu.CSU_PANEL, _abundances(Ruminococcus_gnavus=0.0), frozen_csu)
    row = next(r for r in result.rows if r.feature.source_id == "Ruminococcus_gnavus")
    assert row.state == "below_reporting_limit"
    assert row.contribution is None
    assert result.usable_count == 14
    assert "no validated reporting limit" in row.note


def test_one_usable_feature_is_labelled_limited_not_promoted(frozen_csu) -> None:
    """CU36: a one-of-fifteen result says so, with the bounds left wide."""
    abundances = {"Coprococcus_catus": 1e-3}
    result = cu.score_panel(cu.CSU_PANEL, abundances, frozen_csu)
    assert result.usable_count == 1
    assert result.limited_feature
    assert result.status == "computed_partial"
    lo, hi = result.full_panel_bounds
    assert hi - lo > 80.0
    assert any("limited-feature" in n for n in result.notes)
    assert result.to_json()["limited_feature_result"] is True


def test_every_panel_feature_gets_a_row_even_when_unavailable(frozen_csu) -> None:
    """CU67: all fifteen rows are visible, with state, and absence is not deficiency."""
    result = cu.score_panel(cu.CSU_PANEL, {"Coprococcus_catus": 1e-3}, frozen_csu)
    assert len(result.rows) == 15
    assert {r.feature.source_id for r in result.rows} == set(cu.CSU_PANEL.ids)
    assert all(r.state in cu.MEASUREMENT_STATES for r in result.rows)
    for row in result.rows:
        if row.note:
            assert "deficiency" not in row.note.lower()


def test_bounds_are_labelled_as_coverage_not_confidence(frozen_csu) -> None:
    """§6.3: the bounds must never be read as a confidence interval."""
    payload = cu.score_panel(cu.CSU_PANEL, _abundances(), frozen_csu).to_json()
    assert payload["bounds_type"] == "missingness_and_censoring_not_confidence_interval"


def test_scoring_is_independent_of_any_batch(frozen_csu) -> None:
    """CU30: a sample scored alone equals the same sample scored in a batch.

    Nothing about the query enters the fit, so there is no batch to depend on.
    """
    alone = cu.score_panel(cu.CSU_PANEL, _abundances(), frozen_csu).to_json()
    for _ in range(3):
        cu.score_panel(cu.CSU_PANEL, _abundances(Escherichia_coli=0.2), frozen_csu)
    again = cu.score_panel(cu.CSU_PANEL, _abundances(), frozen_csu).to_json()
    assert alone == again


def test_results_are_only_comparable_on_a_common_mask(frozen_csu) -> None:
    """CU38: a different feature mask means the two dates are not comparable."""
    full = cu.score_panel(cu.CSU_PANEL, _abundances(), frozen_csu)
    partial_abundances = _abundances()
    del partial_abundances["Coprococcus_catus"]
    partial = cu.score_panel(cu.CSU_PANEL, partial_abundances, frozen_csu)
    assert cu.comparable(full, full)
    assert not cu.comparable(full, partial)


def test_degenerate_reference_slots_are_ineligible_not_saturated() -> None:
    """A median-nondetection reference cannot produce a graded position.

    Rather than let such a feature clip to ±6 and act as a full-weight
    presence switch, the slot is reported as reference-ineligible: it lowers
    coverage and widens the bounds.
    """
    stats = {
        # a centre that is itself a nondetection: unusable
        "Subdoligranulum": {"mu": math.log10(cu.EPS), "scale": cu.SCALE_FLOOR, "n_controls": 30.0},
        # a real distribution among carriers: usable
        "Ruminococcus_bromii": {"mu": -3.0, "scale": 0.9, "n_controls": 1881.0},
    }
    frozen = cu.FrozenReference(
        bundle_id="fixture", scope="compatible_independent_controls", stats=stats,
        median_abundance={"Subdoligranulum": 0.0, "Ruminococcus_bromii": 1e-3},
        prevalence={"Subdoligranulum": 0.003, "Ruminococcus_bromii": 0.62}, n_participants=30,
    )
    result = cu.score_panel(cu.SD_PANEL, {"Subdoligranulum_variabile": 1e-4, "Ruminococcus_bromii": 1e-3}, frozen)
    row = next(r for r in result.rows if r.feature.source_id == "Subdoligranulum")
    assert row.state == "reference_ineligible"
    assert row.contribution is None
    assert "presence/absence switch" in row.note
    assert result.usable_count == 1 and result.limited_feature


def test_too_few_carriers_leaves_a_feature_unfitted() -> None:
    """A handful of carriers is not a distribution to rank against."""
    rows = [{"a": 0.01 if i < 10 else 0.0, "b": 0.02} for i in range(40)]
    stats = cu.fit_reference_detected(rows, ["a", "b"])
    assert "a" not in stats, "10 carriers must not produce a reference"
    assert stats["b"]["n_controls"] == 40.0


def test_control_nondetections_are_not_fitted_as_observed_zeros() -> None:
    """§6.1: the control side gets the same treatment as the query side.

    Including nondetections at the floor drags the centre down and collapses
    the spread. The production fit uses the detected values, so the centre is
    a real position among carriers.
    """
    # A quarter of the cohort carries the organism, which is typical of these
    # fifteen species and is where the literal formula breaks down.
    rows = [{"a": 0.01 if i % 4 == 0 else 0.0} for i in range(100)]
    literal = cu.fit_reference(rows, ["a"])["a"]
    detected = cu.fit_reference_detected(rows, ["a"])["a"]
    assert literal["mu"] == pytest.approx(math.log10(cu.EPS))
    assert literal["scale"] == cu.SCALE_FLOOR  # no spread left to standardise by
    assert cu.z_value(0.01, literal) == cu.Z_CLIP  # every carrier clips
    assert detected["mu"] == pytest.approx(math.log10(0.01 + cu.EPS))
    assert detected["n_controls"] == 25.0
    assert abs(cu.z_value(0.01, detected)) < 1e-9  # a carrier at the centre


def test_the_real_catalogue_produces_a_typed_result() -> None:
    """End to end against the shipped reference cohort, if it is present."""
    import json

    path = REPO / "refs" / "taxonomic_cohort.json"
    if not path.is_file():
        pytest.skip("no reference cohort built in this checkout")
    cohort = json.loads(path.read_text())
    catalogue, matrix = cohort["taxa"], cohort["abundance"]
    sample = {t: matrix[i][0] for i, t in enumerate(catalogue)}  # percent scale

    class _Cohort:
        taxa = catalogue
        abundance = matrix
        snapshot = "test"
        n_samples = len(matrix[0])

    results = cu.urticaria_indices(sample, _Cohort())
    assert set(results) == {"CSU_STOOL_WGS_V1", "SD_STOOL_TRANSLATED_V1"}
    csu_result = results["CSU_STOOL_WGS_V1"]
    assert csu_result.panel_count == 15
    assert csu_result.sensitivity is not None
    assert csu_result.sensitivity.panel_count == 8
    if csu_result.computed:
        lo, hi = csu_result.full_panel_bounds
        assert 0.0 <= lo <= csu_result.score <= hi <= 100.0
    payload = csu_result.to_json()
    assert payload["score_semantics"] == "unvalidated_directional_pattern_index"
    assert payload["disease_specificity"] == "not_established"
    assert payload["treatment_recommendation"] is None
    assert payload["donor_eligibility"] is None
    assert results["SD_STOOL_TRANSLATED_V1"].panel.assay_transport == "amplicon_taxon_to_shotgun"


# --------------------------------------------------------------------------- #
# intervention evidence (CU63, CU64, CU65, CU71, CU78, CU79, CU80, CU81)
# --------------------------------------------------------------------------- #


def _cu_mappings(registry) -> list:
    return [
        m for m in registry.mappings
        if m.trigger_selector.target in URTICARIA_PROFILES
    ]


def test_hives_cards_exist_for_both_profiles(registry) -> None:
    targets = {m.trigger_selector.target for m in _cu_mappings(registry)}
    assert targets == set(URTICARIA_PROFILES)


def test_retrieval_is_balanced_not_positive_only(registry) -> None:
    """CU64: the null and mostly-non-response studies are retrieved too."""
    csu_cards = [m for m in _cu_mappings(registry) if m.trigger_selector.target == "csu"]
    refs = {r for m in csu_cards for r in (*m.evidence_assertion_refs, *m.against_evidence_assertion_refs)}
    # the favourable non-randomized study and the null trial both present
    assert "IE_CU_PROBIOTIC_CAN_2024_V1" in refs
    assert "IE_CU_PROBIOTIC_GODSE_2025_X_V1" in refs
    assert "IE_CU_PROBIOTIC_NETTIS_2016_V1" in refs
    lanes = {registry.assertions[r].lane for r in refs}
    assert "X" in lanes, "no against-evidence retrieved for hives"
    assert lanes & {"A", "C", "D"}, "no supporting evidence retrieved for hives"


def test_nettis_is_recorded_as_mostly_non_response(registry) -> None:
    a = registry.assertions["IE_CU_PROBIOTIC_NETTIS_2016_V1"]
    assert a.lane == "D"
    summary = a.summary.lower()
    assert "27" in summary and "non-response" in summary
    assert "no control arm" in summary or "no comparator" in summary


def test_bi_2021_formulation_and_population_are_exact(registry) -> None:
    """CU78: B. animalis LK011, not B. breve, and paediatric evidence only."""
    protocol = registry.protocols["PROTOCOL_CU_YIMINGJIA_4W_PAEDIATRIC_V1"]
    intervention = registry.interventions[protocol.intervention_id]
    taxa = {c["taxon"] for c in intervention.identity["components"]}
    assert "Bifidobacterium animalis" in taxa
    assert "Bifidobacterium breve" not in taxa
    strains = {c.get("strain_designation") for c in intervention.identity["components"]}
    assert "LK011" in strains
    assert protocol.population["age_range_years"] == [6, 12]
    assert any("6–12" in w or "children" in w.lower() for w in intervention.label_warnings)
    assert "children" in registry.assertions["IE_CU_PROBIOTIC_BI_2021_V1"].summary.lower()


def test_can_2024_strain_and_design_are_exact(registry) -> None:
    """CU79: ATCC 55730 and non-randomized — not DSM 17938, not an RCT."""
    a = registry.assertions["IE_CU_PROBIOTIC_CAN_2024_V1"]
    assert a.source["source_type"] == "nonrandomized_study"
    assert "not randomized" in a.summary.lower() or "non-randomized" in a.summary.lower()
    intervention = registry.interventions["PROBIOTIC_LREUTERI_ATCC55730_V1"]
    assert "ATCC 55730" in intervention.display_name
    assert any("DSM 17938" in w for w in intervention.label_warnings)


def test_source_conflicts_are_displayed_not_smoothed(registry) -> None:
    """CU80: Godse and Dabaghzadeh keep their internal inconsistencies."""
    godse = registry.assertions["IE_CU_PROBIOTIC_GODSE_2025_X_V1"]
    assert godse.lane == "X" and godse.against_kind == "null_primary"
    assert "baseline" in godse.summary.lower()
    daba = registry.assertions["IE_CU_PROBIOTIC_DABAGHZADEH_2023_V1"]
    assert "main effect" in daba.summary.lower()
    assert "conflict" in daba.summary.lower() or "inconsistency" in daba.summary.lower()


def test_missing_strain_ids_are_declared_not_invented(registry) -> None:
    """CU63: an unreported strain stays unreported."""
    for iid in ("PROBIOTIC_LACTOCARE_SEVEN_ORGANISM_V1", "PROBIOTIC_FEMILACT_SEVEN_ORGANISM_V1"):
        intervention = registry.interventions[iid]
        strains = {c.get("strain_designation") for c in intervention.identity["components"] if c["entity_type"] == "probiotic_strain"}
        assert strains == {"not_reported"}
        assert intervention.exact_protocol_reproducibility == "incomplete"
        assert intervention.label_warnings


def test_the_fmt_case_report_is_present_and_unverified(registry) -> None:
    """CU81: a real case report, with its protocol fields left empty."""
    a = registry.assertions["IE_CU_FMT_WU_2023_V1"]
    assert a.source["pmid"] == 37077050
    assert a.lane == "X"
    protocol = registry.protocols["PROTOCOL_FMT_CU_CASE_WU_V1"]
    assert protocol.population["age_range_years"] is None
    assert protocol.administration["route"] == "not_reported"
    assert "not verified" in protocol.administration["dose_source_text"]


def test_no_card_says_fmt_is_the_only_solution(registry) -> None:
    """§10.4: neither 'no human report exists' nor 'FMT is the only option'."""
    fmt = next(m for m in registry.mappings if m.action_mapping_id == "MAP_PROFILE_CSU_FMT_V1")
    assert "FMT is the only solution" in fmt.never_say
    assert fmt.fmt_status == "investigational_non_cdi"
    assert "not the only option" in fmt.context_note
    assert "no controlled human trial" in fmt.context_note.lower()


def test_antimicrobial_stacks_are_explicitly_unsupported(registry) -> None:
    """The user's own antimicrobial experience is context, not a template."""
    a = registry.assertions["IE_CU_ANTIMICROBIAL_STACK_X_V1"]
    assert a.lane == "X"
    summary = a.summary.lower()
    for agent in ("berberine", "allicin", "oregano", "black seed"):
        assert agent in summary
    assert "no human study" in summary


def test_cards_never_claim_a_cure_or_a_microbiome_selected_treatment(registry) -> None:
    """CU71: affirmative cure, eradication and selection claims are blocked."""
    banned = (
        "cures hives", "cure your hives", "treats urticaria", "will clear your hives",
        "eradicate the bacteria causing", "microbiome-selected", "guaranteed",
    )
    for mapping in _cu_mappings(registry):
        text = f"{mapping.title} {mapping.context_note}".lower()
        for phrase in banned:
            assert phrase not in text, f"{mapping.action_mapping_id} claims: {phrase}"
        assert mapping.never_say, f"{mapping.action_mapping_id} declares nothing off-limits"
        assert mapping.render in (
            "evidence_only", "clinician_review", "guideline_route_after_confirmation",
            "inform_only", "not_actionable", "no_supported_targeted_intervention",
        )
    for aid, assertion in registry.assertions.items():
        if not aid.startswith("IE_CU_"):
            continue
        assert assertion.target_binding["microbiome_profile_alone_sufficient"] is False


def test_every_hives_card_requires_a_clinical_confirmation(registry) -> None:
    """CU57/CU58: the profile alone unlocks education, never a prescription."""
    for mapping in _cu_mappings(registry):
        assert mapping.confirmation_rule_refs, f"{mapping.action_mapping_id} needs a confirmation rule"
    donor_words = ("donor pass", "eligible donor", "cleared for donation", "donor rank")
    text = " ".join(
        f"{m.title} {m.context_note}" for m in _cu_mappings(registry)
    ).lower()
    for phrase in donor_words:
        assert phrase not in text


def test_guideline_route_names_antihistamines_and_a_prescriber(registry) -> None:
    a = registry.assertions["IE_CU_GUIDELINE_2026_V1"]
    assert a.lane == "A"
    assert "antihistamine" in a.summary.lower()
    assert "clinician" in a.summary.lower() or "prescriber" in a.summary.lower()


def test_vitamin_d_and_hpylori_need_a_lab_result(registry) -> None:
    """CU54/§10.3: neither route can be triggered from stool DNA."""
    for aid in ("IE_CU_VITD_MONY_2020_V1", "IE_CU_HPYLORI_PAWLOWICZ_2018_V1"):
        a = registry.assertions[aid]
        assert a.target_binding["trigger_type"] == "external_lab_abnormality"
    valsecchi = registry.assertions["IE_CU_HPYLORI_VALSECCHI_1998_X_V1"]
    assert valsecchi.lane == "X"


def test_animal_only_evidence_is_labelled_preclinical(registry) -> None:
    """CU61: LG-1 and short-chain fatty acid treatment were given to mice."""
    a = registry.assertions["IE_CU_SCFA_LG1_PRECLINICAL_E_V1"]
    assert a.lane == "E"
    summary = a.summary.lower()
    assert "mice" in summary or "mouse" in summary
    assert "no human dose" in summary


# --------------------------------------------------------------------------- #
# the profile YAML as a document
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", URTICARIA_PROFILES)
def test_profile_yaml_declares_its_provenance(name) -> None:
    raw = yaml.safe_load((PROFILES / f"{name}.yaml").read_text())
    assert raw["status"] == "research_only"
    assert raw["citation"] and "doi" in raw["citation"].lower()
    assert raw["references"]["primary"]["source"] is None
    assert raw["caveats"]
    assert raw["spec_version"] == 3


def test_no_ellipsis_in_hives_prose(profile_set) -> None:
    """The report never truncates a sentence with an ellipsis."""
    for name in URTICARIA_PROFILES:
        profile = next(p for p in profile_set.profiles if p.name == name)
        blocks = [profile.summary, profile.effect_size_note, profile.citation, *profile.caveats]
        blocks += [f.note or "" for f in profile.features()]
        for block in blocks:
            assert "\u2026" not in block
            assert "..." not in block
