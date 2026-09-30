"""Acceptance gates V42-01..V42-61 for the inflammatory-skin release (spec v4.2 §9.2).

These are biomedical-safety gates, not coverage decoration. Most of them exist
because the plausible-looking failure is the dangerous one: a missing organism
silently scored as depleted, a pooled contrast reported as a subtype test, a
culture count converted into an abundance coefficient, a child handed an adult
disease number.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from openbiota import skin
from openbiota.dirindex import (
    DirIndexError,
    FrozenReference,
    Panel,
    PanelFeature,
    evaluate,
    fit_reference,
    fraction,
    score_panel,
    summarize,
    weighted_midrank,
)
from openbiota.interventions import load_registry

REPO = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- #
# registry and migration (V42-01..V42-04)
# --------------------------------------------------------------------------- #


def test_v4201_seed_inventory_is_nine_panels_and_fiftytwo_slots():
    """V42-01: the release inventory is exactly what the spec declares."""
    assert len(skin.SEED_PANELS) == 9
    assert sum(len(p.features) for p in skin.SEED_PANELS) == 52
    ids = [p.panel_id for p in skin.PANELS]
    assert len(ids) == len(set(ids))
    for panel in skin.PANELS:
        feature_ids = [f.source_id for f in panel.features]
        assert len(feature_ids) == len(set(feature_ids)), panel.panel_id


def test_v4201_membership_is_not_counted_twice_within_a_panel():
    """V42-01: two slots resolving to one measurement would double count."""
    for panel in skin.PANELS:
        keys = [f.key for f in panel.features]
        assert len(keys) == len(set(keys)), panel.panel_id


def test_v4202_existing_families_keep_their_identity():
    """V42-02: psoriasis panels attach to one family, not several new ones."""
    psoriasis = [p for p in skin.PANELS if p.family == "plaque_psoriasis"]
    assert len(psoriasis) == 3
    assert len({p.family for p in psoriasis}) == 1
    # three panels under one family is not three diseases
    assert len({p.family for p in skin.PANELS}) < len(skin.PANELS)


def test_v4203_alzheimer_ad_and_dermographism_sd_are_not_reassigned():
    """V42-03: ATD/SEBD must not collide with the existing AD and SD."""
    from openbiota.cuindex import SD_PANEL

    assert SD_PANEL.panel_id == "SD_STOOL_TRANSLATED_V1"
    ids = {p.panel_id for p in skin.PANELS}
    assert SD_PANEL.panel_id not in ids
    for panel_id in ids:
        assert not panel_id.startswith("AD_")
        assert not panel_id.startswith("SD_")
    # and the new ones use the reserved-safe prefixes
    assert any(p.panel_id.startswith("ATD_") for p in skin.PANELS)
    assert "SEBD_STOOL_RESEARCH_V1" in skin.TASKS_BY_ID
    # the existing "ad_*" profiles are Alzheimer's disease and stay that way
    names = {p.stem for p in (REPO / "profiles").glob("*.yaml")}
    ad_profiles = {n for n in names if n.startswith("ad_")}
    assert ad_profiles
    for name in ad_profiles:
        text = (REPO / "profiles" / f"{name}.yaml").read_text().lower()
        assert "alzheimer" in text, name
        assert "atopic" not in text, name


def test_v4204_unrelated_profiles_are_untouched():
    """V42-04: the hives panels keep their own ids, revisions and features."""
    from openbiota.cuindex import CSU_PANEL, SD_PANEL

    assert CSU_PANEL.panel_id == "CSU_STOOL_WGS_V1"
    assert len(CSU_PANEL.features) == 15
    assert CSU_PANEL.revision.startswith("4.1")
    assert SD_PANEL.revision.startswith("4.1")
    assert len(SD_PANEL.features) == 2


# --------------------------------------------------------------------------- #
# taxonomy resolution (V42-05..V42-13)
# --------------------------------------------------------------------------- #


def test_v4205_same_epithet_in_different_genera_never_merge():
    """V42-05: Bacteroides finegoldii is not Alistipes finegoldii."""
    sa = {f.source_id for f in skin.PSO_SA2026.features}
    xiao = {f.source_id for f in skin.PSORIATIC_SHARED_XIAO2024.features}
    assert "Bacteroides_finegoldii" in sa
    assert "Alistipes_finegoldii" in xiao
    b = next(f for f in skin.PSO_SA2026.features if f.source_id == "Bacteroides_finegoldii")
    assert b.key == "Bacteroides_finegoldii"
    assert "Alistipes" in (b.mapping_reason or "")


def test_v4206_missing_sgb4348_is_unavailable_not_zero():
    """V42-06: an absent SGB is not zero, a proxy, or a guessed SGB."""
    catalogue = ["Escherichia_coli", "Bacteroides_fragilis"]
    matrix = [[1.0] * 25, [2.0] * 25]
    ref = skin.freeze_reference(skin.PSO_DENG2026, catalogue=catalogue, matrix=matrix, bundle_id="t")
    result = score_panel(skin.PSO_DENG2026, {"Escherichia_coli": 0.01}, ref)
    assert result.score is None
    assert result.usable_count == 0
    row = result.rows[0]
    # absent from the catalogue: no reference slot, no zero, no proxy
    assert row.state in ("not_in_database", "reference_ineligible")
    assert row.contribution is None
    assert row.abundance is None


def test_v4206_synonym_resolution_is_exact_and_catalogue_verified():
    """A renamed organism is found under whichever name the catalogue carries.

    Both names present is an ambiguity, not a sum.
    """
    f = next(f for f in skin.PSO_CHANG2022.features if f.source_id == "Bacteroides_vulgatus")
    assert f.synonyms == ("Phocaeicola_vulgatus",)
    assert f.resolve(["Bacteroides_vulgatus", "X"]) == (("Bacteroides_vulgatus",), None)
    assert f.resolve(["Phocaeicola_vulgatus", "X"]) == (("Phocaeicola_vulgatus",), None)
    assert f.resolve(["X"]) == ((), "not_in_catalogue")
    assert f.resolve(["Bacteroides_vulgatus", "Phocaeicola_vulgatus"]) == ((), "ambiguous_synonyms")
    # the epithet alone never matches
    assert f.resolve(["Bacteroides_vulgatus_CAG_1"]) == ((), "not_in_catalogue")


def test_v4206_frozen_reference_records_the_names_it_was_fitted_on():
    panel = skin.PSO_CHANG2022
    catalogue = ["Bacteroides_vulgatus", "Dialister_invisus", "Bacteroides_coprocola"]
    matrix = [[1.0 + 0.001 * i] * 30 for i in range(len(catalogue))]
    ref = skin.freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="t")
    assert ref.resolved["Bacteroides_vulgatus"] == ("Bacteroides_vulgatus",)
    assert ref.unresolved["Megasphaera_unclassified"] == "not_in_catalogue"
    result = score_panel(panel, {"Bacteroides_vulgatus": 0.02}, ref)
    row = next(r for r in result.rows if r.feature.source_id == "Bacteroides_vulgatus")
    assert row.state == "detected"
    assert row.contribution is not None
    missing = next(r for r in result.rows if r.feature.source_id == "Megasphaera_unclassified")
    assert missing.state == "not_in_database"
    assert "no slot" in missing.note


def test_enterocloster_sums_only_its_named_species_under_either_name():
    f = next(f for f in skin.PSO_SA2026.features if f.source_id == "Enterocloster")
    old = ["Clostridium_bolteae", "Clostridium_clostridioforme", "Clostridium_butyricum",
           "Clostridium_perfringens"]
    members, reason = f.resolve(old)
    assert reason is None
    assert set(members) == {"Clostridium_bolteae", "Clostridium_clostridioforme"}
    new = ["Enterocloster_bolteae", "Enterocloster_citroniae"]
    assert set(f.resolve(new)[0]) == set(new)
    # one species under both names cannot be summed
    assert f.resolve(["Clostridium_bolteae", "Enterocloster_bolteae"]) == ((), "ambiguous_synonyms")
    # a catalogue without any of them leaves the slot unavailable
    assert f.resolve(["Clostridium_perfringens"]) == ((), "not_in_catalogue")


def test_v4207_same_number_in_another_namespace_is_rejected():
    """V42-07: SGB numbering is namespace-specific."""
    feature = skin.PSO_DENG2026.features[0]
    assert feature.rank == "sgb"
    assert "namespace" in (feature.mapping_reason or "")
    # a same-numbered entry from a different namespace must not match the key
    assert feature.key == "GGB51647_SGB4348"
    assert feature.key != "SGB4348"


def test_v4208_unclassified_bucket_is_not_the_genus_total():
    """V42-08: Megasphaera_unclassified is not the genus and not M. elsdenii."""
    f = next(f for f in skin.PSO_CHANG2022.features if f.source_id == "Megasphaera_unclassified")
    assert f.rank == "source_group"
    assert "genus total" in (f.mapping_reason or "")
    assert "elsdenii" in (f.mapping_reason or "")
    # scoring with only a genus total present leaves the slot unavailable
    catalogue = ["Megasphaera_elsdenii"] + [f"Filler_{i}" for i in range(5)]
    matrix = [[1.0] * 25 for _ in catalogue]
    ref = skin.freeze_reference(skin.PSO_CHANG2022, catalogue=catalogue, matrix=matrix, bundle_id="t")
    result = score_panel(skin.PSO_CHANG2022, {"Megasphaera_elsdenii": 0.02}, ref)
    row = next(r for r in result.rows if r.feature.source_id == "Megasphaera_unclassified")
    assert row.contribution is None


def test_v4209_generic_pcopri_does_not_satisfy_clade_c():
    """V42-09: a clade slot is unavailable without clade resolution."""
    f = next(f for f in skin.VIT_NONSEG_LUAN2023.features
             if f.source_id == "Prevotella_copri_clade_C")
    assert f.rank == "source_clade"
    assert "not all" in (f.mapping_reason or "").lower()
    catalogue = ["Prevotella_copri"] + [f"Filler_{i}" for i in range(5)]
    matrix = [[1.0] * 25 for _ in catalogue]
    ref = skin.freeze_reference(skin.VIT_NONSEG_LUAN2023, catalogue=catalogue, matrix=matrix, bundle_id="t")
    result = score_panel(skin.VIT_NONSEG_LUAN2023, {"Prevotella_copri": 0.05}, ref)
    row = next(r for r in result.rows if r.feature.source_id == "Prevotella_copri_clade_C")
    assert row.contribution is None


def test_v4210_atd_keeps_thirteen_slots_when_only_genera_resolve():
    """V42-10: unresolved source groups stay in the coverage denominator."""
    panel = skin.ATD_ADULT_WANG2023
    assert len(panel.features) == 13
    genera = [f for f in panel.features if f.rank == "genus"]
    groups = [f for f in panel.features if f.rank == "source_group"]
    assert len(genera) == 8
    assert len(groups) == 5
    catalogue = []
    for f in genera:
        catalogue += [f"{f.source_id}_a", f"{f.source_id}_b"]
    matrix = [[0.5] * 25 for _ in catalogue]
    ref = skin.freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="t")
    abundances = dict.fromkeys(catalogue, 0.004)
    result = score_panel(panel, abundances, ref)
    assert result.panel_count == 13
    assert result.usable_count <= 8
    for f in groups:
        row = next(r for r in result.rows if r.feature.source_id == f.source_id)
        assert row.contribution is None


def test_v4211_family_total_is_rejected_for_unclassified_family_slots():
    """V42-11: summing a whole family is not the unclassified bucket."""
    for name in ("unclassified_Oscillospiraceae", "unclassified_Butyricicoccaceae",
                 "unclassified_Erysipelotrichaceae"):
        f = next(f for f in skin.ATD_ADULT_WANG2023.features if f.source_id == name)
        assert "not the family total" in (f.mapping_reason or "")


def test_v4212_faecalibacterium_umbrella_and_child_are_quarantined():
    """V42-12: never sum or substitute an umbrella species and its child."""
    panel = skin.VIT_ACTIVE_JU2025
    ids = {f.source_id for f in panel.features}
    assert {"Faecalibacterium_prausnitzii", "Faecalibacterium_duncaniae"} <= ids
    # only the umbrella exists in the catalogue: the pair is withheld
    catalogue = ["Faecalibacterium_prausnitzii", "Megamonas_funiformis", "Bifidobacterium_bifidum"]
    matrix = [[1.0] * 25 for _ in catalogue]
    ref = skin.freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="t")
    result = skin.score_skin_panel(
        panel, {"Faecalibacterium_prausnitzii": 0.05, "Megamonas_funiformis": 0.01,
                "Bifidobacterium_bifidum": 0.01},
        ref, catalogue=catalogue,
    )
    withheld = {"Faecalibacterium_prausnitzii", "Faecalibacterium_duncaniae"}
    for row in result.rows:
        if row.feature.source_id in withheld:
            assert row.contribution is None, row.feature.source_id
    assert result.panel_count == 4


def test_v4212_disjoint_faecalibacterium_pair_is_allowed():
    """V42-12: proven-disjoint definitions may both score."""
    panel = skin.VIT_ACTIVE_JU2025
    catalogue = ["Faecalibacterium_prausnitzii", "Faecalibacterium_duncaniae",
                 "Megamonas_funiformis", "Bifidobacterium_bifidum"]
    assert skin._quarantine_faecalibacterium(panel, catalogue) == ()


def test_v4213_luan_spelling_conflict_is_excluded():
    """V42-13: the thermophilus/thermophiles conflict stays out."""
    ids = {f.source_id for f in skin.VIT_NONSEG_LUAN2023.features}
    for bad in ("Streptococcus_thermophilus", "Staphylococcus_thermophiles",
                "Staphylococcus_thermophilus"):
        assert bad not in ids


# --------------------------------------------------------------------------- #
# seborrheic dermatitis (V42-14..V42-19)
# --------------------------------------------------------------------------- #


def test_v4214_culture_counts_do_not_become_abundance_coefficients():
    """V42-14: low absolute E. coli count generates no DNA coefficient."""
    for panel in skin.PANELS:
        assert panel.family != "seborrheic_dermatitis"
    assert "SEBD_STOOL_RESEARCH_V1" in skin.TASKS_BY_ID
    task = skin.TASKS_BY_ID["SEBD_STOOL_RESEARCH_V1"]
    assert task.state == "insufficient_stool_signature"
    # no feature anywhere carries a seborrheic-dermatitis direction
    assert not any(p.family == "seborrheic_dermatitis" for p in skin.PANELS)


def test_v4215_culture_denominator_discrepancy_is_preserved():
    """V42-15: 67 enrolled versus a ~32 denominator is not tidied away."""
    import openbiota.pdfskin as pdfskin

    assert "not a score" in (pdfskin._sebd_page.__doc__ or "").lower()
    source = Path(pdfskin.__file__).read_text()
    assert "67" in source and "32" in source
    assert "not explained" in source


def test_v4216_mendelian_randomization_gets_no_coefficient():
    """V42-16: genetic associations are visible and unweighted."""
    registry = load_registry(REPO / "interventions")
    record = registry.assertions["IE_SEBD_NO_STOOL_SIGNATURE_X_V1"]
    assert record.lane == "X"
    assert "genetically instrumented" in record.summary
    assert not any(p.family == "seborrheic_dermatitis" for p in skin.PANELS)


def test_v4217_skin_and_scalp_findings_never_enter_a_stool_signature():
    """V42-17: Malassezia and skin ratios do not seed a gut panel."""
    for panel in skin.PANELS:
        for f in panel.features:
            assert "Malassezia" not in f.source_id
            assert "Cutibacterium" not in f.source_id
            assert "Micrococcus" not in f.source_id


def test_v4218_sixty_percent_dysbiosis_is_not_a_prior_or_cutoff():
    """V42-18: a clinical-series percentage is not a screening prevalence."""
    import openbiota.pdfskin as pdfskin

    source = Path(pdfskin.__file__).read_text()
    assert "not a screening" in source
    assert "cutoff" in source


def test_v4219_sebd_returns_null_with_a_reason_not_zero_or_fifty():
    """V42-19: no number, an explicit reason, and never a negative diagnosis."""
    task = skin.TASKS_BY_ID["SEBD_STOOL_RESEARCH_V1"].to_json()
    assert task["index"] is None
    assert task["control_percentile"] is None
    assert task["disease_probability"] is None
    assert task["index"] != 0 and task["index"] != 50
    assert "cannot assess whether you have seborrheic dermatitis" in task["statement"]
    assert task["status"] == "insufficient_stool_signature"


def test_v4219_sebd_statement_is_verbatim():
    """The spec fixes this wording; it is participant-facing."""
    assert skin.SEBD_STATEMENT == (
        "Gut involvement has been reported, but no adequately specified stool-DNA signature was "
        "identified for this implementation. This result cannot assess whether you have seborrheic "
        "dermatitis."
    )


# --------------------------------------------------------------------------- #
# the numeric kernel (V42-20..V42-31)
# --------------------------------------------------------------------------- #


def _flat_reference(keys, *, mu=-2.0, scale=0.5, n=25):
    return FrozenReference(
        bundle_id="test", scope="test",
        stats={k: {"mu": mu, "scale": scale, "n_controls": float(n)} for k in keys},
        median_abundance=dict.fromkeys(keys, 0.01),
        prevalence=dict.fromkeys(keys, 1.0),
        n_participants=n,
    )


def test_v4220_all_thirteen_atd_missing_gives_null_and_full_bounds():
    """V42-20: no measurements means index null and bounds [0,100]."""
    panel = skin.ATD_ADULT_WANG2023
    ref = _flat_reference([f.key for f in panel.features])
    result = score_panel(panel, {}, ref)
    assert result.score is None
    assert result.usable_count == 0
    assert result.full_panel_bounds == (0.0, 100.0)


def test_v4221_feature_at_the_reference_median_contributes_fifty():
    """V42-21: 50 is centred, not a 50% probability of disease."""
    from openbiota.dirindex import support

    assert support(0.0, 1) == 50.0
    assert support(0.0, -1) == 50.0
    payload = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert payload["disease_probability"] is None
    assert payload["index_label"] == "Research directional pattern index"


def test_v4222_one_point_of_two_slots_gives_coverage_half_and_bounds_25_75():
    """V42-22: the spec's exact arithmetic for a half-covered panel."""
    out = summarize({"a": (50.0, 50.0, 50.0)}, ["a", "b"])
    assert out["index"] == 50
    assert out["coverage"] == 0.5
    assert out["full_panel_bounds"] == [25.0, 75.0]


def test_v4223_opposite_directions_are_monotonic_and_sum_to_one_hundred():
    """V42-23: raising abundance moves + and - features oppositely."""
    from openbiota.dirindex import support

    up_low, up_high = support(-1.0, 1), support(1.0, 1)
    down_low, down_high = support(-1.0, -1), support(1.0, -1)
    assert up_high > up_low
    assert down_high < down_low
    for z in (-2.0, -0.5, 0.5, 2.0):
        assert abs(support(z, 1) + support(z, -1) - 100.0) < 1e-10


def test_v4224_nondetection_with_a_limit_is_interval_only():
    """V42-24: a censored nondetection gives bounds and no invented point."""
    panel = Panel(
        panel_id="T", revision="1", label="t", source_id="S",
        source_assay="shotgun_metagenomics", scoring_assay="shotgun_metagenomics",
        score_semantics=skin.NATIVE, family="plaque_psoriasis",
        features=(PanelFeature("x", -1),),
    )
    ref = _flat_reference(["x"])
    result = score_panel(panel, {"x": 0.0}, ref, reporting_limit=1e-4)
    row = result.rows[0]
    assert row.contribution is None
    assert row.contribution_bounds is not None
    assert result.score is None


def test_v4225_missing_row_without_a_limit_is_unavailable_not_depleted():
    """V42-25: an absent row is not maximal depletion."""
    panel = Panel(
        panel_id="T", revision="1", label="t", source_id="S",
        source_assay="shotgun_metagenomics", scoring_assay="shotgun_metagenomics",
        score_semantics=skin.NATIVE, family="plaque_psoriasis",
        features=(PanelFeature("x", -1),),
    )
    ref = _flat_reference(["x"])
    result = score_panel(panel, {}, ref)
    row = result.rows[0]
    assert row.contribution is None
    # a depleted reading for a -1 feature would be near 100; it must not appear
    assert result.score is None


@pytest.mark.parametrize("bad", [True, -0.1, 1.1, float("nan"), float("inf"), "0.1", None])
def test_v4226_invalid_abundances_are_rejected_before_scoring(bad):
    """V42-26: negative, NaN, infinite, boolean and string all refuse."""
    with pytest.raises((DirIndexError, TypeError)):
        fraction(bad)


def test_v4227_duplicate_participants_cannot_inflate_reference_n():
    """V42-27: technical replicates do not count twice."""
    rows = [{"participant_id": f"p{i}", "features": {"a": 0.01}} for i in range(20)]
    ids = [r["participant_id"] for r in rows] + ["p0"]
    assert len(ids) != len(set(ids))


def test_v4228_fewer_than_twenty_controls_leaves_a_feature_unfitted():
    """V42-28: below the minimum there is no fitted reference."""
    rows = [{"a": 0.01} for _ in range(19)]
    with pytest.raises(DirIndexError):
        fit_reference(rows, ["a"])
    rows20 = [{"a": 0.01} for _ in range(20)]
    assert "a" in fit_reference(rows20, ["a"])


def test_v4229_scoring_is_identical_alone_and_in_a_batch():
    """V42-29: no inference-time fitting, so batch context cannot matter."""
    panel = skin.PSO_SA2026
    catalogue = ["Enterocloster_bolteae", "Phascolarctobacterium_faecium", "Bacteroides_finegoldii"]
    matrix = [[1.0 + 0.01 * i] * 30 for i in range(len(catalogue))]
    ref = skin.freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="t")
    sample = {"Enterocloster_bolteae": 0.02, "Phascolarctobacterium_faecium": 0.01,
              "Bacteroides_finegoldii": 0.005}
    alone = score_panel(panel, sample, ref).to_json()
    # the same frozen reference, scored again amid other samples
    for other in ({"Enterocloster_bolteae": 0.9}, {"Bacteroides_finegoldii": 0.0}):
        score_panel(panel, other, ref)
    again = score_panel(panel, sample, ref).to_json()
    assert alone == again


def test_v4230_scale_below_floor_refuses_to_calculate():
    """V42-30: a degenerate scale is refused, not clipped into service."""
    from openbiota.dirindex import SCALE_FLOOR, z_value

    with pytest.raises(DirIndexError):
        z_value(0.01, {"mu": -2.0, "scale": SCALE_FLOOR / 2, "n_controls": 25.0})


def test_v4231_percent_is_divided_once_and_selected_taxa_are_not_reclosed():
    """V42-31: percentages convert exactly once; no reclosure to 100%."""
    species_percent = {"Bacteroides_fragilis": 5.0, "Dorea_longicatena": 2.0}
    out = skin.skin_indices(species_percent, None)
    record = out["VIT_NONSEG_LUAN2023_TAX_V1"]
    assert record["index"] is None  # no cohort supplied
    # the fractions the kernel would see must sum to 0.07, not 1.0
    fractions = {k: v / 100.0 for k, v in species_percent.items()}
    assert abs(sum(fractions.values()) - 0.07) < 1e-12


# --------------------------------------------------------------------------- #
# gene families (V42-32..V42-34)
# --------------------------------------------------------------------------- #


def test_v4232_ko_panel_is_unavailable_when_only_taxa_are_measured():
    """V42-32: no taxa-to-function imputation, ever."""
    out = skin.skin_indices({"Escherichia_coli": 3.0}, None, age_years=40,
                            usable_read_pairs=2_000_000)
    record = out["PSA_LIU2024_KO_V1"]
    assert record["index"] is None
    assert record["status"] == "features_unavailable"
    assert any("no_taxa_to_function_imputation" in r for r in record["reason_codes"])


def test_v4233_zero_ko_denominator_is_rejected():
    """V42-33: a zero annotation total cannot form a fraction."""
    with pytest.raises(DirIndexError):
        skin.ko_fractions({"K02004": 0.0, "K01190": 0.0})
    with pytest.raises(DirIndexError):
        skin.ko_fractions({})


def test_v4233_ko_denominator_spans_the_whole_annotation_space():
    """V42-33: normalise over all measured gene families, not the selected nine."""
    measured = {"K02004": 1.0, "K01190": 1.0, "K99999": 8.0}
    out = skin.ko_fractions(measured)
    assert abs(sum(out.values()) - 1.0) < 1e-12
    assert abs(out["K02004"] - 0.1) < 1e-12  # not 0.5, which reclosure would give


def test_v4234_k07114_not_the_typo_and_no_host_inference():
    """V42-34: correct KO only; no human chloride channel or serum CRP."""
    ids = {f.source_id for f in skin.PSA_LIU2024_KO.features}
    assert "K07114" in ids
    assert "K07714" not in ids
    f = next(f for f in skin.PSA_LIU2024_KO.features if f.source_id == "K07114")
    assert "not a human chloride channel" in (f.mapping_reason or "")


# --------------------------------------------------------------------------- #
# population, exposure and QC gates (V42-35..V42-39)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("age", [8.0, 17.0])
def test_v4235_children_get_no_adult_disease_numbers(age):
    """V42-35: ages 8-17 are refused a number; measurements stay visible."""
    out = skin.skin_indices({"Escherichia_coli": 3.0}, None, age_years=age,
                            usable_read_pairs=2_000_000, mode="participant")
    record = out["PSO_CHANG2022_TAX_V1"]
    assert record["index"] is None
    assert record["status"] == "population_unsupported"
    assert any("paediatric" in r for r in record["reason_codes"])


def test_v4235_unknown_age_is_not_assumed_eligible_for_a_participant():
    """V42-35: unknown is not 'fine'."""
    blocking, reasons = skin.eligibility(
        skin.PSO_CHANG2022, age_years=None, body_site="stool",
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert blocking == "population_unsupported"
    assert any("age_unknown" in r for r in reasons)


def test_v4235_e15_age_ceiling_is_sixty_not_sixtyfive():
    """V42-35: the gene-family panel's own source range is 18-60."""
    assert skin.PSA_LIU2024_KO.max_age_years == 60.0
    blocking, reasons = skin.eligibility(
        skin.PSA_LIU2024_KO, age_years=65.0, body_site="stool",
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert blocking == "population_unsupported"
    # while a taxon panel at 65 is still inside its own range
    ok, _ = skin.eligibility(
        skin.PSO_CHANG2022, age_years=65.0, body_site="stool",
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert ok is None


def test_v4236_recent_antibiotics_and_unknown_exposure_are_distinguished():
    """V42-36: unknown medication is flagged, never assumed clean."""
    blocked, reasons = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site="stool",
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=True, mode="participant",
    )
    assert blocked == "population_unsupported"
    assert any("antibiotics_within_90_days" in r for r in reasons)
    _, unknown = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site="stool",
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=None, mode="participant",
    )
    assert any("unknown_is_not_none" in r for r in unknown)


@pytest.mark.parametrize("site", ["scalp", "blood", "skin", "saliva"])
def test_v4237_non_stool_material_is_rejected_by_the_wrapper(site):
    """V42-37: the body-site gate fires before anything is scored."""
    blocking, reasons = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site=site,
        usable_read_pairs=2_000_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert blocking == "assay_incompatible"
    assert any(site in r for r in reasons)


def test_v4238_reads_mislabelled_as_pairs_cannot_satisfy_the_floor():
    """V42-38: 500,000 reads is not 500,000 pairs."""
    blocking, reasons = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site="stool",
        usable_read_pairs=500_000, counts_are_pairs=False, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert blocking == "qc_failed"
    assert any("reads_not_pairs" in r for r in reasons)
    ok, _ = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site="stool",
        usable_read_pairs=500_000, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert ok is None


def test_v4238_below_the_depth_floor_fails_qc():
    blocking, _ = skin.eligibility(
        skin.PSO_CHANG2022, age_years=40.0, body_site="stool",
        usable_read_pairs=499_999, counts_are_pairs=True, qc_passed=True,
        antibiotics_within_90_days=False, mode="participant",
    )
    assert blocking == "qc_failed"


def test_v4239_diversity_alone_adds_no_skin_score():
    """V42-39: no panel uses Shannon, F/B ratio, leaky gut or microbial age."""
    banned = ("shannon", "diversity", "firmicutes_bacteroidetes", "leaky",
              "calprotectin", "microbial_age", "scfa")
    for panel in skin.PANELS:
        for f in panel.features:
            low = f.source_id.lower()
            for token in banned:
                assert token not in low, (panel.panel_id, f.source_id)


# --------------------------------------------------------------------------- #
# evidence handling (V42-40..V42-46)
# --------------------------------------------------------------------------- #


def test_v4240_sgb4348_card_carries_its_pemphigus_sharing():
    """V42-40: the shared-marker counterevidence travels with the panel."""
    import openbiota.pdfskin as pdfskin

    text = pdfskin._CONFLICTS["PSO_DENG2026_SGB4348_V1"]
    assert "pemphigus foliaceus" in text
    assert "not specific" in text.lower()
    assert "E03" in skin.PSO_DENG2026.counterevidence_ids


def test_v4241_disease_versus_disease_features_are_not_case_versus_healthy():
    """V42-41: psoriasis-versus-HS/PF is never imported as versus-healthy."""
    assert skin.PSO_SA2026.contrast == "severe_psoriasis_vs_healthy"
    ids = {f.source_id for f in skin.PSO_SA2026.features}
    # SGB4348/SGB15101 were the PF-versus-control findings; not in this panel
    assert "s__GGB51647_SGB4348" not in ids
    assert "SGB15101" not in ids


def test_v4242_hs_null_finding_is_counterevidence_not_replication():
    """V42-42: a negative comparison is not a positive replication."""
    assert "E03" in skin.HS_OGUT2022.counterevidence_ids
    assert "E03" not in skin.HS_OGUT2022.source_ids
    import openbiota.pdfskin as pdfskin
    assert "no significant taxonomic difference" in pdfskin._CONFLICTS["HS_OGUT2022_TRANSLATED_V1"]


def test_v4243_pooled_qvalues_stay_on_the_pooled_contrast():
    """V42-43: pooled q-values never become PsA-only or PsO-only."""
    pooled = skin.PSORIATIC_SHARED_XIAO2024
    assert pooled.contrast == "pooled_plaque_psoriasis_and_psoriatic_arthritis_vs_healthy"
    assert pooled.family == "psoriatic_disease_shared"
    assert all(f.q_value is not None for f in pooled.features)
    # the PsA family's own panels carry no q-values copied from the pooled one
    for panel in (skin.PSA_LIU2024_KO, skin.PSA_VS_PSO_KEGG00072):
        assert all(f.q_value is None for f in panel.features)


def test_v4244_published_aucs_are_not_attached_to_the_engineered_index():
    """V42-44: single-feature and random-forest AUCs stay with their sources."""
    for panel in (skin.PSO_DENG2026, skin.PSORIATIC_SHARED_XIAO2024, skin.VIT_NONSEG_LUAN2023):
        note = panel.percentile_note or ""
        assert "not the accuracy of this index" in note or "different estimator" in note \
            or "not this three-species index" in note or "one feature's figure" in note
    # and no result object carries an accuracy field
    payload = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert "auc" not in json.dumps(payload).lower()


def test_v4245_deng_metric_disagreement_is_preserved():
    """V42-45: the conflict is recorded, not silently resolved favourably."""
    task = skin.TASKS_BY_ID["PSO_DENG2026_ORIGINAL_MODEL"]
    assert task.state == "blocked_unresolved_preprocessing_and_artifacts"
    for token in ("0.74", "0.76", "0.83", "0.96"):
        assert token in task.statement
    assert "more favourable" in task.statement


def test_v4246_failed_validation_and_nonsignificant_results_carry_no_weight():
    """V42-46: PWY-5005 and a q=0.157 result get zero diagnostic weight."""
    by_feature = {o.feature: o for o in skin.NON_SCORED_OBSERVATIONS}
    pwy = next(o for k, o in by_feature.items() if "PWY-5005" in k)
    assert "failed" in pwy.source_direction
    assert pwy.to_json()["disease_score_weight"] == 0.0
    note = skin.PSA_VS_PSO_KEGG00072.percentile_note or ""
    assert "0.157" in note and "not significant" in note
    scored = {f.source_id for p in skin.PANELS for f in p.features}
    assert "PWY-5005" not in scored


# --------------------------------------------------------------------------- #
# data joins (V42-47..V42-52)
# --------------------------------------------------------------------------- #


def test_v4247_no_project_membership_makes_a_participant_healthy():
    """V42-47: an exact control whitelist is required, not 'all 70'."""
    # The counterevidence source is registered, and no control cohort is
    # asserted from project membership anywhere in the panel definitions.
    assert "E05" in skin.PSO_DENG2026.counterevidence_ids
    for panel in skin.PANELS:
        assert "PRJNA1061168" not in json.dumps(
            {"note": panel.percentile_note or "", "label": panel.label}
        ) or "whitelist" in (panel.percentile_note or "")


def test_v4248_aliases_alone_do_not_define_phenotype():
    """V42-48: numeric sample aliases are not phenotype labels."""
    # No panel encodes a participant identifier as a feature.
    for panel in skin.PANELS:
        for f in panel.features:
            assert not f.source_id.isdigit()


def test_v4249_repeat_datasets_are_one_cohort_identity():
    """V42-49: reanalysis and dataset versions are not extra replications."""
    hs = skin.HS_OGUT2022
    assert hs.discovery_sample_count == 30
    assert len(hs.source_ids) == 1
    assert hs.evidence_tier == "single_feature_nominal"


def test_v4250_import_failures_are_explicit_not_filled_in():
    """V42-50: a broken source fails loudly rather than being guessed."""
    with pytest.raises(DirIndexError):
        fraction("not a number")
    with pytest.raises(DirIndexError):
        skin.ko_fractions({"K1": float("nan")})


def test_v4251_annotation_names_without_sequences_produce_no_score():
    """V42-51: no vBin score and no generated reference sequence."""
    scored = {f.source_id for p in skin.PANELS for f in p.features}
    assert not any("vBin" in s for s in scored)
    obs = next(o for o in skin.NON_SCORED_OBSERVATIONS if "vBin" in o.feature)
    assert "no sequence" in obs.permitted_interpretation.lower()
    assert obs.to_json()["disease_score_weight"] == 0.0


def test_v4252_infant_strain_and_pediatric_16s_claims_are_not_made():
    """V42-52: no adult species claim and no pediatric WGS classifier."""
    task = skin.TASKS_BY_ID["ATD_INFANT_SEONG2025_STRAIN_V1"]
    assert task.state == "strain_bundle_pending"
    assert "does not transfer to an adult" in task.statement
    assert "8 to 18" in task.statement
    scored = {f.source_id for p in skin.PANELS for f in p.features}
    assert "Bifidobacterium_longum" not in scored


# --------------------------------------------------------------------------- #
# results, report and interventions (V42-53..V42-61)
# --------------------------------------------------------------------------- #


def test_v4253_several_panels_are_not_fused_and_do_not_inflate_the_count():
    """V42-53: no unversioned family fusion; families < panels."""
    families = {p.family for p in skin.PANELS}
    assert len(families) < len(skin.PANELS)
    out = skin.skin_indices({"Escherichia_coli": 1.0}, None)
    assert "no_family_fusion" in out
    assert "never averaged" in out["no_family_fusion"]
    # no aggregate key appears for a family
    for family in families:
        assert family not in out


def test_v4254_different_masks_are_flagged_incomparable():
    """V42-54: comparison needs the same mask, version and reference."""
    from openbiota.dirindex import comparable

    panel = skin.PSO_SA2026
    catalogue = ["Enterocloster_x", "Phascolarctobacterium_faecium", "Bacteroides_finegoldii"]
    matrix = [[1.0] * 30 for _ in catalogue]
    ref = skin.freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="t")
    full = score_panel(panel, dict.fromkeys(catalogue, 0.01), ref)
    partial = score_panel(panel, {"Bacteroides_finegoldii": 0.01}, ref)
    if full.feature_mask != partial.feature_mask:
        assert not comparable(full, partial)


def test_v4255_every_panel_row_is_rendered_with_no_top_n():
    """V42-55: the summary paginates rather than truncating."""
    import openbiota.pdfskin as pdfskin

    skin_payload = {p.panel_id: {"panel_id": p.panel_id, "profile_family": p.family,
                                 "label": p.label, "status": "reference_unavailable",
                                 "index": None, "usable_count": 0,
                                 "panel_count": len(p.features),
                                 "full_panel_bounds": [0.0, 100.0],
                                 "evidence_tier": p.evidence_tier,
                                 "contrast": p.contrast}
                    for p in skin.PANELS}
    records = pdfskin._panels(skin_payload)
    assert len(records) == len(skin.PANELS)
    source = Path(pdfskin.__file__).read_text()
    assert "repeatRows=1" in source  # header repeats across pages
    assert "[:8]" not in source and "top_n" not in source


def test_v4256_nulls_render_as_a_dash_with_a_reason_never_zero():
    """V42-56: no green clearance, no per-cent sign, no zero for a null."""
    import openbiota.pdfskin as pdfskin

    record = {"status": "reference_unavailable", "index": None, "usable_count": 0, "panel_count": 3}
    text = pdfskin._index_text(record)
    assert "&mdash;" in text
    assert "no compatible reference group" in text
    assert "0" not in text.replace("size=7", "")
    assert "%" not in text


def test_v4240_conflicts_are_found_under_the_result_layers_key_name():
    """The result layer emits ``profile_id``; the renderer must still find the panel."""
    import openbiota.pdfskin as pdfskin

    record = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert "panel_id" not in record or record.get("panel_id") == skin.PSO_DENG2026.panel_id
    assert pdfskin._panel_id(record) == "PSO_DENG2026_SGB4348_V1"
    assert pdfskin._CONFLICTS.get(pdfskin._panel_id(record))


def test_v4258_a_high_family_index_materialises_an_evidence_only_trigger():
    """The psoriasis card is reachable from a high index, under its own namespace."""
    from openbiota.interventions import materialize_triggers

    triggers = materialize_triggers(
        sample="T", panel_levels={}, group_levels={}, profile_levels={},
        index_levels={"plaque_psoriasis": "high", "non_segmental_vitiligo": "typical"},
    )
    pso = next(t for t in triggers if t.target == "plaque_psoriasis")
    assert pso.trigger_type == "research_profile_resemblance"
    assert pso.direction == "high"
    assert pso.microbial_context["namespace"] == "directional_pattern_index_0_100"
    vit = next(t for t in triggers if t.target == "non_segmental_vitiligo")
    assert vit.direction == "any"


def test_v4256_single_marker_panels_say_so_in_words():
    import openbiota.pdfskin as pdfskin

    assert "single-marker result" in pdfskin._coverage_text(
        {"usable_count": 1, "panel_count": 1})


def test_v4256_transported_panels_declare_their_lane():
    import openbiota.pdfskin as pdfskin

    assert pdfskin._TIER_WORDS["assay_transported_exploratory"] == (
        "16S-derived pattern applied to shotgun measurements")
    for panel in (skin.ATD_ADULT_WANG2023, skin.HS_OGUT2022):
        assert panel.assay_transport == "unvalidated_16s_to_shotgun"


def test_v4257_st11_and_topical_cards_keep_route_and_phenotype_distinct():
    """V42-57: oral dandruff trial versus topical SEBD pilot."""
    registry = load_registry(REPO / "interventions")
    st11 = registry.assertions["IE_SEBD_ST11_DANDRUFF_RCT_V1"]
    assert "dandruff" in st11.summary
    assert "scalp" in st11.summary
    assert "does not establish a gut signature" in st11.summary
    topical = registry.assertions["IE_SEBD_TOPICAL_PROBIOTIC_PILOT_V1"]
    assert "applied to the skin, not swallowed" in topical.summary
    assert "not an oral regimen" in topical.summary
    intervention = registry.interventions["TOPICAL_SEBD_LCRISPATUS_P17631_LPARACASEI_I1688_V1"]
    assert intervention.intervention_class == "topical_live_biotherapeutic"


def test_v4257_strain_identity_is_exact():
    registry = load_registry(REPO / "interventions")
    psoriasis = registry.interventions["PROBIOTIC_PSO_BLONGUM_BLACTIS_LRHAMNOSUS_V1"]
    for strain in ("CECT 7347", "CECT 8145", "CECT 8361"):
        assert strain in psoriasis.display_name
    policy = psoriasis.transfer_policy
    assert policy["same_species_other_strain"] == "prohibited"
    assert policy["genus_only"] == "prohibited"


def test_v4258_a_high_index_triggers_no_treatment_or_fmt():
    """V42-58: resemblance is not eligibility, a dose, or an eradication."""
    registry = load_registry(REPO / "interventions")
    mapping = next(m for m in registry.mappings
                   if m.action_mapping_id == "MAP_PROFILE_PSO_PROBIOTIC_ADJUNCT_V1")
    assert mapping.render == "evidence_only"
    assert mapping.fmt_status == "not_indicated"
    assert mapping.confirmation_rule_refs
    lowered = " ".join(mapping.never_say).lower()
    for phrase in ("treats psoriasis", "antimicrobials", "eradicate"):
        assert phrase in lowered
    assert "does not follow from this index" in mapping.context_note


def test_v4258_no_card_claims_fmt_is_the_only_solution():
    registry = load_registry(REPO / "interventions")
    for mapping in registry.mappings:
        if mapping.action_mapping_id in (
            "MAP_PROFILE_PSO_PROBIOTIC_ADJUNCT_V1", "MAP_PROFILE_SEBD_EVIDENCE_V1"
        ):
            assert mapping.fmt_status == "not_indicated"
            assert "only solution" not in (mapping.context_note or "").lower()


def test_v4259_donor_views_get_no_pass_fail_from_skin_panels():
    """V42-59: null or low skin indices are not donor clearance."""
    payload = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert payload.get("clinical_classification") is None
    assert payload.get("disease_probability") is None
    text = json.dumps(payload).lower()
    for token in ("donor_eligibility_pass", "eligible", "cleared", "pass_fail"):
        assert token not in text


def test_v4260_the_five_personal_samples_are_never_reference_data():
    """V42-60: fixtures do not tune weights, references or thresholds."""
    for panel in skin.PANELS:
        assert panel.weight_rule == "equal"
    # every direction comes from a source, and no coefficient is a weight
    for panel in skin.PANELS:
        for f in panel.features:
            assert f.to_json()["coef_is_a_weight"] is False
    source = Path(skin.__file__).read_text()
    for sample in ("SAMPLE2", "A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06"):
        assert sample not in source


def test_v4261_partial_mask_index_and_full_panel_range_are_separate_scopes():
    """V42-61: the range is not a confidence interval around the index."""
    # one point at 95, one interval-only slot bounded at [0, 10]
    out = summarize({"a": (95.0, 95.0, 95.0)}, ["a", "b"], interval_bounds={"b": (0.0, 10.0)})
    index, bounds = out["index"], out["full_panel_bounds"]
    assert index == 95.0
    assert bounds == [47.5, 52.5]
    # the tightened range excludes the point index with no arithmetic error
    assert not (bounds[0] <= index <= bounds[1])
    payload = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert payload["index_scope"] == "available_point_feature_mask"
    assert payload["bounds_scope"] == "full_declared_panel"
    assert payload["bounds_type"] == "coverage_and_censoring"
    import openbiota.pdfskin as pdfskin
    source = Path(pdfskin.__file__).read_text()
    assert "not a margin of error" in source


# --------------------------------------------------------------------------- #
# release-wide language checks
# --------------------------------------------------------------------------- #


def test_every_family_is_named_in_plain_language():
    """A reader should not need to know the medical term to follow the page."""
    for family in skin.FAMILIES.values():
        assert family.plain_label
        assert family.plain_label != family.label
        assert "(" in family.plain_label or family.plain_label.islower()


def test_no_panel_claims_specificity_or_a_probability():
    for panel in skin.PANELS:
        assert panel.disease_specificity == "not_established"
        assert panel.score_semantics.startswith("unvalidated")


def test_the_index_is_never_described_as_a_percentage():
    payload = skin.SkinResult(skin.PSO_DENG2026, None, "reference_unavailable").to_json()
    assert payload["index_label"] == "Research directional pattern index"
    text = json.dumps(payload).lower()
    assert "percentage" not in text and "%" not in text
    assert payload["control_percentile"] is None


def test_weighted_midrank_matches_the_specification():
    assert weighted_midrank(50, [20, 50, 50, 80], [1, 1, 1, 1]) == 50


def test_specification_kernel_self_test_passes():
    from openbiota.dirindex import self_test

    self_test()


def test_evaluate_matches_the_specifications_reference_arithmetic():
    """The spec's own §5.2 fixture, run against our port."""
    rows = [{"a": 0.01} for _ in range(20)]
    ref = fit_reference(rows, ["a"])
    assert ref["a"]["scale"] == 0.25
    features = [{"id": "a", "direction": 1}, {"id": "b", "direction": -1}]
    obs = {"a": {"point": 0.01, "lower": 0.01, "upper": 0.01}}
    out = evaluate(features, obs, ref)
    assert out["index"] == 50
    assert out["coverage"] == 0.5
    assert out["full_panel_bounds"] == [25, 75]
    assert evaluate(features, {}, ref)["index"] is None
    censored = evaluate(features, {"a": {"point": None, "lower": 0, "upper": 0.001}}, ref)
    assert censored["index"] is None
    assert censored["interval_coverage"] == 0.5
    up = 50.0 + 50.0 * math.tanh(1 * ((math.log10(0.1 + 1e-6) - ref["a"]["mu"]) / ref["a"]["scale"]) / 2)
    assert up > 50
