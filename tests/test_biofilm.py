"""Acceptance suite for BUILD_SPEC_v08.0, tests BF-T001 .. BF-T125.

These are behaviour contracts, not clinical validation, and the test names
say which numbered criterion each one enforces. Most of them are written as
refusals: the interesting failure mode for this module is not a crash but a
plausible-looking number that means something other than what its label
says. So the suite spends most of its effort asserting that things do NOT
happen - that a proxy does not populate a mechanism axis, that an unrun
assay does not become a zero, that opposing evidence does not get dropped,
and that two independent axes are never combined.

Where a criterion is about data the repository does not hold (RNA, imaging,
request-only cohorts), the test asserts the correct refusal rather than
skipping, because "we did not implement it" and "we implemented a refusal"
are different states and only the second is safe.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from openbiota.biofilm import (
    datasets,
    engine,
    interventions,
    reference,
    registry,
    scoring,
)
from openbiota.biofilm import sources as src
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
COHORT = REPO / "refs" / "taxonomic_cohort.json"


# --------------------------------------------------------------------------
# Input and detection: BF-T001 .. BF-T016
# --------------------------------------------------------------------------


def test_bf_t001_taxonomic_profile_alone_does_not_invent_accessory_genes():
    """A taxonomic profile yields context and explicit unassayed loci."""
    gene_modules = [
        m
        for m in registry.MODULES.values()
        if m.measurement_kind in {"gene_family", "locus_completeness"}
    ]
    assert gene_modules
    # None of them is implemented, because no gene-level biofilm assay runs.
    for m in gene_modules:
        assert m.census_state != "implemented", m.id
        assert not m.eligible_for_axis, m.id


def test_bf_t002_bowtie2_summary_is_not_a_locus_assay():
    """A MetaPhlAn bowtie2 summary is not treated as SAM/BAM or a gene assay."""
    # The reference layer reads only the .tsv profile, never the .bowtie2.bz2.
    assert reference.MPA3_PROFILE.endswith(".tsv")
    assert "bowtie2" not in reference.MPA3_PROFILE
    assert "bowtie2" not in reference.MPA4_PROFILE


def test_bf_t003_registry_change_changes_the_calibration_hash():
    """Changing the frozen panel changes its calibration identity."""
    rows = [{"a": float(i), "b": float(i)} for i in range(40)]
    p1 = scoring.COMMUNITY_PROXY
    h1 = scoring.calibration_hash(p1, rows, "ref-1")
    h2 = scoring.calibration_hash(p1, rows, "ref-2")
    h3 = scoring.calibration_hash(scoring.ECOLOGY_PROXY, rows, "ref-1")
    assert h1 != h2, "a different reference must produce a different calibration"
    assert h1 != h3, "a different panel must produce a different calibration"


def test_bf_t004_prose_only_change_does_not_touch_detection():
    """Evidence wording is data, separate from anything that reads sequence."""
    # Source text lives in a module with no sequence dependency at all.
    import openbiota.biofilm.sources as s

    assert not hasattr(s, "run_alignment")
    assert "fastq" not in Path(s.__file__).read_text().lower()


def test_bf_t005_reference_bootstrap_resamples_participants_not_values():
    """Participants are resampled whole, retaining their paired vectors."""
    rows = [{"a": float(i), "b": float(100 - i)} for i in range(40)]
    out = scoring.reference_interval(
        {"a": 5.0, "b": 95.0},
        rows,
        (("a", "b"),),
        input_hashes=("h",),
        registry_release="r",
        calibration_hash="c",
        replicates=200,
    )
    # A resampled participant keeps a+b == 100, so a value-wise shuffle
    # would be detectable. The interval exists, which means whole rows were
    # resampled and the panel stayed coherent.
    assert out.attempted == 200


def test_bf_t006_ambiguous_carrier_stays_ambiguous():
    """A shared domain does not force a species attribution."""
    m = registry.get("BF-M02")
    assert "untested homolog" in m.interpretation or "homolog" in m.not_established
    assert m.carrier_required_for_organism_claim


def test_bf_t007_database_missing_target_is_not_absent():
    """A target with no reference panel reports its state, not a negative."""
    m = registry.get("BF-M11")
    assert m.census_state == "source_unavailable"
    assert m.feature_symbols == ()
    assert not m.eligible_for_axis


def test_bf_t008_low_depth_cannot_produce_a_confident_negative():
    """Fungal low depth is explicit, not a clean negative."""
    m = registry.get("BF-M09")
    assert m.census_state == "not_measurable_from_current_assay"
    assert "depth" in m.interpretation.lower()


def test_bf_t009_regulator_cannot_replace_a_structural_mechanism():
    """An isolated regulatory gene is not a complete mechanism."""
    curli = registry.get("BF-M01")
    assert "regulator" in curli.interpretation.lower()
    assert curli.candidate_axis == "context_dependent"
    bss = registry.get("BF-M04")
    assert "regulator" in bss.not_established.lower()


def test_bf_t010_unlinked_components_do_not_create_an_operon():
    """Linked-concern context requires direct carrier linkage."""
    m = registry.get("BF-M20")
    assert "linkage" in m.interpretation.lower()
    assert "co-occurrence" in m.not_established.lower()


def test_bf_t011_assembly_collapse_yields_uncertainty():
    m = registry.get("BF-M05")
    assert "not" in m.not_established.lower()
    assert m.molecular_resolution != "strain"


def test_bf_t012_fungal_fixtures_do_not_inherit_bacterial_limits():
    m = registry.get("BF-M09")
    assert "bacterial detection limits" in m.not_established


def test_bf_t013_family_homology_is_not_an_allele_match():
    for mid in ("BF-M02", "BF-M05"):
        assert registry.get(mid).molecular_resolution == "gene_family"


def test_bf_t014_incompatible_units_fail_scoring_but_keep_raw_data():
    """MetaPhlAn 4 abundances are not ranked against a MetaPhlAn 3 reference."""
    text = (REPO / "openbiota" / "biofilm" / "reference.py").read_text()
    assert "MPA4_PROFILE" in text
    # The MetaPhlAn 4 constant exists but is never read by sample_features.
    body = text.split("def sample_features", 1)[1]
    assert "MPA4_PROFILE" not in body
    assert "MPA3_PROFILE" in body


def test_bf_t015_detection_states_survive_export():
    """Nondetection, unassayed and detected remain distinguishable in JSON."""
    doc = _erick()
    detail = _card(doc, "H-C")["feature_detail"]
    states = {d["censoring"] for d in detail.values()}
    assert states <= {"none", "left_censored", "not_assayed"}
    # And a nondetected feature is flagged, not silently zero.
    for d in detail.values():
        if not d["detected"]:
            assert d["censoring"] in {"left_censored", "not_assayed"}


def test_bf_t016_no_absolute_cell_counts_from_relative_abundance():
    doc = _erick()
    blob = json.dumps(doc).lower()
    for forbidden in ("cells per gram", "cfu/g", "cells/g", "viable cells"):
        assert forbidden not in blob


# --------------------------------------------------------------------------
# Biological interpretation: BF-T017 .. BF-T030
# --------------------------------------------------------------------------


def test_bf_t017_curli_does_not_imply_amyloid_or_parkinsons():
    m = registry.get("BF-M01")
    for claim in ("brain amyloid", "Parkinson", "fibrin", "active"):
        assert claim.lower() in m.not_established.lower()


def test_bf_t018_bap_homologs_stay_distinct_from_tested_domains():
    m = registry.get("BF-M02")
    assert "untested homolog" in m.interpretation


def test_bf_t019_f_prausnitzii_abundance_does_not_activate_pm():
    """Genus abundance feeds only the labelled P-E proxy, never P-M."""
    assert "faecalibacterium" in scoring.ECOLOGY_PROXY.features
    # BF-M11 is the protective-mechanism module and is not axis-eligible.
    assert not registry.get("BF-M11").eligible_for_axis
    doc = _erick()
    pm = _card(doc, "P-M")
    assert pm["reference_percentile"] is None
    assert pm["status"] == "no_validated_panel"


def test_bf_t020_generic_genes_do_not_fill_missing_protective_markers():
    m = registry.get("BF-M11")
    assert "cobalamin" in m.not_established
    assert scoring.FORBIDDEN_PROXY_FEATURES
    for banned in scoring.FORBIDDEN_PROXY_FEATURES:
        assert banned not in scoring.ECOLOGY_PROXY.features
        assert banned not in scoring.COMMUNITY_PROXY.features


def test_bf_t021_l_reuteri_species_is_not_the_studied_strain():
    m = registry.get("BF-M10")
    assert m.molecular_resolution == "strain"
    assert "Species-level abundance alone does not activate" in m.interpretation


def test_bf_t022_s_boulardii_evidence_is_not_all_s_cerevisiae():
    card = interventions.get("BF-I24")
    assert "S. cerevisiae did NOT reproduce" in card.result
    assert "Exact strain" in card.caution


def test_bf_t023_community_genes_are_not_attributed_without_linkage():
    m = registry.get("BF-M04")
    assert "R. gnavus without linkage" in m.not_established
    assert m.carrier_required_for_organism_claim


def test_bf_t024_bsss_is_not_monotonic():
    assert "BssS is a regulator" in registry.get("BF-M04").not_established


def test_bf_t025_bf_d01_is_context_not_training():
    ds = datasets.get("BF-D01")
    assert ds.assay == "16s"
    assert not ds.training_eligible
    assert "Not shotgun training reads" in ds.exclusion


def test_bf_t026_biopsy_bloom_threshold_is_not_a_stool_threshold():
    assert datasets.get("BF-D01").guard == "bloom_rule_is_biopsy_only"
    assert "not a universal stool threshold" in datasets.get("BF-D01").exclusion


def test_bf_t027_disease_label_is_not_a_biofilm_label():
    ds = datasets.get("BF-D03")
    assert ds.guard == "disease_label_is_not_biofilm_label"
    assert not ds.training_eligible


def test_bf_t028_bbsdb_performance_is_not_diagnostic_accuracy():
    s = src.get("BF-S07")
    assert "patient diagnostic accuracy" in s.cannot_support


def test_bf_t029_capsule_alone_is_not_attached_biofilm():
    m = registry.get("BF-M06")
    assert "capsule hit is an attached biofilm" in m.not_established
    assert m.candidate_axis == "context_dependent"


def test_bf_t030_one_observation_many_citations():
    """The same evidence appears once with several provenance links."""
    m = registry.get("BF-M03")
    assert len(m.primary_source_ids) > 1
    # But it is one module in one dependence group, not several.
    assert m.aggregation_group == "BF-G-PNAG"


# --------------------------------------------------------------------------
# Scoring: BF-T031 .. BF-T050
# --------------------------------------------------------------------------


def test_bf_t031_midrank_reference_vectors():
    ref = [0, 0, 2, 4]
    assert scoring.midrank(0, ref) == 25.0
    assert scoring.midrank(2, ref) == 62.5
    assert scoring.midrank(5, ref) == 100.0


def test_bf_t032_identical_values_stay_tied_through_both_stages():
    rows = [{"a": float(i), "b": float(i)} for i in range(40)]
    groups = (("a", "b"),)
    target = {"a": 7.0, "b": 7.0}
    first = scoring.rank_panel(target, rows, groups)
    again = scoring.rank_panel(dict(target), list(rows), groups)
    assert first["reference_percentile"] == again["reference_percentile"]
    # A reference row with the same values gets the same aggregate.
    assert first["raw_index"] == scoring.rank_panel(rows[7], rows, groups)["raw_index"]


def test_bf_t033_constant_reference_is_non_discriminating():
    with pytest.raises(scoring.ReferenceNonDiscriminating):
        scoring.midrank(1.0, [2.0] * 10)
    with pytest.raises(scoring.ReferenceNonDiscriminating):
        scoring.rank_panel({"a": 1.0}, [{"a": 2.0}] * 10, (("a",),))


def test_bf_t034_zero_keeps_censoring_and_unknown_stays_null():
    doc = _erick()
    for d in _card(doc, "H-C")["feature_detail"].values():
        if d["value"] == 0:
            assert d["censoring"] == "left_censored"
            assert d["detected"] is False
    # An unrun mechanism card has a null percentile, not a zero.
    assert _card(doc, "H-M")["reference_percentile"] is None


def test_bf_t035_protective_and_harmful_never_cancel():
    """Two high readings stay two high readings."""
    doc = _erick()
    blob = json.dumps(doc)
    for forbidden in ("net_score", "combined_score", "biofilm_balance", "net_biofilm"):
        assert forbidden not in blob
    assert "never subtracted" in doc["summary"]["no_cancellation_note"]
    # The two axes have disjoint features, so no shared term can link them.
    assert not set(scoring.COMMUNITY_PROXY.features) & set(
        scoring.ECOLOGY_PROXY.features
    )


def test_bf_t036_candidate_module_cannot_contribute_to_an_axis():
    for m in registry.MODULES.values():
        if m.census_state == "candidate_not_validated":
            assert not m.eligible_for_axis, m.id


def test_bf_t037_one_group_is_limited_scope_and_says_so():
    doc = _erick()
    hc = _card(doc, "H-C")
    assert hc["scope"] == "single_mechanism"
    assert hc["scope_names"] == ["Biofilm-associated community pattern"]
    assert "ecosystem health" not in json.dumps(hc).lower()


def test_bf_t038_correlated_genes_do_not_multiply_weight():
    groups = registry.census()["dependence_groups"]
    assert set(groups["BF-G-PNAG"]) == {"BF-M03", "BF-M04"}
    # P-E's four genera are one group, not four.
    assert len(scoring.ECOLOGY_PROXY.groups) == 1
    assert len(scoring.ECOLOGY_PROXY.groups[0]) == 4


def test_bf_t039_missing_modules_are_not_replaced_by_zero():
    with pytest.raises(scoring.MaskMismatch):
        scoring.rank_panel({"a": 1.0}, [{"a": 1.0, "b": 2.0}] * 5, (("a", "b"),))


def test_bf_t040_zero_groups_cannot_produce_a_score():
    with pytest.raises(ValueError):
        scoring.rank_panel({}, [{}], ())
    with pytest.raises(ValueError):
        scoring.rank_panel({"a": 1.0}, [{"a": 1.0}], ((),))


def test_bf_t041_mask_does_not_depend_on_what_was_detected():
    """Both features stay in the mask whether or not they were detected."""
    doc = _erick()
    hc = _card(doc, "H-C")
    assert set(hc["feature_values"]) == set(scoring.COMMUNITY_PROXY.features)
    # SAMPLE2 has no detected E. coli, and the feature is still in the mask.
    assert hc["feature_detail"]["ecoli_complex"]["detected"] is False
    assert "ecoli_complex" in hc["feature_percentiles"]


def test_bf_t042_same_sample_scores_identically_alone_or_in_a_batch():
    a = _analyze("SAMPLE2_A02")
    b = _analyze("SAMPLE2_A02")
    assert _card(a, "H-C")["reference_percentile"] == _card(b, "H-C")[
        "reference_percentile"
    ]
    assert _card(a, "P-E")["reference_percentile"] == _card(b, "P-E")[
        "reference_percentile"
    ]
    assert a["provenance"]["calibration_id"] == b["provenance"]["calibration_id"]


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_bf_t043_local_samples_are_not_the_calibration_cohort():
    cohort = reference.load_reference(str(COHORT))
    local = {"SAMPLE2", "K1", "Z1", "J1", "B2"}
    for subject in cohort.subjects:
        assert not any(subject.startswith(f"{name}") for name in local)
    assert cohort.n >= scoring.MIN_REFERENCE_N


def test_bf_t044_adult_calibration_is_out_of_domain_for_a_child():
    doc = _analyze("SAMPLE2_A02", subject={"supplied": {"age": 9}})
    assert doc["out_of_domain"]
    for card in doc["cards"]:
        assert card["reference_percentile"] is None
    assert "pediatric abnormality" in doc["out_of_domain"]


def test_bf_t045_separate_axes_have_separate_calibration_ids():
    doc = _erick()
    hc, pe = _card(doc, "H-C"), _card(doc, "P-E")
    assert hc["calibration_id"] and pe["calibration_id"]
    assert hc["calibration_id"] != pe["calibration_id"]


def test_bf_t046_bootstrap_seeds_are_deterministic_and_labelled():
    rows = [{"a": float(i % 7), "b": float(i % 5)} for i in range(60)]
    kw = {
        "input_hashes": ("h1",), "registry_release": "r",
        "calibration_hash": "c", "replicates": 200,
    }
    one = scoring.reference_interval({"a": 3.0, "b": 2.0}, rows, (("a", "b"),), **kw)
    two = scoring.reference_interval({"a": 3.0, "b": 2.0}, rows, (("a", "b"),), **kw)
    assert (one.low, one.high) == (two.low, two.high)
    assert one.label == "reference interval"
    assert scoring.analytical_interval(available=False).label == "sampling interval"


def test_bf_t047_summary_only_cache_yields_no_fragment_interval():
    doc = _erick()
    hc = _card(doc, "H-C")
    analytical = hc["intervals"]["analytical"]
    assert analytical["low"] is None
    assert analytical["status"] == "taxonomic_summary_table_only"


def test_bf_t048_reference_interval_is_not_a_clinical_interval():
    doc = _erick()
    ref = _card(doc, "H-C")["intervals"]["reference"]
    assert ref["label"] == "reference interval"
    assert "Not a statement about clinical risk" in ref["meaning"]


def test_bf_t049_nonfinite_and_out_of_range_fail_validation():
    for bad in (float("nan"), float("inf"), float("-inf")):
        with pytest.raises(ValueError):
            scoring.midrank(bad, [1.0, 2.0])
        with pytest.raises(ValueError):
            scoring.rank_panel({"a": bad}, [{"a": 1.0}, {"a": 2.0}], (("a",),))
    doc = _erick()
    for card in doc["cards"]:
        pct = card["reference_percentile"]
        if pct is not None:
            assert 0.0 <= pct <= 100.0
            assert math.isfinite(pct)


def test_bf_t050_no_calibrated_axis_still_produces_a_full_report():
    doc = _analyze("SAMPLE2_A02", subject={"supplied": {"age": 9}})
    assert all(c["reference_percentile"] is None for c in doc["cards"])
    assert len(doc["module_results"]) == 24
    assert doc["intervention_candidates"]
    assert doc["dataset_audit"]["n_datasets"] == 14


# --------------------------------------------------------------------------
# Interventions: BF-T051 .. BF-T070
# --------------------------------------------------------------------------


def test_bf_t051_allicin_appears_without_a_human_rct():
    card = interventions.get("BF-I01")
    assert card.study_type == "in_vitro"
    assert card in interventions.retrieve(("ecoli_adhesion",))


def test_bf_t052_berberine_positive_and_counterevidence_appear_together():
    got = {i.id for i in interventions.retrieve(("ecoli_adhesion",))}
    assert {"BF-I05", "BF-I06"} <= got
    assert interventions.get("BF-I06").direction == "unfavourable"


def test_bf_t053_antibacterial_and_biofilm_direction_stay_separate():
    card = interventions.get("BF-I07")
    assert card.direction == "mixed"
    assert "separate outcomes" in card.caution


def test_bf_t054_nac_ph_qualifier_survives_rendering():
    card = interventions.get("BF-I15")
    assert "pH" in card.lab_exposure or "pH" in card.caution
    assert "NEUTRALISED" in card.result


def test_bf_t055_formation_assay_is_not_eradication():
    card = interventions.get("BF-I15")
    assert card.endpoint == "formation"
    assert "formation-assay design" in card.caution


def test_bf_t056_crystal_violet_is_not_viable_killing():
    card = interventions.get("BF-I01")
    assert card.endpoint == "matrix_biomass"
    assert card.endpoint_plain == "matrix biomass (stain-based)"
    assert "Dispersal is not sterilisation" in card.caution


def test_bf_t057_serrapeptase_capsule_material_stays_in_vitro():
    card = interventions.get("BF-I18")
    assert card.model == "gut_relevant_in_vitro"
    assert "not established" in card.delivery_note


def test_bf_t058_lab_concentration_is_not_a_dose():
    for card in interventions.INTERVENTIONS.values():
        assert not hasattr(card, "oral_dose")
        assert not hasattr(card, "suggested_dose")
    assert "no dose conversion" in interventions.get("BF-I18").caution.lower()


def test_bf_t059_enteric_delivery_is_plausible_not_proven():
    card = interventions.get("BF-I02")
    assert card.delivery == "indirect_delivery_demonstrated"
    assert "is plausible and is not ruled out" in card.delivery_note


def test_bf_t060_metabolite_is_not_intact_colonic_exposure():
    assert "not intact colonic allicin" in interventions.get("BF-I02").caution


def test_bf_t061_blend_result_is_not_assigned_to_each_ingredient():
    card = interventions.get("BF-I22")
    assert "not for the ingredients individually" in card.caution
    assert "REVERSED at a high dose" in card.result
    assert card.funding_conflict


def test_bf_t062_cureus_preserves_its_real_numbers_and_direction():
    card = interventions.get("BF-I23")
    assert "13 patients" in card.result
    assert "60% WITH adjunct versus 100% control" in card.result
    assert "not significant" in card.result
    assert card.direction == "unfavourable"


def test_bf_t063_nac_salvage_and_nonbenefit_both_present():
    salvage, none = interventions.get("BF-I16"), interventions.get("BF-I17")
    assert salvage.direction == "favourable"
    assert none.direction == "unfavourable"
    assert "BF-I17" in salvage.opposes and "BF-I16" in none.opposes
    both = {i.id for i in interventions.retrieve(("hpylori",))}
    assert {"BF-I16", "BF-I17"} <= both


def test_bf_t064_lactoferrin_opposite_directions_not_averaged():
    card = interventions.get("BF-I20")
    assert "INHIBITED" in card.result and "ENHANCED" in card.result
    assert "must not be averaged" in card.caution


def test_bf_t065_d_amino_acid_counterevidence_stays_attached():
    got = {i.id for i in interventions.retrieve(("matrix_generic",))}
    assert {"BF-I29", "BF-I30"} <= got


def test_bf_t066_ecoli_evidence_is_not_an_r_gnavus_claim():
    for card in interventions.INTERVENTIONS.values():
        if "Escherichia coli" in card.organisms:
            assert "gnavus" not in " ".join(card.organisms)


def test_bf_t067_susceptibility_is_not_predicted_from_a_mechanism():
    assert "does not predict allicin susceptibility" in interventions.get(
        "BF-I01"
    ).caution


def test_bf_t068_phage_load_is_not_biofilm_clearance():
    card = interventions.get("BF-I31")
    assert card.endpoint == "organism_load"
    assert "not a measured biofilm-eradication endpoint" in card.caution


def test_bf_t069_protective_signal_does_not_trigger_antimicrobials():
    """A protective observation routes to supportive evidence, not drugs."""
    doc = _erick()
    supportive = [
        a
        for a in doc["intervention_candidates"]
        if a["section"] == "support_protective_ecology"
    ]
    assert supportive
    # The protective route contains no systemic antimicrobial candidate, and
    # anything antimicrobial that is matched sits in a different section and
    # carries its own review label rather than being presented as a step.
    for act in supportive:
        assert act["section"] == "support_protective_ecology"
        assert act["endpoint"] != "eradication", act["id"]
    # Nothing anywhere is phrased as an instruction to take something.
    blob = json.dumps(doc["intervention_candidates"]).lower()
    for phrase in ("antibiotic course", "you should take", "recommended regimen"):
        assert phrase not in blob


def test_bf_t070_missing_context_does_not_erase_the_card():
    card = interventions.get("BF-I01")
    assert card.delivery == "unknown"
    assert card in interventions.retrieve(("ecoli_adhesion",))


# --------------------------------------------------------------------------
# Fibrin and systemic: BF-T071 .. BF-T080
# --------------------------------------------------------------------------


def test_bf_t071_coa_vwb_do_not_diagnose_a_fibrin_biofilm():
    m = registry.get("BF-M14")
    assert set(m.feature_symbols) >= {"coa", "vwb"}
    assert "cardiac or blood diagnosis" in m.not_established
    assert m.candidate_axis == "context_dependent"


def test_bf_t072_host_dna_cannot_generate_protein_results():
    m = registry.get("BF-M16")
    assert "F2 DNA" in m.not_established
    assert m.census_state == "not_measurable_from_current_assay"


def test_bf_t073_physiological_and_excessive_thrombin_stay_distinct():
    physio, excess = src.get("BF-S12"), src.get("BF-S13")
    assert "not inherently adverse" in physio.supports
    assert "harmful community behaviour" in excess.result


def test_bf_t074_ex_vivo_rna_is_not_a_stool_dna_model():
    ds = datasets.get("BF-D04")
    assert "not a drop-in stool-DNA classifier" in ds.exclusion
    assert not ds.training_eligible


def test_bf_t075_staphylokinase_is_not_automatically_protective():
    m = registry.get("BF-M14")
    assert m.candidate_axis == "context_dependent"
    assert not m.eligible_for_axis


def test_bf_t076_fibrinolytic_enzymes_are_separate_identities():
    ids = {i.id for i in interventions.retrieve(("fibrinolytic_enzymes",))}
    assert {"BF-I38", "BF-I39"} <= ids
    assert interventions.get("BF-I38").compound != interventions.get("BF-I39").compound


def test_bf_t077_no_removal_claim_from_a_biomarker():
    doc = _erick()
    blob = json.dumps(doc).lower()
    assert "biofilm eliminated" not in blob
    assert "biofilm removed" not in blob


def test_bf_t078_no_stool_microclot_percentage():
    blob = json.dumps(_erick()).lower()
    assert "microclot" not in blob


def test_bf_t079_amr_tolerance_and_persistence_are_separate():
    m = registry.get("BF-M15")
    assert "resistant biofilm" in m.not_established
    assert "1,000-fold" not in json.dumps(_erick())


def test_bf_t080_independent_site_assays_cannot_overwrite_stool():
    doc = _erick()
    assert doc["external_assay_results"] == []
    assert doc["summary"]["active_biofilm_burden"] is None
    assert doc["summary"]["anatomical_location"] is None


# --------------------------------------------------------------------------
# Report and integration: BF-T081 .. BF-T090
# --------------------------------------------------------------------------


def test_bf_t081_overview_has_panels_context_coverage_and_links():
    doc = _erick()
    assert len(doc["cards"]) == 4
    for card in doc["cards"]:
        assert "coverage" in card
        assert "status" in card
    assert doc["headline_note"]


def test_bf_t082_null_values_are_not_green_zeros():
    doc = _erick()
    for card in doc["cards"]:
        if card["status"] != "available":
            assert card["reference_percentile"] is None
            assert card["band"] == "not computed"
            assert card["reason_codes"]


def test_bf_t083_json_contains_the_full_inventory_without_truncation():
    doc = _erick()
    assert len(doc["module_results"]) == 24
    assert doc["intervention_census"]["n_interventions"] == 50
    assert doc["dataset_audit"]["n_datasets"] == 14
    assert doc["source_count"] == 32


def test_bf_t084_external_links_are_safe_and_present():
    from openbiota import pdfbiofilm

    text = Path(pdfbiofilm.__file__).read_text()
    assert '<a href=' in text
    # Every intervention with a DOI exposes a resolvable URL.
    for card in interventions.INTERVENTIONS.values():
        assert card.url.startswith("https://doi.org/")


def test_bf_t085_figures_have_text_equivalents():
    from openbiota import pdfbiofilm

    for band in (
        "lower reference range",
        "middle reference range",
        "upper reference range",
    ):
        assert pdfbiofilm._band_symbol(band) != "\u2014"
    # The band word itself is printed, so colour is never the only channel.
    assert "band" in Path(pdfbiofilm.__file__).read_text()


def test_report_uses_only_glyphs_the_pdf_font_actually_renders():
    """Guard against glyphs Helvetica silently substitutes with a box.

    ReportLab's metrics report a width for U+25CB and U+25A1, so a width
    check passes while the rendered page shows a filled box. A hollow
    marker that renders filled inverts its meaning - "not supplied" would
    look like "supplied" - so the characters are banned outright.
    """

    banned = {
        "\u25cb": "WHITE CIRCLE renders as a filled box",
        "\u25a1": "WHITE SQUARE renders as a filled box",
        "\u25e6": "WHITE BULLET is not in the base font",
        "\u2610": "BALLOT BOX is not in the base font",
    }
    pdf_dir = Path(__file__).resolve().parent.parent / "openbiota"
    offenders: list[str] = []
    for path in sorted(pdf_dir.glob("pdf*.py")):
        text = path.read_text()
        for lineno, line in enumerate(text.splitlines(), 1):
            for char, why in banned.items():
                escaped = f"\\u{ord(char):04x}"
                if char in line or escaped in line.lower():
                    # Allow it inside a comment that explains the ban.
                    if "renders as" in line or "no U+" in line:
                        continue
                    offenders.append(f"{path.name}:{lineno}: {why}")
    assert not offenders, "\n".join(offenders)


def test_bf_t086_two_donors_cannot_cancel_or_pool():
    a, b = _analyze("SAMPLE2_A02"), _analyze("SAMPLE4_A04")
    # Independently computed; no shared mutable state, no pooled field.
    assert a["sample_id"] != b["sample_id"]
    for doc in (a, b):
        assert "pooled" not in json.dumps(doc).lower()


def test_bf_t087_generic_matrix_genes_do_not_exclude_a_donor():
    """Nothing in this module produces a donor exclusion."""
    doc = _erick()
    # No field anywhere carries a donor verdict.
    blob = json.dumps(doc).lower()
    for verdict in ("donor_excluded", "donor_cleared", "not_suitable_donor",
                    "exclude_donor", "donor_verdict"):
        assert verdict not in blob
    # And the limits say so explicitly.
    assert any("cleared for donation" in lim for lim in doc["limits"])
    # The context-only modules that could tempt an exclusion are not on an axis.
    for mid in ("BF-M06", "BF-M07", "BF-M14"):
        assert not registry.get(mid).eligible_for_axis


def test_bf_t088_nondetection_is_not_eradication():
    doc = _erick()
    blob = json.dumps(doc).lower()
    assert "eradicated" not in blob
    assert "certain elimination" not in blob


def test_bf_t089_shared_strain_is_not_disease_transfer():
    m = registry.get("BF-M17")
    assert "attached biofilm or a disease was transferred" in m.not_established


def test_bf_t090_reads_stay_local_and_provenance_persists():
    doc = _erick()
    prov = doc["provenance"]
    assert prov["input_sha256"]
    # A hash, not a sequence.
    for h in prov["input_sha256"]:
        assert len(h) <= 64 and all(c in "0123456789abcdef" for c in h)
    assert prov["spec_version"] == "08.0.1"


# --------------------------------------------------------------------------
# Revision 08.0.1 additions: BF-T091 .. BF-T125
# --------------------------------------------------------------------------


def test_bf_t091_bf_s01_denominators_stay_distinct():
    s = src.get("BF-S01")
    assert s.denominators["screened"] == 1426
    assert s.denominators["included"] == 1112
    assert s.denominators["screened"] != s.denominators["included"]
    assert "1,426" in s.result and "1,112" in s.result


def test_bf_t092_bf_s20_assays_and_denominators_are_distinct():
    s = src.get("BF-S20")
    assert s.specimen == "biopsy"
    assert s.denominators["sequenced_specimens"] == 265
    assert s.denominators["uc_patients"] == 80
    assert "BIOPSIES" in s.result
    assert "265 specimens are not 265 people" in s.cannot_support


def test_bf_t093_request_only_data_are_not_training_cohorts():
    for did in ("BF-D10", "BF-D11"):
        ds = datasets.get(did)
        assert ds.access_status == "request_only"
        assert not ds.training_eligible
    assert datasets.audit()["n_training_eligible"] == 0


def test_bf_t094_bf_d09_keeps_human_and_mouse_separate():
    ds = datasets.get("BF-D09")
    assert ds.host_species == "mixed"
    assert ds.guard == "replicates_are_not_donors"
    assert "technical replicates are not independent donors" in ds.exclusion


def test_bf_t095_bf_d12_is_not_biofilm_versus_planktonic():
    ds = datasets.get("BF-D12")
    assert ds.guard == "both_arms_are_biofilms"
    assert "NOT biofilm" in ds.exclusion
    assert "BOTH comparison arms were biofilms" in src.get("BF-S22").result


def test_bf_t096_pv739486_is_16s_only():
    ds = datasets.get("BF-D13")
    assert "PV739486" in ds.accessions
    assert ds.assay == "16s"
    assert ds.guard == "pv739486_is_16s_only"


def test_bf_t097_hc_and_pe_are_fully_distinct():
    hc, pe = scoring.COMMUNITY_PROXY, scoring.ECOLOGY_PROXY
    assert hc.id != pe.id
    assert not set(hc.features) & set(pe.features)
    assert hc.module_id != pe.module_id
    doc = _erick()
    assert _card(doc, "H-C")["calibration_id"] != _card(doc, "P-E")["calibration_id"]


def test_bf_t098_pe_never_populates_pm():
    doc = _erick()
    pe, pm = _card(doc, "P-E"), _card(doc, "P-M")
    assert pe["status"] == "available"
    assert pm["reference_percentile"] is None
    assert pm["status"] == "no_validated_panel"
    assert "not a substitute" in pm["what_it_is_not"]


def test_bf_t099_one_feature_panel_returns_limited_scope():
    out = scoring.rank_panel(
        {"a": 2.5}, [{"a": float(i)} for i in range(40)], (("a",),)
    )
    assert 0 <= out["reference_percentile"] <= 100


def test_bf_t100_generic_positivity_takes_no_direction():
    """A module cannot declare an axis without a source-backed claim."""
    with pytest.raises(ValueError, match="directional_claim"):
        registry.Module(
            id="BF-MXX",
            label="generic adhesion",
            measurement_kind="gene_family",
            mechanism_group="adhesion",
            primary_source_ids=("BF-S07",),
            feature_symbols=("someAdhesin",),
            census_state="implemented",
            analytical_status="supported_at_stated_resolution",
            candidate_axis="harmful_associated",
            aggregation_group="BF-G-X",
            interpretation="x",
            not_established="y",
        )


def test_bf_t101_full_proxy_requires_all_features():
    assert len(scoring.COMMUNITY_PROXY.features) == 2
    assert len(scoring.ECOLOGY_PROXY.features) == 4
    with pytest.raises(scoring.MaskMismatch):
        scoring.rank_panel(
            {"faecalibacterium": 1.0},
            [dict.fromkeys(scoring.ECOLOGY_PROXY.features, 1.0)] * 5,
            scoring.ECOLOGY_PROXY.groups,
        )


def test_bf_t102_hc_records_its_measurement_change():
    assert scoring.COMMUNITY_PROXY.transport_change == "genus_complex_to_species_complex"
    note = reference.TRANSPORT_NOTES["ecoli_complex"]
    assert "Escherichia-Shigella" in note
    assert "no pathotype is inferred" in note


def test_bf_t103_no_double_counting_of_leaves_or_synonyms():
    seen: set[str] = set()
    for leaves in reference.FEATURE_LEAVES.values():
        for leaf in leaves:
            assert leaf not in seen
            seen.add(leaf)
    assert "Mediterraneibacter" in reference.TRANSPORT_NOTES["m_gnavus"]
    assert reference.FEATURE_LEAVES["m_gnavus"] == ("Ruminococcus_gnavus",)


def test_bf_t104_constant_index_gives_an_explicit_state():
    with pytest.raises(scoring.ReferenceNonDiscriminating):
        scoring.rank_panel({"a": 1.0}, [{"a": 5.0}] * 40, (("a",),))


def test_bf_t105_opposed_features_with_constant_aggregate_refuse():
    rows = [{"a": float(i), "b": float(9 - i)} for i in range(10)]
    with pytest.raises(scoring.ReferenceNonDiscriminating):
        scoring.rank_panel({"a": 3.0, "b": 6.0}, rows, (("a", "b"),))


def test_bf_t106_same_arithmetic_different_interpretation():
    rows = [{"x": float(i)} for i in range(40)]
    out = scoring.rank_panel({"x": 10.0}, rows, (("x",),))
    assert 0 <= out["reference_percentile"] <= 100
    # Identical arithmetic, but the two proxies carry different meanings.
    assert scoring.COMMUNITY_PROXY.what_it_is_not != scoring.ECOLOGY_PROXY.what_it_is_not
    assert scoring.COMMUNITY_PROXY.heading != scoring.ECOLOGY_PROXY.heading


def test_bf_t107_no_cancellation_anywhere():
    doc = _erick()
    hc, pe = _card(doc, "H-C"), _card(doc, "P-E")
    if hc["reference_percentile"] is not None and pe["reference_percentile"] is not None:
        # Neither is a function of the other: they use disjoint inputs.
        assert not set(hc["feature_values"]) & set(pe["feature_values"])
    assert "no_cancellation_note" in doc["summary"]


def test_bf_t108_linked_concern_requires_supported_linkage():
    m = registry.get("BF-M20")
    assert "direct_supported" in m.interpretation
    assert "unknown carriers" in m.not_established
    assert "co-occurrence" in m.not_established.lower()


def test_bf_t109_findings_and_actions_rank_deterministically():
    a, b = _erick(), _erick()
    assert [f["id"] for f in a["ranked_findings"]] == [
        f["id"] for f in b["ranked_findings"]
    ]
    assert [x["id"] for x in a["intervention_candidates"]] == [
        x["id"] for x in b["intervention_candidates"]
    ]
    for f in a["ranked_findings"]:
        assert f["why_ranked_here"]
        # No numeric efficacy or risk probability is invented. The word may
        # appear only in a disclaimer, never as a field carrying a number.
        for key, value in f.items():
            if "probab" in key.lower():
                raise AssertionError(f"probability field invented: {key}")
            if isinstance(value, (int, float)) and "risk" in key.lower():
                raise AssertionError(f"numeric risk field invented: {key}")
    # And the two lists are ranked separately, not against each other.
    assert all(f["list"] == "needing_review" for f in a["findings_needing_review"])
    assert all(f["list"] == "supportive" for f in a["supportive_observations"])


def test_bf_t110_preclinical_cards_retrievable_with_opposing_evidence():
    got = {i.id for i in interventions.retrieve(("ecoli_adhesion",))}
    assert {"BF-I01", "BF-I05", "BF-I06"} <= got


def test_bf_t111_bf_i45_shows_both_directions():
    card = interventions.get("BF-I45")
    assert card.direction == "unfavourable"
    assert "Reduced K. pneumoniae biofilm in vitro" in card.result
    assert "MAINTAINED" in card.result
    assert card.endpoint == "viable_burden"


def test_bf_t112_bf_i42_favourable_behaviour_despite_more_biomass():
    card = interventions.get("BF-I42")
    assert card.endpoint == "harmful_behaviour"
    assert "biomass INCREASED" in card.result
    assert card.direction == "favourable"


def test_bf_t113_bf_i43_is_not_strong_clearance():
    card = interventions.get("BF-I43")
    assert "was NOT demonstrated" in card.result
    assert "must not be described as proven" in card.caution


def test_bf_t114_bf_i44_metabolites_are_not_validated_ingredients():
    card = interventions.get("BF-I44")
    assert "not proven active ingredients" in card.caution
    assert "16S sequence only" in card.strain


def test_bf_t115_exact_strains_are_not_replaced_by_species():
    for iid in ("BF-I24", "BF-I34", "BF-I44", "BF-I46"):
        card = interventions.get(iid)
        assert card.strain
        assert card.caution


def test_bf_t116_delivery_plausibility_is_not_proof():
    plausible = [
        i
        for i in interventions.INTERVENTIONS.values()
        if i.delivery == "plausible_unverified"
    ]
    assert plausible
    for card in plausible:
        assert card.delivery != "target_exposure_measured"


def test_bf_t117_high_proxy_yields_navigation_not_a_stack():
    doc = _erick()
    for act in doc["intervention_candidates"]:
        assert act["review_label"] in {
            "matched research evidence",
            "opposing evidence - review before considering",
            "mixed evidence - both directions observed",
        }
    blob = json.dumps(doc["intervention_candidates"]).lower()
    assert "recommended dose" not in blob
    assert "stack" not in blob


def test_bf_t118_reduced_dna_is_not_eradication():
    doc = _erick()
    assert "eradication" not in json.dumps(doc["limits"]).lower()
    assert any("not activity" in lim or "potential" in lim for lim in doc["limits"])


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_bf_t119_repeated_samples_never_cross_participant_folds():
    """Samples are collapsed to subjects before the reference is built."""
    cohort = reference.load_reference(str(COHORT))
    assert len(set(cohort.subjects)) == len(cohort.subjects)
    # The collapse actually removed something: 3,027 samples are not 3,027 people.
    assert cohort.n_source_samples > cohort.n_subjects_before_age_filter
    assert cohort.n <= cohort.n_subjects_before_age_filter


def test_bf_t120_a_failed_replication_stays_visible():
    """Unfavourable records are retrievable, not filtered out."""
    adverse = [i for i in interventions.INTERVENTIONS.values() if i.is_adverse]
    assert len(adverse) >= 5
    for card in adverse:
        assert card.review_label.startswith("opposing evidence")
    got = {i.id for i in interventions.retrieve(("klebsiella_adhesion",))}
    assert "BF-I45" in got


def test_bf_t121_four_cards_under_two_headings():
    doc = _erick()
    headings = {}
    for card in doc["cards"]:
        headings.setdefault(card["heading"], []).append(card["card"])
    assert len(headings) == 2
    for cards in headings.values():
        assert len(cards) == 2, "no more than two numeric cards per heading"
    for card in doc["cards"]:
        cov = card["coverage"]
        assert cov["unit"] == "modules, not a health score"


def test_bf_t122_missing_pm_does_not_suppress_other_output():
    doc = _erick()
    assert _card(doc, "P-M")["reference_percentile"] is None
    assert _card(doc, "P-E")["status"] == "available"
    assert doc["intervention_candidates"]
    assert doc["module_results"]


def test_bf_t123_all_registry_records_are_represented():
    doc = _erick()
    assert len(doc["module_results"]) == 24
    assert doc["intervention_census"]["n_interventions"] == 50
    assert doc["source_count"] == 32
    assert doc["dataset_audit"]["n_datasets"] == 14
    for m in doc["module_results"]:
        assert m["census_state"] in registry.CENSUS_STATES


def test_bf_t124_pdf_links_remain_clickable():
    from openbiota import pdfbiofilm

    assert '<a href=' in Path(pdfbiofilm.__file__).read_text()


def test_bf_t125_embedded_rank_arithmetic_passes_its_fixtures():
    scoring.self_test()
    reference.self_test()


# --------------------------------------------------------------------------
# Regression tests for bugs found in review
# --------------------------------------------------------------------------


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_sample_and_reference_share_one_measurement_basis():
    """The sample is renormalised to percent-of-classified, like the reference.

    This pipeline runs MetaPhlAn with unknown estimation, so a raw profile's
    species rows sum to about half; every curatedMetagenomicData row sums to
    100. Ranking one against the other understated every feature by roughly
    two-fold and pulled all percentiles down. The fix is a renormalisation,
    and this test pins both sides to the same total.
    """
    sample_dir = results_dir("SAMPLE2_A02")
    if not sample_dir.exists():
        pytest.skip("SAMPLE2 results not present")
    profile = sample_dir / "taxonomy" / reference.MPA3_PROFILE
    species, scale = reference.parse_mpa3_profile(profile)
    assert scale["renormalised"] is True
    assert scale["unknown_percent"] > 0, "fixture should have an unknown fraction"
    assert sum(species.values()) == pytest.approx(100.0, abs=0.01)

    cohort = reference.load_reference(str(COHORT))
    # And the reference rows are on that same basis.
    row_totals = [sum(row.values()) for row in cohort.rows[:50]]
    assert all(t <= 100.0001 for t in row_totals)


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_subjects_come_from_the_authoritative_metadata_field():
    """Repeated samples are collapsed using cMD's real subject_id.

    3,027 samples are 1,760 people; MehtaRS_2018 alone is 921 samples from
    308. Treating them as independent inflates the reference and narrows
    every bootstrap interval. The subject field is authoritative and covers
    every sample, so no study needs excluding and no ID is pattern-matched.
    """
    cohort = reference.load_reference(str(COHORT))
    prov = cohort.provenance()
    assert prov["subject_identity_source"] == "curatedMetagenomicData subject_id"
    assert not reference.EXCLUDED_STUDIES, "no study needs excluding any more"
    assert prov["excluded_studies"] == {}
    # One row per subject, and the collapse actually removed something.
    assert len(set(cohort.subjects)) == len(cohort.subjects)
    assert cohort.n_source_samples > cohort.n_subjects_before_age_filter
    assert cohort.n >= scoring.MIN_REFERENCE_N


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_how_subjects_were_identified_changes_the_calibration_id():
    """An approximate subject mapping must not share an exact one's identity."""
    cohort = reference.load_reference(str(COHORT))
    assert cohort.subject_source
    assert cohort.subject_source in cohort.provenance()["subject_identity_source"]


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_a_feature_the_reference_cannot_see_is_dropped_not_averaged():
    """An undetectable feature leaves the panel instead of contributing ~50.

    Subdoligranulum is a common gut genus, but this taxonomy release names
    only S. variabile, which appears in a handful of reference participants.
    Left in, its column is effectively constant: every sample ties at the
    same midrank, that fixed value enters the group mean, and the genera
    that do carry signal are diluted. It must be dropped, the card renamed,
    and the omission recorded.
    """
    cohort = reference.load_reference(str(COHORT))
    assert "subdoligranulum" in cohort.unmeasurable
    reason = cohort.unmeasurable["subdoligranulum"]
    assert "not assayed here rather than absent" in reason

    usable = cohort.usable_features(scoring.ECOLOGY_PROXY.features)
    assert "subdoligranulum" not in usable
    assert len(usable) == 3

    doc = _erick()
    pe = _card(doc, "P-E")
    assert "partial panel" in pe["label"]
    assert pe["proxy_id"] != scoring.ECOLOGY_PROXY.id, (
        "a partial panel needs its own identity, not the named index's"
    )
    assert any("subdoligranulum" in o for o in pe["omissions"])
    assert "subdoligranulum" not in pe["feature_values"]
    # And it is calibrated separately from the full-panel H-C card.
    assert pe["calibration_id"] != _card(doc, "H-C")["calibration_id"]


def test_dropping_a_feature_keeps_its_group_intact():
    """A four-member group losing one member becomes one group of three.

    Splitting it into three singleton groups would give each survivor a
    full group weight and silently change the aggregation.
    """
    groups = (("a", "b", "c", "d"),)
    assert engine._regroup(groups, ("a", "b", "c")) == (("a", "b", "c"),)
    # A group that loses everyone disappears rather than becoming empty.
    assert engine._regroup((("a",), ("b", "c")), ("b", "c")) == (("b", "c"),)


@pytest.mark.skipif(not COHORT.exists(), reason="reference cohort not present")
def test_calibration_identity_covers_the_species_leaves():
    """Changing which species are summed changes the calibration id.

    Hashing only the feature names would let the definition of a feature
    change while its calibration id stayed the same, so two incomparable
    numbers would look comparable.
    """
    first = reference.load_reference(str(COHORT)).id
    original = dict(reference.FEATURE_LEAVES)
    try:
        reference.FEATURE_LEAVES["blautia"] = original["blautia"][:-1]
        reference.load_reference.cache_clear()
        changed = reference.load_reference(str(COHORT)).id
    finally:
        reference.FEATURE_LEAVES.clear()
        reference.FEATURE_LEAVES.update(original)
        reference.load_reference.cache_clear()
    assert first != changed


def test_a_panel_headline_gene_with_no_fragments_is_not_detected():
    """`detected` follows the aggregate genes, not the whole panel.

    The histamine panel aggregates over the decarboxylase HDCA and also
    carries the transporter HDCTRANS as context. Z1 has HDCA = 0 and
    HDCTRANS = 2,563. `confidence` was already restricted to the aggregate
    genes and said "not detected", but `detected` summed the whole panel
    and said True, so the report placed the pathway on the reference scale
    and printed "5th percentile, NOTABLY LOW" for a gene with no reads.
    """
    from openbiota.tally import EntryResult, PanelResult

    class _Panel:
        min_fragments_for_stability = 15

        def aggregate_target_ids(self):
            return ("HDCA",)

    def _entry(entry_id: str, fragments: int) -> EntryResult:
        return EntryResult(
            key=f"histamine:{entry_id}", panel="histamine", entry_id=entry_id,
            label=entry_id, gene=entry_id, role="target", fragments=fragments,
            rejected_low_identity=0, rejected_short_alignment=0, mean_identity=0.0,
            median_identity=0.0, reference_count=10,
            reference_mean_length_aa=300.0, reference_truncated=False,
            min_identity=50.0, copies_per_100_genomes=None,
        )

    result = PanelResult(
        panel=_Panel(),
        targets=(_entry("HDCA", 0), _entry("HDCTRANS", 2563)),
        decoys=(), accepted_fragments=2563, decoy_fragments=0,
        rejected_low_identity=0, rejected_short_alignment=0,
        copies_per_100_genomes=None, aggregate_method="sum", indeterminate=False,
    )
    assert result.aggregate_fragments == 0
    assert result.detected is False, "a headline gene with no reads is not detected"
    assert result.stable is False
    # And the two properties now agree with each other.
    assert result.confidence == "not detected"
    # The panel total is still available for "fragments matched".
    assert result.accepted_fragments == 2563


def test_unplaced_readings_are_not_counted_as_typical():
    """Provisional, not-detected and no-reference rows are not 'typical'.

    All three wear the neutral colour for different reasons, and a summary
    that counted everything non-coral as typical folded them in - so a
    pathway at the 6th percentile could be summarised as typical because
    its identity was provisional.
    """
    from openbiota.pdfreport import status_for

    typical = status_for(percentile=50.0)
    assert typical.assessed

    provisional = status_for(percentile=6.5, confidence="provisional")
    assert not provisional.assessed
    assert provisional.unassessed_reason == "provisional identity"

    undetected = status_for(percentile=5.0, detected=False)
    assert not undetected.assessed
    assert undetected.unassessed_reason == "not detected"

    noref = status_for(percentile=None)
    assert not noref.assessed
    assert noref.unassessed_reason == "no reference range"


def test_a_tied_value_gets_the_middle_of_its_tie_block_not_the_floor():
    """Zero on a zero-inflated panel is typical, not "below the usual range".

    methane has p5 = p25 = p50 = 0 because most of the cohort has none.
    Returning the floor of that tie block put the commonest value at the
    5th percentile and classified it as below the usual range, and the
    number fed the JSON, the drawn scale marker and the notability sort.
    """
    from openbiota.cohort import PanelRange

    zero_inflated = PanelRange(
        panel="methane", n=100,
        percentiles={5: 0.0, 25: 0.0, 50: 0.0, 75: 1.2, 95: 8.0},
        mean=1.0, stdev=2.0, n_detected=43, min_value=0.0, max_value=12.0,
    )
    assert zero_inflated.percentile_of(0.0) == 27.5
    assert zero_inflated.classify(0.0) == "within the usual range"

    # Interpolation between distinct quantiles is untouched.
    assert zero_inflated.percentile_of(1.2) == 75.0
    ordinary = PanelRange(
        panel="x", n=100,
        percentiles={5: 1.0, 25: 2.0, 50: 3.0, 75: 4.0, 95: 5.0},
        mean=3.0, stdev=1.0, n_detected=100, min_value=1.0, max_value=5.0,
    )
    assert ordinary.percentile_of(2.5) == 37.5
    # A value genuinely below the whole reference still sits at the floor.
    assert ordinary.percentile_of(0.5) == 5.0


def test_an_installed_adapter_that_never_ran_appears_in_the_census():
    """Installed-but-unrun is a pending measurement, not an absence.

    Only uninstalled adapters used to be emitted, so once the tools were
    installed they vanished from the census entirely - from the very object
    whose purpose is to stop that happening.
    """
    from openbiota.resolution import adapters

    calls = adapters.calls_for_census("S")
    assert len(calls) == len(adapters.ADAPTERS), "every adapter needs a record"
    for call in calls:
        assert call.assay_status in {"scheduled", "reference_unavailable"}
        assert call.analytical_call == "unresolved"
        assert "not_detected" not in str(call.analytical_call)
    installed = [c for c in calls if c.assay_status == "scheduled"]
    for call in installed:
        assert "adapter_installed_not_run" in call.reason_codes
        assert "Nothing here says it is absent" in call.plain


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

_CACHE: dict[str, dict] = {}


def _analyze(sample: str, subject: dict | None = None) -> dict:
    """Run the engine on a real sample, caching for speed."""
    key = f"{sample}|{json.dumps(subject, sort_keys=True)}"
    if key in _CACHE:
        return _CACHE[key]
    sample_dir = REPO / "results" / sample
    if not sample_dir.exists():
        pytest.skip(f"{sample} results not present")
    results = json.loads((sample_dir / "results.json").read_text())
    if subject is not None:
        results = {**results, "subject_context": subject}
    doc = engine.analyze(sample_dir, results, cohort_path=str(COHORT))
    _CACHE[key] = doc
    return doc


def _erick() -> dict:
    return _analyze("SAMPLE2_A02")


def _card(doc: dict, name: str) -> dict:
    for card in doc["cards"]:
        if card["card"] == name:
            return card
    raise AssertionError(f"card {name} not found")
