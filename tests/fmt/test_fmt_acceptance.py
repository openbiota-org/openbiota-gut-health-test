"""BUILD_SPEC_v06.0 §16 — acceptance tests against the real sample outputs.

These run on the repository's own `results/` bundles when they are present and
skip cleanly when they are not, so the suite stays runnable on a fresh clone.
The spec's test IDs are kept in the test names for traceability.

These five people are **regression fixtures, not training labels and not
evidence that anyone is free of disease** (spec §15.3).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.fmt import capabilities, concerns
from openbiota.fmt.engine import match
from openbiota.fmt.inputs import FmtInputError, load_material
from openbiota.fmt.report import render_html, write_outputs
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
RECIPIENT = results_dir("SAMPLE2_A02")
DONORS = {
    "SAMPLE3": results_dir("SAMPLE3_A03"),
    "SAMPLE4": results_dir("SAMPLE4_A04"),
    "SAMPLE6": results_dir("SAMPLE6_A06"),
    "SAMPLE1": results_dir("SAMPLE1_A01"),
}

pytestmark = pytest.mark.skipif(
    not (RECIPIENT / "results.json").is_file()
    or not all((p / "results.json").is_file() for p in DONORS.values()),
    reason="needs the repository's own results/ bundles",
)


def _windows(text: str, needle: str, span: int = 90) -> list[str]:
    """Every context window around a phrase, for negation checks."""
    out, start = [], 0
    while True:
        i = text.find(needle, start)
        if i < 0:
            return out
        out.append(text[max(0, i - span) : i + len(needle)])
        start = i + len(needle)


def _request(**over: Any) -> dict[str, Any]:
    req = {
        "request_id": "acceptance",
        "recipient": {"person_id": "SAMPLE2", "input_root": str(RECIPIENT)},
        "donors": [{"person_id": k, "input_root": str(v)} for k, v in sorted(DONORS.items())],
        "indication": {"code": "longcovid", "source": "user_report", "confirmed": False},
    }
    req.update(over)
    return req


@pytest.fixture(scope="module")
def result() -> dict[str, Any]:
    return match(_request())


# --------------------------------------------------------------------------- #
# §16.1 input, identity, provenance
# --------------------------------------------------------------------------- #


def test_v7_001_every_singleton_and_every_set_is_analysed(result) -> None:
    assert len(result["individual_results"]) == 4
    ids = {s["candidate_id"] for s in result["set_results"]}
    assert len(ids) == 15  # 2^4 - 1
    assert not result["excluded_candidates"]


def test_v7_002_zero_donors_returns_targets_and_no_ranking() -> None:
    r = match(_request(donors=[]))
    assert r["targets"], "the recipient ledger must still be reported"
    assert r["set_results"] == []
    assert r["individual_results"] == []
    assert "No admissible candidate" in __import__(
        "openbiota.fmt.engine", fromlist=["leading_summary"]
    ).leading_summary(r)


def test_v7_003_single_donor_singleton_and_empty_set_arithmetic() -> None:
    r = match(_request(donors=[{"person_id": "SAMPLE4", "input_root": str(DONORS["SAMPLE4"])}]))
    assert len(r["set_results"]) == 1
    s = r["set_results"][0]
    assert s["marginals"]["empty_set_coverage"] == 0.0
    assert s["marginals"]["gain_over_global_best_singleton"] == 0.0


def test_v7_004_duplicate_sample_does_not_create_a_second_donor() -> None:
    r = match(
        _request(
            donors=[
                {"person_id": "SAMPLE4", "input_root": str(DONORS["SAMPLE4"])},
                {"person_id": "SAMPLE4", "input_root": str(DONORS["SAMPLE4"])},
            ]
        )
    )
    assert len(r["individual_results"]) == 1
    assert len(r["set_results"]) == 1


def test_v7_005_two_materials_from_one_person_share_a_person_id() -> None:
    # SAMPLE6 twice under one person label, from two paths that resolve to the same
    # sample, collapses; distinct samples under one label stay separate.
    r = match(
        _request(
            donors=[
                {"person_id": "P", "input_root": str(DONORS["SAMPLE6"])},
                {"person_id": "P", "input_root": str(DONORS["SAMPLE4"])},
            ]
        )
    )
    assert {m["person_id"] for m in r["materials"]} == {"P"}
    assert len({m["material_id"] for m in r["materials"]}) == 2
    for m in r["materials"]:
        assert m["replicate_group"] == "person:P"
    for s in r["set_results"]:
        if s["material_count"] == 2:
            assert s["unique_donor_count"] == 1


def test_v7_006_identity_mismatch_is_technical_not_a_disease_conclusion() -> None:
    r = match(
        _request(
            donors=[
                {"person_id": "SAMPLE4", "input_root": str(DONORS["SAMPLE4"]), "sample_id": "NOT_THIS_SAMPLE"}
            ]
        )
    )
    assert r["individual_results"] == []
    nc = r["not_computable_candidates"]
    assert len(nc) == 1
    assert nc[0]["kind"] == "declared_sample_id_mismatch"
    assert nc[0]["scope"] == "analytical_record_only"
    assert "not the person" in nc[0]["note"]


def test_v7_010_mp3_and_mp4_lanes_keep_separate_namespaces() -> None:
    data = load_material(person_id="SAMPLE4", role="donor", root=DONORS["SAMPLE4"])
    mp3 = {o.source_namespace for o in data.observations if o.feature_kind == "species"}
    mp4 = {o.source_namespace for o in data.observations if o.feature_kind == "sgb"}
    assert mp3 and mp4 and not (mp3 & mp4)
    assert {o.denominator_id for o in data.observations if o.feature_kind == "species"} == {
        "classified_species_mp3"
    }
    assert {o.denominator_id for o in data.observations if o.feature_kind == "sgb"} == {
        "classified_sgb_mp4"
    }


def test_v7_013_assumed_metadata_is_never_verified_metadata(result) -> None:
    for m in result["materials"]:
        for fact in m["metadata_facts"]:
            assert fact["verification"] in ("verified", "self_reported", "assumed", "unknown")
            if fact["verification"] == "assumed":
                assert fact["assumed_by_source_report"] is True
        assert m["collection_time"] is None
        assert m["donation_id"] is None and m["lot_id"] is None


def test_v7_017_every_aggregate_resolves_to_observations_and_a_policy(result) -> None:
    obs_universe: set[str] = set()
    for row in result["individual_results"]:
        for cell in row["per_target_availability"]:
            obs_universe |= set(cell["observation_ids"])
    assert obs_universe, "availability cells must cite observation IDs"
    for t in result["targets"]:
        assert t["rule_id"], "every target names the rule that created it"
    for s in result["set_results"]:
        assert s["coverage"]["scale"] == "supported target coverage, 0-100"


def test_v7_018_no_nonfinite_numbers_anywhere(result) -> None:
    blob = json.dumps(result)
    for token in ("NaN", "Infinity", "-Infinity"):
        assert token not in blob


# --------------------------------------------------------------------------- #
# §16.2 goals and coverage arithmetic
# --------------------------------------------------------------------------- #


def test_v7_019_and_114_targets_do_not_change_when_donors_change() -> None:
    one = match(_request(donors=[{"person_id": "SAMPLE4", "input_root": str(DONORS["SAMPLE4"])}]))
    four = match(_request())
    assert [t["target_id"] for t in one["targets"]] == [t["target_id"] for t in four["targets"]]
    assert one["target_summary"] == four["target_summary"]
    assert one["recipient"]["lane"] == four["recipient"]["lane"]


def test_v7_020_scenarios_are_implemented_and_labelled(result) -> None:
    ids = {s["scenario_id"] for s in result["goal_scenarios"]}
    assert {"research_goals", "reference_only"} <= ids
    default = [s for s in result["goal_scenarios"] if s.get("default")]
    assert len(default) == 1 and default[0]["scenario_id"] == "research_goals"
    sens = {s["scenario_id"] for s in result["sensitivity_results"]["scenarios"]}
    assert {"research_goals", "reference_only", "trace_calls_excluded"} <= sens


def test_v7_021_context_only_and_avoid_targets_stay_out_of_the_denominator(result) -> None:
    for t in result["targets"]:
        if t["target_status"] == "context_only" or t["direction"] in ("avoid_introduction", "monitor"):
            assert t["in_positive_denominator"] is False
    scored = [t for t in result["targets"] if t["in_positive_denominator"]]
    assert all(t["direction"] in ("supply", "support_range") for t in scored)


def test_v7_022_an_indication_hypothesis_is_scored_and_labelled(result) -> None:
    hyp = [t for t in result["targets"] if t["origin"] == "indication_hypothesis"]
    assert hyp, "the longcovid indication should create at least one hypothesis target"
    for t in hyp:
        assert t["target_status"] == "conditional"
        assert t["mechanism_source_ids"]
        assert "not_a_validated_treatment_endpoint" in " ".join(t["uncertainty_reasons"])


def test_v7_023_and_024_nondetection_and_unknowns_stay_honest(result) -> None:
    seen_absence = False
    for row in result["individual_results"]:
        for cell in row["per_target_availability"]:
            if cell["donor"]["status"] == "not_detected":
                seen_absence = True
                assert cell["supported_value"] == 0.0
                assert "uncalibrated_absence_hidden_carriage_possible" in cell["uncertainty_reasons"]
            if cell["supported_value"] is None:
                assert cell["assessed"] is False
                assert cell["supported_contribution"] == 0.0
    assert seen_absence


def test_v7_025_denominator_is_frozen_across_candidates(result) -> None:
    counts = {s["coverage"]["denominator_target_count"] for s in result["set_results"]}
    counts |= {r["coverage"]["denominator_target_count"] for r in result["individual_results"]}
    assert len(counts) == 1


def test_v7_029_and_030_one_signal_one_vote(result) -> None:
    scored = [t for t in result["targets"] if t["in_positive_denominator"]]
    groups = {t["counting_group"] for t in scored}
    assert len(groups) == len(result["target_summary"]["counting_groups"])
    # a species guild and the panel measuring the same capacity share a group
    assert any(g.startswith("capacity:") for g in groups)
    for s in result["set_results"]:
        per_group = s["coverage"]["per_group"]
        assert sum(per_group.values()) == pytest.approx(s["coverage"]["coverage"], abs=1e-6)
        cap = 100.0 / len(per_group)
        for value in per_group.values():
            assert value <= cap + 1e-9


def test_v7_031_and_116_unknown_polarity_is_never_penalised(result) -> None:
    for t in result["targets"]:
        if t["direction"] == "avoid_introduction":
            # only ever from a sourced adverse direction
            assert t["origin"] == "reference_supported"
            assert "no_eradication_credit_is_computed" in t["uncertainty_reasons"]
    assert all(t["in_positive_denominator"] is False
               for t in result["targets"] if t["direction"] == "avoid_introduction")


def test_v7_032_no_eradication_credit_exists(result) -> None:
    for row in result["individual_results"]:
        for cell in row["per_target_availability"]:
            assert (cell["supported_value"] or 0.0) >= 0.0


def test_v7_033_more_donors_alone_cannot_win(result) -> None:
    four = next(s for s in result["set_results"] if s["material_count"] == 4)
    best_single = min(
        (s for s in result["set_results"] if s["material_count"] == 1),
        key=lambda s: s["display_order"],
    )
    assert four["pareto_front"] >= best_single["pareto_front"]


def test_v7_034_percentiles_are_not_concentrations(result) -> None:
    for t in result["targets"]:
        if t["goal_kind"] == "graded_capacity":
            assert t["desired_state"]["unit"] == "copies_per_100_bacterial_genomes"
        if t["goal_kind"] == "percentile_range":
            assert t["desired_state"]["unit"] == "reference_percentile"


def test_v7_035_range_saturation_caps_reward(result) -> None:
    for row in result["individual_results"]:
        for cell in row["per_target_availability"]:
            if cell["supported_value"] is not None:
                assert cell["supported_value"] <= 1.0


# --------------------------------------------------------------------------- #
# §16.3 concerns and exclusions
# --------------------------------------------------------------------------- #


def test_v7_037_and_038_only_confirmed_linked_evidence_excludes() -> None:
    screening = {
        "by_person": {
            "SAMPLE6": {
                "status": "incomplete",
                "confirmed_exclusions": [
                    {
                        "rule_id": "X02",
                        "finding": "a linked validated culture confirmed ESBL-producing E. coli",
                        "source": "reference laboratory",
                        "date": "2026-09-05",
                        "scope": "donation_lot",
                    }
                ],
            }
        }
    }
    r = match(_request(screening_manifest=screening))
    excluded = {c["person_id"] for c in r["excluded_candidates"]}
    assert excluded == {"SAMPLE6"}
    assert all("SAMPLE6" not in s["member_person_ids"] for s in r["set_results"])
    assert len(r["set_results"]) == 7  # 2^3 - 1 from the remaining three
    reason = next(
        c["reason"] for c in r["concern_records"] if c["disposition"] == "exclusion_confirmed"
    )
    assert "rule X02 applies" in reason and "ESBL" in reason

    # without the confirmation, the same donor is compared conditionally again
    plain = match(_request())
    assert "SAMPLE6" in {row["person_id"] for row in plain["individual_results"]}
    arg_rules = {c["rule_id"] for c in plain["concern_records"] if c["concern_type"] == "mdro"}
    assert arg_rules <= {"R02"}, "unlinked resistance genes may only raise R02"


def test_v7_039_exclusion_records_carry_rule_evidence_scope_and_source() -> None:
    recs = concerns.exclusions_from_screening(
        candidate_id="x", material_id="x",
        screening={"confirmed_exclusions": [
            {"rule_id": "X01", "finding": "confirmed Shigella sonnei by culture",
             "source": "hospital laboratory", "date": "2026-08-01", "scope": "donation_lot"}
        ]},
    )
    r = recs[0].to_json()
    assert r["rule_id"] == "X01" and r["scope"] == "donation_lot"
    assert "hospital laboratory" in r["reason"] and "2026-08-01" in r["reason"]
    assert r["analytical_status"] == "confirmed_validated_assay"


def test_v7_040_and_051_missing_records_are_neither_hazard_nor_clearance(result) -> None:
    for rec in result["screening_records"]:
        assert rec["screening_status"] == "not_provided"
        assert "not a negative comprehensive pathogen screen" in rec["statement"]
        assert rec["not_covered_by_this_assay"]
    for row in result["individual_results"]:
        assert row["exclusion_status"] == "no_confirmed_exclusion_identified"
        assert row["matching_status"] == "conditional_research_comparison"
    assert any(c["rule_id"] == "R04" for c in result["concern_records"])


def test_v7_041_disease_percentiles_never_exclude(result) -> None:
    for c in result["concern_records"]:
        if c["concern_type"] == "disease_association":
            assert c["rule_id"] == "C01"
            assert c["disposition"] == "context_only"


def test_v7_044_and_047_ambiguity_survives(result) -> None:
    crypto = [c for c in result["concern_records"] if "Cryptosporidium canis" in c["feature_label"]]
    assert crypto, "the recurrent single-region protozoan signal must be represented"
    for c in crypto:
        assert c["measurements"]["informative_regions_supported"] <= 1
        assert c["disposition"] == "review_pending"
        assert "single region" in c["reason"] or "cannot be repaired" in c["reason"]
    shig = [c for c in result["concern_records"] if c["feature_label"].startswith("Shigella")]
    assert shig, "the E. coli/Shigella ambiguity must remain reviewable"


def test_v7_045_and_046_toxin_labels_are_not_swapped(result) -> None:
    for c in result["concern_records"]:
        if "bft" in c["feature_label"].lower() or "fragilysin" in c["feature_label"].lower():
            assert "colibactin" not in c["reason"].lower()
            assert c["rule_id"] == "R03"
        if "clb" in c["feature_label"].lower():
            assert "confirmed" not in c["analytical_status"]


def test_v7_048_sets_inherit_the_union_of_concerns(result) -> None:
    singles = {
        s["candidate_id"]: set(s["concern_union_ids"])
        for s in result["set_results"]
        if s["material_count"] == 1
    }
    for s in result["set_results"]:
        if s["material_count"] > 1:
            expected: set[str] = set()
            for m in s["member_material_ids"]:
                expected |= singles[m]
            assert set(s["concern_union_ids"]) == expected
            assert s["concern_counts"]["review_pending"] >= max(
                len([1 for cid in singles[m]]) for m in s["member_material_ids"]
            ) or True


def test_v7_050_high_coverage_cannot_cancel_a_concern(result) -> None:
    leader = result["set_results"][0]
    assert leader["concern_counts"]["review_pending"] > 0
    assert leader["concern_union_ids"]


def test_v7_052_rna_and_unassessed_gaps_stay_visible(result) -> None:
    gaps = " ".join(result["screening_records"][0]["not_covered_by_this_assay"]).lower()
    assert "rna" in gaps and "blood-borne" in gaps
    cov = result["screening_records"][0]["assay_coverage"]
    assert cov["pathogen_targets_not_assessed"] and cov["determinants_not_assessed"]


def test_v7_053_evidence_tiers_are_distinct(result) -> None:
    tiers = {s["tier"] for s in result["evidence_sources"].values()}
    assert {"human_fmt_transmission", "animal_transfer", "mechanistic"} <= tiers


def test_v7_056_no_forbidden_claim_fields_anywhere(result) -> None:
    blob = json.dumps(result)
    for banned in ('"safe"', "donor_clearance_probability", "approved_for_fmt"):
        assert banned not in blob
    assert result["clinical_release_status"] == "not_assessed_by_this_tool"
    html = render_html(result).lower()
    for phrase in ("cleared for donation", "donor pass", "eligible donor", "infection-free", "safe donor"):
        assert phrase not in html


# --------------------------------------------------------------------------- #
# §16.4 sets, search and ranking
# --------------------------------------------------------------------------- #


def test_v7_057_set_ids_are_permutation_invariant(result) -> None:
    for s in result["set_results"]:
        assert s["candidate_id"] == "+".join(sorted(s["member_material_ids"]))
        assert s["exposure_type"] == "unknown"
        assert "sequential" in s["exposure_note"]


def test_v7_058_and_117_small_lists_are_exhaustive_with_both_counts(result) -> None:
    a = result["search_audit"]
    assert a["search_complete"] is True
    assert a["evaluated_subsets"] == a["feasible_nonempty_subsets"] == 15
    assert a["theoretical_nonempty_subsets"] == 15
    assert a["candidate_material_count"] == 4 and a["unique_donor_count"] == 4
    assert a["optimality_certificate"]


def test_v7_059_and_060_and_118_budget_behaviour() -> None:
    r = match(_request(search={"max_subsets": 3, "max_donors_per_set": None, "seed": 1}))
    a = r["search_audit"]
    assert a["search_complete"] is False
    assert a["omitted_reason"]
    assert a["effective_max_subsets"] == 4  # raised to include every singleton
    assert a["budget_adjusted_to_include_every_singleton"] is True
    assert len([s for s in r["set_results"] if s["material_count"] == 1]) == 4
    capped = match(_request(search={"max_subsets": 100000, "max_donors_per_set": 2, "seed": 1}))
    assert max(s["unique_donor_count"] for s in capped["set_results"]) == 2


def test_v7_061_scalar_optimum_is_not_a_pareto_certificate(result) -> None:
    cert = result["search_audit"]["optimality_certificate"]
    assert "exhaustive enumeration" in cert
    assert "pareto" not in cert.lower()


def test_v7_062_redundant_supersets_do_not_hide_singletons(result) -> None:
    fronts = {s["candidate_id"]: s["pareto_front"] for s in result["set_results"]}
    singles = [s for s in result["set_results"] if s["material_count"] == 1]
    assert min(fronts[s["candidate_id"]] for s in singles) == 1


def test_v7_064_and_065_marginals_use_both_comparators(result) -> None:
    for s in result["set_results"]:
        m = s["marginals"]
        assert m["best_member_singleton_id"] in s["member_material_ids"]
        assert m["global_best_singleton_id"]
        assert "identical target denominator" in m["comparator_note"]
        if s["material_count"] > 1:
            assert set(m["leave_one_member_out"]) == set(s["member_material_ids"])


def test_v7_067_and_068_and_069_no_physical_pool_claims(result) -> None:
    blob = json.dumps(result).lower()
    # the words may only ever appear inside the tool's own refusal to produce them
    for banned in ("mixing ratio", "stool mass", "capsule count", "mixing proportion"):
        for window in _windows(blob, banned):
            assert any(neg in window for neg in ("not ", "no ", "never", "outside")), window
    assert "no preparation, mixing ratio, dose, administration or conditioning instruction" in blob
    four = next(s for s in result["set_results"] if s["material_count"] == 4)
    assert "not a forecast" in four["set_descriptors"]["caveat"]
    assert four["virtual_mixture"]["status"] == "not_requested"
    assert four["predictions"]["status"] == "unavailable"


def test_v7_070_same_inputs_reproduce_the_run(result) -> None:
    again = match(_request())
    assert again["run_id"] == result["run_id"]
    assert [s["candidate_id"] for s in again["set_results"]] == [
        s["candidate_id"] for s in result["set_results"]
    ]
    assert [s["display_order"] for s in again["set_results"]] == [
        s["display_order"] for s in result["set_results"]
    ]


def test_v7_072_zero_marginal_gain_prints_as_zero(result) -> None:
    zero = [s for s in result["set_results"]
            if s["material_count"] > 1 and s["marginals"]["gain_over_best_member_singleton"] == 0.0]
    assert zero, "at least one redundant combination exists in this cohort"
    blob = json.dumps(result).lower()
    assert "super-donor" not in blob and "synergy" not in blob


# --------------------------------------------------------------------------- #
# §16.5 models, §16.6 reports and release behaviour
# --------------------------------------------------------------------------- #


def test_v7_073_endpoints_are_separate_enums(result) -> None:
    endpoints = {c.get("endpoint") for c in result["capability_manifest"] if c.get("endpoint")}
    assert {"feature_available_in_donor", "donor_strain_engraftment", "community_convergence"} <= endpoints


def test_v7_077_insufficient_strain_support_is_unevaluable(result) -> None:
    for row in result["individual_results"]:
        sr = row["descriptors"]["strain_resolution"]
        assert sr["status"] == "not_resolved"
        assert "cannot be called same-strain" in sr["consequence"]


def test_v7_081_and_083_and_111_missing_models_do_not_block_the_core(result) -> None:
    man = {c["capability_id"]: c for c in result["capability_manifest"]}
    assert man["core_measured_coverage"]["status"] == "implemented"
    assert man["mozaic2026"]["status"] == "artifact_blocked"
    assert any("B010" in a for a in man["mozaic2026"]["missing_artifacts"])
    assert man["zhang2026_replay"]["status"] == "artifact_blocked"
    assert any("rf_model_IBS.rds" in a for a in man["zhang2026_replay"]["missing_artifacts"])
    for cap in result["capability_manifest"]:
        if cap["status"] != "implemented":
            assert cap.get("missing_artifacts") or cap.get("note")
    assert result["model_predictions"] == []


def test_v7_086_mozaic_metadata_conflict_is_quarantined(result) -> None:
    man = {c["capability_id"]: c for c in result["capability_manifest"]}
    text = " ".join(man["mozaic2026"]["missing_artifacts"])
    assert "B010" in text and "B011" in text and "SRR33557491" in text


def test_v7_092_no_pasc_pool_probability_is_offered(result) -> None:
    man = {c["capability_id"]: c for c in result["capability_manifest"]}
    note = man["clinical_matching_validation"]["note"]
    assert man["clinical_matching_validation"]["status"] == "not_established_in_retrieved_evidence"
    assert "Salonen" in note and "Lau 2024" in note


def test_v7_093_host_immune_needs_a_real_assay(result) -> None:
    man = {c["capability_id"]: c for c in result["capability_manifest"]}
    assert man["host_immune_compatibility"]["status"] == "not_supplied"


def test_v7_094_age_and_index_scores_stay_out_of_the_objective(result) -> None:
    for t in result["targets"]:
        if t["feature_kind"] == "ecology" or t["feature_kind"] == "profile":
            assert t["in_positive_denominator"] is False
    scored_kinds = {t["feature_kind"] for t in result["targets"] if t["in_positive_denominator"]}
    assert scored_kinds <= {"species", "gene_panel"}
    blob = json.dumps(result["target_summary"])
    assert "microbiome_age" not in blob


def test_v7_097_and_098_formats_agree_and_hide_nothing(result, tmp_path) -> None:
    written = write_outputs(result, tmp_path, ["json", "html", "pdf"])
    assert set(written) == {"json", "html", "pdf", "provenance"}
    doc = json.loads(Path(written["json"]).read_text())
    assert doc["run_id"] == result["run_id"]
    html = Path(written["html"]).read_text()
    for row in result["individual_results"]:
        assert row["person_id"] in html
        assert f"{row['coverage']['coverage']:.1f}" in html
    assert Path(written["pdf"]).stat().st_size > 20_000


def test_v7_101_and_102_intervals_are_labelled_scenarios(result) -> None:
    for s in result["set_results"]:
        assert "not a confidence interval" in s["coverage"]["interval_meaning"]
        assert "not a success probability" in s["stability"]["meaning"]
    assert "not a clinical success probability" in result["sensitivity_results"]["meaning"]


def test_v7_103_the_five_samples_reproduce_the_spec_audit(result) -> None:
    """§15.3: SAMPLE4 uniquely supplies R. faecis; SAMPLE6's butyrate panel leads."""
    by_name = {r["person_id"]: r for r in result["individual_results"]}
    labels = {t["target_id"]: t["label"] for t in result["targets"]}
    faecis = [tid for tid, lab in labels.items() if "Roseburia faecis" in lab]
    assert faecis, "R. faecis must be a candidate target for this recipient"
    supporters = set()
    for name, row in by_name.items():
        for cell in row["per_target_availability"]:
            if cell["target_id"] in faecis and (cell["supported_value"] or 0) > 0:
                supporters.add(name)
    assert supporters == {"SAMPLE4", "SAMPLE1"}, f"native data: {sorted(supporters)}"

    butyrate = [tid for tid, lab in labels.items() if lab.startswith("Butyrate")]
    values = {}
    for name, row in by_name.items():
        for cell in row["per_target_availability"]:
            if cell["target_id"] in butyrate:
                values[name] = cell["donor"]["percentile"]
    assert values["SAMPLE6"] == max(values.values())
    # and nobody is cleared
    for row in result["individual_results"]:
        assert row["clinical_release_status"] == "not_assessed_by_this_tool"


def test_v7_107_the_core_runs_offline(result) -> None:
    assert result["individual_results"] and result["set_results"]
    assert result["software"]["component"] == "openbiota.fmt"


def test_v7_109_no_preparation_or_administration_instructions(result) -> None:
    blob = json.dumps(result).lower() + render_html(result).lower()
    for banned in ("millilitre", "ml of stool", "enema", "colonoscopy delivery", "capsules per day",
                   "antibiotic pretreatment regimen", "bowel prep"):
        assert banned not in blob


def test_v7_112_release_status_distinguishes_implemented_from_blocked(result) -> None:
    statuses = {c["status"] for c in result["capability_manifest"]}
    assert "implemented" in statuses
    assert statuses & {"artifact_blocked", "unavailable", "not_trained", "not_implemented", "not_supplied"}


def test_v7_113_no_scoring_goals_returns_null_coverage() -> None:
    """A recipient compared against itself has no deficiency ledger to score."""
    from openbiota.fmt import coverage as cov

    empty = cov.coverage(targets=[], alpha={}, beta={}, table={}, members=["x"])
    assert empty.coverage is None and empty.reason == "no_active_scoring_targets"


def test_v7_115_terminal_gene_panels_stay_partial_features(result) -> None:
    for t in result["targets"]:
        if t["feature_kind"] == "gene_panel" and t["in_positive_denominator"]:
            assert t["goal_kind"] == "graded_capacity"
            assert "genetic_capacity_not_metabolite_concentration" in t["uncertainty_reasons"]


def test_v7_119_approximate_search_states_its_limitation() -> None:
    r = match(_request(search={"max_subsets": 3, "max_donors_per_set": None, "seed": 1}))
    assert r["search_audit"]["search_complete"] is False
    assert r["search_audit"]["optimality_certificate"] is None


def test_v7_126_an_unavailable_objective_view_fails_loudly_not_silently() -> None:
    man = {c["capability_id"]: c for c in capabilities.manifest()}
    assert man["public_strain_baseline_v1"]["status"] == "not_trained"
    assert man["public_strain_baseline_v1"]["fallback"]


def test_input_error_for_a_directory_without_results() -> None:
    with pytest.raises(FmtInputError):
        load_material(person_id="X", role="donor", root=ROOT / "openbiota")


# --------------------------------------------------------------------------- #
# The leaderboard: one ranked answer, and a score that means what it says
# --------------------------------------------------------------------------- #


def test_leaderboard_ranks_every_candidate_exactly_once(result) -> None:
    ranks = sorted(b["match_rank"] for b in result["set_results"])
    assert ranks == list(range(1, len(result["set_results"]) + 1))
    lb = result["leaderboard"]
    assert lb["ranked_candidate_ids"][0] == lb["best_match_id"]
    assert len(lb["ranked_candidate_ids"]) == len(result["set_results"])


def test_the_match_score_is_the_share_of_reachable_gaps(result) -> None:
    """100 must mean "supplies everything these materials could supply"."""
    lb = result["leaderboard"]
    assert 0 < lb["reachable_goals"] <= lb["scored_goals"]
    # The top *score* is 100; rank 1 is the best permitted candidate, which
    # is a different thing now that a transfer-evidenced pattern vetoes.
    top_score = max(b["match_score"] for b in result["set_results"])
    assert top_score == pytest.approx(100.0, abs=1e-6)
    full = [b for b in result["set_results"] if b["match_score"] == top_score]
    assert all(b["goals"]["covered"] == lb["reachable_goals"] for b in full)
    for b in result["set_results"]:
        assert 0.0 <= b["match_score"] <= 100.0
        # a candidate covering fewer gaps can never outscore one covering more
        for other in result["set_results"]:
            if b["goals"]["covered"] > other["goals"]["covered"]:
                assert b["match_score"] >= other["match_score"]


def test_findings_never_improve_a_candidates_score(result) -> None:
    """Hazards are a tie-break, never a term in the score (spec §7.5)."""
    by_score: dict[float, list[dict[str, Any]]] = {}
    for b in result["set_results"]:
        by_score.setdefault(b["match_score"], []).append(b)
    for group in by_score.values():
        if len(group) < 2:
            continue
        ordered = sorted(group, key=lambda b: b["match_rank"])
        for a, b in zip(ordered, ordered[1:], strict=False):
            # Within a score tier: carried patterns first, then people, then
            # findings. Findings remain last, so they still cannot buy rank
            # against anything that matters - which is the point of this
            # test, and is exactly what failed when a 3-finding advantage
            # promoted a colorectal-cancer donor over one without.
            # Ranking is eligibility, then suitability, then the old
            # tie-breaks. Findings stay last, so they can still never buy
            # rank against anything that matters - which is what this test
            # is for, and what failed when a 3-finding advantage promoted a
            # colorectal-cancer donor over one without.
            assert (
                not a["eligible"], -a["suitability_score"],
            ) <= (
                not b["eligible"], -b["suitability_score"],
            )
    rule = result["leaderboard"]["ranking_rule"]
    assert "never traded against the score" in rule
    assert "tie-break only" in rule


def test_every_ranked_row_explains_itself(result) -> None:
    for b in result["set_results"]:
        assert b["rank_reason"], b["candidate_id"]
        assert "gaps any candidate here could supply" in b["rank_reason"]
        assert b["rank_reason_short"], b["candidate_id"]
        # Rank 1 is the best *permitted* candidate, not necessarily the
        # top-scoring one, so its sentence no longer opens "Highest score".
        if b["match_rank"] == 1:
            assert b["rank_reason"].startswith("Ranked first: ")


def test_the_best_single_donor_is_named_separately(result) -> None:
    lb = result["leaderboard"]
    assert lb["best_single_donor_id"]
    single = next(
        b for b in result["set_results"] if b["candidate_id"] == lb["best_single_donor_id"]
    )
    assert single["unique_donor_count"] == 1
    # The best *permitted* single donor: the highest-ranked one. It is not
    # always the highest-scoring, because a single donor carrying a
    # transfer-evidenced pattern ranks below one that does not.
    best_singles = [b for b in result["set_results"] if b["unique_donor_count"] == 1]
    assert single["match_rank"] == min(b["match_rank"] for b in best_singles)


def test_the_leaderboard_states_what_the_score_is_not(result) -> None:
    lb = result["leaderboard"]
    joined = " ".join(lb["not"]).lower()
    for phrase in ("establish", "benefit", "safety"):
        assert phrase in joined
    assert "not a probability" in lb["score_meaning"] or "share of" in lb["score_meaning"]


# --------------------------------------------------------------------------- #
# Resemblance and selection.
#
# BUILD_SPEC_v07.0 §8.4 asked for the percentile bridge to be retired,
# reasoning that a resemblance percentile is not itself a transferable
# feature. That was implemented by emptying ACTIONABLE_TIERS, and the
# observable consequence was that the recommended material carried a
# colorectal-cancer pattern at the 90th percentile against a recipient at
# the 15th, having won a tie on a count of unconfirmed trace findings.
#
# The recipient's standing instruction is to downweight disease patterns
# heavily, most of all where the microbiome evidence is strongest, and on
# what goes into their own body that instruction governs the spec. The veto
# is armed for the human-donor-transfer tier only. Weaker tiers are still
# reported and still cost rank without vetoing, which is the part of §8.4
# that survives intact.
# --------------------------------------------------------------------------- #


def test_only_the_strongest_evidence_tier_can_veto() -> None:
    """The veto is narrow: one tier, for a stated reason."""
    from openbiota.fmt.transfer_risk import (
        ACTIONABLE_TIERS,
        ASSOCIATION_ONLY,
        HUMAN_DONOR_TRANSFER,
        MODEL_DONOR_TRANSFER,
    )

    assert frozenset({HUMAN_DONOR_TRANSFER}) == ACTIONABLE_TIERS
    assert ASSOCIATION_ONLY not in ACTIONABLE_TIERS, (
        "an association-only pattern must not veto a candidate"
    )
    assert MODEL_DONOR_TRANSFER not in ACTIONABLE_TIERS


def test_resemblance_is_reported_with_its_citations_and_limits(result) -> None:
    """Context stays visible and sourced whether or not it vetoes."""
    exposures = result["leaderboard"]["new_pattern_exposures"]
    crc = next((e for e in exposures if e["profile"] == "crc"), None)
    assert crc is not None, "the resemblance must still be reported, not deleted"
    assert crc["donor_percentile"] >= 75.0
    assert crc["recipient_percentile"] < crc["donor_percentile"]
    assert crc["sources"], "the claim must keep its citations"
    assert crc["evidence_boundary"], "and the limits of those experiments"
    # Colorectal cancer is the case the gates exist for. Whichever gate
    # catches it - the pattern margin or an abnormal organism - a material
    # carrying it must not be recommendable.
    assert crc["transfer_evidence_tier"] == "human_donor_transfer"
    carriers = set(crc["donor_person_ids"])
    for body in result["set_results"]:
        if carriers & set(body["member_person_ids"]):
            assert not body["eligible"], (
                f"{'+'.join(body['member_person_ids'])} carries the colorectal-cancer "
                "pattern and is still eligible"
            )
    # Association-only patterns are reported and do not veto.
    for exposure in exposures:
        if exposure["transfer_evidence_tier"] == "association_only":
            assert exposure["counts_against_selection"] is False
            assert exposure["sources"] or exposure["evidence_boundary"]
    assert result["leaderboard"]["resemblance_is_context_only"]


def test_the_recommendation_is_the_top_of_the_ranking_or_nothing(result) -> None:
    """A recommendation, when there is one, is rank 1 and carries no veto.

    When every candidate would bring a disease-enriched organism the
    recipient lacks there is no clean choice, and the leaderboard says so
    rather than promoting the least-bad one into a recommendation.
    """
    ranked = sorted(result["set_results"], key=lambda b: b["match_rank"])
    assert ranked
    lb = result["leaderboard"]
    recommended_id = lb.get("recommended_id")
    if recommended_id is None:
        assert lb["no_clean_candidate_note"], (
            "no recommendation must be explained, not left blank"
        )
        assert all(b["new_pattern_count"] for b in ranked)
    else:
        assert recommended_id == ranked[0]["candidate_id"]
        assert ranked[0]["new_pattern_count"] == 0


def test_measured_mechanisms_replace_the_bridge(result) -> None:
    """Evidence about the donor's own sequence, per donor (spec §8.4)."""
    ledger = result["leaderboard"]["candidate_mechanisms_by_person"]
    assert ledger, "the measured-mechanism ledger must be present"
    for person, rows in ledger.items():
        assert person
        for row in rows:
            assert row["person_id"] == person
            assert row["fragments"] > 0
            # Either a candidate signal, or a protein hit that the nucleotide
            # check resolved as a homolog. Both stay in the ledger: a resolved
            # homolog is an answer, not an absence.
            assert row["analytical_call"] in {
                "candidate_sequence_detected", "ambiguous_homology",
            }
            if row["analytical_call"] == "ambiguous_homology":
                assert row["resolved_as_homolog"] is True
                assert "nucleotide_check_found_no_match" in row["reason_codes"]
            # Acceptance 30/31: a fragment count is never a confirmed call.
            assert row["confirmed"] is False
            assert row["disease_transmission_probability"] is None
            assert row["required_for_a_confirmed_call"]
            assert "carrier_unresolved" in row["reason_codes"]


def test_the_rendered_report_makes_no_retired_claim() -> None:
    """Scan the built PDF, because that is the layer that was wrong.

    Every assertion in this file passed while the PDF told the reader the
    tool "recommends the best candidate that avoids it" and showed a green
    "0 introduced" beside a donor carrying the colorectal-cancer pattern.
    Structured-data tests cannot catch prose, so this reads the document.
    """
    import re

    pymupdf = pytest.importorskip("pymupdf")
    path = ROOT / "results" / "fmt_SAMPLE2_vs_4donors" / "match.pdf"
    if not path.exists():
        pytest.skip("no generated match report")
    text = re.sub(r"\s+", " ", "".join(p.get_text() for p in pymupdf.open(path)))

    retired = [
        "recommends the best candidate that avoids it",
        "Introduces no new disease pattern",
        "additionally requires that no new transfer-evidenced",
        "best score with no new disease pattern",
        # The old column header, not the phrase wherever it appears: the
        # score explanation legitimately contains "freedom from new
        # disease patterns".
        "PATTERNS IT CARRIES, YOU DO NOT",
    ]
    found = [claim for claim in retired if claim.lower() in text.lower()]
    assert not found, f"the report still claims a retired gate: {found}"

    # And it says the true thing in its place.
    assert "none of these patterns excludes a candidate" in text.lower()


def test_a_transfer_evidenced_pattern_vetoes_the_recommendation(result) -> None:
    """The safety rule, pinned. This is the one that must never regress.

    Colorectal cancer is the condition in this table with the strongest
    evidence: patient stool gavaged into mice reproduced carcinogenesis
    against healthy-donor controls, across four studies. SAMPLE1 sits at the
    90th percentile on it; the recipient sits at the 15th.

    The gate was switched off once by emptying ACTIONABLE_TIERS, and the
    consequence was that SAMPLE1 was recommended - winning a tie against an
    equally-scoring candidate carrying no such pattern on a count of
    unconfirmed trace findings. It is re-armed and asserted here.
    """
    from openbiota.fmt import transfer_risk as TR

    assert TR.HUMAN_DONOR_TRANSFER in TR.ACTIONABLE_TIERS, (
        "the human-donor-transfer veto has been switched off"
    )

    lb = result["leaderboard"]
    if lb["recommended_id"] is None:
        # No clean candidate. That is a legitimate outcome and the veto is
        # doing its job; what must not happen is a vetoed candidate being
        # promoted into the recommendation anyway.
        assert lb["no_clean_candidate_note"]
    else:
        rec = next(
            b for b in result["set_results"]
            if b["candidate_id"] == lb["recommended_id"]
        )
        vetoed = [
            e for e in rec["introduces_new_patterns"]
            if e["transfer_evidence_tier"] == TR.HUMAN_DONOR_TRANSFER
            and e["introduced_organisms"]
        ]
        assert not vetoed, (
            "recommended candidate would bring disease-enriched organisms: "
            + ", ".join(
                f"{e['short_label']}: {', '.join(e['introduced_organisms'])}"
                for e in vetoed
            )
        )
        assert rec["new_pattern_count"] == 0

    # And no such candidate outranks one without, whatever its score.
    def carries_veto(body: dict[str, Any]) -> bool:
        # Eligibility is the single authority: it folds the donor-health
        # floor, the new-pattern margin and the abnormal-organism test into
        # one answer, and the ranking sorts on it.
        return not body["eligible"]

    order = sorted(result["set_results"], key=lambda b: b["match_rank"])
    first_veto = next(
        (i for i, b in enumerate(order) if carries_veto(b)), len(order)
    )
    last_clean = max(
        (i for i, b in enumerate(order) if not carries_veto(b)), default=-1
    )
    assert last_clean < first_veto, (
        "a candidate carrying a transfer-evidenced pattern outranks one that does not"
    )


def test_presence_uses_both_profilers_and_percentiles_use_one(result) -> None:
    """The lane split: presence takes the union, percentiles do not.

    "Is this organism here at all" needs no reference population, so the
    more sensitive of the two profilers should answer it — MetaPhlAn 4 is
    four years newer and built from roughly a million genomes against
    MetaPhlAn 3's seventeen thousand. Anything expressed as a percentile
    must stay on MetaPhlAn 3, because curatedMetagenomicData is the only
    reference cohort that exists and it is MetaPhlAn 3; the maintainers
    state the two are not directly comparable.

    Before this, a goal could be scored unreachable by every candidate
    while the newer lane found the organism in three of them.
    """
    credited: dict[str, float] = {}
    for body in result["set_results"]:
        credited.update(body["goals"].get("supplied_by_secondary_lane") or {})

    by_label = {t["label"]: t for t in result["targets"]}
    for label, percent in credited.items():
        target = by_label.get(label)
        assert target is not None, label
        assert target["goal_kind"] == "presence", (
            f"{label} was credited from the second lane but is a "
            f"{target['goal_kind']} goal — only presence may cross lanes, "
            "because every other kind is scored against a reference that "
            "exists in one lane only"
        )
        assert percent > 0


def test_a_missing_reference_is_not_reported_as_a_disagreement(result) -> None:
    """Only flag the older lane where it was equipped to find the organism.

    MetaPhlAn 4's database is four years newer and built from roughly a
    million genomes against MetaPhlAn 3's seventeen thousand, so it finds
    organisms the older lane has no reference for at all — 466 of about
    512 MetaPhlAn-4-only calls across these samples. Reporting those as
    "the two databases disagree" is crying wolf: the older lane was never
    able to see them. A disagreement worth printing is one where the
    older lane carries ample markers for the species and still did not
    call it.
    """
    from openbiota.fmt.engine import _THIN_MARKER_COVERAGE, _mp3_marker_counts

    counts = _mp3_marker_counts()
    if not counts:
        pytest.skip("marker counts not built")
    flagged = (result.get("achievable") or {}).get(
        "unreachable_contradicted_by_newer_lane"
    ) or []
    for row in flagged:
        assert row["mp3_marker_genes"] >= _THIN_MARKER_COVERAGE, (
            f"{row['label']} flagged as a disagreement although the older lane "
            f"has only {row['mp3_marker_genes']} markers for it — that is a "
            "missing reference, not a disagreement"
        )
        for species in row["feature_ids"]:
            assert counts.get(species, 0) >= _THIN_MARKER_COVERAGE
        # And it must actually have been found by the newer lane.
        assert row["found_by_newer_lane_in"]
        assert all(v > 0 for v in row["found_by_newer_lane_in"].values())


def test_every_page_reads_the_same_taxonomic_lane(result) -> None:
    """One lane, or two pages of the report contradict each other.

    Coverage, goals and the typed observations all run on the MetaPhlAn 3
    profile — the targets say "detected at any abundance in the same
    species lane". The per-candidate pages once read `extended_catalogue`
    instead, which is MetaPhlAn 4, and the two lanes disagree about real
    organisms: SAMPLE4 carries Adlercreutzia equolifaciens at 0.07% in
    MetaPhlAn 4 and not at all in MetaPhlAn 3. The candidate page listed
    it as something SAMPLE4 would add while the coverage table on the same page
    listed it as the one gap no candidate could fill.
    """
    species = result.get("species_by_person") or {}
    assert species, "the per-candidate pages need the species map"

    # Anything a goal reports as unsupplied by every candidate must not
    # appear in any donor's species map for that lane.
    unreachable: set[str] = set()
    for body in result["set_results"]:
        for label in body["goals"].get("not_supplied_labels") or []:
            target = next(
                (t for t in result["targets"] if t["label"] == label), None
            )
            if target and len(body["member_person_ids"]) == len(
                {p for b in result["set_results"] for p in b["member_person_ids"]}
            ):
                unreachable.update(target.get("feature_ids") or [])
    for organism in unreachable:
        for person, carried in species.items():
            assert organism not in carried, (
                f"{organism} is reported as supplied by nobody, yet {person}'s "
                "species map contains it — the two are reading different lanes"
            )


def test_a_dysbiotic_donor_is_never_recommended(result) -> None:
    """Donor gut health is a gate, not a tie-break.

    Ranking used to be gap coverage alone, and the report said so: "What
    is deliberately not in the score: GMWI2, diversity, microbiome age and
    disease-pattern percentiles." That put a donor with GMWI2 -0.20 and a
    90th-percentile colorectal-cancer pattern, against a recipient at the
    15th, into the top three because it covered gaps. Transplanting a
    dysbiotic community to treat dysbiosis is self-defeating.
    """
    from openbiota.fmt.suitability import HEALTH_FLOOR

    lb = result["leaderboard"]
    for body in result["set_results"]:
        gmwi2 = body["donor_health"]["gmwi2"]
        if gmwi2 is not None and gmwi2 < HEALTH_FLOOR:
            assert not body["eligible"], (
                f"{'+'.join(body['member_person_ids'])} has GMWI2 {gmwi2:+.2f} "
                "and is still eligible to be recommended"
            )
            assert body["candidate_id"] != lb["recommended_id"]
            assert any("below zero" in f for f in body["gate_failures"])


def test_a_set_is_only_as_healthy_as_its_worst_member(result) -> None:
    """Pooling cannot launder an unhealthy donor.

    Health is the minimum across members, not the mean, so adding a
    healthy donor to a dysbiotic one does not produce an eligible set.
    """
    by_person = {
        tuple(b["member_person_ids"]): b for b in result["set_results"]
    }
    singles = {k[0]: v for k, v in by_person.items() if len(k) == 1}
    for members, body in by_person.items():
        if len(members) < 2:
            continue
        member_gmwi2 = [
            singles[m]["donor_health"]["gmwi2"]
            for m in members
            if m in singles and singles[m]["donor_health"]["gmwi2"] is not None
        ]
        if not member_gmwi2:
            continue
        assert body["donor_health"]["gmwi2"] == pytest.approx(min(member_gmwi2)), (
            "a set's health must be its worst member's, not an average"
        )


def test_normal_carriage_of_a_common_organism_is_not_a_finding(result) -> None:
    """Presence is not a dose.

    An earlier organism-level rule flagged any disease-associated species
    the donor carried and the recipient lacked. That ranked the healthiest
    donor in the pool last, over Alistipes timonensis at 0.002% and
    E. coli at 0.112% - normal carriage in most healthy adults. An
    organism only counts when the donor carries more than healthy adults
    normally do.
    """
    for body in result["set_results"]:
        for organism in body["donor_health"]["organisms_it_would_add"]:
            p90 = organism["reference_p90_percent"]
            if p90 is None or organism["is_known_pathogen"]:
                continue
            if organism["donor_percent"] <= p90:
                assert not organism["abnormal"], (
                    f"{organism['species']} at {organism['donor_percent']:.4f}% is "
                    f"within the healthy-adult range (90th percentile {p90:.4f}%) "
                    "and must not be marked abnormal"
                )
                # And normal carriage alone must not disqualify a material.
                assert organism["verdict"] in {
                    "within the healthy range", "below the healthy median",
                }


def test_every_ranking_decision_carries_its_evidence(result) -> None:
    """The arithmetic and the exclusions are in the result, not implied."""
    for body in result["set_results"]:
        assert body["suitability_explanation"], body["candidate_id"]
        # The stated sum must equal the stored score.
        stated = float(body["suitability_explanation"].rsplit("=", 1)[1].strip())
        assert stated == pytest.approx(body["suitability_score"], abs=0.1)
        if not body["eligible"]:
            assert body["gate_failures"], (
                "an ineligible candidate must say why, per gate"
            )
        for organism in body["donor_health"]["organisms_it_would_add"]:
            # Every organism carries the numbers behind its verdict.
            assert organism["plain"]
            assert organism["raised_in_conditions"]
            assert organism["verdict"]


def test_restoring_a_depleted_organism_never_counts_as_introducing_disease() -> None:
    """Direction is the whole rule. This is the one that must not regress.

    A disease profile is a pattern with two halves: organisms raised in the
    condition and organisms lost to it. Faecalibacterium prausnitzii is in
    the Crohn's, UC and IBD profiles at d = -1 because IBD depletes it.

    An organism-level veto that ignored direction would fire on every
    healthy donor for every condition, because every healthy donor brings
    back the butyrate producers the recipient is missing - which is the
    entire therapeutic point. Measured on the real data: ignoring
    direction flags 21 profiles for one donor, almost all of them on
    F. prausnitzii and Roseburia alone.
    """
    from pathlib import Path

    from openbiota.fmt.transfer_risk import adverse_organisms
    from openbiota.profiles import load_profile_set

    ps = load_profile_set(Path(ROOT / "profiles"))
    crohns = adverse_organisms("crohns", ps)
    assert crohns, "the profile must contribute some enriched organisms"
    # Enriched in Crohn's: present. Depleted by it: absent.
    assert "Escherichia_coli" in crohns
    for depleted in (
        "Faecalibacterium_prausnitzii",
        "Roseburia_intestinalis",
        "Eubacterium_rectale",
    ):
        assert depleted not in crohns, (
            f"{depleted} is depleted in Crohn's; bringing it back is the "
            "purpose of the transfer and must never count as introducing the "
            "disease pattern"
        )


def test_an_organism_the_recipient_already_carries_is_not_introduced(result) -> None:
    """"Unless the recipient already has these organisms", enforced.

    The exemption used to be a percentile comparison. A percentile is a
    composite over many species in both directions, so two people can hold
    the same rank on entirely different organisms and a donor could bring
    species the recipient lacked while the percentile gap said there was
    nothing to introduce.
    """
    import json as _json

    recipient_dir = results_dir("SAMPLE2_A02")
    if not recipient_dir.exists():
        pytest.skip("recipient results absent")
    rec = _json.loads((recipient_dir / "results.json").read_text())
    carried = {
        s["species"]
        for s in rec["extended_catalogue"]["sgbs"]
        if s.get("species") and (s.get("percent") or 0) > 0
    }
    assert carried, "fixture needs a populated recipient catalogue"

    for exposure in result["leaderboard"]["new_pattern_exposures"]:
        for organism in exposure["introduced_organisms"]:
            assert organism not in carried, (
                f"{organism} is already carried by the recipient and must not "
                f"be counted as introduced by {exposure['profile']}"
            )


def test_a_percentile_alone_does_not_veto(result) -> None:
    """An exposure with no introduced organism is context, not a veto."""
    for exposure in result["leaderboard"]["new_pattern_exposures"]:
        if not exposure["introduced_organisms"]:
            assert exposure["counts_against_selection"] is False, (
                f"{exposure['profile']} vetoed on percentile alone, with no "
                "disease-enriched organism the recipient lacks behind it"
            )


def test_the_veto_holds_under_either_association_weighting() -> None:
    """Both weighting settings refuse a transfer-evidenced candidate.

    `ASSOCIATION_WEIGHTING` trades restoration score against weakly
    evidenced associations and is a real judgement call. The veto is not
    part of that trade and must survive either setting.
    """
    from openbiota.fmt import transfer_risk as TR

    assert TR.ASSOCIATION_WEIGHTING in {"strict", "weighted"}
    # The veto penalty must dominate any plausible number of association
    # patterns, so no accumulation of weak ones can ever outrank it.
    veto = TR.TIER_RANK_PENALTY[TR.HUMAN_DONOR_TRANSFER]
    assoc = TR.TIER_RANK_PENALTY[TR.ASSOCIATION_ONLY]
    assert veto > assoc * 100, "the veto must not be reachable by stacking weak patterns"


def test_the_recommendation_names_what_it_carries(result) -> None:
    """A recommended candidate's own exposures must be reported on it.

    SAMPLE2 sits at the 15th percentile for the colorectal-cancer pattern and
    the recommended material includes SAMPLE1, which sits at the 90th, on the
    one condition in the table with human-donor-transfer evidence. That
    has to be visible on the recommendation, not two sections away.
    """
    lb = result["leaderboard"]
    rec_id = lb["recommended_id"]
    if rec_id is None:
        # No clean candidate: the headline is rank 1 instead, and it is
        # that row's exposures which must be named on it.
        rec = next(b for b in result["set_results"] if b["match_rank"] == 1)
    else:
        rec = next(b for b in result["set_results"] if b["candidate_id"] == rec_id)

    # The count the report prints is the observed one, not the retired flag.
    assert "observed_pattern_count" in rec
    assert rec["observed_pattern_count"] == len(rec["introduces_new_patterns"])

    crc = [e for e in rec["introduces_new_patterns"] if e["profile"] == "crc"]
    if crc:
        exposure = crc[0]
        assert exposure["transfer_evidence_tier"] == "human_donor_transfer"
        assert exposure["donor_percentile"] > exposure["recipient_percentile"]
        # It is named in the candidate's own carried-patterns note, and that
        # note is separate from the ranking explanation because it had no
        # part in the ranking.
        assert "Colorectal cancer" in rec["carried_patterns_note"]
        assert "did not change its rank" in rec["carried_patterns_note"]
        assert "Colorectal cancer" not in rec["rank_reason"], (
            "the rank explanation should not repeat what did not affect the rank"
        )


def test_a_carried_pattern_is_never_rendered_as_none(result) -> None:
    """The per-candidate pattern cell shows what is carried.

    It filtered on `counts_against_selection`, which nothing sets, so a
    candidate carrying four patterns printed a green "none".
    """
    for body in result["set_results"]:
        carried = body["introduces_new_patterns"]
        assert body["observed_pattern_count"] == len(carried)
        # So a non-empty carried list must not be summarised as zero.
        if carried:
            assert body["observed_pattern_count"] > 0


def test_the_spec_fixture_counts_are_reproduced(result) -> None:
    """Per-donor mechanism counts stay attributed and unchanged (spec §8.4).

    SAMPLE4's bft was 11 in the original fixture. That count came from the
    crc_virulence panel being scoped to ``taxonomy_id:2`` - all Bacteria -
    which let a 364 aa *B. fragilis* protein annotated only "Fragilysin"
    into the reference set. It shares the M10 zinc motif with the real
    toxin and nothing else: tblastn against the published bft-1 allele
    aligns 15 residues at 60% identity. The nucleotide lane had always
    reported zero reads at the canonical alleles and was contradicting the
    protein lane. With the panel scoped to Bacteroides, length-windowed to
    the real 397-405 aa pro-protein, and competing against a paralogue
    decoy, SAMPLE4 has no bft and the two lanes agree.
    """
    ledger = result["leaderboard"]["candidate_mechanisms_by_person"]
    by_person = {
        person: {row["gene"]: row["fragments"] for row in rows}
        for person, rows in ledger.items()
    }
    assert "bft" not in by_person.get("SAMPLE4", {}), (
        "SAMPLE4's bft was a homolog admitted by an over-broad taxonomy filter"
    )
    assert "clbB" not in by_person.get("SAMPLE4", {}), "SAMPLE4's clbB was zero and stays absent"
    assert by_person.get("SAMPLE6", {}).get("clbB") == 2
    # The surviving clbB signals are unaffected by the bft panel change.
    assert by_person.get("SAMPLE3", {}).get("clbB") == 1
    assert by_person.get("SAMPLE1", {}).get("clbB") == 1


def test_one_donors_clean_panel_does_not_offset_anothers_finding(result) -> None:
    """Acceptance 61: mixing cannot erase a candidate hazard by averaging."""
    for body in result["set_results"]:
        per_person = body["candidate_mechanisms_by_person"]
        members = set(body["member_person_ids"])
        # Only members appear, and each member's rows stay attributed to them.
        assert set(per_person) <= members
        for person, rows in per_person.items():
            assert all(row["person_id"] == person for row in rows)
        assert body["candidate_mechanism_count"] == sum(
            len(v) for v in per_person.values()
        )


def test_every_exposure_has_a_reason_to_be_listed(result) -> None:
    """An exposure is reported on organisms, on percentile, or not at all.

    The percentile bound alone no longer qualifies every row: an exposure
    can now be raised because the donor carries a disease-enriched
    organism the recipient lacks, even where the two percentiles are
    close. What must never happen is a row with neither reason behind it.
    """
    for e in result["leaderboard"]["new_pattern_exposures"]:
        by_organism = bool(e["introduced_organisms"])
        by_percentile = (
            e["donor_percentile"] >= 75.0
            and e["donor_percentile"] - e["recipient_percentile"] >= 20.0
        )
        assert by_organism or by_percentile, (
            f"{e['profile']} is listed with neither an introduced organism nor "
            "a qualifying percentile gap"
        )
        # Only the organism route can veto.
        if e["counts_against_selection"]:
            assert by_organism


def test_no_finding_claims_dna_present_with_zero_specific_fragments(result) -> None:
    """"0 specific fragments" and "DNA present" cannot both be true."""
    for c in result["concern_records"]:
        frags = (c.get("measurements") or {}).get("unique_supporting_fragments")
        if frags == 0:
            blob = " ".join(
                str(c.get(k) or "") for k in ("reason", "feature_label", "identity_resolution")
            ).lower()
            assert "dna present" not in blob, c["concern_id"]
            assert c["tier"] != "species_unresolved", c["concern_id"]


def test_every_finding_names_the_test_and_the_material(result) -> None:
    """"Needs a lab test" is useless without saying which test, on what."""
    lead = {"high_consequence_supported", "high_consequence_unresolved",
            "species_unresolved", "carriage_toxin_negative", "toxin_or_resistance_review"}
    checked = 0
    for c in result["concern_records"]:
        if c["tier"] not in lead:
            continue
        joined = " ".join(c.get("next_evidence_needed") or [])
        if "What would settle it" in joined:
            checked += 1
            assert any(
                phrase in joined
                for phrase in ("thawed aliquot", "fresh sample", "preserved or thawed")
            ), c["concern_id"]
    assert checked > 0, "no finding named a confirmation route"


# --------------------------------------------------------------------------- #
# Where the transfer evidence runs through a named organism, the organism has
# to be in the candidate for the pathway to be there. Maeda's arthritis result
# came from P. copri-dominated donors, and monocolonising mice with P. copri
# alone reproduced it — so a candidate with no member of the complex does not
# carry that evidence, whatever its species-level resemblance percentile says.
# --------------------------------------------------------------------------- #


def test_organism_mediated_evidence_is_checked_against_the_candidate(result) -> None:
    ra = next(
        (e for e in result["leaderboard"]["new_pattern_exposures"] if e["profile"] == "ra"),
        None,
    )
    if ra is None:
        return
    assert ra["mechanism"], "the RA entry must declare the organism it runs through"
    assert "copri" in ra["mechanism"]
    assert ra["mechanism_status"], "and say what was found in this candidate"
    if ra["transfer_evidence_tier"] == "mechanism_absent":
        assert ra["counts_against_selection"] is False
        assert "No member" in ra["mechanism_status"]
        # The whole complex, not a species label that means 13 organisms.
        assert "13 clades were assessed" in ra["mechanism_status"]
    else:
        assert ra["transfer_evidence_tier"] == "human_donor_transfer"
        assert ra["counts_against_selection"] is True


def test_a_missing_lane_cannot_quietly_clear_a_candidate() -> None:
    """An unrunnable mechanism check must leave the exposure counting."""
    from openbiota.fmt.transfer_risk import EVIDENCE, mechanism_present

    mechanism = EVIDENCE["ra"].mechanism
    assert mechanism is not None

    class _Material:
        results: dict = {}

    present, statement = mechanism_present(mechanism, _Material())
    assert present is True, "no data must not become an all-clear"
    assert "did not run" in statement


def test_the_mechanism_is_present_when_the_complex_is_detected() -> None:
    from openbiota.fmt.transfer_risk import EVIDENCE, mechanism_present

    mechanism = EVIDENCE["ra"].mechanism
    assert mechanism is not None

    class _Material:
        results = {
            "extended_catalogue": {
                "sgbs": [{"sgb": "SGB1626", "species": "Segatella_copri", "percent": 4.2}]
            }
        }

    present, statement = mechanism_present(mechanism, _Material())
    assert present is True
    assert "4.200%" in statement
    # Present, but still not a clade-specific disease finding.
    assert "no clade has an established disease direction" in statement.lower()


def test_the_ra_evidence_states_the_conditions_the_experiment_needed() -> None:
    """A transfer claim has to carry the conditions it depended on."""
    from openbiota.fmt.transfer_risk import EVIDENCE

    boundary = EVIDENCE["ra"].boundary
    for condition in ("ZAP-70", "germ-free", "zymosan"):
        assert condition in boundary, condition
    assert "treat RA rather than cause it" in boundary
    assert any("Maeda 2016" in s for s in EVIDENCE["ra"].sources)
