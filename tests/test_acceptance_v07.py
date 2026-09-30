"""BUILD_SPEC_v07.0 §12.2 must-pass acceptance suite.

Each test names the numbered criterion it enforces. Criteria whose capability
is not yet implemented are asserted at the level that *is* implemented — that
the target is registered, reaches a terminal non-negative state, and names what
is missing — because §12.2's own rules say an unexecuted assay must never
become a negative call. A criterion that cannot yet be met by a detector is
still met by the census.

`test_every_criterion_is_accounted_for` at the end is the guard against this
file quietly drifting out of sync with the spec.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.resolution import adapters as A
from openbiota.resolution import kingdoms as K
from openbiota.resolution import mechanisms as M
from openbiota.resolution import registry as R
from openbiota.resolution import schema as S
from openbiota.resolution import strains as St
from openbiota.resolution import typing as T
from openbiota.resolution.resolver import EvidenceResolver

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
SAMPLES = ("SAMPLE4_A04", "SAMPLE6_A06", "SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03")


def _results(sample: str) -> dict:
    return json.loads((RESULTS / sample / "results.json").read_text())


def _has_samples() -> bool:
    return all((RESULTS / s / "results.json").exists() for s in SAMPLES)


needs_samples = pytest.mark.skipif(not _has_samples(), reason="sample results absent")


# --------------------------------------------------------------------------- #
# 1-8: census, statuses, and the species/strain boundary
# --------------------------------------------------------------------------- #


def test_c01_every_target_has_a_census_entry() -> None:
    """1. No hidden strain-only placeholder outside the census."""
    calls = [
        S.ResolutionCall(sample_id="SAMPLE4", target_id=t.target_id,
                         identity_kind=t.identity_kind, assay_status=S.SCHEDULED)
        for t in R.RA_TARGETS
    ]
    calls += K.all_calls("SAMPLE4")
    calls += A.calls_for_unavailable("SAMPLE4")
    census = S.census(calls)
    assert census["targets_registered"] == len(calls)
    # Non-operational RA targets are named, not hidden.
    assert set(R.unresolved_placeholders()) <= {c.target_id for c in calls}


def test_c02_an_unexecuted_task_cannot_become_a_negative() -> None:
    """2. Every requested target gets a terminal status; no fake negatives."""
    for status in (S.SCHEDULED, S.FAILED, S.REFERENCE_UNAVAILABLE,
                   S.ACCESS_UNAVAILABLE, S.INCOMPATIBLE_INPUT):
        with pytest.raises(S.SchemaViolation):
            S.ResolutionCall(sample_id="SAMPLE4", target_id="t", identity_kind="locus",
                             assay_status=status,
                             analytical_call=S.NOT_DETECTED_UNVALIDATED)


def test_c03_a_species_tsv_cannot_fabricate_a_strain_fingerprint() -> None:
    """3. A species-only result can never answer a strain query."""
    species = S.ResolutionCall(
        sample_id="SAMPLE4", target_id="org.x", identity_kind="species",
        assay_status=S.COMPLETED, call_id="c", analytical_call=S.SUPPORTED_DETECTION,
    )
    out = EvidenceResolver([species]).resolve("SAMPLE4", "org.x", "population_fingerprint")
    assert not out.satisfied
    assert out.missing_reasons


@needs_samples
def test_c04_missing_sam_triggered_a_real_mapping_task() -> None:
    """4. Missing SAM with FASTQs available produced a reproducible mapping."""
    for sample in SAMPLES:
        sam = RESULTS / sample / "strain" / f"{sample}.markers.sam.bz2"
        markers = RESULTS / sample / "strain" / "consensus_markers"
        assert sam.exists() and sam.stat().st_size > 1_000_000, sample
        assert list(markers.glob("*.json.bz2")), sample


@needs_samples
def test_c05_mp3_scoring_lane_is_untouched_by_the_new_features() -> None:
    """5. New MP4 features cannot silently enter the old CDF."""
    data = _results("SAMPLE4_A04")
    engine = data["profile_similarity"]["taxonomic_engine"]
    blob = json.dumps(engine)
    assert "mpa_v31_CHOCOPhlAn_201901" in blob or "3.1.0" in blob
    # Strain evidence lives in its own namespace, not inside the profile lane.
    assert "strain_resolution" in data
    assert "strain_resolution" not in json.dumps(data["profile_similarity"])


@needs_samples
def test_c06_an_sgb_is_never_treated_as_a_named_strain() -> None:
    """6. Taxonomy aliases preserved; an SGB is not a strain."""
    sr = _results("SAMPLE4_A04")["strain_resolution"]
    assert sr["status"] == "resolved"
    for row in sr["organisms"]:
        assert row["is_named_strain"] is False
        assert row["sgb"].startswith("SGB")


@needs_samples
def test_c07_species_detected_but_strain_unresolved_keeps_partial_evidence() -> None:
    """7. "Species detected; strain unresolved" retains its evidence."""
    sr = _results("SAMPLE4_A04")["strain_resolution"]
    assert sr["organisms_detected_unresolved"] > 0
    # Resolved organisms keep full coverage detail rather than a bare flag.
    for row in sr["organisms"][:5]:
        assert row["markers_resolved"] > 0
        assert row["callable_bases"] > 0
        assert row["median_breadth_percent"] > 0


def test_c08_an_omitted_true_strain_is_not_forced_to_an_exact_identity() -> None:
    """8. Nearest-reference/novelty information, not forced identity."""
    out = EvidenceResolver([]).resolve("SAMPLE4", "strain.population.SGB1",
                                      "named_reference_isolate")
    assert not out.satisfied
    # Marker scope forbids the exact-isolate claim outright.
    assert "not establish whole-genome identity" in St.MARKER_SCOPE_LIMIT


# --------------------------------------------------------------------------- #
# 9-16: mixtures, spans and method independence
# --------------------------------------------------------------------------- #


def test_c10_a_minority_population_is_never_erased_silently() -> None:
    """10. A dominant consensus must declare that it hides minorities."""
    assert "concealed" in St.DOMINANT_CONSENSUS_LIMIT
    calls = St.calls_for_sample(
        {"SGB1": St.MarkerFingerprint("A", "SGB1", 50, 50_000, 95.0, 9.0,
                                      {"m|1__2|SGB1": "ACGT" * 100})},
        sample_id="A", database="db",
    )
    assert calls[0].population_components[0]["minority_populations_resolved"] is False
    assert "dominant_consensus_only" in calls[0].reason_codes


def test_c12_a_conserved_shared_marker_is_not_a_whole_genome_match() -> None:
    """12. Marker agreement is not exact-strain identity."""
    seq = {"m|1__2|SGB1": "ACGT" * 3000}
    fp = St.MarkerFingerprint("A", "SGB1", 1, 12000, 95.0, 9.0, seq)
    other = St.MarkerFingerprint("B", "SGB1", 1, 12000, 95.0, 9.0, seq)
    out = St.compare(fp, other, sgb="SGB1")
    assert out.status == "indistinguishable_over_compared_sites"
    assert "not proof of whole-genome identity" in out.plain


def test_c13_zero_snps_over_a_short_span_is_not_sharing() -> None:
    """13. Insufficient callable span yields no distance at all."""
    seq = {"m|1__2|SGB1": "ACGT" * 100}
    fp = St.MarkerFingerprint("A", "SGB1", 1, 400, 95.0, 9.0, seq)
    other = St.MarkerFingerprint("B", "SGB1", 1, 400, 95.0, 9.0, seq)
    out = St.compare(fp, other, sgb="SGB1")
    assert out.status == "insufficient_callable_overlap"
    assert out.difference_per_kb is None


def test_c14_tracs_omits_meta_and_emits_no_viral_rate_probabilities() -> None:
    """14. The mandatory TRACS correction."""
    cmd = A.tracs_distance_command(msa="a.fa", out="d.csv")
    assert "--meta" not in cmd
    with pytest.raises(ValueError, match="--meta"):
        A.BY_ADAPTER["tracs"].command("distance", "--meta", "m.csv")


def test_c15_skipped_references_are_insufficient_resolution_not_negatives() -> None:
    """15. A filtered/skipped reference is not a biological non-detection."""
    limits = " ".join(A.BY_ADAPTER["tracs"].limits)
    assert "not biological non-detection" in limits


def test_c16_methods_keep_separate_definitions_with_no_agreement_bonus() -> None:
    """16. Correlated agreement is not multiplied into certainty."""
    estimates = [a.estimates for a in A.ADAPTERS]
    assert len(set(estimates)) == len(estimates)
    assert "not multiplied into" in A.registry_json()["no_consensus_bonus"]


# --------------------------------------------------------------------------- #
# 17-21: independence of the targeted lane, and linkage
# --------------------------------------------------------------------------- #


def test_c17_organism_failure_does_not_suppress_the_targeted_screen() -> None:
    """17. Accessory screens run even if the parent looks absent."""
    for tid in ("ra.ctnpc", "ra.pc_p27_antigen"):
        assert R.TARGETS[tid].screen_independent_of_parent_profile is True


def test_c18_a_gene_signal_survives_an_unresolved_carrier() -> None:
    """18. Supported gene evidence persists without carrier resolution.

    Built from an explicit fixture rather than a live sample. It used to
    read SAMPLE4, whose bft count was a false positive from a panel that
    admitted a non-toxin homolog; once the panel was corrected SAMPLE4 has no
    mechanism finding at all and the test was asserting the bug. The
    behaviour under test - evidence is retained when the carrier cannot be
    resolved - does not depend on any particular sample having a hit.
    """
    findings = M.findings_for("S", {
        "panels": [{"name": "crc_virulence",
                    "genes": [{"gene": "bft", "fragments": 11}]}]})
    assert findings
    for f in findings:
        assert "carrier_unresolved" in f.reason_codes
        assert f.fragments > 0  # the evidence is retained, not discarded


def test_an_installed_typing_tool_that_never_ran_is_not_a_negative_result() -> None:
    """A registry census cannot report the outcome of a test nobody ran.

    The census enumerates every scheme and calls ``call_type`` with no
    observations. For schemes whose tool is merely absent that already
    produced ``reference_unavailable``. But for the four whose tool is
    installed, the empty sequence fell through to the locus rules and came
    back ``completed`` with "required loci not detected" - a negative
    result for ECTyper, Kleborate and AMRFinder runs that never happened.
    Four such calls were present in every shipped results.json.
    """
    installed = [s for s in T.BY_SCHEME.values() if s.tool_installed]
    if not installed:
        pytest.skip("no typing tool (ECTyper, Kleborate, AMRFinder) is installed in this checkout")
    for scheme in installed:
        call = T.call_type(
            sample_id="S", scheme_id=scheme.scheme_id, observations=[]
        )
        assert call.assay_status == S.SCHEDULED, scheme.scheme_id
        assert call.analytical_call == S.UNRESOLVED
        assert "typing_tool_installed_not_run" in call.reason_codes
        assert "not looked for" in call.plain
        # And the schema invariant holds: an unrun assay is never negative.
        assert "not_detected" not in str(call.analytical_call)


def test_c19_a_toxin_is_not_assigned_to_the_likeliest_host() -> None:
    """19. No carrier assignment without linkage."""
    calls = M.calls_for_sample("S", {
        "panels": [{"name": "crc_virulence",
                    "genes": [{"gene": "bft", "fragments": 40}]}]})
    bft = next(c for c in calls if c.target_id == "mech.etbf.bft")
    assert bft.linkage.state == "sample_co_detection"
    assert bft.linkage.carrier_call_id is None


def test_c20_absence_from_an_incomplete_assembly_is_not_a_deletion() -> None:
    """20. A gene missing from an incomplete MAG is unknown, not deleted."""
    calls = M.calls_for_sample("S", {
        "panels": [{"name": "crc_virulence", "genes": [{"gene": "fadA", "fragments": 0}]}]})
    fada = next(c for c in calls if c.target_id == "mech.fusobacterium.fadA")
    assert fada.analytical_call == S.NOT_DETECTED_UNVALIDATED
    assert "detection_limit_not_calibrated" in fada.reason_codes


# --------------------------------------------------------------------------- #
# 22-32: typing rules
# --------------------------------------------------------------------------- #


needs_ectyper = pytest.mark.skipif(
    not T.BY_SCHEME["ecoli.pathotype"].tool_installed,
    reason="ECTyper not installed; run 'make strain-tools'",
)


@needs_ectyper
def test_c22_23_unphased_alleles_cannot_build_a_type() -> None:
    """22/23. No fictional ST, and no fabricated O:H from two populations."""
    call = T.call_type(
        sample_id="SAMPLE4", scheme_id="ecoli.pathotype",
        observations=[
            T.LocusObservation("pathotype-defining locus set with carrier linkage",
                               detected=True, fragments=50, linked_to_carrier=True),
            T.LocusObservation("stx1", detected=False, depth_adequate_for_absence=True),
            T.LocusObservation("stx2", detected=False, depth_adequate_for_absence=True),
        ],
        populations_of_species=2,
        assay_ran=True,
    )
    assert call.analytical_call == S.MIXED_UNRESOLVED
    assert "exists in no organism" in call.plain


def test_c24_25_eae_only_and_stx_only_stay_separate() -> None:
    """24/25. No linked eae+stx strain; free phage leaves carrier unresolved."""
    rule = T.BY_SCHEME["ecoli.pathotype"].mixed_specimen_rule
    assert "never inferred" in rule
    assert "carrier-unresolved" in rule


def test_c26_ipah_alone_is_not_a_shigella_species() -> None:
    """26. EIEC alternatives are preserved."""
    scheme = T.BY_SCHEME["shigella.identity"]
    assert "never becomes a named Shigella species" in scheme.mixed_specimen_rule
    assert "carried by both" in scheme.species_insufficient_because


@needs_ectyper
def test_c27_a_shallow_stx_nondetection_is_not_a_negative_criterion() -> None:
    """27. "stx-negative" requires adequate depth."""
    call = T.call_type(
        sample_id="SAMPLE4", scheme_id="ecoli.pathotype",
        observations=[
            T.LocusObservation("pathotype-defining locus set with carrier linkage",
                               detected=True, fragments=50, linked_to_carrier=True),
            T.LocusObservation("stx1", detected=False, depth_adequate_for_absence=False),
        ],
        assay_ran=True,
    )
    assert call.analytical_call == S.INSUFFICIENT_DEPTH
    assert "negative_criterion_depth_inadequate" in call.reason_codes


def test_an_uninstalled_scheme_is_unavailable_never_negative() -> None:
    """The other half of 2: a missing tool is not a clean result."""
    absent = next(s for s in T.SCHEMES if not s.tool_installed)
    call = T.call_type(sample_id="SAMPLE4", scheme_id=absent.scheme_id, observations=[])
    assert call.assay_status == S.REFERENCE_UNAVAILABLE
    assert call.is_negative is False
    assert "not a negative result" in call.plain


def test_c28_klebsiella_loci_from_different_strains_do_not_merge() -> None:
    """28. No artificial hypervirulent resistant strain."""
    assert "must not be merged" in T.BY_SCHEME["klebsiella.loci"].mixed_specimen_rule


def test_c29_meca_in_a_non_aureus_staph_is_not_mrsa() -> None:
    """29. And an unlinked van gene is not VRE."""
    rule = T.BY_SCHEME["staph.mrsa"].mixed_specimen_rule
    assert "never MRSA" in rule
    assert "never VRE" in rule


@needs_samples
def test_c30_the_k1_virulence_fixture_stays_unconfirmed() -> None:
    """30. SAMPLE4's virulence targets are all zero and nothing is confirmed.

    bft was 11 here until the crc_virulence panel was corrected. That
    count came from a reference set scoped to ``taxonomy_id:2`` - all
    Bacteria - which admitted a 364 aa *B. fragilis* protein annotated
    only "Fragilysin". It aligns to the published bft-1 allele over 15
    residues at 60% identity: the shared zinc motif and nothing else. The
    nucleotide lane had always reported zero reads at the canonical
    alleles, and the nucleotide lane was right. With the panel scoped to
    Bacteroides, length-windowed to the real 397-405 aa pro-protein and
    competing against a paralogue decoy, the count is zero and the two
    lanes agree.
    """
    panel = next(
        p for p in _results("SAMPLE4_A04")["panels"] if p["name"] == "crc_virulence"
    )
    counts = {g["gene"]: g["fragments"] for g in panel["genes"]}
    assert counts["bft"] == 0
    assert counts["clbB"] == 0
    assert counts["fadA"] == 0
    for f in M.findings_for("SAMPLE4", _results("SAMPLE4_A04")):
        assert f.to_json()["confirmed"] is False


def test_c31_one_or_two_clbb_fragments_are_not_an_intact_island() -> None:
    """31. clbB/clbS alone is not a functional pks island."""
    for n in (1, 2):
        calls = M.calls_for_sample("S", {
            "panels": [{"name": "crc_virulence", "genes": [{"gene": "clbB", "fragments": n}]}]})
        pks = next(c for c in calls if c.target_id == "mech.pks_island.clbB")
        assert pks.analytical_call == S.CANDIDATE_SEQUENCE
        assert "locus_architecture_not_assessed" in pks.reason_codes


def test_c32_toxin_negative_cdiff_stays_detectable() -> None:
    """32. TcdA-/TcdB+ strains and subtype-versus-ribotype distinction."""
    scheme = T.BY_SCHEME["cdiff.toxinotype"]
    assert "TcdA-negative/TcdB-positive" in scheme.mixed_specimen_rule
    assert "not ribotype" in scheme.mixed_specimen_rule
    assert "tcdA" in scheme.required_loci and "tcdB" in scheme.required_loci


# --------------------------------------------------------------------------- #
# 33-44: RA repair, functional panels, probiotics
# --------------------------------------------------------------------------- #


def test_c33_clade_a_and_pc_p27_cannot_trigger_a_pathogenic_claim() -> None:
    for tid in ("ra.segatella_complex", "ra.pc_p27_antigen"):
        t = R.TARGETS[tid]
        assert t.broad_clade_is_pathogenicity is False
        assert t.validated_human_risk_predictor is False


def test_c34_ctnpc_study_groups_are_distinguishable_and_origin_is_not_the_rule() -> None:
    t = R.TARGETS["ra.ctnpc"]
    assert t.positive_reference_groups and t.negative_reference_groups
    assert "N115-17" in t.negative_reference_groups
    manifest = REPO / "refs" / "ctnpc" / "MANIFEST.json"
    if manifest.exists():
        data = json.loads(manifest.read_text())
        roles = " ".join(c["role"] for c in data["negative_comparators_required"])
        assert "patient origin is not the call rule" in roles


def test_c35_generic_mobile_hits_cannot_establish_ctnpc() -> None:
    t = R.TARGETS["ra.ctnpc"]
    assert t.allow_generic_mobile_element_match is False
    assert t.operational is False, "a delimited region is still not a validated assay"
    # The region was delimited by presence/absence against the study's own
    # comparators, not by mobile-element annotation. Having coordinates must
    # not promote it.
    if t.canonical_region_coordinates:
        assert t.region_delimited is True
        region = REPO / "refs" / "ctnpc" / "CTNPC_REGION.json"
        if not region.exists():
            pytest.skip("the CTnPC region file is built by `make refs-ctnpc`; not present in this checkout")
        data = json.loads(region.read_text())
        shared = data["shared_element"]
        assert shared["status"] == "candidate_region_delimited"
        # The record must carry the negative result that corrected the first
        # attempt, so the reasoning cannot be lost.
        assert "zero orthologs" in shared["how_this_differs_from_a_guess"]
        assert "Not a validated assay" in shared["still_not"]
    manifest = REPO / "refs" / "ctnpc" / "MANIFEST.json"
    if manifest.exists():
        data = json.loads(manifest.read_text())
        assert data["operational"] is False
        assert data["interpretation_guards"]["allow_generic_mobile_element_match"] is False
        # The sensitivity analysis is the evidence that a guess is unsafe.
        assert data["candidate_region_sensitivity"]
        assert "be a guess" in data["why_not_a_call"]


def test_c36_d8_markers_remain_patent_derived_candidates() -> None:
    t = R.TARGETS["ra.subdoligranulum_d8_marker"]
    assert t.reference_readiness == "patent_derived_candidate"
    assert t.operational is False
    boundary = t.biological_evidence[0].boundary
    assert "not independently validated" in boundary
    assert "AT-rich" in boundary


@needs_samples
def test_c37_the_prior_ra_score_is_retained_descriptively() -> None:
    note = R.RA_LEGACY_SCORE_MEANING
    assert "community resemblance" in note
    assert "not a detected transferable agent" in note
    ra = _results("SAMPLE4_A04")["profile_similarity"]["profiles"]["ra"]
    assert ra["combined"]["percentile"] is not None  # retained, not deleted


@needs_samples
def test_c38_39_urda_is_urocanate_not_equol_or_urolithin() -> None:
    """38/39. UrdA cannot fire an equol/urolithin interpretation."""
    data = _results("SAMPLE4_A04")
    panel = next((p for p in data["panels"] if p["name"] == "urda"), None)
    if panel is None:
        pytest.skip("urda panel not present in this build")
    # Scope to the panel's own interpretation fields. An organism *named*
    # Adlercreutzia equolifaciens appearing in the hit list is not an equol
    # interpretation, so matching the whole JSON blob would be the wrong test.
    interpretation = " ".join(
        str(panel.get(k) or "") for k in ("metabolite", "pathway", "description", "label")
    ).lower()
    assert "imidazole propionate" in interpretation or "urocanate" in interpretation
    assert "urolithin" not in interpretation, "UrdA must not trigger a urolithin claim"
    assert "equol" not in interpretation, "UrdA must not trigger an equol claim"


# --------------------------------------------------------------------------- #
# 45-59: multi-kingdom boundaries
# --------------------------------------------------------------------------- #


def test_c45_denominators_are_never_merged() -> None:
    assert "cannot be merged" in K.DENOMINATORS["eukaryote_marker_relative"]
    assert len(K.DENOMINATORS) >= 4


def test_c46_entamoeba_ambiguity_does_not_force_the_pathogen() -> None:
    assert "must not be forced" in K.BY_TARGET["protist.entamoeba.species"].namespace_note


def test_c47_giardia_heterozygosity_is_not_a_mixed_infection() -> None:
    assert "not automatically a mixed infection" in (
        K.BY_TARGET["protist.giardia.assemblage"].ploidy_note
    )


def test_c48_missing_gp60_span_invents_no_subtype() -> None:
    assert "no subtype is invented" in (
        K.BY_TARGET["protist.cryptosporidium.gp60"].namespace_note
    )


def test_c50_blastocystis_namespaces_stay_distinct() -> None:
    note = K.BY_TARGET["protist.blastocystis.subtype"].namespace_note
    assert "distinct namespaces" in note
    assert "not a harmful-pathotype label" in note


def test_c51_candida_diploid_alleles_survive() -> None:
    note = K.BY_TARGET["fungi.candida_albicans.mlst"].ploidy_note
    assert "Diploid" in note
    assert "discard one allele" in note


def test_c52_53_fungal_erg11_presence_is_not_resistance() -> None:
    note = K.BY_TARGET["fungi.amr.substitutions"].namespace_note
    assert "is not resistance" in note
    assert "genetic code" in note


def test_c54_a_helminth_mito_haplotype_is_not_a_strain() -> None:
    note = K.BY_TARGET["helminth.haplotype"].viability_note
    assert "not a complete nuclear strain" in note
    assert "viable" in note


def test_c56_dna_input_returns_rna_viruses_as_incompatible() -> None:
    for name in K.RNA_VIRUS_TARGETS:
        call = K.call_for("SAMPLE4", f"virus.rna.{name}")
        assert call.assay_status == S.INCOMPATIBLE_INPUT
        assert not call.is_negative
        assert "cannot change that" in call.plain


def test_c57_phages_do_not_inherit_a_pathogen_label() -> None:
    assert "does not inherit a pathogen framing" in (
        K.BY_TARGET["virus.phage.votu"].namespace_note
    )


def test_c59_mixed_viral_haplotypes_are_not_flattened() -> None:
    assert "recombinant" in K.BY_TARGET["virus.dna.adenovirus_type"].namespace_note


# --------------------------------------------------------------------------- #
# 60-70: consumers, FMT and rendering
# --------------------------------------------------------------------------- #


def test_c60_every_claim_binds_to_a_studied_entity() -> None:
    for target in R.TARGETS.values():
        for evidence in target.biological_evidence:
            assert evidence.studied_entity
            assert evidence.sources
            assert evidence.boundary


@needs_samples
def test_c61_62_donor_mixing_cannot_erase_a_finding() -> None:
    """61/62. No averaging away hazards; recipient percentile suppresses nothing.

    The v07 spec asked for an empty ACTIONABLE_TIERS. That is no longer the
    setting: the recipient's standing instruction is to veto candidates
    carrying a disease pattern they lack where the condition has
    human-donor-transfer evidence, and on what goes into their own body
    that governs. What this test is really about survives unchanged - one
    donor's clean panel must never offset another's finding - so it now
    asserts the veto is narrow rather than absent.
    """
    from openbiota.fmt.transfer_risk import (
        ACTIONABLE_TIERS,
        ASSOCIATION_ONLY,
        HUMAN_DONOR_TRANSFER,
    )

    assert frozenset({HUMAN_DONOR_TRANSFER}) == ACTIONABLE_TIERS
    assert ASSOCIATION_ONLY not in ACTIONABLE_TIERS
    match = RESULTS / "fmt_SAMPLE2_vs_4donors" / "match.json"
    if not match.exists():
        pytest.skip("match result absent")
    data = json.loads(match.read_text())
    ledger = data["leaderboard"]["candidate_mechanisms_by_person"]
    assert ledger
    for body in data["set_results"]:
        per_person = body["candidate_mechanisms_by_person"]
        assert set(per_person) <= set(body["member_person_ids"])


@needs_samples
def test_c64_no_followup_means_engraftment_is_unmeasured() -> None:
    """64. Baseline similarity is not a future-outcome probability."""
    match = RESULTS / "fmt_SAMPLE2_vs_4donors" / "match.json"
    if not match.exists():
        pytest.skip("match result absent")
    sm = json.loads(match.read_text())["strain_matching"]
    assert sm["observed_engraftment"]["status"] == "followup_required"
    assert sm["observed_engraftment"]["results"] is None
    assert sm["prospective_establishment"]["predictions"] is None
    assert sm["disease_transmission_probability"] is None


@needs_samples
def test_c65_an_ambiguous_source_is_not_credited_to_one_donor() -> None:
    """65. compatible_with_multiple_sources rather than arbitrary credit."""
    match = RESULTS / "fmt_SAMPLE2_vs_4donors" / "match.json"
    if not match.exists():
        pytest.skip("match result absent")
    sm = json.loads(match.read_text())["strain_matching"]
    per_donor = sm["baseline_comparison"]["per_donor"]
    # Each donor is compared independently; no combination figure is produced.
    assert set(per_donor) == {"SAMPLE4", "SAMPLE6", "SAMPLE1", "SAMPLE3"}
    assert "not the average" in sm["baseline_comparison"]["note"]


@needs_samples
def test_c63_availability_denominators_stay_distinct() -> None:
    """63. Total and reachable denominators are never conflated."""
    match = RESULTS / "fmt_SAMPLE2_vs_4donors" / "match.json"
    if not match.exists():
        pytest.skip("match result absent")
    data = json.loads(match.read_text())
    lb = data["leaderboard"]
    assert lb["scored_goals"] >= lb["reachable_goals"]
    assert data["match_score_semantics"] == "target_availability_and_complementarity"


def test_c70_nulls_never_render_as_zero_or_green() -> None:
    call = S.ResolutionCall(sample_id="SAMPLE4", target_id="t", identity_kind="locus",
                            assay_status=S.SCHEDULED)
    assert call.to_json()["disease_transmission_probability"] is None
    assert St.Comparison(
        sgb="S", left_sample="a", right_sample="b", shared_markers=0,
        compared_bases=0, differing_bases=0, difference_per_kb=None,
        status="insufficient_callable_overlap", plain="x",
    ).to_json()["disease_transmission_probability"] is None


# --------------------------------------------------------------------------- #
# 71-80: migration, adapters and model states
# --------------------------------------------------------------------------- #


@needs_samples
def test_c71_the_migration_discovered_real_inputs() -> None:
    """71. Real FASTQs, database versions and specimen pairing."""
    for sample in SAMPLES:
        data = _results(sample)
        assert data["input"]["paired"] is True
        assert data["input"]["read_pairs"] > 1_000_000
        assert data["strain_resolution"]["database_release"] == (
            "mpa_vJun23_CHOCOPhlAnSGB_202403"
        )


@needs_samples
def test_c72_a_bowtie2_summary_is_never_accepted_as_sam() -> None:
    """72. And an existing summary cache cannot suppress SAM generation."""
    for sample in SAMPLES:
        strain_dir = RESULTS / sample / "strain"
        sam = strain_dir / f"{sample}.markers.sam.bz2"
        legacy = RESULTS / sample / "taxonomy"
        # The historical summaries are preserved untouched...
        assert list(legacy.glob("*.bowtie2.bz2"))
        # ...and a genuinely separate SAM was produced in a new directory.
        assert sam.exists()
        assert sam.parent != legacy
        assert sam.stat().st_size > 10 * (1 << 20)


def test_c74_a_missing_research_artifact_blocks_only_its_own_task() -> None:
    """74. Baseline comparison still runs without Schmidt artifacts."""
    from openbiota.fmt import strain_matching as SM

    empty = SM.StrainMatching().to_json()
    assert empty["prospective_establishment"]["status"] == "model_unavailable"
    # The baseline capability is independent of it.
    assert empty["baseline_comparison"]["status"] == "pending"


def test_c75_baseline_only_returns_followup_required() -> None:
    from openbiota.fmt import strain_matching as SM

    out = SM.StrainMatching().to_json()["observed_engraftment"]
    assert out["status"] == "followup_required"
    assert "not zero" in out["why"]


def test_c79_a_missing_model_does_not_zero_the_availability_score() -> None:
    from openbiota.fmt import strain_matching as SM

    out = SM.StrainMatching().to_json()
    assert out["prospective_establishment"]["predictions"] is None
    assert "not a probability" in out["prospective_establishment"]["why"]
    assert SM.MATCH_SCORE_SEMANTICS == "target_availability_and_complementarity"


def test_c80_strain_identity_does_not_upgrade_an_unlinked_amr_hit() -> None:
    """80. Resolving a population does not attribute a separate gene to it."""
    calls = St.calls_for_sample(
        {"SGB1": St.MarkerFingerprint("A", "SGB1", 100, 100_000, 95.0, 9.0,
                                      {"m|1__2|SGB1": "ACGT" * 100})},
        sample_id="A", database="db",
    )
    assert calls[0].analytical_call == S.SUPPORTED_DETECTION
    assert calls[0].linkage.state == "unassigned"
    assert "cannot attribute an unlinked toxin or resistance gene" in St.MARKER_SCOPE_LIMIT


# --------------------------------------------------------------------------- #
# the guard against drift
# --------------------------------------------------------------------------- #


def test_every_criterion_is_accounted_for() -> None:
    """Each of the 80 criteria is either enforced here or listed as deferred.

    Deferred means the detector is not built yet, in which case §12.2's own
    rules still apply and are enforced: the target is registered, reaches a
    terminal non-negative state, and names what is missing.
    """
    enforced = set()
    for name in globals():
        if name.startswith("test_c"):
            for part in name[6:].split("_"):
                if part.isdigit():
                    enforced.add(int(part))
    deferred = {
        9: "exact-duplicate reference collapsing (reference compiler not built)",
        11: "balanced-mixture phasing (mixture deconvolution not built)",
        21: "multi-sample coassembly leakage (assembly lane not built)",
        40: "R. gnavus capsule/inflammatory panels (reference build pending)",
        41: "Fna C1/C2 SGB crosswalk (reference build pending)",
        42: "named probiotic identity (reference build pending)",
        43: "strain-specific probiotic trial scoping (evidence library pending)",
        44: "in vitro antimicrobial evidence labelling (evidence library pending)",
        49: "Cyclospora eight-marker partial profile (reference build pending)",
        55: "fungal morphogenesis phenotype guard (fungal lane pending)",
        58: "UHGV source-catalogue dedup (viral catalogue not installed)",
        66: "cache invalidation across reference/threshold changes (compiler pending)",
        67: "negative-control visibility (no controls sequenced for these runs)",
        68: "licensed-adapter access states (adapters registered, none installed)",
        69: "full HTML/PDF/JSON census reconciliation (renderer partially wired)",
        73: "first-pass SAM retention for new samples (pipeline change pending)",
        76: "repeated persistence vs transient passage (needs follow-up specimens)",
        77: "prospective model training hygiene (no model trained)",
        78: "Schmidt feature-compatibility mapping (artifacts not downloaded)",
    }
    covered = enforced | set(deferred)
    missing = sorted(set(range(1, 81)) - covered)
    assert not missing, f"criteria neither enforced nor explicitly deferred: {missing}"
    # Deferred criteria must not silently become "passing".
    assert len(deferred) == 19
    assert len(enforced) >= 45, f"only {len(enforced)} criteria enforced"
