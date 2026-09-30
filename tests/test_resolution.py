"""The resolution layer: evidence schema, RA targets, measured mechanisms.

These tests encode the acceptance criteria from BUILD_SPEC_v07.0 §12.2 that the
implemented subset covers. They exist because every one of them describes a way
a report can state something it did not measure.
"""

from __future__ import annotations

import pytest

from openbiota.resolution import mechanisms as M
from openbiota.resolution import registry as R
from openbiota.resolution import schema as S


def _panel(**counts: int) -> dict:
    return {
        "panels": [
            {
                "name": "crc_virulence",
                "genes": [
                    {"gene": g, "fragments": n, "median_identity": 95.0,
                     "residue_check": "pass", "reference_count": 10}
                    for g, n in counts.items()
                ],
            }
        ]
    }


# --------------------------------------------------------------------------- #
# schema: the three status axes (§3.2) and acceptance test 2
# --------------------------------------------------------------------------- #


def test_the_schema_self_checks_pass() -> None:
    assert R.self_test() == 0
    assert M.self_test() == 0


def test_an_unrun_assay_can_never_be_a_negative_result() -> None:
    """Acceptance 2: an unexecuted or failed task cannot become a negative."""
    for status in (S.SCHEDULED, S.FAILED, S.REFERENCE_UNAVAILABLE,
                   S.ACCESS_UNAVAILABLE, S.INCOMPATIBLE_INPUT, S.NOT_REQUESTED):
        for call in (S.NOT_DETECTED_CALIBRATED, S.NOT_DETECTED_UNVALIDATED,
                     S.SUPPORTED_DETECTION):
            with pytest.raises(S.SchemaViolation):
                S.ResolutionCall(
                    sample_id="K1", target_id="t", identity_kind="locus",
                    assay_status=status, analytical_call=call,
                )


def test_a_calibrated_negative_must_name_its_limit() -> None:
    with pytest.raises(S.SchemaViolation, match="limit_of_detection"):
        S.ResolutionCall(
            sample_id="K1", target_id="t", identity_kind="locus", call_id="c",
            assay_status=S.COMPLETED, analytical_call=S.NOT_DETECTED_CALIBRATED,
        )
    ok = S.ResolutionCall(
        sample_id="K1", target_id="t", identity_kind="locus", call_id="c",
        assay_status=S.COMPLETED, analytical_call=S.NOT_DETECTED_CALIBRATED,
        coverage=S.Coverage(limit_of_detection=0.01, limit_scope="relative_abundance"),
    )
    assert ok.is_negative


def test_a_completed_call_must_bind_to_its_evidence() -> None:
    with pytest.raises(S.SchemaViolation, match="call_id"):
        S.ResolutionCall(
            sample_id="K1", target_id="t", identity_kind="locus",
            assay_status=S.COMPLETED, analytical_call=S.SUPPORTED_DETECTION,
        )


def test_transmission_probability_is_always_null() -> None:
    """Acceptance 70: a null must never render as 0% or as green clearance."""
    call = S.ResolutionCall(
        sample_id="K1", target_id="t", identity_kind="locus",
        assay_status=S.SCHEDULED,
    )
    assert call.to_json()["disease_transmission_probability"] is None


def test_the_census_keeps_outstanding_targets_visible() -> None:
    """Acceptance 1: no hidden placeholder outside the census."""
    calls = [
        S.ResolutionCall(sample_id="K1", target_id="a", identity_kind="locus",
                         assay_status=S.SCHEDULED),
        S.ResolutionCall(sample_id="K1", target_id="b", identity_kind="locus",
                         call_id="c", assay_status=S.COMPLETED,
                         analytical_call=S.SUPPORTED_DETECTION),
    ]
    out = S.census(calls)
    assert out["complete"] is False
    assert out["outstanding_target_ids"] == ["a"]
    assert out["targets_registered"] == 2
    assert out["targets_terminal"] == 1


def test_nulls_never_default_to_zero() -> None:
    cov = S.Coverage()
    for value in cov.to_json().values():
        assert value is None


# --------------------------------------------------------------------------- #
# RA targets (§8.1) and acceptance tests 33-37
# --------------------------------------------------------------------------- #


def test_the_ra_placeholder_is_replaced_by_four_real_targets() -> None:
    ids = {t.target_id for t in R.RA_TARGETS}
    assert ids == {
        "ra.segatella_complex", "ra.ctnpc",
        "ra.pc_p27_antigen", "ra.subdoligranulum_d8_marker",
    }


def test_clade_or_antigen_alone_is_not_an_ra_positive() -> None:
    """Acceptance 33."""
    for tid in ("ra.segatella_complex", "ra.pc_p27_antigen"):
        target = R.TARGETS[tid]
        assert target.broad_clade_is_pathogenicity is False
        assert target.validated_human_risk_predictor is False
        assert target.allow_species_fallback is False


def test_ctnpc_keeps_both_study_groups_so_origin_is_not_the_rule() -> None:
    """Acceptance 34: an RA-origin isolate that is CTnPc-negative exists."""
    t = R.TARGETS["ra.ctnpc"]
    assert "RA-N001-13" in t.positive_reference_groups
    assert "N115-17" in t.negative_reference_groups, (
        "N115-17 is RA-origin and CTnPc-negative; patient origin cannot be the call rule"
    )


def test_incomplete_curation_cannot_pass_as_a_complete_region_assay() -> None:
    """Acceptance 35.

    The region has since been delimited by presence/absence against the
    study's own comparators, so coordinates exist. Having them must not
    promote the target to a validated assay — that is the whole point of the
    criterion.
    """
    t = R.TARGETS["ra.ctnpc"]
    assert t.operational is False
    if t.canonical_region_coordinates:
        assert t.region_delimited is True
    assert t.allow_generic_mobile_element_match is False
    assert t.result_when_incomplete == "preserve_partial_evidence_and_reason"
    assert "ra.ctnpc" in R.unresolved_placeholders()


def test_d8_marker_stays_candidate_status() -> None:
    """Acceptance 36."""
    t = R.TARGETS["ra.subdoligranulum_d8_marker"]
    assert t.reference_readiness == "patent_derived_candidate"
    assert t.operational is False
    assert t.validated_human_risk_predictor is False
    assert "H3" in t.negative_reference_groups


def test_the_legacy_ra_score_is_retained_descriptively() -> None:
    """Acceptance 37: a depletion term cannot become a causal agent."""
    note = R.RA_LEGACY_SCORE_MEANING
    assert "community resemblance" in note
    assert "not a detected transferable agent" in note
    assert "vulgatus" in note, "the actual driver must be named"


def test_accessory_screens_run_even_if_the_parent_looks_absent() -> None:
    """Acceptance 17: an incomplete marker catalogue is not a blind spot."""
    for tid in ("ra.ctnpc", "ra.pc_p27_antigen"):
        assert R.TARGETS[tid].screen_independent_of_parent_profile is True


def test_every_target_binds_its_claim_to_what_was_studied() -> None:
    """Acceptance 60: a warning links to the actual measured feature."""
    for target in R.TARGETS.values():
        assert target.biological_evidence
        for evidence in target.biological_evidence:
            assert evidence.studied_entity
            assert evidence.boundary
            assert evidence.sources


# --------------------------------------------------------------------------- #
# measured mechanisms (§8.4) and acceptance tests 30, 31
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("fragments", [1, 2])
def test_a_couple_of_clbb_fragments_are_not_an_intact_pks_island(fragments: int) -> None:
    """Acceptance 31."""
    calls = {c.target_id: c for c in M.calls_for_sample("S", _panel(clbB=fragments))}
    pks = calls["mech.pks_island.clbB"]
    assert pks.analytical_call == S.CANDIDATE_SEQUENCE
    assert "locus_architecture_not_assessed" in pks.reason_codes
    assert "carrier_unresolved" in pks.reason_codes
    assert "cannot establish" in pks.plain
    assert pks.linkage.state == "sample_co_detection"


def test_a_bft_count_is_target_specific_and_unconfirmed() -> None:
    """Acceptance 30: K1's bft=11 stays a fragment count."""
    calls = {c.target_id: c for c in M.calls_for_sample("K1", _panel(bft=11, clbB=0, fadA=0))}
    bft = calls["mech.etbf.bft"]
    assert bft.analytical_call == S.CANDIDATE_SEQUENCE
    assert "subtype_not_discriminated" in bft.reason_codes
    assert bft.validation.phenotype_measured_in_sample is False
    # And the panel's other two genes are independent, not swept along.
    assert calls["mech.pks_island.clbB"].analytical_call == S.NOT_DETECTED_UNVALIDATED
    assert calls["mech.fusobacterium.fadA"].analytical_call == S.NOT_DETECTED_UNVALIDATED


def test_a_zero_count_is_not_a_calibrated_absence() -> None:
    calls = {c.target_id: c for c in M.calls_for_sample("S", _panel(fadA=0))}
    fada = calls["mech.fusobacterium.fadA"]
    assert fada.analytical_call == S.NOT_DETECTED_UNVALIDATED
    assert fada.coverage.limit_of_detection is None
    assert "detection_limit_not_calibrated" in fada.reason_codes


def test_a_missing_panel_row_is_unassessed_not_absent() -> None:
    calls = {c.target_id: c for c in M.calls_for_sample("S", {"panels": []})}
    for call in calls.values():
        assert call.assay_status == S.NOT_REQUESTED
        assert call.is_negative is False
        assert call.analytical_call == S.UNRESOLVED


def test_findings_are_per_person_and_carry_no_probability() -> None:
    findings = M.findings_for("K1", _panel(bft=11, clbB=1))
    assert findings
    for finding in findings:
        blob = finding.to_json()
        assert blob["person_id"] == "K1"
        assert blob["confirmed"] is False
        assert blob["disease_transmission_probability"] is None
        assert blob["required_for_a_confirmed_call"]


# --------------------------------------------------------------------------- #
# Targeted locus mapping (§6): the nucleotide check that can contradict a
# translated protein hit. A gene-family match at high amino-acid identity can
# be to a homolog whose gene shares no nucleotide similarity with the target;
# only mapping the exact reference can tell, and this is the stage that does.
# --------------------------------------------------------------------------- #


def test_the_locus_stage_self_checks_pass() -> None:
    from openbiota.resolution import loci

    assert loci.self_test() == 0


def test_an_unavailable_locus_stage_is_not_a_negative() -> None:
    from openbiota.resolution import loci

    call = loci.call_for(loci.BY_TARGET["locus.bft.subtype"], None, sample_id="S")
    assert call.assay_status == S.REFERENCE_UNAVAILABLE
    assert call.is_negative is False
    assert "Not a negative result" in call.plain


def test_similar_alleles_cannot_be_forced_into_a_subtype() -> None:
    """bft1/2/3 are ~1.5 kb and highly similar; a thin margin names nothing."""
    from openbiota.resolution import loci

    target = loci.BY_TARGET["locus.bft.subtype"]

    def cov(name: str, covered: int) -> loci.RefCoverage:
        n = 1546 // loci.WINDOW_BP + 1
        return loci.RefCoverage(name, 1546, covered, covered * 10, (1,) * n, 300)

    tie = {
        "mech.bft1": cov("mech.bft1", 1400),
        "mech.bft2": cov("mech.bft2", 1380),
        "mech.bft3": cov("mech.bft3", 1370),
    }
    call = loci.call_for(target, tie, sample_id="S")
    assert "subtype_not_discriminated" in call.reason_codes
    assert call.coverage.discriminatory_bases == 0
    assert "not enough to name one" in call.plain


def test_scattered_reads_are_not_a_biosynthetic_island() -> None:
    """Acceptance 31 at nucleotide level: breadth, not read count, decides."""
    from openbiota.resolution import loci

    target = loci.BY_TARGET["locus.pks_island.architecture"]
    n = 55140 // loci.WINDOW_BP + 1
    scattered = {
        "mech.pks_island": loci.RefCoverage(
            "mech.pks_island", 55140, 1200, 1200, (1,) + (0,) * (n - 1), 20
        )
    }
    call = loci.call_for(target, scattered, sample_id="S")
    assert call.analytical_call == S.CANDIDATE_SEQUENCE
    assert "scattered_fragments_only" in call.reason_codes
    assert "not a locus" in call.plain
    # Even broad coverage never claims per-gene integrity.
    broad = {
        "mech.pks_island": loci.RefCoverage(
            "mech.pks_island", 55140, 50000, 400000, (1,) * n, 3000
        )
    }
    call = loci.call_for(target, broad, sample_id="S")
    assert call.analytical_call == S.SUPPORTED_DETECTION
    assert "per_gene_integrity_unassessed" in call.reason_codes


def test_a_contradicted_protein_hit_becomes_ambiguous_homology() -> None:
    """The K1 bft case: 11 protein hits, zero nucleotide match."""
    payload = dict(
        _panel(bft=11),
        locus_resolution={"targets": [{
            "target_id": "locus.bft.subtype",
            "assay_status": S.COMPLETED,
            "analytical_call": S.NOT_DETECTED_UNVALIDATED,
        }]},
    )
    calls = {c.target_id: c for c in M.calls_for_sample("K1", payload)}
    bft = calls["mech.etbf.bft"]
    assert bft.analytical_call == S.AMBIGUOUS_HOMOLOGY
    assert "nucleotide_check_found_no_match" in bft.reason_codes
    assert "match is to a homolog" in bft.plain
    # It stays visible rather than vanishing from the ledger.
    findings = [f for f in M.findings_for("K1", payload) if f.gene == "bft"]
    assert findings
    assert findings[0].to_json()["resolved_as_homolog"] is True


def test_a_corroborated_protein_hit_says_so() -> None:
    payload = dict(
        _panel(clbB=40),
        locus_resolution={"targets": [{
            "target_id": "locus.pks_island.architecture",
            "assay_status": S.COMPLETED,
            "analytical_call": S.SUPPORTED_DETECTION,
        }]},
    )
    calls = {c.target_id: c for c in M.calls_for_sample("K1", payload)}
    pks = calls["mech.pks_island.clbB"]
    assert "nucleotide_check_corroborates" in pks.reason_codes
    assert "nucleotide check agrees" in pks.plain
    assert pks.analytical_call == S.CANDIDATE_SEQUENCE
