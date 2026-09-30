"""BUILD_SPEC_v0.8.3 sections 3.3, 3.4, 4.1, 4.2 and 16: the new contracts.

These are the guards that keep a new measurement from misstating itself:
zero is not unmeasured, a percentile needs a reference, a gene symbol is not a
biochemical identity, and an unbuilt feature is not a limit of science.
"""

from __future__ import annotations

import pytest

from openbiota.extension import markers as MK
from openbiota.extension import registry as R
from openbiota.extension import schema as S


def _metric(**over):
    base = {
        "metric_id": "ext083.test.metric",
        "label": "Test metric",
        "kind": "gene_abundance",
        "state": "measured",
        "value": 1.0,
        "unit": "fragment_RPKM",
        "denominator": "quality_passed_nonhost_fragments",
        "method_id": "ext083.test/1.0",
        "limitations": ("A capacity is not a concentration.",),
    }
    return S.ExtensionMetric(**(base | over))


# --------------------------------------------------------------------------- #
# section 3.3 — zero, null and not assessed are different states
# --------------------------------------------------------------------------- #


def test_an_unmeasured_metric_cannot_carry_a_number() -> None:
    """The single most misleading record this schema could hold.

    A not-assayed gene with 0 in its value field renders as "we looked and
    found none", which is a negative result the run never produced.
    """
    for state in ("not_assayed", "unsupported_by_assay", "requires_external_result",
                  "insufficient_coverage", "not_detected_above_assay_threshold"):
        with pytest.raises(S.ExtensionSchemaError, match="cannot carry a value"):
            _metric(state=state, value=0.0)
        # ...and the same state with no value is fine.
        assert _metric(state=state, value=None, limitations=()).value is None


def test_every_measured_metric_names_its_unit_and_denominator() -> None:
    with pytest.raises(S.ExtensionSchemaError, match="needs a unit"):
        _metric(unit=None)
    with pytest.raises(S.ExtensionSchemaError, match="denominator"):
        _metric(denominator=None)


def test_a_percentile_without_a_reference_is_refused() -> None:
    """AT035/AT116: no fabricated percentile while a cohort is being built."""
    with pytest.raises(S.ExtensionSchemaError, match="reference_state"):
        _metric(reference_percentile=50.0)
    with pytest.raises(S.ExtensionSchemaError, match="needs its reference ID"):
        _metric(reference_percentile=50.0, reference_state="available")
    ok = _metric(reference_percentile=50.0, reference_state="available", reference_id="ref.x/1")
    assert ok.reference_percentile == 50.0
    # A raw value stays visible when no reference exists.
    raw = _metric(reference_state="not_available")
    assert raw.value == 1.0 and raw.reference_percentile is None


def test_a_polarity_must_cite_the_context_that_supports_it() -> None:
    """New abundance readings default to descriptive, not green when high."""
    assert _metric().direction == "descriptive"
    for direction in sorted(S.DIRECTIONS_REQUIRING_CONTEXT):
        with pytest.raises(S.ExtensionSchemaError, match="must cite its context"):
            _metric(direction=direction)
        assert _metric(direction=direction, direction_context="cited study").direction == direction


def test_at022_community_distributed_genes_cannot_claim_genome_linkage() -> None:
    with pytest.raises(S.ExtensionSchemaError, match="linkage evidence"):
        _metric(kind="genetic_capacity", completeness_scope="genome_linked")
    ok = _metric(
        kind="genetic_capacity",
        completeness_scope="genome_linked",
        linkage_evidence=("contig ctg_1 carries both subunits",),
    )
    assert ok.completeness_scope == "genome_linked"
    assert _metric(kind="genetic_capacity", completeness_scope="community_assembled")


def test_a_score_must_disclose_its_definition() -> None:
    with pytest.raises(S.ExtensionSchemaError, match="must name its definition"):
        _metric(score_0_100=70.0)
    assert _metric(score_0_100=70.0, score_definition_id="def.x/1").score_0_100 == 70.0


def test_a_measured_metric_says_what_it_does_not_measure() -> None:
    with pytest.raises(S.ExtensionSchemaError, match="does not measure"):
        _metric(limitations=())


def test_an_imported_lab_value_needs_its_source_record() -> None:
    """A DNA proxy may not occupy a measured-analyte field (AT061)."""
    with pytest.raises(S.ExtensionSchemaError, match="source record"):
        _metric(kind="external_lab_measurement", unit="mmol/kg", denominator=None)
    ok = _metric(
        kind="external_lab_measurement", unit="mmol/kg", denominator=None,
        extra={"lab_record_id": "lab.scfa.butyrate.2026-09-01"},
    )
    assert ok.extra["lab_record_id"]


def test_co_detection_is_not_carriage() -> None:
    """Section 3.3: every contributing taxon states its assignment evidence."""
    with pytest.raises(S.ExtensionSchemaError, match="unknown assignment evidence"):
        S.Contributor("Bacteroides fragilis", "because it was also there")
    unique = S.Contributor("Bacteroides fragilis", "uniquely_mapped_gene")
    assoc = S.Contributor("Bacteroides fragilis", "organism_level_association")
    ambiguous = S.Contributor("Bacteroides fragilis", "ambiguous_reference_hit")
    assert unique.is_carrier
    assert not assoc.is_carrier
    assert not ambiguous.is_carrier
    m = _metric(contributing_taxa=(unique, assoc, ambiguous))
    assert len(m.contributing_taxa) == 3
    assert len(m.carriers) == 1
    assert m.to_json()["n_carrier_supported_taxa"] == 1


# --------------------------------------------------------------------------- #
# section 4.1/4.2 — the arithmetic
# --------------------------------------------------------------------------- #


def test_at012_rpkm_is_depth_invariant_and_protein_lengths_convert() -> None:
    """Doubling depth and hits leaves the abundance unchanged."""
    a = S.fragment_rpkm(100, 1200, 10_000_000)
    b = S.fragment_rpkm(200, 1200, 20_000_000)
    assert a is not None and b is not None
    assert a == pytest.approx(b, rel=1e-12)
    assert S.fragment_rpkm(100, 0, 10) is None
    assert S.fragment_rpkm(100, 1200, 0) is None
    # A 400-aa protein is 1,203 nt with its stop codon, not 400.
    assert S.cds_equivalent_length(400) == 1203
    assert S.cds_equivalent_length(400, include_stop=False) == 1200
    assert S.cds_equivalent_length(0) == 0


def test_at019_a_missing_required_step_makes_the_bottleneck_unknown_not_zero() -> None:
    assert S.bottleneck([2.0, 3.0, 5.0]) == 2.0
    assert S.bottleneck([2.0, None, 5.0]) is None
    assert S.bottleneck([]) is None
    assert S.bottleneck([0.0, 1.0]) == 0.0


def test_at035_at036_midrank_percentile_matches_the_declared_formula() -> None:
    ref = [1.0, 2.0, 3.0, 4.0]
    # count(<3)=2, count(==3)=1 -> 100*(2+0.5)/4
    assert S.midrank_percentile(3.0, ref) == pytest.approx(62.5)
    assert S.midrank_percentile(0.0, ref) == pytest.approx(0.0)
    assert S.midrank_percentile(9.0, ref) == pytest.approx(100.0)
    # An empty reference has no percentile; it does not have percentile 0.
    assert S.midrank_percentile(3.0, []) is None
    # A sparse reference: ties take half their own mass, so a zero value in a
    # mostly-zero cohort does not rank at the bottom.
    sparse = [0.0] * 8 + [5.0, 9.0]
    assert S.midrank_percentile(0.0, sparse) == pytest.approx(40.0)


def test_at014_a_shared_gene_cannot_be_counted_twice_in_an_aggregate() -> None:
    with pytest.raises(S.ExtensionSchemaError, match="counted twice"):
        S.validate_aggregate(10.0, {"inulin": 6.0, "fos": 6.0})
    S.validate_aggregate(10.0, {"inulin": 6.0, "fos": 6.0}, shared={"gh32": 2.0})


def test_route_completeness_honours_alternatives() -> None:
    """AT016: a complete alternative route is accepted on its own."""
    dxp = (
        S.Component("pdxA", "PdxA", "measured", 1.0, "rpkm", alternative_group="b6_route"),
        S.Component("pdxJ", "PdxJ", "measured", 1.0, "rpkm", alternative_group="b6_route"),
    )
    indep = (
        S.Component("pdxS", "PdxS", "not_detected_above_assay_threshold", alternative_group="b6_route"),
        S.Component("pdxT", "PdxT", "not_detected_above_assay_threshold", alternative_group="b6_route"),
    )
    m = _metric(kind="genetic_capacity", components=dxp + indep)
    # One of the two alternative groups is satisfied, and they share a group.
    assert m.completeness_from_components() == pytest.approx(1.0)
    # An AND route with a missing subunit is partial, not complete.
    anded = (
        S.Component("pta", "Pta", "measured", 1.0, "rpkm"),
        S.Component("ackA", "AckA", "insufficient_coverage"),
    )
    assert _metric(kind="genetic_capacity", components=anded).completeness_from_components() == 0.5
    assert _metric().completeness_from_components() is None


def test_fingerprints_are_content_addressed_and_stable() -> None:
    a = S.fingerprint({"db": "dbCAN v15", "k": 31}, ["x", "y"])
    b = S.fingerprint({"k": 31, "db": "dbCAN v15"}, ["x", "y"])
    c = S.fingerprint({"db": "dbCAN v14", "k": 31}, ["x", "y"])
    assert a == b and a != c and a.startswith("sha256:")


# --------------------------------------------------------------------------- #
# section 3.4 — a gene symbol is not a biochemical identity
# --------------------------------------------------------------------------- #


def _reg() -> MK.MarkerRegistry:
    reg = MK.MarkerRegistry()
    reg.add_reaction(MK.Reaction(
        "rxn.urocanate_reductase", "Urocanate reductase", "urocanate",
        "imidazole propionate", "forward", ("F22",),
        must_not_be_confused_with=("rxn.urolithin_c_9_dehydroxylation",),
    ))
    reg.add_reaction(MK.Reaction(
        "rxn.urolithin_c_9_dehydroxylation", "Urolithin C 9-dehydroxylation",
        "urolithin C", "urolithin A", "forward", ("F19",), step_logic="AND",
        must_not_be_confused_with=("rxn.urocanate_reductase",),
    ))
    return reg


def test_at023_urda_belongs_to_urocanate_not_urolithin() -> None:
    """The report has made this exact mistake: urdA is not a urolithin marker."""
    reg = _reg()
    urda = reg.add_marker(MK.Marker(
        "mk.urdA", "urdA", "MT740292.1", "nucleotide", ("rxn.urocanate_reductase",),
        ("F22",), "Eggerthella lenta", "purified_enzyme_assay",
        "validated_against_negatives",
        diagnostic_residues=("Y373",), residue_numbering_reference="UrdA reference alignment",
        negative_homologs=("fumarate reductase FrdA",),
    ))
    assert reg.markers_for("rxn.urolithin_c_9_dehydroxylation") == ()
    assert urda in reg.markers_for("rxn.urocanate_reductase")
    rxn = reg.reactions["rxn.urocanate_reductase"]
    assert rxn.substrate == "urocanate" and rxn.product == "imidazole propionate"
    assert "rxn.urolithin_c_9_dehydroxylation" in rxn.must_not_be_confused_with
    # There is no name-based binding to reach for.
    assert not hasattr(reg, "bind_by_name")


def test_a_marker_cannot_claim_validation_it_has_not_had() -> None:
    reg = _reg()
    with pytest.raises(MK.MarkerRegistryError, match="no negative set"):
        MK.Marker("m", "x", "WP_407418999.1", "protein", ("rxn.urocanate_reductase",),
                  ("F17",), "org", "purified_enzyme_assay", "validated_against_negatives")
    candidate = reg.add_marker(MK.Marker(
        "mk.candidate", "x", "WP_407418999.1", "protein", ("rxn.urocanate_reductase",),
        ("F17",), "org", "orthology_inference", "candidate_unvalidated",
    ))
    assert not candidate.is_reaction_specific
    assert candidate in reg.unvalidated()
    assert reg.audit()["n_reaction_specific_markers"] == 0


def test_diagnostic_residues_need_reference_numbering_and_actual_coverage() -> None:
    """Section 5.6: a raw read offset 373 is not reference position 373."""
    with pytest.raises(MK.MarkerRegistryError, match="numbering reference"):
        MK.Marker("m", "x", "MT740292.1", "nucleotide", ("rxn.urocanate_reductase",),
                  ("F22",), "org", "purified_enzyme_assay", diagnostic_residues=("Y373",))
    marker = MK.Marker(
        "m", "x", "MT740292.1", "nucleotide", ("rxn.urocanate_reductase",), ("F22",),
        "org", "purified_enzyme_assay", "validated_against_negatives",
        diagnostic_residues=("Y373", "M373"),
        residue_numbering_reference="UrdA reference alignment",
        negative_homologs=("FrdA",),
    )
    assert MK.diagnostic_residues_covered(marker, [373]) == (True, ())
    resolved, missing = MK.diagnostic_residues_covered(marker, [10, 20])
    assert not resolved and missing == ("Y373", "M373")


def test_a_reaction_needs_a_named_substrate_and_product() -> None:
    with pytest.raises(MK.MarkerRegistryError, match="substrate and product"):
        MK.Reaction("r", "l", "", "product", "forward", ("F01",))
    with pytest.raises(MK.MarkerRegistryError, match="cite its source"):
        MK.Reaction("r", "l", "s", "p", "forward", ())


def test_markers_join_reactions_by_id_only() -> None:
    reg = _reg()
    with pytest.raises(MK.MarkerRegistryError, match="unknown reaction"):
        reg.add_marker(MK.Marker("m", "x", "WP_1234567.1", "protein", ("rxn.invented",),
                                 ("F01",), "org", "orthology_inference"))


def test_sequence_verification_catches_the_wrong_record() -> None:
    seq = "MKVLQA"
    marker = MK.Marker(
        "m", "x", "WP_407418999.1", "protein", ("rxn.urocanate_reductase",), ("F17",),
        "org", "purified_enzyme_assay", sequence_sha256=MK.sequence_hash(seq),
        amino_acids=len(seq),
    )
    MK.verify_sequence(marker, seq)
    MK.verify_sequence(marker, "mkv\nlqa\n")          # wrapping and case normalise
    with pytest.raises(MK.MarkerRegistryError, match="hash mismatch"):
        MK.verify_sequence(marker, "MKVLQG")
    short = MK.Marker(
        "m2", "x", "WP_407418999.1", "protein", ("rxn.urocanate_reductase",), ("F17",),
        "org", "purified_enzyme_assay", sequence_sha256=MK.sequence_hash(seq), amino_acids=99,
    )
    with pytest.raises(MK.MarkerRegistryError, match="expected 99 aa"):
        MK.verify_sequence(short, seq)


def test_and_or_step_grammar_is_defined_once() -> None:
    supported = {"pdxA": True, "pdxJ": True, "pdxS": False, "pdxT": False}
    assert MK.resolve_step("AND", supported, ["pdxA", "pdxJ"])
    assert not MK.resolve_step("AND", supported, ["pdxS", "pdxT"])
    assert MK.resolve_step("OR", supported, ["pdxA", "pdxS"])
    assert not MK.resolve_step("OR", supported, ["pdxS", "pdxT"])
    assert not MK.resolve_step("AND", supported, [])
    with pytest.raises(MK.MarkerRegistryError):
        MK.resolve_step("MAYBE", supported, ["pdxA"])


# --------------------------------------------------------------------------- #
# section 16 — the completeness manifest
# --------------------------------------------------------------------------- #


def test_the_manifest_carries_every_normative_id() -> None:
    """AT009/AT071: the spec's own index is the completeness contract."""
    reg = R.load()
    reg.require_ids(f"A{i:02d}" for i in range(1, 17))
    reg.require_ids(f"F{i:03d}" for i in range(1, 39))
    # M093 is absent from the specification's table; every other M is required.
    reg.require_ids(f"M{i:03d}" for i in range(1, 120) if i != 93)
    assert "M093" not in reg.capabilities
    assert len(reg.features) == 16


def test_an_unbuilt_feature_is_engineering_debt_not_a_scientific_limit() -> None:
    """Section 15.7: these states may not be collapsed into one flag."""
    assert "not_implemented" in R.ENGINEERING_DEBT
    assert "scientifically_unresolved" not in R.ENGINEERING_DEBT
    assert "missing_input" not in R.ENGINEERING_DEBT
    for reason in ("not_implemented", "missing_input", "access_restricted",
                   "external_assay", "unsupported_by_assay", "scientifically_unresolved"):
        assert reason in R.UNAVAILABLE_REASONS
    reg = R.load()
    coverage = reg.coverage()
    assert coverage["n_engineering_debt"] == len(reg.missing())
    assert set(coverage["engineering_debt"]) == {c.capability_id for c in reg.missing()}


def test_a_capability_cannot_claim_an_implementation_it_does_not_name() -> None:
    with pytest.raises(R.RegistryError, match="must name the legacy metric"):
        R.Capability("M005", "measurement_view", "c", "b", "legacy_surfaced", module="m")
    with pytest.raises(R.RegistryError, match="must name its metric IDs"):
        R.Capability("M096", "measurement_view", "c", "b", "extension_computed", module="m")
    with pytest.raises(R.RegistryError, match="must name its module"):
        R.Capability("M096", "measurement_view", "c", "b", "evidence_card")
    with pytest.raises(R.RegistryError, match="unknown implementation state"):
        R.Capability("M096", "measurement_view", "c", "b", "done")


def test_readiness_axes_stay_independent() -> None:
    """A URL that answered is not an installed, validated database."""
    r = R.Readiness(source_identified=True, retrievable=True)
    assert r.highest_axis == "retrievable"
    assert not r.bytes_verified and not r.application_validated
    assert R.Readiness().highest_axis is None
    assert len(R.READINESS_AXES) == 6


def test_repeated_views_share_one_measurement() -> None:
    """Section 16.1: M064/M090 and M068/M073 are not four assays."""
    reg = R.load()
    repeats = reg.repeats()
    assert repeats.get("M064") == ("M090",)
    assert repeats.get("M090") == ("M064",)
    assert repeats.get("M068") == ("M073",)
    assert repeats.get("M073") == ("M068",)
    views = reg.of_kind("measurement_view")
    assert reg.unique_measurement_count() == len(views) - 2


def test_every_bound_legacy_row_points_at_a_real_panel() -> None:
    """A legacy binding that names a panel which does not exist is a lie."""
    import json

    from openbiota.samples import results_file

    results = results_file("SAMPLE2_A02")
    if not results.is_file():
        pytest.skip("no run available in this checkout")
    panels = {p["name"] for p in json.loads(results.read_text())["panels"]}
    reg = R.load()
    for cap in reg.capabilities.values():
        for legacy in cap.legacy_metric_ids:
            if legacy.startswith("panels."):
                assert legacy.split(".", 1)[1] in panels, f"{cap.capability_id} -> {legacy}"
