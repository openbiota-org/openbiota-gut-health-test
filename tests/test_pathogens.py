"""Pathogen branch — spec v5.0 acceptance tests.

Grouped by what they protect:

* the decision kernel's status contract, run verbatim from the spec;
* eligibility, so an RNA-only target is never called absent from a DNA library;
* the k-mer and alignment engines' own invariants;
* the seed catalogs, including the rule that microsporidia are counted once;
* the report layer, where the dangerous failure is not a crash but a
  *plausible-looking* number — a candidate folded into the supported count,
  or a target that was never assessed rendered as a clean negative.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from reportlab.lib import colors

from openbiota import pdfpathogens
from openbiota.pathogens import align, eligibility, kernel, kmers
from openbiota.pathogens.schema import (
    GROUP_LABELS,
    GROUP_ORDER,
    GROUPS,
    AssayManifest,
    TargetRecord,
)

CATALOG = Path(__file__).resolve().parents[1] / "pathogens" / "catalog"


# --------------------------------------------------------------------------- #
# Decision kernel (spec §5.4)
# --------------------------------------------------------------------------- #


def test_kernel_status_contract() -> None:
    """The spec's own status table, byte-for-byte."""
    assert kernel.reference_contract_checks() >= 17


def test_kernel_self_test() -> None:
    assert kernel.self_test() >= 18


def test_not_assessed_is_never_a_negative() -> None:
    """An ineligible or unusable target must not resolve to `not_detected`."""
    policy = kernel.SupportPolicy()
    for broken in (
        kernel.NormalizedEvidence(assay_eligibility="ineligible"),
        kernel.NormalizedEvidence(assay_eligibility="unknown"),
        kernel.NormalizedEvidence("eligible", "failed", True, True),
        kernel.NormalizedEvidence("eligible", "completed", False, True),
        kernel.NormalizedEvidence("eligible", "completed", True, False),
    ):
        assert kernel.sequence_status(broken, policy) == "not_assessed"


def test_confirmed_artifact_outranks_support() -> None:
    status, note = kernel.display_precedence(
        sequence="supported_sequence", contamination="confirmed_technical_artifact"
    )
    assert status == "technical_artifact"
    assert note


def test_suspected_contamination_qualifies_rather_than_deletes() -> None:
    status, note = kernel.display_precedence(
        sequence="supported_sequence", contamination="suspected"
    )
    assert status == "supported_sequence"
    assert note and "confirmation" in note


def test_support_counts_must_be_nonnegative_integers() -> None:
    policy = kernel.SupportPolicy()
    ready = kernel.NormalizedEvidence("eligible", "completed", True, True)
    with pytest.raises(ValueError, match="nonnegative"):
        kernel.sequence_status(
            kernel.NormalizedEvidence(
                **{**vars(ready), "qualifying_fragments": -1}
            ),
            policy,
        )


# --------------------------------------------------------------------------- #
# Eligibility
# --------------------------------------------------------------------------- #


def test_eligibility_self_test() -> None:
    assert eligibility.self_test() >= 8


def _target(**kw: Any) -> TargetRecord:
    base: dict[str, Any] = {
        "schema_version": "5.0",
        "target_id": "t",
        "group": "bacteria",
        "display_name": "Test organism",
        "aliases": (),
        "ncbi_taxids": (),
        "taxonomic_resolution": "species",
        "interpretation_class": "established_enteric",
        "reference_status": "genome_supported",
        "reference_accessions": (),
        "marker_accessions": (),
        "near_neighbor_target_ids": (),
        "sequence_call_profile_id": "default",
        "validation_id": "v",
        "clinical_validation_status": "computational_research_rule",
        "clinical_source_urls": (),
        "report_template_id": "generic",
        "confirmation_options": (),
        "negative_limitation_codes": (),
        "stool_role": "intestinal_shedding",
        "allowed_nucleic_acids": ("DNA",),
    }
    base.update(kw)
    return TargetRecord(**base)


def _assay(**kw: Any) -> AssayManifest:
    base: dict[str, Any] = {
        "schema_version": "5.0",
        "sample_id": "S1",
        "specimen_type": "stool",
        "collected_at": "2026-01-01T00:00:00Z",
        "participant_age_years": None,
        "nucleic_acid_protocol": "DNA",
        "reverse_transcription": False,
        "library_selection": "shotgun",
        "extraction_protocol_id": "p",
        "input_stool_mass_mg": None,
        "preservative": None,
        "mechanical_lysis": True,
        "wetlab_batch_id": None,
        "library_batch_id": None,
        "sequencing_run_id": "r",
        "read_layout": "paired",
        "sequencer": "illumina",
        "read_length_summary": None,
        "control_sample_ids": (),
        "spike_in_id": None,
        "recent_antimicrobials": None,
        "recent_probiotics": None,
        "recent_live_oral_vaccine": None,
        "symptoms": (),
        "immunocompromise": None,
        "raw_fastq_sha256": "0" * 64,
        "host_removal_bundle_id": "h",
    }
    base.update(kw)
    return AssayManifest(**base)


def test_rna_virus_is_not_assessed_by_a_dna_library() -> None:
    """The whole point of §4.1: DNA cannot answer an RNA genome."""
    rna = _target(
        target_id="rna_viruses.norovirus",
        group="rna_viruses",
        allowed_nucleic_acids=("RNA",),
        host_category="human",
        host_evidence_level="established",
        molecule_type="RNA",
        rna_profile_required=True,
    )
    record = eligibility.check_assay_eligibility(_assay(), rna)
    assert record.eligibility == "ineligible"


def test_eligibility_ignores_age() -> None:
    """Age gates interpretation, never whether the search is capable."""
    target = _target()
    child = eligibility.check_assay_eligibility(_assay(participant_age_years=3), target)
    adult = eligibility.check_assay_eligibility(_assay(participant_age_years=48), target)
    unknown = eligibility.check_assay_eligibility(_assay(), target)
    assert child.eligibility == adult.eligibility == unknown.eligibility == "eligible"


# --------------------------------------------------------------------------- #
# Engines
# --------------------------------------------------------------------------- #


def test_kmer_self_test() -> None:
    assert kmers.self_test() >= 5


def test_align_self_test() -> None:
    assert align.self_test() >= 5


def test_determinants_self_test() -> None:
    from openbiota.pathogens import determinants

    assert determinants.self_test() >= 15


def test_gene_presence_never_names_a_carrier_by_association() -> None:
    """The failure this whole module exists to prevent.

    `blaKPC` present plus *Klebsiella* present must not become
    "carbapenem-resistant Klebsiella". Only physical linkage may name a
    carrier.
    """
    from openbiota.pathogens import determinants as D

    seed = D.DeterminantSeed(
        determinant_id="amr.blakpc",
        source_section="8.3",
        display_name="blaKPC",
        gene_family="blaKPC",
        kind="amr",
        candidate_sources=("NDARO",),
        note="",
        associated_seed_ids=("bacteria.klebsiella_pneumoniae",),
    )
    mapping = D.NdaroFamilyMapping("amr.blakpc", ("blaKPC",), "exact")
    strong = align.TargetEvidence(
        target_id="blaKPC",
        qualifying_fragments=40,
        informative_bases_covered=900,
        median_identity=0.995,
    )
    r = D.adjudicate_determinant(
        seed, sample_id="S", mapping=mapping, evidence=strong, reference_bases=1000
    )
    assert r.determinant_status in D.PRESENT_STATUSES
    assert r.host_linkage == "unlinked"
    assert r.linked_target_ids == ()
    assert "not established" in r.plain_statement
    assert r.clinical_phenotype_inference == "not_inferred"
    # A present call must carry the evidence that justified it, so a reader
    # can see the coverage rather than trusting the label.
    assert r.covered_reference_fraction is not None
    assert r.aligned_identity is not None


def test_overlapping_alignments_do_not_inflate_coverage() -> None:
    """Twenty reads stacked on one window cover that window, not 20x it."""
    from openbiota.pathogens.determinants import _merged_length

    assert _merged_length([(0, 150)] * 20) == 150
    assert _merged_length([(0, 100), (50, 200)]) == 200
    assert _merged_length([(0, 100), (300, 400)]) == 200
    assert _merged_length([]) == 0


def test_installed_ndaro_map_loads_if_present() -> None:
    from openbiota.pathogens.determinants import load_ndaro_map

    path = CATALOG / "ndaro_map.yaml"
    if not path.is_file():
        pytest.skip("ndaro_map.yaml not installed")
    ndaro = load_ndaro_map(path)
    assert ndaro.by_determinant
    for det_id, mapping in ndaro.by_determinant.items():
        assert mapping.match_mode in {"exact", "prefix", "absent"}
        # An absent mapping must be a declared gap, never an empty success.
        if mapping.match_mode == "absent":
            assert mapping.gap_reason, det_id
            assert not mapping.available
        else:
            assert mapping.gene_families or mapping.family_prefixes, det_id
            assert mapping.available, det_id


def test_every_determinant_seed_has_a_mapping_decision(catalog: Any) -> None:
    """A determinant with no mapping row at all would silently vanish."""
    from openbiota.pathogens.determinants import load_ndaro_map

    path = CATALOG / "ndaro_map.yaml"
    if not path.is_file() or not catalog.determinants:
        pytest.skip("mapping or determinants not installed")
    ndaro = load_ndaro_map(path)
    missing = [
        d.determinant_id
        for d in catalog.determinants
        if ndaro.for_determinant(d.determinant_id) is None
    ]
    assert not missing, f"determinants with no mapping row: {missing[:10]}"


def test_unmapped_determinants_are_not_assessed_not_absent(catalog: Any) -> None:
    """Every seed yields a record, and gaps read as gaps."""
    from openbiota.pathogens.determinants import screen_determinants

    if not catalog.determinants:
        pytest.skip("no determinants installed")
    rows = screen_determinants(catalog.determinants, sample_id="S", ndaro=None)
    assert len(rows) == len(catalog.determinants)
    assert {r.determinant_status for r in rows} == {"not_assessed"}
    for r in rows:
        assert "not evidence that" in r.plain_statement


def test_kmer_selection_is_strand_agnostic() -> None:
    """A read and its reverse complement must nominate the same targets."""
    import numpy as np

    seq = b"ACGTTGCAAGGCTTACGGATCCGTTAACGGTTCAGGATTCCAGGTTACGA" * 3
    table = bytes.maketrans(b"ACGT", b"TGCA")
    rc = seq.translate(table)[::-1]
    forward, _ = kmers.select_kmers(seq)
    reverse, _ = kmers.select_kmers(rc)
    assert np.array_equal(np.sort(forward), np.sort(reverse))


# --------------------------------------------------------------------------- #
# Catalogs
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def catalog() -> Any:
    from openbiota.pathogens.catalog import load_catalog

    if not CATALOG.is_dir():
        pytest.skip("seed catalog not installed")
    return load_catalog(CATALOG)


def test_catalog_loads_and_every_group_is_known(catalog: Any) -> None:
    assert catalog.seeds
    for seed in catalog.seeds:
        assert seed.group in GROUPS, f"{seed.seed_id}: unknown group {seed.group}"


def test_seed_ids_are_unique(catalog: Any) -> None:
    assert len(catalog.seeds) == len({s.seed_id for s in catalog.seeds})


def test_microsporidia_are_counted_once(catalog: Any) -> None:
    """Spec P082: microsporidia sit under parasites, never also under fungi."""
    micro = [s for s in catalog.seeds if s.group == "microsporidia"]
    assert micro, "no microsporidia seeds"
    names = {s.requested_name.lower() for s in micro}
    fungal = {
        s.requested_name.lower()
        for s in catalog.seeds
        if s.group == "fungi"
    }
    assert not (names & fungal)


def test_parent_seed_references_resolve(catalog: Any) -> None:
    for seed in catalog.seeds:
        parent = getattr(seed, "parent_seed_id", None)
        if parent:
            assert parent in catalog.by_id, f"{seed.seed_id}: dangling parent {parent}"


def test_every_display_group_is_covered_by_the_group_vocabulary() -> None:
    """The report's six headings must partition the schema's groups."""
    covered: list[str] = []
    for _label, groups, _icon in pdfpathogens.DISPLAY_GROUPS:
        covered.extend(groups)
    assert len(covered) == len(set(covered)), "a group is shown under two headings"
    organism_groups = set(GROUP_ORDER) - {"virulence", "amr"}
    assert set(covered) == organism_groups
    for group in covered:
        assert group in GROUP_LABELS


# --------------------------------------------------------------------------- #
# Report layer
# --------------------------------------------------------------------------- #


def _result(
    target_id: str,
    group: str,
    status: str,
    *,
    name: str | None = None,
    frags: int = 0,
    regions: int = 0,
) -> dict[str, Any]:
    return {
        "target_id": target_id,
        "group": group,
        "display_name": name or target_id.split(".")[-1].replace("_", " ").title(),
        "display_status": status,
        "sequence_status": status,
        "resolution": "species" if status != "not_assessed" else "unresolved",
        "interpretation_class": "established_enteric",
        "reference_status": "genome_supported",
        "clinical_interpretation": f"Test statement for {target_id}.",
        "unique_supporting_fragments": frags,
        "informative_regions_supported": regions,
        "informative_bases_covered": frags * 90,
        "confirmation_options": ("stool_culture",),
    }


@pytest.fixture
def payload() -> dict[str, Any]:
    """A mixed payload: every status, in several groups."""
    return {
        "results": [
            _result("bacteria.salmonella", "bacteria", "supported_sequence",
                    name="Salmonella enterica", frags=42, regions=9),
            _result("bacteria.cdiff", "bacteria", "candidate_signal",
                    name="Clostridioides difficile", frags=3, regions=1),
            _result("bacteria.ecoli", "bacteria", "not_detected",
                    name="Escherichia coli O157"),
            _result("protozoa.giardia", "protozoa", "marker_signal",
                    name="Giardia duodenalis", frags=11, regions=4),
            _result("helminths.ascaris", "helminths", "ambiguous_signal",
                    name="Ascaris lumbricoides", frags=6),
            _result("microsporidia.ebieneusi", "microsporidia", "not_detected",
                    name="Enterocytozoon bieneusi"),
            _result("fungi.candida_auris", "fungi", "not_assessed",
                    name="Candida auris"),
            _result("rna_viruses.norovirus", "rna_viruses", "not_assessed",
                    name="Norovirus GII"),
        ],
        "determinants": [
            {
                "determinant_id": "amr.blakpc",
                "display_name": "blaKPC",
                "gene_family": "blaKPC",
                "kind": "amr",
                "amr_class": "carbapenemase",
                "determinant_status": "supported_intact_sequence",
                "host_linkage": "unlinked",
                "locus_completeness": "complete",
                "supporting_fragments": 17,
            },
            {
                "determinant_id": "toxin.stx2",
                "display_name": "stx2",
                "gene_family": "stx2",
                "kind": "toxin",
                "determinant_status": "not_detected",
                "host_linkage": "unlinked",
                "supporting_fragments": 0,
            },
        ],
        "coverage": {
            "total_targets": 8,
            "assessed": 6,
            "not_assessed": 2,
            "not_assessed_assay": 1,
            "not_assessed_reference": 1,
            "not_assessed_route": 0,
            "determinants_assessed": 2,
            "determinants_not_assessed": 0,
            "reference_gap_reasons": {"no_usable_reference": 1},
        },
    }


def test_summarise_counts_supported_only(payload: dict[str, Any]) -> None:
    """The headline number is supported findings. Candidates stay separate."""
    s = pdfpathogens.summarise(payload)
    assert s["n_supported"] == 2  # supported_sequence + marker_signal
    assert s["n_candidates"] == 1
    assert s["n_ambiguous"] == 1
    # The dangerous bug: blending them.
    assert s["n_supported"] != s["n_supported"] + s["n_candidates"]


def test_marker_signal_counts_as_supported(payload: dict[str, Any]) -> None:
    s = pdfpathogens.summarise(payload)
    names = {pdfpathogens._name(r) for r in s["supported"]}
    assert "Giardia duodenalis" in names


def test_not_assessed_is_excluded_from_every_finding_bucket(
    payload: dict[str, Any]
) -> None:
    s = pdfpathogens.summarise(payload)
    buckets = s["supported"] + s["candidates"] + s["ambiguous"]
    ids = {r["target_id"] for r in buckets}
    assert "fungi.candida_auris" not in ids
    assert "rna_viruses.norovirus" not in ids


def test_not_assessed_does_not_inflate_screened_counts(
    payload: dict[str, Any]
) -> None:
    """A group whose only target was unassessable must not read as searched."""
    s = pdfpathogens.summarise(payload)
    fungal = s["per_group"]["Fungal pathogens"]
    assert fungal["total"] == 1
    assert fungal["screened"] == 0

    viral = s["per_group"]["Viral pathogens"]
    assert viral["screened"] == 0

    bacterial = s["per_group"]["Bacterial pathogens"]
    assert bacterial["screened"] == 3


def test_every_display_group_appears_even_when_empty(
    payload: dict[str, Any]
) -> None:
    s = pdfpathogens.summarise(payload)
    assert len(s["per_group"]) == 6
    for label, _groups, _icon in pdfpathogens.DISPLAY_GROUPS:
        assert label in s["per_group"]


def test_summarise_on_empty_input_reports_that_it_did_not_run() -> None:
    for empty in (None, {}, {"results": []}):
        s = pdfpathogens.summarise(empty)
        assert s["n_supported"] == 0
        assert not s["ran"]


def test_tone_reserves_red_for_supported_findings() -> None:
    supported = _result("x.y", "bacteria", "supported_sequence", frags=9)
    candidate = _result("x.z", "bacteria", "candidate_signal", frags=2)
    absent = _result("x.w", "bacteria", "not_detected")
    unassessed = _result("x.v", "bacteria", "not_assessed")

    assert pdfpathogens._tone(supported)[0] == pdfpathogens.ALERT
    assert pdfpathogens._tone(candidate)[0] == pdfpathogens.WATCH
    # Neither an absence nor a gap may borrow the alert colour.
    assert pdfpathogens._tone(absent)[0] != pdfpathogens.ALERT
    assert pdfpathogens._tone(unassessed)[0] != pdfpathogens.ALERT


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #


def _render(flow: list[Any], path: Path) -> int:
    """Build a real PDF from a story and return its byte size.

    The fragments mention other sections ("Section 7 gives the amount…"),
    and those mentions are links; the stubs give them somewhere to land so
    the fragment can be saved on its own.
    """
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from tests.conftest import section_stubs

    doc = SimpleDocTemplate(str(path), pagesize=A4)
    doc.build([*flow, *section_stubs(_styles())])
    return path.stat().st_size


def test_evil_bug_draws_at_every_size(tmp_path: Path) -> None:
    from reportlab.lib.units import mm
    from reportlab.pdfgen import canvas as pcanvas

    out = tmp_path / "bug.pdf"
    c = pcanvas.Canvas(str(out))
    for size in (4 * mm, 13 * mm, 16 * mm, 60 * mm):
        for menace in (True, False):
            bug = pdfpathogens.EvilBug(size=size, menace=menace)
            bug.canv = c
            assert bug.wrap(0, 0) == (size, size)
            c.saveState()
            bug.draw()
            c.restoreState()
    c.save()
    assert out.stat().st_size > 1000


def _styles() -> dict[str, Any]:
    from openbiota.pdfreport import _styles as make

    return make()


def _flat_text(node: Any) -> str:
    """Every string in a flowable tree, however deeply nested.

    Reaching into `Table._cellvalues` at a fixed depth broke the moment the
    panel nested a row inside its text cell. Walking the tree tests what the
    reader sees rather than how it happens to be assembled.
    """
    if isinstance(node, str):
        return node
    if isinstance(node, (list, tuple)):
        return " ".join(_flat_text(x) for x in node)
    text = getattr(node, "text", None)
    if isinstance(text, str):
        return text
    cells = getattr(node, "_cellvalues", None)
    if cells is not None:
        return _flat_text(cells)
    return ""


def test_alert_badge_renders_with_findings(
    payload: dict[str, Any], tmp_path: Path
) -> None:
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, _styles(), pathogens=payload, section=7)
    assert flow
    assert _render(flow, tmp_path / "alert.pdf") > 1000


def test_alert_badge_renders_the_clear_case(tmp_path: Path) -> None:
    clean = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {"assessed": 1, "not_assessed": 0},
    }
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, _styles(), pathogens=clean, section=7)
    assert flow
    assert _render(flow, tmp_path / "clear.pdf") > 1000


def test_alert_badge_is_silent_when_the_branch_did_not_run() -> None:
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, _styles(), pathogens=None, section=7)
    assert flow == []


def test_glance_and_detail_render(payload: dict[str, Any], tmp_path: Path) -> None:
    st = _styles()
    flow: list[Any] = []
    pdfpathogens.pathogens_glance(
        flow, st, pathogens=payload, section=7, detail_section=15
    )
    pdfpathogens.pathogens_detail_pages(
        flow, st, pathogens=payload, section=15
    )
    assert _render(flow, tmp_path / "section.pdf") > 4000


def test_detail_pages_render_with_no_findings_at_all(tmp_path: Path) -> None:
    clean = {
        "results": [
            _result("bacteria.x", "bacteria", "not_detected"),
            _result("helminths.y", "helminths", "not_assessed"),
        ],
        "coverage": {
            "assessed": 1, "not_assessed": 1,
            "not_assessed_reference": 1, "not_assessed_assay": 0,
            "not_assessed_route": 0,
            "determinants_assessed": 0, "determinants_not_assessed": 0,
        },
    }
    st = _styles()
    flow: list[Any] = []
    pdfpathogens.pathogens_glance(
        flow, st, pathogens=clean, section=7, detail_section=15
    )
    pdfpathogens.pathogens_detail_pages(flow, st, pathogens=clean, section=15)
    assert _render(flow, tmp_path / "clean.pdf") > 2000


def test_engine_output_flows_into_the_report_layer(
    catalog: Any, tmp_path: Path
) -> None:
    """The engine-to-report seam, on real records rather than a fixture.

    With no bundle installed every target must come back `not_assessed`,
    every determinant must still produce a record, and the six report
    headings must account for every organism target exactly once.
    """
    from openbiota.pathogens.detect import run_pathogen_branch

    branch = run_pathogen_branch(
        sample_id="TEST", r1=None, r2=None, catalog=catalog,
        bundle_dir=None, out_dir=tmp_path,
    )
    payload = branch.to_json()

    assert payload["results"], "no target records emitted"
    assert all(r["display_status"] == "not_assessed" for r in payload["results"])
    # A missing bundle is a gap for determinants too, not an empty section.
    assert len(payload["determinants"]) == len(catalog.determinants)
    assert all(
        d["determinant_status"] == "not_assessed" for d in payload["determinants"]
    )

    summary = pdfpathogens.summarise(payload)
    assert summary["ran"]
    assert summary["n_supported"] == 0
    assert summary["n_candidates"] == 0
    assert summary["not_assessed"] == len(payload["results"])
    # Every organism target lands in exactly one of the six headings.
    assert sum(g["total"] for g in summary["per_group"].values()) == len(
        payload["results"]
    )
    assert all(g["screened"] == 0 for g in summary["per_group"].values())

    st = _styles()
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, st, pathogens=payload, section=7)
    pdfpathogens.pathogens_glance(
        flow, st, pathogens=payload, section=7, detail_section=15
    )
    pdfpathogens.pathogens_detail_pages(flow, st, pathogens=payload, section=15)
    assert _render(flow, tmp_path / "engine.pdf") > 4000


def test_report_sections_are_unique_and_contiguous() -> None:
    """Inserting the pathogen sections must not collide or leave a hole."""
    from openbiota.pdfreport import SECTIONS

    numbers = sorted(SECTIONS.values())
    # Zero is the one unnumbered section (the contents page); everything
    # else counts up from 1 without a gap.
    assert numbers.count(0) <= 1
    numbered = [n for n in numbers if n]
    assert numbered == list(range(1, len(numbered) + 1))
    assert "pathogens" in SECTIONS
    assert SECTIONS["pathogens"] < SECTIONS["pathogens_detail"]


def test_bundle_ownership_decoys_and_bins_are_competitive(
    catalog: Any, tmp_path: Path
) -> None:
    """The bundle builder's three promises, on a genome pair we control.

    Sequence two targets share is owned by neither; sequence the host or a
    food animal carries is stripped even when it is unique among targets;
    and a target's informative bins are exactly those left holding
    discriminating k-mers. A region repeated inside one genome is still
    that genome's, credited to its first occurrence, so the repeat's bin
    is empty. Building twice must give the identical index.
    """
    import hashlib

    import numpy as np

    from openbiota.pathogens.kmers import SHARED, KmerIndex, select_kmers
    from openbiota.pathogens.refs import BIN_BASES, Asset, Resolution, build_bundle

    rng = np.random.default_rng(7)

    def random_dna(n: int) -> bytes:
        return rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), n).tobytes()

    nb = BIN_BASES
    # Each segment is its own contig, sized to a whole number of bins, so no
    # k-mer straddles a junction and every segment lands on known bins.
    # A: bins 0-1 shared with B, bin 2 also in the decoy, bins 3-4 unique,
    #    bin 5 an exact repeat of bin 3.
    a_shared, a_decoy = random_dna(2 * nb), random_dna(nb)
    a_unique = random_dna(2 * nb)
    a = [a_shared, a_decoy, a_unique, a_unique[:nb]]
    # B: bins 0-3 unique, bins 4-5 the segment shared with A.
    b = [random_dna(4 * nb), a_shared]
    decoy = [random_dna(3 * nb), a_decoy]

    def fasta(name: str, contigs: list[bytes]) -> Path:
        p = tmp_path / f"{name}.fna"
        with p.open("wb") as out:
            for i, seq in enumerate(contigs):
                out.write(f">{name}_{i}\n".encode() + seq + b"\n")
        return p

    seed_a, seed_b = (s.seed_id for s in catalog.seeds[:2])
    paths = {seed_a: fasta("a", a), seed_b: fasta("b", b)}
    assets = {
        sid: Asset(
            seed_id=sid, path=p, sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            sequence_type="genome", accession_version=f"GCF_TEST_{sid}",
            source_url="test", retrieved_at="now", bases=len(p.read_bytes()),
        )
        for sid, p in paths.items()
    }
    resolutions = {
        sid: Resolution(seed_id=sid, requested_name=sid, resolution_state="resolved",
                        accepted_taxon=sid, ncbi_taxids=(1,), route="genome")
        for sid in paths
    }

    def build(out: Path) -> KmerIndex:
        manifest = build_bundle(
            catalog, resolutions, assets, out_dir=out,
            decoy_fastas=[fasta("decoy", decoy)], build_aligner_index=False,
        )
        assert manifest.entries[seed_a]["reference_status"] == "genome_supported"
        assert seed_b in manifest.entries[seed_a]["near_neighbor_target_ids"]
        assert seed_a in manifest.entries[seed_b]["near_neighbor_target_ids"]
        return KmerIndex.load(out / "kmer_index.npz", mmap=False)

    index = build(tmp_path / "bundle1")
    owner_a, owner_b = index.target_ids.index(seed_a), index.target_ids.index(seed_b)

    def owners_of(seq: bytes) -> np.ndarray:
        kmers, _ = select_kmers(seq)
        found = index.lookup(kmers)
        return found[found != -1]

    # 1. Sequence both genomes carry belongs to neither.
    shared = owners_of(a_shared)
    assert shared.size and np.all(shared == SHARED)

    # 2. Sequence the decoy carries is stripped, though only A had it.
    decoyed = owners_of(a_decoy)
    assert decoyed.size and np.all(decoyed == SHARED)

    # 3. What is left over is the target's, and only that target's.
    assert np.all(owners_of(a_unique) == owner_a)
    assert np.all(owners_of(b[0]) == owner_b)

    # Informative bins follow ownership: A keeps its two unique bins and
    # loses the shared, decoyed and repeated ones.
    assert index.is_informative_bin(owner_a, 3)
    assert index.is_informative_bin(owner_a, 4)
    for empty in (0, 1, 2, 5):
        assert not index.is_informative_bin(owner_a, empty), f"bin {empty}"
    assert not any(index.is_informative_bin(owner_b, i) for i in (4, 5))

    # A rebuild must agree exactly — ownership cannot depend on sort order.
    again = build(tmp_path / "bundle2")
    assert np.array_equal(index.kmers, again.kmers)
    assert np.array_equal(index.owners, again.owners)
    assert index.target_ids == again.target_ids
    for owner in (owner_a, owner_b):
        assert np.array_equal(
            index.informative_bins[owner], again.informative_bins[owner]
        )


def _finding(
    target_id: str,
    group: str,
    klass: str,
    *,
    status: str = "supported_sequence",
    fragments: int = 200,
    fpm: float | None = 25.0,
    stool_role: str = "intestinal_shedding",
    name: str | None = None,
    qualifier: str | None = None,
) -> dict[str, Any]:
    return {
        "target_id": target_id, "group": group, "interpretation_class": klass,
        "stool_role": stool_role,
        "display_name": name or target_id.split(".")[-1].replace("_", " ").title(),
        "display_status": status, "sequence_status": status,
        "unique_supporting_fragments": fragments,
        "informative_regions_supported": 5,
        "normalized_fragments_per_million": fpm,
        "resolution": "species", "display_qualifier": qualifier,
    }


def test_a_normal_gut_or_food_organism_is_never_called_a_pathogen() -> None:
    """The bug that made this section untrustworthy.

    Brewer's yeast and an ordinary gut anaerobe were printed in the same red
    "pathogens found" list as Shigella, and counted in the headline. The
    catalogue already says what each organism is; the report must use it.
    """
    payload = {
        "results": [
            # Brewer's and baker's yeast. Arrives with food. Not a pathogen.
            _finding("fungi.saccharomyces_cerevisiae", "fungi", "background_or_decoy",
                     fragments=9000, fpm=1200.0,
                     stool_role="environmental_or_dietary",
                     name="Saccharomyces cerevisiae"),
            # An abundant, entirely ordinary gut anaerobe.
            _finding("bacteria.mediterraneibacter_gnavus", "bacteria",
                     "background_or_decoy", fragments=30150, fpm=4249.0,
                     stool_role="carriage_common", name="Mediterraneibacter gnavus"),
            # A real one.
            _finding("bacteria.shigella_sonnei", "bacteria", "established_enteric",
                     name="Shigella sonnei"),
        ],
        "coverage": {"assessed": 3, "not_assessed": 0},
    }
    summary = pdfpathogens.summarise(payload)

    assert summary["n_pathogens"] == 1
    assert [r["display_name"] for r in summary["pathogens"]] == ["Shigella sonnei"]
    normal = summary["by_tier"]["normal"]
    assert {r["display_name"] for r in normal["all"]} == {
        "Saccharomyces cerevisiae", "Mediterraneibacter gnavus"
    }
    # Being the most abundant thing in the sample must not make it red.
    for record in normal["all"]:
        colour, _bg = pdfpathogens._tone(record)
        assert colour is not pdfpathogens.ALERT

    # ... and it is still measured and reported, not hidden.
    yeast = normal["all"][0] if normal["all"][0]["display_name"].startswith("S") \
        else normal["all"][1]
    assert pdfpathogens.amount_of(yeast)["band"] == "large signal"


def test_the_headline_counts_pathogens_not_every_named_organism() -> None:
    """Opportunists and unsettled commensals are named, never counted as found."""
    payload = {
        "results": [
            _finding("bacteria.enterococcus_faecalis", "bacteria",
                     "conditional_opportunist", stool_role="carriage_common"),
            _finding("bacteria.parvimonas_micra", "bacteria",
                     "uncertain_enteric_role", stool_role="carriage_common"),
            _finding("protozoa.toxoplasma_gondii", "protozoa",
                     "extraintestinal_watch", stool_role="tissue_restricted"),
        ],
        "coverage": {"assessed": 3, "not_assessed": 0},
    }
    summary = pdfpathogens.summarise(payload)
    assert summary["n_supported"] == 3, "all three are still named findings"
    assert summary["n_pathogens"] == 0, "none of the three is a gut pathogen"
    assert summary["by_tier"]["opportunist"]["n"] == 1
    assert summary["by_tier"]["uncertain"]["n"] == 1
    assert summary["by_tier"]["elsewhere"]["n"] == 1

    st = _styles()
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, st, pathogens=payload, section=7)
    text = _flat_text(flow)
    assert "No disease-causing organisms found" in text


def test_an_unknown_class_is_not_promoted_to_pathogen() -> None:
    """If the catalogue has not said it causes disease, neither may we."""
    payload = {
        "results": [_finding("bacteria.mystery", "bacteria", "uncertain_enteric_role")],
        "coverage": {},
    }
    # A class the tier map has never seen.
    payload["results"][0]["interpretation_class"] = "some_future_class"
    assert pdfpathogens.tier_of(payload["results"][0]) == "uncertain"
    assert pdfpathogens.summarise(payload)["n_pathogens"] == 0


def test_a_trace_with_no_qualifying_fragment_is_not_a_finding() -> None:
    """Zero fragments is no evidence, not weak evidence.

    104 of 145 trace rows across five real samples had nothing behind them,
    and that is how Coccidioides and hookworms reached healthy reports.
    """
    payload = {
        "results": [
            _finding("fungi.coccidioides_posadasii", "fungi", "extraintestinal_watch",
                     status="candidate_signal", fragments=0, fpm=0.0),
            _finding("bacteria.shigella_boydii", "bacteria", "established_enteric",
                     status="candidate_signal", fragments=2, fpm=0.28),
        ],
        "coverage": {},
    }
    summary = pdfpathogens.summarise(payload)
    assert summary["n_candidates"] == 1
    assert summary["n_bare_nominations"] == 1
    named = [r["display_name"] for t in summary["by_tier"].values() for r in t["all"]]
    assert not any("Coccidioides" in n for n in named)


def test_the_amount_is_reported_and_carries_no_clinical_threshold() -> None:
    """"How much" was the missing answer; it must not become a false one."""
    big = pdfpathogens.amount_of(
        _finding("x.y", "bacteria", "established_enteric", fragments=30150, fpm=4249.0)
    )
    small = pdfpathogens.amount_of(
        _finding("x.z", "bacteria", "established_enteric", fragments=5, fpm=0.7)
    )
    none = pdfpathogens.amount_of(
        _finding("x.w", "bacteria", "established_enteric", fragments=0, fpm=None)
    )
    assert big["band"] == "large signal" and small["band"] == "at the detection floor"
    assert big["percent"] == pytest.approx(0.42489, abs=1e-4)
    assert "0.42%" in big["share"]
    assert none["band"] is None and none["figure"] == "—"
    # The bands describe signal size only. No band may imply a diagnosis.
    for _floor, band, words in pdfpathogens._AMOUNT_BANDS:
        for text in (band, words):
            assert not any(
                w in text.lower()
                for w in ("infect", "diagnos", "danger", "toxic", "moderate", "severe")
            ), text


def test_carriage_is_stated_so_normal_can_be_told_from_abnormal() -> None:
    common = pdfpathogens.carriage_words(
        _finding("x.y", "bacteria", "conditional_opportunist",
                 stool_role="carriage_common")
    )
    dietary = pdfpathogens.carriage_words(
        _finding("x.z", "fungi", "background_or_decoy",
                 stool_role="environmental_or_dietary")
    )
    pathogen = pdfpathogens.carriage_words(
        _finding("x.w", "bacteria", "established_enteric",
                 stool_role="intestinal_shedding")
    )
    assert "healthy people" in common
    assert "food" in dietary
    assert "not a normal resident" in pathogen


def test_the_controls_caveat_is_not_printed_on_the_rows() -> None:
    """A per-sample qualifier is not a per-row label.

    The run-wide fine print that used to carry it has gone from the
    overview altogether - it did not help a reader decide anything - so
    the only thing left to check is that the qualifier never leaks into
    the finding rows.
    """
    note = "no extraction or library controls for this sample"
    payload = {
        "results": [
            _finding(f"bacteria.b{i}", "bacteria", "established_enteric",
                     qualifier=note)
            for i in range(6)
        ],
        "coverage": {},
    }
    st = _styles()
    row = pdfpathogens._finding_row(st, payload["results"][0])
    assert not any(note in getattr(cell, "text", "") for cell in row)
    assert not hasattr(pdfpathogens, "_run_caveats"), "the run-wide fine print is back"


def test_the_installed_advice_library_is_valid_and_honest(catalog: Any) -> None:
    """The content that turns a name into something a reader can use.

    Three properties matter more than coverage. Every entry is keyed to a
    real target, or it would never be shown. Every agent states how strong
    its evidence is, because a culture result and a randomised trial are
    not the same claim and supplement content is exactly where that
    distinction gets lost. And an organism with no credible option says so
    outright rather than being padded with a plausible-sounding protocol.
    """
    from openbiota.pathogens.agents import EVIDENCE_TIERS, load_agent_library

    ids = [s.seed_id for s in catalog.seeds]
    library = load_agent_library(CATALOG, known_seed_ids=ids)
    if not library:
        pytest.skip("no agent files installed")

    for seed_id, advice in library.by_seed.items():
        assert seed_id in set(ids)
        # The four reader-facing questions are answered for every organism.
        for field in (
            advice.what_it_is, advice.why_it_matters,
            advice.normal_carriage, advice.overgrowth_signal,
        ):
            assert len(field) > 40, f"{seed_id}: stub text"
        if not advice.has_agents:
            assert len(advice.no_agents_reason) > 40, (
                f"{seed_id}: an empty list must say why"
            )
        for agent in advice.agents:
            assert agent.evidence_tier in EVIDENCE_TIERS
            assert agent.citation_url.startswith("http")
            # A laboratory result must not be worded as a human effect.
            if agent.evidence_tier == "in_vitro":
                assert "no human evidence" in agent.tier_words
        # Strongest evidence first, always.
        ranks = [a.tier_rank for a in advice.agents]
        assert ranks == sorted(ranks), f"{seed_id}: agents out of tier order"

    # Some organisms genuinely have nothing worth taking, and the library
    # must be willing to say so — that is the honest answer for an ordinary
    # resident, and its absence would mean the content had been padded.
    assert any(not a.has_agents for a in library.by_seed.values())
    # ... and some must have real human evidence, or the library is useless.
    assert any(a.human_backed for a in library.by_seed.values())


def test_a_detail_card_longer_than_a_page_still_renders(tmp_path: Path) -> None:
    """A card must be able to run onto the next page.

    ReportLab cannot split a single table cell, so once these cards carried
    what the organism is, whether carrying it is normal and what acts
    against it, a long one exceeded the frame and raised LayoutError. The
    whole report was then written without a PDF while the batch still
    reported success against a stale file.
    """
    from openbiota.pathogens.agents import AgentAdvice, OrganismAdvice

    lots = " ".join(["Sentence about this organism and its behaviour."] * 40)
    advice = OrganismAdvice(
        seed_id="bacteria.x", display_name="Test organism",
        what_it_is=lots, why_it_matters=lots,
        normal_carriage=lots, overgrowth_signal=lots,
        agents=tuple(
            AgentAdvice(
                name=f"Agent {i}", kind="botanical", evidence_tier="in_vitro",
                effect=lots, dose_studied="not established", reaches_gut=lots,
                citation_title="A paper", citation_url="https://example.org/a",
                caution=lots,
            )
            for i in range(6)
        ),
    )
    record = _finding("bacteria.x", "bacteria", "established_enteric")
    st = _styles()
    flow = pdfpathogens._detail_card(st, record, advice)
    # The card links back to its at-a-glance row; alone, give it one to land on.
    from openbiota.pdflinks import Anchor, pathogen_dest

    flow.append(Anchor(pathogen_dest(str(record["target_id"]))))
    assert _render(flow, tmp_path / "long.pdf") > 3000


def test_a_decisively_beaten_relative_is_not_a_candidate() -> None:
    """Losing the competition is an answer, not a weaker kind of finding.

    Conserved and repetitive sequence aligns to half the tree. If every
    reference a fragment merely touched were credited, a healthy stool
    library would report systemic fungi and hookworms as candidates — the
    competitive stage overruled by its own inputs. Only the winner and
    anything inside the ambiguity margin may be credited.
    """
    from openbiota.pathogens.align import _rec, summarise_alignments
    from openbiota.pathogens.catalog import RESEARCH_DNA_V1 as profile

    margin = profile.min_best_minus_other_taxon_alignment_score

    def rec(query: str, target: str, pos: int, score: int) -> Any:
        return _rec(query, target, pos, score=score)

    beaten = [
        r
        for i in range(6)
        for r in (
            rec(f"f{i}", "winner", i * 10_000, 300),
            rec(f"f{i}", "loser", i * 10_000, 300 - margin - 50),
        )
    ]
    out = summarise_alignments(beaten, profile=profile)
    assert "loser" not in out.by_target, "a decisively beaten target is not a candidate"
    assert out.by_target["winner"].any_candidate_signal

    # Inside the margin the two are genuinely indistinguishable, and both
    # must survive — as ambiguous, never silently resolved onto the winner.
    tied = [
        r
        for i in range(6)
        for r in (
            rec(f"g{i}", "winner", i * 10_000, 300),
            rec(f"g{i}", "close", i * 10_000, 300 - max(0, margin - 1)),
        )
    ]
    out = summarise_alignments(tied, profile=profile)
    for name in ("winner", "close"):
        assert out.by_target[name].any_candidate_signal
        assert out.by_target[name].only_unresolved_taxonomic_support
        assert out.by_target[name].qualifying_fragments == 0


def test_kmer_nomination_alone_never_reaches_the_report() -> None:
    """A nomination competitive alignment refuted is a negative, not a hit.

    K-mer classification is candidate *generation*; the competitive stage is
    what adjudicates it (spec §2, Bradford et al. 2024). If a nomination
    could re-enter the kernel as `any_candidate_signal`, the stage that
    exists to check it would be overruled by the stage it checks — measured
    on a real stool library that was ~70 spurious candidates per sample,
    Blastomyces and hookworms among them. The nomination must survive only
    as an audit trail.
    """
    import ast
    import inspect
    from dataclasses import replace

    from openbiota.pathogens import detect as detect_mod

    src = inspect.getsource(detect_mod._result_for)
    assert "nominated_but_not_confirmed_competitively" in src, (
        "the audit trail for a refuted nomination must survive"
    )

    # The kernel call must take its candidate flag from alignment evidence
    # alone: no `or nominated`, however spelled.
    tree = ast.parse(inspect.getsource(detect_mod))
    flags = [
        kw.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "NormalizedEvidence"
        for kw in node.keywords
        if kw.arg == "any_candidate_signal"
    ]
    assert flags, "no NormalizedEvidence(any_candidate_signal=...) call found"
    for flag in flags:
        names = {
            n.id for n in ast.walk(flag) if isinstance(n, ast.Name)
        } | {
            n.attr for n in ast.walk(flag) if isinstance(n, ast.Attribute)
        }
        assert "nominated" not in names, (
            "k-mer nomination must not feed the kernel's candidate flag"
        )
        assert "any_candidate_signal" in names

    # And the kernel itself must still call an unaligned target a negative.
    from openbiota.pathogens.kernel import (
        NormalizedEvidence as NE,
    )
    from openbiota.pathogens.kernel import (
        SupportPolicy,
        sequence_status,
    )

    ready = NE("eligible", "completed", True, True)
    assert sequence_status(ready, SupportPolicy()) == "not_detected"
    assert sequence_status(
        replace(ready, any_candidate_signal=True), SupportPolicy()
    ) == "candidate_signal"


def test_competitive_index_is_shared_and_built_once(tmp_path: Path) -> None:
    """The index every sample competes in: one member set, one build.

    Two samples nominating different targets must land on the same cache
    entry, or each pays a multi-hour rebuild; a second call must not shell
    out again; and two processes arriving together must produce one build,
    not a race into the same directory. A build that dies halfway must
    leave nothing a later run would mistake for a finished index.
    """
    import shutil
    import subprocess

    from openbiota.pathogens import refs as refs_mod
    from openbiota.pathogens.refs import (
        BundleManifest,
        build_competitive_index,
        competitive_members,
    )

    usable = {"a": "genome_supported", "b": "marker_only", "c": "genome_supported"}
    manifest = BundleManifest(
        bundle_id="test", catalog_version="v1",
        entries={
            **{
                t: {"reference_status": s, "near_neighbor_target_ids": []}
                for t, s in usable.items()
            },
            "gap": {"reference_status": "unresolved_taxonomy",
                    "near_neighbor_target_ids": ["a"]},
        },
        lock={}, decoys=(),
    )
    seq_dir = tmp_path / "sequences"
    seq_dir.mkdir()
    for target in (*usable, "gap"):
        (seq_dir / f"{target}.fna").write_text(f">{target}\nACGT\n")

    # Only installed sequence competes: a target with no usable reference
    # contributes its neighbours but never itself.
    assert competitive_members(manifest) == ("a", "b", "c")
    assert competitive_members(manifest, ["gap"]) == ("a",)
    assert set(competitive_members(manifest, ["a", "b"])) <= set(
        competitive_members(manifest)
    )

    calls: list[list[str]] = []

    def fake_build(argv: list[str], **_kw: Any) -> Any:
        calls.append(argv)
        Path(argv[-1] + ".1.bt2").write_text("index")
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkey = pytest.MonkeyPatch()
    monkey.setattr(refs_mod.shutil, "which", lambda _n: "/usr/bin/bowtie2-build")
    monkey.setattr(refs_mod.subprocess, "run", fake_build)
    try:
        first = build_competitive_index(tmp_path, None, manifest=manifest)
        assert first is not None and len(calls) == 1
        # The second sample reuses it rather than shelling out again.
        second = build_competitive_index(tmp_path, None, manifest=manifest)
        assert second == first, "cache key must not vary between samples"
        assert len(calls) == 1, "a cached index must not be rebuilt"
        # Narrowing to one sample's candidates is still possible, but it
        # keys on that sample — which is exactly why detection does not use
        # it, and why the guard below pins the shared call.
        narrowed = build_competitive_index(tmp_path, ["a"], manifest=manifest)
        assert narrowed != first
        import ast

        tree = ast.parse((Path(refs_mod.__file__).parent / "detect.py").read_text())
        sites = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and getattr(node.func, "id", None) == "build_competitive_index"
        ]
        assert len(sites) == 1
        candidates_arg = sites[0].args[1]
        assert isinstance(candidates_arg, ast.Constant) and candidates_arg.value is None, (
            "detection must ask for the shared whole-bundle index"
        )

        # Nothing is published unless bowtie2 succeeded.
        def failing(argv: list[str], **_kw: Any) -> Any:
            return subprocess.CompletedProcess(argv, 1, "", "bowtie2-build: boom")

        monkey.setattr(refs_mod.subprocess, "run", failing)
        shutil.rmtree(tmp_path / "competitive")
        assert build_competitive_index(tmp_path, None, manifest=manifest) is None
        leftovers = list((tmp_path / "competitive").glob("*/index.1.bt2*"))
        assert not leftovers, f"failed build left an index behind: {leftovers}"
    finally:
        monkey.undo()


def test_no_hardcoded_section_numbers_in_report_prose() -> None:
    """Every "see section N" in the report must come from SECTIONS.

    Inserting the pathogen sections renumbered everything after 6; a
    literal number in the prose is a cross-reference that silently
    points at the wrong page. Only f-string lookups are allowed.
    """
    import re

    pdf_dir = Path(pdfpathogens.__file__).parent
    literal = re.compile(r"(?i)\bsections?\s+\d+(?![\d_'\]])")
    offenders: list[str] = []
    for path in sorted(pdf_dir.glob("pdf*.py")):
        if path.name == "pdflinks.py":
            # The module that turns "section N" into links: its examples and
            # self-test are about the phrase itself, not cross-references.
            continue
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if literal.search(line):
                offenders.append(f"{path.name}:{lineno}: {line.strip()}")
    assert not offenders, "\n".join(offenders)


# --------------------------------------------------------------------------- #
# Species resolution and pathotype gating (openbiota.pathogens.resolve)
# --------------------------------------------------------------------------- #


def test_resolve_kernel_contract() -> None:
    """The pass's own contract checks, run as part of the suite."""
    from openbiota.pathogens import resolve

    assert resolve.self_test() == 0


def test_a_shigella_species_is_not_named_without_its_invasion_marker() -> None:
    """The bug that reported ordinary E. coli as three dysentery organisms.

    Shigella is the same genomic species as commensal E. coli (GTDB R07-RS207),
    so genome breadth cannot name it. Without ipaH the species row must become
    an unresolved complex signal, keep every count, and stop counting as a
    disease-causing organism.
    """
    from openbiota.pathogens import resolve

    families = {
        "bacteria.shigella_sonnei": "Shigella/EIEC",
        "bacteria.shigella_eiec_combined_group": "Shigella/EIEC",
    }
    rows = {
        tid: _finding(tid, "bacteria", "established_enteric", fragments=432, fpm=54.0,
                      name="Shigella sonnei" if "sonnei" in tid else "Shigella/EIEC group")
        for tid in families
    }
    rows = {k: dict(v) for k, v in rows.items()}
    for v in rows.values():
        v.update(
            {"sequence_status": "supported_sequence", "informative_regions_supported": 168,
             "informative_bases_covered": 50_000, "reference_breadth_fraction": 0.0103,
             "plain_statement": "base.", "reason_codes": (), "confirmation_options": (),
             "display_qualifier": None, "species_resolution": "resolved",
             "complex_id": None, "complex_label": None, "absorbed_from": (),
             "pathotype_evidence": "not_applicable", "pathotype_markers": (),
             "counts_as_pathogen": True, "report_tier": None,
             "carriage_statement": None, "carriage_source": None}
        )

    class Row(dict):
        def __getattr__(self, item):  # dict rows stand in for the frozen dataclass
            return self[item]

    class Det:
        def __init__(self, did, status):
            self.determinant_id, self.determinant_status = did, status

    import dataclasses

    @dataclasses.dataclass(frozen=True)
    class Frozen:
        target_id: str
        display_name: str
        interpretation_class: str
        sequence_status: str
        display_status: str
        resolution: str
        unique_supporting_fragments: int
        informative_regions_supported: int
        informative_bases_covered: int
        reference_breadth_fraction: float
        plain_statement: str
        reason_codes: tuple = ()
        confirmation_options: tuple = ()
        display_qualifier: str | None = None
        species_resolution: str = "resolved"
        complex_id: str | None = None
        complex_label: str | None = None
        absorbed_from: tuple = ()
        pathotype_evidence: str = "not_applicable"
        pathotype_markers: tuple = ()
        counts_as_pathogen: bool = True
        report_tier: str | None = None
        carriage_statement: str | None = None
        carriage_source: str | None = None

    frozen = {
        tid: Frozen(
            target_id=tid, display_name=v["display_name"],
            interpretation_class="established_enteric",
            sequence_status="supported_sequence", display_status="supported_sequence",
            resolution="species", unique_supporting_fragments=432,
            informative_regions_supported=168, informative_bases_covered=50_000,
            reference_breadth_fraction=0.0103, plain_statement="base.",
        )
        for tid, v in rows.items()
    }

    out = resolve.resolve(
        frozen, [Det("pathotype_marker.ipah", "not_detected")], families=families
    )
    species = out["bacteria.shigella_sonnei"]
    assert species.sequence_status == "ambiguous_signal"
    assert species.counts_as_pathogen is False
    assert species.report_tier == "unresolved_complex"
    assert species.unique_supporting_fragments == 432, "no evidence may be discarded"
    assert "ipah" in species.plain_statement.lower()
    group = out["bacteria.shigella_eiec_combined_group"]
    assert group.counts_as_pathogen is False
    assert "commensal Escherichia coli" in group.plain_statement

    with_marker = resolve.resolve(
        frozen, [Det("pathotype_marker.ipah", "supported_intact_sequence")], families=families
    )
    assert with_marker["bacteria.shigella_sonnei"].counts_as_pathogen is True


def test_the_headline_excludes_carriage_without_the_toxin_gene() -> None:
    """Toxin-negative C. difficile is not a disease-causing organism."""
    payload = {
        "results": [
            _finding("bacteria.clostridioides_difficile", "bacteria",
                     "toxin_or_pathotype_dependent", fragments=220, fpm=27.7,
                     name="Clostridioides difficile")
            | {"pathotype_evidence": "not_detected", "counts_as_pathogen": False,
               "report_tier": "pathotype_negative"},
        ],
        "coverage": {"assessed": 1, "not_assessed": 0},
    }
    summary = pdfpathogens.summarise(payload)
    assert summary["n_pathogens"] == 0
    assert summary["n_pathotype_negative"] == 1
    assert pdfpathogens.tier_of(payload["results"][0]) == "pathotype_negative"
    assert not pdfpathogens.TIERS["pathotype_negative"]["counts_as_pathogen"]


def test_unresolved_species_are_tiered_out_of_the_headline() -> None:
    payload = {
        "results": [
            _finding("bacteria.shigella_sonnei", "bacteria", "established_enteric",
                     fragments=432, fpm=54.0, name="Shigella sonnei")
            | {"report_tier": "unresolved_complex", "counts_as_pathogen": False,
               "species_resolution": "not_resolvable_within_complex"},
        ],
        "coverage": {"assessed": 1, "not_assessed": 0},
    }
    summary = pdfpathogens.summarise(payload)
    assert summary["n_pathogens"] == 0
    assert summary["n_unresolved_complex"] == 1


def test_every_amount_band_describes_a_signal_not_a_clinical_state() -> None:
    for _floor, band, words in pdfpathogens._AMOUNT_BANDS:
        assert "signal" in band or "detection floor" in band, band
        assert any(w in words for w in ("DNA", "detection", "depth")), words


def test_a_row_with_no_specific_fragment_is_never_a_named_species_finding() -> None:
    """"DNA present ... 0 specific fragments" is a self-contradiction.

    The resolution pass must not demote a zero-evidence row into the
    "species cannot be named" tier, and the renderer must not honour such a
    verdict if an older run stored one.
    """
    from openbiota.pdfpathogens import tier_of

    stale = {
        "report_tier": "unresolved_complex",
        "unique_supporting_fragments": 0,
        "informative_regions_supported": 0,
        "interpretation_class": "established_enteric_pathogen",
        "display_name": "Enteroinvasive Escherichia coli (EIEC)",
    }
    assert tier_of(stale) != "unresolved_complex"

    stale_carriage = dict(stale, report_tier="pathotype_negative")
    assert tier_of(stale_carriage) != "pathotype_negative"

    # With real organism-specific evidence the verdict is honoured as before.
    real = dict(stale, unique_supporting_fragments=42, informative_regions_supported=3)
    assert tier_of(real) == "unresolved_complex"


# --------------------------------------------------------------------------- #
# The page-1 band is a budget, not just a graphic: pdfsummary hands the
# disease-signature and metabolite lists whatever vertical space is left once
# this panel has been measured. Height, colour and fit are all load-bearing.
# --------------------------------------------------------------------------- #

#: A negative result is the common case and must stay small: this is the
#: budget the disease-signature and metabolite lists are sized against.
_BAND_CLEAR_MM = 15.5
#: A positive result earns one extra line, because naming the organism on
#: page 1 beats making someone turn the page to find out what was found.
_BAND_ALARM_MM = 19.0


def _band(pathogens: dict[str, Any]) -> list[Any]:
    flow: list[Any] = []
    pdfpathogens.pathogen_alert(flow, _styles(), pathogens=pathogens, section=7)
    return flow


def _band_height_mm(flow: list[Any]) -> float:
    from reportlab.lib.units import mm

    return sum(f.wrap(170 * mm, 400)[1] for f in flow) / mm


def test_the_page_one_band_stays_within_its_height_budget() -> None:
    """Every millimetre here is a row taken off the two lists above it."""
    # Nothing found, but plenty still open: the chips must not push it onto a
    # second line and cost the lists a row each.
    busy = {
        "results": [
            _finding("bacteria.b", "bacteria", "conditional_opportunist"),
            dict(
                _finding("bacteria.c", "bacteria", "established_enteric"),
                report_tier="unresolved_complex",
            ),
            dict(
                _finding("bacteria.d", "bacteria", "established_enteric"),
                report_tier="pathotype_negative",
            ),
        ],
        "coverage": {"assessed": 333, "not_assessed": 153},
    }
    assert pdfpathogens.summarise(busy)["n_pathogens"] == 0
    assert _band_height_mm(_band(busy)) <= _BAND_CLEAR_MM

    clean = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {"assessed": 486, "not_assessed": 0},
    }
    assert _band_height_mm(_band(clean)) <= _BAND_CLEAR_MM


def test_partial_coverage_is_signalled_beside_a_negative_verdict() -> None:
    """"None found" must not sit beside a reassuring green coverage chip.

    A reference genuinely still missing has to show, or page one reads as two
    reassurances rather than one claim and its scope. What must *not* show as
    a gap is a target that was searched and answers at a rank it shares with
    a relative, or one whose genome is RNA and cannot be in a DNA library:
    counting those made a complete screen read "332 of 486 targets searched"
    in amber, as though a third of the catalogue had been skipped.
    """
    partial = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {
            "assessed": 332, "not_assessed": 154, "total_targets": 486,
            "searched": 332, "answerable": 486, "reference_pending_no_public_sequence": 154,
        },
    }
    full = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {
            "assessed": 486, "not_assessed": 0, "total_targets": 486,
            "searched": 486, "answerable": 486,
        },
    }

    def _chips(node: Any) -> list[pdfpathogens.Chip]:
        """Every Chip anywhere in a nested flowable/table structure."""
        if isinstance(node, pdfpathogens.Chip):
            return [node]
        found: list[pdfpathogens.Chip] = []
        for attr in ("_cellvalues", "_argW"):
            if attr == "_cellvalues" and hasattr(node, attr):
                for row in node._cellvalues:
                    for cell in row:
                        found += _chips(cell)
        if isinstance(node, (list, tuple)):
            for item in node:
                found += _chips(item)
        return found

    def _coverage_chip(payload: dict) -> pdfpathogens.Chip:
        chips = [
            c for c in _chips(_band(payload)) if "targets searched" in c.text
        ]
        assert chips, "the coverage chip must always be present"
        return chips[0]

    assert _coverage_chip(partial).colour == pdfpathogens.WATCH
    assert _coverage_chip(full).colour == pdfpathogens.CLEAR
    assert _coverage_chip(full).text == "all 486 targets searched"
    # Both still fit the budget.
    assert _band_height_mm(_band(partial)) <= _BAND_CLEAR_MM

    # A screen that searched everything it can answer for is not partial. The
    # RNA viruses leave the denominator because no depth of DNA sequencing
    # reaches them, and the rows that answer at a shared rank were searched.
    complete = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {
            "assessed": 420, "not_assessed": 66, "total_targets": 486,
            "searched": 454, "answerable": 454, "searched_group_level": 34,
            "out_of_assay_scope": 32,
        },
    }
    chip = _coverage_chip(complete)
    assert chip.colour == pdfpathogens.CLEAR
    assert chip.text == "all 454 targets searched"

    alarm = {
        "results": [
            _finding("bacteria.salmonella_enterica", "bacteria", "established_enteric",
                     name="Salmonella enterica"),
            _finding("bacteria.campylobacter_jejuni", "bacteria", "established_enteric",
                     name="Campylobacter jejuni"),
            _finding("bacteria.kpneu", "bacteria", "conditional_opportunist"),
        ],
        "coverage": {"assessed": 333, "not_assessed": 153},
    }
    assert pdfpathogens.summarise(alarm)["n_pathogens"] == 2
    assert _band_height_mm(_band(alarm)) <= _BAND_ALARM_MM
    text = _flat_text(_band(alarm))
    assert "2" in text and "disease-causing organisms found" in text
    # A positive result names what was found, on page 1.
    assert "Salmonella enterica" in text


def test_a_negative_result_does_not_wear_a_warning_colour() -> None:
    """An amber card over "nothing found" is a warning about nothing."""
    from openbiota.pdfpathogens import CLEAR, CLEAR_BG, WATCH, WATCH_BG

    clean = {
        "results": [
            _result("bacteria.x", "bacteria", "not_detected"),
            dict(
                _finding("bacteria.c", "bacteria", "established_enteric"),
                report_tier="unresolved_complex",
            ),
        ],
        "coverage": {"assessed": 486, "not_assessed": 0},
    }
    panel = _band(clean)[0]
    backgrounds = [c[-1] for c in panel._bkgrndcmds]  # type: ignore[attr-defined]
    # Line commands carry the colour mid-tuple, not last.
    borders = [c for cmd in panel._linecmds for c in cmd if isinstance(c, colors.Color)]  # type: ignore[attr-defined]
    assert CLEAR_BG in backgrounds, "a negative result gets the calm background"
    assert CLEAR in borders, "and the calm border"
    assert WATCH_BG not in backgrounds and WATCH not in borders

    text = _flat_text(panel)
    assert "No disease-causing organisms found" in text


def test_the_count_chips_fit_across_the_band() -> None:
    """Chips measure their own labels, so a row of them cannot overflow."""
    from reportlab.lib.units import mm
    from reportlab.pdfbase.pdfmetrics import stringWidth

    from openbiota.pdfpathogens import CLEAR, CLEAR_BG, CONTENT_WIDTH, Chip

    for label in ("2 traces, unconfirmed", "333 of 486 targets searched", "1 opportunist"):
        chip = Chip(label, CLEAR, CLEAR_BG)
        assert chip.width > 0
        # Wide enough for the text, and not a fixed guess far from it.
        assert chip.width >= stringWidth(label, "Helvetica-Bold", chip.size)
        assert chip.width <= stringWidth(label, "Helvetica-Bold", chip.size) + 6 * mm

    busy = {
        "results": [
            dict(
                _finding("bacteria.c", "bacteria", "established_enteric"),
                report_tier="unresolved_complex",
            ),
            dict(
                _finding("bacteria.d", "bacteria", "established_enteric"),
                report_tier="pathotype_negative",
            ),
        ],
        "coverage": {"assessed": 333, "not_assessed": 153},
    }
    panel = _band(busy)[0]
    text_cell = panel._cellvalues[0][1]  # type: ignore[attr-defined]
    rows = [f for f in text_cell if hasattr(f, "_argW")]
    available = CONTENT_WIDTH - 14 * mm
    for row in rows:
        assert sum(row._argW) <= available, "a chip row ran off the page"


def test_the_band_links_into_the_detail_section() -> None:
    """Brevity on page 1 is only acceptable if the way in is obvious."""
    clean = {
        "results": [_result("bacteria.x", "bacteria", "not_detected")],
        "coverage": {"assessed": 486, "not_assessed": 0},
    }
    text = _flat_text(_band(clean))
    assert "Section 7" in text



def test_strict_cdiff_rule_contract() -> None:
    """The strict rule's contract lives beside the gates it changes."""
    from openbiota.pathogens import resolve

    assert resolve.self_test() == 0


# --------------------------------------------------------------------------- #
# one assembly, one row
# --------------------------------------------------------------------------- #


def test_no_assembly_is_installed_under_two_rows(catalog: Any) -> None:
    """Two rows over one genome blind each other, so only one may hold it.

    The ownership pass decides target-specific sequence by comparing every
    installed reference against every other. An identical second copy makes
    both rows non-specific: each loses the sequence that identified it, and
    each is then reported as unresolvable for every sample. Adding
    genus-level rows once cost twelve targets their species-level call this
    way, *Giardia duodenalis*, *Cryptosporidium parvum* and *Aspergillus
    flavus* among them, with no error anywhere - the genus row's best
    assembly was the species row's assembly.
    """
    import json

    from openbiota.pathogens.refs import Resolution, deduplicate_references

    cache = CATALOG.parents[1] / "refs" / "pathogens" / "resolutions.cache.json"
    if not cache.is_file():
        pytest.skip("no resolution cache in this checkout")
    raw = json.loads(cache.read_text())
    resolutions = {
        seed_id: Resolution(**{
            k: (tuple(v) if isinstance(v, list) else v) for k, v in record.items()
        })
        for seed_id, record in raw.items()
    }
    deduplicate_references(catalog, resolutions)

    installed_by_asset: dict[str, list[str]] = {}
    for seed_id, res in resolutions.items():
        if res.route == "genome" and res.accession:
            installed_by_asset.setdefault(res.accession, []).append(seed_id)
        elif res.route == "marker" and res.marker_accessions:
            installed_by_asset.setdefault(
                ",".join(sorted(res.marker_accessions)), []
            ).append(seed_id)
    clashes = {a: m for a, m in installed_by_asset.items() if len(m) > 1}
    assert not clashes, f"one assembly installed under several rows: {clashes}"

    # Every folded row says which row carries it, so the report can state the
    # rank the answer comes at instead of calling the target unsearched.
    for seed_id, res in resolutions.items():
        if res.resolution_state == "covered_by_relative":
            assert res.covered_by, seed_id
            assert resolutions[res.covered_by].route in {"genome", "marker"}, seed_id
            assert res.covered_by != seed_id


def test_the_row_a_reader_needs_named_keeps_the_reference(catalog: Any) -> None:
    """When rows share a genome, the clinically meaningful one holds it.

    A species outranks a genus umbrella over it; a recognised human pathogen
    outranks an animal relative kept only to absorb its reads; and a row that
    resolved to a named species outranks one that fell back to the genus.
    *Cyclospora cayetanensis* is the notifiable foodborne parasite, and an
    alphabetical tie-break once handed its genome to *Cyclospora ashfordi*, a
    lineage name NCBI does not even carry.
    """
    import json

    from openbiota.pathogens.refs import Resolution, deduplicate_references

    cache = CATALOG.parents[1] / "refs" / "pathogens" / "resolutions.cache.json"
    if not cache.is_file():
        pytest.skip("no resolution cache in this checkout")
    raw = json.loads(cache.read_text())
    resolutions = {
        seed_id: Resolution(**{
            k: (tuple(v) if isinstance(v, list) else v) for k, v in record.items()
        })
        for seed_id, record in raw.items()
    }
    deduplicate_references(catalog, resolutions)

    def holds_its_own(seed_id: str) -> bool:
        return resolutions[seed_id].route in {"genome", "marker"}

    must_keep = (
        "protozoa.giardia_duodenalis",
        "protozoa.cryptosporidium_parvum",
        "protozoa.cyclospora_cayetanensis",
        "fungi.aspergillus_flavus",
        "fungi.candida_krusei",
        "fungi.candida_dubliniensis",
        "bacteria.klebsiella_pneumoniae",
        "bacteria.klebsiella_oxytoca",
        "bacteria.enterobacter_cloacae",
    )
    for seed_id in must_keep:
        if seed_id in resolutions:
            assert holds_its_own(seed_id), (
                f"{seed_id} lost its reference to a row that shares its genome"
            )
