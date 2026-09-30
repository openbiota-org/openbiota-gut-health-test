"""Shared fixtures: a synthetic reference database and synthetic DIAMOND output.

No real FASTQ and no network access. The reference database is constructed
in-memory with hand-chosen reference lengths so the normalisation arithmetic
can be verified against a value computed by hand.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from openbiota.logging_util import Reporter
from openbiota.panels import PanelSet, parse_normalizer, parse_panel
from openbiota.references import EntryStats, ReferenceDatabase, RefSequence
from openbiota.residues import AnchorMapping
from openbiota.search import BASE_FIELDS, RESIDUE_FIELDS

FIELDS = (*BASE_FIELDS, *RESIDUE_FIELDS)


def section_stubs(
    st,  # noqa: ANN001 - reportlab styles dict
    *,
    skip: frozenset[int] | set[int] = frozenset(),
    functions: tuple[str, ...] = (),
    profiles: tuple[str, ...] = (),
) -> list:
    """Destinations so a fragment of the report can be built on its own.

    Page 1 and the at-a-glance pages link to sections, function rows and
    pattern rows, and ReportLab will not save a document whose links point
    nowhere. Tests that build one piece of the report append these so the
    piece can be saved and inspected. `skip` names section numbers the story
    already contains for real; `functions` and `profiles` are the panel and
    profile names whose row destinations are needed.
    """
    from reportlab.platypus import PageBreak

    from openbiota.pdflinks import Anchor, function_dest, profile_dest, section_heading
    from openbiota.pdfreport import SECTIONS

    out: list = [PageBreak()]
    for key, number in SECTIONS.items():
        if number in skip:
            continue
        out.append(PageBreak())
        out.extend(section_heading(number, key.replace("_", " "), st["h1"]))
    for name in functions:
        out.extend((Anchor(function_dest(name)), Anchor(function_dest(name, detail=True))))
    for name in profiles:
        out.extend((Anchor(profile_dest(name)), Anchor(profile_dest(name, detail=True))))
    return out

NORMALIZER_YAML = {
    "description": "rpoB test normaliser",
    "citation": "test",
    "normalizer": {
        "id": "RPOB",
        "label": "rpoB",
        "gene": "rpoB",
        "source": "uniprot",
        "query": "(gene:rpoB) AND (reviewed:true)",
        "min_identity": 50.0,
        "length_range": [1000, 1500],
    },
}

URDA_PANEL_YAML = {
    "name": "urda",
    "description": "Urocanate reductase test panel",
    "metabolite": "Imidazole propionate (ImP)",
    "targets": [
        {
            "id": "URDA",
            "label": "urdA",
            "gene": "urdA",
            "source": "uniprot",
            "query": '(protein_name:"urocanate reductase")',
            "min_identity": 60.0,
        }
    ],
    "decoys": [
        {
            "id": "FRDA",
            "label": "frdA",
            "gene": "frdA",
            "source": "uniprot",
            "query": '(protein_name:"fumarate reductase flavoprotein subunit")',
            "min_identity": 50.0,
        }
    ],
    "residue_check": {
        "enabled": True,
        "anchor_accession": "Q8CVD0",
        "canonical_position": 373,
        "accepted_residues": ["Y", "M"],
        "applies_to": ["URDA"],
    },
    "min_fragments_for_stability": 20,
    "min_alignment_aa": 25,
    "aggregate": "sum",
    "upper_bound_reason": "most reads are too short to span residue 373",
    "interpretation": {
        "what_it_is": "Imidazole propionate, a test metabolite.",
        "made_from": "Dietary histidine",
        "higher_means": "adverse",
        "evidence_strength": "moderate",
        "summary": "A test summary of what the research reports.",
        "cognitive_evidence": "Test evidence describing the reported findings.",
        "citation": "Test et al., Journal (2026)",
        "implications": {
            "higher": (
                "More of your gut bacteria carry this gene than in most people, which is "
                "the direction associated with worse reported outcomes."
            ),
            "typical": (
                "Your capacity for this pathway is in line with most people and nothing "
                "here points toward the reported effects."
            ),
            "lower": (
                "Fewer of your gut bacteria carry this gene than in most people, which is "
                "the reassuring direction for this pathway."
            ),
        },
    },
    "citation": "test citation",
}

BUTYRATE_PANEL_YAML = {
    "name": "butyrate",
    "description": "Butyrate test panel",
    "metabolite": "Butyrate",
    "targets": [
        {
            "id": "BUT",
            "label": "but",
            "gene": "but",
            "source": "uniprot",
            "query": '(protein_name:"butyryl-CoA:acetate CoA-transferase")',
            "min_identity": 55.0,
        }
    ],
    "decoys": [],
    "aggregate": "sum",
    "citation": "test citation",
}


@pytest.fixture
def reporter() -> Reporter:
    """Silent reporter: tests must not print progress to stderr."""
    return Reporter(verbose=False)


@pytest.fixture
def panel_set() -> PanelSet:
    entry, citation, description = parse_normalizer(NORMALIZER_YAML)
    return PanelSet(
        panels=(parse_panel(URDA_PANEL_YAML), parse_panel(BUTYRATE_PANEL_YAML)),
        normalizer=entry,
        normalizer_citation=citation,
        normalizer_description=description,
    )


def _stats(
    key: str,
    panel: str,
    entry_id: str,
    role: str,
    n: int,
    mean_length: float,
    min_identity: float,
    gene: str | None = None,
    n_anchor_mapped: int = 0,
    n_anchor_residue_ok: int = 0,
    mean_length_residue_ok_aa: float = 0.0,
) -> EntryStats:
    return EntryStats(
        key=key,
        panel=panel,
        entry_id=entry_id,
        role=role,
        label=entry_id.lower(),
        gene=gene,
        source="uniprot",
        query="test",
        min_identity=min_identity,
        total_available=n,
        n_fetched=n,
        n_dropped_length=0,
        n_dropped_identical=0,
        n_dropped_claimed=0,
        n_sequences=n,
        mean_length_aa=mean_length,
        median_length_aa=mean_length,
        min_length_aa=int(mean_length),
        max_length_aa=int(mean_length),
        length_window="test",
        n_anchor_mapped=n_anchor_mapped,
        n_anchor_residue_ok=n_anchor_residue_ok,
        mean_length_residue_ok_aa=mean_length_residue_ok_aa,
    )


#: Reference layout, chosen so the arithmetic is easy to verify by hand.
#:
#:   idx 0  urda:URDA      600 aa   anchor position 300, residue OK
#:   idx 1  urda:URDA      600 aa   anchor position 300, residue NOT ok
#:   idx 2  urda:FRDA      600 aa   (decoy)
#:   idx 3  butyrate:BUT   450 aa
#:   idx 4  _normalizer:RPOB 1200 aa
#:   idx 5  _normalizer:RPOB 1200 aa
@pytest.fixture
def database(tmp_path: Path) -> ReferenceDatabase:
    sequences = [
        RefSequence(0, "urda:URDA", "URDA", "target", "urda", "P00001",
                    "Streptococcus pasteurianus", 600, "Bacillota", "Streptococcus",
                    "Streptococcaceae", anchor_pos=300, anchor_residue_ok=True),
        RefSequence(1, "urda:URDA", "URDA", "target", "urda", "P00002",
                    "Eggerthella lenta", 600, "Actinomycetota", "Eggerthella",
                    "Eggerthellaceae", anchor_pos=300, anchor_residue_ok=False),
        RefSequence(2, "urda:FRDA", "FRDA", "decoy", "urda", "P00003",
                    "Escherichia coli", 600, "Pseudomonadota", "Escherichia",
                    "Enterobacteriaceae"),
        RefSequence(3, "butyrate:BUT", "BUT", "target", "butyrate", "P00004",
                    "Faecalibacterium prausnitzii", 450, "Bacillota", "Faecalibacterium",
                    "Oscillospiraceae"),
        RefSequence(4, "_normalizer:RPOB", "RPOB", "normalizer", "_normalizer", "P00005",
                    "Bacteroides fragilis", 1200, "Bacteroidota", "Bacteroides",
                    "Bacteroidaceae"),
        RefSequence(5, "_normalizer:RPOB", "RPOB", "normalizer", "_normalizer", "P00006",
                    "Roseburia intestinalis", 1200, "Bacillota", "Roseburia",
                    "Lachnospiraceae"),
    ]
    entries = {
        "urda:URDA": _stats(
            "urda:URDA", "urda", "URDA", "target", 2, 600.0, 60.0, "urdA",
            n_anchor_mapped=2, n_anchor_residue_ok=1, mean_length_residue_ok_aa=600.0,
        ),
        "urda:FRDA": _stats("urda:FRDA", "urda", "FRDA", "decoy", 1, 600.0, 50.0, "frdA"),
        "butyrate:BUT": _stats("butyrate:BUT", "butyrate", "BUT", "target", 1, 450.0, 55.0, "but"),
        "_normalizer:RPOB": _stats(
            "_normalizer:RPOB", "_normalizer", "RPOB", "normalizer", 2, 1200.0, 50.0, "rpoB"
        ),
    }
    anchors = {
        "urda": AnchorMapping(
            panel="urda",
            anchor_accession="Q8CVD0",
            anchor_length=582,
            canonical_position=373,
            canonical_residue="Y",
            accepted_residues=("M", "Y"),
            positions={0: 300, 1: 300},
            residue_ok={0: True, 1: False},
            n_candidates=2,
            n_aligned=2,
            n_position_spanned=2,
            n_reference_residue_ok=1,
        )
    }
    return ReferenceDatabase(
        fingerprint="testfp",
        dmnd_path=tmp_path / "db.dmnd",
        fasta_path=tmp_path / "db.faa",
        index_path=tmp_path / "index.json",
        sequences=sequences,
        entries=entries,
        diamond_version="diamond version 2.2.6 (test)",
        built_at="2026-01-01T00:00:00+0000",
        anchors=anchors,
    )


def hit_line(
    qseqid: str,
    idx: int,
    entry_id: str,
    accession: str,
    *,
    pident: float = 90.0,
    aln_len: int = 45,
    sstart: int = 100,
    send: int = 144,
    bitscore: float = 95.0,
    qseq_gapped: str = "",
    sseq_gapped: str = "",
) -> str:
    """Build one synthetic ``diamond -f 6`` line matching ``FIELDS``."""
    return "\t".join(
        [
            qseqid,
            f"{idx}~{entry_id}~{accession}",
            f"{pident}",
            f"{aln_len}",
            f"{sstart}",
            f"{send}",
            f"{bitscore}",
            qseq_gapped,
            sseq_gapped,
        ]
    )
