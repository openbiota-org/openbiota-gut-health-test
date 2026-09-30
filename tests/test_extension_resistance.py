"""A11 resistance overview — BUILD_SPEC_v0.8.3 section 7.3, AT054.

The arithmetic that matters here is what must *not* be added together:
ciprofloxacin inside fluoroquinolones, an allele counted twice as a family, a
gene read as a phenotype, or an unresolved homolog counted as a resistance
call.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import resistance as RES

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
SAMPLES = ["SAMPLE2_A02", "SAMPLE4_A04", "SAMPLE6_A06", "SAMPLE1_A01", "SAMPLE3_A03"]


def _pathogens(sample: str) -> dict[str, Any]:
    path = RESULTS / sample / "results.json"
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))["pathogens"]


def _det(**kw) -> RES.Determinant:
    base = {
        "determinant_id": "amr.x", "label": "X", "gene_family": "tet(M)",
        "source_class": "tetracycline", "status": "supported_intact_sequence",
        "fragments": 10, "identity": 1.0, "covered_fraction": 0.9,
        "linked_targets": (), "host_linkage": "unlinked",
        "phenotype_inference": "not_inferred",
    }
    return RES.Determinant(**(base | kw))


# --------------------------------------------------------------------------- #
# AT054 — the three levels, and deduplication
# --------------------------------------------------------------------------- #


def test_at054_a_drug_is_not_counted_again_as_its_own_class() -> None:
    """Ciprofloxacin lives inside fluoroquinolones, as a subordinate view."""
    by_id = {c["class_id"]: c for c in RES.CLASS_DEFINITIONS}
    assert "ciprofloxacin" not in by_id
    fq = by_id["fluoroquinolone"]
    assert fq["subordinate_drug_view"] == "ciprofloxacin"
    assert "ciprofloxacin" in fq["drugs"]
    # Every drug named belongs to exactly one class definition.
    seen: dict[str, str] = {}
    for spec in RES.CLASS_DEFINITIONS:
        for drug in spec.get("drugs") or ():
            assert drug not in seen, f"{drug} claimed by {seen.get(drug)} and {spec['class_id']}"
            seen[drug] = str(spec["class_id"])


def test_at054_enzyme_classes_sit_inside_their_antibiotic_class() -> None:
    """Ambler A/B/C/D are an enzyme axis, not four antibiotic classes."""
    beta = next(c for c in RES.CLASS_DEFINITIONS if c["class_id"] == "beta_lactam")
    enzymes = beta["enzyme_classes"]
    assert "carbapenemase" in enzymes and "esbl_or_ampc" in enzymes
    assert "Ambler" in " ".join(enzymes.values())
    # And they are not separate top-level classes.
    ids = {c["class_id"] for c in RES.CLASS_DEFINITIONS}
    assert "carbapenemase" not in ids and "esbl_or_ampc" not in ids


def test_at054_alleles_of_one_family_count_once_for_richness() -> None:
    cls = RES.ResistanceClass(
        class_id="tetracycline", label="Tetracyclines", mechanism="m", drugs=(),
        determinants=(
            _det(determinant_id="amr.tet_a", gene_family="tet(M)"),
            _det(determinant_id="amr.tet_b", gene_family="tet(M)"),
            _det(determinant_id="amr.tet_c", gene_family="tet(O)"),
        ),
    )
    assert len(cls.detected) == 3
    assert cls.richness == 2, "two copies of tet(M) are one family"


def test_at054_cross_class_totals_do_not_sum_per_class_richness() -> None:
    """A family appearing in two classes counts once overall."""
    shared = _det(gene_family="shared_family", source_class="tetracycline")
    other = _det(determinant_id="amr.y", gene_family="shared_family",
                 source_class="aminoglycoside")
    classes = RES.summarise([shared, other])
    per_class = sum(c.richness for c in classes)
    overall = RES.totals(classes)["total_determinant_richness"]
    assert per_class == 2
    assert overall == 1


# --------------------------------------------------------------------------- #
# evidence kinds stay apart
# --------------------------------------------------------------------------- #


def test_acquired_gene_mutation_and_unresolved_homolog_are_three_states() -> None:
    assert _det(status="supported_intact_sequence").evidence_kind == "acquired_gene"
    assert _det(status="supported_resistance_mutation").evidence_kind == "validated_resistance_mutation"
    assert _det(status="candidate_homolog").evidence_kind == "unresolved_homolog"
    assert _det(status="not_detected").evidence_kind == "not_detected"
    assert _det(status="not_assessed").evidence_kind == "not_assessed"
    assert not _det(status="not_detected").detected
    assert not _det(status="not_assessed").detected
    assert _det(status="candidate_homolog").detected


def test_an_unresolved_homolog_is_visible_as_a_fraction_of_the_detections() -> None:
    cls = RES.ResistanceClass(
        class_id="aminoglycoside", label="Aminoglycosides", mechanism="m", drugs=(),
        determinants=(
            _det(determinant_id="amr.aph", gene_family="aph"),
            _det(determinant_id="amr.aac", gene_family="aac", status="candidate_homolog"),
        ),
    )
    payload = cls.to_json()
    assert payload["n_acquired_genes"] == 1
    assert payload["n_unresolved_homologs"] == 1
    assert payload["unresolved_fraction_of_detected"] == pytest.approx(0.5)


def test_an_unassessed_determinant_is_neither_present_nor_absent() -> None:
    cls = RES.ResistanceClass(
        class_id="colistin", label="Polymyxins", mechanism="m", drugs=(),
        determinants=(_det(status="not_assessed"), _det(determinant_id="amr.z", status="not_detected")),
    )
    payload = cls.to_json()
    assert payload["n_detected"] == 0
    assert payload["n_assessed"] == 1
    assert payload["n_not_assessed"] == 1
    assert payload["unresolved_fraction_of_detected"] is None


# --------------------------------------------------------------------------- #
# a gene is not a phenotype
# --------------------------------------------------------------------------- #


def test_no_class_or_metric_infers_a_clinical_phenotype() -> None:
    cls = RES.ResistanceClass(
        class_id="tetracycline", label="Tetracyclines", mechanism="m", drugs=(),
        determinants=(_det(),),
    )
    joined = " ".join(cls.to_json()["limitations"])
    assert "A gene is not a phenotype" in joined
    assert "does not establish that" in joined
    for metric in RES.metrics([cls]):
        text = " ".join(metric.limitations)
        assert "not a phenotype" in text
        assert metric.direction == "descriptive"
        assert metric.reference_percentile is None


def test_relative_abundance_is_not_converted_to_absolute_load() -> None:
    cls = RES.ResistanceClass(
        class_id="tetracycline", label="Tetracyclines", mechanism="m", drugs=(),
        determinants=(_det(),),
    )
    for metric in RES.metrics([cls]):
        assert metric.denominator
        assert "per gram" not in (metric.unit or "")
        assert any("absolute-load" in x for x in metric.limitations)


def test_carriage_language_does_not_name_an_organism_without_linkage() -> None:
    """Short reads rarely link a gene to its carrier; saying otherwise is the
    commonest overstatement in this part of a report."""
    unlinked = RES.summarise([_det()])
    note = RES.unlinked_note(unlinked)
    assert "not which organism carries it" in note

    linked = RES.summarise([_det(linked_targets=("bacteria.escherichia_coli",))])
    note2 = RES.unlinked_note(linked)
    assert "could be linked" in note2

    none = RES.summarise([_det(status="not_detected")])
    assert "No resistance determinant was detected" in RES.unlinked_note(none)


# --------------------------------------------------------------------------- #
# against real samples
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_summary_runs_and_reconciles_on_a_real_sample(sample: str) -> None:
    view = RES.build(_pathogens(sample))
    totals = view["totals"]
    assert view["feature_id"] == "A11"
    assert totals["n_classes_reported"] == len(RES.CLASS_DEFINITIONS)

    # Every detected determinant appears in exactly one class.
    seen: set[str] = set()
    detected_rows = 0
    for cls in view["classes"]:
        for det in cls["determinants"]:
            assert det["determinant_id"] not in seen, "a determinant is in two classes"
            seen.add(det["determinant_id"])
            detected_rows += 1
    assert detected_rows == totals["n_determinants_detected"]

    # Class fragment counts sum to the total, because no class shares a row.
    assert sum(c["supporting_fragments"] for c in view["classes"]) == \
        totals["total_supporting_fragments"]

    # Overall richness never exceeds the sum of per-class richness.
    assert totals["total_determinant_richness"] <= sum(
        c["determinant_richness"] for c in view["classes"]
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_existing_determinant_calls_are_read_not_recomputed(sample: str) -> None:
    pathogens = _pathogens(sample)
    source = {
        str(d["determinant_id"]): d for d in pathogens["determinants"]
        if str(d.get("kind")) == "amr"
    }
    for determinant in RES.read_determinants(pathogens):
        original = source[determinant.determinant_id]
        assert determinant.status == original["determinant_status"]
        assert determinant.fragments == int(original["supporting_fragments"] or 0)
        assert determinant.gene_family == str(original["gene_family"] or "")
        # The screen's own phenotype stance is carried through untouched.
        assert determinant.phenotype_inference == str(
            original.get("clinical_phenotype_inference") or "not_inferred"
        )


def test_every_required_class_from_the_specification_is_present() -> None:
    required = {
        "beta_lactam", "macrolide_lincosamide_streptogramin", "fluoroquinolone",
        "glycopeptide", "sulfonamide_trimethoprim", "polymyxin", "aminoglycoside",
    }
    assert required <= {c["class_id"] for c in RES.CLASS_DEFINITIONS}
    for spec in RES.CLASS_DEFINITIONS:
        assert spec["mechanism"], f"{spec['class_id']} has no mechanism described"
