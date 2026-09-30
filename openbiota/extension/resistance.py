"""A11 — the resistance overview: classes, mechanisms and carriers.

BUILD_SPEC_v0.8.3 section 7.3. This reads the existing AMR determinant calls
and adds the class-level view they lack. It changes none of them.

Four rules the specification is specific about, and which the arithmetic here
enforces:

* **Enzyme class, antibiotic class and drug are three levels.** Ciprofloxacin
  is a drug inside the fluoroquinolone class; counting it again as its own
  class double-counts the same determinant. Drugs are a subordinate view with
  their own count, never added to the class total.
* **Alleles and parent/child terms deduplicate.** Determinant richness counts
  distinct sequence families, so `tet(M)` and `tet(O)` are two and two copies
  of `tet(M)` are one.
* **A gene is not a phenotype.** Nothing here infers clinical resistance. A
  wild-type housekeeping hit and a species' intrinsic profile are not
  evidence that an antibiotic would fail.
* **An unresolved homolog is its own state.** Acquired gene, validated
  resistance mutation and unresolved homolog are counted separately, because
  a partial hit to an aminoglycoside-modifying family is not a resistance
  gene call.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.extension.schema import ExtensionMetric, fingerprint

METHOD: Final = "ext083.resistance_summary/1.0"

#: Determinant states that count as a detection, and what kind. The existing
#: screen's vocabulary, mapped once here rather than re-interpreted per view.
ACQUIRED_STATES: Final[frozenset[str]] = frozenset({"supported_intact_sequence"})
UNRESOLVED_STATES: Final[frozenset[str]] = frozenset({"candidate_homolog", "partial_sequence"})
MUTATION_STATES: Final[frozenset[str]] = frozenset({"supported_resistance_mutation"})
ABSENT_STATES: Final[frozenset[str]] = frozenset({"not_detected"})
UNASSESSED_STATES: Final[frozenset[str]] = frozenset({"not_assessed", "not_run"})

#: The antibiotic classes section 7.3 requires, mapped from the existing
#: screen's own class vocabulary. Beta-lactamase classes A/B/C/D are the
#: enzyme level and are reported inside their antibiotic class rather than
#: beside it: an enzyme class and a drug class are not the same axis.
CLASS_DEFINITIONS: Final[tuple[dict[str, Any], ...]] = (
    {
        "class_id": "beta_lactam",
        "label": "Beta-lactams",
        "source_classes": ("carbapenemase", "esbl_or_ampc", "methicillin"),
        "enzyme_classes": {
            "carbapenemase": "Ambler A/B/D carbapenemases",
            "esbl_or_ampc": "Ambler A extended-spectrum and C cephalosporinases",
            "methicillin": "altered penicillin-binding protein (mecA/mecC)",
        },
        "mechanism": "enzymatic hydrolysis of the beta-lactam ring, or an altered target",
        "drugs": ("penicillins", "cephalosporins", "carbapenems"),
    },
    {
        "class_id": "macrolide_lincosamide_streptogramin",
        "label": "Macrolides, lincosamides and streptogramin B",
        "source_classes": ("macrolide_lincosamide",),
        "mechanism": "ribosomal methylation (erm), efflux (mef) or enzymatic inactivation (lnu)",
        "drugs": ("erythromycin", "clarithromycin", "azithromycin", "clindamycin"),
    },
    {
        "class_id": "fluoroquinolone",
        "label": "Fluoroquinolones",
        "source_classes": ("quinolone",),
        "mechanism": "target-site mutation, target protection (qnr) or efflux",
        "drugs": ("ciprofloxacin", "levofloxacin", "moxifloxacin"),
        "subordinate_drug_view": "ciprofloxacin",
    },
    {
        "class_id": "glycopeptide",
        "label": "Glycopeptides",
        "source_classes": ("glycopeptide",),
        "mechanism": "remodelled peptidoglycan precursor (van operons)",
        "drugs": ("vancomycin", "teicoplanin"),
    },
    {
        "class_id": "sulfonamide_trimethoprim",
        "label": "Sulfonamides and trimethoprim",
        "source_classes": ("sulfonamide_trimethoprim",),
        "mechanism": "alternative dihydropteroate synthase or dihydrofolate reductase",
        "drugs": ("sulfamethoxazole", "trimethoprim"),
    },
    {
        "class_id": "polymyxin",
        "label": "Polymyxins and colistin",
        "source_classes": ("colistin",),
        "mechanism": "lipid A modification (mcr, chromosomal two-component systems)",
        "drugs": ("colistin", "polymyxin B"),
    },
    {
        "class_id": "aminoglycoside",
        "label": "Aminoglycosides",
        "source_classes": ("aminoglycoside",),
        "mechanism": "enzymatic modification (aac, aph, ant) or 16S methylation (armA)",
        "drugs": ("gentamicin", "amikacin", "tobramycin"),
    },
    {
        "class_id": "tetracycline",
        "label": "Tetracyclines",
        "source_classes": ("tetracycline",),
        "mechanism": "ribosomal protection, efflux or enzymatic inactivation",
        "drugs": ("tetracycline", "doxycycline"),
    },
    {
        "class_id": "oxazolidinone_phenicol",
        "label": "Oxazolidinones and phenicols",
        "source_classes": ("oxazolidinone_phenicol",),
        "mechanism": "23S rRNA methylation (cfr) or enzymatic inactivation (cat)",
        "drugs": ("linezolid", "chloramphenicol"),
    },
    {
        "class_id": "other_mechanism",
        "label": "Other mechanisms",
        "source_classes": ("other_mechanism",),
        "mechanism": "mechanisms outside the classes above",
        "drugs": (),
    },
)


@dataclass(slots=True)
class Determinant:
    """One resistance determinant, read from the existing screen unchanged."""

    determinant_id: str
    label: str
    gene_family: str
    source_class: str
    status: str
    fragments: int
    identity: float | None
    covered_fraction: float | None
    linked_targets: tuple[str, ...]
    host_linkage: str
    phenotype_inference: str
    reason_codes: tuple[str, ...] = ()

    @property
    def detected(self) -> bool:
        return self.status in ACQUIRED_STATES | MUTATION_STATES | UNRESOLVED_STATES

    @property
    def evidence_kind(self) -> str:
        """Acquired gene, resistance mutation or unresolved homolog."""
        if self.status in ACQUIRED_STATES:
            return "acquired_gene"
        if self.status in MUTATION_STATES:
            return "validated_resistance_mutation"
        if self.status in UNRESOLVED_STATES:
            return "unresolved_homolog"
        if self.status in UNASSESSED_STATES:
            return "not_assessed"
        return "not_detected"

    def to_json(self) -> dict[str, Any]:
        return {
            "determinant_id": self.determinant_id,
            "label": self.label,
            "gene_family": self.gene_family,
            "evidence_kind": self.evidence_kind,
            "status": self.status,
            "supporting_fragments": self.fragments,
            "aligned_identity": self.identity,
            "covered_reference_fraction": self.covered_fraction,
            "carrier_evidence": (
                "linked to " + ", ".join(self.linked_targets)
                if self.linked_targets else f"unlinked ({self.host_linkage})"
            ),
            "linked_target_ids": list(self.linked_targets),
            "clinical_phenotype_inference": self.phenotype_inference,
            "reason_codes": list(self.reason_codes),
        }


@dataclass(slots=True)
class ResistanceClass:
    """One antibiotic class, with its determinants and their evidence."""

    class_id: str
    label: str
    mechanism: str
    drugs: tuple[str, ...]
    determinants: tuple[Determinant, ...]
    enzyme_classes: Mapping[str, str] = field(default_factory=dict)
    subordinate_drug_view: str | None = None

    @property
    def detected(self) -> tuple[Determinant, ...]:
        return tuple(d for d in self.determinants if d.detected)

    @property
    def acquired(self) -> tuple[Determinant, ...]:
        return tuple(d for d in self.determinants if d.evidence_kind == "acquired_gene")

    @property
    def unresolved(self) -> tuple[Determinant, ...]:
        return tuple(d for d in self.determinants if d.evidence_kind == "unresolved_homolog")

    @property
    def richness(self) -> int:
        """Distinct sequence families detected, alleles deduplicated.

        Two copies of `tet(M)` are one family; `tet(M)` and `tet(O)` are two.
        Counting alleles here would make a class look broader than it is.
        """
        return len({d.gene_family.lower() for d in self.detected if d.gene_family})

    @property
    def fragments(self) -> int:
        return sum(d.fragments for d in self.detected)

    @property
    def assessed(self) -> int:
        return sum(1 for d in self.determinants if d.status not in UNASSESSED_STATES)

    @property
    def unresolved_fraction(self) -> float | None:
        if not self.detected:
            return None
        return len(self.unresolved) / len(self.detected)

    def to_json(self) -> dict[str, Any]:
        return {
            "class_id": self.class_id,
            "label": self.label,
            "mechanism": self.mechanism,
            "enzyme_classes": dict(self.enzyme_classes),
            "drugs_in_this_class": list(self.drugs),
            "subordinate_drug_view": self.subordinate_drug_view,
            "n_determinants_searched": len(self.determinants),
            "n_assessed": self.assessed,
            "n_not_assessed": len(self.determinants) - self.assessed,
            "n_detected": len(self.detected),
            "determinant_richness": self.richness,
            "richness_definition": (
                "distinct sequence families detected; alleles and duplicate family "
                "names count once"
            ),
            "supporting_fragments": self.fragments,
            "fragment_denominator": "quality_passed_nonhost_fragments",
            "n_acquired_genes": len(self.acquired),
            "n_validated_mutations": sum(
                1 for d in self.determinants if d.evidence_kind == "validated_resistance_mutation"
            ),
            "n_unresolved_homologs": len(self.unresolved),
            "unresolved_fraction_of_detected": self.unresolved_fraction,
            "determinants": [d.to_json() for d in sorted(
                self.determinants, key=lambda d: (-d.fragments, d.determinant_id)
            ) if d.detected],
            "limitations": [
                "A gene is not a phenotype. Carrying a determinant does not establish that "
                "an antibiotic would fail, and its absence here does not establish that one "
                "would work.",
                "An unresolved homolog is a partial match to a family, not a resistance call.",
            ],
        }


def read_determinants(pathogens: Mapping[str, Any] | None) -> tuple[Determinant, ...]:
    """Read the existing AMR determinants. Nothing is recomputed."""
    out: list[Determinant] = []
    for row in (pathogens or {}).get("determinants") or ():
        if str(row.get("kind")) != "amr":
            continue
        out.append(Determinant(
            determinant_id=str(row.get("determinant_id") or ""),
            label=str(row.get("display_name") or row.get("determinant_id") or ""),
            gene_family=str(row.get("gene_family") or ""),
            source_class=str(row.get("amr_class") or "other_mechanism"),
            status=str(row.get("determinant_status") or "not_assessed"),
            fragments=int(row.get("supporting_fragments") or 0),
            identity=row.get("aligned_identity"),
            covered_fraction=row.get("covered_reference_fraction"),
            linked_targets=tuple(row.get("linked_target_ids") or ()),
            host_linkage=str(row.get("host_linkage") or "unlinked"),
            phenotype_inference=str(row.get("clinical_phenotype_inference") or "not_inferred"),
            reason_codes=tuple(row.get("reason_codes") or ()),
        ))
    return tuple(out)


def summarise(determinants: Sequence[Determinant]) -> tuple[ResistanceClass, ...]:
    """Group determinants into the required antibiotic classes."""
    by_source: dict[str, list[Determinant]] = {}
    for d in determinants:
        by_source.setdefault(d.source_class, []).append(d)
    classes: list[ResistanceClass] = []
    for spec in CLASS_DEFINITIONS:
        members: list[Determinant] = []
        for source in spec["source_classes"]:
            members.extend(by_source.get(source, ()))
        classes.append(ResistanceClass(
            class_id=str(spec["class_id"]),
            label=str(spec["label"]),
            mechanism=str(spec["mechanism"]),
            drugs=tuple(spec.get("drugs") or ()),
            determinants=tuple(members),
            enzyme_classes=dict(spec.get("enzyme_classes") or {}),
            subordinate_drug_view=spec.get("subordinate_drug_view"),
        ))
    return tuple(classes)


def totals(classes: Sequence[ResistanceClass]) -> dict[str, Any]:
    """Cross-class totals that do not double-count.

    Richness is over the union of gene families, not the sum of per-class
    richness: a family appearing in two classes is one family. And a drug
    view such as ciprofloxacin is never added to its parent class's total.
    """
    families: set[str] = set()
    detected_ids: set[str] = set()
    fragments = 0
    for c in classes:
        for d in c.detected:
            if d.gene_family:
                families.add(d.gene_family.lower())
            detected_ids.add(d.determinant_id)
            fragments += d.fragments
    return {
        "n_classes_reported": len(classes),
        "n_classes_with_a_detection": sum(1 for c in classes if c.detected),
        "total_determinant_richness": len(families),
        "richness_definition": (
            "distinct gene families across all classes; a family counted in two "
            "classes counts once, and a drug view is never added to its class"
        ),
        "n_determinants_detected": len(detected_ids),
        "total_supporting_fragments": fragments,
        "fragment_denominator": "quality_passed_nonhost_fragments",
    }


def metrics(
    classes: Sequence[ResistanceClass],
    *,
    input_fingerprint: str | None = None,
) -> tuple[ExtensionMetric, ...]:
    """Emit the class-level abundance and richness indices (M012, M013)."""
    fp = input_fingerprint or fingerprint([c.class_id for c in classes])
    summary = totals(classes)
    limits = (
        "Sequence presence, not a phenotype: this does not establish that any "
        "antibiotic would fail, and an absence does not establish that one would work.",
        "Relative sequence abundance. Converting it to organisms or gene copies per gram "
        "needs an absolute-load assay.",
    )
    out = [
        ExtensionMetric(
            metric_id="ext083.resistance.determinant_richness",
            label="Resistance determinant richness",
            kind="gene_abundance",
            state="measured",
            value=float(summary["total_determinant_richness"]),
            unit="distinct_gene_families",
            denominator="determinants searched by the existing screen",
            direction="descriptive",
            method_id=METHOD,
            input_fingerprint=fp,
            group="resistance",
            feature_id="A11",
            view_ids=("M013",),
            analytical_confidence="supported",
            evidence_maturity="characterized_enzyme",
            source_ids=("D04", "D05", "D06"),
            limitations=limits,
            extra=summary,
        ),
        ExtensionMetric(
            metric_id="ext083.resistance.fragment_abundance",
            label="Resistance determinant abundance",
            kind="gene_abundance",
            state="measured",
            value=float(summary["total_supporting_fragments"]),
            unit="fragments",
            denominator="quality_passed_nonhost_fragments",
            direction="descriptive",
            method_id=METHOD,
            input_fingerprint=fp,
            group="resistance",
            feature_id="A11",
            view_ids=("M012",),
            analytical_confidence="supported",
            evidence_maturity="characterized_enzyme",
            source_ids=("D04", "D05", "D06"),
            limitations=limits,
            extra=summary,
        ),
    ]
    for c in classes:
        if not c.detected:
            continue
        out.append(ExtensionMetric(
            metric_id=f"ext083.resistance.class.{c.class_id}",
            label=f"{c.label}: determinant richness",
            kind="gene_abundance",
            state="measured",
            value=float(c.richness),
            unit="distinct_gene_families",
            denominator=f"{len(c.determinants)} determinants searched in this class",
            direction="descriptive",
            method_id=METHOD,
            input_fingerprint=fp,
            group="resistance",
            feature_id="A11",
            analytical_confidence="supported" if c.acquired else "provisional",
            evidence_maturity="characterized_enzyme",
            source_ids=("D04", "D05", "D06"),
            limitations=limits,
            extra=c.to_json(),
        ))
    return tuple(out)


def build(pathogens: Mapping[str, Any] | None) -> dict[str, Any]:
    """The whole A11 view for the report and the extension result."""
    determinants = read_determinants(pathogens)
    classes = summarise(determinants)
    return {
        "schema_version": "openbiota.resistance-summary/1.0",
        "feature_id": "A11",
        "method_id": METHOD,
        "source_ids": ["D04", "D05", "D06"],
        "totals": totals(classes),
        "classes": [c.to_json() for c in classes],
        "n_determinants_read": len(determinants),
        "database_note": (
            "Read from the existing determinant screen and grouped; no call is "
            "recomputed and no phenotype is inferred."
        ),
        "limitations": [
            "Enzyme class, antibiotic class and individual drug are three levels. A drug "
            "view is subordinate to its class and is never added to the class total.",
            "Wild-type housekeeping sequence and a species' intrinsic profile are not "
            "evidence of acquired resistance.",
            "Absolute load per gram needs a validated qPCR, flow-cytometry or spike-in "
            "assay; relative sequence abundance cannot be converted to it.",
        ],
    }


def unlinked_note(classes: Iterable[ResistanceClass]) -> str:
    """One sentence on carriage, for the report.

    Short-read data rarely links a resistance gene to the organism carrying
    it, and saying which organism is resistant without that link is the most
    common overstatement in this part of a report.
    """
    detected = [d for c in classes for d in c.detected]
    if not detected:
        return "No resistance determinant was detected above the screen's threshold."
    linked = [d for d in detected if d.linked_targets]
    if not linked:
        return (
            f"All {len(detected)} detected determinants are unlinked: the reads show the "
            "gene is present in the community, not which organism carries it."
        )
    return (
        f"{len(linked)} of {len(detected)} detected determinants could be linked to a "
        "specific organism; the rest are present in the community without a resolved carrier."
    )


__all__ = [
    "ABSENT_STATES",
    "ACQUIRED_STATES",
    "CLASS_DEFINITIONS",
    "MUTATION_STATES",
    "UNASSESSED_STATES",
    "UNRESOLVED_STATES",
    "Determinant",
    "ResistanceClass",
    "build",
    "metrics",
    "read_determinants",
    "summarise",
    "totals",
    "unlinked_note",
]
