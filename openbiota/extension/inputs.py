"""Your information, inputs and assumptions — A14, BUILD_SPEC_v0.8.3 §10.

One register of everything that went into the report: what the person
supplied, what the report assumed because they did not, what it derived,
what it recorded but does not yet use, and what it never had. Section 25
displays this register; the domain sections show what each input changed.

Three rules drive the design, and each of them exists because the opposite
is the easy mistake.

**An assumption is never a fact.** `input_state` keeps `supplied`,
`derived`, `inherited_default`, `assumed_default`, `not_supplied`,
`not_applicable` and `conflicting` apart, and nothing collapses them.
"Assumed for this calculation" cannot become "patient confirmed" because
they are different states with different provenance fields.

**A laboratory value has no imputed normal.** An absent assay is
`not_supplied` with a null value. It never becomes a normal result, and no
normal/abnormal classification is generated for an assay nobody ran.

**A stored field with no consumer says so.** `affects` is the list of
things that actually read the value. When it is empty the register prints
*Recorded; not used in this report*, because importing a number is not the
same as using it.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

SCHEMA_VERSION: Final = "openbiota.input-register/1.0"

#: Where a value came from, and how much weight it carries.
#:
#: ``supplied``           the person or the laboratory gave it
#: ``derived``            computed from other supplied facts, with the rule named
#: ``inherited_default``  the production application's existing default
#: ``assumed_default``    a new capability's declared fallback for this run
#: ``not_supplied``       nobody gave it and no default applies
#: ``not_applicable``     the question does not arise for this sample
#: ``conflicting``        two supplied records disagree and neither was overwritten
INPUT_STATES: Final[frozenset[str]] = frozenset({
    "supplied", "derived", "inherited_default", "assumed_default",
    "not_supplied", "not_applicable", "conflicting",
})

#: States a reader should read as "this is a fact about you".
FACTUAL_STATES: Final[frozenset[str]] = frozenset({"supplied", "derived"})

#: The categories of §10.2, in the order the register presents them.
CATEGORIES: Final[tuple[tuple[str, str], ...]] = (
    ("participant_specimen", "Participant and specimen"),
    ("symptoms_goals", "Symptoms and goals"),
    ("diet_tolerance", "Diet and tolerance"),
    ("medications_supplements", "Medications and supplements"),
    ("clinical_history", "Clinical history"),
    ("collection_processing", "Collection and processing"),
    ("longitudinal_exposures", "Longitudinal exposures"),
    ("laboratory", "Optional laboratory results"),
    ("analysis_configuration", "Analysis configuration"),
)
CATEGORY_LABELS: Final[Mapping[str, str]] = dict(CATEGORIES)

#: How an input reaches a result. The distinction matters: a value that
#: changes a number is not the same as one that changes a sentence.
INFLUENCE_KINDS: Final[frozenset[str]] = frozenset({
    "changes_calculation", "changes_interpretation", "adds_context", "none",
})


class InputRegisterError(ValueError):
    """A record that would misstate where its value came from."""


# --------------------------------------------------------------------------- #
# influence
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Influence:
    """One thing an input reaches, and how."""

    target: str
    rule: str
    kind: str
    explanation: str

    def __post_init__(self) -> None:
        if self.kind not in INFLUENCE_KINDS:
            raise InputRegisterError(f"{self.target}: unknown influence kind {self.kind!r}")

    def to_json(self) -> dict[str, Any]:
        return {
            "target": self.target, "consuming_rule": self.rule,
            "kind": self.kind, "explanation": self.explanation,
        }


# --------------------------------------------------------------------------- #
# one input
# --------------------------------------------------------------------------- #


@dataclass
class InputRecord:
    """One fact, assumption or absence, with its provenance kept separate."""

    input_id: str
    category: str
    label: str
    input_state: str
    reported_value: Any = None
    effective_value: Any = None
    unit: str | None = None
    source: str | None = None
    recorded_at: str | None = None
    valid_at: str | None = None
    default_id: str | None = None
    default_version: str | None = None
    default_rationale: str | None = None
    derivation: str | None = None
    conflict: tuple[str, ...] = ()
    affects: tuple[Influence, ...] = ()
    linked_results: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    participant_id: str | None = None
    sample_id: str | None = None

    def __post_init__(self) -> None:
        if self.input_state not in INPUT_STATES:
            raise InputRegisterError(
                f"{self.input_id}: unknown input state {self.input_state!r}"
            )
        if self.category not in CATEGORY_LABELS:
            raise InputRegisterError(f"{self.input_id}: unknown category {self.category!r}")
        # An assumption that presents itself as a supplied fact is the single
        # failure this register exists to prevent.
        if self.input_state in {"assumed_default", "inherited_default"}:
            if self.reported_value is not None:
                raise InputRegisterError(
                    f"{self.input_id}: state {self.input_state!r} carries a reported value "
                    f"({self.reported_value!r}); a default is not something the person said"
                )
            if not self.default_rationale:
                raise InputRegisterError(
                    f"{self.input_id}: a default must say why it is the default"
                )
        if self.input_state == "supplied" and self.reported_value is None:
            raise InputRegisterError(
                f"{self.input_id}: 'supplied' with no reported value is 'not_supplied'"
            )
        if self.input_state == "not_supplied" and self.effective_value is not None:
            raise InputRegisterError(
                f"{self.input_id}: 'not_supplied' cannot carry an effective value "
                f"({self.effective_value!r}); that is a default, and it must say so"
            )
        if self.input_state == "derived" and not self.derivation:
            raise InputRegisterError(f"{self.input_id}: a derived value must name its derivation")
        if self.input_state == "conflicting" and len(self.conflict) < 2:
            raise InputRegisterError(
                f"{self.input_id}: 'conflicting' must record the records that disagree"
            )

    @property
    def is_fact(self) -> bool:
        return self.input_state in FACTUAL_STATES

    @property
    def is_used(self) -> bool:
        return any(i.kind != "none" for i in self.affects)

    @property
    def provenance_words(self) -> str:
        """What the register prints in its 'supplied or assumed' column."""
        return {
            "supplied": "Supplied",
            "derived": "Worked out from what you supplied",
            "inherited_default": "Existing default",
            "assumed_default": "Assumed for this run",
            "not_supplied": "Not supplied",
            "not_applicable": "Does not apply",
            "conflicting": "Conflicting records",
        }[self.input_state]

    def to_json(self) -> dict[str, Any]:
        return {
            "input_id": self.input_id,
            "category": self.category,
            "category_label": CATEGORY_LABELS[self.category],
            "label": self.label,
            "participant_id": self.participant_id,
            "sample_id": self.sample_id,
            "reported_value": self.reported_value,
            "effective_value": self.effective_value,
            "unit": self.unit,
            "input_state": self.input_state,
            "source": self.source,
            "recorded_at": self.recorded_at,
            "valid_at": self.valid_at,
            "default_id": self.default_id,
            "default_version": self.default_version,
            "default_rationale": self.default_rationale,
            "derivation": self.derivation,
            "conflict": list(self.conflict),
            "affects": [i.to_json() for i in self.affects],
            "linked_results": list(self.linked_results),
            "limitations": list(self.limitations),
            "is_fact": self.is_fact,
            "used_in_this_report": self.is_used,
            "usage_note": (
                None if self.is_used else "Recorded; not used in this report"
            ),
        }


# --------------------------------------------------------------------------- #
# laboratory results
# --------------------------------------------------------------------------- #

#: The optional panels of §10.4 and where each belongs in the report. A panel
#: nobody supplied is one line saying so, not a page of empty charts.
LAB_PANELS: Final[tuple[dict[str, Any], ...]] = (
    {
        "panel_id": "scfa", "label": "Short-chain fatty acids",
        "analytes": ("acetate", "propionate", "butyrate", "valerate", "total_scfa"),
        "destination": "report.functions_detail",
        "relation": (
            "A measured concentration and a DNA capacity are different quantities. Both are "
            "shown, labelled separately; neither is converted into the other."
        ),
    },
    {
        "panel_id": "inflammation", "label": "Inflammation and immune markers",
        "analytes": ("calprotectin", "lactoferrin", "lysozyme", "secretory_iga"),
        "destination": "report.patterns_detail",
        "relation": (
            "Shown beside the barrier and inflammation context. A raised marker is a "
            "measurement, never an automatic diagnosis."
        ),
    },
    {
        "panel_id": "digestion", "label": "Digestion",
        "analytes": ("fecal_fat", "pancreatic_elastase", "reducing_carbohydrates"),
        "destination": "report.functions_detail",
        "relation": (
            "A direct assay of digestion, which is a different question from what the "
            "microbes can do; the two are kept apart."
        ),
    },
    {
        "panel_id": "stool_chemistry", "label": "Other stool chemistry",
        "analytes": ("occult_blood", "stool_ph", "beta_glucuronidase_activity"),
        "destination": "report.functions_detail",
        "relation": (
            "Measured activity is compared with genetic capacity for the same enzyme, with "
            "the assay's own meaning preserved."
        ),
    },
    {
        "panel_id": "metabolomics", "label": "Targeted metabolomics",
        "analytes": ("bile_acids", "urolithin", "equol", "ethanol", "acetaldehyde"),
        "destination": "report.functions_detail",
        "relation": "Shown on the corresponding metabolic card with its own units and method.",
    },
    {
        "panel_id": "endocrine", "label": "Endocrine and nutrient context",
        "analytes": ("glp1", "tsh", "vitamin_b12", "folate", "vitamin_d"),
        "destination": "report.patterns_detail",
        "relation": (
            "A microbial gene is not a hormone or a nutrient concentration. A supplied "
            "measurement is shown beside the capability, not merged with it."
        ),
    },
    {
        "panel_id": "microbial_load", "label": "Absolute microbial load",
        "analytes": ("total_load_qpcr", "total_load_flow_cytometry"),
        "destination": "report.community",
        "relation": (
            "The one measurement that turns relative composition into absolute amounts. "
            "Without it, every percentage here remains a share."
        ),
    },
    {
        "panel_id": "molecular_pathogen", "label": "Molecular pathogen assays",
        "analytes": ("pcr_target",),
        "destination": "report.pathogens",
        "relation": (
            "Attributed to the external assay and shown beside this report's own coverage, "
            "never merged into it."
        ),
    },
)
LAB_PANEL_BY_ID: Final[Mapping[str, Mapping[str, Any]]] = {
    p["panel_id"]: p for p in LAB_PANELS
}

QUALITATIVE_RESULTS: Final[frozenset[str]] = frozenset({
    "detected", "not_detected", "indeterminate",
})
CENSORING: Final[frozenset[str]] = frozenset({"<", ">", "="})


@dataclass(frozen=True)
class LabResult:
    """One imported laboratory measurement, with everything needed to trust it.

    Censoring is kept rather than resolved: a result reported as `<50` is not
    the number 50 and is not the number 0, and a register that stores it as
    either has lost the measurement.
    """

    analyte_id: str
    panel_id: str
    original_name: str
    specimen: str
    method: str | None = None
    value: float | None = None
    censoring: str = "="
    qualitative: str | None = None
    original_unit: str | None = None
    normalized_unit: str | None = None
    collected_at: str | None = None
    reference_low: float | None = None
    reference_high: float | None = None
    laboratory: str | None = None
    document_hash: str | None = None
    # molecular assays only
    target: str | None = None
    nucleic_acid: str | None = None
    ct: float | None = None
    detection_limit: float | None = None

    def __post_init__(self) -> None:
        if self.panel_id not in LAB_PANEL_BY_ID:
            raise InputRegisterError(f"{self.analyte_id}: unknown lab panel {self.panel_id!r}")
        if self.censoring not in CENSORING:
            raise InputRegisterError(f"{self.analyte_id}: unknown censoring {self.censoring!r}")
        if self.qualitative is not None and self.qualitative not in QUALITATIVE_RESULTS:
            raise InputRegisterError(
                f"{self.analyte_id}: unknown qualitative result {self.qualitative!r}"
            )
        if self.value is None and self.qualitative is None:
            raise InputRegisterError(
                f"{self.analyte_id}: a result needs either a value or a qualitative outcome"
            )
        if self.value is not None:
            if not math.isfinite(float(self.value)):
                raise InputRegisterError(f"{self.analyte_id}: value must be finite")
            if not self.original_unit:
                raise InputRegisterError(f"{self.analyte_id}: a value needs its original unit")
        if not self.specimen:
            raise InputRegisterError(
                f"{self.analyte_id}: a result must name its specimen; a vaginal or blood "
                "measurement cannot be silently read as a stool one"
            )

    @property
    def censored(self) -> bool:
        return self.censoring != "="

    @property
    def display_value(self) -> str:
        if self.qualitative:
            return self.qualitative.replace("_", " ")
        prefix = "" if self.censoring == "=" else self.censoring
        return f"{prefix}{self.value:g} {self.original_unit or ''}".strip()

    @property
    def in_reference_interval(self) -> bool | None:
        """Against the laboratory's own interval, or None when it cannot settle it.

        The laboratory's interval, not one this report invented, and only
        where the censoring actually decides the question. A result of
        `<50` against an interval of 10-200 could be 5 or 30: below the
        interval or inside it. Answering "within" there would be a guess
        presented as a classification, so the answer is None.
        """
        low, high = self.reference_low, self.reference_high
        if self.value is None or (low is None and high is None):
            return None
        value = float(self.value)
        if self.censoring == "<":
            # The true value lies somewhere below `value`.
            if low is not None and value <= low:
                return False        # even the ceiling is at or under the floor
            if low is not None:
                return None         # could be under the floor, could be inside
            return True if high is None or value <= high else None
        if self.censoring == ">":
            # The true value lies somewhere above `value`.
            if high is not None and value >= high:
                return False        # even the floor is at or over the ceiling
            if high is not None:
                return None
            return True if low is None or value >= low else None
        return (low is None or value >= low) and (high is None or value <= high)

    def to_json(self) -> dict[str, Any]:
        panel = LAB_PANEL_BY_ID[self.panel_id]
        return {
            "analyte_id": self.analyte_id,
            "panel_id": self.panel_id,
            "panel_label": panel["label"],
            "original_name": self.original_name,
            "specimen": self.specimen,
            "method": self.method,
            "value": self.value,
            "censoring": self.censoring,
            "censored": self.censored,
            "qualitative_result": self.qualitative,
            "display_value": self.display_value,
            "original_unit": self.original_unit,
            "normalized_unit": self.normalized_unit,
            "collected_at": self.collected_at,
            "reference_interval": {"low": self.reference_low, "high": self.reference_high},
            "within_laboratory_interval": self.in_reference_interval,
            "laboratory": self.laboratory,
            "document_hash": self.document_hash,
            "target": self.target,
            "nucleic_acid": self.nucleic_acid,
            "ct": self.ct,
            "detection_limit": self.detection_limit,
            "destination": panel["destination"],
            "relation_to_this_report": panel["relation"],
        }


def document_hash(payload: bytes | str) -> str:
    """A stable hash of the source document, so an import can be traced back."""
    raw = payload.encode("utf-8") if isinstance(payload, str) else payload
    return "sha256:" + hashlib.sha256(raw).hexdigest()


_IDENTIFIER = re.compile(r"[^a-z0-9]+")


def _slug(text: str) -> str:
    return _IDENTIFIER.sub("_", str(text).lower()).strip("_")


def link_is_allowed(
    result: LabResult, *, sample_specimen: str = "stool",
    sample_collected_at: str | None = None, max_days: int = 90,
) -> tuple[bool, str]:
    """Whether a laboratory result may be shown beside this sample's findings.

    Specimen first: a vaginal or blood measurement is not a stool
    measurement and cannot be quietly treated as one. Then date: a result
    from a year ago describes a different gut.
    """
    if _slug(result.specimen) != _slug(sample_specimen):
        return False, (
            f"this result is from {result.specimen}, and this report analysed "
            f"{sample_specimen}; it is recorded but not compared"
        )
    if result.collected_at and sample_collected_at:
        try:
            from datetime import date  # noqa: PLC0415

            a = date.fromisoformat(str(result.collected_at)[:10])
            b = date.fromisoformat(str(sample_collected_at)[:10])
        except ValueError:
            return True, "dates could not be compared; the result is shown with its own date"
        gap = abs((a - b).days)
        if gap > max_days:
            return False, (
                f"this result was collected {gap} days from the stool sample, beyond the "
                f"{max_days}-day window for a comparable measurement; it is recorded and "
                "labelled historical"
            )
    return True, "same specimen and a comparable collection date"


# --------------------------------------------------------------------------- #
# the register
# --------------------------------------------------------------------------- #


@dataclass
class InputRegister:
    """Every input this report had, in one place."""

    sample_id: str
    records: list[InputRecord] = field(default_factory=list)
    lab_results: list[LabResult] = field(default_factory=list)

    def add(self, record: InputRecord) -> InputRecord:
        self.records.append(record)
        return record

    def by_category(self, category: str) -> list[InputRecord]:
        return [r for r in self.records if r.category == category]

    @property
    def supplied(self) -> list[InputRecord]:
        return [r for r in self.records if r.input_state == "supplied"]

    @property
    def assumed(self) -> list[InputRecord]:
        return [r for r in self.records
                if r.input_state in {"assumed_default", "inherited_default"}]

    @property
    def missing(self) -> list[InputRecord]:
        return [r for r in self.records if r.input_state == "not_supplied"]

    @property
    def conflicting(self) -> list[InputRecord]:
        return [r for r in self.records if r.input_state == "conflicting"]

    @property
    def influential_assumptions(self) -> list[InputRecord]:
        """Assumptions that actually change something. §10.1: every one is visible."""
        return [r for r in self.assumed if r.is_used]

    def summary(self) -> dict[str, Any]:
        return {
            "n_supplied": len(self.supplied),
            "n_assumed": len(self.assumed),
            "n_influential_assumptions": len(self.influential_assumptions),
            "n_not_supplied": len(self.missing),
            "n_conflicting": len(self.conflicting),
            "n_laboratory_results": len(self.lab_results),
            "external_laboratory_results": (
                "Not supplied" if not self.lab_results else f"{len(self.lab_results)} imported"
            ),
        }

    def headline(self) -> str:
        """One sentence a reader can act on."""
        n_supplied, n_assumed = len(self.supplied), len(self.influential_assumptions)
        if not n_supplied and not n_assumed and not self.lab_results:
            return (
                "This report was generated from your sequencing data alone. Nothing was "
                "supplied about you, and nothing was assumed."
            )
        parts = []
        if not n_supplied and not n_assumed:
            parts.append("your sequencing data")
        if n_supplied:
            parts.append(f"{n_supplied} fact{'s' if n_supplied != 1 else ''} you supplied")
        if n_assumed:
            parts.append(
                f"{n_assumed} assumption{'s' if n_assumed != 1 else ''} the report had to make "
                "because they were not"
            )
        tail = (
            " No external laboratory results were supplied, and none are needed for anything "
            "above." if not self.lab_results else
            f" {len(self.lab_results)} external laboratory result"
            f"{'s were' if len(self.lab_results) != 1 else ' was'} imported."
        )
        return f"This report used {' and '.join(parts)}.{tail}"

    def missing_categories(self) -> list[dict[str, Any]]:
        """Optional inputs grouped into concise categories, per §10.1."""
        out: list[dict[str, Any]] = []
        for category, label in CATEGORIES:
            absent = [r for r in self.by_category(category) if r.input_state == "not_supplied"]
            if absent:
                out.append({
                    "category": category, "category_label": label,
                    "n_not_supplied": len(absent),
                    "examples": [r.label for r in absent[:4]],
                    "what_it_would_add": absent[0].limitations[0] if absent[0].limitations else None,
                })
        return out

    def lab_panel_status(self) -> list[dict[str, Any]]:
        """Every supported panel and whether anything was supplied for it."""
        supplied_by_panel: dict[str, list[LabResult]] = {}
        for result in self.lab_results:
            supplied_by_panel.setdefault(result.panel_id, []).append(result)
        return [
            {
                "panel_id": panel["panel_id"], "label": panel["label"],
                "supplied": panel["panel_id"] in supplied_by_panel,
                "n_results": len(supplied_by_panel.get(panel["panel_id"], [])),
                "results": [r.to_json() for r in supplied_by_panel.get(panel["panel_id"], [])],
                "destination": panel["destination"],
                "relation": panel["relation"],
                "analytes_recognised": list(panel["analytes"]),
            }
            for panel in LAB_PANELS
        ]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "feature_id": "A14",
            "sample_id": self.sample_id,
            "headline": self.headline(),
            "summary": self.summary(),
            "records": [r.to_json() for r in self.records],
            "by_category": [
                {
                    "category": category, "category_label": label,
                    "records": [r.to_json() for r in self.by_category(category)],
                }
                for category, label in CATEGORIES if self.by_category(category)
            ],
            "missing_categories": self.missing_categories(),
            "laboratory_panels": self.lab_panel_status(),
            "contract": [
                "An assumption is recorded as an assumption. It never becomes a supplied fact.",
                "An absent laboratory result has no imputed normal value and produces no "
                "normal or abnormal classification.",
                "A stored value with no consumer is labelled 'Recorded; not used in this "
                "report' rather than implied to be integrated.",
                "Conflicting supplied records are shown as conflicting; neither is silently "
                "overwritten.",
            ],
        }


__all__ = [
    "CATEGORIES",
    "CATEGORY_LABELS",
    "CENSORING",
    "FACTUAL_STATES",
    "INFLUENCE_KINDS",
    "INPUT_STATES",
    "LAB_PANELS",
    "LAB_PANEL_BY_ID",
    "QUALITATIVE_RESULTS",
    "SCHEMA_VERSION",
    "Influence",
    "InputRecord",
    "InputRegister",
    "InputRegisterError",
    "LabResult",
    "build_register",
    "document_hash",
    "link_is_allowed",
    "records_from_lab_results",
]


# --------------------------------------------------------------------------- #
# building the register from a real run
# --------------------------------------------------------------------------- #


def _sample_records(results: Mapping[str, Any], sample_id: str) -> list[InputRecord]:
    """Participant, specimen and processing facts this run actually had."""
    out: list[InputRecord] = []
    meta = results.get("run") or results.get("meta") or {}
    out.append(InputRecord(
        input_id="specimen.sample_id", category="participant_specimen",
        label="Sample identifier", input_state="supplied",
        reported_value=sample_id, effective_value=sample_id, source="sample manifest",
        affects=(Influence("report.identity", "report header", "adds_context",
                           "Names this report and links it to your sequencing files."),),
    ))
    out.append(InputRecord(
        input_id="specimen.type", category="participant_specimen",
        label="Sample type", input_state="inherited_default",
        effective_value="stool", default_id="specimen.default_type",
        default_version="1.0",
        default_rationale=(
            "This pipeline analyses stool. No other specimen type is accepted, so the value "
            "is fixed rather than guessed."
        ),
        affects=(Influence("report.all", "specimen gate", "changes_interpretation",
                           "Every reading here describes stool. A measurement from another "
                           "body site cannot be read off this report."),),
    ))
    collected = meta.get("collection_date") or meta.get("collected_at")
    if collected:
        out.append(InputRecord(
            input_id="specimen.collected_at", category="participant_specimen",
            label="Collection date", input_state="supplied",
            reported_value=str(collected), effective_value=str(collected),
            source="sample manifest",
            affects=(Influence("A13.longitudinal", "comparability window", "changes_calculation",
                               "Orders this sample against your others and decides which "
                               "laboratory results are close enough to compare."),),
        ))
    else:
        out.append(InputRecord(
            input_id="specimen.collected_at", category="participant_specimen",
            label="Collection date", input_state="not_supplied",
            limitations=(
                "Without it, a repeat sample cannot be placed on a timeline and an imported "
                "laboratory result cannot be checked for date compatibility.",
            ),
        ))
    # §12.2: the full QC tables live with the sequencing-quality section;
    # the register carries a compact status and points there.
    gates = ((results.get("sequencing_quality") or {}).get("gates")
             or (results.get("quality") or {}).get("gates") or {})
    if gates:
        usable = gates.get("usable_nonhost_pairs")
        out.append(InputRecord(
            input_id="processing.quality_gates", category="collection_processing",
            label="Sequencing quality gates", input_state="derived",
            effective_value=str(gates.get("overall", "unknown")).replace("_", " "),
            derivation="evaluated by the pipeline's own quality gates on this run",
            source="pipeline",
            affects=(Influence(
                "report.sequencing", "quality gate", "changes_interpretation",
                "Decides whether the readings are fit to interpret at all. The full gate "
                "table, with every threshold and what was measured against it, is in the "
                "sequencing-quality section."
                + (f" {usable:,} usable non-host pairs." if isinstance(usable, int) else ""),
            ),),
        ))
    profiler = ((results.get("taxonomy") or {}).get("profiler")
                or (results.get("community_profile") or {}).get("method"))
    if profiler:
        out.append(InputRecord(
            input_id="processing.profiler", category="collection_processing",
            label="Taxonomic profiler", input_state="derived",
            reported_value=None, effective_value=str(profiler),
            derivation="read from the pipeline's own run record",
            source="pipeline",
            affects=(Influence("A09.community_type", "reference compatibility",
                               "changes_interpretation",
                               "Decides which reference cohorts this sample can be compared "
                               "with, and which genus names have to be restated first."),),
        ))
    return out


def _context_records(results: Mapping[str, Any]) -> list[InputRecord]:
    """The existing subject-context ledger, imported rather than redesigned.

    §10.3(2): the production application's defaults already exist and are
    already shown. This brings them into one register without changing what
    any of them means.
    """
    context = results.get("subject_context") or {}
    out: list[InputRecord] = []
    kind_to_category = {
        "medication": "medications_supplements",
        "exposure": "longitudinal_exposures",
        "demographic": "participant_specimen",
        "safety": "clinical_history",
    }
    for key, value in (context.get("supplied") or {}).items():
        out.append(InputRecord(
            input_id=f"context.{key}", category="clinical_history",
            label=str(key).replace("_", " "), input_state="supplied",
            reported_value=value, effective_value=value, source="subject context",
            affects=(Influence("actions.safety_rules", "context ledger",
                               "changes_interpretation",
                               "Read by the safety rules on the action cards."),),
        ))
    for entry in context.get("assumed") or []:
        kind = str(entry.get("kind") or "demographic")
        out.append(InputRecord(
            input_id=f"context.{entry.get('field')}", category=kind_to_category.get(kind, "clinical_history"),
            label=str(entry.get("label") or entry.get("field")),
            input_state="inherited_default",
            effective_value=entry.get("assumed_value"),
            default_id=f"context.{entry.get('field')}",
            default_version=str(context.get("mode") or "default"),
            default_rationale=(
                f"Not supplied. In {context.get('mode', 'this')} mode the report assumes "
                f"{entry.get('assumed_value')!r} so the reading can be produced at all."
            ),
            affects=(Influence(
                str(entry.get("used_by") or "report"), "subject context ledger",
                "changes_interpretation", str(entry.get("if_present") or ""),
            ),),
            limitations=("An assumption made to produce a reading, not a fact about you.",),
        ))
    for entry in context.get("unknown") or []:
        label = entry.get("label") if isinstance(entry, Mapping) else str(entry)
        detail = entry.get("if_present", "") if isinstance(entry, Mapping) else ""
        out.append(InputRecord(
            input_id=f"context.{_slug(str(label))}", category="clinical_history",
            label=str(label), input_state="not_supplied",
            limitations=(str(detail) or "No default applies, so this stays unknown.",),
        ))
    return out


def _analysis_records(results: Mapping[str, Any]) -> list[InputRecord]:
    """Assumptions the new capabilities made, each naming what it changed."""
    out: list[InputRecord] = []
    extension = results.get("extension") or {}
    views = extension.get("views") or {}

    simulation = views.get("simulation") or {}
    scenarios = simulation.get("scenarios") or {}
    if simulation.get("readiness"):
        media = {m.get("id"): m.get("label") for m in simulation.get("media") or []}
        baseline = scenarios.get("baseline_medium") or "european"
        out.append(InputRecord(
            input_id="diet.scenario_medium", category="diet_tolerance",
            label="Diet used for the synbiotic scenarios",
            input_state="assumed_default",
            effective_value=media.get(baseline, baseline),
            default_id="scenario.default_medium",
            default_version=str((simulation.get("protocol") or {}).get("protocol_id") or "1.0"),
            default_rationale=(
                "No record of what you actually eat was supplied, so the scenarios use the "
                "published reference diet the protocol declares, and are additionally run "
                "against every other declared diet so the range is visible."
            ),
            affects=(Influence(
                "A16.scenarios", "scenario medium", "changes_calculation",
                "Sets how much of each nutrient the modelled community has. The sensitivity "
                "run shows how much the answer moves when this assumption changes - here it "
                "moves a great deal, which is why the range is printed.",
            ),),
            limitations=("This describes a modelling assumption, not your actual diet.",),
        ))
        out.append(InputRecord(
            input_id="diet.reported_intake", category="diet_tolerance",
            label="Your actual dietary pattern", input_state="not_supplied",
            limitations=(
                "With it, the scenarios could be run on the diet you actually eat instead of "
                "a reference diet, and the fibre rankings would apply to you directly.",
            ),
        ))
        out.append(InputRecord(
            input_id="diet.intolerances", category="diet_tolerance",
            label="Known fibre intolerances or allergies", input_state="not_supplied",
            limitations=(
                "With it, a fibre you react to would be ranked accordingly and supported "
                "alternatives suggested; its measured genetic capacity would still be shown.",
            ),
        ))

    ecology = views.get("ecology") or {}
    model = ecology.get("community_type_model") or {}
    if model:
        cohort = model.get("cohort") or {}
        out.append(InputRecord(
            input_id="reference.community_type_cohort", category="analysis_configuration",
            label="Reference cohort for the community type",
            input_state="assumed_default",
            effective_value=(
                f"{cohort.get('n_reference_samples', 0):,} samples"
                if isinstance(cohort.get("n_reference_samples"), int) else "frozen cohort"
            ),
            default_id="community_type.frozen_model",
            default_version=str(model.get("model_id")),
            default_rationale=(
                f"No tighter matching information was supplied, so the documented eligible "
                f"reference population is used: {cohort.get('n_reference_samples', 'the')} "
                f"samples from {cohort.get('n_studies', 'several')} studies, frozen as "
                f"{model.get('model_id')}."
            ),
            affects=(Influence(
                "A09.community_type", "nearest-medoid assignment", "changes_calculation",
                "Every distance on the community-type card is measured against this frozen "
                "cohort, and never against the other samples analysed alongside yours.",
            ),),
        ))
    return out


def build_register(
    results: Mapping[str, Any], *, sample_id: str | None = None,
    lab_results: Sequence[LabResult] = (),
) -> InputRegister:
    """The whole register for one run."""
    sample = str(sample_id or results.get("sample") or "sample")
    register = InputRegister(sample_id=sample)
    for record in (
        *_sample_records(results, sample),
        *_context_records(results),
        *_analysis_records(results),
    ):
        record.sample_id = record.sample_id or sample
        register.add(record)
    register.lab_results = list(lab_results)
    if not lab_results:
        register.add(InputRecord(
            input_id="laboratory.any", category="laboratory",
            label="External laboratory results", input_state="not_supplied",
            limitations=(
                "None are needed: every reading in this report is produced from your "
                "sequencing data. A supplied result would be shown beside the finding it "
                "relates to, not used to replace it.",
            ),
        ))
    return register


def records_from_lab_results(results: Iterable[LabResult]) -> list[InputRecord]:
    """One register record per imported laboratory result."""
    out: list[InputRecord] = []
    for result in results:
        panel = LAB_PANEL_BY_ID[result.panel_id]
        out.append(InputRecord(
            input_id=f"laboratory.{result.panel_id}.{result.analyte_id}",
            category="laboratory", label=result.original_name,
            input_state="supplied", reported_value=result.display_value,
            effective_value=result.value if result.value is not None else result.qualitative,
            unit=result.original_unit, source=result.laboratory or "external laboratory",
            recorded_at=result.collected_at, valid_at=result.collected_at,
            affects=(Influence(
                str(panel["destination"]), f"laboratory import: {result.panel_id}",
                "adds_context", str(panel["relation"]),
            ),),
            limitations=(
                "A direct measurement of one analyte. It does not replace, confirm or "
                "invalidate the genetic capacity measured here; the two are different "
                "quantities and are shown separately.",
            ),
        ))
    return out
