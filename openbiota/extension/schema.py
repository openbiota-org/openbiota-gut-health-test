"""The extension metric contract from BUILD_SPEC_v0.8.3 section 3.3.

One record per new measurement, with a stable ID, one unit, one denominator,
and vocabularies that keep apart the things a report must never blur:

* **Zero, null and not-assessed are different states.** A gene present at zero
  abundance, a gene nobody looked for, and a gene the assay cannot see are
  three different facts about a sample, and each has its own `state`.
* **A percentile is not a percentage of anything.** `kind` separates a
  reference percentile from a pathway completeness, an abundance, a model
  prediction and an imported laboratory value.
* **High is not automatically good.** `direction` defaults to descriptive and
  has to name the context that justifies any other reading.
* **Co-detection is not carriage.** Every contributing taxon states how the
  assignment was made, so "this species was present and so was this gene"
  cannot be printed as "this species carries this gene".
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.extension import SPEC_VERSION

SCHEMA_VERSION: Final = "openbiota.extension-metric/1.0"

#: What kind of quantity this is. The report's formatting and every comparison
#: rule keys off this, so a model prediction can never be rendered in the same
#: place as a measured concentration.
KINDS: Final[frozenset[str]] = frozenset({
    "taxon_abundance",
    "gene_abundance",
    "genetic_capacity",
    "community_composition",
    "reference_percentile",
    "experimental_index",
    "model_prediction",
    "external_lab_measurement",
    "evidence_context",
})

#: Why a value is, or is not, there.
#:
#: ``measured``                            quantified from this sample
#: ``partial``                             some required steps supported
#: ``not_detected_above_assay_threshold``  looked for, nothing above the floor
#: ``insufficient_coverage``               looked for, depth cannot answer
#: ``not_assayed``                         not searched in this run
#: ``unsupported_by_assay``                this assay cannot answer, at any depth
#: ``requires_external_result``            needs a laboratory measurement
#: ``not_applicable``                      the question does not apply here
STATES: Final[frozenset[str]] = frozenset({
    "measured",
    "partial",
    "not_detected_above_assay_threshold",
    "insufficient_coverage",
    "not_assayed",
    "unsupported_by_assay",
    "requires_external_result",
    "not_applicable",
})

#: States that carry a number. Anything else must have `value is None`: a
#: not-assayed metric with a 0 in its value field is the single most
#: misleading record this schema could hold.
NUMERIC_STATES: Final[frozenset[str]] = frozenset({"measured", "partial"})

#: Which way is better, and only with a cited context.
DIRECTIONS: Final[frozenset[str]] = frozenset({
    "higher_favourable_in_context",
    "lower_favourable_in_context",
    "interval_favourable_in_context",
    "context_dependent",
    "descriptive",
})

#: A direction other than these two has to say in which context it holds.
DIRECTIONS_REQUIRING_CONTEXT: Final[frozenset[str]] = frozenset(
    DIRECTIONS - {"descriptive", "context_dependent"}
)

#: How strongly the sequence evidence supports the reading.
CONFIDENCE: Final[frozenset[str]] = frozenset({
    "supported", "provisional", "weak", "unresolved", "not_applicable",
})

#: How well established the biology is, independently of this sample's data.
EVIDENCE_MATURITY: Final[frozenset[str]] = frozenset({
    "characterized_enzyme",
    "biochemical_route",
    "genomic_prediction",
    "association_only",
    "mechanistic_hypothesis",
    "not_applicable",
})

#: Whether the steps of a route were found in one genome or across the
#: community. Cross-feeding between organisms is a real biological
#: possibility; it is not the same claim as one organism holding the pathway.
COMPLETENESS_SCOPES: Final[frozenset[str]] = frozenset({
    "genome_linked", "community_assembled", "reference_inferred", "not_applicable",
})

#: How a contributing taxon was tied to the observation. Section 3.3: "Co-
#: detection of a species and a gene is not sufficient to claim that species
#: carries that gene in the sample."
ASSIGNMENT_EVIDENCE: Final[frozenset[str]] = frozenset({
    "direct_strain_specific_sequence",
    "uniquely_mapped_gene",
    "supported_contig_linkage",
    "classified_mag",
    "ambiguous_reference_hit",
    "organism_level_association",
})

#: Assignment levels that support saying this organism carries this function
#: in this sample. The rest are context.
CARRIER_EVIDENCE: Final[frozenset[str]] = frozenset({
    "direct_strain_specific_sequence",
    "uniquely_mapped_gene",
    "supported_contig_linkage",
})

#: The state of a reference comparison, kept apart from the value itself so a
#: raw capacity stays visible when no eligible cohort exists.
REFERENCE_STATES: Final[frozenset[str]] = frozenset({
    "available", "not_available", "incompatible_method", "insufficient_reference_samples",
    "prevalence_only",
})


class ExtensionSchemaError(ValueError):
    """An extension record that would misrepresent its own evidence."""


@dataclass(slots=True)
class Contributor:
    """One organism behind a measurement, with the evidence of assignment."""

    taxon: str
    assignment_evidence: str
    value: float | None = None
    unit: str | None = None
    share_of_metric: float | None = None
    taxon_id: str | None = None
    lane: str | None = None
    note: str | None = None

    def __post_init__(self) -> None:
        if self.assignment_evidence not in ASSIGNMENT_EVIDENCE:
            raise ExtensionSchemaError(
                f"{self.taxon}: unknown assignment evidence {self.assignment_evidence!r}"
            )

    @property
    def is_carrier(self) -> bool:
        """Whether this row supports "this organism carries this function"."""
        return self.assignment_evidence in CARRIER_EVIDENCE

    def to_json(self) -> dict[str, Any]:
        return {
            "taxon": self.taxon,
            "taxon_id": self.taxon_id,
            "assignment_evidence": self.assignment_evidence,
            "is_carrier_evidence": self.is_carrier,
            "value": self.value,
            "unit": self.unit,
            "share_of_metric": self.share_of_metric,
            "lane": self.lane,
            "note": self.note,
        }


@dataclass(slots=True)
class Component:
    """One step or subunit inside a route, with its own support.

    Published separately from the route's headline number so an isolated
    marker cannot masquerade as a complete pathway (section 4.1).
    """

    component_id: str
    label: str
    state: str
    value: float | None = None
    unit: str | None = None
    required: bool = True
    alternative_group: str | None = None
    fragments: int | None = None
    covered_residues: int | None = None
    identity: float | None = None
    reference_ids: tuple[str, ...] = ()
    note: str | None = None

    def __post_init__(self) -> None:
        if self.state not in STATES:
            raise ExtensionSchemaError(f"{self.component_id}: unknown state {self.state!r}")
        if self.state not in NUMERIC_STATES and self.value is not None:
            raise ExtensionSchemaError(
                f"{self.component_id}: state {self.state!r} cannot carry a value"
            )

    @property
    def supported(self) -> bool:
        return self.state in NUMERIC_STATES

    def to_json(self) -> dict[str, Any]:
        return {
            "component_id": self.component_id,
            "label": self.label,
            "state": self.state,
            "value": self.value,
            "unit": self.unit,
            "required": self.required,
            "alternative_group": self.alternative_group,
            "fragments": self.fragments,
            "covered_residues": self.covered_residues,
            "identity": self.identity,
            "reference_ids": list(self.reference_ids),
            "note": self.note,
        }


@dataclass(slots=True)
class ExtensionMetric:
    """One new measurement, in the form section 3.3 fixes.

    Fields may be added by a future release; none may be reused for a second
    meaning. The validation in `__post_init__` refuses the specific
    misstatements the specification calls out, rather than checking types.
    """

    metric_id: str
    label: str
    kind: str
    state: str
    method_id: str
    value: float | None = None
    unit: str | None = None
    denominator: str | None = None
    direction: str = "descriptive"
    direction_context: str | None = None
    reference_percentile: float | None = None
    reference_id: str | None = None
    reference_state: str = "not_available"
    prevalence_in_reference: float | None = None
    positive_carrier_percentile: float | None = None
    score_0_100: float | None = None
    score_definition_id: str | None = None
    analytical_confidence: str = "unresolved"
    evidence_maturity: str = "not_applicable"
    pathway_completeness: float | None = None
    completeness_scope: str = "not_applicable"
    linkage_evidence: tuple[str, ...] = ()
    assessable_fraction: float | None = None
    components: tuple[Component, ...] = ()
    contributing_taxa: tuple[Contributor, ...] = ()
    source_ids: tuple[str, ...] = ()
    input_fingerprint: str | None = None
    limitations: tuple[str, ...] = ()
    action_links: tuple[str, ...] = ()
    legacy_metric_links: tuple[str, ...] = ()
    group: str | None = None
    feature_id: str | None = None
    view_ids: tuple[str, ...] = ()
    extra: dict[str, Any] = field(default_factory=dict)
    schema_version: str = SCHEMA_VERSION
    spec_version: str = SPEC_VERSION

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ExtensionSchemaError(f"{self.metric_id}: unknown kind {self.kind!r}")
        if self.state not in STATES:
            raise ExtensionSchemaError(f"{self.metric_id}: unknown state {self.state!r}")
        if self.direction not in DIRECTIONS:
            raise ExtensionSchemaError(f"{self.metric_id}: unknown direction {self.direction!r}")
        if self.analytical_confidence not in CONFIDENCE:
            raise ExtensionSchemaError(
                f"{self.metric_id}: unknown confidence {self.analytical_confidence!r}"
            )
        if self.evidence_maturity not in EVIDENCE_MATURITY:
            raise ExtensionSchemaError(
                f"{self.metric_id}: unknown evidence maturity {self.evidence_maturity!r}"
            )
        if self.completeness_scope not in COMPLETENESS_SCOPES:
            raise ExtensionSchemaError(
                f"{self.metric_id}: unknown completeness scope {self.completeness_scope!r}"
            )
        if self.reference_state not in REFERENCE_STATES:
            raise ExtensionSchemaError(
                f"{self.metric_id}: unknown reference state {self.reference_state!r}"
            )

        # Zero is not the same as unmeasured.
        if self.state not in NUMERIC_STATES and self.value is not None:
            raise ExtensionSchemaError(
                f"{self.metric_id}: state {self.state!r} cannot carry a value "
                f"({self.value!r}); an unmeasured metric with a number in it reads as a zero result"
            )
        if self.state in NUMERIC_STATES and self.value is None and self.kind != "evidence_context":
            raise ExtensionSchemaError(
                f"{self.metric_id}: state {self.state!r} promises a value and has none"
            )
        if self.value is not None:
            if not isinstance(self.value, (int, float)) or not math.isfinite(float(self.value)):
                raise ExtensionSchemaError(f"{self.metric_id}: value must be finite")
            if not self.unit:
                raise ExtensionSchemaError(f"{self.metric_id}: a value needs a unit")
            if not self.denominator and self.kind in {
                "gene_abundance", "taxon_abundance", "genetic_capacity", "community_composition",
            }:
                raise ExtensionSchemaError(
                    f"{self.metric_id}: an abundance needs its denominator named"
                )

        # A percentile without a reference is a fabrication.
        if self.reference_percentile is not None:
            if self.reference_state != "available":
                raise ExtensionSchemaError(
                    f"{self.metric_id}: a percentile requires reference_state='available'"
                )
            if not self.reference_id:
                raise ExtensionSchemaError(f"{self.metric_id}: a percentile needs its reference ID")
            if not 0.0 <= float(self.reference_percentile) <= 100.0:
                raise ExtensionSchemaError(f"{self.metric_id}: percentile out of range")
        if self.reference_state == "available" and not self.reference_id:
            raise ExtensionSchemaError(f"{self.metric_id}: reference_state 'available' needs an ID")

        # A direction that claims a polarity has to name the context for it.
        if self.direction in DIRECTIONS_REQUIRING_CONTEXT and not self.direction_context:
            raise ExtensionSchemaError(
                f"{self.metric_id}: direction {self.direction!r} must cite its context"
            )

        if self.score_0_100 is not None:
            if not self.score_definition_id:
                raise ExtensionSchemaError(
                    f"{self.metric_id}: a 0-100 score must name its definition"
                )
            if not 0.0 <= float(self.score_0_100) <= 100.0:
                raise ExtensionSchemaError(f"{self.metric_id}: score out of range")

        for name, frac in (
            ("pathway_completeness", self.pathway_completeness),
            ("assessable_fraction", self.assessable_fraction),
            ("prevalence_in_reference", self.prevalence_in_reference),
        ):
            if frac is not None and not 0.0 <= float(frac) <= 1.0:
                raise ExtensionSchemaError(f"{self.metric_id}: {name} must be a fraction 0-1")

        # Genome linkage is a claim about evidence, not a default.
        if self.completeness_scope == "genome_linked" and not self.linkage_evidence:
            raise ExtensionSchemaError(
                f"{self.metric_id}: 'genome_linked' completeness requires linkage evidence"
            )

        # An externally measured analyte cannot be filled from sequence.
        if (
            self.kind == "external_lab_measurement"
            and self.state in NUMERIC_STATES
            and not self.extra.get("lab_record_id")
        ):
            raise ExtensionSchemaError(
                f"{self.metric_id}: a measured laboratory value needs its source record"
            )

        if self.state in NUMERIC_STATES and self.kind != "evidence_context" and not self.limitations:
            raise ExtensionSchemaError(
                f"{self.metric_id}: every measured metric states what it does not measure"
            )

    # -- derived views ----------------------------------------------------- #

    @property
    def carriers(self) -> tuple[Contributor, ...]:
        """Contributors whose evidence supports carriage, not co-detection."""
        return tuple(c for c in self.contributing_taxa if c.is_carrier)

    @property
    def required_components(self) -> tuple[Component, ...]:
        return tuple(c for c in self.components if c.required)

    def completeness_from_components(self) -> float | None:
        """Fraction of required steps supported, honouring alternatives.

        An OR step is satisfied by any one of its alternatives; an AND step
        needs its own support. Returns None when the route has no declared
        components, because an unstated route has no completeness rather than
        a completeness of zero.
        """
        required = self.required_components
        if not required:
            return None
        groups: dict[str, list[Component]] = {}
        for comp in required:
            key = comp.alternative_group or f"__and__{comp.component_id}"
            groups.setdefault(key, []).append(comp)
        satisfied = sum(1 for members in groups.values() if any(m.supported for m in members))
        return satisfied / len(groups)

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "schema_version": self.schema_version,
            "spec_version": self.spec_version,
            "metric_id": self.metric_id,
            "label": self.label,
            "group": self.group,
            "feature_id": self.feature_id,
            "view_ids": list(self.view_ids),
            "kind": self.kind,
            "state": self.state,
            "value": self.value,
            "unit": self.unit,
            "denominator": self.denominator,
            "direction": self.direction,
            "direction_context": self.direction_context,
            "reference_percentile": self.reference_percentile,
            "reference_id": self.reference_id,
            "reference_state": self.reference_state,
            "prevalence_in_reference": self.prevalence_in_reference,
            "positive_carrier_percentile": self.positive_carrier_percentile,
            "score_0_100": self.score_0_100,
            "score_definition_id": self.score_definition_id,
            "analytical_confidence": self.analytical_confidence,
            "evidence_maturity": self.evidence_maturity,
            "pathway_completeness": self.pathway_completeness,
            "completeness_scope": self.completeness_scope,
            "linkage_evidence": list(self.linkage_evidence),
            "assessable_fraction": self.assessable_fraction,
            "components": [c.to_json() for c in self.components],
            "contributing_taxa": [c.to_json() for c in self.contributing_taxa],
            "n_carrier_supported_taxa": len(self.carriers),
            "source_ids": list(self.source_ids),
            "method_id": self.method_id,
            "input_fingerprint": self.input_fingerprint,
            "limitations": list(self.limitations),
            "action_links": list(self.action_links),
            "legacy_metric_links": list(self.legacy_metric_links),
        }
        if self.extra:
            out["detail"] = self.extra
        return out


def fingerprint(*parts: Any) -> str:
    """A stable content fingerprint for an extension record's inputs.

    Used for the dependency manifests of section 3.2: cache reuse requires
    every relevant input, database release, parameter and program version to
    match, and that is only checkable if the identity is derived from the
    content rather than from a timestamp.
    """
    payload = json.dumps(parts, sort_keys=True, default=str, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def midrank_percentile(value: float, reference: Sequence[float]) -> float | None:
    """The documented midrank empirical percentile of section 4.2.

    ``P(x; R) = 100 * (count(R < x) + 0.5 * count(R == x)) / len(R)``

    Ties take half their own mass, so a value sitting exactly on a cluster of
    reference samples lands in the middle of that cluster rather than above or
    below all of it. Returns None for an empty reference, because a percentile
    against nothing is the fabrication section 4.2 forbids.
    """
    if not reference:
        return None
    below = sum(1 for r in reference if r < value)
    equal = sum(1 for r in reference if r == value)
    return 100.0 * (below + 0.5 * equal) / len(reference)


def bottleneck(values: Sequence[float | None]) -> float | None:
    """The minimum across required steps: an engineering bottleneck proxy.

    Section 4.1 is explicit that this is not measured flux or a production
    rate. Returns None when any required step is unsupported, because the
    minimum of a set containing an unknown is unknown, not zero.
    """
    if not values:
        return None
    if any(v is None for v in values):
        return None
    return min(float(v) for v in values if v is not None)


def fragment_rpkm(fragments: int, effective_length_nt: int, total_fragments: int) -> float | None:
    """``1e9 * n_g / (L_g * N)`` from section 4.1.

    The effective length is in nucleotides. A protein reference must convert
    first: using an amino-acid length here would silently inflate every
    protein-referenced abundance threefold.
    """
    if effective_length_nt <= 0 or total_fragments <= 0:
        return None
    return 1e9 * fragments / (effective_length_nt * total_fragments)


def cds_equivalent_length(amino_acids: int, *, include_stop: bool = True) -> int:
    """Nucleotide-equivalent length of a protein reference.

    Section 4.1 requires the stop-codon treatment to be fixed rather than
    left to each call site, so it is fixed here and recorded in the unit.
    """
    if amino_acids <= 0:
        return 0
    return 3 * amino_acids + (3 if include_stop else 0)


def validate_aggregate(
    total: float,
    parts: Mapping[str, float],
    *,
    shared: Mapping[str, float] | None = None,
    rtol: float = 1e-9,
) -> None:
    """Refuse an aggregate that counts one fragment twice.

    A gene shared by two substrate panels appears as linked evidence in both
    (section 5.1), but its fragments may enter an aggregate once. Shared
    contributions are passed separately and subtracted, so the check is on the
    arithmetic rather than on a promise.
    """
    summed = sum(parts.values()) - sum((shared or {}).values())
    if not math.isclose(summed, total, rel_tol=rtol, abs_tol=1e-12):
        raise ExtensionSchemaError(
            f"aggregate {total!r} does not reconcile with its parts {summed!r}: "
            "a shared gene has probably been counted twice"
        )


__all__ = [
    "ASSIGNMENT_EVIDENCE",
    "CARRIER_EVIDENCE",
    "COMPLETENESS_SCOPES",
    "CONFIDENCE",
    "DIRECTIONS",
    "DIRECTIONS_REQUIRING_CONTEXT",
    "EVIDENCE_MATURITY",
    "KINDS",
    "NUMERIC_STATES",
    "REFERENCE_STATES",
    "SCHEMA_VERSION",
    "STATES",
    "Component",
    "Contributor",
    "ExtensionMetric",
    "ExtensionSchemaError",
    "bottleneck",
    "cds_equivalent_length",
    "fingerprint",
    "fragment_rpkm",
    "midrank_percentile",
    "validate_aggregate",
]
