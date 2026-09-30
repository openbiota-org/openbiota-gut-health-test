"""The shared resolution evidence record and its status axes.

Three axes, deliberately independent (spec §3.2), because collapsing them is
how a report comes to state something it never measured:

* **assay_status** — did we run the assay at all?
* **analytical_call** — what did the sequence support?
* **biological_evidence** — what is known about the thing, from the
  literature, as a list rather than a single grade.

Detection accuracy, causal evidence and clinical prediction are different
quantities and must not be averaged into one confidence number. A target whose
assay never ran has no analytical call; a target with a supported detection may
still carry no evidence beyond identity.

The invariants below are enforced in code rather than left to reviewer
discipline, because every one of them is a mistake this report has made or
could make:

* An assay that did not complete cannot carry a negative call (§3.2).
* A null is never rendered as a zero (§10, acceptance 70).
* `disease_transmission_probability` exists only to be null; no computation in
  this codebase is permitted to populate it (§9.1).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from . import RESOLUTION_SCHEMA_VERSION

# --------------------------------------------------------------------------- #
# Axis 1: did the assay run?
# --------------------------------------------------------------------------- #

NOT_REQUESTED: Final = "not_requested"
SCHEDULED: Final = "scheduled"
COMPLETED: Final = "completed"
FAILED: Final = "failed"
REFERENCE_UNAVAILABLE: Final = "reference_unavailable"
ACCESS_UNAVAILABLE: Final = "access_or_license_unavailable"
INCOMPATIBLE_INPUT: Final = "incompatible_input"

ASSAY_STATUSES: Final[frozenset[str]] = frozenset({
    NOT_REQUESTED, SCHEDULED, COMPLETED, FAILED,
    REFERENCE_UNAVAILABLE, ACCESS_UNAVAILABLE, INCOMPATIBLE_INPUT,
})

#: Only a completed assay may carry a call about what is or is not present.
#: Everything else leaves `analytical_call` at `unresolved`.
CONCLUSIVE_ASSAY: Final[frozenset[str]] = frozenset({COMPLETED})

# --------------------------------------------------------------------------- #
# Axis 2: what did the sequence support?
# --------------------------------------------------------------------------- #

UNRESOLVED: Final = "unresolved"
CANDIDATE_SEQUENCE: Final = "candidate_sequence_detected"
SUPPORTED_DETECTION: Final = "supported_detection"
NOT_DETECTED_CALIBRATED: Final = "not_detected_with_calibrated_limit"
NOT_DETECTED_UNVALIDATED: Final = "not_detected_limit_unvalidated"
INSUFFICIENT_DEPTH: Final = "insufficient_informative_depth"
AMBIGUOUS_HOMOLOGY: Final = "ambiguous_homology"
MIXED_UNRESOLVED: Final = "mixed_components_unresolved"
PHENOTYPE_NOT_MEASURED: Final = "phenotype_not_measured"

ANALYTICAL_CALLS: Final[frozenset[str]] = frozenset({
    UNRESOLVED, CANDIDATE_SEQUENCE, SUPPORTED_DETECTION,
    NOT_DETECTED_CALIBRATED, NOT_DETECTED_UNVALIDATED,
    INSUFFICIENT_DEPTH, AMBIGUOUS_HOMOLOGY, MIXED_UNRESOLVED,
    PHENOTYPE_NOT_MEASURED,
})

#: Calls that assert absence. Reachable only from a completed assay, and the
#: calibrated form additionally requires a recorded detection limit.
NEGATIVE_CALLS: Final[frozenset[str]] = frozenset({
    NOT_DETECTED_CALIBRATED, NOT_DETECTED_UNVALIDATED,
})

#: Calls that assert some presence.
POSITIVE_CALLS: Final[frozenset[str]] = frozenset({
    CANDIDATE_SEQUENCE, SUPPORTED_DETECTION,
})

# --------------------------------------------------------------------------- #
# Axis 3: what is known about it, in the literature?
# --------------------------------------------------------------------------- #

IDENTITY_ONLY: Final = "identity_only"
SEQUENCE_PREDICTED_FUNCTION: Final = "sequence_predicted_function"
CHARACTERIZED_FUNCTION: Final = "experimentally_characterized_function"
HUMAN_ASSOCIATION: Final = "human_association"
ANIMAL_MECHANISTIC: Final = "animal_mechanistic_effect"
HUMAN_INTERVENTION: Final = "human_intervention_evidence"
CLINICALLY_VALIDATED: Final = "clinically_validated_for_specific_endpoint"

BIOLOGICAL_EVIDENCE_KINDS: Final[frozenset[str]] = frozenset({
    IDENTITY_ONLY, SEQUENCE_PREDICTED_FUNCTION, CHARACTERIZED_FUNCTION,
    HUMAN_ASSOCIATION, ANIMAL_MECHANISTIC, HUMAN_INTERVENTION,
    CLINICALLY_VALIDATED,
})

# --------------------------------------------------------------------------- #
# What kind of identity a target is about
# --------------------------------------------------------------------------- #

IDENTITY_KINDS: Final[frozenset[str]] = frozenset({
    "lineage", "species", "species_complex", "sgb", "subspecies",
    "population_fingerprint", "strain_component", "named_reference_isolate",
    "sequence_type", "serotype", "pathotype", "biotype",
    "accessory_region", "locus", "allele", "operon", "plasmid_replicon",
    "gene_family", "phenotype",
})

#: Carrier-linkage strength, weakest to strongest (spec §6.1). These describe
#: evidence; they are not scores and must not be arithmetically combined.
LINKAGE_STATES: Final[tuple[str, ...]] = (
    "unassigned",
    "sample_co_detection",
    "read_pair_or_amplicon_linkage",
    "assembled_locus_read_supported",
    "phased_or_long_molecule_linkage",
    "validated_isolate_genome",
)


class SchemaViolation(ValueError):
    """An evidence record that would state more than it measured."""


@dataclass(frozen=True, slots=True)
class BiologicalEvidence:
    """One literature claim about a target, bound to what was studied."""

    kind: str
    statement: str
    #: The organism/strain/region actually studied, not the one detected here.
    studied_entity: str = ""
    endpoint: str = ""
    sources: tuple[str, ...] = ()
    #: Limits of the cited experiment, in its own terms.
    boundary: str = ""

    def __post_init__(self) -> None:
        if self.kind not in BIOLOGICAL_EVIDENCE_KINDS:
            raise SchemaViolation(f"unknown biological evidence kind {self.kind!r}")

    def to_json(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "statement": self.statement,
            "studied_entity": self.studied_entity,
            "endpoint": self.endpoint,
            "sources": list(self.sources),
            "boundary": self.boundary,
        }


@dataclass(frozen=True, slots=True)
class Coverage:
    """Sequence support. Every field is nullable and none defaults to zero.

    `independent_fragments` counts paired fragments, not reads and not
    alignments: overlapping mates and repeated alignments of one fragment are
    not independent confirmations (spec §9.1).
    """

    independent_fragments: int | None = None
    unique_fragments: int | None = None
    target_breadth: float | None = None
    breadth_at_required_depth: float | None = None
    median_depth: float | None = None
    callable_bases: int | None = None
    discriminatory_bases: int | None = None
    minor_component_fraction: float | None = None
    limit_of_detection: float | None = None
    limit_scope: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "independent_fragments": self.independent_fragments,
            "unique_fragments": self.unique_fragments,
            "target_breadth": self.target_breadth,
            "breadth_at_required_depth": self.breadth_at_required_depth,
            "median_depth": self.median_depth,
            "callable_bases": self.callable_bases,
            "discriminatory_bases": self.discriminatory_bases,
            "minor_component_fraction": self.minor_component_fraction,
            "limit_of_detection": self.limit_of_detection,
            "limit_scope": self.limit_scope,
        }


@dataclass(frozen=True, slots=True)
class Linkage:
    """Whether the target's sequence can be attributed to a carrier."""

    state: str = "unassigned"
    carrier_call_id: str | None = None
    evidence_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.state not in LINKAGE_STATES:
            raise SchemaViolation(f"unknown linkage state {self.state!r}")

    def to_json(self) -> dict[str, Any]:
        return {
            "state": self.state,
            "carrier_call_id": self.carrier_call_id,
            "evidence_ids": list(self.evidence_ids),
        }


@dataclass(frozen=True, slots=True)
class Validation:
    """What is established about the assay and the claim, separately."""

    analytical_status: str = "not_validated"
    reference_function_status: str = "not_characterized"
    phenotype_measured_in_sample: bool = False
    clinical_predictive_status: str = "not_established"

    def to_json(self) -> dict[str, Any]:
        return {
            "analytical_status": self.analytical_status,
            "reference_function_status": self.reference_function_status,
            "phenotype_measured_in_sample": self.phenotype_measured_in_sample,
            "clinical_predictive_status": self.clinical_predictive_status,
        }


@dataclass(frozen=True, slots=True)
class Provenance:
    """Enough to rebuild the call, or to invalidate it when an input changes."""

    input_sha256: tuple[str, ...] = ()
    database_release: str | None = None
    reference_accession_versions: tuple[str, ...] = ()
    reference_sha256: tuple[str, ...] = ()
    tool_version: str | None = None
    container_digest: str | None = None
    parameters: Mapping[str, Any] = field(default_factory=dict)
    code_commit: str | None = None
    threshold_version: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "input_sha256": list(self.input_sha256),
            "database_release": self.database_release,
            "reference_accession_versions": list(self.reference_accession_versions),
            "reference_sha256": list(self.reference_sha256),
            "tool_version": self.tool_version,
            "container_digest": self.container_digest,
            "parameters": dict(self.parameters),
            "code_commit": self.code_commit,
            "threshold_version": self.threshold_version,
        }


@dataclass(frozen=True, slots=True)
class ResolutionCall:
    """One target's resolution state for one sample.

    This is the single record every consumer reads — profiles, pathogen
    grading, intervention cards, donor comparison, every renderer. They do not
    each re-derive it from a taxonomy table, which is how the same organism
    came to be described differently in different sections.
    """

    sample_id: str
    target_id: str
    identity_kind: str
    assay_status: str = NOT_REQUESTED
    analytical_call: str = UNRESOLVED
    call_id: str | None = None
    reason_codes: tuple[str, ...] = ()
    taxonomic_assertions: tuple[Mapping[str, Any], ...] = ()
    population_components: tuple[Mapping[str, Any], ...] = ()
    reference_candidates: tuple[Mapping[str, Any], ...] = ()
    biological_evidence: tuple[BiologicalEvidence, ...] = ()
    coverage: Coverage = field(default_factory=Coverage)
    linkage: Linkage = field(default_factory=Linkage)
    validation: Validation = field(default_factory=Validation)
    provenance: Provenance = field(default_factory=Provenance)
    #: Plain-language statement of exactly what this call does and does not say.
    plain: str = ""

    def __post_init__(self) -> None:
        if self.assay_status not in ASSAY_STATUSES:
            raise SchemaViolation(f"unknown assay_status {self.assay_status!r}")
        if self.analytical_call not in ANALYTICAL_CALLS:
            raise SchemaViolation(f"unknown analytical_call {self.analytical_call!r}")
        if self.identity_kind not in IDENTITY_KINDS:
            raise SchemaViolation(f"unknown identity_kind {self.identity_kind!r}")

        # The invariant this whole module exists to enforce: an assay that did
        # not complete cannot report what is or is not there.
        if self.assay_status not in CONCLUSIVE_ASSAY and self.analytical_call != UNRESOLVED:
            raise SchemaViolation(
                f"{self.target_id}: assay_status {self.assay_status!r} cannot carry "
                f"analytical_call {self.analytical_call!r} \u2014 only a completed assay "
                "may state a result, and an unrun assay is never a negative screen"
            )
        # A calibrated negative has to name its limit, or it is not calibrated.
        if self.analytical_call == NOT_DETECTED_CALIBRATED and (
            self.coverage.limit_of_detection is None or not self.coverage.limit_scope
        ):
            raise SchemaViolation(
                f"{self.target_id}: a calibrated non-detection requires "
                "limit_of_detection and limit_scope; use "
                f"{NOT_DETECTED_UNVALIDATED!r} when the limit is not established"
            )
        if self.call_id is None and self.assay_status == COMPLETED:
            raise SchemaViolation(
                f"{self.target_id}: a completed assay must carry a call_id so its "
                "result binds to its reads, reference and thresholds"
            )

    # -- derived views ----------------------------------------------------- #

    @property
    def is_negative(self) -> bool:
        return self.analytical_call in NEGATIVE_CALLS

    @property
    def is_positive(self) -> bool:
        return self.analytical_call in POSITIVE_CALLS

    @property
    def terminal(self) -> bool:
        """Has this target reached an accountable end state (§12.2, test 2)?"""
        return self.assay_status != SCHEDULED and self.assay_status != NOT_REQUESTED

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": RESOLUTION_SCHEMA_VERSION,
            "sample_id": self.sample_id,
            "call_id": self.call_id,
            "target_id": self.target_id,
            "identity_kind": self.identity_kind,
            "assay_status": self.assay_status,
            "analytical_call": self.analytical_call,
            "reason_codes": list(self.reason_codes),
            "taxonomic_assertions": [dict(a) for a in self.taxonomic_assertions],
            "population_components": [dict(c) for c in self.population_components],
            "reference_candidates": [dict(r) for r in self.reference_candidates],
            "biological_evidence": [e.to_json() for e in self.biological_evidence],
            "coverage": self.coverage.to_json(),
            "linkage": self.linkage.to_json(),
            "validation": self.validation.to_json(),
            "provenance": self.provenance.to_json(),
            "plain": self.plain,
            # Present, always null, and never computed. Listed explicitly so
            # that a consumer looking for it finds a null rather than
            # inventing one (§9.1, acceptance 70).
            "disease_transmission_probability": None,
        }


def census(calls: Sequence[ResolutionCall]) -> dict[str, Any]:
    """Completeness census: every requested target's terminal state (§12.2).

    A partial run has to say it is partial. This is the object that makes that
    impossible to hide behind a summary.
    """
    by_assay: dict[str, int] = {}
    by_call: dict[str, int] = {}
    for call in calls:
        by_assay[call.assay_status] = by_assay.get(call.assay_status, 0) + 1
        by_call[call.analytical_call] = by_call.get(call.analytical_call, 0) + 1

    outstanding = [c.target_id for c in calls if not c.terminal]
    return {
        "targets_registered": len(calls),
        "targets_terminal": sum(1 for c in calls if c.terminal),
        "targets_outstanding": len(outstanding),
        "outstanding_target_ids": sorted(outstanding),
        "by_assay_status": dict(sorted(by_assay.items())),
        "by_analytical_call": dict(sorted(by_call.items())),
        "complete": not outstanding,
        "what_this_is": (
            "Every target the run registered, with the state it actually reached. A "
            "target that was never run, failed, or lacked a reference appears here "
            "rather than as a negative result, so a partial run cannot present itself "
            "as a complete screen."
        ),
    }
