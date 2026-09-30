"""The gene and reaction registry from BUILD_SPEC_v0.8.3 section 3.4.

The rule this module exists to enforce: **a gene symbol is not a unique
biochemical identity.** `urdA` is not a urolithin marker, its substrate is
urocanate; `acs` supports acetate assimilation rather than formation; GH13
supports amylolysis in general and no resistant starch in particular. Three
separate report claims have gone wrong that way, so identity here is an
accession plus a sequence hash plus a reaction, never a display name.

Two structures do the work.

`Marker` is one reference sequence with the provenance that makes it usable:
accession *and version*, the hash of the sequence actually installed, the
publication that characterised it, the organism it came from, its substrate
and product, and the homologs it must be told apart from. A marker admitted
because its annotation contained a promising word is an unvalidated candidate
until a specificity fixture has been run against it, and it says so.

`Reaction` is a biochemical step, joined to markers by accession. Metrics bind
to reaction IDs. Nothing in the report may bind by matching a display-name
substring, so `bind_by_name` does not exist.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

#: How a marker's specificity currently stands. A marker may be *installed*
#: and *unvalidated* at the same time; the report has to be able to say so.
SPECIFICITY_STATES: Final[frozenset[str]] = frozenset({
    "validated_against_negatives",   # positives and hard negatives both run
    "positives_only",                # found its targets, negatives not yet run
    "candidate_unvalidated",         # admitted from annotation; nothing checked
    "rejected_nonspecific",          # failed: hits its own negative set
})

#: Evidence that this sequence performs this reaction.
FUNCTION_EVIDENCE: Final[frozenset[str]] = frozenset({
    "purified_enzyme_assay",
    "heterologous_expression",
    "gene_knockout",
    "genetic_transfer",
    "structural_with_activity",
    "operon_transcription",
    "orthology_inference",
})

#: Evidence levels that support naming a specific reaction. `orthology_inference`
#: alone supports a family, which is a different and weaker claim.
REACTION_SPECIFIC_EVIDENCE: Final[frozenset[str]] = frozenset(
    FUNCTION_EVIDENCE - {"orthology_inference"}
)

#: Direction of the catalysed step, because a reversible enzyme running the
#: other way is a different measurement: `acs` assimilating acetate is not
#: `pta`/`ackA` forming it.
REACTION_DIRECTIONS: Final[frozenset[str]] = frozenset({
    "forward", "reverse", "reversible", "direction_unresolved",
})

#: Logic joining a reaction's markers. `AND` needs every subunit; `OR` accepts
#: any one alternative. Applied in `openbiota.extension.schema`.
STEP_LOGIC: Final[frozenset[str]] = frozenset({"AND", "OR", "SINGLE"})

_ACCESSION_RE: Final = re.compile(
    r"^(?:[A-Z]{2}_?\d{6,}(?:\.\d+)?"          # GenBank/RefSeq nucleotide
    r"|[A-Z]{3}\d{5,}(?:\.\d+)?"               # protein
    r"|[A-Z]P_\d{6,}(?:\.\d+)?"                # WP_/NP_/XP_
    r"|GC[AF]_\d{9}\.\d+"                      # assembly
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"             # UniProt
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r"|K\d{5}|COG\d{4}|TIGR\d{5}|GH\d+(?:_\d+)?|PL\d+|CE\d+|CBM\d+|MER\d{7}"
    r"|BT\d{4}|Amuc_\d{4})$"
)


class MarkerRegistryError(ValueError):
    """A marker or reaction record that would misstate its own identity."""


def sequence_hash(sequence: str) -> str:
    """SHA-256 of an uppercase, unwrapped sequence.

    One normalisation, fixed here, so a hash computed at ingestion and a hash
    recomputed at verification can actually be compared. FASTA line wrapping
    and lower-case soft-masking are removed; nothing else is touched.
    """
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise MarkerRegistryError("empty sequence")
    return hashlib.sha256(cleaned.encode("ascii", "strict")).hexdigest()


@dataclass(slots=True, frozen=True)
class Marker:
    """One reference sequence, with the provenance that makes it usable."""

    marker_id: str
    symbol: str
    accession: str
    sequence_type: str                      # protein | nucleotide | hmm | assembly
    reaction_ids: tuple[str, ...]
    source_ids: tuple[str, ...]
    organism: str
    function_evidence: str
    specificity_state: str = "candidate_unvalidated"
    sequence_sha256: str | None = None
    amino_acids: int | None = None
    nucleotides: int | None = None
    orthology_families: tuple[str, ...] = ()
    diagnostic_residues: tuple[str, ...] = ()
    residue_numbering_reference: str | None = None
    negative_homologs: tuple[str, ...] = ()
    qualifying_context: str | None = None
    strain: str | None = None
    note: str | None = None

    def __post_init__(self) -> None:
        if self.sequence_type not in {"protein", "nucleotide", "hmm", "assembly"}:
            raise MarkerRegistryError(f"{self.marker_id}: unknown sequence type")
        if self.function_evidence not in FUNCTION_EVIDENCE:
            raise MarkerRegistryError(
                f"{self.marker_id}: unknown function evidence {self.function_evidence!r}"
            )
        if self.specificity_state not in SPECIFICITY_STATES:
            raise MarkerRegistryError(
                f"{self.marker_id}: unknown specificity state {self.specificity_state!r}"
            )
        if not self.reaction_ids:
            raise MarkerRegistryError(f"{self.marker_id}: a marker must join at least one reaction")
        if not self.source_ids:
            raise MarkerRegistryError(f"{self.marker_id}: a marker must cite its source")
        if self.sequence_type in {"protein", "nucleotide"} and not _ACCESSION_RE.match(self.accession):
            raise MarkerRegistryError(
                f"{self.marker_id}: {self.accession!r} is not a recognisable versioned accession"
            )
        # Diagnostic residues are only meaningful against a stated numbering.
        if self.diagnostic_residues and not self.residue_numbering_reference:
            raise MarkerRegistryError(
                f"{self.marker_id}: diagnostic residues require their numbering reference; "
                "a raw read offset is not a reference position"
            )
        # Claiming validation requires something to have been validated against.
        if self.specificity_state == "validated_against_negatives" and not self.negative_homologs:
            raise MarkerRegistryError(
                f"{self.marker_id}: cannot claim negative validation with no negative set"
            )

    @property
    def is_reaction_specific(self) -> bool:
        """Whether this marker can support a named reaction on its own."""
        return (
            self.function_evidence in REACTION_SPECIFIC_EVIDENCE
            and self.specificity_state == "validated_against_negatives"
        )

    @property
    def effective_length_nt(self) -> int | None:
        """Nucleotide length for RPKM, converting a protein reference once."""
        if self.nucleotides:
            return self.nucleotides
        if self.amino_acids:
            from openbiota.extension.schema import cds_equivalent_length

            return cds_equivalent_length(self.amino_acids)
        return None

    def to_json(self) -> dict[str, Any]:
        return {
            "marker_id": self.marker_id,
            "symbol": self.symbol,
            "accession": self.accession,
            "sequence_type": self.sequence_type,
            "sequence_sha256": self.sequence_sha256,
            "amino_acids": self.amino_acids,
            "nucleotides": self.nucleotides,
            "effective_length_nt": self.effective_length_nt,
            "reaction_ids": list(self.reaction_ids),
            "organism": self.organism,
            "strain": self.strain,
            "function_evidence": self.function_evidence,
            "specificity_state": self.specificity_state,
            "is_reaction_specific": self.is_reaction_specific,
            "orthology_families": list(self.orthology_families),
            "diagnostic_residues": list(self.diagnostic_residues),
            "residue_numbering_reference": self.residue_numbering_reference,
            "negative_homologs": list(self.negative_homologs),
            "qualifying_context": self.qualifying_context,
            "source_ids": list(self.source_ids),
            "note": self.note,
        }


@dataclass(slots=True, frozen=True)
class Reaction:
    """One biochemical step. Metrics bind here, by ID, never by name."""

    reaction_id: str
    label: str
    substrate: str
    product: str
    direction: str
    source_ids: tuple[str, ...]
    step_logic: str = "SINGLE"
    ec_numbers: tuple[str, ...] = ()
    must_not_be_confused_with: tuple[str, ...] = ()
    note: str | None = None

    def __post_init__(self) -> None:
        if self.direction not in REACTION_DIRECTIONS:
            raise MarkerRegistryError(f"{self.reaction_id}: unknown direction {self.direction!r}")
        if self.step_logic not in STEP_LOGIC:
            raise MarkerRegistryError(f"{self.reaction_id}: unknown step logic {self.step_logic!r}")
        if not self.source_ids:
            raise MarkerRegistryError(f"{self.reaction_id}: a reaction must cite its source")
        if not self.substrate or not self.product:
            raise MarkerRegistryError(
                f"{self.reaction_id}: a reaction needs both substrate and product; "
                "a step with an unnamed substrate is a gene family, not a reaction"
            )

    def to_json(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id,
            "label": self.label,
            "substrate": self.substrate,
            "product": self.product,
            "direction": self.direction,
            "step_logic": self.step_logic,
            "ec_numbers": list(self.ec_numbers),
            "must_not_be_confused_with": list(self.must_not_be_confused_with),
            "source_ids": list(self.source_ids),
            "note": self.note,
        }


@dataclass(slots=True)
class MarkerRegistry:
    """Reactions and their markers, joined by ID."""

    reactions: dict[str, Reaction] = field(default_factory=dict)
    markers: dict[str, Marker] = field(default_factory=dict)

    def add_reaction(self, reaction: Reaction) -> Reaction:
        if reaction.reaction_id in self.reactions:
            raise MarkerRegistryError(f"duplicate reaction {reaction.reaction_id}")
        self.reactions[reaction.reaction_id] = reaction
        return reaction

    def add_marker(self, marker: Marker) -> Marker:
        if marker.marker_id in self.markers:
            raise MarkerRegistryError(f"duplicate marker {marker.marker_id}")
        unknown = [r for r in marker.reaction_ids if r not in self.reactions]
        if unknown:
            raise MarkerRegistryError(
                f"{marker.marker_id}: joins unknown reaction(s) {unknown}; "
                "markers bind to reaction IDs, and the reaction must exist first"
            )
        self.markers[marker.marker_id] = marker
        return marker

    def markers_for(self, reaction_id: str) -> tuple[Marker, ...]:
        """Every marker joined to one reaction, by ID."""
        if reaction_id not in self.reactions:
            raise MarkerRegistryError(f"unknown reaction {reaction_id}")
        return tuple(
            m for m in self.markers.values() if reaction_id in m.reaction_ids
        )

    def negative_set_for(self, reaction_id: str) -> tuple[str, ...]:
        """The homologs this reaction's markers must be told apart from."""
        out: set[str] = set()
        for marker in self.markers_for(reaction_id):
            out.update(marker.negative_homologs)
        return tuple(sorted(out))

    def unvalidated(self) -> tuple[Marker, ...]:
        """Markers still awaiting a specificity fixture, named rather than hidden."""
        return tuple(
            m for m in self.markers.values()
            if m.specificity_state in {"candidate_unvalidated", "positives_only"}
        )

    def audit(self) -> dict[str, Any]:
        """A registry-level honesty report for the capability coverage output."""
        by_state: dict[str, int] = {}
        for m in self.markers.values():
            by_state[m.specificity_state] = by_state.get(m.specificity_state, 0) + 1
        orphans = sorted(r for r in self.reactions if not self.markers_for(r))
        return {
            "n_reactions": len(self.reactions),
            "n_markers": len(self.markers),
            "markers_by_specificity_state": dict(sorted(by_state.items())),
            "reactions_without_markers": orphans,
            "n_reaction_specific_markers": sum(
                1 for m in self.markers.values() if m.is_reaction_specific
            ),
        }

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": "openbiota.marker-registry/1.0",
            "reactions": {k: v.to_json() for k, v in sorted(self.reactions.items())},
            "markers": {k: v.to_json() for k, v in sorted(self.markers.items())},
            "audit": self.audit(),
        }


def verify_sequence(marker: Marker, sequence: str) -> None:
    """Check a downloaded sequence against the marker's recorded hash.

    Used by the asset verifier. A retrieval that returns the wrong record - a
    different isoform, a truncated entry, an HTML error page saved as FASTA -
    fails here rather than silently becoming the reference for a reaction.
    """
    if not marker.sequence_sha256:
        raise MarkerRegistryError(f"{marker.marker_id}: no recorded hash to verify against")
    actual = sequence_hash(sequence)
    if actual != marker.sequence_sha256:
        raise MarkerRegistryError(
            f"{marker.marker_id}: sequence hash mismatch for {marker.accession}\n"
            f"  expected {marker.sequence_sha256}\n  actual   {actual}"
        )
    cleaned = "".join(sequence.split())
    if marker.amino_acids and len(cleaned) != marker.amino_acids:
        raise MarkerRegistryError(
            f"{marker.marker_id}: expected {marker.amino_acids} aa, got {len(cleaned)}"
        )
    if marker.nucleotides and len(cleaned) != marker.nucleotides:
        raise MarkerRegistryError(
            f"{marker.marker_id}: expected {marker.nucleotides} nt, got {len(cleaned)}"
        )


def resolve_step(
    logic: str, supported: Mapping[str, bool], members: Sequence[str]
) -> bool:
    """Whether a route step is satisfied, honouring AND/OR.

    Section 4.1: an AND step needs every subunit; an OR step is satisfied by
    any one non-overlapping alternative. Written once so no panel invents its
    own grammar.
    """
    if logic not in STEP_LOGIC:
        raise MarkerRegistryError(f"unknown step logic {logic!r}")
    if not members:
        return False
    if logic == "OR":
        return any(supported.get(m, False) for m in members)
    return all(supported.get(m, False) for m in members)


def diagnostic_residues_covered(
    marker: Marker, covered_positions: Iterable[int]
) -> tuple[bool, tuple[str, ...]]:
    """Whether the positions that distinguish this marker were actually read.

    Section 5.6 on UrdA: discrimination must use aligned reference numbering
    and *covered* diagnostic positions. Incomplete coverage of the
    distinguishing residues means unresolved specificity, not a negative.
    """
    covered = {int(p) for p in covered_positions}
    missing: list[str] = []
    for spec in marker.diagnostic_residues:
        m = re.search(r"(\d+)", spec)
        if not m:
            raise MarkerRegistryError(
                f"{marker.marker_id}: diagnostic residue {spec!r} has no position"
            )
        if int(m.group(1)) not in covered:
            missing.append(spec)
    return (not missing), tuple(missing)


__all__ = [
    "FUNCTION_EVIDENCE",
    "REACTION_DIRECTIONS",
    "REACTION_SPECIFIC_EVIDENCE",
    "SPECIFICITY_STATES",
    "STEP_LOGIC",
    "Marker",
    "MarkerRegistry",
    "MarkerRegistryError",
    "Reaction",
    "diagnostic_residues_covered",
    "resolve_step",
    "sequence_hash",
    "verify_sequence",
]
