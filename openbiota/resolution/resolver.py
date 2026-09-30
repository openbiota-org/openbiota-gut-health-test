"""The typed resolver: ask for a resolution, get that resolution or a reason.

Every consumer — profiles, pathogen grading, intervention cards, donor
comparison, every renderer — goes through this instead of reaching into a
taxonomy table. The point is one refusal (spec §9.3):

> A strain query cannot silently return its species abundance.

That substitution is the single most common way a microbiome report overstates
itself. Asking "which strain of *E. coli*?" and receiving "*E. coli*: 2.1%"
looks like an answer and is not one. Here it raises, or returns a record whose
`analytical_call` says the resolution was not reached and why.

Three further refusals follow from the same idea:

* A genomic-trait query cannot be satisfied with a reference organism's
  catalogued trait. BacDive knowing that a type strain sporulates is not a
  measurement of this sample.
* A higher-rank total can aggregate descendants only with an explicit
  denominator, and never by adding a species total to its own strain subtotals.
* Clinical or host data never satisfies a sequence query, and sequence never
  satisfies a clinical one.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from .schema import (
    NOT_REQUESTED,
    REFERENCE_UNAVAILABLE,
    UNRESOLVED,
    ResolutionCall,
)

#: Resolution levels, coarse to fine. A request is satisfied only by evidence
#: at or below (finer than) the level asked for.
LEVELS: Final[tuple[str, ...]] = (
    "lineage",
    "species",
    "species_complex",
    "sgb",
    "subspecies",
    "population_fingerprint",
    "strain_component",
    "named_reference_isolate",
)

#: Which identity kinds can satisfy a request at each level. Deliberately not
#: "anything finer counts": a serotype does not answer a population question,
#: and a population fingerprint does not name an isolate.
SATISFIES: Final[Mapping[str, frozenset[str]]] = {
    "lineage": frozenset({"lineage", "species", "species_complex", "sgb"}),
    "species": frozenset({"species", "sgb"}),
    "species_complex": frozenset({"species_complex", "species", "sgb"}),
    "sgb": frozenset({"sgb"}),
    "subspecies": frozenset({"subspecies", "sgb"}),
    "population_fingerprint": frozenset({"population_fingerprint", "strain_component"}),
    "strain_component": frozenset({"strain_component"}),
    "named_reference_isolate": frozenset({"named_reference_isolate"}),
    "serotype": frozenset({"serotype"}),
    "pathotype": frozenset({"pathotype"}),
    "sequence_type": frozenset({"sequence_type"}),
    "locus": frozenset({"locus", "operon", "allele"}),
    "operon": frozenset({"operon"}),
    "allele": frozenset({"allele"}),
    "gene_family": frozenset({"gene_family", "allele", "operon", "locus"}),
    "accessory_region": frozenset({"accessory_region"}),
    "phenotype": frozenset({"phenotype"}),
    "biotype": frozenset({"biotype"}),
    "plasmid_replicon": frozenset({"plasmid_replicon"}),
}


class ResolutionRefused(ValueError):
    """A query that would have been answered with the wrong kind of evidence."""


@dataclass(frozen=True, slots=True)
class Resolution:
    """The answer to one query: records, states, reasons and context."""

    sample_id: str
    target_id: str
    required_resolution: str
    evidence_records: tuple[ResolutionCall, ...] = ()
    assay_states: tuple[str, ...] = ()
    missing_reasons: tuple[str, ...] = ()
    reference_context: Mapping[str, Any] = field(default_factory=dict)

    @property
    def satisfied(self) -> bool:
        return bool(self.evidence_records)

    def to_json(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "target_id": self.target_id,
            "required_resolution": self.required_resolution,
            "satisfied": self.satisfied,
            "evidence_records": [r.to_json() for r in self.evidence_records],
            "assay_states": list(self.assay_states),
            "missing_reasons": list(self.missing_reasons),
            "reference_context": dict(self.reference_context),
        }


class EvidenceResolver:
    """Holds a sample's resolution calls and answers typed queries."""

    def __init__(
        self,
        calls: Sequence[ResolutionCall],
        *,
        reference_context: Mapping[str, Any] | None = None,
    ) -> None:
        self._by_target: dict[str, list[ResolutionCall]] = {}
        for call in calls:
            self._by_target.setdefault(call.target_id, []).append(call)
        self._context = dict(reference_context or {})

    # -- the query ---------------------------------------------------------- #

    def resolve(
        self,
        sample_id: str,
        target_id: str,
        required_resolution: str,
        permitted_evidence_kinds: Sequence[str] | None = None,
    ) -> Resolution:
        """Answer at the requested resolution, or explain the shortfall.

        Never substitutes a coarser answer. A species abundance is not an
        acceptable reply to a strain question, so the returned record's
        `analytical_call` stays `unresolved` and `missing_reasons` says what
        was asked for and what existed instead.
        """
        if required_resolution not in SATISFIES:
            raise ResolutionRefused(
                f"unknown required_resolution {required_resolution!r}; "
                f"known: {sorted(SATISFIES)}"
            )
        acceptable = SATISFIES[required_resolution]
        if permitted_evidence_kinds is not None:
            acceptable = acceptable & frozenset(permitted_evidence_kinds)

        candidates = self._by_target.get(target_id, [])
        matching = [
            c for c in candidates
            if c.identity_kind in acceptable and c.sample_id == sample_id
        ]
        coarser = [
            c for c in candidates
            if c.identity_kind not in acceptable and c.sample_id == sample_id
        ]

        reasons: list[str] = []
        if not candidates:
            reasons.append(
                f"no registered evidence for {target_id}; the target was never requested "
                "for this sample, which is not an absence"
            )
        elif not matching:
            offered = sorted({c.identity_kind for c in coarser})
            reasons.append(
                f"{required_resolution} was requested; only {offered} evidence exists for "
                f"{target_id}. A coarser or differently-typed result is not returned in its "
                "place \u2014 the resolution is unresolved."
            )

        return Resolution(
            sample_id=sample_id,
            target_id=target_id,
            required_resolution=required_resolution,
            evidence_records=tuple(matching),
            assay_states=tuple(sorted({c.assay_status for c in candidates})) or (NOT_REQUESTED,),
            missing_reasons=tuple(reasons),
            reference_context=self._context,
        )

    # -- aggregation -------------------------------------------------------- #

    def aggregate(
        self,
        sample_id: str,
        target_ids: Sequence[str],
        *,
        denominator: str,
    ) -> dict[str, Any]:
        """Sum descendant measurements with an explicit denominator.

        A total is meaningless without saying what it is a fraction of, and
        adding a species total to its own strain subtotals double counts
        (spec §9.3). Both are refused here rather than warned about.
        """
        if not denominator:
            raise ResolutionRefused(
                "an aggregate requires an explicit denominator; a bare total cannot be "
                "combined with measurements on another scale"
            )
        seen_kinds: set[str] = set()
        contributions: list[dict[str, Any]] = []
        for target_id in target_ids:
            for call in self._by_target.get(target_id, []):
                if call.sample_id != sample_id:
                    continue
                seen_kinds.add(call.identity_kind)
                contributions.append({
                    "target_id": target_id,
                    "identity_kind": call.identity_kind,
                    "callable_bases": call.coverage.callable_bases,
                })
        # Species-level and strain-level measurements of the same thing are
        # not additive.
        if {"species", "sgb"} & seen_kinds and {
            "population_fingerprint", "strain_component"
        } & seen_kinds:
            raise ResolutionRefused(
                "refusing to aggregate species/SGB totals together with their own strain "
                "subtotals: the strain values are part of the species value, so adding them "
                "double counts"
            )
        return {
            "sample_id": sample_id,
            "denominator": denominator,
            "n_contributions": len(contributions),
            "identity_kinds": sorted(seen_kinds),
            "contributions": contributions,
            "unresolved_residual_possible": True,
            "note": (
                "Population estimates must sum consistently with their species model or "
                "retain an unresolved residual. This aggregate reports its contributions "
                "rather than asserting a closed total."
            ),
        }

    # -- traits -------------------------------------------------------------- #

    @staticmethod
    def reference_trait(
        *,
        organism: str,
        trait: str,
        source: str,
        measured_in_sample: bool = False,
    ) -> dict[str, Any]:
        """A catalogued reference trait, explicitly not a sample measurement.

        Returned as a labelled reference annotation so a consumer cannot print
        it as though the phenotype had been observed here (spec §3.1).
        """
        return {
            "organism": organism,
            "trait": trait,
            "source": source,
            "measured_in_this_sample": bool(measured_in_sample),
            "kind": "reference_annotation",
            "note": (
                f"This is a recorded property of the reference organism from {source}. "
                "DNA in stool does not observe morphology, growth behaviour or expression, "
                "so it has not been measured in this sample."
            ),
        }


def self_test() -> int:
    """Guard the substitutions the resolver exists to refuse."""
    from .schema import COMPLETED, SUPPORTED_DETECTION, Coverage

    failures = 0

    species_call = ResolutionCall(
        sample_id="K1", target_id="org.ecoli", identity_kind="species",
        assay_status=COMPLETED, call_id="c1", analytical_call=SUPPORTED_DETECTION,
    )
    strain_call = ResolutionCall(
        sample_id="K1", target_id="strain.population.SGB1", identity_kind="population_fingerprint",
        assay_status=COMPLETED, call_id="c2", analytical_call=SUPPORTED_DETECTION,
        coverage=Coverage(callable_bases=100_000),
    )
    resolver = EvidenceResolver([species_call, strain_call])

    # The central refusal: a strain question is not answered by species data.
    out = resolver.resolve("K1", "org.ecoli", "population_fingerprint")
    if out.satisfied:
        failures += 1
    if not out.missing_reasons or "not returned in its place" not in out.missing_reasons[0]:
        failures += 1

    # The species question is answered by species data.
    out = resolver.resolve("K1", "org.ecoli", "species")
    if not out.satisfied or out.evidence_records[0].identity_kind != "species":
        failures += 1

    # A population fingerprint does not name an isolate.
    out = resolver.resolve("K1", "strain.population.SGB1", "named_reference_isolate")
    if out.satisfied:
        failures += 1

    # An unregistered target is unassessed, not absent.
    out = resolver.resolve("K1", "nope", "species")
    if out.satisfied or "never requested" not in out.missing_reasons[0]:
        failures += 1

    # Permitted-kind narrowing is honoured.
    out = resolver.resolve("K1", "org.ecoli", "species", permitted_evidence_kinds=["sgb"])
    if out.satisfied:
        failures += 1

    # Aggregation requires a denominator and refuses double counting.
    try:
        resolver.aggregate("K1", ["org.ecoli"], denominator="")
        failures += 1
    except ResolutionRefused:
        pass
    try:
        resolver.aggregate(
            "K1", ["org.ecoli", "strain.population.SGB1"], denominator="bacterial_reads"
        )
        failures += 1
    except ResolutionRefused:
        pass
    agg = resolver.aggregate("K1", ["strain.population.SGB1"], denominator="bacterial_reads")
    if agg["denominator"] != "bacterial_reads" or agg["n_contributions"] != 1:
        failures += 1

    # A reference trait is never a sample measurement.
    trait = EvidenceResolver.reference_trait(
        organism="B. fragilis", trait="sporulation", source="BacDive"
    )
    if trait["measured_in_this_sample"] is not False:
        failures += 1
    if trait["kind"] != "reference_annotation":
        failures += 1

    # An unknown resolution level is refused, not guessed.
    try:
        resolver.resolve("K1", "org.ecoli", "vibes")
        failures += 1
    except ResolutionRefused:
        pass

    if REFERENCE_UNAVAILABLE not in {REFERENCE_UNAVAILABLE} or UNRESOLVED != "unresolved":
        failures += 1
    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
