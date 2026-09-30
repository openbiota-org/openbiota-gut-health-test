"""Species interpretation registry (spec v04 §6.5).

Every species the profiler can name gets an explicit, auditable disposition:
what the literature says a higher or lower level has been *associated* with,
in which direction, with what caveats, and whether any curated human
intervention evidence exists for it. Level and interpretation are separate
axes — a high level may be favourable, concerning, neutral or unknown — and
the registry never carries treatment records; those live in the intervention
registry and are bound through typed triggers.

Entries live in ``taxa/interpretations/*.yaml``; a species with no entry is
reported with ``evidence_state: no_curated_interpretation`` rather than
silently given the nearest neighbour's meaning.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.errors import PanelError

POLARITIES: Final = (
    "health_associated",
    "disease_associated",
    "potential_pathobiont",
    "validated_pathogen_signal",
    "context_dependent",
    "conflicting_evidence",
    "no_curated_interpretation",
)

EVIDENCE_STATES: Final = (
    "curated_human_action_evidence",
    "human_association_only",
    "mechanistic_only",
    "no_curated_interpretation",
)

#: Interpretation summary → the reading-guide label used in the report.
POLARITY_LABEL: Final = {
    "health_associated": "health-associated in the literature",
    "disease_associated": "disease-associated in the literature",
    "potential_pathobiont": "potential pathobiont (context-dependent)",
    "validated_pathogen_signal": "pathogen signal requires clinical confirmation",
    "context_dependent": "context-dependent",
    "conflicting_evidence": "conflicting evidence",
    "no_curated_interpretation": "no curated interpretation",
}


@dataclass(frozen=True, slots=True)
class Assertion:
    """One literature-backed statement about a species."""

    assertion_id: str
    statement: str
    direction: str  # "higher" | "lower" | "presence" | "any"
    citation: str
    population: str = "adults, stool"
    claim: str = "association"  # association | mechanism | intervention_context

    def to_json(self) -> dict[str, Any]:
        return {
            "assertion_id": self.assertion_id, "statement": self.statement,
            "direction": self.direction, "citation": self.citation,
            "population": self.population, "claim": self.claim,
        }


@dataclass(frozen=True, slots=True)
class SpeciesInterpretation:
    species: str
    display_name: str
    polarity: str
    summary: str
    evidence_state: str
    assertions: tuple[Assertion, ...] = ()
    tags: tuple[str, ...] = ()
    #: Names of curated taxon groups (guilds) it belongs to, for the card.
    guilds: tuple[str, ...] = ()
    #: What a concern label would need before it can be more than a flag.
    qualification_needed: str | None = None
    #: Phrases the renderer must never use for this organism.
    prohibited_language: tuple[str, ...] = ()
    last_reviewed: str = ""
    source_path: Path | None = field(default=None, compare=False)

    @property
    def concern(self) -> bool:
        return self.polarity in ("potential_pathobiont", "validated_pathogen_signal")

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "display_name": self.display_name,
            "polarity": self.polarity,
            "polarity_label": POLARITY_LABEL[self.polarity],
            "summary": self.summary,
            "evidence_state": self.evidence_state,
            "assertions": [a.to_json() for a in self.assertions],
            "tags": list(self.tags),
            "guilds": list(self.guilds),
            "qualification_needed": self.qualification_needed,
            "prohibited_language": list(self.prohibited_language),
            "last_reviewed": self.last_reviewed,
        }


def _req(data: dict[str, Any], key: str, where: str) -> Any:
    if key not in data or data[key] in (None, ""):
        raise PanelError(f"{where}: '{key}' is required")
    return data[key]


def parse_species_interpretation(data: dict[str, Any], *, source_path: Path | None = None) -> SpeciesInterpretation:
    where = f"interpretation {data.get('species', '?')}"
    species = str(_req(data, "species", where))
    polarity = str(_req(data, "polarity", where))
    if polarity not in POLARITIES:
        raise PanelError(f"{where}: polarity must be one of {list(POLARITIES)}, got {polarity!r}")
    state = str(_req(data, "evidence_state", where))
    if state not in EVIDENCE_STATES:
        raise PanelError(f"{where}: evidence_state must be one of {list(EVIDENCE_STATES)}, got {state!r}")
    assertions = []
    for i, raw in enumerate(data.get("assertions") or []):
        aw = f"{where} assertion[{i}]"
        direction = str(raw.get("direction", "any"))
        if direction not in ("higher", "lower", "presence", "any"):
            raise PanelError(f"{aw}: direction must be higher/lower/presence/any")
        assertions.append(Assertion(
            assertion_id=str(raw.get("assertion_id") or f"INTERPRET_{species.upper()}_{i + 1}"),
            statement=str(_req(raw, "statement", aw)),
            direction=direction,
            citation=str(_req(raw, "citation", aw)),
            population=str(raw.get("population", "adults, stool")),
            claim=str(raw.get("claim", "association")),
        ))
    if polarity != "no_curated_interpretation" and not assertions:
        raise PanelError(f"{where}: a curated polarity needs at least one cited assertion")
    if polarity in ("potential_pathobiont", "validated_pathogen_signal") and not data.get("qualification_needed"):
        raise PanelError(f"{where}: a concern polarity must say what qualification it needs")
    return SpeciesInterpretation(
        species=species,
        display_name=str(data.get("display_name") or species.replace("_", " ")),
        polarity=polarity,
        summary=str(_req(data, "summary", where)),
        evidence_state=state,
        assertions=tuple(assertions),
        tags=tuple(str(t) for t in data.get("tags") or ()),
        guilds=tuple(str(g) for g in data.get("guilds") or ()),
        qualification_needed=data.get("qualification_needed"),
        prohibited_language=tuple(str(p) for p in data.get("prohibited_language") or ()),
        last_reviewed=str(data.get("last_reviewed", "")),
        source_path=source_path,
    )


@dataclass(frozen=True, slots=True)
class InterpretationRegistry:
    entries: dict[str, SpeciesInterpretation]
    namespace: str = "metaphlan_species"
    namespace_version: str = "mpa_v31_CHOCOPhlAn_201901"

    def get(self, species: str) -> SpeciesInterpretation | None:
        return self.entries.get(species)

    def disposition(self, species: str) -> SpeciesInterpretation:
        """Always an answer: the curated entry or an explicit 'none'."""
        found = self.entries.get(species)
        if found is not None:
            return found
        return SpeciesInterpretation(
            species=species, display_name=species.replace("_", " "),
            polarity="no_curated_interpretation",
            summary="No curated interpretation in this registry version.",
            evidence_state="no_curated_interpretation",
        )

    def __len__(self) -> int:
        return len(self.entries)

    def concern_species(self) -> Iterable[str]:
        return (s for s, e in self.entries.items() if e.concern)


def load_interpretation_registry(directory: Path) -> InterpretationRegistry:
    """Load every ``*.yaml`` in ``directory``; each file may hold one entry or a list."""
    entries: dict[str, SpeciesInterpretation] = {}
    if not directory.is_dir():
        return InterpretationRegistry(entries={})
    for path in sorted(directory.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        items = data if isinstance(data, list) else (data.get("species_interpretations") if isinstance(data, dict) and "species_interpretations" in data else [data])
        for item in items or []:
            entry = parse_species_interpretation(item, source_path=path)
            if entry.species in entries:
                raise PanelError(f"duplicate species interpretation for {entry.species} in {path}")
            entries[entry.species] = entry
    return InterpretationRegistry(entries=entries)


__all__ = [
    "EVIDENCE_STATES",
    "POLARITIES",
    "POLARITY_LABEL",
    "Assertion",
    "InterpretationRegistry",
    "SpeciesInterpretation",
    "load_interpretation_registry",
    "parse_species_interpretation",
]
