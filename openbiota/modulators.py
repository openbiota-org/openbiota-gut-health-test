"""What each organism does, and what is known to move it.

The registry in ``taxa/classification/modulators.yaml`` is keyed by organism
- a species where the evidence is species-specific, a genus where it is not -
and each entry carries the organism's role in the gut, why its level matters,
and every agent the literature has shown to raise or lower it, with the
evidence level of each stated.

Resolution is species, then alias, then genus. A genus-level entry is marked
as such when it is used, so the card can say "known about the genus" rather
than imply a species-specific study exists.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

import yaml

REGISTRY: Final = Path("taxa/classification/modulators.yaml")

#: Reader-facing words for each evidence level, and its rank for sorting.
EVIDENCE: Final[dict[str, tuple[str, int]]] = {
    "human_trial": ("human trial", 4),
    "human_observational": ("human cohorts", 3),
    "animal": ("animal model", 2),
    "in_vitro": ("laboratory culture", 1),
}

#: Reader-facing words for agent kinds.
KIND: Final[dict[str, str]] = {
    "probiotic": "probiotic",
    "prebiotic": "prebiotic fibre",
    "food": "food",
    "compound": "compound",
    "herb": "herb",
    "drug": "medication",
    "diet": "diet pattern",
    "lifestyle": "lifestyle",
}


@dataclass(frozen=True, slots=True)
class Agent:
    """One thing shown to move an organism."""

    agent: str
    kind: str
    evidence: str
    effect: str
    survives_transit: bool
    note: str = ""
    source: str = ""

    @property
    def evidence_words(self) -> str:
        return EVIDENCE.get(self.evidence, (self.evidence.replace("_", " "), 0))[0]

    @property
    def evidence_rank(self) -> int:
        return EVIDENCE.get(self.evidence, ("", 0))[1]

    @property
    def kind_words(self) -> str:
        return KIND.get(self.kind, self.kind)


@dataclass(frozen=True, slots=True)
class Modulators:
    """Everything the registry knows about one organism."""

    organism: str
    role: str
    why_level_matters: str
    decrease: tuple[Agent, ...]
    increase: tuple[Agent, ...]
    competitors: tuple[str, ...]
    #: "species" when the entry names this organism, "genus" when inherited.
    basis: str = "species"

    def for_direction(self, flag: str) -> tuple[Agent, ...]:
        """The agents that push the organism the way the reader would want.

        A high reading wants what lowers it; a low reading wants what raises it.
        Ordered by evidence strength.
        """
        agents = self.decrease if flag == "high" else self.increase if flag == "low" else ()
        return tuple(sorted(agents, key=lambda a: -a.evidence_rank))

    @property
    def has_any(self) -> bool:
        return bool(self.decrease or self.increase or self.competitors)


@lru_cache(maxsize=1)
def _table() -> dict[str, Modulators]:
    if not REGISTRY.is_file():
        return {}
    data = yaml.safe_load(REGISTRY.read_text()) or {}
    out: dict[str, Modulators] = {}
    for rec in data.get("modulators") or []:
        def agents(items: Any) -> tuple[Agent, ...]:
            return tuple(
                Agent(
                    agent=str(i.get("agent") or ""),
                    kind=str(i.get("kind") or "compound"),
                    evidence=str(i.get("evidence") or "in_vitro"),
                    effect=str(i.get("effect") or ""),
                    survives_transit=bool(i.get("survives_transit", True)),
                    note=" ".join(str(i.get("note") or "").split()),
                    source=" ".join(str(i.get("source") or "").split()),
                )
                for i in (items or []) if isinstance(i, dict)
            )
        m = Modulators(
            organism=str(rec["organism"]),
            role=" ".join(str(rec.get("role") or "").split()),
            why_level_matters=" ".join(str(rec.get("why_level_matters") or "").split()),
            decrease=agents(rec.get("decrease")),
            increase=agents(rec.get("increase")),
            competitors=tuple(str(c) for c in (rec.get("competitors") or [])),
        )
        out[m.organism] = m
        for alias in rec.get("aliases") or []:
            out.setdefault(str(alias), m)
    return out


def lookup(species: str, *, gtdb: str | None = None, formerly: str | None = None) -> Modulators | None:
    """The registry entry for an organism, or its genus, or nothing."""
    table = _table()
    for key in (species, (gtdb or "").replace(" ", "_"), formerly or ""):
        if key and key in table:
            m = table[key]
            # A hit through a genus-keyed entry is genus-level even when
            # the alias list named this species.
            return m if m.organism == key or "_" in m.organism else _as_genus(m)
    genus = species.split("_")[0]
    if genus in table:
        return _as_genus(table[genus])
    return None


def _as_genus(m: Modulators) -> Modulators:
    return Modulators(
        organism=m.organism, role=m.role, why_level_matters=m.why_level_matters,
        decrease=m.decrease, increase=m.increase, competitors=m.competitors, basis="genus",
    )


def available() -> bool:
    return bool(_table())


# --------------------------------------------------------------------------- #
# lever families: the same thing under different names is one lever
# --------------------------------------------------------------------------- #

#: Ordered keyword rules. The first rule whose keyword appears in the agent
#: name (case-insensitive) names the family. Agents matching no rule are
#: their own family. Order matters where keywords overlap: "inulin" must
#: beat "dietary fibre", "pasteurised Akkermansia" must beat "probiotic".
_FAMILIES: Final[tuple[tuple[str, tuple[str, ...]], ...]] = (
    # Exposures first, so "Low-fibre Western diet" is never filed under fibre
    # and "Antibiotics; low-fibre diet" is never filed under anything else.
    ("Antibiotics", ("antibiotic", "vancomycin", "cephalosporin")),
    ("Low-fibre / refined-carbohydrate eating", ("low-fibre", "low fibre", "refined carbohydrate", "western diet", "low-carbohydrate", "low-diversity diet")),
    ("Smoking and gum disease", ("smoking", "periodontal disease", "gum disease")),
    ("Inulin / fructo-oligosaccharides", ("inulin", "fructo-oligo", "fructan", "chicory")),
    ("Resistant starch", ("resistant starch",)),
    ("Pectin (apples, citrus, carrots)", ("pectin",)),
    ("Oat / barley beta-glucan", ("beta-glucan", "oat ", "barley")),
    ("Whole grains and legumes", ("whole grain", "whole-grain", "legume")),
    ("Fruit and vegetables", ("fruit and vegetable",)),
    ("Mediterranean diet", ("mediterranean",)),
    ("Polyphenols (berries, pomegranate, cocoa, tea)", ("polyphenol", "pomegranate", "cranberry", "cocoa", "grape", "green tea")),
    ("Pasteurised Akkermansia muciniphila", ("akkermansia",)),
    ("Saccharomyces boulardii", ("saccharomyces",)),
    ("Bifidobacterium probiotics", ("bifidobacterium", "b. infantis", "b. breve")),
    ("Lactobacillus probiotics", ("lactobacillus", "lgg", "rhamnosus", "plantarum", "reuteri", "salivarius")),
    ("Berberine", ("berberine",)),
    ("Metformin (prescription)", ("metformin",)),
    ("DMB (in some olive and grape-seed oils)", ("dimethyl", "dmb")),
    ("Oral hygiene and dental care", ("oral hygiene", "periodontal treatment")),
    ("Review acid-suppressing medication with a prescriber", ("acid-suppress", "proton-pump", "ppi")),
    ("Low-FODMAP diet (dietitian-supervised)", ("fodmap",)),
    ("Exclusive enteral nutrition (supervised)", ("enteral nutrition",)),
    ("Fasting or caloric restriction", ("fasting", "caloric restriction", "weight loss")),
    ("Less red meat, choline and carnitine", ("choline", "carnitine", "red meat", "red-meat")),
    ("Less arginine-rich protein", ("arginine", "high-protein")),
    ("Less saturated fat and animal-based eating", ("saturated fat", "saturated-fat", "lower dietary fat", "lower fat", "plant-based", "animal-protein", "animal-based", "high-fat", "high fat")),
    ("Fermented dairy", ("yoghurt", "kefir", "fermented dairy", "cheese", "buttermilk")),
    ("Dietary fibre (general)", ("dietary fibre", "fibre diversity", "fermentable fibre", "fibre-rich", "fibre,", "fibre ")),
)


def lever_family(agent: Agent | str) -> str:
    """The family an agent belongs to, so one lever is counted once."""
    name = (agent.agent if isinstance(agent, Agent) else agent).lower()
    for family, keys in _FAMILIES:
        if any(k in name for k in keys):
            return family
    return (agent.agent if isinstance(agent, Agent) else agent)


#: Families whose name describes *reducing* an exposure. Their agents are
#: stated either as the exposure ("High-fat diet") or as the reduction
#: ("Lower saturated-fat intake"); see :func:`family_helps`.
INVERSE_FAMILIES: Final[frozenset[str]] = frozenset({
    "Less saturated fat and animal-based eating",
    "Less red meat, choline and carnitine",
    "Less arginine-rich protein",
    "Review acid-suppressing medication with a prescriber",
})

#: Pure exposures. Never offered as a lever; shown as what pushes an
#: organism the wrong way.
AVOID_FAMILIES: Final[frozenset[str]] = frozenset({
    "Antibiotics",
    "Low-fibre / refined-carbohydrate eating",
    "Smoking and gum disease",
})

_REDUCTION_WORDS: Final = ("lower", "less ", "reduced", "reduce ", "restriction", "review", "plant-based")


def _stated_as_reduction(agent: Agent) -> bool:
    name = agent.agent.lower()
    return any(w in name for w in _REDUCTION_WORDS)


def family_helps(agent: Agent, *, helps_as_stated: bool) -> bool | None:
    """Whether doing what the *family* name says moves the organism the way
    the reader wants, given whether the agent as literally stated does.

    Returns None for an avoid-family: those are not levers.

    "High-fat diet" raises R. gnavus; for an overgrown R. gnavus that agent
    hurts as stated - and so the family "Less saturated fat" helps. "Lower
    saturated-fat intake" lowers E. ramosum; that helps as stated and the
    family helps too. The two statements of one lever must agree.
    """
    family = lever_family(agent)
    if family in AVOID_FAMILIES:
        return None
    if family in INVERSE_FAMILIES:
        return helps_as_stated if _stated_as_reduction(agent) else not helps_as_stated
    # Direct family: "Less fermented dairy" is a reduction of the family's
    # thing, so its stated effect is the family's effect inverted.
    if agent.agent.lower().startswith("less "):
        return not helps_as_stated
    return helps_as_stated


__all__ = [
    "AVOID_FAMILIES", "EVIDENCE", "INVERSE_FAMILIES", "KIND", "Agent", "Modulators",
    "available", "family_helps", "lever_family", "lookup",
]
