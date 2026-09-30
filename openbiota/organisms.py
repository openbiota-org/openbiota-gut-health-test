"""One verdict per organism: what it is, which way it cuts, and whether its level is a concern.

Every organism in the inventory gets a :class:`Verdict`: a four-way class in
this report's own vocabulary, a plain-language description, and a flag that
says whether the *amount* found is the kind that matters. That is the step a
category alone does not take. "Opportunist" is a property of the species;
"opportunist at the 97th percentile" is a finding about this sample.

Vocabulary
----------
Four classes. Not borrowed from any other report's labels:

``beneficial``
    Higher levels are associated with better outcomes in the human literature.
``opportunist``
    Normally a low-level resident; expansion is what the literature associates
    with inflammation, infection or a disturbed community. Not a pathogen
    label - a pathogen is a specific organism doing a specific thing, and the
    pathogen screen is where that is decided.
``conditional``
    Cuts both ways: helpful at typical levels or in the right company, and
    associated with problems when expanded or when the community around it is
    disturbed. Most of the gut is this, and calling it "variable" would say
    less than is known.
``unknown``
    Not enough human evidence to say. Said as such rather than guessed.

Resolution order
----------------
1. The species registry (``taxa/interpretations/``), under the organism's
   current name, its former name, or its GTDB name. Cited, per-species.
2. The guild files (``taxa/*.yaml``) where a guild implies a class.
3. The genus table (``taxa/genera.yaml``), marked as genus-level.
4. ``unknown``.

The flag
--------
Set from the percentile against reference adults, in the direction that
matters for the class: high for an opportunist, low for a beneficial
organism, high for a conditional one. An organism with no percentile - one
measured only on the catalogue that has no reference cohort - gets no flag,
because there is nothing to compare it against. That is stated, not hidden.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.inventory import Organism, canonical, former_name

BENEFICIAL: Final = "beneficial"
OPPORTUNIST: Final = "opportunist"
CONDITIONAL: Final = "conditional"
UNKNOWN: Final = "unknown"

CLASSES: Final = (BENEFICIAL, OPPORTUNIST, CONDITIONAL, UNKNOWN)

#: Reader-facing label for each class.
CLASS_LABEL: Final[dict[str, str]] = {
    BENEFICIAL: "Beneficial",
    OPPORTUNIST: "Opportunist",
    CONDITIONAL: "Conditional",
    UNKNOWN: "Unknown",
}

#: One sentence on what each class means, for legends.
CLASS_MEANING: Final[dict[str, str]] = {
    BENEFICIAL: (
        "Higher levels are associated with better outcomes in human studies. Butyrate "
        "producers, mucus-layer keepers and the classic probiotic genera live here."
    ),
    OPPORTUNIST: (
        "Normal at low levels. Expansion is what the research associates with inflammation, "
        "infection or a disturbed community. Not a pathogen label \u2014 that is the pathogen "
        "screen's decision."
    ),
    CONDITIONAL: (
        "Cuts both ways: part of a healthy gut at typical levels, associated with problems when "
        "expanded or when the community around it is disturbed. Most of the gut is this."
    ),
    UNKNOWN: (
        "Not enough human evidence to say. Many are known only from genomes and have never "
        "been grown in a laboratory."
    ),
}

#: Registry polarity -> class.
_POLARITY_TO_CLASS: Final[dict[str, str]] = {
    "health_associated": BENEFICIAL,
    "context_dependent": CONDITIONAL,
    "potential_pathobiont": OPPORTUNIST,
    "disease_associated": OPPORTUNIST,
    "validated_pathogen_signal": OPPORTUNIST,
}

#: Guild file -> class it implies, for species the registry does not name.
_GUILD_TO_CLASS: Final[dict[str, str]] = {
    "beneficial": BENEFICIAL,
    "bifidobacteria": BENEFICIAL,
    "butyrate_producers": BENEFICIAL,
    "probiotics": BENEFICIAL,
    "polyphenol_metabolisers": BENEFICIAL,
    "opportunistic_pathogens": OPPORTUNIST,
    "hexa_lps": OPPORTUNIST,
    "sulfate_reducers": OPPORTUNIST,
    "oral_origin": CONDITIONAL,
    "mucin_degraders": CONDITIONAL,
    "methanogens": CONDITIONAL,
}

#: Percentile at or above which an expansion is flagged.
HIGH_AT: Final = 90.0
#: Percentile at or below which a depletion is flagged.
LOW_AT: Final = 10.0
#: Stronger thresholds, for the "notably" tier.
VERY_HIGH_AT: Final = 97.0
VERY_LOW_AT: Final = 3.0

#: How strong a flag is: 0 none, 1 watch, 2 flag.
FLAG_NONE: Final = 0
FLAG_WATCH: Final = 1
FLAG_ISSUE: Final = 2


@dataclass(frozen=True, slots=True)
class Verdict:
    """What the report says about one organism."""

    organism: Organism
    cls: str
    description: str
    #: "species", "guild", "genus" or "none": where the class came from.
    basis: str
    #: The direction that would be a concern for this class.
    concern_when: str
    #: "high", "low" or "" - the direction actually found, if flagged.
    flag: str
    #: 0 none, 1 watch, 2 issue.
    flag_level: int
    #: Short reason for the flag, for rendering.
    flag_reason: str

    @property
    def label(self) -> str:
        return CLASS_LABEL[self.cls]

    @property
    def flagged(self) -> bool:
        return self.flag_level > FLAG_NONE

    @property
    def is_issue(self) -> bool:
        return self.flag_level >= FLAG_ISSUE

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.organism.species,
            "display": self.organism.display,
            "class": self.cls,
            "basis": self.basis,
            "flag": self.flag,
            "flag_level": self.flag_level,
            "flag_reason": self.flag_reason,
            "percent": round(self.organism.percent, 4),
            "in_primary": self.organism.in_primary,
            "secondary_percent": self.organism.secondary_percent,
            "percentile": self.organism.percentile,
        }


# --------------------------------------------------------------------------- #
# registries
# --------------------------------------------------------------------------- #


@lru_cache(maxsize=1)
def _species_registry() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for path in sorted(Path("taxa/interpretations").glob("*.yaml")):
        data = yaml.safe_load(path.read_text()) or {}
        for row in data.get("species_interpretations") or []:
            if isinstance(row, dict) and row.get("species"):
                out[str(row["species"])] = row
    return out


@lru_cache(maxsize=1)
def _guild_membership() -> dict[str, str]:
    """Species -> the class its strongest guild implies."""
    out: dict[str, str] = {}
    # Beneficial and opportunist guilds are decisive; conditional ones only
    # fill in where nothing stronger has spoken.
    order = [BENEFICIAL, OPPORTUNIST, CONDITIONAL]
    for guild, cls in sorted(_GUILD_TO_CLASS.items(), key=lambda kv: order.index(kv[1])):
        path = Path("taxa") / f"{guild}.yaml"
        if not path.is_file():
            continue
        data = yaml.safe_load(path.read_text()) or {}
        members = data.get("members") or data.get("species") or data.get("taxa") or []
        for m in members:
            name = m.get("species") if isinstance(m, dict) else m
            if not name:
                continue
            key = str(name).strip().replace(" ", "_")
            if key not in out:
                out[key] = cls
    return out


@lru_cache(maxsize=1)
def _genus_table() -> dict[str, dict[str, Any]]:
    path = Path("taxa/classification/genera.yaml")
    if not path.is_file():
        return {}
    data = yaml.safe_load(path.read_text()) or {}
    return {str(r["genus"]): r for r in data.get("genera") or [] if r.get("genus")}


def _genus_key(name: str) -> str:
    """Genus as the tables key it: GTDB's ``Blautia_A`` collapses to ``Blautia``."""
    return name.split(" ", 1)[0].split("_", 1)[0]


# --------------------------------------------------------------------------- #
# resolution
# --------------------------------------------------------------------------- #


def _registry_row(o: Organism) -> dict[str, Any] | None:
    reg = _species_registry()
    for key in (o.species, canonical(o.species), former_name(o.species) or "",
                (o.gtdb or "").replace(" ", "_")):
        if key and key in reg:
            return reg[key]
    return None


def _guild_class(o: Organism) -> str | None:
    g = _guild_membership()
    for key in (o.species, canonical(o.species), former_name(o.species) or ""):
        if key and key in g:
            return g[key]
    return None


def _genus_row(o: Organism) -> dict[str, Any] | None:
    table = _genus_table()
    for candidate in (o.genus, o.gtdb_genus or "", o.display):
        if not candidate:
            continue
        key = _genus_key(candidate)
        if key in table:
            return table[key]
    return None


def _family_row(o: Organism) -> dict[str, Any] | None:
    """The family entry, for a genus the table does not know.

    Dozens of gut genera were described after 2019 from assembled genomes
    - Dysosmobacter, Faecimonas, Frisingicoccus - and almost everything
    known about them is what is known about their family. GTDB places every
    bin in a family, so the family entry is the honest level to speak at.
    """
    from openbiota import gtdb

    hit = gtdb.lookup(o.sgb) if o.sgb else None
    if hit is None:
        return None
    table = _genus_table()
    fam = hit.family
    return table.get(fam) if fam else None


def _flag(cls: str, concern_when: str, pct: float | None) -> tuple[str, int, str]:
    """Direction, level and reason, from class and percentile."""
    if pct is None or concern_when in ("", "none"):
        return "", FLAG_NONE, ""
    if concern_when in ("high", "both") and pct >= HIGH_AT:
        level = FLAG_ISSUE if (pct >= VERY_HIGH_AT or cls == OPPORTUNIST) else FLAG_WATCH
        # A beneficial organism that is abundant is, at most, worth a look:
        # the literature that raises it (very high Akkermansia with low
        # fibre, say) is real but conditional, and it must never sit in the
        # same tier as an opportunist bloom.
        if cls == BENEFICIAL:
            level = FLAG_WATCH
        word = "notably high" if pct >= VERY_HIGH_AT else "high"
        return "high", level, f"{word} \u2014 {_ordinal(pct)} percentile among reference adults who carry it"
    if concern_when in ("low", "both") and pct <= LOW_AT:
        level = FLAG_ISSUE if pct <= VERY_LOW_AT else FLAG_WATCH
        word = "notably low" if pct <= VERY_LOW_AT else "low"
        return "low", level, f"{word} \u2014 {_ordinal(pct)} percentile among reference adults who carry it"
    return "", FLAG_NONE, ""


def _ordinal(p: float) -> str:
    """A percentile as words: 99th, or >99th at the ceiling."""
    if p >= 99.5:
        return ">99th"
    if p < 0.5:
        return "<1st"
    n = int(round(p))
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def verdict(o: Organism) -> Verdict:
    """Classify one organism and judge its level."""
    row = _registry_row(o)
    grow = _genus_row(o)
    if row is not None:
        cls = _POLARITY_TO_CLASS.get(str(row.get("polarity")), CONDITIONAL)
        desc = " ".join(str(row.get("summary") or "").split())
        basis = "species"
        tags = {str(t) for t in (row.get("tags") or [])}
        # The registry files probiotics as context-dependent because they
        # are transient - present when a supplement is being taken, gone
        # when it stops. That is true and worth saying in the prose, but
        # as a class it misleads: a reader sees the best-known probiotic
        # labelled the same way as a mucus degrader. Where the genus table
        # is decisive, the genus class stands and the registry's prose
        # still explains the nuance.
        if cls == CONDITIONAL and grow is not None and grow.get("class") in (BENEFICIAL, OPPORTUNIST):
            cls = str(grow["class"])
        concern = _concern_from_assertions(row, cls, tags)
    else:
        gcls = _guild_class(o)
        frow = _family_row(o) if grow is None else None
        if gcls is not None:
            cls = gcls
            basis = "guild"
            desc = " ".join(str((grow or frow or {}).get("description") or "").split())
            concern = {BENEFICIAL: "low", OPPORTUNIST: "high", CONDITIONAL: "high"}.get(cls, "none")
        elif grow is not None:
            cls = str(grow.get("class") or UNKNOWN)
            basis = "genus"
            desc = " ".join(str(grow.get("description") or "").split())
            concern = str(grow.get("concern_when") or "none")
        elif frow is not None:
            cls = str(frow.get("class") or UNKNOWN)
            basis = "family"
            desc = " ".join(str(frow.get("description") or "").split())
            concern = str(frow.get("concern_when") or "none")
        else:
            cls, basis, desc, concern = UNKNOWN, "none", "", "none"
    if not desc:
        desc = _fallback_description(o)
    flag, level, reason = _flag(cls, concern, getattr(o, "level_percentile", o.percentile))
    # A population counted within a relative's share has no share of its own
    # to be high or low: the composition carries it under the relative, and
    # that relative's level is the one read. Flagging it as well would list
    # one population twice, under two catalogues' names (Phocaeicola vulgatus
    # at the 93rd percentile beside "Phocaeicola SPECIV4_34405" at the 92nd).
    # Its own lane's percentile stays in the full table, marked as its lane's.
    if flag and not o.in_primary and getattr(o, "share_basis", "none") in ("member", "none"):
        flag, level, reason = "", FLAG_NONE, ""
    # An opportunist that few reference adults carry at all is worth a look
    # whatever its level: the level percentile ranks it among carriers, so
    # rarity of carriage is a separate fact, said as such.
    # Only the deep reference population can say how common carriage is: a
    # lane cohort profiled on read prefixes has a detection floor, and "4 of
    # 100 carry it" there is a floor, not a prevalence.
    deep_reference = (getattr(o, "percentile_source", None) or "scoring cohort") == "scoring cohort"
    if not flag and cls == OPPORTUNIST and deep_reference and o.prevalence is not None and o.prevalence < 0.10 \
            and o.in_primary and o.percent > 0:
        flag, level = "uncommon", FLAG_WATCH
        reason = f"uncommon \u2014 carried by {o.prevalence:.0%} of reference adults"
    return Verdict(
        organism=o, cls=cls, description=desc, basis=basis,
        concern_when=concern, flag=flag, flag_level=level, flag_reason=reason,
    )


def _concern_from_assertions(row: dict[str, Any], cls: str, tags: set[str]) -> str:
    """Which direction the cited literature actually reports as the problem.

    Each registry assertion records whether the finding was about a higher
    or a lower level. That is a better guide than the class default: an
    organism whose every citation is about depletion should not be flagged
    for being abundant. Probiotic genera are never flagged high - a
    supplement being taken is not an overgrowth.
    """
    if "probiotic_genus" in tags:
        return "low" if cls == BENEFICIAL else "none"
    dirs = {str(a.get("direction") or "") for a in (row.get("assertions") or [])}
    high = "higher" in dirs
    low = "lower" in dirs
    if high and not low:
        return "high"
    if low and not high:
        return "low"
    if high and low:
        return "both"
    return {BENEFICIAL: "low", OPPORTUNIST: "high", CONDITIONAL: "high"}.get(cls, "none")


def _fallback_description(o: Organism) -> str:
    """A sentence for an organism nothing describes, honest about that."""
    genus = (o.gtdb_genus or o.genus or "").split("_")[0]
    if o.unnamed:
        where = f" of the genus <i>{genus}</i>" if genus else ""
        return (
            f"An organism{where} known only from genomes assembled out of stool samples. "
            "It has never been grown in a laboratory, so its physiology and any effect on "
            "health are not yet described."
        )
    if genus:
        return (
            f"A member of the genus <i>{genus}</i> for which no human evidence has been "
            "curated. Its presence is recorded; nothing here says it is helpful or harmful."
        )
    return "No curated evidence exists for this organism; its presence is recorded only."


def verdicts(organisms: list[Organism]) -> list[Verdict]:
    return [verdict(o) for o in organisms]


def composition(vs: list[Verdict], *, unclassified: float | None = None) -> dict[str, float]:
    """Each class's share of the primary composition, plus the unplaced remainder.

    Only organisms measured by the primary lane contribute; a secondary-only
    detection has no share to contribute and is not invented one. When the
    primary lane reported sequence it could not place, that is returned
    under ``"unclassified"`` so the parts add to the whole.
    """
    out = dict.fromkeys(CLASSES, 0.0)
    for v in vs:
        if v.organism.in_primary:
            out[v.cls] += v.organism.percent
    placed = sum(out.values())
    if unclassified is not None and unclassified > 0:
        out["unclassified"] = float(unclassified)
        total = placed + float(unclassified)
    else:
        total = placed or 1.0
    return {c: 100.0 * s / total for c, s in out.items()}


def importance(v: Verdict) -> float:
    """How much one flagged organism should weigh against another.

    Distance from the median alone gets this wrong: it ranks a trace
    organism at the 99.6th percentile above a 13% bloom at the 98th. The
    bloom is the finding. So the weight is the percentile excess times the
    log of the abundance, with an unambiguous direction (an opportunist
    high, a beneficial organism low) worth more than a conditional one
    expanded.
    """
    import math

    p = v.organism.percentile
    if p is None or not v.flagged:
        return 0.0
    excess = (p - HIGH_AT) if v.flag == "high" else (LOW_AT - p)
    excess = max(excess, 0.5)
    # Abundance on a log scale, floored so a trace organism still counts:
    # 0.01% -> 1, 0.1% -> 2, 1% -> 3, 10% -> 4. An organism seen only by a
    # secondary lane weighs by that lane's reading.
    mass = math.log10(max(v.organism.best_percent, 0.001)) + 4.0
    direction = 1.6 if (v.cls == OPPORTUNIST and v.flag == "high") or (
        v.cls == BENEFICIAL and v.flag == "low") else 1.0
    return v.flag_level * direction * excess * mass


def issues(vs: list[Verdict]) -> list[Verdict]:
    """Flagged organisms, most consequential first."""
    return sorted((v for v in vs if v.flagged), key=lambda v: -importance(v))


__all__ = [
    "BENEFICIAL",
    "CLASSES",
    "CLASS_LABEL",
    "CLASS_MEANING",
    "CONDITIONAL",
    "FLAG_ISSUE",
    "FLAG_NONE",
    "FLAG_WATCH",
    "OPPORTUNIST",
    "UNKNOWN",
    "Verdict",
    "composition",
    "importance",
    "issues",
    "verdict",
    "verdicts",
]
