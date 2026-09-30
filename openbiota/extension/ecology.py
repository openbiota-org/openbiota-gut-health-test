"""A09 — companion diversity, dominance, ratios and aerotolerance.

BUILD_SPEC_v0.8.3 section 6.1-6.4. Everything here is computed from *one*
explicitly declared abundance vector, closed to sum one over a stated
inventory, and none of it replaces the existing Shannon, evenness, richness or
Firmicutes:Bacteroidota results. The point of these companions is that
`exp(H)` is a number of species a reader can picture, where an entropy in nats
is not.

Four rules the specification is emphatic about, and which the code enforces
rather than documents:

* **The vector is declared.** Detection depth, catalogue, abundance filter and
  the unresolved fraction travel with every value, because a diversity that
  does not say what it counted cannot be compared with anything.
* **A ratio is dimensionless.** Numerator, denominator and membership are
  published; a zero denominator is censored, not infinity and not zero.
* **No cutoff is invented to fill a gauge.** Every ratio here is descriptive.
  There is no healthy range for Prevotella:Bacteroides, and no action may
  propose raising *Fusobacterium* to move a ratio.
* **Unknown is not zero.** A trait nobody has measured for an organism reduces
  coverage; it does not move the organism into the other category.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.extension.schema import ExtensionMetric, fingerprint

METHOD_DIVERSITY: Final = "ext083.hill_diversity/1.0"
METHOD_RATIO: Final = "ext083.composition_ratio/1.0"
METHOD_AEROTOLERANCE: Final = "ext083.aerotolerance/1.0"
METHOD_REDUNDANCY: Final = "ext083.functional_redundancy/1.0"

#: The abundance filter below which an organism is not admitted to the closed
#: vector. Declared rather than implied: a different floor gives a different
#: richness, and a diversity compared across two floors is meaningless.
DEFAULT_ABUNDANCE_FLOOR_PERCENT: Final = 0.0


@dataclass(slots=True, frozen=True)
class AbundanceVector:
    """One declared, closed abundance vector and the facts about how it was made.

    `p` sums to one over `labels`. `unresolved_fraction` is the mass the lane
    could not name, kept outside `p` so it can be shown rather than silently
    redistributed over the named organisms.
    """

    labels: tuple[str, ...]
    p: tuple[float, ...]
    lane: str
    catalogue: str
    abundance_floor_percent: float
    unresolved_fraction: float | None = None
    detection_depth_fragments: int | None = None
    genus_of: Mapping[str, str] = field(default_factory=dict)
    phylum_of: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.labels) != len(self.p):
            raise ValueError("labels and abundances differ in length")
        if any(v < 0 for v in self.p):
            raise ValueError("a negative abundance cannot enter a closed vector")
        total = sum(self.p)
        if self.p and not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"vector is not closed to one (sums to {total!r})")

    @property
    def n_taxa(self) -> int:
        return len(self.labels)

    def declaration(self) -> dict[str, Any]:
        """What this vector counted. Travels with every derived value."""
        return {
            "lane": self.lane,
            "catalogue": self.catalogue,
            "abundance_floor_percent": self.abundance_floor_percent,
            "n_taxa_in_vector": self.n_taxa,
            "unresolved_fraction": self.unresolved_fraction,
            "detection_depth_fragments": self.detection_depth_fragments,
        }

    def fingerprint(self) -> str:
        return fingerprint(self.lane, self.catalogue, self.abundance_floor_percent,
                           self.labels, [round(v, 12) for v in self.p])

    def collapse(self, level: str) -> dict[str, float]:
        """Sum the vector to genus or phylum, using the declared mapping.

        A label with no mapping at that level is returned under its own name
        prefixed `unmapped:`, so a rollup always reconciles to one and the
        unmapped mass is visible instead of quietly missing.
        """
        table = {"genus": self.genus_of, "phylum": self.phylum_of}.get(level)
        if table is None:
            raise ValueError(f"unknown level {level!r}")
        out: dict[str, float] = {}
        for label, value in zip(self.labels, self.p, strict=True):
            key = table.get(label) or f"unmapped:{label}"
            out[key] = out.get(key, 0.0) + value
        return out


def vector_from_inventory(
    organisms: Sequence[Mapping[str, Any]],
    *,
    lane: str = "organism_inventory.primary",
    catalogue: str = "unified inventory (marker + whole-genome lanes)",
    floor_percent: float = DEFAULT_ABUNDANCE_FLOOR_PERCENT,
    unresolved_percent: float | None = None,
    depth_fragments: int | None = None,
) -> AbundanceVector:
    """Build the declared vector from the existing organism inventory.

    Reads the inventory rather than recomputing abundance, so these companion
    values describe the same community the rest of the report describes.
    """
    labels: list[str] = []
    raw: list[float] = []
    genus: dict[str, str] = {}
    phylum: dict[str, str] = {}
    for row in organisms:
        pct = row.get("percent")
        if pct is None or float(pct) <= floor_percent:
            continue
        name = str(row.get("species") or "").replace("_", " ").strip()
        if not name:
            continue
        labels.append(name)
        raw.append(float(pct))
        g = row.get("gtdb_genus") or row.get("genus")
        if g:
            genus[name] = str(g).replace("_", " ")
        ph = row.get("phylum") or row.get("gtdb_phylum")
        if ph:
            phylum[name] = str(ph)
    total = sum(raw)
    p = tuple(v / total for v in raw) if total > 0 else ()
    return AbundanceVector(
        labels=tuple(labels),
        p=p,
        lane=lane,
        catalogue=catalogue,
        abundance_floor_percent=floor_percent,
        unresolved_fraction=(None if unresolved_percent is None else unresolved_percent / 100.0),
        detection_depth_fragments=depth_fragments,
        genus_of=genus,
        phylum_of=phylum,
    )


#: Phylum names the GTDB rename split in two. Both spellings are one phylum,
#: and a ratio that matched only one of them would read as an absence.
PHYLUM_ALIASES: Final[Mapping[str, str]] = {
    "Firmicutes": "Bacillota",
    "Bacillota": "Bacillota",
    "Bacteroidetes": "Bacteroidota",
    "Bacteroidota": "Bacteroidota",
    "Actinobacteria": "Actinobacteriota",
    "Actinomycetota": "Actinobacteriota",
    "Actinobacteriota": "Actinobacteriota",
    "Proteobacteria": "Proteobacteria",
    "Pseudomonadota": "Proteobacteria",
    "Verrucomicrobia": "Verrucomicrobiota",
    "Verrucomicrobiota": "Verrucomicrobiota",
    "Fusobacteria": "Fusobacteriota",
    "Fusobacteriota": "Fusobacteriota",
}


def phylum_vector(
    phyla: Mapping[str, float],
    *,
    lane: str,
    catalogue: str,
) -> AbundanceVector:
    """A declared phylum-level vector from a lane's own phylum table.

    Built separately from the species vector because the species inventory
    carries no phylum lineage: a phylum ratio computed by guessing lineage
    from a genus name would be a different measurement wearing these units.
    Unclassified mass is held out as the unresolved fraction rather than
    spread over the named phyla.
    """
    named: dict[str, float] = {}
    unresolved = 0.0
    for raw, pct in phyla.items():
        value = float(pct)
        if value <= 0:
            continue
        key = str(raw).split(" (")[0].strip()
        if "unclassified" in key.lower() or "unknown" in key.lower():
            unresolved += value
            continue
        named[PHYLUM_ALIASES.get(key, key)] = named.get(PHYLUM_ALIASES.get(key, key), 0.0) + value
    total = sum(named.values())
    if total <= 0:
        return AbundanceVector((), (), lane, catalogue, 0.0, unresolved_fraction=None)
    labels = tuple(sorted(named))
    return AbundanceVector(
        labels=labels,
        p=tuple(named[k] / total for k in labels),
        lane=lane,
        catalogue=catalogue,
        abundance_floor_percent=0.0,
        unresolved_fraction=(unresolved / (total + unresolved) if unresolved else None),
    )


# --------------------------------------------------------------------------- #
# 6.1 companion diversity and dominance
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class Diversity:
    """The companion values of section 6.1, with their declaration."""

    shannon_nats: float | None
    effective_species: float | None
    inverse_simpson: float | None
    gini_simpson: float | None
    dominance_top1: float | None
    dominance_top5: float | None
    observed_taxa: int
    declaration: dict[str, Any]

    def to_json(self) -> dict[str, Any]:
        return {
            "shannon_nats": self.shannon_nats,
            "effective_shannon_species": self.effective_species,
            "inverse_simpson": self.inverse_simpson,
            "gini_simpson": self.gini_simpson,
            "dominance_top1": self.dominance_top1,
            "dominance_top5": self.dominance_top5,
            "observed_taxa": self.observed_taxa,
            "vector": self.declaration,
        }


def diversity(vector: AbundanceVector) -> Diversity:
    """Hill numbers and dominance from the declared vector.

    Natural logarithms throughout, stated in the unit, because the existing
    Shannon result may use a different base and the two must not be compared
    as though they were the same number. An empty vector returns nulls: a
    community with nothing in it has no diversity, not a diversity of zero.
    """
    p = [v for v in vector.p if v > 0.0]
    if not p:
        return Diversity(None, None, None, None, None, None, 0, vector.declaration())
    h = -sum(v * math.log(v) for v in p)
    simpson_sum = sum(v * v for v in p)
    ordered = sorted(vector.p, reverse=True)
    return Diversity(
        shannon_nats=h,
        effective_species=math.exp(h),
        inverse_simpson=(1.0 / simpson_sum if simpson_sum > 0 else None),
        gini_simpson=1.0 - simpson_sum,
        dominance_top1=ordered[0],
        dominance_top5=sum(ordered[:5]),
        observed_taxa=len(p),
        declaration=vector.declaration(),
    )


# --------------------------------------------------------------------------- #
# 6.1 descriptive composition ratios
# --------------------------------------------------------------------------- #

#: The ratios section 6.1 requires, with the level they are defined at. The
#: level matters: a genus-level Fusobacterium:Faecalibacterium is not the
#: species-level qPCR ratio of the colorectal literature (E18), and swapping
#: them is called out explicitly by AT040.
RATIO_DEFINITIONS: Final[tuple[dict[str, Any], ...]] = (
    {
        "ratio_id": "ext083.ecology.ratio.proteobacteria_actinobacteriota",
        "label": "Proteobacteria (Pseudomonadota) : Actinobacteriota",
        "level": "phylum",
        "numerator": ("Proteobacteria", "Pseudomonadota"),
        "denominator": ("Actinobacteriota", "Actinobacteria"),
        "source_ids": ("E01", "E02"),
        "note": "Descriptive. No established healthy cut-point in stool.",
    },
    {
        "ratio_id": "ext083.ecology.ratio.prevotella_bacteroides",
        "label": "Prevotella : Bacteroides",
        "level": "genus",
        "numerator": ("Prevotella", "Segatella"),
        "denominator": ("Bacteroides",),
        "source_ids": ("E17", "E02"),
        "note": (
            "Diet-response context only (E17). Genus membership follows the declared "
            "taxonomy: Segatella copri was Prevotella copri."
        ),
    },
    {
        "ratio_id": "ext083.ecology.ratio.fusobacterium_faecalibacterium",
        "label": "Fusobacterium : Faecalibacterium",
        "level": "genus",
        "numerator": ("Fusobacterium",),
        "denominator": ("Faecalibacterium",),
        "source_ids": ("E18",),
        "note": (
            "Genus level from shotgun abundance. Not the species-specific colorectal "
            "qPCR ratio (E18), and not a screening test."
        ),
    },
)


@dataclass(slots=True)
class Ratio:
    """One descriptive ratio, with everything needed to check it."""

    ratio_id: str
    label: str
    level: str
    value: float | None
    state: str
    numerator_value: float
    denominator_value: float
    numerator_members: tuple[str, ...]
    denominator_members: tuple[str, ...]
    source_ids: tuple[str, ...]
    note: str
    vector_lane: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "ratio_id": self.ratio_id,
            "label": self.label,
            "level": self.level,
            "value": self.value,
            "state": self.state,
            "unit": "dimensionless",
            "numerator_value": self.numerator_value,
            "denominator_value": self.denominator_value,
            "vector_lane": self.vector_lane,
            "numerator_members_found": list(self.numerator_members),
            "denominator_members_found": list(self.denominator_members),
            "source_ids": list(self.source_ids),
            "note": self.note,
            "direction": "descriptive",
        }


def ratios(
    vector: AbundanceVector,
    *,
    by_level: Mapping[str, AbundanceVector] | None = None,
) -> tuple[Ratio, ...]:
    """Every required descriptive ratio, each from one declared vector.

    Both sides of a ratio must come from the same vector, or the number is a
    comparison between two lanes wearing the units of one. A ratio defined at
    phylum level therefore uses a phylum vector when one is supplied, and is
    censored rather than estimated when the species vector carries no lineage
    at that level.

    A zero denominator yields `censored_zero_denominator` with a null value.
    Returning zero or infinity there would both be claims the data cannot
    make, and a gauge cannot be drawn from either.
    """
    out: list[Ratio] = []
    for spec in RATIO_DEFINITIONS:
        level = str(spec["level"])
        source = (by_level or {}).get(level, vector)
        table = (
            source.collapse(level) if source is vector
            else dict(zip(source.labels, source.p, strict=True))
        )
        num_found = tuple(sorted(k for k in table if k in spec["numerator"]))
        den_found = tuple(sorted(k for k in table if k in spec["denominator"]))
        num = sum(table[k] for k in num_found)
        den = sum(table[k] for k in den_found)
        if den <= 0.0:
            value, state = None, "censored_zero_denominator"
        elif num <= 0.0:
            value, state = 0.0, "measured_numerator_absent"
        else:
            value, state = num / den, "measured"
        out.append(Ratio(
            ratio_id=str(spec["ratio_id"]),
            label=str(spec["label"]),
            level=level,
            vector_lane=source.lane,
            value=value,
            state=state,
            numerator_value=num,
            denominator_value=den,
            numerator_members=num_found,
            denominator_members=den_found,
            source_ids=tuple(spec["source_ids"]),
            note=str(spec["note"]),
        ))
    return tuple(out)


# --------------------------------------------------------------------------- #
# 6.3 aerotolerance balance
# --------------------------------------------------------------------------- #

#: Oxygen phenotype vocabulary. `unresolved` is a real answer and the most
#: common one; it reduces coverage rather than choosing a side.
OXYGEN_PHENOTYPES: Final[frozenset[str]] = frozenset({
    "strict_anaerobe", "aerotolerant", "unresolved",
})

#: How the phenotype was established, strongest first. Section 6.3: cultured
#: phenotypes from BacDive, then explicitly labelled genomic predictions.
TRAIT_EVIDENCE: Final[tuple[str, ...]] = (
    "cultured_strain", "cultured_species", "genomic_prediction", "genus_generalisation",
)

#: GTDB splits many gut genera with a letter suffix - Blautia_A, Ruminococcus_E,
#: Eubacterium_F. The suffix marks a genus boundary, not a different oxygen
#: requirement, so a trait lookup strips it. An explicitly listed split genus
#: still wins, because the strip is only a fallback.
_GTDB_SUFFIX: Final = re.compile(r"[ _][A-Z]{1,2}$")

TRAITS_FILE: Final = Path(__file__).resolve().parents[2] / "extension" / "oxygen_traits.yaml"


def _base_genus(name: str) -> str:
    """A GTDB genus with its split suffix removed."""
    return _GTDB_SUFFIX.sub("", str(name).strip())


def load_oxygen_traits(path: Path | None = None) -> dict[str, dict[str, str]]:
    """The curated oxygen phenotype table, keyed by species and by genus.

    Species records are kept separate from genus records so `aerotolerance`
    can prefer the more specific one and label the other a generalisation.
    """
    raw = yaml.safe_load((path or TRAITS_FILE).read_text(encoding="utf-8")) or {}
    default_evidence = str(raw.get("default_evidence") or "genus_generalisation")
    default_source = str(raw.get("default_source_id") or "E06")
    out: dict[str, dict[str, str]] = {}
    for phenotype in ("strict_anaerobe", "aerotolerant"):
        for genus in raw.get(phenotype) or ():
            out[str(genus)] = {
                "phenotype": phenotype,
                "evidence": default_evidence,
                "source_id": default_source,
            }
    for species, record in (raw.get("species") or {}).items():
        entry = {
            "phenotype": str(record.get("phenotype") or "unresolved"),
            "evidence": str(record.get("evidence") or "cultured_species"),
            "source_id": str(record.get("source_id") or default_source),
        }
        if record.get("note"):
            entry["note"] = str(record["note"])
        out[str(species)] = entry
    return out


@dataclass(slots=True)
class Aerotolerance:
    """The section 6.3 values, with coverage and the contributor table."""

    aerotolerant_fraction: float | None
    mapi: float | None
    mapi_state: str
    trait_coverage: float | None
    aerotolerant_mass: float
    anaerobe_mass: float
    unresolved_mass: float
    contributors: tuple[dict[str, Any], ...]
    declaration: dict[str, Any]

    def to_json(self) -> dict[str, Any]:
        return {
            "aerotolerant_fraction": self.aerotolerant_fraction,
            "mapi": self.mapi,
            "mapi_state": self.mapi_state,
            "trait_coverage": self.trait_coverage,
            "aerotolerant_mass": self.aerotolerant_mass,
            "strict_anaerobe_mass": self.anaerobe_mass,
            "unresolved_mass": self.unresolved_mass,
            "contributors": list(self.contributors),
            "vector": self.declaration,
            "source_ids": ["E04", "E05", "E06", "E07"],
            "limitations": [
                "An ecological proxy, not a measurement of oxygen in the gut and not a "
                "permeability test.",
                "Organisms with no established oxygen phenotype reduce coverage; they are "
                "not counted as anaerobes.",
            ],
        }


def aerotolerance(
    vector: AbundanceVector,
    traits: Mapping[str, Mapping[str, str]],
) -> Aerotolerance:
    """Aerotolerant fraction, MAPI and trait coverage.

    ``aerotolerant_fraction = A / (A + N)``
    ``trait_coverage = (A + N) / (A + N + U)``
    ``MAPI = ln(A / N)``, only when both are positive.

    With one component zero the fraction is still computable and MAPI is
    censored: section 6.3 forbids inserting a pseudocount, because the log of
    an absent group is not a small number, it is undefined.
    """
    a = n = u = 0.0
    rows: list[dict[str, Any]] = []
    for label, mass in zip(vector.labels, vector.p, strict=True):
        record = traits.get(label)
        if record is None:
            # Genus candidates, most specific first. The GTDB genus and the
            # genus in the organism's own name can disagree after a
            # reclassification - *Pseudoflavonifractor capillosus* sits in
            # GTDB's *Enterenecus* - and an oxygen requirement established
            # under either name still describes the same organism. The split
            # suffix (Blautia_A) is stripped last, as a fallback.
            candidates = [
                vector.genus_of.get(label),
                label.split(" ")[0],
                _base_genus(vector.genus_of.get(label) or ""),
                _base_genus(label.split(" ")[0]),
            ]
            for candidate in candidates:
                if candidate and candidate in traits:
                    record = dict(traits[candidate]) | {"evidence": "genus_generalisation"}
                    break
        phenotype = str((record or {}).get("phenotype") or "unresolved")
        if phenotype not in OXYGEN_PHENOTYPES:
            phenotype = "unresolved"
        if phenotype == "aerotolerant":
            a += mass
        elif phenotype == "strict_anaerobe":
            n += mass
        else:
            u += mass
        rows.append({
            "taxon": label,
            "abundance_fraction": mass,
            "phenotype": phenotype,
            "trait_evidence": (record or {}).get("evidence"),
            "trait_source_id": (record or {}).get("source_id"),
        })

    assessed = a + n
    if assessed <= 0.0:
        return Aerotolerance(
            aerotolerant_fraction=None, mapi=None, mapi_state="no_supported_traits",
            trait_coverage=(0.0 if (assessed + u) > 0 else None),
            aerotolerant_mass=a, anaerobe_mass=n, unresolved_mass=u,
            contributors=tuple(sorted(rows, key=lambda r: -r["abundance_fraction"])),
            declaration=vector.declaration(),
        )
    if a > 0.0 and n > 0.0:
        mapi, mapi_state = math.log(a / n), "measured"
    else:
        mapi, mapi_state = None, "censored_component_absent"
    return Aerotolerance(
        aerotolerant_fraction=a / assessed,
        mapi=mapi,
        mapi_state=mapi_state,
        trait_coverage=assessed / (assessed + u) if (assessed + u) > 0 else None,
        aerotolerant_mass=a,
        anaerobe_mass=n,
        unresolved_mass=u,
        contributors=tuple(sorted(rows, key=lambda r: -r["abundance_fraction"])),
        declaration=vector.declaration(),
    )


# --------------------------------------------------------------------------- #
# 6.4 functional redundancy
# --------------------------------------------------------------------------- #


def effective_carriers(weights: Iterable[float]) -> float | None:
    """``exp(-sum(w ln w))`` over carrier contributions normalised within one
    function.

    Section 6.4: redundancy among *resolved* contributors. Two carriers at
    50/50 give 2.0; two at 99/1 give barely more than 1, which is the point -
    a function held by one dominant organism plus a trace is not redundant.
    Returns None for an empty set rather than 0, because no resolved carrier
    is not the same as no redundancy.
    """
    values = [float(w) for w in weights if float(w) > 0.0]
    if not values:
        return None
    total = sum(values)
    if total <= 0.0:
        return None
    w = [v / total for v in values]
    return math.exp(-sum(x * math.log(x) for x in w))


@dataclass(slots=True)
class Redundancy:
    """One route's redundancy among the carriers actually resolved."""

    function_id: str
    label: str
    n_resolved_carriers: int
    effective_carriers: float | None
    carrier_coverage: float | None
    unassigned_fraction: float | None
    carriers: tuple[dict[str, Any], ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "function_id": self.function_id,
            "label": self.label,
            "n_resolved_carriers": self.n_resolved_carriers,
            "effective_carriers": self.effective_carriers,
            "carrier_coverage": self.carrier_coverage,
            "community_only_unassigned_fraction": self.unassigned_fraction,
            "carriers": list(self.carriers),
            "limitations": [
                "Redundancy among the carriers this assay resolved, not a guarantee of "
                "ecological resilience.",
            ],
        }


def redundancy(
    function_id: str,
    label: str,
    carriers: Sequence[Mapping[str, Any]],
    *,
    unassigned_fraction: float | None = None,
) -> Redundancy:
    """Redundancy for one function from its resolved carrier rows.

    Only rows whose assignment evidence supports carriage are counted, which
    is why `contributing_taxa` carries that field at all: an organism merely
    co-detected with a gene is context, not a second carrier.
    """
    resolved = [c for c in carriers if c.get("is_carrier_evidence")]
    weights = [float(c.get("value") or 0.0) for c in resolved]
    total_all = sum(float(c.get("value") or 0.0) for c in carriers) or None
    return Redundancy(
        function_id=function_id,
        label=label,
        n_resolved_carriers=len(resolved),
        effective_carriers=effective_carriers(weights),
        carrier_coverage=(sum(weights) / total_all if total_all else None),
        unassigned_fraction=unassigned_fraction,
        carriers=tuple(resolved),
    )


# --------------------------------------------------------------------------- #
# metric emission
# --------------------------------------------------------------------------- #


def metrics(
    vector: AbundanceVector,
    div: Diversity,
    ratio_rows: Sequence[Ratio],
    aero: Aerotolerance | None,
) -> tuple[ExtensionMetric, ...]:
    """Emit the A09 values as extension metrics.

    Every one is `descriptive`: none of these has an outcome-linked direction,
    and inventing one to colour a gauge is what section 6.1 forbids.
    """
    fp = vector.fingerprint()
    decl = vector.declaration()
    out: list[ExtensionMetric] = []

    def add(metric_id: str, label: str, value: float | None, unit: str, *,
            view_ids: tuple[str, ...] = (), limits: tuple[str, ...] = (),
            extra: dict[str, Any] | None = None, kind: str = "community_composition",
            denominator: str | None = "closed_named_abundance_vector") -> None:
        out.append(ExtensionMetric(
            metric_id=metric_id,
            label=label,
            kind=kind,
            state="measured" if value is not None else "insufficient_coverage",
            value=value,
            unit=unit if value is not None else None,
            denominator=denominator if value is not None else None,
            direction="descriptive",
            method_id=METHOD_DIVERSITY,
            input_fingerprint=fp,
            group="ecology",
            feature_id="A09",
            view_ids=view_ids,
            analytical_confidence="supported" if value is not None else "unresolved",
            evidence_maturity="not_applicable",
            source_ids=("E01", "E02", "E03"),
            limitations=limits or (
                "Descriptive. Requires a compatible assay and processing fingerprint "
                "before a change between samples is a biological change.",
            ),
            extra={"vector": decl} | (extra or {}),
        ))

    add("ext083.ecology.effective_shannon_species", "Effective number of species",
        div.effective_species, "species_equivalents",
        limits=("exp(Shannon) in natural logs. Not the existing Shannon index in a "
                "different base, and not comparable across detection floors.",))
    add("ext083.ecology.inverse_simpson", "Inverse Simpson diversity",
        div.inverse_simpson, "species_equivalents")
    add("ext083.ecology.gini_simpson", "Gini-Simpson index", div.gini_simpson, "fraction")
    add("ext083.ecology.dominance_top1", "Share of the single most abundant organism",
        div.dominance_top1, "fraction", view_ids=("M089",))
    add("ext083.ecology.dominance_top5", "Share of the five most abundant organisms",
        div.dominance_top5, "fraction", view_ids=("M089",))

    for r in ratio_rows:
        view = {"ext083.ecology.ratio.proteobacteria_actinobacteriota": ("M072",)}.get(r.ratio_id, ())
        # A numerator that was searched for and not found does not make the
        # ratio zero: it bounds it below the detection floor. Recording 0.0
        # would state a measurement nobody made, and the report's own card
        # already says "below the floor".
        if r.state == "measured_numerator_absent":
            metric_state, metric_value = "not_detected_above_assay_threshold", None
        elif r.value is None:
            metric_state, metric_value = "insufficient_coverage", None
        else:
            metric_state, metric_value = "measured", r.value
        out.append(ExtensionMetric(
            metric_id=r.ratio_id,
            label=r.label,
            kind="community_composition",
            state=metric_state,
            value=metric_value,
            unit="ratio" if metric_value is not None else None,
            denominator=(f"{r.level}-level abundance of {', '.join(r.denominator_members) or 'none found'}"
                         if metric_value is not None else None),
            direction="descriptive",
            method_id=METHOD_RATIO,
            input_fingerprint=fp,
            group="ecology",
            feature_id="A09",
            view_ids=view,
            analytical_confidence="supported" if metric_value is not None else "unresolved",
            source_ids=r.source_ids,
            limitations=(
                r.note,
                "A ratio is dimensionless; it is not a share of the community.",
            ),
            extra=r.to_json() | {"vector": decl},
        ))

    if aero is not None:
        out.append(ExtensionMetric(
            metric_id="ext083.ecology.aerotolerant_fraction",
            label="Aerotolerance balance",
            kind="community_composition",
            state="measured" if aero.aerotolerant_fraction is not None else "insufficient_coverage",
            value=aero.aerotolerant_fraction,
            unit="fraction" if aero.aerotolerant_fraction is not None else None,
            denominator="abundance with an established oxygen phenotype" if aero.aerotolerant_fraction is not None else None,
            direction="descriptive",
            method_id=METHOD_AEROTOLERANCE,
            input_fingerprint=fp,
            group="ecology",
            feature_id="A09",
            view_ids=("M119",),
            assessable_fraction=aero.trait_coverage,
            analytical_confidence="supported" if (aero.trait_coverage or 0) >= 0.5 else "provisional",
            source_ids=("E04", "E05", "E06", "E07"),
            limitations=tuple(aero.to_json()["limitations"]),
            extra=aero.to_json(),
        ))
        out.append(ExtensionMetric(
            metric_id="ext083.ecology.mapi",
            label="Microbial aerotolerance profile index (MAPI)",
            kind="experimental_index",
            state="measured" if aero.mapi is not None else "insufficient_coverage",
            value=aero.mapi,
            unit="ln_ratio" if aero.mapi is not None else None,
            denominator=None,
            direction="descriptive",
            method_id=METHOD_AEROTOLERANCE,
            input_fingerprint=fp,
            group="ecology",
            feature_id="A09",
            view_ids=("M119",),
            assessable_fraction=aero.trait_coverage,
            analytical_confidence="provisional",
            evidence_maturity="association_only",
            source_ids=("E04", "E05"),
            limitations=(
                "ln(aerotolerant / strict anaerobe). Censored when either group is absent: "
                "no pseudocount is inserted.",
                "A research ecological association, not a measurement of gut oxygen.",
            ),
            extra={"mapi_state": aero.mapi_state, "vector": decl},
        ))
    return tuple(out)


__all__ = [
    "DEFAULT_ABUNDANCE_FLOOR_PERCENT",
    "PHYLUM_ALIASES",
    "TRAITS_FILE",
    "OXYGEN_PHENOTYPES",
    "RATIO_DEFINITIONS",
    "TRAIT_EVIDENCE",
    "AbundanceVector",
    "Aerotolerance",
    "Diversity",
    "Ratio",
    "Redundancy",
    "aerotolerance",
    "diversity",
    "effective_carriers",
    "load_oxygen_traits",
    "metrics",
    "phylum_vector",
    "ratios",
    "redundancy",
    "vector_from_inventory",
]
