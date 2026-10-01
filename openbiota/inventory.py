"""One organism inventory, assembled from every profiler that ran.

Why this module exists
----------------------
The pipeline runs more than one taxonomic profiler over the same host-filtered
reads. They do not see the same thing. On a typical sample the older catalogue
names about 92 species and the newer one about 159, because the newer one is
built from roughly a million genomes against the older one's seventeen
thousand, and it carries organisms that have never been cultured.

Reporting only one of them throws away real detections. Reporting them as two
separate tables in two separate sections - which is what this report used to do
- makes the reader reconcile them by hand, and invites the reasonable question
of why the older answer is the one on the page.

So detection is unified here, once, and every section that lists organisms
reads the result. The reader never sees a profiler name or a version number;
they see the organisms found in their sample.

What may and may not be merged
------------------------------
A hard line runs through this module, and it is not a matter of taste:

* **Presence and abundance may be merged.** "Is this organism here, and how
  much of it" is answered by the sample alone. Two detectors are strictly
  better than one, and the union is the honest answer.

* **Percentiles may not be merged.** A percentile is a position against a
  reference population, and the only reference population that exists
  (curatedMetagenomicData, 3,027 adults) was processed on one specific
  catalogue. Its maintainers state the catalogues are not directly comparable.
  Ranking a measurement from one catalogue against a cohort built on another
  produces a number that looks authoritative and means nothing.

So an organism seen only by the newer catalogue is *listed, counted, named and
strain-typed*, but its percentile column reads "not in the reference set"
rather than inventing a rank. That is a limit of the available reference data,
stated plainly, and it is the only place the split survives.

Adding a profiler later
-----------------------
Add one :class:`Lane` to the list built by :func:`lanes_from_results` and stop.
Nothing else in the report needs to change: every consumer reads
:class:`Inventory`, not a profiler. Set ``rankable=True`` only if a reference
cohort genuinely exists on that lane's catalogue.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Final

from openbiota.engines.metaphlan4 import RENAMED

#: Reverse of :data:`RENAMED`, for recognising a current name's older form.
#: Identity entries are dropped: the source map lists some names that did not
#: change, to record that they were checked, and rendering "formerly
#: Coprococcus comes" against *Coprococcus comes* would be nonsense.
_FORMERLY: Final[dict[str, str]] = {
    new: old for old, new in RENAMED.items() if old != new
}


def canonical(name: str) -> str:
    """The current accepted name for a species.

    Taxonomy moves: *Bacteroides vulgatus* is now *Phocaeicola vulgatus*,
    *Prevotella copri* is now *Segatella copri*. The older catalogue still uses
    the previous names, so without this the same organism appears twice in a
    merged list under two names and the species count is inflated.
    """
    return RENAMED.get(name, name)


def former_name(name: str) -> str | None:
    """The previous name for a species, if it was renamed, else ``None``."""
    return _FORMERLY.get(name)


@dataclass(frozen=True, slots=True)
class Lane:
    """One profiler's view of the sample.

    ``rankable`` records whether a reference cohort exists on this lane's
    catalogue. It gates percentiles and nothing else: a lane that cannot be
    ranked still contributes detections, abundances and strain joins.
    """

    #: Internal identifier for provenance. Never rendered to a reader.
    name: str
    #: Species name -> percent of classified reads.
    species: Mapping[str, float]
    #: True when a reference population exists on this lane's catalogue.
    rankable: bool
    #: Approximate catalogue size, for the technical record.
    catalogue_size: int = 0
    #: SGB rows, where the lane resolves genome bins.
    rows: Sequence[Mapping[str, Any]] = ()
    #: True for the one lane whose abundances form the composition. Every
    #: other lane contributes detections; its abundances are kept as
    #: secondary values and never added to the primary total.
    primary: bool = False
    #: The lane's own abundance unit. Only percent-like units may be shown
    #: as a secondary reading; anything else is detection-only.
    unit: str = "percent"
    #: The detection method family. Two lanes of one family (MetaPhlAn 3 and
    #: 4, say) are not two independent observations.
    method: str = "marker"
    #: True when the lane's numbers are not percent-like (SingleM coverage).
    detection_only: bool = False
    #: Reference release identifier, for the record.
    release: str = ""


@dataclass(frozen=True, slots=True)
class Organism:
    """One organism in the sample, with everything known about it.

    Whether it came from one profiler or several is deliberately not part of
    the reader-facing surface; :attr:`lanes` exists for the technical record
    and for tests, not for the page.
    """

    #: Current accepted species name, underscored as the catalogues write it.
    species: str
    #: Abundance in the primary composition, as a percentage. Zero when the
    #: primary lane did not see the organism; see `in_primary`.
    percent: float
    genus: str = ""
    #: True when the primary lane measured this organism, so `percent` is
    #: a real share of the composition. False means the organism was
    #: detected only by another lane: it is real, it is listed, and its
    #: share of the composition is unavailable rather than zero.
    in_primary: bool = True
    #: Abundance from a non-primary lane, when one saw the organism.
    #: Reported alongside, never summed into the composition.
    secondary_percent: float | None = None
    #: Genome bin identifier, where one profiler resolved it.
    sgb: str | None = None
    #: True for organisms known only from assembled genomes, with no
    #: cultured isolate and so no formal name.
    unnamed: bool = False
    #: Position against the reference population, where one exists.
    percentile: float | None = None
    #: Share of the reference population carrying this organism at all.
    prevalence: float | None = None
    #: Guild / functional group memberships.
    groups: tuple[str, ...] = ()
    #: Near the detection limit.
    trace: bool = False
    #: Carried by under a tenth of the reference population.
    rare: bool = False
    #: Strain fingerprint for this organism's dominant population.
    strain: Mapping[str, Any] | None = None
    #: Which profilers detected it. Provenance only.
    lanes: tuple[str, ...] = ()
    #: The previous name, when the organism was renamed.
    formerly: str | None = None
    #: The organism's name in the Genome Taxonomy Database, where it has
    #: one. For a bin with no formal name this *is* the name shown; for a
    #: named species it is kept as an alias so a reader holding a GTDB-based
    #: report can still find the organism.
    gtdb: str | None = None
    #: Genus supplied by GTDB when the detection catalogue gave none.
    gtdb_genus: str | None = None
    #: The MetaPhlAn database the SGB identifier belongs to (SGB numbers are
    #: not stable across releases).
    sgb_db: str | None = None
    #: Each lane's own identifier for this organism.
    native_ids: Mapping[str, str] = field(default_factory=dict)
    #: Fetchable reference genomes for this organism (accessions, GlobDB or
    #: panel genome ids), from whichever lanes named one. Provenance for the
    #: competitive confirmation; never shown as a name.
    genome_ids: tuple[str, ...] = ()
    #: Other keys this record was known by before lanes were merged - the
    #: MetaPhlAn 3 name of a species now keyed by its GTDB name, for one.
    #: Kept so the scoring cohort's row is still found, and for the explorer.
    aliases: tuple[str, ...] = ()
    #: supported | provisional | ambiguous. Spec 0.8.4 §6.
    status: str = "supported"
    #: Why the status is what it is, in words.
    confidence_basis: str = ""
    #: named_species | unnamed_species_cluster | unresolved_complex | higher_rank
    count_category: str = "named_species"
    #: baseline (seen by MP3/Jun23/GTDB232) | expansion (only the new lanes)
    incremental_gain: str = "baseline"
    #: Distinct detection methods that saw it.
    methods: tuple[str, ...] = ()
    #: The scoring lane's own reading (percent of its classified reads), when
    #: that lane saw the organism. A percentile from the scoring cohort was
    #: computed on this number, not on the primary lane's share.
    scoring_percent: float | None = None
    #: Which population the percentile was taken in: "scoring cohort" for
    #: the 3,027-adult reference set, "<lane> cohort, n=N" for a lane's own
    #: population. None when there is no percentile. Populations are never
    #: mixed; this names the one that was used.
    percentile_source: str | None = None
    #: The lane whose reading the percentile was taken on: "scoring" for the
    #: scoring cohort, else the inventory lane name of the population used.
    reference_lane: str | None = None
    #: That lane's own raw reading of the organism in this sample (percent of
    #: its classified reads), the number the share is compared with to say
    #: whether the rank describes the share on the page.
    reference_lane_percent: float | None = None
    #: The typical level among reference people who carry the organism (the
    #: median of their readings), in the units of the lane the percentile was
    #: taken in, and the sample's reading on that same lane. Their ratio is
    #: the deviation the report prints ("+400%", "-90%").
    reference_percent: float | None = None
    reference_reading: float | None = None
    #: Rank of `reference_reading` among the reference people who carry the
    #: organism, on their readings (midrank, 0-100). This is the percentile
    #: the report shows and judges levels by: it is taken on the same
    #: distribution as `reference_percent`, so above 50 always means above
    #: the typical carrier and below 50 below. `percentile` is the
    #: prevalence-aware rank (absence as a tied category), kept for the
    #: scoring machinery that was calibrated on it.
    carrier_percentile: float | None = None
    #: How many reference people carry the organism. A rank among four
    #: carriers is not a percentile; below `MIN_REFERENCE_CARRIERS` the level
    #: is not stated and the row says how few there were.
    reference_carriers: int | None = None
    #: Labels earlier catalogues gave this population that name a different
    #: species (Jun23 called SGB4571 "Ruminococcus gnavus"; GTDB and Jan26
    #: place it in Dorea hominis). Kept so an old report can be reconciled;
    #: never a name the organism answers to in a lookup.
    formerly_listed_as: tuple[str, ...] = ()
    #: How `percent` was arrived at (spec 0.8.4 §6, one denominator):
    #:   marker     the primary marker lane's own share
    #:   split      a genus total from the marker lane, divided among the
    #:              populations the whole-genome mapping tells apart
    #:   estimated  an organism the marker catalogue has no entry for; its
    #:              whole-genome reading rescaled to the marker lane's scale
    #:   member     counted inside a relative's share above; no share of its own
    #:   absorbed   no reading of its own in the marker lane; its share is the
    #:              reading the marker lane gave a relative whose call the
    #:              competitive confirmation rejected in this organism's favour
    #:   none       detected without a reading that can be placed
    share_basis: str = "marker"
    #: The marker lane's original reading, kept when `percent` was split.
    marker_percent: float | None = None
    #: For a `member`: the relative whose share it is counted within.
    counted_within: str | None = None
    #: For a `member`: True when the report judges this population on the
    #: member's own lane rank - because its reading is far from the relative's
    #: share (a distinct population the competition told apart), or because
    #: the relative it is counted within has no rank of its own, so this is
    #: the only rank the population has. False when it is the same population
    #: as the relative under another catalogue's name, judged once, as the
    #: relative.
    judged_on_own_rank: bool | None = None
    #: Share the marker lane had read under a relative that the competitive
    #: confirmation rejected; the reads belong here, so the share came here.
    absorbed_percent: float | None = None
    #: The rejected names whose shares were absorbed.
    absorbed_from: tuple[str, ...] = ()
    #: True when the scoring catalogue's reading and the composition share
    #: differ more than fivefold (both above trace): that catalogue measures
    #: this organism differently (its markers, its species boundary), so a
    #: percentile taken on its reading says nothing about the share shown.
    #: The percentile is kept in the record and not compared in the report.
    reference_conflict: bool = False

    @property
    def rankable(self) -> bool:
        """True when this organism can be placed against reference adults."""
        return self.percentile is not None

    @property
    def best_percent(self) -> float:
        """The abundance to sort and speak by: primary, else secondary."""
        if self.in_primary:
            return self.percent
        return self.secondary_percent or 0.0

    @property
    def ranked_percent(self) -> float:
        """The reading a quoted percentile was computed on.

        A scoring-cohort percentile rests on the scoring lane's reading; a
        trace share in the primary lane beside a 98th percentile is not a
        contradiction, it is two lanes' numbers, and the headline must quote
        the one the rank belongs to."""
        if self.percentile is not None and (self.percentile_source or "scoring cohort") == "scoring cohort" \
                and self.scoring_percent is not None:
            return self.scoring_percent
        return self.best_percent

    @property
    def few_reference_carriers(self) -> bool:
        """True when too few reference people carry the organism for a rank among them to mean anything."""
        return self.reference_carriers is not None and self.reference_carriers < MIN_REFERENCE_CARRIERS

    @property
    def deviation_percent(self) -> float | None:
        """Percent above (+) or below (-) the typical carrier's level, on the lane the reference was measured."""
        if self.reference_conflict or self.few_reference_carriers:
            return None
        if not self.reference_percent or not self.reference_reading or self.reference_reading <= 0:
            return None
        return (self.reference_reading / self.reference_percent - 1.0) * 100.0

    @property
    def level_percentile(self) -> float | None:
        """The percentile that describes the level: rank among carriers, on the
        same distribution the deviation is measured against. Falls back to the
        prevalence-aware percentile only for records written before the carrier
        rank existed. None when the reference catalogue reads the organism
        differently from the composition (`reference_conflict`) or when too few
        reference people carry it for a rank among them to mean anything."""
        if self.reference_conflict or self.few_reference_carriers:
            return None
        return self.carrier_percentile if self.carrier_percentile is not None else self.percentile

    @property
    def comparable_percentile(self) -> float | None:
        """Alias of `level_percentile`, kept for callers that used the earlier name."""
        return self.level_percentile

    @property
    def resolved_to_strain(self) -> bool:
        return self.strain is not None

    @property
    def display(self) -> str:
        """The name to print.

        A named species prints as itself. An unnamed bin prints its GTDB
        placement when one exists — ``Gemmiger sp937890665`` rather than
        ``GGB45596_SGB63306`` — because the former carries a genus, sorts
        with its relatives and matches what other reports print for the
        same organism.
        """
        if self.unnamed and self.gtdb:
            head, _, tail = self.gtdb.partition(" ")
            # ``MOTU40_058831 MOTU40_058831``: a cluster whose genus is its own id
            return self.gtdb if (tail and tail != head) else head
        # A GTDB-derived key keeps GTDB's spelling (``Blautia_A wexlerae``):
        # the suffix marks a genus GTDB has split, and "Blautia A" is not a name.
        if self.gtdb and self.gtdb.replace(" ", "_") == self.species:
            return self.gtdb
        return self.species.replace("_", " ")

    @property
    def provisional_name(self) -> bool:
        """True when the printed name is a placeholder, not a formal binomial."""
        return self.unnamed

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "species": self.species,
            "percent": round(self.percent, 4),
            "in_primary": self.in_primary,
            "genus": self.genus,
            "detected_by": list(self.lanes),
        }
        if self.secondary_percent is not None:
            out["secondary_percent"] = round(self.secondary_percent, 4)
        if self.sgb:
            out["sgb"] = self.sgb
        if self.percentile is not None:
            out["percentile"] = round(self.percentile, 1)
        if self.prevalence is not None:
            out["prevalence"] = round(self.prevalence, 4)
        if self.groups:
            out["groups"] = list(self.groups)
        if self.formerly:
            out["formerly"] = self.formerly
        if self.gtdb:
            out["gtdb"] = self.gtdb
        if self.gtdb_genus:
            out["gtdb_genus"] = self.gtdb_genus
        if self.unnamed:
            out["unnamed"] = True
        if self.trace:
            out["trace"] = True
        if self.rare:
            out["rare"] = True
        if self.strain is not None:
            out["strain"] = dict(self.strain)
        if self.sgb_db:
            out["sgb_db"] = self.sgb_db
        if self.native_ids:
            out["native_ids"] = dict(self.native_ids)
        if self.genome_ids:
            out["genome_ids"] = list(self.genome_ids)
        if self.aliases:
            out["aliases"] = list(self.aliases)
        out["status"] = self.status
        if self.confidence_basis:
            out["confidence_basis"] = self.confidence_basis
        out["count_category"] = self.count_category
        out["incremental_gain"] = self.incremental_gain
        if self.methods:
            out["methods"] = list(self.methods)
        if self.percentile_source:
            out["percentile_source"] = self.percentile_source
        if self.reference_lane:
            out["reference_lane"] = self.reference_lane
        if self.reference_lane_percent is not None:
            out["reference_lane_percent"] = round(self.reference_lane_percent, 4)
        if self.reference_percent is not None:
            out["reference_percent"] = round(self.reference_percent, 4)
        if self.reference_reading is not None:
            out["reference_reading"] = round(self.reference_reading, 4)
        if self.carrier_percentile is not None:
            out["carrier_percentile"] = round(self.carrier_percentile, 1)
        if self.reference_carriers is not None:
            out["reference_carriers"] = int(self.reference_carriers)
        if self.formerly_listed_as:
            out["formerly_listed_as"] = list(self.formerly_listed_as)
        if self.deviation_percent is not None:
            out["deviation_percent"] = round(self.deviation_percent, 1)
        out["share_basis"] = self.share_basis
        if self.marker_percent is not None:
            out["marker_percent"] = round(self.marker_percent, 4)
        if self.counted_within:
            out["counted_within"] = self.counted_within
        if self.judged_on_own_rank is not None:
            out["judged_on_own_rank"] = bool(self.judged_on_own_rank)
        if self.absorbed_percent:
            out["absorbed_percent"] = round(self.absorbed_percent, 4)
            out["absorbed_from"] = list(self.absorbed_from)
        if self.reference_conflict:
            out["reference_conflict"] = True
        if self.scoring_percent is not None:
            out["scoring_percent"] = round(self.scoring_percent, 4)
        return out


@dataclass(slots=True)
class Inventory:
    """Every organism detected in one sample, from every profiler that ran."""

    organisms: list[Organism] = field(default_factory=list)
    #: Lane names that contributed, for the technical record.
    lanes: tuple[str, ...] = ()
    #: The lane whose abundances form the composition.
    primary_lane: str = ""
    #: The primary lane's own estimate of sequence it could not place, as a
    #: percentage. Shown as an explicit remainder so the composition adds
    #: up; never silently dropped and never relabelled as organisms.
    unclassified_percent: float | None = None
    #: How shares were unified onto the primary scale (see `unify_shares`).
    share_unification: Mapping[str, Any] | None = None

    @property
    def estimated_percent_total(self) -> float:
        """Shares estimated for organisms the marker catalogue cannot place; drawn from the unclassified band."""
        return sum(o.percent for o in self.organisms if o.share_basis == "estimated")

    @property
    def unplaced_percent(self) -> float | None:
        """The unclassified remainder after the estimated organisms are drawn from it."""
        if self.unclassified_percent is None:
            return None
        return max(0.0, float(self.unclassified_percent) - self.estimated_percent_total)

    def __iter__(self) -> Any:
        return iter(self.organisms)

    def __len__(self) -> int:
        return len(self.organisms)

    @property
    def in_composition(self) -> list[Organism]:
        """Organisms with a share of the primary composition."""
        return [o for o in self.organisms if o.in_primary]

    @property
    def secondary_only(self) -> list[Organism]:
        """Organisms detected only by a lane other than the primary."""
        return [o for o in self.organisms if not o.in_primary]

    @property
    def composition_total(self) -> float:
        """Sum of primary shares. Must not exceed 100."""
        return sum(o.percent for o in self.in_composition)

    @property
    def named(self) -> list[Organism]:
        """Organisms with a formal name, largest first."""
        return [o for o in self.organisms if not o.unnamed]

    @property
    def unnamed(self) -> list[Organism]:
        """Organisms known only from assembled genomes."""
        return [o for o in self.organisms if o.unnamed]

    @property
    def percentile_sources(self) -> dict[str, int]:
        """Population name -> organisms ranked in it. Never merged; listed."""
        out: dict[str, int] = {}
        for o in self.organisms:
            if o.percentile is not None:
                key = o.percentile_source or "scoring cohort"
                out[key] = out.get(key, 0) + 1
        return out

    @property
    def rankable(self) -> list[Organism]:
        """Organisms that can be placed against the reference population."""
        return [o for o in self.organisms if o.rankable]

    @property
    def strain_resolved(self) -> list[Organism]:
        return [o for o in self.organisms if o.resolved_to_strain]

    @property
    def genera(self) -> set[str]:
        return {o.genus for o in self.organisms if o.genus}

    def get(self, species: str) -> Organism | None:
        """Look one organism up under its current or former name."""
        key = canonical(species)
        wanted = species.replace("_", " ")
        for o in self.organisms:
            if o.species == key or o.formerly == species or (o.gtdb or "") == wanted:
                return o
        # the record the reference catalogue itself reported under this name
        for o in self.organisms:
            if (o.native_ids or {}).get("scoring") in (species, key):
                return o
        # a name the organism was listed under before lanes were merged (the
        # same species under another catalogue's spelling; a label that names
        # a different species is in `formerly_listed_as` and does not answer)
        for o in self.organisms:
            if key in o.aliases or species in o.aliases:
                return o
        return None

    def present(self, species: str) -> bool:
        """True when any profiler detected this organism."""
        return self.get(species) is not None

    def counts(self) -> dict[str, int]:
        """Headline counts, for prose that needs to state them.

        Spec 0.8.4 §6: named species, unnamed species-level clusters,
        unresolved complexes and higher-rank findings are counted
        separately, strains separately again, and provisional calls are
        never folded into the supported count.
        """
        by_cat = dict.fromkeys(("named_species", "unnamed_species_cluster", "unresolved_complex", "higher_rank"), 0)
        for o in self.organisms:
            by_cat[o.count_category] = by_cat.get(o.count_category, 0) + 1
        return {
            "organisms": len(self.organisms),
            "named": len(self.named),
            "unnamed": len(self.unnamed),
            "genera": len(self.genera),
            "rankable": len(self.rankable),
            "strain_resolved": len(self.strain_resolved),
            "supported": sum(1 for o in self.organisms if o.status == "supported"),
            "provisional": sum(1 for o in self.organisms if o.status == "provisional"),
            "ambiguous": sum(1 for o in self.organisms if o.status == "ambiguous"),
            "expansion": sum(1 for o in self.organisms if o.incremental_gain == "expansion"),
            "expansion_supported": sum(1 for o in self.organisms
                                       if o.incremental_gain == "expansion" and o.status == "supported"),
            **{f"category_{k}": v for k, v in by_cat.items()},
        }

    #: Rejected confirmation calls, kept for the report (set by from_json).
    rejected_records: list[dict[str, Any]] = field(default_factory=list)

    @property
    def newly_supported(self) -> list[Organism]:
        """Organisms the installed baseline never saw, now supported."""
        return [o for o in self.organisms if o.incremental_gain == "expansion" and o.status == "supported"]

    @property
    def provisional(self) -> list[Organism]:
        return [o for o in self.organisms if o.status == "provisional"]

    def to_json(self) -> dict[str, Any]:
        return {
            "n_organisms": len(self.organisms),
            "n_named": len(self.named),
            "n_unnamed": len(self.unnamed),
            "n_genera": len(self.genera),
            "n_rankable": len(self.rankable),
            "n_strain_resolved": len(self.strain_resolved),
            "n_in_composition": len(self.in_composition),
            "n_secondary_only": len(self.secondary_only),
            "composition_total_percent": round(self.composition_total, 3),
            "unclassified_percent": self.unclassified_percent,
            "unplaced_percent": None if self.unplaced_percent is None else round(self.unplaced_percent, 3),
            "share_unification": dict(self.share_unification) if self.share_unification else None,
            "counts": self.counts(),
            # Spec 0.8.4 §7: the per-sample comparison a reader can act on.
            "comparison": {
                "newly_supported": [
                    {"organism": o.display, "id": o.species, "native_ids": dict(o.native_ids),
                     "reading": o.best_percent, "detected_by": list(o.lanes), "basis": o.confidence_basis}
                    for o in self.newly_supported
                ],
                "provisional": [
                    {"organism": o.display, "id": o.species, "detected_by": list(o.lanes), "reason": o.confidence_basis}
                    for o in self.provisional
                ],
                "unresolved_complexes": [
                    {"organism": o.display, "id": o.species, "detected_by": list(o.lanes)}
                    for o in self.organisms if o.count_category == "unresolved_complex"
                ],
                "renamed": [
                    {"organism": o.display, "formerly": o.formerly} for o in self.organisms if o.formerly
                ],
                "meaning": (
                    "newly_supported: organisms the installed baseline (MetaPhlAn 3, MetaPhlAn 4 Jun23, GTDB R232 "
                    "sylph) never saw and the expanded lanes support; provisional: one method at a marginal reading, "
                    "listed but not counted as confirmed; unresolved_complexes: bins that map to several GTDB R232 "
                    "species; renamed: current name with the previous one."
                ),
            },
            "primary_lane": self.primary_lane,
            "contributing_lanes": list(self.lanes),
            "organisms": [o.to_json() for o in self.organisms],
            "meaning": (
                "every organism detected by any profiler that ran, merged under current "
                "taxonomic names. Percentiles are present only where a reference population "
                "exists on the catalogue that measured the organism"
            ),
        }


def from_json(blob: Mapping[str, Any] | None) -> Inventory | None:
    """Rebuild an inventory from the serialised form in ``results.json``.

    The renderer runs from the results file rather than from live objects, so
    the inventory has to survive the round trip.
    """
    if not blob or blob.get("status") == "stage_error":
        return None
    rows = blob.get("organisms") or []
    if not rows:
        return None
    from openbiota import gtdb as _gtdb

    def _placement(r: Mapping[str, Any]) -> tuple[str | None, str | None]:
        # A results file written before GTDB names existed still carries
        # the bin identifier, which is all the lookup needs.
        if r.get("gtdb") or not r.get("sgb"):
            return r.get("gtdb"), r.get("gtdb_genus")
        hit = _gtdb.lookup(r.get("sgb"), str(r.get("sgb_db") or _gtdb.DEFAULT_DB))
        return (hit.display, hit.genus) if hit is not None else (None, None)

    organisms = [
        Organism(
            species=str(r.get("species") or ""),
            percent=float(r.get("percent") or 0.0),
            in_primary=bool(r.get("in_primary", True)),
            secondary_percent=r.get("secondary_percent"),
            genus=str(r.get("genus") or ""),
            sgb=r.get("sgb"),
            unnamed=bool(r.get("unnamed")),
            percentile=r.get("percentile"),
            prevalence=r.get("prevalence"),
            groups=tuple(r.get("groups") or ()),
            trace=bool(r.get("trace")),
            rare=bool(r.get("rare")),
            strain=r.get("strain"),
            lanes=tuple(r.get("detected_by") or ()),
            formerly=r.get("formerly"),
            gtdb=_placement(r)[0],
            gtdb_genus=_placement(r)[1],
            sgb_db=r.get("sgb_db"),
            native_ids=dict(r.get("native_ids") or {}),
            genome_ids=tuple(r.get("genome_ids") or ()),
            aliases=tuple(r.get("aliases") or ()),
            status=str(r.get("status") or "supported"),
            confidence_basis=str(r.get("confidence_basis") or ""),
            count_category=str(r.get("count_category") or ("unnamed_species_cluster" if r.get("unnamed") else "named_species")),
            incremental_gain=str(r.get("incremental_gain") or "baseline"),
            methods=tuple(r.get("methods") or ()),
            percentile_source=r.get("percentile_source"),
            reference_lane=r.get("reference_lane"),
            reference_lane_percent=r.get("reference_lane_percent"),
            scoring_percent=r.get("scoring_percent"),
            reference_percent=r.get("reference_percent"),
            reference_reading=r.get("reference_reading"),
            carrier_percentile=r.get("carrier_percentile"),
            reference_carriers=r.get("reference_carriers"),
            formerly_listed_as=tuple(r.get("formerly_listed_as") or ()),
            share_basis=str(r.get("share_basis") or ("marker" if r.get("in_primary", True) else "none")),
            marker_percent=r.get("marker_percent"),
            counted_within=r.get("counted_within"),
            judged_on_own_rank=r.get("judged_on_own_rank"),
            absorbed_percent=r.get("absorbed_percent"),
            absorbed_from=tuple(r.get("absorbed_from") or ()),
            reference_conflict=bool(r.get("reference_conflict")),
        )
        for r in rows
    ]
    return Inventory(
        organisms=organisms,
        lanes=tuple(blob.get("contributing_lanes") or ()),
        primary_lane=str(blob.get("primary_lane") or ""),
        rejected_records=list(blob.get("rejected") or []),
        unclassified_percent=blob.get("unclassified_percent"),
        share_unification=blob.get("share_unification"),
    )


def build(
    lanes: Sequence[Lane],
    *,
    ranked_rows: Iterable[Any] = (),
    strain: Mapping[str, Any] | None = None,
    unclassified_percent: float | None = None,
    lane_cohorts: Mapping[str, Any] | None = None,
) -> Inventory:
    """Merge every lane into one inventory.

    ``ranked_rows`` are the reference-scored rows (``community.SpeciesRow``)
    that carry percentile, prevalence and group membership. They are matched
    onto the merged organisms by canonical name. A rank is attached only when
    the organism was measured on a lane that a reference population exists
    for, which is what keeps an unrankable detection from acquiring a
    meaningless percentile.
    """
    merged: dict[str, dict[str, Any]] = {}

    for lane in lanes:
        # SGB rows carry genus, kingdom and bin identity; the flat species map
        # is the fallback for lanes that do not resolve bins.
        seen_in_rows: set[str] = set()
        for row in lane.rows:
            name = str(row.get("species") or "")
            value = row.get("percent")
            if not name or (value is not None and float(value) <= 0) or (value is None and not lane.detection_only):
                continue
            key = canonical(name)
            seen_in_rows.add(key)
            native = {lane.name: str(row.get("native_id") or row.get("sgb") or row.get("accession") or "")}
            _accumulate(
                merged, key, lane,
                percent=float(value) if (value is not None and not lane.detection_only) else 0.0,
                genus=str(row.get("genus") or ""),
                sgb=str(row.get("sgb") or "") or None,
                unnamed=bool(row.get("unnamed")),
                gtdb=str(row.get("gtdb") or "") or None,
                sgb_db=str(row.get("sgb_db") or "") or None,
                native_ids={k: v for k, v in native.items() if v},
                lane_status=str(row.get("status") or "") or None,
                complex_=bool(row.get("complex")),
                genome_id=str(row.get("genome_id") or "") or None,
                aliases=tuple(str(a) for a in (row.get("aliases") or ()) if a),
            )
        for name, percent in lane.species.items():
            key = canonical(name)
            if percent <= 0 or key in seen_in_rows:
                continue
            _accumulate(merged, key, lane, percent=float(percent))

    _reconcile_by_alias(merged)

    by_name = {canonical(getattr(r, "species", "")): r for r in ranked_rows}
    rankable_lanes = {ln.name for ln in lanes if ln.rankable}
    # Lanes with their own reference population (spec: rankable by cohort).
    # Loaded from disk unless the caller supplied them; an empty mapping
    # means "none", so tests can switch the behaviour off.
    if lane_cohorts is None:
        from openbiota.expansion import cohorts as _cohorts

        lane_cohorts = _cohorts.load_all()
    lane_ranks: dict[str, dict[str, Any]] = {}
    for ln in lanes:
        lc = lane_cohorts.get(ln.name)
        if lc is None or not ln.species:
            continue
        from openbiota.expansion import cohorts as _cohorts

        lane_ranks[ln.name] = _cohorts.rank(lc, ln.species)
    primary_lanes = [ln.name for ln in lanes if ln.primary]
    if len(primary_lanes) > 1:
        raise ValueError(f"exactly one lane may be primary, got {primary_lanes}")
    # With no lane declared primary the broadest one (most organisms) is.
    primary = primary_lanes[0] if primary_lanes else max(
        lanes, key=lambda ln: len(ln.rows) or len(ln.species), default=None)
    primary = primary if isinstance(primary, str) else (primary.name if primary else "")
    by_sgb: dict[str, Mapping[str, Any]] = {
        str(row.get("sgb")): row
        for row in ((strain or {}).get("organisms") or [])
        if row.get("sgb")
    }

    from openbiota import gtdb as _gtdb

    lane_by_name = {ln.name: ln for ln in lanes}
    organisms: list[Organism] = []
    for key, rec in merged.items():
        # The scoring cohort's row is keyed by the scoring lane's own name,
        # which a merge may have replaced by the newer catalogue's: look
        # under every name the record was known by.
        candidates_ranked = [by_name[n] for n in (key, *sorted(rec.get("aliases") or ())) if n in by_name]
        ranked = next((r for r in candidates_ranked if getattr(r, "percentile", None) is not None),
                      candidates_ranked[0] if candidates_ranked else None)
        placement = _gtdb.lookup(rec["sgb"], rec.get("sgb_db") or _gtdb.DEFAULT_DB) if rec["sgb"] else None
        gtdb_name = placement.display if placement is not None else rec.get("gtdb")
        if rec.get("sibling_species"):
            # two catalogues placed one population in sibling species: the
            # complex label names both, and neither is asserted alone
            gtdb_name = rec.get("gtdb")
        gtdb_genus = placement.genus if placement is not None else (
            (rec.get("gtdb") or "").split(" ")[0] if rec.get("gtdb") else None)
        # A percentile is only meaningful if a reference population exists on
        # a catalogue that actually measured this organism.
        placeable = bool(rec["lanes"] & rankable_lanes)
        by_lane: dict[str, float] = rec["by_lane"]
        in_primary = primary in by_lane
        percent = float(by_lane.get(primary, 0.0))
        # Secondary reading: the largest percent-like value from another
        # lane. Detection-only lanes (coverage units) never supply it.
        others = [v for k, v in by_lane.items() if k != primary and lane_by_name[k].unit == "percent"
                  and not lane_by_name[k].detection_only]
        secondary = max(others) if others else None
        status, basis, methods = _status_of(rec, lane_by_name, primary, percent)
        gain = "baseline" if (rec["lanes"] & BASELINE_LANES) else "expansion"
        organisms.append(Organism(
            species=key,
            percent=percent,
            in_primary=in_primary,
            secondary_percent=secondary,
            genus=rec["genus"] or gtdb_genus or _genus_of(key),
            sgb=rec["sgb"],
            unnamed=rec["unnamed"],
            percentile=getattr(ranked, "percentile", None) if placeable else None,
            prevalence=getattr(ranked, "cohort_prevalence", None) if placeable else None,
            groups=tuple(getattr(ranked, "groups", ()) or ()),
            trace=bool(getattr(ranked, "trace", False)),
            rare=bool(getattr(ranked, "rare_in_cohort", False)),
            strain=by_sgb.get(str(rec["sgb"])) if rec["sgb"] else None,
            lanes=tuple(sorted(rec["lanes"])),
            formerly=former_name(key),
            gtdb=gtdb_name,
            gtdb_genus=gtdb_genus or None,
            sgb_db=rec.get("sgb_db"),
            native_ids=dict(rec.get("native_ids") or {}),
            genome_ids=tuple(sorted(rec.get("genome_ids") or ())),
            aliases=tuple(sorted(a for a in (rec.get("aliases") or ()) if a != key)),
            status=status,
            confidence_basis=basis,
            count_category=_count_category(rec, placement, key),
            incremental_gain=gain,
            methods=methods,
            scoring_percent=(float(by_lane["scoring"]) if "scoring" in by_lane else None),
            # the reading the rank was taken on - the scoring lane's value
            # rescaled to its classified fraction - so the deviation and the
            # carrier percentile describe the same number
            reference_percent=(getattr(ranked, "reference_median", None) if placeable else None),
            reference_reading=((float(getattr(ranked, "percent", 0.0) or 0.0) or None)
                               if (ranked is not None and placeable) else None),
            carrier_percentile=(getattr(ranked, "carrier_percentile", None) if placeable else None),
            reference_carriers=(getattr(ranked, "reference_carriers", None) if placeable else None),
            formerly_listed_as=tuple(sorted(rec.get("misnamed") or ())),
            share_basis="marker" if in_primary else "none",
            reference_conflict=_reads_differently(percent if in_primary else None, by_lane.get("scoring")),
        ))

    # Organisms the scoring population never measured are ranked in their
    # own lane's population, when one exists. The primary lane's cohort is
    # preferred, then the others in registration order; exactly one is used
    # and it is named on the organism.
    if lane_ranks:
        order = [primary] + [ln.name for ln in lanes if ln.name != primary]
        for i, o in enumerate(organisms):
            if o.percentile is not None and not o.reference_conflict:
                organisms[i] = replace(o, percentile_source="scoring cohort", reference_lane="scoring")
                continue
            # An organism the scoring catalogue reads fivefold differently
            # from the composition (its species boundary, its markers) is
            # ranked in another lane's population instead - the composition
            # lane's own where one exists, else a whole-genome lane that saw
            # it. `unify_shares` then compares that lane's reading with the
            # share and withholds the rank where they too disagree.
            for lane_name in order:
                ranks = lane_ranks.get(lane_name)
                if not ranks or lane_name not in o.lanes:
                    continue
                lane_obj = lane_by_name[lane_name]
                # the lane's own name for the organism: match by canonical key
                hit = None
                known = {o.species, *o.aliases}
                for sp in lane_obj.species:
                    if canonical(sp) in known:
                        hit = ranks.get(sp)
                        if hit is not None and hit.percentile is not None:
                            break
                if hit is None or hit.percentile is None:
                    continue
                organisms[i] = replace(
                    o, percentile=hit.percentile, prevalence=hit.prevalence,
                    rare=bool(hit.rare), percentile_source=lane_cohorts[lane_name].source_label,
                    reference_lane=lane_name,
                    reference_lane_percent=(float(merged[o.species]["by_lane"].get(lane_name) or 0.0) or None),
                    reference_percent=getattr(hit, "reference_median", None),
                    reference_reading=(getattr(hit, "reading", None)
                                       or float(merged[o.species]["by_lane"].get(lane_name) or 0.0) or None),
                    carrier_percentile=getattr(hit, "carrier_percentile", None),
                    reference_carriers=getattr(hit, "reference_carriers", None),
                    reference_conflict=False)
                break
            else:
                if o.percentile is not None:
                    organisms[i] = replace(o, percentile_source="scoring cohort", reference_lane="scoring")
    else:
        organisms = [replace(o, percentile_source="scoring cohort", reference_lane="scoring")
                     if o.percentile is not None else o for o in organisms]

    # Composition members first, by share; then organisms seen only by a
    # secondary lane, by their secondary reading.
    organisms.sort(key=lambda o: (not o.in_primary, -o.best_percent))
    return Inventory(organisms=organisms, lanes=tuple(ln.name for ln in lanes),
                     primary_lane=primary, unclassified_percent=unclassified_percent)


# --------------------------------------------------------------------------- #
# one share of the whole per organism
# --------------------------------------------------------------------------- #


def display_of(rec: Mapping[str, Any]) -> str:
    """The printed name of a serialised record - the same rule as `Organism.display`,
    so a confirmation verdict written against the display finds its record."""
    gtdb, species = str(rec.get("gtdb") or ""), str(rec.get("species") or "")
    if rec.get("unnamed") and gtdb:
        head, _, tail = gtdb.partition(" ")
        return gtdb if (tail and tail != head) else head
    if gtdb and gtdb.replace(" ", "_") == species:
        return gtdb
    return species.replace("_", " ")


_display_of = display_of


def level_percentile_of(rec: Mapping[str, Any]) -> float | None:
    """`Organism.level_percentile` for a serialised record: the rank among reference
    carriers when it can be stated, else the prevalence-aware rank of an older file,
    and None where the reference cannot place the organism (conflict, too few carriers)."""
    if rec.get("reference_conflict"):
        return None
    carriers = rec.get("reference_carriers")
    if carriers is not None and int(carriers) < MIN_REFERENCE_CARRIERS:
        return None
    cp = rec.get("carrier_percentile")
    return float(cp) if cp is not None else (float(rec["percentile"]) if rec.get("percentile") is not None else None)


def _genus_key_of(rec: Mapping[str, Any]) -> str:
    g = str(rec.get("gtdb_genus") or rec.get("genus") or "")
    return g.split("_")[0].lower() if g else ""


def unify_shares(blob: dict[str, Any], confirmation: Mapping[str, Any] | None) -> dict[str, Any]:
    """Give every organism one share of the whole, on the primary lane's scale.

    The primary marker lane measures the composition, but it names
    populations by its own catalogue: a species it has one entry for may be
    two populations that the whole-genome catalogues tell apart
    (*Lachnospira eligens* and *L. eligens_A*), and a cluster it has no
    entry for at all it cannot place. Left as they were, those organisms
    appeared with a reading from another lane - a different denominator -
    beside the composition shares, and a 7% reading sat in the middle of a
    list sorted by share. This step resolves them into the one scale:

    * **split** - a marker species is divided among the populations the
      marker catalogue cannot tell apart from it - its GTDB siblings (*X*,
      *X_A*, *X_B*) that the competitive mapping distinguished - in
      proportion to the reads that mapped uniquely to each genome per
      megabase. The species total is conserved and the marker lane's own
      reading is kept as `marker_percent`; readings for species the marker
      lane does distinguish are never redistributed between them.
    * **estimated** - an organism whose genus the marker lane did not place
      at all takes its whole-genome reading rescaled to the marker lane's
      scale (the median ratio of the two lanes over organisms both measured);
      those shares are drawn from the marker lane's unclassified band.
    * **member** - an organism the whole-genome lanes saw inside a genus the
      marker lane did measure, but which the mapping could not apportion
      (no genome in the competition, or evidence from one method only): it
      is counted within its relatives' shares above and carries none of its
      own.

    Nothing is summed across lanes; the primary lane's total is the total.
    """
    records: list[dict[str, Any]] = list(blob.get("organisms") or [])
    if not records:
        return blob
    # start from the marker lane's own readings, whatever an earlier pass did
    for rec in records:
        basis = rec.get("share_basis") or ("marker" if rec.get("in_primary", True) else "none")
        if rec.get("marker_percent") is not None:
            rec["percent"], rec["in_primary"] = float(rec["marker_percent"]), True
        elif basis in ("split", "estimated", "member", "absorbed", "none"):
            rec["percent"], rec["in_primary"] = 0.0, False
        rec.pop("counted_within", None)
        rec.pop("judged_on_own_rank", None)
        rec.pop("absorbed_percent", None)
        rec.pop("absorbed_from", None)
        rec["share_basis"] = "marker" if rec.get("in_primary") else "none"
        rec.pop("marker_percent", None)

    # unique fragments per genome from the competition, by organism display
    mapping: dict[str, tuple[int, int]] = {}
    for g in ((confirmation or {}).get("genomes") or {}).values():
        u, ln = int(g.get("unique_fragments") or 0), int(g.get("genome_length") or 0)
        if g.get("organism") and ln > 0:
            cur = mapping.get(g["organism"])
            mapping[g["organism"]] = (cur[0] + u, cur[1] + ln) if cur else (u, ln)
    if not mapping:
        for v in (confirmation or {}).get("verdicts") or []:
            if v.get("status") in ("supported", "provisional") and int(v.get("genome_length") or 0) > 0:
                mapping[v["organism"]] = (int(v.get("unique_fragments") or 0), int(v["genome_length"]))

    # the secondary lanes' scale relative to the marker lane
    ratios = sorted(float(r["percent"]) / float(r["secondary_percent"]) for r in records
                    if r.get("in_primary") and float(r.get("percent") or 0) > 0.05
                    and float(r.get("secondary_percent") or 0) > 0.05)
    scale = min(4.0, max(0.25, ratios[len(ratios) // 2])) if len(ratios) >= 10 else 1.0

    groups: dict[str, list[dict[str, Any]]] = {}
    for i, rec in enumerate(records):
        groups.setdefault(_genus_key_of(rec) or f"#{i}", []).append(rec)

    estimated_total = 0.0
    n_split = n_estimated = n_member = 0
    for members in groups.values():
        primaries = [r for r in members if r.get("in_primary")]
        secondary = [r for r in members if not r.get("in_primary")]
        if not secondary:
            continue
        if primaries:
            anchor = max(primaries, key=lambda r: float(r.get("percent") or 0))
            strong = [r for r in secondary if r.get("status") == "supported"]
            # A marker species is divided only among the populations the
            # marker catalogue cannot tell apart from it: its GTDB siblings
            # (X, X_A, X_B ...). The marker lane's readings for species it
            # does distinguish stand; a genus total is never redistributed
            # between them by mapping density.
            def base_of(r: dict[str, Any]) -> str:
                g = str(r.get("gtdb") or "")
                b = _sibling_base(g) if g and " / " not in g and not g.endswith(")") else None
                return b or str(r.get("species") or "")

            by_base: dict[str, list[dict[str, Any]]] = {}
            for r in primaries:
                by_base.setdefault(base_of(r), []).append(r)
            for r in strong:
                b = base_of(r)
                if b in by_base:
                    by_base[b].append(r)
            for members_of_base in by_base.values():
                pool_members = [r for r in members_of_base if _display_of(r) in mapping]
                pool_has_new = any(not r.get("in_primary") for r in pool_members)
                pool_has_marker = any(r.get("in_primary") for r in pool_members)
                if not (pool_has_new and pool_has_marker):
                    continue
                pool = sum(float(r.get("percent") or 0) for r in pool_members if r.get("in_primary"))
                weights = {id(r): mapping[_display_of(r)][0] / (mapping[_display_of(r)][1] / 1e6) for r in pool_members}
                total_w = sum(weights.values())
                if total_w > 0 and pool > 0:
                    for r in pool_members:
                        if r.get("in_primary"):
                            r["marker_percent"] = float(r.get("percent") or 0)
                        r["percent"] = pool * weights[id(r)] / total_w
                        r["in_primary"] = True
                        r["share_basis"] = "split"
                        n_split += 1
            for r in secondary:
                if r.get("share_basis") == "split":
                    continue
                r["share_basis"] = "member"
                r["percent"], r["in_primary"] = 0.0, False
                r["counted_within"] = _display_of(anchor)
                n_member += 1
        else:
            for r in secondary:
                reading = float(r.get("secondary_percent") or 0)
                if reading > 0:
                    r["percent"] = reading * scale
                    r["in_primary"] = True
                    r["share_basis"] = "estimated"
                    estimated_total += r["percent"]
                    n_estimated += 1
                else:
                    r["share_basis"] = "none"
    # A rejected call's reads belong to the relative that took them in the
    # competition; when the marker lane had given that call a share, the
    # share goes with the reads, so the whole still adds up and the relative
    # carries what the mapping says is its. Kept on both records. The
    # competition names a relative by the name its genome was entered under:
    # the record's display name, or its GTDB name for a genome drawn from
    # the species' cluster.
    by_display: dict[str, dict[str, Any]] = {_display_of(r): r for r in records}
    for r in records:
        for name in (str(r.get("gtdb") or ""), str(r.get("species") or "").replace("_", " "),
                     *(str(a).replace("_", " ") for a in (r.get("aliases") or ()))):
            if name:
                by_display.setdefault(name, r)
    n_absorbed = 0
    absorbed_total = 0.0
    for rej in blob.get("rejected") or []:
        share = float(rej.get("marker_percent") if rej.get("marker_percent") is not None else (rej.get("percent") or 0.0))
        target = by_display.get(str(rej.get("reads_belong_to") or ""))
        if not (rej.get("in_primary") and share > 0 and target is not None):
            continue
        if target.get("in_primary") and target.get("share_basis") != "estimated":
            # the marker lane read this organism under two names: both readings are its
            if target.get("marker_percent") is None:
                target["marker_percent"] = float(target.get("percent") or 0.0)
        else:
            # no marker reading of its own until now; a whole-genome estimate
            # of the same population gives way to the marker lane's reading
            # of it, and is no longer drawn from the unclassified band
            if target.get("share_basis") == "estimated":
                estimated_total -= float(target.get("percent") or 0.0)
                n_estimated -= 1
            target["in_primary"], target["percent"] = True, 0.0
            target["share_basis"] = "absorbed"
            target.pop("counted_within", None)
        target["percent"] = float(target.get("percent") or 0.0) + share
        target["absorbed_percent"] = round(float(target.get("absorbed_percent") or 0.0) + share, 4)
        target.setdefault("absorbed_from", []).append(_display_of(rej))
        n_absorbed += 1
        absorbed_total += share
    # Only now, with every share final (absorption included), can a member be
    # told apart from the relative it is counted within: the same population
    # under another catalogue's name reads within a factor of two of that
    # relative's share; a distinct population does not.
    final_by_display = {_display_of(r): r for r in records}
    for rec in records:
        if rec.get("share_basis") != "member":
            continue
        anchor = final_by_display.get(str(rec.get("counted_within") or ""))
        distinct = _distinct_from(
            float(rec.get("secondary_percent") or 0.0),
            float(anchor.get("percent") or 0.0) if anchor is not None else 0.0)
        # the same population as its relative is judged once, as the
        # relative - unless the relative has no rank, in which case the
        # member's lane rank is the only reading of that population there is
        rec["judged_on_own_rank"] = bool(distinct or anchor is None or level_percentile_of(anchor) is None)
    # A percentile was taken on the reference lane's reading of the marker
    # species; once the share on the page is a part of that species (split)
    # or a rescaled reading, compare the two and say "not comparable" where
    # they disagree fivefold - the rank would otherwise describe a different
    # number from the one beside it.
    for rec in records:
        if rec.get("reference_reading") is None:
            rec["reference_conflict"] = False
            continue
        # like with like: the reference lane's *raw* reading (percent of its
        # classified reads, the same footing as the share), not the value
        # rescaled onto the cohort's frame for ranking
        source = str(rec.get("percentile_source") or "scoring cohort")
        lane = rec.get("reference_lane") or ("scoring" if source == "scoring cohort" else "")
        if lane == "scoring":
            raw = rec.get("scoring_percent")
        elif lane and lane == str(blob.get("primary_lane") or ""):
            # the composition lane's own reading: the share before it was
            # split or added to
            raw = rec.get("marker_percent") if rec.get("marker_percent") is not None else rec.get("percent")
            if rec.get("share_basis") == "absorbed":
                raw = 0.0
        elif rec.get("reference_lane_percent") is not None:
            raw = rec.get("reference_lane_percent")
        else:
            raw = rec.get("secondary_percent")
        if raw is None:
            raw = rec.get("reference_reading")
        rec["reference_conflict"] = bool(rec.get("in_primary")) and _reads_differently(
            float(rec.get("percent") or 0.0), float(raw or 0.0))
    blob["organisms"] = records
    blob["share_unification"] = {
        "n_absorbed": n_absorbed, "absorbed_percent_total": round(absorbed_total, 4),
        "secondary_scale_to_primary": round(scale, 3), "n_scale_pairs": len(ratios),
        "n_split": n_split, "n_estimated": n_estimated, "n_member": n_member,
        "estimated_percent_total": round(estimated_total, 4),
        "meaning": ("Every organism carries one share on the primary marker lane's scale: the lane's own reading; "
                    "a genus total divided among the populations the competitive whole-genome mapping tells apart "
                    "(split); a whole-genome reading rescaled for an organism the marker catalogue cannot place "
                    "(estimated, drawn from the unclassified band); or no share, because it is counted within a "
                    "relative's share above (member). A rejected call's share goes to the relative whose reads "
                    "they were (absorbed_percent on that relative)."),
    }
    return blob


def _reconcile_by_alias(merged: dict[str, dict[str, Any]]) -> None:
    """Merge records that are one organism under two names.

    The marker lanes name an organism by their own catalogue and carry its
    GTDB name as an alias; the genome lane names it by GTDB directly. Where
    a genome-lane record's key equals another record's GTDB alias - or a
    marker record's key equals a genome record's alias - they are the same
    organism and are folded into one, lanes and readings combined.

    Placeholder clusters need one more step. GTDB's ``sp`` codes follow the
    representative genome, which can change between releases, so the same
    unnamed organism can be ``Gemmiger sp937890665`` in one release and
    ``Gemmiger sp963557315`` in the next. Where exactly one unmerged
    placeholder from each side sits in the same genus they are treated as
    one organism; where several do, they are left apart and counted as the
    separate detections they might be.
    """
    from openbiota import gtdb as _gtdb

    def alias_of(rec: dict[str, Any]) -> str | None:
        if rec.get("sgb"):
            hit = _gtdb.lookup(rec["sgb"], rec.get("sgb_db") or _gtdb.DEFAULT_DB)
            if hit is not None and hit.resolved:
                return hit.display.replace(" ", "_")
        return (rec.get("gtdb") or "").replace(" ", "_") or None

    # Exact alias matches.
    alias_index: dict[str, str] = {}
    for key, rec in merged.items():
        a = alias_of(rec)
        if a and a != key:
            alias_index.setdefault(a, key)
    for key in list(merged):
        if key not in merged:
            continue
        target = alias_index.get(key)
        if target and target in merged and target != key:
            # the surviving key is the name the current catalogues use
            keep, other = (target, key) if _keep_rank(merged[target]) >= _keep_rank(merged[key]) else (key, target)
            _fold(merged, keep, other)

    # Records that resolve to the same GTDB species are one organism, however
    # many native identifiers they arrived under: two SGBs in a crosswalk
    # merge group, a mOTU and an SGB, a GlobDB representative and a Kraken
    # species. Folded into the record with the largest primary reading, so
    # the composition share is the one the primary lane measured most of.
    by_species: dict[str, list[str]] = {}
    for key, rec in merged.items():
        a = alias_of(rec)
        if a:
            by_species.setdefault(a, []).append(key)
    for keys in by_species.values():
        if len(keys) < 2:
            continue
        keys = [k for k in keys if k in merged]
        if len(keys) < 2:
            continue
        keep = max(keys, key=lambda k: _keep_rank(merged[k]))
        for k in keys:
            if k != keep and k in merged:
                _fold(merged, keep, k)

    # Records that name the same reference genome are one organism, whatever
    # each catalogue calls it: a mOTUs cluster and the GlobDB genome that
    # represents it, a panel genome and its GTDB accession.
    by_genome: dict[str, list[str]] = {}
    for key, rec in merged.items():
        for gid in rec.get("genome_ids") or ():
            by_genome.setdefault(gid, []).append(key)
    for keys in by_genome.values():
        keys = [k for k in keys if k in merged]
        if len(keys) < 2:
            continue
        keep = max(keys, key=lambda k: _keep_rank(merged[k]))
        for k in keys:
            if k != keep and k in merged:
                _fold(merged, keep, k)

    # Same-genus placeholder pairs, one each side.
    def is_placeholder(rec: dict[str, Any], key: str) -> bool:
        a = alias_of(rec) or key
        tail = a.split("_", 1)[-1] if "_" in a else ""
        return bool(rec.get("unnamed")) or (tail.startswith("sp") and tail[2:3].isdigit())

    genome_only = {k: r for k, r in merged.items() if r["lanes"] == {"genome"} and is_placeholder(r, k)}
    marker_only = {k: r for k, r in merged.items() if "genome" not in r["lanes"] and is_placeholder(r, k)}
    by_genus_g: dict[str, list[str]] = {}
    by_genus_m: dict[str, list[str]] = {}
    for k, r in genome_only.items():
        by_genus_g.setdefault(_genus_key(r, k), []).append(k)
    for k, r in marker_only.items():
        by_genus_m.setdefault(_genus_key(r, k), []).append(k)
    for genus, gk in by_genus_g.items():
        mk = by_genus_m.get(genus) or []
        if len(gk) == 1 and len(mk) == 1 and genus:
            _fold(merged, mk[0], gk[0])

    _fold_sibling_species(merged, alias_of)
    _separate_mislabels(merged, alias_of)


def _separate_mislabels(merged: dict[str, dict[str, Any]], alias_of: Any) -> None:  # noqa: ARG001 - same shape as the other reconciliation passes
    """An alias that names a different species is a mislabel, not a synonym.

    When two records fold, the folded record's key is kept as an alias so
    the organism can still be found under it. Most aliases are the same
    species under another catalogue's spelling. Some are not: an older
    catalogue labelled the bin with a species GTDB and the newer catalogue
    place elsewhere. That label must not answer a lookup for the species
    it names - a guild that asks for *Ruminococcus gnavus* must not be
    handed *Dorea hominis* - so it is moved to `misnamed`, kept for the
    record and the explorer, and dropped from the names the organism
    answers to.
    """
    from openbiota import gtdb as _gtdb
    from openbiota.expansion import names as _names

    try:
        bridge = _names.ncbi_species_to_r232()
    except Exception:  # noqa: BLE001 - without the bridge no alias can be judged
        return
    r232: set[str] = set()
    try:
        from openbiota.expansion import gtdb as _g

        r232 = set(_g.species_clusters("r232"))
    except Exception:  # noqa: BLE001
        pass

    def gtdb_spaced(label: str) -> str | None:
        """``Blautia_A_wexlerae_B`` -> ``Blautia_A wexlerae_B`` when that is an R232 species."""
        parts = label.split("_")
        if len(parts) < 2:
            return None
        if len(parts) > 2 and re.fullmatch(r"[A-Z]{1,2}", parts[1]):
            genus, rest = f"{parts[0]}_{parts[1]}", parts[2:]
        else:
            genus, rest = parts[0], parts[1:]
        if not rest:
            return None
        spaced = f"{genus} {'_'.join(rest)}"
        return spaced if spaced in r232 else None

    def species_named(label: str) -> str | None:
        return bridge.get(canonical(label)) or bridge.get(label) or gtdb_spaced(label)

    def own_species_of(rec: dict[str, Any]) -> str | None:
        if rec.get("sgb"):
            hit = _gtdb.lookup(rec["sgb"], rec.get("sgb_db") or _gtdb.DEFAULT_DB)
            if hit is not None and hit.resolved:
                return hit.species
        return rec.get("gtdb") or None

    for rec in merged.values():
        own_species = own_species_of(rec)
        # a complex (several species, or siblings folded as one) and a
        # placeholder species cannot convict a label of naming the wrong
        # species: there is no one species to compare it with
        if not own_species or own_species.endswith(")") or " complex (" in own_species or " / " in own_species \
                or _placeholder(own_species):
            continue
        keep: set[str] = set()
        for a in rec.get("aliases") or ():
            named = species_named(a)
            if named and not _placeholder(named) and named != own_species \
                    and _sibling_base(named) != _sibling_base(own_species):
                rec.setdefault("misnamed", set()).add(a)
            else:
                keep.add(a)
        rec["aliases"] = keep


def _distinct_from(reading: float, anchor_share: float) -> bool:
    """Whether a member's own reading places it as a population distinct from
    the relative it is counted within (see `SAME_POPULATION_FACTOR`)."""
    if reading <= 0 or anchor_share <= 0:
        return reading > 0
    ratio = reading / anchor_share
    return ratio < 1.0 / SAME_POPULATION_FACTOR or ratio > SAME_POPULATION_FACTOR


def _reads_differently(share: float | None, scoring: float | None) -> bool:
    """Fivefold disagreement between the composition share and the scoring
    catalogue's reading, both above trace: the catalogues do not measure the
    same thing under this name (Blautia wexlerae at 7% in Jan26 and 0.006%
    in MetaPhlAn 3), so a rank taken on one is not a statement about the other."""
    if not share or not scoring or share < 0.05 or scoring < 0.05:
        # one lane at half a percent or more, the other below trace: they disagree
        return bool(share and scoring and max(share, scoring) >= 0.5 and min(share, scoring) < 0.05)
    ratio = share / scoring
    return ratio > 5.0 or ratio < 0.2


def _keep_rank(rec: dict[str, Any]) -> tuple[Any, ...]:
    """Which of two records of one organism keeps its key when they fold.

    The name a reader sees is the surviving key, so it should be the one
    the current catalogues use: a record the scoring lane alone created
    carries a 2019 NCBI label and yields to any other; then the record
    with the larger readings, then the named one.
    """
    scoring_only = rec["lanes"] <= {"scoring"}
    # the marker catalogues' name is the one readers know (Agathobacter
    # rectalis, not GTDB's Roseburia rectalis); GTDB stays the alias. Between
    # two marker catalogues the newer one names the bin: Jun23 called SGB4571
    # "Ruminococcus gnavus", Jan26 and GTDB call it Dorea hominis.
    marker_rank = 2 if "jan26" in rec["lanes"] else (1 if "extended" in rec["lanes"] else 0)
    return (not scoring_only, marker_rank, not rec["unnamed"], len(rec["lanes"]), sum(rec["by_lane"].values()))


def _sibling_base(gtdb_name: str) -> str | None:
    """``Blautia_A wexlerae_B`` -> ``Blautia_A wexlerae``: the GTDB sibling stem.

    GTDB splits a classical species into ``X``, ``X_A``, ``X_B`` when its
    genomes fall in several clusters. The suffix is on the epithet; the
    genus keeps its own suffix (``Blautia_A`` is a genus, not a sibling).
    """
    parts = gtdb_name.replace("_", " ", 0).split(" ")
    if len(parts) < 2:
        return None
    genus, epithet = parts[0], " ".join(parts[1:])
    m = re.match(r"^(.+?)_[A-Z]{1,2}$", epithet)
    stem = m.group(1) if m else epithet
    return f"{genus} {stem}"


def _fold_sibling_species(merged: dict[str, dict[str, Any]], alias_of: Any) -> None:  # noqa: ARG001 - same shape as the other reconciliation passes
    """One population that two catalogues place in sibling GTDB species is one organism.

    A marker catalogue's bin is labelled by the majority of its member
    genomes (``wexlerae_B``); a whole-genome sketch places the same reads
    in the type cluster (``wexlerae``). When the two records are seen by
    disjoint sets of lanes - no lane saw both siblings - the sample holds
    one population the catalogues name differently, not two species. They
    are folded into one record and marked a complex, so the split is shown
    and never counted as multiplicity (spec 0.8.4 §5). Where any lane saw
    both siblings they stay apart: that lane could tell them apart.
    """
    from openbiota import gtdb as _gtdb

    def spaced(rec: dict[str, Any]) -> str | None:
        """The GTDB species with its space kept (``Blautia_A wexlerae_B``)."""
        if rec.get("sgb"):
            hit = _gtdb.lookup(rec["sgb"], rec.get("sgb_db") or _gtdb.DEFAULT_DB)
            if hit is not None and hit.resolved:
                return hit.display
        return rec.get("gtdb") or None

    groups: dict[str, list[str]] = {}
    for key, rec in merged.items():
        a = spaced(rec)
        if not a or rec.get("complex"):
            continue
        stem = _sibling_base(a)
        if stem:
            groups.setdefault(stem, []).append(key)
    for keys in groups.values():
        keys = [k for k in keys if k in merged]
        if len(keys) < 2:
            continue
        names = {spaced(merged[k]) for k in keys} - {None}
        if len(names) < 2:
            continue  # same species already handled above
        lane_sets = [merged[k]["lanes"] for k in keys]
        if any(lane_sets[i] & lane_sets[j] for i in range(len(keys)) for j in range(i + 1, len(keys))):
            continue  # a lane saw both: genuinely two populations
        keep = max(keys, key=lambda k: sum(merged[k]["by_lane"].values()))
        siblings = sorted(names)  # type: ignore[type-var]
        for k in keys:
            if k != keep:
                _fold(merged, keep, k)
        rec = merged[keep]
        rec["complex"] = True
        genus = siblings[0].split(" ")[0]
        rec["gtdb"] = f"{genus} " + " / ".join(n.split(" ", 1)[1] for n in siblings)
        rec["sibling_species"] = siblings


def _genus_key(rec: dict[str, Any], key: str) -> str:
    """The genus to pair placeholders on.

    A bare marker-catalogue bin carries a placeholder genus of its own
    (``GGB45596``) that says nothing; its GTDB alias carries the real one.
    So the alias's genus is preferred, then the record's own, then the key.
    GTDB's polyphyletic suffixes (``Blautia_A``) are kept: ``Blautia`` and
    ``Blautia_A`` are different genera in that taxonomy.
    """
    from openbiota import gtdb as _gtdb

    if rec.get("sgb"):
        hit = _gtdb.lookup(rec["sgb"], rec.get("sgb_db") or _gtdb.DEFAULT_DB)
        if hit is not None and hit.genus:
            return hit.genus
    alias = rec.get("gtdb") or ""
    if alias:
        return alias.split(" ")[0]
    g = rec.get("genus") or ""
    if g and not g.startswith("GGB"):
        return g
    return key.split("_")[0] if not key.startswith("GGB") else ""


def _fold(merged: dict[str, dict[str, Any]], into: str, other: str) -> None:
    a, b = merged[into], merged.pop(other)
    for lane, pct in b["by_lane"].items():
        a["by_lane"][lane] = a["by_lane"].get(lane, 0.0) + pct
    a["lanes"] |= b["lanes"]
    a.setdefault("aliases", set()).update({other, *(b.get("aliases") or ())})
    a.setdefault("genome_ids", set()).update(b.get("genome_ids") or ())
    a["genus"] = a["genus"] or b["genus"]
    if not a["sgb"]:
        a["sgb"], a["sgb_db"] = b["sgb"], b.get("sgb_db")
    a["gtdb"] = a["gtdb"] or b.get("gtdb")
    a["unnamed"] = a["unnamed"] and b["unnamed"]
    a["complex"] = a.get("complex") or b.get("complex")
    # Provenance is kept: a second native id in the same lane is appended.
    for lane, nid in (b.get("native_ids") or {}).items():
        cur = a.setdefault("native_ids", {}).get(lane)
        a["native_ids"][lane] = nid if not cur else (cur if nid in cur.split(";") else f"{cur};{nid}")
    for lane, stt in (b.get("lane_status") or {}).items():
        a.setdefault("lane_status", {}).setdefault(lane, stt)


#: Fewer reference carriers than this and a rank among them is not stated.
MIN_REFERENCE_CARRIERS: Final = 10
#: A population counted within a relative's share whose own whole-genome
#: reading is within this factor of that share (either way) is the same
#: population under another catalogue's name - GlobDB's "Phocaeicola
#: SPECIV4_34405" at 11.3% beside Phocaeicola vulgatus at 11.7% - and is not
#: judged a second time. One further from it is a distinct population the
#: competition told apart (Blautia luti at 0.39% inside Blautia wexlerae's
#: 7.2%), and its own lane rank stands.
SAME_POPULATION_FACTOR: Final = 2.0

#: The installed union the expansion is benchmarked against (spec §1):
#: MetaPhlAn 3 scoring, MetaPhlAn 4 Jun23, GTDB R232 sylph.
BASELINE_LANES: frozenset[str] = frozenset({"scoring", "extended", "genome"})


def _status_of(
    rec: Mapping[str, Any], lanes: Mapping[str, Lane], primary: str, percent: float,
) -> tuple[str, str, tuple[str, ...]]:
    """supported / provisional / ambiguous, and why.

    Two independent *methods* (not two lanes of one method) establish
    presence. One method does when its own reading is not marginal: the
    primary marker lane at 0.05% or more, or a whole-genome lane that
    itself called the hit supported. Anything resting on one marginal
    reading is provisional until competitive confirmation looks at it.
    """
    methods = tuple(sorted({lanes[n].method for n in rec["lanes"] if n in lanes}))
    if rec.get("complex"):
        if rec.get("sibling_species"):
            sibs = " and ".join(rec["sibling_species"])
            return ("ambiguous", f"one population that the catalogues place in sibling GTDB species ({sibs}); "
                    "no method saw both, so it is shown as one complex and never counted twice", methods)
        return "ambiguous", "the bin maps to several GTDB R232 species (split); shown as a complex", methods
    # A GlobDB-native cluster (no GTDB accession, no SGB, no mOTU behind it)
    # is a unit of one catalogue's clustering. Two methods that both read
    # GlobDB's genomes are not independent evidence for it: a cluster that
    # duplicates a GTDB species absorbs that species' reads in both. Such a
    # unit is provisional until competitive confirmation has had its genome.
    ids = rec.get("native_ids") or {}
    globdb_only = bool(ids) and set(ids) <= {"globdb", "singlem", "kraken", "rescue"} and not rec.get("sgb")
    if globdb_only:
        return "provisional", ("a GlobDB/UHGG cluster seen only by lanes that read GlobDB-derived references; "
                               "presence rests on competitive confirmation"), methods
    if len(methods) >= 2:
        return "supported", f"seen by {len(methods)} independent methods: {', '.join(methods)}", methods
    if methods == ("read_classification",):
        return "provisional", "read classification alone; a k-mer match is a candidate until competitive confirmation", methods
    if primary in rec["lanes"] and percent >= 0.05:
        return "supported", f"primary marker lane at {percent:.2f}% of classified reads", methods
    # One method that is not the primary marker lane establishes nothing on
    # its own, however strong its reading: a whole-genome sketch at 6% can
    # still be a sister species absorbing a relative's reads. Competitive
    # confirmation is what promotes it. Named by the lane with the largest
    # reading (a set has no first element; picking one made the status, and
    # with it the confirmation cache key, differ from run to run).
    only = max(sorted(rec["lanes"]), key=lambda n: float(rec["by_lane"].get(n) or 0.0))
    reading = rec["by_lane"].get(only) or 0.0
    return "provisional", (f"one method ({only}) at {reading:.2f}; presence rests on competitive confirmation"
                           if reading else f"one method ({only}); presence rests on competitive confirmation"), methods


def _count_category(rec: Mapping[str, Any], placement: Any, key: str) -> str:  # noqa: ARG001 - key kept for symmetry with _status_of
    if rec.get("complex") or (placement is not None and getattr(placement, "relationship", "") == "split"):
        return "unresolved_complex"
    if placement is not None and getattr(placement, "relationship", "") == "parent" and rec.get("unnamed"):
        return "higher_rank"
    if rec.get("unnamed"):
        return "unnamed_species_cluster"
    return "named_species"


def _accumulate(
    merged: dict[str, dict[str, Any]],
    key: str,
    lane: Lane,
    *,
    percent: float,
    genus: str = "",
    sgb: str | None = None,
    unnamed: bool = False,
    gtdb: str | None = None,
    sgb_db: str | None = None,
    native_ids: Mapping[str, str] | None = None,
    lane_status: str | None = None,
    complex_: bool = False,
    genome_id: str | None = None,
    aliases: tuple[str, ...] = (),
) -> None:
    """Fold one lane's reading of one organism into the merged record.

    Readings are kept per lane. The richer catalogue supplies genus and bin
    identity; abundance is resolved at build time against the primary lane.
    An SGB is only meaningful with the database it came from - SGB numbers
    are not stable across MetaPhlAn releases - so `sgb_db` travels with it.
    Every lane's native identifier is kept under the lane's name.
    """
    rec = merged.get(key)
    if rec is None:
        merged[key] = {
            "by_lane": {lane.name: percent}, "genus": genus, "sgb": sgb,
            "unnamed": unnamed, "lanes": {lane.name}, "gtdb": gtdb,
            "sgb_db": sgb_db if sgb else None, "native_ids": dict(native_ids or {}),
            "lane_status": {lane.name: lane_status} if lane_status else {}, "complex": complex_,
            "genome_ids": {genome_id} if genome_id else set(), "aliases": {a for a in aliases if a != key},
        }
        return
    rec["lanes"].add(lane.name)
    if genome_id:
        rec.setdefault("genome_ids", set()).add(genome_id)
    rec.setdefault("aliases", set()).update(a for a in aliases if a != key)
    if lane_status:
        rec.setdefault("lane_status", {})[lane.name] = lane_status
    rec["complex"] = rec.get("complex") or complex_
    for k, v in (native_ids or {}).items():
        rec["native_ids"].setdefault(k, v)
    # The SGB used for GTDB placement is the one from the newest MetaPhlAn
    # release present: SGB numbers are not stable across releases, and the
    # Jan26 bridge is the one walked with the least drift.
    if sgb and (not rec.get("sgb") or _release_rank(sgb_db) > _release_rank(rec.get("sgb_db"))):
        rec["sgb"], rec["sgb_db"] = sgb, sgb_db
        if gtdb:
            rec["gtdb"] = gtdb
        rec["complex"] = complex_
    # Each lane's reading is kept under its own name. They are separate
    # compositions with separate denominators and are never combined: the
    # primary lane's value is the share shown, and the others are kept as
    # secondary readings. Taking the larger of two would let the total
    # drift past 100.
    rec["by_lane"][lane.name] = rec["by_lane"].get(lane.name, 0.0) + percent
    rec["genus"] = rec["genus"] or genus
    rec["gtdb"] = rec["gtdb"] or gtdb
    # A bin named by any catalogue is a named organism.
    rec["unnamed"] = rec["unnamed"] and unnamed


_RELEASE_RANK = {"mpa_vJan26_CHOCOPhlAnSGB_202605": 3, "mpa_vJun23_CHOCOPhlAnSGB_202403": 2}


def _release_rank(db: str | None) -> int:
    return _RELEASE_RANK.get(db or "", 1)


def _genus_of(species: str) -> str:
    """Genus from a binomial, for lanes that do not report it separately."""
    head = species.split("_", 1)[0]
    return head if head and head[:1].isupper() else ""


def lanes_from_results(
    results: Mapping[str, Any],
    *,
    scoring_species: Mapping[str, float] | None = None,
) -> list[Lane]:
    """Assemble the lanes present in one sample's results.

    This is the single place that knows which profilers exist. Adding another
    means adding one :class:`Lane` here.
    """
    out: list[Lane] = []

    if scoring_species:
        # The catalogue the reference population was processed on, and so the
        # only one whose measurements can be given a percentile. Its names
        # are NCBI's of 2019; each is given the GTDB R232 species the
        # genomes of that name sit in, so the record merges with the newer
        # catalogues' record of the same organism instead of standing beside
        # it as a second detection.
        from openbiota.expansion import names as _names

        bridge = _names.ncbi_species_to_r232()
        rows = []
        for name, pct in scoring_species.items():
            gtdb_name = bridge.get(canonical(name)) or bridge.get(name)
            rows.append({"species": name, "percent": pct, "gtdb": gtdb_name,
                         "genus": (gtdb_name or name).replace("_", " ").split(" ")[0],
                         "unnamed": _mp3_unnamed(name), "native_id": name})
        out.append(Lane(
            name="scoring",
            species=dict(scoring_species),
            rankable=True,
            catalogue_size=13_500,
            rows=rows,
            method="marker", release="mpa_v31_CHOCOPhlAn_201901",
        ))

    detection = (results.get("detection") or {}).get("lanes") or {}
    jan26 = detection.get("metaphlan_jan26")
    has_jan26 = bool(jan26 and jan26.get("observations"))

    extended = results.get("extended_catalogue") or {}
    rows = extended.get("sgbs") or []
    if rows:
        species: dict[str, float] = {}
        for row in rows:
            name = str(row.get("species") or "")
            if name:
                species[name] = species.get(name, 0.0) + float(row.get("percent") or 0.0)
        from openbiota import gtdb as _gtdb
        rows = [{**row, "sgb_db": _gtdb.JUN23, "native_id": row.get("sgb")} for row in rows]
        out.append(Lane(
            name="extended",
            species=species,
            rankable=False,
            catalogue_size=int(extended.get("n_sgbs_in_database") or 36_822),
            rows=rows,
            # Jun23 is the composition when the upgraded lane did not run;
            # otherwise it is the preserved baseline and contributes detections.
            primary=not has_jan26,
            method="marker", release="mpa_vJun23_CHOCOPhlAnSGB_202403",
        ))

    if has_jan26:
        out.append(_jan26_lane(jan26))
    for lane_id, maker in (("sylph_globdb", _globdb_lane), ("motus4", _motus_lane),
                           ("kraken_uhgg", _kraken_lane), ("kraken_rescue", _kraken_rescue_lane),
                           ("singlem_globdb", _singlem_lane)):
        blob = detection.get(lane_id)
        if blob and blob.get("observations"):
            lane = maker(blob)
            if lane is not None:
                out.append(lane)

    genome = results.get("genome_profile") or {}
    hits = [h for h in (genome.get("hits") or []) if h.get("confident", True)]
    if hits:
        # Whole-genome containment against every GTDB species. Rows carry
        # the GTDB name as both species (underscored, as the catalogues
        # write names) and alias, so a hit merges with a marker-lane record
        # under either the current name or the GTDB one.
        rows = [{
            "species": str(h["species"]).replace(" ", "_"),
            "gtdb": str(h["species"]),
            "genus": str(h.get("genus") or ""),
            "percent": float(h.get("taxonomic_abundance") or 0.0),
            "unnamed": bool(h.get("placeholder")),
            "accession": h.get("accession"),
        } for h in hits]
        out.append(Lane(
            name="genome",
            species={r["species"]: r["percent"] for r in rows},
            rankable=False,
            catalogue_size=int(genome.get("catalogue_size") or 199_923),
            rows=rows,
            method="genome_sketch", release=str(genome.get("database") or "gtdb-r232"),
        ))

    return out


def _ncbi_style(species: str) -> str:
    return species.strip().replace(" ", "_")


def _mp3_unnamed(name: str) -> bool:
    """MetaPhlAn 3 labels for organisms known only from assemblies:
    ``Blautia_sp_CAG_257``, ``Firmicutes_bacterium_CAG_424``,
    ``Clostridiales_bacterium_S5_A14a``, ``Candidatus_..._bacterium_HUM_18``."""
    parts = name.split("_")
    return "_sp_" in f"_{name}_" or "bacterium" in parts or "unclassified" in name.lower() or name.endswith("_sp")


def _jan26_lane(blob: Mapping[str, Any]) -> Lane:
    """MetaPhlAn Jan26 SGBs as the primary composition, with R232 aliases."""
    from openbiota import gtdb as _gtdb

    rows: list[dict[str, Any]] = []
    species: dict[str, float] = {}
    for o in blob.get("observations") or []:
        if o.get("rank") != "species" or not o.get("abundance_value"):
            continue
        segs = str(o.get("lineage") or "").split(";")
        name = next((x[3:] for x in segs if x.startswith("s__")), "")
        genus = next((x[3:] for x in segs if x.startswith("g__")), "")
        if not name:
            continue
        sgb = str(o.get("native_id") or "")
        hit = _gtdb.lookup(sgb, _gtdb.JAN26)
        unnamed = "_SGB" in name or bool(o.get("support_metrics", {}).get("unnamed_sgb"))
        rows.append({
            "species": name, "genus": genus, "percent": float(o["abundance_value"]), "sgb": sgb, "sgb_db": _gtdb.JAN26,
            "native_id": sgb, "unnamed": unnamed,
            "gtdb": (hit.species if (hit is not None and hit.resolved) else ""),
            "complex": bool(hit is not None and hit.relationship == "split"),
            "status": str(o.get("status") or "supported"),
        })
        species[name] = species.get(name, 0.0) + float(o["abundance_value"])
    return Lane(name="jan26", species=species, rankable=False, catalogue_size=72_000, rows=rows, primary=True,
                method="marker", release=str(blob.get("reference_release_id") or "mpa_vJan26_CHOCOPhlAnSGB_202605"))


def _lineage_parts(lineage: str) -> tuple[str, str]:
    segs = [x.strip() for x in str(lineage).replace("|", ";").split(";")]
    species = next((x[3:] for x in segs if x.startswith("s__")), "").strip()
    genus = next((x[3:] for x in segs if x.startswith("g__")), "").strip()
    return species, genus


_ID_TOKEN = re.compile(r"^(?:[A-Z][A-Z0-9-]{2,}_[A-Za-z0-9.-]*\d[A-Za-z0-9.-]*|[A-Z]{3,}\d{5,}|UBA\d+|CAG-\d+|sp\d{6,})$")


def _cluster_genome(species: str) -> str | None:
    """A GlobDB-native species is named by its representative genome
    (``Phocaeicola SPECIV4_00061``): that genome, when the archive holds it."""
    tail = species.split(" ", 1)[1] if " " in species else species
    if not tail or not _ID_TOKEN.match(tail) or tail.startswith("sp"):
        return None
    from openbiota.expansion import genomes as _genomes

    return tail if _genomes.in_globdb(tail) else None


def _placeholder(species: str) -> bool:
    """Whether a species label is a catalogue identifier rather than a name.

    ``Gemmiger sp937890665``, ``Blautia_A MGYG000001338``, ``Lachnospira
    MOTU40_010844``, ``Phocaeicola SPECIV4_34405``, ``MOTU40_058831``
    (no genus at all) and the profilers' ``Unknown Eggerthella`` are
    placeholders; ``Blautia_A wexlerae_B`` is a name.
    """
    species = species.strip()
    if not species or species.lower().startswith("unknown"):
        return True
    head, _, tail = species.partition(" ")
    if not tail:
        return bool(_ID_TOKEN.match(head))
    return bool(_ID_TOKEN.match(tail))


def _globdb_lane(blob: Mapping[str, Any]) -> Lane | None:
    rows: list[dict[str, Any]] = []
    for o in blob.get("observations") or []:
        species, genus = _lineage_parts(o.get("lineage", ""))
        if o.get("rank") != "species" or not species:
            continue
        accession = str(o.get("native_id") or "").split(":", 1)[-1]
        head, _, tail = species.partition(" ")
        row_aliases: list[str] = []
        if tail and tail == head:
            # ``MOTU40_058831 MOTU40_058831``: GlobDB's genus placeholder is
            # the cluster id itself; the organism is the id, genus unknown
            row_aliases.append(_ncbi_style(species))
            species, genus = head, ""
        rows.append({
            "species": _ncbi_style(species), "gtdb": species, "genus": genus,
            "percent": float(o.get("abundance_value") or 0.0), "unnamed": _placeholder(species),
            "native_id": str(o.get("native_id") or ""), "accession": accession,
            "genome_id": accession or None, "aliases": row_aliases,
            "status": str(o.get("status") or "supported"),
        })
    if not rows:
        return None
    return Lane(name="globdb", species={r["species"]: r["percent"] for r in rows}, rankable=False, catalogue_size=346_233,
                rows=rows, method="genome_sketch", release=str(blob.get("reference_release_id") or "GlobDB r232"))


def _motus_lane(blob: Mapping[str, Any]) -> Lane | None:
    # mOTUs 4.1 names are GTDB R226; walk them to R232 by genome membership.
    try:
        from openbiota.expansion.crosswalk import release_species_map
        r226 = release_species_map("r226")
    except Exception:  # noqa: BLE001 - without the map the R226 name stands, flagged by its release
        r226 = {}
    from openbiota.expansion import names as _names

    rows: list[dict[str, Any]] = []
    for o in blob.get("observations") or []:
        species, genus = _lineage_parts(o.get("lineage", ""))
        if not species:
            continue
        motu = str(o.get("native_id") or "")
        # The profiler's label is a majority vote over the cluster's members
        # and reads ``Unknown <genus>`` when they disagree. The cluster's
        # representative genome has a GTDB species, and GlobDB may hold that
        # very genome: naming by it makes one record where the sylph lane
        # and this one saw the same unit.
        label = species
        resolved = _names.resolve_motu(motu, species, genus)
        species, genus = str(resolved["species"]), str(resolved["genus"] or "")
        walked = r226.get(species)
        if walked:
            species = walked
            genus = walked.split(" ", 1)[0]
        rows.append({
            "species": _ncbi_style(species), "gtdb": species, "genus": genus,
            "percent": float(o.get("abundance_value") or 0.0),
            "unnamed": bool(resolved["unnamed"]) or _placeholder(species) or o.get("rank") != "species",
            "native_id": motu, "genome_id": resolved.get("genome_id"),
            "aliases": [_ncbi_style(label)] if _ncbi_style(label) != _ncbi_style(species) else [],
            "status": str(o.get("status") or "supported"),
        })
    if not rows:
        return None
    return Lane(name="motus", species={r["species"]: r["percent"] for r in rows}, rankable=False, catalogue_size=124_295,
                rows=rows, method="universal_marker", release=str(blob.get("reference_release_id") or "mOTUs DB 4.1"))


def _kraken_lane(blob: Mapping[str, Any]) -> Lane | None:
    """Kraken/UHGG species, resolved through the UHGG reconciliation where it exists."""
    recon: dict[str, tuple[str, str, str]] = {}
    try:
        import csv as _csv
        import gzip as _gzip
        path = Path("refs/expanded/reconciliation/uhgg_v2.0.2.tsv.gz")
        if path.is_file():
            with _gzip.open(path, "rt") as fh:
                for row in _csv.DictReader(fh, delimiter="\t"):
                    recon[row["source_id"]] = (row["category"], row["canonical_species"], row["globdb_id"])
    except OSError:
        recon = {}
    rows: list[dict[str, Any]] = []
    for o in blob.get("observations") or []:
        species, genus = _lineage_parts(o.get("lineage", ""))
        if o.get("rank") != "species" or not species:
            continue
        if (o.get("support_metrics") or {}).get("below_inventory_floor"):
            continue  # retained in the lane record; not an organism
        rep = str((o.get("support_metrics") or {}).get("uhgg_species_rep") or "")
        cat, canonical, gid = recon.get(rep, ("", "", ""))
        if cat == "already_represented" and canonical:
            name, alias, unnamed = _ncbi_style(canonical), canonical, _placeholder(canonical)
        elif cat == "new_cluster" and gid:
            name, alias, unnamed = _ncbi_style(f"{genus} {gid}"), f"{genus} {gid}", True
        else:
            name, alias, unnamed = _ncbi_style(species), species, _placeholder(species)
        rows.append({
            "species": name, "gtdb": alias, "genus": genus, "percent": float(o.get("abundance_value") or 0.0),
            "unnamed": unnamed, "native_id": str(o.get("native_id") or ""),
            # the UHGG species representative (MGYG...) is the genome the
            # reads were classified to; the GlobDB cluster it maps to, if any
            "genome_id": (gid if (cat == "new_cluster" and gid) else rep) or None,
            "status": str(o.get("status") or "supported"),
        })
    if not rows:
        return None
    return Lane(name="kraken", species={r["species"]: r["percent"] for r in rows}, rankable=False, catalogue_size=4_744,
                rows=rows, method="read_classification", release=str(blob.get("reference_release_id") or "UHGG v2.0.2"))


_PANEL_TARGETS: dict[str, dict[str, str]] | None = None


def _panel_targets() -> dict[str, dict[str, str]]:
    """Rescue-panel source genome -> manifest row (GlobDB id, local path)."""
    global _PANEL_TARGETS
    if _PANEL_TARGETS is None:
        _PANEL_TARGETS = {}
        path = Path("refs/kraken2/rescue_panel/panel_manifest.tsv")
        if path.is_file():
            import csv

            with path.open() as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    sid = str(row.get("source_id") or "").rstrip("_")
                    if sid:
                        _PANEL_TARGETS[sid] = dict(row)
                        for suffix in (".fna.gz", ".fa.gz", ".fna", ".fa"):
                            if sid.endswith(suffix):
                                _PANEL_TARGETS.setdefault(sid[: -len(suffix)], dict(row))
    return _PANEL_TARGETS


def _kraken_rescue_lane(blob: Mapping[str, Any]) -> Lane | None:
    """Kraken against the locally built gut rescue panel (spec §4B).

    Lineages are GTDB/GlobDB already, so the species name is the alias.
    Background decoys (``d__Background``) are competitors: a hit on one is
    kept in the lane record and is never an organism. The panel genome
    behind each taxid is the native identifier the confirmation can fetch.
    Same method as the UHGG lane (read classification): the two Kraken
    lanes together count as one method for the independence rule.
    """
    rows: list[dict[str, Any]] = []
    n_decoy_reads = 0
    for o in blob.get("observations") or []:
        lineage = str(o.get("lineage", ""))
        species, genus = _lineage_parts(lineage)
        if o.get("rank") != "species" or not species:
            continue
        metrics = o.get("support_metrics") or {}
        if lineage.startswith("d__Background"):
            n_decoy_reads += int(metrics.get("kraken_clade_reads") or 0)
            continue
        if metrics.get("below_inventory_floor"):
            continue
        genome = str(metrics.get("panel_genome") or "").rstrip("_")
        target = _panel_targets().get(genome) or _panel_targets().get(genome.rsplit(".", 1)[0]) or {}
        genome_id = target.get("globdb_id") or (f"file:{target['local_path']}" if target.get("local_path") else None)
        rows.append({
            "species": _ncbi_style(species), "gtdb": species, "genus": genus,
            "percent": float(o.get("abundance_value") or 0.0), "unnamed": _placeholder(species),
            "native_id": (f"panel:{genome}" if genome else str(o.get("native_id") or "")),
            "accession": genome, "genome_id": genome_id, "status": str(o.get("status") or "supported"),
        })
    if not rows:
        return None
    n_targets = int((blob.get("summary") or {}).get("panel_species") or 0)
    return Lane(name="rescue", species={r["species"]: r["percent"] for r in rows}, rankable=False,
                catalogue_size=n_targets or len(rows), rows=rows, method="read_classification",
                release=str(blob.get("reference_release_id") or "gut rescue panel"))


def _singlem_lane(blob: Mapping[str, Any]) -> Lane | None:
    rows: list[dict[str, Any]] = []
    for o in blob.get("observations") or []:
        species, genus = _lineage_parts(o.get("lineage", ""))
        if o.get("rank") != "species" or not species:
            continue
        rows.append({
            "species": _ncbi_style(species), "gtdb": species, "genus": genus, "percent": None,
            "unnamed": _placeholder(species), "native_id": str(o.get("native_id") or ""),
            "genome_id": _cluster_genome(species),
            "status": str(o.get("status") or "supported"),
        })
    if not rows:
        return None
    return Lane(name="singlem", species={}, rankable=False, catalogue_size=346_233, rows=rows,
                method="marker_window", detection_only=True, unit="coverage",
                release=str(blob.get("reference_release_id") or "GlobDB_r232.metapackage_v4"))


def from_results(
    results: Mapping[str, Any],
    *,
    scoring_species: Mapping[str, float] | None = None,
    ranked_rows: Iterable[Any] = (),
) -> Inventory:
    """Build the inventory for one sample straight from its results."""
    extended = results.get("extended_catalogue") or {}
    unclassified = extended.get("unclassified_percent")
    jan26 = ((results.get("detection") or {}).get("lanes") or {}).get("metaphlan_jan26") or {}
    if jan26.get("observations"):
        unclassified = (jan26.get("summary") or {}).get("unclassified_percent", unclassified)
    return build(
        lanes_from_results(results, scoring_species=scoring_species),
        ranked_rows=ranked_rows,
        strain=results.get("strain_resolution"),
        unclassified_percent=float(unclassified) if unclassified is not None else None,
    )


__all__ = [
    "Inventory",
    "Lane",
    "Organism",
    "build",
    "canonical",
    "from_json",
    "former_name",
    "from_results",
    "lanes_from_results",
]
