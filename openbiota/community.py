"""Community overview: curated taxon groups and the full species table.

Two readings that the disease profiles do not produce but the report needs:

* **Taxon groups** — each curated set in ``taxa/`` (butyrate producers,
  oral-origin organisms, opportunistic pathogens …) summed to one abundance
  and placed against the matched reference cohort. These are the "microbial
  groups" a consumer report leads with, and the same sums feed
  ``CARRIER_ABUNDANCE`` features in the profiles so the two never disagree.
* **Species table** — every species the profiler named, with its abundance,
  its percentile against the matched cohort, and how common it is in that
  cohort. This is the complete inventory; the groups and the profiles are
  views onto it.

Both use the cohort's prevalence-aware percentile so that "not detected" is a
category (spec 5.3), not a value at the zero floor. A species that 3% of the
reference cohort carries is not "at the 50th percentile" when absent; it is
absent, like most people.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from openbiota.refcohort import ReferenceCohort, build_taxon_references, sample_clr
from openbiota.scoring import FeatureReference, summarise_reference
from openbiota.similarity import (
    _frame_mean_log,
    aggregate_log_ratio,
    cohort_group_series,
    genus_of,
    rescale_to_classified,
)
from openbiota.taxongroups import TaxonGroup, TaxonGroupSet

#: Species at or below this relative abundance are listed but flagged as
#: trace: their percentile is dominated by sampling noise at typical depth.
TRACE_PERCENT = 0.01


@dataclass(frozen=True, slots=True)
class MemberReading:
    """One species inside a group: its own abundance and cohort position."""

    species: str
    percent: float
    detected: bool
    percentile: float | None
    cohort_prevalence: float | None
    in_catalogue: bool
    #: Which detection lanes saw this member (from the organism inventory).
    detected_by: tuple[str, ...] = ()
    #: True when only the expanded lanes saw it - the scoring catalogue did not.
    expansion_only: bool = False
    #: Rank among reference people who carry the species, on abundance:
    #: the percentile the report prints beside the level, because it and the
    #: deviation from the typical carrier are taken on the same distribution
    #: and so always point the same way. `percentile` above is prevalence-
    #: aware (absence is a tied category) and is kept for scoring.
    carrier_percentile: float | None = None
    reference_median: float | None = None
    #: how many reference people carry the species: the carrier rank rests on them
    reference_carriers: int | None = None

    @property
    def level_percentile(self) -> float | None:
        if self.detected and self.carrier_percentile is not None:
            return self.carrier_percentile
        return self.percentile

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "percent": round(self.percent, 4),
            "detected": self.detected,
            "detected_by": list(self.detected_by),
            "expansion_only": self.expansion_only,
            "percentile": None if self.percentile is None else round(self.percentile, 1),
            "carrier_percentile": None if self.carrier_percentile is None else round(self.carrier_percentile, 1),
            "reference_median": None if self.reference_median is None else round(self.reference_median, 4),
            "reference_carriers": self.reference_carriers,
            "cohort_prevalence": (
                None if self.cohort_prevalence is None else round(self.cohort_prevalence, 3)
            ),
            "in_catalogue": self.in_catalogue,
        }


@dataclass(frozen=True, slots=True)
class GroupReading:
    """A taxon group's summed abundance placed against the reference cohort."""

    group: TaxonGroup
    percent: float
    detected: bool
    percentile: float | None
    #: Fraction of the reference cohort in which the group sum was non-zero.
    cohort_prevalence: float | None
    members: tuple[MemberReading, ...]
    #: Members the profiler's catalogue can name — the rest cannot be measured.
    n_resolvable: int
    reference_n: int
    reference_source: str
    #: The group total on the reference catalogue's own namespace - the number
    #: the percentile was computed from. `percent` is the composition total
    #: from the inventory's primary lane once enriched; before enrichment the
    #: two are the same number.
    scoring_percent: float | None = None
    #: Members detected only by the expanded lanes.
    n_detected_expansion_only: int = 0
    #: Median of the group total across reference people in whom the group
    #: is present at all, in the scoring lane's units; compared with
    #: `scoring_percent` (the same lane), not with the composition total.
    reference_median: float | None = None
    #: Rank of the group total among reference people in whom the group is
    #: present, on the scoring lane's raw totals (see SpeciesRow).
    carrier_percentile: float | None = None
    reference_carriers: int | None = None

    @property
    def level_percentile(self) -> float | None:
        """The percentile that describes the level: rank among reference people
        who have the group, when the sample has it; the prevalence-aware rank
        otherwise (an absent group is placed among those who also lack it)."""
        if self.detected and self.carrier_percentile is not None:
            return self.carrier_percentile
        return self.percentile

    @property
    def n_detected(self) -> int:
        return sum(1 for m in self.members if m.detected)

    @property
    def deviation_percent(self) -> float | None:
        reading = self.scoring_percent if self.scoring_percent is not None else self.percent
        if not self.reference_median or reading is None or reading <= 0:
            return None
        return (reading / self.reference_median - 1.0) * 100.0

    @property
    def resolved_fraction(self) -> float:
        total = len(self.group.members) + len(self.group.unresolved)
        return self.n_resolvable / total if total else 0.0

    @property
    def dominant(self) -> tuple[MemberReading, ...]:
        """Members contributing to the sum, largest first."""
        return tuple(sorted((m for m in self.members if m.detected), key=lambda m: -m.percent))

    def to_json(self) -> dict[str, Any]:
        return {
            "group": self.group.name,
            "label": self.group.label,
            "category": self.group.category,
            "higher_means": self.group.higher_means,
            "reader_direction": self.group.reader_direction,
            "direction_note": self.group.direction_note,
            "percent": round(self.percent, 4),
            "detected": self.detected,
            "percentile": None if self.percentile is None else round(self.percentile, 1),
            "cohort_prevalence": (
                None if self.cohort_prevalence is None else round(self.cohort_prevalence, 3)
            ),
            "n_members": len(self.group.members),
            "n_resolvable": self.n_resolvable,
            "n_unresolved_in_literature": len(self.group.unresolved),
            "n_detected": self.n_detected,
            "n_detected_expansion_only": self.n_detected_expansion_only,
            "scoring_percent": None if self.scoring_percent is None else round(self.scoring_percent, 4),
            "reference_median": None if self.reference_median is None else round(self.reference_median, 4),
            "deviation_percent": None if self.deviation_percent is None else round(self.deviation_percent, 1),
            "carrier_percentile": None if self.carrier_percentile is None else round(self.carrier_percentile, 1),
            "reference_carriers": self.reference_carriers,
            "members": [m.to_json() for m in self.members],
            "reference_n": self.reference_n,
            "reference_source": self.reference_source,
        }


@dataclass(frozen=True, slots=True)
class SpeciesRow:
    """One row of the full species inventory."""

    species: str
    genus: str
    percent: float
    percentile: float | None
    cohort_prevalence: float | None
    in_catalogue: bool
    #: Which curated groups this species belongs to, for the report's cross-reference.
    groups: tuple[str, ...] = ()
    #: The typical level among reference people who carry the species: the
    #: median of their abundances, in the same units as `percent`. The
    #: percentile says where the sample sits in the distribution; this says
    #: how far, as a ratio a reader can picture ("five times the typical
    #: carrier", "a tenth of it").
    reference_median: float | None = None
    #: Rank of `percent` among the carriers' abundances (midrank, 0-100):
    #: the same distribution `reference_median` is the middle of, so a
    #: carrier percentile above 50 always goes with a positive deviation and
    #: below 50 with a negative one. `percentile` is the prevalence-aware
    #: rank (absence as a tied category), kept for the scoring machinery.
    carrier_percentile: float | None = None
    #: how many reference people carry the species
    reference_carriers: int | None = None

    @property
    def trace(self) -> bool:
        return self.percent <= TRACE_PERCENT

    @property
    def deviation_percent(self) -> float | None:
        """Percent above (+) or below (-) the typical carrier level; None without a reference."""
        if not self.reference_median or self.percent <= 0:
            return None
        return (self.percent / self.reference_median - 1.0) * 100.0

    @property
    def rare_in_cohort(self) -> bool:
        return self.cohort_prevalence is not None and self.cohort_prevalence < 0.10

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "genus": self.genus,
            "percent": round(self.percent, 4),
            "percentile": None if self.percentile is None else round(self.percentile, 1),
            "cohort_prevalence": (
                None if self.cohort_prevalence is None else round(self.cohort_prevalence, 3)
            ),
            "in_catalogue": self.in_catalogue,
            "trace": self.trace,
            "groups": list(self.groups),
            "reference_median": None if self.reference_median is None else round(self.reference_median, 4),
            "deviation_percent": None if self.deviation_percent is None else round(self.deviation_percent, 1),
            "carrier_percentile": None if self.carrier_percentile is None else round(self.carrier_percentile, 1),
            "reference_carriers": self.reference_carriers,
        }


def carrier_median(carriers: Sequence[float]) -> float | None:
    """The typical carrier's level: the median of the sorted carrier readings.

    Interpolated for an even count, so the median and `carrier_rank` describe
    the same distribution. Taking the upper middle value instead put a sample
    that sat at the lower middle at the 50th percentile and 24% *below* "the
    typical carrier", which reads as a contradiction on the page.
    """
    if not carriers:
        return None
    n = len(carriers)
    mid = n // 2
    return carriers[mid] if n % 2 else (carriers[mid - 1] + carriers[mid]) / 2.0


def carrier_rank(carriers: Sequence[float], value: float) -> float | None:
    """Midrank percentile of `value` among the sorted abundances of reference
    people who carry the feature; None when nobody in the reference does."""
    if not carriers or value <= 0:
        return None
    import bisect

    below = bisect.bisect_left(carriers, value)
    ties = bisect.bisect_right(carriers, value) - below
    return (below + 0.5 * ties) / len(carriers) * 100.0


@dataclass(slots=True)
class CommunityOverview:
    """Taxon-group readings and the species table for one sample."""

    groups: list[GroupReading]
    species: list[SpeciesRow]
    reference_n: int
    reference_source: str
    notes: list[str] = field(default_factory=list)

    @property
    def n_species_detected(self) -> int:
        return sum(1 for s in self.species if s.percent > 0)

    @property
    def n_genera_detected(self) -> int:
        return len({s.genus for s in self.species if s.percent > 0})

    @property
    def n_not_in_catalogue(self) -> int:
        return sum(1 for s in self.species if not s.in_catalogue)

    def group(self, name: str) -> GroupReading | None:
        return next((g for g in self.groups if g.group.name == name), None)

    def by_category(self) -> dict[str, list[GroupReading]]:
        out: dict[str, list[GroupReading]] = {}
        for g in self.groups:
            out.setdefault(g.group.category, []).append(g)
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "n_species_detected": self.n_species_detected,
            "n_genera_detected": self.n_genera_detected,
            "n_species_not_in_reference_catalogue": self.n_not_in_catalogue,
            "reference_n": self.reference_n,
            "reference_source": self.reference_source,
            "groups": [g.to_json() for g in self.groups],
            "species": [s.to_json() for s in self.species],
            "notes": list(self.notes),
            "what_this_is": (
                "Every species the taxonomic profiler named in this sample, with its relative "
                "abundance and where that abundance sits against the matched reference cohort. "
                "Groups are curated sums of species that share a function or origin. A "
                "percentile is a position among reference samples, not a health verdict; "
                "the direction the literature reads it in is given per group."
            ),
        }


# --------------------------------------------------------------------------- #
# construction
# --------------------------------------------------------------------------- #


def _group_reference(
    cohort: ReferenceCohort, members: Sequence[str], source: str
) -> FeatureReference | None:
    present = [m for m in members if cohort.index_of(m) is not None]
    if not present:
        return None
    series, detected = cohort_group_series(cohort, present)
    return summarise_reference(
        "group", series, source,
        prevalence=sum(detected) / max(1, len(detected)),
        detected=detected,
    )


def _percentile(reference: FeatureReference | None, value: float, detected: bool) -> float | None:
    """Prevalence-aware where the reference is zero-inflated, rank otherwise."""
    if reference is None or not reference.n:
        return None
    if reference.zero_inflated:
        return reference.prevalence_aware_percentile(value, detected=detected)
    if not detected and reference.n_absent:
        return reference.prevalence_aware_percentile(value, detected=False)
    return reference.percentile_of(value)


def build_community_overview(
    *,
    species_percent: Mapping[str, float],
    cohort: ReferenceCohort | None,
    group_set: TaxonGroupSet,
    source_label: str = "",
) -> CommunityOverview:
    """Place the sample's species and each curated group against the cohort.

    ``species_percent`` is the profiler's species-level output (percent,
    classified fraction). ``cohort`` should be the *matched* cohort — the same
    one the profiles were scored against — so the percentiles here and the
    percentiles in the profile pages come from the same reference samples.
    """
    species = rescale_to_classified(species_percent)
    notes: list[str] = []

    if cohort is None:
        rows = [
            SpeciesRow(s, genus_of(s), v, None, None, False, _groups_of(s, group_set))
            for s, v in sorted(species.items(), key=lambda kv: -kv[1])
            if v > 0
        ]
        groups = [
            GroupReading(
                group=g,
                percent=sum(species.get(m, 0.0) for m in g.members),
                detected=any(species.get(m, 0.0) > 0 for m in g.members),
                percentile=None,
                cohort_prevalence=None,
                members=tuple(
                    MemberReading(m, species.get(m, 0.0), species.get(m, 0.0) > 0, None, None, False)
                    for m in g.members
                ),
                n_resolvable=0,
                reference_n=0,
                reference_source="",
            )
            for g in group_set.groups
        ]
        notes.append("No reference cohort; abundances are reported without percentiles.")
        return CommunityOverview(groups, rows, 0, "", notes)

    source = source_label or f"reference cohort (n={cohort.n_samples:,})"
    taxon_refs = build_taxon_references(cohort)
    clr_values = sample_clr(species, cohort)
    catalogue = set(cohort.taxa)

    # Full species table: everything detected, plus nothing else — an absent
    # species is not a row, it is the default.
    rows: list[SpeciesRow] = []
    for name, value in sorted(species.items(), key=lambda kv: -kv[1]):
        if value <= 0:
            continue
        ref = taxon_refs.get(name)
        if ref is None:
            rows.append(SpeciesRow(name, genus_of(name), value, None, None, False, _groups_of(name, group_set)))
            continue
        raw = cohort.abundance[cohort.index_of(name)]  # type: ignore[index]
        fr = FeatureReference(
            feature=name, median=ref.median, mad=ref.mad, n=ref.n, source=source,
            values=ref.values, prevalence=ref.prevalence,
            present_values=tuple(
                v for v, j in zip(ref.values, range(ref.n), strict=True)
                if raw[j] > 0
            ),
        )
        carriers = sorted(float(v) for v in raw if v > 0)
        rows.append(
            SpeciesRow(
                name, genus_of(name), value,
                _percentile(fr, clr_values[name], detected=True),
                ref.prevalence, True, _groups_of(name, group_set),
                reference_median=carrier_median(carriers),
                carrier_percentile=carrier_rank(carriers, value),
                reference_carriers=len(carriers),
            )
        )
    n_outside = sum(1 for r in rows if not r.in_catalogue)
    if n_outside:
        notes.append(
            f"{n_outside} detected species are not in the reference cohort's catalogue and "
            "are listed without a percentile."
        )

    # Groups: sum of members, log-ratio transformed on the cohort's frame,
    # placed against the same transform of the cohort.
    vector = [max(0.0, species.get(t, 0.0)) for t in cohort.taxa]
    frame_mean = _frame_mean_log(vector, cohort.frame)
    groups: list[GroupReading] = []
    for g in group_set.groups:
        resolvable = [m for m in g.members if m in catalogue]
        total = sum(species.get(m, 0.0) for m in g.members)
        detected = total > 0
        ref = _group_reference(cohort, resolvable, source)
        pct = _percentile(ref, aggregate_log_ratio(total, frame_mean), detected) if ref else None
        members = []
        for m in g.members:
            mv = species.get(m, 0.0)
            tr = taxon_refs.get(m)
            m_pct = None
            if tr is not None:
                row_i = cohort.index_of(m)
                fr = FeatureReference(
                    feature=m, median=tr.median, mad=tr.mad, n=tr.n, source=source,
                    values=tr.values, prevalence=tr.prevalence,
                    present_values=tuple(
                        v for v, j in zip(tr.values, range(tr.n), strict=True)
                        if cohort.abundance[row_i][j] > 0  # type: ignore[index]
                    ),
                )
                m_pct = _percentile(fr, clr_values.get(m, 0.0), detected=mv > 0)
            m_carriers = (sorted(float(v) for v in cohort.abundance[cohort.index_of(m)] if v > 0)  # type: ignore[index]
                          if tr is not None else [])
            members.append(
                MemberReading(
                    m, mv, mv > 0, m_pct,
                    None if tr is None else tr.prevalence,
                    m in catalogue,
                    carrier_percentile=carrier_rank(m_carriers, mv),
                    reference_median=carrier_median(m_carriers),
                    reference_carriers=(len(m_carriers) if tr is not None else None),
                )
            )
        group_median = None
        group_carrier_pct = None
        if resolvable:
            idx = [cohort.index_of(m) for m in resolvable]
            totals = sorted(
                t for t in (sum(float(cohort.abundance[i][j]) for i in idx) for j in range(cohort.n_samples))  # type: ignore[index]
                if t > 0)
            group_median = carrier_median(totals)
            group_carrier_pct = carrier_rank(totals, total)
        groups.append(
            GroupReading(
                group=g,
                percent=total,
                detected=detected,
                percentile=pct,
                cohort_prevalence=None if ref is None else ref.prevalence,
                members=tuple(members),
                n_resolvable=len(resolvable),
                reference_n=0 if ref is None else ref.n,
                reference_source=source,
                reference_median=group_median,
                carrier_percentile=group_carrier_pct,
                reference_carriers=(len(totals) if resolvable else None),
            )
        )
        if not resolvable:
            notes.append(f"Group '{g.name}': none of its members are in the reference catalogue.")

    return CommunityOverview(groups, rows, cohort.n_samples, source, notes)


def _groups_of(species: str, group_set: TaxonGroupSet) -> tuple[str, ...]:
    return tuple(g.name for g in group_set.groups if species in g.members)


def group_carrier_map(group_set: TaxonGroupSet) -> dict[str, tuple[str, ...]]:
    """Taxon groups as carrier sets, so a profile's ``CARRIER_ABUNDANCE``
    feature may name a curated group (``oral_origin``) as well as a panel."""
    return {g.name: g.members for g in group_set.groups}


__all__ = [
    "CommunityOverview",
    "GroupReading",
    "MemberReading",
    "SpeciesRow",
    "build_community_overview",
    "group_carrier_map",
]


def enrich_with_inventory(overview: CommunityOverview, inventory: Any) -> CommunityOverview:
    """Second phase: members' detection and abundance from every lane.

    The percentile of a group is computed on the reference catalogue's own
    namespace and is left exactly as it was. What changes is what the group
    *contains*: a member the scoring catalogue could not see but another
    lane did is now detected, with its composition share from the
    inventory's primary lane, and the group total shown is the sum of those
    shares - so the members listed under a group add up to the number above
    them. The scoring-namespace total moves to `scoring_percent` beside the
    percentile it belongs to.
    """
    from dataclasses import replace as _replace

    from openbiota.taxongroups import species_key

    by_key: dict[str, Any] = {}
    for o in getattr(inventory, "organisms", []) or []:
        by_key[species_key(o.species)] = o
        if getattr(o, "gtdb", None):
            by_key.setdefault(species_key(str(o.gtdb)), o)
        if getattr(o, "formerly", None):
            by_key.setdefault(species_key(str(o.formerly)), o)
        scoring_name = (getattr(o, "native_ids", {}) or {}).get("scoring")
        if scoring_name:
            by_key.setdefault(species_key(str(scoring_name)), o)
        for alias in getattr(o, "aliases", ()) or ():
            by_key.setdefault(species_key(str(alias)), o)

    new_groups: list[GroupReading] = []
    for g in overview.groups:
        members: list[MemberReading] = []
        total = 0.0
        n_exp = 0
        for m in g.members:
            o = by_key.get(species_key(m.species))
            if o is None:
                members.append(m)
                total += m.percent if m.detected else 0.0
                continue
            # The share is the primary lane's composition and nothing else:
            # a member the primary lane did not place carries a share of 0
            # here even when another lane read it, because that reading has
            # a different denominator and may not be added into a group
            # total (spec 0.8.4 §6). It is still marked detected, with the
            # lanes that saw it; its own reading is on its catalogue row.
            share = float(o.percent) if o.in_primary else 0.0
            lanes = tuple(o.lanes)
            exp_only = bool(lanes) and not (set(lanes) & {"scoring", "extended", "genome"}) or (
                bool(lanes) and "scoring" not in lanes and not m.detected)
            detected = bool(lanes) or m.detected
            if exp_only and detected:
                n_exp += 1
            members.append(_replace(m, percent=share if detected else m.percent, detected=detected,
                                    detected_by=lanes, expansion_only=exp_only and detected))
            total += share if detected else 0.0
        new_groups.append(_replace(
            g, members=tuple(members), percent=total, detected=any(m.detected for m in members),
            scoring_percent=g.percent if g.scoring_percent is None else g.scoring_percent,
            n_detected_expansion_only=n_exp,
        ))
    overview.groups[:] = new_groups
    return overview
