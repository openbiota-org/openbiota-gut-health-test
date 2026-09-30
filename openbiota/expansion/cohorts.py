"""Reference populations for the expanded detection lanes.

A percentile is a position inside a population that was measured the same
way. The scoring lane's population is the 3,027-adult curatedMetagenomicData
cohort in the MetaPhlAn 3 namespace; nothing measured by another lane can
be placed in it. `scripts/build_expansion_cohort.py` builds a population
per lane (`refs/expansion_cohort/<lane>.cohort.json`) by running that
lane's own engine on public adult stool metagenomes. This module loads
them and ranks a lane's readings against them with the same
prevalence-aware CLR percentile the scoring lane uses, so that "82nd
percentile" means one thing throughout the report.

Rules, in force wherever a lane cohort is used:

* One population per organism. An organism ranked on the scoring lane
  keeps that rank; a lane cohort only ranks organisms the scoring
  population never measured.
* The population is named. Every percentile taken here carries
  `percentile_source` ("globdb cohort, n=100"), and the report prints
  which population a rank came from.
* Populations are never merged, averaged or fallen back between.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from openbiota.refcohort import ReferenceCohort

COHORT_DIR = Path("refs/expansion_cohort")

#: inventory lane name -> cohort file stem
LANE_FILES: Mapping[str, str] = {"globdb": "globdb", "jan26": "jan26", "motus": "motus", "kraken": "kraken"}


@dataclass(frozen=True, slots=True)
class LaneCohort:
    lane: str
    cohort: ReferenceCohort
    provenance: Mapping[str, Any]

    @property
    def n(self) -> int:
        return self.cohort.n_samples

    @property
    def source_label(self) -> str:
        return f"{self.lane} cohort, n={self.n}"

    @property
    def description(self) -> str:
        p = self.provenance
        study = p.get("study") or next((w for w in str(p.get("source", "")).split() if w.startswith("PRJ")), "")
        return (f"{self.n} adults from ENA study {study} profiled with {p.get('lane_label', self.lane)} "
                f"at {int(p.get('depth_read_pairs') or 0):,} read pairs each ({p.get('sampling', 'prefix')} sampling); "
                "not matched on country, age or sex").replace("study  ", "study ")


@dataclass(frozen=True, slots=True)
class LaneRank:
    percentile: float | None
    prevalence: float | None
    trace: bool
    rare: bool
    #: median level among the lane population's carriers, in the lane's units
    reference_median: float | None = None
    #: rank among those carriers, on the lane's readings
    carrier_percentile: float | None = None
    #: the lane reading the rank and median were compared with (rescaled to the classified fraction)
    reading: float | None = None
    reference_carriers: int | None = None


_CACHE: dict[str, tuple[float, dict[str, LaneCohort]]] = {}


def load_all(root: Path = COHORT_DIR) -> dict[str, LaneCohort]:
    """Every lane cohort present on disk, keyed by inventory lane name."""
    key = str(root)
    stamp = 0.0
    files: list[tuple[str, Path]] = []
    for lane, stem in LANE_FILES.items():
        p = root / f"{stem}.cohort.json"
        if p.is_file():
            files.append((lane, p))
            stamp = max(stamp, p.stat().st_mtime)
    cached = _CACHE.get(key)
    if cached is not None and cached[0] == stamp and len(cached[1]) == len(files):
        return cached[1]
    out: dict[str, LaneCohort] = {}
    for lane, p in files:
        try:
            raw = json.loads(p.read_text(encoding="utf-8"))
            cohort = ReferenceCohort._load_uncached(p)
        except (OSError, ValueError, KeyError):
            continue
        manifest = dict(raw.get("manifest") or {})
        manifest.setdefault("study", "")
        out[lane] = LaneCohort(lane=lane, cohort=cohort, provenance=manifest)
    _CACHE[key] = (stamp, out)
    return out


def rank(lane_cohort: LaneCohort, species_percent: Mapping[str, float]) -> dict[str, LaneRank]:
    """Rank one lane's species readings against that lane's population.

    Uses the community overview machinery so the arithmetic is the scoring
    lane's: prevalence-aware CLR percentile with 'not carried' as its own
    category. Keys are canonical species names as the lane reports them.
    """
    from openbiota.community import build_community_overview
    from openbiota.taxongroups import TaxonGroupSet

    if not species_percent:
        return {}
    overview = build_community_overview(
        species_percent=dict(species_percent), cohort=lane_cohort.cohort,
        group_set=TaxonGroupSet(groups=()), source_label=lane_cohort.source_label)
    out: dict[str, LaneRank] = {}
    for row in overview.species:
        out[row.species] = LaneRank(
            percentile=row.percentile, prevalence=row.cohort_prevalence,
            trace=bool(getattr(row, "trace", False)), rare=bool(getattr(row, "rare_in_cohort", False)),
            reference_median=getattr(row, "reference_median", None),
            carrier_percentile=getattr(row, "carrier_percentile", None), reading=float(row.percent),
            reference_carriers=getattr(row, "reference_carriers", None))
    return out
