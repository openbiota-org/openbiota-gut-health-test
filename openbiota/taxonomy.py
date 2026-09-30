"""Community composition from the rpoB fragments — an extension, not spec.

The rpoB search is already paid for: it is the normalisation denominator. rpoB
is also a classical single-copy phylogenetic marker, so the same fragments
yield a coarse community profile at no extra runtime cost.

Honest limits, restated in every report this feeds:

* The reference set is UniProt-reviewed rpoB only — several hundred sequences.
  Most gut genera have no reviewed rpoB and are therefore attributed to their
  nearest sequenced relative. **Genus labels here are "nearest reference
  organism", not taxonomic assignments.**
* Phylum-level proportions are reliable at rpoB's level of conservation and
  are the numbers worth reading.
* Single-copy marker proportions approximate *cell* fractions, not the
  DNA-mass fractions that a read-count profiler like Kraken or MetaPhlAn
  reports. The two are not interchangeable.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

#: Canonical phylum names, mapping both the historic and the current
#: nomenclature onto one label so the ratios do not silently split.
_PHYLUM_ALIASES: Final = {
    "firmicutes": "Bacillota (Firmicutes)",
    "bacillota": "Bacillota (Firmicutes)",
    "bacteroidetes": "Bacteroidota (Bacteroidetes)",
    "bacteroidota": "Bacteroidota (Bacteroidetes)",
    "actinobacteria": "Actinomycetota (Actinobacteria)",
    "actinomycetota": "Actinomycetota (Actinobacteria)",
    "proteobacteria": "Pseudomonadota (Proteobacteria)",
    "pseudomonadota": "Pseudomonadota (Proteobacteria)",
    "verrucomicrobia": "Verrucomicrobiota",
    "verrucomicrobiota": "Verrucomicrobiota",
    "fusobacteria": "Fusobacteriota",
    "fusobacteriota": "Fusobacteriota",
    "euryarchaeota": "Euryarchaeota",
    "spirochaetes": "Spirochaetota",
    "spirochaetota": "Spirochaetota",
    "synergistetes": "Synergistota",
    "synergistota": "Synergistota",
    "lentisphaerae": "Lentisphaerota",
    "lentisphaerota": "Lentisphaerota",
    "desulfobacterota": "Desulfobacterota",
    "campylobacterota": "Campylobacterota",
}

FIRMICUTES: Final = "Bacillota (Firmicutes)"
BACTEROIDETES: Final = "Bacteroidota (Bacteroidetes)"


def canonical_phylum(name: str) -> str:
    if not name:
        return "unassigned"
    return _PHYLUM_ALIASES.get(name.strip().lower(), name.strip())


@dataclass(frozen=True, slots=True)
class Proportion:
    label: str
    fragments: int
    fraction: float

    def to_json(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "fragments": self.fragments,
            "percent": round(self.fraction * 100.0, 2),
        }


@dataclass(frozen=True, slots=True)
class CommunityProfile:
    """Coarse composition and diversity derived from rpoB best hits."""

    total_fragments: int
    phyla: tuple[Proportion, ...]
    genera: tuple[Proportion, ...]
    distinct_references_hit: int
    shannon: float | None
    simpson: float | None
    evenness: float | None
    firmicutes_percent: float | None
    bacteroidetes_percent: float | None
    fb_ratio: float | None

    @property
    def available(self) -> bool:
        return self.total_fragments > 0

    def to_json(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "basis": "rpoB best-hit attribution (single-copy marker)",
            "total_rpob_fragments": self.total_fragments,
            "distinct_reference_organisms_hit": self.distinct_references_hit,
            "phyla": [p.to_json() for p in self.phyla],
            "nearest_reference_genera": [g.to_json() for g in self.genera],
            "shannon_index": None if self.shannon is None else round(self.shannon, 3),
            "simpson_index": None if self.simpson is None else round(self.simpson, 4),
            "pielou_evenness": None if self.evenness is None else round(self.evenness, 3),
            "firmicutes_percent": (
                None if self.firmicutes_percent is None else round(self.firmicutes_percent, 2)
            ),
            "bacteroidetes_percent": (
                None if self.bacteroidetes_percent is None else round(self.bacteroidetes_percent, 2)
            ),
            "firmicutes_bacteroidetes_ratio": (
                None if self.fb_ratio is None else round(self.fb_ratio, 3)
            ),
            "caveats": [
                "Genus labels are the nearest reviewed rpoB reference organism, not a "
                "taxonomic assignment. Gut genera without a reviewed rpoB are attributed "
                "to a relative.",
                "Phylum-level proportions are the reliable level here; treat genus "
                "proportions as indicative only.",
                "Single-copy-marker proportions approximate cell fractions and are not "
                "comparable to the DNA-mass fractions reported by read-count profilers.",
                "Diversity indices are computed over reference organisms hit, so they "
                "are bounded by the reference set size and understate true diversity.",
            ],
        }


def _proportions(counts: Mapping[str, int], total: int, limit: int | None = None) -> tuple[Proportion, ...]:
    if total <= 0:
        return ()
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if limit is not None:
        items = items[:limit]
    return tuple(Proportion(label=k, fragments=v, fraction=v / total) for k, v in items)


def _shannon(counts: Mapping[str, int]) -> float | None:
    total = sum(counts.values())
    if total <= 0:
        return None
    return -sum((c / total) * math.log(c / total) for c in counts.values() if c > 0)


def _simpson(counts: Mapping[str, int]) -> float | None:
    total = sum(counts.values())
    if total <= 0:
        return None
    return 1.0 - sum((c / total) ** 2 for c in counts.values() if c > 0)


def build_profile(
    *,
    phylum_counts: Mapping[str, int],
    genus_counts: Mapping[str, int],
    organism_counts: Mapping[str, int],
    genus_limit: int = 20,
) -> CommunityProfile:
    canonical: dict[str, int] = {}
    for name, count in phylum_counts.items():
        canonical[canonical_phylum(name)] = canonical.get(canonical_phylum(name), 0) + count

    total = sum(genus_counts.values())
    shannon = _shannon(genus_counts)
    richness = len([g for g, c in genus_counts.items() if c > 0])
    evenness = (
        shannon / math.log(richness) if shannon is not None and richness > 1 else None
    )

    phylum_total = sum(canonical.values())
    firmicutes = canonical.get(FIRMICUTES, 0)
    bacteroidetes = canonical.get(BACTEROIDETES, 0)
    f_pct = (firmicutes / phylum_total * 100.0) if phylum_total else None
    b_pct = (bacteroidetes / phylum_total * 100.0) if phylum_total else None
    fb = (firmicutes / bacteroidetes) if bacteroidetes > 0 else None

    return CommunityProfile(
        total_fragments=total,
        phyla=_proportions(canonical, phylum_total),
        genera=_proportions(genus_counts, total, limit=genus_limit),
        distinct_references_hit=len([o for o, c in organism_counts.items() if c > 0]),
        shannon=shannon,
        simpson=_simpson(genus_counts),
        evenness=evenness,
        firmicutes_percent=f_pct,
        bacteroidetes_percent=b_pct,
        fb_ratio=fb,
    )
