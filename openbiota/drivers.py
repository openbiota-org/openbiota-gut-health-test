"""What is behind each metabolic reading: the organisms the matched DNA came from.

Every panel reading is a count of DNA fragments that matched a reference gene.
Each reference gene came from a named organism, and the alignment step keeps
that name beside every fragment. So the question a reader asks first about a
high or low reading - *what in my gut is doing this?* - has a direct,
measured answer that the report used to leave out. A histamine card said
"whether that matters depends on which organisms carry it, which this
measurement cannot tell you." It could. On the sample that sentence was
written for, half the histidine-decarboxylase fragments matched *Eggerthella
lenta*.

This module reads the fragment table for a panel, tallies fragments by the
organism their best-matching reference came from, and joins each organism to
the sample's own inventory: is it there, how much, what class. The result is
rendered as the first thing on the card after the number.

Two honesty limits, stated where they apply:

* A fragment matching *Blautia wexlerae*'s gene says the sequence is closest
  to that reference among the references searched. When the organism is not
  in the inventory the fragment most likely came from a relative whose gene
  is similar; the card says "closest reference" rather than claiming the
  organism is present.
* Direction is a property of the reading's consequences, and for a few
  pathways it depends on who is producing. Where a panel names organisms
  whose production flips the usual reading, the direction shown follows the
  organisms actually found.
"""

from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

#: Fragments below this identity are matching distant relatives and are not
#: used to name a driver. Same floor the panels themselves use for hits.
MIN_IDENTITY: Final = 60.0
#: Drivers shown per card.
TOP_N: Final = 5


@dataclass(frozen=True, slots=True)
class Driver:
    """One organism's contribution to a reading."""

    #: Organism name as the reference carried it, strain qualifiers removed.
    organism: str
    #: Fragments whose best reference came from this organism.
    fragments: int
    #: Share of the panel's target fragments.
    share: float
    #: The organism's row in the sample inventory, when it is there.
    in_sample: bool
    sample_percent: float | None
    sample_class: str | None
    #: True when only the genus was found in the sample.
    genus_in_sample: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "organism": self.organism,
            "fragments": self.fragments,
            "share": round(self.share, 3),
            "in_sample": self.in_sample,
            "sample_percent": self.sample_percent,
            "sample_class": self.sample_class,
            "genus_in_sample": self.genus_in_sample,
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> Driver:
        """The inverse of `to_json`, deriving nothing."""
        return cls(
            organism=str(payload.get("organism") or ""),
            fragments=int(payload.get("fragments") or 0),
            share=float(payload.get("share") or 0.0),
            in_sample=bool(payload.get("in_sample")),
            sample_percent=payload.get("sample_percent"),
            sample_class=payload.get("sample_class"),
            genus_in_sample=bool(payload.get("genus_in_sample")),
        )


@dataclass(frozen=True, slots=True)
class Drivers:
    """The organisms behind one panel's reading."""

    panel: str
    entry_ids: tuple[str, ...]
    total_fragments: int
    top: tuple[Driver, ...]
    #: Direction resolved from the drivers, or ``None`` to use the panel's own.
    resolved_direction: str | None = None
    resolved_reason: str = ""

    @property
    def named_present(self) -> tuple[Driver, ...]:
        return tuple(d for d in self.top if d.in_sample)

    def to_json(self) -> dict[str, Any]:
        return {
            "panel": self.panel,
            "entries": list(self.entry_ids),
            "total_fragments": self.total_fragments,
            "top": [d.to_json() for d in self.top],
            "resolved_direction": self.resolved_direction,
            "resolved_reason": self.resolved_reason,
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> Drivers:
        """Rebuild the record from the file, including its organisms.

        A row read back from `results.json` has to behave like the one that
        was written, and a bare `dict` here does not: the renderer asks the
        record for `.top` and for `named_present`, and gets an attribute
        error a hundred pages into the build.
        """
        return cls(
            panel=str(payload.get("panel") or ""),
            entry_ids=tuple(str(e) for e in (payload.get("entries") or ())),
            total_fragments=int(payload.get("total_fragments") or 0),
            top=tuple(Driver.from_json(d) for d in (payload.get("top") or ())),
            resolved_direction=payload.get("resolved_direction"),
            resolved_reason=str(payload.get("resolved_reason") or ""),
        )


def clean_name(raw: str) -> str:
    """``Eggerthella lenta (strain ATCC 25559 / DSM ...) (Eubacterium lentum)`` -> ``Eggerthella lenta``.

    Reference headers carry strain designations and synonyms in brackets;
    the reader wants the species.
    """
    name = raw.split("(", 1)[0].strip()
    parts = name.split()
    # Drop trailing strain tokens like "7-10-1-b", "DSM", "ATCC", "sp." tails
    # beyond the binomial, but keep a genuine two-word name.
    if len(parts) >= 2:
        genus, epithet = parts[0], parts[1]
        if epithet.lower() in ("sp.", "sp") and len(parts) >= 3:
            return f"{genus} sp. {parts[2]}"
        return f"{genus} {epithet}"
    return name or raw.strip()


def _iter_hits(sample_dir: Path, panel: str) -> Iterable[dict[str, str]]:
    for mate in ("R1", "R2"):
        path = sample_dir / "alignments" / f"hits_{panel}_{mate}.tsv"
        if not path.is_file():
            continue
        with path.open() as fh:
            reader = csv.DictReader((line.lstrip("#") for line in fh), delimiter="\t")
            yield from reader


def for_panel(
    sample_dir: Path,
    panel: str,
    *,
    entry_ids: Iterable[str] = (),
    inventory: Any = None,
    exceptions: Mapping[str, str] | None = None,
    default_direction: str = "unclear",
) -> Drivers | None:
    """Tally the organisms behind one panel's headline reading.

    ``entry_ids`` restricts to the entries that set the headline (the
    panel's aggregate set); empty means every target entry. ``exceptions``
    maps an organism name prefix to the direction its production implies,
    for panels whose reading depends on the producer.
    """
    wanted = set(entry_ids)
    counts: Counter[str] = Counter()
    total = 0
    for row in _iter_hits(sample_dir, panel):
        if row.get("role") != "target":
            continue
        if wanted and row.get("entry") not in wanted:
            continue
        try:
            if float(row.get("pident") or 0.0) < MIN_IDENTITY:
                continue
        except ValueError:
            continue
        counts[clean_name(row.get("organism") or "")] += 1
        total += 1
    if total == 0:
        return None

    by_species: dict[str, Any] = {}
    by_genus: dict[str, list[Any]] = {}
    if inventory is not None:
        for o in inventory.organisms:
            disp = o.display
            by_species[disp] = o
            if o.formerly:
                by_species[o.formerly.replace("_", " ")] = o
            if o.gtdb:
                by_species[o.gtdb] = o
            by_genus.setdefault(disp.split(" ")[0].split("_")[0], []).append(o)

    top: list[Driver] = []
    for name, n in counts.most_common(TOP_N):
        o = by_species.get(name)
        genus = name.split(" ")[0]
        cls = None
        if o is not None:
            from openbiota import organisms as org

            cls = org.verdict(o).cls
        top.append(Driver(
            organism=name,
            fragments=n,
            share=n / total,
            in_sample=o is not None,
            # the organism's one share of the composition; a population
            # counted within a relative's share has none to quote
            sample_percent=(float(o.percent) if (o is not None and o.in_primary) else None),
            sample_class=cls,
            genus_in_sample=(o is None and genus in by_genus),
        ))

    direction, reason = _resolve_direction(counts, total, exceptions or {}, default_direction)
    return Drivers(
        panel=panel, entry_ids=tuple(sorted(wanted)), total_fragments=total,
        top=tuple(top), resolved_direction=direction, resolved_reason=reason,
    )


def _resolve_direction(
    counts: Counter[str], total: int, exceptions: Mapping[str, str], default: str,
) -> tuple[str | None, str]:
    """Follow the producers where the panel says direction depends on them.

    If organisms named as exceptions account for most of the fragments, the
    exception's direction is used and the reason says so. Otherwise the
    default stands. Nothing is resolved for a panel with no exception list.
    """
    if not exceptions:
        return None, ""
    weight: Counter[str] = Counter()
    for name, n in counts.items():
        for prefix, direction in exceptions.items():
            if name.lower().startswith(prefix.lower()):
                weight[direction] += n
                break
    if not weight:
        return default, "none of the producers that would change the reading were found"
    best, n = weight.most_common(1)[0]
    if n / total >= 0.6:
        names = ", ".join(
            sorted({nm for nm in counts if any(
                nm.lower().startswith(p.lower()) and d == best for p, d in exceptions.items())})[:3])
        return best, f"{n / total:.0%} of the matched fragments came from {names}"
    return default, "the producers found do not point one way"


def summarise_for_json(all_drivers: Mapping[str, Drivers | None]) -> dict[str, Any]:
    return {k: (v.to_json() if v is not None else None) for k, v in all_drivers.items()}


__all__ = ["Driver", "Drivers", "MIN_IDENTITY", "TOP_N", "clean_name", "for_panel"]
