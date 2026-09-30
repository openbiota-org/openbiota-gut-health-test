"""SGB -> GTDB R232 placement, release-aware.

What was wrong before
---------------------
This module used to read whichever ``*SGB2GTDB*.tsv`` sat in ``refs/gtdb``
- the Jan25/R220 file - and apply it to Jun23 profiles that were then read
beside an R232 Sylph lane. Two mismatches: the bridge was built for a
different MetaPhlAn database (SGB numbers are not stable across releases)
and for a different GTDB release (species names are not either). Spec
0.8.4 §5 calls this out and requires it repaired.

What it does now
----------------
Every lookup names the MetaPhlAn database the SGB came from, and the
answer comes from `openbiota.expansion.crosswalk`: the upstream bridge for
*that* database (Jun23->R207, Jan26->R226) walked onward to R232 through
assembly membership, with the relationship typed (exact, synonym, split,
merge, parent, unresolved). A split renders as a complex, not a species.
The old R220 file is kept on disk only so the misuse can be reproduced in
tests; nothing reads it.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Final

#: GTDB release every placement is expressed in.
RELEASE: Final = "r232"
#: Kept for the historical record: what the repo used to apply.
LEGACY_MISAPPLIED_FILE: Final = Path("refs/gtdb/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB_r220.tsv")

JUN23: Final = "mpa_vJun23_CHOCOPhlAnSGB_202403"
JAN26: Final = "mpa_vJan26_CHOCOPhlAnSGB_202605"
DEFAULT_DB: Final = JUN23  # the installed extended-catalogue lane, until Jan26 rows say otherwise


@dataclass(frozen=True, slots=True)
class GtdbName:
    """One SGB's placement in GTDB R232, with how sure the walk was."""

    sgb: str
    database: str
    relationship: str            # exact | synonym | split | merge | parent | unresolved
    species: str                 # R232 species when exact/synonym/merge; '' otherwise
    genus: str
    bridge_species: str          # what the upstream bridge said, in its own release
    complex_members: tuple[str, ...] = ()

    @property
    def is_placeholder(self) -> bool:
        sp = self.species
        return not sp or (" sp" in sp and sp.split(" sp", 1)[1][:1].isdigit())

    @property
    def resolved(self) -> bool:
        return self.relationship in ("exact", "synonym", "merge") and bool(self.species)

    @property
    def display(self) -> str:
        if self.relationship == "split":
            return f"{self.genus} complex ({', '.join(m.split(' ', 1)[-1] for m in self.complex_members[:3])})"
        if self.species:
            return self.species
        if self.genus:
            return f"{self.genus} sp. ({self.sgb})"
        return self.sgb

    @property
    def family(self) -> str:  # retained for callers; the crosswalk keeps genus and up in bridge_taxonomy
        return ""


@lru_cache(maxsize=4)
def table(database: str = DEFAULT_DB) -> dict[str, GtdbName]:
    """SGB -> placement for one MetaPhlAn database. Empty if the crosswalk is unavailable."""
    try:
        from openbiota.expansion import crosswalk
        edges = crosswalk.load(database)
    except (FileNotFoundError, KeyError, OSError):
        return {}
    out: dict[str, GtdbName] = {}
    for sgb, e in edges.items():
        out[sgb] = GtdbName(
            sgb=sgb, database=database, relationship=e.relationship, species=e.canonical_species,
            genus=e.canonical_genus, bridge_species=e.bridge_species,
            complex_members=tuple(m for m in e.complex_members.split(";") if m) if e.complex_members else (),
        )
    return out


def lookup(sgb: str | None, database: str = DEFAULT_DB) -> GtdbName | None:
    if not sgb:
        return None
    return table(database).get(str(sgb).removeprefix("t__"))


def display_name(sgb: str | None, fallback: str, database: str = DEFAULT_DB) -> str:
    hit = lookup(sgb, database)
    return hit.display if hit is not None else fallback


def available(database: str = DEFAULT_DB) -> bool:
    return bool(table(database))


def source(database: str = DEFAULT_DB) -> str:
    """Where the mapping came from, for the technical record."""
    from openbiota.expansion import crosswalk
    path = crosswalk.CROSSWALK_DIR / f"{database}_to_gtdb_{RELEASE}.tsv.gz"
    return f"{path} (upstream bridge walked to R232 by assembly membership)" if path.is_file() else "unavailable"


__all__ = ["RELEASE", "JUN23", "JAN26", "DEFAULT_DB", "GtdbName", "available", "display_name", "lookup", "source", "table"]
