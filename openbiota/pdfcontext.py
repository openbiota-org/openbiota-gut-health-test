"""The one account of what is present, reachable from any renderer.

Most sections draw their numbers from their own block of ``results.json``.
Two of them - the disease-pattern feature tables and the biofilm cards - are
ranked against a cohort profiled with a single reference catalogue, so their
readings must be that catalogue's or the percentile beside them means
nothing. That is correct and stays, but it leaves a row that says a species
was *not detected* next to an organism list that gives the same species a
share of the community, and a reader is owed the reconciliation on the row
rather than left to find the contradiction.

`build_pdf` puts the merged inventory here once; any renderer can then ask
what the pooled methods found for a name, without threading the inventory
through half a dozen signatures.
"""

from __future__ import annotations

import re
from typing import Any

_INVENTORY: Any | None = None


def set_inventory(inv: Any | None) -> None:
    global _INVENTORY  # noqa: PLW0603 - one process, one report
    _INVENTORY = inv


def inventory() -> Any | None:
    return _INVENTORY


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def organism(name: str) -> Any | None:
    """The inventory's record for a reference-catalogue species name, under any
    spelling that names the *same* organism: the current name, the GTDB name,
    or a merge alias.

    `formerly_listed_as` is deliberately not consulted. It records an older
    catalogue's label for a population that turned out to be a different
    species - a bin the 2023 catalogue called *Ruminococcus gnavus* is
    *Dorea hominis* in the current one - so answering a query for the old name
    with that record would assert the wrong organism is present.
    """
    inv = _INVENTORY
    if inv is None or not name:
        return None
    hit = inv.get(str(name).replace(" ", "_"))
    if hit is not None:
        return hit
    wanted = _norm(name)
    for o in inv.organisms:
        keys = {_norm(o.display), _norm(o.species), _norm(o.gtdb or ""),
                _norm(o.formerly or ""), *(_norm(a) for a in (o.aliases or ()))}
        if wanted in keys:
            return o
    return None


def pooled_share(name: str) -> float | None:
    """The organism's share of the community, when the pooled methods found it
    and gave it one. None when they did not, or when it is counted inside a
    relative's share and so has none of its own."""
    o = organism(name)
    if o is None or not o.in_primary or not o.percent:
        return None
    return float(o.percent)


def pooled_note(name: str, *, catalogue_value: float | None) -> str | None:
    """"0.59% pooled" for a row whose own catalogue read the species as absent
    while every method together found it. None when the two agree."""
    if catalogue_value:
        return None
    share = pooled_share(name)
    if share is None:
        return None
    return f"{share:.3f}% pooled" if share < 1 else f"{share:.2f}% pooled"
