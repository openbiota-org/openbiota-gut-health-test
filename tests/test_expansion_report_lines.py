"""Report lines the expansion adds to every organism row (spec 0.8.4 §6).

Pure fixtures. Pinned: the identifier line names the lanes' own identifiers
(SGB number, GlobDB accession, mOTU id) once each, lists the methods that
saw the organism, never repeats an identifier, and is empty for an
organism with neither identifiers nor lanes.
"""

from __future__ import annotations

from openbiota import inventory
from openbiota.pdflibrary import _identifiers_line


def _org(**kw):
    base = {"species": "Phocaeicola_vulgatus", "percent": 4.8, "genus": "Phocaeicola"}
    base.update(kw)
    return inventory.Organism(**base)


def test_identifiers_and_lanes_once_each() -> None:
    o = _org(sgb="SGB1814", lanes=("scoring", "jan26", "globdb", "motus"),
             native_ids={"jan26": "SGB1814", "globdb": "globdb:GCF_000012825.1", "motus": "mOTUv4.0_000123"})
    line = _identifiers_line(o)
    assert line.startswith("SGB 1814 · GlobDB GCF_000012825.1 · mOTU 000123")
    assert line.endswith("4 methods")
    assert line.count("1814") == 1


def test_sgb_alone_is_shown_from_the_organism() -> None:
    o = _org(sgb="SGB4837", lanes=("jan26",))
    assert _identifiers_line(o) == "SGB 4837 · 1 method: Jan26"
    two = _org(sgb="SGB4837", lanes=("jan26", "globdb"))
    assert _identifiers_line(two) == "SGB 4837 · 2 methods: Jan26, GlobDB sketch"


def test_empty_for_nothing_to_say() -> None:
    assert _identifiers_line(_org()) == ""
