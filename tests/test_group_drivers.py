"""The organisms behind each group reading.

A group is a sum, and the sum has named parts. Both the glance row and the
detail card name them through one resolver, :func:`pdflibrary.group_drivers`.
These tests hold what that resolver promises: every member is either present
or absent, the present shares add to one, and the abundance quoted for an
organism is the inventory's figure - the same number every other section
prints for it.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

from openbiota import inventory, pdfattention, pdflibrary

RESULTS = sorted(Path("results").glob("*/results.json"))

pytestmark = pytest.mark.skipif(not RESULTS, reason="no results/ to check against")


def _groups(blob: dict) -> list[NS]:
    out = []
    for gj in ((blob.get("profile_similarity") or {}).get("community") or {}).get("groups") or []:
        members = tuple(
            NS(species=m["species"], percent=float(m.get("percent") or 0.0), detected=bool(m.get("detected")),
               percentile=m.get("percentile"), cohort_prevalence=m.get("cohort_prevalence"), in_catalogue=True)
            for m in gj.get("members") or []
        )
        out.append(NS(members=members, percentile=gj.get("percentile"), group=NS(label=gj["label"])))
    return out


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_every_member_is_named_once_and_shares_sum_to_one(path: Path) -> None:
    blob = json.loads(path.read_text())
    inv = inventory.from_json(blob.get("organism_inventory"))
    carded = frozenset(a.species for a in pdfattention.build(inv, None))
    groups = _groups(blob)
    assert groups, "no group readings in results.json"
    for g in groups:
        present, absent = pdflibrary.group_drivers(g, inv=inv, carded=carded)
        assert len(present) + len(absent) == len(g.members), g.group.label
        assert len({d.species for d in present + absent}) == len(g.members), f"{g.group.label}: a member twice"
        if present:
            assert abs(sum(d.share for d in present) - 1.0) < 1e-6, g.group.label
            assert present == sorted(present, key=lambda d: -d.percent), f"{g.group.label}: not largest first"


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_quoted_abundance_is_the_inventory_figure(path: Path) -> None:
    """Section 6 says 6.48% for an organism; the group card must not say 13%."""
    blob = json.loads(path.read_text())
    inv = inventory.from_json(blob.get("organism_inventory"))
    if inv is None:
        pytest.skip("no inventory")
    for g in _groups(blob):
        present, _ = pdflibrary.group_drivers(g, inv=inv, carded=frozenset())
        for d in present:
            o = inv.get(d.species)
            if o is not None and o.best_percent:
                assert d.sample_percent == pytest.approx(o.best_percent), (g.group.label, d.species)


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_links_only_where_a_card_exists(path: Path) -> None:
    blob = json.loads(path.read_text())
    inv = inventory.from_json(blob.get("organism_inventory"))
    carded = frozenset(a.species for a in pdfattention.build(inv, None))
    for g in _groups(blob):
        present, absent = pdflibrary.group_drivers(g, inv=inv, carded=carded)
        for d in present + absent:
            o = inv.get(d.species) if inv is not None else None
            key = o.species if o is not None else d.species
            assert (d.dest is not None) == (key in carded), (g.group.label, d.species)


def test_glance_line_is_readable_and_names_the_missing() -> None:
    """Not grey, not cut: the line names the carriers with their share and
    the common members that were not found."""
    members = (
        NS(species="Faecalibacterium_prausnitzii", percent=0.0, detected=False, percentile=None,
           cohort_prevalence=0.95, in_catalogue=True),
        NS(species="Coprococcus_comes", percent=0.7, detected=True, percentile=90.0, cohort_prevalence=0.8,
           in_catalogue=True),
        NS(species="Odoribacter_splanchnicus", percent=0.7, detected=True, percentile=92.0, cohort_prevalence=0.7,
           in_catalogue=True),
    )
    g = NS(members=members, percentile=4.0, group=NS(label="Butyrate producers"))
    present, absent = pdflibrary.group_drivers(g, inv=None, carded=frozenset())
    line = re.sub(r"<[^>]+>", "", pdflibrary.group_driver_line(present, absent))
    assert "Coprococcus comes 50%" in line and "Odoribacter splanchnicus 50%" in line
    assert "not detected" in line and "Faecalibacterium prausnitzii" in line
    assert "size='6.2'" not in pdflibrary.group_driver_line(present, absent)
