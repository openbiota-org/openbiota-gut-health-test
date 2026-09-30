"""Release-aware crosswalks: the mismatch is fixed and typed (spec 0.8.4 §5, §7).

The repository had been applying the Jan25/R220 bridge to Jun23 profiles
read beside an R232 Sylph lane. These tests pin the repair: every lookup
names its MetaPhlAn database, the walk to R232 goes through assembly
membership, relationships are typed, a split is a complex and not several
species, and nothing reads the old file.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota import gtdb
from openbiota.expansion import crosswalk

REPO = Path(__file__).resolve().parent.parent
HAVE_REFS = (crosswalk.UPSTREAM_DIR / "mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv").is_file() and (
    gtdb.Path("refs/gtdb/r232/genomes.slim.tsv.gz").is_file() or (REPO / "refs/gtdb/r232/bac120_metadata_r232.tsv.gz").is_file()
)
pytestmark = pytest.mark.skipif(not HAVE_REFS, reason="expanded references not fetched in this checkout")


def test_nothing_reads_the_misapplied_r220_file() -> None:
    src = (REPO / "openbiota" / "gtdb.py").read_text()
    assert "SGB2GTDB" not in src.replace("LEGACY_MISAPPLIED_FILE", "").split("def table")[1], (
        "gtdb.table must read the crosswalk, not a *SGB2GTDB*.tsv file"
    )
    assert "crosswalk.load(" in src


def test_every_lookup_names_its_database() -> None:
    """SGB numbers are not stable across releases; a lookup without a database is a bug."""
    a = gtdb.lookup("SGB4421", gtdb.JUN23)
    b = gtdb.lookup("SGB4421", gtdb.JAN26)
    assert a is not None and a.database == gtdb.JUN23
    assert b is None or b.database == gtdb.JAN26
    # The two tables are different objects with different sizes.
    assert len(gtdb.table(gtdb.JUN23)) != len(gtdb.table(gtdb.JAN26))


@pytest.mark.parametrize("db", list(crosswalk.BRIDGES))
def test_relationships_are_typed_and_complete(db: str) -> None:
    edges = crosswalk.load(db)
    summary = crosswalk.summary(db)
    assert set(summary["relationships"]) <= set(crosswalk.RELATIONSHIPS)
    assert sum(summary["relationships"].values()) == len(edges) == summary["n_sgbs"]
    for e in edges.values():
        if e.relationship == "exact":
            assert e.canonical_species == e.bridge_species and e.canonical_species
        elif e.relationship == "synonym":
            assert e.canonical_species and e.canonical_species != e.bridge_species
        elif e.relationship == "split":
            assert not e.canonical_species, "a split must not name one species"
            assert len([m for m in e.complex_members.split(";") if m]) >= 2
        elif e.relationship == "parent":
            assert not e.bridge_species
        elif e.relationship == "merge":
            assert e.canonical_species


def test_a_merge_is_several_sgbs_on_one_species_and_renders_once() -> None:
    edges = crosswalk.load(gtdb.JAN26)
    merges = [e for e in edges.values() if e.relationship == "merge"]
    assert merges, "the Jan26 crosswalk has merge groups"
    by_species: dict[str, list[str]] = {}
    for e in merges:
        by_species.setdefault(e.canonical_species, []).append(e.sgb)
    assert all(len(v) >= 2 for v in by_species.values())
    # Through the inventory, two SGBs of one merge group are one organism.
    from openbiota import inventory

    species, sgbs = next(iter(by_species.items()))
    rows = [{"species": f"X_{i}", "genus": "X", "percent": 1.0, "sgb": s, "sgb_db": gtdb.JAN26, "native_id": s,
             "gtdb": species, "unnamed": False} for i, s in enumerate(sgbs[:2])]
    lane = inventory.Lane(name="jan26", species={r["species"]: 1.0 for r in rows}, rankable=False, rows=rows, primary=True)
    inv = inventory.build([lane])
    names = [o.gtdb for o in inv.organisms]
    assert names.count(species) == 1, f"a merge group rendered {names.count(species)} times"


def test_a_split_is_a_complex_not_species_multiplicity() -> None:
    edges = crosswalk.load(gtdb.JUN23)
    split = next((e for e in edges.values() if e.relationship == "split"), None)
    assert split is not None
    hit = gtdb.lookup(split.sgb, gtdb.JUN23)
    assert hit is not None and hit.relationship == "split" and not hit.resolved
    assert "complex" in hit.display
    from openbiota import inventory

    lane = inventory.Lane(name="extended", species={"Y_sp": 2.0}, rankable=False, primary=True, rows=[
        {"species": "Y_sp", "genus": "Y", "percent": 2.0, "sgb": split.sgb, "sgb_db": gtdb.JUN23, "native_id": split.sgb,
         "gtdb": "", "unnamed": True, "complex": True}])
    inv = inventory.build([lane])
    assert len(inv.organisms) == 1
    assert inv.organisms[0].count_category == "unresolved_complex"
    assert inv.organisms[0].status == "ambiguous"


def test_the_summaries_record_their_inputs_for_invalidation() -> None:
    for db in crosswalk.BRIDGES:
        s = crosswalk.summary(db)
        assert s["input_fingerprint"] and len(s["inputs"]) == 3
        assert s["target_release"] == "r232"


def test_the_historical_misuse_is_measurable() -> None:
    """How many Jun23 names the R220 file would have got wrong: the renames."""
    s = json.loads((crosswalk.CROSSWALK_DIR / "mpa_vJun23_CHOCOPhlAnSGB_202403_to_gtdb_r232.summary.json").read_text())
    assert s["relationships"].get("synonym", 0) > 1000, "R207->R232 renamed thousands of species; the walk must see them"
