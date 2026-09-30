"""What is behind a reading: the organisms the matched fragments came from."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from openbiota import drivers, inventory

RESULTS = sorted(Path("results").glob("*/results.json"))


def _write_hits(tmp: Path, panel: str, rows: list[tuple[str, str, float, str]]) -> None:
    (tmp / "alignments").mkdir(parents=True, exist_ok=True)
    path = tmp / "alignments" / f"hits_{panel}_R1.tsv"
    with path.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["#qseqid", "sseqid", "pident", "length", "sstart", "send", "bitscore",
                    "qseq_gapped", "sseq_gapped", "role", "entry", "organism"])
        for i, (role, entry, pid, organism) in enumerate(rows):
            w.writerow([f"q{i}", f"s{i}~{entry}~X", pid, 50, 1, 50, 80, "A", "A", role, entry, organism])


def test_strain_qualifiers_are_stripped_from_names() -> None:
    assert drivers.clean_name(
        "Eggerthella lenta (strain ATCC 25559 / DSM 2243) (Eubacterium lentum)") == "Eggerthella lenta"
    assert drivers.clean_name("Gordonibacter pamelaeae 7-10-1-b") == "Gordonibacter pamelaeae"
    assert drivers.clean_name("Blautia sp. CAG:257") == "Blautia sp. CAG:257"


def test_drivers_count_only_target_fragments_above_the_identity_floor(tmp_path: Path) -> None:
    _write_hits(tmp_path, "histamine", [
        ("target", "HDCA", 99.0, "Eggerthella lenta"),
        ("target", "HDCA", 98.0, "Eggerthella lenta"),
        ("target", "HDCA", 97.0, "Jilunia laotingensis"),
        ("target", "HDCA", 40.0, "Klebsiella pneumoniae"),  # below the floor
        ("decoy", "PLPDC3", 99.0, "Escherichia coli"),      # not a target
        ("target", "HDCTRANS", 99.0, "Lactobacillus reuteri"),  # a different entry
    ])
    d = drivers.for_panel(tmp_path, "histamine", entry_ids=["HDCA"])
    assert d is not None and d.total_fragments == 3
    assert d.top[0].organism == "Eggerthella lenta" and d.top[0].fragments == 2
    assert d.top[0].share == pytest.approx(2 / 3)


def test_drivers_are_joined_to_the_sample_inventory(tmp_path: Path) -> None:
    _write_hits(tmp_path, "histamine", [
        ("target", "HDCA", 99.0, "Eggerthella lenta"),
        ("target", "HDCA", 99.0, "Klebsiella pneumoniae"),
    ])
    inv = inventory.Inventory(organisms=[
        inventory.Organism(species="Eggerthella_lenta", percent=0.1, genus="Eggerthella"),
        inventory.Organism(species="Klebsiella_oxytoca", percent=0.01, genus="Klebsiella"),
    ])
    d = drivers.for_panel(tmp_path, "histamine", entry_ids=["HDCA"], inventory=inv)
    by = {x.organism: x for x in d.top}
    assert by["Eggerthella lenta"].in_sample and by["Eggerthella lenta"].sample_percent == pytest.approx(0.1)
    assert by["Eggerthella lenta"].sample_class == "opportunist"
    # The species is absent but its genus is there: said as such, not claimed present.
    assert not by["Klebsiella pneumoniae"].in_sample and by["Klebsiella pneumoniae"].genus_in_sample


def test_direction_follows_the_producers_where_the_panel_says_it_depends(tmp_path: Path) -> None:
    exc = {"Lactobacillus reuteri": "favourable"}
    _write_hits(tmp_path, "histamine", [("target", "HDCA", 99.0, "Lactobacillus reuteri")] * 8
                + [("target", "HDCA", 99.0, "Eggerthella lenta")] * 2)
    d = drivers.for_panel(tmp_path, "histamine", entry_ids=["HDCA"], exceptions=exc,
                          default_direction="adverse")
    assert d.resolved_direction == "favourable", d.resolved_reason
    _write_hits(tmp_path, "h", [("target", "HDCA", 99.0, "Eggerthella lenta")] * 10)
    d2 = drivers.for_panel(tmp_path, "h", entry_ids=["HDCA"], exceptions=exc, default_direction="adverse")
    assert d2.resolved_direction == "adverse"


def test_no_exceptions_means_no_resolution(tmp_path: Path) -> None:
    _write_hits(tmp_path, "x", [("target", "E", 99.0, "Anything at all")])
    d = drivers.for_panel(tmp_path, "x", entry_ids=["E"], default_direction="adverse")
    assert d.resolved_direction is None


def test_no_target_fragments_means_no_drivers(tmp_path: Path) -> None:
    _write_hits(tmp_path, "x", [("decoy", "D", 99.0, "Escherichia coli")])
    assert drivers.for_panel(tmp_path, "x") is None


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_every_detected_panel_in_a_real_report_names_its_drivers() -> None:
    for path in RESULTS:
        r = json.loads(path.read_text())
        drv = r.get("metabolite_drivers")
        if drv is None:
            pytest.skip("results predate drivers")
        # A panel can have fragments on a context entry and none on the gene
        # that sets its headline (Z1's histamine: 2,563 on the antiporter,
        # zero on hdcA). Drivers describe the headline, so the panels that
        # must have them are the ones with fragments on an aggregate entry.
        detected = set()
        for panel in r["panels"]:
            agg = set(panel.get("aggregate_from") or [])
            genes = panel.get("genes") or []
            if any(g.get("entry_id") in agg and g.get("fragments") for g in genes):
                detected.add(panel["name"])
        missing = [n for n in detected if n in drv and drv[n] is None]
        assert not missing, f"{path}: panels with headline fragments but no drivers: {missing}"
        for name, d in drv.items():
            if d:
                assert d["top"], (path, name)
                # Top-N is a subset: shares are each in (0, 1] and never exceed the whole.
                shares = [x["share"] for x in d["top"]]
                assert all(0 < x <= 1 for x in shares), (path, name, shares)
                assert sum(shares) <= 1 + 0.005 * len(shares), (path, name, shares)  # 3-dp rounding
                assert shares == sorted(shares, reverse=True), "largest driver first"
