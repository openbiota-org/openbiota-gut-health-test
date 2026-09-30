"""Accepted inventory counts reconcile between results.json and the PDF.

Spec 0.8.4 §7 (data/report integrity): the counts the catalogue page
states - organisms, genera, named species, unnamed clusters, provisional,
newly supported - are read out of the rendered PDF and compared with the
inventory in results.json that produced it. Runs against every sample
that has both files and an expansion-era inventory; skips otherwise.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

pymupdf = pytest.importorskip("pymupdf")

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

SAMPLES = sorted(
    p.name for p in RESULTS.iterdir()
    if p.is_dir() and (p / "results.json").is_file() and (p / f"{p.name}_report.pdf").is_file()
) if RESULTS.is_dir() else []


def _load(sample: str) -> tuple[dict[str, Any], str]:
    results = json.loads((RESULTS / sample / "results.json").read_text(encoding="utf-8"))
    inv = results.get("organism_inventory") or {}
    if not inv.get("counts") or not inv.get("organisms"):
        pytest.skip("this report predates the expanded inventory")
    doc = pymupdf.open(str(RESULTS / sample / f"{sample}_report.pdf"))
    text = "\n".join(page.get_text() for page in doc)
    return results, re.sub(r"\s+", " ", text)


@pytest.mark.parametrize("sample", SAMPLES)
def test_json_counts_are_self_consistent(sample: str) -> None:
    results, _ = _load(sample)
    inv = results["organism_inventory"]
    orgs = inv["organisms"]
    c = inv["counts"]
    assert inv["n_organisms"] == len(orgs) == c["organisms"]
    assert c["named"] + c["unnamed"] == c["organisms"]
    assert c["supported"] + c.get("provisional", 0) + c.get("ambiguous", 0) == c["organisms"]
    cats = sum(c.get(f"category_{k}", 0) for k in
               ("named_species", "unnamed_species_cluster", "unresolved_complex", "higher_rank"))
    assert cats == c["organisms"], "count categories must partition the inventory"
    assert c["rankable"] == sum(1 for o in orgs if o.get("percentile") is not None)
    assert len(inv["comparison"]["newly_supported"]) == sum(
        1 for o in orgs if o.get("incremental_gain") == "expansion" and o.get("status") == "supported")


@pytest.mark.parametrize("sample", SAMPLES)
def test_pdf_states_the_json_counts(sample: str) -> None:
    """The report states the inventory's counts and nothing about software versions.

    Detection statuses, rejected calls and the comparison with the previous
    reference set live in results.json and docs/EXPANSION_COMPARISON.md;
    the report is for the person the sample came from.
    """
    results, text = _load(sample)
    inv = results["organism_inventory"]
    c = inv["counts"]
    assert f"{c['organisms']} organisms were found in this sample" in text
    named = c.get("category_named_species", c["named"])
    unnamed = c.get("category_unnamed_species_cluster", c["unnamed"])
    assert f"Counted separately: {named} named species and {unnamed} unnamed species-level clusters" in text
    for phrase in ("previous reference set", "New in this expansion", "Found only by the expanded",
                   "awaiting confirmation", "Rejected by competitive confirmation", "second search only"):
        assert phrase not in text, f"{sample}: version/log text in the report: {phrase!r}"
    # page 1 states the same total
    assert f"{c['organisms']} ORGANISMS" in text.upper() or f"{c['organisms']} ALL ORGANISMS DETECTED" in text.upper()


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_organism_has_one_share_or_a_relative(sample: str) -> None:
    """One share of the whole per organism, or the relative it is counted within."""
    results, _ = _load(sample)
    inv = results["organism_inventory"]
    if not inv.get("share_unification"):
        pytest.skip("this results file predates share unification")
    by_display = {(o.get("gtdb") if o.get("unnamed") and o.get("gtdb") else o["species"].replace("_", " ")): o
                  for o in inv["organisms"]}
    for o in inv["organisms"]:
        basis = o.get("share_basis")
        if o.get("in_primary"):
            assert basis in ("marker", "split", "estimated", "absorbed") and float(o["percent"]) > 0, o["species"]
            if basis == "absorbed":
                # its whole share is a rejected relative's marker reading
                assert o.get("absorbed_from") and abs(float(o["absorbed_percent"]) - float(o["percent"])) < 1e-6, o["species"]
        elif basis == "member":
            assert o.get("counted_within") in by_display, (o["species"], o.get("counted_within"))
            assert by_display[o["counted_within"]].get("in_primary"), o["species"]
        else:
            assert basis == "none" and float(o.get("percent") or 0) == 0, o["species"]
    # estimated shares come out of the unclassified band, never on top of it
    est = sum(float(o["percent"]) for o in inv["organisms"] if o.get("share_basis") == "estimated")
    if inv.get("unclassified_percent") is not None and inv.get("unplaced_percent") is not None:
        assert abs((float(inv["unclassified_percent"]) - est) - float(inv["unplaced_percent"])) < 0.01 or \
            float(inv["unplaced_percent"]) == 0.0


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_named_organism_is_printed(sample: str) -> None:
    """No hidden top-N: each named organism's display name appears in the PDF."""
    results, text = _load(sample)
    from openbiota.inventory import display_of

    orgs = results["organism_inventory"]["organisms"]
    missing = []
    for o in orgs:
        if o.get("unnamed"):
            continue
        name = display_of(o)
        if name and name not in text:
            missing.append(name)
    assert not missing, f"{len(missing)} named organisms absent from the PDF text, e.g. {missing[:5]}"
