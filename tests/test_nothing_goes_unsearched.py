"""Every gene a reading requires must be a gene the panels look for.

Three readings reported "none of this reading's genes was searched for in
this run" - L-lactate formation, succinate formation, and the taurine
route to sulfide. They were not failures of the sample. Nothing in the
reference database could have matched them, because no panel carried
*ldh*, *mdh*, *tpa* or *isl* at all, so the search never asked.

A reading that declares a requirement the assay cannot see is a promise
the report cannot keep, and it reads to somebody as an absence. This
fails the build instead.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from openbiota.extension import fermentation as FERM
from openbiota.panels import load_panel_set

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def searched() -> set[str]:
    """Every gene symbol the panels put into the reference database."""
    panels = load_panel_set(REPO / "panels")
    return {str(t.gene) for p in panels.panels for t in p.targets if t.gene}


def test_every_route_gene_is_in_the_database(searched: set[str]) -> None:
    missing: dict[str, list[str]] = {}
    for route in FERM.ROUTES:
        for gene in route.required_genes or ():
            if str(gene) not in searched:
                missing.setdefault(str(gene), []).append(route.route_id)
    assert not missing, (
        "these routes require genes no panel searches for, so they can only ever "
        f"report as unsearched: {missing}"
    )


def test_the_genes_added_for_the_unsearched_routes_are_present(
    searched: set[str],
) -> None:
    """Named individually, so a later edit that drops one is obvious."""
    for gene in ("ldh", "mdh", "tpa", "isl", "lcdB", "mmdA", "lctA", "lctB", "lctC"):
        assert gene in searched, f"{gene} is no longer searched for"


def test_the_new_targets_are_targets_and_not_decoys() -> None:
    """A decoy subtracts signal. Adding a real gene as one would measure it
    backwards, which is the mistake this nearly made."""
    panels = load_panel_set(REPO / "panels")
    added = {"ldh", "mdh", "tpa", "isl", "lcdB", "mmdA", "lctA", "lctB", "lctC"}
    for panel in panels.panels:
        decoy_genes = {str(d.gene) for d in getattr(panel, "decoys", ()) if d.gene}
        overlap = added & decoy_genes
        assert not overlap, f"{panel.name} carries {overlap} as decoys"


def test_the_new_targets_stay_out_of_the_headline_aggregates() -> None:
    """They were added so the routes could be measured, not to move any
    panel's published figure."""
    import yaml

    expected = {
        "fermentation": {"LDH", "MDH"},
        "propionate": {"LCDB", "MMDA"},
        "butyrate": {"LCTA", "LCTB", "LCTC"},
        "h2s": {"TPA", "ISL"},
    }
    for name, added in expected.items():
        data = yaml.safe_load((REPO / "panels" / f"{name}.yaml").read_text()) or {}
        aggregate = set(data.get("aggregate_from") or ())
        assert not (aggregate & added), (
            f"{name} folded {aggregate & added} into its aggregate, which moves a "
            "figure that was not meant to move"
        )
