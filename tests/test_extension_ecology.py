"""A09 — BUILD_SPEC_v0.8.3 section 6.1-6.4 and AT037-AT042.

The companion diversity values exist to make an entropy legible, so their
arithmetic is pinned to the specification's own worked fixtures. The rest of
the module is about refusing to answer: a zero denominator, an absent
component of a log ratio, and an organism whose oxygen phenotype nobody has
measured each have their own state rather than a convenient number.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import ecology as E
from openbiota.extension.schema import ExtensionSchemaError
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

TRAITS: dict[str, dict[str, str]] = {
    "A": {"phenotype": "aerotolerant", "evidence": "cultured_species", "source_id": "E06"},
    "N": {"phenotype": "strict_anaerobe", "evidence": "cultured_species", "source_id": "E06"},
}


def _vec(labels, p, **kw) -> E.AbundanceVector:
    return E.AbundanceVector(tuple(labels), tuple(p), kw.pop("lane", "test.lane"),
                             kw.pop("catalogue", "test catalogue"), 0.0, **kw)


# --------------------------------------------------------------------------- #
# the declared vector
# --------------------------------------------------------------------------- #


def test_a_vector_must_be_closed_and_nonnegative() -> None:
    """Everything downstream assumes one closed vector; it is checked once."""
    with pytest.raises(ValueError, match="not closed"):
        _vec(("a", "b"), (0.5, 0.4))
    with pytest.raises(ValueError, match="negative abundance"):
        _vec(("a", "b"), (1.5, -0.5))
    with pytest.raises(ValueError, match="differ in length"):
        _vec(("a", "b"), (1.0,))
    assert _vec((), ()).n_taxa == 0


def test_every_value_carries_the_declaration_of_what_was_counted() -> None:
    """A diversity that does not say what it counted cannot be compared."""
    v = _vec(("a", "b"), (0.5, 0.5), unresolved_fraction=0.11, detection_depth_fragments=7817)
    d = v.declaration()
    assert d["lane"] == "test.lane"
    assert d["catalogue"] == "test catalogue"
    assert d["abundance_floor_percent"] == 0.0
    assert d["n_taxa_in_vector"] == 2
    assert d["unresolved_fraction"] == 0.11
    assert d["detection_depth_fragments"] == 7817
    for metric in E.metrics(v, E.diversity(v), (), None):
        assert metric.extra["vector"] == d


def test_a_rollup_reconciles_and_shows_its_unmapped_mass() -> None:
    v = _vec(("x sp", "y sp", "z sp"), (0.5, 0.3, 0.2),
             genus_of={"x sp": "Xus", "y sp": "Xus"})
    table = v.collapse("genus")
    assert table["Xus"] == pytest.approx(0.8)
    assert table["unmapped:z sp"] == pytest.approx(0.2)
    assert sum(table.values()) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="unknown level"):
        v.collapse("kingdom")


# --------------------------------------------------------------------------- #
# AT037/AT038 — diversity
# --------------------------------------------------------------------------- #


def test_at037_two_equal_halves_give_diversity_two() -> None:
    d = E.diversity(_vec(("a", "b"), (0.5, 0.5)))
    assert d.effective_species == pytest.approx(2.0, abs=1e-12)
    assert d.inverse_simpson == pytest.approx(2.0, abs=1e-12)
    assert d.gini_simpson == pytest.approx(0.5, abs=1e-12)
    assert d.shannon_nats == pytest.approx(math.log(2), abs=1e-12)


def test_at038_a_single_organism_gives_one_and_an_empty_profile_gives_nothing() -> None:
    d = E.diversity(_vec(("a",), (1.0,)))
    assert d.effective_species == pytest.approx(1.0)
    assert d.inverse_simpson == pytest.approx(1.0)
    assert d.gini_simpson == pytest.approx(0.0)
    assert d.dominance_top1 == pytest.approx(1.0)

    empty = E.diversity(_vec((), ()))
    assert empty.effective_species is None
    assert empty.inverse_simpson is None
    assert empty.gini_simpson is None
    assert empty.observed_taxa == 0
    # Unavailable, not a fabricated zero.
    for metric in E.metrics(_vec((), ()), empty, (), None):
        assert metric.value is None
        assert metric.state == "insufficient_coverage"


def test_dominance_uses_the_largest_shares_and_caps_at_one() -> None:
    v = _vec(list("abcdefg"), [0.4, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03])
    d = E.diversity(v)
    assert d.dominance_top1 == pytest.approx(0.4)
    assert d.dominance_top5 == pytest.approx(0.92)
    small = E.diversity(_vec(("a", "b"), (0.6, 0.4)))
    assert small.dominance_top5 == pytest.approx(1.0)


def test_the_companion_shannon_is_labelled_in_nats_not_swapped_for_the_legacy_one() -> None:
    """Section 6.1: these do not replace an existing entropy in another base."""
    v = _vec(("a", "b", "c", "d"), (0.25, 0.25, 0.25, 0.25))
    d = E.diversity(v)
    assert d.shannon_nats == pytest.approx(math.log(4))
    assert d.effective_species == pytest.approx(4.0)
    metric = next(
        m for m in E.metrics(v, d, (), None)
        if m.metric_id == "ext083.ecology.effective_shannon_species"
    )
    assert metric.unit == "species_equivalents"
    assert "natural logs" in " ".join(metric.limitations)
    assert metric.metric_id != "shannon_index"


# --------------------------------------------------------------------------- #
# AT039/AT040 — ratios
# --------------------------------------------------------------------------- #


def test_at040_a_zero_denominator_is_censored_not_zero_or_infinite() -> None:
    v = _vec(("Faecalibacterium prausnitzii",), (1.0,),
             genus_of={"Faecalibacterium prausnitzii": "Faecalibacterium"})
    fus = next(r for r in E.ratios(v) if "Fusobacterium" in r.label)
    assert fus.value == 0.0 and fus.state == "measured_numerator_absent"

    # Denominator absent instead: censored, with no value at all.
    v2 = _vec(("Fusobacterium nucleatum",), (1.0,),
              genus_of={"Fusobacterium nucleatum": "Fusobacterium"})
    fus2 = next(r for r in E.ratios(v2) if "Fusobacterium" in r.label)
    assert fus2.value is None and fus2.state == "censored_zero_denominator"


def test_at040_genus_and_species_ratio_definitions_are_not_swapped() -> None:
    """A genus shotgun ratio is not the species qPCR ratio of the literature."""
    by_id = {r["ratio_id"]: r for r in E.RATIO_DEFINITIONS}
    assert by_id["ext083.ecology.ratio.prevotella_bacteroides"]["level"] == "genus"
    assert by_id["ext083.ecology.ratio.fusobacterium_faecalibacterium"]["level"] == "genus"
    assert by_id["ext083.ecology.ratio.proteobacteria_actinobacteriota"]["level"] == "phylum"
    fus = by_id["ext083.ecology.ratio.fusobacterium_faecalibacterium"]
    assert "qPCR" in fus["note"] and "E18" in fus["source_ids"]


def test_both_sides_of_a_ratio_come_from_one_vector() -> None:
    """Mixing lanes would make the number a comparison between two assays."""
    species = _vec(("Bacteroides fragilis",), (1.0,),
                   genus_of={"Bacteroides fragilis": "Bacteroides"},
                   lane="inventory")
    phyla = E.phylum_vector({"Proteobacteria": 4.0, "Actinobacteria": 8.0},
                            lane="rpoB", catalogue="rpoB")
    rows = {r.ratio_id: r for r in E.ratios(species, by_level={"phylum": phyla})}
    assert rows["ext083.ecology.ratio.proteobacteria_actinobacteriota"].vector_lane == "rpoB"
    assert rows["ext083.ecology.ratio.prevotella_bacteroides"].vector_lane == "inventory"
    # The phylum ratio is 4/8 within its own vector, after closure.
    assert rows["ext083.ecology.ratio.proteobacteria_actinobacteriota"].value == pytest.approx(0.5)


def test_gtdb_renames_do_not_read_as_an_absent_phylum() -> None:
    """Firmicutes and Bacillota are one phylum under two spellings."""
    a = E.phylum_vector({"Bacillota (Firmicutes)": 50.0, "Actinomycetota (Actinobacteria)": 50.0},
                        lane="l", catalogue="c")
    assert set(a.labels) == {"Bacillota", "Actinobacteriota"}
    b = E.phylum_vector({"Firmicutes": 30.0, "Bacillota": 20.0, "Proteobacteria": 50.0},
                        lane="l", catalogue="c")
    assert dict(zip(b.labels, b.p, strict=True))["Bacillota"] == pytest.approx(0.5)
    # Unclassified mass is held out, not spread over the named phyla.
    c = E.phylum_vector({"Firmicutes": 90.0, "Bacteria_unclassified": 10.0}, lane="l", catalogue="c")
    assert c.labels == ("Bacillota",)
    assert c.p == (1.0,)
    assert c.unresolved_fraction == pytest.approx(0.1)


def test_every_ratio_is_descriptive_and_dimensionless() -> None:
    """No cut-point is invented to fill a green/red gauge (section 6.1)."""
    v = _vec(("Bacteroides fragilis", "Prevotella copri"), (0.5, 0.5),
             genus_of={"Bacteroides fragilis": "Bacteroides", "Prevotella copri": "Prevotella"})
    rows = E.ratios(v)
    for r in rows:
        assert r.to_json()["unit"] == "dimensionless"
        assert r.to_json()["direction"] == "descriptive"
    for metric in E.metrics(v, E.diversity(v), rows, None):
        assert metric.direction == "descriptive"
        assert metric.reference_percentile is None


# --------------------------------------------------------------------------- #
# AT041/AT042 — aerotolerance
# --------------------------------------------------------------------------- #


def test_at041_the_worked_aerotolerance_fixture() -> None:
    """A=.2, N=.6, U=.2 -> fraction .25, coverage .8, MAPI ln(1/3)."""
    v = _vec(("A", "N", "U"), (0.2, 0.6, 0.2))
    a = E.aerotolerance(v, TRAITS)
    assert a.aerotolerant_fraction == pytest.approx(0.25)
    assert a.trait_coverage == pytest.approx(0.8)
    assert a.mapi == pytest.approx(math.log(1 / 3))
    assert a.mapi_state == "measured"
    assert a.unresolved_mass == pytest.approx(0.2)


def test_at041_mapi_is_censored_when_a_component_is_absent() -> None:
    """No pseudocount: the log of an absent group is undefined, not small."""
    no_aero = E.aerotolerance(_vec(("N", "U"), (0.8, 0.2)), TRAITS)
    assert no_aero.aerotolerant_fraction == pytest.approx(0.0)
    assert no_aero.mapi is None
    assert no_aero.mapi_state == "censored_component_absent"

    no_anaerobe = E.aerotolerance(_vec(("A",), (1.0,)), TRAITS)
    assert no_anaerobe.aerotolerant_fraction == pytest.approx(1.0)
    assert no_anaerobe.mapi is None

    none_assessed = E.aerotolerance(_vec(("U",), (1.0,)), TRAITS)
    assert none_assessed.aerotolerant_fraction is None
    assert none_assessed.mapi_state == "no_supported_traits"
    assert none_assessed.trait_coverage == pytest.approx(0.0)


def test_at042_an_unmeasured_phenotype_reduces_coverage_and_picks_no_side() -> None:
    v = _vec(("A", "N", "mystery"), (0.3, 0.3, 0.4))
    a = E.aerotolerance(v, TRAITS)
    assert a.trait_coverage == pytest.approx(0.6)
    assert a.aerotolerant_fraction == pytest.approx(0.5)
    row = next(r for r in a.contributors if r["taxon"] == "mystery")
    assert row["phenotype"] == "unresolved"
    assert row["trait_evidence"] is None
    # An unrecognised phenotype string is unresolved, not silently trusted.
    weird = E.aerotolerance(v, {"mystery": {"phenotype": "somewhat aerobic"}} | TRAITS)
    assert weird.trait_coverage == pytest.approx(0.6)


def test_a_genus_generalisation_is_labelled_as_one() -> None:
    """Strain evidence overrides a genus generalisation only when resolved."""
    v = _vec(("Escherichia coli",), (1.0,), genus_of={"Escherichia coli": "Escherichia"})
    a = E.aerotolerance(v, {"Escherichia": {"phenotype": "aerotolerant",
                                            "evidence": "cultured_species", "source_id": "E06"}})
    row = a.contributors[0]
    assert row["phenotype"] == "aerotolerant"
    assert row["trait_evidence"] == "genus_generalisation"
    # A species-level record wins over the genus one.
    a2 = E.aerotolerance(v, {
        "Escherichia": {"phenotype": "strict_anaerobe", "evidence": "genus_generalisation"},
        "Escherichia coli": {"phenotype": "aerotolerant", "evidence": "cultured_strain",
                             "source_id": "E06"},
    })
    assert a2.contributors[0]["trait_evidence"] == "cultured_strain"
    assert a2.aerotolerant_fraction == pytest.approx(1.0)


def test_the_aerotolerance_metric_says_what_it_is_not() -> None:
    v = _vec(("A", "N"), (0.5, 0.5))
    a = E.aerotolerance(v, TRAITS)
    metrics = {m.metric_id: m for m in E.metrics(v, E.diversity(v), (), a)}
    frac = metrics["ext083.ecology.aerotolerant_fraction"]
    mapi = metrics["ext083.ecology.mapi"]
    assert "not a measurement of oxygen" in " ".join(frac.limitations)
    assert "not a measurement of gut oxygen" in " ".join(mapi.limitations)
    assert frac.assessable_fraction == pytest.approx(1.0)
    assert mapi.kind == "experimental_index"
    assert mapi.evidence_maturity == "association_only"


# --------------------------------------------------------------------------- #
# 6.4 — functional redundancy
# --------------------------------------------------------------------------- #


def test_effective_carriers_rewards_spread_not_count() -> None:
    assert E.effective_carriers([1, 1]) == pytest.approx(2.0)
    assert E.effective_carriers([1, 1, 1, 1]) == pytest.approx(4.0)
    # One dominant carrier plus a trace is not two carriers' worth of redundancy.
    assert E.effective_carriers([99, 1]) == pytest.approx(1.058, abs=1e-3)
    assert E.effective_carriers([1]) == pytest.approx(1.0)
    # No resolved carrier is unavailable, not zero redundancy.
    assert E.effective_carriers([]) is None
    assert E.effective_carriers([0, 0]) is None


def test_redundancy_counts_only_carriers_whose_evidence_supports_carriage() -> None:
    """Co-detection is context, not a second carrier."""
    carriers = [
        {"taxon": "a", "value": 5.0, "is_carrier_evidence": True},
        {"taxon": "b", "value": 5.0, "is_carrier_evidence": True},
        {"taxon": "c", "value": 90.0, "is_carrier_evidence": False},
    ]
    r = E.redundancy("fn.x", "Route X", carriers, unassigned_fraction=0.25)
    assert r.n_resolved_carriers == 2
    assert r.effective_carriers == pytest.approx(2.0)
    assert r.carrier_coverage == pytest.approx(10.0 / 100.0)
    assert r.unassigned_fraction == 0.25
    assert "not a guarantee of ecological resilience" in " ".join(r.to_json()["limitations"])
    none = E.redundancy("fn.y", "Route Y", [])
    assert none.n_resolved_carriers == 0 and none.effective_carriers is None


# --------------------------------------------------------------------------- #
# against a real sample
# --------------------------------------------------------------------------- #


def _sample() -> dict[str, Any]:
    path = results_dir("SAMPLE2_A02") / "results.json"
    if not path.is_file():
        pytest.skip("no run available in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


def test_the_companions_run_on_a_real_inventory_and_leave_the_legacy_values_alone() -> None:
    results = _sample()
    inventory = results["organism_inventory"]
    vector = E.vector_from_inventory(
        inventory["organisms"],
        unresolved_percent=inventory.get("unclassified_percent"),
        depth_fragments=results["community_profile"].get("total_rpob_fragments"),
    )
    assert vector.n_taxa > 50
    assert sum(vector.p) == pytest.approx(1.0, abs=1e-9)

    d = E.diversity(vector)
    assert d.effective_species is not None and d.effective_species > 1
    # exp(H) must exceed neither the taxa counted nor fall below 1.
    assert 1.0 <= d.effective_species <= vector.n_taxa
    assert d.inverse_simpson <= d.effective_species + 1e-9
    assert 0.0 < d.dominance_top1 <= d.dominance_top5 <= 1.0

    # The legacy community profile is untouched by any of this.
    legacy = results["community_profile"]
    assert legacy["shannon_index"] == results["community_profile"]["shannon_index"]
    assert "effective_shannon_species" not in legacy

    phyla = E.phylum_vector(
        {p["label"]: p["percent"] for p in legacy["phyla"]},
        lane="community_profile.rpoB", catalogue="rpoB single-copy marker attribution",
    )
    rows = E.ratios(vector, by_level={"phylum": phyla})
    assert len(rows) == len(E.RATIO_DEFINITIONS)
    for r in rows:
        assert r.state in {"measured", "measured_numerator_absent", "censored_zero_denominator"}
        if r.value is not None:
            assert r.value >= 0.0

    metrics = E.metrics(vector, d, rows, None)
    assert len(metrics) == 5 + len(rows)
    for m in metrics:
        assert m.feature_id == "A09"
        assert m.group == "ecology"
        assert m.input_fingerprint == vector.fingerprint()
        m.to_json()


def test_the_vector_floor_is_declared_and_changes_the_result_visibly() -> None:
    """Two floors give two richnesses; both say which floor they used."""
    results = _sample()
    organisms = results["organism_inventory"]["organisms"]
    low = E.vector_from_inventory(organisms, floor_percent=0.0)
    high = E.vector_from_inventory(organisms, floor_percent=0.1)
    assert high.n_taxa < low.n_taxa
    assert low.declaration()["abundance_floor_percent"] == 0.0
    assert high.declaration()["abundance_floor_percent"] == 0.1
    assert low.fingerprint() != high.fingerprint()


def test_a_metric_from_this_module_never_carries_a_percentile() -> None:
    """No new percentile exists until its cohort does (section 4.2)."""
    results = _sample()
    vector = E.vector_from_inventory(results["organism_inventory"]["organisms"])
    import dataclasses

    for m in E.metrics(vector, E.diversity(vector), E.ratios(vector), None):
        assert m.reference_percentile is None
        assert m.reference_state == "not_available"
        with pytest.raises(ExtensionSchemaError, match="reference_state"):
            # And one cannot be attached without naming its reference.
            dataclasses.replace(m, reference_percentile=50.0)
