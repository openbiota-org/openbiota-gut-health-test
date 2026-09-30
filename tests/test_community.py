"""Curated taxon groups and the community overview built from them."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from openbiota.community import build_community_overview, group_carrier_map
from openbiota.errors import PanelError
from openbiota.refcohort import CohortSample, ReferenceCohort
from openbiota.taxongroups import (
    VALID_GROUP_CATEGORIES,
    TaxonGroupSet,
    load_taxon_group_set,
    parse_taxon_group,
    resolve_against_catalogue,
)

TAXA_DIR = Path(__file__).resolve().parent.parent / "taxa"


def _group(name="g", members=("Alpha one", "Beta two"), **over):
    data = {
        "name": name,
        "label": name.title(),
        "category": "function",
        "higher_means": "favourable",
        "description": "d",
        "summary": "s",
        "citation": "c",
        "members": list(members),
    }
    data.update(over)
    return data


# --------------------------------------------------------------------------- #
# taxon groups: parsing and shipped files
# --------------------------------------------------------------------------- #


def test_member_names_are_normalised_to_profiler_keys():
    g = parse_taxon_group(_group(members=["Roseburia intestinalis", " Faecalibacterium prausnitzii "]))
    assert g.members == ("Roseburia_intestinalis", "Faecalibacterium_prausnitzii")


@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"category": "risk"}, "'category' must be one of"),
        ({"higher_means": "good"}, "'higher_means' must be one of"),
        ({"members": []}, "'members' must not be empty"),
        ({"members": ["A b", "A b"]}, "duplicate"),
        ({"summary": ""}, "'summary' must be a non-empty string"),
    ],
)
def test_invalid_groups_are_refused_with_the_field_named(override, message):
    with pytest.raises(PanelError, match=message):
        parse_taxon_group(_group(**override))


def test_duplicate_group_names_are_refused(tmp_path):
    import yaml

    for fname in ("a.yaml", "b.yaml"):
        (tmp_path / fname).write_text(yaml.safe_dump(_group(name="same")))
    with pytest.raises(PanelError, match="duplicate taxon group names"):
        load_taxon_group_set(tmp_path)


def test_missing_directory_is_an_empty_set(tmp_path):
    assert load_taxon_group_set(tmp_path / "nope").groups == ()


def test_shipped_groups_load_and_cover_the_report_categories():
    group_set = load_taxon_group_set(TAXA_DIR)
    names = {g.name for g in group_set.groups}
    assert {"butyrate", "oral_origin", "methanogens", "sulfate_reducers",
            "mucin_degraders", "opportunistic_pathogens", "hexa_lps"} <= names, names
    assert {g.category for g in group_set.groups} <= VALID_GROUP_CATEGORIES
    # Every group explains itself and cites its source.
    for g in group_set.groups:
        assert len(g.summary) > 80, f"{g.name}: summary too thin"
        assert g.citation
    carriers = group_carrier_map(group_set)
    assert carriers["methanogens"][0] == "Methanobrevibacter_smithii"
    assert carriers == group_set.carrier_map()


def test_resolution_splits_present_from_absent_members():
    gs = TaxonGroupSet(groups=(parse_taxon_group(_group(members=["A b", "C d", "E f"])),))
    present, absent = resolve_against_catalogue(gs, ["A_b", "E_f", "Z_z"])["g"]
    assert present == ("A_b", "E_f")
    assert absent == ("C_d",)


# --------------------------------------------------------------------------- #
# community overview against a synthetic cohort
# --------------------------------------------------------------------------- #


def _cohort(n=60, seed=1) -> ReferenceCohort:
    rng = random.Random(seed)
    taxa = ["Alpha_one", "Beta_two", "Gamma_three", "Delta_four", "Core_a", "Core_b", "Core_c"]
    abundance = []
    for taxon in taxa:
        if taxon == "Delta_four":
            # rare: carried by a tenth of the cohort
            row = [rng.uniform(0.5, 2.0) if rng.random() < 0.1 else 0.0 for _ in range(n)]
        else:
            row = [rng.uniform(1.0, 20.0) for _ in range(n)]
        abundance.append(row)
    ids = [f"s{i}" for i in range(n)]
    meta = {
        sid: CohortSample(sid, "study", "control", "USA", 40.0, "male", 24.0, 10_000_000,
                          "IlluminaHiSeq", None, True)
        for sid in ids
    }
    return ReferenceCohort(snapshot="test", profiler="mpa3", taxa=taxa, sample_ids=ids,
                           abundance=abundance, metadata=meta)


def _groups():
    return TaxonGroupSet(groups=(
        parse_taxon_group(_group(name="ab", members=["Alpha one", "Beta two"])),
        parse_taxon_group(_group(name="rare", members=["Delta four"], category="opportunist",
                                 higher_means="adverse")),
        parse_taxon_group(_group(name="ghost", members=["Not here"], category="origin")),
    ))


def test_overview_places_species_and_groups_against_the_cohort():
    cohort = _cohort()
    sample = {"Alpha_one": 30.0, "Beta_two": 25.0, "Gamma_three": 1.0, "Core_a": 10.0,
              "Core_b": 10.0, "Core_c": 10.0, "Unknown_sp": 14.0}
    overview = build_community_overview(species_percent=sample, cohort=cohort, group_set=_groups())

    assert overview.n_species_detected == 7
    assert overview.n_not_in_catalogue == 1
    rows = {r.species: r for r in overview.species}
    assert rows["Alpha_one"].percentile is not None and rows["Alpha_one"].percentile > 75
    assert rows["Gamma_three"].percentile is not None and rows["Gamma_three"].percentile < 25
    assert rows["Unknown_sp"].percentile is None and not rows["Unknown_sp"].in_catalogue
    assert rows["Alpha_one"].groups == ("ab",)

    ab = overview.group("ab")
    assert ab is not None and ab.detected
    assert ab.percent == pytest.approx(55.0 * (100.0 / 100.0), rel=0.05)
    assert ab.percentile is not None and ab.percentile > 75
    assert ab.n_resolvable == 2
    assert {m.species for m in ab.members} == {"Alpha_one", "Beta_two"}

    rare = overview.group("rare")
    assert rare is not None and not rare.detected
    # Absent, like ~90% of the cohort: placed mid-rank among the absent, not "low".
    assert rare.percentile is not None and 30 < rare.percentile < 60
    assert rare.cohort_prevalence == pytest.approx(0.1, abs=0.06)

    ghost = overview.group("ghost")
    assert ghost is not None and ghost.percentile is None and ghost.n_resolvable == 0
    assert any("ghost" in n for n in overview.notes)

    payload = overview.to_json()
    assert payload["n_species_detected"] == 7
    assert {g["group"] for g in payload["groups"]} == {"ab", "rare", "ghost"}


def test_overview_without_a_cohort_reports_abundance_only():
    overview = build_community_overview(
        species_percent={"Alpha_one": 5.0, "Beta_two": 0.0}, cohort=None, group_set=_groups()
    )
    assert overview.n_species_detected == 1
    assert all(r.percentile is None for r in overview.species)
    ab = overview.group("ab")
    assert ab is not None and ab.detected and ab.percentile is None
    assert any("No reference cohort" in n for n in overview.notes)


def test_by_category_groups_readings():
    overview = build_community_overview(
        species_percent={"Alpha_one": 5.0}, cohort=None, group_set=_groups()
    )
    cats = overview.by_category()
    assert {g.group.name for g in cats["function"]} == {"ab"}
    assert {g.group.name for g in cats["opportunist"]} == {"rare"}


def test_level_percentile_and_deviation_always_point_the_same_way() -> None:
    """A rare organism carried at a low level is not 'high' just because most
    people lack it: the level percentile ranks among carriers, on the same
    distribution the deviation is measured against, so the two agree in sign."""
    cohort = _cohort()
    # Delta_four is carried by a tenth of the cohort at 0.5-2.0; this sample has a little of it
    sample = {"Alpha_one": 40.0, "Beta_two": 20.0, "Gamma_three": 15.0, "Delta_four": 0.6, "Core_a": 8.0,
              "Core_b": 8.0, "Core_c": 8.4}
    overview = build_community_overview(species_percent=sample, cohort=cohort, group_set=_groups())
    rows = {r.species: r for r in overview.species}
    delta = rows["Delta_four"]
    assert delta.percentile is not None and delta.percentile > 90, "prevalence-aware: above the 90% who lack it"
    assert delta.carrier_percentile is not None and delta.carrier_percentile < 50, "among carriers it is low"
    assert delta.deviation_percent is not None and delta.deviation_percent < 0
    for r in overview.species:
        if r.carrier_percentile is None or r.deviation_percent is None:
            continue
        if r.carrier_percentile > 50:
            assert r.deviation_percent >= 0, r.species
        if r.carrier_percentile < 50:
            assert r.deviation_percent <= 0, r.species
    for g in overview.groups:
        if g.detected and g.carrier_percentile is not None and g.deviation_percent is not None:
            if g.carrier_percentile > 50:
                assert g.deviation_percent >= 0, g.group.name
            if g.carrier_percentile < 50:
                assert g.deviation_percent <= 0, g.group.name
