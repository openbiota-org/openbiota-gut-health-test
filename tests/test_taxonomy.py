"""Community profile derived from the rpoB fragments."""

from __future__ import annotations

import math

import pytest

from openbiota.taxonomy import BACTEROIDETES, FIRMICUTES, build_profile, canonical_phylum

# --------------------------------------------------------------------------- #
# nomenclature
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Firmicutes", FIRMICUTES),
        ("Bacillota", FIRMICUTES),
        ("bacillota", FIRMICUTES),
        ("Bacteroidetes", BACTEROIDETES),
        ("Bacteroidota", BACTEROIDETES),
        ("Actinobacteria", "Actinomycetota (Actinobacteria)"),
        ("Proteobacteria", "Pseudomonadota (Proteobacteria)"),
        ("", "unassigned"),
        ("Chlamydiota", "Chlamydiota"),
    ],
)
def test_canonical_phylum_merges_old_and_new_names(name: str, expected: str):
    assert canonical_phylum(name) == expected


def test_old_and_new_names_do_not_split_the_ratio():
    """A reference set mixing 'Firmicutes' and 'Bacillota' must not halve the count."""
    profile = build_profile(
        phylum_counts={"Firmicutes": 300, "Bacillota": 300, "Bacteroidota": 200},
        genus_counts={"Roseburia": 600, "Bacteroides": 200},
        organism_counts={"Roseburia intestinalis": 600, "Bacteroides fragilis": 200},
    )
    firmicutes = next(p for p in profile.phyla if p.label == FIRMICUTES)
    assert firmicutes.fragments == 600
    assert profile.fb_ratio == pytest.approx(3.0)


# --------------------------------------------------------------------------- #
# proportions and diversity
# --------------------------------------------------------------------------- #


def test_proportions_sum_to_one():
    profile = build_profile(
        phylum_counts={"Bacillota": 500, "Bacteroidota": 300, "Actinomycetota": 200},
        genus_counts={"A": 500, "B": 300, "C": 200},
        organism_counts={"A sp": 500, "B sp": 300, "C sp": 200},
    )
    assert sum(p.fraction for p in profile.phyla) == pytest.approx(1.0)
    assert sum(g.fraction for g in profile.genera) == pytest.approx(1.0)
    assert profile.total_fragments == 1000


def test_shannon_of_an_even_community_equals_log_richness():
    counts = {f"G{i}": 100 for i in range(8)}
    profile = build_profile(phylum_counts={}, genus_counts=counts, organism_counts=counts)

    assert profile.shannon == pytest.approx(math.log(8))
    assert profile.evenness == pytest.approx(1.0)
    assert profile.simpson == pytest.approx(1 - 8 * (1 / 8) ** 2)


def test_shannon_of_a_monoculture_is_zero():
    profile = build_profile(
        phylum_counts={"Bacillota": 100},
        genus_counts={"OnlyGenus": 100},
        organism_counts={"Only sp": 100},
    )
    assert profile.shannon == pytest.approx(0.0)
    assert profile.evenness is None  # undefined at richness 1
    assert profile.simpson == pytest.approx(0.0)


def test_genus_list_is_truncated_and_sorted():
    counts = {f"G{i:03d}": i + 1 for i in range(50)}
    profile = build_profile(
        phylum_counts={}, genus_counts=counts, organism_counts=counts, genus_limit=5
    )

    assert len(profile.genera) == 5
    assert [g.label for g in profile.genera] == ["G049", "G048", "G047", "G046", "G045"]
    assert profile.genera[0].fragments == 50


def test_empty_profile_is_unavailable():
    profile = build_profile(phylum_counts={}, genus_counts={}, organism_counts={})

    assert profile.available is False
    assert profile.total_fragments == 0
    assert profile.shannon is None
    assert profile.fb_ratio is None
    assert profile.phyla == ()


def test_zero_bacteroidetes_leaves_the_ratio_undefined():
    profile = build_profile(
        phylum_counts={"Bacillota": 500},
        genus_counts={"Roseburia": 500},
        organism_counts={"Roseburia sp": 500},
    )
    assert profile.fb_ratio is None
    assert profile.firmicutes_percent == pytest.approx(100.0)
    assert profile.bacteroidetes_percent == pytest.approx(0.0)


def test_unassigned_lineage_is_labelled_not_dropped():
    profile = build_profile(
        phylum_counts={"": 40, "Bacillota": 60},
        genus_counts={"unassigned": 40, "Roseburia": 60},
        organism_counts={"unknown": 40, "Roseburia sp": 60},
    )
    labels = {p.label for p in profile.phyla}
    assert "unassigned" in labels
    assert sum(p.fragments for p in profile.phyla) == 100


# --------------------------------------------------------------------------- #
# caveats must survive into the JSON
# --------------------------------------------------------------------------- #


def test_json_carries_the_honest_limits():
    profile = build_profile(
        phylum_counts={"Bacillota": 100},
        genus_counts={"Roseburia": 100},
        organism_counts={"Roseburia sp": 100},
    )
    payload = profile.to_json()

    assert payload["available"] is True
    assert "single-copy marker" in payload["basis"]
    joined = " ".join(payload["caveats"])
    assert "not a taxonomic assignment" in joined
    assert "Phylum-level proportions are the reliable level" in joined
    assert "cell fractions" in joined
