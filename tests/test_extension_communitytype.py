"""The gut community-type classifier and its always-available composition view.

The classifier is a frozen artifact. These tests check the arithmetic, the
refusals, and the one property that matters most: nothing here clusters at
report time, and a sample is never compared with the other samples in its
own run.
"""

from __future__ import annotations

import json
import math

import pytest

from openbiota.extension import communitytype as CT

# --------------------------------------------------------------------------- #
# divergence
# --------------------------------------------------------------------------- #


def test_jensen_shannon_is_zero_for_identical_and_ln2_for_disjoint():
    p = [0.2, 0.3, 0.5]
    assert CT.jensen_shannon(p, p) == pytest.approx(0.0, abs=1e-15)
    assert CT.jensen_shannon([1.0, 0.0], [0.0, 1.0]) == pytest.approx(math.log(2))


def test_jensen_shannon_is_symmetric_and_bounded():
    p, q = [0.7, 0.2, 0.1], [0.1, 0.1, 0.8]
    assert CT.jensen_shannon(p, q) == pytest.approx(CT.jensen_shannon(q, p))
    assert 0 <= CT.jensen_shannon(p, q) <= math.log(2) + 1e-12


def test_root_jsd_satisfies_the_triangle_inequality():
    """The reason PAM runs on `sqrt(JSD)` and not on JSD itself."""
    a, b, c = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
    assert CT.root_jsd(a, c) <= CT.root_jsd(a, b) + CT.root_jsd(b, c) + 1e-12
    # JSD without the root does not: 0.693 > 0.693/2 + 0.693/2 is false here,
    # but the classic violation shows on a near-degenerate triple.
    mid = [0.5, 0.5, 0.0]
    assert CT.root_jsd(a, b) <= CT.root_jsd(a, mid) + CT.root_jsd(mid, b) + 1e-12


def test_mismatched_vector_lengths_are_refused():
    with pytest.raises(CT.CommunityTypeError):
        CT.jensen_shannon([0.5, 0.5], [1.0])


# --------------------------------------------------------------------------- #
# composition
# --------------------------------------------------------------------------- #


def test_composition_sums_to_one_over_non_overlapping_partitions():
    comp = CT.composition({
        "Bacteroides_fragilis": 30.0, "Bacteroides_ovatus": 10.0,
        "Prevotella_copri": 20.0, "Blautia_wexlerae": 40.0,
    })
    assert sum(comp.shares.values()) == pytest.approx(1.0)
    assert comp.shares["Bacteroides"] == pytest.approx(0.4)
    assert comp.n_genera == 3 and comp.n_taxa == 4
    assert comp.unresolved_share == 0.0


def test_unresolved_is_a_partition_not_a_dropped_residual():
    comp = CT.composition(
        {"Bacteroides_fragilis": 40.0, "Weirdbacter_novum": 60.0},
        universe=["Bacteroides"],
    )
    assert comp.shares == {"Bacteroides": pytest.approx(0.4), CT.UNRESOLVED: pytest.approx(0.6)}
    assert sum(comp.shares.values()) == pytest.approx(1.0)
    assert "not a residual" in comp.to_json()["partition_note"]


def test_composition_of_an_empty_sample_is_empty_not_zero_shares():
    comp = CT.composition({"Bacteroides_fragilis": 0.0})
    assert comp.shares == {} and comp.total_percent == 0.0


def test_genus_of_strips_gtdb_suffixes_and_brackets():
    assert CT.genus_of("Prevotella_A copri") == "Prevotella"
    assert CT.genus_of("[Ruminococcus]_gnavus") == "Ruminococcus"
    assert CT.genus_of("Blautia_A_wexlerae") == "Blautia"
    assert CT.genus_of("Escherichia coli") == "Escherichia"


# --------------------------------------------------------------------------- #
# the synonymy table
# --------------------------------------------------------------------------- #


def test_renames_reach_the_reference_taxonomy():
    """Without this the sample's largest genus is absent from the reference."""
    assert CT.canonical_genus("Phocaeicola") == "Bacteroides"
    assert CT.canonical_genus("Mediterraneibacter") == "Ruminococcus"
    assert CT.canonical_genus("Segatella") == "Prevotella"
    assert CT.canonical_genus("Lactiplantibacillus") == "Lactobacillus"
    assert CT.canonical_genus("Agathobacter") == "Eubacterium"


def test_a_merely_related_genus_is_not_a_synonym():
    """*Pseudoflavonifractor* is its own genus; folding it in would be a fabrication."""
    assert CT.canonical_genus("Pseudoflavonifractor") == "Pseudoflavonifractor"
    assert CT.canonical_genus("Faecalibacterium") == "Faecalibacterium"
    assert CT.canonical_genus("GGB45596") == "GGB45596"


def test_under_reference_taxonomy_merges_rather_than_overwrites():
    restated = CT.under_reference_taxonomy(
        {"Phocaeicola": 0.2, "Bacteroides": 0.1, "Blautia": 0.7}
    )
    assert restated["Bacteroides"] == pytest.approx(0.3)
    assert sum(restated.values()) == pytest.approx(1.0)


def test_every_synonym_target_is_itself_canonical():
    """A chain (`A -> B -> C`) would make the restatement order-dependent."""
    for source, target in CT.GENUS_SYNONYMS.items():
        if source == target:
            continue
        assert CT.canonical_genus(target) == target, f"{source} -> {target} -> chained"


# --------------------------------------------------------------------------- #
# assignment against a frozen model
# --------------------------------------------------------------------------- #


def _model(stable: bool = True, k: int = 2) -> CT.CommunityTypeModel:
    features = ("Bacteroides", "Prevotella", "Ruminococcus")
    medoids = ((0.8, 0.1, 0.1), (0.1, 0.8, 0.1))[:k]
    return CT.CommunityTypeModel(
        model_id="test/1.0", features=features, medoids=medoids,
        labels=("Bacteroides-dominant", "Prevotella-dominant")[:k],
        distinguishing_taxa=(("Bacteroides",), ("Prevotella",))[:k],
        k=k, tau=0.35, stable=stable, evaluation={}, cohort={},
        reference_coverage={"median": 0.99},
    )


def test_assignment_reports_every_distance_not_only_the_nearest():
    comp = CT.composition({"Bacteroides_fragilis": 80.0, "Prevotella_copri": 10.0,
                           "Ruminococcus_bromii": 10.0})
    out = CT.assign(comp, _model())
    assert out.state == "assigned"
    assert out.nearest == "Bacteroides-dominant"
    assert len(out.distances) == 2
    assert out.distances[0][1] < out.distances[1][1]
    assert sum(s for _, s in out.similarities) == pytest.approx(1.0)


def test_an_equidistant_sample_is_mixed_and_carries_no_name():
    comp = CT.composition({"Bacteroides_fragilis": 45.0, "Prevotella_copri": 45.0,
                           "Ruminococcus_bromii": 10.0})
    out = CT.assign(comp, _model())
    assert out.state == "mixed"
    assert out.nearest is None
    assert "too close to name one" in out.note
    assert len(out.distances) == 2, "the distances are still reported"


def test_an_unstable_model_declines_to_name_a_type():
    comp = CT.composition({"Bacteroides_fragilis": 90.0, "Prevotella_copri": 10.0})
    out = CT.assign(comp, _model(stable=False))
    assert out.state == "continuous_only"
    assert out.nearest is None
    assert "not stable enough" in out.note
    assert out.distances, "composition and distances are still returned"


def test_a_sample_with_nothing_in_the_universe_is_unavailable_not_assigned():
    comp = CT.composition({"Weirdbacter_novum": 100.0})
    out = CT.assign(comp, _model())
    assert out.state == "unavailable" and out.nearest is None and out.distances == ()


def test_assignment_renormalises_so_coverage_does_not_masquerade_as_distance():
    """Two samples with the same *relative* composition must land identically."""
    pure = CT.composition({"Bacteroides_fragilis": 80.0, "Prevotella_copri": 20.0})
    diluted = CT.composition({"Bacteroides_fragilis": 40.0, "Prevotella_copri": 10.0,
                              "Weirdbacter_novum": 50.0})
    a, b = CT.assign(pure, _model()), CT.assign(diluted, _model())
    assert a.distances[0][1] == pytest.approx(b.distances[0][1])
    assert a.coverage == pytest.approx(1.0)
    assert b.coverage == pytest.approx(0.5)


def test_renamed_genera_are_recorded_so_the_reader_can_see_them():
    comp = CT.composition({"Phocaeicola_vulgatus": 70.0, "Prevotella_copri": 30.0})
    out = CT.assign(comp, _model())
    assert out.nearest == "Bacteroides-dominant"
    assert ("Phocaeicola", "Bacteroides", pytest.approx(0.7)) in out.renamed_genera
    payload = out.to_json()
    assert payload["renamed_for_comparison"][0]["as_reported"] == "Phocaeicola"


def test_assignment_never_claims_to_be_a_diagnosis():
    comp = CT.composition({"Bacteroides_fragilis": 90.0, "Prevotella_copri": 10.0})
    text = " ".join(CT.assign(comp, _model()).to_json()["limitations"]).lower()
    assert "not a diagnosis" in text
    assert "high-protein" in text and "absolute" in text


# --------------------------------------------------------------------------- #
# the model file itself
# --------------------------------------------------------------------------- #


def test_a_malformed_model_is_refused_rather_than_half_loaded(tmp_path):
    bad = tmp_path / "m.json"
    bad.write_text(json.dumps({
        "model_id": "x", "features": ["A", "B"], "medoids": [[1.0, 0.0]],
        "labels": ["one", "two"], "distinguishing_taxa": [["A"], ["B"]],
        "k": 2, "tau": 0.3, "stable": True,
    }))
    with pytest.raises(CT.CommunityTypeError, match="k=2 but carries 1 medoid"):
        CT.CommunityTypeModel.load(bad)


def test_a_medoid_that_does_not_span_the_features_is_refused(tmp_path):
    bad = tmp_path / "m.json"
    bad.write_text(json.dumps({
        "model_id": "x", "features": ["A", "B", "C"], "medoids": [[1.0, 0.0]],
        "labels": ["one"], "distinguishing_taxa": [["A"]], "k": 1, "tau": 0.3, "stable": True,
    }))
    with pytest.raises(CT.CommunityTypeError, match="does not span"):
        CT.CommunityTypeModel.load(bad)


def test_a_missing_model_is_an_error_not_an_empty_model(tmp_path):
    with pytest.raises(CT.CommunityTypeError, match="no frozen community-type model"):
        CT.CommunityTypeModel.load(tmp_path / "absent.json")


def test_the_shipped_model_is_loadable_and_internally_consistent():
    model = CT.CommunityTypeModel.load()
    assert model.k >= 2
    assert len(model.medoids) == model.k == len(model.labels)
    for medoid in model.medoids:
        assert len(medoid) == len(model.features)
        assert sum(medoid) == pytest.approx(1.0, abs=1e-9)
    assert model.tau > 0
    evaluation = model.evaluation
    assert set(evaluation["held_out_silhouette"]) == {"2", "3", "4", "5", "6"}
    assert evaluation["selected_k"] == model.k
    assert model.cohort["n_reference_samples"] > 1000
    assert "never ranked against the other samples in its own run" in model.cohort["note"]


def test_the_shipped_model_recovers_the_bacteroides_prevotella_split():
    """Not a fixture: the well-replicated result the reference should reproduce."""
    model = CT.CommunityTypeModel.load()
    labels = " ".join(model.labels).lower()
    assert "bacteroides" in labels and "prevotella" in labels


# --------------------------------------------------------------------------- #
# clustering primitives
# --------------------------------------------------------------------------- #


def test_pam_finds_the_obvious_two_clusters_and_is_deterministic():
    points = [[1.0, 0.0], [0.95, 0.05], [0.9, 0.1], [0.0, 1.0], [0.05, 0.95], [0.1, 0.9]]
    d = [[CT.root_jsd(a, b) for b in points] for a in points]
    first = CT.pam(d, 2, seed=7)
    assert CT.pam(d, 2, seed=7) == first, "the same seed must give the same medoids"
    labels = [min(range(2), key=lambda j: d[i][first[j]]) for i in range(len(points))]
    assert labels[:3] == [labels[0]] * 3
    assert labels[3:] == [labels[3]] * 3
    assert labels[0] != labels[3]


def test_silhouette_is_high_for_separated_clusters_and_low_for_noise():
    points = [[1.0, 0.0], [0.95, 0.05], [0.0, 1.0], [0.05, 0.95]]
    d = [[CT.root_jsd(a, b) for b in points] for a in points]
    assert CT.silhouette(d, [0, 0, 1, 1], 2) > 0.8
    assert CT.silhouette(d, [0, 1, 0, 1], 2) < 0.0


def test_distinguishing_names_the_genus_that_separates_a_cluster():
    vectors = [[0.9, 0.1], [0.85, 0.15], [0.1, 0.9], [0.15, 0.85]]
    features = ["Bacteroides", "Prevotella"]
    assert CT.distinguishing(vectors, [0, 0, 1, 1], features, 0)[0] == "Bacteroides"
    assert CT.distinguishing(vectors, [0, 0, 1, 1], features, 1)[0] == "Prevotella"


def test_genus_matrix_fixes_a_universe_by_prevalence_and_abundance():
    samples = [
        {"Bacteroides_fragilis": 50.0, "Prevotella_copri": 50.0},
        {"Bacteroides_ovatus": 60.0, "Prevotella_copri": 40.0},
        {"Bacteroides_caccae": 70.0, "Rarebacter_only_here": 30.0},
    ]
    features, vectors = CT.genus_matrix(samples, min_prevalence=0.5, min_mean_share=0.01)
    assert features == ["Bacteroides", "Prevotella"]
    assert len(vectors) == 3
    for v in vectors:
        assert sum(v) == pytest.approx(1.0)
    # The one-sample genus is excluded from the universe, not from the sample.
    assert CT.composition(samples[2]).shares["Rarebacter"] == pytest.approx(0.3)


# --------------------------------------------------------------------------- #
# the page-3 contents (A15 §12.2)
# --------------------------------------------------------------------------- #


def test_every_contents_entry_points_at_a_real_section():
    """A renamed or renumbered section must not leave a dangling row."""
    from openbiota.pdfatlas import CONTENTS_ENTRIES
    from openbiota.pdfreport import SECTIONS

    keys = [key for key, _, _ in CONTENTS_ENTRIES]
    assert len(keys) == len(set(keys)), "a destination listed twice"
    for key, title, blurb in CONTENTS_ENTRIES:
        assert key in SECTIONS, f"{key} is not a section"
        assert title and blurb, f"{key} needs a name and a one-line description"
    # Only overview sections: no detail page, no technical appendix.
    assert not any(k.endswith("_detail") for k in keys)
    assert not {"technical", "accuracy", "sequencing", "limits"} & set(keys)
    # And no self-link to the reading guide.
    assert "guide" not in keys


def test_the_contents_lists_the_new_input_section_and_the_renamed_functions():
    from openbiota.pdfatlas import CONTENTS_ENTRIES

    titles = {key: title for key, title, _ in CONTENTS_ENTRIES}
    assert titles["functions"].endswith("Microbial Functions")
    assert "information" in titles["context"].lower()


def test_sections_renumbered_so_the_register_comes_before_the_limits():
    from openbiota.pdfreport import SECTIONS

    # The register sits immediately before the closing four, whatever the
    # absolute numbers (they moved down by one when the contents page
    # became unnumbered).
    c = SECTIONS["context"]
    assert [SECTIONS[k] for k in ("limits", "sequencing", "accuracy", "technical")] == [c + 1, c + 2, c + 3, c + 4]
    assert len(set(SECTIONS.values())) == len(SECTIONS), "two sections share a number"
