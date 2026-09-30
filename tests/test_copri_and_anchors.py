"""The Prevotella copri complex, and module anchors.

These two features exist for the same reason. "P. copri" is the organism most
often named in the rheumatoid-arthritis microbiome literature, and it was the
least interpretable thing this report could print:

* At species level it is really 13 clades, 13-21% divergent from one another
  (Tett 2019, Manghi 2023), so one abundance number could not say which was
  present. The extended SGB lane already resolves every clade, so the report
  can now name them — or state a measured absence across the whole complex.

* The Scher module's finding *is* P. copri expansion. Without an anchor, that
  module reported the 79th percentile for a sample containing no Prevotella of
  any species, on the strength of one absent Bacteroides species. The accessory
  taxa of a study are not its finding.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from openbiota import copri
from openbiota.profiles import load_profile_set

PROFILES = Path("profiles")


# --------------------------------------------------------------------------- #
# the complex
# --------------------------------------------------------------------------- #


def test_the_resolver_self_test_passes() -> None:
    assert copri.self_test() == 0


def test_every_clade_of_the_complex_is_declared() -> None:
    """13 species-level clades, A through M (Manghi 2023)."""
    letters = [c.clade for c in copri.CLADES]
    assert letters == list("ABCDEFGHIJKLM")
    human = [c for c in copri.CLADES if c.host == "human"]
    assert len(human) == 9, "clades A-I are the human ones; J-M are primate-only"
    # Clade A keeps the name the literature uses.
    assert copri.CLADES[0].species == "Segatella_copri"
    assert "SGB1626" in copri.CLADES[0].sgbs


def test_an_absence_is_reported_as_measured_across_the_whole_complex() -> None:
    out = copri.resolve([{"sgb": "SGB4285", "species": "Ruminococcus_bromii", "percent": 9.1}])
    assert out["status"] == "absent"
    assert out["total_percent"] == 0.0
    assert out["clades_assessed"] == 13
    assert out["human_clades_assessed"] == 9
    # The point of the section: a definite negative, not an unresolved reading.
    assert "measured absence" in out["plain"]


def test_a_clade_is_found_by_sgb_even_after_a_species_rename() -> None:
    """MetaPhlAn renames species between releases; SGB numbers are stable."""
    for label in ("Prevotella_copri", "Segatella_copri", "something_else_entirely"):
        out = copri.resolve([{"sgb": "SGB1626", "species": label, "percent": 3.0}])
        assert out["status"] == "resolved", label
        assert out["clades_present"][0]["clade"] == "A"


def test_the_complex_makes_no_disease_claim() -> None:
    """Manghi tested all 13 clades against 1,635 cases and found nothing."""
    out = copri.resolve([{"sgb": "SGB1626", "species": "Segatella_copri", "percent": 3.0}])
    direction = out["disease_direction"]
    assert "No clade" in direction
    assert "rheumatoid arthritis" in direction
    assert "1,635" in direction
    # No score, no percentile, no direction of any kind in the payload.
    assert "percentile" not in out
    assert "score" not in out


def test_a_missing_extended_lane_is_not_an_absence() -> None:
    out = copri.resolve(None)
    assert out["status"] == "not_assessed"
    assert out["total_percent"] is None


# --------------------------------------------------------------------------- #
# anchors
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def profiles() -> dict[str, Any]:
    return {p.name: p for p in load_profile_set(PROFILES).profiles}


def test_the_ra_profile_anchors_scher_on_the_organism_it_is_about(
    profiles: dict[str, Any],
) -> None:
    module = profiles["ra"].modules["scher_2013"]
    assert module.anchor_features == ("Prevotella_copri",)


def test_an_anchor_must_be_a_feature_of_its_own_module() -> None:
    """Otherwise the gate silently never fires."""
    from openbiota.profiles import parse_profile

    data = {
        "name": "t", "version": "1.0.0", "spec_version": 3, "label": "t",
        "task": "t_vs_healthy", "body_site": "stool",
        "summary": "a test profile", "disease": "t", "family": "test",
        "population": "test", "status": "research_beta", "evidence_maturity": "P3",
        "citation": "none", "effect_size_note": "none",
        "modules": {
            "m": {
                "type": "direction_panel", "evidence_type": "TAX_REL",
                "anchor_features": ["Not_A_Feature"],
                "features": [
                    {"name": "Real_Feature", "level": "species",
                     "direction": "higher_in_case", "w": 1.0},
                ],
            }
        },
    }
    with pytest.raises(Exception, match="anchor_features"):
        parse_profile(data)


def _module_score(*, copri_detected: bool) -> Any:
    """Score the real Scher module against a synthetic sample."""
    from openbiota.scoring import FeatureScore, ModuleScore

    module = {p.name: p for p in load_profile_set(PROFILES).profiles}["ra"].modules["scher_2013"]
    features = []
    for feature in module.features:
        absent = feature.name == "Prevotella_copri" and not copri_detected
        features.append(
            FeatureScore(
                feature=feature,
                measured=True,
                # `raw_value is None` after a successful measurement is the
                # profiler saying it looked and found nothing.
                raw_value=None if absent else 2.0,
                transformed_value=-4.8 if absent else 5.0,
                z=-0.7 if absent else 0.5,
                v=-0.2 if absent else 0.5,
                reference=None,
                percentile=25.0 if absent else 70.0,
            )
        )
    return ModuleScore(
        name=module.name,
        weight_in_combined=1.0,
        score=0.26,
        percentile=78.6,
        features=tuple(features),
        measured_weight=sum(f.feature.factor for f in features),
        planned_weight=sum(f.feature.factor for f in features),
        reference_n=3027,
        reference_source="test",
        module=module,
    )


def test_the_module_withdraws_its_claim_when_the_anchor_is_absent() -> None:
    score = _module_score(copri_detected=False)
    assert score.anchor_absent is True
    assert score.status == "anchor_absent"
    assert score.reportable is False, "no percentile may be shown for a withdrawn claim"
    note = score.anchor_note
    assert "Prevotella copri" in note
    assert "not detected" in note


def test_the_module_scores_normally_when_the_anchor_is_present() -> None:
    score = _module_score(copri_detected=True)
    assert score.anchor_absent is False
    assert score.status == "scored"
    assert score.reportable is True
    assert score.anchor_note == ""


def test_the_measurements_survive_the_withdrawal() -> None:
    """Withdrawing a claim must not delete evidence."""
    score = _module_score(copri_detected=False)
    blob = score.to_json()
    assert blob["anchor_absent"] is True
    assert blob["anchor_features"] == ["Prevotella_copri"]
    assert len(blob["features"]) == 3
    assert all(f["measured"] for f in blob["features"])


def test_a_module_without_anchors_is_unaffected() -> None:
    """The gate is opt-in; no other profile may change behaviour."""
    all_profiles = {p.name: p for p in load_profile_set(PROFILES).profiles}
    anchored = {
        (name, mid)
        for name, prof in all_profiles.items()
        for mid, module in prof.modules.items()
        if module.anchor_features
    }
    assert anchored == {("ra", "scher_2013")}
