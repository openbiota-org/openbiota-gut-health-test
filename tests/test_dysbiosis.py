"""The GMWI2 dysbiosis anchor.

The model is the published one, applied verbatim, so the important tests are
that the arithmetic reproduces the reference implementation and that the
module refuses to score profiles from an incompatible database.
"""

from __future__ import annotations

import pytest

from openbiota.dysbiosis import (
    COMPATIBLE_DATABASES,
    MILDLY_NEGATIVE,
    STRONGLY_NEGATIVE,
    compute_anchor,
    load_model,
)


@pytest.fixture(scope="module")
def model():
    return load_model()


def test_model_ships_with_the_package(model) -> None:
    assert model["coefficients"]
    assert len(model["coefficients"]) == 95, "the published model has 95 non-zero coefficients"
    assert model["presence_cutoff"] == pytest.approx(1e-5)


def test_model_has_both_directions(model) -> None:
    coefficients = model["coefficients"].values()
    assert sum(1 for c in coefficients if c > 0) == 49
    assert sum(1 for c in coefficients if c < 0) == 46


def test_model_records_its_provenance(model) -> None:
    assert "Nat Commun" in model["source"]
    assert "8,069" in model["training"]
    assert model["metaphlan_database"].startswith("mpa_v30")


def _clades(model, *, health: int, disease: int) -> dict[str, float]:
    """Build a synthetic profile carrying N health and N disease taxa."""
    positives = [k for k, v in model["coefficients"].items() if v > 0]
    negatives = [k for k, v in model["coefficients"].items() if v < 0]
    out = dict.fromkeys(positives[:health], 1.0)
    out.update(dict.fromkeys(negatives[:disease], 1.0))
    return out


def test_health_associated_taxa_raise_the_score(model) -> None:
    anchor = compute_anchor(
        _clades(model, health=20, disease=0), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    assert anchor.score > 0
    assert anchor.n_health_taxa_present == 20
    assert anchor.n_disease_taxa_present == 0
    assert anchor.band == "not dysbiotic"


def test_disease_associated_taxa_lower_the_score(model) -> None:
    anchor = compute_anchor(
        _clades(model, health=0, disease=20), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    assert anchor.score < STRONGLY_NEGATIVE
    assert anchor.strongly_negative
    assert anchor.band == "strongly dysbiotic"


def test_score_is_the_sum_of_present_coefficients(model) -> None:
    clades = _clades(model, health=5, disease=3)
    expected = sum(model["coefficients"][k] for k in clades) + model["intercept"]
    anchor = compute_anchor(clades, unknown_percent=0.0, database=COMPATIBLE_DATABASES[0])
    assert anchor.score == pytest.approx(expected)


def test_presence_is_thresholded_not_weighted(model) -> None:
    """The model is presence/absence: abundance above the cutoff is irrelevant."""
    positives = [k for k, v in model["coefficients"].items() if v > 0][:4]
    faint = compute_anchor(
        dict.fromkeys(positives, 0.002), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    abundant = compute_anchor(
        dict.fromkeys(positives, 40.0), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    assert faint.score == pytest.approx(abundant.score)


def test_below_cutoff_does_not_count(model) -> None:
    positives = [k for k, v in model["coefficients"].items() if v > 0][:4]
    anchor = compute_anchor(
        dict.fromkeys(positives, 1e-07), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    assert anchor.n_health_taxa_present == 0
    assert anchor.score == pytest.approx(model["intercept"])


def test_unknown_fraction_normalises_the_abundances(model) -> None:
    """Presence is judged against the classified fraction, not the whole sample.

    Abundances arrive as percentages and are divided by ``100 − UNKNOWN``,
    exactly as the reference implementation does, giving a fraction that is
    compared against the 1e-5 cutoff. A taxon at 7e-4 percent therefore falls
    below the cutoff in a fully classified sample (7e-6) and above it in a
    half-classified one (1.4e-5).
    """
    positive = next(k for k, v in model["coefficients"].items() if v > 0)
    fully_classified = compute_anchor(
        {positive: 7e-4}, unknown_percent=0.0, database=COMPATIBLE_DATABASES[0]
    )
    half_classified = compute_anchor(
        {positive: 7e-4}, unknown_percent=50.0, database=COMPATIBLE_DATABASES[0]
    )
    assert fully_classified.n_health_taxa_present == 0
    assert half_classified.n_health_taxa_present == 1


def test_incompatible_database_refuses_to_score() -> None:
    """SGB identifiers from MetaPhlAn 4 are not what the model was trained on."""
    anchor = compute_anchor(
        {"k__Bacteria": 100.0}, unknown_percent=0.0,
        database="mpa_vJan25_CHOCOPhlAnSGB_202503",
    )
    assert not anchor.valid
    assert anchor.score == 0.0
    assert "not comparable" in anchor.note
    assert "could not be computed" in anchor.caution


def test_caution_text_changes_with_the_band(model) -> None:
    strong = compute_anchor(
        _clades(model, health=0, disease=20), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    healthy = compute_anchor(
        _clades(model, health=20, disease=0), unknown_percent=0.0,
        database=COMPATIBLE_DATABASES[0],
    )
    assert "broadly disturbed" in strong.caution
    assert "not being driven by broad gut disturbance" in healthy.caution


def test_empty_profile_is_indeterminate() -> None:
    anchor = compute_anchor({}, unknown_percent=0.0, database=COMPATIBLE_DATABASES[0])
    assert anchor.valid
    assert anchor.band == "indeterminate"


def test_bands_are_ordered() -> None:
    assert STRONGLY_NEGATIVE < MILDLY_NEGATIVE < 0


def test_json_round_trip_is_serialisable(model) -> None:
    import json

    anchor = compute_anchor(
        _clades(model, health=6, disease=4), unknown_percent=12.0,
        database=COMPATIBLE_DATABASES[0],
    )
    payload = json.loads(json.dumps(anchor.to_json()))
    assert payload["index"] == "GMWI2"
    assert payload["valid"] is True
    assert payload["top_contributions"]
    assert "doi" in payload["citation"]
