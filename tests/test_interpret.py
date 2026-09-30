"""Per-reading explanations and the cohort comparison layer."""

from __future__ import annotations

from pathlib import Path

import pytest

from openbiota.cohort import build_ranges
from openbiota.interpret import build_rows, cohort_note, reference_json
from openbiota.panels import load_panel_set
from openbiota.pdfreport import _glance_sentence
from openbiota.tally import tally

from .conftest import FIELDS, hit_line

PANELS_DIR = Path(__file__).resolve().parent.parent / "panels"


def _ranges(urda_values: list[float], butyrate_values: list[float]):
    return build_ranges(
        study="TEST",
        description="test cohort",
        reads_per_sample=600_000,
        panel_values={"urda": urda_values, "butyrate": butyrate_values},
        gene_values={},
    )


def _result(database, panel_set, urda: int, butyrate: int, rpob: int, pident: float = 95.0):
    lines = [hit_line(f"U{i}", 0, "URDA", "P00001", pident=pident) for i in range(urda)]
    lines += [hit_line(f"B{i}", 3, "BUT", "P00004", pident=pident) for i in range(butyrate)]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(rpob)]
    return tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )


# --------------------------------------------------------------------------- #
# every shipped panel can explain every reading
# --------------------------------------------------------------------------- #


def test_every_panel_explains_every_direction():
    """A status chip must never appear without an explanation behind it."""
    panel_set = load_panel_set(PANELS_DIR)
    for panel in panel_set.panels:
        interpretation = panel.interpretation
        assert interpretation is not None, panel.name
        for percentile in (95, 80, 50, 20, 5):
            text = interpretation.implication_for(percentile)
            assert text, f"{panel.name} has no explanation at the {percentile}th percentile"
            assert len(text.split()) >= 15, f"{panel.name} @{percentile} is too thin"


def test_explanations_switch_on_the_reading():
    panel_set = load_panel_set(PANELS_DIR)
    for panel in panel_set.panels:
        high = panel.interpretation.implication_for(95)
        typical = panel.interpretation.implication_for(50)
        low = panel.interpretation.implication_for(5)
        assert high != typical != low
        assert len({high, typical, low}) == 3, f"{panel.name} reuses an explanation"


def test_no_explanation_when_there_is_no_percentile():
    panel_set = load_panel_set(PANELS_DIR)
    assert panel_set.by_name("urda").interpretation.implication_for(None) == ""


def test_boundaries_are_at_the_quartiles():
    interpretation = load_panel_set(PANELS_DIR).by_name("urda").interpretation
    assert interpretation.implication_for(75) == interpretation.implication_for(95)
    assert interpretation.implication_for(74) == interpretation.implication_for(50)
    assert interpretation.implication_for(25) == interpretation.implication_for(5)
    assert interpretation.implication_for(26) == interpretation.implication_for(50)


# --------------------------------------------------------------------------- #
# the one-line version used on the glance page
# --------------------------------------------------------------------------- #


def test_glance_sentence_is_one_sentence_and_fits():
    panel_set = load_panel_set(PANELS_DIR)
    ranges = _ranges([10.0, 20.0, 30.0, 40.0, 50.0], [10.0, 20.0, 30.0])

    for panel in panel_set.panels:
        for percentile in (95, 50, 5):
            text = panel.interpretation.implication_for(percentile)
            row = type(
                "R", (), {"implication": text, "detected": True, "percentile": percentile}
            )()
            sentence = _glance_sentence(row)  # type: ignore[arg-type]
            assert sentence, f"{panel.name} @{percentile} produced no glance line"
            # never clipped (spec v04: no ellipsis anywhere); the row wraps instead
            assert "\u2026" not in sentence and "..." not in sentence
            assert "\n" not in sentence
    del ranges


def test_glance_sentence_handles_no_detection():
    row = type("R", (), {"implication": "", "detected": False, "percentile": None})()
    assert "No reads matched" in _glance_sentence(row)  # type: ignore[arg-type]


def test_glance_sentence_handles_missing_cohort():
    row = type("R", (), {"implication": "", "detected": True, "percentile": None})()
    assert "reference group" in _glance_sentence(row)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# row assembly
# --------------------------------------------------------------------------- #


def test_rows_carry_percentile_and_explanation(database, panel_set):
    result = _result(database, panel_set, urda=120, butyrate=400, rpob=600)
    # urda comes out at 40.0 copies per 100 with these counts
    ranges = _ranges([5.0, 10.0, 15.0, 20.0, 25.0], [10.0, 20.0, 30.0])

    rows = build_rows(list(result.panels), ranges)
    urda = next(r for r in rows if r.panel == "urda")

    assert urda.value == pytest.approx(40.0)
    assert urda.percentile is not None and urda.percentile > 75
    assert urda.cohort_median == 15.0
    assert urda.confidence == "confirmed"


def test_rows_survive_a_missing_cohort(database, panel_set):
    result = _result(database, panel_set, urda=120, butyrate=400, rpob=600)
    rows = build_rows(list(result.panels), None)

    urda = next(r for r in rows if r.panel == "urda")
    assert urda.percentile is None
    assert urda.cohort_median is None
    assert urda.status.label in ("no result", "provisional")


def test_low_identity_gene_is_called_out(database, panel_set):
    """A confirmed panel built partly on distant homologs must say so."""
    result = _result(database, panel_set, urda=200, butyrate=400, rpob=600, pident=60.0)
    ranges = _ranges([5.0, 10.0, 15.0], [10.0, 20.0, 30.0])

    rows = build_rows(list(result.panels), ranges)
    urda = next(r for r in rows if r.panel == "urda")

    # 60% identity is below the confirmed floor, so this is provisional
    assert urda.confidence == "provisional"
    assert "provisional" in urda.extra_note.lower()


def test_residue_subset_is_mentioned_when_available(database, panel_set):
    result = _result(database, panel_set, urda=120, butyrate=400, rpob=600)
    rows = build_rows(list(result.panels), _ranges([10.0, 20.0], [10.0, 20.0]))
    urda = next(r for r in rows if r.panel == "urda")

    assert "stricter version" in urda.extra_note


# --------------------------------------------------------------------------- #
# JSON surface
# --------------------------------------------------------------------------- #


def test_reference_json_carries_the_explanation(database, panel_set):
    result = _result(database, panel_set, urda=120, butyrate=400, rpob=600)
    ranges = _ranges([5.0, 10.0, 15.0], [10.0, 20.0, 30.0])
    rows = build_rows(list(result.panels), ranges)

    payload = reference_json(rows, ranges)

    assert payload["cohort"]["study"] == "TEST"
    urda = payload["panels"]["urda"]
    assert urda["what_it_means"]
    assert urda["percentile"] is not None
    assert urda["confidence"] == "confirmed"
    assert urda["higher_means"] == "adverse"


def test_reference_json_without_a_cohort(database, panel_set):
    result = _result(database, panel_set, urda=120, butyrate=400, rpob=600)
    payload = reference_json(build_rows(list(result.panels), None), None)
    assert payload["cohort"] is None


def test_cohort_note_wording():
    assert cohort_note(None) == "a reference group (unavailable)"
    assert "33" in cohort_note(_ranges([1.0] * 33, [1.0] * 33))
