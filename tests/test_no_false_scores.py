"""A reading may not wear a score its evidence does not support.

Three contradictions had reached the rendered report and all three said the
same thing in different ways: a position on the scale was drawn for a
reading that had no position.

  * Methane: nought fragments, value 0.0, and a slider sitting at the 28th
    percentile with "NOT DETECTED" printed beside it. Both statements were
    about the same reading and they could not both be true.
  * Colorectal virulence: "above average", 77th percentile, from two
    fragments against a declared floor of five. The sort of sentence that
    frightens somebody for no reason.
  * Histamine: "not detected" on an aggregate of zero while 2,554 fragments
    sat on a target outside the aggregate.

The rules here are deliberately blunt, because the failure was subtle and
the consequence was not. A measured zero is the bottom of the scale. A
count below a panel's own stability floor does not earn a verdict colour.
And anything the report declines to place says so in words.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.pdfreport import status_for

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

SAMPLES = sorted(
    p.name for p in RESULTS.iterdir()
    if p.is_dir() and (p / "results.json").is_file()
) if RESULTS.is_dir() else []


# --------------------------------------------------------------------------- #
# the status vocabulary
# --------------------------------------------------------------------------- #


def test_a_reading_with_nothing_detected_is_not_assessed():
    status = status_for(percentile=27.5, detected=False)
    assert status.label == "not detected"
    assert status.assessed is False
    assert status.unassessed_reason == "not detected"


def test_a_reading_below_its_stability_floor_gets_no_verdict_colour():
    """Two fragments is a rumour, not a measurement."""
    unstable = status_for(percentile=77.5, higher_means="adverse", stable=False)
    assert unstable.assessed is False, (
        "an unstable count must not read as an assessed finding"
    )
    assert "too few reads" in unstable.note
    stable = status_for(percentile=77.5, higher_means="adverse", stable=True)
    assert stable.assessed is True
    assert stable.colour != unstable.colour, (
        "the unstable reading must not carry the same judgement colour as a real one"
    )


def test_the_position_word_is_kept_even_when_the_verdict_is_withheld():
    """Withholding the colour is not the same as hiding where it sits."""
    unstable = status_for(percentile=91.0, stable=False)
    assert unstable.label == status_for(percentile=91.0).label
    assert unstable.level == status_for(percentile=91.0).level


def test_an_unmeasurable_reading_says_so_rather_than_reading_as_typical():
    assert status_for(percentile=None).label == "no result"
    assert status_for(percentile=None).assessed is False


# --------------------------------------------------------------------------- #
# the invariant, over every rendered sample
# --------------------------------------------------------------------------- #


def _results(sample: str) -> dict[str, Any]:
    return json.loads((RESULTS / sample / "results.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_measured_zero_is_reported_at_zero(sample: str) -> None:
    """The headline rule: nothing detected shows zero, and the slider sits
    at zero. A tie-block midpoint is the right convention for a run of
    equal non-zero values and the wrong one here."""
    results = _results(sample)
    offenders = []
    for name, row in (results["reference_comparison"]["panels"] or {}).items():
        value, percentile = row.get("value"), row.get("percentile")
        if value == 0 and percentile not in (0, 0.0, None):
            offenders.append(f"{name}: value {value} but percentile {percentile}")
    assert not offenders, (
        "these readings measured zero and were still given a position on the scale: "
        + "; ".join(offenders)
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_panel_is_called_not_detected_while_showing_a_position(sample: str) -> None:
    results = _results(sample)
    offenders = []
    for name, row in (results["reference_comparison"]["panels"] or {}).items():
        if str(row.get("status")) == "not detected" and row.get("percentile"):
            offenders.append(f"{name}: {row.get('percentile')}th and 'not detected'")
    assert not offenders, "; ".join(offenders)


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_unstable_panel_carries_a_confident_band(sample: str) -> None:
    """A count below the panel's own declared floor must not read as a
    finding, in either direction."""
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"] or {}
    confident = {
        "notably low", "low", "below average", "above average", "high", "notably high",
    }
    offenders = []
    for panel in results["panels"]:
        name = panel["name"]
        row = comparison.get(name) or {}
        if panel.get("stable"):
            continue
        # The label of an unstable reading is still its position - the same
        # vocabulary every other row uses - because the reads are real. What
        # it must not do is present that position as a verdict, and the file
        # says which it is with `assessed` rather than leaving a client to
        # read the note. `histamine` is the case to keep in mind: four
        # thousand fragments in the panel, twelve in the gene the headline
        # is aggregated from, so the count alone looks ample and is not.
        if str(row.get("status")) in confident and row.get("assessed", True):
            offenders.append(
                f"{name}: {panel['accepted_fragments']} fragments, floor "
                f"{panel.get('min_fragments_for_stability')}, reads as "
                f"{row.get('status')!r} and is marked assessed"
            )
    assert not offenders, (
        "these readings are below their own stability floor and still carry a verdict: "
        + "; ".join(offenders)
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_an_unstable_reading_says_so_in_the_file(sample: str) -> None:
    """A client reading the json must not have to parse prose to find out.

    The PDF withholds the colour; a web interface cannot see that, so the
    fact has to be a field.
    """
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"] or {}
    for panel in results["panels"]:
        row = comparison.get(panel["name"]) or {}
        if panel.get("stable") or not row:
            continue
        assert row.get("assessed") is False, (
            f"{panel['name']} is below its stability floor and the file does not say so"
        )
        note = str(row.get("status_note") or "")
        # Three ways a reading can fail to be a verdict, and the stronger one
        # wins: a panel with no reads at all says so, rather than saying it
        # had too few to be sure; and a panel most of the reference cohort
        # lacks says its rank is about rarity, not amount.
        assert (
            "too few reads" in note or "no reads matched" in note
            or "reference cohort carries none" in note
        ), (
            f"{panel['name']} is unstable with no reason recorded (note: {note!r})"
        )


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_zero_aggregate_is_never_described_as_a_quantity(sample: str) -> None:
    """A panel whose aggregate is zero has nothing to be high or low about,
    whatever its other targets found."""
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"] or {}
    for panel in results["panels"]:
        if panel.get("copies_per_100_genomes") not in (0, 0.0):
            continue
        row = comparison.get(panel["name"]) or {}
        assert row.get("percentile") in (0, 0.0, None), (
            f"{panel['name']}: aggregate is zero but the report places it at "
            f"{row.get('percentile')}"
        )


# --------------------------------------------------------------------------- #
# how a zero, and a scored component, are drawn
# --------------------------------------------------------------------------- #


def test_a_measured_zero_is_written_as_zero_and_not_hedged():
    """`<1st` hedges a fact that is not hedged. The value is zero, which is
    the bottom of the scale and not a position near the bottom of it."""
    from openbiota.pdfreport import _ordinal

    assert _ordinal(0) == "0"
    assert _ordinal(0.0) == "0"
    # A genuinely tiny but non-zero position still reads as near the bottom.
    assert _ordinal(0.3) == "<1st"
    assert _ordinal(63) == "63rd"
    assert _ordinal(None) == "\u2014"


def test_a_component_reading_is_coloured_by_the_direction_of_what_it_composes():
    """A grey bar cannot distinguish a good 90th from a bad one, which is
    the whole purpose of the colour. A substrate inherits the direction of
    the carbohydrate reading it helps make."""
    from openbiota import pdfsynbiotic

    rows = [{"reading_id": "carb.cellulose", "group": "Fibre & Dietary Substrates",
             "label": "Cellulose", "found": "GH5", "state": "class_level_only",
             "caveat": ""}]
    out = pdfsynbiotic._attach_scores(
        rows, {"reading_scores": {"carb.cellulose": {"percentile": 90.0,
                                                     "limiting_gene": "celA"}}},
        {"carbohydrates": "favourable"},
    )
    assert out[0]["higher_means"] == "favourable"
    assert out[0]["percentile"] == 90.0
    # And the marker takes the colour a panel row at that position would.
    from openbiota.pdfreport import status_for

    expected = status_for(percentile=90.0, higher_means="favourable").colour
    assert pdfsynbiotic._band_colour(90.0, "favourable") == expected


def test_an_unplaced_component_gets_no_judgement_colour():
    from openbiota import pdfsynbiotic
    from openbiota.pdfreport import CORAL, GREEN

    colour = pdfsynbiotic._band_colour(None, "favourable")
    assert colour not in (GREEN, CORAL), (
        "a reading with no position must not be coloured as though it had one"
    )
