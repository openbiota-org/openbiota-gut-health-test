"""Every functional reading gets a position, or says why it has none.

The readings arrived in the report as adjectives - "class-level only",
"part of it" - in a column beside forty-two rows that each carried a
percentile and a bar. A reader cannot tell high from low in a column of
adjectives, and next to real sliders it reads as missing data.

The distributions were already there: one per panel entry, over the same
91 cohort samples. What these tests pin is the join and the rule for
combining it, because both are places where a plausible-looking average
would quietly mislead.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import readingscore as RS
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parents[1]


LADDER = {"5": 1.0, "10": 2.0, "25": 5.0, "50": 10.0, "75": 20.0, "90": 40.0, "95": 60.0}


def _results(**entries: tuple[str, float | None]) -> dict[str, Any]:
    """A results object holding one panel with the given gene entries."""
    return {"panels": [{
        "name": "p",
        "genes": [
            {"entry_id": eid, "gene": gene, "copies_per_100_genomes": value}
            for eid, (gene, value) in entries.items()
        ],
    }]}


def _ranges(*keys: str) -> dict[str, Any]:
    return {"genes": {k: {"percentiles": dict(LADDER), "detected_in": 80,
                          "cohort_n": 91} for k in keys}}


# --------------------------------------------------------------------------- #
# placing one value
# --------------------------------------------------------------------------- #


def test_a_measured_zero_sits_at_the_bottom_of_the_scale():
    """The same rule as the panel readings: a sample with none of something
    is not a quarter of the way up merely because the reference is mostly
    zeros too."""
    assert RS.percentile_of(LADDER, 0.0) == 0.0


def test_a_value_between_two_quantiles_is_interpolated():
    assert RS.percentile_of(LADDER, 10.0) == pytest.approx(50.0)
    assert RS.percentile_of(LADDER, 15.0) == pytest.approx(62.5)


def test_a_value_beyond_the_ladder_is_clamped_not_extrapolated():
    assert RS.percentile_of(LADDER, 10_000.0) == 95.0
    assert RS.percentile_of(LADDER, 0.5) == 5.0


# --------------------------------------------------------------------------- #
# the combining rule
# --------------------------------------------------------------------------- #


def test_several_genes_are_summed_like_their_panel():
    """The rule, twice corrected.

    It began as the minimum, on the reasoning that a route needs all of
    its parts. These gene sets are not routes in that sense - gshF
    replaces gshA plus gshB, cellulose yields to either GH5 or GH9 - so
    the minimum asked which alternative was rarest and contradicted the
    panel above it, 81st on one line and zero on the next.

    Removing the number was not the answer either: forty readings then
    said "not measured" when the data to measure them existed. The answer
    is the sum, which is what every panel here declares, placed against
    the same sum measured across the cohort.
    """
    results = _results(A=("geneA", 60.0), B=("geneB", 1.0))
    ranges = _ranges("p:A", "p:B")
    ranges["readings"] = {"route": {"percentiles": {
        "5": 1.0, "10": 5.0, "25": 20.0, "50": 40.0, "75": 80.0, "90": 120.0, "95": 200.0,
    }}}
    out = RS.score("route", ["geneA", "geneB"], results=results, ranges=ranges)
    assert out.scored
    assert out.value == pytest.approx(61.0), "the reading is the sum of its genes"
    assert 50 < out.percentile < 75, "placed against the cohort's own sum"
    assert {c.gene for c in out.contributors} == {"geneA", "geneB"}


def test_without_a_distribution_the_parts_are_shown_instead():
    """No cohort sum for this set means no position, and it says so."""
    results = _results(A=("geneA", 60.0), B=("geneB", 1.0))
    out = RS.score("route", ["geneA", "geneB"],
                   results=results, ranges=_ranges("p:A", "p:B"))
    assert not out.scored
    assert out.value == pytest.approx(61.0), "the value is still measured"
    assert "No reference distribution" in (out.unscored_reason or "")


def test_every_contributing_gene_is_listed_with_its_own_position():
    """So the reader can see which gene is the constraint rather than
    taking one number on trust."""
    results = _results(A=("geneA", 60.0), B=("geneB", 1.0))
    out = RS.score("route", ["geneA", "geneB"],
                   results=results, ranges=_ranges("p:A", "p:B"))
    positions = {c.gene: c.percentile for c in out.contributors}
    assert positions == {"geneA": 95.0, "geneB": 5.0}
    assert all(c.cohort_n == 91 for c in out.contributors)


def test_a_single_gene_reading_takes_its_own_position():
    results = _results(A=("geneA", 10.0))
    out = RS.score("step", ["geneA"], results=results, ranges=_ranges("p:A"))
    assert out.percentile == pytest.approx(50.0)
    assert out.limiting == "geneA"


# --------------------------------------------------------------------------- #
# refusing to invent a position
# --------------------------------------------------------------------------- #


def test_a_gene_that_was_never_searched_for_gets_no_position():
    """Not searched for and searched-and-absent are different facts, and
    only one of them is a measurement."""
    results = _results(A=("geneA", 10.0))
    out = RS.score("step", ["geneZ"], results=results,
                   ranges=_ranges("p:A"), searched=["geneA"])
    assert not out.scored
    assert "was searched for" in (out.unscored_reason or "")


def test_a_reading_with_no_cohort_distribution_says_so():
    results = _results(A=("geneA", 10.0))
    out = RS.score("step", ["geneA"], results=results, ranges={"genes": {}})
    assert not out.scored
    assert "no reference distribution" in (out.unscored_reason or "").lower()
    assert out.percentile is None, "an absent distribution must not become a zero"


def test_a_reading_with_no_declared_genes_is_not_scored():
    out = RS.score("step", [], results=_results(), ranges=_ranges())
    assert not out.scored
    assert "no required genes" in (out.unscored_reason or "")


def test_an_undetected_gene_still_places_at_zero_rather_than_vanishing():
    """Zero is a measurement. It belongs on the scale, at the bottom."""
    results = _results(A=("geneA", 0.0))
    out = RS.score("step", ["geneA"], results=results, ranges=_ranges("p:A"))
    assert out.scored
    assert out.percentile == 0.0


def test_score_all_keys_by_reading_and_skips_rows_with_no_id():
    results = _results(A=("geneA", 10.0))
    out = RS.score_all(
        [{"reading_id": "one", "genes": ["geneA"]}, {"genes": ["geneA"]}],
        results=results, ranges=_ranges("p:A"),
    )
    assert set(out) == {"one"}


# --------------------------------------------------------------------------- #
# against the real reports
# --------------------------------------------------------------------------- #

SAMPLES = sorted(
    p.name for p in (REPO / "results").iterdir()
    if p.is_dir() and (p / "results.json").is_file()
) if (REPO / "results").is_dir() else []


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_unplaced_reading_in_a_real_report_says_why(sample: str) -> None:
    """This asked for three quarters of readings to carry a position, and
    that target is what produced the wrong numbers.

    Reaching it meant combining several genes into one figure, and the
    only rule available was a minimum over genes that are alternatives to
    each other. It filled the column and contradicted the panels: 81st on
    one line and zero on the next for the same chemistry. A blank with a
    reason is worth more than a number with none, so the requirement now
    is that every blank explains itself.
    """
    results = json.loads(
        (results_file(sample)).read_text(encoding="utf-8"))
    scores = (results.get("extension") or {}).get("views", {}).get("reading_scores")
    if not scores:
        pytest.skip(f"{sample} predates the reading scores")
    placed = [s for s in scores.values() if s.get("percentile") is not None]
    assert placed, f"{sample}: no reading was placed at all"
    for rid, s in scores.items():
        if s.get("percentile") is None:
            assert s.get("unscored_reason"), f"{rid} is unplaced and does not say why"
        else:
            assert 0.0 <= s["percentile"] <= 100.0, rid
            assert s.get("limiting_gene"), f"{rid} is placed without naming its constraint"
