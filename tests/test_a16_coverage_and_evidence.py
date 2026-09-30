"""A16's two coverage scores and the published evidence behind them.

Section 11.4 defines C and P as weighted shares of a goal set. The ways
those go wrong are arithmetic ways: counting an unsolved model as a failed
one, reporting a pair score for something that has no pair, or quietly
dropping the trials that did not work. Each has a test here.
"""

from __future__ import annotations

from typing import Any

import pytest

from openbiota.extension import synbiotic as SYN
from openbiota.extension.engine import _coverage_metrics, _coverage_scores


def _scenarios(*verdicts: str) -> dict[str, Any]:
    return {
        "baseline_medium": "M1",
        "candidates": [
            {
                "medium_id": "M1",
                "substrate": {"substrate_id": f"S{i}", "label": f"fibre {i}"},
                "verdict": v,
            }
            for i, v in enumerate(verdicts, start=1)
        ],
    }


def test_coverage_follows_the_arithmetic_the_specification_gives() -> None:
    """C counts favourable goals, P only pair-supported ones, both over the weight."""
    scores = _coverage_scores(
        _scenarios("pair_gains", "substrate_suffices", "no_gain", "pair_worse")
    )
    # favourable: pair_gains, substrate_suffices -> 2 of 4
    assert scores["goal_coverage"] == pytest.approx(50.0)
    # pair-supported: pair_gains alone -> 1 of 4
    assert scores["pair_supported_goal_coverage"] == pytest.approx(25.0)
    assert scores["n_goals"] == 4


def test_an_unsolved_candidate_is_not_counted_as_an_uncovered_goal() -> None:
    """A model that did not solve is missing evidence, not evidence of absence.

    Counting it against the score would read as "tried, and it did not
    help", which is the opposite of what an unsolved arm means.
    """
    solved_only = _coverage_scores(_scenarios("pair_gains", "no_gain"))
    with_unsolved = _coverage_scores(_scenarios("pair_gains", "no_gain", "unavailable"))
    assert with_unsolved["n_goals"] == solved_only["n_goals"] == 2
    assert with_unsolved["goal_coverage"] == solved_only["goal_coverage"] == 50.0


def test_a_conflicted_result_still_counts_as_a_goal() -> None:
    """`pair_worse` is a finding; dropping it would flatter the score."""
    scores = _coverage_scores(_scenarios("pair_gains", "pair_worse"))
    assert scores["n_goals"] == 2
    assert scores["goal_coverage"] == pytest.approx(50.0)


def test_nothing_eligible_yields_no_score_rather_than_zero() -> None:
    """Zero would claim a goal set was assessed and found wanting."""
    scores = _coverage_scores({"baseline_medium": "M1", "candidates": []})
    assert scores["goal_coverage"] is None
    assert scores["unavailable_reason"]


def test_a_single_component_option_reports_no_pair_score() -> None:
    """P is null for one component, because there is no pair to support anything."""
    scores = SYN.coverage(
        [{"goal": "g", "weight": 1.0, "favourable": True, "pair_supported": False}],
        n_components=1,
    )
    assert scores["goal_coverage"] == pytest.approx(100.0)
    assert scores["pair_supported_goal_coverage"] is None
    assert scores["pair_coverage_unavailable_reason"]


def test_a_negative_weight_is_refused() -> None:
    with pytest.raises(SYN.CoverageError):
        SYN.coverage([{"goal": "g", "weight": -1.0, "favourable": True}], n_components=2)


def test_both_scores_become_scored_metrics_with_a_position() -> None:
    """The report draws a slider from `score_0_100`, so it has to be set."""
    metrics = {m.metric_id: m for m in _coverage_metrics(
        _coverage_scores(_scenarios("pair_gains", "no_gain"))
    )}
    assert set(metrics) == {"A16.goal_coverage", "A16.pair_supported_goal_coverage"}
    for metric in metrics.values():
        assert metric.state == "measured"
        assert metric.score_0_100 == metric.value
        assert 0.0 <= float(metric.score_0_100) <= 100.0
        assert metric.direction == "higher_favourable_in_context"
        assert metric.direction_context


def test_an_unavailable_coverage_carries_its_reason_and_no_number() -> None:
    metrics = _coverage_metrics(_coverage_scores({"baseline_medium": "M1", "candidates": []}))
    for metric in metrics:
        assert metric.state == "not_applicable"
        assert metric.value is None
        assert metric.score_0_100 is None
        assert any("nothing to cover" in lim or "no pair" in lim for lim in metric.limitations)


def test_every_published_record_is_carried_including_the_negatives() -> None:
    """Section 24 has seventeen records and four of them are negative.

    A summary that reported only the encouraging ones would misrepresent
    the literature the whole feature rests on.
    """
    summary = SYN.seed_summary()
    assert summary["n_records"] == 17
    assert summary["n_no_synergy_demonstrated"] == 4
    assert summary["n_missing_a_component_arm"] == 8
    counted = (
        summary["n_demonstrated_synergy"]
        + summary["n_no_synergy_demonstrated"]
        + summary["n_undetermined"]
    )
    assert counted == summary["n_records"], "every record lands in exactly one verdict"


# --- rendering -------------------------------------------------------------


def _render(view: dict[str, Any]) -> str:
    """Draw the two blocks alone and return the text on the page."""
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota import pdfsynbiotic as PS
    from openbiota.pdfreport import _styles

    story: list[Any] = []
    st = _styles()
    PS._coverage_block(story, st, view)
    PS._evidence_seeds_block(story, st, view)
    if not story:
        return ""
    path = "/tmp/_a16_coverage_render.pdf"
    SimpleDocTemplate(path, pagesize=A4).build(story)
    return "".join(page.extract_text() for page in pypdf.PdfReader(path).pages)


def test_coverage_reaches_the_page_with_both_numbers() -> None:
    text = _render({"coverage": {
        "goal_coverage": 83.3, "pair_supported_goal_coverage": 33.3, "n_goals": 6,
    }})
    assert "Goal coverage (C)" in text
    assert "Pair-supported coverage (P)" in text
    assert "83%" in text and "33%" in text
    # the reader must not mistake a modelled share for a clinical one
    assert "not a response rate" in text


def test_a_zero_coverage_is_drawn_at_zero_rather_than_omitted() -> None:
    """Nothing covered is a result. It gets a number and a slider like any other."""
    from openbiota import pdfsynbiotic as PS

    text = _render({"coverage": {
        "goal_coverage": 0.0, "pair_supported_goal_coverage": 0.0, "n_goals": 3,
    }})
    assert "0%" in text
    assert PS._coverage_marker(0.0) is not PS._coverage_marker(100.0)


def test_an_unavailable_pair_score_shows_its_reason_not_a_number() -> None:
    text = _render({"coverage": {
        "goal_coverage": 40.0,
        "pair_supported_goal_coverage": None,
        "pair_coverage_unavailable_reason": "a single-component option has no pair",
        "n_goals": 5,
    }})
    assert "single-component option has no pair" in text
    assert "40%" in text


def test_the_published_evidence_reaches_the_page_with_its_negatives() -> None:
    text = _render({"evidence_seeds": SYN.seed_summary()})
    assert "no synergy demonstrated" in text
    assert "missing a component-only arm" in text
    assert "of 17" in text


def test_the_blocks_stay_silent_when_there_is_nothing_to_say() -> None:
    assert _render({}) == ""
