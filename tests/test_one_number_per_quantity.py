"""One quantity, one number, everywhere it appears.

The report showed a reading twice under two names and two values. GABA
breakdown read 92nd in the headline list and GABA degradation read 39th in
a table below it; putrescine read 92nd and 38th; glutathione read 81st
beside a zero. The two numbers came from two different statistics - the
panel sums its genes, and an earlier scoring rule took the smallest of
them - presented as though they were the same measurement.

The rule that produced the second number was wrong on its own terms. Every
panel here declares ``aggregate: sum`` and none declares a chain, because
these gene sets are alternatives: gshF replaces gshA plus gshB, putrescine
has three independent routes, cellulose yields to either GH5 or GH9. A
minimum over alternatives answers nothing.

These tests hold the line at the arithmetic rather than at the page, so a
contradiction fails here before anybody reads it in a report.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import readingscore as RS
from openbiota.pdfsynbiotic import owner_of
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06")


def _results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


def _ranges() -> dict[str, Any]:
    from tests._data import ranges

    return ranges()


def _rescored(results: dict[str, Any]) -> dict[str, RS.ReadingScore]:
    """Score every reading again, exactly the way a run scores them.

    Including the step where a reading defers to its panel's curated gene
    choice, or this would test a pipeline nobody runs.
    """
    from openbiota.extension.engine import _defer_to_panel_curation

    stored = (results.get("extension") or {}).get("views", {}).get("reading_scores") or {}
    readings = _defer_to_panel_curation([
        {"reading_id": rid, "genes": [c["gene"] for c in (v.get("contributors") or [])]}
        for rid, v in stored.items()
    ], results)
    return RS.score_all(readings, results=results, ranges=_ranges())


# --- the rule itself -------------------------------------------------------


def test_the_limiting_rule_is_gone() -> None:
    """A minimum over alternative routes is not a measurement of anything."""
    assert RS.RULE == "sum", "readings combine the way their panels declare"
    assert "scarcest" not in RS.RULE_NOTE.lower()


def test_a_reading_is_the_sum_of_the_genes_it_names() -> None:
    """Whatever the position, the value has to be arithmetic anyone can redo."""
    for rid, score in _rescored(_results("SAMPLE2_A02")).items():
        measured = [c for c in score.contributors if c.value is not None]
        if not measured or score.value is None:
            continue
        assert score.value == pytest.approx(sum(c.value for c in measured), rel=1e-4), rid


def test_nearly_every_reading_carries_a_position() -> None:
    """A report full of "not measured" is the failure this replaced.

    Forty readings lost their number when the composite was withdrawn,
    and the cohort had the data to place them the whole time - it was
    reduced to ladders before anybody summed the right columns.
    """
    scores = _rescored(_results("SAMPLE2_A02"))
    placed = [s for s in scores.values() if s.percentile is not None]
    assert len(placed) >= 0.9 * len(scores), (
        f"only {len(placed)} of {len(scores)} readings were placed"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_quantity_is_printed_twice_in_the_overview(sample: str) -> None:
    """The failure a reader can see: one vitamin, two rows, two numbers.

    Vitamin B1 printed 89th beside a thiamine panel printing 76th, and
    Vitamin B5 77th beside pantothenate at 58th. Each pair measured one
    vitamin over two slightly different gene lists, so no comparison of
    gene sets caught them - what catches them is that the owner map
    declares one reading against each of those panels, which makes the
    reading the panel's whole subject.
    """
    from openbiota import pdfextension as PX

    results = _results(sample)
    aggregates = {
        str(p.get("name")): p.get("aggregate_from") or ()
        for p in results.get("panels") or []
    }
    directions = {
        r["panel"]: r.get("higher_means", "unclear")
        for r in results.get("report_rows") or []
    }
    rows = PX.all_functional_rows(results["extension"]["views"], directions, aggregates)
    shown = [r for r in rows if not r.get("duplicates_panel")]

    # Nothing still on the page may name a panel as its sole subject.
    from openbiota.pdfsynbiotic import _sole_reading_of

    sole = _sole_reading_of()
    offenders = [
        r["reading_id"] for r in shown
        if sole.get(str(owner_of(str(r.get("reading_id") or "")) or "")) == r["reading_id"]
    ]
    assert not offenders, (
        "these readings are the whole subject of a panel that is already printed, so "
        f"the quantity appears twice: {offenders}"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_vitamin_appears_once(sample: str) -> None:
    """Named directly, because this is where a reader noticed it."""
    from openbiota import pdfextension as PX

    results = _results(sample)
    rows = PX.all_functional_rows(
        results["extension"]["views"],
        {r["panel"]: r.get("higher_means", "unclear") for r in results["report_rows"]},
        {str(p.get("name")): p.get("aggregate_from") or () for p in results["panels"]},
    )
    shown = {r["reading_id"] for r in rows if not r.get("duplicates_panel")}
    for vitamin in ("b1", "b2", "b3", "b5", "b6", "b7", "b9", "b12", "k2"):
        assert vitamin not in shown, (
            f"{vitamin} is printed as its own row as well as its panel's"
        )


def test_every_reading_distribution_is_built_over_the_whole_cohort() -> None:
    """A ladder from a handful of samples is not a reference."""
    ranges = _ranges()
    readings = ranges.get("readings") or {}
    assert readings, "no reading distributions are installed"
    for rid, spec in readings.items():
        assert spec.get("cohort_n", 0) >= 90, f"{rid} was built over {spec.get('cohort_n')}"
        assert spec.get("rule") == "sum", rid
        assert len(spec.get("percentiles") or {}) >= 5, rid


# --- the invariant the report depends on ----------------------------------


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_reading_of_its_panel_s_quantity_reports_its_panel_s_number(sample: str) -> None:
    """The exact failure: two numbers for one measurement.

    When a reading's genes are precisely the set its panel aggregates, the
    reading and the panel are the same quantity. They may not disagree.
    """
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"]
    panels = {p["name"]: p for p in results["panels"]}
    conflicts: list[str] = []
    for rid, score in _rescored(results).items():
        panel_name = owner_of(rid)
        panel = panels.get(str(panel_name))
        if panel is None:
            continue
        aggregate = set(panel.get("aggregate_from") or ())
        keys = {c.key.split(":", 1)[1] for c in score.contributors if ":" in c.key}
        if not aggregate or keys != aggregate:
            continue  # a narrower reading, which is its own quantity
        reported = (comparison.get(str(panel_name)) or {}).get("percentile")
        if score.percentile is None or reported is None:
            continue
        # A tenth of a percentile is the two ladders rounding their stored
        # values differently, not two answers to one question.
        if abs(score.percentile - reported) > 0.2:
            conflicts.append(
                f"{rid} says {score.percentile} and panel {panel_name} says {reported}"
            )
        # and the underlying value, not only where it lands
        value = (comparison.get(str(panel_name)) or {}).get("value")
        if value is not None and score.value is not None and (
            abs(score.value - value) > max(0.01, abs(value) * 1e-3)
        ):
            conflicts.append(
                f"{rid} measures {score.value} and panel {panel_name} measures {value}"
            )
    assert not conflicts, "one quantity with two numbers: " + "; ".join(conflicts)


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_reading_contradicts_a_measured_absence(sample: str) -> None:
    """A reading may not report a position for something its panel found none of.

    The mirror of the glutathione case, where a panel read 81st and the
    reading beside it read zero.
    """
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"]
    for rid, score in _rescored(results).items():
        panel_name = str(owner_of(rid) or "")
        reported = comparison.get(panel_name) or {}
        if reported.get("percentile") == 0 and score.percentile:
            pytest.fail(
                f"{rid} sits at the {score.percentile}th while its panel "
                f"{panel_name} measured nothing at all"
            )


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_scored_reading_names_the_gene_its_number_came_from(sample: str) -> None:
    """A position with no traceable source cannot be checked by anyone."""
    for rid, score in _rescored(_results(sample)).items():
        if score.percentile is None:
            continue
        assert score.limiting, f"{rid} carries a position and names no gene for it"
        assert any(c.gene == score.limiting for c in score.contributors), rid
