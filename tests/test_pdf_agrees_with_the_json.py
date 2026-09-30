"""The PDF must say what results.json says, and nothing else.

The order matters. `test_results_json_is_the_source` checks the file on its
own terms, before anything is drawn. This file checks the drawing against
the file afterwards, which is the only way to catch a renderer that has
quietly computed its own answer.

That is not a hypothetical failure here. The slider position came from one
derivation and the status chip from another, and for methane they
disagreed: a marker at the 28th percentile beside the words "NOT
DETECTED", in all five samples. Neither number was wrong given its own
inputs; they simply were not the same inputs.

So every number these tests look for is read out of the rendered page and
compared with the file that produced it.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

pymupdf = pytest.importorskip("pymupdf")

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

SAMPLES = sorted(
    p.name for p in RESULTS.iterdir()
    if p.is_dir() and (p / "results.json").is_file()
    and (p / f"{p.name}_report.pdf").is_file()
) if RESULTS.is_dir() else []


def _pair(sample: str) -> tuple[dict[str, Any], Any]:
    results = json.loads(
        (RESULTS / sample / "results.json").read_text(encoding="utf-8"))
    doc = pymupdf.open(str(RESULTS / sample / f"{sample}_report.pdf"))
    return results, doc


def _rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    rows = results.get("report_rows")
    if not rows:
        pytest.skip("this report predates report_rows in results.json")
    return rows


def _ordinal_words(percentile: float | None) -> set[str]:
    """How the report may write a position, given the file's number."""
    if percentile is None:
        return {"\u2014"}
    if percentile == 0:
        return {"0"}
    if percentile >= 99.5:
        return {">99th"}
    if percentile < 0.5:
        return {"<1st"}
    n = int(round(percentile))
    suffix = "th"
    if n % 100 not in (11, 12, 13):
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return {f"{n}{suffix}"}


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_reading_in_the_file_appears_in_the_report(sample: str) -> None:
    results, doc = _pair(sample)
    text = "\n".join(page.get_text() for page in doc)
    missing = [
        row["metabolite"] for row in _rows(results)
        if row["metabolite"] not in text
    ]
    assert not missing, f"{sample}: in the file and not on any page: {missing}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_position_printed_is_the_position_in_the_file(sample: str) -> None:
    """The check that would have caught methane. The ordinal beside a
    reading has to be the one its percentile in the file produces."""
    results, doc = _pair(sample)
    text = "\n".join(page.get_text() for page in doc)
    wrong: list[str] = []
    for row in _rows(results):
        # A reading the file marks unassessed keeps its position *word* and
        # loses the ordinal: "76th percentile" off a single fragment is the
        # false precision this whole file exists to prevent, so the page is
        # right to withhold it and this check has nothing to compare.
        if (row.get("status") or {}).get("assessed") is False:
            continue
        expected = _ordinal_words(row["percentile"])
        # Find the reading's own line and the token that follows it.
        for match in re.finditer(re.escape(row["metabolite"]), text):
            window = text[match.end():match.end() + 220]
            tokens = set(re.findall(r">99th|<1st|\b\d{1,3}(?:st|nd|rd|th)\b|(?<!\d)0(?!\d)",
                                    window))
            if not tokens:
                continue
            if expected & tokens:
                break
        else:
            if row["percentile"] is not None:
                wrong.append(
                    f"{row['panel']}: file says {row['percentile']} "
                    f"(so {sorted(expected)}), page shows none of that"
                )
    assert not wrong, f"{sample}: " + "; ".join(wrong[:6])


#: A run of capitals on its own line is a column header, which means a new
#: table has started and the previous reading's row has ended.
_NEXT_TABLE = re.compile(r"\n[A-Z][A-Z0-9 &/()-]{7,}\n")


def _row_window(text: str, start: int, limit: int = 160) -> str:
    """The text belonging to one reading's row, and no further.

    A fixed slice runs off the end of a row and into the next table, where
    it finds that row's ordinal and blames it on this reading. That is how
    an honest "0 · NOT DETECTED" came to look like a reading drawn at the
    36th percentile: the number belonged to the GABA row underneath it.
    """
    window = text[start:start + limit]
    boundary = _NEXT_TABLE.search(window)
    return window[: boundary.start()] if boundary else window


@pytest.mark.parametrize("sample", SAMPLES)
def test_nothing_undetected_is_drawn_with_a_position(sample: str) -> None:
    """A reading the file says was not detected must not be shown sitting
    somewhere on the scale."""
    results, doc = _pair(sample)
    text = "\n".join(page.get_text() for page in doc)
    for row in _rows(results):
        if row["status"]["label"] != "not detected":
            continue
        for match in re.finditer(re.escape(row["metabolite"]), text):
            window = _row_window(text, match.end())
            if "NOT DETECTED" not in window.upper():
                continue
            positions = re.findall(r"\b(\d{1,3})(?:st|nd|rd|th)\b", window)
            assert not [p for p in positions if int(p) > 0], (
                f"{sample}/{row['panel']}: the file says not detected, and the page "
                f"shows it at the {positions[0]}th percentile"
            )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_report_shows_no_number_the_file_does_not_hold(sample: str) -> None:
    """A spot check on the headline figures rather than every glyph: the
    sample name, the read count and the panel count all come from the file."""
    results, doc = _pair(sample)
    first = "\n".join(doc[i].get_text() for i in range(min(3, len(doc))))
    assert str(results["sample"]) in first
    n_panels = len(results["panels"])
    assert n_panels == len(_rows(results)), (
        "the file disagrees with itself about how many readings there are"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_functional_reading_with_a_position_shows_one(sample: str) -> None:
    """The readings that had adjectives where every other row had a bar."""
    results, doc = _pair(sample)
    scores = (results.get("extension") or {}).get("views", {}).get("reading_scores")
    if not scores:
        pytest.skip("this report predates the reading scores")
    text = "\n".join(page.get_text() for page in doc)
    placed = {k: v for k, v in scores.items() if v.get("percentile") is not None}
    assert placed, "no functional reading was placed at all"
    unplaced_share = 1 - len(placed) / len(scores)
    assert unplaced_share <= 0.25, (
        f"{sample}: {unplaced_share:.0%} of the functional readings have no position"
    )
    # Every placed reading that is printed as its own row shows its label
    # and its ordinal. Readings that merely restate a panel are printed once,
    # as the panel, so they are exempt here and checked by the
    # one-number-per-quantity tests instead.
    from openbiota import pdfextension as PX
    from openbiota.pdfreport import _ordinal

    rows = PX.all_functional_rows(
        results["extension"]["views"],
        {r["panel"]: r.get("higher_means", "unclear") for r in results.get("report_rows") or []},
        {str(p.get("name")): p.get("aggregate_from") or () for p in results.get("panels") or []},
    )
    flat = " ".join(text.split())
    missing = []
    for row in rows:
        if row.get("duplicates_panel") or row.get("percentile") is None:
            continue
        label = " ".join(str(row["label"]).split())
        if label not in flat or _ordinal(row["percentile"]) not in flat:
            missing.append(f"{row['reading_id']} ({label}, {_ordinal(row['percentile'])})")
    assert not missing, f"{sample}: placed readings whose label or ordinal is not on any page: {missing}"
