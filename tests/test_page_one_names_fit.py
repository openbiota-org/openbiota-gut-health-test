"""Every name that can reach page 1 sits on one line there. Measured.

Page 1 is a grid. Two of its rows wrapped - "Dietary fibre and
carbohydrate breakdown" and "Glucosinolate to isothiocyanate conversion" -
because the short-name rule was a 30-character guess and those names were
never in the table. This file makes the rule what it was always meant to
be: a measurement of every candidate name, in the font and size page 1
uses, against the width its column actually has.

"Every candidate" means every panel in `panels/` and every profile in
`profiles/`, not the ten that happened to score highly in one sample. Any
of them can be on page 1 next time.
"""

from __future__ import annotations

import json
from pathlib import Path

import pymupdf
import pytest

from openbiota.panels import load_panel_set
from openbiota.pdfsummary import (
    COLUMN_GUTTER,
    CONTENT_WIDTH,
    _one_line,
    page_one_name_budget,
)
from openbiota.profiles import load_profile_set
from openbiota.shortnames import (
    _METABOLITE_SHORT,
    _PROFILE_SHORT,
    PAGE_ONE_MARGIN,
    fits_page_one,
    short_metabolite,
    short_profile_label,
    text_width,
)

REPO = Path(__file__).resolve().parent.parent
PANELS = load_panel_set(REPO / "panels")
PROFILES = load_profile_set(REPO / "profiles")
RESULTS = sorted(
    p for p in (REPO / "results").glob("*/results.json")
    if not p.parent.name.startswith(("_", "fmt_"))
)


def _panels() -> list:
    return list(getattr(PANELS, "panels", PANELS))


def _profiles() -> list:
    return list(getattr(PROFILES, "profiles", PROFILES))


# --------------------------------------------------------------------------- #
# the rule, at the source
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("panel", _panels(), ids=lambda p: p.name)
def test_every_panel_has_a_one_line_name_for_page_one(panel) -> None:
    """Function rows carry a direction arrow, so the arrow's width is spent
    before the name is measured."""
    budget = page_one_name_budget(arrow=True)
    short = short_metabolite(panel.metabolite)
    assert fits_page_one(short, budget), (
        f"{panel.name}: {short!r} is {text_width(short):.1f}pt and page 1 allows "
        f"{budget - PAGE_ONE_MARGIN:.1f}pt with the arrow. Add a shorter form to "
        f"openbiota.shortnames._METABOLITE_SHORT under {panel.metabolite!r}."
    )


@pytest.mark.parametrize("profile", _profiles(), ids=lambda p: p.name)
def test_every_profile_has_a_one_line_name_for_page_one(profile) -> None:
    budget = page_one_name_budget(arrow=False)
    short = short_profile_label(profile.name, profile.label)
    assert fits_page_one(short, budget), (
        f"{profile.name}: {short!r} is {text_width(short):.1f}pt and page 1 allows "
        f"{budget - PAGE_ONE_MARGIN:.1f}pt. Add a shorter form to "
        f"openbiota.shortnames._PROFILE_SHORT under {profile.name!r}."
    )


def test_the_fixed_page_one_rows_fit_their_column() -> None:
    """The biofilm, mycobiome and skin rows have hand-written labels in a
    column of their own; they are measured against that column."""
    col = (CONTENT_WIDTH - COLUMN_GUTTER) / 2
    budget = col * 0.44  # `_biofilm_block` and its siblings: no cell padding
    for label in (
        "Biofilm-forming (bad)", "Protective gut-lining (good)",
        "Myco-Score (experimental)", "Gut-skin axis (beta)",
    ):
        width = text_width(label, font="Helvetica", size=7.8)
        assert width + PAGE_ONE_MARGIN <= budget, f"{label!r} is {width:.1f}pt of {budget:.1f}"


def test_no_short_name_is_keyed_to_a_name_nothing_uses() -> None:
    """A renamed panel would keep its old key here and silently lose its
    short form. Every key must be a live display name or profile id."""
    displays = {p.metabolite for p in _panels()}
    dead = sorted(k for k in _METABOLITE_SHORT if k not in displays)
    assert not dead, f"short names keyed to no panel: {dead}"
    names = {p.name for p in _profiles()}
    dead = sorted(k for k in _PROFILE_SHORT if k not in names)
    assert not dead, f"short names keyed to no profile: {dead}"


def test_the_budget_is_the_real_geometry() -> None:
    """The number the table is measured against is the page's own column,
    not a constant that could drift from it."""
    col = (CONTENT_WIDTH - COLUMN_GUTTER) / 2
    assert 80 < page_one_name_budget(arrow=False) < col * 0.5
    assert page_one_name_budget(arrow=True) < page_one_name_budget(arrow=False)


# --------------------------------------------------------------------------- #
# the backstop in the renderer
# --------------------------------------------------------------------------- #


def test_an_over_budget_name_is_shrunk_rather_than_wrapped() -> None:
    budget = page_one_name_budget(arrow=False)
    long = "<b>A name that is far too long to sit on one line of the page</b>"
    out = _one_line(long, budget)
    assert out.startswith("<font size="), "an over-budget name must be drawn smaller"
    fine = "<b>Butyrate</b>"
    assert _one_line(fine, budget) == fine, "a name that fits is left alone"


def test_the_row_renderer_goes_through_the_backstop() -> None:
    src = (REPO / "openbiota" / "pdfsummary.py").read_text()
    assert "Paragraph(_one_line(name, name_avail), st[\"cell\"])" in src, (
        "_scale_rows must draw every name through _one_line"
    )


# --------------------------------------------------------------------------- #
# the page itself
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_no_metric_name_on_page_one_spans_two_lines(path: Path) -> None:
    """The check a reader would make: every metric name printed on page 1
    is on one line. Read from the rendered PDF's own line structure."""
    pdf = next(path.parent.glob("*_report.pdf"), None)
    if pdf is None:
        pytest.skip("no report")
    results = json.loads(path.read_text())
    with pymupdf.open(pdf) as doc:
        page = doc[0]
        lines = [
            "".join(span["text"] for span in line["spans"]).strip()
            for block in page.get_text("dict")["blocks"] if block.get("type") == 0
            for line in block["lines"]
        ]
    joined = "\n".join(lines)

    flat = " ".join(joined.split())

    # A long name whose short form differs must never be on page 1 at all:
    # its presence means the short-name lookup was bypassed, which is how
    # the two-line rows got there.
    leaked = []
    for r in results.get("report_rows") or []:
        long, short = r["metabolite"], short_metabolite(r["metabolite"])
        if long != short and " ".join(long.split()) in flat:
            leaked.append(long)
    assert not leaked, f"{path.parent.name}: long names printed on page 1: {leaked}"

    names = [short_metabolite(r["metabolite"]) for r in results.get("report_rows") or []]
    names += [
        short_profile_label(p.get("profile") or p.get("name") or "", p.get("label") or "")
        for p in (results.get("profile_similarity") or {}).get("ranked") or []
    ]
    broken = []
    for name in names:
        if name not in flat:
            continue  # not on page 1 this time
        if not any(name in line for line in lines):
            broken.append(name)
    assert not broken, f"{path.parent.name}: names split across lines on page 1: {broken}"
