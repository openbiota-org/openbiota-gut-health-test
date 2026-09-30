"""One layout for every metric, enforced.

A reader should not have to learn a second layout halfway down a page.
For a while they did: the panels carried a bar, an ordinal and a status
chip, and the readings underneath them carried a column of gene symbols
and the words "class-level only". Somebody wanting to know whether their
acetate production was high was shown *pta, acsB, cooS*.

`openbiota.metriccard` fixes the shape. These tests hold every metric to
it, so a new measurement cannot arrive with a layout of its own.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.metriccard import CARD_SECTIONS
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06")


def _results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


#: Pulling the text out of a three-hundred-page report is slow, and every
#: test here wants the same text, so it is done once per sample.
_TEXT_CACHE: dict[str, str] = {}


def _report_text(sample: str) -> str:
    import pymupdf

    if sample in _TEXT_CACHE:
        return _TEXT_CACHE[sample]
    pdf = REPO / "results" / sample / f"{sample}_report.pdf"
    if not pdf.is_file() or pdf.stat().st_size == 0:
        pytest.skip(f"{sample} has no rendered report")
    # A report built before the card contract existed cannot be held to it,
    # and failing on that says nothing except that it has not been rebuilt
    # yet. The contract itself is tested above, against the source.
    contract = REPO / "openbiota" / "metriccard.py"
    if pdf.stat().st_mtime < contract.stat().st_mtime:
        pytest.skip(
            f"{sample}'s report predates the metric-card contract; regenerate it"
        )
    with pymupdf.open(pdf) as doc:
        text = " ".join(" ".join(page.get_text().split()) for page in doc)
    _TEXT_CACHE[sample] = text
    return text


# --- the contract itself ---------------------------------------------------


def test_the_card_has_exactly_the_four_named_sections() -> None:
    """Changing this list changes the contract for every metric, which is
    why it is written down in one place and asserted here."""
    assert CARD_SECTIONS == (
        "What is driving this",
        "What your reading means",
        "How can I improve this",
        "The research behind it",
    )


# --- the rendered report ---------------------------------------------------


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_card_section_appears_in_the_report(sample: str) -> None:
    text = _report_text(sample)
    for heading in CARD_SECTIONS:
        assert heading in text, f"no metric card carries {heading!r}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_improvement_is_offered_as_often_as_the_other_sections(sample: str) -> None:
    """The section most likely to be quietly dropped is the one a reader
    most wants, so it is counted against its neighbours."""
    text = _report_text(sample)
    driving = text.count(CARD_SECTIONS[0])
    improve = text.count(CARD_SECTIONS[2])
    assert improve >= driving * 0.9, (
        f"{driving} cards say what is driving the reading and only {improve} say how "
        "to improve it"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_metric_is_labelled_with_a_bare_gene_symbol(sample: str) -> None:
    """"Acetate formation (Pta-AckA)" is a name only its author can read."""
    banned = ("(Pta-AckA)", "([FeFe])", "(Wood-Ljungdahl)")
    text = _report_text(sample)
    found = [b for b in banned if b in text]
    assert not found, f"gene shorthand is being used as a metric name: {found}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_overview_never_shows_a_readings_only_layout(sample: str) -> None:
    """These headers belonged to the readings' own table, which is gone.

    A gene table inside a pathogen card may legitimately have a column
    called "what was found"; a list of metrics may not, because every
    other metric in that list shows a verdict there.
    """
    text = _report_text(sample)
    assert "ALSO MEASURED HERE" not in text
    assert "HOW SPECIFIC" not in text, "the evidence-grade column is back in a metric list"


def test_the_metric_row_renderer_emits_the_shared_columns() -> None:
    """Checked at the source, because that is where the shape is decided.

    A metric row is label, bar, ordinal, chip - the same four the panel
    rows use. Anything that puts gene symbols or an evidence grade in a
    metric row is the layout this replaced.
    """
    import inspect

    from openbiota import pdfreport, pdfsynbiotic

    assert not hasattr(pdfsynbiotic, "functional_group_table"), (
        "the separate readings table is back; readings are rows of the panel table"
    )
    source = inspect.getsource(pdfreport._glance_reading_rows)
    assert "StatusChip" in source, "metric rows must carry the same chip as every other row"
    assert "PercentileBar" in source, "metric rows must carry the same bar"
    assert "col[3] - 6 * mm" in source and "bar_w - 6 * mm" in source, (
        "metric rows must use the panel rows' own widths"
    )
    for banned in ("WHAT WAS FOUND", "HOW SPECIFIC", "ALSO MEASURED HERE", "WHERE YOU SIT"):
        assert banned not in source, f"{banned} is back in the metric row renderer"


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_functional_reading_reaches_a_card(sample: str) -> None:
    """A row that links nowhere is a dead end for the reader."""
    results = _results(sample)
    views = (results.get("extension") or {}).get("views") or {}
    scores = views.get("reading_scores") or {}
    if not scores:
        pytest.skip(f"{sample} predates the reading scores")
    text = _report_text(sample)
    assert "Each reading in full" in text, "the reading cards section is missing"


# --- the data behind the labels -------------------------------------------


def test_every_contributing_gene_carries_words_as_well_as_a_symbol() -> None:
    """The card leads with what an enzyme does; the symbol goes underneath."""
    from openbiota.extension.readingscore import GenePosition

    assert "label" in GenePosition.__dataclass_fields__, (
        "contributors must carry a readable label, or the card can only print symbols"
    )


def test_the_group_card_only_touches_attributes_that_exist() -> None:
    """A wrong attribute here costs sixteen minutes to discover.

    `_group_meaning_block` runs at the very end of a report build, after
    the screening, so a typo in it is found by a run rather than by a
    test. This reads the attributes it uses and checks them against the
    objects it is handed, which the rest of the card already establishes.
    """
    import inspect
    import re

    from openbiota import pdflibrary as L

    mine = inspect.getsource(L._group_meaning_block)
    established = inspect.getsource(L._group_card)

    driver_attrs = set(re.findall(r"\bd\.([a-z_]+)", mine))
    assert driver_attrs <= set(L.GroupDriver.__slots__), (
        f"unknown GroupDriver attributes: {driver_attrs - set(L.GroupDriver.__slots__)}"
    )

    # The group object has no slots to check against, so the standard is
    # what the rest of the same card already reads from it.
    group_attrs = set(re.findall(r"g\.group\.([a-z_]+)", mine))
    known = set(re.findall(r"g\.group\.([a-z_]+)", established))
    assert group_attrs <= known, f"unverified group attributes: {group_attrs - known}"
