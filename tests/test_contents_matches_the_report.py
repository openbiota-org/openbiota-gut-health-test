"""The contents must name the sections the report actually prints.

An entry whose title differs from the heading it links to sends a reader
somewhere they did not choose, and the report had several: section 1 was
listed as "What stands out" when that is a different section, and three
overview sections were missing from the list altogether.

Numbers are not written into the contents - they are read from
`SECTIONS` - so these tests are about the titles and the coverage.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from openbiota.pdfatlas import CONTENTS_ENTRIES
from openbiota.pdfreport import SECTIONS

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("SAMPLE5_A05", "SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06")


def test_every_contents_entry_names_a_real_section() -> None:
    for key, title, _ in CONTENTS_ENTRIES:
        assert key in SECTIONS, f"{title!r} points at {key!r}, which is not a section"


def test_the_contents_are_in_page_order() -> None:
    """A list that jumps around is harder to use than no list."""
    numbers = [SECTIONS[key] for key, _, _ in CONTENTS_ENTRIES]
    assert numbers == sorted(numbers), f"out of order: {numbers}"


def test_every_overview_section_is_listed() -> None:
    """The detail and methods sections are reachable from their overview;
    the overviews themselves have to be reachable from here."""
    listed = {key for key, _, _ in CONTENTS_ENTRIES}
    overviews = {
        "summary", "stood_out", "community", "groups", "organisms", "catalogue",
        "pathogens", "age", "functions", "patterns", "biofilm", "mycobiome", "skin",
        "actions", "context",
    }
    assert overviews <= listed, f"missing from the contents: {sorted(overviews - listed)}"


def test_no_entry_is_a_duplicate() -> None:
    keys = [key for key, _, _ in CONTENTS_ENTRIES]
    assert len(keys) == len(set(keys))


def test_the_contents_page_comes_before_what_stood_out() -> None:
    """A reader meets the contents first, then the findings.

    The contents page is unnumbered (0) and is printed between the summary
    and what stood out, so the page order is checked through the renderer's
    own sequence rather than through the numbers.
    """
    assert SECTIONS["guide"] == 0
    assert SECTIONS["summary"] == 1
    assert SECTIONS["stood_out"] == 2
    src = (REPO / "openbiota" / "pdfreport.py").read_text()
    assert src.index('reading_guide_page(') < src.index('insights_page('), (
        "the contents page must be built before what stood out"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_each_listed_title_appears_as_a_heading_in_the_report(sample: str) -> None:
    """The strongest form of the check: the words in the contents are the
    words at the top of the section they lead to."""
    import pymupdf

    pdf = REPO / "results" / sample / f"{sample}_report.pdf"
    if not pdf.is_file() or not pdf.stat().st_size:
        pytest.skip(f"{sample} has no rendered report")
    contract = REPO / "openbiota" / "pdfatlas.py"
    if pdf.stat().st_mtime < contract.stat().st_mtime:
        pytest.skip(f"{sample}'s report predates this contents list; regenerate it")

    with pymupdf.open(pdf) as doc:
        text = " ".join(" ".join(page.get_text().split()) for page in doc)

    missing = []
    for key, title, _ in CONTENTS_ENTRIES:
        number = SECTIONS[key]
        # The heading is printed as "<number> <title>", so look for both
        # together rather than the title anywhere on any page.
        if not re.search(rf"\b{number}\s+{re.escape(title)}", text):
            missing.append(f"{number} {title}")
    assert not missing, f"listed but not printed as a heading: {missing}"
