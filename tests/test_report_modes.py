"""The report in its two modes, and the page break that never leaves a page blank.

A lite report leaves out the sections that need the full install and
numbers the rest without gaps; page one carries no pathogen band and no
"strain typing did not run" tile; the title says which it is. These tests
render the pieces that change and read them back, and check the section
numbering as every module sees it.
"""

from __future__ import annotations

import io
from typing import Any

import pytest
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer

from openbiota import pdfpathogens, pdfreport, reportmode
from openbiota.pdflinks import SoftPageBreak
from tests.test_pathogens import _flat_text, _styles

# --------------------------------------------------------------------------- #
# SoftPageBreak
# --------------------------------------------------------------------------- #


def _pages(story: list[Any]) -> int:
    import pymupdf

    buf = io.BytesIO()
    SimpleDocTemplate(buf, pagesize=A4).build(story)
    return len(pymupdf.open(stream=buf.getvalue(), filetype="pdf"))


def _p(text: str = "x") -> Paragraph:
    return Paragraph(text, getSampleStyleSheet()["Normal"])


def test_soft_break_at_the_top_of_a_page_does_nothing():
    assert _pages([SoftPageBreak(), _p()]) == 1
    assert _pages([_p("a"), SoftPageBreak(), SoftPageBreak(), _p("b")]) == 2
    assert _pages([_p("a"), PageBreak(), SoftPageBreak(), _p("b")]) == 2


def test_soft_break_still_breaks_a_page_with_content():
    assert _pages([_p("a"), SoftPageBreak(), _p("b")]) == 2


def test_the_blank_page_the_hard_break_made_is_gone():
    """Content that ends exactly at the page foot, then a break: a hard
    PageBreak turns a second, empty page; the soft one does not."""
    frame_h = A4[1] - 2 * 72
    assert _pages([Spacer(1, frame_h - 12), SoftPageBreak(), _p("b")]) == 2
    assert _pages([_p("a"), PageBreak(), PageBreak(), _p("b")]) == 3   # the bug, for contrast


def test_the_report_modules_use_no_hard_page_breaks():
    """One slip reintroduces a blank page; the report modules only soft-break."""
    import re
    from pathlib import Path

    root = Path(pdfreport.__file__).parent
    for module in root.glob("pdf*.py"):
        code = "\n".join(line for line in module.read_text().splitlines() if not line.lstrip().startswith("#"))
        assert not re.search(r"(?<![A-Za-z])PageBreak\(\)", code), f"{module.name} uses a hard PageBreak()"


# --------------------------------------------------------------------------- #
# Section numbering
# --------------------------------------------------------------------------- #


@pytest.fixture(autouse=True)
def _full_numbering_afterwards():
    yield
    pdfreport.activate_sections()


def test_lite_numbering_is_contiguous_and_restorable():
    omitted = reportmode.omitted_sections("lite")
    s = pdfreport.activate_sections(omitted)
    assert not (omitted & set(s)), "omitted sections have no number"
    numbers = sorted(n for n in s.values() if n)
    assert numbers == list(range(1, len(numbers) + 1)), numbers
    assert s["guide"] == 0
    # The order is the canonical one, so "actions" still follows "skin".
    assert s["actions"] == s["skin"] + 1
    assert s["technical"] == len(numbers)
    pdfreport.activate_sections()
    assert pdfreport.SECTIONS == pdfreport.FULL_SECTIONS


def test_every_module_sees_the_renumbering():
    """`SECTIONS` is read by name in a dozen modules; the renumbering is in
    place, so what they see is what the headings print."""
    from openbiota import pdfatlas, pdfsummary

    pdfreport.activate_sections(reportmode.omitted_sections("lite"))
    assert pdfsummary.SECTIONS["actions"] == pdfreport.SECTIONS["actions"] == 12
    from openbiota.pdfreport import SECTIONS as seen_by_a_fresh_import
    assert seen_by_a_fresh_import is pdfreport.SECTIONS
    # The contents table drops the absent keys rather than failing on them.
    entries = [key for key, _t, _b in pdfatlas.CONTENTS_ENTRIES if key in pdfreport.SECTIONS]
    assert "pathogens" not in entries and "mycobiome" not in entries
    assert "actions" in entries


def test_full_numbering_is_the_published_one():
    pdfreport.activate_sections()
    assert pdfreport.SECTIONS["pathogens"] == 7
    assert pdfreport.SECTIONS["actions"] == 14
    assert pdfreport.SECTIONS["technical"] == 28


# --------------------------------------------------------------------------- #
# Page one
# --------------------------------------------------------------------------- #


def test_lite_foot_is_the_solutions_panel_alone_at_the_pairs_height():
    from openbiota.pdfsummary import _foot

    st = _styles()
    lite = _foot(st, pathogens=None, pathogen_section=0, solutions_section=12,
                 n_problems=9, n_options=42, report_mode="lite")
    text = _flat_text(lite)
    assert "Possible solutions" in text and "Section 12" in text and "42" in text
    assert "pathogen" not in text.lower()
    height = sum(f.wrap(170 * mm, 400)[1] for f in lite)
    assert abs(height - pdfpathogens.MIN_PAIR_HEIGHT) < 0.5 * mm
    assert abs(lite[0].wrap(170 * mm, 400)[0] - 170 * mm) < 0.5 * mm, "the panel spans the width"


def test_full_foot_without_a_screen_still_says_it_did_not_run():
    """Full mode with a screen that produced nothing is an error state and
    is shown as one - that message is only wrong in a lite report."""
    from openbiota.pdfsummary import _foot

    full = _foot(_styles(), pathogens={"results": [], "coverage": {"assessed": 0, "not_assessed": 10}},
                 pathogen_section=7, solutions_section=14, n_problems=1, n_options=1, report_mode="full")
    assert "did not run" in _flat_text(full)


def test_solutions_panel_needs_a_section_to_point_at():
    story: list[Any] = []
    pdfpathogens.solutions_panel(story, _styles(), section=0)
    assert story == []


def test_the_title_carries_the_mode():
    from openbiota.pdfsummary import _head

    for mode in ("full", "lite"):
        head = _head(_styles(), rows=[], similarity=None, findings=None, meta={"read_pairs": "1"},
                     n_headlines=0, report_mode=mode)
        title = next(f for f in head if isinstance(f, Paragraph) and "Metagenomic Report" in f.text)
        assert f"> {mode}<" in title.text, title.text
        assert "<super>" in title.text


def test_the_lite_tile_grid_has_no_strain_tile():
    """In a lite report the sixth tile is the genus count, which every
    report has; "strain typing did not run" is a full-report statement."""
    from types import SimpleNamespace

    from openbiota.engines.metaphlan import TaxonomicProfile
    from openbiota.pdfsummary import _community_block

    taxonomy = TaxonomicProfile(
        index="mpa_v31", profiler_version="3.1.0",
        clades={
            "k__Bacteria|p__Bacteroidetes": 60.0,
            "k__Bacteria|p__Firmicutes": 40.0,
            "k__Bacteria|p__Bacteroidetes|g__Bacteroides|s__Bacteroides_vulgatus": 60.0,
            "k__Bacteria|p__Firmicutes|g__Faecalibacterium|s__Faecalibacterium_prausnitzii": 40.0,
        },
        unknown_percent=0.0, n_reads_processed=1000, elapsed_s=1.0, cached=True,
    )
    similarity = SimpleNamespace(taxonomy=taxonomy, community=None, references=None, results=[],
                                 matched_cohort=None)

    def tile_labels(node: Any) -> list[str]:
        """The labels of every Tile in a flowable tree."""
        out: list[str] = []
        if isinstance(node, pdfreport.Tile):
            out.append(node.label)
        for child in getattr(node, "_cellvalues", []) or []:
            for cell in child:
                out += tile_labels(cell)
        if isinstance(node, (list, tuple)):
            for item in node:
                out += tile_labels(item)
        return out

    for lite, forbidden, expected in ((True, "Resolved to strain", "Genera detected"),
                                      (False, "Genera detected", "Resolved to strain")):
        story: list[Any] = []
        _community_block(story, _styles(), similarity=similarity, rpob_profile=None, findings=None,
                         width=85 * mm, inv=None, lite_tiles=lite)
        labels = tile_labels(story)
        assert expected in labels, labels
        assert forbidden not in labels
        assert len(labels) == 6, "the grid keeps its shape"


def test_mode_label_words():
    assert reportmode.label("lite") == "lite" and reportmode.label("full") == "full"
    assert reportmode.FULL.name == "full" and reportmode.LITE.name == "lite"
    for m in (reportmode.FULL, reportmode.LITE):
        assert m.time_first and m.time_again and m.sections


def test_the_solution_problems_can_be_built_under_lite_numbering():
    """The first real lite run died here: a pathogen section number was
    read for a report that has no pathogen section."""
    import json
    from pathlib import Path

    from openbiota import pdfattention

    pdfreport.activate_sections(reportmode.omitted_sections("lite"))
    results_dir = Path(pdfreport.__file__).resolve().parents[1] / "results" / "EM1_AMD614"
    results = json.loads((results_dir / "results.json").read_text()) if (results_dir / "results.json").is_file() else {}
    results = dict(results, run={**(results.get("run") or {}), "report_mode": "lite"})
    problems = pdfreport._solution_problems(results, pdfattention.build(None, None), pdfreport.SECTIONS,
                                            draws_group_pages=False)
    assert isinstance(problems, list)
    for problem in problems:
        dest = getattr(problem, "detail_dest", None) or ""
        assert not dest.startswith("pth-") and not dest.startswith("myco-"), dest


def test_whitespace_at_the_top_of_a_page_is_dropped_and_cards_keep_their_rule():
    """The two ways a blank page was still made after the soft break: a
    trailing spacer spilling onto a fresh page so the next KeepTogether
    card moved on without it, and a card's separator rule left alone the
    same way. The report's frame drops the spacer; the rule travels at the
    head of its card and draws nothing at the top of a page."""
    from reportlab.platypus import KeepTogether, PageTemplate, Table

    from openbiota.pdflinks import ReportDocTemplate, ReportFrame, Separator

    def pages(story: list[Any]) -> int:
        import pymupdf

        buf = io.BytesIO()
        doc = ReportDocTemplate(buf, pagesize=A4, leftMargin=72, rightMargin=72, topMargin=72, bottomMargin=72)
        doc.addPageTemplates([PageTemplate(id="p", frames=[ReportFrame(
            72, 72, A4[0] - 144, A4[1] - 144, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)])])
        doc.build(story)
        return len(pymupdf.open(stream=buf.getvalue(), filetype="pdf"))

    frame_h = A4[1] - 144
    full = lambda: Table([["ink"]], rowHeights=[frame_h], colWidths=[100])  # noqa: E731
    card = [_p("x " * 40) for _ in range(40)]
    assert pages([full(), Spacer(1, 10), KeepTogether(card)]) == 2
    assert pages([full(), Spacer(1, 10), SoftPageBreak(), _p("b")]) == 2
    assert pages([full(), KeepTogether([Separator(A4[0] - 144), *card])]) == 2
    assert pages([_p("a"), KeepTogether([Separator(A4[0] - 144), _p("b")])]) == 1
    # A CondPageBreak is a Spacer subclass with a job and is never dropped.
    from reportlab.platypus import CondPageBreak
    assert pages([_p("a"), CondPageBreak(frame_h), _p("b")]) == 2
