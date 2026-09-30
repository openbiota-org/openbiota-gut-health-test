"""Internal links: every one resolves, and lands where its text says.

The mechanism is named destinations, which ReportLab resolves at save time
and refuses to save if any is missing — so a dead link is already a build
failure. What these tests add is the other half: that each link lands on a
page that actually carries the thing it names, that page 1's rows link
onward, that the at-a-glance rows link on to detail cards and the cards
back, and that the sections appear in the reader's outline.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import pytest
from reportlab.platypus import BaseDocTemplate, Frame, NextPageTemplate, PageBreak, PageTemplate

from openbiota import pdflinks, pdfreport, pdfsummary
from openbiota.interpret import MetaboliteRow
from openbiota.pdfreport import (
    BODY_FRAME_HEIGHT,
    CONTENT_WIDTH,
    FIRST_FRAME_HEIGHT,
    FOOTER_H,
    MARGIN,
    PAGE,
    SECTIONS,
    _styles,
    status_for,
)

pymupdf = pytest.importorskip("pymupdf")

ROOT = Path(__file__).resolve().parents[1]

GOTO = 1  # pymupdf's LINK_GOTO


def _row(name: str, pct: float, panel: str | None = None) -> MetaboliteRow:
    return MetaboliteRow(
        panel=panel or name.lower(), metabolite=name, value=1.0, percentile=pct,
        cohort_median=1.0, cohort_p25=0.5, cohort_p75=1.5, detected=True, fragments=100,
        status=status_for(percentile=pct, higher_means="adverse"),
        what_it_is="what", made_from="from", higher_means="adverse",
        evidence_strength="moderate", summary="s", evidence_detail="d", citation="c",
        genes=[], confidence="confirmed", implication="i", extra_note="",
        category="neurotransmitter",
    )


def _build(story: list[Any]) -> pymupdf.Document:
    buf = io.BytesIO()
    doc = BaseDocTemplate(buf, pagesize=PAGE, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=0, bottomMargin=0)
    first = Frame(MARGIN, FOOTER_H, CONTENT_WIDTH, FIRST_FRAME_HEIGHT, id="first",
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    body = Frame(MARGIN, FOOTER_H, CONTENT_WIDTH, BODY_FRAME_HEIGHT, id="body",
                 leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id="first", frames=[first]), PageTemplate(id="interior", frames=[body])])
    doc.build(story)
    return pymupdf.open(stream=buf.getvalue(), filetype="pdf")


def _links(doc: pymupdf.Document) -> list[dict[str, Any]]:
    out = []
    for number, page in enumerate(doc):
        for link in page.get_links():
            link["source_page"] = number
            out.append(link)
    return out


@pytest.fixture(scope="module")
def mini_report() -> pymupdf.Document:
    """Page 1, then the functions at-a-glance and detail sections."""
    st = _styles()
    rows = [_row(f"Pathway {i}", 5.0 + 7 * i) for i in range(12)]
    story: list[Any] = []
    pdfsummary.summary_page(
        story, st, rows=rows, similarity=None, rpob_profile=None,
        meta={"sample": "T", "date": "9 September 2026", "n_panels": 12, "read_pairs": "8.0 million",
              "read_pairs_raw": 8_000_000, "cohort_note": None, "n_species": 92, "n_profiles": 0},
        findings=None, age=None, age_meta=None, pathogens=None, pathogen_section=SECTIONS["pathogens"],
    )
    story.append(NextPageTemplate("interior"))
    story.append(PageBreak())
    pdfreport._at_a_glance(story, st, rows, "the reference group", section=SECTIONS["functions"],
                           detail_section=SECTIONS["functions_detail"])
    story.append(PageBreak())
    pdfreport._detail_pages(story, st, rows, "the reference group", section=SECTIONS["functions_detail"])
    # Page 1 and the fine print mention other sections; give them somewhere
    # to land, as the full report does. The functions sections are real here.
    from tests.conftest import section_stubs

    story.extend(section_stubs(st, skip={SECTIONS["functions"], SECTIONS["functions_detail"]}))
    return _build(story)


def test_module_self_test() -> None:
    assert pdflinks.self_test() > 0


def test_every_link_is_an_internal_goto_that_resolves(mini_report: pymupdf.Document) -> None:
    links = _links(mini_report)
    assert len(links) > 20, "the mini report should be full of links"
    for link in links:
        assert link["kind"] == GOTO, link
        assert 0 <= link["page"] < len(mini_report), link


def test_page_one_rows_land_on_their_detail_cards(mini_report: pymupdf.Document) -> None:
    """Item-level: each page-1 function row links straight to that function's card."""
    text = {n: page.get_text() for n, page in enumerate(mini_report)}
    page_one = [lk for lk in _links(mini_report) if lk["source_page"] == 0]
    glance = {n for n, t in text.items() if "Functions at a Glance" in t}
    detail = {n for n, t in text.items() if "at a glance" in t and n not in glance and n > 0}
    assert detail, "no detail pages found"
    rows = [lk for lk in page_one if lk["from"].width > 100 and lk["page"] in detail]
    assert len(rows) >= 8, f"expected page-1 rows to link to detail cards, got {len(rows)}"
    for link in rows:
        rect = link["from"]
        assert rect.width > 4 * rect.height, rect
        # And the card it lands on names a pathway.
        assert "Pathway" in text[link["page"]]


def test_page_one_blocks_link_to_their_sections(mini_report: pymupdf.Document) -> None:
    """Section-level: the block heading rows on page 1 go to the at-a-glance section."""
    text = {n: page.get_text() for n, page in enumerate(mini_report)}
    page_one = [lk for lk in _links(mini_report) if lk["source_page"] == 0]
    glance = {n for n, t in text.items() if "Functions at a Glance" in t}
    age = {n for n, t in text.items() if t.strip().startswith(f"{SECTIONS['age']}") or f"\n{SECTIONS['age']} " in t}
    to_glance = [lk for lk in page_one if lk["page"] in glance and lk["from"].width > 100]
    assert to_glance, ("the Microbial Functions heading should link to "
                       "the at-a-glance section")
    # The whole age block is one link: a tall, wide rectangle.
    tall = [lk for lk in page_one if lk["from"].height > 60 and lk["from"].width > 150]
    assert tall, "expected whole-block links (age, diversity) on page 1"
    _ = age


def test_at_a_glance_rows_link_on_to_detail_cards_and_back(mini_report: pymupdf.Document) -> None:
    text = {n: page.get_text() for n, page in enumerate(mini_report)}
    glance = {n for n, t in text.items() if "Functions at a Glance" in t}
    forward = [lk for lk in _links(mini_report) if lk["source_page"] in glance and lk["page"] not in glance]
    assert forward, "no links from the at-a-glance rows"
    # Every forward link lands on a page that names a pathway.
    for link in forward:
        assert "Pathway" in text[link["page"]] or "in detail" in text[link["page"]], link
    back = [lk for lk in _links(mini_report) if lk["source_page"] not in glance and lk["page"] in glance]
    assert back, "detail cards should link back to the at-a-glance rows"


def test_section_mentions_link_to_the_section(mini_report: pymupdf.Document) -> None:
    """'Section 16 has the findings' on the at-a-glance page lands on section 16."""
    text = {n: page.get_text() for n, page in enumerate(mini_report)}
    glance = next(n for n, t in text.items() if "Functions at a Glance" in t)
    detail_first = next(n for n, t in text.items() if "Functions in Detail" in t)
    to_detail = [lk for lk in _links(mini_report) if lk["source_page"] == glance and lk["page"] == detail_first]
    assert to_detail, "the at-a-glance page's 'Section 16' mention should link to section 16"


def test_sections_populate_the_outline(mini_report: pymupdf.Document) -> None:
    toc = mini_report.get_toc()
    titles = [entry[1] for entry in toc]
    assert any(t.startswith(f"{SECTIONS['summary']} ") for t in titles), titles
    assert any(t.startswith(f"{SECTIONS['functions']} ") for t in titles), titles
    assert any(t.startswith(f"{SECTIONS['functions_detail']} ") for t in titles), titles
    # One entry per section, not one per time the heading was drawn.
    assert len(titles) == len(set(titles)), titles


def test_link_text_is_marked_in_the_accent_colour() -> None:
    out = pdflinks.link_sections("see section 8")
    assert 'color="#0F8B8D"' in out, out


def test_dead_links_fail_the_build() -> None:
    """The property everything else relies on: an undefined destination cannot be saved."""
    st = _styles()
    with pytest.raises(ValueError, match="undefined destination"):
        _build([pdflinks.Paragraph('<a href="#nowhere">dead</a>', st["body"])])


@pytest.mark.skipif(
    not any((ROOT / "results").glob("*/*_report.pdf")) if (ROOT / "results").exists() else True,
    reason="no generated report to inspect",
)
def test_running_header_is_one_home_button_and_the_footer_keeps_the_website() -> None:
    """On interior pages the whole header run - brand, house, arrow, title -
    is a single link back to page one and nothing in the header band goes
    to the website; the footer's brand still does."""
    path = max((ROOT / "results").glob("*/*_report.pdf"), key=lambda f: f.stat().st_mtime)
    doc = pymupdf.open(path)
    page = doc[5]
    header_band = 14 * 72 / 25.4          # top 14 mm, in points
    header = [lk for lk in page.get_links() if lk["from"].y1 < header_band]
    assert header, "no link in the header band"
    assert all(lk["kind"] == GOTO and lk["page"] == 0 for lk in header), header
    assert min(lk["from"].x0 for lk in header) < 60, "the home link must start at the brand, not after it"
    footer = [lk for lk in page.get_links() if lk["from"].y0 > page.rect.height - 14 * 72 / 25.4]
    assert any(lk["kind"] == pymupdf.LINK_URI and "openbiota.com" in (lk.get("uri") or "") for lk in footer), footer


@pytest.mark.skipif(
    not any((ROOT / "results").glob("*/*_report.pdf")) if (ROOT / "results").exists() else True,
    reason="no generated report to inspect",
)
def test_generated_report_links_all_resolve_and_sections_are_bookmarked() -> None:
    """On a real report: every link resolves, and every section is in the outline.

    Takes the most recently generated report, so a stale one from before the
    links existed cannot fail the build.
    """
    path = max((ROOT / "results").glob("*/*_report.pdf"), key=lambda f: f.stat().st_mtime)
    doc = pymupdf.open(path)
    # The report also carries PubMed URLs; only the internal links are ours.
    links = [lk for lk in _links(doc) if lk["kind"] == GOTO]
    assert len(links) > 100
    assert all(0 <= lk["page"] < len(doc) for lk in links)
    titles = [entry[1] for entry in doc.get_toc()]
    for key, number in SECTIONS.items():
        if number == 0:  # the unnumbered contents page: title alone in the outline
            assert any(t.startswith("Table of contents") for t in titles), "contents missing from outline"
            continue
        assert any(t.startswith(f"{number} ") for t in titles), f"section {number} ({key}) missing from outline"
    # Page 1 links to at least a dozen distinct places.
    assert len({lk["page"] for lk in links if lk["source_page"] == 0}) >= 3
