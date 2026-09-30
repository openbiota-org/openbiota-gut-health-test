"""OpenBiota branding on the report.

The wordmark is a shipped package asset, which makes it the kind of thing
that silently disappears when packaging changes: `package-data` globs
`data/*.json`, and a PNG only travels because that glob was widened for it.
These tests fail loudly if the asset stops being installed, and they pin the
fallback so a missing mark costs the logo rather than the whole report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from reportlab.lib.pagesizes import A4
from reportlab.platypus import Paragraph, SimpleDocTemplate

from openbiota.pdfreport import (
    BRAND,
    LOGO_FILE,
    REPORT_TITLE,
    TOOL_NAME,
    _styles,
    logo_flowable,
    logo_path,
)


def test_the_brand_and_the_title_are_separate_things() -> None:
    """The title says what the report is; the brand says who made it."""
    assert BRAND == "OpenBiota"
    assert REPORT_TITLE == "Gut Health Test Metagenomic Report"
    # The title must not carry the brand, or the page-1 lockup and the
    # running header would both read "OpenBiota OpenBiota ...".
    assert BRAND not in REPORT_TITLE
    # The tool name is what goes in PDF metadata and the footer, and it does
    # carry the brand.
    assert TOOL_NAME.startswith(BRAND)


def _chrome_pdf(tmp_path: Path, *, pages: int = 3) -> Path:
    """A document with the real page furniture on every page, nothing else."""
    from reportlab.platypus import NextPageTemplate, PageBreak

    from openbiota.pdfreport import _Chrome, _document

    st = _styles()
    out = tmp_path / "chrome.pdf"
    doc = _document(out, _Chrome(sample="T", generated="today", meta_line="Sample T"))
    story: list[Any] = [NextPageTemplate("interior"), Paragraph("page one", st["body"])]
    for n in range(2, pages + 1):
        story += [PageBreak(), Paragraph(f"page {n}", st["body"])]
    doc.build(story)
    return out


def test_every_page_links_the_brand_to_the_website(tmp_path: Path) -> None:
    """The brand is the way home on every page, and the footer always has the URL.

    A printed page 40 carries no URL a reader could type, so the footer sets
    "OpenBiota" as a live link to BRAND_URL on every page, and only the brand
    is - the rest of the footer sentence stays plain text.

    The header differs by page on purpose. On page 1 the wordmark is the
    website. On an interior page the whole running header, wordmark and title
    together, jumps to page 1: two hotspots a few millimetres apart, one
    leaving the document and one moving inside it, is a trap for a reader
    aiming at the title, so the header keeps one job and the footer keeps the
    URL.
    """
    pytest.importorskip("pymupdf")
    import pymupdf

    from openbiota.pdfreport import BRAND_URL, MARGIN, PAGE

    assert BRAND_URL == "https://openbiota.com"
    doc = pymupdf.open(_chrome_pdf(tmp_path))
    assert len(doc) == 3
    for page in doc:
        links = page.get_links()
        home = [
            link for link in links
            if link["kind"] == pymupdf.LINK_URI and link["uri"] == BRAND_URL
        ]
        footer = [link for link in home if link["from"].y0 > 0.9 * PAGE[1]]
        assert len(footer) == 1, f"page {page.number + 1} footer: {home}"
        header = [link for link in links if link["from"].y0 < 0.1 * PAGE[1]]
        assert len(header) == 1, f"page {page.number + 1} header: {header}"
        if page.number == 0:
            # The wordmark, straight to the website.
            assert header[0]["kind"] == pymupdf.LINK_URI
            assert header[0]["uri"] == BRAND_URL
            assert header[0]["from"].width < 0.25 * PAGE[0], header[0]["from"]
        else:
            # The running header, back to page 1.
            assert header[0]["kind"] == pymupdf.LINK_GOTO
            assert header[0]["page"] == 0
        for link in [*header, *footer]:
            rect = link["from"]
            # Every hotspot sits at the left margin, where the brand is set.
            assert rect.x0 - MARGIN < 1.0, rect
            # And none of them spans the footer sentence.
            assert rect.width < 0.5 * PAGE[0], rect


def test_the_wordmark_is_installed_as_package_data() -> None:
    path = logo_path()
    assert path is not None, (
        f"{LOGO_FILE} is not installed; check the package-data glob in "
        "pyproject.toml includes data/*.png"
    )
    assert path.is_file() and path.stat().st_size > 0


def test_the_wordmark_keeps_its_aspect_ratio() -> None:
    """Height is fixed by the caller; width must follow from the asset.

    Asserting both dimensions would silently distort the mark if the asset
    is ever replaced at a different aspect.
    """
    from reportlab.lib.utils import ImageReader

    mark = logo_flowable(height=20.0)
    assert mark is not None
    width_px, height_px = ImageReader(str(logo_path())).getSize()
    assert mark.drawHeight == pytest.approx(20.0)
    assert mark.drawWidth == pytest.approx(20.0 * width_px / height_px, rel=1e-6)


def test_a_missing_wordmark_is_survivable(monkeypatch: pytest.MonkeyPatch) -> None:
    """No asset, no logo — but still a report."""
    import openbiota.pdfreport as pdfreport

    monkeypatch.setattr(pdfreport, "logo_path", lambda: None)
    assert pdfreport.logo_flowable() is None


def test_the_wordmark_sits_in_the_header_not_the_title_block() -> None:
    """The mark is drawn on the canvas, opposite the sample line.

    Putting it in the flow cost vertical space and left the title competing
    with it for the same line. On the canvas it costs the page nothing.
    """
    import inspect

    from openbiota import pdfreport, pdfsummary

    chrome = inspect.getsource(pdfreport._Chrome)
    assert "_wordmark" in chrome
    assert "drawImage" in chrome
    # The title block must be the title and nothing else.
    title_src = inspect.getsource(pdfsummary._head)
    assert "REPORT_TITLE" in title_src
    assert "logo" not in title_src.lower()


def test_page_one_fits_at_full_size_and_is_never_scaled_down() -> None:
    """The layout defect this guards against.

    Page 1 sits in a `KeepInFrame(mode="shrink")`. That scales *uniformly*,
    so vertical overflow pulls the right edge in too and leaves a dead strip
    down the side of the page. Adding the pathogen badge to a row budget set
    without it did exactly that, shrinking the page to 94%. The row budget
    now adapts, and this asserts the result reaches the right margin.
    """
    import io

    from reportlab.pdfgen.canvas import Canvas
    from reportlab.platypus import KeepInFrame

    from openbiota import pdfsummary
    from openbiota.interpret import MetaboliteRow
    from openbiota.pdfreport import CONTENT_WIDTH, FIRST_FRAME_HEIGHT, PAGE, status_for

    def row(name: str, pct: float) -> MetaboliteRow:
        return MetaboliteRow(
            panel="p", metabolite=name, value=1.0, percentile=pct,
            cohort_median=1.0, cohort_p25=0.5, cohort_p75=1.5, detected=True,
            fragments=100,
            status=status_for(percentile=pct, higher_means="adverse"),
            what_it_is="x", made_from="y", higher_means="adverse",
            evidence_strength="moderate", summary="s", evidence_detail="d",
            citation="c", genes=[], confidence="confirmed", implication="i",
            extra_note="", category="neurotransmitter",
        )

    st = _styles()
    story: list[Any] = []
    pdfsummary.summary_page(
        story, st,
        rows=[row(f"Pathway {i}", 5.0 + i) for i in range(12)],
        similarity=None, rpob_profile=None,
        meta={
            "sample": "T", "date": "9 September 2026", "n_panels": 12,
            "read_pairs": "8.0 million", "read_pairs_raw": 8_000_000,
            "cohort_note": None, "n_species": 92, "n_profiles": 0,
        },
        findings=None, age=None, age_meta=None, pathogens=None,
        pathogen_section=7,
    )
    # summary_page emits exactly one flowable: the frame holding the page.
    # Ask it what scale it settled on — that is the quantity that decides
    # whether the right edge reaches the margin.
    assert len(story) == 1
    frame = story[0]
    assert isinstance(frame, KeepInFrame)
    frame.canv = Canvas(io.BytesIO(), pagesize=PAGE)
    frame.wrap(CONTENT_WIDTH, FIRST_FRAME_HEIGHT)
    scale = float(getattr(frame, "_scale", 1.0) or 1.0)
    assert scale <= 1.0 + 1e-6, (
        f"page 1 is being scaled to {1 / scale:.1%} of full size; the row "
        "budget should have dropped a row instead of shrinking the page"
    )


def test_the_header_mark_clears_the_first_frame() -> None:
    """Content must start below the mark, or it draws over it."""
    from openbiota.pdfreport import (
        FIRST_FRAME_HEIGHT,
        FOOTER_H,
        HEADER_LOGO_TOP,
        LOGO_HEIGHT,
        PAGE,
    )

    logo_bottom = PAGE[1] - HEADER_LOGO_TOP - LOGO_HEIGHT
    frame_top = FOOTER_H + FIRST_FRAME_HEIGHT
    assert frame_top <= logo_bottom, (
        f"page-1 content starts at {frame_top:.1f} but the wordmark reaches "
        f"down to {logo_bottom:.1f}; they would overlap"
    )


def test_the_title_appears_in_the_rendered_text(tmp_path: Path) -> None:
    """Guards against the lockup table swallowing the title."""
    pytest.importorskip("pymupdf")
    import pymupdf

    from openbiota import pdfsummary
    from openbiota.interpret import MetaboliteRow
    from openbiota.pdfreport import status_for

    st = _styles()
    story: list[Any] = []
    pdfsummary.summary_page(
        story, st,
        rows=[MetaboliteRow(
            panel="p", metabolite="Butyrate", value=1.0, percentile=88.0,
            cohort_median=1.0, cohort_p25=0.5, cohort_p75=1.5, detected=True,
            fragments=100,
            status=status_for(percentile=88.0, higher_means="favourable"),
            what_it_is="x", made_from="y", higher_means="favourable",
            evidence_strength="moderate", summary="s", evidence_detail="d",
            citation="c", genes=[], confidence="confirmed", implication="i",
            extra_note="", category="neurotransmitter",
        )],
        similarity=None, rpob_profile=None,
        meta={
            "sample": "T", "date": "8 September 2026", "n_panels": 1,
            "read_pairs": "8.0 million", "read_pairs_raw": 8_000_000,
            "cohort_note": None, "n_species": 92, "n_profiles": 0,
        },
        findings=None, age=None, age_meta=None, pathogens=None,
        pathogen_section=7,
    )
    from tests.conftest import section_stubs

    out = tmp_path / "title.pdf"
    SimpleDocTemplate(
        str(out), pagesize=A4, topMargin=36, bottomMargin=36,
        leftMargin=52, rightMargin=52,
    ).build([*story, *section_stubs(st, functions=("p",))])
    text = pymupdf.open(out)[0].get_text()
    assert REPORT_TITLE in " ".join(text.split())


def test_every_page_one_block_renders_at_half_width(tmp_path: Path) -> None:
    """Page 1 puts blocks in two columns, so each must survive half width.

    The blocks were written against the full content width with fixed
    millimetre columns, and at half width one of them computed a *negative*
    column and killed the whole PDF — silently, because the report is
    written last. This renders each block on its own at column width.
    """
    from reportlab.lib.units import mm

    from openbiota import pdfsummary
    from openbiota.interpret import MetaboliteRow
    from openbiota.pdfreport import CONTENT_WIDTH, status_for

    col = (CONTENT_WIDTH - 5 * mm) / 2
    st = _styles()
    rows = [
        MetaboliteRow(
            panel="p", metabolite=f"Pathway {i}", value=1.0, percentile=5.0 + i,
            cohort_median=1.0, cohort_p25=0.5, cohort_p75=1.5, detected=True,
            fragments=100,
            status=status_for(percentile=5.0 + i, higher_means="adverse"),
            what_it_is="x", made_from="y", higher_means="adverse",
            evidence_strength="moderate", summary="s", evidence_detail="d",
            citation="c", genes=[], confidence="confirmed", implication="i",
            extra_note="", category="neurotransmitter",
        )
        for i in range(8)
    ]

    cases = {
        "age": lambda out: pdfsummary._age_block(
            out, st, age=None, age_meta=None, width=col),
        "community": lambda out: pdfsummary._community_block(
            out, st, similarity=None, rpob_profile=None, findings=None, width=col),
        "functions": lambda out: pdfsummary._metabolite_block(
            out, st, rows=rows, width=col, chips=False),
        "patterns": lambda out: pdfsummary._profile_block(
            out, st, similarity=None, width=col, chips=False),
    }
    for name, build in cases.items():
        flow: list[Any] = []
        build(flow)
        assert flow, f"{name} produced nothing at half width"
        for item in flow:
            widths = getattr(item, "_argW", None) or getattr(item, "_colWidths", None)
            if widths:
                assert all(w is None or w >= 0 for w in widths), (
                    f"{name} computed a negative column at width {col:.0f}pt: {widths}"
                )
        from tests.conftest import section_stubs

        doc = SimpleDocTemplate(
            str(tmp_path / f"{name}.pdf"), pagesize=A4,
            topMargin=36, bottomMargin=36, leftMargin=52, rightMargin=52,
        )
        doc.build([*flow, *section_stubs(st, functions=tuple(r.panel for r in rows))])
