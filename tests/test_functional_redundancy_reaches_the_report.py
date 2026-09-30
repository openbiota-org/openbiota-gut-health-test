"""§6.4 — how many organisms hold each function, and whether it is shown.

The module that computes this was written, exported and never called, so
the reading existed in the code and in no report. These tests are mostly
about the second half of that sentence: the number has to reach the page,
and it has to mean what the scale says it means.
"""

from __future__ import annotations

from typing import Any

import pytest

from openbiota.extension import ecology as ECO


def _reading(label: str, carriers: list[tuple[str, float]], **kw: Any) -> dict[str, Any]:
    rows = [{"organism": n, "value": v, "is_carrier_evidence": True} for n, v in carriers]
    return ECO.redundancy(f"panels.{label}", label, rows, **kw).to_json()


# --- the measure -----------------------------------------------------------


def test_an_even_split_counts_every_carrier() -> None:
    reading = _reading("butyrate", [("a", 250), ("b", 250), ("c", 250), ("d", 250)])
    assert reading["effective_carriers"] == pytest.approx(4.0, abs=1e-6)


def test_a_lopsided_split_counts_barely_more_than_one() -> None:
    """Two carriers at 99/1 are not two carriers in any sense that matters."""
    reading = _reading("p-cresol", [("dominant", 990), ("trace", 10)])
    assert reading["n_resolved_carriers"] == 2
    assert reading["effective_carriers"] < 1.1


def test_a_single_carrier_is_exactly_one() -> None:
    assert _reading("equol", [("solo", 300)])["effective_carriers"] == pytest.approx(1.0)


def test_an_organism_not_in_the_inventory_is_context_not_a_carrier() -> None:
    """A read matching some organism's copy of a gene does not put it in you."""
    rows = [
        {"organism": "in you", "value": 100, "is_carrier_evidence": True},
        {"organism": "reference only", "value": 100, "is_carrier_evidence": False},
    ]
    reading = ECO.redundancy("panels.x", "x", rows).to_json()
    assert reading["n_resolved_carriers"] == 1
    assert reading["effective_carriers"] == pytest.approx(1.0)
    # the unresolved one still counts against coverage rather than vanishing
    assert reading["carrier_coverage"] == pytest.approx(0.5)


def test_no_carriers_yields_no_number_rather_than_zero() -> None:
    reading = ECO.redundancy("panels.x", "x", []).to_json()
    assert reading["effective_carriers"] is None


# --- the scale used to draw it --------------------------------------------


def test_the_bar_position_is_bounded_and_monotone() -> None:
    from openbiota.pdfextension import _redundancy_position

    assert _redundancy_position(0.0) == 0.0
    assert _redundancy_position(2.0) == pytest.approx(50.0)
    assert _redundancy_position(4.0) == 100.0
    # unbounded above, but the bar is not
    assert _redundancy_position(40.0) == 100.0
    assert _redundancy_position(1.0) < _redundancy_position(3.0)


def test_one_carrier_and_four_carriers_do_not_share_a_colour() -> None:
    from openbiota.pdfextension import _redundancy_tone

    assert _redundancy_tone(1.0) is not _redundancy_tone(4.0)


# --- reaching the page -----------------------------------------------------


def _render(view: dict[str, Any]) -> str:
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota import pdfextension as PX
    from openbiota.pdfreport import _styles

    story: list[Any] = []
    PX.redundancy_block(story, _styles(), view)
    if not story:
        return ""
    path = "/tmp/_redundancy_render.pdf"
    SimpleDocTemplate(path, pagesize=A4).build(story)
    text = "".join(page.extract_text() for page in pypdf.PdfReader(path).pages)
    # A phrase that wraps across a line is still on the page; collapse the
    # layout's newlines so the assertions test the prose, not the column width.
    return " ".join(text.split())


def test_every_function_with_a_count_reaches_the_page() -> None:
    view = {"redundancy": [
        _reading("Butyrate route", [("a", 250), ("b", 250), ("c", 250), ("d", 250)]),
        _reading("Equol production", [("solo", 300)]),
        _reading("Bile salt hydrolase", [("a", 500), ("b", 500)]),
    ]}
    text = _render(view)
    for label in ("Butyrate route", "Equol production", "Bile salt hydrolase"):
        assert label in text, f"{label} never reached the page"
    assert "4.0" in text and "1.0" in text and "2.0" in text


def test_the_fragile_functions_are_listed_first() -> None:
    """A reader scanning the top of the table should meet the risks there."""
    view = {"redundancy": [
        _reading("spread wide", [("a", 250), ("b", 250), ("c", 250), ("d", 250)]),
        _reading("held by one", [("solo", 300)]),
    ]}
    text = _render(view)
    assert text.index("held by one") < text.index("spread wide")


def test_a_function_without_a_count_is_left_out_rather_than_drawn_at_zero() -> None:
    """No carriers resolved is not the same as one carrier, or none."""
    view = {"redundancy": [
        ECO.redundancy("panels.x", "nothing resolved", []).to_json(),
    ]}
    assert _render(view) == ""


def test_a_weakly_confirmed_count_is_marked_as_a_floor() -> None:
    """Most of a function's signal can sit on organisms we cannot place.

    Butyrate in a real sample attributes reads to eight reference organisms
    and confirms one of them in the inventory, so "one carrier" is the
    fewest it could be rather than the number there are. A reader told only
    "1 carrier, 100% share" would take that for demonstrated fragility.
    """
    rows = [
        {"organism": "confirmed", "value": 20, "is_carrier_evidence": True},
        {"organism": "unplaced", "value": 80, "is_carrier_evidence": False},
    ]
    view = {"redundancy": [ECO.redundancy("panels.b", "Butyrate", rows).to_json()]}
    text = _render(view)
    assert "SIGNAL CONFIRMED" in text
    assert "20%" in text, "the confirmed share has to be on the page"
    assert "not as the number there are" in text


def test_a_fully_confirmed_count_carries_no_caveat() -> None:
    rows = [
        {"organism": "a", "value": 50, "is_carrier_evidence": True},
        {"organism": "b", "value": 50, "is_carrier_evidence": True},
    ]
    view = {"redundancy": [ECO.redundancy("panels.b", "Butyrate", rows).to_json()]}
    text = _render(view)
    assert "100%" in text
    assert "not as the number there are" not in text


def test_the_block_says_what_it_is_not() -> None:
    view = {"redundancy": [_reading("x", [("a", 1), ("b", 1)])]}
    text = _render(view)
    assert "not a guarantee of ecological resilience" in text
