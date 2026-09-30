"""Page-1 labels must fit on one line, measured against the real column.

Page 1 is a grid, and a grid only reads as one if every row is the same
height. These tests measure each short label with the font page 1 actually
uses, against the width the name column actually has, so a new profile or
panel whose short form is too long fails here rather than wrapping in the
report of whichever sample happens to score it highly.
"""

from __future__ import annotations

from pathlib import Path

from reportlab.pdfbase.pdfmetrics import stringWidth

from openbiota.pdfreport import CONTENT_WIDTH, _styles
from openbiota.pdfsummary import COLUMN_GUTTER, name_column_width
from openbiota.profiles import load_profile_set
from openbiota.shortnames import (
    _METABOLITE_SHORT,
    _PROFILE_SHORT,
    self_test,
    short_metabolite,
    short_profile_label,
)

ROOT = Path(__file__).resolve().parents[1]


def _fits(text: str, suffix_width: float = 0.0) -> bool:
    st = _styles()
    col = (CONTENT_WIDTH - COLUMN_GUTTER) / 2
    room = name_column_width(col, chips=False) - 2.0  # keep a little air
    return stringWidth(text, "Helvetica-Bold", st["cell"].fontSize) + suffix_width <= room


def test_module_self_test() -> None:
    assert self_test() > 0


def test_every_profile_in_the_library_has_a_short_form_that_fits() -> None:
    """Coverage and width in one pass over the real profile set."""
    profiles = load_profile_set(ROOT / "profiles").profiles
    assert profiles, "no profiles loaded"
    missing = [p.name for p in profiles if p.name not in _PROFILE_SHORT]
    assert not missing, f"profiles with no explicit short label: {missing}"
    for p in profiles:
        short = short_profile_label(p.name, p.label)
        assert _fits(short), f"{p.name}: {short!r} does not fit the page-1 name column"
        # The short form has to be shorter, or it is not doing its job.
        assert len(short) < len(p.label) or short == p.label


def test_every_metabolite_short_form_fits_with_its_direction_arrow() -> None:
    arrow = stringWidth(" \u25b2", "Helvetica", 6)
    for long, short in _METABOLITE_SHORT.items():
        assert short_metabolite(long) == short
        assert _fits(short, arrow), f"{long!r} -> {short!r} does not fit with its arrow"


def test_fallback_trims_the_shared_boilerplate() -> None:
    assert short_profile_label("new", "Some condition — gut pattern") == "Some condition"
    assert short_profile_label("new", "Some condition — gut pattern (opt-in)") == "Some condition"
    assert short_metabolite("Thing (with a long parenthetical)") == "Thing"


def test_page_one_rows_use_the_short_forms() -> None:
    """The dashboard must go through the short-name layer, not the raw label."""
    import inspect

    from openbiota import pdfsummary

    assert "short_profile_label(" in inspect.getsource(pdfsummary._profile_block)
    assert "short_metabolite(" in inspect.getsource(pdfsummary._metabolite_block)
