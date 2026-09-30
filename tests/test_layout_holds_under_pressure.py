"""AT067 and AT084 — the new tables under the inputs that break layouts.

Long taxonomy names, wide evidence tables and empty optional context are
the three ways a page stops being readable. A cell that clips loses a
species name silently, which is worse than an ugly page: the reader cannot
tell that anything is missing.
"""

from __future__ import annotations

from typing import Any

import pytest

from openbiota import pdfextension as PX
from openbiota import pdfsynbiotic as PS
from openbiota.extension import ecology as ECO

#: Real names get this long. *Anaerostipes caccae* carries three collection
#: designations; the fungal and archaeal catalogues are worse.
LONG_NAME = (
    "Anaerostipes caccae (strain DSM 14662 / CCUG 47493 / JCM 13470 / NCIMB 13811 / L1-92)"
)
LONGER_NAME = (
    "Candidatus Methanomethylophilus alvus Mx1201 subsp. intestinalis biovar "
    "thermotolerans strain DSM 100000 / ATCC BAA-0000 / NCTC 00000"
)


def _render(fn: str, module: Any, view: Any) -> str:
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota.pdfreport import _styles

    story: list[Any] = []
    getattr(module, fn)(story, _styles(), view)
    if not story:
        return ""
    path = f"/tmp/_layout_{fn}.pdf"
    SimpleDocTemplate(path, pagesize=A4).build(story)
    text = "".join(page.extract_text() for page in pypdf.PdfReader(path).pages)
    return " ".join(text.split())


# --- AT067: long names must wrap, never clip ------------------------------


def test_a_long_organism_name_survives_the_redundancy_table() -> None:
    rows = [
        {"organism": LONG_NAME, "value": 60, "is_carrier_evidence": True},
        {"organism": LONGER_NAME, "value": 40, "is_carrier_evidence": True},
    ]
    view = {"redundancy": [ECO.redundancy("panels.x", LONG_NAME, rows).to_json()]}
    text = _render("redundancy_block", PX, view)
    # the label is the function, and it is the long string here
    assert "Anaerostipes caccae" in text
    assert "L1-92" in text, "the tail of a long name was clipped"


def test_a_wide_evidence_table_keeps_its_last_column() -> None:
    """The decoy table's explanation column is the widest thing in the report."""
    view = {"biotransformation": {"negative_controls": {
        "polyphenol.a_very_long_reading_identifier_that_keeps_going": {
            "someEnzymeWithAnUnusuallyLongName": (
                "A long explanation that runs well past the width of any single line and "
                "has to wrap several times before it finishes, ending on this word: "
                "terminus."
            ),
        },
    }}}
    text = _render("_decoys_block", PS, view)
    assert "terminus." in text, "the end of a wrapped explanation was lost"


def test_a_long_route_identifier_does_not_push_the_network_table_over() -> None:
    edge = {
        "from": "x", "from_label": "a substrate with a long descriptive label",
        "to": "y", "to_label": "a product with an equally long descriptive label",
        "route_id": "fermentation.some.very.long.route.identifier.that.keeps.going",
        "support": "supported", "style": "solid",
    }
    view = {"fermentation": {"network": {"edges": [edge], "limitations": ["Not flows."]}}}
    text = _render("_network_block", PS, view)
    # A long identifier may be broken across lines to fit its column. That is
    # wrapping, not clipping: what matters is that no character is dropped,
    # so the comparison ignores the break points the layout chose.
    assert "that.keeps.going" in text.replace(" ", "")
    assert "equally long descriptive label" in text


# --- AT084: the report with no optional context ---------------------------


@pytest.mark.parametrize(
    ("module", "fn"),
    [
        (PX, "redundancy_block"),
        (PS, "_network_block"),
        (PS, "_shared_genes_block"),
        (PS, "_paired_readings_block"),
        (PS, "_three_concepts_block"),
        (PS, "_host_steps_block"),
        (PS, "_decoys_block"),
        (PS, "_coverage_block"),
        (PS, "_evidence_seeds_block"),
    ],
)
def test_every_new_block_is_silent_rather_than_broken_without_its_input(
    module: Any, fn: str,
) -> None:
    """A missing optional view must produce no section, not an empty frame.

    An empty heading with nothing under it reads as a failed measurement.
    """
    assert _render(fn, module, {}) == ""
    assert _render(fn, module, None) == ""


@pytest.mark.parametrize(
    ("module", "fn", "view"),
    [
        (PX, "redundancy_block", {"redundancy": []}),
        (PS, "_network_block", {"fermentation": {"network": {"edges": []}}}),
        (PS, "_shared_genes_block", {"substrates": {"shared_genes": {"shared": {}}}}),
        (PS, "_decoys_block", {"biotransformation": {"negative_controls": {}}}),
    ],
)
def test_an_empty_collection_is_also_silent(module: Any, fn: str, view: Any) -> None:
    assert _render(fn, module, view) == ""
