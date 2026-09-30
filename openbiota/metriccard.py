"""The one shape every metric is reported in.

READ THIS BEFORE ADDING A METRIC TO THE REPORT.

A reader should not have to learn a new layout for each measurement. For a
while they did: the panels had a bar, an ordinal and a status chip, and the
functional readings underneath them had a column of gene symbols and the
words "class-level only" instead. Two renderings of the same kind of thing,
one of which said *pta, acsB, cooS* to somebody who wants to know whether
their acetate production is high.

So there is one row and one card, and both are built here.

THE ROW. Every metric, wherever it appears, is drawn as:

    label (with its direction arrow) | scale bar | ordinal | status chip

and one plain sentence under it saying what the reading means. Nothing
else goes in the row: gene symbols, evidence grades and specificity belong
on the card, which is one click away.

THE CARD. Every metric's detail page carries these sections, in this
order, with these headings:

    1. What is driving this      - the organisms and genes behind the
                                   number, in this sample, named
    2. What your reading means   - the reading in plain words
    3. How can I improve this    - which direction is healthier and what
                                   is known about moving it, or plainly
                                   that nothing is established
    4. The research behind it    - the evidence, with citations

`CARD_SECTIONS` is that list, and `tests/test_every_metric_looks_the_same.py`
asserts every card carries all of it. A metric that cannot fill a section
says why in that section; it does not omit the heading, because a missing
heading reads as an oversight and a stated gap reads as a finding.

Do not add a metric with a bespoke layout. If a metric genuinely needs
something this does not cover, add it to every metric or add it as a
subsection inside one of these four.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib.styles import ParagraphStyle

#: The headings every metric card carries, in order. The report is checked
#: against this list; changing it changes the contract for every metric.
CARD_SECTIONS: Final[tuple[str, ...]] = (
    "What is driving this",
    "What your reading means",
    "How can I improve this",
    "The research behind it",
)

#: What "better" means for a direction, in words a reader can act on.
_DIRECTION_ADVICE: Final[Mapping[str, str]] = {
    "adverse": "Lower is the more favourable direction for this reading.",
    "favourable": "Higher is the more favourable direction for this reading.",
    "unclear": (
        "Neither direction is established as better for this reading, so there is no "
        "target to move it towards."
    ),
}


def direction_sentence(higher_means: str) -> str:
    """Which way is healthier, said plainly."""
    return _DIRECTION_ADVICE.get(str(higher_means), _DIRECTION_ADVICE["unclear"])


def improvement_block(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    higher_means: str,
    organisms: Sequence[Mapping[str, Any]] = (),
    levers: Sequence[str] = (),
    nothing_established: str = "",
) -> None:
    """Section 3 of the card: how to move the reading, honestly.

    Most readings have no lever anybody has demonstrated. That is the
    normal case and it is said in as many words, because a reader who is
    offered a plausible-sounding action for a reading nothing can move is
    being misled in the direction of spending money.

    Where organisms behind the reading are present in this sample they are
    named, because those are the only things in the gut this reading is
    actually about.
    """
    from openbiota.pdfparts import hex_of  # noqa: PLC0415
    from openbiota.pdfreport import AMBER, INK_SOFT  # noqa: PLC0415

    story.append(_para(st, f"<b>{CARD_SECTIONS[2]}?</b>", "h3"))
    story.append(_para(st, direction_sentence(higher_means), "body_ink"))

    if organisms:
        named = "; ".join(
            f"<i>{o.get('organism')}</i>"
            + (f" ({o['sample_percent']:.2f}% of your community)"
               if o.get("sample_percent") is not None else "")
            for o in organisms[:4]
        )
        story.append(_para(
            st,
            f"<font size='6.8'>This reading follows the organisms carrying it, and the "
            f"ones found in you are {named}. Anything that changes this reading has to "
            f"change them.</font>", "small"))

    if levers:
        for lever in levers:
            story.append(_para(st, f"<font size='6.8'>\u25aa {lever}</font>", "small"))
        return

    story.append(_para(
        st,
        f"<font size='6.8' color='{hex_of(AMBER)}'>No study on record shows that a diet, "
        "a supplement or a probiotic moves this particular reading.</font>"
        f"<font size='6.8' color='{hex_of(INK_SOFT)}'> "
        + (nothing_established or
           "What is established for the organisms behind it is in the research below. "
           "Treating an unproven lever as a proven one is the failure this report is "
           "built to avoid, so the gap is printed rather than filled.")
        + "</font>", "small"))


def _para(st: Mapping[str, ParagraphStyle], text: str, style: str) -> Any:
    from openbiota.pdflinks import Paragraph  # noqa: PLC0415

    return Paragraph(text, st[style])


__all__ = ["CARD_SECTIONS", "direction_sentence", "improvement_block"]
