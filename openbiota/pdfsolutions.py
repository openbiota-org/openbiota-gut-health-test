"""The actions section: every problem this report found, and what might be done.

One structure, repeated. For each finding, in the order of how much it
matters:

    THE PROBLEM     what was found, with the number that makes it concrete
    WHAT MIGHT HELP the exact food, compound, supplement, probiotic, herb,
                    medication or procedure, what kind of thing it is, where
                    the finding comes from, and what was actually seen
    WHAT WOULD MAKE IT WORSE   the same, in the direction you do not want
    ...and when nothing is on record, the words "nothing on record", which
    is a finding too.

Everything that has been shown is here, whatever setting it came from. A
laboratory result and a randomised trial sit in the same table with the
setting named in the row, a small marker, and one line of fine print at the
foot of the section. The introduction carries no caveats: a reader looking
for what to do should find what to do.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, KeepTogether, Spacer, Table, TableStyle

from openbiota import solutions as sol
from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    linked_cell,
    register_solutions,
    section_heading,
    solution_dest,
)
from openbiota.pdfreport import (
    ACCENT,
    ACCENT_DARK,
    AMBER,
    CONTENT_WIDTH,
    CORAL,
    GREEN,
    INK,
    INK_FAINT,
    INK_SOFT,
    PANEL_BG,
    RULE,
    SLATE,
)

#: How many options to print per problem. Enough to choose from, few enough
#: that the page stays readable; the rest are counted, never silently lost.
PER_PROBLEM: Final = 10
#: And how many of the wrong-direction ones.
PER_PROBLEM_WORSE: Final = 4
#: And how many of the studied-but-undirected ones.
PER_PROBLEM_UNCLEAR: Final = 6

#: How much of an option's limits to print. Long enough for the whole
#: caveat on the rows that carry a real one.
#: And how much of what was actually seen.
EFFECT_CHARS: Final = 340

BOUNDARY_CHARS: Final = 560

#: Used inside f-strings, where a backslash escape is not allowed.
EM_DASH: Final = "\u2014"

#: The colour each problem kind wears, so the eye can sort them.
KIND_COLOUR: Final[dict[str, colors.Color]] = {
    "pathogen": CORAL,
    "organism": AMBER,
    "concern": ACCENT_DARK,
    "panel": ACCENT,
    "pattern": SLATE,
}

#: Fallback destination when a finding has no page of its own.
_SECTION_ANCHOR: Final[dict[str, str]] = {
    "pathogen": "pathogens_detail", "organism": "organisms_detail",
    "concern": "biofilm_detail", "panel": "functions_detail",
    "pattern": "patterns_detail",
}

KIND_WORDS: Final[dict[str, str]] = {
    "pathogen": "PATHOGEN SCREEN",
    "organism": "ORGANISM",
    "concern": "WHOLE-GUT FINDING",
    "panel": "METABOLITE PATHWAY",
    "pattern": "DISEASE PATTERN",
}


def _hex(c: colors.Color) -> str:
    return c.hexval().replace("0x", "#")


def _contents_style() -> TableStyle:
    return TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), PANEL_BG),
        ("LINEBELOW", (0, 0), (-1, 0), 0.4, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor("#F0F3F6")),
        ("LEFTPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 1.1 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.1 * mm),
    ])


def _xml(text: object) -> str:
    """Text safe to drop into a reportlab paragraph.

    Citations are full of ampersands and apostrophes, and an unescaped one
    takes the whole report down with a parser error.
    """
    return (str(text or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _table(rows: list[list[Any]], widths: list[float], *, header_bg: colors.Color = PANEL_BG) -> Table:
    table = Table(rows, colWidths=widths, hAlign="LEFT", repeatRows=1)
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), header_bg),
        ("LINEBELOW", (0, 0), (-1, 0), 0.4, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor("#F0F3F6")),
        ("LEFTPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 1.1 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.1 * mm),
    ]))
    return table


def _option_rows(
    st: Mapping[str, ParagraphStyle], options: Sequence[sol.Option], *, limit: int
) -> list[list[Any]]:
    head = ("WHAT", "KIND", "WHERE THIS COMES FROM", "WHAT WAS SEEN")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    for option in options[:limit]:
        reach = option.reach
        tone = GREEN if option.setting in {"human_trial", "human_observational"} else INK_SOFT
        name = option.label if len(option.label) <= 96 else option.label[:95] + "\u2026"
        rows.append([
            Paragraph(
                f"<font size='6.6'><b>{_xml(name)}</b>"
                f"<super><font size='5.4' color='{_hex(INK_FAINT)}'>{reach.marker_text}</font></super>"
                "</font>", st["cell"]),
            Paragraph(
                f"<font size='6.2' color='{_hex(INK_SOFT)}'>{_xml(option.kind_words)}</font>"
                + (f"<br/><font size='5.8' color='{_hex(CORAL)}'>prescriber decides</font>"
                   if option.needs_prescriber else ""),
                st["cell"]),
            Paragraph(
                f"<font size='6.2' color='{_hex(tone)}'>{reach.setting_words}</font>"
                + (f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>{reach.delivery_words}</font>"
                   if reach.marker_text else "")
                # The paper itself. Eight hundred of these were sitting in the
                # data unprinted: a reader who wants to check a claim, or take
                # it to a clinician, needs the reference, not a promise that
                # one exists.
                + (f'<br/><a href="{_xml(option.source)}">'
                   f"<font size='5.6' color='{_hex(ACCENT_DARK)}'>"
                   f"{_xml(option.citation) or 'source'}</font></a>"
                   if option.source else
                   (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>"
                    f"{_xml(option.citation)}</font>" if option.citation else "")),
                st["cell"]),
            Paragraph(
                (f"<font size='6.2'>{_xml(option.effect[:EFFECT_CHARS])}</font>" if option.effect else "")
                + ("<br/>" if option.effect and option.boundary else "")
                # The boundary is where the limits live - what a material is
                # not, who it was studied in, what it has been shown to
                # harm. Clipping it at a tweet's length cut the sentence a
                # reader most needs off the end of the row, so it gets the
                # room it takes.
                + (f"<font size='5.8' color='{_hex(INK_FAINT)}'>{_xml(option.boundary[:BOUNDARY_CHARS])}</font>"
                   if option.boundary else "")
                + ("" if (option.effect or option.boundary)
                   else f"<font size='6.2' color='{_hex(INK_FAINT)}'>{EM_DASH}</font>"),
                st["cell"]),
        ])
    return rows


def _nav_line(st: Mapping[str, ParagraphStyle], problem: sol.Problem) -> Any:
    """Back to the list, and across to the reading, on one line."""
    parts = [f"<a href='#{CONTENTS_DEST}'><font color='{_hex(ACCENT_DARK)}'>"
             f"\u2039 back to the list</font></a>"]
    if problem.detail_dest:
        parts.append(f"<a href='#{problem.detail_dest}'><font color='{_hex(ACCENT_DARK)}'>"
                     f"the reading this is about \u203a</font></a>")
    joined = NAV_SEPARATOR.join(parts)
    return Paragraph(f"<font size='6.6'>{joined}</font>", st["small"])


def _problem_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], problem: sol.Problem, *, position: int
) -> None:
    colour = KIND_COLOUR.get(problem.kind, SLATE)
    kind_words = KIND_WORDS.get(problem.kind, problem.kind.upper())
    # A disease-pattern resemblance is not a quantity to raise or lower, so
    # it is not offered as one: what is on record for it is what has been
    # tried in people who have the condition.
    if problem.want == sol.WANT_DOWN:
        want_words, heading = "bring it down", "WHAT MIGHT BRING IT DOWN"
    elif problem.want == sol.WANT_UP:
        want_words, heading = "bring it up", "WHAT MIGHT BRING IT UP"
    else:
        want_words, heading = "change it", "WHAT HAS BEEN TRIED IN THIS CONDITION"

    head: list[Any] = [
        Anchor(solution_dest(problem.key), above=5 * mm),
        Paragraph(
            f"<font size='5.8' color='{_hex(colour)}'><b>{kind_words}</b></font>", st["small"]),
        Paragraph(
            f"<font size='9.2' color='{_hex(INK)}'><b>{position}. {problem.title}</b></font>", st["h3"]),
    ]
    # Two ways out, in the order a reader wants them: back to the list they
    # came from, and across to the reading this is about.
    head.insert(1, _nav_line(st, problem))
    if problem.detail:
        head.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>{_xml(problem.detail)}</font>", st["small"]))

    if problem.nothing_known and not problem.unclear:
        head.append(Spacer(1, 1.2 * mm))
        head.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>Nothing on record.</b> No study in any "
            f"setting \u2014 trial, animal or laboratory \u2014 has been recorded here for how to "
            f"{want_words}. That is the state of the research, not a judgement that this does not "
            "matter, and not a reason to do nothing: the sections above explain the finding and a "
            "clinician can act on it.</font>", st["small"]))
        for note in problem.notes:
            head.append(Paragraph(
                f"<font size='6.2' color='{_hex(INK_FAINT)}'>{note}</font>", st["small"]))
        story.append(KeepTogether(head))
        story.append(Spacer(1, 2.5 * mm))
        return

    if problem.nothing_known:
        head.append(Spacer(1, 1.2 * mm))
        head.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>No study on record says which way to "
            f"{want_words}. What has been studied on it is below, with the direction its source "
            "actually reported \u2014 which was none.</font>", st["small"]))
        for note in problem.notes:
            head.append(Paragraph(
                f"<font size='6.2' color='{_hex(INK_FAINT)}'>{note}</font>", st["small"]))
        story.append(KeepTogether(head))
        story.append(Spacer(1, 1.2 * mm))
        story.append(KeepTogether([
            Paragraph(
                f"<font size='6.4' color='{_hex(SLATE)}'><b>ALSO STUDIED ON THIS, DIRECTION NOT "
                f"STATED</b></font>  <font size='5.8' color='{_hex(INK_FAINT)}'>"
                f"{len(problem.unclear)}</font>", st["small"]),
            _table(_option_rows(st, problem.unclear,
                                limit=PER_PROBLEM_UNCLEAR),
                   [CONTENT_WIDTH * w for w in (0.30, 0.16, 0.18, 0.36)]),
        ]))
        story.append(Spacer(1, 3 * mm))
        return

    head.append(Spacer(1, 1.2 * mm))
    shown = min(len(problem.options), PER_PROBLEM)
    more = len(problem.options) - shown
    count = (f"{shown} of {len(problem.options)} on record" if more else str(len(problem.options)))
    head.append(Paragraph(
        f"<font size='6.4' color='{_hex(GREEN)}'><b>{heading}</b></font>"
        f"  <font size='5.8' color='{_hex(INK_FAINT)}'>{count}</font>", st["small"]))
    widths = [CONTENT_WIDTH * w for w in (0.30, 0.16, 0.18, 0.36)]
    head.append(_table(_option_rows(st, problem.options, limit=PER_PROBLEM), widths))
    story.append(KeepTogether(head))

    if problem.worsens:
        tail: list[Any] = [
            Spacer(1, 1.2 * mm),
            Paragraph(
                f"<font size='6.4' color='{_hex(CORAL)}'><b>SHOWN TO PUSH IT THE OTHER WAY</b></font>"
                f"  <font size='5.8' color='{_hex(INK_FAINT)}'>{len(problem.worsens)}</font>",
                st["small"]),
            _table(_option_rows(st, problem.worsens, limit=PER_PROBLEM_WORSE), widths),
        ]
        story.append(KeepTogether(tail))

    # Everything else that has been studied on this finding. These are real
    # materials with real papers whose source never said which way the
    # reading moved, so they are not claimed to help - and they are not
    # hidden either, which is what used to happen to a third of the record.
    if problem.unclear:
        more: list[Any] = [
            Spacer(1, 1.2 * mm),
            Paragraph(
                f"<font size='6.4' color='{_hex(SLATE)}'><b>ALSO STUDIED ON THIS, DIRECTION NOT "
                f"STATED</b></font>  <font size='5.8' color='{_hex(INK_FAINT)}'>"
                f"{len(problem.unclear)}</font>", st["small"]),
            _table(_option_rows(st, problem.unclear, limit=PER_PROBLEM_UNCLEAR), widths),
        ]
        story.append(KeepTogether(more))
    story.append(Spacer(1, 3 * mm))


#: The destination of the section's own list, so every block can go back to
#: the row the reader clicked. Landing on a block with no way back but
#: "the reading this is about" was a dead end.
CONTENTS_DEST: Final = "sol-contents"

#: Between the two ways out of a block.
NAV_SEPARATOR: Final = "  &nbsp;&nbsp;\u00b7&nbsp;&nbsp;  "


def _contents(
    story: list[Any], st: Mapping[str, ParagraphStyle], problems: Sequence[sol.Problem]
) -> None:
    """Everything the section covers, before it covers it."""
    story.append(Anchor(CONTENTS_DEST, above=4 * mm))
    head = ("WHAT WAS FOUND", "WHERE IT IS EXPLAINED", "OPTIONS ON RECORD")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    for i, problem in enumerate(problems, start=1):
        colour = KIND_COLOUR.get(problem.kind, SLATE)
        count = (
            f"<font color='{_hex(INK_FAINT)}'>nothing on record</font>"
            if problem.nothing_known else
            f"<b>{len(problem.options)}</b>"
            + (f" <font color='{_hex(CORAL)}'>\u00b7 {len(problem.worsens)} to avoid</font>"
               if problem.worsens else "")
        )
        # The whole row goes to the block that answers this finding, further
        # down this section. The middle cell is the one exception: it goes to
        # the page that explains the reading, and it goes to the reading
        # itself rather than to the top of the section it lives in.
        row_link = RowLink(href=solution_dest(problem.key))
        where = (
            f"<font size='6' color='{_hex(colour)}'>{KIND_WORDS.get(problem.kind, '')}</font>"
        )
        where_cell: Any
        if problem.detail_dest:
            where_cell = Paragraph(
                where + f"<br/><a href='#{problem.detail_dest}'>"
                f"<font size='6' color='{_hex(ACCENT_DARK)}'>this reading, in full \u203a</font></a>",
                st["cell"])
        else:
            where_cell = Paragraph(
                where + f"<br/><a href='#sec-{_SECTION_ANCHOR.get(problem.kind, 'actions')}'>"
                f"<font size='6' color='{_hex(ACCENT_DARK)}'>section {problem.section}</font></a>",
                st["cell"])
        rows.append([
            # Coloured and arrowed, so a row reads as the link it is.
            linked_cell(row_link, Paragraph(
                f"<font size='6.6'><b>{i}.</b> "
                f"<a href='#{solution_dest(problem.key)}'><font color='{_hex(ACCENT_DARK)}'>"
                f"{_xml(problem.title)}</font></a>"
                f" <font size='6' color='{_hex(ACCENT)}'>\u203a</font></font>", st["cell"])),
            where_cell,
            Paragraph(f"<font size='6.4'>{count}</font>", st["cell"]),
        ])
    table = LinkedTable(rows, colWidths=[CONTENT_WIDTH * w for w in (0.56, 0.28, 0.16)],
                        hAlign="LEFT", repeatRows=1)
    table.setStyle(_contents_style())
    story.append(table)


def section(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    section_no: int,
    problems: Sequence[sol.Problem],
) -> None:
    """The whole of the actions section."""
    # Tell the rest of the report which findings have an answer here, so a
    # detail page offers the way through only when there is one.
    register_solutions(p.key for p in problems)
    story.extend(section_heading(section_no, "What you can do about it", st["h1"]))

    if not problems:
        story.append(Paragraph(
            "This report found nothing that needs changing: no organism outside its usual range, "
            "no metabolite pathway away from the middle of the reference, and nothing on the "
            "pathogen screen. There is accordingly nothing in this section.", st["body"]))
        return

    n_with = sum(1 for p in problems if not p.nothing_known)
    n_options = sum(len(p.options) for p in problems)
    story.append(Paragraph(
        f"This report found <b>{len(problems)}</b> things worth doing something about. "
        f"For <b>{n_with}</b> of them there is published research on changing them \u2014 "
        f"<b>{n_options}</b> findings in all, each one named below with the exact food, "
        "compound, supplement, probiotic, herb, medication or procedure that was studied, "
        "where the finding comes from, and what was actually seen. Each problem is set out the "
        "same way: what was found, what might move it the right way, and what has been shown to "
        "move it the wrong way.",
        st["body"]))
    story.append(Spacer(1, 2 * mm))
    _contents(story, st, problems)
    story.append(Spacer(1, 4 * mm))

    for i, problem in enumerate(problems, start=1):
        story.append(CondPageBreak(55 * mm))
        _problem_block(story, st, problem, position=i)

    notes = sol.footnotes(problems)
    if notes:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            "<font size='6' color='" + _hex(INK_FAINT) + "'>"
            + "  ".join(notes)
            + "</font>", st["small"]))
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font size='6' color='{_hex(INK_FAINT)}'>Every row is a description of published "
        "research, not a recommendation, and not a claim that it is safe or right for you. "
        "Nothing here is a prescription; anything marked as a prescriber's decision is one. "
        "Talk to a qualified clinician before acting on any of it.</font>", st["small"]))


__all__ = ["section"]
