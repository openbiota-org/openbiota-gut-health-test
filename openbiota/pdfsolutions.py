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
from openbiota.pdflinks import Paragraph, section_heading
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

KIND_WORDS: Final[dict[str, str]] = {
    "pathogen": "PATHOGEN SCREEN",
    "organism": "ORGANISM",
    "concern": "WHOLE-GUT FINDING",
    "panel": "METABOLITE PATHWAY",
    "pattern": "DISEASE PATTERN",
}


def _hex(c: colors.Color) -> str:
    return c.hexval().replace("0x", "#")


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
                f"<font size='6.6'><b>{name}</b>"
                f"<super><font size='5.4' color='{_hex(INK_FAINT)}'>{reach.marker_text}</font></super>"
                "</font>", st["cell"]),
            Paragraph(
                f"<font size='6.2' color='{_hex(INK_SOFT)}'>{option.kind_words}</font>"
                + (f"<br/><font size='5.8' color='{_hex(CORAL)}'>prescriber decides</font>"
                   if option.needs_prescriber else ""),
                st["cell"]),
            Paragraph(
                f"<font size='6.2' color='{_hex(tone)}'>{reach.setting_words}</font>"
                + (f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>{reach.delivery_words}</font>"
                   if reach.marker_text else ""),
                st["cell"]),
            Paragraph(
                (f"<font size='6.2'>{option.effect[:220]}</font>" if option.effect else "")
                + ("<br/>" if option.effect and option.boundary else "")
                + (f"<font size='5.8' color='{_hex(INK_FAINT)}'>{option.boundary[:170]}</font>"
                   if option.boundary else "")
                + ("" if (option.effect or option.boundary)
                   else f"<font size='6.2' color='{_hex(INK_FAINT)}'>{EM_DASH}</font>"),
                st["cell"]),
        ])
    return rows


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

    head = [
        Paragraph(
            f"<font size='5.8' color='{_hex(colour)}'><b>{kind_words}</b></font>", st["small"]),
        Paragraph(
            f"<font size='9.2' color='{_hex(INK)}'><b>{position}. {problem.title}</b></font>", st["h3"]),
    ]
    if problem.detail:
        head.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>{problem.detail} \u00b7 "
            f"see section {problem.section}</font>", st["small"]))

    if problem.nothing_known:
        head.append(Spacer(1, 1.2 * mm))
        head.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>Nothing on record.</b> No study in any "
            f"setting \u2014 trial, animal or laboratory \u2014 has been recorded here for how to "
            f"{want_words}. That is the state of the research, not a judgement that this does not "
            "matter, and not a reason to do nothing: the sections above explain the finding and a "
            "clinician can act on it.</font>", st["small"]))
        story.append(KeepTogether(head))
        story.append(Spacer(1, 2.5 * mm))
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
    story.append(Spacer(1, 3 * mm))


def _contents(
    story: list[Any], st: Mapping[str, ParagraphStyle], problems: Sequence[sol.Problem]
) -> None:
    """Everything the section covers, before it covers it."""
    head = ("WHAT WAS FOUND", "WHERE", "OPTIONS ON RECORD")
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
        rows.append([
            Paragraph(f"<font size='6.6'><b>{i}.</b> {problem.title}</font>", st["cell"]),
            Paragraph(f"<font size='6' color='{_hex(colour)}'>"
                      f"{KIND_WORDS.get(problem.kind, '')} \u00b7 section {problem.section}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{count}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.56, 0.28, 0.16)]))


def section(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    section_no: int,
    problems: Sequence[sol.Problem],
) -> None:
    """The whole of the actions section."""
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
