"""The biofilm section: two headings, four cards, and nothing that cancels.

Biofilms are the part of a gut report most likely to be sold badly. The
temptation is a single "biofilm score" with a colour, which is worthless:
protective and harmful matrix communities coexist, and a number that
averages them describes neither.

So the layout here is structural rather than decorative. Two headings,
concerning and protective, side by side. Under each, at most two numeric
cards. The two columns are the same width and the two scales are the same
0-100 reference percentile, so a reader can compare them by eye - and
because they are separate columns, there is nowhere for the arithmetic to
quietly subtract one from the other.

Three rendering rules earn their keep:

* A null percentile is rendered as text explaining why, never as a zero
  and never in green. An empty H-M card says "no validated panel yet",
  which is a different statement from "no harmful mechanisms found".
* Every colour carries a text equivalent. The band words - lower, middle,
  upper reference range - are printed next to the bar, so the figure
  survives being read in greyscale or by someone who cannot distinguish
  teal from amber.
* Proxy cards print "research index; not measured biofilm quantity" on
  their face, at the same size as the rest of the caption rather than in
  a footnote.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    KeepInFrame,
    KeepTogether,
    PageBreak,
    Spacer,
    Table,
    TableStyle,
)

from openbiota.pdflinks import Paragraph, section_heading

CONTENT_WIDTH = 170 * mm
#: Two equal columns. Fixed so the two axes are visually comparable and
#: neither heading can be given more room than the other.
COL_WIDTH = CONTENT_WIDTH / 2


def _hex(c: colors.Color) -> str:
    return f"#{int(c.red * 255):02X}{int(c.green * 255):02X}{int(c.blue * 255):02X}"


def _esc(text: object) -> str:
    """Escape for ReportLab's mini-HTML."""
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _band_symbol(band: str) -> str:
    """An accessible symbol per band, so colour is never the only channel."""
    return {
        "lower reference range": "\u25bc",
        "middle reference range": "\u25c6",
        "upper reference range": "\u25b2",
    }.get(band, "\u2014")


def _bar(
    percentile: float | None,
    colour: colors.Color,
    width: float,
) -> Table:
    """A fixed-width 0-100 scale. Identical geometry on every card."""
    from openbiota.pdfreport import RULE

    track = width
    filled = (
        0.0 if percentile is None
        else max(0.0, min(100.0, percentile)) / 100.0 * track
    )
    rest = max(0.0, track - filled)
    # Two cells: filled and unfilled. Zero-height text, pure geometry.
    row = [["", ""]]
    bar = Table(row, colWidths=[filled or 0.01, rest or 0.01], rowHeights=[3.2 * mm])
    bar.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, 0), colour if percentile is not None else RULE),
                ("BACKGROUND", (1, 0), (1, 0), RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                ("ROUNDEDCORNERS", [2, 2, 2, 2]),
            ]
        )
    )
    return bar


def _card_block(
    card: Mapping[str, Any],
    st: Mapping[str, ParagraphStyle],
    colour: colors.Color,
    width: float,
) -> list[Any]:
    """One quantitative card: number, scale, scope, coverage, reference."""
    from openbiota.pdfreport import (
        ACCENT_DARK,
        INK_FAINT,
        INK_SOFT,
    )

    inner = width - 10 * mm
    out: list[Any] = []

    out.append(
        Paragraph(
            f"<font size='6.4' color='{_hex(ACCENT_DARK)}'><b>"
            f"{_esc(card['card'])} \u00b7 {_esc(card['label']).upper()}</b></font>",
            ParagraphStyle(
                "bfk", parent=st["small"], fontSize=6.4, leading=9,
                spaceAfter=1, spaceBefore=0,
            ),
        )
    )

    pct = card.get("reference_percentile")
    band = str(card.get("band") or "not computed")
    if pct is None:
        # Never a zero, never green. The reason is the content.
        out.append(
            Paragraph(
                f"<font size='13' color='{_hex(INK_SOFT)}'><b>Not scored</b></font>",
                ParagraphStyle(
                    "bfv", parent=st["body"], fontSize=13, leading=17,
                    spaceAfter=1, spaceBefore=1,
                ),
            )
        )
        out.append(_bar(None, colour, inner))
        reason = _status_sentence(card)
        out.append(
            Paragraph(
                f"<font size='7' color='{_hex(INK_SOFT)}'>{_esc(reason)}</font>",
                ParagraphStyle(
                    "bfr", parent=st["small"], fontSize=7, leading=9.8,
                    spaceAfter=0, spaceBefore=1.5,
                ),
            )
        )
        return out

    out.append(
        Paragraph(
            f"<font size='21' color='{_hex(colour)}'><b>{pct:.0f}</b></font>"
            f"<font size='8' color='{_hex(INK_FAINT)}'> / 100 percentile</font>",
            ParagraphStyle(
                "bfv", parent=st["body"], fontSize=21, leading=25,
                spaceAfter=1, spaceBefore=1,
            ),
        )
    )
    out.append(_bar(pct, colour, inner))
    out.append(
        Paragraph(
            f"<font size='7.4' color='{_hex(colour)}'><b>{_band_symbol(band)} "
            f"{_esc(band)}</b></font>",
            ParagraphStyle(
                "bfb", parent=st["small"], fontSize=7.4, leading=10,
                spaceAfter=0, spaceBefore=1.5,
            ),
        )
    )

    n = card.get("reference_n")
    cov = card.get("coverage") or {}
    scope_names = ", ".join(card.get("scope_names") or []) or "\u2014"
    meta = (
        f"Scope: limited to {_esc(scope_names)}. "
        f"Reference: {n:,} adults. "
        f"Panel coverage: {cov.get('supported', 0)} of {cov.get('assayed', 0)} "
        "features detected."
    )
    out.append(
        Paragraph(
            f"<font size='6.8' color='{_hex(INK_SOFT)}'>{meta}</font>",
            ParagraphStyle(
                "bfm", parent=st["small"], fontSize=6.8, leading=9.4,
                spaceAfter=0, spaceBefore=1.5,
            ),
        )
    )

    interval = (card.get("intervals") or {}).get("reference") or {}
    if interval.get("low") is not None:
        out.append(
            Paragraph(
                f"<font size='6.6' color='{_hex(INK_FAINT)}'>Reference interval "
                f"{interval['low']:.0f}\u2013{interval['high']:.0f} "
                f"({interval.get('valid', 0)} resamples of the reference cohort).</font>",
                ParagraphStyle(
                    "bfi", parent=st["small"], fontSize=6.6, leading=9,
                    spaceAfter=0, spaceBefore=1,
                ),
            )
        )
    return out


def _status_sentence(card: Mapping[str, Any]) -> str:
    """Plain words for why a card has no number."""
    status = str(card.get("status") or "")
    return {
        "no_validated_panel": (
            "No mechanism in this axis has a validated sequence panel yet. "
            "This is a gap in the reference data, not a finding of absence."
        ),
        "out_of_domain": (
            "The reference cohort is adults. It cannot rank a child, and an "
            "adult tail would not be a pediatric abnormality."
        ),
        "reference_unavailable": "The reference cohort could not be loaded.",
        "insufficient_reference": "Too few reference participants to rank against.",
        "reference_non_discriminating": (
            "Every reference participant has the same value here, so there is "
            "no rank to report."
        ),
        "not_computable": "The measurement could not be ranked; see the detail below.",
    }.get(status, status.replace("_", " ") or "Not computed.")


def _heading_column(
    title: str,
    cards: Sequence[Mapping[str, Any]],
    st: Mapping[str, ParagraphStyle],
    colour: colors.Color,
) -> list[Any]:
    """One of the two headings with its (at most two) cards."""
    from openbiota.pdfreport import INK, INK_SOFT, RULE

    out: list[Any] = [
        Paragraph(
            f"<font size='9.5' color='{_hex(INK)}'><b>{_esc(title)}</b></font>",
            ParagraphStyle(
                "bfh", parent=st["body"], fontSize=9.5, leading=13,
                spaceAfter=2, spaceBefore=0,
            ),
        )
    ]
    scored = [c for c in cards if c.get("reference_percentile") is not None]
    if not scored:
        out.append(Paragraph(
            f"<font size='7.4' color='{_hex(INK_SOFT)}'>Nothing on this side of the "
            "ledger could be scored for you. The axes are listed under the grid.</font>",
            ParagraphStyle("bfe", parent=st["small"], fontSize=7.4, leading=10),
        ))
        return out
    for i, card in enumerate(scored):
        if i:
            out.append(Spacer(1, 2.2 * mm))
            rule = Table([[""]], colWidths=[COL_WIDTH - 10 * mm], rowHeights=[0.4])
            rule.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, -1), RULE),
                        ("LEFTPADDING", (0, 0), (-1, -1), 0),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                        ("TOPPADDING", (0, 0), (-1, -1), 0),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                    ]
                )
            )
            out.append(rule)
            out.append(Spacer(1, 2.2 * mm))
        out.extend(_card_block(card, st, colour, COL_WIDTH))
    return out


def _score_meaning(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    cards: Sequence[Mapping[str, Any]],
) -> None:
    """What the numbers above actually say, in words.

    A percentile with no interpretation is a number the reader has to guess
    at. This says what was counted, what the position means, what it does not
    mean, and what the unscorable axes are - once, in plain language, instead
    of repeating a disclaimer inside every card.
    """
    from openbiota.pdfreport import ACCENT_DARK, INK, INK_FAINT, INK_SOFT, RULE

    scored = [c for c in cards if c.get("reference_percentile") is not None]
    if not scored:
        return

    n = max((c.get("reference_n") or 0) for c in scored)
    story.append(Paragraph(
        f"<font size='8.6' color='{_hex(INK)}'><b>What these numbers mean</b></font>",
        ParagraphStyle("bfsm", parent=st["body"], fontSize=8.6, leading=12,
                       spaceAfter=2),
    ))

    rows: list[list[Any]] = []
    for c in scored:
        pct = float(c["reference_percentile"])
        label = str(c.get("label") or c.get("card") or "")
        higher = "Concerning" in str(c.get("heading") or "")
        sense = (
            "Higher means more of your community is made up of organisms that "
            "take part in biofilm growth."
            if higher else
            "Higher means more of the community that keeps the gut lining "
            "intact and well fed."
        )
        verdict = (
            f"You are around the middle of the reference range: {pct:.0f} means "
            f"roughly {pct:.0f} in every 100 reference adults measure lower than "
            "you, and the rest measure higher."
        ) if 25 <= pct <= 75 else (
            f"You sit in the upper part of the range: about {100 - pct:.0f} in "
            "every 100 reference adults measure higher than you."
            if pct > 75 else
            f"You sit in the lower part of the range: about {pct:.0f} in every "
            "100 reference adults measure lower than you."
        )
        from openbiota.pdfreport import ACCENT, AMBER

        organisms = _organism_phrase(c, AMBER if higher else ACCENT)
        rows.append([
            Paragraph(
                f"<font size='15' color='{_hex(ACCENT_DARK)}'><b>{pct:.0f}</b></font>"
                f"<font size='7' color='{_hex(INK_FAINT)}'> / 100</font>",
                ParagraphStyle("bfsn", parent=st["body"], fontSize=15, leading=18)),
            Paragraph(
                f"<font size='7.6'><b>{_esc(label)}</b></font><br/>"
                f"<font size='7' color='{_hex(INK_SOFT)}'>{_esc(verdict)} {_esc(sense)}</font>"
                + (f"<br/><font size='7.2'>{organisms}</font>" if organisms else ""),
                ParagraphStyle("bfsb", parent=st["small"], fontSize=7, leading=9.8)),
        ])

    t = Table(rows, colWidths=[16 * mm, CONTENT_WIDTH - 16 * mm], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (0, -1), 0),
        ("RIGHTPADDING", (0, 0), (0, -1), 2 * mm),
        ("LEFTPADDING", (1, 0), (1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
    ]))
    story.append(t)
    story.append(Spacer(1, 1.6 * mm))

    story.append(Paragraph(
        f"<font size='7' color='{_hex(INK_SOFT)}'>A middle reading is the common "
        "one and is not a finding. These two numbers are never added, subtracted "
        "or averaged into a single biofilm score: a high concerning reading and a "
        "high protective reading can both be true at once, and collapsing them "
        "would hide exactly the case worth seeing. Both are read from which "
        f"organisms are present and in what proportion, against {n:,} reference "
        "adults.</font>",
        st["small"]))


def biofilm_overview(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    section: int,
    biofilm: Mapping[str, Any] | None,
    detail_section: int,  # noqa: ARG001 - the caller passes every section number; the cards link by anchor name
) -> None:
    """The front-loaded overview: two headings, four cards, key context.

    One page, always. The section is built into its own list and shrunk to
    fit, so a run with more findings than another does not spill two lines
    onto a second page.
    """
    from openbiota.pdfreport import (
        ACCENT,
        AMBER,
        BODY_FRAME_HEIGHT,
        RULE,
    )

    outer = story
    # The heading carries this section's link destination; a shrinking
    # KeepInFrame drops zero-height flowables, so it stays outside.
    outer.extend(
        section_heading(section, "Biofilm-related potential", st["h1"])
    )
    story = []

    if not biofilm:
        outer.append(
            Paragraph(
                "The biofilm module did not run for this sample.", st["body"]
            )
        )
        return

    story.append(
        Paragraph(
            "Biofilms are organised communities of microbes held together by a "
            "shared matrix. They are a normal way for microbes to live, not a "
            "disease in themselves. Some support stable colonisation and a "
            "protective relationship with the intestinal lining; others help "
            "organisms persist, grow too close to tissue, or take part in "
            "inflammation. <b>This page reports concerning and protective "
            "readings separately, so that one never cancels the other.</b>",
            st["body"],
        )
    )
    story.append(Spacer(1, 3 * mm))

    cards = list(biofilm.get("cards") or [])
    concerning = [c for c in cards if c.get("heading", "").startswith("Concerning")]
    protective = [c for c in cards if c.get("heading", "").startswith("Protective")]

    # Each card carries the organisms that made its number, inside the
    # card and sized to the card. It used to carry a table that named the
    # card again, at full-page column widths, so it clipped and said
    # nothing a reader did not already have.
    inner = COL_WIDTH - 10 * mm
    left = _heading_column("Concerning biofilm-forming (bad)", concerning, st, AMBER)
    for card in concerning:
        if card.get("reference_percentile") is not None:
            _organisms_table(left, st, card, title="What drove this reading",
                             colour=AMBER, width=inner)
    right = _heading_column("Protective gut lining (good)", protective, st, ACCENT)
    for card in protective:
        if card.get("reference_percentile") is not None:
            _organisms_table(right, st, card, title="What contributed to this score",
                             colour=ACCENT, width=inner)

    grid = Table([[left, right]], colWidths=[COL_WIDTH, COL_WIDTH], hAlign="LEFT")
    grid.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BACKGROUND", (0, 0), (0, 0), colors.HexColor("#FDF8EE")),
                ("BACKGROUND", (1, 0), (1, 0), colors.HexColor("#EEF7F7")),
                ("BOX", (0, 0), (-1, -1), 0.5, RULE),
                ("LINEAFTER", (0, 0), (0, -1), 0.5, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 4 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5 * mm),
                ("ROUNDEDCORNERS", [5, 5, 5, 5]),
            ]
        )
    )
    story.append(grid)
    story.append(Spacer(1, 2 * mm))
    _score_meaning(story, st, cards)
    outer.append(KeepInFrame(CONTENT_WIDTH, BODY_FRAME_HEIGHT - 2 * mm, story,
                             mode="shrink", hAlign="LEFT"))


def _findings_table(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    group: Sequence[Mapping[str, Any]],
    *,
    title: str = "",
    colour: Any = None,
) -> None:
    """One ranked list, in full, with the driver behind each row.

    Nothing is truncated: a reader asking what produced a reading is owed
    all of it, not the top few.
    """
    from openbiota.pdfreport import INK_FAINT, PANEL_BG, RULE

    if not group:
        return
    if title:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<font size='7.4' color='{_hex(colour) if colour is not None else _hex(INK_FAINT)}'>"
            f"<b>{title}</b></font>", st["body"]))
    rows: list[list[Any]] = []
    for f in group:
        from openbiota.biofilm.reference import feature_label

        driver = f.get("largest_driver") or {}
        symbol = driver.get("feature")
        driver_name = _esc(feature_label(symbol) if symbol else "\u2014")
        detected = driver.get("detected", True)
        # A feature that was not detected has no abundance to print. "0%"
        # would read as a measured zero, and its rank is a tie among every
        # other participant who also has none - not a measurement of this
        # person. Say what it is instead.
        if detected:
            contribution = (
                f"at {driver.get('value', 0):.3g}% of classified reads "
                f"(rank {driver.get('percentile', 0):.0f})"
            )
        else:
            contribution = (
                "not detected above this assay's limit; its rank is a tie "
                "with every other participant who also has none"
            )
        rows.append(
            [
                Paragraph(
                    f"<font size='7.6'><b>{_esc(f['label'])}</b></font>",
                    st["small"],
                ),
                Paragraph(
                    f"<font size='7.6'>{f['percentile']:.0f}</font>"
                    f"<font size='6.2' color='{_hex(INK_FAINT)}'> "
                    f"{_band_symbol(str(f.get('band', '')))}</font>",
                    st["small"],
                ),
                Paragraph(
                    f"<font size='7'>Largest contributor: <b>{driver_name}</b>, "
                    f"{contribution}.</font>",
                    st["small"],
                ),
            ]
        )
    table = Table(
        rows, colWidths=[62 * mm, 16 * mm, CONTENT_WIDTH - 78 * mm], hAlign="LEFT"
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
                ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                ("INNERGRID", (0, 0), (-1, -1), 0.4, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.8 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2 * mm),
            ]
        )
    )
    story.append(table)


def _card_organisms(card: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The organisms behind a card, from its own feature detail.

    Each feature is an organism or a genus; the file carries its level, its
    rank against the reference, whether it was detected, and which species
    were summed into it. Sorted by rank, highest first, because the
    organism at the top is the one most responsible for the reading.
    """
    from openbiota.biofilm.reference import feature_label

    detail = card.get("feature_detail") or {}
    values = card.get("feature_values") or {}
    ranks = card.get("feature_percentiles") or {}
    out: list[dict[str, Any]] = []
    for symbol in set(values) | set(detail):
        info = detail.get(symbol) or {}
        leaves = info.get("leaves_detected") or {}
        out.append({
            "symbol": symbol,
            "name": feature_label(str(symbol)),
            "value": float(values.get(symbol, info.get("value", 0.0)) or 0.0),
            "percentile": ranks.get(symbol),
            "detected": bool(info.get("detected", bool(values.get(symbol)))),
            "species": sorted(
                ((str(k).replace("_", " "), float(v)) for k, v in leaves.items()),
                key=lambda kv: -kv[1],
            ),
            "inventory_note": info.get("inventory_note"),
        })
    out.sort(key=lambda o: -(o["percentile"] or 0.0))
    return out


def _abbrev(species: str) -> str:
    """``Coprococcus eutactus`` -> ``C. eutactus`` for the species sub-line."""
    parts = species.split(" ", 1)
    return f"{parts[0][0]}. {parts[1]}" if len(parts) == 2 else species


def _organisms_table(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    card: Mapping[str, Any],
    *,
    title: str,
    colour: Any,
    width: float,
) -> None:
    """The organisms that made a card's number, inside the card.

    One row per organism: name, level, rank. An organism that was not
    found is printed at zero, because zero is what was measured; the rank
    column says where zero sits among the reference adults.
    """
    from openbiota.pdfreport import INK_FAINT, INK_SOFT, PANEL_BG, RULE, _ordinal

    organisms = _card_organisms(card)
    if not organisms:
        return
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        f"<font size='7.4' color='{_hex(colour)}'><b>{title}</b></font>", st["body"]))
    rows: list[list[Any]] = [[
        Paragraph(f"<font size='6' color='{_hex(INK_SOFT)}'><b>{h}</b></font>", st["small"])
        for h in ("ORGANISM", "LEVEL", "RANK")
    ]]
    for o in organisms:
        name = f"<font size='7.4'><b><i>{_esc(o['name'])}</i></b></font>"
        if len(o["species"]) > 1:
            parts = ", ".join(f"{_esc(_abbrev(sp))} {v:.2g}%" for sp, v in o["species"])
            name += f"<br/><font size='6.2' color='{_hex(INK_SOFT)}'>{parts}</font>"
        if o["detected"]:
            level = f"<font size='7.4'><b>{o['value']:.2f}%</b></font>"
        else:
            level = (f"<font size='7.4'><b>0%</b></font>"
                     f"<font size='6' color='{_hex(INK_FAINT)}'> not detected"
                     + (" by the reference catalogue" if o.get("inventory_note") else "")
                     + "</font>")
        if o.get("inventory_note"):
            name += f"<br/><font size='6' color='{_hex(INK_FAINT)}'>{_esc(str(o['inventory_note']))}</font>"
        rows.append([
            Paragraph(name, st["small"]),
            Paragraph(level, st["small"]),
            Paragraph(f"<font size='7.4' color='{_hex(colour)}'><b>"
                      f"{_ordinal(o['percentile'])}</b></font>", st["small"]),
        ])
    table = Table(rows, colWidths=[width * 0.56, width * 0.26, width * 0.18], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
        ("BOX", (0, 0), (-1, -1), 0.4, RULE),
        ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("TOPPADDING", (0, 0), (-1, 0), 1 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 1 * mm),
    ]))
    story.append(table)


def _organism_phrase(card: Mapping[str, Any], colour: Any) -> str:
    """The organisms and their levels as one highlighted phrase."""
    from openbiota.pdfreport import _ordinal

    organisms = _card_organisms(card)
    if not organisms:
        return ""
    bits = []
    for o in organisms:
        level = f"{o['value']:.1f}%" if o["detected"] else "0%, not detected"
        bits.append(
            f"<font color='{_hex(colour)}'><b><i>{_esc(o['name'])}</i></b></font> "
            f"{level} ({_ordinal(o['percentile'])})"
        )
    found = [o for o in organisms if o["detected"]]
    if not found:
        tail = (" \u2014 none of these organisms was found in your sample, which is what "
                "places the reading low.")
    elif len(found) < len(organisms):
        tail = " \u2014 the organisms found are what lift the reading."
    else:
        tail = "."
    return "Made up of " + "; ".join(bits) + tail



def _actual_state_strip(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    biofilm: Mapping[str, Any],
) -> None:
    """What would show an actual biofilm, and whether any of it was supplied.

    Required by spec v08.0 §4. Every row here is a measurement DNA cannot
    stand in for, so the strip is never auto-filled: it exists to make the
    difference between "measured and negative" and "never measured"
    impossible to miss.
    """
    from openbiota.pdfreport import ACCENT, INK_FAINT, INK_SOFT, RULE, SLATE

    rows_src = list(biofilm.get("actual_state") or [])
    if not rows_src:
        return
    story.append(
        Paragraph(
            "<font size='8.2'><b>What would show an actual biofilm</b></font>",
            st["body"],
        )
    )
    story.append(
        Paragraph(
            f"<font size='7' color='{_hex(INK_SOFT)}'>None of these can be "
            "derived from DNA. Each row is blank unless the measurement was "
            "actually supplied, so that \u201cnot measured\u201d never reads as "
            "\u201cmeasured and normal\u201d.</font>",
            st["small"],
        )
    )
    story.append(Spacer(1, 1.5 * mm))

    cells: list[list[Any]] = []
    for r in rows_src:
        supplied = r.get("state") == "supplied"
        colour = ACCENT if supplied else SLATE
        # An em dash, not a hollow circle: Helvetica has no U+25CB and
        # substitutes a filled box, which would read as "present" - the
        # exact opposite of what this row means.
        mark = "\u25cf supplied" if supplied else "\u2014 not supplied"
        cells.append(
            [
                Paragraph(
                    f"<font size='7.2'><b>{_esc(r['observation'])}</b></font>",
                    st["small"],
                ),
                Paragraph(
                    f"<font size='6.8' color='{_hex(colour)}'><b>{mark}</b></font>",
                    st["small"],
                ),
                Paragraph(
                    f"<font size='6.8' color='{_hex(INK_FAINT)}'>"
                    f"{_esc(r['would_come_from'])} \u2014 {_esc(r['why_it_matters'])}"
                    "</font>",
                    st["small"],
                ),
            ]
        )
    table = Table(
        cells, colWidths=[48 * mm, 24 * mm, CONTENT_WIDTH - 72 * mm], hAlign="LEFT"
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.6 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.8 * mm),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 3 * mm))


def _headline_uncertainty(biofilm: Mapping[str, Any]) -> str:
    """The single most consequential gap, chosen deterministically."""
    if biofilm.get("out_of_domain"):
        return str(biofilm["out_of_domain"])
    cards = list(biofilm.get("cards") or [])
    if any(c.get("status") == "no_validated_panel" for c in cards):
        return (
            "Neither mechanism card could be computed, because no biofilm "
            "mechanism has a validated sequence panel yet. The two readings "
            "above are community patterns, which is a weaker kind of evidence "
            "than a measured mechanism, and they do not substitute for one."
        )
    if biofilm.get("domain_caveat"):
        return str(biofilm["domain_caveat"])
    return ""


def biofilm_detail(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    section: int,
    biofilm: Mapping[str, Any] | None,
) -> None:
    """The complete inventory: every card, module, source and action."""
    from openbiota.pdfreport import (
        ACCENT,
        AMBER,
        CORAL,
        INK_SOFT,
        PANEL_BG,
        RULE,
        SLATE,
    )

    story.extend(
        section_heading(section, "Biofilm findings in full", st["h1"])
    )
    if not biofilm:
        story.append(Paragraph("The biofilm module did not run.", st["body"]))
        return

    # ---- what was measured -------------------------------------------------
    story.append(
        Paragraph(
            "This section lists every biofilm measurement, every module in the "
            "registry with its actual status, and every matched piece of "
            "evidence \u2014 including the evidence that points the other way. "
            "Nothing here is truncated to a top-N list.",
            st["body"],
        )
    )
    story.append(Spacer(1, 3 * mm))

    ref = biofilm.get("reference") or {}
    if ref:
        story.append(
            Paragraph("<b>The reference population</b>", st["h3"])
            if "h3" in st
            else Paragraph("<b>The reference population</b>", st["body"])
        )
        story.append(
            Paragraph(
                f"Percentiles are ranks against <b>{ref.get('n_participants', 0):,} "
                f"independent adults</b> from {ref.get('n_studies', 0)} public "
                f"control cohorts, aged {ref.get('age_range', ['?', '?'])[0]:.0f} to "
                f"{ref.get('age_range', ['?', '?'])[1]:.0f} (median "
                f"{ref.get('age_median', 0):.0f}), profiled on the same taxonomy "
                "as this sample. "
                f"{ref.get('n_source_samples', 0):,} samples were collapsed to "
                f"{ref.get('n_participants', 0):,} people: several of the source "
                "studies sampled the same person repeatedly, and counting those "
                "as separate people would have made the reference look larger "
                "and steadier than it is.",
                st["small"],
            )
        )
        story.append(Spacer(1, 1.5 * mm))
        story.append(
            Paragraph(
                f"<font color='{_hex(INK_SOFT)}'>{_esc(ref.get('not_established', ''))}"
                "</font>",
                st["small"],
            )
        )
        story.append(Spacer(1, 4 * mm))

    # ---- detail card per measurement --------------------------------------
    # Only measurements that produced a number get a card. An axis with no
    # published panel has nothing to detail: a full card headed "Not scored"
    # takes the space of a result and delivers a disclaimer. Those axes are
    # named once underneath instead.
    all_cards = list(biofilm.get("cards") or [])
    for card in all_cards:
        if card.get("reference_percentile") is None:
            continue
        story.append(_detail_card(card, st))
    skipped = [c for c in all_cards if c.get("reference_percentile") is None]
    if skipped:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            "<b>Axes with no measurement</b>", st["body"]))
        for card in skipped:
            story.append(Paragraph(
                f"<font size='7.4' color='{_hex(INK_SOFT)}'>\u2022 "
                f"<b>{_esc(str(card.get('label') or card.get('card') or ''))}</b> "
                f"\u2014 {_esc(_status_sentence(card))}</font>",
                st["small"]))
        story.append(Spacer(1, 3 * mm))
        story.append(Spacer(1, 3 * mm))

    # ---- agreement table ---------------------------------------------------
    rows_src = [r for r in (biofilm.get("agreement_table") or []) if r["row"] != "_concern"]
    concern_row = next(
        (r for r in (biofilm.get("agreement_table") or []) if r["row"] == "_concern"),
        None,
    )
    if rows_src:
        story.append(Paragraph("<b>Do the different kinds of evidence agree?</b>",
                               st["body"]))
        story.append(
            Paragraph(
                f"<font size='7' color='{_hex(INK_SOFT)}'>Each row is a different "
                "way of looking at the same named concern"
                + (
                    f": {_esc(concern_row['state'])}."
                    if concern_row
                    else "."
                )
                + " There is deliberately no combined confidence number, because "
                "two scores computed from the same reads are not independent "
                "confirmations of each other.</font>",
                st["small"],
            )
        )
        story.append(Spacer(1, 1.5 * mm))
        state_colour = {
            "supports_named_concern": CORAL,
            "opposes_named_concern": ACCENT,
            "mixed": AMBER,
            "not_measured": SLATE,
            "not_applicable": SLATE,
        }
        rows: list[list[Any]] = []
        for r in rows_src:
            colour = state_colour.get(r["state"], SLATE)
            rows.append(
                [
                    Paragraph(f"<font size='7.6'><b>{_esc(r['row'])}</b></font>",
                              st["small"]),
                    Paragraph(
                        f"<font size='7.4' color='{_hex(colour)}'><b>"
                        f"{_esc(r['state'].replace('_', ' '))}</b></font>",
                        st["small"],
                    ),
                    Paragraph(f"<font size='7'>{_esc(r['detail'])}</font>", st["small"]),
                ]
            )
        table = Table(rows, colWidths=[42 * mm, 30 * mm, CONTENT_WIDTH - 72 * mm],
                      hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
                    ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                    ("INNERGRID", (0, 0), (-1, -1), 0.4, RULE),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2.5 * mm),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 2.5 * mm),
                    ("TOPPADDING", (0, 0), (-1, -1), 1.8 * mm),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2 * mm),
                ]
            )
        )
        story.append(table)
        story.append(Spacer(1, 4 * mm))

    # ---- module census -----------------------------------------------------
    story.append(PageBreak())
    story.append(Paragraph("<b>Every mechanism in the registry</b>", st["body"]))
    story.append(
        Paragraph(
            f"<font size='7' color='{_hex(INK_SOFT)}'>All "
            f"{len(biofilm.get('module_results') or [])} modules appear here with "
            "their real status. A module that cannot be measured says so; it is "
            "not left out, because a silently missing mechanism is "
            "indistinguishable from one that was checked and found absent.</font>",
            st["small"],
        )
    )
    story.append(Spacer(1, 1.5 * mm))

    state_label = {
        "implemented": ("implemented", ACCENT),
        "candidate_not_validated": ("candidate, not validated", AMBER),
        "context_only": ("context only", SLATE),
        "source_unavailable": ("no reference panel yet", SLATE),
        "not_measurable_from_current_assay": ("not measurable from stool DNA", SLATE),
    }
    rows = [
        [
            Paragraph("<font size='6.6'><b>ID</b></font>", st["small"]),
            Paragraph("<font size='6.6'><b>MECHANISM</b></font>", st["small"]),
            Paragraph("<font size='6.6'><b>STATUS</b></font>", st["small"]),
            Paragraph("<font size='6.6'><b>AXIS</b></font>", st["small"]),
            Paragraph("<font size='6.6'><b>EVIDENCE</b></font>", st["small"]),
        ]
    ]
    for m in biofilm.get("module_results") or []:
        label, colour = state_label.get(
            str(m.get("census_state")), (str(m.get("census_state")), SLATE)
        )
        axis = str((m.get("score") or {}).get("candidate_axis", "")).replace("_", " ")
        rows.append(
            [
                Paragraph(f"<font size='6.8'>{_esc(m['module_id'])}</font>", st["small"]),
                Paragraph(f"<font size='6.8'>{_esc(m['label'])}</font>", st["small"]),
                Paragraph(
                    f"<font size='6.6' color='{_hex(colour)}'>{_esc(label)}</font>",
                    st["small"],
                ),
                Paragraph(f"<font size='6.6'>{_esc(axis)}</font>", st["small"]),
                Paragraph(
                    f"<font size='6.6'>"
                    f"{_esc(str(m.get('biological_evidence_status', '')).replace('_', ' '))}"
                    "</font>",
                    st["small"],
                ),
            ]
        )
    table = Table(
        rows,
        colWidths=[14 * mm, 62 * mm, 34 * mm, 26 * mm, CONTENT_WIDTH - 136 * mm],
        hAlign="LEFT",
        repeatRows=1,
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EDF1F4")),
                ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                ("INNERGRID", (0, 0), (-1, -1), 0.3, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.3 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5 * mm),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 4 * mm))

    # ---- matched evidence --------------------------------------------------
    actions = list(biofilm.get("intervention_candidates") or [])
    if actions:
        story.append(PageBreak())
        story.append(Paragraph("<b>Matched research evidence</b>", st["body"]))
        story.append(
            Paragraph(
                f"<font size='7' color='{_hex(INK_SOFT)}'>{len(actions)} records "
                "matched the findings above. These are things to read, not things "
                "to take. Where a study points the other way it is listed beside "
                "the one that points favourably, at the same size. Laboratory "
                "concentrations are shown as laboratory concentrations and are "
                "never converted into a dose.</font>",
                st["small"],
            )
        )
        story.append(Spacer(1, 2 * mm))
        for act in actions:
            story.append(_action_card(act, st))
            story.append(Spacer(1, 2.2 * mm))

    # ---- education ---------------------------------------------------------
    story.append(PageBreak())
    story.append(Paragraph("<b>Reading this section</b>", st["body"]))
    for q, a in _EDUCATION:
        story.append(Spacer(1, 2 * mm))
        story.append(
            Paragraph(f"<font size='8'><b>{_esc(q)}</b></font>", st["small"])
        )
        story.append(Paragraph(f"<font size='7.6'>{_esc(a)}</font>", st["small"]))

    limits = list(biofilm.get("limits") or [])
    if limits:
        story.append(Spacer(1, 3.5 * mm))
        story.append(Paragraph("<b>What this section cannot tell you</b>", st["body"]))
        for lim in limits:
            story.append(
                Paragraph(
                    f"<font size='7.4' color='{_hex(INK_SOFT)}'>\u2022 {_esc(lim)}</font>",
                    st["small"],
                )
            )


def _detail_card(card: Mapping[str, Any], st: Mapping[str, ParagraphStyle]) -> Any:
    """One measurement in full: value, units, reference, meaning, limits."""
    from openbiota.biofilm.reference import feature_label
    from openbiota.pdfreport import ACCENT_DARK, CARD_BG, INK_SOFT, RULE

    body: list[Any] = [
        Paragraph(
            f"<font size='8' color='{_hex(ACCENT_DARK)}'><b>{_esc(card['card'])} "
            f"\u00b7 {_esc(card['label'])}</b></font>",
            st["small"],
        )
    ]
    pct = card.get("reference_percentile")
    if pct is not None:
        body.append(
            Paragraph(
                f"<font size='7.6'><b>Percentile {pct:.1f}</b> "
                f"({_esc(card.get('band', ''))}) against "
                f"{card.get('reference_n', 0):,} adults. "
                f"Calibration {_esc(card.get('calibration_id', ''))}.</font>",
                st["small"],
            )
        )
    else:
        body.append(
            Paragraph(
                f"<font size='7.6' color='{_hex(INK_SOFT)}'><b>Not scored.</b> "
                f"{_esc(_status_sentence(card))}</font>",
                st["small"],
            )
        )

    detail = card.get("feature_detail") or {}
    if detail:
        lines = []
        for name, d in detail.items():
            value = d.get("value")
            pctl = (card.get("feature_percentiles") or {}).get(name)
            detected = d.get("detected")
            leaves = d.get("leaves_detected") or {}
            shown = ", ".join(f"{k} {v:.3g}%" for k, v in sorted(leaves.items())) or "none"
            lines.append(
                f"<b>{_esc(feature_label(name))}</b> = "
                + (f"{value:.4g}%" if value is not None else "not assayed")
                + (f", rank {pctl:.0f}" if pctl is not None else "")
                + (
                    ""
                    if detected
                    else " \u2014 <i>not detected above this assay's limit</i>"
                )
                + f". Species summed: {_esc(shown)}."
            )
        body.append(Spacer(1, 1.2 * mm))
        for line in lines:
            body.append(Paragraph(f"<font size='7'>{line}</font>", st["small"]))

    for key, label in (
        ("direction_note", "Why this direction"),
        ("what_it_is_not", "What it is not"),
    ):
        if card.get(key):
            body.append(Spacer(1, 1.2 * mm))
            body.append(
                Paragraph(
                    f"<font size='7' color='{_hex(INK_SOFT)}'><b>{label}:</b> "
                    f"{_esc(card[key])}</font>",
                    st["small"],
                )
            )

    if card.get("omissions"):
        body.append(Spacer(1, 1.2 * mm))
        body.append(
            Paragraph(
                f"<font size='6.8' color='{_hex(INK_SOFT)}'>Recorded omissions: "
                + "; ".join(_esc(o) for o in card["omissions"])
                + "</font>",
                st["small"],
            )
        )

    table = Table([[body]], colWidths=[CONTENT_WIDTH], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BACKGROUND", (0, 0), (-1, -1), CARD_BG),
                ("BOX", (0, 0), (-1, -1), 0.5, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 4 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3.2 * mm),
                ("ROUNDEDCORNERS", [4, 4, 4, 4]),
            ]
        )
    )
    return KeepTogether(table)


def _action_card(act: Mapping[str, Any], st: Mapping[str, ParagraphStyle]) -> Any:
    """One evidence card, with its opposing findings attached."""
    from openbiota.pdfreport import (
        ACCENT,
        AMBER,
        CARD_BG,
        CORAL,
        INK_FAINT,
        INK_SOFT,
        RULE,
    )

    direction = str(act.get("direction", ""))
    colour = {
        "favourable": ACCENT,
        "unfavourable": CORAL,
        "mixed": AMBER,
    }.get(direction, INK_SOFT)

    body: list[Any] = [
        Paragraph(
            f"<font size='7.8'><b>{_esc(act['compound'])}</b></font> "
            f"<font size='6.6' color='{_hex(colour)}'><b>"
            f"{_esc(act.get('review_label', ''))}</b></font>",
            st["small"],
        ),
        Paragraph(
            f"<font size='7'>{_esc(act.get('result', ''))}</font>", st["small"]
        ),
    ]
    facts = (
        f"Measured endpoint: <b>{_esc(act.get('endpoint_plain', ''))}</b>. "
        f"Model: {_esc(str(act.get('model', '')).replace('_', ' '))}. "
        + (f"Strain: {_esc(act['strain'])}. " if act.get("strain") else "")
        + f"Laboratory exposure: {_esc(act.get('lab_exposure', ''))}. "
        f"Delivery to the colon: "
        f"{_esc(str(act.get('delivery', '')).replace('_', ' '))}."
    )
    body.append(Paragraph(f"<font size='6.8' color='{_hex(INK_SOFT)}'>{facts}</font>",
                          st["small"]))
    if act.get("delivery_note"):
        body.append(
            Paragraph(
                f"<font size='6.8' color='{_hex(INK_SOFT)}'>"
                f"{_esc(act['delivery_note'])}</font>",
                st["small"],
            )
        )
    for key, label in (
        ("caution", "Caution"),
        ("collateral", "Effect on protective commensals"),
        ("safety", "Safety"),
        ("funding_conflict", "Funding"),
    ):
        if act.get(key):
            body.append(
                Paragraph(
                    f"<font size='6.8' color='{_hex(INK_SOFT)}'><b>{label}:</b> "
                    f"{_esc(act[key])}</font>",
                    st["small"],
                )
            )
    if act.get("opposes"):
        body.append(
            Paragraph(
                f"<font size='6.8' color='{_hex(CORAL)}'><b>Read together with:</b> "
                + ", ".join(_esc(o) for o in act["opposes"])
                + "</font>",
                st["small"],
            )
        )
    url = act.get("url")
    if url:
        body.append(
            Paragraph(
                f"<font size='6.4' color='{_hex(INK_FAINT)}'>{_esc(act['id'])} \u00b7 "
                f'<a href="{_esc(url)}" color="#0B6466">{_esc(act.get("doi", url))}</a>'
                "</font>",
                st["small"],
            )
        )

    table = Table([[body]], colWidths=[CONTENT_WIDTH], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BACKGROUND", (0, 0), (-1, -1), CARD_BG),
                ("BOX", (0, 0), (-1, -1), 0.4, RULE),
                ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
                ("LEFTPADDING", (0, 0), (-1, -1), 3.4 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3.4 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 2.4 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6 * mm),
            ]
        )
    )
    return KeepTogether(table)


#: The education text from spec v08.0 §13.3, kept as data so the report and
#: the JSON cannot drift apart.
_EDUCATION: tuple[tuple[str, str], ...] = (
    (
        "What are biofilms?",
        "Biofilms are organised communities of microbes held together by a "
        "shared matrix. They are a common microbial way of living, not "
        "automatically a disease. In the gut, some communities can support "
        "stable colonisation and a protective relationship with the intestinal "
        "lining. Others can help pathogens persist, grow too close to tissue, "
        "or participate in inflammation. Their effect depends on the organisms, "
        "their activity, their location and the host.",
    ),
    (
        "What did this test measure?",
        "This test ranks selected biofilm-related genetic capabilities and "
        "community patterns. We show concerning mechanisms and protective "
        "support separately, so one does not cancel the other. The "
        "community-pattern scores are experimental proxies. DNA alone does not "
        "directly measure an attached biofilm, its thickness, location or "
        "current activity.",
    ),
    (
        "Why can harmful biofilms matter?",
        "Some biofilms help organisms survive immune defences or antimicrobial "
        "exposure. That can contribute to persistent or recurrent infection in "
        "the right setting. Biofilm-related tolerance is different from carrying "
        "an antibiotic-resistance gene, although both may occur together. A "
        "genetic finding alone cannot tell us which treatment will work.",
    ),
    (
        "Can harmful biofilms be reduced?",
        "Research includes targeted treatment of specific infections, dietary "
        "and community-supporting approaches, probiotics, compounds and enzymes "
        "that affect matrix or adhesion, and emerging precision therapies. Some "
        "have encouraging laboratory or animal evidence; others have human "
        "evidence in particular conditions. Each option above states exactly "
        "what was tested. Breaking up a matrix does not necessarily kill the "
        "organisms, and removing protective communities or mucus can be "
        "counterproductive.",
    ),
    (
        "What about fibrin and other organs?",
        "Fibrin is a human clotting protein that some microbes can exploit in "
        "blood or tissue infections. It is not the same as bacterial amyloid or "
        "every biofilm matrix. Biofilm-associated problems can occur on teeth, "
        "wounds, airways, urinary devices, implants and heart valves. Stool DNA "
        "does not locate or diagnose those problems or measure circulating clots.",
    ),
    (
        "What can I do with this result?",
        "Use the supported findings to identify which mechanisms deserve "
        "attention, which evidence-linked options are relevant, and what "
        "additional information would make a decision more reliable. Changes "
        "over time can show changes in the measured genetic signals; they do not "
        "by themselves prove that a biofilm has been removed.",
    ),
)


def biofilm_glance(
    st: Mapping[str, ParagraphStyle],
    *,
    biofilm: Mapping[str, Any] | None,
    section: int,
) -> list[Any]:
    """A one-line summary for the summary page."""
    if not biofilm:
        return []
    cards = {c["card"]: c for c in (biofilm.get("cards") or [])}
    hc = cards.get("H-C") or {}
    pe = cards.get("P-E") or {}
    if hc.get("reference_percentile") is None and pe.get("reference_percentile") is None:
        return []
    parts = []
    if hc.get("reference_percentile") is not None:
        parts.append(f"concerning pattern <b>{hc['reference_percentile']:.0f}</b>")
    if pe.get("reference_percentile") is not None:
        parts.append(f"protective support <b>{pe['reference_percentile']:.0f}</b>")
    return [
        Paragraph(
            "Biofilm-related potential: "
            + ", ".join(parts)
            + f" (percentile against adult controls) \u2014 <b>section {section}</b>",
            st["small"],
        )
    ]


def self_test() -> int:
    """Smoke-check both branches and the null-percentile path."""
    from openbiota.pdfreport import _styles

    failures = 0
    st = _styles()

    empty: list[Any] = []
    biofilm_overview(empty, st, section=11, biofilm=None, detail_section=19)
    if not empty:
        failures += 1

    card_null = {
        "card": "H-M", "heading": "Concerning biofilm potential",
        "label": "Harm-associated mechanisms", "proxy_id": "",
        "status": "no_validated_panel", "reference_percentile": None,
        "band": "not computed", "feature_percentiles": {}, "feature_values": {},
        "feature_detail": {}, "reference_n": None, "calibration_id": None,
        "intervals": {"analytical": None, "reference": None}, "scope": "none",
        "scope_names": [], "module_ids": [], "source_ids": [],
        "reason_codes": [], "omissions": [], "what_it_is_not": "x",
        "direction_note": "", "transport_change": "",
        "coverage": {"supported": 0, "assayed": 0}, "largest_driver": None,
        "research_index_note": "",
    }
    card_ok = {
        **card_null,
        "card": "H-C", "label": "Biofilm-associated community pattern",
        "proxy_id": "BF-PROXY-COMMUNITY-1.0", "status": "available",
        "reference_percentile": 61.0, "band": "middle reference range",
        "feature_percentiles": {"ecoli_complex": 20.0, "m_gnavus": 95.0},
        "feature_values": {"ecoli_complex": 0.0, "m_gnavus": 6.7},
        "feature_detail": {
            "ecoli_complex": {"value": 0.0, "detected": False, "leaves_detected": {}},
            "m_gnavus": {"value": 6.7, "detected": True,
                         "leaves_detected": {"Ruminococcus_gnavus": 6.7}},
        },
        "reference_n": 1528, "calibration_id": "abc123",
        "scope_names": ["Biofilm-associated community pattern"],
        "coverage": {"supported": 1, "assayed": 2},
        "largest_driver": {"feature": "m_gnavus", "percentile": 95.0, "value": 6.7},
        "research_index_note": "research index; not measured biofilm quantity",
        "intervals": {"analytical": None,
                      "reference": {"low": 55.0, "high": 66.0, "valid": 200}},
    }
    finding = {
        "id": "BF-PROXY-COMMUNITY-1.0", "label": "x", "percentile": 61.0,
        "band": "middle reference range",
        "largest_driver": {"feature": "m_gnavus", "percentile": 95.0,
                           "value": 6.7, "detected": True},
    }
    doc = {
        "cards": [card_null, card_ok],
        "findings_needing_review": [finding],
        "supportive_observations": [],
        "ranked_findings": [finding],
        "actual_state": [{
            "observation": "Attached biofilm seen directly",
            "would_come_from": "endoscopy", "why_it_matters": "y",
            "state": "not_supplied",
        }],
        "module_census": {"n_modules": 24, "by_census_state": {}},
        "module_results": [{
            "module_id": "BF-M01", "label": "x", "census_state": "implemented",
            "score": {"candidate_axis": "context_dependent"},
            "biological_evidence_status": "human_association",
        }],
        "intervention_candidates": [{
            "id": "BF-I05", "compound": "Berberine", "direction": "favourable",
            "review_label": "matched research evidence", "result": "r",
            "endpoint_plain": "formation", "model": "gut_relevant_in_vitro",
            "strain": "s", "lab_exposure": "sub-MIC", "delivery": "plausible_unverified",
            "delivery_note": "d", "caution": "c", "opposes": ["BF-I06"],
            "doi": "10.x/y", "url": "https://doi.org/10.x/y",
        }],
        "agreement_table": [
            {"row": "Community proxy", "state": "mixed", "detail": "d"},
            {"row": "_concern", "state": "named concern", "detail": "d"},
        ],
        "reference": {
            "n_participants": 1528, "n_studies": 14, "age_range": [19.0, 107.0],
            "age_median": 51.0, "n_source_samples": 3027, "not_established": "n",
        },
        "limits": ["a limit"],
    }
    full: list[Any] = []
    biofilm_overview(full, st, section=11, biofilm=doc, detail_section=19)
    if not full:
        failures += 1
    det: list[Any] = []
    biofilm_detail(det, st, section=19, biofilm=doc)
    if not det:
        failures += 1
    if not biofilm_glance(st, biofilm=doc, section=11):
        failures += 1
    return failures
