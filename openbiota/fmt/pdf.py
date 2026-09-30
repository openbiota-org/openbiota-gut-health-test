"""PDF rendering for the match report (spec §13.1, §13.2).

Reuses the main report's chrome and visual language so the matching document
looks like the rest of the product, and reads the canonical JSON only. Every
status is a word, not just a colour, and graphics stay readable in greyscale
(V7-106).
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Paragraph, Spacer, Table, TableStyle

from openbiota import __version__
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
    _Chrome,
    _document,
    _hex,
    _styles,
)

from . import dossier
from .report import LEAD_TIERS, TIER_LABELS

_GRID = TableStyle(
    [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, INK_FAINT),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor("#E3E7DE")),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]
)


def _p(text: str, style: Any) -> Paragraph:
    return Paragraph(text, style)


def _num(value: Any, digits: int = 1) -> str:
    if value is None:
        return "—"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def _table(rows: list[list[Any]], widths: list[float], st: dict[str, Any]) -> Table:
    body = [
        [
            _p(f"<font color='{_hex(INK_FAINT)}' size='6.4'>{c}</font>", st["cell"])
            if i == 0
            else _p(str(c), st["cell"])
            for c in row
        ]
        for i, row in enumerate(rows)
    ]
    t = Table(body, colWidths=widths, repeatRows=1)
    t.setStyle(_GRID)
    return t


def _table_marked(
    rows: list[list[Any]], widths: list[float], st: dict[str, Any], marked: set[int]
) -> Table:
    """As `_table`, plus a tinted box on the rows named in `marked` (1-based)."""
    t = _table(rows, widths, st)
    extra: list[Any] = []
    for i in marked:
        extra += [
            ("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F3F7F2")),
            ("BOX", (0, i), (-1, i), 0.7, ACCENT),
        ]
    if extra:
        t.setStyle(TableStyle(extra))
    return t


def _table_marked(
    rows: list[list[Any]], widths: list[float], st: dict[str, Any], marked: set[int]
) -> Table:
    """As `_table`, plus a tinted box on the rows named in `marked` (1-based)."""
    t = _table(rows, widths, st)
    extra: list[Any] = []
    for i in marked:
        extra += [
            ("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F3F7F2")),
            ("BOX", (0, i), (-1, i), 0.7, ACCENT),
        ]
    if extra:
        t.setStyle(TableStyle(extra))
    return t


def render_pdf(result: dict[str, Any], path: Path) -> Path:
    st = _styles()
    rec = result["recipient"]
    names = {row["candidate_id"]: row["person_id"] for row in result["individual_results"]}
    meta = " · ".join(
        [
            f"Recipient {rec['person_id']}",
            f"{len(result['materials'])} candidate material(s)",
            dt.datetime.now().strftime("%-d %B %Y"),
            f"openbiota v{__version__}",
        ]
    )
    chrome = _Chrome(sample=f"{rec['person_id']} · FMT matching", generated=meta, meta_line=meta)
    doc = _document(path, chrome)
    story: list[Any] = []

    # ---- 1. summary ------------------------------------------------------- #
    story.append(_p("FMT donor and donor-set matching", st["h1"]))
    story.append(
        _p(
            f"<font color='{_hex(INK_SOFT)}'>Recipient <b>{rec['person_id']}</b> "
            f"(sample {rec['sample_id']}) against {len(result['materials'])} candidate material(s). "
            "Experimental computational comparison: research use only, not a screening result and "
            "not a clinical decision.</font>",
            st["body"],
        )
    )
    ts = result["target_summary"]
    ach = result["achievable"]
    lb = result["leaderboard"]
    ranked = sorted(result["set_results"], key=lambda b: b["match_rank"])

    # ---- 1. the answer: the recommendation, then the highest scorer ------- #
    rec_id = lb.get("recommended_id")
    recommended = next((b for b in ranked if b["candidate_id"] == rec_id), None)
    top = ranked[0] if ranked else None
    headline = recommended or top

    def _cell(label: str, value: str, sub: str = "", colour: Any = INK) -> list[Any]:
        """A label above a value, as separate flowables so nothing overlaps."""
        return [
            Paragraph(
                f"<font size=6.6 color='{_hex(ACCENT_DARK)}'><b>{label}</b></font>",
                ParagraphStyle("cardk", parent=st["small"], fontSize=6.6, leading=9,
                               spaceAfter=0, spaceBefore=0),
            ),
            Paragraph(
                f"<font size=17 color='{_hex(colour)}'><b>{value}</b></font>"
                + (f"<font size=8 color='{_hex(INK_SOFT)}'> {sub}</font>" if sub else ""),
                ParagraphStyle("cardv", parent=st["body"], fontSize=17, leading=21,
                               spaceAfter=0, spaceBefore=1),
            ),
        ]

    if headline is not None:
        story.append(Spacer(1, 2.5 * mm))
        # With no clean candidate the headline is the top of the ranking,
        # which is the one introducing the fewest disease-enriched organisms
        # the recipient lacks - not the top scorer. Calling it "highest
        # score" there would point at a different row than rank 1.
        label = "RECOMMENDED" if recommended is not None else "LEAST INTRODUCED"
        block = [[
            _cell(label, " + ".join(headline["member_person_ids"])),
            _cell("MATCH SCORE", _num(headline["match_score"]), "of 100"),
            _cell("GAPS SUPPLIED", str(headline["goals"]["covered"]),
                  f"of {lb['reachable_goals']}"),
            # The donor's own gut health, on the front page. Its absence
            # from the ranking is what let a donor with a negative index
            # reach the top three on gap coverage alone.
            _cell("DONOR GUT HEALTH",
                  (f"{headline['donor_health']['gmwi2']:+.2f}"
                   if headline["donor_health"]["gmwi2"] is not None else "—"),
                  "GMWI2, 0 is neutral"
                  if headline["donor_health"]["passes_health_floor"]
                  else "below zero — dysbiotic",
                  GREEN if headline["donor_health"]["passes_health_floor"] else CORAL),
            # One source of truth with the detail table below it. These
            # disagreed: the tile counted exposures and the table counted
            # the suitability record, so the front page said nine and the
            # detail said none.
            _cell("ORGANISMS IT WOULD ADD",
                  str(len(
                      (headline.get("donor_health") or {}).get(
                          "organisms_it_would_add") or [])),
                  (f"{(headline.get('donor_health') or {}).get('abnormal_organism_count', 0)}"
                   " above the healthy range"),
                  CORAL if (headline.get("donor_health") or {}).get(
                      "abnormal_organism_count") else GREEN),
        ]]
        table = Table(block, colWidths=[CONTENT_WIDTH * f for f in (0.24, 0.16, 0.16, 0.20, 0.24)])
        table.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F3F7F2")),
                    ("BOX", (0, 0), (-1, -1), 0.9, ACCENT),
                    ("LINEAFTER", (0, 0), (-2, -1), 0.4, colors.HexColor("#DFE7DB")),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4.5 * mm),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
                    ("TOPPADDING", (0, 0), (-1, -1), 3.2 * mm),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3.6 * mm),
                ]
            )
        )
        story.append(table)
        story.append(Spacer(1, 2.2 * mm))
        if lb.get("recommended_note"):
            story.append(
                _p(f"<b>Why this one and not the top of the table.</b> {lb['recommended_note']}",
                   st["small"])
            )
        if lb.get("no_clean_candidate_note"):
            story.append(
                _p(f"<font color='{_hex(AMBER)}'><b>{lb['no_clean_candidate_note']}</b></font>",
                   st["small"])
            )
        story.append(_p(headline["rank_reason"], st["small"]))

        # The arithmetic behind the rank, printed. Every term is a number
        # shown elsewhere on this page, so a reader can check the sum.
        story.append(
            _p(
                f"<b>How this was scored.</b> {headline['suitability_explanation']}. "
                "Restoration is the share of your gaps this material supplies; donor "
                "gut health is GMWI2 on a 0\u2013100 display scale; the third term is "
                "freedom from disease patterns you do not already have. The weights "
                "are fixed and are the same for every candidate.",
                st["small"],
            )
        )

        # The selection rule itself, printed. It existed only in the JSON, so
        # a reader could not see what "RECOMMENDED" was actually chosen on -
        # which is how a rule describing a retired filter went unnoticed.
        story.append(
            _p(
                f"<font color='{_hex(INK_SOFT)}'><b>Recommended on:</b> "
                f"{lb['recommended_rule']}</font>",
                st["small"],
            )
        )

        # Why materials were set aside, named. A reader should not have to
        # infer an exclusion from a low position in a table.
        excluded = [b for b in ranked if not b["eligible"]]
        if excluded:
            reasons: dict[str, list[str]] = {}
            for body in excluded:
                for failure in body["gate_failures"]:
                    reasons.setdefault(failure, []).append(
                        " + ".join(body["member_person_ids"])
                    )
            story.append(Spacer(1, 1.5 * mm))
            story.append(
                _p(
                    f"<font color='{_hex(CORAL)}'><b>{len(excluded)} of "
                    f"{len(ranked)} candidates are not eligible to be "
                    "recommended.</b></font> Each is still ranked and shown, "
                    "with its reason:",
                    st["small"],
                )
            )
            for failure, who in sorted(reasons.items(), key=lambda kv: -len(kv[1])):
                short = sorted(set(who), key=len)[:4]
                story.append(
                    _p(
                        f"<font size=7 color='{_hex(INK_SOFT)}'>\u2022 {failure} "
                        f"\u2014 affects {len(who)} candidate(s): "
                        f"{', '.join(short)}{' …' if len(who) > 4 else ''}</font>",
                        st["small"],
                    )
                )

        # What the headline candidate carries, said next to the headline.
        # Without this the front page showed a green "0 introduced" tile and
        # the patterns appeared two pages later, so the recommended material
        # could sit at the 90th percentile for the colorectal-cancer pattern
        # against a recipient at the 15th with nothing on this page saying so.
        # Lead with the organisms, because they are the criterion now. A
        # list of twelve percentiles is unreadable and is not what decided
        # anything; the species this material would add that the recipient
        # does not already carry is both shorter and actionable.
        carried = list(headline.get("introduces_new_patterns") or [])
        members = set(headline["member_person_ids"])
        verdicts = list((headline.get("donor_health") or {}).get(
            "organisms_it_would_add") or [])
        if verdicts:
            abnormal = [v for v in verdicts if v["abnormal"] or v["is_known_pathogen"]]
            story.append(
                _p(
                    f"<b>Every organism this material would add that you do not "
                    f"already carry ({len(verdicts)}), each judged against healthy "
                    "adults.</b> "
                    + (
                        f"<font color='{_hex(CORAL)}'>{len(abnormal)} above the "
                        "normal range.</font>"
                        if abnormal
                        else f"<font color='{_hex(ACCENT_DARK)}'>None is above the "
                             "normal range for a healthy adult.</font>"
                    ),
                    st["small"],
                )
            )
            rows_o: list[list[Any]] = [[
                "ORGANISM", "IN THIS DONOR", "HEALTHY ADULTS", "VERDICT",
                "RAISED IN",
            ]]
            for v in verdicts:
                bad = v["abnormal"] or v["is_known_pathogen"]
                colour = CORAL if bad else INK_SOFT
                ref = (
                    "not in reference"
                    if v["reference_p90_percent"] is None
                    else (
                        f"median {v['reference_median_percent']:.3f}%, "
                        f"90th {v['reference_p90_percent']:.3f}%<br/>"
                        f"<font size=6>carried by "
                        f"{v['carried_by_percent_of_healthy_adults']:.0f}% of adults</font>"
                    )
                )
                rows_o.append([
                    f"<i>{v['species'].replace('_', ' ')}</i>",
                    f"<b>{v['donor_percent']:.3f}%</b>",
                    f"<font size=6.6>{ref}</font>",
                    f"<font color='{_hex(colour)}'><b>{v['verdict']}</b></font>",
                    f"<font size=6.2>{', '.join(v['raised_in_conditions'])}</font>",
                ])
            story.append(
                _table(
                    rows_o,
                    [38 * mm, 22 * mm, 42 * mm, 30 * mm, CONTENT_WIDTH - 132 * mm],
                    st,
                )
            )
            story.append(
                _p(
                    f"<font size=7 color='{_hex(INK_SOFT)}'>An organism counts against "
                    "a material only when this donor carries <b>more of it than "
                    "healthy adults normally do</b> \u2014 above their 90th "
                    "percentile \u2014 or when the pathogen screen called it. Presence "
                    "alone is not a finding: most of these are ordinary gut "
                    "commensals, and a transplant moving a normal amount of one "
                    "introduces nothing abnormal. Organisms a disease <i>depletes</i> "
                    "are never counted, because restoring those is the point of the "
                    "transfer.</font>",
                    st["small"],
                )
            )
        elif carried:
            listed = ", ".join(
                f"<b>{e['short_label']}</b> ({e['donor_percentile']:.0f}th vs your "
                f"{e['recipient_percentile']:.0f}th)"
                for e in carried
                if members & set(e["donor_person_ids"])
            )
            story.append(
                _p(
                    f"<font color='{_hex(INK_SOFT)}'>This material scores above you on "
                    f"{len(carried)} pattern{'' if len(carried) == 1 else 's'} "
                    f"\u2014 {listed} \u2014 but would add no disease-raised organism "
                    "you do not already carry, so nothing here counts against it.</font>",
                    st["small"],
                )
            )
        if ach["unreachable_goals"]:
            story.append(
                _p(
                    "<b>No candidate supplies:</b> "
                    + "; ".join(u["label"] for u in ach["unreachable_goals"])
                    + ". Those gaps sit outside every score here, which is why the scale runs "
                    "against what is reachable.",
                    st["small"],
                )
            )
        # Where the newer taxonomic database disagrees with that claim, say
        # so. Otherwise the page asserts a biological absence that the
        # better of the two reference databases contradicts.
        for row in ach.get("unreachable_contradicted_by_newer_lane") or []:
            who = ", ".join(
                f"{person} ({percent:.3f}%)"
                for person, percent in sorted(row["found_by_newer_lane_in"].items())
            )
            story.append(
                _p(
                    f"<font color='{_hex(AMBER)}'><b>The two taxonomic databases "
                    "disagree about " + str(row["label"]) + ".</b></font> "
                    + str(row["note"]) + " Found by the newer lane in: " + who + ".",
                    st["small"],
                )
            )
        story.append(
            _p(
                f"<font color='{_hex(INK_SOFT)}'>A match score is restoration evidence only: the "
                "share of your reachable gaps a candidate can supply. It is not a probability "
                "that anything will establish, not a probability of benefit, and not a safety "
                "ranking.</font>",
                st["small"],
            )
        )

    # ---- 2. the full ranked field ----------------------------------------- #
    story.append(_p("Every candidate, ranked", st["h2"]))
    story.append(
        _p(
            "Two groups. <b>Allowed</b> candidates pass every safety gate and could be "
            "used. <b>Ruled out</b> candidates fail at least one and cannot be "
            "recommended whatever they score \u2014 a donor whose own community reads as "
            "dysbiotic, or who would bring you an organism that a transfer-evidenced "
            "disease is raised in. They are still listed, with the reason, because an "
            "exclusion you cannot see is one you cannot check. Within each group the "
            "order is restoration score: the share of your reachable gaps a candidate "
            "supplies. The <b>#</b> column is the single ranking across both.",
            st["small"],
        )
    )
    header = [
        "#", "ALLOWED", "CANDIDATE", "MATCH SCORE", "GAPS SUPPLIED", "PEOPLE",
        "ORGANISMS IT WOULD ADD", "FINDINGS TO TEST", "WHY IT RANKS HERE",
    ]
    widths = [7 * mm, 15 * mm, 23 * mm, 16 * mm, 19 * mm, 12 * mm, 26 * mm,
              20 * mm, CONTENT_WIDTH - 138 * mm]

    def _row(b: dict[str, Any]) -> list[Any]:
        star = bool(b.get("recommended"))
        mark, end_ = ("<b>", "</b>") if star else ("", "")
        # Every pattern the candidate carries above the threshold, not just
        # the ones flagged `counts_against_selection`. Nothing counts against
        # selection since the percentile bridge was retired, so that filter
        # printed a green "none" for a candidate carrying four of them.
        # The organisms, not the condition names. Twelve condition names in
        # a narrow column wraps to a dozen lines per row and says nothing a
        # reader can act on; the species are fewer and are what the veto
        # actually turns on.
        carried = list(b.get("introduces_new_patterns") or [])
        verdicts_row = list((b.get("donor_health") or {}).get(
            "organisms_it_would_add") or [])
        organisms = [v["species"] for v in verdicts_row]
        if organisms:
            # Every one of them, with the abnormal ones marked. A truncated
            # "+3" was unreadable and un-checkable.
            parts = []
            for v in verdicts_row:
                name = v["species"].replace("_", " ")
                bad = v["abnormal"] or v["is_known_pathogen"]
                parts.append(
                    f"<font color='{_hex(CORAL)}'><b>{name}</b></font>" if bad
                    else f"<font color='{_hex(INK_SOFT)}'>{name}</font>"
                )
            newp_cell = "<font size=6>" + ", ".join(parts) + "</font>"
        elif carried:
            newp_cell = (
                f"<font color='{_hex(INK_SOFT)}'>none \u2014 "
                f"{len(carried)} pattern(s) above you, no new organism</font>"
            )
        else:
            newp_cell = f"<font color='{_hex(GREEN)}'>none</font>"
        # The candidate name links to its own page. A row can hold a number;
        # the sentence that makes the number mean something lives there.
        dest = dossier.candidate_dest(b["candidate_id"])
        name = " + ".join(b["member_person_ids"])
        # One column a reader can scan without decoding anything: a green
        # tick for a candidate that may be recommended, a red prohibition
        # sign for one the safety gates rule out. Everything else on the row
        # explains the verdict; this states it.
        ok = bool(b.get("eligible"))
        verdict = (
            f"<font color='{_hex(GREEN)}' size=11><b>\u2714</b></font>"
            f"<br/><font size=5.4 color='{_hex(GREEN)}'>ALLOWED</font>"
            if ok else
            f"<font color='{_hex(CORAL)}' size=11><b>\u2205</b></font>"
            f"<br/><font size=5.4 color='{_hex(CORAL)}'>RULED OUT</font>"
        )
        name_colour = "#0B6466" if ok else _hex(CORAL)
        return [
            f"{mark}{b['match_rank']}{end_}",
            verdict,
            f'<a href="#{dest}" color="{name_colour}">{mark}{name}{end_}</a>'
            f'<br/><font size=5.6 color="#8C99A6">full page \u2192</font>',
            f"{mark}{_num(b['match_score'])}{end_}",
            f"{b['goals']['covered']} of {lb['reachable_goals']}"
            + (f" ({b['goals']['partly_supplied']} partly)"
               if b["goals"]["partly_supplied"] else ""),
            str(b["unique_donor_count"]),
            newp_cell,
            f"<font color='{_hex(AMBER)}'>{b['serious_open_questions']}</font>"
            + f"<font color='{_hex(INK_SOFT)}'> + {b['other_findings']} routine</font>",
            f"<font size=6.2 color='{_hex(INK_SOFT)}'>"
            + (b.get("rank_reason_short") or "")
            + "</font>",
        ]

    def _block(bodies: list[dict[str, Any]], *, allowed: bool) -> None:
        if not bodies:
            return
        story.append(Spacer(1, 3 * mm))
        title = (
            "\u2714  Allowed as a match"
            if allowed else
            "\u2205  Ruled out \u2014 not allowed as a match"
        )
        sub = (
            "Passes every safety gate. Any of these could be used; they are "
            "ordered by how much of your ledger they restore."
            if allowed else
            "Blocked by at least one safety gate. Shown in full, with the reason, "
            "because an exclusion you cannot see is one you cannot check."
        )
        story.append(
            _p(
                f"<font size=8.6 color='{_hex(ACCENT_DARK if allowed else CORAL)}'>"
                f"<b>{title}</b></font>  <font size=7.4 color='{_hex(INK_SOFT)}'>"
                f"{len(bodies)} of {len(ranked)} candidates</font>",
                st["body"],
            )
        )
        story.append(_p(f"<font size=7 color='{_hex(INK_SOFT)}'>{sub}</font>", st["small"]))
        marked = {i for i, b in enumerate(bodies, start=1) if b.get("recommended")}
        story.append(
            _table_marked([header, *[_row(b) for b in bodies]], widths, st, marked)
        )

    # Split on eligibility, because that is the question being asked. The
    # previous split was on whether a candidate carried any pattern at all,
    # which put an allowed candidate and a vetoed one in the same group and
    # left every row looking equally available.
    _block([b for b in ranked if b.get("eligible")], allowed=True)
    _block([b for b in ranked if not b.get("eligible")], allowed=False)

    # ---- 3. what every column means --------------------------------------- #
    # ---- 2b. the evidence behind every flagged pattern -------------------- #
    exposures = lb.get("new_pattern_exposures") or []
    if exposures:
        story.append(_p("Patterns a candidate would bring that you do not have", st["h2"]))
        story.append(
            _p(
                "Each row is a pattern where a candidate sits at or above the 75th percentile and "
                "you sit at least 20 points lower. The evidence column is the strongest published "
                "transfer experiment for that condition, with its limit stated. None of this is a "
                "diagnosis of the donor and none of it is a probability for you.",
                st["small"],
            )
        )
        rows_x: list[list[Any]] = [
            ["PATTERN", "CANDIDATES", "THEM vs YOU", "TRANSFER EVIDENCE",
             "EVIDENCE TIER"]
        ]
        for e in exposures:
            colour = (
                AMBER if e["transfer_evidence_tier"] == "human_donor_transfer"
                else INK_SOFT
            )
            rows_x.append(
                [
                    f"<font color='{_hex(colour)}'><b>{e['short_label']}</b></font>",
                    ", ".join(e["donor_person_ids"]),
                    f"<b>{e['donor_percentile']:.0f}th</b> vs your "
                    f"{e['recipient_percentile']:.0f}th",
                    f"<font size=6.6>{e['evidence']} <font color='{_hex(INK_SOFT)}'>"
                    f"Limit: {e['evidence_boundary']}"
                    + (f" Sources: {'; '.join(e['sources'])}." if e["sources"] else "")
                    + "</font></font>",
                    # The tier, plainly. This cell used to read
                    # "no — human donor transfer", which parses as though the
                    # strongest evidence tier were the reason the pattern did
                    # not count. Nothing counts against selection any more, so
                    # the column now reports the evidence tier and says once,
                    # below the table, that none of it excludes a candidate.
                    (
                        f"<font color='{_hex(AMBER)}'><b>"
                        + e["transfer_evidence_tier"].replace("_", " ")
                        + "</b></font>"
                        if e["transfer_evidence_tier"] == "human_donor_transfer"
                        else f"<font color='{_hex(INK_SOFT)}'>"
                        + e["transfer_evidence_tier"].replace("_", " ")
                        + "</font>"
                    ),
                ]
            )
        story.append(
            _table(
                rows_x,
                [30 * mm, 18 * mm, 24 * mm, CONTENT_WIDTH - 100 * mm, 28 * mm],
                st,
            )
        )
        story.append(
            _p(
                "<b>What this does and does not do to the ranking.</b> A resemblance percentile "
                "is not a diagnosis and is not a transmissible thing, so <b>none of these "
                "patterns excludes a candidate and none of them changes the score</b>. The "
                "ranking above is restoration alone. For a handful of conditions, stool from "
                "patients has been put into animals and reproduced disease features against "
                "healthy-donor controls, and colorectal cancer is the clearest of them: those "
                "are marked as human-donor-transfer evidence. Choosing a donor who is high on "
                "such a pattern when you are not is a decision worth making deliberately. The "
                "tool names it and leaves the decision with you; it does not make it by "
                "reordering the table.",
                st["small"],
            )
        )

    # ---- one dossier page per candidate, linked from the table ------------ #
    # A column can hold a number; it cannot hold the sentence that makes the
    # number mean something. "15 + 12 routine" and "2 patterns above you, no
    # new organism" are unreadable on their own, so each candidate gets a
    # page and its row links to it.
    story.append(PageBreak())
    story.append(_p("Each candidate in full", st["h2"]))
    story.append(
        _p(
            "One page per candidate, in rank order. Each answers the same four "
            "questions: what it puts back that you are missing, what else it "
            "brings that is probably good, what it brings that is neither here "
            "nor there, and what it could be introducing. Then the arithmetic "
            "of its rank and every number on its row, decoded.",
            st["small"],
        )
    )
    dossier.dossier_pages(
        story, st,
        ranked=ranked,
        targets=result.get("targets") or [],
        donor_species=result.get("species_by_person") or {},
        recipient_species=result.get("recipient_species") or {},
        leaderboard=lb,
    )

    # One glossary, not two. "What each column means" and "How to read this
    # report" were separate sections defining several of the same terms in
    # different words, each starting on its own page.
    # ---- strain-level donor comparison (spec §10.1) ----------------------- #
    sm = result.get("strain_matching") or {}
    base = sm.get("baseline_comparison") or {}
    if base.get("status") == "completed":
        story.append(_p("Strain-level comparison with each donor", st["h2"]))
        story.append(
            _p(
                "Species tell you an organism is shared. Strains tell you whether it is the "
                "<b>same population</b>. Where a donor carries a population distinct from "
                "yours, that organism is a candidate for changing something; where it is "
                "indistinguishable, it would add nothing you do not already have. This is "
                "measured per donor, at single-base resolution over the sequence callable in "
                "both specimens.",
                st["small"],
            )
        )
        rows_s: list[list[Any]] = [[
            "DONOR", "ORGANISMS COMPARED", "DISTINCT POPULATIONS",
            "CLOSELY RELATED", "INDISTINGUISHABLE", "STRAIN-RESOLVED IN DONOR ONLY",
        ]]
        for donor, body in sorted((base.get("per_donor") or {}).items()):
            if body.get("status") != "completed":
                rows_s.append([
                    donor,
                    f"<font color='{_hex(INK_SOFT)}'>{body.get('why', 'not compared')}</font>",
                    "\u2014", "\u2014", "\u2014", "\u2014",
                ])
                continue
            rows_s.append([
                f"<b>{donor}</b>",
                str(body.get("organisms_compared", 0)),
                f"<b>{body.get('distinct_populations', 0)}</b>",
                str(body.get("closely_related", 0)),
                str(body.get("indistinguishable", 0)),
                str(body.get("n_organisms_resolved_in_donor_only", 0)),
            ])
        story.append(
            _table(
                rows_s,
                [18 * mm, 28 * mm, 30 * mm, 26 * mm, 30 * mm, CONTENT_WIDTH - 132 * mm],
                st,
            )
        )
        for limit in base.get("limits") or []:
            story.append(_p(f"<font color='{_hex(INK_SOFT)}'>{limit}</font>", st["small"]))
        story.append(
            _p(f"<font color='{_hex(INK_SOFT)}'>{base.get('note', '')}</font>", st["small"])
        )
        # The two capabilities that are not measured, stated as such.
        for key, label in (
            ("observed_engraftment", "Did it establish?"),
            ("prospective_establishment", "Will it establish?"),
        ):
            cap = sm.get(key) or {}
            story.append(
                _p(
                    f"<b>{label}</b> <font color='{_hex(AMBER)}'>"
                    f"{str(cap.get('status', '')).replace('_', ' ')}</font> "
                    f"<font color='{_hex(INK_SOFT)}'>{cap.get('why', '')}</font>",
                    st["small"],
                )
            )

    story.append(_p("What every number and label means", st["h2"]))
    glossary: list[tuple[str, str]] = [
        ("Match score",
         "The share of the recipient's reachable gaps this candidate supplies, each gap counted "
         f"once, with partial credit where supply is graded. {lb['reachable_goals']} of "
         f"{lb['scored_goals']} scored gaps are reachable — the rest are supplied by no candidate "
         "here — so 100 means \u201csupplies everything these materials could supply between "
         "them\u201d."),
        ("Gaps supplied",
         "The plain count behind the score. A gap is an organism or a gene capacity the recipient "
         "is missing or low in, derived from the recipient's own data before any donor was "
         "examined. \u201cPartly\u201d means present but still below the reference range."),
        ("People",
         "How many different donors the candidate involves. At equal scores fewer is ranked higher: "
         "one person means one screening history, one consent and one set of unknowns."),
        ("Serious findings",
         "Organisms that would matter if the evidence held up, and which a laboratory test would "
         "settle. None here is confirmed. This column never moves a candidate up or down the score "
         "\u2014 it breaks ties only, and it is reported in full further down."),
        ("Other findings",
         "The routine residue every stool sample produces: species a close relative makes "
         "unnameable, organisms carried without the genes that cause disease, resistance genes with "
         "no resolved carrier organism."),
        ("Organisms it would add",
         "A pattern this candidate scores high on that you do not — at or above the 75th "
         "percentile, and at least 20 percentile points above your own reading — for a condition "
         "where stool from patients has reproduced disease features in recipient animals. A "
         "pattern you already carry at a similar level is not counted: it is not new to you. "
         "These are resemblance scores, not diagnoses, and every transfer experiment behind them "
         "used an animal, usually a predisposed one. They are listed because \u201cdo not give me "
         "a condition I do not already have\u201d is a reasonable thing to ask of a donor choice."),
        ("Findings to test",
         "Organisms the sequencing found that it cannot settle. The first number is the serious "
         "ones \u2014 an organism that would matter if the evidence held up. The second is the "
         "routine residue every stool sample produces. Each finding names the test that would "
         "settle it and what the test needs: a nucleic-acid test works on a thawed aliquot of the "
         "banked sample, while a culture needs viable organisms and should be run on a fresh "
         "sample from the donor."),
        ("Identity not confirmed",
         "The sequence matched an organism, but not well enough to name the species: either the "
         "support sits in a region the organism shares with a close relative, or the organism is "
         "the same genomic species as a relative most healthy guts carry. No species name is "
         "given because it would be a guess. Nothing has been ruled in or out."),
        ("Ranking rule", lb["ranking_rule"] + " Disease-pattern resemblance does not enter "
         "this rule at all: it neither changes a score nor excludes a candidate. Patterns a "
         "candidate carries that you do not are listed separately so they can be weighed."),
    ]
    seen = {t.lower().rstrip(".") for t, _ in glossary}
    aliases = {"gaps covered": "gaps supplied", "serious organism, identity not confirmed":
               "identity not confirmed"}
    for item in result["how_to_read"]:
        term = str(item["term"])
        key = aliases.get(term.lower().rstrip("."), term.lower().rstrip("."))
        if key in seen:
            continue
        seen.add(key)
        glossary.append((term, str(item["plain"])))
    for term, text in glossary:
        story.append(_p(f"<b>{term}.</b> {text}", st["small"]))
        story.append(Spacer(1, 1.2 * mm))

    # ---- 2. individuals --------------------------------------------------- #
    story.append(_p("Every individual candidate", st["h2"]))
    story.append(
        _p(
            "Read the first column. The two index columns are the same information weighted two ways, "
            "shown together so a candidate covering the same number of gaps as another never appears to "
            "score differently for no visible reason. 'Needs a lab test' counts findings a person should "
            "look at; 'other open' counts rows with no organism-specific evidence. Neither is ever "
            "offset by coverage.",
            st["small"],
        )
    )
    concerns_by_cand: dict[str, list[dict[str, Any]]] = {}
    for c in result["concern_records"]:
        concerns_by_cand.setdefault(c["candidate_id"], []).append(c)
    rows: list[list[Any]] = [
        ["DONOR", "GAPS COVERED", f"INDEX (MAX {_num(ach['index_ceiling'], 0)})", "BY GROUP SIZE",
         "NEEDS A LAB TEST", "OTHER OPEN", "GMWI2", "SCREENING"]
    ]
    for row in sorted(result["individual_results"], key=lambda x: -x["goals"]["covered"]):
        cov, alt, g = row["coverage"], row["coverage_alternatives"], row["goals"]
        mine = concerns_by_cand.get(row["candidate_id"], [])
        serious = sum(1 for c in mine if c["tier"] in LEAD_TIERS and c["tier"] != "screening_gap")
        gm = row["context_not_in_the_score"]["gmwi2"]
        rows.append(
            [
                f"<b>{row['person_id']}</b>",
                f"<b>{g['covered']} of {g['total']}</b>",
                _num(cov["coverage"]),
                _num(alt["goal_groups_weighted_by_size"]),
                f"<font color='{_hex(AMBER)}'>{serious}</font>",
                str(len(mine) - serious),
                _num(gm["score"], 2),
                row["screening_status"].replace("_", " "),
            ]
        )
    story.append(
        _table(rows, [20 * mm, 26 * mm, 24 * mm, 24 * mm, 26 * mm, 20 * mm, 16 * mm, 24 * mm], st)
    )
    for row in sorted(result["individual_results"], key=lambda x: -x["goals"]["covered"]):
        summary = row["plain_summary"]
        if summary.startswith(row["person_id"]):
            summary = summary[len(row["person_id"]) :].lstrip()
        story.append(_p(f"<b>{row['person_id']}</b> {summary}", st["small"]))
    for row in result["excluded_candidates"]:
        story.append(
            _p(
                f"<font color='{_hex(AMBER)}'><b>Excluded on confirmed evidence:</b></font> "
                f"{row['person_id']} — analytics retained in the JSON.",
                st["small"],
            )
        )

    # ---- 3. sets ---------------------------------------------------------- #
    story.append(_p("Candidate sets against the best single material", st["h2"]))
    story.append(
        _p(
            "A set's coverage is an availability ceiling: at least one member has evidence of supplying "
            "the goal. It is not a pooled concentration, a predicted community, a persistence forecast "
            "or a preparation instruction. Zero marginal gain prints as zero.",
            st["small"],
        )
    )
    srows: list[list[Any]] = [
        ["MEMBERS", "GAPS COVERED", "INDEX", "GAIN VS BEST SINGLE", "NEEDS A LAB TEST",
         "FRONT", "STABLE"]
    ]
    for s in result["set_results"]:
        mine = [c for c in result["concern_records"] if c["candidate_id"] in s["member_material_ids"]]
        serious = sum(1 for c in mine if c["tier"] in LEAD_TIERS and c["tier"] != "screening_gap")
        srows.append(
            [
                " + ".join(s["member_person_ids"]),
                f"<b>{s['goals']['covered']} of {s['goals']['total']}</b>",
                _num(s["coverage"]["coverage"]),
                _num(s["marginals"]["gain_over_global_best_singleton"]),
                f"<font color='{_hex(AMBER)}'>{serious}</font>",
                str(s["pareto_front"]),
                f"{s['stability']['leading_front_scenarios']}/{s['stability']['scenario_count']}",
            ]
        )
    story.append(
        _table(srows, [34 * mm, 26 * mm, 18 * mm, 30 * mm, 26 * mm, 14 * mm, 16 * mm], st)
    )
    story.append(
        _p(
            "Front 1 means no other candidate is at least as good on every objective at once: coverage, "
            "how much of it is assessed, how wide the uncertainty is, and how few materials are "
            "involved. A single material on front 1 is there because involving one person is itself an "
            "advantage.",
            st["small"],
        )
    )
    if ach["unreachable_goals"]:
        story.append(
            _p(
                f"<font color='{_hex(AMBER)}'><b>Gaps no candidate can fill:</b></font> "
                + "; ".join(
                    f"{u['label']} (worth {_num(u['weight_points'], 1)} index points)"
                    for u in ach["unreachable_goals"]
                )
                + ".",
                st["small"],
            )
        )

    groups = sorted(ts["counting_groups"])
    story.append(_p("Where each candidate's coverage comes from", st["body"]))
    grows: list[list[Any]] = [["CANDIDATE", *[g.upper()[:18] for g in groups], "TOTAL"]]
    for s_ in [x for x in result["set_results"] if x["material_count"] <= 2][:10]:
        grows.append(
            [
                "+".join(s_["member_person_ids"]),
                *[_num(s_["coverage"]["per_group"].get(g), 1) for g in groups],
                f"<b>{_num(s_['coverage']['coverage'])}</b>",
            ]
        )
    gw = [26 * mm] + [(CONTENT_WIDTH - 46 * mm) / max(1, len(groups))] * len(groups) + [20 * mm]
    story.append(_table(grows, gw, st))
    story.append(
        _p(
            "One goal group, one vote: a guild and the gene panel measuring the same capacity share a "
            "group and their credit is capped at that group's weight. Equal group weights are an "
            "engineering default; the sensitivity table re-runs the ranking weighted by group size.",
            st["small"],
        )
    )

    # ---- 4. target matrix ------------------------------------------------- #
    # The matrix repeats its header row, so it does not need a page to itself;
    # forcing one left a third of the preceding page blank.
    story.append(_p("Every gap, and who has it", st["h2"]))
    story.append(
        _p(
            "<b>1.00</b> the material carries it · <b>0.00</b> the same profiler looked and did not find "
            "it, at its detection limit · <b>n/a</b> not assessed in that material · <b>ctx</b> shown "
            "for orientation and outside the score. Values in between are a graded supply index against "
            "the reference lower quartile. No negative credit exists in this tool.",
            st["small"],
        )
    )
    donor_ids = [row["candidate_id"] for row in result["individual_results"]]
    avail: dict[tuple[str, str], dict[str, Any]] = {}
    for row in result["individual_results"]:
        for cell in row["per_target_availability"]:
            avail[(cell["target_id"], row["candidate_id"])] = cell
    header = ["GOAL", "ORIGIN", "DIR", "RECIPIENT"] + [names[d][:6].upper() for d in donor_ids]
    mrows: list[list[Any]] = [header]
    order = {"supply": 0, "support_range": 1, "avoid_introduction": 2, "monitor": 3}
    for t in sorted(result["targets"], key=lambda x: (order.get(x["direction"], 9), x["label"])):
        rec_txt = (
            "not detected"
            if t["recipient"]["status"] == "not_detected"
            else (f"{_num(t['recipient']['percentile'], 0)}th" if t["recipient"]["percentile"] is not None else "—")
        )
        cells = []
        for did in donor_ids:
            cell = avail.get((t["target_id"], did))
            if cell is None:
                cells.append(f"<font color='{_hex(INK_FAINT)}'>ctx</font>")
            elif cell["supported_value"] is None:
                cells.append(f"<font color='{_hex(INK_FAINT)}'>n/a</font>")
            else:
                v = float(cell["supported_value"])
                col = GREEN if v >= 0.999 else (ACCENT_DARK if v > 0 else INK_FAINT)
                cells.append(f"<font color='{_hex(col)}'>{v:.2f}</font>")
        mrows.append(
            [
                t["label"][:58],
                t["origin"].replace("_supported", "").replace("_", " ")[:12],
                t["direction"].replace("_introduction", "").replace("support_", "")[:7],
                rec_txt,
                *cells,
            ]
        )
    widths = [58 * mm, 20 * mm, 14 * mm, 18 * mm] + [
        (CONTENT_WIDTH - 110 * mm) / max(1, len(donor_ids))
    ] * len(donor_ids)
    story.append(_table(mrows, widths, st))

    # ---- 5. concerns ------------------------------------------------------ #
    # No page break: the matrix routinely ends mid-page and a forced break
    # leaves a near-empty sheet. The heading is the section marker.
    story.append(_p("Transmission concerns and screening gaps", st["h2"]))
    story.append(
        _p(
            f"<font color='{_hex(AMBER)}'><b>No confirmed exclusion identified in the available evidence "
            "is not a negative pathogen screen.</b></font> Nothing below was confirmed by a validated "
            "laboratory assay on the material that would actually be used, no laboratory controls "
            "accompanied this run, and stool DNA sequencing cannot cover much of what a donor programme "
            "must test for.",
            st["small"],
        )
    )
    by_person: dict[str, list[dict[str, Any]]] = {}
    for c in result["concern_records"]:
        by_person.setdefault(names.get(c["candidate_id"], c["candidate_id"]), []).append(c)
    for person in sorted(by_person):
        rows_c = by_person[person]
        lead = [c for c in rows_c if c.get("tier") in LEAD_TIERS]
        rest = [c for c in rows_c if c.get("tier") not in LEAD_TIERS]
        exc = [c for c in rows_c if c["disposition"] == "exclusion_confirmed"]
        story.append(
            _p(
                f"<b>{person}</b> — {len(exc)} confirmed exclusion(s), {len(lead)} finding(s) needing "
                f"review, {len(rest)} technical-ambiguity or context rows",
                st["body"],
            )
        )
        crows: list[list[Any]] = [["FINDING", "WHAT IT MEANS", "EVIDENCE", "RULE"]]
        for c in sorted(lead, key=lambda x: (LEAD_TIERS.index(x["tier"]), x["feature_label"])):
            colour = AMBER if c["disposition"] != "context_only" else INK_FAINT
            m = c["measurements"]
            if "unique_supporting_fragments" in m:
                ev = (
                    f"{m.get('unique_supporting_fragments')} specific fragment(s), "
                    f"{m.get('informative_regions_supported')} region(s)"
                )
            elif "supporting_fragments" in m:
                ev = f"{m.get('supporting_fragments')} fragment(s)"
            else:
                ev = "—"
            crows.append(
                [
                    c["feature_label"][:40],
                    f"<font color='{_hex(colour)}'>{TIER_LABELS.get(c['tier'], c['tier'])}</font>",
                    ev,
                    c["rule_id"],
                ]
            )
        story.append(_table(crows, [50 * mm, 62 * mm, 34 * mm, 12 * mm], st))
        if rest:
            story.append(
                _p(
                    f"{len(rest)} further rows for {person} carry no organism-specific fragment, or are "
                    "background, decoy or disease-association context. They are retained in match.json; "
                    "none is evidence for or against any organism.",
                    st["small"],
                )
            )
        story.append(Spacer(1, 2.5 * mm))

    gaps = (result["screening_records"] or [{}])[0].get("not_covered_by_this_assay") or []
    if gaps:
        story.append(_p("What this assay cannot cover at all", st["h3"] if "h3" in st else st["body"]))
        for g in gaps:
            story.append(_p(f"· {g}", st["small"]))

    # ---- context deliberately outside the score --------------------------- #
    story.append(_p("Facts deliberately kept out of the score", st["h2"]))
    if result["individual_results"]:
        ctx0 = result["individual_results"][0]["context_not_in_the_score"]
        story.append(_p(ctx0["gmwi2"]["excluded_because"] + ".", st["small"]))
        story.append(_p(ctx0["diversity_excluded_because"] + ".", st["small"]))
        story.append(_p(ctx0["patterns_excluded_because"] + ".", st["small"]))
    ctxrows: list[list[Any]] = [["DONOR", "GMWI2", "BAND", "SHANNON", "SPECIES", "HIGHEST PATTERN"]]
    for row in sorted(result["individual_results"], key=lambda x: -x["goals"]["covered"]):
        c = row["context_not_in_the_score"]
        top = c["highest_disease_patterns"][0] if c["highest_disease_patterns"] else {}
        ctxrows.append(
            [
                row["person_id"],
                _num(c["gmwi2"]["score"], 2),
                str(c["gmwi2"]["band"] or "—"),
                _num(c["shannon_diversity"], 2),
                str(c["species_richness"] or "—"),
                f"{top.get('label', '—')} ({_num(top.get('percentile'), 0)}th)" if top else "—",
            ]
        )
    story.append(_table(ctxrows, [20 * mm, 18 * mm, 32 * mm, 20 * mm, 18 * mm, CONTENT_WIDTH - 108 * mm], st))

    # ---- 6-8. ecology, sensitivity, capabilities -------------------------- #
    story.append(_p("Ecology and strain resolution", st["h2"]))
    erows: list[list[Any]] = [["DONOR", "BRAY–CURTIS", "JACCARD", "SHARED SPP.", "DONOR-UNIQUE", "STRAINS"]]
    for row in result["individual_results"]:
        d = row["descriptors"]
        erows.append(
            [
                row["person_id"],
                _num(d["bray_curtis_dissimilarity"], 3),
                _num(d["jaccard_overlap"], 3),
                str(d["shared_species"]),
                str(d["donor_unique_species"]),
                d["strain_resolution"]["status"],
            ]
        )
    story.append(_table(erows, [26 * mm, 28 * mm, 24 * mm, 26 * mm, 30 * mm, 30 * mm], st))
    story.append(
        _p(
            "Distance carries no universal direction: whether similarity helps differs by cohort and "
            "feature type, so these descriptors are displayed and never weighted into coverage.",
            st["small"],
        )
    )

    story.append(_p("Sensitivity", st["h2"]))
    story.append(_p(result["sensitivity_results"]["meaning"], st["small"]))
    srows2: list[list[Any]] = [["SCENARIO", "SCORED GOALS", "LEADING CANDIDATES"]]
    for s in result["sensitivity_results"]["scenarios"]:
        leaders = ", ".join(
            "+".join(names.get(m, m) for m in cid.split("+")) for cid in s.get("leaders") or []
        )
        srows2.append([s["scenario_id"], str(s.get("scoring_target_count") or "—"), leaders or "—"])
    story.append(_table(srows2, [42 * mm, 26 * mm, CONTENT_WIDTH - 68 * mm], st))

    story.append(_p("Capabilities and limitations", st["h2"]))
    krows: list[list[Any]] = [["CAPABILITY", "STATUS", "WHAT IS MISSING"]]
    for c in result["capability_manifest"]:
        miss = "; ".join(c.get("missing_artifacts") or []) or (c.get("note") or "—")
        krows.append([c["capability_id"], c["status"].replace("_", " "), miss[:150]])
    story.append(_table(krows, [46 * mm, 28 * mm, CONTENT_WIDTH - 74 * mm], st))
    for line in result["limitations"]:
        story.append(_p(f"· {line}", st["small"]))
    story.append(
        _p(
            f"<font color='{_hex(INK_FAINT)}' size='6.6'>run {result['run_id']} · "
            f"created {result['created_at']} · every number here is also in match.json</font>",
            st["small"],
        )
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    doc.build(story)
    return path
