"""Pages 1 and 2: the summary dashboard, then what stood out and how to read it.

Page 1 is the report for people who read one page: a gut-health dial, the
headline findings, every metabolite pathway on the shared scale, the community
as a donut, and the disease-pattern resemblance bars. Page 2 says the same
things in sentences and explains the scale once, so nothing later needs a
legend.
"""

from __future__ import annotations

import io
import re
from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Flowable, KeepInFrame, Spacer, Table, TableStyle
from reportlab.platypus.flowables import _listWrapOn

from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    function_dest,
    linked_block,
    linked_cell,
    profile_dest,
    section_dest,
    section_heading,
)
from openbiota.pdfreport import (
    _PASTEL,
    ACCENT,
    AMBER,
    BODY_FRAME_HEIGHT,
    CONTENT_WIDTH,
    CORAL,
    FIRST_FRAME_HEIGHT,
    GREEN,
    INK,
    INK_FAINT,
    INK_SOFT,
    NOT_SCORED,
    ORANGE,
    PAGE,
    PHYLA_PALETTE,
    REPORT_TITLE,
    RULE,
    SECTIONS,
    SLATE,
    SLATE_BG,
    Dial,
    Donut,
    Dot,
    MetaboliteRow,
    PercentileBar,
    Status,
    StatusChip,
    Tile,
    _fit_text,
    _hex,
    _ordinal,
    notability,
    status_for,
)
from openbiota.shortnames import short_metabolite, short_profile_label

#: Page 1 shows at most this many function rows and profile rows; the rest
#: are counted and live in their own sections.
#: Ceilings on the rows searched when filling a page-1 column. The height
#: left in the column is what actually decides the count; these only bound
#: the search so a sample with thirty scored patterns does not try thirty
#: layouts.
#: Page one used to give almost all its room to these two lists, which meant
#: whole sections - the organism catalogue, biofilm, strains - had no way in
#: from the front of the report. Both lists continue in full in their own
#: sections, and a shorter list here buys the strip of links across the foot.
SUMMARY_MAX_FUNCTIONS = 10
SUMMARY_MAX_PROFILES = 11
SUMMARY_MIN_ROWS = 3
SUMMARY_MAX_HEADLINES = 8

#: Vertical gap between the blocks on page 1, and the gutter between its two
#: columns. One value each, so the page has one rhythm.
BLOCK_GAP: Final = 1.6 * mm
COLUMN_GUTTER: Final = 8 * mm

#: Life stages the biological age is described in. A number alone means
#: little to a reader; "reads like a community in its thirties" means
#: something immediately. Upper bound is exclusive.
AGE_BANDS: Final[tuple[tuple[float, str], ...]] = (
    (13.0, "a child's"),
    (20.0, "an adolescent's"),
    (30.0, "a young adult's"),
    (40.0, "an adult in their thirties"),
    (50.0, "an adult in their forties"),
    (60.0, "an adult in their fifties"),
    (75.0, "an older adult's"),
    (200.0, "an elderly adult's"),
)


def age_band(years: float) -> str:
    """Plain-language life stage for a biological age in years."""
    for ceiling, words in AGE_BANDS:
        if years < ceiling:
            return words
    return AGE_BANDS[-1][1]


def age_headline(age: Any) -> str:
    """The biological-age sentence, built from the reported ensemble figure.

    This reads `reported_age_years` — the ensemble estimate after the
    community-health adjustment — because that is the number on the page and
    in the tile. An earlier version quoted the raw model point instead and
    told readers their age "could not be read", which contradicted the
    figure printed beside it and undersold a result the pipeline had
    actually computed.
    """
    reported = getattr(age, "reported_age_years", None)
    if getattr(age, "status", None) != "scored" or reported is None:
        return f"<b>Biological age of biota: {getattr(age, 'status_line', 'not computed')}.</b>"

    adj = getattr(age, "health_adjustment", None) or {}
    band = age_band(float(reported))
    lead = (
        f"<b>Your gut community reads like {band}</b> "
        f"&mdash; <b>{reported:.0f} years</b>"
    )

    # What moved it, in the reader's terms rather than the model's.
    drivers: list[str] = []
    adverse = adj.get("profiles_adverse_at_or_above_90th") or 0
    typical = adj.get("profiles_typical_below_79th") or 0
    gm = adj.get("gmwi2_score")
    if adverse:
        drivers.append(
            f"{adverse} disease pattern{'s' if adverse != 1 else ''} scoring high"
        )
    if typical:
        drivers.append(
            f"{typical} pattern{'s' if typical != 1 else ''} in the typical range or better"
        )
    if isinstance(gm, (int, float)):
        drivers.append(
            "a "
            + ("healthy" if gm > 0.25 else "disturbed" if gm < -0.25 else "middling")
            + f" general gut-health index ({gm:+.2f})"
        )

    if drivers:
        detail = (
            " Built by an ensemble of a transformer neural network and a linear "
            "carriage model reading species composition, then moved by "
            + ", ".join(drivers)
            + "."
        )
    else:
        detail = (
            " Estimated by an ensemble of a transformer neural network and a "
            "linear carriage model reading species composition."
        )

    return (
        lead + "." + detail +
        " Communities carrying fewer disease-associated patterns and a healthier "
        "index read younger; this is an estimate of the community, not of you."
    )

# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _table(rows: list[list[Any]], widths: Sequence[float], extra: Sequence[Any] = ()) -> Table:
    table = LinkedTable(rows, colWidths=list(widths), hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 2.0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.0),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                *extra,
            ]
        )
    )
    return table


def profile_status(result: Any) -> Status:
    """A disease pattern's reading, on the same scale as everything else.

    Higher resemblance is the adverse direction, so the colours read the same
    way as an adverse-direction metabolite.
    """
    if result.abstention.abstained or result.combined_percentile is None:
        return NOT_SCORED
    return status_for(percentile=result.combined_percentile, higher_means="adverse")


def anchor_words(anchor: Any) -> tuple[str, colors.Color]:
    """Plain words for the GMWI2 band, and the colour that goes with them."""
    if anchor is None or not anchor.valid:
        return "not computed", SLATE
    if anchor.strongly_negative:
        return "broadly disturbed", CORAL
    if anchor.band == "mildly dysbiotic":
        return "mildly disturbed", ORANGE
    if anchor.band == "indeterminate":
        return "indeterminate", AMBER
    return "healthy range", GREEN


def organism_headlines(
    findings: Any | None, inv: Any = None, *, limit: int = 3,
) -> list[tuple[Any, str, str]]:
    """What is overgrown, missing or depleted — ranked by consequence, not by category.

    This used to return Missing, Depleted, Expanded in a fixed order and the
    caller took the first two, so an overgrowth was cut whenever anything
    was missing — which is nearly always. A 13% bloom of a Crohn's-associated
    mucus degrader was losing its place to the third example of a missing
    species. Each candidate now carries a weight and the heaviest lead.

    Overgrowths come from the classified inventory, ranked by
    :func:`openbiota.organisms.importance` — percentile excess times
    abundance, with an opportunist worth more than a conditional organism.
    Missing organisms cannot be in an inventory, so they still come from the
    census. Each item links to the organism section, where the full lists
    and the descriptions live.
    """
    from openbiota import organisms as org

    # Every headline here lands on the organisms-that-need-attention section,
    # which lists the overgrown, the missing and the depleted, each linked on
    # to its own page. The full catalogue is the wrong landing: it is every
    # organism, and the reader clicked on a problem.
    where = section_dest("organisms")
    candidates: list[tuple[float, tuple[Any, str, str]]] = []

    def names(ts: list[Any], n: int) -> str:
        shown = ", ".join(f"<i>{t.display_name}</i>" for t in ts[:n])
        rest = len(ts) - n
        return shown + (f" and {rest} more" if rest > 0 else "")

    vs = org.verdicts(list(inv.organisms)) if inv is not None and len(inv) else []
    flagged = org.issues(vs)
    over = [v for v in flagged if v.flag == "high" and v.cls in (org.OPPORTUNIST, org.CONDITIONAL)]
    low_v = [v for v in flagged if v.flag == "low" and v.cls == org.BENEFICIAL]

    if over:
        lead = over[0]
        acute = lead.cls == org.OPPORTUNIST and lead.is_issue
        # Name the two that matter most; count only the rest that are
        # issues, not watches, or the headline reads as two dozen overgrowths.
        serious = [v for v in over if v.is_issue]
        # The share the rest of the report prints, and the level against the
        # typical carrier: one number per quantity, page one included.
        shown = ", ".join(
            f"<i>{v.organism.display}</i> ({_pct_short(v.organism.percent) if v.organism.in_primary else 'detected'}, "
            + (f"{_deviation(v.organism)} vs typical carrier" if _deviation(v.organism) else
               f"{org._ordinal(v.organism.level_percentile or 0)} pct") + ")"
            for v in over[:2]
        )
        rest = len(serious) - 2
        text = (
            f"<b>Overgrown</b> — {shown}"
            + (f" and {rest} more at a level of concern" if rest > 0 else "")
            + (f". <i>{lead.organism.display}</i> is an {lead.label.lower()}: normal at low "
               "levels, associated with inflammation when it expands."
               if acute else ", above the reference range.")
        )
        candidates.append((org.importance(lead), (CORAL if acute else AMBER, text, where)))

    if findings is not None:
        missing = list(findings.bucket(4))
        if missing:
            candidates.append((
                55.0 + 12.0 * len(missing),
                (CORAL, f"<b>Missing</b> — {len(missing)} species most reference adults carry were "
                        f"not detected: {names(missing, 3)}.", where),
            ))
        concern = list(findings.bucket(6))
        if concern:
            candidates.append((
                25.0,
                (SLATE, f"<b>Needs qualification</b> — {names(concern, 2)} detected; species DNA is "
                        "not a strain, toxin or infection.", where),
            ))

    if low_v:
        lead = low_v[0]
        shown = ", ".join(f"<i>{v.organism.display}</i>" for v in low_v[:2])
        rest = len(low_v) - 2
        candidates.append((
            org.importance(lead),
            (ORANGE, f"<b>Depleted</b> — {shown}" + (f" and {rest} more" if rest > 0 else "")
                     + " well below the reference carriers.", where),
        ))
    elif findings is not None and list(findings.bucket(3)):
        low = list(findings.bucket(3))
        candidates.append((
            40.0,
            (ORANGE, f"<b>Depleted</b> — {names(low, 2)} well below the reference carriers.", where),
        ))

    candidates.sort(key=lambda c: -c[0])
    return [item for _, item in candidates[:limit]]


def _pct_short(v: float) -> str:
    return f"{v:.1f}%" if v >= 1 else f"{v:.2f}%"


def _deviation(o: Any) -> str:
    """``+400%`` / ``×110`` against the typical carrier, or '' when there is no reference."""
    from openbiota.pdflibrary import deviation_text

    dev = getattr(o, "deviation_percent", None)
    return deviation_text(dev) if dev is not None else ""


def headline_items(
    rows: Sequence[MetaboliteRow], similarity: Any | None, findings: Any | None = None,
    limit: int = SUMMARY_MAX_HEADLINES, inv: Any = None,
) -> list[tuple[Any, str, str]]:
    """The six or seven things worth knowing: colour, text, and where each links.

    Item-level findings (a metabolite, a disease pattern) link to their own
    detail card; section-level ones (gut health, the organism lists) to the
    section that explains them.
    """
    items: list[tuple[Any, str, str]] = []

    anchor = similarity.anchor if similarity is not None else None
    if anchor is not None and anchor.valid:
        words, colour = anchor_words(anchor)
        items.append(
            (colour, f"<b>General gut health</b> — {words} (index {anchor.score:+.2f}).",
             section_dest("community"))
        )
    # Three organism slots, ranked by consequence. Two was one too few: the
    # overgrowth was always the one to go.
    items.extend(organism_headlines(findings, inv, limit=3))

    judged = [r for r in rows if r.status.colour in (CORAL, ORANGE, AMBER) and r.percentile is not None]
    judged.sort(key=lambda r: notability(r.status, r.percentile))
    for row in judged[:2]:
        items.append(
            (
                row.status.colour,
                f"<b>{short_metabolite(row.metabolite)}</b> — {row.status.label}, "
                f"{_ordinal(row.percentile)} percentile.",
                function_dest(row.panel, detail=True),
            )
        )

    good = [
        r for r in rows
        if r.status.colour is GREEN and r.percentile is not None and abs(r.status.level) >= 1
    ]
    good.sort(key=lambda r: -abs(r.status.level))
    for row in good[:1]:
        items.append(
            (
                GREEN,
                f"<b>{short_metabolite(row.metabolite)}</b> — {row.status.label} in the favourable direction, "
                f"{_ordinal(row.percentile)} percentile.",
                function_dest(row.panel, detail=True),
            )
        )
    if not judged and not good:
        # "Every pathway is typical" may only be said about pathways that
        # were actually placed. Provisional, not-detected and no-reference
        # rows all wear the neutral colour, and counting them as typical
        # turned "we could not judge this" into "this is fine".
        assessed = [r for r in rows if r.status.assessed and r.percentile is not None]
        held_back = len(rows) - len(assessed)
        if assessed:
            items.append((
                GREEN,
                f"<b>{len(assessed)} of {len(rows)} metabolite pathways</b> were "
                "placed against the reference and all came out typical"
                + (
                    f"; the other {held_back} could not be placed and are "
                    "neither typical nor atypical."
                    if held_back else "."
                ),
                section_dest("functions"),
            ))
        else:
            items.append((
                SLATE,
                f"<b>None of the {len(rows)} metabolite pathways</b> could be "
                "placed against the reference, so none is reported as typical.",
                section_dest("functions"),
            ))

    if similarity is not None and similarity.results:
        scored = [
            r for r in similarity.results
            if not r.abstention.abstained and r.combined_percentile is not None
        ]
        high = [r for r in scored if profile_status(r).level >= 1]
        if high:
            top = max(high, key=lambda r: r.combined_percentile)
            status = profile_status(top)
            items.append(
                (
                    status.colour,
                    f"<b>{short_profile_label(top.profile.name, top.profile.label)}</b> — "
                    f"{status.label} resemblance, "
                    f"{_ordinal(top.combined_percentile)} percentile.",
                    profile_dest(top.profile.name, detail=True),
                )
            )
        elif scored:
            items.append(
                (
                    GREEN,
                    f"<b>Disease patterns</b> — none of the {len(scored)} compared showed "
                    "above-average resemblance.",
                    section_dest("patterns"),
                )
            )
        discordant = [
            c for r in similarity.results for c in r.cross_engine if c.verdict == "DISCORDANT"
        ]
        if discordant:
            items.append(
                (
                    SLATE,
                    "<b>Two independent measurements disagree</b> on one point — "
                    f"see section {SECTIONS['patterns_detail']}.",
                    section_dest("patterns_detail"),
                )
            )
    return items[:limit]


# --------------------------------------------------------------------------- #
# page 1 — the dashboard
# --------------------------------------------------------------------------- #


def summary_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    rows: Sequence[MetaboliteRow],
    similarity: Any | None,
    rpob_profile: dict[str, Any] | None,
    meta: dict[str, Any],
    findings: Any | None = None,
    age: Any | None = None,
    age_meta: dict[str, Any] | None = None,
    pathogens: Mapping[str, Any] | None = None,
    pathogen_section: int = 0,
    results: Mapping[str, Any] | None = None,
    inv: Any = None,
) -> None:
    """Page 1 — every measurement, graphics first.

    The page is a title block, the gut-health dial beside what stood out, two
    columns, and the pathogen verdict across the foot as the page's anchor.
    The left column is biological age over disease signatures; the right is
    diversity over metabolite production. The columns flow independently, so
    the short age block does not leave a hole under itself the way it did when
    the two rows were locked to each other's height.

    Each column is filled to the height it actually has. The fixed parts of
    the page are measured first; whatever is left is handed to the two list
    blocks, and each takes as many rows as fit. That is what makes the page
    both full and one page: the count of patterns shown is decided by the
    room, not by a constant, so a sample with a long headline list simply
    shows one row fewer rather than spilling or shrinking.

    The frame stays `mode="shrink"` as a backstop only. It scales uniformly,
    so vertical overflow also pulls the right edge in and leaves a dead strip
    down the page; the fill loop exists so that never has to happen.
    """
    # Measure exactly as `KeepInFrame` will: same `_listWrapOn`, against a
    # real canvas. Wrapping without one reports paragraphs materially taller.
    probe = Canvas(io.BytesIO(), pagesize=PAGE)
    col = (CONTENT_WIDTH - COLUMN_GUTTER) / 2

    def height(flowables: Sequence[Any], width: float) -> float:
        return _listWrapOn(list(flowables), width, probe)[1] if flowables else 0.0

    def scale_of(page: list[Any]) -> float:
        # Ask the frame itself rather than re-deriving its arithmetic; it is
        # the only thing whose opinion decides whether the page gets scaled.
        frame = KeepInFrame(
            CONTENT_WIDTH, FIRST_FRAME_HEIGHT, page, mode="shrink", hAlign="LEFT"
        )
        frame.canv = probe
        frame.wrap(CONTENT_WIDTH, FIRST_FRAME_HEIGHT)
        return float(getattr(frame, "_scale", 1.0) or 1.0)

    def profiles(n: int) -> list[Any]:
        block: list[Any] = []
        _profile_block(block, st, similarity=similarity, limit=n, width=col, chips=False)
        return block

    def functions(n: int) -> list[Any]:
        block: list[Any] = []
        _metabolite_block(block, st, rows=rows, limit=n, width=col, chips=False)
        return block

    def fill(make: Any, fixed_height: float, avail: float, ceiling: int) -> int:
        """Most rows whose block, under its fixed partner, fits in `avail`."""
        for n in range(ceiling, SUMMARY_MIN_ROWS, -1):
            if fixed_height + BLOCK_GAP + height(make(n), col) <= avail:
                return n
        return SUMMARY_MIN_ROWS

    # Parts that do not change with the row budget, measured once.
    foot = _foot(st, pathogens=pathogens, pathogen_section=pathogen_section)
    h_foot = height(foot, CONTENT_WIDTH) + (BLOCK_GAP if foot else 0.0)
    # Section-level elements: the whole block is one link to its section.
    right_top: list[Any] = []
    _community_block(
        right_top, st, similarity=similarity, rpob_profile=rpob_profile,
        findings=findings, width=col, inv=inv,
        primary_phyla=(results or {}).get("extended_catalogue", {}).get("phyla") if results else None,
    )
    right_top = [linked_block(section_dest("community"), right_top, col)]
    h_right_top = height(right_top, col)

    left_top: list[Any] = []
    _age_block(left_top, st, age=age, age_meta=age_meta, width=col)
    left_top = [linked_block(section_dest("age"), left_top, col)]
    # Biofilm sits between the age block and the disease list, in the left
    # column, so the reader meets it on page one rather than twenty pages in.
    bio: list[Any] = []
    _biofilm_block(bio, st, results=results, width=col)
    if bio:
        left_top = [
            *left_top, Spacer(1, BLOCK_GAP),
            linked_block(section_dest("biofilm"), bio, col),
        ]
    # The mycobiome score sits directly under the biofilm score: same grammar,
    # one row, no fine print - the section page carries the explanation.
    myc: list[Any] = []
    _mycobiome_block(myc, st, results=results, width=col)
    if myc:
        left_top = [
            *left_top, Spacer(1, BLOCK_GAP),
            linked_block(section_dest("mycobiome"), myc, col),
        ]
    # And the gut-skin axis directly under it, same grammar again.
    skn: list[Any] = []
    _skin_block(skn, st, results=results, width=col)
    if skn:
        left_top = [
            *left_top, Spacer(1, BLOCK_GAP),
            linked_block(section_dest("skin"), skn, col),
        ]
    h_left_top = height(left_top, col)

    # The age block is shorter than the diversity block beside it, so the
    # disease list starts higher than the metabolite list and shows more
    # rows — that is the point of keeping the age block compact. For the two
    # lists to read as one grid, the head start has to be a whole number of
    # rows: pad under the age block by the remainder, and both lists' rows
    # then sit on the same lines wherever they overlap.
    pitch = height(profiles(2), col) - height(profiles(1), col)
    head_start = h_right_top - h_left_top
    left_pad = 0.0
    rows_ahead = 0
    if pitch > 0 and head_start > 0:
        rows_ahead = int(head_start // pitch)
        left_pad = head_start - rows_ahead * pitch
        if left_pad > 0.35 * pitch:
            # Nearer the next line than this one: pad up to it instead of
            # leaving most of a row's worth of air under the age block.
            left_pad = head_start - (rows_ahead + 1) * pitch
            rows_ahead += 1
        if left_pad < 0:
            # Padding cannot be negative; pull the head start back a row.
            rows_ahead -= 1
            left_pad = head_start - rows_ahead * pitch
    h_left_top += left_pad

    n_profiles_max = min(
        SUMMARY_MAX_PROFILES,
        len([r for r in similarity.ranked if r.reportable]) if similarity is not None else 0,
    )
    n_functions_max = min(SUMMARY_MAX_FUNCTIONS, len(rows))

    page: list[Any] | None = None
    for n_headlines in range(SUMMARY_MAX_HEADLINES, 3, -1):
        head = _head(st, rows=rows, similarity=similarity, findings=findings,
                     meta=meta, n_headlines=n_headlines, inv=inv)
        avail = FIRST_FRAME_HEIGHT - height(head, CONTENT_WIDTH) - BLOCK_GAP - h_foot
        n_pr = fill(profiles, h_left_top, avail, max(n_profiles_max, SUMMARY_MIN_ROWS))
        n_fn = fill(functions, h_right_top, avail, max(n_functions_max, SUMMARY_MIN_ROWS))
        # The disease list starts `rows_ahead` rows earlier; give it exactly
        # that many more rows and the two lists end on the same line, so the
        # foot sits under both at one even gap.
        if n_profiles_max >= SUMMARY_MIN_ROWS + rows_ahead:
            n_fn = max(SUMMARY_MIN_ROWS, min(n_fn, n_pr - rows_ahead))
            n_pr = n_fn + rows_ahead

        # Measurement of the parts and of the whole can disagree by a point
        # or two, so confirm with the frame and give back a row if needed.
        for _ in range(3):
            candidate = _assemble(
                head, left_top, profiles(n_pr), right_top, functions(n_fn), foot, col,
                left_pad=left_pad,
            )
            if scale_of(candidate) <= 1.0 + 1e-6:
                page = candidate
                break
            if n_pr <= SUMMARY_MIN_ROWS and n_fn <= SUMMARY_MIN_ROWS:
                break
            # Take a row from each list, so they keep ending on the same line.
            n_fn -= 1
            n_pr -= 1
        if page is not None:
            break

    if page is None:
        # Nothing fitted even at the tightest budget. Take the tightest and
        # let the frame scale it; a slightly small page beats a lost row.
        head = _head(st, rows=rows, similarity=similarity, findings=findings,
                     meta=meta, n_headlines=4, inv=inv)
        page = _assemble(
            head, left_top, profiles(SUMMARY_MIN_ROWS), right_top,
            functions(SUMMARY_MIN_ROWS), foot, col, left_pad=left_pad,
        )

    story.append(
        KeepInFrame(CONTENT_WIDTH, FIRST_FRAME_HEIGHT, page, mode="shrink", hAlign="LEFT")
    )


def _head(
    st: dict[str, ParagraphStyle],
    *,
    rows: Sequence[MetaboliteRow],
    similarity: Any | None,
    findings: Any | None,
    meta: dict[str, Any],
    n_headlines: int,
    inv: Any = None,
) -> list[Any]:
    """Title block and the hero row: everything above the two columns."""
    # The frame starts just under the wordmark; this puts the kicker at the
    # foot of the header band the wordmark is centred in.
    head: list[Any] = [
        Anchor(section_dest("summary"), above=20 * mm, outline=(f"{SECTIONS['summary']}  Summary", 0)),
        Spacer(1, 1.9 * mm),
    ]
    head.append(
        Paragraph(
            "SHOTGUN METAGENOMICS &nbsp;·&nbsp; STOOL DNA &nbsp;·&nbsp; "
            f"{str(meta.get('read_pairs', '—')).upper()} READ PAIRS, EVERY ONE ANALYSED",
            st["kicker"],
        )
    )
    # The wordmark lives in the page header, drawn on the canvas opposite the
    # sample line, so the title block is the title and nothing else.
    head.append(Paragraph(REPORT_TITLE, st["title"]))
    head.append(
        Paragraph(
            "What your gut bacteria can produce, who is living there, how healthy the "
            "community looks, and which published disease patterns it resembles — read "
            "straight from the DNA in your sample.",
            st["tagline"],
        )
    )
    head.append(Spacer(1, 4.5 * mm))
    _hero(head, st, rows=rows, similarity=similarity, findings=findings,
          n_headlines=n_headlines, inv=inv)
    return head


def _foot(
    st: dict[str, ParagraphStyle],
    *,
    pathogens: Mapping[str, Any] | None,
    pathogen_section: int,
) -> list[Any]:
    """The pathogen verdict, or nothing if the pathogen branch did not run."""
    from openbiota import pdfpathogens

    foot: list[Any] = []
    pdfpathogens.pathogen_alert(foot, st, pathogens=pathogens, section=pathogen_section)
    return foot


#: Bands for the concerning biofilm axis: (upper bound, track pastel, marker
#: colour). Lower is better on this axis, and only a genuinely low reading
#: counts as good. The track is painted in the same pastels every other
#: scale bar uses, so the marker is the one saturated thing on the row.
_CONCERN_BANDS: Final[tuple[tuple[float, colors.Color, colors.Color], ...]] = (
    (20.0, _PASTEL["good"], GREEN),
    (50.0, _PASTEL["amber"], AMBER),
    (70.0, _PASTEL["orange"], ORANGE),
    (85.0, colors.HexColor("#F0C4A8"), colors.HexColor("#D0641E")),  # deep orange, pastel + solid
    (100.0, _PASTEL["coral"], CORAL),
)


def _concern_colour(value: float) -> colors.Color:
    """The saturated colour for the band a reading falls in."""
    for upper, _track, marker in _CONCERN_BANDS:
        if value < upper or upper >= 100.0:
            return marker
    return _CONCERN_BANDS[-1][2]


class GradedBar(Flowable):
    """A bar whose track is painted in the bands the reading is judged on.

    The ordinary percentile bar is symmetric and pastel because most
    readings are two-sided. This one is for a one-sided axis where the
    reader should see at a glance which band the marker sits in.
    """

    def __init__(
        self,
        *,
        width: float,
        value: float,
        bands: tuple[tuple[float, colors.Color, colors.Color], ...],
        marker_colour: colors.Color,
        height: float = 5.0 * mm,
    ) -> None:
        super().__init__()
        self.width = width
        self.height = height
        self.value = max(0.0, min(100.0, value))
        self.bands = bands
        self.marker_colour = marker_colour

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        w = self.width
        t = min(3.4 * mm, self.height * 0.56)
        y = (self.height - t) / 2.0
        c.saveState()
        clip = c.beginPath()
        clip.roundRect(0, y, w, t, t / 2.0)
        c.clipPath(clip, stroke=0, fill=0)
        lower = 0.0
        for upper, track, _marker in self.bands:
            c.setFillColor(track)
            c.rect(w * lower / 100.0, y, w * (upper - lower) / 100.0 + 0.3, t, stroke=0, fill=1)
            lower = upper
        c.restoreState()
        x = w * self.value / 100.0
        r = t * 0.78
        c.setFillColor(colors.white)
        c.circle(x, y + t / 2.0, r + 0.9, stroke=0, fill=1)
        c.setFillColor(self.marker_colour)
        c.circle(x, y + t / 2.0, r, stroke=0, fill=1)


def _biofilm_block(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    results: Mapping[str, Any] | None,
    width: float = CONTENT_WIDTH,
) -> None:
    """The two biofilm readings, compactly, between age and disease patterns.

    Two numbers that mean opposite things, so they are never merged into one
    score: a community that builds biofilm readily and a community that
    maintains the gut lining are separate facts and can both be high. Each
    is shown out of 100 against the same reference, coloured by whether its
    own direction is the unwanted one.
    """
    from openbiota.pdfreport import INK_FAINT, PercentileBar, status_for

    cards = [
        c for c in ((results or {}).get("biofilm") or {}).get("cards") or []
        if c.get("reference_percentile") is not None
    ]
    if not cards:
        return

    _block_heading(story, st, "Biofilm Score", f"section {SECTIONS['biofilm']}",
                   width, href=section_dest("biofilm"))

    rows: list[list[Any]] = []
    for card in cards:
        pct = float(card["reference_percentile"])
        concerning = "Concerning" in str(card.get("heading") or "")
        if concerning:
            # Lower is better, and only a genuinely low reading is good:
            # under 20 green, to 50 yellow, to 70 orange, to 85 deep
            # orange, then red. The number takes the band's colour so it
            # cannot be read against the bar.
            colour = _concern_colour(pct)
            bar: Any = GradedBar(width=width * 0.30, value=pct, height=5.0 * mm,
                                 bands=_CONCERN_BANDS, marker_colour=colour)
            label = "Biofilm-forming <font color='#8C99A6'>(bad)</font>"
        else:
            status = status_for(percentile=pct, higher_means="favourable")
            colour = status.colour
            bar = PercentileBar(width=width * 0.30, percentile=pct,
                                higher_means="favourable",
                                marker_colour=colour, height=5.0 * mm)
            label = "Protective gut-lining <font color='#8C99A6'>(good)</font>"
        rows.append([
            Paragraph(label, st["cell"]),
            bar,
            Paragraph(
                f"<font color='{_hex(colour)}'><b>{pct:.0f}</b></font>"
                f"<font size='6' color='{_hex(INK_FAINT)}'>/100</font>",
                st["pct"]),
        ])

    table = Table(rows, colWidths=[width * 0.44, width * 0.32, width * 0.24],
                  hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("LEFTPADDING", (1, 0), (1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 1.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.0),
    ]))
    story.append(table)


def _mycobiome_block(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    results: Mapping[str, Any] | None,
    width: float = CONTENT_WIDTH,
) -> None:
    """One row: the Myco-Score, in the same grammar as every other row on
    the page - label, the standard percentile bar (higher is healthier),
    number out of 100. No sub-text; the section page explains it."""
    from openbiota.pdfreport import INK_FAINT, PercentileBar, status_for

    myco = (results or {}).get("mycobiome") or {}
    ms = myco.get("myco_score") or {}
    if not myco or myco.get("analysis_status") in (None, "not_assessed") or ms.get("value") is None:
        return
    _block_heading(story, st, "Mycobiome", f"section {SECTIONS['mycobiome']}", width, href=section_dest("mycobiome"))
    value = float(ms["value"])
    colour = status_for(percentile=value, higher_means="favourable").colour
    table = Table([[
        Paragraph("Myco-Score <font color='#8C99A6'>(experimental)</font>", st["cell"]),
        PercentileBar(width=width * 0.30, percentile=value, higher_means="favourable", marker_colour=colour, height=5.0 * mm),
        Paragraph(f"<font color='{_hex(colour)}'><b>{int(value + 0.5)}</b></font>"
                  f"<font size='6' color='{_hex(INK_FAINT)}'>/100</font>", st["pct"]),
    ]], colWidths=[width * 0.44, width * 0.32, width * 0.24], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("LEFTPADDING", (1, 0), (1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 1.0), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.0),
    ]))
    story.append(table)


def _skin_block(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    results: Mapping[str, Any] | None,
    width: float = CONTENT_WIDTH,
) -> None:
    """One row: the gut-skin axis, in the same grammar as the Myco-Score
    above it - label, the standard percentile bar with higher healthier,
    number out of 100. The section page carries the explanation."""
    from openbiota import skinaxis
    from openbiota.pdfreport import INK_FAINT, PercentileBar, status_for

    axis = skinaxis.score((results or {}).get("profile_similarity") if results else None)
    if not axis.scored:
        return
    _block_heading(story, st, "Skin", f"section {SECTIONS['skin']}", width,
                   href=section_dest("skin"))
    value = float(axis.value)
    colour = status_for(percentile=value, higher_means="favourable").colour
    table = Table([[
        Paragraph("Gut-skin axis <font color='#8C99A6'>(beta)</font>", st["cell"]),
        PercentileBar(width=width * 0.30, percentile=value, higher_means="favourable",
                      marker_colour=colour, height=5.0 * mm),
        Paragraph(f"<font color='{_hex(colour)}'><b>{int(value + 0.5)}</b></font>"
                  f"<font size='6' color='{_hex(INK_FAINT)}'>/100</font>", st["pct"]),
    ]], colWidths=[width * 0.44, width * 0.32, width * 0.24], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("LEFTPADDING", (1, 0), (1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 1.0), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.0),
    ]))
    story.append(table)


def _assemble(
    head: list[Any],
    left_top: list[Any],
    left_bottom: list[Any],
    right_top: list[Any],
    right_bottom: list[Any],
    foot: list[Any],
    col: float,
    *,
    left_pad: float = 0.0,
) -> list[Any]:
    """Put the measured parts together in page order.

    `left_pad` is the little extra under the age block that puts the disease
    rows on the same lines as the metabolite rows opposite.
    """
    left = [*left_top, Spacer(1, BLOCK_GAP + left_pad), *left_bottom]
    right = [*right_top, Spacer(1, BLOCK_GAP), *right_bottom]
    page: list[Any] = [*head, Spacer(1, BLOCK_GAP), _two_up(left, right, col, COLUMN_GUTTER)]
    if foot:
        page.extend([Spacer(1, BLOCK_GAP), *foot])
    return page


def _two_up(left: list[Any], right: list[Any], col: float, gap: float) -> Table:
    """Put two blocks side by side, each in its own column.

    Both columns are top-aligned so the two headings share a baseline, and
    the table carries no padding of its own: the gap between them is a real
    column, which keeps the outer edges exactly on the content margins.
    """
    table = Table([[left, right]], colWidths=[col + gap, col], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (0, 0), gap),
                ("RIGHTPADDING", (1, 0), (1, 0), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    return table


def _block_heading(
    story: list[Any], st: dict[str, ParagraphStyle], title: str, where: str,
    width: float = CONTENT_WIDTH,
    href: str | None = None,
) -> None:
    """A block's title with its section reference; the whole row is a link when `href` is given."""
    title_cell: Any = Paragraph(title, st["h2_tight"])
    if href:
        title_cell = linked_cell(RowLink(href=href), title_cell)
    head = LinkedTable(
        [[title_cell, Paragraph(where, st["fine_right"])]],
        colWidths=[width * 0.7, width * 0.3],
        hAlign="LEFT",
    )
    head.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "BOTTOM"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5),
                ("LINEBELOW", (0, 0), (-1, 0), 0.7, ACCENT),
            ]
        )
    )
    story.append(head)
    story.append(Spacer(1, 2 * mm))


def _hero(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    rows: Sequence[MetaboliteRow],
    similarity: Any | None,
    findings: Any | None = None,
    n_headlines: int = SUMMARY_MAX_HEADLINES,
    inv: Any = None,
) -> None:
    anchor = similarity.anchor if similarity is not None else None
    score = anchor.score if anchor is not None and anchor.valid else None
    words, _ = anchor_words(anchor)
    sub = (
        f"{anchor.n_health_taxa_present} healthy · {anchor.n_disease_taxa_present} disease markers"
        if anchor is not None and anchor.valid
        else "species engine did not run"
    )
    dial = Dial(value=score, label=words, sublabel=sub, caption="GMWI2 gut-health index", width=58 * mm)

    # Column arithmetic, stated once. The right-hand cell carries its own
    # gutter as padding, so the table inside it gets the cell width *minus*
    # that gutter. Sizing the inner table against the raw cell width made it
    # overflow by a few points, and ReportLab responded by shrinking every
    # column proportionally — which is what left a dead strip down the right
    # of the page while the blocks below ran full width.
    dial_col = 64 * mm
    gutter = 6 * mm
    bullet_col = 4.5 * mm
    right_col = CONTENT_WIDTH - dial_col
    text_col = right_col - gutter - bullet_col

    items = headline_items(rows, similarity, findings, limit=n_headlines, inv=inv)
    right: list[list[Any]] = [[Paragraph("WHAT STOOD OUT", st["label"]), ""]]
    for colour, text, href in items:
        right.append([linked_cell(RowLink(href=href), Dot(colour)), Paragraph(text, st["headline"])])
    right_table = LinkedTable(right, colWidths=[bullet_col, text_col], hAlign="LEFT")
    right_table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("SPAN", (0, 0), (1, 0)),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 1.7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.7),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1.0),
                ("LINEBELOW", (0, 1), (-1, -2), 0.3, RULE),
            ]
        )
    )

    # The dial is explained where the gut-health index is: the community section.
    hero = Table(
        [[linked_block(section_dest("community"), [dial], dial.width), right_table]],
        colWidths=[dial_col, right_col], hAlign="LEFT",
    )
    hero.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (0, 0), "MIDDLE"),
                ("VALIGN", (1, 0), (1, 0), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                ("LINEAFTER", (0, 0), (0, 0), 0.5, RULE),
                ("LEFTPADDING", (1, 0), (1, 0), gutter),
            ]
        )
    )
    story.append(hero)


_SCALE_COLS = (0.27, 0.44, 0.09, 0.20)
#: The same row at half width, without the chip. The name gets the chip's
#: share and more: every page-1 label has a short form, and the short form
#: is only worth having if the column it sits in is wide enough to hold it.
_SCALE_COLS_NARROW = (0.46, 0.38, 0.16)


def name_column_width(width: float, *, chips: bool) -> float:
    """Usable width for a row label, after the cell's own padding."""
    frac = (_SCALE_COLS if chips else _SCALE_COLS_NARROW)[0]
    return width * frac - 2 * _CELL_PAD


#: ReportLab's default cell padding, which `_table` leaves in place.
_CELL_PAD = 6.0

#: The direction arrow drawn after a function name: a space in the name's
#: font, then the glyph at its own smaller size.
_ARROW_SIZE = 6.0


def _arrow_allowance() -> float:
    from openbiota.shortnames import PAGE_ONE_FONT, PAGE_ONE_SIZE, text_width

    return text_width(" ", font=PAGE_ONE_FONT, size=PAGE_ONE_SIZE) + text_width(
        "\u25b2", font="Helvetica", size=_ARROW_SIZE
    )


def page_one_name_budget(*, arrow: bool) -> float:
    """Points a row name may take on page 1 and still sit on one line.

    Page 1 draws its function and disease rows at half width without a
    chip, so this is the name column of that layout after cell padding,
    less the direction arrow when the row carries one. The short-name table
    and its test measure against this number, not a character count.
    """
    col = (CONTENT_WIDTH - COLUMN_GUTTER) / 2
    budget = name_column_width(col, chips=False)
    return budget - _arrow_allowance() if arrow else budget


_TAGS = re.compile(r"<[^>]+>")


def _one_line(markup: str, avail: float) -> str:
    """A row name that cannot wrap.

    The short-name table is meant to make this a no-op, and the test on it
    is what keeps that true. If a name reaches here over budget anyway, it
    is drawn at the size that fits rather than on a second line: page 1 is
    a grid, and a two-line row is the one thing it must never show.
    """
    from openbiota.shortnames import PAGE_ONE_FONT, PAGE_ONE_MARGIN, PAGE_ONE_SIZE, text_width

    plain = _TAGS.sub("", markup)
    has_arrow = "\u25b2" in plain or "\u25bc" in plain
    name = plain.replace("\u25b2", "").replace("\u25bc", "").rstrip()
    width = text_width(name, font=PAGE_ONE_FONT, size=PAGE_ONE_SIZE)
    room = avail - PAGE_ONE_MARGIN - (_arrow_allowance() if has_arrow else 0.0)
    if width <= room:
        return markup
    size = max(5.0, PAGE_ONE_SIZE * room / width)
    return f"<font size='{size:.1f}'>{markup}</font>"


def _scale_rows(
    rows: Sequence[tuple[str, float | None, str, Status, str | None]],
    st: dict[str, ParagraphStyle],
    width: float = CONTENT_WIDTH,
    *,
    chips: bool = True,
) -> Table:
    """Name · scale bar · percentile · chip — the row every block shares.

    Each row's last element is the destination it links to, or None; a
    linked row is clickable across its whole width.

    At half width the chip has nowhere to go, so `chips=False` drops it and
    redistributes its share to the name, which is the column that actually
    needs the room when two of these sit side by side.
    """
    # Without the chip there are three columns, and the percentile needs
    # more of the freed space than its full-width share: ">99th" has to fit
    # on one line, and at half width the old proportion broke it across two.
    cols = _SCALE_COLS if chips else _SCALE_COLS_NARROW
    widths = [width * w for w in cols]
    grid: list[list[Any]] = []
    name_avail = widths[0] - 2 * _CELL_PAD
    for name, percentile, higher_means, status, href in rows:
        label = Paragraph(_one_line(name, name_avail), st["cell"])
        grid.append(
            [
                linked_cell(RowLink(href=href), label) if href else label,
                PercentileBar(
                    width=widths[1] - 6 * mm,
                    percentile=percentile,
                    higher_means=higher_means,
                    marker_colour=status.colour,
                    height=5.6 * mm,
                ),
                Paragraph(
                    f"<font color='{_hex(status.colour)}'>{_ordinal(percentile)}</font>",
                    st["pct"],
                ),
                *(
                    [StatusChip(status, width=widths[3] - 4 * mm, height=4.8 * mm,
                                font_size=6.0)]
                    if chips else []
                ),
            ]
        )
    return _table(
        grid,
        widths,
        extra=[
            ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
            ("TOPPADDING", (0, 0), (-1, -1), 1.6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
            *([("LEFTPADDING", (3, 0), (3, -1), 4)] if chips else []),
        ],
    )


def _metabolite_block(
    story: list[Any], st: dict[str, ParagraphStyle], *, rows: Sequence[MetaboliteRow],
    limit: int = SUMMARY_MAX_FUNCTIONS,
    width: float = CONTENT_WIDTH,
    chips: bool = True,
) -> None:
    # Spec 0.8.3 §12.1: the umbrella is wider than production. Synthesis,
    # degradation, substrate use, transformation and signalling all sit
    # under it. The section and metric IDs are unchanged; only the label is.
    _block_heading(
        story, st, "Microbial Functions",
        f"sections {SECTIONS['functions']}, {SECTIONS['functions_detail']}",
        width=width, href=section_dest("functions"),
    )
    # `rows` arrive most notable first, so the first `limit` are the ones
    # that matter most; once the notable readings are shown, the typical
    # ones fill whatever room the column has left rather than leaving it
    # blank. The caller sets `limit` from the height available.
    shown = list(rows[:limit])
    rest = len(rows) - len(shown)
    rest_notable = sum(
        1 for r in rows[len(shown):] if r.status.level != 0 and r.percentile is not None
    )
    data = [
        (
            f"<b>{short_metabolite(row.metabolite)}</b>"
            + {
                "adverse": f" <font color='{_hex(CORAL)}' size='6'>\u25b2</font>",
                "favourable": f" <font color='{_hex(GREEN)}' size='6'>\u25bc</font>",
            }.get(row.higher_means, ""),
            row.percentile,
            row.higher_means,
            row.status,
            function_dest(row.panel, detail=True),
        )
        for row in shown
    ]
    story.append(_scale_rows(data, st, width, chips=chips))
    # One compact legend line. The full explanation of arrows, percentiles and
    # the reference band lives beside the data in the detail sections; page 1
    # is a dashboard and every line of small print here costs a row of data.
    if rest > 0:
        more = f"<b>{rest} more</b>"
        more += f" ({rest_notable} notable)" if rest_notable else ", all typical"
        story.append(Spacer(1, 0.8 * mm))
        story.append(
            Paragraph(
                f"<font color='{_hex(CORAL)}'>\u25b2</font> worse "
                f"&nbsp; <font color='{_hex(GREEN)}'>\u25bc</font> better "
                f"&nbsp;·&nbsp; {more} in section {SECTIONS['functions']}",
                st["fine"],
            )
        )


def _community_block(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    similarity: Any | None,
    rpob_profile: dict[str, Any] | None,
    findings: Any | None = None,
    width: float = CONTENT_WIDTH,
    inv: Any = None,
    primary_phyla: dict[str, float] | None = None,
) -> None:
    _block_heading(
        story, st, "Microbiome Diversity",
        f"sections {SECTIONS['community']} to {SECTIONS['organisms']}, "
        f"{SECTIONS['groups_detail']}, {SECTIONS['organisms_detail']}",
        width=width,
    )
    taxonomy = similarity.taxonomy if similarity is not None else None
    if taxonomy is None:
        story.append(
            Paragraph("Species-level composition was not measured for this sample.", st["small"])
        )
        return

    from openbiota.scoring import ecological_metrics
    from openbiota.similarity import rescale_to_classified

    species = {k: v for k, v in rescale_to_classified(taxonomy.species).items() if v > 0}
    metrics = ecological_metrics(species)
    # The phylum ring is the primary composition, the same one the community page and
    # the organism section draw. Falls back to the scoring lane when a run
    # predates the primary catalogue.
    if primary_phyla:
        tot = sum(primary_phyla.values()) or 1.0
        phyla = sorted(((k, 100.0 * v / tot) for k, v in primary_phyla.items()), key=lambda kv: -kv[1])
    else:
        phyla = sorted(rescale_to_classified(taxonomy.phyla).items(), key=lambda kv: -kv[1])
    # Anything that would round to 0% is folded into "Other" rather than listed.
    shown = [(n, v) for n, v in phyla[:6] if v >= 0.5]
    rest = sum(v for _, v in phyla) - sum(v for _, v in shown)
    parts = [(name, value, PHYLA_PALETTE[i % len(PHYLA_PALETTE)]) for i, (name, value) in enumerate(shown)]
    if rest >= 0.5:
        parts.append(("Other", rest, PHYLA_PALETTE[-1]))

    # The headline count is every organism found, pooled across catalogues,
    # not the subset one catalogue could name. The Shannon figure beside it
    # stays on the scored lane because it carries a percentile and a
    # percentile needs a reference population measured the same way.
    n_detected = len(inv) if inv is not None and len(inv) else len(species)
    donut = Donut(parts=parts, size=33 * mm, centre_value=f"{n_detected:,}",
                  centre_label="organisms")
    counts = inv.counts() if inv is not None and len(inv) else {}
    n_all = counts.get("organisms", 0)
    n_genera = counts.get("genera", 0)
    n_strain = counts.get("strain_resolved", 0)

    legend_rows = [
        [
            Dot(colour, size=2.0 * mm),
            Paragraph(
                f"{name.replace('_', ' ')} <font color='{_hex(INK_FAINT)}'>{value:.0f}%</font>",
                st["cell_soft"],
            ),
        ]
        for name, value, colour in parts
    ]
    # Narrow enough and the three-across arrangement stops fitting, so the
    # block folds: donut beside its legend, tiles in a grid underneath.
    narrow = width < 150 * mm
    # The column has to clear the donut plus a real gutter, or the legend
    # sits on the ring.
    donut_col = donut.size + 6 * mm if narrow else donut.size + 10 * mm
    legend_col = width - donut_col if narrow else 50 * mm
    legend = Table(
        legend_rows, colWidths=[4 * mm, max(20 * mm, legend_col - 4 * mm)],
        hAlign="LEFT",
    )
    legend.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 1.1),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.1),
            ]
        )
    )

    # Where the community sits: diversity and evenness against the reference,
    # then the organism census — how many expected species are missing, how
    # many are depleted or expanded, how many need qualification. The census
    # is the answer to "what should be there and is not".
    def _n(bucket: int) -> int:
        return len(findings.bucket(bucket)) if findings is not None else 0

    n_missing, n_low, n_high, n_concern = _n(4), _n(3), _n(2), _n(6)
    # A census that was never taken has no counts. Without this, every tile
    # below printed a green zero - "0 expected species missing", "0 / 0",
    # "none need qualification" - which reads as a clean result when the
    # organism engine did not run at all.
    census_made = findings is not None
    from openbiota.pdfatlas import community_metric_percentiles

    pcts = community_metric_percentiles(similarity)
    div_pct, even_pct = pcts.get("shannon_diversity"), pcts.get("evenness")
    tile_w = (
        (width - 3 * mm) / 2 if narrow
        else (width - 40 * mm - 50 * mm - 3 * mm) / 2
    )
    tiles = Table(
        [
            [
                Tile(value=f"{metrics['shannon_diversity']:.2f}", label="Diversity",
                     note="Shannon index" + (f" · {_ordinal(div_pct)} pct" if div_pct is not None else ""),
                     width=tile_w, height=12 * mm),
                Tile(value=f"{metrics['evenness']:.2f}", label="Evenness",
                     note="0 to 1" + (f" · {_ordinal(even_pct)} pct" if even_pct is not None else ""),
                     width=tile_w, height=12 * mm),
            ],
            [
                Tile(value=str(n_missing) if census_made else "\u2014",
                     label="Expected species missing",
                     note=("in most reference adults, absent here" if census_made
                           else "the organism census did not run"),
                     width=tile_w, height=12 * mm,
                     accent=(CORAL if n_missing else GREEN) if census_made else SLATE),
                Tile(value=f"{n_low} / {n_high}" if census_made else "\u2014",
                     label="Depleted / expanded",
                     note=(
                         (f"{n_concern} need qualification" if n_concern
                          else "none need qualification")
                         if census_made else "not assessed"
                     ),
                     width=tile_w, height=12 * mm,
                     accent=(AMBER if (n_low or n_high) else GREEN) if census_made else SLATE),
            ],
            [
                Tile(value=f"{n_all:,}" if n_all else "\u2014",
                     label="All organisms detected",
                     note=(f"{n_genera} genera, every catalogue pooled"
                           if n_all else "no organism list"),
                     width=tile_w, height=12 * mm,
                     accent=ACCENT if n_all else SLATE),
                Tile(value=f"{n_strain:,}" if n_strain else "\u2014",
                     label="Resolved to strain",
                     note=("read at single-base resolution" if n_strain
                           else "strain typing did not run"),
                     width=tile_w, height=12 * mm,
                     accent=ACCENT if n_strain else SLATE),
            ],
        ],
        colWidths=[tile_w + 3 * mm, tile_w],
        hAlign="LEFT",
    )
    tiles.setStyle(
        TableStyle(
            [
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 1.5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5),
            ]
        )
    )

    if narrow:
        block = Table(
            [[donut, legend], [tiles, ""]],
            colWidths=[donut_col, legend_col],
            hAlign="LEFT",
        )
        block.setStyle(TableStyle([("SPAN", (0, 1), (1, 1))]))
    else:
        block = Table(
            [[donut, legend, tiles]],
            colWidths=[donut_col, legend_col, width - donut_col - legend_col],
            hAlign="LEFT",
        )
    block.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    story.append(block)
    _ = rpob_profile


#: Height of the age tile and chart at half width: enough for the number,
#: its label and the interval with even margins, and no more. The chart
#: beside it takes the same height so the row is one clean band.
AGE_ROW_H: Final = 22 * mm
AGE_TILE_W: Final = 30 * mm


def _age_block(
    story: list[Any], st: dict[str, ParagraphStyle], *, age: Any | None, age_meta: dict[str, Any] | None,
    width: float = CONTENT_WIDTH,
) -> None:
    """Estimated microbiome age as one infographic: number beside chart.

    No small print: the number and its interval are the reading, the chart
    puts them against the reference adults, and the arithmetic behind the
    figure lives in the age section.
    """
    from openbiota.pdfatlas import AgeScale, _json_field

    _block_heading(story, st, "Biological Age Estimate", f"section {SECTIONS['age']}",
                   width=width)
    if age is None:
        story.append(Paragraph("The age model was not available for this run.", st["small"]))
        return
    meta = dict(age_meta or {})
    support = _json_field(meta, "training_age_support") or {}
    model_point = age.predicted_chronological_age_years
    adj = getattr(age, "health_adjustment", None)
    # What the page shows is the health-adjusted figure where it exists: the
    # composition model alone reads every sample as its adult training mean.
    pred = getattr(age, "reported_age_years", None) or model_point
    pi80, pi95 = age.prediction_interval_80, age.prediction_interval_95
    pi90 = getattr(age, "prediction_interval_90", None)
    scored = getattr(age, "status", None) == "scored" and pred is not None
    status_colour = {
        "supported": GREEN, "warning—in distribution": AMBER,
        "abstained—out of distribution": CORAL, "not computable": SLATE,
    }.get(age.status_line, SLATE)
    # The "no age signal" caveat is about the *model's* point estimate. Once
    # the health adjustment has moved the answer, the figure on the page is
    # no longer the reference average, so the flat caveat would be wrong.
    uninformative = str(getattr(age, "information_verdict", "") or "").startswith(
        "indistinguishable_from_population_average") and adj is None

    if adj is not None and pi90:
        note = f"90% interval {pi90[0]:.0f}–{pi90[1]:.0f}"
    elif uninformative:
        note = "= the reference average; no age signal found"
    elif pi80:
        note = f"80% interval {pi80[0]:.0f}–{pi80[1]:.0f}"
    else:
        note = ""
    row_h = AGE_ROW_H
    gutter = 4 * mm
    tile_w = AGE_TILE_W if width < 150 * mm else 40 * mm
    number = Tile(
        value="—" if not scored else f"{pred:.0f}",
        label="Estimated age (yrs)",
        note=age.status_line if not scored else note,
        width=tile_w, height=row_h,
        accent=SLATE if uninformative else status_colour, big=True, centred=True,
    )
    # With the health adjustment the answer is held inside the 90% band, so
    # that is the band drawn; without it, the model's own 80% and 95%.
    scale = AgeScale(
        width=width - tile_w - gutter, height=row_h,
        predicted=pred if scored else None,
        pi80=None if adj is not None else (pi80 if scored else None),
        pi95=None if adj is not None else (pi95 if scored else None),
        pi90=pi90 if (adj is not None and scored) else None,
        actual=None,
        support_edges=support.get("histogram_edges", []),
        support_counts=support.get("histogram_counts", []),
        lo=5.0, hi=85.0, compact=True,
    )
    block = Table([[number, scale]], colWidths=[tile_w + gutter, width - tile_w - gutter], hAlign="LEFT")
    block.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0), ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    story.append(block)


def _profile_block(
    story: list[Any], st: dict[str, ParagraphStyle], *, similarity: Any | None,
    limit: int = SUMMARY_MAX_PROFILES,
    width: float = CONTENT_WIDTH,
    chips: bool = True,
) -> None:
    _block_heading(
        story, st, "Disease Signatures",
        f"sections {SECTIONS['patterns']}, {SECTIONS['patterns_detail']}",
        width=width, href=section_dest("patterns"),
    )
    if similarity is None or not similarity.results:
        story.append(Paragraph("No disease patterns were scored for this sample.", st["small"]))
        return
    ranked = [r for r in similarity.ranked if r.reportable]
    shown = ranked[:limit]
    data = [
        (
            f"<b>{short_profile_label(result.profile.name, result.profile.label)}</b>",
            result.combined_percentile,
            "adverse",
            profile_status(result),
            profile_dest(result.profile.name, detail=True),
        )
        for result in shown
    ]
    if data:
        story.append(_scale_rows(data, st, width, chips=chips))
    n_total = len(similarity.results)
    n_scored = len(ranked)
    story.append(Spacer(1, 0.8 * mm))
    story.append(
        Paragraph(
            f"Top {len(shown)} of {n_scored} scored ({n_total} in the library). "
            f"Resemblance is not disease &mdash; section {SECTIONS['patterns']}",
            st["fine"],
        )
    )


# --------------------------------------------------------------------------- #
# page 2 — what stood out: one spark graphic beside each finding
# --------------------------------------------------------------------------- #

SPARK_W = 46 * mm


class SparkGauge(Flowable):
    """A horizontal gauge for a signed index: red through amber to green, a marker, the value."""

    def __init__(self, *, value: float | None, lo: float = -6.0, hi: float = 6.0, width: float = SPARK_W) -> None:
        super().__init__()
        self.value, self.lo, self.hi = value, lo, hi
        self.width, self.height = width, 12 * mm

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        y, h = 4.2 * mm, 3.2 * mm
        n = 40
        seg = self.width / n
        for i in range(n):
            f = i / (n - 1)
            col = _lerp3(CORAL, AMBER, GREEN, f)
            c.setFillColor(col)
            c.rect(i * seg, y, seg + 0.3, h, stroke=0, fill=1)
        c.setStrokeColor(colors.white)
        c.setLineWidth(0.8)
        mid = self.width * (0 - self.lo) / (self.hi - self.lo)
        c.line(mid, y - 0.6 * mm, mid, y + h + 0.6 * mm)
        c.setFont("Helvetica", 5.4)
        c.setFillColor(INK_FAINT)
        c.drawString(0, 0.6 * mm, "disturbed")
        c.drawRightString(self.width, 0.6 * mm, "healthy")
        if self.value is not None:
            f = (min(max(self.value, self.lo), self.hi) - self.lo) / (self.hi - self.lo)
            x = f * self.width
            c.setFillColor(INK)
            c.setStrokeColor(colors.white)
            c.setLineWidth(1.0)
            c.circle(x, y + h / 2, 2.4, stroke=1, fill=1)
            c.setFont("Helvetica-Bold", 7.6)
            c.drawCentredString(min(max(x, 6 * mm), self.width - 6 * mm), y + h + 1.4 * mm, f"{self.value:+.2f}")


class SparkCensus(Flowable):
    """Stacked count bar: each part a colour and a count, labelled beneath."""

    def __init__(self, parts: Sequence[tuple[str, int, colors.Color]], *, width: float = SPARK_W) -> None:
        super().__init__()
        self.parts = [p for p in parts if p[1] > 0]
        self.total = sum(p[1] for p in self.parts)
        self.width, self.height = width, 12 * mm

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        y, h = 5.4 * mm, 3.6 * mm
        if not self.total:
            c.setFillColor(SLATE_BG)
            c.roundRect(0, y, self.width, h, 1.2, stroke=0, fill=1)
            c.setFont("Helvetica", 5.6)
            c.setFillColor(INK_FAINT)
            c.drawString(0, 0.8 * mm, "none")
            return
        x = 0.0
        for _label, n, col in self.parts:
            w = self.width * n / self.total
            c.setFillColor(col)
            c.rect(x, y, max(w - 0.6, 0.4), h, stroke=0, fill=1)
            if w > 5 * mm:
                c.setFillColor(colors.white)
                c.setFont("Helvetica-Bold", 6.4)
                c.drawCentredString(x + w / 2, y + 1.0 * mm, str(n))
            x += w
        # legend
        c.setFont("Helvetica", 5.4)
        x = 0.0
        for label, n, col in self.parts:
            c.setFillColor(col)
            c.circle(x + 1.1 * mm, 1.7 * mm, 1.1 * mm, stroke=0, fill=1)
            c.setFillColor(INK_SOFT)
            text = f"{n} {label}"
            c.drawString(x + 2.8 * mm, 0.9 * mm, text)
            x += c.stringWidth(text, "Helvetica", 5.4) + 5.2 * mm
            if x > self.width - 8 * mm:
                break


class SparkStrips(Flowable):
    """Up to three named mini percentile strips, marker coloured by status."""

    def __init__(self, rows: Sequence[tuple[str, float | None, str, colors.Color]], *, width: float = SPARK_W) -> None:
        super().__init__()
        self.rows = list(rows)[:3]
        self.width = width
        self.height = max(1, len(self.rows)) * 4.6 * mm + 1.0 * mm

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        label_w = 20 * mm
        bar_w = self.width - label_w
        for i, (name, pct, higher_means, col) in enumerate(self.rows):
            y = self.height - (i + 1) * 4.6 * mm + 0.6 * mm
            c.setFont("Helvetica", 5.3)
            c.setFillColor(INK_SOFT)
            c.drawRightString(label_w - 1.6 * mm, y + 0.6 * mm, _fit_text(_short(name), "Helvetica", 5.3, label_w - 2 * mm))
            bar = PercentileBar(width=bar_w, percentile=pct, higher_means=higher_means, marker_colour=col,
                                height=3.2 * mm, show_scale=False)
            bar.canv = c
            c.saveState()
            c.translate(label_w, y - 0.4 * mm)
            bar.draw()
            c.restoreState()


class SparkDots(Flowable):
    """Every scored pattern as a dot on one 0–100 axis; the typical band shaded."""

    def __init__(self, percentiles: Sequence[tuple[float, colors.Color]], *, width: float = SPARK_W) -> None:
        super().__init__()
        self.pts = list(percentiles)
        self.width, self.height = width, 12 * mm

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        y = 5.6 * mm
        c.setFillColor(SLATE_BG)
        c.rect(0, y - 1.6 * mm, self.width, 3.2 * mm, stroke=0, fill=1)
        c.setFillColor(colors.Color(0.85, 0.90, 0.88))
        c.rect(self.width * 0.25, y - 1.6 * mm, self.width * 0.5, 3.2 * mm, stroke=0, fill=1)
        # jitter identical percentiles vertically so a cluster reads as a cluster
        seen: dict[int, int] = {}
        for pct, col in sorted(self.pts):
            k = int(pct // 3)
            j = seen.get(k, 0)
            seen[k] = j + 1
            x = self.width * min(max(pct, 0), 100) / 100
            c.setFillColor(col)
            c.setStrokeColor(colors.white)
            c.setLineWidth(0.5)
            c.circle(x, y + ((j % 3) - 1) * 1.3 * mm, 1.5, stroke=1, fill=1)
        c.setFont("Helvetica", 5.4)
        c.setFillColor(INK_FAINT)
        for v in (0, 25, 50, 75, 100):
            c.drawCentredString(self.width * v / 100, 0.6 * mm, str(v))


def _short(name: str) -> str:
    """Drop parenthetical qualifiers from a metabolite label for a spark axis."""
    return name.split(" (")[0].replace("_", " ")


def _lerp3(a: colors.Color, b: colors.Color, cc: colors.Color, f: float) -> colors.Color:
    if f < 0.5:
        t = f / 0.5
        return colors.Color(a.red + (b.red - a.red) * t, a.green + (b.green - a.green) * t, a.blue + (b.blue - a.blue) * t)
    t = (f - 0.5) / 0.5
    return colors.Color(b.red + (cc.red - b.red) * t, b.green + (cc.green - b.green) * t, b.blue + (cc.blue - b.blue) * t)


def insights_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    rows: Sequence[MetaboliteRow],
    similarity: Any | None,
    meta: dict[str, Any],
    findings: Any | None = None,
    age: Any | None = None,
    plan: Any | None = None,
    age_meta: dict[str, Any] | None = None,
    inv: Any = None,
) -> None:
    """Page 2: each finding as a spark graphic beside one or two plain sentences.

    Ordered as the report is ordered — the community, who is missing, the
    age, the functions, the patterns, then what can be done and what the two
    engines disagree about.
    """
    from openbiota.pdfatlas import AgeScale, _json_field

    story.extend(section_heading(SECTIONS["stood_out"], "What stood out", st["h1"]))
    story.append(Paragraph(
        "The findings that matter most, each with the graphic it comes from. Everything here is expanded in "
        "the sections named at the end of each line.",
        st["body"],
    ))
    story.append(Spacer(1, 1.5 * mm))

    # Each row carries a priority. The page is filled by priority and
    # measured; whatever does not fit is left out rather than squeezed.
    items: list[tuple[int, Any, str]] = []

    # ---- 1 · general gut health --------------------------------------------- #
    anchor = similarity.anchor if similarity is not None else None
    if anchor is not None and anchor.valid:
        words, _c = anchor_words(anchor)
        if anchor.strongly_negative:
            text = (
                f"<b>General gut health: {anchor.score:+.2f}, {words.lower()}.</b> The most important number here. "
                "A community in this state resembles several disease patterns at once, so read those as reflecting "
                "the broad disturbance rather than any one condition."
            )
        elif anchor.score <= -0.25:
            text = (f"<b>General gut health: {anchor.score:+.2f}, {words.lower()}.</b> Some of any raised pattern "
                    "score may be general disturbance rather than a specific signal.")
        elif anchor.score < 0.25:
            text = f"<b>General gut health: {anchor.score:+.2f}, {words.lower()}.</b> Neither clearly healthy nor clearly disturbed."
        else:
            text = (f"<b>General gut health: {anchor.score:+.2f}, {words.lower()}</b> on an index trained on 8,069 "
                    "labelled stool samples. Raised pattern scores below are therefore not just tracking disturbance.")
        items.append((100, SparkGauge(value=anchor.score), text + f" <font color='#8C99A6'>Section {SECTIONS['community']}.</font>"))

    # ---- 2 · who is missing, who is low or high ----------------------------- #
    if findings is not None:
        missing, low, high = list(findings.bucket(4)), list(findings.bucket(3)), list(findings.bucket(2))
        unusual, concern = list(findings.bucket(5)), list(findings.bucket(6))
        n_typ = len(findings.bucket(1))
        census = SparkCensus([("missing", len(missing), CORAL), ("low", len(low), AMBER),
                              ("high", len(high), ORANGE), ("typical", n_typ, SLATE)])
        if missing or low or high:
            bits = []
            if missing:
                bits.append(f"<b>{len(missing)} expected species missing</b> ("
                            + ", ".join(f"<i>{t.display_name}</i>" for t in missing[:3])
                            + (f", +{len(missing) - 3}" if len(missing) > 3 else "") + ")")
            if low:
                bits.append(f"<b>{len(low)} depleted</b> (" + ", ".join(f"<i>{t.display_name}</i>" for t in low[:2]) + ")")
            if high:
                bits.append(f"<b>{len(high)} expanded</b> (" + ", ".join(f"<i>{t.display_name}</i>" for t in high[:2]) + ")")
            text = "; ".join(bits) + f"; {n_typ} typical."
        else:
            text = (
                f"Nothing commonly carried is missing and all {n_typ} placed "
                "species sit in the reference's central range."
            )
        # Species the reference could not place at all. They are neither
        # missing nor typical - their absence simply cannot be judged here -
        # and leaving them out made a partial census read as a complete one.
        unplaceable = list(findings.bucket(7))
        if unplaceable:
            text += (
                f" A further <b>{len(unplaceable)}</b> expected "
                f"{'species was' if len(unplaceable) == 1 else 'species were'} "
                "not placed against the reference at all ("
                + ", ".join(f"<i>{t.display_name}</i>" for t in unplaceable[:3])
                + (f", +{len(unplaceable) - 3}" if len(unplaceable) > 3 else "")
                + "), so nothing here says whether they are present or absent."
            )
        items.append((90, census, text + f" <font color='#8C99A6'>Sections {SECTIONS['organisms']} and {SECTIONS['organisms_detail']}.</font>"))

    # ---- overgrowths that matter, by consequence ------------------------ #
    # The census row above counts expansions and names two by bucket order,
    # which is how a 13% bloom of a Crohn's-associated pathobiont went
    # unnamed while two probiotics were named. This row is the same
    # question answered by weight: percentile excess times abundance, an
    # opportunist ahead of a conditional organism.
    if inv is not None and len(inv):
        from openbiota import organisms as org

        vs = org.verdicts(list(inv.organisms))
        over = [v for v in org.issues(vs)
                if v.flag == "high" and v.cls in (org.OPPORTUNIST, org.CONDITIONAL)]
        if over:
            lead = over[0]
            spark = SparkStrips([
                (v.organism.display.split(" ")[0][:1] + ". " + " ".join(v.organism.display.split(" ")[1:]),
                 v.organism.level_percentile, "adverse",
                 CORAL if v.is_issue else AMBER)
                for v in over[:3]
            ])
            what = lead.description.split(".")[0].rstrip(".")
            dev = _deviation(lead.organism)
            text = (
                f"<b>Overgrown: <i>{lead.organism.display}</i></b> at "
                f"{_pct_short(lead.organism.percent) if lead.organism.in_primary else 'a detected level'} of your community, "
                + (f"{dev} against the typical carrier, " if dev else "")
                + f"{lead.flag_reason.split(' — ')[-1]} \u2014 {lead.label.lower()}; "
                f"{what[:1].lower() + what[1:]}."
            )
            if len(over) > 1:
                others = ", ".join(f"<i>{v.organism.display}</i>" for v in over[1:3])
                more = len(over) - 3
                text += (f" Also above range: {others}"
                         + (f" and {more} more" if more > 0 else "") + ".")
            text += (f" <font color='#8C99A6'>Section {SECTIONS['organisms']} lists every overgrown "
                     f"organism; section {SECTIONS['organisms_detail']} has a page on each.</font>")
            items.append((95, spark, text))

    # The census's own qualifications (unmatched comparator, detection
    # prior, pediatric gate) are deliberately NOT rows here. Every row on
    # this page is a finding with the graphic it came from; a caveat has
    # no graphic, and rendered as a row it sits beside an empty cell and
    # reads as a defect. They live in the fine print of the organisms
    # section, which is where a caveat belongs.
    if findings is not None and (unusual or concern):
            good = [t for t in unusual if _polarity(t) == "beneficial"]
            bad = [t for t in unusual if _polarity(t) == "adverse"]
            neutral = [t for t in unusual if t not in good and t not in bad]
            census2 = SparkCensus([("beneficial", len(good), GREEN), ("neutral", len(neutral), SLATE),
                                   ("adverse", len(bad), CORAL), ("concern", len(concern), CORAL)])
            bits = []
            if unusual:
                bits.append(f"<b>{len(unusual)} uncommon in the reference</b>"
                            + (f", {len(good)} of them health-associated" if good else "")
                            + (f", {len(bad)} adverse" if bad else ""))
            if concern:
                bits.append(f"<b>{len(concern)} need qualification</b> ("
                            + ", ".join(f"<i>{t.display_name}</i>" for t in concern[:2]) + ") — species DNA is not a strain, toxin or infection")
            items.append((60, census2, "; ".join(bits) + f". <font color='#8C99A6'>Sections {SECTIONS['organisms']} and {SECTIONS['organisms_detail']}.</font>"))

    # ---- 3 · age ------------------------------------------------------------ #
    if age is not None:
        support = _json_field(dict(age_meta or {}), "training_age_support") or {}
        scored = getattr(age, "status", None) == "scored" and age.predicted_chronological_age_years is not None
        # Plot the figure the report states, not the unadjusted model point.
        reported = getattr(age, "reported_age_years", None)
        spark = AgeScale(width=SPARK_W, predicted=reported if scored else None,
                         pi80=age.prediction_interval_80 if scored else None, pi95=age.prediction_interval_95 if scored else None,
                         actual=None, support_edges=support.get("histogram_edges", []),
                         support_counts=support.get("histogram_counts", []), compact=True)
        text = age_headline(age)
        items.append((70, spark, text + f" <font color='#8C99A6'>Section {SECTIONS['age']}.</font>"))

    # ---- 4 · functions ------------------------------------------------------ #
    judged = [r for r in rows if r.status.colour in (CORAL, ORANGE, AMBER) and r.percentile is not None]
    good_rows = [r for r in rows if r.status.colour is GREEN and abs(r.status.level) >= 1]
    if judged:
        spark = SparkStrips([(r.metabolite, r.percentile, r.higher_means, r.status.colour) for r in judged[:3]])
        text = (f"<b>{len(judged)} of {len(rows)} metabolic functions</b> sit outside the typical band in the less "
                "favourable direction: " + ", ".join(f"<b>{r.metabolite}</b> ({_ordinal(r.percentile)})" for r in judged[:3])
                + (f", +{len(judged) - 3} more" if len(judged) > 3 else "") + ".")
    else:
        spark = SparkStrips([(r.metabolite, r.percentile, r.higher_means, r.status.colour) for r in rows[:3]])
        # Same correction as the page-1 line: only count what was placed.
        assessed_rows = [r for r in rows if r.status.assessed and r.percentile is not None]
        held_back = len(rows) - len(assessed_rows)
        if assessed_rows:
            text = (
                f"All {len(assessed_rows)} of {len(rows)} metabolic functions that "
                "could be placed against the reference are typical of it, or "
                "outside it only favourably."
            )
            if held_back:
                reasons = sorted({
                    r.status.unassessed_reason for r in rows
                    if not (r.status.assessed and r.percentile is not None)
                })
                text += (
                    f" The remaining {held_back} were not placed ("
                    + ", ".join(reasons)
                    + ") and are neither typical nor atypical."
                )
        else:
            text = (
                f"None of the {len(rows)} metabolic functions could be placed "
                "against the reference, so none is reported as typical."
            )
    items.append((65, spark, text + f" <font color='#8C99A6'>Sections {SECTIONS['functions']} and {SECTIONS['functions_detail']}.</font>"))
    if good_rows:
        spark = SparkStrips([(r.metabolite, r.percentile, r.higher_means, r.status.colour) for r in good_rows[:3]])
        items.append((40, spark, "<b>Favourably placed:</b> "
                      + ", ".join(f"<b>{r.metabolite}</b> ({_ordinal(r.percentile)})" for r in good_rows[:3])
                      + f". <font color='#8C99A6'>Section {SECTIONS['functions']}.</font>"))

    # ---- 5 · patterns ------------------------------------------------------- #
    if similarity is not None and similarity.results:
        scored = [r for r in similarity.results if not r.abstention.abstained and r.combined_percentile is not None]
        high = [r for r in scored if profile_status(r).level >= 1]
        shape = getattr(similarity, "shape", None)
        dots = SparkDots([(r.combined_percentile, profile_status(r).colour) for r in scored])
        if high and shape is not None and shape.verdict == "diffuse":
            text = (f"<b>{shape.n_above_typical} of {shape.n_scored} disease patterns above typical</b>, {shape.n_notably_high} "
                    "notably high. That many rising together reads as <b>general disturbance</b>, not a specific condition: "
                    "the patterns share features and one disturbed community lights them all.")
        elif high:
            top = max(high, key=lambda r: r.combined_percentile)
            spec = " It stands alone." if shape is not None and shape.verdict == "specific" else ""
            text = (f"Closest disease pattern: <b>{top.profile.label}</b>, {_ordinal(top.combined_percentile)} percentile "
                    f"({profile_status(top).label}).{spec}")
        elif scored:
            text = f"None of the {len(scored)} disease patterns showed above-average resemblance."
        else:
            text = "No disease pattern could be scored."
        abst = [r for r in similarity.results if r.abstention.abstained]
        if abst:
            text += f" {len(abst)} pattern{'s' if len(abst) != 1 else ''} not scored because a required condition was not met."
        items.append((75, dots, text + f" <font color='#8C99A6'>Sections {SECTIONS['patterns']} and {SECTIONS['patterns_detail']}.</font>"))

    # ---- 6 · what can be done ----------------------------------------------- #
    if plan is not None and plan.summaries:
        n_cards = len(plan.summaries)
        census = SparkCensus([("human trials", plan.n_with_human_trials, GREEN),
                              ("clinical route", plan.n_clinical_routes, ACCENT),
                              ("none studied", plan.n_no_supported, SLATE),
                              ("other", max(0, n_cards - plan.n_with_human_trials - plan.n_clinical_routes - plan.n_no_supported), INK_FAINT)])
        items.append((30, census, f"<b>{n_cards} evidence cards</b> for the readings that fired: {plan.n_with_human_trials} with human "
                      f"trial or guideline evidence for an exact product, diet or drug; {plan.n_clinical_routes} pointing to a "
                      f"conventional clinical route; {plan.n_no_supported} stating that nothing studied exists for that reading. "
                      f"<font color='#8C99A6'>Section {SECTIONS['actions']}, and under every reading in Part B.</font>"))

    # ---- 7 · where two measurements disagree ------------------------------- #
    if plan is not None and plan.coherence.divergent:
        c0 = plan.coherence.divergent[0]
        spark = SparkStrips([
            ("gene search", c0.panel_percentile, "favourable", AMBER),
            ("named species", c0.group_percentile, "favourable", AMBER),
        ])
        more = f" (+{len(plan.coherence.divergent) - 1} more)" if len(plan.coherence.divergent) > 1 else ""
        items.append((20, spark, f"<b>Two measurements of {c0.pair.capacity_label} disagree</b>{more}: the gene search puts it at "
                      f"the {_ordinal(c0.panel_percentile)} percentile, the named species at the {_ordinal(c0.group_percentile)}. "
                      "Not an error in either — they count different things, and each card explains the split and which to weight. "
                      f"<font color='#8C99A6'>Sections {SECTIONS['groups']}, {SECTIONS['functions']}, {SECTIONS['groups_detail']} and {SECTIONS['functions_detail']}.</font>"))

    # Invariant: every row is a finding with the graphic it came from. A row
    # with no graphic is a caveat that wandered into a findings list, and
    # it renders as text beside an empty cell. Caveats belong in the fine
    # print of the section they qualify. This raises at build time rather
    # than letting a blank cell reach the page.
    blank = [text[:60] for _, spark, text in items if spark is None or spark == ""]
    if blank:
        raise ValueError(
            "insights page: every row needs a graphic; these have none: " + "; ".join(blank)
        )

    # One page, always. Rows are ranked by how much they matter, then the
    # page is filled from the top and measured against the frame exactly as
    # the frame will measure it. A row that does not fit is left out; the
    # layout is never scaled to make room, because a shrunken page reads as
    # a mistake and a missing low-priority row does not.
    ranked = sorted(items, key=lambda it: -it[0])
    style = TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (0, -1), 5 * mm),
        ("RIGHTPADDING", (1, 0), (1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 3.0 * mm), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.0 * mm),
        ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
    ])

    def build(rows: list[tuple[int, Any, str]]) -> Table:
        # Back in report order once chosen, so the page still reads top to
        # bottom the way the report does.
        chosen = sorted(rows, key=lambda it: items.index(it))
        grid = [[spark, Paragraph(text, st["body_ink"])] for _, spark, text in chosen]
        t = Table(grid, colWidths=[SPARK_W + 5 * mm, CONTENT_WIDTH - SPARK_W - 5 * mm], hAlign="LEFT")
        t.setStyle(style)
        return t

    probe = Canvas(io.BytesIO(), pagesize=PAGE)
    used_above = _listWrapOn(list(story[-3:]), CONTENT_WIDTH, probe)[1]
    budget = BODY_FRAME_HEIGHT - used_above - 4 * mm
    keep = list(ranked)
    while keep:
        table = build(keep)
        if _listWrapOn([table], CONTENT_WIDTH, probe)[1] <= budget:
            break
        keep.pop()  # the least important goes first
    if keep:
        story.append(build(keep))
    _ = meta


def _polarity(t: Any) -> str:
    """Curated polarity of a finding's species: beneficial, adverse, or neutral/unknown."""
    from openbiota.pdfatlas import polarity_word

    word = polarity_word(t)
    return word if word in ("beneficial", "adverse") else "neutral"
