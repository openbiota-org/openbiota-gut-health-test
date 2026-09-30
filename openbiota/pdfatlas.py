"""Front-atlas pages and evidence cards for the PDF (spec v04 §10).

The report is in two parts. **Part A — at a glance** puts every reading in
front of the reader as a graphic before any prose: sample validity, the
global community metrics on neutral range strips, all functional panels, all
microbial groups, the organisms that need attention (what is missing, what is
low, what is high, what is unusual, what needs qualification), every
disease-pattern profile, the evidence overview and the estimated microbiome
age. **Part B — in detail** then breaks every one of those readings down, and
under each one shows *what the human research says about changing it*: the
evidence cards produced by :mod:`openbiota.actions`.

Everything here draws from the same seven-band scale and the same flowables
as the rest of the report, so the eye learns the language once.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from typing import Any, Final

import numpy as np
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    CondPageBreak,
    Flowable,
    KeepInFrame,
    KeepTogether,
    PageBreak,
    Spacer,
    Table,
    TableStyle,
)

from openbiota.findings import TaxonFinding
from openbiota.interventions import LANE_MEANING, EvidenceSummary
from openbiota.pdflinks import (
    LinkedTable,
    Paragraph,
    RowLink,
    group_dest,
    linked_block,
    linked_cell,
    section_dest,
    section_heading,
)
from openbiota.pdfreport import (
    ACCENT,
    ACCENT_DARK,
    AMBER,
    AMBER_BG,
    CONTENT_WIDTH,
    CORAL,
    CORAL_BG,
    GREEN,
    GREEN_BG,
    INK,
    INK_FAINT,
    ORANGE,
    ORANGE_BG,
    PANEL_BG,
    RULE,
    SLATE,
    SLATE_BG,
    PercentileBar,
    Rule,
    Status,
    StatusChip,
    Tile,
    _hex,
    _ordinal,
    status_for,
)
from openbiota.safety import participant_safety_sentence

# --------------------------------------------------------------------------- #
# shared bits
# --------------------------------------------------------------------------- #

LANE_COLOUR: dict[str, tuple[colors.Color, colors.Color]] = {
    "A": (GREEN, GREEN_BG), "B": (GREEN, GREEN_BG), "C": (AMBER, AMBER_BG),
    "D": (ORANGE, ORANGE_BG), "E": (SLATE, SLATE_BG), "X": (CORAL, CORAL_BG), "none": (SLATE, SLATE_BG),
}
LANE_SHORT: dict[str, str] = {
    "A": "guideline / label", "B": "replicated RCT benefit", "C": "single or conflicting RCT",
    "D": "non-randomised human", "E": "animal / in-vitro only", "X": "against-evidence", "none": "no supporting evidence",
}
RENDER_WORDS: dict[str, str] = {
    "guideline_route_after_confirmation": "clinical route",
    "clinician_review": "clinician review",
    "evidence_only": "evidence on record",
    "inform_only": "information",
    "not_actionable": "not actionable",
    "no_supported_targeted_intervention": "no supported intervention",
}
GATE_COLOUR = {"pass": GREEN, "warn": AMBER, "fail": CORAL, "unknown": SLATE, "not_assessable": SLATE}
GATE_BG = {"pass": GREEN_BG, "warn": AMBER_BG, "fail": CORAL_BG, "unknown": SLATE_BG, "not_assessable": SLATE_BG}
BUCKET_COLOUR = {1: GREEN, 2: ORANGE, 3: AMBER, 4: CORAL, 5: ORANGE, 6: CORAL, 7: SLATE}
BUCKET_BG = {1: GREEN_BG, 2: ORANGE_BG, 3: AMBER_BG, 4: CORAL_BG, 5: ORANGE_BG, 6: CORAL_BG, 7: SLATE_BG}
BUCKET_SHORT = {
    1: "typical", 2: "relatively high", 3: "relatively low", 4: "missing",
    5: "unusual to have", 6: "needs qualification", 7: "indeterminate",
}


def _p(text: str, st: ParagraphStyle) -> Paragraph:
    return Paragraph(text, st)


def _plain(rows: list[list[Any]], widths: Sequence[float], *, header: bool = True, zebra: bool = False,
           extra: Sequence[Any] = ()) -> Table:
    t = LinkedTable(rows, colWidths=list(widths), hAlign="LEFT", repeatRows=1 if header else 0)
    style: list[Any] = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2.4), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LINEBELOW", (0, 1 if header else 0), (-1, -1), 0.3, RULE),
    ]
    if header:
        style.append(("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE))
    if zebra:
        for i in range(1 if header else 0, len(rows)):
            if (i % 2) == 0:
                style.append(("BACKGROUND", (0, i), (-1, i), PANEL_BG))
    style.extend(extra)
    t.setStyle(TableStyle(style))
    return t


def _chip(text: str, colour: colors.Color, bg: colors.Color, width: float, *, size: float = 6.6) -> StatusChip:
    return StatusChip(Status(text, colour, bg, ""), width=width, height=4.8 * mm, font_size=size)


def section_divider(story: list[Any], st: dict[str, ParagraphStyle], part: str, title: str, blurb: str,
                    *, contents: Sequence[str] = ()) -> None:
    """A big, bold divider between the two parts of the report."""
    story.append(PageBreak())
    story.append(Spacer(1, 18 * mm))
    story.append(_p(part.upper(), st["kicker"]))
    story.append(_p(title, st["divider"]))
    story.append(Spacer(1, 2 * mm))
    story.append(Rule(CONTENT_WIDTH, colour=ACCENT, thickness=1.6))
    story.append(Spacer(1, 4 * mm))
    story.append(_p(blurb, st["lead"]))
    if contents:
        story.append(Spacer(1, 3 * mm))
        story.append(_p("IN THIS PART", st["label"]))
        for line in contents:
            story.append(_p(f"\u2022&nbsp;&nbsp;{line}", st["body_ink"]))


# --------------------------------------------------------------------------- #
# Part A · 3 — sample validity + global metrics
# --------------------------------------------------------------------------- #


def _cohort_metric_distributions(cohort: Any) -> dict[str, np.ndarray]:
    """Shannon, evenness and richness for every reference sample."""
    A = np.asarray(cohort.abundance, dtype=float)  # taxa × samples (percent)
    if A.size == 0:
        return {}
    with np.errstate(divide="ignore", invalid="ignore"):
        tot = A.sum(axis=0)
        P = np.where(tot > 0, A / np.where(tot > 0, tot, 1.0), 0.0)
        logP = np.where(P > 0, np.log(P), 0.0)
        H = -(P * logP).sum(axis=0)
        S = (A > 0).sum(axis=0)
        E = np.where(S > 1, H / np.log(np.maximum(S, 2)), 0.0)
    return {"shannon_diversity": H, "evenness": E, "richness": S.astype(float)}


def _pct_of(dist: np.ndarray | None, value: float) -> float | None:
    if dist is None or dist.size == 0:
        return None
    return float(100.0 * np.mean(dist < value) + 50.0 * np.mean(dist == value))


def community_metric_percentiles(similarity: Any | None) -> dict[str, float | None]:
    """Shannon diversity, evenness and richness as percentiles of the matched reference."""
    taxonomy = similarity.taxonomy if similarity is not None else None
    if taxonomy is None:
        return {}
    from openbiota.scoring import ecological_metrics
    from openbiota.similarity import rescale_to_classified

    species = {k: v for k, v in rescale_to_classified(taxonomy.species).items() if v > 0}
    metrics = ecological_metrics(species)
    dists = _cohort_metric_distributions(similarity.matched_cohort) if similarity.matched_cohort is not None else {}
    return {k: _pct_of(dists.get(k), metrics[k]) for k in ("shannon_diversity", "evenness", "richness")}


def validity_strip(story: list[Any], st: dict[str, ParagraphStyle], *, gates: Any | None, meta: Mapping[str, Any]) -> None:
    """The QC gates as a row of chips, with the one-line verdict."""
    story.append(_p("IS THIS SAMPLE FIT TO READ?", st["label"]))
    if gates is None:
        story.append(_p("Quality gates were not assembled for this run.", st["small"]))
        return
    overall = gates.overall
    cells: list[Any] = []
    for g in gates.gates:
        cells.append(Table(
            [[_chip(g.status.replace("_", " "), GATE_COLOUR.get(g.status, SLATE), GATE_BG.get(g.status, SLATE_BG), CONTENT_WIDTH / 7 - 4 * mm)],
             [_p(f"<b>{g.label}</b>", st["cell_soft"])],
             [_p(g.display or "—", st["fine"])]],
            colWidths=[CONTENT_WIDTH / 7 - 2 * mm], hAlign="LEFT",
            style=TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                              ("TOPPADDING", (0, 0), (-1, -1), 0.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 0.6)]),
        ))
    per_row = 7
    rows = [cells[i:i + per_row] for i in range(0, len(cells), per_row)]
    for r in rows:
        while len(r) < per_row:
            r.append("")
    grid = Table(rows, colWidths=[CONTENT_WIDTH / per_row] * per_row, hAlign="LEFT")
    grid.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 2), ("TOPPADDING", (0, 0), (-1, -1), 2),
                              ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
    story.append(grid)
    verdict = {
        "pass": "All gates passed: a non-detection can be read as evidence of absence where the detection model says so.",
        "warn": "One or more gates warned: readings are shown, non-detections are qualified.",
        "fail": "A gate failed: readings are shown for completeness, but no non-detection is treated as absence and action language is withheld.",
    }.get(overall, "Gate status unknown.")
    story.append(_p(
        f"<b>Overall: {overall.upper()}.</b> {verdict} {meta.get('read_pairs', '')} read pairs were analysed, every one of them.",
        st["small"],
    ))


def community_metric_strips(
    story: list[Any], st: dict[str, ParagraphStyle], *, similarity: Any | None, inv: Any = None,
    primary_phyla: dict[str, float] | None = None,
) -> None:
    """Diversity, evenness, richness and the F:B ratio against the reference — neutral by design."""
    story.append(_p("COMMUNITY-WIDE METRICS · NEUTRAL RANGE STRIPS", st["label"]))
    taxonomy = similarity.taxonomy if similarity is not None else None
    if taxonomy is None:
        story.append(_p("Species-level composition was not measured, so there are no community-wide metrics.", st["small"]))
        return
    from openbiota.scoring import ecological_metrics
    from openbiota.similarity import rescale_to_classified

    species = {k: v for k, v in rescale_to_classified(taxonomy.species).items() if v > 0}
    metrics = ecological_metrics(species)
    dists = _cohort_metric_distributions(similarity.matched_cohort) if similarity.matched_cohort is not None else {}
    n_ref = similarity.matched_cohort.n_samples if similarity.matched_cohort is not None else 0

    # Same phyla as the tile and the ring above, so the page has one F:B ratio.
    if primary_phyla:
        tot = sum(primary_phyla.values()) or 1.0
        phyla = sorted(((k, 100.0 * v / tot) for k, v in primary_phyla.items()), key=lambda kv: -kv[1])
    else:
        phyla = sorted(rescale_to_classified(taxonomy.phyla).items(), key=lambda kv: -kv[1])
    firm = next((v for k, v in phyla if k.lower().startswith(("firmicutes", "bacillota"))), 0.0)
    bact = next((v for k, v in phyla if k.lower().startswith(("bacteroidetes", "bacteroidota"))), 0.0)
    fb = firm / bact if bact else None

    entries = [
        ("Shannon diversity", f"{metrics['shannon_diversity']:.2f}", metrics["shannon_diversity"], dists.get("shannon_diversity"),
         "How many species there are and how evenly abundance is spread among them. Higher is more diverse; it is a description, not a grade."),
        ("Evenness", f"{metrics['evenness']:.2f}", metrics["evenness"], dists.get("evenness"),
         "0 to 1: whether a few species dominate (low) or many share the community (high)."),
        ("Richness, reference-comparable", f"{int(metrics['richness']):,}", metrics["richness"], dists.get("richness"),
         "The count of organisms on the same catalogue the reference adults were measured on, which is what "
         "makes the percentile valid. It is not the total found: pooling every catalogue finds more, and "
         "that total is the headline count on page one and in the organism section."),
    ]
    bar_w = CONTENT_WIDTH * 0.40
    widths = [CONTENT_WIDTH * 0.22, CONTENT_WIDTH * 0.10, bar_w, CONTENT_WIDTH * 0.08, CONTENT_WIDTH * 0.20]
    rows: list[list[Any]] = [[_p("METRIC", st["label"]), _p("YOU", st["label"]), _p("WHERE YOU SIT AMONG THE REFERENCE", st["label"]),
                              _p("", st["label"]), _p("READING", st["label"])]]
    style_extra: list[Any] = []
    for label, shown, value, dist, note in entries:
        pct = _pct_of(dist, value)
        status = status_for(percentile=pct, higher_means="unclear")
        rows.append([
            _p(f"<b>{label}</b>", st["cell"]),
            _p(f"<b>{shown}</b>", st["cell"]),
            PercentileBar(width=bar_w - 6 * mm, percentile=pct, higher_means="unclear", marker_colour=status.colour, height=6.4 * mm),
            _p(f"<font color='{_hex(status.colour)}'>{_ordinal(pct)}</font>", st["pct"]),
            StatusChip(status, width=widths[4] - 6 * mm, font_size=6.4),
        ])
        n = len(rows)
        rows.append([_p(note, st["glance_note"]), "", "", "", ""])
        style_extra += [("SPAN", (0, n), (-1, n)), ("TOPPADDING", (0, n), (-1, n), 0), ("BOTTOMPADDING", (0, n), (-1, n), 3.5)]
    if inv is not None and len(inv):
        # Diversity over everything detected. No percentile: no reference
        # cohort exists on the pooled catalogues, and a percentile against a
        # cohort measured another way would be a number that means nothing.
        from openbiota.scoring import ecological_metrics as _eco

        pooled = {o.species: o.percent for o in inv.in_composition if o.percent > 0}
        if len(pooled) >= 2:
            m_all = _eco(pooled)
            rows.append([
                _p("<b>Diversity, all organisms</b>", st["cell"]),
                _p(f"<b>{m_all['shannon_diversity']:.2f}</b>", st["cell"]),
                _p(f"<font color='{_hex(INK_FAINT)}'>Shannon index over the {len(pooled):,} organisms in the "
                   "full composition. Not placed on a scale: no reference population has been measured this "
                   "way yet, so there is nothing to rank it against.</font>", st["glance_note"]),
                "", _chip("not ranked", SLATE, SLATE_BG, widths[4] - 6 * mm),
            ])
            n_all = len(rows) - 1
            style_extra += [("SPAN", (2, n_all), (3, n_all))]
    rows.append([
        _p("<b>Firmicutes : Bacteroidetes</b>", st["cell"]),
        _p("<b>—</b>" if fb is None else f"<b>{fb:.2f}</b>", st["cell"]),
        _p(f"<font color='{_hex(INK_FAINT)}'>not placed on a scale — the ratio has no agreed healthy range and "
           "depends strongly on the profiler and the population</font>", st["glance_note"]),
        "", _chip("not judged", SLATE, SLATE_BG, widths[4] - 6 * mm),
    ])
    n = len(rows) - 1
    style_extra += [("SPAN", (2, n), (3, n))]
    story.append(_plain(rows, widths, extra=style_extra))
    story.append(_p(
        f"Reference: {n_ref:,} matched adult stool metagenomes, each profiled with the same method; the bar's centre band is "
        "their middle half. Grey means neutral by design: none of these three numbers has a direction that research reads "
        "as favourable on its own.",
        st["fine"],
    ))


#: The major at-a-glance destinations page 3 offers, in reading order. Only
#: overview sections: §12.2 is explicit that this is not an index of every
#: page, and each entry points at a section's first overview page rather than
#: at its detail. Entries are generated from the stable anchors in `SECTIONS`,
#: so a renamed or renumbered section cannot leave a dangling row here.
#: The contents, in page order, each title matching the heading the section
#: actually prints. A title here that differs from the heading it links to
#: sends a reader somewhere they did not choose, so
#: `tests/test_contents_matches_the_report.py` compares the two.
CONTENTS_ENTRIES: Final[tuple[tuple[str, str, str], ...]] = (
    ("summary", "Your results at a glance", "every score on one page"),
    ("stood_out", "What stood out", "the findings worth your attention first"),
    ("community", "Your gut community", "diversity, composition and community type"),
    ("groups", "Your microbial groups", "curated sets of species, each summed to one reading"),
    ("organisms", "Organisms that need attention", "the ones flagged, one by one"),
    ("catalogue", "Every organism found", "the full inventory, searchable"),
    ("pathogens", "Pathogen screening", "what was searched for, and what was found"),
    ("age", "Estimated biological age of biota", "how your community compares by age"),
    ("functions", "Your Microbial Functions",
     "what your microbes can make, break down and transform"),
    ("patterns", "Resemblance to published disease patterns",
     "how your community compares with published patterns"),
    ("biofilm", "Biofilm-related potential",
     "organisms that build and degrade the mucus layer"),
    ("mycobiome", "Your gut mycobiome", "the fungal side of the community"),
    ("skin", "Gut-skin axis", "one score for how your gut leans across published skin panels"),
    ("actions", "What you can do about it", "options with the evidence behind each one"),
    ("context", "Your information & report context", "what was supplied, and what was assumed"),
)


def _accent_hex() -> str:
    from openbiota.pdfreport import ACCENT  # noqa: PLC0415

    return _hex(ACCENT)


def _contents_table(story: list[Any], st: dict[str, ParagraphStyle]) -> None:
    """A two-column contents whose whole row is clickable.

    Page numbers are deliberately absent: they are not known when this page
    is laid out, and a wrong number is worse than none when the row itself
    is the link. The reader's bookmarks panel carries the pagination.
    """
    from openbiota.pdfreport import SECTIONS as _SECTIONS  # noqa: PLC0415

    colour = _accent_hex()
    half = CONTENT_WIDTH * 0.5
    cells: list[Any] = []
    for key, title, blurb in CONTENTS_ENTRIES:
        number = _SECTIONS.get(key)
        if number is None:  # a section this build does not carry
            continue
        # The whole cell is the target, not just the words in it: a reader
        # aiming at a two-line entry should not have to hit the text.
        cells.append(linked_block(
            section_dest(key),
            [
                _p(f'<a href="#{section_dest(key)}" color="{colour}">'
                   f'<b>{(str(number) + " &nbsp; ") if number else ""}{title}</b></a>', st["cell"]),
                _p(f'<font size="6.6" color="{_hex(INK_FAINT)}">{blurb}</font>', st["cell"]),
            ],
            half,
        ))
    rows = [cells[i:i + 2] for i in range(0, len(cells), 2)]
    if rows and len(rows[-1]) == 1:
        rows[-1].append("")
    # Ample row spacing, per the layout note: the row is the link, so it
    # wants to be a comfortable target rather than a dense line of text.
    story.append(_plain(
        rows, [half, half], header=False,
        extra=[("TOPPADDING", (0, 0), (-1, -1), 5.0),
               ("BOTTOMPADDING", (0, 0), (-1, -1), 5.4)],
    ))


def reading_guide_page(
    story: list[Any], st: dict[str, ParagraphStyle], *, section: int,
) -> None:
    """One page: how to read every number here, then where to find them.

    Spec 0.8.3 \u00a712.2 empties the lower half of this page. The two-engine
    explanation moves to the methods material, the sample-QC strip merges
    with the sequencing-quality section that already prints the same gate
    table, and everything the report knew about the reader moves to its own
    section. What replaces them is a contents table that goes somewhere.
    """
    from openbiota.pdfreport import BODY_FRAME_HEIGHT, ScaleLegend

    outer = story
    story = []
    story.extend(section_heading(section, "Table of contents", st["h1"]))

    # ---- explore your results, and the contents ------------------------- #
    # This page leads with where to go. The explanation of the scale is
    # what a reader needs second, not first, so it sits under the contents
    # in smaller type and takes the lower half of the page.
    story.append(_p(
        "This report is interactive. Click a section for results at a glance, then a reading "
        "for its detailed explanation. Use \u2018Back to overview\u2019 to return.",
        st["body"],
    ))
    story.append(Spacer(1, 1.5 * mm))
    _contents_table(story, st)
    story.append(Spacer(1, 3 * mm))

    # ---- how to read it, as a subheading -------------------------------- #
    story.append(_p("<font size='8.5'><b>How to read this report</b></font>", st["h3"]))
    story.append(_p(
        "<font size='6.8'><b>One scale, everywhere.</b> Every reading \u2014 organism, group, "
        "metabolic function, disease pattern \u2014 is a <b>percentile against a reference "
        "group</b>, translated into the same seven words. The 60th percentile means 60% of the "
        "reference had a lower value than you; the reference\u2019s middle half (25th to 75th) "
        "is <b>typical</b>.</font>",
        st["small"],
    ))
    story.append(Spacer(1, 1 * mm))
    story.append(ScaleLegend(width=CONTENT_WIDTH * 0.80, higher_means="adverse"))
    story.append(Spacer(1, 1 * mm))
    story.append(_p(
        "<font size='6.8'><b>Colour is the judgement; the word is the position.</b> Where "
        "research links higher with worse outcomes (\u25b2) the high end is red; where higher is "
        "favourable (\u25bc) the bar is mirrored and high is green; where the evidence is mixed "
        "the bar is grey and no verdict is given. Disease-pattern resemblance always reads higher "
        "as less favourable. Everything here is <i>capacity</i> or <i>composition</i> read from "
        "DNA \u2014 never a measured amount of a compound in your body, and never a "
        "diagnosis.</font>",
        st["small"],
    ))

    outer.append(KeepInFrame(CONTENT_WIDTH, BODY_FRAME_HEIGHT - 2 * mm, story, mode="shrink", hAlign="LEFT"))


# --------------------------------------------------------------------------- #
# Part A · groups at a glance
# --------------------------------------------------------------------------- #


def groups_glance(story: list[Any], st: dict[str, ParagraphStyle], *, section: int, similarity: Any | None,
                  cohort_note: str, inv: Any | None = None, carded: frozenset[str] = frozenset()) -> None:
    from openbiota.pdflibrary import _arrow, group_driver_line, group_drivers
    from openbiota.taxongroups import GROUP_CATEGORIES

    outer = story
    # The heading carries this section's link destination, so it goes on the
    # real story: a KeepInFrame that shrinks its contents drops zero-height
    # flowables, and the contents page links here.
    outer.extend(section_heading(section, "Your microbial groups at a glance", st["h1"]))
    story = []
    story.append(_p(
        "Curated sets of species that share a job or an origin, each summed to one reading "
        f"and placed against {cohort_note}. \u25bc higher is favourable, \u25b2 higher is "
        "adverse, no arrow: it depends on context. Under each reading are the organisms "
        "carrying it and the common members that were not found.",
        st["body"],
    ))
    community = similarity.community if similarity is not None else None
    if community is None:
        outer.append(_p("The taxonomic engine did not run; no group readings.", st["small"]))
        return
    bar_w = CONTENT_WIDTH * 0.36
    widths = [CONTENT_WIDTH * 0.30, bar_w, CONTENT_WIDTH * 0.08, CONTENT_WIDTH * 0.26]
    rows: list[list[Any]] = [[_p("GROUP", st["label"]), _p("WHERE YOU SIT", st["label"]), _p("", st["label"]), _p("READING", st["label"])]]
    extra: list[Any] = []
    by_cat = community.by_category()
    for key, heading in GROUP_CATEGORIES:
        groups = by_cat.get(key, [])
        if not groups:
            continue
        n = len(rows)
        rows.append([_p(heading.upper(), st["label"]), "", "", ""])
        extra += [("SPAN", (0, n), (-1, n)), ("BACKGROUND", (0, n), (-1, n), PANEL_BG)]
        for g in sorted(groups, key=lambda r: -(abs((r.percentile or 50) - 50))):
            status = status_for(percentile=g.percentile, higher_means=g.group.effective_direction,
                                detected=g.detected or g.percentile is not None)
            if not g.detected and g.percentile is not None:
                status = Status(f"not detected · {status.label}", status.colour, status.background, status.note, status.level)
            rows.append([
                linked_cell(
                    RowLink(href=group_dest(g.group.name, detail=True), anchor=group_dest(g.group.name)),
                    _p(f"<b>{g.group.label}</b>{_arrow(g.group.effective_direction)}<br/>"
                       f"<font color='{_hex(INK_FAINT)}' size='6.4'>{g.percent:.2f}% of reads · {g.n_detected} of {g.n_resolvable} species</font>",
                       st["cell"]),
                ),
                PercentileBar(width=bar_w - 6 * mm, percentile=g.percentile, higher_means=g.group.effective_direction,
                              marker_colour=status.colour, height=6.4 * mm),
                _p(f"<font color='{_hex(status.colour)}'>{_ordinal(g.percentile)}</font>", st["pct"]),
                StatusChip(status, width=widths[3] - 6 * mm, font_size=6.4),
            ])
            # The organisms carrying the reading, on a full-width line of
            # their own directly under it - readable, one line, never cut.
            # The same resolver feeds the detail card, so the names agree.
            present, absent = group_drivers(g, inv=inv, carded=carded)
            n = len(rows)
            rows.append([_p(group_driver_line(present, absent), st["cell"]), "", "", ""])
            extra += [
                ("SPAN", (0, n), (-1, n)),
                ("TOPPADDING", (0, n), (-1, n), 0), ("BOTTOMPADDING", (0, n), (-1, n), 4),
                ("LINEBELOW", (0, n - 1), (-1, n - 1), 0, colors.white),
            ]
    story.append(_plain(rows, widths, extra=extra))
    # One page, always: the group count and the length of the driver lines
    # vary by sample, and this section was spilling two lines onto a second
    # page rather than sitting a shade tighter on one.
    from openbiota.pdfreport import BODY_FRAME_HEIGHT  # noqa: PLC0415

    outer.append(KeepInFrame(CONTENT_WIDTH, BODY_FRAME_HEIGHT - 2 * mm, story,
                             mode="shrink", hAlign="LEFT"))


# --------------------------------------------------------------------------- #
# organism polarity words and colours (shared with the summary page)
# --------------------------------------------------------------------------- #


#: Curated polarity collapsed to the three words a reader needs.
POLARITY_WORD = {
    "health_associated": "beneficial",
    "disease_associated": "adverse",
    "potential_pathobiont": "adverse",
    "validated_pathogen_signal": "adverse",
    "context_dependent": "neutral",
    "conflicting_evidence": "neutral",
    "no_curated_interpretation": "unknown",
}


#: Taxon groups whose membership carries a direction when the species itself
#: has no curated verdict. A butyrate producer or a probiotic species that is
#: uncommon in the reference is good news, not a warning.
#: Group names as declared in taxa/*.yaml (``name:``), not file names. Groups
#: whose ``higher_means`` is favourable, plus the probiotic genera whose
#: presence in stool is at worst neutral and whose species have human trial
#: evidence; and the one group whose ``higher_means`` is adverse.
BENEFICIAL_GUILDS = frozenset({
    "beneficial", "butyrate", "bifidobacteria", "probiotics", "polyphenol_metabolisers",
})
ADVERSE_GUILDS = frozenset({"opportunistic_pathogens"})


def polarity_word(f: TaxonFinding) -> str:
    """beneficial, adverse, neutral or unknown — the curated polarity first, then the guild."""
    word = POLARITY_WORD.get(getattr(f.interpretation, "polarity", ""), "unknown")
    if word in ("beneficial", "adverse"):
        return word
    guilds = set(getattr(f, "guilds", ()) or ())
    if guilds & BENEFICIAL_GUILDS:
        return "beneficial"
    if guilds & ADVERSE_GUILDS:
        return "adverse"
    return word


def polarity_judgement(f: TaxonFinding) -> tuple[colors.Color, colors.Color]:
    """Colour is the judgement: does this organism's *position* read well or badly?

    Position (missing, low, high, unusual) is level; the curated polarity says
    whether that level is good news. A beneficial species missing is red; the
    same species unusually abundant is green. An adverse species unusually
    present is red; missing, green. Neutral or unknown organisms are grey
    whatever their level — the report has no verdict, and says so. Bucket 6
    (needs qualification) is always red-flagged because it exists only for
    concern-polarity species.
    """
    word = polarity_word(f)
    if f.bucket == 6:
        return CORAL, CORAL_BG
    if f.bucket == 1 or f.bucket is None or word in ("neutral", "unknown"):
        return SLATE, SLATE_BG
    lower_side = f.bucket in (3, 4)
    if word == "beneficial":
        if lower_side:
            return (CORAL, CORAL_BG) if f.bucket == 4 else (AMBER, AMBER_BG)
        return GREEN, GREEN_BG
    # adverse
    if lower_side:
        return GREEN, GREEN_BG
    return CORAL, CORAL_BG


# --------------------------------------------------------------------------- #
# Part A · disease patterns at a glance
# --------------------------------------------------------------------------- #


def profiles_glance(story: list[Any], st: dict[str, ParagraphStyle], *, section: int, similarity: Any | None) -> None:
    from openbiota.pdflibrary import _ranked_table, _shape_block

    story.extend(section_heading(section, "Resemblance to published disease patterns at a glance", st["h1"]))
    story.append(_p(
        "For every pattern in the library: how closely this community resembles the group-level pattern published for "
        "people with that condition, as a percentile of the reference adults — the same scale as everything else. "
        "Resemblance is not a diagnosis and not a probability of disease. Every scored pattern has a full page, and the "
        "research on what has been tried for that condition, in Part B.",
        st["body"],
    ))
    if similarity is None or not similarity.results:
        story.append(_p("Profile similarity was not computed for this run.", st["small"]))
        return
    _shape_block(story, st, similarity)
    story.append(Spacer(1, 2 * mm))
    _ranked_table(story, st, similarity)


# --------------------------------------------------------------------------- #
# Part A · evidence overview
# --------------------------------------------------------------------------- #


def _card_status_chip(s: EvidenceSummary, width: float) -> StatusChip:
    if s.display_policy == "suppressed":
        return _chip("not shown", SLATE, SLATE_BG, width)
    if s.render == "no_supported_targeted_intervention":
        return _chip("no supported intervention", SLATE, SLATE_BG, width, size=6.0)
    if s.render == "guideline_route_after_confirmation":
        return _chip("clinical route", ACCENT_DARK, PANEL_BG, width)
    if s.render == "clinician_review":
        return _chip("clinician review", ORANGE, ORANGE_BG, width)
    if s.supporting_evidence_lane in ("A", "B"):
        return _chip("human evidence · strong", GREEN, GREEN_BG, width, size=6.0)
    if s.supporting_evidence_lane == "C":
        return _chip("human trial · limited", AMBER, AMBER_BG, width, size=6.0)
    if s.supporting_evidence_lane == "D":
        return _chip("human · non-randomised", ORANGE, ORANGE_BG, width, size=6.0)
    return _chip("no human trial", SLATE, SLATE_BG, width)


def _lane_chip(lane: str, width: float) -> StatusChip:
    col, bg = LANE_COLOUR.get(lane, (SLATE, SLATE_BG))
    return _chip(f"lane {lane}" if lane != "none" else "none", col, bg, width, size=6.4)


def evidence_overview(story: list[Any], st: dict[str, ParagraphStyle], *, section: int, plan: Any | None,
                      as_subsection: bool = False) -> None:
    """The intervention registry's cards, one row each.

    With ``as_subsection`` the section heading is assumed already written
    (by the levers page that opens the actions section) and this renders under an
    h2 instead.
    """
    if as_subsection:
        story.append(CondPageBreak(70 * mm))
        story.append(Spacer(1, 4 * mm))
        story.append(_p("Where a guideline, clinical route or exact-product trial exists", st["h2"]))
    else:
        story.extend(section_heading(section, "What the research says you can do about it — overview", st["h1"]))
    story.append(_p(
        "Every reading in this report that has a research record about <i>changing</i> it is listed here with the "
        "strength of that record: <b>lane A</b> a current guideline or regulatory label, <b>B</b> replicated randomised "
        "benefit for the exact product, <b>C</b> a single or conflicting randomised trial, <b>D</b> non-randomised human "
        "evidence, <b>E</b> animal or laboratory only. <b>Against</b> counts the null, harm or guideline-against records "
        "shown alongside — they are never averaged away. Part B has the full card for each: the exact product or diet "
        "studied, in whom, the studied protocol, and the safety rules that apply. Nothing here is a treatment score, a "
        "probability of benefit, or a prescription.",
        st["body"],
    ))
    if plan is None or not plan.summaries:
        story.append(_p("No evidence cards were generated for this run.", st["small"]))
        return
    order = {"guideline_route_after_confirmation": 0, "clinician_review": 1, "evidence_only": 2, "inform_only": 3,
             "not_actionable": 4, "no_supported_targeted_intervention": 5}
    lane_rank = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "none": 5, "X": 6}
    summaries = sorted(plan.summaries, key=lambda s: (order.get(s.render, 9), lane_rank.get(s.supporting_evidence_lane, 9), s.title))
    widths = [CONTENT_WIDTH * 0.27, CONTENT_WIDTH * 0.25, CONTENT_WIDTH * 0.08, CONTENT_WIDTH * 0.09, CONTENT_WIDTH * 0.09, CONTENT_WIDTH * 0.22]
    rows: list[list[Any]] = [[_p("READING THAT FIRED", st["label"]), _p("WHAT WAS STUDIED", st["label"]), _p("LANE", st["label"]),
                              _p("FOR / AGAINST", st["label"]), _p("GUIDELINE", st["label"]), _p("STATUS", st["label"])]]
    for s in summaries:
        trig = "; ".join(s.trigger_displays[:2]) + (f"; +{len(s.trigger_displays) - 2} more" if len(s.trigger_displays) > 2 else "")
        what = s.intervention_display or ("clinical route" if s.follow_up_route else "—")
        n_for = len(s.supporting)
        n_against = len(s.against)
        guide = s.guideline_position.replace("_", " ") if s.guideline_position != "unavailable" else "—"
        rows.append([
            _p(f"<b>{s.title}</b><br/><font color='{_hex(INK_FAINT)}' size='6.2'>{trig}</font>", st["cell"]),
            _p(f"{what}" + (f"<br/><font color='{_hex(INK_FAINT)}' size='6.2'>{(s.intervention_class or '').replace('_', ' ')}</font>" if s.intervention_class else ""), st["cell"]),
            _lane_chip(s.supporting_evidence_lane, widths[2] - 5 * mm),
            _p(f"<font color='{_hex(GREEN)}'><b>{n_for}</b></font> / <font color='{_hex(CORAL)}'><b>{n_against}</b></font>", st["cell"]),
            _p(guide, st["cell_soft"]),
            _card_status_chip(s, widths[5] - 6 * mm),
        ])
    story.append(_plain(rows, widths, zebra=True))
    story.append(Spacer(1, 2 * mm))
    story.append(_p(
        f"{len(summaries)} cards from registry {plan.registry_lock.get('registry_version')} "
        f"({plan.registry_lock.get('assertions')} evidence records, {plan.registry_lock.get('interventions')} exact interventions, "
        f"{plan.registry_lock.get('safety_rules')} safety rules). {plan.n_with_human_trials} rest on human trial or guideline "
        f"evidence, {plan.n_clinical_routes} point to a conventional clinical route, and {plan.n_no_supported} say plainly that no "
        "supported targeted intervention exists — which is itself a result. Where fecal microbiota transplantation is the only "
        "route with evidence, the card says so and gives its regulatory status.",
        st["fine"],
    ))


# --------------------------------------------------------------------------- #
# evidence cards (Part B, under each reading)
# --------------------------------------------------------------------------- #


def _assertion_line(a: Mapping[str, Any], st: dict[str, ParagraphStyle], *, against: bool) -> Paragraph:
    lane = str(a.get("lane", ""))
    col, _bg = LANE_COLOUR.get(lane, (SLATE, SLATE_BG))
    design = a.get("design") or a.get("source_type") or ""
    n = a.get("participants")
    eff = a.get("effect_direction", "")
    eff_col = {"benefit": GREEN, "null": CORAL, "harm": CORAL, "mixed": AMBER}.get(eff, INK_FAINT)
    meta_bits = [b for b in (
        f"lane {lane}" if lane else "", str(design).replace("_", " "), f"n = {n}" if n else "",
        f"<font color='{_hex(eff_col)}'>{eff.replace('_', ' ')}</font>" if eff and eff != "not_estimable" else "",
    ) if b]
    url = a.get("url") or ""
    cite = a.get("citation") or url
    cite_html = f'<a href="{url}" color="{_hex(ACCENT_DARK)}"><u>{cite}</u></a>' if url else cite
    prefix = f"<font color='{_hex(CORAL)}'><b>Against or limiting —</b></font> " if against else ""
    return _p(
        f"<font color='{_hex(col)}'><b>\u25a0</b></font>&nbsp; {prefix}{a.get('summary', '')} "
        f"<font color='{_hex(INK_FAINT)}' size='6.6'>[{cite_html}] · {' · '.join(meta_bits)}</font>",
        st["evidence"],
    )


def evidence_cards(story: list[Any], st: dict[str, ParagraphStyle], cards: Sequence[EvidenceSummary], *,
                   heading: str = "What the research says about changing this") -> None:
    """Render every evidence card for one reading."""
    if not cards:
        return
    block: list[Any] = [Spacer(1, 1 * mm), _p(heading, st["h3_accent"])]
    cards_out: list[list[Any]] = []
    for s in cards:
        head = Table([[
            _p(f"<b>{s.title}</b>", st["cell_bold"]),
            _lane_chip(s.supporting_evidence_lane, 17 * mm),
            _card_status_chip(s, 36 * mm),
        ]], colWidths=[CONTENT_WIDTH - 17 * mm - 36 * mm - 4 * mm, 17 * mm + 2 * mm, 36 * mm + 2 * mm], hAlign="LEFT")
        head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                                  ("RIGHTPADDING", (0, 0), (-1, -1), 2), ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
                                  ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                                  ("LEFTPADDING", (0, 0), (0, 0), 4)]))
        card: list[Any] = [head]
        fired = "; ".join(s.trigger_displays[:3]) + (f"; +{len(s.trigger_displays) - 3} more" if len(s.trigger_displays) > 3 else "")
        card.append(_p(f"<font color='{_hex(INK_FAINT)}'>Fired by:</font> {fired}", st["fine"]))
        if s.display_policy == "suppressed":
            for para in s.participant_paragraphs:
                card.append(_p(para, st["small"]))
            cards_out.append(card)
            continue
        if s.context_note:
            card.append(_p(" ".join(s.context_note.split()), st["body_ink"]))
        if s.render == "no_supported_targeted_intervention":
            card.append(_p(
                "<b>No supported targeted intervention</b> — no human study on record shows that changing this "
                "reading changes a health outcome. That is a result, not a gap in the analysis.", st["small_ink"]))
        if s.intervention_display:
            lane_text = LANE_MEANING.get(s.supporting_evidence_lane, "")
            pieces = [f"<b>Studied:</b> {s.intervention_display} ({(s.intervention_class or '').replace('_', ' ')})"]
            if s.population_text:
                pieces.append(f"<b>in:</b> {s.population_text}")
            if s.studied_protocol_text:
                pieces.append(f"<b>studied protocol (a source fact, not an instruction):</b> {s.studied_protocol_text}")
            card.append(_p(" &nbsp;·&nbsp; ".join(pieces), st["small_ink"]))
            card.append(_p(f"<font color='{_hex(INK_FAINT)}'>Evidence lane {s.supporting_evidence_lane}: {lane_text}. "
                           f"Replication: {s.replication.replace('_', ' ')}."
                           + (f" Strongest design: {s.strongest_study_design.replace('_', ' ')}." if s.strongest_study_design else "")
                           + (f" Guideline position: {s.guideline_position.replace('_', ' ')}." if s.guideline_position != 'unavailable' else "")
                           + "</font>", st["fine"]))
        for a in s.supporting:
            card.append(_assertion_line(a, st, against=False))
        for a in s.against:
            card.append(_assertion_line(a, st, against=True))
        if s.eligibility_explanation:
            card.append(_p(f"<font color='{_hex(ORANGE)}'><b>Applicability:</b></font> " + " ".join(s.eligibility_explanation), st["small"]))
        if s.follow_up_route:
            route = re.sub(r"\s*\((?:CLINICAL_ROUTE_)?[A-Z0-9_]+\)\.?$", ".", s.follow_up_route).strip()
            card.append(_p(f"<b>Conventional route:</b> {route}", st["small_ink"]))
        if s.fmt_status:
            fmt_words = {
                "established_recurrent_cdi_only": "established only for recurrent C. difficile infection",
                "investigational_non_cdi": "investigational for this condition (IND-only; no treatment route from this report)",
                "guideline_against_outside_trials": "guidelines recommend against it for this condition outside trials",
            }.get(s.fmt_status, s.fmt_status.replace("_", " "))
            card.append(_p(f"<b>Fecal microbiota transplantation:</b> {fmt_words}.", st["small_ink"]))
        if s.intervention_class == "prescription_drug":
            card.append(_p("Prescription route for a clinician to review — never a recommendation generated from a stool result.", st["small"]))
        card.append(_p(f"<font color='{_hex(INK_FAINT)}'>Safety:</font> {participant_safety_sentence(s.safety)}"
                       + (f" {' '.join(list(s.safety.reasons)[:2])}" if s.safety.reasons else ""), st["fine"]))
        if s.never_say:
            card.append(_p(f"<font color='{_hex(INK_FAINT)}'>This report will not say:</font> "
                           + " · ".join(f"\u201c{w}\u201d" for w in s.never_say), st["fine"]))
        card.append(Spacer(1, 1.2 * mm))
        cards_out.append(card)
    # The closing caveat is one line and belongs to the last card. Left
    # loose it lands alone on the next page, which reads as a defect.
    caveat = _p(
        "Not a treatment score, a probability of benefit or a personalised prescription. "
        "Discuss any change with a clinician or registered dietitian.", st["fine"])
    if cards_out:
        cards_out[-1] = [*cards_out[-1], caveat]
    else:
        block.append(caveat)
    # The heading travels with the first card, and each card is kept whole:
    # a card split across a page break loses the tie between its title, its
    # status and the evidence underneath it.
    story.append(KeepTogether([*block, *(cards_out[0] if cards_out else [])]))
    for card in cards_out[1:]:
        story.append(KeepTogether(card))


# --------------------------------------------------------------------------- #
# age page
# --------------------------------------------------------------------------- #


class AgeScale(Flowable):
    """Predicted age on an 18–80 axis over the training-age support histogram,
    with 80% and 95% intervals and the actual age when supplied."""

    def __init__(self, *, width: float, predicted: float | None, pi80: tuple[float, float] | None,
                 pi95: tuple[float, float] | None, actual: float | None, support_edges: Sequence[float],
                 support_counts: Sequence[int], lo: float = 15.0, hi: float = 85.0,
                 compact: bool = False, height: float | None = None,
                 pi90: tuple[float, float] | None = None) -> None:
        super().__init__()
        self.width = width
        self.compact = compact
        # An explicit height lets the caller size the chart to the room it
        # has; the histogram grows to fill it, everything else stays put.
        self.height = height if height is not None else (24 * mm if compact else 34 * mm)
        self.predicted = predicted
        self.pi80 = pi80
        self.pi95 = pi95
        self.pi90 = pi90
        self.actual = actual
        self.edges = list(support_edges)
        self.counts = list(support_counts)
        self.lo, self.hi = lo, hi

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def _x(self, age: float) -> float:
        f = (min(max(age, self.lo), self.hi) - self.lo) / (self.hi - self.lo)
        return 6 * mm + f * (self.width - 12 * mm)

    def draw(self) -> None:
        c: Canvas = self.canv
        base = 6 * mm if self.compact else 9 * mm
        # Room above the histogram for the interval bar and its number.
        head_room = (3 * mm + 5 * mm) if self.compact else (4 * mm + 6 * mm)
        hist_h = max(6 * mm, self.height - base - head_room)
        # support histogram
        if self.counts:
            peak = max(self.counts) or 1
            c.setFillColor(colors.Color(0.86, 0.90, 0.94))
            for i, n in enumerate(self.counts):
                if i + 1 >= len(self.edges) or n == 0:
                    continue
                x0, x1 = self._x(self.edges[i]), self._x(self.edges[i + 1])
                c.rect(x0 + 0.4, base, max(0.0, x1 - x0 - 0.8), hist_h * n / peak, stroke=0, fill=1)
        # axis
        c.setStrokeColor(RULE)
        c.setLineWidth(0.6)
        c.line(self._x(self.lo), base, self._x(self.hi), base)
        c.setFont("Helvetica", 6.0)
        c.setFillColor(INK_FAINT)
        first_tick = 10 * math.ceil(self.lo / 10.0)
        for age in range(first_tick, int(self.hi) + 1, 10):
            x = self._x(age)
            c.line(x, base, x, base - 1.2 * mm)
            c.drawCentredString(x, base - 4.2 * mm, str(age))
        if not self.compact:
            c.drawString(self._x(self.lo), base - 7.6 * mm, "years · shaded: how many training adults at each age")
        # intervals
        y = base + hist_h + (3 * mm if self.compact else 4 * mm)
        if self.pi95 is not None:
            c.setStrokeColor(colors.Color(0.62, 0.70, 0.78))
            c.setLineWidth(2.2)
            c.line(self._x(self.pi95[0]), y, self._x(self.pi95[1]), y)
        if self.pi80 is not None:
            c.setStrokeColor(ACCENT)
            c.setLineWidth(4.0)
            c.line(self._x(self.pi80[0]), y, self._x(self.pi80[1]), y)
        if self.pi90 is not None:
            # The band the health adjustment is held within: drawn as one
            # bar with round ends, in the accent, so it reads as "the range".
            c.saveState()
            c.setStrokeColor(ACCENT)
            c.setLineCap(1)
            c.setLineWidth(4.0)
            c.line(self._x(self.pi90[0]), y, self._x(self.pi90[1]), y)
            c.restoreState()
        if self.predicted is not None:
            x = self._x(self.predicted)
            c.setFillColor(INK)
            c.setStrokeColor(colors.white)
            c.setLineWidth(1.2)
            c.circle(x, y, 3.2, stroke=1, fill=1)
            c.setFont("Helvetica-Bold", 8.4)
            c.drawCentredString(x, y + 5.2, f"{self.predicted:.0f}")
        if self.actual is not None:
            x = self._x(self.actual)
            c.setStrokeColor(CORAL)
            c.setLineWidth(1.4)
            c.line(x, base, x, y + 3)
            c.setFillColor(CORAL)
            c.setFont("Helvetica-Bold", 6.4)
            c.drawCentredString(x, y + 9.5, f"actual {self.actual:.0f}")


def _json_field(meta: Mapping[str, Any], key: str) -> Any:
    v = meta.get(key)
    if isinstance(v, str):
        try:
            return json.loads(v)
        except ValueError:
            return v
    return v


def age_page(story: list[Any], st: dict[str, ParagraphStyle], *, section: int, age: Any | None,
             age_meta: Mapping[str, Any] | None, subject_age: float | None, mode: str) -> None:
    """The age section, on one page.

    The heading carries the section's link destination and stays outside
    the shrink-to-fit frame, which drops zero-height flowables; everything
    else is built into its own list and shrunk to the page if it needs it.
    """
    from reportlab.platypus import KeepInFrame  # noqa: PLC0415

    from openbiota.pdfreport import BODY_FRAME_HEIGHT  # noqa: PLC0415

    story.extend(section_heading(section, "Estimated biological age of biota", st["h1"]))
    body: list[Any] = []
    _age_page_body(body, st, age=age, age_meta=age_meta, subject_age=subject_age, mode=mode)
    story.append(KeepInFrame(CONTENT_WIDTH, BODY_FRAME_HEIGHT - 14 * mm, body, mode="shrink", hAlign="LEFT"))


def _age_page_body(story: list[Any], st: dict[str, ParagraphStyle], *, age: Any | None,
                   age_meta: Mapping[str, Any] | None, subject_age: float | None, mode: str) -> None:  # noqa: ARG001
    # Read the cohort size and study count from the bundle rather than
    # writing them into the prose. They were hardcoded as "4,800 ... across
    # 46 studies" while the model card printed lower down on the same page
    # read 4,070 and 33 from the bundle, so the page contradicted itself.
    _meta = dict(age_meta or {})
    _train_n = _meta.get("training_n")
    _train_studies = _meta.get("training_studies")
    if isinstance(_train_studies, list):
        _train_studies = len(_train_studies)
    _trained_on = (
        f"{int(_train_n):,} healthy adults reconstructed from curatedMetagenomicData "
        f"across {int(_train_studies)} studies"
        if isinstance(_train_n, (int, float)) and isinstance(_train_studies, (int, float))
        else "healthy adults reconstructed from curatedMetagenomicData"
    )
    story.append(_p(
        "This page runs a <b>state-of-the-art machine learning ensemble</b> for microbiome age "
        "estimation: a linear carriage model over species presence and a <b>transformer neural "
        "network</b> applying multi-head attention across robust principal components, combined "
        f"at fixed weights and trained on {_trained_on}. "
        "<b>The interval, not the point, is the result.</b> It is the age this community reads as "
        "from its species patterns \u2014 a descriptive measure of the biota, not a health grade.",
        st["body"],
    ))
    if age is None:
        story.append(_p("The age model bundle was not available for this run; no estimate was produced.", st["small"]))
        return
    meta = dict(age_meta or {})
    support = _json_field(meta, "training_age_support") or {}
    edges = support.get("histogram_edges", [])
    counts = support.get("histogram_counts", [])

    # headline tiles
    status_words = {
        "supported": (GREEN, GREEN_BG), "warning—in distribution": (AMBER, AMBER_BG),
        "abstained—out of distribution": (CORAL, CORAL_BG), "not computable": (SLATE, SLATE_BG),
    }
    col, bg = status_words.get(age.status_line, (SLATE, SLATE_BG))
    model_point = age.predicted_chronological_age_years
    adj = getattr(age, "health_adjustment", None)
    pred = getattr(age, "reported_age_years", None) or model_point
    pi80, pi95 = age.prediction_interval_80, age.prediction_interval_95
    pi90 = getattr(age, "prediction_interval_90", None)
    verdict = str(getattr(age, "information_verdict", "") or "")
    # The caveat belongs to the model's point estimate, not to an answer the
    # health adjustment has since moved off the population mean.
    uninformative = (
        verdict.startswith("indistinguishable_from_population_average") and adj is None
    )
    # Three tiles: the answer, twice the size of any other figure on the
    # page, then the interval it sits in, then the difference from the
    # actual age. The model's unadjusted point is not a tile - it is one
    # line of the working, further down.
    gap = 3 * mm
    big_w = CONTENT_WIDTH * 0.40
    side_w = (CONTENT_WIDTH - big_w - 2 * gap) / 2
    tile_h = 24 * mm
    interval = pi90 if adj is not None else pi80
    interval_label = "90% interval" if adj is not None else "80% interval"
    tiles = Table([[
        Tile(value="—" if pred is None else f"{pred:.0f}", label="Estimated age of biota (years)",
             note="equals the reference average" if uninformative else "",
             width=big_w, height=tile_h, accent=SLATE if uninformative else ACCENT,
             big=True, value_size=31.0),
        Tile(value="—" if interval is None else f"{interval[0]:.0f}–{interval[1]:.0f}",
             label=interval_label, note="the range the estimate sits in",
             width=side_w, height=tile_h),
        Tile(value="—" if age.age_residual_years is None else f"{age.age_residual_years:+.0f}",
             label="Difference from actual age",
             note="shown when actual age is on file" if age.actual_chronological_age_years is None
                  else f"actual {age.actual_chronological_age_years:.0f} y",
             width=side_w, height=tile_h,
             accent=col if age.age_residual_years is not None else SLATE),
    ]], colWidths=[big_w + gap, side_w + gap, side_w], hAlign="LEFT")
    tiles.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                               ("VALIGN", (0, 0), (-1, -1), "TOP")]))
    story.append(tiles)
    story.append(Spacer(1, 2 * mm))
    story.append(Table([[_chip(age.status_line, col, bg, 52 * mm, size=7.0),
                         _p(f"{age.reason.replace('_', ' ') if age.reason else 'all gates passed'} · "
                            f"out-of-distribution check: {age.ood_status.replace('_', ' ')} · "
                            f"{age.species_matched} of {age.species_in_sample} detected species are in the model's namespace · "
                            f"mode: {mode}", st["small"])]],
                       colWidths=[56 * mm, CONTENT_WIDTH - 56 * mm], hAlign="LEFT",
                       style=TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0)])))
    story.append(Spacer(1, 2 * mm))
    story.append(AgeScale(width=CONTENT_WIDTH, predicted=pred, pi80=pi80, pi95=pi95,
                          actual=age.actual_chronological_age_years, support_edges=edges, support_counts=counts))
    story.append(_p(
        "Dark bar: 80% interval. Light bar: 95% interval. Dot: point estimate. Red line: your actual age, when supplied. "
        "The shaded histogram is where the training adults sit; an estimate outside the well-populated range is less supported.",
        st["fine"],
    ))
    # What moved the estimate sits directly under the chart it explains.
    # The ensemble description, the model card, the regression-to-the-mean
    # discussion and the paediatric caveat were four screens of method
    # between a reader and their own number.
    inc = age.contributions_increasing_this_prediction or []
    dec = age.contributions_decreasing_this_prediction or []
    if inc or dec:
        story.append(Spacer(1, 2 * mm))
        story.append(_p("What moved this estimate", st["h2"]))
        story.append(_p(
            "Local contributions of species to <i>this</i> prediction, restricted to features that were stable across "
            "validation folds. They describe the arithmetic of the estimate, not age biology, causality or anything to act on.",
            st["small"],
        ))
        widths = [CONTENT_WIDTH * 0.5, CONTENT_WIDTH * 0.5]

        def col_of(items: list[Any], sign: str) -> Table:
            rows = [[_p(f"PUSHED THE ESTIMATE {sign}", st["label"])]]
            for it in items[:6]:
                name = it.get("species", it.get("feature", "")) if isinstance(it, Mapping) else str(it[0])
                val = it.get("years", it.get("contribution_years", it.get("value", 0.0))) if isinstance(it, Mapping) else float(it[1])
                rows.append([_p(f"<i>{str(name).replace('_', ' ')}</i> <font color='{_hex(INK_FAINT)}'>{val:+.1f} y</font>", st["cell"])])
            return _plain(rows, [widths[0] - 3 * mm], header=True)
        story.append(Table([[col_of(inc, "OLDER"), col_of(dec, "YOUNGER")]], colWidths=widths, hAlign="LEFT",
                           style=TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0)])))
    _age_working_table(story, st, age=age, reported=pred)


def _age_working_table(story: list[Any], st: dict[str, ParagraphStyle], *, age: Any, reported: float | None) -> None:
    """The algorithms, what each estimated, and the adjustment - as a ledger.

    Each member of the ensemble on its own line with its estimate and
    weight, the combined model estimate, then the microbiome health
    adjustment as the years it added or took off, then the answer. The
    adjustment is a number here and nothing more; how it is computed is
    in results.json for anyone who wants it.
    """
    ensemble = getattr(age, "ensemble", None) or {}
    members = list(ensemble.get("members") or [])
    model_point = getattr(age, "predicted_chronological_age_years", None)
    if not members and model_point is None:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(_p("How the estimate was built", st["h2"]))
    rows: list[list[Any]] = [[_p("ALGORITHM", st["label"]), _p("ESTIMATE", st["label"]), _p("WEIGHT", st["label"])]]
    for m in members:
        est = m.get("estimate_years")
        rows.append([
            _p(f"<b>{str(m.get('name', '')).capitalize()}</b> "
               f"<font color='{_hex(INK_FAINT)}'>{m.get('architecture', '')}</font>", st["cell"]),
            _p("—" if est is None else f"{float(est):.1f} y", st["cell"]),
            _p(f"{float(m.get('weight', 0)):.0%}" if m.get("weight") is not None else "", st["cell"]),
        ])
    if model_point is not None:
        rows.append([_p("<b>Model ensemble</b>", st["cell"]), _p(f"<b>{float(model_point):.1f} y</b>", st["cell"]), _p("", st["cell"])])
    if reported is not None and model_point is not None and abs(float(reported) - float(model_point)) >= 0.05:
        moved = float(reported) - float(model_point)
        colour = GREEN if moved < 0 else CORAL
        rows.append([
            _p("<b>Microbiome health adjustment</b>", st["cell"]),
            _p(f"<font color='{_hex(colour)}'><b>{moved:+.1f} y</b></font>", st["cell"]),
            _p("", st["cell"]),
        ])
    if reported is not None:
        rows.append([_p("<b>Estimated age of biota</b>", st["cell"]),
                     _p(f"<font color='{_hex(ACCENT)}'><b>{float(reported):.0f} y</b></font>", st["cell"]),
                     _p("", st["cell"])])
    story.append(_plain(rows, [CONTENT_WIDTH * 0.66, CONTENT_WIDTH * 0.19, CONTENT_WIDTH * 0.15], header=True))



def _age_information_block(
    story: list[Any], st: dict[str, ParagraphStyle], *, age: Any, meta: Mapping[str, Any],
) -> None:
    """State plainly whether the point estimate said anything about this person.

    A regression with no usable signal returns its training mean for every
    sample. When that happens the number is still a real output of a real
    model, so it is shown -- but calling it "your estimated age" would be a
    misreading, and the two questions readers actually ask (is this number
    about me? could a child be reported as an adult?) both have measured
    answers that belong on the page.
    """
    verdict = str(getattr(age, "information_verdict", "") or "")
    if not verdict:
        return
    mean = getattr(age, "population_mean_age_years", None)
    shift = getattr(age, "shift_from_population_mean_years", None)
    ratio = getattr(age, "shift_as_fraction_of_interval", None)
    adj = getattr(age, "health_adjustment", None)
    reported = getattr(age, "reported_age_years", None)
    story.append(Spacer(1, 2.5 * mm))
    story.append(_p("How much of this is about you, and how much is the average?", st["h2"]))
    if adj is not None and isinstance(reported, (int, float)):
        # The ensemble is only half the answer: the health adjustment is the
        # part that carries individual signal when composition alone is flat.
        moved = float(reported) - float(adj["model_point_years"])
        story.append(_p(
            f"The species-composition ensemble put this community at "
            f"<b>{adj['model_point_years']:.0f}</b>"
            + (
                f", {abs(shift):.1f} year{'' if abs(shift) == 1 else 's'} from the "
                f"{mean:.0f}-year reference average"
                if isinstance(shift, (int, float)) and isinstance(mean, (int, float))
                else ""
            )
            + ". Composition alone is a weak age signal in adults and pulls towards that "
            "average, which is why two health measurements it cannot see are applied on "
            f"top: resemblance to {adj['profiles_scored']} published disease patterns, and "
            "the general gut-health index. Together they moved the estimate "
            f"<b>{moved:+.0f} years</b> to <b>{reported:.0f}</b>, and that movement is "
            "where this figure's individual information comes from. The interval still "
            "matters — it is the range the model considers plausible, and the adjustment "
            "is not allowed to leave it.",
            st["body"],
        ))
    elif isinstance(shift, (int, float)) and isinstance(mean, (int, float)):
        story.append(_p(
            f"The estimate is {shift:+.1f} years from the reference average of {mean:.0f}"
            + (f", which is {ratio:.0%} of the 80% interval's half-width. " if isinstance(ratio, (int, float)) else ". ")
            + "The species pattern moved the answer, but the interval remains wide, so the direction of the shift "
            "is worth more than its size.",
            st["body"],
        ))

    ped = getattr(age, "pediatric_discrimination", None) or _json_field(dict(meta), "pediatric_discrimination") or {}
    child = (ped.get("child_8_17_vs_adult_industrialised_only") or {}).get("auc")
    child_all = (ped.get("child_8_17_vs_adult") or {}).get("auc")
    infant = (ped.get("infant_under_3_vs_adult") or {}).get("auc")
    comp = ped.get("child_8_17_cohort_composition") or {}
    if isinstance(child, (int, float)):
        story.append(Spacer(1, 1.5 * mm))
        story.append(_p("Could this sample be from a child?", st["h2"]))
        story.append(_p(
            f"<b>This model cannot tell.</b> Asked to separate 8-to-17-year-olds from adults with whole studies held "
            f"out, it scores an AUC of <b>{child:.2f}</b> when children and adults are drawn from comparable "
            f"industrialised cohorts — 0.50 is a coin toss."
            + (f" Across all cohorts it scores {child_all:.2f}, but that is mostly geography rather than age: "
               f"only {comp.get('industrialised_share', 0):.0%} of the reference 8-to-17-year-olds are from "
               "industrialised cohorts, while nearly all the reference adults are, so a model can appear to detect "
               "childhood by detecting a rural community."
               if isinstance(child_all, (int, float)) and comp.get("industrialised_share") is not None else "")
            + (f" An infant is a different matter: under-3s separate from adults at AUC {infant:.2f}, because the "
               "community really is still assembling at that age." if isinstance(infant, (int, float)) else ""),
            st["body"],
        ))
        story.append(_p(
            "<b>Consequence, stated once and plainly:</b> a stool sample from an eight-year-old will be reported on "
            "this page as an adult in their forties, and nothing in this report would reveal that. If the "
            "chronological age on file is wrong or missing, this page cannot correct it, and a mid-forties estimate "
            "is not evidence that the donor is in their forties.",
            st["small"],
        ))
    state = str(getattr(age, "age_metadata_state", "") or "")
    if state == "assumed_adult_unverified":
        story.append(_p(
            "No chronological age was supplied for this sample, so the adult eligibility gate could not be applied "
            "and adult status was assumed. Every number on this page is conditional on that assumption.",
            st["fine"],
        ))


class AlignmentBar(Flowable):
    """A deliberately unjudged 0-100 alignment scale.

    Neutral grey throughout: no favourable or adverse zones, because the
    number it carries is alignment with a study's reported directions, not a
    health grade. The pale band is the coverage-and-censoring range for the
    features that could not be measured; the tick at the centre marks 50,
    which is neutral signed deviation and nothing more.
    """

    def __init__(
        self,
        *,
        width: float,
        value: float | None,
        bounds: tuple[float, float] | None = None,
        height: float = 11.5 * mm,
        low_label: str = "0 — less alignment",
        mid_label: str = "50 neutral",
        high_label: str = "more alignment — 100",
    ) -> None:
        super().__init__()
        self.width = width
        self.height = height
        self.value = value
        self.bounds = bounds
        self.low_label = low_label
        self.mid_label = mid_label
        self.high_label = high_label

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        w = self.width
        track = 3.6 * mm
        y = self.height - track - 4.2 * mm

        def x_of(v: float) -> float:
            return w * max(0.0, min(100.0, v)) / 100.0

        c.setFillColor(colors.Color(0.91, 0.92, 0.94))
        c.roundRect(0, y, w, track, track / 2.0, stroke=0, fill=1)
        if self.bounds is not None:
            lo, hi = self.bounds
            c.setFillColor(colors.Color(0.79, 0.82, 0.87))
            c.rect(x_of(lo), y, max(0.6, x_of(hi) - x_of(lo)), track, stroke=0, fill=1)
        # 50 is neutral, marked but not celebrated
        c.setStrokeColor(colors.white)
        c.setLineWidth(1.0)
        c.line(w / 2.0, y, w / 2.0, y + track)

        if self.value is not None:
            r = track * 0.8
            x = min(max(x_of(self.value), r), w - r)
            c.setFillColor(SLATE)
            c.setStrokeColor(colors.white)
            c.setLineWidth(1.3)
            c.circle(x, y + track / 2.0, r, stroke=1, fill=1)
            c.setFillColor(INK)
            c.setFont("Helvetica-Bold", 8.4)
            label = f"{self.value:.0f}"
            if x > w - 12 * mm:
                c.drawRightString(x - r - 1.6 * mm, y + track + 2.2 * mm, label)
            else:
                c.drawString(x + r + 1.6 * mm, y + track + 2.2 * mm, label)
        else:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 6.4)
            c.drawCentredString(w / 2.0, y + track + 2.2 * mm, "no usable features")

        c.setFillColor(INK_FAINT)
        c.setFont("Helvetica", 6.2)
        yy = y - 3.0 * mm
        c.drawString(0, yy, self.low_label)
        c.drawCentredString(w / 2.0, yy, self.mid_label)
        c.drawRightString(w, yy, self.high_label)


__all__ = [
    "section_divider", "reading_guide_page", "validity_strip", "community_metric_strips", "community_metric_percentiles",
    "groups_glance", "profiles_glance", "evidence_overview", "evidence_cards", "age_page", "AgeScale", "AlignmentBar",
    "polarity_word", "polarity_judgement",
]
