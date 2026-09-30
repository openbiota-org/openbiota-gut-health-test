"""The organism section: what is there, what each one is, and which ones matter.

Three parts, in reading order:

1. **Composition.** One bar. How much of the community is beneficial,
   opportunist, conditional, unknown - by share of reads, so a 13% bloom
   shows as 13% of the bar and not as one row among two hundred.
2. **The thirty largest.** Each with a plain-language description, its
   class, its abundance, and - the step a category alone does not take -
   whether the *amount* found is a concern. A reader who stops here has
   the picture.
3. **Everything.** One table, every organism detected by any catalogue,
   colour-coded by class, with flagged rows tinted so that a problem three
   pages in still stands out.

Colour carries meaning throughout and means the same thing everywhere:
green beneficial, red opportunist, amber conditional, grey unknown; a red
row tint marks an issue, an amber tint a watch.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Flowable, KeepTogether, Spacer, Table, TableStyle

from openbiota import organisms as org
from openbiota.pdflinks import Paragraph, section_heading
from openbiota.pdfreport import (
    AMBER,
    AMBER_BG,
    CONTENT_WIDTH,
    CORAL,
    CORAL_BG,
    GREEN,
    GREEN_BG,
    INK,
    INK_FAINT,
    INK_SOFT,
    PANEL_BG,
    RULE,
    SLATE,
    SLATE_BG,
)

TOP_N: Final = 30

CLASS_COLOUR: Final[dict[str, colors.Color]] = {
    org.BENEFICIAL: GREEN,
    org.OPPORTUNIST: CORAL,
    org.CONDITIONAL: AMBER,
    org.UNKNOWN: SLATE,
}
CLASS_BG: Final[dict[str, colors.Color]] = {
    org.BENEFICIAL: GREEN_BG,
    org.OPPORTUNIST: CORAL_BG,
    org.CONDITIONAL: AMBER_BG,
    org.UNKNOWN: SLATE_BG,
}
#: Sequence the primary catalogue could not place: shown, not hidden.
UNCLASSIFIED_COLOUR: Final = colors.HexColor("#D5DCE3")
#: Row tints for flagged organisms.
ISSUE_TINT: Final = colors.HexColor("#FBEAE8")
WATCH_TINT: Final = colors.HexColor("#FDF5E2")


def _hex(c: colors.Color) -> str:
    return c.hexval().replace("0x", "#")


def _pct(v: float) -> str:
    return f"{v:.3f}%" if v < 1 else f"{v:.2f}%"


#: Marks on a share whose basis is not the marker lane's own reading; the
#: footnote under the full table says what each means.
DASH: Final = "\u2014"
SPLIT_MARK: Final = "\u2020"      # dagger: genus total divided by whole-genome mapping
ESTIMATE_MARK: Final = "\u00b0"   # degree sign: estimated from the whole-genome search
ABSORBED_MARK: Final = "\u00a7"   # section sign: carries a rejected relative's marker share


def _share_cell(o: Any) -> str:
    """One share of the whole per organism, on one scale."""
    basis = getattr(o, "share_basis", "marker" if o.in_primary else "none")
    if o.in_primary and basis in ("marker", "split", "estimated", "absorbed"):
        mark = {"split": SPLIT_MARK, "estimated": ESTIMATE_MARK, "absorbed": ABSORBED_MARK}.get(basis, "")
        if basis != "absorbed" and getattr(o, "absorbed_percent", None):
            mark += ABSORBED_MARK
        return (f"<font size='7.4'><b>{_pct(o.percent)}</b></font>"
                + (f"<font size='6' color='{_hex(INK_FAINT)}'>{mark}</font>" if mark else ""))
    within = getattr(o, "counted_within", None)
    if within:
        genus = (getattr(o, "gtdb_genus", None) or o.genus or within.split(" ")[0]).split("_")[0]
        return (f"<font size='5.8' color='{_hex(INK_FAINT)}'>within</font><br/>"
                f"<font size='5.8' color='{_hex(INK_SOFT)}'><i>{genus}</i> above</font>")
    return f"<font size='6.4' color='{_hex(INK_FAINT)}'>\u2014</font>"


def _deviation_text(dev: float | None) -> str:
    """``+420%`` / ``-90%`` against the typical carrier; a tenth of a percent is noise at this scale."""
    if dev is None:
        return ""
    if abs(dev) < 5:
        return "\u2248 typical"
    if dev >= 900:
        # "+10,883%" wraps and is hard to read; ten times and above is said as a multiple
        return f"\u00d7{(dev / 100.0 + 1.0):,.0f}"
    sign = "+" if dev > 0 else "\u2212"
    return f"{sign}{abs(dev):,.0f}%"


def _level_cell(v: org.Verdict, *, with_percentile: bool) -> str:
    """How far from the reference, in words a reader can picture.

    The deviation is the sample's reading against the median level among
    reference adults who carry the organism, on the lane the reference was
    measured with; the percentile says where that sits in the spread. The
    flag (high / low) is the class-aware judgement from the percentile.
    """
    o = v.organism
    dev = getattr(o, "deviation_percent", None)
    if getattr(o, "reference_conflict", False):
        return (f"<font size='6' color='{_hex(INK_FAINT)}'>reference catalogue reads<br/>this species differently</font>")
    if getattr(o, "few_reference_carriers", False):
        n = getattr(o, "reference_carriers", 0)
        return (f"<font size='6' color='{_hex(INK_FAINT)}'>only {n} reference adult{'s' if n != 1 else ''} "
                "carry it:<br/>too few to rank against</font>")
    pct = getattr(o, "level_percentile", o.percentile)
    if pct is None and dev is None:
        return f"<font size='6.4' color='{_hex(INK_FAINT)}'>no reference</font>"
    col = CORAL if v.is_issue else (AMBER if v.flagged else INK)
    parts = []
    if v.flagged:
        arrow = {"high": "\u25b2", "low": "\u25bc"}.get(v.flag, "\u25c7")
        word = v.flag_reason.split(" \u2014 ")[0]
        parts.append(f"<font color='{_hex(col)}' size='7'><b>{arrow} {word}</b></font>")
    dev_text = _deviation_text(dev)
    if dev_text:
        typical = getattr(o, "reference_percent", None)
        parts.append(f"<font color='{_hex(col)}' size='7.4'><b>{dev_text}</b></font>"
                     f"<font size='5.4' color='{_hex(INK_FAINT)}'> vs typical carrier"
                     + (f" ({_pct(typical)})" if typical else "") + "</font>")
    if with_percentile and pct is not None:
        src = getattr(o, "percentile_source", None) or "scoring cohort"
        parts.append(f"<font size='5.8' color='{_hex(INK_SOFT)}'>{org._ordinal(pct)} percentile"
                     + ("" if src == "scoring cohort" else "\u2021") + "</font>")
    return "<br/>".join(parts) if parts else f"<font size='6.4' color='{_hex(INK_FAINT)}'>\u2014</font>"


# --------------------------------------------------------------------------- #
# composition bar
# --------------------------------------------------------------------------- #


class CompositionBar(Flowable):
    """One stacked bar: the community by class, as a share of reads."""

    def __init__(self, shares: dict[str, float], *, width: float, height: float = 9 * mm) -> None:
        super().__init__()
        self.shares = shares
        self.width = width
        self.height = height

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        order = (org.BENEFICIAL, org.CONDITIONAL, org.OPPORTUNIST, org.UNKNOWN, "unclassified")
        total = sum(self.shares.get(k, 0.0) for k in order) or 1.0
        x = 0.0
        r = 1.8 * mm
        # Clip to a rounded rectangle so the segments share the bar's corners.
        p = c.beginPath()
        p.roundRect(0, 0, self.width, self.height, r)
        c.saveState()
        c.clipPath(p, stroke=0, fill=0)
        for k in order:
            w = self.width * self.shares.get(k, 0.0) / total
            if w <= 0:
                continue
            c.setFillColor(CLASS_COLOUR.get(k, UNCLASSIFIED_COLOUR))
            c.rect(x, 0, w, self.height, stroke=0, fill=1)
            # A label inside the segment where there is room for it.
            share = self.shares.get(k, 0.0)
            label = f"{share:.0f}%"
            c.setFont("Helvetica-Bold", 8)
            if w > c.stringWidth(label, "Helvetica-Bold", 8) + 3 * mm:
                c.setFillColor(colors.white)
                c.drawCentredString(x + w / 2, self.height / 2 - 2.8, label)
            x += w
        c.restoreState()
        c.setStrokeColor(RULE)
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, r, stroke=1, fill=0)


def _legend(st: dict[str, ParagraphStyle], shares: dict[str, float], vs: Sequence[org.Verdict]) -> Table:
    """Four cells: class, share, count, and what the class means."""
    cells = []
    for k in (org.BENEFICIAL, org.CONDITIONAL, org.OPPORTUNIST, org.UNKNOWN):
        n = sum(1 for v in vs if v.cls == k)
        flagged = sum(1 for v in vs if v.cls == k and v.flagged)
        head = (
            f"<font color='{_hex(CLASS_COLOUR[k])}' size='9'>\u25cf</font> "
            f"<font size='8.4'><b>{org.CLASS_LABEL[k]}</b></font>"
            f"<font size='8.4' color='{_hex(INK_SOFT)}'>&nbsp;&nbsp;{shares.get(k, 0):.0f}%</font>"
        )
        sub = f"{n} organism{'s' if n != 1 else ''}"
        if flagged:
            sub += (f" \u00b7 <font color='{_hex(CORAL if k == org.OPPORTUNIST else AMBER)}'>"
                    f"{flagged} flagged</font>")
        cells.append([
            Paragraph(head, ParagraphStyle("lgh", parent=st["body"], fontSize=8.4, leading=11)),
            Paragraph(f"<font size='6.6' color='{_hex(INK_FAINT)}'>{sub}</font>",
                      ParagraphStyle("lgs", parent=st["small"], fontSize=6.6, leading=8.6)),
            Paragraph(f"<font size='6.8' color='{_hex(INK_SOFT)}'>{org.CLASS_MEANING[k]}</font>",
                      ParagraphStyle("lgm", parent=st["small"], fontSize=6.8, leading=9.2)),
        ])
    if shares.get("unclassified"):
        # The fifth segment is not an organism class: it is sequence no
        # catalogue could place, drawn so the bar adds to the whole. It must
        # be in the legend or two greys read as one unexplained number.
        cells.append([
            Paragraph(
                f"<font color='{_hex(UNCLASSIFIED_COLOUR)}' size='9'>\u25cf</font> "
                f"<font size='8.4'><b>Unclassified</b></font>"
                f"<font size='8.4' color='{_hex(INK_SOFT)}'>&nbsp;&nbsp;{shares['unclassified']:.0f}%</font>",
                ParagraphStyle("lgh", parent=st["body"], fontSize=8.4, leading=11)),
            Paragraph(f"<font size='6.6' color='{_hex(INK_FAINT)}'>sequence, not organisms</font>",
                      ParagraphStyle("lgs", parent=st["small"], fontSize=6.6, leading=8.6)),
            Paragraph(f"<font size='6.8' color='{_hex(INK_SOFT)}'>DNA that matched nothing in any catalogue. Not counted in the four "
                      "classes; shown so the bar adds to 100%.</font>",
                      ParagraphStyle("lgm", parent=st["small"], fontSize=6.8, leading=9.2)),
        ])
    n = len(cells)
    gap = 4 * mm
    w = (CONTENT_WIDTH - (n - 1) * gap) / n
    t = Table([cells], colWidths=[w + gap] * (n - 1) + [w], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-2, -1), gap),
        ("RIGHTPADDING", (-1, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    return t


# --------------------------------------------------------------------------- #
# shared cells
# --------------------------------------------------------------------------- #


def _class_chip(v: org.Verdict) -> str:
    return (f"<font color='{_hex(CLASS_COLOUR[v.cls])}' size='8'>\u25cf</font> "
            f"<font size='6.8'>{v.label}</font>")


def _name_cell(v: org.Verdict, *, size: float = 7.4) -> str:
    o = v.organism
    name = f"<i>{o.display}</i>" if not o.provisional_name else o.display
    out = f"<font size='{size}'><b>{name}</b></font>"
    if o.formerly:
        out += (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>formerly "
                f"{o.formerly.replace('_', ' ')}</font>")
    elif o.gtdb and o.gtdb.replace("_", " ") != o.display.replace("_", " ") and not o.provisional_name:
        out += f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>also {o.gtdb}</font>"
    if o.provisional_name:
        out += (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>not yet formally named "
                "\u00b7 known from genomes</font>")
    listed = getattr(o, "formerly_listed_as", ()) or ()
    if listed:
        out += (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>an earlier catalogue listed this population as "
                f"{', '.join(x.replace('_', ' ') for x in listed[:2])}</font>")
    absorbed = getattr(o, "absorbed_from", ()) or ()
    if absorbed:
        out += (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>includes reads the marker catalogue read as "
                f"{', '.join(x.replace('_', ' ') for x in absorbed[:2])}</font>")
    # The lanes' own identifiers and how many methods saw it (spec 0.8.4 §6).
    chips = []
    if getattr(o, "status", "supported") == "ambiguous":
        chips.append(f"<font color='{_hex(INK_SOFT)}' size='5.6'>sibling species, one population</font>")
    from openbiota.pdflibrary import _identifiers_line

    ids = _identifiers_line(o)
    if chips or ids:
        out += "<br/>" + " ".join(chips) + (" " if chips and ids else "") + (
            f"<font size='5.2' color='{_hex(INK_FAINT)}'>{ids}</font>" if ids else "")
    return out


def _row_tint(v: org.Verdict) -> colors.Color | None:
    if v.is_issue:
        return ISSUE_TINT
    if v.flagged:
        return WATCH_TINT
    return None


# --------------------------------------------------------------------------- #
# the section
# --------------------------------------------------------------------------- #


def organism_section(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    vs: Sequence[org.Verdict],
    cohort_note: str,
    n_strain_resolved: int = 0,
    unclassified_percent: float | None = None,
    inv: Any = None,
) -> None:
    """Composition, the thirty largest with descriptions, then everything."""
    vs = list(vs)
    if not vs:
        story.extend(section_heading(section, "Every organism detected", st["h1"], opens_section=True))
        story.append(Paragraph("No taxonomic profile was produced for this sample.", st["body"]))
        return

    shares = org.composition(vs, unclassified=unclassified_percent)
    flagged = org.issues(vs)
    n_estimated = sum(1 for v in vs if getattr(v.organism, "share_basis", "") == "estimated")
    est_total = sum(v.organism.percent for v in vs if getattr(v.organism, "share_basis", "") == "estimated")
    n_absorbed = sum(len(getattr(v.organism, "absorbed_from", ()) or ()) for v in vs)
    n_member = sum(1 for v in vs if getattr(v.organism, "share_basis", "") == "member")
    n_issue = sum(1 for v in vs if v.is_issue)
    n_watch = sum(1 for v in vs if v.flagged and not v.is_issue)

    # ---- 1. composition --------------------------------------------------- #
    story.extend(section_heading(section, "Your organisms, classified", st["h1"], sub=1,
                                 opens_section=True))
    story.append(Paragraph(
        f"<b>{len(vs)} organisms</b> were found in this sample, pooling every reference "
        "catalogue available. Each is placed in one of four classes from the human literature, "
        "and \u2014 where a reference population exists for it \u2014 its level is compared with "
        f"{cohort_note}. The bar is the community by share of reads: a dominant organism is a "
        "wide band, a trace one is invisible here and listed below."
        + (f" The grey tail is the {shares.get('unclassified', 0):.0f}% of sequence that matched "
           "nothing on record \u2014 shown so the parts add to the whole."
           if shares.get("unclassified") else ""),
        st["body"]))
    _counts_paragraph(story, st, vs, inv, n_estimated=n_estimated, est_total=est_total, n_member=n_member,
                      n_absorbed=n_absorbed)
    story.append(Spacer(1, 2.5 * mm))
    story.append(CompositionBar(shares, width=CONTENT_WIDTH))
    story.append(Spacer(1, 2.2 * mm))
    story.append(_legend(st, shares, vs))
    story.append(Spacer(1, 3 * mm))

    # The verdict line. This is what the section is for.
    if n_issue or n_watch:
        parts = []
        if n_issue:
            parts.append(f"<font color='{_hex(CORAL)}'><b>{n_issue} organism"
                         f"{'s' if n_issue != 1 else ''} at a level that is a concern</b></font>")
        if n_watch:
            parts.append(f"<font color='{_hex(AMBER)}'><b>{n_watch} worth watching</b></font>")
        lead = flagged[0]
        dev = _deviation_text(getattr(lead.organism, "deviation_percent", None))
        story.append(Paragraph(
            " and ".join(parts) + ". The most consequential is "
            f"<b><i>{lead.organism.display}</i></b> \u2014 {lead.label.lower()}, "
            f"{_pct(lead.organism.percent)} of the community, {lead.flag_reason}"
            + (f" ({dev} against the typical carrier)" if dev else "") + ". "
            "Flagged rows are tinted throughout this section.",
            st["body"]))
    else:
        story.append(Paragraph(
            f"<font color='{_hex(GREEN)}'><b>No organism is at a level that raises a flag.</b>"
            "</font> Every organism with a reference range sits inside it, in the direction that "
            "matters for its class.",
            st["body"]))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_FAINT)}'>Class is a property of the organism, from "
        "the literature; the flag is a property of this sample, from its percentile. "
        "<b>vs typical carrier</b> is how far the level sits from the median of reference adults who "
        "carry the organism at all (that median is the figure in brackets): +400% is five times their level, "
        "\u221290% a tenth of it, \u00d7110 a hundred and ten times. Where fewer than ten reference adults "
        "carry an organism, no rank is stated. An organism "
        "with no reference was measured on a catalogue the reference population was not, so there is "
        "nothing to compare it against; the detection stands on its own.</font>",
        st["small"]))
    _deviation_tables(story, st, vs, cohort_note)

    # ---- 2. the thirty largest ------------------------------------------ #
    story.append(Spacer(1, 5 * mm))
    top = sorted((v for v in vs if v.organism.in_primary), key=lambda v: -v.organism.percent)[:TOP_N]
    story.extend(section_heading(section, f"The {len(top)} most abundant, one by one",
                                 st["h1"], sub=2))
    total = sum(v.organism.percent for v in vs if v.organism.in_primary) or 1.0
    top_share = 100.0 * sum(v.organism.percent for v in top if v.organism.in_primary) / total
    story.append(Paragraph(
        f"Together these make up <b>{min(top_share, 100):.0f}%</b> of the composition. Each has "
        "what it is, what the research says about it, its class, and whether the amount found is "
        "a concern.",
        st["body"]))
    story.append(Spacer(1, 2 * mm))
    _top_table(story, st, top)

    # ---- 3. everything ---------------------------------------------------- #
    story.append(Spacer(1, 5 * mm))
    story.extend(section_heading(section, "Every organism detected", st["h1"], sub=3))
    story.append(Paragraph(
        f"All {len(vs)} organisms, largest share first, colour-coded by class. <b>Level</b> is how far "
        f"the organism sits from the typical carrier in {cohort_note}, with its percentile; "
        "<b>carried by</b> is the share of that group in which the organism is found at all. Rows tinted "
        "red are at a level that is a concern; amber, worth watching. "
        + (f"<b>{n_strain_resolved}</b> were also read to strain level, marked in the last column."
           if n_strain_resolved else ""),
        st["body"]))
    story.append(Spacer(1, 2 * mm))
    _full_table(story, st, vs, strained=n_strain_resolved > 0)


def _counts_paragraph(story: list[Any], st: dict[str, ParagraphStyle], vs: Sequence[org.Verdict], inv: Any,
                      *, n_estimated: int = 0, est_total: float = 0.0, n_member: int = 0,
                      n_absorbed: int = 0) -> None:
    """What the count is made of (spec 0.8.4 §6: named species, unnamed
    clusters and complexes counted separately, never summed into a claim)."""
    orgs = [v.organism for v in vs]
    cat = {k: sum(1 for o in orgs if getattr(o, "count_category", "named_species") == k)
           for k in ("named_species", "unnamed_species_cluster", "unresolved_complex", "higher_rank")}
    n_lanes = len(getattr(inv, "lanes", ()) or ()) if inv is not None else 0
    n_single = sum(1 for o in orgs if len(getattr(o, "methods", ()) or ()) <= 1)
    text = (
        (f"Detection pools <b>{n_lanes} methods</b>, each with its own reference catalogue, so a species one "
         "catalogue has no genome for is still found when another can name it. " if n_lanes else "")
        + f"Counted separately: <b>{cat['named_species']}</b> named species and <b>{cat['unnamed_species_cluster']}</b> "
        "unnamed species-level clusters \u2014 organisms known only from assembled genomes, real and counted, whose "
        "names will change when they are formally described"
        + (f"; <b>{cat['unresolved_complex']}</b> populations that the catalogues place in sibling species and that "
           "are listed once" if cat["unresolved_complex"] else "")
        + (f"; <b>{cat['higher_rank']}</b> placed above species" if cat["higher_rank"] else "")
        + ". Each organism has one share of the whole, on the scale of the marker catalogue that measures the "
        "composition"
        + (f"; for <b>{n_estimated}</b> organisms that catalogue has no entry, so their shares "
           f"({est_total:.1f}% together, marked {ESTIMATE_MARK}) are estimated from the whole-genome search and "
           "drawn from the unclassified band" if n_estimated else "")
        + (f"; <b>{n_member}</b> populations the whole-genome search tells apart within a relative's share are listed "
           "with that relative and carry no share of their own" if n_member else "")
        + (f"; for <b>{n_absorbed}</b> call{'s' if n_absorbed != 1 else ''} the whole-genome competition rejected, "
           f"the marker share went to the relative the reads belong to (marked {ABSORBED_MARK})" if n_absorbed else "")
        + ". "
        + (f"<b>{n_single}</b> organisms rest on a single detection method, said so beside their identifiers."
           if n_single else "")
    )
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(text, st["body"]))


def _deviation_tables(story: list[Any], st: dict[str, ParagraphStyle], vs: Sequence[org.Verdict], cohort_note: str) -> None:
    """The organisms furthest from the reference, in each direction.

    The full table is sorted by share; this is the same data sorted by how
    far each organism sits from the typical carrier, so the largest excesses
    and the deepest shortfalls can be read without scanning hundreds of
    rows. Above: organisms at 0.1% or more (a five-fold excess of a trace
    organism is noise). Below: organisms most reference adults carry (a
    shortfall in something most people lack means nothing).
    """
    def dev(v: org.Verdict) -> float | None:
        return getattr(v.organism, "deviation_percent", None)

    placed = [v for v in vs if v.organism.in_primary and v.organism.percent > 0]
    above = sorted((v for v in placed if (dev(v) or 0) >= 50 and v.organism.percent >= 0.1),
                   key=lambda v: -(dev(v) or 0))[:8]
    below = sorted((v for v in placed if (dev(v) or 0) <= -50 and (v.organism.prevalence or 0) >= 0.5),
                   key=lambda v: dev(v) or 0)[:8]
    if not above and not below:
        return
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Furthest from the reference", st["h3"]))
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_SOFT)}'>Against the typical carrier in {cohort_note}. "
        "Left: the largest excesses among organisms at 0.1% or more. Right: the deepest shortfalls among "
        "organisms carried by at least half of the reference adults.</font>", st["small"]))
    story.append(Spacer(1, 1.5 * mm))

    def column(title: str, items: Sequence[org.Verdict]) -> Table:
        rows: list[list[Any]] = [[Paragraph(title, st["label"]), Paragraph("SHARE", st["label"]),
                                  Paragraph("VS TYPICAL", st["label"])]]
        for v in items:
            o = v.organism
            col = CORAL if v.is_issue else (AMBER if v.flagged else INK)
            dot = f"<font color='{_hex(CLASS_COLOUR[v.cls])}' size='7'>\u25cf</font> "
            level = f"<font size='7' color='{_hex(col)}'><b>{_deviation_text(dev(v))}</b></font>"
            lp = getattr(o, "level_percentile", o.percentile)
            if lp is not None:
                level += f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>{org._ordinal(lp)} pct</font>"
            rows.append([
                Paragraph(dot + f"<font size='6.8'><i>{o.display}</i></font>", st["cell"]),
                Paragraph(f"<font size='6.8'>{_pct(o.percent) if o.in_primary else DASH}</font>", st["cell"]),
                Paragraph(level, st["cell"]),
            ])
        if len(rows) == 1:
            rows.append([Paragraph(f"<font size='6.6' color='{_hex(INK_FAINT)}'>none</font>", st["cell"]), "", ""])
        half = (CONTENT_WIDTH - 4 * mm) / 2
        t = Table(rows, colWidths=[half * 0.56, half * 0.18, half * 0.26], hAlign="LEFT")
        t.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE),
            ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2), ("TOPPADDING", (0, 0), (-1, -1), 1.6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.8),
        ]))
        return t

    outer = Table([[column("MOST ABOVE THE REFERENCE", above), column("MOST BELOW THE REFERENCE", below)]],
                  colWidths=[CONTENT_WIDTH / 2, CONTENT_WIDTH / 2], hAlign="LEFT")
    outer.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                               ("RIGHTPADDING", (0, 0), (0, 0), 4 * mm), ("RIGHTPADDING", (1, 0), (1, 0), 0),
                               ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]))
    story.append(outer)


def _top_table(story: list[Any], st: dict[str, ParagraphStyle], top: Sequence[org.Verdict]) -> None:
    widths = [CONTENT_WIDTH * w for w in (0.05, 0.22, 0.45, 0.12, 0.16)]
    label = st["label"]
    header = [Paragraph(h, label) for h in ("#", "ORGANISM", "WHAT IT IS", "CLASS · SHARE", "LEVEL")]
    body = ParagraphStyle("topd", parent=st["small"], fontSize=6.9, leading=9.2, textColor=INK_SOFT)
    rows: list[list[Any]] = [header]
    tints: list[tuple[int, colors.Color]] = []
    for i, v in enumerate(top, start=1):
        basis = {"species": "", "guild": "", "genus": " (genus-level)", "family": " (family-level)",
                 "none": ""}[v.basis]
        desc = v.description + (f"<font color='{_hex(INK_FAINT)}'>{basis}</font>" if basis else "")
        rows.append([
            Paragraph(f"<font size='7.4' color='{_hex(INK_FAINT)}'>{i}</font>", st["cell"]),
            Paragraph(_name_cell(v), st["cell"]),
            Paragraph(desc, body),
            Paragraph(_class_chip(v) + "<br/>" + _share_cell(v.organism), st["cell"]),
            Paragraph(_level_cell(v, with_percentile=True), st["cell"]),
        ])
        tint = _row_tint(v)
        if tint is not None:
            tints.append((i, tint))
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style: list[Any] = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]
    for r, tint in tints:
        style.append(("BACKGROUND", (0, r), (-1, r), tint))
    t.setStyle(TableStyle(style))
    story.append(t)


def _strain_cell(fp: Mapping[str, Any]) -> str:
    """Marker typing and comparative placement, whichever the organism has."""
    if not fp:
        return f"<font size='6.2' color='{_hex(INK_FAINT)}'>\u2014</font>"
    bits: list[str] = []
    if fp.get("markers_resolved"):
        bits.append(f"{fp['markers_resolved']} markers")
    elif fp.get("n_markers"):
        bits.append(f"{fp['n_markers']} markers")
    near = str(fp.get("nearest_reference") or "").split("|")[0].replace(".fna.gz", "").replace(".fna", "")
    if near:
        dist = fp.get("distance_to_nearest")
        bits.append(f"nearest {near}" + (f" ({float(dist):.3f})" if isinstance(dist, (int, float)) else ""))
    return (f"<font color='{_hex(GREEN)}'>\u25cf</font> <font size='6.2' color='{_hex(INK_SOFT)}'>"
            + "; ".join(bits) + "</font>")


def _full_table(story: list[Any], st: dict[str, ParagraphStyle], vs: Sequence[org.Verdict],
                *, strained: bool) -> None:
    if strained:
        widths = [CONTENT_WIDTH * w for w in (0.28, 0.13, 0.10, 0.15, 0.09, 0.14, 0.11)]
        heads = ("ORGANISM", "CLASS", "SHARE", "LEVEL", "CARRIED BY", "PERCENTILE", "STRAIN")
    else:
        widths = [CONTENT_WIDTH * w for w in (0.32, 0.14, 0.11, 0.17, 0.10, 0.16)]
        heads = ("ORGANISM", "CLASS", "SHARE", "LEVEL", "CARRIED BY", "PERCENTILE")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in heads]]
    tints: list[tuple[int, colors.Color]] = []
    ordered = _ordered(vs)
    for i, v in enumerate(ordered, start=1):
        o = v.organism
        src = getattr(o, "percentile_source", None) or "scoring cohort"
        if getattr(o, "reference_conflict", False):
            pct_text = f"<font size='6.4' color='{_hex(INK_FAINT)}'>not comparable</font>"
        else:
            lp = getattr(o, "level_percentile", o.percentile)
            pct_text = (f"<font size='6.4' color='{_hex(INK_FAINT)}'>not in reference set</font>"
                        if lp is None else f"{org._ordinal(lp)}" + ("" if src == "scoring cohort" else "\u2021"))
        prev = "\u2014" if o.prevalence is None else f"{o.prevalence:.0%}"
        row: list[Any] = [
            Paragraph(_name_cell(v, size=7.0), st["cell"]),
            Paragraph(_class_chip(v), st["cell"]),
            Paragraph(_share_cell(o), st["cell"]),
            Paragraph(_level_cell(v, with_percentile=False), st["cell"]),
            Paragraph(prev, st["cell"]),
            Paragraph(pct_text, st["cell"]),
        ]
        if strained:
            row.append(Paragraph(_strain_cell(o.strain or {}), st["fine"]))
        rows.append(row)
        tint = _row_tint(v)
        if tint is not None:
            tints.append((i, tint))
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style: list[Any] = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 1.8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PANEL_BG]),
    ]
    for r, tint in tints:
        style.append(("BACKGROUND", (0, r), (-1, r), tint))
    t.setStyle(TableStyle(style))
    story.append(t)
    story.append(Spacer(1, 2 * mm))
    n_split = sum(1 for v in vs if getattr(v.organism, "share_basis", "") == "split")
    n_est = sum(1 for v in vs if getattr(v.organism, "share_basis", "") == "estimated")
    n_member = sum(1 for v in vs if getattr(v.organism, "share_basis", "") == "member")
    notes = ["<b>Share</b> is the organism's percentage of the whole community, on one scale."]
    if n_split:
        notes.append(f"{SPLIT_MARK} the marker catalogue measured these populations as one species; its total is "
                     "divided among them in proportion to the reads that mapped uniquely to each genome.")
    if n_est:
        notes.append(f"{ESTIMATE_MARK} the marker catalogue has no entry for this organism; the share is its "
                     "whole-genome reading brought onto the same scale, and comes out of the unclassified band.")
    if any(getattr(v.organism, "absorbed_percent", None) for v in vs):
        notes.append(f"{ABSORBED_MARK} includes a share the marker catalogue had read under a relative; the "
                     "whole-genome competition showed those reads belong to this organism, so the share came "
                     "with them (the relative's name is given beneath). Where that makes the share differ "
                     "fivefold from what the reference catalogue read under this name, the percentile is left out.")
    if n_member:
        notes.append("<i>within</i>: a population the whole-genome search tells apart inside its genus, whose "
                     "reads the marker catalogue counted under the relatives listed above it; it carries no share of "
                     "its own.")
    if any(getattr(v.organism, "reference_conflict", False) for v in vs):
        notes.append("<i>not comparable</i>: the reference catalogue's reading of this species differs more than "
                     "fivefold from the share shown - it draws the species' boundary or its markers differently - so "
                     "its percentile is not a statement about the share and is left out.")
    notes.append("Where an organism has been reclassified, the current name is shown with the previous one "
                 "beneath it; <i>also</i> gives its name in the Genome Taxonomy Database. A name made of a genus "
                 "and a genome identifier marks an organism known only from assembled genomes.")
    story.append(KeepTogether([Paragraph(
        f"<font size='6.4' color='{_hex(INK_FAINT)}'>" + " ".join(notes) + "</font>", st["small"])]))


def _ordered(vs: Sequence[org.Verdict]) -> list[org.Verdict]:
    """Largest share first; a population counted within a relative's share follows that relative."""
    by_display = {v.organism.display: v for v in vs}

    def anchor_percent(v: org.Verdict) -> float:
        o = v.organism
        if o.in_primary:
            return o.percent
        within = getattr(o, "counted_within", None)
        a = by_display.get(within) if within else None
        return a.organism.percent if a is not None else 0.0

    return sorted(vs, key=lambda v: (-anchor_percent(v), 0 if v.organism.in_primary else 1,
                                     -(v.organism.secondary_percent or 0.0), v.organism.display))


__all__ = ["CompositionBar", "TOP_N", "organism_section"]
