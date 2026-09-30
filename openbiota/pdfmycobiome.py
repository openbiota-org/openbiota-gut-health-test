"""Your gut mycobiome: the report pages for the fungal module (spec §3).

Part A is one overview page: the MHS-E1 gauge first, then the fungal
signal, the fungi found, their strains, the within-fungi composition,
beneficial and concerning evidence, context and next steps. Part B has a
card per fungus and the full score trace. Blue is measured composition,
grey is unknown or unassessed, amber is a finding worth reviewing; a green
evidence label means a beneficial property was shown under stated
conditions, not that the organism is good.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, Flowable, KeepTogether, Spacer, Table, TableStyle

from openbiota.pdflinks import Anchor, Paragraph, back_link_line, section_dest, section_heading
from openbiota.pdfreport import (
    ACCENT,
    AMBER,
    AMBER_BG,
    CONTENT_WIDTH,
    CORAL,
    GREEN,
    INK,
    INK_FAINT,
    INK_SOFT,
    PANEL_BG,
    RULE,
    SLATE,
    SLATE_BG,
    Status,
    StatusChip,
    Tile,
)

BLUE = colors.HexColor("#2F6FB7")
BLUE_BG = colors.HexColor("#E3ECF7")
TEAL = ACCENT
GREY = colors.HexColor("#B7C1CA")
DASH = "\u2014"


def _hex(c: colors.Color) -> str:
    return c.hexval().replace("0x", "#")


def _lerp(a: colors.Color, b: colors.Color, t: float) -> colors.Color:
    return colors.Color(a.red + (b.red - a.red) * t, a.green + (b.green - a.green) * t, a.blue + (b.blue - a.blue) * t)


class MhsGauge(Flowable):
    """The 0-100 MHS-E1 gauge: amber toward concerning evidence, neutral grey
    around 50, teal toward favorable evidence. The point is a marker with
    the integer above it; the displayed sensitivity range is a darker band
    under the track. Endpoint words are drawn, so the meaning survives
    greyscale; the score is never conveyed by colour alone."""

    def __init__(self, *, width: float, score: float | None, low: float | None, high: float | None,
                 height: float = 7.5 * mm, compact: bool = False) -> None:
        super().__init__()
        self.width = width
        self.compact = compact
        self.height = height + (1.5 * mm if compact else 9.5 * mm)
        self.score, self.low, self.high = score, low, high
        self.track_h = height

    def wrap(self, _aw: float, _ah: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c = self.canv
        w, th = self.width, self.track_h
        y0 = 0.75 * mm if self.compact else 5.0 * mm
        n = 60
        pale_amber = _lerp(AMBER, colors.white, 0.55)
        pale_teal = _lerp(TEAL, colors.white, 0.55)
        pale_grey = _lerp(GREY, colors.white, 0.45)
        c.saveState()
        p = c.beginPath()
        p.roundRect(0, y0, w, th, th / 2)
        c.clipPath(p, stroke=0, fill=0)
        for i in range(n):
            t = i / (n - 1)
            col = _lerp(pale_amber, pale_grey, t / 0.5) if t < 0.5 else _lerp(pale_grey, pale_teal, (t - 0.5) / 0.5)
            c.setFillColor(col)
            c.rect(w * i / n - 0.5, y0, w / n + 1, th, stroke=0, fill=1)
        c.restoreState()
        # Midpoint tick.
        c.setStrokeColor(colors.white)
        c.setLineWidth(1.2)
        c.line(w / 2, y0 + 1, w / 2, y0 + th - 1)
        if self.score is not None:
            x = w * max(0.0, min(100.0, self.score)) / 100.0
            if self.low is not None and self.high is not None:
                xl, xh = w * self.low / 100.0, w * self.high / 100.0
                c.setFillColor(colors.Color(0.2, 0.25, 0.3, alpha=0.18))
                c.rect(xl, y0 + th * 0.28, max(1.0, xh - xl), th * 0.44, stroke=0, fill=1)
                c.setStrokeColor(INK_SOFT)
                c.setLineWidth(0.8)
                for xx in (xl, xh):
                    c.line(xx, y0 + th * 0.2, xx, y0 + th * 0.8)
            marker = INK
            c.setFillColor(colors.white)
            c.setStrokeColor(marker)
            c.setLineWidth(1.4)
            c.circle(x, y0 + th / 2, th * 0.36, stroke=1, fill=1)
            c.setFillColor(marker)
            c.circle(x, y0 + th / 2, th * 0.17, stroke=0, fill=1)
            if not self.compact:
                c.setFont("Helvetica-Bold", 12)
                c.drawCentredString(x, y0 + th + 2.2 * mm, f"{int(self.score + 0.5)}")
        if self.compact:
            return
        c.setFont("Helvetica", 6.0)
        c.setFillColor(INK_SOFT)
        c.drawString(0, 1.2 * mm, "More concerning evidence")
        c.drawCentredString(w / 2, 1.2 * mm, "Neutral / mixed evidence")
        c.drawRightString(w, 1.2 * mm, "More favorable evidence")


def _pct(v: float | None, digits: int = 4) -> str:
    if v is None:
        return "\u2014"
    return f"{100 * v:.{digits}f}%"


def _state_status(state: str) -> Status:
    return {
        "supported": Status("supported", BLUE, BLUE_BG, "", 0),
        "provisional": Status("provisional", SLATE, SLATE_BG, "", 0),
        "ambiguous_complex": Status("genus only", AMBER, AMBER_BG, "", 0),
        "trace": Status("trace", INK_FAINT, PANEL_BG, "", 0),
        "not_supported": Status("not supported", INK_FAINT, PANEL_BG, "", 0),
    }.get(state, Status(state.replace("_", " "), SLATE, SLATE_BG, "", 0))


def _evidence_dots(t: Mapping[str, Any]) -> str:
    fav = sum(1 for c in t.get("evidence") or [] if c.get("direction") == "favorable")
    con = sum(1 for c in t.get("evidence") or [] if c.get("direction") == "concerning")
    parts = []
    if fav:
        parts.append(f"<font color='{_hex(GREEN)}'>\u25cf</font>" * fav)
    if con:
        parts.append(f"<font color='{_hex(AMBER)}'>\u25cf</font>" * con)
    return " ".join(parts) if parts else f"<font color='{_hex(INK_FAINT)}'>\u2014</font>"


def _searched_table(st, myco: Mapping[str, Any]) -> Table:
    """How many species were searched for, by the groups a reader asks
    about, and what was found in each."""
    sr = myco.get("searched") or {}
    taxa = [t for t in myco.get("taxa") or [] if t["detection_state"] in ("supported", "ambiguous_complex", "provisional", "trace")]
    from openbiota.mycobiome.engine import SEARCH_GROUPS
    genera_of = {label: set(g) for label, g in SEARCH_GROUPS}
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("GROUP", "SPECIES SEARCHED", "FOUND IN YOUR SAMPLE")]]
    for g in sr.get("groups") or []:
        genera = genera_of.get(g["group"], set())
        found = [t for t in taxa if (t.get("genus") or t["accepted_name"].split(" ")[0]) in genera] if genera else \
                [t for t in taxa if not any((t.get("genus") or t["accepted_name"].split(" ")[0]) in gs for gs in genera_of.values())]
        if found:
            parts = []
            for t in found:
                frag = (t.get("genome_lane") or {}).get("fragments") or (t.get("pathogen_screen") or {}).get("fragments")
                state = {"supported": "", "provisional": " (provisional)", "trace": " (trace)", "ambiguous_complex": " (genus)"}[t["detection_state"]]
                parts.append(f"<i>{t['accepted_name']}</i>{state}" + (f" {frag}" if frag else ""))
            txt = f"<font color='{_hex(BLUE)}'><b>{len(found)}</b></font> \u00b7 " + ", ".join(parts)
        else:
            txt = f"<font color='{_hex(INK_FAINT)}'>none</font>"
        rows.append([Paragraph(f"<b>{g['group']}</b>", st["cell"]),
                     Paragraph(f"{g['species']:,} <font size='6' color='{_hex(INK_FAINT)}'>({g['assemblies']:,} genomes)</font>", st["cell"]),
                     Paragraph(f"<font size='6.8'>{txt}</font>", st["cell"])])
    t = Table(rows, colWidths=[CONTENT_WIDTH * w for w in (0.30, 0.18, 0.52)], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                           ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                           ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2)]))
    return t


def _healthy_context(st, myco: Mapping[str, Any]) -> list[Any]:
    """What healthy adults typically carry, and where this sample sits."""
    ctx = myco.get("context") or {}
    out: list[Any] = []
    frac = ctx.get("sample_fungal_fraction")
    hmp = ctx.get("hmp_shotgun_fungal_fraction") or 0.0001
    out.append(Paragraph(
        f"<b>Is this a normal amount of fungal DNA?</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>context, not a target</font>", st["h3"]))
    if frac is not None:
        out.append(Paragraph(
            f"In the largest survey of healthy adults (HMP, 147 volunteers), about <b>{100 * hmp:.2f}%</b> of shotgun reads mapped "
            f"to fungal genomes. Yours: <b>{100 * frac:.4f}%</b> \u2014 {ctx.get('position_versus_hmp', '')}. Fungi are a small "
            "fraction of stool DNA in everyone; there is no validated healthy range, and less is not better or worse. Genome "
            "size, cell walls and yesterday's meal all move this number.", st["small"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("GENUS", "HEALTHY ADULTS", "USUALLY MEANS", "YOU")]]
    state_words = {"supported": ("found", BLUE), "provisional": ("provisional", SLATE), "trace": ("trace", SLATE),
                   "ambiguous_complex": ("found (genus)", BLUE), "not_detected": ("not detected", INK_FAINT)}
    for r in ctx.get("genus_prevalence") or []:
        w, col = state_words.get(r["found_here"], (r["found_here"], SLATE))
        prev = f"in {r['hmp_prevalence_percent']:.0f}% of samples" if r.get("hmp_prevalence_percent") else "commonly reported"
        rows.append([Paragraph(f"<i><b>{r['genus']}</b></i>", st["cell"]), Paragraph(f"<font size='6.8'>{prev}</font>", st["cell"]),
                     Paragraph(f"<font size='6.8'>{r['meaning']}</font>", st["cell"]),
                     Paragraph(f"<font color='{_hex(col)}'><b>{w}</b></font>", st["cell"])])
    t = Table(rows, colWidths=[CONTENT_WIDTH * w for w in (0.20, 0.22, 0.40, 0.18)], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                           ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                           ("TOPPADDING", (0, 0), (-1, -1), 1.8), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    out.append(t)
    out.append(Paragraph(f"<font size='6.2' color='{_hex(INK_FAINT)}'>{ctx.get('note', '')} Source: Nash et al. 2017, Microbiome.</font>", st["small"]))
    return out


class ClassBar(Flowable):
    """One stacked bar: share of fungal DNA by organism class."""

    COLOURS = {"food_and_drink": BLUE, "opportunist": AMBER, "beneficial_preclinical": GREEN,
               "environmental": colors.HexColor("#8FA3B8"), "pathogen_route": colors.HexColor("#C0392B"), "unknown": GREY}

    def __init__(self, classes: list[Mapping[str, Any]], *, width: float) -> None:
        super().__init__()
        self.classes, self.width, self.height = classes, width, 11 * mm

    def wrap(self, _aw, _ah):
        return self.width, self.height

    def draw(self) -> None:
        c = self.canv
        y, h = 5.5 * mm, 4.2 * mm
        x = 0.0
        total = sum(float(k.get("share_of_fungal_fragments") or 0) for k in self.classes) or 1.0
        c.setFillColor(PANEL_BG)
        c.roundRect(0, y, self.width, h, 1.5 * mm, stroke=0, fill=1)
        for k in self.classes:
            share = float(k.get("share_of_fungal_fragments") or 0) / total
            if share <= 0:
                continue
            c.setFillColor(self.COLOURS.get(k["key"], GREY))
            c.rect(x, y, self.width * share, h, stroke=0, fill=1)
            x += self.width * share
        lx = 0.0
        c.setFont("Helvetica", 6.2)
        for k in self.classes:
            share = float(k.get("share_of_fungal_fragments") or 0)
            if not k.get("taxa") and share <= 0:
                continue
            c.setFillColor(self.COLOURS.get(k["key"], GREY))
            c.circle(lx + 1.2 * mm, 1.9 * mm, 1.1 * mm, stroke=0, fill=1)
            c.setFillColor(INK_SOFT)
            label = f"{k['label']} {100 * share:.0f}%"
            c.drawString(lx + 3.0 * mm, 1.2 * mm, label)
            lx += c.stringWidth(label, "Helvetica", 6.2) + 7 * mm


def score_panel(st: Mapping[str, ParagraphStyle], myco: Mapping[str, Any], *, width: float) -> list[Any]:
    """The Myco-Score: the standard percentile bar (higher is healthier), the
    number, and one line saying what moved it."""
    from openbiota.pdfreport import PercentileBar, status_for

    ms = myco.get("myco_score") or {}
    out: list[Any] = []
    value = ms.get("value")
    if value is None:
        out.append(Paragraph(
            f"<font size='11' color='{_hex(INK_SOFT)}'><b>Myco-Score not computed</b></font> &nbsp;"
            f"<font size='7' color='{_hex(INK_FAINT)}'>{(ms.get('reason') or 'analysis incomplete').replace('_', ' ')}</font>", st["h2"]))
        return out
    value = float(value)
    colour = status_for(percentile=value, higher_means="favourable").colour
    comps = ms.get("components") or {}
    words = {
        "no_opportunists": "opportunistic fungi absent or very little", "opportunist_dominance": "opportunistic fungi past a tenth of the fungal DNA",
        "opportunist_supported": "an opportunist confirmed at species level", "saccharomyces_present": "Saccharomyces present",
        "evidence_concern": "a concerning evidence rule", "evidence_benefit": "a beneficial evidence rule",
        "alert_route_ceiling": "an alert-route pathogen",
    }
    minus = "\u2212"
    moved = [f"{'+' if v > 0 else minus}{abs(v):.0f} {words.get(k, k)}" for k, v in comps.items()
             if k not in ("base", "alert_route_ceiling") and abs(v) >= 0.5]
    if "alert_route_ceiling" in comps:
        moved.append("capped at 20: an alert-route pathogen")
    row = Table([[
        Paragraph("<b>Myco-Score</b> <font color='#8C99A6'>(experimental)</font>", st["h2"]),
        PercentileBar(width=width * 0.46, percentile=value, higher_means="favourable", marker_colour=colour, height=7.0 * mm, show_scale=True),
        Paragraph(f"<font size='20' color='{_hex(colour)}'><b>{int(value + 0.5)}</b></font>"
                  f"<font size='8' color='{_hex(INK_FAINT)}'>/100</font>", st["pct"]),
    ]], colWidths=[width * 0.30, width * 0.50, width * 0.20], hAlign="LEFT")
    row.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                             ("RIGHTPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    out.append(row)
    out.append(Paragraph(
        f"<font size='7' color='{_hex(INK_SOFT)}'>Higher is healthier. From a neutral 50: " + ("; ".join(moved) if moved else "nothing moved it")
        + f". <font color='{_hex(INK_FAINT)}'>Policy weights, versioned ({ms.get('model_id', '')}); the evidence rules and every "
        "contribution are traced in the detail section.</font></font>", st["small"]))
    return out


def _pct_same_scale(frac: float | None) -> str:
    """All percentages of read pairs on one scale, 4 decimals, so 0.0096%
    and 0.0001% read against each other."""
    if frac is None:
        return "\u2014"
    if frac == 0:
        return "0%"
    if 100 * frac < 0.00005:
        return "<0.0001%"
    return f"{100 * frac:.4f}%"


def _tiles(myco: Mapping[str, Any]) -> Table:
    m = myco.get("measurement") or {}
    ctx = myco.get("context") or {}
    sr = myco.get("searched") or {}
    taxa = myco.get("taxa") or []
    n_sup = sum(1 for t in taxa if t["detection_state"] == "supported")
    n_gen = sum(1 for t in taxa if t["detection_state"] == "ambiguous_complex")
    n_prov = sum(1 for t in taxa if t["detection_state"] == "provisional")
    n_trace = sum(1 for t in taxa if t["detection_state"] == "trace")
    n_found = n_sup + n_gen + n_prov + n_trace
    w = (CONTENT_WIDTH - 3 * 3 * mm) / 4
    frac = m.get("fungal_fragment_fraction")
    f = m.get("fungal_supported_fragments")
    opp_all = ctx.get("opportunist_fraction_of_all_fragments")
    opp_taxa = ctx.get("opportunist_taxa") or []
    opp_note = (f"of all read pairs, {len(opp_taxa)} found" if opp_taxa else "none found in your sample")
    tiles = [
        Tile(value=_pct_same_scale(frac), label="Fungal DNA",
             note=(f"{f:,} of {(m.get('eligible_fragments') or 0) / 1e6:.1f}M read pairs" if f is not None
                   else m.get("quantification_status", "")),
             width=w, height=16 * mm, accent=BLUE),
        Tile(value=str(n_found), label="Fungi found",
             note=f"{n_sup} confirmed, {n_gen + n_prov} likely, {n_trace} trace",
             width=w, height=16 * mm, accent=BLUE),
        Tile(value=_pct_same_scale(opp_all if opp_all is not None else 0.0), label="Opportunistic fungi",
             note=opp_note, width=w, height=16 * mm, accent=(GREEN if not opp_taxa else AMBER)),
        Tile(value=f"{sr.get('species_total', 0):,}", label="Species searched",
             note=f"{sr.get('assemblies_total', 0):,} genomes, 3 methods", width=w, height=16 * mm, accent=SLATE),
    ]
    t = Table([tiles], colWidths=[w + 3 * mm] * 3 + [w], hAlign="LEFT")
    t.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                           ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    return t


class CompositionBars(Flowable):
    """Horizontal blue bars: share of the quantified fungal community."""

    def __init__(self, rows: list[Mapping[str, Any]], *, width: float, row_h: float = 5.2 * mm) -> None:
        super().__init__()
        self.rows, self.width, self.row_h = rows[:10], width, row_h
        self.height = self.row_h * max(1, len(self.rows))

    def wrap(self, _aw, _ah):
        return self.width, self.height

    def draw(self) -> None:
        c = self.canv
        label_w = self.width * 0.42
        bar_w = self.width - label_w - 14 * mm
        for i, r in enumerate(self.rows):
            y = self.height - (i + 1) * self.row_h + 1.2 * mm
            c.setFont("Helvetica-Oblique", 7)
            c.setFillColor(INK)
            c.drawString(0, y, r["name"][:44])
            c.setFillColor(BLUE_BG)
            c.roundRect(label_w, y - 0.4 * mm, bar_w, 2.6 * mm, 1.3 * mm, stroke=0, fill=1)
            c.setFillColor(BLUE)
            fw = max(1.0, bar_w * float(r["share_of_quantified_fungi"]))
            c.roundRect(label_w, y - 0.4 * mm, fw, 2.6 * mm, 1.3 * mm, stroke=0, fill=1)
            c.setFont("Helvetica", 6.6)
            c.setFillColor(INK_SOFT)
            c.drawRightString(self.width, y, f"{100 * float(r['share_of_quantified_fungi']):.0f}%  ({r['fragments']} frag.)")


def _taxa_table(st, taxa: list[Mapping[str, Any]], comp_share: Mapping[str, float]) -> Table:
    header = ["FUNGUS", "WHAT IT IS", "DETECTION", "FRAGMENTS", "SHARE", "STRAIN", "EVIDENCE"]
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in header]]
    for t in taxa:
        g = t.get("genome_lane") or {}
        mk = t.get("marker_lane") or {}
        frag = g.get("fragments")
        basis = "both lanes" if (g and mk) else ("genome" if g else ("markers" if mk else "pathogen screen"))
        share = comp_share.get(t["finding_id"])
        strain = (t.get("strain") or {}).get("resolution", "not_assessed").replace("_", " ")
        reason = ((t.get("strain") or {}).get("reason_codes") or [""])[0].replace("_", " ")
        rows.append([
            Paragraph(f"<i><b>{t['accepted_name']}</b></i>" + (f"<br/><font size='6' color='{_hex(INK_FAINT)}'>reported as {t['name']}</font>" if t['name'] != t['accepted_name'] else ""), st["cell"]),
            Paragraph(f"<font size='6.6'>{t.get('category_words', '')}</font>", st["cell"]),
            StatusChip(_state_status(t["detection_state"]), width=20 * mm, font_size=5.8),
            Paragraph(f"<b>{frag if frag is not None else (mk.get('total_reads') if mk else DASH)}</b><br/><font size='5.8' color='{_hex(INK_FAINT)}'>{basis}</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(BLUE)}'><b>{100 * share:.0f}%</b></font>" if share is not None else f"<font color='{_hex(INK_FAINT)}'>\u2014</font>", st["cell"]),
            Paragraph(f"<font size='6.6'>{strain}</font>" + (f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>{reason}</font>" if reason else ""), st["cell"]),
            Paragraph(_evidence_dots(t), st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.25, 0.18, 0.12, 0.12, 0.08, 0.16, 0.09)]
    tb = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    tb.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE),
                            ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                            ("RIGHTPADDING", (0, 0), (-1, -1), 3), ("TOPPADDING", (0, 0), (-1, -1), 2.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.8)]))
    return tb


def _findings_columns(st, myco: Mapping[str, Any]) -> Table:
    """Beneficial evidence and potential concerns, side by side."""
    def cell(title: str, items: list[Mapping[str, Any]], colour) -> list[Any]:
        out: list[Any] = [Paragraph(f"<font color='{_hex(colour)}'><b>{title}</b></font>", st["h3"])]
        if not items:
            out.append(Paragraph("<i>None activated under the current policy.</i> Unscored means unscored, not benign.", st["small"]))
        for c in items:
            tag = {"active": "counts", "unresolved": "range only"}.get(c["activation_status"], c["activation_status"])
            out.append(Paragraph(
                f"<b>{c['rule_id']}</b> <font size='6.4' color='{_hex(INK_FAINT)}'>cap {c['base_cap']:.3g} \u00b7 {tag}</font><br/>"
                f"<font size='6.8'>{c.get('reason', '')}</font>"
                + (f"<br/><font size='6' color='{_hex(INK_FAINT)}'>extrapolation: {', '.join(x.replace('_', ' ') for x in c.get('extrapolations') or [])}</font>"
                   if c.get("extrapolations") else ""), st["small"]))
        return out
    left = cell("Potentially beneficial evidence", myco.get("beneficial_evidence_findings") or [], GREEN)
    right = cell("Potential concerns", myco.get("potential_concern_findings") or [], AMBER)
    t = Table([[left, right]], colWidths=[CONTENT_WIDTH / 2 - 2 * mm, CONTENT_WIDTH / 2 - 2 * mm], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                           ("RIGHTPADDING", (0, 0), (-1, -1), 4 * mm), ("LINEBEFORE", (1, 0), (1, 0), 0.4, RULE),
                           ("LEFTPADDING", (1, 0), (1, 0), 4 * mm)]))
    return t


def mycobiome_overview(story: list[Any], st: Mapping[str, ParagraphStyle], *, section: int,
                       myco: Mapping[str, Any] | None, detail_section: int) -> None:
    story.extend(section_heading(section, "Your gut mycobiome", st["h1"]))
    story.append(Paragraph("The fungi detected in your stool \u2014 and what the evidence says about them.", st["body"]))
    if not myco or myco.get("analysis_status") in (None, "not_assessed"):
        story.append(Paragraph(
            "The mycobiome module did not run for this sample" + (f": {(myco or {}).get('limits', [''])[0]}" if myco else "."),
            st["body"]))
        return
    story.extend(score_panel(st, myco, width=CONTENT_WIDTH))
    story.append(Spacer(1, 2 * mm))
    story.append(_tiles(myco))
    m = myco.get("measurement") or {}
    n = m.get("eligible_fragments")
    sr = myco.get("searched") or {}
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_FAINT)}'>Every one of the {n:,} quality-passed read pairs was searched against "
        f"{sr.get('species_total', 0):,} fungal species ({sr.get('assemblies_total', 0):,} genomes, {sr.get('genera_total', 0):,} genera) "
        "plus a single-copy-marker search and the targeted pathogen screen. Fungal signal is the share of all pairs confidently "
        "assigned to fungi after competition against near neighbours and bacterial decoys.</font>" if n else "", st["small"]))

    taxa = [t for t in myco.get("taxa") or [] if t["detection_state"] in ("supported", "ambiguous_complex", "provisional")]
    traces = [t for t in myco.get("taxa") or [] if t["detection_state"] == "trace"]
    comp = myco.get("within_fungi") or {}
    share = {r["key"]: float(r["share_of_quantified_fungi"]) for r in comp.get("composition") or []}

    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("<b>What was searched for, and what was found</b>", st["h2"]))
    story.append(_searched_table(st, myco))

    story.append(Spacer(1, 2.5 * mm))
    story.extend(_healthy_context(st, myco))

    story.append(CondPageBreak(70 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph(f"<b>Fungi detected</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>{len(taxa)} finding{'s' if len(taxa) != 1 else ''}; "
                           f"each has a page in section {detail_section}</font>", st["h2"]))
    if taxa:
        story.append(_taxa_table(st, taxa, share))
    else:
        story.append(Paragraph("No fungus reached the support rule (20 fragments over 3 separated genomic regions).", st["body"]))
    if traces:
        bits = []
        for t in traces:
            frag = (t.get("genome_lane") or {}).get("fragments") or (t.get("pathogen_screen") or {}).get("fragments") or 0
            src = "pathogen screen" if t.get("support_basis") == "pathogen_screen_only" else "whole genome"
            bits.append(f"<i>{t['accepted_name']}</i> {frag} fragment{'s' if frag != 1 else ''} ({src})")
        story.append(Paragraph(
            f"<font size='6.8'><font color='{_hex(INK_SOFT)}'><b>Traces, below the rule and not interpreted:</b></font> "
            + "; ".join(bits) + ". One or two reads consistent with a species is a trace, not a detection.</font>", st["small"]))

    ctx = myco.get("context") or {}
    if ctx.get("class_composition"):
        story.append(Spacer(1, 2.5 * mm))
        story.append(Paragraph(
            f"<b>What kind of fungi</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>share of fungal DNA by organism class</font>", st["h2"]))
        story.append(ClassBar(ctx["class_composition"], width=CONTENT_WIDTH))
    if comp.get("composition"):
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(
            f"<b>Composition of the quantified fungal community</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>"
            "genome-coverage shares among supported species; not cell counts</font>", st["h3"]))
        story.append(CompositionBars(comp["composition"], width=CONTENT_WIDTH))
        if comp.get("unresolved_fragments"):
            story.append(Paragraph(
                f"<font size='6.6' color='{_hex(INK_FAINT)}'>Plus {comp['unresolved_fragments']} fungal fragments placed only at "
                "genus level or below the support rule \u2014 counted in the signal, not inserted into the shares.</font>", st["small"]))
    story.append(Spacer(1, 2.5 * mm))
    story.append(_findings_columns(st, myco))
    cr = myco.get("colonisation_resistance") or {}
    if cr.get("bacteria"):
        bits = []
        for b in cr["bacteria"]:
            if b["present"]:
                # A carrier at 0.0009% must not print as "0.00%", which reads
                # as absent next to a 95th-percentile rank.
                pct = float(b["percent"] or 0.0)
                shown = f"{pct:.2f}%" if pct >= 0.01 else "<0.01%"
                lvl = f" ({shown}" + (f", {b['percentile']:.0f}th pct" if b.get("percentile") is not None else "") + ")"
                bits.append(f"<font color='{_hex(GREEN)}'>\u25cf</font> <i>{b['bacterium']}</i>{lvl}")
            else:
                bits.append(f"<font color='{_hex(INK_FAINT)}'>{DASH}</font> <i>{b['bacterium']}</i> not detected")
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(
            f"<b>Bacteria that hold Candida in check</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>from your bacterial results; "
            "mouse and laboratory evidence</font>", st["h3"]))
        joined = " &nbsp;\u00b7&nbsp; ".join(bits)
        story.append(Paragraph(
            f"<font size='6.8'>{joined}. "
            + ("Relevant here because an opportunistic yeast was found. " if cr.get("relevant") else "")
            + "Depleting these bacteria (antibiotics especially) is the best-described route to Candida expansion; "
            "their presence is context, not protection proven in you.</font>", st["small"]))
    steps = myco.get("next_steps") or []
    if steps:
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph("<b>What to do next</b>", st["h3"]))
        for x in steps[:3]:
            story.append(Paragraph(f"<font size='6.8'>\u2022 {x['text']}</font>", st["small"]))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_FAINT)}'>Section {detail_section} has each fungus in full: lanes, strain loci searched, "
        "evidence with its label, and the complete score trace.</font>", st["small"]))


def _label_chip(label: str) -> str:
    words = {"human_trial": "human trial", "human_association": "human association",
             "human_isolate_plus_experiment": "human isolate + experiment", "animal_experiment": "animal experiment",
             "in_vitro": "laboratory", "genome_inference": "genome inference", "insufficient_evidence": "insufficient"}.get(label, label)
    col = {"human_trial": GREEN, "human_isolate_plus_experiment": INK, "human_association": INK,
           "animal_experiment": AMBER, "in_vitro": SLATE}.get(label, SLATE)
    return f"<font color='{_hex(col)}' size='6.4'><b>{words}</b></font>"


def _kv(st, rows: list[tuple[str, str]]) -> Table:
    t = Table([[Paragraph(k, st["label"]), Paragraph(v, st["cell"])] for k, v in rows],
              colWidths=[34 * mm, CONTENT_WIDTH - 34 * mm], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                           ("TOPPADDING", (0, 0), (-1, -1), 1.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2)]))
    return t


def _taxon_card(story: list[Any], st, t: Mapping[str, Any]) -> None:
    g = t.get("genome_lane") or {}
    mk = t.get("marker_lane") or {}
    sr = t.get("strain") or {}
    block: list[Any] = [Anchor(f"myco-{t['finding_id']}", above=6 * mm), back_link_line(section_dest("mycobiome"), "at a glance")]
    head = Table([[Paragraph(f"<i><b>{t['accepted_name']}</b></i> &nbsp;<font size='7.4'>{t.get('category_words', '')}</font>", st["h3"]),
                   StatusChip(_state_status(t["detection_state"]), width=40 * mm, height=5.6 * mm, font_size=6.6)]],
                 colWidths=[CONTENT_WIDTH - 42 * mm, 42 * mm], hAlign="LEFT")
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("LINEBELOW", (0, 0), (-1, 0), 0.8, BLUE), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
    block.append(head)
    block.append(Paragraph("<b>What it is</b>", st["h3"]))
    block.append(Paragraph(t.get("role") or "No curated description; listed because it was found.", st["body_ink"]))
    block.append(Paragraph("<b>How it was found</b>", st["h3"]))
    kv: list[tuple[str, str]] = []
    if g:
        kv.append(("WHOLE GENOME", f"<b>{g['fragments']}</b> fragments over <b>{g['separated_regions_100kb']}</b> separated regions on "
                                  f"{g['contigs_hit']} contigs \u00b7 identity {100 * g['identity_mean']:.1f}% \u00b7 backbone {g.get('backbone_accession') or DASH}"))
    if mk:
        kv.append(("MARKERS", f"EukDetect2: {mk.get('observed_markers') or '?'} markers, {mk['total_reads']} reads, RPKS {mk['rpks_reads_per_kb_marker']:.3g} "
                              "(marker-length-normalised support, not load)"))
    if t.get("reasons"):
        kv.append(("NOTE", ", ".join(r.replace("_", " ") for r in t["reasons"])))
    kv.append(("ORIGIN", ", ".join(o.replace("_", " ") for o in t.get("origin_interpretation") or ["unknown"]) + " \u2014 a single stool sample cannot prove residence"))
    block.append(_kv(st, kv))
    # Two internal tokens read as one phrase and printed the raw second one:
    # "Strain not assessed - analysis not_assessed". The status only adds
    # something when it differs from the resolution.
    resolution = str(sr.get("resolution") or "not assessed").replace("_", " ")
    status = str(sr.get("analysis_status") or "").replace("_", " ")
    strain_line = resolution if status in ("", resolution) else f"{resolution} \u00b7 {status}"
    block.append(Paragraph(
        f"<b>Strain</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>{strain_line}</font>", st["h3"]))
    skv: list[tuple[str, str]] = []
    if sr.get("eligible_discriminatory_loci") is not None:
        skv.append(("LOCI", f"{sr.get('callable_discriminatory_loci', 0)} of {sr['eligible_discriminatory_loci']} discriminating sites callable "
                            f"in {sr.get('callable_blocks', 0)} linkage blocks \u00b7 panel {sr.get('panel_id')}"))
    if sr.get("mean_depth") is not None:
        skv.append(("COVERAGE", f"mean depth {sr['mean_depth']:.4f}x \u00b7 {sr.get('callable_nuclear_bases', 0):,} bases at \u22653x of "
                                f"{sr.get('backbone_bases', 0):,}" + (f" \u00b7 expected {sr['expected_depth']:.4f}x from the fungal fraction" if sr.get("expected_depth") else "")))
    if sr.get("reference_equivalence_group"):
        grp = sr["reference_equivalence_group"]
        skv.append(("COMPATIBLE", f"{len(grp)} reference genotypes remain compatible" + (f": {', '.join(grp[:4])}" + (" \u2026" if len(grp) > 4 else "") if len(grp) <= 12 else "")))
    if sr.get("mixture_status") and sr["mixture_status"] != "not_assessed":
        skv.append(("MIXTURE", sr["mixture_status"].replace("_", " ")))
    for fu in sr.get("recommended_analytical_followup") or []:
        skv.append(("NEXT", fu))
    if sr.get("error"):
        skv.append(("ERROR", sr["error"]))
    block.append(_kv(st, skv) if skv else Paragraph("Strain analysis was not attempted for this finding.", st["small"]))
    ev = t.get("evidence") or []
    block.append(Paragraph(f"<b>What the evidence says</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>{len(ev)} claim{'s' if len(ev) != 1 else ''}</font>", st["h3"]))
    if ev:
        rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("DIRECTION", "ENDPOINT", "EVIDENCE", "LIMIT")]]
        for c in ev:
            d = c.get("direction", "context")
            dcol = {"favorable": GREEN, "concerning": AMBER}.get(d, SLATE)
            rows.append([Paragraph(f"<font color='{_hex(dcol)}'><b>{d}</b></font>", st["cell"]),
                         Paragraph(f"<font size='6.8'>{c.get('endpoint', '')}<br/><font color='{_hex(INK_FAINT)}' size='6'>{c.get('host', '')} \u00b7 {c.get('site', '')} \u00b7 {c.get('design', '')}</font></font>", st["cell"]),
                         Paragraph(_label_chip(c.get("label", "")) + (f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>{c.get('source', '')}</font>" if c.get("source") else ""), st["cell"]),
                         Paragraph(f"<font size='6.6'>{c.get('limit', '')}</font>", st["cell"])])
        et = Table(rows, colWidths=[CONTENT_WIDTH * w for w in (0.11, 0.34, 0.20, 0.35)], hAlign="LEFT")
        et.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                                ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                                ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4)]))
        block.append(et)
    else:
        block.append(Paragraph("No health evidence is held for this organism. It is listed because it was found.", st["small"]))
    if t.get("must_not_infer"):
        block.append(Paragraph(f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>Not inferred from this finding:</b> {t['must_not_infer']}</font>", st["small"]))
    story.append(KeepTogether(block))
    story.append(Spacer(1, 3 * mm))


MYCO_COMPONENT_WORDS: Final = {
    "base": ("Neutral start", "Every completed analysis starts at 50."),
    "no_opportunists": ("Opportunists absent or very little",
                        "Under a tenth of the fungal DNA and under 0.001% of all read pairs, none confirmed at species level."),
    "opportunist_dominance": ("Opportunist dominance",
                              "Scaled from 0 at a tenth of the fungal DNA to \u221230 when opportunists are all of it: the "
                              "\u201cdominated by Candida\u201d pattern the literature ties to IBD relapse, alcohol-related liver "
                              "disease, cirrhosis and failed FMT."),
    "opportunist_supported": ("Opportunist confirmed", "At least one opportunistic species passed the full support test."),
    "saccharomyces_present": ("Saccharomyces present",
                              "Supported S. cerevisiae: depleted in IBD, liver disease and obesity in cohorts; S. boulardii has trial "
                              "evidence. Largely dietary, so a small weight."),
    "evidence_concern": ("Evidence rules, concerning", "\u221230 \u00d7 C from the MHS-E1 rule table below."),
    "evidence_benefit": ("Evidence rules, beneficial", "+20 \u00d7 B from the MHS-E1 rule table below."),
    "alert_route_ceiling": ("Alert-route ceiling", "A supported alert-route pathogen caps the score at 20 whatever else is true."),
}


def _myco_score_trace(story: list[Any], st, myco: Mapping[str, Any]) -> None:
    ms = myco.get("myco_score") or {}
    story.append(Paragraph("<b>How the Myco-Score was computed</b> &nbsp;<font size='7' color='{}'>{}</font>".format(
        _hex(INK_FAINT), ms.get("model_id", "")), st["h2"]))
    if ms.get("value") is None:
        story.append(Paragraph(f"Not computed: {(ms.get('reason') or 'analysis incomplete').replace('_', ' ')}.", st["small"]))
        return
    inp = ms.get("inputs") or {}
    story.append(Paragraph(
        f"Inputs: opportunistic fungi are {_pct_same_scale(inp.get('opportunist_fraction_of_all_fragments'))} of all read pairs and "
        f"{100 * float(inp.get('opportunist_share_of_fungal_dna') or 0):.1f}% of the fungal DNA; {inp.get('opportunists_supported', 0)} "
        f"opportunist(s) confirmed at species level; {inp.get('alert_route_supported', 0)} alert-route pathogen(s); Saccharomyces "
        f"{'supported' if inp.get('saccharomyces_supported') else 'not supported'}; evidence rules B = {float(inp.get('B') or 0):.2g}, "
        f"C = {float(inp.get('C') or 0):.2g}.", st["small"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("COMPONENT", "POINTS", "WHAT IT MEANS")]]
    minus = "\u2212"
    for k, v in (ms.get("components") or {}).items():
        title, why = MYCO_COMPONENT_WORDS.get(k, (k.replace("_", " "), ""))
        if k == "alert_route_ceiling":
            pts, col = f"cap {v:.0f}", CORAL
        else:
            pts = f"{'+' if v > 0 else (minus if v < 0 else '')}{abs(v):.1f}" if k != "base" else f"{v:.0f}"
            col = GREEN if v > 0 and k != "base" else (AMBER if v < 0 else INK)
        rows.append([Paragraph(f"<b>{title}</b>", st["cell"]), Paragraph(f"<font color='{_hex(col)}'><b>{pts}</b></font>", st["cell"]),
                     Paragraph(f"<font size='6.4'>{why}</font>", st["cell"])])
    rows.append([Paragraph("<b>Myco-Score</b>", st["cell"]), Paragraph(f"<b>{ms.get('display', 0)}</b>/100", st["cell"]),
                 Paragraph("<font size='6.4'>Sum, held within 0\u2013100, then the alert ceiling if any. Higher is healthier.</font>", st["cell"])])
    t = Table(rows, colWidths=[CONTENT_WIDTH * w for w in (0.26, 0.10, 0.64)], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                           ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LINEABOVE", (0, -1), (-1, -1), 0.6, RULE),
                           ("LEFTPADDING", (0, 0), (-1, -1), 3), ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2)]))
    story.append(t)
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_SOFT)}'>Basis: {ms.get('source_basis', '')}. No component rewards fungal load or diversity: "
        "the review describes fungi at 0.01\u20130.1% of the microbiome with no lower bound, and diversity moves in opposite directions "
        "across diseases. Experimental; not validated against outcomes.</font>", st["small"]))
    story.append(Spacer(1, 2.5 * mm))


def _score_trace(story: list[Any], st, myco: Mapping[str, Any]) -> None:
    _myco_score_trace(story, st, myco)
    hs = myco.get("health_score") or {}
    story.append(Paragraph("<b>The evidence rules behind it</b> &nbsp;<font size='7' color='{}'>MHS-E1 \u00b7 policy {}</font>".format(
        _hex(INK_FAINT), hs.get("policy_lock_id", "")), st["h2"]))
    b, c = hs.get("benefit_contribution"), hs.get("concern_contribution")
    if hs.get("score_100") is not None:
        story.append(Paragraph(
            f"B = {b:.3g}, C = {c:.3g}; score = 100 \u00d7 0.5 \u00d7 (1 \u2212 C) \u00d7 (1 + B) = <b>{hs['score_100']:.2f}</b>. "
            "B and C are each the single strongest supported contribution, not sums: more studies of the same organism do not "
            "add, and a benefit cannot erase the dominant concern. Caps are policy, not measured effects.", st["small"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("RULE", "DIR", "CAP", "ACTIVATION", "STATUS", "WHY")]]
    for r in hs.get("contributions") or []:
        col = {"active": (TEAL if r["direction"] == "B" else AMBER), "unresolved": INK_SOFT}.get(r["activation_status"], INK_FAINT)
        rng = r.get("activation_range")
        act = f"{r['activation']:.2g}" + (f" (range {rng[0]:.0f}\u2013{rng[1]:.0f})" if rng and rng[0] != rng[1] else "")
        rows.append([Paragraph(f"<b>{r['rule_id']}</b>", st["cell"]), Paragraph(r["direction"], st["cell"]),
                     Paragraph(f"{r['base_cap']:.3g}", st["cell"]), Paragraph(act, st["cell"]),
                     Paragraph(f"<font color='{_hex(col)}'><b>{r['activation_status']}</b></font>", st["cell"]),
                     Paragraph(f"<font size='6.4'>{r.get('reason', '')}</font>", st["cell"])])
    t = Table(rows, colWidths=[CONTENT_WIDTH * w for w in (0.22, 0.05, 0.06, 0.13, 0.10, 0.44)], hAlign="LEFT")
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                           ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                           ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2)]))
    story.append(t)
    rr, pr, dr = hs.get("resolution_sensitivity_range") or {}, hs.get("policy_sensitivity_range") or {}, hs.get("display_sensitivity_range") or {}
    if dr:
        story.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>Sensitivity: resolution alternatives {rr.get('low', 0):.1f}\u2013{rr.get('high', 0):.1f}; "
            f"policy caps \u00d70.5/\u00d71.5 {pr.get('low', 0):.1f}\u2013{pr.get('high', 0):.1f}; combined (displayed) "
            f"{dr.get('low', 0):.1f}\u2013{dr.get('high', 0):.1f}. Scenario ranges, not confidence intervals. Interpretation confidence "
            f"is {hs.get('interpretation_confidence', 'experimental_limited').replace('_', ' ')} at any depth.</font>", st["small"]))


def mycobiome_detail(story: list[Any], st: Mapping[str, ParagraphStyle], *, section: int, myco: Mapping[str, Any] | None) -> None:
    story.extend(section_heading(section, "Your gut mycobiome, in detail", st["h1"]))
    if not myco or myco.get("analysis_status") in (None, "not_assessed"):
        story.append(Paragraph("The mycobiome module did not run for this sample.", st["body"]))
        return
    m = myco.get("measurement") or {}
    gl = myco.get("genome_lane") or {}
    ml = myco.get("marker_lane") or {}
    rc = myco.get("reference_counts") or {}
    story.append(Paragraph(
        f"Every quality-passed read pair ({(m.get('eligible_fragments') or 0):,}) was aligned competitively against "
        f"{(gl.get('competitive_index') or {}).get('n_fungal_assemblies', 0)} fungal assemblies chosen from a reference of "
        f"{rc.get('assemblies_admitted_fungal', 0):,} ({rc.get('species_admitted', 0):,} species), their near neighbours, and "
        f"{len((gl.get('competitive_index') or {}).get('decoys') or [])} bacterial and control decoys; a separate single-copy-marker "
        f"search ({ml.get('tool', 'EukDetect2')}) ran on the same reads. A fungus is supported when its fragments are spread across "
        "separated regions of its genome at high identity and no decoy competes for them. Below: each fungus found, then the "
        "score trace.", st["body"]))
    led = m
    if led.get("fungal_supported_fragments") is not None:
        story.append(Paragraph(
            f"<font size='6.8'>Ledger: {led['fungal_supported_fragments']} fungal \u00b7 {led.get('cross_kingdom_ambiguous_fragments', 0)} "
            f"cross-kingdom ambiguous (excluded) \u00b7 {(gl.get('ledger') or {}).get('fungal_below_identity_fragments', 0)} fungal hits below "
            f"{100 * 0.97:.0f}% identity (excluded) \u00b7 {(gl.get('ledger') or {}).get('decoy_fragments', 0):,} decoy \u00b7 "
            f"{led.get('unassigned_fragments', 0):,} unassigned. Each pair counted once.</font>", st["small"]))
    taxa = [t for t in myco.get("taxa") or [] if t["detection_state"] in ("supported", "ambiguous_complex", "provisional")]
    for t in taxa:
        story.append(CondPageBreak(70 * mm))
        _taxon_card(story, st, t)
    story.append(CondPageBreak(60 * mm))
    _score_trace(story, st, myco)
    lims = myco.get("limits") or []
    if lims:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph("<b>Limits</b>", st["h3"]))
        for lim in lims:
            story.append(Paragraph(f"<font size='6.6' color='{_hex(INK_SOFT)}'>\u2022 {lim}</font>", st["small"]))


__all__ = ["CompositionBars", "MhsGauge", "mycobiome_detail", "mycobiome_overview", "score_panel"]
