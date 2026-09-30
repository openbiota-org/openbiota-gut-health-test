"""Organisms that need attention: one list, two sections, no gaps.

Part A's glance and Part B's cards are drawn from the same list, built here
once. Every organism the classification flags - an opportunist expanded, a
beneficial organism depleted, a conditional one well past its range - is in
it, and so is every organism the census expected and did not find. Nothing
is truncated. Two sections that used to be built from two different systems
disagreed by twenty-two organisms on one sample; they cannot now.

Each card answers, in this order, what a reader asks when they click:

1. What is this organism and what does it do in a gut.
2. Why is it flagged in *this* sample - the level, the percentile, the class.
3. What has been shown to move it the way the reader would want, with the
   evidence level of every item stated, and what competes with it.
4. What the intervention registry holds for it, where anything does.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, KeepTogether, Spacer, Table, TableStyle

from openbiota import modulators as mods
from openbiota import organisms as org
from openbiota.inventory import Organism
from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    back_link_line,
    linked_cell,
    organism_dest,
    section_heading,
)
from openbiota.pdfreport import (
    AMBER,
    AMBER_BG,
    CONTENT_WIDTH,
    CORAL,
    CORAL_BG,
    GREEN,
    INK,
    INK_FAINT,
    INK_SOFT,
    ORANGE,
    ORANGE_BG,
    PANEL_BG,
    RULE,
    SLATE,
    SLATE_BG,
    PercentileBar,
    Status,
    StatusChip,
    Tile,
)

CLASS_COLOUR: Final = {
    org.BENEFICIAL: GREEN, org.OPPORTUNIST: CORAL, org.CONDITIONAL: AMBER, org.UNKNOWN: SLATE,
}


def _hex(c: colors.Color) -> str:
    return c.hexval().replace("0x", "#")


def _pct(v: float) -> str:
    return f"{v:.3f}%" if v < 1 else f"{v:.2f}%"


def _ordinal(p: float | None) -> str:
    return org._ordinal(p) if p is not None else "\u2014"



# --------------------------------------------------------------------------- #
# the one list
# --------------------------------------------------------------------------- #


#: The share below which an organism that was not found is stated to be.
#: A fixed floor, not a modelled limit: the census's prior detection curve
#: put the 95% limit at 0.04% for samples whose own pages list organisms
#: at 0.0002%. Samples resolve organisms to 0.0001%; a missing one is below that.
MISSING_FLOOR_PERCENT: Final = 0.0001


@dataclass(frozen=True, slots=True)
class Attention:
    """One organism that needs attention, from whichever system saw it."""

    species: str
    display: str
    #: "overgrown" | "missing" | "depleted" | "watch"
    group: str
    cls: str
    percent: float
    percentile: float | None
    prevalence: float | None
    detected: bool
    #: Importance weight, for ordering within a group.
    weight: float
    flag_reason: str
    verdict: org.Verdict | None
    finding: Any | None
    modulators: mods.Modulators | None
    detection_power: float | None = None
    #: Percent above (+) or below (-) the typical carrier's level in the reference.
    deviation: float | None = None
    #: The typical carrier's level (median among reference carriers), on the reference lane's scale.
    typical_percent: float | None = None
    #: False for a population counted within a relative's share in the
    #: composition: it has no share of its own, and `percent` is a lane reading.
    in_primary: bool = True
    #: The relative whose share it is counted within, when `in_primary` is False.
    counted_within: str | None = None
    #: The reading the percentile was taken on, in the units of the lane that
    #: ranked it (the whole-genome lane's own percent for a population the
    #: marker catalogue did not place).
    lane_reading: float | None = None
    #: For an organism listed as missing: what the wider inventory has to say
    #: about it - an unconfirmed single-method hit, or a relative that is present.
    note: str | None = None

    @property
    def level_text(self) -> str:
        """The organism's share as the page states it: its measured share,
        or for an organism that was not found, below the reporting floor."""
        if self.detected:
            return _pct(self.percent)
        return f"&lt;{MISSING_FLOOR_PERCENT:g}%"

    @property
    def status(self) -> Status:
        if self.group == "overgrown":
            return Status("overgrown", CORAL, CORAL_BG, self.flag_reason, 3)
        if self.group == "missing":
            return Status("missing", CORAL if self.cls == org.BENEFICIAL else SLATE,
                          CORAL_BG if self.cls == org.BENEFICIAL else SLATE_BG,
                          "commonly carried; not detected here", -3)
        if self.group == "depleted":
            return Status("low", ORANGE, ORANGE_BG, self.flag_reason, -2)
        if self.flag_reason.startswith("uncommon"):
            return Status("uncommon", AMBER, AMBER_BG, self.flag_reason, 1)
        return Status("high", AMBER, AMBER_BG, self.flag_reason, 1)

    @property
    def colour(self) -> colors.Color:
        return self.status.colour

    @property
    def want(self) -> str | None:
        """Which way the reader would want this organism to move: ``"high"``
        means it is high and lowering it is the aim, ``"low"`` the reverse.

        None for a beneficial organism that is merely high. Akkermansia at
        the 96th percentile is listed so the reader knows, not so they can
        lower it; nothing that raises it is a thing to avoid.
        """
        if self.group == "watch" and (self.cls == org.BENEFICIAL or self.flag_reason.startswith("uncommon")):
            return None
        return "high" if self.group in ("overgrown", "watch") else "low"

    @property
    def counts_against(self) -> bool:
        """Whether a lever pushing this organism the wrong way is worth
        naming. Overgrown, missing and depleted readings yes; a conditional
        resident merely above its range is a note, not a reason to avoid
        yoghurt."""
        return self.group != "watch"


GROUPS: Final[tuple[tuple[str, str, str, colors.Color], ...]] = (
    ("overgrown", "OVERGROWN \u2014 at a level that is a concern",
     "Opportunists above the 90th percentile among reference adults who carry them, and any organism above the 97th. "
     "These are the readings to act on first.", CORAL),
    ("missing", "MISSING \u2014 commonly carried, not detected here",
     "Species most reference adults carry that were not found at a depth where they should have been. "
     "The organisms that should be there.", CORAL),
    ("depleted", "DEPLETED \u2014 well below the reference",
     "Beneficial organisms present but below the 10th percentile of people who carry them.", ORANGE),
    ("watch", "WORTH WATCHING \u2014 above the reference range, or uncommon to carry",
     "Conditional organisms between the 90th and 97th percentile among carriers - normal residents that have "
     "expanded, a signal about diet and balance more than a problem in themselves - and opportunists that fewer "
     "than one in ten reference adults carry at all, whatever their level.", AMBER),
)


def _current_genus(species: str) -> tuple[str, str | None]:
    """(genus word, GTDB genus or None) for a reference-catalogue name today:
    the renamed species where taxonomy moved it (Eubacterium eligens is
    Lachnospira), the GTDB placement where the bridge knows the name, else
    the name's own first word. A relative named under the old genus would be
    no relative; Dorea_A and Dorea_D are two genera."""
    from openbiota.inventory import canonical

    current = canonical(species)
    try:
        from openbiota.expansion import names as _names

        gtdb = _names.ncbi_species_to_r232().get(_names.mpa_style(current)) or ""
    except Exception:  # noqa: BLE001 - a lookup table; its absence must not cost the page
        gtdb = ""
    if gtdb:
        return gtdb.split(" ")[0].split("_")[0], gtdb.split(" ")[0]
    return current.split("_")[0], None


def build(
    inv: Any,
    findings: Any | None,
) -> list[Attention]:
    """Every organism needing attention, from the verdicts and the census."""
    items: list[Attention] = []
    #: Every name under which an already-listed organism is known. The
    #: census writes the reference catalogue's name (Ruminococcus_gnavus);
    #: the inventory writes the current one (Mediterraneibacter_gnavus).
    #: Both must count as seen or the same organism is listed twice.
    seen: set[str] = set()

    def _mark(species: str, *more: str | None) -> None:
        seen.add(species)
        for m in more:
            if m:
                seen.add(m)
                seen.add(m.replace(" ", "_"))

    def _known(species: str) -> bool:
        if species in seen:
            return True
        if inv is not None:
            o = inv.get(species)
            if o is not None and o.species in seen:
                return True
        return False

    if inv is not None and len(inv):
        for v in org.verdicts(list(inv.organisms)):
            if not v.flagged:
                continue
            o = v.organism
            group = (("overgrown" if v.is_issue else "watch") if v.flag == "high"
                     else "watch" if v.flag == "uncommon" else "depleted")
            m = mods.lookup(o.species, gtdb=o.gtdb, formerly=o.formerly)
            flag_reason = v.flag_reason
            if m is None:
                # The expanded catalogues flag organisms the curated registry
                # has never described. The card still says what the organism
                # is - the classification's own description, at whatever rank
                # it was made - and says plainly that no intervention evidence
                # is on record. Without evidence to act on, a high reading is
                # something to watch, not an overgrowth to treat.
                m = mods.Modulators(
                    organism=o.species, role=v.description or "No description on record for this organism.",
                    why_level_matters=("No intervention evidence is on record for this organism; the flag rests on its "
                                       "percentile alone and is listed to watch, not to act on."),
                    decrease=(), increase=(), competitors=(), basis=v.basis or "none")
                if group == "overgrown":
                    group = "watch"
                    flag_reason = f"{flag_reason} (no intervention evidence on record)"
            items.append(Attention(
                species=o.species, display=o.display, group=group, cls=v.cls,
                percent=o.best_percent, percentile=getattr(o, "level_percentile", o.percentile), prevalence=o.prevalence,
                detected=True, weight=org.importance(v), flag_reason=flag_reason,
                verdict=v, finding=None,
                modulators=m, deviation=getattr(o, "deviation_percent", None),
                typical_percent=getattr(o, "reference_percent", None),
                in_primary=bool(o.in_primary), counted_within=getattr(o, "counted_within", None),
                lane_reading=getattr(o, "reference_reading", None),
            ))
            _mark(o.species, o.formerly, o.gtdb, *getattr(o, "aliases", ()))

    if findings is not None:
        # The census knows what the reference catalogue did not find; the
        # inventory, pooling every method, knows what was found under any
        # name. An organism any method found with support is not missing.
        def _found(species: str) -> Any | None:
            o = inv.get(species) if inv is not None else None
            if o is None:
                return None
            if (o.in_primary and o.percent > 0) or (o.status == "supported" and (o.best_percent or 0.0) > 0):
                return o
            return None

        def _missing_note(species: str) -> str | None:
            if inv is None:
                return None
            o = inv.get(species)
            if o is not None and (o.best_percent or 0.0) > 0:
                # a single-method hit the competition did not confirm
                return f"an unconfirmed single-method hit at {_pct(o.best_percent)} is not counted as a find"
            genus, gtdb_genus = _current_genus(species)
            if gtdb_genus:
                # the bridge placed the name in a GTDB genus: kin are that genus exactly
                kin = [x for x in inv.organisms if x.in_primary and x.percent > 0
                       and (x.gtdb_genus or (x.gtdb or "").split(" ")[0]) == gtdb_genus]
            else:
                kin = [x for x in inv.organisms if x.in_primary and x.percent > 0
                       and (x.gtdb_genus or x.genus or x.display.split(" ")[0]).split("_")[0] == genus]
            if kin:
                top = max(kin, key=lambda x: x.percent)
                return f"a relative is present: <i>{top.display}</i> at {_pct(top.percent)}"
            return None

        for t in list(findings.bucket(4)):
            if _known(t.species) or _found(t.species) is not None:
                continue
            # Class the absent organism as if it were present: the registry
            # and guild tables know what it is regardless of the sample.
            v_cls = org.verdict(Organism(species=t.species, percent=0.0, genus=t.species.split("_")[0])).cls
            items.append(Attention(
                species=t.species, display=t.display_name, group="missing", cls=v_cls,
                percent=0.0, percentile=None, prevalence=t.prevalence, detected=False,
                weight=60.0 + 40.0 * float(t.prevalence or 0.0),
                flag_reason="not detected",
                verdict=None, finding=t, modulators=mods.lookup(t.species),
                detection_power=getattr(t, "detection_power", None),
                note=_missing_note(t.species),
            ))
            _mark(t.species)
        # Census high/low that the verdicts did not already flag: organisms
        # the census placed but the merged inventory holds no record of. An
        # organism the inventory does hold is judged there, once, on its one
        # share and its level among carriers; the census's own rank on the
        # reference catalogue's reading is not a second listing of it.
        for bucket, group in ((2, "watch"), (3, "depleted")):
            for t in list(findings.bucket(bucket)):
                if _known(t.species) or not t.detected:
                    continue
                if inv is not None and inv.get(t.species) is not None:
                    continue
                cls = org.verdict(Organism(species=t.species, percent=float(t.percent or 0.0),
                                           genus=t.species.split("_")[0])).cls
                if bucket == 2 and cls == org.OPPORTUNIST:
                    group = "overgrown"
                word = "high" if bucket == 2 else "low"
                items.append(Attention(
                    species=t.species, display=t.display_name, group=group, cls=cls,
                    percent=float(t.percent or 0.0), percentile=t.percentile, prevalence=t.prevalence,
                    detected=True, weight=30.0 + abs((t.percentile or 50.0) - 50.0),
                    flag_reason=f"{word} \u2014 {_ordinal(t.percentile)} percentile among carriers",
                    verdict=None, finding=t, modulators=mods.lookup(t.species),
                ))
                _mark(t.species)

    items.sort(key=lambda a: -a.weight)
    return items


# --------------------------------------------------------------------------- #
# Part A: the glance
# --------------------------------------------------------------------------- #


def glance(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    detail_section: int,
    items: list[Attention],
    findings: Any | None,
    n_all_detected: int,
    cohort_note: str,
) -> None:
    story.extend(section_heading(
        section, "Organisms that need attention: overgrown, missing, low", st["h1"]))
    story.append(Paragraph(
        "Your community, read in both directions. The full inventory of everything detected is the "
        "section before this; this page is the part that inventory cannot show \u2014 <b>what is "
        "expanded past its range</b>, <b>what is not there that usually is</b>, and what is present but "
        "low. Every organism here has its own page in "
        f"section {detail_section}: what it does, why it matters, and what has been shown to move it.",
        st["body"]))

    by = {g: [a for a in items if a.group == g] for g, *_ in GROUPS}
    tile_w = (CONTENT_WIDTH - 3 * 3 * mm) / 4
    tiles = Table([[
        Tile(value=str(len(by["overgrown"])), label="Overgrown", note="opportunists or very high",
             width=tile_w, height=15 * mm, accent=CORAL if by["overgrown"] else GREEN),
        Tile(value=str(len(by["missing"])), label="Missing", note="usually carried, not found",
             width=tile_w, height=15 * mm, accent=CORAL if by["missing"] else GREEN),
        Tile(value=str(len(by["depleted"])), label="Depleted", note="beneficial, well below range",
             width=tile_w, height=15 * mm, accent=ORANGE if by["depleted"] else GREEN),
        Tile(value=str(len(by["watch"])), label="Worth watching", note="expanded, not yet a concern",
             width=tile_w, height=15 * mm, accent=AMBER if by["watch"] else GREEN),
    ]], colWidths=[tile_w + 3 * mm] * 3 + [tile_w], hAlign="LEFT")
    tiles.setStyle(TableStyle([
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(tiles)

    ref = getattr(findings, "reference", None)
    if ref is not None:
        qual = f" {cohort_note.strip()}" if cohort_note and cohort_note.strip() else ""
        story.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_FAINT)}'>Reference: {ref.n_samples:,} adult stool "
            f"metagenomes from {ref.n_studies} studies.{qual} {n_all_detected:,} organisms detected in "
            "all; percentiles are against the adults who carry each organism. An organism with no "
            "percentile is one the reference set has no measurement for, and is flagged on presence "
            "and class alone.</font>", st["small"]))

    for group, heading, blurb, colour in GROUPS:
        rows_g = by[group]
        story.append(Spacer(1, 3 * mm))
        story.append(KeepTogether([
            Paragraph(f"<font color='{_hex(colour)}'><b>{heading}</b></font>"
                      f"  <font size='7' color='{_hex(INK_FAINT)}'>{len(rows_g)}</font>", st["h2"]),
            Paragraph(blurb, st["small"]),
        ]))
        if not rows_g:
            story.append(Paragraph("<i>None in this sample.</i>", st["small"]))
            continue
        _glance_table(story, st, rows_g, group)


_SENTENCE_END: Final = re.compile(r"(?<!\b[A-Z])(?<!\bsp)(?<!\bsubsp)\.\s+(?=[A-Z])")


def _first_sentence(text: str) -> str:
    """The first sentence, not cut at an abbreviation: "A relative of
    R. gnavus that..." keeps the "R." """
    parts = _SENTENCE_END.split(text.strip(), maxsplit=1)
    return parts[0].rstrip(".") if parts else ""


def _glance_table(story: list[Any], st: dict[str, ParagraphStyle], rows_g: list[Attention], group: str) -> None:
    bar_w = CONTENT_WIDTH * 0.22
    widths = [CONTENT_WIDTH * w for w in (0.35, 0.12, 0.11)] + [bar_w] + [CONTENT_WIDTH * 0.08, CONTENT_WIDTH * 0.12]
    header = ["ORGANISM", "CLASS", "SHARE", "YOUR LEVEL", "", "READING"]
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in header]]
    marks: list[RowLink] = []
    for a in rows_g:
        role = _first_sentence(a.modulators.role if a.modulators else "")
        role = (role[:118] + "\u2026") if len(role) > 120 else role
        name = (f"<i><b>{a.display}</b></i><br/><font color='{_hex(INK_SOFT)}' size='6.2'>{role}</font>"
                if role else f"<i><b>{a.display}</b></i>")
        chip = (f"<font color='{_hex(CLASS_COLOUR[a.cls])}' size='8'>\u25cf</font> "
                f"<font size='6.8'>{org.CLASS_LABEL[a.cls]}</font>")
        if a.detected and not a.in_primary:
            # counted within a relative in the composition: no share of its
            # own; the level is the whole-genome lane's reading of it
            share = f"<font size='6.4' color='{_hex(INK_FAINT)}'>no share of its own</font>"
            sub = f"within <i>{a.counted_within}</i>" if a.counted_within else "counted within a relative"
            if a.lane_reading:
                sub += f" \u00b7 whole-genome reading {_pct(a.lane_reading)}"
            if a.prevalence is not None:
                sub += f" \u00b7 carried by {a.prevalence:.0%}"
            if a.typical_percent and a.deviation is not None:
                sub += f" \u00b7 typical {_pct(a.typical_percent)}"
        elif a.detected:
            share = _pct(a.percent)
            sub = f"carried by {a.prevalence:.0%}" if a.prevalence is not None else ""
            if a.typical_percent and a.deviation is not None:
                sub += (" \u00b7 " if sub else "") + f"typical {_pct(a.typical_percent)}"
        else:
            share = a.level_text
            sub = "not detected" + (f" \u00b7 carried by {a.prevalence:.0%}" if a.prevalence is not None else "")
            if a.note:
                sub += f"<br/>{a.note}"
        sub = sub.strip(" \u00b7")
        if sub:
            share += f"<br/><font size='6' color='{_hex(INK_FAINT)}'>{sub}</font>"
        if not a.detected:
            bar: Any = PercentileBar(width=bar_w - 6 * mm, percentile=0.0, higher_means="unclear",
                                     marker_colour=a.colour, height=6 * mm)
            # "~0": below every carrier, and not claimed to be exactly nothing.
            num = f"<font color='{_hex(INK_FAINT)}'>~</font>0"
        elif a.percentile is not None:
            bar = PercentileBar(width=bar_w - 6 * mm, percentile=a.percentile, higher_means="unclear",
                                marker_colour=a.colour, height=6 * mm)
            num = _ordinal(a.percentile)
            if a.deviation is not None:
                from openbiota.pdflibrary import deviation_text
                num = (f"{deviation_text(a.deviation)}<br/><font size='5.6' color='{_hex(INK_FAINT)}'>"
                       f"{_ordinal(a.percentile)} pct</font>")
        else:
            bar = PercentileBar(width=bar_w - 6 * mm, percentile=None, higher_means="unclear",
                                marker_colour=a.colour, height=6 * mm)
            num = "\u2014"
        link = RowLink(href=organism_dest(a.species, detail=True), anchor=organism_dest(a.species))
        marks.append(link)
        rows.append([
            linked_cell(link, Paragraph(name, st["cell"])),
            Paragraph(chip, st["cell"]),
            Paragraph(share, st["cell"]),
            bar,
            Paragraph(f"<font color='{_hex(a.colour)}'><b>{num}</b></font>", st["pct"]),
            StatusChip(a.status, width=widths[5] - 3 * mm, font_size=6.4),
        ])
    t = LinkedTable(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style: list[Any] = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE), ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]
    if group == "overgrown":
        style.append(("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#FDF3F2")))
    t.setStyle(TableStyle(style))
    story.append(t)


# --------------------------------------------------------------------------- #
# Part B: the cards
# --------------------------------------------------------------------------- #


def _card_identifiers(o: Any) -> list[str]:
    """Each lane's identifier, labelled and readable: ``SGB4584``, ``GCF_008121495``,
    ``mOTU 000538``, ``UHGG species 2101``, ``panel BackhedF_2015_..._bin.24`` - not
    the raw tokens with their prefixes and file extensions."""
    from openbiota.expansion import genomes as _genomes

    labels = {"jan26": "SGB", "extended": "SGB (Jun23)", "globdb": "GlobDB", "genome": "GTDB", "motus": "mOTU",
              "kraken": "UHGG species", "rescue": "panel", "singlem": "SingleM"}
    out: list[str] = []
    seen: set[str] = set()
    for lane, raw in (o.native_ids or {}).items():
        for token in str(raw).split(";"):
            val = token.split(":", 1)[-1].strip()
            if not val or val.startswith(("s__", "g__")):
                continue
            if lane == "motus" and val.startswith("mOTUv"):
                val = val.split("_", 1)[-1]
            elif lane in ("rescue", "globdb"):
                val = _genomes.normalise(val)
            elif lane in ("jan26", "extended") and val.startswith("SGB"):
                pass
            base = val.split(".")[0] if val.startswith(("GCA_", "GCF_")) else val
            if base in seen:
                continue
            seen.add(base)
            out.append(val if lane in ("jan26", "extended") or val.startswith(("GCA_", "GCF_"))
                       else f"{labels.get(lane, lane)} {val}")
    return out


def detail(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    items: list[Attention],
    plan: Any | None,
    bacdive: Mapping[str, Any] | None = None,
    strain_analysis: Mapping[str, Any] | None = None,
) -> None:
    from openbiota.pdfatlas import evidence_cards

    story.extend(section_heading(section, "Organisms that need attention, one by one", st["h1"]))
    story.append(Paragraph(
        "One page per organism flagged in Part A. Each says what the organism is and does in a gut, "
        "why it is flagged in your sample, what has been shown to move it in the direction you would "
        "want \u2014 with the strength of that evidence stated on every line \u2014 and what competes "
        "with it. Laboratory and animal findings are listed and labelled as such; they are leads, not "
        "recommendations. Nothing here is a prescription.",
        st["body"]))
    story.append(_evidence_legend(st))

    if not items:
        story.append(Paragraph("<i>No organism needs attention in this sample.</i>", st["body"]))
        return

    for group, heading, _blurb, colour in GROUPS:
        rows_g = [a for a in items if a.group == group]
        if not rows_g:
            continue
        story.append(CondPageBreak(60 * mm))
        story.append(Spacer(1, 3 * mm))
        story.append(Paragraph(f"<font color='{_hex(colour)}'><b>{heading}</b></font>", st["h2"]))
        for a in rows_g:
            _card(story, st, a, plan, evidence_cards, bacdive=bacdive, strain_analysis=strain_analysis)


def _evidence_legend(st: dict[str, ParagraphStyle]) -> Table:
    cells = []
    for key, (words, _rank) in sorted(mods.EVIDENCE.items(), key=lambda kv: -kv[1][1]):
        cells.append(Paragraph(
            f"<font size='6.6'><b>{words.capitalize()}</b></font><br/>"
            f"<font size='6' color='{_hex(INK_SOFT)}'>{_EVIDENCE_MEANING[key]}</font>",
            st["small"]))
    w = CONTENT_WIDTH / 4
    t = Table([cells], colWidths=[w] * 4, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"), ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
        ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("ROUNDEDCORNERS", [3, 3, 3, 3]),
    ]))
    return t


_EVIDENCE_MEANING: Final = {
    "human_trial": "given to people; the organism measured before and after",
    "human_observational": "associations in people; nothing was given",
    "animal": "a mouse or rat gut, not a human one",
    "in_vitro": "a culture dish; no gut at all",
}


def _card(story: list[Any], st: dict[str, ParagraphStyle], a: Attention, plan: Any, evidence_cards: Any,
          *, bacdive: Mapping[str, Any] | None = None, strain_analysis: Mapping[str, Any] | None = None) -> None:
    m = a.modulators
    colour = a.colour
    block: list[Any] = [
        Anchor(organism_dest(a.species, detail=True), above=6 * mm),
        back_link_line(organism_dest(a.species), "at a glance"),
    ]

    # ---- title ------------------------------------------------------------ #
    head = Table([[
        Paragraph(
            f"<i><b>{a.display}</b></i> &nbsp;"
            f"<font color='{_hex(CLASS_COLOUR[a.cls])}' size='8'>\u25cf</font> "
            f"<font size='7.4'>{org.CLASS_LABEL[a.cls]}</font>"
            + (f" <font color='{_hex(INK_FAINT)}' size='6.4'>\u00b7 genus-level knowledge</font>"
               if m is not None and m.basis == "genus" else ""),
            st["h3"]),
        StatusChip(a.status, width=44 * mm, height=5.6 * mm, font_size=6.6),
    ]], colWidths=[CONTENT_WIDTH - 46 * mm, 46 * mm], hAlign="LEFT")
    head.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0), ("LINEBELOW", (0, 0), (-1, 0), 0.8, colour),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    block.append(head)
    block.append(Spacer(1, 1.5 * mm))

    # ---- 1. what it is and does ------------------------------------------- #
    block.append(Paragraph("<b>What it is and what it does</b>", st["h3"]))
    if m is not None and m.role:
        block.append(Paragraph(m.role, st["body_ink"]))
    elif a.verdict is not None and a.verdict.description:
        block.append(Paragraph(a.verdict.description, st["body_ink"]))
    else:
        block.append(Paragraph(
            "No curated description exists for this organism yet. It is listed because it was found "
            "at a level outside the reference range, not because anything is known about its effect.",
            st["body_ink"]))

    # ---- 2. why flagged here --------------------------------------------- #
    block.append(Paragraph("<b>Why it is flagged in your sample</b>", st["h3"]))
    kv: list[tuple[str, str]] = []
    if a.detected and not a.in_primary:
        kv.append(("YOUR LEVEL", "<b>no share of its own</b> \u2014 its reads are counted within "
                   + (f"<i>{a.counted_within}</i>" if a.counted_within else "a relative")
                   + " in the composition"
                   + (f"; the whole-genome lane reads it at <b>{_pct(a.lane_reading)}</b> of its own scale"
                      if a.lane_reading else "")))
        if a.percentile is not None:
            from openbiota.pdflibrary import deviation_text
            kv.append(("POSITION", f"<b>{_ordinal(a.percentile)}</b> percentile among reference adults who carry it, "
                       "on that lane"
                       + (f" \u2014 <b>{deviation_text(a.deviation)}</b> against their typical level"
                          + (f" of {_pct(a.typical_percent)}" if a.typical_percent else "")
                          if a.deviation is not None else "")))
        else:
            kv.append(("POSITION", "no reference percentile \u2014 flagged on class and presence"))
    elif a.detected:
        kv.append(("YOUR LEVEL", f"<b>{_pct(a.percent)}</b> of the community"))
        if a.percentile is not None:
            from openbiota.pdflibrary import deviation_text
            kv.append(("POSITION", f"<b>{_ordinal(a.percentile)}</b> percentile among reference adults who carry it"
                       + (f" \u2014 <b>{deviation_text(a.deviation)}</b> against their typical level"
                          + (f" of {_pct(a.typical_percent)}" if a.typical_percent else "")
                          if a.deviation is not None else "")))
        else:
            kv.append(("POSITION", "no reference percentile \u2014 flagged on class and presence"))
    else:
        kv.append(("YOUR LEVEL", f"<b>{a.level_text}</b> \u2014 not detected; this sample resolves organisms below that level"
                   + (f". In the wider inventory, {a.note}" if a.note else "")))
    if a.prevalence is not None:
        kv.append(("HOW COMMON", f"carried by <b>{a.prevalence:.0%}</b> of reference adults"))
    kv.append(("CLASS", f"{org.CLASS_LABEL[a.cls]} \u2014 {org.CLASS_MEANING[a.cls].split('.')[0]}."))
    kv.extend(_provenance_rows(a, bacdive=bacdive, strain_analysis=strain_analysis))
    kvt = Table([[Paragraph(k, st["label"]), Paragraph(v, st["cell"])] for k, v in kv],
                colWidths=[34 * mm, CONTENT_WIDTH - 34 * mm], hAlign="LEFT")
    kvt.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 1.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2),
    ]))
    block.append(kvt)
    if m is not None and m.why_level_matters:
        block.append(Paragraph(
            f"<font color='{_hex(colour)}'><b>What the level means.</b></font> {m.why_level_matters}",
            st["body_ink"]))

    # ---- 3. what moves it ------------------------------------------------- #
    want = a.want
    if want is None:
        # A beneficial organism that is high: say so, and do not offer to
        # lower it. What sustains it is still worth knowing.
        block.append(Paragraph("<b>Nothing to lower</b>", st["h3"]))
        block.append(Paragraph(
            "This is a beneficial organism at a high level. That is listed so you know, not so you can "
            "change it; nothing that raises it is a thing to avoid. If it is very high alongside a "
            "low-fibre diet the organism may be living on your gut's mucus for want of anything else, "
            "and the answer to that is fibre, not less of this organism.",
            st["body_ink"]))
        keeps = m.for_direction("low") if m is not None else ()
        if keeps:
            block.append(Paragraph(
                f"<font size='6.6' color='{_hex(INK_FAINT)}'>What sustains it: "
                + "; ".join(f"{x.agent} ({x.evidence_words})" for x in keeps[:4]) + ".</font>",
                st["small"]))
    else:
        verb = "lower" if want == "high" else "raise"
        agents = m.for_direction(want) if m is not None else ()
        block.append(Paragraph(
            f"<b>What has been shown to {verb} it</b> &nbsp;<font size='7' color='{_hex(INK_FAINT)}'>"
            f"{len(agents)} item{'s' if len(agents) != 1 else ''}, strongest evidence first</font>",
            st["h3"]))
        if agents:
            block.append(_agents_table(st, agents, colour))
        else:
            block.append(Paragraph(
                f"Nothing in the literature we hold has been shown to {verb} this organism specifically. "
                "The general levers \u2014 fibre diversity, less processed food, avoiding unnecessary "
                "antibiotics \u2014 shift the whole community and are the honest fallback.",
                st["small"]))
        if m is not None and m.competitors and want == "high":
            names = ", ".join(f"<i>{c.replace('_', ' ')}</i>" for c in m.competitors)
            block.append(Paragraph(
                f"<b>Competes with:</b> {names}. Organisms that occupy the same niche or make compounds "
                "that inhibit it; raising them is an indirect way of lowering this one.",
                st["small"]))
        other = m.for_direction("low" if want == "high" else "high") if m is not None else ()
        if other:
            block.append(Paragraph(
                f"<font size='6.6' color='{_hex(INK_FAINT)}'>Known to push it the other way: "
                + "; ".join(f"{x.agent} ({x.evidence_words})" for x in other[:4]) + ".</font>",
                st["small"]))

    story.append(KeepTogether(block))

    # ---- 4. the registry's evidence cards --------------------------------- #
    if plan is not None:
        cards = plan.cards_for("species", a.species)
        if cards:
            evidence_cards(story, st, cards, heading="What the intervention registry holds for this organism")
    story.append(Spacer(1, 3 * mm))


def _agents_table(st: dict[str, ParagraphStyle], agents: tuple[mods.Agent, ...], colour: colors.Color) -> Table:
    ev_colour = {4: GREEN, 3: INK, 2: AMBER, 1: SLATE}
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("WHAT", "KIND", "EVIDENCE", "NOTE")]]
    for ag in agents:
        note = ag.note
        if ag.source:
            note += (f" <font color='{_hex(INK_FAINT)}' size='6'>{ag.source}</font>" if note
                     else f"<font color='{_hex(INK_FAINT)}' size='6'>{ag.source}</font>")
        rows.append([
            Paragraph(f"<b>{ag.agent}</b>", st["cell"]),
            Paragraph(f"<font size='6.6'>{ag.kind_words}</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(ev_colour.get(ag.evidence_rank, SLATE))}' size='6.8'>"
                      f"<b>{ag.evidence_words}</b></font>", st["cell"]),
            Paragraph(f"<font size='6.6'>{note}</font>" if note else "", st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.28, 0.11, 0.13, 0.48)]
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE), ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
        ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2.4), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
    ]))
    return t


# --------------------------------------------------------------------------- #
# actions-section opener: your levers, ranked
# --------------------------------------------------------------------------- #


@dataclass
class _Lever:
    family: str
    helps: dict[str, mods.Agent]      # display name -> the agent that said so
    hurts: dict[str, mods.Agent]
    kinds: list[str]

    @property
    def best_rank(self) -> int:
        return max((a.evidence_rank for a in self.helps.values()), default=0)

    @property
    def n_trials(self) -> int:
        return sum(1 for a in self.helps.values() if a.evidence == "human_trial")

    @property
    def kind(self) -> str:
        from collections import Counter
        return Counter(self.kinds).most_common(1)[0][0] if self.kinds else "compound"


def _aggregate(items: list[Attention]) -> tuple[list[_Lever], list[_Lever]]:
    """Every lever family across every flagged organism: who it helps, who it
    hurts. Returns (levers, exposures-to-avoid)."""
    fam: dict[str, _Lever] = {}
    avoid: dict[str, _Lever] = {}

    def slot(table: dict[str, _Lever], family: str) -> _Lever:
        if family not in table:
            table[family] = _Lever(family=family, helps={}, hurts={}, kinds=[])
        return table[family]

    for a in items:
        m = a.modulators
        want = a.want
        if m is None or want is None:
            continue
        stated_help = m.for_direction(want)
        stated_hurt = m.for_direction("low" if want == "high" else "high")
        for agent, helps_as_stated in [(x, True) for x in stated_help] + [(x, False) for x in stated_hurt]:
            family = mods.lever_family(agent)
            verdict = mods.family_helps(agent, helps_as_stated=helps_as_stated)
            if verdict is None:
                # Pure exposure: record only where it pushes the wrong way,
                # and only for readings that matter.
                if not helps_as_stated and a.counts_against:
                    lv = slot(avoid, family)
                    lv.hurts.setdefault(a.display, agent)
                    lv.kinds.append(agent.kind)
                continue
            if not verdict and not a.counts_against:
                continue
            lv = slot(fam, family)
            (lv.helps if verdict else lv.hurts).setdefault(a.display, agent)
            lv.kinds.append(agent.kind)
    levers = sorted(
        (lv for lv in fam.values() if lv.helps),
        key=lambda lv: (-len(lv.helps), -lv.best_rank, -lv.n_trials, lv.family),
    )
    # Levers with nothing on the helping side are exposures for this sample.
    for lv in fam.values():
        if not lv.helps and lv.hurts:
            avoid[lv.family] = lv
    exposures = sorted(avoid.values(), key=lambda lv: (-len(lv.hurts), lv.family))
    return levers, exposures


def levers(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    items: list[Attention],
    detail_section: int,
) -> None:
    """Actions-section opener: what would move the most of your flagged
    organisms the right way, and what pushes them the wrong way."""
    story.extend(section_heading(section, "What you can do about it", st["h1"]))
    lv, exposures = _aggregate(items)
    n_org = sum(1 for a in items if a.modulators is not None)
    story.append(Paragraph(
        f"Section {detail_section} lists, organism by organism, what has been shown to move each one. "
        f"This page adds them up. Across the <b>{n_org}</b> organisms flagged in your sample, each lever "
        "below is ranked by how many of them it moves the right way, with the strongest evidence for any "
        "of those organisms stated, and with the organisms it would push the <i>wrong</i> way named "
        "beside it rather than hidden. A lever that helps six organisms and hurts one is still a good "
        "lever; one that helps one and hurts four is not, and the table shows which is which. Evidence "
        "words mean what they say: a laboratory result is a lead, a human trial is a result.",
        st["body"]))

    if not lv:
        story.append(Paragraph("<i>No organism in this sample has a documented lever.</i>", st["body"]))
    else:
        story.append(Paragraph(
            f"<font color='{_hex(GREEN)}'><b>YOUR LEVERS, RANKED</b></font>"
            f"  <font size='7' color='{_hex(INK_FAINT)}'>{len(lv)}</font>", st["h2"]))
        story.append(_levers_table(st, lv))

    if exposures:
        story.append(Spacer(1, 3 * mm))
        story.append(KeepTogether([
            Paragraph(
                f"<font color='{_hex(CORAL)}'><b>WHAT PUSHES YOUR FLAGGED ORGANISMS THE WRONG WAY</b></font>"
                f"  <font size='7' color='{_hex(INK_FAINT)}'>{len(exposures)}</font>", st["h2"]),
            Paragraph(
                "Exposures with documented effects on organisms flagged in this sample, in the direction "
                "you would not want. Some are unavoidable (a needed antibiotic course); knowing which "
                "organisms they affect says what to rebuild afterwards.", st["small"]),
            _exposures_table(st, exposures),
        ]))

    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_FAINT)}'>Counts are organisms, not studies: a lever "
        "\u201chelps 5\u201d moves five of your flagged organisms the right way according to at least one "
        "cited study each. The evidence word is the strongest among those studies. Nothing here is a "
        "prescription; medication changes are a prescriber's decision and are marked as such.</font>",
        st["small"]))


def _names(d: dict[str, mods.Agent], *, colour: colors.Color, limit: int = 6) -> str:
    keys = sorted(d, key=lambda k: -d[k].evidence_rank)
    shown = ", ".join(f"<i>{k}</i>" for k in keys[:limit])
    if len(keys) > limit:
        shown += f" <font color='{_hex(INK_FAINT)}'>+{len(keys) - limit} more</font>"
    return f"<font color='{_hex(colour)}'>{shown}</font>" if shown else "\u2014"


def _levers_table(st: dict[str, ParagraphStyle], lv: list[_Lever]) -> Table:
    ev_colour = {4: GREEN, 3: INK, 2: AMBER, 1: SLATE}
    header = ["LEVER", "KIND", "HELPS", "ORGANISMS IT MOVES THE RIGHT WAY", "EVIDENCE", "WORKS AGAINST"]
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in header]]
    for x in lv:
        words = mods.EVIDENCE[
            next(k for k, (_w, r) in mods.EVIDENCE.items() if r == x.best_rank)][0] if x.best_rank else "\u2014"
        trials = f"<br/><font size='6' color='{_hex(INK_FAINT)}'>{x.n_trials} human trial{'s' if x.n_trials != 1 else ''}</font>" if x.n_trials else ""
        against = (f"<font color='{_hex(CORAL)}'><b>{len(x.hurts)}</b></font> "
                   f"<font size='6.4'>{_names(x.hurts, colour=CORAL, limit=3)}</font>") if x.hurts else \
                  f"<font color='{_hex(INK_FAINT)}'>none</font>"
        rows.append([
            Paragraph(f"<b>{x.family}</b>", st["cell"]),
            Paragraph(f"<font size='6.6'>{mods.KIND.get(x.kind, x.kind)}</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(GREEN)}' size='10'><b>{len(x.helps)}</b></font>", st["pct"]),
            Paragraph(f"<font size='6.6'>{_names(x.helps, colour=INK)}</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(ev_colour.get(x.best_rank, SLATE))}' size='6.8'><b>{words}</b></font>{trials}", st["cell"]),
            Paragraph(against, st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.24, 0.09, 0.06, 0.31, 0.11, 0.19)]
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE), ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PANEL_BG]),
    ]))
    return t


def _exposures_table(st: dict[str, ParagraphStyle], ex: list[_Lever]) -> Table:
    header = ["EXPOSURE", "KIND", "PUSHES THE WRONG WAY", "ORGANISMS AFFECTED"]
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in header]]
    for x in ex:
        rows.append([
            Paragraph(f"<b>{x.family}</b>", st["cell"]),
            Paragraph(f"<font size='6.6'>{mods.KIND.get(x.kind, x.kind)}</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(CORAL)}' size='10'><b>{len(x.hurts)}</b></font>", st["pct"]),
            Paragraph(f"<font size='6.6'>{_names(x.hurts, colour=INK, limit=8)}</font>", st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.30, 0.10, 0.12, 0.48)]
    t = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE), ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#FDF3F2")),
    ]))
    return t


__all__ = ["Attention", "GROUPS", "build", "detail", "glance", "levers"]


#: Reader-facing names for the lanes, for the provenance row.
_LANE_WORDS = {
    "scoring": "MetaPhlAn 3", "extended": "MetaPhlAn 4 (Jun23)", "jan26": "MetaPhlAn 4.2 (Jan26)",
    "genome": "GTDB whole-genome sketch", "globdb": "GlobDB whole-genome sketch", "motus": "mOTUs 4",
    "kraken": "Kraken2/UHGG", "singlem": "SingleM",
}


def _provenance_rows(a: Attention, *, bacdive: Mapping[str, Any] | None,
                     strain_analysis: Mapping[str, Any] | None) -> list[tuple[str, str]]:
    """Which methods saw it, how sure, and what the reference strains are like.

    Spec 0.8.4 §6: detection support, stable native identifiers, reference
    release and available interpretation on every organism's detail. A
    BacDive trait is a reference strain's trait and is labelled so.
    """
    o = getattr(a.verdict, "organism", None) if a.verdict is not None else None
    rows: list[tuple[str, str]] = []
    if o is None:
        return rows
    lanes = [_LANE_WORDS.get(x, x) for x in (o.lanes or ())]
    if lanes:
        rows.append(("DETECTED BY", ", ".join(lanes)))
    ids = _card_identifiers(o)
    if ids:
        rows.append(("IDENTIFIERS", ", ".join(ids) + (f"; GTDB R232 {o.gtdb}" if o.gtdb and not o.unnamed else "")))
    status = getattr(o, "status", "supported")
    basis = getattr(o, "confidence_basis", "")
    rows.append(("DETECTION", f"<b>{status}</b>" + (f" \u2014 {basis}" if basis else "")))
    if bacdive and isinstance(bacdive.get("species"), Mapping):
        rec = bacdive["species"].get(o.species) or bacdive["species"].get(o.species.replace(" ", "_"))
        if rec and rec.get("status") == "found":
            bits = []
            if rec.get("culture_collection_ids"):
                bits.append("culture collection " + ", ".join(rec["culture_collection_ids"][:3]))
            if rec.get("oxygen_requirement"):
                bits.append(", ".join(rec["oxygen_requirement"]))
            if rec.get("growth_temperature_c"):
                lo, hi = rec["growth_temperature_c"]
                bits.append(f"grows at {lo:g}\u2013{hi:g} \u00b0C" if lo != hi else f"grows at {lo:g} \u00b0C")
            if rec.get("isolation_sources"):
                bits.append("isolated from " + ", ".join(rec["isolation_sources"][:2]).lower())
            if rec.get("products"):
                bits.append("produces " + ", ".join(rec["products"][:4]))
            if bits:
                rows.append(("REFERENCE STRAINS", f"{rec['n_strains']} in BacDive: " + "; ".join(bits)
                             + ". Traits of reference strains, not measured in your sample."))
    if strain_analysis and strain_analysis.get("status") == "completed" and o.sgb:
        for c in strain_analysis.get("clades") or []:
            if c.get("sgb") == o.sgb and c.get("status") == "placed":
                near = str(c.get("nearest_reference") or "").split("|")[0]
                rows.append(("STRAIN PLACEMENT", f"nearest reference genome {near}, distance {c.get('distance_to_nearest')}; "
                             f"{c.get('mixture_evidence')}. A placement, not proof of an identical strain."))
                break
    return rows
