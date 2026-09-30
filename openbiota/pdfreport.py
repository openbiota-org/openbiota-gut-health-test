"""OpenBiota Gut Health Test — the consumer PDF (the Gut Health Test Metagenomic Report).

The document opens with the answer and works towards the evidence:

  1. Summary dashboard — every measurement, one page, graphics first
  2. What stood out, and how to read this report
  PART A — metabolite pathway gene capacity
    3. Metabolites at a glance — every pathway, one line each
    4. Metabolites in detail — one card per pathway
  PART B — community and disease-pattern resemblance
    5. Your gut community — species-level, with the general gut-health index
    6. Resemblance to published disease patterns — one page per pattern
  7. What this report cannot tell you
  8. How this test performs — measured validation figures for both parts
  9. Technical data

One scale, everywhere. Every percentile in the document — metabolite pathways,
disease patterns, module scores — is translated to the same seven-step
vocabulary by :func:`status_for`, so "high" means the same thing on every
page. Colour carries the judgement (does the literature read that direction
as favourable or adverse); the label carries the position.

Charts are reportlab vector flowables, so the file stays small and sharp.
"""

from __future__ import annotations

import datetime as dt
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    BaseDocTemplate,
    CondPageBreak,
    Flowable,
    Frame,
    Image,
    KeepTogether,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Spacer,
    Table,
    TableStyle,
)

from openbiota import __version__
from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    back_link_line,
    function_dest,
    linked_cell,
    section_heading,
)

#: The brand, and the report's own name. Kept apart: the brand is who made
#: the report, the title is what the report is, and the running header needs
#: both without repeating either.
BRAND: Final = "OpenBiota"
REPORT_TITLE: Final = "Gut Health Test Metagenomic Report"
#: Where the brand points. Every wordmark and every set brand name in the
#: report is a live link here, so a reader holding a printout of page 40 can
#: still find the project without a search engine.
BRAND_URL: Final = "https://openbiota.com"

#: Wordmark shipped with the package. Absent or unreadable, the title block
#: falls back to setting the brand as type, so a missing asset costs the mark
#: rather than the report.
LOGO_FILE: Final = "openbiota_logo.png"
#: Height the wordmark is drawn at. The asset is 300x67, so the width follows
#: from the aspect ratio rather than being asserted here; hard-coding both
#: would distort it if the asset is replaced.
LOGO_HEIGHT: Final = 9.5 * mm
#: The page-1 header band: from the foot of the accent bar to the top of the
#: first line of content (the kicker). The wordmark and the sample line are
#: both centred on this band's midline, so the header reads as one row set
#: comfortably in its space rather than crowded against the page edge.
HEADER_BAND: Final = (2.2 * mm, 21.0 * mm)
HEADER_BAND_MID: Final = (HEADER_BAND[0] + HEADER_BAND[1]) / 2.0
#: Distance from the page top to the top of the wordmark.
HEADER_LOGO_TOP: Final = HEADER_BAND_MID - LOGO_HEIGHT / 2.0
#: Baseline of the sample line on page 1: its cap height (0.72 em of 6.4 pt)
#: centred on the same midline.
HEADER_META_BASELINE: Final = HEADER_BAND_MID + 0.72 * 6.4 / 2.0


#: One place for the section numbers every cross-reference in the report uses.
#: Part A is the reading order the report is built around: who is there and
#: who is missing → how old the community reads → what it can make → what it
#: resembles → what can be done. Part B repeats the order in depth.
SECTIONS: Final[dict[str, int]] = {
    # The contents page carries no number - it is the map, not a stop on
    # the route - so the summary comes first and "what stood out" second.
    # Zero is how "no number" is spelled here: `section_heading` prints the
    # title alone for it, and the contents page lists it without a number.
    "summary": 1, "guide": 0, "stood_out": 2,
    "community": 3, "groups": 4, "organisms": 5,
    # The full organism list sits directly after the organisms that need
    # attention: having just read which ones are flagged, the next question
    # is always "what else is in there", and it should not require a jump to
    # the back of the report to answer.
    "catalogue": 6,
    "pathogens": 7,
    "age": 8, "functions": 9, "patterns": 10, "biofilm": 11, "mycobiome": 12,
    "skin": 13, "actions": 14,
    "groups_detail": 15, "organisms_detail": 16, "strains": 17,
    "pathogens_detail": 18,
    "functions_detail": 19, "patterns_detail": 20, "biofilm_detail": 21, "mycobiome_detail": 22,
    "skin_detail": 23,
    # Spec 0.8.3 §10.1 puts the input register before the closing four
    # sections. Anchors are semantic; only these numbers move.
    "context": 24,
    "limits": 25, "sequencing": 26, "accuracy": 27, "technical": 28,
}
TOOL_NAME: Final = f"{BRAND} Gut Health Test"

#: Named destination for page one. Defined by the page chrome rather than by
#: a section, so the running header's home link can never dangle.
HOME_DEST: Final = "home"

# --------------------------------------------------------------------------- #
# palette
# --------------------------------------------------------------------------- #

INK: Final = colors.HexColor("#15222E")
INK_SOFT: Final = colors.HexColor("#4B5B6B")
INK_FAINT: Final = colors.HexColor("#8C99A6")
ACCENT: Final = colors.HexColor("#0F8B8D")
ACCENT_DARK: Final = colors.HexColor("#0B6466")
ACCENT_SOFT: Final = colors.HexColor("#E3F2F2")
RULE: Final = colors.HexColor("#E1E7EC")
PANEL_BG: Final = colors.HexColor("#F5F8FA")
CARD_BG: Final = colors.HexColor("#FAFBFC")

GREEN: Final = colors.HexColor("#2E9C6A")
GREEN_BG: Final = colors.HexColor("#E4F4EA")
AMBER: Final = colors.HexColor("#D69E1C")
AMBER_BG: Final = colors.HexColor("#FBF1D6")
ORANGE: Final = colors.HexColor("#E07A2B")
ORANGE_BG: Final = colors.HexColor("#FBE8D8")
CORAL: Final = colors.HexColor("#C8453C")
CORAL_BG: Final = colors.HexColor("#F9E4E2")
SLATE: Final = colors.HexColor("#6A7988")
SLATE_BG: Final = colors.HexColor("#EDF1F4")

# Pastel track segments for the scale bar. Saturated colours above mark the
# reading; these mark the bands it can fall in.
_PASTEL: Final = {
    "good": colors.HexColor("#C9E6D4"),
    "typical": colors.HexColor("#DCEFE3"),
    "amber": colors.HexColor("#F5E4B4"),
    "orange": colors.HexColor("#F5D2B4"),
    "coral": colors.HexColor("#EFC3BF"),
    "slate": colors.HexColor("#DEE4EA"),
    "slate_typical": colors.HexColor("#E9EEF2"),
}

# Community composition palette, ordered by abundance rank.
PHYLA_PALETTE: Final = (
    ACCENT,
    colors.HexColor("#4E8FA8"),
    colors.HexColor("#7FB069"),
    colors.HexColor("#E0A93A"),
    colors.HexColor("#9A7AA0"),
    colors.HexColor("#C77B4C"),
    colors.HexColor("#8A98A5"),
)

PAGE: Final = A4
MARGIN: Final = 16 * mm
CONTENT_WIDTH: Final = PAGE[0] - 2 * MARGIN
HEADER_H: Final = 13 * mm
FOOTER_H: Final = 13 * mm
#: Usable height of the page-1 frame; the summary shrinks itself to fit it.
#: Page 1 has no running title, but it does carry the wordmark in the
#: header, so the first frame starts below the mark rather than under it.
FIRST_FRAME_HEIGHT: Final = PAGE[1] - FOOTER_H - 16.6 * mm
#: Usable height of an interior page, for content that must fit on one page.
BODY_FRAME_HEIGHT: Final = PAGE[1] - FOOTER_H - HEADER_H - 2 * mm


def _hex(colour: colors.Color) -> str:
    return "#" + colour.hexval()[2:]


def _lerp(a: colors.Color, b: colors.Color, t: float) -> colors.Color:
    t = max(0.0, min(1.0, t))
    return colors.Color(
        a.red + (b.red - a.red) * t,
        a.green + (b.green - a.green) * t,
        a.blue + (b.blue - a.blue) * t,
    )


def _ramp(stops: Sequence[tuple[float, colors.Color]], t: float) -> colors.Color:
    for (t0, c0), (t1, c1) in zip(stops, stops[1:], strict=False):
        if t0 <= t <= t1:
            return _lerp(c0, c1, (t - t0) / (t1 - t0) if t1 > t0 else 0.0)
    return stops[-1][1] if t > stops[-1][0] else stops[0][1]


# --------------------------------------------------------------------------- #
# the one scale
# --------------------------------------------------------------------------- #

#: Percentile boundaries of the seven bands, low to high.
SCALE_BOUNDS: Final = (0.0, 5.0, 15.0, 25.0, 75.0, 85.0, 95.0, 100.0)
#: Band level for each segment between consecutive bounds.
SCALE_LEVELS: Final = (-3, -2, -1, 0, 1, 2, 3)

LEVEL_LABELS: Final = {
    3: "notably high",
    2: "high",
    1: "above average",
    0: "typical",
    -1: "below average",
    -2: "low",
    -3: "notably low",
}

LEVEL_NOTES: Final = {
    3: "well above most people",
    2: "higher than most people",
    1: "somewhat higher than most people",
    0: "in line with most people",
    -1: "somewhat lower than most people",
    -2: "lower than most people",
    -3: "well below most people",
}


@dataclass(frozen=True, slots=True)
class Status:
    """A reading's headline status: position on the scale, read in the
    direction the literature reports."""

    label: str
    colour: colors.Color
    background: colors.Color
    note: str
    level: int = 0

    @property
    def assessed(self) -> bool:
        """Whether this reading was actually placed against the reference.

        Three states share the neutral colour for three different reasons:
        no reference range, nothing detected, and a provisional identity
        that withholds judgement. None of them is a typical result, but
        summaries that counted "everything not coral, orange or amber" as
        typical folded all three in - so a pathway at the 6th percentile
        could be reported as typical because its identity was provisional.
        """
        if self.label in {"no result", "not detected", "not scored", "present"}:
            return False
        if "too few reads for a stable reading" in self.note:
            return False
        return "below the confirmed identity threshold" not in self.note

    @property
    def unassessed_reason(self) -> str:
        """Why this reading was not placed, in words."""
        if self.assessed:
            return ""
        return {
            "no result": "no reference range",
            "not detected": "not detected",
            "not scored": "not scored",
            "present": "most of the cohort carries none, so a position would overstate it",
        }.get(self.label, "provisional identity")


def logo_path() -> Path | None:
    """Where the wordmark lives, or `None` if it is not installed."""
    try:
        from importlib import resources

        with resources.as_file(
            resources.files("openbiota.data").joinpath(LOGO_FILE)
        ) as path:
            return path if path.is_file() else None
    except (ModuleNotFoundError, FileNotFoundError, OSError):
        return None


def logo_flowable(height: float = LOGO_HEIGHT) -> Flowable | None:
    """The wordmark at `height`, width taken from the asset's own aspect.

    Returns `None` when the asset is missing or unreadable so callers can set
    the brand as type instead. A report that fails to render because a
    decorative PNG moved would be a poor trade.
    """
    path = logo_path()
    if path is None:
        return None
    try:
        reader = ImageReader(str(path))
        width_px, height_px = reader.getSize()
        if not height_px:
            return None
        return Image(str(path), width=height * (width_px / height_px), height=height)
    except Exception:  # noqa: BLE001 - never lose a report over the logo
        return None


def level_for(percentile: float) -> int:
    """Which of the seven bands a percentile falls in."""
    if percentile >= 95:
        return 3
    if percentile >= 85:
        return 2
    if percentile >= 75:
        return 1
    if percentile <= 5:
        return -3
    if percentile <= 15:
        return -2
    if percentile <= 25:
        return -1
    return 0


def _judgement(level: int, higher_means: str) -> tuple[colors.Color, colors.Color]:
    """Colour for a band, given which direction the literature calls adverse."""
    if higher_means == "adverse":
        adverse = level
    elif higher_means == "favourable":
        adverse = -level
    else:
        return SLATE, SLATE_BG
    if adverse >= 3:
        return CORAL, CORAL_BG
    if adverse == 2:
        return ORANGE, ORANGE_BG
    if adverse == 1:
        return AMBER, AMBER_BG
    return GREEN, GREEN_BG


def band_colour(level: int, higher_means: str) -> colors.Color:
    """Pastel track colour for one band of the scale bar."""
    if higher_means not in ("adverse", "favourable"):
        return _PASTEL["slate_typical"] if level == 0 else _PASTEL["slate"]
    adverse = level if higher_means == "adverse" else -level
    if adverse >= 3:
        return _PASTEL["coral"]
    if adverse == 2:
        return _PASTEL["orange"]
    if adverse == 1:
        return _PASTEL["amber"]
    if adverse == 0:
        return _PASTEL["typical"]
    return _PASTEL["good"]


def status_for(
    *,
    percentile: float | None,
    higher_means: str = "adverse",
    detected: bool = True,
    confidence: str = "confirmed",
    stable: bool = True,
    cohort_mostly_absent: bool = False,
) -> Status:
    """Translate a percentile into the report's one shared vocabulary.

    The label is pure position (typical, high, notably low …); the colour is
    the judgement — a high reading is coral when higher is the adverse
    direction and green when higher is the favourable one. Where the
    literature has no settled direction, colour stays neutral.
    """
    if percentile is None:
        return Status("no result", SLATE, SLATE_BG, "could not be measured in this sample")
    if not detected:
        return Status("not detected", SLATE, SLATE_BG, "no reads matched this pathway")
    if cohort_mostly_absent and percentile:
        # More than half the reference has none of this at all, so the
        # median is zero and *any* detection outranks it. The colorectal
        # virulence panel read "above average" on four fragments for that
        # reason alone, with none of the three organisms it looks for
        # present in the sample. The position is arithmetically true and
        # says nothing about the amount, so the word for it is presence.
        return Status(
            "present",
            SLATE,
            SLATE_BG,
            "more than half the reference cohort carries none of this, so any amount "
            "ranks above the middle; that is a statement about how uncommon this is "
            "and not about how much of it you have",
            0,
        )
    if not stable:
        # Above the panel's own minimum the count is a measurement; below it
        # Poisson noise dominates and the band would be an accident of a
        # handful of reads. The colorectal virulence panel read "above
        # average" off two fragments, which is the sort of sentence that
        # frightens somebody for no reason. The position is still shown,
        # because the reads are real, and the colour withholds the verdict.
        return Status(
            LEVEL_LABELS[level_for(percentile)],
            SLATE,
            SLATE_BG,
            "too few reads for a stable reading; Poisson noise dominates at this count",
            level_for(percentile),
        )
    level = level_for(percentile)
    if confidence == "provisional":
        # Below the identity threshold at which validation showed no false
        # positives. The word still says where the reading sits — same
        # vocabulary as everything else — but the colour withholds judgement
        # so it can never be mistaken for a finding, and the reason lives in
        # the fine print of the pathway's own section.
        return Status(
            LEVEL_LABELS[level],
            SLATE,
            SLATE_BG,
            "detected below the confirmed identity threshold; read as indicative only",
            level,
        )
    colour, background = _judgement(level, higher_means)
    return Status(LEVEL_LABELS[level], colour, background, LEVEL_NOTES[level], level)


scale_for = status_for

NOT_SCORED: Final = Status("not scored", SLATE, SLATE_BG, "no score was produced")


def notability(status: Status, percentile: float | None) -> tuple[int, float]:
    """Sort key: most notable first. Judged readings outrank neutral ones."""
    judged = status.colour in (CORAL, ORANGE, AMBER)
    return (
        0 if judged else 1,
        -abs(status.level) * 100 - abs((percentile or 50.0) - 50.0),
    )


DIRECTION_WORDS: Final = {
    "adverse": "Research associates higher levels with worse outcomes",
    "favourable": "Research associates higher levels with better outcomes",
    "context-dependent": "The research picture is mixed",
    "unclear": "No clear direction reported",
}

STRENGTH_WORDS: Final = {
    "strong": "Well established",
    "moderate": "Reasonably supported",
    "limited": "Early evidence",
    "preliminary": "Preliminary",
}


# --------------------------------------------------------------------------- #
# flowables
# --------------------------------------------------------------------------- #


class PercentileBar(Flowable):
    """The scale bar: seven pastel bands, a median tick, and your marker.

    The same drawing is used for every percentile in the report so the eye
    learns it once.
    """

    def __init__(
        self,
        *,
        width: float,
        percentile: float | None,
        higher_means: str = "adverse",
        marker_colour: colors.Color | None = None,
        height: float = 6.0 * mm,
        show_scale: bool = False,
        track: float | None = None,
    ) -> None:
        super().__init__()
        self.width = width
        self.height = height + (3.2 * mm if show_scale else 0)
        self.percentile = percentile
        self.higher_means = higher_means
        self.marker_colour = marker_colour
        self.show_scale = show_scale
        self.track = track or min(3.4 * mm, height * 0.56)

    def wrap(self, *_args: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        w, t = self.width, self.track
        y = (self.height - (3.2 * mm if self.show_scale else 0) - t) / 2.0 + (
            3.2 * mm if self.show_scale else 0
        )

        c.saveState()
        clip = c.beginPath()
        clip.roundRect(0, y, w, t, t / 2.0)
        c.clipPath(clip, stroke=0, fill=0)
        for index, level in enumerate(SCALE_LEVELS):
            x0 = w * SCALE_BOUNDS[index] / 100.0
            x1 = w * SCALE_BOUNDS[index + 1] / 100.0
            c.setFillColor(band_colour(level, self.higher_means))
            c.rect(x0, y, x1 - x0 + 0.3, t, stroke=0, fill=1)
        c.restoreState()

        # median
        c.setStrokeColor(colors.white)
        c.setLineWidth(1.0)
        c.line(w / 2.0, y, w / 2.0, y + t)

        if self.percentile is not None:
            r = t * 0.78
            x = w * max(0.0, min(1.0, self.percentile / 100.0))
            x = min(max(x, r), w - r)
            colour = self.marker_colour or status_for(
                percentile=self.percentile, higher_means=self.higher_means
            ).colour
            c.setFillColor(colour)
            c.setStrokeColor(colors.white)
            c.setLineWidth(1.3)
            c.circle(x, y + t / 2.0, r, stroke=1, fill=1)
        else:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 5.6)
            c.drawCentredString(w / 2.0, y + t + 1.2, "not measured")

        if self.show_scale:
            c.setFont("Helvetica", 5.8)
            c.setFillColor(INK_FAINT)
            yy = y - 2.6 * mm
            c.drawString(0, yy, "lower")
            c.drawCentredString(w * 0.25, yy, "25th")
            c.drawCentredString(w * 0.50, yy, "median")
            c.drawCentredString(w * 0.75, yy, "75th")
            c.drawRightString(w, yy, "higher")


class ScaleLegend(Flowable):
    """The scale bar with every band named — the key to the whole report."""

    def __init__(self, *, width: float, higher_means: str = "adverse") -> None:
        super().__init__()
        self.width = width
        self.height = 19 * mm
        self.higher_means = higher_means

    def wrap(self, *_args: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        w = self.width
        t = 4.2 * mm
        y = self.height - t - 2.0 * mm
        c.saveState()
        clip = c.beginPath()
        clip.roundRect(0, y, w, t, t / 2.0)
        c.clipPath(clip, stroke=0, fill=0)
        for index, level in enumerate(SCALE_LEVELS):
            x0 = w * SCALE_BOUNDS[index] / 100.0
            x1 = w * SCALE_BOUNDS[index + 1] / 100.0
            c.setFillColor(band_colour(level, self.higher_means))
            c.rect(x0, y, x1 - x0 + 0.3, t, stroke=0, fill=1)
        c.restoreState()
        c.setStrokeColor(colors.white)
        c.setLineWidth(1.0)
        for bound in SCALE_BOUNDS[1:-1]:
            x = w * bound / 100.0
            c.line(x, y, x, y + t)

        c.setFont("Helvetica", 5.6)
        c.setFillColor(INK_FAINT)
        for bound in SCALE_BOUNDS[1:-1]:
            c.drawCentredString(w * bound / 100.0, y - 2.4 * mm, f"{bound:.0f}")
        c.drawString(0, y - 2.4 * mm, "0")
        c.drawRightString(w, y - 2.4 * mm, "100")
        c.drawCentredString(w * 0.5, y - 5.6 * mm - 8.0, "percentile of the reference group")

        for index, level in enumerate(SCALE_LEVELS):
            x0 = w * SCALE_BOUNDS[index] / 100.0
            x1 = w * SCALE_BOUNDS[index + 1] / 100.0
            colour, _ = _judgement(level, self.higher_means)
            c.setFillColor(colour if self.higher_means in ("adverse", "favourable") else INK_SOFT)
            c.setFont("Helvetica-Bold", 5.6)
            words = LEVEL_LABELS[level].split()
            if len(words) == 2 and (x1 - x0) < 18 * mm:
                c.drawCentredString((x0 + x1) / 2, y - 5.6 * mm, words[0].upper())
                c.drawCentredString((x0 + x1) / 2, y - 8.0 * mm, words[1].upper())
            else:
                c.drawCentredString((x0 + x1) / 2, y - 5.6 * mm, LEVEL_LABELS[level].upper())


class StatusChip(Flowable):
    """A rounded status pill in the scale's colours."""

    def __init__(
        self, status: Status, *, width: float = 26 * mm, height: float = 5.2 * mm,
        font_size: float = 6.3,
    ) -> None:
        super().__init__()
        self.status = status
        self.width = width
        self.height = height
        self.font_size = font_size

    def wrap(self, *_args: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        c.setFillColor(self.status.background)
        c.roundRect(0, 0, self.width, self.height, self.height / 2, stroke=0, fill=1)
        c.setFillColor(self.status.colour)
        c.setFont("Helvetica-Bold", self.font_size)
        c.drawCentredString(
            self.width / 2, self.height / 2 - self.font_size * 0.36, self.status.label.upper()
        )


class HBar(Flowable):
    """Simple horizontal proportion bar, used for composition lists."""

    def __init__(
        self,
        *,
        width: float,
        fraction: float,
        colour: colors.Color,
        height: float = 3.2 * mm,
    ) -> None:
        super().__init__()
        self.width = width
        self.height = height
        self.fraction = max(0.0, min(1.0, fraction))
        self.colour = colour

    def wrap(self, *_args: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        c.setFillColor(PANEL_BG)
        c.roundRect(0, 0, self.width, self.height, self.height / 2, stroke=0, fill=1)
        if self.fraction > 0:
            c.setFillColor(self.colour)
            c.roundRect(
                0, 0, max(self.width * self.fraction, self.height), self.height,
                self.height / 2, stroke=0, fill=1,
            )


class Rule(Flowable):
    def __init__(self, width: float, colour: colors.Color = RULE, thickness: float = 0.6) -> None:
        super().__init__()
        self.width = width
        self.height = thickness
        self.colour = colour
        self.thickness = thickness

    def wrap(self, *_args: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        self.canv.setStrokeColor(self.colour)
        self.canv.setLineWidth(self.thickness)
        self.canv.line(0, 0, self.width, 0)


_DIAL_STOPS: Final = (
    (0.00, CORAL),
    (0.30, CORAL),
    (0.42, ORANGE),
    (0.50, AMBER),
    (0.58, colors.HexColor("#8FBF6A")),
    (0.70, GREEN),
    (1.00, colors.HexColor("#1F7F55")),
)


class Dial(Flowable):
    """A half-circle meter with a needle.

    Used for the general gut-health index. It is a meter, so it behaves like
    one: the whole spectrum is drawn at full strength from red through amber
    to green, with a scale of ticks across it and its end values labelled,
    and the reading is shown as a *position* on that spectrum — a slim needle
    from the hub and a white-ringed marker on the sweep. Nothing is dimmed:
    dimming the sweep past the reading made a middling score look like "all
    you have is the red part", which is the opposite of what a meter says.

    GMWI2 runs roughly -6..+6 and is clipped to -3..+3 here, which is where
    all the interpretive weight is.
    """

    def __init__(
        self,
        *,
        value: float | None,
        low: float = -3.0,
        high: float = 3.0,
        width: float = 56 * mm,
        label: str = "",
        sublabel: str = "",
        caption: str = "",
    ) -> None:
        super().__init__()
        self.value = value
        self.low = low
        self.high = high
        self.width = width
        self.stroke = width * 0.10
        self.radius = width / 2.0 - self.stroke / 2.0 - 1
        self.text_space = 30.0 + (8.0 if caption else 0.0)
        self.height = self.radius + self.stroke / 2.0 + self.text_space + 2
        self.label = label
        self.sublabel = sublabel
        self.caption = caption

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def _fraction(self, value: float) -> float:
        return (min(max(value, self.low), self.high) - self.low) / (self.high - self.low)

    def _arc(self, c: Canvas, cx: float, cy: float, t0: float, t1: float, colour: colors.Color) -> None:
        c.setStrokeColor(colour)
        path = c.beginPath()
        path.arc(
            cx - self.radius, cy - self.radius, cx + self.radius, cy + self.radius,
            180.0 - 180.0 * t0, -180.0 * (t1 - t0),
        )
        c.drawPath(path, stroke=1, fill=0)

    def _radial(self, c: Canvas, cx: float, cy: float, t: float, r0: float, r1: float) -> None:
        a = math.radians(180.0 - 180.0 * t)
        c.line(cx + math.cos(a) * r0, cy + math.sin(a) * r0, cx + math.cos(a) * r1, cy + math.sin(a) * r1)

    def draw(self) -> None:  # noqa: PLR0915 - one contiguous drawing
        c: Canvas = self.canv
        radius = self.radius
        cx, cy = self.width / 2.0, self.text_space
        segments = 90
        overlap = 0.4 / 180.0  # hide the hairline seams between segments

        # --- the spectrum, full strength end to end ------------------------ #
        c.saveState()
        c.setLineWidth(self.stroke)
        c.setLineCap(0)
        for i in range(segments):
            t0, t1 = i / segments, (i + 1) / segments
            self._arc(c, cx, cy, t0, min(1.0, t1 + overlap), _ramp(_DIAL_STOPS, (t0 + t1) / 2))
        c.setLineCap(1)
        self._arc(c, cx, cy, 0.0, 1e-4, _DIAL_STOPS[0][1])
        self._arc(c, cx, cy, 1.0 - 1e-4, 1.0, _DIAL_STOPS[-1][1])
        c.restoreState()

        # --- the scale: ticks across the sweep, longer at zero ------------ #
        c.saveState()
        c.setStrokeColor(colors.white)
        half = self.stroke / 2.0
        step = 1.0 if self.high - self.low <= 8 else 2.0
        v = math.ceil(self.low / step) * step
        while v <= self.high + 1e-9:
            if self.low < v < self.high:
                t = self._fraction(v)
                if abs(v) < 1e-9:
                    c.setLineWidth(1.4)
                    self._radial(c, cx, cy, t, radius - half - 0.6, radius + half + 0.6)
                else:
                    c.setLineWidth(0.8)
                    self._radial(c, cx, cy, t, radius - half + 1.2, radius + half - 1.2)
            v += step
        c.restoreState()
        # End values, so the reading has a scale to be read against.
        c.setFillColor(INK_FAINT)
        c.setFont("Helvetica", 6.0)
        c.drawCentredString(cx - radius, cy - half - 7.0, f"{self.low:+.0f}")
        c.drawCentredString(cx + radius, cy - half - 7.0, f"{self.high:+.0f}")

        # --- the reading: needle from the hub, marker on the sweep -------- #
        if self.value is not None:
            f = self._fraction(self.value)
            angle = math.radians(180.0 - 180.0 * f)
            colour = _ramp(_DIAL_STOPS, f)
            tip = radius - half - 1.0
            base = 2.4
            c.saveState()
            c.setFillColor(INK)
            needle = c.beginPath()
            needle.moveTo(cx + math.cos(angle) * tip, cy + math.sin(angle) * tip)
            needle.lineTo(cx + math.cos(angle + math.pi / 2) * base, cy + math.sin(angle + math.pi / 2) * base)
            needle.lineTo(cx - math.cos(angle) * 3.2, cy - math.sin(angle) * 3.2)
            needle.lineTo(cx + math.cos(angle - math.pi / 2) * base, cy + math.sin(angle - math.pi / 2) * base)
            needle.close()
            c.drawPath(needle, stroke=0, fill=1)
            c.circle(cx, cy, 3.8, stroke=0, fill=1)
            c.setFillColor(colors.white)
            c.circle(cx, cy, 1.5, stroke=0, fill=1)
            # Marker: white ring on the sweep at the reading, centre in the
            # sweep's own colour there.
            mx, my = cx + math.cos(angle) * radius, cy + math.sin(angle) * radius
            c.setFillColor(colors.white)
            c.setStrokeColor(INK)
            c.setLineWidth(1.1)
            c.circle(mx, my, self.stroke * 0.44, stroke=1, fill=1)
            c.setFillColor(colour)
            c.circle(mx, my, self.stroke * 0.19, stroke=0, fill=1)
            c.restoreState()

        # --- the words ------------------------------------------------------ #
        c.setFillColor(INK)
        c.setFont("Helvetica-Bold", 17)
        c.drawCentredString(cx, cy - 18, "\u2014" if self.value is None else f"{self.value:+.2f}")
        if self.label:
            c.setFillColor(self._band_colour())
            c.setFont("Helvetica-Bold", 6.8)
            c.drawCentredString(cx, cy - 26, self.label.upper())
        if self.sublabel:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 5.9)
            c.drawCentredString(cx, cy - 33, self.sublabel)
        if self.caption:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 5.9)
            c.drawCentredString(cx, cy - 33 - 8, self.caption)

    def _band_colour(self) -> colors.Color:
        if self.value is None:
            return SLATE
        if self.value <= -1.0:
            return CORAL
        if self.value <= -0.25:
            return ORANGE
        if self.value < 0.25:
            return AMBER
        return GREEN


def _tint(colour: colors.Color, amount: float) -> colors.Color:
    """Blend a colour toward white. Opaque, so it prints and renders everywhere."""
    return colors.Color(
        colour.red + (1 - colour.red) * amount,
        colour.green + (1 - colour.green) * amount,
        colour.blue + (1 - colour.blue) * amount,
    )


class Donut(Flowable):
    """A donut chart with a headline number in the hole."""

    def __init__(
        self,
        *,
        parts: Sequence[tuple[str, float, colors.Color]],
        size: float = 36 * mm,
        centre_value: str = "",
        centre_label: str = "",
        thickness: float | None = None,
    ) -> None:
        super().__init__()
        self.parts = parts
        self.size = size
        self.width = size
        self.height = size
        self.centre_value = centre_value
        self.centre_label = centre_label
        self.thickness = thickness or size * 0.19

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        r = self.size / 2.0
        cx = cy = r
        total = sum(max(0.0, v) for _, v, _ in self.parts) or 1.0
        start = 90.0
        c.saveState()
        for _, value, colour in self.parts:
            extent = -360.0 * max(0.0, value) / total
            if extent == 0:
                continue
            c.setFillColor(colour)
            c.setStrokeColor(colors.white)
            c.setLineWidth(0.9)
            c.wedge(cx - r, cy - r, cx + r, cy + r, start, extent, stroke=1, fill=1)
            start += extent
        inner = r - self.thickness
        c.setFillColor(colors.white)
        c.circle(cx, cy, inner, stroke=0, fill=1)
        c.restoreState()
        if self.centre_value:
            c.setFillColor(INK)
            c.setFont("Helvetica-Bold", min(15.0, inner * 0.62))
            c.drawCentredString(cx, cy - 2.5, self.centre_value)
        if self.centre_label:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 5.6)
            c.drawCentredString(cx, cy - 10.5, self.centre_label.upper())


def _fit_text(text: str, font: str, size: float, max_width: float) -> str:
    """Trim ``text`` with an ellipsis so it draws inside ``max_width``.

    Card footnotes are single-line by design; a note that runs past the card
    edge is worse than one that stops short, so it stops short visibly."""
    if pdfmetrics.stringWidth(text, font, size) <= max_width:
        return text
    ell = "\u2026"
    lo, hi = 0, len(text)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if pdfmetrics.stringWidth(text[:mid].rstrip() + ell, font, size) <= max_width:
            lo = mid
        else:
            hi = mid - 1
    return text[:lo].rstrip() + ell


class Tile(Flowable):
    """A key-figure card: big number, small caps label, optional footnote."""

    def __init__(
        self,
        *,
        value: str,
        label: str,
        note: str = "",
        width: float,
        height: float = 15.5 * mm,
        accent: colors.Color = ACCENT,
        value_colour: colors.Color = INK,
        big: bool = False,
        centred: bool = False,
        value_size: float | None = None,
    ) -> None:
        super().__init__()
        self.value_size = value_size
        self.value = value
        self.label = label
        self.note = note
        self.width = width
        self.height = height
        self.accent = accent
        self.value_colour = value_colour
        self.big = big
        self.centred = centred

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c: Canvas = self.canv
        c.setFillColor(CARD_BG)
        c.setStrokeColor(RULE)
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, 1.6 * mm, stroke=1, fill=1)
        c.setFillColor(self.accent)
        c.roundRect(0, 0, 1.4 * mm, self.height, 0.7 * mm, stroke=0, fill=1)
        c.rect(0.7 * mm, 0, 0.7 * mm, self.height, stroke=0, fill=1)
        x = 4.6 * mm
        c.setFillColor(self.value_colour)
        size = 15.5 if len(self.value) <= 6 else 12.5
        if self.big:
            size = 26.0 if len(self.value) <= 3 else 18.0
        if self.value_size is not None:
            size = self.value_size
        # Stack from the bottom: note, label, then the value in whatever room
        # is left, shrinking the value rather than letting it collide.
        note_y = 2.4 * mm
        label_y = note_y + 6.6 if self.note else 3.2 * mm
        label_top = label_y + 4.6
        value_y = label_top + 1.8
        room = self.height - 2.2 - value_y
        size = min(size, room / 0.72)
        if self.centred:
            # A tall tile with its figure sat in the bottom corner looks like
            # a mistake. Measure the stack and lift it to the middle; the text
            # stays left-aligned to the accent bar, which is where the eye
            # enters the card.
            stack_h = (value_y + size * 0.72) - (note_y if self.note else label_y)
            lift = (self.height - stack_h) / 2.0 - (note_y if self.note else label_y)
            note_y, label_y, value_y = note_y + lift, label_y + lift, value_y + lift
        c.setFont("Helvetica-Bold", size)
        c.drawString(x, value_y, self.value)
        c.setFillColor(INK_SOFT)
        c.setFont("Helvetica-Bold", 5.9)
        c.drawString(x, label_y, self.label.upper())
        if self.note:
            c.setFillColor(INK_FAINT)
            c.setFont("Helvetica", 5.7)
            c.drawString(x, note_y, _fit_text(self.note, "Helvetica", 5.7, self.width - x - 2 * mm))


class Dot(Flowable):
    """A coloured bullet for headline lists."""

    def __init__(self, colour: colors.Color, *, size: float = 2.2 * mm) -> None:
        super().__init__()
        self.colour = colour
        self.width = size
        self.height = size

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        self.canv.setFillColor(self.colour)
        self.canv.circle(self.width / 2, self.height / 2, self.width / 2, stroke=0, fill=1)


# --------------------------------------------------------------------------- #
# styles
# --------------------------------------------------------------------------- #


def _styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()["BodyText"]

    def make(name: str, **kw: Any) -> ParagraphStyle:
        return ParagraphStyle(name, parent=base, **kw)

    return {
        "title": make(
            "title", fontName="Helvetica-Bold", fontSize=24, leading=27,
            textColor=INK, spaceAfter=1.5 * mm,
        ),
        "tagline": make(
            "tagline", fontName="Helvetica", fontSize=9.6, leading=13.6,
            textColor=INK_SOFT,
        ),
        "h1": make(
            "h1", fontName="Helvetica-Bold", fontSize=17, leading=21,
            textColor=INK, spaceBefore=0, spaceAfter=1.5 * mm,
        ),
        "h2": make(
            "h2", fontName="Helvetica-Bold", fontSize=11.2, leading=14.5,
            textColor=ACCENT_DARK, spaceBefore=3.5 * mm, spaceAfter=1.2 * mm,
        ),
        "h2_tight": make(
            "h2_tight", fontName="Helvetica-Bold", fontSize=10.6, leading=13,
            textColor=ACCENT_DARK, spaceBefore=0, spaceAfter=0.8 * mm,
        ),
        "h3": make(
            "h3", fontName="Helvetica-Bold", fontSize=9.6, leading=13,
            textColor=INK, spaceAfter=0.8 * mm,
        ),
        "body": make(
            "body", fontName="Helvetica", fontSize=9.0, leading=13.2,
            textColor=INK_SOFT, spaceAfter=2 * mm, alignment=TA_LEFT,
        ),
        "body_ink": make(
            "body_ink", fontName="Helvetica", fontSize=9.0, leading=13.2,
            textColor=INK, spaceAfter=2 * mm,
        ),
        "lead": make(
            "lead", fontName="Helvetica", fontSize=9.6, leading=14,
            textColor=INK, spaceAfter=2 * mm,
        ),
        "small": make(
            "small", fontName="Helvetica", fontSize=7.4, leading=10.2,
            textColor=INK_FAINT, spaceAfter=1.2 * mm,
        ),
        "fine": make(
            "fine", fontName="Helvetica", fontSize=6.6, leading=9.0,
            textColor=INK_FAINT, spaceAfter=0.8 * mm,
        ),
        "fine_right": make(
            "fine_right", fontName="Helvetica", fontSize=6.6, leading=9.0,
            textColor=INK_FAINT, alignment=TA_RIGHT,
        ),
        "label": make(
            "label", fontName="Helvetica-Bold", fontSize=6.8, leading=9,
            textColor=INK_FAINT,
        ),
        "metric": make(
            "metric", fontName="Helvetica-Bold", fontSize=15, leading=17,
            textColor=INK,
        ),
        "kicker": make(
            "kicker", fontName="Helvetica-Bold", fontSize=7.4, leading=10,
            textColor=ACCENT, spaceAfter=1 * mm,
        ),
        "cell": make("cell", fontName="Helvetica", fontSize=7.8, leading=10.4, textColor=INK),
        "cell_soft": make(
            "cell_soft", fontName="Helvetica", fontSize=7.4, leading=9.8, textColor=INK_SOFT
        ),
        "cell_bold": make(
            "cell_bold", fontName="Helvetica-Bold", fontSize=7.8, leading=10.4, textColor=INK
        ),
        "pct": make(
            "pct", fontName="Helvetica-Bold", fontSize=8.6, leading=10.4, textColor=INK,
            alignment=TA_RIGHT,
        ),
        "glance_note": make(
            "glance_note", fontName="Helvetica", fontSize=7.4, leading=10.0,
            textColor=INK_SOFT, spaceAfter=0, leftIndent=1.5 * mm,
        ),
        "centre": make(
            "centre", fontName="Helvetica", fontSize=9, leading=13,
            textColor=INK_SOFT, alignment=TA_CENTER,
        ),
        # A shade tighter than it was (8.6/12), which is what pays for the
        # eighth row on page one - the one that carries the overgrowth.
        "headline": make(
            "headline", fontName="Helvetica", fontSize=8.3, leading=11.0,
            textColor=INK, spaceAfter=0,
        ),
        "divider": make(
            "divider", fontName="Helvetica-Bold", fontSize=30, leading=34,
            textColor=INK, spaceAfter=1 * mm,
        ),
        "h3_accent": make(
            "h3_accent", fontName="Helvetica-Bold", fontSize=9.6, leading=13,
            textColor=ACCENT_DARK, spaceBefore=1.5 * mm, spaceAfter=1 * mm,
        ),
        "evidence": make(
            "evidence", fontName="Helvetica", fontSize=7.8, leading=10.8,
            textColor=INK_SOFT, spaceAfter=1.2 * mm, leftIndent=3 * mm, firstLineIndent=-3 * mm,
        ),
        "small_ink": make(
            "small_ink", fontName="Helvetica", fontSize=7.6, leading=10.4,
            textColor=INK, spaceAfter=1.2 * mm,
        ),
    }


# --------------------------------------------------------------------------- #
# page furniture
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class _Chrome:
    sample: str
    generated: str
    meta_line: str

    def _fine_print(self, canvas: Canvas, baseline_from_top: float = 8.6 * mm) -> None:
        canvas.setFont("Helvetica", 6.4)
        canvas.setFillColor(INK_FAINT)
        canvas.drawRightString(PAGE[0] - MARGIN, PAGE[1] - baseline_from_top, self.meta_line)

    @staticmethod
    def _brand_text(
        canvas: Canvas, x: float, baseline: float, text: str, font: str, size: float
    ) -> float:
        """Set `text` at (x, baseline) and make the brand at its head a link.

        Only the leading `BRAND` (in whatever case `text` uses) becomes the
        hotspot, so the footer sentence and the running header stay plain
        text everywhere except on the word that actually names the project.
        Returns the width of the whole string so callers can continue the
        line after it.
        """
        canvas.setFont(font, size)
        canvas.drawString(x, baseline, text)
        brand = BRAND.upper() if text.startswith(BRAND.upper()) else BRAND
        if text.startswith(brand):
            width = pdfmetrics.stringWidth(brand, font, size)
            canvas.linkURL(
                BRAND_URL,
                (x, baseline - 0.25 * size, x + width, baseline + 0.78 * size),
                relative=0, thickness=0,
            )
        return pdfmetrics.stringWidth(text, font, size)

    def _wordmark(self, canvas: Canvas) -> bool:
        """Draw the brand at the top left of the header. True if it landed.

        Placed on the canvas rather than in the flow so it sits on the same
        line as the sample details, opposite them, and costs the page no
        vertical space at all.
        """
        path = logo_path()
        if path is None:
            return False
        try:
            reader = ImageReader(str(path))
            w_px, h_px = reader.getSize()
            if not h_px:
                return False
            height = LOGO_HEIGHT
            width = height * (w_px / h_px)
            bottom = PAGE[1] - HEADER_LOGO_TOP - height
            canvas.drawImage(
                reader, MARGIN, bottom, width=width, height=height, mask="auto",
            )
            canvas.linkURL(
                BRAND_URL, (MARGIN, bottom, MARGIN + width, bottom + height),
                relative=0, thickness=0,
            )
            return True
        except Exception:  # noqa: BLE001 - the logo never costs a report
            return False

    def first(self, canvas: Canvas, _doc: Any) -> None:
        canvas.saveState()
        # The home destination, defined by the chrome itself rather than by
        # any one section. The running header on every later page links here,
        # and ReportLab refuses to save a document with a dead link, so the
        # target has to be guaranteed: page one always draws this chrome, and
        # a report assembled without the summary page would otherwise fail.
        canvas.bookmarkPage(HOME_DEST)
        canvas.setFillColor(ACCENT)
        canvas.rect(0, PAGE[1] - 2.2 * mm, PAGE[0], 2.2 * mm, stroke=0, fill=1)
        if not self._wordmark(canvas):
            canvas.setFillColor(ACCENT_DARK)
            self._brand_text(
                canvas, MARGIN, PAGE[1] - HEADER_META_BASELINE, BRAND, "Helvetica-Bold", 9.0
            )
        self._fine_print(canvas, HEADER_META_BASELINE)
        self._footer(canvas, _doc)
        canvas.restoreState()

    def interior(self, canvas: Canvas, doc: Any) -> None:
        canvas.saveState()
        canvas.setFillColor(ACCENT)
        canvas.rect(0, PAGE[1] - 2.2 * mm, PAGE[0], 2.2 * mm, stroke=0, fill=1)
        # Running header: brand then title, the brand set bold so the two
        # read as one lockup rather than a repeated sentence.
        canvas.setFillColor(ACCENT_DARK)
        baseline = PAGE[1] - 8.6 * mm
        # The brand in the running header is plain text, not a website link:
        # the whole header run - brand, house, arrow, title - is one button
        # back to page one, so a reader who learns to click here is never
        # sent to two different places from the same spot. The logo on page
        # one and the footer still link to the website.
        canvas.setFont("Helvetica-Bold", 6.6)
        canvas.drawString(MARGIN, baseline, BRAND.upper())
        offset = pdfmetrics.stringWidth(BRAND.upper(), "Helvetica-Bold", 6.6)
        # Home affordance: a drawn house and an up arrow, then the title. The
        # base-14 fonts have no house glyph, so it is drawn rather than set.
        home_x = MARGIN + offset + 1.8 * mm
        icon_w = self._home_icon(canvas, home_x, baseline, 2.5 * mm)
        canvas.setFont("Helvetica", 6.6)
        arrow_x = home_x + icon_w + 1.0 * mm
        canvas.drawString(arrow_x, baseline, "\u2191")
        arrow_w = pdfmetrics.stringWidth("\u2191", "Helvetica", 6.6)
        title_x = arrow_x + arrow_w + 1.2 * mm
        canvas.drawString(title_x, baseline, REPORT_TITLE.upper())
        title_w = pdfmetrics.stringWidth(REPORT_TITLE.upper(), "Helvetica", 6.6)
        # One hotspot from the brand through house, arrow and title: the
        # reader should not have to find which part of it is the button.
        canvas.linkRect(
            "", HOME_DEST,
            (MARGIN - 0.6 * mm, baseline - 1.1 * mm,
             title_x + title_w + 0.6 * mm, baseline + 5.4),
            relative=0, thickness=0,
        )
        self._fine_print(canvas)
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(MARGIN, PAGE[1] - 11.2 * mm, PAGE[0] - MARGIN, PAGE[1] - 11.2 * mm)
        self._footer(canvas, doc)
        canvas.restoreState()

    @staticmethod
    def _home_icon(canvas: Canvas, x: float, y: float, size: float) -> float:
        """Draw a small house sitting on the text baseline.

        None of the base-14 fonts carries a house glyph — U+2302 and U+1F3E0
        both fall back to a .notdef box — so it is drawn as paths, the same
        way the dials and donuts elsewhere in the report are.

        Returns the width consumed, so the caller can lay out after it.
        """
        w = size
        h = size * 0.92
        eaves = y + h * 0.52
        canvas.saveState()
        canvas.setFillColor(ACCENT_DARK)
        canvas.setStrokeColor(ACCENT_DARK)
        canvas.setLineWidth(0.35)
        roof = canvas.beginPath()
        roof.moveTo(x - w * 0.06, eaves)
        roof.lineTo(x + w * 0.5, y + h)
        roof.lineTo(x + w * 1.06, eaves)
        roof.close()
        canvas.drawPath(roof, stroke=0, fill=1)
        # Walls, with a doorway knocked out so the shape reads as a house
        # rather than a solid block at 6.6pt.
        canvas.rect(x + w * 0.14, y, w * 0.72, h * 0.54, stroke=0, fill=1)
        canvas.setFillColor(colors.white)
        canvas.rect(x + w * 0.40, y, w * 0.21, h * 0.32, stroke=0, fill=1)
        canvas.restoreState()
        return w

    def _footer(self, canvas: Canvas, doc: Any) -> None:
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(MARGIN, 11 * mm, PAGE[0] - MARGIN, 11 * mm)
        canvas.setFillColor(INK_FAINT)
        self._brand_text(
            canvas, MARGIN, 7.4 * mm,
            f"{TOOL_NAME} v{__version__} · research use only · not a diagnostic test · "
            "every number in this document is also in results.json",
            "Helvetica", 6.2,
        )
        canvas.setFont("Helvetica", 6.2)
        canvas.drawRightString(PAGE[0] - MARGIN, 7.4 * mm, f"Page {doc.page}")


def _document(path: Path, chrome: _Chrome) -> BaseDocTemplate:
    doc = BaseDocTemplate(
        str(path),
        pagesize=PAGE,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN,
        title=REPORT_TITLE,
        author=f"{TOOL_NAME} (openbiota {__version__})",
        subject=f"Sample {chrome.sample}",
    )
    first_frame = Frame(
        MARGIN, FOOTER_H, CONTENT_WIDTH, FIRST_FRAME_HEIGHT, id="first",
        leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0,
    )
    body_frame = Frame(
        MARGIN, FOOTER_H, CONTENT_WIDTH, BODY_FRAME_HEIGHT, id="body",
        leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0,
    )
    doc.addPageTemplates(
        [
            PageTemplate(id="first", frames=[first_frame], onPage=chrome.first),
            PageTemplate(id="interior", frames=[body_frame], onPage=chrome.interior),
        ]
    )
    return doc


def _coherence_note(coh: Any, st: dict[str, ParagraphStyle]) -> Table:
    """A tinted callout: how this reading sits against its twin measurement.

    Agreement is drawn quietly in slate; a split is drawn in amber, because it
    is the one thing on the card a reader most needs to not skim past. The
    text is the same sentence the twin reading prints, so the two halves of
    the report never say different things about the same disagreement.
    """
    agree = coh.agrees
    colour, bg = (SLATE, SLATE_BG) if agree else (AMBER, AMBER_BG)
    lead = "Two measurements agree" if agree else "Two measurements of the same function disagree"
    cell = [
        Paragraph(f"<font color='{_hex(colour)}'><b>{lead}</b></font>", st["small"]),
        Paragraph(coh.explanation(side="capacity"), st["small"]),
    ]
    if coh.out_of_scope_fraction:
        cell.append(Paragraph(
            f"{coh.out_of_scope_fraction:.0%} of this gene family's matched fragments landed on "
            "reference proteins from organisms that carry the enzyme fold without the pathway; those "
            "are counted as background, not as capacity, and are excluded from the figure above.",
            st["fine"],
        ))
    t = Table([[cell]], colWidths=[CONTENT_WIDTH], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg),
        ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
        ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.5 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6 * mm),
    ]))
    return t


def _kv_table(rows: Sequence[tuple[str, str]], st: dict[str, ParagraphStyle], width: float) -> Table:
    data = [
        [Paragraph(k, st["label"]), Paragraph(v, st["cell"])] for k, v in rows
    ]
    table = Table(data, colWidths=[width * 0.42, width * 0.58], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("TOPPADDING", (0, 0), (-1, -1), 1.6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    return table


def _data_table(
    header: Sequence[str],
    rows: Sequence[Sequence[str]],
    st: dict[str, ParagraphStyle],
    widths: Sequence[float],
) -> Table:
    data = [[Paragraph(f"<b>{h}</b>", st["cell"]) for h in header]]
    data.extend([[Paragraph(str(c), st["cell"]) for c in row] for row in rows])
    table = Table(data, colWidths=list(widths), hAlign="LEFT", repeatRows=1)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, ACCENT),
        ("TOPPADDING", (0, 0), (-1, -1), 2.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]
    for index in range(1, len(data)):
        if index % 2 == 0:
            style.append(("BACKGROUND", (0, index), (-1, index), PANEL_BG))
    table.setStyle(TableStyle(style))
    return table


def fine_print(story: list[Any], st: dict[str, ParagraphStyle], items: Sequence[str]) -> None:
    """A fine-print block: thin rule, small caps label, then the small text.

    Disclaimers, confidence grades and caveats live here — present, honest,
    and out of the way of the result.
    """
    items = [i for i in items if i]
    if not items:
        return
    # Set in two columns once there are enough items to balance. Full width
    # this block is tall enough that it often will not fit in what is left
    # of a page, so the whole of it moves — and a detail page that was
    # otherwise complete acquires a second page holding nothing but its
    # caveats, nine tenths of it blank. Two columns roughly halves the
    # height, which is usually the difference between fitting and not.
    body: list[Any]
    if len(items) >= 3:
        # Balance by rendered length rather than by count: one long note and
        # five short ones is not two equal columns.
        weights = [len(t) for t in items]
        target, run, cut = sum(weights) / 2.0, 0, len(items)
        for i, w in enumerate(weights):
            if run + w / 2.0 >= target:
                cut = i
                break
            run += w
        cut = min(max(cut, 1), len(items) - 1)
        gap = 5 * mm
        col = (CONTENT_WIDTH - gap) / 2
        left = [Paragraph(t, st["fine"]) for t in items[:cut]]
        right = [Paragraph(t, st["fine"]) for t in items[cut:]]
        table = Table([[left, right]], colWidths=[col + gap, col], hAlign="LEFT")
        table.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 0),
            ("RIGHTPADDING", (0, 0), (0, 0), gap),
            ("RIGHTPADDING", (1, 0), (1, 0), 0),
            ("TOPPADDING", (0, 0), (-1, -1), 0),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
        ]))
        body = [table]
    else:
        body = [Paragraph(text, st["fine"]) for text in items]

    # Kept whole. Split across a page break this block strands its rule and
    # its "FINE PRINT" label at the foot of one page with the text overleaf,
    # which reads as a rendering fault rather than a caveat.
    story.append(KeepTogether([
        Spacer(1, 2 * mm),
        Rule(CONTENT_WIDTH, RULE, 0.5),
        Spacer(1, 1 * mm),
        Paragraph("FINE PRINT", st["label"]),
        *body,
    ]))


# --------------------------------------------------------------------------- #
# report input
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class MetaboliteRow:
    """Everything the report needs about one metabolite."""

    panel: str
    metabolite: str
    value: float | None
    percentile: float | None
    cohort_median: float | None
    cohort_p25: float | None
    cohort_p75: float | None
    detected: bool
    fragments: int
    status: Status
    what_it_is: str
    made_from: str
    higher_means: str
    evidence_strength: str
    summary: str
    evidence_detail: str
    citation: str
    genes: list[tuple[str, int, float | None]]
    confidence: str = "confirmed"
    implication: str = ""
    extra_note: str = ""
    #: Report grouping (metabolite, neuroactive, vitamin, gas, detox, virulence).
    category: str = "metabolite"
    #: The organisms behind the reading, from the fragment table. None when
    #: the panel had no target fragments or the table is not on disk.
    drivers: Any = None
    #: Why the reading is coloured the way it is, when the panel's own
    #: direction is "it depends" and a usual case or the producers decided it.
    direction_note: str = ""

    def to_json(self) -> dict[str, Any]:
        """The row as data, for `results.json`.

        Every displayed value, so a consumer - this renderer, or the web
        interface that will read the same file - needs to derive nothing.
        The status is carried as its meaning (label, level, whether it was
        placed at all and why not) rather than as a colour: the colour is a
        presentation choice that belongs to whatever is drawing, and a
        different client will make it differently.
        """
        return {
            "panel": self.panel,
            "metabolite": self.metabolite,
            # Rounded to the precision the rest of the file already uses:
            # the panel block stores the value to four places and the
            # reference comparison the percentile to one. Writing more here
            # would leave the file holding two numbers for one reading, and
            # a client that took the wrong one would disagree with the PDF.
            "value": None if self.value is None else round(self.value, 4),
            "percentile": None if self.percentile is None else round(self.percentile, 1),
            "cohort_median": self.cohort_median,
            "cohort_p25": self.cohort_p25,
            "cohort_p75": self.cohort_p75,
            "detected": self.detected,
            "fragments": self.fragments,
            "status": {
                "label": self.status.label,
                "level": self.status.level,
                "assessed": self.status.assessed,
                "note": self.status.note,
                "unassessed_reason": self.status.unassessed_reason,
            },
            "what_it_is": self.what_it_is,
            "made_from": self.made_from,
            "higher_means": self.higher_means,
            "evidence_strength": self.evidence_strength,
            "summary": self.summary,
            "evidence_detail": self.evidence_detail,
            "citation": self.citation,
            "genes": [list(g) for g in self.genes],
            "confidence": self.confidence,
            "implication": self.implication,
            "extra_note": self.extra_note,
            "category": self.category,
            "drivers": None if self.drivers is None else (
                self.drivers.to_json() if hasattr(self.drivers, "to_json")
                else self.drivers
            ),
            "direction_note": self.direction_note,
        }

    @classmethod
    def from_json(cls, payload: Mapping[str, Any]) -> MetaboliteRow:
        """Rebuild a row from `results.json`, deriving nothing.

        The one thing reconstructed is the status *colour*, looked up from
        the label and the direction already recorded. That is presentation,
        not data: no number and no verdict is recomputed here, so a row
        rendered from the file cannot disagree with the file.
        """
        from openbiota.drivers import Drivers  # noqa: PLC0415 - avoids a cycle

        recorded = dict(payload.get("status") or {})
        label = str(recorded.get("label") or "no result")
        level = int(recorded.get("level") or 0)
        note = str(recorded.get("note") or "")
        direction = str(payload.get("higher_means") or "unclear")
        if label in {"no result", "not detected", "not scored"} or not recorded.get(
            "assessed", True
        ):
            colour, background = SLATE, SLATE_BG
        else:
            colour, background = _judgement(level, direction)
        return cls(
            panel=str(payload["panel"]),
            metabolite=str(payload.get("metabolite") or payload["panel"]),
            value=payload.get("value"),
            percentile=payload.get("percentile"),
            cohort_median=payload.get("cohort_median"),
            cohort_p25=payload.get("cohort_p25"),
            cohort_p75=payload.get("cohort_p75"),
            detected=bool(payload.get("detected")),
            fragments=int(payload.get("fragments") or 0),
            status=Status(label, colour, background, note, level),
            what_it_is=str(payload.get("what_it_is") or ""),
            made_from=str(payload.get("made_from") or ""),
            higher_means=direction,
            evidence_strength=str(payload.get("evidence_strength") or "limited"),
            summary=str(payload.get("summary") or ""),
            evidence_detail=str(payload.get("evidence_detail") or ""),
            citation=str(payload.get("citation") or ""),
            genes=[tuple(g) for g in payload.get("genes") or ()],  # type: ignore[misc]
            confidence=str(payload.get("confidence") or "confirmed"),
            implication=str(payload.get("implication") or ""),
            extra_note=str(payload.get("extra_note") or ""),
            category=str(payload.get("category") or "metabolite"),
            # Rebuilt as the record it was, not left as the dict it was
            # stored as: the renderer asks it for attributes.
            drivers=(
                Drivers.from_json(payload["drivers"])
                if isinstance(payload.get("drivers"), Mapping)
                else None
            ),
            direction_note=str(payload.get("direction_note") or ""),
        )


def rows_from_json(results: Mapping[str, Any]) -> list[MetaboliteRow]:
    """Every report row, read from `results.json` rather than recomputed.

    The contract the whole file exists to keep: the analysis writes what it
    found, and the report draws that. A renderer that recomputes a number
    can disagree with the file it was given, and then two artefacts claim
    different things about the same sample.
    """
    rows = results.get("report_rows")
    if not rows:
        raise ValueError(
            "results.json has no 'report_rows'; it was written before the report "
            "rows were part of the file, so the report cannot be drawn from it. "
            "Re-run the sample."
        )
    return [MetaboliteRow.from_json(row) for row in rows]


def _fmt(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "—"
    if value == 0:
        return "0"
    if abs(value) < 0.1:
        return f"{value:.3f}"
    return f"{value:.{digits}f}"


#: Longest glance-row explanation before it is trimmed at a word boundary.


def _glance_sentence(row: MetaboliteRow) -> str:
    """A one-line version of what this reading means.

    The status chip says *what* the reading is; this says what it *means* for
    this particular metabolite. Taken as the leading sentence of the full
    explanation so the glance page and the detail card cannot drift apart.
    """
    if not row.implication:
        if not row.detected:
            return "No reads matched this pathway in your sample."
        if row.percentile is None:
            return "Could not be placed against the reference group."
        return ""

    text = " ".join(row.implication.split())
    cut = len(text)
    for index, char in enumerate(text):
        if char in ".?!" and index + 1 < len(text) and text[index + 1] == " ":
            candidate = text[: index + 1]
            if len(candidate) >= 60:
                cut = index + 1
                break
    # The full first sentence, never clipped: the glance row wraps instead
    # (spec v04 §10: no ellipsis, no truncation anywhere in the report).
    return text[:cut]


def _ordinal(percentile: float | None) -> str:
    """``63rd``; the extremes read ``>99th`` / ``<1st`` because a percentile is
    a position among the reference samples and "100th" would claim a place
    beyond all of them rather than above the ones there are."""
    if percentile is None:
        return "—"
    if percentile >= 99.5:
        return ">99th"
    if percentile == 0:
        # Nothing was detected. "<1st" hedges a fact that is not hedged: the
        # value is zero, which is the bottom of the scale and not a position
        # near the bottom of it.
        return "0"
    if 0 < percentile < 0.5:
        return "<1st"
    n = int(round(percentile))
    suffix = "th"
    if n % 100 not in (11, 12, 13):
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def _arrow(higher_means: str) -> str:
    return {
        "adverse": f' <font color="{_hex(CORAL)}" size="6.5">\u25b2</font>',
        "favourable": f' <font color="{_hex(GREEN)}" size="6.5">\u25bc</font>',
    }.get(higher_means, "")


# --------------------------------------------------------------------------- #
# PART A pages
# --------------------------------------------------------------------------- #


def _at_a_glance(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    rows: Sequence[MetaboliteRow],
    cohort_note: str,
    *,
    section: int = 4,
    detail_section: int | None = None,
    functional_rows: Sequence[dict[str, Any]] | None = None,
) -> None:
    story.extend(section_heading(
        section, "Your Microbial Functions at a Glance", st["h1"]))
    story.append(Paragraph(
        "<b>What your microbes can make, break down and transform.</b>", st["body"]))
    story.append(
        Paragraph(
            f"Where your community's gene capacity for each function falls among {cohort_note}, "
            "grouped by what the function does and most notable first within each group. "
            "The bar is the same seven-band scale used throughout this report; the centre "
            "band is the middle half of the reference group. Every reading is <i>capacity</i> "
            "— the genes are present — not a measured amount of the compound.",
            st["body"],
        )
    )
    story.append(Spacer(1, 2 * mm))

    bar_w = CONTENT_WIDTH * 0.36
    col = [CONTENT_WIDTH * 0.27, bar_w, CONTENT_WIDTH * 0.10, CONTENT_WIDTH * 0.27]

    header = [
        Paragraph("FUNCTION", st["label"]),
        Paragraph("WHERE YOU SIT", st["label"]),
        Paragraph("", st["label"]),
        Paragraph("READING", st["label"]),
    ]
    data: list[list[Any]] = [header]
    style = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
        ("TOPPADDING", (0, 0), (-1, -1), 3.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.0),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]

    # Spec 0.8.3 §12.1: existing and new readings sit together under one set
    # of reader-facing headings, so the grouping comes from there rather than
    # from the older per-panel category vocabulary.
    groups = _pdfextension().grouped_panels(rows)
    # A reading whose genes are exactly what its panel aggregates is that
    # panel said twice. It used to be printed underneath it in a second
    # table with a second number, and the two disagreed wherever the panel
    # summed more than one gene. There is one row per quantity now, and it
    # is the panel's.
    functional_by_group: dict[str, list[dict[str, Any]]] = {}
    for entry in functional_rows or ():
        if entry.get("duplicates_panel"):
            continue
        functional_by_group.setdefault(str(entry["group"]), []).append(entry)

    # One table per category, so a category heading can never be stranded at
    # the foot of a page: the conditional break in front of each asks for room
    # for the heading and its first entry. The column header is drawn once —
    # the categories are self-labelled and the columns do not change.
    base_style = list(style)
    for index, (heading, group) in enumerate(groups):
        data = [header] if index == 0 else []
        style = list(base_style) if index == 0 else [
            rule for rule in base_style if rule[0] != "LINEBELOW"
        ]
        heading_row = len(data)
        data.append([Paragraph(heading.upper(), st["label"]), "", "", ""])
        style.extend(
            [
                ("SPAN", (0, heading_row), (-1, heading_row)),
                ("BACKGROUND", (0, heading_row), (-1, heading_row), PANEL_BG),
                ("TOPPADDING", (0, heading_row), (-1, heading_row), 3.5),
                ("BOTTOMPADDING", (0, heading_row), (-1, heading_row), 3.5),
            ]
        )
        _glance_rows(data, style, group, st, col, bar_w)
        # The further readings for this chemistry go into the *same* table,
        # in the same four columns, with the same bar width and the same
        # chip. They used to be a second table underneath with its own
        # header row, its own widths and its own type sizes, which read as
        # a child layout rather than as more of the same list.
        _glance_reading_rows(
            data, style, functional_by_group.pop(heading, []), st, col, bar_w)
        table = LinkedTable(data, colWidths=col, hAlign="LEFT", repeatRows=1 if index == 0 else 0)
        table.setStyle(TableStyle(style))
        if index:
            story.append(CondPageBreak(34 * mm))
        story.append(table)

    # A group whose further readings have no panel of their own still gets
    # its heading rather than losing the readings.
    for heading in _pdfextension().glance_group_order():
        remaining = functional_by_group.pop(heading, [])
        if not remaining:
            continue
        data = []
        style = [rule for rule in base_style if rule[0] != "LINEBELOW"]
        data.append([Paragraph(heading.upper(), st["label"]), "", "", ""])
        style.extend([
            ("SPAN", (0, 0), (-1, 0)),
            ("BACKGROUND", (0, 0), (-1, 0), PANEL_BG),
            ("TOPPADDING", (0, 0), (-1, 0), 3.5),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 3.5),
        ])
        _glance_reading_rows(data, style, remaining, st, col, bar_w)
        table = LinkedTable(data, colWidths=col, hAlign="LEFT")
        table.setStyle(TableStyle(style))
        story.append(CondPageBreak(34 * mm))
        story.append(table)

    story.append(Spacer(1, 3 * mm))
    story.append(
        Paragraph(
            f'<font color="{_hex(CORAL)}"><b>\u25b2</b></font> research associates higher '
            "levels with worse outcomes &nbsp;&nbsp;·&nbsp;&nbsp; "
            f'<font color="{_hex(GREEN)}"><b>\u25bc</b></font> research associates higher '
            "levels with better outcomes &nbsp;&nbsp;·&nbsp;&nbsp; no arrow: the evidence is "
            "mixed, so the reading is shown in neutral grey. "
            f"Section {detail_section or section} has the findings, the references and the research "
            "on changing each pathway; every row above is a link to its card.",
            st["small"],
        )
    )


def _glance_reading_rows(
    data: list[list[Any]],
    style: list[Any],
    readings: Sequence[Mapping[str, Any]],
    st: dict[str, ParagraphStyle],
    col: Sequence[float],
    bar_w: float,
) -> None:
    """Append the further readings as rows of the group's own table.

    Deliberately the same four cells `_glance_rows` builds - label with its
    direction arrow, the scale bar at the same width, the ordinal in the
    same style, the same chip at the same width - because a reader should
    not be able to tell which rows arrived in which release.
    """
    for row in readings:
        percentile = row.get("percentile")
        direction = str(row.get("higher_means") or "unclear")
        status = status_for(
            percentile=percentile, higher_means=direction,
            detected=percentile is not None,
        )
        marker = RowLink(
            href=_pdfextension().reading_dest(
                str(row["group"]), str(row["label"]), detail=True),
            anchor=_pdfextension().reading_dest(str(row["group"]), str(row["label"])),
            rows=1,
        )
        data.append([
            linked_cell(marker, Paragraph(
                f"<b>{row['label']}</b>{_arrow(direction)}", st["cell"])),
            PercentileBar(
                width=bar_w - 6 * mm, percentile=percentile,
                higher_means=direction, marker_colour=status.colour,
            ),
            Paragraph(
                f"<font color='{_hex(status.colour)}'>{_ordinal(percentile)}</font>",
                st["pct"],
            ),
            StatusChip(status, width=col[3] - 6 * mm),
        ])
        row_index = len(data) - 1
        style.append(("LINEBELOW", (0, row_index), (-1, row_index), 0.3, RULE))


def _glance_rows(
    data: list[list[Any]],
    style: list[Any],
    group: Sequence[MetaboliteRow],
    st: dict[str, ParagraphStyle],
    col: Sequence[float],
    bar_w: float,
) -> None:
    """Append one value row and one meaning row per panel to a glance table."""
    for row in group:
        # The value row and its meaning row are one link to the detail card,
        # and the destination page-1 lands on.
        marker = RowLink(
            href=function_dest(row.panel, detail=True),
            anchor=function_dest(row.panel),
            rows=2,
        )
        data.append(
            [
                linked_cell(
                    marker,
                    Paragraph(f"<b>{row.metabolite}</b>{_arrow(row.higher_means)}", st["cell"]),
                ),
                PercentileBar(
                    width=bar_w - 6 * mm,
                    percentile=row.percentile,
                    higher_means=row.higher_means,
                    marker_colour=row.status.colour,
                ),
                Paragraph(
                    f"<font color='{_hex(row.status.colour)}'>{_ordinal(row.percentile)}</font>",
                    st["pct"],
                ),
                StatusChip(row.status, width=col[3] - 6 * mm),
            ]
        )
        meaning = _glance_sentence(row)
        note_row = len(data)
        data.append([Paragraph(meaning, st["glance_note"]), "", "", ""])
        style.extend(
            [
                ("SPAN", (0, note_row), (-1, note_row)),
                ("TOPPADDING", (0, note_row), (-1, note_row), 0),
                ("BOTTOMPADDING", (0, note_row), (-1, note_row), 4.0),
                ("VALIGN", (0, note_row), (-1, note_row), "TOP"),
                ("LINEBELOW", (0, note_row), (-1, note_row), 0.3, RULE),
            ]
        )


def _detail_pages(  # noqa: PLR0915 - one long card, built in order
    story: list[Any],
    st: dict[str, ParagraphStyle],
    rows: Sequence[MetaboliteRow],
    cohort_note: str,
    *,
    section: int = 10,
    plan: Any | None = None,
) -> None:
    from openbiota.pdfatlas import evidence_cards

    story.extend(section_heading(
        section, "Your Microbial Functions in Detail", st["h1"]))
    story.append(Paragraph(
        "<b>What your microbes can make, break down and transform.</b>", st["body"]))
    story.append(
        Paragraph(
            "One card per pathway, grouped as on the previous page and most notable first "
            "within each group: what the compound is, where you sit, what that means, and "
            "the research behind it.",
            st["body"],
        )
    )
    story.append(Spacer(1, 1 * mm))

    # The same §12.1 grouping as the overview, because the copy above says
    # "grouped as on the previous page" and that has to be true.
    sequence: list[MetaboliteRow | str] = []
    for heading, group in _pdfextension().grouped_panels(rows):
        sequence.append(heading)
        sequence.extend(group)

    for row in sequence:
        if isinstance(row, str):
            story.append(Spacer(1, 1.5 * mm))
            story.append(Paragraph(row, st["h2"]))
            story.append(Spacer(1, 1 * mm))
            continue
        # The card is where page 1 and the at-a-glance row both lead, and it
        # links back to the row so the reader is never stranded here.
        block: list[Any] = [
            Anchor(function_dest(row.panel, detail=True), above=6 * mm),
            back_link_line(function_dest(row.panel), "at a glance"),
        ]
        header = Table(
            [
                [
                    Paragraph(
                        f"<b>{row.metabolite}</b>{_arrow(row.higher_means)}",
                        st["h3"],
                    ),
                    StatusChip(row.status, width=34 * mm, height=5.6 * mm, font_size=6.8),
                ]
            ],
            colWidths=[CONTENT_WIDTH - 36 * mm, 36 * mm],
            hAlign="LEFT",
        )
        header.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 0),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 1),
                ]
            )
        )
        block.append(header)
        block.append(Paragraph(row.what_it_is, st["body"]))

        numbers = _kv_table(
            [
                ("YOUR RESULT", f"<b>{_fmt(row.value)}</b> copies per 100 genomes"),
                (
                    "PERCENTILE",
                    f"<font color='{_hex(row.status.colour)}'><b>{_ordinal(row.percentile)}</b>"
                    f"</font> of {cohort_note}",
                ),
                (
                    "REFERENCE MIDDLE HALF",
                    f"{_fmt(row.cohort_p25)} – {_fmt(row.cohort_p75)} "
                    f"(median {_fmt(row.cohort_median)})",
                ),
                ("MADE FROM", row.made_from),
                ("DNA FRAGMENTS MATCHED", f"{row.fragments:,}"),
            ],
            st,
            CONTENT_WIDTH * 0.52,
        )
        bar = PercentileBar(
            width=CONTENT_WIDTH * 0.42,
            percentile=row.percentile,
            higher_means=row.higher_means,
            marker_colour=row.status.colour,
            height=9 * mm,
            show_scale=True,
            track=4.0 * mm,
        )
        body = Table(
            [[numbers, bar]],
            colWidths=[CONTENT_WIDTH * 0.55, CONTENT_WIDTH * 0.45],
            hAlign="LEFT",
        )
        body.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (0, 0), "TOP"),
                    ("VALIGN", (1, 0), (1, 0), "MIDDLE"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 0),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ]
            )
        )
        block.append(body)

        # The cause, first. A reading is a count of fragments that matched a
        # gene, and every reference gene came from a named organism; which
        # organisms those were is the answer to the question a reader asks
        # before any other. It is set in the status colour, as a table, not
        # in fine print.
        _drivers_block(block, st, row)
        _organism_evidence_block(block, st, row.drivers)

        if row.implication:
            block.append(
                Paragraph(
                    f"<b>What your reading means</b> &nbsp;"
                    f"<font color='{_hex(INK_FAINT)}' size='7'>{row.status.label} · "
                    f"{row.status.note}</font>",
                    st["h3"],
                )
            )
            block.append(Paragraph(row.implication, st["body_ink"]))
        if row.direction_note:
            block.append(Paragraph(
                f"<font color='{_hex(INK_SOFT)}' size='7'><b>Why this colour.</b> "
                f"{row.direction_note}</font>", st["small"]))

        # Where this function is also measured by counting the organisms that
        # do it, say so right here — agreement as corroboration, disagreement
        # with its mechanism — so the two readings are never met apart.
        coh = plan.coherence.by_panel(row.panel) if plan is not None else None
        if coh is not None:
            block.append(_coherence_note(coh, st))

        # Part three of the shared metric card. See `openbiota.metriccard`:
        # every metric carries the same four headings in the same order, and
        # a metric with no established lever says so here rather than
        # dropping the section.
        _metriccard().improvement_block(
            block, st,
            higher_means=row.higher_means,
            organisms=[
                {"organism": d.organism, "sample_percent": d.sample_percent}
                for d in (row.drivers.top if row.drivers else ()) if d.in_sample
            ],
        )

        strength = STRENGTH_WORDS.get(row.evidence_strength, row.evidence_strength)
        block.append(
            Paragraph(
                f"<b>The research behind it</b> &nbsp;<font color='{_hex(INK_FAINT)}' size='7'>"
                f"{DIRECTION_WORDS.get(row.higher_means, '')} · {strength}</font>",
                st["h3"],
            )
        )
        block.append(Paragraph(row.summary, st["body"]))
        block.append(Paragraph(row.evidence_detail, st["body"]))
        notes = [row.extra_note, f"Reference: {row.citation}"]
        for note in notes:
            if note:
                block.append(Paragraph(note, st["fine"]))
        story.append(KeepTogether(block))
        if plan is not None:
            evidence_cards(story, st, plan.cards_for("panel", row.panel))
        story.append(Spacer(1, 1.5 * mm))
        story.append(Rule(CONTENT_WIDTH))
        story.append(Spacer(1, 3.5 * mm))


# --------------------------------------------------------------------------- #
# validation + technical
# --------------------------------------------------------------------------- #


def _performance_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    validation: dict[str, Any] | None,
    cohort: dict[str, Any] | None,
    depth: dict[str, Any] | None = None,
    *,
    section: int = 16,
) -> None:
    story.extend(section_heading(section, "How this test performs", st["h1"]))
    story.append(
        Paragraph(
            "Every figure below was measured, not estimated. The method was tested against "
            "synthetic samples built from reference genomes whose gene content is known "
            "exactly, so the correct answer was known in advance and could be compared "
            "against what the pipeline reported.",
            st["body"],
        )
    )

    if validation:
        metrics = validation.get("overall_metrics", {})
        counts = validation.get("overall_confusion", {})

        def pct(value: Any) -> str:
            return "—" if value is None else f"{value * 100:.1f}%"

        tiles = [
            ("SENSITIVITY", pct(metrics.get("sensitivity")), "genes present that were found"),
            ("SPECIFICITY", pct(metrics.get("specificity")), "genes absent, correctly absent"),
            ("PRECISION", pct(metrics.get("precision")), "of detections that were real"),
            ("FALSE POSITIVES", pct(metrics.get("false_positive_rate")), "spurious detection rate"),
        ]
        width = (CONTENT_WIDTH - 3 * 3 * mm) / 4
        story.append(
            Table(
                [[Tile(value=v, label=k, note=n, width=width, accent=GREEN) for k, v, n in tiles]],
                colWidths=[width + 3 * mm] * 3 + [width],
                hAlign="LEFT",
                style=TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0),
                                  ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                                  ("TOPPADDING", (0, 0), (-1, -1), 0),
                                  ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]),
            )
        )
        story.append(Spacer(1, 2 * mm))
        story.append(
            Paragraph(
                f"Based on {len(validation.get('genomes', []))} reference genomes across "
                f"{len(validation.get('experiments', []))} synthetic samples "
                f"({counts.get('TP', 0)} correct detections, {counts.get('TN', 0)} correct "
                f"non-detections, {counts.get('FP', 0)} false positives, "
                f"{counts.get('FN', 0)} missed).",
                st["small"],
            )
        )

        calibration = validation.get("calibration") or {}
        if calibration.get("r_squared") is not None:
            story.append(Paragraph("Are the numbers themselves accurate?", st["h2"]))
            story.append(
                Paragraph(
                    "Beyond simply detecting a gene, the reported quantity was compared "
                    "against the known quantity in each synthetic sample. A slope of 1.0 "
                    "would mean the reported figure matches truth exactly, and an "
                    f"R\u00b2 of 1.0 would mean a perfect relationship. Measured: slope "
                    f"<b>{calibration['slope']:.2f}</b>, R\u00b2 "
                    f"<b>{calibration['r_squared']:.2f}</b> across "
                    f"{calibration['n']} comparisons.",
                    st["body"],
                )
            )

        if depth:
            story.append(Paragraph("Was your sample sequenced deeply enough?", st["h2"]))
            story.append(
                Paragraph(
                    "Every fragment in your sample was read, but that raises a fair question: "
                    "was the sequencing itself deep enough that a rarer organism could not "
                    "have been missed? To check, the same sample was re-analysed at "
                    "progressively lower depths. If the results stop changing before full "
                    "depth is reached, then depth was not the limiting factor and more "
                    "sequencing would not alter these numbers.",
                    st["body"],
                )
            )
            depths = sorted({p["read_pairs"] for p in depth.get("points", [])})
            by_key = {(p["panel"], p["read_pairs"]): p for p in depth.get("points", [])}
            panels_seen = sorted({p["panel"] for p in depth.get("points", [])})
            if depths and panels_seen:
                header = ["Pathway"] + [f"{d // 1000:,}k" for d in depths] + ["Settled?"]
                converged = {c["panel"]: c for c in depth.get("convergence", [])}
                rows_out: list[list[str]] = []
                for name in panels_seen:
                    cells = [name]
                    for d in depths:
                        point = by_key.get((name, d))
                        value = point.get("copies_per_100") if point else None
                        cells.append("—" if value is None else f"{value:.1f}")
                    verdict = converged.get(name, {})
                    cells.append(
                        "yes"
                        if verdict.get("converged")
                        else ("too few reads" if "too few" in str(verdict.get("verdict")) else "no")
                    )
                    rows_out.append(cells)
                width = CONTENT_WIDTH / (len(depths) + 2)
                story.append(
                    _data_table(header, rows_out, st, [width * 1.4]
                                + [width * 0.85] * len(depths) + [width * 1.3])
                )
            story.append(Paragraph(depth.get("summary", ""), st["body_ink"]))

        limits = validation.get("detection_limit") or []
        if limits:
            story.append(Paragraph("How rare an organism can this still find?", st["h2"]))
            rows = [
                [
                    f"{row['carrier_cell_fraction']:.2%}",
                    f"{row['true_copies_per_100']:.1f}",
                    "found" if row["detected"] else "missed",
                    f"{row['reported_fragments']:,}",
                ]
                for row in limits
            ]
            story.append(
                _data_table(
                    ["Carrier abundance", "True copies/100", "Result", "Fragments"],
                    rows,
                    st,
                    [
                        CONTENT_WIDTH * 0.28, CONTENT_WIDTH * 0.24,
                        CONTENT_WIDTH * 0.22, CONTENT_WIDTH * 0.26,
                    ],
                )
            )

        story.append(Paragraph("Performance by pathway", st["h2"]))
        rows = []
        for name, m in sorted(validation.get("per_panel_metrics", {}).items()):
            c = m.get("confusion", {})
            rows.append(
                [
                    name,
                    "—" if m.get("sensitivity") is None else f"{m['sensitivity'] * 100:.0f}%",
                    "—" if m.get("specificity") is None else f"{m['specificity'] * 100:.0f}%",
                    "—" if m.get("precision") is None else f"{m['precision'] * 100:.0f}%",
                    f"{c.get('TP', 0)}/{c.get('TN', 0)}/{c.get('FP', 0)}/{c.get('FN', 0)}",
                ]
            )
        story.append(
            _data_table(
                ["Pathway", "Sensitivity", "Specificity", "Precision", "TP/TN/FP/FN"],
                rows,
                st,
                [
                    CONTENT_WIDTH * 0.22, CONTENT_WIDTH * 0.19, CONTENT_WIDTH * 0.19,
                    CONTENT_WIDTH * 0.17, CONTENT_WIDTH * 0.23,
                ],
            )
        )
    else:
        story.append(
            Paragraph(
                "Validation figures were not available when this report was produced. "
                "Run <b>openbiota validate</b> to generate them.",
                st["body"],
            )
        )

    if cohort:
        story.append(Paragraph("The reference group", st["h2"]))
        story.append(Paragraph(cohort.get("description", ""), st["body"]))

    story.append(Paragraph("What is still uncertain", st["h2"]))
    story.append(
        Paragraph(
            "The measured figures above come from synthetic samples, where the answer is "
            "known. They establish that the method finds what is there and does not invent "
            "what is not. They do not establish that any of these pathways predicts a health "
            "outcome — that is a question about the biology, not about the test, and for most "
            f"of these metabolites it remains open. The research summaries in section {SECTIONS['functions_detail']} give "
            "the current state of the evidence, including how strong it is in each case.",
            st["body"],
        )
    )


def _organism_evidence_block(
    block: list[Any], st: dict[str, ParagraphStyle], drv: Any
) -> None:
    """What is known about the organisms actually carrying this reading.

    A reading is high or low because of organisms, and an action that does
    not name them is advice about somebody else. This lists the species
    behind the signal that were found in *this* sample and what the curated
    literature says about each, in the direction it says it.

    Where the record is association only it says so. Inventing a way to
    move an organism that nothing has shown how to move would be worse
    than saying nothing, and this report is read by people making
    decisions.
    """
    from openbiota import organisms as org  # noqa: PLC0415

    present = [d for d in (drv.top if drv else ()) if d.in_sample]
    if not present:
        return
    registry = org._species_registry()  # noqa: SLF001 - the curated record
    rows: list[list[Any]] = []
    for d in present[:6]:
        entry = registry.get(str(d.organism).replace(" ", "_"))
        if not entry:
            continue
        assertions = [a for a in (entry.get("assertions") or []) if a.get("statement")]
        state = str(entry.get("evidence_state") or "")
        lines = [f"<font size='6.6'>{entry.get('summary', '').strip()}</font>"]
        for a in assertions[:2]:
            arrow = {"higher": "when higher", "lower": "when lower",
                     "presence": "when present"}.get(str(a.get("direction")), "")
            lines.append(
                f"<font size='6.2' color='{_hex(INK_SOFT)}'>\u25aa {arrow}: "
                f"{a['statement']} <i>[{a.get('citation', '')}]</i></font>"
            )
        if state == "human_association_only":
            lines.append(
                f"<font size='6.2' color='{_hex(AMBER)}'>Association only \u2014 no study "
                "on record shows that changing this organism changes an outcome.</font>"
            )
        rows.append([
            Paragraph(
                f"<i><b>{d.organism}</b></i><br/>"
                f"<font size='6' color='{_hex(INK_FAINT)}'>{d.sample_percent:.2f}% of your "
                f"community \u00b7 {d.share:.0%} of this signal</font>"
                if d.sample_percent is not None else f"<i><b>{d.organism}</b></i>",
                st["cell"]),
            Paragraph("<br/>".join(lines), st["cell"]),
        ])
    if not rows:
        return
    block.append(Spacer(1, 1.5 * mm))
    block.append(Paragraph(
        "<b>The organisms behind this, and what is known about them</b>", st["h3"]))
    header = [[Paragraph("FOUND IN YOU", st["label"]),
               Paragraph("WHAT THE LITERATURE RECORDS", st["label"])]]
    table = Table(header + rows, colWidths=[CONTENT_WIDTH * 0.26, CONTENT_WIDTH * 0.74],
                  hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
        ("TOPPADDING", (0, 0), (-1, -1), 3.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.0),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
    ]))
    block.append(table)


def _drivers_verdict(block: list[Any], st: dict[str, ParagraphStyle], drv: Any) -> None:
    """Answer the reader's question before showing them the table.

    "What is driving this" followed by a list of organisms that are not in
    the sample is not an answer. The colorectal virulence reading named
    *Escherichia coli* as the whole of its signal, marked "not found", on
    four fragments - a reader is entitled to ask whether they have E. coli,
    and the table alone does not say.

    So the bottom line goes first: how much of the signal belongs to
    organisms found in this person, named, and how much belongs to
    reference sequences nothing here matches.
    """
    total = sum(d.fragments for d in drv.top) or 1
    present = [d for d in drv.top if d.in_sample]
    genus_only = [d for d in drv.top if not d.in_sample and d.genus_in_sample]
    absent = [d for d in drv.top if not d.in_sample and not d.genus_in_sample]
    share_present = sum(d.fragments for d in present) / total

    if not present and not genus_only:
        names = ", ".join(f"<i>{d.organism}</i>" for d in absent[:3])
        block.append(Paragraph(
            f"<font color='{_hex(CORAL)}'><b>No organism behind this reading was found in "
            f"your sample.</b></font> The {drv.total_fragments:,} matched "
            f"fragment{'' if drv.total_fragments == 1 else 's'} came closest to {names}, "
            "which the taxonomic profile did not detect in you. A read resembling an "
            "organism's copy of a gene is not that organism being present \u2014 some other "
            "member of your community carries a similar sequence, or the reference set "
            "lacks whichever one does. Treat this as a trace signal with no identified "
            "source rather than as evidence about those species.", st["body"]))
        return

    parts: list[str] = []
    if present:
        named = "; ".join(
            f"<i>{d.organism}</i>"
            + (f" ({d.sample_percent:.2f}% of your community)"
               if d.sample_percent is not None else "")
            for d in present[:4]
        )
        parts.append(
            f"<b>{share_present:.0%} of the signal comes from organisms found in you:</b> "
            f"{named}."
        )
    if genus_only:
        parts.append(
            f"A further {sum(d.fragments for d in genus_only) / total:.0%} matched "
            + ", ".join(f"<i>{d.organism}</i>" for d in genus_only[:3])
            + ", whose genus is present but whose species was not resolved."
        )
    if absent:
        parts.append(
            f"The remaining {sum(d.fragments for d in absent) / total:.0%} matched reference "
            "sequences from organisms not detected in you, which is a resemblance rather "
            "than a presence."
        )
    block.append(Paragraph(" ".join(parts), st["body"]))


def _drivers_block(block: list[Any], st: dict[str, ParagraphStyle], row: MetaboliteRow) -> None:
    """Which organisms the matched DNA came from, as a table in the status colour."""
    drv = row.drivers
    if drv is None or not drv.top:
        if row.detected:
            block.append(Paragraph(
                "<b>What is driving this</b>", st["h3"]))
            block.append(Paragraph(
                "The fragment table for this pathway was not available when the report was "
                "built, so the organisms behind the reading could not be named here.",
                st["small"]))
        return

    from openbiota import organisms as org

    colour = row.status.colour if row.status.colour not in (SLATE,) else INK_SOFT
    block.append(Paragraph(
        f"<b>What is driving this</b> &nbsp;<font color='{_hex(INK_FAINT)}' size='7'>"
        f"{drv.total_fragments:,} matched fragments, by the organism each came closest to</font>",
        st["h3"]))
    _drivers_verdict(block, st, drv)

    class_colour = {
        org.BENEFICIAL: GREEN, org.OPPORTUNIST: CORAL, org.CONDITIONAL: AMBER, org.UNKNOWN: SLATE,
    }
    rows: list[list[Any]] = [[
        Paragraph("ORGANISM", st["label"]), Paragraph("SHARE OF SIGNAL", st["label"]),
        Paragraph("IN YOUR SAMPLE", st["label"]), Paragraph("CLASS", st["label"]),
    ]]
    for d in drv.top:
        if d.in_sample:
            where = (f"<font color='{_hex(colour)}'><b>yes</b></font> \u00b7 "
                     f"{d.sample_percent:.3f}%" if d.sample_percent is not None else
                     f"<font color='{_hex(colour)}'><b>yes</b></font>")
        elif d.genus_in_sample:
            where = f"<font color='{_hex(INK_SOFT)}'>genus present; closest reference</font>"
        else:
            where = f"<font color='{_hex(INK_FAINT)}'>not found; closest reference</font>"
        cls = d.sample_class
        cls_cell = (
            f"<font color='{_hex(class_colour.get(cls, SLATE))}' size='8'>\u25cf</font> "
            f"<font size='6.8'>{org.CLASS_LABEL.get(cls, '')}</font>" if cls else
            f"<font color='{_hex(INK_FAINT)}'>\u2014</font>"
        )
        rows.append([
            Paragraph(f"<i><b>{d.organism}</b></i>", st["cell"]),
            Paragraph(
                f"<font color='{_hex(colour)}'><b>{d.share:.0%}</b></font> "
                f"<font color='{_hex(INK_FAINT)}' size='6.4'>({d.fragments} fragments)</font>",
                st["cell"]),
            Paragraph(where, st["cell"]),
            Paragraph(cls_cell, st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.36, 0.20, 0.26, 0.18)]
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
        ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2.2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
    ]))
    block.append(t)
    lead = drv.top[0]
    present = [d for d in drv.top if d.in_sample]
    if present:
        names = ", ".join(f"<i>{d.organism}</i>" for d in present[:3])
        block.append(Paragraph(
            f"<font size='7.2'>The organisms carrying this gene in your sample are {names}"
            + (f"; <i>{lead.organism}</i> alone accounts for {lead.share:.0%} of the signal."
               if lead.in_sample else ".")
            + " Changing this reading means changing them.</font>",
            st["small"]))
    else:
        block.append(Paragraph(
            "<font size='7.2'>None of the closest references is itself in your organism list: the "
            "fragments most likely come from close relatives that carry a similar gene. The "
            "genus is named where it was found.</font>",
            st["small"]))
    block.append(Spacer(1, 1.5 * mm))


def _detection_catalogues(results: dict[str, Any]) -> str:
    """Every catalogue the organism inventory drew on, with its release.

    Stated exactly, once, here. The reader-facing sections say "every
    catalogue available" and never name a version; this is where the
    versions live for anyone who needs to reproduce or compare.
    """
    parts: list[str] = []
    ext = results.get("extended_catalogue") or {}
    if ext:
        n = ext.get("n_sgbs_in_database") or 36_822
        parts.append(f"{ext.get('database', 'CHOCOPhlAnSGB')} ({n:,} species-level bins)")
    gen = results.get("genome_profile") or {}
    if gen.get("status") == "resolved":
        parts.append(
            f"GTDB r232 via {gen.get('tool', 'sylph')} "
            f"({gen.get('catalogue_size', 199_923):,} species clusters, "
            f"ANI \u2265 {gen.get('minimum_ani', 95):.0f}%)"
        )
    elif gen.get("status") == "not_assessed":
        parts.append("GTDB whole-genome lane: not assessed (not installed)")
    return "; ".join(parts) if parts else "\u2014"


def _table(rows: list[list[Any]], widths: Sequence[float]) -> Table:
    """A zebra table of prepared Paragraph cells, in the catalogue's style."""
    from openbiota.pdflibrary import _plain_table

    return _plain_table(rows, widths, zebra=True)


def _detection_references(story: list[Any], st: dict[str, ParagraphStyle], results: Mapping[str, Any]) -> None:
    """The reference releases and tool versions every detection lane used.

    Spec 0.8.4 §7: the completed report states measured capability and the
    reference releases used. One row per lane: tool and version, reference
    release, taxonomy release, what it observed here. Then the locked
    reference sources with their licences.
    """
    det = results.get("detection") or {}
    lanes = det.get("lanes") or {}
    status = det.get("lane_status") or {}
    if not lanes and not status:
        return
    story.append(Paragraph("Detection lanes and reference releases", st["h2"]))
    story.append(Paragraph(
        "Every organism in this report was found by one or more of these methods, each run on the same "
        "host-filtered reads against a fixed reference release. A lane recorded as not run says nothing "
        "about which organisms it would have found.",
        st["small"]))
    rows: list[list[Any]] = [[
        Paragraph("LANE", st["label"]), Paragraph("TOOL", st["label"]), Paragraph("REFERENCE", st["label"]),
        Paragraph("TAXONOMY", st["label"]), Paragraph("RESULT", st["label"]),
    ]]
    words = {
        "metaphlan_jan26": "Marker lane", "sylph_globdb": "Broad discovery", "motus4": "Universal markers",
        "kraken_uhgg": "Gut catalogue classifier", "kraken_rescue": "Gut rescue panel", "singlem_globdb": "Lineage rescue",
    }
    for lane_id in ("metaphlan_jan26", "sylph_globdb", "motus4", "kraken_uhgg", "kraken_rescue", "singlem_globdb"):
        blob = lanes.get(lane_id) or {}
        stt = status.get(lane_id) or {}
        if blob:
            result = f"{blob.get('n_species', 0)} species-level observations"
        else:
            result = f"<font color='{_hex(CORAL)}'>{stt.get('state', 'not run')}</font>"
        rows.append([
            Paragraph(f"<b>{words.get(lane_id, lane_id)}</b>", st["cell"]),
            Paragraph(f"{blob.get('tool', '')} {blob.get('tool_version', '')}".strip() or "—", st["fine"]),
            Paragraph(str(blob.get("reference_release_id") or "—"), st["fine"]),
            Paragraph(str(blob.get("taxonomy_release") or "—"), st["fine"]),
            Paragraph(result, st["fine"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.16, 0.22, 0.22, 0.24, 0.16)]))
    releases = det.get("reference_releases") or []
    if releases:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph("Reference sources locked for this run", st["h3"]))
        rrows: list[list[Any]] = [[
            Paragraph("SOURCE", st["label"]), Paragraph("RELEASE", st["label"]), Paragraph("COUNT", st["label"]),
            Paragraph("LICENCE", st["label"]), Paragraph("VERIFIED", st["label"]),
        ]]
        for r in releases:
            count = f"{r['source_count']:,} {r.get('count_unit') or ''}".strip() if r.get("source_count") else (r.get("count_unit") or "—")
            rrows.append([
                Paragraph(str(r.get("source_id") or ""), st["cell"]), Paragraph(str(r.get("release_id") or ""), st["fine"]),
                Paragraph(count, st["fine"]), Paragraph(str(r.get("license") or "")[:60], st["fine"]),
                Paragraph("yes" if r.get("verified") else "size only", st["fine"]),
            ])
        story.append(_table(rrows, [CONTENT_WIDTH * w for w in (0.16, 0.30, 0.18, 0.26, 0.10)]))
        story.append(Paragraph(
            "GlobDB propagates CC BY-SA 4.0; mOTUs data are CC BY 4.0; GTDB is CC BY-SA 4.0; BacDive is CC BY 4.0 "
            "and asks to be contacted for commercial use. Full URLs, checksums and retrieval times are in "
            "refs/expanded/locks and in results.json under detection.reference_releases.",
            st["fine"]))
    cap = det.get("measured_capability") or {}
    if cap.get("n_communities"):
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph("Measured capability", st["h3"]))
        ci = cap.get("pooled_precision_ci95") or [None, None]
        neg = cap.get("negatives") or {}
        text = (
            f"On {cap['n_communities']} simulated communit{'y' if cap['n_communities'] == 1 else 'ies'} of "
            f"{cap.get('species_per_community', '?')} species at "
            f"{int(cap.get('pairs') or 0):,} read pairs with known composition (benchmark <i>{cap.get('name', '')}</i>, "
            f"{str(cap.get('run_at', ''))[:10]}), detection found "
            f"{100 * (cap.get('mean_recall_all_lanes') or 0):.1f}% of the species present, and "
            f"{100 * (cap.get('pooled_precision_all_lanes') or 0):.1f}% of the species it reported as supported were present"
            + (f" (95% CI {100 * ci[0]:.1f}\u2013{100 * ci[1]:.1f}%)" if ci[0] is not None else "")
            + f"; the species reported in error ({cap.get('mean_false_supported_per_community', 0):.1f} per community) were "
            "relatives of species that were present."
        )
        if neg.get("n"):
            text += (f" On {neg['n']} organism-free host and food backgrounds, {neg.get('false_supported', 0)} "
                     "microbial species were reported as supported.")
        text += " These are measurements on simulated data, not a statement about this specimen."
        story.append(Paragraph(text, st["fine"]))
    story.append(Spacer(1, 3 * mm))


def _technical_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    results: dict[str, Any],
    rows: Sequence[MetaboliteRow],
    similarity: Any | None = None,
    *,
    section: int = 17,
) -> None:
    story.extend(section_heading(section, "Technical data", st["h1"]))
    _detection_references(story, st, results)

    story.append(Paragraph("Per-gene results", st["h2"]))
    table_rows: list[list[str]] = []
    for row in rows:
        for gene, fragments, copies in row.genes:
            table_rows.append(
                [row.panel, gene, f"{fragments:,}", _fmt(copies), _fmt(row.value)]
            )
    story.append(
        _data_table(
            ["Pathway", "Gene", "Fragments", "Copies/100", "Pathway total"],
            table_rows,
            st,
            [CONTENT_WIDTH * 0.20] * 5,
        )
    )

    story.append(Paragraph("Sample and run details", st["h2"]))
    run = results.get("run", {})
    sample_qc = results.get("input") or {}
    normalisation = results.get("normalisation", {})
    mates = sample_qc.get("mates") or [{}]
    detail = [
        ("Sample", str(results.get("sample", "—"))),
        ("Read pairs", f"{sample_qc.get('read_pairs') or 0:,}"),
        (
            "Total bases sequenced",
            f"{(sample_qc.get('estimated_total_bases') or 0) / 1e9:.2f} Gbp",
        ),
        ("Read length", f"{mates[0].get('read_length_min', '—')}–{mates[0].get('read_length_max', '—')} bp"),
        ("GC content", f"{mates[0].get('gc_percent', '—')}%"),
        ("Mean quality", f"Phred {mates[0].get('mean_phred', '—')}"),
        ("Calibration gene fragments", f"{normalisation.get('fragments', 0):,} rpoB"),
        ("Reference database", str(run.get("reference proteins", "—"))),
        ("Search tool", str(run.get("diamond", "—"))),
        ("Analysis time", str(run.get("wall clock", "—"))),
        (
            "Pipeline",
            f'<a href="{BRAND_URL}" color="{_hex(ACCENT_DARK)}">{BRAND}</a>'
            f"{TOOL_NAME.removeprefix(BRAND)} (openbiota) v{results.get('version', '—')}",
        ),
        ("Reference set fingerprint", str(run.get("reference fingerprint", "—"))),
    ]
    if similarity is not None:
        tax = similarity.taxonomy
        manifest = similarity.cohort_manifest or {}
        detail += [
            (
                "Organism identification",
                f"{tax.family} {tax.profiler_version}" if tax is not None else "not run",
            ),
            ("Taxonomic database", tax.index if tax is not None else "—"),
            (
                "Detection catalogues",
                _detection_catalogues(results),
            ),

            (
                "Host reads removed",
                (
                    "not assessable — removed by the provider before delivery"
                    if tax.host.upstream_removed
                    else f"{tax.host.host_pairs:,} pairs ({tax.host.host_fraction:.2%})"
                )
                if tax is not None and tax.host is not None else "not filtered",
            ),
            (
                "Species reference cohort",
                f"{manifest.get('n_samples', 0):,} samples, snapshot "
                f"{manifest.get('snapshot', '—')}, id {manifest.get('manifest_id', '—')}",
            ),
            ("Reference matching", str(similarity.references.match_info.get("reason", "—"))),
            (
                "Profiles scored",
                ", ".join(f"{r.profile.name} v{r.profile.version}" for r in similarity.results)
                or "none",
            ),
            (
                "Dysbiosis anchor",
                f"GMWI2 {similarity.anchor.score:+.2f}"
                if similarity.anchor is not None and similarity.anchor.valid else "—",
            ),
        ]
    story.append(_kv_table(detail, st, CONTENT_WIDTH * 0.75))

    story.append(Paragraph("How the calculation works", st["h2"]))
    story.append(
        Paragraph(
            "Each pathway gene's fragment count is divided by the average length of that "
            "gene, then divided again by the same quantity for <i>rpoB</i> — a gene every "
            "bacterium carries exactly once. That ratio gives copies per bacterial genome, "
            "which is then multiplied by 100.",
            st["body"],
        )
    )
    story.append(
        Paragraph(
            "Two consequences are worth knowing. Sequencing depth cancels out, so a deeper "
            "sample does not produce bigger numbers. And human DNA in the sample cancels out "
            "too, because human cells carry neither the pathway genes nor <i>rpoB</i> — which "
            "is why Part A needs no human-DNA removal step.",
            st["body"],
        )
    )
    if similarity is not None:
        story.append(
            Paragraph(
                "Part B works differently. Species abundances are compositional — they sum to "
                "100%, so one going up forces others down — and are therefore compared in "
                "centred log-ratio space, against the median and spread of the reference "
                "cohort, using a robust scale that a single outlier cannot dominate. Each "
                "feature's deviation is clipped, signed by the published direction, and "
                "weighted by four documented factors; module scores are weighted means; and "
                "percentiles come from scoring every reference sample the same way, never from "
                "rescaling. Human DNA does <i>not</i> cancel here, so it is removed first.",
                st["body"],
            )
        )

    story.append(Paragraph("Full method and source code", st["h2"]))
    story.append(
        Paragraph(
            "The complete method, the validation protocol and every gene query used are "
            "documented alongside the source code. The machine-readable form of this report, "
            "containing every number shown here plus a good deal more, is in "
            "<b>results.json</b>.",
            st["small"],
        )
    )


# --------------------------------------------------------------------------- #
# assembly
# --------------------------------------------------------------------------- #




def _metriccard() -> Any:
    """The shared metric-card contract, imported on use.

    `metriccard` reads this module's palette, so importing it at the top
    would be circular. Every metric's detail page is built through it: see
    its docstring before adding a metric anywhere in this report.
    """
    from openbiota import metriccard  # noqa: PLC0415

    return metriccard


def _pdfextension() -> Any:
    """The v0.8.3 rendering module, imported on use.

    `pdfextension` reads this module's palette and flowables, so importing it
    at the top would be circular. It is only needed while building a report.
    """
    from openbiota import pdfextension

    return pdfextension


def _pdfinputs() -> Any:
    """The A14 section renderer, imported lazily like the other extensions."""
    from openbiota import pdfinputs  # noqa: PLC0415

    return pdfinputs


def _extension_view(results: Mapping[str, Any] | None, name: str) -> Mapping[str, Any] | None:
    """One view from the v0.8.3 extension, or None when it did not run.

    The extension is additive: a report renders identically without it, and
    every block that reads through here returns early on None.
    """
    return (((results or {}).get("extension") or {}).get("views") or {}).get(name)


def build_pdf(
    *,
    path: Path,
    results: dict[str, Any],
    rows: Sequence[MetaboliteRow],
    profile: dict[str, Any] | None,
    validation: dict[str, Any] | None,
    cohort: dict[str, Any] | None,
    cohort_note: str,
    depth: dict[str, Any] | None = None,
    similarity: Any | None = None,
    profile_validation: dict[str, Any] | None = None,
    plan: Any | None = None,
    age: Any | None = None,
    age_meta: dict[str, Any] | None = None,
    # Retained in the signature because callers pass it; the QC tables now
    # render from `results["sequencing_quality"]` in their own section, and
    # the compact status reaches the input register through the results object.
    gates: Any | None = None,  # noqa: ARG001
    mode: str = "research",
    subject_age: float | None = None,
) -> Path:
    """Render the Gut Health Test Metagenomic Report.

    Part A (at a glance) puts every reading in front of the reader as a
    graphic; Part B (in detail) breaks each one down and, under each, shows
    what the human research says about changing it.

    The rows are taken from ``results["report_rows"]`` whenever the results
    object carries them, in preference to the ``rows`` argument. The two are
    built from the same panels and should agree; reading the file is what
    guarantees they do, and it is the same path the web interface will take.
    """
    if results and results.get("report_rows"):
        from_file = rows_from_json(results)
        if rows and len(from_file) != len(rows):
            raise ValueError(
                f"results.json holds {len(from_file)} report rows but the renderer was "
                f"passed {len(rows)}. One of them is stale, and drawing either would "
                "produce a report that disagrees with the file beside it."
            )
        rows = from_file
    from openbiota import (
        pdfatlas,
        pdfattention,
        pdfbiofilm,
        pdflibrary,
        pdfmycobiome,
        pdfpathogens,
        pdfprofiles,
        pdfskin,
        pdfstrains,
        pdfsummary,
    )
    st = _styles()
    # The merged organism inventory: every catalogue that ran, under current
    # names. Sections that list organisms read this rather than one lane.
    from openbiota import inventory as _inventory
    inv = _inventory.from_json(results.get("organism_inventory"))
    if inv is None and similarity is not None and similarity.taxonomy is not None:
        # A results file written before the inventory existed, or a run where
        # the merge failed: rebuild from what is on hand rather than lose the
        # listing entirely.
        try:
            inv = _inventory.from_results(
                results,
                scoring_species=similarity.taxonomy.species,
                ranked_rows=getattr(
                    getattr(similarity, "community", None), "species", ()) or (),
            )
        except Exception:  # noqa: BLE001 - the single-lane listing still renders
            inv = None

    sample = str(results.get("sample", "sample"))
    generated = dt.datetime.now().strftime("%-d %B %Y")

    sample_qc = results.get("input") or {}
    read_pairs = sample_qc.get("read_pairs") or 0
    read_text = f"{read_pairs / 1e6:.1f} million" if read_pairs else "—"
    meta_line = " · ".join(
        [
            f"Sample {sample}",
            f"Report {generated}",
            f"{read_text} read pairs" if read_pairs else "",
            f"openbiota v{__version__}",
        ]
    ).replace(" ·  · ", " · ")
    chrome = _Chrome(sample=sample, generated=generated, meta_line=meta_line)

    path.parent.mkdir(parents=True, exist_ok=True)
    doc = _document(path, chrome)

    meta = {
        "sample": sample,
        "date": generated,
        "n_panels": len(rows),
        "read_pairs": read_text,
        "read_pairs_raw": read_pairs,
        "cohort_note": cohort_note,
        "n_species": (
            len(inv) if inv is not None and len(inv)
            else sum(1 for v in similarity.taxonomy.species.values() if v > 0)
            if similarity is not None and similarity.taxonomy is not None
            else None
        ),
        "n_profiles": len(similarity.results) if similarity is not None else 0,
    }

    ordered = sorted(rows, key=lambda r: notability(r.status, r.percentile))

    story: list[Any] = []
    findings = plan.findings if plan is not None else None
    pathogens = results.get("pathogens")
    # What each named organism is, whether carrying it is ordinary, and what
    # has published activity against it. Optional: absent files mean the
    # cards report what they measured and nothing more, rather than failing.
    try:
        from openbiota.pathogens.agents import load_agent_library

        agent_library = load_agent_library(Path("pathogens/catalog")) or None
    except Exception:  # noqa: BLE001 - advice is never worth losing a report over
        agent_library = None
    pdfsummary.summary_page(
        story, st, rows=ordered, similarity=similarity, rpob_profile=profile, meta=meta, findings=findings,
        age=age, age_meta=age_meta, pathogens=pathogens, pathogen_section=SECTIONS["pathogens"],
        results=results, inv=inv,
    )

    taxo_note = (
        f"{similarity.references.taxonomic_n:,} matched reference adults"
        if similarity is not None and similarity.references.taxonomic_n
        else "the reference group"
    )
    S = SECTIONS

    # ------------------------------------------------------------------ #
    # PART A — AT A GLANCE: every reading as a graphic, nothing deep yet.
    # Order: the community and what it is missing → how old it reads →
    # what it can make → what it resembles → what can be done.
    # ------------------------------------------------------------------ #
    story.append(NextPageTemplate("interior"))
    story.append(PageBreak())
    # The contents come before the findings: a reader who has just met the
    # first page wants to know where everything is, and the page that tells
    # them was behind the one that assumes they already know.
    pdfatlas.reading_guide_page(story, st, section=S["guide"])

    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; AT A GLANCE", st["kicker"]))
    pdfsummary.insights_page(
        story, st, rows=ordered, similarity=similarity, meta=meta, findings=findings, age=age, plan=plan,
        age_meta=age_meta, inv=inv,
    )

    # 1 · your microbiome make-up: summary infographic, groups, organisms
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; MICROBIOME DIVERSITY: WHAT'S THERE AND WHAT'S MISSING", st["kicker"]))
    # The neutral range strips that used to follow the species list are gone:
    # this section now runs from the community picture straight to the one
    # index that scores it, and the reader meets a number rather than five
    # unscored strips on the way.
    pdfprofiles.community_page(
        story, st, section=S["community"], similarity=similarity, rpob_profile=profile,
        inv=inv,
        # The composition's phyla come from whichever marker lane is primary:
        # Jan26 when it ran, else the Jun23 catalogue.
        primary_phyla=(
            (((results.get("detection") or {}).get("lanes") or {}).get("metaphlan_jan26") or {}).get("summary") or {}
        ).get("phyla") or (results.get("extended_catalogue") or {}).get("phyla"),
    )

    # A09: the companion diversity, dominance, ratio and aerotolerance views
    # belong beside the existing community readings, not in a chapter of
    # their own (spec 0.8.3 §12.1).
    _pdfextension().composition_block(story, st, _extension_view(results, "ecology"))
    _pdfextension().ecology_block(story, st, _extension_view(results, "ecology"))

    # The groups table follows the community page directly when the index
    # block has left the room; a fresh page otherwise.
    # One list feeds the organisms glance and the cards in Part B: every
    # organism the classification flags plus every one the census expected
    # and did not find. Built once so the two sections cannot disagree. The
    # set of species that get a card is also what the group rows link to.
    attention_items = pdfattention.build(inv, findings)
    carded = frozenset(a.species for a in attention_items)

    story.append(CondPageBreak(120 * mm))
    story.append(Spacer(1, 4 * mm))
    pdfatlas.groups_glance(story, st, section=S["groups"], similarity=similarity, cohort_note=taxo_note,
                           inv=inv, carded=carded)

    story.append(PageBreak())
    pdfattention.glance(story, st, section=S["organisms"], detail_section=S["organisms_detail"],
                        items=attention_items, findings=findings,
                        n_all_detected=len(inv) if inv is not None else 0, cohort_note=taxo_note)

    # The full organism list, directly after the flagged ones. Having just
    # read which organisms need attention, the reader's next question is what
    # else is in there, and the answer should not be at the back of the report.
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; EVERY ORGANISM DETECTED", st["kicker"]))
    pdflibrary.catalogue_page(
        story, st, section=S["catalogue"],
        community=None if similarity is None else similarity.community,
        cohort_note=taxo_note, inv=inv,
    )

    # A10: the searchable organism view, after the full catalogue it extends.
    _pdfextension().explorer_block(story, st, _extension_view(results, "explorer"))

    # 1c · pathogens (spec v5.0 §13). Placed straight after "who is there",
    # because the reader's first question about a stool metagenome is whether
    # anything in it causes disease.
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; PATHOGENS: DISEASE-CAUSING ORGANISMS", st["kicker"]))
    pdfpathogens.pathogens_glance(
        story, st, pathogens=pathogens, section=S["pathogens"],
        detail_section=S["pathogens_detail"],
    )

    # 2 · biological age
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; ESTIMATED BIOLOGICAL AGE OF BIOTA", st["kicker"]))
    pdfatlas.age_page(
        story, st, section=S["age"], age=age, age_meta=age_meta, subject_age=subject_age, mode=mode,
    )

    # 3 · metabolic functions
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; WHAT YOUR BACTERIA CAN MAKE AND BREAK DOWN", st["kicker"]))
    # Spec 0.8.3 §12.1: existing panel totals and the new A01-A08 and §5.8
    # readings sit together under one set of reader-facing headings, in this
    # overview as well as in the detail pages. Without the second argument
    # the panel totals were the whole of the overview and every substrate,
    # route and transformation reading existed only in the detail pages.
    _at_a_glance(
        story, st, ordered, cohort_note,
        section=S["functions"], detail_section=S["functions_detail"],
        functional_rows=_pdfextension().all_functional_rows(
            ((results or {}).get("extension") or {}).get("views"),
            {row.panel: row.higher_means for row in ordered},
            # What each panel aggregates, so a reading that restates one can
            # be recognised and left out rather than printed beside it as a
            # second opinion on the same measurement.
            {
                str(p.get("name")): p.get("aggregate_from") or ()
                for p in ((results or {}).get("panels") or [])
                if isinstance(p, dict)
            }),
    )

    # 4 · disease patterns
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; RESEMBLANCE TO PUBLISHED DISEASE PATTERNS", st["kicker"]))
    pdfatlas.profiles_glance(story, st, section=S["patterns"], similarity=similarity)

    # 4a · biofilm-related potential (spec v8.0 §13.1). Front-loaded among the
    # visual metric pages rather than buried after the prose, because the two
    # axes only make sense next to each other.
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; BIOFILM-RELATED POTENTIAL", st["kicker"]))
    pdfbiofilm.biofilm_overview(
        story, st, section=S["biofilm"],
        biofilm=(results.get("biofilm") if results else None),
        detail_section=S["biofilm_detail"],
    )

    # The mycobiome page follows the biofilm page: one complete overview
    # among the early visual summaries (spec §3), detail in Part B.
    story.append(PageBreak())
    pdfmycobiome.mycobiome_overview(
        story, st, section=S["mycobiome"],
        myco=(results.get("mycobiome") if results else None),
        detail_section=S["mycobiome_detail"],
    )

    # 4b · inflammatory skin research panels (spec v4.2 §7.2). Front-loaded,
    # paginated in full: no top-N truncation.
    story.append(PageBreak())
    story.append(Paragraph("PART A &mdash; SKIN-CONDITION RESEARCH PATTERNS", st["kicker"]))
    from openbiota import skinaxis  # noqa: PLC0415

    _skin_axis = skinaxis.score(
        (results or {}).get("profile_similarity") if results else None
    )
    pdfskin.skin_glance(
        story, st, section=S["skin"],
        skin=None if similarity is None else getattr(similarity, "skin", None),
        skin_axis=_skin_axis,
    )

    # The actions section opens with the levers: every agent the organism registry
    # knows, added up across the organisms flagged in this sample and ranked
    # by how many it moves the right way. The registry's guideline and
    # exact-product cards follow under their own heading.
    story.append(PageBreak())
    # A12: the consolidated view across every evidence card, before the
    # per-organism levers it summarises (spec 0.8.3 §8.2).
    _pdfextension().planner_block(story, st, _extension_view(results, "planner"))
    pdfattention.levers(story, st, section=S["actions"], items=attention_items,
                        detail_section=S["organisms_detail"])
    pdfatlas.evidence_overview(story, st, section=S["actions"], plan=plan, as_subsection=True)

    # A16: the model-assisted comparison of a probiotic, a prebiotic and both
    # together belongs with the actions it informs (spec 0.8.3 §11.1).
    _pdfextension().simulation_block(story, st, _extension_view(results, "simulation"))
    _pdfextension().owed_block(
        story, st, ((results or {}).get("extension") or {}).get("unavailable") or [],
    )

    # ------------------------------------------------------------------ #
    # PART B — IN DETAIL, in the same order as Part A
    # ------------------------------------------------------------------ #
    pdfatlas.section_divider(
        story, st, "Part B", "In detail",
        "Every reading from Part A, one at a time and in the same order: what it is, exactly where you sit, what "
        "that means, the research behind the reading — and, under each one, <b>what the human research says about "
        "changing it</b>: the exact products, diets, drugs or clinical routes that have been studied, in whom, with "
        "the evidence for and against side by side and the safety rules that apply. Where nothing has been shown "
        "to work, the card says so. Where fecal microbiota transplantation is the only route with evidence, "
        "the card says that and gives its regulatory status.",
        contents=(
            f"{S['groups_detail']} · Your microbial groups in detail, and every species detected",
            f"{S['organisms_detail']} · Organisms that need attention, one by one",
            f"{S['pathogens_detail']} · Pathogens, one at a time, and what was not assessed",
            f"{S['functions_detail']} · Your metabolic functions in detail ({len(rows)} pathways)",
            f"{S['patterns_detail']} · Disease-pattern resemblance, every scored pattern",
            f"{S['skin_detail']} · Gut-skin axis, panel by panel",
            f"{S['context']} · Your information, and what this report assumed",
            f"{S['limits']} · What none of this can tell you · {S['sequencing']} · Sequencing quality and reference coverage",
            f"{S['accuracy']} · Measured accuracy · {S['technical']} · Technical record",
        ),
    )

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; WHO IS THERE", st["kicker"]))
    pdflibrary.community_groups_pages(
        story, st, section=S["groups_detail"],
        community=None if similarity is None else similarity.community,
        cohort_note=taxo_note, inv=inv, carded=carded,
        plan=plan,
    )

    story.append(PageBreak())
    pdfattention.detail(story, st, section=S["organisms_detail"], items=attention_items, plan=plan,
                        bacdive=results.get("bacdive"), strain_analysis=results.get("strain_analysis"))

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; PATHOGENS, ONE AT A TIME", st["kicker"]))
    # Strain resolution sits between the organism inventory and the pathogen
    # detail: it is the same organisms, one level finer.
    story.append(PageBreak())
    pdfstrains.strain_pages(
        story, st, section=S["strains"],
        strain=(results.get("strain_resolution") if results else None),
    )
    pdfstrains.placement_pages(
        story, st, analysis=(results.get("strain_analysis") if results else None),
        inventory=(results.get("organism_inventory") if results else None),
    )
    pdfstrains.census_pages(
        story, st, census=(results.get("resolution_census") if results else None),
    )

    pdfpathogens.pathogens_detail_pages(
        story, st, pathogens=pathogens, section=S["pathogens_detail"],
        agents=agent_library,
    )

    # A11: the class, mechanism and carrier view over the existing AMR calls.
    _pdfextension().resistance_block(story, st, _extension_view(results, "resistance"))

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; WHAT YOUR BACTERIA CAN MAKE AND BREAK DOWN", st["kicker"]))
    _detail_pages(story, st, ordered, cohort_note, section=S["functions_detail"], plan=plan)
    # A01-A08: the new substrate, route and transformation readings, grouped
    # as spec 0.8.3 §12.1 asks, inside the section they belong to.
    _pdfextension().functional_lane_block(
        story, st, ((results or {}).get("extension") or {}).get("views"),
        # The same aggregate sets the overview used, so both sides exclude
        # the same readings and no link points at a row that is not there.
        {
            str(p.get("name")): p.get("aggregate_from") or ()
            for p in ((results or {}).get("panels") or [])
            if isinstance(p, dict)
        },
    )
    # Spec 0.8.3 §5.8's drug-metabolism cards, and the beta-glucuronidase
    # substrate detail, which states why it is withheld when its class
    # assignment cannot be reproduced from the published structures.
    _pdfextension().drug_metabolism_block(
        story, st, ((results or {}).get("extension") or {}).get("views"),
    )
    # A16: the mechanism behind the synbiotic summary lives with the metabolism
    # detail it explains (spec 0.8.3 §11.1, AT114).
    _pdfextension().simulation_detail_block(story, st, _extension_view(results, "simulation"))

    # §6.4 closes the functions section: how many organisms hold each of the
    # functions above. It sat in the community section, a chapter away from
    # the readings it qualifies.
    _pdfextension().redundancy_block(story, st, _extension_view(results, "ecology"))

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; DISEASE-PATTERN RESEMBLANCE", st["kicker"]))
    pdflibrary.library_pages(
        story, st, section=S["patterns_detail"], similarity=similarity,
        profile_validation=profile_validation, plan=plan,
        copri_complex=(results.get("copri_complex") if results else None),
    )
    # A13's symptom and condition navigation. Spec 0.8.3 §12.2 places this
    # inside the pattern section rather than in a catch-all of its own, so
    # that a topic sits next to the scores it links to.
    _pdfextension().context_block(
        story, st, _extension_view(results, "contexts"),
    )

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; BIOFILM FINDINGS IN FULL", st["kicker"]))
    pdfbiofilm.biofilm_detail(
        story, st, section=S["biofilm_detail"],
        biofilm=(results.get("biofilm") if results else None),
    )
    story.append(PageBreak())
    pdfmycobiome.mycobiome_detail(
        story, st, section=S["mycobiome_detail"],
        myco=(results.get("mycobiome") if results else None),
    )

    story.append(PageBreak())
    story.append(Paragraph("PART B &mdash; SKIN-CONDITION RESEARCH PATTERNS", st["kicker"]))
    pdfskin.skin_detail_pages(
        story, st, section=S["skin_detail"],
        skin=None if similarity is None else getattr(similarity, "skin", None),
        plan=plan,
    )

    # A14: one combined register of everything supplied, assumed and absent,
    # placed immediately before the limits it explains (spec 0.8.3 §10.1).
    story.append(PageBreak())
    _pdfinputs().inputs_section(
        story, st, _extension_view(results, "input_register"), section=S["context"],
    )

    story.append(PageBreak())
    pdfprofiles.limits_page(story, st, section=S["limits"])

    story.append(PageBreak())
    pdflibrary.sequencing_quality_page(
        story, st, section=S["sequencing"],
        quality=results.get("sequencing_quality"),
        extended=results.get("extended_catalogue"),
    )

    story.append(PageBreak())
    _performance_page(story, st, validation, cohort, depth, section=S["accuracy"])
    pdfprofiles.performance_extension(
        story, st, profile_validation=profile_validation, similarity=similarity
    )

    story.append(PageBreak())
    _technical_page(story, st, results, rows, similarity, section=S["technical"])
    # A15 §12.3: the complete index of the new measurements belongs in the
    # technical area, not in another results section.
    _pdfextension().metric_index_block(story, st, (results or {}).get("extension"))

    doc.build(story)
    return path
