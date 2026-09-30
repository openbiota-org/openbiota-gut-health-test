"""Report pages for the pathogen branch — spec v5.0 §13.

Structure
---------
* `pathogen_alert` — the page-1 badge. A red count and an evil bug when
  something is supported; a calm neutral panel when nothing is.
* `pathogens_glance` — the Part A opener: the six organism groups as cards,
  then every supported and candidate finding in one table.
* `pathogens_detail_pages` — one card per finding, then determinants, then
  the coverage ledger that says what was *not* assessed.

Presentation rules that are not cosmetic
----------------------------------------
* Red is reserved for `supported_sequence` and `marker_signal` findings.
  A candidate signal is amber, an ambiguous one grey. Nothing that failed
  to be assessed is ever coloured as a negative, because "we could not
  look" is not "we looked and it was clean".
* The count in the badge is supported findings only. Candidates are shown
  next to it under their own label, never folded into the same number.
* Groups with nothing found still appear, with the number screened, so an
  empty group reads as evidence rather than as an omission.
* A `not_assessed` target never renders as absent. The coverage ledger is
  part of the result, not an appendix.
* Microsporidia appear once, under their own heading (spec P082); they are
  never also counted as fungi.
* Every colour is paired with a word, and the bug glyph is decorative —
  the count is always written out.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import (
    CondPageBreak,
    Flowable,
    Spacer,
    Table,
    TableStyle,
)

from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    back_link_line,
    linked_cell,
    pathogen_dest,
    section_dest,
    section_heading,
)

CONTENT_WIDTH: Final = 170 * mm

INK: Final = colors.HexColor("#15222E")
INK_SOFT: Final = colors.HexColor("#4B5B6B")
INK_FAINT: Final = colors.HexColor("#8C99A6")
RULE: Final = colors.HexColor("#E1E7EC")
PANEL_BG: Final = colors.HexColor("#F5F8FA")

#: Supported findings. The only states that earn red.
ALERT: Final = colors.HexColor("#C0392B")
ALERT_DARK: Final = colors.HexColor("#8E2A1F")
ALERT_BG: Final = colors.HexColor("#FCEFED")
#: Candidate signals: real sequence, not enough of it to name the organism.
WATCH: Final = colors.HexColor("#D98324")
WATCH_DARK: Final = colors.HexColor("#8A5410")
WATCH_BG: Final = colors.HexColor("#FDF4E8")
#: Nothing found, and the search was actually capable of finding it.
CLEAR: Final = colors.HexColor("#2E8B6F")
CLEAR_DARK: Final = colors.HexColor("#1F6B54")
CLEAR_BG: Final = colors.HexColor("#EDF7F3")
#: Ambiguous, or not assessed. Grey means "no claim", never "clean".
QUIET: Final = colors.HexColor("#64748B")
QUIET_BG: Final = colors.HexColor("#EEF1F4")

#: Determinant statuses that mean the gene's sequence is present. Mirrors
#: `openbiota.pathogens.determinants.PRESENT_STATUSES`; kept as a literal here so
#: the report layer never imports the engine just to render a table.
DETERMINANT_PRESENT: Final[frozenset[str]] = frozenset(
    {"supported_intact_sequence", "supported_partial_sequence", "ambiguous_allele"}
)

#: Statuses that count towards the headline number (spec §13.1).
SUPPORTED: Final[frozenset[str]] = frozenset({"supported_sequence", "marker_signal"})
#: Shown prominently, counted separately.
CANDIDATE: Final[frozenset[str]] = frozenset({"candidate_signal"})
AMBIGUOUS: Final[frozenset[str]] = frozenset({"ambiguous_signal"})

#: The six organism groupings the report presents, in reading order.
#: Viral DNA and RNA lanes merge into one heading for the reader and split
#: again in the detail cards, where the molecule actually changes the caveat.
DISPLAY_GROUPS: Final[tuple[tuple[str, tuple[str, ...], str], ...]] = (
    ("Bacterial pathogens", ("bacteria",), "bacteria"),
    ("Protozoa pathogens", ("protozoa",), "protozoa"),
    ("Parasitic worms", ("helminths",), "worm"),
    ("Microsporidia pathogens", ("microsporidia",), "spore"),
    ("Fungal pathogens", ("fungi", "other_eukaryotes"), "fungus"),
    (
        "Viral pathogens",
        ("dna_viruses", "rna_viruses", "other_viruses"),
        "virus",
    ),
)

_STATUS_WORDS: Final[Mapping[str, str]] = {
    "supported_sequence": "enough sequence to name it",
    "marker_signal": "named from marker sequence",
    "candidate_signal": "trace, not enough to confirm",
    "ambiguous_signal": "shared sequence only — cannot tell apart",
    "not_detected": "not found",
    "not_assessed": "could not be checked",
    "technical_artifact": "excluded as laboratory artifact",
}

# --------------------------------------------------------------------------- #
# What the organism actually is
# --------------------------------------------------------------------------- #
# The catalogue already records this per seed in `interpretation_class`, and
# an earlier version of this section ignored it: every named organism was
# printed under "Pathogens found", which put brewer's yeast and an ordinary
# gut anaerobe in the same red list as Shigella. Naming a normal resident a
# pathogen is not a presentation problem, it is a wrong result.

#: Tier keys in reading order, with the heading and the one-line rule.
TIER_ORDER: Final[tuple[str, ...]] = (
    "pathogen", "pathotype_negative", "unresolved_complex", "opportunist",
    "uncertain", "elsewhere", "normal",
)

TIERS: Final[Mapping[str, Mapping[str, Any]]] = {
    "pathogen": {
        "heading": "Disease-causing organisms",
        "blurb": (
            "Recognised causes of gut illness. A finding here is worth acting on, "
            "and worth showing to a doctor."
        ),
        "colour": ALERT, "background": ALERT_BG, "counts_as_pathogen": True,
    },
    "opportunist": {
        "heading": "Opportunists — harmless in most people, a problem in some",
        "blurb": (
            "These live in many healthy guts and cause disease mainly when "
            "something else has changed: recent antibiotics, a weakened immune "
            "system, surgery, or a large overgrowth. Presence alone is not illness."
        ),
        "colour": WATCH, "background": WATCH_BG, "counts_as_pathogen": False,
    },
    "uncertain": {
        "heading": "Common residents whose role is still unsettled",
        "blurb": (
            "Frequently found in healthy people. Research links some of them to "
            "disease, but whether they cause it is not established. Reported for "
            "completeness, not as a problem to fix."
        ),
        "colour": QUIET, "background": QUIET_BG, "counts_as_pathogen": False,
    },
    "elsewhere": {
        "heading": "Organisms that cause disease outside the gut",
        "blurb": (
            "Stool is the wrong specimen for these. Sequence here does not "
            "indicate infection of another organ, and its absence does not rule "
            "one out."
        ),
        "colour": QUIET, "background": QUIET_BG, "counts_as_pathogen": False,
    },
    "pathotype_negative": {
        "heading": "Carried, but without the genes that cause disease",
        "blurb": (
            "For these organisms the disease is caused by a toxin or a pathotype, not by the "
            "species. Their DNA is here; the genes that would make them dangerous were looked "
            "for and not found. Carrying the harmless form is common, and for C. difficile it "
            "is associated with protection from infection rather than with disease."
        ),
        "colour": QUIET, "background": QUIET_BG, "counts_as_pathogen": False,
    },
    "unresolved_complex": {
        "heading": "DNA that cannot be pinned to a species",
        "blurb": (
            "These organisms are the same genomic species as a relative that lives in most "
            "healthy guts, so alignment cannot tell them apart and a species name would be a "
            "guess. The evidence is shown in full, with the marker that would settle it. In "
            "stool, sequence like this usually comes from the harmless relative."
        ),
        "colour": QUIET, "background": QUIET_BG, "counts_as_pathogen": False,
    },
    "normal": {
        "heading": "Normal residents, food and environmental organisms",
        "blurb": (
            "Not pathogens. These are ordinary members of a healthy gut, or they "
            "arrive with food. They are listed because the screen searches for "
            "them by name, and finding them is expected rather than concerning."
        ),
        "colour": CLEAR, "background": CLEAR_BG, "counts_as_pathogen": False,
    },
}

#: `interpretation_class` from the catalogue to reporting tier.
_CLASS_TO_TIER: Final[Mapping[str, str]] = {
    # Recognised enteric pathogens, and organisms that are pathogens when they
    # carry the right toxin or pathotype genes — the determinant screen reports
    # whether they do, and the card says which way it came out.
    "established_enteric": "pathogen",
    "rare_enteric": "pathogen",
    "toxin_or_pathotype_dependent": "pathogen",
    "conditional_opportunist": "opportunist",
    "uncertain_enteric_role": "uncertain",
    "extraintestinal_watch": "elsewhere",
    # Decoys and background are in the catalogue precisely so that reads which
    # belong to them stop being credited to a pathogen. Finding one is the
    # system working, not a result.
    "background_or_decoy": "normal",
}


def tier_of(record: Mapping[str, Any]) -> str:
    """Which reporting tier one finding belongs in.

    Unknown classes fall to `uncertain` rather than `pathogen`: if the
    catalogue has not said an organism causes disease, this report must not
    say it either.
    """
    # The resolution pass (openbiota.pathogens.resolve) has both the organism
    # and the determinant evidence, so where it has ruled on a row its tier
    # wins: a species the data cannot name, and carriage without the toxin
    # gene, are not pathogen findings however the catalogue classes them.
    explicit = record.get("report_tier")
    # Except when the row has no organism-specific fragment at all. "DNA
    # present, species cannot be named — 0 specific fragments, 0 regions" is a
    # sentence that contradicts itself, and it read as the most alarming line
    # on the page for a row carrying no evidence. Guarding here as well as in
    # the resolution pass means a stored verdict from an older run cannot
    # print it either.
    if explicit in {"unresolved_complex", "pathotype_negative"} and not int(
        record.get("unique_supporting_fragments") or 0
    ):
        explicit = None
    if explicit and str(explicit) in TIERS:
        return str(explicit)
    klass = str(record.get("interpretation_class") or "")
    tier = _CLASS_TO_TIER.get(klass)
    if tier is not None:
        return tier
    # A dietary or environmental role is enough on its own to keep something
    # out of the pathogen list even when its class is unfamiliar.
    if str(record.get("stool_role") or "") == "environmental_or_dietary":
        return "normal"
    return "uncertain"

_RESOLUTION_WORDS: Final[Mapping[str, str]] = {
    "species": "species level",
    "species_complex": "species complex",
    "genus": "genus level",
    "group": "group level",
    "subspecies": "subspecies level",
    "serovar": "serovar level",
    "pathotype": "pathotype level",
    "family": "family level",
    "unresolved": "could not be resolved",
}

_CLASS_WORDS: Final[Mapping[str, str]] = {
    "established_enteric": "established gut pathogen",
    "toxin_or_pathotype_dependent": "depends on toxin or pathotype",
    "conditional_opportunist": "opportunist, context dependent",
    "rare_enteric": "rarely reported in gut disease",
    "uncertain_enteric_role": "role in gut disease uncertain",
    "extraintestinal_watch": "not a gut-disease organism",
    "background_or_decoy": "background, food or environmental",
}


# --------------------------------------------------------------------------- #
# The all-clear mark
# --------------------------------------------------------------------------- #


class AllClearMark(Flowable):
    """A round green checkmark for the page-1 band when nothing was found.

    The bug glyph says "pathogen"; a screen that found none should not open
    with a germ, however green. This is the shape everyone reads as "clear":
    a filled disc on a soft halo ring, with a white tick drawn as one stroke
    with round caps and joins. Drawn as vectors, like everything else on the
    page, and decorative - the words beside it carry the meaning.
    """

    def __init__(self, size: float = 9 * mm, colour: colors.Color = CLEAR) -> None:
        super().__init__()
        self.size = size
        self.colour = colour
        self.width = size
        self.height = size

    def wrap(self, _aw: float, _ah: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c = self.canv
        s = self.size
        cx = cy = s / 2.0
        c.saveState()
        # Halo: a wide, pale ring behind the disc so it sits softly on the
        # band rather than as a hard dot.
        halo = colors.Color(self.colour.red, self.colour.green, self.colour.blue, alpha=0.18)
        c.setFillColor(halo)
        c.setStrokeColor(halo)
        c.circle(cx, cy, s * 0.50, stroke=0, fill=1)
        # Disc.
        c.setFillColor(self.colour)
        c.circle(cx, cy, s * 0.39, stroke=0, fill=1)
        # Tick: short down-stroke, long up-stroke, one path, rounded.
        c.setStrokeColor(colors.white)
        c.setLineWidth(max(0.9, s * 0.085))
        c.setLineCap(1)
        c.setLineJoin(1)
        p = c.beginPath()
        p.moveTo(cx - s * 0.175, cy + s * 0.005)
        p.lineTo(cx - s * 0.045, cy - s * 0.125)
        p.lineTo(cx + s * 0.19, cy + s * 0.135)
        c.drawPath(p, stroke=1, fill=0)
        c.restoreState()


# --------------------------------------------------------------------------- #
# The bug
# --------------------------------------------------------------------------- #


class EvilBug(Flowable):
    """The pathogen glyph: an iconic microbe, drawn as vectors so it scales cleanly.

    A pear-shaped cell with vacuole spots, a gloss highlight near the top and
    pili radiating all the way round — the shape everyone already reads as
    "germ". It has to say pathogen at a glance without saying threat: no
    jaws, no scowl, nothing sharp beyond the pili, which are round-capped.

    Drawn rather than shipped as a bitmap because the report is vector
    throughout. The glyph is decorative — every place it appears also states
    the count in words, so the meaning survives greyscale and screen readers.
    `menace` is kept for callers and no longer changes the drawing.
    """

    #: Pili around the cell as (angle in degrees, length as a fraction of
    #: size, curl). Hand-placed so the spacing looks organic rather than
    #: mechanical; the curl bends each one slightly, alternating sides.
    _PILI: tuple[tuple[float, float, float], ...] = (
        (96, 0.12, 0.35), (66, 0.11, -0.3), (36, 0.13, 0.3), (10, 0.11, -0.35),
        (-18, 0.13, 0.3), (-46, 0.12, -0.3), (-74, 0.11, 0.35), (-104, 0.12, -0.3),
        (-134, 0.13, 0.3), (-162, 0.11, -0.35), (170, 0.13, 0.3), (142, 0.12, -0.3),
        (122, 0.11, 0.35),
    )

    #: Lesion-like patches as (x, y, radius, stretch, tilt in degrees), in
    #: units of size, relative to centre. Each is drawn as an uneven blob,
    #: not a circle, in a lighter shade of the body rather than white — a
    #: diseased, mottled surface rather than polka dots.
    _SPOTS: tuple[tuple[float, float, float, float, float], ...] = (
        (0.05, -0.21, 0.058, 1.35, 20.0), (-0.11, -0.13, 0.040, 1.25, -35.0),
        (0.12, -0.02, 0.030, 1.20, 60.0), (-0.03, 0.02, 0.046, 1.30, -10.0),
        (-0.09, -0.27, 0.026, 1.20, 40.0), (0.06, 0.13, 0.030, 1.40, -25.0),
        (-0.07, 0.15, 0.022, 1.20, 15.0), (0.00, 0.24, 0.018, 1.30, 70.0),
    )

    @staticmethod
    def _wobble(a: float, seed: float) -> float:
        """Smooth radius multiplier around a blob; `seed` varies the shape."""
        return 1.0 + 0.11 * math.sin(3.0 * a + seed) + 0.07 * math.cos(5.0 * a - seed)

    def __init__(
        self,
        size: float = 14 * mm,
        colour: colors.Color = ALERT,
        *,
        menace: bool = True,
    ) -> None:
        super().__init__()
        self.size = size
        self.colour = colour
        self.menace = menace
        self.width = size
        self.height = size

    def wrap(self, _aw: float, _ah: float) -> tuple[float, float]:
        return self.width, self.height

    @staticmethod
    def _edge(theta: float) -> tuple[float, float]:
        """Point on the cell outline at angle `theta`, in units of size.

        A pear: full at the bottom, narrower toward the top. The half-width
        shrinks with height, which is all it takes to stop it being an egg.
        """
        k = (math.sin(theta) + 1.0) / 2.0  # 0 at the bottom, 1 at the top
        half_width = 0.29 - 0.11 * k
        return half_width * math.cos(theta), 0.36 * math.sin(theta)

    def draw(self) -> None:  # noqa: PLR0915 - one contiguous vector drawing
        c = self.canv
        s = self.size
        col = self.colour
        # A softer dark and a lighter patch colour: the earlier version drew
        # the outline and pili at 62% and the spots in pure white, which
        # made the glyph heavier and busier than anything else on the page.
        dark = colors.Color(col.red * 0.72, col.green * 0.72, col.blue * 0.72)
        patch = colors.Color(
            col.red + (1 - col.red) * 0.42,
            col.green + (1 - col.green) * 0.42,
            col.blue + (1 - col.blue) * 0.42,
        )
        cx, cy = s * 0.5, s * 0.49

        def at(x: float, y: float) -> tuple[float, float]:
            return cx + s * x, cy + s * y

        c.saveState()
        c.setLineCap(1)
        c.setLineJoin(1)

        # --- pili: drawn first so their roots tuck under the body ---------- #
        c.setStrokeColor(dark)
        c.setLineWidth(s * 0.038)
        for angle, length, curl in self._PILI:
            theta = math.radians(angle)
            ex, ey = self._edge(theta)
            # Outward direction: from a point a little below centre, so the
            # pili on the narrow top fan out rather than converge.
            nx, ny = ex, ey + 0.04
            norm = math.hypot(nx, ny) or 1.0
            nx, ny = nx / norm, ny / norm
            tx, ty = -ny * curl, nx * curl  # tangential offset for the bend
            root = at(ex - nx * 0.03, ey - ny * 0.03)
            tip = at(ex + nx * length + tx * length * 0.6, ey + ny * length + ty * length * 0.6)
            ctrl = at(ex + nx * length * 0.5, ey + ny * length * 0.5)
            path = c.beginPath()
            path.moveTo(*root)
            path.curveTo(*ctrl, *ctrl, *tip)
            c.drawPath(path, stroke=1, fill=0)

        # --- body: the pear outline as a smooth polygon ------------------- #
        c.setFillColor(col)
        c.setStrokeColor(dark)
        c.setLineWidth(s * 0.016)
        body = c.beginPath()
        steps = 72
        for i in range(steps + 1):
            theta = 2.0 * math.pi * i / steps
            x, y = at(*self._edge(theta))
            if i == 0:
                body.moveTo(x, y)
            else:
                body.lineTo(x, y)
        body.close()
        c.drawPath(body, stroke=1, fill=1)

        # --- mottled patches and a soft highlight -------------------------- #
        c.setFillColor(patch)
        points = 24
        for x, y, r, stretch, tilt in self._SPOTS:
            px, py = at(x, y)
            tilt_r = math.radians(tilt)
            blob = c.beginPath()
            for i in range(points):
                a = 2.0 * math.pi * i / points
                rr = s * r * self._wobble(a, tilt_r)
                # An ellipse, stretched along its tilt, with the wobble on top.
                ex, ey = rr * stretch * math.cos(a), rr * math.sin(a)
                qx = px + ex * math.cos(tilt_r) - ey * math.sin(tilt_r)
                qy = py + ex * math.sin(tilt_r) + ey * math.cos(tilt_r)
                if i == 0:
                    blob.moveTo(qx, qy)
                else:
                    blob.lineTo(qx, qy)
            blob.close()
            c.drawPath(blob, stroke=0, fill=1)
        c.setFillColor(colors.Color(1, 1, 1, alpha=0.55))
        hx, hy = at(-0.065, 0.245)
        c.roundRect(hx, hy, s * 0.13, s * 0.045, s * 0.0225, stroke=0, fill=1)

        c.restoreState()


class Chip(Flowable):
    """A small rounded count chip: number in ink, label beside it."""

    def __init__(
        self,
        text: str,
        colour: colors.Color,
        background: colors.Color,
        *,
        width: float | None = None,
        height: float = 4.9 * mm,
        size: float = 6.7,
    ) -> None:
        super().__init__()
        self.text = text
        self.colour = colour
        self.background = background
        self.height = height
        self.size = size
        # Measure the label rather than trusting a number typed by hand. Fixed
        # widths meant a longer label either overflowed the row or sat in a
        # chip three times its size; a row of these has to fit the page.
        self.width = (
            width if width is not None
            else stringWidth(text, "Helvetica-Bold", size) + 5.0 * mm
        )

    def wrap(self, _aw: float, _ah: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        c = self.canv
        c.saveState()
        c.setFillColor(self.background)
        c.setStrokeColor(self.colour)
        c.setLineWidth(0.5)
        c.roundRect(0, 0, self.width, self.height, self.height * 0.5, stroke=1, fill=1)
        c.setFillColor(self.colour)
        c.setFont("Helvetica-Bold", self.size)
        # Centre the cap height, not the baseline: Helvetica caps are ~0.72 em
        # tall, so this puts the text in the middle of the chip.
        baseline = (self.height - self.size * 0.72) / 2
        c.drawCentredString(self.width / 2, baseline, self.text)
        c.restoreState()


# --------------------------------------------------------------------------- #
# Reading the result payload
# --------------------------------------------------------------------------- #


def _records(pathogens: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Every organism result record, defensively."""
    if not pathogens:
        return []
    found = pathogens.get("results") or pathogens.get("targets") or []
    return [r for r in found if isinstance(r, dict)]


def _determinants(pathogens: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not pathogens:
        return []
    found = pathogens.get("determinants") or []
    return [r for r in found if isinstance(r, dict)]


def _status(record: Mapping[str, Any]) -> str:
    return str(record.get("display_status") or record.get("sequence_status") or "")


def _scope(record: Mapping[str, Any]) -> str:
    """What the reads were searched against for this row. See `search_scope`."""
    return str(record.get("search_scope") or "own_rank")


def _is_supported(record: Mapping[str, Any]) -> bool:
    return _status(record) in SUPPORTED


#: Qualifying fragments a trace must have before the organism is named at
#: all. Zero is not a low number, it is no evidence: the alignment touched
#: the reference and then failed an identity, length or coverage gate, so
#: nothing survived to count. Naming an organism on that basis is what put
#: Coccidioides and hookworms on healthy reports — 104 of 145 trace rows
#: across five samples had no qualifying fragment behind them.
MIN_TRACE_FRAGMENTS: Final = 1


def _fragments(record: Mapping[str, Any]) -> int:
    try:
        return int(record.get("unique_supporting_fragments") or 0)
    except (TypeError, ValueError):
        return 0


def _is_candidate(record: Mapping[str, Any]) -> bool:
    """A trace worth naming: candidate status *and* something behind it."""
    return (
        _status(record) in CANDIDATE
        and _fragments(record) >= MIN_TRACE_FRAGMENTS
    )


def _is_bare_nomination(record: Mapping[str, Any]) -> bool:
    """Candidate status with no qualifying fragment at all."""
    return _status(record) in CANDIDATE and _fragments(record) < MIN_TRACE_FRAGMENTS


def _is_ambiguous(record: Mapping[str, Any]) -> bool:
    return _status(record) in AMBIGUOUS


def summarise(pathogens: Mapping[str, Any] | None) -> dict[str, Any]:
    """The counts every part of the report quotes, derived in exactly one place.

    Two separations do the work here.

    What the organism *is*, from the catalogue's `interpretation_class`: only
    recognised disease-causing organisms reach `n_pathogens`, the number the
    page-1 badge shows. An ordinary gut resident and a food yeast are named
    and quantified like everything else, under their own headings, and never
    counted as pathogens.

    How much sequence is *behind* it: `supported` names the organism,
    `candidates` are traces with at least one qualifying fragment, and bare
    nominations with nothing behind them are counted in the ledger but not
    listed as findings.
    """
    records = _records(pathogens)
    supported = [r for r in records if _is_supported(r)]
    candidates = [r for r in records if _is_candidate(r)]
    ambiguous = [r for r in records if _is_ambiguous(r)]
    bare = [r for r in records if _is_bare_nomination(r)]

    # Named findings, split by what the organism is. A row the resolution pass
    # demoted is deliberately `ambiguous_signal` — the sequence is real, the
    # species name is not supportable — so it is named here too. Leaving it out
    # would hide the very rows a reader saw counted as pathogens before.
    demoted = [
        r for r in records
        if str(r.get("report_tier") or "") == "unresolved_complex"
        and _status(r) not in ("not_assessed", "not_detected")
        and not _is_supported(r)
        and not _is_candidate(r)
    ]
    named = supported + candidates + demoted
    by_tier: dict[str, dict[str, Any]] = {}
    for key in TIER_ORDER:
        members = [r for r in named if tier_of(r) == key]
        by_tier[key] = {
            "meta": TIERS[key],
            "supported": [r for r in members if _is_supported(r)],
            "candidates": [r for r in members if _is_candidate(r)],
            "all": members,
            "n": len(members),
        }

    pathogen_tiers = [k for k in TIER_ORDER if TIERS[k]["counts_as_pathogen"]]

    def _counts(record: Mapping[str, Any]) -> bool:
        """One gate for the headline: tier, and the resolution pass's verdict.

        `counts_as_pathogen` is False when the species could not be named from
        the evidence, or when the organism's disease potential is a toxin the
        determinant screen did not find. Both were previously counted, which
        is how a sample whose only Enterobacteriaceae was ordinary E. coli came
        to be reported as five disease-causing organisms.
        """
        if record.get("counts_as_pathogen") is False:
            return False
        return tier_of(record) in pathogen_tiers

    pathogens_supported = [r for r in supported if _counts(r)]
    pathogens_candidates = [r for r in candidates if _counts(r)]

    per_group: dict[str, dict[str, Any]] = {}
    for label, groups, _icon in DISPLAY_GROUPS:
        members = [r for r in records if str(r.get("group", "")) in groups]
        in_group_named = [r for r in members if _is_supported(r) or _is_candidate(r)]
        per_group[label] = {
            # Searched means the reads were aligned against this row's
            # reference. A row that answers only at the rank it shares with a
            # relative was still searched, and counting it as unsearched is
            # what made a complete screen read as two-thirds finished.
            "screened": sum(
                1 for r in members
                if _status(r) != "not_assessed" or _scope(r) in ("group_rank", "masked_by_host")
            ),
            "answerable": sum(
                1 for r in members if _scope(r) != "out_of_assay_scope"
            ),
            "group_rank": sum(1 for r in members if _scope(r) in ("group_rank", "masked_by_host")),
            "out_of_scope": sum(1 for r in members if _scope(r) == "out_of_assay_scope"),
            "pending": sum(1 for r in members if _scope(r) == "reference_pending"),
            "total": len(members),
            "supported": [r for r in members if _is_supported(r)],
            "candidates": [r for r in members if _is_candidate(r)],
            "ambiguous": [r for r in members if _is_ambiguous(r)],
            # What the group card leads with: real pathogens only.
            "pathogens": [r for r in in_group_named if _counts(r)],
            "pathogens_supported": [
                r for r in members if _is_supported(r) and _counts(r)
            ],
            "other_named": [r for r in in_group_named if not _counts(r)],
        }

    ledger = (pathogens or {}).get("coverage") or {}
    dets = _determinants(pathogens)
    return {
        "supported": supported,
        "candidates": candidates,
        "ambiguous": ambiguous,
        "bare_nominations": bare,
        "n_supported": len(supported),
        "n_candidates": len(candidates),
        "n_ambiguous": len(ambiguous),
        "n_bare_nominations": len(bare),
        # The headline: recognised disease-causing organisms named outright.
        "pathogens": pathogens_supported,
        "n_pathogens": len(pathogens_supported),
        "pathogen_candidates": pathogens_candidates,
        "n_pathogen_candidates": len(pathogens_candidates),
        "by_tier": by_tier,
        "n_pathotype_negative": by_tier["pathotype_negative"]["n"],
        "n_unresolved_complex": by_tier["unresolved_complex"]["n"],
        "pathotype_negative": by_tier["pathotype_negative"]["all"],
        "unresolved_complex": by_tier["unresolved_complex"]["all"],
        "per_group": per_group,
        "ledger": ledger,
        "assessed": ledger.get("assessed"),
        "not_assessed": ledger.get("not_assessed"),
        # The honest fraction: rows searched, out of rows this assay can
        # answer at all. An RNA virus has no DNA genome to find, so it is not
        # a gap in the search and does not sit in the denominator.
        "searched": ledger.get("searched") or sum(
            1 for r in records
            if _status(r) != "not_assessed" or _scope(r) in ("group_rank", "masked_by_host")
        ),
        "answerable": ledger.get("answerable") or sum(
            1 for r in records if _scope(r) != "out_of_assay_scope"
        ),
        "searched_group_level": ledger.get("searched_group_level") or sum(
            1 for r in records if _scope(r) == "group_rank"
        ),
        "searched_masked_by_host": ledger.get("searched_masked_by_host") or sum(
            1 for r in records if _scope(r) == "masked_by_host"
        ),
        "masked_records": [r for r in records if _scope(r) == "masked_by_host"],
        # Every record, so a row that answers at a shared rank can name the
        # row whose reference carried the search.
        "records": list(records),
        "out_of_assay_scope": ledger.get("out_of_assay_scope") or sum(
            1 for r in records if _scope(r) == "out_of_assay_scope"
        ),
        "reference_pending": (
            (ledger.get("reference_pending_no_public_sequence") or 0)
            + (ledger.get("reference_pending_needs_member_list") or 0)
        ) or sum(1 for r in records if _scope(r) == "reference_pending"),
        "pending_records": [r for r in records if _scope(r) == "reference_pending"],
        "group_rank_records": [r for r in records if _scope(r) == "group_rank"],
        "determinants": dets,
        "determinants_present": [
            d for d in dets if str(d.get("determinant_status")) in DETERMINANT_PRESENT
        ],
        "ran": bool(records) or bool(ledger),
        "strict_modes": list((pathogens or {}).get("strict_modes") or []),
    }


def _anything_assessed(summary: Mapping[str, Any]) -> bool:
    """Whether at least one target was actually examined.

    The detector can emit a full set of records that are all
    ``not_assessed`` - a missing database, a failed stage, a depth floor
    nothing cleared. Those records make ``ran`` true and produce no hits,
    which is indistinguishable from a clean screen unless the assessed
    count is read. It is read here.
    """
    assessed = summary.get("assessed")
    if isinstance(assessed, (int, float)):
        return assessed > 0
    if isinstance(assessed, (list, tuple, set)):
        return len(assessed) > 0
    # No ledger at all: fall back to whether any record carries a real state.
    records = summary.get("records") or summary.get("pathogens") or []
    return any(
        str(r.get("analysis_status") or "") not in {"", "not_assessed", "not_run"}
        for r in records
        if isinstance(r, Mapping)
    )


def _name(record: Mapping[str, Any]) -> str:
    return str(record.get("display_name") or record.get("target_id") or "unnamed")


def _tone(record: Mapping[str, Any]) -> tuple[colors.Color, colors.Color]:
    """Colour by what the organism is first, then by evidence strength.

    Red is reserved for a named disease-causing organism. Colouring on
    evidence strength alone painted every confidently detected commensal red,
    which is how an abundant, entirely normal gut anaerobe came to look like
    the most alarming thing on the page.
    """
    tier = tier_of(record)
    meta = TIERS[tier]
    if not meta["counts_as_pathogen"]:
        return meta["colour"], meta["background"]
    status = _status(record)
    if status in SUPPORTED:
        return ALERT, ALERT_BG
    if status in CANDIDATE:
        return WATCH, WATCH_BG
    return QUIET, QUIET_BG


# --------------------------------------------------------------------------- #
# Page 1 badge
# --------------------------------------------------------------------------- #


def pathogen_alert(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    pathogens: Mapping[str, Any] | None,
    section: int,
) -> None:
    """The page-1 pathogen anchor: one compact band, linked into the detail.

    This sits at the foot of page 1, and every millimetre it takes is a row
    taken off the disease-signature and metabolite lists above it — the frame
    budget in `pdfsummary` hands those lists whatever is left once this panel
    has been measured. So it earns its height: a verdict, the counts, and a way
    in. Definitions, evidence and confirmation routes belong in the section
    this band links to, not on the summary page.
    """
    summary = summarise(pathogens)
    if not summary["ran"]:
        return

    # The headline counts organisms that passed every gate: the evidence bar,
    # a species the data can actually name, and — where the disease is caused
    # by a toxin rather than by the species — the toxin gene itself. Counting
    # anything less made this badge read "5 disease-causing organisms found"
    # for a sample whose only Enterobacteriaceae was ordinary E. coli.
    n = summary["n_pathogens"]
    n_trace = summary["n_pathogen_candidates"]
    n_opportunist = summary["by_tier"]["opportunist"]["n"]
    n_carriage = summary.get("n_pathotype_negative") or 0
    n_unresolved = summary.get("n_unresolved_complex") or 0
    hit = n > 0

    # Local styles so the band's height is set here rather than inherited from
    # body text: leading is most of what a two-line panel costs.
    head_style = ParagraphStyle(
        "pathhead", parent=st["fine"], fontSize=11.5, leading=13.0,
        spaceBefore=0, spaceAfter=0,
    )
    line_style = ParagraphStyle(
        "pathline", parent=st["fine"], fontSize=7.4, leading=9.0,
        spaceBefore=0, spaceAfter=0,
    )
    # A cell's ALIGN does not move text inside a Paragraph, so the way-in link
    # needs its own right-aligned style to sit against the band's edge instead
    # of floating in the middle of its column.
    link_style = ParagraphStyle(
        "pathlink", parent=line_style, alignment=TA_RIGHT,
    )

    detail: str | None = None
    if hit:
        bug = EvilBug(size=9 * mm, colour=ALERT)
        headline = (
            f"<font size=15 color='{_hex(ALERT)}'><b>{n}</b></font>"
            f"<font size=11.5 color='{_hex(ALERT_DARK)}'><b> disease-causing "
            f"organism{'' if n == 1 else 's'} found</b></font>"
        )
        # The one case where a name has to be on page 1: you should not have to
        # turn the page to learn *what* was found.
        parts = [f"<i>{_name(r)}</i> detected" for r in summary["pathogens"][:3]]
        detail = ", ".join(parts) + (f", and {n - 3} more" if n > 3 else "")
        border, background = ALERT, ALERT_BG
    elif not _anything_assessed(summary):
        # Nothing was actually examined. "None found" would be a negative
        # result, and there is no result: the screen did not run, or every
        # target came back unassessable. A green all-clear here is the most
        # dangerous output the page can produce, so it is refused outright
        # and the reason is printed instead.
        bug = EvilBug(size=8.5 * mm, colour=QUIET, menace=False)
        headline = (
            f"<font size=11.5 color='{_hex(QUIET)}'><b>The pathogen screen did "
            "not run</b></font>"
        )
        detail = (
            "No organism was assessed, so this is not a negative result. "
            "Nothing on this page says an organism is absent."
        )
        border, background = QUIET, QUIET_BG
    else:
        # Nothing met the bar. That is the finding, so the panel says it in
        # those words and wears the calm colour: an amber card over a negative
        # result reads as a warning about nothing. Anything still open is
        # carried by the chips, which is where a count belongs.
        #
        # The headline carries its own scope when coverage is partial. Here
        # 332 of 486 targets were searched, and an unqualified "none found"
        # is a claim about all 486 - it reads as a clean screen of the whole
        # catalogue rather than of the two thirds that were assayable.
        bug = AllClearMark(size=8.5 * mm, colour=CLEAR)
        # The scope lives in the coverage chip rather than here: page one's
        # band has a 15.5 mm budget and any qualifier long enough to be
        # useful wraps the headline onto a second line. The chip carries the
        # exact counts and turns amber when part of the catalogue could not
        # be assessed, so the claim is not left looking like a clean sweep
        # of all 486 targets when 332 were searched.
        headline = (
            f"<font size=11.5 color='{_hex(CLEAR_DARK)}'><b>No disease-causing "
            "organisms found</b></font>"
        )
        border, background = CLEAR, CLEAR_BG

    # Chips are the metrics: scannable, self-measuring, and few enough that the
    # row fits on one line. Only what a reader would act on or judge the result
    # by — anything still open, then the scope it was judged against. Tier
    # tallies that are merely reassuring stay in the section.
    chips: list[Any] = []
    if n_trace:
        chips.append(Chip(f"{n_trace} trace, unconfirmed" if n_trace == 1
                          else f"{n_trace} traces, unconfirmed", WATCH, WATCH_BG))
    if n_unresolved:
        chips.append(Chip(f"{n_unresolved} species unresolved", QUIET, QUIET_BG))
    if n_carriage:
        chips.append(Chip(f"{n_carriage} toxin-negative", QUIET, QUIET_BG))
    if n_opportunist and hit:
        chips.append(Chip(f"{n_opportunist} opportunist" if n_opportunist == 1
                          else f"{n_opportunist} opportunists", QUIET, QUIET_BG))
    # One chip for the scope. It has to answer "was anything skipped?" in the
    # space of a chip, and the earlier version answered wrongly: it read
    # "332 of 486 targets searched" in amber, which says a third of the
    # catalogue was never looked at. It had been counting as unsearched every
    # row that answers at a shared rank (the diarrhoeagenic E. coli
    # pathotypes are one genomic species, told apart by their toxin genes in
    # the determinant screen) and every RNA virus, which a DNA library cannot
    # contain at any depth. Neither is a skipped target. The chip now counts
    # rows searched against rows this assay can answer, and is only amber when
    # a reference is genuinely still pending.
    searched = summary["searched"]
    answerable = summary["answerable"]
    if answerable:
        pending = summary["reference_pending"]
        label = (
            f"all {answerable:,} targets searched" if not pending
            else f"{searched:,} of {answerable:,} targets searched"
        )
        chips.append(Chip(label, WATCH if pending else CLEAR, WATCH_BG if pending else CLEAR_BG))

    text_width = CONTENT_WIDTH - 6 * mm - 14 * mm
    link = RowLink(href=section_dest("pathogens"))

    # Verdict on the left of the top line, the way in on the right of it.
    top = LinkedTable(
        [[
            Paragraph(headline, head_style),
            linked_cell(
                link,
                Paragraph(
                    f"<font size=7.2 color='{_hex(INK_SOFT)}'>Section {section} "
                    "\u2014 full evidence <b>\u2192</b></font>",
                    link_style,
                ),
            ),
        ]],
        colWidths=[text_width * 0.62, text_width * 0.38],
        hAlign="LEFT",
    )
    top.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "BOTTOM"),
                ("ALIGN", (1, 0), (1, 0), "RIGHT"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    text_cell: list[Any] = [top]

    if detail:
        text_cell.extend([
            Spacer(1, 0.8 * mm),
            Paragraph(f"<font color='{_hex(INK)}'>{detail}</font>", line_style),
        ])

    if chips:
        # Wrap to a second line only if the row genuinely cannot fit, rather
        # than letting the last chip run off the page as it did before.
        lines: list[list[Any]] = [[]]
        used = 0.0
        gap = 1.8 * mm
        for chip in chips:
            need = chip.width + gap
            if used + need > text_width and lines[-1]:
                lines.append([])
                used = 0.0
            lines[-1].append(chip)
            used += need
        for i, line in enumerate(lines):
            row = Table(
                [line],
                colWidths=[c.width + gap for c in line],
                hAlign="LEFT",
            )
            row.setStyle(
                TableStyle(
                    [
                        ("LEFTPADDING", (0, 0), (-1, -1), 0),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                        ("TOPPADDING", (0, 0), (-1, -1), 0),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                    ]
                )
            )
            text_cell.extend([Spacer(1, 1.3 * mm if i == 0 else 1.0 * mm), row])

    # The whole band is one link to the pathogen section.
    panel = LinkedTable(
        [[linked_cell(link, bug), text_cell]],
        colWidths=[14 * mm, CONTENT_WIDTH - 14 * mm],
        hAlign="LEFT",
    )
    panel.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("ALIGN", (0, 0), (0, 0), "CENTRE"),
                ("BACKGROUND", (0, 0), (-1, -1), background),
                # Rounded, thin-bordered and undivided, so it sits with the
                # rest of page 1 rather than looking like a warning sticker.
                ("BOX", (0, 0), (-1, -1), 0.6, border),
                ("ROUNDEDCORNERS", [8, 8, 8, 8]),
                ("LEFTPADDING", (0, 0), (0, 0), 2 * mm),
                ("RIGHTPADDING", (0, 0), (0, 0), 1.2 * mm),
                ("LEFTPADDING", (1, 0), (1, 0), 0),
                ("RIGHTPADDING", (1, 0), (1, 0), 3 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.8 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.0 * mm),
            ]
        )
    )
    story.append(panel)


def _hex(c: colors.Color) -> str:
    return f"#{int(c.red * 255):02X}{int(c.green * 255):02X}{int(c.blue * 255):02X}"


# --------------------------------------------------------------------------- #
# How much
# --------------------------------------------------------------------------- #
# "12 pathogens found" answers the wrong question on its own. The engine
# already measures a per-organism quantity — fragments attributed to that
# organism per million analysed — so the report states it, converts it to a
# share of the library that a reader can picture, and says which band it
# falls in. What it must not do is imply a clinical threshold: there is no
# validated cut-off for stool sequence, and inventing one would be worse
# than saying nothing.

#: Bands for normalised fragments per million. These describe the size of
#: the signal in this library, and nothing more. They are not clinical
#: thresholds and carry no infection or severity claim.
#: The words describe how much *sequence* was seen. "moderate" used to sit in
#: the middle band and read as "a moderate infection", so the bands now name
#: the size of the signal and nothing else.
_AMOUNT_BANDS: Final[tuple[tuple[float, str, str], ...]] = (
    (1000.0, "large signal", "a major component of the sample's DNA"),
    (100.0, "clear signal", "well represented in the sample's DNA"),
    (10.0, "small signal", "a minor component of the sample's DNA"),
    (1.0, "faint signal", "a small fraction of the sample's DNA"),
    (0.0, "at the detection floor", "barely above what this depth can see"),
)


def amount_of(record: Mapping[str, Any]) -> dict[str, Any]:
    """The quantity of one organism, in the several ways a reader needs it."""
    fpm = record.get("normalized_fragments_per_million")
    try:
        fpm = None if fpm is None else float(fpm)
    except (TypeError, ValueError):
        fpm = None
    frags = _fragments(record)
    band = word = None
    if fpm is not None and frags:
        for floor, band_name, band_words in _AMOUNT_BANDS:
            if fpm >= floor:
                band, word = band_name, band_words
                break
    percent = None if fpm is None else fpm / 10_000.0
    if fpm is None or not frags:
        figure = "—"
    elif fpm >= 10:
        figure = f"{fpm:,.0f} per million"
    elif fpm >= 1:
        figure = f"{fpm:.1f} per million"
    else:
        figure = f"{fpm:.2f} per million"
    if percent is None or not frags:
        share = ""
    elif percent >= 0.01:
        share = f"{percent:.2f}% of analysed DNA"
    elif percent >= 0.0001:
        share = f"{percent:.4f}% of analysed DNA"
    else:
        share = "well under a thousandth of a percent"
    return {
        "fpm": fpm,
        "fragments": frags,
        "regions": record.get("informative_regions_supported") or 0,
        "percent": percent,
        "figure": figure,
        "share": share,
        "band": band,
        "band_words": word,
    }


def carriage_words(record: Mapping[str, Any]) -> str:
    """Whether finding this organism at all is ordinary, in one clause.

    This is the "is there a normal amount of it?" question, and the honest
    answer comes from the catalogue's stool role rather than from a number
    we do not have a reference distribution for.
    """
    return {
        "carriage_common": (
            "carried by many healthy people, so presence alone is unremarkable"
        ),
        "environmental_or_dietary": (
            "normally arrives with food or from the environment"
        ),
        "intestinal_shedding": (
            "not a normal resident; healthy stool does not usually contain it"
        ),
        "variable_shedding": (
            "shed intermittently, so the amount seen varies between samples"
        ),
        "extraintestinal": (
            "does not normally live in the gut; stool is not its usual site"
        ),
        "tissue_restricted": (
            "lives in tissue rather than the gut lumen, so stool is a poor window on it"
        ),
    }.get(str(record.get("stool_role") or ""), "")


# --------------------------------------------------------------------------- #
# Part A — the six groups, then every finding
# --------------------------------------------------------------------------- #


def _section_heading(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    title: str,
    blurb: str,
) -> None:
    story.extend(section_heading(section, f"{title}", st["h2"]))
    story.append(Paragraph(blurb, st["small"]))
    story.append(Spacer(1, 3 * mm))


def _group_card(
    st: dict[str, ParagraphStyle],
    label: str,
    stats: Mapping[str, Any],
    width: float,
) -> Table:
    """One of the six group tiles: what was found, out of how many searched."""
    # The card counts disease-causing organisms, not everything the screen
    # can name. A group whose only findings are ordinary residents or food
    # organisms reads as clear, because for the question this section asks
    # — is there a pathogen here — it is.
    n_path = len(stats.get("pathogens_supported") or ())
    n_path_trace = len(stats.get("pathogens") or ()) - n_path
    n_other = len(stats.get("other_named") or ())
    n_amb = len(stats["ambiguous"])
    screened = stats["screened"]
    # Rows this assay can answer. The viral tile counted 64 targets and
    # reported 27 searched, which read as a third of a screen; 29 of those 64
    # are RNA viruses with no DNA genome to find.
    answerable = stats.get("answerable") or stats["total"]

    if n_path:
        colour, background = ALERT, ALERT_BG
        figure, verdict = str(n_path), "found"
    elif n_path_trace:
        colour, background = WATCH, WATCH_BG
        figure = str(n_path_trace)
        verdict = "trace only" if n_path_trace == 1 else "traces only"
    elif screened and not n_amb:
        colour, background = CLEAR, CLEAR_BG
        figure, verdict = "none", "found"
    elif screened:
        # Nothing nameable, but sequence did land here. Green would claim a
        # clean group that this evidence does not support.
        colour, background = QUIET, QUIET_BG
        figure, verdict = "none", "nameable"
    else:
        colour, background = QUIET, QUIET_BG
        figure, verdict = "&mdash;", "not assessed"

    big = 16 if figure.isdigit() else 11
    figure_style = ParagraphStyle(
        "gcardfig", parent=st["body"], fontSize=big, leading=big + 1.5,
        textColor=colour, alignment=1,
    )
    label_style = ParagraphStyle(
        "gcardlab", parent=st["label"], fontSize=6.0, leading=7.6, textColor=INK,
        alignment=1,
    )
    verdict_style = ParagraphStyle(
        "gcardv", parent=st["fine"], fontSize=6.2, leading=7.6, textColor=colour,
        alignment=1,
    )
    fine_style = ParagraphStyle(
        "gcardfine", parent=st["fine"], fontSize=5.6, leading=7.0,
        textColor=INK_FAINT, alignment=1,
    )

    extra = ""
    if n_path and n_path_trace:
        extra = f"+{n_path_trace} trace"
    elif n_other:
        extra = f"{n_other} normal or minor"
    elif n_amb and n_path:
        extra = f"+{n_amb} indistinguishable"

    # A fixed four-row grid keeps all six cards the same height whatever
    # each one has to say, so the strip reads as one comparable row.
    rows: list[list[Any]] = [
        [Paragraph(label.upper(), label_style)],
        [Paragraph(f"<b>{figure}</b>", figure_style)],
        [Paragraph(verdict, verdict_style)],
        [
            Paragraph(
                extra
                or (
                    f"all {answerable:,} searched" if answerable and screened >= answerable
                    else f"{screened:,} of {answerable:,} searched" if answerable
                    else "none in this bundle"
                ),
                fine_style,
            )
        ],
    ]
    card = Table(
        rows,
        colWidths=[width],
        rowHeights=[7.6 * mm, 7.4 * mm, 3.6 * mm, 4.4 * mm],
        hAlign="LEFT",
    )
    card.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), background),
                ("BOX", (0, 0), (-1, -1), 0.6, colour),
                ("VALIGN", (0, 0), (0, 0), "MIDDLE"),
                ("VALIGN", (0, 1), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 1.2 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 1.2 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (0, 0), 1.4 * mm),
            ]
        )
    )
    return card


def _finding_row(
    st: dict[str, ParagraphStyle], record: Mapping[str, Any]
) -> list[Any]:
    colour, _bg = _tone(record)
    status = _status(record)
    name_style = ParagraphStyle(
        "fname", parent=st["small"], fontSize=7.6, leading=9.6, textColor=INK,
    )
    status_style = ParagraphStyle(
        "fstat", parent=st["fine"], fontSize=6.4, leading=8.4, textColor=colour,
    )
    grey_style = ParagraphStyle(
        "fgrey", parent=st["fine"], fontSize=6.4, leading=8.4, textColor=INK_SOFT,
    )

    group = str(record.get("group", ""))
    group_label = next(
        (lab for lab, groups, _i in DISPLAY_GROUPS if group in groups), group
    )
    status_text = _STATUS_WORDS.get(status, status.replace("_", " "))
    # The contamination qualifier is deliberately *not* repeated here. It is
    # a property of the run, not of the organism: printing "no extraction or
    # library controls for this sample" under all thirteen rows said nothing
    # thirteen times and read like unfinished output. It is stated once, for
    # the whole section, in `_run_caveats`.
    qualifier = str(record.get("display_qualifier") or "")
    if qualifier and "control" not in qualifier:
        status_text += f"<br/><font color='{_hex(INK_FAINT)}'>{qualifier}</font>"

    amount = amount_of(record)
    if amount["band"]:
        how_much = (
            f"<b>{amount['band']}</b><br/>"
            f"<font color='{_hex(INK_FAINT)}'>{amount['figure']}"
            + (f" · {amount['share']}" if amount["share"] else "")
            + f" · {amount['fragments']:,} fragments in "
            f"{amount['regions']:,} regions</font>"
        )
    else:
        how_much = (
            f"&mdash;<br/><font color='{_hex(INK_FAINT)}'>no qualifying "
            "sequence to measure</font>"
        )

    # Breadth is the statistic that separates a real organism from reads
    # landing on a conserved slice of its genome, and a reader cannot judge a
    # finding without it. A whole genome covered at 0.5% is not an organism
    # "clearly present", however many fragments piled onto it.
    breadth = record.get("reference_breadth_fraction")
    identity = record.get("median_alignment_identity")
    bits = []
    if isinstance(breadth, (int, float)):
        bits.append(f"{float(breadth) * 100:.2f}% of its genome covered")
    if isinstance(identity, (int, float)):
        bits.append(f"{float(identity) * 100:.1f}% identity")
    evidence = f"<font color='{_hex(INK_FAINT)}'>{' · '.join(bits)}</font>" if bits else "&mdash;"

    marker_words = _marker_words(record)
    if marker_words:
        evidence += f"<br/>{marker_words}"

    return [
        Paragraph(f"<b>{_name(record)}</b>", name_style),
        Paragraph(group_label, grey_style),
        Paragraph(status_text, status_style),
        Paragraph(how_much, grey_style),
        Paragraph(evidence, grey_style),
    ]


def _marker_words(record: Mapping[str, Any]) -> str:
    """Toxin/pathotype status in one phrase, coloured by what it means.

    For an organism whose disease potential *is* the gene, this is the finding.
    Burying it under the species name is how "C. difficile found" came to be
    printed for a sample with no toxin genes at all.
    """
    evidence = str(record.get("pathotype_evidence") or "not_applicable")
    markers = [str(m).split(".", 1)[-1] for m in record.get("pathotype_markers") or ()]
    if evidence == "supported":
        return (
            f"<font color='{_hex(ALERT_DARK)}'><b>disease genes present"
            + (f": {', '.join(markers)}" if markers else "")
            + "</b></font>"
        )
    if evidence == "not_detected":
        return f"<font color='{_hex(CLEAR_DARK)}'>disease genes not found</font>"
    if evidence == "not_assessed":
        return f"<font color='{_hex(INK_FAINT)}'>disease genes not assessed</font>"
    return ""


_FINDING_HEAD: Final = (
    "Organism", "Group", "What the sequence shows", "How much", "Strength of evidence",
)
_FINDING_WIDTHS: Final = (0.23, 0.12, 0.19, 0.24, 0.22)


def _finding_table(
    st: dict[str, ParagraphStyle],
    records: Sequence[Mapping[str, Any]],
    *,
    with_cards: frozenset[str] | set[str] = frozenset(),
    anchored: bool = False,
) -> Table:
    """The findings as rows.

    A row whose organism has a detail card (`with_cards`) is a link to it;
    with `anchored`, each row is also the destination its card links back
    to. The bare-nominations table sets neither.
    """
    head_style = ParagraphStyle(
        "fhead", parent=st["label"], fontSize=6.0, leading=8, textColor=INK_SOFT,
    )
    data: list[list[Any]] = [[Paragraph(h.upper(), head_style) for h in _FINDING_HEAD]]
    for record in records:
        cells = _finding_row(st, record)
        target = str(record.get("target_id") or "")
        marker = RowLink(
            href=pathogen_dest(target, detail=True) if target in with_cards else None,
            anchor=pathogen_dest(target) if anchored and target else None,
        )
        if marker.href or marker.anchor:
            cells[0] = linked_cell(marker, cells[0])
        data.append(cells)

    table = LinkedTable(
        data,
        colWidths=[CONTENT_WIDTH * w for w in _FINDING_WIDTHS],
        hAlign="LEFT",
        repeatRows=1,
    )
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1.4 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 1.5 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5 * mm),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK_FAINT),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
    ]
    for i, record in enumerate(records, start=1):
        colour, background = _tone(record)
        style.append(("BACKGROUND", (0, i), (-1, i), background))
        style.append(("LINEBEFORE", (0, i), (0, i), 1.4, colour))
    table.setStyle(TableStyle(style))
    return table


def pathogens_glance(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    pathogens: Mapping[str, Any] | None,
    section: int,
    detail_section: int,
) -> None:
    """Part A: the six groups as tiles, then every finding worth naming."""
    summary = summarise(pathogens)
    if not summary["ran"]:
        _section_heading(
            story, st, section=section, title="Pathogens",
            blurb="The pathogen branch did not run for this sample.",
        )
        return

    _section_heading(
        story, st, section=section, title="Pathogens",
        blurb=(
            "A direct search of your sequencing reads for known disease-causing "
            "organisms — bacteria, protozoa, parasitic worms, microsporidia, fungi "
            "and viruses — against curated reference genomes. This is a "
            "<b>sequence</b> question, not a diagnosis: it reports what DNA is "
            "present, at what resolution, and how strong the evidence is. "
            f"Section {detail_section} takes every finding one at a time, and lists "
            "what could <i>not</i> be assessed in this sample."
        ),
    )
    # ---- the six group tiles ---------------------------------------------- #
    gap = 2.5 * mm
    card_width = (CONTENT_WIDTH - 5 * gap) / 6
    cards = [
        _group_card(st, label, summary["per_group"][label], card_width)
        for label, _groups, _icon in DISPLAY_GROUPS
    ]
    grid = Table(
        [cards],
        colWidths=[card_width + gap] * 5 + [card_width],
        hAlign="LEFT",
    )
    grid.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), gap),
                ("RIGHTPADDING", (-1, 0), (-1, 0), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    story.append(grid)
    story.append(Spacer(1, 4 * mm))

    # ---- the findings, grouped by what the organism is -------------------- #
    # One flat list headed "pathogens found" was the single worst thing about
    # the old version of this page: it put brewer's yeast and an ordinary gut
    # anaerobe in the same red table as Shigella. Each tier now carries its
    # own heading and its own plain statement of what a finding there means.
    # The organisms that get a card in the detail section: the same lists
    # that section iterates, so a row links onward exactly when there is
    # somewhere to go.
    with_cards = {
        str(r.get("target_id") or "")
        for stats in summary["per_group"].values()
        for r in stats["supported"] + stats["candidates"]
    }
    any_named = False
    for key in TIER_ORDER:
        tier = summary["by_tier"][key]
        if not tier["n"]:
            continue
        any_named = True
        meta = tier["meta"]
        n_sup = len(tier["supported"])
        n_trace = len(tier["candidates"])
        counted = []
        if n_sup:
            counted.append(f"{n_sup} named")
        if n_trace:
            counted.append(f"{n_trace} trace")
        story.append(Spacer(1, 3 * mm))
        story.append(
            Paragraph(
                f"<font color='{_hex(meta['colour'])}'>{meta['heading']}</font>"
                f" &mdash; <font size=8>{', '.join(counted)}</font>",
                st["h3"],
            )
        )
        story.append(Paragraph(meta["blurb"], st["fine"]))
        story.append(Spacer(1, 1.5 * mm))
        rows = sorted(
            tier["all"], key=lambda r: (0 if _is_supported(r) else 1, _name(r))
        )
        story.append(_finding_table(st, rows, with_cards=with_cards, anchored=True))

    if not any_named:
        story.append(
            Paragraph(
                "No organism in the panel returned sequence strong enough to "
                "name in this sample.",
                st["body"],
            )
        )
    else:
        story.append(Spacer(1, 2.5 * mm))
        story.append(
            Paragraph(
                "<b>How to read the two strengths.</b> <b>Enough sequence to "
                "name it</b> means distinct, organism-specific DNA passed every "
                "gate: the organism really is in the sample. <b>Trace</b> means "
                "matching DNA was found but too little of it to be sure of the "
                "species &mdash; it is neither a finding nor a clean negative, "
                "and at these amounts a close relative or a small amount of "
                "carry-over would look the same. <b>How much</b> is that "
                "organism's share of all the DNA analysed here; the bands "
                "describe the size of the signal in this sample and are not "
                "clinical thresholds, because no validated cut-off exists for "
                "stool sequencing.",
                st["fine"],
            )
        )

    if summary["n_ambiguous"]:
        story.append(Spacer(1, 3 * mm))
        story.append(
            Paragraph(
                f"<b>{summary['n_ambiguous']} further target"
                f"{'s' if summary['n_ambiguous'] != 1 else ''}</b> matched only "
                "sequence shared with relatives, so no organism can be named "
                f"from it. Listed in section {detail_section}.",
                st["fine"],
            )
        )


# --------------------------------------------------------------------------- #
# Part B — one card per finding
# --------------------------------------------------------------------------- #


def _advice_block(
    st: dict[str, ParagraphStyle],
    record: Mapping[str, Any],
    advice: Any | None,
    colour: colors.Color,
) -> list[Any]:
    """What it is, whether it is normal, how much, and what acts against it.

    This is the half of the answer the section used to leave out. Naming an
    organism and stopping there tells a reader enough to be frightened and
    not enough to do anything.
    """
    body = ParagraphStyle(
        "adv", parent=st["fine"], fontSize=7.0, leading=9.4, textColor=INK,
    )
    head = ParagraphStyle(
        "advh", parent=st["label"], fontSize=6.2, leading=8.2, textColor=colour,
    )
    flow: list[Any] = []
    amount = amount_of(record)

    def para(title: str, text: str) -> None:
        flow.append(Paragraph(title.upper(), head))
        flow.append(Paragraph(text, body))
        flow.append(Spacer(1, 1.4 * mm))

    strict = str(record.get("strict_mode") or "") == "strict_cdiff"
    if strict:
        # The strict rule was asked for. Its statement leads the card, in
        # the alert colour, and the ordinary prose follows it.
        flow.append(Paragraph("WHAT THIS TEST FOUND", head))
        flow.append(Paragraph(
            f"<font color='{_hex(ALERT_DARK)}'><b>Positive for <i>Clostridioides difficile</i>.</b></font> "
            + str(record.get("plain_statement") or ""),
            body))
        flow.append(Spacer(1, 1.4 * mm))
    if advice is not None:
        para("What it is", advice.what_it_is)
        para("Why it matters", advice.why_it_matters)

    # How much, always — this is measured, so it is stated whether or not
    # any written advice exists for the organism.
    if amount["band"]:
        how = (
            f"<b>{amount['band'].capitalize()}</b> in this sample: "
            f"{amount['figure']}"
            + (f", {amount['share']}" if amount["share"] else "")
            + f", from {amount['fragments']:,} distinct DNA fragments across "
            f"{amount['regions']:,} separate regions of its genome. "
            f"{amount['band_words'].capitalize()}."
        )
    else:
        how = (
            "No qualifying sequence survived the evidence gates, so there is "
            "no amount to report. The organism is listed because matching "
            "sequence appeared, not because a quantity was measured."
        )
    carriage = carriage_words(record)
    if carriage:
        how += f" For this organism, it is {carriage}."
    breadth = record.get("reference_breadth_fraction")
    if isinstance(breadth, (int, float)):
        how += (
            f" The aligned fragments cover <b>{float(breadth) * 100:.2f}%</b> of this "
            "organism's reference genome. A genuinely present organism is covered broadly; a "
            "small percentage means the reads sit on a slice of the genome that relatives may "
            "share, which is why breadth is reported next to the count."
        )
    para("How much is here", how)

    # "Is this normal?" — answered with a published rate where one exists,
    # because a reader who is not told assumes the worst.
    published = str(record.get("carriage_statement") or "")
    if published:
        source = str(record.get("carriage_source") or "")
        para(
            "Is carrying this normal",
            published + (f" <font color='{_hex(INK_FAINT)}'>Source: {source}</font>" if source else ""),
        )

    # For a toxin- or pathotype-defined organism, the gene is the finding.
    evidence = str(record.get("pathotype_evidence") or "not_applicable")
    if evidence != "not_applicable":
        markers = [str(m).split(".", 1)[-1] for m in record.get("pathotype_markers") or ()]
        if evidence == "supported":
            text = (
                "The defining disease genes <b>were found</b>"
                + (f" ({', '.join(markers)})" if markers else "")
                + ". That is what makes this organism a finding rather than carriage, and it is "
                "the part worth showing a doctor."
            )
        elif evidence == "not_detected" and strict:
            text = (
                "The toxin genes <b>were looked for and not found</b> among the fragments in "
                "this sample. That does not clear the organism: toxin genes are a small target "
                "that shotgun sequencing can miss at ordinary depth, so their absence here "
                "lowers but does not remove the possibility of the toxin-producing form."
            )
        elif evidence == "not_detected":
            text = (
                "The defining disease genes <b>were looked for and not found</b>. Without them "
                "this organism is not the disease-causing form, so this is carriage evidence "
                "rather than a disease finding."
            )
        else:
            text = (
                "The defining disease genes were <b>not assessed</b> in this sample, so neither "
                "carriage nor the disease-causing form can be claimed."
            )
        para("The genes that decide whether it matters", text)

    # A species the data cannot name says so here, with the test that settles it.
    resolution_state = str(record.get("species_resolution") or "resolved")
    if resolution_state in ("not_resolvable_within_complex", "group_level_only"):
        label = str(record.get("complex_label") or "this complex")
        para(
            "Why no species name is given",
            f"This organism belongs to the {label}, whose members are the same genomic species "
            "as a relative that lives in most healthy guts. Alignment cannot separate them, so "
            "naming a species here would be a guess. A clinical stool culture or a nucleic-acid "
            "test on a fresh sample is what resolves it.",
        )

    if advice is not None:
        para("What would make the amount concerning", advice.overgrowth_signal)
        if advice.has_agents:
            flow.append(Paragraph("WHAT IS KNOWN TO ACT AGAINST IT", head))
            flow.append(Paragraph(
                "Non-prescription agents with published activity against this "
                "organism. The tier after each name is how strong the evidence "
                "is, and it matters: a culture result is not a demonstration "
                "that swallowing the substance changes anything in a person. "
                "None of this is a prescription or a treatment plan.",
                ParagraphStyle("advi", parent=st["fine"], fontSize=6.2,
                               leading=8.2, textColor=INK_FAINT),
            ))
            flow.append(Spacer(1, 1.2 * mm))
            for agent in advice.agents:
                tier_colour = CLEAR if agent.is_human_evidence else QUIET
                flow.append(Paragraph(
                    f"<b>{agent.name}</b> "
                    f"<font size=6 color='{_hex(INK_FAINT)}'>({agent.kind})</font> "
                    f"&mdash; <font size=6 color='{_hex(tier_colour)}'>"
                    f"{agent.tier_words}</font>",
                    body,
                ))
                detail = (
                    f"{agent.effect} <i>Dose studied:</i> {agent.dose_studied}. "
                    f"<i>Reaching the gut:</i> {agent.reaches_gut} "
                    f"<i>Caution:</i> {agent.caution} "
                    f"<font color='{_hex(INK_FAINT)}'>"
                    f"<link href='{agent.citation_url}'>{agent.citation_title}</link>"
                    "</font>"
                )
                flow.append(Paragraph(
                    detail,
                    ParagraphStyle("advd", parent=st["fine"], fontSize=6.4,
                                   leading=8.6, textColor=INK_SOFT,
                                   leftIndent=3 * mm),
                ))
                flow.append(Spacer(1, 1.0 * mm))
        elif advice.no_agents_reason:
            para("What to do about it", advice.no_agents_reason)
    return flow


def _detail_card(
    st: dict[str, ParagraphStyle],
    record: Mapping[str, Any],
    advice: Any | None = None,
) -> list[Any]:
    """One finding in full: what it is, what the evidence is, what to do."""
    colour, background = _tone(record)
    status = _status(record)

    title_style = ParagraphStyle(
        "dtitle", parent=st["h3"], fontSize=9.4, leading=12, textColor=INK,
    )
    badge_style = ParagraphStyle(
        "dbadge", parent=st["label"], fontSize=6.4, leading=8.4, textColor=colour,
        alignment=2,
    )

    target = str(record.get("target_id") or "")
    header = Table(
        [[
            Paragraph(
                f"{_name(record)}",
                title_style,
            ),
            Paragraph(
                _STATUS_WORDS.get(status, status.replace("_", " ")).upper(),
                badge_style,
            ),
        ]],
        colWidths=[CONTENT_WIDTH * 0.66, CONTENT_WIDTH * 0.34],
        hAlign="LEFT",
    )
    header.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "BOTTOM"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2),
                ("LINEBELOW", (0, 0), (-1, 0), 0.7, colour),
            ]
        )
    )

    flow: list[Any] = [
        Anchor(pathogen_dest(target, detail=True), above=6 * mm),
        back_link_line(pathogen_dest(target), "at a glance"),
        header,
        Spacer(1, 1.6 * mm),
    ]

    # The reader's questions first — what is it, is it normal, how much is
    # there, what can be done — then the technical evidence underneath.
    flow.extend(_advice_block(st, record, advice, colour))

    statement = str(record.get("clinical_interpretation") or "").strip()
    if statement:
        flow.append(Paragraph(statement.replace("_", " "), st["small"]))
        flow.append(Spacer(1, 1.6 * mm))

    # Evidence, as orthogonal facts rather than one blended score.
    facts: list[tuple[str, str]] = [
        ("Resolved to", _RESOLUTION_WORDS.get(
            str(record.get("resolution", "")), str(record.get("resolution", "")) or "—"
        )),
        ("Organism class", _CLASS_WORDS.get(
            str(record.get("interpretation_class", "")),
            str(record.get("interpretation_class", "")) or "—",
        )),
        ("Reference route", str(record.get("reference_status", "—")).replace("_", " ")),
    ]
    frags = record.get("unique_supporting_fragments") or 0
    if frags:
        facts.append(("Distinct fragments", f"{frags:,}"))
        facts.append((
            "Informative regions",
            f"{record.get('informative_regions_supported') or 0:,}",
        ))
        facts.append((
            "Informative bases",
            f"{record.get('informative_bases_covered') or 0:,}",
        ))
    if record.get("ambiguous_fragments"):
        facts.append((
            "Shared-sequence fragments",
            f"{record['ambiguous_fragments']:,} (not counted as support)",
        ))
    identity = record.get("median_alignment_identity")
    if identity is not None:
        facts.append(("Median identity", f"{identity * 100:.1f}%"))
    breadth = record.get("informative_region_breadth_fraction")
    if breadth is not None:
        facts.append(("Informative breadth", f"{breadth * 100:.2f}%"))

    label_style = ParagraphStyle(
        "dlab", parent=st["fine"], fontSize=6.2, leading=8.2, textColor=INK_FAINT,
    )
    value_style = ParagraphStyle(
        "dval", parent=st["fine"], fontSize=6.6, leading=8.6, textColor=INK,
    )
    pairs: list[list[Any]] = []
    for i in range(0, len(facts), 3):
        chunk = facts[i : i + 3]
        row: list[Any] = []
        for name, value in chunk:
            row.append(
                [
                    Paragraph(name.upper(), label_style),
                    Paragraph(value, value_style),
                ]
            )
        while len(row) < 3:
            row.append("")
        pairs.append(row)
    grid = Table(pairs, colWidths=[CONTENT_WIDTH / 3] * 3, hAlign="LEFT")
    grid.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 0.8 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0.8 * mm),
            ]
        )
    )
    flow.append(grid)

    linked = record.get("linked_determinants") or ()
    unlinked = record.get("unlinked_determinants") or ()
    if linked or unlinked:
        bits = []
        if linked:
            bits.append(
                f"<b>Toxin or resistance genes on the same DNA as this organism:</b> "
                f"{', '.join(str(x) for x in linked)}"
            )
        if unlinked:
            bits.append(
                f"<b>Present in the sample but not linked to this organism:</b> "
                f"{', '.join(str(x) for x in unlinked)}"
            )
        flow.append(Spacer(1, 1.2 * mm))
        flow.append(Paragraph(" &nbsp;·&nbsp; ".join(bits), st["fine"]))

    options = record.get("confirmation_options") or ()
    if options:
        flow.append(Spacer(1, 1.2 * mm))
        flow.append(
            Paragraph(
                "<b>What would confirm this:</b> "
                + ", ".join(str(o).replace("_", " ") for o in options)
                + ".",
                st["fine"],
            )
        )

    qualifier = record.get("display_qualifier")
    if qualifier:
        flow.append(Spacer(1, 1.2 * mm))
        flow.append(Paragraph(f"<b>Caveat:</b> {qualifier}.", st["fine"]))

    reasons = record.get("reason_codes") or ()
    if reasons:
        flow.append(Spacer(1, 1 * mm))
        flow.append(
            Paragraph(
                "<font color='"
                + _hex(INK_FAINT)
                + "'>"
                + " · ".join(str(r).replace("_", " ") for r in reasons)
                + "</font>",
                st["fine"],
            )
        )

    # One flowable per row, not one tall cell. ReportLab cannot split a
    # single table cell, so a card taller than the page is a hard
    # LayoutError — which is exactly what happened once these cards started
    # carrying what the organism is and what acts against it. Splitting by
    # row keeps the tinted panel and the coloured spine continuous while
    # letting a long card run onto the next page.
    body = Table(
        [[item] for item in flow], colWidths=[CONTENT_WIDTH], hAlign="LEFT"
    )
    last = len(flow) - 1
    body.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), background),
                ("LINEBEFORE", (0, 0), (0, -1), 1.8, colour),
                ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, 0), 2.4 * mm),
                ("BOTTOMPADDING", (0, last), (-1, last), 2.4 * mm),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return [body, Spacer(1, 3 * mm)]


def pathogens_detail_pages(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    pathogens: Mapping[str, Any] | None,
    section: int,
    agents: Any | None = None,
) -> None:
    """Part B: findings group by group, then determinants, then the ledger.

    `agents` is the optional advice library. When it is absent every card
    still reports what it measured; when it is present each card also says
    what the organism is, whether carrying it is ordinary, and what has
    published activity against it.
    """
    summary = summarise(pathogens)
    if not summary["ran"]:
        return

    def advice_for(record: Mapping[str, Any]) -> Any | None:
        if agents is None:
            return None
        return agents.get(str(record.get("target_id") or ""))

    _section_heading(
        story, st, section=section, title="Pathogens in detail",
        blurb=(
            "Each finding on its own, by group, with the evidence behind it kept "
            "as separate facts rather than one blended score: how much distinct "
            "sequence supported it, how much of the reference that sequence "
            "covered, what resolution it justifies, and what would confirm it. "
            "Groups with nothing found still appear, because a negative from a "
            "search that was capable of finding something is itself a result."
        ),
    )

    for label, _groups, _icon in DISPLAY_GROUPS:
        stats = summary["per_group"][label]
        story.append(CondPageBreak(45 * mm))
        story.append(Spacer(1, 1 * mm))
        story.append(Paragraph(label, st["h3"]))

        members = stats["supported"] + stats["candidates"]
        if members:
            story.append(
                Paragraph(
                    f"{len(stats['supported'])} supported, "
                    f"{len(stats['candidates'])} candidate, out of "
                    f"{stats['screened']:,} targets searched in this group.",
                    st["fine"],
                )
            )
            story.append(Spacer(1, 2 * mm))
            for record in sorted(
                members, key=lambda r: (0 if _is_supported(r) else 1, _name(r))
            ):
                story.extend(_detail_card(st, record, advice_for(record)))
        else:
            answerable = stats.get("answerable") or stats["total"]
            scope = (
                f"all {answerable:,} targets in this group"
                if stats["screened"] >= answerable
                else f"{stats['screened']:,} of {answerable:,} targets in this group"
            )
            note = (
                f"Nothing found. {scope} were searched and none returned supporting or "
                "candidate sequence."
                if stats["screened"]
                else "No target in this group could be assessed in this sample."
            )
            if stats.get("out_of_scope"):
                note += (
                    f" A further {stats['out_of_scope']:,} are RNA viruses, which carry no DNA "
                    "for a DNA test to find."
                )
            if stats["ambiguous"]:
                note += (
                    f" {len(stats['ambiguous'])} matched only sequence shared with "
                    "relatives, which cannot name an organism."
                )
            story.append(Paragraph(note, st["small"]))
            story.append(Spacer(1, 2.5 * mm))

    _ambiguous_block(story, st, summary)
    _determinant_block(story, st, summary)
    _ledger_block(story, st, summary)


def _ambiguous_block(
    story: list[Any], st: dict[str, ParagraphStyle], summary: Mapping[str, Any]
) -> None:
    records = summary["ambiguous"]
    if not records:
        return
    story.append(CondPageBreak(40 * mm))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("Matched only shared sequence", st["h3"]))
    story.append(
        Paragraph(
            "These targets picked up reads, but every read fell in sequence they "
            "share with close relatives. That is not evidence for the named "
            "organism and it is not evidence against it. They are listed so the "
            "signal is not silently discarded.",
            st["fine"],
        )
    )
    story.append(Spacer(1, 2 * mm))
    story.append(_finding_table(st, sorted(records, key=_name)))
    story.append(Spacer(1, 3 * mm))


def _determinant_block(
    story: list[Any], st: dict[str, ParagraphStyle], summary: Mapping[str, Any]
) -> None:
    dets = summary["determinants"]
    if not dets:
        return
    present = summary["determinants_present"]
    story.append(CondPageBreak(50 * mm))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("Toxin, virulence and resistance genes", st["h3"]))
    story.append(
        Paragraph(
            "Whether a gene is present is asked separately from which organism "
            "carries it. A resistance gene found in stool belongs to <i>some</i> "
            "organism in the community, and unless the read or contig physically "
            "links it to a named genome, that organism is unknown. Gene presence "
            "is also not a prediction of treatment failure.",
            st["fine"],
        )
    )
    story.append(Spacer(1, 2 * mm))

    if not present:
        story.append(
            Paragraph(
                f"None of the {len(dets):,} screened determinants returned "
                "supporting sequence.",
                st["small"],
            )
        )
        story.append(Spacer(1, 3 * mm))
        return

    head_style = ParagraphStyle(
        "dethead", parent=st["label"], fontSize=6.0, leading=8, textColor=INK_SOFT,
    )
    name_style = ParagraphStyle(
        "detname", parent=st["small"], fontSize=7.2, leading=9.2, textColor=INK,
    )
    grey_style = ParagraphStyle(
        "detgrey", parent=st["fine"], fontSize=6.4, leading=8.4, textColor=INK_SOFT,
    )
    data: list[list[Any]] = [
        [
            Paragraph(h, head_style)
            for h in ("GENE", "KIND", "WHAT WAS FOUND", "CARRIER ORGANISM", "SUPPORT")
        ]
    ]
    linkage_words = {
        "contig_supported": "linked on the same DNA fragment",
        "read_pair_supported": "linked on the same read pair",
        "genome_supported": "linked to an assembled genome",
        "ambiguous": "linkage ambiguous",
        "unlinked": "carrier not established",
    }
    status_words = {
        "supported_intact_sequence": "present, complete",
        "supported_partial_sequence": "present, partial locus",
        "ambiguous_allele": "family present, variant unresolved",
    }
    for det in sorted(present, key=lambda d: str(d.get("display_name", ""))):
        kind = str(det.get("kind", ""))
        amr_class = det.get("amr_class")
        kind_text = kind.replace("_", " ")
        if amr_class:
            kind_text += f"<br/><font size=5.6>{str(amr_class).replace('_', ' ')}</font>"
        det_status = str(det.get("determinant_status", ""))
        linkage = str(det.get("host_linkage", ""))
        carriers = det.get("linked_target_ids") or ()
        carrier_text = (
            ", ".join(str(c) for c in carriers)
            if carriers
            else linkage_words.get(linkage, linkage.replace("_", " "))
        )
        data.append([
            Paragraph(f"<b>{det.get('display_name', det.get('determinant_id'))}</b>", name_style),
            Paragraph(kind_text, grey_style),
            Paragraph(
                status_words.get(det_status, det_status.replace("_", " ")), grey_style
            ),
            Paragraph(carrier_text, grey_style),
            Paragraph(f"{det.get('supporting_fragments') or 0:,} fragments", grey_style),
        ])
    table = Table(
        data,
        colWidths=[CONTENT_WIDTH * w for w in (0.24, 0.18, 0.26, 0.16, 0.16)],
        hAlign="LEFT",
        repeatRows=1,
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK_FAINT),
                ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
                ("BACKGROUND", (0, 1), (-1, -1), PANEL_BG),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 3 * mm))


def _ledger_block(
    story: list[Any], st: dict[str, ParagraphStyle], summary: Mapping[str, Any]
) -> None:
    """What was not assessed and why — the honest boundary of the screen."""
    ledger = summary["ledger"]
    if not ledger:
        return
    story.append(CondPageBreak(55 * mm))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("What this screen did and did not cover", st["h3"]))
    story.append(
        Paragraph(
            "Every number below is counted from the run itself, not from a fixed "
            "marketing total. Organisms, determinants and subtypes are counted "
            "separately and never summed into a single figure, because they are "
            "different kinds of question.",
            st["fine"],
        )
    )
    story.append(Spacer(1, 2 * mm))

    # Every line names what it is rather than what it is not. The old table
    # had three "Not assessed" rows totalling 154 against 332 assessed, which
    # reads as a screen that gave up on a third of its own catalogue.
    total = ledger.get("total_targets") or (summary["assessed"] + summary["not_assessed"])
    rows: list[tuple[str, str, str]] = [
        (
            "Targets in the catalogue",
            f"{total:,}",
            "every organism this screen carries a definition for",
        ),
        (
            "Searched, answered at its own rank",
            f"{summary['assessed']:,}",
            "own reference sequence, aligned, and separable from its relatives",
        ),
    ]
    if summary["searched_group_level"]:
        rows.append((
            "Searched, answers for a group",
            f"{summary['searched_group_level']:,}",
            "whole genome shared with a catalogue relative; the toxin and virulence "
            "gene screen below is what separates these",
        ))
    if summary["searched_masked_by_host"]:
        rows.append((
            "Searched, masked by the host",
            f"{summary['searched_masked_by_host']:,}",
            "sequence also present in the human reference genome, so a matching read cannot be "
            "assigned to the organism rather than to you",
        ))
    if summary["out_of_assay_scope"]:
        rows.append((
            "RNA genome, outside a DNA test",
            f"{summary['out_of_assay_scope']:,}",
            "no DNA exists in the sample to find; needs a stool RNA panel",
        ))
    if summary["reference_pending"]:
        rows.append((
            "Reference not available yet",
            f"{summary['reference_pending']:,}",
            "no public genome or marker sequence good enough to search; named below",
        ))
    if ledger.get("not_assessed_route"):
        rows.append((
            "Analysis step failed",
            f"{ledger['not_assessed_route']:,}",
            "a stage did not complete for these targets",
        ))
    placed = (
        summary["assessed"] + summary["searched_group_level"] + summary["searched_masked_by_host"]
        + summary["out_of_assay_scope"] + summary["reference_pending"]
        + (ledger.get("not_assessed_route") or 0)
    )
    if placed < total:
        rows.append((
            "Not placed in a category above",
            f"{total - placed:,}",
            "the reason recorded for these does not fall in any category above; each one's own "
            "record carries it",
        ))
    rows.extend([
        (
            "Determinants searched",
            f"{ledger.get('determinants_assessed', 0):,}",
            "toxin, virulence and resistance genes, searched independently of species",
        ),
        (
            "Determinants not searched",
            f"{ledger.get('determinants_not_assessed', 0):,}",
            "no usable reference for the gene family",
        ),
    ])
    label_style = ParagraphStyle(
        "lgl", parent=st["fine"], fontSize=6.8, leading=9, textColor=INK,
    )
    num_style = ParagraphStyle(
        "lgn", parent=st["small"], fontSize=8.2, leading=10, textColor=INK,
        alignment=2,
    )
    note_style = ParagraphStyle(
        "lgd", parent=st["fine"], fontSize=6.2, leading=8.4, textColor=INK_SOFT,
    )
    data = [
        [
            Paragraph(f"<b>{label}</b>", label_style),
            Paragraph(value, num_style),
            Paragraph(note, note_style),
        ]
        for label, value, note in rows
    ]
    table = Table(
        data,
        colWidths=[CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.10, CONTENT_WIDTH * 0.60],
        hAlign="LEFT",
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("RIGHTPADDING", (0, 0), (-1, -1), 1.4 * mm),
                ("TOPPADDING", (0, 0), (-1, -1), 1.3 * mm),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1.3 * mm),
                ("LINEBELOW", (0, 0), (-1, -2), 0.25, RULE),
                ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
            ]
        )
    )
    story.append(table)

    # Name them. A count of what is missing is a claim a reader cannot check;
    # a list is one they can.
    pending = summary.get("pending_records") or []
    if pending:
        names = ", ".join(sorted(f"<i>{_name(r)}</i>" for r in pending))
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<b>The {len(pending):,} targets with no reference yet.</b> {names}. Each is in the "
            "catalogue and each will be searched as soon as a reference exists for it. Until "
            "then this report says nothing either way about them.", st["fine"]))
    grouped = summary.get("group_rank_records") or []
    if grouped:
        by_id = {str(r.get("target_id")): r for r in (summary.get("records") or [])}

        def with_holder(record: Mapping[str, Any]) -> str:
            holder = by_id.get(str(record.get("covered_by_target_id") or ""))
            if holder is None:
                return f"<i>{_name(record)}</i>"
            return f"<i>{_name(record)}</i> (searched as <i>{_name(holder)}</i>)"

        names = ", ".join(sorted(with_holder(r) for r in grouped))
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<b>The {len(grouped):,} targets that answer for a group.</b> {names}. These were "
            "searched. Their genomes are the same as a relative's in this catalogue, so a read "
            "cannot say which of the pair it came from; where the difference between them is a "
            "gene rather than a genome &mdash; the diarrhoea-causing <i>E. coli</i> types, "
            "toxin-producing <i>Bacteroides fragilis</i> &mdash; the gene screen above answers it.",
            st["fine"]))

    story.append(Spacer(1, 2 * mm))
    story.append(
        Paragraph(
            "<b>Limits that apply to every finding above.</b> This is a research "
            "screen on a stool sample, not a diagnostic test. Sequence presence is "
            "not infection, and it is not a measure of how much organism is there: "
            "read counts depend on how much of your sample was that organism's DNA, "
            "not on its abundance in you. A negative here does not rule out an "
            "organism that sheds intermittently, sits in tissue rather than the "
            "gut lumen, or has an RNA genome that a DNA library cannot see. "
            "Nothing here should change treatment without a clinician and, where "
            "listed, a confirmatory test.",
            st["fine"],
        )
    )
    story.append(Spacer(1, 2 * mm))
