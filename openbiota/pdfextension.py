"""Rendering for the v0.8.3 additions.

BUILD_SPEC_v0.8.3 §12.1 places each new reading with the existing
readings it belongs beside, rather than in a chapter of its own: the ecology
companions go into the gut-community section, the resistance classes into the
pathogen detail, and the simulation into the actions and metabolism sections.

The display rules of §12.3, applied here rather than described:
every quantitative graphic shows its value, unit, what the scale means,
coverage and status; a reference band appears only where a real reference
exists; and high, low and unavailable are distinguished by text and shape,
never by colour alone.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, Spacer, Table, TableStyle

from openbiota.pdflinks import Paragraph
from openbiota.pdfreport import (
    AMBER,
    CONTENT_WIDTH,
    CORAL,
    GREEN,
    INK_FAINT,
    INK_SOFT,
    RULE,
    SLATE,
    Tile,
)
from openbiota.pdfsummary import GradedBar

DASH: Final = "\u2014"


def _hex(colour: colors.Color) -> str:
    return f"#{int(colour.red * 255):02X}{int(colour.green * 255):02X}{int(colour.blue * 255):02X}"


def _num(value: float | None, digits: int = 2, *, suffix: str = "") -> str:
    """A number, or an em dash. Never a zero standing in for no answer."""
    if value is None:
        return DASH
    return f"{value:,.{digits}f}{suffix}"


def _table(rows: Sequence[Sequence[Any]], widths: Sequence[float], *, header: bool = True) -> Table:
    table = Table([list(r) for r in rows], colWidths=list(widths), hAlign="LEFT", repeatRows=1 if header else 0)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.2),
        ("TOPPADDING", (0, 0), (-1, -1), 2.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.25, RULE),
    ]
    if header:
        style.append(("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE))
    table.setStyle(TableStyle(style))
    return table


# --------------------------------------------------------------------------- #
# A09 — the ecology companions, inside the gut-community section
# --------------------------------------------------------------------------- #


def ecology_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """Diversity in species-equivalents, dominance, ratios and aerotolerance.

    These sit beside the existing Shannon and evenness values, which they do
    not replace: `exp(H)` is a number of species a reader can picture, where
    an entropy in nats is not, and both are shown with the catalogue and the
    detection floor they were counted over.
    """
    if not view:
        return
    diversity = view.get("diversity") or {}
    vector = view.get("vector") or {}
    if diversity.get("effective_shannon_species") is None:
        return

    story.append(CondPageBreak(58 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "<b>Diversity, in species you can picture</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>companion views; the index above is unchanged</font>",
        st["h3"],
    ))

    width = (CONTENT_WIDTH - 3 * 3 * mm) / 4
    tiles = [
        Tile(value=_num(diversity.get("effective_shannon_species"), 1),
             label="Effective species",
             note="how many equally-common species would give this diversity",
             width=width, height=16 * mm, accent=GREEN),
        Tile(value=_num(diversity.get("inverse_simpson"), 1),
             label="Inverse Simpson",
             note="weighted toward the common organisms",
             width=width, height=16 * mm, accent=GREEN),
        Tile(value=_num(100 * (diversity.get("dominance_top1") or 0), 1, suffix="%"),
             label="Largest single organism",
             note="share of the named community",
             width=width, height=16 * mm, accent=SLATE),
        Tile(value=_num(100 * (diversity.get("dominance_top5") or 0), 1, suffix="%"),
             label="Top five together",
             note=f"of {diversity.get('observed_taxa', 0)} organisms counted",
             width=width, height=16 * mm, accent=SLATE),
    ]
    strip = Table([tiles], colWidths=[width + 3 * mm] * 3 + [width], hAlign="LEFT")
    strip.setStyle(TableStyle([
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    story.append(strip)
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_SOFT)}'>Counted over {vector.get('n_taxa_in_vector', 0)} "
        f"named organisms from the {vector.get('lane', 'inventory')} lane"
        + (f", above {vector.get('abundance_floor_percent')}% abundance" if vector.get("abundance_floor_percent") else "")
        + ". Natural logarithms. These describe the same community as the index above and do not "
        "replace it; comparing them between samples needs the same catalogue and the same floor."
        "</font>", st["small"],
    ))

    # -- ratios ----------------------------------------------------------- #
    ratios = list(view.get("ratios") or [])
    if ratios:
        story.append(Spacer(1, 2.5 * mm))
        story.append(Paragraph(
            "<b>Composition ratios</b> &nbsp;"
            f"<font size='7' color='{_hex(INK_FAINT)}'>descriptive; no healthy cut-point is defined "
            "for any of these</font>", st["h3"],
        ))
        rows: list[list[Any]] = [[
            Paragraph(h, st["label"]) for h in ("RATIO", "VALUE", "WHAT WENT INTO IT", "READING")
        ]]
        for ratio in ratios:
            if ratio["state"] == "censored_zero_denominator":
                value = Paragraph(f"<font color='{_hex(INK_FAINT)}'>not defined</font>", st["cell"])
                reading = "The denominator is absent here, so the ratio has no value."
            elif ratio["state"] == "measured_numerator_absent":
                # Printing "0" here would read as a measured zero. The
                # numerator was looked for and not found above the floor,
                # which bounds the ratio rather than fixing it.
                value = Paragraph(
                    f"<font color='{_hex(INK_FAINT)}'>below the floor</font>", st["cell"])
                reading = ("The numerator was searched for and not detected, so the ratio is "
                           "below what this depth can resolve rather than measured at zero.")
            else:
                value = Paragraph(f"<b>{ratio['value']:,.3f}</b>", st["cell"])
                reading = ratio["note"]
            members = ", ".join(ratio["numerator_members_found"]) or "none found"
            against = ", ".join(ratio["denominator_members_found"]) or "none found"
            rows.append([
                Paragraph(f"<b>{ratio['label']}</b>", st["cell"]),
                value,
                Paragraph(f"<font size='6.4'>{members} &nbsp;against&nbsp; {against} "
                          f"({ratio['level']} level)</font>", st["cell"]),
                Paragraph(f"<font size='6.4'>{reading}</font>", st["cell"]),
            ])
        story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.26, 0.10, 0.30, 0.34)]))

    # -- aerotolerance ---------------------------------------------------- #
    aero = view.get("aerotolerance")
    if aero and aero.get("aerotolerant_fraction") is not None:
        fraction = float(aero["aerotolerant_fraction"])
        coverage = aero.get("trait_coverage")
        story.append(Spacer(1, 2.5 * mm))
        story.append(Paragraph(
            "<b>Aerotolerance balance</b> &nbsp;"
            f"<font size='7' color='{_hex(INK_FAINT)}'>an ecological proxy, not a measurement of "
            "oxygen in your gut</font>", st["h3"],
        ))
        strict = 1.0 - fraction
        story.append(Paragraph(
            f"<font size='7'>Of the organisms whose oxygen requirement is established, "
            f"<b>{100 * strict:.1f}%</b> of the abundance is strictly anaerobic and "
            f"<b>{100 * fraction:.1f}%</b> tolerates oxygen. A gut community is normally "
            f"dominated by strict anaerobes."
            + (f" Oxygen requirements are established for <b>{100 * coverage:.0f}%</b> of the "
               "abundance here; the rest is counted as unknown rather than assigned to either side."
               if coverage is not None else "")
            + "</font>", st["body"],
        ))
        top = [c for c in (aero.get("contributors") or []) if c["phenotype"] == "aerotolerant"][:5]
        if top:
            names = ", ".join(
                f"<i>{c['taxon']}</i> ({100 * c['abundance_fraction']:.2f}%)" for c in top
            )
            story.append(Paragraph(
                f"<font size='6.6' color='{_hex(INK_SOFT)}'>Oxygen-tolerant organisms here: "
                f"{names}.</font>", st["small"],
            ))

    # The redundancy table moved to the end of the functions detail, where the
    # functions it describes are. See `redundancy_block`, which the
    # functions section calls.


#: Where a function's effective carrier count stops being a comfort.
#:
#: One effective carrier means a single organism holds the function whatever
#: the raw count says; the scale is generous above three because the
#: difference between four carriers and five is not a difference a reader
#: should act on.
_REDUNDANCY_BANDS: Final[tuple[tuple[float, colors.Color, colors.Color], ...]] = (
    (25.0, colors.HexColor("#F2CFCB"), CORAL),
    (50.0, colors.HexColor("#F0C4A8"), colors.HexColor("#D0641E")),
    (75.0, colors.HexColor("#F6DDB4"), AMBER),
    (100.0, colors.HexColor("#CFE4CF"), GREEN),
)

#: Effective carriers are unbounded above; the bar is not. Four is the point
#: past which more carriers stop changing the answer, so it anchors the top.
_REDUNDANCY_FULL: Final = 4.0


def _redundancy_position(effective: float) -> float:
    return max(0.0, min(100.0, 100.0 * effective / _REDUNDANCY_FULL))


def _redundancy_tone(effective: float) -> colors.Color:
    position = _redundancy_position(effective)
    for upper, _track, marker in _REDUNDANCY_BANDS:
        if position < upper or upper >= 100.0:
            return marker
    return GREEN


def redundancy_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """§6.4 — how many organisms actually hold each function.

    A different question from how much of the function there is, and the
    one that decides whether a reading is fragile: a route carried by a
    single organism moves with that organism, however healthy the level
    looks today. Two carriers at 50/50 count as two; two at 99/1 count as
    barely more than one, which is why the effective count is used and a
    raw tally is not.
    """
    readings = [
        r for r in ((view or {}).get("redundancy") or [])
        if isinstance(r, dict) and r.get("effective_carriers") is not None
    ]
    if not readings:
        return
    readings.sort(key=lambda r: float(r["effective_carriers"]))

    story.append(CondPageBreak(56 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "<b>How many organisms hold each function</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>fewest confirmed carriers first</font>",
        st["h3"],
    ))
    story.append(Paragraph(
        "Two communities can make the same amount of something with very different security. "
        "If one organism is doing all of it, the function travels with that organism; if four "
        "share it, losing any one changes little. The count below is an <i>effective</i> one: "
        "two organisms splitting a job evenly count as two, but two where one does 99% of the "
        "work count as barely more than one.", st["body"],
    ))

    rows: list[list[Any]] = [[
        Paragraph(h, st["label"])
        for h in ("FUNCTION", "EFFECTIVE CARRIERS", "", "NAMED", "TOP SHARE", "SIGNAL CONFIRMED")
    ]]
    weakly_supported = 0
    for reading in readings:
        effective = float(reading["effective_carriers"])
        tone = _redundancy_tone(effective)
        carriers = list(reading.get("carriers") or [])
        total = sum(float(c.get("value") or 0.0) for c in carriers) or None
        top = max((float(c.get("value") or 0.0) for c in carriers), default=0.0)
        share = f"{100.0 * top / total:.0f}%" if total else DASH

        # What fraction of the function's signal belongs to organisms we can
        # actually place in this sample. When it is low the carrier count is
        # a floor rather than a finding: the rest of the signal sits on
        # reference organisms that may or may not be here.
        confirmed = reading.get("carrier_coverage")
        if confirmed is None:
            confirmed_cell = Paragraph(f"<font color='{_hex(INK_FAINT)}'>{DASH}</font>", st["cell"])
        else:
            pct = 100.0 * float(confirmed)
            weak = pct < 50.0
            weakly_supported += 1 if weak else 0
            colour = AMBER if weak else INK_SOFT
            confirmed_cell = Paragraph(
                f"<font color='{_hex(colour)}'>{pct:.0f}%</font>"
                + ("<font size='5.6'> \u2020</font>" if weak else ""), st["cell"],
            )
        rows.append([
            Paragraph(str(reading.get("label") or reading.get("function_id")), st["cell"]),
            Paragraph(f"<b><font color='{_hex(tone)}'>{effective:.1f}</font></b>", st["cell"]),
            GradedBar(width=CONTENT_WIDTH * 0.22, value=_redundancy_position(effective),
                      bands=_REDUNDANCY_BANDS, marker_colour=tone),
            Paragraph(str(reading.get("n_resolved_carriers") or 0), st["cell"]),
            Paragraph(share, st["cell"]),
            confirmed_cell,
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.28, 0.11, 0.23, 0.08, 0.11, 0.19)]))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'>The bar runs to "
        f"{_REDUNDANCY_FULL:.0f} effective carriers, past which more organisms stop changing "
        "the answer. \u201cNamed\u201d counts only organisms found in your own inventory as well "
        "as on the gene: a reference organism whose copy of a gene matched a read, but which "
        "was not detected in you, is context rather than a carrier. <b>Signal confirmed</b> is "
        "how much of the function's signal those named carriers account for."
        + (f" <font color='{_hex(AMBER)}'>\u2020 On {weakly_supported} row"
           f"{'' if weakly_supported == 1 else 's'} it is under half</font>, meaning most of the "
           "signal sits on organisms this assay could not place in you \u2014 there read the "
           "carrier count as the fewest it could be, not as the number there are."
           if weakly_supported else "")
        + " This is redundancy among the carriers this assay resolved \u2014 not a "
        "guarantee of ecological resilience.</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# A10 — the organism explorer, with the full catalogue
# --------------------------------------------------------------------------- #


def explorer_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """A searchable answer for any organism a reader has heard of.

    The point of this table is the negative answer. Most people look up one
    organism they have read about, and "we searched for it and it is not
    here" is a result they can act on — provided it is not confused with
    "we never looked".
    """
    if not view:
        return
    sentinels = view.get("sentinels") or []
    if not sentinels:
        return

    story.append(CondPageBreak(75 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Organisms people ask about", st["h2"]))
    story.append(Paragraph(
        f"The organisms most often named in gut-health writing, each looked up in your "
        f"results. <b>{view.get('n_sentinels_detected', 0)} of {len(sentinels)}</b> were found "
        "in your sample. A blank here means this test searched and did not find it, which is "
        "not the same as the organism being absent from you — a low-abundance organism can "
        "sit below what one stool sample resolves.", st["body"],
    ))

    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in ("ORGANISM", "IN YOUR SAMPLE", "LEVEL", "NOTE")
    ]]
    for s_ in sentinels:
        state = s_["detection_state"]
        if state == "supported_detection":
            found = f"<font color='{_hex(GREEN)}'><b>found</b></font>"
            level = (f"<b>{s_['percent']:.3f}%</b>" if s_.get("percent") is not None else DASH)
            if s_.get("percentile") is not None:
                level += f" <font size='6' color='{_hex(INK_FAINT)}'>{s_['percentile']:.0f}th pct</font>"
        elif state == "below_quantification_limit":
            # Detected, but under what this sample quantifies. Printing this
            # as "not found" would contradict the catalogue page, which lists
            # the same organism with an abundance.
            found = f"<font color='{_hex(AMBER)}'>trace</font>"
            level = (f"{s_['percent']:.3f}%" if s_.get("percent") else "below the reporting level")
        else:
            found = f"<font color='{_hex(INK_FAINT)}'>not found</font>"
            level = DASH
        note = ""
        if s_["resolution_state"] == "accepted_synonym" and s_.get("accepted_name"):
            note = f"listed here as <i>{s_['accepted_name']}</i>"
        elif s_["resolution_state"] == "rank_rollup":
            note = f"genus total over {len(s_.get('candidates') or [])} species"
        elif s_["resolution_state"] == "not_in_catalogue":
            note = "searched, not detected"
        rows.append([
            Paragraph(f"<i>{s_['query']}</i>", st["cell"]),
            Paragraph(found, st["cell"]),
            Paragraph(level, st["cell"]),
            Paragraph(f"<font size='6.2'>{note}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.34, 0.14, 0.20, 0.32)]))

    role = view.get("role_composition") or {}
    classes = role.get("classes") or []
    if classes:
        story.append(Spacer(1, 2.5 * mm))
        story.append(Paragraph(
            "<b>What kind of organisms make up your community</b>", st["h3"],
        ))
        parts = " · ".join(
            f"<b>{c['percent']:.1f}%</b> {c['role']}" for c in classes
        )
        unassigned = role.get("unassigned_percent") or 0.0
        story.append(Paragraph(
            f"<font size='7'>{parts}"
            + (f" · <font color='{_hex(INK_FAINT)}'>{unassigned:.1f}% not yet classified</font>"
               if unassigned > 0.05 else "")
            + f". Counted over {view.get('n_organisms', 0)} organisms in "
            f"{view.get('n_genera', 0)} genera. Organisms with no assigned class are shown as "
            "unclassified rather than dropped, and “unknown” describes the evidence, "
            "not the organism.</font>", st["body"],
        ))



def resistance_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """Class-level resistance determinants, with mechanism and carrier state."""
    if not view:
        return
    classes = [c for c in (view.get("classes") or []) if c["n_detected"]]
    totals = view.get("totals") or {}

    story.append(CondPageBreak(62 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Antibiotic-resistance genes, by class", st["h2"]))
    story.append(Paragraph(
        "Resistance genes present in the community, grouped by the antibiotics they act "
        "against. <b>A gene is not a phenotype.</b> Carrying one of these does not establish "
        "that an antibiotic would fail for you, and not carrying one does not establish that "
        "it would work. Gut bacteria carry resistance genes routinely, most of them in "
        "harmless residents.", st["body"],
    ))
    if view.get("database_note"):
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{view['database_note']}</font>",
            st["small"],
        ))

    if not classes:
        story.append(Paragraph(
            f"<b>No resistance determinant was detected</b> above the screen's threshold, out of "
            f"{totals.get('n_classes_reported', 0)} classes searched.", st["body"],
        ))
        return

    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in
        ("ANTIBIOTIC CLASS", "GENE FAMILIES", "EVIDENCE", "HOW IT WORKS", "READS")
    ]]
    for cls in classes:
        evidence_bits = []
        if cls["n_acquired_genes"]:
            evidence_bits.append(
                f"<font color='{_hex(AMBER)}'>{cls['n_acquired_genes']} acquired gene"
                f"{'s' if cls['n_acquired_genes'] != 1 else ''}</font>")
        if cls["n_validated_mutations"]:
            evidence_bits.append(f"{cls['n_validated_mutations']} resistance mutation(s)")
        if cls["n_unresolved_homologs"]:
            evidence_bits.append(
                f"<font color='{_hex(INK_FAINT)}'>{cls['n_unresolved_homologs']} unresolved</font>")
        families = ", ".join(
            d["gene_family"] for d in cls["determinants"][:4] if d["gene_family"]
        ) or DASH
        rows.append([
            Paragraph(f"<b>{cls['label']}</b>", st["cell"]),
            Paragraph(f"<b>{cls['determinant_richness']}</b> "
                      f"<font size='6' color='{_hex(INK_FAINT)}'>{families}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{' · '.join(evidence_bits)}</font>", st["cell"]),
            Paragraph(f"<font size='6.2'>{cls['mechanism']}</font>", st["cell"]),
            Paragraph(f"{cls['supporting_fragments']:,}", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.24, 0.20, 0.18, 0.28, 0.10)]))

    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        f"<font size='6.8'><b>{totals.get('total_determinant_richness', 0)} distinct gene families</b> "
        f"across {totals.get('n_classes_with_a_detection', 0)} of "
        f"{totals.get('n_classes_reported', 0)} classes. Families are counted once even when they "
        f"appear in two classes, and a drug is never counted again as its own class. "
        f"{view.get('carriage_note', '')}</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# A16 — the simulation, inside the actions section
# --------------------------------------------------------------------------- #


def context_block(
    story: list[Any], st: Mapping[str, Any], view: Mapping[str, Any] | None
) -> None:
    """A13's symptom and condition navigation cards."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.context_block(story, st, view)


def drug_metabolism_block(
    story: list[Any], st: Mapping[str, Any], views: Mapping[str, Any] | None
) -> None:
    """Characterised microbial drug reactions, as mechanism cards."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.drug_metabolism_block(story, st, views)


def all_functional_rows(
    views: Any,
    directions: Mapping[str, str] | None = None,
    aggregates: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Every A01-A08 and §5.8 reading, as rows, in §12.1 group order."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    return pdfsynbiotic.all_functional_rows(views, directions, aggregates)


def reading_dest(group: str, label: str, *, detail: bool = False) -> str:
    """The destination name for one functional reading."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    return pdfsynbiotic.reading_dest(group, label, detail=detail)


def grouped_panels(rows: Any) -> list[tuple[str, list[Any]]]:
    """Panel rows in the §12.1 reader-facing group order."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    return pdfsynbiotic.grouped_panels(rows)


def glance_group_order() -> tuple[str, ...]:
    """The §12.1 group order, for callers laying out their own headings."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    return pdfsynbiotic.glance_group_order()


def functional_glance_block(
    story: list[Any], st: Mapping[str, Any], views: Mapping[str, Any] | None,
    *, detail_section: int, aggregates: Mapping[str, Sequence[str]] | None = None,
) -> None:
    """A01-A08 and §5.8 readings, in the `report.functions` overview."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.functional_glance_block(story, st, views, detail_section=detail_section, aggregates=aggregates)


def functional_lane_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None,
    aggregates: Mapping[str, Any] | None = None,
) -> None:
    """A01-A08 readings; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.functional_lane_block(story, st, views, aggregates)


def planner_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """A12 consolidated plan; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.planner_block(story, st, view)


def longitudinal_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """A13 change over time; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.longitudinal_block(story, st, view)


def metric_index_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], extension: Mapping[str, Any] | None
) -> None:
    """A15 complete metric index; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.metric_index_block(story, st, extension)


def composition_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """A09 composition and community type; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.composition_block(story, st, view)


def simulation_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """The A16 scenarios; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.simulation_block(story, st, view)


def simulation_detail_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """The A16 mechanism detail; rendered by `openbiota.pdfsynbiotic`."""
    from openbiota import pdfsynbiotic  # noqa: PLC0415

    pdfsynbiotic.simulation_detail_block(story, st, view)


# --------------------------------------------------------------------------- #
# what this release still owes
# --------------------------------------------------------------------------- #


def owed_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], unavailable: Sequence[Mapping[str, Any]]
) -> None:
    """Name what is not here, and say which kind of missing it is.

    An unbuilt feature and a sample that cannot answer are different facts,
    and only one of them is a limit of the science. Printing them under one
    heading would let the first hide behind the second.
    """
    if not unavailable:
        return
    debt = [u for u in unavailable if u.get("is_engineering_debt")]
    data = [u for u in unavailable if not u.get("is_engineering_debt")]
    story.append(CondPageBreak(40 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Measurements not in this edition", st["h3"]))
    if data:
        story.append(Paragraph(
            "<b>Your sample could not answer these.</b> " + "; ".join(
                f"{u['label']} ({u['detail']})" for u in data[:6]
            ) + ".", st["small"],
        ))
    if debt:
        story.append(Paragraph(
            f"<font size='6.8' color='{_hex(INK_SOFT)}'><b>Built but not yet released:</b> "
            + ", ".join(sorted(u["label"] for u in debt))
            + ". These are engineering work in progress, not limits of what the sample can "
            "show.</font>", st["small"],
        ))


__all__ = [
    "ecology_block",
    "explorer_block",
    "redundancy_block",
    "owed_block",
    "resistance_block",
    "simulation_block",
]
