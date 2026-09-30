"""The strain-resolution section: which population, not just which species.

Species profiling says an organism is present. This section says *which one* —
the dominant population's fingerprint across marker genes, for every organism
with enough coverage to characterise. That distinction is most of the
difference between "you have E. coli" and anything useful.

Design rules, because a strain section is easy to make either useless or
overclaiming:

* Lead with what it is for, in one sentence, before any number.
* Never print a strain *name*. A marker fingerprint is this sample's own
  population, not a match to a named isolate (spec §3.1, acceptance 6/12).
* Say on the page that a dominant consensus can hide a minority population
  (acceptance 10), rather than burying it in a technical appendix.
* Organisms detected but not resolved are listed with the reason. Absence of
  resolution is not absence of the organism.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import Spacer, Table, TableStyle

from openbiota.pdflinks import Paragraph, section_heading

CONTENT_WIDTH = 170 * mm
DASH = "\u2014"


def _hex(c: colors.Color) -> str:
    return f"#{int(c.red * 255):02X}{int(c.green * 255):02X}{int(c.blue * 255):02X}"


def strain_pages(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    section: int,
    strain: Mapping[str, Any] | None,
) -> None:
    """Render the strain-resolution section, or say why it is absent."""
    from openbiota.pdfreport import (
        ACCENT,
        CARD_BG,
        INK_FAINT,
        INK_SOFT,
        PANEL_BG,
        RULE,
    )

    # `section_heading` emits the destination and the outline entry itself.
    story.extend(
        section_heading(section, "Which strain, not just which species", st["h1"])
    )

    if not strain or strain.get("status") != "resolved":
        story.append(
            Paragraph(
                str(
                    (strain or {}).get("what_this_is")
                    or "Strain resolution was not run for this sample."
                ),
                st["body"],
            )
        )
        if strain and strain.get("how_to_run"):
            story.append(
                Paragraph(
                    f"<font color='{_hex(INK_FAINT)}'>{strain['how_to_run']}</font>",
                    st["small"],
                )
            )
        return

    resolved = int(strain.get("organisms_resolved") or 0)
    unresolved = int(strain.get("organisms_detected_unresolved") or 0)
    bases = int(strain.get("total_callable_bases") or 0)
    rows = list(strain.get("organisms") or [])

    story.append(
        Paragraph(
            "Two people can carry the same species and a different <b>strain</b> of it, and "
            "the difference often matters more than the species does: whether an "
            "<i>E. coli</i> is a harmless resident or carries a toxin operon is a "
            "strain-level question. This section reports, for each organism with enough "
            "coverage, the fingerprint of the population actually living in this sample \u2014 "
            "read from the marker genes, at single-base resolution.",
            st["body"],
        )
    )
    story.append(Spacer(1, 2 * mm))

    # Headline figures.
    tiles = [[
        _tile("ORGANISMS RESOLVED", f"{resolved}", "population fingerprinted", st, ACCENT),
        _tile("CALLABLE BASES", f"{bases / 1e6:.1f}M", "compared at base level", st, ACCENT),
        _tile("DETECTED, UNRESOLVED", f"{unresolved}", "too little coverage", st, INK_SOFT),
    ]]
    table = Table(tiles, colWidths=[CONTENT_WIDTH / 3] * 3, hAlign="LEFT")
    table.setStyle(
        TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
            ("BOX", (0, 0), (-1, -1), 0.5, RULE),
            ("LINEAFTER", (0, 0), (-2, -1), 0.5, RULE),
            ("LEFTPADDING", (0, 0), (-1, -1), 5 * mm),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
            ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3.4 * mm),
            ("ROUNDEDCORNERS", [4, 4, 4, 4]),
        ])
    )
    story.append(table)
    story.append(Spacer(1, 3 * mm))

    # What this can and cannot say -- on the page, not in an appendix.
    story.append(Paragraph("What a fingerprint can and cannot tell you", st["h2"]))
    for term, text in (
        (
            "It is this sample's own population",
            "Not a match to a named laboratory strain. Naming one would require comparing "
            "the whole genome against that isolate, which marker genes cannot do.",
        ),
        (
            "It is the dominant population",
            "Where two populations of one species are present, this reports the majority "
            "one. A less abundant second population can be hidden by it entirely.",
        ),
        (
            "It enables comparison",
            "Because the fingerprint is read at single-base resolution, two samples can be "
            "compared organism by organism \u2014 which is how donor and recipient "
            "populations are told apart in the matching report.",
        ),
        (
            "It does not attribute genes",
            "Resolving an organism's population does not show which organism carries a "
            "separately detected toxin or resistance gene. That needs physical linkage.",
        ),
    ):
        story.append(Paragraph(f"<b>{term}.</b> {text}", st["small"]))
        story.append(Spacer(1, 1.2 * mm))
    story.append(Spacer(1, 2 * mm))

    # The organisms.
    story.append(Paragraph("Every organism resolved to population level", st["h2"]))
    story.append(
        Paragraph(
            f"<font color='{_hex(INK_FAINT)}'>Ordered by how much sequence was callable, "
            "which is how much of the population could actually be characterised. "
            "<b>Markers</b> is how many marker genes resolved; <b>breadth</b> is how much of "
            "each was covered; <b>depth</b> is the average read depth over them.</font>",
            st["small"],
        )
    )

    header = ["ORGANISM", "SGB", "MARKERS", "CALLABLE BASES", "BREADTH", "DEPTH"]
    body: list[list[Any]] = [[
        Paragraph(f"<font size='6.2' color='{_hex(INK_FAINT)}'><b>{h}</b></font>", st["cell"])
        for h in header
    ]]
    for row in rows:
        species = str(row.get("species") or "").replace("_", " ") or str(row.get("sgb"))
        unnamed = species.startswith("GGB") or not row.get("species")
        label = (
            f"<font color='{_hex(INK_SOFT)}'>{species}</font>"
            if unnamed else f"<i>{species}</i>"
        )
        body.append([
            Paragraph(label, st["cell"]),
            Paragraph(f"<font size='6.4' color='{_hex(INK_FAINT)}'>{row.get('sgb')}</font>",
                      st["cell"]),
            Paragraph(f"{row.get('markers_resolved')}", st["cell"]),
            Paragraph(f"{int(row.get('callable_bases') or 0):,}", st["cell"]),
            Paragraph(f"{row.get('median_breadth_percent')}%", st["cell"]),
            Paragraph(f"{row.get('median_depth')}\u00d7", st["cell"]),
        ])
    widths = [
        CONTENT_WIDTH * 0.36, CONTENT_WIDTH * 0.13, CONTENT_WIDTH * 0.11,
        CONTENT_WIDTH * 0.19, CONTENT_WIDTH * 0.11, CONTENT_WIDTH * 0.10,
    ]
    organisms = Table(body, colWidths=widths, repeatRows=1, hAlign="LEFT")
    organisms.setStyle(
        TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (2, 0), (-1, -1), "RIGHT"),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT),
            ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, CARD_BG]),
            ("TOPPADDING", (0, 0), (-1, -1), 1.7 * mm),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.7 * mm),
            ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
        ])
    )
    story.append(organisms)
    story.append(Spacer(1, 2 * mm))

    if unresolved:
        story.append(
            Paragraph(
                f"<b>{unresolved} further organism(s)</b> were detected but carried too few "
                "resolvable marker genes to characterise a population. They are present; "
                "their strain is simply unresolved at this sequencing depth, which is a "
                "limit of the data rather than a finding about the organism.",
                st["small"],
            )
        )
    for limit in strain.get("limits") or []:
        story.append(
            Paragraph(f"<font color='{_hex(INK_FAINT)}'>{limit}</font>", st["small"])
        )
    story.append(
        Paragraph(
            f"<font color='{_hex(INK_FAINT)}'>Method: "
            f"{strain.get('method')} \u00b7 database {strain.get('database_release')}."
            "</font>",
            st["small"],
        )
    )


def _acc(name: Any) -> str:
    return str(name or "").split("|")[0].replace(".fna.gz", "").replace(".fna", "")


def placement_pages(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    analysis: Mapping[str, Any] | None,
    inventory: Mapping[str, Any] | None = None,
) -> None:
    """Comparative placement (StrainPhlAn 4): where each population sits among
    the reference genomes of its own species.

    Reported separately from the fingerprint above because it answers a
    different question - *which known genome is this population closest to* -
    and carries its own caveat: a nearest reference is a placement, not proof
    of an identical resident strain, and a consensus can hide a second
    population. Unresolved organisms are listed with the reason; a deferred
    tree is said to be deferred.
    """
    from openbiota.pdfreport import ACCENT, CARD_BG, INK_FAINT, INK_SOFT, PANEL_BG, RULE

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("Placed among the reference genomes of its species", st["h2"]))
    if not analysis or analysis.get("status") != "completed":
        reason = str((analysis or {}).get("reason") or "the comparative step did not run for this sample")
        story.append(Paragraph(
            f"Comparative placement was not performed: {reason}. The fingerprints above stand on their own; "
            "nothing below is inferred.", st["small"]))
        return

    names: dict[str, str] = {}
    for rec in (inventory or {}).get("organisms") or []:
        if rec.get("sgb"):
            names[str(rec["sgb"])] = str(rec.get("species") or rec.get("gtdb") or rec["sgb"])
    clades = list(analysis.get("clades") or [])
    placed = [c for c in clades if c.get("status") == "placed"]
    built = [c for c in clades if c.get("status") == "tree_built"]
    failed = [c for c in clades if c.get("status") == "unresolved"]
    too_few = list(analysis.get("unresolved") or [])
    deferred = list(analysis.get("deferred") or [])
    n_eligible = int(analysis.get("n_eligible") or 0)

    story.append(Paragraph(
        "For every organism with enough reconstructed marker sequence, its consensus markers were aligned "
        "with the marker sequences of member genomes of the same species (GTDB R232) and a tree was inferred. "
        "The <b>nearest reference</b> is the genome the sample's population sits closest to in that tree and the "
        "<b>distance</b> is the tree's branch length between them \u2014 comparable within one tree, "
        "not a genome-wide identity. "
        "A placement is not proof of an identical resident strain: it says which known genome is closest, "
        "not that they are the same. <b>Polymorphic sites</b> in the consensus are the mixture evidence — "
        "a high rate suggests more than one population of the species is present, and the consensus then "
        "describes neither exactly. No pathogenicity or resistance is inferred from a neighbour.",
        st["small"]))
    story.append(Spacer(1, 2 * mm))
    tiles = [[
        _tile("PLACED", f"{len(placed)}", f"of {n_eligible} eligible organisms", st, ACCENT),
        _tile("UNRESOLVED", f"{len(failed) + len(too_few)}", "too few markers or references", st, INK_SOFT),
        _tile("DEFERRED", f"{len(deferred)}", "left for the next run", st, INK_SOFT),
    ]]
    table = Table(tiles, colWidths=[CONTENT_WIDTH / 3] * 3, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"), ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
        ("BOX", (0, 0), (-1, -1), 0.5, RULE), ("LINEAFTER", (0, 0), (-2, -1), 0.5, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 5 * mm), ("RIGHTPADDING", (0, 0), (-1, -1), 3 * mm),
        ("TOPPADDING", (0, 0), (-1, -1), 3 * mm), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.4 * mm),
        ("ROUNDEDCORNERS", [4, 4, 4, 4]),
    ]))
    story.append(table)
    story.append(Spacer(1, 3 * mm))

    if placed or built:
        header = ["ORGANISM", "NEAREST REFERENCE", "DISTANCE", "MARKERS", "POLYMORPHIC", "GENOMES"]
        body: list[list[Any]] = [[
            Paragraph(f"<font size='6.2' color='{_hex(INK_FAINT)}'><b>{h}</b></font>", st["cell"]) for h in header]]
        for c in placed + built:
            sgb = str(c.get("sgb") or "")
            species = names.get(sgb, sgb).replace("_", " ")
            unnamed = species.startswith("GGB") or species == sgb
            label = (f"<font color='{_hex(INK_SOFT)}'>{species}</font>" if unnamed else f"<i>{species}</i>") + \
                f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>{sgb}</font>"
            dist = c.get("distance_to_nearest")
            rate = c.get("polymorphic_rate")
            body.append([
                Paragraph(label, st["cell"]),
                Paragraph(f"<font size='6.4'>{_acc(c.get('nearest_reference')) or DASH}</font>", st["cell"]),
                Paragraph(f"{float(dist):.3f}" if isinstance(dist, (int, float)) else DASH, st["cell"]),
                Paragraph(f"{c.get('n_markers') or DASH}", st["cell"]),
                Paragraph(f"{float(rate):.2%}" if isinstance(rate, (int, float)) else "\u2014", st["cell"]),
                Paragraph(f"{len(c.get('references') or [])}", st["cell"]),
            ])
        widths = [CONTENT_WIDTH * 0.34, CONTENT_WIDTH * 0.22, CONTENT_WIDTH * 0.11,
                  CONTENT_WIDTH * 0.10, CONTENT_WIDTH * 0.13, CONTENT_WIDTH * 0.10]
        t = Table(body, colWidths=widths, repeatRows=1, hAlign="LEFT")
        t.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (2, 0), (-1, -1), "RIGHT"),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT), ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, CARD_BG]),
            ("TOPPADDING", (0, 0), (-1, -1), 1.5 * mm), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5 * mm),
            ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm), ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
        ]))
        story.append(t)
        story.append(Spacer(1, 2 * mm))
    else:
        story.append(Paragraph("No organism could be placed in this sample.", st["small"]))

    if failed:
        story.append(Paragraph("<b>Eligible but unresolved</b>, with the reason:", st["small"]))
        for c in failed:
            sgb = str(c.get("sgb") or "")
            species = names.get(sgb, sgb).replace("_", " ")
            reason = str(c.get("reason") or "").split(":")[0]
            story.append(Paragraph(
                f"<font color='{_hex(INK_FAINT)}'>\u2022 <i>{species}</i> ({sgb}): {reason}</font>", st["small"]))
        story.append(Spacer(1, 1.5 * mm))
    if too_few:
        story.append(Paragraph(
            f"<font color='{_hex(INK_FAINT)}'>{len(too_few)} further marker clades had fewer than "
            f"{(analysis.get('operating_points') or {}).get('min_markers', 20)} reconstructed markers and were not "
            "placed; they are present as organisms where the inventory lists them, and unresolved only at strain "
            "level.</font>", st["small"]))
    if deferred:
        story.append(Paragraph(
            f"<font color='{_hex(INK_FAINT)}'>{len(deferred)} eligible organisms were deferred by this run's "
            "time budget and will be placed on the next run; they are not reported as unresolved.</font>",
            st["small"]))
    story.append(Paragraph(
        f"<font color='{_hex(INK_FAINT)}'>Method: StrainPhlAn 4 on the {analysis.get('database')} markers, "
        "reference genomes from the organism's GTDB R232 species cluster (representative first, at most "
        f"{(analysis.get('operating_points') or {}).get('max_references', 6)}), PhyloPhlAn fast mode. "
        "Strain results never add to species counts.</font>", st["small"]))


def census_pages(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    *,
    census: Mapping[str, Any] | None,
) -> None:
    """What was asked for, and what each request actually reached.

    This exists so an incomplete analysis cannot present itself as a complete
    screen. Every registered target appears with its terminal state, including
    the ones that could not run and why — a target that was never assayed is
    shown as never assayed, not as a clean result.
    """
    from openbiota.pdfreport import ACCENT, CARD_BG, INK_FAINT, RULE

    if not census:
        return

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("What was asked for, and what came back", st["h2"]))
    story.append(
        Paragraph(
            "A screen is only meaningful alongside its own coverage. Every resolution "
            "target registered for this sample is listed below with the state it reached. "
            "<b>An unrun, failed or unavailable assay is never reported as a negative "
            "result</b> \u2014 the distinction between \u201clooked and found nothing\u201d "
            "and \u201cnever looked\u201d is preserved throughout.",
            st["small"],
        )
    )

    labels = {
        "completed": "Completed",
        "not_requested": "Not requested",
        "scheduled": "Registered, awaiting reference curation",
        "failed": "Failed",
        "reference_unavailable": "Reference or tool unavailable",
        "access_or_license_unavailable": "Access or licence unavailable",
        "incompatible_input": "Incompatible with DNA input",
    }
    meanings = {
        "completed": "the assay ran and produced a result",
        "not_requested": "not part of this run's target set",
        "scheduled": "registered, with reference curation outstanding",
        "failed": "attempted and did not complete",
        "reference_unavailable": "the tool or reference is not installed",
        "access_or_license_unavailable": "blocked by access or licence terms",
        "incompatible_input": "an RNA genome cannot be assayed from DNA sequencing",
    }

    rows: list[list[Any]] = [[
        Paragraph(f"<font size='6.2' color='{_hex(INK_FAINT)}'><b>{h}</b></font>", st["cell"])
        for h in ("STATE", "TARGETS", "WHAT IT MEANS")
    ]]
    for state, count in (census.get("by_assay_status") or {}).items():
        rows.append([
            Paragraph(f"<b>{labels.get(state, state)}</b>", st["cell"]),
            Paragraph(f"{count}", st["cell"]),
            Paragraph(
                f"<font color='{_hex(INK_FAINT)}'>{meanings.get(state, '')}</font>",
                st["cell"],
            ),
        ])
    table = Table(
        rows,
        colWidths=[CONTENT_WIDTH * 0.42, CONTENT_WIDTH * 0.12, CONTENT_WIDTH * 0.46],
        hAlign="LEFT",
    )
    table.setStyle(
        TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (1, 0), (1, -1), "RIGHT"),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, ACCENT),
            ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, CARD_BG]),
            ("TOPPADDING", (0, 0), (-1, -1), 1.8 * mm),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.8 * mm),
            ("LEFTPADDING", (0, 0), (-1, -1), 2 * mm),
        ])
    )
    story.append(table)
    story.append(Spacer(1, 2 * mm))

    registered = int(census.get("targets_registered") or 0)
    outstanding = int(census.get("targets_outstanding") or 0)
    story.append(
        Paragraph(
            f"<b>{registered} targets registered</b>, of which <b>{outstanding}</b> remain "
            "outstanding. "
            + (
                "Every registered target reached a terminal state."
                if census.get("complete")
                else "The outstanding targets are listed in results.json under "
                "<font face='Courier'>resolution_census</font>, each with the reason it has "
                "not yet been determined."
            ),
            st["small"],
        )
    )

    schemes = census.get("typing_schemes") or {}
    if schemes:
        story.append(
            Paragraph(
                f"<font color='{_hex(INK_FAINT)}'>Organism-specific typing: "
                f"{schemes.get('n_schemes')} schemes registered "
                f"({schemes.get('n_installed')} with their tool installed). A scheme without "
                "its tool reports an unavailable capability for every sample; it never "
                "reports a negative type, and a species name is never allowed to stand in "
                "for a type.</font>",
                st["small"],
            )
        )
    adapters = census.get("tool_adapters") or {}
    if adapters:
        story.append(
            Paragraph(
                f"<font color='{_hex(INK_FAINT)}'>Strain-comparison tools: "
                f"{adapters.get('n_installed')} of {adapters.get('n_adapters')} installed. "
                f"{adapters.get('no_consensus_bonus', '')}</font>",
                st["small"],
            )
        )


def _tile(
    label: str,
    value: str,
    sub: str,
    st: Mapping[str, ParagraphStyle],
    colour: colors.Color,
) -> list[Any]:
    from openbiota.pdfreport import ACCENT_DARK, INK_FAINT

    return [
        Paragraph(
            f"<font size='6.4' color='{_hex(ACCENT_DARK)}'><b>{label}</b></font>",
            ParagraphStyle("tk", parent=st["small"], fontSize=6.4, leading=9,
                           spaceAfter=0, spaceBefore=0),
        ),
        Paragraph(
            f"<font size='19' color='{_hex(colour)}'><b>{value}</b></font>",
            ParagraphStyle("tv", parent=st["body"], fontSize=19, leading=23,
                           spaceAfter=0, spaceBefore=1),
        ),
        Paragraph(
            f"<font size='7' color='{_hex(INK_FAINT)}'>{sub}</font>",
            ParagraphStyle("ts", parent=st["small"], fontSize=7, leading=9.5,
                           spaceAfter=0, spaceBefore=0),
        ),
    ]


def strain_glance(
    st: Mapping[str, ParagraphStyle],
    *,
    strain: Mapping[str, Any] | None,
    section: int,
) -> list[Any]:
    """A one-line summary for the summary page, or nothing."""
    if not strain or strain.get("status") != "resolved":
        return []
    resolved = int(strain.get("organisms_resolved") or 0)
    if not resolved:
        return []
    bases = int(strain.get("total_callable_bases") or 0)
    return [
        Paragraph(
            f"<b>{resolved} organisms</b> resolved to strain level "
            f"({bases / 1e6:.1f}M callable bases) \u2014 <b>section {section}</b>",
            st["small"],
        )
    ]


def self_test() -> int:
    """Smoke-check that both branches build flowables. Returns failures."""
    from openbiota.pdfreport import _styles

    failures = 0
    st = _styles()

    absent: list[Any] = []
    strain_pages(absent, st, section=15, strain={"status": "not_run",
                                                "what_this_is": "x", "how_to_run": "y"})
    if not absent:
        failures += 1

    placed: list[Any] = []
    placement_pages(placed, st, analysis={
        "status": "completed", "database": "db", "n_eligible": 3,
        "clades": [{"sgb": "SGB1", "status": "placed", "nearest_reference": "GCA_1.fna", "distance_to_nearest": 0.01,
                    "n_markers": 50, "polymorphic_rate": 0.001, "references": ["a", "b", "c"]},
                   {"sgb": "SGB2", "status": "unresolved", "reason": "only 2 reference genomes could be fetched"}],
        "unresolved": [{"sgb": "SGB3", "n_markers": 4}], "deferred": ["SGB4"],
    }, inventory={"organisms": [{"sgb": "SGB1", "species": "Escherichia_coli"}]})
    if not placed:
        failures += 1

    present: list[Any] = []
    strain_pages(
        present, st, section=15,
        strain={
            "status": "resolved",
            "database_release": "db",
            "method": "m",
            "organisms_resolved": 1,
            "organisms_detected_unresolved": 2,
            "total_callable_bases": 12345,
            "organisms": [{
                "sgb": "SGB1", "species": "Escherichia_coli", "markers_resolved": 42,
                "callable_bases": 12345, "median_breadth_percent": 95.0,
                "median_depth": 8.0, "is_named_strain": False,
            }],
            "limits": ["a limit"],
        },
    )
    if not present:
        failures += 1
    return failures
