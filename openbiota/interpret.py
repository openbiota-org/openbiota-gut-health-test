"""Turn tallied results plus reference ranges into report-ready rows.

Kept separate from both the tally (which does the science) and the PDF (which
does the layout), so the comparison logic is testable on its own.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from openbiota.cohort import ReferenceRanges
from openbiota.confidence import (
    CONFIRMED,
    MIN_MEAN_IDENTITY,
    PROVISIONAL,
    reasons_below_confirmed,
)
from openbiota.pdfreport import MetaboliteRow, status_for
from openbiota.tally import PanelResult


def load_optional_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def cohort_note(ranges: ReferenceRanges | None) -> str:
    """The reference group, named. A percentile against an unnamed group is
    a number without a meaning; against a named one the reader can judge the
    comparison — and for a child compared against an older-adult cohort, that
    judgement matters."""
    if ranges is None:
        return "a reference group (unavailable)"
    who = f"{ranges.n_samples} adult stool samples"
    if ranges.study_title:
        who += f" ({ranges.study_title})"
    return who


def build_rows(
    panels: list[PanelResult],
    ranges: ReferenceRanges | None,
    *,
    sample_dir: Path | None = None,
    inventory: Any = None,
) -> list[MetaboliteRow]:
    """Assemble one report row per panel, with cohort position where available.

    With ``sample_dir`` the fragment tables are read and each row learns which
    organisms its reading came from; with ``inventory`` those organisms are
    joined to the sample's own list. Direction follows the panel's reader
    direction, overridden per sample where the producers found say so.
    """
    from openbiota import drivers as drivers_mod

    rows: list[MetaboliteRow] = []

    for panel in panels:
        interpretation = panel.panel.interpretation
        direction = interpretation.effective_direction if interpretation else "unclear"
        direction_note = interpretation.direction_note if interpretation else ""
        drv = None
        if sample_dir is not None:
            drv = drivers_mod.for_panel(
                sample_dir, panel.panel.name,
                entry_ids=panel.panel.aggregate_target_ids(),
                inventory=inventory,
                exceptions=(interpretation.driver_exceptions if interpretation else None),
                default_direction=direction,
            )
            if drv is not None and drv.resolved_direction and drv.resolved_direction != direction:
                direction = drv.resolved_direction
                direction_note = (
                    f"Coloured as {direction} for this sample because {drv.resolved_reason}."
                )
        value = panel.copies_per_100_genomes
        panel_range = ranges.panels.get(panel.panel.name) if ranges else None

        percentile = (
            panel_range.percentile_of(value)
            if panel_range is not None and value is not None
            else None
        )
        # Nothing detected is zero, and zero belongs at the bottom of the
        # scale. `percentile_of` shares a tie block's midpoint, which is
        # right for a run of equal non-zero values and wrong here: on a
        # zero-inflated panel it put a sample with no reads at the 28th
        # percentile while the chip beside it read "not detected". Both
        # statements were about the same reading and they could not both be
        # true. A measured zero is the minimum of the scale, and that is
        # where it is now drawn and what the number now says.
        if percentile is not None and not value:
            percentile = 0.0
        # A percentile is only a comparison when there is a spread to
        # compare against. Where more than half the cohort carries none of
        # something, the median is zero and any detection at all lands
        # above it: crc_virulence read "above average" off four fragments
        # for exactly that reason, with none of its three organisms in the
        # sample. Recorded here rather than in the renderer so the file
        # says it too.
        cohort_mostly_absent = bool(
            panel_range is not None and not panel_range.median
        )
        status = status_for(
            percentile=percentile,
            higher_means=direction,
            detected=panel.detected,
            confidence=panel.confidence,
            stable=panel.stable,
            cohort_mostly_absent=cohort_mostly_absent,
        )

        extra: list[str] = []
        if panel.confidence == PROVISIONAL:
            why = "; ".join(
                reasons_below_confirmed(
                    max((g.fragments for g in panel.targets), default=0),
                    max(
                        (g.mean_identity for g in panel.targets if g.mean_identity),
                        default=None,
                    ),
                )
            )
            extra.append(
                "This is a provisional result rather than a confirmed one"
                + (f" — {why}." if why else ".")
                + " Read it as a possible signal."
            )
        elif not panel.stable and panel.detected:
            extra.append(
                f"Based on {panel.accepted_fragments} matching DNA fragments, few enough that "
                "the exact figure carries real uncertainty; the broad position is more "
                "reliable than the precise number."
            )
        # A panel can be "confirmed" overall while one of the genes setting its
        # headline is matching distant relatives. That halves the reliability
        # of the figure and has to be said.
        aggregate_ids = set(panel.panel.aggregate_target_ids())
        weak = [
            g
            for g in panel.targets
            if g.entry_id in aggregate_ids
            and g.fragments > 0
            and (g.mean_identity or 0) < MIN_MEAN_IDENTITY
        ]
        if weak and panel.confidence == CONFIRMED:
            names = ", ".join(g.gene or g.entry_id for g in weak)
            extra.append(
                f"One of the genes behind this figure ({names}) is matching distant relatives "
                f"rather than close ones, so treat the number as an upper bound. The pathway is "
                "still detected — other genes in it match well."
            )

        if panel.copies_per_100_genomes_residue_consistent is not None:
            extra.append(
                "A stricter version of this measurement, which additionally requires the "
                "matched gene to carry the correct active-site amino acid, gives "
                f"{panel.copies_per_100_genomes_residue_consistent:.1f} copies per 100 "
                "genomes. The true value lies between the two."
            )

        rows.append(
            MetaboliteRow(
                panel=panel.panel.name,
                metabolite=panel.panel.metabolite,
                value=value,
                percentile=percentile,
                cohort_median=panel_range.median if panel_range else None,
                cohort_p25=panel_range.percentiles.get(25) if panel_range else None,
                cohort_p75=panel_range.percentiles.get(75) if panel_range else None,
                detected=panel.detected,
                fragments=panel.accepted_fragments,
                status=status,
                what_it_is=interpretation.what_it_is if interpretation else panel.panel.description,
                made_from=interpretation.made_from if interpretation else "—",
                higher_means=direction,
                evidence_strength=interpretation.evidence_strength if interpretation else "limited",
                summary=interpretation.summary if interpretation else "",
                evidence_detail=interpretation.evidence_detail if interpretation else "",
                citation=interpretation.citation if interpretation else panel.panel.citation,
                confidence=panel.confidence,
                # The band sentences describe where a measured reading sits
                # ("less capacity than most people"). None of them is true
                # of a pathway whose headline genes produced no reads, and
                # printing one beside a NOT DETECTED chip reads as a
                # measured-low result. Withheld when nothing was detected.
                implication=(
                    interpretation.implication_for(percentile)
                    if interpretation and panel.detected
                    else ""
                ),
                genes=[
                    (g.gene or g.entry_id, g.fragments, g.copies_per_100_genomes)
                    for g in panel.targets
                ],
                extra_note=" ".join(extra),
                category=panel.panel.category,
                drivers=drv,
                direction_note=direction_note,
            )
        )

    return rows


def reference_json(
    rows: list[MetaboliteRow], ranges: ReferenceRanges | None
) -> dict[str, Any]:
    """Machine-readable form of the cohort comparison, for results.json."""
    return {
        "cohort": None if ranges is None else {
            "study": ranges.study,
            "description": ranges.description,
            "n_samples": ranges.n_samples,
            "reads_per_sample": ranges.reads_per_sample,
            "built_at": ranges.built_at,
        },
        "panels": {
            row.panel: {
                "value": row.value,
                "percentile": None if row.percentile is None else round(row.percentile, 1),
                "cohort_p25": row.cohort_p25,
                "cohort_median": row.cohort_median,
                "cohort_p75": row.cohort_p75,
                "status": row.status.label,
                "status_note": row.status.note,
                # The label of an unstable reading is its position and not a
                # verdict: the PDF says so by withholding the colour, which a
                # client reading this file cannot see. Said here in a field
                # rather than left in the prose of `status_note`, so nobody
                # has to parse a sentence to find out whether "above average"
                # was a finding or a handful of reads.
                "assessed": row.status.assessed,
                "confidence": row.confidence,
                "what_it_means": row.implication,
                "higher_means": row.higher_means,
                "evidence_strength": row.evidence_strength,
            }
            for row in rows
        },
    }
