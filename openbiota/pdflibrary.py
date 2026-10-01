"""Report pages for the spec-3 library: shape, ranked table, notable profiles,
microbial groups, species inventory, sequencing quality.

The old report had three profiles and gave each a page. With thirty-plus, a
page each would bury the reader; the spec (7.8) asks instead for the *shape*
of the result across the library — is one pattern standing out, or is the
sample high on everything because it is generally disturbed? — followed by a
ranked table and full pages only for the patterns worth reading about.

Everything here draws on the same flowables and the same seven-band scale as
the rest of the report.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import KeepTogether, PageBreak, Spacer, Table, TableStyle

from openbiota.pdflinks import (
    Anchor,
    LinkedTable,
    Paragraph,
    RowLink,
    back_link_line,
    group_dest,
    linked_cell,
    profile_dest,
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
    INK_FAINT,
    INK_SOFT,
    ORANGE,
    ORANGE_BG,
    PANEL_BG,
    RULE,
    SECTIONS,
    SLATE,
    SLATE_BG,
    Dot,
    PercentileBar,
    Rule,
    Status,
    StatusChip,
    Tile,
    _hex,
    _ordinal,
    fine_print,
    status_for,
)
from openbiota.pdfsummary import anchor_words, profile_status

#: Profiles at or above this percentile get a full detail page.
NOTABLE_PERCENTILE = 75.0
#: Never more than this many detail pages, most notable first.
MAX_DETAIL_PAGES = 8

MATURITY_SHORT = {
    "P0": "preclinical",
    "P1": "one small study",
    "P2": "human cohort",
    "P3": "replicated",
    "P4": "externally validated",
    "P5": "clinical-grade",
}

EVIDENCE_SHORT = {
    "TAX_REL": "species",
    "TAX_PRES": "presence",
    "GENE_ABUND": "genes",
    "CARRIER_ABUNDANCE": "carriers",
    "ECO": "diversity",
    "STRAIN": "strain",
    "PANCNV": "pangenome",
    "SEQ_HMM": "HMM",
    "PATH_ABUND": "pathways",
    "PATH_COMPL": "completeness",
    "CAZY": "CAZymes",
    "AMR": "resistome",
    "FLUX": "flux",
    "METAB": "metabolites",
    "MODEL": "model",
}

VERDICT_STATUS = {
    "specific": Status("one pattern stands out", ORANGE, ORANGE_BG, ""),
    "diffuse": Status("broadly raised", CORAL, CORAL_BG, ""),
    "mixed": Status("mixed picture", AMBER, AMBER_BG, ""),
    "quiet": Status("nothing stands out", GREEN, GREEN_BG, ""),
}


def _tiles(metrics: Sequence[tuple[str, str, str]], accent: colors.Color = ACCENT) -> Table:
    gap = 3 * mm
    width = (CONTENT_WIDTH - gap * (len(metrics) - 1)) / len(metrics)
    table = Table(
        [[Tile(value=v, label=k, note=n, width=width, accent=accent) for k, v, n in metrics]],
        colWidths=[width + gap] * (len(metrics) - 1) + [width],
        hAlign="LEFT",
    )
    table.setStyle(
        TableStyle(
            [
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    return table


def _plain_table(rows: list[list[Any]], widths: Sequence[float], *, header: bool = True,
                 zebra: bool = False) -> Table:
    table = LinkedTable(rows, colWidths=list(widths), hAlign="LEFT", repeatRows=1 if header else 0)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
    ]
    if header:
        style.append(("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE))
        style.append(("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE))
    if zebra:
        for i in range(1, len(rows)):
            if i % 2 == 0:
                style.append(("BACKGROUND", (0, i), (-1, i), PANEL_BG))
    table.setStyle(TableStyle(style))
    return table


def _badge(text: str, colour: colors.Color, st: dict[str, ParagraphStyle]) -> Paragraph:
    return Paragraph(
        f"<font color='{_hex(colour)}' size='6.2'><b>{text}</b></font>", st["cell"]
    )


def _maturity_colour(maturity: str) -> tuple[colors.Color, colors.Color]:
    return {
        "P0": (SLATE, SLATE_BG), "P1": (SLATE, SLATE_BG), "P2": (AMBER, AMBER_BG),
        "P3": (GREEN, GREEN_BG), "P4": (GREEN, GREEN_BG), "P5": (GREEN, GREEN_BG),
    }.get(maturity, (SLATE, SLATE_BG))


def _evidence_types(result: Any) -> list[str]:
    seen: list[str] = []
    for m in result.modules.values():
        et = m.module.evidence_type if m.module is not None else {
            "taxonomic": "TAX_REL", "functional": "GENE_ABUND", "ecological": "ECO",
            "phenotype": "phenotype",
        }.get(m.name, m.name)
        if et not in seen:
            seen.append(et)
    return seen


# --------------------------------------------------------------------------- #
# shape + ranked library
# --------------------------------------------------------------------------- #


def library_pages(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    similarity: Any | None,
    profile_validation: dict[str, Any] | None,
    plan: Any | None = None,
    copri_complex: dict[str, Any] | None = None,
) -> None:
    from openbiota.pdfatlas import evidence_cards

    story.extend(section_heading(section, "Resemblance to published disease patterns, in detail", st["h1"]))
    story.append(
        Paragraph(
            "Researchers have described how the gut communities of people with certain "
            "conditions differ, on average, from those of healthy people. This section asks, "
            "for each of the patterns in the library: <b>how closely does this sample resemble "
            "that published group-level pattern?</b> The answer is a percentile of a matched "
            "reference group, on the same scale as every other reading in this report. It is "
            "not a diagnosis and not a probability of disease.",
            st["body"],
        )
    )
    if similarity is None or not similarity.results:
        story.append(Spacer(1, 4 * mm))
        story.append(Paragraph("Profile similarity was not computed for this run.", st["body"]))
        return

    story.append(Spacer(1, 3 * mm))
    _ledger_block(story, st, similarity)
    story.append(Spacer(1, 3 * mm))
    _how_to_read(story, st)
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        "Every scored pattern follows, most resemblant first, each with the human research on what has been "
        "tried for that condition. Patterns the engine declined to score are listed with the reason at the end.",
        st["small"],
    ))

    for result in _detailed(similarity):
        story.append(PageBreak())
        _profile_page(story, st, result, profile_validation,
                      urticaria=getattr(similarity, "urticaria", None),
                      copri_complex=copri_complex)
        if plan is not None:
            evidence_cards(
                story, st, plan.cards_for("profile", result.profile.name),
                heading="What the research says about this condition and this pattern",
            )
    declined = [r for r in similarity.ranked if not r.reportable]
    if declined:
        story.append(PageBreak())
        story.extend(section_heading(section, "Patterns not scored, and why", st["h1"], sub=1))
        story.append(Paragraph(
            "A pattern is not scored when its own conditions are not met. Each line says "
            "what would need to change for the pattern to be scored.", st["body"]))
        for r in declined:
            if r.abstention.abstained:
                reason = "; ".join(r.abstention.triggered) or "declined"
                if r.abstention.remedies:
                    reason += ". Would need: " + "; ".join(r.abstention.remedies)
            else:
                reason = "needs a measurement engine this pipeline does not run"
            story.append(Paragraph(f"<b>{r.profile.label}</b> — {' '.join(str(reason).split())}", st["small_ink"]))


def _detailed(similarity: Any) -> list[Any]:
    """Every scored pattern gets a detail page (spec v04: no cap)."""
    return [r for r in similarity.ranked if r.reportable and r.combined_percentile is not None]


def _notable(similarity: Any) -> list[Any]:
    ranked = similarity.ranked
    picked = [
        r for r in ranked
        if r.reportable and r.combined_percentile is not None
        and r.combined_percentile >= NOTABLE_PERCENTILE
    ]
    for r in ranked:
        if r.profile.status == "validation_control" and r.reportable and r not in picked:
            picked.append(r)
    return picked


def _shape_block(story: list[Any], st: dict[str, ParagraphStyle], similarity: Any) -> None:
    shape = similarity.shape
    if shape is None:
        return
    status = VERDICT_STATUS.get(shape.verdict, Status(shape.verdict, SLATE, SLATE_BG, ""))
    head = Table(
        [[Paragraph("The shape of the result", st["h2"]),
          StatusChip(status, width=44 * mm, height=6.2 * mm, font_size=7.2)]],
        colWidths=[CONTENT_WIDTH - 46 * mm, 46 * mm], hAlign="LEFT",
    )
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                              ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    story.append(head)
    story.append(Paragraph(f"<b>{shape.headline}</b>", st["body_ink"]))
    story.append(Paragraph(" ".join(shape.explanation.split()), st["body"]))
    story.append(Spacer(1, 2 * mm))

    anchor = similarity.anchor
    anchor_text = "—"
    anchor_note = "not computed"
    if anchor is not None and anchor.valid:
        words, _ = anchor_words(anchor)
        anchor_text = f"{anchor.score:+.2f}"
        anchor_note = words
    spec = (
        f"{shape.mean_specificity_of_high:.2f}" if shape.mean_specificity_of_high is not None else "—"
    )
    story.append(
        _tiles(
            [
                ("PATTERNS SCORED", f"{shape.n_scored}", f"of {shape.n_profiles} in the library"),
                ("ABOVE TYPICAL", f"{shape.n_above_typical}",
                 f"above 75th pct · {shape.n_notably_high} notably high"),
                ("SPECIFICITY OF THOSE", spec,
                 "no high scorers" if not shape.high_scorers
                 else "part from challenge-cohort AUCs" if "AUROC" in shape.specificity_basis
                 else "literature prior only"),
                ("NOT SCORED", f"{shape.n_abstained + shape.n_not_computable}",
                 f"{shape.n_abstained} declined · {shape.n_not_computable} unmeasurable"),
                ("GENERAL GUT-HEALTH INDEX", anchor_text, anchor_note),
            ]
        )
    )


def _ranked_table(story: list[Any], st: dict[str, ParagraphStyle], similarity: Any) -> None:
    story.append(Paragraph("Every pattern in the library, most resemblant first", st["h2"]))
    story.append(
        Paragraph(
            "<b>Maturity</b> is how far the underlying research has come (P0 preclinical → P5 "
            "clinical-grade); <b>evidence</b> lists which kinds of measurement the pattern uses. "
            "A pattern marked <i>declined</i> did not meet its own conditions for a meaningful "
            "score; <i>not measurable</i> means it needs an engine this pipeline does not yet run. "
            "Every scored pattern has a full page in Part B, and its row here is a link to it.",
            st["small"],
        )
    )
    story.append(Spacer(1, 2 * mm))
    widths = [CONTENT_WIDTH * w for w in (0.30, 0.13, 0.28, 0.07, 0.22)]
    rows: list[list[Any]] = [[
        Paragraph("PATTERN", st["label"]),
        Paragraph("MATURITY · EVIDENCE", st["label"]),
        Paragraph("RESEMBLANCE", st["label"]),
        Paragraph("", st["label"]),
        Paragraph("READING", st["label"]),
    ]]
    with_pages = {r.profile.name for r in _detailed(similarity)}
    for result in similarity.ranked:
        profile = result.profile
        pct = None if not result.reportable else result.combined_percentile
        status = profile_status(result)
        if result.abstention.abstained:
            status = Status("declined", SLATE, SLATE_BG, "")
        elif result.combined_percentile is None:
            status = Status("not measurable", SLATE, SLATE_BG, "")
        elif result.status == "low_coverage":
            status = Status(f"{status.label} · low coverage", status.colour, status.background, "")
        m_col, m_bg = _maturity_colour(profile.evidence_maturity)
        types = " · ".join(EVIDENCE_SHORT.get(t, t) for t in _evidence_types(result))
        extra = ""
        if result.competing:
            rivals = ", ".join(n for n, _ in result.competing[:2]).replace("_", " ")
            extra = f" &nbsp;<font color='{_hex(INK_FAINT)}' size='6'>competes with {rivals}</font>"
        name_cell = Paragraph(
            f"<b>{profile.label}</b>"
            + (f"<br/><font color='{_hex(INK_FAINT)}' size='6.2'>{profile.disease or ''}"
               f"{' · opt-in' if profile.opt_in else ''}"
               f"{' · validation control' if profile.status == 'validation_control' else ''}</font>"
               if (profile.disease or profile.opt_in or profile.status == 'validation_control') else "")
            + extra,
            st["cell"],
        )
        # Where page 1 lands, and the way on to the pattern's own page.
        marker = RowLink(
            href=profile_dest(profile.name, detail=True) if profile.name in with_pages else None,
            anchor=profile_dest(profile.name),
        )
        rows.append([
            linked_cell(marker, name_cell),
            Paragraph(
                f"<font color='{_hex(m_col)}'><b>{profile.evidence_maturity}</b></font> "
                f"<font size='6.2'>{MATURITY_SHORT.get(profile.evidence_maturity, '')}</font>"
                f"<br/><font color='{_hex(INK_FAINT)}' size='6.2'>{types}</font>",
                st["cell"],
            ),
            PercentileBar(width=widths[2] - 6 * mm, percentile=pct, higher_means="adverse",
                          marker_colour=status.colour, height=6.0 * mm),
            Paragraph(f"<font color='{_hex(status.colour)}'>{_ordinal(pct)}</font>", st["pct"]),
            StatusChip(status, width=widths[4] - 4 * mm, font_size=6.4),
        ])
    story.append(_plain_table(rows, widths))
    if similarity.withheld:
        story.append(Spacer(1, 1.5 * mm))
        story.append(
            Paragraph(
                f"Not scored unless requested: {', '.join(p.label for p in similarity.withheld)}. "
                "These cover psychiatric and other sensitive conditions and are opt-in.",
                st["fine"],
            )
        )


def _ledger_block(story: list[Any], st: dict[str, ParagraphStyle], similarity: Any) -> None:
    ledger = similarity.ledger
    if ledger is None:
        return
    story.append(Paragraph("What we knew about you when scoring", st["h2"]))
    story.append(
        Paragraph(
            "Several patterns depend on facts the sequencing cannot see — medication, how long "
            "you have been ill, recent antibiotics. Those supplied are listed; those not supplied "
            "are listed with what they would have changed. Where a pattern needs one of them, the "
            "ordinary case was assumed (no medication, no recent antibiotics) and the pattern was "
            f"scored on that stated assumption rather than withheld; section {SECTIONS['guide']} lists every assumption.",
            st["small"],
        )
    )
    story.append(Spacer(1, 1.5 * mm))
    rows: list[list[Any]] = [[Paragraph("FACTOR", st["label"]), Paragraph("STATUS", st["label"]),
                              Paragraph("EFFECT ON SCORING", st["label"])]]
    for _key, label, value in ledger.recorded:
        rows.append([Paragraph(label, st["cell"]),
                     Paragraph(f"<font color='{_hex(GREEN)}'>recorded</font>: {value}", st["cell"]),
                     Paragraph("used where a pattern declares it", st["small"])])
    for _key, label, effect in ledger.assumed:
        rows.append([Paragraph(label, st["cell"]),
                     Paragraph(f"<font color='{_hex(AMBER)}'>not supplied · assumed absent</font>", st["cell"]),
                     Paragraph(effect, st["small"])])
    for _key, label, effect in ledger.missing:
        rows.append([Paragraph(label, st["cell"]),
                     Paragraph(f"<font color='{_hex(AMBER)}'>not supplied</font>", st["cell"]),
                     Paragraph(effect, st["small"])])
    if len(rows) > 1:
        story.append(_plain_table(rows, [CONTENT_WIDTH * 0.26, CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.44]))
    for note in ledger.notes:
        story.append(Paragraph(note, st["fine"]))


#: Profile name -> the urticaria (hives) panel scored on the source study's
#: own fixed scale. Both panels measure the same species this report already
#: ranks; the second scale exists for its coverage and censoring bounds.
URTICARIA_PANELS: Final = {
    "csu": "CSU_STOOL_WGS_V1",
    "symptomatic_dermographism": "SD_STOOL_TRANSLATED_V1",
}


def _urticaria_index_block(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    result: Any,
    urticaria: Mapping[str, Any] | None,
) -> None:
    """The study-scale index for a hives panel, beside the cohort percentile.

    Two summaries of one set of measurements, reconciled out loud: the
    percentile above is this sample's rank in a matched reference group, and
    this is the source study's own fixed 0-100 scale, which carries coverage
    and censoring bounds that a rank cannot.
    """
    panel_id = URTICARIA_PANELS.get(result.profile.name)
    cu = (urticaria or {}).get(panel_id) if panel_id else None
    if cu is None:
        return

    from openbiota.pdfatlas import AlignmentBar

    story.append(Spacer(1, 3.5 * mm))
    story.append(Paragraph("The same species on the source study's own scale", st["h3"]))

    if not cu.computed:
        states = "; ".join(
            f"<i>{r.feature.source_id.replace('_', ' ')}</i> — {r.state.replace('_', ' ')}"
            for r in cu.rows
        )
        story.append(Paragraph(
            f"<b>Not computed on this scale</b> ({cu.status.replace('_', ' ')}). No number is "
            "shown rather than a zero, because a zero here would be a measurement that was never "
            f"made. Panel features and why: {states}.",
            st["small"],
        ))
        story.append(Paragraph(
            "The percentile above still has a value because the two calculations treat a "
            "nondetection differently, and the difference is the point of showing both. The "
            "percentile places an undetected organism at a fixed floor and ranks the sample from "
            "there, which is how every reading in this report is built and what makes them "
            "comparable with each other. This scale declines to convert an undetected organism "
            "into a number at all unless a validated detection limit exists for it, and none "
            "does at this sequencing depth. Read the percentile as the position of this "
            "community, and read the blank here as the reason not to put much weight on it.",
            st["fine"],
        ))
        return

    lo, hi = cu.full_panel_bounds
    story.append(AlignmentBar(width=CONTENT_WIDTH, value=cu.score, bounds=(lo, hi)))
    story.append(Paragraph(
        f"<b>{cu.score:.0f} out of 100</b> — with no percent sign, because this is not a "
        f"percentage, a probability or a match score. Coverage {cu.usable_count} of "
        f"{cu.panel_count} panel features; the pale band, {lo:.0f} to {hi:.0f}, is where the "
        "figure could sit once the unmeasured features are accounted for. Fifty is neutral "
        "signed deviation under this transformation, not a definition of health, and higher "
        "means stronger alignment with the directions the source study reported.",
        st["small"],
    ))
    if cu.sensitivity is not None and cu.sensitivity.computed:
        s = cu.sensitivity
        story.append(Paragraph(
            f"Strict subset, the {s.panel_count} species that survive multiple-testing "
            f"correction in the source table: <b>{s.score:.0f}</b> at coverage "
            f"{s.usable_count} of {s.panel_count}. This is a sensitivity analysis on the same "
            "measurements, not a second independent test — the two panels share every one of "
            "these features, so a gap between them means the wider figure is being carried by "
            "its exploratory members.",
            st["small"],
        ))
    story.append(Paragraph(
        "These are coverage and censoring bounds, not a 95% confidence interval: they widen "
        "when a feature cannot be measured and say nothing about sampling error, about how "
        "well a single-centre discovery cohort transports to this population, or about how "
        "specific the pattern is to this condition. The percentile above and the figure here read "
        "the same measurements two ways, and are not two opinions about one quantity: the "
        "percentile is a rank, so it says how this sample compares with a matched reference "
        "group and it moves when that group's spread changes; this figure is a fixed "
        "transformation, so it says how far each organism sits from its reference centre in the "
        "direction the study reported, and it stays put. A high rank can therefore sit beside a "
        "middling figure, or the reverse, depending on how tightly the reference group clusters "
        "and on which features could be measured at all.",
        st["fine"],
    ))
    unusable = [r for r in cu.rows if not r.usable]
    if unusable:
        story.append(Paragraph(
            "Features outside the figure: "
            + "; ".join(
                f"<i>{r.feature.source_id.replace('_', ' ')}</i> ({r.state.replace('_', ' ')})"
                for r in unusable
            )
            + ". None of them was scored as zero.",
            st["fine"],
        ))


def _how_to_read(story: list[Any], st: dict[str, ParagraphStyle]) -> None:
    story.append(Rule(CONTENT_WIDTH))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("How the pages that follow are built", st["h2"]))
    story.append(
        Paragraph(
            "Each pattern is made of one or more <b>study modules</b> — one per independent "
            "published cohort, so a finding replicated by three groups counts three times and a "
            "finding from one group once. Modules are scored separately so you can see which "
            "line of evidence is driving the result, then fused with equal weight per "
            "independent study. Below each result is a table of every contributing feature with "
            "its measured value, its position in the comparison group, the direction the "
            "research reports, and its weight, so the score can be recomputed by hand.",
            st["body"],
        )
    )


# --------------------------------------------------------------------------- #
# one profile
# --------------------------------------------------------------------------- #


def _module_hint(module_score: Any) -> str:
    m = module_score.module
    if m is None:
        return {
            "taxonomic": "which species are present",
            "functional": "gene capacity from Part A",
            "ecological": "diversity and richness",
            "phenotype": "symptom-domain features",
        }.get(module_score.name, "")
    parts = [EVIDENCE_SHORT.get(m.evidence_type, m.evidence_type)]
    if m.claim_level and m.claim_level != "observed":
        parts.append(m.claim_level.replace("_", " "))
    if m.assay_transport and m.assay_transport != "native":
        parts.append("16S study" if m.assay_transport == "amplicon_taxon_to_shotgun" else m.assay_transport.replace("_", " "))
    if m.group:
        parts.append(f"study {m.group}")
    # When a module has withdrawn its claim because the organism the finding is
    # about is absent, that is the most important thing on the row — a reader
    # seeing "not scored" needs to know it means "the organism is not there",
    # not "we could not measure it".
    if getattr(module_score, "anchor_absent", False):
        anchors = ", ".join(n.replace("_", " ") for n in m.anchor_features)
        parts.append(f"<b>{anchors} not detected, so this pattern is not asserted</b>")
    return " · ".join(parts)


def _copri_block(
    story: list[Any], st: dict[str, ParagraphStyle], resolved: dict[str, Any]
) -> None:
    """Clade-level composition of the Segatella (Prevotella) copri complex."""
    story.append(Paragraph("The Prevotella copri complex, clade by clade", st["h2"]))
    status = str(resolved.get("status") or "")
    if status == "not_assessed":
        story.append(Paragraph(str(resolved.get("reason") or ""), st["fine"]))
        return

    story.append(Paragraph(str(resolved.get("plain") or ""), st["small"]))
    present = resolved.get("clades_present") or []
    if present:
        rows: list[list[Any]] = [[
            Paragraph("<b>CLADE</b>", st["cell"]), Paragraph("<b>SPECIES</b>", st["cell"]),
            Paragraph("<b>SGB</b>", st["cell"]), Paragraph("<b>OF COMMUNITY</b>", st["cell"]),
        ]]
        for clade in present:
            rows.append([
                Paragraph(f"<b>{clade['clade']}</b>", st["cell"]),
                Paragraph(f"<i>{clade['species']}</i>", st["cell"]),
                Paragraph(str(clade["sgb"]), st["cell"]),
                Paragraph(f"{clade['percent']:.3f}%", st["cell"]),
            ])
        story.append(_plain_table(
            rows,
            [CONTENT_WIDTH * 0.12, CONTENT_WIDTH * 0.46, CONTENT_WIDTH * 0.18,
             CONTENT_WIDTH * 0.24],
        ))
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font color='{_hex(INK_FAINT)}'>{resolved.get('disease_direction') or ''}</font>",
        st["fine"],
    ))
    if sources := resolved.get("sources"):
        story.append(Paragraph(
            f"<font color='{_hex(INK_FAINT)}' size='6.6'>{'; '.join(sources)}.</font>",
            st["fine"],
        ))
    story.append(Spacer(1, 3 * mm))


def _profile_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    result: Any,
    profile_validation: dict[str, Any] | None,
    urticaria: Any | None = None,
    copri_complex: dict[str, Any] | None = None,
) -> None:
    profile = result.profile
    status = profile_status(result)
    combined_pct = None if result.abstention.abstained else result.combined_percentile

    # Where the at-a-glance row leads; the small link takes the reader back.
    story.append(Anchor(profile_dest(profile.name, detail=True), above=6 * mm))
    story.append(back_link_line(profile_dest(profile.name), "at a glance"))
    head = Table(
        [[Paragraph(
            f"{profile.label}",
            st["h1"],
          ),
          StatusChip(status, width=36 * mm, height=6.2 * mm, font_size=7.2)]],
        colWidths=[CONTENT_WIDTH - 38 * mm, 38 * mm], hAlign="LEFT",
    )
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                              ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    story.append(head)
    m_col, _ = _maturity_colour(profile.evidence_maturity)
    story.append(
        Paragraph(
            f"Pattern version <b>{profile.version}</b> <font color='#8C99A6'>({profile.name})</font> &middot; "
            f"maturity <font color='{_hex(m_col)}'><b>{profile.evidence_maturity}</b></font> "
            f"({profile.maturity_meaning}) &middot; {profile.status.replace('_', ' ')}"
            + (f" &middot; {profile.population}" if profile.population else ""),
            st["small"],
        )
    )
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(" ".join(profile.summary.split()), st["body"]))
    story.append(Spacer(1, 2 * mm))

    # headline
    interval = ""
    if result.vector is not None and result.vector.interval is not None:
        lo, hi = result.vector.interval
        interval = f"<br/><font color='{_hex(INK_FAINT)}' size='6.5'>likely range {lo:.0f}–{hi:.0f}</font>"
    headline = Table(
        [[
            Paragraph("<b>COMBINED RESEMBLANCE</b><br/>"
                      f"<font color='{_hex(INK_FAINT)}' size='7'>percentile of the matched "
                      "reference group</font>", st["cell"]),
            PercentileBar(width=CONTENT_WIDTH * 0.46, percentile=combined_pct, higher_means="adverse",
                          marker_colour=status.colour, height=9 * mm, show_scale=True, track=4.2 * mm),
            Paragraph((f"<font color='{_hex(status.colour)}'>{_ordinal(combined_pct)}</font>{interval}"
                       if combined_pct is not None else "—"), st["metric"]),
        ]],
        colWidths=[CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.52, CONTENT_WIDTH * 0.18], hAlign="LEFT",
    )
    headline.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
        ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("ROUNDEDCORNERS", [4, 4, 4, 4]),
    ]))
    story.append(headline)

    # A withdrawn module changes what the number above means, so the
    # qualification belongs against the number rather than further down the
    # page. Without this, a profile whose defining organism is absent still
    # prints a high percentile and a "HIGH" chip with nothing beside them.
    scored_modules = list(result.modules.values())
    withdrawn = [m for m in scored_modules if getattr(m, "anchor_absent", False)]
    if withdrawn and combined_pct is not None:
        remaining = [
            m for m in scored_modules
            if m.reportable and result.weights.get(m.name, 0.0) > 0
        ]
        anchors = ", ".join(
            n.replace("_", " ")
            for m in withdrawn
            for n in (m.module.anchor_features if m.module is not None else ())
        )
        rests_on = (
            f"one study module ({(remaining[0].module.label if remaining[0].module else '') or remaining[0].name})"
            if len(remaining) == 1
            else f"{len(remaining)} study modules"
        )
        # Imported here because this function rebinds AMBER locally further
        # down, which would otherwise shadow the module-level import.
        from openbiota.pdfreport import AMBER as _AMBER

        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(
            f"<font color='{_hex(_AMBER)}'><b>Read this percentile narrowly.</b></font> "
            f"<i>{anchors}</i> — the organism this profile's central finding is about — was "
            f"looked for and not detected, so the module carrying that finding is not scored. "
            f"The percentile above therefore rests on {rests_on}, and reflects accessory taxa "
            f"rather than the organism the condition is associated with.",
            st["small"],
        ))
    story.append(Spacer(1, 3 * mm))

    if result.abstention.abstained:
        story.append(Paragraph("No score was produced", st["h2"]))
        story.append(Paragraph("The engine <b>declined to score this pattern</b>: one or more "
                               "conditions for a meaningful score were not met.", st["body"]))
        for reason in result.abstention.triggered:
            story.append(Paragraph(f"•&nbsp;&nbsp;{reason}", st["body"]))
        story.append(Paragraph("<b>What would change this</b>", st["small"]))
        for remedy in result.abstention.remedies:
            story.append(Paragraph(f"•&nbsp;&nbsp;{remedy}", st["body"]))
        story.append(Spacer(1, 2 * mm))
    elif result.abstention.assumptions:
        # Scored on a stated assumption: the fact the pattern depends on was
        # not supplied, so the ordinary case was assumed and the score stands
        # with the assumption beside it and the consequence of its being wrong.
        from openbiota.context import FIELD_BY_KEY
        from openbiota.pdfreport import AMBER, AMBER_BG

        lines: list[Any] = [Paragraph(
            f"<font color='{_hex(AMBER)}'><b>Scored on an assumption.</b></font> "
            "This pattern depends on a fact about you that was not supplied, so the report assumed the "
            "ordinary case rather than withhold the score:", st["small"])]
        for key in result.abstention.assumptions:
            f = FIELD_BY_KEY.get(key)
            label = f.label if f is not None else key.replace("_", " ")
            assumed = {"medication": "not taken", "exposure": "none", "safety": "no"}.get(
                f.kind if f is not None else "", "absent")
            consequence = f.if_present if f is not None else "the score would need to be re-read."
            lines.append(Paragraph(
                f"•&nbsp;&nbsp;<b>{label[:1].upper()}{label[1:]}: assumed {assumed}.</b> "
                f"If that is wrong: {consequence}", st["small"]))
        lines.append(Paragraph(
            "Supplying the fact (for example <font face='Courier'>--medications</font>) replaces the "
            "assumption and re-scores this pattern.", st["fine"]))
        box = Table([[lines]], colWidths=[CONTENT_WIDTH], hAlign="LEFT")
        box.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), AMBER_BG), ("LINEBEFORE", (0, 0), (0, -1), 1.6, AMBER),
            ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm), ("RIGHTPADDING", (0, 0), (-1, -1), 2.5 * mm),
            ("TOPPADDING", (0, 0), (-1, -1), 1.6 * mm), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6 * mm),
        ]))
        story.append(box)
        story.append(Spacer(1, 2 * mm))

    # modules
    _urticaria_index_block(story, st, result, urticaria)

    story.append(Paragraph("What is driving it", st["h2"]))
    rows: list[list[Any]] = []
    modules = sorted(result.modules.values(),
                     key=lambda m: -result.weights.get(m.name, m.weight_in_combined))
    for module in modules:
        label = (module.module.label if module.module is not None and module.module.label
                 else module.name.replace("_", " ").title())
        pct = module.percentile if module.reportable else None
        w = result.weights.get(module.name, module.weight_in_combined)
        share = f"{w:.0%} of combined" if w else "not in combined"
        if module.module is not None and not module.module.bound:
            share = f"not measurable — {module.module.would_bind_if or 'engine missing'}"
        m_status = status_for(percentile=pct, higher_means="adverse") if pct is not None else Status(
            module.status.replace("_", " "), SLATE, SLATE_BG, "")
        rows.append([
            Paragraph(f"<b>{label}</b><br/><font color='{_hex(INK_FAINT)}' size='7'>"
                      f"{_module_hint(module)} &middot; {module.n_measured}/{len(module.features)} "
                      f"features &middot; {share}</font>", st["cell"]),
            PercentileBar(width=CONTENT_WIDTH * 0.34, percentile=pct, higher_means="adverse",
                          marker_colour=m_status.colour, height=6.4 * mm),
            Paragraph(f"<font color='{_hex(m_status.colour)}'>{_ordinal(pct)}</font>" if pct is not None else "—",
                      st["pct"]),
            StatusChip(m_status, width=CONTENT_WIDTH * 0.17 - 4 * mm, font_size=6.4),
        ])
    if rows:
        table = Table(rows, colWidths=[CONTENT_WIDTH * 0.40, CONTENT_WIDTH * 0.36, CONTENT_WIDTH * 0.07,
                                       CONTENT_WIDTH * 0.17], hAlign="LEFT")
        table.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3), ("LEFTPADDING", (0, 0), (-1, -1), 0),
            ("LEFTPADDING", (3, 0), (3, -1), 4), ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
        ]))
        story.append(table)
    story.append(Spacer(1, 3 * mm))

    # The Prevotella/Segatella copri complex, where it is relevant. "P. copri"
    # is 13 species-level clades 13-21% apart, so a species-level number could
    # never say which was present. The extended lane resolves all of them, and
    # this is the profile where that matters.
    if copri_complex and profile.name in {"ra", "t2d", "obesity"}:
        _copri_block(story, st, copri_complex)

    # confidence vector + competition
    _confidence_block(story, st, result)

    # decomposition
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("Every feature that produced this score", st["h2"]))
    story.append(Paragraph(
        "Direction is what the research reports for the condition (▲ higher, ▼ lower). "
        "<b>v</b> is how far this sample moves in that direction, from −1 (opposite) to +1 "
        "(strongly concordant). <b>w</b> is the feature's weight in its module; q, s, c "
        "(cohort independence, specificity, cross-study agreement) are recorded for audit and "
        "not multiplied into the score.", st["small"]))
    story.append(Paragraph(
        f"Every species value here is the reference catalogue's own measurement (MetaPhlAn 3, the catalogue the reference cohort was profiled with), because that is what makes a percentile against that cohort mean anything. It is not the organism's share of the community in section {SECTIONS['catalogue']}, which pools nine detection methods and the current catalogues: a species that catalogue has no genome for reads as <i>absent</i> here and still has a share there; where that is so, the pooled share is printed beneath. Section {SECTIONS['catalogue']} is the account of what is present; this table is the account of what produced this score.",
        st["small"]))
    story.append(Spacer(1, 1.5 * mm))
    _decomposition_table(story, st, modules)

    notes: list[str] = []
    if profile.status == "validation_control":
        notes.append("<b>Validation control.</b> Included to check that the scoring engine recovers a "
                     "known signal on labelled data. Not a cancer screen and not a risk estimate; "
                     "validated clinical tests exist for this condition.")
    if result.abstention.abstained:
        notes.append("<b>Confidence:</b> not graded — no score was produced.")
    else:
        reasons = "; ".join(result.confidence.reasons[:3])
        more = max(0, len(result.confidence.reasons) - 3)
        notes.append(f"<b>Confidence grade {result.confidence.grade.lower()}.</b> "
                     + (f"Why not higher: {reasons}"
                        + (f"; +{more} more, all in results.json." if more else ".")
                        if reasons else "All confidence checks passed."))
    if profile.fusion_note:
        notes.append(" ".join(profile.fusion_note.split()))
    if profile_validation:
        notes.append(_validation_text(profile.name, profile_validation))
    notes.extend(" ".join(c.split()) for c in profile.caveats)
    notes.append(f"Sources: {' '.join(profile.citation.split())}")
    fine_print(story, st, notes)


def _confidence_block(story: list[Any], st: dict[str, ParagraphStyle], result: Any) -> None:
    vec = result.vector
    if vec is None:
        return
    story.append(Paragraph("How much to trust it", st["h2"]))
    story.append(Paragraph(
        "Confidence is reported beside the score, never multiplied into it. Each row answers a "
        "different question; a weak row means 'less certain', not 'lower'.", st["small"]))

    def frac(v: float | None) -> str:
        return "—" if v is None else f"{v:.2f}"

    def dot_for(v: float | None, good: float = 0.7, ok: float = 0.4) -> colors.Color:
        if v is None:
            return SLATE
        return GREEN if v >= good else AMBER if v >= ok else CORAL

    rows = [
        ("Technical reliability", frac(vec.technical_reliability), dot_for(vec.technical_reliability),
         "sequencing depth, host fraction, reference match"),
        ("Independent evidence coverage", frac(vec.independent_feature_coverage),
         dot_for(vec.independent_feature_coverage),
         "share of the pattern's independent evidence clusters this sample could measure"),
        ("Evidence maturity", vec.evidence_maturity, _maturity_colour(vec.evidence_maturity)[0], vec.maturity_meaning),
        ("Assay transportability", vec.assay_transportability.replace("_", " "),
         GREEN if vec.assay_transportability == "native" else AMBER,
         "whether the studies measured what this sample measured"),
        ("Population match", "unknown" if vec.population_transportability is None else frac(vec.population_transportability),
         dot_for(vec.population_transportability), "how well the reference group matches you"),
        ("Specificity", "unknown" if vec.empirical_specificity is None else frac(vec.empirical_specificity),
         dot_for(vec.empirical_specificity),
         "does this pattern separate its condition from general disturbance? measured where labelled data exist"),
        ("Direction conflict", vec.direction_conflict,
         GREEN if vec.direction_conflict == "none" else AMBER if vec.direction_conflict == "minor" else CORAL,
         "studies disagreeing on direction for one or more features"),
    ]
    table_rows: list[list[Any]] = []
    for label, value, colour, note in rows:
        table_rows.append([
            Table([[Dot(colour, size=2.2 * mm), Paragraph(f"<b>{label}</b>", st["cell"])]],
                  colWidths=[4 * mm, CONTENT_WIDTH * 0.30 - 4 * mm], hAlign="LEFT",
                  style=TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                    ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0)])),
            Paragraph(value, st["cell"]),
            Paragraph(note, st["small"]),
        ])
    story.append(_plain_table(table_rows, [CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.14, CONTENT_WIDTH * 0.56],
                              header=False))
    if result.competing:
        comp = ", ".join(f"{n.replace('_', ' ')} ({_ordinal(p)})" for n, p in result.competing[:4])
        margin = f" Margin over the next pattern: {result.top_two_margin:.0f} points." if result.top_two_margin is not None else ""
        story.append(Spacer(1, 1 * mm))
        story.append(Paragraph(f"<b>Competing patterns</b> scoring nearly as high: {comp}.{margin} "
                               "Where several patterns share features, a high score on one is not "
                               "evidence for it over the others.", st["small"]))
    for note in vec.notes[:4]:
        story.append(Paragraph(note, st["fine"]))


def _decomposition_table(story: list[Any], st: dict[str, ParagraphStyle], modules: Sequence[Any]) -> None:
    from openbiota import pdfcontext as _pdfcontext

    head = [Paragraph(x, st["label"]) for x in ("FEATURE", "MEASURED", "REF PCT", "DIR", "w", "q", "s", "c", "v", "SOURCES")]
    rows: list[list[Any]] = [head]
    spans: list[int] = []
    for module in modules:
        label = (module.module.label if module.module is not None and module.module.label
                 else module.name.replace("_", " ").title())
        spans.append(len(rows))
        rows.append([Paragraph(f"<b>{label}</b> <font color='{_hex(INK_FAINT)}'>· {_module_hint(module)}</font>", st["cell"])] + [""] * 9)
        for f in module.features:
            engine = f.feature.engine
            if f.raw_value is None:
                # "absent" asserts the organism is not there. A CLR-imputed
                # non-detection supports only "not detected at this depth":
                # sample_clr imputes every cohort taxon at the floor, so an
                # undetected species arrives here measured with no value,
                # at any depth and whatever the QC gate said.
                measured = (
                    "missing" if getattr(f, "missing", False) else "not detected"
                )
                pooled = (_pdfcontext.pooled_note(f.feature.name, catalogue_value=None)
                          if engine in ("metaphlan", "carrier") else None)
                if pooled:
                    # this catalogue read the species as absent; pooling every
                    # method found it, and the organism list gives it that share
                    measured += f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>{pooled}</font>"
            elif engine in ("metaphlan", "carrier"):
                measured = f"{f.raw_value:.3f}%"
            elif engine == "diamond":
                measured = f"{f.raw_value:.1f}/100"
            else:
                measured = f"{f.raw_value:.2f}"
            v = f.v
            v_colour = (_hex(CORAL) if v is not None and v > 0.33 else _hex(GREEN) if v is not None and v < -0.33 else _hex(SLATE))
            sources = ", ".join(f.feature.sources)
            name = f.feature.name.replace("_", " ")
            if f.feature.level == "genus":
                name += " (genus)"
            if f.feature.measured_as:
                # The source paper's name is the identifier; this pipeline's
                # profiler reports it under another. Print both, always.
                name += f", measured as {f.feature.measured_as.replace('_', ' ')}"
            rows.append([
                Paragraph(f"<i>{name}</i>" if engine in ("metaphlan",) else name, st["cell"]),
                Paragraph(measured, st["cell"]),
                Paragraph(f"{f.percentile:.0f}" if f.percentile is not None else "—", st["cell"]),
                Paragraph("▲" if f.feature.direction == "increased" else "▼", st["cell"]),
                Paragraph(f"{f.feature.w:.2f}", st["cell"]),
                Paragraph(f"{f.feature.q:.2f}", st["cell"]),
                Paragraph(f"{f.feature.s:.2f}", st["cell"]),
                Paragraph(f"{f.feature.c:.2f}", st["cell"]),
                Paragraph(f"<font color='{v_colour}'><b>{v:+.2f}</b></font>" if v is not None else "—", st["cell"]),
                Paragraph(sources, st["fine"]),
            ])
    widths = [0.24, 0.10, 0.07, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.27]
    table = Table(rows, colWidths=[CONTENT_WIDTH * w for w in widths], hAlign="LEFT", repeatRows=1)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("TOPPADDING", (0, 0), (-1, -1), 1.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6), ("LEFTPADDING", (0, 0), (-1, -1), 1.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1.5), ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
    ]
    for i in spans:
        style.append(("SPAN", (0, i), (-1, i)))
        style.append(("BACKGROUND", (0, i), (-1, i), PANEL_BG))
    table.setStyle(TableStyle(style))
    story.append(table)


def _validation_text(name: str, validation: dict[str, Any]) -> str:
    results = [r for r in validation.get("results", []) if r.get("profile") == name]
    if not results:
        return ("<b>Measured discrimination: not available.</b> No labelled cohort exists for "
                "this pattern, so its ability to separate cases from controls is unmeasured.")
    own = next((r for r in results if r["condition"].lower() == name.lower()), None)
    parts = []
    for r in sorted(results, key=lambda x: -(x.get("auc") or 0)):
        ci = r.get("auc_ci95")
        parts.append(f"{r['condition']} {r['auc']:.2f}" + (f" ({ci[0]:.2f}–{ci[1]:.2f})" if ci else ""))
    verdict = validation.get("specificity_verdicts", {}).get(name, "")
    return ("<b>Measured on labelled cohorts</b> (AUC, 0.5 = chance, 1.0 = perfect): " + "; ".join(parts) + ". "
            + (f"{verdict[0].upper() + verdict[1:]}." if verdict else "")
            + (" The published 0.80+ for this condition comes from learned models over hundreds of features; "
               "a handful of prespecified markers cannot match it, and is not meant to."
               if own and (own.get("auc") or 0) < 0.75 else ""))


# --------------------------------------------------------------------------- #
# microbial groups + species inventory
# --------------------------------------------------------------------------- #


def community_groups_pages(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    community: Any | None,
    cohort_note: str,
    plan: Any | None = None,
    inv: Any | None = None,
    carded: frozenset[str] = frozenset(),
) -> None:
    from openbiota.taxongroups import GROUP_CATEGORIES

    story.extend(section_heading(section, "Your microbial groups in detail", st["h1"]))
    if community is None:
        story.append(Paragraph("The taxonomic engine did not run; no group readings.", st["body"]))
        return
    story.append(
        Paragraph(
            "Curated sets of species that share a job or an origin — the fibre fermenters, the "
            "mouth bacteria that should not be in stool in quantity, the organisms that make "
            "methane or sulfide — each summed to one reading and placed against "
            f"{cohort_note}. The arrow says which direction the research reads as favourable; "
            "no arrow means it depends on context. Under each group are the species that made "
            "up your reading.",
            st["body"],
        )
    )
    story.append(Spacer(1, 2 * mm))

    by_cat = community.by_category()
    for key, heading in GROUP_CATEGORIES:
        groups = by_cat.get(key, [])
        if not groups:
            continue
        ordered = sorted(groups, key=lambda r: -(abs((r.level_percentile or 50) - 50)))
        for i, g in enumerate(ordered):
            # The category heading travels with its first card so it is never
            # left alone at the foot of a page.
            lead = [Spacer(1, 1.5 * mm), Paragraph(heading, st["h2"])] if i == 0 else None
            _group_card(
                story, st, g, lead=lead,
                cards=None if plan is None else plan.cards_for("group", g.group.name),
                coherence=None if plan is None else plan.coherence.by_group(g.group.name),
                inv=inv, carded=carded,
            )

    story.append(PageBreak())



def _arrow(higher_means: str) -> str:
    return {
        "adverse": f' <font color="{_hex(CORAL)}" size="6.5">\u25b2</font>',
        "favourable": f' <font color="{_hex(GREEN)}" size="6.5">\u25bc</font>',
    }.get(higher_means, "")


# --------------------------------------------------------------------------- #
# the organisms behind a group reading
# --------------------------------------------------------------------------- #


class GroupDriver:
    """One member species of a group, with everything the page needs: its
    share of the group's reads, its own abundance and cohort position, its
    class, and the destination of its own card if the report has one.

    Built once per group by :func:`group_drivers` and used by both the
    glance row and the detail card, so the two cannot name different
    organisms for the same reading.
    """

    __slots__ = ("cls", "counted_within", "dest", "detected", "deviation", "display", "expansion_only", "in_primary",
                 "percent", "percentile", "prevalence", "sample_percent", "share", "species")

    def __init__(self, m: Any, *, total: float, inv: Any | None, carded: frozenset[str]) -> None:
        from openbiota import organisms as org
        from openbiota.inventory import Organism
        from openbiota.pdflinks import organism_dest

        self.species: str = m.species
        self.detected: bool = bool(m.detected and m.percent > 0)
        #: Abundance in the lane the group was summed in; the share below
        #: is computed from it so the parts add up to the group's reading.
        self.percent: float = float(m.percent or 0.0)
        self.share: float = (self.percent / total) if (total > 0 and self.detected) else 0.0
        o = inv.get(m.species) if inv is not None else None
        #: Abundance as the rest of the report states it - the inventory's
        #: one share of the composition - so an organism has one number
        #: everywhere. A population counted within a relative's share has
        #: none of its own, and is said so rather than given a lane reading.
        self.in_primary: bool = bool(o.in_primary) if o is not None else True
        self.counted_within: str | None = getattr(o, "counted_within", None) if o is not None else None
        self.sample_percent: float = (
            float(o.percent) if (o is not None and o.in_primary)
            else 0.0 if o is not None else self.percent
        )
        self.display: str = o.display if o is not None else m.species.replace("_", " ")
        #: Seen only by the expanded lanes: the reference catalogue could not name it.
        self.expansion_only: bool = bool(getattr(m, "expansion_only", False)) or (
            o is not None and getattr(o, "incremental_gain", "baseline") == "expansion")
        #: The level percentile: rank among reference adults who carry it, on
        #: the same distribution as the deviation below, so the two agree.
        self.percentile: float | None = (
            o.level_percentile if (o is not None and getattr(o, "level_percentile", None) is not None)
            else getattr(m, "level_percentile", m.percentile))
        #: Against the typical carrier in the reference (the organism's own lane and cohort).
        self.deviation: float | None = getattr(o, "deviation_percent", None) if o is not None else None
        self.prevalence: float | None = (
            o.prevalence if (o is not None and o.prevalence is not None) else m.cohort_prevalence
        )
        probe = o if o is not None else Organism(species=m.species, percent=self.percent, genus=m.species.split("_")[0])
        self.cls: str = org.verdict(probe).cls
        # Link to the organism's own card only where one exists; the card
        # set is the attention list, and a dead link is worse than none.
        key = o.species if o is not None else m.species
        self.dest: str | None = organism_dest(key, detail=True) if key in carded else None


def _sample_percent_cell(d: GroupDriver) -> str:
    """The member's share of the composition, or where it is counted when it has none."""
    if not d.in_primary:
        within = f"<i>{d.counted_within}</i>" if d.counted_within else "a relative"
        return f"<font size='5.8' color='{_hex(INK_FAINT)}'>within {within}</font>"
    return f"{d.sample_percent:.3f}%" if d.sample_percent < 1 else f"{d.sample_percent:.2f}%"


def group_drivers(g: Any, *, inv: Any | None, carded: frozenset[str]) -> tuple[list[GroupDriver], list[GroupDriver]]:
    """(detected members largest first, undetected members most-common first)."""
    total = sum(m.percent for m in g.members if m.detected and m.percent > 0)
    all_ = [GroupDriver(m, total=total, inv=inv, carded=carded) for m in g.members]
    present = sorted((d for d in all_ if d.detected), key=lambda d: -d.percent)
    absent = sorted((d for d in all_ if not d.detected), key=lambda d: -(d.prevalence or 0.0))
    return present, absent


def deviation_text(dev: float | None) -> str:
    """``+420%`` / ``\u221290%`` against the typical carrier; within 5% reads as typical."""
    if dev is None:
        return ""
    if abs(dev) < 5:
        return "\u2248 typical"
    if dev >= 900:
        return f"\u00d7{(dev / 100.0 + 1.0):,.0f}"
    return f"{'+' if dev > 0 else chr(0x2212)}{abs(dev):,.0f}%"


def _class_dot(cls: str) -> str:
    from openbiota import organisms as org
    colour = {org.BENEFICIAL: GREEN, org.OPPORTUNIST: CORAL, org.CONDITIONAL: AMBER, org.UNKNOWN: SLATE}.get(cls, SLATE)
    return f"<font color='{_hex(colour)}' size='7'>\u25cf</font>"


def _driver_name(d: GroupDriver) -> str:
    name = f"<i>{d.display}</i>"
    return f'<a href="#{d.dest}" color="{_hex(ACCENT_DARK)}">{name}</a>' if d.dest else name


def group_driver_line(present: list[GroupDriver], absent: list[GroupDriver], *, limit: int = 4) -> str:
    """One readable line for the glance row: who carries the reading."""
    from openbiota.pdfreport import INK

    if not present:
        common = [d for d in absent if (d.prevalence or 0) >= 0.5][:3]
        if common:
            names = ", ".join(f"<i>{d.display}</i>" for d in common)
            return (f"<font size='6.8' color='{_hex(CORAL)}'><b>none detected</b></font> "
                    f"<font size='6.8' color='{_hex(INK)}'>\u2014 usually present: {names}</font>")
        return f"<font size='6.8' color='{_hex(INK_FAINT)}'>none of this group's species were detected</font>"
    parts = [f"{_class_dot(d.cls)} {_driver_name(d)} <b>{d.share:.0%}</b>" for d in present[:limit]]
    rest = len(present) - min(len(present), limit)
    tail = f" <font color='{_hex(INK_FAINT)}'>+{rest} more</font>" if rest > 0 else ""
    missing = ""
    if absent:
        common = [d for d in absent if (d.prevalence or 0) >= 0.5]
        if common:
            missing = (f" &nbsp;<font color='{_hex(CORAL)}'>not detected:</font> "
                       + ", ".join(f"<i>{d.display}</i>" for d in common[:3])
                       + (f" <font color='{_hex(INK_FAINT)}'>+{len(common) - 3}</font>" if len(common) > 3 else ""))
    sep = " \u00b7 "
    return f"<font size='6.8' color='{_hex(INK)}'>{sep.join(parts)}{tail}{missing}</font>"


def _group_card(story: list[Any], st: dict[str, ParagraphStyle], g: Any,
                *, lead: list[Any] | None = None, cards: Any | None = None,
                coherence: Any | None = None, inv: Any | None = None,
                carded: frozenset[str] = frozenset()) -> None:
    from openbiota.pdfatlas import evidence_cards

    status = status_for(percentile=g.level_percentile, higher_means=g.group.effective_direction,
                        detected=g.detected or g.level_percentile is not None)
    if not g.detected and g.level_percentile is not None:
        # Absent, but placed: most of the reference lacks it too, or does not.
        status = Status(f"not detected · {status.label}", status.colour, status.background, status.note, status.level)
    present, absent = group_drivers(g, inv=inv, carded=carded)
    block: list[Any] = [
        *(lead or []),
        Anchor(group_dest(g.group.name, detail=True), above=6 * mm),
        back_link_line(group_dest(g.group.name), "at a glance"),
    ]
    head = Table(
        [[
            Paragraph(f"<b>{g.group.label}</b>{_arrow(g.group.effective_direction)}<br/>"
                      f"<font color='{_hex(INK_FAINT)}' size='6.5'>{g.percent:.2f}% of classified reads"
                      + (f" · <font color='{_hex(status.colour)}'><b>{deviation_text(getattr(g, 'deviation_percent', None))}</b></font> vs typical"
                         if getattr(g, "deviation_percent", None) is not None else "")
                      + f" · {g.n_detected} of {g.n_resolvable} species detected"
                      + (f" · {len(g.group.unresolved)} in the literature not nameable here" if g.group.unresolved else "")
                      + "</font>", st["cell"]),
            PercentileBar(width=CONTENT_WIDTH * 0.30, percentile=g.level_percentile, higher_means=g.group.effective_direction,
                          marker_colour=status.colour, height=6.4 * mm),
            Paragraph(f"<font color='{_hex(status.colour)}'>{_ordinal(g.level_percentile)}</font>", st["pct"]),
            StatusChip(status, width=CONTENT_WIDTH * 0.20 - 4 * mm, font_size=6.4),
        ]],
        colWidths=[CONTENT_WIDTH * 0.43, CONTENT_WIDTH * 0.31, CONTENT_WIDTH * 0.06, CONTENT_WIDTH * 0.20],
        hAlign="LEFT",
    )
    head.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 2), ("LEFTPADDING", (3, 0), (3, 0), 4)]))
    block.append(head)
    block.append(Paragraph(" ".join(g.group.summary.split()), st["small"]))
    if g.group.caveats:
        block.append(Paragraph(" ".join(g.group.caveats[0].split()), st["fine"]))
    # The organisms behind the number: every member, in the status colour,
    # the same table the metabolite cards carry. A group is a sum and the
    # sum has named parts; this is where they are named.
    _group_drivers_block(block, st, g, present, absent, status)
    # The same function measured by gene search: the identical sentence the
    # panel card prints, so the two readings explain each other in both places.
    if coherence is not None:
        block.append(Spacer(1, 1 * mm))
        block.append(_group_coherence_note(coherence, st))
    # Parts two, three and four of the shared metric card. Every metric in this
    # report carries the same four headings in the same order; see
    # `openbiota.metriccard` before changing what a card contains.
    _group_meaning_block(block, st, g, status, present)
    story.append(KeepTogether(block))
    if cards:
        evidence_cards(story, st, cards)
    story.append(Spacer(1, 1 * mm))
    story.append(Rule(CONTENT_WIDTH))
    story.append(Spacer(1, 2 * mm))


def _group_meaning_block(block: list[Any], st: dict[str, ParagraphStyle], g: Any,
                         status: Status, present: list[GroupDriver]) -> None:
    """The last three sections of the shared metric card, for a group.

    A group reading is a sum over species, so what would move it is which
    species are there - and for a low group, the ones that are not. That
    makes the improvement section answerable in a way it often is not
    elsewhere: the named absentees are the gap.
    """
    from openbiota import metriccard as MC
    from openbiota.pdfreport import INK_SOFT

    # 2 · What your reading means
    block.append(Paragraph(
        f"<b>{MC.CARD_SECTIONS[1]}</b> &nbsp;<font color='{_hex(INK_FAINT)}' size='7'>"
        f"{status.label}</font>", st["h3"]))
    detected = f"{g.n_detected} of the {g.n_resolvable} species this group is defined by"
    block.append(Paragraph(
        f"This group is {g.percent:.2f}% of your classified reads, counted over {detected}. "
        + (f"That places you at the {_ordinal(g.level_percentile)} percentile among reference adults who have the group. "
           if g.level_percentile is not None else "There is no reference position for it. ")
        + "A group reading is the sum of its species, so it moves when they do.",
        st["body_ink"]))

    # 3 · How can I improve this
    MC.improvement_block(
        block, st,
        higher_means=g.group.effective_direction,
        organisms=[
            {"organism": d.display, "sample_percent": d.sample_percent if d.in_primary else None}
            for d in present[:4]
        ],
        nothing_established=(
            "What is established for each species is on its own page, and the members "
            "not detected are listed above: for a low reading those absentees are the "
            "gap, and they are the specific thing any change would have to supply."
        ),
    )

    # 4 · The research behind it
    block.append(Paragraph(f"<b>{MC.CARD_SECTIONS[3]}</b>", st["h3"]))
    block.append(Paragraph(" ".join(g.group.summary.split()), st["body"]))
    for caveat in (g.group.caveats or ())[:2]:
        block.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>{' '.join(caveat.split())}</font>",
            st["small"]))


def _group_drivers_block(block: list[Any], st: dict[str, ParagraphStyle], g: Any,
                         present: list[GroupDriver], absent: list[GroupDriver], status: Status) -> None:
    """'What is driving this reading': every member species of the group.

    Detected members largest first, each with its share of the group, its own
    abundance, its class and its position among reference carriers, linked
    to its own card where the report has one. Then the members that were not
    detected, most-common first - for a low group, those *are* the cause.
    """
    from openbiota import organisms as org
    from openbiota.pdfreport import INK, INK_SOFT

    colour = status.colour if status.colour != SLATE else INK_SOFT
    n_total = len(present) + len(absent)
    block.append(Paragraph(
        f"<b>What is driving this reading</b> &nbsp;<font color='{_hex(INK_FAINT)}' size='7'>"
        f"{len(present)} of {n_total} species detected; the reading is their sum</font>",
        st["h3"]))

    if present:
        rows: list[list[Any]] = [[
            Paragraph("ORGANISM", st["label"]), Paragraph("CLASS", st["label"]),
            Paragraph("SHARE OF GROUP", st["label"]), Paragraph("IN YOUR SAMPLE", st["label"]),
            Paragraph("YOUR LEVEL", st["label"]),
        ]]
        for d in present:
            level = (f"<b>{_ordinal(d.percentile)}</b> percentile"
                     if d.percentile is not None else f"<font color='{_hex(INK_FAINT)}'>no reference</font>")
            if d.deviation is not None:
                level = f"<b>{deviation_text(d.deviation)}</b> <font size='5.6' color='{_hex(INK_FAINT)}'>vs typical carrier</font><br/>" + level
            if d.prevalence is not None:
                level += f"<br/><font size='6' color='{_hex(INK_FAINT)}'>carried by {d.prevalence:.0%} of adults</font>"
            name = f"<b>{_driver_name(d)}</b>"
            if d.dest:
                name += f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>full page \u2192</font>"
            rows.append([
                Paragraph(name, st["cell"]),
                Paragraph(f"{_class_dot(d.cls)} <font size='6.8'>{org.CLASS_LABEL[d.cls]}</font>", st["cell"]),
                Paragraph(f"<font color='{_hex(colour)}'><b>{d.share:.0%}</b></font>", st["cell"]),
                Paragraph(_sample_percent_cell(d), st["cell"]),
                Paragraph(level, st["cell"]),
            ])
        widths = [CONTENT_WIDTH * w for w in (0.34, 0.15, 0.14, 0.13, 0.24)]
        t = Table(rows, colWidths=widths, hAlign="LEFT")
        t.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE), ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
            ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
            ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ]))
        block.append(t)

    # The root cause, in one sentence, in the status colour.
    lead = present[0] if present else None
    common_absent = [d for d in absent if (d.prevalence or 0) >= 0.5]
    high = g.level_percentile is not None and g.level_percentile >= 75
    low = g.level_percentile is not None and g.level_percentile <= 25
    if lead is not None and high:
        why = (f"This reading is carried by {_driver_name(lead)}, which alone is <b>{lead.share:.0%}</b> of the "
               f"group's reads"
               + (f" and sits at the <b>{_ordinal(lead.percentile)}</b> percentile among adults who carry it"
                  if lead.percentile is not None else "") + ".")
        if len(present) > 1 and present[1].share >= 0.2:
            why += f" {_driver_name(present[1])} contributes another <b>{present[1].share:.0%}</b>."
        why += " Changing this reading means changing them."
    elif low and common_absent:
        names = ", ".join(f"{_driver_name(d)} <font color='{_hex(INK_FAINT)}' size='6.4'>({d.prevalence:.0%})</font>"
                          for d in common_absent[:4])
        why = (f"This reading is low because <b>{len(common_absent)}</b> species most reference adults carry "
               f"were not detected: {names}"
               + (f", and {len(common_absent) - 4} more" if len(common_absent) > 4 else "") + "."
               + (f" What is present, led by {_driver_name(lead)}, is not enough to make up the sum." if lead else ""))
    elif low and lead is not None:
        why = (f"Every common member is present; the reading is low because they are scarce. {_driver_name(lead)} "
               + (f"is the largest at {lead.sample_percent:.3f}% of your community" if lead.in_primary else
                  "is the largest, counted within a relative's share of your community")
               + (f", the <b>{_ordinal(lead.percentile)}</b> percentile among carriers" if lead.percentile is not None else "")
               + ".")
    elif lead is not None:
        why = (f"Led by {_driver_name(lead)} at <b>{lead.share:.0%}</b> of the group"
               + (f", the {_ordinal(lead.percentile)} percentile among carriers" if lead.percentile is not None else "")
               + ". The group as a whole sits within the reference range.")
    else:
        why = "None of this group's species were detected in your sample."
    block.append(Paragraph(f"<font color='{_hex(colour)}'><b>Why.</b></font> <font color='{_hex(INK)}'>{why}</font>", st["small"]))

    if absent:
        names = ", ".join(
            f"<i>{d.display}</i>"
            + (f" <font color='{_hex(INK_FAINT)}' size='6'>({d.prevalence:.0%})</font>" if d.prevalence is not None else "")
            for d in absent
        )
        block.append(Paragraph(
            f"<font size='6.8'><font color='{_hex(CORAL if common_absent else INK_FAINT)}'><b>Not detected "
            f"({len(absent)} of {n_total}):</b></font> {names}. "
            f"<font color='{_hex(INK_FAINT)}'>Percentages are how many reference adults carry each.</font></font>",
            st["small"]))
    block.append(Spacer(1, 1 * mm))


def _group_coherence_note(coh: Any, st: dict[str, ParagraphStyle]) -> Any:
    """The group-side twin of the panel card's note: same facts, this side's wording."""
    from openbiota.pdfreport import AMBER, AMBER_BG, SLATE, SLATE_BG, Table, TableStyle, _hex

    agree = coh.agrees
    colour, bg = (SLATE, SLATE_BG) if agree else (AMBER, AMBER_BG)
    lead = "Two measurements agree" if agree else "Two measurements of the same function disagree"
    cell = [
        Paragraph(f"<font color='{_hex(colour)}'><b>{lead}</b></font>", st["small"]),
        Paragraph(coh.explanation(side="carriers"), st["small"]),
    ]
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


def _species_table(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    community: Any,
    *,
    section: int,
    cohort_note: str,
    inv: Any = None,
    sub: int | None = None,
) -> None:
    """Every organism detected, merged across every catalogue that ran.

    This used to list one catalogue's 92 species. The merged inventory
    typically carries around 187, because the newer catalogue is built from
    roughly a million genomes and sees organisms the older one has no
    reference for. Everything detected is listed here; the percentile column
    is the one thing that cannot cross catalogues, and says so per row.
    """
    if inv is None or not len(inv):
        _species_table_single(
            story, st, community, section=section, cohort_note=cohort_note, sub=sub)
        return

    named, unnamed = inv.named, inv.unnamed
    c = inv.counts()
    story.extend(section_heading(
        section, "Every organism detected", st["h1"], sub=sub,
        opens_section=sub is None))
    n_lanes = len(getattr(inv, "lanes", ()) or ())
    story.append(Paragraph(
        f"All <b>{c['organisms']} organisms</b> found in this sample across "
        f"{c['genera']} genera, largest first. Detection pools <b>{n_lanes} independent methods and "
        "reference catalogues</b>, so this is the widest view the sequencing supports — a species one "
        "catalogue has no reference genome for is still reported if another can name it. "
        f"Counted separately: <b>{c.get('category_named_species', c['named'])}</b> named species, "
        f"<b>{c.get('category_unnamed_species_cluster', c['unnamed'])}</b> unnamed species-level clusters"
        + (f", <b>{c['category_unresolved_complex']}</b> unresolved complexes" if c.get("category_unresolved_complex") else "")
        + (f", <b>{c['category_higher_rank']}</b> placed above species" if c.get("category_higher_rank") else "")
        + f". <b>{c['strain_resolved']}</b> were characterised further, to the specific population "
        "you carry rather than the species alone.",
        st["body"]))
    story.append(Paragraph(
        f"<b>Percentile</b> is your abundance against {cohort_note}, treating 'not carried' as "
        "its own category; <b>carried by</b> is the share of that reference group carrying the "
        "organism at all — one carried by 8% of people is unusual to have, whatever the amount. "
        f"{c['organisms'] - c['rankable']} of these organisms have no percentile: the reference "
        "group predates their description, so they can be measured in you but not yet placed "
        "against other people. That is a gap in the reference data, and the detection stands "
        "on its own.",
        st["body"]))
    lane_pops = {k: v for k, v in (getattr(inv, "percentile_sources", None) or {}).items() if k != "scoring cohort"}
    if lane_pops:
        from openbiota.expansion import cohorts as _cohorts

        described = []
        for key, n in lane_pops.items():
            lane = key.split(" cohort", 1)[0]
            lc = _cohorts.load_all().get(lane)
            described.append(f"<b>{n}</b> organisms marked ‡ are ranked in the {key}: "
                             + (lc.description if lc is not None else "a population profiled by that lane alone") + ".")
        story.append(Paragraph(
            " ".join(described) + " A rank from one population is never compared with a rank from another; "
            "each says only where you sit among the people it contains.",
            st["body"]))
    story.append(Spacer(1, 2 * mm))

    strained = c["strain_resolved"] > 0
    widths = [CONTENT_WIDTH * w for w in (
        (0.25, 0.11, 0.15, 0.09, 0.18, 0.22) if strained else (0.30, 0.12, 0.19, 0.10, 0.29)
    )]
    header = [
        Paragraph("ORGANISM", st["label"]), Paragraph("ABUNDANCE", st["label"]),
        Paragraph("PERCENTILE", st["label"]), Paragraph("CARRIED BY", st["label"]),
    ]
    if strained:
        header.append(Paragraph("STRAIN DETAIL", st["label"]))
    header.append(Paragraph("GROUPS", st["label"]))
    rows: list[list[Any]] = [header]
    for o in named:
        rows.append(_organism_row(o, st, strained=strained))
    story.append(_plain_table(rows, widths, zebra=True))
    story.append(Spacer(1, 2 * mm))

    if unnamed:
        story.append(Paragraph(
            f"<b>A further {len(unnamed)} organisms have no name yet.</b> They are real and "
            "measured — assembled from metagenomes like yours, given a stable identifier, but "
            "never grown in a laboratory, so nobody has described or named them. They are "
            "listed separately because a name is what links an organism to published research: "
            "without one there is nothing to look up. Between them they account for "
            f"{sum(o.percent for o in unnamed):.1f}% of your classified reads, so this is not a "
            "rounding error — it is the part of your gut that science has not caught up with.",
            st["body"]))
        story.append(Spacer(1, 1.5 * mm))
        uw = [CONTENT_WIDTH * w for w in (0.30, 0.14, 0.26, 0.30)]
        urows: list[list[Any]] = [[
            Paragraph("IDENTIFIER", st["label"]), Paragraph("ABUNDANCE", st["label"]),
            Paragraph("NEAREST NAMED GROUP", st["label"]), Paragraph("STRAIN DETAIL", st["label"]),
        ]]
        for o in unnamed:
            fp = o.strain or {}
            ident = o.sgb or next(iter((getattr(o, "native_ids", {}) or {}).values()), "") or o.display
            ident = str(ident).split(":", 1)[-1]
            flags = ""
            urows.append([
                Paragraph(f"<font color='{_hex(INK_FAINT)}'>{ident}</font>{flags}", st["cell"]),
                Paragraph(_pct(o.percent) if (o.in_primary and o.percent) else
                          (f"<font color='{_hex(INK_FAINT)}'>within <i>{o.counted_within}</i></font>"
                           if getattr(o, "counted_within", None) else
                           f"<font color='{_hex(INK_FAINT)}'>detected</font>"), st["cell"]),
                Paragraph(f"<i>{o.genus.replace('_', ' ')}</i>" if o.genus else "—", st["fine"]),
                Paragraph(
                    f"<font color='{_hex(ACCENT_DARK)}'>resolved</font> "
                    f"<font color='{_hex(INK_FAINT)}'>{fp.get('markers_resolved')} markers</font>"
                    if fp else f"<font color='{_hex(INK_FAINT)}'>—</font>", st["fine"]),
            ])
        story.append(_plain_table(urows, uw, zebra=True))
        story.append(Spacer(1, 2 * mm))

    _detection_coverage(story, st, inv)
    for note in community.notes:
        story.append(Paragraph(note, st["fine"]))
    story.append(Paragraph(
        "Abundance is the share of your classified reads. Where an organism has been "
        "reclassified, the current name is shown with its previous name beneath, so a name you "
        "recognise from older results is still findable. Orange 'carried by' figures mark "
        "organisms found in fewer than 10% of the reference group.",
        st["fine"]))


def _pct(v: float) -> str:
    return f"{v:.3f}%" if v < 1 else f"{v:.2f}%"


#: Short words for the identifier line under an organism: lane -> label.
_ID_LABEL = {"jan26": "SGB", "extended": "SGB", "globdb": "GlobDB", "motus": "mOTU", "kraken": "UHGG",
             "rescue": "panel", "singlem": "SingleM", "genome": "GTDB"}
_LANE_SHORT = {"scoring": "MP3", "extended": "Jun23", "jan26": "Jan26", "genome": "GTDB sketch",
               "globdb": "GlobDB sketch", "motus": "mOTUs", "kraken": "Kraken", "rescue": "Kraken panel", "singlem": "SingleM"}


def _identifiers_line(o: Any) -> str:
    """Stable native identifiers and detection support, one faint line.

    Spec 0.8.4 §6: every accepted organism shows its native identifiers and
    which methods saw it. Identifiers are the lanes' own (SGB number,
    GlobDB/GTDB accession, mOTU id); the lane list is the provenance.
    """
    parts: list[str] = []
    seen_ids: set[str] = set()
    ids = dict(getattr(o, "native_ids", {}) or {})
    if o.sgb and "jan26" not in ids and "extended" not in ids:
        ids = {"jan26": o.sgb, **ids}
    for lane in ("jan26", "extended", "globdb", "genome", "motus", "kraken", "singlem"):
        raw = ids.get(lane)
        if not raw:
            continue
        val = str(raw).split(":", 1)[-1].split(";")[0]
        if val.startswith(("s__", "g__")):
            continue  # a lineage name, not an identifier; the organism's name already says it
        if lane in ("jan26", "extended") and val.startswith("SGB"):
            val = val[3:]
        if lane == "motus" and val.startswith("mOTUv"):
            val = val.split("_", 1)[-1]
        base = val.split(".")[0] if val.startswith(("GCA_", "GCF_")) else val
        if base in seen_ids or not val:
            continue
        seen_ids.add(base)
        parts.append(f"{_ID_LABEL.get(lane, lane)} {val}")
        if len(parts) >= 3:
            break
    lanes = [_LANE_SHORT.get(x, x) for x in getattr(o, "lanes", ()) or ()]
    if lanes:
        # The methods are named in full in "How this list was detected";
        # here the count keeps the row to one faint line.
        parts.append(f"{len(lanes)} method" + ("s" if len(lanes) != 1 else "") + f": {', '.join(lanes)}"
                     if len(lanes) <= 2 else f"{len(lanes)} methods")
    return " · ".join(parts)


def _organism_row(o: Any, st: dict[str, ParagraphStyle], *, strained: bool) -> list[Any]:
    """One organism's row in the merged inventory."""
    pct_text = f"<font color='{_hex(INK_FAINT)}' size='6'>not in reference set</font>"
    if o.percentile is not None:
        lvl = status_for(percentile=o.percentile, higher_means="unclear")
        src = getattr(o, "percentile_source", None) or "scoring cohort"
        mark = "" if src == "scoring cohort" else "‡"
        pct_text = (f"{_ordinal(o.percentile)}{mark} <font color='{_hex(INK_FAINT)}' "
                    f"size='6'>{lvl.label}</font>")
    prev = "—" if o.prevalence is None else f"{o.prevalence:.0%}"
    if o.rare and o.prevalence is not None:
        prev = f"<font color='{_hex(ORANGE)}'>{prev}</font>"
    name = f"<i>{o.display}</i>"
    if o.trace:
        name += f" <font color='{_hex(INK_FAINT)}' size='6'>trace</font>"
    if getattr(o, "status", "") == "ambiguous":
        name += f" <font color='{_hex(INK_SOFT)}' size='6'>sibling species, one population</font>"
    if o.formerly:
        name += (f"<br/><font color='{_hex(INK_FAINT)}' size='5.6'>formerly "
                 f"{o.formerly.replace('_', ' ')}</font>")
    ids_line = _identifiers_line(o)
    if ids_line:
        name += f"<br/><font color='{_hex(INK_FAINT)}' size='5.4'>{ids_line}</font>"
    row: list[Any] = [
        Paragraph(name, st["cell"]),
        Paragraph(_pct(o.percent), st["cell"]),
        Paragraph(pct_text, st["cell"]),
        Paragraph(prev, st["cell"]),
    ]
    if strained:
        fp = o.strain or {}
        if fp and fp.get("nearest_reference"):
            near = str(fp["nearest_reference"]).split("|")[0].replace(".fna.gz", "").replace(".fna", "")
            label = "resolved, placed" if fp.get("markers_resolved") else "placed"
            row.append(Paragraph(
                f"<font color='{_hex(ACCENT_DARK)}'><b>{label}</b></font> "
                f"<font color='{_hex(INK_FAINT)}'>nearest {near}, "
                f"{fp.get('markers_resolved') or fp.get('n_markers')} markers</font>", st["fine"]))
        else:
            row.append(Paragraph(
                f"<font color='{_hex(ACCENT_DARK)}'><b>resolved</b></font> "
                f"<font color='{_hex(INK_FAINT)}'>{fp.get('markers_resolved')} markers, "
                f"{int(fp.get('callable_bases') or 0) / 1000:.0f}k bases</font>"
                if fp else f"<font color='{_hex(INK_FAINT)}'>—</font>", st["fine"]))
    row.append(Paragraph(", ".join(g.replace("_", " ") for g in o.groups), st["fine"]))
    return row


def _species_table_single(
    story: list[Any], st: dict[str, ParagraphStyle], community: Any,
    *, section: int, cohort_note: str, sub: int | None = None,
) -> None:
    """Fallback listing when only one catalogue produced a result."""
    story.extend(section_heading(
        section, "Every organism detected", st["h1"], sub=sub,
        opens_section=sub is None))
    story.append(Paragraph(
        f"All {community.n_species_detected} species named in this sample across "
        f"{community.n_genera_detected} genera, largest first. <b>Percentile</b> is your "
        f"abundance against {cohort_note}, treating 'not carried' as its own category; "
        "<b>carried by</b> is the share of that reference group in which the species is found "
        "at all. Species without a percentile are not in the reference catalogue and cannot "
        "be placed.",
        st["body"]))
    story.append(Spacer(1, 2 * mm))
    widths = [CONTENT_WIDTH * w for w in (0.34, 0.10, 0.19, 0.10, 0.27)]
    rows: list[list[Any]] = [[
        Paragraph("SPECIES", st["label"]), Paragraph("ABUNDANCE", st["label"]),
        Paragraph("PERCENTILE", st["label"]), Paragraph("CARRIED BY", st["label"]),
        Paragraph("GROUPS", st["label"]),
    ]]
    for s in community.species:
        pct_text = "—"
        if s.percentile is not None:
            lvl_status = status_for(percentile=s.percentile, higher_means="unclear")
            pct_text = (f"{_ordinal(s.percentile)} <font color='{_hex(INK_FAINT)}' "
                        f"size='6'>{lvl_status.label}</font>")
        prev = "—" if s.cohort_prevalence is None else f"{s.cohort_prevalence:.0%}"
        if s.rare_in_cohort:
            prev = f"<font color='{_hex(ORANGE)}'>{prev}</font>"
        rows.append([
            Paragraph(f"<i>{s.species.replace('_', ' ')}</i>"
                      + (f" <font color='{_hex(INK_FAINT)}' size='6'>trace</font>"
                         if s.trace else ""), st["cell"]),
            Paragraph(_pct(s.percent), st["cell"]),
            Paragraph(pct_text, st["cell"]),
            Paragraph(prev, st["cell"]),
            Paragraph(", ".join(g.replace("_", " ") for g in s.groups), st["fine"]),
        ])
    story.append(_plain_table(rows, widths, zebra=True))
    story.append(Spacer(1, 2 * mm))
    for note in community.notes:
        story.append(Paragraph(note, st["fine"]))
    story.append(Paragraph(
        "Abundance is the share of your classified reads. Orange 'carried by' figures mark "
        "species found in fewer than 10% of the reference group.",
        st["fine"]))


# --------------------------------------------------------------------------- #
# sequencing quality + reference coverage
# --------------------------------------------------------------------------- #


GATE_COLOUR = {"pass": GREEN, "warn": AMBER, "fail": CORAL, "unknown": SLATE, "not_assessable": SLATE}


def sequencing_quality_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    quality: dict[str, Any] | None,
    extended: dict[str, Any] | None,
) -> None:
    story.extend(section_heading(section, "Sequencing quality and reference coverage", st["h1"]))
    story.append(
        Paragraph(
            "Whether the sample was sequenced well enough for the readings above to mean what they "
            "say. Each gate has a prespecified threshold; a warning means the reading is shown but "
            "should be read with that in mind, and a failure means the affected section says so.",
            st["body"],
        )
    )
    story.append(Spacer(1, 2 * mm))
    if not quality:
        story.append(Paragraph("Quality gates were not assembled for this run.", st["body"]))
    else:
        gates = quality.get("gates") or {}
        fastp = quality.get("fastp")
        overall = gates.get("overall", "unknown")
        colour = GATE_COLOUR.get(overall, SLATE)
        tiles = [
            ("OVERALL", overall.upper(), "across all gates"),
            ("USABLE NON-HOST PAIRS",
             f"{gates.get('usable_nonhost_pairs'):,}" if gates.get("usable_nonhost_pairs") is not None else "—",
             f"minimum {gates.get('minimum_usable_pairs', 500000):,}"),
        ]
        if fastp:
            tiles += [
                ("BASES AT Q30", f"{fastp['q30_rate']:.1%}", "≤ 0.1% error probability"),
                ("DUPLICATE RATE", f"{fastp['duplication_rate']:.2%}", "whole-file, fastp"),
                ("READ LENGTH", f"{fastp['read1_mean_length']:.0f} bp", fastp.get("sequencing", "")),
            ]
        story.append(_tiles(tiles, accent=colour))
        story.append(Spacer(1, 3 * mm))
        rows: list[list[Any]] = [[Paragraph("GATE", st["label"]), Paragraph("VALUE", st["label"]),
                                  Paragraph("", st["label"]), Paragraph("NOTE", st["label"])]]
        for g in gates.get("gates", []):
            c = GATE_COLOUR.get(g["status"], SLATE)
            rows.append([
                Paragraph(f"<b>{g['label']}</b>", st["cell"]),
                Paragraph(g["display"], st["cell"]),
                Table([[Dot(c, size=2.4 * mm), Paragraph(g["status"], st["cell"])]], colWidths=[4 * mm, 14 * mm],
                      hAlign="LEFT", style=TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                                       ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0)])),
                Paragraph(g["note"], st["small"]),
            ])
        story.append(_plain_table(rows, [CONTENT_WIDTH * 0.26, CONTENT_WIDTH * 0.22, CONTENT_WIDTH * 0.10, CONTENT_WIDTH * 0.42]))
        prov = gates.get("provenance") or {}
        if prov.get("header_style") == "illumina":
            story.append(Spacer(1, 1.5 * mm))
            story.append(Paragraph(
                f"Sequenced on instrument {prov.get('instrument')}, run {prov.get('run')}, flowcell "
                f"{'+'.join(prov.get('flowcells', []))}, lane(s) {', '.join(str(x) for x in prov.get('lanes', []))}. "
                "Samples from different flowcells carry a batch difference that can exceed the biology; "
                "compare across samples only with that in mind.", st["fine"]))
        if fastp:
            story.append(Paragraph(
                f"Read preprocessing: fastp {fastp['tool'].split()[-1]} ({fastp['mode']} mode). "
                f"{fastp['insert_note'][0].upper() + fastp['insert_note'][1:]}.", st["fine"]))

    story.append(Spacer(1, 4 * mm))
    story.extend(section_heading(section, "How much of your community has a name", st["h2"], sub=1))
    if not extended:
        story.append(Paragraph(
            "Only one reference catalogue was searched for this sample, covering about 13,500 "
            "species \u2014 the one the reference population and the disease patterns are built "
            "on. The organism list is complete with respect to that catalogue.", st["body"]))
        return
    story.append(Paragraph(
        "Your reads were searched against every reference catalogue available, and the organism "
        f"list in section {SECTIONS['catalogue']} is the merged result. The largest of them holds "
        "about 26,000 organisms "
        "assembled from roughly a million genomes, including some 4,900 that have never been grown "
        "in a laboratory and so have no formal name. The figures below say how much of your "
        "community each search could account for. Percentiles are the one thing that cannot be "
        "pooled: a position is only meaningful against a reference population measured the same "
        "way, and that population exists for one catalogue.", st["body"]))
    story.append(Spacer(1, 2 * mm))
    kc = extended.get("kingdom_counts", {})
    story.append(_tiles([
        ("GENOME BINS (SGBs)", f"{extended['n_sgbs_detected']}", f"{extended['n_species_detected']} species · {extended['n_genera_detected']} genera"),
        ("WITHOUT A NAME", f"{extended['n_unnamed_sgbs']}", f"{extended['unnamed_sgb_percent']:.1f}% of classified reads"),
        ("ARCHAEA", f"{kc.get('Archaea', 0)}", "methanogens and relatives"),
        ("EUKARYOTES", f"{kc.get('Eukaryota', 0)}", "fungi and protists in the catalogue"),
        ("UNCLASSIFIED", f"{extended['unclassified_percent']:.1f}%", "reads matching nothing on record"),
    ]))
    bridge = extended.get("bridge_to_scoring_lane")
    if bridge:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<b>How much the catalogues agree.</b> {bridge['named_by_both']} organisms were found "
            f"by all of them. {bridge['named_only_by_metaphlan4']} were found only by the larger "
            f"catalogue \u2014 recently described organisms, renamed genera, and organisms with no "
            f"name \u2014 and {bridge['named_only_by_metaphlan3']} only by the smaller one, mostly "
            "names that have since been split or merged. All of them appear in your organism list: "
            "pooling the catalogues is what makes that list complete, and a detection found by one "
            "search and not another is a difference in reference coverage, not a contradiction.",
            st["small"]))
        if bridge.get("metaphlan4_only_examples"):
            ex = ", ".join(f"<i>{n.replace('_', ' ')}</i>" for n in bridge["metaphlan4_only_examples"][:10])
            story.append(Paragraph(
                f"Largest found only by the wider search: {ex}.", st["fine"]))
    unnamed = [r for r in extended.get("sgbs", []) if r.get("unnamed")][:8]
    if unnamed:
        story.append(Spacer(1, 1.5 * mm))
        rows = [[Paragraph("IDENTIFIER", st["label"]), Paragraph("NEAREST NAMED LINEAGE", st["label"]),
                 Paragraph("ABUNDANCE", st["label"])]]
        for r in unnamed:
            rows.append([Paragraph(r["sgb"], st["cell"]),
                         Paragraph(f"{r['kingdom']} · <i>{r['genus'].replace('_', ' ') or '—'}</i>", st["cell"]),
                         Paragraph(f"{r['percent']:.3f}%", st["cell"])])
        story.append(_plain_table(rows, [CONTENT_WIDTH * 0.24, CONTENT_WIDTH * 0.56, CONTENT_WIDTH * 0.20]))
        story.append(Paragraph(
            "Unnamed organisms are known only from assembled genomes. They are real and "
            "often abundant; what is missing is a cultured isolate and therefore any physiology to "
            "report. Their share is one measure of how much of your community current science has not "
            "yet characterised.", st["fine"]))

def catalogue_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    community: Any | None,
    cohort_note: str,
    inv: Any = None,
) -> None:
    """Every organism detected, as its own section.

    With a merged inventory this is the classified section: composition
    bar, the thirty largest with descriptions, then everything colour-coded.
    Without one (a results file from before the inventory existed) it
    falls back to the single-catalogue table.
    """
    if community is None and inv is None:
        story.extend(section_heading(
            section, "Every organism detected", st["h1"], opens_section=True))
        story.append(Paragraph(
            "No taxonomic profile was produced for this sample, so there is no "
            "organism list.", st["body"]))
        return
    if inv is not None and len(inv):
        from openbiota import organisms as org
        from openbiota.pdforganisms import organism_section

        organism_section(
            story, st, section=section,
            vs=org.verdicts(list(inv.organisms)),
            cohort_note=cohort_note,
            n_strain_resolved=len(inv.strain_resolved),
            unclassified_percent=(inv.unplaced_percent if inv.unplaced_percent is not None else inv.unclassified_percent),
            inv=inv,
        )
        return
    _species_table(
        story, st, community, section=section, cohort_note=cohort_note, inv=inv)


__all__ = [
    "catalogue_page",
    "community_groups_pages",
    "library_pages",
    "sequencing_quality_page",
]


#: Reader-facing names for the detection lanes, and what each one is.
LANE_WORDS: dict[str, tuple[str, str]] = {
    "scoring": ("MetaPhlAn 3 markers", "the calibrated catalogue every percentile is measured on"),
    "extended": ("MetaPhlAn 4 markers (Jun23)", "36,822 species-level genome bins; the previous detection catalogue"),
    "jan26": ("MetaPhlAn 4.2 markers (Jan26)", "72,000 species-level genome bins; the composition shown"),
    "genome": ("Whole-genome sketch, GTDB R232", "199,923 species clusters"),
    "globdb": ("Whole-genome sketch, GlobDB r232", "346,233 representatives from 26 catalogues"),
    "motus": ("Universal marker genes, mOTUs 4.1", "about 124,300 species-level units"),
    "kraken": ("Read classification, UHGG v2.0.2", "4,744 gut species; Bracken abundance"),
    "rescue": ("Read classification, gut rescue panel", "reconciliation targets with near neighbours and decoys; built locally"),
    "singlem": ("Marker windows, SingleM / GlobDB", "places lineages no genome catalogue holds"),
}


def _detection_coverage(story: list[Any], st: dict[str, ParagraphStyle], inv: Any) -> None:
    """How this list was detected: every lane that contributed, and the counts.

    Spec 0.8.4 §6 asks for a compact coverage summary beside the inventory.
    Each lane is one row: what it is, how many of the organisms above it
    saw, how many it alone saw. The releases are the ones the run used.
    """
    lanes = list(getattr(inv, "lanes", ()) or ())
    if not lanes:
        return
    organisms = list(inv.organisms)
    block: list[Any] = [Spacer(1, 2 * mm), Paragraph("<b>How this list was detected</b>", st["h3"]), Paragraph(
        "Each method reads the same DNA a different way. An organism seen by two independent methods is "
        "supported; one seen by a single method at a marginal level is provisional. Composition shares come "
        "from one method only and are never added across methods.",
        st["small"])]
    rows: list[list[Any]] = [[
        Paragraph("METHOD", st["label"]), Paragraph("REFERENCE", st["label"]),
        Paragraph("ORGANISMS SEEN", st["label"]), Paragraph("SEEN BY THIS METHOD ONLY", st["label"]),
    ]]
    for lane in lanes:
        words, ref = LANE_WORDS.get(lane, (lane, ""))
        seen = sum(1 for o in organisms if lane in o.lanes)
        only = sum(1 for o in organisms if o.lanes == (lane,))
        rows.append([
            Paragraph(f"<b>{words}</b>", st["cell"]), Paragraph(ref, st["fine"]),
            Paragraph(str(seen), st["cell"]), Paragraph(str(only) if only else "—", st["cell"]),
        ])
    widths = [CONTENT_WIDTH * w for w in (0.32, 0.36, 0.14, 0.18)]
    block.append(_plain_table(rows, widths, zebra=True))
    story.append(KeepTogether(block))  # heading and its table on one page
    story.append(Spacer(1, 2 * mm))

