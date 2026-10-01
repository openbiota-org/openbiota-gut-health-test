"""PART B pages of the report: species composition and disease-pattern resemblance.

Kept apart from `pdfreport` so the two halves can evolve independently, but
drawn with the same flowables, palette and — above all — the same scale, so a
"high" here means what "high" meant on the metabolite pages.

Rendering rules that come from the spec rather than taste:

* Profile pages show the module bars on the shared scale, the combined
  percentile, the cross-engine check, and the dysbiosis anchor — and below all
  of that a decomposition table listing every feature with its measured value,
  reference percentile, direction, all four weight factors and sources, so the
  score is always reducible to what produced it.
* The confidence grade, its reasons, the validation figures and every caveat
  are printed — as fine print at the foot of the page, where they belong.
* An abstaining profile shows the reasons and what would need to change, not
  a blank or a zero.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    CondPageBreak,
    KeepTogether,
    PageBreak,
    Spacer,
    Table,
    TableStyle,
)

from openbiota.pdflinks import Paragraph, section_heading
from openbiota.pdfreport import (
    ACCENT,
    AMBER,
    AMBER_BG,
    CONTENT_WIDTH,
    CORAL,
    CORAL_BG,
    GREEN,
    GREEN_BG,
    INK_FAINT,
    NOT_SCORED,
    ORANGE,
    ORANGE_BG,
    PANEL_BG,
    PHYLA_PALETTE,
    RULE,
    SECTIONS,
    SLATE,
    SLATE_BG,
    Dial,
    Donut,
    Dot,
    HBar,
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

MODULE_LABELS = {
    "taxonomic": ("Taxonomic", "which species are present"),
    "functional": ("Functional", "gene capacity from Part A"),
    "ecological": ("Ecological", "diversity and richness"),
    "phenotype": ("Phenotype", "symptom-domain features"),
}


def _pct_status(percentile: float | None) -> Status:
    """A resemblance percentile on the shared scale. Higher = more like the pattern."""
    if percentile is None:
        return NOT_SCORED
    return status_for(percentile=percentile, higher_means="adverse")


def _anchor_status(anchor: Any) -> Status:
    if anchor is None or not anchor.valid:
        return Status("not computed", SLATE, SLATE_BG, "")
    if anchor.strongly_negative:
        return Status("broadly disturbed", CORAL, CORAL_BG, "")
    if anchor.band == "mildly dysbiotic":
        return Status("mildly disturbed", ORANGE, ORANGE_BG, "")
    if anchor.band == "indeterminate":
        return Status("indeterminate", AMBER, AMBER_BG, "")
    return Status("healthy range", GREEN, GREEN_BG, "")


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


# --------------------------------------------------------------------------- #
# community composition
# --------------------------------------------------------------------------- #


def community_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    similarity: Any | None,
    rpob_profile: dict[str, Any] | None,
    after_species: Callable[[list[Any]], None] | None = None,
    inv: Any = None,
    primary_phyla: dict[str, float] | None = None,
) -> None:
    """Species-level composition from the taxonomic engine, rpoB as fallback.

    ``after_species`` lets the caller drop extra flowables (the community-wide
    metric strips) between the species table and the gut-health index so the
    infographic page stays one continuous page.
    """
    story.extend(section_heading(section, "Your gut community", st["h1"]))
    taxonomy = similarity.taxonomy if similarity is not None else None

    if taxonomy is None:
        _rpob_fallback(story, st, rpob_profile)
        return

    from openbiota.scoring import ecological_metrics
    from openbiota.similarity import rescale_to_classified

    species = {k: v for k, v in rescale_to_classified(taxonomy.species).items() if v > 0}
    metrics = ecological_metrics(species)
    # Headline count and phylum picture come from the merged inventory and
    # the primary composition - the same numbers page one and the organism
    # section print - so the reader is never told 251 on one page and 92 on
    # the next. The diversity index beside them stays on the scoring lane,
    # because that is the only lane a reference percentile exists for, and
    # the tile says so.
    n_all = len(inv) if inv is not None and len(inv) else len(species)
    # The count on the reference's own catalogue - the same figure the
    # richness row below reports, so the page carries two numbers, not three.
    n_on_reference = len(species)
    if primary_phyla:
        total = sum(primary_phyla.values()) or 1.0
        phyla = sorted(((k, 100.0 * v / total) for k, v in primary_phyla.items()), key=lambda kv: -kv[1])
    else:
        phyla = sorted(rescale_to_classified(taxonomy.phyla).items(), key=lambda kv: -kv[1])
    firmicutes = next((v for k, v in phyla if k.lower().startswith(("firmicutes", "bacillota"))), 0.0)
    bacteroidetes = next(
        (v for k, v in phyla if k.lower().startswith(("bacteroidetes", "bacteroidota"))), 0.0
    )
    fb = firmicutes / bacteroidetes if bacteroidetes else None

    story.append(
        Paragraph(
            "Every reference catalogue available was searched, so an organism one "
            "catalogue cannot name is still found by another. These are "
            "<b>measured identifications</b>, not closest matches."
            + (
                (
                    " Human DNA had already been removed by the sequencing provider."
                    if taxonomy.host.upstream_removed
                    else f" Human DNA ({taxonomy.host.host_fraction:.2%} of fragments) "
                         "was removed first."
                )
                if taxonomy.host is not None
                else ""
            ),
            st["body"],
        )
    )
    story.append(Spacer(1, 2 * mm))
    story.append(
        _tiles(
            [
                ("Organisms detected", f"{n_all:,}", f"{n_on_reference:,} on the reference catalogue"),
                ("Diversity", f"{metrics['shannon_diversity']:.2f}", "Shannon index, reference-comparable"),
                ("Evenness", f"{metrics['evenness']:.2f}", "0 to 1, reference-comparable"),
                ("Firmicutes : Bacteroidetes", "—" if fb is None else f"{fb:.2f}", "two dominant groups"),
            ]
        )
    )
    story.append(Spacer(1, 5 * mm))

    # Donut + phylum bars side by side.
    shown = [(n, v) for n, v in phyla[:7] if v >= 0.1]
    parts = [(n, v, PHYLA_PALETTE[i % len(PHYLA_PALETTE)]) for i, (n, v) in enumerate(shown)]
    rest = sum(v for _, v in phyla) - sum(v for _, v in shown)
    if rest >= 0.1:
        parts.append(("Other", rest, PHYLA_PALETTE[-1]))
    donut = Donut(
        parts=parts, size=44 * mm,
        centre_value=f"{shown[0][1]:.0f}%" if shown else "—",
        centre_label=shown[0][0][:14] if shown else "",
    )
    rows: list[list[Any]] = [[Paragraph("COMPOSITION BY PHYLUM", st["label"]), "", ""]]
    for name, value, colour in parts:
        rows.append(
            [
                Paragraph(name, st["cell"]),
                HBar(width=CONTENT_WIDTH * 0.40, fraction=value / 100.0, colour=colour),
                Paragraph(f"<b>{value:.1f}%</b>", st["cell"]),
            ]
        )
    bars = Table(
        rows,
        colWidths=[CONTENT_WIDTH * 0.20, CONTENT_WIDTH * 0.42, CONTENT_WIDTH * 0.09],
        hAlign="LEFT",
    )
    bars.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("SPAN", (0, 0), (-1, 0)),
                ("TOPPADDING", (0, 0), (-1, -1), 2.0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.0),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    side = Table([[donut, bars]], colWidths=[CONTENT_WIDTH * 0.27, CONTENT_WIDTH * 0.73], hAlign="LEFT")
    side.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    story.append(side)

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("Most abundant organisms", st["h2"]))
    # From the primary composition, under current names - the same numbers
    # the organism section prints for the same organisms. The scoring-lane
    # readings differ slightly and use older names; showing those here put
    # "Bacteroides vulgatus 22%" on this page against "Phocaeicola vulgatus
    # 12%" four pages on.
    if inv is not None and len(inv.in_composition):
        top = [(o.display, o.percent) for o in inv.in_composition[:16]]
    else:
        top = [(k.replace("_", " "), v) for k, v in sorted(species.items(), key=lambda kv: -kv[1])[:16]]
    top_value = top[0][1] if top else 1.0
    half = (len(top) + 1) // 2
    left, right = top[:half], top[half:]
    pairs: list[list[Any]] = []
    for index in range(half):
        row: list[Any] = []
        for source in (left, right):
            if index < len(source):
                name, value = source[index]
                row.extend(
                    [
                        Paragraph(f"<i>{name}</i>", st["cell"]),
                        HBar(width=CONTENT_WIDTH * 0.12, fraction=value / top_value, colour=ACCENT,
                             height=2.4 * mm),
                        Paragraph(f"{value:.1f}%", st["cell_soft"]),
                    ]
                )
            else:
                row.extend([Paragraph("", st["cell"])] * 3)
        pairs.append(row)
    if pairs:
        gt = Table(
            pairs,
            colWidths=[CONTENT_WIDTH * w for w in (0.28, 0.13, 0.08, 0.28, 0.13, 0.08)],
            hAlign="LEFT",
        )
        gt.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("TOPPADDING", (0, 0), (-1, -1), 1.6),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
                    ("LEFTPADDING", (0, 0), (-1, -1), 0),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
                ]
            )
        )
        story.append(gt)

    story.append(Spacer(1, 3 * mm))
    if after_species is not None:
        after_species(story)
        story.append(Spacer(1, 3 * mm))
    # The section ends on the index that scores this community and the
    # organisms driving it. The marker-gene cross-check and the note about
    # proportions were fine print between the reader and that number.
    _anchor_block(story, st, similarity.anchor)


def _rpob_fallback(story: list[Any], st: dict[str, ParagraphStyle], profile: dict[str, Any] | None) -> None:
    if not profile or not profile.get("available"):
        story.append(
            Paragraph("Community composition was not available for this run.", st["body"])
        )
        return
    story.append(
        Paragraph(
            "The species-level profiler did not run for this sample, so composition below "
            "comes from the universal marker gene used to calibrate Part A. It is reliable at "
            "the phylum level and indicative at genus; names are closest matches rather than "
            "measured identifications.",
            st["body"],
        )
    )
    fb = profile.get("firmicutes_bacteroidetes_ratio")
    story.append(Spacer(1, 2 * mm))
    story.append(
        _tiles(
            [
                ("Diversity", f"{profile.get('shannon_index', 0):.2f}", "Shannon index"),
                ("Evenness", f"{profile.get('pielou_evenness') or 0:.2f}", "0 to 1"),
                ("Firmicutes : Bacteroidetes", f"{fb:.2f}" if fb else "—", "two dominant groups"),
                ("Distinct organisms", f"{profile.get('distinct_reference_organisms_hit', 0):,}",
                 "reference matches"),
            ]
        )
    )


def _anchor_block(story: list[Any], st: dict[str, ParagraphStyle], anchor: Any) -> None:
    story.append(CondPageBreak(62 * mm))  # heading, dial and prose travel together
    story.append(Paragraph("General gut-health index", st["h2"]))
    if anchor is None or not anchor.valid:
        story.append(
            Paragraph(
                anchor.note if anchor is not None else
                "Not computed: the taxonomic engine did not run.",
                st["body"],
            )
        )
        return
    words, _ = anchor_words(anchor)
    dial = Dial(value=anchor.score, label=words, width=40 * mm,
                sublabel=f"{anchor.n_health_taxa_present} healthy · {anchor.n_disease_taxa_present} disease markers")
    text = [
        Paragraph(
            "<b>GMWI2</b> is an independent, published index of general gut health trained on "
            "8,069 stool metagenomes across 54 studies. Positive values lean healthy, negative "
            "lean disturbed. It is repeated beside every pattern score in the next section "
            "because <b>a community that is generally disturbed will resemble many disease "
            "patterns at once</b>; the index tells you whether that is what is happening here.",
            st["body"],
        ),
        Paragraph(anchor.caution, st["small"]),
    ]
    block = Table([[dial, text]], colWidths=[46 * mm, CONTENT_WIDTH - 46 * mm], hAlign="LEFT")
    block.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                ("LEFTPADDING", (1, 0), (1, 0), 5 * mm),
            ]
        )
    )
    story.append(block)
    _anchor_drivers(story, st, anchor)


def _anchor_drivers(story: list[Any], st: dict[str, ParagraphStyle], anchor: Any) -> None:
    """The species moving the index, in two columns: pulling it down, holding it up.

    The index is a weighted sum over the organisms present, and every weight
    is known. This is the number that leads page one, so the organisms behind
    it belong on the page it links to - in the colour of their effect, not in
    a footnote.
    """
    contributions = list(getattr(anchor, "contributions", ()) or ())
    if not contributions:
        return
    down = [(c, w) for c, w in contributions if w < 0][:6]
    up = [(c, w) for c, w in contributions if w > 0][:6]

    def name(clade: str) -> str:
        last = clade.rsplit("|", 1)[-1]
        rank, _, label = last.partition("__")
        label = label.replace("_", " ")
        return f"<i>{label}</i>" if rank == "s" else f"<i>{label}</i> <font size='6' color='{_hex(INK_FAINT)}'>(genus)</font>"

    def column(title: str, rows: list[tuple[str, float]], colour: Any) -> Table:
        cells: list[list[Any]] = [[Paragraph(
            f"<font color='{_hex(colour)}'><b>{title}</b></font>", st["small"]), ""]]
        for clade, w in rows:
            cells.append([
                Paragraph(name(clade), st["cell"]),
                Paragraph(f"<font color='{_hex(colour)}'><b>{w:+.2f}</b></font>", st["cell"]),
            ])
        if not rows:
            cells.append([Paragraph(f"<font color='{_hex(INK_FAINT)}'>none</font>", st["cell"]), ""])
        t = Table(cells, colWidths=[(CONTENT_WIDTH / 2 - 6 * mm) * 0.78, (CONTENT_WIDTH / 2 - 6 * mm) * 0.22],
                  hAlign="LEFT")
        t.setStyle(TableStyle([
            ("SPAN", (0, 0), (1, 0)),
            ("LINEBELOW", (0, 0), (-1, 0), 0.7, colour),
            ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
            ("LINEBEFORE", (0, 0), (0, -1), 1.6, colour),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2),
        ]))
        return t

    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph(
        "<b>What is driving this index</b> &nbsp;<font color='#8C99A6' size='7'>the organisms "
        "present, weighted as the published model weights them</font>", st["h3"]))
    grid = Table(
        [[column("Pulling the index down", down, CORAL), column("Holding it up", up, GREEN)]],
        colWidths=[CONTENT_WIDTH / 2, CONTENT_WIDTH / 2], hAlign="LEFT",
    )
    grid.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (0, 0), 6 * mm),
        ("RIGHTPADDING", (1, 0), (1, 0), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    story.append(grid)
    lead = down[0] if down else None
    if lead is not None:
        story.append(Paragraph(
            f"<font size='7.2'>The single largest downward pull is {name(lead[0])} at "
            f"{lead[1]:+.2f}. Each weight is the model's coefficient for that organism when "
            "present; the index is their sum. These are the organisms whose change would move "
            "the number.</font>", st["small"]))


# --------------------------------------------------------------------------- #
# profile similarity pages
# --------------------------------------------------------------------------- #


def profile_pages(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    similarity: Any | None,
    profile_validation: dict[str, Any] | None,
) -> None:
    story.extend(section_heading(section, "Resemblance to published disease patterns", st["h1"]))
    story.append(
        Paragraph(
            "Researchers have described how the gut communities of people with certain "
            "conditions differ, on average, from those of healthy people. Each page below "
            "asks one question: <b>how closely does this sample resemble one of those "
            "published group-level patterns?</b> The answer is a percentile of a matched "
            "reference group, on the same scale as every other reading in this report.",
            st["body"],
        )
    )
    if similarity is None or not similarity.results:
        story.append(Spacer(1, 4 * mm))
        story.append(Paragraph("Profile similarity was not computed for this run.", st["body"]))
        return

    story.append(Spacer(1, 3 * mm))
    _overview_table(story, st, similarity)
    story.append(Spacer(1, 5 * mm))
    _how_to_read_a_profile(story, st, similarity)

    for result in similarity.results:
        story.append(PageBreak())
        _profile_page(story, st, result, similarity, profile_validation)


def _how_to_read_a_profile(
    story: list[Any], st: dict[str, ParagraphStyle], similarity: Any
) -> None:
    """Guidance that applies to every profile page, said once here."""
    story.append(Rule(CONTENT_WIDTH))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("How the pages that follow are built", st["h2"]))
    story.append(
        Paragraph(
            "Each pattern is scored in up to four separate parts, so you can see which kind "
            "of evidence is driving it rather than only a single number.",
            st["body"],
        )
    )
    rows = [
        [Paragraph(f"<b>{label}</b>", st["cell"]), Paragraph(text, st["small"])]
        for label, text in (
            (
                "Taxonomic",
                "Which species are present, from marker-gene profiling. Species level only "
                "— aggregating to genus reverses direction for several of these organisms.",
            ),
            (
                "Functional",
                "The pathway gene capacity from Part A, reused directly. An independent "
                "measurement, which is what makes the cross-engine check possible.",
            ),
            (
                "Ecological",
                "Diversity and richness. Deliberately given the least weight: these shift "
                "with diet, travel and antibiotics far more readily than with any disease.",
            ),
            (
                "Phenotype",
                "Symptom-linked features. Reported separately and never folded into the "
                "combined score, so a symptom pattern cannot be picked after seeing the data.",
            ),
        )
    ]
    table = Table(rows, colWidths=[CONTENT_WIDTH * 0.16, CONTENT_WIDTH * 0.84], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("TOPPADDING", (0, 0), (-1, -1), 2.4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 3 * mm))
    story.append(
        Paragraph(
            "Below each result is a table of every contributing feature with its measured "
            "value, its position in the comparison group, the direction the research reports, "
            "and the four factors that weight it. The score can be recomputed by hand from that "
            "table — nothing is hidden in the arithmetic. Where a meaningful score is not "
            "possible at all, the page says so and lists what would need to change.",
            st["body"],
        )
    )
    abstained = [r for r in similarity.results if r.abstention.abstained]
    if abstained:
        story.append(
            Paragraph(
                f"<b>{len(abstained)} of {len(similarity.results)} patterns were not scored "
                "for this sample.</b>",
                st["small"],
            )
        )


def _overview_table(story: list[Any], st: dict[str, ParagraphStyle], similarity: Any) -> None:
    widths = [CONTENT_WIDTH * w for w in (0.34, 0.38, 0.09, 0.19)]
    rows: list[list[Any]] = [
        [
            Paragraph("PATTERN", st["label"]),
            Paragraph("RESEMBLANCE", st["label"]),
            Paragraph("", st["label"]),
            Paragraph("READING", st["label"]),
        ]
    ]
    for result in similarity.results:
        pct = None if result.abstention.abstained else result.combined_percentile
        status = profile_status(result)
        rows.append(
            [
                Paragraph(f"<b>{result.profile.label}</b>", st["cell"]),
                PercentileBar(width=widths[1] - 6 * mm, percentile=pct, higher_means="adverse",
                              marker_colour=status.colour, height=6.4 * mm),
                Paragraph(f"<font color='{_hex(status.colour)}'>{_ordinal(pct)}</font>", st["pct"]),
                StatusChip(status, width=widths[3] - 4 * mm),
            ]
        )
    table = Table(rows, colWidths=widths, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("LEFTPADDING", (3, 0), (3, -1), 4),
                ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                ("LINEBELOW", (0, 1), (-1, -1), 0.4, RULE),
            ]
        )
    )
    story.append(table)
    anchor = similarity.anchor
    if anchor is not None and anchor.valid:
        words, colour = anchor_words(anchor)
        story.append(Spacer(1, 3 * mm))
        story.append(
            Paragraph(
                f"General gut-health index <b>{anchor.score:+.2f}</b> — "
                f"<font color='{_hex(colour)}'><b>{words}</b></font>. {anchor.caution}",
                st["small"],
            )
        )


def _profile_page(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    result: Any,
    similarity: Any,
    profile_validation: dict[str, Any] | None,
) -> None:
    profile = result.profile
    status = profile_status(result)
    combined_pct = None if result.abstention.abstained else result.combined_percentile

    head = Table(
        [
            [
                Paragraph(profile.label, st["h1"]),
                StatusChip(status, width=36 * mm, height=6.2 * mm, font_size=7.2),
            ]
        ],
        colWidths=[CONTENT_WIDTH - 38 * mm, 38 * mm],
        hAlign="LEFT",
    )
    head.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    story.append(head)
    story.append(
        Paragraph(
            f"Pattern version <b>{profile.version}</b> <font color='#8C99A6'>({profile.name})</font> &middot; "
            f"{profile.status.replace('_', ' ')}",
            st["small"],
        )
    )
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(" ".join(profile.summary.split()), st["body"]))
    story.append(Spacer(1, 2 * mm))

    # -- headline bar ---------------------------------------------------- #
    headline = Table(
        [
            [
                Paragraph(
                    "<b>COMBINED RESEMBLANCE</b><br/>"
                    f"<font color='{_hex(INK_FAINT)}' size='7'>percentile of the matched "
                    "reference group</font>",
                    st["cell"],
                ),
                PercentileBar(width=CONTENT_WIDTH * 0.46, percentile=combined_pct,
                              higher_means="adverse", marker_colour=status.colour,
                              height=9 * mm, show_scale=True, track=4.2 * mm),
                Paragraph(
                    f"<font color='{_hex(status.colour)}'>{_ordinal(combined_pct)}</font>"
                    if combined_pct is not None else "—",
                    st["metric"],
                ),
            ]
        ],
        colWidths=[CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.52, CONTENT_WIDTH * 0.18],
        hAlign="LEFT",
    )
    headline.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("BACKGROUND", (0, 0), (-1, -1), PANEL_BG),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("ROUNDEDCORNERS", [4, 4, 4, 4]),
            ]
        )
    )
    story.append(headline)
    story.append(Spacer(1, 3 * mm))

    # -- abstention ---------------------------------------------------- #
    if result.abstention.abstained:
        story.append(Paragraph("No score was produced", st["h2"]))
        story.append(
            Paragraph(
                "The engine <b>declined to score this pattern</b>: one or more conditions for "
                "a meaningful score were not met. The module values below are shown for "
                "information only.",
                st["body"],
            )
        )
        story.append(Paragraph("<b>Why</b>", st["small"]))
        for reason in result.abstention.triggered:
            story.append(Paragraph(f"•&nbsp;&nbsp;{reason}", st["body"]))
        story.append(Paragraph("<b>What would change this</b>", st["small"]))
        for remedy in result.abstention.remedies:
            story.append(Paragraph(f"•&nbsp;&nbsp;{remedy}", st["body"]))
        story.append(Spacer(1, 2 * mm))

    # -- module bars --------------------------------------------------- #
    story.append(Paragraph("What is driving it", st["h2"]))
    rows: list[list[Any]] = []
    for name in ("taxonomic", "functional", "ecological", "phenotype"):
        module = result.modules.get(name)
        if module is None:
            continue
        label, hint = MODULE_LABELS[name]
        pct = module.percentile if module.reportable else None
        share = (
            f"{module.weight_in_combined:.0%} of combined" if module.weight_in_combined
            else "not in combined"
        )
        m_status = _pct_status(pct)
        rows.append(
            [
                Paragraph(
                    f"<b>{label}</b><br/><font color='{_hex(INK_FAINT)}' size='7'>{hint} "
                    f"&middot; {module.n_measured}/{len(module.features)} features &middot; "
                    f"{share}</font>",
                    st["cell"],
                ),
                PercentileBar(width=CONTENT_WIDTH * 0.36, percentile=pct, higher_means="adverse",
                              marker_colour=m_status.colour, height=6.4 * mm),
                Paragraph(
                    f"<font color='{_hex(m_status.colour)}'>{_ordinal(pct)}</font>"
                    if pct is not None else "—",
                    st["pct"],
                ),
                StatusChip(m_status, width=CONTENT_WIDTH * 0.17 - 4 * mm),
            ]
        )
    table = Table(
        rows,
        colWidths=[CONTENT_WIDTH * 0.38, CONTENT_WIDTH * 0.38, CONTENT_WIDTH * 0.07,
                   CONTENT_WIDTH * 0.17],
        hAlign="LEFT",
    )
    table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("LEFTPADDING", (3, 0), (3, -1), 4),
                ("LINEBELOW", (0, 0), (-1, -2), 0.4, RULE),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 3 * mm))

    # -- anchor + cross-engine ----------------------------------------- #
    anchor = similarity.anchor
    words, colour = anchor_words(anchor)
    side_by_side: list[Any] = [
        [
            Paragraph("GENERAL GUT-HEALTH INDEX", st["label"]),
            Paragraph(
                (f"<b>{anchor.score:+.2f}</b> &nbsp;<font color='{_hex(colour)}'><b>"
                 f"{words.upper()}</b></font>") if anchor is not None and anchor.valid else "—",
                st["cell"],
            ),
            Paragraph(
                anchor.caution if anchor is not None else
                "Not computed: the taxonomic engine did not run.",
                st["small"],
            ),
        ]
    ]
    checks: list[Any] = []
    for check in result.cross_engine:
        c_colour = GREEN if check.concordant else (SLATE if check.verdict == "UNAVAILABLE" else CORAL)
        checks.extend(
            [
                Paragraph("CROSS-ENGINE CHECK", st["label"]),
                Table(
                    [[Dot(c_colour, size=2.4 * mm), Paragraph(f"<b>{check.verdict.title()}</b>", st["cell"])]],
                    colWidths=[4 * mm, 40 * mm], hAlign="LEFT",
                    style=TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0),
                                      ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                      ("TOPPADDING", (0, 0), (-1, -1), 0),
                                      ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]),
                ),
                Paragraph(f"Species: {check.taxonomic_detail}", st["small"]),
                Paragraph(f"Gene panel: {check.functional_detail}", st["small"]),
                Paragraph(check.message, st["small"]),
            ]
        )
    if checks:
        two = Table(
            [[side_by_side[0], checks]],
            colWidths=[CONTENT_WIDTH * 0.5, CONTENT_WIDTH * 0.5],
            hAlign="LEFT",
        )
        two.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                                 ("LEFTPADDING", (0, 0), (0, 0), 0),
                                 ("LEFTPADDING", (1, 0), (1, 0), 6),
                                 ("LINEBEFORE", (1, 0), (1, 0), 0.5, RULE)]))
        story.append(two)
    else:
        for item in side_by_side[0]:
            story.append(item)

    if profile.duration_dependent:
        story.append(Spacer(1, 2.5 * mm))
        if result.stratum:
            stratum = next((s for s in profile.strata if s.name == result.stratum), None)
            text = (
                f"<b>Illness duration stratum: {stratum.label if stratum else result.stratum}.</b> "
                + (" ".join(stratum.note.split()) if stratum else "")
            )
        else:
            text = (
                "<b>Illness duration was not supplied.</b> The published pattern differs "
                "between short-term and long-term illness — marked microbial changes early, "
                "largely resolved microbial changes after a decade with persisting metabolic "
                "ones — and the two must not be averaged. Both readings therefore apply. A low "
                "taxonomic score in long-standing illness is consistent with the published "
                "pattern, not evidence against it."
            )
        story.append(Paragraph(text, st["body"]))

    # -- decomposition table ------------------------------------------- #
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("Every feature that produced this score", st["h2"]))
    story.append(
        Paragraph(
            "Direction is what the research reports for the condition (▲ higher, ▼ lower). "
            "<b>v</b> is how far this sample moves in that direction, from −1 (opposite) to "
            "+1 (strongly concordant). The four factors weight each feature by evidence "
            "strength (w), cohort independence (q), condition specificity (s) and direction "
            "agreement across studies (c); all four are in the versioned profile.",
            st["small"],
        )
    )
    story.append(Paragraph(
        f"Every species value here is the reference catalogue's own measurement (MetaPhlAn 3, the catalogue the reference cohort was profiled with), because that is what makes a percentile against that cohort mean anything. It is not the organism's share of the community in section {SECTIONS['catalogue']}, which pools nine detection methods and the current catalogues: a species that catalogue has no genome for reads as <i>absent</i> here and still has a share there; where that is so, the pooled share is printed beneath. Section {SECTIONS['catalogue']} is the account of what is present; this table is the account of what produced this score.",
        st["small"]))
    story.append(Spacer(1, 1.5 * mm))
    _decomposition_table(story, st, result)

    # -- fine print: confidence, validation, caveats, sources ---------- #
    notes: list[str] = []
    if profile.status == "validation_control":
        notes.append(
            "<b>Validation control.</b> This pattern is included to check that the scoring "
            "engine recovers a known signal on labelled data. It is not a cancer screen and not "
            "a risk estimate; validated clinical tests exist for this condition (faecal "
            "immunochemical testing, stool DNA testing, colonoscopy)."
        )
    if result.abstention.abstained:
        notes.append("<b>Confidence:</b> not graded — no score was produced.")
    else:
        reasons = "; ".join(result.confidence.reasons[:3])
        more = max(0, len(result.confidence.reasons) - 3)
        notes.append(
            f"<b>Confidence grade {result.confidence.grade.lower()}.</b> "
            + (f"Why not higher: {reasons}"
               + (f"; +{more} more, all in results.json." if more else ".")
               if reasons else "All confidence checks passed.")
        )
    notes.append(
        "The combined score weights taxonomic 50 / functional 35 / ecological 15 — an "
        "engineering starting point, not a validated weighting."
    )
    if profile_validation:
        notes.append(_validation_text(profile.name, profile_validation))
    notes.extend(" ".join(caveat.split()) for caveat in profile.caveats)
    notes.append(f"Sources: {' '.join(profile.citation.split())}")
    fine_print(story, st, notes)


def _validation_text(name: str, validation: dict[str, Any]) -> str:
    results = [r for r in validation.get("results", []) if r.get("profile") == name]
    if not results:
        return (
            "<b>Measured discrimination: not available.</b> No labelled cohort exists for "
            "this pattern, so its ability to separate cases from controls is unmeasured."
        )
    own = next((r for r in results if r["condition"].lower() == name.lower()), None)
    parts = []
    for r in sorted(results, key=lambda x: -(x.get("auc") or 0)):
        ci = r.get("auc_ci95")
        parts.append(
            f"{r['condition']} {r['auc']:.2f}" + (f" ({ci[0]:.2f}–{ci[1]:.2f})" if ci else "")
        )
    verdict = validation.get("specificity_verdicts", {}).get(name, "")
    return (
        "<b>Measured on labelled cohorts</b> (AUC, 0.5 = chance, 1.0 = perfect): "
        + "; ".join(parts) + ". "
        + (f"{verdict[0].upper() + verdict[1:]}." if verdict else "")
        + (
            " The published 0.80+ for this condition comes from learned models over hundreds "
            "of features; a handful of prespecified markers cannot match it, and is not meant to."
            if own and (own.get("auc") or 0) < 0.75
            else ""
        )
    )


def _decomposition_table(story: list[Any], st: dict[str, ParagraphStyle], result: Any) -> None:
    from openbiota import pdfcontext as _pdfcontext

    head = [
        Paragraph(x, st["label"])
        for x in ("FEATURE", "MEASURED", "REF PCT", "DIR", "w", "q", "s", "c", "v", "SOURCES")
    ]
    rows: list[list[Any]] = [head]
    spans: list[tuple[int, str]] = []
    for name in ("taxonomic", "functional", "ecological", "phenotype"):
        module = result.modules.get(name)
        if module is None:
            continue
        spans.append((len(rows), MODULE_LABELS[name][0]))
        rows.append([Paragraph(f"<b>{MODULE_LABELS[name][0]} module</b>", st["cell"])] + [""] * 9)
        for f in module.features:
            if f.raw_value is None:
                measured = "absent"
                pooled = (_pdfcontext.pooled_note(f.feature.name, catalogue_value=None)
                          if f.feature.engine == "metaphlan" else None)
                if pooled:
                    measured += f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>{pooled}</font>"
            elif f.feature.engine == "metaphlan":
                measured = f"{f.raw_value:.3f}%"
            elif f.feature.engine == "diamond":
                measured = f"{f.raw_value:.1f}/100"
            else:
                measured = f"{f.raw_value:.2f}"
            v = f.v
            v_colour = (
                _hex(CORAL) if v is not None and v > 0.33 else
                _hex(GREEN) if v is not None and v < -0.33 else _hex(SLATE)
            )
            sources = ", ".join(f.feature.sources)
            rows.append(
                [
                    Paragraph(
                        f"<i>{f.feature.name.replace('_', ' ')}</i>"
                        if f.feature.engine == "metaphlan"
                        else f.feature.name.replace("_", " "),
                        st["cell"],
                    ) if not f.feature.measured_as else Paragraph(
                        f"<i>{f.feature.name.replace('_', ' ')}</i>, measured as "
                        f"<i>{f.feature.measured_as.replace('_', ' ')}</i>",
                        st["cell"],
                    ),
                    Paragraph(measured, st["cell"]),
                    Paragraph(f"{f.percentile:.0f}" if f.percentile is not None else "—", st["cell"]),
                    Paragraph("▲" if f.feature.direction == "increased" else "▼", st["cell"]),
                    Paragraph(f"{f.feature.w:.2f}", st["cell"]),
                    Paragraph(f"{f.feature.q:.2f}", st["cell"]),
                    Paragraph(f"{f.feature.s:.2f}", st["cell"]),
                    Paragraph(f"{f.feature.c:.2f}", st["cell"]),
                    Paragraph(
                        f"<font color='{v_colour}'><b>{v:+.2f}</b></font>" if v is not None else "—",
                        st["cell"],
                    ),
                    Paragraph(sources, st["fine"]),
                ]
            )
    widths = [0.24, 0.10, 0.07, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.27]
    table = Table(rows, colWidths=[CONTENT_WIDTH * w for w in widths], hAlign="LEFT", repeatRows=1)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 1.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 1.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 1.5),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
    ]
    for row_index, _ in spans:
        style.append(("SPAN", (0, row_index), (-1, row_index)))
        style.append(("BACKGROUND", (0, row_index), (-1, row_index), PANEL_BG))
    table.setStyle(TableStyle(style))
    story.append(table)


# --------------------------------------------------------------------------- #
# what this cannot tell you
# --------------------------------------------------------------------------- #

CANNOT_TELL_YOU: tuple[str, ...] = (
    "<b>It does not measure metabolite levels.</b> Carrying the gene means the bacteria "
    "<i>can</i> make the metabolite. How much they actually make depends on your diet, your "
    "gut transit time and much else. Measuring the metabolites themselves requires a blood "
    "or stool chemistry test.",
    "<b>It does not diagnose or predict disease.</b> No threshold in this report corresponds "
    "to a clinical diagnosis, because no such threshold has been published for any of these "
    "pathways.",
    "<b>It is a snapshot.</b> Gut composition shifts with diet over days to weeks. A single "
    "measurement is far less informative than the same test repeated over time, which is why "
    "the trend matters more than any one value.",
    "<b>Pattern resemblance is resemblance, not diagnosis.</b> A resemblance percentile says how "
    "closely this community resembles the <i>average</i> pattern published for a group of "
    "patients. It is not a probability of having the condition. The underlying associations "
    "largely do not distinguish one condition from general gut disturbance — the same species "
    "shift in inflammatory bowel disease, after antibiotics, and with a low-fibre diet — which "
    "is why every pattern score is shown beside the general gut-health index.",
)


def limits_page(story: list[Any], st: dict[str, ParagraphStyle], *, section: int) -> None:
    story.extend(section_heading(section, "What this report cannot tell you", st["h1"]))
    story.append(
        Paragraph(
            "The four boundaries of what a stool metagenome can show.",
            st["body"],
        )
    )
    for text in CANNOT_TELL_YOU:
        story.append(Paragraph(f"•&nbsp;&nbsp;{text}", st["body"]))
        story.append(Spacer(1, 1.5 * mm))
    fine_print(
        story,
        st,
        [
            "Research use only. Nothing in this report is a diagnosis, a risk estimate or a "
            "treatment recommendation. If you have symptoms or concerns, the right next step is "
            "a clinician and the validated tests they can order."
        ],
    )


# --------------------------------------------------------------------------- #
# performance extension: labelled-cohort validation and specificity matrix
# --------------------------------------------------------------------------- #


def _worst_off_target(matrix: Mapping[str, Any]) -> dict[str, Any] | None:
    """The strongest cross-condition signal actually present in the matrix.

    Returns the profile/condition pair with the highest AUC against a
    condition that is not the profile's own, or None when the matrix has
    nothing to say. Used so the worked example under the specificity table
    describes a row the reader can see rather than a hardcoded claim about
    a profile that may not be scored at all.
    """
    best: dict[str, Any] | None = None
    for profile_name, row in (matrix or {}).items():
        if not isinstance(row, Mapping):
            continue
        for condition, entry in row.items():
            if not isinstance(entry, Mapping):
                continue
            auc = entry.get("auc")
            if not isinstance(auc, (int, float)):
                continue
            # Skip a profile scored against its own condition.
            if str(condition).lower() in str(profile_name).lower():
                continue
            if best is None or auc > best["auc"]:
                best = {
                    "profile": profile_name,
                    "condition": condition,
                    "auc": float(auc),
                }
    # Only worth printing when it is meaningfully above chance.
    return best if best and best["auc"] >= 0.55 else None


def performance_extension(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    profile_validation: dict[str, Any] | None,
    similarity: Any | None,
) -> None:
    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("Part B — how the resemblance scoring was tested", st["h2"]))
    manifest = (similarity.cohort_manifest if similarity is not None else {}) or {}
    n = manifest.get("n_samples")
    comp = manifest.get("composition", {})
    if n:
        story.append(
            Paragraph(
                f"Species abundances are compared against <b>{n:,} healthy adult stool "
                f"metagenomes</b> from {comp.get('n_studies', '?')} published studies "
                f"(curatedMetagenomicData snapshot {manifest.get('snapshot')}), measured the "
                "same way as your sample so the two are directly comparable. "
                f"Reference identifier {manifest.get('manifest_id')}.",
                st["body"],
            )
        )
    if not profile_validation:
        story.append(
            Paragraph("Labelled-cohort validation results were not available.", st["body"])
        )
        return

    results = profile_validation.get("results", [])
    if results:
        story.append(Spacer(1, 2 * mm))
        story.append(
            Paragraph(
                "The scoring engine was run end to end on public cohorts that carry disease "
                "labels, and asked to separate cases from controls. AUC is the standard "
                "measure: 0.5 is a coin flip, 1.0 is perfect. <b>The same pattern was also "
                "scored against unrelated conditions</b> — a pattern that scores high on "
                "everything is measuring general gut disturbance, not its own condition.",
                st["body"],
            )
        )
        story.append(Spacer(1, 2 * mm))
        matrix = profile_validation.get("specificity_matrix", {}).get("matrix", {})
        conditions = sorted({c for row in matrix.values() for c in row})
        head = [Paragraph("PATTERN", st["label"])] + [
            Paragraph(f"vs {c}", st["label"]) for c in conditions
        ]
        rows: list[list[Any]] = [head]
        for profile_name in sorted(matrix):
            row = [Paragraph(f"<b>{profile_name}</b>", st["cell"])]
            for c in conditions:
                entry = matrix[profile_name].get(c) or {}
                auc = entry.get("auc")
                ci = entry.get("auc_ci95")
                own = c.lower() == profile_name.lower()
                text = (
                    f"<b>{auc:.2f}</b>" if own else f"{auc:.2f}"
                ) if auc is not None else "—"
                if ci and auc is not None:
                    text += f"<br/><font size='6.5' color='{_hex(INK_FAINT)}'>{ci[0]:.2f}–{ci[1]:.2f}</font>"
                row.append(Paragraph(text, st["cell"]))
            rows.append(row)
        widths = [CONTENT_WIDTH * 0.28] + [CONTENT_WIDTH * 0.72 / max(len(conditions), 1)] * len(conditions)
        table = Table(rows, colWidths=widths, hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
                ]
            )
        )
        story.append(KeepTogether([table]))
        story.append(Spacer(1, 2 * mm))
        for name, verdict in profile_validation.get("specificity_verdicts", {}).items():
            story.append(Paragraph(f"<b>{name}:</b> {verdict}.", st["small"]))
        unmeasured = [
            r.profile.name for r in (similarity.results if similarity else [])
            if not any(x.get("profile") == r.profile.name for x in results)
        ]
        if unmeasured:
            story.append(Spacer(1, 1.5 * mm))
            story.append(
                Paragraph(
                    f"<b>Unmeasured:</b> {', '.join(unmeasured)}. No labelled cohort with "
                    "outcome keys was available for these patterns, so their discrimination is "
                    "unknown and no accuracy figure is printed for them.",
                    st["small"],
                )
            )
    story.append(Spacer(1, 2 * mm))
    # The worked example is derived from the matrix that was actually
    # printed, not asserted. It used to name the colorectal pattern
    # unconditionally, and no colorectal row exists in the shipped
    # validation file - only longcovid and mecfs are scored - so the
    # report was drawing a conclusion from a row the reader could not see.
    _off_target = _worst_off_target(
        (profile_validation or {}).get("specificity_matrix", {}).get("matrix", {})
    )
    story.append(
        Paragraph(
            "Reading the table: a pattern should score highest against its own condition and "
            "near 0.5 against the others."
            + (
                f" Here the <b>{_off_target['profile']}</b> pattern scores <i>higher</i> "
                f"against {_off_target['condition']} (AUC {_off_target['auc']:.2f}) than the "
                "0.5 a pattern with no cross-condition signal would give, which means part of "
                "what it detects is general disturbance rather than that one condition."
                if _off_target
                else ""
            )
            + " That is why the general gut-health index accompanies every resemblance score.",
            st["small"],
        )
    )
