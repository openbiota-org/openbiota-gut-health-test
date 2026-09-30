"""Report pages for the inflammatory-skin research panels — spec v4.2 §7.2.

Two parts. A front-loaded section carries *every* skin reading before any
detailed explanation, paginating rather than truncating to a top-N. Detail
pages then give each panel its declared features one by one, including the
slots that could not be measured, because a missing slot is a fact about
coverage and hiding it would make a partial panel look complete.

Presentation rules that are not cosmetic
----------------------------------------
* The scale is labelled **pattern concordance**, low to high, in neutral
  greys. No green "healthy" and no red "diseased": the number is agreement
  with one study's reported directions, not a verdict.
* Nulls render as an em dash with a reason. Never zero, never a green tick.
  An unavailable score is not a negative result.
* Where point coverage is incomplete, two separately labelled tracks appear:
  the measured-feature index, and the possible full-panel range. They have
  different scopes, and the second is not a confidence interval around the
  first, so they are never drawn as one bar with whiskers.
* One-feature panels say "single-marker result" in words.
* Transported panels say the pattern came from 16S and was applied to
  shotgun measurements, so a reader cannot mistake the lane for a native one.
* Everything stays legible in greyscale and has a text equivalent.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

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
    linked_cell,
    section_heading,
    skin_dest,
)

CONTENT_WIDTH = 170 * mm

INK = colors.HexColor("#15222E")
INK_SOFT = colors.HexColor("#4B5B6B")
INK_FAINT = colors.HexColor("#8C99A6")
RULE = colors.HexColor("#E1E7EC")
PANEL_BG = colors.HexColor("#F5F8FA")
SLATE = colors.HexColor("#64748B")
SLATE_BG = colors.HexColor("#EEF1F4")

#: Statuses that carry a number.
_SCORED = frozenset({"computed", "computed_partial", "computed_research", "computed_limited_features"})

#: Plain-English rendering of every reason a panel produced no number.
_STATUS_WORDS: Mapping[str, str] = {
    "reference_unavailable": "no compatible reference group",
    "population_unsupported": "outside the supported population",
    "assay_incompatible": "measurement not compatible",
    "qc_failed": "sequencing quality gate not met",
    "features_unavailable": "none of the organisms could be measured",
    "interval_only": "only bounded measurements available",
    "insufficient_stool_signature": "stool-DNA signature not established",
    "strain_bundle_pending": "strain-level reference not available",
    "blocked_unresolved_preprocessing_and_artifacts": "source method unresolved",
}

_TIER_WORDS: Mapping[str, str] = {
    "same_center_tested_candidate": "tested in the same centre, not externally validated",
    "single_cohort_selected_subset": "one cohort, selected subset of features",
    "small_confounded_cohort": "small cohort with heavy treatment confounding",
    "shared_psoriasis_signal": "pooled psoriasis and psoriatic arthritis contrast",
    "single_cohort_nominal": "one cohort, nominal significance only",
    "assay_transported_exploratory": "16S-derived pattern applied to shotgun measurements",
    "cross_study_control_confounding": "cases and controls from different studies",
    "single_feature_nominal": "one feature at nominal significance",
    "subtype_contrast_capacity_observation": "comparison between two disease subtypes",
}


def _p(text: str, style: ParagraphStyle) -> Paragraph:
    return Paragraph(text, style)


def _hex(c: colors.Color) -> str:
    return f"#{int(c.red * 255):02X}{int(c.green * 255):02X}{int(c.blue * 255):02X}"


def _panels(skin: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Every numeric panel record, in the registry's report order."""
    if not skin:
        return []
    from openbiota.skin import PANELS

    out = []
    for panel in PANELS:
        record = skin.get(panel.panel_id)
        if isinstance(record, dict):
            out.append(record)
    return out


def _panel_id(record: Mapping[str, Any]) -> str:
    """The panel identifier, under either key the result layer has used."""
    return str(record.get("panel_id") or record.get("profile_id") or "")


def _index_text(record: Mapping[str, Any]) -> str:
    """The index as words, so the graphic is never the only carrier."""
    index = record.get("index", record.get("score_0_100"))
    if index is None:
        status = str(record.get("status", ""))
        return f"&mdash; <font size=7>({_STATUS_WORDS.get(status, status.replace('_', ' '))})</font>"
    return f"<b>{index:.0f}</b> <font size=7 color='{_hex(INK_FAINT)}'>of 100</font>"


def _coverage_text(record: Mapping[str, Any]) -> str:
    usable = record.get("usable_count", 0)
    total = record.get("panel_count", 0)
    if not total:
        return "&mdash;"
    label = f"{usable}/{total}"
    if total == 1:
        return f"{label} <font size=6.5>single-marker result</font>"
    if usable == 1:
        return f"{label} <font size=6.5>one feature only</font>"
    return label


def _bounds_text(record: Mapping[str, Any]) -> str:
    bounds = record.get("full_panel_bounds")
    if not bounds or len(bounds) != 2:
        return "&mdash;"
    return f"{bounds[0]:.0f}&ndash;{bounds[1]:.0f}"


def skin_glance(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    skin: Mapping[str, Any] | None,
    skin_axis: Any | None = None,
) -> None:
    """Front-loaded section: every skin reading, paginated, nothing dropped."""
    story.extend(section_heading(section, "Gut-skin axis", st["h1"]))
    axis_panel(story, st, skin_axis)
    story.append(_p(
        "Researchers have compared the stool communities of people with several inflammatory skin "
        "diseases against healthy controls. This section asks, for each published list of organisms: "
        "<b>does this sample move in the direction that study reported?</b> "
        "There is no established, uniquely diagnostic stool signature for any of these conditions. "
        "A high reading is not a diagnosis and a low one does not rule anything out.",
        st["body"],
    ))
    records = _panels(skin)
    if not records:
        story.append(Spacer(1, 3 * mm))
        story.append(_p("Skin research panels were not computed for this run.", st["body"]))
        return

    story.append(Spacer(1, 2 * mm))
    story.append(_p(
        "The scale is <b>pattern concordance</b> from low to high, drawn in neutral greys because it "
        "is agreement with a study's reported directions and not a health grade. It carries no per "
        "cent sign, because it is not a percentage or a probability. Where some organisms could not "
        "be measured, two separate figures appear: the index over the features that <i>were</i> "
        "measured, and the range the whole declared panel could occupy. The second is a coverage "
        "range, not a margin of error on the first.",
        st["small"],
    ))
    story.append(Spacer(1, 2.5 * mm))

    from openbiota.skin import FAMILIES

    header = [
        _p("CONDITION AND COMPARISON", st["label"]),
        _p("PANEL", st["label"]),
        _p("INDEX", st["label"]),
        _p("FEATURES", st["label"]),
        _p("PANEL RANGE", st["label"]),
        _p("EVIDENCE", st["label"]),
    ]
    widths = [
        CONTENT_WIDTH * 0.235, CONTENT_WIDTH * 0.215, CONTENT_WIDTH * 0.105,
        CONTENT_WIDTH * 0.115, CONTENT_WIDTH * 0.105, CONTENT_WIDTH * 0.225,
    ]
    rows: list[list[Any]] = [header]
    for record in records:
        family = FAMILIES.get(str(record.get("profile_family") or ""))
        contrast = str(record.get("contrast") or "").replace("_", " ")
        panel_id = _panel_id(record)
        rows.append([
            linked_cell(
                RowLink(href=skin_dest(panel_id, detail=True), anchor=skin_dest(panel_id)),
                _p(f"<b>{family.label if family else record.get('profile_family', '')}</b><br/>"
                   f"<font size=6.5 color='{_hex(INK_FAINT)}'>{contrast}</font>", st["cell"]),
            ),
            _p(str(record.get("label") or _panel_id(record)), st["cell"]),
            _p(_index_text(record), st["cell"]),
            _p(_coverage_text(record), st["cell"]),
            _p(_bounds_text(record), st["cell"]),
            _p(_TIER_WORDS.get(str(record.get("evidence_tier", "")),
                               str(record.get("evidence_tier", "")).replace("_", " ")), st["cell"]),
        ])
    table = LinkedTable(rows, colWidths=widths, hAlign="LEFT", repeatRows=1)
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), PANEL_BG),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    story.append(table)

    # Registered skin tasks that produce no number are recorded in
    # results.json (``skin.unavailable_tasks``) for the audit trail; they are
    # not printed here. A list of things that were not scored is not a
    # result a reader can use, and it read as one.

    story.append(KeepTogether([
        Spacer(1, 2.5 * mm),
        _p(
            "Panels for the same condition are shown side by side and never averaged into "
            "one condition score: they come from different cohorts, they disagree in "
            "places, and averaging would hide that. One organism appearing in several "
            "panels is one measurement, not several agreeing results.",
            st["fine"],
        ),
    ]))


def _feature_table(
    story: list[Any], st: dict[str, ParagraphStyle], record: Mapping[str, Any],
) -> None:
    """Every declared slot, measured or not."""
    observations: Sequence[Mapping[str, Any]] = (
        record.get("observations") or record.get("feature_results") or []
    )
    if not observations:
        story.append(_p(
            "No feature-level measurements were produced for this panel, so there is nothing to "
            "break down. The declared organisms are listed in the panel definition and all of them "
            "remain in the coverage denominator.",
            st["small"],
        ))
        return
    header = [
        _p("ORGANISM OR GENE FAMILY", st["label"]),
        _p("MEASURED", st["label"]),
        _p("STUDY DIRECTION", st["label"]),
        _p("CONTRIBUTION", st["label"]),
        _p("NOTE", st["label"]),
    ]
    widths = [CONTENT_WIDTH * 0.25, CONTENT_WIDTH * 0.14, CONTENT_WIDTH * 0.13,
              CONTENT_WIDTH * 0.13, CONTENT_WIDTH * 0.35]
    rows: list[list[Any]] = [header]
    for obs in observations:
        name = str(obs.get("measured_as") or obs.get("id") or "")
        state = str(obs.get("state", ""))
        abundance = obs.get("abundance_fraction")
        measured = (
            f"{abundance * 100:.3g}%" if isinstance(abundance, (int, float)) and abundance > 0
            else state.replace("_", " ")
        )
        direction = "higher in cases" if obs.get("direction") == 1 else "lower in cases"
        contribution = obs.get("contribution")
        bounds = obs.get("contribution_bounds")
        if isinstance(contribution, (int, float)):
            contrib = f"{contribution:.0f}"
        elif bounds:
            contrib = f"{bounds[0]:.0f}&ndash;{bounds[1]:.0f}<br/><font size=6>range only</font>"
        else:
            contrib = "&mdash;"
        # For a slot the catalogue cannot resolve, the declared mapping reason
        # (a SILVA group, an SGB namespace, a clade) explains more than the
        # generic "no slot" note does.
        mapping_reason = str(obs.get("mapping_reason") or "")
        note = str(obs.get("note") or mapping_reason)
        if state == "not_in_database" and mapping_reason:
            note = mapping_reason
        rows.append([
            _p(f"<i>{name.replace('_', ' ')}</i>", st["cell"]),
            _p(measured, st["cell"]),
            _p(direction, st["cell"]),
            _p(contrib, st["cell"]),
            _p(note, st["cell"]),
        ])
    table = Table(rows, colWidths=widths, hAlign="LEFT", repeatRows=1)
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), PANEL_BG),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, RULE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    story.append(table)


def _panel_detail(
    story: list[Any], st: dict[str, ParagraphStyle], record: Mapping[str, Any], *, section: int, n: int,
) -> None:
    from openbiota.pdfatlas import AlignmentBar
    from openbiota.skin import BY_ID, FAMILIES

    panel = BY_ID.get(_panel_id(record))
    family = FAMILIES.get(str(record.get("profile_family") or ""))
    story.append(Anchor(skin_dest(_panel_id(record), detail=True), above=6 * mm))
    story.append(back_link_line(skin_dest(_panel_id(record)), "at a glance"))
    story.extend(section_heading(
        section, f"{record.get('label', '')}", st["h1"], sub=n,
    ))

    if family is not None:
        story.append(_p(f"<b>{family.plain_label.capitalize()}.</b> {family.summary}", st["body"]))
        story.append(Spacer(1, 1.5 * mm))

    index = record.get("index", record.get("score_0_100"))
    status = str(record.get("status", ""))
    if index is not None:
        story.append(AlignmentBar(
            width=CONTENT_WIDTH, value=float(index),
            bounds=tuple(record["full_panel_bounds"]) if record.get("full_panel_bounds") else None,
            low_label="0 — low pattern concordance",
            mid_label="50 — at the reference centre",
            high_label="high pattern concordance — 100",
        ))
        usable, total = record.get("usable_count", 0), record.get("panel_count", 0)
        bounds = record.get("full_panel_bounds") or [0.0, 100.0]
        if total and usable < total:
            story.append(_p(
                f"<b>Measured-feature index {index:.0f}</b>, over the {usable} of {total} declared "
                f"features that could be measured. <b>Possible full-panel range "
                f"{bounds[0]:.0f}&ndash;{bounds[1]:.0f}</b>, which is where the index could sit once "
                f"the {total - usable} unmeasured feature{'s' if total - usable != 1 else ''} are "
                "accounted for. These two describe different things: the range is not a margin of "
                "error around the index, and the index can fall outside it without any arithmetic "
                "being wrong.",
                st["small"],
            ))
        elif total == 1:
            story.append(_p(
                f"<b>Index {index:.0f} of 100 &mdash; a single-marker result.</b> This panel declares "
                "one feature, so the figure is one organism's position against the reference "
                "group, not a signature. It cannot be more than that however far from 50 it sits.",
                st["small"],
            ))
        elif bounds[1] - bounds[0] < 0.5:
            story.append(_p(
                f"<b>Index {index:.0f} of 100</b> over all {total} declared features, each measured as "
                "a point, so the coverage range collapses onto the index itself. That says every "
                "declared organism was measured; it says nothing about sampling error, which is not "
                "estimated here.",
                st["small"],
            ))
        else:
            story.append(_p(
                f"<b>Index {index:.0f} of 100</b> over all {total} declared features. The pale band is "
                f"the censoring range {bounds[0]:.0f}&ndash;{bounds[1]:.0f}, which reflects how "
                "precisely each organism could be measured — not sampling error.",
                st["small"],
            ))
    else:
        story.append(_p(
            f"<b>No index &mdash; {_STATUS_WORDS.get(status, status.replace('_', ' '))}.</b> "
            "This is not a negative result and not clearance. It means the number could not be "
            "produced, for the reasons below.",
            st["body"],
        ))
        reasons = record.get("reason_codes") or []
        if reasons:
            story.append(_p(
                "Reasons recorded: " + "; ".join(str(r).replace("_", " ") for r in reasons) + ".",
                st["small"],
            ))
        for note in record.get("notes") or []:
            story.append(_p(str(note), st["small"]))

    # what was measured, and against what
    story.append(Spacer(1, 2 * mm))
    story.append(_p("What was measured, and against whom", st["h2"]))
    transport = str(record.get("assay_transport", ""))
    bits = [
        f"Source comparison: <b>{str(record.get('contrast') or '').replace('_', ' ')}</b>.",
        f"Measurement: {str(record.get('measurement') or '').replace('_', ' ')}.",
    ]
    if transport.startswith("unvalidated_16s"):
        bits.append(
            "<b>This is a 16S-derived pattern applied to shotgun measurements.</b> The source study "
            "named organisms with a different method from the one used here, so the translation is "
            "unvalidated and inherits none of the source's sensitivity or specificity."
        )
    if panel is not None and panel.percentile_note:
        bits.append(panel.percentile_note)
    story.append(_p(" ".join(bits), st["small"]))
    story.append(_p(
        "No percentile is offered on this scale: that would need a held-out reference score "
        "distribution on this same set of features, and none exists.",
        st["fine"],
    ))

    # every feature
    story.append(Spacer(1, 2 * mm))
    story.append(_p("Every declared feature, including the ones that could not be measured", st["h2"]))
    _feature_table(story, st, record)

    # opposing evidence
    counter = record.get("counterevidence_ids") or []
    conflicts = _CONFLICTS.get(_panel_id(record))
    if counter or conflicts:
        story.append(Spacer(1, 2 * mm))
        story.append(_p("Evidence pointing the other way", st["h2"]))
        if conflicts:
            story.append(_p(conflicts, st["small"]))
        if counter:
            story.append(_p(
                "Sources recorded as counterevidence for this panel: "
                + ", ".join(str(c) for c in counter) + ".",
                st["fine"],
            ))
    not_independent = record.get("not_independent_of") or []
    if not_independent:
        story.append(_p(
            "Shares participants with " + ", ".join(str(x) for x in not_independent)
            + ", so the two are one evidence group rather than independent agreement.",
            st["fine"],
        ))


#: Named, specific conflicts the specification requires be prominent.
_CONFLICTS: Mapping[str, str] = {
    "PSO_DENG2026_SGB4348_V1": (
        "<b>This marker is not specific to psoriasis.</b> The same SGB4348 was reported as increased "
        "in pemphigus foliaceus, a different blistering skin disease, in a separate cohort. A high "
        "reading on a single shared marker therefore does not point at psoriasis in particular. A "
        "further study of 53 people with psoriasis and 47 controls found no broad diversity "
        "difference at all, with inflammation confounding the apparent profiles."
    ),
    "PSO_CHANG2022_TAX_V1": (
        "<b>Two of these studies disagree on direction.</b> <i>Bacteroides coprocola</i> is recorded "
        "here as increased in psoriasis, following its source, but a separate psoriasis cohort "
        "reported that organism as enriched in <i>controls</i>. The slot is kept as its own source "
        "described it and the disagreement is not resolved by preferring one paper."
    ),
    "PSORIATIC_SHARED_XIAO2024_TAX_V1": (
        "<b>This is a pooled contrast.</b> Skin-only psoriasis and psoriatic arthritis were combined "
        "against healthy controls, and no species difference separated the two from each other in "
        "that cohort. A high reading here cannot be read as psoriatic arthritis rather than "
        "psoriasis, and the source's q-values belong to the pooled comparison only."
    ),
    "HS_OGUT2022_TRANSLATED_V1": (
        "<b>A separate shotgun study found nothing.</b> Comparing hidradenitis suppurativa against "
        "controls, it reported no significant taxonomic difference and no significant diversity "
        "difference. With one translated feature here against an explicit negative comparison there, "
        "the honest summary is that stool evidence for this condition is very weak."
    ),
    "VIT_ACTIVE_JU2025_TAX_V1": (
        "<b>Cases and controls came from different studies.</b> The ten cases and twenty controls were "
        "deposited by separate projects, so sequencing batch and disease cannot be separated. This "
        "panel is also not a replication of the other vitiligo panel: they share almost no features."
    ),
    "PSO_SA2026_TAX_V1": (
        "<b>Treatment is not separable from disease here.</b> Eighteen of the twenty-four psoriasis "
        "participants were on immunosuppression and five had parasitic infection, against ten "
        "controls. Differences attributed to psoriasis may belong to its treatment."
    ),
}


def _sebd_page(story: list[Any], st: dict[str, ParagraphStyle], *, section: int, n: int, task: Mapping[str, Any]) -> None:
    """Seborrheic dermatitis: the real evidence, and why it is not a score."""
    story.extend(section_heading(section, "Seborrheic dermatitis — the evidence, and why there is no score", st["h1"], sub=n))
    story.append(_p(f"<b>{task.get('statement', '')}</b>", st["body"]))
    story.append(_p(
        "Seborrheic dermatitis is a scaly, greasy rash of the scalp, face and upper trunk; dandruff "
        "is its mildest form. Gut involvement has been reported, and the three blocks below are what "
        "that reporting actually consists of. None of them defines a stool-DNA signature, and "
        "publishing a number anyway would be inventing one.",
        st["small"],
    ))
    story.append(Spacer(1, 2 * mm))

    blocks = (
        (
            "Stool bacterial culture",
            "A 2019 series enrolled 67 people aged 18 to 57 and describes 30 controls, using stool "
            "<b>bacteriology</b> rather than sequencing. It reported how often counts were low: "
            "bifidobacteria in 59.4%, enterococci 62.5%, <i>E. coli</i> 65.7%, lactobacilli 90.6%, "
            "bacteroides 71.9%. Those percentages imply a denominator near 32, which is not explained "
            "against the 67 people enrolled, and the denominator discrepancy is preserved here rather "
            "than tidied away. Per-subject values, control distributions and usable thresholds are "
            "absent. Counts of living organisms grown on a plate are a different quantity from the "
            "share of DNA in a sequencing run, so these cannot be converted into abundance "
            "coefficients — and a study reporting <i>low</i> E. coli does not become a report of high "
            "E. coli because a generic dysbiosis template expects it.",
        ),
        (
            "Genetic hypotheses",
            "A 2024 study used genetic instruments and national disease records rather than measuring "
            "stool from people with the condition. Its nominal associations point in various "
            "directions across broad, overlapping taxonomic ranks. These concern genetically "
            "instrumented exposures, not an observed diagnostic stool profile, so no organism here "
            "receives a weight. The printed significance threshold in that paper is also internally "
            "inconsistent. Human genetic risk is not recovered from stool sequencing in this report.",
        ),
        (
            "Scalp and skin findings",
            "Studies of <i>Malassezia</i>, of Staphylococcus-to-Cutibacterium ratios, of scalp lipases "
            "and of a topical probiotic all measure <b>skin</b>, not stool. They cannot seed a gut "
            "classifier. Fungal organisms seen in stool remain exactly that: measured stool fungal "
            "observations, which neither demonstrate scalp overgrowth nor prove contamination.",
        ),
    )
    for title, text in blocks:
        story.append(_p(title, st["h2"]))
        story.append(_p(text, st["small"]))
        story.append(Spacer(1, 1.2 * mm))

    story.append(Spacer(1, 1 * mm))
    story.append(_p("What would have to exist for this to become a score", st["h2"]))
    story.append(_p(str(task.get("promotion_requirement", "")), st["small"]))
    story.append(_p(
        "A reported figure such as \"60% had dysbiosis\" in a clinical series is not a screening "
        "prevalence, a disease prior, or a cutoff, and is not used as one.",
        st["fine"],
    ))


def skin_detail_pages(
    story: list[Any],
    st: dict[str, ParagraphStyle],
    *,
    section: int,
    skin: Mapping[str, Any] | None,
    plan: Any | None = None,
) -> None:
    """One page per panel, then the tasks that produce no number."""
    story.extend(section_heading(section, "Gut-skin axis, in detail", st["h1"]))
    records = _panels(skin)
    if not records:
        story.append(_p("Skin research panels were not computed for this run.", st["body"]))
        return
    story.append(_p(
        "Each panel below gets its own page: what the source study compared, every organism or gene "
        "family it declared, what this sample measured for each, and the evidence that argues "
        "against reading too much into it. Slots that could not be measured are listed as such, "
        "because leaving them out would make a partly-measured panel look complete.",
        st["body"],
    ))
    n = 0
    for record in records:
        n += 1
        if n > 1:
            story.append(PageBreak())
        else:
            # The first panel shares the page with the paragraph that
            # introduces it: a lone intro followed by a forced break is a
            # page with four lines on it.
            story.append(Spacer(1, 3 * mm))
        _panel_detail(story, st, record, section=section, n=n)
        if plan is not None:
            from openbiota.pdfatlas import evidence_cards

            family = str(record.get("profile_family") or "")
            cards = plan.cards_for("profile", family)
            if cards:
                evidence_cards(
                    story, st, cards,
                    heading="What the human research has actually tested for this condition",
                )
            else:
                story.append(Spacer(1, 2 * mm))
                story.append(_p("Matched intervention evidence", st["h2"]))
                story.append(_p(
                    "No supported signature-matched intervention evidence was identified for this "
                    "pattern. Resembling a research pattern is not a diagnosis and does not establish "
                    "treatment eligibility, predicted benefit, or any need to remove organisms. "
                    "Nothing on this page triggers antibiotics, antifungals, botanical antimicrobials, "
                    "a prescription change, or a stool transplant.",
                    st["small"],
                ))

    tasks = (skin or {}).get("unavailable_tasks") or {}
    sebd = tasks.get("SEBD_STOOL_RESEARCH_V1")
    if sebd:
        n += 1
        story.append(PageBreak())
        _sebd_page(story, st, section=section, n=n, task=sebd)
    for task_id, task in tasks.items():
        if task_id == "SEBD_STOOL_RESEARCH_V1":
            continue
        n += 1
        story.append(Spacer(1, 3 * mm))
        story.extend(section_heading(section, f"{task.get('label') or task_id.replace('_', ' ')}", st["h2"], sub=n))
        story.append(_p(str(task.get("statement", "")), st["small"]))
        story.append(_p(
            "<b>What would have to exist:</b> " + str(task.get("promotion_requirement", "")),
            st["fine"],
        ))

    observations = (skin or {}).get("non_scored_observations") or []
    if observations:
        story.append(Spacer(1, 3 * mm))
        story.append(_p("Measured alongside, with no disease-score weight", st["h2"]))
        story.append(_p(
            "These were reported by the source studies and are shown because they were measured, not "
            "because they contribute to any score. DNA capacity for a pathway is not a measurement of "
            "a metabolite, and an annotation is not a validated gene call.",
            st["small"],
        ))
        rows = [[
            _p("SOURCE FINDING", st["label"]), _p("WHAT IT MAY BE READ AS", st["label"]),
        ]]
        for obs in observations:
            rows.append([
                _p(f"<b>{obs.get('feature', '')}</b><br/>"
                   f"<font size=6.5 color='{_hex(INK_FAINT)}'>{obs.get('source_direction', '')}</font>",
                   st["cell"]),
                _p(str(obs.get("permitted_interpretation", "")), st["cell"]),
            ])
        table = Table(rows, colWidths=[CONTENT_WIDTH * 0.34, CONTENT_WIDTH * 0.66], hAlign="LEFT", repeatRows=1)
        table.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("BACKGROUND", (0, 0), (-1, 0), PANEL_BG),
            ("LINEBELOW", (0, 0), (-1, 0), 0.5, RULE),
            ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        story.append(table)


def axis_panel(
    story: list[Any], st: dict[str, ParagraphStyle], axis: Any | None
) -> None:
    """The one headline number for this section, drawn as the Myco-Score is.

    The panels below each answer one study's question; this is their
    coverage-weighted middle, turned the healthy way up so it reads like
    every other favourable number in the report.
    """
    from openbiota.pdfreport import CONTENT_WIDTH, PercentileBar, status_for

    if axis is None or not getattr(axis, "scored", False):
        return
    value = float(axis.value)
    colour = status_for(percentile=value, higher_means="favourable").colour
    row = Table([[
        _p("<b>Gut-skin axis</b> <font color='#8C99A6'>(beta)</font>", st["h2"]),
        PercentileBar(width=CONTENT_WIDTH * 0.40, percentile=value,
                      higher_means="favourable", marker_colour=colour,
                      height=7.0 * mm, show_scale=True),
        _p(f"<font size='20' color='{_hex(colour)}'><b>{int(value + 0.5)}</b></font>"
           f"<font size='8' color='#8C99A6'>/100</font>", st["h2"]),
    ]], colWidths=[CONTENT_WIDTH * 0.30, CONTENT_WIDTH * 0.44, CONTENT_WIDTH * 0.26],
        hAlign="LEFT")
    row.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    story.append(row)
    top = ", ".join(
        f"{c['label'].split('—')[0].strip()} ({c['index']:.0f})"
        for c in axis.contributors[:3]
    )
    story.append(_p(
        f"<font size='6.6' color='#6B7785'>Higher is the better direction. The weighted "
        f"middle of {axis.n_scored} published panels, each weighted by how much of its own "
        f"marker list this sample could measure. Carried mostly by {top}.</font>",
        st["small"],
    ))
