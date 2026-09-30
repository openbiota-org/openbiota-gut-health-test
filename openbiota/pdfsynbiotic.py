"""Rendering for the A16 synbiotic scenarios (BUILD_SPEC_v0.8.3 §11).

Two blocks. `simulation_block` sits in the actions section: the four arms of
every candidate, the interaction contrast, who each fibre would feed -
including the organisms this report has already flagged - and who would make
the butyrate. `simulation_detail_block` sits in the metabolism detail: the
full member tables behind those summaries.

Every number is the optimum of a linear programme and the block says how it
was reproduced. Nothing here is a measurement, a concentration or a dose,
and the text says so in every place a reader could take it for one.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from functools import lru_cache
from typing import Any, Final

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, Flowable, Spacer, Table, TableStyle

from openbiota.extension.communitytype import UNRESOLVED
from openbiota.pdflinks import (
    Anchor,
    Paragraph,
    back_link_line,
    section_dest,
)
from openbiota.pdfparts import DASH, num
from openbiota.pdfparts import hex_of as _hex
from openbiota.pdfparts import pct as _pct
from openbiota.pdfparts import table as _table
from openbiota.pdfreport import (
    AMBER,
    CONTENT_WIDTH,
    CORAL,
    GREEN,
    INK_FAINT,
    INK_SOFT,
    RULE,
    SLATE,
)
from openbiota.pdfsummary import GradedBar

#: verdict -> (words a reader sees, tone)
VERDICT_WORDS: Final[Mapping[str, tuple[str, str]]] = {
    "pair_gains": ("pair gains", "green"),
    "pair_additive": ("pair adds up", "green"),
    "substrate_suffices": ("fibre alone suffices", "slate"),
    "probiotic_suffices": ("consortium alone suffices", "slate"),
    "no_gain": ("no gain", "slate"),
    "pair_worse": ("parts interfere", "amber"),
    "unavailable": ("not solved", "amber"),
}

ARM_WORDS: Final[Mapping[str, str]] = {
    "neither": "as you are",
    "probiotic_only": "+ consortium",
    "substrate_only": "+ fibre",
    "both": "+ both",
}

TONES: Final[Mapping[str, colors.Color]] = {"green": GREEN, "amber": AMBER, "slate": SLATE}


def species_label(name: str) -> str:
    """A model member id rendered as a species name a reader can read.

    `Clostridium_sp_7_2_43FAA` is a real AGORA2 id; printed raw it reads as a
    typo, and title-cased it reads as a different organism. The genus keeps
    its capital, the epithet keeps its case, and an unnamed-species tail is
    joined back up rather than split into words.
    """
    parts = [p for p in str(name).replace("_", " ").split() if p]
    if not parts:
        return str(name)
    genus = parts[0][:1].upper() + parts[0][1:]
    rest = parts[1:]
    if rest and rest[0].lower() in {"sp", "sp."}:
        tail = "_".join(rest[1:])
        return f"{genus} sp. {tail}" if tail else f"{genus} sp."
    return " ".join([genus, *rest])


def _short(display: str) -> str:
    """*Genus species* -> *G. species* for a dense cell; other forms untouched."""
    parts = species_label(display).split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return parts[0] if parts else display


def _italic(display: str) -> str:
    return f"<i>{_short(display)}</i>"


def _verdict_cell(verdict: str, st: Mapping[str, ParagraphStyle]) -> Paragraph:
    label, tone = VERDICT_WORDS.get(verdict, (verdict, "slate"))
    return Paragraph(f"<font color='{_hex(TONES[tone])}'><b>{label}</b></font>", st["cell"])


def _baseline_candidates(scenarios: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    baseline = scenarios.get("baseline_medium")
    cands = [c for c in scenarios.get("candidates") or [] if c.get("medium_id") == baseline]
    ranking = (scenarios.get("ranking") or {}).get("by_raw_butyrate_flux") or {}
    order = {r["candidate_id"]: (r["rank"] is None, r["rank"] or 0)
             for r in ranking.get(baseline) or []}
    cands.sort(key=lambda c: order.get(c["candidate_id"], (True, 0)))
    return cands


# --------------------------------------------------------------------------- #
# the actions-section block
# --------------------------------------------------------------------------- #


def simulation_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """What the metabolic model predicts for each fibre with the consortium."""
    if not view:
        return
    readiness = view.get("readiness") or {}
    protocol = view.get("protocol") or {}
    scenarios = view.get("scenarios") or {}

    story.append(CondPageBreak(90 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Personalised synbiotic options: what the model predicts", st["h2"]))
    story.append(Paragraph(
        "Your own community is rebuilt as a metabolic model from the organisms found in your "
        "sample and run four ways for each fibre \u2014 as it is, with a five-species probiotic "
        "consortium added, with the fibre added, and with both \u2014 so the question \u201cis the "
        "combination worth more than either part alone?\u201d can be asked directly, and the "
        "question \u201cwho would this fibre actually feed in <i>my</i> gut?\u201d can be answered "
        "organism by organism. These are <b>experimental predictions from a model</b>, not "
        "measurements, and not a guarantee of how you would respond.", st["body"],
    ))
    _status_strip(story, st, readiness, protocol, scenarios)
    _consortium_note(story, st, view)

    state = scenarios.get("execution_state")
    if not scenarios or state in (None, "not_run"):
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<font size='6.8' color='{_hex(INK_SOFT)}'><b>Scenarios not solved in this run.</b> "
            + str(scenarios.get("execution_state_meaning")
                  or "The scenario solver was not started for this sample.")
            + "</font>", st["small"],
        ))
        _replay_note(story, st, view)
        return

    _community_note(story, st, scenarios)
    _four_arm_table(story, st, scenarios)
    _feeding_table(story, st, scenarios)
    _makers_note(story, st, scenarios)
    _coverage_block(story, st, view)
    _stability_note(story, st, scenarios)
    _replay_note(story, st, view)
    _evidence_seeds_block(story, st, view)
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>Units are mmol per gram of community dry "
        "weight per hour in the model: a production potential, not a stool concentration, an "
        "absorbed amount or a dose. Fibre caps in the model are exchange limits, not grams to "
        "eat. The consortium is a species-level perturbation following the published protocol, "
        "not a named product. Organisms the model library lacks are counted as unmapped, never "
        "substituted. Feeding of a flagged organism is a specific model edge \u2014 it carries the "
        "transporter and the community state permits the uptake \u2014 not a rule that fibre feeds "
        "pathogens. In ISAPP terms every row is a <i>candidate combination</i>; a positive "
        "contrast in the model is evidence about the model, not a demonstrated host benefit."
        "</font>", st["small"],
    ))


#: Coverage runs the good way up, so the bands are the concern bands reversed:
#: a high share of goals covered is the green end.
_COVERAGE_BANDS: Final[tuple[tuple[float, colors.Color, colors.Color], ...]] = (
    (25.0, colors.HexColor("#F2CFCB"), CORAL),
    (50.0, colors.HexColor("#F0C4A8"), colors.HexColor("#D0641E")),
    (75.0, colors.HexColor("#F6DDB4"), AMBER),
    (100.0, colors.HexColor("#CFE4CF"), GREEN),
)


def _coverage_marker(value: float) -> colors.Color:
    for upper, _track, marker in _COVERAGE_BANDS:
        if value < upper or upper >= 100.0:
            return marker
    return GREEN


def _coverage_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """§11.4's C and P, each on its own scale.

    Two different questions, so two bars rather than one blended number.
    C asks how many of the modelled goals have a favourable edge at all; P
    asks the stricter question of whether the *combination* is what carries
    the endpoint, rather than the fibre or the consortium alone. P can be
    much the lower of the two without anything being wrong: it usually means
    the fibre works on its own.
    """
    scores = (view or {}).get("coverage") or {}
    if not scores:
        return
    story.append(CondPageBreak(30 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("How much of the goal set these options cover", st["h3"]))

    rows: list[list[Any]] = []
    for key, label, gloss in (
        ("goal_coverage", "Goal coverage (C)",
         "share of modelled goals with a favourable edge"),
        ("pair_supported_goal_coverage", "Pair-supported coverage (P)",
         "share where the pair, not either part alone, carries it"),
    ):
        value = scores.get(key)
        if value is None:
            reason = (
                scores.get("pair_coverage_unavailable_reason")
                or scores.get("unavailable_reason")
                or "not applicable"
            )
            rows.append([
                Paragraph(f"<b>{label}</b>", st["cell"]),
                Paragraph(f"<font color='{_hex(INK_FAINT)}'>{DASH}</font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{reason}</font>",
                          st["cell"]),
            ])
            continue
        number = float(value)
        rows.append([
            Paragraph(f"<b>{label}</b>", st["cell"]),
            Paragraph(f"<b><font color='{_hex(_coverage_marker(number))}'>"
                      f"{number:.0f}%</font></b>", st["cell"]),
            GradedBar(width=CONTENT_WIDTH * 0.40, value=number,
                      bands=_COVERAGE_BANDS, marker_colour=_coverage_marker(number)),
        ])
        rows.append([
            Paragraph("", st["cell"]),
            Paragraph("", st["cell"]),
            Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{gloss}</font>", st["small"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.30, 0.10, 0.60)]))
    n_goals = scores.get("n_goals")
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>Over {n_goals} modelled "
        f"goal{'' if n_goals == 1 else 's'}, each weighted equally. A candidate whose model did "
        "not solve is left out of the set rather than counted against it. These are percentages "
        "of a modelled goal set \u2014 not a response rate, and not a share of people helped."
        "</font>", st["small"],
    ))


def _evidence_seeds_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """§24's published synbiotic records, counted honestly.

    The negatives and the missing arms are the point of showing this at all:
    a list of only the encouraging trials would misrepresent the literature
    this feature is built on.
    """
    seeds = (view or {}).get("evidence_seeds") or {}
    total = seeds.get("n_records")
    if not total:
        return
    story.append(CondPageBreak(34 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("What the published synbiotic trials actually show", st["h3"]))

    counts = (
        ("demonstrated synergy", seeds.get("n_demonstrated_synergy"), GREEN),
        ("no synergy demonstrated", seeds.get("n_no_synergy_demonstrated"), CORAL),
        ("undetermined", seeds.get("n_undetermined"), SLATE),
        ("missing a component-only arm", seeds.get("n_missing_a_component_arm"), AMBER),
    )
    rows: list[list[Any]] = [[
        Paragraph(f"<b><font color='{_hex(tone)}'>{value}</font></b>"
                  f"<font size='6.2' color='{_hex(INK_SOFT)}'> of {total}</font><br/>"
                  f"<font size='6.2'>{label}</font>", st["cell"])
        for label, value, tone in counts
    ]]
    story.append(_table(rows, [CONTENT_WIDTH * 0.25] * 4))
    note = seeds.get("note")
    if note:
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_SOFT)}'>{note}</font>", st["small"],
        ))


def _status_strip(story: list[Any], st: Mapping[str, ParagraphStyle], readiness: Mapping[str, Any],
                  protocol: Mapping[str, Any], scenarios: Mapping[str, Any]) -> None:
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("COMPONENT", "STATE", "DETAIL")]]
    rows.append([
        Paragraph("<b>Simulation engine</b>", st["cell"]),
        Paragraph(f"<font color='{_hex(GREEN)}'>ready</font>"
                  if readiness.get("micom", {}).get("installed") else "not installed", st["cell"]),
        Paragraph(f"<font size='6.4'>MICOM {protocol.get('micom_version', '')}; community growth "
                  f"held at {protocol.get('tradeoff')} of its maximum while each output is bounded"
                  "</font>", st["cell"]),
    ])
    rows.append([
        Paragraph("<b>Metabolic model library</b>", st["cell"]),
        Paragraph(f"<font color='{_hex(GREEN)}'>ready</font>"
                  if readiness.get("model_library", {}).get("installed") else "not installed",
                  st["cell"]),
        Paragraph("<font size='6.4'>AGORA2 species models, checksum verified</font>", st["cell"]),
    ])
    repro = scenarios.get("reproducibility") or {}
    agree = repro.get("all_arms_agree")
    if agree:
        state_cell = f"<font color='{_hex(GREEN)}'>exact</font>"
        words = "agreed to every digit"
    elif agree is False:
        state_cell = f"<font color='{_hex(AMBER)}'>disagreement</font>"
        words = "did not agree \u2014 that arm is withheld"
    else:
        state_cell = "ready"
        words = "are compared"
    solved = repro.get("linear_programmes_solved")
    cached = repro.get("linear_programmes_from_cache")
    counts = (f" ({solved:,} programmes solved, {cached:,} served from cache)"
              if isinstance(solved, int) and isinstance(cached, int) else "")
    rows.append([
        Paragraph("<b>Numerical solver</b>", st["cell"]),
        Paragraph(state_cell, st["cell"]),
        Paragraph(
            "<font size='6.4'>every value is the optimum of a linear programme; each arm's growth "
            f"maximum was recomputed by independent worker processes and {words}{counts}</font>",
            st["cell"],
        ),
    ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.26, 0.14, 0.60)]))


def _consortium_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """Which five organisms the model adds.

    The text above says \u201ca five-species consortium\u201d, which tells the
    reader nothing they could check. These are species added to a model, not
    a product on a shelf, and naming them is what makes that difference
    visible.
    """
    consortium = [str(x) for x in ((view or {}).get("consortium") or [])]
    if not consortium:
        return
    names = ", ".join(f"<i>{n}</i>" for n in consortium[:-1])
    names = f"{names} and <i>{consortium[-1]}</i>" if len(consortium) > 1 else (
        f"<i>{consortium[0]}</i>"
    )
    story.append(Paragraph(
        f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>The consortium added in the model "
        f"is {names}.</b> These are species added to a metabolic model following the "
        "published protocol \u2014 not a named product, not a dose, and not a recommendation "
        "to take any of them.</font>", st["small"],
    ))


def _community_note(story: list[Any], st: Mapping[str, ParagraphStyle],
                    scenarios: Mapping[str, Any]) -> None:
    community = scenarios.get("community") or {}
    coverage = community.get("coverage") or {}
    n_members = community.get("n_members")
    mapped = coverage.get("mapped_fraction_of_input")
    n_input = coverage.get("n_input_taxa")
    aliases = community.get("aliases_resolved") or {}
    limits = community.get("growth_limits") or {}
    parts: list[str] = []
    if n_members is not None and n_input:
        parts.append(
            f"<b>{n_members} of the {n_input} organisms found in your sample have a metabolic "
            "model in the library</b>"
            + (f", covering {mapped * 100:.0f}% of the bacterial abundance they represent"
               if isinstance(mapped, (int, float)) else "") + "."
        )
    if aliases:
        shown = "; ".join(
            f"{_italic(sample)} under its earlier name {_italic(model)}"
            for model, sample in list(aliases.items())[:3]
        )
        parts.append(
            "The model library predates some renames, so two of your organisms were matched "
            f"{shown}. They are the same organism, not an extra one."
        )
    if limits:
        nutrient = next(iter(limits))
        parts.append(
            f"In this diet definition the model's total growth is capped by "
            f"<b>{_nutrient_word(nutrient)}</b> supply, so growth is identical in every arm and "
            "only what the community <i>makes</i> is compared."
        )
    if parts:
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(" ".join(parts), st["small"]))


NUTRIENT_WORDS: Final[Mapping[str, str]] = {
    "cobalt2": "cobalt", "fe2": "iron", "fe3": "iron", "zn2": "zinc", "mn2": "manganese",
    "cu2": "copper", "mg2": "magnesium", "ca2": "calcium", "k": "potassium", "pi": "phosphate",
    "so4": "sulfate", "nh4": "ammonium", "cl": "chloride", "na1": "sodium",
}


def _nutrient_word(exchange: str) -> str:
    key = exchange.removeprefix("EX_").removesuffix("_m")
    return NUTRIENT_WORDS.get(key, key.replace("_", " "))


def _four_arm_table(story: list[Any], st: Mapping[str, ParagraphStyle],
                    scenarios: Mapping[str, Any]) -> None:
    """One row per fibre: the four arms, the interaction and the verdict."""
    cands = _baseline_candidates(scenarios)
    if not cands:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>Butyrate potential under each option</b> "
                           f"<font size='6.4' color='{_hex(INK_SOFT)}'>mmol/gDW/h in the model; "
                           "ranked by the combined option</font>", st["h3"]))
    head = ("FIBRE", "AS YOU ARE", "+ CONSORTIUM", "+ FIBRE", "+ BOTH", "INTERACTION", "VERDICT")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    proxies: list[str] = []
    for cand in cands:
        sub = cand.get("substrate") or {}
        arms = (cand.get("contrasts") or {}).get("arms") or {}
        interaction = (cand.get("contrasts") or {}).get("interaction")
        label = str(sub.get("label") or cand.get("candidate_id"))
        mark = ""
        if sub.get("is_model_proxy"):
            proxies.append(f"{label}: {sub.get('proxy_note') or 'a model stand-in for the food'}")
            mark = "*"
        food = sub.get("food_context")
        rows.append([
            Paragraph(f"<b>{label}{mark}</b>"
                      + (f"<br/><font size='6.2' color='{_hex(INK_SOFT)}'>{food}</font>" if food else ""),
                      st["cell"]),
            Paragraph(num(arms.get("neither"), 1), st["cell"]),
            Paragraph(num(arms.get("probiotic_only"), 1), st["cell"]),
            Paragraph(num(arms.get("substrate_only"), 1), st["cell"]),
            Paragraph(f"<b>{num(arms.get('both'), 1)}</b>", st["cell"]),
            Paragraph(("+" if (interaction or 0) > 0 else "") + num(interaction, 1), st["cell"]),
            _verdict_cell(str(cand.get("verdict")), st),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.11, 0.13, 0.11, 0.11, 0.12, 0.20)]))
    ranking = scenarios.get("ranking") or {}
    first_raw, first_norm = ranking.get("first_by_raw_flux"), ranking.get("first_by_flux_per_growth")
    note = (
        "<b>Interaction</b> is the combined option minus both single options plus the baseline "
        "(F11 \u2212 F10 \u2212 F01 + F00): above zero, the pair does more than the sum of its parts "
        "in the model; below zero, less. A fibre that adds nothing on its own is still shown, and "
        "a verdict of \u201cfibre alone suffices\u201d means the consortium adds little on top of it."
    )
    if first_raw and first_norm and first_raw != first_norm:
        note += (f" Ranked by output per unit of community growth instead, <b>{first_norm}</b> "
                 f"comes first rather than <b>{first_raw}</b>; the two objectives are shown "
                 "separately in the data file and neither is a recommendation on its own.")
    def _fibre_alone_is_flat(cand: Mapping[str, Any]) -> bool:
        arms = (cand.get("contrasts") or {}).get("arms") or {}
        alone, base = arms.get("substrate_only"), arms.get("neither")
        return alone is not None and base is not None and abs(alone - base) < 1e-6

    flat = [c for c in cands if _fibre_alone_is_flat(c)]
    if flat and len(flat) == len(cands):
        note += (
            " <b>Every fibre alone leaves butyrate unchanged in your community</b>, while the "
            "same fibre with the consortium raises it: your present butyrate producers are "
            "already at capacity in the model, so extra fibre has nowhere to go until producers "
            "are added. That is the interaction these rows are measuring, and it is the case "
            "where a combination is worth more than either part."
        )
    if proxies:
        note += " * " + " ".join(proxies)
    story.append(Paragraph(f"<font size='6.4' color='{_hex(INK_SOFT)}'>{note}</font>", st["small"]))


def _feeding_table(story: list[Any], st: Mapping[str, ParagraphStyle],
                   scenarios: Mapping[str, Any]) -> None:
    """Who each fibre would feed in this community, with the report's own flags.

    The question a reader has about a fibre is not only "does it raise
    butyrate?" but "does it feed the organism you told me is overgrown?".
    Both are answered from the same arms: each member's possible share of the
    fibre in the state where the community makes the most butyrate it can.
    """
    cands = _baseline_candidates(scenarios)
    rows_out: list[list[Any]] = []
    for cand in cands:
        eaters = (cand.get("eaters") or {}).get("both") or (cand.get("eaters") or {}).get("substrate_only") or []
        if not eaters and cand.get("solved"):
            eaters = []
        possible = [e for e in eaters if e.get("possible")]
        flagged = [e for e in possible if e.get("flagged")]
        beneficial = [e for e in possible if e.get("beneficial") or e.get("added_by_consortium")]
        other = [e for e in possible if e not in flagged and e not in beneficial]
        label = str((cand.get("substrate") or {}).get("label") or cand.get("candidate_id"))

        def _cell(items: Sequence[Mapping[str, Any]], *, tone: colors.Color | None) -> str:
            if not items:
                return f"<font color='{_hex(INK_FAINT)}'>none in the model</font>"
            shown = items[:5]
            text = "; ".join(
                f"{_italic(str(e.get('display')))} \u2264{_pct(e.get('can_share'))}"
                + (" <b>required</b>" if e.get("required") else "")
                + (f" <font size='5.8'>(p{e['percentile']:.0f})</font>"
                   if tone is AMBER and isinstance(e.get("percentile"), (int, float)) else "")
                for e in shown
            )
            if len(items) > len(shown):
                text += f"; +{len(items) - len(shown)} more"
            return f"<font color='{_hex(tone)}'>{text}</font>" if tone else text

        if not possible:
            summary = (f"<font color='{_hex(INK_FAINT)}'>no organism in your modelled community "
                       "carries a transporter for it</font>" if cand.get("solved")
                       else f"<font color='{_hex(AMBER)}'>not solved</font>")
            rows_out.append([Paragraph(f"<b>{label}</b>", st["cell"]),
                             Paragraph(summary, st["cell"]), Paragraph(DASH, st["cell"]),
                             Paragraph(DASH, st["cell"])])
            continue
        rows_out.append([
            Paragraph(f"<b>{label}</b><br/><font size='6.2' color='{_hex(INK_SOFT)}'>"
                      f"{len(possible)} of your organisms can use it</font>", st["cell"]),
            Paragraph(_cell(beneficial, tone=GREEN), st["cell"]),
            Paragraph(_cell(flagged, tone=AMBER), st["cell"]),
            Paragraph(_cell(other, tone=None), st["cell"]),
        ])
    if not rows_out:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>Who each fibre would feed in your community</b> "
                           f"<font size='6.4' color='{_hex(INK_SOFT)}'>share of the fibre each "
                           "organism can take up when the community makes the most butyrate it "
                           "can</font>", st["h3"]))
    head = ("FIBRE", "BENEFICIAL OR ADDED ORGANISMS", "ORGANISMS THIS REPORT FLAGGED", "OTHER ORGANISMS")
    rows = [[Paragraph(h, st["label"]) for h in head], *rows_out]
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.16, 0.30, 0.30, 0.24)]))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'>\u201c\u2264 40%\u201d means that in some "
        "butyrate-optimal state of your community this organism takes up to 40% of the fibre "
        "supplied; <b>required</b> means no such state exists without it. A flagged organism "
        "listed here has the transporter and the community state permits the uptake \u2014 a "
        "specific, testable hypothesis about your sample, which is why it is shown rather than "
        "assumed either way.</font>", st["small"],
    ))
    _no_route_note(story, st, cands)


def _no_route_note(story: list[Any], st: Mapping[str, ParagraphStyle],
                   cands: Sequence[Mapping[str, Any]]) -> None:
    """Which flagged organisms a fibre cannot feed, by name.

    A table of who a fibre *can* feed leaves the more important question
    unanswered: the organism a reader is worried about might simply be
    missing from the row. An organism that carries no transporter for the
    substrate cannot take it up in the model at all, and saying so by name is
    the strongest statement available.
    """
    without: dict[str, set[str]] = {}
    for cand in cands:
        label = str((cand.get("substrate") or {}).get("label") or "")
        for organism in cand.get("flagged_without_route") or []:
            without.setdefault(str(organism.get("display") or organism.get("taxon")), set()).add(label)
    if not without:
        return
    all_labels = {str((c.get("substrate") or {}).get("label") or "") for c in cands}
    universal = sorted(name for name, labels in without.items() if labels >= all_labels)
    partial = sorted(
        (name, sorted(all_labels - labels)) for name, labels in without.items()
        if labels and not labels >= all_labels
    )
    parts: list[str] = []
    if universal:
        shown = ", ".join(_italic(n) for n in universal[:8])
        more = f" and {len(universal) - 8} other flagged organisms" if len(universal) > 8 else ""
        parts.append(
            f"<b>None of these fibres can feed {shown}{more}</b> \u2014 they carry no transporter "
            "for any of them in the model, so no amount of any of these fibres reaches them "
            "directly."
        )
    for name, can_eat in partial[:4]:
        parts.append(f"{_italic(name)} can use only {', '.join(can_eat)} of the fibres tested.")
    if not parts:
        return
    story.append(Spacer(1, 1.2 * mm))
    story.append(Paragraph(
        f"<font size='6.6'>{' '.join(parts)} This concerns direct uptake only: an organism can "
        "still gain from products another organism releases, which the model does allow and "
        "which the detail section shows.</font>", st["small"],
    ))


def _makers_note(story: list[Any], st: Mapping[str, ParagraphStyle],
                 scenarios: Mapping[str, Any]) -> None:
    """Who would make the butyrate: as you are, and under the best-ranked option."""
    cands = _baseline_candidates(scenarios)
    if not cands:
        return
    top = cands[0]
    makers_now = (top.get("makers") or {}).get("neither") or []
    makers_both = (top.get("makers") or {}).get("both") or []
    if not makers_now and not makers_both:
        return

    def _list(rows: Sequence[Mapping[str, Any]], limit: int = 7) -> str:
        required = [r for r in rows if r.get("required")]
        optional = [r for r in rows if not r.get("required") and r.get("possible")]
        parts: list[str] = []
        for r in required[:limit]:
            mark = " <font size='5.8'>(added)</font>" if r.get("added_by_consortium") else ""
            must, can = _pct(r.get("must_share")), _pct(r.get("can_share"))
            span = must if must == can else f"{must}\u2013{can}"
            parts.append(f"{_italic(str(r.get('display')))} {span}{mark}")
        text = "; ".join(parts) if parts else "no organism is required"
        if len(required) > limit:
            text += f"; +{len(required) - limit} more"
        if optional:
            text += (f" <font color='{_hex(INK_SOFT)}'>(optional: "
                     + ", ".join(_italic(str(r.get("display"))) for r in optional[:4])
                     + (f", +{len(optional) - 4} more" if len(optional) > 4 else "") + ")</font>")
        return text

    label = str((top.get("substrate") or {}).get("label") or top.get("candidate_id"))
    arms = (top.get("contrasts") or {}).get("arms") or {}
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>Who would make the butyrate</b> "
                           f"<font size='6.4' color='{_hex(INK_SOFT)}'>each organism's share of the "
                           "community's potential; a range is what it must and what it can "
                           "contribute</font>", st["h3"]))
    rows = [[Paragraph(h, st["label"]) for h in ("SCENARIO", "REQUIRED PRODUCERS")]]
    rows.append([
        Paragraph(f"<b>As you are</b><br/><font size='6.2' color='{_hex(INK_SOFT)}'>"
                  f"{num(arms.get('neither'), 1)} mmol/gDW/h</font>", st["cell"]),
        Paragraph(_list(makers_now), st["cell"]),
    ])
    rows.append([
        Paragraph(f"<b>{label} + consortium</b><br/><font size='6.2' color='{_hex(INK_SOFT)}'>"
                  f"{num(arms.get('both'), 1)} mmol/gDW/h</font>", st["cell"]),
        Paragraph(_list(makers_both), st["cell"]),
    ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.78)]))
    tight = [r for r in makers_now if r.get("required")
             and r.get("must_share") is not None and r.get("can_share") is not None
             and (r["can_share"] - r["must_share"]) < 0.03]
    if makers_now and len(tight) >= max(1, len([r for r in makers_now if r.get("required")]) - 1):
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>Your present producers' ranges are narrow "
            "\u2014 each must contribute almost exactly what it can \u2014 which means they are "
            "already at capacity in the model. That is why a fibre on its own may raise acetate "
            "and propionate without raising butyrate: the extra carbon has no spare "
            "butyrate-making capacity to flow into, until producers are added.</font>", st["small"],
        ))


def _arm_values(scenarios: Mapping[str, Any], medium_id: str, arm: str) -> list[float]:
    return [
        v for c in scenarios.get("candidates") or []
        if c.get("medium_id") == medium_id
        and (v := ((c.get("contrasts") or {}).get("arms") or {}).get(arm)) is not None
    ]


def _stability_note(story: list[Any], st: Mapping[str, ParagraphStyle],
                    scenarios: Mapping[str, Any]) -> None:
    """What the other declared diets do to the answer.

    A rank that moves is worth saying. But the interesting case in practice
    is not a reshuffle: it is a diet on which every candidate ties, because
    that means the background diet already supplies what the supplement
    would have added. Reporting only "unstable" would bury the finding.
    """
    stability = (scenarios.get("ranking") or {}).get("stability_across_media")
    media = {str(m.get("id")): str(m.get("label", m.get("id"))) for m in scenarios.get("media") or []}
    if not stability or len(media) < 2:
        return
    baseline = str(scenarios.get("baseline_medium"))
    parts: list[str] = []
    for medium_id, label in media.items():
        if medium_id == baseline:
            continue
        both = _arm_values(scenarios, medium_id, "both")
        probiotic = _arm_values(scenarios, medium_id, "probiotic_only")
        if not both or not probiotic:
            continue
        spread = max(both) - min(both)
        scale = max(abs(v) for v in both) or 1.0
        adds_nothing = spread <= 0.01 * scale and abs(max(both) - max(probiotic)) <= 0.01 * scale
        if adds_nothing:
            base_both = _arm_values(scenarios, baseline, "both")
            gain = (max(both) / max(base_both)) if base_both and max(base_both) > 0 else None
            parts.append(
                f"<b>On {label.lower()}, no fibre supplement changes the prediction.</b> Every "
                "candidate reaches the same output as the consortium alone, because that diet "
                "already supplies what the supplement would have added"
                + (f" \u2014 and it reaches {gain:.0%} of the best combined result on the "
                   "baseline diet without any supplement at all." if gain and gain > 1.01 else ".")
                + " If your own diet is already high in fibre, the model puts the benefit in the "
                "organisms rather than in the fibre."
            )
        else:
            moved = [sid for sid, ok in (stability.get("stable") or {}).items() if not ok]
            if moved:
                positions = stability.get("positions") or {}
                changes = "; ".join(
                    f"<b>{sid}</b> " + " \u2192 ".join(
                        f"{media.get(m, m)} #{r if r is not None else DASH}"
                        for m, r in (positions.get(sid) or {}).items()
                    ) for sid in moved[:3]
                )
                parts.append(f"<b>Ranking shifts under {label.lower()}:</b> {changes}.")
    if not parts and stability.get("all_stable"):
        others = ", ".join(v.lower() for k, v in media.items() if k != baseline)
        parts.append(f"<b>Ranking is stable across diets.</b> Rerun under {others}, every fibre "
                     "keeps its place.")
    if not parts:
        return
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font size='6.6'>{' '.join(parts)} This is a range across predeclared dietary "
        "assumptions, not a confidence interval and not a probability of response.</font>",
        st["small"],
    ))


def _replay_note(story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any]) -> None:
    replay = view.get("published_replay") or {}
    if not replay.get("passed"):
        return
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font size='6.6'><b>The method is validated against published predictions.</b> Run "
        "against the released predictions of the study it follows, this pipeline reproduces every "
        f"value to {replay.get('tolerance', 1e-8):g}, including the result that a fibre alone can "
        "outrank a fibre-plus-probiotic combination, and that ranking by total output and by "
        "output per unit of growth picks different winners. That replay and your solved model are "
        "different things, and the data file labels them differently.</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# the metabolism-detail block
# --------------------------------------------------------------------------- #


def simulation_detail_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """The mechanism behind each candidate: every arm's outputs and member tables.

    Placed in the metabolism detail so the actions section can stay a summary
    with the full accounting one section away.
    """
    scenarios = (view or {}).get("scenarios") or {}
    cands = _baseline_candidates(scenarios)
    if not cands or scenarios.get("execution_state") in (None, "not_run"):
        return
    story.append(CondPageBreak(70 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Synbiotic scenarios in detail: every arm, every output, every organism",
                           st["h2"]))
    story.append(Paragraph(
        "The full accounting behind the summary in the actions section. For each fibre: what "
        "the community could release of four fermentation products in each of the four arms, then "
        "which organisms would take up the fibre and which would make the butyrate. Each is a "
        "range over every state in which the community makes at least 99% of its butyrate "
        "maximum: the least an organism contributes in any such state, and the most.", st["body"],
    ))
    for cand in cands:
        _candidate_detail(story, st, cand)


def _candidate_detail(story: list[Any], st: Mapping[str, ParagraphStyle],
                      cand: Mapping[str, Any]) -> None:
    sub = cand.get("substrate") or {}
    label = str(sub.get("label") or cand.get("candidate_id"))
    arms = cand.get("arms") or {}
    story.append(CondPageBreak(55 * mm))
    story.append(Spacer(1, 2.5 * mm))
    verdict_label, tone = VERDICT_WORDS.get(str(cand.get("verdict")), (str(cand.get("verdict")), "slate"))
    story.append(Paragraph(
        f"<b>{label}</b> + consortium \u2014 <font color='{_hex(TONES[tone])}'>{verdict_label}</font>"
        f"<font size='6.4' color='{_hex(INK_SOFT)}'> · model cap {sub.get('cap_mmol_per_gDW_per_h')} "
        f"mmol/gDW/h on {sub.get('exchange_reaction')}</font>", st["h3"],
    ))
    # outputs per arm
    head = ("OUTPUT (mmol/gDW/h)", *[ARM_WORDS[a].upper() for a in ARM_WORDS])
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    for name in ("butyrate", "propionate", "acetate", "lactate"):
        cells = [Paragraph(f"<b>{name}</b>", st["cell"])]
        for arm_name in ARM_WORDS:
            arm = arms.get(arm_name) or {}
            value = arm.get("raw_flux") if name == "butyrate" else (arm.get("other_exchanges") or {}).get(name)
            if name == "butyrate" and arm.get("raw_flux") is None and arm.get("unavailable_reason"):
                cells.append(Paragraph(f"<font size='6' color='{_hex(AMBER)}'>unavailable</font>", st["cell"]))
            else:
                cells.append(Paragraph(num(value, 1), st["cell"]))
        rows.append(cells)
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.24, 0.19, 0.19, 0.19, 0.19)]))
    growth = (arms.get("neither") or {}).get("community_growth")
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_SOFT)}'>Each output is its own exchange bounded "
        "separately; they are never summed into one short-chain-fatty-acid number. Community "
        f"growth {num(growth, 4)} per hour in every arm (diet-capped).</font>", st["small"],
    ))
    # members
    for title, table_key, arm_name in (
        ("Organisms that would take up the fibre (+ both)", "eaters", "both"),
        ("Organisms that would make the butyrate (+ both)", "makers", "both"),
        ("Organisms that make the butyrate as you are", "makers", "neither"),
    ):
        members = (cand.get(table_key) or {}).get(arm_name) or []
        members = [m for m in members if m.get("possible") or m.get("required")]
        if not members:
            continue
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(f"<b>{title}</b>", st["small"]))
        rows = [[Paragraph(h, st["label"]) for h in ("ORGANISM", "MUST", "CAN", "STATUS IN THIS REPORT")]]
        for m in members[:14]:
            status = []
            if m.get("added_by_consortium"):
                status.append("added by the consortium")
            if m.get("flagged"):
                status.append(f"<font color='{_hex(AMBER)}'>flagged {m.get('flag_level') or ''}"
                              + (f" (p{m['percentile']:.0f})" if isinstance(m.get("percentile"), (int, float)) else "")
                              + "</font>")
            elif m.get("beneficial"):
                status.append(f"<font color='{_hex(GREEN)}'>beneficial</font>")
            elif m.get("class"):
                status.append(str(m["class"]))
            rows.append([
                Paragraph(_italic(str(m.get("display"))) + (" <b>required</b>" if m.get("required") else ""),
                          st["cell"]),
                Paragraph(_pct(m.get("must_share")), st["cell"]),
                Paragraph(_pct(m.get("can_share")), st["cell"]),
                Paragraph("; ".join(status) or DASH, st["cell"]),
            ])
        if len(members) > 14:
            rows.append([Paragraph(f"<font color='{_hex(INK_FAINT)}'>+{len(members) - 14} more with "
                                   "smaller possible shares</font>", st["cell"]), "", "", ""])
        story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.40, 0.12, 0.12, 0.36)]))


__all__ = [
    "ARM_WORDS",
    "VERDICT_WORDS",
    "composition_block",
    "glance_group_order",
    "grouped_panels",
    "group_for_panel",
    "GLANCE_GROUP_ORDER",
    "PANEL_GROUP",
    "CheckBox",
    "all_functional_rows",
    "context_block",
    "drug_metabolism_block",
    "functional_glance_block",
    "functional_lane_block",
    "longitudinal_block",
    "planner_block",
    "metric_index_block",
    "simulation_block",
    "simulation_detail_block",
]


# --------------------------------------------------------------------------- #
# A09 — composition and community type, inside the gut-community section
# --------------------------------------------------------------------------- #


def _bar(share: float, width_mm: float, colour: colors.Color) -> Table:
    """A proportional bar. Width carries the value; the number is printed too."""
    filled = max(0.4, width_mm * min(1.0, max(0.0, share)))
    bar = Table([[""]], colWidths=[filled * mm], rowHeights=[2.4 * mm])
    bar.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colour),
        ("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
    ]))
    return bar


def _genus_cell(genus: str) -> str:
    """A genus name, or an unnamed genome bin said plainly.

    `GGB45596` is a MetaPhlAn bin with no described species in it. Set in
    italics beside *Bacteroides* it reads as a Latin name a reader might try
    to look up; it has to be marked as what it is.
    """
    if genus.startswith(("GGB", "SGB", "CAG")) or genus == UNRESOLVED:
        label = "organisms not yet named" if genus == UNRESOLVED else genus
        return f"<font color='{_hex(INK_SOFT)}'>{label}</font>" + (
            "<br/><font size='5.8'>an unnamed genome bin</font>" if genus != UNRESOLVED else ""
        )
    return f"<i>{genus}</i>"


def composition_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """The genus composition and, where the frozen model allows, a community type.

    The composition is always available and needs no model. The type is an
    assignment to the nearest medoid of a model frozen on an external cohort,
    with every distance shown - never a cluster computed from whoever else
    happens to be in this run.
    """
    if not view:
        return
    comp = view.get("composition") or {}
    top = comp.get("top_genera") or []
    if not top:
        return

    story.append(CondPageBreak(66 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "<b>What your community is made of</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>genus shares of "
        f"{comp.get('n_taxa', 0)} organisms, as non-overlapping partitions</font>", st["h3"],
    ))
    rows: list[list[Any]] = []
    for entry in top[:10]:
        share = float(entry.get("share") or 0.0)
        rows.append([
            Paragraph(_genus_cell(str(entry.get("genus") or "")), st["cell"]),
            Paragraph(f"<b>{share * 100:.1f}%</b>", st["cell"]),
            _bar(share / max(0.01, float(top[0].get("share") or 1.0)), 62.0, SLATE),
        ])
    other = 1.0 - sum(float(e.get("share") or 0.0) for e in top[:10])
    if other > 0.001:
        rows.append([
            Paragraph(f"<font color='{_hex(INK_SOFT)}'>{comp.get('n_genera', 0) - len(top[:10])} "
                      "other genera</font>", st["cell"]),
            Paragraph(f"<font color='{_hex(INK_SOFT)}'>{other * 100:.1f}%</font>", st["cell"]),
            _bar(other / max(0.01, float(top[0].get("share") or 1.0)), 62.0, RULE),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.26, 0.10, 0.64)], header=False))
    unnamed = sum(
        float(e.get("share") or 0.0) for e in (comp.get("top_genera") or [])
        if str(e.get("genus", "")).startswith(("GGB", "SGB", "CAG"))
    )
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'>Shares of {comp.get('denominator')}. "
        "These are partitions of one total, so they add to 100%: a genus is not counted twice "
        "and nothing is dropped."
        + (f" {unnamed * 100:.0f}% of your community sits in genome bins that have been "
           "assembled and counted but not yet given a name \u2014 a real part of your "
           "community, not a gap in the measurement." if unnamed > 0.01 else "")
        + "</font>", st["small"],
    ))
    _community_type_note(story, st, view)


def _community_type_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    ct = view.get("community_type")
    model = view.get("community_type_model") or {}
    if not ct:
        return
    state = ct.get("state")
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph(
        "<b>Which research community type yours resembles</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>a resemblance, not a diagnosis</font>",
        st["h3"],
    ))
    if state == "assigned":
        taxa = ", ".join(f"<i>{t}</i>" for t in ct.get("distinguishing_taxa") or [])
        lead = (
            f"Your community is closest to the <b>{ct.get('nearest_group')}</b> reference group"
            + (f", the group distinguished in that cohort by {taxa}" if taxa else "")
            + f". It is {ct.get('margin_to_next', 0) * 100:.0f}% closer to that group than to the "
            "next, which is a clear rather than a borderline placement."
        )
    elif state == "mixed":
        lead = f"<b>Your community sits between the reference groups.</b> {ct.get('note')}"
    elif state == "continuous_only":
        lead = f"<b>No named type is reported.</b> {ct.get('note')}"
    else:
        lead = f"<b>No comparison is possible.</b> {ct.get('note')}"
    story.append(Paragraph(f"<font size='7'>{lead}</font>", st["small"]))

    rows: list[list[Any]] = [[Paragraph(h, st["label"])
                              for h in ("REFERENCE GROUP", "DISTANCE", "RESEMBLANCE", "")]]
    distances = ct.get("distances") or []
    sims = {s["group"]: s["similarity"] for s in ct.get("similarities") or []}
    nearest = ct.get("nearest_group")
    for entry in distances:
        group = str(entry["group"])
        similarity = float(sims.get(group, 0.0))
        rows.append([
            Paragraph((f"<b>{group}</b>" if group == nearest else group), st["cell"]),
            Paragraph(f"{entry['root_jsd']:.3f}", st["cell"]),
            Paragraph(f"{similarity * 100:.0f}%", st["cell"]),
            _bar(similarity, 46.0, GREEN if group == nearest else SLATE),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.30, 0.12, 0.13, 0.45)]))

    coverage = ct.get("feature_coverage")
    reference = ct.get("reference_feature_coverage") or {}
    renamed = ct.get("renamed_for_comparison") or []
    notes: list[str] = []
    evaluation = model.get("evaluation") or {}
    if model:
        n_ref = (model.get("cohort") or {}).get("n_reference_samples")
        n_studies = (model.get("cohort") or {}).get("n_studies")
        notes.append(
            f"Groups were fixed once on {n_ref:,} reference samples from {n_studies} studies, "
            f"{model.get('k')} groups chosen by held-out silhouette among those that survived "
            "study-by-study resampling. Your sample is compared with that frozen reference and "
            "never with the other samples analysed alongside it."
            if isinstance(n_ref, int) else ""
        )
    if renamed:
        shown = "; ".join(
            f"<i>{r['as_reported']}</i> as <i>{r['as_in_reference']}</i> ({r['share'] * 100:.0f}%)"
            for r in renamed[:3]
        )
        notes.append(
            "The reference cohort predates several genus reclassifications, so your organisms "
            f"were restated under its names to be comparable: {shown}. Same organisms, older "
            "names."
        )
    if isinstance(coverage, (int, float)) and reference.get("median"):
        median = float(reference["median"])
        notes.append(
            f"The comparison uses the {coverage * 100:.0f}% of your community that the "
            f"reference's genus list names, against {median * 100:.0f}% for a typical reference "
            "sample."
            + (" The gap is newer unnamed genome bins that the older reference catalogue has no "
               "entry for; they are excluded from the comparison rather than forced into a "
               "neighbouring genus." if coverage < median - 0.1 else "")
        )
    if evaluation.get("stable") is False:
        notes.append(f"<b>{evaluation.get('reason')}</b>")
    if notes:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{' '.join(n for n in notes if n)}</font>",
            st["small"],
        ))
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>A community type is a resemblance to a "
        "published composition, not a diagnosis and not a probability of any disease. A "
        "<i>Bacteroides</i>-dominant result does not establish a high-protein or high-fat diet. "
        "Relative composition cannot establish how much microbial mass you carry; that needs an "
        "absolute-load assay.</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# A15 — the complete metric index, in the technical section
# --------------------------------------------------------------------------- #

#: Which section each feature's readings live in, so the index can say where
#: to go rather than only what exists.
FEATURE_HOME: Final[Mapping[str, str]] = {
    "A09": "Your gut community",
    "A10": "Every organism found",
    "A11": "Pathogen screening",
    "A12": "What you can do about it",
    "A13": "Your gut community",
    "A14": "Your information & report context",
    "A16": "What you can do about it",
}

STATE_WORDS: Final[Mapping[str, str]] = {
    "measured": "measured",
    "partial": "partly measured",
    "not_detected_above_assay_threshold": "searched, nothing above the floor",
    "insufficient_coverage": "depth cannot answer",
    "not_assayed": "not searched in this run",
    "unsupported_by_assay": "this assay cannot answer",
    "requires_external_result": "needs a laboratory measurement",
    "not_applicable": "does not apply here",
}


def metric_index_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], extension: Mapping[str, Any] | None
) -> None:
    """Every new measurement, with its ID, value, unit, status and home.

    Spec 0.8.3 §12.3: the complete index belongs in the technical area, not
    in another results section. Per-feature counts separate unique
    measurements from the ones that could not be produced, because a long
    list is not the same as a complete one.
    """
    metrics = list((extension or {}).get("metrics") or [])
    unavailable = list((extension or {}).get("unavailable") or [])
    if not metrics and not unavailable:
        return
    story.append(CondPageBreak(60 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Complete index of the additional measurements", st["h2"]))

    by_feature: dict[str, list[Mapping[str, Any]]] = {}
    for metric in metrics:
        by_feature.setdefault(str(metric.get("feature_id") or "other"), []).append(metric)
    measured = sum(1 for m in metrics if m.get("state") == "measured")
    debt = sum(1 for u in unavailable if u.get("is_engineering_debt"))
    story.append(Paragraph(
        f"<b>{len(metrics)} additional measurement{'s' if len(metrics) != 1 else ''}</b> in this "
        f"report, {measured} of them quantified for your sample, across "
        f"{len(by_feature)} capabilit{'ies' if len(by_feature) != 1 else 'y'}. "
        f"{len(unavailable) - debt} could not be produced from this sample and "
        f"{debt} {'are' if debt != 1 else 'is'} not built yet; the two are listed separately "
        "because only one of them is a limit of the science. Each row below is one unique "
        "measurement \u2014 a reading shown in more than one place is counted once.",
        st["small"],
    ))
    head = ("METRIC ID", "READING", "VALUE", "STATUS", "WHERE IT APPEARS")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    for feature in sorted(by_feature):
        entries = sorted(by_feature[feature], key=lambda m: str(m.get("metric_id")))
        rows.append([
            Paragraph(f"<font size='6.4' color='{_hex(INK_SOFT)}'><b>{feature} \u00b7 "
                      f"{len(entries)} measurement{'s' if len(entries) != 1 else ''}</b></font>",
                      st["cell"]), "", "", "", "",
        ])
        for metric in entries:
            value, unit = metric.get("value"), metric.get("unit") or ""
            shown = DASH if value is None else (
                f"{value:,.4g} {unit}".strip() if isinstance(value, (int, float)) else str(value)
            )
            state = str(metric.get("state") or "")
            rows.append([
                Paragraph(f"<font size='5.8' color='{_hex(INK_SOFT)}'>"
                          f"{metric.get('metric_id')}</font>", st["cell"]),
                Paragraph(f"<font size='6.4'>{metric.get('label')}</font>", st["cell"]),
                Paragraph(f"<font size='6.4'><b>{shown}</b></font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='"
                          f"{_hex(GREEN if state == 'measured' else SLATE)}'>"
                          f"{STATE_WORDS.get(state, state)}</font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                          f"{FEATURE_HOME.get(feature, DASH)}</font>", st["cell"]),
            ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.30, 0.26, 0.16, 0.16, 0.12)]))
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>Every value here also appears in "
        "<b>results.json</b> and <b>all_metrics.tsv</b> beside this report, with its full "
        "precision, its method identifier and the fingerprint of the inputs that produced it. "
        "Rounding on this page is presentation only. An optional laboratory field nobody "
        "supplied is not counted as a reading.</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# A13 — change over time, inside the gut-community section
# --------------------------------------------------------------------------- #


def longitudinal_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """What changed since last time, or what a second sample would add.

    §6.4 is explicit that this section stays visible when recovery cannot
    be estimated. A first sample gets the honest version: what a comparison
    would tell you, and what has to match for it to be a comparison at all.
    """
    if not view:
        return
    if view.get("state") == "single_timepoint":
        _first_sample_block(story, st, view)
        return
    _series_block(story, st, view)


def _first_sample_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    story.append(CondPageBreak(46 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "<b>Change over time</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>your first sample</font>", st["h3"],
    ))
    story.append(Paragraph(f"<font size='7'>{view.get('headline')}</font>", st["small"]))
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in
        ("WHAT A SECOND SAMPLE WOULD ADD", "WHAT MAKES TWO SAMPLES COMPARABLE")
    ]]
    adds = view.get("what_a_second_sample_adds") or []
    needs = view.get("what_makes_them_comparable") or []
    rows.append([
        Paragraph("".join(f"<font size='6.6'>\u00b7 {a}</font><br/>" for a in adds), st["cell"]),
        Paragraph("".join(f"<font size='6.6'>\u00b7 {n}</font><br/>" for n in needs), st["cell"]),
    ])
    story.append(_table(rows, [CONTENT_WIDTH * 0.5, CONTENT_WIDTH * 0.5]))
    _follow_up_note(story, st, view)
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>A single sample describes one day. It "
        "cannot separate what is usual for you from what was true that week, which is the "
        "main thing a second sample settles.</font>", st["small"],
    ))


def _series_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    comparisons = view.get("comparisons") or []
    story.append(CondPageBreak(60 * mm))
    story.append(Spacer(1, 3 * mm))
    n = view.get("n_timepoints", 0)
    story.append(Paragraph(
        "<b>Change over time</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>{n} samples on file</font>", st["h3"],
    ))
    head = ("BETWEEN", "APART", "HOW MUCH THE COMMUNITY MOVED", "COMPARABLE?")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    for row in comparisons:
        bc, jd = row.get("bray_curtis"), row.get("jaccard_distance")
        state = str(row.get("comparability"))
        if bc is None:
            moved = f"<font color='{_hex(INK_FAINT)}'>not compared</font>"
        else:
            # A Jaccard of zero has a plain meaning worth stating: the same
            # organisms in different amounts. Printing "0.000" makes a
            # reader hunt for what the number is telling them.
            if jd is None:
                membership = ""
            elif jd <= 1e-9:
                membership = "; the same organisms throughout, in different amounts"
            else:
                membership = f"; {jd:.3f} by which organisms are present"
            moved = f"<b>{bc:.3f}</b> of the community by abundance{membership}"
        days = row.get("days_apart")
        rows.append([
            Paragraph(f"<font size='6.4'>{row.get('from')} \u2192 {row.get('to')}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{days} days</font>" if days is not None else DASH, st["cell"]),
            Paragraph(f"<font size='6.4'>{moved}</font>", st["cell"]),
            Paragraph(
                f"<font size='6.2' color='"
                f"{_hex(GREEN if state == 'native' else AMBER)}'>{state.replace('_', ' ')}</font>"
                # The reason earns its space only when something differs.
                + (f"<br/><font size='5.8' color='{_hex(INK_SOFT)}'>"
                   f"{row.get('comparability_reason')}</font>" if state != "native" else ""),
                st["cell"],
            ),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.20, 0.10, 0.32, 0.38)]))
    if all(str(r.get("comparability")) == "native" for r in comparisons):
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_SOFT)}'>Every pair is <b>native</b>: the same "
            "lane, the same tool and database versions throughout, so a difference between two "
            "of these samples is a difference in you rather than in how they were measured."
            "</font>", st["small"],
        ))
    _recovery_note(story, st, view)
    _follow_up_note(story, st, view)
    if view.get("repeat_processing_excluded"):
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
            f"{len(view['repeat_processing_excluded'])} result(s) re-analyse a sample already "
            f"on this timeline and are not charted as a change. "
            f"{view.get('repeat_processing_note')}</font>", st["small"],
        ))


def _recovery_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    recovery = view.get("recovery")
    if not recovery:
        return
    state = str(recovery.get("state"))
    if state == "measured":
        last = (recovery.get("trajectory") or [])[-1]
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(
            f"<font size='6.8'><b>Return toward the pre-event community.</b> Your community "
            f"moved {recovery['peak_distance']:.3f} away from its own baseline after the "
            f"recorded event \u2014 well beyond the {recovery['baseline_envelope']:.3f} these "
            f"samples normally differ from each other \u2014 and has since returned "
            f"<b>{last['return_percent']:.0f}%</b> of the way back. This is a retrospective "
            "trajectory: a later sample can revise it, and returning toward the earlier "
            "community does not establish that the earlier community was healthy.</font>",
            st["small"],
        ))
    elif state in {"insufficient_baseline", "no_displacement"}:
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>No return percentage is given.</b> "
            f"{recovery.get('reason')}</font>", st["small"],
        ))


def _follow_up_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    templates = view.get("follow_up") or []
    if not templates:
        return
    story.append(Spacer(1, 1.2 * mm))
    parts = []
    for template in templates:
        parts.append(
            f"<b>{template.get('question')}</b> Repeat in <b>{template.get('interval')}</b>. "
            f"{template.get('rationale')}"
        )
    story.append(Paragraph(
        f"<font size='6.6'><b>When a repeat would tell you something.</b> {' '.join(parts)}"
        "</font>", st["small"],
    ))


# --------------------------------------------------------------------------- #
# A12 — the consolidated plan, at the head of the actions section
# --------------------------------------------------------------------------- #


def planner_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """Start here, then everything else, then the shopping list.

    The per-finding cards below this are unchanged; this is the view across
    them. Every option names the findings it addresses, the setting its
    evidence came from, and the question that would tell you whether it did
    anything.
    """
    if not view:
        return
    start = view.get("start_here") or []
    coverage = view.get("coverage") or {}
    if not start:
        return

    story.append(CondPageBreak(80 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Start here: the options that address the most of your findings",
                           st["h2"]))
    story.append(Paragraph(
        f"Across {coverage.get('n_options', 0)} studied options that between them address "
        f"{coverage.get('n_findings_addressed', 0)} of your readings, these are the ones that "
        "reach the most of them on the best-matched evidence. This is a <b>prioritisation with "
        "its reasons printed</b>, not a prediction that any of them will work for you, and the "
        "full list follows.", st["body"],
    ))
    for i, option in enumerate(start, start=1):
        _option_card(story, st, option, position=i)

    _ranking_policy_note(story, st, view)
    _catalogue_summary(story, st, view)
    _checklist_block(story, st, view)


def _ranking_policy_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """In what order, and why, so the ordering can be argued with.

    A ranked list whose rule is hidden asks to be taken on trust. The last
    tie-break is the option's own identifier, which is what stops two
    equally-placed options from swapping between runs and looking like a
    change in the evidence.
    """
    policy = (view or {}).get("ranking_policy") or {}
    order = [str(x) for x in (policy.get("order") or [])]
    note = (view or {}).get("start_here_note")
    if not order and not note:
        return
    story.append(Spacer(1, 1.5 * mm))
    if note:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{note}</font>", st["small"],
        ))
    if order:
        steps = "; then ".join(order)
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_FAINT)}'><b>Ordered by</b> {steps}."
            "</font>", st["small"],
        ))


def _option_card(
    story: list[Any], st: Mapping[str, ParagraphStyle], option: Mapping[str, Any], *, position: int
) -> None:
    lane = str(option.get("best_evidence_lane"))
    tone = GREEN if lane in {"A", "B"} else SLATE
    flagged = option.get("addresses_flagged") or []
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        f"<b>{position}. {option.get('label')}</b> &nbsp;"
        f"<font size='6.6' color='{_hex(INK_SOFT)}'>{option.get('category_label')}</font>",
        st["h3"],
    ))
    rows: list[list[Any]] = [
        [Paragraph("<font size='6.4'><b>What it addresses</b></font>", st["cell"]),
         Paragraph(
             f"<font size='6.4'>{'; '.join(option.get('addresses') or [])}</font>"
             + (f"<br/><font size='6.2' color='{_hex(AMBER)}'>including "
                f"{len(flagged)} reading(s) this report flagged: {'; '.join(flagged)}</font>"
                if flagged else ""),
             st["cell"]),
         ],
        [Paragraph("<font size='6.4'><b>Evidence</b></font>", st["cell"]),
         Paragraph(
             f"<font size='6.4' color='{_hex(tone)}'>{option.get('best_evidence_words')}</font>"
             f"<font size='6.4'> \u2014 {option.get('n_supporting', 0)} supporting and "
             f"{option.get('n_against', 0)} contrary finding(s) on record. "
             f"{option.get('counts_note')}</font>", st["cell"]),
         ],
    ]
    if option.get("population"):
        rows.append([
            Paragraph("<font size='6.4'><b>Studied in</b></font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{option['population']}</font>", st["cell"]),
        ])
    if option.get("studied_exposure"):
        rows.append([
            Paragraph("<font size='6.4'><b>What the study did</b></font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{option['studied_exposure']}</font>"
                      f"<br/><font size='6' color='{_hex(INK_SOFT)}'>"
                      f"{option.get('exposure_note')}</font>", st["cell"]),
        ])
    if option.get("monitoring_question"):
        rows.append([
            Paragraph("<font size='6.4'><b>How you would know</b></font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{option['monitoring_question']}</font>", st["cell"]),
        ])
    reasons = option.get("rank_reasons") or []
    if reasons:
        rows.append([
            Paragraph("<font size='6.4'><b>Why it is here</b></font>", st["cell"]),
            Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                      + ", ".join(r.replace("_", " ") for r in reasons) + "</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * 0.20, CONTENT_WIDTH * 0.80], header=False))


def _clip(text: str, limit: int) -> str:
    """Shorten on a word boundary. A name cut mid-word reads as a typo."""
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0].rstrip(" ,;(")
    return f"{cut}\u2026"


def _catalogue_summary(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    catalogue = view.get("catalogue") or []
    coverage = view.get("coverage") or {}
    if not catalogue:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>Everything else that has been studied</b> &nbsp;"
                           f"<font size='7' color='{_hex(INK_FAINT)}'>each one's card follows "
                           "under the reading it belongs to</font>", st["h3"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"])
                              for h in ("CATEGORY", "OPTIONS", "BEST-MATCHED EVIDENCE")]]
    for group in catalogue:
        options = group.get("options") or []
        names = "; ".join(_clip(str(o.get("label")), 56) for o in options[:4])
        if len(options) > 4:
            names += f"; +{len(options) - 4} more"
        best = min(
            (str(o.get("best_evidence_lane", "X")) for o in options),
            key=lambda lane: "ABCDEX".find(lane),
            default="X",
        )
        rows.append([
            Paragraph(f"<b>{group.get('category_label')}</b><br/>"
                      f"<font size='6.2' color='{_hex(INK_SOFT)}'>{group.get('n_options')}"
                      "</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{names}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{VERDICT_LANE_WORDS.get(best, best)}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.18, 0.58, 0.24)]))
    without = coverage.get("n_without_human_evidence") or 0
    if without:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{without} of these have been studied "
            "only in animals, cells or by mechanism. They are listed here with their actual "
            "setting named rather than hidden: nobody having run the trial is not the same as "
            "nothing being worth trying, and it is not the same as evidence that it works."
            "</font>", st["small"],
        ))
    excluded = view.get("excluded_for_this_person") or []
    if excluded:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(AMBER)}'>{len(excluded)} option(s) are kept out of "
            "the start-here list for you specifically: "
            + "; ".join(
                f"<b>{o.get('label')}</b> ({(o.get('exclusions') or [{}])[0].get('detail')})"
                for o in excluded[:3]
            )
            + ". Their evidence is unchanged and their cards still follow.</font>", st["small"],
        ))


VERDICT_LANE_WORDS: Final[Mapping[str, str]] = {
    "A": "human trial in this condition",
    "B": "human trial in a related group",
    "C": "human study, other design",
    "D": "animal or ex-vivo study",
    "E": "laboratory study",
    "X": "mechanistic or indirect",
}


class CheckBox(Flowable):
    """An empty box to tick, drawn rather than typed.

    The ballot-box character is not in the report's font, and the
    substituted glyph is a *filled* square - which reads as already done
    and is worse than no box at all. Two millimetres of stroked rectangle
    has no such opinion.
    """

    def __init__(self, size: float = 2.4 * mm) -> None:
        super().__init__()
        self.width = self.height = size

    def wrap(self, *_: float) -> tuple[float, float]:
        return self.width, self.height

    def draw(self) -> None:
        self.canv.setStrokeColor(SLATE)
        self.canv.setLineWidth(0.5)
        self.canv.rect(0, 0, self.width, self.height, stroke=1, fill=0)


def _checklist_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    checklist = view.get("food_checklist") or []
    patterns = view.get("dietary_patterns") or []
    if not checklist and not patterns:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>What to put on the list</b> &nbsp;"
                           f"<font size='7' color='{_hex(INK_FAINT)}'>deduplicated across every "
                           "option above</font>", st["h3"]))
    if checklist:
        # One line per item rather than per group, each with a box to tick:
        # §12.3 asks for a printable checklist, and a comma-separated list
        # of items on one row is a table, not something anybody can take to
        # a shop and mark off.
        rows: list[list[Any]] = [[Paragraph(h, st["label"])
                                  for h in ("", "ITEM", "GROUP", "WHAT IT IS FOR")]]
        for group in checklist:
            for item in group.get("items") or []:
                targets = sorted(set(item.get("addresses") or []))
                rows.append([
                    CheckBox(),
                    Paragraph(f"<b>{str(item.get('item')).capitalize()}</b>", st["cell"]),
                    Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                              f"{group.get('group_label')}</font>", st["cell"]),
                    Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                              + "; ".join(targets[:3])
                              + (f"; +{len(targets) - 3} more" if len(targets) > 3 else "")
                              + "</font>", st["cell"]),
                ])
        story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.05, 0.25, 0.28, 0.42)]))
    if patterns:
        story.append(Paragraph(
            "<font size='6.6'><b>Ways of eating, rather than items to buy:</b> "
            + "; ".join(str(p.get("pattern")) for p in patterns[:6])
            + (f"; +{len(patterns) - 6} more" if len(patterns) > 6 else "")
            + ". Each is a whole pattern that was studied as a pattern; the card for each one "
            "gives the protocol that was actually followed.</font>", st["small"],
        ))
    story.append(Paragraph(
        f"<font size='6.2' color='{_hex(INK_FAINT)}'>{view.get('checklist_note', '')}</font>",
        st["small"],
    ))


# --------------------------------------------------------------------------- #
# A01-A08 — the functional lane, in the metabolism detail section
# --------------------------------------------------------------------------- #

#: The reader-facing groupings of §12.1, in order, and which views feed each.
FUNCTION_GROUPS: Final[tuple[tuple[str, str, tuple[str, ...]], ...]] = (
    ("Fibre & Dietary Substrates", "substrates", ("substrates",)),
    ("Fermentation & Cross-feeding", "fermentation", ("fermentation",)),
    ("Vitamins & Nutrients", "vitamins", ("vitamins",)),
    ("Protein & Nitrogen Metabolism", "nitrogen", ("nitrogen",)),
    ("Gut-Brain & Metabolic Signalling", "biotransformation", ("A05", "A07")),
    ("Plant-Compound Conversion", "biotransformation", ("A06",)),
    ("Mucus, Bile & Other Transformations", "biotransformation", ("A08",)),
    # Spec 0.8.3 §12.1 asks for the urate, polyamine, glutathione, ethanol and
    # drug-metabolism detail in this part of the report, and asks for smaller
    # subsections rather than one overfilled card. These are those.
    ("Urate & Purines", "additional_capacities", ("capacity.urate_anaerobic",)),
    ("Polyamines", "additional_capacities",
     ("capacity.putrescine", "capacity.agmatine", "capacity.spermidine")),
    ("Glutathione & Ethanol", "additional_capacities",
     ("capacity.glutathione", "capacity.ethanol")),
)

#: Spec 0.8.3 §12.1: the reader-facing groups, in order, for both the
#: functions overview and the functions detail. The first seven are the
#: specification's named groupings; the three after them are §5.8's
#: subsections, kept separate because §5.8 asks for smaller subsections
#: rather than one overfilled card; and virulence is last because a toxin is
#: not one of the seven chemistries and has no sensible home among them.
GLANCE_GROUP_ORDER: Final[tuple[str, ...]] = (
    "Fibre & Dietary Substrates",
    "Fermentation & Cross-feeding",
    "Vitamins & Nutrients",
    "Protein & Nitrogen Metabolism",
    "Gut-Brain & Metabolic Signalling",
    "Plant-Compound Conversion",
    "Mucus, Bile & Other Transformations",
    "Urate & Purines",
    "Polyamines",
    "Glutathione & Ethanol",
    "Toxins & Virulence",
)

#: Which of those groups each existing panel belongs to. Assigned by the
#: chemistry, not by the panel's own category field: the categories are an
#: older display vocabulary that splits vitamins from nutrients and lumps
#: bile acids in with oxalate, and §12.1 asks for existing and new readings
#: to sit together under one set of headings.
PANEL_GROUP: Final[Mapping[str, str]] = {
    "carbohydrates": "Fibre & Dietary Substrates",
    # §12.1: existing SCFAs plus acetate, lactate, succinate and hydrogen
    # production and consumption.
    "butyrate": "Fermentation & Cross-feeding",
    "propionate": "Fermentation & Cross-feeding",
    "fermentation": "Fermentation & Cross-feeding",
    "h2s": "Fermentation & Cross-feeding",
    "methane": "Fermentation & Cross-feeding",
    "b12": "Vitamins & Nutrients",
    "b6": "Vitamins & Nutrients",
    "biotin": "Vitamins & Nutrients",
    "folate": "Vitamins & Nutrients",
    "k2": "Vitamins & Nutrients",
    "niacin": "Vitamins & Nutrients",
    "pantothenate": "Vitamins & Nutrients",
    "riboflavin": "Vitamins & Nutrients",
    "thiamine": "Vitamins & Nutrients",
    # §12.1 puts the aromatic metabolites here with the amino-acid and
    # ammonia functions, which is also where the A04 module links them.
    "bcaa": "Protein & Nitrogen Metabolism",
    "nitrogen": "Protein & Nitrogen Metabolism",
    "urease": "Protein & Nitrogen Metabolism",
    "peptides": "Protein & Nitrogen Metabolism",
    "indole": "Protein & Nitrogen Metabolism",
    "ipa": "Protein & Nitrogen Metabolism",
    "pcresol": "Protein & Nitrogen Metabolism",
    "cutc": "Protein & Nitrogen Metabolism",
    # Imidazole propionate is made from histidine, so it belongs with the
    # amino-acid chemistry and not with the plant compounds - urdA is not a
    # urolithin marker, which §5.6 is emphatic about.
    "urda": "Protein & Nitrogen Metabolism",
    "gaba": "Gut-Brain & Metabolic Signalling",
    "gababreak": "Gut-Brain & Metabolic Signalling",
    "dopamine": "Gut-Brain & Metabolic Signalling",
    "histamine": "Gut-Brain & Metabolic Signalling",
    "tryptamine": "Gut-Brain & Metabolic Signalling",
    # A07's experimental protein sits with the GLP-1 mechanisms it belongs to.
    "p9": "Gut-Brain & Metabolic Signalling",
    "glucosinolate": "Plant-Compound Conversion",
    "urolithin": "Plant-Compound Conversion",
    "equol": "Plant-Compound Conversion",
    "bsh": "Mucus, Bile & Other Transformations",
    "bai": "Mucus, Bile & Other Transformations",
    "mucinlps": "Mucus, Bile & Other Transformations",
    "bglucuronidase": "Mucus, Bile & Other Transformations",
    "oxalate": "Mucus, Bile & Other Transformations",
    # Each of these sits with the capacity computed over it, so a panel row
    # and its route reading are never in different parts of the section.
    "urate": "Urate & Purines",
    "putrescine": "Polyamines",
    "agmatine": "Polyamines",
    "spermidine": "Polyamines",
    "glutathione": "Glutathione & Ethanol",
    "ethanol": "Glutathione & Ethanol",
    "crc_virulence": "Toxins & Virulence",
}


def group_for_panel(name: str) -> str:
    """The reader-facing group a panel is shown under.

    An unmapped panel goes to the transformations group rather than
    disappearing, and the test over this map is what stops a new panel
    quietly landing there.
    """
    return PANEL_GROUP.get(str(name), "Mucus, Bile & Other Transformations")


def grouped_panels(rows: Sequence[Any]) -> list[tuple[str, list[Any]]]:
    """Panel rows in the §12.1 group order, empty groups omitted."""
    buckets: dict[str, list[Any]] = {}
    for row in rows:
        buckets.setdefault(group_for_panel(getattr(row, "panel", "")), []).append(row)
    out = [(label, buckets.pop(label)) for label in GLANCE_GROUP_ORDER if label in buckets]
    # Anything left is a group name nobody declared; shown rather than lost.
    out.extend(sorted(buckets.items()))
    return out


#: Which subsection each additional capacity is reported under.
_CAPACITY_GROUP: Final[Mapping[str, str]] = {
    capacity_id: label
    for label, view, keys in FUNCTION_GROUPS if view == "additional_capacities"
    for capacity_id in keys
}

STATE_TONE: Final[Mapping[str, Any]] = {
    "substrate_specific": GREEN, "supported": GREEN, "complete": GREEN, "present": GREEN,
    "class_level_only": AMBER, "partial": AMBER, "negative_control_matched": AMBER,
    "insufficient": SLATE, "absent": SLATE, "not_assayed": INK_FAINT,
    "not_applicable": INK_FAINT,
}

STATE_WORD: Final[Mapping[str, str]] = {
    "substrate_specific": "specific evidence",
    "class_level_only": "class-level only",
    "insufficient": "not enough to say",
    "supported": "supported",
    "partial": "part of it",
    "absent": "searched, not found",
    "not_assayed": "not searched for",
    "not_applicable": "no such route",
    "complete": "complete route",
    "present": "present",
    "negative_control_matched": "withheld",
}


def _state_cell(state: str, st: Mapping[str, ParagraphStyle]) -> Paragraph:
    tone = STATE_TONE.get(state, SLATE)
    return Paragraph(
        f"<font size='6.4' color='{_hex(tone)}'>{STATE_WORD.get(state, state)}</font>",
        st["cell"],
    )


def functional_lane_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None,
    aggregates: Mapping[str, Sequence[str]] | None = None,
) -> None:
    """The A01-A08 readings, grouped as §12.1 asks and stating their limits.

    Each row carries what it found, how specific that is, and the one
    conclusion it must not support. A row that was never searched for says
    so rather than reading as a negative.
    """
    if not views:
        return
    # The same exclusion the overview makes. A reading that restates its
    # panel has no row up there to link back to, and a link to a row that
    # was removed is a destination nothing defines.
    rows = [r for r in all_functional_rows(views, None, aggregates)
            if not r.get("duplicates_panel")]
    if not rows:
        return
    story.append(CondPageBreak(70 * mm))
    story.append(Spacer(1, 3 * mm))
    # §12.3: a detail card carries Back to overview and Contents.
    story.append(back_link_line(section_dest("functions"), "the functions overview"))
    story.append(Paragraph(
        "Substrates, routes and transformations, one at a time", st["h2"]))
    story.append(Paragraph(
        "What your microbes carry the genes to do, grouped by the kind of chemistry. "
        "Each row says how specific the evidence is, because a gene family that acts on "
        "a whole class of molecules is not evidence about one member of it \u2014 and each "
        "says the one conclusion it does not support.", st["body"],
    ))
    # This used to be a flat table of every reading with columns for the
    # genes found and an evidence grade. Each reading has a full card below
    # now, in the shape every other metric uses, so the table said the same
    # things twice - once in a layout of its own, with gene symbols where
    # the rest of the report shows a verdict.


    _evidence_level_note(story, st, views)
    _reading_detail_cards(story, st, views, rows)
    _network_block(story, st, views)
    _shared_genes_block(story, st, views)


def _reading_detail_cards(
    story: list[Any], st: Mapping[str, ParagraphStyle],
    views: Mapping[str, Any] | None, rows: Sequence[Mapping[str, Any]],
) -> None:
    """One card per reading: what it is made of, and what it rests on.

    The overview line says where a reading sits. It does not say which
    genes produced that, how each of them sits on its own, or which scored
    panel carries the headline for the same chemistry - and without those
    a reader cannot check the number or act on it.

    The components are the point. A reading is a claim about a set of
    genes, and showing the set with each gene's own position is what turns
    a bar into something arguable.
    """
    scores = ((views or {}).get("reading_scores") or {})
    if not rows:
        return
    story.append(CondPageBreak(60 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Each reading in full", st["h2"]))
    story.append(Paragraph(
        "One card per reading, all in the same four parts: what is driving it, what it "
        "means, how it might be moved, and the research behind it. Each card names the "
        "genes the reading is measured from and where each of those sits in the reference "
        "cohort, so the number above can be taken apart rather than taken on trust.",
        st["body"],
    ))
    for row in rows:
        _reading_card(story, st, row, scores.get(str(row.get("reading_id"))) or {})


def _reading_card(
    story: list[Any], st: Mapping[str, ParagraphStyle],
    row: Mapping[str, Any], score: Mapping[str, Any],
) -> None:
    from openbiota.pdfreport import PercentileBar, _ordinal  # noqa: PLC0415

    story.append(CondPageBreak(30 * mm))
    story.append(Spacer(1, 2.5 * mm))
    # §12.3: the card is what the overview row links to, and links back.
    dest = reading_dest(str(row["group"]), str(row["label"]), detail=True)
    story.append(Anchor(dest))
    # §12.3 and the shared card: every detail page offers the way back to
    # the line the reader came from.
    story.append(back_link_line(
        reading_dest(str(row["group"]), str(row["label"])), "at a glance"))
    owner = owner_of(str(row.get("reading_id") or ""))
    story.append(Paragraph(
        f"<b>{row['label']}</b>"
        + (f"<font size='6.6' color='{_hex(INK_SOFT)}'> \u00b7 headline for this chemistry is "
           f"the <b>{owner}</b> panel</font>" if owner else ""),
        st["h3"],
    ))

    pct = row.get("percentile")
    if pct is not None:
        story.append(_table([[
            PercentileBar(width=CONTENT_WIDTH * 0.30, percentile=float(pct),
                          higher_means=str(row.get("higher_means") or "unclear"),
                          marker_colour=_band_colour(float(pct),
                                                     str(row.get("higher_means") or "unclear"))),
            Paragraph(f"<b>{_ordinal(pct)}</b>"
                      f"<font size='6.2' color='{_hex(INK_SOFT)}'> of the reference cohort"
                      "</font>", st["cell"]),
        ]], [CONTENT_WIDTH * 0.34, CONTENT_WIDTH * 0.66], header=False))
    elif score.get("unscored_reason"):
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{score['unscored_reason']}</font>",
            st["small"]))

    # --- the shared metric card, in the order `openbiota.metriccard` fixes.
    # Do not reorder these or add a fifth heading to this card alone.
    from openbiota import metriccard as MC  # noqa: PLC0415

    # 1 · What is driving this. The organisms come first: a reading is
    # caused by organisms, and a card that lists only genes cannot tell
    # anybody what is behind their number or what to act on.
    story.append(Paragraph(f"<b>{MC.CARD_SECTIONS[0]}</b>", st["h3"]))
    contributors = [c for c in (score.get("contributors") or []) if isinstance(c, dict)]
    _reading_organisms(story, st, contributors)
    if contributors:
        cells: list[list[Any]] = [[
            Paragraph(h, st["label"])
            for h in ("THE GENE BEHIND IT", "HOW MUCH IS IN YOU", "WHERE THAT SITS")
        ]]
        for c in contributors:
            value = c.get("value")
            position = c.get("percentile")
            cells.append([
                Paragraph(_gene_words(c), st["cell"]),
                Paragraph(
                    f"{value:.1f}<font size='5.8' color='{_hex(INK_SOFT)}'> copies per 100 "
                    "genomes</font>" if isinstance(value, int | float) else
                    f"<font color='{_hex(INK_FAINT)}'>{DASH}</font>", st["cell"]),
                Paragraph(
                    f"{_ordinal(position)}"
                    + (f"<font size='5.8' color='{_hex(INK_SOFT)}'> of "
                       f"{c['cohort_n']} people</font>" if c.get("cohort_n") else "")
                    if position is not None else
                    f"<font size='6' color='{_hex(INK_FAINT)}'>no reference range</font>",
                    st["cell"]),
            ])
        story.append(_table(cells, [CONTENT_WIDTH * w for w in (0.34, 0.30, 0.36)]))
    else:
        story.append(Paragraph(
            "No gene behind this reading was measured in this sample, so nothing here "
            "identifies a source for it.", st["body_ink"]))

    # 2 · What your reading means
    story.append(Paragraph(f"<b>{MC.CARD_SECTIONS[1]}</b>", st["h3"]))
    story.append(Paragraph(_reading_meaning(row, score), st["body_ink"]))

    # 3 · How can I improve this
    MC.improvement_block(story, st, higher_means=str(row.get("higher_means") or "unclear"))

    # 4 · The research behind it
    story.append(Paragraph(f"<b>{MC.CARD_SECTIONS[3]}</b>", st["h3"]))
    story.append(Paragraph(
        f"<font size='6.8'>{row.get('note') or row.get('found') or ''}</font>", st["small"]))
    if row.get("caveat"):
        story.append(Paragraph(
            f"<font size='6.6' color='{_hex(INK_SOFT)}'><b>What it does not mean.</b> "
            f"{row['caveat']}</font>", st["small"]))


def _reading_organisms(
    story: list[Any], st: Mapping[str, ParagraphStyle],
    contributors: Sequence[Mapping[str, Any]],
) -> None:
    """The organisms carrying this reading, in this sample, largest first.

    Summed across the reading's genes, because a route carried by one
    organism and a route spread over six are different findings and the
    gene counts alone do not distinguish them. This is the part of the
    card a reader can act on.
    """
    totals: dict[str, int] = {}
    for contributor in contributors:
        for organism in contributor.get("organisms") or ():
            name = str(organism.get("organism") or "").strip()
            if name:
                totals[name] = totals.get(name, 0) + int(organism.get("fragments") or 0)
    if not totals:
        story.append(Paragraph(
            "<font size='6.6'>No organism in this sample could be matched to the genes "
            "behind this reading, so nothing here identifies a source for it.</font>",
            st["small"]))
        return
    grand = sum(totals.values()) or 1
    ranked = sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))[:6]
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"])
        for h in ("ORGANISM CARRYING IT", "SHARE OF THIS READING", "FRAGMENTS")
    ]]
    for name, fragments in ranked:
        rows.append([
            Paragraph(f"<i>{species_label(name)}</i>", st["cell"]),
            Paragraph(f"<b>{fragments / grand:.0%}</b>", st["cell"]),
            Paragraph(f"<font size='6.4' color='{_hex(INK_SOFT)}'>{fragments:,}</font>",
                      st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.52, 0.26, 0.22)]))
    if len(totals) > len(ranked):
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_SOFT)}'>and {len(totals) - len(ranked)} "
            "more organisms with smaller shares.</font>", st["small"]))


def _gene_words(contributor: Mapping[str, Any]) -> str:
    """A gene a reader can do something with.

    "pta, acsB, cooS" tells somebody nothing about whether their acetate
    production is high. The panel entry already carries a description -
    "pta (phosphate acetyltransferase)" - so the row leads with what the
    enzyme does and keeps the symbol underneath, small, because that is
    what makes the reading checkable against a database.
    """
    gene = str(contributor.get("gene") or "")
    label = str(contributor.get("label") or gene)
    # Panel labels read "pta (phosphate acetyltransferase)"; the part in
    # brackets is the English and the part before it is the symbol again.
    words = label.split("(", 1)[1].rstrip(")").strip() if "(" in label else label
    words = words or gene
    return (
        f"<b>{words[:1].upper()}{words[1:]}</b>"
        + (f"<br/><font size='5.6' color='{_hex(INK_FAINT)}'>{gene}</font>" if gene else "")
    )


def _reading_meaning(row: Mapping[str, Any], score: Mapping[str, Any]) -> str:
    """The reading in a sentence, whether or not it has a position."""
    from openbiota.pdfreport import _ordinal  # noqa: PLC0415

    label = str(row.get("label") or "This reading")
    pct = row.get("percentile")
    if pct is None:
        return (
            f"{label} has no position in the reference cohort. "
            + str(score.get("unscored_reason") or
                  "No distribution has been built for it yet.")
        )
    where = ("higher than most people's" if pct >= 75 else
             "lower than most people's" if pct <= 25 else
             "in the middle of the reference range")
    return (
        f"Your capacity for {label.lower()} sits {where}, at the {_ordinal(pct)} of "
        f"{score.get('cohort_n') or 99} reference samples. This counts the genes for the "
        "step, which is a capacity and not a measurement of how much is made."
    )


#: What an evidence level means for a reader deciding how much weight to put
#: on a reading. Spelled out because the bare word is a category name.
_EVIDENCE_LEVEL_WORDS: Final[Mapping[str, str]] = {
    "research": (
        "These are research-grade readings. The biology behind them is published and the "
        "genes are real, but none of them has the weight of a clinical measurement, and "
        "none should be acted on as though it had."
    ),
    "established": (
        "The biology behind these readings is well established, though a gene count is "
        "still not a measurement of the molecule."
    ),
    "exploratory": (
        "These readings are exploratory. They are shown so the evidence is visible, not "
        "because it is settled."
    ),
}


def _evidence_level_note(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """How much weight the additional capacities can carry."""
    level = ((views or {}).get("additional_capacities") or {}).get("evidence_level")
    if not level:
        return
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'><b>Evidence level: {level}.</b> "
        + _EVIDENCE_LEVEL_WORDS.get(str(level), "")
        + "</font>", st["small"],
    ))
    _paired_readings_block(story, st, views)
    _three_concepts_block(story, st, views)
    _host_steps_block(story, st, views)
    _decoys_block(story, st, views)


def _network_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Which fermentation steps this community has the genes for.

    Drawn as a list of steps rather than a graph, because the honest thing
    to show is which of them are supported and which are only compatible
    with the literature \u2014 a picture makes every arrow look equally real.
    None of them is a flow: an edge is a capacity, and carries no rate.
    """
    network = ((views or {}).get("fermentation") or {}).get("network") or {}
    edges = [e for e in (network.get("edges") or []) if isinstance(e, dict)]
    unsupported = [e for e in (network.get("unsupported_edges") or []) if isinstance(e, dict)]
    if not edges and not unsupported:
        return
    story.append(CondPageBreak(44 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("The fermentation steps your community can run", st["h3"]))
    story.append(Paragraph(
        "Each line is one step, from what goes in to what comes out. A <b>supported</b> step "
        "is one this community carries the genes for. A step marked <i>literature only</i> is "
        "consistent with what is known but is not established by these genes, and is shown so "
        "the gap is visible rather than filled in. None of these is a flow: they are things "
        "the community could do, with no rate attached, and a step that uses a compound is not "
        "evidence the compound is there.", st["body"],
    ))
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in ("FROM", "TO", "ROUTE", "EVIDENCE")
    ]]
    for edge in [*edges, *unsupported]:
        supported = str(edge.get("support")) == "supported"
        rows.append([
            Paragraph(f"<font size='6.6'>{edge.get('from_label', '')}</font>", st["cell"]),
            Paragraph(f"<b><font size='6.6'>{edge.get('to_label', '')}</font></b>", st["cell"]),
            Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                      f"{edge.get('route_id', '')}</font>", st["cell"]),
            Paragraph(
                f"<font size='6.4' color='{_hex(GREEN if supported else INK_FAINT)}'>"
                + ("genes present" if supported else "literature only")
                + "</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.24, 0.24, 0.30, 0.22)]))
    limits = [str(x) for x in (network.get("limitations") or [])]
    if limits:
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_SOFT)}'>{' '.join(limits)}</font>", st["small"],
        ))


def _shared_genes_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Why the substrate readings must not be added up.

    One gene can serve several substrates \u2014 amyA counts toward starch and
    resistant starch and IMO alike \u2014 so a reader who totals the cards
    counts it once per card. Saying which genes are shared is what makes
    that visible instead of leaving a plausible sum to be made.
    """
    shared_view = ((views or {}).get("substrates") or {}).get("shared_genes") or {}
    shared = shared_view.get("shared") or {}
    if not shared:
        return
    story.append(CondPageBreak(34 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("Why these fibre readings cannot be added together", st["h3"]))
    story.append(Paragraph(
        f"<b>{shared_view.get('n_shared_genes', 0)}</b> of the "
        f"<b>{shared_view.get('n_genes', 0)}</b> genes behind these fibre readings serve more "
        "than one fibre. The same stretch of DNA that helps digest starch also helps digest "
        "resistant starch, so each reading counts it and a total over the readings counts it "
        "several times.", st["body"],
    ))
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in ("GENE", "COUNTS TOWARD")
    ]]
    for gene, substrates in sorted(shared.items()):
        names = ", ".join(str(s).split(".")[-1].replace("_", " ") for s in substrates)
        rows.append([
            Paragraph(f"<b><i>{gene}</i></b>", st["cell"]),
            Paragraph(f"<font size='6.4'>{names}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.78)]))
    if shared_view.get("note"):
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{shared_view['note']}</font>",
            st["small"],
        ))


def _paired_readings_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Two gene counts that are not a flux, shown side by side and not divided.

    Synthesis and degradation genes for GABA are separate measurements. A
    ratio of them would look like a net balance and is nothing of the kind:
    both count genes, neither counts molecules.
    """
    balance = ((views or {}).get("biotransformation") or {}).get("gaba_balance") or {}
    if not balance or balance.get("synthesis") is None and balance.get("degradation") is None:
        return
    story.append(CondPageBreak(30 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("GABA: two counts, deliberately not combined", st["h3"]))

    def cell(value: Any) -> Paragraph:
        if value is None:
            return Paragraph(
                f"<font color='{_hex(INK_FAINT)}'>{DASH}</font>"
                f"<font size='6.2' color='{_hex(INK_SOFT)}'> no supported reading</font>",
                st["cell"],
            )
        return Paragraph(f"<b>{float(value):.1f}</b>", st["cell"])

    rows = [
        [Paragraph(h, st["label"]) for h in ("READING", "VALUE", "UNIT")],
        [Paragraph("Genes for making GABA", st["cell"]), cell(balance.get("synthesis")),
         Paragraph(f"<font size='6.4'>{balance.get('unit', '')}</font>", st["cell"])],
        [Paragraph("Genes for breaking GABA down", st["cell"]), cell(balance.get("degradation")),
         Paragraph(f"<font size='6.4'>{balance.get('unit', '')}</font>", st["cell"])],
    ]
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.44, 0.22, 0.34)]))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'>{balance.get('note', '')}</font>",
        st["small"],
    ))


#: Acronyms that a sentence-casing rule would spoil.
_CONCEPT_ACRONYMS: Final[Mapping[str, str]] = {
    "bcaa": "BCAA", "scfa": "SCFA", "gaba": "GABA", "tma": "TMA", "hmo": "HMO",
}


def _concept_label(key: str) -> str:
    words = key.replace("_", " ").split()
    return " ".join(
        _CONCEPT_ACRONYMS.get(w.casefold(), w if i else w.capitalize())
        for i, w in enumerate(words)
    )


def _three_concepts_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Three nitrogen readings that get confused for one another."""
    concepts = ((views or {}).get("nitrogen") or {}).get("three_different_concepts") or {}
    named = [(k, v) for k, v in concepts.items() if k != "note" and isinstance(v, str)]
    if not named:
        return
    story.append(CondPageBreak(30 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("Three protein readings that are not the same reading", st["h3"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in ("READING", "WHAT IT MEANS")]]
    for key, meaning in named:
        rows.append([
            Paragraph(f"<b>{_concept_label(key)}</b>", st["cell"]),
            Paragraph(f"<font size='6.4'>{meaning}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.30, 0.70)]))
    if concepts.get("note"):
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{concepts['note']}</font>", st["small"],
        ))


def _host_steps_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Where your own liver finishes what the microbes started.

    Worth stating plainly: these end products are not microbial readings,
    and a stool gene count cannot predict how much of one you make.
    """
    steps = [s for s in (((views or {}).get("nitrogen") or {}).get("host_steps") or [])
             if isinstance(s, dict)]
    if not steps:
        return
    story.append(CondPageBreak(30 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("Where your own body finishes the job", st["h3"]))
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"])
        for h in ("MICROBES MAKE", "YOUR ENZYME", "WHICH BECOMES")
    ]]
    for step in steps:
        rows.append([
            Paragraph(f"<b>{step.get('microbial_product', '')}</b>", st["cell"]),
            Paragraph(f"<font size='6.4'>{step.get('host_enzyme', '')}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{step.get('host_product', '')}</font>", st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.28, 0.38, 0.34)]))
    notes = [str(s.get("note")) for s in steps if s.get("note")]
    if notes:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{notes[0]}</font>", st["small"],
        ))


def _decoys_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """What each hard panel is built to tell itself apart from.

    A panel for a rare reaction is only as good as the near-miss it rejects.
    Naming the decoy is what lets a reader see that \u201cnot found\u201d here
    means the look-alike was excluded, not that nothing was looked at.
    """
    controls = ((views or {}).get("biotransformation") or {}).get("negative_controls") or {}
    named = [(k, v) for k, v in controls.items() if isinstance(v, dict) and v]
    if not named:
        return
    story.append(CondPageBreak(34 * mm))
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("What the difficult panels are told apart from", st["h3"]))
    story.append(Paragraph(
        "Some reactions are carried by enzymes that look very like commoner ones. For those, "
        "the panel is built around the look-alike it has to reject: if a read matches the "
        "decoy, the panel has failed rather than found something. These are the decoys.",
        st["body"],
    ))
    rows: list[list[Any]] = [[
        Paragraph(h, st["label"]) for h in ("READING", "THE LOOK-ALIKE", "WHY IT IS HARD")
    ]]
    for reading_id, decoys in named:
        label = reading_id.split(".")[-1].replace("_", " ")
        for decoy, why in decoys.items():
            rows.append([
                Paragraph(f"<b>{label}</b>", st["cell"]),
                Paragraph(f"<font size='6.4'><i>{decoy}</i></font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{why}</font>", st["cell"]),
            ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.18, 0.60)]))


#: Which scored panel each functional reading is a component of.
#:
#: None of these readings is independent. Every §5.8 capacity is computed
#: over a panel that already has a percentile and a slider; the five reused
#: vitamins *are* existing panels; and the substrates, routes, modules and
#: steps are graded over the targets of a panel that is scored. Shown as
#: rows of their own in the overview they had no score to show, which made
#: a third of that section unreadable and repeated five vitamins that were
#: already there with their sliders a few inches above.
#:
#: So each one declares the reading it helps make, and is shown on that
#: reading's detail card instead. A value of ``None`` means the chemistry is
#: defined and no panel measures it yet - urolithin and equol - and those
#: are reported as not searched for rather than as readings.
READING_OWNER: Final[Mapping[str, str | None]] = {
    # A01 - every substrate is graded over the carbohydrate panel's CAZy targets.
    "carb.cellulose": "carbohydrates", "carb.starch_general": "carbohydrates",
    "carb.resistant_starch": "carbohydrates", "carb.chitin": "carbohydrates",
    "carb.pectin": "carbohydrates", "carb.inulin": "carbohydrates",
    "carb.fos": "carbohydrates", "carb.gos": "carbohydrates",
    "carb.xos": "carbohydrates", "carb.imo": "carbohydrates",
    "carb.lactose": "carbohydrates", "carb.beta_glucan": "carbohydrates",
    "carb.arabinoxylan": "carbohydrates", "carb.galactomannan": "carbohydrates",
    # A02 - each route sits with the product it makes, where that is measured.
    "acetate.pta_ack": "fermentation", "acetate.wood_ljungdahl": "fermentation",
    "lactate.d_formation": "fermentation", "lactate.l_formation": "fermentation",
    "lactate.utilisation_butyrate": "butyrate",
    "lactate.utilisation_propionate": "propionate",
    "succinate.formation": "propionate", "succinate.to_propionate": "propionate",
    "propionate.propanediol": "propionate",
    "hydrogen.production_fefe": "fermentation",
    "hydrogen.consumption_methanogenesis": "methane",
    "hydrogen.consumption_sulfate": "h2s", "sulfur.taurine": "h2s",
    # A03 - the four new vitamins have their own panels; the five reused ones
    # are those panels, and stop being repeated.
    "b1": "thiamine", "b3": "niacin", "b5": "pantothenate", "b6": "b6",
    "b2": "riboflavin", "b7": "biotin", "b9": "folate", "b12": "b12", "k2": "k2",
    # A04
    "proteolysis": "nitrogen", "peptide_transport": "peptides",
    "amino_acid_fermentation": "nitrogen", "bcfa": "bcaa",
    "ammonia": "urease", "aromatic": "nitrogen",
    # A05-A08
    "neuro.gaba_synthesis": "gaba", "neuro.gaba_degradation": "gababreak",
    "neuro.gaba_to_butyrate": "gababreak",
    "polyphenol.urolithin_9_dehydroxylation": "urolithin",
    "polyphenol.equol_daidzein_conversion": "equol",
    "diet.glucosinolate_isothiocyanate_conversion": "glucosinolate",
    "mucin.glycan_foraging": "mucinlps", "lps.lipid_a_modification": "mucinlps",
    "bile.deconjugation": "bsh", "bile.dehydroxylation": "bai",
    "bile.hsdh_transformation": "mucinlps",
    # A07
    "glp1.scfa_signalling": "butyrate", "glp1.bile_acid_signalling": "bsh",
    "glp1.indole_signalling": "indole", "glp1.microbial_protein": "p9",
    # §5.8 - each computed over its own scored panel.
    "capacity.urate_anaerobic": "urate", "capacity.putrescine": "putrescine",
    "capacity.agmatine": "agmatine", "capacity.spermidine": "spermidine",
    "capacity.glutathione": "glutathione", "capacity.ethanol": "ethanol",
}


def owner_of(reading_id: str) -> str | None:
    """The scored panel a functional reading is a component of."""
    return READING_OWNER.get(str(reading_id))


def components_by_panel(rows: Sequence[Mapping[str, Any]]) -> dict[str, list[Any]]:
    """Group functional readings under the panel each one helps make."""
    out: dict[str, list[Any]] = {}
    for row in rows:
        panel = owner_of(str(row.get("reading_id") or ""))
        if panel:
            out.setdefault(panel, []).append(row)
    return out


def reading_dest(group: str, label: str, *, detail: bool = False) -> str:
    """The destination name for one functional reading.

    Derived from the group and the label rather than from a running index,
    so a link keeps working when a reading is added above it. §12.3 wants
    every overview reading to reach its own detail card, which needs a name
    per reading and not only per section.
    """
    slug = re.sub(r"[^a-z0-9]+", "-", f"{group}-{label}".lower()).strip("-")
    return f"fr-{slug}" + ("-detail" if detail else "")


def _row(group: str, label: str, found: str, state: str, caveat: str,
         reading_id: str = "") -> dict[str, Any]:
    """One functional reading, as a row.

    ``reading_id`` is what lets the row find the scored panel it is a
    component of, so it can be shown on that panel's card instead of as a
    row with no score of its own.
    """
    return {"group": group, "label": label, "found": found, "state": state,
            "caveat": caveat, "reading_id": reading_id}


def _substrate_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not view:
        return []
    out = []
    for entry in view.get("substrates") or []:
        families = entry.get("required_present") or []
        found = ", ".join(families) if families else DASH
        out.append(_row(
            "Fibre & Dietary Substrates", str(entry.get("label")), found,
            str(entry.get("specificity")), str(entry.get("must_not_conclude") or ""),
            reading_id=str(entry.get("substrate_id") or ""),
        ))
    return out


def _fermentation_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not view:
        return []
    out = []
    for entry in view.get("routes") or []:
        present = entry.get("required_present") or []
        missing = entry.get("required_missing") or []
        found = ", ".join(present) if present else DASH
        if missing and present:
            found += f" <font color='{_hex(INK_FAINT)}'>(missing {', '.join(missing)})</font>"
        out.append(_row(
            "Fermentation & Cross-feeding", str(entry.get("label")), found,
            str(entry.get("state")), str(entry.get("note") or ""),
            reading_id=str(entry.get("route_id") or ""),
        ))
    return out


def _vitamin_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not view:
        return []
    out = []
    for entry in view.get("new") or []:
        columns = entry.get("columns") or {}
        found = "; ".join(
            f"{name}: {STATE_WORD.get(data.get('state'), data.get('state'))}"
            for name, data in columns.items() if data.get("state") != "not_applicable"
        )
        out.append(_row(
            "Vitamins & Nutrients", str(entry.get("label")), found,
            str(columns.get("synthesis", {}).get("state")),
            str(entry.get("must_not_conclude") or ""),
            reading_id=str(entry.get("vitamin_id") or ""),
        ))
    for entry in view.get("reused") or []:
        value = entry.get("existing_value")
        row = _row(
            "Vitamins & Nutrients", str(entry.get("label")),
            (f"{value:,.1f} copies per 100 genomes" if isinstance(value, (int, float))
             else DASH) + " (existing reading, unchanged)",
            "present" if value else "not_assayed", str(entry.get("route_detail") or ""),
            reading_id=str(entry.get("vitamin_id") or ""),
        )
        # This module says of these rows that they are the existing panel
        # reused unchanged, which makes them that panel's row written
        # again. They were printing as "not measured" beside the panel that
        # holds their number, so they are marked duplicates here from the
        # declaration rather than inferred from a gene set they do not have.
        if entry.get("reused_unchanged") and entry.get("panel"):
            row["duplicates_panel"] = str(entry["panel"])
        out.append(row)
    return out


def _nitrogen_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not view:
        return []
    out = []
    for entry in view.get("modules") or []:
        informative = entry.get("informative_genes_present") or []
        found = ", ".join(informative) if informative else DASH
        caveat = "; ".join((entry.get("distinct_from") or {}).values())
        out.append(_row(
            "Protein & Nitrogen Metabolism", str(entry.get("label")), found,
            str(entry.get("state")), caveat,
            reading_id=str(entry.get("module_id") or ""),
        ))
    return out


def _biotransform_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not view:
        return []
    group_for = {
        "A05": "Gut-Brain & Metabolic Signalling",
        "A06": "Plant-Compound Conversion",
        "A07": "Gut-Brain & Metabolic Signalling",
        "A08": "Mucus, Bile & Other Transformations",
    }
    out = []
    for feature in view.get("by_feature") or []:
        group = group_for.get(str(feature.get("feature_id")))
        if group is None:
            continue
        for step in feature.get("steps") or []:
            requirement = step.get("requirement") or {}
            present = requirement.get("present") or []
            found = ", ".join(present) if present else DASH
            out.append(_row(
                group, str(step.get("label")), found, str(step.get("state")),
                str(step.get("does_not_establish") or ""),
                reading_id=str(step.get("step_id") or ""),
            ))
    for channel in (view.get("glp1") or {}).get("channels") or []:
        requirement = channel.get("requirement") or {}
        present = requirement.get("present") or []
        out.append(_row(
            "Gut-Brain & Metabolic Signalling", str(channel.get("label")),
            ", ".join(present) if present else DASH, str(channel.get("state")),
            str(channel.get("does_not_establish") or ""),
            reading_id=str(channel.get("step_id") or ""),
        ))
    return out


def _capacity_rows(view: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Spec 0.8.3 §5.8's additional capacities, as rows of the same table."""
    if not view:
        return []
    out: list[dict[str, Any]] = []
    for capacity in view.get("capacities") or []:
        group = _CAPACITY_GROUP.get(str(capacity.get("capacity_id")))
        if group is None:
            continue
        present = (capacity.get("requirement") or {}).get("present") or []
        copies = capacity.get("copies_per_100_genomes")
        found = ", ".join(present) if present else DASH
        if isinstance(copies, int | float):
            found += f" \u00b7 {copies:,.0f} per 100 genomes"
        out.append(_row(
            group, str(capacity.get("label")), found, str(capacity.get("state")),
            str(capacity.get("does_not_establish") or ""),
            reading_id=str(capacity.get("capacity_id") or ""),
        ))
    return out


def _attach_scores(rows: list[dict[str, Any]],
                   views: Mapping[str, Any] | None,
                   directions: Mapping[str, str] | None = None,
                   aggregates: Mapping[str, Sequence[str]] | None = None,
                   ) -> list[dict[str, Any]]:
    """Give each row its cohort position and the direction to colour it by.

    The direction is inherited from the panel the reading composes, so a
    substrate bar means what the carbohydrate bar above it means. Without
    that the bars are grey and a reader cannot tell a good 90th from a bad
    one - which is the whole purpose of the colour.
    """
    scores = (views or {}).get("reading_scores") or {}
    directions = directions or {}
    for row in rows:
        rid = str(row.get("reading_id") or "")
        score = scores.get(rid) or {}
        row["percentile"] = score.get("percentile")
        row["limiting_gene"] = score.get("limiting_gene")
        row["unscored_reason"] = score.get("unscored_reason")
        owner = owner_of(rid)
        row["higher_means"] = directions.get(str(owner), "unclear") if owner else "unclear"
        # A row may already know it restates a panel, because the module
        # that built it said so. Only work it out where nobody has.
        if not row.get("duplicates_panel"):
            row["duplicates_panel"] = _is_same_quantity_as_panel(score, owner, aggregates)
    return rows


@lru_cache(maxsize=1)
def _sole_reading_of() -> dict[str, str]:
    """Panels whose whole subject is one reading, from the owner map.

    `READING_OWNER` already declares which panel each reading belongs to,
    and the shape of that map answers the question directly: a panel owned
    by exactly one reading *is* that reading's subject, while a panel owned
    by fourteen has fourteen narrower questions under it.

    Vitamin B1 and the thiamine panel are the same vitamin; chitin and the
    carbohydrate panel are not the same thing.
    """
    counts: dict[str, list[str]] = {}
    for reading, panel in READING_OWNER.items():
        counts.setdefault(str(panel), []).append(str(reading))
    return {panel: rs[0] for panel, rs in counts.items() if len(rs) == 1}


def _is_same_quantity_as_panel(
    score: Mapping[str, Any],
    owner: str | None,
    aggregates: Mapping[str, Sequence[str]] | None,
) -> str | None:
    """The panel this reading merely restates, if it restates one.

    Two ways to be the same quantity, and the report needs both.

    The first is declared: a panel with exactly one reading against it has
    that reading as its whole subject. Vitamin B1 read 89th beside a
    thiamine panel reading 76th, and Vitamin B5 77th beside pantothenate at
    58th, because each pair measured one vitamin over two slightly
    different gene lists. Gene-set algebra could not catch those - the
    reading's genes overlapped the panel's aggregate without containing it
    - and the declaration could.

    The second is arithmetic: a reading whose genes are exactly what its
    panel aggregates is that panel written again, whatever the owner map
    says about how many readings share it.
    """
    if not owner:
        return None
    reading_id = str(score.get("reading_id") or "")
    if _sole_reading_of().get(owner) == reading_id and reading_id:
        return owner
    aggregate = set((aggregates or {}).get(owner) or ())
    if not aggregate:
        return None
    keys = {
        str(c.get("key", "")).split(":", 1)[1]
        for c in (score.get("contributors") or [])
        if ":" in str(c.get("key", ""))
    }
    return owner if keys and keys == aggregate else None


def all_functional_rows(
    views: Mapping[str, Any] | None,
    directions: Mapping[str, str] | None = None,
    aggregates: Mapping[str, Sequence[str]] | None = None,
) -> list[dict[str, Any]]:
    """Every A01-A08 and §5.8 reading, as rows, in the §12.1 group order.

    One source for both homes. The at-a-glance section and the detail
    section render the same rows, so a reading cannot appear in one and go
    missing from the other - which is exactly what happened when the glance
    was left to the panel table alone.
    """
    if not views:
        return []
    return _attach_scores([
        *_substrate_rows(views.get("substrates")),
        *_fermentation_rows(views.get("fermentation")),
        *_vitamin_rows(views.get("vitamins")),
        *_nitrogen_rows(views.get("nitrogen")),
        *_biotransform_rows(views.get("biotransformation")),
        *_capacity_rows(views.get("additional_capacities")),
    ], views, directions, aggregates)


def _band_colour(percentile: float | None, higher_means: str) -> Any:
    """The marker colour a panel row would use for this position."""
    from openbiota.pdfreport import status_for  # noqa: PLC0415

    if percentile is None:
        return INK_FAINT
    return status_for(percentile=percentile, higher_means=higher_means).colour


def glance_group_order() -> tuple[str, ...]:
    """The group order, for callers that lay out their own headings."""
    return GLANCE_GROUP_ORDER


def _direction_arrow(higher_means: str) -> str:
    """The same arrow the panel rows use, so the two read alike."""
    from openbiota.pdfreport import CORAL, GREEN  # noqa: PLC0415

    if higher_means == "adverse":
        return f" <font color='{_hex(CORAL)}' size='7'><b>\u25b2</b></font>"
    if higher_means == "favourable":
        return f" <font color='{_hex(GREEN)}' size='7'><b>\u25bc</b></font>"
    return ""


def functional_glance_block(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    views: Mapping[str, Any] | None,
    *,
    detail_section: int,
    aggregates: Mapping[str, Sequence[str]] | None = None,
) -> None:
    """The A01-A08 and §5.8 readings, at a glance.

    Spec 0.8.3 §12.1 gives every one of A01 to A08 two homes: the
    `report.functions` overview and the `report.functions_detail` pages,
    grouped under seven reader-facing subheadings. This is the first of
    those two: one line per reading, what was found and how specific it is.
    What a reading does *not* establish belongs to the detail pages and is
    not repeated here. The number of each section is taken from the caller,
    because those numbers move and these anchors do not.

    The section is allowed to run onto further pages. §12.1 is explicit
    that results are not truncated to fit a page count.
    """
    rows = [r for r in all_functional_rows(views, None, aggregates)
            if not r.get("duplicates_panel")]
    if not rows:
        return
    story.append(CondPageBreak(60 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "Substrates, routes, transformations and capacities", st["h2"]))
    measured = sum(1 for r in rows if r["state"] not in {"not_assayed", "not_applicable"})
    story.append(Paragraph(
        f"{len(rows)} further readings, beyond the panel totals above: what your community "
        "can break down, which routes it carries end to end, and which transformations it "
        f"can perform. {measured} of them were measured in this sample. Each is grouped by "
        "the kind of chemistry, and each says how specific its evidence is, because a gene "
        "family that acts on a whole class of molecules is not evidence about one member of "
        f"it. Section {detail_section} gives each one its full card, including the one "
        "conclusion it does not support.", st["body"],
    ))

    widths = [CONTENT_WIDTH * w for w in (0.34, 0.44, 0.22)]
    for group_label, view_key, _keys in FUNCTION_GROUPS:
        del view_key
        group_rows = [r for r in rows if r["group"] == group_label]
        if not group_rows:
            continue
        story.append(CondPageBreak(30 * mm))
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(f"<b>{group_label}</b>", st["h3"]))
        table_rows: list[list[Any]] = [[
            Paragraph(h, st["label"])
            for h in ("READING", "WHAT WAS FOUND", "HOW SPECIFIC")
        ]]
        for row in group_rows:
            table_rows.append([
                Paragraph(f"<b>{row['label']}</b>", st["cell"]),
                Paragraph(f"<font size='6.4'>{row['found']}</font>", st["cell"]),
                _state_cell(row["state"], st),
            ])
        story.append(_table(table_rows, widths))



def context_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], view: Mapping[str, Any] | None
) -> None:
    """A13's fifteen symptom and condition cards.

    Each card names every reading that feeds it, because a contributor the
    reader cannot see is one they cannot check, and states plainly that the
    number of contributors is a count of things to read rather than a
    likelihood of anything.
    """
    if not view or not view.get("contexts"):
        return
    story.append(CondPageBreak(70 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Topics people ask about, and where the answers are", st["h2"]))
    story.append(Paragraph(str(view.get("what_this_is") or ""), st["body"]))
    story.append(Paragraph(
        f"<font size='6.4' color='{_hex(INK_SOFT)}'><b>{view.get('standing_note')}</b></font>",
        st["small"],
    ))
    for card in view.get("contexts") or []:
        contributing = card.get("contributing") or []
        absent = card.get("not_in_this_run") or []
        story.append(CondPageBreak(42 * mm))
        story.append(Spacer(1, 2.5 * mm))
        story.append(Paragraph(f"<b>{card.get('label')}</b>", st["h3"]))
        story.append(Paragraph(
            f"<font size='6.4'><b>What this keeps separate.</b> {card.get('distinction')}</font>",
            st["small"],
        ))
        for mechanism in card.get("mechanisms") or []:
            story.append(Paragraph(f"<font size='6.3'>{mechanism}</font>", st["small"]))

        rows: list[list[Any]] = [[
            Paragraph(h, st["label"])
            for h in ("IN YOUR REPORT", "WHY IT IS HERE", "WHAT IT SHOWS")
        ]]
        for link in contributing:
            copies = link.get("copies_per_100_genomes")
            state = link.get("state")
            shows = (
                f"{copies:,.0f} per 100 genomes" if isinstance(copies, int | float)
                else (str(state) if state else "linked as calculated")
            )
            rows.append([
                Paragraph(f"<b>{link.get('key')}</b>", st["cell"]),
                Paragraph(f"<font size='6.2'>{link.get('why')}</font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{shows}</font>", st["cell"]),
            ])
        if len(rows) > 1:
            story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.56, 0.22)]))
        story.append(Paragraph(
            f"<font size='6' color='{_hex(INK_FAINT)}'>{card.get('navigation_count_note')}</font>",
            st["small"],
        ))
        if absent:
            story.append(Paragraph(
                f"<font size='6' color='{_hex(INK_FAINT)}'>Not produced in this run, and so "
                "not part of the card: " + ", ".join(str(a.get("key")) for a in absent)
                + ".</font>", st["small"],
            ))
        missing = card.get("missing_information") or []
        if missing:
            story.append(Paragraph(
                "<font size='6.2'><b>What this cannot tell you:</b> "
                + " ".join(missing) + "</font>", st["small"],
            ))
        follow = card.get("follow_up") or []
        if follow:
            story.append(Paragraph(
                f"<font size='6.2' color='{_hex(INK_SOFT)}'><b>If you want to take this "
                f"further:</b> " + " ".join(follow) + "</font>", st["small"],
            ))
    story.append(Spacer(1, 2 * mm))


def drug_metabolism_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], views: Mapping[str, Any] | None
) -> None:
    """Characterised drug reactions, and the substrate detail that is withheld.

    These are cards rather than readings. Whether a drug is transformed
    depends on the strain, the dose and the timing, and none of that is in a
    stool gene count, so every card states what it does not mean and none of
    them supports changing a medication.
    """
    if not views:
        return
    view = views.get("additional_capacities") or {}
    reactions = view.get("drug_reactions") or []
    gus = views.get("gus_substrate_classes") or {}
    if not reactions and not gus:
        return
    story.append(CondPageBreak(60 * mm))
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph("Medicines your gut bacteria are known to alter", st["h2"]))
    story.append(Paragraph(
        "Each of these transformations has been demonstrated with purified enzymes or "
        "defined strains. They are here as mechanism, not as a result about you: none "
        "of them is a measurement of what is happening to a medicine you take, and "
        "nothing in this section is a reason to change a dose or stop a drug. That is a "
        "conversation with the person who prescribed it.", st["body"],
    ))
    table_rows: list[list[Any]] = [[
        Paragraph(h, st["label"])
        for h in ("MEDICINE", "WHAT THE BACTERIA DO", "IN THIS REPORT", "WHAT IT DOES NOT MEAN")
    ]]
    for card in reactions:
        state = str(card.get("state") or "")
        if card.get("measured_by"):
            word = {
                "gene_detected": f"{card.get('gene')} detected",
                "gene_not_detected": f"{card.get('gene')} not detected",
            }.get(state, "measured elsewhere")
            tone = GREEN if state == "gene_detected" else SLATE
        else:
            word, tone = "not measured", INK_FAINT
        table_rows.append([
            Paragraph(f"<b>{card.get('drug')}</b><br/>"
                      f"<font size='5.8' color='{_hex(INK_SOFT)}'>{card.get('organism')}</font>",
                      st["cell"]),
            Paragraph(f"<font size='6.2'>{card.get('mechanism')}</font>", st["cell"]),
            Paragraph(f"<font size='6.2' color='{_hex(tone)}'>{word}</font>", st["cell"]),
            Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{card.get('does_not_mean')}</font>",
                      st["cell"]),
        ])
    story.append(_table(table_rows, [CONTENT_WIDTH * w for w in (0.17, 0.32, 0.13, 0.38)]))

    cards = view.get("neuroactive_cards") or []
    if cards:
        story.append(CondPageBreak(60 * mm))
        story.append(Spacer(1, 3 * mm))
        story.append(Paragraph("Neurotransmitters, and what can actually be measured",
                               st["h2"]))
        story.append(Paragraph(
            "This is the part of the field where the distance between what a stool "
            "sample can show and what people want to conclude is widest, so each of "
            "these says where its evidence stops. Three of the five have no specific "
            "test behind them at all, and say so rather than borrowing a number from "
            "a related molecule.", st["body"],
        ))
        rows: list[list[Any]] = [[
            Paragraph(h, st["label"])
            for h in ("MOLECULE", "WHAT THE EVIDENCE SHOWS", "IN THIS REPORT",
                      "WHERE IT STOPS")
        ]]
        for card in cards:
            state = str(card.get("state") or "")
            if card.get("sequence_panel"):
                word = {
                    "gene_detected": f"{card['sequence_panel']} detected",
                    "gene_not_detected": f"{card['sequence_panel']} not detected",
                }.get(state, f"measured as {card['sequence_panel']}")
                tone = GREEN if state == "gene_detected" else SLATE
                detail = word
            else:
                tone = INK_FAINT
                detail = "no specific test exists"
            rows.append([
                Paragraph(f"<b>{card.get('label')}</b>", st["cell"]),
                Paragraph(f"<font size='6.2'>{card.get('evidence')}</font>", st["cell"]),
                Paragraph(
                    f"<font size='6.2' color='{_hex(tone)}'>{detail}</font>"
                    + (f"<br/><font size='5.8' color='{_hex(INK_FAINT)}'>"
                       f"{card.get('no_panel_because')}</font>"
                       if not card.get("sequence_panel") else ""),
                    st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>"
                          f"{card.get('endpoint_boundary')}</font>", st["cell"]),
            ])
        story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.15, 0.31, 0.24, 0.30)]))

    if gus and not gus.get("available"):
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph(
            f"<font size='6.3'><b>Why there is no breakdown of your beta-glucuronidase "
            f"by substrate.</b> {gus.get('reason')}</font>", st["small"],
        ))
        checks = (gus.get("validation") or {}).get("checks") or []
        failed = [c for c in checks if not c.get("agrees")]
        if failed:
            story.append(Paragraph(
                f"<font size='6' color='{_hex(INK_FAINT)}'>The test that withheld it: "
                + "; ".join(
                    f"{c['organism']} should be {c['expected']}, came out as "
                    + (", ".join(c["observed"]) if c["observed"] else "absent")
                    for c in failed[:4]
                )
                + f". {gus.get('what_would_unlock_it', '')}</font>", st["small"],
            ))
    elif gus.get("available"):
        story.append(Spacer(1, 2 * mm))
        story.append(Paragraph("Your beta-glucuronidase enzymes, by structural class",
                               st["h3"]))
        class_rows: list[list[Any]] = [[
            Paragraph(h, st["label"])
            for h in ("CLASS", "SHARE", "WHAT THIS CLASS ACTS ON", "WHAT IT DOES NOT MEAN")
        ]]
        for row in gus.get("classes") or []:
            if not row.get("fragments"):
                continue
            class_rows.append([
                Paragraph(f"<b>{row.get('label')}</b>", st["cell"]),
                Paragraph(f"<font size='6.4'>{row.get('percent_of_gus')}%</font>", st["cell"]),
                Paragraph(f"<font size='6.2'>{row.get('substrates')}</font>", st["cell"]),
                Paragraph(f"<font size='6.2' color='{_hex(INK_SOFT)}'>{row.get('does_not_mean')}</font>",
                          st["cell"]),
            ])
        story.append(_table(class_rows, [CONTENT_WIDTH * w for w in (0.16, 0.10, 0.36, 0.38)]))
