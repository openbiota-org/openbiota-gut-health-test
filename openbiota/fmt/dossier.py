"""One full page per candidate: what it gives, what it adds, why it ranks there.

The ranked table can only hold a number per column, and numbers like
"15 + 12 routine" or "2 patterns above you, no new organism" are unreadable
without the sentence behind them. This module writes that sentence — a
dossier for every candidate, linked from its row, structured around the
four questions a person actually has to answer before accepting material:

1. **What does it put back that I am missing?** The gaps it fills, each
   with the organism, why that organism is on the list, and how far the
   donor sits above the bar.
2. **What else does it bring that is probably good?** Organisms it carries
   that the recipient lacks and that no disease profile marks as raised.
3. **What does it bring that is neither here nor there?** Normal carriage
   of common organisms: present, unremarkable, not worth a decision.
4. **What could it be introducing?** Organisms a transfer-evidenced disease
   is raised in, each with the donor's amount against the healthy-adult
   range, so "risk" is a number and not an adjective.

Then the arithmetic of its rank, and every gate it failed, in words.

The section deliberately repeats numbers that appear on the summary page.
A reader arriving here from a link should not have to hold the front page
in their head, and a figure that only exists in one place cannot be
checked against anything.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Spacer

from openbiota.pdflinks import Anchor


#: Destination name for a candidate's dossier, derived from its stable ID.
def candidate_dest(candidate_id: str) -> str:
    return f"fmtcand-{candidate_id}"


def _e(text: object) -> str:
    return (
        str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def _organism_buckets(
    body: Mapping[str, Any],
    donor_species: Mapping[str, float],
    recipient_species: Mapping[str, float],
    target_species: frozenset[str],
) -> dict[str, list[dict[str, Any]]]:
    """Sort what this material would add into the four decision buckets.

    The flagged organisms come with their own verdicts already. Everything
    else the donor carries and the recipient does not is split on whether
    it is something the recipient was short of (good), or simply more of
    the ordinary gut (unremarkable).
    """
    flagged = {
        o["species"]: o
        for o in (body.get("donor_health") or {}).get("organisms_it_would_add") or []
    }
    good: list[dict[str, Any]] = []
    neutral: list[dict[str, Any]] = []
    for species, percent in sorted(donor_species.items()):
        if percent <= 0 or recipient_species.get(species, 0.0) > 0:
            continue
        if species in flagged:
            continue
        row = {"species": species, "percent": percent}
        (good if species in target_species else neutral).append(row)
    return {
        "replenishes": [],  # filled by the caller from the goal labels
        "good": good,
        "neutral": neutral,
        "risk": list(flagged.values()),
    }


def dossier_pages(
    story: list[Any],
    st: Mapping[str, Any],
    *,
    ranked: Sequence[Mapping[str, Any]],
    targets: Sequence[Mapping[str, Any]],
    donor_species: Mapping[str, Mapping[str, float]],
    recipient_species: Mapping[str, float],
    leaderboard: Mapping[str, Any],
) -> None:
    """A dossier per candidate, in rank order — every one of them.

    Not a top-N: every row in the ranked table links here, so every row
    needs a page to land on. A reader comparing two options should not
    find that one of them was the cutoff.
    """
    from openbiota.fmt.pdf import _p, _table
    from openbiota.pdfreport import (
        ACCENT_DARK,
        AMBER,
        CORAL,
        GREEN,
        INK_SOFT,
        _hex,
    )

    target_species = frozenset(
        f for t in targets for f in (t.get("feature_ids") or [])
    )
    by_label = {t["label"]: t for t in targets}

    for body in ranked:
        members = " + ".join(body["member_person_ids"])
        health = body.get("donor_health") or {}
        story.append(PageBreak())
        story.append(
            Anchor(
                candidate_dest(body["candidate_id"]),
                # Level 0: this report has no bookmark hierarchy to nest under.
                outline=(f"Candidate {body['match_rank']}: {members}", 0),
            )
        )
        story.append(
            _p(f"Candidate {body['match_rank']}: {_e(members)}", st["h2"])
        )

        # ---- verdict line -------------------------------------------------
        if body["eligible"]:
            verdict = (
                f"<font color='{_hex(ACCENT_DARK)}'><b>Eligible to be "
                "recommended.</b></font> It passes every safety gate: the donor's "
                "own gut-health index is at or above zero, it brings no disease "
                "pattern you do not already have, and it carries no organism above "
                "the healthy-adult range."
            )
        else:
            verdict = (
                f"<font color='{_hex(CORAL)}'><b>Not eligible to be "
                f"recommended.</b></font> {len(body['gate_failures'])} gate"
                f"{'' if len(body['gate_failures']) == 1 else 's'} failed; each is "
                "listed below. It is still ranked and shown so you can see what it "
                "would have offered."
            )
        story.append(_p(verdict, st["small"]))
        story.append(Spacer(1, 2 * mm))

        # ---- the arithmetic ----------------------------------------------
        story.append(
            _p(
                f"<b>Why it ranks {body['match_rank']}.</b> "
                f"{_e(body['suitability_explanation'])}. Candidates are ordered "
                "first by whether they pass the safety gates, then by this score. "
                "Nothing else moves a candidate up.",
                st["small"],
            )
        )
        if body["gate_failures"]:
            for failure in body["gate_failures"]:
                story.append(
                    _p(
                        f"<font size=7 color='{_hex(CORAL)}'>\u2022 "
                        f"{_e(failure)}</font>",
                        st["small"],
                    )
                )
        story.append(Spacer(1, 3 * mm))

        # ---- 1. what it replenishes ---------------------------------------
        filled = list(body["goals"].get("fully_supplied_labels") or [])
        partly = list(body["goals"].get("partly_supplied_labels") or [])
        missed = list(body["goals"].get("not_supplied_labels") or [])
        story.append(
            _p(
                f"<b>1. What it puts back that you are missing</b> "
                f"<font color='{_hex(INK_SOFT)}'>\u2014 {len(filled)} of "
                f"{leaderboard['reachable_goals']} reachable gaps</font>",
                st["body"],
            )
        )
        if filled:
            rows: list[list[Any]] = [["GAP IT FILLS", "WHY THIS IS ON YOUR LIST",
                                      "FOUND BY"]]
            secondary = body["goals"].get("supplied_by_secondary_lane") or {}
            for label in filled:
                t = by_label.get(label, {})
                why = (
                    (t.get("desired_state") or {}).get("statement")
                    or t.get("notes")
                    or "a reading below the reference range in your sample"
                )
                # Which profiler found it. Worth naming: where only the
                # newer lane saw it, the reader should know the older one
                # disagreed rather than assume both agreed.
                found = secondary.get(label)
                lane = (
                    f"<font size=6.2 color='{_hex(AMBER)}'>MetaPhlAn 4 only "
                    f"({found:.3f}%)<br/>MetaPhlAn 3 did not see it</font>"
                    if found is not None
                    else "<font size=6.2>both profilers</font>"
                )
                rows.append([
                    f"<font size=7><b>{_e(label)}</b></font>",
                    f"<font size=6.6>{_e(why)}</font>",
                    lane,
                ])
            story.append(_table(rows, [54 * mm, 78 * mm, 38 * mm], st))
        else:
            story.append(_p("It fills none of your gaps.", st["small"]))
        if partly:
            story.append(
                _p(
                    f"<font size=7 color='{_hex(AMBER)}'>Partly supplied "
                    f"({len(partly)}): {_e(', '.join(partly))}. The donor carries "
                    "the organism but below the level the goal asks for.</font>",
                    st["small"],
                )
            )
        if missed:
            story.append(
                _p(
                    f"<font size=7 color='{_hex(INK_SOFT)}'>Not supplied "
                    f"({len(missed)}): {_e(', '.join(missed))}.</font>",
                    st["small"],
                )
            )
        story.append(Spacer(1, 3 * mm))

        # ---- buckets 2-4 ---------------------------------------------------
        buckets = _organism_buckets(
            body,
            {
                k: v
                for person in body["member_person_ids"]
                for k, v in (donor_species.get(person) or {}).items()
            },
            recipient_species,
            target_species,
        )

        story.append(
            _p(
                f"<b>2. What else it brings that is probably good</b> "
                f"<font color='{_hex(INK_SOFT)}'>\u2014 {len(buckets['good'])} "
                "organisms</font>",
                st["body"],
            )
        )
        if buckets["good"]:
            story.append(
                _p(
                    "<font size=7>Organisms you do not carry that are on your own "
                    "goal list as things to restore: "
                    + ", ".join(
                        f"<i>{_e(r['species'].replace('_', ' '))}</i> "
                        f"({r['percent']:.3f}%)"
                        for r in buckets["good"][:24]
                    )
                    + (f", and {len(buckets['good']) - 24} more"
                       if len(buckets["good"]) > 24 else "")
                    + ".</font>",
                    st["small"],
                )
            )
        else:
            story.append(
                _p(
                    "<font size=7>Nothing beyond the gaps listed above.</font>",
                    st["small"],
                )
            )
        story.append(Spacer(1, 2.5 * mm))

        story.append(
            _p(
                f"<b>3. What it brings that is neither here nor there</b> "
                f"<font color='{_hex(INK_SOFT)}'>\u2014 {len(buckets['neutral'])} "
                "organisms</font>",
                st["body"],
            )
        )
        story.append(
            _p(
                f"<font size=7 color='{_hex(INK_SOFT)}'>Ordinary gut organisms the "
                "donor has and you do not, which are not on your goal list and are "
                "not raised in any disease with transfer evidence. Two people never "
                "carry the same species list, and most of this difference is "
                "unremarkable. "
                + (
                    ", ".join(
                        f"<i>{_e(r['species'].replace('_', ' '))}</i>"
                        for r in buckets["neutral"][:20]
                    )
                    + (f", and {len(buckets['neutral']) - 20} more."
                       if len(buckets["neutral"]) > 20 else ".")
                    if buckets["neutral"] else "None."
                )
                + "</font>",
                st["small"],
            )
        )
        story.append(Spacer(1, 2.5 * mm))

        # ---- 4. risk --------------------------------------------------------
        risk = buckets["risk"]
        abnormal = [r for r in risk if r["abnormal"] or r["is_known_pathogen"]]
        story.append(
            _p(
                f"<b>4. What it could be introducing</b> "
                f"<font color='{_hex(CORAL if abnormal else INK_SOFT)}'>\u2014 "
                f"{len(risk)} organism{'' if len(risk) == 1 else 's'} flagged, "
                f"{len(abnormal)} above the healthy range</font>",
                st["body"],
            )
        )
        if risk:
            rows = [["ORGANISM", "IN THIS DONOR", "HEALTHY ADULTS WHO CARRY IT",
                     "VERDICT"]]
            for r in risk:
                bad = r["abnormal"] or r["is_known_pathogen"]
                ref = (
                    "not in the reference"
                    if r["reference_p90_percent"] is None
                    else (
                        f"median {r['reference_median_percent']:.3f}%, "
                        f"90th {r['reference_p90_percent']:.3f}%, "
                        f"{r['carried_by_percent_of_healthy_adults']:.0f}% carry it"
                    )
                )
                rows.append([
                    f"<font size=7><i>{_e(r['species'].replace('_', ' '))}</i><br/>"
                    f"<font size=6 color='{_hex(INK_SOFT)}'>raised in "
                    f"{_e(', '.join(r['raised_in_conditions']))}</font></font>",
                    f"<font size=7><b>{r['donor_percent']:.3f}%</b></font>",
                    f"<font size=6.4>{_e(ref)}</font>",
                    f"<font size=6.8 color='{_hex(CORAL if bad else GREEN)}'>"
                    f"<b>{_e(r['verdict'])}</b></font>",
                ])
            story.append(_table(rows, [54 * mm, 22 * mm, 60 * mm, 34 * mm], st))
            for r in risk:
                story.append(
                    _p(f"<font size=6.6 color='{_hex(INK_SOFT)}'>{_e(r['plain'])}"
                       "</font>", st["small"])
                )
        else:
            story.append(
                _p(
                    f"<font size=7 color='{_hex(GREEN)}'>Nothing. This material "
                    "carries no organism that a transfer-evidenced disease is "
                    "raised in and that you do not already have.</font>",
                    st["small"],
                )
            )
        story.append(Spacer(1, 3 * mm))

        # ---- the row's own numbers, decoded --------------------------------
        story.append(_p("<b>The numbers on this candidate's row, decoded</b>",
                        st["body"]))
        decoded = [
            (
                f"Match score {body['match_score']:.1f}",
                f"It supplies {body['goals']['covered']} of the "
                f"{leaderboard['reachable_goals']} gaps that any material in this "
                "comparison could supply. Not a probability of anything working.",
            ),
            (
                f"Findings to test: {body['serious_open_questions']} + "
                f"{body['other_findings']} routine",
                f"{body['serious_open_questions']} organism(s) the sequencing found "
                "that it cannot settle from DNA alone and that would matter if the "
                "evidence held up \u2014 each needs a laboratory test on the banked "
                f"sample to resolve. The other {body['other_findings']} are the "
                "routine residue every stool sample produces: species a close "
                "relative makes unnameable, organisms carried without the genes that "
                "cause disease, resistance genes with no resolved host. Neither "
                "number moves the score; they are a tie-break only.",
            ),
            (
                f"People involved: {body['unique_donor_count']}",
                "Each additional donor is another screening history, another "
                "consent, and another set of unknowns. At equal scores fewer is "
                "ranked higher.",
            ),
            (
                f"Donor gut health {health.get('gmwi2'):+.2f}"
                if health.get("gmwi2") is not None else "Donor gut health —",
                "GMWI2, a published index trained on 8,069 labelled stool samples. "
                "Zero is its neutral point; below zero the community reads as "
                "dysbiotic. For a set this is the worst member's score, not an "
                "average \u2014 pooling does not make an unwell community well.",
            ),
        ]
        rows = [["WHAT THE ROW SAYS", "WHAT IT MEANS"]]
        for term, meaning in decoded:
            rows.append([
                f"<font size=7><b>{_e(term)}</b></font>",
                f"<font size=6.6>{_e(meaning)}</font>",
            ])
        story.append(_table(rows, [44 * mm, 126 * mm], st))
