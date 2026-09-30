"""Front-loaded HTML and PDF renderers (spec §13).

Both formats read the canonical JSON only. Neither recomputes a score, so the
three outputs cannot disagree (V7-097). Section order follows §13.1, with one
addition the spec asks for in spirit and a first reader needs in practice: the
report defines its own terms before it uses them, leads with a plain count of
how many of the recipient's gaps a candidate covers, states the highest index
anything could reach, and shows the facts it deliberately keeps out of the
score next to every candidate with the reason.

Status is always carried in words as well as colour, and no green shield, risk
percentage or combined benefit/safety dial exists anywhere (§13.2).
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

#: Plain-language display labels for the concern tiers. The machine-readable
#: tier stays in the JSON; a reader gets a sentence.
TIER_LABELS: dict[str, str] = {
    "confirmed_exclusion": "confirmed by a laboratory test — excluded",
    "high_consequence_supported": "serious organism, sequence support — needs a lab test",
    "high_consequence_unresolved": "serious organism, identity not confirmed",
    "toxin_or_resistance_review": "toxin or resistance gene, carrier organism unknown",
    "screening_gap": "no screening records supplied",
    "species_unresolved": "DNA present, species cannot be named",
    "carriage_toxin_negative": "carried, without the genes that cause disease",
    "technical_ambiguity": "no organism-specific evidence",
    "context": "context only",
}

LEAD_TIERS = (
    "confirmed_exclusion",
    "high_consequence_supported",
    "high_consequence_unresolved",
    "species_unresolved",
    "carriage_toxin_negative",
    "toxin_or_resistance_review",
    "screening_gap",
)

_CSS = """
:root { --paper:#f6f6ef; --ink:#193e30; --muted:#61765b; --line:#d7dfd0; --forest:#092b22;
        --leaf:#668857; --amber:#9a6b12; --slate:#5b6b76; }
* { box-sizing:border-box; }
body { margin:0; background:var(--paper); color:var(--ink);
       font:400 15px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif; }
.wrap { max-width:1180px; margin:0 auto; padding:28px 24px 80px; }
h1 { font:400 2.2rem/1.15 Georgia,serif; letter-spacing:-.02em; margin:0 0 6px; }
h2 { font:400 1.55rem/1.2 Georgia,serif; margin:42px 0 6px; padding-top:18px; border-top:1px solid var(--line); }
h3 { font-size:1.02rem; margin:22px 0 6px; }
p { margin:0 0 12px; max-width:76ch; }
.eyebrow { font:600 11px/1.5 inherit; letter-spacing:.14em; text-transform:uppercase; color:var(--muted); }
.lede { font:400 1.15rem/1.5 Georgia,serif; max-width:78ch; }
.answer { background:#fff; border:1px solid var(--line); border-left:4px solid var(--leaf);
          border-radius:6px; padding:16px 20px; margin:18px 0 6px; max-width:88ch; }
.winner { display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:0;
          background:#f3f7f2; border:1px solid var(--leaf); border-radius:8px; margin:20px 0 10px; }
.winner > div { padding:14px 18px; border-right:1px solid #dfe7db; }
.winner > div:last-child { border-right:0; }
.winner .k { font:600 10.5px/1.5 inherit; letter-spacing:.12em; text-transform:uppercase;
             color:#3f6b50; }
.winner .v { font:400 1.85rem/1.15 Georgia,serif; letter-spacing:-.02em; margin-top:2px; }
.winner .v small { font:400 .82rem/1 inherit; color:var(--muted); letter-spacing:0; }
.winner .amber .v { color:var(--amber); }
.winner .green .v { color:var(--leaf); }
h3.ok, h3.warn { font:600 .95rem/1.3 inherit; letter-spacing:.01em;
                 margin:26px 0 8px; }
h3.ok, h3.warn { padding-left:10px; border-left:3px solid currentColor; }
h3.ok { color:var(--forest); }
h3.warn { color:var(--amber); }
h3.ok small, h3.warn small { font:400 .8rem/1 inherit; color:var(--muted);
                             margin-left:6px; }
p.alert { border-left:3px solid var(--amber); background:#fbf4e2; padding:10px 14px;
          margin:12px 0; }
.tag.ok { border-color:#c3d6b8; background:#eef4ea; color:var(--forest); }
tr.rank1 td { background:#f3f7f2; font-weight:600; }
td.why { font-size:12px; color:var(--muted); }
.answer strong { font-weight:600; }
table { border-collapse:collapse; width:100%; margin:14px 0 6px; font-size:13.5px; }
th,td { text-align:left; padding:7px 10px; border-bottom:1px solid var(--line); vertical-align:top; }
th { font-size:11px; letter-spacing:.08em; text-transform:uppercase; color:var(--muted); font-weight:600; }
td.num,th.num { text-align:right; font-variant-numeric:tabular-nums; }
tr.lead td { background:#fbfcf8; }
.bar { position:relative; height:9px; background:#e6e9de; border-radius:2px; overflow:hidden; min-width:120px; }
.bar i { position:absolute; top:0; bottom:0; left:0; background:var(--leaf); }
.bar u { position:absolute; top:0; bottom:0; width:2px; background:#9aa891; }
.tag { display:inline-block; font-size:11px; letter-spacing:.04em; padding:2px 7px; border:1px solid var(--line);
       border-radius:3px; background:#fff; color:var(--slate); white-space:nowrap; }
.tag.review { border-color:#e0cda0; background:#fbf4e2; color:var(--amber); }
.tag.excluded { border-color:#d8b3ad; background:#f9ecea; color:#8d3b2f; }
.tag.context { border-color:var(--line); background:#fff; color:var(--muted); }
.card { background:#fff; border:1px solid var(--line); border-radius:8px; padding:16px 18px; margin:12px 0; }
.grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:12px; }
.kv { font-size:13px; color:var(--muted); }
.kv strong { color:var(--ink); font-weight:600; display:block; font-size:15px; margin-top:2px; }
dl.terms { max-width:88ch; margin:10px 0 0; }
dl.terms dt { font-weight:600; margin-top:14px; }
dl.terms dd { margin:3px 0 0; color:var(--muted); font-size:14px; }
code { font:12.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace; background:#eef0e8; padding:1px 5px; border-radius:3px; }
.note { font-size:12.5px; color:var(--muted); max-width:84ch; }
.warn { border-left:3px solid #c9a227; padding-left:12px; }
footer { margin-top:46px; border-top:1px solid var(--line); padding-top:14px; font-size:12.5px; color:var(--muted); }
"""


def _e(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _num(value: Any, digits: int = 1, dash: str = "—") -> str:
    if value is None:
        return dash
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return _e(value)


def _bar(point: Any, ceiling: Any) -> str:
    """Filled to the value, with a tick at the achievable ceiling."""
    if point is None:
        return '<span class="note">no coverage index</span>'
    tick = ""
    if isinstance(ceiling, (int, float)) and ceiling > 0:
        tick = f'<u style="left:{min(100.0, float(ceiling)):.2f}%"></u>'
    return f'<div class="bar">{tick}<i style="width:{float(point):.2f}%"></i></div>'


def render_html(result: dict[str, Any]) -> str:  # noqa: PLR0915 - one long document builder
    r = result
    rec = r["recipient"]
    ach = r["achievable"]
    ceiling = ach["index_ceiling"]
    donor_names = {row["candidate_id"]: row["person_id"] for row in r["individual_results"]}
    ind = (rec.get("indication") or {}).get("code") if isinstance(rec.get("indication"), dict) else None
    ts = r["target_summary"]

    out: list[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>FMT donor matching — recipient {_e(rec['person_id'])}</title>",
        f"<style>{_CSS}</style></head><body><div class='wrap'>",
        "<p class='eyebrow'>OpenBiota · experimental computational donor matching</p>",
        f"<h1>Recipient {_e(rec['person_id'])} against {len(r['materials'])} candidate material(s)</h1>",
    ]

    # ---- 1. the answer: one winner, then the ranked field ----------------- #
    lb = r["leaderboard"]
    ranked = sorted(r["set_results"], key=lambda b: b["match_rank"])
    top = ranked[0] if ranked else None
    recommended = next(
        (b for b in ranked if b["candidate_id"] == lb.get("recommended_id")), None
    )
    winner = recommended or top
    if winner is not None:
        obs = winner.get("observed_pattern_count", 0)
        out.append("<div class='winner'>")
        for cls, key, value, sub in (
            ("", "Recommended" if recommended is not None else "Least introduced",
             " + ".join(winner["member_person_ids"]), ""),
            ("", "Match score", _num(winner["match_score"]), "of 100"),
            ("", "Gaps supplied", str(winner["goals"]["covered"]),
             f"of {lb['reachable_goals']} reachable"),
            # Was `new_pattern_count`, which is structurally zero since the
            # percentile bridge was retired, so it printed a green "0
            # introduced" for every candidate whatever it carried.
            ("amber" if obs else "green", "Organisms it would add",
             str(obs), "context, not a score" if obs else "none at the threshold"),
        ):
            out.append(
                f"<div class='{cls}'><div class='k'>{_e(key)}</div>"
                f"<div class='v'>{_e(value)}" + (f" <small>{_e(sub)}</small>" if sub else "")
                + "</div></div>"
            )
        out.append("</div>")
        if lb.get("recommended_note"):
            out.append(
                "<p><strong>Why this one and not the top of the table.</strong> "
                f"{_e(lb['recommended_note'])}</p>"
            )
        if lb.get("no_clean_candidate_note"):
            out.append(
                f"<p class='alert'><strong>{_e(lb['no_clean_candidate_note'])}</strong></p>"
            )
        out.append(f"<p>{_e(winner['rank_reason'])}</p>")
        if top is not None and recommended is not None and top["candidate_id"] != winner["candidate_id"]:
            intro = ", ".join(
                f"{e['short_label']} ({e['donor_percentile']:.0f}th vs your "
                f"{e['recipient_percentile']:.0f}th)"
                for e in top["introduces_new_patterns"]
            )
            out.append(
                f"<p><strong>Highest score:</strong> {_e(' + '.join(top['member_person_ids']))} "
                f"at {_num(top['match_score'])}, supplying {top['goals']['covered']} of "
                f"{lb['reachable_goals']} gaps — but it would introduce {_e(intro)}.</p>"
            )
        if lb.get("best_single_donor_note"):
            out.append(f"<p class='note'>{_e(lb['best_single_donor_note'])}</p>")
        if lb.get("tie_note"):
            out.append(f"<p class='note'>{_e(lb['tie_note'])}</p>")
        if ach["unreachable_goals"]:
            out.append(
                "<p class='note'><strong>No candidate supplies:</strong> "
                + _e("; ".join(u["label"] for u in ach["unreachable_goals"]))
                + f". Those {len(ach['unreachable_goals'])} gap(s) sit outside every score here, "
                "which is why the scale runs against what is reachable.</p>"
            )
        out.append(
            "<p class='note'>A match score is restoration evidence only: the share of the "
            "recipient's reachable gaps a candidate can supply. It is not a probability that "
            "anything will establish, not a probability of benefit, and not a safety ranking — "
            "findings are listed beside it and never traded against it.</p>"
        )

    out.append("<h2>Every candidate, ranked</h2>")
    out.append(
        "<p class='note'>Ranked on restoration score, after the patterns a candidate "
        "would bring. A candidate carrying a pattern you do not have, for a condition "
        "where patient stool has reproduced disease features in recipient animals, is "
        "ranked below every candidate carrying none &mdash; whatever it scores &mdash; "
        "and is never the recommendation. Weaker associations cost rank without "
        "vetoing. The <strong>#</strong> column is the single ranking after all of "
        "that.</p>"
    )

    def _rows(bodies: list[dict[str, Any]], *, clean: bool) -> None:
        if not bodies:
            return
        title = (
            "Carries no pattern you do not have"
            if clean else "Carries a pattern you do not have — context, not an exclusion"
        )
        out.append(
            f"<h3 class=\"{'ok' if clean else 'warn'}\">{title} "
            f"<small>{len(bodies)} of {len(ranked)} candidates</small></h3>"
        )
        out.append(
            "<table><thead><tr><th>#</th><th>Candidate</th><th class='num'>Match score</th>"
            "<th class='num'>Gaps supplied</th><th class='num'>People</th>"
            "<th>Organisms it would add</th><th class='num'>Findings to test</th>"
            "<th>Why it ranks here</th></tr></thead><tbody>"
        )
        for b in bodies:
            partly = (
                f" <span class='note'>({b['goals']['partly_supplied']} partly)</span>"
                if b["goals"]["partly_supplied"] else ""
            )
            names = [x["short_label"] for x in b["introduces_new_patterns"]]
            cell = (
                "<span class='tag ok'>none</span>" if not names
                else "".join(f"<span class='tag review'>{_e(n)}</span> " for n in names)
            )
            out.append(
                f"<tr class=\"{'rank1' if b.get('recommended') else ''}\">"
                f"<td>{b['match_rank']}</td>"
                f"<td>{_e(' + '.join(b['member_person_ids']))}</td>"
                f"<td class='num'>{_num(b['match_score'])}</td>"
                f"<td class='num'>{b['goals']['covered']} of {lb['reachable_goals']}{partly}</td>"
                f"<td class='num'>{b['unique_donor_count']}</td>"
                f"<td>{cell}</td>"
                f"<td class='num'><span class='tag review'>"
                f"{b['serious_open_questions']}</span>"
                f" <span class='note'>+ {b['other_findings']} routine</span></td>"
                f"<td class='why'>{_e(b.get('rank_reason_short') or '')}</td></tr>"
            )
        out.append("</tbody></table>")

    # Split on what a candidate actually carries. `new_pattern_count` is
    # zero for everyone while the percentile bridge stays retired, so the
    # old split filed all candidates under "introduces no new disease
    # pattern" and then listed the patterns they carry two sections later.
    # Clean group first: rank 1 is in it.
    _rows([b for b in ranked if not b.get("observed_pattern_count")], clean=True)
    _rows([b for b in ranked if b.get("observed_pattern_count")], clean=False)

    # ---- the evidence behind every flagged pattern ------------------------ #
    exposures = lb.get("new_pattern_exposures") or []
    if exposures:
        out.append("<h2>Patterns a candidate would bring that you do not have</h2>")
        out.append(
            "<p class='note'>Each row is a pattern where a candidate sits at or above the 75th "
            "percentile and you sit at least 20 points lower. The evidence column is the strongest "
            "published transfer experiment for that condition, with its limit stated. None of this "
            "is a diagnosis of the donor and none of it is a probability for you.</p>"
        )
        out.append(
            "<table><thead><tr><th>Pattern</th><th>Candidates</th><th class='num'>Them vs you"
            "</th><th>Transfer evidence</th><th>Evidence tier</th></tr></thead><tbody>"
        )
        for e in exposures:
            # The tier itself. This used to read "no — human donor
            # transfer", which parses as though the strongest evidence tier
            # were the reason the pattern did not count.
            tier = _e(e["transfer_evidence_tier"].replace("_", " "))
            against = (
                f"<span class='tag review'>{tier}</span>"
                if e["transfer_evidence_tier"] == "human_donor_transfer"
                else f"<span class='note'>{tier}</span>"
            )
            src = (
                f" <span class='note'>Sources: {_e('; '.join(e['sources']))}.</span>"
                if e["sources"] else ""
            )
            out.append(
                f"<tr><td><strong>{_e(e['short_label'])}</strong></td>"
                f"<td>{_e(', '.join(e['donor_person_ids']))}</td>"
                f"<td class='num'><strong>{e['donor_percentile']:.0f}th</strong> vs your "
                f"{e['recipient_percentile']:.0f}th</td>"
                f"<td>{_e(e['evidence'])} <span class='note'>Limit: "
                f"{_e(e['evidence_boundary'])}</span>{src}</td>"
                f"<td>{against}</td></tr>"
            )
        out.append("</tbody></table>")
        out.append(
            "<p class='note'><strong>What this does to the ranking.</strong> A "
            "resemblance percentile is not a diagnosis and is not itself a transmissible "
            "thing, and it never changes the match score. It does change the order. For "
            "a handful of conditions, stool from patients has been put into animals and "
            "reproduced disease features against healthy-donor controls &mdash; "
            "colorectal cancer most clearly of all. A candidate carrying one of those, "
            "where you do not, is <strong>never recommended</strong> and is ranked below "
            "every candidate that carries none, whatever it scores. Patterns with "
            "weaker, association-only evidence do not veto; they cost rank and are "
            "listed here so you can see what you are accepting.</p>"
        )

    out.append("<h2>What each column means</h2><dl class='terms'>")
    for term, text in (
        ("Match score",
         "The share of the recipient's reachable gaps this candidate supplies, each gap counted "
         f"once, with partial credit where supply is graded. {lb['reachable_goals']} of "
         f"{lb['scored_goals']} scored gaps are reachable — the rest are supplied by no candidate "
         "here — so 100 means \u201csupplies everything these materials could supply between "
         "them\u201d."),
        ("Gaps supplied",
         "The plain count behind the score. A gap is an organism or gene capacity the recipient is "
         "missing or low in, derived from the recipient's own data before any donor was examined. "
         "\u201cPartly\u201d means present but still below the reference range."),
        ("Which profiler answered",
         "Two taxonomic profilers run on every sample. Presence goals \u2014 \u201cis this organism "
         "here at all\u201d \u2014 are answered by either, because that question needs no reference "
         "population and MetaPhlAn 4 is four years newer and built from roughly a million genomes "
         "against MetaPhlAn 3's seventeen thousand. Every percentile in this report stays on "
         "MetaPhlAn 3, because curatedMetagenomicData is the only reference cohort that exists and "
         "it is MetaPhlAn 3; its maintainers state the two are not directly comparable. The PDF "
         "names the profiler for each gap, and match.json records it under "
         "\u201csupplied_by_secondary_lane\u201d."),
        ("People",
         "How many different donors the candidate involves. At equal scores fewer ranks higher: one "
         "person means one screening history, one consent and one set of unknowns."),
        ("Serious findings",
         "Organisms that would matter if the evidence held up, and that a laboratory test would "
         "settle. None here is confirmed. This column never moves a candidate up or down the score "
         "— it breaks ties only, and every finding is reported in full below."),
        ("Other findings",
         "The routine residue every stool sample produces: species a close relative makes "
         "unnameable, organisms carried without the genes that cause disease, resistance genes with "
         "no resolved carrier."),
        ("Organisms it would add",
         "A pattern this candidate scores high on that you do not — at or above the 75th "
         "percentile, and at least 20 percentile points above your own reading — for a condition "
         "where stool from patients has reproduced disease features in recipient animals. A "
         "pattern you already carry at a similar level is not counted: it is not new to you. "
         "These are resemblance scores, not diagnoses, and every transfer experiment behind them "
         "used an animal, usually a predisposed one. They are listed because \u201cdo not give me "
         "a condition I do not already have\u201d is a reasonable thing to ask of a donor choice."),
        ("Findings to test",
         "Organisms the sequencing found that it cannot settle. The first number is the serious "
         "ones \u2014 an organism that would matter if the evidence held up. The second is the "
         "routine residue every stool sample produces. Each finding names the test that would "
         "settle it and what the test needs: a nucleic-acid test works on a thawed aliquot of the "
         "banked sample, while a culture needs viable organisms and should be run on a fresh "
         "sample from the donor."),
        ("Identity not confirmed",
         "The sequence matched an organism, but not well enough to name the species: either the "
         "support sits in a region the organism shares with a close relative, or the organism is "
         "the same genomic species as a relative most healthy guts carry. No species name is "
         "given because it would be a guess. Nothing has been ruled in or out."),
        ("Ranking rule", lb["ranking_rule"] + " Disease-pattern resemblance does not enter "
         "this rule at all: it neither changes a score nor excludes a candidate. Patterns a "
         "candidate carries that you do not are listed separately so they can be weighed."),
    ):
        out.append(f"<dt>{_e(term)}</dt><dd>{_e(text)}</dd>")
    out.append("</dl>")

    out.append("<div class='grid'>")
    for label, value in (
        ("Recipient sample", rec["sample_id"]),
        ("Indication", ind or "not declared"),
        ("Scored gaps", f"{ts['positive_scoring']} of {ts['total']} ledger entries"),
        ("Reachable by these materials", f"{lb['reachable_goals']} gaps"),
        ("Candidates evaluated", r["search_audit"]["evaluated_subsets"]),
        ("Confirmed exclusions", len(r["excluded_candidates"])),
        ("Clinical release status", r["clinical_release_status"].replace("_", " ")),
        ("Software", f"openbiota {r['software']['version']}"),
    ):
        out.append(f"<div class='card kv'><div>{_e(label)}</div><strong>{_e(value)}</strong></div>")
    out.append("</div>")

    # ---- 2. how to read this --------------------------------------------- #
    out.append("<h2>How to read this report</h2>")
    out.append("<dl class='terms'>")
    for item in r["how_to_read"]:
        out.append(f"<dt>{_e(item['term'])}</dt><dd>{_e(item['plain'])}</dd>")
    out.append("</dl>")

    # ---- 3. individual candidates ---------------------------------------- #
    out.append("<h2>Every individual candidate</h2>")
    out.append(
        "<p class='note'>Read the first column. The two index columns are the same information "
        "weighted two ways, and are shown together so a candidate covering the same number of gaps "
        "as another never appears to score differently for no visible reason.</p>"
    )
    out.append(
        "<table><thead><tr><th>Donor</th><th class='num'>Gaps covered</th>"
        f"<th class='num'>Index (equal groups, max {_num(ceiling, 0)})</th>"
        "<th class='num'>Index (weighted by group size)</th><th>Scenario range</th>"
        "<th class='num'>Needs a lab test</th><th class='num'>Other open questions</th>"
        "<th class='num'>GMWI2</th><th>Screening</th></tr></thead><tbody>"
    )
    concerns_by_cand: dict[str, list[dict[str, Any]]] = {}
    for c in r["concern_records"]:
        concerns_by_cand.setdefault(c["candidate_id"], []).append(c)
    for row in sorted(r["individual_results"], key=lambda x: -x["goals"]["covered"]):
        cov, alt, g = row["coverage"], row["coverage_alternatives"], row["goals"]
        mine = concerns_by_cand.get(row["candidate_id"], [])
        serious = sum(1 for c in mine if c["tier"] in LEAD_TIERS and c["tier"] != "screening_gap")
        other = len(mine) - serious
        gm = row["context_not_in_the_score"]["gmwi2"]
        out.append(
            f"<tr><td><strong>{_e(row['person_id'])}</strong>"
            f"<div class='note'><code>{_e(row['sample_id'])}</code></div></td>"
            f"<td class='num'><strong>{g['covered']} of {g['total']}</strong></td>"
            f"<td class='num'>{_num(cov['coverage'])}</td>"
            f"<td class='num'>{_num(alt['goal_groups_weighted_by_size'])}</td>"
            f"<td>{_bar(cov['coverage'], ceiling)}"
            f"<span class='note'>{_num(cov['scenario_lower'])}–{_num(cov['scenario_upper'])}</span></td>"
            f"<td class='num'><span class='tag review'>{serious}</span></td>"
            f"<td class='num'>{other}</td>"
            f"<td class='num'>{_num(gm['score'], 2)}</td>"
            f"<td><span class='tag context'>{_e(row['screening_status'].replace('_', ' '))}</span></td></tr>"
        )
    out.append("</tbody></table>")

    for row in sorted(r["individual_results"], key=lambda x: -x["goals"]["covered"]):
        g = row["goals"]
        singleton = next(
            (s for s in r["set_results"] if s["candidate_id"] == row["candidate_id"]), None
        )
        summary = row["plain_summary"]
        if summary.startswith(row["person_id"]):
            summary = summary[len(row["person_id"]) :].lstrip()
        out.append(
            f"<div class='card'><strong>{_e(row['person_id'])}</strong> {_e(summary)}"
        )
        if g["partly_supplied_labels"]:
            out.append(
                f"<div class='note'>Partly covered (present but below the reference range, or a gene "
                f"capacity short of it): {_e('; '.join(g['partly_supplied_labels']))}</div>"
            )
        if singleton:
            uniq = singleton["unique_target_ids_by_member"].get(row["candidate_id"]) or []
            if uniq:
                labels = {t["target_id"]: t["label"] for t in r["targets"]}
                out.append(
                    f"<div class='note'>Covers: {_e('; '.join(labels[t] for t in uniq))}</div>"
                )
        only_here = [
            t["label"]
            for t in r["targets"]
            if t["in_positive_denominator"]
            and sum(
                1
                for other in r["individual_results"]
                for cell in other["per_target_availability"]
                if cell["target_id"] == t["target_id"] and (cell["supported_value"] or 0) > 0
            )
            == 1
            and any(
                cell["target_id"] == t["target_id"] and (cell["supported_value"] or 0) > 0
                for cell in row["per_target_availability"]
            )
        ]
        if only_here:
            out.append(
                f"<div class='note'><strong>No other candidate has:</strong> "
                f"{_e('; '.join(only_here))}</div>"
            )
        out.append("</div>")

    for row in r["excluded_candidates"]:
        out.append(
            f"<div class='card warn'><span class='tag excluded'>excluded on confirmed evidence</span> "
            f"<strong>{_e(row['person_id'])}</strong> — full analytics retained in match.json.</div>"
        )

    # ---- 4. combinations -------------------------------------------------- #
    out.append("<h2>Does a combination help?</h2>")
    out.append(
        "<p class='note'>A combination's coverage is a ceiling: at least one member has evidence of "
        "supplying the gap. It is not a pooled concentration, not a prediction of the community a "
        "physical combination would produce, and not a preparation instruction. A combination also "
        "inherits every member's open pathogen questions — they never average out.</p>"
    )
    out.append(
        "<table><thead><tr><th>Members</th><th class='num'>Gaps covered</th>"
        "<th class='num'>Index</th><th class='num'>Gain vs best single</th>"
        "<th class='num'>Needs a lab test</th><th>Front</th><th class='num'>Stable across scenarios</th>"
        "</tr></thead><tbody>"
    )
    for s in r["set_results"]:
        mine = [c for c in r["concern_records"] if c["candidate_id"] in s["member_material_ids"]]
        serious = sum(1 for c in mine if c["tier"] in LEAD_TIERS and c["tier"] != "screening_gap")
        cls = " class='lead'" if s["pareto_front"] == 1 else ""
        out.append(
            f"<tr{cls}><td>{_e(' + '.join(s['member_person_ids']))}</td>"
            f"<td class='num'><strong>{s['goals']['covered']} of {s['goals']['total']}</strong></td>"
            f"<td class='num'>{_num(s['coverage']['coverage'])}</td>"
            f"<td class='num'>{_num(s['marginals']['gain_over_global_best_singleton'])}</td>"
            f"<td class='num'><span class='tag review'>{serious}</span></td>"
            f"<td>{_e(s['pareto_front'])}</td>"
            f"<td class='num'>{s['stability']['leading_front_scenarios']}/{s['stability']['scenario_count']}</td></tr>"
        )
    out.append("</tbody></table>")
    indispensable = {
        t["target_id"]: [
            other["person_id"]
            for other in r["individual_results"]
            for cell in other["per_target_availability"]
            if cell["target_id"] == t["target_id"] and (cell["supported_value"] or 0) > 0
        ]
        for t in r["targets"]
        if t["in_positive_denominator"]
    }
    solo = {tid: who for tid, who in indispensable.items() if len(who) == 1}
    labels_all = {t["target_id"]: t["label"] for t in r["targets"]}
    if solo:
        out.append(
            "<p class='note'><strong>Gaps that depend on one candidate:</strong> "
            + _e("; ".join(f"{labels_all[tid]} — only {who[0]}" for tid, who in solo.items()))
            + ".</p>"
        )
    else:
        out.append(
            "<p class='note'><strong>No gap here depends on a single candidate:</strong> every gap that "
            "can be covered at all is covered by at least two of these materials, so there is no "
            "organism that only one of them could contribute.</p>"
        )
    out.append(
        "<p class='note'>Front 1 means no other candidate is at least as good on every objective at "
        "once (coverage, assessed evidence, uncertainty, and how few materials are involved). "
        "A single material on front 1 is there because using one person is itself an advantage.</p>"
    )

    if ach["unreachable_goals"]:
        rows = "; ".join(
            f"{u['label']} (worth {_num(u['weight_points'], 1)} index points)"
            for u in ach["unreachable_goals"]
        )
        out.append(
            f"<div class='card warn'><strong>Gaps no candidate can fill:</strong> {_e(rows)}. "
            f"That is why the highest reachable index is {_num(ceiling)} and not 100.</div>"
        )

    # ---- 5. goal-group breakdown ----------------------------------------- #
    groups = sorted(ts["counting_groups"])
    out.append("<h3>Where each candidate's index comes from</h3>")
    out.append(
        "<p class='note'>One goal group, one vote: a guild of organisms and the gene panel measuring "
        "the same capacity share a group, and their combined credit is capped at that group's weight, "
        "so redundant organisms cannot multiply a capacity's importance. Equal group weights are the "
        "default; the group-size weighting in the table above is the alternative.</p>"
    )
    out.append("<table><thead><tr><th>Candidate</th>")
    for g in groups:
        out.append(f"<th class='num'>{_e(g)}</th>")
    out.append("<th class='num'>Total</th></tr></thead><tbody>")
    for s_ in [x for x in r["set_results"] if x["material_count"] <= 2][:10]:
        out.append(f"<tr><td>{_e(' + '.join(s_['member_person_ids']))}</td>")
        for g in groups:
            out.append(f"<td class='num'>{_num(s_['coverage']['per_group'].get(g), 1)}</td>")
        out.append(f"<td class='num'><strong>{_num(s_['coverage']['coverage'])}</strong></td></tr>")
    out.append("</tbody></table>")

    # ---- 6. target matrix ------------------------------------------------ #
    out.append("<h2>Every gap, and who has it</h2>")
    out.append(
        "<p class='note'><code>1.00</code> the material carries it · <code>0.00</code> the same "
        "profiler looked and did not find it, at its detection limit · <code>n/a</code> not assessed "
        "in that material · <code>context</code> shown for orientation and outside the score. "
        "Values between 0 and 1 are a graded supply index against the reference lower quartile.</p>"
    )
    donor_ids = [row["candidate_id"] for row in r["individual_results"]]
    avail: dict[tuple[str, str], dict[str, Any]] = {}
    for row in r["individual_results"]:
        for cell in row["per_target_availability"]:
            avail[(cell["target_id"], row["candidate_id"])] = cell
    out.append("<table><thead><tr><th>Gap</th><th>Where it came from</th><th>Recipient</th>")
    for did in donor_ids:
        out.append(f"<th class='num'>{_e(donor_names[did])}</th>")
    out.append("</tr></thead><tbody>")
    order = {"supply": 0, "support_range": 1, "avoid_introduction": 2, "monitor": 3}
    for t in sorted(r["targets"], key=lambda x: (order.get(x["direction"], 9), x["label"])):
        rec_txt = (
            "not detected"
            if t["recipient"]["status"] == "not_detected"
            else (
                f"{_num(t['recipient']['percentile'], 0)}th pct"
                if t["recipient"]["percentile"] is not None
                else "—"
            )
        )
        out.append(
            f"<tr><td>{_e(t['label'])}<div class='note'>{_e(t['desired_state'].get('statement'))}</div></td>"
            f"<td><span class='tag context'>{_e(t['origin'].replace('_', ' '))}</span>"
            f"<div class='note'>{_e(t['direction'].replace('_', ' '))}</div></td>"
            f"<td class='num'>{_e(rec_txt)}</td>"
        )
        for did in donor_ids:
            cell = avail.get((t["target_id"], did))
            if cell is None:
                out.append("<td class='num note'>context</td>")
            elif cell["supported_value"] is None:
                out.append("<td class='num'><span class='tag'>n/a</span></td>")
            else:
                out.append(f"<td class='num'>{_num(cell['supported_value'], 2)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")

    # ---- 7. concerns ----------------------------------------------------- #
    out.append("<h2>Open pathogen questions and screening gaps</h2>")
    out.append(
        "<p class='note warn'>Nothing below is a diagnosis and nothing below has been confirmed by a "
        "laboratory test. \"No confirmed exclusion identified\" is not the same statement as a negative "
        "pathogen screen: this run had no laboratory controls, and stool DNA sequencing cannot cover "
        "much of what a donor programme must test for. The recipient's own baseline carries findings of "
        "the same kind.</p>"
    )
    by_person: dict[str, list[dict[str, Any]]] = {}
    for c in r["concern_records"]:
        by_person.setdefault(donor_names.get(c["candidate_id"], c["candidate_id"]), []).append(c)
    for name in sorted(by_person):
        rows_c = by_person[name]
        lead = [c for c in rows_c if c["tier"] in LEAD_TIERS]
        rest = [c for c in rows_c if c["tier"] not in LEAD_TIERS]
        exc = [c for c in rows_c if c["disposition"] == "exclusion_confirmed"]
        out.append(f"<h3>{_e(name)}</h3>")
        out.append(
            f"<p class='kv'>{len(exc)} confirmed exclusion(s) · {len(lead)} finding(s) a person should "
            f"look at · {len(rest)} rows with no organism-specific evidence or context only</p>"
        )
        out.append(
            "<table><thead><tr><th>Finding</th><th>What it means</th><th>Evidence</th>"
            "<th>Rule</th><th>Why it is not settled</th></tr></thead><tbody>"
        )
        for c in sorted(lead, key=lambda x: (LEAD_TIERS.index(x["tier"]), x["feature_label"])):
            cls = {
                "exclusion_confirmed": "excluded",
                "review_pending": "review",
                "context_only": "context",
            }[c["disposition"]]
            m = c["measurements"]
            ev = (
                f"{m.get('unique_supporting_fragments')} specific fragment(s), "
                f"{m.get('informative_regions_supported')} region(s)"
                if "unique_supporting_fragments" in m
                else (
                    f"{m.get('supporting_fragments')} fragment(s)"
                    if "supporting_fragments" in m
                    else "—"
                )
            )
            out.append(
                f"<tr><td>{_e(c['feature_label'])}</td>"
                f"<td><span class='tag {cls}'>{_e(TIER_LABELS.get(c['tier'], c['tier']))}</span></td>"
                f"<td class='note'>{_e(ev)}</td>"
                f"<td><code>{_e(c['rule_id'])}</code></td>"
                f"<td class='note'>{_e(c['reason'])}</td></tr>"
            )
        out.append("</tbody></table>")
        if rest:
            out.append(
                f"<p class='note'>{len(rest)} further rows for {_e(name)} carry no organism-specific "
                "fragment at all, or are background, decoy or disease-association context. They are "
                "kept in full in <code>match.json</code>; none is evidence for or against any "
                "organism.</p>"
            )
    if r["screening_records"]:
        gaps = r["screening_records"][0]["not_covered_by_this_assay"]
        out.append("<h3>What this assay cannot cover at all</h3><ul class='note'>")
        out.extend(f"<li>{_e(g)}</li>" for g in gaps)
        out.append("</ul>")

    # ---- 8. context deliberately outside the score ----------------------- #
    out.append("<h2>Facts deliberately kept out of the score</h2>")
    first = r["individual_results"][0]["context_not_in_the_score"] if r["individual_results"] else {}
    if first:
        out.append(f"<p class='note'>{_e(first['gmwi2']['excluded_because'])}.</p>")
        out.append(f"<p class='note'>{_e(first['diversity_excluded_because'])}.</p>")
        out.append(f"<p class='note'>{_e(first['patterns_excluded_because'])}.</p>")
    out.append(
        "<table><thead><tr><th>Donor</th><th class='num'>GMWI2</th><th>Band</th>"
        "<th class='num'>Shannon</th><th class='num'>Species</th>"
        "<th>Highest disease-pattern resemblance</th></tr></thead><tbody>"
    )
    for row in sorted(r["individual_results"], key=lambda x: -x["goals"]["covered"]):
        c = row["context_not_in_the_score"]
        pats = "; ".join(
            f"{p['label']} ({_num(p['percentile'], 0)}th)" for p in c["highest_disease_patterns"]
        )
        out.append(
            f"<tr><td>{_e(row['person_id'])}</td><td class='num'>{_num(c['gmwi2']['score'], 2)}</td>"
            f"<td>{_e(c['gmwi2']['band'])}</td><td class='num'>{_num(c['shannon_diversity'], 2)}</td>"
            f"<td class='num'>{_e(c['species_richness'])}</td><td class='note'>{_e(pats)}</td></tr>"
        )
    out.append("</tbody></table>")

    # ---- 9. ecology, sensitivity, capabilities --------------------------- #
    out.append("<h2>Ecology and strain resolution</h2>")
    for row in r["individual_results"]:
        d = row["descriptors"]
        out.append(
            f"<div class='card'><strong>{_e(row['person_id'])}</strong> · Bray–Curtis "
            f"{_num(d['bray_curtis_dissimilarity'], 3)} · Jaccard {_num(d['jaccard_overlap'], 3)} · "
            f"{d['shared_species']} species shared with the recipient · {d['donor_unique_species']} "
            f"the recipient does not have"
            f"<div class='note'>strain resolution: {_e(d['strain_resolution']['status'])} — "
            f"{_e(d['strain_resolution']['reason'])}. {_e(d['strain_resolution']['consequence'])}.</div>"
            "</div>"
        )
    out.append(f"<p class='note'>{_e(r['individual_results'][0]['descriptors']['interpretation'])}</p>"
               if r["individual_results"] else "")

    out.append("<h2>Does the answer depend on my choices?</h2>")
    out.append(f"<p class='note'>{_e(r['sensitivity_results']['meaning'])}</p>")
    out.append(
        "<table><thead><tr><th>If we...</th><th class='num'>Scored gaps</th>"
        "<th>Leading candidates</th></tr></thead><tbody>"
    )
    for s in r["sensitivity_results"]["scenarios"]:
        leaders = ", ".join(
            " + ".join(donor_names.get(m, m) for m in cid.split("+")) for cid in s.get("leaders") or []
        )
        out.append(
            f"<tr><td>{_e(s['note'])}<div class='note'><code>{_e(s['scenario_id'])}</code></div></td>"
            f"<td class='num'>{_e(s.get('scoring_target_count'))}</td><td>{_e(leaders) or '—'}</td></tr>"
        )
    out.append("</tbody></table>")

    out.append("<h2>What this tool could not do</h2>")
    out.append("<table><thead><tr><th>Capability</th><th>Status</th><th>What is missing</th></tr></thead><tbody>")
    for c in r["capability_manifest"]:
        miss = "; ".join(c.get("missing_artifacts") or []) or c.get("note") or "—"
        out.append(
            f"<tr><td><code>{_e(c['capability_id'])}</code></td>"
            f"<td>{_e(c['status'].replace('_', ' '))}</td><td class='note'>{_e(miss)}</td></tr>"
        )
    out.append("</tbody></table><ul class='note'>")
    out.extend(f"<li>{_e(x)}</li>" for x in r["limitations"])
    out.append("</ul>")
    out.append(
        f"<footer>run <code>{_e(r['run_id'])}</code> · created {_e(r['created_at'])} · "
        "experimental computational matching. Not a diagnostic test, not a screening result and not a "
        "clinical decision. Physical preparation, dosing and administration are outside this tool."
        "</footer></div></body></html>"
    )
    return "".join(out)


def write_outputs(
    result: dict[str, Any], out_dir: Path, formats: list[str]
) -> dict[str, str]:
    """Write match.json / match.html / match.pdf as requested."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}
    if "json" in formats:
        p = out_dir / "match.json"
        p.write_text(json.dumps(result, indent=2, sort_keys=False, default=str))
        written["json"] = str(p)
    if "html" in formats:
        p = out_dir / "match.html"
        p.write_text(render_html(result))
        written["html"] = str(p)
    if "pdf" in formats:
        from .pdf import render_pdf

        p = out_dir / "match.pdf"
        render_pdf(result, p)
        written["pdf"] = str(p)
    manifest = out_dir / "provenance.json"
    manifest.write_text(
        json.dumps(
            {
                "run_id": result["run_id"],
                "created_at": result["created_at"],
                "software": result["software"],
                "request": result["request"],
                "input_audit": result["input_audit"],
                "search_audit": result["search_audit"],
                "capability_manifest": result["capability_manifest"],
                "outputs": written,
            },
            indent=2,
            default=str,
        )
    )
    written["provenance"] = str(manifest)
    return written
