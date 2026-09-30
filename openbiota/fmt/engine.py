"""The matching engine: one recipient, any number of donor materials.

Follows the deterministic core of spec §9.4. Everything the report and the
HTML/PDF renderers show is read from the single canonical result document this
module returns, so no format can compute a different number (V7-097).
"""

from __future__ import annotations

import datetime as dt
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

from openbiota import __version__
from openbiota.resolution import mechanisms as mech_mod

from . import capabilities, compatibility, strain_matching, suitability, transfer_risk
from . import concerns as concerns_mod
from . import coverage as cov_mod
from . import sets as sets_mod
from .identity import IdentityConflict, content_id, dedupe_materials
from .inputs import FmtInputError, MaterialData, lane_compatibility, load_material
from .targets import build_targets, goal_weights

SCHEMA = "fmt-match-result/6.0.0"


def _request_document(req: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "fmt-match-request/6.0.0",
        "request_id": req.get("request_id") or "fmt-match",
        "recipient": req["recipient"],
        "donors": req["donors"],
        "indication": req.get("indication"),
        "goal_policy": req.get("goal_policy") or "recipient-targets-v1",
        "custom_goals": req.get("custom_goals") or [],
        "screening_manifest": req.get("screening_manifest"),
        "model_policy": req.get("model_policy") or "standard",
        "objective_views": req.get("objective_views")
        or [{"type": "measured_coverage", "model_id": None, "endpoint": None, "horizon": None,
             "scenario_set_id": None}],
        "search": req.get("search") or {"max_subsets": 100000, "max_donors_per_set": None,
                                        "max_materials_per_set": None, "seed": 1701},
        "output_formats": req.get("output_formats") or ["json", "html", "pdf"],
    }


def match(request: dict[str, Any]) -> dict[str, Any]:
    """Run a complete comparison and return the canonical result document."""
    req = _request_document(request)
    started = dt.datetime.now().astimezone()

    # ---- ingest ---------------------------------------------------------- #
    conflicts: list[IdentityConflict] = []
    recipient_spec = req["recipient"]
    recipient = load_material(
        person_id=str(recipient_spec["person_id"]),
        role="recipient",
        root=Path(recipient_spec["input_root"]),
    )
    declared = recipient_spec.get("sample_id")
    if declared and declared != recipient.material.sample_id:
        raise FmtInputError(
            f"identity conflict: {recipient_spec['input_root']} holds sample "
            f"{recipient.material.sample_id!r}, not the declared {declared!r}"
        )

    donors: list[MaterialData] = []
    for spec in req["donors"]:
        try:
            data = load_material(
                person_id=str(spec["person_id"]), role="donor", root=Path(spec["input_root"])
            )
        except FmtInputError as exc:
            conflicts.append(
                IdentityConflict(
                    material_id=None,
                    sample_id=spec.get("sample_id"),
                    person_id=str(spec.get("person_id")),
                    kind="no_usable_analytical_input",
                    detail=str(exc),
                )
            )
            continue
        declared = spec.get("sample_id")
        if declared and declared != data.material.sample_id:
            conflicts.append(
                IdentityConflict(
                    material_id=data.material.material_id,
                    sample_id=data.material.sample_id,
                    person_id=data.material.person_id,
                    kind="declared_sample_id_mismatch",
                    detail=(
                        f"input root holds sample {data.material.sample_id!r} but the request declared "
                        f"{declared!r}; rule X03 excludes this analytical record, not the person"
                    ),
                )
            )
            continue
        if data.material.person_id == recipient.material.person_id:
            conflicts.append(
                IdentityConflict(
                    material_id=data.material.material_id,
                    sample_id=data.material.sample_id,
                    person_id=data.material.person_id,
                    kind="recipient_listed_as_donor",
                    detail="the recipient's own material cannot be a donor candidate",
                )
            )
            continue
        donors.append(data)

    materials, dedupe_notes = dedupe_materials([d.material for d in donors])
    kept_ids = {m.material_id for m in materials}
    unique: list[MaterialData] = []
    for d in donors:
        if d.material.material_id in kept_ids and all(
            d.material.material_id != u.material.material_id for u in unique
        ):
            unique.append(d)
    donors = unique

    # ---- goals: built once, from the recipient only ---------------------- #
    indication = (req.get("indication") or {}).get("code") if isinstance(req.get("indication"), dict) else None
    targets = build_targets(recipient, indication=indication, custom_goals=req["custom_goals"])
    alpha, beta = goal_weights(targets)

    # ---- per-material evidence ------------------------------------------- #
    screening_by_person = {
        str(k): v for k, v in ((req.get("screening_manifest") or {}).get("by_person") or {}).items()
    }
    lanes: dict[str, dict[str, Any]] = {}
    concern_records: list[concerns_mod.ConcernRecord] = []
    screening_records: list[dict[str, Any]] = []
    excluded_material_ids: set[str] = set()
    for d in donors:
        cid = d.material.material_id
        lanes[cid] = lane_compatibility(recipient, d)
        screening = screening_by_person.get(d.material.person_id)
        recs, screen = concerns_mod.evaluate(d, candidate_id=cid, screening=screening)
        hard = concerns_mod.exclusions_from_screening(
            candidate_id=cid, material_id=cid, screening=screening
        )
        concern_records.extend(recs + hard)
        screening_records.append(screen)
        if hard:
            excluded_material_ids.add(cid)

    scoring_targets = [t for t in targets if t.is_positive_scoring]
    table = cov_mod.build_table(targets, donors)
    person_of = {d.material.material_id: d.material.person_id for d in donors}
    label_of = {d.material.material_id: d.material.person_id for d in donors}
    sample_of = {d.material.material_id: d.material.sample_id for d in donors}

    # ---- individual results ---------------------------------------------- #
    individual: list[dict[str, Any]] = []
    excluded_view: list[dict[str, Any]] = []
    searchable_ids: list[str] = []
    for d in donors:
        mid = d.material.material_id
        cov = cov_mod.coverage(targets=targets, alpha=alpha, beta=beta, table=table, members=[mid])
        my_concerns = [c for c in concern_records if c.candidate_id == mid]
        is_excluded = mid in excluded_material_ids
        record = {
            "candidate_id": mid,
            "kind": "singleton",
            "person_id": d.material.person_id,
            "sample_id": d.material.sample_id,
            "member_material_ids": [mid],
            "material_count": 1,
            "unique_donor_count": 1,
            "analysis_status": "complete" if lanes[mid]["compatible"] else "partial",
            "exclusion_status": (
                "excluded_on_confirmed_evidence" if is_excluded else "no_confirmed_exclusion_identified"
            ),
            "screening_status": next(
                (s["screening_status"] for s in screening_records if s["candidate_id"] == mid), "not_provided"
            ),
            "review_status": (
                "high_consequence_review_pending"
                if any(c.severity == "high_consequence" and c.disposition == "review_pending" for c in my_concerns)
                else ("other_review_pending" if any(c.disposition == "review_pending" for c in my_concerns) else "none_identified")
            ),
            "clinical_release_status": "not_assessed_by_this_tool",
            "matching_status": "evidence_excluded" if is_excluded else "conditional_research_comparison",
            "lane_compatibility": lanes[mid],
            "coverage": cov.to_json(),
            "goals": _goal_counts(scoring_targets, table, [mid]),
            "context_not_in_the_score": _donor_context(d),
            "coverage_alternatives": {
                "equal_goal_groups": cov.coverage,
                "goal_groups_weighted_by_size": None,  # filled in after the sensitivity run
                "achievable_ceiling_equal_groups": None,
            },
            "per_target_availability": [
                table[(t.target_id, mid)].to_json() for t in targets if t.is_positive_scoring
            ],
            "descriptors": compatibility.descriptors(recipient, d),
            "concern_ids": [c.concern_id for c in my_concerns],
            "concern_counts": {
                "confirmed_exclusions": sum(1 for c in my_concerns if c.disposition == "exclusion_confirmed"),
                "review_pending": sum(1 for c in my_concerns if c.disposition == "review_pending"),
                "high_consequence_review_pending": sum(
                    1 for c in my_concerns
                    if c.disposition == "review_pending" and c.severity == "high_consequence"
                ),
                "context_only": sum(1 for c in my_concerns if c.disposition == "context_only"),
            },
            "model_prediction_ids": [],
            "unavailable_capability_reasons": [
                c["capability_id"] for c in capabilities.manifest() if c["status"] != "implemented"
            ],
        }
        if is_excluded:
            excluded_view.append(record)
        else:
            individual.append(record)
            searchable_ids.append(mid)

    # ---- set search ------------------------------------------------------- #
    search_cfg = req["search"]
    subsets, audit = sets_mod.feasible_subsets(
        searchable_ids,
        person_of=person_of,
        max_donors=search_cfg.get("max_donors_per_set"),
        max_materials=search_cfg.get("max_materials_per_set"),
        max_subsets=int(search_cfg.get("max_subsets") or 100000),
    )
    candidates = sets_mod.evaluate_candidates(
        targets=targets, alpha=alpha, beta=beta, table=table, subsets=subsets,
        person_of=person_of, label_of=label_of,
    )
    sets_mod.pareto_fronts(candidates)
    ordered = sets_mod.display_order(candidates)

    sens = sets_mod.sensitivity(
        targets=targets, table=table, subsets=subsets, person_of=person_of, label_of=label_of
    )
    freq = sens["leading_front_frequency"]
    n_scen = max(1, sens["scenario_count"])
    # The alternative weighting is computed anyway for sensitivity; surfacing it
    # beside the default explains why a candidate covering the same number of
    # gaps can score differently, instead of leaving that looking arbitrary.
    by_size = next(
        (s2.get("coverage_by_candidate") or {} for s2 in sens["scenarios"]
         if s2["scenario_id"] == "group_weights_by_size"),
        {},
    )

    # What is the best any combination could reach? Without this the index
    # lies about its own scale: if nothing supplies a goal group, its weight is
    # unreachable and a "75 of 100" is really "75 of 75".
    ceiling = cov_mod.coverage(
        targets=targets, alpha=alpha, beta=beta, table=table, members=searchable_ids
    )
    reachable_ids = {
        t.target_id for t in scoring_targets
        if (ceiling.per_target.get(t.target_id) or 0.0) > 0.0
    }
    unreachable = [
        {
            "target_id": t.target_id,
            "label": t.label,
            "counting_group": t.counting_group,
            "weight_points": round(100.0 * alpha.get(t.counting_group, 0.0) * beta.get(t.target_id, 0.0), 2),
            "why": "no candidate material has evidence of supplying it",
        }
        for t in scoring_targets
        if (ceiling.per_target.get(t.target_id) or 0.0) <= 0.0
    ]

    # ---- patterns a donor would introduce that the recipient does not have -- #
    # A resemblance percentile is still not a diagnosis. But a candidate at the
    # 90th percentile of a pattern the recipient is at the 15th percentile of,
    # for a condition whose patient stool has reproduced disease features in
    # recipient animals, is a different proposition from one that merely shares
    # a pattern the recipient already carries. That distinction now reaches the
    # recommendation instead of sitting in a context table nobody reads.
    # Measured candidate mechanisms, per donor. This is what replaced the
    # retired percentile bridge: evidence about the donor's own sequence,
    # reported per person and never pooled (spec §8.4, §10).
    mechanisms: dict[str, list[dict[str, Any]]] = {}
    for donor in donors:
        person = person_of.get(donor.material.material_id, donor.material.person_id)
        found = mech_mod.findings_for(person, donor.results)
        if found:
            mechanisms.setdefault(person, []).extend(f.to_json() for f in found)

    admissible = set(searchable_ids)
    # The profile registry is what makes the organism-level test possible:
    # it carries each profile's per-species direction, so "would this donor
    # bring a disease-enriched organism the recipient lacks" can be asked
    # instead of comparing two composite percentiles.
    try:
        from openbiota.profiles import load_profile_set

        profile_set = load_profile_set(Path("profiles"))
    except Exception:  # noqa: BLE001 - fall back to the percentile rule
        profile_set = None
    exposures = transfer_risk.exposures_for(
        recipient=recipient,
        donors=[d for d in donors if d.material.material_id in admissible],
        person_of=person_of,
        profile_set=profile_set,
    )

    donor_by_id = {d.material.material_id: d for d in donors}
    set_results: list[dict[str, Any]] = []
    for c in ordered:
        members = c.member_material_ids
        inherited = [x for x in concern_records if x.candidate_id in members]
        counts = _goal_counts(scoring_targets, table, members)
        score = match_score(scoring_targets, table, members, reachable_ids)
        mine_exposed = transfer_risk.candidate_exposures(exposures, c.member_person_ids)
        actionable_exposed = [e for e in mine_exposed if e.actionable]
        # Kept per person: a clean panel from one member does not offset a
        # candidate signal from another (acceptance 61).
        mine_mechanisms = {
            person: list(rows)
            for person, rows in mechanisms.items()
            if person in set(c.member_person_ids)
        }
        body = {
            "candidate_id": c.candidate_id,
            "kind": "singleton" if not c.is_set else "set",
            "member_material_ids": members,
            "member_person_ids": c.member_person_ids,
            "member_sample_ids": [sample_of[m] for m in members],
            "material_count": c.material_count,
            "unique_donor_count": c.unique_donor_count,
            "exposure_type": "unknown",
            "exposure_note": (
                "a computational hypothesis. Whether these materials would ever be combined, and how, "
                "is outside this tool: physical pooling, non-pooled multidonor and sequential exposure "
                "are different interventions and only the first is permutation-invariant"
            ),
            "candidate_mechanisms_by_person": mine_mechanisms,
            "candidate_mechanism_count": sum(len(v) for v in mine_mechanisms.values()),
            "introduces_new_patterns": [e.to_json() for e in mine_exposed],
            "introduces_with_transfer_evidence": [
                e.profile for e in mine_exposed
                if e.tier == transfer_risk.HUMAN_DONOR_TRANSFER
            ],
            # Patterns this candidate carries that the recipient does not, at
            # the exposure threshold.
            "observed_pattern_count": len(mine_exposed),
            # The species this candidate would actually add: disease-raised
            # organisms the recipient does not already carry. This is the
            # number the report headlines, because it is what the veto turns
            # on and what a reader can do something about.
            "introduced_organisms": sorted({
                organism
                for e in mine_exposed
                if e.tier == transfer_risk.HUMAN_DONOR_TRANSFER
                for organism in e.introduced_organisms
            }),
            "introduced_organism_count": len({
                organism
                for e in mine_exposed
                if e.tier == transfer_risk.HUMAN_DONOR_TRANSFER
                for organism in e.introduced_organisms
            }),
            # Of those, the ones whose condition has the strongest transfer
            # evidence tier. Reported so the recommendation can name them;
            # it does not exclude anybody.
            "transfer_evidenced_pattern_count": sum(
                1 for e in mine_exposed
                if e.tier == transfer_risk.HUMAN_DONOR_TRANSFER
            ),
            # Always zero while the percentile bridge stays retired. Kept so
            # the retirement is visible in the data rather than implied.
            "new_pattern_count": len(actionable_exposed),
            "match_score": score,
            "match_score_scale": (
                "0-100: the share of the recipient's reachable gaps this candidate supplies, "
                "each gap counted once, partial credit where supply is graded. Not a probability "
                "of engraftment and not a probability of benefit."
            ),
            "coverage": c.coverage.to_json(),
            "coverage_alternatives": {
                "equal_goal_groups": c.coverage.coverage,
                "goal_groups_weighted_by_size": by_size.get(c.candidate_id),
                "achievable_ceiling_equal_groups": ceiling.coverage,
                "why_they_differ": (
                    "the default gives every goal group an equal share, so a group holding one organism "
                    "counts as much as a group holding twelve. The alternative weights each group by how "
                    "many goals it contains. Both cap a group's credit, so neither lets one signal vote "
                    "twice; the ranking is reported under both"
                ),
            },
            "goals": counts,
            "plain_summary": _plain_summary(
                labels=[label_of[m] for m in members], counts=counts, marginals=c.marginals
            ),
            "marginals": c.marginals,
            "unique_target_ids_by_member": c.unique_targets,
            "redundant_target_ids": c.redundant_target_ids,
            "pareto_front": c.pareto_front,
            "display_order": c.display_order,
            "concern_union_ids": sorted({x.concern_id for x in inherited}),
            "concern_counts": {
                "review_pending": sum(1 for x in inherited if x.disposition == "review_pending"),
                "high_consequence_review_pending": sum(
                    1 for x in inherited
                    if x.disposition == "review_pending" and x.severity == "high_consequence"
                ),
                "context_only": sum(1 for x in inherited if x.disposition == "context_only"),
            },
            "concern_note": (
                "a set inherits every member's unresolved concern. Adding a donor cannot dilute or "
                "average one away"
            ),
            "stability": {
                "leading_front_scenarios": freq.get(c.candidate_id, 0),
                "scenario_count": n_scen,
                "fraction": round(freq.get(c.candidate_id, 0) / n_scen, 4),
                "meaning": "share of declared sensitivity scenarios on the leading front, not a success probability",
            },
            "predictions": {
                "status": "unavailable",
                "reason": "no engraftment, convergence or clinical-outcome model artifact is available",
            },
        }
        if c.is_set:
            body["set_descriptors"] = compatibility.set_descriptors(
                recipient, [donor_by_id[m] for m in members]
            )
            body["virtual_mixture"] = compatibility.virtual_mixture(
                [donor_by_id[m] for m in members], request.get("mixture_weights")
            )
        set_results.append(body)

    # ---- leaderboard ------------------------------------------------------ #
    # One ranked list over every evaluated candidate, singletons and
    # combinations together, so there is exactly one answer to "which is best"
    # and the order is reproducible.
    # Two counts, because they mean different things to a reader. "Serious"
    # is an organism that would matter if the evidence held up. "Other" is the
    # routine residue every stool sample produces: species a relative makes
    # unnameable, carriage without the toxin gene, resistance genes with no
    # resolved carrier.
    SERIOUS_TIERS = {
        "confirmed_exclusion", "high_consequence_supported", "high_consequence_unresolved",
    }
    OTHER_TIERS = {
        "species_unresolved", "carriage_toxin_negative", "toxin_or_resistance_review",
    }
    serious_by_candidate: dict[str, int] = {}
    other_by_candidate: dict[str, int] = {}
    for body in set_results:
        mine = [c for c in concern_records if c.candidate_id in body["member_material_ids"]]
        serious_by_candidate[body["candidate_id"]] = sum(
            1 for c in mine if c.tier in SERIOUS_TIERS
        )
        other_by_candidate[body["candidate_id"]] = sum(1 for c in mine if c.tier in OTHER_TIERS)

    # Donor suitability: the donor's own gut health, and what their material
    # would bring that the recipient does not have. Ranking used to ignore
    # both, which is how a donor with a negative gut-health index and a
    # 90th-percentile colorectal-cancer pattern reached the top three on
    # gap coverage alone.
    person_health: dict[str, suitability.DonorHealth] = {}
    for donor in donors:
        if donor.material.material_id not in admissible:
            continue
        person = person_of.get(donor.material.material_id, donor.material.person_id)
        person_health[person] = suitability.assess_donor(
            donor=donor,
            recipient=recipient,
            person_id=person,
            evidence=transfer_risk.EVIDENCE,
            actionable_tiers=transfer_risk.ACTIONABLE_TIERS,
            adverse_organisms=transfer_risk.adverse_organisms,
            profile_set=profile_set,
            high_percentile=transfer_risk.HIGH_PERCENTILE,
            new_margin=transfer_risk.NEW_MARGIN,
        )
    # Species abundances, carried into the result so the per-candidate pages
    # can sort what a material would add into "restores a gap", "probably
    # good", "unremarkable" and "could be introducing" without re-reading
    # every donor's results.json.
    species_by_person: dict[str, dict[str, float]] = {
        person_of.get(d.material.material_id, d.material.person_id):
            suitability._species_abundance(d)
        for d in donors
        if d.material.material_id in admissible
    }
    recipient_species = suitability._species_abundance(recipient)

    for body in set_results:
        health = suitability.combine(person_health, body["member_person_ids"])
        body["donor_health"] = health.to_json()
        body["suitability_score"] = round(
            suitability.suitability_score(
                restoration=body["match_score"] or 0.0, health=health
            ),
            1,
        )
        body["suitability_explanation"] = suitability.score_explanation(
            restoration=body["match_score"] or 0.0, health=health
        )
        body["eligible"] = health.eligible
        body["gate_failures"] = list(health.gate_failures)

    def _pattern_penalty(body: dict[str, Any]) -> int:
        """Rank cost of the patterns this candidate would bring.

        Sorted on before the score, so a candidate carrying a
        transfer-evidenced pattern cannot outrank one that does not,
        whatever their scores. Eight candidates previously tied at 100.0
        here and the tie fell to `serious_open_questions` - a count of
        unconfirmed trace findings - which handed the top slot to a donor
        carrying colorectal cancer at the 90th percentile against a
        recipient at the 15th, three trace findings ahead of an identical
        candidate that carried none.
        """
        # The veto weight is charged only where the donor would actually
        # bring a disease-enriched organism the recipient lacks. An
        # exposure flagged by the percentile alone, with no such organism
        # behind it, is context and costs a single point - otherwise the
        # penalty and `Exposure.actionable` would disagree about the same
        # exposure.
        total = 0
        for e in body["introduces_new_patterns"]:
            tier = e["transfer_evidence_tier"]
            brings = bool(e.get("introduced_organisms"))
            if tier in transfer_risk.ACTIONABLE_TIERS and not brings:
                total += transfer_risk.TIER_RANK_PENALTY[transfer_risk.ASSOCIATION_ONLY]
            else:
                total += transfer_risk.TIER_RANK_PENALTY.get(tier, 1)
        if transfer_risk.ASSOCIATION_WEIGHTING == "weighted":
            # Keep only the veto. Association-only patterns fall through to
            # the tie-break below, so a materially better restoration score
            # is not surrendered to avoid a weakly-evidenced association.
            veto = transfer_risk.TIER_RANK_PENALTY[transfer_risk.HUMAN_DONOR_TRANSFER]
            return total - (total % veto)
        return total

    ranked = sorted(
        set_results,
        key=lambda b: (
            # 1. Eligibility. A material that fails a safety gate - a
            #    dysbiotic donor, a genuinely new transfer-evidenced
            #    pattern, an organism above the healthy range - ranks below
            #    every material that passes, whatever it would supply.
            not b["eligible"],
            # 2. Suitability: restoration and donor health together, with
            #    the arithmetic printed beside it in the report.
            -b["suitability_score"],
            # 3. Then the old tie-breaks, which now only separate materials
            #    already equal on safety and suitability.
            -(b["match_score"] or 0.0),
            -b["goals"]["fully_supplied"],
            b["unique_donor_count"],
            serious_by_candidate.get(b["candidate_id"], 0),
            b["candidate_id"],
        ),
    )
    # The former constraint -- no newly introduced high-resemblance disease
    # pattern -- is retired with the bridge that produced it (spec §8.4): a
    # resemblance percentile is not a transferable feature and may not exclude
    # a donor. Measured candidate mechanisms are reported per donor beside the
    # ranking for a person to weigh; they are not folded into the score, and
    # one donor's clean panel never offsets another's finding.
    # The recommendation may not carry a transfer-evidenced new pattern.
    # `new_pattern_count` counts exposures in ACTIONABLE_TIERS, which now
    # holds human_donor_transfer, so this filter has teeth again.
    #
    # `leader` stays the top of the ranking so the report can say plainly
    # when the two differ and why. They differ exactly when the best
    # restoration score belongs to a candidate carrying a veto pattern.
    clean = [b for b in ranked if b["eligible"]]
    recommended = clean[0] if clean else None
    leader = ranked[0] if ranked else None
    for position, body in enumerate(ranked, start=1):
        body["match_rank"] = position
        body["recommended"] = bool(
            recommended and body["candidate_id"] == recommended["candidate_id"]
        )
        body["serious_open_questions"] = serious_by_candidate.get(body["candidate_id"], 0)
        body["other_findings"] = other_by_candidate.get(body["candidate_id"], 0)
        body["rank_reason"] = _rank_reason(
            rank=position,
            counts=body["goals"],
            reachable=len(reachable_ids),
            people=body["unique_donor_count"],
            serious=body["serious_open_questions"],
            leader=leader,
        )
        short: list[str] = []
        if body["recommended"]:
            # Was "best score with no new disease pattern", which asserted
            # something the ranking never checked and which was false here:
            # the top candidate carried three.
            short.append("highest score on restoration")
        reachable_missing = len(reachable_ids) - body["goals"]["covered"]
        if reachable_missing > 0:
            short.append(f"misses {reachable_missing} reachable gap(s)")
        if leader and body["unique_donor_count"] > leader["unique_donor_count"] and \
                body["match_score"] == leader["match_score"]:
            short.append(f"same score using {body['unique_donor_count']} people")
        body["rank_reason_short"] = (
            "; ".join(short)
            or ("supplies every reachable gap" if not reachable_missing else "—")
        )
        # Patterns this candidate carries that the recipient does not, kept
        # in their own field rather than appended to `rank_reason`. They did
        # not affect the rank, so they do not belong in the sentence that
        # explains the rank, and a renderer that shows both would say it
        # twice. This used to be gated on `new_pattern_count` (always zero,
        # so it never printed) and to open "Not recommended:", which was
        # never true of anything.
        body["carried_patterns_note"] = (
            (
                "Carries, and you do not: "
                + "; ".join(
                    f"{e['short_label']} ({e['donor_percentile']:.0f}th vs your "
                    f"{e['recipient_percentile']:.0f}th, carried by "
                    f"{', '.join(e['donor_person_ids'])})"
                    for e in body["introduces_new_patterns"]
                )
                + ". This did not change its rank \u2014 a resemblance percentile is not "
                "a transferable feature \u2014 and is listed so it can be weighed."
            )
            if body["introduces_new_patterns"]
            else ""
        )
    tied = [
        b["candidate_id"] for b in ranked[1:]
        if leader and b["match_score"] == leader["match_score"]
    ]
    best_single = next((b for b in ranked if b["unique_donor_count"] == 1), None)

    for row in individual:
        row["match_score"] = next(
            (b["match_score"] for b in set_results if b["candidate_id"] == row["candidate_id"]), None
        )
        row["match_rank"] = next(
            (b["match_rank"] for b in set_results if b["candidate_id"] == row["candidate_id"]), None
        )
        row["coverage_alternatives"]["goal_groups_weighted_by_size"] = by_size.get(row["candidate_id"])
        row["coverage_alternatives"]["achievable_ceiling_equal_groups"] = ceiling.coverage
        row["plain_summary"] = _plain_summary(
            labels=[row["person_id"]],
            counts=row["goals"],
            marginals=next(
                (s2["marginals"] for s2 in set_results if s2["candidate_id"] == row["candidate_id"]), {}
            ),
        )

    # ---- canonical document ---------------------------------------------- #
    scoring = scoring_targets
    result: dict[str, Any] = {
        "schema_version": SCHEMA,
        "run_id": content_id(
            "fmtrun",
            {
                "request": req,
                "recipient": recipient.material.material_id,
                "donors": sorted(person_of),
                "targets": sorted(t.target_id for t in targets),
                "version": __version__,
            },
        ),
        "created_at": started.isoformat(timespec="seconds"),
        "software": {"name": "openbiota", "version": __version__, "component": "openbiota.fmt"},
        "request": req,
        "recipient": {
            "person_id": recipient.material.person_id,
            "sample_id": recipient.material.sample_id,
            "material_id": recipient.material.material_id,
            "material": recipient.material.to_json(),
            "lane": recipient.lane,
            "read_pairs": (recipient.results.get("input") or {}).get("read_pairs"),
            "indication": req.get("indication"),
        },
        "materials": [d.material.to_json() for d in donors],
        "input_audit": [
            {
                "person_id": d.material.person_id,
                "sample_id": d.material.sample_id,
                "input_root": d.material.input_root,
                "artifact_ids": d.material.input_artifact_ids,
                "evidence_level": "native_measurement",
                "chosen": "results.json",
                "rejected": [],
            }
            for d in [recipient, *donors]
        ] + [{"note": n} for n in dedupe_notes],
        "goal_scenarios": [
            {
                "scenario_id": "research_goals",
                "default": True,
                "inclusion": "every active or conditional target",
                "scoring_target_count": len(scoring),
                "goal_group_weights": alpha,
                "within_group_weights": beta,
                "weight_note": "equal group weights then equal weights within a group: engineering defaults, not efficacy coefficients",
            },
            {
                "scenario_id": "reference_only",
                "default": False,
                "inclusion": "active or conditional targets whose origin is reference_supported",
                "scoring_target_count": len(
                    [t for t in scoring if t.origin == "reference_supported"]
                ),
            },
        ],
        "targets": [t.to_json() for t in targets],
        "target_summary": {
            "total": len(targets),
            "positive_scoring": len(scoring),
            "by_direction": _count(targets, "direction"),
            "by_origin": _count(targets, "origin"),
            "by_status": _count(targets, "target_status"),
            "counting_groups": sorted({t.counting_group for t in scoring}),
            "note": (
                "context-only and avoid-introduction targets are displayed but never enter the benefit "
                "denominator, and no negative abundance credit exists anywhere in this tool"
            ),
        },
        "individual_results": individual,
        "set_results": set_results,
        "excluded_candidates": excluded_view,
        "not_computable_candidates": [c.to_json() for c in conflicts],
        "concern_records": [c.to_json() for c in concern_records],
        "screening_records": screening_records,
        "concern_rule_registry": concerns_mod.RULES,
        "evidence_sources": concerns_mod.SOURCES,
        "capability_manifest": capabilities.manifest(),
        "model_predictions": [],
        "sensitivity_results": sens,
        "search_audit": audit,
        # Strain matching, as three separate capabilities (spec §10.1). The
        # baseline comparison is computable from the specimens in hand; the
        # other two report what they are missing rather than returning zero.
        "strain_matching": strain_matching.baseline(
            recipient=recipient.material.person_id,
            donors=[d.material.person_id for d in donors],
            materials={
                m.material.person_id: Path(m.material.input_root)
                for m in [recipient, *donors]
            },
        ).to_json(),
        "match_score_semantics": strain_matching.MATCH_SCORE_SEMANTICS,
        # Per-person species abundances, for the per-candidate pages.
        "species_by_person": {
            person: {k: round(v, 5) for k, v in sorted(m.items()) if v > 0}
            for person, m in species_by_person.items()
        },
        "recipient_species": {
            k: round(v, 5) for k, v in sorted(recipient_species.items()) if v > 0
        },
        "leaderboard": {
            "ranked_candidate_ids": [b["candidate_id"] for b in ranked],
            "recommended_id": recommended["candidate_id"] if recommended else None,
            "recommended_members": recommended["member_person_ids"] if recommended else [],
            "recommended_score": recommended["match_score"] if recommended else None,
            "recommended_rank": recommended["match_rank"] if recommended else None,
            # What actually happens. Verified against behaviour by
            # test_no_output_claims_a_gate_that_cannot_fire, which exists
            # because this string once described a gate that had been
            # switched off.
            "recommended_rule": (
                "the highest-scoring candidate that would bring you no organism a "
                "transfer-evidenced disease is enriched for and that you do not "
                "already carry. Transfer-evidenced means patient stool reproduced "
                "disease features in recipient animals against healthy-donor "
                "controls. The test is at organism level, not percentile: a "
                "resemblance percentile mixes organisms raised by a disease with "
                "organisms lost to it, and restoring the lost ones is the point of "
                "the transfer, so only the raised half can count as introducing "
                "anything. Organisms you already carry never count. Candidates that "
                "would bring one are ranked below those that would not, whatever "
                "they score."
            ),
            "recommended_note": (
                (
                    f"{'+'.join(recommended['member_person_ids'])} is recommended over the "
                    f"highest-scoring candidate: it scores {recommended['match_score']:.1f} at "
                    f"rank {recommended['match_rank']} on restoration alone, and it is the best "
                    "candidate that introduces no new transfer-evidenced pattern."
                )
                if recommended and leader and recommended["candidate_id"] != leader["candidate_id"]
                else None
            ),
            "no_clean_candidate_note": (
                None if recommended else
                "Every candidate would introduce at least one new transfer-evidenced pattern. "
                "There is no clean choice here; the exposures are listed so the trade-off is "
                "explicit."
            ),
            "highest_scoring_id": leader["candidate_id"] if leader else None,
            "highest_scoring_members": leader["member_person_ids"] if leader else [],
            "new_pattern_exposures": [e.to_json() for e in exposures],
            "resemblance_is_context_only": transfer_risk.RETIREMENT_NOTE,
            "candidate_mechanisms_by_person": mechanisms,
            "mechanism_ledger_note": (
                "Candidate mechanisms measured in each donor's own sequence data. A "
                "gene-family fragment count is the weakest evidence there is: none of these "
                "is a confirmed pathotype, the carrier organism is unresolved, and each row "
                "states what a confirmed call would require. Reported per person and never "
                "pooled across a donor set."
            ),
            "exposure_rule": (
                f"a donor pattern at or above the {transfer_risk.HIGH_PERCENTILE:.0f}th percentile "
                "that exceeds the recipient's own reading by at least "
                f"{transfer_risk.NEW_MARGIN:.0f} percentile points"
            ),
            "best_match_id": leader["candidate_id"] if leader else None,
            "best_match_members": leader["member_person_ids"] if leader else [],
            "best_match_score": leader["match_score"] if leader else None,
            "best_match_rank_reason": leader["rank_reason"] if leader else None,
            "tied_with_best": tied,
            "tie_note": (
                (
                    f"{len(tied) + 1} candidates supply every reachable gap and therefore share the "
                    "top score. They are separated by the number of people involved, then by how "
                    "many findings need a laboratory test — never by trading one against the other."
                )
                if tied
                else "No other candidate matches the leader's score."
            ),
            "best_single_donor_id": best_single["candidate_id"] if best_single else None,
            "best_single_donor_members": best_single["member_person_ids"] if best_single else [],
            "best_single_donor_score": best_single["match_score"] if best_single else None,
            "best_single_donor_note": (
                (
                    f"{'+'.join(best_single['member_person_ids'])} is the strongest single material "
                    f"at {best_single['match_score']:.1f}, involving one person rather than "
                    f"{leader['unique_donor_count']}."
                    + (
                        " It is not recommended: it would introduce "
                        + ", ".join(
                            e["short_label"]
                            for e in best_single["introduces_new_patterns"]
                            if e["counts_against_selection"]
                        )
                        + ", which you do not have."
                        if best_single["new_pattern_count"]
                        else ""
                    )
                )
                if best_single and leader and best_single["candidate_id"] != leader["candidate_id"]
                else None
            ),
            "reachable_goals": len(reachable_ids),
            "scored_goals": len(scoring),
            "ranking_rule": (
                "match score, then how many gaps are fully rather than partly supplied, then the "
                "fewest people involved, then the fewest findings needing a lab test, then a stable "
                "candidate ID. Findings are a tie-break only: they are never traded against the "
                "score, and a candidate is never promoted for having fewer of them."
            ),
            "score_meaning": (
                "the share of the recipient's reachable gaps a candidate supplies. 100 means it "
                "supplies everything these materials could supply between them"
            ),
            "not": [
                "a probability that any organism will establish",
                "a probability of clinical benefit",
                "a safety ranking or a screening result",
            ],
        },
        "achievable": {
            "index_ceiling": ceiling.coverage,
            "index_scale": 100.0,
            "goals_reachable": len(scoring) - len(unreachable),
            "goals_total": len(scoring),
            "unreachable_goals": unreachable,
            # A goal "no candidate supplies" is a statement about the lane
            # the matching runs on, which is MetaPhlAn 3 because that is
            # the namespace of the 3,027-sample reference cohort. The
            # MetaPhlAn 4 lane is four years newer and built from far more
            # genomes, and where it disagrees the reader should see it
            # rather than read a database limit as a biological fact.
            "unreachable_contradicted_by_newer_lane": _newer_lane_check(
                unreachable, targets, donors, person_of
            ),
            "note": (
                "the highest index any combination of the supplied materials could reach. Points for a "
                "goal no material supplies are unreachable, so read coverage against this ceiling, not "
                "against 100"
            ),
        },
        "how_to_read": _how_to_read(ceiling.coverage, len(scoring), len(unreachable)),
        "clinical_release_status": "not_assessed_by_this_tool",
        "claim_audit": [
            {
                "quantity": "supported target coverage",
                "scale": "index 0-100",
                "claim_class": "measured_comparison",
                "means": "share of the recipient's frozen target weight that at least one member has evidence of supplying",
                "does_not_mean": [
                    "percent of microbiome restored",
                    "engraftment or persistence probability",
                    "compatibility probability",
                    "clinical match accuracy",
                ],
            },
            {
                "quantity": "scenario bounds",
                "claim_class": "missing_data_bounds",
                "means": "lower and upper analytical scenarios over unassessed and trace features",
                "does_not_mean": ["a confidence interval", "a statistical prediction interval"],
            },
            {
                "quantity": "leading-front frequency",
                "claim_class": "sensitivity_summary",
                "means": "how often a candidate leads across declared scenarios",
                "does_not_mean": ["probability of clinical benefit"],
            },
        ],
        "limitations": [
            "Sequence evidence is not an infection diagnosis and the absence of a confirmed exclusion "
            "is not a negative pathogen screen.",
            "No donor here is screened, cleared or released for any clinical use by this tool.",
            "Coverage compares supply evidence. Whether any organism would establish in this recipient "
            "is a different endpoint, and no model for it is available in this environment.",
            "A donor set is a computational hypothesis; the behaviour of physically combined material "
            "cannot be inferred from separate sequencing.",
            "No preparation, mixing ratio, dose, administration or conditioning instruction is produced.",
        ],
    }
    return result


def match_score(
    scoring: list[Any],
    table: dict[tuple[str, str], Any],
    members: list[str],
    reachable_ids: set[str],
) -> float | None:
    """The single number the leaderboard ranks on, 0-100.

    **The share of the recipient's reachable gaps this candidate supplies**,
    counting each gap once and giving partial credit where supply is graded.
    Two deliberate choices make it readable:

    * the denominator is what *any* candidate here could supply, so 100 means
      "supplies everything available from these materials" rather than an
      unreachable ideal;
    * every gap weighs the same, so the number tracks the plain count a reader
      can check. The group-capped index stays beside it as the guarantee that
      no single biological signal is counted twice.

    Hazards are deliberately absent. A pathogen finding is not tradeable
    against restoration, so it is never folded into this score; it sits in its
    own column and in its own section.
    """
    pool = [t for t in scoring if t.target_id in reachable_ids]
    if not pool:
        return None
    total = 0.0
    for t in pool:
        total += max(
            (table[(t.target_id, m)].contribution for m in members if (t.target_id, m) in table),
            default=0.0,
        )
    return round(100.0 * total / len(pool), 1)


#: Below this many MetaPhlAn 3 marker genes, the older lane is too thinly
#: equipped for a species for its silence to mean anything. Chosen against
#: the database's own distribution: the median species has 100 markers and
#: the 10th percentile has 5, so 20 sits comfortably in the thin tail.
_THIN_MARKER_COVERAGE: Final = 20


@lru_cache(maxsize=1)
def _mp3_marker_counts() -> dict[str, int]:
    """Marker genes per species in the MetaPhlAn 3 database.

    Precomputed to a small JSON because the database pickle is large and
    takes about ten seconds to load. Absent file means no check, which
    degrades to reporting nothing rather than to reporting noise.
    """
    path = Path("refs/metaphlan_db/mpa_v31_marker_counts.json")
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text()).get("markers_per_species", {})
    except (OSError, json.JSONDecodeError):
        return {}


def _newer_lane_check(
    unreachable: list[dict[str, Any]],
    targets: Any,
    donors: Any,
    person_of: Any,
) -> list[dict[str, Any]]:
    """Goals called unsupplied where the newer lane genuinely disagrees.

    Matching runs on MetaPhlAn 3 because the reference cohort is in that
    namespace and a percentile needs a same-namespace reference.

    The two lanes differ constantly, and almost all of it is uninteresting:
    MetaPhlAn 4's database is four years newer and built from roughly a
    million genomes against MetaPhlAn 3's seventeen thousand, so it finds
    organisms the older lane has no reference for at all. Across these
    samples that accounts for 466 of about 512 MetaPhlAn-4-only calls.
    Flagging those as "the databases disagree" would be crying wolf: the
    older lane was never able to see them.

    A disagreement worth printing is one where the older lane is **well
    equipped for that species and still did not call it**. Adlercreutzia
    equolifaciens is the case here: 150 markers in the MetaPhlAn 3
    database, the 87th percentile of species coverage, and the newer lane
    puts it at 0.20% in a donor the older lane calls clean. That is not a
    missing reference and it is not a depth limit.
    """
    marker_counts = _mp3_marker_counts()
    by_label = {t.label: t for t in targets}
    out: list[dict[str, Any]] = []
    for goal in unreachable:
        target = by_label.get(goal.get("label"))
        if target is None:
            continue
        features = set(getattr(target, "feature_ids", ()) or ())
        # Only species the older lane was actually equipped to find.
        covered = {
            f for f in features
            if marker_counts.get(f, 0) >= _THIN_MARKER_COVERAGE
        }
        if not covered:
            continue
        found: dict[str, float] = {}
        for donor in donors:
            person = person_of.get(
                donor.material.material_id, donor.material.person_id
            )
            catalogue = (
                donor.results.get("extended_catalogue") or {}
            ).get("sgbs") or []
            total = sum(
                float(row.get("percent") or 0.0)
                for row in catalogue
                if row.get("species") in covered
            )
            if total > 0:
                found[person] = round(total, 4)
        if not found:
            continue
        markers = max(marker_counts.get(f, 0) for f in covered)
        out.append({
            "label": goal.get("label"),
            "feature_ids": sorted(covered),
            "found_by_newer_lane_in": found,
            "matching_lane": "metaphlan3:mpa_v31_CHOCOPhlAn_201901",
            "newer_lane": "metaphlan4:mpa_vJun23_CHOCOPhlAnSGB_202403",
            "mp3_marker_genes": markers,
            "note": (
                "The lane this matching runs on does not detect this organism "
                "in any candidate, so the goal is scored as unreachable. The "
                f"newer MetaPhlAn 4 lane does. The older lane carries {markers} "
                "marker genes for this species, so this is not a case of it "
                "lacking a reference - it was equipped to find it and did not. "
                "Treat this gap as open rather than impossible."
            ),
        })
    return out


def _rank_reason(
    *,
    rank: int,
    counts: dict[str, Any],
    reachable: int,
    people: int,
    serious: int,
    leader: dict[str, Any] | None,
) -> str:
    """One sentence saying why this candidate sits where it does."""
    covered = counts["covered"]
    bits = [f"supplies {covered} of the {reachable} gaps any candidate here could supply"]
    if counts["partly_supplied"]:
        bits.append(f"{counts['partly_supplied']} of them only partly")
    bits.append(f"{people} {'person' if people == 1 else 'people'} involved")
    if serious:
        bits.append(f"{serious} finding{'s' if serious != 1 else ''} needing a lab test")
    # Rank 1 is the best *permitted* candidate, which is not always the
    # highest-scoring one: a candidate carrying a transfer-evidenced pattern
    # is ranked below every candidate without one whatever it scores.
    head = "Ranked first: " if rank == 1 else "This candidate "
    tail = ""
    if rank > 1 and leader is not None:
        lead_cov = leader["goals"]["covered"]
        lead_people = leader["unique_donor_count"]
        if covered == lead_cov and people > lead_people:
            tail = (
                f" Ranked below the leader because it reaches the same {covered} gaps "
                f"using {people} people instead of {lead_people}."
            )
        elif covered < lead_cov:
            missing = [
                lab for lab in counts["not_supplied_labels"]
                if lab in leader["goals"]["fully_supplied_labels"]
            ]
            tail = (
                f" Ranked below the leader because it misses {lead_cov - covered} gap(s) the "
                f"leader supplies"
                + (f": {'; '.join(missing[:2])}" if missing else "")
                + "."
            )
    return head + "; ".join(bits) + "." + tail


def _goal_counts(
    scoring: list[Any],
    table: dict[tuple[str, str], Any],
    members: list[str],
) -> dict[str, Any]:
    """How many of the recipient's gaps this candidate covers, by count.

    The weighted index answers "how much of the goal *weight*"; a plain count
    answers "how many of the things I am missing", which is the question a
    reader actually asks first.
    """
    full: list[str] = []
    partial: list[str] = []
    none: list[str] = []
    # Goals credited only because the newer profiler saw the organism. The
    # reader should be able to tell those apart from goals both profilers
    # agreed on, so the provenance travels with the label.
    secondary: dict[str, float] = {}
    for t in scoring:
        cells = [table[(t.target_id, m)] for m in members if (t.target_id, m) in table]
        best = max((c.contribution for c in cells), default=0.0)
        winner = max(cells, key=lambda c: c.contribution, default=None)
        if (
            winner is not None
            and getattr(winner, "rule_id", "") == "A-PRESENCE-SECONDARY-LANE"
            and winner.contribution >= 0.999
        ):
            secondary[t.label] = float(getattr(winner, "donor_value", 0.0) or 0.0)
        if best >= 0.999:
            full.append(t.label)
        elif best > 0:
            partial.append(t.label)
        else:
            none.append(t.label)
    return {
        "total": len(scoring),
        "covered": len(full) + len(partial),
        "fully_supplied": len(full),
        "partly_supplied": len(partial),
        "not_supplied": len(none),
        "fully_supplied_labels": full,
        "partly_supplied_labels": partial,
        "not_supplied_labels": none,
        "supplied_by_secondary_lane": secondary,
        "meaning": (
            "counted goals, not weighted points: a goal is covered when at least one member has "
            "evidence of supplying it at all"
        ),
    }


def _plain_summary(
    *, labels: list[str], counts: dict[str, Any], marginals: dict[str, Any]
) -> str:
    """One sentence a reader can act on, in words rather than indices."""
    who = " + ".join(labels)
    missing = counts["not_supplied_labels"]
    gain = marginals.get("gain_over_global_best_singleton")
    head = (
        f"{who} covers {counts['covered']} of the {counts['total']} gaps in the recipient's ledger"
    )
    if counts["partly_supplied"]:
        head += f" ({counts['partly_supplied']} of them only partly)"
    tail = ""
    if missing:
        shown = "; ".join(missing[:3])
        more = f" and {len(missing) - 3} more" if len(missing) > 3 else ""
        tail = f". Not covered: {shown}{more}"
    if len(labels) > 1 and isinstance(gain, (int, float)):
        if abs(gain) < 1e-9:
            tail += ". Adding these materials together covers nothing the best single material does not"
        elif gain > 0:
            tail += f". Together they cover {gain:.1f} index points more than the best single material"
    return head + tail + "."


def _donor_context(data: Any) -> dict[str, Any]:
    """Facts a reader will look for that are deliberately not in the score."""
    def val(kind: str, key: str) -> Any:
        obs = data.observation(kind, key)
        return obs.value if obs else None

    gmwi2 = data.observation("ecology", "gmwi2")
    patterns = sorted(
        (o for o in data.observations if o.feature_kind == "profile" and o.percentile is not None),
        key=lambda o: -(o.percentile or 0),
    )[:3]
    return {
        "gmwi2": {
            "score": gmwi2.value if gmwi2 else None,
            "band": (gmwi2.extra.get("band") if gmwi2 else None),
            "excluded_because": (
                "GMWI2 asks whether this donor's own community looks healthy. The match asks whether "
                "this donor has what the recipient lacks. They are different questions, and GMWI2 is "
                "computed from the same species table the coverage already uses, so scoring it too "
                "would count one signal twice (rule C01)"
            ),
        },
        "shannon_diversity": val("ecology", "shannon_index"),
        "species_richness": val("ecology", "species_richness"),
        "diversity_excluded_because": (
            "a high-diversity donor can carry unwanted organisms and a low-diversity donor can carry "
            "the one missing strain, so diversity has no universal direction"
        ),
        "highest_disease_patterns": [
            {"label": o.extra.get("label"), "percentile": o.percentile} for o in patterns
        ],
        "patterns_excluded_because": (
            "pattern resemblance is not a diagnosis and not a transmissible entity; it is shown as "
            "context and can never exclude a material"
        ),
    }


def _how_to_read(ceiling: float | None, n_goals: int, n_unreachable: int) -> list[dict[str, str]]:
    """The definitions a first-time reader needs, in the report itself."""
    top = f"{ceiling:.0f}" if isinstance(ceiling, (int, float)) else "—"
    return [
        {
            "term": "Gaps covered",
            "plain": (
                f"The recipient's ledger has {n_goals} scored gaps: organisms or gene capacities the "
                "recipient is missing or low in, each one derived from the recipient's own data before "
                "any donor was looked at. 'Covered' means at least one candidate material has evidence "
                "of carrying it. Read this number first."
            ),
        },
        {
            "term": "Supported target coverage",
            "plain": (
                f"The same information as a weighted index from 0 to 100, where gaps are grouped so one "
                f"biological signal cannot be counted twice. The most any combination here could reach "
                f"is {top}, because {n_unreachable} gap(s) are supplied by no candidate at all. Compare "
                f"candidates against {top}, not against 100."
            ),
        },
        {
            "term": "Scenario range",
            "plain": (
                "How the index moves if borderline calls are read the other way (a trace-level "
                "detection dropped, an unmeasured feature counted). It is a range of readings, not a "
                "statistical confidence interval."
            ),
        },
        {
            "term": "Open pathogen question",
            "plain": (
                "Sequencing found something that could be a pathogen, toxin gene or resistance gene, "
                "and the data cannot confirm or dismiss it. None of these is a diagnosis, and none of "
                "them has been confirmed by a laboratory test. They are listed so a person can decide "
                "what to test."
            ),
        },
        {
            "term": "Serious organism, identity not confirmed",
            "plain": (
                "An open question about an organism that would matter if it were really there. Most are "
                "a handful of DNA fragments, or many fragments piled on a single genome region, which is "
                "what a shared or conserved sequence looks like. A laboratory test settles it; this "
                "report cannot."
            ),
        },
        {
            "term": "What is deliberately not in the score",
            "plain": (
                "GMWI2, diversity, microbiome age and disease-pattern percentiles. Each is shown as "
                "context beside every candidate with the reason it is excluded. A donor with the best "
                "GMWI2 is not necessarily the donor who carries what this recipient lacks."
            ),
        },
    ]


def _count(items: list[Any], attr: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for item in items:
        key = str(getattr(item, attr))
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


def leading_summary(result: dict[str, Any]) -> str:
    """One honest sentence for the top of the report (spec §13.1 item 1).

    Leads with a count of the recipient's gaps, because that is the question a
    reader asks first and an index cannot answer without its scale.
    """
    sets_ = result["set_results"]
    if not sets_:
        return (
            "No admissible candidate remained, so no ranking is offered. The recipient's target ledger "
            "and every concern record are still reported in full."
        )
    ach = result.get("achievable") or {}
    ceiling = ach.get("index_ceiling")
    best = max(result["individual_results"], key=lambda x: x["goals"]["covered"], default=None)
    best_set = next((s for s in sets_ if s["material_count"] > 1), None)
    if best is None:
        return "No individual candidate could be compared."
    g = best["goals"]
    if g["total"] == 0:
        return (
            "No eligible scoring goal was found for this recipient, so no coverage index is reported. "
            "The descriptive donor comparison and the concern ledger stand on their own."
        )
    parts = [
        f"{best['person_id']} covers {g['covered']} of the {g['total']} gaps in the recipient's "
        "ledger, more than any other single material"
    ]
    if best_set and best_set["goals"]["covered"] > g["covered"]:
        parts.append(
            f"the best combination reaches {best_set['goals']['covered']} of {g['total']} "
            f"({' + '.join(best_set['member_person_ids'])})"
        )
    unreachable = ach.get("unreachable_goals") or []
    if unreachable:
        parts.append(
            f"{len(unreachable)} gap(s) are supplied by nobody, so the highest reachable index is "
            f"{ceiling:.1f} of 100 rather than 100"
        )
    return (
        "; ".join(parts)
        + ". These are counts of supply evidence in sequencing data: not engraftment, not benefit, and "
        "not a screening or clinical decision."
    )
