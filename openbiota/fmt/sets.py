"""Donor-set enumeration, marginal value, Pareto fronts and sensitivity.

Spec §9 and §7.5. The search is over *material units*, set identity is the
sorted member list, and both `material_count` and `unique_donor_count` enter
the comparison so replicate lots cannot evade the parsimony objective.

What the set arithmetic means: the availability ceiling A_j(S) = max_{d∈S} a_dj
answers "does at least one member have evidence of supplying this target". It
is not a pooled concentration, not a final community composition, not an
engraftment probability and not a persistence forecast (§9.2).
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

from .coverage import Availability, CoverageResult, coverage
from .targets import TargetRecord


@dataclass(slots=True)
class Candidate:
    """A singleton or a set, with everything needed to rank and explain it."""

    candidate_id: str
    member_material_ids: list[str]
    member_person_ids: list[str]
    labels: list[str]
    coverage: CoverageResult
    marginals: dict[str, Any] = field(default_factory=dict)
    unique_targets: dict[str, list[str]] = field(default_factory=dict)
    redundant_target_ids: list[str] = field(default_factory=list)
    pareto_front: int | None = None
    display_order: int | None = None
    stability: dict[str, Any] = field(default_factory=dict)

    @property
    def material_count(self) -> int:
        return len(self.member_material_ids)

    @property
    def unique_donor_count(self) -> int:
        return len(set(self.member_person_ids))

    @property
    def is_set(self) -> bool:
        return self.material_count > 1


def set_id(member_material_ids: Iterable[str]) -> str:
    """Permutation-invariant identity (V7-057). Order lives in a separate object."""
    return "+".join(sorted(member_material_ids))


def feasible_subsets(
    materials: list[str],
    *,
    person_of: dict[str, str],
    max_donors: int | None,
    max_materials: int | None,
    max_subsets: int,
) -> tuple[list[tuple[str, ...]], dict[str, Any]]:
    """Enumerate admissible non-empty subsets, singletons always included.

    The budget is raised to at least the candidate count so every individual is
    always evaluated, and that adjustment is reported (V7-118).
    """
    n = len(materials)
    theoretical = (2 ** n - 1) if n else 0
    limit_k = min(max_materials or n, n)
    feasible_all: list[tuple[str, ...]] = []
    for k in range(1, limit_k + 1):
        for combo in itertools.combinations(sorted(materials), k):
            if max_donors is not None and len({person_of[m] for m in combo}) > max_donors:
                continue
            feasible_all.append(combo)
    effective = max(max_subsets, n)
    adjusted = effective != max_subsets
    complete = len(feasible_all) <= effective
    if complete:
        chosen = feasible_all
        omitted = None
    else:
        # Singletons first, then smallest sets, so the archive always contains
        # every individual and the most parsimonious combinations.
        chosen = feasible_all[:effective]
        omitted = (
            f"{len(feasible_all) - len(chosen)} larger feasible subsets were not evaluated because the "
            f"effective budget is {effective}"
        )
    return chosen, {
        "candidate_material_count": n,
        "unique_donor_count": len(set(person_of.values())),
        "theoretical_nonempty_subsets": theoretical,
        "feasible_nonempty_subsets": len(feasible_all),
        "max_subsets": max_subsets,
        "effective_max_subsets": effective,
        "budget_adjusted_to_include_every_singleton": adjusted,
        "evaluated_subsets": len(chosen),
        "search_complete": complete,
        "algorithm": "exhaustive" if complete else "size_ordered_truncated",
        "optimality_certificate": (
            "exhaustive enumeration of every feasible non-empty subset under the stated constraints"
            if complete
            else None
        ),
        "omitted_reason": omitted,
    }


def evaluate_candidates(
    *,
    targets: list[TargetRecord],
    alpha: dict[str, float],
    beta: dict[str, float],
    table: dict[tuple[str, str], Availability],
    subsets: list[tuple[str, ...]],
    person_of: dict[str, str],
    label_of: dict[str, str],
) -> list[Candidate]:
    """Coverage, marginal gains, unique contributions and redundancy per subset."""
    out: list[Candidate] = []
    singleton_cov: dict[str, float] = {}
    for combo in subsets:
        if len(combo) == 1:
            cov = coverage(targets=targets, alpha=alpha, beta=beta, table=table, members=list(combo))
            singleton_cov[combo[0]] = cov.coverage or 0.0

    global_best = max(singleton_cov, key=lambda m: singleton_cov[m], default=None)

    for combo in subsets:
        members = list(combo)
        cov = coverage(targets=targets, alpha=alpha, beta=beta, table=table, members=members)
        cand = Candidate(
            candidate_id=set_id(members),
            member_material_ids=sorted(members),
            member_person_ids=[person_of[m] for m in sorted(members)],
            labels=[label_of[m] for m in sorted(members)],
            coverage=cov,
        )
        # Leave-one-member-out marginals; shared targets are counted once
        # because coverage itself takes a max per target.
        loo: dict[str, Any] = {}
        for m in members:
            rest = [x for x in members if x != m]
            rest_cov = (
                coverage(targets=targets, alpha=alpha, beta=beta, table=table, members=rest).coverage or 0.0
                if rest
                else 0.0
            )
            loo[m] = round((cov.coverage or 0.0) - rest_cov, 10)
        best_member = max(members, key=lambda m: singleton_cov.get(m, 0.0)) if members else None
        cand.marginals = {
            "leave_one_member_out": loo,
            "best_member_singleton_id": best_member,
            "best_member_singleton_coverage": singleton_cov.get(best_member or "", 0.0),
            "gain_over_best_member_singleton": round(
                (cov.coverage or 0.0) - singleton_cov.get(best_member or "", 0.0), 10
            ),
            "global_best_singleton_id": global_best,
            "global_best_singleton_coverage": singleton_cov.get(global_best or "", 0.0),
            "gain_over_global_best_singleton": round(
                (cov.coverage or 0.0) - singleton_cov.get(global_best or "", 0.0), 10
            ),
            "empty_set_coverage": 0.0,
            "comparator_note": (
                "both comparators use the identical target denominator and goal policy"
            ),
        }
        unique: dict[str, list[str]] = {m: [] for m in members}
        redundant: list[str] = []
        for t in targets:
            if not t.is_positive_scoring:
                continue
            supporters = [m for m in members if table[(t.target_id, m)].contribution > 0]
            if len(supporters) == 1:
                unique[supporters[0]].append(t.target_id)
            elif len(supporters) > 1:
                redundant.append(t.target_id)
        cand.unique_targets = unique
        cand.redundant_target_ids = redundant
        out.append(cand)
    return out


#: The objective vector of §7.5. Hazards are deliberately absent: they are
#: shown alongside, never traded against benefit.
OBJECTIVES = (
    ("scenario_lower", "max"),
    ("supported_coverage", "max"),
    ("assessed_goal_weight", "max"),
    ("scenario_width", "min"),
    ("unique_donor_count", "min"),
    ("material_count", "min"),
)

TOLERANCE = 1e-9


def _vector(c: Candidate) -> dict[str, float]:
    cov = c.coverage
    return {
        "scenario_lower": cov.lower or 0.0,
        "supported_coverage": cov.coverage or 0.0,
        "assessed_goal_weight": cov.assessed_goal_weight,
        "scenario_width": ((cov.upper or 0.0) - (cov.lower or 0.0)),
        "unique_donor_count": float(c.unique_donor_count),
        "material_count": float(c.material_count),
    }


def dominates(a: Candidate, b: Candidate) -> bool:
    """True when a is no worse on every objective and better on at least one."""
    va, vb = _vector(a), _vector(b)
    better = False
    for key, sense in OBJECTIVES:
        x, y = va[key], vb[key]
        if sense == "max":
            if x < y - TOLERANCE:
                return False
            if x > y + TOLERANCE:
                better = True
        else:
            if x > y + TOLERANCE:
                return False
            if x < y - TOLERANCE:
                better = True
    return better


def pareto_fronts(candidates: list[Candidate]) -> None:
    """Assign 1-based front indices in place."""
    remaining = list(candidates)
    front = 1
    while remaining:
        current = [c for c in remaining if not any(dominates(o, c) for o in remaining if o is not c)]
        if not current:  # pragma: no cover - cycles are impossible with this order
            current = list(remaining)
        for c in current:
            c.pareto_front = front
        remaining = [c for c in remaining if c.pareto_front is None]
        front += 1


def display_order(candidates: list[Candidate]) -> list[Candidate]:
    """The deterministic order of §7.5, with stable candidate ID as final key."""
    ordered = sorted(
        candidates,
        key=lambda c: (
            c.pareto_front or 10**6,
            -(c.coverage.lower or 0.0),
            -(c.coverage.coverage or 0.0),
            -c.coverage.assessed_goal_weight,
            c.unique_donor_count,
            c.material_count,
            c.candidate_id,
        ),
    )
    for i, c in enumerate(ordered, start=1):
        c.display_order = i
    return ordered


def sensitivity(
    *,
    targets: list[TargetRecord],
    table: dict[tuple[str, str], Availability],
    subsets: list[tuple[str, ...]],
    person_of: dict[str, str],
    label_of: dict[str, str],
) -> dict[str, Any]:
    """Declared scenarios, and how often each candidate leads (spec §12.2, §12.3).

    A rank frequency over these scenarios is not a probability that a donor
    will work: it says how sensitive the ordering is to the stated choices.
    """
    from .targets import goal_weights

    scenarios: list[dict[str, Any]] = []

    def run(name: str, subset_targets: list[TargetRecord], note: str, mutate: Any = None) -> None:
        if not subset_targets:
            scenarios.append({"scenario_id": name, "note": note, "coverage": None,
                              "reason": "no_active_scoring_targets", "leaders": []})
            return
        local_table = dict(table)
        if mutate is not None:
            local_table = mutate(local_table)
        a, b = goal_weights(subset_targets)
        cands = evaluate_candidates(
            targets=subset_targets, alpha=a, beta=b, table=local_table,
            subsets=subsets, person_of=person_of, label_of=label_of,
        )
        pareto_fronts(cands)
        leaders = sorted(c.candidate_id for c in cands if c.pareto_front == 1)
        scenarios.append(
            {
                "scenario_id": name,
                "note": note,
                "scoring_target_count": len([t for t in subset_targets if t.is_positive_scoring]),
                "leaders": leaders,
                "coverage_by_candidate": {
                    c.candidate_id: c.coverage.coverage for c in sorted(cands, key=lambda x: x.candidate_id)
                },
            }
        )

    run("research_goals", targets, "every active or conditional goal, the default view")
    run(
        "reference_only",
        [t for t in targets if t.origin == "reference_supported"],
        "reference-supported goals only",
    )
    run(
        "no_indication_hypotheses",
        [t for t in targets if t.origin != "indication_hypothesis"],
        "drops indication hypotheses, keeping reference-supported and user goals",
    )

    def drop_trace(tbl: dict[tuple[str, str], Availability]) -> dict[tuple[str, str], Availability]:
        out = {}
        for key, av in tbl.items():
            if "trace_call_near_detection_limit" in av.uncertainty_reasons:
                out[key] = Availability(
                    target_id=av.target_id, material_id=av.material_id,
                    supported_value=av.scenario_lower, scenario_lower=av.scenario_lower,
                    scenario_upper=av.scenario_upper, rule_id=av.rule_id + "-NO-TRACE",
                    observation_ids=av.observation_ids, assessed=av.assessed,
                    donor_value=av.donor_value, donor_percentile=av.donor_percentile,
                    donor_status=av.donor_status,
                    uncertainty_reasons=av.uncertainty_reasons,
                    note="trace calls excluded in this scenario",
                )
            else:
                out[key] = av
        return out

    run("trace_calls_excluded", targets, "trace-level detections do not count as supply", drop_trace)

    # Equal group weights give a one-target group the same say as a
    # twelve-target group. The alternative is proportional weighting, so the
    # report shows whether the ordering depends on that choice (§12.2).
    scoring = [t for t in targets if t.is_positive_scoring]
    if scoring:
        sizes: dict[str, int] = {}
        for t in scoring:
            sizes[t.counting_group] = sizes.get(t.counting_group, 0) + 1
        total = sum(sizes.values())
        alpha_prop = {g: n / total for g, n in sizes.items()}
        beta_eq = {t.target_id: 1.0 / sizes[t.counting_group] for t in scoring}
        cands = evaluate_candidates(
            targets=targets, alpha=alpha_prop, beta=beta_eq, table=table,
            subsets=subsets, person_of=person_of, label_of=label_of,
        )
        pareto_fronts(cands)
        scenarios.append(
            {
                "scenario_id": "group_weights_by_size",
                "note": "goal-group weights proportional to the number of targets in each group",
                "scoring_target_count": len(scoring),
                "leaders": sorted(c.candidate_id for c in cands if c.pareto_front == 1),
                "coverage_by_candidate": {
                    c.candidate_id: c.coverage.coverage for c in sorted(cands, key=lambda x: x.candidate_id)
                },
            }
        )

    counted = [s for s in scenarios if s.get("leaders")]
    freq: dict[str, int] = {}
    for s in counted:
        for cid in s["leaders"]:
            freq[cid] = freq.get(cid, 0) + 1
    return {
        "scenarios": scenarios,
        "scenario_count": len(counted),
        "leading_front_frequency": dict(sorted(freq.items())),
        "meaning": (
            "the share of declared sensitivity scenarios in which a candidate is on the leading "
            "objective front. This is not a clinical success probability and not a bootstrap."
        ),
    }
