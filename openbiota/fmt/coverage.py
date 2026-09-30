"""Supported availability and target coverage (spec §7.2).

The number this module produces is called **supported target coverage, 0–100**
and it means one thing: how much of the recipient's target ledger at least one
candidate material has evidence of being able to supply. It is not a
percentage of a microbiome restored, not an engraftment probability, not a
compatibility probability and not a clinical match score.

Three arithmetic rules carry the spec's honesty requirements:

* An unassessed target contributes 0 to the *supported* figure while keeping
  its measurement null and scenario bounds [0, 1]: no evidence earns no
  credit, but no evidence is also not proof of absence (F3, V7-024).
* The denominator is frozen across every candidate, so a sparsely measured
  donor cannot score well by having fewer targets (V7-025).
* Supply is capped at 1 per target, so excess earns nothing and ten redundant
  taxa cannot multiply one goal's weight (V7-029, V7-035).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .inputs import MaterialData
from .targets import TargetRecord


@dataclass(slots=True)
class Availability:
    """One material's evidence for supplying one target."""

    target_id: str
    material_id: str
    supported_value: float | None
    scenario_lower: float
    scenario_upper: float
    rule_id: str
    observation_ids: list[str]
    assessed: bool
    donor_value: float | None = None
    donor_percentile: float | None = None
    donor_status: str = "unassessed"
    uncertainty_reasons: list[str] = field(default_factory=list)
    note: str = ""

    @property
    def contribution(self) -> float:
        """Numeric credit: the supported value when known, otherwise zero."""
        return float(self.supported_value) if self.supported_value is not None else 0.0

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "material_id": self.material_id,
            "supported_value": self.supported_value,
            "supported_contribution": self.contribution,
            "scenario_lower": self.scenario_lower,
            "scenario_upper": self.scenario_upper,
            "assessed": self.assessed,
            "rule_id": self.rule_id,
            "observation_ids": list(self.observation_ids),
            "donor": {
                "value": self.donor_value,
                "percentile": self.donor_percentile,
                "status": self.donor_status,
            },
            "uncertainty_reasons": list(self.uncertainty_reasons),
            "note": self.note,
        }


def availability(target: TargetRecord, donor: MaterialData) -> Availability:
    """Compute a_dj under the target's own declared rule."""
    mid = donor.material.material_id
    feature = target.feature_ids[0] if target.feature_ids else ""
    obs = donor.observation(target.feature_kind, feature)

    if obs is None:
        # A feature absent from a lane that ran successfully was looked for and
        # not found: an assay-qualified non-detection, not an unassessed
        # feature. A feature from a lane that did not run is unassessed.
        lane_ran = any(o.feature_kind == target.feature_kind for o in donor.observations)
        if lane_ran and target.feature_kind == "species":
            # Before recording a non-detection, ask the second lane. This
            # applies to presence goals only: they ask "is this organism
            # here at all", which needs no reference population, so the
            # more sensitive detector should answer. Every other goal kind
            # compares against a percentile and must stay on the lane the
            # reference cohort was built with.
            #
            # Without this, a goal could be scored unreachable by every
            # candidate while the newer lane found the organism in three of
            # them - which is what happened to Adlercreutzia equolifaciens,
            # the single gap capping the achievable ceiling.
            secondary = (
                donor.secondary_presence(feature)
                if target.goal_kind == "presence" else None
            )
            if secondary is not None:
                return Availability(
                    target_id=target.target_id,
                    material_id=mid,
                    supported_value=1.0,
                    scenario_lower=1.0,
                    scenario_upper=1.0,
                    rule_id="A-PRESENCE-SECONDARY-LANE",
                    observation_ids=[],
                    assessed=True,
                    donor_status="detected_secondary_lane",
                    uncertainty_reasons=["detected_by_secondary_lane_only"],
                    note=(
                        f"not seen by the MetaPhlAn 3 lane, but the MetaPhlAn 4 lane "
                        f"reports it at {secondary:.3f}%. Presence needs no reference "
                        "population, so the more sensitive of the two detectors "
                        "answers it: MetaPhlAn 4 is four years newer and built from "
                        "roughly a million genomes against MetaPhlAn 3's seventeen "
                        "thousand. Percentile comparisons stay on the MetaPhlAn 3 "
                        "lane, which is the one the reference cohort was built with."
                    ),
                )
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=0.0,
                scenario_lower=0.0,
                scenario_upper=0.0,
                rule_id="A-PRESENCE",
                observation_ids=[],
                assessed=True,
                donor_status="not_detected",
                uncertainty_reasons=["uncalibrated_absence_hidden_carriage_possible"],
                note=(
                    "not reported by either profiler in this material: an "
                    "assay-qualified non-detection at their detection limits, not "
                    "certain absence"
                    if target.goal_kind == "presence"
                    else "not reported by the same profiler on the same database in "
                    "this material: an assay-qualified non-detection at its detection "
                    "limit, not certain absence"
                ),
            )
        return Availability(
            target_id=target.target_id,
            material_id=mid,
            supported_value=None,
            scenario_lower=0.0,
            scenario_upper=1.0,
            rule_id="A-UNASSESSED",
            observation_ids=[],
            assessed=False,
            uncertainty_reasons=["feature_outside_this_material_analytical_scope"],
            note="not assessed in this material: no supported credit, and no claim of absence",
        )

    ids = [obs.observation_id]
    trace = bool(obs.extra.get("trace"))

    if obs.status == "unassessed":
        return Availability(
            target_id=target.target_id,
            material_id=mid,
            supported_value=None,
            scenario_lower=0.0,
            scenario_upper=1.0,
            rule_id="A-UNASSESSED",
            observation_ids=ids,
            assessed=False,
            donor_status=obs.status,
            uncertainty_reasons=["unassessed_feature"],
            note="unassessed: bounds span the full range",
        )

    if target.goal_kind == "presence":
        if obs.status == "not_detected":
            # Ask the second lane before recording a non-detection. A
            # presence goal asks "is this organism here at all", which
            # needs no reference population, so the more sensitive of the
            # two profilers should answer it. Percentile goals never reach
            # this branch and stay on the MetaPhlAn 3 lane, which is the
            # one curatedMetagenomicData was built with.
            secondary = donor.secondary_presence(feature)
            if secondary is not None:
                return Availability(
                    target_id=target.target_id,
                    material_id=mid,
                    supported_value=1.0,
                    scenario_lower=1.0,
                    scenario_upper=1.0,
                    rule_id="A-PRESENCE-SECONDARY-LANE",
                    observation_ids=ids,
                    assessed=True,
                    donor_value=secondary,
                    donor_status="detected_secondary_lane",
                    uncertainty_reasons=["detected_by_secondary_lane_only"],
                    note=(
                        f"not seen by the scoring catalogue; a wider search reports it "
                        f"at {secondary:.3f}%. Presence needs no reference population, "
                        "so the more sensitive search answers it. Every percentile "
                        "stays on the scoring catalogue, because that is the one the "
                        "reference cohort exists in."
                    ),
                )
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=0.0,
                scenario_lower=0.0,
                scenario_upper=0.0,
                rule_id="A-PRESENCE",
                observation_ids=ids,
                assessed=True,
                donor_value=obs.value,
                donor_percentile=obs.percentile,
                donor_status=obs.status,
                uncertainty_reasons=["uncalibrated_absence_hidden_carriage_possible"],
                note=(
                    "assessed and not detected by either profiler: no supported supply. "
                    "Detection is uncalibrated, so this is an assay-qualified "
                    "non-detection and not certain biological absence"
                ),
            )
        return Availability(
            target_id=target.target_id,
            material_id=mid,
            supported_value=1.0,
            scenario_lower=0.0 if trace else 1.0,
            scenario_upper=1.0,
            rule_id="A-PRESENCE",
            observation_ids=ids,
            assessed=True,
            donor_value=obs.value,
            donor_percentile=obs.percentile,
            donor_status=obs.status,
            uncertainty_reasons=["trace_call_near_detection_limit"] if trace else [],
            note=(
                "trace-level call: the lower scenario drops it entirely"
                if trace
                else "detected in the same species lane"
            ),
        )

    if target.goal_kind == "percentile_range":
        lower = float(target.desired_state.get("lower") or 0.0)
        if obs.status == "not_detected":
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=0.0,
                scenario_lower=0.0,
                scenario_upper=0.0,
                rule_id="A-RANGE-PCT",
                observation_ids=ids,
                assessed=True,
                donor_value=obs.value,
                donor_percentile=obs.percentile,
                donor_status=obs.status,
                uncertainty_reasons=["uncalibrated_absence_hidden_carriage_possible"],
                note=(
                    "assessed and not detected, so no supply evidence. A zero abundance still carries a "
                    "cohort percentile; that percentile is not evidence of supply"
                ),
            )
        pct = obs.percentile
        if pct is None:
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=None,
                scenario_lower=0.0,
                scenario_upper=1.0,
                rule_id="A-RANGE-PCT",
                observation_ids=ids,
                assessed=False,
                donor_value=obs.value,
                donor_status=obs.status,
                uncertainty_reasons=["no_reference_percentile_for_this_feature"],
                note="no reference percentile available in this material",
            )
        value = min(1.0, pct / lower) if lower > 0 else (1.0 if pct > 0 else 0.0)
        return Availability(
            target_id=target.target_id,
            material_id=mid,
            supported_value=value,
            scenario_lower=0.0 if trace else value,
            scenario_upper=1.0 if trace else value,
            rule_id="A-RANGE-PCT",
            observation_ids=ids,
            assessed=True,
            donor_value=obs.value,
            donor_percentile=pct,
            donor_status=obs.status,
            uncertainty_reasons=["trace_call_near_detection_limit"] if trace else [],
            note=(
                f"donor reads at the {pct:.0f}th percentile against the same reference; credit is "
                f"capped at the {lower:.0f}th"
            ),
        )

    if target.goal_kind == "graded_capacity":
        lower = target.desired_state.get("lower")
        if not isinstance(lower, (int, float)) or lower <= 0:
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=None,
                scenario_lower=0.0,
                scenario_upper=1.0,
                rule_id="A-INVALID-THRESHOLD",
                observation_ids=ids,
                assessed=False,
                donor_value=obs.value,
                donor_status=obs.status,
                uncertainty_reasons=["threshold_not_finite_and_positive"],
                note="graded supply needs a finite positive threshold; none was available",
            )
        if obs.value is None:
            return Availability(
                target_id=target.target_id,
                material_id=mid,
                supported_value=None,
                scenario_lower=0.0,
                scenario_upper=1.0,
                rule_id="A-GRADED",
                observation_ids=ids,
                assessed=False,
                donor_status=obs.status,
                uncertainty_reasons=["no_value_for_this_panel"],
                note="panel not quantified in this material",
            )
        value = min(1.0, float(obs.value) / float(lower))
        return Availability(
            target_id=target.target_id,
            material_id=mid,
            supported_value=value,
            scenario_lower=value,
            scenario_upper=value,
            rule_id="A-GRADED",
            observation_ids=ids,
            assessed=True,
            donor_value=obs.value,
            donor_percentile=obs.percentile,
            donor_status=obs.status,
            note=(
                f"{obs.value:.3g} {obs.unit} against a reference lower quartile of {float(lower):.3g}; "
                "a graded supply index, not a dose"
            ),
        )

    return Availability(
        target_id=target.target_id,
        material_id=mid,
        supported_value=None,
        scenario_lower=0.0,
        scenario_upper=1.0,
        rule_id="A-NO-CORE-RULE",
        observation_ids=ids,
        assessed=False,
        donor_value=obs.value,
        donor_percentile=obs.percentile,
        donor_status=obs.status,
        uncertainty_reasons=["no_implemented_core_scoring_rule_for_this_goal_kind"],
        note="shown for comparison; outside the scored denominator",
    )


@dataclass(slots=True)
class CoverageResult:
    """C(S) with its frozen denominator and honest weight accounting."""

    coverage: float | None
    lower: float | None
    upper: float | None
    reason: str | None
    assessed_goal_weight: float
    supported_goal_weight: float
    unassessed_goal_weight: float
    per_target: dict[str, float]
    per_group: dict[str, float]
    denominator_target_ids: list[str]

    def to_json(self) -> dict[str, Any]:
        return {
            "coverage": self.coverage,
            "scenario_lower": self.lower,
            "scenario_upper": self.upper,
            "scenario_width": (
                None if self.coverage is None or self.lower is None or self.upper is None
                else round(self.upper - self.lower, 10)
            ),
            "interval_meaning": "scenario bounds from analytical uncertainty, not a confidence interval",
            "reason": self.reason,
            "assessed_goal_weight": self.assessed_goal_weight,
            "supported_goal_weight": self.supported_goal_weight,
            "unassessed_goal_weight": self.unassessed_goal_weight,
            "per_group": dict(sorted(self.per_group.items())),
            "denominator_target_count": len(self.denominator_target_ids),
            "scale": "supported target coverage, 0-100",
        }


def coverage(
    *,
    targets: list[TargetRecord],
    alpha: dict[str, float],
    beta: dict[str, float],
    table: dict[tuple[str, str], Availability],
    members: list[str],
) -> CoverageResult:
    """C(S) = 100 · Σ_g α_g Σ_j β_j max_{d∈S} a_dj (spec §7.2)."""
    scoring = [t for t in targets if t.is_positive_scoring]
    if not scoring or not alpha:
        return CoverageResult(
            coverage=None,
            lower=None,
            upper=None,
            reason="no_active_scoring_targets",
            assessed_goal_weight=0.0,
            supported_goal_weight=0.0,
            unassessed_goal_weight=0.0,
            per_target={},
            per_group={},
            denominator_target_ids=[],
        )

    per_target: dict[str, float] = {}
    per_group: dict[str, float] = dict.fromkeys(alpha, 0.0)
    total = low = high = 0.0
    assessed_w = supported_w = unassessed_w = 0.0

    for t in scoring:
        w = alpha[t.counting_group] * beta[t.target_id]
        cells = [table[(t.target_id, m)] for m in members if (t.target_id, m) in table]
        point = max((c.contribution for c in cells), default=0.0)
        lo = max((c.scenario_lower for c in cells), default=0.0)
        hi = max((c.scenario_upper for c in cells), default=0.0) if cells else 0.0
        fully_assessed = bool(cells) and (
            all(c.assessed for c in cells) or any(c.contribution >= 1.0 for c in cells)
        )
        per_target[t.target_id] = point
        per_group[t.counting_group] = per_group.get(t.counting_group, 0.0) + beta[t.target_id] * point
        total += w * point
        low += w * lo
        high += w * hi
        if fully_assessed:
            assessed_w += w
        else:
            unassessed_w += w
        if point > 0:
            supported_w += w

    return CoverageResult(
        coverage=round(100.0 * total, 10),
        lower=round(100.0 * low, 10),
        upper=round(100.0 * high, 10),
        reason=None,
        assessed_goal_weight=round(assessed_w, 10),
        supported_goal_weight=round(supported_w, 10),
        unassessed_goal_weight=round(unassessed_w, 10),
        per_target=per_target,
        per_group={g: round(100.0 * alpha[g] * v, 10) for g, v in per_group.items()},
        denominator_target_ids=[t.target_id for t in scoring],
    )


def build_table(
    targets: list[TargetRecord], donors: list[MaterialData]
) -> dict[tuple[str, str], Availability]:
    """Every (target, material) availability cell, computed once."""
    table: dict[tuple[str, str], Availability] = {}
    for t in targets:
        for d in donors:
            table[(t.target_id, d.material.material_id)] = availability(t, d)
    return table
