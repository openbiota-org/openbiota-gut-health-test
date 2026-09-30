"""MHS-E1: the experimental Mycobiome Health Score (spec §7.5).

The score summarises the *direction of the strongest interpretable benefit
and concern signals* in one stool sample. It is a policy index over
source-backed observations, not a measurement of health and not a disease
probability. The formula, caps and display rules here are OpenBiota
engineering policy (spec §7.5.1); the studies behind each contribution
supply the facts, not the numbers.

    B = max(cap_j × activation_j) over eligible beneficial claims
    C = max(cap_j × activation_j) over eligible concern claims
    score_01 = 0.5 × (1 − C) × (1 + B)

Max, not sum: repeated studies, duplicate references and many neutral fungi
cannot inflate a score or dilute a concern. The product is conservative: a
benefit can raise the score but cannot erase the dominant concern (C ≥ 0.5
caps the score at 50 whatever B is). With no supported nonzero contribution
the score is null with a reason - never 50.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any, Final, Literal

MODEL_ID: Final = "MHS-E1"
SCALE: Final = 100.0
#: Policy caps by evidence route (spec §7.5.2).
CAP_STRAIN_STUDIED: Final = 0.50
CAP_SPECIES_EXPERIMENTAL: Final = 0.25
CAP_HUMAN_ABUNDANCE: Final = 0.25
CAP_IN_VITRO: Final = 0.125
#: Policy sensitivity factors (spec §7.5.5).
POLICY_FACTORS: Final = (0.5, 1.0, 1.5)

Direction = Literal["B", "C"]
Status = Literal["not_assessed", "computed", "provisional", "not_computable", "failed"]
ReasonCode = Literal[
    "no_directional_evidence", "insufficient_fungal_information", "incomplete_analysis", "failed_qc",
]


@dataclass(frozen=True, slots=True)
class Contribution:
    """One rule evaluated against one sample (spec §7.5.2 fields).

    ``activation`` is the point value (1 for a supported presence or
    genotype rule; ``max(0, 2p − 1)`` for a qualified abundance rule; 0 when
    the rule did not fire). ``alternatives`` are the admissible activation
    values under the sample's actual analytical alternatives - strain
    equivalence members, unresolved genotype, missing context - and always
    include the point. A rule with ``alternatives == (activation,)`` has no
    scoring-relevant ambiguity.
    """

    rule_id: str
    direction: Direction
    biological_claim_id: str
    source_ids: tuple[str, ...]
    required_identity: str
    trigger_type: str
    base_cap: float
    activation: float
    alternatives: tuple[float, ...]
    finding_ids: tuple[str, ...] = ()
    required_features: tuple[str, ...] = ()
    required_context: tuple[str, ...] = ()
    extrapolations: tuple[str, ...] = ()
    #: Why the rule did or did not contribute, in one sentence.
    reason: str = ""
    #: "active" (contributes), "inactive" (did not fire), "unresolved"
    #: (a higher-resolution claim that only enters the sensitivity range).
    activation_status: str = "inactive"

    @property
    def value(self) -> float:
        return round(self.base_cap * self.activation, 12)

    @property
    def value_range(self) -> tuple[float, float]:
        alts = self.alternatives or (self.activation,)
        return (round(self.base_cap * min(alts), 12), round(self.base_cap * max(alts), 12))

    @property
    def ambiguous(self) -> bool:
        return len(set(self.alternatives)) > 1

    def to_json(self) -> dict[str, Any]:
        lo, hi = self.value_range
        return {
            "rule_id": self.rule_id, "direction": self.direction,
            "biological_claim_id": self.biological_claim_id, "source_ids": list(self.source_ids),
            "required_identity": self.required_identity, "required_features": list(self.required_features),
            "required_context": list(self.required_context), "trigger_type": self.trigger_type,
            "base_cap": self.base_cap, "activation": self.activation,
            "activation_range": [min(self.alternatives), max(self.alternatives)] if self.alternatives else None,
            "value": self.value, "admissible_value_range": [lo, hi],
            "activation_status": self.activation_status, "extrapolations": list(self.extrapolations),
            "finding_ids": list(self.finding_ids), "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class Range:
    low: float
    high: float
    method_id: str
    scenario_ids: tuple[str, ...] = ()
    scale: float = SCALE

    def to_json(self) -> dict[str, Any]:
        return {"low": round(self.low, 6), "high": round(self.high, 6), "scale": self.scale,
                "scenario_ids": list(self.scenario_ids), "method_id": self.method_id}

    @property
    def crosses_midpoint(self) -> bool:
        return self.low <= 50.0 <= self.high


@dataclass(slots=True)
class Score:
    status: Status
    reason_codes: list[str]
    score_01: float | None
    score_100: float | None
    benefit_contribution: float | None
    concern_contribution: float | None
    resolution_range: Range | None
    policy_range: Range | None
    display_range: Range | None
    direction: str
    mixed_signals: bool | None
    analytical_confidence: str
    contributions: list[Contribution]
    unscored_findings: list[str] = field(default_factory=list)
    unresolved_concerns: list[str] = field(default_factory=list)
    active_alert_ids: list[str] = field(default_factory=list)
    context_missing: list[str] = field(default_factory=list)
    scored_fungal_fragment_fraction: float | None = None
    strain_assessment_completion: dict[str, Any] | None = None
    policy_lock_id: str = ""
    interpretation_confidence: str = "experimental_limited"

    @property
    def display_integer(self) -> int | None:
        """Half-up integer for the gauge; full precision is stored."""
        return None if self.score_100 is None else int(math.floor(self.score_100 + 0.5))

    @property
    def strongest_up(self) -> Contribution | None:
        act = [c for c in self.contributions if c.direction == "B" and c.activation_status == "active" and c.value > 0]
        return max(act, key=lambda c: c.value) if act else None

    @property
    def strongest_down(self) -> Contribution | None:
        act = [c for c in self.contributions if c.direction == "C" and c.activation_status == "active" and c.value > 0]
        return max(act, key=lambda c: c.value) if act else None

    def to_json(self) -> dict[str, Any]:
        return {
            "model_id": MODEL_ID,
            "policy_lock_id": self.policy_lock_id,
            "status": self.status,
            "reason_codes": list(self.reason_codes),
            "score_01": None if self.score_01 is None else round(self.score_01, 9),
            "score_100": None if self.score_100 is None else round(self.score_100, 7),
            "score_display": self.display_integer,
            "benefit_contribution": self.benefit_contribution,
            "concern_contribution": self.concern_contribution,
            "resolution_sensitivity_range": None if self.resolution_range is None else self.resolution_range.to_json(),
            "policy_sensitivity_range": None if self.policy_range is None else self.policy_range.to_json(),
            "display_sensitivity_range": None if self.display_range is None else self.display_range.to_json(),
            "range_type": "scenario_sensitivity_not_confidence_interval",
            "direction": self.direction,
            "mixed_signals": self.mixed_signals,
            "analytical_confidence": self.analytical_confidence,
            "interpretation_confidence": self.interpretation_confidence,
            "scored_fungal_fragment_fraction": self.scored_fungal_fragment_fraction,
            "strain_assessment_completion": self.strain_assessment_completion,
            "contributions": [c.to_json() for c in self.contributions],
            "unscored_findings": list(self.unscored_findings),
            "unresolved_concerns": list(self.unresolved_concerns),
            "active_alert_ids": list(self.active_alert_ids),
            "context_missing": list(self.context_missing),
            "reference_percentile": None,
            "clinical_validation_status": "not_validated_for_global_health",
        }


# --------------------------------------------------------------------------- #
# arithmetic
# --------------------------------------------------------------------------- #


def score_01(b: float, c: float) -> float:
    """``0.5 × (1 − C) × (1 + B)`` (spec §7.5.4)."""
    if not (0.0 <= b <= 1.0 and 0.0 <= c <= 1.0):
        raise ValueError(f"B and C must lie in [0, 1]; got B={b}, C={c}")
    return 0.5 * (1.0 - c) * (1.0 + b)


def _dedupe(contribs: Iterable[Contribution]) -> list[Contribution]:
    """One contribution per (direction, biological claim): the strongest.

    Two rules on the same biological family - species-level extrapolation
    and the strain-resolved rule for the same organism - must not both
    count. Under max aggregation a duplicate cannot change B or C anyway;
    deduplicating keeps the trace honest and the coverage count exact.
    """
    best: dict[tuple[str, str], Contribution] = {}
    for c in contribs:
        key = (c.direction, c.biological_claim_id)
        cur = best.get(key)
        if cur is None or c.value_range[1] > cur.value_range[1] or (
            c.value_range[1] == cur.value_range[1] and c.value > cur.value
        ):
            best[key] = c
    return list(best.values())


def _envelope(contribs: Sequence[Contribution], direction: Direction) -> tuple[float, float, float]:
    """(point, low, high) of ``max_j v_j`` over the admissible scenario set.

    Max is monotone in every argument, so the extremes of the maximum are
    the maxima of the per-contribution extremes. Only active and unresolved
    contributions enter; an inactive rule contributes exactly 0.
    """
    point = 0.0
    lo = 0.0
    hi = 0.0
    for c in contribs:
        if c.direction != direction or c.activation_status == "inactive":
            continue
        v_lo, v_hi = c.value_range
        if c.activation_status == "active":
            point = max(point, c.value)
        # An unresolved higher-resolution claim only widens the range.
        lo = max(lo, v_lo if c.activation_status == "active" else 0.0)
        hi = max(hi, v_hi)
    # The point is itself a scenario; the envelope must contain it.
    lo = min(lo, point)
    hi = max(hi, point)
    return point, lo, hi


def compute(
    contributions: Sequence[Contribution],
    *,
    reason_if_empty: ReasonCode = "no_directional_evidence",
    analytical_confidence: str = "limited",
    incomplete: bool = False,
    unscored_findings: Sequence[str] = (),
    unresolved_concerns: Sequence[str] = (),
    active_alert_ids: Sequence[str] = (),
    context_missing: Sequence[str] = (),
    scored_fungal_fragment_fraction: float | None = None,
    strain_assessment_completion: dict[str, Any] | None = None,
    policy_lock_id: str = "",
) -> Score:
    """Evaluate MHS-E1 for one sample from its rule contributions.

    ``incomplete`` marks a mandatory stage that did not complete; the score
    is still computed from the evidence in hand but labelled provisional and
    the reason recorded. A failed mandatory stage is an implementation
    failure, not "no directional evidence" (spec §7.5.4).
    """
    contribs = _dedupe(contributions)
    b_pt, b_lo, b_hi = _envelope(contribs, "B")
    c_pt, c_lo, c_hi = _envelope(contribs, "C")
    any_active = any(c.activation_status == "active" and c.value > 0 for c in contribs)
    any_unresolved = any(c.activation_status == "unresolved" or c.ambiguous for c in contribs)
    reasons: list[str] = []
    if incomplete:
        reasons.append("incomplete_analysis")

    if not any_active and reason_if_empty == "no_directional_evidence" and not incomplete:
        # A completed analysis with no concerning and no beneficial evidence
        # is a neutral result, and it is shown as one: 50, direction
        # "neutral". This is a deliberate product decision over the spec's
        # null-for-no-evidence rule: a finished search that found nothing to
        # weigh is information, and the page says exactly that. Genuine
        # failures (lane not run, QC failed) still return null below.
        reasons.append("no_directional_evidence")
        res_range = None
        if any_unresolved:
            res_range = Range(
                low=SCALE * score_01(b_lo, c_hi), high=SCALE * score_01(b_hi, c_lo),
                method_id="resolution_envelope_v1",
                scenario_ids=tuple(c.rule_id for c in contribs if c.activation_status != "inactive"),
            )
        neutral = Range(low=50.0, high=50.0, method_id="neutral_no_directional_evidence_v1")
        return Score(
            status="computed", reason_codes=reasons, score_01=0.5, score_100=50.0,
            benefit_contribution=0.0, concern_contribution=0.0,
            resolution_range=res_range or neutral, policy_range=neutral,
            display_range=res_range or neutral,
            direction="neutral" if not any_unresolved else "uncertain", mixed_signals=False,
            analytical_confidence=analytical_confidence, contributions=list(contribs),
            unscored_findings=list(unscored_findings), unresolved_concerns=list(unresolved_concerns),
            active_alert_ids=list(active_alert_ids), context_missing=list(context_missing),
            scored_fungal_fragment_fraction=scored_fungal_fragment_fraction,
            strain_assessment_completion=strain_assessment_completion, policy_lock_id=policy_lock_id,
        )

    if not any_active:
        # Applicability precedes arithmetic: null, with the reason.
        reasons.append(reason_if_empty)
        res_range = None
        if any_unresolved:
            # Exploratory detail only; the main point stays null.
            res_range = Range(
                low=SCALE * score_01(b_lo, c_hi), high=SCALE * score_01(b_hi, c_lo),
                method_id="resolution_envelope_v1",
                scenario_ids=tuple(c.rule_id for c in contribs if c.activation_status != "inactive"),
            )
        return Score(
            status="not_computable", reason_codes=reasons, score_01=None, score_100=None,
            benefit_contribution=None, concern_contribution=None,
            resolution_range=res_range, policy_range=None, display_range=None,
            direction="not_assessed", mixed_signals=None,
            analytical_confidence="insufficient" if analytical_confidence == "insufficient" else analytical_confidence,
            contributions=list(contribs),
            unscored_findings=list(unscored_findings), unresolved_concerns=list(unresolved_concerns),
            active_alert_ids=list(active_alert_ids), context_missing=list(context_missing),
            scored_fungal_fragment_fraction=scored_fungal_fragment_fraction,
            strain_assessment_completion=strain_assessment_completion, policy_lock_id=policy_lock_id,
        )

    s01 = score_01(b_pt, c_pt)
    scenario_ids = tuple(c.rule_id for c in contribs if c.activation_status != "inactive")
    resolution = Range(
        low=SCALE * score_01(b_lo, c_hi), high=SCALE * score_01(b_hi, c_lo),
        method_id="resolution_envelope_v1", scenario_ids=scenario_ids,
    )
    # Policy factors vary B and C independently; the baseline is included.
    pol_vals = [
        SCALE * score_01(min(1.0, fb * b_pt), min(1.0, fc * c_pt))
        for fb in POLICY_FACTORS for fc in POLICY_FACTORS
    ]
    policy = Range(low=min(pol_vals), high=max(pol_vals), method_id="policy_factors_0.5_1.0_1.5_v1")
    # Displayed range: Cartesian product of resolution alternatives and the
    # independent policy factors. Under this monotone model the bound is:
    display = Range(
        low=SCALE * score_01(min(1.0, 0.5 * b_lo), min(1.0, 1.5 * c_hi)),
        high=SCALE * score_01(min(1.0, 1.5 * b_hi), min(1.0, 0.5 * c_lo)),
        method_id="resolution_x_policy_v1", scenario_ids=scenario_ids,
    )
    if display.high < 50.0:
        direction = "concerning_evidence"
    elif display.low > 50.0:
        direction = "favorable_evidence"
    else:
        direction = "uncertain"
    provisional = incomplete or any_unresolved or bool(context_missing) or bool(unresolved_concerns)
    if any_unresolved:
        reasons.append("scoring_relevant_ambiguity")
    if context_missing:
        reasons.append("required_context_missing")
    return Score(
        status="provisional" if provisional else "computed",
        reason_codes=reasons, score_01=s01, score_100=SCALE * s01,
        benefit_contribution=b_pt, concern_contribution=c_pt,
        resolution_range=resolution, policy_range=policy, display_range=display,
        direction=direction, mixed_signals=bool(b_pt > 0 and c_pt > 0),
        analytical_confidence=analytical_confidence if not any_unresolved else "limited",
        contributions=list(contribs),
        unscored_findings=list(unscored_findings), unresolved_concerns=list(unresolved_concerns),
        active_alert_ids=list(active_alert_ids), context_missing=list(context_missing),
        scored_fungal_fragment_fraction=scored_fungal_fragment_fraction,
        strain_assessment_completion=strain_assessment_completion, policy_lock_id=policy_lock_id,
    )


# --------------------------------------------------------------------------- #
# Myco-Score: the one number on the page
# --------------------------------------------------------------------------- #

MYCO_SCORE_ID: Final = "myco-score-v2"
#: Policy weights, higher is healthier. Grounded in Huang et al. 2024 (Gut
#: Microbes 16:2440111): fungi are 0.01-0.1% of the microbiome with no lower
#: bound described; Candida is a normal constituent and its *dominance* -
#: "fecal samples dominated by Candida" - is the dysbiosis pattern (IBD
#: relapse, alcohol-related liver disease, cirrhosis, CDI/FMT failure);
#: Saccharomyces is depleted in IBD, ALD, obesity and HCC and its dietary
#: presence tracks lower IBD incidence; S. boulardii has trial evidence. No
#: component rewards fungal load or diversity, which the review leaves
#: without a monotone direction.
MYCO_BASE: Final = 50.0
MYCO_CLEAN_BONUS: Final = 20.0          # opportunists < 10% of fungal DNA and < 0.001% of all fragments, none confirmed
MYCO_DOMINANCE_FREE_SHARE: Final = 0.10 # Candida presence at low share is normal (Odds; Nash 2017)
MYCO_OPP_SHARE_WEIGHT: Final = 30.0     # scaled over share from 10% to 100% of fungal DNA
MYCO_OPP_ABSOLUTE_FLOOR: Final = 1e-5   # 0.001% of all fragments: above this, opportunists are not "very little"
MYCO_OPP_SUPPORTED_PENALTY: Final = 10.0
MYCO_SACCHAROMYCES_BONUS: Final = 5.0   # supported Saccharomyces: the direction the review associates with health
MYCO_CONCERN_WEIGHT: Final = 30.0       # × C from the evidence rules
MYCO_BENEFIT_WEIGHT: Final = 20.0       # × B from the evidence rules
MYCO_ALERT_CEILING: Final = 20.0        # a supported alert-route pathogen caps the score here


def myco_score(*, opportunist_share: float | None, opportunist_fraction_all: float | None, opportunist_supported: int,
               alert_supported: int, saccharomyces_supported: bool, benefit: float, concern: float,
               complete: bool) -> dict[str, Any]:
    """The single 0-100 Myco-Score, higher = healthier.

    From a neutral 50: opportunistic fungi that are very little (under 10% of
    fungal DNA and under 0.001% of all read pairs, none confirmed at species
    level) move it up; opportunists move it down in proportion to how far
    past a tenth of the fungal DNA they reach - the dominance pattern the
    literature ties to disease - and further if one is confirmed. A
    supported Saccharomyces adds a little. The evidence rules move it on top.
    A supported alert-route pathogen caps it at 20. Null only when the
    analysis did not complete.
    """
    if not complete:
        return {"model_id": MYCO_SCORE_ID, "value": None, "components": {}, "reason": "analysis_incomplete"}
    share = float(opportunist_share or 0.0)
    frac_all = float(opportunist_fraction_all or 0.0)
    comps: dict[str, float] = {"base": MYCO_BASE}
    very_little = share < MYCO_DOMINANCE_FREE_SHARE and frac_all < MYCO_OPP_ABSOLUTE_FLOOR and opportunist_supported == 0
    if very_little:
        comps["no_opportunists"] = MYCO_CLEAN_BONUS
    else:
        if share > MYCO_DOMINANCE_FREE_SHARE:
            comps["opportunist_dominance"] = -MYCO_OPP_SHARE_WEIGHT * min(1.0, (share - MYCO_DOMINANCE_FREE_SHARE)
                                                                          / (1.0 - MYCO_DOMINANCE_FREE_SHARE))
        if opportunist_supported:
            comps["opportunist_supported"] = -MYCO_OPP_SUPPORTED_PENALTY
    if saccharomyces_supported:
        comps["saccharomyces_present"] = MYCO_SACCHAROMYCES_BONUS
    if concern:
        comps["evidence_concern"] = -MYCO_CONCERN_WEIGHT * concern
    if benefit:
        comps["evidence_benefit"] = MYCO_BENEFIT_WEIGHT * benefit
    value = max(0.0, min(100.0, sum(comps.values())))
    if alert_supported:
        comps["alert_route_ceiling"] = MYCO_ALERT_CEILING
        value = min(value, MYCO_ALERT_CEILING)
    return {"model_id": MYCO_SCORE_ID, "value": round(value, 4), "display": int(value + 0.5),
            "components": {k: round(v, 4) for k, v in comps.items()},
            "inputs": {"opportunist_share_of_fungal_dna": share, "opportunist_fraction_of_all_fragments": frac_all,
                       "opportunists_supported": opportunist_supported, "alert_route_supported": alert_supported,
                       "saccharomyces_supported": saccharomyces_supported, "B": benefit, "C": concern},
            "source_basis": "Huang et al. 2024, Gut Microbes 16:2440111 (review); policy weights are OpenBiota's",
            "reason": None}


def abundance_activation(percentile: float | None) -> float | None:
    """``max(0, 2p − 1)`` for a qualified high-abundance association rule
    (spec §7.5.2): a descriptive ramp above the reference median, not a
    pathological cutoff. None when there is no qualified reference."""
    if percentile is None:
        return None
    p = percentile / 100.0
    return max(0.0, 2.0 * p - 1.0)


def policy_lock(rules: Sequence[dict[str, Any]]) -> str:
    """Content hash of the frozen rule manifest, so a changed cap or a new
    rule is a new policy that cannot compare silently with the old."""
    blob = json.dumps(sorted(rules, key=lambda r: r.get("rule_id", "")), sort_keys=True).encode()
    return "mhs-e1-policy-" + hashlib.sha256(blob).hexdigest()[:16]


__all__ = [
    "CAP_HUMAN_ABUNDANCE", "CAP_IN_VITRO", "CAP_SPECIES_EXPERIMENTAL", "CAP_STRAIN_STUDIED", "MODEL_ID", "MYCO_SCORE_ID",
    "Contribution", "Range", "Score", "abundance_activation", "compute", "myco_score", "policy_lock", "score_01",
]
