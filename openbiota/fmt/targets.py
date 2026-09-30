"""The recipient target ledger (spec §7.1, §5.3).

Built **once, from the recipient alone**, before any donor is examined. If the
donor list could change what the recipient is said to need, every ranking
would be circular (V7-019, V7-114).

Four origins are kept apart: a reference-supported low or absent feature, an
explicit research goal, an indication hypothesis, and a report candidate.
Neutral or unknown-polarity features are *context only* by default and stay
out of the benefit denominator until a user activates them (§7.1 rule 5).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Final, Literal

from openbiota import coherence

from .identity import content_id
from .inputs import MaterialData, Observation

Direction = Literal["supply", "support_range", "avoid_introduction", "monitor"]
Origin = Literal["reference_supported", "explicit_user_goal", "indication_hypothesis", "report_candidate"]
Status = Literal["active", "conditional", "context_only", "invalid"]

#: A species must be this common in the reference cohort before its absence in
#: the recipient is treated as a candidate restoration target. Below it, the
#: species is ordinary variation between healthy people, not a gap.
PREVALENCE_FOR_PRESENCE_TARGET: Final = 0.50

#: Percentile below which a present feature counts as low against the
#: reference distribution. The reference's own lower quartile, not a
#: biological threshold.
LOW_PERCENTILE: Final = 25.0

#: Percentile at or above which a feature counts as high.
HIGH_PERCENTILE: Final = 75.0


@dataclass(slots=True)
class TargetRecord:
    """One thing the recipient's data suggests a donor might supply."""

    target_id: str
    feature_group_id: str
    feature_kind: str
    feature_ids: list[str]
    label: str
    direction: Direction
    goal_kind: str
    origin: Origin
    target_status: Status
    recipient_observation_ids: list[str]
    reference_id: str | None
    desired_state: dict[str, Any]
    clinical_endpoint: None
    mechanism_source_ids: list[str]
    priority_class: str
    counting_group: str
    within_domain_weight: float = 1.0
    uncertainty_reasons: list[str] = field(default_factory=list)
    rule_id: str = ""
    recipient_value: float | None = None
    recipient_percentile: float | None = None
    recipient_status: str = "unassessed"
    notes: list[str] = field(default_factory=list)

    @property
    def is_positive_scoring(self) -> bool:
        """Whether this target enters the positive-coverage denominator (§7.3)."""
        return self.direction in ("supply", "support_range") and self.target_status in ("active", "conditional")

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "feature_group_id": self.feature_group_id,
            "feature_kind": self.feature_kind,
            "feature_ids": list(self.feature_ids),
            "label": self.label,
            "direction": self.direction,
            "goal_kind": self.goal_kind,
            "origin": self.origin,
            "target_status": self.target_status,
            "recipient_observation_ids": list(self.recipient_observation_ids),
            "reference_id": self.reference_id,
            "desired_state": dict(self.desired_state),
            "clinical_endpoint": None,
            "mechanism_source_ids": list(self.mechanism_source_ids),
            "priority_class": self.priority_class,
            "counting_group": self.counting_group,
            "within_domain_weight": self.within_domain_weight,
            "uncertainty_reasons": list(self.uncertainty_reasons),
            "rule_id": self.rule_id,
            "recipient": {
                "value": self.recipient_value,
                "percentile": self.recipient_percentile,
                "status": self.recipient_status,
            },
            "notes": list(self.notes),
            "in_positive_denominator": self.is_positive_scoring,
        }


def _counting_group_for_guild(guild: str) -> str:
    """One underlying signal, one vote (§7.1 rule 8).

    A guild and the gene panel that measures the same capacity share a
    counting group, using the pipeline's own capacity pairs. Ten butyrate taxa
    plus a butyrate gene panel therefore cannot multiply butyrate's importance.
    """
    pair = coherence.PAIR_BY_GROUP.get(guild)
    return f"capacity:{pair.panel}" if pair is not None else f"guild:{guild}"


def _counting_group_for_panel(panel: str) -> str:
    pair = coherence.PAIR_BY_PANEL.get(panel)
    return f"capacity:{panel}" if pair is not None else f"panel:{panel}"


def _priority(percentile: float | None, prevalence: float | None, detected: bool) -> str:
    if not detected and (prevalence or 0) >= 0.9:
        return "high"
    if percentile is not None and percentile < 5:
        return "high"
    return "moderate"


def build_targets(
    recipient: MaterialData,
    *,
    indication: str | None = None,
    custom_goals: list[dict[str, Any]] | None = None,
) -> list[TargetRecord]:
    """Compile the ledger. Deterministic order, stable IDs."""
    targets: list[TargetRecord] = []
    seen: set[str] = set()

    def push(t: TargetRecord) -> None:
        if t.target_id in seen:
            return
        seen.add(t.target_id)
        targets.append(t)

    _guild_species_targets(recipient, push)
    _panel_targets(recipient, push)
    if indication:
        _indication_targets(recipient, indication, push)
    _context_targets(recipient, push)
    for goal in custom_goals or []:
        _custom_goal(recipient, goal, push)

    # A target may not sit in two active counting groups (§7.2). Presence and
    # range targets for the same feature collapse to the stronger one.
    by_feature: dict[tuple[str, str], TargetRecord] = {}
    ordered: list[TargetRecord] = []
    for t in targets:
        key = (t.feature_kind, t.feature_ids[0] if t.feature_ids else t.target_id)
        if t.is_positive_scoring and key in by_feature:
            by_feature[key].notes.append(
                f"a second candidate goal for this feature ({t.goal_kind}, {t.origin}) was folded in "
                "so the feature is counted once"
            )
            by_feature[key].mechanism_source_ids.extend(t.mechanism_source_ids)
            continue
        if t.is_positive_scoring:
            by_feature[key] = t
        ordered.append(t)
    return ordered


def _guild_species_targets(recipient: MaterialData, push: Any) -> None:
    """Species targets from curated guilds, which are the only curated polarity.

    Species outside every curated guild have no sourced direction, so they are
    context only however striking the number looks.
    """
    for guild_obs in sorted(
        (o for o in recipient.observations if o.feature_kind == "guild"), key=lambda o: o.source_id
    ):
        polarity = str(guild_obs.extra.get("higher_means") or "unknown")
        group_label = str(guild_obs.extra.get("label") or guild_obs.source_id)
        counting = _counting_group_for_guild(guild_obs.source_id)
        for member in guild_obs.extra.get("members") or []:
            name = str(member.get("species"))
            sp = recipient.observation("species", name)
            detected = bool(member.get("detected"))
            pct = member.get("percentile")
            prev = member.get("cohort_prevalence")
            obs_ids = [sp.observation_id] if sp else [guild_obs.observation_id]
            if polarity != "favourable":
                # Adverse or unknown polarity: never an automatic supply goal,
                # and a high value is not automatically an avoidance goal
                # (V7-031, V7-116).
                continue
            if not member.get("in_catalogue", True):
                continue
            if not detected and (prev or 0) >= PREVALENCE_FOR_PRESENCE_TARGET:
                push(
                    TargetRecord(
                        target_id=content_id("target", {"k": "presence", "f": name}),
                        feature_group_id=counting,
                        feature_kind="species",
                        feature_ids=[name],
                        label=f"{name.replace('_', ' ')} present",
                        direction="supply",
                        goal_kind="presence",
                        origin="reference_supported",
                        target_status="conditional",
                        recipient_observation_ids=obs_ids,
                        reference_id=guild_obs.reference_bundle_id,
                        desired_state={
                            "kind": "presence",
                            "statement": "detected at any abundance in the same species lane",
                            "reference_prevalence": prev,
                        },
                        clinical_endpoint=None,
                        mechanism_source_ids=[f"guild:{guild_obs.source_id}"],
                        priority_class=_priority(pct, prev, detected),
                        counting_group=counting,
                        uncertainty_reasons=[
                            "uncalibrated_absence",
                            "population_reference_not_matched_to_recipient",
                        ],
                        rule_id="T-PRESENCE-1",
                        recipient_value=sp.value if sp else None,
                        recipient_percentile=pct,
                        recipient_status="not_detected",
                        notes=[
                            f"member of the curated guild '{group_label}' (favourable direction), "
                            f"present in {(prev or 0) * 100:.0f}% of the reference cohort and not "
                            "detected here"
                        ],
                    )
                )
            elif detected and pct is not None and pct < LOW_PERCENTILE:
                push(
                    TargetRecord(
                        target_id=content_id("target", {"k": "range", "f": name}),
                        feature_group_id=counting,
                        feature_kind="species",
                        feature_ids=[name],
                        label=f"{name.replace('_', ' ')} toward the typical range",
                        direction="support_range",
                        goal_kind="percentile_range",
                        origin="reference_supported",
                        target_status="conditional",
                        recipient_observation_ids=obs_ids,
                        reference_id=guild_obs.reference_bundle_id,
                        desired_state={
                            "kind": "percentile_range",
                            "lower": LOW_PERCENTILE,
                            "upper": 100.0,
                            "unit": "reference_percentile",
                            "statement": (
                                "a donor reading at or above the reference lower quartile for this "
                                "species, in the same lane and against the same reference"
                            ),
                        },
                        clinical_endpoint=None,
                        mechanism_source_ids=[f"guild:{guild_obs.source_id}"],
                        priority_class=_priority(pct, prev, detected),
                        counting_group=counting,
                        uncertainty_reasons=[
                            "relative_abundance_not_absolute_load",
                            "population_reference_not_matched_to_recipient",
                        ],
                        rule_id="T-RANGE-1",
                        recipient_value=sp.value if sp else member.get("percent"),
                        recipient_percentile=pct,
                        recipient_status="quantified",
                        notes=[
                            f"member of the curated guild '{group_label}' reading at the "
                            f"{pct:.0f}th percentile of the reference cohort"
                        ],
                    )
                )


def _panel_targets(recipient: MaterialData, push: Any) -> None:
    """Gene-capacity targets, the only lane with sourced quartiles."""
    for obs in sorted(
        (o for o in recipient.observations if o.feature_kind == "gene_panel"), key=lambda o: o.source_id
    ):
        polarity = str(obs.extra.get("higher_means") or "unknown")
        pct = obs.percentile
        counting = _counting_group_for_panel(obs.source_id)
        label = str(obs.extra.get("metabolite") or obs.source_id)
        if pct is None or obs.reference_p25 is None:
            continue
        if polarity == "favourable" and pct < LOW_PERCENTILE and obs.reference_p25 > 0:
            push(
                TargetRecord(
                    target_id=content_id("target", {"k": "panel", "f": obs.source_id}),
                    feature_group_id=counting,
                    feature_kind="gene_panel",
                    feature_ids=[obs.source_id],
                    label=f"{label} gene capacity toward the typical range",
                    direction="support_range",
                    goal_kind="graded_capacity",
                    origin="reference_supported",
                    target_status="conditional",
                    recipient_observation_ids=[obs.observation_id],
                    reference_id=obs.reference_bundle_id,
                    desired_state={
                        "kind": "graded_supply",
                        "lower": obs.reference_p25,
                        "upper": obs.reference_p75,
                        "unit": obs.unit,
                        "statement": (
                            "donor capacity relative to the reference cohort's lower quartile "
                            f"({obs.reference_p25:.3g} {obs.unit}), capped at 1"
                        ),
                    },
                    clinical_endpoint=None,
                    mechanism_source_ids=[f"panel:{obs.source_id}"],
                    priority_class=_priority(pct, None, True),
                    counting_group=counting,
                    uncertainty_reasons=[
                        "genetic_capacity_not_metabolite_concentration",
                        "population_reference_not_matched_to_recipient",
                    ],
                    rule_id="T-PANEL-1",
                    recipient_value=obs.value,
                    recipient_percentile=pct,
                    recipient_status=obs.status,
                    notes=[
                        f"recipient reads at the {pct:.0f}th percentile of {obs.extra.get('reference_n')} "
                        "reference stools for this panel"
                    ],
                )
            )
        elif polarity == "adverse" and pct >= HIGH_PERCENTILE:
            push(
                TargetRecord(
                    target_id=content_id("target", {"k": "avoid", "f": obs.source_id}),
                    feature_group_id=counting,
                    feature_kind="gene_panel",
                    feature_ids=[obs.source_id],
                    label=f"{label} gene capacity already high — avoid adding more",
                    direction="avoid_introduction",
                    goal_kind="excess_context",
                    origin="reference_supported",
                    target_status="conditional",
                    recipient_observation_ids=[obs.observation_id],
                    reference_id=obs.reference_bundle_id,
                    desired_state={
                        "kind": "avoid_excess",
                        "statement": (
                            "shown for comparison only; a donor lacking this capacity is not "
                            "evidence that it would be reduced in the recipient"
                        ),
                    },
                    clinical_endpoint=None,
                    mechanism_source_ids=[f"panel:{obs.source_id}"],
                    priority_class="moderate",
                    counting_group=counting,
                    uncertainty_reasons=["no_eradication_credit_is_computed"],
                    rule_id="T-AVOID-1",
                    recipient_value=obs.value,
                    recipient_percentile=pct,
                    recipient_status=obs.status,
                    notes=["adverse direction is sourced from the panel's own reference direction"],
                )
            )
        elif polarity == "favourable" and pct >= HIGH_PERCENTILE:
            push(
                _context(
                    obs,
                    counting,
                    f"{label} gene capacity already at the {pct:.0f}th percentile",
                    "already high: no benefit is awarded for adding more (§7.1 rule 7)",
                )
            )


def _indication_targets(recipient: MaterialData, indication: str, push: Any) -> None:
    """Mechanistic targets from the recipient's own scored pattern for the indication.

    Uses the pipeline's resolved feature rows, so every direction comes from
    the study module that published it. A pattern weight is not a therapeutic
    weight, so these stay `conditional` and carry their module as the source
    (§7.1 rule 4).
    """
    profiles = ((recipient.results.get("profile_similarity") or {}).get("profiles") or {})
    prof = profiles.get(indication)
    if not isinstance(prof, dict):
        return
    for mod_name, module in sorted((prof.get("modules") or {}).items()):
        if not module.get("bound"):
            continue
        for feat in module.get("features") or []:
            if not feat.get("measured"):
                continue
            direction = str(feat.get("direction") or "")
            engine = str(feat.get("engine") or "")
            name = str(feat.get("feature") or feat.get("name") or "")
            ref_pct = feat.get("reference_percentile")
            low_in_case = direction in ("decreased", "lower_in_case")
            if not (low_in_case and isinstance(ref_pct, (int, float)) and ref_pct < LOW_PERCENTILE):
                continue
            kind = "gene_panel" if engine == "diamond" else "species"
            counting = (
                _counting_group_for_panel(name)
                if kind == "gene_panel"
                else _guild_counting_for_species(recipient, name)
            )
            obs = recipient.observation(kind, name)
            push(
                TargetRecord(
                    target_id=content_id("target", {"k": "indication", "f": name, "i": indication}),
                    feature_group_id=counting,
                    feature_kind=kind,
                    feature_ids=[name],
                    label=f"{name.replace('_', ' ')} (published {indication} direction)",
                    direction="support_range",
                    goal_kind="percentile_range" if kind == "species" else "graded_capacity",
                    origin="indication_hypothesis",
                    target_status="conditional",
                    recipient_observation_ids=[obs.observation_id] if obs else [],
                    reference_id=obs.reference_bundle_id if obs else None,
                    desired_state={
                        "kind": "percentile_range" if kind == "species" else "graded_supply",
                        "lower": LOW_PERCENTILE if kind == "species" else (obs.reference_p25 if obs else None),
                        "upper": 100.0 if kind == "species" else (obs.reference_p75 if obs else None),
                        "unit": "reference_percentile" if kind == "species" else (obs.unit if obs else None),
                        "statement": "toward the reference range for a feature the literature reports as depleted in this condition",
                    },
                    clinical_endpoint=None,
                    mechanism_source_ids=[
                        f"profile:{indication}/{mod_name}",
                        *[str(s) for s in (module.get("sources") or [])],
                    ],
                    priority_class="moderate",
                    counting_group=counting,
                    uncertainty_reasons=[
                        "indication_hypothesis_not_a_validated_treatment_endpoint",
                        "pattern_weight_is_not_a_therapeutic_weight",
                    ],
                    rule_id="T-INDICATION-1",
                    recipient_value=obs.value if obs else feat.get("raw_value"),
                    recipient_percentile=ref_pct,
                    recipient_status=obs.status if obs else "quantified",
                    notes=[f"module {module.get('label') or mod_name}"],
                )
            )


def _guild_counting_for_species(recipient: MaterialData, species: str) -> str:
    obs = recipient.observation("species", species)
    guilds = list((obs.extra.get("groups") if obs else None) or [])
    if guilds:
        return _counting_group_for_guild(str(guilds[0]))
    return f"species:{species}"


def _context(obs: Observation, counting: str, label: str, why: str) -> TargetRecord:
    return TargetRecord(
        target_id=content_id("target", {"k": "context", "kind": obs.feature_kind, "f": obs.source_id}),
        feature_group_id=counting,
        feature_kind=obs.feature_kind,
        feature_ids=[obs.source_id],
        label=label,
        direction="monitor",
        goal_kind="context",
        origin="report_candidate",
        target_status="context_only",
        recipient_observation_ids=[obs.observation_id],
        reference_id=obs.reference_bundle_id,
        desired_state={"kind": "context", "statement": why},
        clinical_endpoint=None,
        mechanism_source_ids=[],
        priority_class="context",
        counting_group=counting,
        uncertainty_reasons=["no_sourced_direction" if "polarity" in why else "already_at_or_above_reference"],
        rule_id="T-CONTEXT-1",
        recipient_value=obs.value,
        recipient_percentile=obs.percentile,
        recipient_status=obs.status,
        notes=[why],
    )


def _context_targets(recipient: MaterialData, push: Any) -> None:
    """Everything shown but deliberately outside the benefit denominator."""
    for obs in recipient.observations:
        if obs.feature_kind == "ecology" and obs.source_id in ("shannon_index", "species_richness", "gmwi2"):
            push(
                _context(
                    obs,
                    f"ecology:{obs.source_id}",
                    f"{obs.source_id.replace('_', ' ')} (descriptor)",
                    "community descriptor with no universal polarity; never a matching objective "
                    "(rule C01)",
                )
            )
        elif obs.feature_kind == "guild" and str(obs.extra.get("higher_means")) != "favourable":
            push(
                _context(
                    obs,
                    _counting_group_for_guild(obs.source_id),
                    f"{obs.extra.get('label') or obs.source_id} (uncurated or adverse polarity)",
                    "guild polarity is not favourable, so neither supply nor avoidance is a sourced "
                    "goal here",
                )
            )
        elif obs.feature_kind == "profile" and (obs.percentile or 0) >= HIGH_PERCENTILE:
            push(
                _context(
                    obs,
                    f"pattern:{obs.source_id}",
                    f"{obs.extra.get('label') or obs.source_id} resemblance ({obs.percentile:.0f}th)",
                    "disease-pattern resemblance is context: it is not a transmissible entity and "
                    "cannot exclude a donor (rule C01)",
                )
            )


def _custom_goal(recipient: MaterialData, goal: dict[str, Any], push: Any) -> None:
    """An explicit research goal: allowed without trial evidence, clearly labelled."""
    kind = str(goal.get("feature_kind") or "species")
    name = str(goal.get("feature_id") or "")
    if not name:
        return
    obs = recipient.observation(kind, name)
    direction = str(goal.get("direction") or "supply")
    counting = (
        _counting_group_for_panel(name)
        if kind == "gene_panel"
        else _guild_counting_for_species(recipient, name)
    )
    push(
        TargetRecord(
            target_id=content_id("target", {"k": "user", "f": name, "kind": kind}),
            feature_group_id=counting,
            feature_kind=kind,
            feature_ids=[name],
            label=str(goal.get("label") or f"{name.replace('_', ' ')} (user goal)"),
            direction=direction,  # type: ignore[arg-type]
            goal_kind=str(goal.get("goal_kind") or ("presence" if kind == "species" else "graded_capacity")),
            origin="explicit_user_goal",
            target_status="active",
            recipient_observation_ids=[obs.observation_id] if obs else [],
            reference_id=obs.reference_bundle_id if obs else None,
            desired_state=dict(goal.get("desired_state") or {"kind": "presence", "statement": "requested by the user"}),
            clinical_endpoint=None,
            mechanism_source_ids=[str(s) for s in (goal.get("source_ids") or [])],
            priority_class=str(goal.get("priority_class") or "user"),
            counting_group=counting,
            uncertainty_reasons=["explicit_user_goal_without_clinical_efficacy_evidence"],
            rule_id="T-USER-1",
            recipient_value=obs.value if obs else None,
            recipient_percentile=obs.percentile if obs else None,
            recipient_status=obs.status if obs else "unassessed",
        )
    )


def goal_weights(targets: list[TargetRecord]) -> tuple[dict[str, float], dict[str, float]]:
    """Equal group weights, then equal weights inside each group (§7.2).

    Transparent engineering defaults, not published efficacy coefficients.
    """
    groups: dict[str, list[TargetRecord]] = {}
    for t in targets:
        if t.is_positive_scoring:
            groups.setdefault(t.counting_group, []).append(t)
    if not groups:
        return {}, {}
    alpha = {g: 1.0 / len(groups) for g in sorted(groups)}
    beta: dict[str, float] = {}
    for members in groups.values():
        for t in members:
            beta[t.target_id] = 1.0 / len(members)
    return alpha, beta
