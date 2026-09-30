"""Restricted safety-rule DSL (spec v04 §8.1).

Rules are declarative YAML, never Python. A rule names the intervention
classes it applies to, a severity, a boolean condition over typed context
fields, and — mandatorily — what to do when a referenced field is unknown.
The compiler rejects undeclared field paths, unknown operators, invalid units
and nullable fields lacking ``on_unknown``; the evaluator fails closed.

Three propositions are kept apart everywhere downstream (§8.1):

* ``known_safety_signal`` — what is medically known about the intervention
  in this context;
* ``display_policy`` — how it may be rendered;
* ``action_status`` — whether it can be discussed as an action at all.

Precedence is deterministic: labeled contraindication > missing-safety-data
suppression > clinician review > caution > none known. No positive study
overrides a higher safety state.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.errors import OpenBiotaError

# --------------------------------------------------------------------------- #
# the context schema — every field a rule may reference
# --------------------------------------------------------------------------- #

#: field path -> (type, unit or None, ontology/note, nullable)
CONTEXT_SCHEMA: Final[dict[str, tuple[str, str | None, str, bool]]] = {
    "subject.age_years": ("number", "years", "years since birth; None when unknown", True),
    "subject.age_known": ("bool", None, "age assurance present", False),
    "subject.is_minor": ("bool", None, "age_years < 18; None when age unknown", True),
    "subject.adult_status": (
        "enum:minor,adult,assumed_adult,unknown", None,
        "age posture: from a supplied age, assumed adult in an assuming mode, or unknown (fails closed)", False,
    ),
    "subject.sex": ("enum:male,female", None, "sex at birth", True),
    "subject.mode": ("enum:participant,clinician,research", None, "report mode", False),
    "context.pregnant": ("bool", None, "HP:0001197-like context flag", True),
    "context.breastfeeding": ("bool", None, "", True),
    "context.immunocompromise.status": ("bool", None, "any immunocompromising condition or therapy", True),
    "context.immunocompromise.severe": ("bool", None, "severe (e.g. neutropenia, transplant)", True),
    "context.central_venous_catheter": ("bool", None, "", True),
    "context.critical_illness": ("bool", None, "ICU-level illness", True),
    "context.recent_hospitalization": ("bool", None, "within 90 days", True),
    "context.preterm_infant": ("bool", None, "", True),
    "context.short_bowel_or_structural_gi": ("bool", None, "short bowel, stricture, obstruction history", True),
    "context.active_ibd_flare": ("bool", None, "", True),
    "context.eating_disorder_or_malnutrition": ("bool", None, "", True),
    "context.liver_disease": ("bool", None, "", True),
    "context.kidney_disease": ("bool", None, "any CKD stage or dialysis", True),
    "context.egfr": ("number", "mL/min/1.73m2", "", True),
    "context.gi_bleeding_history": ("bool", None, "", True),
    "context.surgery_planned_within_14_days": ("bool", None, "", True),
    "context.bleeding_disorder": ("bool", None, "", True),
    "context.kidney_stone_history": ("bool", None, "", True),
    "context.diabetes_on_glucose_lowering": ("bool", None, "", True),
    "context.medications": ("list:medication_class", None, "RxNorm-class labels, see MEDICATION_CLASSES", False),
    "context.conditions": ("list:condition", None, "clinician-confirmed conditions", False),
    "context.self_reported_conditions": ("list:condition", None, "self-reported, unconfirmed", False),
    "context.allergies": ("list:allergen", None, "", False),
    "context.symptoms": ("list:symptom", None, "", False),
    "context.external_labs": ("list:lab", None, "external lab results supplied", False),
    "context.breath_methane_ppm": ("number", "ppm", "standardized breath test", True),
}

MEDICATION_CLASSES: Final = frozenset({
    "anticoagulant", "antiplatelet", "glucose_lowering", "insulin", "metformin", "levodopa",
    "mao_inhibitor", "immunosuppressant", "chemotherapy", "irinotecan", "warfarin", "ppi",
    "antibiotic", "ssri", "snri", "biologic", "jak_inhibitor", "corticosteroid", "laxative",
    "rifaximin", "lactulose", "statin", "antihypertensive", "salicylate", "nsaid",
    "thyroid_hormone", "lithium", "opioid", "gluten_free_diet_started",
})

OPERATORS: Final = frozenset({
    "eq", "ne", "lt", "lte", "gt", "gte", "in", "not_in", "exists",
    "all", "any", "none", "medication_class_present", "condition_present", "allergy_present",
})

SEVERITIES: Final = ("labeled_contraindication", "suppress", "clinician_only", "caution")
#: known_safety_signal enum, in precedence order (highest first)
SIGNALS: Final = ("labeled_contraindication", "unknown_incomplete_context", "clinician_review", "caution", "none_known")
SIGNAL_RANK: Final = {s: i for i, s in enumerate(SIGNALS)}

INTERVENTION_CLASSES: Final = frozenset({
    "probiotic", "live_biotherapeutic", "prebiotic", "fermented_food", "botanical_supplement",
    "vitamin_supplement", "diet", "prescription_drug", "fmt_microbiota_product", "multi_ingredient_product",
    "lifestyle", "clinical_route", "otc_drug",
    # Live organisms applied to skin. A separate class from the oral ones on
    # purpose: route changes both what the evidence supports and which safety
    # rules apply, and a topical study must never read as an oral regimen.
    "topical_live_biotherapeutic",
})


class RuleCompileError(OpenBiotaError):
    pass


@dataclass(frozen=True, slots=True)
class Condition:
    op: str
    field: str | None = None
    value: Any = None
    children: tuple[Condition, ...] = ()

    def referenced_fields(self) -> set[str]:
        out = {self.field} if self.field else set()
        for c in self.children:
            out |= c.referenced_fields()
        return out


@dataclass(frozen=True, slots=True)
class SafetyRule:
    rule_id: str
    intervention_classes: frozenset[str]
    severity: str
    when: Condition
    on_unknown_fields: frozenset[str]
    on_unknown_action: str
    suppress_self_directed_option: bool
    display_safety_reason: bool
    reason: str
    source_refs: tuple[str, ...]
    participant_text: str | None = None

    def applies_to(self, intervention_class: str) -> bool:
        return intervention_class in self.intervention_classes


@dataclass(frozen=True, slots=True)
class RuleOutcome:
    rule_id: str
    fired: bool
    unknown_fields: tuple[str, ...]
    signal: str                 # one of SIGNALS
    reason: str
    participant_text: str | None
    source_refs: tuple[str, ...]


# --------------------------------------------------------------------------- #
# compiler
# --------------------------------------------------------------------------- #


def _compile_condition(node: Any, rule_id: str) -> Condition:
    if not isinstance(node, Mapping):
        raise RuleCompileError(f"{rule_id}: condition must be a mapping, got {type(node).__name__}")
    if len(node) == 1 and next(iter(node)) in ("all", "any", "none"):
        op = next(iter(node))
        kids = node[op]
        if not isinstance(kids, list) or not kids:
            raise RuleCompileError(f"{rule_id}: '{op}' needs a non-empty list")
        return Condition(op=op, children=tuple(_compile_condition(k, rule_id) for k in kids))
    op = node.get("op")
    if op not in OPERATORS or op in ("all", "any", "none"):
        raise RuleCompileError(f"{rule_id}: unknown operator {op!r}")
    field_path = node.get("field")
    if op in ("medication_class_present", "condition_present", "allergy_present"):
        field_path = field_path or {
            "medication_class_present": "context.medications",
            "condition_present": "context.conditions",
            "allergy_present": "context.allergies",
        }[op]
        if op == "medication_class_present" and node.get("value") not in MEDICATION_CLASSES:
            raise RuleCompileError(f"{rule_id}: undeclared medication class {node.get('value')!r}")
    if field_path not in CONTEXT_SCHEMA:
        raise RuleCompileError(f"{rule_id}: undeclared field path {field_path!r}")
    ftype = CONTEXT_SCHEMA[field_path][0]
    value = node.get("value")
    if op in ("lt", "lte", "gt", "gte") and (ftype != "number" or not isinstance(value, (int, float))):
        raise RuleCompileError(f"{rule_id}: {op} needs a numeric field and value ({field_path})")
    if op in ("eq", "ne") and ftype == "bool" and not isinstance(value, bool):
        raise RuleCompileError(f"{rule_id}: {field_path} is boolean; value must be true/false")
    if op in ("eq", "ne") and ftype.startswith("enum:") and value not in ftype[5:].split(","):
        raise RuleCompileError(f"{rule_id}: {value!r} is not in the enum for {field_path}")
    if op in ("in", "not_in") and not isinstance(value, list):
        raise RuleCompileError(f"{rule_id}: {op} needs a list value")
    if "unit" in node and node["unit"] != CONTEXT_SCHEMA[field_path][1]:
        raise RuleCompileError(f"{rule_id}: unit {node['unit']!r} does not match {field_path}")
    return Condition(op=op, field=field_path, value=value)


def compile_rule(raw: Mapping[str, Any]) -> SafetyRule:
    rule_id = str(raw.get("rule_id") or "")
    if not rule_id:
        raise RuleCompileError("rule without rule_id")
    classes = frozenset(raw.get("applies_to", {}).get("intervention_classes", []))
    bad = classes - INTERVENTION_CLASSES
    if bad or not classes:
        raise RuleCompileError(f"{rule_id}: unknown intervention classes {sorted(bad)}" if bad else f"{rule_id}: applies_to is empty")
    severity = raw.get("severity")
    if severity not in SEVERITIES:
        raise RuleCompileError(f"{rule_id}: severity must be one of {SEVERITIES}")
    when = _compile_condition(raw.get("when"), rule_id)
    referenced = when.referenced_fields()
    nullable = {f for f in referenced if CONTEXT_SCHEMA[f][3]}
    on_unknown = raw.get("on_unknown") or {}
    declared = frozenset(on_unknown.get("fields", []))
    if nullable and not on_unknown:
        raise RuleCompileError(f"{rule_id}: references nullable fields {sorted(nullable)} without on_unknown")
    undeclared = nullable - declared
    if undeclared:
        raise RuleCompileError(f"{rule_id}: nullable fields {sorted(undeclared)} lack on_unknown coverage")
    unknown_action = on_unknown.get("action", "suppress") if nullable else "none"
    if nullable and unknown_action not in ("suppress", "clinician_only", "caution"):
        raise RuleCompileError(f"{rule_id}: on_unknown.action must fail closed (suppress|clinician_only|caution)")
    action = raw.get("action") or {}
    return SafetyRule(
        rule_id=rule_id, intervention_classes=classes, severity=severity, when=when,
        on_unknown_fields=declared, on_unknown_action=unknown_action,
        suppress_self_directed_option=bool(action.get("suppress_self_directed_option", True)),
        display_safety_reason=bool(action.get("display_safety_reason", True)),
        reason=str(raw.get("reason") or ""), source_refs=tuple(str(s) for s in raw.get("source_refs", [])),
        participant_text=raw.get("participant_text"),
    )


def load_rules(path: Path) -> list[SafetyRule]:
    data = yaml.safe_load(path.read_text()) or {}
    rules = [compile_rule(r) for r in data.get("safety_rules", [])]
    ids = [r.rule_id for r in rules]
    if len(set(ids)) != len(ids):
        raise RuleCompileError("duplicate rule_id in safety rules")
    return rules


# --------------------------------------------------------------------------- #
# evaluation
# --------------------------------------------------------------------------- #


_UNKNOWN = object()


def lookup(context: Mapping[str, Any], path: str) -> Any:
    """Dotted lookup; missing or None -> _UNKNOWN sentinel."""
    node: Any = context
    for part in path.split("."):
        if not isinstance(node, Mapping) or part not in node:
            return _UNKNOWN
        node = node[part]
    return _UNKNOWN if node is None else node


def _eval(cond: Condition, context: Mapping[str, Any]) -> tuple[bool | None, set[str]]:
    """Returns (truth or None if unknown, unknown fields encountered)."""
    if cond.op in ("all", "any", "none"):
        results = [_eval(c, context) for c in cond.children]
        unknown = set().union(*(u for _, u in results))
        truths = [t for t, _ in results]
        if cond.op == "all":
            if any(t is False for t in truths):
                return False, unknown
            return (None if any(t is None for t in truths) else True), unknown
        if cond.op == "any":
            if any(t is True for t in truths):
                return True, unknown
            return (None if any(t is None for t in truths) else False), unknown
        if any(t is True for t in truths):
            return False, unknown
        return (None if any(t is None for t in truths) else True), unknown

    assert cond.field is not None
    value = lookup(context, cond.field)
    if cond.op == "exists":
        return value is not _UNKNOWN, set()
    if value is _UNKNOWN:
        return None, {cond.field}
    if cond.op in ("medication_class_present", "condition_present", "allergy_present", "in", "not_in"):
        if cond.op in ("in", "not_in"):
            hit = value in cond.value
            return (hit if cond.op == "in" else not hit), set()
        items = value if isinstance(value, (list, tuple, set)) else [value]
        return cond.value in items, set()
    try:
        if cond.op == "eq":
            return value == cond.value, set()
        if cond.op == "ne":
            return value != cond.value, set()
        if cond.op == "lt":
            return float(value) < float(cond.value), set()
        if cond.op == "lte":
            return float(value) <= float(cond.value), set()
        if cond.op == "gt":
            return float(value) > float(cond.value), set()
        if cond.op == "gte":
            return float(value) >= float(cond.value), set()
    except (TypeError, ValueError):
        return None, {cond.field}
    return None, {cond.field}


def _signal_for(severity_or_action: str) -> str:
    return {
        "labeled_contraindication": "labeled_contraindication",
        "suppress": "unknown_incomplete_context",
        "clinician_only": "clinician_review",
        "caution": "caution",
    }[severity_or_action]


def evaluate_rule(rule: SafetyRule, context: Mapping[str, Any]) -> RuleOutcome:
    truth, unknown = _eval(rule.when, context)
    if truth is True:
        return RuleOutcome(rule.rule_id, True, tuple(sorted(unknown)), _signal_for(rule.severity),
                           rule.reason, rule.participant_text, rule.source_refs)
    if truth is None:
        # fail closed with the declared unknown action
        return RuleOutcome(rule.rule_id, False, tuple(sorted(unknown)), _signal_for(rule.on_unknown_action),
                           f"{rule.reason} (not on file: {', '.join(sorted(field_label(u) for u in unknown))})",
                           rule.participant_text, rule.source_refs)
    return RuleOutcome(rule.rule_id, False, (), "none_known", "", None, rule.source_refs)


@dataclass(slots=True)
class SafetyVerdict:
    known_safety_signal: str
    display_policy: str        # participant_summary | clinician_only | evidence_only | suppressed
    action_status: str         # discussable_after_confirmation | evidence_only | suppressed
    outcomes: list[RuleOutcome] = field(default_factory=list)

    @property
    def reasons(self) -> list[str]:
        return [o.reason for o in self.outcomes if o.signal != "none_known" and o.reason]

    def to_json(self) -> dict[str, Any]:
        return {
            "known_safety_signal": self.known_safety_signal,
            "display_policy": self.display_policy,
            "action_status": self.action_status,
            "rules": [
                {"rule_id": o.rule_id, "fired": o.fired, "unknown_fields": list(o.unknown_fields),
                 "signal": o.signal, "reason": o.reason, "source_refs": list(o.source_refs)}
                for o in self.outcomes if o.signal != "none_known"
            ],
        }


def apply_rules(
    rules: Iterable[SafetyRule], intervention_class: str, context: Mapping[str, Any],
    *, base_display: str, base_action: str,
) -> SafetyVerdict:
    """Evaluate every applicable rule; combine by precedence; never let a
    positive study lift a higher safety state."""
    outcomes = [evaluate_rule(r, context) for r in rules if r.applies_to(intervention_class)]
    signal = "none_known"
    for o in outcomes:
        if SIGNAL_RANK[o.signal] < SIGNAL_RANK[signal]:
            signal = o.signal
    display, action = base_display, base_action
    if signal == "labeled_contraindication":
        display, action = "suppressed", "suppressed"
    elif signal == "unknown_incomplete_context":
        display = "evidence_only" if base_display != "suppressed" else "suppressed"
        action = "suppressed"
    elif signal == "clinician_review":
        display = "clinician_only" if base_display in ("participant_summary",) else base_display
        action = "evidence_only" if base_action == "discussable_after_confirmation" else base_action
    return SafetyVerdict(signal, display, action, outcomes)


FIELD_LABELS: Final[dict[str, str]] = {
    "subject.age_years": "age", "subject.age_known": "age", "subject.is_minor": "age",
    "subject.sex": "sex", "context.pregnant": "pregnancy", "context.breastfeeding": "breastfeeding",
    "context.immunocompromise.status": "immune status", "context.immunocompromise.severe": "immune status",
    "context.central_venous_catheter": "central line", "context.critical_illness": "critical illness",
    "context.recent_hospitalization": "recent hospital stay", "context.preterm_infant": "prematurity",
    "context.short_bowel_or_structural_gi": "bowel structure (short bowel, stricture, obstruction)",
    "context.active_ibd_flare": "active IBD flare", "context.eating_disorder_or_malnutrition": "eating disorder or malnutrition",
    "context.liver_disease": "liver disease", "context.kidney_disease": "kidney disease", "context.egfr": "kidney function (eGFR)",
    "context.gi_bleeding_history": "GI bleeding history", "context.surgery_planned_within_14_days": "planned surgery",
    "context.bleeding_disorder": "bleeding disorder", "context.kidney_stone_history": "kidney stone history",
    "context.diabetes_on_glucose_lowering": "glucose-lowering treatment", "context.medications": "medication list",
    "context.conditions": "confirmed conditions", "context.self_reported_conditions": "self-reported conditions",
    "context.allergies": "allergies", "context.symptoms": "symptoms", "context.external_labs": "external lab results",
    "context.breath_methane_ppm": "breath methane",
}


def field_label(field: str) -> str:
    return FIELD_LABELS.get(field, field.split(".")[-1].replace("_", " "))


def participant_safety_sentence(verdict: SafetyVerdict) -> str:
    """Never 'safe' or 'contraindicated' from incomplete metadata (§8.1)."""
    if verdict.known_safety_signal == "none_known":
        return "No known safety signal is on record for this context; that is not a statement that it is safe for you."
    if verdict.known_safety_signal == "unknown_incomplete_context":
        fields = sorted({field_label(f) for o in verdict.outcomes for f in o.unknown_fields})
        return (
            "The safety context needed to discuss this is not on file"
            + (f" — {', '.join(fields)} — " if fields else ", ")
            + "so it is shown as evidence only. Supplying a sample manifest with these fields unlocks the full card."
        )
    if verdict.known_safety_signal == "labeled_contraindication":
        return "A labeled contraindication applies in this context; this is not shown as an option."
    if verdict.known_safety_signal == "clinician_review":
        return "A known safety consideration applies; review with a clinician before any change."
    return "A caution is on record for this context."
