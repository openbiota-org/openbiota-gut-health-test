"""Findings → triggers → evidence cards: the glue between the engines and the report.

This is the one place where the run's readings (functional panels, taxon
groups, disease-pattern profiles, organism findings) are turned into the
typed triggers of spec v04 §7.7 and matched against the intervention
registry. The result, an :class:`ActionPlan`, is what the PDF renders as
"what the human research says about this reading" under every detail
section and as the evidence overview in the front atlas.

Nothing here scores an intervention. Every card is the registry's evidence
(supporting and against), the studied population and protocol, the safety
verdict and the fixed participant sentences. If a reading has no mapping,
the plan says so explicitly — that is a result, not a gap.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from openbiota.coherence import CoherenceReport, check_coherence
from openbiota.context import ContextLedger, resolve_context
from openbiota.detection import DetectionModel
from openbiota.findings import (
    BUCKETS,
    FindingsResult,
    TaxonFinding,
    build_findings,
    find_contradictions,
    reference_context,
)
from openbiota.interpretations import load_interpretation_registry
from openbiota.interventions import (
    EvidenceSummary,
    Registry,
    TriggerRecord,
    default_registry_dir,
    generate,
    load_registry,
    materialize_triggers,
)
from openbiota.logging_util import Reporter
from openbiota.taxongroups import load_taxon_group_set

LEVEL_FOR_BAND: dict[int, str] = {
    -3: "notably_low", -2: "low", -1: "somewhat_low", 0: "typical",
    1: "somewhat_high", 2: "high", 3: "notably_high",
}

#: Findings bucket → trigger bucket string (drives selector direction).
BUCKET_TRIGGER: dict[int, str] = {
    2: "high", 3: "low", 4: "absent_expected", 5: "unexpected_present", 6: "concern",
}


def band_of(percentile: float | None) -> int | None:
    if percentile is None:
        return None
    if percentile >= 95:
        return 3
    if percentile >= 85:
        return 2
    if percentile >= 75:
        return 1
    if percentile <= 5:
        return -3
    if percentile <= 15:
        return -2
    if percentile <= 25:
        return -1
    return 0


def level_of(percentile: float | None) -> str | None:
    b = band_of(percentile)
    return None if b is None else LEVEL_FOR_BAND[b]


@dataclass(slots=True)
class ActionPlan:
    """Everything the report needs from the intervention layer."""

    findings: FindingsResult | None
    triggers: list[TriggerRecord]
    summaries: list[EvidenceSummary]
    registry_lock: dict[str, Any]
    mode: str
    notes: list[str] = field(default_factory=list)
    #: Every function measured both by gene panel and by carrier group, graded
    #: for agreement. Computed whenever both engines ran, whether or not the
    #: organism findings succeeded, because the two percentile maps come from
    #: the engines directly.
    coherence: CoherenceReport = field(default_factory=CoherenceReport)
    #: What was supplied about the person, what was assumed, what is unknown.
    context_ledger: ContextLedger | None = None
    # indexes ------------------------------------------------------------- #
    by_panel: dict[str, list[EvidenceSummary]] = field(default_factory=dict)
    by_group: dict[str, list[EvidenceSummary]] = field(default_factory=dict)
    by_profile: dict[str, list[EvidenceSummary]] = field(default_factory=dict)
    by_species: dict[str, list[EvidenceSummary]] = field(default_factory=dict)

    def cards_for(self, kind: str, target: str) -> list[EvidenceSummary]:
        table = {"panel": self.by_panel, "group": self.by_group,
                 "profile": self.by_profile, "species": self.by_species}[kind]
        return table.get(target, [])

    @property
    def n_with_human_trials(self) -> int:
        return sum(1 for s in self.summaries if s.supporting_evidence_lane in ("A", "B", "C"))

    @property
    def n_no_supported(self) -> int:
        return sum(1 for s in self.summaries if s.render == "no_supported_targeted_intervention")

    @property
    def n_clinical_routes(self) -> int:
        return sum(1 for s in self.summaries if s.render == "guideline_route_after_confirmation")

    def to_json(self) -> dict[str, Any]:
        return {
            "registry": self.registry_lock,
            "mode": self.mode,
            "n_triggers": len(self.triggers),
            "n_evidence_summaries": len(self.summaries),
            "triggers": [t.to_json() for t in self.triggers],
            "evidence_summaries": [s.to_json() for s in self.summaries],
            "findings": None if self.findings is None else self.findings.to_json(include_unremarkable=False),
            "notes": list(self.notes),
            "what_this_is": (
                "For every reading in this report, the human research on record about changing it: "
                "supporting and against evidence side by side, the exact studied product or diet, the "
                "population it was studied in, and the safety rules that apply. It is never a treatment "
                "score, a probability of benefit or a personalised prescription."
            ),
        }


def _trigger_index(summaries: Sequence[EvidenceSummary], triggers: Sequence[TriggerRecord]) -> dict[str, dict[str, list[EvidenceSummary]]]:
    by_id = {t.trigger_id: t for t in triggers}
    out: dict[str, dict[str, list[EvidenceSummary]]] = {
        "panel": defaultdict(list), "group": defaultdict(list), "profile": defaultdict(list), "species": defaultdict(list),
    }
    kind_of = {
        "functional_capacity_finding": "panel", "ecological_finding": "group",
        "research_profile_resemblance": "profile", "exact_microbiome_feature": "species",
    }
    for s in summaries:
        seen: set[tuple[str, str]] = set()
        for tid in s.trigger_ids:
            t = by_id.get(tid)
            if t is None:
                continue
            kind = kind_of.get(t.trigger_type)
            if kind is None or (kind, t.target) in seen:
                continue
            seen.add((kind, t.target))
            out[kind][t.target].append(s)
    return out


def build_action_plan(
    *,
    sample: str,
    rows: Sequence[Any],                       # MetaboliteRow (pdfreport)
    similarity: Any | None,                    # SimilarityStage
    gates: Any | None,                         # QCGates
    subject_age: float | None,
    subject_sex: str | None,
    subject_country: str | None,
    medications: Sequence[str] | None,
    mode: str,
    taxa_dir: Path,
    registry_dir: Path | None = None,
    reporter: Reporter | None = None,
    out_of_scope_fractions: Mapping[str, float] | None = None,
    context_ledger: ContextLedger | None = None,
) -> ActionPlan:
    """Build findings, triggers and evidence cards for one sample."""
    notes: list[str] = []
    registry: Registry = load_registry(registry_dir or default_registry_dir())

    # ---- organism findings (spec §5–6) ------------------------------------ #
    findings: FindingsResult | None = None
    taxa_buckets: dict[str, str] = {}
    if similarity is not None and similarity.taxonomy is not None and similarity.matched_cohort is not None:
        cohort = similarity.matched_cohort
        interp = load_interpretation_registry(taxa_dir / "interpretations")
        ref = reference_context(
            cohort, similarity.references.match_info or {},
            subject_age=subject_age, age_known=subject_age is not None,
        )
        model = DetectionModel.prior(
            profiler=similarity.taxonomy.family, database=similarity.taxonomy.index,
        )
        usable = gates.usable_pairs if gates is not None else None
        classified = max(0.0, 1.0 - similarity.taxonomy.unknown_percent / 100.0)
        qc_pass = bool(gates is not None and gates.overall in ("pass", "warn"))
        mentions: dict[str, list[str]] = defaultdict(list)
        ranked_profiles: list[dict[str, Any]] = []
        for r in similarity.results:
            feats = []
            for feature in r.profile.features_for_engine("metaphlan"):
                if feature.level not in ("species", "strain", "sgb"):
                    continue
                mentions[feature.name].append(r.profile.name)
                feats.append({"species": feature.name, "direction": feature.direction, "weight": float(feature.w)})
            if r.reportable:
                ranked_profiles.append({
                    "name": r.profile.name, "label": getattr(r.profile, "label", "") or r.profile.name,
                    "percentile": r.combined_percentile, "features": feats,
                })
        try:
            findings = build_findings(
                species_percent=similarity.taxonomy.species, cohort=cohort, reference=ref,
                detection_model=model, registry=interp,
                group_set=load_taxon_group_set(taxa_dir),
                usable_reads=usable, classified_fraction=classified, qc_pass=qc_pass,
                profile_mentions=mentions,
            )
            panel_pct = {r.panel: r.percentile for r in rows}
            panel_det = {r.panel: bool(r.detected) for r in rows}
            group_pct = {
                g.group.name: g.percentile for g in (similarity.community.groups if similarity.community else [])
            }
            findings.contradictions = find_contradictions(
                panel_percentiles=panel_pct, panel_detected=panel_det, group_percentiles=group_pct,
                findings=findings, ranked_profiles=ranked_profiles,
            )
        except Exception as exc:  # noqa: BLE001 — findings are a report layer, never fatal
            notes.append(f"organism findings unavailable: {type(exc).__name__}: {exc}")
            if reporter is not None:
                reporter.warn(f"organism findings failed: {type(exc).__name__}: {exc}")
        if findings is not None:
            for t in findings.taxa:
                if t.bucket in BUCKET_TRIGGER and (t.headline or t.bucket in (2, 3)):
                    taxa_buckets[t.species] = BUCKET_TRIGGER[t.bucket]
    else:
        notes.append("No matched reference cohort: organism findings and missing-species calls are unavailable.")

    # ---- double-measured functions ----------------------------------------- #
    coherence = CoherenceReport()
    if similarity is not None and similarity.community is not None:
        coherence = check_coherence(
            panel_percentiles={r.panel: r.percentile for r in rows},
            group_percentiles={g.group.name: g.percentile for g in similarity.community.groups},
            out_of_scope_fractions=out_of_scope_fractions or {},
        )

    # ---- triggers --------------------------------------------------------- #
    panel_levels = {r.panel: level_of(r.percentile) for r in rows if r.detected or r.percentile is not None}
    group_levels: dict[str, str | None] = {}
    if similarity is not None and similarity.community is not None:
        group_levels = {g.group.name: level_of(g.percentile) for g in similarity.community.groups}
    profile_levels: dict[str, str | None] = {}
    abstained: list[str] = []
    if similarity is not None:
        for r in similarity.results:
            if r.reportable:
                profile_levels[r.profile.name] = level_of(r.combined_percentile)
            else:
                abstained.append(r.profile.name)
    # v4.2 skin families: banded from the highest computed panel index in the
    # family. Several panels never fuse into one number; the band is only what
    # decides whether an evidence-only card is shown at all.
    index_levels: dict[str, str | None] = {}
    skin = getattr(similarity, "skin", None) if similarity is not None else None
    if isinstance(skin, Mapping):
        best: dict[str, float] = {}
        for record in skin.values():
            if not isinstance(record, Mapping) or record.get("index") is None:
                continue
            family = str(record.get("profile_family") or "")
            if family:
                best[family] = max(best.get(family, float("-inf")), float(record["index"]))
        index_levels = {family: level_of(value) for family, value in best.items()}
    qc_ok = gates is None or gates.overall != "fail"
    # The safety rules read one nested context. It comes from the same ledger
    # the profile engine used, so a medication assumed absent for scoring is
    # assumed absent for safety too, and both say so in the same words.
    if context_ledger is None:
        context_ledger = resolve_context(
            mode=mode, medications=medications, subject_age=subject_age,
            subject_sex=subject_sex, subject_country=subject_country,
        )
    context: dict[str, Any] = context_ledger.safety_context(
        subject_age=subject_age, subject_sex=subject_sex, subject_country=subject_country,
    )
    triggers = materialize_triggers(
        sample=sample, panel_levels=panel_levels, group_levels=group_levels,
        profile_levels=profile_levels, profile_abstained=abstained,
        taxa_buckets=taxa_buckets, context=context, qc_qualified=qc_ok,
        index_levels=index_levels,
    )
    summaries = generate(registry, triggers, sample=sample, context=context, mode=mode)
    idx = _trigger_index(summaries, triggers)
    plan = ActionPlan(
        findings=findings, triggers=triggers, summaries=summaries, registry_lock=registry.to_lock(),
        mode=mode, notes=notes, coherence=coherence, context_ledger=context_ledger,
        by_panel=dict(idx["panel"]), by_group=dict(idx["group"]),
        by_profile=dict(idx["profile"]), by_species=dict(idx["species"]),
    )
    if reporter is not None:
        reporter.record(
            f"    evidence: {len(triggers)} triggers → {len(summaries)} evidence cards "
            f"({plan.n_with_human_trials} with human trial/guideline evidence, "
            f"{plan.n_clinical_routes} clinical routes, {plan.n_no_supported} 'no supported targeted intervention')"
        )
        if findings is not None:
            c = findings.counts
            reporter.record(
                f"    organisms: {c['bucket_4']} expected-but-not-detected, {c['bucket_3']} relatively low, "
                f"{c['bucket_2']} relatively high, {c['bucket_5']} uncommon-but-present, "
                f"{c['bucket_6']} potential-concern signals"
            )
    return plan


__all__ = ["ActionPlan", "build_action_plan", "level_of", "band_of", "BUCKETS", "TaxonFinding"]
