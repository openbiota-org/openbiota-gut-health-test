"""Matched-reference findings: what is here, what is high or low, what is
commonly there but not here, what is unusual, and what needs qualifying
(spec v04 §5–6).

Three independent axes for every reference-eligible species:

* **detection** — detected / qualified not detected / indeterminate, decided
  by the feature-specific depth-power model in :mod:`openbiota.detection`;
* **reference position** — where the abundance sits among the reference's
  *carriers* (a study-balanced percentile on the frozen fixed-frame
  log-ratio), plus a study-balanced prevalence class
  (core / common / variable / rare) with leave-one-study-out stability;
* **interpretation** — what the literature has associated a higher or lower
  level with, from the curated registry, or an explicit "no curated
  interpretation".

Level and interpretation never collapse: a high level may be favourable,
concerning, neutral or unknown, and only the interpretation axis says which.

Every organism lands in exactly one headline bucket (spec §6.1):

    1  detected in the reference's central range
    2  relatively high versus the reference           (promoted: outside the
    3  relatively low versus the reference             central 90%, BH q ≤ 0.10,
                                                       LOSO-stable, QC pass)
    4  commonly detected in the reference, not detected here (qualified)
    5  uncommon in the reference, detected here
    6  potential-concern signal requiring strain/toxin/clinical qualification
    7  indeterminate or reference unavailable

Multiplicity is controlled once, within a frozen family that records the
evidence type, namespace, feature universe and tail construction; display
filtering never changes the family. A relative call is worded as *relative
expansion/depletion*; "overgrowth" is reserved for an absolute-load assay.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

import numpy as np

from openbiota.detection import (
    DETECTED,
    INDETERMINATE,
    NOT_DETECTED,
    QUALIFIED_POWER,
    DetectionModel,
)
from openbiota.interpretations import InterpretationRegistry, SpeciesInterpretation
from openbiota.refcohort import (
    CLR_FLOOR_PERCENT,
    ReferenceCohort,
    build_taxon_references,
    sample_clr,
)
from openbiota.similarity import rescale_to_classified
from openbiota.taxongroups import TaxonGroupSet

# --------------------------------------------------------------------------- #
# frozen policy
# --------------------------------------------------------------------------- #

FINDING_SCHEMA_VERSION: Final = "4.0.0"
NAMESPACE: Final = "metaphlan_species"
NAMESPACE_VERSION: Final = "mpa_v31_CHOCOPhlAn_201901"
TRANSFORM: Final = "fixed_reference_logratio_v1"

#: A study contributes to a study-balanced estimate only with this many samples.
MIN_STUDY_N: Final = 10
#: ...and to an abundance percentile only with this many carriers.
MIN_POS_PER_STUDY: Final = 5
#: Abundance positioning needs this many carriers from this many studies.
MIN_POSITIVE_N: Final = 50
MIN_INDEPENDENT_STUDIES: Final = 3
#: Prevalence-class certainty and thresholds (spec §5.4).
CLASS_CERTAINTY: Final = 0.95
CORE_AT: Final = 0.80
COMMON_AT: Final = 0.50
RARE_AT: Final = 0.05
#: Fraction of contributing studies that must individually meet the class threshold.
MIN_STUDY_FRACTION: Final = 0.50
#: Headline promotion (spec §5.5).
CENTRAL_LOW: Final = 5.0
CENTRAL_HIGH: Final = 95.0
BH_Q: Final = 0.10
MIN_DIRECTION_STABILITY: Final = 0.80
#: Study-block bootstrap draws.
BOOTSTRAP_DRAWS: Final = 300
#: Detected below this relative abundance (percent) is a trace call.
TRACE_PERCENT: Final = 0.01
#: Strong expected-but-not-detected label needs this much reference.
STRONG_MISSING_MIN_N: Final = 100

BUCKETS: Final = {
    1: "Detected in the typical reference range",
    2: "Relatively high versus the reference",
    3: "Relatively low versus the reference",
    4: "Commonly detected in the reference, not detected here",
    5: "Uncommon in the reference, detected here",
    6: "Potential-concern signal requiring strain/toxin/clinical qualification",
    7: "Indeterminate or reference unavailable",
}

PREVALENCE_CLASSES: Final = ("reference_core", "reference_common", "reference_variable", "reference_rare")

MATCH_EXACT: Final = "exact"
MATCH_FALLBACK: Final = "validated_fallback"
MATCH_UNMATCHED: Final = "unmatched_context_only"
MATCH_UNAVAILABLE: Final = "unavailable"


# --------------------------------------------------------------------------- #
# records
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ReferenceContext:
    """Which comparator was used and how well it fits the subject."""

    reference_model_id: str
    match_status: str
    criteria: dict[str, Any]
    n_samples: int
    n_studies: int
    reason: str
    transform: str = TRANSFORM
    transform_digest: str = ""
    population_transportability: float | None = None
    applicability_reasons: tuple[str, ...] = ()
    pediatric_gate: bool = False

    @property
    def numeric_allowed(self) -> bool:
        """Adult references cannot label a minor; nothing numeric then."""
        return not self.pediatric_gate

    def to_json(self) -> dict[str, Any]:
        return {
            "reference_model_id": self.reference_model_id,
            "match_status": self.match_status,
            "criteria": self.criteria,
            "n_samples": self.n_samples,
            "n_studies": self.n_studies,
            "reason": self.reason,
            "transform": self.transform,
            "transform_digest": self.transform_digest,
            "population_transportability": self.population_transportability,
            "applicability_reasons": list(self.applicability_reasons),
            "pediatric_gate": self.pediatric_gate,
            "namespace": NAMESPACE,
            "namespace_version": NAMESPACE_VERSION,
        }


@dataclass(slots=True)
class TaxonFinding:
    species: str
    display_name: str
    percent: float
    detection_status: str
    detection_power: float | None
    lod95_percent: float | None
    # reference prevalence
    prevalence_class: str | None
    prevalence: float | None
    prevalence_interval: tuple[float, float] | None
    prevalence_pooled: float | None
    prevalence_loso_stable: bool | None
    n_studies_contributing: int
    n_reference_carriers: int
    # reference position (among carriers)
    percentile: float | None
    percentile_interval: tuple[float, float] | None
    tail_p: float | None
    fdr_q: float | None
    direction_stability: float | None
    abundance_comparator_available: bool
    position: str
    level: str
    bucket: int | None
    headline: bool
    label_strength: str  # "strong" | "qualified" | "none"
    in_catalogue: bool
    trace: bool
    interpretation: SpeciesInterpretation
    guilds: tuple[str, ...] = ()
    flags: tuple[str, ...] = ()
    profile_mentions: tuple[str, ...] = ()

    @property
    def detected(self) -> bool:
        return self.detection_status == DETECTED

    @property
    def report_language(self) -> str:
        name = self.display_name
        if self.bucket == 4:
            qual = "the matched reference" if self.label_strength == "strong" else "the reference (context not matched)"
            return (
                f"{name} is commonly detected in {qual} but was not detected in this specimen at the "
                f"stated sequencing depth (detection power {self.detection_power:.0%})."
            )
        if self.bucket == 5:
            return (
                f"{name} was detected here; it is uncommon in the reference "
                f"(carried by {self.prevalence:.0%} of reference samples)."
            )
        if self.bucket == 2:
            return f"{name}: relative expansion versus the reference carriers ({_ordinal(self.percentile)} percentile)."
        if self.bucket == 3:
            return f"{name}: relative depletion versus the reference carriers ({_ordinal(self.percentile)} percentile)."
        if self.bucket == 6:
            return (
                f"{name} was detected. Species-level detection is not a pathogen or infection finding; "
                f"{self.interpretation.qualification_needed or 'strain, toxin or clinical confirmation'} would be needed."
            )
        if self.bucket == 7 and self.detected:
            return f"{name} was detected; a matched comparator was unavailable."
        if self.bucket == 7:
            return f"{name}: absence cannot be judged at this depth (detection power {'unknown' if self.detection_power is None else f'{self.detection_power:.0%}'})."
        if self.bucket == 1:
            return f"{name} was detected within the reference's central range ({_ordinal(self.percentile)} percentile among carriers)."
        return f"{name} was not detected; it is {'' if self.prevalence is None else f'carried by {self.prevalence:.0%} of the reference and '}not expected in every specimen."

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "display_name": self.display_name,
            "percent": round(self.percent, 5),
            "detection": {
                "status": self.detection_status,
                "power_if_present": None if self.detection_power is None else round(self.detection_power, 4),
                "lod95_percent": None if self.lod95_percent is None else round(self.lod95_percent, 5),
                "trace": self.trace,
            },
            "reference_prevalence": {
                "class": self.prevalence_class,
                "study_balanced": None if self.prevalence is None else round(self.prevalence, 4),
                "interval_90": None if self.prevalence_interval is None else [round(v, 4) for v in self.prevalence_interval],
                "pooled": None if self.prevalence_pooled is None else round(self.prevalence_pooled, 4),
                "leave_one_study_out_stable": self.prevalence_loso_stable,
                "n_studies_contributing": self.n_studies_contributing,
                "n_reference_carriers": self.n_reference_carriers,
            },
            "reference_position": {
                "position": self.position,
                "level": self.level,
                "percentile_among_carriers": None if self.percentile is None else round(self.percentile, 1),
                "interval_90": None if self.percentile_interval is None else [round(v, 1) for v in self.percentile_interval],
                "tail_p": None if self.tail_p is None else round(self.tail_p, 5),
                "fdr_q": None if self.fdr_q is None else round(self.fdr_q, 5),
                "leave_one_study_out_direction_stability": None if self.direction_stability is None else round(self.direction_stability, 3),
                "abundance_comparator_available": self.abundance_comparator_available,
                "measurement_basis": "relative_abundance",
                "absolute_load_available": False,
            },
            "bucket": self.bucket,
            "bucket_label": None if self.bucket is None else BUCKETS[self.bucket],
            "headline": self.headline,
            "label_strength": self.label_strength,
            "in_catalogue": self.in_catalogue,
            "interpretation": self.interpretation.to_json(),
            "guilds": list(self.guilds),
            "flags": list(self.flags),
            "profile_mentions": list(self.profile_mentions),
            "report_language": self.report_language,
        }


@dataclass(frozen=True, slots=True)
class Contradiction:
    kind: str
    title: str
    statement: str
    evidence: dict[str, Any]
    consequence: str

    def to_json(self) -> dict[str, Any]:
        return {"kind": self.kind, "title": self.title, "statement": self.statement,
                "evidence": self.evidence, "consequence": self.consequence}


@dataclass(slots=True)
class FindingsResult:
    reference: ReferenceContext
    taxa: list[TaxonFinding]
    contradictions: list[Contradiction]
    detection_model: DetectionModel
    multiplicity_family: dict[str, Any]
    usable_reads: int | None
    qc_pass: bool
    notes: list[str] = field(default_factory=list)

    def bucket(self, n: int) -> list[TaxonFinding]:
        return [t for t in self.taxa if t.bucket == n]

    @property
    def counts(self) -> dict[str, int]:
        out = {f"bucket_{n}": len(self.bucket(n)) for n in BUCKETS}
        out["detected"] = sum(1 for t in self.taxa if t.detected)
        out["not_detected_unremarkable"] = sum(
            1 for t in self.taxa if t.bucket is None and t.detection_status == NOT_DETECTED
        )
        out["indeterminate_unremarkable"] = sum(
            1 for t in self.taxa if t.bucket is None and t.detection_status == INDETERMINATE
        )
        out["reference_eligible"] = len(self.taxa)
        return out

    def to_json(self, *, include_unremarkable: bool = True) -> dict[str, Any]:
        taxa = self.taxa if include_unremarkable else [t for t in self.taxa if t.bucket is not None]
        return {
            "finding_schema_version": FINDING_SCHEMA_VERSION,
            "reference": self.reference.to_json(),
            "detection_model": self.detection_model.to_json(),
            "usable_microbial_reads": self.usable_reads,
            "technical_qc_pass": self.qc_pass,
            "multiplicity_family": self.multiplicity_family,
            "buckets": {str(k): v for k, v in BUCKETS.items()},
            "counts": self.counts,
            "headline_policy": {
                "outside_central_90_percent": True, "bh_q_lte": BH_Q,
                "direction_stability_gte": MIN_DIRECTION_STABILITY, "technical_qc": "pass",
            },
            "taxa": [t.to_json() for t in taxa],
            "contradictions": [c.to_json() for c in self.contradictions],
            "notes": list(self.notes),
            "user_facing_phrase_for_bucket_4": (
                "Commonly detected in the matched reference, but not detected in this specimen "
                "at the stated sequencing depth."
            ),
            "prohibited_language": ["required organism", "deficient species", "must be replaced",
                                    "overgrowth (without absolute load)", "eradication", "infection (from abundance)"],
        }


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _ordinal(p: float | None) -> str:
    if p is None:
        return "—"
    if p >= 99.5:
        return ">99th"
    if p < 0.5:
        return "<1st"
    n = int(round(p))
    suffix = "th" if n % 100 in (11, 12, 13) else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def _bh(pvalues: Sequence[float]) -> list[float]:
    """Benjamini–Hochberg adjusted q-values (monotone)."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        q[i] = min(1.0, running)
    return q


def _level(pct: float | None) -> str:
    if pct is None:
        return "not_positioned"
    if pct < CENTRAL_LOW:
        return "very_low"
    if pct < 25:
        return "low"
    if pct <= 75:
        return "typical"
    if pct <= CENTRAL_HIGH:
        return "high"
    return "very_high"


def transform_digest(cohort: ReferenceCohort) -> str:
    frame = cohort.frame or list(range(cohort.n_taxa))
    payload = json.dumps({
        "transform": TRANSFORM, "floor_percent": CLR_FLOOR_PERCENT,
        "frame_taxa": sorted(cohort.taxa[i] for i in frame), "weights": "equal",
        "namespace_version": NAMESPACE_VERSION,
    }, sort_keys=True)
    return "sha256:" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def reference_context(
    cohort: ReferenceCohort,
    match_info: Mapping[str, Any],
    *,
    subject_age: float | None,
    age_known: bool,
    reference_model_id: str = "OPENBIOTA_HEALTHY_STOOL_MPA3_ADULT_CMD2021_V1",
) -> ReferenceContext:
    """Turn the similarity stage's match record into a typed context."""
    studies = {cohort.metadata[s].study_name for s in cohort.sample_ids if s in cohort.metadata}
    criteria = dict(match_info.get("criteria") or {})
    reasons: list[str] = []
    pediatric = subject_age is not None and subject_age < 18
    if pediatric:
        status = MATCH_UNAVAILABLE
        reasons.append("pediatric_reference_insufficient")
        transport: float | None = None
    elif not criteria:
        status = MATCH_UNMATCHED
        reasons.append("no_match_criteria_supplied" if age_known else "age_sex_country_not_supplied")
        transport = None
    elif match_info.get("matched"):
        status = MATCH_EXACT
        transport = 1.0
    else:
        status = MATCH_FALLBACK
        reasons.append("stratum_too_small_fell_back_to_full_reference")
        transport = 0.7
    return ReferenceContext(
        reference_model_id=reference_model_id,
        match_status=status,
        criteria=criteria,
        n_samples=cohort.n_samples,
        n_studies=len(studies),
        reason=str(match_info.get("reason", "")),
        transform_digest=transform_digest(cohort),
        population_transportability=transport,
        applicability_reasons=tuple(reasons),
        pediatric_gate=pediatric,
    )


# --------------------------------------------------------------------------- #
# the engine
# --------------------------------------------------------------------------- #


def _prevalence_class(p_core: float, p_common: float, p_rare: float) -> str:
    if p_core >= CLASS_CERTAINTY:
        return "reference_core"
    if p_common >= CLASS_CERTAINTY:
        return "reference_common"
    if p_rare >= CLASS_CERTAINTY:
        return "reference_rare"
    return "reference_variable"


def build_findings(
    *,
    species_percent: Mapping[str, float],
    cohort: ReferenceCohort,
    reference: ReferenceContext,
    detection_model: DetectionModel,
    registry: InterpretationRegistry,
    group_set: TaxonGroupSet | None,
    usable_reads: int | None,
    classified_fraction: float | None,
    qc_pass: bool,
    profile_mentions: Mapping[str, Sequence[str]] | None = None,
    seed: int = 20260907,
) -> FindingsResult:
    """Compute every axis for every reference-eligible or detected species."""
    rng = np.random.default_rng(seed)
    species = rescale_to_classified(species_percent)
    notes: list[str] = []

    taxa = list(cohort.taxa)
    index = {t: i for i, t in enumerate(taxa)}
    A = np.asarray(cohort.abundance, dtype=float)  # taxa x samples (percent)
    P = A > 0
    refs = build_taxon_references(cohort)
    Z = np.asarray([refs[t].values for t in taxa], dtype=float)
    z_sample = sample_clr(species, cohort)

    studies = np.asarray([cohort.metadata[s].study_name if s in cohort.metadata else "?" for s in cohort.sample_ids])
    study_names = sorted(set(studies))
    masks = [studies == s for s in study_names]
    study_n = np.asarray([m.sum() for m in masks])
    eligible_studies = [k for k, n in enumerate(study_n) if n >= MIN_STUDY_N]
    if len(eligible_studies) < MIN_INDEPENDENT_STUDIES:
        notes.append(
            f"Only {len(eligible_studies)} reference studies have ≥{MIN_STUDY_N} samples; "
            "study-balanced estimates fall back to pooled ones."
        )

    # ---- prevalence: per study, balanced, bootstrap, LOSO ---------------- #
    prev = np.zeros((len(taxa), len(study_names)))
    for k, m in enumerate(masks):
        prev[:, k] = P[:, m].mean(axis=1) if m.any() else 0.0
    E = np.asarray(eligible_studies, dtype=int) if eligible_studies else np.arange(len(study_names))
    prev_e = prev[:, E]
    balanced_prev = prev_e.mean(axis=1)
    boot_idx = rng.integers(0, len(E), size=(BOOTSTRAP_DRAWS, len(E)))
    boot_prev = prev_e[:, boot_idx].mean(axis=2)  # taxa x B
    p_core = (boot_prev >= CORE_AT).mean(axis=1)
    p_common = (boot_prev >= COMMON_AT).mean(axis=1)
    p_rare = (boot_prev <= RARE_AT).mean(axis=1)
    prev_lo = np.percentile(boot_prev, 5, axis=1)
    prev_hi = np.percentile(boot_prev, 95, axis=1)
    # LOSO means
    total = prev_e.sum(axis=1, keepdims=True)
    loso_prev = (total - prev_e) / max(1, len(E) - 1)  # taxa x studies
    frac_core = (prev_e >= CORE_AT).mean(axis=1)
    frac_common = (prev_e >= COMMON_AT).mean(axis=1)
    pooled_prev = P.mean(axis=1)
    n_carriers = P.sum(axis=1)

    # ---- abundance position among carriers: sample vs each study --------- #
    z_vec = np.asarray([z_sample.get(t, math.log(CLR_FLOOR_PERCENT)) for t in taxa])
    pct_study = np.full((len(taxa), len(study_names)), np.nan)
    for k, m in enumerate(masks):
        if k not in set(E.tolist()):
            continue
        Zs = Z[:, m]
        Ps = P[:, m]
        n_pos = Ps.sum(axis=1)
        below = ((Zs < z_vec[:, None]) & Ps).sum(axis=1)
        ties = ((Zs == z_vec[:, None]) & Ps).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = (below + 0.5 * ties) / n_pos * 100.0
        pct[n_pos < MIN_POS_PER_STUDY] = np.nan
        pct_study[:, k] = pct
    valid = ~np.isnan(pct_study)
    n_valid_studies = valid.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        balanced_pct = np.where(valid, pct_study, 0.0).sum(axis=1) / np.where(n_valid_studies > 0, n_valid_studies, np.nan)
    # bootstrap over valid studies per taxon (vectorised by resampling all
    # studies and re-normalising over valid ones)
    pct_filled = np.where(valid, pct_study, 0.0)
    boot_all = rng.integers(0, len(study_names), size=(BOOTSTRAP_DRAWS, len(study_names)))
    boot_sum = pct_filled[:, boot_all].sum(axis=2)  # taxa x B
    boot_cnt = valid.astype(float)[:, boot_all].sum(axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        boot_pct = boot_sum / boot_cnt
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        pct_lo = np.nanpercentile(np.where(boot_cnt > 0, boot_pct, np.nan), 5, axis=1)
        pct_hi = np.nanpercentile(np.where(boot_cnt > 0, boot_pct, np.nan), 95, axis=1)
    # LOSO percentiles
    sum_valid = pct_filled.sum(axis=1, keepdims=True)
    cnt_valid = valid.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        loso_pct = (sum_valid - pct_filled) / (cnt_valid - valid)  # taxa x studies
    loso_pct = np.where(valid & (cnt_valid - valid > 0), loso_pct, np.nan)

    # ---- detection power for non-detections -------------------------------- #
    n_eff = None
    if usable_reads:
        n_eff = usable_reads * (classified_fraction if classified_fraction else 1.0)
    lod = None if n_eff is None else detection_model.lod95(n_eff)

    # ---- assemble per-taxon records --------------------------------------- #
    catalogue = set(taxa)
    detected_names = {s for s, v in species.items() if v > 0}
    all_names = sorted(catalogue | detected_names)
    guild_of: dict[str, tuple[str, ...]] = {}
    if group_set is not None:
        for g in group_set.groups:
            for m in g.members:
                guild_of.setdefault(m, ())
                guild_of[m] = (*guild_of[m], g.name)
    mentions = profile_mentions or {}

    findings: list[TaxonFinding] = []
    tail_ps: list[tuple[int, float]] = []  # (finding index, p) for the BH family
    for name in all_names:
        value = species.get(name, 0.0)
        detected = value > 0
        interp = registry.disposition(name)
        i = index.get(name)
        flags: list[str] = []
        if i is None:
            # detected but outside the reference catalogue
            findings.append(TaxonFinding(
                species=name, display_name=interp.display_name, percent=value,
                detection_status=DETECTED, detection_power=None, lod95_percent=None,
                prevalence_class=None, prevalence=None, prevalence_interval=None, prevalence_pooled=None,
                prevalence_loso_stable=None, n_studies_contributing=0, n_reference_carriers=0,
                percentile=None, percentile_interval=None, tail_p=None, fdr_q=None,
                direction_stability=None, abundance_comparator_available=False,
                position="comparator_unavailable", level="not_positioned",
                bucket=6 if interp.concern else 7, headline=interp.concern, label_strength="none",
                in_catalogue=False, trace=value <= TRACE_PERCENT, interpretation=interp,
                guilds=guild_of.get(name, ()), flags=("not_in_reference_catalogue",),
                profile_mentions=tuple(mentions.get(name, ())),
            ))
            continue

        pclass = _prevalence_class(p_core[i], p_common[i], p_rare[i])
        # LOSO stability of the class promotion
        if pclass == "reference_core":
            loso_ok = bool(np.all(loso_prev[i] >= CORE_AT)) and bool(frac_core[i] >= MIN_STUDY_FRACTION)
        elif pclass == "reference_common":
            loso_ok = bool(np.all(loso_prev[i] >= COMMON_AT)) and bool(frac_common[i] >= MIN_STUDY_FRACTION)
        elif pclass == "reference_rare":
            loso_ok = bool(np.all(loso_prev[i] <= RARE_AT))
        else:
            loso_ok = True
        if not loso_ok and pclass in ("reference_core", "reference_common"):
            pclass = "reference_variable"  # demoted: promotion did not survive LOSO
            flags.append("class_demoted_leave_one_study_out")
        if not loso_ok and pclass == "reference_rare":
            pclass = "reference_variable"
            flags.append("class_demoted_leave_one_study_out")

        carriers = int(n_carriers[i])
        comparator_ok = bool(carriers >= MIN_POSITIVE_N and n_valid_studies[i] >= MIN_INDEPENDENT_STUDIES)
        power = None
        if n_eff is not None and carriers:
            positives = A[i][P[i]] / 100.0
            power = detection_model.power(n_eff, positives.tolist())

        percentile = None
        interval = None
        tail_p = None
        stability = None
        position = "comparator_unavailable"
        level = "not_positioned"
        if detected and comparator_ok and reference.numeric_allowed:
            percentile = float(balanced_pct[i])
            interval = (float(pct_lo[i]), float(pct_hi[i])) if not math.isnan(pct_lo[i]) else None
            tail_p = max(2.0 * min(percentile, 100.0 - percentile) / 100.0, 1.0 / max(carriers, 1))
            level = _level(percentile)
            loso_vals = loso_pct[i][~np.isnan(loso_pct[i])]
            if len(loso_vals):
                if percentile < CENTRAL_LOW:
                    stability = float((loso_vals < CENTRAL_LOW).mean())
                elif percentile > CENTRAL_HIGH:
                    stability = float((loso_vals > CENTRAL_HIGH).mean())
                else:
                    stability = float(((loso_vals < 50) == (percentile < 50)).mean())
            position = level
        elif detected and not reference.numeric_allowed:
            flags.append("pediatric_reference_insufficient")
        elif detected:
            flags.append("abundance_comparator_unavailable")

        if detected:
            status = DETECTED
        elif not reference.numeric_allowed:
            status = INDETERMINATE
        elif power is not None and power >= QUALIFIED_POWER and qc_pass:
            status = NOT_DETECTED
        else:
            status = INDETERMINATE
            if power is not None and power < QUALIFIED_POWER:
                flags.append("low_depth_for_this_organism")
            if not qc_pass:
                flags.append("technical_qc_not_passed")

        trace = detected and value <= TRACE_PERCENT
        identity_ok = detected and not trace and (lod is None or value / 100.0 >= 2.0 * lod)

        finding = TaxonFinding(
            species=name, display_name=interp.display_name, percent=value,
            detection_status=status, detection_power=power,
            lod95_percent=None if lod is None else lod * 100.0,
            prevalence_class=pclass, prevalence=float(balanced_prev[i]),
            prevalence_interval=(float(prev_lo[i]), float(prev_hi[i])),
            prevalence_pooled=float(pooled_prev[i]), prevalence_loso_stable=bool(loso_ok),
            n_studies_contributing=int(len(E)), n_reference_carriers=carriers,
            percentile=percentile, percentile_interval=interval, tail_p=tail_p, fdr_q=None,
            direction_stability=stability, abundance_comparator_available=comparator_ok,
            position=position, level=level, bucket=None, headline=False, label_strength="none",
            in_catalogue=True, trace=trace, interpretation=interp,
            guilds=guild_of.get(name, ()), flags=tuple(flags),
            profile_mentions=tuple(mentions.get(name, ())),
        )
        findings.append(finding)
        if tail_p is not None:
            tail_ps.append((len(findings) - 1, tail_p))
        # stash identity gate on the record via flags
        if detected and not identity_ok:
            finding.flags = (*finding.flags, "identity_gate_trace_or_near_lod")

    # ---- multiplicity within the frozen family ---------------------------- #
    qs = _bh([p for _, p in tail_ps])
    for (fi, _), q in zip(tail_ps, qs, strict=True):
        findings[fi].fdr_q = q
    universe_digest = hashlib.sha256("|".join(taxa).encode()).hexdigest()[:16]
    family = {
        "multiplicity_family_id": f"MULTIPLICITY_TAX_MPA3_V1:{universe_digest}",
        "evidence_type": "TAX_REL", "namespace": NAMESPACE, "namespace_version": NAMESPACE_VERSION,
        "feature_universe_digest": "sha256:" + universe_digest,
        "n_tests": len(tail_ps), "method": "Benjamini-Hochberg",
        "tail": "two-sided empirical percentile among reference carriers, study-balanced",
        "inclusion": f"detected, ≥{MIN_POSITIVE_N} reference carriers from ≥{MIN_INDEPENDENT_STUDIES} studies",
        "nondetection_policy": "nondetections are classified through TAX_PRES and are not tested here",
    }

    # ---- buckets ----------------------------------------------------------- #
    strong_ok = reference.match_status in (MATCH_EXACT, MATCH_FALLBACK) and reference.n_samples >= STRONG_MISSING_MIN_N and reference.n_studies >= MIN_INDEPENDENT_STUDIES
    for f in findings:
        if not f.in_catalogue:
            continue
        if not reference.numeric_allowed:
            f.bucket = 7 if f.detected else None
            f.position = "comparator_unavailable"
            continue
        if f.detected:
            if f.interpretation.concern:
                f.bucket, f.headline, f.label_strength = 6, True, "qualified"
                f.flags = (*f.flags, "clinical_concern_signal:modality_not_assessed")
            elif f.prevalence_class == "reference_rare" and "identity_gate_trace_or_near_lod" not in f.flags:
                f.bucket, f.headline, f.label_strength = 5, True, ("strong" if strong_ok else "qualified")
                f.position = "uncommon_detected"
            elif f.abundance_comparator_available and f.percentile is not None:
                # A single specimen's percentile is a placement, not a test, so
                # LOW/HIGH are descriptive: the point sits in the outer 5%, the
                # study-block bootstrap keeps it outside the central 80%, and the
                # direction survives dropping any one reference study. The BH
                # q-value across the whole family only decides how strongly the
                # label may be worded (a q that clears BH_Q with a matched
                # reference is 'strong'; otherwise 'qualified').
                lo, hi = f.percentile_interval if f.percentile_interval else (f.percentile, f.percentile)
                extreme = f.percentile < CENTRAL_LOW or f.percentile > CENTRAL_HIGH
                interval_clear = (hi < 10.0) if f.percentile < CENTRAL_LOW else (lo > 90.0)
                promoted = (
                    extreme and interval_clear
                    and f.direction_stability is not None and f.direction_stability >= MIN_DIRECTION_STABILITY
                    and qc_pass
                )
                if promoted:
                    f.bucket = 2 if f.percentile > CENTRAL_HIGH else 3
                    f.headline = True
                    f.label_strength = (
                        "strong" if strong_ok and f.fdr_q is not None and f.fdr_q <= BH_Q else "qualified"
                    )
                    if f.fdr_q is None or f.fdr_q > BH_Q:
                        f.flags = (*f.flags, "descriptive_placement_not_fdr_significant")
                else:
                    f.bucket = 1
                    if extreme:
                        f.flags = (*f.flags, "extreme_but_not_promoted")
            else:
                f.bucket = 7
                f.position = "comparator_unavailable"
        elif f.detection_status == NOT_DETECTED and f.prevalence_class in ("reference_core", "reference_common"):
            f.bucket, f.headline = 4, True
            f.label_strength = "strong" if strong_ok else "qualified"
            f.position = "commonly_detected_not_detected"
        elif f.detection_status == INDETERMINATE and f.prevalence_class in ("reference_core", "reference_common"):
            f.bucket = 7  # would have been expected; absence cannot be judged
        else:
            f.bucket = None  # rare/variable nondetection: the default, counted only

    # Reader-facing prose. These end up printed in the report, so they carry
    # no internal status tokens and no instructions to run anything: a person
    # reading their own results can act on neither, and a bare identifier in
    # the middle of a sentence reads as a defect.
    if reference.match_status == MATCH_UNMATCHED:
        notes.append(
            "No age, sex or country was recorded for this sample, so every comparison here "
            "is against healthy adults in general rather than against people matched to you. "
            "Readings called high, low or missing should be read with that in mind."
        )
    if reference.pediatric_gate:
        notes.append(
            "The reference group is adults and this sample is from someone under 18, so no "
            "high, low or missing labels are applied. Organisms found are listed without a "
            "position, because an adult range would not describe a child."
        )
    if not detection_model.calibrated:
        notes.append(
            "How reliably a low-abundance organism can be detected at this sequencing depth "
            "is estimated from a conservative general curve rather than measured on this "
            "laboratory's own data, so detection limits are approximate."
        )

    return FindingsResult(
        reference=reference, taxa=findings, contradictions=[], detection_model=detection_model,
        multiplicity_family=family, usable_reads=usable_reads, qc_pass=qc_pass, notes=notes,
    )


# --------------------------------------------------------------------------- #
# cross-engine contradictions (spec §15.3)
# --------------------------------------------------------------------------- #

#: Gene panel → curated carrier group that should move with it.
PANEL_GROUP_PAIRS: Final = (
    ("butyrate", "butyrate", "butyrate production capacity", "butyrate-producing taxa"),
    ("h2s", "sulfate_reducers", "hydrogen-sulfide capacity", "sulfate-reducing taxa"),
    ("methane", "methanogens", "methane capacity (mcrA)", "methanogen abundance"),
)

#: Virulence gene → the carrier species the fragments would belong to.
VIRULENCE_CARRIERS: Final = {
    "crc_virulence": ("Fusobacterium_nucleatum", "Bacteroides_fragilis", "Escherichia_coli"),
}

#: Panels whose product is a host exposure nothing in a FASTQ can measure.
UNMEASURED_HOST_EXPOSURE: Final = {
    "cutc": "plasma TMAO", "b12": "serum vitamin B12", "k2": "vitamin K status",
    "folate": "blood folate", "riboflavin": "riboflavin status", "biotin": "biotin status",
    "gaba": "luminal or CNS GABA", "dopamine": "dopamine/tyramine exposure",
    "histamine": "histamine exposure", "tryptamine": "tryptamine exposure",
    "indole": "indole/indoxyl sulfate", "ipa": "indole-3-propionate", "pcresol": "p-cresyl sulfate",
    "imidazole": "imidazole propionate", "bai": "fecal or serum bile acids",
}


def find_contradictions(
    *,
    panel_percentiles: Mapping[str, float | None],
    panel_detected: Mapping[str, bool],
    group_percentiles: Mapping[str, float | None],
    findings: FindingsResult,
    ranked_profiles: Sequence[Mapping[str, Any]],
    symptoms: Sequence[str] = (),
    external_assays: Mapping[str, Any] | None = None,
) -> list[Contradiction]:
    """Search the engines' outputs for disagreements and record each one."""
    out: list[Contradiction] = []
    external = external_assays or {}

    for panel, group, panel_label, group_label in PANEL_GROUP_PAIRS:
        pp, gp = panel_percentiles.get(panel), group_percentiles.get(group)
        if pp is None or gp is None:
            continue
        if pp >= 75 and gp <= 25:
            out.append(Contradiction(
                kind="capacity_vs_carriers", title=f"High {panel_label}, low {group_label}",
                statement=(
                    f"Gene capacity for {panel_label} sits at the {_ordinal(pp)} percentile while the curated "
                    f"{group_label} sit at the {_ordinal(gp)} percentile. The genes are being carried by "
                    "organisms outside the curated set, or the set's members are present below detection."
                ),
                evidence={"panel": panel, "panel_percentile": pp, "group": group, "group_percentile": gp},
                consequence="Neither reading is discarded; action language for this function is blocked and uncertainty is widened.",
            ))
        elif gp >= 75 and pp <= 25:
            out.append(Contradiction(
                kind="carriers_vs_capacity", title=f"High {group_label}, low {panel_label}",
                statement=(
                    f"The curated {group_label} sit at the {_ordinal(gp)} percentile while gene capacity for "
                    f"{panel_label} sits at the {_ordinal(pp)} percentile: the carriers are present but the "
                    "pathway genes were found less than expected — incomplete pathway capacity, or a panel "
                    "that does not cover these organisms' gene variants."
                ),
                evidence={"panel": panel, "panel_percentile": pp, "group": group, "group_percentile": gp},
                consequence="Neither reading is discarded; action language for this function is blocked.",
            ))

    for panel, carriers in VIRULENCE_CARRIERS.items():
        if not panel_detected.get(panel):
            continue
        detected = {t.species for t in findings.taxa if t.detected}
        if not any(c in detected for c in carriers):
            out.append(Contradiction(
                kind="virulence_without_carrier", title="Virulence gene fragments without a resolved carrier",
                statement=(
                    "Fragments matching the CRC-virulence panel were found, but none of the species that carry "
                    "these loci (F. nucleatum, enterotoxigenic B. fragilis, pks+ E. coli) was detected. The "
                    "fragments cannot be attributed to an organism."
                ),
                evidence={"panel": panel, "expected_carriers": list(carriers)},
                consequence="Reported as unattributed fragments; no pathotype or concern label follows from them.",
            ))

    for panel, exposure in UNMEASURED_HOST_EXPOSURE.items():
        if panel in panel_percentiles and panel not in external:
            out.append(Contradiction(
                kind="capacity_without_measurement", title=f"{panel}: capacity measured, {exposure} not",
                statement=(
                    f"The {panel} panel reports genomic capacity only. {exposure[0].upper() + exposure[1:]} was not "
                    "measured in this run, so no level, deficiency or exposure can be stated."
                ),
                evidence={"panel": panel, "not_measured": exposure},
                consequence="Cards for this panel say 'capacity', never production, concentration or deficiency.",
            ))

    if "constipation" in {s.lower() for s in symptoms}:
        mp = panel_percentiles.get("methane")
        if mp is not None and not panel_detected.get("methane", True):
            out.append(Contradiction(
                kind="methane_nondetection_with_constipation",
                title="Methane genes not detected, constipation reported",
                statement="Reported constipation with no mcrA detected in stool: a standardized breath methane test, not stool DNA, is the relevant measurement.",
                evidence={"panel": "methane", "symptom": "constipation"},
                consequence="Stool nondetection is not evidence against intestinal methanogen overgrowth.",
            ))

    by_species = {t.species: t for t in findings.taxa}
    for prof in ranked_profiles:
        pct = prof.get("percentile")
        if pct is None or pct < 75:
            continue
        opposed: list[str] = []
        weight_total = 0.0
        weight_opposed = 0.0
        for feat in prof.get("features", ()):
            sp = feat.get("species")
            direction = feat.get("direction")
            w = float(feat.get("weight", 1.0))
            f = by_species.get(sp)
            if f is None or direction not in ("up", "down"):
                continue
            weight_total += w
            opposes = False
            if direction == "up" and (f.bucket in (3, 4) or (f.percentile is not None and f.percentile < 25)):
                opposes = True
            if direction == "down" and (f.bucket in (2, 5) or (f.percentile is not None and f.percentile > 75)):
                opposes = True
            if opposes:
                weight_opposed += w
                opposed.append(f"{f.display_name} ({'higher' if direction == 'up' else 'lower'} in the signature; {f.level.replace('_', ' ')} here)")
        if weight_total and weight_opposed / weight_total >= 0.30 and opposed:
            out.append(Contradiction(
                kind="profile_direction_opposed",
                title=f"{prof.get('label', prof.get('name'))}: resemblance opposed by matched-reference findings",
                statement=(
                    f"The pattern scores at the {_ordinal(pct)} percentile, yet {len(opposed)} of its taxa point the "
                    f"other way in this specimen ({weight_opposed / weight_total:.0%} of the evidence weight): "
                    + "; ".join(opposed[:6])
                    + (f"; and {len(opposed) - 6} more" if len(opposed) > 6 else "") + "."
                ),
                evidence={"profile": prof.get("name"), "percentile": pct, "opposed_weight_fraction": round(weight_opposed / weight_total, 3), "opposed": opposed},
                consequence="The resemblance is driven by a subset of the signature; it is not read as a coherent match.",
            ))
    return out


__all__ = [
    "BUCKETS",
    "FINDING_SCHEMA_VERSION",
    "Contradiction",
    "FindingsResult",
    "ReferenceContext",
    "TaxonFinding",
    "build_findings",
    "find_contradictions",
    "reference_context",
    "transform_digest",
]
