"""Validation of the profile scoring engine against labelled cohorts.

Three things are measured here, in the order the spec puts them.

**7.2 — CRC end-to-end.** curatedMetagenomicData carries colorectal cancer
cohorts with disease labels. Scoring cases and controls through the identical
engine and computing ROC AUC is the only test of the whole machine against
clinical ground truth. It is the build gate: if a known-good CRC signal cannot
be recovered on labelled data, no number the engine produces for any other
profile means anything.

**7.4 — Cross-profile specificity.** Every profile scored against every
labelled condition. A profile that scores high on everything is measuring
general dysbiosis rather than its own condition, and this matrix is what makes
the dysbiosis anchor interpretable.

**7.5 — Reference stability.** Handled in `refcohort.stability_curve`.

Study-level holdout, never random splitting: pooling cohorts and splitting
randomly leaks recruitment-batch signal and inflates every number.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.logging_util import Reporter
from openbiota.profiles import Profile
from openbiota.refcohort import ReferenceCohort, build_taxon_references, sample_clr
from openbiota.scoring import (
    Z_CLIP,
    FeatureReference,
    summarise_reference,
)

#: Minimum cases and controls for an AUC to be worth reporting.
MIN_PER_ARM: Final = 20


# --------------------------------------------------------------------------- #
# ROC
# --------------------------------------------------------------------------- #


def roc_auc(scores: Sequence[float], labels: Sequence[int]) -> float | None:
    """AUC via the Mann-Whitney U statistic, with tie correction."""
    pairs = [(s, y) for s, y in zip(scores, labels, strict=True) if s is not None]
    positives = [s for s, y in pairs if y == 1]
    negatives = [s for s, y in pairs if y == 0]
    if not positives or not negatives:
        return None

    ordered = sorted(pairs, key=lambda p: p[0])
    ranks: list[float] = [0.0] * len(ordered)
    index = 0
    while index < len(ordered):
        end = index
        while end + 1 < len(ordered) and ordered[end + 1][0] == ordered[index][0]:
            end += 1
        average = (index + end) / 2.0 + 1.0
        for position in range(index, end + 1):
            ranks[position] = average
        index = end + 1

    rank_sum = sum(r for r, (_, y) in zip(ranks, ordered, strict=True) if y == 1)
    n_pos, n_neg = len(positives), len(negatives)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def auc_confidence_interval(
    scores: Sequence[float], labels: Sequence[int], *, replicates: int = 2000, seed: int = 20260904
) -> tuple[float | None, float | None]:
    """Bootstrap percentile interval for the AUC."""
    import random

    rng = random.Random(seed)
    pairs = [(s, y) for s, y in zip(scores, labels, strict=True) if s is not None]
    if len(pairs) < 10:
        return None, None
    values: list[float] = []
    for _ in range(replicates):
        sample = [pairs[rng.randrange(len(pairs))] for _ in range(len(pairs))]
        auc = roc_auc([s for s, _ in sample], [y for _, y in sample])
        if auc is not None:
            values.append(auc)
    if len(values) < 100:
        return None, None
    values.sort()
    return (
        values[int(0.025 * len(values))],
        values[int(0.975 * len(values))],
    )


# --------------------------------------------------------------------------- #
# scoring a labelled cohort
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class LabelledScores:
    """Taxonomic-module scores for a labelled case/control set."""

    profile: str
    condition: str
    scores: list[float] = field(default_factory=list)
    labels: list[int] = field(default_factory=list)
    studies: list[str] = field(default_factory=list)
    sample_ids: list[str] = field(default_factory=list)

    @property
    def n_cases(self) -> int:
        return sum(self.labels)

    @property
    def n_controls(self) -> int:
        return len(self.labels) - self.n_cases


def taxonomic_module_score(
    profile: Profile,
    clr_by_taxon: Mapping[str, float],
    references: Mapping[str, FeatureReference],
) -> float | None:
    """Taxonomic module score for one sample, no side effects.

    Duplicated deliberately from `scoring.score_profile` rather than reused:
    validation must score thousands of samples in a tight loop and does not
    need the full decomposition machinery.
    """
    numerator = 0.0
    denominator = 0.0
    for feature in profile.features(module="taxonomic"):
        reference = references.get(feature.name)
        value = clr_by_taxon.get(feature.name)
        if reference is None or value is None or reference.scale <= 0:
            continue
        z = reference.robust_z(value)
        if z is None:
            continue
        v = feature.d * max(-Z_CLIP, min(Z_CLIP, z)) / Z_CLIP
        numerator += feature.factor * v
        denominator += feature.factor
    return numerator / denominator if denominator > 0 else None


def score_labelled_cohort(
    *,
    profile: Profile,
    case_cohort: ReferenceCohort,
    control_cohort: ReferenceCohort,
    condition: str,
    reporter: Reporter,
) -> LabelledScores:
    """Score cases and controls against a reference built from the controls.

    The reference distribution comes from the **control arm**, never from the
    cases — otherwise the cases help define the distribution they are being
    measured against.

    Using the control arm rather than the global healthy cohort matters more
    than it looks. Oral-origin taxa are three to four times more prevalent in
    one study's controls than across the pooled healthy cohort, because
    prevalence of a rare taxon tracks sequencing depth and population. Scoring
    a study's cases against a pooled reference therefore shifts cases and
    controls up together and destroys the contrast. This is the same
    population-matching requirement the reference cohort applies at query
    time, enforced here.
    """
    taxon_refs = build_taxon_references(control_cohort)
    references = {
        name: summarise_reference(
            name,
            taxon_refs[name].values,
            f"within-study controls (n={control_cohort.n_samples})",
            prevalence=taxon_refs[name].prevalence,
        )
        for name in taxon_refs
    }
    out = LabelledScores(profile=profile.name, condition=condition)

    for cohort, label in ((case_cohort, 1), (control_cohort, 0)):
        for j, sample_id in enumerate(cohort.sample_ids):
            abundances = {
                cohort.taxa[i]: cohort.abundance[i][j]
                for i in range(cohort.n_taxa)
                if cohort.abundance[i][j] > 0
            }
            transformed = sample_clr(abundances, control_cohort)
            score = taxonomic_module_score(profile, transformed, references)
            if score is None:
                continue
            out.scores.append(score)
            out.labels.append(label)
            meta = cohort.metadata.get(sample_id)
            out.studies.append(meta.study_name if meta else "unknown")
            out.sample_ids.append(sample_id)

    reporter.record(
        f"    {profile.name} vs {condition}: {out.n_cases} cases, {out.n_controls} controls"
    )
    return out


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class AucResult:
    profile: str
    condition: str
    n_cases: int
    n_controls: int
    auc: float | None
    ci_low: float | None
    ci_high: float | None
    per_study: dict[str, Any] = field(default_factory=dict)
    holdout_auc: float | None = None
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "profile": self.profile,
            "condition": self.condition,
            "n_cases": self.n_cases,
            "n_controls": self.n_controls,
            "auc": None if self.auc is None else round(self.auc, 4),
            "auc_ci95": (
                None
                if self.ci_low is None
                else [round(self.ci_low, 4), round(self.ci_high or 0.0, 4)]
            ),
            "study_holdout_mean_auc": (
                None if self.holdout_auc is None else round(self.holdout_auc, 4)
            ),
            "per_study": self.per_study,
            "note": " ".join(self.note.split()),
        }


def evaluate_auc(scores: LabelledScores, *, note: str = "") -> AucResult:
    """Pooled AUC plus leave-one-study-out, which is the honest number."""
    auc = roc_auc(scores.scores, scores.labels)
    low, high = auc_confidence_interval(scores.scores, scores.labels)

    per_study: dict[str, Any] = {}
    holdouts: list[float] = []
    for study in sorted(set(scores.studies)):
        indices = [i for i, s in enumerate(scores.studies) if s == study]
        subset_scores = [scores.scores[i] for i in indices]
        subset_labels = [scores.labels[i] for i in indices]
        if sum(subset_labels) == 0 or sum(subset_labels) == len(subset_labels):
            per_study[study] = {
                "n": len(indices),
                "auc": None,
                "note": "single-class study; AUC undefined",
            }
            continue
        study_auc = roc_auc(subset_scores, subset_labels)
        per_study[study] = {
            "n": len(indices),
            "n_cases": sum(subset_labels),
            "auc": None if study_auc is None else round(study_auc, 4),
        }
        if study_auc is not None:
            holdouts.append(study_auc)

    return AucResult(
        profile=scores.profile,
        condition=scores.condition,
        n_cases=scores.n_cases,
        n_controls=scores.n_controls,
        auc=auc,
        ci_low=low,
        ci_high=high,
        per_study=per_study,
        holdout_auc=statistics.fmean(holdouts) if holdouts else None,
        note=note,
    )


@dataclass(slots=True)
class SpecificityMatrix:
    """Every profile scored against every labelled condition."""

    rows: dict[str, dict[str, Any]] = field(default_factory=dict)

    def add(self, profile: str, condition: str, result: AucResult) -> None:
        self.rows.setdefault(profile, {})[condition] = result.to_json()

    def to_json(self) -> dict[str, Any]:
        return {
            "matrix": self.rows,
            "interpretation": (
                "A profile with a high AUC against its own condition and near 0.5 against "
                "the others is measuring that condition. A profile with similar AUCs across "
                "many conditions is measuring general gut disturbance, and its score should "
                "be read alongside the dysbiosis anchor rather than as condition-specific."
            ),
        }

    def print_table(self) -> None:
        conditions = sorted({c for row in self.rows.values() for c in row})
        if not conditions:
            print("  (no labelled conditions scored)")
            return
        print(f"  {'profile':<14}" + "".join(f"{c[:11]:>13}" for c in conditions))
        print("  " + "-" * (14 + 13 * len(conditions)))
        for profile in sorted(self.rows):
            line = f"  {profile:<14}"
            for condition in conditions:
                entry = self.rows[profile].get(condition)
                auc = entry.get("auc") if entry else None
                line += f"{('—' if auc is None else f'{auc:.2f}'):>13}"
            print(line)


def spread_verdict(matrix: SpecificityMatrix, profile: str) -> str:
    """Is this profile condition-specific, or measuring general dysbiosis?"""
    row = matrix.rows.get(profile, {})
    aucs = {c: e.get("auc") for c, e in row.items() if e.get("auc") is not None}
    if len(aucs) < 2:
        return "not enough labelled conditions to judge specificity"
    # cMD condition labels are upper case ("CRC"); profile names are lower.
    key = profile.strip().lower()
    own = next((v for c, v in aucs.items() if c.strip().lower() == key), None)
    others = [v for c, v in aucs.items() if c.strip().lower() != key]
    if own is None:
        best = max(others)
        worst = min(others)
        if best - worst <= 0.10 and best >= 0.65:
            return (
                f"no labelled cohort exists for this condition, so its own discrimination is "
                f"UNMEASURED. It does score {best:.2f} against another condition, which means "
                "its features are not specific to it"
            )
        return (
            f"no labelled cohort exists for this condition, so its own discrimination is "
            f"UNMEASURED (it scores {worst:.2f}-{best:.2f} against other conditions)"
        )
    if not others:
        return "no other labelled conditions to compare against"
    best_other = max(others)
    best_other_name = max(
        (c for c in aucs if c.strip().lower() != key), key=lambda c: aucs[c]
    )
    if own - best_other >= 0.10:
        return (
            f"condition-specific: AUC {own:.2f} against its own condition versus "
            f"{best_other:.2f} for the next highest ({best_other_name})"
        )
    if own < best_other - 0.05:
        return (
            f"NOT condition-specific: it separates {best_other_name} better "
            f"({best_other:.2f}) than the condition it was built for ({own:.2f}). These "
            "features are a general-disturbance signal, and the score must be read beside "
            "the dysbiosis anchor rather than as condition-specific"
        )
    if math.isclose(own, best_other, abs_tol=0.05):
        return (
            f"NOT condition-specific: AUC {own:.2f} against its own condition is "
            f"indistinguishable from {best_other:.2f} against {best_other_name}. This "
            "profile is probably measuring general dysbiosis"
        )
    return (
        f"weakly specific: AUC {own:.2f} own versus {best_other:.2f} for "
        f"{best_other_name} — the margin is small"
    )
