"""Shape analysis across the profile library (spec 7.8) and its supporting
statistics: empirical specificity, the literature-uniqueness prior, competing
profiles, and the confounder ledger.

With twenty-odd profiles scoring, the individual numbers stop being the
headline. Most disease-associated microbiome changes are a shared response to
illness rather than disease-specific markers (Duvallet 2017), so a sample that
resembles *many* patterns is telling you about general disturbance, not about
many diseases. The shape — how many profiles score high, how specific those
profiles are known to be, and what the general dysbiosis index says — is what
separates a specific resemblance from a diffuse one, and the report leads with
it.
"""

from __future__ import annotations

import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from openbiota.profiles import Profile, ProfileSet
from openbiota.scoring import ProfileResult

#: Percentile thresholds the shape counts against. These match the report's
#: seven-band scale: 75 is the top of the typical band (the reference group's
#: middle half), 95 is where "notably high" begins.
HIGH_AT: float = 75.0
VERY_HIGH_AT: float = 95.0

#: Fraction of scored profiles above HIGH_AT beyond which the resemblance is
#: called diffuse rather than specific.
DIFFUSE_FRACTION: float = 0.4

#: Mean specificity of the high scorers below which they are treated as
#: general-disturbance markers.
LOW_SPECIFICITY: float = 0.5


# --------------------------------------------------------------------------- #
# specificity
# --------------------------------------------------------------------------- #


def literature_uniqueness_prior(profile: Profile, library: ProfileSet) -> float | None:
    """How unique this profile's features are across the registry.

    Mean over fused, bound features of 1 / (number of profiles using the same
    feature). A labelled fallback prior, not empirical specificity: it
    changes whenever the registry is curated, and a feature that many profiles
    share may still be discriminating in practice.
    """
    counts: dict[str, int] = {}
    for other in library.profiles:
        seen = {f.name for m in other.fused_modules for f in m.features if f.bound}
        for name in seen:
            counts[name] = counts.get(name, 0) + 1
    names = [f.name for m in profile.fused_modules for f in m.features if f.bound]
    if not names:
        return None
    return statistics.fmean(1.0 / max(1, counts.get(n, 1)) for n in names)


def empirical_specificity(
    profile_name: str, validation: Mapping[str, Any] | None
) -> tuple[float | None, str]:
    """Specificity from the challenge-cohort matrix, if the profile is in it.

    Uses the worst-case off-target discrimination: a profile that separates
    another disease from health as well as it separates its own is not
    specific, however good its own AUC looks. Returns ``(value, basis)``.
    """
    if not validation:
        return None, "no challenge-cohort matrix available"
    matrix = (validation.get("specificity_matrix") or {}).get("matrix") or {}
    row = matrix.get(profile_name)
    if not isinstance(row, dict) or not row:
        return None, "profile not yet scored against the labelled challenge cohorts"
    target = None
    off_target: list[float] = []
    for condition, entry in row.items():
        auc = entry.get("auc") if isinstance(entry, dict) else None
        if auc is None:
            continue
        if _is_own_condition(profile_name, condition):
            target = auc
            continue
        off_target.append(float(auc))
    if not off_target:
        return None, "no off-target challenge cohort scored"
    worst = max(off_target)
    excess = max(0.0, worst - 0.5) / 0.5
    value = max(0.0, min(1.0, 1.0 - excess))
    basis = (
        f"worst off-target AUROC {worst:.2f}"
        + (f"; own condition {target:.2f}" if target is not None else "; own condition unmeasured")
    )
    return value, basis


_CONDITION_ALIASES: dict[str, tuple[str, ...]] = {
    "crc": ("CRC",),
    "adenoma": ("adenoma",),
    "ibd": ("IBD",),
    "crohns": ("IBD",),
    "uc": ("IBD",),
    "t2d": ("T2D",),
}


def _is_own_condition(profile_name: str, condition: str) -> bool:
    return condition in _CONDITION_ALIASES.get(profile_name, ())


# --------------------------------------------------------------------------- #
# shape
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ShapeAnalysis:
    n_profiles: int
    n_scored: int
    n_abstained: int
    n_not_computable: int
    n_above_typical: int
    n_notably_high: int
    high_scorers: tuple[str, ...]
    mean_specificity_of_high: float | None
    specificity_basis: str
    dysbiosis_score: float | None
    dysbiosis_band: str | None
    verdict: str
    headline: str
    explanation: str
    top: tuple[tuple[str, float], ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "n_profiles": self.n_profiles,
            "n_scored": self.n_scored,
            "n_abstained": self.n_abstained,
            "n_not_computable": self.n_not_computable,
            "n_above_typical_band": self.n_above_typical,
            "n_notably_high": self.n_notably_high,
            "high_scorers": list(self.high_scorers),
            "mean_specificity_of_high_scorers": (
                None if self.mean_specificity_of_high is None
                else round(self.mean_specificity_of_high, 2)
            ),
            "specificity_basis": self.specificity_basis,
            "gmwi2_score": None if self.dysbiosis_score is None else round(self.dysbiosis_score, 3),
            "gmwi2_band": self.dysbiosis_band,
            "verdict": self.verdict,
            "headline": self.headline,
            "explanation": self.explanation,
            "top": [{"profile": n, "percentile": round(p, 1)} for n, p in self.top],
        }


def _specificity_for(result: ProfileResult) -> tuple[float | None, bool]:
    """(value, empirical?) — empirical where measured, else the labelled prior."""
    vector = result.vector
    if vector is None:
        return None, False
    if vector.empirical_specificity is not None:
        return vector.empirical_specificity, True
    return vector.literature_uniqueness_prior, False


def shape_analysis(results: Sequence[ProfileResult], anchor: Any | None) -> ShapeAnalysis:
    scored = [r for r in results if r.reportable and r.combined_percentile is not None]
    abstained = [r for r in results if r.abstention.abstained]
    not_computable = [
        r for r in results if not r.abstention.abstained and r.combined_percentile is None
    ]
    high = [r for r in scored if r.combined_percentile >= HIGH_AT]
    very_high = [r for r in scored if r.combined_percentile >= VERY_HIGH_AT]

    specs = [_specificity_for(r) for r in high]
    values = [v for v, _ in specs if v is not None]
    mean_spec = statistics.fmean(values) if values else None
    n_empirical = sum(1 for v, emp in specs if v is not None and emp)
    basis = (
        "no high scorers"
        if not high
        else (
            f"{n_empirical} of {len(high)} from challenge-cohort AUROCs, the rest from the "
            "literature-uniqueness prior"
            if n_empirical
            else "literature-uniqueness prior only (no challenge cohort covers these profiles)"
        )
    )

    dys_score = getattr(anchor, "score", None) if anchor is not None else None
    dys_band = getattr(anchor, "band", None) if anchor is not None else None

    fraction_high = len(high) / len(scored) if scored else 0.0
    if not scored:
        verdict = "unscored"
        headline = "No disease pattern could be scored."
        explanation = (
            "Every profile abstained or lacked a bound engine; see the reasons listed "
            "with each profile."
        )
    elif not high:
        verdict = "quiet"
        headline = "No published disease pattern stands out."
        explanation = (
            f"All {len(scored)} scored profiles sit within or below the typical band "
            f"(at or under the {HIGH_AT:.0f}th percentile of the reference group). The "
            "community does not resemble any of the registered patterns more than most "
            "healthy people do."
        )
    elif fraction_high >= DIFFUSE_FRACTION and (mean_spec is None or mean_spec < LOW_SPECIFICITY):
        verdict = "diffuse"
        headline = "A general-disturbance picture, not a specific pattern."
        explanation = (
            f"{len(high)} of {len(scored)} scored profiles sit above the typical band, and "
            "the features they rest on are shared across many conditions "
            + (f"(mean specificity {mean_spec:.2f})" if mean_spec is not None else "")
            + ". When most profiles score high at once, the honest reading is that the "
            "community is disturbed in the way many illnesses disturb it — the general "
            "dysbiosis index is the headline, not any one disease."
        )
    elif len(high) <= 2 and (mean_spec is None or mean_spec >= LOW_SPECIFICITY):
        verdict = "specific"
        headline = f"Resemblance concentrated in {len(high)} pattern{'s' if len(high) > 1 else ''}."
        explanation = (
            f"Only {len(high)} of {len(scored)} scored profiles sit above the typical band "
            "and their features are comparatively distinctive. That is the shape a specific "
            "resemblance takes — though resemblance is still not risk."
        )
    else:
        verdict = "mixed"
        headline = (
            f"{len(high)} of {len(scored)} patterns sit above the typical band; "
            "specificity is partial."
        )
        explanation = (
            "Several profiles score high and they share some but not all of their "
            "features. Read the ranked table with the specificity column: the profiles "
            "whose distinctive features drive the score matter more than those riding "
            "on shared butyrate-producer depletion or oral-taxa enrichment."
        )
    if dys_band is not None and verdict in ("diffuse", "mixed"):
        explanation += f" GMWI2 general dysbiosis index: {dys_score:+.2f} ({dys_band})."

    ranked = sorted(scored, key=lambda r: -(r.combined_percentile or 0))
    return ShapeAnalysis(
        n_profiles=len(results),
        n_scored=len(scored),
        n_abstained=len(abstained),
        n_not_computable=len(not_computable),
        n_above_typical=len(high),
        n_notably_high=len(very_high),
        high_scorers=tuple(r.profile.name for r in ranked if r in high),
        mean_specificity_of_high=mean_spec,
        specificity_basis=basis,
        dysbiosis_score=dys_score,
        dysbiosis_band=dys_band,
        verdict=verdict,
        headline=headline,
        explanation=explanation,
        top=tuple((r.profile.name, r.combined_percentile) for r in ranked[:5]),
    )


def annotate_competition(results: Sequence[ProfileResult]) -> None:
    """Fill in competing profiles and the top-two margin on every result.

    The margin is the gap between the library's highest percentile and its
    second highest. A large margin means one pattern stands apart; a small
    one means the top profile is not distinguishable from its neighbour on
    this sample, whatever its absolute percentile.
    """
    scored = sorted(
        (r for r in results if r.reportable and r.combined_percentile is not None),
        key=lambda r: -(r.combined_percentile or 0),
    )
    if not scored:
        return
    margin = (
        scored[0].combined_percentile - scored[1].combined_percentile
        if len(scored) > 1
        else None
    )
    for result in scored:
        others = [
            (r.profile.name, r.combined_percentile)
            for r in scored
            if r is not result and r.profile.family != result.profile.family
        ][:3]
        result.competing = tuple(others)
        result.top_two_margin = margin


# --------------------------------------------------------------------------- #
# confounder ledger (spec 12)
# --------------------------------------------------------------------------- #

CONFOUNDER_FIELDS: tuple[tuple[str, str, str], ...] = (
    ("antibiotics", "Antibiotics in the last 6 months", "dominates the community for months; profiles abstain within 90 days"),
    ("metformin", "Metformin", "produces much of the published type 2 diabetes signature; assumed absent when no medication list is supplied (participant mode abstains)"),
    ("ppi", "Proton-pump inhibitor", "large documented shift; recorded and lowers technical reliability"),
    ("gluten_free_diet", "Gluten-free diet", "large effect in prevalent celiac disease; required for the celiac profile"),
    ("country", "Country", "geography often exceeds the disease effect; reference matched on it"),
    ("age", "Age", "large across the lifespan; mandatory for age-specific profiles"),
    ("sex", "Sex", "androgenetic alopecia signals are sex-specific; mandatory there"),
    ("stool_form", "Stool form", "correlates with diversity; recorded"),
    ("hospitalised", "Hospitalisation / ICU", "drives E. coli, Klebsiella and resistance genes; recorded"),
    ("kidney_function", "Kidney function", "decides whether uremic-toxin capacity matters; context for CKD"),
    ("medications", "Other medications", "recorded"),
    ("sequencing_run", "Sequencing run / flowcell", "batch effects can exceed biology; recorded from FASTQ headers"),
)


@dataclass(frozen=True, slots=True)
class ConfounderLedger:
    recorded: tuple[tuple[str, str, str], ...]
    missing: tuple[tuple[str, str, str], ...]
    abstentions: tuple[tuple[str, str], ...]
    notes: tuple[str, ...] = field(default=())
    #: Confounders that were not supplied and were assumed absent so that the
    #: patterns depending on them could still be scored (see openbiota.context).
    assumed: tuple[tuple[str, str, str], ...] = field(default=())

    def to_json(self) -> dict[str, Any]:
        return {
            "recorded": [{"field": k, "label": lab, "value": v} for k, lab, v in self.recorded],
            "assumed_absent": [{"field": k, "label": lab, "effect": e} for k, lab, e in self.assumed],
            "missing": [{"field": k, "label": lab, "effect": e} for k, lab, e in self.missing],
            "profiles_abstained": [{"profile": p, "reason": r} for p, r in self.abstentions],
            "notes": list(self.notes),
        }


def confounder_ledger(
    metadata: Mapping[str, Any], results: Sequence[ProfileResult],
    assumed_keys: frozenset[str] = frozenset(),
) -> ConfounderLedger:
    recorded: list[tuple[str, str, str]] = []
    missing: list[tuple[str, str, str]] = []
    assumed: list[tuple[str, str, str]] = []
    for key, label, effect in CONFOUNDER_FIELDS:
        value = metadata.get(key)
        if key in assumed_keys:
            assumed.append((key, label, effect))
        elif value in (None, "", False):
            missing.append((key, label, effect))
        else:
            recorded.append((key, label, str(value)))
    abstentions = tuple(
        (r.profile.label, "; ".join(r.abstention.triggered))
        for r in results
        if r.abstention.abstained
    )
    return ConfounderLedger(
        recorded=tuple(recorded),
        missing=tuple(missing),
        assumed=tuple(assumed),
        abstentions=abstentions,
        notes=(
            "What was recorded, what was missing, and which profiles refused to score "
            "because of it. A missing confounder does not lower a score; it widens what "
            "the score could mean.",
        ),
    )
