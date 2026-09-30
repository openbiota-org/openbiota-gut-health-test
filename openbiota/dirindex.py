"""Directional pattern index: the numeric kernel shared by the research panels.

One transformation, used by every study-derived panel in the library that
scores a sample against a published list of "higher in cases" / "lower in
cases" organisms. Chronic hives (:mod:`openbiota.cuindex`) and the inflammatory
skin diseases (:mod:`openbiota.skin`) both call it, and both get the same
guarantees, because the arithmetic and the missingness policy are the part
that must not drift between releases.

What the number is
------------------
A 0-100 index. Fifty is approximately centred on the reference after
transformation. Higher means movement in the direction the source study
reported for its case group; lower means the opposite. It is not a
percentage, not a probability of disease, and not a calibrated boundary
between health and illness, so it never carries a ``%`` sign and never
acquires a "positive" threshold.

For an abundance fraction ``a``::

    x            = log10(a + 1e-6)
    z            = clip((x - mu) / scale, -6, +6)
    contribution = 50 + 50 * tanh(direction * z / 2)

with ``mu`` and ``scale`` fitted offline from reference rows only. The panel
index is the equal-weight mean over usable point contributions; every
declared slot stays in the coverage denominator, so an unmeasurable organism
lowers coverage and widens the bounds rather than scoring zero.

Equal weights are deliberate. Source coefficients are recorded as provenance
and never used as weights: a coefficient fitted on a discovery cohort is a
description of that cohort, not a deployable classifier.

Bounds are coverage and censoring, not confidence
-------------------------------------------------
``full_panel_bounds`` is where the index could sit once the slots that could
not be measured are accounted for. It says nothing about sampling error, and
it is not a confidence interval around the point. When coverage is partial
the two have genuinely different scopes -- the point is over the measured
mask, the bounds are over the whole declared panel -- so a point can sit
outside a tightened bound without any arithmetic being wrong. Callers must
label the two separately.

Fitting is offline in the sense that matters: no query sample ever enters the
fit, so a sample scored alone and the same sample scored in a batch give
identical numbers.

One deliberate departure from the source formulas
-------------------------------------------------
Both specifications take each feature's median over every control row, with
nondetections sitting at ``log10(EPS)``. For an organism carried by only part
of a cohort that is the same thing as treating a control nondetection as an
observed zero, which the specifications forbid on the query side and require
be handled identically on the control side. It also degenerates: the centre
lands on the floor, the median absolute deviation collapses, and every
carrier clips to the rail, turning a graded feature into a full-weight
presence switch. :func:`fit_reference_detected` is therefore the production
estimator -- it fits on the detected control values, so the question becomes
where this sample sits among the participants who carry the organism.
:func:`fit_reference` is kept verbatim for the specifications' own numerical
fixtures.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable, Mapping, Sequence, Set
from dataclasses import dataclass, field, replace
from typing import Any, Final

from openbiota.errors import OpenBiotaError

# --------------------------------------------------------------------------- #
# constants — every one an engineering choice
# --------------------------------------------------------------------------- #

#: Added before the log so a zero abundance is representable.
EPS: Final = 1e-6
#: Floor on the robust scale, so a feature that is constant across controls
#: cannot divide by zero or turn rounding noise into a large deviation.
SCALE_FLOOR: Final = 0.25
#: Numerical stabilisation of the standardised value. Not a normal range and
#: not an out-of-distribution test.
Z_CLIP: Final = 6.0
#: Softness of the squash from standard deviations onto 0–100.
TEMPERATURE: Final = 2.0
#: Engineering minimum of distinct reference participants per fitted feature.
#: Twenty controls does not confer clinical validity.
MIN_CONTROLS: Final = 20
#: Minimum share of reference participants in whom a feature was detected
#: before its statistics may be used at all. Deliberately permissive: the
#: binding constraint is :data:`MIN_CONTROLS` applied to the *detected*
#: values, which for a large cohort bites first.
MIN_REFERENCE_PREVALENCE: Final = 0.0
#: 1 / Phi^-1(3/4): scales the median absolute deviation to a standard
#: deviation for a normal distribution.
MAD_TO_SIGMA: Final = 1.4826

SCORE_SEMANTICS: Final = "unvalidated_directional_pattern_index"
SCORE_SEMANTICS_TRANSPORTED: Final = "unvalidated_16S_to_shotgun_directional_pattern_index"
BOUNDS_TYPE: Final = "missingness_and_censoring_not_confidence_interval"

#: States a single feature measurement can be in. Zero is not one of them:
#: a profiler that does not report a species has not measured its absence.
MEASUREMENT_STATES: Final = (
    "detected",
    "below_reporting_limit",
    "not_in_database",
    "unresolved_mapping",
    "reference_ineligible",
    "assay_not_run",
    "failed_qc",
)

UNAVAILABLE_STATES: Final = (
    "no_usable_features",
    "qc_failed",
    "assay_incompatible",
    "population_not_supported",
    "reference_unavailable",
    "analytic_validation_required",
)


class DirIndexError(OpenBiotaError):
    """A malformed abundance, direction, panel or frozen reference."""


# --------------------------------------------------------------------------- #
# the kernel (spec §7), dependency-free
# --------------------------------------------------------------------------- #


def fraction(value: Any) -> float:
    """Validate an abundance as a finite fraction in [0, 1].

    A boolean is rejected rather than silently read as 0 or 1, and so is a
    percentage above one: units are declared by the caller, never guessed.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        # A boolean would read as 0 or 1 and a string would parse silently.
        # Schema validation belongs before the kernel, so this is a hard stop.
        raise DirIndexError(f"abundance must be a real number, got {value!r}")
    value = float(value)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise DirIndexError(f"expected a finite fraction in [0, 1], got {value!r}")
    return value


def fit_reference(
    rows: Sequence[Mapping[str, float | None]],
    feature_ids: Iterable[str],
) -> dict[str, dict[str, float]]:
    """Freeze ``mu`` and ``scale`` per feature from control rows only.

    One row is one distinct participant's first eligible baseline sample. The
    caller is responsible for that identity check: repeated rows would pass
    the count here while carrying no independent information.
    """
    if len(rows) < MIN_CONTROLS:
        raise DirIndexError(f"at least {MIN_CONTROLS} distinct eligible controls required")
    out: dict[str, dict[str, float]] = {}
    for fid in feature_ids:
        values = [
            math.log10(fraction(row[fid]) + EPS)
            for row in rows
            if fid in row and row[fid] is not None
        ]
        if len(values) < MIN_CONTROLS:
            continue
        mu = statistics.median(values)
        mad = statistics.median([abs(v - mu) for v in values])
        out[fid] = {
            "mu": mu,
            "scale": max(MAD_TO_SIGMA * mad, SCALE_FLOOR),
            "n_controls": float(len(values)),
        }
    return out


def fit_reference_detected(
    rows: Sequence[Mapping[str, float | None]],
    feature_ids: Iterable[str],
) -> dict[str, dict[str, float]]:
    """Fit ``mu`` and ``scale`` from the control values that were *detected*.

    This is the production fit, and it differs from :func:`fit_reference`
    deliberately. That function follows the specification's formula literally,
    taking the median over every control row with the nondetections sitting at
    ``log10(EPS)``. For a gut species carried by only some people that is the
    same thing as treating a nondetection as an observed zero, which §6.1 of
    the specification forbids on the query side and requires be handled
    identically on the control side.

    It also does not work. Where most controls are nondetections the fitted
    centre lands on the floor, the median absolute deviation collapses to
    zero, the scale drops to :data:`SCALE_FLOOR`, and every standardised value
    clips to ±6 — so the feature stops being a graded position and becomes a
    presence/absence switch carrying full weight. Six of the fifteen chronic
    hives species behave that way against this pipeline's reference cohort,
    including the one with the largest reported effect.

    Fitting on the detected values instead asks a question that has an answer:
    among the reference participants who carry this organism, where does this
    sample sit? A sample that does not carry it is not scored against that
    distribution at all — it goes down the censoring path and keeps only its
    bounds. A feature with fewer than :data:`MIN_CONTROLS` detected controls is
    left unfitted, so it lowers coverage rather than resting on a handful of
    carriers.
    """
    if len(rows) < MIN_CONTROLS:
        raise DirIndexError(f"at least {MIN_CONTROLS} distinct eligible controls required")
    out: dict[str, dict[str, float]] = {}
    for fid in feature_ids:
        values = [
            math.log10(fraction(row[fid]) + EPS)
            for row in rows
            if fid in row and row[fid] is not None and float(row[fid]) > 0.0
        ]
        if len(values) < MIN_CONTROLS:
            continue
        mu = statistics.median(values)
        mad = statistics.median([abs(v - mu) for v in values])
        out[fid] = {
            "mu": mu,
            "scale": max(MAD_TO_SIGMA * mad, SCALE_FLOOR),
            "n_controls": float(len(values)),
            "fitted_on": 1.0,  # 1 = detected values only
        }
    return out


def z_value(abundance: Any, ref: Mapping[str, float]) -> float:
    """Standardise one abundance against a frozen reference, clipped."""
    mu, scale = float(ref["mu"]), float(ref["scale"])
    if not math.isfinite(mu) or not math.isfinite(scale) or scale < SCALE_FLOOR:
        raise DirIndexError("invalid frozen reference")
    z = (math.log10(fraction(abundance) + EPS) - mu) / scale
    return max(-Z_CLIP, min(Z_CLIP, z))


def support(z: float, direction: Any) -> float:
    """Map a signed standardised value onto 0–100 alignment."""
    if isinstance(direction, bool) or direction not in (-1, 1):
        raise DirIndexError("direction must be -1 or +1")
    if not math.isfinite(z):
        raise DirIndexError("nonfinite standardised feature")
    return 50.0 + 50.0 * math.tanh(direction * z / TEMPERATURE)


def summarize(
    contributions: Mapping[str, tuple[float, float, float]],
    panel_ids: Sequence[str],
    *,
    interval_bounds: Mapping[str, tuple[float, float]] | None = None,
) -> dict[str, Any]:
    """Average the usable contributions and bound the unusable ones.

    ``interval_bounds`` carries interval-only features: they tighten the
    full-panel bounds without adding a point. The point index and the range
    therefore have different scopes, and the index may legitimately lie
    outside the range (spec v4.2 V42-61).
    """
    if not panel_ids or len(set(panel_ids)) != len(panel_ids):
        raise DirIndexError("panel requires unique IDs")
    intervals = interval_bounds or {}
    if set(intervals) & set(contributions):
        raise DirIndexError("a feature cannot be both a point and interval-only")
    usable = [contributions[f] for f in panel_ids if f in contributions]
    censored = [intervals[f] for f in panel_ids if f in intervals]
    k, m = len(panel_ids), len(usable)
    unavailable = k - m - len(censored)
    point = sum(v[0] for v in usable) / m if m else None
    lower = (sum(v[1] for v in usable) + sum(b[0] for b in censored)) / k
    upper = (sum(v[2] for v in usable) + sum(b[1] for b in censored) + 100.0 * unavailable) / k
    return {
        "index": point,
        "usable_count": m,
        "panel_count": k,
        "coverage": m / k,
        "full_panel_bounds": [lower, upper],
        "mask": [f for f in panel_ids if f in contributions],
    }


def evaluate(
    features: Sequence[Mapping[str, Any]],
    observations: Mapping[str, Mapping[str, float]],
    reference: Mapping[str, Mapping[str, float]],
) -> dict[str, Any]:
    """Score one sample against one panel and one frozen reference.

    ``observations`` holds only accepted measurements, each an interval
    ``{point, lower, upper}``. A nondetection with a validated reporting limit
    but no validated point estimator is ``point=None`` over ``[0, L]``: it is
    *interval-only*, tightens the full-panel bounds and contributes no point
    (spec v4.2 V42-24). A nondetection without a limit is not present here at
    all, because even its bounds would be a guess.
    """
    ids = [str(f["id"]) for f in features]
    if not ids or len(set(ids)) != len(ids):
        raise DirIndexError("panel requires unique IDs")
    if any(type(f["direction"]) is not int or f["direction"] not in (-1, 1) for f in features):
        raise DirIndexError("invalid direction")
    contributions: dict[str, tuple[float, float, float]] = {}
    interval_bounds: dict[str, tuple[float, float]] = {}
    details: dict[str, dict[str, Any]] = {}
    for feature in features:
        fid = str(feature["id"])
        if fid not in observations or fid not in reference:
            details[fid] = {"state": "unavailable", "point": None, "bounds": (0.0, 100.0)}
            continue
        obs = observations[fid]
        low, high = fraction(obs["lower"]), fraction(obs["upper"])
        if low > high:
            raise DirIndexError(f"{fid}: reversed interval")
        direction = feature["direction"]
        # Direction flips the interval, so the endpoints are reordered after
        # transformation rather than before it.
        ends = [support(z_value(v, reference[fid]), direction) for v in (low, high)]
        bounds = (min(ends), max(ends))
        raw_point = obs.get("point")
        if raw_point is None:
            interval_bounds[fid] = bounds
            details[fid] = {"state": "interval_only", "point": None, "bounds": bounds}
            continue
        point = fraction(raw_point)
        if not low <= point <= high:
            raise DirIndexError(f"{fid}: point must lie inside its interval")
        centre = support(z_value(point, reference[fid]), direction)
        contributions[fid] = (centre, *bounds)
        details[fid] = {"state": "point", "point": centre, "bounds": bounds}
    result = summarize(contributions, ids, interval_bounds=interval_bounds)
    result["contributions"] = contributions
    result["interval_bounds"] = interval_bounds
    result["interval_coverage"] = sum(d["state"] != "unavailable" for d in details.values()) / len(ids)
    result["details"] = details
    return result


def weighted_midrank(score: float, scores: Sequence[float], weights: Sequence[float]) -> float:
    """Mid-rank percentile of ``score`` within one named, weighted cohort."""
    if not scores or len(scores) != len(weights):
        raise DirIndexError("nonempty equal-length arrays required")
    if not math.isfinite(score) or any(not math.isfinite(s) for s in scores):
        raise DirIndexError("nonfinite score")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise DirIndexError("positive finite weights required")
    numerator = sum(
        w * (1.0 if s < score else 0.5 if s == score else 0.0)
        for s, w in zip(scores, weights, strict=True)
    )
    return 100.0 * numerator / sum(weights)


# --------------------------------------------------------------------------- #
# panels
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class PanelFeature:
    """One panel slot: the source study's name, and how it is measured here."""

    #: The immutable identifier the source paper used.
    source_id: str
    #: +1 for higher in cases, −1 for lower in cases.
    direction: int
    #: The rank the source reported this feature at.
    rank: str = "species"
    #: The production catalogue name, when it differs from ``source_id``.
    measured_as: str | None = None
    #: Why that equivalence holds.
    mapping_reason: str | None = None
    #: Non-overlapping children summed for a genus slot, if any. A member may
    #: be written ``"Old_name|New_name"`` when one organism has been renamed:
    #: whichever spelling the pinned catalogue carries is used, and a
    #: catalogue carrying both makes the slot unresolvable rather than summed.
    aggregate_of: tuple[str, ...] = ()
    #: Verified nomenclatural synonyms for the *same* organism (a formal
    #: reclassification, never a fuzzy match). Resolution picks the one the
    #: pinned catalogue carries; if it carries more than one, the slot is
    #: withheld, because two entries cannot be proven disjoint from here.
    synonyms: tuple[str, ...] = ()
    #: BH-adjusted q from the source table, where the source reports one.
    q_value: float | None = None
    #: The source's association coefficient. Provenance only: it is never a
    #: weight here, because a coefficient fitted on a discovery cohort is not
    #: a deployable classifier.
    coef: float | None = None

    @property
    def key(self) -> str:
        return self.measured_as or self.source_id

    def resolve(self, catalogue: Sequence[str] | Set[str]) -> tuple[tuple[str, ...], str | None]:
        """Which catalogue entries this slot is measured from, or why none.

        Returns ``(members, reason)``. ``members`` is empty when the slot
        cannot be measured from this catalogue, and ``reason`` then says why:
        ``not_in_catalogue`` or ``ambiguous_synonyms``. Matching is exact and
        by declared name only; there is no fuzzy or epithet matching.
        """
        names = set(catalogue)
        if self.aggregate_of:
            members: list[str] = []
            for spec in self.aggregate_of:
                alternatives = [n for n in spec.split("|") if n in names]
                if len(alternatives) > 1:
                    return (), "ambiguous_synonyms"
                members.extend(alternatives)
            return (tuple(members), None) if members else ((), "not_in_catalogue")
        if self.rank == "genus":
            prefix = f"{self.source_id}_"
            members_g = tuple(sorted(t for t in names if t.startswith(prefix)))
            return (members_g, None) if members_g else ((), "not_in_catalogue")
        candidates = [self.key, *self.synonyms]
        present = [n for n in candidates if n in names]
        if len(present) > 1:
            return (), "ambiguous_synonyms"
        return (tuple(present), None) if present else ((), "not_in_catalogue")

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.source_id,
            "measured_as": self.measured_as,
            "measurement_key": self.key,
            "synonyms": list(self.synonyms),
            "mapping_reason": self.mapping_reason,
            "aggregate_of": list(self.aggregate_of),
            "rank": self.rank,
            "direction": self.direction,
            "q_value": self.q_value,
            "coef": self.coef,
            "coef_is_a_weight": False,
        }


@dataclass(frozen=True, slots=True)
class Panel:
    """A named feature manifest with its own score semantics."""

    panel_id: str
    revision: str
    label: str
    source_id: str
    source_assay: str
    scoring_assay: str
    features: tuple[PanelFeature, ...]
    score_semantics: str
    discovery_sample_count: int | None = None
    assay_transport: str = "native"
    weight_rule: str = "equal"
    #: The clinical family this panel reads for. Several panels may share one
    #: family without the family becoming several diseases.
    family: str | None = None
    #: The comparison the source study actually made. A pooled contrast can
    #: never be reported as a subtype-specific headline.
    contrast: str | None = None
    #: How the production measurement is defined for this panel.
    measurement: str = "shotgun_taxonomic_fraction"
    #: How strong the source evidence is, in the source's own terms.
    evidence_tier: str = "exploratory_single_center_discovery"
    #: Whether this pattern has been shown to separate its condition from
    #: other diseases, medications and populations. Almost never yes.
    disease_specificity: str = "not_established"
    #: Sources supporting the panel and sources arguing against it.
    source_ids: tuple[str, ...] = ()
    counterevidence_ids: tuple[str, ...] = ()
    #: Inclusive participant-facing age bounds for a numeric result.
    min_age_years: float = 18.0
    max_age_years: float = 65.0
    #: Why no percentile is offered on this scale.
    percentile_note: str | None = None
    #: Panels whose features come from the same participants as this one, and
    #: which therefore are not independent replications of it.
    same_participants_as: tuple[str, ...] = ()

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(f.source_id for f in self.features)

    @property
    def all_source_ids(self) -> tuple[str, ...]:
        return self.source_ids or (self.source_id,)

    def subset(self, panel_id: str, label: str, *, max_q: float) -> Panel:
        """The features whose source q-value is below ``max_q``."""
        kept = tuple(f for f in self.features if f.q_value is not None and f.q_value < max_q)
        return replace(self, panel_id=panel_id, label=label, features=kept)



# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class FeatureRow:
    """One panel feature as it came out for this sample."""

    feature: PanelFeature
    state: str
    abundance: float | None = None
    reporting_limit: float | None = None
    contribution: float | None = None
    contribution_bounds: tuple[float, float] | None = None
    reference_median_abundance: float | None = None
    reference_prevalence: float | None = None
    n_controls: int | None = None
    note: str | None = None

    @property
    def usable(self) -> bool:
        return self.contribution is not None

    def to_json(self) -> dict[str, Any]:
        return {
            **self.feature.to_json(),
            "state": self.state,
            "abundance_fraction": self.abundance,
            "units": "relative abundance as a fraction of the total community",
            "reporting_limit": self.reporting_limit,
            "contribution": self.contribution,
            "contribution_bounds": list(self.contribution_bounds) if self.contribution_bounds else None,
            "reference_median_abundance": self.reference_median_abundance,
            "reference_prevalence": self.reference_prevalence,
            "n_controls": self.n_controls,
            "note": self.note,
        }


@dataclass(slots=True)
class PanelResult:
    """A typed result for one registered panel task, computed or not."""

    panel: Panel
    status: str
    score: float | None = None
    coverage: float = 0.0
    usable_count: int = 0
    full_panel_bounds: tuple[float, float] = (0.0, 100.0)
    feature_mask: tuple[str, ...] = ()
    rows: tuple[FeatureRow, ...] = ()
    reference_bundle_id: str | None = None
    reference_scope: str | None = None
    reason_codes: tuple[str, ...] = ()
    sensitivity: PanelResult | None = None
    notes: tuple[str, ...] = field(default=())

    @property
    def computed(self) -> bool:
        return self.status in ("computed", "computed_partial")

    @property
    def limited_feature(self) -> bool:
        """One usable feature out of many is labelled, never promoted."""
        return self.usable_count == 1 and self.panel_count > 1

    @property
    def panel_count(self) -> int:
        return len(self.panel.features)

    @property
    def interval_coverage(self) -> float:
        """Share of slots that produced at least a bounded contribution.

        Distinct from :attr:`coverage`, which counts point estimates. A
        censored nondetection tightens the bounds without contributing a
        point, and that is worth reporting rather than hiding.
        """
        if not self.rows:
            return 0.0
        bounded = sum(1 for r in self.rows if r.contribution_bounds is not None)
        return bounded / len(self.panel.features)

    def to_json(self) -> dict[str, Any]:
        return {
            "profile_id": self.panel.panel_id,
            "profile_revision": self.panel.revision,
            "label": self.panel.label,
            "status": self.status,
            "score_semantics": self.panel.score_semantics,
            "score_0_100": self.score,
            "index": self.score,
            "score_is_not_a_percentage": True,
            "coverage": self.coverage,
            "usable_count": self.usable_count,
            "panel_count": self.panel_count,
            "limited_feature_result": self.limited_feature,
            "full_panel_bounds": list(self.full_panel_bounds),
            "bounds_type": BOUNDS_TYPE,
            "strict_subset_result": self.sensitivity.to_json() if self.sensitivity else None,
            "control_percentile": None,
            "case_percentile": None,
            "percentile_note": self.panel.percentile_note or (
                "No held-out compatible case or control score distribution exists for this panel, "
                "so no percentile is reported on this scale. A panel selected in the same "
                "participants that produced it can only be ranked against that cohort, which is "
                "reconstruction rather than external validation."
            ),
            "reference_bundle_id": self.reference_bundle_id,
            "reference_scope": self.reference_scope,
            "feature_mask": list(self.feature_mask),
            "assay_transport": self.panel.assay_transport,
            "source_assay": self.panel.source_assay,
            "scoring_assay": self.panel.scoring_assay,
            "evidence_maturity": self.panel.evidence_tier,
            "evidence_tier": self.panel.evidence_tier,
            "disease_specificity": self.panel.disease_specificity,
            "profile_family": self.panel.family,
            "contrast": self.panel.contrast,
            "measurement": self.panel.measurement,
            "index_label": "Research directional pattern index",
            "index_scope": "available_point_feature_mask",
            "bounds_scope": "full_declared_panel",
            "interval_coverage": self.interval_coverage,
            "weight_rule": self.panel.weight_rule,
            "reason_codes": list(self.reason_codes),
            "source_ids": list(self.panel.all_source_ids),
            "counterevidence_ids": list(self.panel.counterevidence_ids),
            "not_independent_of": list(self.panel.same_participants_as),
            "disease_probability": None,
            "clinical_classification": None,
            "observations": [r.to_json() for r in self.rows],
            "notes": list(self.notes),
            "treatment_recommendation": None,
            "donor_eligibility": None,
        }


# --------------------------------------------------------------------------- #
# fitting against this pipeline's reference cohort
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class FrozenReference:
    """``mu``/``scale`` per feature plus the provenance of the rows they came from."""

    bundle_id: str
    scope: str
    stats: dict[str, dict[str, float]]
    median_abundance: dict[str, float] = field(default_factory=dict)
    prevalence: dict[str, float] = field(default_factory=dict)
    n_participants: int = 0
    #: source_id -> the exact catalogue entries each slot was fitted from.
    #: Provenance of the name resolution, frozen with the statistics.
    resolved: dict[str, tuple[str, ...]] = field(default_factory=dict)
    #: source_id -> why a slot could not be resolved against the catalogue.
    unresolved: dict[str, str] = field(default_factory=dict)

    def __contains__(self, key: str) -> bool:
        return key in self.stats

    def members(self, feature: PanelFeature) -> tuple[str, ...] | None:
        """Catalogue entries for ``feature``; ``None`` when never resolved.

        References frozen before resolution was recorded fall back to the
        feature's declared key so that older bundles still score.
        """
        if feature.source_id in self.resolved:
            return self.resolved[feature.source_id]
        if feature.source_id in self.unresolved:
            return ()
        return None

    def eligible(self, key: str) -> bool:
        """Whether this feature's fitted statistics may be used at all.

        Two ways a slot fails: too few reference participants carried the
        organism for a distribution to exist, or the fitted centre is itself
        a nondetection, which would make the feature a presence/absence switch
        rather than a graded position.
        """
        stats = self.stats.get(key)
        if stats is None:
            return False
        if stats.get("n_controls", 0.0) < MIN_CONTROLS:
            return False
        if stats["mu"] <= math.log10(EPS) + 1e-9:
            return False
        return self.prevalence.get(key, 0.0) >= MIN_REFERENCE_PREVALENCE

    def carriers(self, key: str) -> int | None:
        stats = self.stats.get(key)
        return None if stats is None else int(stats.get("n_controls", 0.0))

    def ineligible_reason(self, key: str) -> str:
        prevalence = self.prevalence.get(key)
        share = "" if prevalence is None else f" It was detected in {prevalence:.1%} of them."
        carriers = self.carriers(key)
        if carriers is None or carriers < MIN_CONTROLS:
            return (
                f"fewer than {MIN_CONTROLS} reference participants carried this organism, so there "
                "is no distribution to place this sample in and the slot lowers coverage rather "
                "than contributing a number." + share
            )
        return (
            "the reference centre for this organism is itself a nondetection, so standardising "
            "against it would turn the feature into a presence/absence switch carrying full "
            "weight. The slot is reported as reference-ineligible: it lowers coverage and widens "
            "the bounds instead." + share
        )


def _abundance_rows(
    catalogue: Sequence[str],
    matrix: Sequence[Sequence[float]],
    keys: Mapping[str, tuple[str, ...]],
) -> list[dict[str, float]]:
    """Per-participant fractions for each wanted key, summing genus members.

    ``matrix`` is taxa-by-samples in percent, the shape the reference cohort
    ships in. Percentages become fractions here, once, by the declared unit
    rather than by inspecting the numbers.
    """
    index = {name: i for i, name in enumerate(catalogue)}
    n_samples = len(matrix[0]) if matrix else 0
    rows: list[dict[str, float]] = []
    for s in range(n_samples):
        row: dict[str, float] = {}
        for key, members in keys.items():
            positions = [index[m] for m in members if m in index]
            if not positions:
                continue
            total = sum(float(matrix[p][s]) for p in positions) / 100.0
            row[key] = min(max(total, 0.0), 1.0)
        rows.append(row)
    return rows


def freeze_reference(
    panel: Panel,
    *,
    catalogue: Sequence[str],
    matrix: Sequence[Sequence[float]],
    bundle_id: str,
    scope: str = "compatible_independent_controls",
) -> FrozenReference:
    """Fit the panel's reference statistics from cohort rows only.

    No query sample is passed here, ever. That is what makes a sample scored
    alone identical to the same sample scored in a batch.
    """
    keys: dict[str, tuple[str, ...]] = {}
    unresolved: dict[str, str] = {}
    for feature in panel.features:
        members, reason = feature.resolve(catalogue)
        if members:
            keys[feature.source_id] = members
        else:
            unresolved[feature.source_id] = reason or "not_in_catalogue"

    rows = _abundance_rows(catalogue, matrix, keys)
    stats = fit_reference_detected(rows, keys)
    median_abundance: dict[str, float] = {}
    prevalence: dict[str, float] = {}
    for key in keys:
        values = [r[key] for r in rows if key in r]
        if values:
            median_abundance[key] = statistics.median(values)
            prevalence[key] = sum(1 for v in values if v > 0.0) / len(values)
    return FrozenReference(
        bundle_id=bundle_id,
        scope=scope,
        stats=stats,
        median_abundance=median_abundance,
        prevalence=prevalence,
        n_participants=len(rows),
        resolved=dict(keys),
        unresolved=unresolved,
    )


def score_panel(
    panel: Panel,
    abundances: Mapping[str, float],
    reference: FrozenReference | None,
    *,
    reporting_limit: float | None = None,
    detected: Mapping[str, bool] | None = None,
    aggregates: Mapping[str, tuple[str, ...]] | None = None,
) -> PanelResult:
    """Score one sample's abundances against one panel.

    ``abundances`` are fractions of the total community, keyed by production
    catalogue name. A key that is absent means the profiler did not report the
    organism: that is ``not_in_database``, not zero. A key present at zero is
    a nondetection, censored over ``[0, reporting_limit]`` when a limit is
    supplied and dropped from the point estimate when it is not.
    """
    if reference is None:
        return PanelResult(
            panel=panel,
            status="reference_unavailable",
            reason_codes=("no_compatible_frozen_reference",),
            rows=tuple(
                FeatureRow(feature=f, state="assay_not_run", note="no frozen reference for this panel")
                for f in panel.features
            ),
        )

    rows: list[FeatureRow] = []
    features: list[dict[str, Any]] = []
    observations: dict[str, dict[str, float]] = {}

    for feature in panel.features:
        # Which catalogue entries this slot reads from. The frozen reference
        # carries the resolution it was fitted on, so the sample is measured
        # on exactly the same names as the controls were.
        members = (aggregates or {}).get(feature.source_id)
        if members is None:
            members = reference.members(feature)
        if members is None and feature.rank == "genus":
            members = tuple(k for k in abundances if k.startswith(f"{feature.source_id}_"))
        row = FeatureRow(
            feature=feature,
            state="not_in_database",
            reference_median_abundance=reference.median_abundance.get(feature.source_id),
            reference_prevalence=reference.prevalence.get(feature.source_id),
            n_controls=int(reference.stats.get(feature.source_id, {}).get("n_controls", 0)) or None,
        )

        if feature.source_id in reference.unresolved:
            # Never fitted, because the catalogue has no entry for it (or has
            # two that cannot be proven disjoint). Not zero, not a proxy.
            why = reference.unresolved[feature.source_id]
            row.state = "not_in_database"
            row.note = (
                "the pinned reference catalogue carries more than one entry that could be this "
                "organism, and they cannot be proven disjoint from here, so the slot is withheld "
                "rather than summed or guessed; it lowers coverage instead of scoring"
                if why == "ambiguous_synonyms" else
                "this organism has no slot in the pinned reference catalogue used here"
                + (
                    f" (checked under {', '.join(n.replace('_', ' ') for n in (feature.key, *feature.synonyms))})"
                    if feature.synonyms else ""
                )
                + ", so its absence has not been measured; the slot lowers coverage instead of scoring zero"
            )
            rows.append(row)
            continue

        if not reference.eligible(feature.source_id):
            row.state = "reference_ineligible"
            row.note = reference.ineligible_reason(feature.source_id)
            lookup = members if members else (feature.key,)
            present_now = [k for k in lookup if k in abundances]
            if present_now:
                row.abundance = min(max(sum(float(abundances[k]) for k in present_now), 0.0), 1.0)
            rows.append(row)
            continue

        if members is not None:
            present = [k for k in members if k in abundances]
            value: float | None = sum(abundances[k] for k in present) if present else None
            seen = any((detected or {}).get(k, abundances.get(k, 0.0) > 0.0) for k in present)
        elif feature.key in abundances:
            value = abundances[feature.key]
            seen = (detected or {}).get(feature.key, value > 0.0)
        else:
            value = None
            seen = False

        if value is None:
            # The profiler lists only what it detected. If the organism is in
            # the reference catalogue then the database can measure it and
            # this is a nondetection; if it is not, the database has no slot
            # for it and its absence was never measurable at all.
            in_database = feature.source_id in reference.median_abundance
            if in_database:
                row.state = "below_reporting_limit"
                prevalence = reference.prevalence.get(feature.source_id)
                share = "" if prevalence is None else f" It was detected in {prevalence:.0%} of the reference participants."
                row.note = (
                    "the profiler reports only the organisms it detects, and this one is not in "
                    "this sample's table. No validated reporting limit exists for it at this "
                    "sequencing depth, so its point is omitted and only its bounds are kept: a "
                    "nondetection, not a measured absence." + share
                )
            else:
                row.note = (
                    "this organism has no slot in the reference catalogue used here, so its "
                    "absence has not been measured; the slot lowers coverage instead of scoring zero"
                )
            rows.append(row)
            continue

        row.abundance = min(max(float(value), 0.0), 1.0)
        row.n_controls = reference.carriers(feature.source_id)
        if row.abundance > 0.0 and seen:
            row.state = "detected"
            observations[feature.source_id] = {
                "point": row.abundance,
                "lower": row.abundance,
                "upper": row.abundance,
            }
        elif reporting_limit is not None:
            row.state = "below_reporting_limit"
            row.reporting_limit = reporting_limit
            row.note = (
                f"not detected at a reporting limit of {reporting_limit:.2e} of the community; "
                "this is a detection statement, not biological absence. No validated point "
                "estimator exists for a censored value, so this slot tightens the panel bounds "
                "and contributes no point"
            )
            observations[feature.source_id] = {"point": None, "lower": 0.0, "upper": reporting_limit}
        else:
            row.state = "below_reporting_limit"
            row.note = (
                "not detected, and no validated reporting limit is available for this feature, so "
                "its point is omitted and only its bounds are kept"
            )
            rows.append(row)
            features.append({"id": feature.source_id, "direction": feature.direction})
            continue

        features.append({"id": feature.source_id, "direction": feature.direction})
        rows.append(row)

    # Features whose reference exists but which produced no observation still
    # belong in the panel: they are what widens the bounds.
    all_features = [{"id": f.source_id, "direction": f.direction} for f in panel.features]
    usable_stats = {k: v for k, v in reference.stats.items() if reference.eligible(k)}
    summary = evaluate(all_features, observations, usable_stats)

    by_id = {r.feature.source_id: r for r in rows}
    for fid, (point, low, high) in summary["contributions"].items():
        row = by_id.get(fid)
        if row is not None:
            row.contribution = point
            row.contribution_bounds = (low, high)
    for fid, (low, high) in summary["interval_bounds"].items():
        row = by_id.get(fid)
        if row is not None:
            row.contribution = None
            row.contribution_bounds = (low, high)

    usable = int(summary["usable_count"])
    status = "computed" if usable == len(panel.features) else "computed_partial"
    reasons: list[str] = []
    if usable == 0 and summary["interval_bounds"]:
        # Bounded measurements exist but no point does: index null, bounds
        # possibly informative (spec v4.2 §7.1 ``interval_only``).
        status = "interval_only"
        reasons.append("only_bounded_measurements_available")
    elif usable == 0:
        status = "no_usable_features"
        reasons.append("no_usable_features")
    elif usable < len(panel.features):
        reasons.append("partial_panel_coverage")

    notes: list[str] = []
    if usable == 1 and len(panel.features) > 1:
        notes.append(
            f"A limited-feature result: one of {len(panel.features)} panel features was usable, so "
            "the figure is one measurement rather than a panel and the bounds are correspondingly wide."
        )
    if panel.assay_transport != "native":
        notes.append(
            "Every figure for this panel crosses assay types: the source study named these "
            "organisms by 16S amplicon and this report measures them by shotgun sequencing. The "
            "transport is unvalidated."
        )

    return PanelResult(
        panel=panel,
        status=status,
        score=summary["index"],
        coverage=float(summary["coverage"]),
        usable_count=usable,
        full_panel_bounds=(summary["full_panel_bounds"][0], summary["full_panel_bounds"][1]),
        feature_mask=tuple(summary["mask"]),
        rows=tuple(rows),
        reference_bundle_id=reference.bundle_id,
        reference_scope=reference.scope,
        reason_codes=tuple(reasons),
        notes=tuple(notes),
    )



def comparable(a: PanelResult, b: PanelResult) -> bool:
    """Whether two results may be compared directly across dates.

    Same panel, same reference bundle and the same feature mask. Anything else
    needs an explicit recomputation on a named common mask.
    """
    return (
        a.panel.panel_id == b.panel.panel_id
        and a.reference_bundle_id == b.reference_bundle_id
        and a.feature_mask == b.feature_mask
    )


# --------------------------------------------------------------------------- #
# self-test (spec §7)
# --------------------------------------------------------------------------- #


def self_test() -> None:
    """The specification's numerical assertions, verbatim.

    The repeated rows below are synthetic numerical fixtures. Production
    enforces distinct participant identity before fitting; nothing here is
    reference data.
    """
    assert support(0, 1) == 50
    assert abs(support(2, 1) - 88.07970779778825) < 1e-10
    assert abs(support(2, -1) - 11.92029220221176) < 1e-10
    r = summarize({"a": (80.0, 70.0, 90.0)}, ["a", "b"])
    assert r["index"] == 80 and r["coverage"] == 0.5
    assert r["full_panel_bounds"] == [35.0, 95.0]
    assert summarize({}, ["a"])["index"] is None
    assert weighted_midrank(50, [20, 50, 50, 80], [1, 1, 1, 1]) == 50
    refs = fit_reference([{"a": 0.01}] * 20, ["a"])
    assert refs["a"]["scale"] == SCALE_FLOOR
    f = [{"id": "a", "direction": -1}]
    obs = {"a": {"point": 0.0, "lower": 0.0, "upper": 0.001}}
    result = evaluate(f, obs, refs)
    assert result["full_panel_bounds"][0] <= result["index"] <= result["full_panel_bounds"][1]
    assert result == evaluate(f, obs, refs)
    # V42-24: a censored nondetection without a point estimator is interval-only.
    censored = evaluate(f, {"a": {"point": None, "lower": 0.0, "upper": 0.001}}, refs)
    assert censored["index"] is None and censored["coverage"] == 0.0
    assert censored["interval_coverage"] == 1.0
    assert censored["details"]["a"]["state"] == "interval_only"
    # V42-61: a tightened range can exclude a partial-mask point index.
    two = [{"id": "a", "direction": -1}, {"id": "b", "direction": 1}]
    refs2 = {**refs, "b": refs["a"]}
    mixed = evaluate(two, {"a": {"point": 0.01, "lower": 0.01, "upper": 0.01},
                           "b": {"point": None, "lower": 0.0, "upper": 1e-6}}, refs2)
    assert mixed["coverage"] == 0.5 and mixed["interval_coverage"] == 1.0
    assert mixed["index"] == 50.0
    assert not mixed["full_panel_bounds"][0] <= mixed["index"] <= mixed["full_panel_bounds"][1]
    for bad in (True, -0.1, 1.1, float("nan"), float("inf")):
        try:
            fraction(bad)
        except DirIndexError:
            pass
        else:  # pragma: no cover - the guard is the point of the test
            raise AssertionError(f"invalid abundance accepted: {bad!r}")


if __name__ == "__main__":  # pragma: no cover
    self_test()
    print("directional index numerical self-test passed")
