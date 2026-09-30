"""Generic profile scoring. Knows nothing about which disease it is scoring.

Consumes a profile plus a sample's measurements and produces per-feature
values, module scores, percentiles against an untouched reference, a
confidence grade, and — where warranted — a refusal to score at all.

The arithmetic, from the spec:

    z_i = (CLR_i(sample) − median_i(ref)) / (1.4826 × MAD_i(ref))
    v_i = d_i × clip(z_i, −3, +3) / 3
    module = Σ(w q s c v) / Σ(w q s c)

Median and MAD rather than mean and SD, because microbiome abundance
distributions are skewed enough that one outlier otherwise sets the scale.
Clipping at ±3 SD then dividing by 3 maps each feature onto −1..+1, so no
single extreme feature can dominate a module.

Percentiles come from scoring every reference sample through the identical
path and ranking the sample against that null — never from linearly rescaling
the raw score, which would silently assume a distribution shape.
"""

from __future__ import annotations

import bisect
import math
import random
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from typing import Any, Final

from openbiota.profiles import (
    COMBINED_MODULES,
    MODULE_FUNCTIONAL,
    MODULE_TAXONOMIC,
    Feature,
    Module,
    Profile,
)

#: Clip bound for the robust z, in MAD-scaled standard deviations.
Z_CLIP: Final = 3.0

#: Below this many usable non-host read pairs a taxon's non-detection is
#: ``missing`` rather than absence (spec 4.2 / 7.2). Enforced by the
#: presence/absence adapter, not left to interpretation.
MIN_PAIRS_FOR_ABSENCE: Final = 500_000

#: Independent-cluster coverage below which a score is returned but flagged
#: ``scored_low_coverage`` with a widened interval.
LOW_COVERAGE_BELOW: Final = 0.5

#: Bootstrap draws for the uncertainty interval.
INTERVAL_DRAWS: Final = 300

STATUS_SCORED: Final = "scored"
STATUS_LOW_COVERAGE: Final = "scored_low_coverage"
STATUS_NOT_COMPUTABLE: Final = "not_computable"
STATUS_ABSTAINED: Final = "abstained"
STATUS_INVALID: Final = "invalid_sample"
STATUS_ANCHOR_ABSENT: Final = "anchor_absent"

#: Draws used to build the combined-score null distribution.
COMBINED_NULL_DRAWS: Final = 4000

#: Confidence bands.
HIGH: Final = "HIGH"
MODERATE: Final = "MODERATE"
LOW: Final = "LOW"
INSUFFICIENT: Final = "INSUFFICIENT"

#: A module needs this fraction of its planned weight measured to be reported.
MIN_MEASURED_WEIGHT_FRACTION: Final = 0.5

#: Below this the reference is too small for a percentile to mean anything,
#: and the engine abstains rather than reporting one.
MIN_REFERENCE_N: Final = 100

#: Above this a percentile estimate is stable enough not to cost confidence.
#: Derived from the cohort bootstrap, not chosen: see `reference_size` below.
STABLE_REFERENCE_N: Final = 200

#: |taxonomic − functional| above this counts as the modules disagreeing.
MODULE_DISAGREEMENT: Final = 0.5


def _probit(p: float) -> float:
    """Inverse standard normal CDF, Acklam's rational approximation.

    Accurate to about 1.15e-9 over the open unit interval, which is far more
    than needed here, and avoids a SciPy dependency for one function.
    """
    if not 0.0 < p < 1.0:
        return 0.0
    a = (-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00)
    b = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00)
    low, high = 0.02425, 1 - 0.02425

    if p < low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
        )
    if p > high:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
        )
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (
        ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1
    )


# --------------------------------------------------------------------------- #
# reference distributions
# --------------------------------------------------------------------------- #


#: Below this prevalence a feature is treated as zero-inflated and scored by
#: the prevalence-aware rank method rather than the MAD-based z.
#:
#: The distinction is not stylistic. For a taxon absent from both the sample
#: and most of the reference, the CLR values still differ — they carry the
#: geometric-mean denominator, which varies with how many other taxa each
#: sample happened to detect. Comparing those values compares sequencing
#: incidentals. Measured on this pipeline it turned five absent colorectal
#: markers into a spurious 74th-percentile similarity score.
#:
#: So absence is treated as a category rather than a number: a sample with the
#: taxon undetected sits at the mid-rank of the reference samples that also
#: lack it, and a sample with it detected ranks above all of them.
ZERO_INFLATED_BELOW: Final = 0.80


def _rank_counts(ordered: Sequence[float], value: float) -> tuple[int, int]:
    """(number strictly below ``value``, number equal to it) in a sorted sequence."""
    below = bisect.bisect_left(ordered, value)
    ties = bisect.bisect_right(ordered, value) - below
    return below, ties


@dataclass(frozen=True, slots=True)
class FeatureReference:
    """Robust location and scale for one feature, plus the raw reference values.

    ``values`` is retained so the null distribution can be rebuilt exactly and
    percentiles taken by ranking rather than by assuming a shape.
    ``n_absent`` records how many reference samples lacked the feature
    entirely, which is what makes the zero-inflated path possible.
    """

    feature: str
    median: float
    mad: float
    n: int
    source: str
    values: tuple[float, ...] = ()
    prevalence: float | None = None
    #: Reference values from samples where the feature *was* detected.
    present_values: tuple[float, ...] = ()
    #: Sorted copies, built once. Ranking a value against the reference is
    #: then a bisection rather than a pass over every reference sample — the
    #: null distribution asks this question once per reference sample per
    #: feature, so it dominated the similarity stage before this.
    _sorted: tuple[float, ...] = field(default=(), init=False, repr=False, compare=False)
    _sorted_present: tuple[float, ...] = field(default=(), init=False, repr=False, compare=False)
    _spread: float | None = field(default=None, init=False, repr=False, compare=False)
    present_set: frozenset[float] = field(default=frozenset(), init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        ordered = tuple(sorted(self.values))
        object.__setattr__(self, "_sorted", ordered)
        object.__setattr__(self, "_sorted_present", tuple(sorted(self.present_values)))
        object.__setattr__(self, "present_set", frozenset(self.present_values))
        spread = None
        if len(ordered) >= 20:
            low = ordered[int(0.10 * (len(ordered) - 1))]
            high = ordered[int(0.90 * (len(ordered) - 1))]
            candidate = (high - low) / 2.563
            spread = candidate if candidate > 0 else None
        object.__setattr__(self, "_spread", spread)

    @property
    def n_absent(self) -> int:
        return self.n - len(self.present_values)

    @property
    def zero_inflated(self) -> bool:
        """True when absence should be treated as a category, not a number.

        Note this holds even when *no* reference sample carried the feature.
        That case is not an edge case to fall through: with a matched stratum
        of 74, four of the five colorectal markers appear in none of them, and
        letting those fall back to a value comparison put an absent taxon at
        the 74th percentile on nothing but reference-frame noise.
        """
        return self.prevalence is not None and self.prevalence < ZERO_INFLATED_BELOW

    def prevalence_aware_percentile(self, value: float, *, detected: bool) -> float | None:
        """Percentile treating "not detected" as one tied category.

        A sample without the feature sits at the mid-rank of the reference
        samples that also lack it. A sample with it ranks above every absent
        reference sample, then by value among the present ones.
        """
        if not self.n:
            return None
        absent = self.n_absent
        if not detected:
            return (absent / 2.0) / self.n * 100.0
        below, ties = _rank_counts(self._sorted_present, value)
        return (absent + below + 0.5 * ties) / self.n * 100.0

    @property
    def scale(self) -> float:
        return 1.4826 * self.mad

    def robust_z(self, value: float, *, detected: bool | None = None) -> float | None:
        """Robust deviation for one measurement, on a z-like scale.

        Three paths, in order of preference:

        1. **Zero-inflated features** (prevalence below `ZERO_INFLATED_BELOW`)
           use the prevalence-aware rank, converted to a normal deviate. The
           best-replicated colorectal cancer taxa are exactly this shape —
           present in a few percent of healthy people — and the MAD-based z is
           not defined for them, because the median and most of the mass sit
           on the zero-replacement floor.
        2. **Well-populated features** use the spec's MAD-based z,
           ``(value − median) / (1.4826 × MAD)``.
        3. **Degenerate spread** falls back to the plain empirical percentile,
           for the case where a feature is prevalent but almost constant.
        """
        if detected is not None and self.zero_inflated:
            percentile = self.prevalence_aware_percentile(value, detected=detected)
            if percentile is not None:
                return _probit(min(max(percentile / 100.0, 0.5 / self.n), 1.0 - 0.5 / self.n))

        if self.scale > 0:
            spread = self._robust_spread()
            # A MAD far smaller than the distribution's overall spread means
            # the centre is a spike rather than a location.
            if spread is None or self.scale >= 0.2 * spread:
                return (value - self.median) / self.scale
        return self._percentile_z(value)

    def _robust_spread(self) -> float | None:
        """Percentile-based scale estimate: (P90 − P10) / 2.563 equals sigma."""
        return self._spread

    def _percentile_z(self, value: float) -> float | None:
        """Empirical percentile mapped through the inverse normal CDF."""
        if not self.values:
            return None
        n = len(self.values)
        below, ties = _rank_counts(self._sorted, value)
        # Mid-rank, then clamped away from 0 and 1 so the inverse is finite.
        fraction = (below + 0.5 * ties) / n
        fraction = min(max(fraction, 0.5 / n), 1.0 - 0.5 / n)
        return _probit(fraction)

    def percentile_of(self, value: float) -> float | None:
        if not self.values:
            return None
        below, ties = _rank_counts(self._sorted, value)
        return (below + 0.5 * ties) / len(self.values) * 100.0

    def to_json(self) -> dict[str, Any]:
        return {
            "feature": self.feature,
            "reference_median": round(self.median, 4),
            "reference_mad": round(self.mad, 4),
            "reference_n": self.n,
            "reference_source": self.source,
            "prevalence": None if self.prevalence is None else round(self.prevalence, 4),
        }


def summarise_reference(
    feature: str,
    values: Sequence[float],
    source: str,
    *,
    prevalence: float | None = None,
    detected: Sequence[bool] | None = None,
) -> FeatureReference:
    """Summarise a reference distribution.

    ``detected`` marks, per reference sample, whether the feature was actually
    observed rather than imputed at the zero floor. Supplying it enables the
    prevalence-aware path; without it a feature is scored as well-populated.
    """
    series = [v for v in values if v is not None and math.isfinite(v)]
    if not series:
        return FeatureReference(feature=feature, median=0.0, mad=0.0, n=0, source=source)
    median = statistics.median(series)
    mad = statistics.median([abs(v - median) for v in series])
    if detected is not None and len(detected) == len(series):
        present = tuple(v for v, seen in zip(series, detected, strict=True) if seen)
    else:
        present = tuple(series)
    return FeatureReference(
        feature=feature,
        median=median,
        mad=mad,
        n=len(series),
        source=source,
        values=tuple(series),
        prevalence=prevalence,
        present_values=present,
    )


# --------------------------------------------------------------------------- #
# per-feature scoring
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class FeatureScore:
    """One feature's contribution, fully decomposed."""

    feature: Feature
    measured: bool
    raw_value: float | None
    transformed_value: float | None
    z: float | None
    v: float | None
    reference: FeatureReference | None
    percentile: float | None
    unmeasured_reason: str = ""

    @property
    def weighted(self) -> float:
        return (self.v or 0.0) * self.feature.factor

    @property
    def missing(self) -> bool:
        """Unmeasured because the data cannot say, as opposed to genuinely absent."""
        return not self.measured and self.unmeasured_reason.startswith("missing")

    def to_json(self) -> dict[str, Any]:
        return {
            **self.feature.to_json(),
            "measured": self.measured,
            "missing": self.missing,
            "raw_value": None if self.raw_value is None else round(self.raw_value, 5),
            "transformed_value": (
                None if self.transformed_value is None else round(self.transformed_value, 4)
            ),
            "robust_z": None if self.z is None else round(self.z, 3),
            "v": None if self.v is None else round(self.v, 4),
            "reference_percentile": None if self.percentile is None else round(self.percentile, 1),
            "weighted_contribution": round(self.weighted, 4),
            "unmeasured_reason": self.unmeasured_reason,
            **({} if self.reference is None else self.reference.to_json()),
        }


def score_feature(
    feature: Feature,
    *,
    transformed_value: float | None,
    raw_value: float | None,
    reference: FeatureReference | None,
    detected: bool | None = None,
    usable_pairs: int | None = None,
) -> FeatureScore:
    """Robust z, direction-signed and clipped onto −1..+1.

    Typed adapter behaviour lives here, keyed on the feature's evidence type:

    * ``TAX_PRES`` — presence/absence. A non-detection at inadequate depth is
      ``missing``, not absence, and does not enter the score.
    * ``TAX_REL`` — zero-aware log-ratio; zero-inflated features use the
      prevalence-aware rank (see `FeatureReference.robust_z`).
    * ``GENE_ABUND`` / ``CARRIER_ABUNDANCE`` / ``ECO`` — continuous values
      against the matching reference distribution.
    """
    if not feature.bound:
        return FeatureScore(
            feature=feature,
            measured=False,
            raw_value=None,
            transformed_value=None,
            z=None,
            v=None,
            reference=None,
            percentile=None,
            unmeasured_reason=f"missing: no engine in this pipeline measures {feature.evidence_type}",
        )
    if (
        feature.evidence_type == "TAX_PRES"
        and detected is False
        and usable_pairs is not None
        and usable_pairs < MIN_PAIRS_FOR_ABSENCE
    ):
        return FeatureScore(
            feature=feature,
            measured=False,
            raw_value=raw_value,
            transformed_value=transformed_value,
            z=None,
            v=None,
            reference=reference,
            percentile=None,
            unmeasured_reason=(
                f"missing: not detected, but only {usable_pairs:,} usable pairs — below the "
                f"{MIN_PAIRS_FOR_ABSENCE:,} at which non-detection can be read as absence"
            ),
        )
    if transformed_value is None:
        return FeatureScore(
            feature=feature,
            measured=False,
            raw_value=raw_value,
            transformed_value=None,
            z=None,
            v=None,
            reference=reference,
            percentile=None,
            unmeasured_reason="not detected in the sample and not present in the reference axis",
        )
    if reference is None or reference.n == 0:
        return FeatureScore(
            feature=feature,
            measured=False,
            raw_value=raw_value,
            transformed_value=transformed_value,
            z=None,
            v=None,
            reference=None,
            percentile=None,
            unmeasured_reason="no reference distribution available for this feature",
        )
    z = reference.robust_z(transformed_value, detected=detected)
    if z is None:
        return FeatureScore(
            feature=feature,
            measured=False,
            raw_value=raw_value,
            transformed_value=transformed_value,
            z=None,
            v=None,
            reference=reference,
            percentile=reference.percentile_of(transformed_value),
            unmeasured_reason=(
                "the reference distribution has zero spread for this feature, so a "
                "z-score is undefined"
            ),
        )
    v = feature.d * max(-Z_CLIP, min(Z_CLIP, z)) / Z_CLIP
    percentile = (
        reference.prevalence_aware_percentile(transformed_value, detected=detected)
        if detected is not None and reference.zero_inflated
        else reference.percentile_of(transformed_value)
    )
    return FeatureScore(
        feature=feature,
        measured=True,
        raw_value=raw_value,
        transformed_value=transformed_value,
        z=z,
        v=v,
        reference=reference,
        percentile=percentile,
    )


# --------------------------------------------------------------------------- #
# module scoring
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ModuleScore:
    name: str
    weight_in_combined: float
    score: float | None
    percentile: float | None
    features: tuple[FeatureScore, ...]
    measured_weight: float
    planned_weight: float
    reference_n: int
    reference_source: str
    module: Module | None = None

    @property
    def n_measured(self) -> int:
        return sum(1 for f in self.features if f.measured)

    @property
    def coverage(self) -> float:
        """Raw feature-weight coverage."""
        return self.measured_weight / self.planned_weight if self.planned_weight else 0.0

    @property
    def cluster_coverage(self) -> float:
        """Coverage over independent evidence clusters (spec 7.4).

        Correlated features share a cluster and count once, so a module that
        measured five species from one lineage cluster is not "five out of
        six" covered.
        """
        expected: dict[str, float] = {}
        observed: dict[str, float] = {}
        for score in self.features:
            key = score.feature.cluster_key
            expected[key] = max(expected.get(key, 0.0), score.feature.factor)
            if score.measured:
                observed[key] = max(observed.get(key, 0.0), score.feature.factor)
        total = sum(expected.values())
        return sum(observed.values()) / total if total else 0.0

    @property
    def n_clusters(self) -> int:
        return len({f.feature.cluster_key for f in self.features})

    @property
    def concordance(self) -> float | None:
        """Pattern concordance on 0–100: ``50 × (1 + Σ(w·r) / Σ(w))`` (spec 7.1)."""
        return None if self.score is None else 50.0 * (1.0 + self.score)

    @property
    def anchor_absent(self) -> bool:
        """Is the organism this module's finding is *about* simply not there?

        A study result of the form "X expanded in cases" is a claim about X.
        The accessory taxa that moved alongside X in that cohort are evidence
        for the pattern only once X is present; on their own they are a
        different sample that happens to share a feature or two.

        So when every declared anchor was looked for and not detected, the
        module reports `anchor_absent` rather than a resemblance score. The
        features keep their measurements and stay in the report — only the
        claim is withdrawn.
        """
        if self.module is None or not self.module.anchor_features:
            return False
        anchors = [
            f for f in self.features if f.feature.name in set(self.module.anchor_features)
        ]
        if not anchors:
            return False
        # `raw_value is None` after a successful measurement is the profiler
        # saying it looked and found nothing, as distinct from `measured`
        # being False, which means it could not look.
        return all(f.measured and f.raw_value is None for f in anchors)

    @property
    def reportable(self) -> bool:
        return (
            self.score is not None
            and self.coverage >= MIN_MEASURED_WEIGHT_FRACTION
            and not self.anchor_absent
        )

    @property
    def status(self) -> str:
        if self.module is not None and not self.module.bound:
            return STATUS_NOT_COMPUTABLE
        if self.score is None:
            return STATUS_NOT_COMPUTABLE
        if self.anchor_absent:
            return STATUS_ANCHOR_ABSENT
        if self.cluster_coverage < LOW_COVERAGE_BELOW:
            return STATUS_LOW_COVERAGE
        return STATUS_SCORED

    @property
    def anchor_note(self) -> str:
        """Why the module withdrew its claim, in words, for the report."""
        if not self.anchor_absent or self.module is None:
            return ""
        names = ", ".join(
            n.replace("_", " ") for n in self.module.anchor_features
        )
        return (
            f"Not scored: {names} — the organism this finding is about — was looked for and "
            "not detected. The remaining features of this module describe taxa that moved "
            "alongside it in the source cohort; without it they do not constitute the pattern."
        )

    @property
    def n_missing(self) -> int:
        return sum(1 for f in self.features if f.missing)

    def to_json(self) -> dict[str, Any]:
        module_meta = {} if self.module is None else {
            "type": self.module.type,
            "evidence_type": self.module.evidence_type,
            "claim_level": self.module.claim_level,
            "assay_transport": self.module.assay_transport,
            "study_group": self.module.group,
            "module_status": self.module.status,
            "maturity": self.module.maturity,
            "bound": self.module.bound,
            "would_bind_if": self.module.would_bind_if,
            "fuses": self.module.fuses,
            "label": self.module.label or self.module.name,
        }
        return {
            "module": self.name,
            **module_meta,
            "status": self.status,
            "weight_in_combined": self.weight_in_combined,
            "score": None if self.score is None else round(self.score, 4),
            "pattern_concordance_percent": (
                None if self.concordance is None else round(self.concordance, 1)
            ),
            "percentile": None if self.percentile is None else round(self.percentile, 1),
            "features_measured": self.n_measured,
            "features_missing": self.n_missing,
            "features_planned": len(self.features),
            "raw_feature_coverage": round(self.coverage, 3),
            "independent_cluster_coverage": round(self.cluster_coverage, 3),
            "n_clusters": self.n_clusters,
            "weight_coverage": round(self.coverage, 3),
            "reportable": self.reportable,
            "anchor_features": list(
                self.module.anchor_features if self.module is not None else ()
            ),
            "anchor_absent": self.anchor_absent,
            "anchor_note": self.anchor_note,
            "reference_n": self.reference_n,
            "reference_source": self.reference_source,
            "features": [f.to_json() for f in self.features],
        }


def aggregate_module(scores: Sequence[FeatureScore]) -> float | None:
    """Σ(w q s c v) / Σ(w q s c) over measured features only."""
    numerator = 0.0
    denominator = 0.0
    for score in scores:
        if not score.measured or score.v is None:
            continue
        weight = score.feature.factor
        numerator += weight * score.v
        denominator += weight
    return numerator / denominator if denominator > 0 else None


# --------------------------------------------------------------------------- #
# sample measurements
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class SampleMeasurements:
    """Everything measured on one sample, keyed by engine."""

    #: species -> CLR value, on the reference cohort's taxon axis
    taxonomic_clr: dict[str, float] = field(default_factory=dict)
    #: species -> relative abundance percent, retained for the human report
    taxonomic_abundance: dict[str, float] = field(default_factory=dict)
    #: genus -> aggregate log-ratio (sum of member species), for transported
    #: 16S features; and genus -> summed relative abundance percent
    genus_clr: dict[str, float] = field(default_factory=dict)
    genus_abundance: dict[str, float] = field(default_factory=dict)
    #: panel name -> aggregate log-ratio of the panel's curated carrier
    #: species, and the summed relative abundance percent behind it
    carrier_clr: dict[str, float] = field(default_factory=dict)
    carrier_abundance: dict[str, float] = field(default_factory=dict)
    #: panel or gene name -> copies per 100 genomes
    functional: dict[str, float] = field(default_factory=dict)
    #: shannon_diversity, richness, ...
    ecological: dict[str, float] = field(default_factory=dict)
    #: qc and metadata used by the abstention rules
    usable_nonhost_reads: int | None = None
    antibiotics_days_ago: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    #: Metadata keys whose value in ``metadata`` was assumed rather than
    #: supplied (see :mod:`openbiota.context`). A rule that would abstain on a
    #: missing key scores instead and records the assumption.
    assumed_keys: frozenset[str] = frozenset()

    @property
    def usable_pairs(self) -> int | None:
        return None if self.usable_nonhost_reads is None else self.usable_nonhost_reads // 2

    def value_for(self, feature: Feature) -> tuple[float | None, float | None]:
        """Return ``(transformed, raw)`` for a feature, or ``(None, None)``."""
        key = feature.key
        if feature.engine == "metaphlan":
            if feature.level == "genus":
                return self.genus_clr.get(key), self.genus_abundance.get(key)
            transformed = self.taxonomic_clr.get(key)
            raw = self.taxonomic_abundance.get(key)
            return transformed, raw
        if feature.engine == "carrier":
            return self.carrier_clr.get(key), self.carrier_abundance.get(key)
        if feature.engine == "diamond":
            raw = self.functional.get(key)
            return raw, raw
        if feature.engine == "unbound":
            return None, None
        raw = self.ecological.get(key)
        return raw, raw

    def detected(self, feature: Feature) -> bool | None:
        """Whether the feature was observed, as opposed to imputed at the floor.

        Only meaningful for taxonomic and carrier features; gene-capacity and
        ecological values are continuous and always "observed".
        """
        if feature.engine == "metaphlan":
            if feature.level == "genus":
                return self.genus_abundance.get(feature.key, 0.0) > 0
            return self.taxonomic_abundance.get(feature.key, 0.0) > 0
        if feature.engine == "carrier":
            return self.carrier_abundance.get(feature.key, 0.0) > 0
        return None


@dataclass(slots=True)
class ReferenceBundle:
    """Reference distributions, per engine, with provenance."""

    taxonomic: dict[str, FeatureReference] = field(default_factory=dict)
    genus: dict[str, FeatureReference] = field(default_factory=dict)
    carrier: dict[str, FeatureReference] = field(default_factory=dict)
    functional: dict[str, FeatureReference] = field(default_factory=dict)
    ecological: dict[str, FeatureReference] = field(default_factory=dict)
    taxonomic_source: str = "unavailable"
    functional_source: str = "unavailable"
    ecological_source: str = "unavailable"
    taxonomic_n: int = 0
    functional_n: int = 0
    ecological_n: int = 0
    #: Size of the *unmatched* cohort. Abstention keys off this, because a
    #: small matched stratum is a reason to lower confidence, not to refuse:
    #: the spec's own matching threshold is 30, and refusing everything below
    #: 100 would make matching self-defeating.
    taxonomic_full_n: int = 0
    match_info: dict[str, Any] = field(default_factory=dict)
    profiler_matches_reference: bool | None = None

    def for_feature(self, feature: Feature) -> FeatureReference | None:
        if feature.engine == "metaphlan":
            if feature.level == "genus":
                return self.genus.get(feature.key)
            return self.taxonomic.get(feature.key)
        if feature.engine == "carrier":
            return self.carrier.get(feature.key)
        if feature.engine == "diamond":
            return self.functional.get(feature.key)
        if feature.engine == "unbound":
            return None
        return self.ecological.get(feature.key)

    def source_for(self, module: str, engine: str | None = None) -> tuple[str, int]:
        engine = engine or {
            MODULE_TAXONOMIC: "metaphlan",
            MODULE_FUNCTIONAL: "diamond",
        }.get(module, "ecological")
        if engine in ("metaphlan", "carrier"):
            return self.taxonomic_source, self.taxonomic_n
        if engine == "diamond":
            return self.functional_source, self.functional_n
        if engine == "unbound":
            return "no engine bound", 0
        return self.ecological_source, self.ecological_n

    def per_sample_values(self, feature: Feature) -> tuple[float, ...]:
        reference = self.for_feature(feature)
        return reference.values if reference else ()


# --------------------------------------------------------------------------- #
# null distributions and percentiles
# --------------------------------------------------------------------------- #


def module_null(
    features: Sequence[Feature], references: ReferenceBundle
) -> list[float]:
    """Score every reference sample through the identical path.

    Requires the features in a module to share a reference sample axis, which
    holds by construction: a module's features all use one engine.
    """
    usable = [f for f in features if references.for_feature(f) is not None]
    if not usable:
        return []
    lengths = {len(references.per_sample_values(f)) for f in usable}
    lengths.discard(0)
    if len(lengths) != 1:
        # Features measured on different sample sets cannot form a joint null.
        return []
    n = lengths.pop()

    numerators = [0.0] * n
    denominators = [0.0] * n
    for feature in usable:
        reference = references.for_feature(feature)
        if reference is None or not reference.values:
            continue
        present = reference.present_set if reference.zero_inflated else None
        factor = feature.factor
        scale = feature.d / Z_CLIP
        for index, value in enumerate(reference.values):
            # A reference sample counts as "detected" when its value is among
            # the observed ones, which keeps the null on the same footing as
            # the sample being scored.
            seen = (value in present) if present is not None else None
            z = reference.robust_z(value, detected=seen)
            if z is None:
                continue
            numerators[index] += factor * scale * max(-Z_CLIP, min(Z_CLIP, z))
            denominators[index] += factor
    return [num / den for num, den in zip(numerators, denominators, strict=True) if den > 0]


def percentile_against(null: Sequence[float], value: float | None) -> float | None:
    if value is None or not null:
        return None
    below = sum(1 for v in null if v < value)
    ties = sum(1 for v in null if v == value)
    return (below + 0.5 * ties) / len(null) * 100.0


def combined_null(
    module_nulls: Mapping[str, Sequence[float]],
    weights: Mapping[str, float],
    *,
    draws: int = COMBINED_NULL_DRAWS,
    seed: int = 20260904,
) -> list[float]:
    """Null for the weighted combination of modules.

    The modules do not share a reference sample axis — taxonomic features come
    from curatedMetagenomicData, functional ones from the DIAMOND cohort — so
    the joint null is built by resampling each module independently. That
    assumes the modules are independent across reference samples, which is
    stated rather than assumed silently: it makes the combined null slightly
    wider than reality if the modules are positively correlated, which is the
    conservative direction for a percentile.
    """
    usable = {k: list(v) for k, v in module_nulls.items() if v and weights.get(k, 0) > 0}
    if not usable:
        return []
    rng = random.Random(seed)
    total_weight = sum(weights[k] for k in usable)
    if total_weight <= 0:
        return []

    out: list[float] = []
    for _ in range(draws):
        value = 0.0
        for name, series in usable.items():
            value += weights[name] * series[rng.randrange(len(series))]
        out.append(value / total_weight)
    return out


# --------------------------------------------------------------------------- #
# cross-engine concordance
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class CrossEngineResult:
    taxonomic_feature: str
    functional_feature: str
    taxonomic_v: float | None
    functional_v: float | None
    taxonomic_detail: str
    functional_detail: str
    verdict: str
    message: str

    @property
    def concordant(self) -> bool:
        return self.verdict == "CONCORDANT"

    def to_json(self) -> dict[str, Any]:
        return asdict(self) | {"concordant": self.concordant}


def evaluate_cross_engine(
    profile: Profile, module_scores: Mapping[str, ModuleScore]
) -> list[CrossEngineResult]:
    """Compare two independent measurements of the same biology."""
    by_name: dict[str, FeatureScore] = {
        f.feature.name: f for m in module_scores.values() for f in m.features
    }
    out: list[CrossEngineResult] = []

    for check in profile.cross_engine_checks:
        tax = by_name.get(check.taxonomic_feature)
        fun = by_name.get(check.functional_feature)
        if tax is None or fun is None:
            continue
        tax_detail = (
            f"{tax.feature.name}: {tax.z:+.1f} SD vs reference "
            f"({'depleted' if tax.z < 0 else 'enriched'}, "
            f"{'as' if tax.v > 0 else 'against'} the published direction)"
            if tax.z is not None and tax.v is not None
            else f"{tax.feature.name}: not measured"
        )
        fun_detail = (
            f"{fun.feature.name}: {fun.raw_value:.1f} copies/100 genomes"
            + (f", {fun.percentile:.0f}th pct" if fun.percentile is not None else "")
            if fun.raw_value is not None
            else f"{fun.feature.name}: not measured"
        )

        if tax.v is None or fun.v is None:
            verdict, message = (
                "UNAVAILABLE",
                "One side of the comparison was not measured, so the engines cannot be "
                "cross-checked.",
            )
        elif tax.v * fun.v >= 0 or abs(tax.v - fun.v) <= MODULE_DISAGREEMENT:
            verdict = "CONCORDANT"
            message = (
                f"Both engines agree on {check.functional_feature}. Two independent "
                "measurements of the same biology point the same way."
            )
        else:
            # v is direction-signed: positive means concordant with the
            # direction the profile declares, whichever way that points.
            verdict = "DISCORDANT"
            tax_says = "matches" if tax.v > 0 else "contradicts"
            fun_says = "matches" if fun.v > 0 else "contradicts"
            message = (
                f"The two engines disagree about {check.functional_feature}. The species "
                f"measurement {tax_says} the published pattern while the gene-capacity "
                f"measurement {fun_says} it. These are independent measurements of the same "
                "biology, so at least one is misleading here — confidence is reduced and the "
                "combined score should not be read on its own."
            )
        out.append(
            CrossEngineResult(
                taxonomic_feature=check.taxonomic_feature,
                functional_feature=check.functional_feature,
                taxonomic_v=None if tax.v is None else round(tax.v, 3),
                functional_v=None if fun.v is None else round(fun.v, 3),
                taxonomic_detail=tax_detail,
                functional_detail=fun_detail,
                verdict=verdict,
                message=message,
            )
        )
    return out


# --------------------------------------------------------------------------- #
# confidence and abstention
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ConfidenceComponent:
    name: str
    passed: bool
    detail: str
    penalty: float = 0.0

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class Confidence:
    grade: str
    components: tuple[ConfidenceComponent, ...]

    @property
    def reasons(self) -> list[str]:
        return [c.detail for c in self.components if not c.passed]

    def to_json(self) -> dict[str, Any]:
        return {
            "grade": self.grade,
            "reasons_for_downgrade": self.reasons,
            "components": [c.to_json() for c in self.components],
        }


def grade_confidence(
    *,
    profile: Profile,
    module_scores: Mapping[str, ModuleScore],
    references: ReferenceBundle,
    measurements: SampleMeasurements,
    cross_engine: Sequence[CrossEngineResult],
) -> Confidence:
    """Grade on coverage, reference quality, versions, depth, metadata, agreement."""
    components: list[ConfidenceComponent] = []

    coverage = [
        m.coverage for m in module_scores.values()
        if (m.module.fuses if m.module is not None else m.name in COMBINED_MODULES)
    ]
    mean_coverage = statistics.fmean(coverage) if coverage else 0.0
    components.append(
        ConfidenceComponent(
            "feature_coverage",
            mean_coverage >= 0.75,
            f"{mean_coverage:.0%} of planned feature weight was measurable"
            + ("" if mean_coverage >= 0.75 else " (below 75%)"),
            penalty=0.0 if mean_coverage >= 0.75 else (1.0 if mean_coverage < 0.5 else 0.5),
        )
    )

    # Bootstrapping the cohort showed the spread of the median estimate falling
    # from 0.81 at n=25 to 0.32 at n=100 and 0.10 at n=800, so a few hundred
    # samples is where a tail percentile becomes stable.
    tax_n = references.taxonomic_n
    stable = tax_n >= STABLE_REFERENCE_N
    components.append(
        ConfidenceComponent(
            "reference_size",
            stable,
            f"comparison group holds {tax_n:,} samples"
            + (
                ""
                if stable
                else f" — below {STABLE_REFERENCE_N}, where bootstrapping shows percentile "
                "estimates are still noticeably unstable"
            ),
            penalty=0.0 if stable else (1.0 if tax_n < MIN_REFERENCE_N else 0.5),
        )
    )

    matched = bool(references.match_info.get("matched"))
    components.append(
        ConfidenceComponent(
            "reference_matching",
            matched,
            str(references.match_info.get("reason", "no matching information")),
            penalty=0.0 if matched else 0.5,
        )
    )

    versions_ok = references.profiler_matches_reference is not False
    components.append(
        ConfidenceComponent(
            "profiler_version_match",
            versions_ok,
            (
                "the sample and the reference cohort were profiled with the same tool family"
                if versions_ok
                else "the sample and the reference cohort were profiled with different tool "
                "versions, so species abundances are only approximately comparable"
            ),
            penalty=0.0 if versions_ok else 0.5,
        )
    )

    reads = measurements.usable_nonhost_reads
    depth_ok = reads is None or reads >= 5_000_000
    components.append(
        ConfidenceComponent(
            "usable_depth",
            depth_ok,
            (
                f"{reads:,} usable non-host reads" if reads is not None
                else "usable non-host read count unknown"
            )
            + ("" if depth_ok else " (thin for species-level profiling)"),
            penalty=0.0 if depth_ok else 0.5,
        )
    )

    missing = [
        key
        for key in ("antibiotics", "stool_form", "onset_date", "sampling_date", "medications")
        if not measurements.metadata.get(key)
    ]
    components.append(
        ConfidenceComponent(
            "metadata_completeness",
            not missing,
            (
                "all confounder metadata present"
                if not missing
                else f"missing metadata: {', '.join(missing)} — confounders cannot be excluded"
            ),
            penalty=0.0 if not missing else min(1.0, 0.25 * len(missing)),
        )
    )

    primary_ok = profile.references.has_primary
    components.append(
        ConfidenceComponent(
            "exposure_matched_contrast",
            primary_ok,
            (
                "an exposure-matched contrast group is available"
                if primary_ok
                else "no exposure-matched contrast exists for this profile, so "
                "syndrome-specific signal cannot be separated from post-exposure signal"
            ),
            penalty=0.0 if primary_ok else 0.5,
        )
    )

    discordant = [c for c in cross_engine if c.verdict == "DISCORDANT"]
    components.append(
        ConfidenceComponent(
            "cross_engine_agreement",
            not discordant,
            (
                "taxonomic and functional engines agree"
                if not discordant
                else (
                    "a cross-engine check disagrees: two independent measurements of the "
                    "same biology point different ways"
                    if len(discordant) == 1
                    else f"{len(discordant)} cross-engine checks disagree"
                )
            ),
            penalty=0.0 if not discordant else 1.0,
        )
    )

    tax = module_scores.get(MODULE_TAXONOMIC) or _first_fused(module_scores, "metaphlan")
    fun = module_scores.get(MODULE_FUNCTIONAL) or _first_fused(module_scores, "diamond")
    modules_agree = True
    if tax and fun and tax.score is not None and fun.score is not None:
        modules_agree = abs(tax.score - fun.score) <= MODULE_DISAGREEMENT or (
            tax.score * fun.score >= 0
        )
    components.append(
        ConfidenceComponent(
            "module_agreement",
            modules_agree,
            (
                "taxonomic and functional modules point the same way"
                if modules_agree
                else "taxonomic and functional modules point opposite ways"
            ),
            penalty=0.0 if modules_agree else 1.0,
        )
    )

    penalty = sum(c.penalty for c in components)
    if any(not c.passed and c.penalty >= 1.0 for c in components) and penalty >= 2.0:
        grade = INSUFFICIENT
    elif penalty >= 2.0:
        grade = LOW
    elif penalty >= 0.75:
        grade = MODERATE
    else:
        grade = HIGH
    return Confidence(grade=grade, components=tuple(components))


def _first_fused(module_scores: Mapping[str, ModuleScore], engine: str) -> ModuleScore | None:
    for m in module_scores.values():
        if m.module is not None and m.module.fuses and m.module.engine == engine:
            return m
    return None


@dataclass(frozen=True, slots=True)
class Abstention:
    """A refusal to score, with what would need to change.

    ``assumptions`` records the facts this profile depends on that were not
    supplied and were assumed instead of triggering an abstention. A profile
    scored on an assumption is a real score, printed with the assumption
    beside it, never a silent default.
    """

    abstained: bool
    triggered: tuple[str, ...] = ()
    remedies: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "abstained": self.abstained,
            "reasons": list(self.triggered),
            "what_would_change_this": list(self.remedies),
            "scored_on_assumptions": list(self.assumptions),
        }


def check_abstention(
    profile: Profile, measurements: SampleMeasurements, references: ReferenceBundle
) -> Abstention:
    """Abstention is a feature. A refused score beats a confident wrong one."""
    triggered: list[str] = []
    remedies: list[str] = []
    assumptions: list[str] = []

    def _absent(keys: Sequence[str]) -> tuple[list[str], list[str]]:
        """Split required keys into truly absent and assumed-absent."""
        missing, assumed = [], []
        for k in keys:
            if k in measurements.assumed_keys:
                assumed.append(k)
            elif measurements.metadata.get(k) in (None, "", False):
                missing.append(k)
        return missing, assumed

    for rule in profile.abstain_if:
        if rule.kind == "recent_antibiotics_days_lt":
            days = measurements.antibiotics_days_ago
            if days is not None and days < int(rule.value):
                triggered.append(rule.reason)
                remedies.append(
                    f"re-sample at least {rule.value} days after finishing antibiotics"
                )
        elif rule.kind == "usable_nonhost_reads_lt":
            reads = measurements.usable_nonhost_reads
            if reads is not None and reads < int(rule.value):
                triggered.append(rule.reason)
                remedies.append(f"sequence to at least {int(rule.value):,} usable non-host reads")
        elif rule.kind == "usable_nonhost_pairs_lt":
            pairs = measurements.usable_pairs
            if pairs is not None and pairs < int(rule.value):
                triggered.append(rule.reason)
                remedies.append(f"sequence to at least {int(rule.value):,} usable non-host read pairs")
        elif rule.kind == "missing_metadata":
            required = rule.value if isinstance(rule.value, list) else [rule.value]
            absent, assumed = _absent(required)
            if absent:
                triggered.append(f"{rule.reason} (missing: {', '.join(absent)})")
                remedies.append(f"supply {', '.join(absent)}")
            assumptions.extend(assumed)

    # Mandatory matching variables (spec 4.3 / 12): a profile whose signal is
    # sex- or age-specific cannot be placed without them, and a drug that
    # produces the published signature must be known before the signature is
    # read as disease. A drug that was *assumed* absent (research and
    # clinician modes, nothing supplied) does not block the score: the profile
    # is placed on the untreated assumption and says so.
    absent, assumed = _absent(profile.mandatory_match)
    if absent:
        triggered.append(
            f"this profile requires {', '.join(absent)} before it can be placed; the published "
            "signal is specific to it"
        )
        remedies.append(f"supply {', '.join(absent)} on the command line")
    assumptions.extend(assumed)

    full_n = references.taxonomic_full_n or references.taxonomic_n
    if full_n and full_n < MIN_REFERENCE_N:
        triggered.append(
            f"the taxonomic reference cohort holds only {full_n} samples, below the "
            f"{MIN_REFERENCE_N} needed for a percentile to be meaningful"
        )
        remedies.append("rebuild the reference cohort with `openbiota taxonomic-cohort`")

    return Abstention(
        abstained=bool(triggered),
        triggered=tuple(dict.fromkeys(triggered)),
        remedies=tuple(dict.fromkeys(remedies)),
        assumptions=tuple(dict.fromkeys(assumptions)),
    )


# --------------------------------------------------------------------------- #
# confidence as a vector (spec 7.5)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ConfidenceVector:
    """Confidence reported beside the score, never multiplied into it.

    Each dimension answers a different question. Folding them into one number
    would hide which one is the problem, and multiplying them into the score
    would make weak evidence read as "typical" instead of "unknown".
    """

    technical_reliability: float
    independent_feature_coverage: float
    evidence_maturity: str
    assay_transportability: str
    population_transportability: float | None
    empirical_specificity: float | None
    literature_uniqueness_prior: float | None
    out_of_distribution: bool | None
    direction_conflict: str
    interval_width: float | None
    interval: tuple[float, float] | None = None
    notes: tuple[str, ...] = ()

    @property
    def maturity_meaning(self) -> str:
        from openbiota.profiles import MATURITY_MEANING

        return MATURITY_MEANING.get(self.evidence_maturity, "")

    def to_json(self) -> dict[str, Any]:
        return {
            "technical_reliability": round(self.technical_reliability, 2),
            "independent_feature_coverage": round(self.independent_feature_coverage, 2),
            "evidence_maturity": self.evidence_maturity,
            "evidence_maturity_meaning": self.maturity_meaning,
            "assay_transportability": self.assay_transportability,
            "population_transportability": (
                "unknown" if self.population_transportability is None
                else round(self.population_transportability, 2)
            ),
            "empirical_specificity": (
                "unknown" if self.empirical_specificity is None
                else round(self.empirical_specificity, 2)
            ),
            "literature_uniqueness_prior": (
                None if self.literature_uniqueness_prior is None
                else round(self.literature_uniqueness_prior, 2)
            ),
            "out_of_distribution": (
                "unknown" if self.out_of_distribution is None else self.out_of_distribution
            ),
            "direction_conflict": self.direction_conflict,
            "interval_width": None if self.interval_width is None else round(self.interval_width, 1),
            "uncertainty_interval": (
                None if self.interval is None else [round(self.interval[0], 1), round(self.interval[1], 1)]
            ),
            "notes": list(self.notes),
        }


def bootstrap_interval(
    *,
    module_scores: Mapping[str, ModuleScore],
    weights: Mapping[str, float],
    nulls: Mapping[str, Sequence[float]],
    draws: int = INTERVAL_DRAWS,
    seed: int = 20260905,
) -> tuple[float, float] | None:
    """Percentile interval from resampling the measured features.

    Within each fused module the measured features are resampled with
    replacement and the module re-aggregated; the modules are then fused with
    the same weights and the fused value mapped through the same null. The
    5th–95th percentiles of that distribution are the interval. A module
    resting on two features produces a wide interval by construction, which
    is the behaviour the spec asks for: a score built on little evidence
    returns, but returns visibly uncertain.
    """
    usable = {
        name: [f for f in module_scores[name].features if f.measured and f.v is not None]
        for name in weights
        if name in module_scores and module_scores[name].score is not None
    }
    usable = {k: v for k, v in usable.items() if v}
    if not usable:
        return None
    combined = combined_null({k: v for k, v in nulls.items() if k in usable}, weights)
    if not combined:
        return None
    rng = random.Random(seed)
    total = sum(weights[k] for k in usable)
    values: list[float] = []
    for _ in range(draws):
        fused = 0.0
        for name, scores in usable.items():
            picks = [scores[rng.randrange(len(scores))] for _ in scores]
            numerator = sum(p.feature.factor * (p.v or 0.0) for p in picks)
            denominator = sum(p.feature.factor for p in picks)
            fused += weights[name] * (numerator / denominator if denominator else 0.0)
        values.append(fused / total)
    percentiles = sorted(
        p for p in (percentile_against(combined, v) for v in values) if p is not None
    )
    if not percentiles:
        return None
    low = percentiles[int(0.05 * (len(percentiles) - 1))]
    high = percentiles[int(0.95 * (len(percentiles) - 1))]
    return low, high


def build_confidence_vector(
    *,
    profile: Profile,
    module_scores: Mapping[str, ModuleScore],
    weights: Mapping[str, float],
    references: ReferenceBundle,
    measurements: SampleMeasurements,
    cross_engine: Sequence[CrossEngineResult],
    interval: tuple[float, float] | None,
    empirical_specificity: float | None,
    literature_uniqueness_prior: float | None,
) -> ConfidenceVector:
    notes: list[str] = []

    # Technical reliability: starts at 1 and loses ground for each thing
    # that makes the measurement itself less trustworthy.
    reliability = 1.0
    discordant = [c for c in cross_engine if c.verdict == "DISCORDANT"]
    if discordant:
        reliability -= 0.3
        notes.append("independent measurements of the same biology disagree")
    if references.profiler_matches_reference is False:
        reliability -= 0.3
        notes.append("sample and reference profiled with different tool versions")
    reads = measurements.usable_nonhost_reads
    if reads is not None and reads < 5_000_000:
        reliability -= 0.2
        notes.append("thin sequencing depth for species-level profiling")
    if measurements.metadata.get("ppi"):
        reliability -= 0.1
        notes.append("proton-pump inhibitor use recorded; large documented community shift")
    if measurements.metadata.get("hospitalised"):
        reliability -= 0.1
        notes.append("hospitalisation recorded; drives E. coli, Klebsiella and resistance genes")
    reliability = max(0.0, min(1.0, reliability))

    fused = [module_scores[n] for n in weights if n in module_scores]
    if fused:
        total = sum(weights[m.name] for m in fused)
        coverage = sum(weights[m.name] * m.cluster_coverage for m in fused) / total if total else 0.0
    else:
        coverage = 0.0

    transported = any(
        m.module is not None and m.module.assay_transport != "native" for m in fused
    )
    transport = "labeled_transport" if transported else "native"

    matched = references.match_info.get("matched")
    if not references.match_info:
        population: float | None = None
    elif matched:
        population = 1.0
    else:
        population = 0.6
        notes.append("no matched reference stratum; full cohort used")

    ood: bool | None = None
    age = measurements.metadata.get("age")
    pediatric = any(k in profile.population.lower() for k in ("pediatric", "child", "infant"))
    if pediatric:
        if age is None:
            ood = None
        elif float(age) >= 18:
            ood = True
            notes.append("profile derived from children; adult sample is out of distribution")
        else:
            ood = False
    elif references.match_info:
        ood = False

    conflict = "none"
    conflicts = [f for f in profile.features() if f.direction_conflict != "none"]
    if any(f.direction_conflict == "material" for f in conflicts):
        conflict = "material"
    elif conflicts:
        conflict = "low"

    return ConfidenceVector(
        technical_reliability=reliability,
        independent_feature_coverage=coverage,
        evidence_maturity=profile.evidence_maturity,
        assay_transportability=transport,
        population_transportability=population,
        empirical_specificity=empirical_specificity,
        literature_uniqueness_prior=literature_uniqueness_prior,
        out_of_distribution=ood,
        direction_conflict=conflict,
        interval_width=None if interval is None else interval[1] - interval[0],
        interval=interval,
        notes=tuple(notes),
    )


# --------------------------------------------------------------------------- #
# the whole result
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ProfileResult:
    profile: Profile
    modules: dict[str, ModuleScore]
    combined_score: float | None
    combined_percentile: float | None
    confidence: Confidence
    abstention: Abstention
    cross_engine: tuple[CrossEngineResult, ...]
    stratum: str | None = None
    dysbiosis: Any | None = None
    reference_manifest: dict[str, Any] = field(default_factory=dict)
    vector: ConfidenceVector | None = None
    weights: dict[str, float] = field(default_factory=dict)
    #: Filled in by the shape analysis once every profile has scored.
    competing: tuple[tuple[str, float], ...] = ()
    top_two_margin: float | None = None

    @property
    def reportable(self) -> bool:
        return not self.abstention.abstained and self.combined_percentile is not None

    @property
    def concordance(self) -> float | None:
        """``pattern_concordance_percent`` — resemblance on 0–100, not risk."""
        return None if self.combined_score is None else 50.0 * (1.0 + self.combined_score)

    @property
    def status(self) -> str:
        if self.abstention.abstained:
            return STATUS_ABSTAINED
        if self.combined_percentile is None:
            return STATUS_NOT_COMPUTABLE
        fused = [self.modules[n] for n in self.weights if n in self.modules]
        if fused and self.vector is not None and self.vector.independent_feature_coverage < LOW_COVERAGE_BELOW:
            return STATUS_LOW_COVERAGE
        return STATUS_SCORED

    def module(self, name: str) -> ModuleScore | None:
        return self.modules.get(name)

    def primary_module(self, engine: str) -> ModuleScore | None:
        """The highest-weighted fused module measured by ``engine``.

        Legacy profiles have exactly one module per engine; evidence-typed
        profiles may have several study modules, and the report wants the one
        that carries the most weight.
        """
        candidates = [
            m for m in self.modules.values()
            if (m.module.engine if m.module is not None else _legacy_engine(m.name)) == engine
        ]
        if not candidates:
            return None
        candidates.sort(key=lambda m: -self.weights.get(m.name, m.weight_in_combined))
        return candidates[0]

    @property
    def fused_modules(self) -> list[ModuleScore]:
        return [self.modules[n] for n in self.weights if n in self.modules]

    @property
    def unbound_modules(self) -> list[ModuleScore]:
        return [m for m in self.modules.values() if m.module is not None and not m.module.bound]

    def to_json(self) -> dict[str, Any]:
        return {
            "profile": self.profile.to_json(),
            "status": self.status,
            "stratum": self.stratum,
            "modules": {k: v.to_json() for k, v in self.modules.items()},
            "combined": {
                "score": None if self.combined_score is None else round(self.combined_score, 4),
                "pattern_concordance_percent": (
                    None if self.concordance is None else round(self.concordance, 1)
                ),
                "percentile": (
                    None
                    if self.combined_percentile is None
                    else round(self.combined_percentile, 1)
                ),
                "control_tail_percentile": (
                    None
                    if self.combined_percentile is None
                    else round(self.combined_percentile, 1)
                ),
                "weights": {k: round(v, 4) for k, v in self.weights.items()},
                "fusion": self.profile.fusion,
                "note": " ".join(self.profile.fusion_note.split()),
                "excluded_modules": [
                    n for n, m in self.modules.items() if n not in self.weights
                ],
                "calibration": "uncalibrated_pattern_concordance",
            },
            "confidence": self.confidence.to_json(),
            "confidence_vector": None if self.vector is None else self.vector.to_json(),
            "competing_profiles": [
                {"profile": name, "percentile": round(p, 1)} for name, p in self.competing
            ],
            "top_two_margin": None if self.top_two_margin is None else round(self.top_two_margin, 1),
            "abstention": self.abstention.to_json(),
            "cross_engine_checks": [c.to_json() for c in self.cross_engine],
            "dysbiosis_anchor": None if self.dysbiosis is None else self.dysbiosis.to_json(),
            "reference_manifest": self.reference_manifest,
            "interpretation_note": (
                "A percentile here means the sample is more concordant with this profile's "
                "prespecified pattern than that fraction of the matched reference set. "
                "Pattern concordance is resemblance to a published group-level pattern on the "
                "profile's frozen reference scale. Neither is a probability of having the "
                "condition."
            ),
        }


def _legacy_engine(module_name: str) -> str:
    return {
        MODULE_TAXONOMIC: "metaphlan",
        MODULE_FUNCTIONAL: "diamond",
    }.get(module_name, "ecological")


def _features_for_stratum(module: Module, stratum: str | None) -> tuple[Feature, ...]:
    """Features that apply to the selected stratum (or to every stratum)."""
    return tuple(f for f in module.features if f.stratum is None or f.stratum == stratum)


def score_profile(
    *,
    profile: Profile,
    measurements: SampleMeasurements,
    references: ReferenceBundle,
    dysbiosis: Any | None = None,
    stratum: str | None = None,
    reference_manifest: Mapping[str, Any] | None = None,
    empirical_specificity: float | None = None,
    literature_uniqueness_prior: float | None = None,
) -> ProfileResult:
    """Score one profile against one sample."""
    module_scores: dict[str, ModuleScore] = {}
    nulls: dict[str, list[float]] = {}

    for name, module in profile.modules.items():
        features = _features_for_stratum(module, stratum)
        if module.features and not features and stratum is not None:
            # Every feature belongs to another stratum; the module has nothing
            # to say about this one and must not be scored on someone else's
            # pattern.
            features = ()
        scores: list[FeatureScore] = []
        for feature in features:
            transformed, raw = measurements.value_for(feature)
            scores.append(
                score_feature(
                    feature,
                    transformed_value=transformed,
                    raw_value=raw,
                    reference=references.for_feature(feature),
                    detected=measurements.detected(feature),
                    usable_pairs=measurements.usable_pairs,
                )
            )
        score = aggregate_module(scores)
        null = module_null(features, references) if module.bound else []
        nulls[name] = null
        source, n = references.source_for(name, module.engine)
        module_scores[name] = ModuleScore(
            name=name,
            weight_in_combined=module.weight_in_combined,
            score=score,
            percentile=percentile_against(null, score),
            features=tuple(scores),
            measured_weight=sum(f.feature.factor for f in scores if f.measured),
            planned_weight=sum(f.feature.factor for f in scores),
            reference_n=n,
            reference_source=source,
            module=module,
        )

    cross_engine = evaluate_cross_engine(profile, module_scores)

    # Fusion. Legacy profiles use their declared split; evidence-typed
    # profiles give each independent study group one vote. Either way only
    # modules that actually scored take part, and the remaining weight is
    # renormalised over them.
    declared = profile.combined_weights
    weights = {
        name: w for name, w in declared.items()
        # A module whose anchor organism is absent has withdrawn its claim, so
        # it must not carry weight in the fused score either. Leaving it in
        # would let a study's accessory taxa vote for a pattern whose defining
        # organism this sample does not have.
        if name in module_scores
        and module_scores[name].score is not None
        and not module_scores[name].anchor_absent
    }
    total = sum(weights.values())
    combined_score = (
        sum(weights[name] * module_scores[name].score for name in weights) / total
        if total > 0
        else None
    )
    combined = combined_null(
        {k: v for k, v in nulls.items() if k in weights}, weights
    )
    combined_percentile = percentile_against(combined, combined_score)
    interval = bootstrap_interval(
        module_scores=module_scores, weights=weights, nulls=nulls
    )

    confidence = grade_confidence(
        profile=profile,
        module_scores=module_scores,
        references=references,
        measurements=measurements,
        cross_engine=cross_engine,
    )
    vector = build_confidence_vector(
        profile=profile,
        module_scores=module_scores,
        weights=weights,
        references=references,
        measurements=measurements,
        cross_engine=cross_engine,
        interval=interval,
        empirical_specificity=empirical_specificity,
        literature_uniqueness_prior=literature_uniqueness_prior,
    )
    abstention = check_abstention(profile, measurements, references)

    return ProfileResult(
        profile=profile,
        modules=module_scores,
        combined_score=combined_score,
        combined_percentile=combined_percentile,
        confidence=confidence,
        abstention=abstention,
        cross_engine=tuple(cross_engine),
        stratum=stratum,
        dysbiosis=dysbiosis,
        reference_manifest=dict(reference_manifest or {}),
        vector=vector,
        weights={k: v / total for k, v in weights.items()} if total > 0 else {},
    )


# --------------------------------------------------------------------------- #
# ecological features
# --------------------------------------------------------------------------- #


def shannon(abundances: Sequence[float]) -> float:
    total = sum(v for v in abundances if v > 0)
    if total <= 0:
        return 0.0
    return -sum((v / total) * math.log(v / total) for v in abundances if v > 0)


def richness(abundances: Sequence[float], *, threshold: float = 0.0) -> int:
    return sum(1 for v in abundances if v > threshold)


def ecological_metrics(abundances: Mapping[str, float]) -> dict[str, float]:
    values = list(abundances.values())
    diversity = shannon(values)
    observed = richness(values)
    return {
        "shannon_diversity": diversity,
        "richness": float(observed),
        "evenness": diversity / math.log(observed) if observed > 1 else 0.0,
    }
