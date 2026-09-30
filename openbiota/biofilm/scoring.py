"""Deterministic rank arithmetic, the two named proxies, and bootstrap intervals.

The arithmetic in :func:`midrank` and :func:`rank_panel` is normative: it is
transcribed from section 8.11 of the build spec and must not be "improved".
Its behaviour on edge cases is the point. In particular a constant reference
raises rather than returning 50, because a reference in which everyone has
the same value cannot rank anybody.

Two stages, both using the same midrank operation:

1. Each feature is ranked against that feature's column in the reference.
2. Those percentiles are averaged within fixed dependence groups, then
   across groups, and the resulting aggregate is itself ranked against the
   aggregates of every reference participant computed the same way.

Stage 2 is what makes the output a percentile rather than an average of
percentiles. It also forces the reference to be transformed with the same
frozen full-reference transforms as the sample, which is why
:func:`rank_panel` computes ``aggregate_reference`` from ``reference_rows``
rather than accepting a precomputed distribution.

Nothing here knows what a biofilm is. The labels, the evidence and the
refusal to cancel one axis against another live in the engine; this module
only guarantees that identical inputs produce identical, tied, reproducible
numbers.
"""

from __future__ import annotations

import hashlib
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from math import isfinite
from statistics import mean
from typing import Final

#: Engineering default from section 8.5.
BOOTSTRAP_REPLICATES: Final = 200
#: Below this many valid replicates out of the attempted count, the interval
#: is withheld entirely rather than reported from a thin sample.
MIN_VALID_REPLICATES: Final = 180
#: Minimum independent reference participants before a percentile is shown.
#: An exploratory display minimum, not a claim of clinical adequacy.
MIN_REFERENCE_N: Final = 30
#: Stored numerical precision. Ties must survive both ranking stages, so
#: values are rounded once, here, before any comparison happens.
VALUE_PRECISION: Final = 10

#: Versioned seed string. Changing it changes every interval, so it is
#: pinned and included in the calibration hash.
SEED_STRING: Final = "openbiota.biofilm.bootstrap/1"


class ReferenceNonDiscriminating(ValueError):
    """Every reference participant has the same value, so no rank exists."""


class MaskMismatch(ValueError):
    """The sample's feature set differs from the calibration's frozen mask."""


def midrank(value: float, reference: list[float] | tuple[float, ...]) -> float:
    """Percentile of ``value`` within ``reference``, counting ties at half weight.

    ``P(x) = 100 * (count(r < x) + 0.5 * count(r == x)) / n``

    The half-weight on ties is what makes a zero-inflated reference behave
    sensibly: if half the reference is zero and the sample is zero, the
    result is 25 rather than 0 or 50. That number is a rank, not a
    detection probability, and the caller is responsible for printing
    "not detected above this assay's limit" beside it.
    """
    if not reference:
        raise ValueError("empty reference")
    if not isfinite(value):
        raise ValueError("missing or nonfinite input")
    if any(not isfinite(x) for x in reference):
        raise ValueError("nonfinite reference")
    if min(reference) == max(reference):
        raise ReferenceNonDiscriminating("reference_non_discriminating")
    below = sum(1 for x in reference if x < value)
    equal = sum(1 for x in reference if x == value)
    return 100.0 * (below + 0.5 * equal) / len(reference)


def _round(value: float) -> float:
    """Round to stored precision so equal values compare equal."""
    return round(float(value), VALUE_PRECISION)


def _midrank_sorted(value: float, ordered: list[float]) -> float:
    """:func:`midrank` against a pre-sorted column, in O(log n).

    Exactly equivalent to the linear form: ``bisect_left`` counts strictly
    smaller values and ``bisect_right - bisect_left`` counts ties. The
    reference bootstrap ranks every reference participant inside every
    replicate, which is quadratic with the linear scan and takes hours at
    n=1528; this makes it seconds without changing a single output value.
    :func:`self_test` asserts the two agree.
    """
    n = len(ordered)
    if n == 0:
        raise ValueError("empty reference")
    if not isfinite(value):
        raise ValueError("missing or nonfinite input")
    if ordered[0] == ordered[-1]:
        raise ReferenceNonDiscriminating("reference_non_discriminating")
    below = bisect_left(ordered, value)
    equal = bisect_right(ordered, value) - below
    return 100.0 * (below + 0.5 * equal) / n


def rank_panel(
    sample: dict[str, float],
    reference_rows: list[dict[str, float]],
    groups: tuple[tuple[str, ...], ...] | list[list[str]],
) -> dict[str, object]:
    """Two-stage grouped ranking of one sample against a reference cohort.

    ``groups`` is the frozen dependence structure. Features inside one group
    are averaged before groups are averaged, so a mechanism represented by
    six correlated genes does not outweigh a mechanism represented by one.

    Raises rather than guessing on every degenerate input: an empty group, a
    duplicated feature, a mask that differs from the calibration, a missing
    or nonfinite value, or a constant reference column. Each of those has a
    way to be silently wrong, and silence is what the spec forbids.
    """
    groups_t = tuple(tuple(g) for g in groups)
    names = [name for group in groups_t for name in group]
    if not groups_t or any(not group for group in groups_t):
        raise ValueError("empty panel/group")
    if len(names) != len(set(names)):
        raise ValueError("duplicate feature in dependence groups")
    if set(sample) != set(names):
        raise MaskMismatch(
            f"sample mask differs from calibration: "
            f"sample has {sorted(set(sample) - set(names)) or 'no extras'}, "
            f"missing {sorted(set(names) - set(sample)) or 'nothing'}"
        )
    if not reference_rows:
        raise ValueError("empty reference")
    if any(set(row) != set(names) for row in reference_rows):
        raise MaskMismatch("reference mask differs from calibration")
    rows = [sample, *reference_rows]
    if any(
        row.get(name) is None or not isfinite(float(row[name]))
        for row in rows
        for name in names
    ):
        raise ValueError("missing/nonfinite value")

    sample_r = {n: _round(sample[n]) for n in names}
    reference_r = [{n: _round(row[n]) for n in names} for row in reference_rows]
    columns = {n: [row[n] for row in reference_r] for n in names}
    # Sorted once, reused for every row. Both stages rank against the same
    # fixed reference columns, so the sort is loop-invariant.
    ordered = {n: sorted(col) for n, col in columns.items()}
    for n, col in ordered.items():
        if any(not isfinite(x) for x in col):
            raise ValueError("nonfinite reference")
        if col[0] == col[-1]:
            raise ReferenceNonDiscriminating(
                f"reference column {n!r} is constant; no informative percentile"
            )

    def aggregate(row: dict[str, float]) -> float:
        return mean(
            mean(_midrank_sorted(row[n], ordered[n]) for n in group)
            for group in groups_t
        )

    aggregate_reference = [aggregate(row) for row in reference_r]
    raw = aggregate(sample_r)
    if min(aggregate_reference) == max(aggregate_reference):
        raise ReferenceNonDiscriminating(
            "reference aggregate is constant; no informative axis percentile"
        )
    return {
        "raw_index": raw,
        "reference_percentile": midrank(raw, aggregate_reference),
        "feature_percentiles": {
            n: _midrank_sorted(sample_r[n], ordered[n]) for n in names
        },
        "reference_n": len(reference_rows),
        "aggregate_reference": aggregate_reference,
    }


# --------------------------------------------------------------------------
# Bootstrap intervals
# --------------------------------------------------------------------------


def _seed(
    *,
    input_hashes: tuple[str, ...],
    registry_release: str,
    calibration_hash: str,
    interval_kind: str,
) -> int:
    """Deterministic seed from the inputs that legitimately change an interval.

    Derived from SHA-256 so the same sample, registry and calibration always
    produce the same interval, and so an interval cannot be quietly
    re-rolled until it looks tighter.
    """
    material = "|".join(
        (SEED_STRING, registry_release, calibration_hash, interval_kind, *input_hashes)
    )
    digest = hashlib.sha256(material.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


def _quantile(sorted_values: list[float], q: float) -> float:
    """Linear-interpolation quantile, documented method per section 8.5."""
    if not sorted_values:
        raise ValueError("empty sample")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = q * (len(sorted_values) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


@dataclass(frozen=True, slots=True)
class Interval:
    """One bootstrap interval, with its validity accounting.

    ``kind`` is either ``analytical`` (resampling read fragments, holding
    the calibration fixed) or ``reference`` (resampling reference
    participants, holding the sample measurement fixed). They answer
    different questions and section 8.5 forbids pooling them into one
    number, so they are separate objects with separate labels.
    """

    kind: str
    low: float | None
    high: float | None
    attempted: int
    valid: int
    status: str
    failure_reasons: dict[str, int] = field(default_factory=dict)

    @property
    def label(self) -> str:
        return {
            "analytical": "sampling interval",
            "reference": "reference interval",
        }.get(self.kind, self.kind)

    @property
    def invalid_fraction(self) -> float:
        return 0.0 if not self.attempted else 1.0 - (self.valid / self.attempted)

    def to_json(self) -> dict[str, object] | None:
        if self.low is None or self.high is None:
            return {
                "kind": self.kind,
                "label": self.label,
                "low": None,
                "high": None,
                "status": self.status,
                "attempted": self.attempted,
                "valid": self.valid,
                "failure_reasons": dict(self.failure_reasons),
            }
        return {
            "kind": self.kind,
            "label": self.label,
            "low": round(self.low, 2),
            "high": round(self.high, 2),
            "status": self.status,
            "attempted": self.attempted,
            "valid": self.valid,
            "invalid_fraction": round(self.invalid_fraction, 4),
            "failure_reasons": dict(self.failure_reasons),
            "meaning": (
                "Where this rank would land if the "
                + ("reference cohort" if self.kind == "reference" else "sequencing")
                + " were repeated. Not a statement about clinical risk."
            ),
        }


def reference_interval(
    sample: dict[str, float],
    reference_rows: list[dict[str, float]],
    groups: tuple[tuple[str, ...], ...],
    *,
    input_hashes: tuple[str, ...],
    registry_release: str,
    calibration_hash: str,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> Interval:
    """Resample reference participants, holding the sample measurement fixed.

    Participants are resampled whole, retaining their paired module vectors,
    because a participant is the independent unit. Resampling feature values
    independently would break their correlation and produce a falsely tight
    interval.

    Both ranking stages are rebuilt inside each replicate under the identical
    algorithm. A replicate whose resampled reference happens to be constant
    is invalid, and is counted rather than repaired: substituting zero or
    dropping a module to rescue it would bias the interval.
    """
    import random

    rng = random.Random(
        _seed(
            input_hashes=input_hashes,
            registry_release=registry_release,
            calibration_hash=calibration_hash,
            interval_kind="reference",
        )
    )
    n = len(reference_rows)
    values: list[float] = []
    failures: dict[str, int] = {}
    for _ in range(replicates):
        resampled = [reference_rows[rng.randrange(n)] for _ in range(n)]
        try:
            result = rank_panel(sample, resampled, groups)
        except ReferenceNonDiscriminating:
            failures["reference_non_discriminating"] = (
                failures.get("reference_non_discriminating", 0) + 1
            )
            continue
        except (ValueError, MaskMismatch) as exc:
            key = type(exc).__name__
            failures[key] = failures.get(key, 0) + 1
            continue
        values.append(float(result["reference_percentile"]))

    return _finish(values, replicates, failures, kind="reference")


def analytical_interval(
    *,
    available: bool,
    reason: str = "summary_only_cache",
) -> Interval:
    """Fragment-level resampling, or an explicit refusal to invent one.

    Analytical resampling requires the read fragments themselves so pairs
    can be kept together. When only a summary table is available the
    interval is unavailable, not estimated: section 8.5 is explicit that a
    summary-only cache cannot yield a fragment-resampling interval, and
    BF-T047 tests exactly this.
    """
    if not available:
        return Interval(
            kind="analytical",
            low=None,
            high=None,
            attempted=0,
            valid=0,
            status=reason,
            failure_reasons={reason: 1},
        )
    return Interval(
        kind="analytical",
        low=None,
        high=None,
        attempted=0,
        valid=0,
        status="not_implemented_for_taxonomic_proxy",
    )


def _finish(
    values: list[float],
    attempted: int,
    failures: dict[str, int],
    *,
    kind: str,
) -> Interval:
    """Apply the validity threshold and build the interval."""
    valid = len(values)
    if valid < MIN_VALID_REPLICATES:
        return Interval(
            kind=kind,
            low=None,
            high=None,
            attempted=attempted,
            valid=valid,
            status="bootstrap_unstable",
            failure_reasons=failures,
        )
    ordered = sorted(values)
    return Interval(
        kind=kind,
        low=_quantile(ordered, 0.025),
        high=_quantile(ordered, 0.975),
        attempted=attempted,
        valid=valid,
        status="conditional_on_valid_replicates",
        failure_reasons=failures,
    )


# --------------------------------------------------------------------------
# The two named proxies
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Proxy:
    """A fixed, versioned, transparent index over taxonomic features.

    These are engineering hypotheses, not published predictive models. The
    ``id`` carries a version because changing the feature set requires a new
    index and a new calibration rather than an edit in place.
    """

    id: str
    card: str
    heading: str
    label: str
    #: Fixed feature set. All of them are required for the named proxy.
    features: tuple[str, ...]
    #: Fixed dependence structure passed to rank_panel.
    groups: tuple[tuple[str, ...], ...]
    source_ids: tuple[str, ...]
    module_id: str
    direction_note: str
    transport_change: str
    what_it_is_not: str

    @property
    def n_features(self) -> int:
        return len(self.features)


#: H-C. Two features, one group: they are one association pattern, not two
#: independent findings.
COMMUNITY_PROXY: Final = Proxy(
    id="BF-PROXY-COMMUNITY-1.0",
    card="H-C",
    heading="Concerning biofilm potential",
    label="Biofilm-associated community pattern",
    features=("ecoli_complex", "m_gnavus"),
    groups=(("ecoli_complex", "m_gnavus"),),
    source_ids=("BF-S01",),
    module_id="BF-M18",
    direction_note=(
        "Higher means more of this two-feature association pattern, which "
        "BF-S01 reported enriched in biofilm-positive participants."
    ),
    transport_change="genus_complex_to_species_complex",
    what_it_is_not=(
        "Not measured harmful biofilm, not a stool biofilm probability, and "
        "not a pathotype. BF-S01's 16S Escherichia/Shigella association cannot "
        "be resolved identically by a WGS species-complex measurement, and not "
        "every Escherichia/Shigella read belongs to pathogenic E. coli."
    ),
)

#: P-E. Four genera, one group: section 8.7 requires them treated as one
#: correlated ecological group, not four independent protective mechanisms.
ECOLOGY_PROXY: Final = Proxy(
    id="BF-PROXY-ECOLOGY-1.0",
    card="P-E",
    heading="Protective biofilm/ecosystem support",
    label="Protective ecological support",
    features=("faecalibacterium", "coprococcus", "subdoligranulum", "blautia"),
    groups=(("faecalibacterium", "coprococcus", "subdoligranulum", "blautia"),),
    source_ids=("BF-S01",),
    module_id="BF-M19",
    direction_note=(
        "Higher means more of this narrow community pattern, which BF-S01 "
        "reported relatively depleted in biofilm-positive biopsies."
    ),
    transport_change="biopsy_association_to_stool_measurement",
    what_it_is_not=(
        "Not a measurement of beneficial biofilm amount. Low values do not "
        "diagnose a missing protective biofilm, a species deficiency or a need "
        "for transplant, and this proxy does not instruct supplementation with "
        "these organisms."
    ),
)

PROXIES: Final = (COMMUNITY_PROXY, ECOLOGY_PROXY)

#: Measurements explicitly forbidden from entering the fixed proxies.
#: They remain useful elsewhere in the report; they are just not these indices.
FORBIDDEN_PROXY_FEATURES: Final = frozenset(
    {
        "shannon_diversity",
        "firmicutes_bacteroidetes_ratio",
        "butyrate_genes",
        "all_lactobacillus",
    }
)


def calibration_hash(
    proxy: Proxy, reference_rows: list[dict[str, float]], reference_id: str
) -> str:
    """Stable hash of a frozen calibration.

    Covers the proxy version, its exact feature order, its dependence
    structure, the reference identity and the reference size. Two axes, or
    one axis under two missingness masks, therefore get different hashes and
    cannot be compared as though they were the same scale.
    """
    material = "|".join(
        (
            proxy.id,
            ",".join(proxy.features),
            ";".join(",".join(g) for g in proxy.groups),
            reference_id,
            str(len(reference_rows)),
        )
    )
    return hashlib.sha256(material.encode("utf-8")).hexdigest()[:16]


def self_test() -> None:
    """Fixtures from section 8.6 and BF-T031/T032/T125.

    Run at import time by the acceptance suite and by ``biofilm validate``.
    """
    # Section 8.6 midrank vectors.
    ref = [0.0, 0.0, 2.0, 4.0]
    assert midrank(0, ref) == 25.0, midrank(0, ref)
    assert midrank(2, ref) == 62.5, midrank(2, ref)
    assert midrank(5, ref) == 100.0, midrank(5, ref)

    # A constant reference has no ranks to give.
    try:
        midrank(1.0, [3.0, 3.0, 3.0])
    except ReferenceNonDiscriminating:
        pass
    else:  # pragma: no cover
        raise AssertionError("constant reference must raise")

    # Ties survive both stages: two identical rows get identical output.
    rows = [
        {"a": 1.0, "b": 1.0},
        {"a": 2.0, "b": 2.0},
        {"a": 3.0, "b": 3.0},
        {"a": 4.0, "b": 9.0},
    ]
    groups = (("a", "b"),)
    first = rank_panel({"a": 2.0, "b": 2.0}, rows, groups)
    second = rank_panel({"a": 2.0, "b": 2.0}, rows, groups)
    assert first == second
    assert first["reference_percentile"] == rank_panel(rows[1], rows, groups)[
        "reference_percentile"
    ]

    # A single group is supported and is not an error.
    single = rank_panel({"a": 2.5}, [{"a": 1.0}, {"a": 2.0}, {"a": 4.0}], (("a",),))
    assert 0.0 <= float(single["reference_percentile"]) <= 100.0

    # Nonfinite input is refused rather than coerced.
    for bad in (float("nan"), float("inf")):
        try:
            rank_panel({"a": bad}, [{"a": 1.0}, {"a": 2.0}], (("a",),))
        except ValueError:
            pass
        else:  # pragma: no cover
            raise AssertionError(f"{bad} must be refused")

    # A mask mismatch is an error, never a silent drop-and-rescale.
    try:
        rank_panel({"a": 1.0}, [{"a": 1.0, "b": 2.0}], (("a", "b"),))
    except MaskMismatch:
        pass
    else:  # pragma: no cover
        raise AssertionError("mask mismatch must raise")

    # Duplicate feature across groups would double-count.
    try:
        rank_panel({"a": 1.0}, [{"a": 1.0}], (("a",), ("a",)))
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("duplicate feature must raise")

    # The binary-search fast path must agree with the linear definition on
    # every value, including ties, below-minimum and above-maximum.
    import random as _random

    _rng = _random.Random(20260917)
    for _ in range(200):
        col = [float(_rng.randrange(0, 6)) for _ in range(_rng.randrange(2, 40))]
        if min(col) == max(col):
            continue
        ordered_col = sorted(col)
        for probe in {*col, min(col) - 1.0, max(col) + 1.0, 2.5}:
            assert midrank(probe, col) == _midrank_sorted(probe, ordered_col), (
                probe,
                col,
            )

    # Oppositely ordered features with a constant aggregate cannot produce a
    # meaningful composite rank (BF-T105).
    opposed = [
        {"a": 1.0, "b": 4.0},
        {"a": 2.0, "b": 3.0},
        {"a": 3.0, "b": 2.0},
        {"a": 4.0, "b": 1.0},
    ]
    try:
        rank_panel({"a": 2.0, "b": 3.0}, opposed, (("a", "b"),))
    except ReferenceNonDiscriminating:
        pass
    else:  # pragma: no cover
        raise AssertionError("constant aggregate must raise")
