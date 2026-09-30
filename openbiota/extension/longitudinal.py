"""Longitudinal views — A13, BUILD_SPEC_v0.8.3 §6.5, and the recovery
envelope of §6.4.

The hard part of comparing two samples is not the distance. It is deciding
whether the two are comparable at all, and saying so when they are not.

A **comparability fingerprint** records everything that could make two
results differ for a reason that is not biology: the specimen and
participant, the collection date, the extraction and library method, host
filtering, read depth, and every tool, database and model version in the
lane. Two samples whose fingerprints differ in a version are a *labelled
segment boundary*, not a change in the person, and the comparison says
`descriptive_only` rather than quietly drawing a trend line across it.

Reprocessing the same FASTQs is not a new timepoint. The specimen identity
travels with the specimen, so a second analysis of one stool sample cannot
become two points on a chart.

Recovery is only computed where it means something: a recorded
perturbation, at least three eligible baseline samples, and an observed
displacement larger than the baseline's own spread. Anything less shows
displacement without claiming return.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import date
from typing import Any, Final

from openbiota.extension.schema import fingerprint

METHOD_DISTANCE: Final = "ext083.longitudinal.distance/1.0"
METHOD_RECOVERY: Final = "ext083.longitudinal.recovery/1.0"

#: How far two results may be compared.
#:
#: ``native``                 same lane, same versions: compare directly
#: ``uniformly_reprocessed``  re-run together in one independent lane
#: ``bridge_validated``       different lanes with a measured conversion
#: ``descriptive_only``       shown side by side; no distance is quoted
COMPARABILITY: Final[tuple[str, ...]] = (
    "native", "uniformly_reprocessed", "bridge_validated", "descriptive_only",
)

#: Fingerprint fields that, if they differ, break native comparability.
#: Read depth is deliberately not among them: it varies between runs and is
#: handled by the detection rule rather than by refusing the comparison.
BLOCKING_FIELDS: Final[tuple[str, ...]] = (
    "profiler", "profiler_version", "database", "database_version",
    "extraction_method", "library_method", "host_filter", "specimen_type",
)

#: Fields that are recorded and shown but do not block a comparison.
CONTEXT_FIELDS: Final[tuple[str, ...]] = (
    "collection_date", "stool_form", "storage", "read_depth",
    "diet_events", "antibiotic_events", "supplement_events",
)

#: A taxon counts as present at or above this share of the sample. Fixed and
#: documented, per §6.5: a presence set that moves with the sequencing depth
#: turns a deeper run into a "richer" gut.
PRESENCE_FLOOR: Final = 0.001  # 0.1% of the abundance vector

#: The four follow-up templates of §6.5. No universal retesting interval:
#: each says what the reader wants to learn and what would have changed.
FOLLOW_UP_TEMPLATES: Final[tuple[dict[str, Any], ...]] = (
    {
        "template_id": "routine_review",
        "trigger": "no recorded intervention or exposure",
        "question": "Has anything drifted on its own?",
        "interval": "6 to 12 months",
        "rationale": (
            "An adult gut community is fairly stable month to month without an exposure to "
            "move it, so a sooner repeat mostly measures sampling variation."
        ),
    },
    {
        "template_id": "new_probiotic_or_prebiotic",
        "trigger": "a probiotic or prebiotic started since the last sample",
        "question": (
            "Did the organism establish, and did the functions it was taken for move with it?"
        ),
        "interval": "4 to 8 weeks on the product, and again 4 weeks after stopping",
        "rationale": (
            "Most probiotic strains are detectable while being taken and fade within weeks of "
            "stopping. Sampling only while taking it cannot separate passage from "
            "colonisation, which is why the second point matters."
        ),
    },
    {
        "template_id": "major_dietary_change",
        "trigger": "a substantial, sustained change in what is eaten",
        "question": "Did the fibre-using and fermentation capacity follow the diet?",
        "interval": "4 to 6 weeks after the change is established",
        "rationale": (
            "Composition responds to a sustained dietary change within days, but the first "
            "week reflects the transition rather than the new steady state."
        ),
    },
    {
        "template_id": "recent_antibiotics",
        "trigger": "an antibiotic course since the last sample",
        "question": "Has the community returned toward where it was before the course?",
        "interval": "at least 8 weeks after the last dose, then again at 6 months",
        "rationale": (
            "Diversity and the butyrate producers fall during a course and recover over "
            "weeks to months, incompletely in some people. Sampling during or just after a "
            "course measures the disturbance, not the outcome."
        ),
    },
)
TEMPLATE_BY_ID: Final[Mapping[str, Mapping[str, Any]]] = {
    t["template_id"]: t for t in FOLLOW_UP_TEMPLATES
}


class LongitudinalError(ValueError):
    """A comparison that would compare things that are not comparable."""


# --------------------------------------------------------------------------- #
# distances
# --------------------------------------------------------------------------- #


def bray_curtis(p: Mapping[str, float], q: Mapping[str, float]) -> float | None:
    """`sum|p_i - q_i| / sum(p_i + q_i)` over the union of both vectors.

    None when both vectors are empty: there is no distance between two
    absences, and zero would read as perfect agreement.
    """
    keys = set(p) | set(q)
    if not keys:
        return None
    numerator = sum(abs(float(p.get(k, 0.0)) - float(q.get(k, 0.0))) for k in keys)
    denominator = sum(float(p.get(k, 0.0)) + float(q.get(k, 0.0)) for k in keys)
    if denominator <= 0:
        return None
    return numerator / denominator


def presence_set(vector: Mapping[str, float], *, floor: float = PRESENCE_FLOOR) -> frozenset[str]:
    """Everything at or above the documented detection floor."""
    total = sum(float(v) for v in vector.values() if float(v) > 0)
    if total <= 0:
        return frozenset()
    return frozenset(k for k, v in vector.items() if float(v) / total >= floor)


def jaccard_distance(
    p: Mapping[str, float], q: Mapping[str, float], *, floor: float = PRESENCE_FLOOR
) -> float | None:
    """`1 - |P and Q| / |P or Q|` over the presence sets.

    An empty union is undefined, not perfect agreement: two samples in which
    nothing was detected have not been shown to be identical.
    """
    a, b = presence_set(p, floor=floor), presence_set(q, floor=floor)
    union = a | b
    if not union:
        return None
    return 1.0 - len(a & b) / len(union)


def closure(vector: Mapping[str, float]) -> dict[str, float]:
    """A vector rescaled to sum to one, dropping non-positive entries."""
    positive = {k: float(v) for k, v in vector.items() if float(v) > 0}
    total = sum(positive.values())
    if total <= 0:
        return {}
    return {k: v / total for k, v in positive.items()}


def mean_composition(vectors: Sequence[Mapping[str, float]]) -> dict[str, float]:
    """The arithmetic mean composition, followed by closure, per §6.4."""
    if not vectors:
        return {}
    keys: set[str] = set()
    for vector in vectors:
        keys |= set(vector)
    mean = {k: sum(float(v.get(k, 0.0)) for v in vectors) / len(vectors) for k in keys}
    return closure(mean)


# --------------------------------------------------------------------------- #
# comparability
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Fingerprint:
    """Everything that could make two results differ for a non-biological reason."""

    sample_id: str
    participant_id: str | None = None
    specimen_id: str | None = None
    fields: Mapping[str, Any] = field(default_factory=dict)

    @property
    def identity(self) -> str:
        """The specimen this result describes.

        Two analyses of one stool sample share this, which is what stops a
        reprocessing run from becoming a second biological timepoint.
        """
        return str(self.specimen_id or self.sample_id)

    @property
    def digest(self) -> str:
        return fingerprint(
            "a13.fingerprint/1", self.identity, self.participant_id,
            sorted((k, str(v)) for k, v in self.fields.items()),
        )

    @property
    def blocking_digest(self) -> str:
        """Only the fields whose change breaks a direct comparison."""
        return fingerprint(
            "a13.blocking/1",
            sorted((k, str(self.fields.get(k))) for k in BLOCKING_FIELDS),
        )

    def differences(self, other: Fingerprint) -> dict[str, tuple[Any, Any]]:
        out: dict[str, tuple[Any, Any]] = {}
        for key in (*BLOCKING_FIELDS, *CONTEXT_FIELDS):
            mine, theirs = self.fields.get(key), other.fields.get(key)
            if mine != theirs:
                out[key] = (mine, theirs)
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "participant_id": self.participant_id,
            "specimen_id": self.specimen_id,
            "specimen_identity": self.identity,
            "fields": dict(self.fields),
            "digest": self.digest,
            "blocking_digest": self.blocking_digest,
        }


def comparability(
    a: Fingerprint, b: Fingerprint, *, reprocessed_together: bool = False,
    bridge_validated: bool = False,
) -> tuple[str, str]:
    """How far these two may be compared, and why.

    `uniformly_reprocessed` is a claim about an independent analysis lane
    and only the caller knows whether one was run; it is not inferred.
    """
    differences = {k: v for k, v in a.differences(b).items() if k in BLOCKING_FIELDS}
    if not differences:
        return "native", "same lane, same tool and database versions"
    named = ", ".join(
        f"{k.replace('_', ' ')} {old!r} then {new!r}" for k, (old, new) in sorted(differences.items())
    )
    if reprocessed_together:
        return "uniformly_reprocessed", (
            f"these differ in {named}, and both samples were re-run together in one lane, "
            "so the comparison is made on the reprocessed results"
        )
    if bridge_validated:
        return "bridge_validated", (
            f"these differ in {named}, with a measured conversion between the two lanes"
        )
    return "descriptive_only", (
        f"these differ in {named}. A number computed across that change would measure the "
        "change of method as much as any change in you, so the results are shown side by "
        "side and no distance is quoted"
    )


def is_repeat_processing(a: Fingerprint, b: Fingerprint) -> bool:
    """Whether these two results describe the same physical specimen."""
    return a.identity == b.identity


# --------------------------------------------------------------------------- #
# a series
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TimePoint:
    """One result in a series."""

    fingerprint: Fingerprint
    collected_at: date | None
    abundances: Mapping[str, float]
    metrics: Mapping[str, float] = field(default_factory=dict)
    events: tuple[str, ...] = ()

    @property
    def sample_id(self) -> str:
        return self.fingerprint.sample_id


def as_date(value: Any) -> date | None:
    """A date from whatever the manifest recorded, or None if it is not one."""
    if isinstance(value, date):
        return value
    try:
        return date.fromisoformat(str(value)[:10])
    except (TypeError, ValueError):
        return None


@dataclass
class Series:
    """One participant's timepoints, deduplicated by specimen and ordered."""

    participant_id: str
    points: list[TimePoint] = field(default_factory=list)
    #: Results dropped because they re-analyse a specimen already present.
    repeat_processing: list[str] = field(default_factory=list)

    @classmethod
    def build(cls, participant_id: str, points: Iterable[TimePoint]) -> Series:
        series = cls(participant_id=participant_id)
        seen: dict[str, TimePoint] = {}
        for point in points:
            identity = point.fingerprint.identity
            if identity in seen:
                series.repeat_processing.append(point.sample_id)
                continue
            seen[identity] = point
        series.points = sorted(
            seen.values(), key=lambda p: (p.collected_at or date.min, p.sample_id),
        )
        return series

    @property
    def n_timepoints(self) -> int:
        return len(self.points)

    @property
    def latest(self) -> TimePoint | None:
        return self.points[-1] if self.points else None

    def last(self, n: int = 3) -> list[TimePoint]:
        """The last `n` measurements, for the compact cards of §6.5."""
        return self.points[-n:]

    def segments(self) -> list[dict[str, Any]]:
        """Runs of consecutive timepoints that share a blocking fingerprint.

        A tool or database upgrade starts a new segment. Charts draw within
        a segment and break across one, rather than running a line through
        a change of method.
        """
        out: list[dict[str, Any]] = []
        for point in self.points:
            digest = point.fingerprint.blocking_digest
            if out and out[-1]["blocking_digest"] == digest:
                out[-1]["sample_ids"].append(point.sample_id)
            else:
                out.append({
                    "blocking_digest": digest,
                    "sample_ids": [point.sample_id],
                    "lane": {k: point.fingerprint.fields.get(k) for k in BLOCKING_FIELDS},
                })
        return out

    def pairwise(self) -> list[dict[str, Any]]:
        """Consecutive comparisons, each with its own comparability verdict."""
        out: list[dict[str, Any]] = []
        for earlier, later in zip(self.points, self.points[1:], strict=False):
            state, why = comparability(earlier.fingerprint, later.fingerprint)
            days = (
                (later.collected_at - earlier.collected_at).days
                if earlier.collected_at and later.collected_at else None
            )
            record: dict[str, Any] = {
                "from": earlier.sample_id, "to": later.sample_id,
                "days_apart": days, "comparability": state, "comparability_reason": why,
                "events_between": list(later.events),
                "bray_curtis": None, "jaccard_distance": None,
            }
            if state != "descriptive_only":
                record["bray_curtis"] = bray_curtis(earlier.abundances, later.abundances)
                record["jaccard_distance"] = jaccard_distance(
                    earlier.abundances, later.abundances,
                )
                record["metric_changes"] = {
                    key: later.metrics[key] - earlier.metrics[key]
                    for key in set(earlier.metrics) & set(later.metrics)
                }
            out.append(record)
        return out


# --------------------------------------------------------------------------- #
# recovery (§6.4)
# --------------------------------------------------------------------------- #

#: A recovery estimate needs at least this many eligible baseline samples.
MIN_BASELINE_SAMPLES: Final = 3


@dataclass(frozen=True)
class Recovery:
    """Return toward a pre-event community, or exactly why there is no number."""

    state: str
    baseline_samples: tuple[str, ...] = ()
    envelope: float | None = None
    peak_distance: float | None = None
    peak_sample: str | None = None
    trajectory: tuple[tuple[str, float, float | None], ...] = ()
    reason: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "method_id": METHOD_RECOVERY,
            "state": self.state,
            "baseline_samples": list(self.baseline_samples),
            "baseline_envelope": self.envelope,
            "peak_distance": self.peak_distance,
            "peak_sample": self.peak_sample,
            "trajectory": [
                {"sample_id": s, "distance_from_baseline": d, "return_percent": r}
                for s, d, r in self.trajectory
            ],
            "reason": self.reason,
            "label": "Return toward the pre-event community",
            "limitations": [
                "An empirical envelope from the baseline samples themselves, not a 95% "
                "confidence interval.",
                "A retrospective trajectory: a later sample can revise the observed peak, "
                "and reports already issued are not rewritten.",
                "Returning toward the earlier community does not establish that the earlier "
                "community was healthy.",
            ],
        }


def recovery(
    series: Series, *, event_sample_id: str, technical_noise_bound: float | None = None,
) -> Recovery:
    """Return toward the pre-event community, where the data support it.

    Baseline is the mean composition of the eligible pre-event samples,
    closed. The envelope is the largest distance any baseline sample sits
    from that mean, widened to a separately validated technical bound if one
    exists. A displacement inside the envelope is not a displacement.
    """
    ids = [p.sample_id for p in series.points]
    if event_sample_id not in ids:
        return Recovery("unavailable", reason=f"{event_sample_id} is not in this series")
    index = ids.index(event_sample_id)
    before = series.points[:index]
    after = series.points[index:]
    eligible = [p for p in before if p.abundances]
    if len(eligible) < MIN_BASELINE_SAMPLES:
        return Recovery(
            "insufficient_baseline",
            baseline_samples=tuple(p.sample_id for p in eligible),
            reason=(
                f"{len(eligible)} sample(s) before the event; an envelope needs at least "
                f"{MIN_BASELINE_SAMPLES}, because a spread cannot be estimated from fewer. "
                "Changes and similarities are still shown; only the return percentage is not."
            ),
        )
    baseline = mean_composition([closure(p.abundances) for p in eligible])
    spreads = [bray_curtis(closure(p.abundances), baseline) for p in eligible]
    finite = [s for s in spreads if s is not None]
    if not finite:
        return Recovery("unavailable", reason="the baseline samples carry no comparable vector")
    envelope = max(finite)
    if technical_noise_bound is not None:
        envelope = max(envelope, float(technical_noise_bound))

    distances: list[tuple[str, float]] = []
    for point in after:
        distance = bray_curtis(closure(point.abundances), baseline)
        if distance is not None:
            distances.append((point.sample_id, distance))
    if not distances:
        return Recovery(
            "unavailable", baseline_samples=tuple(p.sample_id for p in eligible),
            envelope=envelope, reason="no comparable sample after the event",
        )
    peak_distance = max(d for _, d in distances)
    peak_sample = next(s for s, d in distances if d == peak_distance)
    if peak_distance <= envelope or peak_distance <= 0:
        return Recovery(
            "no_displacement", baseline_samples=tuple(p.sample_id for p in eligible),
            envelope=envelope, peak_distance=peak_distance, peak_sample=peak_sample,
            trajectory=tuple((s, d, None) for s, d in distances),
            reason=(
                f"the largest post-event distance ({peak_distance:.3f}) is within the "
                f"baseline's own spread ({envelope:.3f}), so nothing moved further than "
                "these samples normally differ from each other"
            ),
        )
    peak_index = next(i for i, (_, d) in enumerate(distances) if d == peak_distance)
    trajectory: list[tuple[str, float, float | None]] = []
    for i, (sample_id, distance) in enumerate(distances):
        if i < peak_index:
            trajectory.append((sample_id, distance, None))  # displacement only
        else:
            fraction = (peak_distance - distance) / peak_distance
            trajectory.append((sample_id, distance, 100.0 * min(1.0, max(0.0, fraction))))
    return Recovery(
        "measured", baseline_samples=tuple(p.sample_id for p in eligible),
        envelope=envelope, peak_distance=peak_distance, peak_sample=peak_sample,
        trajectory=tuple(trajectory),
        reason=(
            f"displacement peaked at {peak_distance:.3f} in {peak_sample}, beyond the "
            f"baseline envelope of {envelope:.3f}; return is measured from that peak onward"
        ),
    )


# --------------------------------------------------------------------------- #
# follow-up
# --------------------------------------------------------------------------- #


def follow_up(events: Sequence[str]) -> list[dict[str, Any]]:
    """Which follow-up templates this person's recorded exposures call for.

    No universal interval: the templates differ because the questions do,
    and a routine review is not the same question as a post-antibiotic one.
    """
    text = " ".join(events).lower()
    chosen: list[str] = []
    if any(w in text for w in ("antibiotic", "amoxicillin", "azithromycin", "ciprofloxacin")):
        chosen.append("recent_antibiotics")
    if any(w in text for w in ("probiotic", "prebiotic", "inulin", "synbiotic")):
        chosen.append("new_probiotic_or_prebiotic")
    if any(w in text for w in ("diet", "vegan", "vegetarian", "fibre", "fiber", "carnivore")):
        chosen.append("major_dietary_change")
    if not chosen:
        chosen.append("routine_review")
    return [dict(TEMPLATE_BY_ID[t]) for t in chosen]


def summarise(series: Series, *, event_sample_id: str | None = None) -> dict[str, Any]:
    """The whole longitudinal view for one participant.

    `event_sample_id` names a recorded perturbation. Without one there is
    no recovery estimate, because "recovery" needs something to recover
    from; the trend and the distances are reported either way.
    """
    points = series.points
    events = [e for p in points for e in p.events]
    if event_sample_id is None:
        event_sample_id = next(
            (p.sample_id for p in points if p.events), None,
        )
    return {
        "feature_id": "A13",
        "method_id": METHOD_DISTANCE,
        "participant_id": series.participant_id,
        "n_timepoints": series.n_timepoints,
        "sample_ids": [p.sample_id for p in points],
        "collection_dates": [
            p.collected_at.isoformat() if p.collected_at else None for p in points
        ],
        "repeat_processing_excluded": list(series.repeat_processing),
        "repeat_processing_note": (
            "Re-analysing the same stool sample does not create a second timepoint; those "
            "results keep the specimen's identity and are not charted as a change."
        ),
        "segments": series.segments(),
        "comparisons": series.pairwise(),
        "last_three": [p.sample_id for p in series.last(3)],
        "presence_floor": PRESENCE_FLOOR,
        "presence_floor_note": (
            f"A taxon counts as present at or above {PRESENCE_FLOOR:.1%} of the sample. The "
            "rule is fixed so that a deeper run does not appear to be a richer gut."
        ),
        "recovery": (
            recovery(series, event_sample_id=event_sample_id).to_json()
            if event_sample_id else None
        ),
        "recorded_exposures": events,
        "follow_up": follow_up(events),
        "limitations": [
            "A change between two samples is a change in what was measured, which includes "
            "how it was measured; the comparability verdict on each pair says which.",
            "Raw changes and reference-percentile changes are different quantities and are "
            "reported separately.",
            "Earlier reports are never rewritten when a later sample arrives.",
        ],
    }


def not_enough_history(sample_id: str, *, n_available: int = 1) -> dict[str, Any]:
    """The honest view when there is only one sample.

    §6.4 is explicit: do not suppress the longitudinal section merely
    because recovery cannot be estimated. This says what a second sample
    would add and what would make it comparable.
    """
    return {
        "feature_id": "A13",
        "method_id": METHOD_DISTANCE,
        "state": "single_timepoint",
        "n_timepoints": n_available,
        "sample_ids": [sample_id],
        "headline": (
            "This is your first sample, so there is nothing to compare it with yet."
            if n_available <= 1 else
            f"{n_available} samples are on file, but none of them can be compared with this one."
        ),
        "what_a_second_sample_adds": [
            "Which readings moved and which held, separately from where they sit against "
            "the reference group.",
            "Whether an organism that is high or low now is usually that way for you.",
            "Whether a change you made reached the community at all.",
        ],
        "what_makes_them_comparable": [
            f"the same collection and storage method ({', '.join(BLOCKING_FIELDS[:3])} and "
            "the rest of the lane are recorded automatically)",
            "processing in the same lane, or re-running both together if the tools have "
            "changed in between",
            "a recorded collection date, so the interval is known",
        ],
        "follow_up": follow_up([]),
        "limitations": [
            "A single sample describes one day. It cannot separate what is usual for you "
            "from what was true that week.",
        ],
    }


def percentile_change(before: float | None, after: float | None) -> dict[str, Any]:
    """A percentile change, kept apart from a raw change.

    A move from the 40th to the 60th percentile is not a doubling and is
    not a 20% increase; it is a move through the reference distribution.
    """
    if before is None or after is None:
        return {"state": "unavailable", "delta_percentile": None}
    delta = float(after) - float(before)
    return {
        "state": "measured",
        "from_percentile": float(before),
        "to_percentile": float(after),
        "delta_percentile": delta,
        "unit": "percentile points",
        "note": (
            "A move through the reference distribution, not a change in the underlying "
            "amount. The raw change is reported separately."
        ),
    }


def stability(values: Sequence[float]) -> dict[str, Any]:
    """How much a reading moves for this person, across their own samples."""
    finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if len(finite) < 2:
        return {"state": "insufficient_history", "n": len(finite)}
    return {
        "state": "measured",
        "n": len(finite),
        "mean": statistics.fmean(finite),
        "spread": max(finite) - min(finite),
        "stdev": statistics.stdev(finite) if len(finite) > 2 else None,
        "note": (
            "The spread of your own results, which is the only baseline that says whether a "
            "new value is unusual for you rather than unusual for the reference group."
        ),
    }


__all__ = [
    "BLOCKING_FIELDS",
    "COMPARABILITY",
    "CONTEXT_FIELDS",
    "FOLLOW_UP_TEMPLATES",
    "METHOD_DISTANCE",
    "METHOD_RECOVERY",
    "MIN_BASELINE_SAMPLES",
    "PRESENCE_FLOOR",
    "TEMPLATE_BY_ID",
    "Fingerprint",
    "LongitudinalError",
    "Recovery",
    "Series",
    "TimePoint",
    "as_date",
    "bray_curtis",
    "closure",
    "comparability",
    "follow_up",
    "is_repeat_processing",
    "jaccard_distance",
    "mean_composition",
    "not_enough_history",
    "percentile_change",
    "presence_set",
    "recovery",
    "stability",
    "summarise",
]
