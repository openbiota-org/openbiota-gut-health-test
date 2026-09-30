"""One score for the gut-skin axis, from the skin pattern panels.

Each published skin panel asks the same question in a different study's
words: does this stool community lean the way that study's patients leaned?
A panel's index runs 0-100 with higher meaning *more* like the disease
cohort, so the axis score is the inverse - higher is the better direction,
the way every other favourable reading in this report runs.

WHY A WEIGHTED MEAN, AND WEIGHTED BY WHAT

The panels do not carry equal evidence. One resolved three of its three
marker organisms in this sample; another resolved one of seven, because
the rest were below the reporting floor. Averaging those equally would let
a panel that saw a seventh of its own evidence pull the number as hard as
one that saw all of it, so each panel is weighted by its coverage - the
fraction of its markers this sample could actually measure.

That makes the score a statement about the evidence there is, and it makes
the weight visible: the number is reported with the coverage behind it, so
a score resting on thin panels can be read as such.

WHAT IT IS NOT

Not a diagnosis, not a probability, and not validated against skin
outcomes. No stool signature is established as diagnostic for any of these
conditions; these are directional pattern indices from published case-
control comparisons, and this score is their weighted middle. It carries
"(beta)" wherever it is printed for that reason.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

METHOD: Final = "openbiota.gut_skin_axis/1.0"

#: Where the panels live in the similarity stage's output.
_SOURCES: Final = ("skin_pattern_indices", "urticaria_pattern_indices")

#: Keys in those mappings that are not panels.
_NOT_PANELS: Final = frozenset({
    "unavailable_tasks", "non_scored_observations", "families",
})


@dataclass(frozen=True)
class AxisScore:
    """The gut-skin axis reading, with everything behind it."""

    value: float | None
    n_panels: int
    n_scored: int
    total_coverage: float
    contributors: tuple[dict[str, Any], ...]
    reason: str = ""

    @property
    def scored(self) -> bool:
        return self.value is not None

    def to_json(self) -> dict[str, Any]:
        return {
            "method": METHOD,
            "value": self.value,
            "n_panels": self.n_panels,
            "n_scored": self.n_scored,
            "total_coverage": round(self.total_coverage, 3),
            "higher_means": "favourable",
            "contributors": list(self.contributors),
            "reason": self.reason,
            "what_this_is_not": (
                "A weighted middle of published directional pattern indices, not a "
                "diagnosis, not a probability, and not validated against skin outcomes."
            ),
        }


def _panels(similarity: Mapping[str, Any] | None) -> list[tuple[str, dict[str, Any]]]:
    out: list[tuple[str, dict[str, Any]]] = []
    for source in _SOURCES:
        for panel_id, record in ((similarity or {}).get(source) or {}).items():
            if panel_id in _NOT_PANELS:
                continue
            if isinstance(record, list):
                record = record[0] if record and isinstance(record[0], dict) else {}
            if isinstance(record, dict) and record.get("profile_id"):
                out.append((str(panel_id), record))
    return out


def score(similarity: Mapping[str, Any] | None) -> AxisScore:
    """The axis score, or a stated reason there is none."""
    panels = _panels(similarity)
    if not panels:
        return AxisScore(None, 0, 0, 0.0, (), "the skin panels were not computed for this run")

    contributors: list[dict[str, Any]] = []
    weighted = 0.0
    weight = 0.0
    for panel_id, record in panels:
        index = record.get("score_0_100")
        coverage = record.get("coverage")
        if index is None or not coverage:
            continue
        # A panel that resolved none of its markers carries no weight, and a
        # panel that resolved all of them carries a full one.
        w = float(coverage)
        weighted += float(index) * w
        weight += w
        contributors.append({
            "panel_id": panel_id,
            "label": str(record.get("label") or panel_id),
            "index": round(float(index), 1),
            "coverage": round(w, 3),
            "markers_used": record.get("usable_count"),
            "markers_in_panel": record.get("panel_count"),
        })

    if weight <= 0:
        return AxisScore(
            None, len(panels), 0, 0.0, tuple(contributors),
            "none of the skin panels resolved enough of its marker organisms in this "
            "sample to carry a reading",
        )

    # The panels score "more like the disease cohort" upward; the axis runs
    # the healthy way up, like every other favourable reading here.
    leaning = weighted / weight
    contributors.sort(key=lambda c: -c["coverage"])
    return AxisScore(
        value=round(100.0 - leaning, 1),
        n_panels=len(panels),
        n_scored=len(contributors),
        total_coverage=weight,
        contributors=tuple(contributors),
    )


__all__ = ["METHOD", "AxisScore", "score"]
