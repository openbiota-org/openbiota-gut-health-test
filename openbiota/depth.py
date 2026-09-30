"""Saturation analysis: is the sequencing depth sufficient?

The pipeline always reads the whole input — every fragment in both mate files.
The open question is whether the *sequencer* sampled the stool deeply enough
that a rarer organism could not have been missed.

That is answerable rather than arguable. Screen the same sample at increasing
depths and watch the normalised values. If they stop moving well before full
depth, depth is not the limiting factor and more sequencing would buy nothing.
If they are still drifting at full depth, deeper sequencing genuinely would
change the answer.

Two quantities are tracked per depth:

* the normalised value, which should converge to a plateau;
* the fragment count, which grows linearly with depth by construction and is
  what sets the counting precision.

Convergence is judged on the relative change between the last two depths,
compared against the Poisson error expected from the fragment count at full
depth. A value that moves by less than its own counting noise has converged —
there is no signal left to recover.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from openbiota.logging_util import Reporter, human_duration

#: Depths screened, as a fraction of the full input. The full depth is always
#: appended so the series ends at the real answer.
DEFAULT_FRACTIONS: tuple[float, ...] = (0.0625, 0.125, 0.25, 0.5, 1.0)

#: A panel needs at least this many fragments at full depth for convergence to
#: be meaningful. Below it, the value is noise-dominated at every depth and the
#: question "has it converged" does not have a useful answer.
MIN_FRAGMENTS_FOR_VERDICT = 30


@dataclass(frozen=True, slots=True)
class DepthPoint:
    """One panel measured at one depth."""

    panel: str
    read_pairs: int
    fragments: int
    rpob_fragments: int
    copies_per_100: float | None

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PanelConvergence:
    """Whether one panel's value has stopped moving with depth."""

    panel: str
    full_depth_value: float | None
    full_depth_fragments: int
    half_depth_value: float | None
    relative_change: float | None
    poisson_error: float | None
    converged: bool
    verdict: str

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class DepthReport:
    sample: str
    full_read_pairs: int
    fractions: list[float]
    points: list[DepthPoint] = field(default_factory=list)
    convergence: list[PanelConvergence] = field(default_factory=list)
    generated: str = ""
    elapsed_s: float = 0.0

    @property
    def all_converged(self) -> bool:
        evaluable = [c for c in self.convergence if c.verdict != "too few fragments to judge"]
        return bool(evaluable) and all(c.converged for c in evaluable)

    def summary_sentence(self) -> str:
        evaluable = [c for c in self.convergence if c.verdict != "too few fragments to judge"]
        converged = [c for c in evaluable if c.converged]
        if not evaluable:
            return "Too few fragments at every depth to judge convergence."
        if len(converged) == len(evaluable):
            return (
                f"All {len(evaluable)} measurable pathways had stopped changing before full "
                "depth was reached, so the sequencing was deep enough — more reads would not "
                "change these results."
            )
        drifting = [c.panel for c in evaluable if not c.converged]
        return (
            f"{len(converged)} of {len(evaluable)} measurable pathways had converged by full "
            f"depth. Still moving: {', '.join(drifting)}. Deeper sequencing would sharpen "
            "those."
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "sample": self.sample,
            "full_read_pairs": self.full_read_pairs,
            "fractions": self.fractions,
            "generated": self.generated,
            "elapsed_s": round(self.elapsed_s, 1),
            "all_converged": self.all_converged,
            "summary": self.summary_sentence(),
            "method": (
                "The same sample screened at increasing depths. A pathway counts as converged "
                "when its value changes between half and full depth by less than the Poisson "
                "counting error implied by its own fragment count — i.e. by less than the "
                "noise floor, leaving no recoverable signal."
            ),
            "points": [p.to_json() for p in self.points],
            "convergence": [c.to_json() for c in self.convergence],
        }

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> dict[str, Any]:
        return json.loads(path.read_text(encoding="utf-8"))


def poisson_relative_error(count: int) -> float | None:
    """Relative standard error of a Poisson count: 1/sqrt(n)."""
    return 1.0 / math.sqrt(count) if count > 0 else None


def assess_convergence(points: Sequence[DepthPoint], panel: str) -> PanelConvergence:
    """Compare the last two depths against the counting noise at full depth."""
    series = sorted(
        (p for p in points if p.panel == panel), key=lambda p: p.read_pairs
    )
    if len(series) < 2:
        return PanelConvergence(
            panel=panel,
            full_depth_value=series[-1].copies_per_100 if series else None,
            full_depth_fragments=series[-1].fragments if series else 0,
            half_depth_value=None,
            relative_change=None,
            poisson_error=None,
            converged=False,
            verdict="not enough depths screened",
        )

    full, half = series[-1], series[-2]
    error = poisson_relative_error(full.fragments)

    if full.fragments < MIN_FRAGMENTS_FOR_VERDICT or error is None:
        return PanelConvergence(
            panel=panel,
            full_depth_value=full.copies_per_100,
            full_depth_fragments=full.fragments,
            half_depth_value=half.copies_per_100,
            relative_change=None,
            poisson_error=error,
            converged=False,
            verdict="too few fragments to judge",
        )

    if full.copies_per_100 is None or half.copies_per_100 is None or full.copies_per_100 == 0:
        return PanelConvergence(
            panel=panel,
            full_depth_value=full.copies_per_100,
            full_depth_fragments=full.fragments,
            half_depth_value=half.copies_per_100,
            relative_change=None,
            poisson_error=error,
            converged=False,
            verdict="not measurable at one or both depths",
        )

    change = abs(full.copies_per_100 - half.copies_per_100) / full.copies_per_100
    # Comparing two independent counts, so the noise on their difference is
    # the quadrature sum; the half-depth point carries sqrt(2) more error.
    threshold = error * math.sqrt(3.0)
    converged = change <= threshold

    return PanelConvergence(
        panel=panel,
        full_depth_value=full.copies_per_100,
        full_depth_fragments=full.fragments,
        half_depth_value=half.copies_per_100,
        relative_change=change,
        poisson_error=error,
        converged=converged,
        verdict=(
            "converged — change is within counting noise"
            if converged
            else "still changing beyond counting noise"
        ),
    )


def print_depth_summary(report: DepthReport, reporter: Reporter | None = None) -> None:
    depths = sorted({p.read_pairs for p in report.points})
    panels = sorted({p.panel for p in report.points})
    by_key = {(p.panel, p.read_pairs): p for p in report.points}

    print()
    print("=" * 88)
    print(" SEQUENCING DEPTH CHECK — does reading more of the sample change the answer?")
    print("=" * 88)
    print()
    print(f"  sample                {report.sample}")
    print(f"  full depth            {report.full_read_pairs:,} read pairs (100% of the input)")
    print(f"  depths screened       {', '.join(f'{d:,}' for d in depths)}")
    print(f"  elapsed               {human_duration(report.elapsed_s)}")
    print()

    header = "  " + f"{'pathway':<12}" + "".join(f"{d // 1000:>9}k" for d in depths)
    print("  copies per 100 genomes at each depth")
    print(header)
    print("  " + "-" * (12 + 10 * len(depths)))
    for panel in panels:
        row = f"  {panel:<12}"
        for depth in depths:
            point = by_key.get((panel, depth))
            value = point.copies_per_100 if point else None
            row += f"{('—' if value is None else f'{value:.2f}'):>10}"
        print(row)

    print()
    print("  fragments matched at each depth")
    print(header)
    print("  " + "-" * (12 + 10 * len(depths)))
    for panel in panels:
        row = f"  {panel:<12}"
        for depth in depths:
            point = by_key.get((panel, depth))
            row += f"{(point.fragments if point else 0):>10,}"
        print(row)

    print()
    print(f"  {'pathway':<12}{'value':>9}{'change':>10}{'noise':>9}  verdict")
    print("  " + "-" * 74)
    for c in sorted(report.convergence, key=lambda c: c.panel):
        change = "—" if c.relative_change is None else f"{c.relative_change:.1%}"
        noise = "—" if c.poisson_error is None else f"±{c.poisson_error:.1%}"
        value = "—" if c.full_depth_value is None else f"{c.full_depth_value:.2f}"
        print(f"  {c.panel:<12}{value:>9}{change:>10}{noise:>9}  {c.verdict}")

    print()
    print("  " + report.summary_sentence())
    print()
    if reporter is not None:
        reporter.record("    " + report.summary_sentence())
