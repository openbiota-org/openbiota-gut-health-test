"""Feature-specific detection power (spec v04 §4.2–4.3).

"Not detected" is only a finding when the sample had a real chance of
detecting the organism. That chance depends on how many usable reads there
were and on how abundant the organism is when it is present. A species that
carriers hold at 2% of the community is found in a shallow sample; one they
hold at 0.005% is not, and its absence from that sample means nothing.

The model here is deliberately simple and is calibrated against the profiler
that is actually run rather than assumed:

    P(detected | N, a) = Φ((log(N·a) − log τ) / σ)

where ``N`` is the number of usable microbial reads, ``a`` the organism's
relative abundance (fraction) and ``τ`` the effective "marker reads needed"
at which detection is a coin toss. ``τ`` and ``σ`` are fitted by profiling
one deep sample at a ladder of depths and recording which species drop out
(``openbiota calibrate-detection``). The fitted values live in
``refs/detection_model.json`` with their calibration record; until that
file exists a conservative prior is used and every output says so.

The detection power for a feature at this sample's depth is then

    D_f(N) = mean over the reference's positive abundances a_i of P(detected | N, a_i)

— the probability the organism would have been seen here if it were present
at a typical carrier's abundance. A nondetection is *qualified* only when
``D_f(N) >= QUALIFIED_POWER``; below that it is ``indeterminate``.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

#: Minimum detection power for a nondetection to count as evidence of absence.
QUALIFIED_POWER: Final = 0.95

#: Detection states every feature receives (spec §4.2).
DETECTED: Final = "detected"
NOT_DETECTED: Final = "not_detected"
INDETERMINATE: Final = "indeterminate"
NOT_MEASURED: Final = "not_measured"

MODEL_ID: Final = "feature_depth_detection_curve_v1"

#: Prior used before calibration. MetaPhlAn 3 keeps ~100 kb of markers per
#: species (~2.5–5% of a genome) and needs reads across a good share of them,
#: so a few hundred reads-per-unit-abundance is where detection turns over.
PRIOR_LOG_TAU: Final = math.log(500.0)
PRIOR_SIGMA: Final = 0.7


def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


@dataclass(frozen=True, slots=True)
class DetectionModel:
    """Frozen depth-detection curve for one profiler/database."""

    log_tau: float
    sigma: float
    profiler: str
    database: str
    calibrated: bool
    calibration: dict[str, Any] = field(default_factory=dict)
    model_id: str = MODEL_ID

    @property
    def digest(self) -> str:
        payload = json.dumps(
            {"model_id": self.model_id, "log_tau": self.log_tau, "sigma": self.sigma,
             "profiler": self.profiler, "database": self.database},
            sort_keys=True,
        )
        return "sha256:" + hashlib.sha256(payload.encode()).hexdigest()[:16]

    def p_detect(self, n_reads: float, abundance_fraction: float) -> float:
        """Probability of detecting an organism at ``abundance_fraction`` with ``n_reads``."""
        if n_reads <= 0 or abundance_fraction <= 0:
            return 0.0
        return _phi((math.log(n_reads * abundance_fraction) - self.log_tau) / self.sigma)

    def power(self, n_reads: float, positive_abundances: Sequence[float]) -> float | None:
        """``D_f(N)``: mean detection probability over the carriers' abundances (fractions)."""
        if not positive_abundances or n_reads <= 0:
            return None
        return sum(self.p_detect(n_reads, a) for a in positive_abundances) / len(positive_abundances)

    def lod95(self, n_reads: float) -> float | None:
        """Abundance fraction detected with 95% probability at this depth."""
        if n_reads <= 0:
            return None
        return math.exp(self.log_tau + 1.6449 * self.sigma) / n_reads

    def to_json(self) -> dict[str, Any]:
        return {
            "model_id": self.model_id,
            "digest": self.digest,
            "log_tau": round(self.log_tau, 4),
            "tau_reads_per_unit_abundance": round(math.exp(self.log_tau), 1),
            "sigma": round(self.sigma, 4),
            "profiler": self.profiler,
            "database": self.database,
            "calibrated": self.calibrated,
            "calibration": self.calibration,
            "qualified_power": QUALIFIED_POWER,
            "form": "P(detected | N, a) = Phi((ln(N*a) - ln tau) / sigma); D_f(N) = mean_i P(detected | N, a_i) over reference carriers",
        }

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_json(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> DetectionModel:
        data = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            log_tau=float(data["log_tau"]), sigma=float(data["sigma"]),
            profiler=str(data.get("profiler", "")), database=str(data.get("database", "")),
            calibrated=bool(data.get("calibrated", False)),
            calibration=dict(data.get("calibration", {})),
            model_id=str(data.get("model_id", MODEL_ID)),
        )

    @classmethod
    def prior(cls, profiler: str = "MetaPhlAn 3", database: str = "") -> DetectionModel:
        # `note` is printed in the report, so it reads as a sentence to a
        # person and not as an instruction to an operator. The machine-
        # readable state is `calibrated`, and the command to fix it belongs
        # in the documentation, not on someone's results.
        return cls(log_tau=PRIOR_LOG_TAU, sigma=PRIOR_SIGMA, profiler=profiler,
                   database=database, calibrated=False,
                   calibration={
                       "note": (
                           "estimated from a conservative general curve rather than "
                           "measured on this laboratory's own data, so detection "
                           "limits are approximate"
                       ),
                       "status": "prior",
                   })


def load_detection_model(refs_dir: Path, *, profiler: str = "MetaPhlAn 3", database: str = "") -> DetectionModel:
    path = refs_dir / "detection_model.json"
    if path.is_file():
        try:
            return DetectionModel.load(path)
        except (OSError, ValueError, KeyError):
            pass
    return DetectionModel.prior(profiler, database)


# --------------------------------------------------------------------------- #
# calibration
# --------------------------------------------------------------------------- #


def fit_detection_curve(
    observations: Sequence[tuple[float, float, bool]],
) -> tuple[float, float, dict[str, Any]]:
    """Fit ``(log_tau, sigma)`` by maximum likelihood.

    ``observations`` are ``(n_reads, abundance_fraction, detected)`` triples
    from a depth ladder: for every species called at full depth, whether it
    was still called at each shallower depth. Grid search over a generous
    box, refined once; the likelihood is smooth and unimodal in practice.
    """
    xs = [math.log(n * a) for n, a, _ in observations if n > 0 and a > 0]
    ys = [d for n, a, d in observations if n > 0 and a > 0]
    if len(xs) < 10 or all(ys) or not any(ys):
        raise ValueError("calibration needs mixed detected/undetected observations across depths")

    def nll(mu: float, sigma: float) -> float:
        total = 0.0
        for x, y in zip(xs, ys, strict=True):
            p = min(max(_phi((x - mu) / sigma), 1e-9), 1 - 1e-9)
            total -= math.log(p if y else 1 - p)
        return total

    best = (float("inf"), PRIOR_LOG_TAU, PRIOR_SIGMA)
    lo, hi = min(xs), max(xs)
    for step in (0.25, 0.05, 0.01):
        mu_grid = [best[1] + k * step for k in range(-20, 21)] if step < 0.25 else [
            lo + (hi - lo) * k / 60 for k in range(61)
        ]
        sig_grid = [max(0.05, best[2] + k * step) for k in range(-20, 21)] if step < 0.25 else [
            0.1 + 0.1 * k for k in range(30)
        ]
        for mu in mu_grid:
            for sigma in sig_grid:
                value = nll(mu, sigma)
                if value < best[0]:
                    best = (value, mu, sigma)
    _, mu, sigma = best
    n_det = sum(1 for y in ys if y)
    return mu, sigma, {
        "n_observations": len(xs), "n_detected": n_det, "n_undetected": len(xs) - n_det,
        "negative_log_likelihood": round(best[0], 2),
    }


def calibrate_from_depth_ladder(
    *,
    full_profile: dict[str, float],
    ladder: Sequence[tuple[int, dict[str, float]]],
    full_reads: int,
    profiler: str,
    database: str,
    sample: str,
) -> DetectionModel:
    """Build a model from a full-depth species profile and shallower re-profiles.

    ``full_profile`` maps species -> percent at full depth; each ladder entry
    is ``(n_reads, species -> percent)`` for one subsample. Species called at
    full depth are the truth set; their full-depth abundance is the ``a``.
    """
    truth = {s: v / 100.0 for s, v in full_profile.items() if v > 0}
    observations: list[tuple[float, float, bool]] = []
    for n_reads, profile in ladder:
        seen = {s for s, v in profile.items() if v > 0}
        for species, a in truth.items():
            observations.append((float(n_reads), a, species in seen))
    log_tau, sigma, fit = fit_detection_curve(observations)
    return DetectionModel(
        log_tau=log_tau, sigma=sigma, profiler=profiler, database=database, calibrated=True,
        calibration={
            "method": "empirical depth ladder, maximum-likelihood probit in ln(N*a)",
            "sample": sample,
            "full_depth_reads": full_reads,
            "n_species_at_full_depth": len(truth),
            "ladder_reads": [n for n, _ in ladder],
            "fitted_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            **fit,
        },
    )


__all__ = [
    "DETECTED",
    "INDETERMINATE",
    "MODEL_ID",
    "NOT_DETECTED",
    "NOT_MEASURED",
    "QUALIFIED_POWER",
    "DetectionModel",
    "calibrate_from_depth_ladder",
    "fit_detection_curve",
    "load_detection_model",
]
