"""Detection confidence tiers, calibrated against measured validation data.

A gene either clears the bar for a confirmed call or it does not, and the two
cases are reported differently. The thresholds below are not guesses — they are
the operating point at which validation against known ground truth produced
**zero false positives**.

Measured on 24 reference genomes across 5 synthetic communities
(``docs/validation_results.json``):

    operating point                     sens    spec    prec    FP
    any fragment at all                100.0%   93.5%   89.2%    4
    >= 2 fragments                     100.0%   96.8%   94.3%    2
    >= 2 fragments AND >= 75% identity  97.0%  100.0%  100.0%    0   <- CONFIRMED

The single call lost at the confirmed tier is a carrier present at 0.2% of
cells, which is the hardest case in the whole validation set. Trading that for
the elimination of every false positive is the right way round: a detection
that is reported is one you can rely on, and the borderline cases are still
shown rather than discarded.

Why these two criteria and not others:

* **Fragment count.** A single read is not evidence. Both single-fragment
  false positives in validation disappear at two, and no true call is affected
  — the lowest true positive had 64 fragments.
* **Mean identity.** Every false positive in the validation set fell below
  73% mean translated identity, while genuine calls clustered at 81-99% with a
  median of 98%. A read from an organism actually carrying the gene matches
  its reference closely; matches in the 60s are landing on distant relatives
  within the same protein fold.
"""

from __future__ import annotations

from typing import Final

#: A detection needs more than one read.
MIN_FRAGMENTS: Final = 2

#: Mean translated identity floor for a confirmed call.
MIN_MEAN_IDENTITY: Final = 75.0

CONFIRMED: Final = "confirmed"
PROVISIONAL: Final = "provisional"
NOT_DETECTED: Final = "not detected"
INDETERMINATE: Final = "indeterminate"

#: Human-readable explanation per tier, used in reports.
EXPLANATION: Final = {
    CONFIRMED: (
        "Meets the criteria at which validation against known samples produced no false "
        "positives."
    ),
    PROVISIONAL: (
        "Detected, but below the threshold at which detections were error-free in "
        "validation. Treat as a possible signal rather than a confirmed one."
    ),
    NOT_DETECTED: "No reads matched this pathway.",
    INDETERMINATE: "Could not be measured — no calibration gene was detected.",
}


def classify_confidence(
    *,
    fragments: int,
    mean_identity: float | None,
    indeterminate: bool = False,
) -> str:
    """Assign a detection confidence tier."""
    if indeterminate:
        return INDETERMINATE
    if fragments <= 0:
        return NOT_DETECTED
    if fragments >= MIN_FRAGMENTS and (mean_identity or 0.0) >= MIN_MEAN_IDENTITY:
        return CONFIRMED
    return PROVISIONAL


def reasons_below_confirmed(fragments: int, mean_identity: float | None) -> list[str]:
    """Why a detection did not reach the confirmed tier."""
    out: list[str] = []
    if 0 < fragments < MIN_FRAGMENTS:
        out.append(f"only {fragments} matching fragment (at least {MIN_FRAGMENTS} required)")
    if fragments > 0 and (mean_identity or 0.0) < MIN_MEAN_IDENTITY:
        out.append(
            f"matches average {mean_identity:.0f}% identity "
            f"(at least {MIN_MEAN_IDENTITY:.0f}% required), so the reads are matching "
            "distant relatives rather than the gene itself"
            if mean_identity is not None
            else "identity could not be determined"
        )
    return out
