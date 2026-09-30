#!/usr/bin/env python3
"""Show where a running report has got to, and what is left.

`openbiota run --quiet` writes its log at the end, so while a run is going
there is nothing to read: the only evidence is which stage directories
have been written and when. This reads those, in the order the run uses
them, and prints a line per stage with the time it took.

    scripts/progress.py SAMPLE5_A05          once
    scripts/progress.py SAMPLE5_A05 --watch  every 20 seconds

Stage durations are taken from the last completed run of any sample, so
the estimate is this machine's own timing rather than a guess. A stage
that is slower from cold - the translated search, the fungal alignment,
the metabolic simulation - says so, because a cached run is not evidence
about an uncached one.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: The stages in the order the run performs them, with the directory each
#: writes into. A stage with no directory is only visible in the log.
STAGES: tuple[tuple[str, str], ...] = (
    ("preprocess", "read preprocessing (fastp)"),
    ("alignments", "translated search (diamond blastx)"),
    ("taxonomy", "taxonomic engine (MetaPhlAn)"),
    ("genome", "whole-genome profile"),
    ("pathogens", "pathogen screen"),
    ("mycobiome", "mycobiome (fungal alignment)"),
    ("strain", "strain resolution"),
)

#: Stages that are much slower the first time, and why.
COLD = {
    "alignments": "no cached alignment for this sample",
    "mycobiome": "fungal alignment runs from scratch",
    "taxonomy": "MetaPhlAn maps every read",
}


def _mtime(path: Path) -> float | None:
    try:
        return max(
            (p.stat().st_mtime for p in path.rglob("*") if p.is_file()), default=None
        )
    except OSError:
        return None


def _reference_timings() -> dict[str, float]:
    """How long each stage took on the most recent completed run."""
    best: dict[str, float] = {}
    for log in sorted(REPO.glob("results/*/run.log"), key=lambda p: -p.stat().st_mtime):
        marks = re.findall(r"^\[ *([0-9.]+)s\] . (.+?)(?: done in| reused|$)",
                           log.read_text(errors="ignore"), re.M)
        if not marks:
            continue
        previous = 0.0
        for at, name in marks:
            best.setdefault(name.strip(), max(0.0, float(at) - previous))
            previous = float(at)
        break
    return best


def report(sample: str) -> bool:
    """Print the current state. Returns True when the report is finished."""
    out = REPO / "results" / sample
    pdf = out / f"{sample}_report.pdf"
    if not out.is_dir():
        print(f"no run directory for {sample}")
        return True

    started = min(
        (p.stat().st_mtime for p in out.rglob("*") if p.is_file()), default=time.time()
    )
    done: list[tuple[str, float]] = []
    for directory, label in STAGES:
        at = _mtime(out / directory)
        if at is not None:
            done.append((label, at))
    done.sort(key=lambda row: row[1])

    print(f"\n{sample} — {time.strftime('%H:%M:%S')}, "
          f"{(time.time() - started) / 60:.0f} min elapsed")
    previous = started
    for label, at in done:
        print(f"  ✓ {label:38} {(at - previous) / 60:5.1f} min")
        previous = at

    if pdf.is_file() and pdf.stat().st_size and pdf.stat().st_mtime >= previous:
        print(f"  ✓ {'report written':38} {pdf.stat().st_size / 1e6:5.1f} MB")
        print("\nfinished.")
        return True

    last = done[-1][0] if done else "starting"
    ago = (time.time() - previous) / 60
    print(f"  · in progress: after {last}, {ago:.0f} min so far")
    remaining = [lab for _, lab in STAGES if lab not in {d[0] for d in done}]
    if remaining:
        print(f"    still to come: {', '.join(remaining)}")
    print("    then: findings, report extension (metabolic simulation), PDF")
    if not list((REPO / "refs" / "micom" / "lp").glob("*")) if (
        REPO / "refs" / "micom" / "lp"
    ).is_dir() else False:
        print("    note: the simulation cache is cold, so that stage solves from scratch")
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample")
    parser.add_argument("--watch", action="store_true", help="refresh every 20 seconds")
    args = parser.parse_args()
    while True:
        if report(args.sample) or not args.watch:
            return 0
        time.sleep(20)


if __name__ == "__main__":
    sys.exit(main())
