"""Build a reference distribution for every functional reading.

A reading is a set of genes, and its value is their sum - the same rule
every panel in this reference declares. What was missing was something to
compare that sum against, so readings over more than one gene had no
position and the report printed "not measured" beside forty of them.

The cohort's own per-sample numbers make one. They are recovered from the
alignments still on disk, summed per reading per sample, and reduced with
the cohort's own percentile function over the cohort's own sample set, so
a reading's ladder is built exactly the way a panel's was.

Verifies itself first: the gene ladders rebuilt this way must reproduce
the shipped ones, or the arithmetic is not the arithmetic that made them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.cohort import _percentiles  # noqa: E402,PLC2701 - the reference's own rule

RANGES = REPO / "refs" / "reference_ranges.json"

#: Samples to leave out, so the reading ladders are built over exactly the
#: cohort the panel ladders were. Kept empty unless the two disagree: the
#: self-check below is what decides whether they do, and it refused the
#: build outright when a stale exclusion made this 99 against the shipped
#: 100.
EXCLUDE: set[str] = set()


def main() -> int:
    harvest = json.loads(Path("/tmp/cohort_gene_values.json").read_text())
    ranges = json.loads(RANGES.read_text())

    keep = [i for i, s in enumerate(harvest["samples"]) if s not in EXCLUDE]
    samples = [harvest["samples"][i] for i in keep]
    genes = {k: [v[i] for i in keep] for k, v in harvest["genes"].items()
             if len(v) == len(harvest["samples"])}
    print(f"{len(samples)} cohort samples, {len(genes)} genes")

    # --- self-check: rebuild the gene ladders and compare with the shipped
    exact = close = 0
    for key, values in genes.items():
        shipped = ranges["genes"].get(key)
        if not shipped:
            continue
        mine = _percentiles(values)
        theirs = {int(k): v for k, v in shipped["percentiles"].items()}
        gaps = [abs(mine[p] - theirs[p]) for p in mine if p in theirs]
        scale = max([abs(v) for v in theirs.values()] + [1e-9])
        if max(gaps) < 1e-6:
            exact += 1
        elif max(gaps) / scale < 0.01:
            close += 1
    print(f"  gene ladders reproduced exactly: {exact}, within 1%: {close}, "
          f"of {len(genes)}")
    if exact + close < 0.95 * len(genes):
        print("  REFUSING: this arithmetic does not reproduce the published ladders")
        return 1

    # --- the readings, from the gene sets each module declares
    readings = json.loads(Path("/tmp/reading_genes.json").read_text())
    built: dict[str, dict] = {}
    for reading_id, keys in readings.items():
        usable = [k for k in keys if k in genes]
        if not usable:
            continue
        per_sample = [sum(genes[k][i] for k in usable) for i in range(len(samples))]
        built[reading_id] = {
            "reading": reading_id,
            "entries": sorted(usable),
            "rule": "sum",
            "cohort_n": len(samples),
            "detected_in": sum(1 for v in per_sample if v > 0),
            "percentiles": {str(p): round(v, 6) for p, v in _percentiles(per_sample).items()},
            "mean": round(sum(per_sample) / len(per_sample), 6),
            "max": round(max(per_sample), 6),
        }
    print(f"  built {len(built)} reading distributions")

    ranges["readings"] = built
    RANGES.write_text(json.dumps(ranges, indent=2))
    print(f"wrote {len(built)} reading ladders into {RANGES.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
