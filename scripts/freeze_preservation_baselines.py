#!/usr/bin/env python3
"""Re-freeze the preservation baselines, with the reason on the record.

The baselines exist so that an unintended change to a legacy reading
fails a test. That makes re-freezing them a deliberate act, never a way
to make a red test green: run this only when every difference has been
traced to a change that was asked for, and say so in `--reason`, which
goes into the manifest beside the new hashes.

    scripts/freeze_preservation_baselines.py --reason "..." [--dry-run]

Without `--force` it refuses to freeze a sample whose run was
incomplete, because a baseline built from a run with a stage missing
would freeze the gap in.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from openbiota.completeness import audit  # noqa: E402
from openbiota.extension import preservation as P  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures" / "preservation"
RESULTS = REPO / "results"


def protected_subset(results: dict[str, Any]) -> dict[str, Any]:
    """Exactly the keys `preservation.compare` walks, and nothing else."""
    return {k: results[k] for k in P.PROTECTED_KEYS if k in results}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reason", required=True, help="why these baselines are being replaced")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="freeze even an incomplete run")
    ap.add_argument("samples", nargs="*", help="default: every sample already in the manifest")
    args = ap.parse_args()

    manifest_path = FIXTURES / "protected_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    samples = args.samples or sorted(manifest["samples"])

    changed: list[str] = []
    from openbiota.samples import local_name, published_name, results_file

    # frozen under the published identifier, read from the local results directory
    samples = sorted({published_name(s) for s in samples})
    for sample in samples:
        results_path = results_file(sample)
        if not results_path.is_file():
            print(f"  SKIP  {sample}: no results.json")
            continue
        results = json.loads(results_path.read_text())

        # A run that carries no completeness record was produced before the
        # gate existed, which means it also predates whatever change is
        # being frozen. Freezing it would pin the old numbers under a new
        # reason, and the next real run would fail against them.
        if not (results.get("run") or {}).get("completeness") and not args.force:
            print(f"  REFUSED  {sample}: this run predates the completeness gate; re-run it first")
            continue

        report = audit(results)
        if not report.complete and not args.force:
            print(f"  REFUSED  {sample}: run is incomplete ({', '.join(report.missing)})")
            continue

        results = P.anonymise(results, local=local_name(sample), published=sample)
        subset = protected_subset(results)
        blob = json.dumps(subset, indent=2, sort_keys=True).encode()
        digest = hashlib.sha256(blob).hexdigest()
        record = P.build_manifest(results, sample=sample)

        old = manifest["samples"].get(sample, {})
        if old.get("frozen_protected_sha256") == digest:
            print(f"  same  {sample}: unchanged")
            continue

        n_before = old.get("metric_id_count", 0)
        n_after = record["metric_id_count"]
        print(f"  FREEZE {sample}: {n_before} -> {n_after} metric ids, sha {digest[:12]}")
        changed.append(sample)
        if args.dry_run:
            continue

        (FIXTURES / f"{sample}.protected.json.gz").write_bytes(gzip.compress(blob, 9))
        manifest["samples"][sample] = {
            "baseline_version": record["baseline_version"],
            "frozen_file": f"{sample}.protected.json.gz",
            "frozen_protected_sha256": digest,
            "metric_id_count": n_after,
            "metric_ids": record["metric_ids"],
            "frozen_on": dt.date.today().isoformat(),
            "frozen_because": args.reason,
            "previous_protected_sha256": old.get("frozen_protected_sha256"),
        }

    if changed and not args.dry_run:
        manifest["last_reason"] = args.reason
        manifest["last_frozen_on"] = dt.date.today().isoformat()
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
        print(f"\nre-froze {len(changed)}: {', '.join(changed)}")
        print(f"reason recorded: {args.reason}")
    elif not changed:
        print("\nnothing to do")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
