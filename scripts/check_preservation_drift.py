#!/usr/bin/env python3
"""Is every protected difference one the detection expansion is allowed to make?

The preservation baselines (spec 0.8.3) freeze every leaf under the
protected result keys. Spec 0.8.4 deliberately changes some of those
leaves: the organism inventory grows, the group members' detection
provenance and composition shares are read from the merged inventory, the
metabolite-driver rows learn which drivers are actually in the sample. The
calibrated numbers - percentiles, scores, indices, bands, prevalences, the
pathogen and functional panels - must not move.

This script compares every sample with its baseline and sorts each
difference into `allowed` (a detection-derived leaf, or a new leaf) or
`violation` (anything else). Exit 0 only when there are no violations;
that is the condition under which `freeze_preservation_baselines.py` may
be run, with this script's summary as the reason.

    scripts/check_preservation_drift.py [--json out.json]
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.extension import preservation as P  # noqa: E402

#: Leaves the expansion may change, by regex on the difference path.
ALLOWED = [
    r"^organism_inventory(\.|\[|$)",
    r"^organism_verdicts(\.|\[|$)",
    r"^profile_similarity\.community\.groups\[\d+\]\.members\[\d+\]\.(detected_by|expansion_only|percent|detected|sample_percent)$",
    r"^profile_similarity\.community\.groups\[\d+\]\.(n_detected_expansion_only|n_detected|scoring_percent|percent|detected|n_resolvable|members)$",
    r"^profile_similarity\.community\.(species|n_species_detected|n_genera_detected|notes)(\.|\[|$)",
    r"^metabolite_drivers\.[^.]+\.top\[\d+\]\.(in_sample|genus_in_sample|sample_class|sample_percent|percent)$",
    r"^metabolite_drivers\.[^.]+\.(n_in_sample|in_sample|sample_percent|coverage)$",
    r"^findings_and_evidence(\.|\[|$)",  # cards follow the inventory
    # which Candida-suppressing bacteria the inventory holds, and at what share: detection-derived context
    r"^mycobiome\.colonisation_resistance\.bacteria\[\d+\]\.(percent|present|percentile)$",
    r"^mycobiome\.colonisation_resistance\.relevant$",
    r"^extended_catalogue\.(sam|mapout|elapsed_s|command)$",
]
#: Words that mark a calibrated number; a changed leaf whose path holds one
#: is a violation even inside an otherwise allowed object.
NEVER = re.compile(r"percentile|score|z_|robust_z|index|band|prevalence|combined|threshold|median|mad\b", re.I)


#: Detection-derived context that happens to be called a percentile: the
#: inventory's own placement of an organism, copied into another block.
DERIVED_PERCENTILES = re.compile(r"^mycobiome\.colonisation_resistance\.bacteria\[\d+\]\.percentile$")

#: Rule changes made deliberately in the v0.8.4 release, each documented in
#: docs/EXPANDED_DETECTION.md and docs/OUTPUT.md and covered by tests:
#:   - the typical carrier's level is the interpolated median of the carrier
#:     readings (community.carrier_median), so it and the midrank percentile
#:     describe one distribution; the deviation follows it;
#:   - a scoring-cohort rank that conflicts with the composition share falls
#:     back to a lane cohort whose reading agrees with it (inventory.build);
#:   - the pathogen screen's bacterial calls are reconciled with the organism
#:     inventory (pathogens.reconcile): statuses, counts and statements move.
#: Nothing else calibrated may change; a path outside this list still fails.
RELEASE_RULE_CHANGES = [
    r"^profile_similarity\.community\.(groups\[\d+\]\.(members\[\d+\]\.)?|species\[\d+\]\.)(reference_median|deviation_percent|carrier_percentile)$",
    r"^organism_verdicts\.verdicts\[\d+\]\.percentile$",
    r"^pathogens\.results\[\d+\]\.(plain_statement|reason_codes|display_qualifier|display_status|sequence_status|species_resolution|counts_as_pathogen|report_tier|inventory_\w+)$",
    r"^pathogens\.counts\.(pathogen_count|supported|attention|opportunists|uncertain)$",
    r"^pathogens\.inventory_reconciliation(\.|\[|$)",
]


def classify(path: str, baseline, candidate) -> str:  # noqa: ARG001 - the candidate value is reported by the caller
    if baseline is None:
        return "allowed"  # an addition (or a new leaf that is null): no calibrated number changed
    if DERIVED_PERCENTILES.search(path):
        return "allowed"
    if any(re.search(rx, path) for rx in RELEASE_RULE_CHANGES):
        return "allowed"
    if NEVER.search(path.rsplit(".", 1)[-1]) and not path.startswith("organism_inventory"):
        return "violation"
    return "allowed" if any(re.search(rx, path) for rx in ALLOWED) else "violation"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    fixtures = REPO / "tests" / "fixtures" / "preservation"
    manifest_path = fixtures / "protected_manifest.json"
    if not manifest_path.is_file():
        print("no preservation manifest found")
        return 2
    manifest = json.loads(manifest_path.read_text())
    out = {"samples": {}, "violations": 0}
    for sample in sorted(manifest["samples"]):
        from openbiota.samples import local_name, results_file

        base = fixtures / f"{sample}.protected.json.gz"
        current = results_file(sample)
        if not base.is_file() or not current.is_file():
            out["samples"][sample] = {"skipped": "missing baseline or results"}
            continue
        baseline = json.loads(gzip.decompress(base.read_bytes()).decode("utf-8"))
        live = P.anonymise(json.loads(current.read_text(encoding="utf-8")), local=local_name(sample), published=sample)
        report = P.compare(baseline, live)
        kinds = Counter()
        violations = []
        for d in report.differences:
            k = classify(d.path, d.baseline, d.candidate)
            kinds[k] += 1
            if k == "violation":
                violations.append({"path": d.path, "baseline": str(d.baseline)[:80], "candidate": str(d.candidate)[:80]})
        out["samples"][sample] = {"compared_leaves": report.compared_leaves, "differences": len(report.differences),
                                  "allowed": kinds["allowed"], "violations": violations}
        out["violations"] += len(violations)
        print(f"{sample:16} leaves {report.compared_leaves:>7,}  differences {len(report.differences):>5}  "
              f"allowed {kinds['allowed']:>5}  violations {len(violations)}")
        for v in violations[:10]:
            print(f"    VIOLATION {v['path']}: {v['baseline']} -> {v['candidate']}")
    if args.json:
        args.json.write_text(json.dumps(out, indent=1))
    print("OK: every protected difference is detection-derived" if out["violations"] == 0
          else f"{out['violations']} calibrated value(s) changed; do not re-freeze")
    return 0 if out["violations"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
