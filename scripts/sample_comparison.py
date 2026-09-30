#!/usr/bin/env python3
"""Per-sample comparison of what the expanded detection changed (spec 0.8.4 §7).

For every regenerated sample this reads results.json and reports, side by
side: organisms the installed baseline saw, organisms now supported that
it never saw, renamed or deduplicated existing organisms, strain
placements added, calls rejected by competitive confirmation with their
reasons, and unresolved candidates (provisional, pending a genome, or
complexes). Samples of unknown truth show practical change, not accuracy;
the benchmark measures accuracy.

Writes results/EXPANSION_COMPARISON.json and results/EXPANSION_COMPARISON.md.

    scripts/sample_comparison.py [--results results]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

BASELINE_LANES = {"scoring", "extended", "genome"}


def one(results_path: Path) -> dict | None:
    try:
        r = json.loads(results_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    inv = r.get("organism_inventory") or {}
    orgs = inv.get("organisms") or []
    if not orgs:
        return None
    det = r.get("detection") or {}
    conf = det.get("confirmation") or {}
    verdicts = conf.get("verdicts") or []
    strain = r.get("strain_analysis") or {}
    lane_status = det.get("lane_status") or {}
    baseline_seen = [o for o in orgs if set(o.get("detected_by") or []) & BASELINE_LANES]
    supported = [o for o in orgs if o.get("status") == "supported"]
    newly = [o for o in orgs if o.get("incremental_gain") == "expansion" and o.get("status") == "supported"]
    provisional = [o for o in orgs if o.get("status") == "provisional"]
    complexes = [o for o in orgs if o.get("count_category") == "unresolved_complex"]
    renamed = [o for o in orgs if o.get("formerly")]
    placed = [c for c in (strain.get("clades") or []) if c.get("status") == "placed"]
    rejected = inv.get("rejected") or []
    pending = [v for v in verdicts if v.get("status") == "pending"]
    counts = inv.get("counts") or {}
    return {
        "sample": r.get("sample") or results_path.parent.name,
        "lanes_ran": sorted(k for k, v in lane_status.items() if v.get("state") == "ran"),
        "lanes_not_run": sorted(k for k, v in lane_status.items() if v.get("state") != "ran"),
        "organisms_total": len(orgs),
        "organisms_supported": len(supported),
        "organisms_seen_by_installed_baseline": len(baseline_seen),
        "newly_supported": [
            {"organism": o.get("gtdb") or o.get("species"), "reading_percent": o.get("secondary_percent") or o.get("percent"),
             "detected_by": o.get("detected_by"), "basis": o.get("confidence_basis"), "native_ids": o.get("native_ids")}
            for o in newly],
        "renamed_or_deduplicated": [{"organism": o.get("gtdb") or o.get("species"), "formerly": o.get("formerly")} for o in renamed],
        "strain_placements_added": [
            {"sgb": c.get("sgb"), "nearest_reference": c.get("nearest_reference"), "distance": c.get("distance_to_nearest"),
             "n_markers": c.get("n_markers")} for c in placed],
        "strain_unresolved": strain.get("n_unresolved"),
        "strain_deferred": len(strain.get("deferred") or []),
        "rejected_calls": [{"organism": x.get("gtdb") or x.get("species") or x.get("organism"),
                            "reason": x.get("confidence_basis") or x.get("reason"), "detected_by": x.get("detected_by")}
                           for x in rejected],
        "unresolved_candidates": {
            "provisional": [{"organism": o.get("gtdb") or o.get("species"), "detected_by": o.get("detected_by"), "reason": o.get("confidence_basis")}
                            for o in provisional],
            "pending_confirmation_genome": [{"organism": v.get("organism"), "reason": v.get("reason")} for v in pending],
            "complexes": [{"organism": o.get("gtdb") or o.get("species")} for o in complexes],
        },
        "count_categories": {k: v for k, v in counts.items() if k.startswith("category_") or k in ("supported", "provisional", "ambiguous", "rejected")},
        "percentile_populations": _populations(orgs),
        "assembly": {k: (r.get("assembly") or {}).get(k) for k in ("status", "reads") if (r.get("assembly") or {}).get(k) is not None},
    }


def _populations(orgs: list[dict]) -> dict[str, int]:
    out: dict[str, int] = {}
    for o in orgs:
        if o.get("percentile") is not None:
            k = o.get("percentile_source") or "scoring cohort"
            out[k] = out.get(k, 0) + 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=REPO / "results")
    args = ap.parse_args()
    rows = []
    for p in sorted(args.results.glob("*/results.json")):
        rec = one(p)
        if rec:
            rows.append(rec)
    out = {"written_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"), "n_samples": len(rows), "samples": rows,
           "meaning": ("Unknown-truth specimens: this shows what the expansion changed in each report, not sensitivity or "
                       "accuracy. A specimen need not gain species for the implementation to be correct.")}
    (args.results / "EXPANSION_COMPARISON.json").write_text(json.dumps(out, indent=1) + "\n")
    md = [f"# Expanded detection: per-sample comparison ({len(rows)} samples)", "",
          "| sample | organisms | supported | seen by installed baseline | newly supported | renamed | strains placed | rejected | provisional | pending genome | complexes |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s in rows:
        u = s["unresolved_candidates"]
        md.append(f"| {s['sample']} | {s['organisms_total']} | {s['organisms_supported']} | {s['organisms_seen_by_installed_baseline']} | "
                  f"{len(s['newly_supported'])} | {len(s['renamed_or_deduplicated'])} | {len(s['strain_placements_added'])} | "
                  f"{len(s['rejected_calls'])} | {len(u['provisional'])} | {len(u['pending_confirmation_genome'])} | {len(u['complexes'])} |")
    for s in rows:
        md += ["", f"## {s['sample']}", "",
               f"lanes ran: {', '.join(s['lanes_ran']) or 'none'}" + (f"; not run: {', '.join(s['lanes_not_run'])}" if s["lanes_not_run"] else ""),
               f"percentile populations: {s['percentile_populations']}", ""]
        if s["newly_supported"]:
            md += ["**Newly supported** (installed baseline never saw them):", ""]
            md += [f"- {n['organism']} — {n['reading_percent']}% by {', '.join(n['detected_by'] or [])}; {n['basis']}" for n in s["newly_supported"]]
            md.append("")
        if s["renamed_or_deduplicated"]:
            md += ["**Renamed / deduplicated:** " + "; ".join(f"{x['organism']} (formerly {x['formerly']})" for x in s["renamed_or_deduplicated"][:40]), ""]
        if s["strain_placements_added"]:
            md += ["**Strain placements added:** " + "; ".join(f"{x['sgb']} → {str(x['nearest_reference']).split('.fna')[0]} (d={x['distance']})" for x in s["strain_placements_added"]), ""]
        if s["rejected_calls"]:
            md += ["**Rejected by competitive confirmation:**", ""] + [f"- {x['organism']}: {x['reason']}" for x in s["rejected_calls"]] + [""]
        u = s["unresolved_candidates"]
        if u["provisional"] or u["pending_confirmation_genome"] or u["complexes"]:
            md += [f"**Unresolved:** {len(u['provisional'])} provisional, {len(u['pending_confirmation_genome'])} pending a reference genome, "
                   f"{len(u['complexes'])} complexes", ""]
    (args.results / "EXPANSION_COMPARISON.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:len(rows) + 4]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
