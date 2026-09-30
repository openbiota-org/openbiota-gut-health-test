#!/usr/bin/env python3
"""Write docs/MEASURED_CAPABILITY.md from the benchmark and the per-sample comparison.

Spec 0.8.4 §7: the completed report states measured capability. This
document is generated, never hand-edited: it reads the newest benchmark
summary (results/benchmarks/*/summary.json), its organism-free background
controls, and results/EXPANSION_COMPARISON.json, and states what was
measured beside what the spec asks for, gate by gate.

    scripts/measured_capability_doc.py
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.expansion import lanes as _lanes  # noqa: E402

OUT = REPO / "docs" / "MEASURED_CAPABILITY.md"


def _pct(x) -> str:
    return "—" if x is None else f"{100 * float(x):.1f}%"


def main() -> int:
    cap = _lanes.measured_capability(REPO / "results" / "benchmarks")
    comp_path = REPO / "results" / "EXPANSION_COMPARISON.json"
    comp = json.loads(comp_path.read_text()) if comp_path.is_file() else {"samples": []}
    lines = [
        "# Measured capability (generated)",
        "",
        f"Generated {dt.datetime.now().astimezone().isoformat(timespec='minutes')} by "
        "`scripts/measured_capability_doc.py` from the newest benchmark summary and the per-sample comparison. "
        "Numbers here are measurements; the spec's targets are listed beside them and are not claimed where they were not met.",
        "",
    ]
    if cap:
        gates = cap.get("gates") or {}
        ci = cap.get("pooled_precision_ci95") or [None, None]
        gain = cap.get("paired_recall_gain_ci95_bootstrap") or [None, None]
        neg = cap.get("negatives") or {}
        lines += [
            f"## Benchmark `{cap.get('name')}` ({str(cap.get('run_at', ''))[:16]})",
            "",
            f"{cap.get('n_communities')} simulated communit{'y' if cap.get('n_communities') == 1 else 'ies'} × "
            f"{cap.get('species_per_community')} species × "
            f"{int(cap.get('pairs') or 0):,} read pairs (InSilicoSeq, HiSeq error model; truth = HRGM2 representatives frozen to GTDB R232).",
            "",
            "| measure | installed baseline (MP3 + Jun23 + GTDB sylph) | all lanes | spec target |",
            "|---|---|---|---|",
            f"| mean species recall | {_pct(cap.get('mean_recall_baseline'))} | **{_pct(cap.get('mean_recall_all_lanes'))}** | ≥ 95% at 8M pairs |",
            f"| pooled precision (95% CI) | — | **{_pct(cap.get('pooled_precision_all_lanes'))}** ({_pct(ci[0])}–{_pct(ci[1])}) | ≥ 99% |",
            f"| false supported species / community | — | {cap.get('mean_false_supported_per_community')} "
            f"({cap.get('mean_fp_sibling_split_per_community', '—')} sibling-split, {cap.get('mean_fp_other_per_community', '—')} same-genus other) | ≤ 5 |",
            f"| precision at species-complex level (sibling splits credited) | — | {_pct(cap.get('pooled_precision_complex_level'))} | — |",
            f"| paired recall gain (bootstrap 95% CI) | — | {_pct(cap.get('paired_recall_gain_mean'))} ({_pct(gain[0])} to {_pct(gain[1])}) | CI above zero |",
            "",
            "Gates: " + ", ".join(f"{k} = {'met' if v else ('not met' if v is False else 'n/a')}" for k, v in gates.items()),
            "",
        ]
        if neg:
            lines += [
                f"Organism-free backgrounds (human chr21/22 + cattle + chicken DNA, {neg.get('n')} communities): "
                f"**{neg.get('false_supported', 0)} false supported microbial species**, {neg.get('provisional', 0)} provisional "
                "(read classification alone, never supported by rule). Real blanks: none available; disclosed.",
                "",
            ]
        lines += [
            f"Communities run: {cap.get('n_communities')} of the 20 the spec asks for; seeds and depth tiers not yet exhausted. "
            "The confidence intervals are correspondingly wide and the gates above are reported as measured, not as achieved at scale.",
            "",
            "What the false supported species are: on every community measured so far, all of them sit in a genus that is "
            "present in the truth - a relative of a true species called beside it, most often by read classification and "
            "the marker lanes agreeing on a shared-sequence neighbour - and none is an unrelated organism. Competitive "
            "confirmation rejects the shadows it can test; widening it to every same-genus co-call, and calibrating the "
            "Kraken operating point per genus, is the next precision step and is not done.",
            "",
        ]
    else:
        lines += ["No completed benchmark summary found.", ""]
    if comp.get("samples"):
        # Published identifiers only: the documentation is tracked by git and
        # read by anyone; the local sample names stay in results/, which is not.
        from openbiota.samples import published_name

        lines += [
            f"## Per-sample change on the project's {len(comp['samples'])} specimens (unknown truth)",
            "",
            "These are real stool samples with no known truth; the table shows what the expansion changed for "
            "each, not sensitivity or accuracy. Samples are identified by neutral codes.",
            "",
            "| sample | organisms | supported | seen by installed baseline | newly supported | strains placed | rejected | provisional |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for s in sorted(comp["samples"], key=lambda s: published_name(s["sample"])):
            u = s["unresolved_candidates"]
            lines.append(f"| {published_name(s['sample'])} | {s['organisms_total']} | {s['organisms_supported']} | "
                         f"{s['organisms_seen_by_installed_baseline']} | {len(s['newly_supported'])} | "
                         f"{len(s['strain_placements_added'])} | {len(s['rejected_calls'])} | {len(u['provisional'])} |")
        lines += ["", "The per-organism lists behind this table are written to `results/EXPANSION_COMPARISON.md` by "
                  "`scripts/sample_comparison.py` and stay with the results; they name organisms sample by sample and "
                  "are not part of the documentation.", ""]
    OUT.write_text("\n".join(lines))
    print("\n".join(lines[:30]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
