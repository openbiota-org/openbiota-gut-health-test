"""`openbiota mycobiome` command surface (spec §13).

    openbiota mycobiome references prepare   lock + sketch DB + species index + strain panels
    openbiota mycobiome references status    what is built, with counts
    openbiota mycobiome analyze --sample S --r1 ... --r2 ... --output results/S
    openbiota mycobiome score --sample-results results/S
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from openbiota.mycobiome import engine, panel_build, references, registry, score


def _say(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def run(args: argparse.Namespace) -> int:
    cmd = getattr(args, "mycobiome_cmd", None)
    if cmd == "references":
        return _references(args)
    if cmd == "analyze":
        return _analyze(args)
    if cmd == "score":
        return _score(args)
    _say("usage: openbiota mycobiome {references,analyze,score} ...")
    return 2


def _references(args: argparse.Namespace) -> int:
    base = Path(args.refs_dir) / "mycobiome"
    if args.verb == "status":
        lock = references.load_lock(base)
        if lock is None:
            _say("no reference lock; run `openbiota mycobiome references prepare`")
            return 1
        _say(f"lock {lock.lock_id} built {lock.built_at}")
        for k, v in lock.counts.items():
            _say(f"  {k}: {v:,}")
        _say(f"  sketch db: {lock.sketch_db}")
        _say(f"  species index: {'yes' if (base / 'species_index' / 'fungi_species.mmi').is_file() else 'no'}")
        panels = sorted((base / "panels").glob("*.panel.json"))
        _say(f"  strain panels: {len(panels)} " + ", ".join(p.name for p in panels))
        _say(f"  deposits: {len(lock.deposits)} files recorded")
        return 0
    _say("building reference lock and sketch database")
    lock = references.build_lock(base, sketch=True, threads=args.threads, log=_say)
    _say(f"lock {lock.lock_id}: {lock.counts}")
    _say("building species-representative index")
    references.build_species_index(lock, base, threads=args.threads, log=_say)
    if not args.no_panels:
        _say("deriving strain panels")
        for p in panel_build.build_all(lock, base, threads=args.threads, log=_say):
            _say(f"  {p.name}")
    return 0


def _analyze(args: argparse.Namespace) -> int:
    out_dir = Path(args.output)
    results_path = out_dir / "results.json"
    if not results_path.is_file():
        _say(f"no results.json under {out_dir}; run the main pipeline first")
        return 1
    results = json.loads(results_path.read_text())
    sketch = out_dir / "genome" / "sketches" / f"{args.sample}.paired.sylsp"
    rec = engine.analyze(sample=args.sample, r1=Path(args.r1), r2=Path(args.r2) if args.r2 else None, results=results,
                         out_dir=out_dir, refs_dir=Path(args.refs_dir), threads=args.threads,
                         read_sketch=sketch if sketch.is_file() else None, log=_say)
    results["mycobiome"] = rec
    results_path.write_text(json.dumps(results, indent=1, default=str))
    hs = rec.get("health_score") or {}
    _say(f"{args.sample}: {rec['analysis_status']}; score {hs.get('score_display')} ({hs.get('status')}, {hs.get('direction')})")
    return 0


def _score(args: argparse.Namespace) -> int:
    """Recompute the score from the stored record (same rule manifest)."""
    results_path = Path(args.sample_results) / "results.json"
    results = json.loads(results_path.read_text())
    rec = results.get("mycobiome") or {}
    if not rec.get("taxa"):
        _say("no mycobiome record to score")
        return 1
    reg = registry.load()
    supported = [
        registry.SupportedTaxon(
            finding_id=t["finding_id"], name=t["accepted_name"], rank=t["rank"], detection_state=t["detection_state"],
            strain_resolution=t["strain"]["resolution"], reference_accessions=tuple(t["strain"]["reference_accessions"]),
            compatible_genotypes=tuple(t["strain"]["compatible_genotypes"]),
            compatible_accessions=tuple(t["strain"]["reference_accessions"]) if t["strain"]["compatible_genotypes"] else (),
            strain_analysis_status=t["strain"]["analysis_status"],
        )
        for t in rec["taxa"]
    ]
    bact = [str(o.get("species", "")).replace("_", " ") for o in (results.get("organism_inventory") or {}).get("organisms") or []]
    contribs = registry.evaluate(reg, supported, bacteria_present=bact)
    sc = score.compute(contribs, policy_lock_id=reg.policy_lock_id)
    rec["health_score"] = sc.to_json()
    results_path.write_text(json.dumps(results, indent=1, default=str))
    _say(f"score {sc.display_integer} ({sc.status}, {sc.direction}); policy {reg.policy_lock_id}")
    return 0


__all__ = ["run"]
