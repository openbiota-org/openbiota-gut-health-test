#!/usr/bin/env python
"""Build reference cohorts in the expanded detection lanes' own namespaces.

A percentile is a position inside a population that was measured the same
way. The scoring lane has one (3,027 adults, curatedMetagenomicData,
MetaPhlAn 3). The expanded lanes - GlobDB via sylph, MetaPhlAn Jan26,
mOTUs 4, Kraken/UHGG - had none, so an organism they alone detect could be
reported as present but never as high or low for a person.

This script gives each lane its own population: it fetches public adult
stool metagenomes, runs the lane's engine on each exactly as the pipeline
does, and reduces the profiles to a `ReferenceCohort` matrix
(`refs/expansion_cohort/<lane>.cohort.json`) that `openbiota.inventory`
reads to rank that lane's organisms. Lanes are never mixed: one cohort per
lane, and an organism's percentile always names the population it came
from.

The default population is the same 100 adults of ENA study PRJNA1271016
that the pathway reference ranges rest on, so both sets of ranges describe
one group of people. Depth is a prefix of each run (the head of the file);
the manifest records the depth, the sampling, the engine and the database
release, because a percentile is only as honest as its provenance.

    scripts/build_expansion_cohort.py --lanes globdb --pairs 3000000
    scripts/build_expansion_cohort.py --status
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import shutil
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from openbiota import cohort as cohort_mod  # noqa: E402
from openbiota import inventory  # noqa: E402
from openbiota.engines import kraken, metaphlan_jan26, motus, sylph_globdb  # noqa: E402
from openbiota.logging_util import Reporter  # noqa: E402

OUT_DIR = ROOT / "refs" / "expansion_cohort"
RUNS_FILE = ROOT / "refs" / "reference_cohort_samples.json"
STUDY = cohort_mod.DEFAULT_STUDY

#: lane -> (engine callable, inventory lane builder, human label)
LANES: dict[str, dict[str, Any]] = {
    "globdb": {
        "engine": lambda **kw: sylph_globdb.run_sylph_globdb(**kw),
        "work": "globdb", "builder": inventory._globdb_lane,
        "label": "sylph 1.0.0 against GlobDB r232",
    },
    "jan26": {
        "engine": lambda **kw: metaphlan_jan26.run_metaphlan_jan26(
            sample=kw["sample"], r1=kw["r1"], r2=kw["r2"], db_dir=kw["refs_dir"] / "metaphlan4_db",
            work_dir=kw["work_dir"], strain_dir=kw["work_dir"] / "strain", threads=kw["threads"]),
        "work": "taxonomy", "builder": inventory._jan26_lane,
        "label": "MetaPhlAn 4.2.6 against mpa_vJan26_CHOCOPhlAnSGB_202605",
    },
    "motus": {
        "engine": lambda **kw: motus.run_motus(**kw),
        "work": "motus", "builder": inventory._motus_lane,
        "label": "mOTUs 4.1.0 against mOTUs DB 4.1",
    },
    "kraken": {
        "engine": lambda **kw: kraken.run_kraken(**kw),
        "work": "kraken", "builder": inventory._kraken_lane,
        "label": "Kraken2 + Bracken against UHGG v2.0.2",
    },
}


def _runs(reporter: Reporter, limit: int) -> list[cohort_mod.CohortRun]:
    wanted = [r["run"] for r in json.loads(RUNS_FILE.read_text())][:limit]
    listed = {r.run: r for r in cohort_mod.list_runs(STUDY, reporter)}
    missing = [r for r in wanted if r not in listed]
    if missing:
        reporter.warn(f"{len(missing)} runs from {RUNS_FILE.name} are not in ENA's listing of {STUDY}: {missing[:3]}")
    return [listed[r] for r in wanted if r in listed]


def _profile_path(lane: str, run: str) -> Path:
    return OUT_DIR / "profiles" / lane / f"{run}.json"


def _status_line(lanes: list[str], runs: list[str]) -> str:
    parts = []
    for lane in lanes:
        n = sum(1 for r in runs if _profile_path(lane, r).is_file())
        parts.append(f"{lane} {n}/{len(runs)}")
    return "; ".join(parts)


def profile_run(run: cohort_mod.CohortRun, r1: Path, r2: Path, lanes: list[str], threads: int,
                reporter: Reporter) -> None:
    for lane in lanes:
        dest = _profile_path(lane, run.run)
        if dest.is_file():
            continue
        spec = LANES[lane]
        work_dir = OUT_DIR / "work" / run.run / spec["work"]
        work_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.monotonic()
        try:
            result = spec["engine"](sample=run.run, r1=r1, r2=r2, refs_dir=ROOT / "refs", work_dir=work_dir, threads=threads)
        except Exception as exc:  # noqa: BLE001 - one failed run is one missing column, recorded
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.with_suffix(".failed").write_text(f"{type(exc).__name__}: {exc}"[:2000])
            reporter.warn(f"  {run.run} {lane}: {type(exc).__name__}: {str(exc)[:160]}")
            continue
        if result is None:
            reporter.warn(f"  {run.run} {lane}: engine or database not installed")
            continue
        blob = result.to_json()
        lane_obj = spec["builder"](blob)
        payload = {
            "run": run.run, "sample": run.sample, "lane": lane, "tool": blob.get("tool"),
            "tool_version": blob.get("tool_version"), "reference_release_id": blob.get("reference_release_id"),
            "taxonomy_release": blob.get("taxonomy_release"), "elapsed_s": round(time.monotonic() - t0, 1),
            "summary": blob.get("summary") or {},
            "species": dict(lane_obj.species) if lane_obj is not None else {},
        }
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(payload, separators=(",", ":")))
        reporter.record(f"  {run.run} {lane}: {len(payload['species'])} species in {payload['elapsed_s']:.0f}s")
    shutil.rmtree(OUT_DIR / "work" / run.run, ignore_errors=True)


def assemble(lane: str, runs: list[str], pairs: int, sampling: str, reporter: Reporter) -> Path | None:
    """One `ReferenceCohort`-shaped JSON from the lane's per-run profiles."""
    profiles = []
    for run in runs:
        p = _profile_path(lane, run)
        if p.is_file():
            profiles.append(json.loads(p.read_text()))
    if len(profiles) < 30:
        reporter.warn(f"{lane}: only {len(profiles)} profiles; a cohort needs at least 30 to say anything about tails")
        return None
    taxa = sorted({s for prof in profiles for s in prof["species"]})
    index = {t: i for i, t in enumerate(taxa)}
    abundance = [[0.0] * len(profiles) for _ in taxa]
    for j, prof in enumerate(profiles):
        for s, v in prof["species"].items():
            abundance[index[s]][j] = round(float(v), 6)
    sample_ids = [prof["run"] for prof in profiles]
    versions = sorted({(prof.get("tool_version") or "", prof.get("reference_release_id") or "") for prof in profiles})
    if len(versions) != 1:
        raise SystemExit(f"{lane}: profiles come from more than one engine/database version: {versions}; "
                         "a cohort may not mix them. Delete the stale profiles and rerun.")
    tool_version, release = versions[0]
    metadata = {
        prof["run"]: {
            "sample_id": prof["run"], "study_name": STUDY, "condition": "adult stool", "country": None,
            "age": None, "sex": None, "bmi": None, "n_reads": pairs, "platform": "Illumina NovaSeq X",
            "antibiotics": None, "westernized": None,
        }
        for prof in profiles
    }
    manifest = {
        "manifest_id": f"{lane}-{len(profiles)}-{pairs}",
        "snapshot": release,
        "profiler": f"{profiles[0].get('tool')} {tool_version}",
        "built_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "n_samples": len(profiles), "n_taxa": len(taxa),
        "source": f"ENA study {STUDY} (adult stool, NovaSeq X), profiled locally by the same engine the pipeline uses",
        "study": STUDY,
        "lane": lane, "lane_label": LANES[lane]["label"],
        "taxonomy_release": profiles[0].get("taxonomy_release"),
        "depth_read_pairs": pairs, "sampling": sampling,
        "matched": False,
        "caveat": (f"{len(profiles)} unmatched adults at {pairs:,} read pairs each: the 5th and 95th percentiles "
                   f"rest on about {max(1, len(profiles) // 20)} people, and organisms below roughly "
                   f"{100.0 / (pairs * 2 * 150 / 3_000_000):.3f}% of a 3 Mb genome's coverage are near the "
                   "detection floor of this depth, so their prevalence here is a floor, not a fact"),
        "citation": "profiles computed by OpenBiota; reads from ENA study " + STUDY,
    }
    out = OUT_DIR / f"{lane}.cohort.json"
    out.write_text(json.dumps({"manifest": manifest, "taxa": taxa, "sample_ids": sample_ids,
                               "abundance": abundance, "metadata": metadata}, separators=(",", ":")))
    reporter.info(f"{lane}: cohort of {len(profiles)} samples x {len(taxa)} taxa -> {out}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lanes", default="globdb", help="comma-separated: " + ",".join(LANES))
    ap.add_argument("--pairs", type=int, default=3_000_000, help="read pairs per run (prefix)")
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--workers", type=int, default=4, help="parallel downloads")
    ap.add_argument("--parallel", type=int, default=2, help="runs profiled at once")
    ap.add_argument("--threads", type=int, default=6, help="threads per engine run")
    ap.add_argument("--transport", default="ena", choices=cohort_mod.TRANSPORTS)
    ap.add_argument("--sampling", default="prefix", choices=cohort_mod.SAMPLING_MODES)
    ap.add_argument("--status", action="store_true", help="print progress and exit")
    ap.add_argument("--assemble-only", action="store_true", help="skip fetching/profiling; rebuild the cohort files")
    ap.add_argument("--keep-fastq", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    lanes = [x.strip() for x in args.lanes.split(",") if x.strip()]
    unknown = [x for x in lanes if x not in LANES]
    if unknown:
        raise SystemExit(f"unknown lanes {unknown}; choose from {list(LANES)}")
    reporter = Reporter(verbose=not args.quiet)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wanted = [r["run"] for r in json.loads(RUNS_FILE.read_text())][:args.samples]
    if args.status:
        print(_status_line(lanes, wanted))
        return 0
    if not args.assemble_only:
        runs = _runs(reporter, args.samples)
        todo = [r for r in runs if any(not _profile_path(ln, r.run).is_file() for ln in lanes)]
        reporter.info(f"{len(runs)} runs; {len(todo)} still to profile for {lanes}")
        fastq_dir = OUT_DIR / "fastq"
        fastq_dir.mkdir(exist_ok=True)
        progress = OUT_DIR / "progress.txt"
        # Fetch in batches so profiling overlaps downloading and the disk
        # never holds more than a batch of prefixes.
        batch = max(args.parallel * 2, 4)
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel)
        pending: list[concurrent.futures.Future[None]] = []
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            fetched = cohort_mod.download_cohort(
                chunk, dest_dir=fastq_dir, read_pairs=args.pairs, workers=args.workers, reporter=reporter,
                sampling=args.sampling, transport=args.transport)

            def job(item: tuple[cohort_mod.CohortRun, Path, Path]) -> None:
                run, r1, r2 = item
                try:
                    profile_run(run, r1, r2, lanes, args.threads, reporter)
                finally:
                    if not args.keep_fastq:
                        for p in (r1, r2, r1.with_suffix(r1.suffix + ".done"), r2.with_suffix(r2.suffix + ".done")):
                            with_suppress_unlink(p)
                    progress.write_text(f"{datetime.now().strftime('%H:%M:%S')} {_status_line(lanes, wanted)}\n")

            pending.extend(pool.submit(job, item) for item in fetched)
            # Keep at most one batch in flight ahead of the profilers.
            while sum(1 for f in pending if not f.done()) > batch:
                time.sleep(5)
        for f in pending:
            f.result()
        pool.shutdown(wait=True)
    for lane in lanes:
        assemble(lane, wanted, args.pairs, args.sampling, reporter)
    return 0


def with_suppress_unlink(p: Path) -> None:
    try:
        p.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        pass


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    raise SystemExit(main())
