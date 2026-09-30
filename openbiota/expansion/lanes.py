"""Run every expanded detection lane over one sample and record the outcome.

One entry point for the pipeline. Each lane runs on the same host-filtered
reads, in its own try/except: a lane that is not installed or that fails
is recorded as `not_assessed` with the reason, and never as a set of
absent organisms. The block written to results.json carries every lane's
native observations (spec §6 `DetectionObservation`), the reference
releases each lane used (from the locks), and the lane status table the
completeness gate reads.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from openbiota.engines import kraken, metaphlan_jan26, motus, singlem, sylph_globdb
from openbiota.engines.lanebase import LaneResult

LANE_ORDER = ("metaphlan_jan26", "sylph_globdb", "motus4", "kraken_uhgg", "kraken_rescue", "singlem_globdb")


def _reference_releases() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    lock_dir = Path("refs/expanded/locks")
    if not lock_dir.is_dir():
        return out
    for path in sorted(lock_dir.glob("*.lock.json")):
        try:
            lock = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        files = lock.get("files") or {}
        out.append({
            "source_id": lock.get("source_id"), "release_id": lock.get("release_id"),
            "taxonomy_release": lock.get("taxonomy_release"), "license": lock.get("license"),
            "count_unit": lock.get("count_unit"), "source_count": lock.get("source_count"),
            "n_files": len(files), "verified": all(f.get("verified") for f in files.values()) if files else False,
            "sha256": {k: v.get("sha256") for k, v in files.items()},
            "urls": [v.get("url") for v in files.values()],
            "retrieved_at": max((v.get("retrieved_at") or "" for v in files.values()), default=""),
        })
    return out


def measured_capability(root: Path = Path("results/benchmarks")) -> dict[str, Any]:
    """The newest completed benchmark summary, reduced to what the report states.

    Spec §7: the completed report states measured capability. Only summaries
    with at least one truth community count; the negative-control block is
    carried when present. Nothing here is per sample: it is what the
    installed lanes measured on simulated communities with known truth.
    """
    newest: tuple[str, Path] | None = None
    if not root.is_dir():
        return {}
    for p in root.glob("*/summary.json"):
        try:
            data = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not data.get("n_communities"):
            continue
        stamp = str(data.get("run_at") or "")
        if newest is None or stamp > newest[0]:
            newest = (stamp, p)
    if newest is None:
        return {}
    data = json.loads(newest[1].read_text())
    negatives = data.get("negatives") or []
    keys = ("name", "run_at", "n_communities", "species_per_community", "pairs", "mean_recall_baseline",
            "mean_recall_all_lanes", "pooled_precision_all_lanes", "pooled_precision_ci95",
            "mean_false_supported_per_community", "paired_recall_gain_mean", "paired_recall_gain_ci95_bootstrap", "gates",
            "pooled_precision_complex_level", "mean_fp_sibling_split_per_community", "mean_fp_other_per_community")
    out = {k: data.get(k) for k in keys if k in data}
    if negatives:
        out["negatives"] = {"n": len(negatives), "false_supported": sum(int(n.get("n_supported") or 0) for n in negatives),
                            "provisional": sum(int(n.get("n_provisional") or 0) for n in negatives)}
    out["source"] = str(newest[1])
    return out


def run_all(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    out_dir: Path,
    threads: int,
    reporter: Any,
    only: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Every lane, each recorded whether it ran or not."""
    lanes: dict[str, dict[str, Any]] = {}
    status: dict[str, dict[str, Any]] = {}

    def attempt(lane_id: str, label: str, fn: Callable[[], LaneResult | None]) -> None:
        if only is not None and lane_id not in only:
            status[lane_id] = {"state": "skipped", "reason": "not requested"}
            return
        t0 = time.monotonic()
        try:
            with reporter.stage(label):
                result = fn()
                if result is None:
                    reason = ("panel not built (scripts/build_rescue_panel_kraken.py)" if lane_id == "kraken_rescue"
                              else "tool or reference not installed")
                    status[lane_id] = {"state": "not_assessed", "reason": reason}
                    reporter.warn(f"{label}: {reason}; recorded as not assessed")
                    return
                lanes[lane_id] = result.to_json()
                status[lane_id] = {
                    "state": "ran", "cached": result.cached, "elapsed_s": round(result.elapsed_s, 1),
                    "n_observations": len(result.observations),
                    "n_species": sum(1 for o in result.observations if o.rank == "species"),
                }
                reporter.record(
                    f"    {status[lane_id]['n_species']} species-level observations"
                    + (" (cached)" if result.cached else f" in {result.elapsed_s:.0f}s")
                )
        except Exception as exc:  # noqa: BLE001 - a failed lane is a recorded absence of assessment
            status[lane_id] = {"state": "failed", "reason": f"{type(exc).__name__}: {str(exc)[:300]}",
                               "error_detail": str(exc)[-2500:], "elapsed_s": round(time.monotonic() - t0, 1)}
            reporter.warn(f"{label} failed: {type(exc).__name__}: {str(exc)[:200]}")

    attempt("metaphlan_jan26", "marker lane (MetaPhlAn 4.2.6 / Jan26)", lambda: metaphlan_jan26.run_metaphlan_jan26(
        sample=sample, r1=r1, r2=r2, db_dir=refs_dir / "metaphlan4_db", work_dir=out_dir / "taxonomy",
        strain_dir=out_dir / "strain", threads=threads))
    attempt("sylph_globdb", "broad discovery (sylph / GlobDB r232)", lambda: sylph_globdb.run_sylph_globdb(
        sample=sample, r1=r1, r2=r2, refs_dir=refs_dir, work_dir=out_dir / "globdb", threads=threads))
    attempt("motus4", "universal-marker lane (mOTUs 4.1)", lambda: motus.run_motus(
        sample=sample, r1=r1, r2=r2, refs_dir=refs_dir, work_dir=out_dir / "motus", threads=threads))
    attempt("kraken_uhgg", "gut rescue classifier (Kraken2 + Bracken / UHGG v2.0.2)", lambda: kraken.run_kraken(
        sample=sample, r1=r1, r2=r2, refs_dir=refs_dir, work_dir=out_dir / "kraken", threads=threads))
    attempt("kraken_rescue", "gut rescue panel (Kraken2 + Bracken / reconciliation panel)", lambda: kraken.run_kraken(
        sample=sample, r1=r1, r2=r2, refs_dir=refs_dir, work_dir=out_dir / "kraken_rescue", threads=threads, panel="rescue"))
    attempt("singlem_globdb", "unrepresented-lineage rescue (SingleM / GlobDB)", lambda: singlem.run_singlem(
        sample=sample, r1=r1, r2=r2, refs_dir=refs_dir, work_dir=out_dir / "singlem", threads=threads))

    return {
        "schema_version": "openbiota.detection/1.0",
        "lanes": lanes,
        "lane_status": status,
        "reference_releases": _reference_releases(),
        "what_this_is": (
            "Every expanded detection lane's own observations on this sample, in its own units and identifiers. "
            "Nothing here is summed across lanes; the organism inventory merges identities and takes composition "
            "from one primary lane. A lane recorded as not_assessed or failed says nothing about which organisms "
            "it would have found."
        ),
    }
