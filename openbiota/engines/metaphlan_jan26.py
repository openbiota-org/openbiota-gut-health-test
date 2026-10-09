"""MetaPhlAn 4.2.6 / Jan26 marker lane (spec 0.8.4 R1).

The upgraded detection lane: 72,000 SGBs (68,266 bacterial, 3,245 archaeal,
489 eukaryotic) against the installed Jun23 index's 36,822. The Jun23
adapter (`openbiota.engines.metaphlan4`, MetaPhlAn 4.1.1 in `.venv-mpa4`)
is left exactly as it is: it is the reproducible baseline the spec asks to
preserve, and its raw results stay on disk.

This lane runs from `.venv-mpa42` (MetaPhlAn 4.2.6, tag c55b299; the tag
self-reports 4.2.5) with `--offline`, so a sample run can never resolve a
newer database. It keeps the mapout and the SAM: StrainPhlAn's
`sample2markers` needs the Jan26 SAM, and Jun23 consensus objects are not
compatible inputs for it.

Units: MetaPhlAn relative abundance, percent of classified reads
(unclassified estimation on, so the UNCLASSIFIED row is recorded too).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Final

from openbiota.engines.lanebase import LaneResult, Observation, cache_key, input_sha256, read_cached
from openbiota.engines.metaphlan import _env_for

LANE_ID: Final = "metaphlan_jan26"
INDEX: Final = "mpa_vJan26_CHOCOPhlAnSGB_202605"
TOOL_VERSION: Final = "4.2.6 (tag c55b299; self-reports 4.2.5)"
SGB_TOTAL: Final = 72_000
CACHE_VERSION: Final = 1


def _tool() -> Path | None:
    root = Path(__file__).resolve().parents[2]
    explicit = os.environ.get("OPENBIOTA_METAPHLAN42")
    if explicit and Path(explicit).is_file():
        return Path(explicit)
    exe = root / ".venv-mpa42" / "bin" / "metaphlan"
    return exe if exe.is_file() else None


def database_ready(db_dir: Path) -> bool:
    """The index is installed when every bowtie2 shard and the pickle exist."""
    shards = list(db_dir.glob(f"{INDEX}.*.bt2l")) + list(db_dir.glob(f"{INDEX}.*.bt2"))
    return len(shards) >= 6 and (db_dir / f"{INDEX}.pkl").is_file()


def available(db_dir: Path) -> bool:
    return _tool() is not None and database_ready(db_dir) and shutil.which("bowtie2") is not None


def cache_file(work_dir: Path, r1: Path, r2: Path | None) -> Path:
    """Where this lane's result for these reads is cached. One place, so a tool
    that renames inputs can re-key the cache instead of recomputing the lane."""
    sha = input_sha256([r1, r2])
    return work_dir / f"metaphlan_jan26.{cache_key(str(CACHE_VERSION), TOOL_VERSION, INDEX, sha)}.json"


def run_metaphlan_jan26(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    db_dir: Path,
    work_dir: Path,
    strain_dir: Path,
    threads: int = 8,
) -> LaneResult | None:
    exe = _tool()
    if exe is None or not database_ready(db_dir):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    strain_dir.mkdir(parents=True, exist_ok=True)
    sha = input_sha256([r1, r2])
    cached = cache_file(work_dir, r1, r2)
    hit = read_cached(cached)
    if hit is not None and Path(hit.raw_result_uri).is_file():
        return hit

    profile = work_dir / f"metaphlan_jan26.{INDEX}.tsv"
    mapout = work_dir / f"metaphlan_jan26.{INDEX}.mapout.bz2"
    sam = strain_dir / f"{sample}.jan26.markers.sam.bz2"
    inputs = str(r1) if r2 is None else f"{r1},{r2}"
    cmd = [str(exe), inputs, "--input_type", "fastq", "--offline", "--db_dir", str(db_dir), "-x", INDEX,
           "--nproc", str(threads), "--mapout", str(mapout), "-s", str(sam), "-o", str(profile)]
    # An interrupted run leaves empty or truncated intermediates behind; a
    # mapout or SAM under 1 MB cannot be a real alignment of a sample and
    # is discarded so the lane aligns afresh rather than failing on it.
    for stale in (mapout, sam, profile):
        if stale.is_file() and stale.stat().st_size < (1_000_000 if stale is not profile else 200):
            stale.unlink()
    if mapout.is_file() and not sam.is_file():
        mapout.unlink()  # MetaPhlAn refuses to overwrite a mapout; the SAM only comes from a fresh alignment
    if mapout.is_file() and sam.is_file() and not profile.is_file():
        cmd = [str(exe), str(mapout), "--input_type", "mapout", "--offline", "--db_dir", str(db_dir), "-x", INDEX,
               "--nproc", str(threads), "-o", str(profile)]
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=_env_for(str(exe)))
    if proc.returncode != 0 or not profile.is_file():
        raise RuntimeError(f"MetaPhlAn Jan26 failed (exit {proc.returncode}): {proc.stderr[-3000:]}")
    observations, summary = _parse(profile)
    summary["sam"] = str(sam) if sam.is_file() else None
    summary["mapout"] = str(mapout) if mapout.is_file() else None
    result = LaneResult(
        lane_id=LANE_ID, tool="metaphlan", tool_version=TOOL_VERSION, reference_release_id=INDEX,
        taxonomy_release="SGB Jan26; GTDB R232 via refs/crosswalk (assembly-membership walk from the R226 bridge)",
        observations=observations, elapsed_s=time.monotonic() - t0, cached=False, command=tuple(cmd),
        input_sha256=sha, raw_result_uri=str(profile), summary=summary,
        notes=[f"{SGB_TOTAL:,} SGBs in the index; relative abundance is percent of classified reads",
               "SAM and mapout retained for StrainPhlAn"],
    )
    cached.write_text(json.dumps(result.to_json(), indent=1))
    return result


def _parse(profile: Path) -> tuple[list[Observation], dict[str, Any]]:
    obs: list[Observation] = []
    unclassified = 0.0
    n_reads = 0
    with profile.open() as fh:
        for line in fh:
            if line.startswith("#"):
                if "reads processed" in line:
                    digits = "".join(ch for ch in line if ch.isdigit())
                    n_reads = int(digits) if digits else 0
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            clade, taxid, rel = parts[0], parts[1], parts[2]
            try:
                value = float(rel)
            except ValueError:
                continue
            if clade == "UNCLASSIFIED":
                unclassified = value
                continue
            if "|t__" not in clade:
                continue  # SGB rows carry the species and the native id
            segs = clade.split("|")
            sgb = segs[-1][3:]
            lineage = ";".join(s for s in segs[:-1])
            species = next((s[3:] for s in segs if s.startswith("s__")), "")
            obs.append(Observation(
                native_id=sgb, lineage=lineage, rank="species" if species else "genus", status="supported",
                abundance_value=value, abundance_unit="percent of classified reads (MetaPhlAn relative abundance)",
                denominator="classified reads", support_metrics={"ncbi_taxid": taxid, "unnamed_sgb": species.endswith("_SGB" + sgb[3:]) or "_SGB" in species},
            ))
    obs.sort(key=lambda o: -(o.abundance_value or 0.0))
    phyla: dict[str, float] = {}
    kingdoms: dict[str, float] = {}
    for o in obs:
        segs = o.lineage.split(";")
        ph = next((x[3:] for x in segs if x.startswith("p__")), "")
        kd = next((x[3:] for x in segs if x.startswith("k__")), "")
        if ph:
            phyla[ph] = phyla.get(ph, 0.0) + float(o.abundance_value or 0.0)
        if kd:
            kingdoms[kd] = kingdoms.get(kd, 0.0) + float(o.abundance_value or 0.0)
    return obs, {"unclassified_percent": unclassified, "reads_processed": n_reads, "n_sgbs": len(obs),
                 "n_unnamed_sgbs": sum(1 for o in obs if o.support_metrics.get("unnamed_sgb")),
                 "phyla": dict(sorted(phyla.items(), key=lambda kv: -kv[1])),
                 "kingdoms": dict(sorted(kingdoms.items(), key=lambda kv: -kv[1]))}
