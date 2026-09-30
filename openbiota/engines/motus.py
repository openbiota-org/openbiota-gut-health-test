"""mOTUs 4 lane: marker-gene profiling against ~124k species-level units.

Spec 0.8.4 R4. A different detection method from MetaPhlAn's clade-
specific markers: mOTUs aligns reads (bwa) to ten universal single-copy
marker gene families and resolves them into marker-gene clusters and
mOTUs. Its references overlap GlobDB/GTDB heavily - the point is the
method, not new sequence.

Pinned: tool 4.1.0 (`.venv-motus`), database 4.1 (Zenodo 20322482,
md5-verified by the fetcher). The mOTU count is read from the installed
database rather than from either published figure.

Units: `INSERT_SCALED` counts, and relative abundance computed here over
assigned inserts (the `unassigned` mOTU excluded from the denominator but
recorded). Neither is a MetaPhlAn percentage; the unit travels with the
number.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Final

from openbiota.engines.lanebase import LaneResult, Observation, cache_key, input_sha256, read_cached

LANE_ID: Final = "motus4"
TOOL_VERSION: Final = "4.1.0"
DB_RELEASE: Final = "mOTUs DB 4.1"
#: The GTDB release mOTUs 4.1 names its taxonomy in (from the artifact's GTDB_VERSION column).
TAXONOMY_RELEASE: Final = "GTDB R226"
DB_DIRNAME: Final = "db"
MIN_MARKER_GENES: Final = 3          # mOTUs default: precision/recall balance
MIN_ALIGNMENT_LENGTH: Final = 75
COUNT_MODE: Final = "INSERT_SCALED"
CACHE_VERSION: Final = 1


def _tool() -> Path | None:
    root = Path(__file__).resolve().parents[2]
    exe = root / ".venv-motus" / "bin" / "motus"
    if exe.is_file():
        return exe
    found = shutil.which("motus")
    return Path(found) if found else None


def database_ready(refs_dir: Path) -> bool:
    inner = refs_dir / "motus" / DB_DIRNAME / "db_mOTU"
    return inner.is_dir() and (inner / "db_mOTU.downloaded").is_file() and any(inner.glob("mOTUsv*.db"))


def _db_dir(refs_dir: Path) -> Path:
    """The directory mOTUs wants: the one *containing* `db_mOTU/`."""
    db = refs_dir / "motus" / DB_DIRNAME
    return db if (db / "db_mOTU").is_dir() else db


def available(refs_dir: Path) -> bool:
    return _tool() is not None and database_ready(refs_dir) and shutil.which("bwa") is not None


def _inner(refs_dir: Path) -> Path:
    return refs_dir / "motus" / DB_DIRNAME / "db_mOTU"


def database_motu_count(refs_dir: Path) -> int | None:
    """Number of mOTUs in the installed database, from its own taxonomy file.

    Recorded from the artifact, as the spec asks, rather than from either
    published figure (124,295 on the website, 124,300 in the tutorial).
    """
    import gzip

    db = _inner(refs_dir)
    for cand in sorted(db.glob("*gtdb.taxonomy.rep.tsv*")) + sorted(db.glob("*taxonomy*")):
        opener = gzip.open if cand.suffix == ".gz" else open
        try:
            with opener(cand, "rt") as fh:  # type: ignore[arg-type]
                return sum(1 for line in fh if line.strip() and not line.startswith(("#", "mOTU\t", "MOTU\t")))
        except (OSError, UnicodeDecodeError):
            continue
    return None


def run_motus(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    work_dir: Path,
    threads: int = 8,
) -> LaneResult | None:
    exe = _tool()
    if exe is None or not database_ready(refs_dir):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    sha = input_sha256([r1, r2])
    key = cache_key(str(CACHE_VERSION), TOOL_VERSION, DB_RELEASE, str(MIN_MARKER_GENES),
                    str(MIN_ALIGNMENT_LENGTH), COUNT_MODE, sha)
    cached = work_dir / f"motus.{key}.json"
    hit = read_cached(cached)
    if hit is not None:
        return hit

    out_tsv = work_dir / f"{sample}.motus.tsv"
    cmd = [str(exe), "profile", "-f", str(r1)]
    if r2 is not None:
        cmd += ["-r", str(r2)]
    cmd += ["-n", sample, "-db", str(_db_dir(refs_dir)), "-t", str(threads), "-o", str(out_tsv),
            "-g", str(MIN_MARKER_GENES), "-l", str(MIN_ALIGNMENT_LENGTH), "-y", COUNT_MODE]
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not out_tsv.is_file():
        raise RuntimeError(f"mOTUs failed (exit {proc.returncode}): {proc.stderr[-2000:]}")
    observations, summary = _parse(out_tsv)
    result = LaneResult(
        lane_id=LANE_ID, tool="motus", tool_version=TOOL_VERSION, reference_release_id=DB_RELEASE,
        taxonomy_release=f"mOTUs 4.1 taxonomy, {TAXONOMY_RELEASE} names (walked to R232 by the crosswalk)", observations=observations,
        elapsed_s=time.monotonic() - t0, cached=False, command=tuple(cmd), input_sha256=sha,
        raw_result_uri=str(out_tsv), summary=summary,
        notes=[f"mOTU called present with >= {MIN_MARKER_GENES} marker genes; abundances are {COUNT_MODE} "
               "inserts and a relative abundance over assigned inserts"],
    )
    cached.write_text(json.dumps(result.to_json(), indent=1))
    return result


def _parse(path: Path) -> tuple[list[Observation], dict[str, Any]]:
    rows: list[tuple[str, str, float]] = []
    unassigned = 0.0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("mOTU\t") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            motu, tax, val = parts[0], parts[1], parts[2]
            try:
                v = float(val)
            except ValueError:
                continue
            if v <= 0:
                continue
            if "unassigned" in motu.lower() or "unassigned" in tax.lower():
                unassigned += v
                continue
            rows.append((motu, tax, v))
    total = sum(v for _, _, v in rows)
    obs: list[Observation] = []
    for motu, tax, v in rows:
        lineage = tax.replace("|", ";")
        rank = "species" if "s__" in lineage and lineage.rsplit("s__", 1)[-1].strip() not in ("", "unassigned") else "genus"
        obs.append(Observation(
            native_id=motu, lineage=lineage, rank=rank, status="supported",
            abundance_value=(100.0 * v / total) if total else None,
            abundance_unit="percent of mOTU-assigned inserts", denominator="mOTU-assigned inserts",
            support_metrics={"insert_scaled": round(v, 4), "counting_mode": COUNT_MODE},
        ))
    obs.sort(key=lambda o: -(o.abundance_value or 0.0))
    summary = {
        "assigned_inserts_scaled": round(total, 2), "unassigned_inserts_scaled": round(unassigned, 2),
        "unassigned_fraction": round(unassigned / (total + unassigned), 4) if (total + unassigned) else None,
        "n_motus": len(obs),
    }
    return obs, summary
