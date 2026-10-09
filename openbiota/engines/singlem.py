"""SingleM lane: unrepresented-lineage rescue (spec 0.8.4 §4D).

SingleM reads single-copy marker windows straight out of the reads and
places each against the GlobDB r232 metapackage, so a lineage that no
genome catalogue represents still shows up - at the deepest rank its
marker windows support. That is the lane's job here: every *unresolved*
lineage (placed at genus or above with real coverage) is a candidate for
targeted assembly, because coverage that high with no species-level home
is either a novel organism or a catalogue gap, and only sequence decides
which.

Coverage is SingleM's own unit (mean marker coverage, roughly genome
coverage). A marker fragment is not an organism; an unresolved lineage
with coverage under `ASSEMBLY_COVERAGE` is reported at its supported rank
and nothing more.

Assembly is triggered when the lane reports an unresolved lineage at or
above `ASSEMBLY_COVERAGE`; `openbiota.expansion.assembly` runs it. The
trigger is a data condition, not a feature flag.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tarfile
import time
from pathlib import Path
from typing import Any, Final

from openbiota.engines.lanebase import LaneResult, Observation, cache_key, input_sha256, read_cached

LANE_ID: Final = "singlem_globdb"
METAPACKAGE_TAR: Final = "GlobDB_r232.metapackage_v4.smpkg.tar.gz"
METAPACKAGE_DIR: Final = "GlobDB_r232.metapackage_v4.smpkg"
DB_RELEASE: Final = "GlobDB_r232.metapackage_v4"
#: Mean marker coverage at which an unresolved lineage is worth assembling:
#: below ~3x, an assembler returns fragments that cannot be binned or judged.
ASSEMBLY_COVERAGE: Final = 3.0
CACHE_VERSION: Final = 1

_RANKS = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


def _tool() -> Path | None:
    root = Path(__file__).resolve().parents[2]
    exe = root / ".venv-singlem" / "bin" / "singlem"
    if exe.is_file():
        return exe
    found = shutil.which("singlem")
    return Path(found) if found else None


def _pkg_dir(refs_dir: Path) -> Path:
    return refs_dir / "singlem" / METAPACKAGE_DIR


def ensure_extracted(refs_dir: Path) -> bool:
    pkg = _pkg_dir(refs_dir)
    tar = refs_dir / "singlem" / METAPACKAGE_TAR
    if pkg.is_dir() and any(pkg.iterdir()):
        return True
    if not tar.is_file():
        return False
    lock = Path("refs/expanded/locks/singlem_globdb.lock.json")
    if lock.is_file():
        rec = json.loads(lock.read_text()).get("files", {}).get(f"singlem/{METAPACKAGE_TAR}")
        if not (rec and rec.get("verified")):
            return False
    with tarfile.open(tar) as tf:
        tf.extractall(refs_dir / "singlem")  # noqa: S202 - publisher archive, md5-verified by the fetcher
    return pkg.is_dir()


def _version(exe: Path) -> str:
    out = subprocess.run([str(exe), "--version"], capture_output=True, text=True, check=False)
    return (out.stdout or out.stderr).strip().splitlines()[-1] if (out.stdout or out.stderr) else ""


def available(refs_dir: Path) -> bool:
    exe = _tool()
    return exe is not None and ensure_extracted(refs_dir)


def cache_file(work_dir: Path, r1: Path, r2: Path | None, *, version: str | None = None) -> Path:
    """Where this lane's result for these reads is cached (see metaphlan_jan26.cache_file)."""
    if version is None:
        exe = _tool()
        version = _version(exe) if exe is not None else ""
    sha = input_sha256([r1, r2])
    return work_dir / f"singlem.{cache_key(str(CACHE_VERSION), version, DB_RELEASE, str(ASSEMBLY_COVERAGE), sha)}.json"


def run_singlem(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    work_dir: Path,
    threads: int = 8,
) -> LaneResult | None:
    exe = _tool()
    if exe is None or not ensure_extracted(refs_dir):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    version = _version(exe)
    sha = input_sha256([r1, r2])
    cached = cache_file(work_dir, r1, r2, version=version)
    hit = read_cached(cached)
    if hit is not None:
        return hit

    profile = work_dir / f"{sample}.singlem.profile.tsv"
    otu = work_dir / f"{sample}.singlem.otu_table.tsv"
    cmd = [str(exe), "pipe", "-1", str(r1)]
    if r2 is not None:
        cmd += ["-2", str(r2)]
    cmd += ["--metapackage", str(_pkg_dir(refs_dir)), "-p", str(profile), "--otu-table", str(otu),
            "--threads", str(threads)]
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not profile.is_file():
        raise RuntimeError(f"singlem pipe failed (exit {proc.returncode}): {proc.stderr[-3000:]}")
    observations, summary = _parse(profile)
    result = LaneResult(
        lane_id=LANE_ID, tool="singlem", tool_version=version, reference_release_id=DB_RELEASE,
        taxonomy_release="GlobDB r232 metapackage", observations=observations,
        elapsed_s=time.monotonic() - t0, cached=False, command=tuple(cmd), input_sha256=sha,
        raw_result_uri=str(profile), summary=summary,
        notes=["coverage is SingleM's mean marker coverage; a lineage placed above species is unresolved, not absent",
               f"unresolved lineages at >= {ASSEMBLY_COVERAGE}x trigger targeted assembly"],
    )
    cached.write_text(json.dumps(result.to_json(), indent=1))
    return result


def _deepest_rank(lineage: str) -> str:
    segs = [s.strip() for s in lineage.split(";") if s.strip()]
    if not segs:
        return ""
    last = segs[-1]
    return {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}.get(
        last[:1], "")


def _parse(profile: Path) -> tuple[list[Observation], dict[str, Any]]:
    obs: list[Observation] = []
    unresolved_cov = 0.0
    total_cov = 0.0
    triggers: list[dict[str, Any]] = []
    with profile.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        cov_i = idx.get("coverage", 1)
        tax_i = idx.get("taxonomy", 2)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(cov_i, tax_i):
                continue
            try:
                cov = float(parts[cov_i])
            except ValueError:
                continue
            if cov <= 0:
                continue
            lineage = parts[tax_i].replace("Root; ", "").replace("Root;", "").strip()
            lineage = ";".join(s.strip() for s in lineage.split(";") if s.strip())
            rank = _deepest_rank(lineage)
            total_cov += cov
            resolved = rank == "species"
            if not resolved:
                unresolved_cov += cov
                if cov >= ASSEMBLY_COVERAGE:
                    triggers.append({"lineage": lineage, "rank": rank, "coverage": round(cov, 2)})
            obs.append(Observation(
                native_id=f"singlem:{lineage.rsplit(';', 1)[-1] or 'root'}", lineage=lineage,
                rank=rank or "root", status="supported" if resolved else "ambiguous",
                abundance_value=cov, abundance_unit="SingleM mean marker coverage (x)", denominator="marker windows",
                support_metrics={"resolved_to_species": resolved, "assembly_trigger": (not resolved and cov >= ASSEMBLY_COVERAGE)},
            ))
    obs.sort(key=lambda o: -(o.abundance_value or 0.0))
    return obs, {
        "total_coverage": round(total_cov, 2), "unresolved_coverage": round(unresolved_cov, 2),
        "unresolved_fraction": round(unresolved_cov / total_cov, 4) if total_cov else None,
        "n_lineages": len(obs), "n_species_resolved": sum(1 for o in obs if o.rank == "species"),
        "assembly_triggers": triggers,
    }
