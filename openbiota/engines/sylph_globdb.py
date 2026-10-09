"""GlobDB r232 discovery lane on sylph (spec 0.8.4 R3).

The installed GTDB R232 sylph lane is the baseline and is untouched
(`openbiota.engines.sylph`). This lane profiles the same sample sketch
against GlobDB r232's two-stage index: 346,233 species representatives
from 26 source datasets, of which 199,923 are GTDB and 146,310 are
non-GTDB clusters under GlobDB's own rules (96% ANI / 50% AF, GTDB
priority). Those 146,310 are not "146,310 more gut species"; the lane
records GlobDB's identifiers and lineage as published and the crosswalk
decides what each one is against GTDB.

Provenance is kept separate from the GTDB lane: different database,
different taxonomy file, different lock. The sample sketch (`-c 200`) is
identical for both and is reused, as the spec allows.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Final

from openbiota.engines.lanebase import LaneResult, Observation, cache_key, input_sha256, read_cached

LANE_ID: Final = "sylph_globdb"
SYLPH_VERSION: Final = "1.0.0"
DB_RELEASE: Final = "GlobDB r232"
DATABASE_FILE: Final = "globdb_r232_sylph_v2_c200.syl2db"
TAXONOMY_FILE: Final = "globdb_r232_taxonomy_sylph.tsv.gz"
CATALOGUE_SIZE: Final = 346_233
NON_GTDB_REPRESENTATIVES: Final = 146_310
MINIMUM_ANI: Final = 95.0
SKETCH_C: Final = 200
CACHE_VERSION: Final = 1

_ACC = re.compile(r"(GC[AF]_\d+(?:\.\d+)?)")


def _sylph() -> Path | None:
    root = Path(__file__).resolve().parents[2]
    exe = root / "vendor" / "sylph" / "bin" / "sylph"
    if exe.is_file():
        return exe
    found = shutil.which("sylph")
    return Path(found) if found else None


def _db_dir(refs_dir: Path) -> Path:
    return refs_dir / "sylph" / "globdb_r232"


def database_ready(refs_dir: Path) -> bool:
    d = _db_dir(refs_dir)
    if not ((d / DATABASE_FILE).is_file() and (d / TAXONOMY_FILE).is_file()):
        return False
    lock = Path("refs/expanded/locks/globdb_r232.lock.json")
    if lock.is_file():
        rec = json.loads(lock.read_text()).get("files", {}).get(f"sylph/globdb_r232/{DATABASE_FILE}")
        return bool(rec and rec.get("verified"))
    return True


def available(refs_dir: Path) -> bool:
    return _sylph() is not None and database_ready(refs_dir)


def genome_key(name: str) -> str:
    """`GCF_003697165.fa.gz` / `MGYG000001234.fa.gz` / `.../HRGMV2_00123.fa` -> bare id."""
    head = (name.split()[0] if name else "").rsplit("/", 1)[-1]
    for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta", "_genomic"):
        if head.endswith(suffix):
            head = head[: -len(suffix)]
    m = _ACC.search(head)
    return m.group(1) if m else head


def taxonomy(refs_dir: Path) -> dict[str, str]:
    """GlobDB genome id -> lineage, from the publisher's sylph taxonomy file."""
    out: dict[str, str] = {}
    with gzip.open(_db_dir(refs_dir) / TAXONOMY_FILE, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and "d__" in parts[1]:
                out[genome_key(parts[0])] = parts[1]
    return out


def cache_file(work_dir: Path, r1: Path, r2: Path | None) -> Path:
    """Where this lane's result for these reads is cached (see metaphlan_jan26.cache_file)."""
    sha = input_sha256([r1, r2])
    key = cache_key(str(CACHE_VERSION), SYLPH_VERSION, DATABASE_FILE, TAXONOMY_FILE, str(MINIMUM_ANI), str(SKETCH_C), sha)
    return work_dir / f"sylph_globdb.{key}.json"


def run_sylph_globdb(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    work_dir: Path,
    threads: int = 8,
) -> LaneResult | None:
    exe = _sylph()
    if exe is None or not database_ready(refs_dir):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    sha = input_sha256([r1, r2])
    cached = cache_file(work_dir, r1, r2)
    hit = read_cached(cached)
    if hit is not None:
        return hit

    t0 = time.monotonic()
    # The GTDB lane's sketch lives in genome/sketches; reuse it when present.
    sketch_dir = work_dir.parent / "genome" / "sketches"
    sketch = sketch_dir / f"{sample}.paired.sylsp"
    if not sketch.is_file():
        sketch_dir.mkdir(parents=True, exist_ok=True)
        cmd = [str(exe), "sketch", "-1", str(r1)]
        if r2 is not None:
            cmd += ["-2", str(r2)]
        cmd += ["-S", sample, "-c", str(SKETCH_C), "-t", str(threads), "-d", str(sketch_dir)]
        subprocess.run(cmd, check=True, capture_output=True)
    profile_tsv = work_dir / f"{sample}.sylph_globdb.tsv"
    cmd = [str(exe), "profile", "-d", str(_db_dir(refs_dir) / DATABASE_FILE), str(sketch),
           "--minimum-ani", str(MINIMUM_ANI), "-u", "-t", str(threads), "-o", str(profile_tsv)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not profile_tsv.is_file():
        raise RuntimeError(f"sylph (GlobDB) failed (exit {proc.returncode}): {proc.stderr[-2000:]}")
    lineages = taxonomy(refs_dir)
    observations, summary = _parse(profile_tsv, lineages)
    result = LaneResult(
        lane_id=LANE_ID, tool="sylph", tool_version=SYLPH_VERSION, reference_release_id=DB_RELEASE,
        taxonomy_release="GlobDB r232 (GTDB R232 for GTDB members; GlobDB-native clusters otherwise)",
        observations=observations, elapsed_s=time.monotonic() - t0, cached=False, command=tuple(cmd),
        input_sha256=sha, raw_result_uri=str(profile_tsv), summary=summary,
        notes=[f"{CATALOGUE_SIZE:,} representatives, {NON_GTDB_REPRESENTATIVES:,} of them non-GTDB clusters under "
               "GlobDB's 96% ANI / 50% AF rule; a non-GTDB hit is a GlobDB unit until the crosswalk says otherwise",
               f"minimum adjusted ANI {MINIMUM_ANI}; sequence and taxonomic abundance are sylph's own units"],
    )
    cached.write_text(json.dumps(result.to_json(), indent=1))
    return result


def _f(v: Any) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _parse(profile_tsv: Path, lineages: dict[str, str]) -> tuple[list[Observation], dict[str, Any]]:
    obs: list[Observation] = []
    n_gtdb = n_native = 0
    with profile_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = genome_key(row.get("Genome_file") or row.get("Contig_name") or "")
            lineage = lineages.get(gid, "")
            if not lineage:
                continue
            is_gtdb = bool(_ACC.match(gid))
            n_gtdb += is_gtdb
            n_native += not is_gtdb
            species = ""
            for seg in lineage.split(";"):
                if seg.startswith("s__"):
                    species = seg[3:].strip()
            lam = str(row.get("Eff_lambda") or "").strip()
            ani = _f(row.get("Adjusted_ANI"))
            # sylph's own uncertainty: LOW lambda means too little signal for
            # a confident ANI; those are provisional until confirmation.
            status = "supported" if (lam not in ("LOW", "") and (ani or 0) >= MINIMUM_ANI) else "provisional"
            obs.append(Observation(
                native_id=f"globdb:{gid}", lineage=lineage, rank="species" if species else "genus", status=status,
                abundance_value=_f(row.get("Taxonomic_abundance")), abundance_unit="sylph taxonomic abundance (percent)",
                denominator="sylph-profiled genomes", support_metrics={
                    "adjusted_ani": ani, "sequence_abundance": _f(row.get("Sequence_abundance")),
                    "true_cov": _f(row.get("True_cov")), "eff_lambda": lam, "gtdb_member": is_gtdb,
                    "kmers_reassigned": _f(row.get("kmers_reassigned")),
                },
            ))
    obs.sort(key=lambda o: -(o.abundance_value or 0.0))
    return obs, {"n_hits": len(obs), "n_gtdb_representatives": n_gtdb, "n_globdb_native_clusters": n_native}
