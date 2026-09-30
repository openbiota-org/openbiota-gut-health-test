#!/usr/bin/env python3
"""Inventory of everything under refs/, for packaging and installation elsewhere.

Every reference the pipeline reads lives under `refs/`. Some of it is
downloaded from a publisher and covered by a lock (`refs/expanded/locks`,
with URL, size and SHA-256), some is built here from those downloads
(indexes, crosswalks, reconciliation, marker caches), some is a cohort
profiled locally, and some is a cache of genomes fetched on demand. To
archive the repo and install it on another machine, one needs to know
which is which, how big each part is, and how each part is regenerated.

This script writes `refs/MANIFEST.json` and `refs/MANIFEST.md` with one
entry per top-level directory: bytes, file count, category, the lock or
manifest that covers it, and the command that rebuilds it. It reads file
metadata only - no hashing; the locks already hold the checksums.

    scripts/refs_manifest.py            # write the manifest
    scripts/refs_manifest.py --check    # exit 1 if a lock-covered file is missing or the wrong size
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REFS = REPO / "refs"
LOCKS = REFS / "expanded" / "locks"

#: directory -> (category, how it is rebuilt)
KNOWN: dict[str, tuple[str, str]] = {
    "expanded": ("locks + reconciliation (derived)", "make refs-expanded; make expansion"),
    "crosswalk": ("derived from GTDB/GlobDB/MetaPhlAn bridges", "make expansion-crosswalks"),
    "gtdb": ("downloaded (locked)", "make refs-expanded"),
    "globdb_r232": ("downloaded (locked); genomes/ unpacked from the archive", "make refs-expanded; make refs-globdb-genomes"),
    "metaphlan4_db": ("downloaded by MetaPhlAn (Jan26 + Jun23 indexes); marker dump and clade caches derived",
                      "make refs-expanded (metaphlan --install); caches rebuild on first use"),
    "metaphlan_db": ("downloaded by MetaPhlAn 3 (scoring lane)", "make setup"),
    "motus": ("downloaded (locked), unpacked", "make refs-expanded"),
    "kraken2": ("downloaded (locked) publisher Kraken2/Bracken database", "make refs-expanded"),
    "singlem": ("downloaded (locked) GlobDB metapackage", "make refs-expanded"),
    "sylph": ("downloaded GTDB R232 sylph database", "make genome-lane"),
    "genomes": ("cache of reference genomes fetched on demand (NCBI/UHGG/GlobDB); each has a .source sidecar",
                "refilled on demand by openbiota.expansion.genomes"),
    "expansion_cohort": ("profiled locally from public reads; manifest inside each cohort file",
                         "make expansion-cohorts"),
    "cohort": ("pathway reference cohort work dir (reads deleted after use)", "openbiota cohort"),
    "cohort_random": ("pathway reference cohort, uniform-random sampling work dir (retained reads; regenerable)", "openbiota cohort --sampling random"),
    "age.bak_rclr": ("superseded age model artefacts (backup; not read)", "none"),
    "reference_cohort_samples.json": ("pathway cohort sample table", "openbiota cohort"),
    "reference_ranges.json": ("pathway reference ranges", "openbiota cohort"),
    "reference_ranges.prefix600k.json": ("superseded pathway ranges (prefix sampling; kept for audit)", "none"),
    "taxonomic_cohort.json": ("scoring cohort matrix (curatedMetagenomicData, MetaPhlAn 3)", "openbiota taxonomic-cohort"),
    "taxonomic_cohort_manifest.json": ("scoring cohort manifest", "openbiota taxonomic-cohort"),
    "cmd": ("curatedMetagenomicData profiles + ExperimentHub index (scoring cohort)", "openbiota taxonomic-cohort"),
    "bacdive": ("BacDive API v2 response cache", "refilled on demand"),
    "micom": ("AGORA2 model library + community/LP caches", "make extension-refs"),
    "host": ("GRCh38 + PhiX bowtie2 host index", "openbiota build-host-index"),
    "decoys": ("food/background genomes", "make substrate-refs"),
    "pathogens": ("pathogen target references", "make refs-strain"),
    "strain_refs.lock.json": ("lock", "make strain-tools-lock"),
    "panels": ("built DIAMOND panel databases", "make build-db"),
    "db": ("built databases", "make build-db"),
    "age": ("age model artefacts", "scripts/train_*"),
    "mycobiome": ("mycobiome references", "make mycobiome-refs"),
    "validation": ("validation fixtures", "openbiota validate"),
    "cache": ("fast-load caches (regenerable)", "automatic"),
    "ctnpc": ("references", "make refs-ctnpc"),
    "supplements": ("supplement references", "make substrate-refs"),
    "nii": ("references", "make setup"),
}


def _du(path: Path) -> tuple[int, int]:
    total = count = 0
    if path.is_file():
        return path.stat().st_size, 1
    for p in path.rglob("*"):
        try:
            if p.is_file() and not p.is_symlink():
                total += p.stat().st_size
                count += 1
        except OSError:
            continue
    return total, count


def _locks() -> dict[str, dict]:
    covered: dict[str, dict] = {}
    if LOCKS.is_dir():
        for lock in sorted(LOCKS.glob("*.lock.json")):
            try:
                data = json.loads(lock.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            for dest, rec in (data.get("files") or {}).items():
                covered[dest] = {**rec, "source_id": data.get("source_id"), "lock": lock.name}
    return covered


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    covered = _locks()
    entries = []
    grand = 0
    for child in sorted(REFS.iterdir()):
        if child.name.startswith(".") or child.name.startswith("MANIFEST"):
            continue
        size, n = _du(child)
        grand += size
        cat, rebuild = KNOWN.get(child.name, ("unclassified", "see docs/EXPANDED_DETECTION.md"))
        locked = [d for d in covered if d.split("/", 1)[0] == child.name]
        entries.append({"path": f"refs/{child.name}", "bytes": size, "gb": round(size / 1e9, 2), "files": n,
                        "category": cat, "rebuild": rebuild, "locked_files": len(locked)})
    problems = []
    consumed = []
    for dest, rec in covered.items():
        p = REFS / dest
        if not p.is_file():
            # An installer that unpacks and removes its archive (MetaPhlAn's
            # --install) leaves the index where the archive was.
            stem = p.name.split(".")[0]
            if any(q.name.startswith(stem) and q.suffix in (".bt2l", ".pkl") for q in p.parent.glob(f"{stem}*")):
                consumed.append(f"{dest}: archive consumed by its installer; index present")
                continue
            problems.append(f"{dest}: missing (lock {rec['lock']})")
        elif rec.get("size") and p.stat().st_size != rec["size"]:
            problems.append(f"{dest}: size {p.stat().st_size} != lock {rec['size']}")
        elif not rec.get("verified"):
            problems.append(f"{dest}: not verified in lock")
    manifest = {
        "written_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "total_bytes": grand, "total_gb": round(grand / 1e9, 2),
        "n_lock_covered_files": len(covered), "lock_problems": problems, "consumed_archives": consumed,
        "entries": entries,
        "install": ["make tools-expanded", "make refs-expanded", "make expansion", "make refs-globdb-genomes",
                    "make expansion-cohorts", "scripts/fetch_expanded_refs.py --status", "scripts/refs_manifest.py --check"],
        "note": ("sizes and counts from file metadata; checksums are in refs/expanded/locks. Directories marked derived or "
                 "cache are regenerable and need not be archived; everything else is either locked-downloadable or profiled "
                 "from public reads with the recipe recorded in its own manifest."),
    }
    (REFS / "MANIFEST.json").write_text(json.dumps(manifest, indent=1) + "\n")
    md = [f"# refs/ inventory ({manifest['total_gb']} GB, {len(covered)} lock-covered files)", "",
          "| path | GB | files | category | rebuild |", "|---|---|---|---|---|"]
    for e in entries:
        md.append(f"| {e['path']} | {e['gb']} | {e['files']:,} | {e['category']} | `{e['rebuild']}` |")
    if problems:
        md += ["", "## lock problems", *[f"- {p}" for p in problems]]
    (REFS / "MANIFEST.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    if args.check and problems:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
