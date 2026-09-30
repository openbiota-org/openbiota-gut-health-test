#!/usr/bin/env python3
"""Place the UHGG species representatives no accession or name could, by sequence.

Spec 0.8.4 §3: a supplement representative that matches nothing by
accession or name is compared to the reference catalogue by ANI, and is
`already_represented` when it is the same species as a catalogue genome.
The R6 reconciliation left 470 UHGG v2.0.2 representatives `unresolved_
mapping`: MGnify MAGs GTDB never took in, with no species name of their
own. On the simulated communities these turned into "new clusters" that
were in fact truth species under another catalogue's identifier - counted
as false positives in the benchmark and listed as unnamed clusters in the
reports.

For each such representative the genome is sketched as a sample and
queried against the GlobDB r232 sylph database (346,232 genomes), which
shortlists candidates by k-mer containment; the best hit's naive ANI and
containment decide. Same species: naive ANI >= 95.5 with at least a quarter
of the genome's k-mers contained. Boundary hits (94-95.5) are recorded and
left unresolved, not forced.

    scripts/reconcile_uhgg_ani.py [--threads N] [--limit N]

Rewrites refs/expanded/reconciliation/uhgg_v2.0.2.tsv.gz in place (a copy
of the previous file is kept beside it) and writes uhgg_v2.0.2.ani.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from openbiota.expansion import genomes  # noqa: E402

RECON = ROOT / "refs" / "expanded" / "reconciliation" / "uhgg_v2.0.2.tsv.gz"
SYLPH = ROOT / "vendor" / "sylph" / "bin" / "sylph"
DB = ROOT / "refs" / "sylph" / "globdb_r232" / "globdb_r232_sylph_v2_c200.syl2db"
TAXONOMY = ROOT / "refs" / "sylph" / "globdb_r232" / "globdb_r232_taxonomy_sylph.tsv.gz"
ANI_SAME: float = 95.5
ANI_BOUNDARY: float = 94.0
MIN_CONTAINMENT: float = 0.25


def _taxonomy() -> dict[str, str]:
    out: dict[str, str] = {}
    with gzip.open(TAXONOMY, "rt") as fh:
        for line in fh:
            name, _, tax = line.rstrip("\n").partition("\t")
            out[name.rsplit("/", 1)[-1]] = tax
    return out


def _species(tax: str) -> str:
    for seg in tax.split(";"):
        if seg.startswith("s__"):
            return seg[3:].strip()
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--categories", default="unresolved_mapping")
    a = ap.parse_args()
    with gzip.open(RECON, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    wanted = [r for r in rows if r["category"] in set(a.categories.split(","))]
    if a.limit:
        wanted = wanted[: a.limit]
    print(f"{len(wanted)} representatives to place by sequence", flush=True)

    paths: dict[str, Path] = {}

    def _get(r: dict[str, str]) -> None:
        try:
            paths[r["source_id"]] = genomes.fetch(r["source_id"])
        except Exception as exc:  # noqa: BLE001
            print(f"  {r['source_id']}: {type(exc).__name__}: {str(exc)[:80]}", flush=True)

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_get, wanted))
    print(f"{len(paths)} genomes on disk", flush=True)
    if not paths:
        return 1

    with tempfile.TemporaryDirectory(prefix="uhgg_ani_") as tmp:
        sketch_dir = Path(tmp) / "sk"
        subprocess.run([str(SYLPH), "sketch", "-r", *[str(p) for p in paths.values()], "-d", str(sketch_dir),
                        "-t", str(a.threads)], check=True, capture_output=True)
        sketches = sorted(sketch_dir.glob("*.sylsp"))
        out = Path(tmp) / "query.tsv"
        subprocess.run([str(SYLPH), "query", str(DB), *[str(s) for s in sketches], "-t", str(a.threads),
                        "-o", str(out)], check=True, capture_output=True)
        best: dict[str, tuple[float, float, str]] = {}
        with out.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for h in reader:
                sid = Path(h["Sample_file"]).name.split(".")[0]
                naive = float(h["Naive_ANI"])
                num, _, den = h["Containment_ind"].partition("/")
                cont = int(num) / int(den) if den and int(den) else 0.0
                genome = h["Genome_file"].rsplit("/", 1)[-1]
                if sid not in best or (naive, cont) > best[sid][:2]:
                    best[sid] = (naive, cont, genome)

    tax = _taxonomy()
    n_same = n_boundary = n_none = 0
    by_id = {r["source_id"]: r for r in rows}
    placed: dict[str, dict] = {}
    for sid in paths:
        hit = best.get(sid)
        r = by_id[sid]
        if hit is None:
            n_none += 1
            placed[sid] = {"result": "no_hit"}
            continue
        naive, cont, genome = hit
        species = _species(tax.get(genome, ""))
        gid = genomes.normalise(genome)
        rec = {"globdb_genome": gid, "naive_ani": round(naive, 2), "containment": round(cont, 4), "species": species}
        if naive >= ANI_SAME and cont >= MIN_CONTAINMENT and species:
            r["category"], r["method"] = "already_represented", "ani_sylph"
            r["canonical_species"], r["globdb_id"] = species, gid
            r["ani"], r["af_shorter"] = f"{naive:.2f}", f"{100 * cont:.1f}"
            rec["result"] = "same_species"
            n_same += 1
        elif naive >= ANI_BOUNDARY:
            r["method"] = "ani_sylph_boundary"
            r["ani"], r["af_shorter"] = f"{naive:.2f}", f"{100 * cont:.1f}"
            rec["result"] = "boundary"
            n_boundary += 1
        else:
            rec["result"] = "distinct"
            n_none += 1
        placed[sid] = rec

    backup = RECON.with_suffix(".before_ani.tsv.gz")
    if not backup.is_file():
        shutil.copyfile(RECON, backup)
    with gzip.open(RECON, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summary = {
        "run_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"), "tool": "sylph 1.0.0 (naive ANI, containment)",
        "database": DB.name, "n_queried": len(paths), "n_same_species": n_same, "n_boundary": n_boundary,
        "n_distinct_or_no_hit": n_none, "thresholds": {"naive_ani_same": ANI_SAME, "naive_ani_boundary": ANI_BOUNDARY,
                                                        "min_containment": MIN_CONTAINMENT},
        "placements": placed,
    }
    (RECON.parent / "uhgg_v2.0.2.ani.json").write_text(json.dumps(summary, indent=1))
    print(f"same species {n_same}; boundary {n_boundary}; distinct or no hit {n_none}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
