#!/usr/bin/env python3
"""Build the gut rescue panel as a Kraken2 + Bracken database (spec 0.8.4 §4B).

The panel is what the reconciliation admitted beyond GTDB/UHGG: GlobDB
non-GTDB clusters from gut sources (UHGG, HRGM2), ELGG member genomes
that add strain diversity to known species, and the food/background
competitors GlobDB carries (cFMD) once their genomes are on disk. One
explicit custom taxonomy, built here from each genome's lineage, with the
mapping back to source identifiers written beside it - reversible, never
an integer taxonomy borrowed from another database.

    scripts/build_rescue_panel_kraken.py [--threads N] [--max-genomes N]

Writes refs/kraken2/rescue_panel/ with hash.k2d etc., taxonomy/, a
`panel_manifest.tsv` (taxid <-> source id <-> lineage <-> category) and a
`build.json` lock. Bracken distributions are built for 150 and 100 bp.
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
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from openbiota.expansion import genomes, globdb  # noqa: E402

MANIFEST = REPO / "refs" / "expanded" / "reconciliation" / "rescue_panel.tsv"
OUT = REPO / "refs" / "kraken2" / "rescue_panel"
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


class Taxonomy:
    """A minimal NCBI-style taxonomy built from lineage strings."""

    def __init__(self) -> None:
        self.nodes: dict[int, tuple[int, str, str]] = {1: (1, "no rank", "root")}
        self.by_key: dict[tuple[str, ...], int] = {(): 1}
        self.next_id = 2

    def add(self, lineage: str) -> int:
        segs = [x.strip() for x in lineage.split(";") if x.strip()]
        key: tuple[str, ...] = ()
        parent = 1
        for i, seg in enumerate(segs[:7]):
            key = key + (seg,)
            tid = self.by_key.get(key)
            if tid is None:
                tid = self.next_id
                self.next_id += 1
                self.by_key[key] = tid
                self.nodes[tid] = (parent, RANKS[i] if i < len(RANKS) else "no rank", seg.split("__", 1)[-1] or seg)
            parent = tid
        return parent

    def write(self, tax_dir: Path) -> None:
        tax_dir.mkdir(parents=True, exist_ok=True)
        with (tax_dir / "nodes.dmp").open("w") as n, (tax_dir / "names.dmp").open("w") as m:
            for tid, (parent, rank, name) in sorted(self.nodes.items()):
                n.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\t-\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\n")
                m.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")


def _retag(src: Path, dest: Path, taxid: int, gid: str) -> int:
    """Copy a genome with `kraken:taxid` headers; returns the number of contigs."""
    opener = gzip.open if src.suffix == ".gz" else open
    n = 0
    with opener(src, "rt") as fh, dest.open("w") as out:  # type: ignore[arg-type]
        for line in fh:
            if line.startswith(">"):
                n += 1
                out.write(f">{gid}__{n}|kraken:taxid|{taxid}\n")
            else:
                out.write(line)
    return n


def _retag_limited(src: Path, dest: Path, taxid: int, gid: str, max_bases: int) -> int:
    """Like _retag but stops after max_bases (0 = whole file); for large decoy genomes."""
    opener = gzip.open if src.suffix == ".gz" else open
    n = taken = 0
    with opener(src, "rt") as fh, dest.open("w") as out:  # type: ignore[arg-type]
        for line in fh:
            if line.startswith(">"):
                if max_bases and taken >= max_bases:
                    break
                n += 1
                out.write(f">{gid}__{n}|kraken:taxid|{taxid}\n")
            else:
                out.write(line)
                taken += len(line) - 1
    return n


_ACC_SPECIES: dict[str, str] | None = None


def _lineage_for_accession(acc: str, reps: dict) -> str:
    """GTDB R232 lineage for an NCBI accession or GlobDB id.

    GlobDB representatives carry their lineage; any other accession is
    looked up in the R232 species clusters (every member genome is listed
    there) and given a species-only lineage under its genus."""
    global _ACC_SPECIES
    if acc in reps:
        return reps[acc].lineage
    base = acc.split(".")[0]
    if _ACC_SPECIES is None:
        from openbiota.expansion import gtdb as _gtdb

        _ACC_SPECIES = {}
        for sp, cluster in _gtdb.species_clusters("r232").items():
            for member in (cluster.representative, *cluster.members):
                _ACC_SPECIES[str(member).split(".")[0].removeprefix("RS_").removeprefix("GB_")] = sp
    sp = _ACC_SPECIES.get(base.removeprefix("RS_").removeprefix("GB_"))
    if not sp:
        return ""
    # The same species must map to one taxid: reuse a representative's full
    # lineage when GlobDB has the species, else a species-only lineage.
    return _species_lineage(sp, reps)


_SPECIES_LINEAGE: dict[str, str] | None = None


def _species_lineage(sp: str, reps: dict) -> str:
    """One lineage per species: a GlobDB representative's, else GTDB R232's own."""
    global _SPECIES_LINEAGE
    if _SPECIES_LINEAGE is None:
        _SPECIES_LINEAGE = {}
        for rep in reps.values():
            if rep.species and rep.species not in _SPECIES_LINEAGE:
                _SPECIES_LINEAGE[rep.species] = rep.lineage
        from openbiota.expansion import gtdb as _gtdb

        for name, cluster in _gtdb.species_clusters("r232").items():
            if name not in _SPECIES_LINEAGE and cluster.taxonomy:
                _SPECIES_LINEAGE[name] = cluster.taxonomy
    return _SPECIES_LINEAGE.get(sp) or f"d__Bacteria;p__;c__;o__;f__;g__{sp.split(' ')[0]};s__{sp}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--max-genomes", type=int, default=None)
    ap.add_argument("--no-neighbours", action="store_true", help="targets only (not recommended: misassigns relatives)")
    ap.add_argument("--decoy-bases", type=int, default=200_000_000,
                    help="bases taken from each food/background genome (0 = whole genome)")
    ap.add_argument("--nice", type=int, default=0, help="run kraken2-build/bracken-build at this niceness")
    args = ap.parse_args()
    if not MANIFEST.is_file():
        print("no rescue panel manifest; run `make expansion` first")
        return 1
    kraken_build = shutil.which("kraken2-build")
    bracken_build = shutil.which("bracken-build")
    if not (kraken_build and bracken_build):
        print("kraken2-build / bracken-build not installed")
        return 1

    reps = globdb.representatives()
    rows = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
    if args.max_genomes:
        rows = rows[: args.max_genomes]
    OUT.mkdir(parents=True, exist_ok=True)
    lib = OUT / "library_src"
    lib.mkdir(exist_ok=True)
    tax = Taxonomy()
    manifest_rows: list[dict] = []
    n_ok = n_pending = n_failed = 0
    t0 = time.monotonic()
    for r in rows:
        gid = r["globdb_id"] or r["source_id"]
        lineage = ""
        if gid in reps:
            lineage = reps[gid].lineage
        elif r.get("canonical_species"):
            sp = r["canonical_species"]
            lineage = _species_lineage(sp, reps)
        if not lineage:
            manifest_rows.append({**r, "taxid": "", "status": "no_lineage"})
            continue
        path: Path | None = None
        if r.get("local_path") and Path(r["local_path"]).is_file():
            path = Path(r["local_path"])
        else:
            fetch_id = gid if (gid in reps or gid.startswith(("GC", "MGYG", "HRGM", "HROM"))) else r["source_id"]
            try:
                path = genomes.fetch(fetch_id)
            except genomes.Pending:
                n_pending += 1
                manifest_rows.append({**r, "taxid": "", "status": "pending_archive"})
                continue
            except Exception as exc:  # noqa: BLE001 - one genome must not stop the build
                n_failed += 1
                manifest_rows.append({**r, "taxid": "", "status": f"fetch_failed: {type(exc).__name__}"})
                continue
        taxid = tax.add(lineage)
        dest = lib / f"{gid}.fna"
        if not dest.is_file():
            _retag(path, dest, taxid, gid)
        manifest_rows.append({**r, "taxid": str(taxid), "lineage": lineage, "status": "in_panel"})
        n_ok += 1
    # Near neighbours and already-detected alternatives (spec §4 item 5): a
    # panel holding only the additions would hand every read of a common
    # gut species to whatever novel cluster of that genus is present. Every
    # GlobDB representative from a gut source in a target's genus, and every
    # reference genome any sample's confirmation has fetched, compete too.
    target_genera = {reps[g].genus for g in (r["globdb_id"] or r["source_id"] for r in rows) if g in reps}
    # Gut-relevant species: every canonical species the supplement
    # reconciliation touched (UHGG, HRGM2, ELGG, HumGut2, HROM) plus every
    # GTDB species any sample's inventory has held. GTDB representatives of
    # those species in a target genus are the near neighbours that keep a
    # common species' reads from landing on a novel cluster of its genus;
    # all 19k GTDB species in those genera would triple the build for
    # organisms that never occur in a gut.
    gut_species: set[str] = set()
    for path in (REPO / "refs" / "expanded" / "reconciliation").glob("*.tsv.gz"):
        with gzip.open(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("canonical_species"):
                    gut_species.add(row["canonical_species"])
    for res in (REPO / "results").glob("*/results.json"):
        try:
            for o in (json.loads(res.read_text()).get("organism_inventory") or {}).get("organisms") or []:
                if o.get("gtdb"):
                    gut_species.add(str(o["gtdb"]).split(" / ")[0])
        except (OSError, json.JSONDecodeError):
            continue
    gut_sources = {"UHGG", "HRGM", "HRGM2", "MGnify"}
    n_neighbours = n_alternatives = 0
    if not args.no_neighbours:
        for gid, rep in reps.items():
            if rep.genus not in target_genera:
                continue
            if not (rep.source in gut_sources or (rep.source == "GTDB" and rep.species in gut_species)):
                continue
            if True:
                dest = lib / f"{gid}.fna"
                if dest.is_file():
                    n_neighbours += 1
                    continue
                try:
                    path = genomes.fetch(gid)
                except Exception:  # noqa: BLE001 - archive still unpacking, or not a member
                    continue
                _retag(path, dest, tax.add(rep.lineage), gid)
                manifest_rows.append({"source": rep.source, "source_id": gid, "globdb_id": gid, "category": "near_neighbour",
                                      "canonical_species": rep.species, "taxid": str(tax.add(rep.lineage)),
                                      "lineage": rep.lineage, "status": "in_panel", "role": "competitor"})
                n_neighbours += 1
        for fna in sorted((REPO / "refs" / "genomes").glob("*.fna.gz")):
            acc = fna.name.replace(".fna.gz", "")
            dest = lib / f"{acc}.fna"
            if dest.is_file():
                continue
            lineage = _lineage_for_accession(acc, reps)
            if not lineage:
                continue
            _retag(fna, dest, tax.add(lineage), acc)
            manifest_rows.append({"source": "detected_alternative", "source_id": acc, "globdb_id": "", "category": "detected_alternative",
                                  "canonical_species": lineage.split("s__")[-1], "taxid": str(tax.add(lineage)),
                                  "lineage": lineage, "status": "in_panel", "role": "competitor"})
            n_alternatives += 1
    # Decoys: PhiX and the food/background genomes on disk. Host reads are
    # removed upstream by the bowtie2 host filter, so the human genome is
    # not repeated here.
    n_decoys = 0
    for name, src in (("phiX174", REPO / "refs" / "host" / "phix.fna.gz"),
                      ("Bos_taurus", REPO / "refs" / "decoys" / "bos_taurus" / "genome.fna"),
                      ("Gallus_gallus", REPO / "refs" / "decoys" / "gallus_gallus" / "genome.fna"),
                      ("Sus_scrofa", REPO / "refs" / "decoys" / "sus_scrofa" / "genome.fna")):
        if not src.is_file():
            continue
        lineage = f"d__Background;p__;c__;o__;f__;g__{name.split('_')[0]};s__{name}"
        dest = lib / f"decoy_{name}.fna"
        if not dest.is_file():
            _retag_limited(src, dest, tax.add(lineage), f"decoy_{name}", max_bases=args.decoy_bases)
        manifest_rows.append({"source": "decoy", "source_id": name, "globdb_id": "", "category": "background_decoy",
                              "canonical_species": name, "taxid": str(tax.add(lineage)), "lineage": lineage,
                              "status": "in_panel", "role": "competitor",
                              "note": f"first {args.decoy_bases:,} bases" if args.decoy_bases else "whole genome"})
        n_decoys += 1
    print(f"panel genomes: {n_ok} targets in, {n_neighbours} near neighbours, {n_alternatives} detected alternatives, "
          f"{n_decoys} decoys; {n_pending} pending the GlobDB archive, {n_failed} failed to fetch "
          f"({time.monotonic() - t0:.0f}s)")
    tax.write(OUT / "taxonomy")
    with (OUT / "panel_manifest.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in manifest_rows for k in r}), delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)
    if n_ok == 0:
        print("nothing to build")
        return 1

    # Kraken2 build over the retagged library, then Bracken distributions.
    # Everything in one library file: add-to-library per genome forks a
    # process per file and copies it; a single concatenation is one pass.
    t1 = time.monotonic()
    nice = ["nice", "-n", str(args.nice)] if args.nice else []
    combined = OUT / "library_all.fna"
    if not (OUT / "library" / "added").is_dir():
        with combined.open("w") as out:
            for fna in sorted(lib.glob("*.fna")):
                with fna.open() as fh:
                    shutil.copyfileobj(fh, out)
        subprocess.run([*nice, kraken_build, "--add-to-library", str(combined), "--db", str(OUT), "--no-masking"],
                       check=True, capture_output=True, text=True)
        combined.unlink(missing_ok=True)
    if not (OUT / "hash.k2d").is_file():
        proc = subprocess.run([*nice, kraken_build, "--build", "--db", str(OUT), "--threads", str(args.threads)],
                              capture_output=True, text=True, check=False)
        if "OMP only wants you to use 1 threads" in (proc.stderr or "") + (proc.stdout or ""):
            # This platform's kraken2 was built without OpenMP; the hash build
            # is single-threaded here and refuses anything else.
            print("kraken2 build_db has no OpenMP on this platform; building with one thread")
            proc = subprocess.run([*nice, kraken_build, "--build", "--db", str(OUT), "--threads", "1"],
                                  capture_output=True, text=True, check=False)
        if proc.returncode != 0 or not (OUT / "hash.k2d").is_file():
            print("kraken2-build --build failed:\n" + (proc.stderr or proc.stdout)[-3000:])
            return 1
    for rl in (150, 100):
        if (OUT / f"database{rl}mers.kmer_distrib").is_file():
            continue
        bp = subprocess.run([*nice, bracken_build, "-d", str(OUT), "-t", str(args.threads), "-k", "35", "-l", str(rl)],
                            check=False, capture_output=True, text=True)
        if bp.returncode != 0:
            print(f"bracken-build -l {rl} failed:\n" + (bp.stderr or bp.stdout)[-1500:])
    build_s = time.monotonic() - t1
    lock = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_genomes": n_ok, "n_near_neighbours": n_neighbours, "n_detected_alternatives": n_alternatives, "n_decoys": n_decoys,
        "n_pending_archive": n_pending, "n_fetch_failed": n_failed,
        "n_taxa": len(tax.nodes), "kraken2": subprocess.run(["kraken2", "--version"], capture_output=True, text=True).stdout.strip().splitlines()[0],
        "build_seconds": round(build_s, 1), "manifest": str(OUT / "panel_manifest.tsv"),
        "taxonomy": "custom, built from GlobDB/GTDB R232 lineages; taxid <-> source id in panel_manifest.tsv",
        "bracken_read_lengths": [150, 100],
        "sources_in_panel": sorted({r["source"] for r in manifest_rows if r.get("status") == "in_panel"}),
    }
    (OUT / "build.json").write_text(json.dumps(lock, indent=2) + "\n")
    print(json.dumps(lock, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
