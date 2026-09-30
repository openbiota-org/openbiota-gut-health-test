"""Fungal reference preparation and the content-addressed reference lock.

Spec §4 and §6.2. The lock freezes what the lane searched against: every
admitted assembly with its accession version, sequence hash, organism,
taxid, lineage and genetic code; exact-duplicate alias groups; the decoy
set; the two deposited databases (EukDetect2, FungiGutDB) with their
hashes and licences; and the admitted counts kept separate (assemblies,
species, genera, marker-covered taxa). Published catalogue sizes are
provenance, not sensitivity.

Layout under ``refs/mycobiome/``::

    genomes/ncbi_dataset/data/<ACC>/<ACC>_*_genomic.fna   fungal assemblies
    decoys/ncbi_dataset/data/<ACC>/...                    gut bacterial decoys
    eukdetect2/                                           EukDetect2 database
    fungigut/FungiGut_db.tar.gz                           FungiGutDB v1.0
    sylph/fungi-c200.syldb                                whole-genome sketch DB
    manifests/                                            NCBI summaries used
    reference_lock.json                                   this module's output
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

PRESET: Final = "v08.2"
LOCK_VERSION: Final = 3
SKETCH_C: Final = 200          # must match the read sketches made for the GTDB lane
FUNGI_TAXID: Final = 4751
CGF_BIOPROJECT: Final = "PRJNA833221"

#: Audited upstream snapshots (spec §4).
UPSTREAM_COMMITS: Final = {
    "FungiGut": "dcadf124c7910ebda2cc0f7f9544f6709e44e3df",
    "EukDetect2": "8d69014727b2c5956de30b0811644b6c3a7bd4f3",
    "CGF": "353df4bf0fa31ba1146cd7e729cae8a90ceb5bd5",
    "FunOMIC": "d8ef2438a4e44ed6839e434778d1ca8e32ea3569",
}
LICENSES: Final = {
    "ncbi_genomes": "NCBI public sequence data; no licence asserted, cite the submitters",
    "cgf_genomes": "NCBI public sequence data (Cultivated Gut Fungi, Cell 2024); cite the paper",
    "eukdetect2_database": "CC BY 4.0 (Zenodo 19056625); code MIT",
    "fungigut_database": "CC BY 4.0 (Zenodo 17581472); workflow code GPL-3.0",
}


def root() -> Path:
    return Path(__file__).resolve().parents[2] / "refs" / "mycobiome"


@dataclass(frozen=True, slots=True)
class Assembly:
    accession: str
    organism: str
    taxid: int
    species_taxid: int | None
    species: str
    genus: str
    family: str
    kingdom: str
    level: str
    source: str                   # cgf | refseq_reference | panel | both
    path: str
    bases: int
    n_contigs: int
    sha256: str
    genetic_code: int | None
    mito_code: int | None
    strain: str = ""
    duplicate_of: str | None = None
    is_fungus: bool = True

    def to_json(self) -> dict[str, Any]:
        return {
            "accession": self.accession, "organism": self.organism, "taxid": self.taxid,
            "species_taxid": self.species_taxid, "species": self.species, "genus": self.genus,
            "family": self.family, "kingdom": self.kingdom, "level": self.level, "source": self.source,
            "path": self.path, "bases": self.bases, "n_contigs": self.n_contigs, "sha256": self.sha256,
            "genetic_code": self.genetic_code, "mito_code": self.mito_code, "strain": self.strain,
            "duplicate_of": self.duplicate_of, "is_fungus": self.is_fungus,
        }


@dataclass(slots=True)
class Lock:
    lock_id: str
    built_at: str
    assemblies: list[Assembly]
    decoys: list[dict[str, Any]]
    deposits: list[dict[str, Any]]
    taxonomy: dict[int, dict[str, Any]]
    counts: dict[str, int] = field(default_factory=dict)
    sketch_db: str | None = None
    licenses: dict[str, str] = field(default_factory=lambda: dict(LICENSES))
    upstream_commits: dict[str, str] = field(default_factory=lambda: dict(UPSTREAM_COMMITS))

    def by_accession(self) -> dict[str, Assembly]:
        return {a.accession: a for a in self.assemblies}

    def species_assemblies(self) -> dict[int, list[Assembly]]:
        out: dict[int, list[Assembly]] = defaultdict(list)
        for a in self.assemblies:
            if a.is_fungus and a.duplicate_of is None and a.species_taxid:
                out[a.species_taxid].append(a)
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "lock_version": LOCK_VERSION, "preset": PRESET, "lock_id": self.lock_id, "built_at": self.built_at,
            "counts": self.counts, "sketch_db": self.sketch_db, "sketch_c": SKETCH_C,
            "licenses": self.licenses, "upstream_commits": self.upstream_commits,
            "deposits": self.deposits, "decoys": self.decoys,
            "taxonomy": {str(k): v for k, v in self.taxonomy.items()},
            "assemblies": [a.to_json() for a in self.assemblies],
        }


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _sha256(path: Path, cache: dict[str, Any]) -> str:
    st = path.stat()
    key = str(path)
    hit = cache.get(key)
    if hit and hit.get("size") == st.st_size and hit.get("mtime") == int(st.st_mtime):
        return str(hit["sha256"])
    h = hashlib.sha256()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:  # type: ignore[operator]
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    digest = h.hexdigest()
    cache[key] = {"size": st.st_size, "mtime": int(st.st_mtime), "sha256": digest}
    return digest


def _fasta_stats(path: Path) -> tuple[int, int]:
    """(total bases, number of records) without loading the file."""
    bases = 0
    n = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith(b">"):
                n += 1
            else:
                bases += len(line.strip())
    return bases, n


def _accessions(path: Path) -> set[str]:
    with path.open() as fh:
        return {json.loads(line)["accession"] for line in fh if line.strip()}


def _datasets_bin() -> Path | None:
    p = Path(__file__).resolve().parents[2] / "vendor" / "ncbi-datasets" / "datasets"
    if p.is_file():
        return p
    found = shutil.which("datasets")
    return Path(found) if found else None


def _sylph_bin() -> Path | None:
    p = Path(__file__).resolve().parents[2] / "vendor" / "sylph" / "bin" / "sylph"
    if p.is_file():
        return p
    found = shutil.which("sylph")
    return Path(found) if found else None


def taxonomy_table(taxids: set[int], cache_path: Path) -> dict[int, dict[str, Any]]:
    """Lineage, current name and genetic codes for every taxid, via the NCBI
    datasets CLI, cached on disk. Missing taxids are recorded as unknown."""
    table: dict[int, dict[str, Any]] = {}
    if cache_path.is_file():
        table = {int(k): v for k, v in json.loads(cache_path.read_text()).items()}
    missing = sorted(t for t in taxids if t not in table)
    ds = _datasets_bin()
    if missing and ds is not None:
        ids = cache_path.with_suffix(".query.txt")
        ids.write_text("\n".join(str(t) for t in missing) + "\n")
        proc = subprocess.run([str(ds), "summary", "taxonomy", "taxon", "--inputfile", str(ids), "--as-json-lines"],
                              capture_output=True, text=True, check=False)
        for line in proc.stdout.splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            tax = rec.get("taxonomy") or {}
            cls = tax.get("classification") or {}
            def _name(rank: str, cls: dict[str, Any] = cls) -> str:
                return str((cls.get(rank) or {}).get("name") or "")

            def _id(rank: str, cls: dict[str, Any] = cls) -> int | None:
                v = (cls.get(rank) or {}).get("id")
                return int(v) if v is not None else None
            tid = int(tax.get("tax_id") or (rec.get("query") or [0])[0])
            table[tid] = {
                "name": str((tax.get("current_scientific_name") or {}).get("name") or ""),
                "rank": str(tax.get("rank") or "").lower(),
                "kingdom": _name("kingdom"), "phylum": _name("phylum"), "class": _name("class"),
                "order": _name("order"), "family": _name("family"), "genus": _name("genus"),
                "species": _name("species"), "species_taxid": _id("species"),
                "genetic_code": ((tax.get("genetic_code") or {}).get("primary") or {}).get("id"),
                "mito_code": ((tax.get("genetic_code") or {}).get("mitochondrial") or {}).get("id"),
                "parents": [int(p) for p in (tax.get("parents") or [])],
            }
        for t in missing:
            table.setdefault(t, {"name": "", "rank": "", "kingdom": "", "genus": "", "family": "", "species": "",
                                 "species_taxid": None, "genetic_code": None, "mito_code": None, "parents": []})
        cache_path.write_text(json.dumps({str(k): v for k, v in table.items()}, indent=0))
    return table


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #


def build_lock(base: Path | None = None, *, sketch: bool = True, threads: int = 32, log=print) -> Lock:
    base = base or root()
    gdir = base / "genomes" / "ncbi_dataset" / "data"
    report = gdir / "assembly_data_report.jsonl"
    if not report.is_file():
        raise FileNotFoundError(f"no assembly report at {report}; run the genome fetch first")
    hash_cache_path = base / "hash_cache.json"
    hash_cache: dict[str, Any] = json.loads(hash_cache_path.read_text()) if hash_cache_path.is_file() else {}

    cgf = _accessions(base / "manifests" / "cgf_summary.jsonl")
    refseq = _accessions(base / "manifests" / "refseq_fungi_reference.jsonl")
    set((base / "manifests" / "panel_accessions.txt").read_text().split())

    records: dict[str, dict[str, Any]] = {}
    with report.open() as fh:
        for line in fh:
            rec = json.loads(line)
            records[rec["accession"]] = rec
    taxids = {int(r["organism"]["taxId"]) for r in records.values()}
    tax = taxonomy_table(taxids, base / "taxonomy_cache.json")

    assemblies: list[Assembly] = []
    seen_hash: dict[str, str] = {}
    n = 0
    for acc, rec in sorted(records.items()):
        d = gdir / acc
        fnas = sorted(d.glob("*_genomic.fna")) + sorted(d.glob("*_genomic.fna.gz"))
        if not fnas:
            continue  # not rehydrated yet; the lock records what is present
        fna = fnas[0]
        n += 1
        if n % 100 == 0:
            log(f"  hashed {n} assemblies")
        digest = _sha256(fna, hash_cache)
        bases, contigs = _fasta_stats(fna)
        tid = int(rec["organism"]["taxId"])
        t = tax.get(tid, {})
        src = "cgf" if acc in cgf else ("refseq_reference" if acc in refseq else "panel")
        if acc in cgf and acc in refseq:
            src = "both"
        parents = set(t.get("parents") or [])
        is_fungus = FUNGI_TAXID in parents or tid == FUNGI_TAXID
        dup = seen_hash.get(digest)
        if dup is None:
            seen_hash[digest] = acc
        assemblies.append(Assembly(
            accession=acc, organism=str(rec["organism"].get("organismName") or t.get("name") or ""), taxid=tid,
            species_taxid=t.get("species_taxid") or (tid if t.get("rank") == "species" else None),
            species=str(t.get("species") or rec["organism"].get("organismName") or ""), genus=str(t.get("genus") or ""),
            family=str(t.get("family") or ""), kingdom=str(t.get("kingdom") or ""),
            level=str((rec.get("assemblyInfo") or {}).get("assemblyLevel") or ""), source=src,
            path=str(fna.relative_to(base)), bases=bases, n_contigs=contigs, sha256=digest,
            genetic_code=t.get("genetic_code"), mito_code=t.get("mito_code"),
            strain=str(((rec["organism"].get("infraspecificNames") or {}).get("strain")) or ""),
            duplicate_of=dup, is_fungus=is_fungus,
        ))
    hash_cache_path.write_text(json.dumps(hash_cache))

    # Decoys: gut bacteria fetched for this lane, plus PhiX from the host set.
    decoys: list[dict[str, Any]] = []
    ddir = base / "decoys" / "ncbi_dataset" / "data"
    if (ddir / "assembly_data_report.jsonl").is_file():
        for line in (ddir / "assembly_data_report.jsonl").read_text().splitlines():
            rec = json.loads(line)
            fnas = sorted((ddir / rec["accession"]).glob("*_genomic.fna"))
            if fnas:
                decoys.append({"kind": "bacterium", "accession": rec["accession"],
                               "organism": rec["organism"]["organismName"], "taxid": int(rec["organism"]["taxId"]),
                               "path": str(fnas[0].relative_to(base)), "sha256": _sha256(fnas[0], hash_cache)})
    phix = base.parent / "host" / "phix.fna.gz"
    if phix.is_file():
        decoys.append({"kind": "control_virus", "accession": "PhiX174", "organism": "Escherichia phage phiX174",
                       "taxid": 10847, "path": str(Path("..") / "host" / "phix.fna.gz"), "sha256": _sha256(phix, hash_cache)})
    hash_cache_path.write_text(json.dumps(hash_cache))

    # Deposits, from the fetch manifest.
    deposits: list[dict[str, Any]] = []
    man = base / "deposits.manifest.tsv"
    if man.is_file():
        for row in man.read_text().splitlines()[1:]:
            dep, fname, url, size, sha, lic = row.split("\t")
            deposits.append({"deposit": dep, "file": fname, "url": url, "bytes": int(size), "sha256": sha, "license": lic})

    fungi = [a for a in assemblies if a.is_fungus and a.duplicate_of is None]
    counts = {
        "assemblies_present": len(assemblies),
        "assemblies_admitted_fungal": len(fungi),
        "assemblies_exact_duplicates": sum(1 for a in assemblies if a.duplicate_of),
        "assemblies_nonfungal_excluded": sum(1 for a in assemblies if not a.is_fungus),
        "species_admitted": len({a.species_taxid for a in fungi if a.species_taxid}),
        "genera_admitted": len({a.genus for a in fungi if a.genus}),
        "cgf_assemblies": sum(1 for a in fungi if a.source in ("cgf", "both")),
        "refseq_reference_assemblies": sum(1 for a in fungi if a.source in ("refseq_reference", "both")),
        "panel_assemblies": sum(1 for a in fungi if a.source == "panel"),
        "decoy_genomes": len(decoys),
        "admitted_bases": sum(a.bases for a in fungi),
    }
    # The id covers what the analysis aligns against - assemblies and
    # decoys - so recording a deposit (FungiGutDB is a baseline asset, not an
    # alignment target) does not invalidate every sample's cached work.
    blob = json.dumps({
        "v": LOCK_VERSION, "a": sorted((a.accession, a.sha256) for a in assemblies),
        "d": sorted((d["accession"], d["sha256"]) for d in decoys), "c": SKETCH_C,
    }, sort_keys=True).encode()
    lock = Lock(
        lock_id="myco-ref-" + hashlib.sha256(blob).hexdigest()[:16],
        built_at=time.strftime("%Y-%m-%dT%H:%M:%S"), assemblies=assemblies, decoys=decoys, deposits=deposits,
        taxonomy={t: tax[t] for t in taxids if t in tax}, counts=counts,
    )
    if sketch:
        lock.sketch_db = build_sketch(base, lock, threads=threads, log=log)
    elif (base / "sylph" / f"fungi-c{SKETCH_C}.syldb").is_file():
        lock.sketch_db = f"sylph/fungi-c{SKETCH_C}.syldb"
    save_lock(lock, base)
    return lock


def build_sketch(base: Path, lock: Lock, *, threads: int = 32, log=print) -> str | None:
    """Sketch every admitted fungal assembly into one sylph database at the
    same ``c`` as the read sketches, so a sample's existing read sketch can
    be profiled against fungi without re-reading the FASTQs."""
    sylph = _sylph_bin()
    if sylph is None:
        log("  sylph not available; no whole-genome fungal sketch built")
        return None
    out_dir = base / "sylph"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / f"fungi-c{SKETCH_C}"
    db = prefix.with_suffix(".syldb")
    stamp = out_dir / "sketch.lock.json"
    fungi = [a for a in lock.assemblies if a.is_fungus and a.duplicate_of is None]
    want = {"lock_id": lock.lock_id, "n": len(fungi), "c": SKETCH_C}
    if db.is_file() and stamp.is_file() and json.loads(stamp.read_text()) == want:
        return str(db.relative_to(base))
    listing = out_dir / "genomes.txt"
    listing.write_text("\n".join(str(base / a.path) for a in fungi) + "\n")
    log(f"  sketching {len(fungi)} fungal assemblies at c={SKETCH_C}")
    subprocess.run([str(sylph), "sketch", "-t", str(threads), "-c", str(SKETCH_C), "--gl", str(listing),
                    "-o", str(prefix)], check=True, capture_output=True, text=True)
    stamp.write_text(json.dumps(want))
    return str(db.relative_to(base))


def representative(assemblies: list[Assembly]) -> Assembly:
    """One assembly per species for discovery: a RefSeq reference first,
    then the most complete, then the longest."""
    order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    return sorted(assemblies, key=lambda a: (0 if a.source in ("refseq_reference", "both") else 1,
                                             order.get(a.level, 4), -a.bases))[0]


def build_species_index(lock: Lock, base: Path | None = None, *, threads: int = 32, log=print) -> Path | None:
    """One representative assembly per admitted species plus the decoys, as a
    single minimap2 short-read index for candidate discovery. Every read pair
    is aligned to it, so a fungus is found wherever it sits in the catalogue,
    not only if it is on a fixed panel. Built once; keyed on the lock."""
    base = base or root()
    mm = shutil.which("minimap2")
    if mm is None:
        log("  minimap2 not available; no species index built")
        return None
    out_dir = base / "species_index"
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "fungi_species.fa"
    mmi = out_dir / "fungi_species.mmi"
    meta = out_dir / "fungi_species.json"
    stamp = out_dir / "index.lock.json"
    reps = [representative(asms) for asms in lock.species_assemblies().values()]
    want = {"lock_id": lock.lock_id, "n": len(reps), "decoys": len(lock.decoys)}
    if mmi.is_file() and meta.is_file() and stamp.is_file() and json.loads(stamp.read_text()) == want:
        return mmi
    log(f"  writing species index FASTA: {len(reps)} representatives + {len(lock.decoys)} decoys")
    contig_acc: dict[str, str] = {}
    contig_len: dict[str, int] = {}
    with fasta.open("w") as out:
        for a in sorted(reps, key=lambda a: a.accession):
            _append_fasta(base / a.path, a.accession, out, contig_acc, contig_len)
        for d in lock.decoys:
            _append_fasta((base / d["path"]).resolve(), d["accession"], out, contig_acc, contig_len)
    meta.write_text(json.dumps({"contig_acc": contig_acc, "contig_len": contig_len,
                                "representatives": [a.accession for a in reps]}))
    log("  minimap2 -x sr index over the species FASTA (single part)")
    subprocess.run([mm, "-x", "sr", "-I", "64G", "-t", str(threads), "-d", str(mmi), str(fasta)],
                   check=True, capture_output=True, text=True)
    stamp.write_text(json.dumps(want))
    return mmi


def _append_fasta(src: Path, acc: str, out, contig_acc: dict[str, str], contig_len: dict[str, int]) -> None:
    opener = gzip.open if src.suffix == ".gz" else open
    name = None
    n = 0
    with opener(src, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contig_len[name] = n
                name = f"{acc}|{line[1:].split()[0]}"
                contig_acc[name] = acc
                n = 0
                out.write(f">{name}\n")
            else:
                n += len(line.strip())
                out.write(line)
        if name is not None:
            contig_len[name] = n


def save_lock(lock: Lock, base: Path | None = None) -> Path:
    base = base or root()
    out = base / "reference_lock.json"
    out.write_text(json.dumps(lock.to_json(), indent=0))
    return out


def load_lock(base: Path | None = None) -> Lock | None:
    base = base or root()
    path = base / "reference_lock.json"
    if not path.is_file():
        return None
    d = json.loads(path.read_text())
    if d.get("lock_version") != LOCK_VERSION:
        return None
    assemblies = [Assembly(**dict(a.items())) for a in d["assemblies"]]
    return Lock(
        lock_id=d["lock_id"], built_at=d["built_at"], assemblies=assemblies, decoys=d.get("decoys") or [],
        deposits=d.get("deposits") or [], taxonomy={int(k): v for k, v in (d.get("taxonomy") or {}).items()},
        counts=d.get("counts") or {}, sketch_db=d.get("sketch_db"),
        licenses=d.get("licenses") or dict(LICENSES), upstream_commits=d.get("upstream_commits") or dict(UPSTREAM_COMMITS),
    )


def available(base: Path | None = None) -> bool:
    base = base or root()
    lock = load_lock(base)
    return lock is not None and bool(lock.assemblies)


__all__ = ["Assembly", "Lock", "PRESET", "SKETCH_C", "available", "build_lock", "build_sketch", "load_lock", "root", "save_lock",
           "taxonomy_table"]
