"""R5-R10 supplements against GTDB R232 and GlobDB r232 (spec 0.8.4 §3).

For every supplement, each species representative is placed in exactly one
category, with the evidence that put it there:

    already_represented       the genome, or its species, is a GTDB R232 /
                              GlobDB unit - matched by accession where the
                              source carries one, else by a genome-backed
                              walk (representative accession -> R232 species),
                              else by ANI+AF against representatives
    new_cluster               a GlobDB non-GTDB representative from this
                              source (GlobDB's clustering already admitted it)
    additional_strain_reference
                              a member genome of a known species that this
                              source contributes (kept for discrimination)
    better_reference          a higher-quality genome for a known species than
                              the one in hand
    unresolved_mapping        no accession, no genome-backed walk, no ANI run
                              yet - it will be compared when its genome is
                              fetched; never counted as new
    excluded_quality          fails the quality policy (completeness < 50 or
                              contamination > 10, or the source's own flag)

Counts are per source and never summed across sources: GlobDB already
consolidates several of them, and "146,310 non-GTDB representatives" is
not "146,310 gut species".

Outputs under refs/expanded/reconciliation/: one TSV per source, one
summary JSON with counts and the method used per row, and the rescue-panel
manifest (`rescue_panel.tsv`) naming every genome the gut panel should
carry with its source, category and where to fetch it.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import re
import tarfile
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from openbiota.expansion import ani, globdb, gtdb

OUT_DIR: Final = Path("refs/expanded/reconciliation")
SUPP: Final = Path("refs/supplements")
CATEGORIES: Final = ("already_represented", "new_cluster", "additional_strain_reference", "better_reference",
                     "unresolved_mapping", "excluded_quality")
MIN_COMPLETENESS: Final = 50.0
MAX_CONTAMINATION: Final = 10.0

_ACC = re.compile(r"(GC[AF]_\d+(?:\.\d+)?)")


@dataclass(frozen=True)
class Row:
    source: str
    source_id: str
    source_species_name: str      # the species name the source gives (its own release)
    source_taxonomy_release: str
    category: str
    method: str                    # accession | assembly_walk | globdb_rep | ani | name_match | none
    canonical_species: str         # GTDB R232 species when represented
    globdb_id: str                 # GlobDB representative when known
    ani: float | None
    af_shorter: float | None
    completeness: float | None
    contamination: float | None
    local_path: str                # genome on disk, when we have it
    fetch_hint: str                # where to get it otherwise


def _num(v: object) -> float | None:
    try:
        return float(str(v)) if str(v).strip() not in ("", "NA", "None", "nan") else None
    except ValueError:
        return None


def _quality_ok(comp: float | None, cont: float | None) -> bool:
    if comp is not None and comp < MIN_COMPLETENESS:
        return False
    return not (cont is not None and cont > MAX_CONTAMINATION)


def _write(source: str, rows: list[Row]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{source}.tsv.gz"
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(list(Row.__dataclass_fields__))
        for r in rows:
            w.writerow([getattr(r, f) for f in Row.__dataclass_fields__])
    return path


def _summary(source: str, rows: list[Row], extra: dict | None = None) -> dict:
    counts = Counter(r.category for r in rows)
    methods = Counter(r.method for r in rows)
    return {"source": source, "n_representatives": len(rows), "categories": {c: counts.get(c, 0) for c in CATEGORIES},
            "methods": dict(methods), **(extra or {})}


# --------------------------------------------------------------------------- #
# shared lookups
# --------------------------------------------------------------------------- #

class Backbone:
    """Everything a supplement is compared against, loaded once."""

    def __init__(self) -> None:
        self.r232 = gtdb.genome_table("r232")
        self.r232_species = {g.species for g in self.r232.values() if g.species}
        self.r232_by_unversioned = {gtdb.unversioned(a): g for a, g in self.r232.items()}
        self.reps = globdb.representatives()
        self.globdb_quality = globdb.quality()
        # UHGG v1 assembly names inside GTDB: MGYG-HGUT-00653 -> GCA accession
        self.uhgg_v1_in_gtdb: dict[str, gtdb.Genome] = {}
        with gzip.open(gtdb.GTDB_DIR / "r232" / "bac120_metadata_r232.tsv.gz", "rt") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                name = row.get("ncbi_assembly_name") or ""
                m = re.search(r"MGYG-HGUT-(\d{5})", name)
                if m:
                    g = self.r232.get(gtdb.bare(row["accession"]))
                    if g is not None:
                        self.uhgg_v1_in_gtdb[m.group(1)] = g
        self._clusters: dict[str, dict[str, gtdb.SpeciesCluster]] = {}

    def clusters(self, release: str) -> dict[str, gtdb.SpeciesCluster]:
        if release not in self._clusters:
            try:
                self._clusters[release] = gtdb.species_clusters(release)
            except FileNotFoundError:
                self._clusters[release] = {}
        return self._clusters[release]

    def walk_species(self, release: str, species: str) -> str:
        """A species name in an older release -> its R232 species, by the genomes."""
        cluster = self.clusters(release).get(species)
        if cluster is None:
            return ""
        counts: Counter[str] = Counter()
        for acc in cluster.members[:2000]:
            g = self.r232.get(acc) or self.r232_by_unversioned.get(gtdb.unversioned(acc))
            if g is not None and g.species:
                counts[g.species] += 1
        if not counts:
            return ""
        top, n = counts.most_common(1)[0]
        return top if n >= 0.5 * sum(counts.values()) else ""

    def accession_species(self, accession: str) -> str:
        g = self.r232.get(accession) or self.r232_by_unversioned.get(gtdb.unversioned(accession))
        return g.species if g else ""


# --------------------------------------------------------------------------- #
# R6 UHGG v2.0.2
# --------------------------------------------------------------------------- #

def reconcile_uhgg(bb: Backbone) -> tuple[Path, dict]:
    path = SUPP / "uhgg_v2.0.2" / "genomes-all_metadata.tsv"
    rows: list[Row] = []
    species_seen: set[str] = set()
    n_genomes = 0
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_genomes += 1
            if row["Genome"] != row["Species_rep"]:
                continue
            gid = row["Genome"]
            species_seen.add(gid)
            comp = _num(row.get("Completeness"))
            cont = _num(row.get("Contamination"))
            name = gtdb.species_of(row["Lineage"])
            fetch = f"https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/species_catalogue/{gid[:-2]}/{gid}/genome/{gid}.fna"
            if not _quality_ok(comp, cont):
                rows.append(Row("uhgg_v2.0.2", gid, name, "GTDB r95-era (MGnify)", "excluded_quality", "none", "", "", None, None, comp, cont, "", fetch))
                continue
            rep = bb.reps.get(gid)
            if rep is not None:
                rows.append(Row("uhgg_v2.0.2", gid, name, "GTDB r95-era (MGnify)", "new_cluster", "globdb_rep", "", gid, None, None, comp, cont, "", fetch))
                continue
            digits = re.sub(r"\D", "", gid.split(".", 1)[0].replace("MGYG", ""))
            num = int(digits) if digits else 0
            g = bb.uhgg_v1_in_gtdb.get(str(num).zfill(5)) if 0 < num <= 4644 else None
            if g is not None and g.species:
                rows.append(Row("uhgg_v2.0.2", gid, name, "GTDB r95-era (MGnify)", "already_represented", "accession",
                                g.species, g.accession, None, None, comp, cont, "", fetch))
                continue
            if name and name in bb.r232_species:
                rows.append(Row("uhgg_v2.0.2", gid, name, "GTDB r95-era (MGnify)", "already_represented", "name_match",
                                name, "", None, None, comp, cont, "", fetch))
                continue
            rows.append(Row("uhgg_v2.0.2", gid, name, "GTDB r95-era (MGnify)", "unresolved_mapping", "none", "", "", None, None, comp, cont, "", fetch))
    out = _write("uhgg_v2.0.2", rows)
    return out, _summary("uhgg_v2.0.2", rows, {"n_genomes_in_catalogue": n_genomes, "n_species": len(species_seen)})


# --------------------------------------------------------------------------- #
# R5 HRGM2
# --------------------------------------------------------------------------- #

def _hrgm2_rep_paths() -> dict[str, Path]:
    root = SUPP / "hrgm2" / "HRGMv2_Rep_Genome"
    tar = SUPP / "hrgm2" / "HRGMv2_Rep_Genome.tar.gz"
    if not root.is_dir() and tar.is_file():
        with tarfile.open(tar) as tf:
            tf.extractall(SUPP / "hrgm2")  # noqa: S202 - Zenodo archive, size-verified
    return {p.stem: p for p in root.rglob("*.fna")} if root.is_dir() else {}


def reconcile_hrgm2(bb: Backbone) -> tuple[Path, dict]:
    meta = SUPP / "hrgm2" / "HRGMv2_Cluster_metadata.tsv"
    r220 = SUPP / "hrgm2" / "HRGMv2_gtdbr220_results.tsv"
    names220: dict[str, str] = {}
    with r220.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            names220[row["HRGMv2"]] = gtdb.species_of(row.get("GTDBr220_results", ""))
    paths = _hrgm2_rep_paths()
    rows: list[Row] = []
    with meta.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = row["HRGMv2 Cluster"]
            sp220 = names220.get(cid, "")
            gid = f"HRGMV2_{cid.split('_', 1)[1]}"
            local = str(paths.get(cid, ""))
            fetch = "Zenodo 19482781 (HRGMv2_Rep_Genome.tar.gz)"
            if gid in bb.reps:
                rows.append(Row("hrgm2", cid, sp220, "GTDB r220", "new_cluster", "globdb_rep", "", gid, None, None, None, None, local, fetch))
                continue
            walked = bb.walk_species("r220", sp220) if sp220 else ""
            if walked:
                rows.append(Row("hrgm2", cid, sp220, "GTDB r220", "already_represented", "assembly_walk", walked, "", None, None, None, None, local, fetch))
                continue
            if sp220 and sp220 in bb.r232_species:
                rows.append(Row("hrgm2", cid, sp220, "GTDB r220", "already_represented", "name_match", sp220, "", None, None, None, None, local, fetch))
                continue
            rows.append(Row("hrgm2", cid, sp220, "GTDB r220", "unresolved_mapping", "none", "", "", None, None, None, None, local, fetch))
    out = _write("hrgm2", rows)
    return out, _summary("hrgm2", rows, {"n_rep_genomes_on_disk": len(paths)})


# --------------------------------------------------------------------------- #
# R7 ELGG  (accessions exact; MAGs by ANI against the gut representatives)
# --------------------------------------------------------------------------- #

def _elgg_paths() -> dict[str, Path]:
    root = SUPP / "elgg" / "ELGG_representatives_2172"
    z = SUPP / "elgg" / "ELGG_representatives_2172.zip"
    if not root.is_dir() and z.is_file():
        with zipfile.ZipFile(z) as zf:
            zf.extractall(SUPP / "elgg")  # noqa: S202
    return {p.name: p for p in root.glob("*.fna*")} if root.is_dir() else {}


def reconcile_elgg(bb: Backbone, *, threads: int = 16, run_ani: bool = True) -> tuple[Path, dict]:
    paths = _elgg_paths()
    hrgm = _hrgm2_rep_paths()
    rows: list[Row] = []
    mags: list[Path] = []
    for name, p in sorted(paths.items()):
        m = _ACC.match(name)
        if m:
            sp = bb.accession_species(m.group(1))
            if sp:
                rows.append(Row("elgg", name, "", "GTDB (accession)", "already_represented", "accession", sp, "", None, None, None, None, str(p), "NCBI"))
            else:
                rows.append(Row("elgg", name, "", "", "unresolved_mapping", "none", "", "", None, None, None, None, str(p), "NCBI"))
        else:
            mags.append(p)
    ani_summary: dict = {"n_mags": len(mags), "ani_run": False}
    if run_ani and mags and hrgm and ani.available():
        pairs = ani.dist(mags, list(hrgm.values()), work_dir=OUT_DIR / "ani_elgg_vs_hrgm2", threads=threads)
        best = ani.best_by_query(pairs)
        # HRGM2 reps carry their own reconciliation; a MAG matching one is
        # represented by whatever that representative is.
        hrgm_rows = {r.source_id: r for r in _read("hrgm2")} if (OUT_DIR / "hrgm2.tsv.gz").is_file() else {}
        n_same = n_boundary = 0
        for p in mags:
            b = best.get(p.name)
            if b is not None and b.verdict == "same_species":
                n_same += 1
                cid = Path(b.reference).stem
                target = hrgm_rows.get(cid)
                canonical = target.canonical_species if target else ""
                cat = "additional_strain_reference" if (target and target.category in ("already_represented", "new_cluster")) else "unresolved_mapping"
                rows.append(Row("elgg", p.name, "", "", cat, "ani", canonical, target.globdb_id if target else "", b.ani, b.af_shorter, None, None, str(p), "Zenodo 6969520"))
            elif b is not None and b.verdict == "boundary":
                n_boundary += 1
                rows.append(Row("elgg", p.name, "", "", "unresolved_mapping", "ani", "", "", b.ani, b.af_shorter, None, None, str(p), "Zenodo 6969520"))
            else:
                rows.append(Row("elgg", p.name, "", "", "unresolved_mapping", "ani", "", "", b.ani if b else None, b.af_shorter if b else None, None, None, str(p), "Zenodo 6969520"))
        ani_summary.update({"ani_run": True, "n_same_species_as_hrgm2_rep": n_same, "n_boundary": n_boundary,
                            "rule": f"ANI>={ani.ANI_SAME_SPECIES} and AF(shorter)>={ani.AF_SAME_SPECIES}; skani {ani.skani_version()}"})
    else:
        for p in mags:
            rows.append(Row("elgg", p.name, "", "", "unresolved_mapping", "none", "", "", None, None, None, None, str(p), "Zenodo 6969520"))
    out = _write("elgg", rows)
    return out, _summary("elgg", rows, ani_summary)


# --------------------------------------------------------------------------- #
# R8 HumGut2  (genome level, accession where present)
# --------------------------------------------------------------------------- #

def reconcile_humgut2(bb: Backbone) -> tuple[Path, dict]:
    path = SUPP / "humgut2" / "HumGut2.tsv"
    rows: list[Row] = []
    clusters95: set[str] = set()
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            clusters95.add(row.get("cluster95", ""))
            name = row.get("gtdb_name", "") or gtdb.species_of(row.get("gtdb_taxonomy", ""))
            comp = _num(row.get("completeness"))
            cont = _num(row.get("contamination"))
            fname = (row.get(r.fieldnames[-1]) or "").strip()
            m = _ACC.search(fname) or _ACC.search(row.get("HumGut_name", ""))
            if not _quality_ok(comp, cont):
                rows.append(Row("humgut2", row["HumGut_name"], name, "GTDB (authors' release)", "excluded_quality", "none", "", "", None, None, comp, cont, "", fname))
                continue
            if m:
                sp = bb.accession_species(m.group(1))
                if sp:
                    cat = "already_represented"
                    rows.append(Row("humgut2", row["HumGut_name"], name, "GTDB (authors' release)", cat, "accession", sp, "", None, None, comp, cont, "", f"NCBI {m.group(1)}"))
                    continue
            walked = bb.walk_species("r220", name) or bb.walk_species("r207", name) if name else ""
            if walked:
                rows.append(Row("humgut2", row["HumGut_name"], name, "GTDB (authors' release)", "already_represented", "assembly_walk", walked, "", None, None, comp, cont, "", fname))
            elif name in bb.r232_species:
                rows.append(Row("humgut2", row["HumGut_name"], name, "GTDB (authors' release)", "already_represented", "name_match", name, "", None, None, comp, cont, "", fname))
            else:
                rows.append(Row("humgut2", row["HumGut_name"], name, "GTDB (authors' release)", "unresolved_mapping", "none", "", "", None, None, comp, cont, "", fname))
    out = _write("humgut2", rows)
    return out, _summary("humgut2", rows, {"n_clusters95": len(clusters95), "count_unit": "genome (not species)"})


# --------------------------------------------------------------------------- #
# R9 HROM
# --------------------------------------------------------------------------- #

def reconcile_hrom(bb: Backbone) -> tuple[Path, dict]:
    path = SUPP / "hrom" / "HROM-Species-metadata.tsv"
    rows: list[Row] = []
    shared = 0
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            sid = row["species"]
            name = gtdb.species_of(row.get("gtdb_taxonomy", ""))
            comp = _num(row.get("completeness"))
            cont = _num(row.get("contamination"))
            if "shared" in (row.get("gut-oral shared status") or "").lower():
                shared += 1
            fetch = f"https://www.decodebiome.org/HROM/data/genome_catalog/HROM_nonredundant_genomes/{sid}/{sid}_1.fna"
            if not _quality_ok(comp, cont):
                rows.append(Row("hrom", sid, name, "GTDB (authors' release)", "excluded_quality", "none", "", "", None, None, comp, cont, "", fetch))
                continue
            walked = bb.walk_species("r220", name) if name else ""
            if walked:
                rows.append(Row("hrom", sid, name, "GTDB (authors' release)", "already_represented", "assembly_walk", walked, "", None, None, comp, cont, "", fetch))
            elif name and name in bb.r232_species:
                rows.append(Row("hrom", sid, name, "GTDB (authors' release)", "already_represented", "name_match", name, "", None, None, comp, cont, "", fetch))
            else:
                rows.append(Row("hrom", sid, name, "GTDB (authors' release)", "unresolved_mapping", "none", "", "", None, None, comp, cont, "", fetch))
    out = _write("hrom", rows)
    return out, _summary("hrom", rows, {"n_gut_oral_shared": shared})


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def _read(source: str) -> list[Row]:
    path = OUT_DIR / f"{source}.tsv.gz"
    out: list[Row] = []
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for k in ("ani", "af_shorter", "completeness", "contamination"):
                row[k] = float(row[k]) if row[k] not in ("", "None") else None
            out.append(Row(**row))
    return out


def rescue_panel_manifest() -> Path:
    """Every genome the gut rescue panel should carry, and why.

    New clusters from gut sources, plus additional strain references, plus
    the food/background competitors GlobDB carries (cFMD). Host and PhiX
    come from the existing host filter reference.
    """
    rows: list[dict] = []
    for source in ("uhgg_v2.0.2", "hrgm2", "elgg", "humgut2", "hrom"):
        if not (OUT_DIR / f"{source}.tsv.gz").is_file():
            continue
        for r in _read(source):
            if r.category in ("new_cluster", "additional_strain_reference", "better_reference"):
                rows.append({"source": r.source, "source_id": r.source_id, "category": r.category,
                             "canonical_species": r.canonical_species, "globdb_id": r.globdb_id,
                             "local_path": r.local_path, "fetch_hint": r.fetch_hint, "role": "target"})
    for gid, rep in globdb.representatives().items():
        if rep.source == "cFMD":
            rows.append({"source": "globdb_cfmd", "source_id": gid, "category": "background_competitor",
                         "canonical_species": rep.species, "globdb_id": gid, "local_path": "", "fetch_hint": "GlobDB genome_fasta", "role": "competitor"})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "rescue_panel.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["source"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return path


def run_all(*, threads: int = 16, run_ani: bool = True) -> dict:
    bb = Backbone()
    summaries: dict[str, dict] = {}
    for fn in (reconcile_uhgg, reconcile_hrgm2, reconcile_humgut2, reconcile_hrom):
        _, s = fn(bb)
        summaries[s["source"]] = s
    _, s = reconcile_elgg(bb, threads=threads, run_ani=run_ani)
    summaries[s["source"]] = s
    panel = rescue_panel_manifest()
    result = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "backbone": {"gtdb": "R232", "globdb": "r232", "n_globdb_representatives": len(bb.reps),
                     "n_uhgg_v1_genomes_in_gtdb": len(bb.uhgg_v1_in_gtdb)},
        "quality_policy": {"min_completeness": MIN_COMPLETENESS, "max_contamination": MAX_CONTAMINATION},
        "clustering_rule": {"ani": ani.ANI_SAME_SPECIES, "af_shorter": ani.AF_SAME_SPECIES, "chaining": "none"},
        "sources": summaries,
        "rescue_panel": str(panel),
        "note": "counts are per source and are not summed; GlobDB already consolidates UHGG, HRGM2 and mOTUs",
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":  # pragma: no cover
    import sys
    print(json.dumps(run_all(threads=int(sys.argv[1]) if len(sys.argv) > 1 else 16), indent=1))
