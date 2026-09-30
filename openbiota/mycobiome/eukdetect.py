"""EukDetect2: the single-copy-marker lane (spec §4, §5.2, §6.3).

Runs the pinned EukDetect2 (commit 8d69014, database Zenodo 19056625) in
its own virtual environment and parses its native tables. The fields keep
their native meaning: RPKS is reads per kilobase of *marker* sequence, a
marker-length-normalised support figure, not a sequencing-depth-normalised
fungal load; ``Relative_abundance`` is relative to the eukaryotic universe
EukDetect saw, not to the sample. Non-fungal eukaryotes (protists,
helminths, plants, animals) are returned separately and never counted as
fungi.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

EUKDETECT_COMMIT: Final = "8d69014727b2c5956de30b0811644b6c3a7bd4f3"
DATABASE_DEPOSIT: Final = "zenodo.19056625"
CACHE_VERSION: Final = 2
FUNGI_TAXID: Final = 4751
#: EukDetect2's own filter: at least this many markers and reads per taxon.
NATIVE_MIN_MARKERS: Final = 2
NATIVE_MIN_READS: Final = 4

REQUIRED_DB_FILES: Final = (
    "eukdb.1.bt2l", "eukdb.2.bt2l", "eukdb.3.bt2l", "eukdb.4.bt2l", "eukdb.rev.1.bt2l", "eukdb.rev.2.bt2l",
    "eukdb.fasta", "taxa.sqlite", "taxa.sqlite.traverse.pkl", "busco_taxid_genome_link.txt",
    "specific_and_inherited_markers_per_taxid.txt", "taxid_and_genome_cumulativelength.txt", "all_genomes.txt",
)


@dataclass(frozen=True, slots=True)
class MarkerHit:
    """One taxon EukDetect2 reported after its native filter."""

    name: str
    taxid: int
    rank: str
    lineage: str
    total_reads: int
    marker_length_bp: int
    rpks: float
    reads_aligned: int
    pid_aligned: float | None
    genomes: str
    reads_reassigned: int
    observed_markers: int | None
    percent_identity: float | None
    relative_abundance_eukaryotic: float | None
    is_fungus: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name, "taxid": self.taxid, "rank": self.rank, "lineage": self.lineage,
            "total_reads": self.total_reads, "marker_length_bp": self.marker_length_bp,
            "rpks_reads_per_kb_marker": self.rpks, "reads_aligned": self.reads_aligned,
            "percent_identity_aligned": self.pid_aligned, "genomes": self.genomes,
            "reads_reassigned": self.reads_reassigned, "observed_markers": self.observed_markers,
            "percent_identity": self.percent_identity,
            "relative_abundance_within_eukdetect_eukaryotes": self.relative_abundance_eukaryotic,
            "is_fungus": self.is_fungus,
            "units_note": "RPKS is marker-length-normalised support; relative abundance is within EukDetect's "
                          "eukaryotic universe. Neither is a fraction of the sample.",
        }


@dataclass(slots=True)
class MarkerProfile:
    status: str                       # resolved | not_assessed | failed
    tool: str
    database: str
    hits: list[MarkerHit]
    nonfungal: list[MarkerHit]
    all_hits_table: str | None
    cached: bool
    error: str | None = None

    @property
    def fungal(self) -> list[MarkerHit]:
        return [h for h in self.hits if h.is_fungus]

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status, "tool": self.tool, "database": self.database, "cached": self.cached,
            "error": self.error, "native_filter": {"min_markers": NATIVE_MIN_MARKERS, "min_reads": NATIVE_MIN_READS},
            "n_fungal": len(self.fungal), "n_nonfungal_eukaryotes": len(self.nonfungal),
            "fungal": [h.to_json() for h in self.fungal], "nonfungal_eukaryotes": [h.to_json() for h in self.nonfungal],
        }


def _tools() -> Path | None:
    root = Path(__file__).resolve().parents[2]
    exe = root / ".venv-eukdetect" / "bin" / "eukdetect"
    if exe.is_file():
        return exe
    found = shutil.which("eukdetect")
    return Path(found) if found else None


def available(db_dir: Path) -> bool:
    return _tools() is not None and all((db_dir / f).is_file() for f in REQUIRED_DB_FILES)


def _stamp(p: Path | None) -> str:
    if p is None or not p.is_file():
        return ""
    st = p.stat()
    return f"{p.name}:{st.st_size}:{int(st.st_mtime)}"


def _fp(*parts: str) -> str:
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def _f(v: str) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _i(v: str) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


def _read_len(fastq: Path, n: int = 2000) -> int:
    import gzip
    opener = gzip.open if fastq.suffix == ".gz" else open
    lengths: list[int] = []
    with opener(fastq, "rt") as fh:  # type: ignore[operator]
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lengths.append(len(line.strip()))
                if len(lengths) >= n:
                    break
    lengths.sort()
    return lengths[len(lengths) // 2] if lengths else 150


def run_eukdetect(*, sample: str, r1: Path, r2: Path | None, db_dir: Path, work_dir: Path,
                  threads: int = 8) -> MarkerProfile | None:
    """Run EukDetect2 on the QC non-host reads, or return the cached parse.

    Returns None when the tool or database is absent: the lane is then
    ``not_assessed``, never "no fungi".
    """
    exe = _tools()
    if exe is None or not available(db_dir):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    key = _fp(str(CACHE_VERSION), EUKDETECT_COMMIT, DATABASE_DEPOSIT, _stamp(r1), _stamp(r2))
    cached = work_dir / f"eukdetect.{key}.json"
    if cached.is_file():
        d = json.loads(cached.read_text())
        return _from_json(d, cached=True)

    out = work_dir / "eukdetect_out"
    out.mkdir(parents=True, exist_ok=True)
    table = out / f"{sample}_filtered_hits_table.txt"
    eukfrac = out / f"{sample}_filtered_hits_eukfrac.txt"
    all_hits = out / "filtering" / f"{sample}_all_hits_table.txt"
    if table.is_file():
        # The tool already ran for these reads; only the parse changed.
        prof = parse(table, eukfrac if eukfrac.is_file() else None, all_hits if all_hits.is_file() else None)
        cached.write_text(json.dumps(prof.to_json()))
        return prof
    cmd = [str(exe), "single", "-1", str(r1)]
    if r2 is not None:
        cmd += ["-2", str(r2)]
    cmd += ["-n", sample, "--outdir", str(out), "--database", str(db_dir), "--cores", str(max(1, threads)),
            "--readlen", str(_read_len(r1))]
    import os
    env = dict(os.environ)
    env["PATH"] = f"{exe.parent}{os.pathsep}" + env.get("PATH", "")
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
    (work_dir / "eukdetect.log").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0 or not table.is_file():
        prof = MarkerProfile(status="failed", tool=f"EukDetect2 {EUKDETECT_COMMIT[:7]}", database=DATABASE_DEPOSIT,
                             hits=[], nonfungal=[], all_hits_table=None, cached=False,
                             error=(proc.stderr or proc.stdout)[-2000:])
        cached.write_text(json.dumps(prof.to_json()))
        return prof
    prof = parse(table, eukfrac if eukfrac.is_file() else None, all_hits if all_hits.is_file() else None)
    cached.write_text(json.dumps(prof.to_json()))
    return prof


#: EukDetect2 lineages carry no kingdom token ("phylum-Ascomycota|class-...");
#: fungi are recognised by phylum.
FUNGAL_PHYLA: Final = frozenset({
    "Ascomycota", "Basidiomycota", "Mucoromycota", "Zoopagomycota", "Chytridiomycota", "Blastocladiomycota",
    "Microsporidia", "Cryptomycota", "Rozellomycota", "Glomeromycota", "Mortierellomycota", "Kickxellomycota",
    "Entomophthoromycota", "Olpidiomycota", "Neocallimastigomycota", "Aphelidiomycota", "Monoblepharomycota",
    "Sanchytriomycota", "Basidiobolomycota", "Entorrhizomycota",
})


def _is_fungal_lineage(lineage: str) -> bool:
    if "Fungi" in lineage:
        return True
    for tok in lineage.replace(";", "|").split("|"):
        tok = tok.strip()
        if tok.lower().startswith("phylum-"):
            return tok.split("-", 1)[1] in FUNGAL_PHYLA
        if tok.startswith("p__"):
            return tok[3:] in FUNGAL_PHYLA
    return False


def parse(table: Path, eukfrac: Path | None, all_hits: Path | None) -> MarkerProfile:
    rel: dict[int, float] = {}
    if eukfrac is not None:
        rows = eukfrac.read_text().splitlines()
        if rows:
            hdr = rows[0].split("\t")
            for line in rows[1:]:
                parts = line.split("\t")
                if len(parts) < len(hdr):
                    continue
                rec = dict(zip(hdr, parts, strict=False))
                if rec.get("TaxID"):
                    rel[_i(rec["TaxID"])] = _f(rec.get("Relative_abundance", "")) or 0.0
    obs: dict[int, tuple[int | None, float | None]] = {}
    if all_hits is not None:
        rows = all_hits.read_text().splitlines()
        if rows:
            hdr = rows[0].split("\t")
            for line in rows[1:]:
                parts = line.split("\t")
                if len(parts) < len(hdr):
                    continue
                rec = dict(zip(hdr, parts, strict=False))
                obs[_i(rec.get("Taxid", "0"))] = (_i(rec.get("Observed_markers", "0")) or None,
                                                   _f(rec.get("Percent_identity", "")))
    hits: list[MarkerHit] = []
    non: list[MarkerHit] = []
    rows = table.read_text().splitlines()
    hdr = rows[0].split("\t") if rows else []
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        rec = dict(zip(hdr, parts, strict=False))
        tid = _i(rec.get("Taxid", "0"))
        lineage = rec.get("Lineage", "")
        is_fungus = _is_fungal_lineage(lineage)
        h = MarkerHit(
            name=rec.get("Name", ""), taxid=tid, rank=rec.get("Rank", ""), lineage=lineage,
            total_reads=_i(rec.get("Total_reads", "0")), marker_length_bp=_i(rec.get("Total_marker_length", "0")),
            rpks=_f(rec.get("RPKS", "")) or 0.0, reads_aligned=_i(rec.get("Reads_aligned", "0")),
            pid_aligned=_f(rec.get("PID_aligned", "")), genomes=rec.get("Genomes", ""),
            reads_reassigned=_i(rec.get("Reads_reassigned", "0")),
            observed_markers=obs.get(tid, (None, None))[0], percent_identity=obs.get(tid, (None, None))[1],
            relative_abundance_eukaryotic=rel.get(tid), is_fungus=is_fungus,
        )
        (hits if is_fungus else non).append(h)
    return MarkerProfile(status="resolved", tool=f"EukDetect2 {EUKDETECT_COMMIT[:7]}", database=DATABASE_DEPOSIT,
                         hits=hits, nonfungal=non, all_hits_table=str(all_hits) if all_hits else None, cached=False)


def _from_json(d: dict[str, Any], *, cached: bool) -> MarkerProfile:
    def hit(x: dict[str, Any]) -> MarkerHit:
        return MarkerHit(
            name=x["name"], taxid=int(x["taxid"]), rank=x["rank"], lineage=x["lineage"], total_reads=int(x["total_reads"]),
            marker_length_bp=int(x["marker_length_bp"]), rpks=float(x["rpks_reads_per_kb_marker"]),
            reads_aligned=int(x["reads_aligned"]), pid_aligned=x.get("percent_identity_aligned"), genomes=x.get("genomes", ""),
            reads_reassigned=int(x.get("reads_reassigned", 0)), observed_markers=x.get("observed_markers"),
            percent_identity=x.get("percent_identity"),
            relative_abundance_eukaryotic=x.get("relative_abundance_within_eukdetect_eukaryotes"),
            is_fungus=bool(x["is_fungus"]),
        )
    return MarkerProfile(status=d["status"], tool=d["tool"], database=d["database"],
                         hits=[hit(x) for x in d.get("fungal", [])], nonfungal=[hit(x) for x in d.get("nonfungal_eukaryotes", [])],
                         all_hits_table=None, cached=cached, error=d.get("error"))


__all__ = ["DATABASE_DEPOSIT", "EUKDETECT_COMMIT", "MarkerHit", "MarkerProfile", "available", "parse", "run_eukdetect"]
