"""Whole-genome profiling with Sylph against a GTDB species universe.

What this lane adds
-------------------
The marker lanes identify an organism from a few dozen marker genes drawn
from a catalogue of reference genomes. Sylph works differently: it sketches
the whole read set and asks, for each of ~200,000 species-representative
genomes, how much of that genome's k-mer content is contained in the sample,
correcting for coverage. Its universe is GTDB - every bacterial and archaeal
species cluster with a genome - and it resolves organisms the marker
catalogues have no markers for.

It is a detection lane. Its abundances are a coverage-normalised share of
the genomes it recognised, a different quantity with a different denominator
from the marker lanes, and no reference cohort exists on it. So its
detections join the inventory, its abundances are kept as secondary
readings, and nothing here is given a percentile.

What a hit means, and does not
------------------------------
Sylph reports a containment ANI for each genome: how similar the sample's
copy of that species is to the reference. That is a similarity to a genome,
not a strain identity, and it is only a reliable estimate when coverage is
high enough for the correction to work. Where the correction fails Sylph
says so (``Contig_name`` rows with ``LOW`` status and no interval), and
those rows are kept but never promoted to a confident species call. The
minimum ANI is 95%, GTDB's own species boundary, and is not lowered to find
more.

Every column Sylph writes is retained. The parser reads the profile through
``sylph-tax``, which attaches the GTDB lineage to each genome accession.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

#: Pinned tool and database. Changing either invalidates every cached profile.
SYLPH_VERSION: Final = "1.0.0"
DATABASE_FILE: Final = "gtdb-r232-c200-dbv2.syl2db"
TAXONOMY_ID: Final = "GTDB_r232"
#: Species clusters in the pinned GTDB release, for the technical record.
CATALOGUE_SIZE: Final = 199_923
#: GTDB's species boundary. The starting point, not a knob.
MINIMUM_ANI: Final = 95.0
#: Sketch subsampling rate; c200 is Sylph's recommended sensitivity setting.
SKETCH_C: Final = 200
CACHE_VERSION: Final = 1


@dataclass(frozen=True, slots=True)
class GenomeHit:
    """One species-representative genome Sylph found in the sample."""

    #: Full GTDB lineage, ``d__;p__;c__;o__;f__;g__;s__``.
    lineage: str
    #: GTDB species name, e.g. ``Phocaeicola vulgatus`` or ``Gemmiger sp937890665``.
    species: str
    genus: str
    family: str
    domain: str
    #: Coverage-normalised share of recognised genomes, in percent.
    taxonomic_abundance: float
    #: Share of sequence content, in percent.
    sequence_abundance: float
    #: Containment ANI to the representative genome, in percent.
    ani: float
    #: Effective coverage of the representative.
    coverage: float
    #: Whether the ANI correction succeeded. ``LOW`` means it did not and the
    #: ANI printed is not a confident estimate.
    correction: str
    #: GTDB genome accession of the representative.
    accession: str
    #: Every column Sylph wrote, untouched.
    raw: dict[str, str] = field(default_factory=dict)

    @property
    def placeholder(self) -> bool:
        """``Genus sp012345678``: a species cluster with no formal name."""
        parts = self.species.split(" ", 1)
        return len(parts) == 2 and parts[1].startswith("sp") and parts[1][2:3].isdigit()

    @property
    def confident(self) -> bool:
        """The coverage model fitted and the call clears the species boundary.

        ``LOW`` means Sylph could not fit its coverage correction - the
        printed ANI is a lower bound, not an estimate - so the row is kept
        but not treated as a confident species call. ``HIGH`` is the
        opposite regime (coverage saturated the model) and is fine.
        """
        return self.correction.upper() != "LOW" and self.ani >= MINIMUM_ANI

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "genus": self.genus,
            "family": self.family,
            "domain": self.domain,
            "lineage": self.lineage,
            "accession": self.accession,
            "taxonomic_abundance": round(self.taxonomic_abundance, 4),
            "sequence_abundance": round(self.sequence_abundance, 4),
            "ani": round(self.ani, 2),
            "coverage": round(self.coverage, 3),
            "correction": self.correction,
            "confident": self.confident,
            "placeholder": self.placeholder,
        }


@dataclass(frozen=True, slots=True)
class GenomeProfile:
    """One sample's whole-genome profile."""

    sample: str
    database: str
    tool_version: str
    taxonomy: str
    hits: tuple[GenomeHit, ...]
    elapsed_s: float
    cached: bool
    command: tuple[str, ...] = ()

    @property
    def confident_hits(self) -> tuple[GenomeHit, ...]:
        return tuple(h for h in self.hits if h.confident)

    def to_json(self) -> dict[str, Any]:
        conf = self.confident_hits
        return {
            "status": "resolved",
            "tool": f"sylph {self.tool_version}",
            "database": self.database,
            "taxonomy": self.taxonomy,
            "catalogue_size": CATALOGUE_SIZE,
            "minimum_ani": MINIMUM_ANI,
            "n_hits": len(self.hits),
            "n_confident": len(conf),
            "n_named": sum(1 for h in conf if not h.placeholder),
            "n_placeholder": sum(1 for h in conf if h.placeholder),
            "n_low_correction": sum(1 for h in self.hits if not h.confident),
            "elapsed_s": round(self.elapsed_s, 1),
            "cached": self.cached,
            "hits": [h.to_json() for h in self.hits],
            "what_this_is": (
                "Whole-genome containment profiling against every bacterial and archaeal "
                "species cluster in GTDB. Detections join the organism inventory; abundances "
                "are a coverage-normalised share of recognised genomes and are kept separate "
                "from the marker-lane composition. No reference cohort exists on this lane, "
                "so nothing here carries a percentile."
            ),
        }


def _tools() -> tuple[Path | None, Path | None]:
    root = Path(__file__).resolve().parents[2]
    sylph = root / "vendor" / "sylph" / "bin" / "sylph"
    if not sylph.is_file():
        found = shutil.which("sylph")
        sylph = Path(found) if found else None
    tax = root / ".venv" / "bin" / "sylph-tax"
    if not tax.is_file():
        found = shutil.which("sylph-tax")
        tax = Path(found) if found else None
    return sylph, tax


def available(db_dir: Path) -> bool:
    sylph, tax = _tools()
    return sylph is not None and tax is not None and (db_dir / DATABASE_FILE).is_file()


def _fingerprint(*parts: str) -> str:
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()[:16]


def run_sylph(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    db_dir: Path,
    work_dir: Path,
    threads: int = 8,
) -> GenomeProfile | None:
    """Profile one sample, caching by input content and pinned versions.

    Returns ``None`` when the tool or database is not installed; a missing
    lane is recorded as not assessed, never as zero organisms.
    """
    sylph, tax = _tools()
    db = db_dir / DATABASE_FILE
    if sylph is None or tax is None or not db.is_file():
        return None

    work_dir.mkdir(parents=True, exist_ok=True)
    key = _fingerprint(
        str(CACHE_VERSION), SYLPH_VERSION, DATABASE_FILE, TAXONOMY_ID, str(MINIMUM_ANI),
        str(SKETCH_C), _file_stamp(r1), _file_stamp(r2) if r2 else "",
    )
    cached = work_dir / f"sylph.{key}.json"
    if cached.is_file():
        data = json.loads(cached.read_text())
        return _profile_from_json(data, cached=True)

    t0 = time.monotonic()
    sketch_dir = work_dir / "sketches"
    sketch_dir.mkdir(exist_ok=True)
    sketch = sketch_dir / f"{sample}.paired.sylsp"
    if not sketch.is_file():
        cmd = [str(sylph), "sketch", "-1", str(r1)]
        if r2 is not None:
            cmd += ["-2", str(r2)]
        cmd += ["-S", sample, "-c", str(SKETCH_C), "-t", str(threads), "-d", str(sketch_dir)]
        subprocess.run(cmd, check=True, capture_output=True)

    profile_tsv = work_dir / f"{sample}.sylph.tsv"
    cmd = [
        str(sylph), "profile", "-d", str(db), str(sketch),
        "--minimum-ani", str(MINIMUM_ANI), "-u", "-t", str(threads), "-o", str(profile_tsv),
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    # sylph-tax writes <prefix><sample>.sylphmpa with the lineage per row.
    prefix = str(work_dir / f"{sample}.")
    subprocess.run(
        [str(tax), "taxprof", str(profile_tsv), "-t", TAXONOMY_ID, "-o", prefix],
        check=True, capture_output=True,
    )
    hits = _parse(profile_tsv)
    profile = GenomeProfile(
        sample=sample, database=DATABASE_FILE, tool_version=SYLPH_VERSION,
        taxonomy=TAXONOMY_ID, hits=tuple(hits), elapsed_s=time.monotonic() - t0,
        cached=False, command=tuple(cmd),
    )
    cached.write_text(json.dumps(profile.to_json(), indent=1))
    return profile


def _file_stamp(path: Path | None) -> str:
    if path is None:
        return ""
    st = path.stat()
    return f"{path.name}:{st.st_size}:{int(st.st_mtime)}"


def _parse(profile_tsv: Path) -> list[GenomeHit]:
    """Join Sylph's per-genome rows to the GTDB lineage of each genome.

    The raw profile's ``Contig_name`` carries the representative's GTDB
    accession; the lineage comes from the taxonomy table sylph-tax
    installed, 199,923 rows of ``accession<TAB>lineage``.
    """
    lineage_by_acc = _lineages()
    hits: list[GenomeHit] = []
    with profile_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = _accession(row.get("Genome_file") or row.get("Contig_name") or "")
            lineage = lineage_by_acc.get(acc, "")
            ranks = {seg[:3]: seg[3:] for seg in lineage.split(";") if "__" in seg}
            species = ranks.get("s__", "").strip()
            if not species:
                continue
            hits.append(GenomeHit(
                lineage=lineage,
                species=species,
                genus=ranks.get("g__", "").strip(),
                family=ranks.get("f__", "").strip(),
                domain=ranks.get("d__", "").strip(),
                taxonomic_abundance=_f(row.get("Taxonomic_abundance")),
                sequence_abundance=_f(row.get("Sequence_abundance")),
                ani=_f(row.get("Adjusted_ANI")),
                coverage=_f(row.get("True_cov")),
                # Eff_lambda is numeric when the coverage model fitted, LOW
                # when it could not (too little signal: the ANI is then not
                # a confident estimate), HIGH when coverage saturated it.
                correction=str(row.get("Eff_lambda") or "").strip() or "n/a",
                accession=acc,
                raw=dict(row),
            ))
    hits.sort(key=lambda h: -h.taxonomic_abundance)
    return hits


def _accession(genome_file: str) -> str:
    """GTDB genome accession from Sylph's genome path.

    ``.../GCF/964/248/265/GCF_964248265.1_genomic.fna.gz`` -> ``GCF_964248265.1``.
    """
    import re

    m = re.search(r"(GC[AF]_\d+\.\d+)", genome_file or "")
    if m:
        return m.group(1)
    head = (genome_file.split()[0] if genome_file else "").rsplit("/", 1)[-1]
    for suffix in (".fna.gz", ".fa.gz", ".fna", ".fa", "_genomic"):
        if head.endswith(suffix):
            head = head[: -len(suffix)]
    return head


def _lineages() -> dict[str, str]:
    """Accession -> lineage, from the taxonomy metadata sylph-tax installed."""
    out: dict[str, str] = {}
    cfg = Path.home() / ".config" / "sylph-tax" / "config.json"
    tax_dir: Path | None = None
    if cfg.is_file():
        try:
            tax_dir = Path(json.loads(cfg.read_text()).get("taxonomy_dir", ""))
        except (json.JSONDecodeError, OSError, TypeError):
            tax_dir = None
    if tax_dir is None or not tax_dir.is_dir():
        tax_dir = Path("refs/sylph/tax")
    for path in sorted(tax_dir.glob(f"*{TAXONOMY_ID.lower()}*")) + sorted(tax_dir.glob(f"*{TAXONOMY_ID}*")):
        opener = _opener(path)
        try:
            with opener(path) as fh:
                for line in fh:
                    if isinstance(line, bytes):
                        line = line.decode("utf-8", "replace")
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and "d__" in parts[1]:
                        out[_accession(parts[0])] = parts[1]
        except OSError:
            continue
        if out:
            break
    return out


def _opener(path: Path):
    if path.suffix == ".gz":
        import gzip
        return gzip.open
    if path.suffix == ".bz2":
        import bz2
        return bz2.open
    return open


def _f(v: Any) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _profile_from_json(data: dict[str, Any], *, cached: bool) -> GenomeProfile:
    hits = tuple(
        GenomeHit(
            lineage=str(h.get("lineage") or ""),
            species=str(h.get("species") or ""),
            genus=str(h.get("genus") or ""),
            family=str(h.get("family") or ""),
            domain=str(h.get("domain") or ""),
            taxonomic_abundance=_f(h.get("taxonomic_abundance")),
            sequence_abundance=_f(h.get("sequence_abundance")),
            ani=_f(h.get("ani")),
            coverage=_f(h.get("coverage")),
            correction=str(h.get("correction") or "n/a"),
            accession=str(h.get("accession") or ""),
        )
        for h in data.get("hits") or []
    )
    return GenomeProfile(
        sample="", database=str(data.get("database") or DATABASE_FILE),
        tool_version=str(data.get("tool") or "").replace("sylph ", "") or SYLPH_VERSION,
        taxonomy=str(data.get("taxonomy") or TAXONOMY_ID), hits=hits,
        elapsed_s=_f(data.get("elapsed_s")), cached=cached,
    )


__all__ = [
    "CATALOGUE_SIZE",
    "DATABASE_FILE",
    "GenomeHit",
    "GenomeProfile",
    "MINIMUM_ANI",
    "available",
    "run_sylph",
]
