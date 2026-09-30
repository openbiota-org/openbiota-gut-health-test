"""GTDB release backbone: genomes, species clusters, circumscription radii.

The R232 bacterial metadata is 900k rows and 270 MB compressed. Sample
runs and crosswalk builds need six of its 110 columns, so the first call
writes a slim table beside it (`genomes.slim.tsv.gz`) and everything after
reads that. The slim table is keyed by the bare accession (`GCF_000006945.2`)
with the `RS_`/`GB_` prefix stripped, because every other source in this
project uses the bare form.
"""

from __future__ import annotations

import csv
import functools
import gzip
import json
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

#: The E. coli cluster alone lists 44,640 member accessions in one field.
csv.field_size_limit(1 << 30)

GTDB_DIR: Final = Path("refs/gtdb")
SLIM_COLUMNS: Final = (
    "accession", "gtdb_taxonomy", "gtdb_representative", "gtdb_genome_representative",
    "checkm2_completeness", "checkm2_contamination", "ncbi_genbank_assembly_accession", "genome_size",
)


def bare(accession: str) -> str:
    """`RS_GCF_000006945.2` / `GB_GCA_...` -> `GCF_000006945.2`."""
    acc = accession.strip()
    if acc[:3] in ("RS_", "GB_"):
        acc = acc[3:]
    return acc


def unversioned(accession: str) -> str:
    return bare(accession).split(".", 1)[0]


def species_of(taxonomy: str) -> str:
    """`s__Escherichia coli` from a full GTDB string, or '' when unassigned."""
    for part in taxonomy.split(";"):
        part = part.strip()
        if part.startswith("s__"):
            return part[3:].strip()
    return ""


def genus_of(taxonomy: str) -> str:
    for part in taxonomy.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            return part[3:].strip()
    return ""


@dataclass(frozen=True)
class Genome:
    accession: str
    taxonomy: str
    is_representative: bool
    representative: str
    completeness: float | None
    contamination: float | None
    genbank_accession: str
    genome_size: int | None

    @property
    def species(self) -> str:
        return species_of(self.taxonomy)


def _slim_path(release: str) -> Path:
    return GTDB_DIR / release / "genomes.slim.tsv.gz"


def _metadata_files(release: str) -> list[Path]:
    d = GTDB_DIR / release
    return sorted(p for p in d.glob(f"*_metadata_{release}.tsv.gz"))


def build_slim(release: str = "r232") -> Path:
    """Write the six-column genome table once; idempotent."""
    out = _slim_path(release)
    if out.is_file():
        return out
    files = _metadata_files(release)
    if not files:
        raise FileNotFoundError(f"no GTDB {release} metadata under {GTDB_DIR / release}; run scripts/fetch_expanded_refs.py")
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SLIM_COLUMNS)
        for path in files:
            with gzip.open(path, "rt") as src:
                r = csv.DictReader(src, delimiter="\t")
                for row in r:
                    w.writerow([
                        bare(row["accession"]), row.get("gtdb_taxonomy", ""),
                        "t" if row.get("gtdb_representative") == "t" else "f",
                        bare(row.get("gtdb_genome_representative", "")),
                        row.get("checkm2_completeness", ""), row.get("checkm2_contamination", ""),
                        row.get("ncbi_genbank_assembly_accession", ""), row.get("genome_size", ""),
                    ])
    tmp.replace(out)
    return out


def iter_genomes(release: str = "r232") -> Iterator[Genome]:
    path = build_slim(release)
    with gzip.open(path, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            def _f(v: str) -> float | None:
                try:
                    return float(v) if v not in ("", "none", "None", "NA") else None
                except ValueError:
                    return None
            yield Genome(
                accession=row["accession"], taxonomy=row["gtdb_taxonomy"],
                is_representative=row["gtdb_representative"] == "t",
                representative=row["gtdb_genome_representative"],
                completeness=_f(row["checkm2_completeness"]), contamination=_f(row["checkm2_contamination"]),
                genbank_accession=row["ncbi_genbank_assembly_accession"],
                genome_size=int(row["genome_size"]) if row["genome_size"].isdigit() else None,
            )


def genome_table(release: str = "r232") -> dict[str, Genome]:
    """accession -> Genome, for the whole release (about 900k rows)."""
    return {g.accession: g for g in iter_genomes(release)}


def species_by_unversioned(release: str = "r232") -> dict[str, str]:
    """Unversioned accession -> species, for matching sources that drop versions."""
    return {unversioned(g.accession): g.species for g in iter_genomes(release) if g.species}


@dataclass(frozen=True)
class SpeciesCluster:
    species: str
    representative: str
    taxonomy: str
    ani_radius: float
    members: tuple[str, ...]


@functools.lru_cache(maxsize=4)
def species_clusters(release: str) -> dict[str, SpeciesCluster]:
    """`sp_clusters_rNNN.tsv` -> species -> cluster (representative + members).

    Parsed once per process: the R232 table lists 900,000 member genomes,
    and the strain stage asks for it once per clade.
    """
    path = GTDB_DIR / release / f"sp_clusters_{release}.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    out: dict[str, SpeciesCluster] = {}
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rep = bare(row["Representative genome"])
            members = tuple(bare(m) for m in row["Clustered genomes"].split(",") if m.strip())
            sp = row["GTDB species"].replace("s__", "", 1).strip()
            try:
                radius = float(row["ANI circumscription radius"])
            except (ValueError, KeyError):
                radius = 95.0
            out[sp] = SpeciesCluster(sp, rep, row["GTDB taxonomy"], radius, (rep, *members))
    return out


def representative_accessions(release: str = "r232") -> dict[str, str]:
    """species -> representative accession."""
    return {sp: c.representative for sp, c in species_clusters(release).items()}


def summary(release: str = "r232") -> Mapping[str, int | str]:
    clusters = species_clusters(release)
    return {
        "release": release, "species_clusters": len(clusters),
        "genomes": sum(len(c.members) for c in clusters.values()),
    }


if __name__ == "__main__":  # pragma: no cover
    build_slim("r232")
    print(json.dumps(summary("r232")))
