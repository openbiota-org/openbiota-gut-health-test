"""Name bridges the inventory needs to recognise one organism under two labels.

Two catalogues name the same population differently, and a merge that
works only on identical strings leaves duplicates: a MetaPhlAn 3 species
carrying its 2019 NCBI name beside the same organism's GTDB R232 record, or
a mOTUs cluster the profiler labels ``Unknown Eggerthella`` beside the
GlobDB genome that *is* that cluster. Each bridge here is built once from
the reference files on disk, cached under ``refs/crosswalk``, and answers
in constant time.

``ncbi_species_to_r232``
    NCBI organism names, as MetaPhlAn writes them (``Blautia_sp_CAG_257``),
    to the GTDB R232 species the genomes carrying that name were placed
    in. Majority rule over the release's genomes; a name whose genomes
    scatter maps to nothing and keeps its own label.

``motus_representative``
    mOTUs 4.1 cluster -> its representative genome and that genome's GTDB
    species (R226, walked to R232 by the caller). The profiler's own label
    is an 80 % majority vote over members and reads ``Unknown <genus>`` when
    members disagree; the representative genome is what GlobDB and sylph
    would call the cluster.

``motus_globdb_id``
    mOTUs 4.1 cluster -> the GlobDB r232 genome that represents it, through
    GlobDB's own dictionary of source genome names. A cluster GlobDB kept is
    the same unit under both catalogues; one it dereplicated away has no
    GlobDB id and is named by its representative's species instead.
"""

from __future__ import annotations

import functools
import gzip
import re
from pathlib import Path
from typing import Final

CROSSWALK_DIR: Final = Path("refs/crosswalk")
GTDB_DIR: Final = Path("refs/gtdb/r232")
MOTUS_DB: Final = Path("refs/motus/db/db_mOTU")
GLOBDB_DICTS: Final = Path("refs/globdb_r232/globdb_r232_dictionaries")

_NCBI_MAP_FILE: Final = CROSSWALK_DIR / "ncbi_names_to_gtdb_r232.tsv.gz"
_MOTU_ID = re.compile(r"(?:mOTUv4\.0_|MOTU40_)(\d+)")
_RSGB_ACC = re.compile(r"(GC[AF])-(\d{9})-V(\d+)")


def mpa_style(name: str) -> str:
    """``Blautia sp. CAG:257`` -> ``Blautia_sp_CAG_257`` (MetaPhlAn's spelling)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def _species_of(taxonomy: str) -> str:
    for seg in taxonomy.split(";"):
        seg = seg.strip()
        if seg.startswith("s__"):
            return seg[3:].strip()
    return ""


def build_ncbi_map(*, force: bool = False) -> Path:
    """Write the NCBI-name bridge from the R232 metadata (once; minutes)."""
    if _NCBI_MAP_FILE.is_file() and not force:
        return _NCBI_MAP_FILE
    from collections import Counter, defaultdict

    votes: dict[str, Counter[str]] = defaultdict(Counter)
    for fname in ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz"):
        path = GTDB_DIR / fname
        if not path.is_file():
            continue
        with gzip.open(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {h: i for i, h in enumerate(header)}
            i_tax, i_org = col["gtdb_taxonomy"], col["ncbi_organism_name"]
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(i_tax, i_org):
                    continue
                species = _species_of(parts[i_tax])
                organism = parts[i_org].strip()
                if not species or not organism or organism.lower() in ("none", "na", ""):
                    continue
                organism = organism.replace("[", "").replace("]", "")
                votes[mpa_style(organism)][species] += 1
                words = organism.split()
                # the binomial without strain designation, when it is one
                if len(words) >= 2 and words[1] not in ("sp.", "bacterium", "cf.") and words[0][:1].isupper():
                    votes[mpa_style(" ".join(words[:2]))][species] += 1
    CROSSWALK_DIR.mkdir(parents=True, exist_ok=True)
    tmp = _NCBI_MAP_FILE.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as out:
        for name, counter in sorted(votes.items()):
            total = sum(counter.values())
            species, n = counter.most_common(1)[0]
            if n * 2 >= total:
                out.write(f"{name}\t{species}\t{n}\t{total}\n")
    tmp.replace(_NCBI_MAP_FILE)
    return _NCBI_MAP_FILE


@functools.lru_cache(maxsize=1)
def ncbi_species_to_r232() -> dict[str, str]:
    """MetaPhlAn-style NCBI name -> GTDB R232 species (majority of that name's genomes)."""
    try:
        path = build_ncbi_map()
    except (FileNotFoundError, KeyError, OSError):
        return {}
    out: dict[str, str] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            name, _, rest = line.rstrip("\n").partition("\t")
            species = rest.split("\t", 1)[0]
            if name and species:
                out[name] = species
    return out


def motu_number(identifier: str) -> str | None:
    """``mOTUv4.0_078089`` or ``MOTU40_078089`` -> ``078089``."""
    m = _MOTU_ID.search(identifier or "")
    return m.group(1) if m else None


@functools.lru_cache(maxsize=1)
def motus_representatives() -> dict[str, tuple[str, str]]:
    """mOTU id -> (representative genome name, GTDB species of that genome, R226)."""
    path = MOTUS_DB / "mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz"
    if not path.is_file():
        return {}
    out: dict[str, tuple[str, str]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            motu, genome, tax = parts[col.get("MOTU", 0)], parts[col.get("GENOME", 1)], parts[col.get("GTDB", 2)]
            out[motu] = (genome, _species_of(tax))
    return out


@functools.lru_cache(maxsize=1)
def _globdb_motu_dictionary() -> dict[str, str]:
    """Source genome name -> GlobDB id, for the mOTUs genomes GlobDB kept."""
    path = GLOBDB_DICTS / "MOTU4_dictionary.tsv"
    if not path.is_file():
        return {}
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            name, _, gid = line.rstrip("\n").partition("\t")
            if name and gid:
                out[name] = gid
    return out


def motus_globdb_id(motu: str) -> str | None:
    """The GlobDB r232 genome representing this mOTUs cluster, if GlobDB kept it."""
    from openbiota.expansion import genomes

    rep = motus_representatives().get(motu)
    if rep is not None:
        gid = _globdb_motu_dictionary().get(rep[0])
        if gid and genomes.in_globdb(gid):
            return gid
    n = motu_number(motu)
    if n is None:
        return None
    # GlobDB numbers the clusters it kept by the mOTU number itself.
    candidate = f"MOTU40_{n}"
    return candidate if genomes.in_globdb(candidate) else None


def motus_genome_accession(motu: str) -> str | None:
    """A fetchable genome for the cluster: its GlobDB copy, or the NCBI assembly
    behind an ``RSGB`` representative (``RSGB23-1_GCF-016902095-V1_...`` ->
    ``GCF_016902095.1``)."""
    gid = motus_globdb_id(motu)
    if gid:
        return gid
    rep = motus_representatives().get(motu)
    if rep is None:
        return None
    m = _RSGB_ACC.search(rep[0])
    if m:
        return f"{m.group(1)}_{m.group(2)}.{m.group(3)}"
    return None


def resolve_motu(motu: str, label_species: str, label_genus: str) -> dict[str, str | bool | None]:
    """The inventory's name for a mOTUs observation.

    A cluster the profiler named (a GTDB species) keeps that name. One it
    labelled ``Unknown <genus>`` is named by its representative genome's
    species when that is a species; otherwise it is an unnamed cluster,
    labelled by its GlobDB id when GlobDB kept it (so it is one record with
    the sylph lane) and by its mOTU id when not.
    """
    unknown = (not label_species) or label_species.lower().startswith("unknown")
    genome = motus_genome_accession(motu)
    if not unknown:
        return {"species": label_species, "genus": label_genus, "unnamed": False, "genome_id": genome}
    rep = motus_representatives().get(motu)
    rep_species = rep[1] if rep else ""
    if rep_species and not rep_species.lower().startswith("unknown"):
        genus = rep_species.split(" ", 1)[0]
        return {"species": rep_species, "genus": genus, "unnamed": False, "genome_id": genome}
    # genuinely unnamed: the genus the vote could settle on, then the cluster id
    genus = label_genus or (rep_species.split(" ")[-1] if rep_species else "")
    if genus.lower().startswith("unknown"):
        genus = genus[len("unknown"):].strip()
    gid = motus_globdb_id(motu)
    cluster = gid or motu
    return {"species": f"{genus} {cluster}".strip(), "genus": genus, "unnamed": True, "genome_id": genome}
