"""The expanded reference sources, declared once, with everything a lock needs.

Spec 0.8.4 §3: every reference a sample run touches is acquired once,
resumably, with its release, URL, publisher checksum where one is
published, local SHA-256, retrieval time, size, taxonomy release and
licence recorded. A sample run never resolves a moving alias: the URLs here
name fixed releases, and `refs/expanded/locks/*.lock.json` is what a run
reads to know what it is looking at.

Adding a source: append a `RefFile` here. `scripts/fetch_expanded_refs.py`
does the rest, and `make refs-expanded` calls it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

#: Where every expanded reference lives, under `refs/`.
EXPANDED_DIR: Final = Path("refs/expanded")
LOCK_DIR: Final = EXPANDED_DIR / "locks"
PROGRESS_FILE: Final = EXPANDED_DIR / "fetch.progress"


@dataclass(frozen=True)
class RefFile:
    """One file of one reference release."""

    source_id: str          # R1..R10 family, e.g. "mpa_jan26"
    release_id: str         # fixed release string, e.g. "mpa_vJan26_CHOCOPhlAnSGB_202605"
    url: str
    dest: str               # path under refs/
    license: str
    taxonomy_release: str = ""
    #: Publisher checksum: "md5:<hex>" or "md5-list:<url>" (looked up by basename).
    publisher_checksum: str = ""
    count_unit: str = ""
    source_count: int | None = None
    note: str = ""
    #: Skip in the default fetch (very large; fetched by a later, targeted step).
    deferred: bool = False
    tags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def path(self) -> Path:
        return Path("refs") / self.dest


_MPA = "https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/"
_GLOB = "https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/"
_GLOB_MD5 = _GLOB + "globdb_r232_md5sum.txt"
_GLOB_TAX_MD5 = _GLOB + "taxonomic_profiling/globdb_r232_tax_profile_md5sums.txt"
_UHGG = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/"
_GTDB = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
_GTDB_MD5 = _GTDB + "MD5SUM.txt"
_HROM = "https://www.decodebiome.org/HROM/data/genome_catalog/"

JAN26: Final = "mpa_vJan26_CHOCOPhlAnSGB_202605"

SOURCES: Final[tuple[RefFile, ...]] = (
    # ----------------------------------------------------------- R1 MetaPhlAn Jan26
    RefFile("mpa_jan26", JAN26, _MPA + f"{JAN26}.tar", f"metaphlan4_db/{JAN26}.tar",
            license="MIT (MetaPhlAn); database CC BY-NC-SA per bioBakery", taxonomy_release="SGB Jan26 (GTDB R226 bridge upstream)",
            publisher_checksum=f"md5-list:{_MPA}{JAN26}.md5", count_unit="SGB", source_count=72_000,
            note="marker database; bowtie2 index is built locally by `metaphlan --install`"),
    RefFile("mpa_jan26", JAN26, _MPA + f"{JAN26}.md5", f"metaphlan4_db/{JAN26}.md5", license="as above"),
    RefFile("mpa_jan26", JAN26, _MPA + f"bowtie2_indexes/{JAN26}_bt2.tar", f"metaphlan4_db/{JAN26}_bt2.tar",
            license="as above", publisher_checksum=f"md5-list:{_MPA}bowtie2_indexes/{JAN26}_bt2.md5",
            note="publisher-built bowtie2 index; `metaphlan --install` fetches it and the lock is written post hoc by "
                 "scripts/lock_installed_refs.py, or by this fetcher when run with --only mpa_jan26 --include-deferred",
            deferred=True),
    RefFile("mpa_jan26", JAN26, _MPA + f"{JAN26}_species.txt.bz2", f"metaphlan4_db/{JAN26}_species.txt.bz2",
            license="as above", count_unit="SGB", source_count=72_000,
            publisher_checksum="sha256:441f86af40f587b969b40ebdb79ffb02b311ffec4c81024d8001e529c6ea5572"),
    RefFile("mpa_jan26", JAN26, _MPA + f"{JAN26}_marker_info.txt.bz2", f"metaphlan4_db/{JAN26}_marker_info.txt.bz2", license="as above"),
    RefFile("mpa_jan26", JAN26, _MPA + f"{JAN26}.nwk", f"metaphlan4_db/{JAN26}.nwk", license="as above"),
    # ----------------------------------------------------------- R3 GlobDB r232 (sylph two-stage index + publisher taxonomy)
    RefFile("globdb_r232", "globdb_r232", "http://faust.compbio.cs.cmu.edu/sylph-stuff/globdb_r232_sylph_v2_c200.syl2db",
            "sylph/globdb_r232/globdb_r232_sylph_v2_c200.syl2db", license="CC BY-SA 4.0 (GlobDB); index by sylph authors",
            taxonomy_release="GlobDB r232", count_unit="species representative", source_count=346_233,
            note="official sylph pre-built list; no publisher hash - integrity established by loading and by matching the publisher taxonomy"),
    RefFile("globdb_r232", "globdb_r232", _GLOB + "taxonomic_profiling/globdb_r232_taxonomy_sylph.tsv.gz",
            "sylph/globdb_r232/globdb_r232_taxonomy_sylph.tsv.gz", license="CC BY-SA 4.0", taxonomy_release="GlobDB r232",
            publisher_checksum=f"md5-list:{_GLOB_TAX_MD5}"),
    RefFile("globdb_r232", "globdb_r232", _GLOB + "globdb_r232_dictionaries.tar.gz", "globdb_r232/globdb_r232_dictionaries.tar.gz",
            license="CC BY-SA 4.0", publisher_checksum=f"md5-list:{_GLOB_MD5}", note="source-identifier crosswalk dictionaries"),
    RefFile("globdb_r232", "globdb_r232", _GLOB + "globdb_r232_taxonomy.tsv.gz", "globdb_r232/globdb_r232_taxonomy.tsv.gz",
            license="CC BY-SA 4.0", publisher_checksum=f"md5-list:{_GLOB_MD5}"),
    RefFile("globdb_r232", "globdb_r232", _GLOB + "globdb_r232_checkm2.tsv.gz", "globdb_r232/globdb_r232_checkm2.tsv.gz",
            license="CC BY-SA 4.0", publisher_checksum=f"md5-list:{_GLOB_MD5}", note="genome quality for the quality policy"),
    RefFile("globdb_r232", "globdb_r232", _GLOB + "globdb_r232_dataset_list.tsv", "globdb_r232/globdb_r232_dataset_list.tsv",
            license="CC BY-SA 4.0", publisher_checksum=f"md5-list:{_GLOB_MD5}", note="the 26 source datasets"),
    RefFile("globdb_genomes", "globdb_r232", _GLOB + "globdb_r232_genome_fasta.tar.gz", "globdb_r232/globdb_r232_genome_fasta.tar.gz",
            license="CC BY-SA 4.0", publisher_checksum=f"md5-list:{_GLOB_MD5}", count_unit="representative genome",
            source_count=346_233, deferred=True,
            note="272 GB one-time archive: GlobDB publishes no per-genome endpoint, and competitive confirmation "
                 "of a GlobDB-native cluster needs its genome. Subsets are extracted on demand into globdb_r232/genomes/"),
    # ----------------------------------------------------------- SingleM + GlobDB metapackage (§4D)
    RefFile("singlem_globdb", "GlobDB_r232.metapackage_v4", _GLOB + "taxonomic_profiling/GlobDB_r232.metapackage_v4.smpkg.tar.gz",
            "singlem/GlobDB_r232.metapackage_v4.smpkg.tar.gz", license="CC BY-SA 4.0", taxonomy_release="GlobDB r232",
            publisher_checksum=f"md5-list:{_GLOB_TAX_MD5}"),
    # ----------------------------------------------------------- R2 GTDB R232 metadata (reconciliation backbone)
    RefFile("gtdb_r232", "R232", _GTDB + "bac120_metadata_r232.tsv.gz", "gtdb/r232/bac120_metadata_r232.tsv.gz",
            license="CC BY-SA 4.0", taxonomy_release="GTDB R232", publisher_checksum=f"md5-list:{_GTDB_MD5}",
            count_unit="genome", source_count=901_341),
    RefFile("gtdb_r232", "R232", _GTDB + "ar53_metadata_r232.tsv.gz", "gtdb/r232/ar53_metadata_r232.tsv.gz",
            license="CC BY-SA 4.0", taxonomy_release="GTDB R232", publisher_checksum=f"md5-list:{_GTDB_MD5}"),
    RefFile("gtdb_r232", "R232", _GTDB + "bac120_taxonomy_r232.tsv.gz", "gtdb/r232/bac120_taxonomy_r232.tsv.gz",
            license="CC BY-SA 4.0", taxonomy_release="GTDB R232", publisher_checksum=f"md5-list:{_GTDB_MD5}"),
    RefFile("gtdb_r232", "R232", _GTDB + "ar53_taxonomy_r232.tsv.gz", "gtdb/r232/ar53_taxonomy_r232.tsv.gz",
            license="CC BY-SA 4.0", taxonomy_release="GTDB R232", publisher_checksum=f"md5-list:{_GTDB_MD5}"),
    RefFile("gtdb_r232", "R232", _GTDB + "auxillary_files/sp_clusters_r232.tsv", "gtdb/r232/sp_clusters_r232.tsv",
            license="CC BY-SA 4.0", taxonomy_release="GTDB R232", note="species clusters with ANI radii"),
    RefFile("gtdb_r232", "R232", _GTDB + "VERSION.txt", "gtdb/r232/VERSION.txt", license="CC BY-SA 4.0"),
    RefFile("gtdb_r232", "R232", _GTDB + "MD5SUM.txt", "gtdb/r232/MD5SUM.txt", license="CC BY-SA 4.0"),
    # ----------------------------------------------------------- release-aware bridges (§5)
    RefFile("sgb_bridges", "MetaPhlAn utils @ master 2026-09", "https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv",
            "crosswalk/upstream/mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv", license="MIT (MetaPhlAn)", taxonomy_release="GTDB R226",
            note="upstream SGB->GTDB bridge for the Jan26 database; an intermediate step, never a label"),
    RefFile("sgb_bridges", "MetaPhlAn utils @ master 2026-09", "https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB_r207.tsv",
            "crosswalk/upstream/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB_r207.tsv", license="MIT (MetaPhlAn)", taxonomy_release="GTDB R207",
            note="the correct bridge for the installed Jun23 database (the repo had been applying the Jan25/R220 file)"),
    RefFile("sgb_bridges", "MetaPhlAn utils @ master 2026-09", "https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB_r220.tsv",
            "crosswalk/upstream/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB_r220.tsv", license="MIT (MetaPhlAn)", taxonomy_release="GTDB R220",
            note="kept only to reproduce and audit the historical misuse"),
    RefFile("gtdb_r226", "R226", "https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/sp_clusters_r226.tsv",
            "gtdb/r226/sp_clusters_r226.tsv", license="CC BY-SA 4.0", taxonomy_release="GTDB R226",
            note="representative accession per R226 species: the assembly membership that resolves R226 names onward to R232"),
    RefFile("gtdb_r220", "R220", "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/sp_clusters_r220.tsv",
            "gtdb/r220/sp_clusters_r220.tsv", license="CC BY-SA 4.0", taxonomy_release="GTDB R220",
            note="HRGM2 and HumGut2 carry r220-era species names; their representative accessions walk to R232 through this"),
    RefFile("gtdb_r207", "R207", "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/sp_clusters_r207.tsv",
            "gtdb/r207/sp_clusters_r207.tsv", license="CC BY-SA 4.0", taxonomy_release="GTDB R207",
            note="for the historical Jun23 reconciliation repair"),
    # ----------------------------------------------------------- R4 mOTUs 4.1 marker database
    RefFile("motus_db", "mOTUs DB 4.1 (Zenodo 20322482)", "https://zenodo.org/records/20322482/files/db_mOTU.tar.gz?download=1",
            "motus/db_mOTU.tar.gz", license="CC BY 4.0 (mOTUs data)", taxonomy_release="mOTUs 4.1 (GTDB-anchored)",
            publisher_checksum="md5:471ea128f0c0839f5c4629b949ea5f8a", count_unit="mOTU (species-level unit)",
            source_count=None, note="count recorded from the downloaded artifact, not from either published figure"),
    RefFile("motus_db", "mOTUs DB 4.1 (Zenodo 20343612)", "https://zenodo.org/records/20343612/files/mOTUsv4.1.annotation.db?download=1",
            "motus/mOTUsv4.1.annotation.db", license="CC BY 4.0", publisher_checksum="md5:5cb1a220f28b427bd0d7b907a415eccc",
            note="18 GB annotation database (genome context); not needed for profiling", deferred=True),
    # ----------------------------------------------------------- R6 UHGG v2.0.2 (Kraken2 DB + metadata)
    RefFile("uhgg_v2.0.2", "v2.0.2", _UHGG + "genomes-all_metadata.tsv", "supplements/uhgg_v2.0.2/genomes-all_metadata.tsv",
            license="EMBL-EBI terms of use (open)", taxonomy_release="GTDB (as annotated by MGnify)",
            count_unit="genome", source_count=289_232),
    RefFile("uhgg_v2.0.2", "v2.0.2", _UHGG + "README_v2.0.2.txt", "supplements/uhgg_v2.0.2/README_v2.0.2.txt", license="EMBL-EBI terms of use"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/hash.k2d", "kraken2/uhgg_v2.0.2/hash.k2d",
            license="EMBL-EBI terms of use", count_unit="species", source_count=4_744,
            note="publisher-built Kraken2 database over the UHGG v2.0.2 species representatives"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/opts.k2d", "kraken2/uhgg_v2.0.2/opts.k2d", license="EMBL-EBI terms of use"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/taxo.k2d", "kraken2/uhgg_v2.0.2/taxo.k2d", license="EMBL-EBI terms of use"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/seqid2taxid.map", "kraken2/uhgg_v2.0.2/seqid2taxid.map", license="EMBL-EBI terms of use"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/database150mers.kmer_distrib", "kraken2/uhgg_v2.0.2/database150mers.kmer_distrib",
            license="EMBL-EBI terms of use", note="Bracken distribution matching 150 bp reads"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/database100mers.kmer_distrib", "kraken2/uhgg_v2.0.2/database100mers.kmer_distrib",
            license="EMBL-EBI terms of use", note="Bracken distribution for shorter reads"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/taxonomy/names.dmp", "kraken2/uhgg_v2.0.2/taxonomy/names.dmp", license="EMBL-EBI terms of use"),
    RefFile("uhgg_kraken2", "v2.0.2", _UHGG + "kraken2_db_uhgg_v2.0.2/taxonomy/nodes.dmp", "kraken2/uhgg_v2.0.2/taxonomy/nodes.dmp", license="EMBL-EBI terms of use"),
    # ----------------------------------------------------------- R5 HRGM2
    RefFile("hrgm2", "HRGMv2", "https://zenodo.org/records/19482781/files/HRGMv2_Rep_Genome.tar.gz?download=1",
            "supplements/hrgm2/HRGMv2_Rep_Genome.tar.gz", license="CC BY 4.0 (Zenodo record 19482781)",
            taxonomy_release="GTDB r220 (authors' assignment)", count_unit="species representative", source_count=4_824),
    RefFile("hrgm2", "HRGMv2", "https://zenodo.org/records/19480672/files/HRGMv2_Cluster_metadata.tsv?download=1",
            "supplements/hrgm2/HRGMv2_Cluster_metadata.tsv", license="CC BY 4.0 (Zenodo record 19480672)"),
    RefFile("hrgm2", "HRGMv2", "https://zenodo.org/records/19480672/files/HRGMv2_gtdbr220_results.tsv?download=1",
            "supplements/hrgm2/HRGMv2_gtdbr220_results.tsv", license="CC BY 4.0"),
    RefFile("hrgm2", "HRGMv2", "https://zenodo.org/records/19480672/files/Dereplication_genomes_metadata.tsv?download=1",
            "supplements/hrgm2/Dereplication_genomes_metadata.tsv", license="CC BY 4.0", count_unit="genome", source_count=230_632),
    # ----------------------------------------------------------- R7 ELGG
    RefFile("elgg", "ELGG (Zenodo 6969520)", "https://zenodo.org/records/6969520/files/ELGG_representatives_2172.zip?download=1",
            "supplements/elgg/ELGG_representatives_2172.zip", license="CC BY 4.0 (Zenodo record 6969520)",
            count_unit="species representative", source_count=2_172,
            note="early-life discovery cohort; sequence coverage only, not a paediatric reference interval"),
    # ----------------------------------------------------------- R8 HumGut2
    RefFile("humgut2", "HumGut2", "https://arken.nmbu.no/~larssn/humgut/HumGut2.tsv", "supplements/humgut2/HumGut2.tsv",
            license="CC BY 4.0 (authors)", count_unit="genome (97.5% ANI dereplicated)", source_count=31_225,
            note="metadata first; genomes are fetched per accession from R10 after reconciliation - the 18 GB archive is not mirrored"),
    # ----------------------------------------------------------- R9 HROM
    RefFile("hrom", "HROM", _HROM + "HROM-Species-metadata.tsv", "supplements/hrom/HROM-Species-metadata.tsv",
            license="see HROM README (PRJNA1206836)", count_unit="oral species", source_count=3_426),
    RefFile("hrom", "HROM", _HROM + "HROM_Conspecific-genomes-metadata.tsv", "supplements/hrom/HROM_Conspecific-genomes-metadata.tsv",
            license="see HROM README", count_unit="genome", source_count=72_641),
    RefFile("hrom", "HROM", _HROM + "README.txt", "supplements/hrom/README.txt", license="see file"),
)

#: Licence terms the report must carry beside the results (spec §6).
LICENCE_NOTES: Final[dict[str, str]] = {
    "globdb_r232": "GlobDB propagates CC BY-SA 4.0.",
    "motus": "mOTUs data are CC BY 4.0.",
    "bacdive": "BacDive is CC BY 4.0; DSMZ asks to be contacted for commercial use.",
    "gtdb_r232": "GTDB is CC BY-SA 4.0.",
}


def by_source() -> dict[str, list[RefFile]]:
    out: dict[str, list[RefFile]] = {}
    for f in SOURCES:
        out.setdefault(f.source_id, []).append(f)
    return out
