# Attribution

Generated 2026-09-29 by `scripts/attribution_doc.py`. OpenBiota redistributes none of the
software or data below: each is fetched from its publisher at install time and pinned by version and
checksum under `refs/expanded/`. Every run records the exact versions and release identifiers it used in
`results.json` (`detection.reference_releases`, `run`) and in `run.log`. Licences are as each project publishes
them at the time of writing; check the source before any commercial use, and note that the MetaPhlAn
databases and BacDive carry noncommercial or contact-first terms of their own.

## Software

| tool | used for | licence | home |
|---|---|---|---|
| DIAMOND | translated search of every read against the metabolite gene panels | GPL-3.0 | <https://github.com/bbuchfink/diamond> |
| Bowtie2 | host-read removal; competitive confirmation alignments | GPL-3.0 | <https://github.com/BenLangmead/bowtie2> |
| samtools | alignment sorting, indexing and depth | MIT/Expat | <https://github.com/samtools/samtools> |
| MetaPhlAn 3.1 | the scoring lane the reference cohort and GMWI2 were built with | MIT | <https://github.com/biobakery/MetaPhlAn> |
| MetaPhlAn 4 (Jun23, Jan26), StrainPhlAn 4, sample2markers | species-level genome bins; comparative strain placement | MIT (software); CHOCOPhlAn databases CC BY-NC-SA | <https://github.com/biobakery/MetaPhlAn> |
| PhyloPhlAn, RAxML | marker alignment and tree inference inside StrainPhlAn | MIT; GPL-3.0 | <https://github.com/biobakery/phylophlan> |
| sylph | whole-genome containment against GTDB R232 and GlobDB r232 | MIT | <https://github.com/bluenote-1577/sylph> |
| mOTUs 4.1 | universal single-copy marker profiling | GPL-3.0 | <https://github.com/motu-tool/mOTUs> |
| Kraken 2, Bracken | read classification against UHGG v2.0.2 and the gut rescue panel | MIT; GPL-3.0 | <https://github.com/DerrickWood/kraken2> |
| SingleM (with OrfM, smafa, mfqe, HMMER) | marker-window profiling of lineages no catalogue names | GPL-3.0 (SingleM, HMMER); MIT (OrfM, smafa, mfqe) | <https://github.com/wwood/singlem> |
| metaSPAdes | targeted assembly of reads no method could place | GPL-2.0 | <https://github.com/ablab/spades> |
| skani, FastANI | genome-wide ANI for supplement reconciliation | MIT; Apache-2.0 | <https://github.com/bluenote-1577/skani> |
| seqkit | read selection for targeted assembly; marker extraction | MIT | <https://github.com/shenwei356/seqkit> |
| aria2 | parallel, resumable reference downloads | GPL-2.0-or-later | <https://aria2.github.io> |
| fastp | read QC and trimming | MIT | <https://github.com/OpenGene/fastp> |
| GMWI2 | the published gut-health index | MIT | <https://github.com/danielchang2002/GMWI2> |
| MICOM (with HiGHS, OSQP) | model-assisted metabolic scenarios | Apache-2.0; MIT; Apache-2.0 | <https://github.com/micom-dev/micom> |
| ReportLab | the PDF | BSD-3-Clause | <https://www.reportlab.com/opensource/> |
| NumPy, PyYAML, scikit-learn | arithmetic, configuration, the age model | BSD-3-Clause; MIT; BSD-3-Clause | <https://numpy.org> |
| curatedMetagenomicData | the 3,027-adult reference cohort's profiles and metadata | Artistic-2.0 | <https://waldronlab.io/curatedMetagenomicData/> |
| AGORA2 / VMH | genome-scale metabolic models behind the scenarios | CC BY 4.0 | <https://www.vmh.life> |
| BacDive (DSMZ) | culture-collection identities and reference-strain traits, via API v2 | CC BY 4.0; DSMZ asks to be contacted for commercial use | <https://bacdive.dsmz.de> |
| UniProt | reference proteins for the metabolite gene panels | CC BY 4.0 | <https://www.uniprot.org> |

## Reference releases

Each lock file under `refs/expanded/locks/` carries the URLs, sizes, SHA-256 digests and retrieval date.

| source | release | licence | home |
|---|---|---|---|
| elgg | ELGG (Zenodo 6969520) (2,172 species representative) | CC BY 4.0 (Zenodo record 6969520) | <https://zenodo.org/records/6969520/files> |
| globdb_genomes | globdb_r232 (346,233 representative genome) | CC BY-SA 4.0 | <https://fileshare.lisc.univie.ac.at/globdb/globdb_r232> |
| globdb_r232 | globdb_r232 (346,233 species representative) | CC BY-SA 4.0 (GlobDB); index by sylph authors | <https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/taxonomic_profiling> |
| gtdb_r207 | R207 | CC BY-SA 4.0 | <https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files> |
| gtdb_r220 | R220 | CC BY-SA 4.0 | <https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files> |
| gtdb_r226 | R226 | CC BY-SA 4.0 | <https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files> |
| gtdb_r232 | R232 (901,341 genome) | CC BY-SA 4.0 | <https://data.gtdb.ecogenomic.org/releases/release232/232.0> |
| hrgm2 | HRGMv2 (4,824 species representative) | CC BY 4.0 (Zenodo record 19482781) | <https://zenodo.org/records/19480672/files> |
| hrom | HROM | see file | <https://www.decodebiome.org/HROM/data/genome_catalog> |
| humgut2 | HumGut2 (31,225 genome (97.5% ANI dereplicated)) | CC BY 4.0 (authors) | <https://arken.nmbu.no/~larssn/humgut> |
| motus_db | mOTUs DB 4.1 (Zenodo 20322482) | CC BY 4.0 (mOTUs data) | <https://zenodo.org/records/20322482/files> |
| mpa_jan26 | mpa_vJan26_CHOCOPhlAnSGB_202605 (72,000 SGB) | MIT (MetaPhlAn); database CC BY-NC-SA per bioBakery | <https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases> |
| sgb_bridges | MetaPhlAn utils @ master 2026-09 | MIT (MetaPhlAn) | <https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils> |
| singlem_globdb | GlobDB_r232.metapackage_v4 | CC BY-SA 4.0 | <https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/taxonomic_profiling> |
| uhgg_kraken2 | v2.0.2 (4,744 species) | EMBL-EBI terms of use | <https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2> |
| uhgg_v2.0.2 | v2.0.2 (289,232 genome) | EMBL-EBI terms of use (open) | <https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2> |

## Published research

Every metabolite pathway, disease-pattern profile, microbial group, organism description and intervention
assertion names the study it rests on (DOI or PubMed identifier) in its definition file and on the page of
the report where it is used. The complete catalogue of sources is on the website's research page and in
`specs/research/`.
