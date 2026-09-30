# BUILD_SPEC_v08.1 — Expanded organism discovery, trustworthy abundance, and complete interpretation

**Project:** OpenBiota. **Research cutoff:** 2026-09-19. **Status:** implementation specification; performance targets below have not yet been demonstrated on the project's FASTQs.

**Handoff:** This document is the complete handoff for this increment. All implementation requirements, source locations, evidence boundaries, data contracts, migration rules, and acceptance criteria are here. No companion research report, source ledger, or earlier build-spec document is required. The coding agent must inspect the existing repository and adapt the proposed module names to its conventions. The repository and FASTQs were not supplied for this research audit; paths and commands introduced here are contracts to implement, not claims about code already present.

## 1. Required outcome

Build a substantially more sensitive and better explained organism inventory from the existing stool shotgun FASTQs. Expand beyond the current union of legacy MetaPhlAn profiles using current marker and whole-genome references, credible low-abundance confirmation, explicit taxonomy reconciliation, and complete organism cards.

The product must:

1. Search the current, quality-controlled bacterial/archaeal reference universe and a substantially broader complementary genome universe. Include gut-specific, early-life, geographically diverse, fungal, viral, and eukaryotic resources where relevant.
2. Demonstrate recovery of **at least 570 of 600 resolvable bacterial/archaeal species in the defined 8-million-pair benchmark, and at least 588 of 600 at 30 million pairs**, while meeting the precision gates in §14. These are analytical capability targets, not a requirement to report 600 organisms in every person.
3. Count independent organisms at an evidence-supported resolution. A name change, an additional strain, a synonymous reference, a gene, a plasmid, or another database hit must not inflate species richness.
4. Eliminate cross-profiler abundance arithmetic errors. A composition must have one estimator, one denominator, one coherent taxon partition, and an explicit treatment of unclassified/unresolved mass.
5. Give **100% of reported organisms** a readable identity card. Launch a curated priority registry covering **at least 600 canonical species/cluster concepts**, with explicit unknowns where evidence is absent. Six hundred cards does not mean six hundred proven health effects.
6. Preserve the calibrated legacy scoring lane, existing pathogen targets, strain analyses, functional panels, biofilm modules, and donor/recipient analyses while adding the new inventory. Revalidate any downstream calculation that changes its inputs.
7. Explain what the assay searched, supported, could not distinguish, or could not assess. Do not imply that more database entries automatically mean better measured sensitivity.

There may be trillions of microbial cells in a gut, but cell count is not species richness. The recoverable species count depends on biology, sampling, extraction, sequencing depth, reference representation, and the operating point for false positives. Improve all computationally addressable limits; do not manufacture a desired count.

## 2. Confirmed defects in the supplied OpenBiota report

The audited artifact is the 219-page SAMPLE2_A02 report generated on 2026-09-19. These are PDF observations, not a repository audit.

| Observed result | Required correction |
|---|---|
| Pages 8–12 list 186 inventory rows: 142 named/non-bare-SGB rows plus 44 bare SGB identifiers. | Reconcile identities before declaring a biological species count; replace bare IDs with supported readable taxonomy while preserving IDs. |
| The displayed inventory percentages total **130.777%**: 112.837% in the named list plus 17.940% in the bare-SGB list. | Find the generating merge logic. Never add, maximize, or average abundances from independent compositions. Do not merely divide all 186 values by 130.777. |
| Page 12 calls the percentages a share of classified reads; page 218 says the composition sums to 100%. | Use the estimator's actual quantity. Marker abundance is not a directly counted read fraction. Enforce mass conservation before rendering. |
| There are 104 rows without a percentile, while a footer refers to only 10. | Generate all coverage counts and footers from the rendered records, with compatible-reference status. |
| Summary richness is 92 in one location; another section describes 159 MP4 species/164 SGBs; the merged inventory has 186 rows. | Show separate legacy-scoring richness, native marker units, and reconciled inventory concepts. Use one clearly defined headline count. |
| Page 212 reports 65 overlapping, 94 MP4-only and 27 MP3-only entries. | Recompute overlap using accession/cluster mappings, including splits, merges, and unresolved complexes; names alone cannot establish overlap. |
| The report calls the Jun23 marker catalog approximately 26,000 units. | Its published size is 36,822 SGBs. Update catalog metadata automatically; the newer Jan26 catalog contains 72,000 SGBs. |
| Pages 3–4 and 212 present different assigned/unrepresented estimates without a clear common denominator. | Separate actual fragment accounting, marker-based estimates, and whole-genome modeled sequence coverage. |
| The report has 52 strain characterizations. | Preserve this capability; explain whether each is marker consensus, comparative strain placement, named reference matching, or a validated pathotype call. |
| A Paraprevotella clara entry has 0.006% abundance while its strain section reports substantial marker support, including 153 markers and 4.16× marker depth. | Trace lane, reference, units, and sample identity. This is an investigation trigger, not proof that either number alone is wrong. |
| Input, passing, and downstream read counts differ; the report describes a fast assessment mode using original reads downstream. | Record assessment versus filtering explicitly, including reads versus pairs and the exact file used by every lane. Do not label assessed reads as filtered reads. |
| Claims include effectively searching every reference catalog and treating single-catalog calls as necessarily benign coverage differences. | Describe installed, admitted catalogs only. Single-lane calls can reflect reference differences, thresholds, ambiguity, or false positives and require evidence-specific adjudication. |

The coding agent's stated 113 interpretation records, 3,027 reference adults, 22 studies, and 486 pathogen targets are starting inventory claims to verify against the repository. Do not overwrite a newer local manifest to force these historical counts.

## 3. Architecture and invariants

### 3.1 Separate six concerns

1. **Observation:** raw lane output, exact input files, database, parameters, and supporting reads/markers.
2. **Identity:** versioned reference clusters, assemblies, nomenclature, taxonomy mappings, and supported rank.
3. **Detection:** accepted, provisional, ambiguous, unsupported, or not assessed under a validated operating point.
4. **Quantification:** abundance estimate, estimator, denominator, uncertainty, and unresolved mass.
5. **Interpretation:** organism traits, population associations, measured functions, strain-specific findings, and intervention evidence with provenance.
6. **Reference scoring:** percentiles/models requiring the same feature space and processing as their training/reference cohort.

An accepted detection need not have a compatible percentile. A taxon can have a readable name without having a known health role. A species can be present without its pathogenic strain or virulence determinant being present. A detected DNA sequence does not establish viability or persistent colonization.

### 3.2 Execution modes

| Mode | Mandatory work | Purpose |
|---|---|---|
| `comprehensive` — new default | Frozen scoring lane; Jan26 marker inventory; validated broad whole-genome profile; complementary search; read-level adjudication; current pathogen/strain/function modules; kingdom-specific search; complete rendering | Maximize supported recovery from the available data. |
| `standard` | Frozen scoring lane; Jan26 marker inventory; GTDB whole-genome profile; adjudication of new/ambiguous calls; preserved existing modules | Explicit resource-limited alternative, never silently substituted for comprehensive. |
| `deep-discovery` | Comprehensive plus assembly/binning and protein-homology rescue where depth supports them | Investigate unrepresented organisms and sequences. It must not promise strain assemblies or hundreds of MAGs from shallow libraries. |
| `benchmark` | All selected alternatives on truth sets, fixed thresholds, complete performance export | Choose and qualify production engines without running every redundant engine on every sample. |

During development, compare GTDB and GlobDB profiles on the same samples. After qualification, use **one** selected whole-genome abundance universe per report. Comprehensive should select the reconciled GlobDB-based universe if it passes the precision, incremental-recall, and quantification gates; otherwise retain the qualified GTDB universe and expose GlobDB-supported additions separately until the broader abundance model qualifies. This is an explicit staged release, not permission to omit the broader detection work.

### 3.3 Frozen scoring lane

Keep the existing MetaPhlAn 3.1/database/reference-cohort environment immutable for existing scores. Verify its exact database string, cohort participants, preprocessing, model versions, and file checksums. Inventory upgrades cannot change its feature vector or denominator.

Do not feed Jan26, GTDB, GlobDB, Kraken/Bracken, renamed placeholders, or merged abundances into MP3-trained percentiles, GMWI2, disease coefficients, or donor metrics without an explicit validated migration. Build a future reference lane by reprocessing eligible reference FASTQs uniformly. Taxonomic string conversion alone does not harmonize measurement methods.

### 3.4 Additive migration

Keep legacy raw outputs, cached FASTQs, marker mapping files, and strain results. New references get new cache keys and output directories. Reuse a cache only when input content hashes, preprocessing, tool/container, database files, parameters, and schema match. Preserve lineage between old and new results.

Existing `*.bowtie2.bz2`/mapping outputs are not substitutes for SAM. If the correct SAM is missing, realign the clean FASTQs to the selected marker database. Retain both mapping output and SAM for future strain work.

## 4. Reference catalog registry: core and complementary genomes

All counts below describe specific versions, not unique organisms after combining catalogs. “Current” means checked by 2026-09-19. Resolve fixed artifacts and hashes during reference preparation; never use `latest` or a moving branch during a sample run.

### 4.1 Production priorities

| Resource | Verified size and release | Required use, access, and important qualifications |
|---|---|---|
| **GTDB R11-RS232 / 232.0**, 2026-04-15 | **901,341 genomes; 199,923 species clusters**: 189,801 bacterial and 10,122 archaeal clusters | Main stable bacterial/archaeal taxonomic backbone. All environments, not 199,923 gut species. Use species representatives for profiling; retain members for strain/near-neighbor confirmation. [Release](https://data.gtdb.ecogenomic.org/releases/release232/232.0/), [project](https://gtdb.ecogenomic.org/), [terms](https://gtdb.ecogenomic.org/downloads). Data CC BY-SA 4.0. |
| **MetaPhlAn Jan26** `mpa_vJan26_CHOCOPhlAnSGB_202605` | **72,000 unique SGB records**, directly counted from the species metadata | Required marker inventory upgrade, paired with MetaPhlAn 4.2.6. SGBs are not necessarily named species. [Metadata](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan26_CHOCOPhlAnSGB_202605_species.txt.bz2), [database directory](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/), [software releases](https://github.com/biobakery/MetaPhlAn/releases). |
| **GlobDB r232**, 2026-06-25 | **346,233 representative clusters from 26 source datasets**; includes all 199,923 GTDB representatives and 146,310 additional representatives | Broad candidate and comprehensive profiling universe, subject to qualification. Its between-source clustering uses 96% ANI/50% alignment fraction, with GTDB priority; do not relabel every representative as an independent GTDB species. [Downloads](https://globdb.org/downloads), [methods](https://globdb.org/methods), [fixed release](https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/), [source manifest](https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_dataset_list.tsv). CC BY-SA 4.0 plus retained upstream provenance. |
| **HRGM2** | **155,211 nonredundant near-complete genomes / 4,824 species**, 41 countries | Required gut-reference gap audit; include absent clusters and useful close-relative/strain members. Near-complete criteria ≥90% completeness/≤5% contamination do not alone establish full MIMAG high quality. GTDB r220 labels require versioned mapping. [Paper](https://doi.org/10.1038/s41564-025-02206-1), [author download index](https://github.com/netbiolab/HRGM2/blob/main/DATA.md), [representatives](https://zenodo.org/records/19482781), [metadata](https://zenodo.org/records/19480672). Representative deposit CC0; preserve per-artifact terms. |
| **UHGG v2.0.2**, 2024-01-23 | **289,232 genomes / 4,744 species representatives** | Required complementary gut/strain panel and provenance bridge. Source taxonomy GTDB r202. Do not profile against all redundant genomes as independent species. [API metadata](https://www.ebi.ac.uk/metagenomics/api/v1/genome-catalogues/human-gut-v2-0-2), [fixed files](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/), [foundational paper](https://doi.org/10.1038/s41587-020-0603-3). Retain original reuse conditions. |
| **HumGut2**, metadata last modified 2025-06-11; inspected 2026-09-19 | **31,225 genomes at the 97.5% dereplication tier; 6,914 distinct 95% clusters** in the inspected metadata | Healthy-gut prevalence-informed complement. Its custom numeric IDs are not necessarily NCBI taxids. Use species clusters for counting, genome members for discrimination. [Project](https://github.com/larssnip/HumGut), [metadata](https://arken.nmbu.no/~larssn/humgut/HumGut2.tsv), [archive](https://arken.nmbu.no/~larssn/humgut/HumGut2.tar). Confirm dataset redistribution terms before bundling. |
| **ELGG** | **32,277 genomes / 2,172 species**, children under three in the published resource | Include an early-life representation audit, especially for younger participants; not an age-matched reference interval for children aged 8–18. [Paper](https://doi.org/10.1038/s41467-022-32805-z), [representatives/data](https://zenodo.org/records/6969520). Deposit CC BY 4.0. |
| **IMGG**, long-read gut collection | **6,729 genomes: 802 closed genomes and 5,927 high-quality MAGs**, not 6,729 species | Strong structural/reference completeness and strain validation complement. [Paper](https://doi.org/10.1038/s41564-022-01270-1), [data](https://doi.org/10.6084/m9.figshare.19661523), [reads](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA763692). Verify selected deposit terms. |
| **Human gut archaeal catalog** | **1,167 genomes** in the published catalog, including 15 then-new species | Audit archaeal representation and extraction/detection limitations independently of bacteria. [Paper](https://doi.org/10.1038/s41564-021-01020-9). Retrieve its accession/supplement manifest and archive it; do not add genome count to species richness. |
| **Hadza/industrialization-diverse collection** | **91,662 mixed-kingdom genomes from 351 fecal samples** in the 2023 study | Audit geographic and lifestyle gaps; 44% novelty was relative to catalogs at publication, not today's novelty. [Paper and data links](https://doi.org/10.1016/j.cell.2023.05.046). Reconcile with GlobDB/GTDB before admitting additions. |
| **AWI-Gen African gut data** | Published cohort includes 1,801 women and 19 men; reported 1,005 new MAGs | Representation audit and optional read reprocessing; “new MAG” is not “new species.” [Primary study](https://doi.org/10.1038/s41586-024-08485-8). Retain cohort/consent/access conditions. |
| **Expanded gut genome collection, 2022** | **241,118 genomes / 3,594 species**, 310 then-novel species | Historical complement, mostly overlapping newer sources. [Paper](https://doi.org/10.1038/s41467-022-31502-1), [representative deposit](https://doi.org/10.6084/m9.figshare.16885261). Some full data may require request; do not block the open core on them. |
| **NCBI RefSeq/GenBank assemblies** | Rolling, accession-versioned collection; no fixed gut-species count claimed | Fill identified sequence gaps, supply near-neighbors and type/reference strains, and preserve assembly quality and source metadata. [Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/), [assembly](https://www.ncbi.nlm.nih.gov/assembly/), [BioSample](https://www.ncbi.nlm.nih.gov/biosample/docs/). |

GlobDB already incorporates substantial contributions from mOTU, SPIRE, gcMeta, NGDC, GEM, MGnify, HRGM and other collections. Its actual 26-source manifest is authoritative; older articles/docs describing 14 sources are not the r232 manifest. Produce a **marginal coverage table** for every proposed addition: admitted assemblies, already represented clusters, genuinely additional clusters, better-quality replacements, and new strain/near-neighbor references. Do not concatenate all FASTAs and call their summed size the detection universe.

### 4.2 Integrity anchors verified in this audit

| Artifact | Integrity anchor |
|---|---|
| Jan26 `species.txt.bz2` linked above | 1,238,616 compressed bytes; SHA-256 `441f86af40f587b969b40ebdb79ffb02b311ffec4c81024d8001e529c6ea5572`; 72,000 unique first-column SGB IDs |
| HRGM2 representative archive `HRGMv2_Rep_Genome.tar.gz` in record 19482781 | 3,772,220,057 bytes; publisher MD5 `f473e5ceb67b433c887bc690ad20c01e`; also calculate SHA-256 locally |
| ELGG `ELGG_representatives_2172.zip` in record 6969520 | 1,587,031,082 bytes; publisher MD5 `3eb992c0adc6aeb0dd0f61dc1512be9f`; also calculate SHA-256 locally |

Other full database files were not downloaded for this research. The reference preparation command must obtain their official checksums where available, compute SHA-256, count records, check archive contents safely, and write an immutable local release lock before use. A missing download is an operational state, not a biological absence.

### 4.3 Reference admission and taxonomy reconciliation

For every genome retain accession.version, source ID, source release, sequence hash, length, completeness/contamination methods, contig count/N50, host/body site, type-strain status, taxonomy, cluster representative, license, and provenance links. Inspect suspect or chimeric genomes with appropriate quality methods; absence of a quality field is not a fabricated passing value.

Separate the complete **taxonomy-anchor registry** from the **analytical sequence-admission manifest**. GTDB/GlobDB retain documented quality exceptions. Keep their identifiers and exception annotations to preserve taxonomy/count integrity; select qualified additional member references for analytical confirmation where needed. Report source count, admitted analytical count, and unsupported anchors separately. An excluded low-quality sequence must not make the source inventory appear to have a different published size, and retaining its taxonomy must not falsely certify its sequence quality.

Use GTDB representatives as the stable initial backbone. GTDB species assignment uses representative-specific ANI radii and alignment-fraction criteria; a universal transitive 95% connected-component graph is not equivalent. GTDB's documented methods use approximately 95% ANI and a 0.5 alignment fraction, with exceptions/radii recorded in release metadata. GlobDB's 96% between-source threshold is a different operation. Preserve both memberships. [GTDB-Tk FAQ](https://ecogenomics.github.io/GTDBTk/faq.html), [GTDB release methods](https://data.gtdb.ecogenomic.org/releases/release232/232.0/), [GlobDB methods](https://globdb.org/methods).

For supplements outside existing GTDB boundaries, implement this deterministic operational clustering policy:

1. Remove exact sequence duplicates, retaining all provenance. Test membership against existing GTDB representatives using that release's own radii/method, retaining ambiguous multi-match cases rather than forcing a new species.
2. Sort remaining qualified genomes by completeness descending, contamination ascending, contiguity descending, then accession.version/sequence SHA-256 lexically. Freeze the quality-method inputs and ANI-tool binary in the release lock.
3. Compare each genome directly with accepted supplemental representatives. The initial OpenBiota supplemental membership rule is ANI ≥95% and aligned fraction ≥0.50 in **both** directions; record directional values and validate this operational policy on near-neighbor/partial-genome truth. This is not a replacement definition of every source catalog's species.
4. If one representative qualifies, join it. If more than one qualifies, preserve ambiguous/overlap membership for read-level adjudication; do not count all matching representatives as detected. If none qualifies and sequence quality is sufficient, create a supplemental representative; otherwise retain unresolved status. Never use transitive chaining through member genomes.
5. Mint stable `ob:cluster:` identifiers from the initial representative sequence hash plus policy version, and preserve them in the identity registry across representative replacement. Record splits/merges through versioned edges, not silent ID reuse. Label these as operational species-level clusters, retaining original GlobDB/UHGG/other memberships separately.

Publish how this reconciliation changes the raw GlobDB representative count. The canonical organism universe may be smaller than 346,233 even though all source records remain searchable/provenanced.

Implement typed mapping edges: `same_accession`, `exact_sequence`, `member_of_cluster`, `representative_classified_as`, `nomenclatural_synonym`, `renamed`, `split_into`, `merged_into`, `overlaps`, `ambiguous`, `unresolved`. Only justified equivalence edges collapse observations. A representative's annotation is not proof that two complete cluster concepts are identical.

If a group cannot be distinguished at the read-supported rank, report one unresolved complex, not all its candidate species. Keep a stable OpenBiota concept ID and versioned external mappings; display names must never be database primary keys.

## 5. Profiling engines and exact migration

### 5.1 Required marker lane

Pin **MetaPhlAn 4.2.6 + Jan26 database** in its own environment. The official release list is newer than some changelog text. [Releases](https://github.com/biobakery/MetaPhlAn/releases), [4.2.6 parser](https://github.com/biobakery/MetaPhlAn/blob/4.2.6/metaphlan/metaphlan.py).

Database preparation and illustrative per-sample commands:

```bash
metaphlan --install \
  --db_dir /db/metaphlan/jan26 \
  --index mpa_vJan26_CHOCOPhlAnSGB_202605

metaphlan sample.R1.clean.fastq.gz,sample.R2.clean.fastq.gz \
  --input_type fastq --nproc 16 \
  --db_dir /db/metaphlan/jan26 \
  --index mpa_vJan26_CHOCOPhlAnSGB_202605 --offline \
  --mapout sample.inventory.mapout.bz2 \
  --samout sample.inventory.sam.bz2 \
  --skip_unclassified_estimation \
  -o sample.inventory.conditional.tsv
```

These are documented interfaces, not executed analyses of the user's samples. Thread counts are examples. Version 4.2 uses `--db_dir` and `--mapout`; do not copy old `--bowtie2db`/`--bowtie2out` invocations into this environment. MP3 compatibility is a separate installation. Preserve a separate native unclassified estimate if useful, but never combine its denominator with the conditional-composition output.

New strain consensus and StrainPhlAn references must use the same Jan26 marker namespace. Keep old strain histories separately identifiable; realign/rebuild rather than relabeling old marker IDs.

Current MetaPhlAn also has optional native viral-cluster profiling through `--profile_vsc`. Evaluate it as an additional viral evidence lane with its exact compatible database; do not reuse the old `--add_viruses` interface or repeat the blanket historical claim that current MP4 cannot profile viruses. It does not replace the broader kingdom-specific workflow or establish human-pathogen status.

### 5.2 Required whole-genome lane

Pin **Sylph 1.0.0**. Use c200 sensitivity as the starting configuration and validate its new two-stage `.syl2db` behavior against the ordinary `.syldb` engine before promotion. The official GTDB r232 profile sketch covers 199,923 species representatives. [Databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/), [quick start](https://sylph-docs.github.io/install%2Bquickstart/), [source/changelog](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md).

```bash
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r232-c200-dbv2.syl2db

sylph sketch -1 sample.R1.clean.fastq.gz -2 sample.R2.clean.fastq.gz \
  -S sample -c 200 -t 16 -d sketches

sylph profile -d gtdb-r232-c200-dbv2.syl2db \
  sketches/sample.paired.sylsp --minimum-ani 95 -u -t 16 \
  -o sample.sylph.tsv

sylph-tax taxprof sample.sylph.tsv -t GTDB_r232 -o sample_
```

Download once in reference preparation, capture the resolved URL and integrity, then run offline. Pin the sylph-tax package and taxonomy files as well. Where HTTP is the official endpoint, prefer an official secure mirror if available; verify the published checksum or an independently authenticated artifact before redistribution.

Preserve every original output column, including taxonomic and sequence abundance, corrected/naive ANI, coverage, containment counts, k-mer depth, lambda state, and interval fields. `--minimum-ani 95` is the profiling starting point; do not lower it to force more species. Current source defaults for minimum k-mers differ from older cookbook passages; archive `--help` and parser schema for the actual binary rather than silently using stale thresholds.

Sylph's ANI is a containment-based reference similarity estimate, not a reconstructed genome's ANI and not strain identity. `LOW` correction failures with unavailable intervals do not become confident species calls merely because a numerical ANI is printed. `HIGH` coverage failures have different semantics and must not be rejected automatically. [Output definitions](https://sylph-docs.github.io/Output-format/), [capabilities and limits](https://sylph-docs.github.io/sylph-can-cannot/).

With `-u`, `Sequence_abundance` is scaled by estimated captured sequence content and coverage semantics change; `Taxonomic_abundance` remains a coverage-normalized conditional estimate. Do not use `--estimate-read-counts` under the percentage schema. Sylph does not provide per-read assignments; mapping/read-classification supplies those separately.

The prebuilt all-genome UHGG sketch is explicitly **not dereplicated for profiling**. Do not use its 289,232 genomes as a species abundance universe. Similarly, do not add GTDB and GlobDB profiles to make one abundance vector.

### 5.3 Read-assignment and alternative-engine benchmark

Implement a reference-accessible read-assignment/competitive-mapping component. Benchmark at least one broad read classifier against Sylph and use the qualified engine for read accounting and candidate rescue. Alternatives are comparison tools, not independent votes when they share the same evidence/reference.

| Tool | Verified version and role | Required safeguards/source |
|---|---|---|
| Kraken2 + Bracken | 2.17.1 + 3.1; broad fragment classification plus abundance estimation | Kraken confidence is a k-mer support fraction, not posterior probability. Bracken uses Kraken assignments and is not independent confirmation. Match Bracken read-length distributions to the database and library. [Kraken manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown), [releases](https://github.com/DerrickWood/kraken2/releases), [Bracken](https://github.com/jenniferlu717/Bracken). |
| Struo2 | 2.3.0; reproducible database building | Do not assume old prebuilt R207 indexes are current R232. Build from the admitted manifest or verify an exact current artifact. [Repository](https://github.com/leylabmpi/Struo2). |
| ganon2 | 2.4.3; alternative whole-genome classifier and read redistribution | Preserve report type; redistribution and genome-size normalization change quantities. Calibrate default low-abundance filters. [Documentation](https://pirovc.github.io/ganon/), [releases](https://github.com/pirovc/ganon/releases). |
| mOTUs | 4.1.0; independent marker audit, 124,295 marker units in current database | 30,256 isolate-derived and 94,039 MAG-derived units; ten-marker system. Unknown-marker mass is not a species count. [Profiler](https://www.motus-tool.org/profiler/), [releases](https://github.com/motu-tool/mOTUs/releases). |
| Metabuli | 1.2.0; combined nucleotide/amino-acid rescue | Protein homology can support broader ancestry without resolving a named species. Its GTDB232 database is listed at approximately 744 GB on disk, not 744 GB RAM. [Repository](https://github.com/steineggerlab/Metabuli), [databases](https://jaebeom-kim.github.io/metabuli-doc/databases/new-database/). |

Lock classifier thresholds using calibration data and then evaluate blind held-out data. Do not hardcode unvalidated Kraken confidence, Bracken count thresholds, mapping identity cutoffs, or a mandatory two-engine agreement rule as universal truth.

## 6. Candidate discovery and read-level adjudication

### 6.1 Candidate union, not automatic result union

Generate candidates from marker detections, whole-genome profiles, targeted pathogen modules, broad read classification, and optional assemblies. A single well-supported whole-genome call may be accepted even if the marker lane lacks that organism. Conversely, three tools calling the same shared sequence do not provide three independent confirmations.

Track which evidence is genuinely distinct: unique genome regions, independent marker families, assembly linkage, read-pair linkage, and different reference coverage. Same reads/same reference/same algorithm family are correlated.

### 6.2 Competitive mapping procedure

For every newly added, discordant, low-abundance, medically material, or taxonomically ambiguous candidate:

1. Retrieve its admitted representative(s), plausible conspecific members, nearest species/complex neighbors, and relevant decoys. Include host/food/background sequences where they can plausibly explain the signal. Do not remap reads only to the desired organism.
2. Use the same audited nonhost paired FASTQs. Scan the complete dataset or a demonstrably lossless candidate-read index; do not let a primary classifier's missed assignment prevent rescue.
3. Use an aligner/configuration validated for short metagenomic reads and close-relative competition. Record supplementary/secondary alignments, mismatch and identity distributions, pairing, edit distances, informative bases, and mapping scores.
4. Mask or separately label rRNA, repeats, low-complexity sequence, plasmids, conserved mobile elements, and other regions that cannot uniquely identify the chromosomal host. Preserve mobile-element findings separately.
5. Count independent fragments and taxon-discriminating bases, not mates twice, optical duplicates, identical alignment records, or multi-mappers as separate support.
6. Measure coverage breadth, observed-versus-expected breadth, window occupancy, coverage evenness, unique informative support, and competition with near neighbors. Require distributed evidence rather than one conserved island.
7. Compare with calibrated positives, absent close relatives, extraction blanks if available, and plausible background. If blanks do not exist, say contamination was not empirically assessed; do not fabricate a negative control.
8. Produce a structured reason for acceptance, provisional status, rank downgrade, or rejection. Keep rejected evidence in the audit output, not in the main species count.

**Development starting rule, not a validated assay limit:** a new trace whole-genome species candidate should initially require at least 20 nonduplicate informative paired fragments across at least three separated informative genomic regions, with acceptable near-neighbor discrimination. Benchmark alternative thresholds by organism/genome class, then freeze validated settings. This initial rule must not automatically reject an independently validated marker or targeted assay with a different support model.

Twenty fragments can still be contamination; three conserved regions can still be ambiguous. Do not promote a candidate using these counts alone. A MAPQ threshold alone is also insufficient when the reference panel is incomplete.

### 6.3 Detection states

| State | Meaning and report treatment |
|---|---|
| `supported` | Meets the locked validated evidence rule at the stated resolution; included in the appropriate supported count. |
| `provisional` | Credible signal but insufficient validation/discrimination; visible in a clearly separate section, excluded from headline confirmed richness. |
| `ambiguous_complex` | Signal supports one unresolved group, not each candidate member; count the complex separately from resolved species. |
| `not_supported` | Candidate failed adjudication; retain machine-readable evidence and reason. |
| `not_detected` | Assessed by a specified method; no supported call. Include assay scope and limit information where validated. |
| `not_assessed` | Missing compatible assay, database, depth, input, or execution; never replace with zero. |

“Supported” is analytical support, not proof of harmlessness, causality, or infection. Do not label analytical confidence as a calibrated probability unless it has actually been calibrated against held-out truth.

## 7. One coherent abundance system

### 7.1 Required quantity types

| Quantity type | Meaning | Forbidden interpretation |
|---|---|---|
| `marker_organismal_fraction` | Marker-depth normalized composition within that marker lane's modeled universe | Direct fraction of raw reads or absolute cells/g |
| `genome_coverage_fraction` | Coverage/genome-copy proxy normalized within the whole-genome lane | Measured viable cell fraction |
| `estimated_sequence_fraction` | Modeled fraction of sequence content under a specified algorithm | Exact classified-fragment count |
| `assigned_fragment_fraction` | Read-pair/fragment assignments divided by an explicit stage denominator | Direct cell abundance |
| `estimated_read_fraction` | Classifier-derived redistribution such as Bracken | Independent raw alignments or independent evidence |
| `absolute_quantity` | Only with a separately validated quantitative assay/spike-in and units | Available from relative shotgun data by default |

Fraction-type values are stored in [0,1]; rendering converts them to percentages. Absolute quantities instead require explicit units and nonnegative values without a universal upper bound. Keep method-native values and units intact as well. Every abundance record includes lane, database, estimator, input stage, denominator ID, conditional/full fraction status, and uncertainty method.

### 7.2 Primary report composition

Select one qualified primary whole-genome lane and immutable abundance universe for a report. Its native taxonomic vector is partitioned into accepted identities, unresolved groups, and unresolved/unaccepted modeled mass. A display restricted to accepted taxa may be renormalized only when explicitly labeled **“composition within accepted taxa”**, with the retained fraction shown and original values preserved.

Additional organisms supported only outside this primary universe may have `primary_abundance: null`. Show their detection evidence and a clearly labeled secondary-lane abundance, if available. Null means unavailable, not absent. Do not insert secondary percentages into the primary pie, stack, diversity calculation, or total.

Preserve alternate profiles in an inspectable comparison view. A profile difference is a measured methodological difference requiring investigation; neither profile is automatically correct, and one should not be called inflated merely because its value is higher.

### 7.3 Mass conservation and mapping

For mutually exclusive categories in a conditional primary composition:

`sum(accepted_leaf_fractions) + unresolved_fraction + other_modeled_fraction = 1`

Tolerance before rounding: 1e-6, unless a documented estimator requires a different validated tolerance. For estimated sequence fractions with modeled unknown content, store that content as an explicitly **modeled** complement; do not relabel it as directly counted unclassified reads.

Many-to-one mappings can aggregate mutually exclusive source cells. A one-to-many taxonomy split cannot duplicate its parent abundance into each child. Keep a complex-level value until read evidence resolves the split. Ancestor rows are rollups and never summed with descendants. An SGB alias and a named species representing the same accepted concept get one row.

Shannon diversity/evenness/richness must state the included domains, estimator, detection threshold, and universe. A combined inventory with some null primary abundances is not a valid composition for Shannon. Preserve legacy diversity when its reference interval uses the legacy lane; offer the new lane's diversity separately until new references qualify.

### 7.4 Fragment accounting

Write a stage ledger: raw pairs, raw singleton reads, quality-assessed pairs, quality-retained pairs, host-removed pairs, nonhost pairs used, discarded/unpaired reads, and duplicate handling. Account for overlapping mates in base/coverage estimates. Store content hashes of the exact files each lane consumed.

If a fast QC mode only assesses reads, label it `assessment_only`; do not claim the downstream input was filtered. Absence of human hits in an initial subsample does not prove all host DNA was removed or establish the sequencing provider's preprocessing history. Actual human-read screening/filtering remains a separately recorded operation.

## 8. Taxonomy, readable names, and complete organism cards

### 8.1 Mapping files are database-specific

The current MetaPhlAn utility directory provides:

| Input marker database | Verified mapping target and artifact |
|---|---|
| Jan25 `mpa_vJan25_CHOCOPhlAnSGB_202503` | GTDB r220: [exact TSV](https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB_r220.tsv) |
| Jan26 `mpa_vJan26_CHOCOPhlAnSGB_202605` | GTDB r226: [exact TSV](https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv) |
| Older Jan21/Oct22/Jun23 | Select the actual database-specific mapping in the [official directory](https://api.github.com/repos/biobakery/MetaPhlAn/contents/metaphlan/utils); do not substitute Jan25/Jan26 mappings because the installed filename looks similar. |

Pin a commit/tag and SHA-256 for these moving-branch URLs before use. A generic `SGB2GTDB_r220.tsv` URL is not a reliable contract. Jan26's official crosswalk does not target r232; a further mapping requires accession/representative evidence and must retain both edges.

Crosswalk membership may yield a named species, a stable unnamed cluster, a group, or an empty species rank. Audit the historical “43 of 44 resolved” claim as separate counts: ID hits, species assignments, named species, placeholders, ambiguous groups, unresolved. A lookup hit is not proof that an unknown organism has been scientifically named.

### 8.2 Display policy

Prefer a supported accepted name; show a familiar synonym when it aids recognition and does not collapse distinct concepts. Otherwise show a stable placeholder such as `Gemmiger sp937890665` with “unnamed species cluster.” If only higher rank is supported, use “unresolved species within [rank]” and the stable internal ID.

Preserve meaningful suffixes such as `_A`/`_B`, `spNNN` codes, and alternative cluster identifiers. Explain them in a concise tooltip. Never silently strip a suffix, invent a species epithet, or equate a GTDB placeholder with a clinically characterized species.

### 8.3 Four independent interpretation dimensions

Every card has:

1. **Identity confidence:** supported taxon, unknown cluster, complex, or higher rank.
2. **Biological knowledge:** measured traits, strain-limited traits, predicted traits, and unknowns.
3. **Health evidence:** context-dependent beneficial functions, opportunism/established pathogenicity, observational disease associations, or insufficient characterization.
4. **Sample findings:** abundance, compatible reference position, direct gene/pathotype evidence, strain support, and limitations.

Do not derive a single good/bad flag from a case-control association. If the interface requires a summary badge, allow `commensal_context`, `potential_concern_context`, `mixed_or_context_dependent`, and `insufficient_characterization`, with supporting claims and an explicit separation from detection confidence. A species-wide badge must not overrule a resolved pathogenic strain or specific toxin finding.

### 8.4 Interpretation and trait sources

| Resource | Implementation use | Retrieval/terms and limits |
|---|---|---|
| NCBI Taxonomy | Original/current taxids, name classes, ranks, merged/deleted IDs | [Taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz), [readme](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt), [policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). Pin date/hash. Government information policy does not make every submitted or linked record CC0. |
| GTDB metadata | Genome-based taxonomy, representative membership, genome quality and type-material links | [r232 files](https://data.gtdb.ecogenomic.org/releases/release232/232.0/). CC BY-SA 4.0. Accessions, not display-name joins. |
| LPSN | Nomenclature, basonyms, valid publication, type strains, taxonomic opinions | [API](https://api.lpsn.dsmz.de/), [API fields](https://lpsn.dsmz.de/text/lpsn-api), [terms](https://lpsn.dsmz.de/text/copyright). Free API registration; automated access via API/download, not page scraping. CC BY-SA 4.0; attribute/link taxon records. |
| BacDive v2 | Culture/strain physiology, oxygen tolerance, Gram reaction, substrates, products, isolation metadata | [API](https://api.bacdive.dsmz.de/), [terms](https://bacdive.dsmz.de/about). Registration removed February 2026; v1 frozen April 2025. CC BY 4.0 with a commercial-contact request also stated by the provider. Keep strain specificity and source references. |
| BugSigDB | Structured study/experiment/signature associations with direction and context | [Primary paper](https://www.nature.com/articles/s41587-023-01872-y), [export guide](https://bugsigdb.org/Help:Export), [license](https://bugsigdb.org/Project:About#Data_licensing). ODC-BY 1.0. Snapshot studies/experiments/signatures CSV; direction is not health valence. |
| GMrepo v3 | Cross-study prevalence/phenotype context and association discovery | [2026 paper](https://academic.oup.com/nar/article/54/D1/D734/8340991), [site](https://gmrepo.humangut.info/), [API examples](https://github.com/evolgeniusteam/GMrepoProgrammableAccess). 118,965 runs/samples in the paper, not independent healthy participants. Original data terms CC BY-NC 3.0; broader current commercial rights unverified. |
| Disbiome | Curated disease/organism/comparison/publication records | [Export](https://disbiome.ugent.be/export), [experiment JSON](https://disbiome.ugent.be:8080/experiment), [paper](https://pubmed.ncbi.nlm.nih.gov/29866037/). Noncommercial/personal use; commercial permission required. |
| gutMDisorder v3, 2026-06-10 | Human/mouse, gut/oral, phenotype/intervention and recomputed-data associations | [Resource](https://bio-computing.hrbmu.edu.cn/gutMDisorder/resource.dhtml), [API](https://bio-computing.hrbmu.edu.cn/gutMDisorder_api/api/resource/table-data), [literature download](https://bio-computing.hrbmu.edu.cn/gutMDisorder_api/api/resource/download?fileName=gutMDisorder_v3_Literature-based_All.xlsx). Current bulk schema/redistribution license require verification; do not conflate its v2 paper with a v3 methods paper. |
| MiMeDB 2.0 | Microbe/metabolite/reaction/trait discovery and literature links | [About](https://mimedb.org/about), [downloads](https://mimedb.org/downloads), [2026 paper](https://pubmed.ncbi.nlm.nih.gov/41273085/). Bulk release dates 2025-10-01; CC BY-NC 4.0. Model/inference edges are not measured sample metabolite concentrations. |
| Madin et al. trait synthesis | Broader trait coverage with source-level provenance | [Repository](https://github.com/bacteria-archaea-traits/bacteria-archaea-traits), [paper](https://doi.org/10.1038/s41597-020-0497-4), [archive](https://doi.org/10.6084/m9.figshare.c.4843290). Combines 26 resources; deduplicate original evidence and check inherited data rights. |
| ProTraits | Explicitly predicted bacterial/archaeal traits | [Resource](http://protraits.irb.hr/), [paper](https://pubmed.ncbi.nlm.nih.gov/27915291/). 424 traits/3,046 species; inference, not direct phenotyping. Separate data grant unverified; optional. |
| Rhea + UniProt + ChEBI | Reaction/chemical vocabulary and curated protein links | [Rhea downloads](https://www.rhea-db.org/help/download), [files](https://ftp.expasy.org/databases/rhea/), [terms](https://www.rhea-db.org/help/license-disclaimer). Rhea CC BY 4.0. Require a gene/protein link before describing sample metabolic capacity. |
| VMH / AGORA2 | Optional strain-specific metabolic models and constrained simulation | [VMH paper](https://pubmed.ncbi.nlm.nih.gov/30371894/), [AGORA2 paper](https://pubmed.ncbi.nlm.nih.gov/36658342/), [API](https://www.vmh.life/_api). Verify exact model release/access/license. Model flux is not measured metabolism. |
| HMP / MetaHIT | Eligible prevalence cohorts after uniform reprocessing | [NIH HMP](https://commonfund.nih.gov/hmp), [PRJNA43017](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA43017), [MetaHIT ERP000108](https://www.ebi.ac.uk/ena/browser/view/ERP000108). Use accession-based sources; do not rely on an unverified legacy domain. |
| curatedMetagenomicData | Existing standardized cohort data and metadata | [Documentation](https://waldronlab.io/curatedMetagenomicData/). Existing MP3/HUMAnN3 profiles are not automatically compatible with the new inventory. Preserve source cohort conditions and participant grouping. |

Public BugSigDB snapshot endpoints to pin:

```text
https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/studies.csv
https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/experiments.csv
https://bugsigdb-csv.s3.us-east-va.perf.cloud.ovh.us/signatures.csv
```

Use ECO for evidence type, OBI for assays, UBERON/ENVO for body site/material, MONDO for disease, and ChEBI for chemicals. Ontology terms standardize vocabulary; they are not evidence that an organism has a trait. Sources: [ECO](https://obofoundry.org/ontology/eco.html), [OBI](https://obofoundry.org/ontology/obi.html), [UBERON](https://obofoundry.org/ontology/uberon.html), [ENVO](https://obofoundry.org/ontology/envo.html), [MONDO](https://obofoundry.org/ontology/mondo.html), [ChEBI](https://obofoundry.org/ontology/chebi.html).

### 8.5 Build 600 useful cards without inventing claims

Prioritize the union of actual project organisms, prevalent taxa in diverse compatible stool cohorts, existing pathogen/strain targets, and newly admitted gut species. Keep this priority set versioned. Include both named organisms and stable unnamed species clusters.

Each card must state a supported identity, observed detection evidence, the meaning of its abundance, and what is known/unknown. Add primary-cited traits and disease/intervention evidence only when supported. For a poorly characterized cluster, a concise honest card is preferable to copied genus-level claims.

Count coverage separately: identity coverage; card coverage; direct trait evidence; predicted trait evidence; human association evidence; intervention evidence; compatible percentile availability. Publishing 600 empty templates does not meet the curated-priority-registry requirement, but a reviewed record concluding that biological characterization is unavailable is valid and must cite the searched source snapshots.

Deduplicate evidence by original paper, cohort, experiment, population, comparison, and outcome. The same paper appearing in three aggregators is one study. Track direction contradictions, covariate adjustment, medication, geography, age, body site, 16S versus shotgun, and taxonomic resolution. Do not inherit species-specific findings across an unresolved genus.

Intervention evidence remains reachable from these cards, including mechanistic, in-vitro, animal, observational, and human-intervention studies under their correct evidence labels. This increment must not silently remove preclinical options or upgrade them into demonstrated human benefit. A detection or association alone cannot establish a treatment response.

### 8.6 Healthy-reference prevalence

A “common in healthy participants” statement requires eligible independent-participant denominator, age/body-site population, health definition, assay, detection threshold, depth, profiler/database, repeated-sample handling, geography, and uncertainty interval. Report the actual compatible cohort size, not all records in an aggregator.

Pediatric samples cannot inherit adult normal ranges by default. Unknown eligibility produces “reference not compatible,” not a misleading low/high percentile. An organism absent from a small stool aliquot is “not detected,” not necessarily biologically missing. Additional detected species do not by themselves establish a healthier ecosystem or greater donor suitability.

## 9. Beyond bacteria: separate, rigorous organism lanes

The bacterial/archaeal benchmark must not be padded with phages, fungi, parasites, strains, or genes. Add kingdom-specific results with their own count units and validation.

| Resource | Verified release/count | Required integration |
|---|---|---|
| **UHGV** | v1.0, 2025-10-24; **873,995 viral genome records / 168,536 species-level vOTUs** | Primary human-gut viral reference. [Deposit](https://zenodo.org/records/17402089), [maintainer](https://github.com/snayfach/UHGV). CC BY 4.0. Select quality tiers, retain uncertainty/completeness and cluster membership. |
| **UHGV-classifier** | Assembly classification workflow | [Repository](https://github.com/snayfach/UHGV-classifier). Classifies viral contigs; not a read profiler. Retain its separate software license. Use validated read mapping for read-level support. |
| **Cultivated Gut Fungi (CGF)** | Published 2024 catalog: **760 genomes / 206 species / 48 families**, including 69 previously unidentified species | [Paper](https://doi.org/10.1016/j.cell.2024.04.043), [repository](https://github.com/yexianingyue/Cultivated-Gut-Fungi). Repository later expanded NCBI content without a verified replacement count; pin accessions and terms instead of guessing. |
| **EukDetect2** | Database deposited 2026-03-16; **6,948 reference entries**, not species | [Code](https://github.com/allind/EukDetect), [database](https://zenodo.org/records/19056625). MIT code/CC BY 4.0 database. New database incompatible with EukDetect v1. Orthogonal eukaryote markers plus targeted whole-genome confirmation. |
| **WormBase ParaSite** | WBPS19, March 2024; **274 genomes / 208 species** | [Official resource](https://parasite.wormbase.org/). Includes free-living and nonhuman parasites; select relevant references. Development is on hold; retain release/access metadata. |
| **NCBI viral/eukaryotic references** | Accession-versioned, date-frozen subset | [Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/). Fill fungal/protist/helminth and human-pathogen gaps with appropriate close relatives and quality checks. |
| **MetaVR** | 2026 paper: **24,435,662 uncultivated viral genome records / 12,705,385 vOTUs**, all environments | Successor to IMG/VR, includes UHGV. Optional broad discovery; exact deployable release/access/terms must be verified. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807716/), [site](https://www.meta-virome.org/). Do not add its raw count to UHGV. |
| **MicroEuk v3** | 79,920,431 proteins in MicroEuk100; 51,767,730 in MicroEuk90; 29,898,853 in MicroEuk50 | [Deposit](https://zenodo.org/records/10139451), CC BY 4.0. Optional protein/assembly rescue. Proteins are not organism detections or whole-genome references. |
| **FungiDB / VEuPathDB** | Current release/count not verified in this audit | [FungiDB](https://fungidb.org/fungidb/), [VEuPathDB](https://veupathdb.org/veupathdb/). Optional annotations once access/version/terms are verified; use accession-pinned NCBI fallback for the operational core. |

UHGV already integrates MGV, GPD, CHVD, IMG/VR and other sources. Keep them for provenance/reproduction, not naive concatenation: [MGV](https://portal.nersc.gov/MGV/) has 189,680 draft genomes/54,118 vOTUs; [GPD](https://www.sanger.ac.uk/data/gut-phage-database/) lists 142,809 nonredundant phage genomes; [CHVD v1.1](https://zenodo.org/records/4498884) has 45,033 representative vOTUs across several human body sites. These count units and scopes differ.

For viruses retain predicted host, host-assignment method/confidence, taxonomy release, and unknown host. A phage is not automatically a human pathogen; phage abundance is not bacterial-host abundance. DNA-only libraries cannot establish absence of RNA-only viruses. Separate viral cluster detection, named taxon assignment, prophage sequence, and active infection.

For fungi/protists/helminths, validate extraction/library sensitivity, human/food cross-mapping, near-neighbor discrimination, and coverage distribution. A larger reference database cannot repair absent DNA extraction or a sample containing no shed material. Keep “not assessed” distinct from “not detected.”

Preserve and independently regression-test the existing pathogen catalog and strain-level workflow. No new broad inventory lane is a substitute for determinant-specific pathogenicity/AMR interpretation or resolved carrier attribution.

## 10. Unknown sequences, functional data, and optional assembly

### 10.1 Report distinct unknown categories

Separate: known unnamed reference cluster; known organism below naming resolution; unresolved close-relative complex; modeled unrepresented sequence content; directly unassigned fragments; candidate novel assembly; host/food/nonmicrobial sequence; and low-quality/unusable sequence. “Unknown” is not one biological class.

Explain whether unclassified sequence estimates arise from a marker model, a genome sketch model, or actual read accounting. Marker mapping covers selected genomic regions and its unmapped fraction is not the fraction of unknown organisms.

### 10.2 Functional catalogs are useful but not species counters

Use [GMGC](https://gmgc.embl.de/download.cgi), UHGP/IGC and existing curated gene panels to explain otherwise unannotated functions. GMGC v1's approximately 302.7 million 95%-identity unigenes are not 302.7 million organisms. Shared genes and horizontally transferred sequences cannot uniquely assign a host without genomic/read linkage.

Preserve species, strain, pathotype/serotype where validated, and gene-level calls in separate fields. Do not infer every strain's accessory genome from a species representative. AMR or toxin reads can be reported with unresolved carrier rather than assigned to whichever abundant species is nearby.

### 10.3 Deep-discovery assembly

Use a reproducible assembler/binning workflow consistent with project conventions, then CheckM2/GUNC or qualified alternatives, independent coverage support, taxonomic assignment, and chimera review. Record tool versions and model databases. Sources: [CheckM2](https://pubmed.ncbi.nlm.nih.gov/37500759/), [GUNC](https://doi.org/10.1186/s13059-021-02393-0), [MIMAG standards](https://doi.org/10.1038/nbt.3893).

A proposed new reportable MAG should initially require ≥70% estimated completeness, ≤5% contamination, distributed read support, acceptable chimera checks, and stable cluster assignment; validate this policy on truth sets. This is an engineering admission threshold, not the published MIMAG definition. Published MIMAG high quality requires **>90% completeness, <5% contamination, 23S/16S/5S rRNA genes, and tRNAs for at least 18 amino acids**. Completeness/contamination estimates alone are insufficient. In simulated communities, require admitted bins to demonstrate ≥95% truth-based purity and ≥70% truth-based genome completeness; CheckM2/GUNC predictions must not serve as their own ground truth.

Short contigs, multiple bins from one genome, phages/plasmids, and contaminants must not inflate species count. A novel MAG belongs in a research discovery section until its identity and specificity qualify. Assembly failure at low coverage is not evidence that an organism is absent.

## 11. Machine-readable contracts

Implement strict versioned schemas in the repository. The following fields and semantics are normative; serialization may follow existing project conventions. Files generated by the implementation are normal runtime/code artifacts, not additional handoff documents required by this specification.

### 11.1 Reference lock

```yaml
schema_version: openbiota.reference_lock.v1
release_id: openbiota-v08.1-YYYYMMDD
research_cutoff: '2026-09-19'
artifacts:
  - artifact_id: metaphlan-jan26-species-metadata
    source_release: mpa_vJan26_CHOCOPhlAnSGB_202605
    url: https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan26_CHOCOPhlAnSGB_202605_species.txt.bz2
    sha256: 441f86af40f587b969b40ebdb79ffb02b311ffec4c81024d8001e529c6ea5572
    record_type: species_genome_bin
    observed_records: 72000
    sequence_type: metadata
required_artifact_fields:
  - source_release
  - resolved_url
  - retrieved_at
  - sha256
  - bytes
  - license_id_or_review_status
  - count_unit
  - admitted_records
  - source_provenance
  - quality_policy_id
  - taxonomy_namespace
required_tool_fields:
  - version
  - container_digest_or_binary_sha256
  - dependency_lock_hash
  - command_schema_hash
```

The YAML shows an integrity anchor and schema requirements, not a complete production lock. The preparation step must populate every actual artifact/tool field; unresolved placeholders are forbidden in a released lock. Record license-unverified optional sources as unavailable with a reason, not falsely admitted. Core reference preparation must continue using available qualified sources.

### 11.2 Observation and accepted result

```json
{
  "schema_version": "openbiota.organism_observation.v1",
  "sample_id": "example",
  "observation_id": "content-addressed-id",
  "lane_id": "inventory_whole_genome",
  "input_manifest_id": "sha256-reference",
  "reference_lock_id": "release-reference",
  "native_taxon_id": "source-release:cluster-id",
  "canonical_concept_id": "ob:taxon:stable-id",
  "supported_rank": "species_cluster",
  "status": "supported",
  "detection_rule_id": "validated-rule-id",
  "taxonomic_resolution": "reference_cluster",
  "native_values": {},
  "evidence": {
    "informative_fragments": null,
    "informative_regions": null,
    "genome_breadth": null,
    "coverage": null,
    "marker_count": null,
    "near_neighbor_resolved": null,
    "negative_control_status": "not_available"
  },
  "abundances": [],
  "mapping_edge_ids": [],
  "reason_codes": [],
  "provenance_record_ids": []
}
```

Null values above are schema illustrations, not permission for a production supported call with no evidence. Validate the evidence requirements of its stated rule. Reject NaN/Infinity and distinguish missing, censored, below-limit, and measured zero.

An abundance entry requires `value`, `quantity_type`, `estimator`, `lane_id`, `denominator_id`, `conditional`, `interval`, `interval_method`, `native_unit`, and `reference_lock_id`. Missing intervals are explicit. Do not fabricate bootstrap uncertainty by treating correlated reads, genomes, or samples as independent biological replicates.

### 11.3 Identity and interpretation records

An identity record requires stable concept ID, rank, kingdom/domain, source cluster memberships, accession.version members, display name, nomenclatural names/synonyms, naming status, mapping edges, representative, release lineage, and merge/split history.

An evidence claim requires:

```yaml
claim_required_fields:
  - claim_id
  - subject_concept_id
  - subject_rank
  - predicate
  - object
  - evidence_kind
  - primary_source_url
  - source_record_locator
  - original_study_id
  - cohort_or_experiment_id
  - host_species
  - body_site
  - assay
  - direction_or_effect
  - population_and_context
  - measured_inferred_or_predicted
  - strain_scope
  - contradiction_group
  - extraction_method
  - reviewer_status
  - license_provenance
allowed_evidence_kinds:
  - cultured_trait
  - genome_inference
  - computational_prediction
  - in_vitro_experiment
  - animal_experiment
  - human_observational
  - human_intervention
  - systematic_synthesis
  - nomenclature_or_identity
```

Fields can contain explicitly documented `unknown`/`not_applicable`; never guess a study's effect size, dose, population, or causal design. Primary source links must support the actual statement. Aggregator URLs alone can support a database fact, not a stronger unverified biomedical claim.

### 11.4 Output summary

The summary contains separate counts for supported named species, supported unnamed species clusters, unresolved complexes, provisional taxa, supported strains, other kingdoms, and failed/not-assessed lanes. Also include composition-universe ID, coverage-card statistics, reference-compatible percentile count, catalog-admission statistics, and validation-release ID.

For a combined species/cluster headline, define `supported_species_level_concepts = supported_named_species + supported_unnamed_species_clusters`, after deduplication. Unresolved complexes are displayed separately and not represented as a known number of species. Include all kingdom-specific totals without adding them to the bacterial/archaeal benchmark result.

## 12. CLI, modules, caching, and failure handling

Implement equivalent commands within the existing CLI; the following command names define required behavior:

```bash
openbiota references prepare --preset v08.1 --output /db/openbiota/v08.1
openbiota references audit --lock /db/openbiota/v08.1/reference-lock.json

openbiota inventory run --sample SAMPLE2_A02 \
  --r1 SAMPLE2_A02.R1.fastq.gz --r2 SAMPLE2_A02.R2.fastq.gz \
  --mode comprehensive --reference-lock /db/openbiota/v08.1/reference-lock.json \
  --reuse-compatible-cache --output results/SAMPLE2_A02

openbiota inventory reconcile --sample-results results/SAMPLE2_A02
openbiota evidence build --preset v08.1 --minimum-reviewed-concepts 600
openbiota report render --sample-results results/SAMPLE2_A02 --complete-inventory

openbiota benchmark run --suite v08.1-release --split heldout
openbiota migration compare --sample SAMPLE2_A02 --from legacy --to v08.1
```

`--preset v08.1` is an embedded registry implemented from this document's source tables, not a missing supplemental file. Runtime locks, models, catalogs, caches, JSON and benchmark artifacts are generated by the code.

Proposed module boundaries: reference registry/admission; preprocessing ledger; profiler adapters; candidate generator; competitive adjudicator; taxonomy graph; abundance model; interpretation evidence store; report inventory; benchmark/migration runner. Keep heavy reference downloads separate from sample execution. Permit resume and retry after partial failure without re-running valid expensive work.

Exit status must distinguish successful complete run, successful explicitly partial run, validation failure, unavailable reference, incompatible cache, and invalid input. Do not silently downgrade comprehensive to standard or render a missing lane as zero organisms. Do not discard a whole report because an optional enrichment service is offline; display source availability and cached evidence freshness.

### Resource planning

The official Jan26 metadata archive is approximately 5.6 GB and its Bowtie2 index archive approximately 39 GB; these are archive sizes, not RAM requirements. Sylph r232 c200 `.syl2db` is approximately 19 GB and GlobDB's converted sketch approximately 30 GB. Full reference genomes, mapping indexes, SAMs and assembly workspace add substantial disk requirements. Author-reported Sylph memory examples do not bound the entire pipeline.

Implement disk-space preflight, resumable downloads, checksum verification, bounded parallelism, reference-shard scheduling, and measured CPU/RAM/disk telemetry. Accept user-provided resource limits. Benchmark wall time and peak memory on the target environment before setting defaults. Do not claim exact costs or runtimes from this paper review.

## 13. Report specification

### 13.1 Front-loaded overview

Show a readable “Organisms detected” panel before long explanations:

- Supported bacterial/archaeal species-level concepts, divided into named and unnamed clusters.
- Unresolved complexes and provisional findings, separately.
- Fungi, protists/parasites and viral clusters, using their correct units.
- Number with compatible reference percentiles and number without them.
- Identity/card/evidence coverage and the exact profiling release.
- A compact sensitivity/depth statement and a link to the full inventory.

The headline is evidence-based richness, not a health score. Do not color a higher count automatically green or a lower count automatically red. Biofilm/disease/FMT modules retain their own meaning and validation.

### 13.2 Full inventory

Render every accepted organism/cluster without top-N truncation. Provide a stable table with readable name, aliases/ID, domain, detection status, primary abundance or unavailable state, reference percentile if compatible, interpretation badge, and concise evidence summary. Use sorting, filters and search in HTML; paginate cleanly in PDF with repeated headers and complete continuation counts.

Make secondary-lane values inspectable, clearly separated from primary composition. Very small trace values should use sensible significant digits or scientific notation; do not round a supported call to an unexplained 0.00%. Group provisional/complex findings in their own visible section. Retain an export of all observations and reasons.

### 13.3 Organism detail card

Each card answers: What is it? How confidently was it detected? How much was estimated and by which method? What is known about this specific organism/strain? What relevant genes/functions were measured in this sample? Is a compatible reference comparison available? What evidence or intervention research relates to this finding, and how directly?

Show study links near the claims. Generated HTML external links must use `target="_blank" rel="noopener noreferrer"`. Clearly identify observational, mechanistic, animal, and human-intervention evidence. Do not hide limited evidence; do not promote it beyond its design.

### 13.4 Patient-facing language

Use: “Detected with supporting sequence evidence,” “unnamed species cluster,” “reference comparison unavailable,” “estimated composition within [universe],” and “not detected under this assay.” Avoid “all organisms present,” “complete ecosystem restored,” “every database searched,” and definitive organism absence from a negative shotgun result.

A donor pair covering more catalogued species does not establish that all recipient needs will be restored, that strains will engraft, or that disease transmission risk is eliminated. This increment supplies better measured inputs to existing donor analyses; it does not convert detection breadth into a clinical clearance rule.

## 14. Validation: real sensitivity gains with controlled false positives

### 14.1 Split design and truth

Use separate development/calibration and locked held-out test sets. Split by community, source genome/strain, and where possible study. Do not fit thresholds on SAMPLE2 or the five project samples. Benchmark database-present exact representatives as a smoke test only; the decisive tests include held-out strains, close relatives, mixed strains, reference-absent species, low abundance, contamination, and extraction biases.

Audit held-out membership against **every** admitted catalog. Removing a genome from GTDB while retaining it in GlobDB/UHGG/HRGM2 is not a valid open-set test. Freeze taxonomy/cluster definitions before truth generation. Evaluate named species, cluster-level, complex-level and higher-rank outputs separately.

### 14.2 Public benchmark sources and executable accessions

| Benchmark | Use and verified access |
|---|---|
| CAMI II | Community/strain complexity and standard profiling evaluation. [Paper](https://doi.org/10.1038/s41592-022-01431-4). Strain Madness [data DOI](https://doi.org/10.4126/FRL01-006425521), [truth](https://zenodo.org/records/5006866); Toy HMP [data](https://doi.org/10.4126/FRL01-006425518). |
| CAMI III public Toy Human Gut | 20 paired-end 150-bp samples, about 100 Gbp total; use corrected sample/subject mapping announced August 2026. [Public dataset](https://cami-challenge.org/datasets/toy-human-gut/). Full challenge data have separate terms and are not a redistributable default CI fixture. |
| OPAL | Rank-specific profiling precision, recall, abundance distance and related metrics. [Paper](https://doi.org/10.1186/s13059-019-1646-y). Preserve native truth and canonical mapping audit. |
| LEMMI v2 | Reproducible benchmarking framework, not proof of this pipeline's accuracy. [2026 paper](https://doi.org/10.1186/s13059-026-04089-9), [v2.2.0 deposit](https://zenodo.org/records/18293506). |
| MBARC-26 physical mock | 23 bacterial and 3 archaeal strains; **SRR3656745**, PRJNA324704, approximately 173,981,994 paired spots before QC. Supports 8M/30M pair subsets if sufficient retained reads. [Study](https://www.nature.com/articles/sdata201681), [run](https://www.ebi.ac.uk/ena/browser/view/SRR3656745). |
| Zymo D6300 physical mock | **SRR12324253**, PRJNA648136, approximately 51,466,358 paired spots. [Run](https://www.ebi.ac.uk/ena/browser/view/SRR12324253). Validate exact lot/composition and library metadata. |
| Zymo logarithmic mock | **ERR2935805**, PRJEB29504, approximately 47,832,553 paired spots, roughly 101-bp mates. [Run](https://www.ebi.ac.uk/ena/browser/view/ERR2935805), [study](https://doi.org/10.1093/gigascience/giz043). |
| Zymo even mock | **ERR2984773**, PRJEB29504, approximately 8,785,731 paired spots before QC. [Run](https://www.ebi.ac.uk/ena/browser/view/ERR2984773). May fall below 8M after QC; never duplicate reads to manufacture depth. |
| ATCC MSA-1003 physical mock | **SRR8359173**, PRJNA510527, approximately 5,019,157 paired spots, 125-bp mates. [Run](https://www.ebi.ac.uk/ena/browser/view/SRR8359173), [benchmark context](https://doi.org/10.1186/s12859-022-05103-0). Native-depth validation only; cannot become an 8M/30M observed library by resampling duplicates. |

Use ENA/SRA metadata and file checksums to resolve files at execution. Distinguish spots/read pairs from individual reads, raw from retained depth, and DNA-mix proportions from cell-mix proportions. Physical mocks test extraction/preparation effects only to the extent their actual experimental design supports.

### 14.3 Six-hundred-species capability test

Generate at least 20 independent ground-truth communities: 10 approximately even and 10 long-tailed. Each contains 600 distinct, reference-resolvable bacterial/archaeal species-level truth concepts, primarily gut-relevant with archaea included. Use source genomes of 1–8 Mb and five read/error seeds per community. At least the decisive held-out-strain tier uses different strain genomes from the profiling representatives, while remaining within validated species boundaries.

For this capability tier, every truth species has **at least 0.02% of microbial DNA mass**; the 600 minima therefore consume 12% of total microbial DNA, leaving the remaining mass for the selected distribution. Simulate paired reads using empirical read-length, overlap, fragment, error, and GC-bias distributions; do not rely only on perfect 150-bp reads. Create 8M and 30M retained nonhost-pair tiers. Preserve separately input DNA-mass/genome-copy truth and the realized generated/retained fragment- and base-origin truth. Bias can make these fractions differ. Keep the capability gate's declared input-mixture tiers and report performance additionally by realized support. Add separate host/food/contamination challenges rather than silently altering the defined microbial denominator.

Required held-out targets — engineering release criteria, not published clinical standards:

| Metric | Target |
|---|---|
| Mean species recall, separately for even and long-tail 600-species tiers | ≥95% at 8M pairs; ≥98% at 30M pairs |
| Mean supported species recovered | ≥570/600 at 8M; ≥588/600 at 30M |
| Macro species precision across communities | ≥99%; one-sided 95% community-bootstrap lower bound ≥98.5% |
| Per 600-species sample false positives | ≤5 supported incorrect species-level calls |
| Abundant taxa ≥0.1% DNA fraction | Recall ≥99% in applicable reference-present tests |
| Intermediate taxa 0.01–0.1% DNA fraction | Recall ≥90% at 8M and ≥95% at 30M in the separate gradient test |
| Rare taxa 0.001–0.01% | Report abundance/depth-specific detection curves and precision; no blanket complete-recovery claim |
| Incremental calls absent from the legacy accepted inventory | Precision ≥99% on truth-known validation; audit this subset separately from aggregate precision |
| Broad catalog versus current baseline, fixed test/precision | Positive paired recall gain with 95% confidence interval above zero; target ≥10 percentage points in catalog-gap test sets |
| Shared baseline taxa | Recall loss ≤0.5 percentage points unless a documented false-positive correction accounts for the difference |
| Physical simple mocks | Zero unadjudicated high-confidence unexpected calls; explain contaminant evidence separately rather than changing truth post hoc |
| Open-set species, absent from every admitted reference and with ≥0.05× true coverage | Incorrect supported assignment to any reference species-level concept, **named or unnamed**, ≤1% of those truth organisms; also count every false reported taxon in ordinary precision. Correct higher-rank abstention remains separate. |

Do not average away a failure in the long-tail or low-depth tier. If targets fail, identify whether reference coverage, classifier settings, adjudication, ambiguity, or library sensitivity is limiting. Improve the component and retest on a new untouched held-out set where tuning has consumed the previous set. Do not lower the precision threshold or relabel provisional results simply to hit 600.

### 14.4 Quantification and depth

Benchmark each quantity against matching truth: realized retained sequence-origin truth for sequence fraction, realized retained fragment-origin truth for read fractions, and appropriate genome-copy truth for coverage/cell proxies. Report additional bias against the original DNA/cell mixture separately, because extraction, GC and library effects intervene. Report Bray–Curtis/L1 error, bias by abundance/genome size/GC/relatedness, absolute fraction error, and intervals where implemented.

For all truth taxa with true abundance ≥0.1%, target median absolute relative abundance error `median(abs(p_hat - p) / p) ≤0.20`; a missed truth taxon has `p_hat = 0` and stays in the calculation. Target Bray–Curtis ≤0.10 over the **complete matched union** of truth and predicted concepts, including false-call mass and explicit unknown categories, with matching full or conditional denominators. Do not evaluate only successfully recovered taxa or renormalize away false mass. Evaluate physical mocks separately under their correct mixture assumptions. These are proposed engineering gates, not measured performance.

Evaluate modeled unknown sequence fraction at 0%, 10%, 30%, and 50% deliberately unrepresented microbial DNA; target mean absolute error ≤5 percentage points. If unqualified, display it as an unvalidated model estimate or withhold that numerical estimator, without suppressing actual read-accounting results.

Expected fragment support is `N_pairs × retained_fragment_fraction`. Approximate coverage is `N_pairs × retained_fragment_fraction × taxon_effective_unique_aligned_bases_per_pair / genome_length`. Using input DNA fraction in place of retained fragment fraction requires an explicit approximately unbiased, comparable-fragment assumption. The effective bases term is smaller than two full read lengths when mates overlap; it also depends on trimming/alignment. Ideal Poisson breadth `1 - exp(-coverage)` is a planning approximation, not an assay sensitivity guarantee.

Rarefy observed libraries at 1, 2, 4, 8, 15, 30 and 60 million retained pairs only where sufficient original reads exist, with five random subsets per available depth. Show accepted richness, newly confirmed taxa, false-positive estimates from controls/truth, and saturation. Do not extrapolate a real sample's species count to a precise unseen total without an explicitly validated model.

### 14.5 Project-sample migration audit

Reprocess all available project samples, including SAMPLE2 and the donor candidates, through locked lanes. Produce an internal per-sample change table: supported gains, deduplications, losses and reasons; named versus unnamed changes; primary abundance differences; old/new reference compatibility; strain namespace changes; pathway/pathogen consistency; and all downstream score changes.

Every reported additional organism needs its identity/provenance/evidence, not just a larger total. Independently inspect the highest-abundance discrepancies and the most medically material newly detected taxa. This audit estimates practical behavior, not truth-based sensitivity or clinical efficacy.

## 15. Acceptance tests

Implement these as meaningful automated fixtures/integration checks where possible, with explicit manual review artifacts for evidence and layout. IDs are stable release requirements.

| ID | Required result |
|---|---|
| REF-01 | Jan26 metadata hash and 72,000 unique SGB records match the pinned artifact. |
| REF-02 | GTDB r232 source/taxonomy representatives count 199,923; domain counts, quality exceptions and separately admitted analytical counts are recorded. |
| REF-03 | GlobDB r232 manifest records 346,233 representatives and 26 sources, without adding GTDB's count again. |
| REF-04 | Every admitted sequence has source identity, version/hash, quality policy, and license provenance. |
| REF-05 | All-genome UHGG sketches cannot be selected as an undereplicated species-profile universe. |
| REF-06 | New complementary sources emit marginal cluster/strain/quality contributions. |
| REF-07 | Missing/changed/checksum-failing files prevent that artifact's admission, not silent substitution. |
| REF-08 | Runtime execution uses immutable locks and no automatic latest-version downloads. |
| REF-09 | Optional source unavailability does not disable the qualified open core. |
| REF-10 | Disk/resource preflight distinguishes archive storage from peak RAM. |
| TAX-01 | One biological concept with MP3, MP4 and GTDB aliases renders once. |
| TAX-02 | A cluster split does not duplicate parent abundance or count all children without evidence. |
| TAX-03 | A complex remains one unresolved complex, separate from resolved species count. |
| TAX-04 | Jan26→r226 crosswalk cannot masquerade as direct Jan26→r232 equivalence. |
| TAX-05 | A mapping hit, placeholder, named species, and ambiguous assignment have separate statistics. |
| TAX-06 | Meaningful GTDB suffixes and stable placeholder IDs survive display conversion. |
| TAX-07 | Renamed/merged/deleted NCBI IDs retain historical provenance. |
| TAX-08 | Two different concepts sharing a familiar name cannot be merged by string alone. |
| TAX-09 | HumGut custom numeric identifiers cannot be treated as NCBI taxids without mapping. |
| TAX-10 | Full taxonomy mappings survive reference migration, including unresolved edges; supplemental clustering is deterministic and nontransitive. |
| DET-01 | Strong validated whole-genome support is not rejected solely because MP4 lacks a marker call. |
| DET-02 | A shared conserved island cannot support multiple species calls. |
| DET-03 | Candidate mapping includes near neighbors and relevant decoys. |
| DET-04 | Paired/duplicate/multimapped reads cannot inflate independent-fragment support. |
| DET-05 | Low-complexity, rRNA and mobile-only evidence cannot establish a chromosomal species. |
| DET-06 | Sylph LOW/NA correction output cannot be interpreted as calibrated high-confidence ANI. |
| DET-07 | Kraken plus Bracken are not counted as independent corroboration. |
| DET-08 | Unsupported/provisional candidates remain inspectable but outside supported richness. |
| DET-09 | A missing lane renders not assessed, not zero detected organisms. |
| DET-10 | Full-data rescue can recover a truth organism missed by the initial candidate-read classifier. |
| DET-11 | A novel MAG cannot duplicate an existing bin/cluster, plasmid, or fragmented genome in species counts. |
| DET-12 | Small-genome/viral validation does not inherit bacterial two-stage sketch assumptions. |
| ABU-01 | The historical 130.777% merged fixture is rejected as a single composition. |
| ABU-02 | Maximum/mean/sum merging of independently normalized profiles is prohibited by schema/API. |
| ABU-03 | Primary abundance partitions conserve mass before rendering. |
| ABU-04 | Additional secondary-only detections retain null primary abundance rather than stealing mass. |
| ABU-05 | Display normalization states the retained denominator and preserves native values. |
| ABU-06 | Parent and descendant rollups are not added together. |
| ABU-07 | Estimated sequence fraction, read fraction and organismal fraction have different quantity types. |
| ABU-08 | Zero, null, below-limit, NA and nonfinite values have distinct validated handling; absolute quantities are not constrained to [0,1]. |
| ABU-09 | Shannon/evenness only consume a coherent eligible composition with a documented universe. |
| ABU-10 | Rounded table totals and raw fractions are both checked; tiny positives remain legible. |
| QC-01 | Each profiler input hash matches the preprocessing ledger. |
| QC-02 | Assessment-only QC is never labeled applied filtering. |
| QC-03 | Reads, pairs and effective aligned bases are separately counted. |
| QC-04 | A 250,000-read host-negative subsample cannot certify whole-library host removal. |
| QC-05 | Reference/tool/parameter/input changes invalidate only affected caches. |
| QC-06 | Interrupted execution resumes without replacing valid outputs with partial files. |
| EVD-01 | Every rendered accepted/provisional/complex record has an identity card. |
| EVD-02 | At least 600 priority concepts have reviewed records; explicit knowledge gaps are documented. |
| EVD-03 | Literature-card coverage is distinct from identity/card-template coverage. |
| EVD-04 | A paper imported from three databases counts as one original study. |
| EVD-05 | Genus associations cannot become unsupported claims about all species. |
| EVD-06 | A strain-specific trait cannot become a universal species trait. |
| EVD-07 | Predicted traits and metabolic models remain labeled predictions/models. |
| EVD-08 | Association direction cannot automatically become beneficial/harmful polarity. |
| EVD-09 | Contradictory findings and relevant host/site/assay contexts survive summarization. |
| EVD-10 | Noncommercial/permission-dependent data are not silently redistributed under the app's code license. |
| EVD-11 | Preclinical intervention evidence remains accessible with correct labels, without fabricated human efficacy. |
| EVD-12 | Primary claim citations support their exact subject, rank and evidence level. |
| SCORE-01 | Legacy scoring feature vectors and outputs remain unchanged for unchanged inputs/locks. |
| SCORE-02 | New inventory abundances cannot enter old model/reference APIs without an explicit validated adapter. |
| SCORE-03 | Adult reference percentiles do not silently appear for ineligible pediatric samples. |
| SCORE-04 | Percentile-unavailable counts and report footers are derived from actual rows. |
| SCORE-05 | New richness alone cannot generate donor clearance, engraftment probability, or complete-restoration claims. |
| STR-01 | New SAM/mapping outputs are retained and identified by exact marker namespace. |
| STR-02 | Jun23 and Jan26 marker profiles cannot be mixed in one strain comparison without validated reconstruction. |
| STR-03 | Marker consensus, named strain identity, pathotype and carrier linkage are distinguishable outputs. |
| PATH-01 | Existing pathogen targets and accepted near-neighbor specificity tests remain covered. |
| PATH-02 | AMR/toxin without carrier linkage retains unresolved carrier. |
| KING-01 | Viral vOTUs, bacterial species and eukaryote reference entries use different count units. |
| KING-02 | UHGV and its integrated legacy sources are not independently added to richness. |
| KING-03 | Unknown/predicted viral host remains unknown/predicted; phage detection is not human infection. |
| KING-04 | RNA-only targets are not assessed by an unvalidated DNA-only assay. |
| KING-05 | Fungal/parasite positives pass group-specific cross-mapping/background checks. |
| BENCH-01 | Six-hundred-species tests meet both distribution tiers, depths and precision/recall gates. |
| BENCH-02 | Incremental-call precision is reported separately and meets its target. |
| BENCH-03 | Open-set truth is absent from every admitted catalog, including supplements; wrong named and unnamed species-cluster assignments are both penalized. |
| BENCH-04 | Real mock libraries are never padded with duplicate reads to reach a depth tier. |
| BENCH-05 | Quantification is evaluated against the matching DNA/genome-copy truth. |
| BENCH-06 | Unknown-fraction estimates remain marked unvalidated until their gate passes. |
| BENCH-07 | Confidence intervals resample independent communities/appropriate units, not correlated rows. |
| BENCH-08 | All project samples receive a legacy-versus-new identity/abundance/score audit. |
| UI-01 | Every inventory row appears in HTML and PDF; no hidden top-N truncation. |
| UI-02 | Headline counts, subsection totals and pagination continuations agree with JSON. |
| UI-03 | Species, unnamed clusters, complexes and provisional findings are visibly distinct. |
| UI-04 | Every external HTML link carries the specified target/rel attributes. |
| UI-05 | The report cannot render a composition above 100% without a validation error. |
| UI-06 | Complete and partial execution states are immediately visible and specific. |
| UI-07 | PDF visual inspection confirms readable names, small abundances, repeated headers and working links. |

## 16. Delivery sequence and definition of done

### Phase A — repair measurement integrity and naming

Inspect repository manifests, confirm the current baseline, locate the 130.777% composition construction, implement typed abundance/denominator contracts, deduplicate supported identity mappings, replace hardcoded counters/footers, and render readable IDs with explicit unknowns. Preserve scoring outputs. Run the relevant arithmetic, taxonomy, cache and report tests before reissuing reports.

### Phase B — qualify the new detection universe

Prepare the current marker and whole-genome references, quantify catalog overlap/marginal coverage, implement candidate adjudication, and run held-out benchmarks including the 600-species tests. Compare current Jun23, Jan26, GTDB and GlobDB behavior at matched precision. Activate the qualified comprehensive mode; keep any unqualified additions visible only as provisional/research output.

### Phase C — complete interpretation and kingdoms

Build the reviewed 600-concept priority registry, render a card for every observation, integrate source-specific terms/provenance, and qualify fungal/eukaryotic/viral extensions independently. Retain existing intervention, pathogen, function and strain integrations. Do not make restricted optional enrichment sources a dependency for basic detection or naming.

### Phase D — migrate samples and publish measured results

Reprocess all project samples, inspect every material discrepancy, validate downstream reference isolation, and regenerate complete reports. The release must state actual achieved sensitivity/precision, database versions, computation requirements, naming/evidence coverage, and sample-specific gains or losses. Publish enough locked benchmark/configuration detail for an independent open-source user to reproduce the claims.

**Done means:** the software produces a broader, evidence-adjudicated organism inventory; every reported organism is understandable; abundance arithmetic and count units are correct; new detections have measured specificity; existing calibrated analyses remain reproducible; and the 600-species capability claim is backed by the defined held-out tests. No particular human sample is required to contain a predetermined number of species.
