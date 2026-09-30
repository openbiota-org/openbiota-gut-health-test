# OpenBiota — Expanded gut bacterial recognition

**Research checked:** 2026-09-27. **Revision:** reconciled against the supplied installed inventory. **Handoff:** this single document contains the implementation requirements, reference downloads, reconciliation rules and acceptance tests. Download sequence data from the linked publishers; do not embed datasets in this document.

## 1. Required outcome and confirmed starting point

Expand supported **gut bacterial and archaeal detection** beyond the current combined inventory. Integrate into the existing multi-source abstraction, taxon registry, caches and organism report sections. Implement the requirements below; this is an additive detection improvement, not a replacement application or a new version-named directory hierarchy.

The objective is higher recall at measured high precision, including unnamed gut species and divergent strains. Catalog genomes, species clusters, strains and organisms actually supported in a sample are different counts. Do not force a sample to contain a target number of organisms.

| Existing component, supplied by the coding agent | Required treatment |
|---|---|
| MetaPhlAn 3.1.0 / `mpa_v31_CHOCOPhlAn_201901` | Retain existing calibrated scoring inputs and environment unchanged. |
| MetaPhlAn 4.1.1 / `mpa_vJun23_CHOCOPhlAnSGB_202403` | Upgrade the detection adapter to MetaPhlAn **4.2.6 / Jan26**, below; retain historical raw results. |
| Sylph 1.0.0 / `gtdb-r232-c200-dbv2.syl2db` | **Already installed.** Reuse its 199,923-cluster reference, validated outputs and caches. Adding this same database is not an expansion. |
| Jan25-SGB → GTDB R220 bridge applied alongside Jun23 profiles and R232 Sylph | Repair both source-release and destination-release mismatches; invalidate affected reconciliation caches. |
| Jun23 `sample2markers` consensus reconstruction | Existing capability is consensus-marker reconstruction, **not full comparative StrainPhlAn**. Complete the comparative workflow on eligible organisms. |
| Pathogen bundle `pathogens-dna-d12c28da90aa0f80`, catalog `143a64c812e4d77e`, 486 targets; AMRFinderPlus/StxTyper | Preserve and reconcile its organism evidence with the expanded inventory. |
| Fungal lock `myco-ref-0fba621a9babdd35`, 1,696 assemblies / 961 species; Sylph + EukDetect | Preserve. New fungal, viral or parasite catalog projects are outside this focused increment. |
| DIAMOND 2.2.6 / 462,075 proteins; rpoB community summaries; GRCh38 no-alt + PhiX174 filtering | Preserve. Protein-panel or rpoB matches are not substitutes for whole-genome species confirmation. |

Verify actual manifests at implementation start; record differences from this supplied inventory. Do not reimplement components already present and validated. Benchmark improvements against **the existing MP3 + Jun23 MP4 + GTDB R232 union**, not against MP3 alone.

### The MetaPhlAn upgrade is larger than the supplied estimate

Counts below were checked against official species metadata, rather than inferred from the installed software version.

| Marker database | SGB records | Increase over installed Jun23 |
|---|---:|---:|
| `mpa_vJun23_CHOCOPhlAnSGB_202403` | **36,822** | baseline; the supplied approximately 26k count is inaccurate for this index |
| `mpa_vJan25_CHOCOPhlAnSGB_202503` | 58,331 | 58.4% |
| **`mpa_vJan26_CHOCOPhlAnSGB_202605`** | **72,000** | **95.5%** |

Upgrade directly to Jan26. These are species-level genome bins, not 72,000 individually identifiable strains, and a 95.5% catalog increase does not imply twice as many detections in a specimen. Jan26 contains 68,266 bacterial, 3,245 archaeal and 489 eukaryotic SGB records. Its additional bacterial marker coverage is the relevant gain here.

Sources: [official database files](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/), [software releases](https://github.com/biobakery/MetaPhlAn/releases). Jan26 compressed species-file size: 1,238,616 bytes; independently checked SHA-256: `441f86af40f587b969b40ebdb79ffb02b311ffec4c81024d8001e529c6ea5572`.

## 2. Reference sources and their specific contribution

Integrate the following ten reference families in the roles specified. This does **not** require ten redundant full-sample profilers. Genome supplements must pass the overlap audit in section 3. The large catalogs cover many habitats; their sizes are not counts of gut species.

| ID | Resource / scale | Implementation role and primary downloads |
|---|---|---|
| R1 | **MetaPhlAn Jan26: 72,000 SGBs** | Required upgrade of the existing marker lane; retain alignments for strain analysis. [Release files](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/), [installation and profiling](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4). |
| R2 | **GTDB R232: 199,923 species clusters / 901,341 genomes** | Existing taxonomic backbone and baseline. Reuse the installed index; acquire member genomes for missing strain diversity and discrimination when needed. [Release](https://data.gtdb.ecogenomic.org/releases/release232/232.0/), [statistics](https://gtdb.ecogenomic.org/stats/r232). |
| R3 | **GlobDB r232: 346,233 representatives, 26 source datasets** | Add broad Sylph discovery. Includes GTDB plus **146,310 non-GTDB representatives** under its clustering rules; these are not 146,310 additional gut species. [Downloads](https://globdb.org/downloads), [methods](https://globdb.org/methods), [source manifest](https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_dataset_list.tsv). |
| R4 | **mOTUs4: approximately 124,300 species-level units**, derived from about 3.9 million genomes | Add a complementary marker-based profiling method. Its references overlap GlobDB/GTDB; its distinct detection method is useful even where sequences overlap. Use tool 4.1.0 / database 4.1. [Database](https://motus-db.org/), [code and downloads](https://github.com/motu-tool/mOTUs), [database study](https://doi.org/10.1093/nar/gkae1004). |
| R5 | **HRGM2: 155,211 near-complete genomes / 4,824 species** | Priority gut-specific supplement: species representatives plus useful divergent member genomes. Add missing clusters and better-matching references after reconciliation. [Study](https://doi.org/10.1038/s41564-025-02206-1), [download manifest](https://github.com/netbiolab/HRGM2/blob/main/DATA.md), [representative FASTAs](https://zenodo.org/records/19482781). |
| R6 | **UHGG v2.0.2: 289,232 genomes / 4,744 species** | Gut-specific residual coverage and within-species diversity. Do not assume these are absent from GTDB232. **UHGP is a gene/protein catalog, not a species-genome database.** [Metadata](https://www.ebi.ac.uk/metagenomics/api/v1/genome-catalogues/human-gut-v2-0-2), [fixed release files](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/), [foundational study](https://doi.org/10.1038/s41587-020-0603-3). |
| R7 | **ELGG: 32,277 genomes / 2,172 species** | Audit early-life gut genomes for residual clusters and useful strain references. The discovery cohort was children under three; this is sequence coverage, not a pediatric reference interval. [Study](https://doi.org/10.1038/s41467-022-32805-z), [representative sequences](https://zenodo.org/records/6969520). |
| R8 | **HumGut2: 31,225 genomes at 97.5% ANI dereplication** | Gut-associated reference diversity supplement after overlap removal. These are genomes, not 31,225 species. Its local integer IDs are not automatically NCBI taxids. [Author repository](https://github.com/larssnip/HumGut), [metadata](https://arken.nmbu.no/~larssn/humgut/HumGut2.tsv), [FASTA archive](https://arken.nmbu.no/~larssn/humgut/HumGut2.tar). |
| R9 | **HROM: 72,641 near-complete genomes / 3,426 oral species** | Targeted supplement for oral-origin organisms detectable in stool. Admit residuals and useful member genomes; oral provenance does not establish disease or exclude a valid stool detection. [Catalog](https://www.decodebiome.org/HROM/), [genome files](https://www.decodebiome.org/HROM/listdir.php?directory=data%2Fgenome_catalog), [project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1206836). |
| R10 | **NCBI RefSeq + selected GenBank bacterial/archaeal assemblies** | Obtain accession-versioned missing genomes, near relatives, named strains and linked BacDive genomes. Extend beyond the existing 486-target pathogen bundle where coverage warrants it. [Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/), [assembly metadata](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/). |

The mOTUs site lists 124,295 units while the current database 4.1 tutorial reports 124,300: record the exact downloaded artifact's count rather than forcing either figure. [Current metadata example](https://www.motus-tool.org/profiler/tutorials/motus_classify.html).

GlobDB already consolidates resources including mOTUs, SPIRE, GEM, gcMeta, MGnify and HRGM2. Do not sum their catalog counts. Older GTDB Kraken builds are not new sequence coverage beyond the installed R232 database. A Kraken2 adapter is valuable because it searches reads differently, not because an older GTDB label makes its catalog larger.

## 3. Prepare references once and measure actual added coverage

Extend existing reference preparation commands; make acquisition resumable. Lock tool/container versions, reference release, source URLs, retrieval timestamp, sizes, publisher checksums where supplied, local SHA-256, taxonomy release, build parameters and licenses. Sample runs never resolve moving `latest` aliases.

For each R5–R10 supplement:

1. Download metadata and species representatives first. Match exact accession versions and sequence checksums against GTDB232, GlobDB and other supplements.
2. Reconcile known source cluster memberships. Compare unmatched representatives using genome-wide ANI **and alignment fraction**; apply GTDB species-specific radii and existing representative assignments where applicable. For provisional non-GTDB clusters, start with 95% ANI and at least 65% alignment fraction of the shorter genome as a project clustering rule, preserve both directional fractions, and flag boundary/incomplete-genome cases instead of forcing equivalence. Do not use transitive single-link chaining to collapse distinct species.
3. For a known species, retain additional member genomes when they add previously unrepresented sequence diversity or improve discrimination. Dereplicate redundant members; do not use an arbitrary one-genome-per-species cap. Select diverse members from metadata/sketch distances and fetch more automatically for poorly matching candidates.
4. Record `already_represented`, `new_cluster`, `additional_strain_reference`, `better_reference`, `unresolved_mapping`, and `excluded_quality` counts. Retain quality and source provenance. Reference contamination or incompleteness must influence usable regions and confirmation, not silently create new organisms.
5. Build the **combined gut bacterial/archaeal rescue panel** from these representatives and admitted members. Include relevant GTDB near-neighbors, already-detected alternatives, host/PhiX and food/background competitors. These are one mapped reference universe; do not concatenate incompatible integer taxonomies.

GlobDB uses different operational cross-source clustering rules, including 96% ANI / 50% aligned fraction and GTDB priority. Preserve its source identities; its units are not automatically interchangeable with GTDB species. The project crosswalk records the relationship rather than rewriting a publisher's clustering result.

### Broad reference downloads

Reuse installed Sylph **1.0.0** and GTDB232. Add the published approximately **30-GB** GlobDB two-stage index:

`http://faust.compbio.cs.cmu.edu/sylph-stuff/globdb_r232_sylph_v2_c200.syl2db`

Obtain its matching taxonomy from:

`https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/taxonomic_profiling/globdb_r232_taxonomy_sylph.tsv.gz`

[Authoritative Sylph database list](https://sylph-docs.github.io/pre%E2%80%90built-databases/), [GlobDB fixed release root](https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/), [Sylph changelog](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md). Pin hashes and matching metadata; prefer an authenticated HTTPS publisher mirror where available or rebuild from verified source assets if integrity cannot be established. `.syl2db` requires Sylph 1.0+.

Use the full broad sketch without first mirroring hundreds of GB of GlobDB FASTAs. Retrieve candidate and near-neighbor genomes using GlobDB dictionaries and source manifests, including publisher-local identifiers: many additions have no NCBI/ENA accession. Use per-genome downloads where available; otherwise perform a one-time resumable, checksum-verified source/archive download and cache extracted subsets for confirmation. Existing GTDB and GlobDB overlap heavily; reuse sample sketches where supported, but preserve separate run provenance. Index size is not a peak-RAM guarantee.

## 4. Required detection and confirmation workflow

All added discovery methods below are part of the completed expansion. Data-dependent strain reconstruction and assembly run when their coverage conditions are met; those conditions are not feature flags for deferring implementation.

**A. Inputs and discovery.** Feed all lanes the same quality-controlled, host-filtered FASTQs, including supported singleton handling. Record retained fragment counts and input hashes. Run upgraded Jan26, existing GTDB232/Sylph, new GlobDB/Sylph, mOTUs4 and the gut rescue classifier. Do not restrict discovery to previously unclassified reads: incorrect earlier assignments could otherwise conceal new taxa. Reports/abundance JSON cannot replace reads for new sequence searches.

```bash
# Upstream invocation patterns; bind these variables in the existing runner.
metaphlan --install --index mpa_vJan26_CHOCOPhlAnSGB_202605 --db_dir "$OB_MPA_DB"
metaphlan "$OB_R1,$OB_R2" --input_type fastq --index mpa_vJan26_CHOCOPhlAnSGB_202605 --db_dir "$OB_MPA_DB" --mapout "$OB_MAPOUT" -s "$OB_SAM" --nproc "$OB_THREADS" -o "$OB_MPA_TSV"
sylph profile -d "$OB_GLOBDB_INDEX" -1 "$OB_R1" -2 "$OB_R2" -t "$OB_THREADS" > "$OB_GLOBDB_TSV"
motus downloadMGDB -db "$OB_MOTUS_DB"
motus profile -db "$OB_MOTUS_DB" -f "$OB_R1" -r "$OB_R2" -n "$OB_SAMPLE" -t "$OB_THREADS" -o "$OB_MOTUS_TSV"
```

Install/download steps run once through reference preparation. Check help against pinned binaries and smoke-test invocations. Sources: [MetaPhlAn](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4), [Sylph](https://sylph-docs.github.io/install%2Bquickstart/), [mOTUs](https://github.com/motu-tool/mOTUs).

**B. Whole-genome read rescue.** Implement **Kraken2** against the combined gut panel, using [Struo2](https://github.com/leylabmpi/Struo2) or equivalent existing build infrastructure. Use a single explicit custom taxonomy with reversible mappings to source accessions/GTDB/native IDs. Calibrate `--confidence` and `--minimum-hit-groups`; retain `--report-minimizer-data` and fragment assignments. Kraken confidence is not a probability. Build and run **Bracken** with matching database, k-mer and read-length settings where its native abundance estimates are reported; it is an abundance redistribution step, not independent confirmation. [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual), [Bracken](https://github.com/jenniferlu717/Bracken).

Do not build an enormous second all-GlobDB Kraken index merely to duplicate the broad sketch. The required Kraken lane searches the focused gut panel and informative member genomes; broad unexpected candidates are handled by competitive confirmation. Use uncompressed or appropriately decompressed inputs as required by each pinned tool.

**C. Competitive confirmation.** Validate every newly added species and every discordant/confusable candidate using the full eligible reads against candidate genomes, closest alternatives and relevant decoys together. Reuse the existing alignment infrastructure, extended beyond the pathogen target list. Measure independent fragments, uniquely discriminating regions, alignment identity, genome breadth/depth, locus distribution, coverage evenness, and alternative assignments. Conserved rRNA, mobile elements, repeats, or one suspect MAG contig cannot independently establish a species. Multiple strains of one species must compete as well.

Train operating points on separate tuning truth sets by depth, divergence and genome completeness; lock them before final testing. A single high-quality discovery lane plus discriminating sequence evidence can establish presence: requiring agreement of all profilers would discard their complementary gains. Conversely, three tools matching the same ambiguous sequence are not three independent biological observations. If no representative genome exists, validated species-specific markers may support a native cluster; otherwise retain the best-supported higher rank.

**D. Unrepresented-lineage rescue.** Add [SingleM](https://wwood.github.io/singlem/) with the [GlobDB metapackage](https://globdb.org/downloads) on retained reads. Use its unresolved marker lineages to trigger targeted assembly/binning when coverage supports reconstruction; reuse existing assemblers, or integrate a pinned metagenomic assembler/binning workflow if absent. Apply completeness, contamination and chimerism checks, competitively remap reads, and dereplicate recovered bins against all reference concepts. Report insufficient-depth lineages at their supported rank. A marker fragment or several contigs are not automatically several organisms. [Primary method](https://doi.org/10.1038/s41587-025-02738-1).

## 5. Repair taxonomy merging and complete strain analysis

### Release-aware taxonomy

Upstream bridge pairs currently include:

| Native MetaPhlAn database | Supplied GTDB mapping |
|---|---|
| Jun23_202403 | R207 |
| Jan25_202503 | R220 |
| Jan26_202605 | R226 |

There is no shipped R232 bridge in the checked upstream utility inventory. [Upstream utilities](https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils). The current Jan25/R220 file is not a valid general bridge for Jun23→R232.

Preserve native observations and create versioned, accession/genome-backed crosswalks. The proper Jan26/R226 mapping can be an intermediate step; resolve onward to R232 using assembly memberships, representative/member comparisons and release taxonomy. Classify or compare missing genomes against the pinned R232 reference. Never relabel R220/R226 names as R232 or join SGBs from different releases solely by numeric suffix.

Store `exact`, `synonym`, `parent`, `overlap`, `split`, `merge`, and `unresolved` relationships with supporting accessions and comparison evidence. Ambiguous splits remain complexes, not multiple certain species or forced one-to-one matches. Correct historical Jun23 reconciliation separately. Invalidate only affected names, merges, report counts and dependent caches; preserve raw detections and calibrated MP3 inputs. Bridge repair is a naming/reconciliation gain, not new sequence detection.

### Strains

`sample2markers` supplies consensus-marker inputs. Complete **StrainPhlAn**: obtain matching Jan26 SAMs, reconstruct consensus markers, extract clade markers, add suitable reference genomes, filter eligible samples/references, align and infer comparative phylogenies. Do not reuse Jun23 consensus objects as Jan26-compatible inputs. [Official workflow](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4).

Report comparative placement, closest-reference similarity, mixture evidence, and validated named subtype/pathotype separately. Insufficient marker coverage yields `unresolved`. A nearest reference is not proof of an exact resident strain; a consensus may conceal coexisting strains. Use competitive whole-genome evidence for candidate strain distinctions where feasible. Do not infer pathogenicity or resistance-carrier linkage from a phylogenetic neighbor alone. Strain detections never inflate species richness.

## 6. BacDive enrichment and merged report contract

**BacDive is downloadable and useful**, but it is primarily strain/phenotype metadata, not another FASTQ classifier. Its dashboard currently lists **102,187 strains / 22,480 species**. Use it for culture-collection identities, linked sequence accessions, oxygen requirements, substrates, products and cited reference-strain traits. [Dashboard](https://bacdive.dsmz.de/dashboard), [API](https://api.bacdive.dsmz.de/), [export documentation](https://hub.dsmz.de/wiki/bacdive/export/).

Use API v2, freely accessible without registration since February 2026; v1 is frozen at April 2025. Resolve `/v2/sequence_genome/{accession}`, `/v2/taxon/{genus}/{species}` or culture identifiers; follow pagination. Fetch `/v2/fetch/{id1;id2;...}` in batches of at most 100. Cache raw responses, citations and timestamps with retries/backoff. Record `no_matching_record` where appropriate. Download linked genome sequences through R10 only when they add coverage/discrimination. A trait measured in a reference strain is not automatically a trait of the detected patient strain.

Preserve data licenses and attribution. GlobDB propagates CC BY-SA 4.0; mOTUs data use CC BY 4.0. [BacDive terms](https://bacdive.dsmz.de/about) state CC BY 4.0 and request contact for commercial use; retain these terms alongside linked sequence-source licenses.

Extend existing equivalents of these records:

```text
ReferenceRelease:
  source_id, release_id, taxonomy_release, tool_version, sha256,
  urls[], license, count_unit, source_count, admitted_count, quality_policy
DetectionObservation:
  sample_id, input_sha256, lane_id, reference_release_id, native_id,
  canonical_id?, rank, status, support_metrics{}, alternatives[],
  abundance_value?, abundance_unit?, denominator?, raw_result_uri
CanonicalOrganism:
  stable_id, display_name, rank, domain, aliases[], identity_edges[],
  observation_ids[], strain_results[], traits[], confidence_basis,
  count_category, unresolved_reason?, incremental_gain_category
```

Statuses: `supported`, `provisional`, `ambiguous`, `not_detected`, `not_assessed`. A failed lane/missing reference is not a negative. Count named species, unnamed species-level clusters, unresolved complexes and higher-rank findings separately; display strain counts separately. Collapse aliases and overlapping evidence without losing provenance.

Add rows to the existing organism inventory and detail pages, with no hidden top-N truncation in PDF/HTML/JSON. Display all accepted organisms, readable names, stable native identifiers, detection support, reference release and available interpretation. Unknown biological significance is a valid label, not a reason to hide a detected organism. Show a compact coverage summary and link to methods.

Preserve the existing calibrated abundance/scoring inputs. Store each lane's native units and denominator; **do not sum, maximize or average independently normalized percentages**. An additional supported organism can have `abundance=null` rather than disappearing. New detections must be available through the existing downstream data interface without silently treating detection-only values as calibrated percentiles. DNA detection does not establish viability, persistent residence or disease.

## 7. Acceptance tests and completion

These are engineering targets, not claimed achieved performance. Reuse stronger existing gates. Separate tuning and held-out test communities/strains; publish uncertainty by abundance, depth, divergence and reference availability. Freeze truth-to-canonical mappings first. Score old and new native detector outputs through the same corrected canonical mapping; compare the as-deployed faulty bridge separately as a reporting repair. This separates added sequence detection from taxonomy cleanup. Parent/complex calls are not exact-species true positives.

| Gate | Required result |
|---|---|
| Installed baseline | Preserve a reproducible MP3 + Jun23 + GTDB232 baseline. New manifests show Jan26, GlobDB, mOTUs4, gut rescue panel and SingleM; no silent fallback to older references. |
| Reference contribution | Each supplement reports new clusters, extra strain references, better references, overlap and unresolved mappings. No summed catalog inflation or claim that all GlobDB additions are gut organisms. |
| Correct merging | Fixtures cover the supplied Jun23/Jan25/R220/R232 mismatch, aliases, taxonomic splits/merges and overlapping clusters. One organism renders once; an unresolved complex is not fabricated species multiplicity. |
| Precision | Held-out species precision **≥99% overall and ≥99% for newly added calls**, with confidence intervals and false positives/sample reported. Include close relatives, contaminated MAGs, shared genes, food/host backgrounds and strain mixtures. |
| Broad capability | At least 20 independent 600-species bacterial/archaeal communities, half even and half long-tailed; truth organisms each ≥0.02% microbial DNA mass. Test 8M and 30M retained microbial paired fragments, five read/error seeds. Target mean recall ≥95% and ≥98%, respectively, within each community class, while meeting precision and ≤5 false supported species per community. Report exact-reference and held-out-strain tiers separately. |
| Actual incremental value | Against the **full installed union**, demonstrate positive paired recall gain at matched precision, with a 95% confidence interval above zero on prespecified catalog-gap/strain-divergence fixtures. Publish ordinary gut-community performance and each lane's ablation separately. Do not tune final fixtures after inspecting results. |
| Rare/open-set behavior | Publish sensitivity below 0.02%; challenge with species excluded from every admitted catalog. Unsupported nearest-species matches must remain unresolved. Test simulated organism-free host/defined-food backgrounds and real blanks; zero false supported microbial calls on clean simulated negatives. Missing real blanks are disclosed. |
| Strain capability | Exercise full comparative outputs, mixtures, insufficient coverage and wrongly named nearest references. Consensus-only output cannot satisfy this gate. No strain counts enter species richness. |
| Data/report integrity | Accepted inventory counts reconcile with JSON, HTML and PDF; null abundance, unnamed species and long names render. Existing calibrated inputs remain unchanged. External HTML links use `target="_blank" rel="noopener noreferrer"`. |
| Reproducibility | Input/tool/reference/crosswalk changes invalidate the right caches. Interrupted preparation resumes; record actual peak RAM, disk and runtime. No unbounded repeated whole-reference downloads. |

For confidence intervals, resample independent communities/strain sets; seeds and depth tiers are technical replicates, not independent biological samples. Benchmark and physical-mock starting points: [CAMI II](https://doi.org/10.1038/s41592-022-01431-4), [MBARC-26](https://doi.org/10.1038/sdata.2016.81), [MBARC reads](https://www.ebi.ac.uk/ena/browser/view/SRR3656745). Resolve exact compositions and library metadata; do not duplicate reads to manufacture depth.

Regenerate every available project sample and deliver a per-sample comparison of **newly supported species**, renamed/deduplicated existing species, added strain information, rejected calls with reasons, and unresolved candidates. These unknown-truth samples demonstrate practical changes, not sensitivity or accuracy. A specimen need not gain species for an implementation to be correct, but the expansion must demonstrate genuine benchmark gains.

Complete all adapters, source acquisition, reference locks, crosswalks, competitive confirmation, strain continuation, report integration and tests. Failed gates require correction and rerun, not silently lowering thresholds or declaring unsupported hits detected. The completed report states measured capability and the reference releases used.
