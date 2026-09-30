# BUILD_SPEC_v05.0 — Broad stool pathogen and opportunist sequence screening

**Executable implementation specification · Research cut-off: 7 September 2026**

## 0. Implement this release

Add a dedicated **Pathogens, parasites and opportunistic organisms** system to the existing stool metagenomics application. Screen bacterial pathogens and pathotypes, protozoa, helminths, microsporidia, fungi, other relevant eukaryotes and viruses. Add separate virulence/toxin-gene and antimicrobial-resistance evidence. Search broadly; verify candidate assignments; show what was assessed, what was found, what was ambiguous and what this particular assay could not assess.

This document is the entire v05.0 handoff. It contains the target catalogs, data sources, acquisition rules, schemas, algorithms, report contract, build sequence and acceptance criteria. **No companion research report, source ledger or supplemental specification is required.** The coding agent should create normal application source files, tests and versioned reference manifests in the repository while implementing it. Large public sequence databases are acquired from the sources below; their actual bytes cannot sensibly be embedded in a Markdown specification.

The coding agent's supplied description of the current repository is a **reported audit, not independently verified repository access**. Begin by inspecting the real code. Preserve all existing disease profiles, organism/ecology panels, interventions, microbiome-age features and report sections that are actually implemented. This release adds a pathogen surveillance branch; it does not replace community profiling or convert disease-concordance scores into infection diagnoses. Implementation of v05.0 does not require reading the earlier research handoffs.

### 0.1 Required outcomes

1. Reprocess eligible FASTQs through a dedicated pathogen branch; an abundance table alone cannot support the new analysis.
2. Search broad genomic references, including harmless relatives and background organisms, and a curated human-relevance catalog. The displayed catalog is not the search engine's entire universe.
3. Explicitly close major worm, protozoan, fungal and human-virus reference gaps. A tool name or genus in a database is insufficient evidence of target coverage.
4. Retain supported low-abundance sequence findings without requiring a complete assembled genome. Preserve ambiguous and marker-only findings with their true limitations.
5. Separate organism evidence, pathotype/toxin evidence, AMR evidence, physical host linkage and clinical meaning.
6. Report DNA-virus capability now. Implement an optional validated RNA/cDNA input branch. A conventional DNA-only library must show RNA-only targets as **not assessed**, not negative.
7. Show a front-loaded pathogen summary, detailed evidence cards and the complete target-coverage ledger within the report. No green “all clear,” infection-free score, donor pass or arbitrary “badness percentage.”
8. Every installed target has a source, a reference-coverage status, an assay-eligibility rule, a technical detection rule and a validation status. A missing reference produces an explicit gap, never a fabricated negative.
9. Existing samples can be analyzed without waiting for new clinical trials. Their outputs are **research sequence findings**, with the applicable validation limitations. Promoting a finding to a validated clinical assay claim is a separate, evidence-based gate.

### 0.2 What “all” means operationally

Implement an extensible, versioned catalog covering established human enteric pathogens, rare documented enteric organisms, opportunists, toxin/pathotype determinants and unexpected human-pathogen sequence findings. Include broad discovery against additional reference taxa. Every release publishes the exact searchable reference inventory and explicit gaps.

No finite database contains every strain or unknown organism. More importantly, inclusion in a database does not establish that the organism sheds DNA into this specimen, survives extraction, appears in the sequenced aliquot or can be distinguished from relatives. “All known targets searched where the assay and references permit” is a defensible engineering objective. “All infections reliably excluded from any stool DNA sample” is not an available capability.

### 0.3 Corrections to the reported current-state assumptions

| Reported assumption | Verified correction and implementation consequence |
|---|---|
| Many thousands of MetaPhlAn bacterial species means comprehensive pathogen screening | Community marker coverage is not a validated clinical target list. Inspect actual package/database/flags and per-target markers. Add a separate pathogen workflow. |
| MetaPhlAn 4 has no viruses | Version-dependent: 4.1 introduced `--profile_vsc`; subsequent releases expanded viral outputs. These predominantly microbial-virus catalogs are not comprehensive human enteric-virus screens. Do not remove working phage outputs, but do not count them as human-pathogen coverage. [Official announcement](https://forum.biobakery.org/t/announcing-metaphlan-4-1-release-new-virome-database-and-sgb-database-update/6668), [changelog](https://github.com/biobakery/MetaPhlAn/blob/master/CHANGELOG.md). |
| EukDetect plus Kraken2 PlusPF automatically solves parasites and fungi | Neither guarantees the requested species. Current EukDetect2 lacks the major target worms and several important intestinal protozoa; PlusPF adds RefSeq fungi/protozoa, not a comprehensive helminth library. Use actual accession manifests and custom references. |
| All bad organisms should appear in DNA data | RNA-only genomes require an appropriate RNA/RT workflow. Tissue-restricted infection, intermittent shedding, low biomass, extraction bias, reference gaps and insufficient sequencing can also prevent detection. |
| A species is inherently bad wherever detected | Many requested bacteria and fungi are colonizers, food-associated organisms or conditional opportunists. Pathogenic strain/trait, body site and host context matter. Preserve findings without manufacturing infection. |

## 1. Scientific basis and decisions that follow from it

The following are the pivotal implementation sources, not borrowed performance guarantees.

| Evidence | Finding relevant to this build | Required implementation decision |
|---|---|---|
| Parks et al., February 2026, [176-target stool assay](https://www.frontiersin.org/journals/cellular-and-infection-microbiology/articles/10.3389/fcimb.2026.1759322/full) | Clinical comparisons covered 22 targets in 510 stools; contrived-sample experiments covered 19; all 176 had computational evaluation. The assay was DNA-based. Clinical performance varied substantially by target, with especially poor H. pylori sensitivity. | Store target-specific validation scope. Do not label all 176 clinically validated or copy this commercial assay's performance to our implementation. Public clinical reads: [PRJNA1200893](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1200893). |
| Cunningham-Oakes et al., 2025, [INTEGRATE](https://link.springer.com/article/10.1186/s13073-025-01478-w) | Large gastroenteritis study with paired metagenomic/metatranscriptomic data; 1,067 participants and 985 paired datasets. | Benchmark DNA and RNA separately; pair samples by participant/specimen. [PRJEB62473](https://www.ebi.ac.uk/ena/browser/view/PRJEB62473). |
| Haugum et al., 2025, [clinical/spiked stool study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0331288) | Twelve clinical stools, 36 spiked specimens and seven controls; assembly/MAG-based detection missed findings available from reads or the full assembly. | Assembly is additional context, not a mandatory presence gate. [PRJNA1218764](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1218764). |
| Peterson et al., 2022, [304-stool benchmark](https://www.mdpi.com/2076-2607/10/2/441) | Tool/database choice materially changed bacterial pathogen detection performance. | Include a published targeted-marker/k-mer comparator and orthogonal diagnostic truth. Raw accession not verified for this spec: resolve from the paper rather than inventing one. |
| [ParaRef, 2025](https://doi.org/10.1186/s13059-025-03818-w) | A decontaminated endoparasite genome resource; contaminating bacterial/host sequences in parasite references can produce false assignments. | Use screened parasite references, near-neighbor/background competition and reference-region masks. Do not treat a large count of matches to a contaminated reference as proof of a worm. |
| Avershina et al., 2025, [fungal-tool benchmark](https://link.springer.com/article/10.1186/s40168-025-02048-3) | Major variation in fungal reference coverage and detection; the tested historical PlusPF snapshot was limited, and even the strongest method missed many simulated fungi. | Supplement whole genomes with fungal markers and validate current snapshots. [Reproducible mocks](https://github.com/Rounge-lab/mock_mycobiome). |
| Bradford et al., 2024, [false-positive analysis](https://link.springer.com/article/10.1186/s12859-024-05952-x) | Database composition, close relatives and classifier thresholds affect pathogen false positives. | K-mer classification is candidate generation. Require competitive alignment and organism-specific resolution. |
| Goulet et al., May 2026, [CroCoDeEL](https://www.nature.com/articles/s41467-026-72637-9) | Batch cross-contamination can create low-abundance findings; detection depends on depth and profiling method. | Add batch context and optional contamination screening, but do not substitute it for blanks or automatically erase shared organisms. |

**Evidence types must remain distinct:** a reference genome establishes searchability; an in-silico mixture tests software; a DNA spike tests downstream library steps; a whole-organism stool spike tests more of extraction; a clinical comparison evaluates the matched assay/population. None automatically supplies the others.

## 2. Scope, clinical interpretation classes and assay eligibility

### 2.1 Target classes

Use these exact machine values:

| `interpretation_class` | Meaning | Default participant-facing placement |
|---|---|---|
| `established_enteric` | Established human GI pathogen or a defined enteric pathogen group | Supported findings in main Pathogens section; confirmation context |
| `toxin_or_pathotype_dependent` | Species-level identification alone cannot establish pathogenic trait | Organism and determinant evidence shown separately |
| `conditional_opportunist` | Can cause disease in particular hosts/sites; carriage alone is common or possible | Opportunistic organisms subsection, neutral detection wording |
| `rare_enteric` | Documented but uncommon human intestinal/hepatobiliary infection | Main section when supported, with specimen/clinical limitations |
| `uncertain_enteric_role` | Human stool association or disputed disease attribution | “Detected; clinical significance uncertain” |
| `extraintestinal_watch` | Important organism for which stool is not a general diagnostic specimen | Unexpected findings requiring appropriate confirmation; never negative-screen reassurance |
| `background_or_decoy` | Commensal, dietary, environmental, harmless relative or nonhuman host sequence | Coverage/other-organisms view; no pathogen alarm |

These classes describe clinical relevance, not technical confidence. A very strong Candida sequence finding can still be a colonization finding; a weak Giardia signal remains technically uncertain despite established pathogenic relevance.

### 2.2 Required assay manifest

```json
{
  "schema_version": "pathogens.assay.v1",
  "sample_id": "existing-internal-sample-id",
  "specimen_type": "stool",
  "collected_at": null,
  "participant_age_years": null,
  "nucleic_acid_protocol": "DNA",
  "reverse_transcription": false,
  "library_selection": "shotgun",
  "extraction_protocol_id": null,
  "input_stool_mass_mg": null,
  "preservative": null,
  "mechanical_lysis": null,
  "wetlab_batch_id": null,
  "library_batch_id": null,
  "sequencing_run_id": null,
  "read_layout": "paired",
  "sequencer": null,
  "read_length_summary": null,
  "control_sample_ids": [],
  "spike_in_id": null,
  "recent_antimicrobials": null,
  "recent_probiotics": null,
  "recent_live_oral_vaccine": null,
  "symptoms": [],
  "immunocompromise": null,
  "raw_fastq_sha256": [],
  "host_removal_bundle_id": null
}
```

Allowed protocol values: `DNA`, `RNA_with_RT`, `total_nucleic_acid_with_RT`, `unknown`. `null` means unknown, never false. Verify the upstream laboratory protocol; the FASTQ alphabet cannot distinguish DNA from reverse-transcribed RNA. `library_selection` must distinguish unbiased shotgun, poly-A selected RNA, rRNA-depleted RNA, amplicon and targeted capture. An amplicon result cannot inherit shotgun target coverage.

| Target nucleic acid | DNA-only | RNA with RT | Total nucleic acid with RT | Unknown protocol |
|---|---|---|---|---|
| Bacterial/fungal/parasite genomic DNA and DNA viruses | Eligible for DNA sequence search, subject to extraction/reference limits | Transcript findings possible; not equivalent to validated genomic-DNA screening | Eligible only for components the actual preparation preserves and validates | Findings may be explored; no assay-wide negative claim |
| RNA-only enteric viruses | **Not assessed** | Eligible if library chemistry captures their RNA | Eligible if RT and RNA preservation are documented | Not assessed until protocol is resolved |
| Retroviruses | Proviral DNA may be searchable, but stool is not an HIV/HTLV diagnostic assay | RNA may be searchable with appropriate protocol | Component-dependent | No diagnostic inference |
| RNA-derived expression/viability | DNA cannot measure it | Transcript detection does not independently establish infectious viability | Same limitation | Not assessed |

For a DNA-only sample, a read that resembles an RNA virus is retained as an **unexpected assay-inconsistent signal**, with possible contamination, synthetic material, homology or metadata error considered. It must not unlock RNA-virus negatives or silently become a supported RNA infection.

**Age:** presence screening for pathogen sequence is not an adult microbiome-reference percentile. Do not suppress all pathogen searches in minors under the existing adult disease-profile abstention rule. Keep age-specific clinical interpretation, assay validation and confirmations explicit. For example, C. difficile carriage in very young children requires age-aware interpretation. No antimicrobial or antiparasitic dosing is generated.

## 3. Reference architecture and exact acquisition sources

### 3.1 Required bundles

Implement five independently versioned reference bundles under one release lock:

1. `broad_taxonomy`: broad bacteria/archaea, fungi, protists, DNA/RNA viral sequences and background/decoy genomes. Clinical labels do not determine inclusion in competitive search.
2. `human_relevance_targets`: the seed catalog in sections 8–11, curated expansion rules and clinical interpretation metadata.
3. `parasite_fungal_supplement`: screened helminth/protozoan/fungal nuclear genomes, organellar genomes and targeted loci missing or inadequately represented in the broad bundle.
4. `determinants`: curated AMR, virulence, toxin and specialist locus references, with exact alleles and source-specific rules.
5. `verification`: target-specific informative regions, near-neighbor sets, repeat/contamination masks, assay-specific calling profiles and validation records.

A target is not “covered” simply because its name exists in taxonomy. Build coverage from actual sequence accessions, usable unmasked sequence and validated taxonomic resolution. Store `genome`, `organelle_only`, `marker_only`, `protein_only`, `no_usable_reference` separately. A genus-level hit cannot satisfy species-level coverage.

### 3.2 Primary sequence and taxonomy resources

| Resource | Acquisition entry point | Use and constraints |
|---|---|---|
| NCBI RefSeq / GenBank | [Datasets](https://www.ncbi.nlm.nih.gov/datasets/), [genome download guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/), [RefSeq FTP](https://ftp.ncbi.nlm.nih.gov/refseq/), [GenBank assemblies](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/) | Primary whole-genome and accession metadata. RefSeq is a preferred baseline, not a rule that excludes valid rare GenBank-only targets. Retain assembly quality, source, status and contamination information. |
| NCBI Taxonomy | [Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz), [accession-to-taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/) | Lock names, nodes, merged and deleted IDs. Preserve scientific names and historical aliases. Do not mix NCBI IDs with GTDB/SGB identifiers without explicit mappings. |
| NCBI Virus | [Portal](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/), [viral downloads](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/virus-download/) | Nucleotide/protein records, host, molecule and segments; add representative diversity beyond one exemplar. Isolation from a human sample is not proof of human viral tropism. |
| ICTV VMR / MSL | [VMR](https://ictv.global/vmr), [July 2026 version archive](https://zenodo.org/records/21694279) | Current verified workbook: `VMR_MSL41.v1.20260729.xlsx`. Use exemplar accessions, molecule/host and taxonomy mappings. The earlier July 21 workbook was superseded. VMR is not a pathogenicity list. |
| NCBI targeted loci | [Targeted loci](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/), [fungal ITS](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/ITS_process) | Curated marker references, including partial-coverage fallback. Distinguish marker evidence from whole-genome evidence. |
| Broad Kraken2 indexes | [Official index table](https://benlangmead.github.io/aws-indexes/k2) | Current table has June 2026 PlusPF/PlusPFP snapshots. Preserve `library_report.tsv`, `names.dmp`, `nodes.dmp`, index checksums and source date. PlusPF adds fungi/protozoa; PlusPFP also adds plant sequences. Neither guarantees all worms. |
| NCBI Pathogen Detection | [Isolate browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates), [documentation](https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/) | Diverse pathogen assemblies and linked isolate/AST metadata. Isolate performance is not stool-mixture performance; deduplicate outbreak clones for validation. |
| BV-BRC | [Data documentation](https://www.bv-brc.org/docs/system_documentation/data.html), [downloads](https://www.bv-brc.org/docs/quick_references/ftp.html), [ViPR/IRD migration](https://www.bv-brc.org/docs/quick_start/ird-vipr_bv-brc_mapping.html) | Successor to PATRIC/ViPR/IRD. Respect source annotations and rights; current transfer documentation uses FTPS. Do not depend on retired legacy sites. |

For broad searches, an uncapped current PlusPFP snapshot is a practical starting point when its resource footprint is affordable. Alternatively build an equivalent reproducible NCBI-taxonomy bundle. Add the missing eukaryotes separately and perform final competition across all relevant bundles. Do not report a reduced 8/16-GB index as equivalent sensitivity to the uncapped reference. Never allow a memory failure to silently substitute a smaller database or a subset of reads.

### 3.3 Parasite and fungal resources

| Resource | Acquisition | Required behavior |
|---|---|---|
| ParaRef | [Paper](https://doi.org/10.1186/s13059-025-03818-w), [code](https://github.com/Schroeder-Group/ParaRef), [versioned initial archive](https://doi.org/10.5281/zenodo.13744644) | Strong base for screened endoparasite genomes. Read actual FASTA/taxonomy/mask manifest and species coverage; supplement human targets absent from it. Its research thresholds are not universal clinical thresholds. |
| WormBase ParaSite | [Portal](https://parasite.wormbase.org/), [WBPS19 downloads](https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/) | Worm nuclear genomes, annotations and proteins. Verified release remains WBPS19, March 2024; development funding pause means use NCBI and newer primary genomes to fill gaps. Preserve BioProject and assembly version. |
| VEuPathDB family | [Portal](https://veupathdb.org/), [CryptoDB](https://cryptodb.org/), [GiardiaDB](https://giardiadb.org/), [AmoebaDB](https://amoebadb.org/), [ToxoDB](https://toxodb.org/) | Curated protist genomes and annotations. Resolve current download links and upstream projects. Do not infer human intestinal coverage from the resource's genus or disease name. |
| EukDetect2 | [Official repository](https://github.com/allind/EukDetect), [March 2026 data](https://doi.org/10.5281/zenodo.19056625) | Additional fungi/protist marker evidence. Pin a tested release/commit and matching v2 database; v1/v2 assets are not interchangeable. Actual manifest lacks major target worms and several intestinal protozoa. |
| CORRAL | [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10084625/) | Additional eukaryotic marker-analysis approach and method comparator. Not a universal replacement for missing reference sequences. |
| UNITE | [Versioned repository](https://unite.ut.ee/repository.php) | ITS references and species-hypothesis identifiers. Pin fungal-only versus all-eukaryote release. ITS copy number and variable resolution prohibit cell-load inference. |
| FungiDB | [Downloads](https://fungidb.org/fungidb/app/downloads), [resource paper](https://academic.oup.com/genetics/article/227/1/iyae035/7634838) | Fungal/oomycete genomes and annotations; preserve per-dataset source and rights. |
| Candida Genome Database | [CGD](https://www.candidagenome.org/) | Curated Candida genes and nomenclature; supplement whole-genome data. |
| MycoBank | [MycoBank](https://www.mycobank.org/) | Synonyms/type information. Names and renamed genera must map to stable targets without double counting. |
| FunOMIC | [Resource](https://manichanh.vhir.org/funomic/), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9293737/) | Optional fungal marker/protein comparator; benchmark current accessible assets before enabling. |

EukDetect2 source snapshot inspected during research: `8d69014727b2c5956de30b0811644b6c3a7bd4f3` (July 16, 2026). The March 2026 database manifest contained 6,948 records; this is a database-record count, not a number of validated human pathogens. It lacked **Cystoisospora belli, Balantioides/Balantidium coli, Sarcocystis hominis/suihominis**, and major helminths. Related animal species must not silently substitute for these targets. Record this as an observed snapshot, then regenerate coverage for the installed snapshot.

ParaRef repository snapshot `6f20a948d70afd55dbef7bce67d124652563135f` pointed to [taxonomic mapping files](https://doi.org/10.5281/zenodo.21413835), distinct from the genome archive. Those mapping-file contents were not downloaded/verified during this research. The builder may adopt them only after resolving metadata, checksums, license, contents and compatibility with the selected genome archive; do not present their compatibility as already tested here.

### 3.4 Viral supplemental resources

| Resource | Acquisition | Use |
|---|---|---|
| RVDB nucleotide | [Archive](https://rvdb.dbi.udel.edu/previous-release), [U-RVDBv32.0](https://rvdb.dbi.udel.edu/download/U-RVDBv32.0.fasta.gz), [C-RVDBv32.1](https://rvdb.dbi.udel.edu/download/C-RVDBv32.1.fasta.gz) | Broad eukaryotic-virus discovery. May 2026 release banner and clustered filename differ; preserve the actual artifact identifier and sidecars. Supplement NCBI/ICTV; do not call all matches human pathogens. |
| RVDB-prot | [Pasteur resource](https://rvdb-prot.pasteur.fr/), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7492780/) | Optional protein/HMM rescue for divergent viral candidates; homology-only evidence does not support exact species or infection. |
| geNomad | [Code](https://github.com/apcamargo/genomad) | Optional assembled viral/plasmid discovery and classification. Not a human-pathogen classifier. |
| CheckV | [Primary method](https://www.nature.com/articles/s41587-020-00774-7) | Viral-contig completeness/host-contamination assessment; not a clinical presence validator. |

### 3.5 AMR, virulence and toxin resources

| Resource | Acquisition | Default decision |
|---|---|---|
| NCBI NDARO / AMRFinderPlus | [Database root](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/), [formats](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database), [program](https://github.com/ncbi/amr) | Mandatory default. Lock dated database plus compatible software. Retain `ReferenceGeneCatalog.txt`, hierarchy, `AMR_CDS.fa`, `AMRProt.fa`, HMMs and organism-specific mutation assets. Selected virulence coverage is not all virulence. |
| ResFinder / PointFinder | [Program](https://github.com/genomicepidemiology/resfinder), [ResFinder DB](https://bitbucket.org/genomicepidemiology/resfinder_db), [PointFinder DB](https://bitbucket.org/genomicepidemiology/pointfinder_db) | Optional corroboration and gap fill when database terms permit. Organism-specific mutation models require appropriate organism assignment. |
| CARD / RGI | [Downloads](https://card.mcmaster.ca/download), [analysis modes](https://card.mcmaster.ca/analyze), [terms](https://card.mcmaster.ca/about) | Optional licensed adapter. Current commercial-use restrictions require authorization under provider terms. Do not silently include restricted data or disable the open-resource core because this adapter is unavailable. |
| VFDB | [Database](https://www.mgc.ac.cn/VFs/main.htm), [downloads](https://www.mgc.ac.cn/VFs/download.htm) | Optional licensed adapter; current CC BY-NC 4.0. Preserve experimentally supported core A versus expanded B evidence. |
| VirulenceFinder | [Official code](https://bitbucket.org/genomicepidemiology/virulencefinder), [primary study](https://journals.asm.org/doi/10.1128/jcm.03617-13) | Selected bacterial genera; inspect database rights independently of program license. |
| Victors | [Resource](https://phidias.us/victors/), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6324020/) | Curated virulence experiments across hosts. Filter human relevance and experimental context; confirm reuse terms. |
| PHI-base | [PHI5](https://phi5.phi-base.org/), [2025 paper](https://academic.oup.com/nar/article/53/D1/D826/7908791), [versioned data](https://zenodo.org/records/16738930) | Gene/host/experimental-phenotype evidence; much concerns nonhuman hosts. A mutant phenotype is not a diagnostic gene rule. |
| StxTyper | [NCBI code](https://github.com/ncbi/stxtyper), [2026 paper](https://www.mdpi.com/2076-2607/14/8/1607) | Specialist assembled Shiga-toxin operon typing; not a raw-read taxonomic classifier. |
| Kleborate | [Code](https://github.com/klebgenomics/Kleborate), [documentation](https://kleborate.readthedocs.io/) | Specialist Klebsiella assembly interpretation, with host/assembly quality checks. |
| MIBiG | [Catalog](https://mibig.secondarymetabolites.org/) | Curated biosynthetic loci for targeted toxin-capacity annotations. Preserve quality/completeness annotations. |
| FungAMR | [Code/data](https://github.com/Landrylab/FungAMR), [interface](https://card.mcmaster.ca/fungamrhome), [2025 paper](https://doi.org/10.1038/s41564-025-02084-7) | Optional fungal resistance-associated variant evidence. Includes susceptibility observations, experimental/agricultural records and non-SNP mechanisms; not a ready-made all-fungi clinical resistance caller. |

Open-resource fallback for mandatory bacterial determinants: use the NDARO catalog and experimentally characterized nucleotide/protein records linked from the primary studies/FDA methods in section 8. If a determinant is absent from NDARO, acquire its explicitly characterized GenBank/RefSeq accession and record its identity/function evidence. **Do not substitute a gene-name search or an unlicensed copied VFDB/CARD record.** Unresolved exact allele provenance remains an explicit determinant reference gap; the rest of the system runs.

### 3.6 Rights and reproducibility

Record software license and data license independently. NCBI's [reuse policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) does not impose NCBI restrictions on molecular data but notes possible original-source rights. UNITE, the referenced ICTV archive and EukDetect2 archive provide attribution-based terms; retain the exact downloaded release license. RVDB has a [specific agreement](https://rvdb.dbi.udel.edu/RVDBonlinelicenseterms040921.pdf) with redistribution, disclaimer and nonendorsement conditions. Aggregating annotations through another site does not erase upstream rights.

Do not expose an approval prompt for routine open-data downloads. If a required optional licensed adapter is unavailable, mark it unavailable and complete the executable default path. Database rights inspection is a release-build step, not an excuse to leave all pathogen screening unfinished.

### 3.7 Deterministic acquisition procedure

The coding agent must implement `pathogens refs resolve`, `pathogens refs build`, `pathogens refs audit` and `pathogens refs lock` as application commands or equivalent repository entry points. These are **new interfaces to implement**, not commands claimed to exist today.

1. Resolve catalog names/aliases against a frozen taxonomy. Require explicit human review of ambiguous name-to-taxon mappings; continue other resolved targets. Never choose the first substring match.
2. Fetch assembly/sequence metadata first. Select reference/type representatives plus diverse clinically relevant strains/subspecies, without excluding rare targets merely because only scaffold/contig assemblies exist.
3. Download sequence records by accession.version and archive their metadata. A fallback species query is allowed during resolution, but the released build must contain explicit accession.version lists.
4. For records absent from genome packages, retrieve verified nuclear markers, ITS/18S/28S, mitochondrial or plastid references from GenBank/RefSeq/primary genome studies. Organellar evidence has its own resolution rules; an apicoplast is not a helminth mitochondrion.
5. Add host, common food, environmental eukaryote and near-neighbor competitors. Include pig/cow/sheep/chicken/fish and common dietary plants/fungi where relevant, plus laboratory vectors and common bacterial contaminants.
6. Deduplicate identical sequences while preserving all origin/alias records. Screen assemblies for foreign sequence contamination, low complexity, repeats and ambiguous taxonomy. Retain masks and rejection reasons. Shared regions remain eligible for higher-rank assignment, never species-specific evidence.
7. Generate informative regions against both within-genus relatives and broader cross-kingdom/background references. Within-genus uniqueness alone does not exclude bacterial sequence contaminating a worm assembly.
8. Build indexes with frozen parameters. Preserve complete source FASTAs and taxonomy for rebuilding. Never rebuild from mutable `latest` at sample-analysis time.
9. Audit every seed and dynamically added target: intended name/taxid, actual accessions, usable base count, route coverage, near neighbors, specificity limitations and reference gap reason.
10. Lock source URLs, retrieval times, accession versions, SHA-256, software/container digests, taxonomy release, masking algorithms and the complete target inventory. Old reports retain their old lock IDs.

Genome acquisition uses the documented NCBI Datasets CLI, for example `datasets download genome accession ACCESSION --include genome --filename PACKAGE.zip` after resolving an actual accession. For high-volume download use documented dehydrated packages and rehydration, with resumable transfers. Do not embed fabricated assembly IDs or hashes in this specification.

WormBase downloads are discoverable under WBPS19 by species and BioProject; many follow `species_lower_underscores/BIOPROJECT/species_lower_underscores.BIOPROJECT.WBPS19.genomic.fa.gz`. Parse the release listing and manifest rather than assuming every file follows this template. Keep the source's assembly identity when merging with NCBI.

## 4. Data contracts

### 4.1 Target record

```yaml
schema_version: pathogens.target.v1
target_id: parasite.giardia_duodenalis
display_name: Giardia duodenalis
aliases: [Giardia intestinalis, Giardia lamblia]
group: protozoa
interpretation_class: established_enteric
ncbi_taxids: []                 # resolved by refs resolve; empty is not search-ready
taxonomic_resolution: species_complex
allowed_nucleic_acids: [DNA]
stool_role: intestinal_shedding
reference_status: unresolved
reference_accessions: []
marker_accessions: []
near_neighbor_target_ids: []
sequence_call_profile_id: research_dna_v1
validation_id: null
clinical_validation_status: not_established_for_this_pipeline
clinical_source_urls:
  - https://www.cdc.gov/dpdx/giardiasis/index.html
report_template_id: intestinal_parasite
confirmation_options: [clinical_stool_NAAT, clinical_stool_antigen]
negative_limitation_codes: [sampling, shedding, extraction, depth, reference_diversity]
```

This example illustrates a **pre-resolution seed**, not a fully installed target. Compiler rule: a target with unresolved identifiers/accessions cannot receive `not_detected` or species-supported status. Resolved catalog records must contain valid accessions or an explicit `no_usable_reference` reason.

`group` enum: `bacteria`, `protozoa`, `helminths`, `microsporidia`, `fungi`, `other_eukaryotes`, `dna_viruses`, `rna_viruses`, `other_viruses`, `virulence`, `amr`. Microsporidia appear once under Parasites, with fungal relationship metadata; never counted twice.

For viral target and result records additionally require `host_category` (`human`, `bacterial`, `archaeal`, `fungal`, `protist`, `plant`, `other_animal`, `unresolved`), `host_evidence_level`, `host_evidence_source_ids`, `sample_isolation_host`, `molecule_type` and `segment_manifest_id`. Isolation from human stool does not set the biological host to human. Map uncertain viral attribution to `uncertain_enteric_role` or `background_or_decoy` as appropriate; do not create an undeclared interpretation-class enum.

### 4.2 Reference lock

Required fields: `bundle_id`, `created_at`, `catalog_version`, `taxonomy_snapshot`, `taxonomy_sha256`, `software[{name,version,source_commit,container_digest}]`, `assets[{source_url,retrieved_at,accession_version,sha256,license_id,sequence_type,taxid,mask_sha256}]`, `index_parameters`, `target_coverage`, `validation_profile_ids`, `build_command_log_sha256`.

`target_coverage` must distinguish `genome_supported`, `marker_only`, `organelle_only`, `protein_only`, `unresolved_taxonomy`, `no_usable_reference`, `license_unavailable`, `not_in_this_bundle`. A later bundle can fill another bundle's gap. Report combined coverage from eligible routes, not the worst individual-tool gap or a falsely optimistic union of incompatible assays.

### 4.3 Result object

```json
{
  "schema_version": "pathogens.result.v1",
  "sample_id": "sample-001",
  "target_id": "bacteria.escherichia_shigella_complex",
  "group": "bacteria",
  "display_name": "Escherichia coli / Shigella complex",
  "interpretation_class": "toxin_or_pathotype_dependent",
  "assay_eligibility": "eligible",
  "analysis_status": "completed",
  "reference_status": "genome_supported",
  "sequence_status": "supported_sequence",
  "contamination_status": "controls_unavailable",
  "resolution": "species_complex",
  "clinical_interpretation": "organism_sequence_not_infection_diagnosis",
  "validation_scope": "computational_research_rule",
  "calling_profile_id": "research_dna_v1",
  "reference_bundle_id": "resolved-release-lock-id",
  "unique_supporting_fragments": 12,
  "ambiguous_fragments": 31,
  "informative_regions_supported": 4,
  "informative_bases_covered": 810,
  "reference_breadth_fraction": 0.0002,
  "informative_region_breadth_fraction": 0.013,
  "median_alignment_identity": 0.997,
  "normalized_fragments_per_million": 1.2,
  "normalization_denominator": "10000000_qc_pass_nonhost_fragments",
  "absolute_load": null,
  "clinical_lod": null,
  "infection_probability": null,
  "linked_determinants": [],
  "unlinked_determinants": ["example-determinant-result-id"],
  "reason_codes": ["species_resolution_limited", "controls_unavailable"],
  "evidence_artifact_ids": [],
  "confirmation_options": ["clinician_selected_stool_NAAT_or_culture"]
}
```

Numbers above are **synthetic schema examples**, not measurements from SAMPLE2 or any other sample. Production output must derive all values from real evidence. `clinical_lod`, absolute load and infection probability remain null unless independently supported; do not fill them from literature values or read fractions.

### 4.4 Sequence status, precedence and missingness

Store orthogonal fields rather than a single overloaded red/green status:

- `assay_eligibility`: `eligible`, `ineligible`, `unknown`.
- `analysis_status`: `completed`, `not_run`, `failed`, `qc_failed`.
- `sequence_status`: `supported_sequence`, `marker_signal`, `candidate_signal`, `ambiguous_signal`, `not_detected`, `not_assessed`.
- `contamination_status`: `no_flag_under_applied_checks`, `suspected`, `confirmed_technical_artifact`, `controls_unavailable`, `not_evaluated`.
- `validation_scope`: `computational_research_rule`, `analytical_assay_validated`, `clinical_assay_validated` plus the exact validation record; no inherited whole-panel validation.

Precedence for participant display: failed/not-run/assay-ineligible/reference-unusable becomes **Not assessed** with reason; sequence evidence remains available in an assay-inconsistent or technical-review record. Suspected contamination is displayed beside the finding and downgrades its clinical prominence; it does not silently turn it into a negative. Ambiguous related-species evidence is displayed at the supported group/rank. `not_detected` requires an eligible assay, usable reference, completed searches and passed applicable QC; it means no signal meeting this documented analytical rule, not organism absence.

A child species can be unresolved while its parent group has supported evidence. Count the parent once. Do not count every candidate reference, synonymous name, viral segment or subtype as a separate infection.

## 5. Pipeline implementation

### 5.1 Repository audit and integration seam

Inspect the FASTQ ingestion, QC, human-read removal, MetaPhlAn invocation/database, DIAMOND panels, existing opportunist registry, result schema and report renderer. Record actual findings in the PR and machine-readable build manifest. Locate the earliest reusable QC-passed reads. Do not use only reads classified as microbes by MetaPhlAn: those would discard the very gaps being fixed.

Add a resumable task branch from eligible reads, caching by raw-input hash + assay manifest + software + reference lock + parameters. Add feature flag `pathogens_v05`; run migrations without deleting old results. If FASTQs are unavailable, emit “Pathogen screening not run: raw reads required,” retaining existing reports.

A reference orchestrator is [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler); [GMS metaval](https://github.com/genomic-medicine-sweden/metaval) provides useful taxid-driven read verification and coverage reporting. Reuse compatible components or implement the same contracts in the existing framework. Do not force a wholesale Nextflow migration. Metaval's custom-taxid route extends beyond viruses, but its default workflows and thresholds are not a universal clinical pathogen assay.

### 5.2 Processing order

1. **Input validation:** sample identity, paired-read synchronization, compressed-file integrity, assay metadata and checksums. Keep DNA/RNA libraries separate even if from one stool.
2. **QC:** adapter/quality processing compatible with the actual platform. Record input fragments, retained fragments, length distribution, ambiguous bases and low-complexity fractions. Do not abundance-filter rare reads or downsample the production pathogen branch.
3. **Human-read handling:** use the approved, pinned host-removal workflow. Record losses and assess its effect on known pathogen positives. Competitively evaluate host-homologous candidates locally; do not upload raw reads to public BLAST or external classifiers. Host-homologous evidence alone cannot support a pathogen. Never simply disable host filtering to improve sensitivity.
4. **Broad discovery:** run a frozen Kraken2/custom equivalent against broad references, retaining per-read assignments including ambiguous/unclassified reads. No `--quick` or unrecorded memory-capped fallback. Candidate generation may be permissive; it is not a positive call.
5. **Independent sensitive target search:** screen all eligible reads against curated informative-target/marker/locus references, including target families the broad classifier misses. This route must not depend on a prior Kraken/MetaPhlAn hit.
6. **Eukaryote route:** EukDetect2 or a tested equivalent for its represented fungi/protists; additional ParaRef/WormBase/NCBI/marker search for uncovered targets. A negative or absent EukDetect target is not a veto on other valid evidence.
7. **Competitive confirmation:** align reads against target references plus close relatives, host/food/vector/background decoys. Retrieve mates and relevant ambiguous reads, not only the first classifier's species-assigned reads. Confirm rank, unique informative support, breadth, depth, contamination and alternative explanations.
8. **Determinant search:** screen curated nucleotide genes in parallel; optionally use protein/HMM search for divergent candidates. Preserve partial gene signals. Assemble where feasible, then run AMRFinderPlus and specialist typing on appropriate assemblies.
9. **Optional discovery/assembly:** assemble candidate and unclassified microbial reads where resources permit. geNomad/DIAMOND/BLAST may reveal novel/divergent candidates. Homology-only or novel sequences enter a review queue, never acquire invented species/pathogenicity labels.
10. **Evidence adjudication:** apply the frozen target-specific or provisional research rules below. Store all underlying evidence and reasons. Generate the new report from canonical results, not raw classifier tables.

### 5.3 Competitive alignment and taxonomic resolution

Use an appropriate nucleotide aligner such as Bowtie2 for Illumina; minimap2 presets must match long-read or other supported data. Freeze software and scoring parameters. Candidate verification must evaluate competing alignments, not only the best hit in a pathogen-only reference.

Required per-target metrics:

- Number of **distinct supporting fragments**; overlapping mates are one fragment. For UMI libraries collapse by validated UMI rules. For non-UMI libraries count conservative duplicate families once for the provisional support gate: equivalent canonical template start/end/orientation after conspecific-reference normalization. Different read names alone do not establish independent support. Retain raw and duplicate-like counts; coordinate deduplication is not exact molecule counting and may conservatively merge genuine molecules.
- Aligned query fraction, alignment identity, base quality, mapping quality, best competing taxon/score and ambiguity.
- Distinct nonoverlapping informative regions and covered informative bases, separately for nuclear genome, organelle, marker, gene or viral segment.
- Whole-reference breadth and depth with explicit denominators; informative-region breadth is not whole-genome breadth.
- Distribution of support: repeat-only, single conserved locus, contaminant contig, organelle-only, strand/end anomalies and coverage hotspots.
- Fraction of support shared with close relatives and alternative taxonomic rank justified by the data.

Build all informative regions from actual sequences and frozen background comparisons. Do not enforce “one species per genus”: true coinfections and multiple colonizing relatives must survive if each has distinct evidence. Avoid using abundance to break unresolved ties. If all reads distinguish a species complex but not its members, the complex is the result.

### 5.4 Executable research detection policy

Clinical thresholds must eventually be target/assay validated. To make the first build executable now, implement the following **provisional research policy**, visibly labeled `computational_research_rule`. These values are engineering starting points to test, not thresholds borrowed from a clinical guideline and not a claim of diagnostic sensitivity.

```yaml
research_dna_v1:
  candidate_generation:
    kraken_confidence: 0.0
    minimum_hit_groups: 2
    quick_mode: false
    preserve_unclassified: true
  informative_alignment:
    min_aligned_query_fraction: 0.90
    min_identity: 0.95
    min_aligned_bases: 75
    min_target_level_mapq: 30
    conspecific_multimapping_alone_is_not_a_veto: true
    min_best_minus_other_taxon_alignment_score: 10
    exclude_low_complexity_or_masked_only_support: true
  genome_support:
    min_distinct_fragments: 5
    min_nonoverlapping_informative_regions: 3
    min_informative_bases_covered: 300
    require_resolution_specific_regions: true
  marker_or_organelle_support:
    min_distinct_fragments: 5
    min_nonoverlapping_supported_windows: 3
    min_informative_bases_covered: 300
    report_as: marker_signal
  absence_claim: no_supported_sequence_under_this_rule
  clinical_lod: null
  clinical_sensitivity: null
  clinical_specificity: null
```

Additional mandatory semantics:

1. MAPQ/score requirements apply to the frozen aligner and competitive reference. They are not comparable across arbitrary programs. Use **target-level** uniqueness against other taxa: equivalent placements in multiple strains of the same target must not be rejected merely because raw mapper MAPQ is low. Collapse equivalent references/placements or implement and benchmark a target-aware scoring adapter. Preserve unresolved coordinates/alleles separately. A reference lacking meaningful competitors is not species-reportable regardless of MAPQ; an ad hoc score must not be mislabeled a calibrated probability.
2. Genomic support requires three frozen independent region IDs, not three adjacent windows within one marker, repeat copies or three references containing the same homolog. Annotated overlapping/homologous loci count once. For large unannotated genomes, default to nonoverlapping 10-kb canonical bins, merging bins intersecting the same recognized marker/repeat locus; for small viral genomes use three nonoverlapping genomic thirds, excluding shared/repetitive-only support. The builder stores this region policy per target. Marker/organelle references may use nonoverlapping windows but remain `marker_signal`. Segments have their own coordinates and cannot be double counted as organisms.
3. Qualifying genomic fragments must include resolution-specific support. Five reads in a conserved gene cannot establish a species. If no distinguishing sequence exists, use the validated group/complex target.
4. Signals below these support gates remain `candidate_signal` or `ambiguous_signal`, available in the report. They are not discarded as biologically absent.
5. Single-locus/mitochondrial/ITS-only findings remain `marker_signal` under this default, even when their species attribution looks persuasive. Promote them to a target-specific supported marker assay only after appropriate specificity testing.
6. A novel/divergent hit below 95% identity may be retained by the rescue route at a broader rank; it cannot silently bypass the exact-species rule.
7. Reads shorter than the policy's minimum, unusual platforms or incompatible assay preparations require a separate tested policy. Do not call those targets negative using an inapplicable rule.
8. Severe sample/route QC failure produces `not_assessed`. A low-depth but technically completed research analysis can report candidate/supported findings and “not detected under the research rule,” with **sensitivity unvalidated**; it cannot claim validated negatives.
9. Clinical validation may justify different thresholds. Version any change and re-run the regression set. Never tune on SAMPLE2 versus the healthy young donors to make expected findings appear.

Required dependency-free decision kernel, to be implemented directly or with behaviorally identical code. `NormalizedEvidence` is the explicit adapter between alignment evidence and the final result schema: `assay_eligibility` and `analysis_status` use the same enums as the result; counts have already undergone the alignment, duplicate-family, independent-region and reference-specificity rules above. Protein-only homology never populates these nucleotide-support counts. The two policy instances map directly to the corresponding nested YAML sections.

```python
from dataclasses import dataclass, replace

@dataclass(frozen=True)
class SupportPolicy:
    min_distinct_fragments: int = 5
    min_regions: int = 3
    min_informative_bases: int = 300

@dataclass(frozen=True)
class NormalizedEvidence:
    assay_eligibility: str = "unknown"
    analysis_status: str = "not_run"
    reference_usable: bool = False
    qc_pass: bool = False
    qualifying_fragments: int = 0
    nonoverlapping_informative_regions: int = 0
    informative_bases_covered: int = 0
    resolution_specific_support: bool = False
    only_unresolved_taxonomic_support: bool = False
    any_candidate_signal: bool = False
    marker_or_organelle_only: bool = False
    protein_only: bool = False

def sequence_status(e: NormalizedEvidence, policy: SupportPolicy) -> str:
    if e.assay_eligibility not in {"eligible", "ineligible", "unknown"}:
        raise ValueError("Invalid assay eligibility")
    if e.analysis_status not in {"completed", "not_run", "failed", "qc_failed"}:
        raise ValueError("Invalid analysis status")
    counts = (e.qualifying_fragments, e.nonoverlapping_informative_regions,
              e.informative_bases_covered)
    if any(type(x) is not int or x < 0 for x in counts):
        raise ValueError("Support counts must be nonnegative integers")
    if (e.assay_eligibility != "eligible" or e.analysis_status != "completed"
            or not e.reference_usable or not e.qc_pass):
        return "not_assessed"
    if e.protein_only:
        return "candidate_signal"  # no nucleotide-negative or species call
    if e.only_unresolved_taxonomic_support:
        return "ambiguous_signal"
    if e.qualifying_fragments == 0:
        return "candidate_signal" if e.any_candidate_signal else "not_detected"
    supported = (
        e.qualifying_fragments >= policy.min_distinct_fragments
        and e.nonoverlapping_informative_regions >= policy.min_regions
        and e.informative_bases_covered >= policy.min_informative_bases
        and e.resolution_specific_support
    )
    if not supported:
        return "candidate_signal"
    return "marker_signal" if e.marker_or_organelle_only else "supported_sequence"

def reference_contract_checks():
    p = SupportPolicy()
    ready = NormalizedEvidence("eligible", "completed", True, True)
    strong = replace(ready, qualifying_fragments=5,
                     nonoverlapping_informative_regions=3,
                     informative_bases_covered=300,
                     resolution_specific_support=True)
    cases = [
        (NormalizedEvidence(), "not_assessed"),
        (ready, "not_detected"),
        (strong, "supported_sequence"),
        (replace(strong, assay_eligibility="unknown"), "not_assessed"),
        (replace(strong, assay_eligibility="ineligible"), "not_assessed"),
        (replace(strong, analysis_status="failed"), "not_assessed"),
        (replace(strong, reference_usable=False), "not_assessed"),
        (replace(strong, qc_pass=False), "not_assessed"),
        (replace(strong, qualifying_fragments=4), "candidate_signal"),
        (replace(strong, nonoverlapping_informative_regions=1), "candidate_signal"),
        (replace(strong, informative_bases_covered=299), "candidate_signal"),
        (replace(strong, resolution_specific_support=False), "candidate_signal"),
        (replace(strong, marker_or_organelle_only=True), "marker_signal"),
        (replace(strong, only_unresolved_taxonomic_support=True), "ambiguous_signal"),
        (replace(strong, protein_only=True), "candidate_signal"),
        (replace(ready, any_candidate_signal=True), "candidate_signal"),
    ]
    for evidence, expected in cases:
        assert sequence_status(evidence, p) == expected
    try:
        sequence_status(replace(ready, qualifying_fragments=-1), p)
    except ValueError:
        pass
    else:
        raise AssertionError("Negative counts accepted")
    return len(cases) + 1

if __name__ == "__main__":
    print(f"{reference_contract_checks()} reference contract checks passed")
```

This kernel intentionally does not output infection or treatment decisions. Contamination flags, validation scope and clinical interpretation are independent mandatory fields. A confirmed laboratory artifact remains visible as an artifact record and must be excluded from the participant's supported-organism count.

### 5.5 Quantification

Report observed fragment counts, informative support and normalized fragments per million QC-passed nonhost fragments. Also retain total QC-passed fragments including host for audit and cross-run normalization checks. An RNA library has its own denominator; do not pool DNA/RNA relative abundances.

Do not convert read fractions, Bracken relative abundance, marker coverage or EukDetect normalized abundance into cells/g, worms/g, viable organisms, infection severity or a percent probability of infection. Absolute load requires an appropriate calibrated spike-in/extraction model or independent quantification. A mitochondrial/rDNA count cannot be compared directly with single-copy bacterial genes as an organism-load ratio.

### 5.6 Contamination and cross-sample leakage

For newly generated libraries require extraction blanks, library blanks, appropriate positive controls and batch/run metadata. Use unique dual-index information where available. Evaluate candidate sequences in controls and other samples with normalized counts and coverage, not an arbitrary rule that any blank hit deletes every sample hit.

Implement a validated batch-background rule per assay. Until available, show control comparisons and flag plausible contamination for review; do not create false mathematical certainty with an unvalidated universal sample:blank ratio. Optional CroCoDeEL can supplement batch screening, especially on compatible profiles, but cannot reconstruct missing blanks or distinguish every natural shared strain from contamination.

For historical samples lacking controls: run the analysis and show `controls_unavailable`. Do not mark contamination checks passed. Reproducible low-level findings may warrant a new independent specimen/library and clinical confirmation; do not claim a technical rerun of the same reads is biological confirmation.

## 6. Organism, determinant and host-linkage rules

### 6.1 Keep evidence layers separate

Each determinant record includes source accession/allele, gene family, aligned identity, covered fraction, discriminating positions, intact/partial status, contig ID, host-linkage status and model provenance. A sequence match to a ubiquitous housekeeping, adhesion, iron-acquisition or secretion gene is not sufficient to label an organism pathogenic.

`host_linkage` enum:

- `unlinked`: found in the same community only.
- `read_pair_supported`: a high-quality pair spans determinant/context with validated uniqueness; weaker than a reconstructed locus.
- `contig_supported`: determinant on a well-supported nonchimeric contig with taxonomically informative flanks.
- `genome_supported`: coherent validated assembly/bin plus locus support; plasmids and mixed strains remain separate concerns.
- `ambiguous`: competing host assignments or mobile-element context prevents attribution.

Only appropriate physical evidence can support “organism carrying determinant.” Same-sample co-occurrence, correlated abundance, a generic ARG-host k-mer prediction or membership in an uncertain MAG is insufficient. Do not hide unlinked toxin/AMR genes: report them as independent genetic findings.

### 6.2 Gene screening and assembly

Run sensitive curated nucleotide-reference screening directly on reads; retain partial candidates. Confirm full genes using source-specific identity/coverage/allele models. For a provisional generic homolog screen use ≥90% covered reference length and ≥95% identity as a **candidate-family** filter only; known AMR/function models override this and may require exact substitutions or different validated thresholds. Never label all similar proteins as resistance or all partial matches as intact toxins.

Assemble sufficiently supported reads without requiring a complete MAG. Run AMRFinderPlus on assembled nucleotide/protein inputs; it is not a FASTQ program. Use `--plus` for its additional selected virulence/stress catalog. Run mixed assemblies without falsely declaring one organism. Apply organism-specific point-mutation models only to confidently assigned appropriate sequence. Preserve partial, frameshifted, internal-stop and contig-end evidence.

Kleborate and StxTyper are specialist assembly tools. Their successful execution does not prove a mixed stool assembly is a single strain. Record mixed alleles and unresolved strain composition. Toxin-gene absence is only assessed when locus/strain coverage supports it; lack of reads at low abundance cannot create a “nontoxigenic” classification.

### 6.3 Fungal resistance

Optional FungAMR integration requires exact species, allele, mutation/CNV/promoter mechanism, drug and experimental evidence. Its signed confidence values must be interpreted correctly: +1 is strongest resistance evidence and +8 weakest; negative values encode susceptibility observations. Human-relevant, experimentally supported records must be separated from agricultural and indirect data. Many mechanisms cannot be inferred from a short conserved protein match. Do not output fungal drug susceptibility from pooled stool.

Preserve contradictory resistance/susceptibility observations for the same variant with their organisms, experiments and conditions. A negative confidence value is a source observation, not permission to report clinical susceptibility from stool.

## 7. Catalog compilation rules

Sections 8–11 are **mandatory seed catalogs**. Each semicolon-separated organism is an individual seed, except explicitly named species complexes, pathotypes, viral groups and genus-level expansion rules. Preserve every alias but create one canonical target. Target IDs are stable application identifiers; names/taxids can change via versioned mappings.

1. Expand a species complex only when current references and distinguishing evidence support members; otherwise retain a complex-level target plus member coverage gaps.
2. Expand named genera/viral groups using documented human-infection or human-host relevance and the source references below. Animal-only and environmental relatives remain competitors/background, not automatically human pathogens.
3. Add newly recognized clinically relevant targets by evidence review plus a tested reference update. Broad classification must continue to expose confidently assigned organisms outside the curated list as “other sequence findings” instead of silently discarding them.
4. Every seed produces either a searchable target or an explicit unresolved/reference-gap record. No absent row is allowed. Broad-genome species counts, clinical-target counts, pathotype counts and determinant counts are separate; do not inflate coverage by adding them together as “organisms.”
5. An unsequenced organism remains a documented gap. Use genus/marker resolution when defensible; do not fabricate species-specific detection instructions or substitute a related animal species.
6. Clinical negative display is limited by stool suitability even when genomic references exist. Tissue-only targets remain unexpected-detection watch items, never part of the “infections excluded” denominator.

7. Expand abbreviated genus names within their explicit row context before machine compilation. Store fully spelled accepted names, explicit alias arrays and explicit complex/member relationships. Slash-separated aliases, species complexes and semicolon-separated species are different constructs and must not share a naive splitting rule. In particular, **Entamoeba coli and Escherichia coli are unrelated target names**, not aliases.

## 8. Bacterial target catalog and pathotype rules

### 8.1 Required bacterial seeds

The first column specifies a stable target family; the builder resolves individual names and current synonyms. All taxa remain searchable irrespective of whether their default report class is alarming. Useful clinical anchors are the [CDC diarrheagenic E. coli resource](https://www.cdc.gov/ecoli/php/technical-info/index.html), [CDC travelers' diarrhea chapter](https://www.cdc.gov/yellow-book/hcp/preparing-international-travelers/travelers-diarrhea.html) and [FDA Bad Bug Book](https://www.fda.gov/food/foodborne-pathogens/bad-bug-book-second-edition). Use the additional primary sources below for uncommon organisms and traits.

| Family | Mandatory organisms/targets | Default interpretation and resolution |
|---|---|---|
| Salmonella | Salmonella enterica, all six subspecies and Typhi/Paratyphi/nontyphoidal serovar diversity; S. bongori | Established enteric. Species/group before serovar; no serotyping from insufficient coverage. |
| Shigella/EIEC | Shigella sonnei; S. flexneri; S. dysenteriae; S. boydii; enteroinvasive E. coli | Established enteric/pathotype-dependent. Preserve combined Shigella/EIEC or E. coli/Shigella result when unresolved. |
| Diarrheagenic E. coli | STEC; ETEC; EPEC; EAEC; EIEC; DAEC; hybrid pathotypes | Trait-defined targets, not separate species. DAEC remains heterogeneous/uncertain; require appropriate linked marker sets for strain-level naming. |
| Other Escherichia | E. albertii; E. fergusonii | E. albertii emerging enteropathogen; E. fergusonii conditional opportunist. |
| Core Campylobacter | C. jejuni; C. coli; C. upsaliensis; C. lari | Established enteric; preserve species complexes if unresolved. |
| Expanded Campylobacter | C. fetus; C. hyointestinalis; C. concisus; C. curvus; C. ureolyticus; C. showae | Rare/conditional or uncertain enteric role, species-specific context. |
| Arcobacter aliases | Arcobacter/Aliarcobacter butzleri; A. cryaerophilus; A. skirrowii | Emerging/rare enteric. Resolve current taxonomy rather than duplicate renamed genera. |
| Yersinia | Y. enterocolitica; Y. pseudotuberculosis | Pathogenic lineage/virulence context; not every Y. enterocolitica strain is equivalent. Other Yersinia retained as competitors. |
| Vibrio | V. cholerae; V. parahaemolyticus; V. mimicus; V. fluvialis; V. furnissii; V. vulnificus; V. alginolyticus | Enteric/pathotype-dependent first five; last two have important extraintestinal/seafood context. V. cholerae species alone is not cholera. |
| Aeromonas | A. caviae; A. veronii; A. hydrophila; A. dhakensis; A. jandaei; A. trota; A. schubertii | Enteric/conditional, with variable causality and species discrimination. |
| Other enteric | Plesiomonas shigelloides; Edwardsiella tarda; Providencia alcalifaciens; Laribacter hongkongensis | Rare recognized enteric targets. Use Laribacter spelling and current taxonomy. |
| Toxin-dependent anaerobes | Clostridioides difficile; Clostridium perfringens; enterotoxigenic Bacteroides fragilis | Organism and toxin evidence separate; no disease from abundance alone. |
| Antibiotic-associated Klebsiella | Klebsiella oxytoca complex, including K. oxytoca; K. michiganensis; K. grimontii | Conditional organism detection; toxin-associated hemorrhagic colitis requires additional trait/clinical evidence. |
| Food intoxication/toxicoinfection | Bacillus cereus sensu lato; Staphylococcus aureus; Clostridium botulinum; neurotoxigenic C. baratii; neurotoxigenic C. butyricum | Genetic trait assessment only. Preformed toxin illness may occur without detectable organism DNA; generic B. cereus-group reads cannot establish a particular high-consequence species. |
| Listeria | Listeria monocytogenes | Enteric/conditional; stool DNA is not invasive listeriosis. |
| Intestinal spirochetes | Brachyspira aalborgi; B. pilosicoli | Conditional intestinal association; carriage can be asymptomatic. |
| Gastric Helicobacter | H. pylori; H. suis; H. heilmannii sensu stricto; H. felis; H. bizzozeronii; H. salomonis | Gastric targets with variable stool shedding; preserve non-H. pylori group if unresolved. |
| Enterohepatic Helicobacter | H. cinaedi; H. fennelliae | Conditional/rare; stool cannot establish bloodstream infection. |
| Whipple organism | Tropheryma whipplei | Conditional carriage versus disease. A positive stool sequence is not Whipple disease. |
| Mycobacterial watch | Mycobacterium tuberculosis complex; M. avium complex; curated clinically relevant nontuberculous mycobacteria | Extraintestinal/rare intestinal watch. Swallowed respiratory material is possible; no intestinal-TB or disseminated-MAC diagnosis from stool. |
| K. pneumoniae complex | K. pneumoniae; K. quasipneumoniae; K. variicola; plus K. aerogenes | Opportunistic colonization/reservoir findings. Traits can be added; detection does not require eradication. |
| Enterobacter complex | Enterobacter cloacae; E. hormaechei; E. asburiae; E. kobei; E. ludwigii; E. roggenkampii and resolved E. cloacae-complex members | Conditional opportunists. |
| Other Enterobacterales | Citrobacter freundii complex; C. koseri; Serratia marcescens; Proteus mirabilis; P. vulgaris; P. penneri; Morganella morganii; Providencia rettgeri; P. stuartii; Hafnia alvei; H. paralvei | Conditional opportunists; Hafnia enteric attribution uncertain. |
| Cronobacter | C. sakazakii; C. malonaticus; C. turicensis | Opportunist with neonatal/invasive context; no universal stool infection label. |
| Nonfermenters | Pseudomonas aeruginosa; Acinetobacter baumannii–calcoaceticus complex; Stenotrophomonas maltophilia; Burkholderia cepacia complex; Achromobacter xylosoxidans | Colonization, environmental contamination and infection distinguished. |
| Enterococci | Enterococcus faecalis; E. faecium; E. gallinarum; E. casseliflavus | Conditional opportunists; intrinsic/acquired AMR separated. |
| Other Gram-positive opportunists | Streptococcus agalactiae; S. anginosus group; S. gallolyticus group; Clostridium innocuum; C. tertium; C. septicum | Conditional stool findings, not invasive-disease diagnoses. |
| Anaerobic pathobionts | Bacteroides fragilis group; Fusobacterium nucleatum; F. necrophorum; Parvimonas micra; Peptostreptococcus anaerobius; Ruminococcus/Mediterraneibacter gnavus | Existing ecological/disease associations can cross-link; ordinary detection is not a pathogen red flag. |

Additional source anchors for the uncommon entries: [human Aeromonas study](https://wwwnc.cdc.gov/eid/article/9/5/02-0451_article), [Laribacter genome study](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000416), [human Brachyspira genomics](https://journals.asm.org/doi/10.1128/jb.00272-19), [T. whipplei diagnostic study](https://academic.oup.com/cid/article-abstract/47/5/659/296177), [2025 non-H. pylori surveillance](https://wwwnc.cdc.gov/eid/article/31/6/24-1315_article).

**Archaea:** include them in broad taxonomy and ecological/methane modules. Do not invent a catalog of established human archaeal enteric pathogens or label methanogens infectious pathogens because of a methane/constipation association.

### 8.2 Mandatory pathotype/toxin logic

Genes below refer to curated experimentally characterized families/alleles with exact reference records, not a substring match to an annotation. Negative components require evaluability. Community-level co-occurrence never supplies strain linkage.

| Target | Required evidence and allowed interpretation |
|---|---|
| STEC | stx1/stx2 family evidence, A/B subunit and subtype where resolved. Do not require eae or O157: eae-negative and non-O157 STEC exist. Unlinked stx becomes “Shiga-toxin gene DNA; host unresolved,” including possible phage context. Only credible E. coli linkage supports the organism-trait statement. |
| ETEC | Human-relevant eltA/eltB and/or estA/STa alleles. Do not treat every est-like annotation as enterotoxin or infer which coexisting E. coli carries it. |
| EPEC | eae-associated evidence; stx assessed independently. Typical EPEC requires appropriate bfpA/EAF context. If stx absence is unassessable, output EPEC-associated markers, not a definitive EPEC or atypical-EPEC strain. |
| EAEC | aggR, aatA, aaiC and supported adhesin/regulon context. astA/EAST1 alone is insufficient. Preserve heterogeneous atypical phenotypes and mixed-strain uncertainty. |
| Shigella/EIEC | ipaH family and invasion-context evidence; default combined group unless validated discriminatory markers support narrower resolution. ipaH is not uniquely Shigella. |
| DAEC | afa/daa/dr adhesin evidence is supportive/exploratory, not universally specific. |
| Hybrid E. coli | Multiple marker sets in one community are not a hybrid strain. Require within-strain linkage. |
| C. difficile | Organism DNA; tcdB and tcdA independently; cdtA/cdtB separately. Major-toxin DNA implies genetic potential, not expressed toxin or CDI. Binary toxin does not substitute for major toxins. |
| C. perfringens | cpe for enterotoxin-associated potential; curated other toxin families separately. plc/cpa alone does not establish food-poisoning potential. |
| Enterotoxigenic B. fragilis | bft-1/bft-2/bft-3 family and B. fragilis linkage. B. fragilis abundance alone is not ETBF. |
| Toxigenic V. cholerae | ctxA/ctxB independently from species; tcpA supportive, not toxin replacement. O1/O139 or other lineage claims require validated typing evidence. |
| V. parahaemolyticus | tdh/trh disease-associated markers; tlh identification-associated, not proof of toxigenicity. Lack of marker reads cannot establish a harmless strain. |
| Pathogenic Yersinia traits | Species/lineage-specific ail, ystA, inv, virF/yadA and plasmid context. Do not use one generic homolog as a universal pathogenicity test. |
| B. cereus group toxins | ces cluster for cereulide capacity; nheABC, hbl operon and cytK alleles for diarrheal-toxin capacity. Report locus completeness and taxonomic ambiguity. |
| Staphylococcal enterotoxins | Curated sea/seb/sec/sed/see and experimentally supported additional enterotoxin families. Distinguish enterotoxin from broad superantigen/adhesion annotations. |
| Botulinum neurotoxin-associated DNA | Curated bont family/allele and relevant cluster context with strict near-neighbor/false-positive checks. Species detection alone and unrelated protease homology are insufficient. Any supported signal needs appropriate clinical confirmation; do not infer expressed toxin, viability or absence of intoxication from DNA results. |
| K. oxytoca-associated toxin cluster | npsA/npsB/thdA and surrounding kleboxymycin/tilimycin/tilivalline locus context. Generic NRPS homology is insufficient. |
| K. pneumoniae virulence-associated loci | iuc, iro, ybt, clb and rmp/rmpA2 context, using appropriate Kleborate analysis when possible. A siderophore gene alone does not establish hypervirulence or invasive disease. |
| H. pylori traits | Identification independent from cagA/vacA. Resolve relevant loci/alleles only with adequate coverage. Absence of cagA does not imply absence of H. pylori or benignity. |
| Colibactin capacity | clb/pks island completeness and host linkage; a clbP-like hit alone is insufficient. Link to existing CRC mechanisms without diagnosing cancer or measuring toxin exposure. |

Implementation references: [FDA BAM E. coli methods](https://www.fda.gov/food/laboratory-methods-food/bam-chapter-4a-diarrheagenic-escherichia-coli), [pathotype marker study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5852171/), [Shigella/EIEC discrimination](https://pmc.ncbi.nlm.nih.gov/articles/PMC8767346/), [ETBF genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC4922554/), [K. oxytoca toxin mechanism](https://www.pnas.org/doi/10.1073/pnas.1819154116), [CDC C. difficile test interpretation](https://www.cdc.gov/c-diff/hcp/diagnosis-testing/index.html).

For the K. oxytoca locus, [MIBiG BGC0000446.5](https://mibig.secondarymetabolites.org/go/BGC0000446.5) links GenBank `HG425356.1`; the record labels completeness “complete” but quality “questionable.” Preserve that qualification and validate against primary evidence; do not silently treat it as a perfect reference.

### 8.3 AMR reporting panel

Screen the complete installed, rights-compatible acquired-resistance catalog, with special display priorities below. These are gene/allele families, not diagnoses or antibiotic recommendations.

| Group | Prioritized families and key constraint |
|---|---|
| Carbapenemases | blaKPC, blaNDM, blaVIM, blaIMP, blaOXA-48-like and additional curated carbapenemase alleles. Not every OXA allele is a carbapenemase. |
| ESBL/AmpC | blaCTX-M; appropriate ESBL subsets of blaSHV/blaTEM; blaCMY/blaDHA and other curated AmpC. Not every TEM/SHV allele is ESBL. |
| Glycopeptides | vanA/vanB and supported others, distinguishing intrinsic species-associated systems. |
| Methicillin | mecA/mecC. S. aureus plus unlinked mecA cannot become MRSA: another staphylococcus may carry the determinant. |
| Colistin | mcr families and appropriate species-specific mechanisms; gene family alone does not establish expression. |
| Quinolones | qnr, aac(6′)-Ib-cr and organism-assigned gyrA/parC mechanisms when position/allele coverage permits. |
| Macrolide/lincosamide | erm, mef, mph, lnu curated alleles. |
| Tetracyclines | tet mechanisms separated into curated efflux, ribosomal protection and enzymatic families. |
| Aminoglycosides | aac/aph/ant alleles; armA/rmt methylases. |
| Sulfonamide/trimethoprim | sul/dfr families. |
| Oxazolidinone/phenicol | cfr, optrA, poxtA and appropriate cat/floR families. |
| Other supported mechanisms | Curated rifamycin/nitroimidazole and other determinants; rpoB and other substitutions interpreted only in the correct organism and reference coordinate system. |

AMR cards state “resistance-associated genetic determinant detected,” plus completeness, host linkage and validation. They must not label the person infected, the entire microbiome resistant, a specific organism resistant without linkage, or a drug clinically effective/ineffective without appropriate phenotypic evidence.

## 9. Parasite catalog, genomic routes and specimen limitations

Reference-route notation: **E2** means exact species/reference representation observed in the audited March 2026 EukDetect2 manifest, not clinical validation. **WB:BioProject** means a genome source verified in the WBPS19 list; these are project IDs, not assembly accessions. **N** means acquire and verify NCBI/primary genome or marker records. ParaRef membership must be computed from its actual accession manifest. No route code makes a target reportable before reference checks.

### 9.1 Intestinal protists and microsporidia

| Target seeds | Route | Interpretation |
|---|---|---|
| Giardia duodenalis (= G. intestinalis, G. lamblia), human assemblages A/B | E2; GiardiaDB; N | Established enteric; aliases count once, assemblage only if resolved. |
| Entamoeba histolytica | E2; AmoebaDB; N | Established enteric; distinguish E. dispar/moshkovskii and other relatives. |
| Cryptosporidium parvum; C. hominis | E2; CryptoDB | Established enteric. |
| C. meleagridis; C. felis; C. canis; C. ubiquitum; C. viatorum; C. muris | E2 | Less common human-infecting species; do not require one of the two common species. |
| C. cuniculus; Cryptosporidium mink genotype | N | Rare human-infecting/associated targets; verify exact reference and host evidence. |
| Cryptosporidium chipmunk genotype I | E2 | Preserve genotype identity; do not invent a species assignment. |
| Cyclospora cayetanensis sensu lato umbrella; reviewed C. cayetanensis/lineage A, C. ashfordi/lineage B, C. henanensis/lineage C | E2 for represented legacy references; N for verified lineage genomes | Established enteric umbrella; current species-boundary uncertainty and marker resolution must remain explicit. |
| Cystoisospora belli (= Isospora belli) | N | Established enteric; E2's C. suis is not coverage. |
| Balantioides coli (= Balantidium coli, Neobalantidium coli) | N; reviewed SSU/ITS and other loci | Established enteric; fish Balantidium is a competitor, not a substitute. |
| Sarcocystis hominis; S. suihominis | N | Intestinal sarcocystosis; distinguish muscle-associated Sarcocystis. |
| Sarcocystis heydorni; S. sigmoideus | N; curated COI markers | Documented human intestinal infections; clinical pathogenicity remains uncertain. No invented WGS accession. |
| Enterocytozoon bieneusi | E2; MicrosporidiaDB; N | Intestinal opportunist; genotype does not alone identify transmission source. |
| Encephalitozoon intestinalis (= Septata intestinalis) | E2 | Intestinal opportunist. |
| Encephalitozoon cuniculi; E. hellem | E2 | Contextual opportunists; stool does not establish systemic involvement. |
| Dientamoeba fragilis | E2 | Detection with variable/disputed clinical significance. |
| Blastocystis spp., curated human subtypes | E2; N | Uncertain enteric role; no automatic parasite-treatment alarm. |
| Entamoeba moshkovskii | E2 | Emerging/uncertain pathogenicity; separate from E. histolytica. |
| Entamoeba bangladeshi | N | Emerging/uncertain significance. |
| Entamoeba nuttalli | E2 | Rare/zoonotic context; review human-relevance evidence before clinical tier promotion. |

Clinical/identity anchors: [Giardia](https://www.cdc.gov/dpdx/giardiasis/index.html), [Cryptosporidium](https://www.cdc.gov/dpdx/cryptosporidiosis/index.html), [amebiasis](https://www.cdc.gov/dpdx/amebiasis/index.html), [balantidiasis](https://www.cdc.gov/dpdx/balantidiasis/index.html), [sarcocystosis](https://www.cdc.gov/dpdx/sarcocystosis/index.html), [microsporidia](https://www.cdc.gov/dpdx/microsporidiosis/index.html), [Dientamoeba](https://www.cdc.gov/dpdx/dientamoeba/index.html), [Blastocystis](https://www.cdc.gov/dpdx/blastocystis/index.html).

**Updated human-species evidence:** a [2025 Sarcocystis study](https://wwwnc.cdc.gov/eid/article/31/3/24-1640_article) identified S. heydorni and S. sigmoideus in microscopy-positive human stools, often with other species. It used a 332-bp COI assay; generated datasets are available on request, not a verified public shotgun cohort. Use reviewed reference-marker accessions and retain species-specific clinical uncertainty. Animal-associated S. cruzi remains a competitor/uncertain finding.

For Cyclospora, the [CDC lineage account](https://www.cdc.gov/advanced-molecular-detection/php/success-stories/cyclospora.html) recognizes the A/B/C names, while a [2024 genome study](https://link.springer.com/article/10.1186/s12864-024-10163-y) discusses uncertainty in species boundaries. Acquire [PRJNA1045665](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1045665) metadata and verify assembly-to-lineage mapping. Shared markers yield the sensu-lato umbrella; never count umbrella and resolved members as separate infections.

### 9.2 Intestinal nematodes

| Target | Verified route / required follow-up | Important distinction |
|---|---|---|
| Ascaris lumbricoides | WB:PRJEB4950 | A. suum/related or hybrid ambiguity retained. |
| Ascaris suum | WB:PRJNA62057; PRJNA80881 | Human zoonotic infection possible; no forced species if complex unresolved. |
| Trichuris trichiura | WB:PRJEB535 | Bacterial-contaminant masks essential. |
| Necator americanus | WB:PRJNA1007425; PRJNA72135 | Core hookworm. |
| Ancylostoma duodenale | WB:PRJNA72581 | Core hookworm. |
| Ancylostoma ceylanicum | WB:PRJNA72583; PRJNA231479 | Human zoonotic hookworm. |
| Ancylostoma caninum | WB:PRJNA72585 | Contextual human intestinal disease, usually no patent human egg shedding. |
| Strongyloides stercoralis | WB:PRJNA930454; PRJEB528 | Low/intermittent larval output; negative stool DNA cannot exclude it. |
| S. fuelleborni subsp. fuelleborni; S. fuelleborni subsp. kellyi | N | Rare human intestinal infections; do not substitute another Strongyloides. |
| Enterobius vermicularis | WB:PRJEB503 | Stool has limited sensitivity; perianal tape is the relevant confirmation approach. |
| Capillaria/Paracapillaria philippinensis | N | Intestinal capillariasis; distinguish hepatic Capillaria/Calodium. |
| Trichostrongylus orientalis; T. colubriformis; T. axei | N | Zoonotic intestinal targets. |
| Oesophagostomum bifurcum; O. aculeatum | N | Human/rare intestinal targets; pig O. dentatum cannot substitute. |
| Ternidens deminutus | N | Rare intestinal target with reference-sufficiency gate. |
| Syphacia obvelata | N | Exceptional human reports; manual-review tier until human relevance/reference specificity is resolved. |

Sources: [Ascaris](https://www.cdc.gov/dpdx/ascariasis/index.html), [hookworm](https://www.cdc.gov/dpdx/hookworm/index.html), [Strongyloides](https://www.cdc.gov/dpdx/strongyloidiasis/index.html), [pinworm](https://www.cdc.gov/dpdx/enterobiasis/index.html), [intestinal capillariasis](https://www.cdc.gov/dpdx/intestinalcapillariasis/index.html), [Trichostrongylus](https://www.cdc.gov/dpdx/trichostrongylosis/index.html), [Oesophagostomum](https://www.cdc.gov/dpdx/oesophagostomiasis/index.html). Do not create an independent Enterobius gregorii disease count without resolving its disputed identity.

### 9.3 Intestinal tapeworms

| Target | Route | Distinction |
|---|---|---|
| Taenia solium | WB:PRJNA170813 | Intestinal taeniasis; stool does not diagnose tissue cysticercosis. |
| Taenia saginata | WB:PRJNA71493 | Intestinal. |
| Taenia asiatica | WB:PRJNA299871; PRJEB532 | Intestinal; related-species ambiguity retained. |
| Hymenolepis/Rodentolepis nana | WB:PRJEB508 | One canonical target. |
| Hymenolepis diminuta | WB:PRJEB30942; PRJEB507 | Intestinal. |
| Dipylidium caninum | N | Intestinal zoonosis. |
| Dibothriocephalus latus (= Diphyllobothrium latum) | WB:PRJEB1206 | Fish tapeworm. |
| Dibothriocephalus nihonkaiense; D. dendriticus | N | Fish tapeworms; preserve historical Diphyllobothrium aliases. |
| Adenocephalus pacificus (= Diphyllobothrium pacificum) | N | Fish tapeworm. |
| Diphyllobothrium stemmacephalum (= D. yonagoense); D. balaenopterae (= Diplogonoporus grandis) | N | Rare human fish-tapeworm records; current synonym review required. |
| Bertiella studeri; B. mucronata | N | Rare intestinal. |
| Inermicapsifer madagascariensis | N | Rare intestinal. |
| Raillietina celebensis | N | Rare intestinal; genus result if species resolution unsupported. |
| Mesocestoides spp. | N | Rare intestinal; require human-species evidence for expansion, no invented species assignment. |

Sources: [taeniasis](https://www.cdc.gov/dpdx/taeniasis/index.html), [hymenolepiasis](https://www.cdc.gov/dpdx/hymenolepiasis/index.html), [fish tapeworms](https://www.cdc.gov/dpdx/diphyllobothriasis/index.html), [CDC parasite index](https://www.cdc.gov/dpdx/az.html).

### 9.4 Intestinal/hepatobiliary flukes and acanthocephalans

| Target | Route | Interpretation |
|---|---|---|
| Schistosoma mansoni | WB:PRJEA36577 | Intestinal egg shedding. |
| S. japonicum | WB:PRJEA34885; PRJNA520774; PRJNA724792 | Intestinal egg shedding. |
| S. mekongi | N | Intestinal egg shedding. |
| S. intercalatum; S. guineensis | WB:PRJEB44434 | Distinguish current taxonomy and shared project records. |
| S. haematobium | WB:PRJNA78265; PRJEB44434 | Primarily urinary; stool-negative result does not screen out urinary schistosomiasis. |
| Clonorchis sinensis | WB:PRJNA386618; PRJDA72781 | Biliary infection with stool eggs. |
| Opisthorchis viverrini | WB:PRJNA222628 | Biliary infection. |
| Opisthorchis felineus | WB:PRJNA413383 | Biliary infection. |
| Fasciola hepatica | WB:PRJEB58756; PRJNA179522 | Eggs mainly in patent/chronic infection; early infection can be stool-negative. |
| Fasciola gigantica | WB:PRJNA230515 | Related species/hybrid ambiguity retained. |
| Fasciolopsis buski | WB:PRJNA284521 | Intestinal fluke. |
| Dicrocoelium dendriticum | WB:PRJEB44434 | Biliary; ingestion of infected liver can cause spurious passage. |
| Dicrocoelium hospes | N | Rare biliary target. |
| Heterophyes heterophyes; H. nocens | N | Intestinal minute flukes. |
| Metagonimus yokogawai; M. miyatai; M. takahashii | N | Intestinal minute flukes. |
| Haplorchis taichui; H. pumilio; H. yokogawai | N; published mitochondrial/ribosomal loci | Intestinal; marker-only route likely for some targets. |
| Centrocestus formosanus; Procerovum varium; Stellantchasmus falcatus | N; published ribosomal loci | Intestinal. |
| Heterophyopsis continua; Pygidiopsis summa | N | Rare intestinal; current taxonomy/sequence verification. |
| Gastrodiscoides hominis | N; published ITS/ribosomal/mitochondrial loci | Whole-genome reference not verified; explicit marker-only gap if applicable. |
| Echinostoma revolutum complex; E. ilocanum; E. cinetorchis; E. macrorchis; E. trivolvis; E. echinatum (= E. lindoense); E. fujianensis | N | Rare intestinal flukes; historical names and species complexes require taxonomic resolution. |
| Isthmiophora hortensis (= Echinostoma hortense); Artyfechinostomum malayanum | N | Rare intestinal flukes; preserve accepted aliases. |
| Echinoparyphium; Acanthoparyphium; Episthmium; Himasthla; Hypoderaeum | N; review-only genus expansion | Sporadic human echinostomid records; require species-level human evidence and usable sequences before adding a clinical species target. |
| Moniliformis moniliformis | N | Rare intestinal acanthocephalan; humans may not shed eggs. |
| Macracanthorhynchus hirudinaceus; M. ingens | N | Rare intestinal acanthocephalans. |
| Acanthocephalus rauschi; Pseudoacanthocephalus bufonis; Corynosoma strumosum; Bolbosoma spp. | N | Exceptional/manual-review candidates; not clinically validated stool targets by inclusion here. |

Sources and sequence leads: [schistosomiasis](https://www.cdc.gov/dpdx/schistosomiasis/index.html), [dicrocoeliasis](https://www.cdc.gov/dpdx/dicrocoeliasis/index.html), [echinostomiasis](https://www.cdc.gov/dpdx/echinostomiasis/index.html), [acanthocephaliasis](https://www.cdc.gov/dpdx/acanthocephaliasis/index.html), [heterophyid primary sequence study](https://doi.org/10.1186/s13071-017-1968-0), [H. taichui mitogenome](https://pubmed.ncbi.nlm.nih.gov/24516279/), [H. pumilio molecular detection](https://wwwnc.cdc.gov/eid/article/28/11/22-0653_article), [Gastrodiscoides molecular characterization](https://pubmed.ncbi.nlm.nih.gov/19198879/), [2026 molecular human case](https://pmc.ncbi.nlm.nih.gov/articles/PMC12926597/).

### 9.5 Tissue-associated and nonstandard stool targets

Search these when references permit, but **never include them in a claim that negative stool sequencing excludes their infections**. Report unexpected supported sequences with source/specimen ambiguity.

| Targets | Why a stool-negative result is not exclusion |
|---|---|
| Toxoplasma gondii | Humans do not ordinarily shed its oocysts; food/environmental sequences can occur. |
| Trichinella spiralis; T. britovi; T. nativa; T. nelsoni; T. murrelli; T. pseudospiralis; T. papuae; T. patagoniensis; T. zimbabwensis; recognized genotypes | Tissue infection is not established by ingested-meat DNA or excluded by stool negativity. |
| Toxocara canis; T. cati; Baylisascaris procyonis | Humans usually harbor larvae rather than egg-shedding adult worms. |
| Echinococcus granulosus sensu lato; E. multilocularis; E. vogeli; E. oligarthrus | Humans are intermediate hosts; stool is not a hydatid/alveolar disease screen. |
| Anisakis simplex complex; A. pegreffii; Pseudoterranova spp. | Foodborne larval GI disease generally lacks a patent human stool-egg cycle. |
| Angiostrongylus cantonensis; A. costaricensis | Human infections generally do not yield ordinary stool larvae/eggs; abdominal disease can be stool-negative. |
| Gnathostoma spp.; Spirometra spp.; Taenia multiceps; T. serialis; Calodium hepaticum (= Capillaria hepatica) | Tissue-associated disease; appropriate non-stool confirmation needed. |
| Paragonimus westermani; P. kellicotti; P. heterotremus and clinically relevant relatives | Pulmonary disease; swallowed sputum can yield stool eggs. Stool detection does not imply intestinal adult worms. |
| Plasmodium; Babesia; Trypanosoma; Leishmania; human filarial parasites | Broad unexpected sequence watch only; no stool diagnostic-negative claim. |
| Anncaliia algerae; A. connori; A. vesicularum; Trachipleistophora hominis; T. anthropophthera; Pleistophora ronneafiei; Vittaforma corneae; Tubulinosema acridophagus; Nosema ocularum; Microsporidium africanum; M. ceylonensis | Extraintestinal/contextual microsporidia; retain site-specific clinical relevance and unresolved taxonomy. |

Life-cycle sources: [Toxoplasma](https://www.cdc.gov/dpdx/toxoplasmosis/index.html), [Trichinella](https://www.cdc.gov/dpdx/trichinellosis/index.html), [Toxocara](https://www.cdc.gov/dpdx/toxocariasis/index.html), [Echinococcus](https://www.cdc.gov/dpdx/echinococcosis/index.html), [Anisakis](https://www.cdc.gov/dpdx/anisakiasis/index.html), [Angiostrongylus](https://www.cdc.gov/dpdx/angiostrongyliasis/index.html), [Paragonimus](https://www.cdc.gov/dpdx/paragonimiasis/index.html).

### 9.6 Mandatory nonpathogenic/animal competitors

Include Entamoeba dispar; Entamoeba coli; Entamoeba hartmanni; Entamoeba polecki; Endolimax nana; Iodamoeba buetschlii; Chilomastix mesnili; Enteromonas hominis; Retortamonas intestinalis; Pentatrichomonas hominis. Detection does not automatically imply a pathogenic parasite. [CDC intestinal amebae](https://www.cdc.gov/dpdx/intestinalamebae/index.html), [nonpathogenic flagellates](https://www.cdc.gov/dpdx/nonpathogenic_flagellates/index.html).

Also include relevant animal Giardia/Cryptosporidium/Cystoisospora/Eimeria, fish ciliates and nonhuman Sarcocystis; Strongyloides ratti/venezuelensis/papillosus, Trichuris suis/muris, Oesophagostomum dentatum, free-living/plant nematodes, neighboring flukes/tapeworms and common food-animal genomes. Use food sequences as competitors, not an indiscriminate pre-filter that discards potential true parasite reads.

### 9.7 Parasite-specific implementation rules

1. Do not copy ParaRef's research rule selecting only the highest-k-mer species per genus. It can suppress genuine coinfection. Use uniquely supported multi-species assignments and a shared genus pool.
2. ParaRef's published k-mer/read/entropy thresholds are reproduction settings for that study, not universal clinical cutoffs. Our provisional research policy and later assay-specific validation remain separately versioned.
3. Nuclear and mitochondrial/rDNA evidence are separate tracks. A marker-only target never reports whole-genome breadth. Do not demand mitochondrial confirmation for microsporidia or other organisms without the relevant conventional organelle genome.
4. EukDetect and CORRAL sharing markers is computational corroboration, not independent specimen/clinical confirmation.
5. Validate egg/oocyst/spore/larval recovery with the matched extraction protocol. Naked DNA spikes do not establish those lysis capabilities.
6. Do not infer infection intensity, worm count or egg count from read abundance without a validated organism/protocol-specific quantitative assay.
7. Confirmation options must match the organism: appropriate stool NAAT/antigen/microscopy; pinworm perianal tape; Strongyloides specialist-selected stool methods and/or serology; tissue parasites require appropriate additional specimens/tests.
8. Nonintestinal or rare uncertain findings must remain visible for technical review without becoming an automatically alarming clinical diagnosis.

## 10. Fungal and other opportunistic-eukaryote catalog

Every organism below is a reference/search target or a documented expansion group. **This is not a list of organisms that normally cause intestinal infection.** Colonization, food passage, skin/oral contamination and invasive disease are separate interpretations. Clinical prioritization can use the [WHO fungal priority list](https://www.who.int/publications/i/item/9789240060241) and [CDC fungal disease catalog](https://www.cdc.gov/fungal/about/types-of-fungal-diseases.html); these across-body-site resources do not validate stool screening for all their members.

| Group | Mandatory seeds / expansion | Default stool interpretation |
|---|---|---|
| Common Candida | Candida albicans; C. dubliniensis; C. tropicalis; C. parapsilosis; C. orthopsilosis; C. metapsilosis | Conditional opportunists/colonization. |
| Glabrata relatives | C. glabrata; C. nivariensis; C. bracarensis and corresponding Nakaseomyces names | Keep these three species distinct; deduplicate only each species' legacy/current synonyms. Use an unresolved related-species group when discrimination is insufficient. |
| Other Candida-related yeasts | C. krusei (= Pichia kudriavzevii); C. lusitaniae (= Clavispora lusitaniae); C. guilliermondii; C. fermentati; C. kefyr (= Kluyveromyces marxianus); C. inconspicua; C. norvegensis; C. rugosa complex | Opportunists; preserve Meyerozyma and other current synonyms. |
| Auris/haemulonii group | Candida/Candidozyma auris; C. haemulonii complex, including clinically recognized closely related taxa | Strict near-neighbor resolution. Supported C. auris findings warrant confirmation/infection-control review, not an invasive-candidiasis diagnosis. |
| Saccharomyces | S. cerevisiae, including boulardii-associated lineage | Food/probiotic/colonization context; species reads do not identify a specific probiotic strain. |
| Rhodotorula | R. mucilaginosa; R. glutinis; R. minuta | Conditional opportunists/environmental context; explicitly close prior omission. |
| Trichosporon relatives | T. asahii; T. inkin; T. asteroides; T. mucoides; curated human-infecting Trichosporon/Cutaneotrichosporon members | Conditional opportunists, with accepted taxonomy and complex resolution. |
| Cryptococcus | C. neoformans complex; C. gattii complex; clinically documented reclassified opportunistic cryptococcal yeasts | Extraintestinal/conditional watch; stool is not a cryptococcosis exclusion test. |
| Geotrichum-related | Geotrichum candidum; Magnusiomyces/Saprochaete/Geotrichum capitatus; Saprochaete clavata | Conditional opportunists, aliases deduplicated. |
| Other yeasts | Wickerhamomyces anomalus; Cyberlindnera fabianii; Kodamaea ohmeri; Debaryomyces hansenii | Conditional/food-associated context. |
| Malassezia | M. restricta; M. globosa; M. furfur; M. sympodialis; M. pachydermatis | Skin-associated/ecological stool finding. Does not diagnose seborrheic dermatitis or fungal gut disease. |
| Additional food/ecological yeasts | Curated Pichia, Kazachstania and other common dietary yeasts | Background/ecology unless a separately supported clinical organism-specific context applies. |
| Aspergillus | A. fumigatus; A. flavus; A. niger complex; A. terreus; A. nidulans; clinically relevant cryptic relatives | Opportunistic mould; stringent food/environmental contamination review. |
| Fusarium | Clinically relevant Fusarium species complexes and current reclassified members | Opportunistic mould; preserve complex-level results. |
| Scedosporium/Lomentospora | Scedosporium apiospermum/boydii complex; Lomentospora prolificans | Opportunistic mould watch. |
| Other hyaline moulds | Purpureocillium lilacinum; Paecilomyces variotii; clinically relevant Acremonium/Sarocladium | Conditional/extraintestinal watch. |
| Mucorales | Rhizopus arrhizus/oryzae; R. microsporus; Mucor circinelloides complex; Lichtheimia corymbifera; L. ramosa; Rhizomucor pusillus; Cunninghamella bertholletiae | Rare invasive GI involvement possible; stool does not establish mucormycosis. |
| Extended Mucorales | Human-infecting Apophysomyces, Saksenaea, Syncephalastrum members | Conditional watch; expand only with documented human relevance. |
| Basidiobolus | Basidiobolus ranarum and verified relevant relatives | Rare GI-invasive fungal context; specimen-specific confirmation. |
| Endemic/systemic fungi | Histoplasma species complex; Talaromyces marneffei; Coccidioides immitis; C. posadasii; Blastomyces dermatitidis; B. gilchristii and relevant relatives; Paracoccidioides brasiliensis/lutzii complexes; pathogenic Sporothrix complex; Emergomyces | Extraintestinal/rare GI watch. Genome inclusion does not imply validated stool sensitivity. |
| Other extraintestinal fungi | Pneumocystis jirovecii; clinically relevant Exophiala, Cladophialophora, Fonsecaea, Alternaria, Curvularia | Unexpected findings; no stool-negative reassurance for pulmonary/skin/systemic infection. |
| Dermatophytes | Human-infecting Trichophyton/Microsporum/Epidermophyton, including T. indotineae | Searchable background/watch; stool is not a skin-fungal diagnostic assay. |
| Nonfungal look-alikes | Pythium insidiosum; clinically documented Prototheca spp. | Pythium is an oomycete, Prototheca are algae. Use `other_eukaryotes`, not fungi or bacteria. |

Specific C. auris interpretation anchor: [CDC clinical overview](https://www.cdc.gov/candida-auris/hcp/clinical-overview/index.html). The card can say “sequence finding warrants confirmation”; it cannot equate colonization with invasive infection or automatically recommend an antifungal.

Fungal read counts and ITS findings must carry extraction/reference limitations. [Controlled-diet primary work](https://pmc.ncbi.nlm.nih.gov/articles/PMC5874442/) illustrates that food and oral sources can account for detected stool fungi in a studied cohort; it does not establish that all gut fungi are transient in everyone. Store recent probiotic/fermented-food exposure as optional context, never as a reason to erase a well-supported result.

**Do not report:** “Candida overgrowth” from relative abundance alone, fungal viability, invasive mycosis, a mycotoxin measurement, or a resistance phenotype from an unlinked homolog. If the existing system contains fungal ecological scores, keep them separate from this sequence-detection interpretation.

## 11. Viral target catalog and DNA/RNA separation

Resolve clinical aliases against the frozen NCBI taxonomy and ICTV VMR. Do not hand-invent new ICTV binomials. Include representative strain/genogroup diversity, avoid redundant millions of near-identical records, and keep human-host evidence distinct from human-sample isolation.

### 11.1 Required viral targets

| Class | Mandatory groups/targets | Assay and interpretation |
|---|---|---|
| Established enteric DNA | Human adenovirus 40; human adenovirus 41 | DNA-compatible; distinguish these from other adenoviruses and report type only if resolved. |
| Other human adenoviruses | Human adenovirus species/types beyond F40/F41 | Shedding and conditional disease context; not automatically equivalent to enteric F40/F41 gastroenteritis. |
| Established enteric RNA | Human noroviruses, including documented human-infecting genogroups/variants; human sapoviruses; rotavirus A/B/C and other documented human-associated lineages; classical human astroviruses | Requires compatible RNA/RT assay. DNA-only results explicitly not assessed. |
| Additional enteric RNA | Aichi virus | RNA/RT eligible; clinical attribution remains contextual. |
| Enteroviruses | Human enteroviruses including polioviruses, coxsackieviruses, echoviruses and numbered enteroviruses | RNA/RT; stool shedding does not establish systemic disease. Polio-type unexpected findings need appropriate expert/laboratory confirmation. |
| Parechoviruses | Human parechoviruses | RNA/RT; age and clinical syndrome matter. |
| Hepatitis viruses shed in stool | Hepatitis A; human-infecting hepatitis E lineages | RNA/RT; stool result is not a complete hepatitis diagnosis. |
| Respiratory/systemic RNA shed in stool | SARS-CoV-2 and other supported human-virus sequences surfaced by broad search | RNA/RT and specimen limitations; no inference of infectiousness or cause of GI symptoms. |
| Herpesvirus DNA | CMV; HSV-1; HSV-2; EBV; VZV; HHV-6A; HHV-6B; HHV-7; HHV-8 | Conditional/extraintestinal findings; no automatic reactivation, colitis or antiviral recommendation. |
| Other human DNA viruses | BK/JC and other human polyomaviruses; human papillomaviruses; human bocaviruses; parvovirus B19 and curated human parvoviruses | Conditional/incidental/watch; not all established causes of diarrhea. |
| Uncertain stool-disease attribution | Astrovirus MLB/VA/HMO lineages; saliviruses; cosaviruses; Saffold cardioviruses; bufaviruses; tusaviruses; cutaviruses; other human stool-associated parvoviruses; anelloviruses | Molecule-specific eligibility, uncertain-enteric class. No alarm from stool association alone. |
| Human torovirus-like findings | Curated candidate references only when host/taxonomy can be verified | Uncertain-enteric/host class; reference/attribution gap if unresolved. |
| Uncertain host/background | Picobirnaviruses; smacoviruses; CRESS-DNA/circovirus-like findings; fungal/protist viruses | Sequence finding does not establish human infection. |
| Phage/background | Bacterial/archaeal viruses, including MetaPhlAn VSC outputs; plant/food viruses such as pepper mild mottle virus | Ecological/background section. Never count all viral contigs as human pathogens. |
| Non-enteric systemic watch | HIV/HTLV; HBV/HCV and additional curated human viruses | Broad unexpected-sequence search where molecule permits, not routine stool diagnostic-negative targets. Proviral DNA does not make stool a validated HIV test. |

Viral anchors: [CDC adenovirus](https://www.cdc.gov/adenovirus/hcp/outbreaks/index.html), [CDC HAV diagnostics](https://www.cdc.gov/hepatitis-a/hcp/diagnosis-testing/index.html), [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/), [ICTV VMR](https://ictv.global/vmr). The broad registry must use actual molecule/segment annotations rather than infer them from the row heading.

### 11.2 Viral implementation details

1. Adenovirus benchmark references include `NC_001454.1` for F40 and `DQ315364.2` for F41 in INTEGRATE. Use them as validated acquisition seeds, then include relevant diversity and non-F competitors; two references alone are not complete adenovirus coverage.
2. Segment-aware bundles must connect all segments of a virus. Detecting three rotavirus segments is one virus-group finding, not three pathogens. Partial segment recovery can support a limited sequence result but cannot imply a complete genome or resolved strain/reassortant.
3. Preserve recent live oral vaccine context. Rotavirus vaccine shedding is documented; infer vaccine versus wild-type only with adequate discriminatory evidence. [Primary 2024 cohort](https://academic.oup.com/jid/article/230/3/754/7603816).
4. EBV/CMV DNA alone cannot be labeled reactivation or end-organ disease. [CDC EBV testing](https://www.cdc.gov/epstein-barr/php/laboratories/index.html), [NIH CMV guidance](https://clinicalinfo.hiv.gov/en/guidelines/hiv-clinical-guidelines-adult-and-adolescent-opportunistic-infections/cytomegalovirus).
5. Picobirnavirus is not a default human-pathogen result: bacterial-lysin evidence supports microbial hosts for at least some members. [Primary study](https://www.pnas.org/doi/10.1073/pnas.2309647120).
6. A virus protein/HMM match, geNomad classification or CheckV completeness estimate does not prove human tropism. Keep host certainty and clinical attribution fields mandatory.
7. Distinguish library absence from bioinformatic absence. Additional software cannot recover norovirus RNA that the original DNA-only extraction/library did not sequence.
8. Provide a future RNA-input path using documented RNA preservation, extraction, RT and library selection appropriate to the target panel. Software support can ship before that assay is operational; actual reports must show its availability honestly.
9. Do not equate viral DNA/RNA, assembled genomes or high read counts with infectious virus, transmission risk or need for treatment.

### 11.3 RNA branch implementation boundary

Implement a separate `research_rna_virus_v1` calling profile for RNA-virus nucleotide sequences. Initial engineering support gates may use the same five conservative fragments, three independent viral genomic regions and 300 informative covered bases as the DNA research template, but **must be benchmarked and locked separately** on the RNA portion of INTEGRATE before activation. Region/segment coordinates, controls and library selection must be RNA appropriate. No RNA route may inherit `research_dna_v1` by filename or silently claim validated genomic coverage.

Until that profile is installed and its computational tests pass, uploaded RNA libraries return `not_assessed` for supported/negative RNA-virus calls; discovery candidates remain reviewable. For RNA-only libraries, bacterial/fungal/parasite and DNA-virus transcripts are optional transcript-associated candidates, not genomic-DNA presence/absence screens under this profile. A separate transcript-specific assay/profile would be required to expand those claims. An RNA-preserving total-nucleic-acid protocol can use both branches only where its documented preparation supports both and each branch's computational profile is installed; this still does not establish clinical sensitivity.

This separation is an implementation requirement, not a request to postpone the DNA release. Ship DNA-compatible screening and complete RNA capability/status plumbing; activate the RNA research route once its defined benchmark is passed. New wet-lab RNA preparation is required for DNA-only legacy samples to gain actual RNA-virus coverage.

## 12. Coverage compiler, commands and internal interfaces

### 12.1 Normalized seed records

The coding agent must convert catalog rows into machine-readable records in the repository. The seed compiler must be deterministic and must not require a language model during sample analysis. For groups such as “human noroviruses” or a species complex, retain the explicit group record and separately resolve reviewed members; do not split prose arbitrarily.

Required compilation output for every seed:

```json
{
  "seed_id": "helminth.ascaris_lumbricoides",
  "source_section": "9.2",
  "requested_name": "Ascaris lumbricoides",
  "accepted_taxon": null,
  "aliases": [],
  "resolution_state": "pending_resolution",
  "candidate_sources": ["WBPS19:PRJEB4950", "ParaRef", "NCBI"],
  "resolved_accessions": [],
  "usable_reference_types": [],
  "eligible_assays": ["DNA"],
  "interpretation_class": "established_enteric",
  "stool_suitability": "intestinal_shedding",
  "negative_claim": "not_detected_under_applied_rule_only"
}
```

Resolution and acquisition happen during reference builds; sample runs consume a frozen installed manifest. The example's pending fields must become concrete accessions/IDs or explicit gaps before release. Do not leave them as invisible TODOs. An unavailable rare organism does not block usable targets, but its missing state must appear in the coverage ledger and release notes.

### 12.2 Required application commands

Implement these semantics within the existing language/framework; names can be adapted to established repository conventions. Document exact implemented equivalents in the repository README and PR, with no separate handoff needed.

| Command | Behavior and exit condition |
|---|---|
| `pathogens audit-existing --repo PATH` | Inspect pipeline and report seams; emit actual versions/database/flags and input availability. Does not claim to have run sequencing. |
| `pathogens refs resolve --catalog VERSION` | Resolve aliases, taxids, assemblies and rights-compatible sources; emit all coverage gaps. Fail only malformed/ambiguous build inputs; retain explicit unresolved target records. |
| `pathogens refs build --lock LOCK` | Download/checksum/mask/index fixed assets; resumable; no mutable network lookup during later sample calls. |
| `pathogens refs audit --lock LOCK` | Validate every seed against actual sequence/index coverage, near-neighbor sets and artifact hashes; fail false coverage claims. |
| `pathogens run --sample MANIFEST --refs LOCK` | Execute all eligible routes, preserve per-route errors and evidence, produce canonical result JSON in the application's normal output tree. |
| `pathogens validate --suite SUITE --refs LOCK` | Run held-out computational and available clinical-reference fixtures; keep truth unknown unless independently supplied. |
| `pathogens report --result RESULT` | Render the integrated Pathogens section from canonical data, including full target coverage and assay exclusions. |
| `pathogens compare --old RESULT --new RESULT` | Explain differences from inputs, database, software or rule changes; never silently relabel a historical report. |

Command arguments are data, not shell snippets. Use subprocess argument arrays, validate paths, and reject malformed accession/URL inputs. Sequence databases are untrusted scientific inputs; downloaded scripts must not execute merely because they arrived with a FASTA archive.

### 12.3 Required interfaces

```text
resolve_target(seed, taxonomy_snapshot, source_metadata) -> TargetRecord | CoverageGap
build_reference_bundle(resolved_targets, source_lock) -> ReferenceBundle
check_assay_eligibility(assay, target) -> EligibilityRecord
discover_candidates(reads, broad_bundle, sensitive_bundle) -> CandidateSet
verify_candidates(reads, candidates, competitive_bundle) -> AlignmentEvidence[]
screen_determinants(reads, assemblies, determinant_bundle) -> DeterminantEvidence[]
infer_host_linkage(determinants, alignments, assemblies) -> LinkageRecord[]
evaluate_sequence(evidence, calling_profile) -> SequenceResult
evaluate_contamination(result, controls, batch_context) -> ContaminationRecord
interpret_finding(sequence_result, target, determinants, clinical_context) -> ReportCard
compile_coverage(all_targets, assay, reference_lock, completed_routes) -> CoverageLedger
```

Do not implement a single `is_bad(taxon)` boolean. Clinical interpretation is a typed record with context, evidence and limitations. Likewise, `screened=true` cannot replace assay eligibility + actual reference coverage + successful route completion.

### 12.4 Determinant result states

Use `determinant_status`: `supported_intact_sequence`, `supported_partial_sequence`, `candidate_homolog`, `ambiguous_allele`, `not_detected`, `not_assessed`. Separately record `required_discriminating_positions_assessed`, `locus_completeness`, `host_linkage`, `gene_function_evidence` and `clinical_phenotype_inference`.

Missing/disabled licensed databases, no curated allele, inadequate discriminating-position coverage, a mixed-strain assembly or missing locus context must never yield “toxin absent,” “nontoxigenic,” “nonpathogenic” or “susceptible.” An organism may be supported while its toxin genotype is not assessed. An unlinked toxin/AMR determinant may be supported while its organism is unresolved.

## 13. Report and API requirements

### 13.1 Front-loaded visual summary

Add the Pathogens summary to the report's existing early visual-results pages, before long explanatory sections. Show:

1. **Assay capability:** “Stool DNA sequence screen” or the actual supported DNA/RNA combination; sample date, analysis date and reference version.
2. **Findings requiring attention:** supported sequence findings of established/rare enteric pathogens or important determinants, with technical qualification and confirmation context.
3. **Opportunistic organisms:** supported conditional organisms without labeling colonization as infection.
4. **Uncertain/low-support findings:** count and readable link to all candidate, marker-only, ambiguous and contamination-qualified records. Do not hide them behind a top-N cutoff.
5. **Coverage:** counts of unique catalog targets assessed, not assessed by this assay, missing usable references, and failed/not-run routes. Keep organism targets separate from determinant and subtype counts.
6. **Prominent limit:** “Not detected does not rule out infection. This result depends on the specimen, assay, reference coverage and sequencing sensitivity.” For DNA-only samples additionally show “RNA-virus screening was not performed by this DNA assay.”

Use neutral blue/gray for assessment states, amber for uncertain/confirmation findings and text/icons in addition to color. Do not use a green shield, a “100% safe” circle, an infection-free score or a donor rank. The exact counts come from the coverage ledger, never a hardcoded marketing number.

### 13.2 Group ordering and details

Detailed section order: established bacterial pathogens/pathotypes; protozoa; worms; microsporidia; fungi/opportunists; DNA viruses; RNA viruses; other/uncertain viral findings; toxin/virulence determinants; AMR determinants; complete assessed-target ledger and technical coverage.

Each finding card includes:

- Accepted name, aliases where useful, supported taxonomic rank and clinical class.
- Sequence status with a plain-language statement: for example, “DNA sequences consistent with Giardia detected under the research rule,” or “low-level signal; species assignment unresolved.”
- What the finding can mean: infection-compatible organism, conditional colonizer, uncertain association, dietary/background or nonstandard specimen.
- Evidence: distinct fragment count, supported independent loci, coverage type, near-neighbor ambiguity, reference and calling-policy versions. Detailed alignments belong in an expandable technical panel, not the main patient narrative.
- Separate toxin/pathotype/AMR results and physical linkage status.
- Validation scope, controls status and relevant limitations, including sample/assay/reference gaps.
- Organism-appropriate confirmation options for a qualified clinician/laboratory to choose; primary source links.

Provide all findings and all catalog coverage rows across as many pages as needed. PDF output must repeat table headers, wrap names without clipping, keep status labels readable and preserve live sources. HTML can offer filters/search, but the PDF cannot silently omit records hidden by an interactive filter.

### 13.3 Required wording examples

| Situation | Approved meaning |
|---|---|
| Supported Giardia sequence | “Giardia-associated DNA detected. This sequence finding may warrant confirmation with an appropriate clinical stool test; it does not by itself establish the cause of symptoms.” |
| Strongyloides no supported signal | “No supported Strongyloides sequence detected in this sample. Stool sampling and low/intermittent shedding can miss infection; clinical suspicion may require other testing.” |
| DNA-only norovirus row | “Not assessed: this library was prepared for DNA, while norovirus has an RNA genome.” |
| Candida detected | “Candida DNA detected. Stool carriage or food/oral sources may occur; this result does not diagnose invasive fungal infection or prove a need for antifungal treatment.” |
| stx with unknown host | “Shiga-toxin-associated gene DNA detected; the carrying organism/strain is unresolved.” |
| mecA with S. aureus but no linkage | “mecA DNA detected; host unresolved. This does not establish MRSA.” |
| Parasite reference gap | “Not assessed at species level: no usable validated-discrimination reference is installed. Broader/marker findings, if any, are shown separately.” |
| Species resolved, clinical LOD absent | “Sequence support meets the documented research rule. Clinical sensitivity and limit of detection have not been established for this target and assay.” |
| No supported catalog findings | “No supported sequence findings among the assessed targets. This is not clearance of infection; see unassessed targets and assay limitations.” |

Do not state that the amount of stool is the only limitation. Extraction, shedding, sequencing depth, library chemistry, divergent strains, reference contamination and taxonomic ambiguity also matter. Do not imply that sequencing more reads always solves these problems.

### 13.4 Existing intervention and FMT integration

Link appropriate general evidence/confirmation context into the existing intervention architecture. A pathogen sequence, Candida abundance, AMR gene or unconfirmed toxin locus must not automatically trigger antibiotics, antifungals, antiparasitics, antimicrobial supplements or “eradication” advice. This release does not add dosing or treatment selection.

FMT research views may expose qualified findings and outstanding confirmation needs. They cannot emit donor clearance, “safe donor,” pass/fail, or a combined pathogen-negative safety score. Existing required clinical donor testing remains separate; this broad research screen is supplementary.

Historical sample labels are not truth: SAMPLE2's reported long-COVID status does not establish a particular pathogen, and visibly healthy young donors are not certified pathogen-negative controls. Never train or tune the screen to separate those five samples.

## 14. Validation data, experiments and release evidence

### 14.1 Public data acquisition inventory

| Dataset/resource | Exact entry point | Validation purpose and limitations |
|---|---|---|
| Parks 2026 clinical stool | [PRJNA1200893](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1200893), [paper](https://doi.org/10.3389/fcimb.2026.1759322) | Clinical comparison for the actually tested targets and sample-level reference results. Do not label all targets in every specimen negative merely because they were not tested. |
| INTEGRATE 2025 | [PRJEB62473](https://www.ebi.ac.uk/ena/browser/view/PRJEB62473), [paper](https://doi.org/10.1186/s13073-025-01478-w) | Paired DNA/RNA with diagnostic metadata; principal viral assay-eligibility benchmark. |
| Haugum 2025 | [PRJNA1218764](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1218764), [paper](https://doi.org/10.1371/journal.pone.0331288) | Clinical, spiked and control stool; rare-read and assembly-loss regression. |
| Fungal mocks | [Rounge-lab mock_mycobiome](https://github.com/Rounge-lab/mock_mycobiome), [paper](https://doi.org/10.1186/s40168-025-02048-3) | Reproducible synthetic communities; independent database inclusivity/specificity challenges, not clinical/extraction LOD. |
| Fungal extraction/diet study | [PRJNA925700](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA925700), [paper](https://doi.org/10.1186/s40168-023-01693-w), [code](https://github.com/ManichanhLab/LongitudinalMycobiomeWithDiet) | 143 preparations across repeated/enriched/unenriched stools, not 143 independent patients. Test preparation bias and duplicate participant handling. |
| Global soil-transmitted-helminth diversity | [PRJEB90452](https://www.ebi.ac.uk/ena/browser/view/PRJEB90452), [2025 study](https://doi.org/10.1038/s41467-025-61687-0) | Geographic inclusivity/marker-diversity challenge; not automatically a matched stool clinical-sensitivity cohort. |
| Parasite genome skimming | [2023 primary paper](https://doi.org/10.1016/j.ijpara.2022.12.002), [author repository](https://www.repository.cam.ac.uk/handle/1810/344263) | Feasibility and reference/marker comparisons. Resolve exact raw accessions and specimen truth from paper metadata; accession not verified here. |
| Parasite capture | [2025 hybrid-capture paper](https://doi.org/10.1111/1755-0998.70005) | Optional future wet-lab sensitivity extension, not a transform that can add missing molecules to existing FASTQs. |
| Hunting for Helminths | [March 2026 preprint](https://doi.org/10.64898/2026.03.09.710549) | Emerging mitochondrial/shotgun approach. Full methods/accession not independently recovered in this research; do not copy unverified sensitivity claims. |
| NCBI isolate/AST collections | [Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/isolates) | Gene/allele and species-resolution truth; not a substitute for stool-mixture and extraction validation. |
| ParaRef modern reanalysis | [Paper/data links](https://doi.org/10.1186/s13059-025-03818-w) | Contaminated-reference/animal-food challenges; computational calls alone are not independent infection ground truth. |

For a BioProject linked to SRA, resolve its ENA study alias before assuming the exact project string is accepted by every API. Use the official [ENA portal API](https://www.ebi.ac.uk/ena/portal/api/) to obtain read-run metadata. Request fields supported by the current `read_run` schema, including study/sample/run accessions, library strategy/source, layout, FASTQ URLs and MD5s. Preserve raw responses, selected runs and sample-to-diagnostic mappings. Never infer positive/negative infection status from filenames alone.

Example acquisition pattern for the verified ENA study (URL parameters are ordinary data):

```text
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB62473&result=read_run&fields=run_accession,sample_accession,study_accession,library_strategy,library_source,library_layout,fastq_ftp,fastq_md5&format=tsv&download=true
```

The agent must verify the current response schema, map publication participant/specimen IDs to archive samples and retain unknown labels as unknown. Controlled/unreleased metadata are an explicit gap; do not manufacture clinical labels to finish a benchmark.

### 14.2 Computational validation matrix

Create reusable fixtures from held-out reference strains and public reads, with source accession and random seed recorded. Use negative/background stool matrices, multiple pathogen concentrations, different read lengths/error profiles, human fractions and mixtures. Hold out references by strain/lineage and preferably study/geography, not merely random reads from the same genome used to build the database.

Required challenges:

- Major bacterial enteric taxa/pathotypes and near neighbors; toxin-positive and toxin-negative strains; unlinked mobile genes.
- Major protozoa, intestinal worms and fungal representatives with exact-reference, divergent-reference and reference-absent conditions.
- Two species within the same genus and several kingdoms in one sample, to detect winner-take-all suppression.
- DNA versus RNA viruses, segment mixtures, phages/food viruses, host-uncertain references and assay-inconsistent signals.
- Low-complexity/rDNA-only reads, contaminated parasite genomes, food-animal DNA, human homologs, PCR duplicates and shuffled/mislabeled controls.
- Complete and partial genes, disruptive variants, closely related nonresistance alleles, multiple candidate hosts and assembly chimeras.

Measure per-target/per-depth sensitivity, specificity, precision, false findings per specimen and resolution accuracy. Also measure family-level false-positive burden across the full panel: hundreds of weak independent calls can make a misleading report even if individual classifiers look accurate. Report confidence intervals and test denominators. No global AUROC or average sensitivity can substitute for weak target-specific results.

Default research-rule values can be revised only against held-out results and a documented trade-off, with a new policy version. Do not optimize on the final validation cohort, the same reference genome used to construct informative regions, or the five user samples.

### 14.3 End-to-end analytical validation

Before publishing target-specific clinical LOD or sensitivity claims, evaluate the actual collection/preservative/extraction/library/sequencer/pipeline combination. Include stool-matrix positives and controls, difficult-to-lyse organisms/reference materials where relevant, dilution series, multiple donors/matrices, between-run/lot/operator variation and relevant interferents. Whole-organism recovery cannot be inferred from naked-DNA or in-silico spikes.

Design target-specific dilution/replicate experiments capable of estimating the intended detection probability with uncertainty; store observations and the fitted LOD model. A commonly used 95% detection point is a statistical target, not something established by one successful dilution or a few replicates. Do not label “19/20 detected” as a precisely established 95% lower-bound sensitivity. Replication and confidence intervals must accompany the estimate.

Absolute quantification needs separate calibration for input mass, recovery, genome/marker copy number and measurement error. A target's cell/egg/spore biology may make genome-equivalent quantities different from viable organism counts.

### 14.4 Clinical comparison and scope

Use appropriate independent clinical reference tests for each organism/syndrome. Preserve discordance; do not assume either sequencing or one comparator is perfect. Interpret discrepant low-load findings with planned adjudication rather than selectively confirming only favorable cases. Separate children/adults, symptomatic/asymptomatic populations, immunocompromise and relevant sample protocols when evaluating claims.

Clinical validation is target + assay + population + intended use. A donor-screening intended use cannot inherit performance measured only in symptomatic diarrhea patients. A positive blood/tissue reference does not guarantee stool sensitivity for the organism.

### 14.5 What ships at each maturity level

| Level | Required evidence | Permitted output |
|---|---|---|
| Search candidate | Correctly installed reference and completed candidate route | Reviewable signal, with unresolved evidence explicit |
| Research sequence support | Tested software contracts, held-out computational specificity/inclusivity and frozen research rule | Supported/marker/ambiguous research sequence findings; clinical LOD/sensitivity absent |
| Analytical assay validation | Matched end-to-end validation, targets and uncertainty documented | Target-specific analytical claims within validated assay scope |
| Clinical assay validation | Appropriate clinical comparison and intended-use review | Only the claims supported for that target, assay and population |

The first software release must include the complete implementation and coverage ledger, not wait indefinitely for every rare organism to obtain clinical validation. Equally, incomplete validation must never be hidden behind a generic “high-confidence pathogen test” label.

## 15. Implementation order and operational requirements

### 15.1 Build sequence

1. Audit real repository and add schemas/statuses, stable IDs, feature flag and migration without breaking existing reports.
2. Compile all seeds into reviewed target/gap records; implement source resolution, acquisition, rights metadata, immutable locks and coverage auditing.
3. Build broad DNA references plus parasite/fungal supplements and target-specific sensitive routes. Add positive/negative reference fixtures immediately.
4. Implement competitive confirmation, duplicate-family counting, target-level ambiguity, independent-locus evidence and provisional research policy.
5. Add determinant screening, assembly context, host linkage and specialist pathotype/AMR modules. Complete the open-resource route first; optional licensed adapters are explicit extensions.
6. Add RNA input/capability plumbing and the separately tested RNA-virus research profile. DNA-only legacy samples retain unassessed RNA-virus rows.
7. Implement report/API integration and complete coverage/technical detail pages; test a specimen with many findings and many gaps, not only an empty report.
8. Run public/computational regression suites, fix false assignments and freeze a release lock. Reprocess available existing samples and label their validation/controls status accurately.
9. Publish the implementation PR with actual versions, measured resource requirements, acceptance results and remaining target-specific gaps. Do not claim this specification itself downloaded all references or clinically validated the pipeline.

### 15.2 Resource management

Database builds and sample analyses are separate jobs. Cache immutable references by content hash, use resumable downloads and preflight actual compressed/expanded/index sizes. Large reference builds require adequate scratch and RAM; measure them on the selected snapshot instead of copying old memory estimates from a tool paper.

Allow resource-tier configurations only as explicitly different analytical profiles. A constrained machine may queue work, use supported disk-backed/chunked processing or execute a declared smaller research profile; it must not silently omit fungi/worms/viruses. Report every skipped target/route. When using multiple indexes, handle taxonomy compatibility and cross-index competing evidence; do not union positive classifications as if each was independently specific.

Assembly is resource-limited and optional for organism presence; determinant/strain linkage may remain unresolved when assembly fails. A task timeout produces a route failure record and resumable work, not a negative result.

### 15.3 Privacy and data integrity

Retain host-containing data only under the existing authorized controlled-storage policy. Do not send raw reads or participant metadata to public sequence services. Public APIs are for reference acquisition and public validation data. Local evidence artifacts should use internal sample IDs, access control and retention rules. Reports must not embed human read sequences.

Source assets, metadata, alignments and API inputs may contain malformed text; validate parsing and escape report content. Never execute downloaded reference-package scripts automatically. Checksum failures invalidate the affected bundle and dependent results, not unrelated existing reports.

## 16. Acceptance tests — required executable repository tests

These are implementation acceptance requirements. The coding agent must implement fixtures/assertions in the real repository; this document does not claim the current application already passes them. The small reference kernel in section 5.4 supplies an independently runnable status-contract baseline.

| ID | Fixture/action | Required assertion |
|---|---|---|
| P001 | Actual repository audit | Installed profiler/database/flags are recorded; supplied audit is not copied as fact. |
| P002 | FASTQs missing, abundance table present | Pathogen route is not run; no inferred negative panel. |
| P003 | DNA-only stool | RNA-only targets are not assessed, not negative. |
| P004 | Unknown library protocol | Unknown does not evaluate truthy as eligible. |
| P005 | RNA reads without installed RNA policy | Candidate review allowed; supported/negative RNA calls not assessed. |
| P006 | Paired DNA/RNA specimen | Separate libraries/denominators, linked specimen identity. |
| P007 | Amplicon library | No inherited shotgun coverage claims. |
| P008 | RNA-virus-like hit in DNA-only library | Assay-inconsistent finding retained; RNA panel remains unassessed. |
| P009 | Minor's eligible stool DNA | Pathogen sequence search runs without adult-reference scoring; clinical age context retained. |
| P010 | Five user samples | Never used as pathogen-positive/negative truth from health labels. |
| P011 | Cystoisospora suis reference only | C. belli remains a reference gap. |
| P012 | Fish Balantidium reference only | No Balantioides coli coverage/call. |
| P013 | Related Sarcocystis reference only | No false S. hominis/suihominis coverage. |
| P014 | EukDetect2 current manifest | Major worm gaps trigger supplemental route, not coverage success. |
| P015 | PlusPF installed | Worm coverage not inferred from database name. |
| P016 | RefSeq missing, valid reviewed GenBank genome present | Target can use the qualified GenBank reference. |
| P017 | No usable reference | Explicit gap row; no not-detected species result. |
| P018 | Taxonomy synonym update | One canonical target; old report identity preserved. |
| P019 | Ambiguous name resolver | No first-match selection; target gap/review state. |
| P020 | Source SHA mismatch | Affected reference rejected; no silent use. |
| P021 | Memory limit | No unreported reduced database or downsampled read fallback. |
| P022 | Multi-index analysis | Competing taxa across indexes are considered before supported calls. |
| P023 | E. coli-only mock | No Trichuris call from contaminated reference regions. |
| P024 | Pig-food mock | No Taenia solium call from pig sequence contamination. |
| P025 | Entamoeba dispar mock | No E. histolytica species call. |
| P026 | C. glabrata/nivariensis/bracarensis references | Three distinct species retained; synonyms deduplicated within species. |
| P027 | Multiple conspecific reference placements | Valid target-level unique evidence not rejected solely for low raw MAPQ. |
| P028 | Equally good other-species placements | Ambiguous/group result, not forced species. |
| P029 | High-count conserved rDNA only | Marker/ambiguous status, not genome-distributed species support. |
| P030 | Three adjacent windows of one locus | Counts as one independent genomic locus. |
| P031 | PCR duplicates with different read names | Cannot satisfy five independent-fragment gate. |
| P032 | Overlapping paired reads | Count one supporting fragment. |
| P033 | Valid UMI library | Use tested UMI collapse; retain molecule-rule provenance. |
| P034 | Five fragments/three independent regions/300 informative bases | Supported research sequence only if specificity/eligibility/QC gates pass. |
| P035 | Same evidence, marker-only reference | Marker signal, not whole-genome support. |
| P036 | Four qualifying fragments | Candidate retained, not supported under default rule. |
| P037 | Protein-only homolog | Candidate; no nucleotide-supported species or negative call. |
| P038 | Divergent hit below exact-species identity gate | Broader candidate can survive; no forced species. |
| P039 | Two uniquely supported species within one genus | Both retained; no winner-per-genus suppression. |
| P040 | Giardia plus Cryptosporidium plus bacterium | All independently supported targets retained. |
| P041 | Assembly failure with supported read evidence | Organism evidence preserved; linkage remains limited. |
| P042 | Empty/failed analysis route | Not assessed, never negative. |
| P043 | Historical controls absent | `controls_unavailable`, not contamination-check passed. |
| P044 | Low blank signal plus strong specimen signal | No automatic blanket deletion; qualified background comparison. |
| P045 | Confirmed technical artifact | Excluded from supported-organism count, artifact evidence retained. |
| P046 | Batch shared strains | No automatic contamination assumption from sharing alone. |
| P047 | High relative E. coli without resolved traits | No pathogenic-E. coli or infection label from abundance. |
| P048 | ipaH evidence | No automatic Shigella species identity. |
| P049 | Unlinked stx plus E. coli | Independent toxin-gene finding, no invented within-strain linkage. |
| P050 | eae-negative, linked stx-positive E. coli | STEC-associated evidence retained. |
| P051 | eae with unevaluable stx absence | No definitive atypical-EPEC classification. |
| P052 | astA alone | No definitive EAEC call. |
| P053 | Two pathotype sets in different strains | No hybrid-strain call. |
| P054 | C. difficile with major-toxin genes | Toxigenic potential distinct from diagnosed CDI. |
| P055 | C. difficile binary toxin only | Does not substitute for major-toxin assessment. |
| P056 | B. fragilis abundance, bft unassessed | No ETBF call. |
| P057 | V. cholerae without assessed ctx | Species result separate from cholera/toxigenicity. |
| P058 | V. parahaemolyticus tlh alone | No toxigenic-strain assertion. |
| P059 | C. perfringens plc/cpa alone | No enterotoxin/food-poisoning confirmation. |
| P060 | Complete toxin locus | Genetic capacity, not measured toxin or disease. |
| P061 | Single nonspecific NRPS hit | No K. oxytoca toxin-cluster claim. |
| P062 | Low-resolution B. cereus-group reads | No forced high-consequence species identity. |
| P063 | S. aureus plus unlinked mecA | No MRSA call. |
| P064 | Generic blaTEM/blaSHV allele | No automatic ESBL classification. |
| P065 | Generic OXA allele | No automatic carbapenemase classification. |
| P066 | Partial/disrupted ARG | Partial/ambiguous state retained; no intact determinant assertion. |
| P067 | Wrong-species mutation model | Cannot assign resistance from an unrelated coordinate/model. |
| P068 | Mobile ARG in uncertain MAG | Host linkage remains unresolved/qualified. |
| P069 | Missing licensed CARD/VFDB adapter | Open-resource analysis completes; unavailable determinant coverage explicit. |
| P070 | AMR evidence in commensal | No infection diagnosis or antibiotic recommendation. |
| P071 | Fungal variant has conflicting evidence | Both directions/context retained; no stool susceptibility verdict. |
| P072 | FungAMR signed confidence | +1 strongest resistance; negative is source susceptibility evidence, not patient phenotype. |
| P073 | C. auris-like sequence | Close relatives assessed; qualified confirmation context, no invasive disease claim. |
| P074 | Saccharomyces species reads | No unsupported boulardii strain identification. |
| P075 | Malassezia stool finding | No seborrheic-dermatitis diagnosis. |
| P076 | Mould DNA | No mycotoxin measurement or inferred toxin exposure. |
| P077 | Strongyloides not detected | Explicit sampling/shedding limitation and appropriate confirmation context. |
| P078 | Enterobius not detected | No parasite clearance; tape-test context preserved. |
| P079 | Taenia solium stool sequence | Does not diagnose cysticercosis. |
| P080 | Tissue-only parasite reference covered, no stool signal | No excluded-infection or ordinary negative-screen claim. |
| P081 | Blastocystis/Dientamoeba finding | Uncertain/variable role retained; no automatic treatment alarm. |
| P082 | Microsporidia detected | One parasite count, not an additional fungal count. |
| P083 | Microsporidian evidence | No mandatory mitochondrial-confirmation gate. |
| P084 | Phage/geNomad viral contig | Not automatically a human pathogen. |
| P085 | Picobirnavirus/smacovirus uncertain host | Host field unresolved/contextual; no automatic human infection. |
| P086 | Virus isolated from human stool in metadata | Isolation host does not establish biological host. |
| P087 | CMV/EBV DNA | No reactivation/end-organ disease/antiviral inference. |
| P088 | Multiple rotavirus segments | One grouped finding with per-segment evidence. |
| P089 | Vaccine-like rotavirus signal | Vaccine/wild-type claim only with validated distinguishing sequence. |
| P090 | Non-F adenovirus | Not relabeled F40/F41 enteric adenovirus. |
| P091 | Normalized fungal/parasite counts | No cells/g, worms/g or viability without independent calibration. |
| P092 | Low-depth completed research analysis | Sensitivity/LOD remain unvalidated; no clinical-negative claim. |
| P093 | Clinical comparison lacked target testing | Target truth remains unknown, not negative. |
| P094 | Repeated participant/multiple preparations | Grouped validation split; no leakage from replicates. |
| P095 | Database/reference strain used for training | Held-out inclusivity tests use independent strains where available. |
| P096 | Many findings and gaps in report | All detail/coverage rows render; no top-N truncation. |
| P097 | No supported findings | No green clearance/donor pass or infection-free score. |
| P098 | Database or threshold update | New versioned results; historical reports reproducible and unchanged. |
| P099 | Existing disease/age/intervention reports | Regression preserves implemented features and separate interpretation. |
| P100 | Entire v05 handoff | Build can proceed from this document plus existing repo/public sources; no supplemental handoff file required. |

## 17. Completion gate and explicit remaining scientific limits

The implementation is complete when the new branch runs end-to-end, every seed has an installed reference or an explicit gap, supported targets pass the applicable computational checks, all result states and report pages work, available public benchmarks are reproducibly evaluated, and the acceptance suite passes. A clinical validation claim is complete only for the targets and assay conditions actually validated.

Known gaps must remain explicit rather than being “filled” by invented data:

- Not every rare human parasite has a complete, uncontaminated, discriminatory reference genome. Marker-only and no-reference rows are legitimate, visible outcomes.
- Newly linked ParaRef artifacts and current portal downloads require acquisition/checksum verification during the build; their full contents were not downloaded in this research.
- Public clinical truth and raw-read availability vary by study/target. Follow the accession/metadata rules; do not invent labels or accession IDs.
- CARD/VFDB and other data have specific reuse terms. Optional access restrictions cannot silently become broad negative findings.
- Clinical LOD, sensitivity, specificity and quantitative load for this application's assay are not established merely by selecting published databases or implementing these provisional rules.
- DNA-only legacy samples cannot gain RNA-virus coverage from software alone. Tissue-restricted disease and low/intermittent shedding also remain outside reliable exclusion by one stool sequence sample.

**Required product result:** a broad, inspectable pathogen-sequence screen that can reveal supported organisms and determinants missed by the existing community profiler, shows the evidence behind every finding, and makes every unassessed target and assay limitation visible. No organism is hidden merely because it is absent from a small handpicked opportunist list; no weak sequence association is upgraded into an infection, treatment instruction or donor clearance.
