# BUILD_SPEC_v08.2 — Gut mycobiome measurement and interpretation

**Project:** OpenBiota. **Research cutoff:** 2026-09-21. **Deliverable:** one self-contained implementation specification for a dedicated mycobiome report page and supporting analysis. No supplemental research file is required.

**Revision: mandatory fungal strain analysis.** This revision supersedes the earlier optional strain extension. Species detection is the entry point; every fungal finding must enter the strain-resolution workflow below. The initial release includes strain/lineage comparison, mixed-population analysis, relevant genetic determinants and donor/recipient comparison. A sample may have inadequate evidence for a particular strain claim; the software may not omit the analysis merely because a species-level profiler does not supply it.

**Revision: required single mycobiome score.** Implement the **Mycobiome Health Score — experimental**, `MHS-E1`, on a 0–100 scale (also serialize 0–1). Section 7.5 supplies the complete initial scoring policy, seed rules, uncertainty and examples. This replaces the earlier prohibition on a composite. It is a transparent summary of the direction of observed fungal evidence, not a validated measure of overall intestinal health or a disease probability. Species, strain and mechanistic evidence feed it at their actual supported resolution; confidence is displayed separately. All implementation instructions remain in this file.

## 1. Decision and scope

**Implement this feature.** Stool shotgun FASTQs can contain useful fungal sequence information. Dedicated fungal references and methods can identify fungi that a predominantly bacterial pipeline misses, quantify the fungal signal in the sequenced specimen, and provide evidence-based interpretations of potentially beneficial, concerning, dietary, or unresolved findings.

The report should answer:

1. How much fungal sequence signal was detected in this sample?
2. Which fungi were supported, at what taxonomic resolution and confidence?
   Include lineage/clade, study/reference strain, sample-specific genotype and supported variants; show the actual resolution achieved and why finer resolution could not be established.
3. Which findings have beneficial-function evidence, opportunistic/pathogenic potential, or uncertain/food-associated interpretation?
4. Is a finding unusual compared with an eligible, method-compatible reference population or the person's previous samples?
5. What sensible follow-up is supported, and what remains unmeasured?
6. What is the single experimental mycobiome score, what moves it up or down, and how much of its interpretation remains uncertain?

This is an additive module. The current v08.1 design already names CGF, EukDetect2 and fungal confirmation; this increment supplies the missing fungal quantification, software corrections, reference/cohort strategy, interpretation registry, and dedicated report interface. All requirements needed for this increment are repeated here. Inspect the actual repository before choosing module paths or adapting CLI names. Neither the repository nor the participant FASTQs were analyzed during this research task; no individual fungal result is asserted.

Preserve existing bacterial scoring, pathogen, strain, function, biofilm and donor-matching outputs. New fungal findings must use the shared identity/provenance system and remain distinguishable from those existing scores.

## 2. What “less than 1%” actually means

Fungi are usually a low-abundance component of stool microbial sequencing, often far below 1%. That does **not** mean everyone should have a specific percentage of “good fungi.” In the HMP mycobiome study, approximately 0.01% of more than 27 billion shotgun reads mapped to fungal genomes. That was a method- and cohort-specific read fraction, not a healthy cutoff or the proportion of living fungal cells. [Nash et al., 2017](https://link.springer.com/article/10.1186/s40168-017-0373-4).

Keep these quantities separate:

| Quantity | Available from existing FASTQs? | Meaning |
|---|---|---|
| Fraction of analyzed fragments confidently assigned to fungi | Yes, after validated competitive classification | Sequence signal under a specified assay/reference; can underestimate unrepresented or poorly recovered fungi. |
| Fraction of analyzed sequence bases assigned to fungi | Yes, with explicit base accounting | A sequence-content estimate; not identical to fragment fraction when lengths differ. |
| Composition within the detected fungal community | Yes, subject to sufficient support | Distribution among detected fungi, normalized within fungi only. |
| Fungal genome-equivalent abundance | A conditional estimate with a qualified model | Depends on genome length, ploidy, copy number, extraction and mapping; not automatically cell abundance. |
| Fungal cells, biomass, viable organisms or copies per gram of stool | Not from ordinary relative FASTQs alone | Requires appropriate calibrated quantitative/viability measurements. |
| Fungal growth phase, gene expression, tissue invasion or active toxin production | Not established by DNA alone | Requires other observations/assays. |

Genome size, fungal cell walls, variable lysis efficiency, ploidy, repeated sequences and nonviable dietary DNA all affect interpretation. Do not set a green zone at “0–1%,” label >1% infection, label zero healthy, or label low fungal diversity deficient.

Stool also does not fully represent the intestinal mucosa. A 2026 study found *Cladosporium sphaerospermum* depleted in Crohn ileal mucosa while fecal abundance was unchanged; protective mechanisms were investigated experimentally. This is a concrete example of biology that stool abundance alone can miss. [Huang et al., 2026](https://pubmed.ncbi.nlm.nih.gov/41501528/).

## 3. Required report page

Title: **Your gut mycobiome**. Subtitle: **The fungi detected in your stool—and what the evidence says about them.**

Place one complete overview page among the early visual summaries. Detailed organism/evidence pages can follow later. A sample with no supported fungal calls still gets this page, explaining analytical coverage and limits.

### 3.1 Page composition

| Panel | Required content |
|---|---|
| Mycobiome Health Score — experimental | Prominent 0–100 horizontal gauge, higher = more favorable observed evidence under the frozen policy. Show the §7.5 sensitivity range, confidence, strongest contributors and completeness next to the number. This is the first panel, also included in the report's main summary. |
| Fungal signal | Percentage of QC nonhost fragments assigned to fungi, matching fragments per million, informative fragment count, and sampling uncertainty. State denominator in the label. |
| Fungi detected | Supported named species, unresolved fungal complexes/higher-rank findings, and provisional findings as separate counts. |
| Strains and functional variants | For every fungal finding: analysis status, achieved resolution, nearest reference/equivalence group, discriminatory evidence, mixed-strain status and relevant measured variants. Show phenotype evidence separately from genetic identity. |
| Composition | Horizontal bars for quantifiable fungi using one qualified within-fungi estimator. Show unresolved read support separately unless its mass can be estimated in the same units. No unlabeled pie implying fungal cell counts. |
| Potentially beneficial evidence | Specific detected organism/strain, studied property, human/animal/laboratory evidence tier, and whether the relevant strain/function was actually resolved. |
| Potential concerns | Specific supported findings, host/clinical context, reference-relative expansion where qualified, measured determinants where supported, and the strength of evidence. |
| Context and confidence | Diet/probiotic/antibiotic/antifungal exposures, data adequacy, relevant uncertainty, compatible reference availability and assay version. |
| What to do next | Result-specific follow-up from §11, with links to evidence and deeper cards. |

Use blue for measured composition, neutral gray for unknown/unassessed data, and amber for contextual findings worth reviewing. A green evidence label means a beneficial property was demonstrated under the stated conditions; it does not certify the organism or the sample as universally beneficial. Red clinical alerts belong only to the existing validated alert policy, not to generic detection of Candida or mold DNA.

Show total nonhost fungal sequence percentage and within-fungal percentage in different columns. Example layout text, **illustrative only and never substituted into a real report**:

> Fungal sequence signal: 0.020% of analyzed nonhost fragments (200 fragments per million). Of the measured fungal community, 70% was assigned to taxon A. These percentages use different denominators.

The within-fungal value can be large while total fungal signal is small. Increasing fungal richness is not inherently an improvement. Do not imply that all fungi should be eliminated or that every healthy gut requires a fixed list of fungal species.

### 3.2 One health summary score, with measurement ranks kept distinct

Provide **reference percentiles for fungal signal and selected taxa only when the reference cohort is eligible and identically processed**. Label them as comparison ranks, not probabilities of disease or scores of goodness. Before a compatible reference is qualified, still show measured signal, taxa, confidence and evidence cards; do not suppress the page.

Implement the single score in §7.5 now; an unavailable population reference must not block its evidence-based computation. Keep `beneficial_evidence_findings` and `potential_concern_findings` as the auditable inputs and explanatory detail. A taxon can appear in both. Do not count papers, organisms or detected genes as independent health points. The score's 0–100 normalization is a mathematical scale, not a population percentile. An independently validated outcome model is a future, separately versioned improvement, not a hidden prerequisite for this experimental index.

The gauge endpoints read **More concerning evidence** and **More favorable evidence**. The midpoint reads **Neutral / mixed evidence**, never “normal.” Show the integer score without a percent sign, its sensitivity range, and the word **Experimental** in the same panel. Use amber toward the concerning end, neutral gray around 50 and teal toward the favorable end; avoid unvalidated red/green diagnostic zones or labels such as “excellent health.” Display “Mixed signals” whenever both benefit and concern contributions are positive. A numerical midpoint can reflect conflicting findings, not reassurance. Where the range crosses 50, label direction uncertain even when the point estimate lies on one side. Confirmed alerts remain prominent at every score.

No fungal evidence or a failed analysis displays **Not enough information to score**, with a reason, rather than 0, 50 or 100. This is a data limitation, not an option to skip computation: all eligible observations must be evaluated and the score returned whenever at least one directional contribution is supported. Browser-rendered source links use `target="_blank" rel="noopener noreferrer"`.

## 4. Reference databases and reusable methods

Freeze assembly accession versions, sequence hashes, taxonomy snapshots, code commits/container digests, licenses and source URLs. Published genome/taxon counts are provenance checks, not measured sensitivity.

| Resource | Evidence, access and role | Integration requirement |
|---|---|---|
| **Cultivated Gut Fungi (CGF), Cell 2024** | 760 cultured isolate genomes, 206 species, 48 families; 69 previously unidentified species. [Paper](https://pubmed.ncbi.nlm.nih.gov/38776919/), [author repository](https://github.com/yexianingyue/Cultivated-Gut-Fungi), [PRJNA833221](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833221). | Required gut-isolate genome complement. These are isolates, not 760 stool metagenomes. At audit, the BioProject exposed 760 BioSamples but 709 assembly records; reconcile the paper's accession table and alternate author sources rather than claiming all 760 assemblies were retrieved. |
| **FungiGutDB v1.0 / FungiGut** | December 2025 preprint and public Nextflow workflow; curated list of 304 taxa associated with culture-based human-gut literature. [Preprint](https://doi.org/10.64898/2025.12.03.691829), [code](https://github.com/diegocoleto7/FungiGut), [analysis code](https://github.com/diegocoleto7/FungiGut-Paper), [frozen database](https://zenodo.org/records/17581472). | Required candidate reference/baseline comparison. Code GPL-3.0, deposited database CC BY 4.0. The approximately 12.47-GB fungal archive is distinct from optional host/UHGG archives. Do not interpret its fungal-normalized output as total sample fungal percentage. Apply §5 software corrections before production use. |
| **EukDetect2** | June 2026 preprint; broad eukaryotic single-copy-marker detection. [Code](https://github.com/allind/EukDetect), [preprint](https://doi.org/10.64898/2026.06.24.734308), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13320746/), [database deposit](https://zenodo.org/records/19056625). | Required orthogonal marker lane. MIT code/CC BY 4.0 database. March 2026 deposit, approximately 7.1 GB, incompatible with old EukDetect databases. Filter outputs by fungal taxonomy; other eukaryotes are not fungi. Pin a working code commit and its documented DB, not an ambiguous unpinned package installation. |
| **NCBI fungal RefSeq/GenBank** | Accession-versioned reference and clinical/food/environmental neighbors. [Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/), [taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). | Required completeness and near-neighbor supplement. Exclude known contamination segments, audit genome quality, and keep nonfungal decoys. Public availability is not permission to relicense all upstream annotations. |
| **FunOMIC** | The original 2022 paper describes approximately 1.69 million taxonomic markers from 4,816 genomes representing 1,916 species, plus approximately 3.41 million functional protein sequences. The 2025 benchmark describes a different, 3,060-species downloaded inventory; retain release-specific counts. [Paper](https://doi.org/10.1016/j.csbj.2022.07.010), [code](https://github.com/ManichanhLab/FunOMIC), [database portal](https://manichanh.vhir.org/funomic/). | Benchmark and optional functional expansion. GPL-3.0 code; verify separate database terms. Validate actual artifact counts/install. The 2025 comparison found required repairs; do not make it the sole critical runtime dependency. |
| **Metax, Cell 2026** | Recent cross-domain profiler using genome-coverage modeling. [Paper](https://doi.org/10.1016/j.cell.2026.08.024), [code](https://github.com/hzi-bifo/Metax), [manual](https://github.com/hzi-bifo/Metax/blob/main/MANUAL.md), [archived v0.9.23-beta code](https://zenodo.org/records/20140652), [benchmark datasets](https://research.bifo.helmholtz-hzi.de/downloads/metax/benchmark_datasets/). | Include in the benchmark comparison; promote to production if its fungal precision/recovery improves the validated workflow. AGPL-3.0-or-later code. The authors' cross-domain benchmark is not an independent comparison with EukDetect2. Fractional-index counts require special handling in §5.3. |
| **MiCoP** | Whole-genome fungal/viral community profiling. [Primary paper](https://doi.org/10.1186/s12864-019-5699-9). | Useful algorithmic baseline with identical admitted references; its old default reference is not a current comprehensive fungal database. |
| **UNITE** | Fungal ITS reference and Species Hypothesis system. [Resource](https://unite.ut.ee/), [downloads](https://unite.ut.ee/repository.php). | ITS annotation/orthogonal testing and taxonomy cross-checks. ITS copy number/primer bias prevents direct conversion of ITS percentages into whole-microbiome fractions. Pin SH release and clustering threshold; do not treat every SH as a validated named species. |
| **MycoCosm / JGI** | Broad fungal genome and annotation source. [Portal](https://mycocosm.jgi.doe.gov/), [1000 Fungal Genomes](https://mycocosm.jgi.doe.gov/mycocosm/home/1000-fungal-genomes). | Supplement identified gaps and functional annotations with accession/permission provenance. Environmental reference breadth is not proof of human gut residence. |
| **FungiDB / VEuPathDB** | Specialized fungal annotations and comparative resources. [FungiDB](https://fungidb.org/fungidb/), [VEuPathDB](https://veupathdb.org/veupathdb/). | Optional source after verifying the exact accessible release and terms; use public accession-based sequences as fallback. |
| **UFCG** | Universal fungal core genes for reference/phylogeny auditing. [Paper](https://doi.org/10.1093/nar/gkac894). | Use fungal-specific quality/identity assessment. Do not apply bacterial MAG rules or a universal bacterial 95%-ANI species boundary to fungi. |

The CGF repository's later human-associated NCBI expansion is not a clearly tagged replacement release with a verified new total; do not describe it as a new fixed species count. Assemble a reproducible accession manifest from accessible sources. Deduplicate exact assemblies, synonyms, source-overlapping strains and species complexes before counting.

The ordinary GTDB/GlobDB bacterial/archaeal profile does not provide a fungal percentage. MetaPhlAn fungal calls remain useful candidate evidence, but its bacteria-focused inventory cannot be assumed sufficient for this page.

Audited source snapshots: FungiGut `dcadf124c7910ebda2cc0f7f9544f6709e44e3df`; EukDetect2 `8d69014727b2c5956de30b0811644b6c3a7bd4f3`; CGF `353df4bf0fa31ba1146cd7e729cae8a90ceb5bd5`; FunOMIC `d8ef2438a4e44ed6839e434778d1ca8e32ea3569`. Start reproducibility checks from these commits and the cited immutable deposits; a newer selected release requires a recorded diff, compatible database and fresh qualification. Record reference-database licenses independently from software licenses. Do not copy unlicensed repository code merely because it is public.

## 5. Corrections required before wrapping public code

Public research code is a starting point, not a validation certificate. The 2025 comparison used 18 fungal mock communities, found substantial inter-tool disagreement and implementation problems, and favored different tools for precision versus recovery. Its lowest tested fungal fraction in bacterial-background experiments was about 1%; it does not establish a detection limit for 0.001% fungal signal. [Avershina et al., 2025](https://link.springer.com/article/10.1186/s40168-025-02048-3), [reproducible workflow](https://github.com/Rounge-lab/mock_mycobiome).

### 5.1 FungiGut code audit

The inspected [abundance script](https://github.com/diegocoleto7/FungiGut/blob/dcadf124c7910ebda2cc0f7f9544f6709e44e3df/bin/compute-abundances.py) has the following semantics/issues. Freeze the inspected revision and recheck these statements against the selected release because upstream may change.

| Finding | Required implementation response |
|---|---|
| `PERCENTAGE` normalizes genome-length-adjusted fungal counts to 100%. | Store as a within-fungal estimate. Compute sample fungal signal separately using the original eligible nonhost denominator. |
| Multimappers are allocated with `random.random()` without a seed in that script. | Replace with a deterministic, tested ambiguity/quantification policy; seeded random allocation is acceptable only for reproduction, never as proof of species identity. |
| `min_map` is a matched-alignment-length threshold, not mapping quality. | Rename the wrapper parameter to its true meaning; separately handle MAPQ. |
| `pct_id` uses CIGAR matched-length/total-length; `M` can include mismatches. | Compute actual identity from validated alignment/edit information; do not label this ratio 99% nucleotide identity. |
| Edit distance is read from a fixed optional-field position. | Parse `NM`/`MD` and alignment tags by name with a standards-compliant parser. |
| Read handling assumes adjacent read names. | Explicitly collate paired-fragment alignments; account for mates, secondary/supplementary alignments, duplicate records and input order. |
| Genome-length normalization uses reference information that can depend on the number of included assemblies. | Ensure adding redundant strains does not artificially change species abundance. Use validated effective reference length/genome coverage instead of summing all conspecific assembly lengths as one organism. |

Run FungiGut unmodified only as a documented reproduction baseline. Production can use the reference data and qualified alignment steps with an OpenBiota evidence parser, rather than preserving these failure modes. Keep required GPL notices/obligations if adapting code.

### 5.2 EukDetect2 quantity semantics

Keep its native fields. RPKS is marker-length-normalized support, not sequencing-depth-normalized fungal load. RelEuk/EukFrac-style outputs refer to their documented eukaryotic universe. The separate `eukdetect-normalize` workflow produces RPKSB using sequencing-library size; that is still a sequence-based proxy, not independently measured cells/g. Validate names, units and formula against the pinned source.

Never relabel a 100%-normalized eukaryotic composition as “100% fungus.” Do not use an old database with new software or call protist/helminth marker output fungal. A close-relative cluster that cannot separate, for example, toxin-associated and food-associated Aspergillus relatives must remain unresolved; a species label must not manufacture mycotoxin evidence.

### 5.3 Metax fractional-index semantics

Its fractional index searches selected genome segments, scales reported counts for the reduced reference size, and disables read-level classification. Those scaled counts are model estimates, not the observed fungal-fragment numerator in §7.1. Use a qualified full-index read assignment or the independent competitive-confirmation ledger for that numerator. Preserve observed/expected breadth and presence-likelihood fields as algorithm outputs; their values are not clinical disease probabilities. [Official documentation](https://github.com/hzi-bifo/Metax).

## 6. Analysis workflow

### 6.1 Inputs and eligibility

Accept paired or single-end **shotgun** FASTQs and existing compatible caches. A bacteria-only profile TSV or report JSON is insufficient to recover fungal reads that were never analyzed. 16S sequencing cannot support this module; ITS sequencing is a separate assay adapter with different quantities.

Required metadata: sample ID, stool body site, dates, raw/cleaned/host-filtered input status, read count/length, extraction/library method where known, database/tool hashes, and age/reference eligibility. Collect optional diet, fermented-food/yeast-product/probiotic use, antibiotics, antifungals, bowel prep and relevant clinical context; missing metadata remain unknown.

Maintain read-stage counts **before any bacterial depletion**. A valid original denominator is necessary but not sufficient: filtering may also have removed fungal or ambiguous reads. For comparable total-sample fungal signal, require the original eligible pre-depletion data or a validated, demonstrably lossless fungal candidate subset plus the original denominator. Otherwise permit detection and clearly labeled observed lower-bound support where calculable, but set total-sample fungal signal to `not_quantifiable_from_available_input`; never normalize the remaining fungi-enriched reads as though they were the original sample.

### 6.2 Reference preparation

1. Combine admissible FungiGutDB, CGF and selected NCBI fungal references through accession/version and sequence-hash deduplication. Preserve alternative strains for confirmation.
2. Audit suspicious bacterial/human/adapter sequence embedded in fungal assemblies. Retain mask coordinates and origin. High reference contamination can create apparent fungi from abundant bacterial reads.
3. Add human, bacteria/archaea, plants/food, other eukaryotes and relevant viral/low-complexity decoys for competitive adjudication. Decoy references are not reportable fungal species.
4. Keep taxonomy aliases, complexes and nomenclatural updates. *Candida glabrata*/*Nakaseomyces glabratus* is one concept, not two detected fungi. Do not substitute genus-level evidence for a named species.
5. Track the actual number of admitted fungal species/complexes, reference assemblies, marker-covered taxa and validation-qualified targets separately.

### 6.3 Detection and confirmation

Run EukDetect2 and a whole-genome fungal candidate lane. Reuse existing MetaPhlAn/pathogen evidence as additional candidates. Do not require both primary lanes to agree for every accepted call: well-supported whole-genome calls outside the marker catalog can qualify after confirmation.

For each candidate, competitively align the complete eligible data or a demonstrably lossless candidate subset to the fungus, its near neighbors and decoys. Use fungal-specific thresholds qualified in §12. Maintain the complete read ledger even when optional bacterial filtering accelerates search.

Require taxon-discriminating evidence across separated nuclear genomic regions; audit rDNA, mitochondrial sequences, repeats, conserved genes and mobile/shared sequences separately. Organellar/rDNA-only findings may be visible as provisional evidence but must not silently receive whole-genome quantification. Resolve multi-mapping at the lowest supported common taxonomic rank.

Starting calibration grid for new genome-based calls: 10/20/50 independent informative fragments, 2/3/5 separated informative regions, and several identity/competition thresholds. These are values to test, not clinical cutoffs. Select rules on calibration/training positives and near-neighbor/contamination negatives, lock them, then evaluate independent held-out strains, communities and studies. Never tune on the final evaluation set. Existing validated marker rules can have their own support model.

Output states: `supported`, `provisional`, `ambiguous_complex`, `not_supported`, `not_detected`, `not_assessed`. Separately record quantification state: `quantifiable`, `detected_below_quantification_limit`, `no_valid_denominator`, `not_quantifiable_from_available_input`, or `not_assessed`. Low quantitative precision must not erase a valid qualitative detection.

### 6.4 Required fungal strain-resolution subsystem

**Implement and run this subsystem by default. It is part of release completeness.** Each supported species/complex receives an assessment; provisional fungal evidence receives an assessment status without being promoted to a named strain. Missing software/reference artifacts are execution or reference-coverage failures, not biological negatives. Do not use a generic sentence such as “this is a species-level pipeline” to explain away an unimplemented component.

There is direct published precedent: West et al. reconstructed Candida genomes from stool metagenomes, analyzed strain relationships/heterozygosity and investigated copy-number variation. That establishes feasibility under adequate sample coverage, not guaranteed recovery from every low-abundance adult stool sample. [Primary study](https://link.springer.com/article/10.1186/s40168-021-01085-y), [analysis repository](https://github.com/patrickwest/c_parapsilosis_analysis).

#### 6.4.1 Separate resolution from phenotype

Return the strongest supported combination of these outputs, retaining underlying evidence:

| Output | Meaning and required distinction |
|---|---|
| Species/complex | Established by §6.3; never relabel a species as a strain. |
| Lineage/clade | Supported membership of a versioned reference-defined population; not an individual isolate identity. |
| Reference strain or equivalence group | A particular reference genotype is supported against appropriate near neighbors, or several deposited strains remain indistinguishable in the observed regions. A closest hit alone is insufficient. |
| Sample genotype fingerprint | A reproducible masked allele/variant profile even when no deposited strain matches. Preserve it for comparisons; do not invent a named clinical strain. |
| Mixed-population evidence | Supported distinct genotypes, an identifiable mixture estimate, or mixture ambiguity. Diploid heterozygosity is not automatically two strains. |
| Genetic determinant | A supported allele, haplotype, gene-content difference or structural event with exact evidence and carrier attribution. |
| Phenotype evidence | A separate claim about what a study measured, a validated genotype–phenotype association, or an explicitly experimental inference. Genetic similarity does not prove the phenotype occurs in this gut. |

A named reference-strain call is not proof that a sample originated from a particular commercial product. Where commercial isolates are genomically indistinguishable, report their reference-equivalence group and allow the user's recorded exposure to provide separate context.

#### 6.4.2 Reference and marker construction

Use the required panels in §6.5. Preserve all distinct high-quality intraspecies genomes while collapsing exact sequence duplicates into explicit alias groups. Build a species pangenome/variant representation with a stable coordinate backbone, one-to-one orthologs, accessible nuclear masks, repeat/paralog masks and alternative reference sequences. Keep mitochondrial/ITS/MLST information as additional evidence with its lower or different resolution.

Separate the all-isolate candidate-search panel from the coordinate-mapping reference. Naively concatenating hundreds of near-identical conspecific genomes can split read support and collapse MAPQ. Use a selected species/population backbone plus appropriate cross-species decoys and explicit strain discriminators, or a separately qualified graph approach. Preserve all candidate hypotheses and check sensitivity to backbone choice. Record assembly representation separately from biological ploidy: a haploid/pseudohaploid reference assembly does not establish that the sampled fungus is biologically haploid.

For each study strain, retrieve assembly versions or reconstruct from its public isolate FASTQs, link BioSample/run/strain names to the paper's phenotype table and record mapping confidence. A phenotype label must never be inferred from an accession order or similar spelling. Maintain accession, specimen site, host, study conditions, control strains, measured endpoint, assay result and source/correction links.

Generate candidate discriminatory SNP/indel, haplotype, gene-content and structural-junction markers from the admitted panel. Screen them against broad conspecific genomes, near-neighbor species and nonfungal decoys. A feature seen once in one assembly is a candidate, not a validated strain barcode. Recheck raw isolate reads, repeat structure and assembly quality; use independent held-out strains and unseen relatives to measure false identity assignment. Keep “none of the panel” as a possible outcome.

Each target definition must contain reference accession/version and checksum; coordinate convention; REF/ALT or junction sequence; strain/clade membership; ambiguity/equivalence group; unique flanks; positive/negative genomes; evidence provenance; benchmark-supported read requirements; and claim scope (`lineage`, `reference_genotype`, `determinant`, `phenotype`). Generate these machine-readable assets during implementation from the specified sources; no guessed nucleotide coordinates are authorized.

#### 6.4.3 Required computational steps

1. **Recover reads and align.** Reuse compatible fungal alignments or realign retained QC nonhost FASTQs. A bacterial marker profile/Bowtie summary is insufficient. Use a pinned whole-genome aligner such as BWA-MEM2 or Bowtie2 against the selected fungal references and competitor/decoy set; retain BAM/CRAM, reference hashes, mate/read-group information and rejected/ambiguous evidence. Avoid deleting potentially informative fungal reads through an irreversible bacteria-first filter.
2. **Measure callability.** Store nuclear breadth/depth, marker callability, MAPQ/base quality, strand/read-position bias, duplicate/overlapping-mate policy and mapping competition. Mask repeats/paralogs and represent uncovered sites as unknown rather than reference alleles. A genome-wide average must not conceal zero coverage of the actual distinguishing loci.
3. **Count alleles with a fungal model.** Use a pinned SAMtools/BCFtools or FreeBayes-based adapter plus an explicit allele-evidence table. FreeBayes supports ploidy settings and pooled-frequency calls; neither setting automatically deconvolves a mixed stool population. Retain allele counts, genotype likelihoods, quality, phase sets and callable reference positions. [BCFtools](https://github.com/samtools/bcftools), [FreeBayes](https://github.com/freebayes/freebayes).
4. **Handle biological complexity.** Configure species/region ploidy, heterozygosity and copy-number uncertainty; investigate loss of heterozygosity, aneuploidy, hybrids and conspecific mixtures. Do not run a haploid default across diploid Candida/Saccharomyces and interpret heterozygotes as contaminating strains. Snippy's haploid workflow is not a general solution for this task. [Snippy scope](https://github.com/tseemann/snippy).
5. **Compare hypotheses.** Score observed alleles/junctions against each admissible reference lineage/genotype and explicit mixture/unknown hypotheses. Use ploidy-aware genotype likelihoods or a calibrated overdispersed allele-count likelihood, accounting for sequencing/mapping errors. Do not count linked nearby SNPs as independent proof; assess uncertainty by genomic blocks and validate the decision margin on held-out strains. Return an equivalence set if the observed loci cannot distinguish the leaders.
6. **Resolve mixtures where identifiable.** Fit candidate genotype contributions using dispersed discriminating loci with nonnegative weights summing to one; include an unknown component and test residual mismatch. Compare a single heterozygous genotype with alternative multi-genotype models. Validate minor-genotype detection by abundance and coverage. Report only partial/unphased genotype evidence when several explanations remain compatible; do not merge their alleles into one fictitious strain.
7. **Assemble when warranted.** For sufficiently covered populations that disagree with references or contain informative unmapped sequence, run a validated assembly/targeted-assembly branch and fungal contig assessment. Retain metagenome assembly identities and check contamination, haplotype collapse and completeness. Assembly is an escalation path, not a prerequisite that prevents low-depth marker/variant assessment.
8. **Call determinants and compare samples.** Apply §§6.6–6.7, store per-operation status and propagate supported strain-level evidence into the report, interventions, pathogen and FMT modules. A missing clinical classifier must not suppress a valid genotype, lineage or variant finding.

The implementation must retain both species discovery and intraspecies variation; EukDetect2/FungiGut/Metax species abundance outputs alone do not satisfy these steps. Existing bacterial strain tools may contribute reusable components only after fungal-specific validation; they do not establish appropriate ploidy or clinical meaning by their presence in the pipeline.

Mixture weights are modeled DNA/genome-copy contributions under the chosen estimator, not automatically proportions of fungal cells. Core phylogeny masks must not permanently erase clinically relevant repetitive loci: recover such loci through dedicated, validated repeat/junction/CNV assays. For CNV, use GC/mappability-adjusted local depth against suitable single-copy regions, report relative copy ratios when ploidy is unresolved, and test whether strain mixtures could explain the signal. A mixture shift is not automatically an acquired LOH/aneuploidy event.

An [inStrain](https://github.com/MrOlm/inStrain) adapter may supply microdiversity/pairwise comparison statistics after fungal qualification; shared-allele population ANI must not be relabeled exact strain identity, particularly for heterozygous or mixed populations. [WhatsHap](https://whatshap.readthedocs.io/) may provide local read-backed phasing under a justified genotype/ploidy model; it is not automatically a solver for unknown mixtures of diploid fungi. Keep phase-block boundaries. Do not coassemble different people or time points for individual strain identity unless every conclusion is independently supported by the original per-sample reads.

#### 6.4.4 Coverage-aware execution, not blanket omission

Assess strain-discriminating loci even when full fungal genome reconstruction is impossible. Report which loci were searched, how many were callable, competing identities and the reason for unresolved results. Do not hard-code a universal minimum fungal abundance that skips all strain analysis.

Use `expected_reference_depth ≈ eligible_sequenced_bases × taxon_sequence_fraction / reference_genome_length` only as a planning estimate. For example, 14 million 2×150-bp pairs provide 4.2 billion bases; if one 12-Mb fungus contributes 0.01% of sequence bases, expected average depth is approximately 0.035×. This can support some detection while leaving most strain loci unobserved. Additional databases cannot recover reads that were never sequenced. Measure actual coverage before deciding whether further sequencing, targeted capture, or isolate whole-genome sequencing would address the limitation.

When estimating additional depth, state the assumed unchanged fungal fraction/extraction and the target callability model; identify capture/culture enrichment as a different assay that cannot replace original relative-abundance measurement. Store unresolved findings with specific reasons: `insufficient_discriminatory_coverage`, `reference_equivalence`, `novel_or_unrepresented_genotype`, `mixed_population_unresolved`, `reference_panel_missing`, `failed_qc` or `execution_error`. Execution/reference failures must be visible in the run's completeness summary.

### 6.5 Required strain panels and source registry

The panels below are initial release assets. Apply the general strain workflow to every other supported fungus; the initial curated panels do not restrict which species can be analyzed. Extend available intraspecies references through CGF/NCBI and emit specific reference-gap records for organisms without adequate genomes.

| Required panel | Verified starting data | Implementation and claim scope |
|---|---|---|
| **Saccharomyces / boulardii lineages** | [2026 genomic comparison](https://doi.org/10.1016/j.isci.2026.115706), [full text and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC13141078/), [PRJNA1312645](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1312645): eight public Nanopore isolate WGS runs. | Build whole-genome and candidate junction/allele panels; validate against broader conspecific diversity, not only the study's three non-boulardii controls. Paper assemblies collapse haplotypes; verify candidate variants against raw reads and independent high-quality truth. |
| **Broad Saccharomyces near neighbors** | [Peter et al., 2018: 1,011 isolate genomes](https://doi.org/10.1038/s41586-018-0030-5), [ERP014555](https://www.ebi.ac.uk/ena/browser/view/ERP014555), [author genome/variant downloads](http://1002genomes.u-strasbg.fr/files/). | Add diverse food, environmental and clinical backgrounds for population placement and false product-strain assignment challenges. Preserve project metadata and biological ploidy. |
| **C. albicans damaging/low-damaging research isolates** | [Li et al., 2022](https://www.nature.com/articles/s41586-022-04502-w), [accessible full text/supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC9166917/), [correction](https://doi.org/10.1038/s41586-022-05102-4), [PRJNA702809](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA702809): 34 public paired Illumina isolate WGS runs; comparative [PRJNA432884](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA432884). | Reconstruct/genotype the panel and join measured phenotypes to exact isolates. Separate engineered ECE1/EFG1 mutants from naturally sampled isolates; 34 runs are not 34 independent natural-phenotype training cases. The four candidalysin sequence variants did not distinguish high- from low-damaging isolates, so they cannot serve as that classifier. |
| **Experimentally protective Clavispora P4013B** | [2025 primary study](https://www.nature.com/articles/s41467-025-64914-w); exact assembly [GCA_023627835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023627835.1/), BioSample SAMN27964042, WGS JAMAIB01, within PRJNA833221. | Compare against other C. lusitaniae strains and validate unique markers. Link a supported P4013B-compatible genotype to its experimental evidence, not proven human benefit. A generic pyruvate-decarboxylase/PDC gene is not an identifier for this strain. |
| **Stool Candida genome recovery** | [West 2021](https://doi.org/10.1186/s40168-021-01085-y), [code](https://github.com/patrickwest/c_parapsilosis_analysis), [PRJNA717139](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA717139), [PRJNA471744](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA471744). | Reproduce a published stool-based positive benchmark for fungal genotype/coverage and relevant structural evidence. Select exact fungal samples/modalities from the study, then challenge dilution and independent near neighbors; infant hospital biology is not an adult normal reference. |
| **Paired gut/blood strain comparison** | [Zhai 2020](https://doi.org/10.1038/s41591-019-0709-7), [PRJNA579121](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA579121): current manifest has 61 isolate WGS and 128 amplicon runs. | Use the 54 fecal/seven blood cultured-isolate genomes for comparison truth; preserve specimen provenance. Archive aliases now include Lodderomyces parapsilosis/metapsilosis/orthopsilosis. These isolates validate genomic relationships; they do not provide a general stool invasion-risk probability. |
| **Other CGF/clinical/food fungal strains** | All distinct qualified genomes from [CGF PRJNA833221](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833221), accession-versioned NCBI genomes and matching isolate WGS. [Candida Genome Database](https://www.candidagenome.org/), [2025 CGD publication](https://pubmed.ncbi.nlm.nih.gov/39776186/), [Candida PubMLST](https://pubmlst.org/organisms/candida-albicans). | Build generic intraspecies fingerprint/variant panels for supported taxa, prioritizing Candida/Nakaseomyces, Malassezia, Debaryomyces, Rhodotorula and relevant Aspergillus complexes. Include all accessible C. auris/Candidozyma auris clades and related species when this pathogen is a candidate; do not assume an old four-clade list is complete. PubMLST type/allele profiles aid placement but are not unique whole-strain identities; mixed/incomplete loci must not generate a synthetic sequence type. |

**Exact Saccharomyces starting runs:** I-745 = SRR35186293/SAMN50873286; DBVPG6763 = SRR35186292/SAMN50873287; MYA796 = SRR35186291/SAMN50873288; I-3799 = SRR35186290/SAMN50873289; I-1079 = SRR35186289/SAMN50873290; Swedish baker = SRR35186288; Danish baker = SRR35186287; CEN.PK113-5D = SRR35186286. Resolve FASTQ URLs/checksums through the official archive. Candidate traits include the reported approximately 68-kb chromosome-XVI inversion; unique breakpoint-spanning support is required because sequence inside an inversion does not establish orientation. Derive deployed coordinates/junctions from the exact aligned assemblies, not an unversioned paper chromosome label.

The 2026 comparison's Table S12 contains candidate isolate-discriminating sites, including sites where I-745 carries the reference allele. Search REF as well as ALT evidence. Do not count tightly clustered sites as independent markers or deploy nanopore assembly differences before error/heterozygosity verification. Retrieve Table S12 from the paper's supplements or the [Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13141078/supplementaryFiles), which contains `mmc2.zip` and the workbooks. Verify supplement permissions before redistribution; compute deployment markers from accessioned sequences. Do not depend on the paper's draft-like Mendeley link.

Verified Li-study starting isolate/run joins: high-damaging IDB311 = SRR13741117 and IDB101 = SRR13741122; low-damaging IDC561 = SRR13741104 and IDB891 = SRR13741110. Preserve the original assay labels and reconcile the correction/source measurements before model development; these four examples do not constitute a validated phenotype classifier.

For C. auris, [PRJNA328792](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA328792) provides an established assembly-discovery starting point. B8441 is accessioned as GCA_002759435.3 / corresponding GCF_002759435.1, BioSample SAMN05379624; preserve the exact selected assembly namespace/version and coordinate hashes. Extend with current genomes and source-verified clades/near neighbors. [2026 six-clade research](https://doi.org/10.1038/s43856-026-01642-2) is relevant, but its PRJNA1447732 yielded no ENA read-run rows during this audit. Treat inaccessible genomes as a specific coverage gap, not proof they do not exist. A clade assignment does not establish the resistance genotype or clinical infection.

**Additional current validation source:** [Pan et al., September 2026](https://doi.org/10.1038/s41591-026-04616-y) compares 30 paired oral/fecal C. albicans isolates using SNP/LOH evidence. The paper mentions shotgun data under PRJNA1449485, but the audited public manifest contains 180 amplicon runs, with no verified isolate-WGS runs. [Published code/archive](https://zenodo.org/records/21072154) supplies statistical analysis rather than a released strain caller. Track this as an explicit availability gap and proceed with the verified panels above. Do not call an amplicon-only archive a usable WGS strain panel.

### 6.6 Required determinant and phenotype evidence layer

Evaluate relevant loci for every supported fungal taxon with an applicable curated panel, including supported/provisional status and reasons when not callable. Keep five fields separate: sequence/allele observed, carrier species, carrier strain, measured/associated phenotype and phenotype prediction status. Do not require complete named-strain identification before reporting an independently supported determinant.

Use locus-specific breadth/depth, species-specific reference coordinates, unique mapping and verified allele definitions. Require linkage/phasing or another validated joint attribution method before joining separated alleles into one strain or multi-locus resistance genotype. Repeats, copy-number changes and low depth need their own uncertainty. Variant absence is reportable only when the relevant allele/region was adequately interrogated; an uncovered locus is not wild type or susceptibility.

Annotate with the species/reference-specific nuclear genetic code and exon structure. Several yeasts use the alternative yeast nuclear code (NCBI table 12); do not translate all fungi with bacterial table 11 or assume every Candida-named organism uses the same code. Keep mitochondrial translation rules separate. [NCBI genetic-code definitions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG12).

Start with ECE1/candidalysin-related sequence variation and relevant studied Candida virulence/biofilm loci; species-specific antifungal determinant panels such as ERG11 and FKS hotspots where primary genotype–phenotype evidence exists; and supported structural/CNV mechanisms. Broad efflux-gene presence is not resistance. The exact allele/coordinate/source and species applicability are mandatory; similarly named fungal genes are not interchangeable.

**FungAMR integration:** the 2025 publication describes 35,792 curated entries across 95 fungal species, with experimental support and contradictory/non-resistance observations. Import exact species, protein reference, mutation, drug, experiment and evidence grade rather than treating every entry as proven resistance. [Paper](https://www.nature.com/articles/s41564-025-02084-7), [database](https://card.mcmaster.ca/fungamrhome), [versioned source repository](https://github.com/Landrylab/FungAMR). Positive evidence grades run from 1 (stronger) to 8 (weaker); negative grades encode evidence that a mutation did not cause resistance. Never sort greater absolute numbers as stronger proof or drop negative evidence. Collapse duplicate study observations separately from distinct alleles.

[ChroQueTas](https://github.com/nmquijada/ChroQueTas) is a useful assembly/proteome annotation adapter, not a raw mixed-FASTQ strain caller; its documented mutation/species subset does not cover all FungAMR mechanisms. Confirm assembly calls against participant reads. Its README and FungAMR source specify CC BY-NC-ND 4.0: do not bundle these as unrestricted open-source/commercial assets or distribute modified derivatives without appropriate permission. Keep a terms-compatible optional external-data adapter, and implement the core detector with independently curated primary-study facts and appropriately licensed sequences. Licensing of this particular resource must not remove the required strain/determinant capability. Record exact versions and permitted use.

The biofilm report receives measured genetic capacity and supported carrier attribution. It does not receive an asserted active biofilm, hyphal transition or toxin production from DNA alone. The intervention engine can retrieve evidence matching a determinant/strain and retain preclinical experiments; it cannot pretend a study phenotype was measured in the participant. Each unresolved strain question must still leave species-level evidence and supported variant findings visible.

### 6.7 Mandatory longitudinal and donor–recipient strain comparison

When multiple samples are supplied, compare the fungal genotypes under a shared reference/mask/version. Evaluate pairwise callable intersections, report their size and discriminatory information, allele concordance/conflicts, mixture ambiguity and reference equivalence. Never treat zero observed differences across a handful of covered bases as proof of an identical strain. Calibrate relationship thresholds by species using repeated isolates/samples and unrelated near neighbors; no universal bacterial ANI/SNP threshold is allowed.

Return `compatible_genotype`, `distinct_genotypes`, `shared_lineage_only` or `insufficient_comparable_evidence`, with the uncertainty and measured loci. For post-FMT samples, identify genotypes consistent with donor origin, recipient persistence, or unresolved/shared origin. Use every available pre-FMT donor/recipient sample; if both donor and recipient carried an indistinguishable genotype, source is unresolved. Baseline non-detection at poor depth is not proof of absence, and a post-transfer stool DNA hit alone is not durable engraftment. Report persistence only across appropriate longitudinal observations.

Donor matching must show known strain/determinant findings and unresolved questions separately. Do not infer fungal compatibility, absence of transmissible risk or engraftment probability from species overlap alone. Do not discard useful strain comparisons merely because a validated disease-risk percentage is unavailable.

## 7. Quantification: formulas and denominator rules

### 7.1 Primary observed signal

Define `N` as all eligible QC nonhost fragments entering comprehensive classification, **before bacterial filtering**. Count each read pair once; retain singletons with a separately documented policy. Define `F` as fragments confidently assigned to fungal ancestry after competitive adjudication, including fungal higher-rank/complex assignments. Exclude cross-kingdom ambiguous assignments from `F` and report them separately.

`fungal_fragment_fraction = F / N`

`fungal_fragments_per_million = 1_000_000 × F / N`

These two numbers must agree exactly before display rounding. The label is **fungal sequence signal (% of analyzed nonhost fragments)**. It is not “percentage of living microbes that are fungi.” Supply counts and a sampling interval, clearly distinguished from extraction/classification uncertainty.

For base-level signal, sum nonoverlapping eligible fungal-assigned read bases divided by eligible nonhost bases with explicit paired-overlap treatment. Do not interchange this with fragment fraction. Any estimated correction for mappability/reference coverage is a separate model output with its own validated uncertainty.

Keep a separate `fungal_fraction_among_assigned_cellular_microbes` if desired, defined as fungal assignments divided by fungi + bacteria + archaea + other cellular-microbe assignments. Exclude plant/host/viral/ambiguous sequence according to a frozen policy. It is conditional on the assigned microbial universe and must never replace the primary denominator without a label.

No integer read can contribute to two mutually exclusive categories. Unassigned sequence remains in the primary `N` denominator; do not assume it is bacterial. Report its fraction so database changes cannot invisibly alter the meaning of the fungal percentage.

### 7.2 Within-fungi composition

Choose one qualified estimator and a coherent species/complex partition. If using coverage-normalized abundance, define `q_i = c_i / sum(c_j)` where `c_i` is a qualified fungal nuclear genome-coverage or marker-depth estimate. If using assigned-fragment composition, label it explicitly as such. Do not combine the two estimators into one bar chart.

Complex-level mass stays at complex level. A secondary-lane organism can be supported while having unavailable primary quantitative abundance. No maximum/mean/sum merging across independently normalized outputs. Adding redundant reference strains or synonyms must not change biological mass or organism count.

An unresolved fungal-fragment fraction uses an explicitly stated read denominator; do not insert it into a genome-equivalent composition bar. Show an unresolved composition segment only when the same estimator can quantify that segment. Higher-rank evidence may overlap named descendants: a Candida-genus assignment plus C. albicans is not automatically two organisms. Keep unresolved assignment counts separate from species richness unless discriminatory evidence establishes an additional distinct taxon.

When only one species is quantifiable, “100% of the quantified fungal community” is permitted with the total fungal signal shown beside it; it does not mean the gut is 100% fungal or that all other fungi are absent. If fungal support is below the qualified quantification limit, show detections and counts but omit precise composition percentages.

### 7.3 Diversity, confidence and trends

Report observed fungal richness and optionally Shannon/evenness only when rarefaction/coverage demonstrates sufficient stability. These are ecological descriptors, not universally directional health measures. Zero fungal detections do not yield a reassuring diversity score of zero.

Longitudinal comparisons require compatible extraction, sequencing, pipeline and database—or an explicit harmonized reanalysis. Give effect size and uncertainty; annotate yeast products, antifungals and other exposures. Changes in bacterial DNA alone can change the fungal percentage without any increase in absolute fungal cells.

At 8–14 million eligible pairs, a hypothetical fungal sequence fraction of 0.001% yields only about 80–140 fungal pairs, before marker and specificity filtering. At 0.0001%, it yields about 8–14. These arithmetic examples illustrate sensitivity limits; they are not validated detection thresholds or claims about the user's samples.

### 7.4 Optional quantitative follow-up

Import fungal qPCR/ddPCR, culture or spike-in calibrated results as separate assays. Store specimen, target, units, calibration, extraction, target-copy assumptions and uncertainty. ITS/rDNA copies cannot be silently converted into cells because copy number varies. A spike-in added after extraction cannot measure fungal lysis losses that already occurred.

### 7.5 Required single Mycobiome Health Score — experimental (`MHS-E1`)

#### 7.5.1 What is being scored

The initial target is **the direction of the strongest interpretable benefit and concern signals detected in this stool sample, conditional on the available evidence**. It is not a measurement of total fungal wellness, cumulative infection risk, gut-wall activity or treatment response. This explicit target allows an immediate, reproducible summary while preserving the substantial limitations of current fungal science.

The rationale is grounded in the primary studies already linked in §§2, 8 and 9: common healthy fungi can also be opportunists; fungal phenotypes can differ within a species; food/probiotic exposure affects stool signals; and experimental protection does not establish protection in a particular human. None of those publications supplies the formula below. **The formula, caps and display rules are OpenBiota engineering policy, not published clinical coefficients.** Version and expose them accordingly. Do not describe the index itself as clinically validated because its component observations have published support.

Separate three questions throughout: is the sequence assignment supported; is its connection to a studied phenotype supported; and does that phenotype predict a human outcome? A successful strain call answers the first question more specifically, not all three.

#### 7.5.2 Frozen contribution rules

Each eligible rule has `rule_id`, `direction`, `biological_claim_id`, `source_ids`, `required_identity`, `required_features`, `required_context`, `trigger_type`, `base_cap`, `activation`, `activation_range`, `extrapolations`, and `finding_ids`. The source-backed *fact* and the policy *cap* are different fields. No LLM may assign a polarity, weight or dose ad hoc during report rendering. Resolve synonyms, multiple databases, dependent genes and repeated publications to the same biological claim before scoring.

| Initial evidence route | Default cap on a contribution | Required handling |
|---|---:|---|
| Resolved genotype matching a specifically studied beneficial or damaging strain | 0.50 | Require the phenotype-linked isolate/genotype rule in §6.5, adequate discriminating evidence, and the study's actual direction. Mark as studied potential; DNA does not confirm activity or benefit/harm in the participant. |
| Resolved species with relevant experimental gut evidence, without matching the tested strain | 0.25 | Explicitly label species-level extrapolation. Never claim the experimental strain was detected. Direction must be tied to a named endpoint, not generic “good” or “bad” taxonomy. |
| Human stool abundance–outcome association | 0.25 | Require an applicable, quantitatively comparable measurement and qualified reference contrast in the published direction. Presence alone does not activate an “increased abundance” rule. |
| In-vitro mechanism demonstrated for the observed identity/feature combination | 0.125 | Require relevant assay and identity; no transfer of a tissue-invasion or expression phenotype from a generic housekeeping/adhesion gene. Keep the laboratory endpoint explicit. |
| Generic opportunist label, gene capacity without directional phenotype evidence, uncharacterized organism, mere ecological unusualness | 0 | Still display the finding, context, alerts and evidence. Zero contribution means unscored, not benign. |

These numbers deliberately limit the influence of weaker or less specific evidence; they are not probabilities, severity estimates or measured effect sizes. The initial sequence-only preset uses no cap greater than 0.50. The scale is still normalized to 0–100; do not stretch the realized scores to fill that range or force every sample to one extreme. Any future increase in caps or addition of an outcome-trained component requires a new `model_id` and documented validation. Host context can determine applicability but must not invent a disease diagnosis.

For a supported presence/genotype rule, `activation = 1`. If unresolved, use a supported lower-resolution rule for the point estimate if one exists; retain the unresolved higher-resolution claim only in the sensitivity scenarios. Do not convert a classifier confidence, read fraction or abundance into an activation multiplier. There is no established universal dose–benefit or dose–harm function for fungal abundance. Features below their validated calling threshold are provisional rather than weakly positive points.

For a supported **high-abundance association** rule only, define `activation = max(0, 2p - 1)`, where `p` is that taxon's qualified same-pipeline midrank reference percentile divided by 100 (§8). This is a declared descriptive ramp above the reference median, not a pathological cutoff. Its maximum contribution remains 0.25. It requires the original association's direction and specimen/population compatibility; an ITS or mucosal rank cannot activate a stool-shotgun rule. If the reference is missing, leave this activation unassessed and still score other eligible rules. The initial preset contains no rule awarding health points for a missing organism, low fungal burden or a low-abundance depletion contrast.

For known contradictory or context-dependent evidence, retain both applicable directions. A missing required host-context field does not silently mean “healthy”: mark that rule unresolved and evaluate explicit applicable/nonapplicable sensitivity scenarios. Experimental model context (for example, mouse colitis) must remain disclosed as extrapolation; it does not require assigning that disease to the participant. A known contrary patient context can make a human-indication-specific rule inapplicable.

#### 7.5.3 Required initial rule manifest

Implement these rows as the complete initial preset; sources and exact accessions are embedded in §§6.5 and 9. `B` means beneficial-potential contribution; `C` means concern contribution. Additional literature stays in evidence cards until a new explicit rule with the required fields is reviewed and versioned. This table is enough to build the first scorer without an absent companion file.

| Rule ID | Point trigger and contribution | Sensitivity / interpretation requirement |
|---|---|---|
| `P4013B_B` | Resolved P4013B reference genotype (`GCA_023627835.1`) → B=0.50. | Link to the [2025 experimental protection study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12615640/). Presence of generic PDC or species identity alone does not activate this row. |
| `CL_LUSITANIAE_B` | Supported C. lusitaniae species without a resolved P4013B match → B=0.25. | Explicit species-to-tested-strain extrapolation; also applies to a resolved different genotype unless genotype-specific contrary evidence invalidates the extrapolation. Same biological family as `P4013B_B`, so no double counting. Upper B=0.50 only when P4013B remains genuinely compatible with observed discriminating evidence, not merely because it belongs to that species. Excluding P4013B removes that alternative, not the underlying species evidence. Human benefit is not established. |
| `CA_HD_C` | Resolved natural HD isolate IDB311 or IDB101 from §6.5 → C=0.50. | Link to [Li 2022 and its correction](https://pmc.ncbi.nlm.nih.gov/articles/PMC9166917/). Exact study-isolate genotype evidence is needed; engineered mutants are separate. Generic ECE1, candidalysin variants, species identity, or a nearest-neighbor result alone do not activate this rule. |
| `CA_LD_CONTEXT` | A resolved studied LD isolate contributes no automatic B points. | Lower damage in one assay is not a demonstrated health benefit. If an observed genotype equivalence set contains HD and LD references, use C range 0–0.50 for that concrete ambiguity; do not choose the preferred member. |
| `CPAR_GUT_C` | Supported C. parapsilosis species → C=0.25. | [Experimental diet-related metabolic phenotype](https://pmc.ncbi.nlm.nih.gov/articles/PMC8546080/); clearly label transfer from tested isolates/models. An unresolved parapsilosis complex does not identify this species. No human obesity probability follows. |
| `DH_REPAIR_C` | Supported D. hansenii species → C=0.25. | [Experimental intestinal wound-repair evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC10114606/). Stool-to-injured-tissue and strain extrapolations remain visible; no claim of impaired repair in this person. |
| `MR_COLITIS_C` | Supported M. restricta species → C=0.25. | [Experimental colitis and host-context evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC6417942/). Keep susceptibility/model/site limitations; unknown CARD9 genotype is not imputed. |
| `CT_EC_SM_INTERACTION_C` | Supported C. tropicalis plus E. coli plus S. marcescens → C=0.125. | [Laboratory interaction study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5030358/). Species co-detection is experimental interaction potential, not a detected gut biofilm. Within one sample, multiple biofilm genes/papers cannot multiply the contribution. |
| `SB_STUDIED_B` | Resolved probiotic genotype with an exact, verified genotype-to-beneficial-study join and applicable indication context → B=0.50. | A product label, boulardii-like lineage or the 2026 isolate panel alone is insufficient to identify the strain from an older trial. Leave the rule unbound until that exact source join is verified. The bound manifest must name the indication and required context: an antibiotic-associated-diarrhea trial requires documented relevant antibiotic exposure, not an unrestricted baseline-health claim. Missing required context leaves only an explicit applicable/nonapplicable sensitivity alternative, not a point contribution. No points for ordinary S. cerevisiae or mere probiotic use; those findings remain visible. |
| `FUNGAL_ABUNDANCE_ASSOCIATION` | A reviewed human stool rule satisfying §7.5.2 → B or C = 0.25 × activation. | An unbound generic adapter, not permission to label all high Candida or all low Saccharomyces unhealthy. Requires its own exact study/feature/reference manifest before activation. |

Do not add all protective consortium members as independent beneficial fungi, assign a mucosal protective phenotype from fecal absence, use mold DNA as measured mycotoxin, or use resistance-marker presence as current intestinal harm. Existing confirmed fungal pathogen/resistance alerts remain visible independently of whether the health index can interpret them. A suspected pathogen must not disappear behind a favorable score.

#### 7.5.4 Formula and exact behavior

For every supported, active directional claim, compute `v_j = base_cap_j × activation_j`. Deduplicate claims; retain all evidence provenance. Define:

```text
B = max(v_j for eligible beneficial claims, default=0)
C = max(v_j for eligible concern claims, default=0)
score_01  = 0.5 × (1 − C) × (1 + B)
score_100 = 100 × score_01
```

**Applicability precedes arithmetic.** If there is no supported nonzero directional contribution, both public score fields are null and the reason is `no_directional_evidence`, `insufficient_fungal_information`, `incomplete_analysis` or `failed_qc`, as applicable. A mathematical empty default of 50 must never leak into the report. A failed or omitted mandatory scoring/strain stage is an implementation failure, not no directional evidence. A complete strain attempt can legitimately remain unresolved and still permit a species-level contribution.

The maximum is a deliberate **dominant-evidence** aggregation: repeated studies, reference duplicates, or many neutral fungi cannot inflate a score or dilute a concern. It does not quantify accumulated biological risk; show the full number/list of independent concerns next to it. Adding another concern below the current maximum leaves the number unchanged but still adds its visible finding and follow-up. Do not present unchanged score as unchanged biology.

The product is conservative: a benefit can raise the score, but cannot erase the dominant concern. For C≥0.50, even B=1 cannot produce a score above 50. Absence of a concern contribution means no eligible concern was established under this policy, not proof that no harmful fungus exists. Confidence must not be multiplied into the score; poor sequencing is not poor health, and missing evidence is not good health.

Illustrative arithmetic, not participant results:

| Supported B | Supported C | Score /100 | Meaning under this policy |
|---:|---:|---:|---|
| 0.50 | 0 | 75 | A specific beneficial-potential finding with no activated concern rule; incomplete knowledge still applies. |
| 0 | 0.50 | 25 | A specific concerning-potential finding. |
| 0.50 | 0.50 | 37.5 | Mixed signals; the benefit does not cancel the concern. |
| 0.25 | 0 | 62.5 | Less specific beneficial evidence. |
| 0 | 0.25 | 37.5 | Less specific concerning evidence. |
| 0 | 0 | null | No supported direction; not a healthy midpoint. |

Report whole-number gauge values while storing full precision; compare unrounded values and use deterministic half-up rounding. Do not print “75% healthy.” Whole-fungal abundance, diversity, number of named strains and population normality do not enter this formula by themselves.

#### 7.5.5 Uncertainty, coverage and confidence

For each relevant ambiguous finding, build an admissible scenario set from actual analytical alternatives: strain equivalence members, supported variant uncertainty, contamination adjudication, conflicting study directions, or missing required context. Retain the unrepresented-genotype hypothesis where appropriate. Do not assign an uncharacterized fungus hypothetical pathogenicity 0–1 solely to force a range.

For the resulting envelope B∈[B_lo,B_hi] and C∈[C_lo,C_hi], including the point scenario:

```text
range_low  = 50 × (1 − C_hi) × (1 + B_lo)
range_high = 50 × (1 − C_lo) × (1 + B_hi)
```

This may be a conservative outer envelope if endpoints involve incompatible scenarios; label it **resolution/evidence sensitivity range**, never “95% confidence interval.” Keep the contributing scenarios inspectable. If all directions are unresolved, provide the range as exploratory detail but leave the main point null. For example, a supported B=0.50 plus a detected C. albicans equivalence group containing HD/LD alternatives yields point 75 from the supported evidence and range 37.5–75; the panel must show **Direction uncertain — strain unresolved**, not a favorable overall verdict. New genuinely unrepresented organisms remain outside this range and lower interpretation coverage; the range is not an exhaustive bound on actual health.

Additionally calculate a **policy sensitivity envelope** by rerunning with all nonzero caps multiplied by 0.5 and 1.5 (clip to 1), allowing B and C scaling to vary independently. Include the baseline in the envelope. This exposes the initial policy's arbitrary numerical precision. Serialize the resolution-only and policy-only ranges separately, but derive the **displayed combined range from the Cartesian product of resolution alternatives and independent policy factors**, not the union of the separate ranges. An equivalent conservative bound under this monotone model is:

```text
display_low  = 50 × (1 − min(1, 1.5 × C_hi)) × (1 + 0.5 × B_lo)
display_high = 50 × (1 − 0.5 × C_lo) × (1 + min(1, 1.5 × B_hi))
```

For B=0.50/C=0 without resolution ambiguity, policy sensitivity alone is 62.5–87.5. For B=0.50 and unresolved C∈[0,0.50], combined sensitivity is 15.625–87.5; a simple union of separate ranges would miss the lower extreme. These are scenario values, not inferred population variability. No bootstrap can turn literature-policy weights into validated clinical effect sizes.

Show these separate fields:

- `analytical_confidence`: supported by the detection/strain benchmark and sample-specific callability; use `supported | limited | insufficient`. It is not a probability of good health. Supported requires the dominant contributions to pass their validated assay gates and no unresolved scoring-relevant analytical alternative.
- `interpretation_confidence`: `experimental_limited` for every initial MHS-E1 result. Do not label clinical reliability high because the DNA coverage is high. Preserve the strongest study design, extrapolations and conflicts in detail.
- `scored_fungal_fragment_fraction`: unique supported fungal fragments attributable to at least one active scored finding, divided by all supported fungal fragments, counting each fragment once. Report null where attribution is unavailable; never mix this with genome-equivalent fractions. This is interpretation coverage of detected sequence, not the fraction of health understood.
- `strain_assessment_completion`: successfully completed mandatory taxon strain attempts divided by required attempts, plus the distinct proportion that actually resolved. Unresolved and unattempted are different states.
- `unscored_findings`, `unresolved_concerns`, `active_alert_ids`, and `context_missing`: always visible/retrievable; none contributes reassuring points.

Label the score `provisional` when a scoring-relevant ambiguity, incomplete assessment, unassessed compatible concern, unavailable interpretation coverage, or out-of-domain population prevents a complete reading of its inputs. The observed-evidence number may still render with the exact limitation; QC failure invalidating the contributing evidence makes it null. Age, immune status and exposure determine context, not invented numerical multipliers. Adult reference percentiles must not be transported to children. For children or unvalidated host strata, retain research observations and visibly mark extrapolation rather than claim pediatric health calibration.

#### 7.5.6 Explanation, action and future calibration

Every point has a trace: source observation → taxon/genotype/determinant → evidence claim → rule/cap/activation → B or C → score. Show the strongest up/down drivers and all concerns. Explain that improvement in this *index* can result from new information rather than improved biology. Serial or donor comparisons require identical model, sources, assay compatibility and adequate matched information; rerun both samples under the same lock when these differ. Loss of coverage must never be labeled improvement, and a pooled donor union must not inherit the highest component score.

Actions attach to the actual finding using §11, not the gauge number. A low score can surface its organism-specific human or preclinical intervention evidence; it does not select an antifungal or infer infection. A high score does not clear an FMT donor. The algorithm must never recommend ingesting a detected opportunistic fungus just to raise B.

Implement an optional separate `reference_percentile` of MHS-E1 only after eligible same-pipeline cohorts exist; label it as a score rank, not the health score itself. For a later outcome-trained model, preregister the particular clinical endpoint and population, reprocess independent stool-WGS cohorts, split by participant/study, preserve held-out external evaluation, and report discrimination, calibration, precision/recall and sensitivity to host/exposure confounding. Stool samples from healthy people alone cannot validate a global health scale. Neither six longitudinal donors nor isolate phenotype experiments establish population clinical performance. Version the replacement model and retain MHS-E1 for reproducibility; do not silently rebrand an IBD classifier as universal mycobiome health.

## 8. Open cohorts and validation data

Reference **genomes** identify sequences; reference **cohorts** establish population comparisons. These are different assets.

| Dataset | Modality and verified role | Access |
|---|---|---|
| HMP mycobiome, Nash 2017 | 317 analyzed stool ITS2 samples from 147 volunteers, with linked shotgun/18S analysis. Useful adult prevalence context and cross-assay comparison. Repeated samples are not independent people. | [Study](https://link.springer.com/article/10.1186/s40168-017-0373-4), [HMP shotgun project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA43017), [ITS2 PRJNA356769](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA356769), [analysis code](https://github.com/cmmr/NashMicrobiome2017). Resolve exact sample links; neither all HMP runs nor all 390 archived amplicon runs are the 317 analyzed stool samples. |
| Longitudinal mycobiome/diet, Xie 2023 | 143 shotgun stool samples from six healthy adults over two months, with diet/method work. Useful technical/temporal/diet sensitivity dataset; far too few independent participants for a general healthy percentile. | [Paper](https://link.springer.com/article/10.1186/s40168-023-01693-w), [PRJNA925700](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA925700), [analysis code](https://github.com/ManichanhLab/LongitudinalMycobiomeWithDiet). Keep enriched and unenriched protocol arms distinct. |
| CHILD, Nature Communications 2026 | 2,256 longitudinal samples from 1,409 infants, paired ITS2/shotgun study with later allergy outcomes. Published fungal analysis uses ITS2; shotgun analysis primarily profiles bacteria. Reanalyzing matched WGS for fungi is a new validation task. Not an adult/adolescent normal reference. | [Paper](https://www.nature.com/articles/s41467-026-74418-w), [WGS PRJNA838575](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA838575), [ITS2 PRJNA1368998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1368998), [metabolomics MTBLS7919](https://www.ebi.ac.uk/metabolights/MTBLS7919). Clinical metadata access is governed separately. |
| iHMP / IBDMDB | Longitudinal IBD and control multi-omics; useful shotgun fungal method validation and disease-context comparison. PRJNA398089 currently contains 2,041 WGS, 762 RNA-seq and 176 amplicon runs; select the DNA stool subset and exact participant metadata. | [Project PRJNA398089](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA398089), [published EukDetect application](https://link.springer.com/article/10.1186/s40168-021-01015-y). The mixed project is not a ready-made healthy shotgun reference cohort. |
| Twin mycobiome/bacterial interactions, 2026 | 212 participants, ITS + 16S. Useful ecological and host-context evidence; cannot directly calibrate shotgun fungal fractions. | [Paper](https://doi.org/10.1016/j.isci.2026.115786), [PRJNA1254145](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1254145). |
| Human gut-derived C. albicans isolates, Nature 2022 | Strain genomes linked to measured immune-damaging phenotypes and human/mouse studies. A resource for research strain comparison, not an off-the-shelf risk classifier. | [Paper](https://www.nature.com/articles/s41586-022-04502-w), [isolate WGS PRJNA702809](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA702809), [ITS PRJNA610042](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA610042), [reference-strain PRJNA432884](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA432884). Read the linked author correction before reproducing phenotype mapping. |
| CGF human metagenome reanalyses | More than 11,000 Chinese/non-Chinese metagenomes analyzed in the publication. Useful source discovery and external disease/healthy contexts. | [CGF paper/repository](https://github.com/yexianingyue/Cultivated-Gut-Fungi). Retrieve the paper's project/sample manifest; do not assume all samples are healthy, independent, public or method-compatible. |
| Avershina 2025 fungal mocks | Known genomes, community definitions and reproducible simulation. Baseline taxonomy/abundance testing; expand to much lower fungal fractions and additional contaminants. | [Workflow/data definitions](https://github.com/Rounge-lab/mock_mycobiome). |
| Physical metagenomic mock with fungi | Example Zymo D6300 sequencing run SRR12324253; verify the exact lot, composition and extraction/library protocol. | [Run](https://www.ebi.ac.uk/ena/browser/view/SRR12324253). It tests only the fungal members actually included, not the entire mycobiome catalog. |

At reference construction, query official SRA/ENA metadata to produce an accession-versioned run manifest with file checksums, paired/single status, assay, read length, enrichment, extraction, host age, participant ID and study ID. Do not identify a cohort from a project accession alone. Public reads and restricted clinical metadata have different permissions; exclude unavailable phenotype fields from claims.

Verified starting run examples: PRJNA925700 contains 143 WGS runs; SRR23378988 supplies paired FASTQs. SRR5935746 is an iHMP WGS example. Do not trust a single archive field blindly: CHILD run SRR19421665 is labeled SINGLE in retrieved metadata while supplying `_1`/`_2` FASTQs; validate mate identifiers and file contents before assigning layout. The twin project has 424 amplicon runs, not 424 shotgun metagenomes. Freeze the queried metadata and checksums because archive counts/annotations can change.

### Reference percentile procedure

1. Define the intended population and exclusion/stratification variables before model fitting: age, site, protocol, medications, major disease criteria and repeat sampling.
2. Reprocess raw eligible reads using the identical locked module. Never compare a new genome-normalized shotgun quantity with published ITS percentages.
3. Give each participant appropriate weight. Use participant/study-grouped resampling and validation; do not treat 143 samples from six people as n=143 independent controls.
4. Report the effective eligible participant count and detection-censoring policy. Zero/undetected values need explicit handling, not a pseudocount masquerading as measurement.
5. Use a midpoint empirical percentile, `100 × (n_less + 0.5*n_equal)/n_eligible`, with uncertainty and censoring-aware interpretation. Do not label percentile ≥95 “disease” or <5 “fungal deficiency.”
6. Qualify stability and external transport before displaying an “unusually high” flag. Otherwise label it a descriptive cohort rank. Adult, infant and adolescent references are separate.

## 9. Initial health-evidence registry

Implement the following seed interpretations with citations. Add every detected fungus to the page even if its health evidence is unknown. Source-level facts are not universal species prescriptions.

| Organism/finding | Evidence-supported interpretation | What the report must not infer |
|---|---|---|
| **Saccharomyces cerevisiae; studied S. boulardii strains** | Food-associated yeast and species containing studied probiotic strains. A primary antibiotic-associated-diarrhea trial reported 22% diarrhea with placebo versus 9.5% with yeast in 180 completing participants. [Trial](https://pubmed.ncbi.nlm.nih.gov/2494098/). | Species-level reads do not identify the probiotic strain, prove benefit, or establish deficiency when absent. |
| **Saccharomyces safety/context** | Some S. cerevisiae strains worsened colitis in mice; hospital studies document probiotic-associated fungemia in susceptible settings. [Mouse study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5994919/), [hospital study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8314839/). | Do not mark all Saccharomyces harmful or all probiotic yeast risk-free. Hospital proportions are not risk estimates for healthy users. |
| **Candida albicans** | Commonly observed in healthy people but also an opportunist. Human isolates differed in immune-cell damage; high-damaging strains aggravated experimental inflammation through candidalysin-related mechanisms. [Strain study](https://www.nature.com/articles/s41586-022-04502-w). | Presence is not candidiasis; species abundance alone does not identify the high-damaging phenotype. |
| **Mucosa-associated fungal consortium** | A defined consortium promoted barrier protection/type-17 responses in mouse experiments, illustrating potentially beneficial fungal-host interactions. [Primary study](https://pubmed.ncbi.nlm.nih.gov/35176228/). | A consortium effect cannot be assigned to each species independently or treated as demonstrated human benefit. |
| **Candida tropicalis** | Associated with Crohn disease in a family study; interaction with E. coli and Serratia marcescens enhanced biofilm formation in vitro. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5030358/). | Co-detection is not an observed intestinal biofilm or proof of disease causation in this person. |
| **Candida parapsilosis complex** | Opportunistic potential; C. parapsilosis promoted diet-induced obesity through lipid-related mechanisms in mice. [Experimental study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8546080/). | Do not translate it into an individual human obesity probability, or assign the mechanism to unresolved C. orthopsilosis/metapsilosis. |
| **Nakaseomyces glabratus / Candida glabrata** | Clinically important opportunistic yeast. [Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5478&mode=info), [clinical epidemiology](https://academic.oup.com/cid/article-abstract/83/1/e109/8653119). | No validated healthy stool abundance threshold follows from bloodstream-infection data. |
| **Debaryomyces hansenii** | Food-associated yeast enriched in injured/inflamed intestinal tissue; impaired wound repair in mice. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10114606/). | Stool presence cannot establish tissue localization, impaired healing, or the need for eradication. |
| **Malassezia restricta** | Found in healthy stool and associated with Crohn mucosa/host CARD9 context; worsened colitis in susceptible experimental models. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6417942/). | Do not equate a stool hit with Crohn disease, skin disease or invasive fungal infection. |
| **Malassezia globosa and related species** | Interpret the resolved species, sample site and study context; retain potential skin/handling contribution and opportunistic context. | Tumor-tissue/mouse research is not a stool cancer-risk threshold. [Tissue/model study](https://pubmed.ncbi.nlm.nih.gov/31578522/). |
| **Clavispora lusitaniae, especially studied P4013B** | Human IBD depletion association and mouse colitis protection involving indole-3-ethanol/AHR; a promising preclinical beneficial strain. [2025 study](https://www.nature.com/articles/s41467-025-64914-w). | Species-level detection is not P4013B identification, proven human benefit or a supplementation recommendation; opportunistic behavior can coexist. |
| **Cladosporium sphaerospermum** | Depleted in Crohn ileal mucosa but unchanged in feces; protective mechanisms studied in mice/cells. [2026 study](https://pubmed.ncbi.nlm.nih.gov/41501528/). | Do not classify all Cladosporium as bad mold or use fecal non-detection to infer missing mucosal protection. |
| **Aspergillus / Penicillium** | Food/environmental passage is a plausible explanation for some stool signals; species/complex specificity matters. [Controlled-diet study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5874442/), [longitudinal culture study](https://doi.org/10.3389/fmicb.2019.01575). | DNA does not establish mycotoxin production/exposure, viability, colonization or invasive disease. Closely related food/toxin-capable organisms must not be conflated. |
| **Rhodotorula mucilaginosa** | Early-life experimental effects on metabolic outcomes differ from other fungi; retain opportunistic and host-context evidence. [2025 mouse study](https://www.nature.com/articles/s41467-025-56743-8). | Do not use a mouse phenotype as a universal human good/bad label. |
| **Other cultured or unnamed fungal taxa** | Show identity, provenance, exposure plausibility and explicit evidence gaps. | Isolation from a healthy stool sample does not prove a species is universally beneficial or a necessary resident. |

Additional 2025–26 work on fungal metabolism, host genetics and early-life maturation should be indexed as emerging evidence, not a ready-made clinical score: [host genetic associations/Mendelian randomization](https://doi.org/10.1371/journal.pbio.3003339), [infant maturation/allergy](https://www.nature.com/articles/s41467-026-74418-w), [Fusarium experimental metabolic protection](https://pubmed.ncbi.nlm.nih.gov/40310917/), [filamentous-fungus radioprotection in experimental models](https://doi.org/10.1073/pnas.2608386123). A Mendelian-randomization inference retains its assumptions and is not equivalent to experimentally transferring disease or a validated stool risk prediction.

### 9.1 Residency and dietary context

Use `origin_interpretation = unknown | compatible_with_food_or_supplement | compatible_with_oral_or_skin_source | recurrent_detection | viable_isolate_confirmed | mucosal_localization_confirmed` as nonexclusive evidence fields. A single stool sample cannot prove long-term residence.

The controlled-diet study supporting food/oral passage was small; it does not establish that all intestinal fungi are transient. Longitudinal culture and mucosal studies provide other evidence. Preserve both findings rather than choosing a universal colonization narrative. [Auchtung 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5874442/), [Raimondi 2019](https://doi.org/10.3389/fmicb.2019.01575).

### 9.2 Evidence labels

Each claim states the organism/strain, host species, anatomical site, outcome, study design, direction, replication, primary citation and inference limit. Allowed short labels: `human_trial`, `human_association`, `human_isolate_plus_experiment`, `animal_experiment`, `in_vitro`, `genome_inference`, `insufficient_evidence`.

The same fungus may have protective and harmful effects across endpoints or contexts. Neither missing clinical trials nor positive in-vitro findings justify deleting the evidence or overstating it. Keep mechanistic/preclinical options visible with the correct label.

## 10. Integration with other OpenBiota results

**Bacterial ecosystem:** show relevant co-occurrence and literature context, but do not infer a microbial interaction network from one person's single sample. Shared presence of Candida/E. coli/Serratia may prompt a research note; it is not an actual measured biofilm.

**Disease profiles:** route observed fungal features only into models trained and validated with compatible fungal measurements. Do not add fungal percentages to existing MP3-trained bacterial disease scores. A research association can be displayed independently without manufacturing an adjusted diagnostic score.

**FMT matching:** run §6.7 for supplied donor/recipient samples and include fungal strains, supported determinants and unresolved concerns. A donor's higher fungal diversity or lower Candida percentage is not by itself proof of superior suitability. Do not use the new page to issue donor clearance or to assert that a fungal community will engraft. Provide sample-specific analytical resolution alongside every comparison.

**Biofilms:** link fungal adhesion/hyphal/biofilm-capacity findings to the existing genetic-potential module while retaining its measurement limits. Stool DNA does not visualize adherent fungal biofilms or establish their location/volume/activity.

**Interventions:** use the existing evidence store and exact taxon/strain/indication matching. Do not recommend antifungals automatically because fungi were detected, or a probiotic automatically because a potentially beneficial fungus was absent. Preserve preclinical evidence as research options for review, without converting it into proven human efficacy.

## 11. Result-specific actions

| Result | Appropriate next-step text |
|---|---|
| Low fungal support/no supported calls | Explain achieved sensitivity and what was searched. If fungal symptoms/clinical concerns persist, a targeted fungal assay or appropriately collected follow-up sample may answer questions that low-depth shotgun cannot. |
| Saccharomyces or food-associated signal with relevant exposure | Record yeast probiotic/fermented-food exposure and collection timing. Consider a standardized follow-up when the clinical question is persistent colonization; do not prescribe stopping a beneficial therapy from the sequence result alone. |
| Technically supported high reference-relative opportunist signal | Explain the magnitude, compatible reference and context; recommend clinical correlation/targeted confirmation where symptoms, immune status or persistence make it relevant. “Higher than reference” is not “infection confirmed.” |
| Fungal signal dominated by one conserved locus, near-neighbor ambiguity or batch contamination | Say the result needs analytical confirmation; do not generate an organism-specific treatment suggestion from the uncertain species label. |
| A studied beneficial strain is genuinely resolved | Display its actual human/preclinical evidence and indication; stool presence is not evidence that the host is receiving the studied therapeutic effect. |
| Preclinical beneficial species detected, strain unresolved | Show promising research with “studied strain not established in this sample.” Do not turn it into a commercial supplementation recommendation. |
| Persistent/recurrent supported finding across comparable samples | Show the time series and exposure changes. Recurrence strengthens evidence of persistent detection, but repeated ingestion can still contribute. |
| Suspected small-intestinal or invasive fungal disease | State that this stool module does not establish SIFO, mucosal invasion, fungemia or systemic fungal infection. Direct the finding to the appropriate clinical evaluation rather than inventing a diagnosis. |

Provide no eradication target for all fungi and no fixed “anti-Candida” regimen generated solely from relative abundance. If therapeutic evidence is retrieved, include strain/compound specificity, experiment type, delivery relevance, human evidence where available and clinically relevant risks; do not hide nonclinical studies or misrepresent them as successful treatment in this individual.

## 12. Validation and release qualification

### 12.1 Truth sets

Reproduce the public 2025 fungal mock workflow, then expand it. Include independent/held-out strains, related species, mixed strains, hybrid/heterozygous references, incomplete assemblies, bacterial/human contamination in reference genomes, food/plant DNA, synthetic-negative libraries and physical mocks. Match truth to quantity: retained read-origin truth for sequence fractions, input mass/cell composition separately for extraction/library bias.

Test 8M, 14M and 30M paired-fragment libraries where available or explicitly simulated, with total fungal sequence fractions of 0%, 0.0001%, 0.001%, 0.01%, 0.1% and 1%. Include one-, three-, ten- and higher-diversity fungal mixtures and unequal species abundance. Do not promise that all species are recoverable at every tier. Do not duplicate real reads to manufacture depth.

Split calibration from evaluation by strain/community and study. Compare old MetaPhlAn-only fungal calls, FungiGut reproduction, corrected competitive genome mapping, EukDetect2 and Metax on identical truth, with optional additional profilers. A publication's claims require replication in the intended workflow before becoming product performance claims.

### 12.2 Engineering gates

These are proposed analytical release targets, not existing clinical standards:

- Species/complex precision ≥99% within the claimed validated operating range, with uncertainty and per-taxon near-neighbor results. Where a species cannot meet specificity, report the supported complex/rank rather than falsely naming it.
- Define the claimed LoD per organism/complex and protocol as a signal range with at least 95% detection in appropriate repeated challenges, accompanied by its confidence interval; do not extrapolate a shallow simulation into a clinical LoD.
- Define a separate LoQ based on measured reproducibility and abundance error. A starting engineering criterion is median absolute relative error ≤25% and replicate CV ≤25% within the claimed range; validate bias from genome size, extraction and mixture context separately.
- No unexplained supported fungal call in no-fungus computational controls. Physical blanks/controls must be investigated for real contamination as well as algorithmic false positives, without retroactively changing truth to make a test pass.
- Repeated execution and reference/input ordering changes produce identical supported taxa, counts and deterministic primary estimates. Stochastic research reproductions store their seed and are not production truth.
- Fungal signal and within-fungi mass conservation pass exactly under the declared denominators; unrelated bacterial abundance changes cannot alter raw fungal fragment counts.
- New genes/variants require independent locus-level validation and species attribution; taxonomic sensitivity does not validate resistance or virulence phenotype prediction.

Report precision/recall against all eligible truth taxa, including misses. Evaluate new calls separately, not hidden within overall bacterial accuracy. False assignments to unnamed fungal clusters count as false positives too. The five project samples are migration examples, not a sufficient validation or healthy-reference cohort.

### Required strain validation

Benchmark known isolate genomes/raw reads and real stool strain examples independently of the species benchmark. Include reference-present and reference-absent strains, near-identical commercial/food isolates, haploid and heterozygous diploid populations, LOH/CNV/aneuploidy, and single versus mixed strains. Test 100:0, 95:5, 80:20 and 50:50 mixtures across a fungal coverage ladder (for example 0.1×, 1×, 5×, 10×, 20× and 50×); these are engineering challenges, not guaranteed operating ranges or hard reporting cutoffs.

Report correct clade assignment, correct reference/equivalence assignment, unsupported named-strain assignment, unknown-strain rejection, genotype error, callable genome fraction, minor-strain detection, determinant sensitivity/specificity and pairwise relationship accuracy separately. A starting named-reference precision target is ≥99% within the declared range, with confidence intervals and near-neighbor results; do not turn that target into “99% accurate” marketing without adequate independent validation. Define per-panel callability/minor-strain limits and evaluate error versus unresolved-rate so an implementation cannot pass by refusing every call.

Use both downsampled real isolate reads and realistic simulated communities with errors, contamination and mapping competitors. Keep the source individual/strain and close related isolates out of conflicting training/test partitions. Also evaluate multiple independent isolates from the same named lineage rather than memorizing one assembly. Test genotype comparison over the actual shared callable mask; “same strain” thresholds must be learned/qualified separately for each organism and claim type.

### 12.3 Essential tests

| ID | Required check |
|---|---|
| MYC-01 | 16S-only input returns incompatible assay, not no fungi. |
| MYC-02 | Total-sample fungal percentage requires an original denominator and pre-depletion reads or validated lossless fungal retention; the denominator alone cannot repair lost fungal reads. |
| MYC-03 | 100% FungiGut composition is not displayed as 100% fungus in the microbiome. |
| MYC-04 | EukDetect eukaryotic output excludes nonfungal taxa from fungal totals. |
| MYC-05 | RPKS/RPKSB/RelEuk/native fields retain correct denominators and units. |
| MYC-06 | `F/N` and fungal fragments per million agree; unknown sequence remains in primary N. |
| MYC-07 | Each pair contributes at most one fragment to mutually exclusive categories. |
| MYC-08 | Cross-kingdom ambiguous reads cannot become confident fungal reads. |
| MYC-09 | Mitochondrial/rDNA-only support is distinguished from nuclear genome-wide evidence. |
| MYC-10 | Human/bacterial contamination in a fungal reference cannot create a supported fungal call in the negative fixture. |
| MYC-11 | Removing an organism from one database while retaining it in another does not qualify as an open-set test. |
| MYC-12 | Strain duplicates, synonyms and unresolved ancestor/descendant overlap do not inflate richness or alter effective abundance denominator. |
| MYC-13 | Food-associated/clinical near-neighbor complexes are not forced to a named species. |
| MYC-14 | FungiGut min_map is parsed as alignment length, not MAPQ. |
| MYC-15 | Mismatches encoded within CIGAR M are included in real identity calculations. |
| MYC-16 | Permuting SAM optional-tag order does not change parsing/results. |
| MYC-17 | Shuffling/collating read records and paired mates does not change fragment results. |
| MYC-18 | Unseeded multimapper randomness is absent from production estimates. |
| MYC-19 | Supported qualitative detections survive low quantitative precision with the proper label. |
| MYC-20 | Missing database/execution becomes not assessed, not zero fungi. |
| MYC-21 | Fraction types are [0,1]; calibrated absolute values use their own nonnegative units. |
| MYC-22 | Generic S. cerevisiae does not become a verified S. boulardii product strain. |
| MYC-23 | Generic C. albicans does not become a verified high-damaging isolate. |
| MYC-24 | Gene presence does not become expression, biofilm activity, toxin exposure or drug susceptibility. |
| MYC-25 | Same fungus may retain both beneficial and concerning evidence without forced cancellation. |
| MYC-26 | Mouse/cell/association/clinical-trial evidence labels remain distinct in rendered text. |
| MYC-27 | Less than 1% is never used as a universal green reference interval. |
| MYC-28 | No detection does not become a fungal deficiency, clearance or reassurance score. |
| MYC-29 | Repeated cohort samples retain participant grouping. |
| MYC-30 | ITS percentages cannot calibrate WGS total fungal signal. |
| MYC-31 | Infant/adult references cannot silently supply adolescent percentiles. |
| MYC-32 | Enriched fungal libraries cannot calibrate unenriched whole-stool fungal fractions without a validated adjustment. |
| MYC-33 | Physical mock depth/mixture assumptions are recorded and not padded with duplicated reads. |
| MYC-34 | LoD, LoQ and low-signal uncertainty are independently handled. |
| MYC-35 | Old bacterial/model outputs remain unchanged for unchanged inputs/reference locks. |
| MYC-36 | Fungal taxa appear in the overall inventory without duplicate count or duplicated abundance mass. |
| MYC-37 | FMT/biofilm integration does not convert detection into engraftment, donor clearance or a visualized biofilm. |
| MYC-38 | Every supported fungus and unresolved complex has a readable evidence card. |
| MYC-39 | Page renders a useful insufficient-signal state instead of disappearing. |
| MYC-40 | PDF/HTML labels distinguish total signal from fungal composition and fit without truncation. |
| MYC-41 | External HTML links use target="_blank" rel="noopener noreferrer". |
| MYC-42 | All database sequences/code have version, source, hash and applicable license provenance. |
| MYC-43 | Citation subject, strain, host, body site and conclusion match the underlying study. |
| MYC-44 | Absolute qPCR/culture values remain separate from sequence fractions with calibration metadata. |
| MYC-45 | Project migration audit shows actual gained/lost fungal calls and evidence, not a promised richness increase. |
| MYC-46 | Metax fractional-index modeled counts cannot replace observed fungal-fragment counts. |
| MYC-47 | Unresolved fragment percentages cannot be inserted into genome-equivalent composition without a compatible estimator. |
| MYC-48 | Unphased mixed-strain alleles cannot manufacture a named strain or combined resistance genotype. |
| MYC-49 | Threshold selection cannot access final held-out benchmark outcomes; locked evaluation remains reproducible. |
| MYC-50 | Default analysis schedules strain assessment for every fungal finding; species-only success cannot satisfy full release completeness. |
| MYC-51 | Missing strain software/reference execution is an explicit incomplete/error state, never no harmful strain detected. |
| MYC-52 | Adequately covered known-strain fixtures resolve correctly; the test suite cannot pass by making every strain unresolved. |
| MYC-53 | Held-out unrepresented strains cannot become a named reference merely because it is the nearest available genome. |
| MYC-54 | Indistinguishable reference strains produce an equivalence group, not arbitrary product/isolate attribution. |
| MYC-55 | One heterozygous diploid isolate is not automatically split into two strains. |
| MYC-56 | Distinct conspecific mixtures retain supported minor alleles and are not collapsed into a chimeric named consensus strain. |
| MYC-57 | Zero coverage is missing data, not a reference allele, gene deletion or absence of a resistance variant. |
| MYC-58 | Strain markers undergo cross-reference near-neighbor/decoy tests and assembly/raw-read verification. |
| MYC-59 | A few matching covered bases or 100% shared-allele ANI cannot establish exact strain identity. |
| MYC-60 | Reference orientation/liftover preserves REF/ALT, chromosome and phase; a paper's descriptive coordinates cannot directly become an unchecked VCF. |
| MYC-61 | Diploid allele dosage, LOH, CNV and multiple-strain explanations are tested independently. |
| MYC-62 | Correct fungal genetic code/exon structure is used for every protein variant; nuclear/mitochondrial annotations cannot mix. |
| MYC-63 | Supported determinants remain reportable with carrier-strain unresolved, without inventing linkage or a combined genotype. |
| MYC-64 | Donor/recipient shared genotypes yield unresolved source where both were present; post-FMT detection alone is not durable engraftment. |
| MYC-65 | Variant/lineage identification does not silently become demonstrated candidalysin activity or a high-damaging phenotype. |
| MYC-66 | FungAMR negative evidence and evidence-grade direction are preserved; novel changes at a hotspot remain unvalidated. |
| MYC-67 | Probiotic-lineage marker support does not identify a commercial product or prove its studied benefit occurred. |
| MYC-68 | Strain cache invalidates on reference panel, mask, marker, ploidy, annotation/genetic-code or algorithm changes. |
| MYC-69 | Report shows callable discriminatory loci, reason for unresolved strain calls and a technically relevant next measurement. |
| MYC-70 | Closely linked marker sites are not counted as independent evidence, and phase cannot extend beyond supported blocks. |
| MYC-71 | Low-depth locus assessment runs even when whole-genome assembly or a genome-wide strain claim is impossible. |
| MYC-72 | Public archive assay fields and actual reads distinguish isolate WGS, stool WGS and ITS; papers' availability statements alone cannot label a strain truth set verified. |
| MYC-73 | The overview and main report summary render MHS-E1 on one 0–100 gauge; JSON also serializes the exactly corresponding 0–1 value. |
| MYC-74 | At least one adequate positive and one concern fixture produce non-null scores through the full observation-to-rule path; an always-null implementation fails. |
| MYC-75 | Property tests over admissible B/C inputs show score bounded by 0–100, nondecreasing in B, nonincreasing in C and deterministic. |
| MYC-76 | No fungal calls, unknown-only findings or zero active directional contributions yield null with a specific reason; they do not yield 0, 50 or 100. |
| MYC-77 | B=0.50/C=0 yields 75; B=0/C=0.50 yields 25; B=C=0.50 yields 37.5 before rounding. |
| MYC-78 | Duplicate taxa aliases, reference genomes, genes or publications do not change contributions or interpretation-coverage fragment counts. |
| MYC-79 | Adding unknown organisms or unscored neutral evidence leaves the point score unchanged; unknowns remain visible and coverage is recomputed. |
| MYC-80 | For C≥0.50, no B≤1 can raise the score above 50; adding a probiotic cannot hide a concern card or alert. |
| MYC-81 | Generic Candida/ECE1/candidalysin variants do not activate the HD-isolate rule; the specific positive HD fixture does. |
| MYC-82 | A compatible HD/LD equivalence set widens C and the score range; an unresolved higher-resolution claim is never silently scored as a known safe strain. |
| MYC-83 | P4013B reference support activates B=0.50; species-only support activates only the explicit B=0.25 extrapolation; generic PDC alone does neither. |
| MYC-84 | An old S. boulardii trial without a verified genotype-to-trial join does not activate the named probiotic rule; species identity and supplement history are not substitute joins. |
| MYC-85 | High fungal fraction, high diversity, normal-range abundance and non-detection of a fungus cannot independently award beneficial points. |
| MYC-86 | A high-abundance rule needs the qualified source/reference contract; at p=0.50 activation=0 and at p=0.90 activation=0.80. Its absence does not block other rules. |
| MYC-87 | Sensitivity envelopes contain the point; resolution × B-policy × C-policy combinations are evaluated jointly. B=0.50/C=0 gives 62.5–87.5; with C unresolved over 0–0.50 the combined range is 15.625–87.5. |
| MYC-88 | Range labels are sensitivity ranges, not statistical CIs. Midpoint labels never say normal, and a range crossing 50 displays direction uncertain. |
| MYC-89 | Existing pathogen/resistance alerts survive a favorable score and remain visible; a max aggregation does not suppress additional non-dominant concerns. |
| MYC-90 | Downsampling that loses a concerning call must worsen information status or range and cannot be presented as improved health; version/depth-incompatible time series are not ranked as improvement. |
| MYC-91 | Failed mandatory analysis is an incomplete/error state, whereas completed unresolved strain analysis can support a lower-resolution score. No cache silently skips mandatory work. |
| MYC-92 | Every contribution exposes source, identity/context predicate, cap, activation and deduplication family; all score-relevant changes invalidate the cache. |
| MYC-93 | Interpretation confidence remains experimental_limited regardless of sequencing depth; coverage is not represented as the percentage of health understood. |
| MYC-94 | Adult-reference abundance rules cannot activate on minors; donor union/pooling does not inherit the best donor's score or generate clearance. |
| MYC-95 | Renderer links use target="_blank" rel="noopener noreferrer" and the gauge has a text alternative; score is not conveyed by color alone. |
| MYC-96 | Removing the population-reference dataset still permits all bound presence/genotype/mechanistic seed rules to execute and score qualifying fixtures. |

## 13. Data contract and executable handoff

Implement a strict schema with these required fields; use existing project conventions for serialization/storage. The sample below is a **schema illustration**, not a measured participant result.

```json
{
  "schema_version": "openbiota.mycobiome.v3",
  "sample_id": "example",
  "assay": "stool_shotgun_dna",
  "reference_lock_id": "content-addressed-lock",
  "input_manifest_id": "content-addressed-input",
  "analysis_status": "not_assessed",
  "measurement": {
    "denominator_stage": "qc_nonhost_before_bacterial_filter",
    "eligible_fragments": null,
    "fungal_supported_fragments": null,
    "cross_kingdom_ambiguous_fragments": null,
    "unassigned_fragments": null,
    "fungal_fragment_fraction": null,
    "fungal_fragments_per_million": null,
    "sampling_interval": null,
    "quantification_status": "not_assessed"
  },
  "within_fungi": {
    "estimator_id": null,
    "quantity_type": null,
    "denominator_id": null,
    "unresolved_composition_fraction": null,
    "unresolved_fragment_fraction": null,
    "unresolved_fragment_denominator_id": null,
    "composition": []
  },
  "taxa": [],
  "strain_analysis": {
    "required": true,
    "status": "not_assessed",
    "panel_lock_id": null,
    "taxa_eligible": null,
    "taxa_assessed": null,
    "taxa_with_resolved_reference_or_lineage": null,
    "execution_failures": [],
    "pairwise_comparisons": []
  },
  "beneficial_evidence_findings": [],
  "potential_concern_findings": [],
  "health_score": {
    "model_id": "MHS-E1",
    "policy_lock_id": null,
    "status": "not_assessed",
    "reason_codes": [],
    "score_01": null,
    "score_100": null,
    "benefit_contribution": null,
    "concern_contribution": null,
    "resolution_sensitivity_range": null,
    "policy_sensitivity_range": null,
    "display_sensitivity_range": null,
    "range_type": "scenario_sensitivity_not_confidence_interval",
    "direction": "not_assessed",
    "mixed_signals": null,
    "analytical_confidence": "not_assessed",
    "interpretation_confidence": "experimental_limited",
    "scored_fungal_fragment_fraction": null,
    "strain_assessment_completion": null,
    "contributions": [],
    "unscored_findings": [],
    "unresolved_concerns": [],
    "active_alert_ids": [],
    "context_missing": [],
    "reference_percentile": null,
    "clinical_validation_status": "not_validated_for_global_health"
  },
  "exposure_context": [],
  "reference_comparisons": [],
  "next_steps": [],
  "limits": [],
  "validation_release_id": null
}
```

Every taxon record includes canonical identity, accepted/synonym names, supported rank, all source observations, read/marker/coverage evidence, detection and quantification states, the required strain record below, gene/variant assessment results, primary abundance with units/denominator, mapping/ambiguity reasons, source URLs and interpretation claims. Use null for unavailable numeric data; no fabricated zeros or patient results.

For `health_score`, enforce `status = not_assessed | computed | provisional | not_computable | failed`; `direction = not_assessed | favorable_evidence | concerning_evidence | uncertain`. Status describes computation/completeness, not validation: a `computed` MHS-E1 remains experimental. Favorable/concerning requires the displayed sensitivity interval wholly above/below 50; otherwise direction is uncertain. Before computation the confidence field may be `not_assessed`; afterward use §7.5.5 values. Every range is an object with `low`, `high`, `scale` (100), `scenario_ids` and `method_id`. `contributions` contains the rule fields in §7.5.2, plus actual value, admissible alternatives, activation status and the reason a rule did or did not contribute. `reference_percentile` stays null until separately qualified. Migrate existing v2 objects to v3 with a not-assessed score, then compute from retained evidence; never backfill a synthetic midpoint. No participant result appears in the schema example.

Required per-taxon strain object, again an unmeasured schema illustration:

```json
{
  "analysis_status": "not_assessed",
  "resolution": "not_assessed",
  "reason_codes": [],
  "panel_lock_id": null,
  "reference_accessions": [],
  "reference_equivalence_group": [],
  "lineage_id": null,
  "genotype_fingerprint_id": null,
  "ploidy_model": null,
  "genetic_code_id": null,
  "callable_nuclear_bases": null,
  "eligible_discriminatory_loci": null,
  "callable_discriminatory_loci": null,
  "supporting_loci": [],
  "contradictory_loci": [],
  "mixture_status": "not_assessed",
  "mixture_components": [],
  "phase_blocks": [],
  "determinants": [],
  "phenotype_claims": [],
  "evidence_artifact_ids": [],
  "validation_release_id": null,
  "recommended_analytical_followup": []
}
```

Enforce `analysis_status = not_assessed | complete | partial | failed`; `resolution = not_assessed | unresolved | species_only | lineage | reference_equivalence_group | reference_genotype | sample_genotype`; `mixture_status = not_assessed | no_mixture_evidence_at_achieved_sensitivity | mixture_supported | ambiguous`. These are different axes: a complete analysis may correctly yield an unresolved strain, whereas an execution failure is not a complete analysis. Validate the reason codes against §6.4.4. Every determinant stores its own callability, carrier resolution, alleles/copy evidence, reference coordinates, applicable drug/trait and source-backed claim status.

The reference lock includes code commit/container digest, source release, sequence/file hashes, accession versions, reference masks, fungal/nonfungal taxonomy and decoy membership, actual admitted counts, licenses, parameters, validated operating range and calibration version. Strain caches additionally key on pangenome/coordinate backbone, marker panel, callable mask, ploidy and genetic-code settings, genotype/mixture model, phenotype mapping and sample inputs. Native output snapshots remain inspectable.

Score caches additionally include model/policy/rule/source hashes, all phenotype joins, input observations and callability, assay/host-context fields, reference-cohort lock where used, unresolved hypotheses and confidence/range algorithms. A score-only recomputation can reuse compatible expensive alignments. Different models or policies must not compare silently.

Implement CLI equivalents:

```bash
openbiota mycobiome references prepare --preset v08.2 --output /db/openbiota/mycobiome-v08.2
openbiota mycobiome benchmark --suite v08.2 --split heldout
openbiota mycobiome strains benchmark --suite v08.2 --split heldout
openbiota mycobiome score benchmark --model MHS-E1 --suite v08.2
openbiota mycobiome analyze --sample SAMPLE2_A02 \
  --r1 SAMPLE2_A02.R1.fastq.gz --r2 SAMPLE2_A02.R2.fastq.gz \
  --reuse-compatible-cache --output results/SAMPLE2_A02
openbiota mycobiome strains compare --recipient results/SAMPLE2_A02 \
  --donor results/SAMPLE4_A04 --donor results/SAMPLE6_A06 \
  --output results/fungal_strain_comparison
openbiota mycobiome score --sample-results results/SAMPLE2_A02 --model MHS-E1
openbiota report render --sample-results results/SAMPLE2_A02 --include-mycobiome
```

These are commands to implement within the existing CLI, not assertions that they currently exist. The embedded preset must contain the references and policies in this document; there is no missing companion configuration handoff. Runtime reference locks, data files, model artifacts and tests are generated implementation assets.

The ordinary `mycobiome analyze` command runs the required strain, applicable determinant and MHS-E1 scoring stages automatically; users do not need a hidden opt-in flag. The separate score command supports recomputation from existing results. The compare command consumes those results and performs compatible reanalysis when needed. A progress/partial report can render actual findings, but it must show incomplete stages rather than declare a successful full analysis.

Release in this order: (1) reference, strain-panel and software-adapter audit; (2) fungal-specific detection/denominator parser; (3) mandatory strain, determinant and pairwise-comparison implementation; (4) species and strain benchmark qualification; (5) MHS-E1 scoring, traceability, sensitivity and tests, complete evidence cards and summary page; (6) matched-cohort percentiles and quantitative-assay extensions where qualified. A missing optional population cohort must not block reporting actual fungal detections, strain evidence or eligible MHS-E1 contributions.

**Definition of done:** the report shows what fungal sequence signal was measured, which species and strains/lineages were supported, which relevant genetic determinants were assessed, the strain resolution achieved for every finding, why a finding may be beneficial or concerning, and what follow-up would resolve remaining analytical questions. It includes the single MHS-E1 gauge with exact score, sensitivity, confidence, supporting contributions and visible limits; adequate benefit/concern fixtures must actually score. The initial release must demonstrate successful strain resolution on adequate positive fixtures and correct handling of inadequate/ambiguous data. It must not confuse fungal DNA with living biomass, within-fungal percentages with whole-sample percentages, genetic identity with proven biological activity, a policy index with measured health, or study associations with an individual diagnosis.
