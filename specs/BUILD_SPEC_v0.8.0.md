# BUILD_SPEC_v08.0 — Gut biofilm mechanisms, protective potential, harmful potential, and intervention evidence

**Project:** OpenBiota  
**Specification version:** 08.0.1 — research and scoring revision  
**Research cutoff:** 18 September 2026  
**Deliverable:** This single document is the complete implementation handoff. All decisions, evidence, data sources, scoring rules, contracts, and acceptance criteria required for this increment are included here. External research and database downloads are source inputs, not additional handoff documents.

**Revision objective:** Produce useful, reproducible experimental scores and ranked actions from the available data. This revision replaces the previous two-group minimum with explicitly scoped single- or multi-mechanism rankings; adds separately labeled community-pattern and protective-ecology scores; adds human-derived functional-phenotype evidence; and ranks intervention relevance while preserving contrary results. It does not require human clinical validation before displaying an analytically supported research measurement. It does require that the label identify what the number actually measures.

**What changes for the user:** A DNA-only report can now show up to four quantitative research cards under two headings—concerning potential and protective support—plus a ranked list of findings and evidence-linked actions. A sparse mechanism panel remains useful and visible. Actual attached biofilm, host injury, dispersal and treatment response remain distinct outcomes; where additional measured data exist, the report shows them alongside the DNA findings. This entire file supersedes the previous 08.0.0 text; no separate revision or research file is required by the coding agent.

## 1. Build objective and scientific conclusion

Add a complete biofilm section to the existing shotgun-metagenomics analysis and reports. Analyze existing FASTQs, compatible cached gene calls, strain calls, reference distributions, clinical metadata, and optional independently measured assays. Preserve existing disease, organism, intervention, and donor-matching features.

The product must distinguish **protective-associated biofilm potential**, **harmful-associated biofilm potential**, and **context-dependent matrix mechanisms**. A generic total-biofilm score has no reliable health interpretation. Protective and concerning mechanisms can coexist; one must never cancel the other.

Published research supports measuring selected biofilm-related genes in stool and explaining their organism-specific mechanisms. It does **not** yet establish a broadly validated stool-DNA test that quantifies all protective versus harmful intestinal biofilms, their mucosal location, or their current activity. No validated universal good/bad biofilm classifier or transferable clinical coefficient set was identified in this review. This specification creates an explicit research measurement layer and a validation path, not a claim that the missing validation already exists.

Deliver numerical research rankings whenever an analytically supported measurement and an applicable comparison exist. A single supported mechanism is sufficient for a **named, limited-scope score**; it is not sufficient to claim comprehensive gut health. Compute the additional community-pattern/ecological-support indices in §8.7 from compatible existing taxonomic profiles, with their explicit proxy labels. If no compatible reference is available, show the raw measurement and measured/assayed panel coverage while building the reference; do not fabricate percentiles. Lack of clinical trials must not suppress these experimental outputs or directly supported intervention candidates.

The optimization target is **less damaging behavior with preservation of beneficial community function**. Lower total matrix biomass is not a universal improvement. Reduction in viable harmful organisms, epithelial injury, excessive dispersal or a source-specific damaging mechanism can matter more than reducing all biofilm.

### Required result layers

| Layer | What it can establish | What it cannot establish alone |
|---|---|---|
| Stool shotgun DNA | Detected genetic capabilities; selected accessory loci; carrier/strain linkage when supported; relative feature abundance | Active expression, attached architecture, viable biomass, location, host fibrin, or treatment response |
| Stool/biopsy RNA | Expression of selected mechanisms at sampling, with RNA-specific controls | Proof of assembled matrix or an attached mucosal biofilm |
| Protein, matrix, or metabolite assay | The analyte actually measured, with assay-specific specificity | All biofilms; anatomical location from a nonspecific stool analyte |
| Endoscopy plus appropriately preserved spatial imaging | Location and architecture in the sampled tissue; proximity/encroachment when assessed | Entire-gut absence from one negative site; species identity without a supporting assay |
| Culture/phenotyping | Organism-specific behavior under the tested conditions; viable killing versus biomass changes | The same phenotype in the patient's gut without transport evidence |

Use the headline **“Biofilm-related potential”** for DNA-only results. “Biofilm state” is reserved for an explicitly described multimodal observation; it must name its specimen, site, method, and limitations.

## 2. Scope, existing application, and implementation responsibilities

The known application already has taxonomic profiles, functional panels, disease profiles, intervention evidence, pathogen/AMR features, and donor matching. The coding agent also reports locally retained FASTQs and MetaPhlAn databases; retained `*.bowtie2.bz2` summaries are not equivalent to SAM/BAM alignments. Biofilm genes frequently lie outside MetaPhlAn marker sets. Analyze the full host-filtered microbial reads for this increment rather than assuming marker alignments cover them.

This research task did not execute the application repository or the user's FASTQs. Do not interpret example objects below as measured SAMPLE2/SAMPLE4/SAMPLE6 results. Repository paths below are proposed logical responsibilities; bind them to the actual repository after inspection.

1. Inventory CLI entry points, cached data schemas, read filtering, gene search, taxonomy versions, the shared resolution service, intervention retrieval, report rendering, and donor matcher.
2. Reuse the existing strain-resolution contract where present: `resolve(sample_id, target_id, required_resolution, permitted_evidence_kinds)`. Add typed adapters instead of maintaining conflicting gene calls in two pipelines.
3. Implement a versioned biofilm registry, analysis engine, evidence-card retrieval, independent score calibration, report panels, and longitudinal comparisons.
4. Add explicit assay availability records. Missing RNA, microscopy, carrier linkage, or a reference cohort are separate limitations, not failed DNA analysis.
5. Preserve raw findings when a downstream panel cannot be scored. No “everything unavailable” fallback from one missing optional assay.
6. Never claim a clinical sample is biofilm-free, infection-free, a cleared FMT donor, or free of disease-transfer risk from this module.

### Proposed code responsibilities

```text
biofilm/registry        source-specific modules, ontology, reviewed claims, intervention links
biofilm/acquisition     versioned public reference and cohort acquisition manifests
biofilm/detection       adapters to existing gene, assembly, and strain resolution services
biofilm/interpretation  biological context and evidence-grade rules
biofilm/scoring         frozen reference transformations and independent axes
biofilm/assays          optional RNA, protein/activity, microscopy and phenotype inputs
biofilm/report         overview, full findings, explanation and intervention cards
biofilm/validation     analytical benchmarks, cohort evaluation, regression fixtures
```

These are implementation responsibilities, not additional documents the user must supply. Generated registries, database indexes, JSON, and test fixtures belong inside the application repository and caches.

## 3. The evidence that changes the design

The following findings are mandatory implementation constraints, not optional background reading.

| ID | Source and actual result | Implementation consequence |
|---|---|---|
| BF-S01 | [Baumgartner et al., 2021](https://doi.org/10.1053/j.gastro.2021.06.024): visible ileocolonic biofilms in 57% of IBS, 34% of UC and 6% of controls; 1,426 were screened and 1,112 met inclusion criteria, with 212 biofilm-positive participants overall. Molecular analyses used smaller subsets. | Human mucosal biofilms matter, but this is not a validated stool-shotgun classifier. Preserve site and assay. Association is not universal causation. |
| BF-S02 | [Buret and Allain, 2023](https://doi.org/10.1084/jem.20221743): synthesis of host regulation and protective-community disruption. | Dissolving every biofilm is not a desirable treatment objective. |
| BF-S03 | [Jandl et al., 2024](https://doi.org/10.1128/cmr.00133-23): intestinal biofilms, host defense and therapeutic opportunities. | Useful evidence map; its review statements do not become independently validated diagnostic weights. |
| BF-S04 | [Biophysical determinants, 2021](https://doi.org/10.1016/j.cobme.2021.100275): flow, mucus and spatial conditions shape gut biofilms. | DNA alone cannot reconstruct thickness, attachment or location. |
| BF-S05 | [Hillman et al., 2025](https://doi.org/10.1016/j.gastha.2025.100712): shotgun stool study of 26 bile-acid-diarrhea cases and 21 controls, with higher biofilm-associated gene abundance and altered bile-acid metabolism. | Add the exact exploratory BAD-associated gene pattern. It is a disease-associated DNA pattern, not a direct measurement of mucosal biofilm or proof of causation. |
| BF-S06 | [BAP/amyloid study, 2024](https://doi.org/10.1038/s41467-024-48309-x): accessory biofilm-associated proteins and amyloid-related experiments, stool molecular measurements, and Parkinson's metagenome reanalysis. | Add separate BAP genetic potential and measured-amyloid evidence. Neither proves brain amyloid, Parkinson's transmission, or harmful activity in every carrier. |
| BF-S07 | [BBSdb, 2024](https://doi.org/10.3389/fcimb.2024.1428784): curated biofilm-associated proteins plus transcriptomic/proteomic evidence across organisms. | Useful candidate annotations; not a person-level biofilm or good/bad classifier. |
| BF-S08 | [Klebsiella isolate study, 2024](https://doi.org/10.1038/s41522-024-00629-y): substantial phenotype variation across 100 clinical isolates. | Species abundance and one adhesion-gene hit do not establish a strain's biofilm phenotype. |
| BF-S09 | [Klebsiella biofilm evolution, 2026](https://doi.org/10.1038/s41467-026-71505-w): background- and environment-dependent genomic changes in non-gut infection models. | Curated variants can be research context, not unvalidated universal gut-risk weights. |
| BF-S10 | [Tomkovich et al., 2019](https://doi.org/10.1172/JCI124196): human biofilm-positive mucosal communities promoted tumors in susceptible mice. | Mechanistic concern is real; do not turn it into a donor-specific human transmission probability. |
| BF-S11 | [Dejea et al., 2014](https://doi.org/10.1073/pnas.1406199111): spatially organized communities associated with colorectal cancer. | Tissue location and community behavior matter; a common-genus list is insufficient. |
| BF-S12 | [Motta et al., 2019](https://doi.org/10.1038/s41467-019-11140-w): physiological epithelial thrombin contributes to mucosal biofilm control. | Host protease activity is not inherently adverse and is not inferred from DNA. |
| BF-S13 | [Excess thrombin/Crohn's study, 2026](https://doi.org/10.1080/19490976.2026.2687903): excessive thrombin induced harmful community behavior with limited taxonomic change. | Static gene counts miss important state changes; add optional protein/activity and RNA inputs without fabricating them. |
| BF-S14 | [F. prausnitzii–syndecan-1 study, 2026](https://doi.org/10.1080/19490976.2026.2665870): host glycans supported protective bacterial biofilm-related behavior in experimental colitis, with human tissue correlations. | Protective biology exists. F. prausnitzii abundance, generic vitamin genes, or residual host SDC1 DNA are not a validated protective-biofilm score. |
| BF-S15 | [L. reuteri microsphere study, 2017](https://doi.org/10.3389/fmicb.2017.00489): a specified strain/formulation showed improved probiotic properties in biofilm form. | Preserve strain and formulation; do not assume all L. reuteri or all biofilms have the same benefit. |
| BF-S16 | [Wheat-fiber community study, 2024](https://doi.org/10.1039/D4FO01294A): a nine-species consortium formed a stable biofilm on fiber in vitro. | Fiber-associated biofilm formation can be supportive; more matrix is not automatically harmful. |
| BF-S17 | [Staphylococcal fibrin model, 2025](https://doi.org/10.1016/j.bioflm.2025.100261): organism/strain/substrate-specific host fibrin recruitment. | Separate microbial fibrin-interaction capability from measured host fibrin and systemic infection. |
| BF-S18 | [Staphylococcal context study, 2026](https://doi.org/10.1016/j.bioflm.2026.100389): genetic effects differed between surface biofilms and abscess communities. | No universal positive/negative weight for a regulatory gene across sites. |

### Additional primary evidence incorporated in revision 08.0.1

| ID | Primary source and observed result | Implementation consequence |
|---|---|---|
| BF-S19 | [Prefrail human gut communities, 2025](https://doi.org/10.1038/s41522-025-00716-8): 42 donors, including 13 adults aged 30–70, 15 robust adults over 70 and 14 prefrail adults over 70. Feces-derived laboratory biofilms differed in dispersal and inflammatory effects. Grape polyphenols reduced selected adverse effects while increasing biomass. | Add community-behavior measurements and stabilization-oriented actions. Human-derived ex-vivo behavior is not observed in-vivo intestinal architecture. Age/frailty patterns cannot calibrate minors. The transfer experiment used mouse microbiota, not a human treatment trial. |
| BF-S20 | [UC biofilms and oncotraits, 2023](https://doi.org/10.1093/ecco-jcc/jjad092): 80 UC patients and 35 controls; 873 longitudinal biopsies, 265 shotgun-sequenced specimens. Sequencing concerned **biopsies**; stool oncotraits used qPCR. Biofilm presence was not significantly associated with dysplasia, whereas colibactin-marker carriage was. FadA had an inverse association in this UC cohort. | Prioritize specific damage mechanisms and context. Do not use all biofilm, all Fusobacterium or all FadA as a universal risk rule. This is valuable paired tissue evidence, not an openly reusable stool-WGS classifier. |
| BF-S21 | [Cirrhosis virulence-factor study, 2021](https://doi.org/10.1080/19490976.2021.1993584): 233 people, including 40 controls. Stool adhesion/biofilm-associated factors tracked decompensation and outcomes; changes were also studied around FMT. | Supports disease-specific functional patterns beyond taxonomy. It does not validate the same risk weights for otherwise healthy donors or prove that reducing a score changes outcomes. |
| BF-S22 | [Nine-species gut-community RNA/metabolite study, 2025](https://doi.org/10.3390/microorganisms13020234): mixed versus single-species laboratory biofilms differed in expression and metabolites. | RNA research dataset and interaction hypothesis resource. Both comparison arms were biofilms: this is **not** an active-biofilm-versus-planktonic training set and cannot be relabeled as one. |
| BF-S23 | [Opposing Lactobacillus/Klebsiella effects, 2017](https://doi.org/10.3920/bm2017.0002): L. plantarum CIRM653 reduced K. pneumoniae biofilm in vitro, but treated mice maintained intestinal carriage while controls declined. | Require viable burden/dispersal and in-vivo counterevidence on action cards. A favorable biomass assay alone must not rank this strain as a proven helpful intestinal treatment. |
| BF-S24 | [Pomegranate/barrier repair, 2023](https://doi.org/10.3390/nu15071771): improved experimental-colitis recovery and barrier-related outcomes; broad Enterobacterales biofilm inhibition was not demonstrated, with only a slight effect on C. freundii. | Include a barrier-support candidate; do not market the study as strong Klebsiella/E. coli biofilm clearance. |
| BF-S25 | [L. johnsonii against attaching/effacing pathogens, July 2026](https://doi.org/10.3389/fimmu.2026.1749001): EPEC biofilm effects in vitro and improved pathogen/inflammation outcomes in a C. rodentium mouse model. | Add a strain-qualified protective-antagonism candidate. PV739486 is a **16S sequence**, not a complete reference genome or validated strain marker panel. Putative metabolite annotations do not establish the activity of the named compounds. |
| BF-S26 | [Resveratrol/L. paracasei ATCC334, 2020](https://doi.org/10.3390/ijms21155423): enhanced adhesion and biofilm formation in vitro, with strain dependence. | Protective-support research card; not evidence that every Lactobacillus or every resveratrol exposure benefits the gut. |
| BF-S27 | [L. plantarum Y42, 2024](https://doi.org/10.1021/acs.jafc.4c00460): biofilm-state organisms and associated EPS/surface proteins supported barrier/immune outcomes in a mouse challenge model. | Broaden protective mechanisms while keeping strain, formulation and host model intact. EPS genes alone are not a unique protective signature. |
| BF-S28 | [Lactiplantibacillus LR-1, 2023](https://doi.org/10.1039/d3fo02733c): biofilm-state cells outperformed planktonic cells on selected outcomes in experimental colitis. | Supports testing functional state, not universal probiotic efficacy or a DNA activity claim. |
| BF-S29 | [Multi-strain restorative biofilm consortium, June 2026](https://doi.org/10.1007/s11274-026-05071-0): mouse studies examined biofilm formation and recovery after antibiotics. | Consortia may restore functions, but taxonomic co-presence does not establish the administered formulation, interaction or effect in a patient. Abstract-level evidence in this audit; full methods needed before importing quantitative effects. |
| BF-S30 | [Colon-targeted E. coli biofilm inhibitor, March 2026](https://doi.org/10.1021/acsmedchemlett.5c00675): experimental prodrug improved colonic delivery and selected mouse-colitis outcomes. | Demonstrates that delivery can be engineered; retain as an investigational approach, not an available supplement or proven human treatment. |
| BF-S31 | [B. subtilis protective matrix/delivery study, June 2026](https://doi.org/10.1080/19490976.2026.2684066): experimental systems support matrix-mediated stress protection and delivery. | Context for beneficial biofilm formulations; neither ordinary species carriage nor stool gene abundance establishes the studied formulation or clinical benefit. |
| BF-S32 | [Indian CRC tissue biofilms, 2025](https://doi.org/10.1007/s00253-025-13537-8): imaging of 15 tumors and 15 paired adjacent samples found heterogeneous species distributions. | Adds geographic/tissue validation context; adjacent tissue is not an independent healthy person, and imaging results are not a validated fecal signature. |

### Findings that materially improve interpretation

1. **Architecture, function and consequence are separate targets.** BF-S20 contains biopsy imaging plus tissue sequencing, with stool qPCR in a subset. Its 265 sequencing specimens must not be described as 265 stool samples or 265 independent people. Biofilm–dysplasia aOR 1.45 (95% CI 0.63–3.40) and the different oncotrait associations demonstrate why a generic presence score cannot substitute for a damaging-function assessment.
2. **Dispersal is not automatically success.** BF-S19 found greater adverse community behavior despite lower biomass in aging comparisons. BF-S23 provides a direct counterexample to assuming in-vitro disruption improves intestinal clearance. Record biomass, viable burden, dispersal and host effects separately.
3. **Protective support is experimentally real, but its sequence specificity varies.** BF-S14/S15/S25–S31 justify protective modules and action candidates. Where an exclusive functional marker or exact strain reference is absent, report ecological/strain-family support as a proxy instead of inventing a protective-biofilm gene panel.
4. **Endoscopic visibility is not microscopic presence.** Different patient selection, tissue methods, sites and phenotype definitions contribute to different prevalence estimates. Store the observed endpoint verbatim in metadata, with a controlled normalized endpoint; do not pool all labels as one truth standard.

### Especially useful new stool evidence: the BAD study

BF-S05 used MetaPhlAn 3.1 and HUMAnN 3.0.3. It examined `bssS` and `pgaA/B/C/D`; the abstract emphasizes `bssS`, `pgaA`, and `pgaB`, while the results discuss the five-gene set. Store both statements and their locations; import supplementary identifiers before selecting a production feature mapping. Record *Mediterraneibacter gnavus* as the current-name synonym of *Ruminococcus gnavus*.

Important limitations: small cross-sectional cohort; cases and controls differed in age, sex balance, recruitment site and collection period; some cases had cholecystectomy, metformin exposure, or Crohn's disease; no SeHCAT-negative IBS comparison group; no direct mucosal-biofilm measurement for all participants. Community-level gene abundance does not attribute those genes to R. gnavus. BssS is a regulator, not a one-direction “more biofilm” switch. The paper's association does not establish that increasing or decreasing any one gene causes BAD.

Implement this as **“BAD-associated biofilm-gene resemblance: exploratory”**, with a separate bile-acid context panel and no diagnosis. Do not copy published fold differences into the new profiler's reference distribution. Do not treat a species read cutoff or gene CPM from that pipeline as a validated threshold in another pipeline. The paper also repeats a KEGG transporter identifier in its text; verify supplementary mappings rather than silently assigning two functions to the duplicated identifier.

## 4. Vocabulary and the two-axis product

### Patient-facing labels and four quantitative cards

Use two main headings: **Concerning biofilm potential** and **Protective biofilm/ecosystem support**. Under them show the following distinct cards, when computable:

| Card | Default label | Measurement and interpretation |
|---|---|---|
| H-M | Harm-associated mechanisms | Named DNA-mechanism reference percentile, with single- or multi-mechanism scope and adverse-context evidence. |
| H-C | Biofilm-associated community pattern | Experimental resemblance to a specified published community association. A transferred tissue-derived signature stays visibly labeled as such. |
| P-M | Protective-associated mechanisms | Named DNA-mechanism reference percentile, with the exact strain/functional evidence and limited scope. |
| P-E | Protective ecological support | Experimental ranking of specified taxa associated with less biofilm encroachment or supportive ecology; **not a measurement of beneficial biofilm amount**. |

H-M and P-M use the mechanism specification in §8.1–8.6. H-C and P-E use the fixed experimental proxy definitions in §8.7. They must not be collapsed into a single number, mixed with each other, or silently substituted for one another. This allows useful numbers from current profiles while distinguishing a proxy from a directly supported mechanism. A positive P-E does not activate P-M.

For every item distinguish molecular confidence, carrier resolution, mechanism evidence, human-gut relevance, clinical predictive validation and sample-domain compatibility. High-confidence DNA detection may still have low predictive value for actual biofilm behavior.

### Mandatory visual rules

Display the four cards in a compact two-column arrangement under the two headings; only two columns, no more than two numeric cards in each. Use equal-size 0–100 **reference-percentile** scales when calibrated. Show scope, panel coverage and reference population alongside each number. Both high concerning potential and high protective support can coexist; do not subtract, average, divide or offset them.

Use restrained teal and amber/purple, text equivalents, and accessible symbols. Percentiles <25, 25–75 and >75 may be called “lower,” “middle” and “upper reference range,” explicitly descriptive quartiles rather than clinical normal/abnormal thresholds. Print “research index; not measured biofilm quantity” on proxy cards. Never label a score “83% bad biofilm” or “83% chance of harm.”

A one-group score is allowed and must be labeled “limited scope: [named mechanism].” Show supported/assayed and assayed/planned counts; these are coverage metrics, not health scores. An unavailable mechanism card does not hide a computable community proxy or other module. No applicable reference means a raw measurement and explicit comparison status, not an invented zero. Do not use the same numeric scale for coverage and percentiles without visible different units/labels.

Higher protective support is not an established optimum. Lower concerning potential does not exclude harmful biofilms. The actual-state strip shows separate measured imaging, behavior or inflammation observations if provided; it never gets auto-filled from DNA.

## 5. Registries and public data acquisition

Every download must record source URL, publication/accession, retrieval date, license/access conditions, release/version, file SHA-256, sequence provenance and parser version. Do not hard-code a current database release without checking the downloadable release manifest. Use openly accessible downloads; do not bypass unavailable or restricted services. An unavailable optional source must not stop the rest of the module.

### Reference resources

| Resource | Source | Use and limitation |
|---|---|---|
| BBSdb | [primary resource paper](https://doi.org/10.3389/fcimb.2024.1428784) | Candidate protein/function and experimental-context import. Count actual imported records: prose/tables report inconsistent totals. Expression changes require experiment metadata and are not DNA signs of active biofilm. |
| UniProtKB | [official API](https://www.uniprot.org/help/api_queries) | Curated protein identity, sequence versions, organism and experimental annotation; reviewed status and evidence codes must persist. Gene-name homonyms are not enough. |
| NCBI Datasets/RefSeq/GenBank | [official documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) | Genome and protein references, strain aliases, annotation coordinates and accession versions. Use a versioned broad gut panel plus appropriate taxonomic decoys. |
| InterPro/Pfam | [InterPro](https://www.ebi.ac.uk/interpro/) | Domain architecture support and homology disambiguation, not automatic clinical phenotype. |
| dbCAN | [resource](https://bcb.unl.edu/dbCAN2/) | Carbohydrate-active enzyme context; broad polysaccharide synthesis/degradation functions do not identify a good or bad biofilm by themselves. |
| VFDB | [resource](http://www.mgc.ac.cn/VFs/) | Reuse existing virulence calls, retaining organism/site context and database evidence. Do not count the same virulence locus twice. |
| AMRFinderPlus / CARD | [AMRFinderPlus](https://github.com/ncbi/amr), [CARD](https://card.mcmaster.ca/) | Keep genetic AMR distinct from biofilm tolerance. Reuse compatible existing calls. Respect resource licensing. |
| Candida Genome Database / FungiDB | [CGD](http://www.candidagenome.org/), [FungiDB](https://fungidb.org/) | Fungal adhesion/matrix context with species/orthology checks. Common fungal regulatory genes do not diagnose candidiasis or an active fungal biofilm. |
| ENA / SRA | [ENA browser](https://www.ebi.ac.uk/ena/browser/home), [SRA](https://www.ncbi.nlm.nih.gov/sra) | Public study reads and explicit assay metadata. Accessions are study inputs, not prevalidated clinical panels. |
| BAP study supplements | [2024 paper and supplementary tables](https://doi.org/10.1038/s41467-024-48309-x) | Table S2 provides the authors' BAP reference identifiers; retain full protein/domain identities and distinguish gene survey from biochemical validation. Do not replace it with every annotation containing “bap.” |
| Existing OpenBiota reference cohorts | Local existing pipeline | Reuse only if assay, taxonomic namespace, extraction/read characteristics and processing are compatible; otherwise uniformly reprocess eligible raw cohorts. |

Use all these as complementary annotation sources, not votes that make a protein harmful merely because it appears in several databases. Database membership is not experimental validation.

### Public research datasets and correct uses

| ID | Public input | Assay / intended use / exclusion |
|---|---|---|
| BF-D01 | [PRJNA644520](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA644520), BF-S01 | Stool/biopsy 16S and endoscopic phenotype context. **Not shotgun training reads.** Stool collection included the first stool after bowel preparation began; preserve that covariate. A reported >3% bloom rule concerned biopsy analyses, not a universal stool threshold. |
| BF-D02 | [PRJEB66343](https://www.ebi.ac.uk/ena/browser/view/PRJEB66343), BF-S06 | Original BAP study stool metagenomes; 49-person analysis described 35 IBS and 14 controls. Verify run/sample mapping before use. |
| BF-D03 | [PRJNA834801](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA834801), BF-S06 reanalysis | Parkinson's cohort, 490 PD and 234 controls in the cited analysis. Disease labels are not biofilm imaging labels. |
| BF-D04 | [PRJNA1279062](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279062), BF-S13 | Paper reports 16S, metagenomic and metatranscriptomic data. Separate specimen, treatment and assay before importing. Ex-vivo treated-community RNA is not a drop-in stool-DNA classifier. |
| BF-D05 | [PRJNA1048869](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1048869), [PRJNA473315](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473315), [PRJNA473316](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473316), [PRJNA857654](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA857654), BF-S09 | Isolate/reference and outbreak contexts. Useful to evaluate annotation/strain specificity; these are not a healthy-stool reference or a validated gut-biofilm outcome dataset. |
| BF-D06 | [PRJNA1228183](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1228183), [2025 CRC tissue study](https://doi.org/10.1038/s41698-025-00873-1) | Tissue 16S context, not whole-gut shotgun truth. |
| BF-D07 | BF-S05 paper and supplements | Retrieve exact gene-family mappings and available abundance tables. A reusable public FASTQ accession was not confirmed in this audit: mark `raw_access_status=unconfirmed`, do not invent one or block the independent gene module. |
| BF-D08 | [Baumgartner image-analysis code](https://github.com/MaximilianBaumgartner/U_Net_bacteria_detection) | Bacterial image segmentation; not a FASTQ biofilm predictor. Optional imaging research adapter only after validation on compatible images. |

### Additional phenotype-linked data and acquisition requirements

| ID | Input | Correct use and availability |
|---|---|---|
| BF-D09 | [PRJNA1076354](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1076354), BF-S19 | Public 16S data from experimental biofilms and mouse work. Audit runs against the article and supplements; do not assume all original human fecal sequences or clinical covariates are publicly deposited. Additional INSPIRE clinical metadata require a governed request. Supports community/phenotype research, not shotgun molecular calibration. |
| BF-D10 | BF-S20 [paper](https://doi.org/10.1093/ecco-jcc/jjad092) and [anonymized supplementary data](https://www.medrxiv.org/content/10.1101/2022.09.09.22279675v1.supplementary-material) | Paired tissue architecture/tissue metagenomes and subset stool qPCR. Raw reads with human sequence are not public; processed sequence data are available on request. Use openly available tables immediately; keep any unavailable sample mapping explicit. No claim of an open stool-WGS training cohort. |
| BF-D11 | BF-S21 [paper and supplements](https://doi.org/10.1080/19490976.2021.1993584) | Cirrhosis stool functions/outcomes. Paper restricts VA metadata and states intent to release deidentified sequences; no specific raw accession was verified here. Source-specific model/feature replication may proceed from available tables; patient-level outcome validation cannot be invented. |
| BF-D12 | [PRJNA1188810](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1188810), BF-S22 | RNA from experimental mixed/single-species communities. Preserve assay, species and experimental condition. Does not directly distinguish biofilm from planktonic state or validate a protective human-stool score. |
| BF-D13 | [Zenodo 16886812](https://doi.org/10.5281/zenodo.16886812), [PV739486](https://www.ncbi.nlm.nih.gov/nuccore/PV739486), BF-S25 | Article links figure/analysis data and a 16S sequence. Validate actual files/version/license. Do not use the 16S accession to claim strain-resolved protective-gene detection or to fill a missing whole genome. |
| BF-D14 | BF-S32 [paper and image/supplementary data](https://doi.org/10.1007/s00253-025-13537-8) | Tissue imaging/context validation; not a public stool-WGS calibration set. Paired samples remain grouped by patient. |

Every dataset gets `host_species`, `specimen`, `assay`, `participant_id`, `sample_id`, `collection_timing`, `phenotype_method`, `phenotype_definition`, `body_site`, `disease_context`, `age`, `medications`, `access_status`, and `label_mapping_status`. An accession alone does not prove compatible labels. Source data without enough mapping can inform a registry but cannot silently enter supervised training.

Add acquisition commands that resolve/download public inputs, inventory their contents, and report missing dependencies individually. Do not make a request-only cohort a prerequisite for the functioning DNA report. Reference construction uses independent assay-compatible controls, whose “healthy” status must not be conflated with microscopy-confirmed biofilm negativity.

No healthy-reference percentile may be calibrated on the five local samples, on SAMPLE2 as the single disease exemplar, or on the two favored donors. That would make the score circular.

## 6. Required biological module census

Each module must be represented in the registry and report as `implemented`, `candidate_not_validated`, `context_only`, `source_unavailable`, or `not_measurable_from_current_assay`. This prevents silent omissions without pretending every candidate is score-ready.

| Module ID | Detectable/context features | Interpretation and scoring policy |
|---|---|---|
| BF-M01 | Curli-family structural/assembly gene evidence with carrier and locus context | Bacterial amyloid potential. A complete supported locus is stronger than an isolated regulator. Healthy carriage also occurs; harmful-axis eligibility requires an explicit experimentally supported adverse-context claim, not merely csgA detection. Never equate with host fibrin or brain plaques. |
| BF-M02 | BAP-family targets from BF-S06 Table S2; exact target/domain annotations | Amyloid-associated adhesive potential. Keep validated domain, untested homolog, organism and study model separate. Domain hits in repetitive surface proteins are especially vulnerable to ambiguity. |
| BF-M03 | PNAG/PGA synthesis/export locus, including the pga gene family; organism-specific ica systems where relevant | Matrix-production potential, initially context-dependent. Do not force universal harm or claim all systems are interchangeable. |
| BF-M04 | BAD study gene pattern: bssS and pgaA/B/C/D; separate R. gnavus/E. coli/bile-acid context | Exploratory human disease-associated pattern. One grouped signal, not six independent proofs of harm. Community correlation is not carrier linkage. |
| BF-M05 | Klebsiella type-3/type-1 adhesion and other supported adhesin loci with species/strain attribution | Persistence/adhesion potential with clinical-isolate evidence. Not active gut biofilm or universal virulence. BF-S08 anchors include MrkD WP_004149659.1, FimH BAH65076.1, EcpD WP_002890060.1; these are annotation anchors, not sufficient single-reference diagnostic panels. |
| BF-M06 | Organism-specific extracellular polysaccharide systems and matrix-associated proteins | Broad discovery inventory. Capsule, cellulose and other matrix capabilities must remain separately typed; a capsule hit alone is not attached biofilm. |
| BF-M07 | DNABII/eDNA-related matrix associations, quorum and cyclic-di-GMP regulatory context | Many genes are broadly conserved/pleiotropic. Report annotation context; no harmful-score inflation from housekeeping counts. Total extracted DNA cannot distinguish intracellular from matrix DNA. |
| BF-M08 | C. difficile adhesion/matrix/spore context from existing validated organism and toxin assays | Separate biofilm potential, viable spores and toxin potential. Stool DNA cannot measure spore viability or prove a biofilm reservoir. |
| BF-M09 | Candida species-specific adhesion, matrix and regulatory context | Research fungal biofilm potential; presence alone is neither candidiasis nor invasive disease. Low eukaryotic depth and paralogs must be explicit. |
| BF-M10 | L. reuteri ATCC 23272-associated GtfW mechanism from BF-S15 | Candidate protective-associated mechanism, conditioned on identity/context. The published microsphere formulation cannot be inferred from stool DNA. Species-level abundance alone does not activate the score. |
| BF-M11 | F. prausnitzii/host-glycan mechanism from BF-S14 | Protective evidence card and optional multimodal context. No validated specific DNA panel established here; generic F. prausnitzii/cobalamin abundance is not substituted. |
| BF-M12 | Fiber-associated mixed-community mechanisms from BF-S16 | Supportive ecological evidence; no exclusive DNA signature established. Do not count fiber-degradation genes as proof of protective biofilm. |
| BF-M13 | Other verified probiotic-strain adhesion/EPS mechanisms | Expand only with exact strain/gene/phenotype evidence and a biological claim audit. Generic “Lactobacillus equals good” rules are forbidden. |
| BF-M14 | S. aureus host fibrin-interaction capability, with organism/strain linkage | Extraintestinal/context panel, not automatic harmful-gut-axis input. Fibrin is a host substrate, not detected by these genes. |
| BF-M15 | Co-occurring existing pathogen, toxin, AMR or damaging-function calls | Cross-reference as independent clinical/mechanistic context. Co-occurrence is not same-cell linkage or a measured resistant biofilm. |
| BF-M16 | Measured host thrombin/protein activity and tissue architecture | Optional assay layer. Never infer from F2 DNA, taxon abundance or generic fliC expression. |
| BF-M17 | Longitudinal donor/recipient shared strain plus supported biofilm-associated loci | Track potential and strain persistence. No conclusion that an attached biofilm or a disease was transferred. |

### Additional scoring, behavior and validation modules

| Module ID | Measurement | Required implementation |
|---|---|---|
| BF-M18 | Two-feature human biofilm-associated community proxy | BF-S01 tissue-derived E. coli/Escherichia-Shigella and R. gnavus association translated to a deliberately narrower stool-WGS hypothesis in §8.7. Explicit transport label, no stool biofilm probability. |
| BF-M19 | Four-genus biofilm-negative-associated ecological proxy | Faecalibacterium, Coprococcus, Subdoligranulum and Blautia, from BF-S01 biopsy associations. §8.7 specifies a separately named proxy; these do not become exclusive beneficial-biofilm organisms. |
| BF-M20 | Linked matrix/adhesion plus independently supported damaging capability | Join existing molecular calls only with supported carrier linkage. List the independent damage claim and source; co-occurrence on different unknown carriers cannot satisfy the conjunction. It is a concern-context overlay, not an added copy of the same feature in the axis. |
| BF-M21 | Feces-derived community behavior | Optional externally measured biomass, dispersal, viable burden, epithelial/barrier effect and intervention response, anchored to BF-S19/S23. Preserve lab-condition scope; no stool-DNA reconstruction of these values. |
| BF-M22 | Tissue architecture and host-interface observations | Optional image/pathology-derived coverage, thickness, epithelial proximity and mucus/host response. Use validated measurement adapters and per-site sampling, without creating a universal whole-gut amount. |
| BF-M23 | Matched expression/activity support | Optional RNA/protein support for a named DNA mechanism. DNA and RNA must match sampling/context before a joint claim; generic expression is not physical attachment. |
| BF-M24 | Additional source-specific protective-antagonism evidence | BF-S25–S31 expand the protective evidence registry. Implement exact strain resolution where public reference support permits; otherwise retain species/family context and action cards. No assumption that a named probiotic species implies the tested strain/formulation. |

BF-M18/19 are intentionally experimental proxies, analytically computable from compatible taxonomy. Their inclusion in a separate proxy display is not a relaxation of BF-M10/11/13's requirements for **direct protective-mechanism** scoring. BF-M20–24 add observations and interpretation; they are not five automatic extra weights.

### Module manifest contract

```yaml
schema_version: openbiota.biofilm-module/1.1
module_id: BF-M04
label: BAD-associated biofilm-gene pattern
measurement_kind: community_gene_pattern
mechanism_group: pnag_and_regulation
primary_source_ids: [BF-S05]
feature_symbols: [bssS, pgaA, pgaB, pgaC, pgaD]
sequence_mapping_status: requires_source_identifier_resolution
analytical_status: pending_validation
biological_evidence_status: human_association
clinical_validation_status: not_established
molecular_resolution: gene_family
carrier_required_for_organism_claim: true
interpretation:
  default_valence: context_dependent
  evidence_relation: human_cross_sectional_association
  phenotype_directly_measured: false
  disease_context: bile_acid_diarrhea
  clinical_predictive_validation: not_established
score:
  candidate_axis: harmful_associated
  status: candidate_not_validated
  aggregation_group: BF-G-PNAG
  independent_of: []
  activation_requires:
    - reviewed_sequence_manifest
    - analytical_validation
    - compatible_reference_transform
    - explicit_association_only_label
report_even_if_unscored: true
```

`feature_symbols` are semantic inventory, **not a sequence classifier**. Resolve public accession versions, taxonomic scope, gene-family mapping, length/domain architecture and evidence before the assay becomes active. An unavailable mapping results in a visible target-specific status, not an invented match or a fabricated accession. Preserve source coordinates and reference sequences in the application's pinned registry build.

## 7. Detection and evidence pipeline

### 7.1 Resolve inputs and caches

Prefer original host-filtered paired reads and compatible gene/assembly caches. Record sample provenance, extraction protocol, read length, library type, depth, duplicate handling and host filtering. Keep paired fragments paired. Do not upload patient reads to public services.

A taxonomic profile alone may support context cards but cannot reconstruct missing accessory loci. A HUMAnN gene-family table may support family abundance if its identifiers and units are compatible; it does not automatically establish strain, locus integrity, carrier or operon completeness. A MetaPhlAn marker alignment may support existing strain work but is not a substitute for a whole-read biofilm gene screen.

Cache key: input hashes + preprocessing version + assay configuration + reference manifest hash + tool/container versions + registry release + scoring/calibration versions. Source-evidence wording changes invalidate interpretation/report caches; changed sequences invalidate detection caches. Do not re-run all taxonomic profiling for a prose-only change.

### 7.2 Molecular identification

Use existing validated nucleotide/protein search and competitive mapping adapters. Broad candidate discovery and confirmatory identification are separate. Include close homologs and taxonomic decoys so a shared domain does not force a species-specific attribution. Preserve ambiguous evidence instead of assigning it to the best-known pathogen.

For every call store independent fragment support, target/reference coverage, breadth at the required depth, uniquely informative positions, identity/alignment quality, read-pair support, alternative matches, low-complexity/repeat fraction, and actual assay detection limits. Protein homology supports family identity; exact strain/allele claims need corresponding discriminatory evidence.

No universal identity, read-count, or breadth threshold is scientifically appropriate for all these loci. For each target family, establish and freeze analytical criteria using positives, close negatives, mixtures and realistic read simulations from public reference sequences. Until that target is validated, label it `candidate_match`; it can appear in research details but cannot silently become a definitive allele/strain or clinical claim. Never copy the BAP publication's short-match/presence rule as a clinically validated threshold.

### 7.3 Locus completeness and carrier linkage

An operon is not proven intact merely because its component genes appear somewhere in a mixed community. Use contig/assembly linkage, read-pair constraints, supported taxon-stratified calls and strain evidence, each with its own confidence. Mark weak co-abundance-based inference as inference. Short-read assembly can collapse strains, separate plasmids from hosts or miss repeats; retain alternatives.

`linkage_state` values: `direct_supported`, `probabilistic`, `taxon_stratified_only`, `unassigned`, `conflicting`. Only the first state can satisfy a rule requiring direct same-carrier linkage; probabilistic rules must name their calibration and error rate. AMR and biofilm genes on separate unlinked contigs cannot be declared part of one resistant strain.

Distinguish:

- `detected_supported`: validated molecular evidence for the stated target/resolution.
- `candidate_match`: suggestive but not validated at requested resolution.
- `not_detected_above_assay_limit`: adequate assay, no detection; not biological absence.
- `below_resolution`: species/family signal but insufficient finer resolution.
- `not_assayed`: no applicable assay run.
- `unresolved`: contradictory/ambiguous evidence.
- `failed_qc`: data unsuitable for the requested measurement.

A target absent from the database is `not_assayed`, not negative. Nondetection is left-censored, not proof of a true numeric zero. Store `value`, `censoring`, `detection_limit` and `quantification_limit` separately. A frozen registry may encode adequate-assay nondetection as zero for ranking, but must apply exactly the same encoding to references and retain the censoring flag. Unresolved, failed-QC and unassayed measurements stay null. A low-depth organism does not get a “biofilm genes absent” badge.

### 7.4 Quantification

Keep units explicit: gene-family CPM, length-normalized fragment abundance, estimated copies per microbial genome equivalent, carrier-normalized copies, or taxonomic relative abundance are different measurements. Never blend them without a documented transformation.

Preferred new locus quantification is length-normalized, confidently attributed fragment support, with both per-sample microbial-library normalization and per-carrier normalization when reliable carrier genome coverage exists. Report each separately. Use validated single-copy microbial reference markers for genome-equivalent normalization; exclude unreliable denominators. No invented absolute cells per gram without an appropriate external measurement/spike-in.

For a supported multi-gene mechanism, retain all components and compute module abundance only under the module's predeclared rule. An essential-component minimum may be used as an engineering abundance estimate only when that locus biology and same-carrier linkage justify it; otherwise keep a vector or community-pattern score. A regulatory gene alone cannot stand in for the physical matrix.

### 7.5 Optional RNA, proteins and imaging

Import them as separate observations with specimen/date/site, units, assay validation and matched-sample status. RNA abundance must not be renamed DNA abundance, and DNA-normalized transcription requires matched validated measurements. Amyloid dye binding is not uniquely bacterial, not uniquely curli, and not a quantitative fibrin test. Preserve antibody/protein specificity and orthogonal controls where reported.

Clinical endoscopy is not automatically ordered from a research score. If an independently indicated procedure provides observations, accept structured findings with site, bowel preparation, fixation, staining, microscopy and sampling limitations. Standard specimen preparation can disturb mucus/architecture; absence in one sample does not establish whole-gut absence.

## 8. Reproducible scoring specification

### 8.1 Separate molecular observations from engineered ranking

There are no literature-derived universal good/bad coefficients to copy. The following algorithm is a **versioned engineering summary of registered evidence**, not a validated clinical diagnostic model. All weights are transparent, and clinical validation is a separate status.

Each scoring module needs a fixed feature transform, evidence context, eligible population, fixed dependence group and active registry release. Assign the biological direction from the reviewed claim, not from a model's accidental case/control correlation. Generic adhesion, capsule, quorum, BssS regulation, or a beneficial species name alone cannot establish directional health value.

Keep two forms of output:

1. **Measured module view:** feature abundance, call confidence, mechanism context and optional module percentile. Always available to the degree the data permit.
2. **Axis reference summary:** computed only for a frozen supported module set and a compatible population. It is explicitly “research reference percentile.”

### 8.2 Reference construction

Use independent participants with compatible stool-shotgun processing and appropriate metadata. Freeze discovery, reference and validation cohorts separately. Maintain age, medication, recent antibiotics, bowel preparation, geography, stool consistency, extraction batch and relevant disease metadata. Subject-level splits prevent repeated samples crossing folds.

For the existing SAMPLE4/SAMPLE6 minors, never silently describe adult-reference tails as pediatric abnormality. If an age-compatible reference is unavailable, show molecular findings and optionally a separately labeled adult research comparison; do not use that comparison to claim healthy, unhealthy or donor clearance. Do not infer chronological age from file names.

Default engineering requirements for a displayed reference percentile: at least 30 independent compatible reference participants and nonconstant reference values. This is a minimum for an exploratory display, not a clinical adequacy claim. Display n and uncertainty. Tail interpretation below the 5th or above the 95th percentile should remain descriptive; do not make tail-specific claims when reference sample size is too small to resolve them. Reference-health questionnaires alone do not prove absence of mucosal biofilm.

### 8.3 Exact reference-rank algorithm

For module m, define `x_m` as its frozen nonnegative quantitative measurement after the exact registry transform. Transform definitions and denominator must match the reference. A log transformation can be stored for modeling, but the empirical rank does not require choosing a pseudocount.

For reference observations `r_1 ... r_n`, define midrank percentile:

`P_m(x) = 100 * (count(r < x) + 0.5 * count(r == x)) / n`.

Use stored numerical precision/tie rules. If the entire reference is constant, a module percentile is not informative: return `reference_non_discriminating`, keep the raw measurement and detection state. In a zero-inflated but nonconstant reference, a detected absence may have a nonzero midrank. Always display “not detected above this assay's limit” beside it; a percentile is not a detection probability. Never relabel a tied zero as measurable positive potential.

Within each prespecified biological dependence group, average the eligible module percentiles with equal weights, normalized by the fixed group's module count. Across groups, average with equal group weights. This prevents a long operon, many homologs, or many papers on one mechanism from dominating the axis. Assign weights before viewing SAMPLE2/SAMPLE4/SAMPLE6. No high-confidence mapping multiplier that turns molecular certainty into biological harm.

For v08 descriptive reference summaries, compute every reference participant's aggregate with the **same frozen full-reference module transforms** used for a new sample, then store the aggregate distribution. The reported axis is the empirical midrank percentile of the sample aggregate against that distribution. Identical observations must remain tied through both stages. Do not mix leave-one-out reference transforms with full-reference sample transforms. This is a descriptive empirical ranking, not out-of-sample predictive validation. A future independent fit/calibration design may use separate cohorts, but calibration and new samples must still use the identical fit-reference transform. Fit and version transformations once; do not refit on the recipient/donor batch. Freeze a separate calibration object for each supported missingness mask if a partial axis is offered. Never silently drop missing modules and rescale against the full-panel reference.

Do not include BF-M03 and BF-M04 as independent evidence for the same PNAG signal. Shared features get one dependence group. Amyloid genes and their annotated matrix function likewise cannot count repeatedly just because several databases report them.

The exact aggregation is `G_g = sum(P_m for m in fixed_group_g) / len(fixed_group_g)` and `A = sum(G_g for g in fixed_axis_groups) / len(fixed_axis_groups)`. Every input in the frozen calibration mask must be eligible. Unresolved, constant-reference, incompatible-unit, inactive or failed inputs cannot contribute zero. A partial mask must retain at least one eligible biological group, and a one-group result must retain a limited-scope label. Masks are selected only by prespecified assay availability/QC rules, never by removing detected, unfavorable or inconvenient findings. Record a reason for every omitted module.

### 8.4 Availability, research activation and sparse panels

Molecular assay validity and clinical validity are separate. A supported gene-family assay plus a source-specific biological claim can be activated for **experimental reference ranking** without a human outcome trial. The registry must distinguish `analytical_status`, `biological_evidence_status` and `clinical_validation_status`; one generic “validated” flag is insufficient.

A single active group can receive a percentile under its exact named scope. Two or more groups can receive a multi-mechanism summary, still limited to those groups. Remove the old `insufficient_active_groups` stop at two: the correct stop is zero usable groups. The display shows `scope=single_mechanism` or `scope=multiple_mechanisms`, their names, and coverage. A single GtfW-related comparison cannot be titled “overall protective ecosystem health.”

Preset masks are selected by available assays/QC, including valid nondetections; **not by whether a favorable or unfavorable target was detected**. A nondetected but adequately assayed target stays in the mask. Use a separately calibrated mask when necessary. Compare samples only within the same mask, reference and measurement scope.

All candidate discoveries remain visible, but ambiguous sequence identity does not become a confident molecular call because the feature is interesting. No clinical-validation gate may hide otherwise supported research scores. Conversely, no proxy may silently replace a missing DNA-mechanism result. The four-card layout in §4 allows the complementary evidence to remain useful.

### 8.5 Uncertainty

Report uncertainty sources separately:

- Analytical uncertainty: fragment ambiguity, depth, locus resolution and sequencing resampling.
- Reference uncertainty: bootstrap independent reference participants, not repeated samples as independent people.
- Mechanistic transport: non-gut model, species/strain mismatch, unmeasured expression, formulation/site mismatch.
- Clinical uncertainty: no validated relationship between the index and benefit, harm or treatment response.

Use 200 bootstrap replicates as an initial engineering default, with a convergence check on representative data. Report separate 95% analytical and reference intervals using 2.5th/97.5th percentiles with a documented linear quantile method. Derive a deterministic seed from SHA-256 of sample input hashes, registry release, calibration hash, interval kind and a versioned seed string. Analytical resampling uses independent read fragments, keeps pairs together, and reruns the supported quantitative transform while holding the calibration fixed; if only a summary table is available, analytical resampling is unavailable rather than invented. Reference resampling uses independent participants, retains their paired module vectors, rebuilds both ranking stages in each replicate under the identical algorithm, and holds the sample measurement fixed. Do not pool these intervals into one unless a separately specified joint bootstrap is implemented. Label them “sampling interval” and “reference interval,” not “95% certainty about biofilm health.” A bootstrap replicate with constant module/aggregate references, failed analytical evidence or an incompatible frozen panel is invalid; do not change the panel, substitute zero, or drop that module to rescue it. Record attempted/valid counts and failure reasons. With fewer than 180 valid replicates out of the initial 200, return a null interval with `bootstrap_unstable`; otherwise label the interval as conditional on valid replicates and disclose the invalid fraction. The 90% rule is an engineering stability threshold, not a coverage guarantee. A constant original axis aggregate likewise returns `reference_non_discriminating`. Resampling cannot capture unknown biology. Do not synthesize a clinical confidence interval from a mapping quality score.

### 8.6 Test vectors

- Reference [0, 0, 2, 4], x=0: midrank25; raw result remains nondetection if the assay did not detect the target. x=2:62.5. x=5:100.
- Reference [0, 0, 0, 0]: no informative module percentile; no fabricated 50/100 biological potential.
- Protective 90, harmful 90: display both; never compute zero concern by cancellation.
- Full panel A/B/C, C unresolved: no full-panel axis. Use a separately frozen A/B calibration if available, otherwise module results only.
- Same recipient processed alone or beside ten donors: identical calls and scores when inputs/versions are unchanged.
- One gene represented in three databases and four papers: one molecular observation with several citations, not twelve evidence units.

### 8.7 Explicit starter proxy indices from existing profiles

These are **new engineering hypotheses**, not published predictive models. BF-S01 supplies the associations; this specification supplies the following fixed, transparent arithmetic. They are deliberately separate from H-M/P-M. Validate whether they predict anything beyond general dysbiosis before promoting their interpretation.

**H-C: BF-PROXY-COMMUNITY-1.0**

- Feature A: E. coli species-complex relative abundance from the pinned WGS profiler. BF-S01's 16S Escherichia/Shigella association cannot resolve this identically; record `transport_change=genus_complex_to_species_complex`. Do not assign every Escherichia/Shigella read to pathogenic E. coli, and do not infer a pathotype.
- Feature B: Ruminococcus/Mediterraneibacter gnavus relative abundance, with a fixed synonym/taxonomy mapping and no duplicated counts.
- For each feature obtain its ordinary abundance percentile P against the same compatible WGS reference. Compute `raw_index = (P_A + P_B)/2`; then re-rank this index against reference individuals transformed by the same frozen reference, using §8.3. Higher means more of this **two-feature association pattern**, not more measured harmful biofilm.
- Both features must be assayed. If one is unavailable, show the other as a named feature percentile; do not retain the two-feature proxy title. Nondetection above a qualified assay limit is permitted and retains its censoring flag.

**P-E: BF-PROXY-ECOLOGY-1.0**

- Four inputs: genus-level relative abundances of **Faecalibacterium, Coprococcus, Subdoligranulum, Blautia**. BF-S01 reported these as relatively depleted in biofilm-positive **biopsies**. Their stool-WGS use is a hypothesis, not validated stool transport or proof they build protective biofilms.
- Aggregate each genus from the pinned taxonomy without overlapping leaves, preserve unclassified assignments, and never combine incompatible MetaPhlAn/GTDB taxonomies without an audited mapping. Current-name changes are versioned.
- Compute the four abundance percentiles against the same WGS reference, take their equal mean, and re-rank that mean as in §8.3. Treat the entire set as **one correlated ecological group**; it is not four independent protective mechanisms. All four are required for the named proxy.
- Low values indicate relatively less of this narrow community pattern. They do not diagnose a missing beneficial biofilm, species deficiency or need for FMT. The proxy does not instruct supplementation with these organisms.

Do not add Shannon diversity, a Firmicutes/Bacteroidetes ratio, generic butyrate genes or all Lactobacillus species to these fixed indices. Such measurements remain useful separately in the existing report. Changes to feature sets require a new index version and calibration. These proxies can be tested immediately against matched phenotype datasets at the **assay level actually available**: a 16S reproduction is not validation of a WGS deployment.

### 8.8 Ranked findings and action priorities

Produce two ranked lists, **findings needing review** and **protective/supportive observations**. Ranking is an evidence-navigation order, not a disease-risk probability. Use a deterministic lexicographic key:

1. Relevant measured adverse phenotype or independently confirmed clinical concern, when supplied; then supported same-carrier damaging-function context; then source-qualified molecular association; then community proxy; then ambiguous candidate/context.
2. Applicable human-gut evidence; then relevant human-derived community evidence; then animal-gut evidence; then other in-vitro/site evidence; then in-silico hypothesis.
3. Within an otherwise identical class and compatible assay/panel, descending measured concerning-feature rank, or descending supportive-feature rank on the supportive list. Missing percentiles sort after available ones, but cannot suppress a confirmed clinical concern.
4. Stable module/call ID for deterministic ties.

Do not describe this ordering as a clinically validated severity score. Do not rank incompatible percentile masks as though a 90 in one measured scope means more biological harm than an 80 in another. Display why the finding is placed where it is, the source, and what would resolve the most important uncertainty.

#### Required explanation of every high or low result

Beside every score expose the exact contributing feature names, raw values/units, individual reference percentiles, mean before final reranking, and evidence context. Give the single largest driver and whether its abundance or carrier-normalized capability changed. Where the aggregate is dominated by one common taxon, say so plainly. A percentile of 84 means the selected measurement ranks around that position in the named reference; it does not mean an 84% harmful biofilm, an 84% disease match, or an 84% chance of transferring disease.

Show an **agreement table**, not a fabricated combined confidence number: molecular mechanism; community proxy; expression/activity; directly measured architecture; measured host consequence. Each has `supports_named_concern`, `opposes_named_concern`, `mixed`, `not_measured`, or `not_applicable`, with the named concern and source-specific rule. Generic inflammation cannot fill the architecture row. Two scores computed from the same reads or taxa are not independent confirmations. Contradictory direct observations remain visible; no majority vote turns three correlated proxies into proof.

### 8.9 Behavior-based validation and optional scores

BF-M21 supports a research laboratory's returned quantitative results. Store the measured vector: matrix biomass, viable attached burden, viable dispersed burden or dispersal rate, epithelial injury/permeability/inflammatory readouts, and beneficial-community function where actually tested. Some assays may be unavailable. Do not infer the entire vector from one stain or cytokine.

For every assay with a compatible external laboratory reference, a named phenotype percentile can use §8.3. Without such a reference, retain the raw value and, if supplied, its actual comparator effect. Separate matrix amount from **harm-related behavior** and **protective function**; their scales need not share physical units. A laboratory phenotype remains “feces-derived community, ex vivo,” not the patient's directly measured intestinal state. Optional assays must not prevent the DNA report or experimental proxies from completing.

The software ingests qualified laboratory results; this spec does not prescribe culturing, manipulating or exposing people to organisms. Do not recommend an invasive investigation solely to fill this panel. Use observations from independently indicated clinical assessment or appropriate research arrangements.

### 8.10 Longitudinal interpretation and criterion for improvement

Freeze the assay, proxy definitions, references and specimen-handling metadata across time. If versions change, reprocess all compared time points under the same version or show a break in the series. Report change in raw gene signal, carrier abundance, carrier-normalized signal and each independent index. A higher rank for one taxon can result from another taxon declining; it is not automatically absolute expansion.

Use labels that match endpoints: “lower measured gene signal,” “less of the study-associated community pattern,” “lower ex-vivo epithelial injury,” or “lower viable pathogen burden” only when that respective outcome was measured. Reserve “biofilm eliminated” for evidence that genuinely supports eradication at a specified sampled site; DNA nondetection alone cannot meet it. A temporary symptom flare is not evidence of successful biofilm killing. Biomass reduction accompanied by increased viable dispersal or host injury is an unfavorable/mixed result, not an automatic success.

### 8.11 Reference implementation for deterministic rank arithmetic

The following standard-library Python is normative for the arithmetic only. The surrounding application must enforce source, assay, age, scope, frozen-mask and minimum-reference-size rules before calling it. Small arrays here are mathematical fixtures, not real calibration cohorts.

```python
from math import isfinite
from statistics import mean


def midrank(value, reference):
    if not reference or not isfinite(value):
        raise ValueError("missing or nonfinite input")
    if any(not isfinite(x) for x in reference):
        raise ValueError("nonfinite reference")
    if min(reference) == max(reference):
        raise ValueError("reference_non_discriminating")
    return 100.0 * (sum(x < value for x in reference)
                    + 0.5 * sum(x == value for x in reference)) / len(reference)


def rank_panel(sample, reference_rows, groups):
    # Inputs have already been rounded to the calibration's stored precision.
    names = [name for group in groups for name in group]
    if not groups or any(not group for group in groups):
        raise ValueError("empty panel/group")
    if len(names) != len(set(names)):
        raise ValueError("duplicate feature in dependence groups")
    if set(sample) != set(names):
        raise ValueError("sample mask differs from calibration")
    if not reference_rows or any(set(row) != set(names) for row in reference_rows):
        raise ValueError("reference mask differs from calibration")
    values = [sample, *reference_rows]
    if any(value is None or not isfinite(value)
           for row in values for value in row.values()):
        raise ValueError("missing/nonfinite value")
    columns = {name: [row[name] for row in reference_rows] for name in names}

    def aggregate(row):
        return mean(mean(midrank(row[name], columns[name]) for name in group)
                    for group in groups)

    aggregate_reference = [aggregate(row) for row in reference_rows]
    raw = aggregate(sample)
    return {"raw_index": raw,
            "reference_percentile": midrank(raw, aggregate_reference),
            "feature_percentiles": {n: midrank(sample[n], columns[n]) for n in names}}
```

`groups=[["ecoli_complex", "m_gnavus"]]` defines the H-C arithmetic; `groups=[["faecalibacterium", "coprococcus", "subdoligranulum", "blautia"]]` defines P-E. H-M/P-M use their frozen dependence groups. The grouped means and the final reranking are identical operations across these cases, while labels and evidence differ. Constant references raise an explicit state; they are not repaired by adding pseudodata. One feature/group is supported. Preserve equal input/output ordering in the actual serialization.

## 9. Unified result contracts

```json
{
  "schema_version": "openbiota.biofilm/1.1",
  "sample_id": "EXAMPLE_NOT_A_MEASUREMENT",
  "spec_version": "08.0.1",
  "input_assays": ["stool_shotgun_dna"],
  "registry_release": "resolved_at_build",
  "summary": {
    "measurement_label": "Biofilm-related genetic potential",
    "protective_associated": {
      "status": "not_analyzed",
      "reference_percentile": null,
      "intervals": {"analytical": null, "reference": null},
      "calibration": null,
      "omissions": [],
      "module_ids": [],
      "reason_codes": []
    },
    "harmful_associated": {
      "status": "not_analyzed",
      "reference_percentile": null,
      "intervals": {"analytical": null, "reference": null},
      "calibration": null,
      "omissions": [],
      "module_ids": [],
      "reason_codes": []
    },
    "community_pattern_proxy": {"status": "not_analyzed", "reference_percentile": null},
    "protective_ecology_proxy": {"status": "not_analyzed", "reference_percentile": null},
    "active_biofilm_burden": null,
    "anatomical_location": null,
    "clinical_disease_probability": null
  },
  "module_results": [],
  "context_results": [],
  "external_assay_results": [],
  "intervention_evidence_ids": [],
  "provenance": {
    "input_sha256": [],
    "reference_manifest_sha256": null,
    "tool_versions": {},
    "parameters": {},
    "code_commit": null,
    "calibration_id": null
  }
}
```

Each `module_result` must include: stable module/call IDs, requested and achieved resolution, quantitative value/unit, component gene calls, coverage, linkage, biological-context claim ID, source IDs, active/candidate status, score contribution/dependence group, alternatives, uncertainties, and explanatory text. Reuse `openbiota.resolution/1.0` evidence IDs when present.

Each `external_assay_result` needs: specimen, body site, collection date, laboratory, method, analyte, unit, value, assay limit, reference population, interpretation scope and linked sample. A report renderer must never create a value that is absent from the structured result.

### 9.1 Normative type definitions and invariants

The coding agent must implement these as actual schema validation (the repository's existing schema library is preferred), not just prose comments. Unknown object keys should be rejected in core measurement records unless placed in an explicitly versioned extension field.

| Type / field | Required type and rule |
|---|---|
| `ModuleResult.module_id`, `call_ids` | Stable string; array of referenced immutable call IDs |
| `analytical_status` | `pending_validation`, `supported_at_stated_resolution`, `failed_validation`; applies to molecular identity/quantification, not clinical utility |
| `biological_evidence_status` | `human_association`, `human_derived_ex_vivo`, `animal_mechanism`, `in_vitro_mechanism`, `in_silico_hypothesis`, `mixed`; claim-specific sources remain mandatory |
| `registry_status` | `implemented`, `candidate_not_validated`, `context_only`, `source_unavailable`, `not_measurable_from_current_assay` |
| `measurement_status` | `detected_supported`, `candidate_match`, `not_detected_above_assay_limit`, `below_resolution`, `not_assayed`, `unresolved`, `failed_qc` |
| `measurement` | Object with nullable finite numeric `value`, nonempty `unit`, `censoring` (`none`, `left`, `interval`, `unknown`), nullable finite nonnegative `detection_limit` and `quantification_limit`, and `normalization_id` |
| `components` | Array of component IDs, measurements, validity and linkage evidence; empty if not assayed |
| `requested_resolution`, `achieved_resolution` | Controlled strings from shared resolver; achieved may be null |
| `linkage` | Object with state from §7.3, nullable carrier ID, evidence IDs and alternatives |
| `valence` | `protective_associated`, `harmful_associated`, `context_dependent`; a module may have different source-specific claims, not a silently averaged biological direction |
| `claim_ids`, `source_ids` | Arrays of resolvable registry IDs; no unresolvable ID in a published card |
| `clinical_validation_status` | `not_established`, `research_association`, `externally_validated_for_named_endpoint`; the last requires an endpoint/population/model validation record |
| `module_percentile` | Nullable finite number in [0,100], with its own calibration provenance |
| `AxisResult.status` | `available_research`, `available_limited_scope`, `not_analyzed`, `no_eligible_measurements`, `reference_unavailable`, `reference_non_discriminating`, `out_of_domain`, `missing_required_measurements` |
| `AxisResult.reference_percentile` | Number in [0,100] only with a valid axis calibration; otherwise null. An explicitly requested out-of-domain comparison lives in `research_comparison`, not the primary value |
| `AxisResult.intervals` | Separate nullable `analytical` and `reference` objects: `{level:0.95,lower:number,upper:number,method:string,replicates:integer,attempted:integer,valid:integer,invalid_reasons:object,seed:string}`; bounds ordered and within [0,100] |
| `AxisResult.calibration` | Null for unavailable scores; otherwise an object containing calibration ID/hash, transform version, fixed module IDs, fixed group mapping, mask ID, selection-rule version, participant n and applicability status |
| `AxisResult.omissions` | Array of module ID and controlled reason; not an unexplained denominator reduction |
| `AssayResult` | Specimen ID, sample link, site, collection date, lab/method, analyte, nullable `value` and explicit `unit`, measurement status, assay limits, provenance and separately applicable reference |
| `CalibrationManifest` | Immutable ID/hash; axis/module name; registry/assay versions; feature transforms and units; censoring encoding; group membership; weights; mask rule; participant-level metadata summary; reference data hash; algorithm/version; tie precision; applicability predicate; validation status |

Additional mandatory fields for every available score: `score_kind` (`mechanism_percentile`, `community_proxy_percentile`, `ecological_proxy_percentile`, `measured_phenotype_percentile`), `scope`, `feature_definition_id`, `assayed_feature_ids`, `planned_feature_ids`, `measured_scope_label`, `biological_evidence_status`, `transport_assumptions`, `clinical_validation_status`, and `reference_percentile` with the full §9.1 calibration/interval structure. The short null proxy objects above are permitted only for `not_analyzed`; available proxy objects require the same provenance and coverage fields as mechanism axes.

`RankedFinding` contains `finding_id`, `call_ids`, `rank`, `ranking_rule_version`, explicit `sort_key`, `reason_text`, `evidence_ids`, and `next_information_needed`. `ActionCandidate` adds the underlying intervention/evidence IDs, matched target, intended endpoint, evidence tier, delivery/formulation context, contradictory effects, practical access status and a non-probabilistic relevance rank. No numerical efficacy field is populated without a valid study-specific estimate and context.

Invariants: all cited calls/sources exist; NaN/infinity are rejected; null differs from zero; module IDs in an axis exactly match its calibration mask; group weights sum to one after the prespecified equal-weight construction; all reference comparisons preserve age/domain status in JSON/HTML/PDF. A top-level calibration ID is optional convenience only: every scored module/axis carries its own provenance.

Migration: preserve existing 1.0 measurement calls and provenance. Recalculate interpretation/score objects under 1.1, map the old `insufficient_active_groups` state to an explicit re-evaluation task, and retain the prior result in its original report history. Do not rewrite old output as if it had originally used these new proxy definitions. Candidate-module status alone must not obscure whether the remaining blocker is sequence identity, reference availability or clinical validity.

A partial mask cannot be chosen because an unfavorable or favorable feature was present. Calibration activation and sample availability are independent states. An empty calibrated mechanism panel produces a complete report with null mechanism values, any separately computable proxy values, visible molecular candidates and explicit reasons; the application must not crash or activate an unvalidated candidate merely to populate a graphic.

## 10. Intervention evidence engine

**A human clinical trial is not required to list a directly supported candidate.** Include laboratory, animal, ex-vivo and human evidence; identify each clearly. Relevant allicin, berberine and enzyme evidence must be retrievable even when the clinical effect is unknown.

The rule is: relevant experimental effect + a plausible or investigational delivery route => an evidence-linked candidate card. It is not automatically a proven treatment, a personalized dose, or a prediction that the patient's strain will respond. Unknown oral exposure is displayed as unknown, rather than assumed impossible or successful.

Distinguish prevention of formation, reduction in matrix biomass, dispersal of mature biofilm, killing of viable organisms, regrowth-free eradication, symptom change, and pathogen-load change. Crystal-violet decline alone is not eradication. Planktonic antibacterial activity alone is not antibiofilm activity. A mixture's result cannot be assigned to every ingredient.

### Evidence selection and ordering

Retrieve by exact molecule/formulation, organism or supported mechanism, body site, model and endpoint. Allow species-level matching to an experimentally tested strain, but label the strain transport uncertainty. Never assign E. coli-specific biofilm evidence to R. gnavus without supporting data.

Display four groups: dietary/ecological support; compound/enzyme candidates; organism-specific medical options; investigational precision approaches. Within each group sort by directness to the measured target and site, then evidence maturity and study quality. Do not manufacture efficacy probabilities or combine unlike endpoints into a benefit percentage. Strong contradictory results appear on the same card, not hidden in a bibliography.

### Deterministic retrieval order and meaningful action categories

Create four action sections: **support protective ecology**, **reduce a matched concerning mechanism**, **clarify or treat an independently established condition**, and **investigational precision approaches**. List relevant preclinical candidates even without clinical trials. Surface the most relevant options in the overview and retain all matched entries in the expanded table.

For each section, group each intervention's positive and negative records before ranking. Use the following explicit, lexicographic key (lower index first):

1. Target match: `exact_tested_strain`, `same_species`, `supported_mechanism_other_carrier`, `community_model_match`, `context_only`.
2. Endpoint match to the selected goal: measured harmful-behavior/viable-burden reduction or protective-function improvement as applicable; then direct matrix/adhesion change; then organism-load-only; then clinical outcome without a measured biofilm endpoint; then hypothesis. Do not call biomass reduction an eradication endpoint.
3. Site/model applicability: applicable controlled human-gut result; human-derived gut-community result; applicable animal-gut result; other gut-relevant in-vitro result; extraintestinal result; in-silico only. This differs from a blanket hierarchy in which an unrelated human trial automatically wins.
4. Delivery evidence for the **tested formulation**: active target exposure measured; relevant delivery/activity demonstrated indirectly; plausible but unverified; route/site mismatch. Unknown exposure remains eligible, with its uncertainty prominent.
5. Relevant opposing findings: reproducible favorable findings without identified direct conflict; mixed/uncertain; direct unfavorable in-vivo evidence. This is a relevance/review ordering, not a claim of quantified net benefit. A direct adverse result must be visible even when a positive paper ranks highly.
6. Assay-quality flag completeness and stable intervention ID. Funding does not alone decide rank.

An adverse-result card must not masquerade as an option recommended for use. Keep it in the matched list with “opposing evidence—review before considering.” Independently established individual contraindications move an otherwise relevant candidate to “not suitable with the recorded clinical context,” with the explicit reason; missing medical history stays unknown. Do not infer a safe supplement stack from rankings.

Quality flags: appropriate comparator, biological replication, growth confounding addressed, established versus developing biofilm distinguished, viable burden measured, dispersal measured, exposure/pH controls, relevant host/barrier outcome, commensal collateral effects, funding/conflicts, and human randomization/blinding when applicable. Values are `yes/no/unknown/not_applicable`, not an invented efficacy score.

Every option answers: **What finding does this address? What changed in the study? Was this in a gut model? Is oral/colonic delivery supported or only plausible? What opposing findings exist? What could be followed to know whether the intended outcome improved?**

### Required intervention evidence record

`intervention_id`, exact compound, formulation/composition, DOI, publication status, study type, organism/taxid, strain, body site, host/model, biofilm maturity, model type, endpoint, effect direction, effect size/comparator, growth normalization, medium pH/control, experimental exposure units/duration, contradictory sources, delivery evidence, local target exposure, collateral commensal effects, funding/conflicts, safety context, and sample-match strength.

Record experimental concentrations internally and in expanded evidence details as **laboratory exposures**, not suggested oral doses. Never convert a well concentration into a capsule dose. Any human study regimen requires an actual source, its context and a separate medication/formulation review; this module does not derive treatment regimens from gene abundance.

### Delivery evidence dimensions

1. Is the tested molecule the same as the product ingredient?
2. Is activity retained after manufacturing/storage and relevant digestion?
3. Is release gastric, small-intestinal, colonic, local/topical, or systemic?
4. Was intact active compound/enzyme measured at the target, or only a metabolite?
5. Are mucus, pH, food binding, proteases and contact time represented?
6. Was the target a monospecies laboratory biofilm, mixed community, animal mucosa or human biofilm?
7. Were protective-community collateral effects or viable dispersal measured?

Enteric coating is a legitimate technology and may improve delivery. It does not establish human colonic exposure without product-specific evidence. Ordinary garlic, purified allicin and a nanoemulsion are different interventions; black seed oil and purified thymoquinone are different; a researched enzyme blend is not evidence for every commercial “biofilm breaker.”

### Mandatory intervention seed evidence

The following records are required seed content, including their limitations and contrary results. Their presence is an evidence-retrieval requirement, not an instruction to take every listed substance. Source DOI/PMID links identify the actual study.

| ID | Intervention and direct source | Model/endpoint and interpretation | Gut delivery/implementation note |
|---|---|---|---|
| BF-I01 | Allicin: [2016 UPEC study](https://doi.org/10.3390/ijms17070979) | UPEC CFT073 and J96; 12–50 µg/mL tested at growth-sparing concentrations. At 50 µg/mL, pre-exposure reduced formation about 33% and 17%; 1 h post-treatment dispersed about 40% and 30% respectively. Crystal violet biomass, microscopy, adhesion/motility, RT-qPCR. This includes both prevention and established-biofilm dispersal, not universal sterilization. | Real positive allicin evidence. Urinary isolates on plates, not human intestinal biofilms. Do not equate an oral dose to these concentrations. Target-linked card can reference fimH/adhesion/curli context, but presence of fimH does not prove allicin susceptibility. |
| BF-I02 | Allicin delivery: [human supplement study](https://doi.org/10.3390/nu10070812) | 13 subjects; garlic products compared using breath allyl methyl sulfide metabolite. Enteric-tablet bioequivalence 36–104%, reduced 22–57% with high-protein meal; marked product variability. | Supports oral generation/delivery of allicin-derived compounds, not measured intact-allicin concentration at colonic biofilms. Enteric formulation is plausible and must not be ruled out; record disintegration site, alliinase activity, product-specific release evidence, and meal effects. Systemic metabolite exposure is not colonic target exposure. |
| BF-I03 | Allicin nanoemulsion + epsilon-polylysine [2025 study](https://doi.org/10.1016/j.foodchem.2025.142949), [PubMed](https://pubmed.ncbi.nlm.nih.gov/39842203/) | In vitro E. coli planktonic and mature biofilm disruption; food-preservation/raw-beef application. | Formulation research card only; does not validate oral allicin+polylysine treatment or gut efficacy. |
| BF-I04 | Berberine, Klebsiella: [2013 primary paper](https://pubmed.ncbi.nlm.nih.gov/24377137/) | Seven strong-biofilm clinical isolates from 35 screened. Berberine MBIC 0.0635 mg/mL; study also reports chitosan/eugenol the same MBIC, linoleic acid 0.0312 mg/mL, curcumin 0.25 mg/mL. Formation endpoint, not mature eradication. | Positive research options for these compounds can be included. Do not label berberine unsupported. Isolate context and culture concentration must travel with evidence. |
| BF-I05 | Berberine, E. coli [2019 study](https://doi.org/10.3389/fmicb.2019.02584) | Antimicrobial-resistant E. coli; sub-MIC berberine reduced formation and downregulated quorum-related genes including luxS/pfs/sdiA. In vitro. | Relevant formation/QS evidence; not proof mature intestinal removal. |
| BF-I06 | Berberine contradictory E. coli evidence [2025 study](https://doi.org/10.3389/fcimb.2025.1565714) | Repeated half-MIC exposure selected increased MIC (>32-fold), increased biofilm and csgD expression; csgD-overexpression experiment supports mechanism. | Mandatory linked counterevidence on E. coli berberine cards. Do not extrapolate it to “berberine never works”; do not ignore adaptation risk or imply prolonged low exposure is invariably helpful. |
| BF-I07 | Berberine + vancomycin/C. difficile [2020 study](https://doi.org/10.1007/s10096-020-03857-0) | Twelve strains; planktonic antibacterial/synergistic effects but half-MIC berberine or combination enhanced biofilm in particular strains. | Store planktonic efficacy and biofilm direction as separate outcomes. An antibacterial hit must not automatically become an antibiofilm hit. |
| BF-I08 | Berberine/Candida [2020 primary study](https://doi.org/10.2147/DDDT.S230857) | Candida isolates assessed separately in planktonic and biofilm conditions; antifungal and biofilm effects support experimental candidal-biofilm option. | Match exact Candida species and assay; oral/intestinal colonization is not candidiasis. Do not promote treatment solely from stool DNA. |
| BF-I09 | Carvacrol/thymol/geraniol [2022 K. pneumoniae study](https://doi.org/10.3390/antibiotics11020147) | Uropathogenic NDM-1-producing strains; 15 essential-oil components compared, thymol/carvacrol/geraniol performed best for antibacterial and biofilm activity. Formation assay. | Include promising options; concentrated purified compound differs from whole oregano/thyme oil. No validated human gut MBEC or local exposure. Essential oils can irritate mucosa; avoid oral use of products manufactured only for fragrance/topical use. |
| BF-I10 | Carvacrol/E. coli [2023 study](https://doi.org/10.1186/s12866-023-02797-x) | Carvacrol alone and with cefixime examined for antibacterial and antibiofilm effects. | Experimental option, antibiotic combination belongs to clinician-review context. Need preserve exact assay rather than infer eradication. |
| BF-I11 | Thyme/oregano mixed-species biofilms [2026 study](https://doi.org/10.1007/s10482-026-02284-z), [PubMed](https://pubmed.ncbi.nlm.nih.gov/41854771/) | K. pneumoniae and A. baumannii polymicrobial biofilm in vitro, oils alone/combination. | Broadens from single-species evidence; not intact human gut-community evidence. Import full methods before quantitative effect reuse. |
| BF-I12 | Thymoquinone [2024 CR-UPEC study](https://doi.org/10.1007/s12088-024-01231-8), [PubMed](https://pubmed.ncbi.nlm.nih.gov/39678958/) | Carbapenem-resistant uropathogenic E. coli: formation inhibition and eradication assays; 256 µg/mL antibacterial; 128 µg/mL reduced motility; changes in blaKPC/efflux/motility gene expression. | Purified thymoquinone not equivalent to arbitrary black-seed-oil dose/content. Gut formulation plausible; active local concentration unknown. |
| BF-I13 | Carvacrol-thymoquinone nanocarrier/Candida [2024 study](https://doi.org/10.1016/j.diagmicrobio.2024.116606), [PubMed](https://pubmed.ncbi.nlm.nih.gov/39586149/) | Vaginal C. albicans/C. glabrata isolates; formation prevention at half-MIC/MIC using 50 nm carriers. | Delivery-specific research; vaginal context, not human gut eradication. |
| BF-I14 | NAC/Klebsiella [2023 study](https://doi.org/10.1186/s12866-023-02969-9) | NAC, metformin and secnidazole at one-eighth MIC reduced formation/virulence expression. Mouse protection tested after experimental manipulation; not a human gut biofilm trial. | Candidate adjuncts, no assumption that stool genes establish need for a drug. |
| BF-I15 | NAC crucial pH qualifier [2026 study](https://doi.org/10.3390/microorganisms14020512) |34 K. pneumoniae strains for susceptibility, selected strains for crystal-violet biofilms. Acidic NAC + polymyxin B additive; neutralized NAC much weaker, significant activity generally only at 32 mg/mL. | Capture media pH and pH-matched controls. Colon/host buffering can erase activity achieved by acidification; do not propose acidifying patient's gut to laboratory pH. Study wording “disrupts” does not replace its actual formation-assay design. |
| BF-I16 | NAC human gastric salvage [2010 RCT](https://doi.org/10.1016/j.cgh.2010.05.006), [PubMed](https://pubmed.ncbi.nlm.nih.gov/20478402/) |40 refractory H. pylori patients with ≥4 failures, NAC before culture-guided antibiotics: eradication = 13/20 versus 4/20. Human gastric biofilm rationale and biopsy evaluation. | Stronger clinical context, but gastric H. pylori-specific salvage strategy; not NAC-alone or colonic dysbiosis cure. |
| BF-I17 | NAC human nonbenefit [2020 multicenter RCT](https://doi.org/10.1177/1756284820927306), [PubMed](https://pubmed.ncbi.nlm.nih.gov/32821287/) |680 treatment-naive H. pylori patients; adjunct NAC ITT eradication 81.7% versus 84.3%, no improvement. | Mandatory linked contrary study. Regimen and prior-treatment context differ from BF-I16; present both. |
| BF-I18 | Serrapeptase [2025 E. coli primary study](https://doi.org/10.3390/microorganisms13081875) |E. coli ATCC 25922; capsule-derived preparation tested for formation and 24 h treatment of preformed biofilms, biomass/viability/amyloid stains. Formation IC50 = 14.2 ng/mL; lower curli-associated signal. One strain; proposed direct CsgA/CsgB binding from docking, not biochemical proof. | Credible enzyme research card including mature-biofilm assay. Oral enteric delivery possible but active enzyme reaching human colon at effective concentration was not established. No dose conversion from ng/mL, no assumption fibrinolysis equals antibiofilm action. |
| BF-I19 | Lactoferrin/lysozyme/dextranase [primary study](https://doi.org/10.2436/20.1501.01.171), [PubMed](https://pubmed.ncbi.nlm.nih.gov/23844477/) |Single-species E. coli/K. pneumoniae MBEC-device biofilms, ATP-based viable-biomass readout; reductions varied; no agent destroyed both completely. | Matrix-targeted enzyme/lactoferrin candidates; protein digestion and intact local exposure matter. ATP reduction not histologic colon clearance. |
| BF-I20 | Lactoferrin opposite responses [2026 study](https://doi.org/10.1007/s11259-026-11512-w) |24 E. coli and 20 K. pneumoniae bovine/environment isolates. At 200/1,000 µg/mL growth decreased, K. pneumoniae biofilm inhibited but E. coli biofilm enhanced. | Mandatory organism/strain-direction counterevidence. Cannot call lactoferrin a universally selective beneficial biofilm remover. |
| BF-I21 | Lactoferrin commensal effects [2020 study](https://doi.org/10.1016/j.anaerobe.2020.102232), [PubMed](https://pubmed.ncbi.nlm.nih.gov/32634470/) |B. fragilis/B. thetaiotaomicron formation/laminin binding inhibited at 12.5 µg/mL, growth not inhibited by physiological concentration 2 mg/mL; antibiotic synergy negative. | Demonstrates matrix/colonization effects can extend to commensals. Record protective-community collateral effects. |
| BF-I22 | Enzyme/botanical mix BioDisrupt [2023 study](https://doi.org/10.4014/jmb.2212.12010) |Manufacturer-sponsored in vitro established biofilms: Candida and two Staphylococcus strains improved; P. aeruginosa response reversed at a high dose with increased biomass/activity; B. burgdorferi residual metabolic activity increased. Combination enzymes, NAC, cranberry, berberine, rosemary and peppermint. | Study supports exact-mixture research option, not all marketed “biofilm breakers”, all ingredients individually, or human systemic eradication. Product 1–40 mg/mL tested, no gut exposure proven. Funding and dose-direction reversal must display. |
| BF-I23 | Human SIBO/IMO exploratory add-on: [2025 study](https://doi.org/10.7759/cureus.99116) | Retrospective 13 patients: herbs alone (n=5) versus herbs plus adjunct (n=8; NAC n=2, enzyme blend n=6). SIBO eradication was 60% with adjunct versus 100% control, nonsignificant; IMO eradication was zero in both. Within-group gas changes; no symptom or direct biofilm measurement. Formulation-owner conflict disclosed. | Preliminary human evidence. Within-group significance is not between-group superiority; breath-gas changes do not demonstrate biofilm removal. |
| BF-I24 | S. boulardii CNCM I-745 [2022 study](https://doi.org/10.3390/microorganisms10061082) |Three C. difficile strains co-cultured with live yeast: reduced formation, thickness and eDNA. Genetically close S. cerevisiae did not reproduce effect; direct contact appeared required. | Useful exact-strain supportive-antagonism card, prevention rather than human mature-biofilm eradication; clinical CDI management separate. Yeast infection risk relevant for severely immunocompromised/critically ill/central-line patients. |
| BF-I25 | Polydextrose/Klebsiella: [2025 study](https://doi.org/10.1128/spectrum.01017-25) | Mouse prevention and laboratory work associated polydextrose with reduced bacterial loads and changes in TamA-related adhesion/biofilm mechanisms. | Relevant food/prebiotic research; not demonstrated removal of established human biofilms. Preserve organism, model and tolerability context. |
| BF-I26 | Fiber-supported communities: [2024 study](https://doi.org/10.1039/D4FO01294A), [PubMed](https://pubmed.ncbi.nlm.nih.gov/39082112/) | A nine-species gut consortium formed a stable biofilm on wheat fiber in vitro and withstood relevant environmental challenges. | More biofilm is not automatically bad. This supports a protective ecological hypothesis, not established human clinical benefit. |
| BF-I27 | EGCG/myricetin: [2018 study](https://doi.org/10.1038/s41598-018-26748-z) | E. coli K-12 curli-related formation and regulation assays; EGCG activity supports a targeted experimental option. | Beverage and concentrated extract differ. Human benefit and effective gut exposure are unknown; concentrated-extract safety context must remain visible. |
| BF-I28 | EGCG/Salmonella: [2024 study](https://doi.org/10.3389/fcimb.2024.1432111) | Laboratory quorum/biofilm inhibition and a mouse enteritis model. | Organism-specific research, not proof that green tea removes all human intestinal biofilms. |
| BF-I29 | D-amino acids: [2010 study](https://doi.org/10.1126/science.1188628) | Biofilm matrix-disassembly/formation experiments in defined laboratory organisms. | Experimental candidate; link the important strain-context and reproduction findings in BF-I30. |
| BF-I30 | D-amino-acid counterevidence: [2013 mechanism](https://doi.org/10.1128/JB.00975-13), [2015 reproduction study](https://doi.org/10.1371/journal.pone.0117613) | A laboratory-strain defect explained important susceptibility; another study did not reproduce inhibition in tested Bacillus/Staphylococcus contexts. | Preserve isomer, organism, strain and exposure. No universal D-amino-acid efficacy claim. |
| BF-I31 | Precision phage delivery: [2022 Cell study](https://doi.org/10.1016/j.cell.2022.07.003), [PubMed](https://pubmed.ncbi.nlm.nih.gov/35931020/) | A phage consortium suppressed IBD-associated K. pneumoniae and reduced inflammation in mice; artificial-gut and healthy-volunteer work supported delivery/safety research. | Pathobiont suppression and delivery are not direct human biofilm-eradication endpoints. Host range and treatment context require separate evidence. |
| BF-I32 | Targeted phages: [PSC model, 2023](https://doi.org/10.1038/s41467-023-39029-9), [high-alcohol-producing Klebsiella model, 2023](https://doi.org/10.1038/s41467-023-39028-w) | Targeted gut Klebsiella reduction improved experimental hepatobiliary/metabolic disease outcomes. | Promising organism-specific approaches; do not claim a measured biofilm endpoint where none was assessed. |
| BF-I33 | Natural products/C. difficile: [2024 study](https://doi.org/10.3390/pharmaceutics16040478) | Garlic oil, peppermint oil, curcumin and cinnamaldehyde-related compounds produced strain-dependent antibacterial/adhesion changes; biofilm architecture could also become denser/thicker. | Positive growth-inhibition findings cannot be converted into uniformly beneficial antibiofilm effects. |


### Additional required protective and systemic candidate cards

| ID | Evidence | Implementation |
|---|---|---|
| BF-I34 | [L. reuteri biofilm delivery, 2017](https://doi.org/10.3389/fmicb.2017.00489), BF-S15 | Protective/formulation research. The ATCC 23272/GtfW context is not a generic probiotic-species effect. Do not imply the microsphere formulation is naturally present in a donor. |
| BF-I35 | [SB-121 phase Ib human trial, 2023](https://doi.org/10.1038/s41598-023-30909-0) | Small randomized, placebo-controlled crossover study in 15 participants with autism; L. reuteri, dextran microparticles and maltose, with safety/tolerability and exploratory outcomes. It is relevant evidence that supportive biofilm-oriented formulations have reached human research, not validation of a stool protective score or proof of broad clinical efficacy. |
| BF-I36 | [F. prausnitzii/host-glycan work, 2026](https://doi.org/10.1080/19490976.2026.2665870) | Protective mechanism in experimental colitis and human tissue correlations. Glycan intervention/model findings are not a validated oral treatment for the current user or a unique DNA biomarker. |
| BF-I37 | [Serrapeptase/S. aureus, 2023](https://doi.org/10.1007/s00253-022-12356-5), [P. aeruginosa work, 2023](https://doi.org/10.1007/s00253-023-12772-1) | Strain-specific laboratory enzyme evidence; separate from the E. coli study. No automatic systemic or gut clearance claim. |
| BF-I38 | [Nattokinase/S. mutans, 2024](https://doi.org/10.3390/pathogens13040286) | Oral-bacterial research in 46 dental-plaque isolates. This is direct organism-specific research but not evidence of orally administered nattokinase clearing colonic biofilms. |
| BF-I39 | [Lumbrokinase/emodin biomaterial, 2023](https://doi.org/10.3390/bioengineering10080906) | Local combination/material context. Do not assign the result to oral lumbrokinase alone. |
| BF-I40 | [Anti-DNABII antibody work, 2023](https://doi.org/10.3389/fmicb.2023.1202215), [host HMGB1 work, 2021](https://doi.org/10.1172/JCI140527) | Investigational matrix-directed approaches. Evidence for increased susceptibility after release or model clearance is not proof that unsupervised dispersal is harmless. Clinical availability and route are separate from mechanistic relevance. |
| BF-I41 | [Antibiotics versus biofilm-targeted antibody microbiome effects, 2024](https://doi.org/10.1038/s41522-024-00481-0) | Duff et al.: oral/middle-ear antibiotic versus targeted-antibody microbiome effects in chinchillas. Preserve the animal model and local treatment context; this is not a human gut therapy or proof that all antibody approaches spare every commensal. |

The BF-I34–BF-I41 cards expand supportive and site-specific research coverage. They must retain the same organism, delivery and endpoint restrictions as the other seed cards.

### New gut-specific and protective intervention seeds

| ID | Candidate and direct source | Required interpretation and action target |
|---|---|---|
| BF-I42 | Grape-pomace polyphenols, [BF-S19](https://doi.org/10.1038/s41522-025-00716-8) | Human-derived prefrail-community experiment, not an oral clinical trial. The tested extract reduced dispersal and IL-6-inducing effects while biomass increased; epithelial adhesion, IL-8 and TNF did not significantly improve. Small donor subset. Show as **community stabilization / reduced adverse behavior**, not proven eradication or an ordinary grape-food equivalent. Oral exposure at the tested community conditions remains unestablished. |
| BF-I43 | Pomegranate extract, [BF-S24](https://doi.org/10.3390/nu15071771) | Barrier repair/experimental-colitis support. Biofilm effects on Enterobacterales were weak or absent in the tested system; do not infer robust Klebsiella/E. coli clearance. Product composition and gut-delivery uncertainty remain explicit. |
| BF-I44 | L. johnsonii studied isolate, [BF-S25](https://doi.org/10.3389/fimmu.2026.1749001) | Direct EPEC biofilm experiments plus mouse enteric-infection outcomes support a protective-antagonism research option. Do not substitute any commercial L. johnsonii strain without identity/equivalence evidence. Public 16S data do not provide an exact strain genome. Do not recommend putatively annotated metabolites as proven active ingredients. |
| BF-I45 | L. plantarum CIRM653, [BF-S23](https://doi.org/10.3920/bm2017.0002) | Include both promising in-vitro Klebsiella results and prolonged mouse intestinal carriage. Display prominently as **opposing in-vivo evidence**, not a selected decolonization treatment. Essential regression case for the action-ranking engine. |
| BF-I46 | Resveratrol plus L. paracasei ATCC334, [BF-S26](https://doi.org/10.3390/ijms21155423) | Strain-dependent promotion of adhesion/biofilm in vitro. Possible supportive formulation hypothesis, not proven human engraftment, beneficial biofilm quantity or an automatic reason to take concentrated resveratrol. |
| BF-I47 | Biofilm-state L. plantarum Y42 / associated matrix products, [BF-S27](https://doi.org/10.1021/acs.jafc.4c00460) | Protective barrier evidence in mice; exact strain, state and components matter. Ordinary planktonic products or generic EPS genes are not equivalent. |
| BF-I48 | Biofilm-state Lactiplantibacillus LR-1, [BF-S28](https://doi.org/10.1039/d3fo02733c) | Selected mouse-colitis benefits exceeded planktonic comparison. Include as experimental protective support, without a claimed human treatment effect. |
| BF-I49 | Multi-strain restorative consortium, [BF-S29](https://doi.org/10.1007/s11274-026-05071-0) | Experimental mouse recovery after antibiotics. Full tested composition/formulation and study methods must be imported before any exact-product matching. Co-detection of taxa in stool or donor material does not recreate the consortium. |
| BF-I50 | Colon-targeted anti-E. coli biofilm prodrug approach, [BF-S30](https://doi.org/10.1021/acsmedchemlett.5c00675) | Investigational delivery/therapy evidence. No validated human regimen or generally available supplement is established here. Link the paper; do not provide synthesis or unsupervised experimental treatment instructions. |

### Finding-to-action routing table

| Finding | Options that must be retrievable | Next useful information / success endpoint |
|---|---|---|
| Supported E. coli adhesion/amyloid-related mechanism, with relevant concerning context | Allicin BF-I01/02, berberine BF-I05/06 with the opposing study, relevant EGCG/myricetin BF-I27, serrapeptase BF-I18, and any other exact matched records. BF-I50 belongs in investigational options. | Confirm carrier/pathotype and clinical context as needed; follow raw signal/carrier abundance and symptoms separately. No promise that these remove a human intestinal biofilm. |
| Supported Klebsiella-associated mechanism | Exact Klebsiella records already in BF-I01–41; include essential-oil compound, NAC and other matched evidence with formulation/pH limitations. Always include BF-I45's counterexample where probiotic decolonization is discussed. | Distinguish colonization from infection; actual viable burden and relevant clinical assessment are more informative than a generic matrix reduction. |
| High H-C community-pattern proxy alone | Explain the two contributing taxa; retrieve source-matched dietary/community options, related measured functional findings and evidence that may clarify an established GI condition. | Gene/carrier resolution or an independently measured phenotype can clarify relevance. Proxy elevation alone is not an eradication instruction. |
| Low P-E ecological proxy | Existing individualized dietary-fiber/prebiotic support options with tolerability context; BF-I26, BF-I36 and relevant BF-I42/43/46–49 as explicitly experimental research. | Track the measured ecological pattern and tolerance/clinical outcomes. Do not diagnose absence of protective biofilm or promise stable engraftment. |
| Measured excessive community dispersal or host inflammatory effect | BF-I42 is directly relevant to the studied community behavior; BF-I43 is barrier-support context. Other options need endpoint-specific matching. | Repeat the **same measured endpoint** under comparable conditions, not just total biomass. |
| Protective strain/mechanism present alongside concerning signals | Show preservation-oriented options, exact probiotic evidence and collateral-effect information alongside target-specific candidates. | Improvement is not defined as disappearance of every biofilm-related gene. |
| Only generic matrix/housekeeping genes detected | Explain context and show the measured inventory; no organism-eradication plan follows from these genes alone. | More specific carrier/function or phenotype evidence, if independently appropriate. |

No relevant positive candidate is excluded merely because it lacks an RCT. No candidate is called effective in this person solely because a compound reaches the gut and affects one strain in vitro. Both rules are required. The report is a practical evidence navigator with ranked possibilities, rather than an unsupported treatment-response oracle.

### Clinical and procedural context

- The [UQ IBS page supplied by the user](https://imb.uq.edu.au/preventing-biofilms-ibs) describes developing research and observations around visible biofilms. BF-S01 supports the underlying association. It does not establish endoscopic lavage as a controlled, validated cure for IBS; do not provide home irrigation instructions.
- Established organism-specific treatment can matter more than a generic matrix-directed product. A confirmed gastric H. pylori infection has a different clinical context from stool E. coli carriage or an experimental amyloid gene signal.
- Bismuth-containing H. pylori therapy, novel bismuth materials, and ordinary over-the-counter bismuth are not interchangeable evidence. Add specific cards only with the actual compound, formulation, indication and study.
- DNase, dispersin-type enzymes and phage depolymerases are substrate-specific research strategies. Matrix reduction, viable-cell release and organism killing must be separate outcomes. No genetic readout predicts that an untested enzyme formulation will work in the gut.
- [A randomized rhDNase study in non-CF bronchiectasis](https://doi.org/10.1378/chest.113.5.1329) found worse outcomes, illustrating the danger of transferring a successful mucus/biofilm strategy across diseases. It does not imply all DNase approaches are harmful.
- [C. difficile gut-model work](https://doi.org/10.1038/s41522-021-00184-w) found biofilm-associated spores could persist despite an intervention context involving FMT. This is not proof that clinical FMT fails, nor that DNA detects viable spores.
- Do not generate an automatic supplement “stack,” antibiotic course, anticoagulant, or FMT instruction from these scores. Supply evidence-linked options and explicit missing information needed to choose between them.

## 11. Fibrin, clotting, and systemic biofilms

### 11.1 Substrates must remain distinct

| Term | Meaning | Relevant measurement |
|---|---|---|
| Microbial polysaccharide matrix | Organism-produced structural carbohydrates | Genes indicate capacity; actual matrix requires an appropriate material assay |
| Bacterial amyloid, including curli/BAP-associated structures | Particular microbial protein assemblies | Sequence potential, expression, or separately validated protein/structural assay |
| Fibrin | Host clotting protein formed from fibrinogen | Host protein/clotting assays or tissue characterization |
| Fibronectin | Host extracellular matrix/adhesion protein | Host-substrate interaction, distinct from fibrin |
| Extracellular DNA | DNA outside cells that may contribute to matrix | A compartment-specific assay; total stool DNA cannot distinguish its origin |
| Mucus | Host barrier material supporting spatial separation | Barrier/mucin/tissue measurements; not waste to strip away |
| Fibrosis | Host tissue remodeling/scarring | Organ-specific clinical assessment; not a microbial matrix-gene result |
| NETs | Host neutrophil extracellular structures | Host-specific immunologic/structural assessment, not generic microbial DNA |

A DNA hit in `csgA` cannot produce “fibrin detected.” A hit in a microbial fibrin-binding gene cannot produce “blood clots present.” Residual human `F2`, fibrinogen or `SDC1` reads cannot measure protein abundance or activity and should not be repurposed from incidental host contamination.

### 11.2 Required host-interaction context registry

Detecting a supported microbial mechanism is legitimate; its anatomical interpretation is constrained to the specimen and evidence.

| Target context | Interpretation | Sources |
|---|---|---|
| S. aureus Coa/vWbp | Host-coagulation recruitment capability when appropriate substrates/context exist | [2015 device model](https://doi.org/10.1093/infdis/jiv319), [2024 pseudocapsule study](https://doi.org/10.1016/j.bioflm.2024.100233), BF-S17 |
| S. aureus ClfA/ClfB, FnbA/FnbB | Host-protein binding/adhesion; preserve substrate distinctions | BF-S17 and curated protein annotations |
| S. aureus staphylokinase | Host fibrinolytic interaction, potentially related to dissemination | Context-only; fibrin degradation is not automatically a beneficial microbial trait |
| Broad anchoring/regulatory systems | Effects can differ by strain and experimental site | BF-S18; no monotonic harmful-score weight from one experiment |
| E. faecalis fibrinogen-associated adhesion | Evidence from particular urinary/device contexts | Retain site-specific evidence; stool detection does not diagnose those infections |

Reviewed protein identity examples for cross-checking curated annotations include ClfA Q2G015/Q53653, ClfB O86476/Q2FUY2, and FnbA P14738/Q2FE03/Q7A3J7. These are reference examples, not an exhaustive diagnostic panel or a claim that all sequence variants behave identically. Resolve versions and organisms through the pinned UniProt import.

### 11.3 Physiological versus excessive thrombin

BF-S12 and BF-S13 must be linked together. Physiological epithelial thrombin contributes to containment; excessive activity can disrupt community architecture and increase harmful behavior. This is a contextual host–microbe mechanism, not evidence that suppressing all thrombin or dissolving all fibrin improves the gut.

The 2026 study measured fecal thrombin as a protein and used additional experimental assays. It did not validate host activity inference from microbial genes. Increased flagellin accessibility/impact may occur without a corresponding increase in DNA or transcription; do not turn `fliC` abundance into a surrogate thrombin assay.

Add an optional `fecal_thrombin` external-assay type with method, units, activity-versus-protein distinction, experimental reference population and status. Until a laboratory provides a validated assay, the patient-facing value is “not measured,” not zero. Experimental intracolonic pharmacology in an animal model is not an oral human treatment protocol.

### 11.4 Long COVID and circulating material

[Long-COVID fibrinaloid/microclot proteomic work](https://doi.org/10.1186/s12933-022-01623-4) concerns host blood/plasma material. It does not establish that every microclot is a microbial biofilm or that stool sequencing can quantify circulating clots. Research hypotheses investigating such connections should be tagged `hypothesis`, not reported as validated diagnostic equivalence.

Do not create a “systemic biofilm burden,” “fibrin toxicity,” “microclot percentage,” or predicted clot location from stool reads. If independent blood/tissue results are imported, show them as separate measurements with their own evidentiary status.

### 11.5 Enzyme delivery and medication context

Serratiopeptidase is a zinc metalloprotease; it is not the same enzyme class or evidence base as nattokinase. A capsule-derived enzyme dissolved directly into a laboratory biofilm assay remains in-vitro evidence. Nattokinase oral studies suggest biological exposure/effects, but this does not establish a measured intact enzyme concentration in colonic or cardiac biofilms: [pharmacokinetic pilot](https://pubmed.ncbi.nlm.nih.gov/23709455/), [small pharmacodynamic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4479826/).

Keep possible bleeding/medication interactions attached to fibrinolytic candidates. A [reported severe bleeding case](https://doi.org/10.7759/cureus.20074) is a safety signal with confounding, not an incidence estimate or proof of causation in every user. Use actual current medication metadata; do not infer the user's present regimen from older conversation examples. Missing medication information is visible.

## 12. Antibiotic resistance, tolerance, and persistence

The report must explain three different ideas:

1. **Genetic resistance:** inherited or acquired determinants that can affect drug susceptibility; interpreted by the existing AMR system with its organism/linkage limits.
2. **Biofilm-associated tolerance:** survival under particular exposure conditions through matrix effects, growth gradients, physiological state or related mechanisms; not necessarily a heritable change in MIC.
3. **Persister behavior:** a surviving subpopulation under particular conditions; gene presence alone is not a measured persister fraction.

Do not display a universal “1,000 times more resistant” multiplier from biofilm genes or repeat such a number from an explainer as if it were the sample's measurement. Effects differ by organism, matrix, drug and assay. Matrix disruption may improve susceptibility but can also release viable cells; reduced staining does not prove that organisms were killed.

Co-detection of an AMR gene and matrix genes is clinically interesting context. It is not a measured minimum biofilm eradication concentration and does not identify a drug that will work. A carrier link requires actual supporting evidence. Clinical susceptibility testing and site-specific assessment remain separate inputs.

## 13. Report layout and exact explanatory copy

### 13.1 Front-loaded overview

Place the biofilm overview among the early visual metric pages, with cross-links to complete detail pages. Do not hide it after many pages of explanatory prose.

Required overview content:

- The two headings with up to four quantitative cards from §4, each with scope, percentile, reference n, coverage and clinical-validation status.
- “DNA measures potential; activity and location are not measured” when applicable.
- Three most informative supported findings, plus a count and link to **all** findings.
- Context-dependent mechanisms count, supported/assayed module counts and unresolved-target count.
- Most consequential uncertainty, e.g. carrier unresolved, sparse protective evidence or unmatched reference.
- Whether any direct assay was supplied, with date/site.
- The highest-ranked matched action candidates, their intended endpoint and strongest supporting/opposing evidence, with a link to the full list. No automatic multi-supplement plan.

Do not show only the first eight or first N detailed results. JSON, HTML and PDF must agree on complete inventory. A compact summary can prioritize findings; the detail section cannot silently omit them.

### 13.2 Detail-card fields

For every module show: finding, measured value and units, reference if valid, sequence/assay support, carrier/resolution, biological meaning, protective/harmful/context category, evidence maturity, what is not measured, related existing profile findings, and relevant intervention evidence. Users must be able to trace a conclusion back to the measurement and source.

For example, a curli-family finding reads “Genes associated with bacterial amyloid production detected” followed by carrier/confidence and supporting detail. It must not read “You have a harmful biofilm” or “You have brain amyloid.” For a supported BAD-associated pattern, say “Resembles a gene pattern reported in one bile-acid-diarrhea study”; the diagnosis and physical-biofilm status remain unestablished.

### 13.3 Ready-to-use education text

**What are biofilms?**  
“Biofilms are organized communities of microbes held together by a shared matrix. They are a common microbial way of living, not automatically a disease. In the gut, some communities can support stable colonization and a protective relationship with the intestinal lining. Others can help pathogens persist, grow too close to tissue, or participate in inflammation. Their effect depends on the organisms, their activity, their location and the host.”

**What did this test measure?**  
“This test ranks selected biofilm-related genetic capabilities and community patterns. We show concerning mechanisms and protective support separately, so one does not cancel the other. The community-pattern scores are experimental proxies. DNA alone does not directly measure an attached biofilm, its thickness, location or current activity.”

**Why can harmful biofilms matter?**  
“Some biofilms help organisms survive immune defenses or antimicrobial exposure. That can contribute to persistent or recurrent infection in the right setting. Biofilm-related tolerance is different from carrying an antibiotic-resistance gene, although both may occur together. A genetic finding alone cannot tell us which treatment will work.”

**Can harmful biofilms be reduced?**  
“Research includes targeted treatment of specific infections, dietary and community-supporting approaches, probiotics, compounds and enzymes that affect matrix or adhesion, and emerging precision therapies. Some have encouraging laboratory or animal evidence; others have human evidence in particular conditions. Each option below states exactly what was tested. Breaking up a matrix does not necessarily kill the organisms, and removing protective communities or mucus can be counterproductive.”

**What about fibrin and other organs?**  
“Fibrin is a human clotting protein that some microbes can exploit in blood or tissue infections. It is not the same as bacterial amyloid or every biofilm matrix. Biofilm-associated problems can occur on teeth, wounds, airways, urinary devices, implants and heart valves. Stool DNA does not locate or diagnose those problems or measure circulating clots.”

**What can I do with this result?**  
“Use the supported findings to identify which mechanisms deserve attention, which evidence-linked options are relevant, and what additional information would make a decision more reliable. Changes over time can show changes in the measured genetic signals; they do not by themselves prove that a biofilm has been removed.”

### 13.4 Intervention card example

**Allicin — experimental evidence relevant to E. coli biofilms**  
“Laboratory studies in particular E. coli strains found reduced biofilm formation and dispersal of preformed biofilms. These were not human intestinal clearance trials. Product formulation affects delivery, and the active concentration at a person's colonic biofilm is not established. This is a relevant candidate to review, rather than a predicted cure based on DNA.” Link BF-I01 and BF-I02; attach current medication/formulation context without hiding the positive findings.

**Berberine — positive and conflicting E. coli evidence**  
“Some experiments found less biofilm formation; another found increased biofilm and resistance after repeated sub-inhibitory exposure. The organism, strain and exposure conditions matter. These findings support an evidence-linked candidate, not a universal claim that berberine always removes or always promotes biofilms.” Link BF-I05 and BF-I06 together.

### 13.5 Accessibility and external links

Use readable labels, tooltips only as supplements, text equivalents for colors, and unclipped footnotes. All external HTML links must include `target="_blank" rel="noopener noreferrer"`. In generated PDF, retain clickable source links. Include source title/year and actual model on the card, not a bare DOI without explanation.

## 14. FMT matching and longitudinal integration

This module adds information to matching; it does not validate donor clearance, disease causation, or mixing benefits.

- Show each donor and recipient's supported molecular mechanisms, strain/locus resolution, assay gaps and independent potential displays.
- A protective mechanism in one donor cannot mathematically cancel a concerning mechanism in another. Do not blend scores and call the pooled product safer.
- A union of genes predicts only a possible combined inventory. Relative abundances are compositional; pooled material cannot be inferred by adding percentages. Engraftment, competition, dose, viability and host conditions remain unmeasured unless independently modeled/observed.
- A low harmful-associated score cannot override a confirmed pathogen or another validated exclusion reason already handled by donor screening.
- A common matrix gene, high research percentile, isolated curli hit or context-only fibrin-binding capability cannot by itself exclude a donor as disease-transmitting. Explain the finding and its evidence.
- Post-transfer strain/locus tracking can measure which genetic populations persist when sufficient depth exists. It does not prove attached biofilm establishment, clinical benefit or disease transfer.
- Report current sample age and assay date; a donor's old DNA profile does not establish the present composition or infectious status of material collected later.

Example required matching language: “These candidates differ in detected biofilm-associated mechanisms. The measurements do not establish whether their communities will form protective or harmful biofilms in this recipient. No validated numerical probability of transferring a biofilm-related disease can be calculated from these results.”

Longitudinal plots require comparable extraction/pipeline/reference versions. Show genes, carrier abundance and carrier-normalized signal separately when available; falling community relative abundance alone is not eradication. Antibiotics, supplements, bowel preparation and stool consistency should be visible timeline covariates. The module must not reinterpret a post-treatment nondetection as certain elimination.

## 15. CLI and machine-readable behavior

The following are **new interface requirements**, not claims that commands already exist. Adapt the executable prefix to the repository's established CLI while preserving these semantics.

```bash
openbiota biofilm inventory --sample SAMPLE2 --output-format json
openbiota biofilm analyze --sample SAMPLE2 --reuse-cache --registry biofilm-v08
openbiota biofilm analyze --sample SAMPLE4 --reuse-cache --registry biofilm-v08
openbiota biofilm analyze --sample SAMPLE6 --reuse-cache --registry biofilm-v08
openbiota biofilm explain --sample SAMPLE2 --include-candidate-evidence
openbiota biofilm compare --recipient SAMPLE2 --donors SAMPLE4 SAMPLE6
openbiota biofilm validate --registry biofilm-v08 --suite analytical-and-report
openbiota biofilm reference build --manifest approved-reference-manifest.json
openbiota biofilm rank-actions --sample SAMPLE2 --include-preclinical --include-opposing
openbiota biofilm datasets audit --registry biofilm-v08
openbiota biofilm longitudinal --samples SAMPLE2_BASELINE SAMPLE2_FOLLOWUP
```

`reference build` consumes an application-generated or user-supplied manifest listing independent sample IDs, source accessions, metadata, pinned assay versions, feature definitions, QC rules and hashes; this manifest is a runtime input, not a missing supplemental specification. It reuses compatible reference feature tables or reprocesses eligible available reads, enforces the frozen panel/domain rules and writes a versioned calibration with a reproducible audit. `rank-actions` reads existing structured results without rerunning sequencing. `datasets audit` reports assay/label/access compatibility for every BF-D source. `longitudinal` compares compatible time points and emits explicit version/domain differences.

`inventory` enumerates available reads/caches/assays/reference compatibility before expensive work. `analyze` produces full structured results even if some modules are unresolved. `explain` reads source-linked evidence and does not rerun alignment. `compare` preserves separate profiles without inventing a pooled-biofilm outcome. `validate` runs implementation checks and reports separately which scientific validations remain unavailable.

Add options for manifest path, output directory, compute budget, resume, explicit reference selection and optional external assay input using existing application conventions. Invalid references or incompatible caches must produce useful errors before computation. Process exit codes: 0 for completed analysis with explicit partial results; nonzero for corrupt inputs, failed mandatory schema validation or execution failure. A missing optional assay is not a process crash.

### Processing pseudocode

```python
def analyze_biofilm(sample, registry, reference=None):
    inventory = resolve_available_inputs(sample)
    evidence = reuse_or_run_compatible_assays(inventory, registry)
    calls = validate_and_link_molecular_calls(evidence, registry)
    modules = interpret_supported_modules(calls, registry)
    context = collect_context_without_inventing_measurements(sample, calls)
    axes = score_independent_axes(modules, reference, registry)
    proxies = score_named_experimental_proxies(inventory, reference, registry)
    options = retrieve_positive_and_contrary_evidence(modules, context, registry)
    findings = rank_supported_findings(modules, proxies, context, registry)
    ranked_options = rank_action_relevance(options, findings, registry)
    result = assemble_typed_result(sample, modules, axes, context, ranked_options,
                                  proxies=proxies, ranked_findings=findings)
    validate_result_schema_and_claim_limits(result)
    return result
```

No scoring or treatment branch may read free-form report prose as its source of truth when structured data exist. Legacy report-only imports are explicitly marked lower-resolution; never reverse-engineer missing gene calls from narrative text.

## 16. Validation and release design

### 16.1 Analytical validation

Use public reference-derived read fixtures, difficult close homologs, mixed communities, low-depth targets, fragmented assemblies, repeat-rich proteins, multiple strains, and gene-negative relatives. Evaluate correct-resolution precision/recall, cross-assignment, abundance error, locus-completeness error and unresolved frequency. Short reads with no discriminating bases must remain unresolved even if a classifier offers a best hit.

Publish target-specific analytical performance in the application release metadata. Gene-family validation is not strain validation. A single laboratory isolate cannot establish performance across every gut carrier. Include fungal lower-depth cases and protein repeats; do not hide those under a bacteria-only benchmark.

### 16.2 Biological validation

Separate three questions:

1. Can the software identify the DNA feature accurately?
2. Does that feature track an experimentally measured biofilm phenotype in the specified organism/context?
3. Does a stool measurement predict a clinically relevant mucosal phenotype or outcome in independent people?

Success at one question does not answer the next. Human disease labels alone are not gold-standard biofilm labels. Predicting PD, IBS or BAD is not necessarily predicting biofilm architecture.

A future clinical validation dataset should pair stool WGS with prospectively described site-specific biofilm assessment, relevant host measurements, and covariates. Include healthy participants and disease controls, pre-bowel-preparation stool where feasible, and blinded outcome assessment. Avoid treating all healthy controls as biofilm-negative. Account for extraction/center and medications, and hold out entire cohorts. Prospective collection is a research dependency, not something the coding agent can simulate by inventing labels.

### 16.3 Model evaluation

For any future trained predictor, report external calibration, discrimination, uncertainty, missingness and age/domain transport. Use nested subject-grouped splits; fit scalers/feature selection only on training data; serialize all transforms. Compare against taxonomy-only, organism abundance, basic covariates and simple module baselines to determine whether new genes add useful information. Do not report training accuracy as validation.

If tissue imaging and stool genes disagree, preserve both rather than forcing agreement. Evaluate whether a model predicts medication/bowel preparation rather than biofilm. Biofilm-positive tissue in one region is not an unbiased label for the entire intestine.

### 16.3.1 Executable retrospective validation work packages

The coding agent must attempt the available work, not leave “validate later” as the entire implementation:

1. **BF-D01 community associations:** acquire open 16S data and available metadata; reproduce genus-level feature associations separately in stool and biopsy, with participant grouping and disease/age/bowel-preparation covariates. The reported molecular subsets included 51 biofilm-positive and 54 biofilm-negative stool samples, and 35 positive/38 negative biopsy samples. Verify the deposited usable mapping rather than forcing these counts. BF-M18/19 feature directions come from biopsy associations; test whether direction and discrimination persist in stool. Do not assume they do.
2. **BF-D09 community behavior:** distinguish original feces, cultured communities and mouse samples; relate available taxonomic tables to dispersal/host-effect endpoints only where participant/condition mapping exists. Stratify adult/robust-aged/prefrail groups; do not let age alone masquerade as an adverse-biofilm signature. Restrict any intervention analysis to actual independent donor counts.
3. **BF-D10 tissue/oncotraits:** import available anonymous supplementary records. Evaluate biofilm and damaging-function associations as different endpoints. Raw/processed sequencing access limitations may block a subanalysis; emit its dependency record, while completing available table-based analyses. No bypass of access restrictions or patient privacy.
4. **BF-D02/03 BAP and BF-D07 BAD:** reproduce molecular feature prevalence/association when usable data exist. Disease association remains distinct from physical-biofilm validation. Do not report classification of IBS/PD/BAD as a validated biofilm classifier.
5. **Model comparisons:** compare fixed proxies, taxonomy-only models, mechanism-only models, basic clinical covariates, and combined models. For binary phenotype endpoints use participant-grouped nested evaluation and report AUROC, AUPRC with prevalence, calibration if probabilities are claimed, sensitivity/specificity at prespecified operating points, and uncertainty. For continuous behavior use MAE/rank association with intervals and assay reproducibility. No single cohort's training performance qualifies as external validation.
6. **Transport and falsification:** test disease/site/medication leakage, batch effects, unrelated dysbiosis controls, and whether a proxy adds value beyond its dominant taxon. Hold out cohorts where possible. A failed replication remains in the report/evidence registry and changes the relevant interpretation; it must not be deleted to preserve a desirable score.

A reference percentile can ship as a descriptive experimental ranking before phenotype validation. A prediction of **actual biofilm state** requires evidence on that endpoint and population, with prespecified acceptable error and calibration established by the intended-use validation plan. No arbitrary AUROC alone is treated as clinical adequacy. After retrospective work, propose only the additional paired observations needed to resolve specific failures; do not promise that more computation alone creates physical-state truth.

### 16.4 Release levels

- **Release A — useful research report:** complete module census; compatible DNA analyses; named single-/multi-mechanism rankings wherever calibrated; H-C/P-E proxies from compatible existing profiles; ranked findings/actions; explicit candidates; optional assay intake. Missing optional sources must not suppress computable results.
- **Release B — expanded and tested research scores:** independent frozen reference panels, retrospective replication results, calibrated partial masks and additional direct measurements. Report failed or limited transport alongside successful findings.
- **Release C — phenotype prediction:** only for models with actual assay-matched external validation. Name the predicted phenotype and population; do not rename it generic “biofilm health.”

These levels allow useful functionality to ship without hiding research evidence or pretending an unavailable clinical validation is already complete.

## 17. Acceptance tests

The following tests are mandatory behavior contracts. They test the software and report; they do not constitute clinical validation.

| Test | Area | Required outcome |
|---|---|---|
| BF-T001 | Input and detection | A taxonomic species profile without reads/gene calls yields taxonomic context and explicit unassayed loci, never invented accessory genes. |
| BF-T002 | Input and detection | MetaPhlAn bowtie2 summary files are not treated as SAM/BAM or a complete biofilm locus assay. |
| BF-T003 | Input and detection | A compatible existing gene cache is reused; changing registry sequence hashes invalidates the relevant detection cache. |
| BF-T004 | Input and detection | Changing only an evidence-card explanation does not trigger a full FASTQ realignment. |
| BF-T005 | Input and detection | Read pairs remain together during mapping, filtering and bootstrap resampling. |
| BF-T006 | Input and detection | A common short domain with multiple plausible carriers stays ambiguous. |
| BF-T007 | Input and detection | A database-missing target reports not_assayed, not absent. |
| BF-T008 | Input and detection | Low species depth cannot produce a confident strain-negative or intact-locus-negative result. |
| BF-T009 | Input and detection | An isolated regulatory-gene match cannot stand in for a complete structural mechanism. |
| BF-T010 | Input and detection | Different locus components on unlinked carriers do not create a fictional complete operon. |
| BF-T011 | Input and detection | Assembly collapse of multiple strains produces uncertainty rather than a manufactured single strain. |
| BF-T012 | Input and detection | Fungal low-depth/paralog fixtures preserve unresolved results and do not inherit bacterial detection limits. |
| BF-T013 | Input and detection | Gene-family homology cannot automatically become an exact allele/strain match. |
| BF-T014 | Input and detection | Gene abundance units and normalization versions must match the reference; incompatible units fail scoring but retain raw data. |
| BF-T015 | Input and detection | Nondetection, below-resolution, unassayed and failed-QC states survive all exports as different states. |
| BF-T016 | Input and detection | No absolute cells-per-gram or viable-cell result is derived from relative DNA abundance alone. |
| BF-T017 | Biological interpretation | A curli-related hit does not produce brain amyloid, Parkinson disease, fibrin or active-biofilm claims. |
| BF-T018 | Biological interpretation | BAP untested homologs remain distinct from the experimentally evaluated reference domains. |
| BF-T019 | Biological interpretation | F. prausnitzii abundance alone does not activate P-M; genus abundance may contribute only to the explicitly labeled fixed P-E proxy. |
| BF-T020 | Biological interpretation | Generic cobalamin, butyrate or fiber-metabolism genes are not used to fill missing protective-biofilm markers. |
| BF-T021 | Biological interpretation | L. reuteri species presence alone does not establish the studied GtfW strain/formulation phenotype. |
| BF-T022 | Biological interpretation | S. boulardii CNCM I-745 evidence is not automatically assigned to all S. cerevisiae strains. |
| BF-T023 | Biological interpretation | Community bssS/pga enrichment is not attributed to R. gnavus without linkage. |
| BF-T024 | Biological interpretation | BssS presence is not treated as a universally monotonic biofilm-forming effect. |
| BF-T025 | Biological interpretation | BF-D01 is imported as 16S/phenotype context, never as WGS classifier training. |
| BF-T026 | Biological interpretation | The BF-S01 biopsy bloom threshold is not reused as a stool decision threshold. |
| BF-T027 | Biological interpretation | A disease-labelled cohort is not relabelled biofilm-positive without physical-phenotype evidence. |
| BF-T028 | Biological interpretation | BBSdb protein classification performance is not advertised as patient diagnostic accuracy. |
| BF-T029 | Biological interpretation | A capsule annotation alone cannot create an attached-biofilm finding. |
| BF-T030 | Biological interpretation | The same sequence evidence in multiple databases appears once with multiple provenance links. |
| BF-T031 | Scoring | The midrank examples [0,0,2,4] produce 25, 62.5 and 100 for 0, 2 and 5, respectively. |
| BF-T032 | Scoring | Identical sample/reference values remain tied through both module and axis transforms. |
| BF-T033 | Scoring | A constant reference produces reference_non_discriminating and no invented meaningful percentile. |
| BF-T034 | Scoring | Nondetection encoded as zero retains left-censoring; unknown/unassayed observations remain null. |
| BF-T035 | Scoring | Protective 90 and harmful 90 remain separate; no net or subtractive score is computed. |
| BF-T036 | Scoring | An inactive/candidate-only module cannot contribute to an active axis. |
| BF-T037 | Scoring | One biological group receives an available_limited_scope score with the mechanism named; it is never titled comprehensive ecosystem health. |
| BF-T038 | Scoring | Correlated genes within one locus or PNAG dependence group do not multiply evidence weight. |
| BF-T039 | Scoring | Missing modules cannot be silently replaced by zero or dropped from the denominator. |
| BF-T040 | Scoring | A partial panel requires its own frozen mask/calibration and at least one group; zero groups cannot produce a mechanism score. |
| BF-T041 | Scoring | Mask selection never depends on keeping only favorable, unfavorable or detected findings. |
| BF-T042 | Scoring | The same sample receives identical outputs when analyzed alone or in a donor batch. |
| BF-T043 | Scoring | SAMPLE2/SAMPLE4/SAMPLE6 cannot be used as the complete healthy calibration cohort. |
| BF-T044 | Scoring | Adult-only calibration is flagged out_of_domain for a child and cannot generate pediatric normal/abnormal claims. |
| BF-T045 | Scoring | Separate axes and partial masks retain separate calibration IDs/hashes. |
| BF-T046 | Scoring | Paired-fragment and participant bootstraps have deterministic seeds and separate interval labels. |
| BF-T047 | Scoring | A summary-only cache cannot yield an invented fragment-resampling interval. |
| BF-T048 | Scoring | Reference intervals are not renamed clinical risk or treatment-response confidence intervals. |
| BF-T049 | Scoring | NaN, infinity, out-of-range percentiles and unresolved source IDs fail schema validation. |
| BF-T050 | Scoring | A valid analysis with no calibrated axes still produces the full measured/candidate report. |
| BF-T051 | Interventions | Relevant E. coli/allicin laboratory evidence appears without requiring a human RCT. |
| BF-T052 | Interventions | Berberine positive E. coli evidence and the 2025 adaptation/biofilm counterevidence appear together. |
| BF-T053 | Interventions | Berberine antibacterial activity plus increased C. difficile biofilm remain separate outcomes. |
| BF-T054 | Interventions | The NAC 2026 pH qualification survives both summary and detailed card rendering. |
| BF-T055 | Interventions | A formation-only assay is not labelled mature-biofilm eradication because of its title. |
| BF-T056 | Interventions | Crystal-violet decline alone cannot establish viable killing or regrowth-free eradication. |
| BF-T057 | Interventions | Serrapeptase capsule material added directly to culture remains in-vitro evidence with unverified intestinal exposure. |
| BF-T058 | Interventions | Laboratory concentration is not automatically converted into an oral dose. |
| BF-T059 | Interventions | Enteric delivery is accepted as plausible while missing measured target exposure stays unknown. |
| BF-T060 | Interventions | Human allicin-derived metabolite detection is not renamed intact colonic allicin exposure. |
| BF-T061 | Interventions | A blend result is not assigned to each ingredient, and high-dose adverse-direction findings remain visible. |
| BF-T062 | Interventions | Cureus 2025 preserves n = 13, no direct biofilm/symptom measurement, nonsignificant comparison and correct 60% versus 100% SIBO direction. |
| BF-T063 | Interventions | NAC H. pylori salvage benefit and larger first-line nonbenefit are both present with different contexts. |
| BF-T064 | Interventions | Lactoferrin organism-specific opposite biofilm directions are not averaged into universal benefit. |
| BF-T065 | Interventions | D-amino-acid contradictory/reproduction evidence remains attached to the candidate. |
| BF-T066 | Interventions | E. coli antibiofilm experiments cannot create an R. gnavus-specific eradication claim. |
| BF-T067 | Interventions | Antimicrobial susceptibility is not predicted merely from the presence of a target mechanism. |
| BF-T068 | Interventions | Phage pathogen-load improvement is not labelled direct biofilm clearance when biofilm was not measured. |
| BF-T069 | Interventions | A protective-community signal cannot automatically trigger broad antimicrobial options. |
| BF-T070 | Interventions | Missing medication/formulation information is displayed and does not erase the underlying evidence card. |
| BF-T071 | Fibrin and systemic | Stool coa/vwb capability cannot produce a heart, blood or intestinal fibrin-biofilm diagnosis. |
| BF-T072 | Fibrin and systemic | Host F2/SDC1 DNA cannot generate host protein/activity results. |
| BF-T073 | Fibrin and systemic | Physiological and excessive thrombin findings preserve their different contexts. |
| BF-T074 | Fibrin and systemic | Ex-vivo RNA cannot be deployed as a stool-DNA model without a validated assay-transport record. |
| BF-T075 | Fibrin and systemic | Staphylokinase does not become a protective signal simply because it interacts with fibrinolysis. |
| BF-T076 | Fibrin and systemic | Serrapeptase, nattokinase and lumbrokinase remain separate evidence identities. |
| BF-T077 | Fibrin and systemic | No gut/systemic biofilm removal claim comes solely from a fibrinolytic biomarker change. |
| BF-T078 | Fibrin and systemic | Long-COVID blood microclot research does not create a stool microclot percentage. |
| BF-T079 | Fibrin and systemic | AMR, biofilm tolerance and persister state are separate fields; no universal 1,000-fold multiplier appears. |
| BF-T080 | Fibrin and systemic | Independent-site clinical assays retain body site/date and cannot overwrite stool molecular findings. |
| BF-T081 | Report and integration | Overview contains independent panels, context, coverage and links to all detail findings. |
| BF-T082 | Report and integration | Null/unavailable values cannot appear as green zero scores. |
| BF-T083 | Report and integration | JSON, HTML and PDF contain the same full module and evidence inventory without fixed-N truncation. |
| BF-T084 | Report and integration | Every external HTML source link includes target=_blank and rel=noopener noreferrer. |
| BF-T085 | Report and integration | All figures have readable text/color-independent meaning, visible reference population and units. |
| BF-T086 | Report and integration | Two donors cannot cancel each other’s concerns or receive an invented pooled safety score. |
| BF-T087 | Report and integration | Generic matrix genes do not by themselves exclude donors; confirmed pathogen handling remains in its existing validated workflow. |
| BF-T088 | Report and integration | Post-treatment nondetection is not labelled certain eradication; relative decreases are not absolute depletion. |
| BF-T089 | Report and integration | Post-FMT shared strain/locus detection does not establish transfer of a disease or physical biofilm. |
| BF-T090 | Report and integration | Default processing keeps patient reads local and persists reproducible provenance without exposing host sequence data. |

### Revision 08.0.1 additional acceptance tests

| Test | Area | Required outcome |
|---|---|---|
| BF-T091 | Evidence | BF-S01 stores 1,426 screened versus 1,112 included; the denominators are not conflated. |
| BF-T092 | Evidence | BF-S20's biopsy WGS and stool qPCR are distinct assays; 265 specimens are not 265 people or stool-WGS samples. |
| BF-T093 | Evidence | BF-D10/11 request-only or unmapped data do not become fabricated public training cohorts. |
| BF-T094 | Evidence | BF-D09 human-derived communities, fecal inputs and mouse samples remain separate; technical replicates are not independent donors. |
| BF-T095 | Evidence | BF-D12 compares multi-species with single-species biofilms, not biofilm with planktonic state. |
| BF-T096 | Evidence | PV739486 is typed as 16S; it cannot satisfy a complete-genome or exact protective-strain requirement. |
| BF-T097 | Scoring | H-C and P-E have distinct IDs, full feature definitions, research labels and independent calibrations. |
| BF-T098 | Scoring | P-E never populates P-M or claims measured beneficial-biofilm quantity. |
| BF-T099 | Scoring | A qualified one-feature mechanism panel returns a limited-scope percentile without requiring clinical validation. |
| BF-T100 | Scoring | Generic adhesion/matrix positivity does not acquire a harmful direction without an explicit source-context claim. |
| BF-T101 | Scoring | Both H-C features and all four P-E genera are required for their named full proxy; missing assay is not zero. |
| BF-T102 | Scoring | H-C preserves the Escherichia/Shigella-to-E. coli measurement change and does not infer pathogenicity. |
| BF-T103 | Scoring | Taxonomic synonyms and genus sums do not double-count reads or species leaves. |
| BF-T104 | Scoring | Constant reference or constant combined index produces an explicit nondiscriminating state. |
| BF-T105 | Scoring | Oppositely ordered features with a constant aggregate cannot produce a fabricated meaningful composite rank. |
| BF-T106 | Scoring | Same panel values with different biological labels share arithmetic but never interchangeable interpretation. |
| BF-T107 | Scoring | Higher H-M does not decrease because P-M/P-E rises; no cancellation anywhere in reporting or donor matching. |
| BF-T108 | Findings | Linked concern context requires supported carrier linkage; unlinked community co-presence stays unlinked. |
| BF-T109 | Findings | Findings and actions have deterministic ranks, stored reasons and source IDs; no efficacy probability is invented. |
| BF-T110 | Actions | Allicin/berberine preclinical cards remain retrievable without human trials; opposing results accompany the positive ones. |
| BF-T111 | Actions | BF-I45 shows improved in-vitro biofilm findings and unfavorable mouse-carriage findings together. |
| BF-T112 | Actions | BF-I42 can show a favorable behavior endpoint despite increased matrix biomass. |
| BF-T113 | Actions | BF-I43 is not described as proven strong Enterobacterales biofilm clearance. |
| BF-T114 | Actions | BF-I44's putative metabolite annotations cannot become validated active-ingredient recommendations. |
| BF-T115 | Actions | Exact probiotic strains/formulations are not silently replaced by any product with the same species. |
| BF-T116 | Actions | Delivery plausibility allows an experimental candidate, but never becomes proof of effective human colonic exposure. |
| BF-T117 | Actions | A high proxy alone yields evidence/navigation options, not an automatic antimicrobial stack or eradication instruction. |
| BF-T118 | Longitudinal | Reduced DNA abundance is not reported as eradication; increased viable dispersal prevents a biomass-only success claim. |
| BF-T119 | Validation | Paired tissue samples/repeated observations never cross participant-level evaluation folds. |
| BF-T120 | Validation | A failed stool replication or domain shift remains visible and cannot be hidden by a successful tissue result. |
| BF-T121 | Report | Four cards remain separately labeled beneath two headings; coverage is not formatted as health probability. |
| BF-T122 | Report | Missing RNA, microscopy, request-only data or P-M cannot suppress a computable DNA module, proxy or matched action card. |
| BF-T123 | Report | All 50 intervention seeds, 24 modules, 32 core sources and 14 dataset records are represented with actual status. |
| BF-T124 | Report | External links in HTML include target=_blank and rel=noopener noreferrer; PDF links remain clickable. |
| BF-T125 | Validation | The embedded rank arithmetic passes fixtures for ties, single groups, missing masks, nonfinite values and constant aggregates. |


## 18. Implementation sequence and definition of done

### Increment 1 — inventory, provenance and truthful output

Inspect the repository, map existing schemas and sample IDs, implement the complete module census and typed result schema, and produce a report from compatible existing measurements. All unavailable analyses are explicit. Import source and intervention seed records from this specification. Correct prior report language that treats all biofilms as harmful or lack of trials as absence of evidence.

### Increment 2 — targeted observational analysis

Build versioned reference manifests from the named public resources, implement compatible whole-read gene analysis and existing strain/carrier adapters, validate each activated molecular assay, and preserve candidate results. Add gene-family, carrier-normalized and module measurements without inventing physical state. Integration is additive to existing caches; do not discard past sample analysis.

### Increment 3 — references and independent displays

Acquire/reprocess eligible public or existing reference cohorts, freeze transforms and dependence groups, run the scoring fixtures, implement dual panels and partial-calibration behavior. Implement the four score cards and the §8.7 starter proxies, with named limited-scope mechanism results rather than a two-group stop. Build an applicable reference as an explicit task; do not substitute the local recipient/donor set. Record scientific validation dependencies in application metadata.

### Increment 4 — complete evidence and report experience

Implement positive/negative intervention cards, oral-delivery distinctions, fibrin/systemic explanations, optional laboratory/imaging imports, longitudinal view, and donor comparison. Verify all sources and counts, accessible graphics and complete detail pages. Use actual source-linked claims rather than AI-generated citations or extrapolated efficacy estimates.

### Increment 5 — independent validation and maintenance

Run analytical benchmarks and all BF-T tests; evaluate held-out data where labels exist. Separate computational test results from clinical validity. Pin databases, code, schemas, evidence and calibration manifests. Provide a command to show what changed between releases and whether a changed source invalidates only report wording, interpretation, reference calibration, or molecular calls.

Definition of done: the CLI runs end to end on the existing sample inputs with complete per-module results/statuses; deterministic schemas and score fixtures pass; all 50 seed intervention IDs, 24 mechanism/proxy/assay IDs, 32 core source IDs and 14 dataset IDs are accounted for; every primary claim links to its source and retains model/endpoint; the four score cards under two headings, ranked findings, action options and context panel render without truncation; donor and longitudinal views retain the correct limits; the application explains why a value is unavailable rather than pretending it is zero. Actual molecular results, runtime and cohort performance must be reported by the coding agent after execution, not preclaimed by this specification.

## 19. Research provenance, source scope and remaining unknowns

Research included the user's five starting sources, their primary references, human endoscopic and stool studies, biofilm-associated genetic/protein resources, strain-level genotype–phenotype work, protective-community experiments, host thrombin/fibrin mechanisms, oral-delivery studies, direct compound/enzyme assays, contradictory findings and published precision approaches available by the cutoff. Source identities and selected full texts were checked using publisher/PubMed/PMC/Europe PMC records. Reviews and institutional explainers were used to locate and interpret primary work, not as substitutes for tested patient-level classifiers.

The user's starting links are explicitly incorporated:

- [PMC9884580](https://pmc.ncbi.nlm.nih.gov/articles/PMC9884580/): BF-S02, protective regulation/dispersal and therapy context.
- [UQ: preventing biofilms in IBS](https://imb.uq.edu.au/preventing-biofilms-ibs): research explainer and procedural-evidence limits.
- [PMC11391705](https://pmc.ncbi.nlm.nih.gov/articles/PMC11391705/): BF-S03, intestinal biofilm review.
- [ScienceDirect S2468451121000167](https://www.sciencedirect.com/science/article/pii/S2468451121000167): BF-S04, biophysical determinants.
- [Biofilms: The Good, the Bad & the Groundbreaking](https://pediatricsnationwide.org/2024/09/23/biofilms-the-good-the-bad-the-groundbreaking/): protective probiotic delivery and matrix-targeted research, with primary sources in BF-I34–BF-I41. Company development statements in this 2024 article must not be presented as verified current trial status.

This is a broad targeted evidence synthesis, not a claim that every biofilm publication worldwide was retrieved or that a clinical systematic review/meta-analysis was completed. No clinical good/bad stool-DNA score was validated on SAMPLE2, SAMPLE4 or SAMPLE6 in this work. The newly specified community proxies and ranking rules are explicit engineering hypotheses, not coefficients attributed to their source papers. Literature review through 18 September 2026 added BF-S19–32, including published 2026 protective/delivery research and opposing intervention findings. Direct active-state information remains limited by the input assay. The executable response to that limitation is to measure and explain supported mechanisms now, preserve all relevant intervention evidence, and add actual activity/spatial measurements when available.

### Required ongoing evidence-update behavior

Every new study is ingested as a source-specific claim with model, organism, site, assay, phenotype and direction. New publications cannot automatically change user scores. Sequence/annotation changes require analytical checks; changed module membership or weights require a new calibration; retractions/corrections invalidate affected claims and regenerate reports. Record negative results and conflicts alongside positive evidence. Deduplicate preprints against their later peer-reviewed paper so one study is not counted twice. Keep archived report provenance reproducible.

**End of specification. No supplemental handoff files are required.**
