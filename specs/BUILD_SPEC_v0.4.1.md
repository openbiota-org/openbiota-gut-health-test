# BUILD SPEC v4.1 — Chronic Urticaria Gut-Profile System

Revision: 4.1.1 · Evidence checked through 7 September 2026 · Research-use implementation specification

## 0. Coding-agent contract

Implement this specification in the existing microbiome-report repository. **This is the complete handoff for the chronic-urticaria increment. No supplemental research report, claim ledger, or previous build-spec document is required.** All necessary evidence, feature definitions, formulas, acquisition instructions, implementation steps, and tests are included here. Normal application source files, registries, datasets, model artifacts, and tests are expected implementation outputs—not additional handoff documents.

This revision replaces the earlier v4.1 specification. Preserve existing disease profiles, functional panels, reference logic, graphics, and sample results except for the chronic-urticaria changes explicitly specified here. Adapt the contracts to the actual repository; do not silently implement an entire assumed v04 architecture.

Deliver operational research profiles, not only literature cards. Compute the defined index when an eligible sample has usable features and a compatible frozen reference. Display partial coverage and bounds. Missing analytical prerequisites produce explicit unavailable states, never fabricated references or model weights. Implementation completion and independent clinical validation are different milestones.

### 0.1 Tasks

| Stable task ID | Required output | Evidence boundary |
|---|---|---|
| `CSU_STOOL_WGS_V1` | CSU-associated directional pattern index, 15-species panel plus eight-species sensitivity result | Small human shotgun discovery cohort; not an externally validated diagnosis |
| `SD_STOOL_TRANSLATED_V1` | Symptomatic-dermographism-associated directional pattern index | Small human 16S study translated to measured shotgun features; transport unvalidated |

Register by set union with the installed task IDs. A 33-task baseline becomes 35; a 30-task baseline becomes 32. Do not add unrelated profiles to force a count. Migration is idempotent. Existing CU IDs are revised in place; invalidate only their obsolete score/model caches and retain historical provenance.

Add non-task evidence sections for CSU+SD coexistence, antihistamine response, omalizumab response, functional potential, and interventions. These are not extra diagnostic or response probabilities.

### 0.2 Corrections to the previous v4.1

1. **Bacteroides stercoris is higher; Bilophila wadsworthia is lower** in the primary CSU differential analysis. The earlier manifest had these wrong. B. stercoris is not a verified figure/text direction conflict.
2. Only **8/15 species have BH q<0.05** in the numerical source table. Retain all 15 as exploratory; expose the strict subset and exact q-values.
3. MaAsLin2 association coefficients are not deployable classifier weights. Do not apply them to new abundances and call that the published model.
4. Remove the invented 60% taxonomy/25% function/15% consensus composite. Keep correlated taxonomy, pathways, metabolites, and clinical observations separate.
5. Pure SD data are **OEP001229**; **OEP001891** is CSU+SD. Cohort independence must be checked rather than assumed.
6. Add the overlooked 2026 human 16S cohort: **PRJNA1096288** is human amplicon data, **PRJNA1101664** mouse amplicons, **PRJNA1359023** a probiotic isolate genome—not patient shotgun.
7. The 2026 JACI Global MR paper is not new patient multi-omics; its associations did not survive FDR<0.05.
8. A genuine human CU FMT case exists, but does not establish efficacy. The eczema FMT report previously considered is not CU evidence.
9. Symptom duration must not block research scoring of asymptomatic eligible adults. Microbiome scores cannot establish clinical chronicity.
10. Replace ambiguous percentiles, reference fitting, missingness, and uncertainty with the exact contract below.

## 1. Evidence decision and subtype coverage

Direct human evidence is sufficient for **exploratory CSU-associated stool-pattern scoring**, particularly [Zhu et al., Nature Communications 2024](https://www.nature.com/articles/s41467-023-44373-x). Weaker direct evidence supports an explicitly translated SD profile from [Liu et al., Experimental Dermatology 2021](https://doi.org/10.1111/exd.14326). Neither establishes a routine clinical microbiome diagnostic test. The index measures directional alignment with selected associations, not disease probability, causation, severity, treatment response, or donor safety.

| Clinical category | Required behavior |
|---|---|
| CSU; historical chronic idiopathic urticaria terminology where used | Score CSU when analytically eligible; preserve each study's original phenotype |
| Symptomatic dermographism, a chronic inducible subtype | Separate SD index with permanent transport warning until independently validated |
| Coexisting CSU+SD | Show both eligible indices plus coexistence evidence; no third score or multiplication of indices |
| Cold, delayed-pressure, solar, heat, vibratory angioedema, cholinergic, contact, aquagenic urticaria | Each gets `evidence_not_established`, no subtype score or CSU substitution; no sufficient direct human stool-signature cohort was located by this audit |
| CSU with angioedema | Metadata modifier, not a separate stool classifier |
| Type I autoallergic/type IIb autoimmune CSU | Clinical/laboratory metadata only; never assigned from taxa |
| Acute urticaria, urticarial vasculitis, systemic mastocytosis, MCAS, hereditary/bradykinin-mediated angioedema | Not interchangeable with CU; no borrowed profile labels |

Use the [2026 international urticaria guideline](https://doi.org/10.1111/all.70210) for clinical classification context. “Not established” records a search result, not proof that no association exists. The [CIndU phenotyping study NCT07359430](https://clinicaltrials.gov/study/NCT07359430) is not a published validated stool signature. Wheat-dependent exercise-induced anaphylaxis and skin/sweat Malassezia mechanisms are not gut profiles for cholinergic urticaria.

## 2. Integration and input contracts

### 2.1 Repository inspection and components

Inspect the actual profile loader, feature namespace, scored profiler/database, reference loader, result schema, renderer, and test runner before editing. Historical reports used MetaPhlAn 3 and newer extended profiles; do not assume the scored lane is MetaPhlAn 4. Pin the actual lane for CU, or explicitly build a separate compatible lane without changing other scores.

Implement these logical components using repository-native paths:

| Component | Responsibility |
|---|---|
| CU evidence registry | Study and intervention records embedded below |
| CU profile compiler | Exact feature manifests, taxonomy mapping, version/digest |
| Reference adapter | Offline fitting, frozen statistics, population/assay compatibility |
| CU engine | Numerical contract and censoring/missingness |
| Result adapter | Typed computed/unavailable records |
| Report extension | Front-loaded cards, all feature rows, evidence/safety text |
| Acquisition/reproduction command | Digest-checked source download and optional raw reprocessing, never during inference |
| Tests | Numerical, provenance, integration, regression, and visual checks |

### 2.2 Metadata

Use numeric nulls, not strings, for unknown numbers. Minimum extension:

```json
{
  "age_years": null,
  "body_site": "stool",
  "urticaria": {
    "clinical_status": "unknown",
    "diagnosis_source": null,
    "duration_weeks": null,
    "spontaneous_wheals": null,
    "angioedema": null,
    "inducible_subtypes": [],
    "uas7": null,
    "uct": null,
    "antihistamine_response": "unknown",
    "omalizumab_exposure": "unknown",
    "omalizumab_response": "unknown",
    "autoimmune_endotype_confirmed": null
  },
  "exposures": {
    "antibiotic_last_date": null,
    "probiotic_last_date": null,
    "bowel_preparation_last_date": null,
    "acute_gastroenteritis_last_date": null,
    "diet_restriction": [],
    "current_medications": [],
    "history_complete": false
  }
}
```

`clinical_status`: clinician_diagnosed/self_reported/suspected/absent/unknown. Validate nonnegative age/duration; UAS7 integer0–42, UCT integer0–16; preserve measurement dates. Missing questionnaires do not block stool scoring. Symptoms persisting >6 weeks concern clinical chronicity, not index computation. Diagnosis/duration contradictions get a review flag rather than a silent rewrite.

Never use diagnosis, questionnaires, medication response, filename, SAMPLE2's known state, or donor-health descriptions as microbiome predictors. A healthy adult with duration0 remains eligible for research computation. Do not infer raw-study phenotype from an alias unless independently audited.

### 2.3 Measurement records

Each feature must include canonical ID, original database ID, rank/namespace, body site, assay, profiler/database digests, value, units/denominator, mapping/QC provenance, and one state:

`detected | below_reporting_limit | not_in_database | unresolved_mapping | assay_not_run | failed_qc`.

Include an experimentally established reporting limit when available. Zero does not represent missing data. A profiler lacking a species has not measured its absence. Do not repurpose host-filtered reads for host genetics or immune endotyping.

## 3. Embedded research registry

These are implementation seed records, not instructions to commission another research document. Preserve phenotype, assay, contrast, population, source, access status and limitations. Unverified metadata are explicit blocked fields, not guessed values. A paper can support context without qualifying as a training dataset.

| Study / primary source | Verified design/contribution | Implementation use and limits |
|---|---|---|
| `ZHU_2024` — [10.1038/s41467-023-44373-x](https://www.nature.com/articles/s41467-023-44373-x) | 26 adult CSU+26 HC human stool shotgun; separate plasma SCFA29+38 and LPS16+16 comparisons | Main taxonomy seed. Separate subsets; do not add counts into one paired cohort. Mouse interventions do not establish human benefit |
| `NABIZADEH_2017` — [10.1016/j.anai.2017.05.006](https://pubmed.ncbi.nlm.nih.gov/28668239/) | Targeted qPCR20+20; lower patient A. muciniphila; C. leptum/F. prausnitzii trends not conventionally significant | Context only; no WGS validation or confirmed deficiency from nonsignificant trends |
| `LU_2019` — [primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6881578/) | Small CSU16S10+10 cohort | Early association context, not classifier validation |
| `WANG_2020` — [10.3389/fcimb.2020.00024](https://doi.org/10.3389/fcimb.2020.00024) | Human16S+metabolomics | Association context; no serum concentration imputation from DNA |
| `WANG_2021` — [10.3389/fimmu.2021.691304](https://doi.org/10.3389/fimmu.2021.691304) | 39+40, metabolomics12+12; no clear alpha-diversity difference; Lactobacillus/Turicibacter/Lachnobacterium enrichment | Conflict evidence against universal low-diversity or low-Lactobacillus rules |
| `LIU_SD_2021` — [10.1111/exd.14326](https://doi.org/10.1111/exd.14326) | Pure SD22+22,16S; OEP001229 | Two curated SD candidates, not an entire inferred species/pathway list |
| `LIU_CSD_2021` — [10.3389/fcimb.2021.703126](https://doi.org/10.3389/fcimb.2021.703126) | CSU+SD25+25,16S+qPCR; OEP001891 | Separate coexistence phenotype; not pure-SD or pure-CSU validation |
| `SONG_2022` — [10.3389/fcimb.2022.831489](https://doi.org/10.3389/fcimb.2022.831489) | Resistant CSU25/responsive19/HC19,16S; PRJNA809140 | Response research only, no treatment selector |
| `LIU_RESPONSE_2022` — [10.1111/exd.14460](https://doi.org/10.1111/exd.14460) | Response discovery15+15, qPCR30+30; Lachnospira association | Treatment-specific evidence, not automatic resistance inference |
| `LUO_2023` — [10.3389/fcimb.2022.1094737](https://doi.org/10.3389/fcimb.2022.1094737) | CSU15+HC15,16S/metabolomics; PRJNA901136 | Diversity-direction conflict; not external WGS performance |
| `CESIC_2023` — [10.3390/life13061280](https://doi.org/10.3390/life13061280) | Croatian CSU22+HC23,16S; lower evenness/Lachnospiraceae, higher Lactobacillus | Conflicts retained. [Krišto Life13(1):152](https://doi.org/10.3390/life13010152) is a review, not another independent cohort |
| `WANG_OMALIZUMAB_2023` — [10.2147/CCID.S393406](https://doi.org/10.2147/CCID.S393406) | Adolescent longitudinal16S before/after treatment | Not pediatric disease-classifier validation |
| `PODDER_2025` — [10.1002/clt2.70027](https://doi.org/10.1002/clt2.70027) | Refractory CSU20+HC15,16S; four CSU also SD; treatment exposure including omalizumab; increased diversity | Retain treatment/subtype confounding and diversity disagreement |
| `PARK_2025` — [10.1111/all.16601](https://doi.org/10.1111/all.16601) | CU84+HC30,16S; plasma LPS/LL-37 and activity associations | Broad CU, not necessarily untreated pure CSU; data by request. No Firmicutes:Bacteroidetes health dial |
| `CIFTCI_2025` — [Pathogens14(11):1140](https://www.mdpi.com/2076-0817/14/11/1140) | Urticaria33+HC34; Blastocystis strata; full-length Nanopore16S despite loose metagenomic terminology; no significant diversity difference | Chronicity/CSU criterion not established; age means29 vs41 contradict matching claim. Exclude from clean CSU training pending clarification; no Blastocystis causality/eradication rule |
| `CHO_2026` — [10.1111/exd.70208](https://doi.org/10.1111/exd.70208) | 14 refractory CSU, longitudinal16S and omalizumab response; [OTU data](https://doi.org/10.6084/m9.figshare.30286705) | Small response study, metadata partly by request; no independent prediction. Industry relationships recorded |
| `WANG_LG1_2026` — [10.1111/1751-7915.70316](https://doi.org/10.1111/1751-7915.70316) | Human35 untreated CSU+21 HC,18–65;16S+blood metabolomics; Lactobacillus/Blautia/Bifidobacterium lower, Ruminococcus/Bacteroides higher, Proteobacteria lower; circulating urate/hypoxanthine higher | Human amplicon replication/conflict branch. LG-1 intervention only mice. Human, mouse, and isolate accessions must remain separate |
| `DENG_MR_2026` — [10.1016/j.jacig.2026.100686](https://doi.org/10.1016/j.jacig.2026.100686) | MR430 taxa/1400 metabolites/731 immune phenotypes/91 proteins; no associations FDR<0.05; no colocalization | Hypothesis only; no patient FASTQs, fitted weights, or causal confirmation |

### 3.1 Conflict and scope rules

- Diversity is lower, unchanged, or higher across cohorts. Neither disease index contains Shannon/evenness.
- Genus Bacteroides/Lactobacillus/Bifidobacterium and phylum Proteobacteria are inconsistent. Do not overwrite a species sign with a genus summary.
- Ruminococcus species have opposite signs. Never substitute genus abundance for R. gnavus, R. obeum, or R. bromii.
- Xiangya-associated OEP001229/OEP001891/OEP002960 may overlap. Audit participant/control identity before declaring independent replication.
- In the 2026 LG-1 study, paper region/primer descriptions and some archive metadata disagree; reconcile before reanalysis (§5).
- Shared diet, medication, inflammation, constipation, and other disease signals limit specificity. A review, MR, qPCR assay, isolate genome, and patient WGS are different evidence types and cannot be counted as equivalent replications.
- The 2025–cutoff update identified the five direct patient reports above (Podder, Park, Ciftci, Cho, Wang); no newer patient-shotgun CU cohort was located. This is a search limitation, not a guarantee of exhaustive global coverage.

## 4. Exact feature manifests

### 4.1 CSU canonical seed JSON

The following values were extracted from the publisher's source workbook, sheet `Fig 1b`, and directions checked against Figure1b. `coef` and `q_value` are provenance, **not model weights**. All computational weights equal1.

```json
{
  "profile_id":"CSU_STOOL_WGS_V1",
  "profile_revision":"4.1.1",
  "source_id":"ZHU_2024",
  "assay":"shotgun_metagenomics",
  "rank":"species",
  "features":[
    {"id":"Alistipes_onderdonkii","direction":-1,"coef":-2.469669346,"q_value":0.015384911,"n_nonzero":35},
    {"id":"Alistipes_putredinis","direction":-1,"coef":-2.459051270,"q_value":0.027690473,"n_nonzero":45},
    {"id":"Alistipes_shahii","direction":-1,"coef":-2.039284071,"q_value":0.027690473,"n_nonzero":41},
    {"id":"Bacteroides_cellulosilyticus","direction":-1,"coef":-1.642682855,"q_value":0.100549642,"n_nonzero":38},
    {"id":"Bacteroides_intestinalis","direction":-1,"coef":-1.888324131,"q_value":0.042168812,"n_nonzero":27},
    {"id":"Coprococcus_catus","direction":-1,"coef":-2.167381929,"q_value":0.002208388,"n_nonzero":34},
    {"id":"Odoribacter_splanchnicus","direction":-1,"coef":-1.955020950,"q_value":0.098787634,"n_nonzero":31},
    {"id":"Roseburia_hominis","direction":-1,"coef":-1.359086378,"q_value":0.128357101,"n_nonzero":40},
    {"id":"Ruminococcus_obeum","direction":-1,"coef":-0.831660425,"q_value":0.022195148,"n_nonzero":52},
    {"id":"Veillonella_parvula","direction":1,"coef":1.350313455,"q_value":0.047789302,"n_nonzero":31},
    {"id":"Bacteroides_stercoris","direction":1,"coef":1.827477087,"q_value":0.119496077,"n_nonzero":48},
    {"id":"Bilophila_wadsworthia","direction":-1,"coef":-1.083777733,"q_value":0.119496077,"n_nonzero":37},
    {"id":"Escherichia_coli","direction":1,"coef":1.389339140,"q_value":0.107643794,"n_nonzero":52},
    {"id":"Klebsiella_pneumoniae","direction":1,"coef":1.713724620,"q_value":0.124562782,"n_nonzero":31},
    {"id":"Ruminococcus_gnavus","direction":1,"coef":3.057857054,"q_value":0.002208388,"n_nonzero":40}
  ],
  "weight_rule":"equal",
  "sensitivity_subset_rule":"q_value < 0.05",
  "discovery_sample_count":52,
  "score_semantics":"unvalidated_directional_pattern_index"
}
```

Strict subset: Alistipes onderdonkii, A. putredinis, A. shahii, Bacteroides intestinalis, Coprococcus catus, Ruminococcus obeum, Veillonella parvula, Ruminococcus gnavus. Compile from q-values and assert exactly eight. E. coli/K. pneumoniae remain exploratory, not adjusted-significant simply because the mechanistic narrative emphasizes them. None of these detections establishes infection or a need for eradication.

### 4.2 Taxonomy resolution

Keep source names as immutable IDs. Resolve to production stable taxon/SGB IDs through a versioned explicit mapping. Historical synonyms are accepted only when the database establishes equivalent biological scope; preserve both names. No fuzzy genus substitution.

If one historical species maps to several SGBs, aggregate only a documented non-overlapping set. Never sum a parent total plus its children. An unresolved mapping lowers coverage, not abundance. An E. coli/Shigella complex is not automatically an exact E. coli measurement; require equivalent resolution in query and reference or exclude that slot.

### 4.3 SD canonical seed JSON

```json
{
  "profile_id":"SD_STOOL_TRANSLATED_V1",
  "profile_revision":"4.1.1",
  "source_id":"LIU_SD_2021",
  "source_assay":"16S_amplicon",
  "scoring_assay":"shotgun_metagenomics",
  "features":[
    {"id":"Subdoligranulum","rank":"genus","direction":-1},
    {"id":"Ruminococcus_bromii","rank":"species","direction":-1}
  ],
  "weight_rule":"equal",
  "score_semantics":"unvalidated_16S_to_shotgun_directional_pattern_index"
}
```

Both were enriched in the healthy-control core in [SD supplementary FigureS5](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111/exd.14326&file=exd14326-sup-0002-FigS1-S6.pdf), supporting depletion-oriented seeds. That core used ≥90% prevalence and mean abundance>0.3%; these discovery definitions are **not patient healthy cutoffs**. Only these curated candidates are imported; technically implausible species calls elsewhere in the 16S table are not copied into shotgun/pathogen rules.

Aggregate non-overlapping children for Subdoligranulum. Ensure the two mapped slots are not ancestor/descendant or otherwise overlapping. Broad Verrucomicrobia/Ruminococcaceae/Enterobacteriales findings are context only, not extra terms. Prevotella stercorea's duration association is not a case-control sign. CSD qPCR must not be falsely described as validation conducted in the pure-SD cohort.

SD FigureS6 is Tax4Fun-predicted **histidine metabolism**, not directly measured histamine production. Do not generate a histamine-gene module from the abstract's wording.

### 4.4 Separate functional observations

Exact values from source sheet `Fig 1h left panel`; `log2_ratio` is the group-mean CSU/control contrast. These are unadjusted `t_pvalue` entries, not q-values, individual cutoffs, or classifier coefficients.

| MetaCyc ID | Meaning | log2_ratio | p-value |
|---|---|---:|---:|
| `COLANSYN-PWY` | Colanic-acid building blocks | 1.3417623472 | 0.0277745409 |
| `ECASYN-PWY` | Enterobacterial common antigen | 1.9514882754 | 0.0069812917 |
| `ENTBACSYN-PWY` | Enterobactin biosynthesis | 1.5776287273 | 0.0125725900 |
| `LPSSYN-PWY` | LPS-biosynthesis superpathway | 1.5078246893 | 0.0255302956 |
| `PWY0-1338` | Polymyxin-resistance pathway annotation | 1.4994420163 | 0.0177212517 |
| `KDO-NAGLIPASYN-PWY` | (Kdo)2-lipid-A biosynthesis | 1.7960406518 | 0.0081928234 |

Display directly measured DNA-based potential and a compatible reference position if available. Missing assays stay unmeasured, not inferred from taxa. These overlapping pathways are not six independent mechanisms: never add parent/child or stratified/unstratified totals, and never sum them into the taxonomy index. A resistance pathway annotation alone is not a clinically validated antimicrobial-susceptibility result.

Existing validated SCFA/histamine/bile-acid/tryptophan modules may be linked as separate generic functions. If absent, label the function unmeasured; do not invent sequences, concentrations, or CSU-specific weights. Specific boundaries:

- Butyrate-producer taxa, production genes, fecal butyrate and plasma butyrate are distinct measurements. Zhu's acetate/propionate/caproate results were plasma measurements.
- R. gnavus is Gram-positive; its inflammatory polysaccharide evidence does not prove canonical LPS production. [Primary molecular study](https://doi.org/10.1073/pnas.1904099116).
- Histidine-decarboxylase capacity is strain/context dependent and does not diagnose histamine intolerance, DAO deficiency, MCAS or CSU. [Mou et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8465708/) is not CSU diagnostic validation.
- The 2026 urate/hypoxanthine results are circulating metabolomics, not concentrations measurable in a FASTQ-only report.

## 5. Reproducible acquisition and reference data

### 5.1 Publisher source artifact: verified download

Download the [Zhu source-data ZIP](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-44373-x/MediaObjects/41467_2023_44373_MOESM4_ESM.zip). It was downloaded and numerically inspected for this revision. This public research dataset is not a separately required handoff document.

```json
{
  "source_id":"ZHU_2024_SOURCE_DATA",
  "zip_sha256":"a1d4ce95eb94caa6a05eaf03f4b27a8bf907e4ff242d2180adb118a728497c24",
  "workbook":"Source data-Fig1.xlsx",
  "workbook_sha256":"e660e1f33dc5ca2a86d6c1c9ef0b7f1912ad0335b25e7f6592b456be8c0bfde3",
  "feature_sheet":"Fig 1b",
  "feature_header_row_1_based":2,
  "feature_rows_1_based":[3,17],
  "feature_columns":["Species","coef","stderr","N","N.not.0","p-value","q-value"],
  "matrix_sheet":"Fig1a left panel",
  "matrix_header_row_1_based":2,
  "matrix_sample_columns_1_based":[2,53],
  "label_sheet":"Fig 1c",
  "label_header_row_1_based":2,
  "label_columns_1_based":[1,2],
  "pathway_sheet":"Fig 1h left panel",
  "pathway_header_row_1_based":2
}
```

Validate HTTP status, ZIP integrity, safe archive paths, and digest. Changed content becomes `source_revision_requires_review`; do not silently overwrite the locked source. Store URL, retrieval date, digests, parser version, and licensing/data-use status.

Verified parsing constraints:

1. `Fig 1b` contains15 rows, all N=52. Compare coefficient signs and q-values with §4.1 within1e-9.
2. The matrix contains **308 rows mixing species and strains**,52 samples, and an extra statistics column. Retain terminal `s__` rows only for species analysis, not terminal `t__` rows.
3. This yields **190 species**, with sample-column sums **73.25779–98.79070**, not100. It is a filtered published matrix, not a complete composition. Do not silently reclose it, use it as a full CLR denominator, or mix it with production references.
4. The label sheet supplies explicit sample IDs and26 HC/26 CSU labels. Join by ID, never by column order: source sheets differ in their ordering.
5. Matrix abundances are percentage-scale; some plotted sums are fractions. Read per-sheet units rather than blindly rescaling all sheets.
6. Some non-species figure/workbook p-values differ, including plasma LPS. Preserve their source-specific values; do not turn them into patient cutoffs. Verified species signs agree with Figure1b.
7. Relapse data cover22 participants. The exploratory survival analysis is not a deployable, independently validated relapse model.

The source artifact supports exact provenance/parser reconstruction. It does not establish production database equivalence or independent diagnostic performance.

### 5.2 Accession routing

| Record | Material and routing | Access/metadata boundary |
|---|---|---|
| [OEP002960](https://www.biosino.org/node/project/detail/OEP002960); [publication review link](http://www.biosino.org/node/review/detail/OEV000435?code=K4MN2WY7) | Human CSU/HC stool shotgun; primary reprocessing target | Publication-declared; raw download and participant join not completed in this audit |
| [OEP001229](https://www.biosino.org/node/project/detail/OEP001229) | Pure SD/HC16S | Publication-declared; no shotgun substitution |
| [OEP001891](https://www.biosino.org/node/project/detail/OEP001891) | CSU+SD/HC16S | Separate phenotype; raw join must be verified |
| [PRJNA809140](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA809140) | CSU antihistamine-response16S | Response research branch only |
| [PRJNA901136](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA901136) | CSU/HC16S | Cross-cohort direction audit, not WGS validation |
| [PRJNA1096288](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1096288) | 56 human AMPLICON runs,2026 study | Run metadata verified; body-site/phenotype/date inconsistencies must be reconciled before importing |
| [PRJNA1101664](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1101664) | 30 mouse-gut AMPLICON runs | Animal research only; excluded from human reference |
| [PRJNA1359023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1359023) | LG-1 isolate genome | Not patient stool metagenomics |
| [OEP002997](https://www.biosino.org/node/project/detail/OEP002997) | Zhu experimental mouse16S | Animal mechanism only |
| [CNP0005021](https://db.cngb.org/search/project/CNP0005021/) | Zhu mass spectrometry | External metabolomics only; exact downloadable matrix/join not verified |
| [Figshare30286705](https://doi.org/10.6084/m9.figshare.30286705) | Cho omalizumab OTU data | Not proof of public FASTQs, clinical labels, or a pretrained model |

For PRJNA1096288, [SAMN40750031](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN40750031) reports `isolate=feces`, `tissue=plasma`, and June12 collection while the paper describes March–May recruitment. Alias `FDMP_H32` does not establish healthy-control status. The paper's V3–V4 wording and515F/806R primers also require reconciliation. Block a phenotype/body-site-dependent import until resolved; do not make assumptions from a matching run count.

Use distinct access states: `download_verified`, `metadata_verified`, `publication_declared`, `author_request_required`, `access_blocked`, `not_reported`. Track license and participant mapping separately. “Public” is not sufficient evidence of usable labelled FASTQs.

An acquisition command can query the official [ENA run-report endpoint](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1096288&result=read_run&fields=run_accession,sample_accession,scientific_name,library_strategy,library_source,fastq_ftp,fastq_md5&format=tsv), changing accession for each study. Validate organism, library strategy, body site, paired-end layout, checksums and participant join. `library_source=GENOMIC` does not override `library_strategy=AMPLICON`. One run is not necessarily one participant. Respect access/data-use conditions; blockers are explicit states, not fabricated data or permission workarounds.

### 5.3 Uniform reprocessing

Zhu reports KneadData0.10.0, Trimmomatic0.33, Bowtie2 2.2 and HUMAnN2. These names do not establish its exact taxonomic database. Do not declare equivalence to the current MetaPhlAn lane.

For production reference/model work, reprocess permitted human controls/cases under identical pinned QC, host filtering, profiler/database, normalization, denominator, and detection policy. Preserve container/database digests and all parameters, read counts, QC and host-removal status. Audit extraction/platform/geography differences. Do not train or calibrate on the current five personal samples; do not treat healthy minors as adult references.

## 6. Exact scoring contract

### 6.1 Meaning and reference

This is an **unvalidated directional pattern index**,0–100. It is not a probability, percentage match, default percentile, or geometric distance to a disease centroid. Fifty means neutral signed deviation in this transformation, not a diagnostic threshold or definition of health. Higher values mean stronger alignment with selected study directions. Every constant below is an engineering choice, not a clinically derived cutoff.

Use a frozen bundle containing ID/digest, source/participant provenance, code/profile digests, body site, assay, profiler/database, rank mappings, abundance denominator/units, detection policy, population/exposure eligibility, and per-feature log median/scale/count/detection prevalence. Set `reference_scope=compatible_independent_controls` or `discovery_controls_only` accurately. Separate any held-out calibration subjects from fitting and feature discovery.

Require ≥20 distinct eligible adult control participants per fitted feature as an engineering minimum. Use the prespecified first eligible baseline sample per participant; technical replicates do not raise N. Apply the same measurement-state and reporting-limit policy to controls as queries before fitting: an unresolved taxon or unquantifiable nondetection is not an observed zero. Twenty controls does not confer clinical validity. Existing compatible production controls may support the directional index without waiting for CU author data, but not a disease-specific validation claim. The filtered published matrix alone cannot become a production reference without established analytic equivalence.

Default supported population: adults18–65, intersected with the bundle's actual supported age range. Missing/out-of-range age blocks participant-facing numeric indices; measurements and evidence remain visible. A pediatric intervention study does not validate this adult classifier. Reference unavailable means an explicit data blocker, not a fabricated “healthy” baseline.

### 6.2 Transformation and aggregation

Normalize only by declared units: percentage/100 becomes fraction. Require finite a∈[0,1]. Use the reference's total-community denominator; never reclose only the selected CU panel. Handle unclassified/non-bacterial mass identically in controls and queries.

For controls i and feature j:

`x_ij = log10(a_ij + 0.000001)`

`mu_j = median_i(x_ij)`

`sigma_j = max(1.4826 * median_i(abs(x_ij - mu_j)), 0.25)`

For the query:

`z_j = clip((log10(a_j + 0.000001) - mu_j) / sigma_j, -6, 6)`

`contribution_j = 50 + 50 * tanh(direction_j * z_j / 2)`

Index = arithmetic mean of usable contributions. Equal weights only: no q-value, publication-count or mechanistic-importance weights. Fit/serialize mu and sigma offline on controls only; never refit on inference batches. This explicit log-relative-abundance algorithm replaces the prior incomplete CLR formula and is not presented as the published classifier.

Compute CSU full15 and strict8 results independently under the same transformation. The latter is sensitivity analysis, not a second independent test. SD uses its two candidates and always carries `16S_to_shotgun_unvalidated`.

### 6.3 Missingness, censoring and comparisons

Detected measurement: use its value; propagate a validated quantitative interval if available. Without such an interval, the calculation uses a point, but must not imply measurement is error-free.

Below a validated reporting limit L: computational point a=0, censoring interval[0,L]. Transform both endpoints and reorder after applying direction. This means “not detected at this reporting limit,” not biological absence. If L is unavailable, treat the feature as unquantifiable for the index: omit its point and preserve missingness bounds.

Missing/not-in-database/unresolved/failed feature: no point contribution; possible contribution[0,100]. With locked panel size K and usable set M:

`coverage = |M|/K`

`index = sum(point_j for j in M)/|M|` when M is nonempty; otherwise null.

`full_panel_lower = sum(lower_j for j in M)/K`

`full_panel_upper = (sum(upper_j for j in M) + 100*(K-|M|))/K`

Zero usable features yields null with bounds[0,100]. One usable feature may produce a **limited-feature directional index**, explicitly labelled1/K, with the full-panel bounds. Do not silently promote it to a full-panel result. Compare different dates/samples only with the same reference and feature mask, or explicitly recompute on a named common mask.

These are **coverage/censoring bounds, not95% confidence intervals**. They omit unknown sampling/population transport and disease-specificity uncertainty. Evidence maturity, analytical coverage, reference compatibility and disease specificity are separate fields; score magnitude is not confidence. Do not display zero-width propagated bounds as proof of certainty.

### 6.4 Optional empirical percentiles

Default release does not need percentiles. If genuinely held-out compatible scores exist, report control and case percentiles separately using:

`P(s)=100 * sum_i(w_i * (I(s_i<s)+0.5*I(s_i==s))) / sum_i(w_i)`

Equal participant weights within cohort; equal total weights across cohorts. Store cohort IDs, weights, score revision and feature mask. Recompute the reference distribution for the query mask or leave percentile null. Never mix case/control scores into one unnamed population. An in-sample rank is not held-out; feature selection on the original52 subjects prevents an external-validation claim. A case percentile is not posterior disease probability. Do not print a percentage sign next to the uncalibrated index.

### 6.5 QC/domain gates

Reuse a validated existing QC gate with its provenance. A new lane without established analytical QC is `analytic_validation_required`, not approved using an invented universal read-count threshold. Non-stool/non-shotgun input is incompatible.

Gate order: structural validity → body site/assay → QC → age/reference population → reference analytical compatibility → feature usability → calculation. Preserve all material reasons. Exposures produce specific flags and assessment against reference inclusion criteria; unknown history is not clean. Without a fitted/validated OOD model use `domain_assessment=limited_metadata_based`, never assert proven in-distribution status. Clipping |z| at6 is numerical stabilization, not a biological normal range or OOD test.

## 7. Executable reference kernel

Port this dependency-free behavior into the repository language. Eligibility/provenance checks above run before the kernel. `observations` contains only accepted feature intervals; unusable observations are omitted. Reference fitting runs offline.

```python
import math
from statistics import median

EPS = 1e-6
SCALE_FLOOR = 0.25
Z_CLIP = 6.0
TEMPERATURE = 2.0

def fraction(value):
    if isinstance(value, bool):
        raise ValueError("boolean is not abundance")
    value = float(value)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("expected finite fraction in [0,1]")
    return value

def fit_reference(rows, feature_ids):
    if len(rows) < 20:
        raise ValueError("at least 20 distinct eligible controls required")
    out = {}
    for fid in feature_ids:
        vals = [math.log10(fraction(r[fid]) + EPS)
                for r in rows if fid in r and r[fid] is not None]
        if len(vals) < 20:
            continue
        mu = median(vals)
        mad = median(abs(v - mu) for v in vals)
        out[fid] = {"mu": mu, "scale": max(1.4826 * mad, SCALE_FLOOR),
                    "n_controls": len(vals)}
    return out

def z_value(abundance, ref):
    mu, scale = float(ref["mu"]), float(ref["scale"])
    if not math.isfinite(mu) or not math.isfinite(scale) or scale < SCALE_FLOOR:
        raise ValueError("invalid frozen reference")
    z = (math.log10(fraction(abundance) + EPS) - mu) / scale
    return max(-Z_CLIP, min(Z_CLIP, z))

def support(z, direction):
    if isinstance(direction, bool) or direction not in (-1, 1):
        raise ValueError("direction must be -1 or +1")
    if not math.isfinite(z):
        raise ValueError("nonfinite standardized feature")
    return 50.0 + 50.0 * math.tanh(direction * z / TEMPERATURE)

def summarize(contributions, panel_ids):
    if not panel_ids or len(set(panel_ids)) != len(panel_ids):
        raise ValueError("panel requires unique IDs")
    usable = [contributions[f] for f in panel_ids if f in contributions]
    k, m = len(panel_ids), len(usable)
    point = sum(v[0] for v in usable) / m if m else None
    lower = sum(v[1] for v in usable) / k
    upper = (sum(v[2] for v in usable) + 100.0 * (k - m)) / k
    return {"index": point, "usable_count": m, "panel_count": k,
            "coverage": m / k, "full_panel_bounds": [lower, upper],
            "mask": [f for f in panel_ids if f in contributions]}

def evaluate(features, observations, reference):
    ids = [f["id"] for f in features]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate feature")
    contributions = {}
    for f in features:
        fid = f["id"]
        if fid not in observations or fid not in reference:
            continue
        o = observations[fid]
        a, lo, hi = map(fraction, (o["point"], o["lower"], o["upper"]))
        if not lo <= a <= hi:
            raise ValueError("point must lie in interval")
        p = support(z_value(a, reference[fid]), f["direction"])
        ends = [support(z_value(v, reference[fid]), f["direction"])
                for v in (lo, hi)]
        contributions[fid] = (p, min(ends), max(ends))
    result = summarize(contributions, ids)
    result["contributions"] = contributions
    return result

def weighted_midrank(score, scores, weights):
    if not scores or len(scores) != len(weights):
        raise ValueError("nonempty equal-length arrays required")
    if not math.isfinite(score) or any(not math.isfinite(s) for s in scores):
        raise ValueError("nonfinite score")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("positive finite weights required")
    return 100.0 * sum(w * (1.0 if s < score else 0.5 if s == score else 0.0)
                       for s, w in zip(scores, weights)) / sum(weights)

def self_test():
    assert support(0, 1) == 50
    assert abs(support(2, 1) - 88.07970779778825) < 1e-10
    assert abs(support(2, -1) - 11.92029220221176) < 1e-10
    r = summarize({"a": (80.0, 70.0, 90.0)}, ["a", "b"])
    assert r["index"] == 80 and r["coverage"] == 0.5
    assert r["full_panel_bounds"] == [35.0, 95.0]
    assert summarize({}, ["a"])["index"] is None
    assert weighted_midrank(50, [20, 50, 50, 80], [1, 1, 1, 1]) == 50
    refs = fit_reference([{"a": 0.01}] * 20, ["a"])
    assert refs["a"]["scale"] == SCALE_FLOOR
    f = [{"id": "a", "direction": -1}]
    obs = {"a": {"point": 0.0, "lower": 0.0, "upper": 0.001}}
    r = evaluate(f, obs, refs)
    assert r["full_panel_bounds"][0] <= r["index"] <= r["full_panel_bounds"][1]
    assert r == evaluate(f, obs, refs)
    for bad in (True, -0.1, 1.1, float("nan"), float("inf")):
        try:
            fraction(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid abundance accepted")
    print("CU numerical self-test passed")

if __name__ == "__main__":
    self_test()
```

Repeated rows in `self_test` are synthetic numerical fixtures only. Production must enforce distinct participant identities before fitting. Never use these fixtures as patient reference data. Validate schema types before passing data to this kernel.

## 8. Result schema and states

Emit a record for every registered task, including unavailable tasks. Example unavailable CSU result:

```json
{
  "profile_id":"CSU_STOOL_WGS_V1",
  "profile_revision":"4.1.1",
  "status":"reference_unavailable",
  "score_semantics":"unvalidated_directional_pattern_index",
  "score_0_100":null,
  "coverage":0.0,
  "usable_count":0,
  "panel_count":15,
  "full_panel_bounds":[0.0,100.0],
  "bounds_type":"missingness_and_censoring_not_confidence_interval",
  "strict_subset_result":null,
  "control_percentile":null,
  "case_percentile":null,
  "reference_bundle_id":null,
  "reference_scope":null,
  "feature_mask":[],
  "assay_transport":"human_shotgun_discovery_to_production_requires_compatibility",
  "evidence_maturity":"exploratory_single_center_discovery",
  "disease_specificity":"not_established",
  "clinical_status":"unknown",
  "reason_codes":["no_compatible_frozen_reference"],
  "confounder_flags":[],
  "source_ids":["ZHU_2024"],
  "observations":[],
  "functional_observations":[],
  "intervention_evidence_ids":[],
  "treatment_recommendation":null,
  "donor_eligibility":null
}
```

For eligible data set `computed` or `computed_partial` and fill real values. Unavailable states: `no_usable_features`, `qc_failed`, `assay_incompatible`, `population_not_supported`, `reference_unavailable`, `analytic_validation_required`. A null/failed result never becomes a numeric zero. Preserve all material reasons; the first failed gate supplies the primary status.

Use separate clinical-observation records for diagnoses, symptoms and external labs, with dates and sources. Scores can retrieve educational evidence but cannot create confirmed diagnoses or prescriptions. Unsupported subtype entries use `evidence_not_established` and are not additional score-bearing tasks.

## 9. Report and visual implementation

Place CU/SD in the complete at-a-glance disease-profile section near the beginning, after the existing summary/“what stood out” content and before explanations. Do not impose an eight-profile or one-page cap. Preserve all old cards and later detail sections.

Each calculated card shows the0–100 **index without a percent sign**, “research pattern—not a diagnosis,” a neutral directional bar labelled less/more alignment, coverage N/K, full-panel bounds, reference identity and evidence maturity. CSU also shows its strict-eight sensitivity index and coverage. SD shows the16S-to-shotgun warning. No green healthy/red diseased zones or clinical score thresholds.

Show the numerical full-versus-strict difference without turning it into a validated disagreement threshold. The analyses share features and are not independent tests. Missing/partial masks must remain visible. Accessible text labels supplement color; bars and bounds must survive PDF/narrow layouts. Propagated point-only bounds must not imply confidence.

### 9.1 Complete detail pages

Every panel feature, including unavailable ones, gets a row with: source/current taxon name, study direction, measured abundance and units, measurement state, compatible-reference comparison, score contribution, source q-value if applicable, and limitations. Distinguish detection, reference-relative high/low, nondetection with a limit, unmeasured/unrepresentable, and reference-ineligible.

Never label a taxon “missing that should be there” merely because it appears in this panel. Existing healthy-reference logic may show commonly detected-but-not-detected only with matched reference prevalence and reporting limit. Relative depletion is not deficiency. Do not label E. coli/Klebsiella abundance as infection or an eradication target.

Follow the feature table with separate functional observations, conflicting studies, coexistence/response context, applicable intervention-evidence cards and direct citations. State whether each intervention endpoint was symptoms, quality of life, a lab marker, microbiome composition or only an animal mechanism. No “solution” is selected merely because it sounds opposite to a taxonomic finding.

### 9.2 Research screening/FMT boundary

An access-controlled FMT research prescreen may show nonspecific patterns, but never a donor pass/rank/clearance, pathogen-transmission exclusion, or reassurance based on a low CU index. Adult profiles do not qualify minors as donors. CU/FMT remains an investigational question subject to clinical and applicable regulatory oversight; this build does not issue an FMT protocol. No microbiome score authorizes antimicrobial eradication, FMT or stopping current urticaria treatment.

## 10. Intervention evidence: exact seed cards

This layer answers **what has been studied**, not what an individual should take. The user's antimicrobial experience is context, not an efficacy study or a treatment template.

### 10.1 Common card contract

Store `evidence_id`, source/registration, design, human/animal status, population, exact intervention identity, comparator, duration, endpoint/timepoint, between-group result, uncertainty, adverse-event limitations, source conflicts and trigger type. All current CU cards have:

```json
{
  "microbiome_selection_biomarker_validated":false,
  "microbiome_score_can_select_intervention":false,
  "microbiome_normalization_as_response_mechanism":"not_established",
  "protocol_details_are_recommendations":false,
  "participant_action":"none",
  "purpose":"evidence_for_discussion"
}
```

Allowed triggers: `profile_evidence_interest`, `clinician_confirmed_diagnosis`, `documented_symptom_context`, `external_lab_confirmed`, `general_education`. Profile interest permits educational retrieval only. Prescription/infection/deficiency discussion requires corresponding clinical or lab context and clinician review; the software never diagnoses it.

Exact probiotic identity matters. Missing strain IDs remain `not_reported`, not invented. A named mixture is not interchangeable with another brand or any organism of the same genus. Changes in fecal Lactobacillus/Bifidobacterium are not a validated treatment-success surrogate.

Any doses below are **study metadata for the internal evidence registry**, not regimens to recommend. Participant reports omit dose/procedure instructions, antimicrobial stacks and high-dose protocols. Unknown safety fields say not established; absence of reported adverse events is not proof of safety. Clinician-facing evidence may show source protocols with an explicit study-only label, never a generated prescription.

### 10.2 Probiotic/synbiotic records

| Evidence ID / source | Exact design/identity | Outcome and required interpretation |
|---|---|---|
| `CU_PROBIOTIC_NETTIS_2016` — [PMID27608474](https://pubmed.ncbi.nlm.nih.gov/27608474/), [primary report](https://www.eurannallergyimm.com/wp-content/uploads/2016/09/volume-probiotics-refractory-chronic-spontaneous-urticaria-1173allasp1.pdf) | Uncontrolled refractory-CSU adjunct;52 enrolled,38 evaluated. Bifiderm: L. salivarius LS01+B. breve BR03, each≥10^9 CFU/g; twice daily8weeks | Completers:27 no improvement,9 mild,1 significant,1 complete. Fourteen withdrawals,11 reporting no improvement. `mostly_nonresponse`, high attrition bias; do not highlight the complete responder while omitting nonresponders |
| `CU_PROBIOTIC_BI_2021` — [10.1186/s13223-021-00544-3](https://doi.org/10.1186/s13223-021-00544-3) | Pediatric CU age6–12;213 randomized,206 analyzed104/102;4-week placebo-controlled adjunct to desloratadine. Yimingjia: L. gasseri LK00140%, L. salivarius LK00220%, L. johnsonii LK00315%, L. paracasei LK0045%, L. reuteri LK0055%, **B. animalis LK01115%**. Product5×10^9CFU/g,1.5g twice daily. NCT03328897 | Response counts84/104 vs63/102; article percentage reporting differs. Wheal-size/frequency endpoints favored intervention; several other symptom endpoints did not. Not B. breve, not adult efficacy, and not pediatric stool-profile validation |
| `CU_PROBIOTIC_ATEFI_2022` — [10.23822/EurAnnACI.1764-1489.200](https://www.eurannallergyimm.com/wp-content/uploads/2022/05/volume-probiotic-adjuvant-therapy-chronic-urticaria-4842allasp1.pdf) | Adults18–45;52 enrolled,42 completed21/arm. LactoCare twice daily8weeks+antihistamines vs antihistamines; analyst blinded, no placebo/patient blinding. L. rhamnosus/L. casei/L. acidophilus/B. breve/L. bulgaricus/B. longum/S. thermophilus+FOS; strain IDs/CFU not verified. IRCT20190825044613N1 | Week8 UAS7 11±11.41 vs16.86±13.54; between-group activity advantage not significant. Secondary DLQI favored treatment. Code activity `null_result`, QoL `favorable_secondary`, high bias; not a positive placebo-controlled UAS7 trial |
| `CU_PROBIOTIC_DABAGHZADEH_2023` — [10.22088/cjim.14.2.192](https://caspjim.com/article-1-2885-en.pdf) | Blinded placebo-controlled trial; 38 analyzed (20/18), mixed age cohort. Femilact/ZistTakhmir twice daily for 8 weeks+cetirizine vs placebo+cetirizine; same seven listed species+FOS, 10^9 CFU formulation reported without strain-specific quantities. IRCT20110531006660N9 | Week 8 UAS7 9.6±6.4 vs12.7±8.1. P=.036 is reported main between-arm effect across measurements, not a verified standalone final-visit contrast. QoL P=.805 null. Conclusion conflicts with favorable abstract/main-effect wording; preserve `source_internal_inconsistency=true` |
| `CU_PROBIOTIC_CAN_2024` — [10.4103/ds.DS-D-24-00136](https://doi.org/10.4103/ds.DS-D-24-00136), [institutional study record](https://acikerisim.bahcesehir.edu.tr/server/oai/request?identifier=oai%3Aacikerisim.bahcesehir.edu.tr%3A20.500.14719%2F19401&metadataPrefix=oai_dc&verb=GetRecord) | **Nonrandomized** prospective30+30;4weeks. Limosilactobacillus reuteri **ATCC55730**,10^8CFU+ebastine vs ebastine | Response25/30 vs15/30,P=.006; changes in UAS7/QoL favored intervention. `favorable_nonrandomized`; do not call RCT or substitute DSM17938/unspecified L. reuteri |
| `CU_PROBIOTIC_GODSE_2025` — [10.18231/j.ijced.2025.015](https://doi.org/10.18231/j.ijced.2025.015) | Adult double-blind placebo-controlled97 randomized/87 analyzed;12weeks Lactogut twice daily+antihistamine. Combined5×10^9CFU/capsule: B. coagulans Unique IS-2, L. rhamnosus UBLR-58, B. longum UBBL-64, B. bifidum UBBB-55, **S. boulardii Unique-28 yeast**, S. thermophilus UBST-50; FOS20mg+lactitol10mg. CTRI/2023/06/054428 | Abstract favorable itch claim conflicts with body: between-group itch differences at baseline, not follow-up; UCT advantage absent. `null_or_uncertain`, baseline imbalance, source conflict, manufacturer sponsorship. Not an all-bacterial or interchangeable generic probiotic |

No pooled “probiotic cure rate” is permitted from these heterogeneous populations, products, endpoints, designs and conflicts. Treatment cards are retrieved as a balanced group, including null/mostly nonresponse studies.

### 10.3 Diet, nutrients, infection and mechanistic records

| Evidence ID / primary source | Evidence | Permitted behavior |
|---|---|---|
| `CU_DIET_WAGNER_2017` — [10.1111/jdv.13966](https://pubmed.ncbi.nlm.nih.gov/27624921/) | Uncontrolled56 CSU patients with GI symptoms;≥3-week low-histamine diet;34/56 met UAS4-improvement endpoint; DAO unchanged | Limited human symptom evidence, not microbiome correction or DAO-deficiency proof. Clinical-history-based, dietitian-supported time-limited assessment, not universal restriction |
| `CU_DIET_SON_2018` — [10.5021/ad.2018.30.2.164](https://pmc.ncbi.nlm.nih.gov/articles/PMC5839887/) | Uncontrolled22-adult,4-week dietary study; symptom improvement | Small uncontrolled evidence; no randomized-efficacy, permanent avoidance or individual response claim |
| `CU_VITD_MONY_2020` — [primary trial](https://pubmed.ncbi.nlm.nih.gov/31926152/) | Randomized placebo-controlled120 vitamin-D-deficient CU patients | Discussion requires external25-OH-vitamin-D assessment/clinical review; no deficiency inference or high-dose regimen from metagenomics |
| `CU_HPYLORI_PAWLOWICZ_2018` — [primary placebo-controlled study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5949544/) | Reported temporary benefit from infection eradication as add-on therapy | Require clinically appropriate infection confirmation; stool taxon abundance is not a treatment indication or guaranteed cure |
| `CU_HPYLORI_VALSECCHI_1998` — [primary negative study](https://medicaljournalssweden.se/actadv/article/view/14676) | Eradication did not clearly improve urticaria course | Retrieve alongside positive findings; no selective efficacy narrative |
| `CU_FIBER_FERMENTED_GENERAL` | No validated CU-specific response-selection rule established here for generic fiber/prebiotics/fermented foods | General diet context only. Do not automatically recommend fermented foods to someone reporting food/histamine-related reactions |
| `CU_ANTIMICROBIAL_STACK_UNSUPPORTED` | Adequate CU-specific human evidence not established in this audit for taxon-selected berberine/allicin/oregano/black-seed-oil/commercial multi-herbal stacks | No automatic recommendation, eradication protocol or promised restoration of commensals |
| `CU_SCFA_LG1_PRECLINICAL` — [Zhu](https://www.nature.com/articles/s41467-023-44373-x), [Wang2026](https://doi.org/10.1111/1751-7915.70316) | SCFA/bacterial mechanistic interventions and mouse LG-1 treatment | Preclinical-only card; no human-dose extrapolation or symptom-response prediction |

The [2026 guideline](https://doi.org/10.1111/all.70210) supports a separate clinician-care context, not a stool-profile-directed treatment algorithm. It may name second-generation H1 antihistamines and specialist escalation. For omalizumab/dupilumab/remibrutinib or other agents, current local indication/age/label review is required; this document does not promise universal approval or generate an individualized drug ladder. Autoimmune endotypes and treatment resistance require clinical assessment.

### 10.4 FMT: include the real report, exclude false matches

`CU_FMT_WU_2023`: [Wu et al., Faecal microbiota transplantation for treatment of chronic urticaria with recurrent abdominal pain and food allergy](https://pubmed.ncbi.nlm.nih.gov/37077050/), [10.4103/singaporemedj.SMJ-2021-423](https://doi.org/10.4103/singaporemedj.SMJ-2021-423). Bibliographic identity is verified; this is an uncontrolled case report, not proof of CSU efficacy. Full protocol was not accessible in this audit: age, dose, donor, delivery and durability remain unverified and must not be populated from secondary summaries. `research_only=true`, `participant_action=none`.

No controlled human FMT study establishing CU/CSU benefit was identified in this search. Do not say no human report exists; do not say FMT is the only solution. The [2026 report of allergic reactions after FMT in two children with ASD](https://doi.org/10.3389/fped.2026.1847568) supplies adverse-event context, not an incidence estimate or proof of the authors' proposed immune mechanisms. Do not use it to validate food-IgG testing.

Exclude [CCID.S443542](https://pmc.ncbi.nlm.nih.gov/articles/PMC10826708/) from CU efficacy: it concerns generalized eczema after vaccination. Atopic-dermatitis/peanut-allergy FMT studies and a patient's incidental CU comorbidity without CU-specific outcomes are not independent CU efficacy confirmations. A [late-August2026 FMT review](https://pubmed.ncbi.nlm.nih.gov/42618999/) is not itself a treatment trial.

### 10.5 Quarantined/unverified intervention records

- [Bux/Laique,10.53350/pjmhs02023171138](https://pjmhsonline.com/pjmhs/article/download/5372/5240/5897): methods specify age6–12 while results describe means near38; product identity/CFU missing. No actionable formulation record or pooled efficacy weight pending resolution.
- [Sharifi/Atefi/Fallahpour2025,10.61882/smmr.202501.08.03](https://www.simmr.info/article_234055.html): discovered but full methods/product/registration and independence not verified. Not a second independent positive replication by default.
- Secondary-only cold-urticaria probiotic citations do not supply a gut-signature manifest. Without verified primary design/product/outcomes, no quantitative treatment or profile seed.

Store quarantine reasons so a future curator can resolve them; do not silently erase contrary or incomplete research, and do not pretend it is usable evidence.

## 11. Implementation order, validation and definition of done

1. Capture installed task IDs and existing score/report regression fixtures. Preserve user edits. Resolve actual module paths rather than assuming a repository layout.
2. Register the two tasks, study/card records, subtype map and typed unavailable states. If upgrading an old CU implementation, migrate its IDs and invalidate incompatible CU-only caches.
3. Compile exact feature directions/statistics and explicit taxonomy mappings, digest all inputs, and implement the offline frozen-reference adapter.
4. Implement §6/§7, analytic/population gates, missingness/censoring and strict-subset comparison. Connect genuine existing compatible controls if available; otherwise implement the explicit data-blocked path.
5. Implement the digest-checked source parser and optional raw acquisition/reprocessing. Unit tests use cached fixtures; production inference does not download studies or fit models.
6. Connect measured functional/external-assay observations without adding them to the taxonomy index. Add balanced evidence retrieval and participant/clinician audience filters.
7. Render front-loaded cards and complete detail sections; preserve every existing profile and all previous report content.
8. Run the embedded numerical test, source-statistics reconciliation, repository unit/integration/regression tests, result-schema validation and visual checks.
9. With a compatible real reference, generate a real computed research result. Without one, produce a synthetic test fixture explicitly labelled synthetic plus a real `reference_unavailable` report. Report the exact blocker; do not claim production numeric activation or model validation merely because code passes tests.

### 11.1 Future learned-model branch, not a release shortcut

A learned classifier requires participant-grouped splits, feature selection and preprocessing fitted inside training only, nested tuning, serialized transformations and independent same-assay validation. Audit overlap across recruiting centers/cohorts. Include other inflammatory/dermatologic/GI diseases, not just healthy controls, to evaluate specificity. Preserve medication/diet/age/geography strata and report uncertainty and coverage.

The15-species panel was selected using the original52 subjects. Evaluation on those same subjects is discovery/reconstruction even with a later cross-validation wrapper. Do not quote internal AUC as expected external or clinical performance. Do not use SAMPLE2's known label or the other four personal samples to tune features, thresholds or weights.

SD needs independent shotgun validation before removing its transport warning. Antihistamine/omalizumab response needs a separate treatment-response validation task; existing association cards cannot be silently promoted to predictors. A future probabilistic classifier must use a new score-semantic type and its own calibration/validation rather than relabelling this index.

No validated CU-specific stool virome, mycobiome, MAG, host-variant, or DNA-imputed-metabolite diagnostic panel was established here. Keep the generic observation architecture extensible, but do not invent seed sequences or coefficients to fill these modalities.

### 11.2 Handoff from the coding agent

Report implemented modules, inventory before/after, profile/source/reference digests, passed/failed tests, real reference availability, and access/data blockers. Distinguish software complete, exploratory numeric profile activated, and independent clinical validation not established. No additional supplemental research handoff is necessary; ordinary application code, tests and machine-readable runtime artifacts are expected.

## 12. Acceptance tests

Implement deterministic assertions as automated tests. Cache the verified source artifact rather than requiring live network for unit tests. Visual checks are explicitly marked. A blocked external-data test must report blocked, not pass or silently skip the limitation.

| ID | Fixture/action | Required result |
|---|---|---|
| CU01 | Migrate registry twice | Two unique new task IDs first time; no duplicate second time |
| CU02 | Baseline30 or33 | New count32 or35; no unrelated profiles added |
| CU03 | Existing non-CU score fixtures | Numerical profile payloads unchanged; intentional top-level version/report changes excluded |
| CU04 | Compile CSU JSON | 15 unique species; five positive/ten negative signs |
| CU05 | Filter q<0.05 | Exactly the eight listed IDs |
| CU06 | B. stercoris rises | Its contribution rises; no fabricated direction-conflict flag |
| CU07 | B. wadsworthia rises | Its contribution falls |
| CU08 | Source coef/q reconciliation | Match within1e-9; coefficient not used as classifier weight |
| CU09 | Changed ZIP digest | Source-revision review blocks ingestion |
| CU10 | Matrix parse | 52 sample columns; statistics column excluded |
| CU11 | Terminal-species filter | 190 species; no parent+strain double counting |
| CU12 | Species sums<100 | No silent panel closure or claim of full production composition |
| CU13 | Permuted source columns | Explicit ID join still returns26 HC/26 CSU |
| CU14 | AMPLICON+GENOMIC archive row | Amplicon, not shotgun |
| CU15 | PRJNA1101664/OEP002997 | Excluded from human training/reference |
| CU16 | PRJNA1359023 | Isolate genome, not patient cohort |
| CU17 | SD versus CSD | OEP001229/OEP001891 retain separate phenotypes |
| CU18 | Possible cohort overlap | No assumed independent replication |
| CU19 | Ruminococcus genus only | Does not fill three different species slots |
| CU20 | Parent abundance plus SGB/strain children | Duplicate counting prevented |
| CU21 | E. coli/Shigella unresolved complex | Unavailable species slot unless equivalent bundle resolution established |
| CU22 | SD genus aggregate | Non-overlapping children; no overlap with R. bromii |
| CU23 | 1% versus0.01 fraction | Same transformed measurement under explicit units |
| CU24 | NaN/inf/negative/>1/boolean | Validation error, never finite patient score |
| CU25 | 19 distinct controls | Reference cannot fit |
| CU26 | 20 rows from one person | Rejected before reference fitting |
| CU27 | Constant control feature | Finite scale-floor transform, no divide-by-zero |
| CU28 | Query at log reference median | Contribution50 |
| CU29 | z=2, direction±1 | Contributions88.0797077978/11.9202922022 |
| CU30 | Sample alone versus arbitrary batch | Identical result; no inference fitting |
| CU31 | Taxon not in database | Not zero-imputed; lowers coverage |
| CU32 | Nondetection with known limit | Correct signed censoring endpoints |
| CU33 | Nondetection without limit | Omitted point, missingness bound retained |
| CU34 | One of two features,point80,bounds70–90 | Index80,coverage0.5,full-panel bounds35–95 |
| CU35 | No usable features | Null score, reason, bounds0–100; not healthy zero |
| CU36 | One of15 usable | Partial index explicitly limited-feature, broad bounds visible |
| CU37 | Full/strict different masks | Separate coverage; no independent-test claim |
| CU38 | Longitudinal mask/reference differs | Noncomparable flag or named common-mask recomputation |
| CU39 | No held-out calibration | Percentiles null; index has no percentage sign |
| CU40 | Equal-weight scores20,50,50,80,query50 | Midrank percentile50 |
| CU41 | Case/control calibration | Separate labelled distributions; no posterior probability |
| CU42 | Unknown age or8/17/70 with adult18–65 bundle | No participant-facing adult index; measurements/evidence retained |
| CU43 | Healthy adult,status absent,duration0 | Research computation allowed; no diagnosis assigned |
| CU44 | Diagnosis/duration conflict | Review flag, no silent rewrite |
| CU45 | Missing UAS7/UCT/response | No imputation or use as stool features |
| CU46 | Recent/unknown exposures | Correct flags; not presumed clean/in-distribution |
| CU47 | Profiler/database/denominator differs | Reference incompatibility blocks number |
| CU48 | New lane lacks validated QC | Analytic-validation-required state |
| CU49 | Shannon high/low | No CU index term; ecology displayed separately |
| CU50 | Conflicting genus/phylum study | Retained; does not overwrite species signs |
| CU51 | Functional assay not run | Not measured; no taxonomic imputation |
| CU52 | Overlapping LPS pathways | Not counted as independent taxonomy-score terms |
| CU53 | R. gnavus detected | No canonical-LPS-production claim |
| CU54 | SCFA/histamine DNA potential | No fabricated concentration,DAO status,MCAS or nutrient deficiency |
| CU55 | PWY0-1338 only | No clinical drug-resistance claim |
| CU56 | Unsupported inducible subtype | Evidence-not-established; no CSU substitution |
| CU57 | High index without diagnosis | Educational retrieval only; no prescription or eradication |
| CU58 | Low index in proposed donor | No pass/rank/clearance/transmission reassurance |
| CU59 | Healthy minors as controls | Not admitted to adult reference/donor clearance |
| CU60 | MR/review/isolate data | No patient WGS rows or invented feature weights |
| CU61 | LG-1 mouse intervention | Preclinical label, no human dosing/response claim |
| CU62 | Eczema FMT article | Not included as CU efficacy |
| CU63 | Missing strain ID | Explicit identity limitation, no substitution |
| CU64 | Positive and null probiotic studies | Balanced retrieval, not positive-only solutions |
| CU65 | Internal study protocols | No participant prescription/FMT dose or procedural leak |
| CU66 | >8 high-profile results | All cards/details retained; no truncation |
| CU67 | 15 CSU taxa incl.unavailable | All rows visible with source/status; no absence=deficiency |
| CU68 | Unavailable task | Reason/null display, not zero bar |
| CU69 | Non-CU report sections | Existing content retained after layout expansion |
| CU70 | PDF/narrow-screen visual check | Readable labels,bounds,rows,links; no color-only meaning |
| CU71 | Semantic language lint | Blocks affirmative diagnosis/probability/cure/clearance claims; permits disclaimers containing those words |
| CU72 | Run §7 self_test | Numerical assertions pass |
| CU73 | Discovery-subject evaluation | Label reconstruction, not external validation |
| CU74 | Offline inference | Frozen local artifacts suffice; no runtime research/download |
| CU75 | Handoff dependency audit | No required supplemental report/ledger/previous spec |
| CU76 | SD inferred histidine pathway | Not relabelled measured histamine production |
| CU77 | 2026 archive tissue/alias/date conflicts | Ingestion blocked until reconciled; H32 not labelled healthy by guess |
| CU78 | Bi2021 formulation | B. animalis LK011, not B. breve; pediatric evidence label |
| CU79 | Can2024 formulation/design | ATCC55730/nonrandomized; not DSM17938/RCT |
| CU80 | Godse2025/Dabaghzadeh2023 | Endpoint/source conflicts displayed; not unequivocal efficacy |
| CU81 | FMT case fields not verified | Null protocol/age/durability; no secondary-source filling |
| CU82 | Duplicate feature IDs or malformed intervals | Compile/validation failure, not silent acceptance |

End of normative specification. New evidence changes a reviewed source/profile revision and tests; it must not silently mutate frozen sample scores.
