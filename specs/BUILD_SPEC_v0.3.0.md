# BUILD SPEC 3 — Evidence-typed universal disease signature system

**Version:** 3.0.0 research specification
**Primary input:** human stool shotgun-metagenomic FASTQ
**Optional inputs:** measured stool/plasma metabolomics, metatranscriptomics, qPCR or spike-in microbial load, consented clinical metadata
**Output:** disease-associated microbiome **pattern concordance** — resemblance, not diagnosis or disease probability

Follow-on to `BUILD_SPEC.md` (gene capacity engine, v1.1.0) and `BUILD_SPEC_2.md`
(taxonomic engine and profile system). Self-contained; references nothing else.

---

## 1. Executive build decision

Build a versioned profile compiler that scores any documented disease-associated
microbiome signature without pretending every signature is equally mature or measured the
same way.

Three requirements are non-negotiable.

**1. Every registered condition with an observable feature returns a score.** Alzheimer's
disease, alopecia areata, and androgenetic alopecia are scored profiles, not
literature-only cards. When at least one supported feature passes analytical QC, the
profile returns a 0–100 pattern-concordance value. Evidence limitations appear as separate
coverage, maturity, transportability, conflict, out-of-distribution, and uncertainty
outputs alongside that value.

**2. Heterogeneous evidence stays heterogeneous.** A taxon, a sequence allele, a gene
family, a pathway, a measured metabolite, a predicted flux, an ecological metric, and a
frozen machine-learning model do not share one normalization function. Each evidence type
gets a typed observation contract and a typed scorer (Section 7).

**3. The reported percentage is resemblance, not risk.** `pattern_concordance = 80` means
the observed evidence resembles the published case pattern on that profile's frozen
reference scale. It does not mean an 80% chance of disease.

Standard report label:

> **Research-use disease-associated microbiome pattern concordance. Not a diagnosis,
> disease probability, screening result, treatment recommendation, or substitute for
> clinical evaluation.**

### 1.1 Runtime rule

For every registered profile:

- one or more supported features pass QC → return a score;
- only part of the signature is observable → return the score with status
  `scored_low_coverage` and a widened interval;
- the profile was derived from 16S and the sample is shotgun → score it and emit
  `assay_transport: amplicon_taxon_to_shotgun`;
- published function was PICRUSt-predicted and the sample has measured functions → score
  the explicitly transported module and emit `assay_transport: predicted_to_measured`;
- an optional measured modality is absent → omit that module and select the
  availability-mask calibrator; never impute it from DNA;
- zero supported features survive QC → `not_computable`, never a fabricated zero;
- the whole sample fails QC → `invalid_sample`.

### 1.2 First-wave profiles

| Family | Required task separation | Initial state |
|---|---|---|
| Alzheimer's continuum | `AD_SCS`, `AD_SCD`, `AD_MCI`, `AD_CLINICAL`, `AD_PRECLINICAL_AMYLOID` | scored research beta |
| Alopecia areata | pediatric shotgun, human taxonomic ensemble, optional measured VOC, mechanism | scored experimental |
| Androgenetic alopecia | female stool, male stool, sex-unspecified two-stratum range; scalp separate | scored experimental |
| ME/CFS | BioMapAI-compatible species and KO tasks plus production-native retrain | scored beta after repair |
| Long COVID / PASC | chronic state, symptom domain, and future risk kept separate | scored beta |
| IBS | overall plus IBS-C, IBS-D, IBS-M subtypes | scored experimental |
| High-maturity set | CRC, IBD composite, Crohn's, UC, T2D, cirrhosis, PDAC, Parkinson's, RRMS | production candidate after uniform validation |

---

## 2. Claim levels

Every observation and every report sentence carries one `claim_level`.

| Claim level | Meaning | Permitted example | Forbidden inference |
|---|---|---|---|
| `observed` | Directly measured in this specimen and assay | stool relative abundance of *F. prausnitzii* | absolute cell count without load measurement |
| `genomic_capacity` | DNA supports presence or completeness of genes or a pathway | butyrate-pathway gene coverage | pathway expression or metabolite concentration |
| `carrier_proxy` | Abundance of organisms with curated potential | abundance of butyrate-capable carriers | actual production rate |
| `flux_predicted` | Model-derived flux under declared diet constraints | MICOM ensemble flux | measured in-vivo flux |
| `external_measured` | Direct measurement from a separate assay | fecal imidazole propionate by LC-MS | value inferred from FASTQ |
| `model_output` | Output of a frozen model and calibrator bundle | BioMapAI-compatible species score | calibrated clinical probability |
| `association_only` | Literature-linked resemblance without adequate clinical validation | transported AGA pathway direction | diagnosis, risk, causality |

**Body site is part of identity.** `stool`, `scalp`, `oral`, `skin`, and `tumor`
observations never merge because a taxon name matches. Alopecia areata and androgenetic
alopecia are different diseases. Pediatric and adult AA are different populations. AD
stages are different prediction tasks. ME/CFS and PASC are different labels.

**Host reads are removed by default and not retained.** Host-genetic analysis from stool
is disabled unless a separate consent, governance, and validation path is explicitly
enabled, and never contributes to a microbial profile score.

---

## 3. Architecture

### 3.1 Core entities

| Entity | Role |
|---|---|
| `Assay` | What was measured: platform, chemistry, depth, body site, extraction, batch |
| `Observation` | One typed measurement with value, unit, QC state, and `claim_level` |
| `EvidenceAssertion` | One published claim: feature, direction, effect, source, cohort, assay of origin |
| `SignatureModule` | A set of assertions from one study or one coherent model, scored together |
| `Profile` | A disease task composed of modules, with fusion rules and calibrators |
| `CompatibilityContainer` | Pinned toolchain reproducing a published feature space exactly |
| `registry.lock` | Immutable manifest of every profile, model artifact, database, and digest |

### 3.2 Repository layout

```
gutscreen/
├── engines/          native production pipeline
├── compat/           pinned containers per published feature space
│   ├── metaad_2025/
│   ├── biomapai_2025/
│   └── gmwi2/
├── adapters/         one typed scorer per evidence type
├── profiles/         declarative profile YAML
├── namespaces/       controlled vocabularies + mapping files
├── calibrators/      per availability-mask calibration bundles
├── registry.lock
└── validation/
```

### 3.3 Migration from the existing application

- Wrap current DIAMOND/`rpoB` outputs as `GENE_ABUND` observations, preserving the
  existing copies-per-100-genomes unit, the exact protein reference-set fingerprint,
  mapping thresholds, and saturation QC. Keep `claim_level: genomic_capacity` — the nine
  pathways are not measured metabolites.
- Wrap current MetaPhlAn output as `TAX_REL` with the exact database identifier. A
  MetaPhlAn 3 `mpa_v31_CHOCOPhlAn_201901` space and a MetaPhlAn 4 SGB space are
  incompatible unless reprocessed or explicitly bridged.
- Convert each `BUILD_SPEC_2` feature into an `EvidenceAssertion` plus `SignatureModule`.
  Preserve source, direction, and original weight for audit, but recompute concordance —
  do not carry forward multiplied `w × q × s × c` values (Section 7.1).
- Preserve the existing PDF and `results.json` as immutable legacy artifacts. Generate v3
  results alongside during a shadow period with a field-level migration table.
- A legacy feature lacking an exact namespace and version becomes
  `unresolved_legacy_feature` and cannot enter a v3 score until curated.

---

## 4. Sample contract

### 4.1 Required manifest

Platform and chemistry; read length and layout; body site; extraction kit; library kit;
sequencing run, flowcell, and lane; collection date; preservative; freeze delay. Optional
but recorded when present: age, sex, BMI, country, antibiotics within 6 months, PPI,
metformin, other microbiome-active drugs, stool form, diet pattern, symptom onset date,
disease duration, treatment status.

### 4.2 QC gates

Input pairs, post-trim pairs, host fraction, non-host pairs, classified fraction, duplicate
rate, read-length profile, GC distribution, contamination flags, and `rpoB` fragment count.

Prespecified minimum: **500,000 usable non-host read pairs.** Below this, taxonomic
non-detection is `missing`, not absence — this distinction is enforced by the
presence/absence adapter (Section 7.2), not left to interpretation.

Run and flowcell identifiers are extracted from FASTQ headers and recorded. Batch effects
in metagenomics routinely exceed the biological differences being measured, so
cross-sample comparisons carry a batch-match flag.

### 4.3 Reference matching

Match on country, age band, and sex where the profile declares them mandatory. When a
matching stratum holds fewer than 30 samples, fall back to the full reference cohort and
set `population_transportability` accordingly. Never silently substitute an unmatched
reference.

---

## 5. Evidence types

| Code | Type | Engine | Native claim level |
|---|---|---|---|
| `TAX_REL` | Species/genus relative abundance | MetaPhlAn 4 | `observed` |
| `TAX_PRES` | Presence/absence | MetaPhlAn 4 + depth model | `observed` |
| `STRAIN` | Strain, SNV, allele | StrainPhlAn / MIDAS2 | `observed` |
| `PANCNV` | Pangenome copy number | MIDAS2 / PanPhlAn | `observed` |
| `GENE_ABUND` | Gene family abundance | DIAMOND (v1.1.0) | `genomic_capacity` |
| `SEQ_HMM` | Sequence/HMM with active-site logic | HMMER + decoys | `genomic_capacity` |
| `PATH_ABUND` | Pathway abundance | HUMAnN 3 | `genomic_capacity` |
| `PATH_COMPL` | Pathway completeness | frozen reaction rules | `genomic_capacity` |
| `CAZY` | Carbohydrate-active enzymes | dbCAN3 | `genomic_capacity` |
| `AMR` | Resistome | CARD / RGI | `genomic_capacity` |
| `CARRIER_ABUNDANCE` | Abundance of curated carriers of a pathway | MetaPhlAn + gene-to-taxon map | `carrier_proxy` |
| `ECO` | Diversity, richness, evenness | MetaPhlAn-derived | `observed` |
| `FLUX` | Predicted metabolite flux | AGORA2 / MICOM | `flux_predicted` |
| `METAB` | Measured metabolite | external LC-MS/GC-MS | `external_measured` |
| `MODEL` | Frozen published model output | compatibility container | `model_output` |

### 5.1 `CARRIER_ABUNDANCE`

Names a metabolite and resolves to the summed abundance of taxa whose reference genomes
carry the pathway's genes. Build the map from the existing DIAMOND reference sets; cache
it keyed by reference fingerprint.

**Carriage is not production.** The claim level is `carrier_proxy`, never `observed` or
`external_measured`. Its value is as an independent line of evidence alongside
`GENE_ABUND`: genes counted one way, organisms carrying them counted another. Agreement
between the two is a cross-evidence check; disagreement indicates a broken map or an
unusual community, and downgrades `technical_reliability`.

### 5.2 Assay transport

Evidence derived from a different assay than the sample binds explicitly, with the
transport labeled and a resolution penalty recorded.

| Transport | Permitted | Constraint |
|---|---|---|
| `amplicon_taxon_to_shotgun` | genus and species names with a validated mapping | never transfer ASV or OTU sequence IDs without a validated mapping; unresolvable groups become `missing`, not guessed |
| `predicted_to_measured` | PICRUSt-inferred function → measured shotgun function | separate module, separate sensitivity analysis, direction is a hypothesis test |
| `scalp_to_stool` | not permitted | different body site; requires a validated paired-site fusion model |

This replaces any notion that 16S or predicted-function evidence is unusable. It is
usable, labeled, and penalized.

---

## 6. Engines and compatibility containers

**Production-native core:** fastp → KneadData/Bowtie2 host removal (GRCh38 + CHM13 +
PhiX) → MetaPhlAn 4 (pinned database) → HUMAnN 3 → DIAMOND panels (v1.1.0, unchanged) →
optional MEGAHIT/Prodigal assembly for `CAZY`, `AMR`, `PATH_COMPL`.

**Compatibility containers** reproduce a published feature space exactly, because
published models are not portable across database releases. Each container pins every tool
version, every database build, and emits a digest recorded in `registry.lock`.

Required containers: `metaad_2025` (Section 10.1), `biomapai_2025` (Section 10.4),
`gmwi2` — the last must reproduce GMWI2's official MetaPhlAn 3-compatible feature space
rather than feeding MetaPhlAn 4 output into an incompatible formula.

**Provenance discrepancies are blocking.** When a source's methods, README, and
supplementary tables disagree on feature counts, record all claims, compute the imported
manifest size and SHA-256, and identify which list actually enters the model before
freezing an artifact.

**Audit source code for leakage.** Published performance is provisional until a
leakage-corrected reproduction. Check specifically for feature selection performed outside
cross-validation folds and for reused random seeds across repeated CV. Where leakage is
found, record it and treat the published metric as optimistic.

---

## 7. Typed scoring

### 7.1 Common output, not common likelihood

Every module maps its native evidence to a signed resemblance `r ∈ [-1, 1]`, where `+1` is
strongly case-concordant. How it obtains `r` is type-specific.

For a directional panel with valid features `V`:

```
S_m = 50 × (1 + Σ(w_i · r_i) / Σ(w_i))
```

Weights are frozen biological or statistical weights derived inside training data —
preferably shrinkage estimates based on effect precision — and clustered so correlated
features do not count as independent replications.

**Evidence quality, assay transport, and confidence never multiply `S_m`.** Multiplying
them changes the quantity being measured: a score built on weak evidence would drift toward
the midpoint and read as *biologically neutral* rather than *unknown*. Confidence is
reported as a separate vector (Section 7.5).

### 7.2 Typed adapters

| Module type | Required calculation |
|---|---|
| Relative abundance | Fixed-reference zero-aware log-ratio; transform by frozen case/control density or robust control scale. Never use the sample's own changing feature set as denominator |
| Presence/absence | Depth-dependent detection model; beta-binomial prevalence or frozen likelihood ratio. Nondetection at inadequate depth is `missing`, not absence |
| Gene/pathway abundance | Correct for feature length, breadth, identity, multimapping, and the source model's normalization; exact release |
| Pathway completeness | Score required and alternative reactions with frozen rules; report completeness and abundance separately |
| SNV/allele | Coverage-qualified allele likelihood with base and strand filters; `missing` below eligibility |
| Pangenome CNV | Species-coverage-qualified copy estimate in the exact pangenome build |
| Sequence/HMM | Calibrated bit-score, coverage, and active-site logic against positives and decoys; preserve homolog uncertainty |
| Measured metabolite | Assay-, matrix-, unit-, and batch-specific reference distribution. No DNA substitution |
| Predicted flux | Ensemble across declared diet scenarios; output median, scenario interval, sensitivity |
| Ecology | Metric-specific empirical case/control null. A cohort-level network is not computable from one sample |
| Frozen ML model | Exact training-fitted transformer and ordered schema, then estimator, OOD detector, frozen calibrator |

**Never infer direction from feature importance.** Tree-model importance has magnitude but
no case/control sign. Store `model_importance` and `direction_evidence` as independent
fields; a feature may have high importance and unknown direction.

### 7.3 Fallback and preferred transforms

For a prespecified direction feature with a compatible control reference:

```
z_i = (x_i − median(x_control)) / (1.4826 · MAD(x_control) + ε)
r_i = tanh(d_i · z_i / τ_i)
```

`d_i` is the published case direction; `τ_i` is frozen from derivation data. Label this
`direction_panel`, never a learned diagnostic model.

When individual-level data exist, prefer the empirical density ratio:

```
r_i(x) = 2σ[ log( (f_case,i(x) + ε) / (f_control,i(x) + ε) ) ] − 1
```

### 7.4 Coverage

Coverage is computed over independent evidence clusters, not raw feature count:

```
coverage_m = Σ(W_c over observed clusters) / Σ(W_c over expected clusters)
```

Return both `raw_feature_coverage` and `independent_cluster_coverage`. The latter governs
`scored_low_coverage` and interval expansion. A score built on 1 of 100 expected features
returns — but it returns visibly incomplete with a wide interval.

### 7.5 Confidence is a vector

```yaml
confidence:
  technical_reliability: 0.00-1.00
  independent_feature_coverage: 0.00-1.00
  evidence_maturity: P0-P5
  assay_transportability: native | labeled_transport
  population_transportability: 0.00-1.00 | unknown
  empirical_specificity: 0.00-1.00 | unknown
  out_of_distribution: false | true | unknown
  direction_conflict: none | low | material
  interval_width: 0-100
```

| Maturity | Meaning |
|---|---|
| `P0` | animal, case report, or mechanistic hypothesis only |
| `P1` | one small human cohort, amplicon-only, or inaccessible individual data |
| `P2` | usable human shotgun cohort, or multiple compatible human association cohorts |
| `P3` | independent human replication with coherent task and assay |
| `P4` | external population validation with frozen model and proper participant split |
| `P5` | prospective intended-use validation with calibrated clinical performance |

Every grade may support a research pattern score when a feature is measurable. No grade
authorizes a clinical claim.

### 7.6 Late fusion

- Large profiles: fit a stacker inside participant-grouped, study-aware nested validation.
  Inputs are module outputs plus an availability mask.
- Small literature profiles: average **independent study modules** with equal prespecified
  study weight, publishing every module score. A study contributing 20 correlated taxa must
  not overwhelm a study contributing one replicated feature.
- Fit or freeze a calibrator per permitted availability mask. Without one, return the
  transparent weighted module mean labeled `uncalibrated_pattern_concordance`.
- Never impute a measured metabolite, transcript, protein, or absolute load from DNA.

### 7.7 Profile result

`pattern_concordance_percent` plus:

- `case_typicality_percentile` — position within the frozen case distribution
- `control_tail_percentile` — extremeness relative to matched controls
- `uncertainty_interval` — bootstrap, measurement, and model interval
- `general_dysbiosis` — separate official GMWI2-compatible result
- `disease_specific_residual` — target signal after frozen general-dysbiosis adjustment
- `competing_profiles` and top-two margin
- empirical specificity statistics (Section 8)

### 7.8 Shape analysis across profiles

With the full library scoring, report above the individual results: count above the 60th
and 85th percentiles, mean empirical specificity of the high scorers, and correlation with
GMWI2. When most profiles score high with low mean specificity, the general-disturbance
conclusion is the headline. A report of this size is that many opportunities to over-read
noise, and the shape is what distinguishes a specific resemblance from a diffuse one.

---

## 8. Empirical specificity

Estimate specificity from challenge cohorts, not from how often a feature appears across
profile files — the latter is circular and changes whenever the registry is curated.

For each target profile evaluate: target versus matched healthy controls; target versus
every available non-target disease; target versus diseases sharing an
inflammation/oralization/diarrhea/constipation pattern; target false-positive rate within
every challenge cohort; worst-case non-target false-positive rate; target-versus-all-others
AUROC and AUPRC; top-two margin and error at each abstention threshold; and performance
after residualizing the general dysbiosis component.

Registry feature frequency may be displayed as `literature_uniqueness_prior`. It is a
labeled fallback prior, not empirical specificity.

---

## 9. Profile contract

```yaml
profile_id: AD_CLINICAL_METAAD_V1
version: 1.0.0
label: "Clinical Alzheimer's disease — stage-specific gut pattern concordance"
disease: alzheimers_disease
task: clinical_ad_vs_normal_control
body_site: stool
population: adult_chinese_cohort_derived
status: research_beta
evidence_maturity: P2

reference:
  cohort: metaad_2025_china
  control_class: NC
  match_on: [age_band, sex, country]

compatibility_container: compat/metaad_2025

modules:

  - module_id: METAAD_STRICT_DIRECTION_PANEL_V1
    type: direction_panel
    claim_level: observed
    evidence_type: TAX_REL
    assay_transport: native
    study_group: metaad_2025          # shares participants with the KO module
    features:
      - feature: Odoribacter_splanchnicus
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Barnesiella_viscericola
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Alistipes_communis
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Alistipes_finegoldii
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Alistipes_megaguti
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Alistipes_senegalensis
        direction: higher_in_case
        q_value: "<0.05"
        simulator_supported: true
      - feature: Alistipes_ihumii
        direction: higher_in_case
        q_value: "<0.05"
      - feature: Prevotella_intermedia
        direction: higher_in_case
        q_value: "<0.05"
      - feature: Streptococcus_ilei
        direction: lower_in_case
        q_value: "<0.05"
      - feature: Streptococcus_rubneri
        direction: lower_in_case
        q_value: "<0.05"

  - module_id: METAAD_STRICT_KO_PANEL_V1
    type: direction_panel
    claim_level: genomic_capacity
    evidence_type: PATH_ABUND
    study_group: metaad_2025          # correlated with the taxon module
    features:
      - {feature: K00857, direction: higher_in_case}
      - {feature: K00876, direction: higher_in_case}
      - {feature: K01031, direction: higher_in_case}
      - {feature: K03179, direction: higher_in_case}
      - {feature: K11782, direction: higher_in_case}
      - {feature: K18285, direction: higher_in_case}
      - {feature: K20036, direction: higher_in_case}
      - {feature: K01754, direction: lower_in_case}

  - module_id: METAAD_RELAXED_Q25_PANEL_V1
    type: sensitivity_panel
    status: exploratory
    note: "q < 0.25 exploratory set; reported as sensitivity, never called validated"

  - module_id: METAAD_LOCKED_MODEL_V2
    type: frozen_model
    claim_level: model_output
    status: pending_leakage_corrected_reproduction

fusion:
  method: independent_study_module_average
  study_group_covariance: true        # taxon and KO modules are one group, not two
  availability_masks: [tax_only, tax_plus_ko, full]

abstain_if:
  - usable_nonhost_pairs_lt: 500000
  - recent_antibiotics_days_lt: 90

output_notes:
  - >
    Derived from a Chinese staged cohort. Directions are source-cohort
    findings, not universal AD biology. Conflicting directions in other
    cohorts are held as separate assertions.
```

### 9.1 Assertion contract

Every `EvidenceAssertion` records: source DOI or accession; cohort size and composition;
originating assay and body site; feature ID in a declared namespace; mapping confidence if
transported; direction; effect size and dispersion where published; significance and
correction method; and the study group it belongs to for covariance adjustment.

**Same-participant modules are one evidence group.** A taxon panel and a KO panel from the
same stool DNA are not two independent replications.

### 9.2 `registry.lock` and compiler failures

The loader fails and refuses to score when: a feature has no resolvable namespace ID; a
transported feature has no validated mapping; a module declares a direction its source
does not support; a frozen model artifact digest does not match; a compatibility container
version drifts; or two modules from the same participants are declared independent.

---

## 10. Profile catalog

### 10.1 Alzheimer's continuum — five stage-specific profiles

**Primary source.** Jia L, Ke Y, Zhao S, et al. Metagenomic analysis characterizes
stage-specific gut microbiota in Alzheimer's disease. *Mol Psychiatry* 2025;30:3951–3962.
https://doi.org/10.1038/s41380-025-02973-7 — fecal shotgun, Chinese cohort of **476
participants across five stages**, reporting average AUC 0.80 in cross-validation and 0.75
in independent validation, with over 10% of microbial species and gene families
significantly altered across progression.

| Stage | N | Profile ID |
|---|---:|---|
| Normal control | 63 | reference class |
| Subjective cognitive symptoms | 82 | `AD_SCS_METAAD_V1` |
| Subjective cognitive decline | 90 | `AD_SCD_METAAD_V1` |
| Mild cognitive impairment | 119 | `AD_MCI_METAAD_V1` |
| Clinical AD | 122 | `AD_CLINICAL_METAAD_V1` |

Data: `PRJCA022804` / `CRA014435` (NGDC); analysis repository
`github.com/JiaLonghao1997/metaAD_analysis`. Verify the commit digest and every table row
count at import — record any discrepancy between methods text, README, and supplementary
tables as blocking before freezing an artifact.

**At strict `q < 0.05`, no species qualified for SCS, SCD, or MCI.** Only clinical AD has a
species-level strict panel. Stage profiles below clinical AD score from KO panels and the
relaxed sensitivity set, and their reports must say so.

**Three operators per stage:** `STRICT_DIRECTION_PANEL` (primary, transparent),
`RELAXED_Q25_PANEL` (sensitivity only, never called validated), and `LOCKED_MODEL`
(preferred once leakage-corrected reproduction completes, and it does not erase the panel
results).

**Complementary cohorts, each its own module:**

- **AlzBiom** — https://doi.org/10.3389/fnins.2022.792996, `PRJEB47976`; 75 amyloid-positive
  AD and 100 controls with same-site holdout.
- **Ferreiro et al.** preclinical amyloid — PMC10680783, `PRJNA798058`; 164 cognitively
  unimpaired (49 Aβ+, 115 Aβ−). Creates `AD_PRECLINICAL_AMYLOID_V1`, which is **not** a
  clinical-AD label.
- **Kang et al. 2025** MARS/AGMP — PMC12221809, AD Knowledge Portal `syn53071655`; MARS
  n=232 with independent AGMP n=448 replication of a taxon log-ratio. The strongest
  cross-cohort replicated taxonomic ratio available; weight accordingly.
- **Fan et al. 2025** Taiwan MCI — PMC12123878, `PRJNA1258384`; 119 MCI, 320 normal.
  Direction and contradiction testing.
- Longitudinal: normal-to-MCI PMC12112180 / `PRJEB85717`; MCI-to-AD PMC10855790.
  Hypothesis-only until external validation.

**Separate metabolite-capacity lane.** The nine existing v1.1.0 panels carry AD-relevant
mechanisms — secondary bile acids (elevated DCA with reduced primary bile acid, ratio
tracking cognitive decline; DCA interacts with nicastrin in the gamma-secretase complex),
TMA/TMAO (elevated in cerebrospinal fluid in AD and MCI), butyrate (drives GSK3β acetylation
at lysine 15, regulating tau phosphorylation), propionate (Mendelian randomization
associates higher production with increased AD likelihood, with documented protective
mitochondrial effects creating direction conflict), indole-3-propionate (interferes with
amyloid-beta fibril formation), and imidazole propionate (impairs insulin signalling via
p38-gamma and mTORC1).

These bind as `GENE_ABUND` at `claim_level: genomic_capacity` in a **separate module** from
the taxonomic panels. They are mechanism lanes, not stage classifiers, and do not fuse into
the stage score. The imidazole-propionate/`urdA` lane stays independent — Vemuganti et al.
found `urdA` carriage related only weakly to plasma ImP across 294 samples.

### 10.2 Alopecia areata

Task separation: pediatric AA, adult AA, patchy AA, alopecia totalis, alopecia universalis
where labels permit. An adult sample still scores against the pediatric module but is
flagged OOD with a wide population-transport interval.

**Native shotgun module — `AA_PEDIATRIC_STOOL_V1`.** Rangu S, Lee J-J, Hu W, Bittinger K,
Castelo-Soccio L. *JID Innov* 2021;1:100051. https://doi.org/10.1016/j.xjidi.2021.100051 —
41 children aged 4–17 with AA and 41 unaffected siblings, shotgun metagenomics, median ~2.7
million QC read pairs per sample. *Ruminococcus bicirculans* (NCBI Taxonomy 1160721) was
decreased in AA. Alpha diversity and Bray-Curtis were null; a small Jaccard difference was
significant. Twenty gene orthologs differed between groups.

**Encode the KO direction conflict explicitly.** The publication's prose reverses `gerKA`
(K06295), `gerKC` (K06297), and `fbpA` (K02012) relative to its figure, and creates
additional ambiguity around K07552. Compile `RANGU_FIGURE_ENCODING_V1` and
`RANGU_TEXT_ENCODING_V1`, score both as sensitivity modules, expose the range and delta, and
mark all 20 KO directions `direction_uncertain` for production fusion pending author
clarification or raw reanalysis. The taxon result remains independently scoreable. Never
silently pick a direction.

The taxon and KO signals come from the same participants and stool DNA — one evidence group
with covariance adjustment, not two replications.

**Human taxonomic modules**, each independent, mapped with `amplicon_taxon_to_shotgun`:

| Module | Case-higher | Case-lower | Note |
|---|---|---|---|
| Moreno-Arrones 2020 (15 AU / 15 control) | *Holdemania filiformis*, Erysipelotrichaceae, Lachnospiraceae, *Parabacteroides johnsonii*, Clostridiales vadinBB60, *Bacteroides eggerthii*, *Parabacteroides distasonis*, Eggerthellaceae | *Phascolarctobacterium succinatutens*, *Dorea longicatena* | internal two-marker AUROC 0.804 (CI 0.633–0.976), no external validation |
| Lu 2021 (33 AA / 35 control) | *Blautia*, *Anaerostipes*, Erysipelotrichaceae, *Dorea*, *Collinsella*, *Megasphaera*, *Achromobacter* | Bacteroidetes, Fusobacteria | internal RF AUROC 0.935 train / 0.867 holdout; low-abundance environmental genera need contamination review |
| Bain 2022 (41 AA / 19 control) | high-SALT subgroup: *Alistipes*, *Bacteroides*, *Barnesiella*, Lachnospiraceae NK4A136, *Tyzzerella*, *Ruminiclostridium*, *Erysipelatoclostridium* | *Lachnospira*, *Lachnoclostridium*, *Ruminococcus* | overall tests largely null after correction; exploratory subgroup only |
| Lee 2024 (19 AA / 20 control) | *Blautia*, *Eubacterium* g5 | *Bacteroides* | alpha diversity null |
| Nikoloudaki 2024 (24 AA / 18 control) | Firmicutes, Clostridia, Lachnospirales, Lachnospiraceae, *Blautia*, *Holdemania* | *Coprococcus* | `PRJNA857562`; deposited run count requires reconciliation; internal DIABLO error ~35%; collapse correlated ranks into one lineage cluster |
| Juhasz 2020 (25 / 25) | Bacilli, Lactobacillales | — | fecal fungal and diversity results largely null |

**Conflicts do not resolve into a pooled consensus.** *Bacteroides* is case-enriched in the
severe Moreno and Bain branches but lower in Lee. *Dorea* is higher at genus level in Lu
while *D. longicatena* is control-associated in Moreno. Return the study-module vector,
heterogeneity, leave-one-study-out score range, and whether any direction survives in at
least two independent cohorts.

**Mechanism module**, `claim_level: genomic_capacity`, maturity P0, scored separately:
biotin biosynthesis (`bioA`, `bioB`, `bioC`, `bioD`, `bioF`, `bioH`) — Hayashi A, et al.
*Cell Rep* 2017;20:1513, germ-free mice under biotin deprivation developed alopecia with
*Lactobacillus murinus* overgrowth, reversed by supplementation; and propionate capacity
(`pct`, `pduP`, `lcdA`, existing panel) — Borde and Åstrand observed regrowth in five of five
C3H/HeJ mice then did not reproduce it, so this feature carries `direction_conflict:
material`.

**Optional measured VOC module** when fecal volatile metabolomics is supplied
(`claim_level: external_measured`). Never imputed from DNA.

### 10.3 Androgenetic alopecia

Biologically distinct from AA. Never merge their evidence or labels.

Source: Jung et al. 2022, *Front Microbiol*
https://doi.org/10.3389/fmicb.2022.1076242 — 95 Korean AGA participants (49 women, 46 men)
and 46 controls (25 women, 21 men), stool and scalp 16S V4–V5. Gut: `PRJNA891926`; scalp:
`PRJNA891901`.

No strong overall alpha-diversity difference or clear PCoA separation. Signals were
sex-specific, and male age distributions differed materially (control mean ~32.05 years,
AGA ~43.61). **Age and sex are mandatory matching variables.** Maturity P1.

| Stratum | Direction in AGA | Taxon |
|---|---|---|
| Female | lower | *Bifidobacterium* |
| Female | higher | *Lachnoclostridium*, Oscillospiraceae NK4A214, *Escherichia/Shigella* |
| Male | lower | *Butyricicoccus*, Ruminococcaceae incertae sedis, *Parabacteroides*, *Alistipes*, unclassified Oscillospiraceae |
| Male | lower | Christensenellaceae R7 — `nominal_status: not_below_0.05` (source figure ~p=0.058); include and exclude in a sensitivity result |

Functional directions were PICRUSt2-inferred, so they bind as `predicted_to_measured` in a
separate hypothesis-test module: female lower amino-sugar/nucleotide-sugar metabolism,
primary and secondary bile-acid biosynthesis; female higher siderophore NRP biosynthesis
and bacterial invasion of epithelial cells.

*Escherichia/Shigella* and several family-level groups cannot be resolved cleanly by
marker-based shotgun taxonomy. The mapping file records resolution loss and alternatives;
unresolved groups become `missing`.

**Scalp stays a separate profile.** Do not fuse scalp and stool without a validated
paired-body-site model.

### 10.4 ME/CFS

Primary: Xiong R, Aiken E, Caldwell R, et al. *Nat Med* 2025;31:2991–3001.
https://doi.org/10.1038/s41591-025-03788-3 — 249 participants (153 ME/CFS, 96 healthy),
four years, 515 timepoints, 479 with stool. Data `PRJNA1125469`; code, trained model, and
processed dataset at `github.com/ohlab/BioMapAI`.

**Bind the Species and KEGG sub-models only.** BioMapAI trains separate models per data
type and reports each independently; only Species and KEGG are reproducible from stool.
The combined figure integrates plasma metabolomics, immune profiling, and blood labs and
must not be cited beside a stool-only score.

**Two source-stated constraints.** The reported biomarkers were calculated on the entire
dataset and not validated on held-out data — maturity P2, `direction_evidence` separate
from `model_importance`. And in ordination, except for clinical scores, controls are
indistinguishable from patients. A direction panel is a distance-based construct, so an
uninformative panel result here is the expected outcome; the frozen `MODEL` module is the
operator that can recover signal, and its performance must be re-derived under nested
participant-grouped validation before use.

Named feature: *Dysosmobacter welbionis*, higher in case — top-ranked species feature.

**Duration stratification.** Xiong R, et al. *Cell Host Microbe* 2023;31:273 (`PRJNA878603`,
`github.com/ohlab/MECFS_2021`) found short-term patients (<4 yr, n=75) showed significant
microbial dysbiosis while long-term patients (>10 yr, n=79) had largely resolved microbial
dysbiosis while retaining metabolic abnormalities. Separate strata, never averaged; low
concordance in long-duration illness matches the published pattern.

**Cross-evidence check.** Guo C, et al. *Cell Host Microbe* 2023;31:288 (106 cases, 91
controls, multi-site US) found reduced *F. prausnitzii* and *E. rectale* driving butyrate
deficiency, confirmed by qPCR and measured SCFAs. Check the `TAX_REL` claim against the
existing `GENE_ABUND` butyrate panel — 82.8 copies per 100 genomes on sample A02 from
1,917 fragments — and against `CARRIER_ABUNDANCE`. Disagreement downgrades
`technical_reliability`.

### 10.5 Long COVID / PASC

Chronic-state, symptom-domain, and future-risk tasks are separate profiles. Features from
the acute-risk task never read as chronic biomarkers.

Calibration: in the largest independent prospective outpatient study (799 participants),
species AUROC was 0.62 and genus 0.59, against 0.72 for clinical variables alone. A
community study of 2,561 participants found no association with illness duration. The
frequently quoted 0.96 classifies a study's own internally defined enterotypes.

Datasets: `PRJNA714459` (Liu, *Gut* 2022), `PRJNA1004982` (Su, *Cell Host Microbe* 2024),
`PRJNA1167328` (Comba, *Gut Microbes* 2026), `PRJNA1181655` (Blankestijn).

**Verify label availability before building validation.** Several of these publish reads
without a sample-to-outcome key. Where labels are absent, restrict to unsupervised work and
record that discrimination is unmeasured.

**Cohort independence.** `PRJNA714459`, `PRJNA876804`, `PRJNA841786`, `PRJNA1004982`,
`PRJNA995807`, and `PRJNA1056888` come from potentially overlapping programs at one
institution — one study group with covariance adjustment, and held out whole during
validation splits.

### 10.6 IBS

Overall profile plus IBS-C, IBS-D, and IBS-M subtype tasks. Also serves as **the specificity
challenge** for Long COVID and ME/CFS: roughly a quarter of COVID microbiome biomarkers also
replicate in IBS, so a sample scoring high on both is the informative comparison.

### 10.7 High-maturity set

| Profile | Primary source | Notes |
|---|---|---|
| Colorectal cancer | Wirbel J, et al. *Nat Med* 2019;25:679 — 8-study meta-analysis, n=768, 29 core species at FDR<1e-5, signatures retained accuracy across studies | Binds `TAX_REL`, `GENE_ABUND` (*clbB*, *bft*, *fadA*), `CAZY`. **The validation gate** — the only profile with labeled public data. Not a cancer screen; validated clinical tests exist |
| Colorectal adenoma | curatedMetagenomicData labels | Lower-effect natural test of graded signal |
| Crohn's disease | Lewis JD, et al. *Cell Host Microbe* 2015 (`SRP057027`); HMP2; Franzosa 2019 | **Dysbiotic-subgroup stratification required** — only 28 of 85 newly diagnosed patients were distinctly dysbiotic |
| Ulcerative colitis | Machiels K, et al. *Gut* 2014 | Decreased *R. hominis* and *F. prausnitzii* |
| Liver cirrhosis | Qin N, et al. *Nature* 2014 | Oral taxa translocating to gut; *Veillonella* and *Streptococcus* higher — profile-local direction, see 10.9 |
| Type 2 diabetes | Qin 2012; Karlsson 2013; **Forslund K, et al. *Nature* 2015;528:262** | **Requires metformin status; abstains when unknown.** Much of the published signature is the drug |
| Parkinson's disease | Wallen ZD, et al. *Nat Commun* 2022;13:6958 (n=490/234) | *Akkermansia* higher here, beneficial elsewhere — `direction_conflict: material` |
| Multiple sclerosis | iMSMS Consortium *Cell* 2022 — 576 household-control pairs | Best confounder control in the library |
| ACVD | Jie Z, et al. *Nat Commun* 2017;8:845 | Cross-evidence check against existing `cutC`/`cntA` |
| Rheumatoid arthritis | Scher 2013; Zhang 2015 | Binds `STRAIN` — *P. copri* at strain level, since genus *Prevotella* has no stable direction |
| MASLD | Loomba 2017; Hoyles 2018 | Heavy T2D/obesity confounding; expect low empirical specificity |
| MDD | Valles-Colomer 2019; Radjabzadeh 2022 | Opt-in flag required; separate section; excluded from headline shape analysis |
| AMR burden | CARD/RGI | Binds `AMR`; relevant to donor comparison work |

Additional registered profiles: T1D, celiac (binds `CAZY`), ankylosing spondylitis,
hypertension, chronic kidney disease (functional-first via existing p-cresol and indole
panels), obesity (binds `ECO` only — Sze and Schloss showed taxonomic effect sizes too small
for classification and F:B does not replicate).

### 10.8 Cross-disease resources

GMrepo v3 (`gmrepo.humangut.info`) — 118,965 curated samples, 31,917 shotgun, 302
phenotypes, processed with MetaPhlAn 4.1.0. **Database presence is cohort discovery, not
validation.** Every phenotype used requires label audit, uniform reprocessing, grouped
validation, and cross-disease challenge before a profile leaves experimental status.

curatedMetagenomicData for the healthy reference cohort (minimum 300 samples, target 1,000+,
full depth, versioned manifest). MicrobiomeHD for the shared-response prior only, kept in
its own 16S space. Duvallet C, et al. *Nat Commun* 2017;8:1784 — the finding that many
disease-associated changes are a shared response to illness rather than disease-specific
markers, which motivates Section 8.

### 10.9 Direction conflicts

At genus level these carry per-profile directions only: *Faecalibacterium*, *Streptococcus*,
*Veillonella*, *Prevotella*, *Blautia*, *Lachnospiraceae*, *Alistipes*, *Bacteroides*,
*Dorea*. The loader requires species or strain resolution for any feature on this list, or
an explicit `direction_conflict` declaration.

---

## 11. New gene panels

**Biotin biosynthesis** — `bioA`, `bioB`, `bioC`, `bioD`, `bioF`, `bioH`. Tenth panel for
the v1.1.0 engine, required by 10.2. Same reference-plus-decoy discipline as existing
panels. Genes are widespread, so expect a compressed reference distribution — inspect the
shape before assigning percentile bands and report without bands if the interquartile range
is narrow.

**CRC virulence** — `clbB` (pks island, colibactin), `bft` (*B. fragilis* toxin), `fadA`.
Required by 10.7.

---

## 12. Confounders

| Confounder | Effect | Handling |
|---|---|---|
| Metformin | Produces much of the published T2D signature | Required for T2D; abstain if unknown |
| Antibiotics | Dominates for months | Abstain within 90 days |
| PPIs | Large documented shift | Record; reduce `technical_reliability` |
| Gluten-free diet | Large in prevalent celiac | Required match for celiac |
| Geography | Often exceeds disease effect | Match on country; set `population_transportability` |
| Age | Large across lifespan; critical for T1D and AGA | Mandatory match where declared |
| Sex | AGA signals are sex-specific | Mandatory match where declared |
| Diet and fiber | Days-to-weeks variation | Record; report |
| Stool consistency | Correlates with diversity | Record if available |
| Hospitalization / ICU | Drives *E. coli*, *Klebsiella*, AMR | Record; penalize those features |
| Kidney function | Determines whether uremic-toxin capacity matters | Required context for CKD |
| Sequencing batch | Can exceed biological signal | Record run and flowcell |

The confounder ledger is a required report section: what was recorded, what was missing,
which profiles abstained.

---

## 13. Report contract

**Shape analysis first** (7.8), then the ranked profile table with concordance, evidence
badges, maturity grade, and confidence vector summary.

Every score renders with: bound evidence types, `claim_level` per module, coverage,
uncertainty interval, assay transport labels, and direction-conflict status. Module scores
are always individually visible — a fused number that cannot be decomposed is a black box.

**Abstention renders with reasons** and what would need to change. Never a blank or a zero.

**Unbound sources render with `would_bind_if`**, so the report documents its own gaps.

**Opt-in sections** (MDD) require an explicit flag and render separately.

The existing metabolite section continues, gaining the biotin panel as the tenth pathway.
`results.json` carries every module, feature, adapter output, coverage figure, confidence
dimension, calibrator ID, and every version identifier: profile version, container digests,
MetaPhlAn database, DIAMOND reference fingerprint, reference cohort manifest, pipeline
version.

---

## 14. Validation

| # | Step | Passing criterion |
|---|---|---|
| V1 | CRC end-to-end on labeled curatedMetagenomicData | AUROC with CI, compared against published cross-cohort performance |
| V2 | Per-feature direction reproduction against GMrepo | Reported per feature, not as an aggregate; failures lose weight or drop |
| V3 | Cross-disease specificity matrix (Section 8) | Diagonal dominates; matrix published including off-diagonal signal |
| V4 | Leakage-corrected reproduction of every imported model | Nested, participant-grouped, feature selection inside folds |
| V5 | Cross-evidence agreement, `GENE_ABUND` vs `CARRIER_ABUNDANCE` | Correlation reported across the reference cohort |
| V6 | Synthetic communities with known composition | Sensitivity, specificity, quantitative slope and R², as the existing harness reports |
| V7 | Confounder sensitivity | Score movement per confounder published; if metformin moves T2D more than T2D does, that goes in the README |
| V8 | Reference stability bootstrap | Percentile stabilization curve published |
| V9 | Availability-mask calibration | Each permitted mask has a fitted or frozen calibrator, or output is labeled uncalibrated |

---

## 15. Adding a profile

1. Locate evidence; record the originating assay and body site.
2. Assign maturity P0–P5 per source, not per profile.
3. Choose evidence types per Section 5; declare transports explicitly.
4. Write one `SignatureModule` per study. Same-participant modules share a `study_group`.
5. Where a source's figure and prose disagree, compile both encodings and score them as
   sensitivity modules with the delta exposed.
6. Declare fusion, availability masks, and abstention rules.
7. Register in `registry.lock` with digests.
8. Validate: V2, then add to V3 and confirm the diagonal still dominates.

---

## 16. Build order

| # | Step | Done when |
|---|---|---|
| 1 | Typed adapter interface and evidence-type registry (Sections 5, 7.2) | `TAX_REL`, `TAX_PRES`, `GENE_ABUND`, `ECO` adapters pass unit tests on synthetic input |
| 2 | Profile contract loader with `BUILD_SPEC_2` shim and compiler failures (9.2) | Every existing profile loads; malformed profiles fail loudly |
| 3 | Confidence vector and coverage (7.4, 7.5) | Score and confidence computed independently; verified non-multiplicative |
| 4 | Reference cohort build, MetaPhlAn 4.1.0 pinned (4.3) | ≥300 samples, stratified query, versioned manifest |
| 5 | `gmwi2` compatibility container | Reproduces official MetaPhlAn 3-compatible output |
| 6 | High-maturity profiles: CRC, adenoma, Crohn's (subgroups), UC, cirrhosis, T2D | All score; T2D abstains without metformin status |
| **7** | **GATE — V1 CRC validation and V3 specificity matrix** | **AUROC with CI published; diagonal dominates** |
| 8 | `metaad_2025` container; five AD stage profiles (10.1) | Strict panels score; relaxed set labeled sensitivity; locked model pending V4 |
| 9 | Biotin and CRC virulence panels (Section 11) | Pass V6 like the existing nine |
| 10 | AA profile with dual KO encodings (10.2) | Both encodings score; delta exposed; taxon module independent |
| 11 | AGA profile, sex-stratified (10.3) | Sex and age matching enforced; PICRUSt module labeled `predicted_to_measured` |
| 12 | `biomapai_2025` container; ME/CFS profile (10.4) | Sub-models bound; V4 completed before the MODEL operator ships |
| 13 | `STRAIN`, `AMR`, `CAZY`, `PATH_COMPL` adapters | Each binds in ≥1 shipped profile |
| 14 | Remaining profiles: Parkinson's, MS, ACVD, RA, MASLD, MDD, T1D, celiac, AS, hypertension, CKD, obesity, Long COVID, IBS | Each passes V2 |
| 15 | `FLUX` adapter, AGORA2/MICOM | Separate section; direction-only comparison |
| 16 | Report contract (Section 13) | Shape analysis renders above individual scores |
| 17 | Re-run V3 with the full library | Diagonal dominates; matrix republished |

**Step 7 is a hard gate.** CRC is the only profile with labeled public data, so it is the
only place the engine is checkable against ground truth.

**Step 17 repeats step 7** because specificity estimates depend on the challenge set;
adding profiles changes them.

---

## 17. Sources

**Framework** — Duvallet C, et al. *Nat Commun* 2017;8:1784 · Gloor GB, et al. *Front
Microbiol* 2017;8:2224 · GMWI2: Chang D, et al. *Nat Commun* 2024;15:7447, code
`github.com/danielchang2002/GMWI2` · SIAMCAT: Wirbel J, et al. *Genome Biol* 2021;22:93 ·
curatedMetagenomicData: Pasolli E, et al. *Nat Methods* 2017;14:1023 · GMrepo v3: *Nucleic
Acids Res* 2026;54(D1):D734 · MetaPhlAn 4: Blanco-Míguez A, et al. *Nat Biotechnol*
2023;41:1633 · HUMAnN 3: Beghini F, et al. *eLife* 2021;10:e65088 · StrainPhlAn: Truong DT,
et al. *Genome Res* 2017;27:626 · MIDAS2: Zhao C, et al. *Bioinformatics* 2023 · CARD 2023:
Alcock BP, et al. *Nucleic Acids Res* 51:D690 · dbCAN3 · AGORA2: Heinken A, et al. *Nat
Biotechnol* 2023;41:1320 · MICOM: Diener C, et al. *mSystems* 2020

**Alzheimer's** — Jia L, et al. *Mol Psychiatry* 2025;30:3951.
https://doi.org/10.1038/s41380-025-02973-7 · AlzBiom: *Front Neurosci* 2022,
https://doi.org/10.3389/fnins.2022.792996 · Ferreiro et al. PMC10680783 · Kang JW, et al.
*Alzheimers Dement* 2025;21:e70417 · Fan et al. PMC12123878 · MahmoudianDehkordi S, et al.
*Alzheimers Dement* 2019;15:76 · Vogt NM, et al. *Alzheimers Res Ther* 2018;10:124 ·
Zhang Y, et al. *Mol Psychiatry* 2023;28:4421 · Sanna S, et al. *Nat Genet* 2019;51:600 ·
Pappolla MA, et al. *Neurobiol Dis* 2021;156:105403 · Vemuganti V, et al. *Nat Commun*
2026;17:8017

**Alopecia areata** — Rangu S, et al. *JID Innov* 2021;1:100051.
https://doi.org/10.1016/j.xjidi.2021.100051 · Moreno-Arrones et al.
https://doi.org/10.1111/jdv.15885 · Lu et al. https://doi.org/10.1016/j.jdermsci.2021.04.003
· Bain et al. https://doi.org/10.1093/cei/uxac088 · Lee JH, et al. *IJMS* 2024;25:4256 ·
Nikoloudaki O, et al. *Nutrients* 2024;16:858, `PRJNA857562` · Juhasz et al.
https://doi.org/10.25251/skin.4.1.4 · Hayashi A, et al. *Cell Rep* 2017;20:1513 · Borde A,
Åstrand A. *Expert Opin Ther Targets* 2018;22:503

**Androgenetic alopecia** — Jung et al. *Front Microbiol* 2022,
https://doi.org/10.3389/fmicb.2022.1076242, `PRJNA891926` (gut), `PRJNA891901` (scalp)

**ME/CFS** — Xiong R, et al. *Nat Med* 2025;31:2991, `PRJNA1125469`,
`github.com/ohlab/BioMapAI` · Xiong R, et al. *Cell Host Microbe* 2023;31:273,
`PRJNA878603`, `github.com/ohlab/MECFS_2021` · Guo C, et al. *Cell Host Microbe* 2023;31:288

**High-maturity set** — Wirbel J, et al. *Nat Med* 2019;25:679 · Qin N, et al. *Nature*
2014;513:59 · Lewis JD, et al. *Cell Host Microbe* 2015;18:489 · Machiels K, et al. *Gut*
2014;63:1275 · Qin J, et al. *Nature* 2012;490:55 · Karlsson FH, et al. *Nature*
2013;498:99 · Forslund K, et al. *Nature* 2015;528:262 · Wallen ZD, et al. *Nat Commun*
2022;13:6958 · iMSMS Consortium *Cell* 2022;185:3467 · Jie Z, et al. *Nat Commun*
2017;8:845 · Scher JU, et al. *eLife* 2013;2:e01202 · Zhang X, et al. *Nat Med*
2015;21:895 · Valles-Colomer M, et al. *Nat Microbiol* 2019;4:623 · Loomba R, et al. *Cell
Metab* 2017;25:1054 · Sze MA, Schloss PD. *mBio* 2016;7:e01018-16

**Long COVID** — Liu Q, et al. *Gut* 2022;71:544 · Su Q, et al. *Cell Host Microbe*
2024;32:651 · Comba IY, et al. *Gut Microbes* 2026 · Österdahl MF, et al. *Sci Rep* 2023
