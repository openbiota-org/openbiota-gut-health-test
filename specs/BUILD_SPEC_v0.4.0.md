# BUILD SPEC 4 — Evidence-gated microbiome findings, interventions, age, and visual reporting

**Specification version:** 4.0.0  
**Date:** 2026-09-07  
**Input:** human stool shotgun-metagenomic FASTQ, plus optional validated external assays and consented metadata  
**Primary outputs:** reproducible microbial observations, matched-reference findings, research-use disease-pattern concordance, independently curated research-evidence summaries, microbiome-predicted chronological age, and an experimental FMT research prescreen  
**Clinical status:** research use only; not a diagnosis, disease probability, prescription, donor clearance, regulatory authorization, or substitute for clinical testing

This specification is the next build after the live system described by `BUILD_SPEC_v03.md`. It preserves the v03 disease-profile engine and adds a separate findings-to-evidence-to-intervention layer. It also fixes reference matching, expected-but-not-detected organisms, microbiome age, report ordering, and visual/PDF accessibility.

**Normative dependency:** v04 supersedes v03 only where this document explicitly changes behavior. The v03 typed-observation system, profile definitions, source assertions, score semantics, and 30-task live disease-signature registry remain normative and must be frozen beside this file. V04's target registry contains exactly 33 profile tasks: the 30 live tasks plus `AD_SCS_METAAD_V1`, `AD_SCD_METAAD_V1`, and `PDAC_STOOL_MULTINATIONAL_V1`. No implementer may reinterpret or silently rebuild disease signatures from the intervention literature below.

---

## 0. Executive build decision

Build three systems that exchange typed records but never collapse into one score:

1. **Measurement and reference interpretation:** what was observed, whether detection was technically possible, and where it sits in a compatible healthy/control reference.
2. **Disease-pattern research scoring:** how closely the sample resembles a frozen published signature, with coverage, specificity, uncertainty, transportability, and out-of-distribution status.
3. **Intervention evidence retrieval:** what exact intervention protocols have evidence for the exact finding, symptom phenotype, external lab result, or clinician-confirmed condition, after safety and applicability gates.

There is no direct `microbe -> treatment` shortcut. A disease-pattern percentile cannot by itself unlock a drug, antimicrobial, probiotic, or FMT instruction. It can retrieve an unranked research-evidence summary and an independently justified clinical-follow-up route. “Option to discuss” language requires claim-specific clinical, regulatory, and substantiation approval; it is not granted automatically by an evidence lane.

The universal report statement is:

> This is a research-use analysis of stool shotgun-metagenomic data. It reports measurements, matched-reference comparisons, and resemblance to published microbiome patterns. It does not diagnose or rule out disease, establish infection or eradication, measure metabolites unless an external assay was supplied, predict that a treatment will work, prescribe treatment, or clear an FMT donor.

### 0.1 What v04 must accomplish

- Retain all 30 live disease-profile tasks, including alopecia areata, androgenetic alopecia, the Alzheimer continuum, ME/CFS, and Long COVID; restore the PDAC task promised by v03 but absent from the live registry.
- Render every scored, declined, low-coverage, abstained, or non-computable profile; remove the hidden maximum-eight detail-page limit.
- Show what is detected, relatively low, relatively high, commonly detected in matched references but not detected here, uncommon here, potentially concerning only after qualification, and analytically indeterminate.
- Add a research-use microbiome-predicted chronological-age page based on a frozen, independently reproduced model.
- Add independently curated research-evidence summaries to every detail page, including an explicit `no_supported_targeted_intervention` result when that is the honest answer.
- Keep exact probiotic strain/formulation, diet protocol, compound identity, population, outcome, dose, duration, and adverse evidence attached to the study that supports them.
- Put all complete visual readings before the long explanations.
- Make every metric understandable at a glance without falsely treating every high value as bad or every low value as good.
- Fail closed for minors and for missing safety-critical metadata.
- Separate an experimental FMT research prescreen from clinical donor eligibility and from FMT treatment evidence.

### 0.2 What v04 must not claim

- No “disease detected,” “risk percent,” “you have,” “rules out,” or “diagnostic” language from a pattern-concordance result.
- No “healthy gut checklist” or claim that every healthy person must carry a named species.
- No unqualified absolute “overgrowth,” “depletion,” or “eradication” from relative abundance alone. “Relative expansion/depletion versus the named reference” is permitted.
- No pathogenic label for species-level *Escherichia coli*, *Klebsiella pneumoniae*, or another context-dependent organism without a validated pathotype, toxin, strain, virulence constellation, or clinical test.
- No vitamin deficiency, neurotransmitter level, metabolite concentration, or pathway activity inferred from DNA capacity.
- No generic-species transfer of strain-specific probiotic evidence.
- No ingredient-level claims split out of a studied multi-ingredient product.
- No personalized dose. Store exact studied protocols in the evidence registry; prescription and FMT doses never appear in participant mode. Other exact doses appear to participants only when a claim-specific release flag confirms that the net impression is not treatment instruction.
- No non-CDI FMT recommendation outside a specific clinical investigation conducted under an effective IND, IRB oversight, and informed consent, unless FDA has determined that an IND is not required, or a future FDA-approved indication. A guideline statement or trial registration alone does not authorize treatment.
- No minor labeled an FMT donor candidate.
- No red/green health judgment for descriptive metrics such as Firmicutes:Bacteroidetes ratio, Shannon diversity, or chronological-age residual unless a validated directional clinical interpretation exists.

---

## 1. Current-system audit and mandatory corrections

The supplied reports were generated by `openbiota v2.0.0`, despite the source specification being named v03. The live pipeline contains 37 Python modules, 71 YAML files, 512 tests, 25 functional panels, 10 microbial groups, 30 profile tasks, 133 study modules, and 109 cited studies.

### 1.1 Version every layer independently

Add these required fields to every run:

```yaml
build_spec_version: 4.0.0
pipeline_version: 2.x.y
observation_schema_version: 3.x.y
finding_schema_version: 4.0.0
disease_profile_schema_version: 3.x.y
intervention_schema_version: 4.0.0
results_schema_version: 4.0.0
report_schema_version: 4.0.0
profile_registry_version: ...
reference_registry_version: ...
intervention_registry_version: ...
registry_lock_digest: sha256:...
```

The pipeline version must not be inferred from a specification filename.

### 1.2 Findings from the three supplied reports

| Finding | SAMPLE2 | A04 | A06 | Required correction |
|---|---:|---:|---:|---|
| GMWI2 | -1.35 | +1.65 | +2.63 | Keep, but show compatible-reference and validation scope |
| Profiles above 75th percentile | 17/28 | 5/28 | 3/28 | Render all; do not cap detailed profiles at eight |
| Shannon | 3.17 | 3.66 | 3.49 | Descriptive reference strip, not good/bad dial |
| Evenness | 0.70 | 0.77 | 0.75 | Same |
| F:B ratio | 0.54 | 3.37 | 2.03 | Neutral/de-emphasized; no universal optimum |
| MetaPhlAn 3 species | 92 | 114 | 105 | Clearly label scoring lane |
| MetaPhlAn 4 SGBs | 164 | 203 | 191 | Clearly label extended, currently unscored lane |

The two reportedly healthy young samples were rendered with adult references while age, sex, and country were recorded as absent. A04 nevertheless showed clinical-Alzheimer 86th, rheumatoid-arthritis 84th, celiac 84th, and alopecia-areata 82nd percentiles. A06 showed ACVD 98th, ankylosing-spondylitis 79th, and IBS-D 78th. These are not valid age-matched screening conclusions. Participant reports for pediatric samples must abstain from adult normative labels until a validated pediatric comparator is available; visibly qualified experimental/OOD numbers may exist only in access-controlled internal research output.

SAMPLE2's Long-COVID profile was 82nd percentile while sixteen other profiles were equal or higher, including several above the 99th percentile. This is a classic `diffuse` profile shape and demonstrates why a high pattern result must not produce disease-specific treatment.

### 1.3 Correct current contradictions and defects

- Generate score-method prose from the same canonical formula/configuration used by code. Current pages conflict on whether `q/s/c` factors multiply the disease score.
- AGA (androgenetic alopecia) must abstain without sex or return a declared two-stratum range; the live reports currently score it despite mandatory sex being missing.
- Implement all five v03 Alzheimer stage tasks. The live reports expose only preclinical amyloid, MCI, and clinical AD; `AD_SCS` and `AD_SCD` are missing.
- Restore pancreatic ductal adenocarcinoma as a separately versioned profile. V03 promised PDAC in its first-wave table, but the live 30-profile inventory and supplied reports omit it.
- Replace “matched reference” when no match criteria were supplied with `unmatched_context_only`.
- If host reads were removed upstream, report `host_fraction: not_assessable_upstream_removed`, never `0.00% pass`.
- Keep MetaPhlAn 3 scoring counts and MetaPhlAn 4 extended-catalog counts visually distinct.
- Remove text truncation and literal ellipses, eliminate orphan pages, embed fonts, add PDF tags/language/bookmarks/clickable citations, and enforce readable type sizes.
- Front-load the complete visual results atlas; do not force readers through function prose before seeing community, taxa, profiles, and interventions.

Add `profiles/expected_profiles_v04.yaml` as the canonical set manifest. It imports the exact IDs and digests of all 30 live v03 tasks and adds only these three tasks:

```yaml
profile_set_id: OPENBIOTA_PROFILE_SET_V04
expected_count: 33
inherits:
  profile_set_id: OPENBIOTA_LIVE_V03_30
  literal_profile_ids_file: profiles/live_v03_30.lock.yaml
  exact_profile_ids_digest: sha256:...
additions:
  - AD_SCS_METAAD_V1
  - AD_SCD_METAAD_V1
  - PDAC_STOOL_MULTINATIONAL_V1
```

`profiles/live_v03_30.lock.yaml` contains the sorted literal IDs, versions, and content digests exported from the currently compiled live registry; the digest alone is insufficient. Build and report generation require set equality among this manifest, compiled profile files, actionability mappings, validation fixtures, `results.json`, and PDF rows. A study module cannot appear as a standalone profile unless separately declared. Every v03-mentioned task must be marked `shipped`, `pending_blocked`, or `retired_with_reason`; silent omission is a compiler failure.

---

## 2. Architecture

### 2.1 Required data flow

```text
FASTQ + manifest + optional external assays
        |
        v
QC -> typed observations -> compatible reference calls -> findings
                    |                         |
                    v                         v
         v03 disease profiles        interpretation assertions
                    |                         |
                    +------------+------------+
                                 v
                    intervention evidence retrieval
                                 v
                  applicability + safety-rule gates
                                 v
          evidence summaries / no-supported-evidence state
                                 |
                                 v
                    results.json -> PDF/text report
```

Measurement, interpretation, disease resemblance, and intervention evidence remain individually inspectable in `results.json` and in the report.

### 2.2 Preserve v03 and add v04 entities

Keep v03 `Assay`, `Observation`, `EvidenceAssertion`, `SignatureModule`, `Profile`, `CompatibilityContainer`, typed adapters, availability masks, module scoring, fusion, coverage, confidence vector, and shape analysis unchanged.

Add:

| Entity | Purpose |
|---|---|
| `ReferenceDefinition` | Frozen compatible comparator, population, assay space, estimators, and fallback policy |
| `ReferenceCall` | Detection state and reference position for one observation |
| `Finding` | Reportable high, low, expected-but-not-detected, uncommon, or concern-qualified result |
| `InterpretationAssertion` | Literature-backed meaning of a finding, independent of treatment |
| `InterventionIdentity` | Stable composition, strain/compound/product identity, formulation, and product version—never the administered regimen |
| `InterventionProtocol` | Exact administered amount, serving, route, frequency, duration, indication/population, and co-interventions studied |
| `InterventionEvidenceAssertion` | One PICO-specific result from one study or guideline, including null/harm |
| `SafetyRule` | Deterministic contraindication, restriction, warning, or suppression rule |
| `InterventionEvidenceResult` | Evidence and applicability for this participant/context; never predicted benefit |
| `ExposureEvent` | Time-stamped medication, supplement, diet, infection, or treatment exposure |
| `LongitudinalResult` | Compatible within-person change with uncertainty and confounders |
| `AgeModelBundle` | Frozen preprocessing, feature namespace, model, calibration, and validation artifacts |
| `ResearchDonorPrescreen` | Access-controlled, experimental, non-clearing FMT research summary for researcher or credentialed-clinician mode only |

### 2.3 Repository additions

```text
openbiota/
├── schemas/v4/
├── references/
│   ├── definitions/
│   ├── manifests/
│   ├── detection_models/
│   └── distributions/
├── findings/
├── interpretations/
├── interventions/
│   ├── identities/
│   ├── protocols/
│   ├── assertions/
│   ├── guidelines/
│   └── source_snapshots/
├── safety_rules/
├── outcome_ontology/
├── age_models/
├── longitudinal/
├── donor_prescreen/
├── reports/v4/
└── registry.lock
```

No YAML is trusted at runtime unless its schema validates and its digest is present in `registry.lock`.

### 2.4 Live-repository change map

Keep the current deterministic Python/YAML design and add modules rather than rewriting the working FASTQ pipeline.

| Current area | v04 change |
|---|---|
| `openbiota/preprocess.py` | Add explicit upstream-host-removal state, specimen-level QC artifact, control/batch fields |
| `openbiota/refcohort.py` | Split reference definition, stratum resolver, prevalence posterior, detection model, study bootstrap, OOD |
| `openbiota/scoring.py` | Preserve v03 scores; add typed reference-call engine without changing disease math |
| `openbiota/community.py` | Add neutral ecology semantics and global multivariate unusualness |
| `openbiota/taxongroups.py` | Add group-specific context; split mucin degraders; rename health-associated core |
| `openbiota/profiles.py`, `similarity.py`, `shape.py` | Preserve v03; enforce missing-metadata abstention/fallback and full detail count |
| `openbiota/engines/metaphlan.py` | Preserve MetaPhlAn 3 scoring lane and exact namespace |
| `openbiota/engines/metaphlan4.py` | Preserve extended inventory; score only after its own rebuilt reference passes validation |
| New `openbiota/detection.py` | Feature-specific downsampling curves, LOD/LOQ states, nondetection qualification |
| New `openbiota/findings.py` | Independent detection/reference/interpretation axes and FDR/stability promotion |
| New `openbiota/interventions.py` | Exact identity/protocol/PICO retrieval and evidence-result records |
| New `openbiota/safety.py` | Restricted rule compiler and fail-closed evaluation |
| New `openbiota/age.py` | Frozen age-bundle inference, interval, OOD, attribution, abstention |
| New `openbiota/longitudinal.py` | Compatibility, exposure timeline, paired delta/uncertainty |
| New `openbiota/donor_prescreen.py` | Access-controlled, non-clearing research summary and assay-gap inventory; participant serializer has no prescreen access |
| `pdfreport.py`, `pdfsummary.py`, `pdflibrary.py`, `pdfprofiles.py` | Recompose front atlas, typed visuals, unlimited detail pages, evidence cards, accessible PDF |

Required commands:

```text
openbiota run --manifest sample.yaml --report-schema 4 --mode research
openbiota validate-registry --lock registry.lock
openbiota validate-reference REFERENCE_ID
openbiota validate-interventions --include-negative-evidence
openbiota age-model audit AGE_MODEL_ID
openbiota report-audit results/SAMPLE/report.pdf results/SAMPLE/results.json
make references-v4
make interventions-v4
make age-model-v4
make test-v4
```

No language model runs during production analysis or report generation. All report text is deterministic from schema-validated, curator-approved claim records.

### 2.5 Heterogeneous evidence graph

Retain v03 typed adapters and add a declarative expression DAG with these operators:

```yaml
operators:
  - directional_abundance
  - presence
  - qualified_nondetection
  - percentile_range
  - sequence_allele
  - strain_match
  - virulence_constellation
  - gene_abundance
  - pathway_abundance
  - pathway_completeness
  - ecological_metric
  - external_measurement
  - predicted_flux
  - frozen_model_output
  - all
  - any
  - at_least_k
  - independent_study_fusion
```

Every leaf declares exact evidence type, feature namespace/version, body site, unit, tool/database, and adapter. Every transport is explicit. Missing evidence propagates as missing, never zero. Same-participant modalities remain one evidence cluster rather than becoming false replication.

---

## 3. Manifest, metadata, and safety intake

### 3.1 Required specimen manifest

```yaml
sample:
  sample_id: ...
  subject_id: ...
  collection_datetime: ...
  age_years: ...
  date_of_birth_precision: exact | year_only | unknown
  minor_status: true | false | unknown
  sex_at_birth: female | male | intersex | unknown
  country: ...
  body_site: stool
  stool_form_bristol: 1-7 | unknown
  preservative: ...
  freeze_delay_hours: ...
  extraction_kit: ...
  library_kit: ...
  platform: ...
  read_layout: paired
  read_length: ...
  sequencing_run: ...
  flowcell: ...
  lane: ...
```

Chronological age is not an input feature to the microbiome-age predictor. It is required for age-dependent references, age residuals, pediatric/safety logic, and validated-stratum applicability. Missing age never defaults to adult; participant-mode numeric age prediction also fails closed until age assurance can rule in the applicable minor/adult release policy. Sex is required for sex-stratified profiles such as androgenetic alopecia. Country, stool form, medications, diet, and collection factors are required when the selected reference declares them mandatory.

### 3.2 Required clinical and safety-context fields

```yaml
context:
  age_years: decimal | null
  age_source: date_of_birth | verified_record | self_report | guardian_report | unknown
  minor_status: true | false | unknown
  preterm_infant: true | false | not_applicable | unknown
  pregnancy: true | false | unknown
  breastfeeding: true | false | unknown
  immunocompromise:
    status: true | false | unknown
    severity: mild | moderate | severe | unknown | not_applicable
    causes: [neutropenia | chemotherapy | hematologic_malignancy | primary_immunodeficiency | immunosuppressive_drug | other]
    absolute_neutrophil_count: {value: decimal | null, unit: cells_per_uL, measured_at: date | null}
  central_venous_catheter: true | false | unknown
  critical_illness: true | false | unknown
  recent_hospitalization: true | false | unknown
  fever: true | false | unknown
  suspected_or_confirmed_acute_infection: true | false | unknown
  gi_bleeding: true | false | unknown
  renal:
    disease: true | false | unknown
    stage: string | null
    egfr: {value: decimal | null, unit: mL_min_1_73m2, measured_at: date | null}
    potassium: {value: decimal | null, unit: mmol_L, measured_at: date | null}
    phosphorus: {value: decimal | null, unit: mg_dL, measured_at: date | null}
  hepatic:
    disease: true | false | unknown
    cirrhosis: true | false | unknown
    decompensated: true | false | unknown
  transplant:
    solid_organ_or_stem_cell: true | false | unknown
    type: string | null
    date: date | null
  active_ibd:
    status: true | false | unknown
    severity: mild | moderate | severe | unknown | not_applicable
  obstruction_stricture_or_short_bowel: true | false | unknown
  eating_disorder_or_malnutrition_risk: true | false | unknown
  allergies: []
  medications:
    - {rxnorm_id: string | null, name: string, dose: decimal | null, unit: string | null, frequency: string | null, route: string | null}
  supplements: []
  antibiotics_last_180_days: []
  ppi: true | false | unknown
  metformin: true | false | unknown
  anticoagulant_or_antiplatelet: true | false | unknown
  diabetes_medication: true | false | unknown
  planned_surgery: {status: true | false | unknown, date: date | null}
  clinician_confirmed_conditions:
    - condition_code: SNOMEDCT:...
      confirmation_source_type: clinician_attestation | ehr_record | accredited_lab_result
      source_record_ref: ...
      confirming_organization_or_clinician_ref: ...
      confirmation_date: date
      credential_or_attestation_ref: ...
      verification_status: verified | pending | rejected
  self_reported_conditions: []
  symptoms: []
```

Unknown safety-critical fields are not treated as false. Under the versioned rule they set `known_safety_signal: unknown_incomplete_context`, then set `display_policy: clinician_only | evidence_only | suppressed` and `action_status: evidence_only | suppressed` according to intervention class and the missing field. Retired combined statuses such as `clinician_review_required` and `evidence_only_not_personalized` are schema errors and cannot compile.

Participant-entered or otherwise self-reported conditions remain `self_reported`; they cannot be promoted to `verified`, satisfy a confirmed-indication selector, or unlock prescription or FMT routes. Confirmation provenance is immutable and independently auditable.

`sample.age_years` and its provenance are canonicalized once into `context.age_years`; `context` is the only safety-rule runtime namespace. The safety compiler rejects `subject.*`, unresolved aliases, unit mismatches, and rules without explicit `on_unknown`. `minor_status` is derived from a verified age when available; a conflicting supplied flag is a hard intake error, never resolved silently.

### 3.3 Exposure history

The SAMPLE2 supplement history is valuable longitudinal context but is not treatment evidence. Store every exposure as an event:

```yaml
exposure_event:
  event_id: EXP_...
  class: antibiotic | prescription | botanical | probiotic | fermented_food | diet | infection | hospitalization | FMT
  exact_identity_or_label: ...
  product_lot_if_known: ...
  start_date: ...
  stop_date: ...
  dose_as_reported: ...
  adherence: unknown | low | moderate | high
  source: subject_report | clinician_record | medication_record
```

Berberine, allicin, black-seed oil, oregano oil, GI Microb-X, and high-dose probiotics must be recorded as separate or exact-combination exposures. A post-exposure nondetection may be reported as temporal association only. It must never be called eradication or proof of efficacy without an appropriate absolute/clinical test and a causal design.

---

## 4. Analytical QC and detection contract

### 4.1 Sample-level gates

Retain current fastp, pairing, host filtering, classifier, and `rpoB` controls. Record at minimum:

- input and post-trim pairs;
- usable non-host pairs;
- Q20/Q30, GC, duplication, read-length distribution, insert size, adapter content;
- host fraction status, including `not_assessable_upstream_removed`;
- taxonomic classified fraction by profiler lane;
- `rpoB` fragment count and functional saturation;
- negative-control, contamination, and batch status when available;
- database/tool/reference fingerprints.

The existing minimum of 500,000 usable non-host pairs remains a sample-level floor, but taxon nondetection also requires a feature-specific detection model.

### 4.2 Detection states

Every feature returns one:

```yaml
detection_status:
  detected | not_detected | indeterminate | not_measured
```

`not_detected` means the feature was searchable and the sample had adequate feature-specific power. `indeterminate` means an absence claim is not analytically supported. `not_measured` means the modality was never assayed, such as a plasma metabolite in a FASTQ-only run.

### 4.3 Feature-specific detection power

For feature `f`, estimate detection probability as a function of usable microbial reads `N`, read length, profiler/database, ambiguity, and the positive-abundance distribution:

\[
D_f(N)=\int P(\text{detected}\mid N,a)\,dF_{f,+}(a)
\]

Build curves by empirical downsampling of positive reference specimens and synthetic-spike controls. A nondetection is qualified only when:

- specimen QC passes;
- `f` exists in the exact database/namespace;
- expected sensitivity at this sample depth is at least 0.95;
- mapping ambiguity and contamination gates pass;
- the matching reference itself had comparable searchability.

The report must display the estimated detection power and limit context for any “commonly detected but not detected” call.

---

## 5. Matched healthy/control reference system

### 5.1 Why there is no universal healthy-species checklist

Healthy gut taxa vary strongly across age, geography, stool consistency, diet, medication, extraction, sequencing, and database. Function can be more conserved than taxonomy. Therefore “missing” is a statistical statement about a frozen compatible reference—not a claim that a person lacks an essential organism.

Use this user-facing phrase:

> Commonly detected in the matched reference, but not detected in this specimen at the stated sequencing depth.

Never use “required organism,” “deficient species,” or “this organism must be replaced.”

### 5.2 Independent axes for every feature

```yaml
detection_status:
  detected | not_detected | indeterminate | not_measured

match_status:
  exact | validated_fallback | unmatched_context_only | unavailable

reference_position:
  very_low | low | typical | high | very_high |
  uncommon_detected | commonly_detected_not_detected |
  comparator_unavailable

interpretation:
  summary: health_associated | disease_associated | potential_pathobiont |
    validated_pathogen_signal | context_dependent |
    conflicting_evidence | no_curated_interpretation
  assertion_refs: []
  tags: []

analytical_qualification:
  overall: pass | warn | fail | not_applicable
  flags: [low_depth | ambiguous_mapping | strain_unresolved |
    compositional_only | batch_mismatch | out_of_distribution]

population_transportability: 0.0-1.0 | unknown
applicability_reasons: []
```

Level and interpretation are separate. A high level may be favorable, concerning, neutral, or unknown.

### 5.3 Reference definition

```yaml
reference_model_id: EXAMPLE_HEALTHY_STOOL_MPA3_ADULT_V1
version: 1.0.0
status: design_example_not_releaseable
body_site: stool
compatible_observation:
  evidence_type: TAX_REL
  namespace: metaphlan_species
  namespace_version: mpa_v31_CHOCOPhlAn_201901
  unit: relative_abundance
  transform: fixed_reference_logratio_v1
  transform_state_digest: sha256:...
cohort:
  manifest_digest: sha256:...
  release_and_query_digest: sha256:...
  inclusion_definition: versioned_healthy_control_predicate
  exclusions:
    - antibiotics_within_90_days
    - acute_gastrointestinal_infection
    - duplicate_or_nonbaseline_specimen
  n_participants: 3027
  n_studies: 22
strata:
  mandatory: [age_band]
  preferred: [country, sex, stool_form]
  minimum_n_for_percentile: 100
  minimum_positive_n_for_abundance: 50
  minimum_independent_studies: 3
  maximum_percentile_interval_width: 30
fallback_ladder:
  - [age_band, country, sex]
  - [age_band, country]
  - [age_band, validated_transport_only]
  - unmatched_context_only
estimators:
  prevalence: study_balanced_hierarchical_prevalence_v1
  exploratory_prevalence_fallback: jeffreys_beta_posterior
  positive_abundance: study_balanced_empirical_cdf
  detection: feature_depth_detection_curve_v1
  uncertainty: study_block_bootstrap
  multiplicity_policy_ref: MULTIPLICITY_TAX_MPA3_V1
```

Disease-profile fallback and expected/missing rules are intentionally different. A disease profile may return a visibly low-transportability internal research score against a broader reference. A strong expected-but-not-detected label requires at least 100 matched participants, three independent study blocks, feature-specific held-out-country/stool-form transportability, and a stable detection model; otherwise show only `unmatched_context_only` or `comparator_unavailable`.

The 3,027-participant/22-study and 34-adult counts describe the current legacy seed, not a releaseable reference. Replace `EXAMPLE_*` only when a checked-in participant/sample/run manifest, exact cMD release and query, disease/BMI/pregnancy/hospitalization/medication eligibility, antibiotic window, baseline and deduplication rules, preprocessing/database hashes, and expected counts all validate. A count/manifest mismatch fails the build.

The frozen univariate taxon transform is not RCLR. For a reference-learned stable denominator set \(D\), frozen weights \(w_j\), and feature-specific reference/LOD replacement \(\delta_j\):

\[
g_i=\exp\!\left(\sum_{j\in D} w_j\log(x_{ij}+\delta_j)\right),\qquad
z_{if}=\log\frac{x_{if}+\delta_f}{g_i}.
\]

The bundle freezes the eligible feature universe, \(D\), weights, every \(\delta\), minimum observed denominator fraction, training-only fit manifest, and state digest. If denominator support is below its frozen gate, the comparison is unavailable. RCLR/RPCA may be used separately for multivariate distance/OOD only and must be fitted within training folds.

### 5.4 Prevalence classes

The primary prevalence estimator is a study-balanced hierarchical binomial/logistic model with a study random effect, detection-power adjustment, recorded study-specific prevalence/heterogeneity, and leave-one-study-out stability. No single large cohort may dominate. Where that model is unsupported, a pooled Jeffreys posterior may be shown as exploratory context only:

\[
\pi_f\mid k,n\sim\operatorname{Beta}(k+0.5,n-k+0.5)
\]

- first, `reference_core` if \(P(\pi_f\ge0.80)\ge0.95\);
- else `reference_common` if \(P(\pi_f\ge0.50)\ge0.95\);
- else `reference_rare` if \(P(\pi_f\le0.05)\ge0.95\);
- else `reference_variable`.

These are mutually exclusive by ordered precedence. Core/common promotion must also survive leave-one-study-out and the declared minimum fraction of contributing studies.

A qualified nondetection of a core/common feature becomes `commonly_detected_not_detected`. It does not become a treatment target automatically.

### 5.5 High/low calls and multiplicity

For detected taxa that pass the frozen quantitative gate at or above LOQ, use the frozen fixed-reference zero-aware log-ratio above, not raw relative abundance alone. A nondetected taxon is classified through `TAX_PRES` first and cannot also receive a low/depleted abundance headline; a detected value below LOQ or with no established LOQ may retain a secondary descriptive value but receives no low/high percentile. Gene/pathway records use a separate two-part hurdle: detection/prevalence plus positive relative encoded-potential abundance in a frozen compatible HUMAnN/gene-family namespace and frozen compositional denominator. `PATH_COMPL` reaction/gene-set completeness is separate from HUMAnN pathway abundance/coverage. Never call either absolute activity or metabolite production. For every eligible feature:

1. calculate a study-balanced empirical percentile and interval;
2. calculate a two-sided empirical tail probability;
3. adjust using Benjamini-Hochberg within a frozen `multiplicity_family_id` that records evidence type, namespace, eligible feature-universe digest, inclusion/exclusion rules, tail construction, and nondetection policy;
4. repeat with study-block bootstrap and leave-one-study-out references;
5. report all values, but promote a headline only if:

```yaml
outside_central_90_percent: true
bh_q_lte: 0.10
direction_stability_gte: 0.80
technical_qc: pass
```

Abundance positioning additionally requires at least 50 positive reference observations from at least three study blocks, sufficient discrete percentile resolution, and an interval-width pass; otherwise prevalence may render but the abundance comparator is unavailable. Validate both per-feature q values and per-report false-flag burden on held-out controls. `fdr_q` is reference-extremeness control, not clinical significance or disease evidence; display filtering must never change the testing family.

Call taxa `relative_expansion` or `relative_depletion`. Use “overgrowth” only if a compatible absolute-load assay supports it.

### 5.6 Pediatric reference gate

The adult 3,027-sample taxonomic reference and 34-adult functional reference are not valid pediatric norms. Build separate pediatric strata from public data only after harmonized reprocessing and governance checks. Initial anchors include:

- Hollister et al., ages 7–12, 37 16S and 22 WGS specimens: https://doi.org/10.1186/s40168-015-0101-x
- Laue et al., 247 shotgun specimens including 167 healthy adolescents ages 11–14: https://doi.org/10.1186/s13073-025-01481-1

These are seed cohorts, not sufficient proof of a universal pediatric reference. Required release gates are at least 100 participants and three independent studies per displayed age stratum, held-out-study stability, compatible preprocessing, and prespecified metadata coverage. Until then:

```yaml
reference_position: comparator_unavailable
population_transportability: unknown
applicability_reasons: [pediatric_reference_insufficient]
```

No adult `healthy range`, low/high, expected/missing, or intervention applicability label may be silently applied to an 8–17-year-old.

---

## 6. Finding ontology and schemas

### 6.1 User-facing organism buckets

Every reference-eligible organism belongs to exactly one primary display bucket:

1. Detected in the typical reference range
2. Relatively high versus matched reference
3. Relatively low versus matched reference
4. Commonly detected in matched references but not detected here
5. Uncommon in matched references and detected here
6. Potential-concern signal requiring strain/toxin/clinical qualification
7. Indeterminate or reference unavailable

An organism may carry secondary evidence tags, but it cannot be simultaneously counted in multiple headline buckets.

### 6.2 Taxon finding schema

```yaml
finding_id: FINDING_SAMPLE123_TAX_0042
finding_schema_version: 4.0.0
sample_id: SAMPLE123
feature:
  evidence_type: TAX_REL
  namespace: metaphlan_species
  namespace_version: mpa_v31_CHOCOPhlAn_201901
  feature_id: metaphlan_mpa3_clade_stable_id
  display_name: Escherichia coli
  full_lineage: k__Bacteria|...
  body_site: stool
detection:
  status: detected
  identity_confidence: null
  identity_gate_status: pass | warn | fail | not_established
  identity_gate_method: frozen_marker_breadth_depth_ambiguity_v1
  identity_calibration_bundle_id: IDENTITY_MPA3_MOCKSPIKE_V1 | null
  identity_calibration_digest: sha256:... | null
  marker_policy_digest: sha256:...
  identity_thresholds: {minimum_breadth: ..., minimum_depth: ..., maximum_ambiguity: ...}
  detection_power_if_present: null
  marker_or_read_breadth: 0.91
  marker_or_read_depth: 22.4
  ambiguity_fraction: ...
  lod95: {value: 0.002, unit: relative_abundance_percent, method: empirical_downsampling_v1}
  loq: {value: null, unit: relative_abundance_percent, method: not_established}
  quantitation_gate: loq_not_established
  quantitative_status: relative_only
  quantitative_status_reason: no_absolute_load_assay
reference_comparison:
  reference_model_id: EXAMPLE_HEALTHY_STOOL_MPA3_ADULT_V1
  match_status: unavailable
  percentile: null
  interval: null
  position: comparator_unavailable
  fdr_q: null
  leave_one_study_out_direction_stability: null
  measurement_basis: relative_abundance
  absolute_load_available: false
interpretation:
  label: detected_context_dependent
  polarity: context_dependent
  pathogen_status: unresolved_at_species_level
  assertion_refs: [INTERPRET_ECOLI_CONTEXT_V1]
report_language:
  preferred: "E. coli was detected; a releaseable matched comparator was unavailable."
  prohibited:
    - "E. coli overgrowth"
    - "pathogenic E. coli"
    - "E. coli infection"
    - "E. coli eradication"
```

`identity_confidence` is nullable and may be populated only by a locked, independently calibrated probability model. MetaPhlAn output is not itself such a probability. The deterministic identity gate freezes database/marker policy, mock/spike calibration artifacts, ambiguity handling, and breadth/depth thresholds. Unknown or unsupported identity blocks uncommon-present, potential-concern, strain, virulence, AMR, and intervention promotion. Every `reference_model_id` is a compile-time foreign key, and a missing or `design_example_not_releaseable` record cannot produce `match_status: exact` or numeric comparison fields.

### 6.3 Feature-type semantics

| Evidence | Permitted finding |
|---|---|
| `TAX_REL` | Relatively low/high versus compatible reference |
| `TAX_PRES` | Qualified commonly-detected nondetection or uncommon detection |
| `GENE_ABUND` | Lower/higher genomic capacity |
| `PATH_COMPL` | Pathway capacity incomplete or not detected |
| `STRAIN` | Strain not detected/detected at qualifying coverage |
| `AMR` | Resistance determinant detected; phenotype not established |
| `METAB` | Measured low/high only when the external assay exists |
| `FLUX` | Predicted lower/higher under the declared diet scenario |
| `ECO` | Below/above reference; clinical polarity may remain neutral |
| `MODEL` | Frozen model output with OOD and calibration status |

### 6.4 Guild encoded potential and redundancy before single-species missingness

Add encoded-potential/redundancy summaries for butyrate production, resistant-starch degradation, mucin utilization, bile-acid transformations, methanogenesis, sulfate reduction, oxalate degradation, and other supported guilds. Freeze required/optional gene and reaction sets, copy-number normalization, taxon-stratified attribution, coverage/ambiguity gates, missing-data behavior, and database versions. This is not metabolic completeness, expression, production, or flux. A low or absent single carrier is secondary to whether alternate encoded carriers are present. Contradictions are first-class findings—for example SAMPLE2's high butyrate gene capacity but 4th-percentile butyrate-producer abundance must render as an explicit cross-engine contradiction rather than a single conclusion.

### 6.5 Every species gets an explicit evidence state

The software does not need a fabricated treatment for every species; it needs a complete, auditable disposition for every species/SGB in both profiler lanes. Build a `SpeciesInterpretation` registry keyed by exact namespace and version:

```yaml
species_interpretation:
  namespace: metaphlan_species
  namespace_version: mpa_v31_CHOCOPhlAn_201901
  feature_id: ...
  display_name: ...
  taxonomy_mappings:
    - source_namespace: metaphlan_species
      source_release: ...
      source_id: ...
      target_namespace: NCBITaxon | GTDB | UHGG | SGB
      target_release: ...
      target_ids: [...]
      relation: exact | synonym | merged | split | ambiguous
      confidence: 0.0-1.0 | unknown
      curator_status: reviewed | pending
  body_site: stool
  reference_models: [...]
  guild_memberships: [...]
  interpretation_assertions: [...]
  disease_assertions: [...]
  virulence_pathotype_requirements: [...]
  intervention_assertions: [...]
  evidence_state:
    curated_human_action_evidence |
    human_association_only |
    mechanistic_only |
    no_curated_interpretation |
  last_reviewed: ...
```

For each detected or reference-eligible taxon the report shows abundance/detection, matched prevalence/percentile, analytical confidence, guild role, association assertions with direction conflicts, and either exact intervention evidence or `no_supported_targeted_intervention`. Unclassified SGBs remain named by stable SGB ID and never inherit an interpretation from a loose name match.

Namespace completeness applies to machine-readable JSON, not one PDF page per database feature. The PDF includes every detected taxon, qualified core/common nondetection, uncommon detection, concern signal, and explicitly requested taxon. Variable/rare nondetections with no interpretation are summarized by counts and remain available in a machine-readable appendix/export. Evidence never transfers across a split/ambiguous taxonomy mapping unless the assertion's required resolution is demonstrably preserved.

Seed discovery from curatedMetagenomicData/GMrepo, the Human Microbiome Project, UHGG/GTDB/NCBI taxonomy, primary disease studies, and the curated v03 assertions. Do not automatically turn text-mined associations into human-facing claims. Every assertion needs curator review, source digest, body site, assay, population, and direction.

### 6.6 Clinical-concern and FMT-exclusion signals need a separate lane

Community profiling is not a clinical pathogen test. Add a separately validated lane for strain/pathotype, toxin, virulence, and AMR signals, with assembly or read-level breadth/depth rules and explicit host-linkage uncertainty. Examples include Shiga-toxin/pathotype genes for *E. coli*, toxin genes for *C. difficile*, enterotoxigenic *B. fragilis*, pks/colibactin loci, hypervirulence/AMR context for *K. pneumoniae*, and clinically relevant enteric viruses/parasites.

Required states:

```yaml
clinical_concern_signal:
  detected_confirmatory_testing_required |
  no_signal_detected_at_this_depth |
  indeterminate |
  modality_not_assessed
```

Short-read AMR detection does not necessarily identify the microbial host or phenotypic susceptibility. Antibiotic choice requires infection site, isolate/organism, susceptibility, allergies, organ function, and clinician judgment. Follow [IDSA infectious-diarrhea guidance](https://www.idsociety.org/practice-guideline/infectious-diarrhea/) and [CDC *E. coli* guidance](https://www.cdc.gov/ecoli/hcp/guidance/index.html); possible STEC must not generate an antimicrobial or antimotility suggestion.

For FMT research prescreening, the report lists bacterial, viral, eukaryotic/parasitic, toxin, and AMR/MDRO lanes as `assessed`, `partially_assessed`, or `not_assessed`. No collection of negative shotgun results becomes clinical clearance.

---

## 7. Intervention evidence system

### 7.1 Intervention evidence is not disease evidence

Do not auto-convert any v03 disease association into a treatment assertion. Association, causal target, and intervention benefit are three different propositions. Every intervention card needs its own source and population-intervention-comparator-outcome record.

For each detailed finding, the engine must render one or more of:

- Guideline/label context after independent clinical confirmation
- Human research evidence summary
- Mechanistic evidence only
- Evidence against or a safety warning
- No supported targeted intervention found in registry version X

The last state is a successful result, not an engine failure.

### 7.2 Exact intervention identity

| Class | Required identity fields |
|---|---|
| Probiotic | Every strain, stable accession/culture collection ID, exact combination, delivery/formulation, and product version. A marketed fixed strength may be stored as label metadata, but the amount actually administered belongs only to the protocol. |
| Prebiotic | Molecule, chain length/degree of polymerization, purity, formulation |
| Fermented food | Food/product identity, fermentation method, starter/live-culture status, pasteurization, and version; serving belongs to the protocol |
| Botanical/supplement | Latin species, plant part, extraction solvent/ratio, standardized constituents, chemical form |
| Diet | Named diet identity/version; food targets, exclusions, energy constraints, and duration belong to the protocol |
| Prescription drug | RxNorm identity, stable formulation, and product version; route and confirmed indication belong to the protocol/PICO |
| FMT/microbiota product | Product/donor-source identity, processing/formulation, donor-screening standard, and version; route and schedule belong to the protocol |
| Multi-ingredient product | Exact complete stable composition and formulation; evidence cannot be split across ingredients. Administered amounts belong to the protocol. |

```yaml
intervention_id: PROBIOTIC_EXACT_COMBINATION_001
version: 1.0.0
class: probiotic
identity:
  evidence_transfer_key: exact_strain_combination
  components:
    - entity_type: probiotic_strain
      taxon_id: NCBITaxon:...
      strain_designation: ...
      culture_collection_ids: [...]
  formulation:
    delivery: capsule
    coating: specified | unknown
  product_version: ...
  label_declared_strength_metadata: {value: ..., unit: ..., basis: ...} | null
transfer_policy:
  exact_protocol: permitted
  same_strain_different_dose: sensitivity_only
  same_species_other_strain: prohibited
  genus_only: prohibited
```

### 7.3 Protocol and PICO assertion

```yaml
protocol_id: PROTOCOL_...
intervention_id: ...
population:
  age_range_years: [...]
  confirmed_condition: ...
  geography: ...
  inclusion_criteria: [...]
  exclusions: [...]
administration:
  dose_source_text: ...
  normalized_dose: ...
  serving_definition: ...
  route: ...
  frequency: ...
  duration_days: ...
  cointerventions: [...]
protocol_scope:
  studied_protocol_not_personalized_instruction: true
```

```yaml
evidence_assertion_id: IE_...
source:
  source_type: randomized_trial | nonrandomized_study | observational_study |
    case_report | systematic_review | meta_analysis | guideline |
    regulatory_label | regulatory_safety_communication | protocol | preprint
  doi: ...
  pmid: ...
  trial_registration: ...
  publication_status: peer_reviewed | preprint | abstract | issued_guidance |
    current_label | withdrawn | corrected | retracted
  correction_or_retraction_status: ...
  snapshot_date: ...
  digest: sha256:...
study:
  design: ... | null
  participants_randomized: ... | null
  analysis_population: ... | null
  risk_of_bias: low | some_concerns | high | not_assessed
  funding_and_conflicts: ...
pico:
  population_ref: ...
  protocol_ref: ...
  comparator: ...
  outcome:
    ontology_id: ...
    class: patient_clinical | patient_reported | validated_biomarker | microbiome_surrogate | mechanistic
    timepoint_days: ...
  effect:
    direction: benefit | null | harm | mixed | not_estimable
    measure: ...
    estimate: ...
    confidence_interval: [...]
    multiplicity_handling: ...
target_binding:
  trigger_type: ...
  required_condition: ...
  microbiome_profile_alone_sufficient: false
  feature_selectors: [...]
```

Null, contradictory, and harmful studies are mandatory registry records. Publication status and conflicts must be visible.

Schema conditionals make study-only fields nullable for guidelines, labels, and safety communications while requiring issuing body, jurisdiction, effective date, recommendation/indication text, and supersession status for those source types.

### 7.4 Evidence lanes

| Lane | Meaning |
|---|---|
| A | Current authoritative guideline or regulatory label for an independently confirmed indication |
| B | Replicated human randomized clinical benefit for the exact intervention/population |
| C | Single/pilot/conflicting human randomized evidence, including clinical or surrogate endpoints |
| D | Nonrandomized human, pilot, case series, or observational evidence |
| E | Animal, in-vitro, or mechanistic evidence only |
| X | Against-evidence record: guideline against, null/conflicting evidence, or harm signal; personal safety is stored separately |

Also retain separate, non-collapsed dimensions:

```yaml
guideline_position:
  strong_for | conditional_for | no_recommendation |
  trial_only | conditional_against | strong_against | unavailable
outcome_relevance:
  clinical | patient_reported | validated_biomarker |
  microbiome_surrogate | mechanistic
replication:
  replicated | single_study | conflicting | none
applicability:
  exact | partial | outside_population | unknown
supporting_evidence_lane: A | B | C | D | E | none
against_evidence: [] | [guideline_against, null_primary, harm_signal, conflicting]
```

`against_evidence` is a deterministic set: union every applicable active assertion, preserve all assertion references, sort by the frozen enum order above, and never discard a null, harm, guideline-against, or conflict record because a supportive record exists. `supporting_evidence_lane` uses only the A–E lane enum; a separately typed `strongest_study_design` may describe design without overloading the lane field.

Never calculate a “treatment score” or percent chance that an intervention will help. Participant evidence summaries are unranked. Internal/clinician evidence review may sort lexicographically by safety, guideline position, clinical versus surrogate endpoint, exact target, exact population, replication, protocol reproducibility, then recency; this ordering must not appear as personalized benefit ranking.

### 7.5 Trigger types

```yaml
trigger_type:
  - clinician_confirmed_condition
  - symptom_phenotype
  - external_lab_abnormality
  - exact_microbiome_feature
  - functional_capacity_finding
  - ecological_finding
  - research_profile_resemblance
```

Materialize all triggers as a typed union rather than pretending every trigger is a microbial finding:

```yaml
trigger_record:
  trigger_id: ...
  trigger_type: reference_finding | disease_profile_result | external_lab_result |
    symptom_phenotype | clinician_confirmed_condition
  source_record_ref: ...
  microbial_context: {body_site: ..., namespace: ..., direction: ...} | null
  confirmation_status: verified | self_reported | profile_only | unknown
```

Body-site, feature-namespace, and microbial-direction checks apply only when `microbial_context` is non-null.

Examples:

- A trial enrolling diagnosed IBS-D requires confirmed IBS-D; an IBS-D resemblance percentile alone is insufficient.
- A trial that changed *Bifidobacterium* but had no clinical benefit belongs in the surrogate lane.
- A probiotic-strain trial cannot be triggered by low stool abundance of that species.
- An antimicrobial cannot be proposed solely because a taxon's relative abundance is high.
- An animal FMT Alzheimer experiment belongs in the mechanistic lane, not a human-treatment lane.
- A disease-profile resemblance may offer clinician follow-up appropriate to symptoms, but not disease medication.

### 7.6 Applicability vector

```yaml
applicability:
  target_match: exact | related | indirect | none
  population_match: exact | partial | outside | unknown
  age_match: exact | outside | unknown
  sex_match: exact | not_required | outside | unknown
  disease_confirmation: confirmed | self_reported | profile_only | absent
  medication_context_match: exact | partial | unknown
  assay_match: native | transported | not_applicable
  endpoint_relevance: clinical | surrogate | mechanistic
  protocol_identity: exact | partial | class_only
  known_dimension_coverage: 0.0-1.0
```

`known_dimension_coverage` is metadata completeness, not strength of evidence or benefit probability.

### 7.7 Candidate-generation algorithm

```text
INPUT: observations, profiles, external measurements, reference calls,
       findings, subject metadata, intervention registry

1. Build Findings only from QC-qualified evidence.
2. Materialize typed `TriggerRecord` values from Findings, disease profiles, external labs, symptoms, and clinician-confirmed conditions; retrieve matching evidence selectors.
3. For microbial triggers only, require exact body site, namespace, claim level, and direction.
4. Bind evidence by exact intervention/protocol equivalence key.
5. Never pool strains, products, formulations, doses, or combinations unless a
   prespecified meta-analysis explicitly permits it.
6. Evaluate population, indication, and confirmatory-test requirements.
7. Retain positive, null, harm, and conflicting assertions.
8. Apply every matching safety rule.
9. If a required diagnosis is absent, downgrade to evidence-only.
10. If safety-critical metadata are unknown, fail closed.
11. Assign evidence lane and applicability vector.
12. Sort transparently; do not predict treatment benefit.
13. If nothing exact survives, emit no_supported_targeted_intervention.
14. Emit every provenance path and rejection reason.
```

### 7.8 Card content

Every internal evidence card includes:

- trigger record(s), and—only for microbial triggers—the intended microbial direction;
- intervention class and exact identity;
- studied protocol as a source fact, never personalized dosing;
- study population and all applicability mismatches;
- clinical and microbiome endpoints shown separately;
- evidence lane, guideline position, replication, risk of bias, and conflicts;
- contraindications, interactions, age/pregnancy/immunocompromise limits;
- conflicts with other findings or options;
- monitoring variables and research retest interval when supported;
- clinician-review requirement;
- direct DOI/PubMed/guideline/trial links;
- source version and retrieval date.

Participant mode renders an unranked `research_evidence_summary`, not an “option” by default. It omits all prescription and FMT doses. It may show an exact food/probiotic/supplement dose only when an immutable claim-level release record documents legal/clinical review, age/population applicability, safety display policy, and why the report's net impression is educational rather than treatment instruction. Otherwise it shows intervention class, exact study identity, population, endpoint, result, uncertainty, and clinician/dietitian follow-up without dosing.

### 7.9 Compiled research-evidence summary output

```json
{
  "summary_id": "EVIDENCE_SUMMARY_SAMPLE123_004",
  "intervention_id": "PROBIOTIC_EXACT_COMBINATION_001",
  "trigger_ids": ["TRIGGER_SAMPLE123_PROFILE_LONGCOVID"],
  "supporting_evidence_lane": "C",
  "against_evidence": [],
  "display_policy": "evidence_only",
  "known_safety_signal": "unknown_incomplete_context",
  "action_status": "suppressed",
  "guideline_position": "no_recommendation",
  "applicability": {
    "target_match": "related",
    "population_match": "unknown",
    "age_match": "exact",
    "disease_confirmation": "profile_only",
    "endpoint_relevance": "clinical",
    "protocol_identity": "exact",
    "known_dimension_coverage": 0.64
  },
  "evidence_summary": {
    "benefit_assertions": 1,
    "null_assertions": 0,
    "harm_assertions": 0,
    "conflicting": false,
    "strongest_study_design": "randomized_controlled_trial"
  },
  "eligibility_explanation": [
    "The study required clinically defined PACS.",
    "This report provides microbiome resemblance only.",
    "Clinical eligibility is therefore unconfirmed."
  ],
  "studied_protocol_refs": ["PROTOCOL_SIM01_PACS_6M_V1"],
  "evidence_assertion_refs": ["IE_SIM01_PACS_RCT_V1"],
  "not_an_effect_probability": true,
  "not_a_personalized_prescription": true
}
```

### 7.10 Evidence curation and update workflow

For every new finding/profile/intervention pair:

1. Run reproducible searches in PubMed, ClinicalTrials.gov/WHO ICTRP as applicable, regulator/product labels, professional guidelines, and systematic-review reference lists.
2. Store the exact query, date, result count, and source snapshot.
3. Deduplicate by DOI, PMID, registry ID, correction, and underlying participant cohort.
4. Screen title/abstract, then full text; record exclusion reasons.
5. Extract PICO, exact identity/protocol, primary versus secondary outcomes, multiplicity, adverse events, conflicts, and risk of bias.
6. Pair positive records with relevant null/harm studies and current guideline position.
7. Require independent scientific and clinical review before activation.
8. Compile only schema-valid assertions with immutable digests.
9. Recheck quarterly and immediately after guideline changes, regulator alerts, corrections, retractions, or formulation changes.
10. Expired or unreviewed assertions remain visible in the audit log but cannot render action language.

An LLM may assist candidate discovery or extraction review, but cannot directly generate a production treatment sentence. Production prose remains deterministic from approved claim records.

---

## 8. Safety, prescription, and FMT firewall

### 8.1 Safety-rule DSL

Use a restricted declarative language; never arbitrary Python or `eval`.

```yaml
rule_id: PROBIOTIC_IMMUNOCOMPROMISE_V1
applies_to:
  intervention_classes: [probiotic, live_biotherapeutic]
severity: clinician_only
when:
  any:
    - {field: context.immunocompromise.status, op: eq, value: true}
    - {field: context.central_venous_catheter, op: eq, value: true}
    - {field: context.critical_illness, op: eq, value: true}
    - {field: context.recent_hospitalization, op: eq, value: true}
on_unknown:
  fields: [context.immunocompromise.status, context.central_venous_catheter, context.critical_illness, context.recent_hospitalization]
  action: clinician_only
action:
  suppress_self_directed_option: true
  display_safety_reason: true
source_refs: [...]
```

Allowed operators are limited to `eq`, `ne`, `lt`, `lte`, `gt`, `gte`, `in`, `not_in`, `exists`, `all`, `any`, `none`, `medication_class_present`, `condition_present`, and `allergy_present`. Every field has a type, unit, and ontology. Every rule specifies unknown-field behavior.

Keep evidence, known medical-safety information, and rendering policy separate:

```yaml
known_safety_signal: none_known | caution | clinician_review | labeled_contraindication | unknown_incomplete_context
display_policy: participant_summary | clinician_only | evidence_only | suppressed
action_status: discussable_after_confirmation | evidence_only | suppressed
```

Deterministic precedence is: labeled contraindication > missing-safety-data suppression > clinician review > caution > no known signal. No positive study overrides a higher safety state. Participant output never says “safe” or “contraindicated” solely from incomplete/self-reported metadata; it states the known signal and why clinician review or suppression applies.

The compiler rejects every undeclared or misspelled field path, unresolved alias, invalid unit, and nullable referenced field lacking `on_unknown`. Runtime unknown values execute the declared fail-closed action. A preterm-infant rule always suppresses live-microbe participant content and presents the FDA invasive/fatal-infection warning and the fact that no probiotic is FDA-approved as a drug or biologic for infants: https://www.fda.gov/safety/medical-product-safety-information/risk-invasive-disease-preterm-infants-given-probiotics-formulated-contain-live-bacteria-or-yeast

### 8.2 Prescription medications

Prescription cards are clinician-only and require an exact confirmed indication. A microbiome pattern never establishes that indication. Participant mode never displays a prescription rank or dose. Each clinician/internal card shows the guideline/label source, diagnostic prerequisite, monitored risks, and the sentence:

> This is a clinician-review evidence route, not a recommendation generated from the microbiome result.

Do not reproduce a full disease-treatment guideline inside the microbiome report. Link to the current authoritative guideline and show only why the source is relevant.

### 8.3 Natural antimicrobials

Do not characterize berberine, allicin, oregano oil, black-seed oil, or multi-ingredient botanical products as proven eradication therapies for stool *E. coli* or *K. pneumoniae*. The often-cited 2014 herbal SIBO study was nonrandomized, choice-based, used a breath-test endpoint, and does not establish stool-pathogen eradication: https://pubmed.ncbi.nlm.nih.gov/24891990/

Botanical evidence based only on in-vitro susceptibility or animal studies remains lane E. Multi-ingredient products retain exact-product identity. Required safety rules include anticoagulant/antiplatelet interactions, glucose-lowering therapy, pregnancy, liver/kidney disease, allergy, surgery, GI bleeding, and product-quality uncertainty. No antimicrobial sequence is suggested from relative-abundance data.

### 8.4 FMT treatment evidence

The March 2024 AGA recommendations are conditional, low-certainty, indication- and immune-status-specific. They support fecal microbiota-based therapy primarily for prevention of recurrent *C. difficile* infection after standard-of-care antibiotics in selected immunocompetent or mildly/moderately immunocompromised adults, recommend against it in severely immunocompromised adults, and recommend against conventional FMT for Crohn's disease, ulcerative colitis, IBS, and pouchitis outside clinical trials: https://pubmed.ncbi.nlm.nih.gov/38395525/

FDA-approved microbiota products and conventional FMT are not interchangeable. REBYOTA and VOWST are approved only for adults age 18 or older to prevent recurrent CDI after antibacterial treatment; neither treats active CDI. Store exact infectious-agent/allergen warnings and current labeling, but never display participant doses or use a donor/profile result to select either product. Registry records must name the exact product/procedure, indication, route, and current regulatory source: https://www.fda.gov/vaccines-blood-biologics/fecal-microbiota-products

FDA's November 2022 enforcement-discretion policy is narrower: it concerns non-stool-bank FMT for *C. difficile* infection not responding to standard therapy, with treating-licensed-clinician donor/stool qualification, informed consent describing investigational status and risks, and no identified safety concern. It does not authorize FMT for other diseases: https://www.fda.gov/media/86440/download

Non-CDI FMT may be investigated only in a specific clinical investigation conducted under an effective IND, IRB oversight, and informed consent, unless FDA has determined that an IND is not required, or under a future FDA-approved indication/enforcement policy. An IND is not a treatment authorization. A favorable paper, professional guideline, or trial registration does not authorize treatment.

FDA safety communications document transmission of multidrug-resistant organisms, pathogenic *E. coli*, and other infectious risks, including a death after ESBL-producing *E. coli* transmission. Encode and date the full alert set, required informed-consent text, donor testing, donation-window/bookend testing, and quarantine implications:

- MDRO/ESBL Enterobacteriaceae, VRE, CRE, and MRSA — June 13/18, 2019: https://www.fda.gov/safety/medical-product-safety-information/fecal-microbiota-transplantation-safety-communication-risk-serious-adverse-reactions-due
- EPEC/STEC — March 12/13 and April 6, 2020: https://www.fda.gov/vaccines-blood-biologics/safety-availability-biologics/information-pertaining-additional-safety-protections-regarding-use-fecal-microbiota-transplantation
- SARS-CoV-2 — March 23 and April 9, 2020: https://www.fda.gov/vaccines-blood-biologics/fecal-microbiota-products
- Mpox — August 22, 2022: https://www.fda.gov/vaccines-blood-biologics/safety-availability-biologics/safety-alert-regarding-use-fecal-microbiota-transplantation-and-additional-safety-protections-0

For Alzheimer disease, alopecia, Long COVID, ME/CFS, Parkinson disease, metabolic disease, cancer, IBD, and IBS, FMT remains investigational and may be studied only in a specific clinical investigation under an effective IND, IRB oversight, and informed consent unless FDA has determined that an IND is not required, or unless a future FDA action changes the status. Two incidental alopecia-universalis regrowth cases after FMT are hypothesis-generating only: https://pubmed.ncbi.nlm.nih.gov/28932754/

### 8.5 Research donor prescreen

The only permitted output name is `research_donor_prescreen`, and it exists only in access-controlled `researcher` or authenticated `credentialed_clinician` mode. It renders per-lane findings and unassessed gaps only. It must never create an aggregate score, rank, green/pass state, “no exclusions” conclusion, or the words `eligible`, `safe`, `healthy donor`, `approved`, `cleared`, or `candidate`. Diversity and disease-pattern shape may appear only as non-weighted experimental observations; they are not donor-qualification criteria.

It must explicitly list clinical history, physical assessment, blood tests, validated stool pathogen tests, MDRO testing, viral/parasitic screening, repeat-window requirements, jurisdictional rules, and any unavailable lanes. A negative shotgun screen does not clear a donor. All subjects under 18 receive:

```yaml
research_donor_prescreen:
  status: not_supported_under_product_policy_minor
  policy_basis: categorical_product_policy_not_claimed_as_universal_federal_prohibition
```

Participant reports never expose those per-lane negative findings because they could create false reassurance or facilitate unsupervised FMT. They show only: **“Clinical donor qualification was not performed.”** No participant API payload may include the access-controlled prescreen object.

### 8.6 Regulatory release gate

Because the software analyzes NGS-derived patterns for disease-related purposes and displays treatment information to patients, the team must obtain a documented function-by-function classification and premarket-path determination before public release. This gate covers patient-facing disease concordance, intervention retrieval/selection, donor prescreening, and any health interpretation of microbiome age. A Q-Submission is a feedback process, not authorization, and a “research use only” disclaimer does not by itself remove device or advertising obligations. Review at minimum:

- FDA Clinical Decision Support Software guidance, January 2026: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/clinical-decision-support-software
- FDA General Wellness guidance, January 2026: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/general-wellness-policy-low-risk-devices
- FTC Health Products Compliance Guidance: https://www.ftc.gov/business-guidance/resources/health-products-compliance-guidance

Release modes such as `researcher`, `clinician`, and `participant` control presentation but do not override legal classification. Claim language, marketing pages, charts, and implied net impression all require review.

If any function is regulated as a device, implement QMSR design/development controls, risk management, analytical/clinical/software validation, change control, cybersecurity, complaint handling, applicable medical-device reporting, and corrections/removals. QMSR became effective February 2, 2026: https://www.fda.gov/medical-devices/postmarket-requirements-devices/quality-management-system-regulation-qmsr

Individual health-assessment reporting must use an appropriate CLIA-certified laboratory unless counsel documents a valid nonclinical research workflow; map state direct-access ordering/reporting rules. Do not assume a blanket LDT exemption or automatic FDA coverage: the 2024 LDT rule was vacated March 31, 2025, and prior regulatory text was restored September 19, 2025: https://www.fda.gov/medical-devices/in-vitro-diagnostics/laboratory-developed-tests

Operationalize the FTC Health Breach Notification Rule rather than treating it as a citation. Maintain a documented applicability determination for the product and each data flow; instrument detection and triage of unauthorized acquisition or disclosure; retain incident facts and affected-population counts; and assign controller/service-provider duties. When the rule applies, affected individuals must be notified without unreasonable delay and no later than 60 calendar days after discovery; the FTC must be notified within 10 business days for a breach involving 500 or more people; media notice is required where applicable; smaller breaches require the prescribed annual report; and service providers must notify the covered vendor or PHR-related entity. Version these deadlines and duties against the current rule: https://www.ftc.gov/business-guidance/resources/complying-ftcs-health-breach-notification-rule

### 8.7 Minor governance

Participant mode abstains completely from numeric adult-derived disease percentiles for a minor when population transport is unvalidated; such values may exist only in access-controlled internal research output. Under-18 users require guardian-controlled accounts, pediatric-specific human evidence, and clinician review before any probiotic/supplement content. Implement age assurance, conflict detection between age and minor status, verified parental consent when COPPA applies to an under-13 user, state-law consent mapping, assent/parental permission when covered research requires it, deletion/retention controls, and separate consent for secondary genomic use. Unknown or conflicting age fails closed.

---

## 9. Microbiome-predicted chronological age

### 9.1 Claim boundary

Call the output **estimated stool-microbiome chronological age (research use)**. Do not call it biological age, gut age, health age, aging rate, rejuvenation, or longevity score. An older- or younger-than-chronological prediction is not inherently good or bad.

Primary source: Myers et al., *Communications Biology* (2025), https://www.nature.com/articles/s42003-025-08590-y

Code, supplement, and data routes:

- Supplement PDF: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-025-08590-y/MediaObjects/42003_2025_8590_MOESM1_ESM.pdf
- Supplement workbook: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-025-08590-y/MediaObjects/42003_2025_8590_MOESM2_ESM.xlsx
- Analysis repository: https://github.com/tydymy/TRPCA_analysis
- Model repository: https://github.com/tydymy/TRPCA
- Frozen code DOI: https://doi.org/10.5281/zenodo.15801638
- WGS retrieval recipe: https://waldronlab.io/curatedMetagenomicDataAnalyses/articles/MLdatasets.html
- Earlier comparison data/code: https://github.com/shihuang047/age-prediction

Benchmark, but do not automatically import, earlier clocks on the same leakage-free held-out folds: [Galkin 2020](https://pubmed.ncbi.nlm.nih.gov/32534441/), [Chen 2022](https://pubmed.ncbi.nlm.nih.gov/35040752/), and [gAge 2024](https://pubmed.ncbi.nlm.nih.gov/38289284/). Claims that an age residual represents biological aging, frailty, or health require separate prospective validation; chronological prediction alone is insufficient.

The study aggregates 8,959 16S samples from 10 studies and 9,356 WGS samples from 56 studies; the stool WGS pool contains 8,641 samples spanning ages 0–91. Reported stool-WGS performance is about 8.83 years MAE for all samples and 9.72 years for a one-sample-per-subject analysis, versus roughly 9.98 and 10.86 years for reported KNN baselines. Global error is not pediatric accuracy.

### 9.2 Critical reproducibility finding

The published repositories do not supply a complete inference bundle: there is no verified production checkpoint with fitted RCLR/PCA/scalers, frozen feature order/schema, and exact training manifest suitable for direct new-FASTQ inference. Code audit also indicates preprocessing was performed globally before cross-validation in key notebooks, including fitting transformations and target scaling, which creates leakage risk. Methods, notebooks, and supplementary tables contain additional inconsistencies in filtering, pseudo-count construction, fold count, subject totals, and a country-accuracy value.

Therefore the report page is approved only after clean retraining and independent validation. Do not copy the published CV score into a production “age test.”

### 9.3 Reproduction and validation gate

1. Freeze one input namespace. Either reproduce cMD3 MetaPhlAn 3/CHOCOPhlAn-201901 or reprocess every training and inference FASTQ with one pinned current profiler/database; never mix MetaPhlAn 3 training with MetaPhlAn 4 inference.
2. Freeze exact sample/run/participant/project identifiers from the supplementary manifest, one baseline per participant, and remove duplicates/overlapping cohorts.
3. Fit prevalence filtering, RCLR/RPCA, PCA, all scalers, feature selection, and hyperparameters inside each training fold.
4. Use leave-one-study/project-out as the primary split; also participant-grouped, country-held-out, and prespecified pediatric holdout evaluations.
5. Benchmark KNN, regularized linear, random forest, SVR, and TRPCA on identical folds. Ship the simplest model whose external performance is not inferior.
6. Report MAE, median absolute error, RMSE, R², calibration slope/intercept, 80%/95% interval coverage, and error by age band, sex, BMI, region, study, and sequencing depth.
7. Validate on external raw FASTQs that were never used in model selection.
8. Freeze checkpoint, RCLR/PCA/scalers, feature order/IDs, taxonomy database, container, training manifest, checksums, model card, and inference CLI.
9. Verify all code/data licenses or use a documented clean-room reimplementation.
10. Predeclare adult-only/teen-specific models or fit age-bin weighting/sampling inside each training fold; the infant-heavy pooled dataset must not dominate model selection. Freeze the training age histogram/support density, select using age-stratified—not global-only—metrics, and prohibit extrapolation outside the validated range.
11. Use a frozen study-aware split-conformal interval method: fit on training folds, reserve an untouched calibration manifest, apply the finite-sample residual quantile (optionally prespecified normalized/heteroscedastic conformal), and validate 80%/95% coverage within every supported age stratum.
12. Freeze three OOD components—analytical/QC, latent-feature distance, and metadata/training-support density—with thresholds learned without the test specimen. Hard OOD must abstain; it cannot render an unqualified experimental number.

The supplementary workbook contains only about 332 samples ages 8–18 in the all-sample stool WGS sheet and about 143 in the one-sample-per-subject sheet. That is inadequate to infer reliable performance for the current young donors. Numeric pediatric output requires at least 200 independent 8–18-year-old participants from at least three external cohorts, laboratories, and regions with no training overlap and preregistered age-band MAE, calibration-slope/intercept, and 80%/95% coverage bounds; otherwise show `insufficient_pediatric_reference_support`.

### 9.4 Model bundle

```yaml
age_model_bundle:
  model_id: MPA_WGS_AGE_V1
  task: predicted_chronological_age
  input_namespace: ...
  database_version: ...
  preprocessing_digest: sha256:...
  feature_digest: sha256:...
  estimator_digest: sha256:...
  training_manifest_digest: sha256:...
  eligible_age_range: ...
  eligible_countries: [...]
  validation:
    participant_grouped: true
    study_held_out: true
    country_held_out: true
    pediatric_external: ...
    mae_years: ...
    calibration_slope: ...
    interval_coverage_95: ...
    conformal_method: ...
    calibration_manifest_digest: sha256:...
    validation_acceptance_policy_ref: ...
  training_age_support:
    histogram_digest: sha256:...
    support_density_digest: sha256:...
  ood_policy:
    analytical_threshold: ...
    latent_distance_threshold: ...
    metadata_support_threshold: ...
  license_review: pass
```

### 9.5 Required age output

```yaml
microbiome_age:
  status: scored | abstained_ood | not_computable
  reason: ... | null
  predicted_chronological_age_years: ... | null
  prediction_interval_95: [...] | null
  actual_chronological_age_years: ... | null
  age_residual_years: ... | null
  matched_stratum_mae_years: ...
  calibration_population: ...
  ood_status: in_distribution | warning | out_of_distribution
  model_bundle_id: ...
  contributions_increasing_this_prediction: [...]
  contributions_decreasing_this_prediction: [...]
  attribution_method_and_background_digest: sha256:...
```

Chronological age is not a predictor feature. It is required for residual calculation, age-stratified applicability, and minor/adult release control; missing actual age never defaults to adult. An estimate may be generated in access-controlled internal research mode without actual age, but participant mode with unknown age assurance returns `status: not_computable` and `reason: unknown_age_assurance`. Display the residual only when actual age is supplied. Participant numeric estimates require every OOD gate to pass and sufficient support density around the predicted age; hard OOD returns `abstained_ood` with null estimate/interval. Attribution values are local contributions to this prediction, retain only features stable across folds/bootstrap, and are not general age biology, causality, or treatment targets. No intervention is generated from age or residual alone.

### 9.6 Age-page visual

The age page contains:

1. A neutral chronological axis with actual age, predicted-age point, and 95% prediction interval.
2. A status line mapped exactly to the output contract: `supported` for `status: scored` plus `ood_status: in_distribution`; `warning—in distribution` for `status: scored` plus `ood_status: warning`; `abstained—out of distribution` for `status: abstained_ood`; or `not computable` for `status: not_computable`. Hard-OOD participant output has null estimate and interval. Any qualified numeric OOD research result lives only in a separately named access-controlled internal object and never on this page.
3. Matched-stratum error and sample size beside the estimate—not hidden in a footnote.
4. Calibration and held-out-study performance mini-table.
5. The nearest reference age/country/study regions and the sample's OOD distance.
6. Stable local contributions increasing/decreasing this prediction shown in blue/gray, with attribution method/background, never labeled universal older/younger organisms or good/bad.
7. A permanent statement: “This estimates chronological age from stool community patterns. It is not biological age, a health grade, or a treatment target.”

If the prediction interval spans more than one life-stage band, the headline must visually emphasize the interval over the point estimate.

---

## 10. Report information architecture

### 10.1 Front results atlas

The report is not limited to one summary page. All clear visual readings come first; detailed prose comes afterward.

1. **Executive snapshot:** overall state, reference applicability, top favorable/concerning findings, counts, and disease-profile shape verdict.
2. **What stood out:** concise prose, urgent safety/clinical-follow-up flags, and strongest limitations.
3. **Data validity and reference fit:** QC, age/sex/country, medications/antibiotics, stool form, matched-stratum `n`, fallback, OOD, and upstream host-removal status.
4. **Complete global metrics:** GMWI2, richness, Shannon, evenness, neutral F:B, classified fractions, and microbiome-age output.
5. **All 25 functional measurements:** complete visual matrix, no selective subset.
6. **All 10 microbial groups:** complete matrix plus cross-engine contradictions.
7. **Taxa needing attention:** relatively high, relatively low, qualified expected-but-not-detected, uncommon present, potential-concern qualified, and indeterminate.
8. **Complete detected-species inventory:** as many pages as required; scoring and extended lanes distinct.
9. **All disease profiles:** every score/status plus shape, intervals, coverage, maturity, specificity, OOD, and residualized/sensitivity results.
10. **Research-evidence overview:** deduplicated trigger/evidence mappings, evidence lanes, conflicts, safety gates, and clinician-only routes; no participant benefit ranking.
11. **Experimental FMT research prescreen:** access-controlled researcher/credentialed-clinician mode only and separated from clinical eligibility. Participant mode shows only “Clinical donor qualification was not performed.”
12. **Functional detail pages:** measurement, meaning, contradiction checks, and evidence cards.
13. **Taxon/group detail pages:** reference distribution, detection qualification, context, and evidence cards.
14. **Disease-profile detail pages:** every notable or user-requested profile, including AA, AGA, AD continuum, ME/CFS, and Long COVID.
15. **Limitations, methods, QC, extended catalog, validation, technical provenance, and references.**

Every section begins with a full-width titled divider and a plain-language reading guide.

Detail selection is deterministic. All 33 registered profile tasks appear in the front table. Every `scored` or `scored_low_coverage` profile receives a full detail record/page; every abstained/declined/not-computable profile receives a compact reason card with missing prerequisites. An explicitly requested profile always receives detail. The renderer asserts `detail_count == count(render_detail=true)` and has no maximum cap. Taxon JSON remains namespace-complete under Section 6.5, while PDF inclusion follows its explicit detected/expected/uncommon/concern/requested rules.

### 10.2 Universal infographic grammar

| Metric type | Required visual |
|---|---|
| Bounded continuous | 5th–95th reference strip, IQR, median, sample marker, exact value, percentile, interval |
| Unbounded continuous | Reference-density or quantile strip, never an unlabeled dial |
| Two-sided optimum | Stable axis and evidence-backed target corridor; low/high adverse sides explicit |
| One-sided direction | Separate `level` from `interpretation` |
| Presence/bimodal | Detected / qualified not detected / indeterminate tile with prevalence and detection power |
| Descriptive ratio | Neutral distribution; no red/green judgment without validation |
| Disease pattern | Concordance with interval, coverage, maturity, specificity, OOD, and shape context |
| Expected nondetection | Reference prevalence, expected abundance, detection power, and reason |
| Uncommon detection | Abundance percentile, prevalence, confirmation/strain qualification |
| Composition | Neutral stacked bar/donut; no intrinsic healthy phylum claim |
| QC | Threshold strip with pass/warn/fail/not-assessable and exact rule |
| Model prediction | Estimate, prediction interval, matched-stratum error, OOD, model/version |

Use green only for an evidence-supported favorable interpretation, red only for an evidence-supported unfavorable interpretation, amber for uncertainty/conflict/transport, and blue/gray for descriptive/neutral/not assessable. Never rely on color; every status also has text and an icon/shape.

### 10.3 Metric-card fields

Every metric card shows:

- title and typed claim level;
- exact value and unit;
- reference band and sample marker;
- reference cohort, `n`, age/sex/country stratum, database, and version;
- percentile and uncertainty;
- level (`low/high/etc.`) separate from interpretation (`favorable/concerning/neutral/unknown`);
- OOD/reference fallback;
- longitudinal delta when compatible;
- direct link or reference marker.

Shannon and evenness do not receive universal “good” labels. Firmicutes:Bacteroidetes is descriptive and visually de-emphasized. Gene-capacity cards must say “capacity,” not production or concentration.

### 10.4 PDF and accessibility contract

- Main body at least 9.5 pt; tables at least 8.5 pt; footnotes at least 7.5 pt; never shrink below minima.
- WCAG-AA contrast for digital output; grayscale and deuteranopia QA.
- No dynamic-content truncation or unexplained ellipsis.
- No orphan continuation page below 30% content area unless intentional.
- Embed all fonts; set document language; generate tagged reading order, bookmarks, repeated table headers, alt text, and clickable DOI/PubMed/guideline links.
- Header/footer and reference markers remain legible in print.
- All result counts in the PDF equal `results.json`; no hidden renderer cap.
- Run automated `pdffonts`, `pdfinfo`, link-annotation, text-overflow, and page-density tests plus representative raster/contact-sheet review.

---

## 11. Longitudinal analysis

### 11.1 Compatibility before change

Two specimens may be compared quantitatively only when body site, collection method, preservative, extraction, assay, databases, and pipeline are compatible, or when raw reads are reprocessed together through one frozen bundle. Subject identity and collection/exposure dates must be known. Otherwise emit `longitudinal_status: incompatible_measurement`.

### 11.2 Change model

For transformed quantitative feature `f`:

\[
\Delta_f=T_f(x_{t2})-T_f(x_{t1})
\]

Combine read-sampling uncertainty, technical variation, and reference-study bootstrap. Call a reliable change only if its interval excludes zero, it exceeds a frozen minimal detectable change, and direction is stable under sensitivity analyses. Disease-profile deltas use a paired bootstrap over independent evidence clusters.

Permitted language: “Relative abundance decreased after the reported exposure.”  
Prohibited language: “The supplement eradicated the organism,” “caused recovery,” or “treated the disease.”

Repeat specimens are especially important because day-to-day abundance can vary substantially. Display symptoms and clinical labs on a parallel timeline, never merged into the microbial score.

---

## 12. Top-level results schema

```json
{
  "schema_version": "4.0.0",
  "provenance": {},
  "sample_manifest": {},
  "qc": {},
  "observations": [],
  "reference_calls": [],
  "findings": [],
  "functional_panels": [],
  "microbial_groups": [],
  "taxonomic_inventories": {},
  "disease_profiles": [],
  "disease_shape": {},
  "microbiome_age": {},
  "research_evidence_summaries": [],
  "safety_evaluation": {},
  "longitudinal_results": [],
  "access_controlled_outputs": {
    "research_donor_prescreen": {}
  },
  "subsystem_status": {
    "intervention_retrieval": "complete|unavailable|failed",
    "age": "scored|abstained_ood|not_computable|not_deployed",
    "donor_prescreen": "complete|not_supported_product_policy|mode_redacted|unavailable"
  },
  "source_ledger": [],
  "legacy_v3": {}
}
```

`access_controlled_outputs` is omitted—not merely emptied—from participant payloads. The participant renderer emits the fixed non-qualification statement and cannot retrieve or serialize per-lane donor findings. The intervention registry is optional at runtime. If it is unavailable, all measurement, reference, and disease-profile stages still complete and the report declares that evidence retrieval was unavailable.

Examples in this document do not replace normative schemas. Ship JSON Schema 2020-12 contracts for `Observation`, `ReferenceDefinition`, `ReferenceCall`, `Finding`, `InterpretationAssertion`, `TriggerRecord`, intervention identities/protocols/evidence, `LongitudinalResult`, age output, safety evaluation, and the complete donor-prescreen record. Every schema has `$id`, required/nullability rules, enums, units, `additionalProperties: false`, foreign-key checks, and registry-lock references. Every model/reference bundle carries a preregistered `validation_acceptance` block; a failed gate forces experimental/not-deployable status, never post-hoc threshold adjustment.

---

## 13. Disease-profile intervention actionability registry

All 30 live profile tasks remain disease-associated microbiome **pattern concordance** under v03. V04 adds two missing Alzheimer stages and restores the v03-promised PDAC task for an exact target of 33. The actionability-family rows below cover all of them; they seed evidence assertions and do not redefine disease signatures. Every medication or disease-management route requires independent clinical confirmation. Every non-CDI FMT status here is `specific_effective_IND_investigation_or_FDA_determined_IND_not_required_only`; this records a research-route constraint, not treatment authorization.

### 13.1 Gastrointestinal, metabolic, renal, and hepatic profiles

| Profile | Evidence cards permitted | Required suppression/gate | Primary sources |
|---|---|---|---|
| Colorectal cancer | Population prevention evidence for fiber/whole grains, weight, alcohol, and red/processed meat; established screening routes | No microbiome-directed prevention/treatment; no “cancer detected”; no FIT/colonoscopy substitution; do not treat *F. nucleatum* from a profile | [USPSTF screening](https://www.uspreventiveservicestaskforce.org/uspstf/recommendation/colorectal-cancer-screening); [WCRF report](https://www.wcrf.org/wp-content/uploads/2024/10/Colorectal-cancer-report.pdf) |
| Pancreatic ductal adenocarcinoma (restored v04 task) | Two independent 2022 multinational stool-shotgun classifier studies as disease-signature evidence only; conventional evaluation route for independently concerning symptoms/risk | No “pancreatic cancer detected,” population screening, CA19-9 imputation, supplement/probiotic/antimicrobial/FMT treatment, or oncology decision from stool | [Kartal et al.](https://pubmed.ncbi.nlm.nih.gov/35260444/); [Nagata et al.](https://pubmed.ncbi.nlm.nih.gov/35398347/); [NCI clinical route](https://www.cancer.gov/types/pancreatic/patient/pancreatic-treatment-pdq) |
| Colorectal adenoma | Exact CBM588 human trial as early/conflicting evidence; meta-analysis context | No adenoma diagnosis or surveillance change; no generic *Clostridium butyricum* substitution | [CBM588 RCT](https://pubmed.ncbi.nlm.nih.gov/41425719/); [2024 meta-analysis](https://pubmed.ncbi.nlm.nih.gov/38126945/) |
| Crohn disease | Guideline position: probiotics and conventional FMT only in trials; clinician-confirmed disease may link to current guideline care | No diagnosis, activity/severity claim, steroid/biologic/JAK/surgery decision from stool | [AGA probiotics](https://pubmed.ncbi.nlm.nih.gov/32531291/); [AGA FMT](https://pubmed.ncbi.nlm.nih.gov/38395525/); [CDED+PEN RCT](https://pubmed.ncbi.nlm.nih.gov/31170412/) |
| Ulcerative colitis | Positive FMT trial signal displayed beside current trial-only guideline; probiotic evidence remains formulation-specific/trial context | No diagnosis, activity, mesalamine/steroid/biologic/JAK/colectomy decision | [Paramsothy FMT RCT](https://pubmed.ncbi.nlm.nih.gov/28214091/); [AGA FMT](https://pubmed.ncbi.nlm.nih.gov/38395525/); [AGA probiotics](https://pubmed.ncbi.nlm.nih.gov/32531291/) |
| Composite IBD | No targeted intervention for a generic composite; route to symptoms, fecal calprotectin, and gastroenterology evaluation if clinically relevant | Never collapse Crohn and UC into one treatment; do not infer mucosal inflammation | [AGA FMT](https://pubmed.ncbi.nlm.nih.gov/38395525/) |
| Type 2 diabetes | High-fiber dietary RCT affecting HbA1c/SCFA-producer ecology; lean-donor FMT and exact *Akkermansia* products as research-stage evidence | Metformin status remains mandatory; no diabetes diagnosis or medication change; endogenous low *Akkermansia* does not authorize an OTC product | [high-fiber RCT](https://pubmed.ncbi.nlm.nih.gov/29590046/); [FMT 2012](https://pubmed.ncbi.nlm.nih.gov/22728514/); [FMT 2017](https://pubmed.ncbi.nlm.nih.gov/28978426/); [AKK-WST01 RCT](https://pubmed.ncbi.nlm.nih.gov/39879980/) |
| Cirrhosis / hepatic encephalopathy | Clinician-confirmed HE guideline route for lactulose/rifaximin; probiotic and FMT studies as adjunct/research evidence | No cirrhosis/HE diagnosis, drug start, transplant decision, or hepatotoxic botanical card from profile | [EASL 2022](https://pubmed.ncbi.nlm.nih.gov/35724930/); [probiotic RCT](https://pubmed.ncbi.nlm.nih.gov/24246768/); [meta-analysis](https://pubmed.ncbi.nlm.nih.gov/26561214/); [FMT pilot](https://pubmed.ncbi.nlm.nih.gov/28586116/) |
| Obesity | Exact pasteurized *A. muciniphila* product evidence and null-primary adolescent FMT evidence | No obesity diagnosis, GLP-1/bariatric decision, or inference that endogenous abundance predicts response | [pasteurized *Akkermansia* RCT](https://doi.org/10.1038/s41591-026-04394-7); [adolescent FMT RCT](https://pubmed.ncbi.nlm.nih.gov/33346848/) |
| MASLD | Clinician pathway based on metabolic risk/fibrosis; null small FMT trial | No MASLD/MASH diagnosis, fibrosis stage, drug/biopsy/supplement decision from stool | [EASL-EASD-EASO 2024](https://pubmed.ncbi.nlm.nih.gov/38851997/); [AASLD guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC10735173/); [FMT RCT](https://pubmed.ncbi.nlm.nih.gov/32618656/) |
| Type 1 diabetes | Small new-onset adult FMT trial as hypothesis-generating evidence | No diagnosis, insulin initiation/reduction, immune therapy, or ketoacidosis reassurance | [2021 FMT RCT](https://pubmed.ncbi.nlm.nih.gov/33106354/) |
| Celiac disease | Gluten-free diet only after proper clinical diagnosis; route to serology/biopsy guidance | Do not initiate gluten avoidance from a microbiome pattern because it may impair diagnostic testing | [ACG 2023 guideline](https://pubmed.ncbi.nlm.nih.gov/36602836/) |
| CKD | Small FMT RCT and resistant-starch/probiotic/prebiotic/synbiotic surrogate evidence; conventional CKD guideline route | No CKD diagnosis/stage/dialysis/drug decision; diet options require eGFR, potassium, phosphorus, and renal-diet review | [FMT RCT](https://pubmed.ncbi.nlm.nih.gov/38674803/); [resistant-starch review](https://pubmed.ncbi.nlm.nih.gov/39444299/); [biotic meta-analysis](https://pubmed.ncbi.nlm.nih.gov/34640474/); [KDIGO 2024](https://pubmed.ncbi.nlm.nih.gov/38490803/) |
| IBS composite | Dietitian-guided limited low-FODMAP trial after clinical evaluation, including restriction, reintroduction, and personalization; soluble fiber/peppermint evidence in guideline | No diagnosis without criteria/red-flag review; ACG suggests against probiotics for global IBS; no taxon-as-infection language | [ACG 2021](https://pubmed.ncbi.nlm.nih.gov/33315591/); [AGA diet update](https://www.gastrojournal.org/article/S0016-5085%2821%2904084-1/fulltext) |
| IBS-C | Same diet/fiber route; clinician-confirmed constipation therapies as guideline links | No secretagogue/laxative/motility decision from the profile | [ACG 2021](https://pubmed.ncbi.nlm.nih.gov/33315591/) |
| IBS-D | Clinician-confirmed guideline route may mention rifaximin for global symptoms | Rifaximin is not treatment for organisms seen in shotgun data; no infection/SIBO diagnosis | [ACG 2021](https://pubmed.ncbi.nlm.nih.gov/33315591/) |
| IBS-M | Dietitian-guided diet/soluble-fiber evidence after clinical diagnosis | No antimicrobial or secretagogue selection from fluctuating microbial results | [AGA diet update](https://www.gastrojournal.org/article/S0016-5085%2821%2904084-1/fulltext) |

### 13.2 Neurologic, immune, cardiovascular, psychiatric, hair, and post-infectious profiles

| Profile | Evidence cards permitted | Required suppression/gate | Primary sources |
|---|---|---|---|
| Parkinson disease | Exact multistrain probiotic for constipation endpoint only; conflicting FMT trials displayed together | No diagnosis, prognosis, motor/neuroprotection claim, or dopaminergic medication change | [constipation RCT](https://pubmed.ncbi.nlm.nih.gov/33046607/); [FMT RCT, no meaningful benefit](https://pubmed.ncbi.nlm.nih.gov/39073834/); [FMT RCT, motor signal](https://pubmed.ncbi.nlm.nih.gov/38686220/) |
| Multiple sclerosis | Small probiotic immune/microbiome studies and nine-person terminated FMT pilot as research-only | No diagnosis/activity/MRI/relapse or disease-modifying-therapy decision | [FMT pilot](https://pubmed.ncbi.nlm.nih.gov/35571974/); [2018 probiotic trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6181139/); [2024 crossover](https://pmc.ncbi.nlm.nih.gov/articles/PMC11470793/) |
| Atherosclerotic CVD | Mediterranean dietary-pattern event evidence; red-meat feeding effect on measured circulating TMAO | No ASCVD diagnosis/risk or statin/antiplatelet/anticoagulant decision; TMA gene capacity is not plasma TMAO | [PREDIMED reanalysis](https://pubmed.ncbi.nlm.nih.gov/29897866/); [red-meat/TMAO trial](https://pubmed.ncbi.nlm.nih.gov/30535398/) |
| Rheumatoid arthritis | Small exact-formulation probiotic trials, including positive and null results, presented as conflicting early evidence | No diagnosis/activity or steroid/DMARD/JAK/biologic decision | [2014 trial](https://pubmed.ncbi.nlm.nih.gov/24673738/); [2016 trial](https://pubmed.ncbi.nlm.nih.gov/27135916/); [negative trial](https://pubmed.ncbi.nlm.nih.gov/21629190/); [FMT case](https://pmc.ncbi.nlm.nih.gov/articles/PMC7869316/) |
| Ankylosing spondylitis | Negative probiotic RCT and case-level FMT literature | No diagnosis/imaging/activity or NSAID/biologic/JAK decision | [negative probiotic RCT](https://pubmed.ncbi.nlm.nih.gov/20716665/); [FMT case](https://pubmed.ncbi.nlm.nih.gov/36911747/) |
| Hypertension | DASH-style lifestyle evidence independent of stool; transient/null-duration FMT trial context | No diagnosis or medication titration; require repeated validated blood-pressure measurement | [PREMIER RCT](https://pubmed.ncbi.nlm.nih.gov/15385946/); [FMT RCT](https://pubmed.ncbi.nlm.nih.gov/40410854/) |
| Alzheimer-associated, preclinical | Multidomain lifestyle evidence independent of microbiome; no established microbiome treatment | Never label the subject preclinical AD or infer amyloid/tau/future dementia; no preventive drug | [FINGER RCT](https://pubmed.ncbi.nlm.nih.gov/25771249/) |
| Mild cognitive impairment continuum | Small/inconsistent probiotic evidence and null comparative MIND-diet trial | No MCI diagnosis/prognosis, biomarker, capacity/driving, or dementia medication decision | [2024 evidence review](https://pmc.ncbi.nlm.nih.gov/articles/PMC11863739/); [MIND trial](https://pubmed.ncbi.nlm.nih.gov/37466280/) |
| Clinical Alzheimer disease | Small probiotic-milk trial shown beside inconsistent pooled evidence; human FMT case literature only | No AD diagnosis/stage or cholinesterase inhibitor, memantine, or anti-amyloid route from stool | [probiotic RCT](https://pubmed.ncbi.nlm.nih.gov/27891089/); [FMT review](https://pubmed.ncbi.nlm.nih.gov/36381829/) |
| ME/CFS | BioMapAI as classifier evidence only; small probiotic pilot; randomized FMT pilot with no symptom/QoL benefit | No diagnosis, fatigue attribution, graded-exercise, pacing, or medication decision from a taxon/profile | [BioMapAI study](https://www.nature.com/articles/s41591-025-03788-3); [probiotic pilot](https://pmc.ncbi.nlm.nih.gov/articles/PMC2664325/); [FMT pilot RCT](https://pubmed.ncbi.nlm.nih.gov/37516837/) |
| Major depressive disorder | Exact high-dose multistrain add-on trial, null-remission FMT-plus-escitalopram trial, and diet adjunct evidence | Opt-in only; no diagnosis/antidepressant decision/suicide inference; crisis logic uses symptoms, never stool | [probiotic RCT](https://pubmed.ncbi.nlm.nih.gov/35654766/); [FMT RCT](https://pubmed.ncbi.nlm.nih.gov/42309058/); [SMILES RCT](https://pubmed.ncbi.nlm.nih.gov/28137247/) |
| Alopecia areata | Two incidental post-FMT regrowth cases and a tiny probiotic RCT with nonsignificant primary hair outcomes, visibly labeled hypothesis-generating/null-primary | No diagnosis/subtype/severity/JAK/steroid/immunotherapy/FMT recommendation | [FMT case report](https://pmc.ncbi.nlm.nih.gov/articles/PMC5599691/); [additional case](https://pubmed.ncbi.nlm.nih.gov/31624757/); [probiotic pilot RCT](https://www.mdpi.com/2079-9284/11/4/119) |
| Androgenetic alopecia | Exact three-strain *Lactiplantibacillus* 16-week trial, with significant and nonsignificant endpoints separated | Sex/age matching required; no diagnosis, minoxidil/finasteride decision, generic probiotic substitution, or pregnancy-safety inference | [2024 RCT](https://pubmed.ncbi.nlm.nih.gov/39275216/) |
| Long COVID / PASC | Exact SIM01 six-month RCT and smaller STOP-FATIGUE trial as formulation-specific evidence | No diagnosis, anticoagulant/antiviral/immunomodulator/antibiotic/supplement stack; classifier does not predict response | [SIM01 RCT](https://pubmed.ncbi.nlm.nih.gov/38071990/); [STOP-FATIGUE](https://pubmed.ncbi.nlm.nih.gov/39592468/) |

### 13.3 Alzheimer and alopecia are mandatory scored systems

The presence of weak or early intervention evidence does not remove these profiles. Canonical profile-to-module mappings prevent one cohort's multiple modalities from becoming duplicate profiles:

- Alzheimer profiles: `AD_SCS_METAAD_V1`, `AD_SCD_METAAD_V1`, `AD_MCI_METAAD_V1`, `AD_CLINICAL_METAAD_V1`, and `AD_PRECLINICAL_AMYLOID_V1`.
- `ALOPECIA_AREATA_V1` contains `AA_PEDIATRIC_STOOL_V1`, adult/human taxonomic cohort modules, the direction-conflict KO sensitivity modules, and optional measured-VOC module. These modules contribute once by study group; they are not standalone profiles.
- `ANDROGENETIC_ALOPECIA_V1` contains female- and male-stool modules. Sex-unknown two-stratum range is permitted only when the frozen fallback explicitly enables it and never becomes two contributions to profile shape.

Every module declares `parent_profile_id` and `study_group`. The compiler rejects module IDs in the expected-profile set and rejects duplicate contributions from the same participants.

Their cards must say “Similarity to published [condition]-associated microbiome signatures.” They must never say “[condition] probability,” “disease detected,” or “you have [condition].” Conventional confirmation routes are neurologic/cognitive evaluation and accepted biomarkers for the Alzheimer continuum, and dermatologist/scalp examination for alopecia.

### 13.4 New 2026 two-year Long-COVID mechanistic module

Register the newly published Zhang et al. source as a **candidate association/mechanism module**, not a production score or human-treatment assertion: https://doi.org/10.1007/s12602-026-11197-2

The reported human cohort was only 11 two-year Long-COVID participants and 11 SARS-CoV-2-unexposed controls from Hainan, China. It combined stool shotgun metagenomics with fecal LC-MS, then transferred donor material into antibiotic-treated mice. Reported human directions included higher *Streptococcus salivarius* and *S. parasanguinis*, lower `Faecalibacterium_SGB15346` and *Alistipes onderdonkii*, lower diversity, and impaired carbohydrate degradation, indole, SCFA, and fatty-acid-degradation signals. Mouse recipients showed lung/intestinal inflammation, anxiety-like outcomes, and worse *K. pneumoniae* infection outcomes; *S. salivarius* experiments were also performed.

Implementation rules:

- bind only taxonomic and DNA-observable features to FASTQ;
- keep LC-MS metabolites `not_measured` unless the exact external assay is supplied;
- keep mouse-transfer and isolate experiments in lane E/mechanistic;
- do not call *S. salivarius* universally pro-inflammatory—the paper itself cites strain/context-dependent prior work;
- do not generate antimicrobial, probiotic, or FMT treatment from these directions;
- require independent human cohort replication, participant-level data, and study-held-out validation before score contribution.

The article's public page says “No datasets were generated or analysed during the current study,” despite describing new shotgun/LC-MS data, and exposes no accession. Treat this as a reproducibility blocker. Until authors provide participant-level feature tables/raw accession and exact methods, compile the module as `unbound_source_candidate` with `would_bind_if: public_or_auditable_data_and_replication`.

### 13.5 Restore the missing pancreatic ductal adenocarcinoma task

Create `PDAC_STOOL_MULTINATIONAL_V1` from two independent 2022 stool-shotgun programs, not from a hand-written union of favorable taxa:

1. Kartal et al. used 57 treatment-naive PDAC cases, 50 controls, and 29 chronic-pancreatitis participants in Spanish discovery, plus a 76-person German validation set. Its 27-species fecal model reached reported AUROC up to 0.84, was tested across disease stages and 25 external disease datasets, and improved when combined with measured serum CA19-9. Source, supplements, conflict/patent disclosure, and public-data route: https://pubmed.ncbi.nlm.nih.gov/35260444/
2. Nagata et al. derived gut/oral signatures in treatment-naive Japanese participants and externally evaluated Spanish/German cohorts; it reported 30 gut species and model AUROCs about 0.78–0.82. Keep stool and oral models separate and use only stool-observable features here: https://pubmed.ncbi.nlm.nih.gov/35398347/

Implementation contract:

- retrieve exact supplementary feature coefficients, preprocessing, outcome labels, accession/run manifests, and licenses; never reconstruct the classifier from abstracts, figures, or a patent list;
- uniformly reprocess every accessible FASTQ in a single frozen namespace or reproduce the papers' exact frozen feature space;
- group overlapping Spanish/German participants across the papers before any split and hold out whole country/study groups;
- preserve chronic pancreatitis, other-cancer/GI-disease challenges, treatment status, age, country, diabetes, pancreatitis, antibiotics, PPI, jaundice/biliary obstruction, and sequencing protocol as covariates/challenge strata;
- encode each published model as a locked module and add a transparent direction panel only as sensitivity analysis; do not average models before external validation;
- measured CA19-9 is an optional external-lab module and remains `not_measured` for FASTQ-only input; never impute it;
- require participant-grouped, country-held-out AUROC/AUPRC, calibration, fixed-specificity sensitivity, decision-curve context, cross-disease false-positive matrix, and prospective asymptomatic/high-risk validation before any claim beyond uncalibrated pattern concordance;
- display `Similarity to published PDAC-associated stool microbiome signatures`, never “PDAC screen,” “early detection,” or cancer probability;
- no microbe-directed PDAC intervention exists for this report. Symptoms such as jaundice, unexplained weight loss, or persistent upper abdominal/back pain trigger ordinary urgent clinical evaluation independently of the score: https://www.cancer.gov/types/pancreatic/patient/pancreatic-treatment-pdq

---

## 14. Reference, training, and validation data registry

### 14.1 Primary healthy/control substrates

| Resource | Role | Strength | Bounded limitation | Direct access |
|---|---|---|---|---|
| curatedMetagenomicData v3 | Primary uniformly processed taxonomic/function discovery and reference substrate | 22,710 metagenomes, 94 cohorts, 42 countries, extensive metadata; MetaPhlAn 3 and HUMAnN 3 matrices | Moving resource, cohort imbalance, MetaPhlAn 3 namespace, some controlled raw data; pin exact release/sample IDs | [2026 paper](https://www.nature.com/articles/s41467-025-66888-1); [pipeline](https://waldronlab.io/curatedMetagenomicData/articles/our-pipeline.html); [CLI snapshot](https://zenodo.org/records/17498288); [metadata](https://zenodo.org/records/17498348) |
| HMP1 | Adult analytical/reference anchor | About 300 healthy US adults and >1,200 WGS specimens across body sites | Narrow adult US population; repeated body sites; not pediatric/global | [HMPDACC](https://www.hmpdacc.org/hmp/HMASM); [NIH HMP](https://commonfund.nih.gov/hmp) |
| PREDICT1 | Adult diet/metabolic context | 1,098 UK/US volunteers; 1,203 deep stool metagenomes | Adult, geographically limited; not a clinical reference interval | [paper/data description](https://pmc.ncbi.nlm.nih.gov/articles/PMC8353542/); `PRJEB39223` |
| Dutch Microbiome Project | Host/exposure matching research | 8,208 Dutch individuals, family structure, 241 host/exposure variables | Controlled access and geographically narrow; explicitly does not establish a universal healthy signature | [paper](https://www.nature.com/articles/s41586-022-04567-7); [EGA](https://ega-archive.org/studies/EGAS00001005027) |
| LifeLines-DEEP | Adult Dutch multi-omic/reference context | 1,135 specimens in the age-study pool | Controlled access; regional cohort | [EGA study](https://ega-archive.org/studies/EGAS00001001704); [dataset](https://ega-archive.org/datasets/EGAD00001001991) |
| GMrepo | Discovery and cohort-stratification catalog | Tens of thousands of samples and many phenotypes | Not a ready-made normative matrix; raw data require uniform reprocessing | [GMrepo paper](https://academic.oup.com/nar/article/50/D1/D777/6426060); [documentation](https://evolgeniusteam.github.io/gmrepodocumentation/) |
| GutFeelingKB | Historical candidate-organism catalog | Curated candidate list with public code/data | Small/self-reported cohort and older alignment; never a universal core checklist | [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6738582/); [repository](https://github.com/GW-HIVE/GFKB) |
| NIST RM 8048 | Batch, repeatability, and reproducibility QC | Standardized human-stool reference material | Analytical control, not a healthy-population interval | [NIST program](https://www.nist.gov/programs-projects/human-gut-microbiome-reference-material); [2026 characterization](https://doi.org/10.6028/NIST.SP.260-252) |
| UHGG | Taxonomy/genome namespace and mapping resource | Large human-gut genome catalog | Catalog, not prevalence or healthy range | [UHGG paper](https://www.nature.com/articles/s41587-020-0603-3) |

Every normative reference model must occupy one exact frozen processing/database space: either use one pinned curatedMetagenomicData processing snapshot and reproduce that exact inference pipeline, or reprocess every included raw FASTQ together. Incompatible, unavailable, or controlled-only cohorts without a compatible frozen profile are excluded or remain separate reference models; “when possible” mixing is prohibited. Use one baseline per participant, model study/lab/batch effects, and preserve controlled-access status. A portal sample count does not equal a validated reference cohort.

### 14.2 Age-training project map

The Myers et al. supplement can be joined to curatedMetagenomicData accessions, but the paper does not provide a frozen raw-run manifest. Create `refs/age/MYERS_WGS_STOOL_RECONSTRUCTION_V1.parquet` with study, participant, BioSample, run list, time point, accession, age, sex, country, eligibility, overlap cluster, access status/license, file checksum, and inclusion/exclusion reason. Assert the reported 8,641-row pre-dedup reconstruction and record the post-dedup count. The 52 parenthesized project-seed counts below do arithmetically sum to 8,641, but aggregate totals are not a row-level identity check. Regenerate the project map directly from a workbook-to-curatedMetagenomicData accession join, enumerate every missing/changed/duplicated row and final study label, and store `reported_pre_dedup_n: 8641`, `listed_seed_n: 8641`, workbook/release checksums, and join diagnostics. Training is blocked until the row-level manifest, count, deduplication assertions, and project totals reconcile. The following map is therefore a download seed only; biological-sample-to-run mappings and overlaps must be resolved before training:

```text
AsnicarF_2017 PRJNA339914 (8); AsnicarF_2021 PRJEB39223 (1098)
BackhedF_2015 PRJEB6456 (287); Bengtsson-PalmeJ_2015 PRJEB7369 (70)
BritoIL_2016 PRJNA217052 (172); BrooksB_2017 PRJNA376566 (5)
ChuDM_2017 PRJNA322188 (65); DeFilippisF_2019 PRJNA340216 (97)
DhakanDB_2019 PRJNA397112 (110); FengQ_2015 PRJEB7774 (61)
GuptaA_2019 PRJNA397112 (30); HMP_2012 PRJNA48479 (147)
HMP_2019_ibdmdb PRJNA398089 (426); HMP_2019_t2d PRJNA497499 (46)
HanniganGD_2017 PRJNA389927 (28); HansenLBS_2018 PRJNA491335 (207)
Heitz-BuschartA_2016 PRJNA289586 (26); IjazUZ_2017 PRJEB18780 (37)
KarlssonFH_2013 PRJEB1786 (43); KaurK_2020 PRJNA531203 (31)
KeohaneDM_2020 PRJEB36820 (117); KosticAD_2015 PRJNA231909 (89)
LifeLinesDeep_2016 EGAS00001001704/EGAD00001001991 (1135)
LokmerA_2019 PRJEB27005 (57); NagySzakalD_2017 PRJNA379741 (50)
Obregon-TitoAJ_2015 PRJNA268964 (57); PasolliE_2019 PRJNA485056 (112)
PehrssonE_2016 PRJNA300541 (191); QinN_2014 PRJEB6337 (114)
RampelliS_2015 PRJNA278393 (38); RaymondF_2016 PRJEB8094 (36)
RubelMA_2020 PRJNA547591 (86); SankaranarayananK_2015 PRJNA299502 (18)
ShaoY_2019 PRJEB32631 (1619); SmitsSA_2017 PRJNA392180 (27)
TettAJ_2019_a PRJNA529400 (68); TettAJ_2019_b PRJNA529124 (43)
TettAJ_2019_c PRJNA504891 (49); ThomasAM_2018a PRJNA447983 (24)
ThomasAM_2018b PRJEB27928 (27); ThomasAM_2019_c DRA006684 (40)
VatanenT_2016 PRJNA290380 (615); VincentC_2016 PRJNA297252 (196)
VogtmannE_2016 PRJEB12449 (52); WampachL_2018 PRJNA379120 (53)
WirbelJ_2018 PRJEB27928 (65); XieH_2016 PRJEB9576 (177)
YachidaS_2019 DRA006684 (251); YeZ_2018 PRJNA431482 (45)
YuJ_2015 PRJEB10878 (54); ZellerG_2014 PRJEB6070 (61)
ZhuF_2020 PRJEB29127 (81)
```

Known project overlaps include Gupta/Dhakan, Thomas/Wirbel, and Thomas/Yachida. Deduplicate at participant and biological-sample levels, not by accession string alone.

### 14.3 Required methods/reporting standards

- STORMS human-microbiome reporting checklist: https://www.nature.com/articles/s41591-021-01552-x and https://www.stormsmicrobiome.org/
- MIxS/MIMS contextual metadata: https://genomicsstandardsconsortium.github.io/mixs/
- STREAMS technical reporting checklist and frozen artifact: https://www.nature.com/articles/s41564-025-02186-2 and https://doi.org/10.5281/zenodo.15014818
- Compositional analysis rationale: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full
- Protocol bias and multiplicative error: https://elifesciences.org/articles/46923
- Extraction-method effect: https://pubmed.ncbi.nlm.nih.gov/28967887/
- Contamination-control method: https://pmc.ncbi.nlm.nih.gov/articles/PMC6298009/

Each report/model artifact records collection/storage, extraction kit and lot, library prep, sequencer/read length, blank/mock/spike controls, raw and post-QC depth, host-depletion reference, profiler/database hashes, filters/zero handling, model checksum, training manifest/split, code/data DOI, and build date.

### 14.4 Absolute abundance option

Healthy microbial load varies widely and relative trade-offs can be artifacts. Add an optional spike-in, flow-cytometry/cells-per-gram, or validated qPCR/ddPCR lane before using absolute terms. Source: https://www.nature.com/articles/nature24460

The relative-only report remains useful, but its vocabulary is constrained to relative expansion/depletion and qualified detection.

---

## 15. Functional panels, microbial groups, and counterfactual models

### 15.1 Integrate the 2018 personalized in-silico microbiota concept as a research scenario engine

The motivating paper—Bauer & Thiele, “From metagenomic data to personalized in silico microbiotas: predicting dietary supplements for Crohn's disease”—is directly relevant, but it is computational hypothesis generation rather than a validated treatment selector: https://www.nature.com/articles/s41540-018-0063-2

The study used 26 controls and 28 selected dysbiotic, newly diagnosed pediatric Crohn specimens from `SRP057027`, mapped them into a 773-strain metabolic reconstruction/BacArena environment, captured on average 73.5% of input relative abundance, simulated SCFA behavior, and predicted individualized nutrient additions. Suggested sets ranged from 1 to 55 metabolites, with a median of 19; pectin appeared for 17 of 24 model-responsive patients. Four modeled patients had no predicted solution. The authors explicitly said the dietary predictions require nutritional-trial validation.

Implement this as a separate `COUNTERFACTUAL_FLUX_RESEARCH` lane:

```yaml
scenario_model:
  model_id: ...
  claim_level: flux_predicted
  input_observation_refs: [...]
  compatibility:
    body_site: stool
    taxonomy_namespace: ...
    mapping_coverage: ...
  reconstruction:
    resource: AGORA2 | pinned_BacArena_reproduction
    version: ...
    digest: sha256:...
  environmental_constraints:
    baseline_diet_source: measured | questionnaire | generic_unknown
    nutrient_bounds: ...
  counterfactual:
    exact_substrate: ...
    modeled_change: ...
    uncertainty: ...
  validation_status: mechanistic_unvalidated
  intervention_binding: prohibited
```

Requirements:

1. Reproduce the original cohort/model before extending it.
2. Report community mapping coverage and excluded taxa/reactions.
3. Keep predicted flux separate from gene capacity and measured LC-MS/GC-MS metabolites.
4. Include diet constraints; a generic “rich diet” is a sensitivity scenario, not personal baseline.
5. Run ensembles across plausible abundance, diet, and model uncertainty.
6. Compare predictions with measured stool/plasma metabolites when supplied, without calling stool and plasma equivalent.
7. Never convert a predicted substrate into a participant dosing instruction.
8. A modeled scenario remains lane E until the exact protocol has prospective human validation.
9. Publish when the model returns no solution or mutually conflicting scenarios.
10. Do not use the model to optimize FMT donors or select added organisms without separate safety/regulatory/clinical validation.

This same typed scenario interface can later host AGORA2/MICOM. The computation may explore “what if” diet constraints, but the report must label every result **modeled, not measured or proven**.

### 15.2 Universal function-to-evidence rule

Each of the 25 current panels reports **genomic capacity**. For every panel, the intervention engine must ask four questions in order:

1. Was this molecule/activity directly measured? If no, the report cannot state its concentration or host exposure.
2. Is direction clinically interpretable in this population? If no, display reference position neutrally.
3. Is there human intervention evidence for this exact measured/capacity finding, rather than merely for a related disease? If no, use mechanistic context or `no_supported_targeted_intervention`.
4. Does the subject match the intervention population and safety requirements? If no, downgrade or suppress.

Gene capacity may retrieve an educational research summary whose internal source record contains a studied protocol. It cannot promote that protocol to a personalized action, clinician route, or recommendation and cannot select a drug, antimicrobial, probiotic strain, or dose.

### 15.3 Cross-engine contradictions

The report must actively search for and render contradictions such as:

- high butyrate gene capacity with low abundance of curated butyrate-producing taxa;
- high carrier abundance but incomplete pathway capacity;
- TMA capacity with no measured plasma TMAO;
- vitamin-synthesis capacity with unknown host vitamin status;
- virulence gene fragments without a resolved carrier/pathotype;
- methane gene nondetection with constipation symptoms or a positive clinical breath test;
- disease-profile association direction opposed by the matched-reference finding.

Contradictions widen uncertainty and block simplified action language. They are not averaged away.

### 15.4 Actionability map for all 25 functional panels

Evidence letters refer to Section 7 lanes. They describe source type, not predicted benefit. `General food pattern` means educational human evidence that still requires population/safety matching; it is not an automatic prescription from the percentile. Doses below are research-extraction facts for the internal registry and clinician audit; the participant renderer applies Section 7.8 dose suppression.

| # | Panel | Default render and permitted evidence | Confirmation, counterevidence, and safety |
|---:|---|---|---|
| 1 | Butyrate production | Low capacity may retrieve a gradual, varied fermentable-plant-food/resistant-starch evidence card; high capacity needs no action. RS2 potato starch was titrated 12→24→48 g and then 48 g/day for seven days; butyrate rose only in baseline responders, with *R. bromii* predicting response: [Venkataraman 2016](https://pubmed.ncbi.nlm.nih.gov/27357127/). AXOS 10.4 g/day for four weeks altered taxa without metabolic-marker improvement: [Kjølbæk](https://pubmed.ncbi.nlm.nih.gov/30827722/). | Capacity is not fecal/systemic butyrate. Responses are heterogeneous: [Wastyk](https://pubmed.ncbi.nlm.nih.gov/34256014/). Suppress generic fiber escalation for obstruction/stricture, severe flare, or prescribed restriction. Do not promise species restoration. |
| 2 | Propionate production | `inform_only`. In 60 overweight adults, 10 g/day inulin-propionate ester for 24 weeks reduced weight gain/adiposity versus inulin: [Chambers 2015](https://pubmed.ncbi.nlm.nih.gov/25500202/). | A 270-person 12-month trial did not reduce weight gain; adjusted difference favored control: [iPREVENT](https://pubmed.ncbi.nlm.nih.gov/40663640/). Experimental delivery ingredient, not correction of stool gene capacity. |
| 3 | BCAA biosynthesis | `not_actionable`; link routine metabolic confirmation when context warrants. | Microbial BCAA pathways correlated with serum BCAA/insulin resistance, while causal *P. copri* work was in mice: [Pedersen 2016](https://pubmed.ncbi.nlm.nih.gov/27409811/). Never restrict protein/BCAAs from DNA; use HbA1c/glucose/lipids/BP/waist clinically. |
| 4 | TMA/TMAO potential | `inform_only`; may ask about diet and renal function. Red-meat feeding raised measured plasma TMAO and withdrawal lowered it: [Wang 2019](https://pubmed.ncbi.nlm.nih.gov/30535398/). | Comparator-dependent feeding evidence is inconsistent: [2025 review](https://pubmed.ncbi.nlm.nih.gov/40419218/). Probiotic meta-analysis was null: [eight RCTs](https://pubmed.ncbi.nlm.nih.gov/35871952/); a very-high-dose dialysis trial was also null: [trial](https://pubmed.ncbi.nlm.nih.gov/29651635/). DMB inhibition is animal-only: [study](https://pubmed.ncbi.nlm.nih.gov/26687352/). Do not restrict essential choline or treat capacity as plasma TMAO/ASCVD risk. |
| 5 | Imidazole propionate | `not_actionable`; optional targeted LC-MS research measurement. | Association and metformin-interference work is mechanistic/observational, not a treatment trial: [Koh et al.](https://pubmed.ncbi.nlm.nih.gov/32783890/). Never alter histidine or metformin from the panel. |
| 6 | Secondary bile-acid production | `inform_only`; high/low `bai` capacity is not independently actionable. | Short controlled animal- versus plant-based diet changed bile-tolerant organisms and fecal deoxycholate but did not establish disease treatment: [David 2014](https://pubmed.ncbi.nlm.nih.gov/24336217/). Targeted bile-acid LC-MS is required to measure levels. |
| 7 | Bile salt hydrolase | Contextual/U-shaped; `inform_only`. Exact *L. reuteri* NCIMB 30242 studies in hypercholesterolemic adults may be shown only for that phenotype/product: [yogurt trial](https://pubmed.ncbi.nlm.nih.gov/22067612/); [2.9×10^9 CFU twice-daily capsule trial](https://pubmed.ncbi.nlm.nih.gov/22990854/); [safety](https://pubmed.ncbi.nlm.nih.gov/22561556/). | Total BSH genes do not predict lipid benefit or justify another strain/product. Standard lipid evaluation/care remains external. |
| 8 | Indole production | `not_actionable`; neutral/contextual explanation. | Indoles include potentially favorable AhR ligands and precursors of indoxyl sulfate; renal clearance shapes host exposure. Gene potential cannot resolve the net effect. |
| 9 | Indole-3-propionate | `not_actionable`; optional targeted LC-MS research measurement. | No replicated human intervention selectively raises IPA on the basis of stool `fld` capacity; *C. sporogenes* intervention work is mainly preclinical. |
| 10 | p-Cresol production | General population: `not_actionable`. In established CKD plus measured serum p-cresyl sulfate: clinician-review research evidence. A six-week SYNERGY crossover used 15 g/day mixed prebiotic plus a nine-strain probiotic escalating 45→90 billion CFU/day and reduced PCS but not indoxyl sulfate: [trial](https://pubmed.ncbi.nlm.nih.gov/26772193/). | A 16-week resistant-starch CKD trial lowered PCS but not other biomarkers/vascular function and reduced several alpha-diversity metrics: [Headley](https://pubmed.ncbi.nlm.nih.gov/39362281/). Surrogate evidence does not authorize action for stool genes or without renal review. |
| 11 | GABA production | `not_actionable`. | No validated human trial maps total stool GABA-biosynthesis capacity to luminal output, CNS exposure, symptoms, or “psychobiotic” selection. |
| 12 | Dopamine/tyramine production | Default `inform_only`; levodopa or monoamine-oxidase-inhibitor use can trigger a medication-safety clinician card. | Bacterial tyrosine decarboxylase abundance was associated with levodopa exposure/dose: [Maini Rekdal](https://pubmed.ncbi.nlm.nih.gov/30659181/). Inhibition evidence is mouse/preclinical: [Rekdal 2019](https://pubmed.ncbi.nlm.nih.gov/31196984/). No approved microbiome inhibitor/eradication protocol. |
| 13 | Histamine production | DNA alone: `not_actionable`. Diagnosed IBS may retrieve a clinician/dietitian-guided limited low-FODMAP card, not a chronic “low-histamine” prescription. | A high-histamine IBS subgroup and *K. aerogenes hdc* variant were identified, with causal treatment work mainly in mice: [De Palma](https://pubmed.ncbi.nlm.nih.gov/35895832/). A three-week low-FODMAP randomized study improved symptoms and lowered urinary histamine: [McIntosh](https://pubmed.ncbi.nlm.nih.gov/26976734/). Restriction needs reintroduction/personalization and may lower bifidobacteria: [ACG](https://pubmed.ncbi.nlm.nih.gov/33315591/). |
| 14 | Tryptamine production | `not_actionable`. | No human intervention validates metagenomic tryptamine capacity as a target. Abundance does not determine luminal concentration, receptor exposure, motility, or symptoms. |
| 15 | Vitamin B12 synthesis | `inform_only`; explicitly state that host B12 status was not measured. | Confirm suspected deficiency with serum B12 and, when appropriate, methylmalonic acid; supplementation follows clinical/diet/malabsorption context: [NIH ODS](https://ods.od.nih.gov/factsheets/VitaminB12-HealthProfessional/). Colonic microbial synthesis is not host bioavailability. |
| 16 | Vitamin K2 synthesis | `inform_only`; no K2 action from microbial genes. | Gene capacity does not measure coagulation or sufficiency. Warfarin requires stable intake and clinician management: [NIH ODS](https://ods.od.nih.gov/factsheets/VitaminK-HealthProfessional/). |
| 17 | Folate synthesis | `inform_only`; no folate dose from microbial genes. | Assess diet, medication, pregnancy, and blood status clinically. High folic-acid intake can obscure B12 deficiency: [NIH ODS](https://ods.od.nih.gov/factsheets/Folate-HealthProfessional/). |
| 18 | Riboflavin synthesis | `inform_only`. RIBOGUT tested 50 or 100 mg/day for two weeks in 105 healthy adults; individual doses did not significantly increase *F. prausnitzii*, though pooled groups showed higher butyrate: [trial](https://pubmed.ncbi.nlm.nih.gov/35943883/). | A 100 mg/day, three-week Crohn single-arm study did not materially alter MGS diversity/taxonomy/pathways: [RISE-UP](https://pubmed.ncbi.nlm.nih.gov/31873717/). Do not convert into a microbiome-directed vitamin recommendation. |
| 19 | Biotin synthesis | `inform_only`; specifically block hair-loss/alopecia biotin suggestions from microbial genes. | Biotin deficiency is rare: [NIH ODS](https://ods.od.nih.gov/factsheets/Biotin-HealthProfessional/). High-dose supplements can interfere with assays including troponin: [FDA](https://www.fda.gov/medical-devices/in-vitro-diagnostics/biotin-interference-troponin-lab-tests-assays-subject-biotin-interference). |
| 20 | Hydrogen-sulfide production | DNA alone: `not_actionable`. A bismuth study in ten healthy adults used 524 mg four times daily for 3–7 days and reduced ex-vivo fecal H₂S, with no disease endpoint: [Suarez](https://pubmed.ncbi.nlm.nih.gov/9558280/). | A small open-label eight-week UC 4-SURE diet had both response and worsening: [trial](https://pubmed.ncbi.nlm.nih.gov/35451489/). No generic sulfur restriction. Bismuth cards require salicylate/bleeding, anticoagulant, pregnancy, pediatric viral-illness, and toxicity gates. |
| 21 | Methane production | Stool `mcrA`/methanogen results are `inform_only`. Constipation/bloating plus standardized breath methane may unlock clinician review. | North American consensus uses methane ≥10 ppm: [consensus](https://pubmed.ncbi.nlm.nih.gov/28323273/); [ACG SIBO guidance](https://pubmed.ncbi.nlm.nih.gov/32023228/). A 31-person methane-positive IBS-C study used neomycin 500 mg twice daily plus rifaximin 550 mg three times daily for 14 days under an older threshold and improved some constipation symptoms, not pain: [trial](https://pubmed.ncbi.nlm.nih.gov/24788320/). Never automate prescription advice or diagnose IMO from stool. |
| 22 | β-Glucuronidase | `not_actionable`; block “estrogen detox” and calcium-D-glucarate claims. | Selective bacterial GUS inhibition for irinotecan toxicity is preclinical: [Wallace](https://pubmed.ncbi.nlm.nih.gov/21051639/); [Bhatt](https://pubmed.ncbi.nlm.nih.gov/32170007/). Oncology management is required for irinotecan toxicity. |
| 23 | Oxalate degradation | Low capacity/taxon nondetection is `inform_only`. Stone history plus 24-hour urine abnormalities routes to standard clinician evaluation; no *O. formigenes* probiotic recommendation. | Oxabact colonized but missed its plasma-oxalate endpoint after 52 weeks: [ePHex](https://pubmed.ncbi.nlm.nih.gov/35552824/). Lactobacillus/Bifidobacterium trial was null: [trial](https://pubmed.ncbi.nlm.nih.gov/34129232/); controlled diet outperformed probiotics: [Lieske](https://pubmed.ncbi.nlm.nih.gov/20736987/); [AUA guideline](https://www.auanet.org/guidelines-and-quality/guidelines/medical-management-of-kidney-stones). |
| 24 | Urease | `not_actionable`; genes do not diagnose hyperammonemia or hepatic encephalopathy. | In diagnosed recurrent HE, rifaximin 550 mg twice daily for six months reduced breakthrough events, with >90% also on lactulose: [Bass](https://pubmed.ncbi.nlm.nih.gov/20335583/). This is diagnosis-dependent guideline care, never stool-urease treatment: [AASLD](https://www.aasld.org/practice-guidelines/hepatic-encephalopathy). |
| 25 | Virulence (`clb/pks`, `bft`, `fadA`) | `clinician_review` for evidence interpretation and independent CRC-screening reminder; no eradication card. | Require full-locus breadth/depth, strain attribution, contamination control, technical replicate, and preferably absolute quantification. Carriage potential is not toxin expression or cancer. No human RCT proves eliminating these organisms prevents CRC: [review](https://pmc.ncbi.nlm.nih.gov/articles/PMC10308883/); screening remains conventional: [USPSTF](https://www.uspreventiveservicestaskforce.org/uspstf/recommendation/colorectal-cancer-screening). |

### 15.5 Actionability map for all 10 microbial groups

| Group | Required v04 behavior | Human evidence/context |
|---|---|---|
| Butyrate producers | Low relative abundance may retrieve the same gradual varied-fiber/resistant-starch evidence as panel 1, with all safety gates. Never promise species replacement. A high-pathway/low-carrier result is discordance, not failure. | Baseline degraders predict heterogeneous resistant-starch response: [Venkataraman](https://pubmed.ncbi.nlm.nih.gov/27357127/); responses to fiber/fermented foods differ: [Wastyk](https://pubmed.ncbi.nlm.nih.gov/34256014/). |
| Health-associated core | Rename `reference-associated taxa`. Low is not deficiency; default `inform_only`. | Geography, age, stool form, medications, protocol, and database determine the distribution; there is no universal core: [Dutch Microbiome Project](https://www.nature.com/articles/s41586-022-04567-7). |
| Probiotic-genus species | Low is not deficiency; high can reflect recent food/supplement transit, not durable colonization. Never select a product at genus/species level. | An 11-strain product showed person-, strain-, and body-region-specific colonization: [Zmora](https://pubmed.ncbi.nlm.nih.gov/30193112/). After antibiotics it delayed native recovery relative to spontaneous recovery, while autologous FMT recovered fastest: [Suez](https://pubmed.ncbi.nlm.nih.gov/30193113/). AGA evidence is formulation/indication specific: [guideline](https://pubmed.ncbi.nlm.nih.gov/32531291/). |
| Mucin degraders | Remove one favorable/unfavorable aggregate. Separate *A. muciniphila*, *R. gnavus*, *Bacteroides*, *Barnesiella*, and other taxa; interpretation is context-specific. | Live/pasteurized *A. muciniphila* MucT pilot: [Depommier](https://pubmed.ncbi.nlm.nih.gov/31263284/); later AKK-WST01 had no overall weight/HbA1c benefit, only exploratory baseline subgroup signals: [Zhang](https://pubmed.ncbi.nlm.nih.gov/39879980/); a 2026 study missed its primary insulin-sensitivity endpoint: [Suenaert](https://pubmed.ncbi.nlm.nih.gov/42343233/); another exact post-weight-loss protocol reduced weight regain: [Mount](https://pubmed.ncbi.nlm.nih.gov/42120725/). Show together; do not infer response from endogenous abundance without independent replication. |
| Oral-origin organisms | Do not equate relative enrichment with oral invasion. Require oral-health context, contamination checks, and preferably absolute load. Ordinary dental evaluation only if overdue/symptomatic. | Relative oral taxa can rise because resident gut biomass falls: [Liao 2024](https://www.nature.com/articles/s41564-024-01680-3). Periodontal treatment lowered stool *F. nucleatum* in CRC patients but did not prove cancer/community benefit: [trial](https://pubmed.ncbi.nlm.nih.gov/34887459/). |
| Opportunists/pathobionts | Split ecological pathobionts from clinically virulent strains. Species-level relative abundance cannot trigger infection, eradication, antibiotic, botanical, or FMT action. | Require symptoms, absolute evidence when relevant, strain/virulence/AMR resolution, and appropriate culture/clinical testing. *R. gnavus* must not share an infection rule with Enterobacterales. |
| Hexa-acylated LPS carriers | `not_actionable`; do not call “endotoxemia” or suggest “LPS detox.” | Taxonomy cannot reliably establish lipid-A structure/expression, translocation, or circulating endotoxin. No validated intervention is selected by this group score. |
| Sulfate reducers | `inform_only`; no eradication or generic sulfur restriction. | *Bilophila* and *Desulfovibrio* can be commensal at low abundance. Apply H₂S-panel cautions; current human evidence does not support abundance-triggered antimicrobials. |
| Methanogens | Use only methane-panel logic. | Stool abundance is not intestinal methanogen overgrowth; require constipation phenotype and standardized breath methane before clinician review. |
| Polyphenol metabolizers | Presence is not a functional result. `inform_only` unless a standardized substrate challenge plus urine/plasma metabolomics is supplied. | Equol-producer classification used 2×250 mL soy milk/day for three days and urinary S-equol:daidzein ratio: [Setchell](https://pubmed.ncbi.nlm.nih.gov/16857839/). Urolithin-producing strains remain in-vitro candidates: [study](https://pubmed.ncbi.nlm.nih.gov/36840624/). Direct urolithin A 1,000 mg/day for four months in adults 65–90 improved a secondary endurance measure but missed primary walking/ATP endpoints and bypasses rather than restores microbes: [RCT](https://pubmed.ncbi.nlm.nih.gov/35050355/). |

### 15.6 Natural-antimicrobial evidence block

The system explicitly suppresses an anecdote-to-protocol path for berberine, allicin, oregano oil, black-seed oil, GI Microb-X, or a multi-product “kill phase.”

- The herbal SIBO report was retrospective, nonrandomized, patient-choice chart review. It compared four weeks of rifaximin 400 mg three times daily with exact two-product herbal pairs; breath-test normalization was 34% versus 46%, not statistically significant, with no controlled stool-eradication or symptom endpoint: https://pubmed.ncbi.nlm.nih.gov/24891990/
- BRIEF-SIBO is a protocol, not outcome evidence: https://pubmed.ncbi.nlm.nih.gov/36873985/
- A historical 400 mg berberine trial concerned acute enterotoxigenic-*E. coli*/cholera diarrhea, not chronic decolonization: https://pubmed.ncbi.nlm.nih.gov/3549923/
- Oregano/*Klebsiella*, allicin, and black-seed findings are predominantly in-vitro or animal work and remain lane E.
- The subject's experience is an `ExposureEvent`, never a source assertion or causal conclusion.

### 15.7 Compiled action mapping

```yaml
action_mapping_id: ...
trigger_selector_ref: ...
intervention_protocol_ref: ...
evidence_assertion_refs: [...]
against_evidence_assertion_refs: [...]
safety_rule_refs: [...]
confirmation_rule_refs: [...]
participant_claim_release_ref: ... | null
```

This is a compiled foreign-key mapping only. It cannot duplicate or override protocol identity, source facts, evidence grade, safety, regulatory status, or rendering policy held in their canonical records.

Store three independent propositions:

1. `evidence_signal_associated_with_state`
2. `evidence_intervention_changes_signal`
3. `evidence_intervention_changes_clinical_outcome`

No proposition inherits the truth of another.

### 15.8 Priority exact-protocol records

These are the first internal evidence-registry records to encode and full-text verify. Values shown are source-reported studied protocols, not personalized instructions. Prescription and FMT doses are stripped from participant output; every other participant dose is suppressed unless a claim-specific release record passes the rule in Section 7.8.

| Condition/context | Exact protocol to capture | Primary result scope | Mandatory limitation/source |
|---|---|---|---|
| Clinically defined PACS | SIM01, 10 billion CFU in sachets twice daily for six months; capture the complete proprietary strain/prebiotic formula from supplement/label rather than guessing missing components | 463 randomized adults; improved several six-month symptom-alleviation outcomes; adverse-event frequency similar | Single-center Hong Kong exact product; immune-compromised/pregnant/breastfeeding patients excluded; commercial/IP conflicts; no stool-profile responder validation: [trial](https://pubmed.ncbi.nlm.nih.gov/38071990/) |
| Alopecia areata | *L. rhamnosus* Bths-08 CECT30580 + *B. longum* Bths-06 CECT30616, 1:1, total 10^9 CFU/day for 24 weeks, adjunct to intralesional triamcinolone every four weeks | 26 randomized/19 completed; primary plaque/SALT and other between-group hair outcomes were not statistically significant; no gut-microbiome change | Tiny adjunct pilot, not evidence for monotherapy or FMT: [trial](https://www.mdpi.com/2079-9284/11/4/119) |
| Androgenetic alopecia | *L. plantarum* DCn_07 CECT30102 + DCn_06 CECT30103 + *L. pentosus* DCn1_05 CECT30104; capture source-reported dose/formulation; 16 weeks | 136 randomized/115 completed; mixed hair-cycle endpoints, several nonsignificant variables | Industry/product-specific; no stool responder signature or class transfer: [trial](https://pubmed.ncbi.nlm.nih.gov/39275216/) |
| Alzheimer disease | 200 mL/day exact four-organism probiotic milk for 12 weeks; source-level identity/CFU must be captured before activation | One 60-person trial reported MMSE/metabolic changes | Later severe-AD study was negative and conservative pooled evidence is not robust: [positive trial](https://pubmed.ncbi.nlm.nih.gov/27891089/); [negative trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6104449/); [meta-analysis](https://pubmed.ncbi.nlm.nih.gov/34957177/) |
| Diagnosed IBS | Heat-inactivated *B. bifidum* MIMBb75, 10^9 nonviable cells daily for eight weeks | 443 adults; composite response 34% versus 19% | Exact nonviable product; not interchangeable with live *B. bifidum* or selected by abundance: [trial](https://pubmed.ncbi.nlm.nih.gov/32277872/) |
| Diagnosed IBS | *B. infantis* 35624; three dose arms for four weeks, with signal at 10^8 CFU/day but not lower/higher doses | 362 women; dose-specific symptom result | Illustrates nonmonotonic dose and strain specificity; guidelines still do not endorse generic probiotics: [trial](https://pubmed.ncbi.nlm.nih.gov/16863564/); [AGA](https://pubmed.ncbi.nlm.nih.gov/32531291/) |
| Parkinson disease with constipation | Exact multistrain product for four weeks; full formula/dose from source | 72 patients; spontaneous bowel movements improved; endpoint was constipation | No motor/neuroprotection claim and no Parkinson-score trigger: [trial](https://pubmed.ncbi.nlm.nih.gov/33046607/) |
| Insulin resistance/overweight | Pasteurized *A. muciniphila* MucT, 10^10 cells/day for three months | Small pilot metabolic-surrogate effects | Later exact-product/strain trials have null primary or subgroup-dependent results; show together: [pilot](https://pubmed.ncbi.nlm.nih.gov/31263284/); [2025 trial](https://pubmed.ncbi.nlm.nih.gov/39879980/); [2026 primary-null trial](https://pubmed.ncbi.nlm.nih.gov/42343233/) |
| Diagnosed IBS diet | Limited low-FODMAP restriction, ordinarily 4–6 weeks, then reintroduction and personalization; dietitian support | Symptom-management guideline pathway | Not microbiome normalization; suppress with eating-disorder risk/malnutrition and personalize: [AGA diet guidance](https://pubmed.ncbi.nlm.nih.gov/35337654/) |
| Recurrent CDI in adults | VOWST label: four capsules daily for three days after antibacterial treatment; REBYOTA label: single 150 mL rectal dose 24–72 hours after antibiotics | Prevention of recurrent CDI, not treatment of active infection | Exact FDA indications/products, age ≥18, clinician-only; never triggered by dysbiosis/profile: [VOWST](https://www.fda.gov/news-events/press-announcements/fda-approves-first-orally-administered-fecal-microbiota-product-prevention-recurrence-clostridioides); [REBYOTA label](https://www.fda.gov/media/163587/download) |

Every source undergoes exact formula/strain accession verification before activation. If a paper reports only a brand name or incomplete composition, mark `exact_protocol_reproducibility: incomplete` rather than filling the gaps from retail listings.


---

## 16. Migration and backward compatibility

1. Every valid v03 `Observation`, assertion, module, profile, score, confidence vector, coverage field, availability mask, and shape result retains its meaning. `Observation` remains schema 3.x; v04 adds separate records rather than silently changing it.
2. Every valid v03 profile file compiles unchanged under v04.
3. v04 adds records; it does not silently rename or reinterpret v03 fields.
4. Preserve every existing v03 JSON path at its original location for at least two release cycles and dual-write an immutable audit copy under `legacy_v3`; `legacy_v3` does not replace the original path. Check in `compatibility/legacy_v3_path_map.yaml` with every old path, canonical v04 path, dual-write type/coercion, deprecation release, and fixture. An unmapped legacy path fails compatibility tests.
5. For identical inputs and registry locks, score equality within a frozen floating-point tolerance of \(10^{-12}\) applies only when v03 analytical and metadata preconditions were actually satisfied. Corrected AGA-sex, pediatric-reference, host-read, and other explicitly listed conformance defects intentionally change invalid legacy statuses and must appear in the migration ledger.
6. No disease-association assertion may be auto-converted into an intervention assertion.
7. Existing “beneficial,” “opportunist,” and microbial-group labels may seed a human curation queue, but never treatment records.
8. Adult references support adult context only. They cannot silently produce pediatric normative labels.
9. Run v03 and v04 side-by-side for all five existing specimens and publish a field-level migration ledger.
10. Rebuild a report from raw input, manifest, source snapshots, and `registry.lock`; output digests must match.
11. The intervention subsystem can fail without preventing the measurement/profile report from completing.
12. Any profiler/database upgrade requires complete reference reprocessing or an independently validated bridge; name crosswalks alone are insufficient.

---

## 17. Validation and acceptance tests

### 17.1 Analytical and biological validation gates

- Analytical repeatability/reproducibility across extraction, library, lane, and run.
- Taxon/workflow-specific LOD95, LOQ, linearity, specificity, ambiguity, and contamination using mock, dilution, spike-in, blanks, and replicates.
- Held-out study, country, age band, sex, stool form, medication, and batch performance.
- Reference false-flag rate and interval coverage on healthy/control holdouts.
- Disease-profile AUROC/AUPRC only on labeled independent cohorts, with calibration and confidence intervals.
- Cross-disease specificity matrix including off-diagonal signal and `diffuse`/`specific` behavior.
- Longitudinal minimal detectable change and same-person repeat variation.
- Independent clinical review of every human-facing interpretation and intervention assertion.

Without labeled participant-grouped study-held-out validation, a disease result remains `uncalibrated_pattern_concordance`, never probability.

Check in a versioned `validation/validation_policy_v04.yaml` before opening any final holdout. It must contain numeric, assay/profile/age-stratum-specific limits for repeatability, LOD95/LOQ, ambiguity and contamination, matched-reference per-feature and per-report false-flag burden, interval coverage/tolerance, maximum age MAE, calibration slope/intercept, subgroup disparity, and OOD false-accept/abstention. Distinguish study/country/batch external validation from subgroup evaluation by age, sex, stool form, and medication. A missing limit, post-hoc threshold change, or failed holdout gate makes the relevant function experimental/not deployable.

### 17.2 Core machine acceptance tests

| ID | Given | Required result |
|---|---|---|
| A01 | Valid v03 profile whose analytical/metadata preconditions were satisfied, identical input and lock | v03 score unchanged |
| A02 | Feature without namespace/version/body site | Compiler failure |
| A03 | Nondetection below feature-depth gate | `indeterminate`, never absent |
| A04 | Core/common reference feature, adequate power, exact/prevalidated fallback, stable across studies, nondetected | `commonly_detected_not_detected` |
| A05 | Rare feature detected with high-confidence identity and blank/contamination pass | `uncommon_detected`, not automatically harmful |
| A06 | Adult reference applied to a 12-year-old | Comparator unavailable; no adult low/high/missing label |
| A07 | *E. coli* detected, no pathotype/toxin evidence | No pathogen/infection label |
| A08 | High relative *Klebsiella*, no absolute load | `relative_expansion`; “overgrowth” prohibited |
| A09 | FASTQ-only sample plus metabolite assertion | `not_measured`; no imputation |
| A10 | Low probiotic species, strain-specific trial | No strain intervention match |
| A11 | Multi-ingredient trial | Evidence cannot be split across components |
| A12 | Different product/dose/formulation from source | Partial/class-only applicability, not exact protocol |
| A13 | Disease resemblance without diagnosis | Drug/FMT remains evidence-only |
| A14 | Confirmed indication with guideline-backed drug | Clinician-review lane; never self-treatment |
| A15 | Minor and adult-only supplement study | Outside population or suppressed |
| A16 | Immunocompromise/central-line status unknown | Live probiotic/FMT fails closed |
| A17 | Positive and null/harm trials | All render; conflict status set |
| A18 | Guideline recommends against FMT | Lane X visible; option not promoted |
| A19 | “Dysbiosis” or high profile alone | Cannot trigger FMT |
| A20 | Minor donor sample in authorized research/clinician mode | `not_supported_under_product_policy_minor`; no per-lane result or eligibility output |
| A21 | Longitudinal samples on different databases | No quantitative delta unless reprocessed together |
| A22 | Relative abundance falls below detection after supplement | “Not detected at this depth,” never “eradicated” |
| A23 | Profile score falls after exposure | Temporal association only; no causal claim |
| A24 | Hundreds of comparisons | Headline requires FDR and stability gates |
| A25 | No exact intervention evidence | `no_supported_targeted_intervention` |
| A26 | Harm-only evidence | Warning renders; no promotion |
| A27 | Safety rule references undeclared/misspelled field | Unconditional compiler failure |
| A27b | Safety rule references declared nullable field without `on_unknown` | Compiler failure |
| A28 | Undeclared rule operator/unit | Compiler failure |
| A29 | Registry artifact digest changes | Reproducibility failure |
| A30 | Intervention registry unavailable | Screening/profile report still completes |
| A31 | Profile has partial evidence mask | Low-coverage score/status and widened interval |
| A32 | LC-MS result has wrong matrix/unit | No binding until validated conversion |
| A33 | Taxa/pathways from same participants | One evidence cluster, not two replications |
| A34 | Botanical evidence is class-level only | Exact product not recommended |
| A35 | Direction unstable leave-one-study-out | Not promoted as standout |
| A36 | “Common” taxon supported by one study | Context only; expected/missing suppressed |
| A37 | Retraction/correction appears | Registry blocks or flags by source policy |
| A38 | Prescription assertion lacks exact indication | Evidence-only/clinician lane, no action language |
| A39 | Probiotic CFU basis unspecified | Protocol reproducibility incomplete |
| A40 | Disease model is OOD | Visible OOD warning or abstention by frozen policy |
| A41 | Age preprocessing fitted before CV split | Test fails model build |
| A42 | Age model lacks frozen scaler/PCA/feature order | `not_deployable` |
| A43 | Minor lacks validated pediatric age calibration | No numeric age; insufficient support message |
| A44 | Actual chronological age/age assurance absent | Internal research estimate may exist without residual; participant output `not_computable` |
| A45 | Age estimate older than actual | No red/failure interpretation |
| A46 | Host reads removed by provider | `not_assessable_upstream_removed`, not pass |
| A47 | AGA score with sex missing | Abstention or declared two-stratum interval |
| A48 | SAMPLE2 has 17 profiles above the highlight threshold | All 17 highlighted and every other scored profile detailed under Section 10; never capped at eight |
| A49 | MetaPhlAn 3 and 4 inventories | Visibly separate counts/namespaces |
| A50 | Conflicting butyrate capacity and carrier abundance | Explicit contradiction card; no collapsed answer |
| A51 | Feature meets core and common posterior rules | One primary class: `core` by ordered precedence |
| A52 | One large study dominates pooled prevalence | Hierarchical/study-balanced result and heterogeneity shown; no pooled headline |
| A53 | Total reference `n` passes but positive `n`/study count fails | Prevalence may render; abundance percentile unavailable |
| A54 | Geography/stool-form absent and fallback not feature-validated | `unmatched_context_only` or comparator unavailable; no missing call |
| A55 | Denominator/zero-replacement digest changes | New reference-model version; reproducibility test fails old lock |
| A56 | HUMAnN capacity input | Labeled relative encoded potential; no absolute production/activity claim |
| A57 | Ambiguous/split taxonomy mapping | Assertion transfer blocked unless required resolution preserved |
| A58 | Age prediction outside validated range/support | `abstained_ood` with null numeric participant output |
| A59 | Infant-heavy training produces poor teen/adult stratum metrics | Model fails even if global MAE passes |
| A60 | Untouched conformal calibration fixture | Frozen 80%/95% coverage tolerances pass by supported stratum |
| A61 | Rare detection fails blank/ambiguity gate | Indeterminate/possible contaminant; no uncommon-present headline |
| A62 | Display filter changes but testing universe/data do not | Identical p/q values and `multiplicity_family_id` |
| A63 | Registry/profile/report sets differ from 33-ID manifest | Compiler or report build failure |
| A64 | Study module appears as independent profile | Compiler failure; no duplicate shape contribution |
| A65 | Participant is a minor with unvalidated disease transport | Adult-derived numeric percentile suppressed; internal research record only |
| A66 | Participant renderer receives prescription/FMT studied dose | Dose removed and clinician route rendered |
| A67 | Evidence summary is triggered by a profile, external lab, symptom, or verified condition | Every `trigger_ids[]` value resolves to a typed `TriggerRecord.trigger_id`; no finding-only foreign key |
| A68 | Self-reported condition matches a prescription/FMT indication selector | It remains unconfirmed and cannot unlock the route |
| A69 | Confirmed-condition record lacks source/date/credential or verified status | Compiler failure or selector does not match |
| A70 | Legacy combined safety status enters compiler | Compiler failure; canonical safety-signal, display-policy, and action-status fields required |
| A71 | Positive, null, harm, and guideline-against assertions coexist | `against_evidence[]` preserves every applicable category and reference; supportive evidence remains separate |
| A72 | Participant report/API requests donor prescreen | Per-lane object is absent and only “Clinical donor qualification was not performed” renders |
| A73 | Finding references a missing or nonreleaseable reference model as `exact` | Compiler failure; numeric comparator fields cannot render |
| A74 | Core/common taxon is qualified not detected | One primary `commonly_detected_not_detected` state; no low/depleted abundance state |
| A75 | Taxon is detected below LOQ or LOQ is not established | No low/high abundance percentile; identity/detection context may remain secondary |
| A76 | Identity probability is populated without a calibrated locked bundle, or identity gate is unsupported | Compiler failure or fail-closed nonpromotion |
| A77 | Age reconstruction project totals or row-level join do not equal the locked 8,641-row assertion | Training blocked with enumerated missing/changed/duplicated rows |

### 17.3 Report acceptance tests

1. Count parity: 25/25 functions, 10/10 groups, every profile status, every detected taxon, and every expectedness bucket agree with JSON.
2. No adult `healthy range` for minors without validated pediatric transport.
3. Every normative visual names reference, `n`, studies, stratum, version, and fallback.
4. No dynamic text truncation or literal ellipsis.
5. No unintended page with less than 30% content area.
6. Main text at least 9.5 pt, tables 8.5 pt, footnotes 7.5 pt.
7. WCAG-AA contrast plus grayscale/deuteranopia snapshot tests.
8. Embedded fonts, tagged structure, language, reading order, bookmarks, repeated headers, alt text, and clickable links.
9. Red cannot be assigned merely because a value is high; level and interpretation fields are independent.
10. Disease resemblance cannot render medication/FMT action language without a separate confirmed-indication gate.
11. Score-method prose is generated from the exact canonical implementation/configuration.
12. Every research-evidence summary exposes source, population, exact identity, primary endpoint/result, null/harm evidence, known safety state, and applicability; participant mode obeys the dose/redaction policy.

### 17.4 Shadow-run expectations for the supplied samples

- **SAMPLE2:** adult reference only if age/sex/country metadata satisfy the reference; show diffuse disease shape; highlight the butyrate capacity-versus-carrier contradiction; record prior supplement exposures without efficacy/eradication claims; no automatic antimicrobial or FMT recommendation.
- **A04 and A06:** if under 18 or age remains missing, suppress adult normative high/low/missing language, numeric pediatric age, adult treatment applicability, and participant-facing adult-derived disease percentiles. Unvalidated numeric OOD results may remain access-controlled internal research records only.
- **All five samples:** no FMT donor clearance. In access-controlled researcher/credentialed-clinician mode, subjects under 18 receive `not_supported_under_product_policy_minor` and adults receive only non-clearing, per-lane research-prescreen states. Participant reports receive no per-lane states and show only “Clinical donor qualification was not performed.”

---

## 18. Implementation plan and release gates

### Phase 0 — Freeze and govern

1. Freeze current v03 outputs and golden fixtures for all five samples.
2. Create independent version fields and `registry.lock` digests.
3. Adopt STORMS/MIxS/STREAMS metadata and source-claim ledger requirements.
4. Obtain documented FDA function-by-function classification/premarket-path determination (using Q-Submission feedback where appropriate), CLIA/state DTC analysis, FTC claims review, FMT legal/clinical review, pediatric/COPPA/privacy/IRB determinations.

**Gate:** no patient-facing disease concordance, intervention retrieval/selection, donor prescreen, or health interpretation of age until regulatory owners sign the function-by-function release inventory and required authorization/compliance path.

### Phase 1 — Matched reference and findings

5. Implement `ReferenceDefinition`, frozen manifests, and hierarchical stratum resolver.
6. Build empirical detection curves and analytical QC controls.
7. Implement prevalence posterior, log-ratio percentiles, study-block bootstrap, FDR, and leave-one-study-out stability.
8. Implement `ReferenceCall`, `Finding`, interpretation registry, and the seven organism buckets.
9. Build pediatric references only where the cohort/study gates pass; otherwise abstain.
10. Add optional absolute-load lane.

**Gate:** no expected-but-not-detected output until LOD, prevalence, matched-reference, and false-flag validations pass.

### Phase 2 — Intervention evidence

11. Implement exact `InterventionIdentity`, `InterventionProtocol`, PICO assertion, outcome ontology, and null/harm records.
12. Implement restricted safety-rule compiler and applicability matcher.
13. Seed evidence mappings for the exact 33-profile target set plus all 25 functional panels and 10 microbial groups in Sections 13 and 15.
14. Add retraction/correction monitoring, review owners, next-review dates, and immutable source snapshots.
15. Implement `no_supported_targeted_intervention` and evidence-only fallback.

**Gate:** clinical, regulatory, pharmacy, nutrition, and scientific review of every human-facing card; no free-text uncurated treatment generation.

### Phase 3 — Microbiome age

16. Reconstruct/freeze the age training manifest and eliminate participant/project overlap.
17. Reimplement leakage-free preprocessing and benchmark simple models versus TRPCA.
18. Validate by held-out study/country/age band and an independent external cohort.
19. Release adult model only for validated strata; release pediatric numeric output only after its separate gate.

**Gate:** deployable artifact completeness, external interval coverage, OOD behavior, and model/data license review.

### Phase 4 — Longitudinal and FMT firewall

20. Implement exposure timeline, compatibility, minimal-detectable-change, and paired uncertainty.
21. Implement the access-controlled per-lane `research_donor_prescreen` with prohibited language, no aggregate/rank/pass state, clinical-test gap list, and a participant serializer that cannot expose it.
22. Encode current FDA alerts and indication-specific FMT guideline status.

**Gate:** no donor ranking/clearance; no disease-profile-to-FMT edge; minors hard-stopped.

### Phase 5 — Report rebuild

23. Build the front visual atlas in the order specified in Section 10.
24. Add complete function/group/taxon/profile/intervention detail sections without renderer caps.
25. Implement universal typed visual grammar and reference context.
26. Add tagged accessible PDF, embedded fonts, bookmarks, clickable sources, and QA tests.

**Gate:** all JSON/PDF parity, accessibility, readability, overflow, and three-supplied-report regression tests pass.

### Phase 6 — Validation and staged release

27. Shadow-run all five current samples; document every difference from v03.
28. Publish analytical validation, reference false-positive rates, disease cross-specificity, age validation, and bounded limitations.
29. Release first to internal research mode, then clinician-review mode only after external review; participant mode requires the strictest claims/safety approval.
30. Schedule evidence-source review at least quarterly and immediately on FDA safety alerts, guideline changes, retractions, or product/formulation changes.

---

## 19. Launch blockers

The build is not complete until all are resolved:

1. Documented FDA classification and premarket/compliance path—not merely a Q-Submission strategy—for patient-facing disease matching, intervention retrieval/selection, donor prescreening, and health interpretation of age.
2. An appropriate CLIA-certified laboratory for individual health-assessment reporting unless counsel documents a valid nonclinical research workflow, plus state direct-access ordering/reporting and medical-practice review. Do not assume a blanket LDT exemption or automatic FDA coverage after the March 31, 2025 vacatur and September 19, 2025 text restoration.
3. FTC claim-by-claim substantiation review of report and marketing together.
4. FMT legal/clinical review confirming no output qualifies donors or promotes non-CDI FMT.
5. Pediatric, COPPA, consumer-health/genetic-privacy, retention/deletion, and research-IRB determinations.
6. Independent analytical validation and locked references before “normal,” “good,” “bad,” diagnostic, predictive, or treatment-responsive claims.
7. Security controls for raw FASTQ and derived health profiles: host-read quarantine/removal, encryption, least privilege, access audit, retention limits, deletion workflow, prohibition on secondary use without consent/IRB basis, a documented Health Breach Notification Rule applicability decision, and tested incident-notification clocks/duties.
8. If any function is a device, operational QMSR lifecycle controls—design/development, risk, validation, changes, cybersecurity, complaints, applicable reporting, and corrections/removals—not a documentation link alone.

Relevant current sources:

- [FDA CDS guidance PDF, January 2026](https://www.fda.gov/media/109618/download)
- [FDA CDS FAQ, June 2026](https://www.fda.gov/medical-devices/software-medical-device-samd/clinical-decision-support-software-frequently-asked-questions-faqs)
- [FDA Q-Submission guidance](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/requests-feedback-and-meetings-medical-device-submissions-q-submission-program)
- [FDA Quality Management System Regulation](https://www.fda.gov/medical-devices/postmarket-requirements-devices/quality-management-system-regulation-qmsr)
- [42 CFR 493.2](https://www.ecfr.gov/current/title-42/chapter-IV/subchapter-G/part-493/subpart-A/section-493.2)
- [CMS CLIA](https://www.cms.gov/medicare/quality/clinical-laboratory-improvement-amendments)
- [FTC Health Products Compliance Guidance](https://www.ftc.gov/business-guidance/resources/health-products-compliance-guidance)
- [COPPA final rule](https://www.federalregister.gov/documents/2025/04/22/2025-05904/childrens-online-privacy-protection-rule)
- [FTC Health Breach Notification Rule](https://www.federalregister.gov/documents/2024/05/30/2024-10855/health-breach-notification-rule)
- [45 CFR Part 46](https://www.hhs.gov/ohrp/regulations-and-policy/regulations/45-cfr-46/index.html)
