# OpenBiota gut intervention evidence atlas — implementation handoff

Compact edition · 24 September 2026 · Research integration for the existing recommendation engine

**Read this file as the coding instructions. Query the separate data package with code. Do not load the data package into the agent's conversation context.** This is a packaging and handoff correction to the existing research atlas; it is not a replacement application build specification.

## Supplied assets and how to use them

| Asset | Use |
|---|---|
| `OPENBIOTA_GUT_INTERVENTION_EVIDENCE_ATLAS.md` | This compact handoff: engine decisions, evidence rules, integration sequence, acceptance criteria and source acquisition instructions. |
| [OPENBIOTA_INTERVENTION_DATA.zip](sandbox:/workspace/scratch/1d3f677dbe34/OPENBIOTA_INTERVENTION_DATA.zip) | Machine-readable research asset. Extract locally; do not paste its contents into a prompt. |
| `OPENBIOTA_INTERVENTION_DATA/FORMAT.md` inside the archive | Complete file-format description, column definitions, interpretation rules and SQL/helper examples. Read this before writing the importer. |

Give the coding agent this Markdown and filesystem access to the downloaded/extracted archive. The attachment link is for downloading in this conversation; an agent running elsewhere needs the archive placed in its workspace or an accessible URL. No other earlier research document or giant Markdown file is required.

The archive contains `openbiota_intervention_data.sqlite`, `FORMAT.md`, `query.py`, and `manifest.json`. It preserves **87,010 data rows from 110 tables, 18 structured JSON objects, all original research notes and 115 source-catalogue entries**. These are heterogeneous records, not counts of independent studies or proven treatments. Public snapshots remain in the machine dataset to preserve the previously added identifiers and annotations; future refreshes use the original download links below.

The conversion is lossless: database records and archived narrative reconstruct the original 23,400,377-byte source exactly, verified by SHA-256. This checks preservation, not the scientific correctness of every source assertion. The large document no longer serves as mandatory agent context.

## Agent workflow

1. Read this handoff and the archive's `FORMAT.md`. Reuse the repository's existing architecture, identifiers and data locations; do not invent a new application directory hierarchy or rewrite unrelated report algorithms.
2. Extract the ZIP and open SQLite read-only. Read `metadata`, inspect `data_tables`, and use bounded record searches. Data import can process entire tables in code without sending their contents to an LLM.
3. Import the curated observations, ecological records, source-quality overlays and exact-material/strain identities into the existing evidence store. Import report-specific mappings as private regression fixtures. Treat ingredient aliases, genome indexes and publication metadata as supporting data, not efficacy observations.
4. Add separate adapters for the upstream data sources below. Preserve their original identities, study contexts and known source anomalies. Do not discard enrichment just because an original download is available.
5. Apply the engine decisions and acceptance checks in this handoff. Query the archived research notes and full source-catalogue entries only when implementing the relevant component.

Example bounded access, from the extracted directory:

```bash
python3 query.py tables --match 'Curated intervention'
python3 query.py search 'rosemary' --limit 5
python3 query.py search '"calcium" AND "glucarate"' --limit 5
python3 query.py record PI09914
python3 query.py sources --match 'MICOM'
python3 query.py notes 'Source-quality corrections' --limit 3
```

The helper runs read-only and emits JSON Lines. Its search results identify records to inspect; the search excerpt is not a substitute for the record and its source-quality/context notes. See `FORMAT.md` for direct SQL and decoded cell handling.

## Common data contract

`data_tables` identifies each source table and its original column labels. `records` provides one common envelope: unique package record ID, source table ID, original evidence ID when present, source row/line and `fields_json`. Values retain the original source strings and citations; this deliberately avoids converting unknown values or incompatible measurements into a fabricated universal biological schema. `structured_objects` contains original structured pathways and manifests. `source_blocks` retains the research rationale and overlays. `source_catalog` stores source access, scope and license notes. `records_fts` is a lexical retrieval index.

The original raw row is retained only for audit. Production adapters use `fields_json` plus source-specific column mapping. A record's evidence type, experimental system, source quality, exact material and direction must be resolved before it becomes an application action assertion. Existing source facts, locally derived annotations and new predictions remain distinguishable.


## Integration decisions

| Layer | Selection and required behavior |
|---|---|
| Inputs | Reuse versioned sample JSON and validated taxonomy/strain/function outputs, preserving method, reference version and uncertainty. Reprocess reads only for justified missing features; FASTQ belongs in bioinformatics, not a chat prompt. |
| Canonical evidence | PostgreSQL with typed relationships, full-text search and optional pgvector discovery. Exact identifiers precede semantic matches. A graph is a data model; a separate graph server needs a demonstrated query requirement. |
| Research maintenance | DeepEvidence-style primary-source research and BioCypher/BioChatter-style normalization. Proposed assertions retain source spans, experimental context and review status before becoming authoritative. |
| Metabolic scenarios | MICOM with compatible AGORA2 models and qualified extensions. Persist assumptions, coverage, solver identity and sensitivity. Outputs remain scenario predictions. |
| Learned responses | Evaluate McMLP and Venturelli miRNN within their supported intervention/endpoint domains. Separate response prediction, generation, imputation and causal estimation. |
| Synthetic data | MIDASim and MB-GAN adapters for robustness and distributional tests. Generated profiles never increase empirical evidence counts. |
| Planning | Deterministic conditional decision graph; OR-Tools CP-SAT when combination/order constraints justify it. The solver optimizes declared objectives, not biological truth. |
| Execution/explanation | Checkpointed workflow such as LangGraph, with a provider-neutral reasoning LLM selected by project-specific evaluation. Record retrieval, assumptions, model runs, contradictions and plan revisions. |
| Feedback | Versioned longitudinal exposure, context and measured outcomes. Update after observations; simulated follow-up never replaces observed state. |

Adapt interfaces to the existing repository before introducing dependencies. Preserve a usable evidence-only path when compatible model inputs are unavailable. Culture, enzyme, cell, animal, ex vivo and human findings remain available; their context governs the claim. Lack of a human trial alone does not erase a candidate.

## Interpretation-layer contracts

The companion SQLite file is the preserved source corpus, with its physical format in FORMAT.md. These logical objects describe the engine's interpretation layer; missing source fields remain missing.

| Object | Required fields/semantics |
|---|---|
| Source | Stable ID, URL/DOI/PMID, version/commit, retrieval date, checksum where available, access/license, primary/secondary status, correction/retraction status. |
| Entity | Type, exact original identity, canonical IDs, aliases with scope, formulation/strain/chemical attributes and identity-match state. Preserve unresolved identity. |
| Observation | Sample/date, assay/pipeline/database, entity, value/unit, relative/absolute/viable/activity distinction, measured/inferred status, uncertainty and provenance. |
| Experiment | Study/experiment/arm IDs, model/population, exact intervention/comparator, context, compartment, time window, outcome, measurement method and source spans. |
| Assertion | Stable evidence ID, typed subject/relation/object, observed direction, endpoint, experiment/source IDs, evidence state, causal support, applicability, quality flags, counterevidence and review state. |
| Pathway | Ordered assertion IDs, compatibility, prerequisites, opposing sources/sinks, end-to-end support, earliest unsupported bridge, feedback and contradictions. |
| Combination | Exact members/formulations, context/order, observed endpoint, experiment ID, interaction modifier and combination support. Whole-mixture outcomes are not independent effects of every member. |
| ModelRun | Adapter/model/version, input hashes, feature mapping, assumptions, coverage, parameter provenance, seed where relevant, numerical status, outputs, uncertainty and validation domain. |
| Assessment | Candidate/dimension, assessment status, evidence IDs/type, context match, observed direction, separate prediction references, counterevidence and unknowns. |
| Plan | Immutable version, goals, baseline observation IDs, exposures/unknowns, stages, dependencies/conflicts, alternatives, rank basis, conditional reassessment and source IDs. |

Keep distinct enums:

- `evidence_state`: observed, database_curated, predicted, hypothesis, identity_only.
- `assessment_status`: not_assessed, assessed_with_evidence, assessed_no_relevant_evidence, conflicting, not_applicable with reason.
- `observed_direction`: increase, decrease, mixed, null_result, unknown.
- `path_support`: measured_same_experiment, linked_across_experiments, hypothesis; causal mediation is a separate attribute.
- `sequence_support`: tested_sequence, tested_combination, assembled_from_separate_studies, hypothesis.

Relations include consumes, produces, releases, inhibits, competes, transforms, alters_habitat, host_response and associated_with. Production, secretion and host absorption are separate edges. Gene carriage, expression, activity, metabolite concentration, relative abundance and symptoms are separate endpoints. Absence of evidence is not a measured null; an empty concern list is not a safety result.

Preserve stable IDs, raw values and all normalization/quality overlays when acquiring upstream snapshots. Aliases, copied database rows, reanalyses, multiple endpoints and multiple versions of one experiment are not independent replication. Do not silently merge strains, related species, GTDB suffixes, metagenomic species labels or grouped chemicals. Public access does not itself establish redistribution permission.

## Universal assessment vocabulary

Every candidate and combination inherits all applicable dimensions with explicit state. Inheritance initializes unknowns; it does not supply scientific evidence. A finding resolves only its actual identity, setting and endpoint. `not_applicable` requires a reason.

| Original dimensions 1–12 | Original dimensions 13–24 |
|---|---|
| exact_material | redox_and_electron_acceptors |
| baseline_exposure | mucus_and_attachment |
| host_nutrient_need | biofilms |
| delivery | bile_acid_ecology |
| direct_growth_or_inhibition | nitrogen_and_protein |
| collateral_inhibition | host_response |
| primary_degradation | drug_and_compound_metabolism |
| public_breakdown_products | resistance_and_persistence |
| metabolic_cross_feeding | ecological_loss_and_rebound |
| partner_competition | combined_exposures |
| host_nutrient_competition | time_and_location |
| pH_and_buffering | contrary_and_null_evidence |

| Twelve refinements | Distinction to preserve |
|---|---|
| absolute_load_and_viability | Relative/absolute, viable/total/active quantities |
| compartment_exposure_and_host_clearance | Local exposure, uptake, retention and clearance |
| bioaccumulation_and_release | Reversible sequestration versus transformation |
| functional_redundancy_and_backup_paths | Alternative contributors; inferred versus measured function |
| regulatory_and_phenotypic_state | Genetic capacity versus active, conditional phenotype |
| priority_history_and_residency | Residents, incoming organisms, prior exposure and order |
| hysteresis_and_reversibility | Persistence, withdrawal and recovery evidence |
| evolution_and_trait_stability | Time-stamped phenotype and stability |
| resistance_ecology_and_mobile_elements | Selection, attribution, mobility and measured transfer |
| nonbacterial_and_nonmetabolic_interactions | Phages, fungi, protists, contact effects and host immunity |
| sampling_phase_and_baseline_variation | Collection timing, recent context and within-person variation |
| observation_and_identifiability | Observed quantities, assumptions and unmeasured state |

Direct, second-order and third-order effects are pathway views relative to an action. Preserve supported feedback and longer paths. A higher-order interaction, where another component modifies an interaction, also needs a combination/modifier record. Connected edges do not prove the complete path. Include opposing consumers, host sinks, compartment differences and unmeasured bridges. Never multiply edge signs into a personal success probability.

## Model comparison and integration conditions

Integrate an adapter for a defined question with compatible inputs and an evaluation target. This comparison does not require installing every framework. The companion corpus retains the full source/model catalogues.

| System | Role | Entry condition and boundary |
|---|---|---|
| MICOM + AGORA2 | Primary metabolic scenarios | Matched taxa/models and explicit diet assumptions. Steady-state metabolism does not supply immunity, delivery, time-to-response or durable engraftment. |
| COMETS | Dynamic/spatial metabolic scenarios | Environmental, kinetic and spatial parameters beyond a stool profile. Model geometry is not measured gut structure. |
| iDynoMiCS 2 | Individual-cell 2D/3D scenarios, mechanics, adhesion/EPS and reaction–diffusion | Requires specified geometry, rules and calibrated parameters. FASTQ cannot identify these. Calibrated scenarios remain distinct from measured gut biofilm and clinical outcomes. |
| BacArena | Spatial individual-agent FBA comparator | Compatible models and spatial/environmental rules; same identification limits. |
| MDSINE2; BEEM comparators | Longitudinal ecosystem dynamics | Adequate repeated observations, perturbation context and load data where required. One stool sample does not identify personal dynamics. |
| SMETANA; GutCP | Interaction potential / proposed missing exchanges | Dependencies and inferred edges remain hypotheses. |
| MIMOSA2 | Relate metabolic potential to measured chemistry | Independent paired microbiome/metabolome observations; imputed chemistry cannot validate itself. |
| McMLP | Human dietary response comparator | Match preprocessing, intervention domain and endpoint; performance varies by endpoint. |
| Venturelli miRNN | Experimental fiber/community dynamics and combinations | Match species, fibers and measured conditions. Existing 1,193 endpoint records are this same study. |
| MIDASim; MB-GAN | Synthetic distribution and robustness tests | Development-only fitting. Distributional fidelity does not establish an intervention effect or causal transition. |
| Conditional generators/foundation models | Evaluate response heads, representations or desired-state proposals | Audit conditioning, task-specific validation and leakage. A generated desired community still needs evidence of reachability. |
| Whole-body/mgPipe; Nutrition Toolbox | Host exchange and food-to-model inputs | Actual host/intake data. Computational feasibility additions are not observed consumption. |

Primary references: [MICOM](https://doi.org/10.1128/mSystems.00606-19), [2026 intervention analysis](https://doi.org/10.1371/journal.pbio.3003638), [AGORA2](https://doi.org/10.1038/s41587-022-01628-0), [COMETS](https://doi.org/10.1038/s41596-021-00593-3), [iDynoMiCS 2](https://doi.org/10.1371/journal.pcbi.1011303) and [code](https://github.com/kreft/iDynoMiCS-2), [McMLP](https://doi.org/10.1038/s41467-025-56165-6), [miRNN](https://doi.org/10.1038/s41589-026-02272-4), [MIDASim](https://doi.org/10.1186/s40168-024-01822-z).

For MICOM, explicitly verify the selected LP/QP backend: the reviewed open path uses HiGHS + OSQP; GLPK alone cannot perform cooperative-tradeoff QP. Record solver and completion status. Pin code/model versions and separate licenses. Report abundance coverage and relevant organism/function coverage; missing low-abundance decision targets matter. Retain strain-mapping and gap-filling uncertainty. Numerical feasibility is separate from biological validity.

## Conditional sequence logic

1. Establish goals and measured state: dates, resolution, quantities, symptoms, exposures, prior response and unknowns. Compare addition, replacement or withdrawal with actual baseline; unknown consumption is not zero.
2. Retrieve exact identities and favorable, null, adverse and contradictory findings across eligible study types.
3. Assess direct/downstream effects, all applicable dimensions, exact combinations and current coexposures. Separate documented conflicts from hypothetical concerns.
4. Compare options with separate benefit, concern, applicability, uncertainty, reversibility and burden dimensions. Include the current effective approach, no additional intervention and a discriminating measurement. Ranking weights are declared preferences, not clinical probabilities.
5. Propose the next supported action and conditional later stages. Each stage carries action/formulation ID, purpose, prerequisites, evidence/counterevidence, sequence-support status, observable endpoints, unknowns and alternatives.
6. Define continue_if, change_if, stop_if and next_stage_conditions against observed response. Reassessment timing has a source basis or stays unknown. Follow intended and collateral outcomes; record changed coexposures and competing explanations.

Do not hardcode suppression → biofilm removal → reseeding. Do not invent administration, dosage or duration to complete a schema. Keep tested concurrent combinations intact without inventing synergy from separate component findings. Nutritional needs stay in the objective; microbial substrate use alone does not justify host nutrient restriction. Hypothetical adverse paths affect conditional ranking rather than automatically excluding an option. Missing paths cannot establish harmlessness.

Use short-horizon observation-updated planning now. Offline RL, POMDP policies and longitudinal causal estimators need appropriate action/outcome trajectories and identification assumptions. Literature rows and the supplied reports are not patient treatment trajectories. Collect exact exposure identity, adherence, coexposures, decision reasons, prespecified endpoints, adverse effects and assay metadata for later evaluation.

## Priority corrections

Carry the complete source-quality overlay from the companion corpus into imports, retrieval and explanations. Regression cases must include:

- `FIX-GUS-01`–`06`: enzyme activity, gene abundance, blood and fecal endpoints differ. Keep experimental model, preliminary/trend status, null results and lack of extra mixture benefit. Aliases do not create independent studies.
- `FIX-RGN-01`–`03`: aqueous extract, isolated constituent and host/community findings have distinct identities/settings. A target-specific culture result establishes neither broad selectivity nor human response.
- Arginine/digoxin findings must not become a reversed restriction claim; growth and drug metabolism are different endpoints ([PMID 23869020](https://pubmed.ncbi.nlm.nih.gov/23869020/)).
- The Zamani trial does not support the claimed Eggerthella endpoint and has an [Expression of Concern](https://doi.org/10.1111/1756-185X.15299); this is distinct from retraction.
- Collinsella/fiber evidence is an observational genus association, not exact-species feeding evidence ([PMID 29144833](https://pubmed.ncbi.nlm.nih.gov/29144833/)).
- Exact Dorea-species restoration/reduction is not established by the cited inulin/low-FODMAP studies (PMIDs 28213610, 25016597). The Parasutterella physiology source does not establish the proposed exact-species fat-reduction intervention (PMID 30742017).
- Preserve the oil-tea/tea-tree-oil material correction, all retraction flags, unresolved mismatches and disputed values. Related-species and genus/family results do not become exact-target claims. Retracted records remain traceable but cannot support positive action claims.

## Implementation sequence and acceptance

Extend in order: import/provenance → identity/assertion semantics → assessments/contradictions → conditional planning/explanation → scoped model adapters → prospective feedback. Added models must show incremental value against the evidence-only implementation.

1. Imports retain every row/object, stable IDs, exact values and quality/citation fields. Counts and hashes reconcile with manifest.json; repeated ingestion is idempotent. References resolve or remain explicitly flagged historical gaps.
2. Strain/formulation mismatch, genus-to-species transfer and alias collisions cannot create an exact-match claim. Unknown, blank, tested negative and conflicting results remain distinguishable.
3. Retractions/material corrections affect retrieval and explanation. Shared-study rows do not inflate support; mixture outcomes do not become independent component successes.
4. All 36 assessment dimensions have explicit states. Missing collateral assays cannot create selectivity claims. Relative abundance and gene presence do not become absolute growth or measured activity.
5. Cross-study pathways expose context mismatch and unsupported bridges. Opposing sinks and contradictory evidence remain retrievable. Predictions never overwrite observations.
6. Plans include alternatives and response-dependent reassessment. No unsupported fixed sequence, invented administration or fabricated success probability appears. New observations create traceable plan versions.
7. Adapters reject incompatible/missing mandatory inputs and expose coverage, assumptions, failures and out-of-domain status. Synthetic results remain labeled and cannot serve as independent validation.
8. Evaluate independently adjudicated holdouts with participant/experiment/study isolation and appropriate chronology. Compare evidence retrieval, text RAG, graph planning and added prediction layers. Measure citation entailment, applicability, contradictions, harmful omissions, calibration and sequence consistency separately.
9. Compare predictive models with simple baselines and class-specific metrics. Independent real holdouts determine biological prediction performance; synthetic tests alone do not. Prospective outcome improvement remains a separate evaluation from faithful evidence retrieval and coherent plans.

## Upstream data acquisition and refresh

Extracted from the existing atlas and its research files, whose stated retrieval date is **2026-09-24**. This is an acquisition plan, not a new network check. Download data into separate files; ingest them with code. Do not paste source records into implementation instructions or an agent's prompt. Counts below reconcile the old snapshot, not every future release.

### Source downloads and existing snapshot

| Source / old embedded scope | Exact acquisition route and snapshot | Import requirements / interpretation | Reuse status recorded in atlas |
|---|---|---|---|
| **gutMDisorder 3.0** — 10,378 gut intervention–taxon rows plus 832 comparison records | [Literature XLSX](https://bio-computing.hrbmu.edu.cn/gutMDisorder_api/api/resource/download?fileName=gutMDisorder_v3_Literature-based_All.xlsx); [reanalyzed-data XLSX](https://bio-computing.hrbmu.edu.cn/gutMDisorder_api/api/resource/download?fileName=gutMDisorder_v3_Rawdata-based_All.xlsx); [resource page](https://bio-computing.hrbmu.edu.cn/gutMDisorder/resource.dhtml). Version 3.0, retrieved 2026-09-24; endpoints are mutable. | Full files have 35,643 association rows: 13,402 literature + 22,241 reanalyzed. Filter intervention class and gut site: 4,921 literature + 5,457 reanalyzed = 10,378. Keep 873 oral and 21 mixed/unresolved intervention rows separately. Join comparison and sample metadata; retain source row, condition_1, condition_2, direction, host label, PMID/accession, statistical values and origin. Exact source workbook worksheet names and a runnable extraction script were not recovered in this audit; inspect the downloaded workbooks before mapping. | Database bulk reuse terms not established. Attribution does not resolve redistribution rights. |
| **Probio-Ichnos** — 12,962 property rows | [Pinned JSON](https://raw.githubusercontent.com/Mtsif/probio-ichnos/e74b82f433aa026c7db6fb4b899938236270d55d/probioIchnos.json); [repository](https://github.com/Mtsif/probio-ichnos). Commit `e74b82f433aa026c7db6fb4b899938236270d55d`. | Preserve source species, exact strain and PMID. 11,165 raw species–strain pairs; 2,255 PMIDs; six exact duplicates. Decode Yes/No/None as reported positive/reported negative/not reported. Unexpected adhesion `A` remains unrecognized (old PI09914). Assay details and target organisms are absent from this source; do not invent them. Keep deduplication and source assertions separate from independent primary evidence. | Atlas records MIT, copyright 2026 Tsifintaris Margaritis; retain [license](https://github.com/Mtsif/probio-ichnos/blob/main/LICENSE) and confirm its scope for redistributed contents. |
| **Glycobif / carbohydrate compendium** — 1,290 strain–substrate outcomes | [Pinned growth matrix](https://raw.githubusercontent.com/Arzamasov/compendium_manuscript/63b5216d891cdfd7a4a0fc8a3419bc6a92e61fab/data/growth/summary/Bif_30_strains_growth_data_full.txt); [compendium repository](https://github.com/Arzamasov/compendium_manuscript); [Glycobif repository](https://github.com/Arzamasov/glycobif); [paper](https://doi.org/10.1038/s41564-025-02056-x). Compendium commit `63b5216d891cdfd7a4a0fc8a3419bc6a92e61fab`. | 30 strains × 43 substrates; 1,247 tested and 43 not tested. Preserve exact strain/substrate identifiers and source classifications; growth, weak growth, no growth, not tested are distinct. Do not fold the larger genome-predicted pathway collection into measured outcomes. Related accessions: GSE239955 and PRJNA1126848. | Compendium MIT, copyright 2022 Aleksandr Arzamasov; separate Glycobif software GPL-3.0. |
| **Connors–Thompson / Venturelli 2026** — 1,193 processed condition endpoints | [Data archive](https://doi.org/10.5281/zenodo.19210709); [code archive](https://doi.org/10.5281/zenodo.20336398), code v1.0.1; [repository](https://github.com/VenturelliLab/Connors-Thompson_et_al_2026); [paper](https://doi.org/10.1038/s41589-026-02272-4). See file list below. | Seven CSVs; 1,193 file/condition endpoints, 1,185 unique condition labels. Group by source file + condition; choose highest numeric Time. Time is source passage/index, not hours. Keep decimal strings and nulls; exp0 organic acids are unmeasured. Keep initial design separately from measured output. Codes require the source strain map, not guessed species aliases. This dataset is already represented in the atlas: do not import it twice. | Both named Zenodo archives CC BY 4.0. |
| **IPDB** — 52 genome annotation profiles | [IPDB_latest.zip](https://probiogenomics.unipr.it/files/IPDB_latest.zip); [extended resource paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11684986/). Mutable archive retrieved 2026-09-24. | Contains `IPDB_info_V2.xlsx`, 52 GenBank + 52 corresponding nucleotide files and comparative assets. Two formats are not 104 strains; preserve aliases and potentially redundant entries. Workbook annotation counts are predicted capabilities. Keep numeric 0 distinct from source dash. Cell B4 for ADO11 contains literal `time`; preserve raw value and normalize numeric value to null + quality flag. | Separate archive redistribution license not verified. |
| **iProbiotics** — 650 training labels | [Download page](https://bioinfor.imu.edu.cn/iprobiotics/public/download.html); main [positive file](https://bioinfor.imu.edu.cn/iprobiotics/public/static/data/Ptobiotics_239.csv) and [negative file](https://bioinfor.imu.edu.cn/iprobiotics/public/static/data/Non-probiotics_412.csv). Preserve spelling `Ptobiotics_239.csv`. | Actual downloaded files contain 239 positive + 411 negative records. Negative filename/paper say 412, but accessible file has 411 sequence-header lines and no header: do not fabricate record 651. Preserve assembly, BioSample, BioProject or sequence accession as supplied. Subgroup/balanced files overlap and are not additional observations. These are training labels, not treatment evidence. | Download-data redistribution terms not established. |
| **gutMGene 2.0** — 4,860 relations | [Portal](https://bio-computing.hrbmu.edu.cn/gutmgene/); [paper](https://doi.org/10.1093/nar/gkae1002). Three direct CSVs below. Curation cutoff recorded as 2023-10-31. | 2,488 microbe–metabolite + 1,323 microbe–host-gene + 1,049 metabolite–host-gene. Use correct encoding per file; retain original strings and flag replacement characters. Preserve relation type, exact strain, substrate, chemical/gene identifiers, experimental method and source host/sample context. Source `human` and `causally` do not automatically mean a human clinical intervention. | Bulk redistribution terms unverified. |
| **cFMD** — 4,125 samples, 109 datasets, 18,180 MAG index rows | [Repository](https://github.com/SegataLab/cFMD), v1.3.2 commit `7c09082463a88a9c71dc9c38ba0f9786ac7c50f6`. Download pinned TSVs listed below plus README, metadata rules and license. | 3,293 fermented + 832 other samples. Join dataset/sample/run/study/project identifiers; preserve food/process categories and specimen stage. Full FASTQs, MAG FASTAs and functional profiles were not downloaded or embedded; the indexes point to them. Sequence presence is not a viable-organism dose or health-effect claim. | Downloaded repository license recorded as CC BY 4.0; preserve attribution and any individual dataset terms. |
| **NJS16** — curated references, assertions, nodes and 4,483 transport/degradation relationships | [Dryad](https://doi.org/10.5061/dryad.mc1j9); [paper](https://doi.org/10.1038/ncomms15393); original retrieval used [Europe PMC supplementaryFiles](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5467172/supplementaryFiles), **Supplementary Data 2 XLSX**. | Retain reference numbers linking assertions to workbook citations; preserve source row, grouped organism/chemical labels, metabolic activity and direction. Includes 570 organism/host nodes and 244 compound groups. Commas inside grouped synonyms are not delimiters for inventing strains. Broad chemical groupings are not exact stereochemical identities. Workbook sheet mapping/extraction code not recovered here. | Atlas attributes source-derived network under paper CC BY 4.0; retain original Sung et al. attribution and dataset terms. |
| **GutCP** — 293 consensus predictions + 7,348 genome-model exchange rows | [Repository](https://github.com/maslov-group/ML_human_gut), commit `5e5451ef27f4df12b3c9192793f068c1d8c57906`; [paper](https://doi.org/10.1038/s41467-021-21586-6). Download `data/names_ID.txt`, `data/pruned_chia_network.csv`, `cluster/SI_table1.csv`, `GSMM/SI_table2.csv` from pinned raw base below. | SI_table1 is prediction layer; SI_table2 is genome-model capability layer. Keep separate from NJS16 curated assertions and avoid counting overlap as independent evidence. Preserve p-value and simulation prevalence as source simulation measures, not patient probabilities. | Retrieved repository tree had no repository-level license; article openness does not settle software/data reuse scope. |
| **2026 Gibbons model-diet inputs** — two published constraint tables | [Code](https://github.com/Gibbons-Lab/2024_probiotic_engraftment), commit `a494f4e21c8510cf02618af1861f8649d92fe698`; [archive](https://doi.org/10.5281/zenodo.18037976); [paper](https://doi.org/10.1371/journal.pbio.3003638). Download `european_medium.csv` and `high_fiber_medium.csv` from pinned raw base below. | Preserve duplicate `reaction` columns with positional suffixes instead of overwriting them. Flux values are model constraints, not dietary quantities or doses. Some raw trial-A sequence data require investigator access; trial B accession PRJNA755324. | Retain archive/publication terms; no blanket software license inferred. |
| **NIH DSLD ingredient groups** — 6,467 groups and 78,788 aliases | [API docs](https://api.ods.od.nih.gov/dsld/v9/); [guide](https://dsld.od.nih.gov/api-guide). Request `https://api.ods.od.nih.gov/dsld/v9/ingredient-groups/?method=by_letter&term=A&from=0&size=1000`; iterate A–Z and Other, paginating from offsets until reported totals reconcile. | Preserve group ID, name, source categories and aliases. API 9.5.0/May 2026; backend index names reported 2025-09-25. Capture API/index version, retrieval date and counts. Aliases support search; they do not prove chemically/biologically identical ingredients. The source is a US label catalogue, not an efficacy database or global product inventory. | API documentation declared CC0-1.0. |

### File-level details needed by adapters

**gutMGene direct CSVs:**

- [Gut Microbe–Microbial metabolite](https://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/dow/downloadFolder?filepath=allfile/Gut%20Microbe-Microbial%20metabolite.csv): `utf-8-sig`, 2,488 rows.
- [Gut Microbe–Host Gene](https://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/dow/downloadFolder?filepath=allfile/Gut%20Microbe-Host%20Gene.csv): `gb18030`, 1,323 rows.
- [Microbial metabolite–Host Gene](https://bio-computing.hrbmu.edu.cn/gutMGene2.0_api/dow/downloadFolder?filepath=allfile/Microbial%20metabolite-Host%20Gene.csv): `utf-8-sig`, 1,049 rows.

**cFMD pinned raw base:** `https://raw.githubusercontent.com/SegataLab/cFMD/7c09082463a88a9c71dc9c38ba0f9786ac7c50f6/`. Files: `cFMD_metadata.tsv`, `cFMD_datasets.tsv`, `cFMD_mags_list.tsv`, `cFMD_metadata_rules.tsv`, `README.md`, `LICENSE.txt`.

**GutCP pinned raw base:** `https://raw.githubusercontent.com/maslov-group/ML_human_gut/5e5451ef27f4df12b3c9192793f068c1d8c57906/`.

**Gibbons model inputs pinned raw base:** `https://raw.githubusercontent.com/Gibbons-Lab/2024_probiotic_engraftment/a494f4e21c8510cf02618af1861f8649d92fe698/`.

**Venturelli CSV paths:** use [the archived release](https://doi.org/10.5281/zenodo.19210709) to reproduce the snapshot; raw URLs in the old atlas used mutable `main`, with hashes to detect change. Base: `https://raw.githubusercontent.com/VenturelliLab/Connors-Thompson_et_al_2026/main/`.

| Path | Source rows | File/condition endpoints |
|---|---:|---:|
| `data/exp0/exp0_metabolites.csv` | 774 | 387 |
| `data/exp1/exp1_metabolites.csv` | 1,052 | 263 |
| `data/exp2/exp2_metabolites.csv` | 948 | 237 |
| `data/exp3/exp3_metabolites.csv` | 1,128 | 282 |
| `data/exp4/exp4_metabolites_best_reps.csv` | 32 | 8 |
| `data/exp4/exp4_metabolites_new_best.csv` | 32 | 8 |
| `data/exp4/exp4_metabolites_new_worst.csv` | 32 | 8 |

Its [strain mapping workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-026-02272-4/MediaObjects/41589_2026_2272_MOESM4_ESM.xlsx) provides exact identities; the [supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-026-02272-4/MediaObjects/41589_2026_2272_MOESM1_ESM.pdf) documents units. Existing audit verifies abundance as an OD600-based proxy, pH, and acetate/butyrate mM. Lactate unit was not independently verified: preserve source scale without inventing units.

### Custom data and enrichment preservation

The five report-derived inventories and lookup mappings are specific to supplied reports. The supplied SQLite preserves these derived records, including report/page identity, original labels and the recorded PDF hashes. Their underlying reports are private inputs, not upstream public datasets. Reports: SAMPLE1_A01 (236 pages), SAMPLE2_A02 (274), SAMPLE6_A06 (230), SAMPLE4_A04 (213), SAMPLE3_A03 (248). Their 1,309 printed organism rows are report observations, not 1,309 independent intervention studies.

Keep bespoke curated observations, quality corrections, evidence links and ecological chains as project data too: an upstream database download does not recreate that work. The publication metadata index can be refreshed through [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), but retain its PMID/DOI keys, known retraction/concern flags and the curated evidence links before replacing it. A refreshed index does not prove full-text review.

For reproducibility, persist each downloaded file's URL, immutable release/commit where available, retrieval timestamp, SHA-256, encoding, source counts, transformation version and license/access status in machine-readable provenance. Keep original row identifiers and a mapping to any local IDs: simply downloading originals does not recreate locally assigned PI/GMD/GMG/NJS IDs or derived joins unless the adapter implements that mapping.
