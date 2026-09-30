# BUILD_SPEC v4.2 — Gut microbiome profiles for inflammatory skin disease

**Research cutoff:** 7 September 2026. **Artifact:** one complete coding-agent handoff. **Target:** extend the existing live microbiome application without changing unrelated profile results.

## 0. Execute this document

Implement the registry, observation adapters, directional scoring, explicit unavailable states, evidence retrieval, report integration, and acceptance tests below. This document contains the research needed for these changes. Do not require a separate research report, claim ledger, previous build spec, or supplementary handoff document. Study datasets and software assets fetched during implementation are runtime/build inputs, not additional handoff documents.

Inspect the actual repository before editing: profile loader, scored profiler/database, reference bundles, results schema, report renderer, and test runner. Historical output used MetaPhlAn 3; some newer observations used MetaPhlAn 4. Do not assume that the extended lane already supports scoring. Adapt the contracts below to the repository language and architecture. Preserve existing alopecia, Alzheimer’s, chronic urticaria, disease, functional, organism, age, and intervention modules.

**Deliver working research-pattern calculations when their analytical prerequisites exist.** Missing reference data or an unresolved feature produces an explicit reason and partial coverage; it must not be replaced by invented healthy abundances. Independent clinical validation is a separate milestone from implementing this specification.

### 0.1 Product decisions

There is no established, uniquely diagnostic stool signature for the skin diseases reviewed here. There are useful, measurable research patterns. “Unique” would require testing against other relevant diseases, medications, and populations; healthy-versus-case AUC alone does not establish this.

| Condition | v4.2 implementation | Meaning of the result |
|---|---|---|
| Plaque psoriasis | Upgrade the existing family with separately versioned native-shotgun study panels, especially the 2026 SGB4348 signal; preserve the legacy result | Research pattern concordance; a high result does not establish psoriasis |
| Psoriatic arthritis, PsA | Separate nine-KO shotgun panel; shared psoriatic-disease taxonomy context and optional subtype functional observation | Research resemblance; not a validated test distinguishing PsA from skin-only psoriasis |
| Adult atopic dermatitis | Thirteen-slot 16S-derived genus/group panel with explicit shotgun transport and mapping coverage | Exploratory translated pattern; source groups without an exact crosswalk remain unavailable |
| Advanced non-segmental vitiligo | Eleven-slot native-shotgun panel plus a separate four-species 2025 exploratory panel | Cohort-specific resemblance; no diagnostic probability |
| Seborrheic dermatitis | Complete evidence/measurement page and registered unavailable stool-DNA task; preserve culture, genetic-association, and scalp findings as different evidence types | No defensible stool-DNA score at this cutoff; this is not a negative test |
| Hidradenitis suppurativa, HS | One-feature translated stool research panel; display the negative 2026 shotgun comparison prominently | Very limited-feature resemblance; no claim of uniqueness |
| Infant atopic dermatitis | Register a separate, non-default strain-reproduction task | Six-month-old infant strain evidence; never applied to adults or children aged 8–18 |
| Pemphigus foliaceus, PF | Cross-disease counterevidence and research dataset; no new default disease score | Particularly important because it shares the proposed psoriasis SGB4348 marker |

The strongest additional **native-shotgun** candidate identified is non-segmental vitiligo; PsA adds a clinically meaningful, distinct phenotype. Adult atopic dermatitis has substantial human gut evidence, with the exact openly recoverable panel here derived from 16S. HS evidence is considerably weaker. Acne and rosacea do not enter the default scored additions in this release. Existing stronger non-skin families such as CRC, IBD, cirrhosis, and T2D remain intact rather than being duplicated under a skin release.

### 0.2 Naming and migration

- Use `SEBD` for seborrheic dermatitis. Existing `SD` can mean symptomatic dermographism and must not be reassigned.
- Use `ATD` for atopic dermatitis. Existing `AD` belongs to Alzheimer’s and must not be reassigned.
- Distinguish `PSO` plaque psoriasis from `PSA` psoriatic arthritis, including coexisting disease.
- A clinical family, a study panel, and a frozen classifier are different registry entities. Several panels do not inflate the disease count.
- Match existing families using explicit phenotype, body site, task, and population identity. Export the existing literal IDs and digests before migration. Add by set union; never use a guessed historical total as the migration rule.
- Add these panels under the existing psoriasis family when present. Keep its old score, label, and version. New study panels appear as additional readings until a separately validated replacement is deliberately versioned. Do not average them into an existing score silently.
- This handoff specifies no invented file paths in the existing repo. The agent creates normal code/config/test assets within its actual structure and uses the existing build/test entry points.

## 1. Evidence records and exact limits

All features below are **source-cohort associations**, not universal beneficial/harmful assignments. Store one assertion per `(source, participant cohort, contrast, feature, assay)`. A missing numerical effect or q-value is `null`, not zero. The sign is sufficient for the engineering index defined later; it is not a published model coefficient.

### E01 — Deng 2026: psoriasis SGB4348

[Primary study, DOI 10.1186/s12967-026-08013-4](https://doi.org/10.1186/s12967-026-08013-4). Stool shotgun, MetaPhlAn 4/HUMAnN 3; 98 psoriasis participants, 28 cohabiting and 17 unrelated controls. Discovery used 28 cases/28 cohabitants; the same-center test used 70 cases/17 unrelated controls. **Fimenecus sp000432435**, source `s__GGB51647_SGB4348`, was increased. LinDA discovery used FDR<0.1. Single-feature test ROC AUC was 0.84, 95% CI 0.74–0.94. This was not geographical external validation; selecting the best of three models on the test group weakens its independence. Household matching changed diversity findings. Biotin pathway `PWY-5005` did not validate. [Raw project CNP0000322](https://db.cngb.org/search/project/CNP0000322/) is publication-declared; the complete run/phenotype join has not been verified.

**Implementation:** one-feature directional index, visibly labeled as such. The published AUC does not validate our index and must never be attached to the individual result as its accuracy.

### E02 — Deng supplementary audit and model blocker

[Published supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12967-026-08013-4/MediaObjects/12967_2026_8013_MOESM3_ESM.docx) disagree with the main text: three-feature ROC AUC appears as 0.74 versus 0.76; original/SMOTE PR-AUC values 0.83 and 0.96 are reversed between accounts. The methods mix relative abundance and “+1 to read counts.” Repeated five-fold CV is not documented as household-grouped. The clinical workbook contains missing control BMI despite a complete-BMI statement. The [declared repository](https://github.com/qqwxp1987/psoriasis_microbiome) returned 404 during this audit.

Do not resolve these by guessing. Register the original classifier as `blocked_unresolved_preprocessing_and_artifacts`. Recompute its metrics only from a reproducible, locked implementation and source labels. SMOTE must be restricted to inner training folds; original untouched test prevalence must be retained. The separate directional index remains implementable.

### E03 — Sá 2026: multi-disease specificity check

[Primary study, DOI 10.3390/ijms27020838](https://www.mdpi.com/1422-0067/27/2/838): **55 total people**, comprising 24 severe psoriasis, 10 Hurley-III HS, 11 severe PF, and 10 controls. Stool shotgun used MetaPhlAn 4.0.6/CHOCOPhlAnSGB 202212 and HUMAnN 3; approximately one million reads/sample. Psoriasis-versus-control decreases: **Enterocloster**, **Phascolarctobacterium faecium**, **Bacteroides finegoldii**. There was no significant HS-versus-control taxonomic difference and no significant group alpha/beta diversity difference. PF showed increased **SGB4348** and **SGB15101** versus controls/HS. Confounding is substantial: immunosuppression in 18/24 psoriasis, 5/10 HS, 10/11 PF; parasites in 5/24, 2/10, 6/11 respectively. [Versioned analysis deposit](https://doi.org/10.17632/yjz7v8rknx.2), CC BY 4.0, describes processed analyses; downloadable FASTQs/participant joins were not verified.

**Implementation inference:** SGB4348 cannot be advertised as psoriasis-specific. Preserve its PF counterassociation. Do not convert psoriasis-versus-HS or psoriasis-versus-PF features into case-versus-healthy features. Enzyme annotations labeled “virulence” or “resistance” in this paper are not validated resistance/virulence gene calls; do not import those clinical interpretations.

### E04 — Xiao 2024: shared psoriatic disease versus healthy controls

[Primary study, DOI 10.1128/spectrum.01154-23](https://doi.org/10.1128/spectrum.01154-23): 44 skin-only plaque psoriasis, 26 CASPAR-classified PsA, 25 controls, ages 18–65. Stool shotgun used MetaPhlAn 2/HUMAnN 2.8.1. **Pooled psoriatic disease versus controls** showed decreases in **Eubacterium rectale** (synonym Agathobacter rectalis), **Alistipes finegoldii**, and **Alistipes shahii**, with reported q-values <0.01, <0.001, <0.001 respectively. These are not PsA-only q-values. E. rectale decreased in both subgroups; internal discrimination was approximately AUC 0.74–0.77. Same-subject PCR/assembly is not external validation. No significant species difference separated PsA from skin-only psoriasis. `KEGG:00072` decreased in PsA versus skin-only psoriasis (OR 0.22, q=0.044); butanoate-metabolism q=0.157 was not significant. [PRJNA938297](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA938297) has 95 paired-end run records; [declared code](https://github.com/ETaSky/xbiome_psoriasis_ms) returned 404 in this audit.

**Implementation:** one shared three-species psoriatic-disease-versus-healthy panel, linked under both psoriasis and PsA detail pages with its pooled contrast stated. It cannot supply a PsA-specific headline. Keep KEGG00072 as a separately gated, assay-compatible functional observation for a PsA-versus-PsO research comparison; never infer it from the three taxa.

The E04 AUC range concerns **E. rectale abundance alone**, not the new three-species index. Alistipes subgroup-specific estimates remain unknown; do not duplicate pooled q-values into PsA-only or PsO-only assertions.

### E05 — Psoriasis counterevidence and 2026 PsA subtyping

[Yunusbayev, DOI 10.1128/spectrum.01382-24](https://doi.org/10.1128/spectrum.01382-24), online December 2024/2025 issue, studied 53 psoriasis participants and 47 controls with shotgun sequencing. It did not find broad alpha/beta dysbiosis; inflammation/calprotectin confounded apparent profiles. [PRJNA1102742](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1102742) contains 53 psoriasis runs. [PRJNA1061168](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1061168) contains **70 runs from a prior TB case-control project**, not 47 uniformly healthy controls. Import only the paper’s explicitly identified 47 controls. Do not learn “psoriasis = low Shannon.”

[Boix-Amorós 2026, DOI 10.1016/j.ard.2026.01.018](https://doi.org/10.1016/j.ard.2026.01.018) describes 192 axial/peripheral PsA participants and five newly diagnosed IBD cases. It concerns **within-PsA subtype and IBD comorbidity**, not healthy-person PsA screening. Exact assay feature definitions were not recoverable here; no numeric seed is invented from the abstract. Its earlier preprint and final paper count as one cohort.

### E06 — Wang 2023: adult atopic dermatitis

[Primary study, DOI 10.3390/ijms241612856](https://www.mdpi.com/1422-0067/24/16/12856): 104 cases/130 controls, ages 18–68, **16S V4–V5** (515F/907R), QIIME2 2020.11/DADA2/SILVA138. Thirteen genus/group findings are encoded in §3. LEfSe used p<0.05 and LDA>2; only **Clostridium_sensu_stricto_1** was also significant with ANCOM (W=218). No significant Shannon or Firmicutes/Bacteroidetes ratio difference. “Mild/severe” was a median EASI split, not standard severity classification. [PRJNA778863](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA778863) and [metadata/analysis repository](https://github.com/evy-yiweiwang999/AD_gutmicrobiome_case_control), including [metainfo.csv](https://github.com/evy-yiweiwang999/AD_gutmicrobiome_case_control/blob/main/metainfo.csv), are publication-declared; complete raw-to-clinical joining was not verified.

**Implementation:** translated research panel, not a native-shotgun classifier. Do not claim that all butyrate-producing genera decrease: several increase in this cohort. Preserve the stronger ANCOM finding even if its group mapping is unavailable; never substitute the entire genus Clostridium.

### E07 — Additional atopic dermatitis evidence that must not be misencoded

[Spanish adult metagenomic study, DOI 10.1089/derm.2024.0536](https://journals.sagepub.com/doi/10.1089/derm.2024.0536), 2025, includes 38 cases/32 controls. The accessible abstract describes depleted species belonging to several families. Exact species tables and an accession were not verified. A statement about “species belonging to Akkermansiaceae” does not establish an Akkermansia muciniphila coefficient or a family-total direction. Register `source_features_pending`; do not fill its slots from general knowledge.

[Seong 2025, DOI 10.1038/s41522-025-00714-w](https://www.nature.com/articles/s41522-025-00714-w), 31 cases/29 controls at six months, found **no significant species abundance differences** (p>0.1). The signal was *B. longum subsp. infantis* **subclade II** associated with atopic dermatitis versus **subclade I** with controls. [PRJNA979436](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA979436), [PRJEB45443](https://www.ebi.ac.uk/ena/browser/view/PRJEB45443), and [released MAGs](https://doi.org/10.6084/m9.figshare.27367887) support a strain-reproduction branch. The first project has 13 experiments/60 BioSamples and alone does not supply all 60 metagenomes. Supplementary Table 2 contains mapping information. Do not encode B. longum abundance↑, treat strain-specific KOs as community KOs, or transfer this result to an adult.

[Peng 2026, DOI 10.3389/fimmu.2026.1836716](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1836716/full), 53 children aged 2–12 and 16 controls, used 16S and reported AUC 0.941. Feature preselection preceded five-fold CV, with no external validation. CRA071715/OMIX017284 are declared source identifiers; 97%-OTU species labels and this AUC do not justify pediatric shotgun disease scores.

### E08 — Luan 2023: non-segmental vitiligo

[Primary study, DOI 10.1186/s12866-023-03020-7](https://link.springer.com/article/10.1186/s12866-023-03020-7): 25 advanced non-segmental cases/25 matched controls, Xi’an, mean age approximately 31, stool shotgun/MetaPhlAn 3.0.5. The eleven selected directional features are in §3. Methods describe Wilcoxon/BH FDR<0.05 and LEfSe p<0.05/LDA>2; do not fabricate individual q-values. Published 28-feature random-forest nested-CV AUC was 0.786, without external validation; beta diversity p=0.495. [PRJNA1013076](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1013076) is verified as 50 experiments/50 BioSamples. [Published workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12866-023-03020-7/MediaObjects/12866_2023_3020_MOESM1_ESM.xlsx): Tables 3/5 contain taxonomy evidence; Table 8 lists model features. Our eleven-slot directional index is a new engineering implementation, not that random forest.

The paper’s “Staphylococcus thermophiles” versus “Streptococcus thermophilus” discrepancy is quarantined. Do not infer enterotoxigenic *B. fragilis* from species abundance or human NOD activation from a microbial KEGG annotation.

### E09 — Ju 2025: active-spreading vitiligo

[Primary study, DOI 10.3390/ijms26072939](https://www.mdpi.com/1422-0067/26/7/2939): 10 cases versus 20 selected external controls, Kraken2/custom RefSeq and HUMAnN3. Reported directions: *F. prausnitzii*, *F. duncaniae*, *Megamonas funiformis* decreased; *Bifidobacterium bifidum* increased. [Case data PRJNA1222148](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1222148); [control source PRJEB33013](https://www.ebi.ac.uk/ena/browser/view/PRJEB33013). Cases and controls originate in different studies; control selection from 37 people used an unusual combinations procedure. Results/methods disagree on case sex counts. Alpha diversity was not significantly lower. Fungal, plant-virus, and thermophilic-archaeal findings require contamination/mapping review.

**Implementation:** a separate exploratory four-species panel, not a replication of E08’s feature set. In a modern taxonomy that splits historical *F. prausnitzii*, verify distinct source definitions before scoring *F. prausnitzii* and *F. duncaniae* separately; never double-count an umbrella species and one of its children.

### E10 — Hidradenitis suppurativa

[Öğüt 2022, DOI 10.5826/dpc.1204a191](https://dpcj.org/index.php/dpc/article/view/2257): 15 cases/15 matched controls, stool 16S V3–V4. *Fusicatenibacter* decreased (p=0.046); unclassified Clostridiales/Firmicutes findings cannot become independent species. No verified raw accession or external classifier was recovered. Implement the one-genus translated research panel with a prominent specificity/coverage limitation.

[McCarthy 2022, DOI 10.1016/j.jid.2021.05.036](https://doi.org/10.1016/j.jid.2021.05.036) and [Cronin 2023 reanalysis, DOI 10.3389/fmicb.2023.1289374](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1289374/full) are related evidence, not independent replications. Reanalysis includes HS/Crohn’s/controls and repeated samples. [PRJEB43835](https://www.ebi.ac.uk/ena/browser/view/PRJEB43835) and [PRJNA414072](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA414072) require participant and phenotype deduplication. Do not import reported *R. gnavus*/*C. ramosum* directions without confirming the original feature tables. E03 provides an explicit negative shotgun HS-versus-control comparison; keep it visible.

## 2. Seborrheic dermatitis: retain the actual evidence, not a fabricated stool signature

### E11 — Direct human stool culture

[Odintsova/Dyudyun 2019, DOI 10.37321/dermatology.2019.1-2-05](https://repo.dma.dp.ua/5151/), [full primary paper](https://repo.dma.dp.ua/5151/1/05_%D0%9E%D0%B4%D0%B8%D0%BD%D1%86%D0%BE%D0%SAMPLE1%D0%B0.pdf), enrolled 67 cases aged 18–57 and describes 30 controls. It used **stool bacteriology**, not mucosal biopsies or sequencing. Reported low-count frequencies: Bifidobacterium 19/59.4%, Enterococcus 20/62.5%, E. coli 21/65.7%, lactobacilli 29/90.6%, Bacteroides 23/71.9%. These imply a denominator near 32, unexplained against 67 enrolled. Actual per-subject culture values, comparative control distributions, and usable CFU thresholds are absent. H. pylori testing mixed stool antigen, breath testing, and IgG. Detection of Candida, Staphylococcus, or Klebsiella is not a demonstrated relative-abundance enrichment. No reusable sequencing accession was reported.

**Engineering decision:** encode these as culture-source evidence assertions with `relative_abundance_direction=null`. Absolute viable counts cannot be converted into DNA-relative-abundance coefficients. Show available corresponding stool observations as **related measurements whose relationship to this culture finding is unvalidated**, without numeric SEBD concordance. Do not reclassify low E. coli in this paper as high E. coli because of a generic dysbiosis template.

### E12 — Genetic associations

[Zhu 2024, DOI 10.3389/fimmu.2024.1427276](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1427276/full) uses microbiome GWAS and FinnGen disease outcomes. It is not a measured case-control SEBD stool cohort. Nominal forward associations include Tenericutes/Firmicutes/Mollicutes/Senegalimassilia/Victivallis toward increased risk and Butyrivibrio/Eubacterium eligens group/Howardella/Lachnospiraceae NC2004 group/Ruminiclostridium5 toward reduced risk. The printed genus Bonferroni calculation is inconsistent: 0.05/131≈0.0003817, not 0.00382. Broad ranks overlap. These directions concern genetically instrumented exposures, not an observed diagnostic stool profile.

Store `evidence_type=mendelian_randomization`, `observed_case_control_direction=null`, `scoring_enabled=false`. Neither nominal associations nor host SNP effects become taxon weights. Do not recover human genetic risk from host reads in stool FASTQs.

### E13 — Latest clinical series and scalp studies

[Solomakha 2026, DOI 10.22141/ogh.7.2.2026.301](https://oralhealth-journal.com/index.php/journal/article/view/301), published June 27, includes a 40-patient clinical series aged 25–63. It reports dysbiosis/GI comorbidities but supplies no reproducible sequencing feature matrix, accession, or disease classifier. The full methods describe clinical history/examination rather than a defined metagenomic case-control assay. Its “60% dysbiosis” cannot become a disease prior or an index cutoff.

[Skin study, DOI 10.1111/1348-0421.12398](https://doi.org/10.1111/1348-0421.12398) and [2024 topical probiotic study, DOI 10.1038/s41598-024-53016-0](https://www.nature.com/articles/s41598-024-53016-0) measure **skin**, not stool. Malassezia, Staphylococcus/Cutibacterium ratios, scalp lipases, and skin treatment responses cannot seed a gut SEBD classifier. Stool fungal observations remain measured fungal observations; they neither identify scalp overgrowth nor prove contamination automatically.

### 2.1 Required SEBD page and result

Always retain a SEBD row when this release’s skin section is enabled. Status is `insufficient_stool_signature`, score/percentile/probability are null, and text is:

> Gut involvement has been reported, but no adequately specified stool-DNA signature was identified for this implementation. This result cannot assess whether you have seborrheic dermatitis.

The detail page includes three evidence tabs/blocks: stool culture findings, genetic hypotheses, scalp findings. It may link to actual observed Bifidobacterium, Enterococcus, E. coli, Bacteroides, and fungal rows. Do not infer “missing lactobacilli” by collapsing an undefined historical culture group into modern taxa. It also shows the relevant intervention evidence in §8, even though no gut signature score is computed.

The promotion task requires a measured SEBD case-control stool dataset with exact phenotype, body site, feature definitions, assay, clinical metadata, and reproducible analysis. Test against psoriasis/sebopsoriasis, atopic dermatitis, rosacea, and healthy controls. A future culture-based module must have compatible quantitative culture inputs and reference distributions; it is not implemented by converting percentages or literature frequency counts into CFUs.

## 3. Embedded feature manifests

### E14 — Chang 2022: independent psoriasis shotgun panel

[Primary study, DOI 10.1016/j.xjidi.2022.100115](https://doi.org/10.1016/j.xjidi.2022.100115), [author manuscript with supplements](https://escholarship.org/content/qt06t0v09n/qt06t0v09n.pdf): 33 psoriasis/15 controls, San Francisco, MetaPhlAn2 version 2.6.0 and HUMAnN2 version 0.11.1. Seven selected taxon findings are encoded below, including *B. coprocola*↑, which conflicts with the control-enrichment reported by E04. Supplementary Table S1 contains the larger 62-taxon analysis; Table S12 maps participants. Large sparse DESeq2 effects must not become weights. Host sigmoid RNA/PBMC measurements are separate assays. Exploratory PSO1/2/3 clusters are not established clinical subtypes. [PRJNA634145](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA634145) contains 48 paired-end run entries; [SAMN14986030](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN14986030) has explicit `host_disease=PSO` and participant 7305, agreeing with Table S12.

Selected functional findings from Table S3—`P163-PWY`↓, `PWY-5304`↓, `PWY-7200`↑—can be displayed as DNA capacity observations when directly measured. They are not added to the taxon score or treated as independent replication.

### E15 — Liu 2024: native-shotgun PsA functional panel

[Primary study, DOI 10.1177/1759720X241266720](https://doi.org/10.1177/1759720X241266720): 20 treatment-naive CASPAR PsA/10 controls, ages 18–60, stool shotgun. Nine source KOs are encoded below. Selection was LEfSe p<0.05/LDA>2, not external classifier validation. No public accession was recovered. `K07114` is a microbial CaCC-homolog annotation; a conflicting discussion spelling `K07714` is not imported. Neither this KO nor CARD `crp` measures a human chloride channel or serum C-reactive protein.

**Implementation:** a separately scored, directly measured microbial KO panel. No taxa-to-KO imputation. Ten discovery controls do not meet this build’s reference minimum; use an independently compatible frozen functional reference or return `reference_unavailable`. Source feature directions are usable without a pretrained model, but require a declared production measurement definition.

### 3.1 Canonical seed JSON

This is the complete initial directional seed inventory. `source_feature` is an exact source concept, **not an assertion that every production database has that identifier**. Resolve it through §4 and store the resolved native ID/digest. `direction=+1` means higher in the stated case group; `-1` means lower. Unless given, `effect_size`, `q_value`, and disease specificity are unknown. Never fill them from the engineered index. These seed panels do not reconstruct published classifiers.

```json
{
  "release": "4.2",
  "panels": [
    {
      "id": "PSO_DENG2026_SGB4348_V1",
      "family": "plaque_psoriasis",
      "source": "E01",
      "cohort_group": "deng2026",
      "contrast": "psoriasis_vs_controls",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "same_center_tested_candidate",
      "features": [
        {"id":"SGB4348","source_feature":"s__GGB51647_SGB4348","rank":"sgb","direction":1}
      ]
    },
    {
      "id": "PSO_CHANG2022_TAX_V1",
      "family": "plaque_psoriasis",
      "source": "E14",
      "cohort_group": "chang2022",
      "contrast": "psoriasis_vs_healthy",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "single_cohort_selected_subset",
      "features": [
        {"id":"megasphaera_unclassified","source_feature":"Megasphaera_unclassified","rank":"source_group","direction":1,"q_value":0.004186369},
        {"id":"dialister_invisus","source_feature":"Dialister_invisus","rank":"species","direction":1,"q_value":0.001738535},
        {"id":"bacteroides_vulgatus","source_feature":"Bacteroides_vulgatus","rank":"species","direction":1,"q_value":0.01176506},
        {"id":"phascolarctobacterium_succinatutens","source_feature":"Phascolarctobacterium_succinatutens","rank":"species","direction":-1,"q_value":0.040279303},
        {"id":"bifidobacterium_pseudocatenulatum","source_feature":"Bifidobacterium_pseudocatenulatum","rank":"species","direction":-1,"q_value":0.002393849},
        {"id":"coprococcus_art55_1","source_feature":"Coprococcus_sp_ART55_1","rank":"source_species","direction":-1,"q_value":0.000961438},
        {"id":"bacteroides_coprocola","source_feature":"Bacteroides_coprocola","rank":"species","direction":1,"q_value":1.33958e-10}
      ]
    },
    {
      "id": "PSO_SA2026_TAX_V1",
      "family": "plaque_psoriasis",
      "source": "E03",
      "cohort_group": "sa2026",
      "contrast": "severe_psoriasis_vs_healthy",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "small_confounded_cohort",
      "features": [
        {"id":"enterocloster","source_feature":"Enterocloster","rank":"genus","direction":-1},
        {"id":"phascolarctobacterium_faecium","source_feature":"Phascolarctobacterium_faecium","rank":"species","direction":-1},
        {"id":"bacteroides_finegoldii","source_feature":"Bacteroides_finegoldii","rank":"species","direction":-1}
      ]
    },
    {
      "id": "PSORIATIC_SHARED_XIAO2024_TAX_V1",
      "family": "psoriatic_disease_shared",
      "source": "E04",
      "cohort_group": "xiao2024",
      "contrast": "pooled_plaque_psoriasis_and_psoriatic_arthritis_vs_healthy",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "shared_psoriasis_signal",
      "features": [
        {"id":"eubacterium_rectale","source_feature":"Eubacterium_rectale","rank":"species","direction":-1},
        {"id":"alistipes_finegoldii","source_feature":"Alistipes_finegoldii","rank":"species","direction":-1},
        {"id":"alistipes_shahii","source_feature":"Alistipes_shahii","rank":"species","direction":-1}
      ]
    },
    {
      "id": "PSA_LIU2024_KO_V1",
      "family": "psoriatic_arthritis",
      "source": "E15",
      "cohort_group": "liu2024",
      "contrast": "treatment_naive_psoriatic_arthritis_vs_healthy",
      "measurement": "shotgun_ko_fraction",
      "evidence": "single_cohort_nominal",
      "features": [
        {"id":"K02004","source_feature":"K02004","rank":"ko","direction":1},
        {"id":"K01190","source_feature":"K01190","rank":"ko","direction":1},
        {"id":"K05349","source_feature":"K05349","rank":"ko","direction":1},
        {"id":"K01897","source_feature":"K01897","rank":"ko","direction":1},
        {"id":"K01187","source_feature":"K01187","rank":"ko","direction":1},
        {"id":"K07114","source_feature":"K07114","rank":"ko","direction":1},
        {"id":"K06131","source_feature":"K06131","rank":"ko","direction":-1},
        {"id":"K07729","source_feature":"K07729","rank":"ko","direction":-1},
        {"id":"K01005","source_feature":"K01005","rank":"ko","direction":-1}
      ]
    },
    {
      "id": "ATD_ADULT_WANG2023_TRANSLATED_V1",
      "family": "adult_atopic_dermatitis",
      "source": "E06",
      "cohort_group": "wang2023_atopic",
      "contrast": "adult_atopic_dermatitis_vs_healthy",
      "measurement": "shotgun_genus_fraction_translated_from_16s",
      "evidence": "assay_transported_exploratory",
      "features": [
        {"id":"blautia","source_feature":"Blautia","rank":"genus","direction":1},
        {"id":"butyricicoccus","source_feature":"Butyricicoccus","rank":"genus","direction":1},
        {"id":"lachnoclostridium","source_feature":"Lachnoclostridium","rank":"genus","direction":1},
        {"id":"eubacterium_hallii_group","source_feature":"Eubacterium_hallii_group","rank":"source_group","direction":1},
        {"id":"erysipelatoclostridium","source_feature":"Erysipelatoclostridium","rank":"genus","direction":1},
        {"id":"megasphaera","source_feature":"Megasphaera","rank":"genus","direction":1},
        {"id":"oscillibacter","source_feature":"Oscillibacter","rank":"genus","direction":1},
        {"id":"flavonifractor","source_feature":"Flavonifractor","rank":"genus","direction":1},
        {"id":"unclassified_oscillospiraceae","source_feature":"unclassified_Oscillospiraceae","rank":"source_group","direction":1},
        {"id":"romboutsia","source_feature":"Romboutsia","rank":"genus","direction":-1},
        {"id":"clostridium_sensu_stricto_1","source_feature":"Clostridium_sensu_stricto_1","rank":"source_group","direction":-1,"ancom_W":218},
        {"id":"unclassified_butyricicoccaceae","source_feature":"unclassified_Butyricicoccaceae","rank":"source_group","direction":-1},
        {"id":"unclassified_erysipelotrichaceae","source_feature":"unclassified_Erysipelotrichaceae","rank":"source_group","direction":-1}
      ]
    },
    {
      "id": "VIT_NONSEG_LUAN2023_TAX_V1",
      "family": "nonsegmental_vitiligo",
      "source": "E08",
      "cohort_group": "luan2023",
      "contrast": "advanced_nonsegmental_vitiligo_vs_healthy",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "single_cohort_selected_subset",
      "features": [
        {"id":"bacteroides_fragilis","source_feature":"Bacteroides_fragilis","rank":"species","direction":1},
        {"id":"massilioclostridium_coli","source_feature":"Massilioclostridium_coli","rank":"species","direction":1},
        {"id":"lachnospiraceae_bx3","source_feature":"Lachnospiraceae_bacterium_BX3","rank":"source_species","direction":1},
        {"id":"tm7_oral_348","source_feature":"TM7_phylum_sp_oral_taxon_348","rank":"source_species","direction":1},
        {"id":"prevotella_copri_clade_c","source_feature":"Prevotella_copri_clade_C","rank":"source_clade","direction":-1},
        {"id":"dorea_longicatena","source_feature":"Dorea_longicatena","rank":"species","direction":-1},
        {"id":"allisonella_histaminiformans","source_feature":"Allisonella_histaminiformans","rank":"species","direction":-1},
        {"id":"coprobacter_secundus","source_feature":"Coprobacter_secundus","rank":"species","direction":-1},
        {"id":"coprococcus_comes","source_feature":"Coprococcus_comes","rank":"species","direction":-1},
        {"id":"sellimonas_intestinalis","source_feature":"Sellimonas_intestinalis","rank":"species","direction":-1},
        {"id":"bacteroides_bouchesdurhonensis","source_feature":"Bacteroides_bouchesdurhonensis","rank":"species","direction":-1}
      ]
    },
    {
      "id": "VIT_ACTIVE_JU2025_TAX_V1",
      "family": "nonsegmental_vitiligo",
      "source": "E09",
      "cohort_group": "ju2025",
      "contrast": "active_spreading_vitiligo_vs_external_controls",
      "measurement": "shotgun_taxonomic_fraction",
      "evidence": "cross_study_control_confounding",
      "features": [
        {"id":"faecalibacterium_prausnitzii","source_feature":"Faecalibacterium_prausnitzii","rank":"species","direction":-1},
        {"id":"faecalibacterium_duncaniae","source_feature":"Faecalibacterium_duncaniae","rank":"species","direction":-1},
        {"id":"megamonas_funiformis","source_feature":"Megamonas_funiformis","rank":"species","direction":-1},
        {"id":"bifidobacterium_bifidum","source_feature":"Bifidobacterium_bifidum","rank":"species","direction":1}
      ]
    },
    {
      "id": "HS_OGUT2022_TRANSLATED_V1",
      "family": "hidradenitis_suppurativa",
      "source": "E10",
      "cohort_group": "ogut2022",
      "contrast": "hidradenitis_suppurativa_vs_healthy",
      "measurement": "shotgun_genus_fraction_translated_from_16s",
      "evidence": "single_feature_nominal",
      "features": [
        {"id":"fusicatenibacter","source_feature":"Fusicatenibacter","rank":"genus","direction":-1,"p_value":0.046}
      ]
    }
  ],
  "unavailable_tasks": [
    {"id":"SEBD_STOOL_RESEARCH_V1","family":"seborrheic_dermatitis","state":"insufficient_stool_signature","sources":["E11","E12","E13"]},
    {"id":"ATD_INFANT_SEONG2025_STRAIN_V1","family":"infant_atopic_dermatitis","state":"strain_bundle_pending","sources":["E07"]},
    {"id":"PSO_DENG2026_ORIGINAL_MODEL","family":"plaque_psoriasis","state":"blocked_unresolved_preprocessing_and_artifacts","sources":["E01","E02"]}
  ]
}
```

Nine initial numeric panel tasks contain **52 feature slots**, including the group/clade slots that can be unavailable. This is a compiler check on this release inventory, not a count of unique diseases, unique organisms, independent studies, or the whole application’s profiles. Canonicalized feature membership cannot be counted twice within a panel.

### 3.2 Non-scored functional and mechanistic observations

| Source | Feature / source direction | Exact permitted interpretation |
|---|---|---|
| E01 | `PWY-5005`↑ in discovery, failed validation | Measured biotin-biosynthesis DNA potential; zero disease-score weight |
| E14 | `P163-PWY`↓; `PWY-5304`↓; `PWY-7200`↑ | Three separately reported pathway-capacity observations; no inference of stool/serum metabolite levels |
| E04 | `KEGG:00500`, `02030`, `02040`↓ in all psoriasis versus HC | Source-associated functional potential; not skin cytokines or clinical motility |
| E04 | `KEGG:00072`↓ in PsA versus PsO | Subtype-contrast capacity observation; not a healthy-screening signal |
| E15 | CAZy `GH43`↑, `CE4`↓, `CBM50`↓ | Source-specific gene-family potential; not independent replication of KOs from the same samples |
| E01 | vBin_422 / three published contig names | Annotation-only lead; no sequence score without actual FASTA and a validated mapping/reference bundle |

No current numeric skin panel uses Shannon, the Firmicutes/Bacteroidetes ratio, inferred “leaky gut,” inflammation, microbial age, calprotectin, or generic SCFA-producing taxonomy as substitute disease features. Existing measured functional and ecological panels still render independently.

## 4. Observation, taxonomy, assay, and population contracts

### 4.1 Required observation envelope

Every accepted observation has: `sample_id`, `participant_id`, collection time, `body_site`, `assay`, raw/QC provenance, profiler name/version, database name/version/digest, source/native feature ID, mapping version, rank, abundance definition, units, denominator, detection state, value or censoring interval, and a measurement-quality flag. `body_site=stool` is mandatory for the numeric panels above. Observations from skin, scalp, blood, diet records, cultures, metabolomics, and host RNA remain typed external observations.

Accepted numeric modes:

- `shotgun_taxonomic_fraction`: production taxonomic abundance expressed as a fraction of the declared profiled community. Percentages may be divided by 100 exactly once. Do not reclose the selected skin taxa to 100%.
- `shotgun_genus_fraction_translated_from_16s`: sum a frozen, non-overlapping set of native species/SGBs belonging to the exact mapped source genus/group. Query and reference both use this same shotgun computation. Mark `assay_transport=unvalidated_16s_to_shotgun`; it never inherits a 16S model’s sensitivity or specificity.
- `shotgun_ko_fraction`: directly annotated microbial KO abundances from a pinned protein/gene pipeline, after a declared length-normalized abundance step. Normalize each unstratified KO’s nonnegative abundance by the sum across **all measured microbial KOs in that same pinned annotation space**, not just the nine selected KOs. Do not count both stratified and unstratified totals. If the denominator is zero or the upstream quantity is not comparable, all KO slots are unavailable. Version the mapping of multi-KO assignments. Query and reference use identical processing. Source-to-production normalization differences remain an explicit transport limitation.

MetaCyc pathway values, RPKs, copies/cell, absolute CFUs, Ct values, and plasma concentrations do not enter the fraction kernel without their own defined compatible reference adapter. Display them separately. An organism’s presence is not proof of a gene, strain, pathway, metabolite, or treatment response.

### 4.2 Taxonomy resolver

The resolver returns `exact_native`, `verified_synonym`, `verified_aggregate`, or `unresolved`. Store the source definition and actual native member IDs, not only a display name. Verify aliases against the pinned profiler/database taxonomy; downloading a modern label from a search engine is not a crosswalk.

Required cases:

1. SGB4348 requires the **same biological SGB definition** as E01/E03. Display Fimenecus sp000432435 only when verified. A database without that SGB returns unavailable; never replace it with all Firmicutes, an arbitrary unclassified taxon, or a same-number SGB from another namespace.
2. *Eubacterium rectale* / *Agathobacter rectalis* and *Bacteroides vulgatus* / *Phocaeicola vulgatus* are candidate synonym mappings requiring database verification. *Bacteroides finegoldii* is distinct from *Alistipes finegoldii*.
3. `Megasphaera_unclassified` is the source unclassified bucket, not the genus total. Do not substitute *M. elsdenii*.
4. `Prevotella_copri_clade_C` is not all *P. copri*, nor an unverified Segatella replacement. If clade resolution is absent, that slot is unavailable.
5. `TM7_phylum_sp_oral_taxon_348`, `Lachnospiraceae_bacterium_BX3`, and `Coprococcus_sp_ART55_1` require exact source concept mappings; they are not their enclosing phylum/family/genus.
6. Unclassified Oscillospiraceae/Butyricicoccaceae/Erysipelotrichaceae do not equal entire families. SILVA138 group membership must be resolved before a shotgun aggregation can be used. Preserve all 13 ATD slots in the coverage denominator when only eight genera are mapped.
7. Before using the Ju panel, prove that historical *F. prausnitzii* and *F. duncaniae* inputs are disjoint in the source-to-production mapping. If this cannot be established, quarantine the ambiguous pair; do not silently sum or replace either.
8. Maintain spelling discrepancies as source flags; the Luan thermophilus/thermophiles conflict remains excluded. No fuzzy matching across genera.

An ambiguous mapping lowers coverage. It does not imply low abundance or biological absence. Database changes require a new mapping digest, reference bundle, and panel result version; old results remain reproducible.

### 4.3 Population and exposure eligibility

Default participant-facing numeric skin panels support adults **18–65**, intersected with each reference bundle’s verified age range and, for E15, ages **18–60**. This is an engineering release restriction, not a claim that the sources validate every adult. Missing age or an unsupported age produces `population_unsupported` for numbers, while measurements and evidence remain visible. Ages 8–17 do not get adult disease numbers; age 18 alone does not establish comparability with an adult clinical cohort.

Require metadata fields with explicit unknown states: sex, country/region, antibiotic/antifungal/probiotic exposure and dates, systemic immune therapy/biologic class, topical therapy, diet pattern, BMI, stool consistency, smoking, IBD, other skin diagnoses, and known parasitic/GI infection. Unknown is not “none.” Preserve phenotype source: clinician-confirmed, self-reported, or absent/unknown.

Exposure compatibility follows the frozen reference’s declared inclusion rules. Apply existing stricter project gates. Recent antibiotics within 90 days are a conservative default exclusion for participant-facing numeric panels unless a separately validated exposure-specific reference is supplied; store the reason. Do not assume a patient receiving immunotherapy matches a treatment-naive reference. Allow a research-only confounded-domain output when explicitly requested by the research application, labeled separately; do not silently bypass participant gates.

QC requires the existing validated stool pipeline to pass and at least **500,000 usable nonhost read pairs** as an initial engineering floor for ordinary panels. A higher profiler/feature-specific validated depth requirement overrides this floor. Record whether counts refer to reads or pairs. Passing this floor does not establish an organism’s detection limit or sufficient strain depth. No automatic threshold is specified for infant strain reconstruction; that task stays pending until an actual marker/coverage validation bundle exists.

The five current personal samples are regression fixtures, not a reference cohort or disease-validation set. One participant's known Long-COVID status must not tune skin-profile weights; visibly healthy younger donors must not be labeled as microbiome-negative for every skin disease.

### 4.4 Detection and missingness

Distinguish `quantified`, `detected_below_quantification`, `not_detected_with_limit`, `not_assayed`, `mapping_unresolved`, and `qc_unusable`.

- Quantified abundance: use its validated interval if available; otherwise a point interval with `measurement_uncertainty_unavailable=true`.
- Nondetection with a validated reporting limit L: use an interval [0,L]. If the pipeline explicitly defines zero as its nondetection estimator, point=0 with this interval; otherwise no point index for that slot. The interval may still widen full-panel bounds.
- Nondetection without a validated limit is not an observed biological zero; omit the slot from the point calculation and assign the full [0,100] contribution interval.
- Never use a missing table row as a depleted disease feature unless the profiler’s declared complete-output and detection contract justifies the interpretation.
- “Not detected” is not “should be present,” and a low reference position is not a taxonomic deficiency or eradication target.

## 5. Exact directional index and executable reference behavior

### 5.1 Interpretation and reference fitting

The output is a **0–100 directional pattern index**. Fifty means approximately centered on the selected reference after transformation. It is neither 50% disease probability nor an empirically calibrated “normal” boundary. Higher values indicate movement in the selected study’s reported direction; lower values indicate the opposite. Do not append `%`, set “disease positive” thresholds, or reuse published classifier AUCs for this index.

Fit reference statistics offline, using at least **20 distinct eligible control participants per feature** as an engineering minimum. Select one prespecified baseline sample per participant; deduplicate technical replicates and related archive submissions. The fitting kernel enforces unique participant IDs; phenotype/population/assay compatibility is checked before it runs. Existing compatible project controls may support the directional index without waiting for all source-cohort data. This does not validate disease discrimination.

For abundance fraction a, define:

`x = log10(a + 0.000001)`

`mu_j = median(control x_j)`

`scale_j = max(1.4826 × median(abs(control x_j − mu_j)), 0.25)`

`z_j = clip((x_j − mu_j)/scale_j, −6, +6)`

`contribution_j = 50 + 50 × tanh(direction_j × z_j / 2)`

The constants are engineering choices, not learned disease coefficients or biological cutoffs. Controls use the same detection estimator and denominator as queries. Unknown/mapping-failed values cannot be encoded as zeros during fitting. Record detection prevalence, usable control N, estimator/censoring policy, cohort, eligibility rules, native feature definitions, and reference digest.

The panel index is the equal-weight mean over usable point contributions. All declared slots remain in the coverage denominator. One measured feature can produce a limited-feature index, explicitly labeled `1/K`. No point features gives `index=null`.

For K slots and per-slot contribution bounds `[l_j,u_j]`, unavailable slots use `[0,100]`. Full-panel bounds are `sum(l_j)/K` to `sum(u_j)/K`. Valid interval-only measurements may tighten bounds without contributing a point. These are **coverage/censoring bounds, not confidence intervals**. A point interval is not certainty about disease biology or sampling error.

The available-feature index and possible full-panel range have **different scopes when the point mask is incomplete**. An interval-only feature can make the full-panel range exclude the available-feature index without a numerical error. Store `index_scope=available_point_feature_mask` and `bounds_scope=full_declared_panel`; never present full-panel bounds as a confidence interval around a partial-mask index. Render separate labeled tracks when coverage is incomplete.

Do not average study panels into a new family-level number in this release. Display them together with their source contrasts, distinct evidence tiers, and overlap/conflict notes. Shared taxonomy plus functions from the same participants are one evidence group. The same feature appearing in several panels is not several independent observations. Existing family fusion remains versioned and unchanged.

Optional percentiles require a separate, held-out, identically processed reference score distribution on the **same feature mask**. Compute weighted midrank `100 × sum[w_i × (I(S_i<S)+0.5 I(S_i=S))]/sum(w_i)`, equal participant weights within cohort and equal total cohort weights. Name the comparator (`held_out_controls` or `held_out_cases`). In-sample ranks must be labeled as such and not marketed as calibration. Percentile is null when mask-specific calibration is absent.

### 5.2 Dependency-free executable kernel

Port this behavior to the repository language or use it directly. The wrapper must enforce §4 and frozen-bundle integrity before calling it. `reference_rows` are prequalified controls with participant IDs. `observations` are accepted fraction observations only; a null point with finite limits is allowed. All IDs inside one panel are unique. The code is intentionally independent of article-specific parsing.

```python
import math
from statistics import median

EPS = 1e-6
SCALE_FLOOR = 0.25
Z_CLIP = 6.0

def fraction(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ValueError("abundance must be numeric")
    if not math.isfinite(v) or not 0 <= v <= 1:
        raise ValueError("expected finite fraction in [0,1]")
    return float(v)

def fit_reference(rows, feature_ids):
    ids = [r["participant_id"] for r in rows]
    if any(not isinstance(x, str) or not x for x in ids):
        raise ValueError("participant IDs required")
    if len(ids) != len(set(ids)) or len(ids) < 20:
        raise ValueError("at least 20 distinct eligible controls required")
    if not feature_ids or len(feature_ids) != len(set(feature_ids)):
        raise ValueError("unique nonempty feature list required")
    result = {}
    for fid in feature_ids:
        values = [fraction(r["features"][fid]) for r in rows
                  if r["features"].get(fid) is not None]
        if len(values) < 20:
            continue
        xs = [math.log10(v + EPS) for v in values]
        mu = median(xs)
        mad = median(abs(v-mu) for v in xs)
        result[fid] = {"mu": mu, "scale": max(1.4826*mad, SCALE_FLOOR),
                       "n_controls": len(values)}
    return result

def contribution(a, direction, ref):
    if type(direction) is not int or direction not in (-1, 1):
        raise ValueError("direction must be integer -1 or +1")
    mu, scale = ref["mu"], ref["scale"]
    if (isinstance(mu, bool) or isinstance(scale, bool)
            or not math.isfinite(mu) or not math.isfinite(scale)
            or scale < SCALE_FLOOR
            or type(ref["n_controls"]) is not int or ref["n_controls"] < 20):
        raise ValueError("invalid frozen feature reference")
    z = (math.log10(fraction(a)+EPS)-mu)/scale
    z = max(-Z_CLIP, min(Z_CLIP, z))
    return 50.0+50.0*math.tanh(direction*z/2.0)

def evaluate(features, observations, reference):
    fids = [f["id"] for f in features]
    if not fids or len(fids) != len(set(fids)):
        raise ValueError("panel needs unique feature IDs")
    if any(type(f["direction"]) is not int or f["direction"] not in (-1, 1)
           for f in features):
        raise ValueError("invalid direction")
    points, bounds, details = {}, {}, {}
    for f in features:
        fid, d = f["id"], f["direction"]
        if fid not in observations or fid not in reference:
            bounds[fid] = (0.0, 100.0)
            details[fid] = {"state": "unavailable"}
            continue
        o, r = observations[fid], reference[fid]
        lo, hi = fraction(o["lower"]), fraction(o["upper"])
        if lo > hi:
            raise ValueError("reversed interval")
        ends = [contribution(v, d, r) for v in (lo, hi)]
        bounds[fid] = (min(ends), max(ends))
        if o.get("point") is not None:
            point = fraction(o["point"])
            if not lo <= point <= hi:
                raise ValueError("point outside interval")
            points[fid] = contribution(point, d, r)
        details[fid] = {"state": "point" if fid in points else "interval_only",
                        "point": points.get(fid), "bounds": bounds[fid]}
    k, m = len(fids), len(points)
    return {"index": sum(points.values())/m if m else None,
            "usable_count": m, "panel_count": k, "coverage": m/k,
            "interval_coverage": sum(x["state"] != "unavailable"
                                     for x in details.values())/k,
            "full_panel_bounds": [sum(v[0] for v in bounds.values())/k,
                                  sum(v[1] for v in bounds.values())/k],
            "mask": [fid for fid in fids if fid in points],
            "details": details}

def weighted_midrank(score, scores, weights):
    if not scores or len(scores) != len(weights):
        raise ValueError("nonempty aligned scores/weights required")
    if not math.isfinite(score) or any(not math.isfinite(v) for v in scores):
        raise ValueError("invalid score")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("positive finite weights required")
    return 100*sum(w*(1 if s < score else 0.5 if s == score else 0)
                   for s, w in zip(scores, weights))/sum(weights)

def self_test():
    rows = [{"participant_id": f"synthetic-{i}", "features": {"a": .01}}
            for i in range(20)]
    r = fit_reference(rows, ["a"])
    assert r["a"]["scale"] == SCALE_FLOOR
    f = [{"id":"a", "direction":1}, {"id":"b", "direction":-1}]
    o = {"a":{"point":.01, "lower":.01, "upper":.01}}
    x = evaluate(f, o, r)
    assert x["index"] == 50 and x["coverage"] == .5
    assert x["full_panel_bounds"] == [25,75]
    assert evaluate(f, {}, r)["index"] is None
    c = evaluate(f, {"a":{"point":None,"lower":0,"upper":.001}}, r)
    assert c["index"] is None and c["interval_coverage"] == .5
    up = contribution(.1,1,r["a"])
    down = contribution(.1,-1,r["a"])
    assert up > 50 > down and abs(up+down-100) < 1e-10
    assert weighted_midrank(50,[20,50,50,80],[1,1,1,1]) == 50
    assert x == evaluate(f, o, r)
    for bad in (True,-.1,1.1,float("nan"),float("inf"),"0.1"):
        try:
            fraction(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid abundance accepted")
    try:
        fit_reference(rows+[rows[0]], ["a"])
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate participant accepted")
    print("v4.2 numerical contract passed")

if __name__ == "__main__":
    self_test()
```

Synthetic rows above test arithmetic only. They must never enter a patient reference bundle. Reference fitting and inference run as separate stages; query samples, test cases, and other samples in the inference batch cannot update medians, scales, feature selection, imputation, calibration, or weights.

### 5.3 Frozen reference/model bundle

Store `bundle_id`, digest, creation date, participant/cohort manifests, body site, age/population eligibility, exposure policy, profiler/container/database digests, feature-map digest, abundance denominator and units, QC/detection settings, fitted per-feature statistics, source versus independent-reference scope, calibration subjects, feature masks, and validation report embedded in the bundle metadata. This is a runtime data object, not a supplemental handoff.

Never make “20 healthy controls” a global validation claim. If E03, E10, E14, or E15’s small control group cannot meet the feature-level minimum, the panel can use existing compatible independent controls or remain explicitly unavailable. Combining controls across geography/assays without a frozen justified reference policy is prohibited. No holdout participant or relative may leak into feature discovery or fitting.

Retraining a classifier requires participant/household grouping; leave-cohort-out outer validation where multiple cohorts exist; all feature selection, transforms, imputation, scaling, hyperparameters, and SMOTE inside inner training folds; a serialized training transform; an untouched final test; and explicit clinical-variable-only and microbiome-only baselines. Prefer a transparent regularized logistic baseline before attempting a larger model. These requirements do not claim that a deployable original classifier was recovered.

## 6. Data acquisition and reproducible joins

### 6.1 Dataset inventory and priorities

| Dataset | Planned use | Verified access state / constraint |
|---|---|---|
| PRJNA634145, E14 | Independent psoriasis training/reprocessing | 48 paired-end run entries; explicit phenotype join example verified |
| PRJNA938297, E04 | Psoriasis/PsA/healthy reprocessing and discrimination audit | 95 paired-end entries; already trimmed/host-depleted; healthy label example verified, full subtype join required |
| CNP0000322, E01 | Reproduce SGB4348 and household effects | Public landing/review page reachable; run/subject split join not verified |
| PRJNA1102742 + exact healthy subset of PRJNA1061168, E05 | Independent psoriasis direction/specificity check | 53 psoriasis runs; second project has 70 mixed-phenotype runs; exact 47-control whitelist required |
| PRJNA1013076, E08 | Vitiligo reprocessing | 50 experiments/50 BioSamples verified; per-participant phenotype join must be checked |
| PRJNA1222148 + selected PRJEB33013 controls, E09 | Separate exploratory vitiligo reproduction | Publication-declared; independent-study batch effects cannot be removed by simply combining tables |
| PRJNA778863, E06 | Reproduce adult ATD source 16S result | Publication-declared; does not validate the shotgun transport by itself |
| PRJEB43835 + PRJNA414072, E10 | HS versus other inflammatory conditions | Publication-declared; distinguish original cohorts, reanalysis, and repeat samples |
| [Mendeley yjz7v8rknx version 2](https://data.mendeley.com/datasets/yjz7v8rknx/2), E03 | Psoriasis/HS/PF analysis and counterevidence | Dataset metadata verified, CC BY4.0; actual file inventory and raw-read availability not verified |
| PRJNA979436 + PRJEB45443 + Seong MAG deposit, E07 | Optional infant strain reproduction | Multi-project assembly/join necessary; not an adult reference |
| E15 PsA KO cohort | Independent biological source for nine KOs | No public accession recovered; source-cohort reproduction requires author data |

Do not download all human FASTQs indiscriminately before inspecting sizes, licenses, source labels, and the existing project’s data policy. Metadata inventory and code implementation can proceed immediately. Use the project’s configured data budget for bulk retrieval. No author messaging is authorized by this spec; report a data-request blocker if required.

### 6.2 Exact ENA query and metadata procedure

For each PRJNA/PRJEB project use ENA’s read-run report with `run_accession,sample_accession,sample_alias,library_strategy,library_source,library_layout,scientific_name,fastq_ftp,fastq_md5,fastq_bytes`. Follow project-specific redirects and validate header names. An example [psoriasis run manifest](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA634145&result=read_run&fields=run_accession,sample_accession,sample_alias,fastq_ftp&format=tsv) is directly accessible. Fetch sample XML from `https://www.ebi.ac.uk/ena/browser/api/xml/{sample_accession}` and parse its named `SAMPLE_ATTRIBUTES`.

Verified join examples:

- Chang: [SAMN14986030](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN14986030), sample alias/participant 7305, `host_disease=PSO`. Match Table S12’s PID and metagenomic availability.
- Xiao: [SAMN33451581](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN33451581), alias 9455 but `host_subject_id=HC-36`, `host_disease=healthy`, `host_phenotype=healthy`. Numeric sample aliases alone do not define phenotype.

Require one-to-one run→sample→participant links, explicit case/control/subtype labels, baseline/follow-up status, and consistency with source counts. Enumerate discrepant/missing entries and exclude them from label-dependent training until resolved. A run total matching the paper does not prove a correct join. Preserve biological replicates; technical read splits may be concatenated only after sample identity and library compatibility checks.

For E05, obtain the exact healthy whitelist from the [working analysis repository](https://github.com/Annanielle/Gut-microbiota-in-Psoriasis), audited tree `dc3d91ea299057bc02d9ef2a4e505aa394523e50`, including `phyloseq-psor-control-metagenome.RData`. Read only trusted, reviewed data through an isolated data importer; do not execute downloaded notebook code merely to extract labels. Verify the 47 controls against the archive. No patient with TB or unknown phenotype becomes healthy by project membership.

### 6.3 Publisher assets

- E01 [clinical/annotation workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12967-026-08013-4/MediaObjects/12967_2026_8013_MOESM2_ESM.xlsx): `Table S1` has 143 phenotype rows and family identifiers; `TableS2` describes 28 pairs. `Table S3` is viral annotation, not a FASTA; `Table S4` is Fimenecus annotation, not measured blood metabolites. Preserve missing BMI.
- E08 workbook Tables 3/5 can fill exact per-feature statistics after identifying their headers; Table 8 is a different, 28-feature classifier. Do not expand the eleven-feature panel without a new version.
- E14 source Table S12 joins 48 stool participants; gene expression accession GSE150851 is not another stool cohort. Tables S1/S3 contain taxon/pathway findings; their sparse fold changes do not become weights.
- E03’s deposit version 1 and version 2 are versions of the same cohort. Downloadable “raw analysis data” does not necessarily mean FASTQ. Detect actual MIME/file type and inspect inventory before assigning `raw_reads_verified`.

Pin every downloaded asset with source URL, retrieval date, actual byte size, SHA-256, license, and parsing version. Preserve archive-provided checksums as well. Treat HTML error pages, corrupt workbooks, missing columns, checksum mismatch, ambiguous phenotype joins, and unexpected taxa as explicit import failures. Do not invent expected checksums or pin mutable `latest` database identifiers.

States are separate: `publication_declared`, `metadata_verified`, `download_verified`, `phenotype_join_verified`, `reprocessed_compatible`, `author_data_required`, `access_blocked`. Code availability, data availability, license permission, analytical compatibility, and clinical validity are independent fields. The coding agent must not translate one into another.

### 6.4 Compatibility and source reproduction

Reprocess usable source FASTQs using the actual pinned production lane for new reference/model work. A named species crosswalk does not make MP2, MP3, MP4, Kraken2, and 16S relative abundances interchangeable. Keep a separate source-reproduction lane when reproducing published analysis, and a production-equivalent lane for sample scoring.

Recheck sample identity, cross-project duplication, household membership, assay, and medication/IBD status before fitting. Report direction agreement after reprocessing rather than forcing the production result to match the publication. A negative replication remains evidence. Do not tune the source panel on the user’s five samples.

## 7. Results API and report layout

### 7.1 Canonical result object

Extend the existing result schema with these fields or exact equivalents. Serialize nonfinite numbers as validation errors, never NaN/Infinity. Example below is an **unavailable fixture**, not an analyzed participant:

```json
{
  "profile_family":"plaque_psoriasis",
  "panel_id":"PSO_DENG2026_SGB4348_V1",
  "panel_version":"1.0.0",
  "score_revision":"directional_fraction_v42",
  "status":"reference_unavailable",
  "index":null,
  "index_label":"Research directional pattern index",
  "usable_count":0,
  "panel_count":1,
  "coverage":0.0,
  "interval_coverage":0.0,
  "feature_mask":[],
  "full_panel_bounds":[0.0,100.0],
  "bounds_type":"coverage_and_censoring",
  "index_scope":"available_point_feature_mask",
  "bounds_scope":"full_declared_panel",
  "control_percentile":null,
  "case_percentile":null,
  "disease_probability":null,
  "clinical_classification":null,
  "reference_bundle_id":null,
  "reference_scope":null,
  "assay_transport":"native_shotgun_source",
  "source_population_transport":"not_established",
  "disease_specificity":"not_established",
  "evidence_tier":"same_center_tested_candidate",
  "source_ids":["E01","E02","E03"],
  "counterevidence_ids":["E03","E05"],
  "reason_codes":["no_compatible_frozen_reference"],
  "feature_results":[],
  "intervention_evidence_ids":[],
  "domain_assessment":"limited_metadata_based"
}
```

`native_shotgun_source` describes the evidence assay, not proof of exact production transport. `source_population_transport` and `reference_compatibility` must be separate. A trained OOD model is not implied by clipping z or checking age; absent such a model, keep `domain_assessment=limited_metadata_based`.

Result statuses:

| State | Numeric behavior | UI behavior |
|---|---|---|
| `computed_research` | Index, feature contributions, coverage, bounds | Research badge and source contrast always visible |
| `computed_limited_features` | Index from fewer than all declared slots or a one-feature panel | Explicit `m/K features`; never imply comprehensive signature |
| `interval_only` | Index null; possibly informative bounds | Show “only bounded measurements available” |
| `reference_unavailable` | Index null | Explain which compatible reference is absent |
| `population_unsupported` | Participant-facing index null | State population limitation; retain observations |
| `assay_incompatible` | Index null | Identify body-site/assay/denominator mismatch |
| `qc_failed` | Index null | Show failed QC reason; no disease interpretation |
| `features_unavailable` | Index null | Distinguish mapping, detection, depth, and not-assayed causes |
| `insufficient_stool_signature` | Index, percentile, probability null | SEBD evidence remains visible; never “no disease” |
| `strain_bundle_pending` / `blocked_unresolved_preprocessing_and_artifacts` | No purported model result | Exact missing asset/preprocessing reason |

Apply gates in order: structural integrity → body site/assay → QC → population/exposure → reference compatibility → feature mapping/detection → calculation. Preserve all relevant reason codes, not only the first one. For an unavailable nonnumeric task, do not invoke a zero-feature kernel or manufacture coverage from hypothesis taxa.

### 7.2 Front-loaded report

Preserve the existing summary and “what stood out,” then put **all skin profile readings** in the front visual section before detailed explanations. Allow multiple pages; never truncate to fit one page, “top eight,” or an arbitrary top-N list.

Each numeric panel row shows: condition/contrast, panel name, index on a clearly labeled 0–100 scale, `m/K` features, bounds, source evidence tier, and direct link to its complete detail section. For one-feature panels, write “single-marker result.” For transported ATD/HS panels, write “16S-derived pattern applied to shotgun measurements.” A new assay lane must not look equivalent to an existing score solely because both use 0–100.

Use a neutral sequential palette for concordance and a marker for the point estimate. With complete point coverage, a horizontal interval can show the same panel’s censoring bounds. With incomplete point coverage, show two labeled tracks: **measured-feature index** and **possible full-panel range**; do not imply they are the same estimand or that the second is a confidence interval around the first. Label low→high **pattern concordance**; do not use green “healthy” and red “diseased.” Show null values as em dash plus a short explanation, never as zero or a green check. The SEBD row states “stool-DNA signature not established.” Graphics must remain readable in grayscale and have text equivalents.

Detail pages include:

1. What was actually measured, source assay/population, contrast, and the meaning of the index.
2. Every declared feature, observed value/state, reference position, study direction, contribution, and mapping/transport limitation—including missing slots.
3. Opposing evidence and known shared markers; SGB4348/PF and pooled PsO/PsA are prominent examples.
4. Existing general gut function and organism readings linked as context, with no extra disease votes.
5. Relevant intervention evidence or an explicit “no signature-matched intervention evidence identified” state.

Do not turn a high shared psoriasis index into PsA, infer skin location from stool, infer a clinical stage from abundance, or hide the absence of a KO assay behind a numeric taxon proxy. Do not describe taxon depletion as something that must be “replaced.”

### 7.3 Longitudinal comparisons and donor views

Compare indices only with the same panel version, reference, database, denominator, exposure policy, and feature mask. If masks differ, recompute both using a named common mask and its own bounds, or state that direct comparison is unavailable. Index changes are not validated treatment response.

Preserve existing access controls for FMT research views. These skin modules cannot issue donor pass/fail, eligibility rank, treatment selection, or a reassurance that FMT is safe. A low score cannot exclude the disease; an unavailable score is not clearance. This is a product output constraint, not a claim that a particular donor or treatment is unsafe.

## 8. Intervention-evidence integration

Reuse the existing intervention registry. Seed the specific records below, with exact strain/formulation identity and studied indication. Evidence can be retrieved for a **confirmed diagnosis or explicitly selected educational topic** independently of a score. Profile resemblance may retrieve a general research explanation but cannot establish treatment eligibility, predicted benefit, or the need to eradicate organisms. These additions do not replace the existing disease-management evidence library.

### I01 — Psoriasis adjunct probiotic trial

[Navarro-López 2019, DOI 10.2340/00015555-3305](https://www.medicaljournals.se/acta/content/html/10.2340/00015555-3305): 90 adults, plaque psoriasis, 12-week randomized placebo-controlled adjunct trial. Product contained **B. longum CECT7347, B. lactis CECT8145, L. rhamnosus CECT8361**, 1:1:1, total 10^9 CFU/day. Both groups received topical treatment. PASI75 response was 66.7% versus 41.9%; the endpoint is **at least 75% improvement**, despite misleading abstract wording. The protocol retained achieved success at later visits. Post-trial relapse follow-up was outside the original protocol. This is formulation-specific adjunct evidence, not proof that a stool profile predicts response or that any probiotic blend works.

Store the exact studied regimen in the internal evidence record; participant cards summarize the formulation, study duration, outcome and limitations without generating a personalized dose. Do not promise that treatment normalizes SGB4348 or any v4.2 feature. Trial adverse-event reporting does not establish universal safety.

### I02 — Oral ST11 and dandruff

[Reygagne 2017, DOI 10.3920/BM2016.0144](https://pubmed.ncbi.nlm.nih.gov/28789559/): 60 adult men with moderate/severe dandruff, ages 18–60, ST11 (**Lactobacillus paracasei NCC2461**, modern genus Lacticaseibacillus) versus placebo for 56 days, 10^9 CFU/day in the study. Dandruff/erythema and scalp measures improved. The phenotype was dandruff and the microbiome outcome was scalp; it does not establish a gut SEBD signature or validate stool-based probiotic selection.

Link as `related_phenotype=dandruff`, `applicability_to_sebd=indirect`, never an exact proven SEBD treatment match. Strain identity matters; neither other paracasei strains nor generic fermented foods inherit the trial effect.

### I03 — Topical SEBD probiotic pilot

[2024 oily suspension study, DOI 10.1038/s41598-024-53016-0](https://www.nature.com/articles/s41598-024-53016-0) studied 25 SEBD patients using **topical** L. crispatus P17631/Lacticaseibacillus paracasei I1688 suspension for one week. Its route and skin outcome must remain explicit. It is not an oral regimen, not a stool intervention-response classifier, and does not justify gut antifungal treatment.

For ATD, vitiligo, HS and PsA, retrieve applicable existing condition-specific evidence. If no curated record matches, render the explicit no-supported-signature-matched-intervention state. Do not fabricate a recommendation to supply every taxon row. Missing bacteria, elevated skin-pattern indices, KO annotations, or inferred dysbiosis do not trigger antibiotics, botanical antimicrobial regimens, systemic antifungals, prescription changes, or FMT. Do not claim FMT is the only solution. Maintain the project’s clinician-review and pediatric rules.

Minimal intervention record fields: `id`, `source_url`, `phenotype`, `confirmed_diagnosis_required`, `study_population`, `route`, `exact_product_or_strains`, `study_regimen_internal`, `comparator`, `cointerventions`, `clinical_endpoint`, `microbiome_endpoint`, `result_direction`, `limitations`, `signature_response_validation=false`, and `participant_copy`. Exact doses are evidence metadata, not an instruction to the participant.

## 9. Build order, validation, and completion

### 9.1 Implementation sequence

1. Inspect and snapshot the actual registry/results; establish semantic family mappings and ID collision tests.
2. Add the embedded evidence records, nine numeric-panel manifests, three unavailable tasks, and non-scored counterevidence/functional observations.
3. Implement or adapt typed observation and taxonomy resolution, eligibility gates, frozen reference construction, and the exact numeric kernel. Fit compatible existing references if available; record real blockers otherwise.
4. Build metadata importers for the prioritized public cohorts and verify joins before scheduling bulk reads. Store deterministic provenance and quality failures. Do not execute arbitrary downloaded analysis scripts.
5. Integrate the new result objects and complete front-loaded summary/detail pages, including null/partial states and intervention cards.
6. Run the numerical self-test, the relevant regression/source-integrity tests below, and render representative complete/partial/unsupported/SEBD reports. Inspect every added page for missing rows, cut-off labels, broken links, and misleading color/percent signs.
7. Deliver the code changes with a concise implementation status: which modules calculate on current assay/reference inputs, which are unavailable and why, which public joins are verified, and which disease-validation claims remain unsupported. Do not claim completion by adding only prose cards for supported numeric panels.

### 9.2 Acceptance tests

These are meaningful gates for biomedical scoring and data integrity. Translate them into the actual test framework; do not make mocked author data pass a real-data gate.

| ID | Fixture / action | Required outcome |
|---|---|---|
| V42-01 | Parse seed JSON | Nine numeric panels; 52 slots; unique panel IDs and per-panel feature IDs |
| V42-02 | Migrate a registry already containing psoriasis/ATD | Semantic deduplication; existing family identity preserved |
| V42-03 | Registry includes Alzheimer AD and dermographism SD | Neither is renamed or overwritten by ATD/SEBD |
| V42-04 | Existing unrelated sample report before/after | Existing score values/versions unchanged unless explicitly migrated |
| V42-05 | A species with identical epithet in different genera | Bacteroides finegoldii never maps to Alistipes finegoldii |
| V42-06 | Missing SGB4348 in MP3 reference | Unavailable feature, not zero, Firmicutes proxy, or guessed SGB |
| V42-07 | Same numeric SGB in another namespace | Reject absent verified biological crosswalk |
| V42-08 | Source Megasphaera_unclassified but only genus total | That slot unavailable; do not use genus total |
| V42-09 | Generic P. copri observation for vitiligo clade C | Clade slot unavailable |
| V42-10 | Resolve eight adult ATD genera but no five source groups | Point coverage 8/13; unresolved slots remain in bounds |
| V42-11 | Unclassified-family source feature | Entire family sum rejected |
| V42-12 | Historical F. prausnitzii umbrella includes F. duncaniae | Ambiguous pair quarantined; no double-counting |
| V42-13 | Luan thermophiles/thermophilus typo | No fuzzy genus substitution or active marker |
| V42-14 | Stool culture shows lower absolute E. coli count | No negative DNA-relative coefficient automatically generated |
| V42-15 | E11 67 enrolled/printed percentages near32 denominator | Preserve denominator discrepancy; do not relabel 32 as controls |
| V42-16 | SEBD MR association loaded | Evidence visible; numeric coefficient forbidden |
| V42-17 | Micrococcus skin association or Malassezia scalp abundance | Does not enter any stool signature |
| V42-18 | SEBD clinical-series “60% dysbiosis” | Never becomes a prior, prevalence estimate for screening, or cutoff |
| V42-19 | SEBD task runs without a supported stool signature | Null with explicit reason, not0/50 or a negative diagnosis |
| V42-20 | All 13 ATD observations missing | Index null; full-panel bounds [0,100] |
| V42-21 | One feature at its frozen reference median | Contribution50; never interpreted as50% disease probability |
| V42-22 | One point50 of two slots, other unavailable | Index50, coverage.5, full-panel bounds[25,75] |
| V42-23 | Abundance rises for + versus − features | Opposite monotonic contributions; complementary sum100 |
| V42-24 | Nondetection with limit but no estimator | Interval-only contribution; no invented point |
| V42-25 | Missing row without a detection limit | Unavailable rather than maximally depleted |
| V42-26 | Negative/NaN/Infinity/boolean/string abundance | Reject before scoring |
| V42-27 | Duplicate participant/control technical replicates | Fitting rejects inflated reference N |
| V42-28 | Nineteen eligible controls or <20 per feature | No fitted reference for that feature |
| V42-29 | Query alone versus same query in a larger inference batch | Exactly identical output; no inference-time fitting |
| V42-30 | Corrupt reference digest or scale below floor | Refuse calculation |
| V42-31 | Selected-feature reclosure or percent/fraction mismatch | Reject incompatible denominator or convert declared percent once |
| V42-32 | KO assay absent but taxa measured | PsA KO index unavailable; no taxa-to-function imputation |
| V42-33 | KO total denominator zero / mixed stratification | Reject functional calculation |
| V42-34 | K07114 versus typo K07714; microbial crp | Correct KO only; no host chloride-channel/serum-CRP inference |
| V42-35 | Age8,17,unknown; E15 age65 | Adult/E15 numeric gates fire; measurements remain visible |
| V42-36 | Recent antibiotics or unknown medication metadata | Apply explicit compatibility rule/unknown flag; never assume clean |
| V42-37 | Scalp, blood, or culture sample supplied to stool kernel wrapper | Assay/body-site gate rejects |
| V42-38 | 500,000 reads mislabeled as pairs | Cannot satisfy pair-count floor |
| V42-39 | Low Shannon alone with no signature taxa | No added skin-profile score or disease-positive flag |
| V42-40 | E01 SGB4348 evidence card | Includes E03 PF sharing and limited specificity |
| V42-41 | Sá psoriasis-versus-HS taxon differences | Not imported as psoriasis-versus-healthy directions |
| V42-42 | Sá HS-versus-HC null finding | Retained as counterevidence, not counted as positive replication |
| V42-43 | E04 pooled PsO/PsA q-values | Attached only to pooled contrast, never PsA-only or PsO-only |
| V42-44 | E04 E.rectale AUC; E08 random-forest AUC | Not attributed to newly engineered multi-feature indices |
| V42-45 | Deng text/supplement PR-AUC disagreement | Preserve conflict; do not silently select favorable number |
| V42-46 | PWY-5005 failure to validate; KEGG00650 q=.157 | Zero diagnostic weight; no significant replicated label |
| V42-47 | PRJNA1061168 project imported | Exact healthy whitelist required; no all70=healthy assumption |
| V42-48 | Chang alias/PID and Xiao alias/host_subject_id differ | Explicit attributes and source tables determine joins |
| V42-49 | HS reanalysis or dataset versions imported twice | One participant/cohort identity; no duplicated replication |
| V42-50 | Source metadata conflict, HTML404 disguised as XLSX | Import fails with specific reason; no fabricated fallback |
| V42-51 | Viral annotation names without sequences | No vBin score, no generated reference FASTA |
| V42-52 | Seong infant strain findings or Peng leaky-CV AUC | No generic adult species or pediatric WGS classifier claim |
| V42-53 | Whole family has several study indices | No unversioned fusion or inflated disease count |
| V42-54 | Different masks/references at two dates | Common-mask recomputation or explicitly incomparable |
| V42-55 | Full summary exceeds one page | Paginate all rows; no top-N truncation |
| V42-56 | Null/partial/single-marker report visual | Correct labels/bounds; no green clearance or probability sign |
| V42-57 | ST11 and topical SEBD trial cards | Correct strain, route, dandruff-versus-SEBD distinction |
| V42-58 | High profile index with no confirmed diagnosis | No treatment eligibility, dose, eradication, or FMT trigger |
| V42-59 | Donor research view with low/all-null skin indices | No donor pass/fail, rank, or safety clearance |
| V42-60 | Current five personal samples used as fixtures | Never enter training, reference fitting, or threshold tuning |
| V42-61 | Partial-mask point index lies outside a tightened full-panel range | Both scopes labeled separately; renderer does not present the range as that index’s confidence interval |

### 9.3 Research promotion gates

Before claiming disease specificity or validated screening, evaluate independent geographic/site cohorts using the locked pipeline; compare against relevant skin diseases, IBD, medication exposures and healthy controls; report sample/participant counts, calibration, uncertainty, ROC/PR metrics at actual test prevalence, clinical-only baselines, and errors in subgroup transport. Do not use healthy-control discrimination alone to claim that psoriasis and PsA or SEBD can be distinguished.

For infant strains, require the exact MAG/marker reference, source subclade labels, sequence quality, unique-marker alignment criteria, strain coverage/mixture behavior, and independent age-matched validation. A genus or species observation is insufficient. Preserve the pending task rather than generating an ad hoc phylogenetic cutoff.

For SEBD, the present research establishes gut-related observations and hypotheses, not an operational stool-DNA signature. New measured evidence can activate a versioned panel through the same typed contracts; no new architectural rewrite is needed.

**Definition of done:** the specified supported panels execute deterministically when compatible measurements and references exist; unavailable tasks remain visible with truthful reasons; every source, feature, data join, evidence limitation, reporting rule, and acceptance gate necessary for this release is contained in this document.
