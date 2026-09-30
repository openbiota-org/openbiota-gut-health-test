# BUILD SPEC 2 — Disease similarity profiles

Follow-on to `BUILD_SPEC.md`. Read that first; this assumes the v1.1.0 pipeline is
working and extends it.

---

## Overview

### Where the tool is now

v1.1.0 measures **gene capacity**: DIAMOND translated search of every read against
curated protein references, normalized to copies per 100 bacterial genomes via `rpoB`,
across nine pathways and nineteen genes. It reports percentiles against a reference
cohort, validates against synthetic communities built from known genomes, confirms
depth saturation by rarefaction, and renders a PDF plus `results.json`.

That answers: *which metabolite-producing genes does this gut community carry?*

### What this spec adds

A second axis: **which organisms are present, and how closely does this community
resemble published disease-associated patterns?**

Four components, in dependency order:

1. **An expanded reference cohort** — the prerequisite for everything else.
2. **A taxonomic engine** (MetaPhlAn 4) alongside the existing DIAMOND engine.
3. **A profile system** — declarative configs describing a published
   disease-associated microbiome pattern.
4. **A scoring engine** that turns a profile plus a sample into module scores,
   percentiles, and a confidence grade.

Three profiles ship: **colorectal cancer** (positive control), **ME/CFS**, and
**Long COVID**.

### What the new output looks like

```
ME/CFS-LIKE COMMUNITY SIMILARITY                      Sample A02
═══════════════════════════════════════════════════════════════════
  Taxonomic module        +0.41    68th pct    12/14 features
  Functional module       +0.22    59th pct    butyrate, tryptophan
  Ecological module       −0.10    44th pct    low weight by design
  ───────────────────────────────────────────────────────────────
  Combined similarity              63rd percentile
  Confidence                       MODERATE

  CROSS-ENGINE CHECK                                    CONCORDANT
    Taxonomic butyrate producers    −0.3 SD vs reference
    Gene panel (but + buk)          82.8 copies/100 genomes, 71st pct
    Both engines agree on butyrate capacity.

  GENERAL DYSBIOSIS ANCHOR (GMWI2)              +0.4  not dysbiotic
    A strongly negative anchor would mean the score above reflects
    general gut disturbance rather than anything profile-specific.

  Research use only. This is resemblance to a published group-level
  pattern. It is not a diagnosis and not a probability of disease.
═══════════════════════════════════════════════════════════════════
```

### What "63rd percentile" means

The sample is more concordant with this profile's prespecified pattern than 63% of the
matched reference set. It is not a 63% chance of having the condition. Every piece of
the design below follows from holding that distinction.

---

## PART 1 — Reference cohort expansion

**Build this first. Nothing downstream is meaningful without it.**

### 1.1 Why it comes first

v1.1.0's percentiles rest on 33 samples subsampled to 600,000 read pairs. That is
adequate for ranking a single pathway value and inadequate for what follows: a
percentile at the tails is unstable at n=33, taxonomic profiling needs full depth
rather than a 600k subsample, and disease profiles require references *matched* on
country, age band, and sex — which means several hundred samples so that each matching
stratum still has enough members.

### 1.2 Requirements

Source healthy adult stool metagenomes from **curatedMetagenomicData**, which supplies
uniformly processed profiles with subject metadata attached.

- Minimum 300 samples; target 1,000+.
- Record per sample: country, age, sex, BMI where available, study accession,
  sequencing platform, read depth.
- Process at full depth through the identical pipeline — same DIAMOND reference set
  fingerprint, same MetaPhlAn database version.
- Cache the computed reference distribution as a versioned artifact with a manifest
  recording cohort composition, sample count, all tool versions, and both the DIAMOND
  reference fingerprint and the MetaPhlAn database version.
- Support stratified subsetting at query time: `--match country,age_band,sex`. When a
  stratum holds fewer than 30 samples, fall back to the full cohort and mark the
  result as unmatched in the confidence grade.

### 1.3 Reference contrasts for profiles

Profiles need two contrasts, not one:

- **Primary:** where available, people with the same exposure who did *not* develop the
  syndrome — recovered-without-PACS, for instance. This isolates syndrome-specific
  signal from post-exposure signal.
- **Secondary:** general healthy adults, matched on the non-exposure variables.

For most published cohorts the primary contrast is unavailable. When it is missing,
say so explicitly in the report and downgrade confidence rather than substituting
healthy controls silently.

---

## PART 2 — Taxonomic engine

### 2.1 Why a second engine

The DIAMOND engine answers "how much of gene X is present." Disease profiles need
"how much of species Y is present," which requires marker-gene profiling. The
`rpoB`-derived composition already in section 4 of the report is a useful by-product
but is marker-limited — reliable at phylum level, indicative at genus. Profile scoring
needs species and strain resolution.

### 2.2 Implementation

```
openbiota/engines/
├── diamond.py       # existing
└── metaphlan.py     # new
```

Add a required `engine:` field to the panel and profile schemas; existing panels
declare `engine: diamond`. Panels and profiles of different engine types must run in a
single invocation and merge into one report.

Requirements:

- **Pin the MetaPhlAn database version and record it in `results.json` and the PDF.**
  SGB identifiers are database-version dependent — `Faecalibacterium_SGB15346` is
  meaningless without the database that produced it, and profiles computed against
  different databases are not comparable.
- Cache the `--bowtie2out` intermediate so re-profiling does not re-align.
- Run on the full sample, not a subsample. The 600k subsampling used for the reference
  cohort in v1.1.0 must not carry over here.
- Roughly 15 GB for the Bowtie2 index — within budget, but it is the largest memory
  consumer in the pipeline.
- Same native-build consideration as DIAMOND on Apple Silicon.

### 2.3 Host read handling

The existing pipeline needs no host removal, because human reads carry neither pathway
genes nor `rpoB` and cancel in the ratio. **MetaPhlAn has no such property.** Add
KneadData/Bowtie2 host filtering against GRCh38 + CHM13 + PhiX ahead of the taxonomic
engine only. Report host fraction as a QC field. Discard host-aligned files per the
project's data protocol.

### 2.4 Species-level discipline

Profiles specify species or strain. **Never substitute a genus for its constituent
species.** Aggregation reverses direction in practice — genus-level *Faecalibacterium*
has been reported higher in a cohort where species-level *F. prausnitzii* was lower,
and qPCR has found one *F. prausnitzii* strain lower and another higher against the
same outcome. The loader must reject a profile feature that declares
`level: genus` for any taxon on the conflicted list in 4.4.

---

## PART 3 — Profile system

### 3.1 Schema

```yaml
name: mecfs
version: 0.1.0
label: "ME/CFS-like gut microbiome similarity"
status: research_only
duration_dependent: true

references:
  primary:
    source: PRJNA878603
    note: "healthy controls n=79; no recovered-ME/CFS group exists"
  secondary:
    source: curatedMetagenomicData
    match_on: [country, age_band, sex]

modules:
  taxonomic:
    weight_in_combined: 0.50
    features:
      - name: Faecalibacterium_prausnitzii
        level: species
        direction: decreased
        w: 1.0      # evidence weight
        q: 1.0      # cohort independence
        s: 0.5      # disease specificity
        c: 1.0      # direction agreement
        conflict_note: "genus Faecalibacterium conflicts; species only"
        sources: ["Guo 2023 CHM", "Xiong 2023 CHM"]
      - name: Eubacterium_rectale
        level: species
        direction: decreased
        w: 1.0
        q: 0.7
        s: 0.6
        c: 1.0
        sources: ["Guo 2023 CHM"]

  functional:
    weight_in_combined: 0.35
    features:
      - name: butyrate
        engine: diamond          # reuses the existing v1.1.0 panel
        direction: decreased
        w: 1.0
        q: 1.0
        s: 0.5
        c: 1.0

  ecological:
    weight_in_combined: 0.15
    features:
      - name: shannon_diversity
        direction: decreased
        w: 0.5
        s: 0.2

cross_engine_checks:
  - taxonomic_feature: Faecalibacterium_prausnitzii
    functional_feature: butyrate
    expect: concordant
    on_disagreement: downgrade_confidence

abstain_if:
  - recent_antibiotics_days_lt: 90
  - usable_nonhost_reads_lt: 500000
  - missing_metadata: [onset_date, sampling_date]

combined_split_note: >
  The 50/35/15 module split is an engineering starting point, not
  evidence-derived. Do not present it as validated weighting.
```

### 3.2 Versioning

Every profile carries a semantic version. Changing any weight bumps it. `results.json`
records the profile version, MetaPhlAn database version, DIAMOND reference
fingerprint, reference cohort manifest ID, and pipeline version — so any score can be
recomputed exactly, and two results are only comparable when all five match.

---

## PART 4 — Scoring engine

Generic. It consumes a profile and a sample and has no knowledge of which disease it
is scoring.

### 4.1 Per-feature robust value

Using only the matched reference distribution:

```
z_i = (CLR_i(sample) − median_i(reference)) / (1.4826 × MAD_i(reference))
v_i = d_i × clip(z_i, −3, +3) / 3
```

`d_i` is +1 for a reported increase, −1 for a decrease. `v_i` runs from −1 (opposite
the pattern) to +1 (strongly concordant).

Median and MAD rather than mean and standard deviation, because microbiome abundance
distributions are skewed enough that one outlier otherwise dominates the scale.
Centered log-ratio because relative abundances are compositional — they sum to 1, so
an apparent increase in one taxon can be entirely produced by decreases elsewhere, and
sums or differences of raw abundances are invalid. Retain raw relative abundance
alongside CLR for the human-readable report.

The existing `copies per 100 genomes` values feed the functional module directly; they
are already depth-normalized and need only the same median/MAD treatment against the
reference cohort.

### 4.2 Module aggregation

```
module = Σ(wᵢ qᵢ sᵢ cᵢ vᵢ) / Σ(wᵢ qᵢ sᵢ cᵢ)
```

| Factor | Meaning | Low when |
|---|---|---|
| `w` | Evidence weight | Small cohort, single study |
| `q` | Cohort independence | Replicated only within one recruitment family |
| `s` | Disease specificity | Feature also marks IBD, IBS, antibiotics, hospitalization |
| `c` | Direction agreement | Cohorts disagree on direction |

All 0–1, all in the versioned profile manifest, all reproducible by hand from the
manifest plus the measured values.

**The `q` factor matters more than it looks.** Several accessions in the Long COVID
literature — `PRJNA714459`, `PRJNA876804`, `PRJNA841786`, `PRJNA1004982`,
`PRJNA995807`, `PRJNA1056888` — come from potentially overlapping recruitment programs
at one institution. A taxon reported across four of those papers carries roughly one
independent vote. Encode that as `q`, and hold out entire studies rather than pooling
and randomly splitting them during any validation.

### 4.3 Four modules

1. **Taxonomic** — species and strain features only.
2. **Functional** — existing DIAMOND pathway values, plus HUMAnN 3 MetaCyc pathways
   where a profile needs them.
3. **Ecological** — Shannon, richness, butyrate-producer fraction, distance to matched
   controls. **Deliberately low weight**: diversity findings conflict across cohorts,
   and diversity shifts with diet, travel, and antibiotics far more readily than any
   disease signal.
4. **Phenotype** — symptom-domain submodules, used only when the symptom label was
   declared before analysis. Excluded from the combined number, so a symptom pattern
   cannot be chosen after seeing the data.

Report all four separately. If a combined number is produced, initialize at 50%
taxonomic, 35% functional, 15% ecological, and state in both output and README that
this split is an engineering starting point rather than an evidence-derived weighting.

Convert to a 0–100 percentile against an untouched reference set. Never by linear
rescaling of the raw value.

### 4.4 Conflicted features

These have no stable direction across cohorts and must never carry a fixed global
direction at genus level: *Faecalibacterium*, *Streptococcus*, *Veillonella*,
*Prevotella*, *Blautia*, *Lachnospiraceae*, *Alistipes*.

Documented conflicts to encode as `c` penalties and surface in the report:
*Streptococcus* and *Veillonella* are enriched in Dutch and Hainan cohorts but depleted
in GI-predominant long COVID; *Blautia* and *Lachnospiraceae* were enriched in a US
long COVID cohort while treated as beneficial elsewhere; *Alistipes* direction varies
by sex and neurological phenotype in recent work.

### 4.5 General dysbiosis anchor

Run **GMWI2** on every sample and report it beside every profile score. It is a
published, externally validated general gut-health index trained on 8,069 metagenomes
across 54 studies.

Its role is specificity. If GMWI2 is strongly negative, the community is generally
disturbed, and a high score on any profile probably reflects that rather than anything
disease-specific. The report must say so in that case rather than presenting the
profile score alone.

### 4.6 Confidence and abstention

Grade on: fraction of planned weighted features measured; reference cohort size and
match quality; whether pipeline, profiler, and database versions match the reference;
host-free usable depth; missing antibiotic, medication, stool-form, and timing
metadata; out-of-distribution distance; and agreement between the taxonomic and
functional modules.

**Abstain rather than emit a number** when: broad-spectrum antibiotics within the
prespecified window; usable non-host reads below the minimum (start at 500,000 and
sensitivity-test it); the sample falls outside the multivariate reference envelope;
onset or sampling date missing; or the taxonomic and functional modules point opposite
ways with no phenotype explanation.

Abstention is a feature, not a failure path. A refused score is more useful than a
confident wrong one, and the existing depth-saturation machinery already provides most
of the QC signal needed to decide.

---

## PART 5 — The three profiles

### 5.1 Colorectal cancer — positive control

**This is a pipeline validation profile, not a cancer screen.** Validated clinical
tests exist — FIT, multi-target stool DNA, colonoscopy — and this is not one of them.
The report must state that plainly and must not present the output as a risk estimate.

It ships because CRC is the only microbiome disease association with genuinely
cross-cohort-validated performance, and because curatedMetagenomicData carries CRC
cohorts **with disease labels attached**. That makes it the one profile on which the
entire scoring engine can be validated end to end against ground truth.

If the engine cannot recover a known-good CRC signal on labeled data, no number it
produces for any other profile means anything.

**Taxonomic features:** *Fusobacterium nucleatum*, *Parvimonas micra*,
*Peptostreptococcus stomatis*, *Gemella morbillorum*, *Solobacterium moorei* — all
enriched, replicated across independent international cohorts.

**Functional features, via the existing DIAMOND engine:** the *pks* island (*clbB*,
colibactin) and *Bacteroides fragilis* toxin (*bft*). These are direct gene targets
needing no taxonomic engine, so this profile exercises both paths and can be partially
run before Part 2 is finished.

### 5.2 ME/CFS

The best-evidenced profile in the set. Two cohorts from different groups, published
simultaneously, with functional confirmation by qPCR and measured short-chain fatty
acids rather than DNA inference alone.

**Guo et al. 2023** — 106 cases, 91 controls, geographically diverse US sites. Reduced
*F. prausnitzii* and *E. rectale* contributing to butyrate deficiency; low
*F. prausnitzii* correlates with more severe fatigue.

**Xiong et al. 2023** — short-term (<4 years, n=75) and long-term (>10 years, n=79)
patients against 79 controls. Reduced microbial butyrate biosynthesis; reduced plasma
butyrate, bile acids, and benzoate. Data at `PRJNA878603`, **code published at
github.com/ohlab/MECFS_2021**.

**Implement the duration dependence.** Xiong found short-term patients had significant
microbial dysbiosis while long-term patients had largely *resolved* microbial
dysbiosis, retaining metabolic and clinical abnormalities. Consequence: a low
taxonomic score in someone ill more than ten years is consistent with the published
pattern, not evidence against it. Carry `duration_dependent: true`, keep separate
references per stratum, never average across them, and state this in the report. If
duration is unknown, report both strata and downgrade confidence.

**Do not cite the 90% accuracy figure** from the 2025 BioMapAI paper next to this
profile. That figure comes from integrating metagenomics with plasma metabolomics,
immune cell profiling, and blood tests; microbiome data alone best predicted
gastrointestinal, emotional, and sleep symptoms specifically. A stool-only tool cannot
reproduce it.

**This profile carries the most important cross-engine check in the tool.** The
taxonomic module claims butyrate-producer depletion; the existing `but`/`buk` panel
measures butyrate gene capacity directly — 82.8 copies per 100 genomes on sample
A02, from 1,917 fragments. Two independent measurements of the same biology. If
they disagree, something is wrong in the pipeline, the sample, or both. Print the
comparison and trigger reduced confidence on disagreement. Build the pattern here and
reuse it.

### 5.3 Long COVID

Weakest of the three. Implement last, and let the report reflect the evidence.

**Two prediction tasks, never mixed.** This profile scores *chronic-state
resemblance* — stool at least three months post-infection with ongoing symptoms.
*Acute future-risk* prediction from stool at 0–2 months is a different task; features
from it must never be read as biomarkers in a chronic sample. If built, it is a
separately named and separately validated profile.

**Calibrate expectations to these numbers.** In the largest independent prospective
outpatient study (799 participants), species-level AUROC was 0.62 and genus 0.59,
against 0.72 for clinical variables alone — adding microbiome features did not
materially improve on a symptom questionnaire. A community study of 2,561 participants
found no association between gut microbiome profile and illness duration. A separate
cohort found lower richness but no species surviving multiple-testing correction. The
frequently quoted AUC of 0.96 is for classifying a study's own internally defined
enterotypes, not an externally validated diagnosis.

**Feature tiers.**

*Tier A, weight 1.0 before penalties:* higher species-level *R. gnavus*, *C. bolteae*,
*C. innocuum*; lower species-level *F. prausnitzii*, *B. adolescentis* and
*B. pseudocatenulatum*; lower community butyrate and SCFA capacity. Most originate in
one institutional family — apply the `q` penalty.

*Tier B, weight 0.5:* higher *E. ramosum*, *B. vulgatus*, *F. plautii*; lower
*Gemmiger formicilis*; lower indole/tryptophan and L-isoleucine biosynthesis; higher
urea and arginine degradation.

*Tier C, 0.25 or phenotype module only:* GI-predominant oral–gut pattern;
respiratory-associated *Streptococcus* species; *A. naeslundii* and *C. innocuum* in
fatigue and neuropsychiatric symptoms.

Apply specificity penalties: *B. vulgatus* is common across many conditions; flagellin
potential is shared with rheumatoid arthritis, ankylosing spondylitis, and IBD;
*E. coli*, *Klebsiella*, and antimicrobial resistance genes track ICU exposure and
antibiotics rather than the syndrome.

**Weight sources by what they can support.** The two-year Hainan cohort
(Zhang et al. 2026) is n=11 per arm with uninfected rather than recovered controls,
and its data availability statement reports no deposited datasets, so its SGB-level
calls are neither reusable nor database-version-portable. Set `w: 0.25` and flag it as
a hypothesis source.

---

## PART 6 — Report extension

Extend the existing PDF rather than producing a second document. Add after the
existing metabolite sections:

**Section: Community composition.** Replace the current `rpoB`-derived genus list with
MetaPhlAn species-level output. Keep the phylum summary, Shannon, evenness, and
Firmicutes:Bacteroidetes ratio. Label species calls as measured rather than
"closest match," which the marker-based approach now supports.

**Section: Profile similarity.** One page per profile. Four module bars in the visual
language already used for pathway percentiles, the combined percentile, the confidence
grade with its specific reasons, the cross-engine check result, and the GMWI2 anchor.
Below each, a decomposition table listing every contributing feature with its measured
value, reference percentile, direction, all four weight factors, and its source
citations — so the score is always reducible to what produced it.

**Section: What this cannot tell you.** Extend the existing three-item list with a
fourth: profile similarity is resemblance to a group-level published pattern, not a
diagnosis or a probability, and the underlying associations largely do not distinguish
one condition from general gut disturbance.

**Abstention rendering.** When a profile abstains, the page must show the reasons and
what would need to change, not a blank or a zero.

`results.json` gains the full profile block: every feature value, factor, module score,
percentile, confidence component, and all version identifiers.

---

## PART 7 — Validation

Ordered by what is achievable, following the pattern already established with synthetic
communities.

**7.1 Extend the synthetic harness to taxonomy.** The existing framework builds
synthetic samples from reference genomes with known gene content. Extend it to build
synthetic communities with known *taxonomic* composition, spiking profile taxa at
defined abundances. Verify MetaPhlAn recovers those abundances, and that the scoring
engine recovers the correct direction and approximate magnitude. This reuses machinery
that already works and tests the new engine against known truth. Report sensitivity,
specificity, and quantitative slope and R² as the existing report does.

**7.2 CRC end-to-end on labeled data.** curatedMetagenomicData CRC cohorts carry
disease labels. Run the full engine and report ROC AUC with confidence intervals
against published cross-cohort CRC performance. **This is the only validation that
tests the whole machine against clinical ground truth. Do it before building the other
two profiles.**

**7.3 ME/CFS per-feature reproduction.** `PRJNA878603` has published code and a key
resources table. Reproduce the reported per-feature directions, and separately check
whether the duration-stratified pattern appears. Report per feature, not as a pass/fail
aggregate — a feature that fails to replicate should have its weight reconsidered.

**7.4 Cross-profile specificity matrix.** Score every profile against IBD, IBS, CRC,
recovered post-COVID, recent-antibiotic, and healthy cohorts, and publish the matrix.
A profile scoring high on everything is measuring general dysbiosis. This single table
is the most informative output the project can produce, and it is what makes the GMWI2
anchor interpretable.

**7.5 Reference stability.** Bootstrap the reference cohort to show how percentile
estimates stabilize with cohort size, and publish the curve. This justifies the Part 1
sample-count target empirically rather than by assertion.

**A note on labels.** Several Long COVID accessions publish reads without a
sample-to-outcome key. Verify label availability before building any validation step
around a dataset, and fail loudly when labels are absent rather than substituting a
proxy. Where labels cannot be obtained, restrict use to unsupervised work and state in
the README that the profile's discrimination is unmeasured — do not print an AUC for
it.

---

## PART 8 — Sources

Read the primary source before implementing the section it supports. Cite all in the
README with working links, and verify each resolves.

### 8.1 ME/CFS

- Guo C, Che X, Briese T, et al. Deficient butyrate-producing capacity in the gut
  microbiome is associated with bacterial network disturbances and fatigue symptoms in
  ME/CFS. *Cell Host Microbe* 2023;31:288–304.
  https://doi.org/10.1016/j.chom.2023.01.004
- Xiong R, Gunter C, Fleming E, et al. Multi-omics of gut microbiome-host interactions
  in short- and long-term ME/CFS patients. *Cell Host Microbe* 2023;31:273–287.
  https://doi.org/10.1016/j.chom.2023.01.001 — data `PRJNA878603`, code
  https://github.com/ohlab/MECFS_2021
- Xiong R, Aiken E, Caldwell R, et al. BioMapAI: AI multi-omics modeling of ME/CFS.
  *Nature Medicine* 2025. https://doi.org/10.1038/s41591-025-03788-3

### 8.2 Long COVID — primary

- Liu Q, Mak JWY, Su Q, et al. *Gut* 2022;71:544–552.
  https://doi.org/10.1136/gutjnl-2021-325989 — data `PRJNA714459`
- Liu Q, Su Q, Zhang F, et al. *Nat Commun* 2022;13:6806.
  https://doi.org/10.1038/s41467-022-34535-8 — data `PRJNA876804`
- Su Q, Lau RI, Liu Q, et al. *Gut* 2023;72:1230–1232.
- Su Q, et al. The gut microbiome associates with phenotypic manifestations of
  post-acute COVID-19 syndrome. *Cell Host Microbe* 2024;32:651–660.
  https://pubmed.ncbi.nlm.nih.gov/38657605/ — data `PRJNA1004982`
- Zhang F, et al. *Gastroenterology* 2022;162:548–561.
  https://pubmed.ncbi.nlm.nih.gov/34687739/ — functional backbone: reduced SCFA and
  L-isoleucine biosynthesis, elevated urea production. Data `PRJNA689961`
- Zhang D, et al. *Probiotics Antimicrob Proteins* 2026.
  https://doi.org/10.1007/s12602-026-11197-2 — hypothesis source, `w: 0.25`

### 8.3 Long COVID — effect-size calibration

**Read these before writing any interpretive text.**

- Comba IY, et al. *Gut Microbes* 2026 — 799 outpatients; species AUROC 0.62, genus
  0.59, clinical 0.72. https://www.biorxiv.org/content/10.1101/2024.12.10.626852v1 —
  data `PRJNA1167328`
- Österdahl MF, et al. *Scientific Reports* 2023 — 2,561 participants, no association
  with illness duration
- Hamrefors V, et al. *Scientific Reports* 2024 — no species survived correction
- Upadhyay V, et al. *mBio* 2023 — associations vanished after longitudinal
  adjustment. Data `PRJNA885137`
- Robust cross-cohort gut microbiome associations with COVID-19 severity.
  *Gut Microbes* 2023;15:2242615. https://doi.org/10.1080/19490976.2023.2242615
- Blankestijn JM, et al. *IJMS* 2025;26:1781 — within-cohort heterogeneity, confounded
  by ICU exposure. Data `PRJNA1181655`
- Porcari S, et al. International consensus on microbiome testing, *Lancet Gastroenterol
  Hepatol* — evidence insufficient for routine clinical microbiome testing

### 8.4 Methods and tools

- MetaPhlAn 4: Blanco-Míguez A, et al. *Nat Biotechnol* 2023.
  https://doi.org/10.1038/s41587-023-01688-w — pin the SGB database version
- HUMAnN 3: Beghini F, et al. *eLife* 2021;10:e65088.
  https://doi.org/10.7554/eLife.65088
- KneadData: https://github.com/biobakery/kneaddata
- curatedMetagenomicData: Pasolli E, et al. *Nat Methods* 2017;14:1023–1024.
  https://doi.org/10.1038/nmeth.4468 — reference cohort and CRC validation source
- Compositional data: Gloor GB, et al. *Front Microbiol* 2017;8:2224.
  https://doi.org/10.3389/fmicb.2017.02224 — **read before implementing 4.1**
- GMWI2: Chang D, et al. *Nat Commun* 2024;15:7447.
  https://doi.org/10.1038/s41467-024-51651-9 — the dysbiosis anchor; bioconda package
  available
- GMHI: Gupta VK, et al. *Nat Commun* 2020;11:4635.
  https://doi.org/10.1038/s41467-020-18476-8 — the closest published analogue to this
  scoring design
- SIAMCAT: Wirbel J, et al. *Genome Biol* 2021;22:93.
  https://doi.org/10.1186/s13059-021-02306-1 — required if learned coefficients are
  ever fitted; naively transferred models lose accuracy and disease specificity
- MMUPHin, for cohort batch adjustment:
  https://huttenhower.sph.harvard.edu/mmuphin
- STORMS reporting checklist for microbiome studies

---

## PART 9 — Build order

1. **Reference cohort expansion** (Part 1). Publish the stability curve from 7.5.
2. **MetaPhlAn engine** with host filtering and version pinning (Part 2).
3. **Synthetic taxonomic validation** (7.1), reusing the existing harness.
4. **Scoring engine** (Part 4) with unit tests on synthetic profiles.
5. **CRC profile and end-to-end validation on labeled data** (5.1, 7.2).
   **Stop here and report the AUC before continuing.**
6. **ME/CFS profile**, including the cross-engine concordance check (5.2, 7.3).
7. **Long COVID profile** (5.3).
8. **Cross-profile specificity matrix** (7.4).
9. **Report extension** (Part 6) and README.

Step 5 is the gate. The CRC number tells you whether the scoring engine works at all,
and it is cheap to get wrong quietly if the other profiles are built first.
