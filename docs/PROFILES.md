# Disease similarity profiles

A profile is a declarative description of a published, group-level
disease-associated microbiome pattern. It says which features were reported,
which way each moved, how the finding was measured, and how far the evidence
has come. It contains no scoring logic — that lives in `openbiota/scoring.py`,
which has no knowledge of which condition it is scoring.

Thirty-three profiles ship. List them with `openbiota profiles`, inspect one
with `openbiota profiles --show crohns`. Every profile in the library is scored
on every sample, and the result is read as a *shape* across the library rather
than as thirty-three independent verdicts.

## What a percentile means

> The sample is more concordant with this profile's prespecified pattern than
> that fraction of the matched reference set.

It is **not** a probability of having the condition. Every design decision
below follows from holding that distinction.

## The library

| Profile | Family | Maturity | Study modules (fused) | Features | Evidence | Challenge for |
|---|---|---|---|---|---|---|
| `crc` | colorectal | P3 | 5 (2) | 22 | taxa · genes · carriers | adenoma, crohns, cirrhosis |
| `adenoma` | colorectal | P3 | 3 (2) | 11 | taxa · genes | crc |
| `crohns` | ibd | P3 | 7 (5) | 41 | taxa · genes · carriers · ecology | crc, uc |
| `uc` | ibd | P3 | 6 (4) | 19 | taxa · genes · carriers · ecology | crohns, crc |
| `ibd` | ibd | P3 | 6 (4) | 27 | taxa · genes · carriers · ecology | — |
| `ibs`, `ibs_c`, `ibs_d`, `ibs_m` | ibs | P2 / P1 | 5 (4) / 3 (2) / 3 (2) / 2 (1) | 19 / 4 / 11 / 6 | taxa · genes (· ecology) | — |
| `t2d` | cardiometabolic | P3 | 5 (3) | 23 | taxa · genes · carriers | — |
| `acvd` | cardiometabolic | P3 | 4 (2) | 22 | taxa · genes · carriers | cirrhosis, t2d |
| `hypertension` | cardiometabolic | P2 | 3 (2) | 13 | taxa · genes · ecology | acvd, obesity, t2d |
| `ckd` | cardiometabolic | P2 | 4 (2) | 21 | taxa · genes · carriers | t2d, hypertension, cirrhosis |
| `obesity` | cardiometabolic | P3 | 6 (3) | 9 | taxa · genes · carriers · ecology | — |
| `cirrhosis` | liver | P3 | 5 (3) | 23 | taxa · genes · carriers · ecology | crc, crohns |
| `masld` | liver | P2 | 6 (4) | 16 | taxa · genes · carriers · ecology | cirrhosis, obesity, t2d |
| `parkinsons` | neurodegeneration | P3 | 5 (3) | 30 | taxa · genes · carriers | — |
| `ad_clinical`, `ad_mci`, `ad_preclinical_amyloid` | alzheimers_continuum | P2 | 4 (1) / 3 (1) / 2 (1) | 20 / 10 / 5 | taxa · genes (· pathways, unbound) | — |
| `ms` | autoimmune | P3 | 5 (3) | 16 | taxa · genes · carriers | — |
| `ra` | autoimmune | P3 | 4 (2) | 11 | taxa · genes (· strain, unbound) | — |
| `t1d` | autoimmune | P2 | 5 (3) | 13 | taxa · genes · carriers · ecology | — |
| `celiac` | autoimmune | P2 | 4 (2) | 14 | taxa · genes (· CAZy, unbound) | — |
| `ankylosing_spondylitis` | autoimmune | P2 | 4 (2) | 16 | taxa · genes · carriers | — |
| `mdd` | psychiatric | P2 | 4 (3) | 16 | taxa · genes | — |
| `alopecia_areata` | alopecia_areata | P1 | 9 (5) | 29 | taxa · genes (· pathways, unbound) | — |
| `androgenetic_alopecia` | androgenetic_alopecia | P1 | 4 (3) | 11 | taxa · genes | — (sex-stratified) |
| `mecfs` | post_infectious | P2 | 3 (3) | 7 | taxa · genes · ecology | longcovid, ibs |
| `longcovid` | post_infectious | P1 | 4 (3) | 15 | taxa · genes · ecology | mecfs, ibs |

"Study modules" counts every module in the profile; "fused" is how many enter
the combined score (the rest are sensitivity panels, phenotype panels, or
modules whose evidence type has no engine yet). Measured accuracy where a
labelled cohort exists: [VALIDATION.md](VALIDATION.md).

## The contract (spec 3)

Every profile declares what kind of claim it is making, and the loader refuses
one that does not.

**Evidence type** — what was measured. Bound types have an engine in this
pipeline; unbound types load, are shown, and say what would bind them.

| Type | Claim level | Engine | Status |
|---|---|---|---|
| `TAX_REL`, `TAX_PRES` | observed | MetaPhlAn 3 species relative abundance / presence | bound |
| `GENE_ABUND` | genomic_capacity | DIAMOND gene panels, rpoB-normalised | bound |
| `CARRIER_ABUNDANCE` | carrier_proxy | sum of a curated carrier set (panel `organisms` or a taxon group) | bound |
| `ECO` | observed | Shannon, richness, evenness | bound |
| `STRAIN`, `PANCNV` | observed | StrainPhlAn / MIDAS2 / PanPhlAn | unbound |
| `SEQ_HMM`, `PATH_ABUND`, `PATH_COMPL`, `CAZY`, `AMR` | genomic_capacity | HMMER / HUMAnN 3 / assembly + dbCAN3 / CARD | unbound |
| `FLUX`, `METAB`, `MODEL` | flux_predicted / external_measured / model_output | AGORA2 / metabolomics file / frozen container | unbound |

**Evidence maturity** — how far the research has come, P0 to P5:

| | |
|---|---|
| P0 | animal, case report, or mechanistic hypothesis only |
| P1 | one small human cohort, amplicon-only, or inaccessible individual data |
| P2 | usable human shotgun cohort, or multiple compatible human association cohorts |
| P3 | independent human replication with coherent task and assay |
| P4 | external population validation with a frozen model and proper participant split |
| P5 | prospective intended-use validation with calibrated clinical performance |

Nothing in the library is above P3. The report prints the grade beside every
score.

**Assay transport** — whether the source study measured what this pipeline
measures. `native` means shotgun species or genes. `amplicon_taxon_to_shotgun`
means a 16S genus-level finding carried onto shotgun species, and is the only
way a genus-level feature is accepted (below). `predicted_to_measured` means a
PICRUSt inference tested against genes actually counted.

**Study modules and fusion.** A v3 profile is a set of modules, one per
independent published cohort, each carrying a `study_group`. Modules that share
participants share one group. Fusion is `study_average`: independent groups are
averaged with equal weight, so a study contributing twenty correlated taxa
cannot outweigh a study contributing one replicated feature. `sensitivity_panel`
and `phenotype_panel` modules are scored and shown but never fused. The three
legacy profiles (`crc`, `mecfs`, `longcovid`) keep their taxonomic /
functional / ecological split with the documented 50/35/15 weights.

## Scoring

**1 — Per-feature robust value**, against the matched reference only:

```
z_i = (CLR_i(sample) − median_i(ref)) / (1.4826 × MAD_i(ref))
v_i = d_i × clip(z_i, −3, +3) / 3
```

`d_i` is +1 for a reported increase and −1 for a decrease, so `v_i` runs from
−1 (opposite) to +1 (strongly concordant). Zero-inflated features (carried by
under 80% of the reference) use a prevalence-aware rank instead: absence is a
category, not a small number, because the CLR of an absent taxon carries only
the denominator's noise. Median and MAD, not mean and SD, because abundance
distributions are skewed enough that one outlier otherwise sets the scale.
Centred log-ratio, because relative abundances are compositional.

Genus-level and carrier features are the log-ratio of a *sum* of species
against the same prevalent-core frame, computed identically for the sample and
for every reference sample.

**2 — Module aggregation.**

```
module = Σ(wᵢ vᵢ) / Σ(wᵢ)
```

`w` is the evidence weight. The three other factors every feature carries —
`q` cohort independence, `s` disease specificity, `c` cross-study direction
agreement — are **recorded for audit and printed in the decomposition table,
not multiplied into the score**. Multiplying them made every profile's number
depend on four curator judgements at once; the v3 contract keeps them visible
and keeps the score simple.

**3 — Percentile.** Every reference sample is scored through the identical
path and the sample is ranked against that null. Never by rescaling the raw
value. A bootstrap over the reference gives a likely range (the report's
"likely 71–100").

**4 — Fusion**, as declared by the profile (`study_average` or the legacy
weighted split).

## Two rules the loader enforces

**Genus level is conditional.** Aggregating a genus reverses direction in
practice: genus *Faecalibacterium* has been reported higher in a cohort where
*F. prausnitzii* was lower. So a `genus` feature is accepted only when the
module declares `assay_transport: amplicon_taxon_to_shotgun` — the source could
not have named the species — and, for the conflicted genera
(*Faecalibacterium*, *Streptococcus*, *Veillonella*, *Prevotella*, *Blautia*,
*Lachnospiraceae*, *Alistipes*, *Bacteroides*, *Dorea*), only when it also
declares `direction_conflict: material`. Shotgun-derived genus evidence is
never accepted: the study could have named the species and did not.

**Every feature carries all four factors and a source.** A feature without
them cannot be audited, so it is not accepted.

## Confidence is a vector, not a multiplier

Beside every score the report prints, and never folds into the number:

| Dimension | Answers |
|---|---|
| technical reliability | depth, host fraction, reference match, profiler-version agreement |
| independent feature coverage | share of the pattern's independent evidence clusters this sample could measure |
| evidence maturity | the P-grade |
| assay transportability | native / transported |
| population match | how well the matched reference resembles the subject |
| specificity | empirical (challenge-cohort AUROC, where measured) or the literature-uniqueness prior |
| out-of-distribution | whether the sample sits outside the reference for the features used |
| direction conflict | none / documented / material |
| interval width | the bootstrap range |

A weak row means "less certain", not "lower".

## Shape of the result (spec 7.8)

Thirty percentiles are not thirty findings. After every profile is scored the
library is read as a whole:

- **quiet** — nothing above the 75th percentile (the report's typical band)
- **specific** — one or two profiles above it, with distinctive features
- **diffuse** — a large fraction above it, resting on shared features
  (butyrate-producer loss, oral-taxa gain); the general dysbiosis index is the
  headline, not any disease
- **mixed** — several above it with partial specificity

Each scored profile also lists its **competing** profiles (the next-highest
from other families) and the **top-two margin**, so a 90th-percentile result
two points ahead of its neighbour reads as what it is.

**Specificity** is empirical where the challenge-cohort matrix in
`docs/profile_validation.json` covers the profile (worst off-target AUROC
against its own), and otherwise the labelled *literature-uniqueness prior*:
the mean over the profile's features of 1 / (number of profiles sharing that
feature). The report always says which one it used.

## The dysbiosis anchor

GMWI2 runs on every sample and is reported beside every profile score. If the
anchor is strongly negative the community is generally disturbed and a high
score on any profile probably reflects that. The 95 published coefficients are
stored verbatim in `openbiota/data/gmwi2_model.json`; `tests/test_dysbiosis.py`
reproduces the authors' worked example to four decimal places.

## Abstention and the confounder ledger (spec 12)

A refused score is more useful than a confident wrong one. Profiles abstain
when:

- broad-spectrum antibiotics fall inside the declared window (90 days)
- usable non-host pairs fall below 500,000
- required metadata is missing — `t2d` needs metformin status, `celiac` needs
  gluten-free-diet status, `androgenetic_alopecia` needs sex
- a `mandatory_match` field (age band, sex) cannot be matched in the reference

An abstaining profile renders the reasons and what would change them. Beside
the library the report prints a **confounder ledger**: what was known about the
subject when scoring (antibiotics, metformin, PPI, gluten-free diet, country,
age, sex, stool form, hospitalisation, kidney function, medications,
sequencing batch from the FASTQ headers), what was missing and what it would
have changed, and which profiles abstained because of it.

**Opt-in profiles.** The mechanism exists — a profile marked `opt_in: true`
is scored only with `--include-opt-in` — but no profile in the shipped library
uses it. `mdd` did until the transfer evidence was reviewed: patient stool has
induced depression-like behaviour in germ-free mice and in microbiota-depleted
rats, so the pattern is measured and reported like every other one, under the
same caveat that a percentile is resemblance and not a diagnosis.

## Reference contrasts

Profiles need two contrasts: people with the same exposure who did *not*
develop the syndrome (primary), and general healthy adults matched on the
non-exposure variables (secondary). For ME/CFS and Long COVID the primary
contrast **does not exist** in the published cohorts; both profiles declare
`primary: source: null`, say so in the report, and take a confidence penalty.

## Adding a profile

Drop a YAML file in `profiles/`. The loader checks:

- `disease`, `task`, `body_site`, `population`, `evidence_maturity`
- semantic `version` — changing any weight bumps it
- every module has an `evidence_type`, a `study_group`, an `assay_transport`
  and a `module_type`
- every feature has `direction`, `w q s c` in [0, 1], a `level`, a `source`
- genus features obey the transport rule above
- `stratify_on` names a subject field; `mandatory_match` fields exist
- `challenge_for` names profiles that exist
- cross-engine checks name features that exist

`profiles/crohns.yaml` is the fullest v3 example (seven study modules,
carrier features, an amplicon-transported module, challenge declarations);
`profiles/mecfs.yaml` the fullest legacy one (duration strata, cross-engine
check, absent primary contrast).

## Sources

Every profile carries its citations inline, per feature and at the profile
level; `openbiota profiles --show <name> --json` prints them all. Method sources:

- Compositional data: Gloor et al., *Front Microbiol* 2017;8:2224,
  [doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224)
- GMWI2: Chang & Gupta et al., *Nat Commun* 2024;15:7447,
  [doi:10.1038/s41467-024-51651-9](https://doi.org/10.1038/s41467-024-51651-9)
- MetaPhlAn: Blanco-Míguez et al., *Nat Biotechnol* 2023,
  [doi:10.1038/s41587-023-01688-w](https://doi.org/10.1038/s41587-023-01688-w)
- curatedMetagenomicData: Pasolli et al., *Nat Methods* 2017;14:1023–1024,
  [doi:10.1038/nmeth.4468](https://doi.org/10.1038/nmeth.4468)
- Colorectal cancer: Wirbel et al., *Nat Med* 2019;25:679–689; Thomas et al.,
  *Nat Med* 2019;25:667–678
- ME/CFS: Guo et al. and Xiong et al., *Cell Host Microbe* 2023;31
- Long COVID: Liu et al., *Gut* 2022;71:544–552; Su et al., *Cell Host
  Microbe* 2024;32:651–660

Measured performance and limitations: [VALIDATION.md](VALIDATION.md).
