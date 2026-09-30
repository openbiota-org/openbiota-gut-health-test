# Reading the output

Every field, and the order to read them in.

---

## Files

```
results/
├── comparison.json              side-by-side across samples, after `openbiota run-all`
└── <sample>/
    ├── <sample>_report.pdf      the Gut Health Test Metagenomic Report, ~330 pages
    ├── summary.txt              the same content as text, one section per subject
    ├── results.json             every number, machine-readable
    ├── run.log                  versions, flags, timings, reference fingerprints
    ├── depth_check.json         rarefaction, after `openbiota depth-check`
    ├── taxonomy/                species profile and cached Bowtie2 alignment
    ├── alignments/              raw and per-panel DIAMOND hits
    └── depth/                   rarefaction scratch, regenerable
```

| File | Contents |
|---|---|
| `<sample>_report.pdf` | The report, in two parts (see below). Part A puts every reading in front of the reader as a graphic in fourteen short sections; Part B breaks each one down and, under each, shows what the human research says about changing it. Every percentile uses one seven-band scale (notably low → notably high) |
| `summary.txt` | Full technical summary, one section per subject |
| `results.json` | Everything in the text report, plus more, machine-readable |
| `run.log` | Versions, flags, timings, reference set sizes, stage progress |
| `depth_check.json` | Saturation analysis, if `openbiota depth-check` has been run |
| `comparison.json` | Cross-sample table, if `openbiota run-all` was used |
| `alignments/hits_all_<mate>.tsv` | Raw DIAMOND output |
| `alignments/hits_<panel>_<mate>.tsv` | Per-panel view, with role and organism appended |
| `taxonomy/metaphlan.*.tsv` | Species profile, plus the cached alignment so re-profiling is instant |

The directory name is the sample name, which comes from the FASTQ filename —
so `SAMPLE2_A02_1.fastq.gz` produces `results/SAMPLE2_A02/SAMPLE2_A02_report.pdf`.
Renaming the input renames the output; `--sample` overrides it.

`alignments/` and `depth/` together account for most of the disk a run uses and
are fully regenerable; `openbiota prune --alignments --yes` clears them.

Subsampled runs write to `results/<sample>.subsampleN/` so smoke-test and
full-run caches do not invalidate each other.

Exit status is `1` if any `[ERROR]` diagnostic fired, so the tool composes into
a pipeline.

---

## The report's two parts

**Part A — at a glance** (sections 1–14, about 55 pages). Every reading as a
graphic, nothing deep; the full organism list sits here because "what else is
in there" is the question that follows "what is flagged".

| § | What it shows |
|---|---|
| 1 | Summary: gut-health dial, what stood out (index, missing/depleted/expanded organisms, most notable functions, top pattern), every function as a bar, community donut |
| – | Table of contents |
| 2 | What stood out, in prose: the flagged organisms, the functions, the evidence-card tally, the age estimate |
| 3 | Your gut community: composition by class, diversity and evenness as *neutral* range strips (none has a favourable direction on its own) |
| 4 | Your microbial groups at a glance: the curated groups (butyrate producers, mucin degraders, oral-origin, opportunists …) on reference bars, each with its level against the typical carrier |
| 5 | Organisms that need attention: **overgrown** (above the 90th percentile among carriers for an opportunist, above the 97th for anything, or — for an opportunist most people do not carry — several times above the level of the people who do), **missing** (commonly carried, not detected where the detection model says it should have been), **depleted**, **worth watching** (expanded, or uncommon to carry). Each with its share, the typical carrier's level, its deviation and its rank among carriers |
| 6 | Your organisms, classified: every organism detected, one share of the whole each, largest first; the thirty most abundant described one by one; the organisms furthest above and below the reference; then the full table with class, share, level against the typical carrier, how many reference adults carry it, percentile, strain |
| 7 | Pathogens: every target searched, graded (see below) |
| 8 | Estimated biological age of biota: point estimate on the training-support histogram, conformal intervals, what moved it, the model card |
| 9 | Your microbial functions at a glance: every gene-capacity pathway on the seven-band bar |
| 10 | Resemblance to published disease patterns at a glance: shape verdict (diffuse vs specific), ranked table of every scored pattern |
| 11 | Biofilm-related potential |
| 12 | Your gut mycobiome |
| 13 | Gut-skin axis |
| 14 | What you can do about it: one line per evidence card — what fired, what was studied, evidence lane (A–E, X), for/against count, status |

**Part B — in detail** (sections 15–28). Each reading broken down; under each
one an evidence card block.

| § | What it shows |
|---|---|
| 15 | Your microbial groups in detail: every group card with its members, each member's share, level and percentile, and its evidence |
| 16 | Organisms that need attention, one by one: detection power, reference carriage, level, curated interpretation with citations, disease signatures it appears in, evidence cards; then the cross-checks where two measurements disagree |
| 17 | Which strain, not just which species: the marker fingerprint of each organism's dominant population, and its comparative placement among the reference genomes of its species (nearest genome, distance, mixture evidence) |
| 18 | Pathogens in detail |
| 19 | Your microbial functions in detail: result, percentile, reference band, what it means, the research behind it, **what the research says about changing it** |
| 20 | Resemblance to published disease patterns, in detail: every scored pattern, most resemblant first, module scores, feature table, confounders, evidence cards; patterns not scored, and why |
| 21 | Biofilm findings in full |
| 22 | Your gut mycobiome, in detail |
| 23 | Gut-skin axis, in detail |
| 24 | Your information and report context: the input register — what was supplied about the person and the sample, and what was missing |
| 25 | What this report cannot tell you |
| 26 | Sequencing quality and reference coverage: the reads, the QC gates, every detection method that ran and the reference release each used |
| 27 | How this test performs: measured accuracy of the engines on simulated communities with known composition |
| 28 | Technical data |

**Organisms, on one scale.** Section 6 gives every organism one share of the
whole, on the scale of the marker catalogue that measures the composition. A
† marks a share divided among populations the whole-genome mapping tells
apart within one marker species; a ° marks an organism the marker catalogue
has no entry for, whose share is estimated from the whole-genome search and
drawn from the unclassified band; *within* marks a population counted inside
its relatives' shares above it; a § marks a share that includes what the
marker catalogue had read under a relative the whole-genome competition
rejected — the reads belong to the marked organism, so its share came with
them, and the rejected name is printed beneath. **Level** is the deviation from the typical
carrier — the median level among reference adults who carry the organism at
all, printed in brackets — with the percentile among those carriers; the two
are taken on the same distribution, so above the 50th percentile is always a
positive deviation. Where fewer than ten reference adults carry an organism no
rank is stated. Where the reference catalogue reads a species fivefold
differently from the composition, the organism is ranked instead in the
population of a lane whose reading agrees with the share (marked ‡, the
population named); where none does, the row says *not comparable*.

**Evidence cards.** Each card names the exact intervention (strain and dose as

studied, diet protocol, drug, or clinical route), who it was studied in, the
evidence *for* and *against* side by side with PMIDs/DOIs, the evidence lane
(A guideline/label · B replicated RCT · C single or conflicting RCT · D
non-randomised human · E animal/in vitro · X against), the applicability of
the studied population to this sample, the safety rules that fired and which
context fields would be needed to say more, and the phrases the report will
not use. Where no human study shows that changing a reading changes an
outcome, the card says *no supported targeted intervention* — that is a
result, not a gap. Where FMT is the only route with evidence, the card says so
and gives its regulatory status. Nothing on a card is a recommendation; the
disclaimer is on every card.

---

## Read it in this order

**1. Diagnostics.** An `[ERROR]` means the run is invalid. Zero `rpoB` makes
every figure indeterminate; a dead `butyrate` positive control means something
is broken upstream. Stop there.

**2. The `rpoB` denominator.** Everything is divided by it. The summary reports
its yield per million read pairs and flags implausible values. An
under-detected denominator inflates every pathway at once.

**3. Stability flags.** `UNSTABLE` means the fragment count is low enough that
counting noise dominates, however precise the number looks.

**4. The numbers.** Per-gene values are more interpretable than the panel
aggregate. Where a residue-consistent subset is reported, the headline is the
upper bound and the subset is the tighter estimate.

**5. The cross-pathway pattern.** Everything elevated at once points at the
denominator, not at biology.

---

## Key fields

### Copies per 100 bacterial genomes

The primary measurement. For a single-copy gene: 100 means essentially every
bacterium carries it, 10 means roughly one in ten, 1 means roughly one in a
hundred.

### Percentile

Position within the reference cohort. The 60th percentile means 60% of cohort
samples had a lower value. Requires `refs/reference_ranges.json` from
`openbiota cohort`.

### Bands

`trace` / `low` / `moderate` / `high` / `very high` are scan aids computed from
the absolute value. The **percentile** is the meaningful comparison; the band is
just a label.

### Flags

| Flag | Meaning |
|---|---|
| `ok` | Nothing to note |
| `UNSTABLE` | Below the panel's minimum fragment count; counting noise dominates |
| `LOW-ID` | Fragments average under 75% translated identity — matching distant relatives rather than the gene itself. Upper bound only. |
| `few-refs` | Fewer than 25 reference sequences; a low value is weak evidence of absence |
| `refs-capped` | Reference set hit its size cap; sensitivity is bounded by it |
| `upper-bound` | Panel has an active-site filter that reads are mostly too short to apply |
| `INDETERMINATE` | No `rpoB` denominator |
| `extension` | Pathway beyond the original scope |
| `none` | Nothing detected |

### Diagnostic codes

| Code | Level | Meaning |
|---|---|---|
| `rpob-zero` | error | No denominator; every figure indeterminate |
| `positive-control-zero` | error | `butyrate` returned nothing; the run is invalid |
| `rpob-low-yield` / `rpob-high-yield` | warn | Denominator outside the plausible band |
| `all-panels-elevated` | warn | Everything high at once — a normalisation artefact |
| `low-identity-<gene>` | warn | Distant-homolog matching; see `LOW-ID` |
| `implausible-<panel>` | warn | More gene copies than there are genomes to carry them |
| `residue-concordant-<gene>` | info | The two independent residue checks agree |
| `residue-discordant-<gene>` | warn | They disagree; trust neither over the headline |
| `unstable-<panel>` | info | Below the stability threshold |
| `refs-small-<gene>` / `refs-truncated-<gene>` | info | Reference set size notes |
| `unresolved-subjects` | warn | Cached hits came from a different database; re-run with `--force` |

### Residue check

Reported in two halves for `urdA`.

*Reference side* — how many accepted fragments matched a reference that itself
carries an acceptable active-site residue. Applies to every fragment.

*Read side* — of the fragments whose alignment reached the site, how many
carried an acceptable residue. This is the published filter, applied where
available. Typically a small percentage of fragments at 151 bp.

When both are available the report cross-checks them against each other, since
they share no evidence.

### Closest reference organisms

The reference sequences the reads matched best. This is not taxonomic
identification — a read matching a *Clostridium* reference means the nearest
sequence in the reference set came from one, nothing more.

### Confidence tier

`confirmed` means the reading meets the criteria at which validation produced
no false positives: at least 2 matching fragments and at least 75% mean
translated identity. `provisional` means detected but below that bar — a
possible signal rather than a reliable one. See
[VALIDATION.md](VALIDATION.md#operating-points-and-why-the-thresholds-are-what-they-are).

### What it means

Each pathway carries an explanation matching where the reading fell — above the
75th percentile, below the 25th, or in between. The text is specific to that
metabolite and comes from `interpretation.implications` in its panel file, so
the report and the panel definition cannot drift apart.

### Community composition

Derived from the `rpoB` fragments, which are already computed as the
denominator. Phylum-level proportions are reliable. Genus labels are
nearest-reference matches. Proportions approximate cell fractions, not the
DNA-mass fractions a read-count profiler like Kraken reports.

---

## Profile similarity fields

| Field | Meaning |
|---|---|
| `modules.<name>.score` | Weighted mean of feature values, −1 to +1 |
| `modules.<name>.percentile` | Rank of that score against the reference null |
| `modules.<name>.weight_coverage` | Fraction of planned feature weight actually measured |
| `combined.percentile` | The headline figure. Fraction of the matched reference less concordant with the pattern |
| `confidence.grade` | HIGH / MODERATE / LOW / INSUFFICIENT |
| `confidence.reasons_for_downgrade` | Why it is not higher, in plain language |
| `abstention.abstained` | Whether the engine refused to score |
| `abstention.what_would_change_this` | What to fix, if it did |
| `cross_engine_checks[].verdict` | CONCORDANT / DISCORDANT / UNAVAILABLE |
| `dysbiosis_anchor.score` | GMWI2. Positive leans healthy, negative leans disturbed |
| `features[].v` | That feature's signed contribution, −1 to +1 |
| `features[].w/q/s/c` | Evidence weight, cohort independence, specificity, direction agreement |
| `features[].reference_percentile` | Where the sample sits for that feature |

A `combined.percentile` of 63 means the sample is more concordant with the
profile's prespecified pattern than 63% of the matched reference set. It is
**not** a 63% chance of having the condition.

Added by the v3 contract, on each profile result:

| Field | Meaning |
|---|---|
| `status` | `scored`, `low_coverage`, `not_computable` or `abstained` |
| `confidence_vector.*` | The dimensions behind the grade, kept separate rather than multiplied: `technical_reliability`, `independent_feature_coverage`, `evidence_maturity`, `assay_transportability`, `population_transportability`, `empirical_specificity`, `literature_uniqueness_prior`, `out_of_distribution`, `direction_conflict`, `interval_width` |
| `confidence_vector.uncertainty_interval` | Bootstrap interval on the percentile |
| `competing_profiles` | Other profile families that score within reach of this one |
| `top_two_margin` | Percentile points between this profile and the next-highest family |
| `modules.<name>.module_status` / `independent_cluster_coverage` | Whether a study module could be scored, and what fraction of its independent feature clusters were measured |
| `modules.<name>.bound` / `would_bind_if` | Whether the module's evidence type has a measurement in this pipeline, and what would bind it if not |
| `features[].evidence_type` / `claim_level` | What kind of measurement the feature is and what the source study could legitimately claim |
| `features[].scoring_weight` / `audit_factor_wqsc` | Only `w` enters the score; `w·q·s·c` is recorded for audit |

And alongside the profiles, under `profile_similarity`:

| Field | Meaning |
|---|---|
| `ranked[]` | Every profile with status, percentile, concordance and maturity, highest first |
| `shape_analysis.verdict` | `specific`, `diffuse`, `quiet`, `mixed` or `unscored` — the library-wide reading |
| `shape_analysis.n_above_typical_band` / `n_notably_high` | How many profiles sit above the 75th / 95th percentile |
| `shape_analysis.mean_specificity_of_high_scorers` | Whether the high scorers discriminate their own condition or track general disturbance |
| `confounder_ledger.recorded` / `missing` / `profiles_abstained` | Which subject facts were supplied, which were not, and which profiles refused to score for want of them |
| `withheld_opt_in_profiles` | Profiles skipped without `--include-opt-in` |
| `community.groups[]` | Ten curated microbial groups: summed abundance, cohort percentile, prevalence, each member's own placement |
| `community.species[]` | Every species named by the scoring lane: abundance, cohort percentile, prevalence, and which groups it belongs to |

Full detail on the scoring: [PROFILES.md](PROFILES.md).

## Sequencing quality fields

`sequencing_quality.gates` is the spec 4.2 gate table.

| Field | Meaning |
|---|---|
| `overall` | `pass`, `warn`, `fail` or `unknown` — the worst individual gate |
| `usable_nonhost_pairs` / `depth_adequate` / `minimum_usable_pairs` | The prespecified depth floor (500,000 pairs) and whether it was met |
| `gates[]` | One row per gate: `gate`, `label`, `value`, `display`, `status`, `note`. A gate that could not be measured is `unknown`, never `pass` |
| `provenance` | Instrument, run, flowcell, lane and index read from the FASTQ headers, and the batch key they form |

`sequencing_quality.fastp` is the parsed fastp report: pairs and bases in and
out, Q20/Q30, GC, duplication rate, read-length profile, adapter-trimmed reads,
insert-size peak, and whether the reads arrived pre-trimmed.

## Extended catalogue fields

`extended_catalogue` is the MetaPhlAn 4 SGB inventory. It is reported, never
scored — `what_this_is` says so in the file.

| Field | Meaning |
|---|---|
| `n_sgbs_detected` / `n_species_detected` / `n_genera_detected` | Catalogue breadth in this sample |
| `n_unnamed_sgbs` / `unnamed_sgb_percent` | Genome bins with no cultured representative, and their share of the sample |
| `unclassified_percent` | What even this catalogue could not name |
| `sgbs[]` | Every bin: `sgb`, `species`, `genus`, `kingdom`, `percent`, `unnamed` |
| `bridge_to_scoring_lane` | Species named by both profilers, by MetaPhlAn 4 alone, and by MetaPhlAn 3 alone, with examples |

## Version identifiers

Two profile results are comparable only when all five of these match. They are
all in `results.json` and on the technical page of the PDF.

| Identifier | Where |
|---|---|
| Pipeline version | `version` |
| DIAMOND reference fingerprint | `run["reference fingerprint"]` |
| MetaPhlAn database | `profile_similarity.taxonomic_engine.database` |
| Reference cohort manifest | `profile_similarity.reference.taxonomic.manifest.manifest_id` |
| Profile version | `profile_similarity.profiles.<name>.profile.version` |

## results.json structure

```
tool, version, sample
measures, does_not_measure
run                     flags, timings, versions, fingerprints
input                   read counts, lengths, GC, quality, mate structure
normalisation           rpoB counts, the formula, indeterminate flag
search                  hit lines, mate collapse accounting
panels[]                per pathway:
  accepted_fragments, decoy_fragments, copies_per_100_genomes
  copies_per_100_genomes_residue_consistent
  aggregate_method, aggregate_from, stable, caveats
  genes[]               per gene: fragments, identity, references, residue_check
  decoys[]
reference_comparison    percentile, cohort quartiles, status per pathway
community_profile       phyla, genera, diversity indices
diagnostics[]           level, code, message
reference_entries       provenance for every reference set
profile_similarity      disease-pattern engine (see below)
microbiome_age          status, estimate, 80%/95% intervals, OOD detail,
                        contributions, model bundle id
findings_and_evidence
  findings              reference context, detection model, every species with
                        detection status/power, reference prevalence class,
                        percentile among carriers with 90% bootstrap interval,
                        BH q, bucket 1-7, curated interpretation; contradictions
  triggers[]            what fired (panel/group/profile/species, direction)
  evidence_summaries[]  per card: intervention identity, studied protocol,
                        supporting/against assertions with sources, lane,
                        applicability, safety verdict, display policy,
                        never_say, follow-up route, FMT status
  registry              digests of the intervention/safety/interpretation registries
```

---

## Comparing over time

The same sample analysed twice gives the same answer; two samples from the same
person weeks apart give a trend. That trend is more informative than any single
value, because the method's systematic biases — reference set composition,
identity floors, read length — are held constant across runs and largely cancel
in a difference.

Keep results directories per timepoint and compare
`reference_comparison.panels.<name>.value` across them.

## Strain resolution (section 17)

Species profiling says an organism is present. Strain resolution says *which
one* — and for much of gut biology that is the question that matters: whether
an *E. coli* is a harmless resident or carries an intact Shiga-toxin operon is
a strain-level question, not a species-level one.

The report resolves, for every organism with enough coverage, the **dominant
population's consensus** across marker genes, at single-base resolution. In
the supplied samples that is 49–78 organisms each, over 5.5–7.5 million
callable bases.

### What a fingerprint is, and is not

| It is | It is not |
|---|---|
| This sample's own population of that organism | A match to a named laboratory strain |
| The dominant population | Proof of a single strain — a minority population can be concealed |
| Comparable between specimens at base level | A whole-genome identity claim from conserved markers |
| Evidence of an organism's identity | Attribution of a separately detected toxin or resistance gene |

An organism detected but not resolved is reported as exactly that. Absence of
resolution is a limit of sequencing depth, not a finding about the organism.

### How it is produced

The marker lane is a second pass over the reads, pinned to MetaPhlAn **4.1.1**
and the `mpa_vJun23_CHOCOPhlAnSGB_202403` database — the same pair that
produced the sample's abundance profile. Mixing the 4.1 and 4.2 workflows is
not supported and is actively prevented: a 4.2 `strainphlan` on `PATH` will not
be used against 4.1 alignments.

```bash
# emit SAM alongside the usual outputs, in a fresh directory so an existing
# bowtie2 summary cannot cause the alignment to be skipped
metaphlan reads_1.fastq.gz,reads_2.fastq.gz --input_type fastq \
  --bowtie2db refs/metaphlan4_db --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
  --nproc 32 -s sample.markers.sam.bz2 \
  --bowtie2out sample.markers.bowtie2.bz2 -o sample.markers.species.tsv

sample2markers.py -i sample.markers.sam.bz2 -o consensus_markers \
  -d refs/metaphlan4_db/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl -n 32
```

`samtools` is a hard dependency of `sample2markers.py`.

## What was asked for, and what came back

Every report now carries a **completeness census**: each registered resolution
target with the state it actually reached. The distinction it preserves is the
one that matters most in a screening report:

> An unrun, failed or unavailable assay is never reported as a negative result.

States include `completed`, `reference_unavailable` (the tool or reference is
not installed), `incompatible_input` (an RNA virus cannot be assayed from DNA
sequencing, and adding reference sequences cannot change that) and `scheduled`
(registered, with reference curation outstanding). The full per-target detail,
with a reason for every non-terminal target, is in `results.json` under
`resolution_census`.
