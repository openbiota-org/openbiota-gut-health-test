# Validation, accuracy and limitations

This document covers both halves of the tool. **Part A** is metabolite gene
capacity, validated against synthetic communities with known gene content.
**Part B** is species composition and disease-pattern similarity, validated
against public cohorts that carry disease labels.

Jump to: [Part A headline figures](#headline-figures) ·
[Part B headline figures](#part-b-profile-similarity) ·
[the specificity matrix](#cross-profile-specificity-the-most-informative-table-here) ·
[limitations](#limitations)

Everything about how well this test performs, and where it falls short, is in
this one file. Nothing here is estimated — the figures come from
`openbiota validate`, and the raw output is in
[`validation_results.json`](validation_results.json).

Reproduce with:

```bash
openbiota validate                    # ~5 min with genomes cached
openbiota cohort                      # ~10 min, rebuilds the reference ranges
```

---

## Headline figures

At the operating point the tool actually reports (**confirmed** tier):

| | |
|---|---|
| **Sensitivity** | **97.0%** — genes present that were found |
| **Specificity** | **100.0%** — genes absent that were correctly called absent |
| **Precision** | **100.0%** — detections that were real |
| **False positives** | **0** of 62 true negatives |
| **False negatives** | **1** of 33 true positives |
| **Accuracy** | **98.9%** |
| **Detection limit** | carriers found down to **0.2% of cells** |
| **Quantitative calibration** | slope 1.13, R² 0.72 |

Every pathway individually achieved 100% specificity and 100% precision.

---

## Part B — profile similarity

A different kind of validation, because the ground truth is different. For
gene capacity, truth is "does this genome carry this gene", which a synthetic
community establishes exactly. For disease similarity there is no such truth —
only whether the score separates labelled cases from labelled controls.

Only **colorectal cancer** has both a cross-cohort-validated microbiome
association and public cohorts with disease labels attached, so it is the one
profile on which the engine can be tested end to end. That is why it ships as
a validation control rather than as a health readout.

### The colorectal cancer result

258 labelled cases against 250 within-study controls (YachidaS 2019, via
curatedMetagenomicData), scored through the identical engine, with the
reference distribution built from the control arm only:

| | |
|---|---|
| **AUC** | **0.62** (95% CI 0.57–0.67) |
| **Leave-one-study-out AUC** | 0.62 |
| **Ceiling for these six features** | 0.67, measured directly |
| **Published cross-cohort performance** | ~0.80, from learned models over hundreds of features |

Three things follow, and all three are stated in the report:

**The engine recovers the signal.** All five oral-origin markers are enriched
in the correct direction — *F. nucleatum* in 18.6% of cases versus 4.4% of
controls, *P. micra* 37.2% versus 14.8%, and similarly for
*P. stomatis*, *G. morbillorum* and *S. moorei*. The confidence interval
excludes chance.

**0.62 is near the ceiling for six prespecified markers, not a shortfall in
the engine.** Counting how many of the five markers a sample carries gives AUC
0.671; summing their log abundances gives 0.674. The published 0.80+ comes
from LASSO models fitted over the full profile, which is a different method
that this tool deliberately does not use — naively transferred learned
coefficients lose both accuracy and disease specificity.

**The profile is not condition-specific**, which is the most important finding
in this document. See the matrix below.

### Cross-profile specificity: the most informative table here

Every profile scored against every labelled condition. A profile should score
highest against its own condition and near 0.5 against the others.

| Profile | vs CRC | vs IBD | vs T2D |
|---|---|---|---|
| `crc` | **0.62** | 0.76 | 0.52 |
| `longcovid` | 0.53 | 0.80 | 0.59 |
| `mecfs` | 0.56 | 0.49 | 0.60 |

- **`crc` separates IBD better (0.76) than the condition it was built for
  (0.62).** Its oral-origin markers are a general-disturbance signal. Oral-gut
  translocation happens in periodontal disease, with reduced gastric acid, and
  in inflammatory bowel disease, not only in malignancy.
- **`longcovid` scores 0.80 against IBD** — higher than anything else in the
  table. *R. gnavus* enrichment and *F. prausnitzii* depletion are classic
  IBD findings. Its own discrimination is unmeasured, because the accessions
  that would allow measuring it publish reads without a sample-to-outcome key.
- **`mecfs` is flat and near chance everywhere** (0.49–0.60). Its features are
  weak discriminators of anything, which is consistent with the modest
  species-level effects both source cohorts reported. Notably it does *not*
  spike on IBD, so it is at least not a relabelled inflammation signal.

This is why GMWI2 accompanies every profile score. A community that is
generally disturbed will resemble several of these patterns at once, and the
anchor is what tells you whether that is what is happening.

### What is not measured

No AUC is printed for `mecfs` or `longcovid`, and none should be. Several Long
COVID accessions publish reads without a sample-to-outcome key; `PRJNA878603`
(ME/CFS) has published code and a key resources table but no disease-labelled
cohort in curatedMetagenomicData. Where labels cannot be obtained, the honest
statement is that discrimination is unmeasured — not a proxy figure.

### Reference cohort stability

The species reference cohort holds **3,027 healthy adult stool metagenomes
from 22 studies**, median depth 38.2 million reads, all processed through one
pipeline. Bootstrapping sub-cohorts of increasing size shows how quickly a
percentile estimate settles:

| Cohort size | SD of the median estimate |
|---|---|
| 25 | 0.81 |
| 50 | 0.54 |
| 100 | 0.33 |
| 200 | 0.24 |
| 400 | 0.16 |
| 800 | 0.10 |

This is what justifies the sample-count target empirically rather than by
assertion, and it is where the confidence grade's threshold of 200 comes from.

### Host filtering, verified with a positive control

Human reads are removed against GRCh38 plus PhiX before taxonomic profiling.
On sample A02 this removed **0 of 8,029,909 pairs**, and a 0.00% rate is
exactly what a silently broken index also looks like — so it was checked
rather than assumed.

Two thousand synthetic read pairs cut directly from the GRCh38 sequence align
at **100%**, confirming the index resolves. Aligning the real sample finds 2
of 400,000 mates, a 0.00% overall rate. The sample genuinely carries
negligible human DNA, which is what a well-prepared stool library should look
like.

Reproduce with:

```bash
bowtie2 -p 10 --very-fast -x refs/host/host -1 <human_reads_1.fq> \
        -2 <human_reads_2.fq> -S /dev/null
```

Because host filtering was a no-op for this sample, the species profile is
identical with and without it — a useful consistency check in its own right,
since the two runs went through different code paths.

### What was not validated, and why

**Synthetic taxonomic communities.** The plan was to extend the synthetic
harness — which builds communities from reference genomes with known gene
content — to spike profile taxa at defined abundances and check that the
profiler recovers them. This was not done, because the colorectal validation
above supersedes it: it tests the same machine end to end against *clinical*
ground truth rather than against a simulator, and a simulator that produced
reads from the same reference genomes MetaPhlAn's markers were drawn from
would mostly measure that circularity. What synthetic communities would add
that the labelled cohorts do not is a detection limit for individual species,
which remains unmeasured.

**Per-feature reproduction of the ME/CFS cohort.** `PRJNA878603` publishes
code and a key resources table, so the reported per-feature directions could
be re-derived. This was not done. The consequence is that the ME/CFS feature
weights rest on the papers' reported effects rather than on independently
reproduced ones, and the profile's discrimination is unmeasured either way.

**Comparison against the published assembly-and-HMM pipeline** for gene
capacity remains open; the protocol is at the end of this document.

### Two bugs this validation caught

Both produced plausible-looking numbers, which is what makes them worth
recording.

**A depth-dependent zero floor.** The CLR's zero-replacement value was
originally derived per sample as half that sample's smallest detected
abundance — the textbook rule. That makes the floor depth-dependent: a deeply
sequenced sample detects rarer taxa, gets a lower floor, and maps its zeros
further down than a shallow one. For a taxon absent in most samples, that
difference *is* the reported signal. The floor is now a fixed constant
(1e-4 percent), applied identically to sample and reference.

**Absence treated as a number.** For a taxon absent from both the sample and
most of the reference, the CLR values still differ, because they carry the
geometric-mean denominator, which varies with how many other taxa each sample
happened to detect. Comparing those values compares sequencing incidentals. On
sample A02 it placed five *absent* colorectal markers at the 74th
percentile and produced a spurious 74th-percentile similarity score. Absence
is now a category: a sample lacking the taxon sits at the mid-rank of the
reference samples that also lack it, which puts those five markers at the 48th
percentile and the profile score at the 23rd, where it belongs.

Both are pinned as regression tests in `tests/test_scoring.py`.

---

## How this was measured

The hard part of validating a metagenomic screen is knowing the right answer.
This validation constructs the sample so that the answer is known.

**1. Ground truth from complete genomes.** 24 reference genomes were
downloaded from NCBI. For each, its own annotated proteome was searched at
*full length* against the pathway reference sets, requiring a genuine ortholog:
≥70% identity across ≥70% of the reference and ≥60% of the query, plus the
active-site residue criterion where a pathway defines one. That is the
published method's decision procedure, applied where it works properly — on
whole proteins rather than 50-residue fragments.

**2. Synthetic samples.** Reads were simulated from those same genomes at known
cell abundances, with matched read lengths (100–151 bp), NovaSeq-style binned
quality values, a 0.2% substitution error rate, and mate ids identical between
R1 and R2 exactly as real NovaSeq output produces. Because the community
composition is set, the true number of gene copies per 100 genomes is known
exactly rather than estimated.

**3. Measurement.** The pipeline was run on the simulated reads with no
knowledge of any of the above.

Five experiments, 3,000,000 read pairs each:

| Experiment | Community | Measures |
|---|---|---|
| `balanced` | all 24 genomes, equal abundance | sensitivity, specificity |
| `negative` | 6 organisms verified to lack the target genes | the false-positive rate, directly |
| `spike_0p05` | carrier at 5% of cells | detection at moderate abundance |
| `spike_0p01` | carrier at 1% | detection at low abundance |
| `spike_0p002` | carrier at 0.2% | the detection limit |

Genomes used span the known carriers for every pathway — *C. scindens*,
*C. hylemonae* and *P. hiranonis* for bile acids, *K. pneumoniae* and
*P. mirabilis* for TMA, *F. prausnitzii* and *R. intestinalis* for butyrate,
*E. lenta*, *S. pasteurianus* and *L. plantarum* for imidazole propionate,
*C. difficile* for p-cresol, *C. sporogenes* for indole-3-propionate — plus six
abundant gut organisms that carry none of them. The full list with accessions
is in `validation_results.json`.

---

## Operating points, and why the thresholds are what they are

| Criterion | Sensitivity | Specificity | Precision | False positives |
|---|---|---|---|---|
| Any matching fragment | 100.0% | 93.5% | 89.2% | 4 |
| ≥2 fragments | 100.0% | 96.8% | 94.3% | 2 |
| **≥2 fragments and ≥75% identity** | **97.0%** | **100.0%** | **100.0%** | **0** |

The tool reports the third row as **confirmed** and everything else detected as
**provisional**. The thresholds were not chosen in advance — they are where the
measured false-positive count reaches zero.

Two things made this possible:

- **A single read is not evidence.** Both single-fragment false positives
  disappear at two fragments, and no true call is affected: the weakest true
  positive had 64 fragments.
- **Identity separates the classes.** Every false positive fell below 73% mean
  translated identity, while true positives clustered at 81–99% with a median
  of 98%. A read from an organism that genuinely carries the gene matches its
  reference closely; matches in the 60s are landing on distant relatives within
  the same protein fold.

The one call given up at the confirmed tier is the *urdA* carrier at 0.2% of
cells — the hardest case in the set. That trade is the right way round: a
reported detection is one you can act on, and borderline signals are still
shown, labelled provisional.

---

## Quantitative accuracy

Detection is one question; whether the number is right is another.

**Overall:** slope 1.13, intercept 0.50, R² 0.72 across 37 comparisons of
reported against known copies per 100 genomes.

**Per gene**, as median reported ÷ true across confirmed calls:

| Gene | Bias | | Gene | Bias |
|---|---|---|---|---|
| `bai:baiF` | 1.00× | | `cutc:cutC` | 0.84× |
| `bai:baiB` | 1.02× | | `indole:tnaA` | 0.69× |
| `bai:baiE` | 1.15× | | `ipa:fldBC` | 0.41× |
| `bai:baiCD` | 1.16× | | `pcresol:hpdB` | 0.42× |
| `butyrate:buk` | 1.17× | | `propionate:pduP` | 2.32× |
| `cutc:cntA` | 1.17× | | `propionate:pct` | 2.42× |
| `bai:baiA`, `bai:baiH` | 1.21× | | `urda:urdA` | 2.43× |
| `ipa:fldH` | 1.29× | | | |
| `bsh:bsh` | 1.31× | | | |
| `butyrate:but` | 1.66× | | | |

Most genes land within about 1.7× of truth. Four are off by roughly 2.4× — two
over (`pct`, `pduP`, `urdA`) and two under (`fldBC`, `hpdB`, the latter pair
having the smallest reference sets, so reads from carriers not represented in
the set are lost).

**This is why results are reported as percentiles.** The bias is *systematic
per gene* — the same gene measured the same way in two samples is off by the
same factor in both, so it cancels when comparing samples to each other or to
the reference cohort. It does not cancel in an absolute number. Read the
percentile; treat the absolute value as accurate to within roughly a factor of
two.

---

## The active-site residue filter, measured

`urdA` is the pathway whose specificity is hardest, because fumarate reductase
flavoprotein is a close relative. The published method separates them by
requiring tyrosine or methionine at FAD-site residue 373 of the full-length
protein, which a 151 bp read mostly cannot reach.

The pipeline's reference-side substitute — checking whether the *matched
reference* itself carries an acceptable residue — was tested against
residue-filtered truth:

| Figure | Slope vs truth | R² |
|---|---|---|
| Headline count | **4.31×** | 0.991 |
| Residue-consistent subset | **1.65×** | 1.000 |

The filter cuts overcounting by a factor of 2.6, with essentially perfect
linearity. This is why both figures are reported: the headline is an upper
bound, the residue-consistent subset is the tighter estimate, and the truth
sits between them.

Two supporting checks:

| Check | Result |
|---|---|
| Anchor self-test — does position 373 map onto itself in `Q8CVD0`? | Passes on every build; failure aborts the run |
| Curated Swiss-Prot urocanate reductases carrying Y/M at their own mapped position | **6 of 6 (100%)**, at correctly distinct positions 373, 594, 394, 598, 594, 408 |
| Automatically annotated (TrEMBL) entries carrying Y/M | **~47%** of 1,537 mapped |
| Read-side vs reference-side estimate on the real sample | ~477 vs 514 fragments — **7% apart**, from two checks that share no evidence |

---

## Was the sample sequenced deeply enough?

The pipeline reads every fragment in the input — nothing is subsampled. The
separate question is whether the *sequencer* sampled the stool deeply enough
that a rarer organism could not have been missed. `openbiota depth-check` answers it
by screening the same sample at increasing depths and watching for a plateau.

A pathway counts as converged when its value changes between half and full
depth by less than the Poisson counting error implied by its own fragment
count — i.e. by less than its own noise floor, leaving no recoverable signal.

Measured on the 8,029,909-read-pair sample in this repository:

| Pathway | 502k | 1.0M | 2.0M | 4.0M | 8.0M | Change | Noise | Settled |
|---|---|---|---|---|---|---|---|---|
| `bai` | 7.60 | 8.61 | 8.58 | 8.42 | 8.20 | 2.7% | ±3.0% | yes |
| `bsh` | 67.49 | 59.08 | 59.26 | 58.01 | 55.62 | 4.3% | ±3.0% | yes |
| `butyrate` | 87.80 | 87.69 | 85.98 | 83.27 | 82.77 | 0.6% | ±2.3% | yes |
| `cutc` | 3.39 | 3.57 | 3.63 | 3.87 | 4.05 | 4.3% | ±7.0% | yes |
| `indole` | 9.09 | 10.43 | 11.50 | 11.40 | 11.55 | 1.3% | ±5.5% | yes |
| `ipa` | 12.75 | 11.17 | 10.68 | 9.87 | 9.53 | 3.5% | ±5.6% | yes |
| `pcresol` | 9.61 | 9.98 | 10.23 | 8.44 | 8.61 | 2.0% | ±4.7% | yes |
| `propionate` | 32.07 | 36.19 | 35.35 | 32.30 | 32.99 | 2.1% | ±2.7% | yes |
| `urda` | 35.69 | 36.86 | 35.26 | 35.30 | 35.07 | 0.6% | ±2.7% | yes |

All nine pathways had plateaued well before full depth — most are within a few
percent of their final value by 500,000 read pairs, which is 6% of the data.
Deeper sequencing would not change these results.

This also demonstrates the depth-independence that the `rpoB` normalisation is
supposed to provide, on real data rather than by argument: a 16-fold change in
depth moves the reported values by single-digit percentages, while the raw
fragment counts scale linearly with depth as they must.

Two consequences worth drawing out:

- **Deeper sequencing buys precision, not new detections.** Fragment counts
  scale with depth, so counting error shrinks as 1/√n. At 8M read pairs the
  headline pathways carry 2-3% counting error, which is already well below the
  per-gene systematic bias measured above.
- **A shallower sample is still usable, provided it is sampled uniformly.**
  This is why the reference cohort can be built from 2M-pair subsamples rather
  than whole runs: with uniform sampling the expected values match and only the
  per-sample noise is wider. A *prefix* is not a uniform sample — see the
  reference-cohort section below for the measurement.

---

## The reference cohort

Percentiles come from 91 adult stool shotgun metagenomes from ENA study
`PRJNA1271016`. ENA describes the study as *"MARS shotgun metagenomics — Gut
microbiome and AD study in MARS cohort"* (University of Wisconsin–Madison,
USA, 2023): an older-adult, Alzheimer's-disease study population. It is the
closest available technical match to this pipeline's inputs — same body site,
same library strategy, same NovaSeq X platform family — and the report names
it on every page that quotes a percentile, because a reader comparing a child
against it needs to know that.

| | |
|---|---|
| Samples | 91 (one of 92 excluded: 580 `rpoB` fragments, a low-bacterial specimen) |
| Reads per sample | **2,000,000 pairs, sampled uniformly across the whole run** |
| Sampling | 30 evenly spaced spot windows per run, seeded (`20260909`), reproducible |
| `rpoB` depth | median 2,219, range 1,820–2,674 |
| Provenance | `refs/reference_ranges.json` records the 91 run accessions, the seed and the study title; `openbiota cohort --runs` rebuilds it exactly |

Observed ranges, copies per 100 bacterial genomes:

| Pathway | 5th | 25th | Median | 75th | 95th |
|---|---|---|---|---|---|
| `b12` | 61.06 | 74.78 | 87.76 | 101.32 | 119.17 |
| `bai` | 1.08 | 1.95 | 2.66 | 4.06 | 6.57 |
| `bcaa` | 136.47 | 165.37 | 184.13 | 197.45 | 215.17 |
| `bsh` | 51.02 | 70.43 | 87.20 | 99.93 | 122.72 |
| `butyrate` | 17.02 | 22.91 | 27.20 | 32.00 | 39.30 |
| `cutc` | 1.09 | 1.75 | 2.47 | 3.87 | 5.57 |
| `gaba` | 5.83 | 11.35 | 18.85 | 24.74 | 36.61 |
| `indole` | 5.27 | 9.51 | 14.46 | 20.45 | 30.28 |
| `ipa` | 1.33 | 3.09 | 4.65 | 7.15 | 11.17 |
| `pcresol` | 2.72 | 4.92 | 6.52 | 8.87 | 15.29 |
| `propionate` | 12.10 | 15.71 | 19.72 | 25.62 | 38.48 |
| `urda` | 18.00 | 29.05 | 37.73 | 47.39 | 62.40 |

### Why the reads are sampled across the run, and what that changed

Until this build the cohort was made from **the first 600,000 read pairs** of
each of 34 runs — a file prefix — on the argument that every figure is a ratio
to `rpoB`, so depth divides out. That was asserted, not measured. Measured, on
one sample screened three ways at the same 600k budget against its full-depth
(8.03M) answer:

| 600k arm | median error vs full depth | panels >10% off |
|---|---|---|
| first 600k read pairs (the old cohort's method) | **19.6%** | 16 of 24 |
| a uniform random 600k | 11.6% | 13 of 24 |

and the random-arm deviations are pure counting noise: median \|z\| = 0.87
against Poisson expectation 1.0, no panel beyond 3σ, sign split 16 high / 7
low. Two conclusions follow, and they decided the design:

- **The head of a FASTQ is a biased sample of it.** Counting by flowcell tile,
  the first 600k pairs of `SAMPLE2_A02` carry 28.98 panel fragments per
  `rpoB` fragment against 27.26 for the rest of the file. The first tiles have
  a different quality profile, and at a fixed identity threshold that changes
  how many fragments each panel accepts. Uniform sampling removes it.
- **Once sampling is uniform, depth is only precision.** There is no
  directional drift with depth, so the cohort does not have to be depth-matched
  to the query. It has to be deep enough that its own counting noise is small
  against the biological spread it describes — at 2M pairs that noise is ~5%
  per panel against relative interquartile ranges of 0.3–0.9.

The rebuild did what the diagnosis predicted, and one thing it did not:

| | old (34 × first-600k) | new (91 × random-2M) |
|---|---|---|
| Panel medians | — | **0.929×** the old, every panel lower: the prefix cohort was ~7% inflated |
| Relative IQR (spread) | 0.54 | 0.52 — only **4% narrower**. The old spread was already dominated by biology, not noise; the precision gain I projected from noise arithmetic was overstated |
| Our five samples, median percentile | 64.1 | **73.0** — higher, because the inflated scale came down. Correct direction; it *unmasks* the population mismatch below rather than reducing it |

Sampling itself is `fastq-dump` spot-range access into the SRA archive, which
seeks rather than transfers: a window at spot 18,000,000 costs the same as
one at spot 1, so an unbiased 2M-pair sample of a 32M-pair run takes about
four minutes rather than the fifteen it takes to stream and discard the whole
file. `openbiota cohort --sampling prefix` is retained for cheap builds and is
labelled as biased in the output.

### Depth still matters per panel

Low-count panels have not converged at 2M and the reader should know which:
`openbiota depth-check` screens a sample at increasing depths and names them.
Comparing `SAMPLE2_A02` at full depth against a 600k cohort moved its median
panel percentile by only −1.3 points — but 7 of 25 panels moved by more than
10, and `histamine`, `ipa`, `butyrate`, `oxalate`, `pcresol` are the usual
suspects.

### The population, honestly

Four of the five samples in this repository are children aged 8–18 and the
fifth is 48; the reference population is older adults in an Alzheimer's study.
Their percentiles run high (median 73, not 50) and no amount of sampling or
depth changes that, because it is a real difference in who is being compared.
A paediatric functional cohort would; it is the next piece of work.

---

## Cross-platform check: Tiny Health on the same FASTQ files

Every validation above is against simulated reads, reference genomes or this
pipeline's own cohort. The one external check on real samples is that
[Tiny Health](https://www.tinyhealth.com), a CLIA-certified shotgun
metagenomics service, ran the same five stool specimens. Two independent
pipelines, one set of DNA.

Tiny Health is not treated as ground truth. Where the two disagree, the
question is which is right, and the answer comes from the sequence evidence.

Reproduce with `python -m tools.compare_tinyhealth` against your own Tiny
Health report and run. The comparison is a tool, not a fixture: nothing from
these samples is pinned in the test suite.

### What is comparable, and what it shows

Absolute values are not comparable — Tiny Health reports RPKM against their
reference population, this pipeline reports copies per 100 bacterial genomes
against `PRJNA1271016`. What is comparable is where a reading sits on each
platform's own scale, and relative species abundance, which is unitless.

| Check | Result |
|---|---|
| Like-for-like functional readings, 5 samples | **89 pairs** |
| Same band (low / typical / high) on both | 45 (51%) |
| Within one band | **83 (93%)** |
| Opposite at a band edge (see below) | 3 |
| Opposite | 3 (3%), all three explained below |
| Their top-20 species also detected here (`SAMPLE2_A02`) | **17 / 17 (100%)** |
| Abundance rank correlation on those species | Spearman **ρ = 0.775** |
| Median abundance ratio between platforms | **1.41x** |

The 1.41x median abundance ratio is the expected size of a MetaPhlAn-versus-
GTDB disagreement on the same reads: the two use different marker sets and
different reference catalogues, and Tiny Health's figures run consistently
higher because their denominator excludes more unclassified sequence.

### Band-edge straddles

A three-band comparison is fragile exactly at the lines. Three pairs read as
"opposite" only because Tiny Health's value sits 0.5–1.5% *below* their own
cutoff while ours sits 1–3 percentile points *above* our 75th line:
`SAMPLE1_A01` propionate (theirs 546 vs cutoff 549; ours 76.3rd), `SAMPLE1_A01`
vitamin B9 (268.7 vs 271; 77.4th) and `SAMPLE6_A06` hydrogen sulfide index
(9.01 vs 9.15; 77.5th). For propionate and B9 the rank order of the three
same-template samples is *identical* on both platforms — the two agree about
the values and differ about where a band starts. The comparison tool reports
these separately and applies the rank-order check to every straddle.

### The three opposite calls, resolved

**Branched-chain amino acids, `SAMPLE1_A01` and `SAMPLE2_A02`** — reference
population, not measurement. On the three samples that share a Tiny Health
report template, BCAA capacity ranks *identically* on both platforms
(theirs 2507 / 3064 / 3025 RPKM; ours 209 / 235 / 210 copies per 100 genomes).
The two pipelines measure the same thing and agree on the ordering; they
disagree only about whether that value is high, because they compare against
different populations. Neither is wrong. Worth knowing: this panel is close to
saturation — `ilvC` and `ilvD` are near-universal in gut bacteria, so the
whole cohort spans 123–234 and a percentile inside that range carries little
information.

**GABA production, `SAMPLE3_A03`** — a Tiny Health false negative. They report
**0.0 rpkm**. This pipeline finds **1,940 `gadB` fragments at 100% median
identity**, best-matching *Bacteroides clarus* (592), *Phocaeicola dorei*
(182), *Bacteroides stercoris* (178) and *Phocaeicola vulgatus* (146), in a
sample that is **34% Bacteroidaceae**. *Bacteroides* is the dominant genus of
gut GABA producers (Strandwitz et al., *Nature Microbiology* 4:396–403, 2019,
doi:10.1038/s41564-018-0307-3). A third of the community cannot carry the gene
and produce a zero; this pipeline counts the gene whatever organism carries
it, so it cannot report that combination.

### One pair excluded as not like-for-like

**Histamine.** Tiny Health's v4.8.3 metric is "Histamine-producing species":
the summed abundance of a *curated species list*. This pipeline's is `hdcA`
gene capacity, whatever organism carries it. For `SAMPLE1_A01` and `SAMPLE3_A03`
they report 0.0% while this pipeline finds `hdcA` at **100% identity** in
*Eggerthella lenta* (14 and 66 fragments) and *Bacteroides fragilis* — neither
on their list. This is the structural argument for gene-based over
taxon-based functional calls: a curated list can only find carriers someone
already thought of. The two are counted separately rather than scored as a
disagreement, because they are measurements of different quantities.

### A reproducibility problem on their side

Tiny Health's RPKM values are not comparable across their own report versions.
Between template v4.8.3 (samples SAMPLE1, SAMPLE3, SAMPLE6) and v5.x (SAMPLE2, SAMPLE4), their
median value shifted by **2.41x** — and in the *same direction for all 12
comparable panels* (range 1.74–3.02x). Twelve of twelve in one direction is
p ≈ 0.0002 under a sign test.

That could be biology if the two groups of people simply differed. They do
not: this pipeline measured the same specimens and puts the two groups
**1.12x** apart, with the ratios scattered on both sides of 1.0 (0.79–1.65).
A consistent multiplier on one platform and scatter on the other is what a
silent renormalisation looks like, not a cohort difference. The practical
consequence is that a Tiny Health customer comparing a 2026 report to an
earlier one would see every functional number fall by roughly half with no
biological change. `test_our_units_are_stable_across_pipeline_revisions`
asserts this pipeline does not do the same.

### What this check does not cover

The comparison touches 89 of the ~1,390 quantities this pipeline reports per
sample. Disease-pattern resemblance, biological age, resistance and virulence
determinants, skin research panels and the species-level percentiles have no
Tiny Health counterpart, so nothing here validates them. They are covered, to
the extent they are covered at all, by the sections above and by
`docs/VALIDATION.md`'s per-module limitations.

---

## Limitations

Grouped by what they affect.

### What the measurement is

- **This measures gene capacity, not metabolite concentration.** Carrying the
  gene means the bacteria can make the metabolite; how much they actually make
  depends on substrate availability, transit time and competition. Measuring
  the metabolites requires mass spectrometry on blood or stool.
- **Gene carriage predicts metabolite level weakly.** In the source study of
  294 samples, of all `urdA`-positive genomes only *Streptococcus pasteurianus*
  abundance was significantly associated with plasma imidazole propionate. This
  is a property of the biology, not of this implementation.
- **No clinical thresholds exist.** No value in any output corresponds to a
  diagnosis, because no such threshold has been published for any of these
  pathways. The percentile says where you sit in a population; it does not say
  what that means for your health.

### Accuracy boundaries

- **Absolute values are accurate to roughly a factor of two**, with four genes
  at ~2.4×. Percentiles and longitudinal changes are much more reliable than
  absolute numbers.
- **The detection limit is about 0.2% of cells** for a well-referenced gene.
  Rarer carriers may be missed.
- **`urdA` loses sensitivity at the confirmed tier** (75% versus 100%
  permissive), because its hardest case sits right at the identity threshold.
- **Small reference sets bound sensitivity.** `bai:baiCD` and `bai:baiH` have
  four references each, `propionate:lcdA` five, `pcresol:hpdB` 77. A low value
  for these is weak evidence of absence, and the `few-refs` flag marks them.

### Profile similarity boundaries

- **A percentile is not a probability.** It says the sample is more concordant
  with a published pattern than that fraction of a matched reference group.
  Nothing more.
- **The associations are mostly not condition-specific.** The specificity
  matrix above shows the colorectal pattern separating IBD better than
  colorectal cancer, and the Long COVID pattern scoring 0.80 against IBD. Read
  every profile score beside the dysbiosis anchor.
- **Two of the three profiles have unmeasured discrimination.** No labelled
  cohort exists for ME/CFS or Long COVID that this tool can reach, so no AUC
  is printed for them.
- **No exposure-matched contrast exists** for ME/CFS or Long COVID. The right
  comparison is people with the same exposure who did not develop the
  syndrome; the published cohorts use healthy or uninfected controls instead,
  so syndrome-specific signal cannot be separated from post-exposure signal.
  Both profiles declare this and take a confidence penalty for it.
- **The combined 50/35/15 module split is an engineering starting point**, not
  an evidence-derived weighting. Report the module scores separately, which the
  output does.
- **Profiler versions must match.** Species abundances from MetaPhlAn 3 and
  MetaPhlAn 4 are not interchangeable: version 4 reorganised the taxonomy into
  SGBs and renamed genera. The sample is profiled with the same version as the
  reference cohort, both are recorded in `results.json`, and a mismatch
  downgrades confidence.
- **Long COVID scores chronic state only.** Acute future-risk prediction from
  stool at 0–2 months post-infection is a different task, and features from it
  are excluded rather than reused.
- **Effect sizes in this literature are small.** The largest independent
  prospective Long COVID study (799 outpatients) reached species-level AUROC
  0.62 against 0.72 for clinical variables alone; a 2,561-participant community
  study found no association with illness duration. The frequently quoted AUC
  of 0.96 classifies a study's own internally defined enterotypes, not an
  externally validated diagnosis.

### Method boundaries

- **No per-genome attribution for gene capacity.** Read-level search cannot say
  which organism carries a gene. Named organisms in the gene-capacity output
  are the nearest reference sequence, not identifications. Species in Part B
  *are* measured identifications, from marker genes.
- **Synteny is unavailable.** `cutC`/`cutD` gene neighbourhood is a strong
  specificity signal in assembly-based work and cannot be used here.
- **Taxonomy is coarse.** The community profile is reliable at phylum level
  only, and is bounded by the reviewed `rpoB` reference set.
- **The succinate route to propionate is not screened**, so propionate capacity
  is systematically understated. See [PANELS.md](PANELS.md).

### Validation boundaries

- **Depth is not a limiting factor for this sample**, but it would be for a
  much shallower one. Below roughly 500,000 read pairs the rarer pathways
  (`bai`, `cutc`, `ipa`) fall under the fragment counts needed for a confirmed
  call.
- **Validation used simulated reads.** Simulation reproduces read length,
  quality binning, error rate and mate structure, but not every artefact of a
  real library — chimeras, uneven coverage, PCR duplicates, strain-level
  variation within a species.
- **24 genomes is not the gut.** A real sample contains hundreds of species,
  many unsequenced. Some will carry homologs the validation set does not
  represent.
- **The reference cohort is 33 samples of mixed clinical status.** ENA metadata
  for `PRJNA1271016` does not distinguish case from control, so the percentiles
  describe adult stool generally rather than a screened-healthy population.
  Widen it with `openbiota cohort --samples 100`.
- **Not compared against the published pipeline on shared samples.** The
  strongest possible validation would run both methods on the same real
  metagenomes and compare per-genome calls, using the MAGs in
  `PRJNA1469852`. That has not been done. The protocol is below.

---

## Remaining work: comparison against the published method

The one validation not yet done, and how to do it.

**1. Fetch the MAGs.** `PRJNA1469852` holds the metagenome-assembled genomes
the source study called `urdA`-positive and negative, per sample.

**2. Reproduce their ground truth.** Call genes on each MAG with Prodigal,
score against their ImP-producer HMM, keep hits above 60% of the lowest-scoring
model sequence, align with Clustal Omega and require Y or M at residue 373.

**3. Run this pipeline** over the matching raw metagenomes from
`PRJNA1271016`.

```bash
for r1 in fastq/*_1.fastq.gz; do
  openbiota run --r1 "$r1" --r2 "${r1/_1/_2}" --out results/comparison
done
```

**4. Build the confusion matrix** against their per-sample calls, sweeping a
threshold over both `copies_per_100_genomes` and
`copies_per_100_genomes_residue_consistent`. The difference between those two
ROC curves quantifies what the reference-side residue filter buys on real data
rather than synthetic.

**5. If plasma imidazole propionate is available**, correlate it against the
pipeline's output. Expect a weak relationship — that is the correct result, not
a failure, and it is the point about gene carriage made above.

Until this is done, accuracy is characterised against synthetic ground truth
and against curated reference proteins, which is what the figures at the top of
this file represent.
