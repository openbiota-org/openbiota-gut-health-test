# Method

How the measurement works, and why each design choice was made.

---

## 1. Read-level translated search

Reads are searched directly with `diamond blastx` against a curated protein
database. Nothing is assembled.

The reason is coverage. Assembling a single sample only recovers organisms
abundant enough to reach roughly 8–10× coverage. At 2.4 Gbp, a 3 Mb genome must
be about 1% of the community to clear that bar. Several target organisms sit
well below it — 7α-dehydroxylating Clostridia are typically under 0.1% — so an
assembly-first pipeline would not see them at all. Read-level search finds the
gene wherever it occurs.

The trade-off is that a 151 bp read translates to ~50 residues, so full-length
protein context is lost. Section 3 is how that is handled.

## 2. Normalisation to genome equivalents

Raw hit counts scale with sequencing depth and with gene length, so they are
not comparable between samples or between genes. Everything is divided by
`rpoB`, the universal single-copy bacterial RNA polymerase beta subunit:

```
copies_per_genome = (target_fragments / mean_target_ref_length_aa)
                  / (rpoB_fragments   / mean_rpoB_ref_length_aa)
```

reported ×100 as **copies per 100 bacterial genomes**. For a single-copy gene,
100 means essentially every genome carries it, 10 means roughly one in ten.

Three properties follow directly:

- **Depth independence.** Two samples at different depths are comparable. This
  is measured, not assumed — see the depth-independence check in
  [VALIDATION.md](VALIDATION.md).
- **Gene-length correction.** A 1,200 aa target is comparable to a 400 aa one.
  Longer genes intercept proportionally more reads; dividing by mean reference
  length removes that.
- **Host DNA cancels.** Human reads carry no bacterial `rpoB`, so they inflate
  neither numerator nor denominator. This is why there is no host-filtering
  step — it would be redundant.

If `rpoB` comes out at zero the result is reported as indeterminate. No
per-million-reads substitute is offered in its place.

## 3. Specificity

Every target here has close homologs. Four mechanisms handle them.

### One combined database

All panels' targets, all decoys and `rpoB` live in a single DIAMOND database,
searched in one pass with `--max-target-seqs 1`. A read is credited to a target
only when that target beats every decoy *and* every other panel's target.

This makes cross-panel competition automatic — most usefully for `cutC` versus
`hpdB`, which are both glycyl-radical enzymes and both targets of different
panels. It also means the FASTQ is translated once rather than once per panel,
which is where essentially all the runtime goes.

### Role precedence

A protein claimed by `rpoB` or by a target can never also appear in a decoy
set. Without that guard, competition between identical sequences would be a
coin flip.

### Reference length windows

Each reference set is trimmed around its median length, removing database
fragments and multi-domain fusions. Both would otherwise corrupt the length
term in the denominator.

### Active-site residue checks

Where a panel defines one — currently `urdA` — two independent checks run.

**Reference side.** At database build time, every target reference is aligned
to a designated anchor protein with `diamond blastp`, and the position
corresponding to the canonical active-site position is recorded, along with
whether that reference's own residue there is acceptable. Fragments matching a
reference that fails are reported as a separate, weaker subset. This applies to
*every* fragment.

**Read side.** For the minority of fragments whose alignment spans the mapped
position, the aligned residue is read off and tested directly. This is the
published filter, applied where it is available.

The reference-side check is calibrated against curated data: all 6 Swiss-Prot
urocanate reductases carry Y or M at their own mapped position, at correctly
*distinct* positions (373, 594, 394, 598, 594, 408), while only ~47% of
automatically annotated entries do. An anchor self-test runs on every build —
the anchor must map position 373 onto itself, or the build aborts.

## 4. Mate collapse

The two mates of one DNA fragment are not independent observations; counting
both would inflate every result by up to 2×.

Hits are keyed by query id across both mate files and the higher-scoring hit
wins. On NovaSeq output the read ids are byte-identical between R1 and R2 with
no `/1` or `/2` suffix, so a set union suffices — but suffix-stripping is
implemented for input that does carry mate markers.

Organism attribution is derived from mate-collapsed fragments, so attribution
counts reconcile with the fragment total rather than the read total. This is
asserted at runtime.

## 5. Reference cohort

A single number has no meaning without a population to compare it against.
`openbiota cohort` downloads read prefixes from a public ENA study, screens each
sample through the identical pipeline, and reduces the results to percentiles.

The default cohort is ENA `PRJNA1271016` — the MARS cohort, 294 adult stool
shotgun metagenomes from an older-adult Alzheimer's-disease study — same body
site, same library strategy, and the same NovaSeq X platform family. The
active build uses 91 of them. Samples below a minimum `rpoB` depth are excluded
so that imprecise samples do not widen the percentiles artificially.

Each cohort sample contributes **2,000,000 read pairs sampled uniformly across
its whole run** — 30 evenly spaced spot windows, seeded — rather than the
first 2M. Every reported figure is a ratio to `rpoB`, so depth divides out to
first order; but the head of a FASTQ is not a random sample of it (flowcell
position changes which fragments pass a fixed identity threshold), and a
prefix cohort was measured to be ~7% inflated on every panel. The measurement
and the before/after are in `docs/VALIDATION.md`. `reference_ranges.json`
records the run accessions, seed and study title, so `openbiota cohort --runs`
rebuilds the cohort exactly.

## 6. Per-panel decisions

Panel definitions are declarative YAML, and each records why it is built the
way it is. The decisions that most affect the numbers:

| Panel | Decision | Reason |
|---|---|---|
| `cutc` | Query on protein name, never `gene:cutC` | `gene:cutC` returns ~6,300 entries that are mostly the *copper homeostasis* protein CutC. Same trap for `gene:cntA`, which returns a staphylopine transporter. |
| `bai` | Headline is the median of `baiCD` and `baiH` only | The other four genes' UniProt names are heavily over-applied by automated annotation. `baiE` hits landed on *Phocaeicola*; `baiF` on an alkaliphilic soil *Bacillus*. Neither performs the chemistry. |
| `bai` | `baiCD` = 3-oxocholoyl-CoA 4-desaturase | The obvious query, on `"3-oxo-Delta(4,5)-steroid 5-beta-reductase"`, is doubly wrong: the comma inside the parentheses silently breaks UniProt's boolean parser, and the hits are mostly *plant* progesterone 5β-reductases. |
| `butyrate` | Positive control | Butyrate producers are abundant in any stool sample. A near-zero result means the run is broken, and is raised as an error. |
| `butyrate` | `ackA` is a mandatory decoy | Acetate kinase is homologous to butyrate kinase and far more widespread. |
| `indole` | `tpl` is a mandatory decoy | Tyrosine phenol-lyase is `tnaA`'s closest relative — same fold, same length, different product. |
| `ipa` | Headline is `fldBC` only | `fldH` shares a fold with the near-universal D-lactate dehydrogenases. |
| `propionate` | Succinate route excluded | Its marker genes sit in the core methylmalonyl-CoA pathway, so nearly every genome would score positive. The panel covers the acrylate and propanediol routes only and therefore understates total capacity. |
| `pcresol` | `hpdB` large subunit only | The unqualified name query returns a mix of large subunit, an 85 aa small subunit and an activase, whose combined median length is meaningless. |

## 7. The taxonomic engine

The gene-capacity engine answers "how much of gene X is present". Disease
profiles need "how much of species Y is present", which requires marker-gene
profiling. The `rpoB`-derived composition in section 5 of the report is a
useful by-product of normalisation, but it is marker-limited: reliable at
phylum, indicative at genus, and not usable at species.

**MetaPhlAn 3, pinned.** The species reference cohort comes from
curatedMetagenomicData, whose profiles were computed with MetaPhlAn 3 against
the CHOCOPhlAn 201901 markers. Species abundances from different profiler
versions are not comparable — version 4 reorganised the taxonomy into
species-level genome bins, renamed genera (*Bacteroides vulgatus* became
*Phocaeicola vulgatus*) and changed detection behaviour. A sample profiled
with one version and scored against a reference built with another measures
the version difference. So the sample is profiled with the same version as the
reference, both are recorded in `results.json`, and a mismatch downgrades
confidence. The GMWI2 dysbiosis anchor was also trained on MetaPhlAn 3 output,
which settles the choice.

**Host reads are removed here, and only here.** The gene-capacity pipeline
needs no host removal: human reads carry neither the pathway genes nor `rpoB`,
so they cancel in the ratio. Marker-gene profiling has no such property —
human reads inflate the denominator and can occasionally hit markers. Reads
are filtered against GRCh38 plus PhiX with Bowtie2 before profiling, the host
fraction is reported as a QC field, and only the non-host reads are written to
disk. Consumer-kit FASTQs often arrive with human reads already stripped by the
provider — the five kits this was developed on contain not one human pair in
eight million, which real stool never does — so an exact zero in a deep sample
is reported as "removed upstream", not as a measurement of the sample.

**Scale must match.** MetaPhlAn with unknown estimation reports species as a
share of the whole sample, so they sum to 100 minus the unclassified fraction
— about 50% for a typical stool sample. The reference profiles carry no
unknown row and sum to 100. Because the compositional transform's zero floor
is an absolute constant, that mismatch shifts every value, and an *absent*
taxon ends up scoring near the 80th percentile of a reference where it is also
absent. Sample abundances are therefore rescaled to the classified fraction
first, which is the same normalisation the published GMWI2 implementation
applies.

**Two lanes, one scored.** The MetaPhlAn 3 lane above is the *scoring* lane:
sample and reference are named the same way, so percentiles mean something.
It is also the older catalogue, and it leaves 30–55% of a typical stool sample
unclassified. A second, *extended* lane profiles the same host-filtered reads
with MetaPhlAn 4 against CHOCOPhlAnSGB (Jun23) — ~26,000 species-level genome
bins built from a million genomes, including ~4,900 organisms with no cultured
representative. On the first sample run through both, the extended lane named
191 genome bins across 177 species where the scoring lane named 105, and left
14% unclassified where the scoring lane left 48%. It also resolves a species
into its constituent bins — *F. prausnitzii* into six — which is the closest
a read-level profiler comes to strain resolution.

The extended lane is reported in full (section 9.2 of the report lists every
bin) and **never scored**: there is no reference cohort profiled on that
catalogue, so a bin can be counted but not placed. The two inventories are
bridged by name, tolerating the renames the SGB catalogue applies to some
genera (*Bacteroides* → *Phocaeicola*) and not others, and the report says
which species each lane found alone. When a MetaPhlAn 4-profiled reference
cohort of adequate size exists, the scoring lane can move to it; until then
the older catalogue is the honest one to score against.

**Bowtie2's gzip writer is checked.** Host filtering writes only the unaligned
pairs, compressed, via `--un-conc-gz`. That writer has intermittently produced
files with trailing non-gzip bytes under multi-threaded runs, which MetaPhlAn's
reader rejects. Every output is decompressed end-to-end before use — on a cache
hit as well — and a corrupt file triggers one re-run with plain output
compressed in-process.

## 7a. Read preprocessing and QC gates

Before any alignment, `fastp` makes one pass over the reads. By default it
assesses and does not alter them: duplication rate, Q20/Q30, GC, read-length
distribution, adapter content by read-pair overlap, insert-size peak, polyG
tails. `--trim` writes trimmed reads and uses them downstream; adapter
auto-detection (`--detect-adapters`) is opt-in because it is single-threaded
and the overlap method already handles paired-end adapters. The result is
cached by input fingerprint, so a re-run costs nothing.

The FASTQ headers are read for provenance — instrument, run, flowcell, lane and
index for Illumina headers; accession for SRA — and the flowcell:lane set forms
a *batch key*. Two samples sharing it were sequenced together; two that differ
carry a batch difference that can exceed the biology, and the report says so.

Every stage's measurements are assembled into the spec 4.2 gate table: input
pairs, pairs passing filters, host fraction, usable non-host pairs (the
prespecified minimum is 500,000 and a profile abstains below it), classified
fraction, duplicate rate, Q30, adapter content, read length, GC, contamination
flags, `rpoB` fragment count, and sequencing batch. A gate that could not be
measured is `unknown`, never `pass`.

`fastp` itself is invoked as a binary. Its source is C++ with SIMD kernels
and cannot be imported into Python; what it offers this pipeline is the QC
pass and, optionally, the trimming, not a faster FASTQ reader for the stages
that follow — those are DIAMOND, Bowtie2 and MetaPhlAn, which read the files
themselves.

## 7b. Group features: genus sums and carrier sets

Some published evidence names a genus (16S studies) or a functional guild
("butyrate producers", "oral-origin taxa") rather than a species. Both are
scored as the log-ratio of a *sum* of species abundances against the same
prevalent-core frame the species CLR uses, computed identically for the sample
and for every reference sample, so the null is built the same way as the
value. The curated guilds live in `taxa/*.yaml` and are also what the report's
"microbial groups" section reads from; a panel's `organisms` list serves the
same role for a metabolite pathway, giving an independent carriage-based line
of evidence beside the gene count.

## 8. Compositional comparison

Relative abundances are compositional: they sum to a constant, so an apparent
increase in one taxon can be produced entirely by decreases elsewhere. Sums
and differences of raw abundances are not valid, and the centred log-ratio is
what makes them valid.

Three details matter more than they look, and each was found by validation
rather than reasoning:

**The zero floor is a fixed constant.** Zeros need replacing before taking
logs. The textbook rule is half the smallest observed value, which makes the
floor depend on the sample: a deeper sample detects rarer taxa, gets a lower
floor, and maps its zeros further down. For a taxon absent in most samples
that difference *is* the signal, so the floor is a constant (1e-4 percent)
applied identically to sample and reference.

**The denominator uses a stable core.** The textbook CLR divides by the
geometric mean of every component, which makes a near-floor taxon's value
depend on how many other rare taxa the sample happened to detect. Dividing
instead by the geometric mean of taxa prevalent in the reference removes that
noise. On the colorectal validation set it recovered a third of the gap
between whole-composition CLR and the theoretical ceiling for those features.

**Absence is a category, not a number.** Even with a fixed floor and a stable
frame, two samples that both lack a taxon get slightly different values,
because the denominator still varies. Comparing those values compares
sequencing incidentals. For a zero-inflated feature — and the best-replicated
colorectal markers are present in 1–4% of healthy people — a sample lacking
the taxon is instead placed at the mid-rank of the reference samples that also
lack it, and a sample carrying it ranks above all of them. Before this change
five absent markers scored at the 74th percentile; after it they score at the
48th, and the profile score moved from a spurious 74th to 23rd.

## 9. Comparison with the published assembly method

Vemuganti et al. assembled and binned. This pipeline does not.

| Step | Published method | This pipeline | Consequence |
|---|---|---|---|
| Assembly | SPAdes `-k 21,33,55,77` | none | Finds low-abundance organisms the assembly misses; loses per-genome resolution |
| Binning | MetaBAT2 + CheckM + dRep | none | No per-genome calls |
| Taxonomy | GTDB-Tk on bins | coarse `rpoB` profile, phylum level | Much weaker |
| Gene calling | Prodigal on contigs | reads searched directly | No full-length context |
| Homology model | HMM via MUSCLE + hmmbuild | curated reference set + blastx | Comparable curation, weaker scoring model |
| Threshold | 60% of lowest-scoring model sequence | identity floor + minimum alignment length | A length-blind proxy |
| Y/M at FAD site 373 | on full-length proteins | reference-side on all fragments, read-side where the read reaches the site | See VALIDATION.md |
| Abundance | kallisto TPM + CLR | fragments / gene length, over `rpoB` | Different scale, same intent |

### Why the read-level approach is defensible

- **Detection at low abundance** is the real argument. Genes in organisms below
  assembly coverage are invisible to the published pipeline and visible here.
- **Normalisation needs no contigs.** `rpoB` per genome equivalent works
  directly on reads and cancels host contamination for free.
- **Longitudinal use.** Systematic biases are held constant across runs of the
  same pipeline, so a change over time is far more interpretable than any
  single absolute value.
- **Cost.** Minutes on one workstation, versus hours of assembly and binning.

### Where it is weaker

- No per-genome or per-organism attribution.
- Weaker wherever the discriminating feature is a specific residue or a gene
  neighbourhood — `urdA` (residue 373) and `cutC` (`cutC`/`cutD` synteny).
- Taxonomy is coarse and bounded by the reference set.

## 10. Deviations from the source specification

Two, both documented where they occur in the code.

**Rath et al.'s cutC/cntA protein FASTA is not published.** All eleven
supplementary files were checked: they contain genome and taxonomy lists (1,107
`cutC`-positive and 6,738 `cntA`-positive genomes), phylogenetic trees and
figures. There is no sequence database to fetch. UniProt's curated
recommended-name assignment is used instead, which is itself a curation layer
built on family membership and signature residues, and yields a comparable set
size (~370 versus their 454).

**Their HMM cutoffs are not transferable.** The published thresholds (bit score
906.4 for `cutC`, 440.6 for `cntA`) are defined on full-length predicted
proteins. A 33–50 residue read fragment cannot reach that score regardless of
match quality. Applying them would reject everything; rescaling them by length
would be inventing a threshold and attributing it to them. Neither is done.

**One correction to the specification.** The `urdA` residue-check anchor was
given as `Q8EAP8`. That accession is a 180 aa YdbS-like PH-domain protein and
cannot carry a residue 373. The urocanate reductase of *S. oneidensis* MR-1 is
`Q8CVD0` (582 aa, gene `urdA`), whose residue 373 is tyrosine — the expected
Y/M identity, which independently confirms both the anchor and the position.

---

## References

- Vital M., Howe A.C., Tiedje J.M., *mBio* 5:e00889-14 (2014) — butyrate pathways
- Rath S. et al., *Microbiome* 5:54 (2017) — TMA-producing bacteria
- Ridlon J.M. et al., *J Lipid Res* 47:241–259 (2006) — bile acid transformations
- Reichardt N. et al., *ISME J* 8:1323–1335 (2014) — propionate pathways
- Dodd D. et al., *Nature* 551:648–652 (2017) — aromatic amino acid metabolism
- Koh A. et al., *Cell* 175:947–961 (2018) — imidazole propionate
- Case R.J. et al., *Appl Environ Microbiol* 73:278–288 (2007) — `rpoB` as a marker
- Buchfink B. et al., *Nat Methods* 18:366–368 (2021) — DIAMOND

Per-panel citations are in each `panels/*.yaml` and are reproduced in the
report.
