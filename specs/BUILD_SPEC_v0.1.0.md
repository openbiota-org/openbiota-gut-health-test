# BUILD SPEC — openbiota-gut-health-test

You are building a complete, production-quality bioinformatics tool in this
repository. Read this entire document before writing any code. Ask me before
deviating from the scientific method described in Part 3; those choices are
deliberate and several of them are counterintuitive.

---

## PART 1 — GOAL AND SCOPE

### What this tool does

Quantifies the abundance of specific bacterial metabolic pathway genes in a shotgun
stool metagenome, normalized to bacterial genome equivalents, so the result is
comparable across samples and sequencing depths.

The initial target is `urdA` (urocanate reductase), the enzyme producing imidazole
propionate (ImP). Additional gene panels follow. The architecture must make adding a
gene panel a config change, not a code change.

### What this tool must never claim to do

- **It does not measure metabolites.** ImP, TMAO, butyrate and the rest are small
  molecules quantified from plasma or stool by mass spectrometry. They are not
  recoverable from DNA sequence by any method. The tool measures genetic *capacity*.
- **It does not estimate disease risk.** No output, log line, docstring, or README
  sentence may imply otherwise.
- **It has no validated reference range.** If the report prints qualitative bands
  ("low", "moderate"), each must be explicitly labeled in the output text as an
  arbitrary heuristic with no clinical meaning. Do not invent thresholds and present
  them as if derived from literature.

Treat these three as hard constraints, not stylistic preferences.

---

## PART 2 — INPUT DATA (already on disk)

```
openbiota-gut-health-test/
└── fastq/
    ├── A02_1.fastq      (~2.6 GB, uncompressed)
    ├── A02_1.fastq.gz
    ├── A02_2.fastq      (~2.6 GB, uncompressed)
    └── A02_2.fastq.gz
```

Verified characteristics — do not re-derive these, but do validate them at runtime:

| Property | Value |
|---|---|
| Read pairs | 8,029,909 |
| Total bases | ~2.38 Gbp (1.19 Gbp per mate file) |
| Read length | 100–151 bp, mean 148.4 (variable = already trimmed) |
| Platform | Illumina NovaSeq X (`LH00469` instrument prefix) |
| Quality encoding | Phred+33, binned (3–4 distinct characters — normal, not corruption) |
| Data type | Shotgun metagenomic, paired-end |

**Critical detail that affects correctness:** read IDs are *identical* between the R1
and R2 files. There is no `/1` or `/2` suffix and no space-delimited `1:N:0:INDEX`
field. Example ID present in both files: `LH00469:636:23VMTGLT4:8:1101:1003:20056`.

Consequences you must handle:

1. Mates from one DNA fragment are **not independent observations**. Counting both
   inflates results by up to 2x. Deduplicate hits by query ID across both files
   before counting. Since IDs are already identical, a set union on `qseqid` is
   sufficient — but write the suffix-stripping logic anyway so the tool works on data
   that does use `/1` and `/2`.
2. Files are order-synchronized (last record in both is
   `LH00469:636:23VMTGLT4:8:2498:52140:8019`). Do not assume this for arbitrary
   input; do not rely on it either.

**Preprocessing:** none required. The variable read lengths confirm the provider
already quality- and adapter-trimmed. Do not add fastp, Trimmomatic, or cutadapt.
Do not filter human reads — see Part 3 for why that is unnecessary.

---

## PART 3 — SCIENTIFIC METHOD

### 3.1 Core design: read-level translated search

Do **not** assemble. Search the reads directly with `diamond blastx`.

Rationale, which you should preserve in a comment: assembling a single sample only
recovers organisms abundant enough to reach roughly 8–10x coverage. At 2.38 Gbp, a
3 Mb genome needs to be about 1% of the community to clear that bar. Several target
organisms sit well below it. Read-level search detects the gene wherever it occurs.
The cost is losing full-length protein context, which matters for specificity (3.3).

### 3.2 Normalization — this is the most important design decision

Raw hit counts are meaningless. They scale with sequencing depth and with the length
of each gene. Normalize against `rpoB` (RNA polymerase beta subunit), a universal
single-copy bacterial gene:

```
copies_per_genome = (target_fragments / mean_target_ref_length_aa)
                  / (rpoB_fragments   / mean_rpoB_ref_length_aa)
```

Report as **copies per 100 bacterial genomes**.

This gives three properties for free:
- Depth-independent, so samples are comparable.
- Gene-length-corrected, so a 1200 aa target is comparable to a 400 aa one.
- **Human read contamination cancels automatically.** Human reads carry no bacterial
  `rpoB`, so they inflate neither numerator nor denominator. This is why no host
  filtering step is needed.

If `rpoB` fragments are zero, the result is *indeterminate*. Say so explicitly. Do
not fall back to a per-million-reads figure and present it as equivalent.

### 3.3 Specificity — the hard part, handled per gene family

Every target here has close homologs that will produce false positives. Handling
differs by family. **Do not apply a single generic approach.**

**`urdA` (urocanate reductase).** Close homolog: fumarate reductase flavoprotein
(FrdA), same flavoprotein family. The source paper discriminated them by requiring
tyrosine or methionine at FAD-binding-site residue 373 on full-length predicted
proteins. On 151 bp reads (~50 aa) that filter is mostly unavailable.

Implement two layers:
- *Required:* include FrdA sequences in the search database as decoys and accept a
  read only when its best hit is a UrdA reference. Best-hit competition against
  decoys suppresses obvious false positives.
- *Stretch, implement if feasible:* a partial residue check. At database build time,
  align all UrdA references to the *Shewanella oneidensis* UrdA reference and record,
  for each reference sequence, which of its own residue positions corresponds to
  canonical position 373. At search time, request subject alignment coordinates in
  the DIAMOND output; for the subset of reads whose alignment spans that mapped
  position, verify the aligned residue is Y or M. Report the checked subset
  separately — e.g. "412 accepted, of which 31 spanned residue 373; 29 passed."
  This partially recovers the paper's filter. Verify which output fields your DIAMOND
  version supports (`diamond help`) before relying on translated query sequence.

Until the stretch filter is in place, **the urdA count is an upper bound** and the
report must say so.

**`cutC` (choline TMA-lyase) — hardest target, do not improvise.** CutC is a glycyl
radical enzyme. That family also contains pyruvate formate-lyase (`pflB`), glycerol/
propanediol dehydratase, and 4-hydroxyphenylacetate decarboxylase (`hpdB`) — which
is itself a target in another panel. A naive UniProt keyword fetch plus identity
threshold will produce a confidently wrong number.

Use the curated reference set from Rath et al. 2017 (*Microbiome* 5:54,
doi:10.1186/s40168-017-0271-9), which screened 67,134 genomes on three criteria: HMM
similarity, conservation of signature amino acid residues, and phylogenetic distance
to top-scoring sequences, with synteny to partner genes as a further check. They
published explicit discrimination cutoffs: HMM score 906.4 for `cutC`, 440.6 for
`cntA`. Their final sets contain 454 unique cutC proteins and 491 unique cntA
proteins.

Fetch their supplementary databases rather than constructing your own. If HMMER is
available, prefer `hmmsearch` with their cutoffs over DIAMOND identity thresholds for
this panel specifically. Note in code and README that `cutC`/`cutD` synteny — a
strong specificity signal in the original work — requires contigs and is unavailable
in read-level search.

**`bai` operon (secondary bile acids).** No severe homolog problem, but 7α-dehydroxylating
Clostridia typically sit under 0.1% of the community. Expect single-digit fragment
counts where Poisson noise dominates. The report must flag any panel whose fragment
count falls below a configurable minimum (default 20) as statistically unstable
rather than printing a precise-looking normalized figure.

**`but` / `buk` (butyrate).** High abundance, well-behaved. Use this as the positive
control — if this panel returns near zero, the pipeline is broken.

### 3.4 Gene panels

Implement `urdA` first and verify end-to-end. Then add, in this order: `cutC`/`cntA`,
`but`/`buk`, `hpdBCA`, `bai`/`bsh`, `tnaA`/`fldH`.

**Deliberately excluded — do not add these, and document why in the README:**
`ldh` (lactate), `gadB` (GABA), `speE` (spermidine), `murI` (D-glutamate). These
genes are near-universal across bacteria, so essentially every genome scores
positive and the resulting number carries no information. Urolithin A cannot be
screened at all — the responsible gut bacteria have not been identified.

### 3.5 Interpretation guardrail

A high value for one panel against normal values for others is mildly interesting. A
high value across all panels indicates a systematic normalization artifact, not a
finding. Build this comparison into the multi-panel report.

---

## PART 4 — PERFORMANCE REQUIREMENTS

Target: full 16M-read run in under 30 minutes on an Apple Silicon Mac Pro
(24 cores, 192 GB unified memory). The reference databases are small (thousands of
proteins), so the bottleneck is query translation and I/O, not index size.

**Install DIAMOND natively.** Bioconda's osx-64 builds run under Rosetta and cost
roughly 2x. Try, in order: Homebrew (`brew install diamond`, or the `brewsci/bio`
tap), compiling from source, then `CONDA_SUBDIR=osx-64 conda` as the fallback.
Detect and warn at runtime if running translated under Rosetta.

**Read the `.gz` files directly.** DIAMOND accepts gzipped FASTQ. This halves disk
I/O and lets the user delete the 5.2 GB of uncompressed files. Default to `.gz` when
both exist.

**DIAMOND flags:**

| Flag | Value | Reason |
|---|---|---|
| `-p` | core count | All cores |
| `-k 1` / `--max-target-seqs 1` | 1 | Best hit only; large speedup |
| `-c 1` / `--index-chunks 1` | 1 | Single in-memory index chunk; ample RAM |
| `-b` | 8–12 | Block size; roughly 6x this in GB RAM |
| `-f 6` | minimal fields | Do not request fields you will not parse |
| sensitivity | default | `--very-sensitive` is ~10x slower and unnecessary at 50–60% identity |
| `--quiet` | | Parse-friendly output |

Make threads, block size, and sensitivity CLI-configurable with these as defaults.

**Other requirements:**
- Build the DIAMOND index once and cache it. Re-runs must not rebuild.
- Cache the fetched reference FASTA and its metadata; `--refresh-refs` forces
  re-download.
- Cache per-file DIAMOND hit output; an interrupted run must resume, not restart.
- Stream the hits TSV line by line. Never load it fully into memory.
- Emit progress to stderr with elapsed time per stage. A 30-minute silent run is a
  bug report waiting to happen.
- Provide `--subsample N` to run against the first N reads for a fast smoke test.

---

## PART 5 — ARCHITECTURE

Panels must be **declarative**. Adding a gene family is a config file edit.

```
openbiota-gut-health-test/
├── README.md
├── BUILD_SPEC.md
├── pyproject.toml
├── panels/
│   ├── urda.yaml
│   ├── cutc.yaml
│   └── ...
├── openbiota/
│   ├── __init__.py
│   ├── cli.py           # argparse; orchestration only
│   ├── panels.py        # load + validate panel YAML
│   ├── references.py    # UniProt/curated fetch, FASTA build, length stats
│   ├── search.py        # DIAMOND wrapper, caching, resume
│   ├── tally.py         # parse, filter, mate-collapse, normalize
│   ├── residues.py      # position-373-style checks (stretch)
│   └── report.py        # text + JSON output
├── tests/
│   ├── test_tally.py
│   ├── test_panels.py
│   └── fixtures/
└── fastq/
```

Panel schema:

```yaml
name: urda
description: Urocanate reductase — imidazole propionate production
targets:
  - id: URDA
    source: uniprot
    query: '(protein_name:"urocanate reductase" OR gene:urdA) AND (taxonomy_id:2)'
    min_identity: 60.0
decoys:
  - id: FRDA
    source: uniprot
    query: '(protein_name:"fumarate reductase flavoprotein subunit") AND (reviewed:true)'
    min_identity: 50.0
residue_check:
  enabled: true
  anchor_accession: Q8EAP8      # S. oneidensis UrdA
  canonical_position: 373
  accepted_residues: [Y, M]
min_fragments_for_stability: 20
citation: "Vemuganti et al., Nat Commun 17:8017 (2026)"
```

`rpoB` is a built-in normalizer included in every run, not a per-panel entry.

**Engineering standards:** Python 3.11+, type hints throughout, no bare `except`,
`pathlib` over string paths, `subprocess.run(check=True)` with real error messages.
Fail loudly on a missing dependency or an empty reference fetch. Never silently
produce a number from partial data — a zero-sequence reference fetch must abort, not
yield a count of zero.

---

## PART 6 — OUTPUT

Write to `results/<sample>/`:

- `summary.txt` — human-readable, one section per panel
- `results.json` — machine-readable; every number in the text report
- `hits_<panel>_<mate>.tsv` — raw DIAMOND output
- `run.log` — versions, flags, timings, reference set sizes

Each panel section must report: fragments accepted, fragments rejected as decoys,
`rpoB` fragments, the normalized figure, a stability flag if below the minimum
fragment count, residue-check results when applicable, and the closest reference
organisms with an explicit note that best-hit reference similarity is **not**
taxonomic identification.

Include the standing caveats: this measures capacity not metabolite; no validated
reference range exists; the count is an upper bound where the residue filter is
unavailable.

---

## PART 7 — TESTING

Unit tests, no real FASTQ required, using synthetic DIAMOND TSV fixtures:

1. **Mate collapse** — 40 fragments appearing as 80 read hits must count as 40.
2. **Decoy rejection** — hits whose best match is a decoy count as rejected, never
   as target.
3. **Identity threshold** — a target hit below `min_identity` is dropped.
4. **Normalization arithmetic** — assert against a hand-computed expected value.
5. **Zero rpoB** — returns indeterminate, does not divide by zero.
6. **Panel validation** — malformed YAML fails with a clear message.
7. **Organism counts reconcile** — attribution counts must sum to the fragment total,
   not the read total. (This was a real bug in the prototype.)

Integration: `--subsample 100000` against the real files, asserting the run completes
and `rpoB` fragments exceed zero.

Include a `make validate` target documenting how to benchmark against the source
paper's public data: SRA `PRJNA1271016` (raw metagenomes), `PRJNA1469852` (MAGs).
Comparing this tool's calls against their published `urdA`-positive results is the
only honest way to characterize its false-positive rate. State in the README that
until this is done, accuracy relative to the published method is unmeasured.

---

## PART 8 — README

Write a full `README.md` covering: what the tool measures and explicitly what it does
not; both source papers with links; the method; a table comparing this pipeline to
Vemuganti et al.'s published method step by step, with the missing residue-373 filter
called out as the main weakness; where the read-level approach is defensible and why;
usage and all CLI flags; how to interpret output, including that longitudinal
comparison of one person over time is more meaningful than any single value; how to
add a panel; the excluded genes and why; and the validation gap.

Reference method for the comparison table — Vemuganti et al. used: SPAdes assembly
(`-k 21,33,55,77`), contigs under 500 bp dropped, MetaBAT2 binning, CheckM filtering
at >90% completeness and <5% contamination, dRep dereplication, GTDB-Tk taxonomy,
Prodigal gene calling, an HMM built from known ImP-producer references via MUSCLE and
hmmbuild, a threshold of 60% of the lowest-scoring model sequence, Clustal Omega
alignment, the Y/M-at-373 FAD-site check, and kallisto TPM abundance with CLR
transformation.

Note in the README that of all `urdA`-positive genomes in that study, only
*Streptococcus pasteurianus* abundance was significantly positively associated with
plasma ImP — gene carriage was a weak predictor of the metabolite even across 294
samples.

Licence: MIT.

---

## PART 9 — BUILD ORDER

1. Skeleton, `pyproject.toml`, CLI stub, panel loader, tests for the loader.
2. Reference fetch and DIAMOND database build for the `urda` panel.
3. Search wrapper with caching, resume, and progress reporting.
4. Tally and normalization, with the full unit test suite passing on fixtures.
5. Report writers.
6. End-to-end `--subsample 100000` run, then the full run.
7. README.
8. `cutC`/`cntA` panel using the Rath curated references.
9. Residue-373 check.
10. Remaining panels.

Stop after step 6 and show me the output before continuing.
