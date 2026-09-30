# Usage

Every command and flag.

---

## Commands

```
openbiota [run] [options]     screen a sample (default command)
openbiota panels              list or inspect pathways
openbiota build-db            fetch reference sets, build the search database
openbiota cohort              build percentile reference ranges
openbiota validate            measure sensitivity, specificity and accuracy
openbiota depth-check         test whether the sequencing depth was sufficient
openbiota doctor              check dependencies, hardware and inputs
```

`openbiota` with no subcommand means `openbiota run`.

---

## `openbiota run`

### Inputs

| Flag | Default | Purpose |
|---|---|---|
| `--fastq-dir DIR` | `./fastq` | Directory holding the FASTQ pair |
| `--r1 PATH`, `--r2 PATH` | — | Explicit mate paths, overriding discovery |
| `--sample NAME` | inferred | Sample name; required if the directory holds several |
| `--no-prefer-gzip` | off | Use uncompressed FASTQ when both forms exist |

Filenames like `NAME_1.fastq.gz` / `NAME_2.fastq.gz` or `NAME_R1.fq` /
`NAME_R2.fq` are recognised automatically.

`.gz` is preferred by default because DIAMOND decompresses natively, halving
disk I/O and letting you delete the uncompressed copies. Record counting and
subsampling independently prefer the plain twin when present, since for those
stages sequential reads beat inflating through zlib.

### Pathways

| Flag | Default | Purpose |
|---|---|---|
| `--panels LIST` | `all` | Pathways to **report** |
| `--db-panels LIST` | `all` | Pathways to include in the search **database** |
| `--min-fragments N` | per-panel | Override every stability threshold |
| `--no-residue-check` | off | Skip active-site checks; drops two output columns |

Keep `--db-panels` at `all`. Cross-pathway best-hit competition is a
specificity mechanism, and a wider database costs almost nothing at runtime.

### Performance

| Flag | Default | Purpose |
|---|---|---|
| `--threads`, `-p N` | all cores | DIAMOND threads |
| `--block-size`, `-b F` | `8` | Query block size in billions of letters; RAM ≈ 6× this in GB |
| `--index-chunks`, `-c N` | `1` | Single in-memory index chunk |
| `--sensitivity MODE` | `default` | `--very-sensitive` is ~10× slower and unnecessary at 50–60% identity |
| `--evalue F` | `1e-5` | Stricter than DIAMOND's `1e-3` default |
| `--subsample N` | off | First N read pairs only |
| `--diamond PATH` | `diamond` | Path to the executable |

### Caching

| Flag | Effect |
|---|---|
| `--refresh-refs` | Re-download reference sets from UniProt |
| `--rebuild-db` | Rebuild the DIAMOND database |
| `--force` | Ignore cached DIAMOND output and input profile |

Hit output is cached per mate file against a manifest covering the database
fingerprint, the exact flags, the field list, and the input's size and mtime.
`--threads` is excluded from the manifest — it changes how fast DIAMOND runs,
not what it writes, so changing the thread count never invalidates a cache.
An interrupted run resumes at the mate boundary. Output goes to a `.partial`
path and is renamed on success, so a killed run never leaves a truncated file
that a later run mistakes for complete.

The reference database is keyed by a fingerprint over every query, cap and
filter, so it rebuilds only when a pathway definition actually changes.

### Input profiling

| Flag | Effect |
|---|---|
| `--no-qc` | Skip input validation and read statistics |
| `--no-read-count` | Skip the exact record count (a full pass over the input) |
| `--qc-sample-reads N` | Reads examined for sampled metrics (default 400,000) |
| `--no-taxonomy` | Skip the community profile |

### Inputs by path, and the sample name

| Flag | Effect |
|---|---|
| `--r1 PATH`, `--r2 PATH` | Explicit mate files, when they are not in `./fastq` or not named `<sample>_1/_2` |
| `--sample NAME` | The sample name (default: inferred from the file names); it names the results directory and appears in the report header |

### Organism detection and confirmation

| Flag | Effect |
|---|---|
| `--lanes LANE ...` | Run only these expanded lanes (`metaphlan_jan26 sylph_globdb motus4 kraken_uhgg singlem_globdb`); default all that are installed |
| `--metaphlan4 PATH` | A MetaPhlAn 4 executable for the Jun23 extended catalogue |
| `--no-assembly` | Skip targeted assembly even when SingleM triggers it; the completeness gate records the skip |
| `--no-pathogens` | Skip the pathogen screen |
| `--strict-cdiff` | Report a *Clostridioides difficile* species detection as a positive finding even without the toxin genes (`tcdA`/`tcdB`), instead of filing it as carriage; the report states the mode |
| `--no-age` | Skip the estimated microbiome age |
| `--mode {participant,clinician,research}` | Report mode; `research` (the default) states that the age estimate carries no residual against an actual age |
| `--profile-validation PATH` | Labelled-cohort AUCs and specificity matrix from `openbiota validate-profiles` |
| `--panels-dir DIR` | An alternative panel definition directory |

### Model-assisted simulation

| Flag | Effect |
|---|---|
| `--simulate {auto,always,never}` | `auto` solves the synbiotic scenarios and caches every linear programme, so a regenerated report costs seconds; `never` reports readiness only (the fast path while iterating on the report: `--simulate never --no-pdf`) |
| `--lp-workers N` | Worker processes for the scenario solver (default: cores minus two, at most 16); each programme has a ten-minute limit and a batch gives stragglers up after thirty |

### Output

| Flag | Effect |
|---|---|
| `--compact` | One-screen answer only |
| `--json` | Print `results.json` to stdout |
| `--no-pdf` | Skip the PDF report |
| `--no-panel-hits` | Skip per-panel hit views |
| `--reference-ranges PATH` | Percentile ranges (default `refs/reference_ranges.json`) |
| `--validation PATH` | Measured performance (default `docs/validation_results.json`) |
| `--out DIR` | Results root (default `./results`) |
| `--refs-dir DIR` | Reference cache root (default `./refs`) |
| `--quiet`, `-q` | Suppress progress on stderr |

---

## `openbiota run-all`

Screens every mate pair in the FASTQ directory, one after another, and prints a
side-by-side comparison. Each sample gets its own directory and its own named
report, exactly as a single run would.

```bash
openbiota run-all                                   # everything in ./fastq
openbiota run-all --only SAMPLE2,A03               # a subset
openbiota run-all --skip-existing                   # resume an interrupted batch
openbiota run-all --stop-on-error                   # abort on the first failure
```

Accepts every `openbiota run` flag; they apply to all samples. Defaults to carrying
on when one sample fails, so a single bad input cannot lose a long batch — the
failure is reported in the comparison table and in `comparison.json`.

Writes `results/comparison.json` alongside the per-sample directories:

```
results/
├── comparison.json
├── SAMPLE2/
│   ├── SAMPLE2_report.pdf
│   └── ...
└── A03/
    ├── A03_report.pdf
    └── ...
```

Percentiles in the comparison are against the reference cohort, not against the
other samples in the batch.

## `openbiota prune`

Each distinct combination of pipeline version, DIAMOND version, reference set
and residue-check configuration gets its own database directory keyed by a
fingerprint, so that a past result can always be reproduced against the exact
references that produced it. They accumulate at about 54 MB each.

```bash
openbiota prune              # list what is superseded; deletes nothing
openbiota prune --yes        # delete it
openbiota prune --alignments --yes   # also drop per-sample alignment TSVs
```

The database the current configuration would use is never deleted. Alignment
TSVs are large and fully regenerable from the cached DIAMOND output, so
`--alignments` is safe if you only need the reports.

## Input handling

Inputs stay gzipped. Every consumer — DIAMOND, Bowtie2, MetaPhlAn — reads gzip
natively, and where a full pass is needed the file is streamed through the
system `gzip -dc` rather than Python's `gzip` module, which is about five times
faster (3.9 s against 19.8 s for a 594 MB mate). `pigz` was measured and is
slower: a single gzip stream cannot be inflated in parallel.

**Nothing is ever decompressed to a file.** Keeping only the `.gz` is both the
efficient and the supported arrangement; there is no plain-FASTQ twin to
manage.

Mate pairs are recognised from `NAME_1`/`NAME_2` or `NAME_R1`/`NAME_R2`, with
`.fastq`, `.fq`, `.fastq.gz` or `.fq.gz`. The sample name is whatever precedes
the mate number, so renaming a file renames its report and its output
directory.

## Read QC flags

`fastp` makes one pass over the reads before any alignment. By default it only
measures; nothing downstream changes.

| Flag | Default | Meaning |
|---|---|---|
| `--no-fastp` | off | Skip the fastp pass; the dependent QC gates read `unknown` |
| `--fastp PATH` | auto | Path to a fastp executable |
| `--trim` | off | Write adapter/polyG-trimmed reads and run every later stage on those |
| `--detect-adapters` | off | fastp's adapter auto-detection pre-pass (single-threaded and slow; paired-end overlap trimming already handles adapters) |

Results are cached by input fingerprint under `<out>/<sample>/preprocess/`.

## Profile similarity flags

Part of `openbiota run`. Skip the whole stage with `--no-profiles`.

| Flag | Default | Meaning |
|---|---|---|
| `--no-profiles` | off | Skip the species engine and profile scoring entirely |
| `--profile NAMES` | `all` | Comma-separated profiles to score |
| `--profiles-dir DIR` | `./profiles` | Where profile YAML lives |
| `--metaphlan PATH` | auto | MetaPhlAn 3 executable; also honours `$OPENBIOTA_METAPHLAN` |
| `--metaphlan4 PATH` | auto (`.venv-mpa4`) | MetaPhlAn 4 executable for the extended SGB catalogue |
| `--no-extended-catalogue` | off | Skip the MetaPhlAn 4 inventory (reported, never scored) |
| `--include-opt-in` | off | Also score profiles marked `opt_in` in their spec. No shipped profile is marked, so this changes nothing today |

## FMT donor matching

`openbiota fmt match` compares one recipient against candidate donor materials
and evaluates donor sets. It reads the native `results.json` of each sample, so
run the screen first. See [FMT donor matching](FMT.md) for the full contract —
including what the coverage index does and does not mean, and why no donor is
ever described as cleared.

```bash
openbiota fmt match --recipient results/RECIPIENT --donor A=results/DONOR_A \
  --donor B=results/DONOR_B --indication longcovid --out results/fmt_run
openbiota fmt models list          # capability and artifact status
openbiota fmt validate --run DIR   # integrity of a completed run
```
| `--taxa-dir DIR` | `./taxa` | Curated taxon groups for the microbial-groups section and carrier features |
| `--skip-host-filter` | off | Profile without removing human reads |
| `--taxonomic-cohort PATH` | `refs/taxonomic_cohort.json` | Species reference cohort |
| `--functional-cohort PATH` | `refs/reference_cohort_samples.json` | Per-sample gene-capacity values |
| `--match VARS` | `country,age_band,sex` | Reference matching variables |

Subject metadata. Supplying it narrows the comparison group and raises
confidence; omitting it costs confidence but never blocks a run.

| Flag | Effect |
|---|---|
| `--subject-country USA` | Match the reference on country (ISO3) |
| `--subject-age 46` | Match on age band |
| `--subject-sex male` | Match on sex |
| `--illness-duration-years 12` | Select the stratum for duration-dependent profiles |
| `--antibiotics-days-ago 400` | Under 90 triggers abstention |
| `--onset-date`, `--sampling-date` | Required by some profiles; absence can trigger abstention |
| `--stool-form`, `--medications` | Recorded as confounder metadata |

A stratum with fewer than 30 members falls back to the full cohort and is
flagged as unmatched. Below 200 the confidence grade takes a penalty, because
bootstrapping shows percentile estimates are still unstable there.

```bash
# gene capacity only — no MetaPhlAn needed
openbiota run --no-profiles

# everything, with subject metadata
openbiota run --subject-country USA --subject-age 46 --subject-sex male \
         --antibiotics-days-ago 400 --stool-form "Bristol 4"

# one profile, unmatched reference
openbiota run --profile mecfs --match ""
```

## `openbiota taxonomic-cohort`

Builds the species reference cohort from curatedMetagenomicData. Downloads
uniformly processed profiles plus curated subject metadata, and assembles a
reference distribution. One snapshot only, so one profiler version.

```bash
openbiota taxonomic-cohort --stability
```

| Flag | Default | Meaning |
|---|---|---|
| `--snapshot` | `2021-10-14` | cMD snapshot; determines the profiler version |
| `--condition` | `control` | cMD `study_condition` to select |
| `--min-reads` | 1,000,000 | Minimum sequencing depth per reference sample |
| `--max-samples` | none | Cap the cohort size |
| `--stability` | off | Also compute the cohort-size stability curve |

About 8 minutes and 15 MB of download for 3,027 samples. Writes
`refs/taxonomic_cohort.json` plus a manifest recording cohort composition,
sample count, profiler version and a content hash.

## `openbiota build-host-index`

One-time. Downloads GRCh38 and PhiX and builds the Bowtie2 index used to
remove human reads before taxonomic profiling. Roughly an hour, mostly index
construction. Until it exists, `openbiota run` warns and profiles without host
filtering.

## `openbiota validate-profiles`

Scores every profile against public cohorts that carry disease labels, and
prints the cross-profile specificity matrix.

```bash
openbiota validate-profiles --conditions CRC,IBD,T2D
```

The reference distribution is built from the control arm, never from the cases
being scored. Leave-one-study-out AUC is reported alongside the pooled figure,
because pooling cohorts and splitting randomly leaks recruitment-batch signal.
Writes `docs/profile_validation.json`, which the PDF reads.

## `openbiota profiles`

```bash
openbiota profiles                  # list
openbiota profiles --show mecfs     # one profile in full, with all weight factors
openbiota profiles --show mecfs --json
```

## `openbiota cohort`

Builds percentile ranges from a public ENA study.

| Flag | Default | Purpose |
|---|---|---|
| `--study ACC` | `PRJNA1271016` | ENA study accession |
| `--samples N` | 30 | Runs to use |
| `--reads N` | 600,000 | Read pairs per cohort sample |
| `--min-rpob N` | 300 | Exclude samples below this denominator depth |
| `--workers N` | 10 | Parallel downloads |
| `--keep-fastq` | off | Do not delete downloaded prefixes |
| `--out PATH` | `refs/reference_ranges.json` | Where to write the ranges |

Uses `fastq-dump -X` when the SRA toolkit is installed, falling back to HTTP
range requests against the ENA gzip files.

---

## `openbiota validate`

Measures performance against known ground truth.

| Flag | Default | Purpose |
|---|---|---|
| `--reads N` | 4,000,000 | Read pairs per synthetic community |
| `--error-rate F` | 0.002 | Simulated substitution rate |
| `--seed N` | 20260904 | RNG seed; the whole run is deterministic |
| `--refresh` | off | Re-download reference genomes |
| `--out PATH` | `docs/validation_results.json` | Where to write results |

---

## `openbiota depth-check`

Screens the same sample at increasing depths and reports whether the values
have stopped changing — answering "would deeper sequencing change this?" with a
measurement instead of an assertion. The result is written to
`results/<sample>/depth_check.json` and appears in the PDF report.

| Flag | Default | Purpose |
|---|---|---|
| `--fractions LIST` | `0.0625,0.125,0.25,0.5,1.0` | Fractions of full depth to screen |
| `--fastq-dir DIR` | `./fastq` | Input directory |
| `--out PATH` | `results/<sample>/depth_check.json` | Where to write the result |

Full depth is always included. Cached DIAMOND output is reused, so re-running
after a full screen only pays for the shallower depths.

---

## Make targets

```
make setup       create .venv and install
make doctor      check dependencies, hardware and inputs
make panels      list pathways
make build-db    fetch reference sets and build the database
make cohort      build percentile reference ranges
make validate    measure performance
make depth-check test whether the sequencing depth was sufficient
make all         build-db, cohort, validate, then run
make smoke       fast end-to-end test on 100k read pairs
make run         full screen                        [default]
make compact     full screen, one-screen output
make json        full screen, JSON to stdout
make test        unit tests
make lint        ruff
make clean       remove results/
make clean-all   remove results/, refs/ and .venv/
```

Variables: `THREADS`, `FASTQ_DIR`, `RESULTS`, `SUBSAMPLE`, `COHORT_SAMPLES`.

```bash
make run THREADS=8
make cohort COHORT_SAMPLES=50
```

---

## Performance

Measured on a 16-core / 32-thread Xeon W-3245, 192 GB RAM, DIAMOND 2.2.6 built
natively, 8,029,909 read pairs / 2.38 Gbp:

| Stage | Cold | Cached |
|---|---|---|
| Input validation, incl. exact record count | 20 s | < 0.1 s |
| Reference database (36 UniProt fetches + index) | ~16 min | 0.2 s |
| Translated search, both mates | 6 min 11 s | < 0.1 s |
| Tally and normalisation | 0.7 s | 0.7 s |
| PDF report | ~1 s | ~1 s |
| **Full run, warm database** | **~6 min** | **~2 s** |

One-off setup: `build-db` ~16 min, `cohort` ~10 min, `validate` ~5 min.

Progress goes to stderr with elapsed time per stage, plus a heartbeat during
the long DIAMOND stages.

### Tuning

- `--threads` defaults to all cores; DIAMOND scales well to 32.
- `--block-size 8` loads the whole query file in one block. Reduce it on a
  memory-constrained machine — RAM is roughly 6× the value in GB.
- Leave `--sensitivity` at `default`. At the 50–60% identity these searches
  operate at, `--very-sensitive` costs ~10× for no useful gain.
- Deleting uncompressed FASTQ is safe. DIAMOND reads `.gz` natively; the only
  cost is a slower exact record count, which is cached anyway.
