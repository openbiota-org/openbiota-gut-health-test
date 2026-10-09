<p class="eyebrow">Get started</p>

# Quickstart

From a fresh clone to your first Gut Health Test Metagenomic Report. The short
version:

```bash
git clone https://github.com/openbiota-org/openbiota-gut-health-test
cd openbiota-gut-health-test
brew install diamond bowtie2 fastp   # or apt / conda
make setup && make doctor            # Python env, dependency check
make profilers build-db host-index   # one-off reference builds
cp ~/SAMPLE_1.fastq.gz ~/SAMPLE_2.fastq.gz fastq/
make run                             # → results/SAMPLE/SAMPLE_report.pdf
```

The rest of this page explains each step, what it costs, and where the report
ends up. To see what you are working towards first, two finished reports from
real samples are in
[`sample-reports/`](https://github.com/openbiota-org/openbiota-gut-health-test/tree/main/sample-reports):
one community the index places as broadly disturbed, one in the healthy range.

## 1. What you need

| | Minimum | Comfortable |
|---|---|---|
| Operating system | macOS 13+ or Linux | either |
| Python | 3.11 | 3.13 |
| CPU | 8 cores | 16–32 cores (the search and alignment stages scale linearly) |
| Memory | 16 GB | 32 GB or more (the MetaPhlAn 4 index alone is ~20 GB) |
| Disk for reference data | 60 GB | 100 GB — see [Reference data](#3-reference-data) |
| Command-line tools | [DIAMOND](https://github.com/bbuchfink/diamond), [Bowtie2](https://github.com/BenLangmead/bowtie2) | plus [fastp](https://github.com/OpenGene/fastp) for read-quality gates |

And a **shotgun-sequenced stool sample as paired FASTQ files** — the raw reads,
not a 16S amplicon run and not a PDF. If you tested with
[Tiny Health](https://www.tinyhealth.com), their support team will send download
links for your raw FASTQ files on request. Any Illumina paired-end stool
metagenome with roughly five million read pairs or more works; the report states
its own depth and flags the readings that depth limits.

## 2. Install the software

```bash
git clone https://github.com/openbiota-org/openbiota-gut-health-test
cd openbiota-gut-health-test
brew install diamond bowtie2 fastp   # Linux: apt or conda, see Installation
make setup                           # creates .venv, installs `openbiota`
make doctor                          # checks tools, hardware, ./fastq
```

`make setup` is the only Python step; it installs the package in editable mode
with its three runtime dependencies (PyYAML, NumPy, ReportLab). Everything
`openbiota` needs beyond that is a command-line tool found on `PATH`, and
`doctor` tells you which are missing. See [Installation](INSTALL.md) for
platform notes, Apple Silicon, and building DIAMOND from source.

## 3. Reference data

The report compares your sample against published references, and those
references are built once and cached under `refs/` (git-ignored, so a fresh
clone starts empty). Each build is resumable and reused by every later run.

| Build | Command | What it fetches / builds | One-off cost |
|---|---|---|---|
| MetaPhlAn 3 + 4 environments | `make profilers` | Two Python environments (the two versions cannot share one); their marker databases (~28 GB) download on first use | 20–40 min, mostly download |
| Metabolite pathway proteins | `make build-db` | Curated protein references from UniProt → DIAMOND database (184,607 proteins) | ~16 min |
| Human-read filter | `make host-index` | GRCh38 + PhiX → Bowtie2 index (~6 GB) | one to a few hours of CPU |
| Species reference cohort | `make taxonomic-cohort` | 3,027 healthy-control stool profiles from curatedMetagenomicData | minutes |
| Metabolite reference cohort | `make cohort` | Percentile ranges from public stool metagenomes, sampled uniformly from each run | several hours and ~140 GB of temporary FASTQ; see [Validation](VALIDATION.md#the-reference-cohort) |
| Microbiome-age model | `openbiota age-model train` | Trains and freezes the age ensemble on curatedMetagenomicData | an hour or more; needs the TRPCA environment described in [Method](METHOD.md) |

Every build is optional in the sense that `openbiota run` degrades gracefully:
a missing reference makes the report *omit* that section with a note saying
why, rather than fail. `make profilers build-db host-index` is the set that
gives you the community, pathway and pathogen sections; the cohorts and the age
model add percentiles, disease-pattern scores and the age page.

!!! info "The pathogen reference bundle"
    The 485-target pathogen screen aligns against a curated bundle of reference
    genomes and k-mer indexes under `refs/pathogens/` (about 40 GB). The bundle
    is built by the maintainers from NCBI assemblies and is not yet distributed
    with the repository. Until it is installed, the report marks every pathogen
    target as *not assessed*, states the reason, and never lets an empty screen
    read as "nothing found".

## 4. Add your sample

Put the paired FASTQ files in `./fastq`, gzipped or not:

```
fastq/
├── SAMPLE_1.fastq.gz
└── SAMPLE_2.fastq.gz
```

`SAMPLE_R1` / `SAMPLE_R2` naming also works. Do **not** trim, filter or
subsample the reads first — the pipeline does its own quality assessment and
host filtering, and it needs every read to normalise correctly. Several samples
can sit side by side; `make run-all` screens them all and adds a comparison.

## 5. Generate the report

```bash
make run     # full run on ./fastq → results/<sample>/
make smoke   # same pipeline, first 100,000 read pairs only
```

A full run on an eight-million-pair sample takes about 20 minutes cold on a
32-core workstation and about 5 minutes when the cached stages are reused;
the last line of `run.log` records the total. The output directory holds:

| File | What it is |
|---|---|
| `<sample>_report.pdf` | The report. Page 1 is the summary; every row on it links to its detail card. |
| `results.json` | Every number in the PDF, machine-readable, with provenance. |
| `summary.txt` | The one-screen text summary (`--compact` prints the same thing). |
| `run.log` | The full log of the run, ending with its total wall-clock time. |
| `preprocess/`, `alignments/`, `taxonomy/`, `pathogens/` | Per-stage outputs and caches, so a re-run only redoes what changed. |

Open the PDF and start at page 1. [Reading the report](OUTPUT.md) explains the
seven-band scale every reading is placed on, and what each section means.

## 6. Useful next commands

```bash
openbiota run --help             # every flag: threads, panels, subject data
openbiota run --subject-age 42   # the age model also reports its error
openbiota run --compact          # one-screen answer in the terminal
openbiota panels                 # the 25 metabolite pathways and their genes
openbiota profiles               # the 33 disease-pattern profiles
openbiota depth-check            # was the sequencing deep enough?
make test                        # unit tests; no network or FASTQ needed
```

The full command reference is in [Usage](USAGE.md).
