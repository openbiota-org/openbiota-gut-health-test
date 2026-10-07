# Installation

---

## Requirements

| | |
|---|---|
| Python | 3.11 or newer |
| DIAMOND | 2.x, installed natively |
| Disk | ~25 GB for the core reference data, cohort downloads and results; ~1.3 TB more for the full organism detection (see below) |
| RAM | 8 GB minimum; `--block-size 8` wants ~48 GB, lower it if you have less |
| Part B | Bowtie2, MetaPhlAn 3.1 in `.venv-mpa3` (species scoring lane) |
| Optional | fastp (read QC and trimming); MetaPhlAn 4.1 in `.venv-mpa4` (extended SGB catalogue, ~25 GB database); SRA toolkit (`openbiota cohort`); reportlab + matplotlib (PDF report) |

---

## macOS

```bash
brew install diamond
brew install sratoolkit          # optional, for openbiota cohort
make setup
make doctor
```

### Apple Silicon

Install DIAMOND **natively**. Bioconda's `osx-64` package runs under Rosetta
and costs roughly 2×.

```bash
brew install diamond             # arm64 native
# or, if the formula is unavailable:
brew install brewsci/bio/diamond
```

`openbiota doctor` reports the DIAMOND binary's architecture and warns when it does
not match the host. The runtime also warns if the Python process itself is
running translated.

Last resort only:

```bash
CONDA_SUBDIR=osx-64 conda install -c bioconda diamond
```

---

## Linux

```bash
# Debian / Ubuntu
sudo apt install diamond-aligner
# or conda
conda install -c bioconda diamond sra-tools

python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/openbiota doctor
```

---

## From source

If no package is available:

```bash
git clone https://github.com/bbuchfink/diamond
cd diamond && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j"$(nproc)" && sudo make install
```

---

## Python package

```bash
make setup                                  # creates .venv, installs editable
# or manually
python3.13 -m venv .venv
.venv/bin/pip install -e '.[dev]'
```

Dependencies are deliberately minimal: **PyYAML** for the pathway files, plus
**reportlab** and **matplotlib** for the PDF. Everything else is standard
library. `pytest` and `ruff` for development.

If reportlab is missing, the run still completes and simply skips the PDF.

---

## First-run setup

Three one-off builds, all cached afterwards:

```bash
make build-db      # ~16 min — fetches reference protein sets from UniProt
make cohort        # ~10 min — builds percentile ranges from public metagenomes
make validate      # ~5 min  — measures sensitivity, specificity, accuracy
```

Or all of it plus a screen of `./fastq`:

```bash
make all
```

`build-db` needs network access to `rest.uniprot.org`. `cohort` needs access to
NCBI SRA or `ftp.sra.ebi.ac.uk`. `validate` needs
`api.ncbi.nlm.nih.gov/datasets`. Once cached, runs are fully offline.

### The two profiler environments

MetaPhlAn 3 and 4 cannot share a Python environment, so each gets its own,
found automatically by `openbiota run`:

```bash
python3 -m venv .venv-mpa3 && .venv-mpa3/bin/pip install 'metaphlan==3.1.0'
python3 -m venv .venv-mpa4 && .venv-mpa4/bin/pip install 'metaphlan>=4.1'
.venv-mpa4/bin/metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
    --bowtie2db refs/metaphlan4_db        # ~25 GB, one-off
```

Version 3 is the *scoring* lane — the reference cohort and the GMWI2 index were
built with it. Version 4 is the *extended catalogue*: it names roughly twice as
many organisms and is reported in full, but is never scored against a cohort
profiled with a different catalogue. If it is not installed the report simply
omits that section; `--no-extended-catalogue` skips it explicitly.

`fastp` comes from Homebrew (`brew install fastp`) or conda. Without it the
QC gates that depend on it (duplication, Q30, adapter content, insert size)
read `unknown` and the run continues.

### The full organism detection

The report's organism list pools nine detection methods (MetaPhlAn 4 Jan26,
GlobDB r232 by sylph, mOTUs 4, Kraken 2 against UHGG and against a locally
built gut rescue panel, SingleM, targeted assembly, competitive read mapping,
StrainPhlAn 4). Their tools and references are installed and locked in four
steps, all resumable:

```bash
make tools-expanded       # every pinned tool (sylph, mOTUs, Kraken 2, SingleM, metaSPAdes,
                          # bowtie2, samtools, aria2, seqkit, MetaPhlAn 4.2 in .venv-mpa42)
make refs-expanded        # every reference release, checksummed and locked (~200 GB)
make refs-globdb-genomes  # the GlobDB r232 genome archive, 292 GB compressed, ~1 TB unpacked
make expansion            # crosswalks to GTDB R232, supplement reconciliation, the rescue panel
```

Budget about 1.3 TB of disk for all of it and a fast connection for the
first fetch (files over 1 GB are downloaded with parallel range requests).
Every tool version and every reference file is pinned by checksum under
`refs/expanded/`, and `make refs-expanded` never downloads a file whose lock
already verifies. `refs/MANIFEST.md` lists what is on disk and the command
that rebuilds each part. Without these installed a run still produces the
report from the core lanes and states which methods did not run.
[docs/EXPANDED_DETECTION.md](EXPANDED_DETECTION.md) describes the lanes.

### Hardware and runtime

Measured on a 32-core Apple Silicon workstation with 192 GB of memory and
an NVMe drive, one 8-million-pair sample at a time:

| stage | time | notes |
|---|---|---|
| core report (Part A and B, MetaPhlAn 3 and 4) | ~35 min | 30 threads |
| the nine detection lanes | +35 min | most of it MetaPhlAn Jan26, SingleM and the two Kraken panels; every lane is cached on the reads |
| competitive confirmation | +15-20 min | one Bowtie2 alignment against the candidates and their relatives; cached on reads, genomes and operating points |
| strain placement | +5-40 min | per clade, cached; a first run on a new sample places every eligible organism |
| model-assisted simulation | +45-60 min | linear programmes on the community model; cached on the community |
| PDF | ~15 s | |

Two samples run comfortably side by side on that machine. Every run records
its own peak memory in `results.json` (`run.resources`): on these samples
the pipeline process peaks at 12-14 GB and the largest child process
(the Kraken 2 lane with the UHGG database loaded) at ~19 GB, so 32 GB of
memory is a comfortable minimum for the full detection and 8 GB is enough
for the core report. Disk: ~25 GB for the core references, ~1.3 TB with the
full organism detection installed, and 4-9 GB of results per sample
(`openbiota prune` removes the regenerable intermediates).

---

## Verifying the install

```bash
make doctor
```

```
openbiota 1.1.0
python              3.13.2
platform            Darwin/x86_64, Intel(R) Xeon(R) W-3245 CPU @ 3.20GHz
cpus                32 logical
memory              192.0GB
diamond             diamond version 2.2.6
                    /usr/local/bin/diamond [x86_64]
PyYAML              6.0.3
panels              9 panels, 36 reference entries
inputs              sample 'A02'
                    R1: A02_1.fastq.gz (566.2MB, gzip)
                    R2: A02_2.fastq.gz (584.5MB, gzip)

OK
```

Then:

```bash
make test          # unit tests, no network or FASTQ needed
make smoke         # 17 s end-to-end on 100k read pairs
```

---

## Input data

Paired FASTQ in `./fastq`, gzipped or not:

```
fastq/
├── SAMPLE_1.fastq.gz
└── SAMPLE_2.fastq.gz
```

`NAME_R1` / `NAME_R2` naming works too. No preprocessing is needed or wanted —
no trimming, no adapter removal, no host filtering. See
[METHOD.md](METHOD.md) for why host filtering is unnecessary rather than
merely optional.
