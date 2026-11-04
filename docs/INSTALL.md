# Installation

One script installs everything. This page is what it does, what it needs,
what to do when it stops, and how to do any part of it by hand.

```bash
git clone https://github.com/openbiota-org/openbiota-gut-health-test
cd openbiota-gut-health-test
./install.sh            # or: make install
```

If you have never used a terminal, start with the [Quickstart](QUICKSTART.md),
which explains the three lines above from the beginning.

---

## Requirements

| | Standard install | Full install (`--full`) |
|---|---|---|
| Operating system | macOS 13 or newer (Intel or Apple silicon), or Linux (Debian, Ubuntu, Fedora; anything else via the conda fallback) | the same; the full lanes are built and tested on macOS, best-effort on Linux |
| Disk | ~20 GB (~10 GB with `--skip-host-index`) | ~1.5 TB |
| Memory | 8 GB; 16 GB or more is faster | 32 GB (the Kraken 2 lane loads a 19 GB database) |
| Time | about an hour, mostly downloads | most of a day |
| Network | once, for the reference downloads; runs are offline afterwards | a fast connection; files over 1 GB are fetched with parallel range requests |
| Privileges | your own user. The installer asks for your password once, to install the command-line tools through Homebrew or `apt`/`dnf`. It never needs root itself and refuses to run as root on macOS. | the same |

Nothing is installed outside the project folder except those command-line
tools, a link named `openbiota` in `~/.local/bin`, and one line on your
shell's start-up file that puts that folder on your `PATH` (`--no-path`
skips both). Deleting the project folder removes everything else.

---

## What the installer does

`./install.sh` runs eleven steps, each of which first checks whether its
work is already done and skips it if so. That is what makes it safe to run
again after an interruption, and what makes `./install.sh --list` a true
report of the state of the machine.

| # | Step | What it produces | Skipped when |
|---|---|---|---|
| 1 | `preflight` | nothing; it checks the OS, free disk, memory, network and write access, and on a Mac the Xcode command-line tools and Homebrew | never |
| 2 | `system` | `git`, `curl`, Python 3.11+, `diamond`, `bowtie2`, `samtools`, `fastp` (and `pigz` if available) | all are on the `PATH` |
| 3 | `python` | `.venv/` with `openbiota` installed in editable mode | `.venv/bin/openbiota` imports |
| 4 | `launcher` | `bin/openbiota` executable; `~/.local/bin/openbiota` linked to it; `~/.local/bin` on the `PATH` | the link and the `PATH` line exist |
| 5 | `profiler` | `.venv-mpa3/` with MetaPhlAn 3.1.0 — the scoring lane, pinned because the reference cohort was profiled with it | it reports that version |
| 6 | `data` | the two bundles under `data/` unpacked into `refs/`: the project's own reference data and the UniProt responses behind the protein database (below) | every item each bundle carries is present |
| 7 | `metaphlan_db` | `refs/metaphlan_db/` — the `mpa_v31_CHOCOPhlAn_201901` marker database: a 0.4 GB download from the MetaPhlAn servers, indexed by Bowtie2 to 3.5 GB (20–40 min of CPU) | MetaPhlAn's own readiness check passes |
| 8 | `diamond_db` | `refs/db/<fingerprint>/` — every panel's reference proteins, from the cached UniProt responses, indexed by DIAMOND (about two minutes). The fingerprint includes the DIAMOND version, so this is always built on the machine that will use it | the fingerprinted index exists |
| 9 | `models` | MICOM, an open LP solver and the AGORA2 species-model library under `refs/micom/` (~300 MB), for the synbiotic scenarios | the library and the pinned fixtures are present |
| 10 | `host_index` | `refs/host/` — GRCh38 + PhiX as a Bowtie2 index. Downloads the human genome (0.9 GB) and indexes it: about an hour on a desktop, longer on a laptop | the index is complete, or `--skip-host-index` |
| 11 | `check` | nothing; it runs `openbiota doctor`, imports the report modules and re-asks every readiness question | never |

Readiness is decided by the pipeline, not guessed by the script:
`python -m openbiota.installcheck --all` prints the same table the installer
consults, and `openbiota doctor` the same dependency check.

Everything each step prints goes to `install.log` in the project folder; the
terminal shows the progress bar, a live line with the elapsed time and the
step's latest message, and — if a step fails — the last lines of the log.

### How the tools are installed

The installer tries, in order:

1. **The system's package manager.** Homebrew on macOS (offering to install
   it if it is missing); `apt-get` on Debian and Ubuntu; `dnf` on Fedora.
   Everything the standard install needs is in Homebrew core and in the
   Debian/Ubuntu *universe* repositories (`diamond-aligner`, `bowtie2`,
   `samtools`, `fastp`).
2. **conda-forge / bioconda**, for anything the package manager could not
   supply — a Linux without `sudo`, a distribution without a `fastp`
   package, a Python older than 3.11. The installer fetches `micromamba`
   (a single static binary) into `vendor/bin/` and creates `vendor/conda/`
   with the missing tools. Nothing is written outside the project.
3. **A clear message** naming what is missing and how to install it, if
   both routes fail.

`bin/openbiota` and the `Makefile` put `vendor/bin` and `vendor/conda/bin`
first on the `PATH`, so a tool that came from conda is found exactly as it
was during the install.

### Options

```
./install.sh --help
```

| Option | Effect |
|---|---|
| `--full` | also run the full-install steps below |
| `--skip-host-index` | leave out step 10 (saves ~7 GB and an hour; the report then cannot state the human DNA fraction) |
| `--no-path` | do not link `openbiota` into `~/.local/bin` or touch the shell start-up file |
| `--only STEP[,STEP]` / `--skip STEP[,STEP]` | run a subset of the steps |
| `--threads N` | CPU threads for the builds (default: all) |
| `--force` | redo steps even when they look complete |
| `--yes` | do not stop to ask before installing system packages (for scripts and CI) |
| `--offline` | fail fast instead of downloading; verifies what is already there |
| `--dry-run` | print what would run, run nothing |
| `--list` | show every step and whether it is done |
| `--no-color` | plain output (also honours `NO_COLOR`) |

`make install` runs the script; `make install INSTALL_ARGS="--skip-host-index"`
passes options through; `make install-full` is `--full`.

---

## The reference data that travels with the repository

`refs/` on a complete installation is about 1.3 TB, and almost all of it is
someone else's — GlobDB, GTDB, the MetaPhlAn indexes, Kraken 2, AGORA2 —
with a home of its own that the installer fetches from. The project never
hosts those, so nobody pays to serve them.

Two bundles are committed instead, under `data/`:

| Bundle | Size | What it is |
|---|---|---|
| `derived-refs` | 50 MB | What this project computed itself: the percentile ranges every pathway is scored against, the scoring cohort matrix, the trained microbiome-age models, the GTDB crosswalks, the reconciliation locks, the resolution lanes' reference sequences and the BacDive cache. Rebuilding them is not a download — it is hours of compute, R and ExperimentHub for the cohort, an API account for BacDive. |
| `uniprot-cache` | 200 MB | The UniProt REST responses behind the protein search database: one gzipped table per reference entry, 244 of them, each with a sidecar recording its query. UniProt answers these queries at roughly ninety seconds each, which would make step 8 a six-hour fetch on a fresh machine; with the responses present it is a two-minute local build. UniProt's data is [CC BY 4.0](https://www.uniprot.org/help/license). |

Each bundle is a plain `tar.gz` — nothing extra to install — cut into parts
under GitHub's 50 MB per-file warning, with a manifest holding every part's
size and SHA-256 and a line on every item:

```
data/derived-refs.tar.gz.part00 … part01      data/derived-refs.manifest.json
data/uniprot-cache.tar.gz.part00 … part04     data/uniprot-cache.manifest.json
```

Step 6 verifies every part against its manifest, joins them and unpacks into
`refs/`, leaving any file that is already there alone. The script behind it:

```bash
scripts/derived_refs.py list                   # what is in them, and why
scripts/derived_refs.py verify [BUNDLE]        # do the parts match the manifest?
scripts/derived_refs.py check  [BUNDLE]        # is the content present under refs/?
scripts/derived_refs.py restore [BUNDLE]       # unpack into refs/  (--force to overwrite)
scripts/derived_refs.py pack [BUNDLE]          # maintainers: rebuild from refs/
```

The UniProt bundle is derived from the panel files: `pack uniprot-cache`
asks the pipeline which cached responses the current panels need and packs
exactly those, and a test fails when the panels change and the bundle has
not been repacked.

---

## When something goes wrong

The installer stops at the first problem, prints a red ✗ with the reason and
what to do, and the last 25 lines of `install.log`. Fix the cause and run
`./install.sh` again; finished steps are kept. The situations it knows:

| It says | What it means |
|---|---|
| *the Xcode command line tools are not installed* | macOS needs them for `git` and compilers. The installer offers to open the dialog; finish it, then re-run. |
| *Homebrew is not installed* | The installer offers to run Homebrew's installer, which asks for your password. Or install it from [brew.sh](https://brew.sh) and re-run. |
| *no sudo: cannot use apt* | On Linux without `sudo`, the tools come from conda instead, automatically. Nothing to do. |
| *only N GB free … needs about M GB* | Free space on the drive holding the project folder, or clone onto a larger drive. `--skip-host-index` needs 10 GB rather than 20. |
| *cannot reach the internet* | The install downloads references. Check the connection; a corporate proxy needs `https_proxy` set. `--offline` verifies what is already there without downloading. |
| *another ./install.sh (pid N) is running* | Two installers building the same index would corrupt it. Wait for the other one, or delete `.install.lock` if it is stale. |
| *python -m venv failed* | On Debian/Ubuntu: `sudo apt-get install python3-venv`, then re-run. |
| *MetaPhlAn could not install its database* | Almost always a dropped connection. Re-run; the download resumes. |
| *openbiota build-db failed* | With the UniProt bundle unpacked this step needs no network; see `install.log` for DIAMOND's message. Without it, UniProt's REST service was unreachable or rate-limited: re-run later, everything fetched so far is cached under `refs/cache/`. |
| *could not install X by any route* | Install it yourself with your package manager or conda, make sure it is on the `PATH`, and re-run. |
| *do not run the installer as root on macOS* | Homebrew refuses root. Run it as your own user. |

Three commands answer most other questions:

```bash
./install.sh --list                         # which steps are done
bin/openbiota doctor                        # tools, hardware, inputs
.venv/bin/python -m openbiota.installcheck --all   # every readiness question
```

`install.log` has everything every step printed, with timestamps.

---

## The full install

`./install.sh --full` continues past the standard eleven steps with the
optional detection lanes. Each is a `Makefile` target, so every pinned
version lives in one place (`make help` lists them), and each is resumable.
These are large: budget about 1.5 TB and most of a day.

| Step | `make` target | Adds | Size |
|---|---|---|---|
| `profiler4` | `profilers` | MetaPhlAn 4.1.1 in `.venv-mpa4` — the extended SGB catalogue | the Jun23 index is fetched on first use, ~25 GB |
| `vendor_tools` | `vendor-tools` | BLAST+, AMRFinderPlus and StxTyper under `vendor/` | 2 GB |
| `strain_tools` | `strain-tools` | ECTyper, StrainGE, MIDAS, SameStr, TRACS, inStrain, Kleborate in their own environments | 5 GB |
| `refs_strain` | `refs-strain` | the strain and typing reference panels, lock-verified | 1 GB |
| `pathogens` | `tools/build_pathogen_bundle.py` | the 485-target pathogen screen: 377 NCBI genomes, k-mer ownership against human and food-animal decoys, a competitive Bowtie2 index | 130 GB, the heaviest build |
| `substrate_refs` | `substrate-refs` | the dbCAN snapshot behind the substrate panel | 2 GB |
| `tools_expanded` | `tools-expanded` | sylph, mOTUs 4, Kraken 2, Bracken, SingleM, SPAdes, MetaPhlAn 4.2 | 3 GB |
| `refs_expanded` | `refs-expanded` | the expanded references: MetaPhlAn Jan26, mOTUs, Kraken 2 UHGG, SingleM, GTDB, GlobDB r232 | ~800 GB |
| `expansion` | `expansion` | GTDB crosswalks and the supplement reconciliation | minutes |
| `genome_lane` | `genome-lane` | sylph and the GTDB R232 species sketch | 50 GB |
| `mycobiome_refs` | `mycobiome-refs` | the fungal references: EukDetect2, FungiGutDB, RefSeq fungal genomes, strain panels | 230 GB |

Until a lane's references are present the report marks it *not assessed* and
says why — never as a negative result. [Organism detection](EXPANDED_DETECTION.md)
describes the lanes and [`refs/MANIFEST.md`](https://github.com/openbiota-org/openbiota-gut-health-test/blob/main/refs/MANIFEST.md)
lists what each one puts on disk.

!!! note "Platform"
    The full-install lanes were built and validated on macOS. The `Makefile`
    fetches macOS binaries for BLAST+, SPAdes and the NCBI datasets CLI;
    on Linux, install those with your package manager or conda and the
    targets will find them on the `PATH`.

---

## Doing it by hand

Every step of the installer is an ordinary command, and the `Makefile` has
had a target for each since before the script existed. For a machine that
is managed some other way, or to redo one piece:

```bash
# tools
brew install diamond bowtie2 samtools fastp python@3.13        # macOS
sudo apt-get install diamond-aligner bowtie2 samtools fastp \
     python3 python3-venv python3-pip                           # Debian / Ubuntu
micromamba create -p vendor/conda -c conda-forge -c bioconda \
     python=3.12 diamond bowtie2 samtools fastp                 # anywhere

# environments
make setup                                   # .venv with openbiota (editable)
python3 -m venv .venv-mpa3 && .venv-mpa3/bin/pip install 'metaphlan==3.1.0'

# reference data
scripts/derived_refs.py restore              # the committed bundles → refs/
.venv-mpa3/bin/metaphlan --install --index mpa_v31_CHOCOPhlAn_201901 \
     --bowtie2db refs/metaphlan_db           # 0.4 GB download, 3.5 GB built
make build-db                                # cached UniProt tables → DIAMOND, ~2 min
make extension-refs                          # MICOM + AGORA2, ~300 MB
make host-index                              # GRCh38 + PhiX → Bowtie2, ~1 h

# verify
bin/openbiota doctor
make test                                    # unit tests; no network or FASTQ
```

### Apple silicon

Install DIAMOND **natively**, which `brew install diamond` does. Bioconda's
`osx-64` package runs under Rosetta and costs roughly 2×; `openbiota doctor`
reports the binary's architecture and warns when it does not match the host.

### DIAMOND from source

If no package is available:

```bash
git clone https://github.com/bbuchfink/diamond
cd diamond && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j"$(nproc)" && sudo make install
```

---

## Hardware and runtime

Measured on a 32-core workstation with 192 GB of memory and an NVMe drive,
one eight-million-pair sample at a time:

| stage | time | notes |
|---|---|---|
| read search, profiling, patterns, age, attention, solutions | 10–20 min cold, ~2 min cached | 30 threads; the DIAMOND search scales with cores |
| model-assisted simulation | 35–40 min cold, seconds cached | measured on a fresh clone: ~4 min to build each community model from the AGORA2 library, then 82 linear programmes per arm, single-threaded. Cached on the community, so a new *sample* pays it again; the same sample re-run does not |
| the expanded detection lanes (full install) | +35 min | most of it MetaPhlAn Jan26, SingleM and the two Kraken panels; every lane is cached on the reads |
| competitive confirmation (full install) | +15–20 min | one Bowtie2 alignment against the candidates and their relatives |
| strain placement (full install) | +5–40 min | per clade, cached |
| PDF | ~15 s | |

So the first run of a sample on the standard install is about an hour on
this machine and longer on a laptop, almost all of it the simulation; a
repeat run is about five minutes. The standard install's run ends with the
completeness census listing the full-install lanes as *not installed*, in
red, and exits non-zero — see the [Quickstart](QUICKSTART.md#step-3-run).

Every run records its own peak memory in `results.json` (`run.resources`):
the pipeline process peaks at 12–14 GB and the largest child process (the
Kraken 2 lane with the UHGG database loaded) at ~19 GB, so 32 GB is a
comfortable minimum for the full detection and 8 GB is enough for the
standard report. Results take 4–9 GB per sample; `openbiota prune` removes
the regenerable intermediates.

---

## Input data

One folder per sample under `fastq/`, the two mates named by the sample:

```
fastq/
└── my-sample/
    ├── my-sample_1.fastq.gz
    └── my-sample_2.fastq.gz
```

`NAME_R1` / `NAME_R2` naming works too, and files placed directly in `fastq/`
are found as well. No preprocessing is needed or wanted — no trimming, no
adapter removal, no host filtering. See [Method](METHOD.md) for why host
filtering is unnecessary for the pathway half rather than merely optional.
Nothing under `fastq/` is tracked by git.
