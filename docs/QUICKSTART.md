<p class="eyebrow">Get started</p>

# Quickstart

From nothing to your first Gut Health Test report in three steps: **install**
(one command), **add your sample** (copy two files), **run** (one command).
This page assumes you have never used a terminal. If you have, the whole
thing is:

```bash
git clone https://github.com/openbiota-org/openbiota-gut-health-test
cd openbiota-gut-health-test
./install.sh                          # ~1 hour, mostly downloads; safe to re-run
# copy SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz into fastq/SAMPLE/
bin/openbiota run                     # → results/SAMPLE/SAMPLE_report.pdf
```

To see what you are working towards, two finished reports from real samples
are in
[`sample-reports/`](https://github.com/openbiota-org/openbiota-gut-health-test/tree/main/sample-reports):
one community the index places as broadly disturbed, one in the healthy range.

## Before you start: open a terminal

Everything below is typed into a **terminal** — a window where you type a
command, press <kbd>Return</kbd>, and the computer runs it and prints what happened.
It looks like this:

![A terminal window with the three install commands typed into it](images/terminal.svg)

=== "macOS"

    Press <kbd>&#8984; Cmd</kbd> + <kbd>Space</kbd>, type `Terminal`, press <kbd>Return</kbd>. (It is also in
    *Applications → Utilities → Terminal*.)

=== "Ubuntu / Linux"

    Press <kbd>Ctrl</kbd> + <kbd>Alt</kbd> + <kbd>T</kbd>, or search your applications for *Terminal*.

Three things to know, and then you know enough:

- **Commands are typed exactly as shown**, then <kbd>Return</kbd>. The `$` at the
  start of a line in a picture is the terminal's prompt — do not type it.
  Copying a line from this page and pasting it in is fine, and the little
  copy button at the right of each grey box does that for you.
- **`git clone <address>` downloads a copy of a project.** It creates a
  folder named after the project, inside whatever folder the terminal is
  currently in (your home folder, when it opens). **`cd <folder>`** moves the
  terminal into that folder, which is where every later command is typed.
  That folder is a normal folder: you can open it in the Finder or the file
  manager, and you will, to drop your files in.
- **A command that is still working has not frozen.** The installer shows a
  progress bar and what it is doing; a long download can take a while. If
  you close the window or lose power, run the same command again and it
  continues from where it stopped.

## What you need

| | |
|---|---|
| A computer | A Mac (Intel or Apple silicon) or a Linux machine. 8 GB of memory is enough; more is faster. |
| Free disk | About **20 GB** for the standard install, **on the computer's own disk** — your home folder is the right place. Not a network share and not an external drive formatted for Windows (exFAT/NTFS): `git clone` itself fails there with *unable to write file … Permission denied*, and the installer refuses to start. (The optional full install is about 1.5 TB; see the end of this page.) |
| Internet | The installer downloads the reference databases once, roughly 3 GB. |
| About an hour | Mostly waiting for downloads and one index build. You can leave it. |
| Your sample | **Two FASTQ files** from a shotgun-sequenced stool sample, usually named something like `SAMPLE_1.fastq.gz` and `SAMPLE_2.fastq.gz`. |

The FASTQ files are the **raw reads**, not a PDF and not a 16S amplicon run.
If you tested with [Tiny Health](https://www.tinyhealth.com), their support
team will send download links for your raw FASTQ files on request. Any
Illumina paired-end stool metagenome of roughly five million read pairs or
more works; the report states its own depth and flags what depth limits.

## Step 1 — Install

Open a terminal and type (or paste) these three lines, pressing <kbd>Return</kbd>
after each:

```bash
git clone https://github.com/openbiota-org/openbiota-gut-health-test
cd openbiota-gut-health-test
./install.sh
```

What each line does:

1. **`git clone …`** downloads the project into a new folder called
   `openbiota-gut-health-test`, inside the folder the terminal is in — your
   home folder when it opens, which is where it should be. On a Mac that has
   never run a developer tool, macOS may first offer to install its *command
   line developer tools* — accept, wait for it to finish, then run the line
   again. If the clone ends with `unable to write file … Permission denied`,
   the terminal was in a network share or an external drive; type `cd ~`
   and run the line again.
2. **`cd openbiota-gut-health-test`** moves into that folder.
3. **`./install.sh`** does everything else, in eleven numbered steps with a
   progress bar:

   | | |
   |---|---|
   | Checks the machine | operating system, disk space, memory, internet |
   | Installs the tools | `diamond`, `bowtie2`, `samtools`, `fastp`, a recent Python. On a Mac it uses [Homebrew](https://brew.sh) and offers to install it if you do not have it; on Linux it uses `apt` or `dnf` and asks before using `sudo`. **It will ask for your computer password once**, for that. |
   | Builds the environments | the program itself and MetaPhlAn 3, each in its own folder inside the project; nothing else on your computer is changed |
   | Unpacks the reference data | 250 MB that came with the download: the ranges every reading is scored against, the reference cohort, the age model, and the protein references the search database is built from |
   | Downloads and builds the databases | the MetaPhlAn marker database (a 0.4 GB download, indexed to 3.5 GB), the protein search database (built locally), the metabolic model library |
   | Builds the human-read filter | downloads the human genome and indexes it; this is the long step, an hour or so on a laptop. `./install.sh --skip-host-index` leaves it out; the report then cannot state the human DNA fraction, which for stool is usually under 1 %. |
   | Checks everything | and prints the two steps left |

   Every step is skipped if it was already done, so if anything interrupts
   the install — a dropped connection, a closed laptop — just run
   `./install.sh` again and it continues. `./install.sh --list` shows what
   is done and what is not. The full detail of every step is in
   `install.log`, in the project folder.

!!! tip "It ended with a red ✗"
    Read the line after the ✗ — it says what went wrong and what to do — then
    run `./install.sh` again. The last lines of `install.log` are printed
    with it. [Installation](INSTALL.md#when-something-goes-wrong) lists the
    situations the installer knows how to explain.

## Step 2 — Add your sample

Inside the project folder there is a folder called **`fastq`**. Make a folder
in it for your sample, and copy your two FASTQ files into that folder:

```
openbiota-gut-health-test/
└── fastq/
    └── my-sample/
        ├── my-sample_1.fastq.gz
        └── my-sample_2.fastq.gz
```

That is, the two files end up at `fastq/my-sample/my-sample_1.fastq.gz` and
`fastq/my-sample/my-sample_2.fastq.gz`. You can do this in the Finder or the
file manager — drag the two files into the new folder — or in the terminal:

```bash
mkdir fastq/my-sample
cp ~/Downloads/my-sample_1.fastq.gz ~/Downloads/my-sample_2.fastq.gz fastq/my-sample/
```

A few rules, all of them simple:

- **Keep the file names you were given.** The sample is named by its files:
  `NAME_1.fastq.gz` and `NAME_2.fastq.gz` — or `NAME_R1` and `NAME_R2` —
  are the sample `NAME`. The folder's name is up to you.
- **Do not unzip, trim, filter or shorten the files.** The pipeline does its
  own quality assessment and host filtering and needs every read to
  normalise correctly. `.fastq.gz` and plain `.fastq` are both fine.
- **One folder per sample.** Several samples can sit side by side, each in
  its own folder; the run command below then asks which one you mean, and
  `bin/openbiota run-all` does them all and adds a comparison.
- Nothing you put in `fastq/` is ever added to git; the whole folder is
  ignored. There is a `README.md` in it that says the same thing.

## Step 3 — Run

In the terminal, still in the project folder:

```bash
bin/openbiota run
```

That is the actual program; `bin/openbiota --help` lists everything it can
do. (The installer also puts `openbiota` on your path, so in any *new*
terminal window the same command is just `openbiota run`.) With more than one
sample in `fastq/`, say which: `bin/openbiota run --sample my-sample`.

The run prints each stage as it goes. Expect **about an hour** for a sample
on a desktop, longer on a laptop: the read search and profiling take ten to
twenty minutes, and the metabolic-model simulation near the end takes
thirty to forty minutes on its own, because it solves several thousand
linear programmes for your particular community. Every stage caches its
work, so running the same sample again takes about five minutes. The report
is written to:

```
results/my-sample/my-sample_report.pdf
```

!!! note "The last lines say INCOMPLETE RUN"
    On the standard install they will. The run ends with a census of every
    stage, and the detection lanes that only the full install provides —
    the pathogen screen, the mycobiome, strain typing, the expanded
    catalogues — are listed as *not installed* or *not assessed*, in red,
    and the command exits non-zero so a script cannot mistake a partial
    report for a complete one. The report itself is written and says the
    same thing on its pages: each of those lanes is *not assessed*, never a
    negative result. `./install.sh --full` is what turns them on.

Open it:

=== "macOS"

    ```bash
    open results/my-sample/my-sample_report.pdf
    ```

=== "Linux"

    ```bash
    xdg-open results/my-sample/my-sample_report.pdf
    ```

Start at page 1. [Reading the report](OUTPUT.md) explains the seven-band
scale every reading is placed on and what each section means. Next to the
PDF, `results.json` holds every number in it, machine-readable, and
`run.log` the full log of the run.

## Optional: the full install

The standard install produces the complete report: community, pathways,
disease patterns, microbiome age, the organisms that need attention and what
the research says about them. Six further detection lanes — the 485-target
pathogen screen, the mycobiome, strain placement, whole-genome and expanded
catalogue detection — need about **1.5 TB** of reference data and most of a
day to fetch and build. When you have the disk:

```bash
./install.sh --full
```

It is resumable like everything else, and each lane is reported in the PDF
as *not assessed* until its references are present — never as a negative
result. [Installation](INSTALL.md) describes every step, and
[Organism detection](EXPANDED_DETECTION.md) the lanes.

## Useful next commands

```bash
bin/openbiota run --help             # every flag: threads, panels, subject data
bin/openbiota run --subject-age 42   # the age model also reports its error
bin/openbiota run --compact          # one-screen answer in the terminal
bin/openbiota run-all                # every sample in fastq/, plus a comparison
bin/openbiota doctor                 # is everything installed? what inputs are there?
bin/openbiota depth-check            # was the sequencing deep enough?
bin/openbiota panels                 # the metabolite pathways and their genes
bin/openbiota profiles               # the disease-pattern profiles
make test                            # the unit tests; no network or FASTQ needed
```

The full command reference is in [Usage](USAGE.md).
