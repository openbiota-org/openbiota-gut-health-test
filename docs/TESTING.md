# Testing and development

---

## Running the tests

```bash
make test            # or: .venv/bin/python -m pytest
make lint            # ruff
```

No network access and no FASTQ required. Fixtures are synthetic DIAMOND TSV
lines plus an in-memory reference database with hand-chosen sequence lengths,
so the normalisation arithmetic can be checked against values computed by hand.

---

## What is covered

| Area | Tests |
|---|---|
| Mate collapse | 40 fragments seen as 80 read hits count as 40; `/1` `/2` suffix stripping; higher-scoring mate wins |
| Decoy rejection | A decoy best-hit counts as rejected, never as a target; weak decoy hits are neither |
| Thresholds | Identity floor and minimum alignment length, applied per gene |
| Normalisation | Against hand-computed values; depth independence; gene-length comparability |
| Zero `rpoB` | Returns indeterminate rather than dividing by zero |
| Aggregation | `sum` / `median` / `max` / `mean`, and `aggregate_from` narrowing |
| Residue mapping | Position mapping through gapped alignments, gaps in either sequence, reversed coordinates |
| Residue checks | Reference-side split, read-side classification, concordance and discordance detection |
| Organism attribution | Sums to the fragment total, not the read total — asserted at runtime too |
| Pathway validation | Every malformed-YAML case fails with a specific message |
| Query hazards | `gene:cutC` never used; no comma-inside-parentheses in any shipped query; excluded genes stay excluded |
| Input handling | Pair discovery in `fastq/` and one folder down, the file naming the sample, hidden files ignored, the same sample in two folders refused, gzip preference, subsampling, truncated and malformed FASTQ |
| Installation | `install.sh` under shellcheck and bash 3.2 rules; the first step says when the machine's memory is short of what the MetaPhlAn 4 lane holds; every advertised step has its functions, label and help text; options parse and bad ones are refused; `--list`, `--dry-run`, the lock, the log; the `launcher` step on a throwaway `$HOME` (link, `PATH` line, idempotence, a foreign file left alone); `bin/openbiota` with and without a `.venv` and through a symlink; every readiness check on an empty folder; `make install` |
| Reference-data archive | `data/derived-refs.tar.gz` matches its manifest, fits GitHub's file limit, holds every promised item and no reads or results; pack → restore round trip, existing files kept unless `--force`, tampering and path traversal refused, reproducible bytes, split-and-join when the part size is exceeded |
| The run command | `openbiota` alone prints help and starts nothing; flags without a command are refused; positional inputs resolve (a pair in either order, one file finding its mate, a folder, three files refused); inputs that are not FASTQ (FASTA, PDF, HTML, empty, a gzip cut short, two different samples as a pair) are refused by name; several samples without a terminal list the ways to choose; full/lite in every non-interactive case, full refused without the full install; the census asks a lite run only about lite stages; the report path is the last line; the viewer opens once, at a terminal only, never under `CI` or `OPENBIOTA_NO_OPEN`; the two prompts driven on a real pty |
| What a run prints | The short form is one line per stage with no timestamps, and a stage too quick to matter leaves none; on a terminal the running stage is one line redrawn in place that never wraps, a step inside it shows beside it, and work between stages shows as *working*; into a log a long stage is announced once; a failed stage says why once, in words, and the handler's explanation is not printed a second time; a warning is one line, said once; `--quiet` prints nothing; the log keeps every line as plain text; the ending is one line when complete and names only the missing stages, each with a reason a person can act on, when not; the text summary is printed only with `--summary`; `Report: <path>` is the only thing on stdout and the last line |
| Self-test | The reads are a gzipped pair with a distinct name per read; a lane passes only if it ran *and* found the genera the reads were drawn from, however the catalogue spells them; a failed, missing or crashing lane is reported with its reason, never raised; each lane is measured in a process of its own and a measuring process that dies is a failed lane, not a hang; exit status 0 / 1 / 2; the work folder is removed unless kept; the installer's final check runs it and a failure stops the install |
| Lane robustness | Kraken 2 is told the reads are gzipped only when they are, and reading nothing is a failure, not an empty gut; sylph-tax is given the taxonomy folder and no per-user config, and its failure carries what it said; SingleM names a missing program before it starts, and the readiness check and the installer both count HMMER and OrfM; an optional lane that is not installed is recorded, not announced as a finished stage; the project's own tool folders are on the `PATH` whichever way the program is started |
| Report modes | `SoftPageBreak` never yields an empty page and no report module uses a hard `PageBreak()`; lite numbering is contiguous and every module sees it; the lite foot is the solutions panel alone at the pair's height; the title carries the mode; the lite tile grid has the genus count where the strain tile was |
| Installation docs | Every `install.sh` option the docs mention exists; every step is described; the quickstart tells the three steps in order and explains the terminal first; the run step is the program, not a make target; the `fastq/` layout drawn in the docs is the one the code reads; the installer's closing message matches the quickstart; CI exercises the installer on Ubuntu and macOS |
| Cohort ranges | Percentile interpolation, classification, old/new phylum name merging |
| QC | Read statistics, caching and invalidation, binned-quality detection |
| Reporting | Band boundaries, diagnostic firing, JSON completeness |
| CLI | Argument defaults, command dispatch, clean error presentation |

Two invariants are asserted at runtime, not only in tests:

- organism attribution must reconcile with the fragment total;
- the residue anchor must map its own canonical position onto itself, or the
  database build aborts.

---

## Profile and scoring coverage

| File | Covers |
|---|---|
| `test_profiles.py` | Profile schema, the conflicted-genus refusal, weight-factor validation, duration strata, the shipped profiles |
| `test_scoring.py` | The compositional transform, robust z, module aggregation, percentiles, abstention, confidence |
| `test_dysbiosis.py` | GMWI2 against the authors' own worked example, and its refusal on an incompatible database |
| `test_profilevalidate.py` | ROC/AUC with tie correction, per-study holdout, specificity verdicts |

Two regression tests are worth knowing about, because both bugs produced
plausible-looking numbers:

- `test_clr_zero_floor_does_not_depend_on_the_samples_minimum` — the zero
  floor must be a constant, not derived per sample, or the transform becomes
  depth-dependent.
- `test_absent_sample_lands_mid_tie_block` and
  `test_absent_sample_contributes_almost_nothing` — an undetected taxon must
  not read as evidence. Before the fix, five absent colorectal markers scored
  at the 74th percentile on reference-frame noise alone.

## Integration checks

```bash
bin/openbiota run --subsample 100000 # end to end on real reads, ~17 s
bin/openbiota validate               # synthetic communities with known truth
```

The installer has an integration check of its own in CI: a fresh Ubuntu and
a fresh macOS runner clone the repository and run `./install.sh --only
preflight,system,python,launcher,profiler,data --yes --no-path`, then run it
again and require every step to report *already done*. The steps left out
are the ones that download gigabytes of reference data.

`openbiota validate` is the substantive one — it builds samples whose correct answer
is known and measures what the pipeline reports. See
[VALIDATION.md](VALIDATION.md).

---

## Layout

```
openbiota/
├── cli.py           argparse and orchestration only
├── panels.py        pathway YAML loading and validation
├── references.py    UniProt fetch, filtering, DIAMOND database
├── residues.py      anchor mapping and active-site checks
├── search.py        DIAMOND wrapper, caching, resume
├── tally.py         parsing, mate collapse, normalisation
├── interpret.py     cohort comparison and status logic
├── cohort.py        reference range construction
├── validate.py      validation experiments
├── simulate.py      read simulator
├── taxonomy.py      community profile
├── qc.py            input validation and read statistics
├── report.py        text and JSON output
├── pdfreport.py     consumer PDF
├── net.py           HTTP with retry
├── seqio.py         FASTA helpers
├── logging_util.py  stderr progress and run log
└── errors.py        exception hierarchy
```

The science lives in `tally.py`, `references.py` and `residues.py`. `cli.py`
orchestrates and contains no scientific logic.

---

## Conventions

- Python 3.11+, type hints throughout, no bare `except`.
- `pathlib` over string paths.
- `subprocess.run(check=True)` with the failing command and stderr in the
  error message.
- Fail loudly on a missing dependency or an empty reference fetch. A
  zero-sequence fetch is a hard error, never a count of zero.
- Every expensive stage is cached with a fingerprint over everything that
  affects its output.
- Long-running stages print progress to stderr.

Line length 100. `make lint` must pass clean.

---

## Adding a pathway

See [PANELS.md](PANELS.md). After adding one:

```bash
openbiota panels --show mypanel
openbiota build-db
openbiota run --subsample 100000
openbiota validate                   # confirm accuracy has not degraded
make test
```

If the pathway needs a reference genome that carries its gene for validation,
add it to `VALIDATION_GENOMES` in `openbiota/validate.py`.

---

## Debugging a suspicious result

1. `results/<sample>/run.log` — flags, timings, reference set sizes.
2. Diagnostics in `summary.txt` — most surprises are already flagged there.
3. `hits_<panel>_<mate>.tsv` — per-panel alignments with organism appended.
   The identity distribution and which organisms matched usually explain it.
4. `openbiota build-db` — kept counts, mean lengths and length windows per entry.
5. `openbiota panels --show <panel>` — the exact queries used.

A gene with a low median identity is matching distant relatives rather than
itself; a well-covered gene shows a median near 100%, because a read from a
sequenced carrier matches its own reference exactly.
