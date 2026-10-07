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
| Input handling | Pair discovery, gzip preference, subsampling, truncated and malformed FASTQ |
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
openbiota run --subsample 100000     # end to end on real reads, ~17 s
openbiota validate                   # synthetic communities with known truth
```

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
