# Contributing

OpenBiota turns raw shotgun stool metagenome reads into a report a person can
read. Every contribution that brings more published research into that report,
makes a finding easier to understand, or makes a number more accurate is
welcome.

## Before you start

- Read [docs/METHOD.md](docs/METHOD.md) for how the measurements are made and
  why, and [docs/TESTING.md](docs/TESTING.md) for how the suite is run.
- Open an issue for anything larger than a typo so the design can be agreed
  before the work is done. Bug reports need the `run.log` of the run and the
  `openbiota --version` output; never attach `results.json`, a report PDF or
  reads from a real person.
- The licence is the [PolyForm Noncommercial License 1.0.0](LICENSE). By
  submitting a contribution you agree that it is licensed under the same terms
  and that OpenBiota may also offer it under its commercial licence
  (a Developer Certificate of Origin sign-off - `git commit -s` - records
  this on each commit).

## What a change needs

1. **A test that fails before and passes after.** The suite is the
   specification: `make test` (about 2,500 tests, no network). A change to a
   number that a report prints needs the test that pins the number.
2. **Provenance for any new reference or claim.** A panel, profile, group or
   intervention entry carries its source (DOI or accession), the release it
   came from and the date it was retrieved, in the same YAML/JSON the existing
   entries use. A tool or reference release is pinned by version and checksum
   under `refs/expanded/`.
3. **One number per quantity.** If a quantity is printed in two places it is
   computed once and both places read it. The tests in
   `tests/test_one_number_per_quantity.py` and
   `tests/test_inventory_counts_reconcile.py` enforce this for the report.
4. **Nothing personal.** No sample from a real person, no report, no sequencing
   read, no name, in code, tests, fixtures or documentation. Test fixtures use
   the published sample identifiers (`SAMPLE1_A01`, …; see
   `openbiota/samples.py`); frozen baselines are anonymised by
   `scripts/freeze_preservation_baselines.py` before they are written.
5. **Plain language on the page.** The report is for the person the sample came
   from. A sentence a reader has to look up is a sentence to rewrite; a claim
   the literature does not support is a claim to remove.

## Style

- Python 3.11+, type hints, `ruff` clean (`make lint`). Docstrings say what
  the function is for and why it is the way it is, not what each line does.
- Report text: British spelling, no marketing, no exclamation marks, numbers
  with units, every reference a working link.
- Commit messages describe the change and its reason in one sentence; the
  first line is the summary.

## Running the pieces

```bash
make setup                 # virtualenv and the package
make test                  # the suite
make lint                  # ruff
make docs                  # documentation preview at http://127.0.0.1:8000
scripts/build_site.py      # rebuild web/openbiota.com (docs, sitemap, icons)
```

A full report needs the reference data and the profilers described in
[docs/INSTALL.md](docs/INSTALL.md); the suite does not.

## Adding a pathway, profile, group or intervention

Each has its own guide: [docs/PANELS.md](docs/PANELS.md) (metabolite
pathways), [docs/PROFILES.md](docs/PROFILES.md) (disease-pattern profiles),
`openbiota/taxongroups.py`, `taxa/` and `interventions/` (microbial groups, organism descriptions and the
intervention evidence). Every entry needs its citation, its evidence grade
and a sentence a reader can act on.

## Reporting a security or privacy problem

See [SECURITY.md](SECURITY.md). Do not open a public issue for a problem that
could expose a person's data.
