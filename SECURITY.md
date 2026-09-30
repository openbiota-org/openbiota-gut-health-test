# Security and privacy

OpenBiota runs entirely on the machine that holds the reads. It makes network
requests only to fetch reference databases and tools during installation
(`make setup`, `make build-db`, `make refs-expanded`), and to the BacDive API
during a run for culture-collection metadata about named species; no sample
data leaves the machine. The test suite makes no network requests.

## What is sensitive

A stool metagenome identifies a person as surely as a fingerprint, and the
report describes their health. Treat these as personal data at every step:

- the FASTQ files (`fastq/`), which git ignores;
- everything under `results/` (`results.json`, the report, logs, alignments),
  which git ignores;
- third-party reports (`tinyhealth/`), which git ignores.

The repository tracks no reads, no results and no report. Frozen test
baselines under `tests/fixtures/preservation/` are derived from real runs but
carry only published sample identifiers (`SAMPLE1_A01`, …) and no paths or
names; `scripts/freeze_preservation_baselines.py` anonymises them before they
are written.

## Reporting a problem

If you find a way the software could expose a person's data, a dependency
with a known vulnerability, or an unsafe default, email
**security@openbiota.com** with the details and a way to reproduce it. Please
do not open a public issue until a fix is available. You will get an
acknowledgement within three days and a fix or a plan within fourteen.

## Supported versions

The latest release on the `main` branch. Reference and tool versions are
pinned by checksum under `refs/expanded/`; a run records every version it
used in `results.json` and `run.log`.
