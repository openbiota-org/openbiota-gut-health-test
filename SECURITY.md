# Security and privacy

OpenBiota runs entirely on the machine that holds the reads. It makes network
requests only to fetch reference databases and tools during installation
(`./install.sh`, and the `make` targets it runs: `build-db`, `extension-refs`,
`host-index`, `refs-expanded` and so on), and to the BacDive API during a run
for culture-collection metadata about named species; no sample data leaves
the machine. The test suite makes no network requests. Nothing under `fastq/`
or `results/` is tracked by git, at any depth.

## What is sensitive

A stool metagenome identifies a person as surely as a fingerprint, and the
report describes their health. Treat these as personal data at every step:

- the FASTQ files (`fastq/`), which git ignores;
- everything under `results/` (`results.json`, the report, logs, alignments),
  which git ignores;

The repository tracks no reads, no results, no report and no fixture derived
from a real person's sample. Tests that need a finished run read it from the
local `results/` directory and skip when it is absent, so the suite runs
anywhere; every assertion about a report is a property any sample must
satisfy, never a frozen copy of one person's numbers.

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
