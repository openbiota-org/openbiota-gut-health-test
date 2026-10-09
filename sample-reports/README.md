# Sample reports

Two complete OpenBiota Gut Health Test reports from real stool metagenomes,
produced by the software in this repository. They are here so a reader can see
exactly what a report contains before sequencing anything, and so the
documentation's screenshots come from real output.

| Report | What the index found | Reads |
|---|---|---|
| [`EM1_AMD614_report.pdf`](EM1_AMD614_report.pdf) | Broadly disturbed: an opportunist overgrown far above its typical carrier, several core species missing, butyrate capacity low | 8.0 million read pairs |
| [`MM1_FXX745_report.pdf`](MM1_FXX745_report.pdf) | Healthy range: health markers well ahead of disease markers, a protective gut lining, two organisms modestly above their typical carrier | 5.7 million read pairs |

Both are research-use reports, not diagnostic tests; see the note on page 1 of
either.

The reports and the pictures under `docs/images/em1` and `docs/images/mm1`
are regenerated together with `make sample-reports`
(`scripts/generate_sample_reports.py`) whenever the report changes, so what is
published here always shows the current software.
