---
name: Research addition
about: A published study, pathway, profile, group or intervention the report should use
labels: research
---

**Source**
DOI or accession, and the exact claim it supports.

**Where it belongs**
Metabolite pathway (`panels/`), disease-pattern profile (`profiles/`),
microbial group (`openbiota/taxongroups.py`), intervention evidence
(`interventions/`), organism description (`taxa/interpretations/`, `taxa/classification/`), or somewhere else.

**Evidence grade**
Study design, population, effect size, replication. See the grading in
`docs/PROFILES.md`, and the `lane` (A–E, X) and `study.risk_of_bias` fields of the entries in `interventions/assertions.yaml`.

**What the report would say**
One sentence a reader can act on, and what it must not claim.
