"""Strain, subtype and genomic-function resolution (BUILD_SPEC_v07.0).

Species profiling answers "what is here". It cannot answer "which one", and
for a great deal of gut biology that is the question that matters: whether an
*E. coli* is a commensal or carries an intact Shiga-toxin operon, whether a
*C. difficile* carries the toxin locus at all, whether a donor's population of
a species is the same population the recipient already has.

This package adds the layer that answers "which one", and — just as
importantly — says precisely when it cannot. Its rules come from the spec and
are all versions of one idea: **preserve the resolution the sequence supports,
never invent the resolution it does not.**

* A species-only abundance can never produce a strain fingerprint.
* An unrun, failed or unavailable assay is never a negative result.
* A nearest reference is not the identical cultured isolate.
* Two genes in one stool sample are not automatically in one organism.
* Nulls never become zeros, and a null transmission probability never renders
  as 0% or as a green all-clear.

Modules:

* `schema` — the evidence record and its three independent status axes.
* `registry` — target definitions, the compiler and the completeness census.
"""

from __future__ import annotations

RESOLUTION_SCHEMA_VERSION = "openbiota.resolution/1.0"

__all__ = ["RESOLUTION_SCHEMA_VERSION"]
