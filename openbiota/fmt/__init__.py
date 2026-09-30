"""openbiota.fmt — recipient-specific FMT donor and donor-set matching.

Implements BUILD_SPEC_v06.0: compare one recipient against any number of
candidate donor materials, evaluate donor sets for complementary restoration
potential, and report exactly what each donor could supply, what remains
unresolved and what a combination adds.

What this package computes is a **research comparison**. Nothing here is a
clinical release decision, a safety clearance, a preparation instruction or a
statement that any donor will help anyone. The spec's separations are carried
through the code: measured coverage, predicted engraftment, simulated
metabolism and clinical outcome are four different things and never share a
number.
"""

from __future__ import annotations

#: Schema family version for request/result documents (spec §5.1, §13.3).
FMT_SCHEMA_VERSION = "6.0.0"

__all__ = ["FMT_SCHEMA_VERSION"]
