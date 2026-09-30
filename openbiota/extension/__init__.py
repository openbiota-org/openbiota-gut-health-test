"""The v0.8.3 report extension: additional measurements over the same run.

Everything in this package is *additive*. The rule the whole package is built
around, from BUILD_SPEC_v0.8.3 section 1: an existing metric's value, status
label, component values and calculation fingerprint must be identical for
identical resolved inputs, with or without the extension. New measurements
live in their own objects with their own IDs, and a new measurement may never
enter a protected model's input graph.

Layout:

``schema``
    The extension metric contract - one stable ID, one unit, one denominator,
    and states that keep zero, null and not-assessed apart.
``preservation``
    The protected-output manifest and the comparison that proves nothing
    existing moved.
``registry``
    Which capabilities this release owes, what implements each one, and the
    coverage report that says so honestly.
``sources``
    The embedded research register: every claim resolves to a primary source
    with its actual access terms.
"""

from __future__ import annotations

#: The specification this package implements. Recorded in every emitted
#: object so a stored result can be traced to the contract it was built to.
SPEC_VERSION: str = "0.8.3"

__all__ = ["SPEC_VERSION"]
