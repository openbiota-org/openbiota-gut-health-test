"""Measurement engines.

Two engines answer two different questions about the same reads:

* **diamond** — *how much of gene X is present?* Translated search of every
  read against curated protein references, normalised to copies per 100
  bacterial genomes. Lives in `openbiota.search` and `openbiota.tally`; the module here
  is a thin facade so both engines present the same shape.
* **metaphlan** — *how much of species Y is present?* Marker-gene profiling
  with MetaPhlAn against a pinned database, after host-read removal.

Panels declare `engine: diamond`; profile features declare their engine per
feature. Both run in one invocation and merge into one report.
"""

from __future__ import annotations

from typing import Final

ENGINE_DIAMOND: Final = "diamond"
ENGINE_METAPHLAN: Final = "metaphlan"
ENGINES: Final = (ENGINE_DIAMOND, ENGINE_METAPHLAN)
