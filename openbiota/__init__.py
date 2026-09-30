"""openbiota — the OpenBiota Gut Health Test pipeline.

Turns raw shotgun stool metagenome reads (FASTQ) into a readable gut-health
report: metabolite pathway gene capacity normalised to bacterial genome
equivalents via ``rpoB``, species-level community composition, a published
gut-health index, and resemblance to published disease-associated patterns.

This package measures *genetic capacity* and *community resemblance*. It does
not measure metabolites and it does not diagnose disease. See ``README.md``.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["BUILD_VERSION_FILE", "__version__"]

#: The build version lives in one place, ``specs/build_version.txt``, and is
#: read at import. The report header and footer, ``results.json`` and
#: ``openbiota --version`` all print this value, so bumping the file bumps
#: them together. A missing or empty file yields a visibly unknown version
#: rather than a stale hard-coded one.
BUILD_VERSION_FILE = Path(__file__).resolve().parent.parent / "specs" / "build_version.txt"


def _build_version() -> str:
    try:
        raw = BUILD_VERSION_FILE.read_text(encoding="utf-8").strip()
    except OSError:
        return "0.0.0+unknown"
    return raw.lstrip("vV").strip() or "0.0.0+unknown"


__version__ = _build_version()
