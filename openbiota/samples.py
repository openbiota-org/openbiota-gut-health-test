"""Sample identifiers that can be published, and the local names they stand for.

Everything the repository tracks - documentation, tests, frozen baselines -
refers to the project's specimens by neutral identifiers (``SAMPLE1_A01``,
short form ``SAMPLE1``). The directories under ``results/`` and the FASTQ
files may carry whatever names the sequencing provider gave them; the
mapping between the two lives in one local file that git ignores:

    results/.sample_aliases.json
    {"SAMPLE1_A01": "B2_XXXXXX", "SAMPLE2_A02": "..."}

Without the file, or for a name not in it, an identifier is its own local
name - a fresh checkout with samples named ``SAMPLE1_A01`` needs nothing.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
ALIAS_FILE = RESULTS / ".sample_aliases.json"


@lru_cache(maxsize=1)
def aliases() -> dict[str, str]:
    """Published identifier -> local sample name, from the local alias file."""
    if not ALIAS_FILE.is_file():
        return {}
    try:
        data = json.loads(ALIAS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {str(k): str(v) for k, v in data.items() if isinstance(v, str)}


def local_name(sample: str) -> str:
    """The local sample name for a published identifier (identity when unmapped)."""
    return aliases().get(sample, sample)


def published_name(local: str) -> str:
    """The published identifier for a local sample name (identity when unmapped)."""
    for alias, name in aliases().items():
        if name == local:
            return alias
    return local


def short_name(sample: str) -> str:
    """``SAMPLE1_A01`` -> ``SAMPLE1``: the part before the first underscore, as the FMT tool labels people."""
    return published_name(sample).split("_")[0]


def results_dir(sample: str) -> Path:
    """The results directory for a published identifier or a local name."""
    return RESULTS / local_name(sample)


def results_file(sample: str) -> Path:
    return results_dir(sample) / "results.json"


def reference_ranges() -> Path:
    """The cohort's reference ranges (`make cohort`), used by the pipeline and by the tests that read it."""
    return REPO / "refs" / "reference_ranges.json"


def report_pdf(sample: str) -> Path:
    local = local_name(sample)
    return RESULTS / local / f"{local}_report.pdf"


__all__ = ["ALIAS_FILE", "RESULTS", "aliases", "local_name", "published_name", "reference_ranges", "report_pdf",
           "results_dir", "results_file", "short_name"]
