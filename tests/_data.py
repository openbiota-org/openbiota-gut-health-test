"""Helpers for tests that read a finished run or the built reference data.

A fresh checkout has neither: `results/` is written by a run and
`refs/reference_ranges.json` by `make cohort`. Tests that need them skip
rather than fail, so the suite is green on a clean machine and exercises
everything it can there.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.samples import reference_ranges, results_file


def results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


def results_path(sample: str) -> Path:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return path


def ranges() -> dict[str, Any]:
    path = reference_ranges()
    if not path.is_file():
        pytest.skip("refs/reference_ranges.json is not built in this checkout (make cohort)")
    return json.loads(path.read_text(encoding="utf-8"))
