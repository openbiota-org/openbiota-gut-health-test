"""The build version has one source: ``specs/build_version.txt``."""

from __future__ import annotations

import re
from pathlib import Path

import openbiota

ROOT = Path(__file__).resolve().parents[1]


def test_version_is_read_from_the_build_version_file() -> None:
    raw = (ROOT / "specs" / "build_version.txt").read_text().strip().lstrip("vV")
    assert openbiota.__version__ == raw
    assert re.fullmatch(r"\d+\.\d+\.\d+([-+.][0-9A-Za-z.-]+)?", raw), raw


def test_pyproject_reads_the_same_file() -> None:
    text = (ROOT / "pyproject.toml").read_text()
    assert 'dynamic = ["version"]' in text
    assert 'version = { file = "specs/build_version.txt" }' in text
    assert not re.search(r'^version = "\d', text, re.M), "a static version in pyproject would drift from the file"
