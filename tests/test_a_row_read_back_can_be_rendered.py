"""A row rebuilt from results.json must behave like the row that was written.

Round-tripping the *values* is not enough. The renderer asks a row's
drivers for `.top` and for `named_present`, and a record stored as a dict
and read back as a dict answers neither. That failure surfaced a hundred
pages into a build, after eight minutes of screening, with the json
written and no PDF beside it - so these tests go at the object, not at
the page.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.drivers import Driver, Drivers
from openbiota.pdfreport import MetaboliteRow, rows_from_json
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06")


def _results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has no results on disk")
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(sample: str) -> list[MetaboliteRow]:
    results = _results(sample)
    if not results.get("report_rows"):
        pytest.skip(f"{sample} predates report_rows in results.json")
    return rows_from_json(results)


# --- the records, in isolation --------------------------------------------


def test_a_driver_survives_the_round_trip() -> None:
    original = Driver(
        organism="Phocaeicola vulgatus", fragments=412, share=0.37,
        in_sample=True, sample_percent=6.7, sample_class="commensal",
        genus_in_sample=True,
    )
    assert Driver.from_json(original.to_json()) == original


def test_a_drivers_record_survives_the_round_trip() -> None:
    original = Drivers(
        panel="agmatine",
        entry_ids=("ADIA", "SPEB"),
        total_fragments=1100,
        top=(
            Driver("A", 600, 0.55, True, 3.1, "commensal", True),
            Driver("B", 500, 0.45, False, None, None, False),
        ),
        resolved_direction="favourable",
        resolved_reason="the carriers agree",
    )
    assert Drivers.from_json(original.to_json()) == original


def test_an_empty_drivers_record_round_trips() -> None:
    original = Drivers(panel="x", entry_ids=(), total_fragments=0, top=())
    assert Drivers.from_json(original.to_json()) == original


# --- the rows the report is actually built from ---------------------------


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_row_carries_its_drivers_as_a_bare_mapping(sample: str) -> None:
    """The exact shape that broke the build."""
    for row in _rows(sample):
        if row.drivers is None:
            continue
        assert isinstance(row.drivers, Drivers), (
            f"{row.panel}: drivers came back as {type(row.drivers).__name__}, "
            "which the renderer cannot ask for .top"
        )


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_row_answers_what_the_renderer_asks_it(sample: str) -> None:
    """Touch each attribute the detail page reaches for, on every row.

    Cheaper than rendering, and it fails on the row rather than on page
    one hundred.
    """
    for row in _rows(sample):
        assert isinstance(row.status.label, str)
        assert isinstance(row.detected, bool)
        drv = row.drivers
        if drv is None:
            continue
        assert isinstance(drv.top, tuple)
        assert isinstance(drv.named_present, tuple)
        for driver in drv.top:
            assert isinstance(driver, Driver), (
                f"{row.panel}: a driver came back as {type(driver).__name__}"
            )
            assert isinstance(driver.organism, str)
            assert isinstance(driver.in_sample, bool)


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_rebuilt_rows_match_the_file_they_came_from(sample: str) -> None:
    rebuilt = _rows(sample)  # skips a file written before the contract existed
    stored = _results(sample)["report_rows"]
    assert len(rebuilt) == len(stored)
    for row, payload in zip(rebuilt, stored, strict=True):
        assert row.panel == payload["panel"]
        assert row.value == payload.get("value")
        assert row.percentile == payload.get("percentile")
        if payload.get("drivers") is None:
            assert row.drivers is None
        else:
            assert row.drivers is not None
            assert len(row.drivers.top) == len(payload["drivers"].get("top") or [])
