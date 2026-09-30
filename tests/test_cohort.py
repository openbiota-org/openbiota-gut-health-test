"""Reference ranges, percentile placement and status interpretation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.cohort import PERCENTILES, PanelRange, ReferenceRanges, build_ranges
from openbiota.panels import load_panel_set
from openbiota.pdfreport import AMBER, CORAL, GREEN, SLATE, status_for

RANGES = Path(__file__).resolve().parent.parent / "refs" / "reference_ranges.json"


def _range(values: list[float], name: str = "test") -> PanelRange:
    return build_ranges(
        study="TEST",
        description="test",
        reads_per_sample=1000,
        panel_values={name: values},
        gene_values={},
    ).panels[name]


# --------------------------------------------------------------------------- #
# percentile construction
# --------------------------------------------------------------------------- #


def test_percentiles_of_a_uniform_distribution():
    pr = _range([float(v) for v in range(1, 101)])

    assert pr.n == 100
    assert pr.percentiles[50] == pytest.approx(50.5, abs=0.6)
    assert pr.percentiles[25] == pytest.approx(25.75, abs=0.6)
    assert pr.percentiles[75] == pytest.approx(75.25, abs=0.6)
    assert pr.min_value == 1.0
    assert pr.max_value == 100.0
    assert set(pr.percentiles) == set(PERCENTILES)


def test_percentiles_are_monotonic():
    pr = _range([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0])
    values = [pr.percentiles[p] for p in sorted(PERCENTILES)]
    assert values == sorted(values)


def test_detection_rate():
    pr = _range([0.0, 0.0, 1.0, 2.0, 3.0])
    assert pr.n_detected == 3
    assert pr.detection_rate == pytest.approx(0.6)


def test_single_value_distribution_does_not_crash():
    pr = _range([7.0])
    assert pr.median == 7.0
    assert pr.stdev == 0.0


# --------------------------------------------------------------------------- #
# placing a value
# --------------------------------------------------------------------------- #


def test_percentile_of_places_values_sensibly():
    pr = _range([float(v) for v in range(0, 101)])

    assert pr.percentile_of(pr.percentiles[50]) == pytest.approx(50, abs=2)
    assert pr.percentile_of(pr.percentiles[25]) == pytest.approx(25, abs=2)
    assert pr.percentile_of(pr.percentiles[95]) == pytest.approx(95, abs=2)


def test_values_outside_the_range_are_clamped():
    pr = _range([10.0, 20.0, 30.0, 40.0, 50.0])

    assert pr.percentile_of(-100.0) == float(min(PERCENTILES))
    assert pr.percentile_of(1e9) == float(max(PERCENTILES))


def test_classify_wording():
    pr = _range([float(v) for v in range(0, 101)])

    assert pr.classify(None) == "not measurable"
    assert "below" in pr.classify(pr.percentiles[5] - 1)
    assert pr.classify(pr.percentiles[50]) == "within the usual range"
    assert "above" in pr.classify(pr.percentiles[95] + 1)


# --------------------------------------------------------------------------- #
# round trip
# --------------------------------------------------------------------------- #


def test_reference_ranges_round_trip(tmp_path: Path):
    original = build_ranges(
        study="TEST",
        description="a test cohort",
        reads_per_sample=600_000,
        panel_values={"urda": [10.0, 20.0, 30.0], "cutc": [1.0, 2.0, 3.0]},
        gene_values={"urda:URDA": [10.0, 20.0, 30.0]},
    )
    path = tmp_path / "ranges.json"
    original.save(path)
    restored = ReferenceRanges.load(path)

    assert restored.study == original.study
    assert restored.n_samples == original.n_samples
    assert restored.panels["urda"].percentiles == original.panels["urda"].percentiles
    assert restored.genes["urda:URDA"].median == original.genes["urda:URDA"].median


# --------------------------------------------------------------------------- #
# status = position read in the direction the literature reports
# --------------------------------------------------------------------------- #


def test_adverse_direction_flags_high_values():
    high = status_for(percentile=95, higher_means="adverse", detected=True)
    mid = status_for(percentile=50, higher_means="adverse", detected=True)
    low = status_for(percentile=5, higher_means="adverse", detected=True)

    assert high.colour is CORAL
    assert mid.colour is GREEN
    assert low.colour is GREEN          # low is good when higher is adverse


def test_favourable_direction_flags_low_values():
    high = status_for(percentile=95, higher_means="favourable", detected=True)
    mid = status_for(percentile=50, higher_means="favourable", detected=True)
    low = status_for(percentile=5, higher_means="favourable", detected=True)

    assert high.colour is GREEN
    assert mid.colour is GREEN
    assert low.colour is CORAL          # low is bad when higher is favourable


def test_intermediate_values_are_amber():
    assert status_for(percentile=80, higher_means="adverse", detected=True).colour is AMBER
    assert status_for(percentile=20, higher_means="favourable", detected=True).colour is AMBER


def test_context_dependent_gives_position_without_a_verdict():
    for percentile in (5, 50, 95):
        status = status_for(
            percentile=percentile, higher_means="context-dependent", detected=True
        )
        assert status.colour is SLATE


def test_provisional_never_reads_as_a_finding():
    """A result below the confirmed threshold must not be presented as one.

    The word still places the reading on the shared scale — the report never
    switches vocabulary — but the colour withholds judgement, which is what
    keeps it out of the headline findings.
    """
    status = status_for(
        percentile=99, higher_means="adverse", detected=True, confidence="provisional"
    )
    assert status.label == "notably high"
    assert status.colour is SLATE
    assert "threshold" in status.note


def test_undetected_and_unmeasurable():
    assert status_for(percentile=None, higher_means="adverse", detected=True).label == "no result"
    assert status_for(percentile=50, higher_means="adverse", detected=False).label == "not detected"


# --------------------------------------------------------------------------- #
# the shipped ranges
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not RANGES.is_file(), reason="run `openbiota cohort` first")
def test_shipped_ranges_cover_every_panel():

    ranges = ReferenceRanges.load(RANGES)
    panels = load_panel_set(Path(__file__).resolve().parent.parent / "panels")

    assert ranges.n_samples >= 20, "a cohort this small makes percentiles unreliable"
    for panel in panels.panels:
        assert panel.name in ranges.panels, f"no reference range for {panel.name}"


@pytest.mark.skipif(not RANGES.is_file(), reason="run `openbiota cohort` first")
def test_shipped_ranges_are_plausible():
    ranges = ReferenceRanges.load(RANGES)
    panels = load_panel_set(Path(__file__).resolve().parent.parent / "panels")

    for name, pr in ranges.panels.items():
        assert pr.percentiles[5] <= pr.median <= pr.percentiles[95], name
        assert pr.median >= 0, name
        # A single-copy gene above ~150 copies per 100 genomes is implausible;
        # the same ceiling the run-time diagnostic uses, scaled for `sum`
        # panels that add several near-universal genes (bcaa, riboflavin) and
        # for declared multi-copy families (GUS).
        limit = panels.by_name(name).plausible_ceiling(200.0)
        assert pr.percentiles[95] < limit, f"{name} 95th percentile is implausible"

    # butyrate producers are abundant in every healthy adult stool sample
    assert ranges.panels["butyrate"].detection_rate == 1.0
    assert ranges.panels["butyrate"].median > 20


@pytest.mark.skipif(not RANGES.is_file(), reason="run `openbiota cohort` first")
def test_cohort_provenance_is_recorded():
    payload = json.loads(RANGES.read_text(encoding="utf-8"))
    assert payload["study"]
    assert payload["description"]
    assert payload["reads_per_sample"] > 0
    assert payload["built_at"]
