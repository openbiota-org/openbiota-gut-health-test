"""Saturation analysis: convergence judgement and its arithmetic."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from openbiota.depth import (
    MIN_FRAGMENTS_FOR_VERDICT,
    DepthPoint,
    DepthReport,
    assess_convergence,
    poisson_relative_error,
)

DEPTH_JSON = (
    Path(__file__).resolve().parent.parent / "results" / "A02" / "depth_check.json"
)


def _points(
    values: list[tuple[int, float | None, int]], panel: str = "urda"
) -> list[DepthPoint]:
    return [
        DepthPoint(
            panel=panel,
            read_pairs=depth,
            fragments=fragments,
            rpob_fragments=depth // 1000,
            copies_per_100=value,
        )
        for depth, value, fragments in values
    ]


# --------------------------------------------------------------------------- #
# Poisson error
# --------------------------------------------------------------------------- #


def test_poisson_relative_error():
    assert poisson_relative_error(100) == pytest.approx(0.10)
    assert poisson_relative_error(10_000) == pytest.approx(0.01)
    assert poisson_relative_error(0) is None


def test_error_shrinks_with_count():
    assert poisson_relative_error(1000) < poisson_relative_error(100)


# --------------------------------------------------------------------------- #
# convergence judgement
# --------------------------------------------------------------------------- #


def test_a_flat_series_converges():
    points = _points([(1_000_000, 35.0, 170), (2_000_000, 35.1, 340)])
    result = assess_convergence(points, "urda")

    assert result.converged
    assert "within counting noise" in result.verdict
    assert result.relative_change is not None
    assert result.relative_change < result.poisson_error * math.sqrt(3)


def test_a_still_rising_series_does_not_converge():
    # doubling depth doubles the value — nowhere near a plateau
    points = _points([(1_000_000, 10.0, 1000), (2_000_000, 20.0, 2000)])
    result = assess_convergence(points, "urda")

    assert not result.converged
    assert "still changing" in result.verdict


def test_noise_threshold_scales_with_fragment_count():
    """The same relative change is noise at low counts and signal at high ones."""
    sparse = assess_convergence(
        _points([(1_000_000, 10.0, 40), (2_000_000, 11.0, 80)]), "urda"
    )
    dense = assess_convergence(
        _points([(1_000_000, 10.0, 40_000), (2_000_000, 11.0, 80_000)]), "urda"
    )

    assert sparse.converged, "10% change on 80 fragments is within counting noise"
    assert not dense.converged, "10% change on 80,000 fragments is real"


def test_too_few_fragments_gives_no_verdict():
    points = _points(
        [(1_000_000, 1.0, 5), (2_000_000, 2.0, MIN_FRAGMENTS_FOR_VERDICT - 1)]
    )
    result = assess_convergence(points, "urda")

    assert not result.converged
    assert result.verdict == "too few fragments to judge"
    assert result.relative_change is None


def test_indeterminate_values_give_no_verdict():
    points = _points([(1_000_000, None, 500), (2_000_000, None, 1000)])
    result = assess_convergence(points, "urda")

    assert not result.converged
    assert "not measurable" in result.verdict


def test_a_single_depth_cannot_be_judged():
    result = assess_convergence(_points([(1_000_000, 35.0, 500)]), "urda")
    assert not result.converged
    assert "not enough depths" in result.verdict


def test_zero_value_at_full_depth_is_not_measurable():
    points = _points([(1_000_000, 0.0, 100), (2_000_000, 0.0, 200)])
    result = assess_convergence(points, "urda")
    assert "not measurable" in result.verdict


# --------------------------------------------------------------------------- #
# report assembly
# --------------------------------------------------------------------------- #


def test_report_summary_when_everything_converged():
    points = _points([(1_000_000, 35.0, 500), (2_000_000, 35.1, 1000)])
    report = DepthReport(
        sample="T",
        full_read_pairs=2_000_000,
        fractions=[0.5, 1.0],
        points=points,
        convergence=[assess_convergence(points, "urda")],
    )

    assert report.all_converged
    assert "deep enough" in report.summary_sentence()


def test_report_summary_names_what_is_still_drifting():
    points = _points([(1_000_000, 10.0, 5000), (2_000_000, 20.0, 10_000)])
    report = DepthReport(
        sample="T",
        full_read_pairs=2_000_000,
        fractions=[0.5, 1.0],
        points=points,
        convergence=[assess_convergence(points, "urda")],
    )

    assert not report.all_converged
    assert "urda" in report.summary_sentence()


def test_report_ignores_unjudgeable_panels_in_the_verdict():
    good = _points([(1_000_000, 35.0, 500), (2_000_000, 35.1, 1000)], panel="urda")
    thin = _points([(1_000_000, 0.1, 2), (2_000_000, 0.2, 4)], panel="bai")
    report = DepthReport(
        sample="T",
        full_read_pairs=2_000_000,
        fractions=[0.5, 1.0],
        points=good + thin,
        convergence=[
            assess_convergence(good + thin, "urda"),
            assess_convergence(good + thin, "bai"),
        ],
    )

    # bai cannot be judged, so it must not veto the overall verdict
    assert report.all_converged


def test_report_round_trip(tmp_path: Path):
    points = _points([(1_000_000, 35.0, 500), (2_000_000, 35.1, 1000)])
    report = DepthReport(
        sample="T",
        full_read_pairs=2_000_000,
        fractions=[0.5, 1.0],
        points=points,
        convergence=[assess_convergence(points, "urda")],
        generated="now",
    )
    path = tmp_path / "depth.json"
    report.save(path)
    payload = DepthReport.load(path)

    assert payload["sample"] == "T"
    assert payload["all_converged"] is True
    assert len(payload["points"]) == 2
    assert "method" in payload


# --------------------------------------------------------------------------- #
# the real result for the shipped sample
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not DEPTH_JSON.is_file(), reason="run `openbiota depth-check` first")
def test_shipped_sample_was_sequenced_deeply_enough():
    payload = json.loads(DEPTH_JSON.read_text(encoding="utf-8"))

    assert payload["all_converged"], payload["summary"]
    assert len(payload["fractions"]) >= 3, "need several depths to judge convergence"
    assert payload["full_read_pairs"] > 1_000_000


@pytest.mark.skipif(not DEPTH_JSON.is_file(), reason="run `openbiota depth-check` first")
def test_fragment_counts_scale_with_depth():
    """A sanity check on the subsampling itself: counts must grow with depth."""
    payload = json.loads(DEPTH_JSON.read_text(encoding="utf-8"))
    by_panel: dict[str, list[tuple[int, int]]] = {}
    for point in payload["points"]:
        by_panel.setdefault(point["panel"], []).append(
            (point["read_pairs"], point["fragments"])
        )

    for panel, series in by_panel.items():
        series.sort()
        counts = [fragments for _, fragments in series]
        assert counts == sorted(counts), f"{panel} fragment count fell as depth rose"
        if counts[0] >= 20:
            # roughly linear: 16x the depth should give 8-32x the fragments
            ratio = counts[-1] / counts[0]
            depth_ratio = series[-1][0] / series[0][0]
            assert 0.5 < ratio / depth_ratio < 2.0, f"{panel} scaled non-linearly"
