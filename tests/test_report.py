"""Diagnostics, band labelling and report rendering."""

from __future__ import annotations

import dataclasses
import json

import pytest

from openbiota.qc import MateQC, SampleQC
from openbiota.report import (
    IMPLAUSIBLE_COPIES_PER_100,
    STANDING_CAVEATS,
    band_for,
    build_results_json,
    render_compact,
    render_summary,
    run_diagnostics,
)
from openbiota.tally import tally
from openbiota.taxonomy import build_profile

from .conftest import FIELDS, hit_line


def _qc(read_pairs: int = 1_000_000) -> SampleQC:
    mate = MateQC(
        label="R1",
        path="/tmp/x_1.fastq.gz",
        file_bytes=1234,
        gzipped=True,
        records=read_pairs,
        sampled_reads=1000,
        min_length=100,
        max_length=151,
        mean_length=148.4,
        total_bases_sampled=148_400,
        gc_percent=48.5,
        mean_quality=35.6,
        quality_alphabet="#,9J",
        n_percent=0.01,
        duplicate_fraction=0.02,
        first_read_id="LH00469:636:23VMTGLT4:8:1101:1003:20056",
        last_sampled_read_id="LH00469:636:23VMTGLT4:8:2498:52140:8019",
        estimated_bases=int(read_pairs * 148.4),
    )
    return SampleQC(
        mates=(mate, dataclasses.replace(mate, label="R2")),
        paired=True,
        ids_identical_between_mates=True,
        mate_suffixes_present=False,
        record_counts_match=True,
        notes=("Read ids are identical between mates.",),
    )


def _flat(text: str) -> str:
    """Collapse whitespace so assertions survive the report's line wrapping."""
    return " ".join(text.split())


def _tally(database, panel_set, urda: int, butyrate: int, rpob: int):
    lines = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(urda)]
    lines += [hit_line(f"B{i}", 3, "BUT", "P00004") for i in range(butyrate)]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(rpob)]
    return tally(
        streams=[("R1", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )


# --------------------------------------------------------------------------- #
# bands
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("value", "label"),
    [
        (None, "indeterminate"),
        (0.0, "not detected"),
        (0.001, "trace"),
        (0.99, "trace"),
        (1.0, "low"),
        (9.99, "low"),
        (10.0, "moderate"),
        (49.9, "moderate"),
        (50.0, "high"),
        (99.9, "high"),
        (100.0, "very high"),
        (10_000.0, "very high"),
    ],
)
def test_band_boundaries(value, label):
    assert band_for(value) == label


# --------------------------------------------------------------------------- #
# diagnostics
# --------------------------------------------------------------------------- #


def test_zero_rpob_raises_an_error_diagnostic(database, panel_set):
    result = _tally(database, panel_set, urda=50, butyrate=50, rpob=0)
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    codes = {d.code: d for d in diagnostics}
    assert codes["rpob-zero"].level == "error"
    assert "INDETERMINATE" in codes["rpob-zero"].message


def test_zero_positive_control_raises_an_error_diagnostic(database, panel_set):
    result = _tally(database, panel_set, urda=100, butyrate=0, rpob=2000)
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    codes = {d.code: d for d in diagnostics}
    assert codes["positive-control-zero"].level == "error"
    assert "PIPELINE IS BROKEN" in codes["positive-control-zero"].message


def test_healthy_run_has_no_error_diagnostics(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    assert not [d for d in diagnostics if d.level == "error"]
    assert {"rpob-yield-ok", "positive-control-ok"} <= {d.code for d in diagnostics}


def test_low_rpob_yield_warns(database, panel_set):
    result = _tally(database, panel_set, urda=10, butyrate=100, rpob=50)
    # 50 rpoB fragments across 1e6 pairs = 50 per million, below the floor
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    codes = {d.code for d in diagnostics}
    assert "rpob-low-yield" in codes


def test_implausible_copy_number_warns(database, panel_set):
    # 300 urdA fragments / 600 aa vs 20 rpoB / 1200 aa -> 3000 copies per 100
    result = _tally(database, panel_set, urda=300, butyrate=300, rpob=20)
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    implausible = [d for d in diagnostics if d.code.startswith("implausible-")]
    assert implausible
    urda = result.panel("urda")
    assert urda.copies_per_100_genomes is not None
    assert urda.copies_per_100_genomes > IMPLAUSIBLE_COPIES_PER_100


def test_unstable_panel_is_flagged(database, panel_set):
    result = _tally(database, panel_set, urda=5, butyrate=400, rpob=2000)
    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))

    assert "unstable-urda" in {d.code for d in diagnostics}
    assert result.panel("urda").stable is False
    assert result.panel("urda").detected is True


def test_unresolved_subject_ids_warn(database, panel_set):
    lines = [hit_line("A", 999, "GHOST", "P9") for _ in range(3)]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]
    lines += [hit_line(f"B{i}", 3, "BUT", "P00004") for i in range(400)]
    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )

    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))
    assert "unresolved-subjects" in {d.code for d in diagnostics}


def test_low_identity_capture_warns(database, panel_set):
    """Many fragments at low identity means distant homologs, not carriage."""
    from openbiota.report import LOW_IDENTITY_CAPTURE, low_identity_capture

    lines = [
        hit_line(f"U{i}", 0, "URDA", "P00001", pident=67.0) for i in range(200)
    ]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]
    lines += [hit_line(f"B{i}", 3, "BUT", "P00004") for i in range(400)]
    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )

    gene = result.panel("urda").targets[0]
    assert gene.mean_identity is not None
    assert gene.mean_identity < LOW_IDENTITY_CAPTURE
    assert low_identity_capture(gene) is True

    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))
    hit = next(d for d in diagnostics if d.code == "low-identity-urda:URDA")
    assert hit.level == "warn"
    assert "distant-homolog capture" in hit.message
    # urdA is in the panel aggregate, so the doubt must be propagated
    assert "inherits the same doubt" in hit.message


def test_high_identity_does_not_warn(database, panel_set):
    from openbiota.report import low_identity_capture

    result = _tally(database, panel_set, urda=200, butyrate=400, rpob=2000)

    gene = result.panel("urda").targets[0]
    assert gene.mean_identity == pytest.approx(90.0)
    assert low_identity_capture(gene) is False

    diagnostics = run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels))
    assert not [d for d in diagnostics if d.code.startswith("low-identity-")]


def test_low_identity_needs_enough_fragments(database, panel_set):
    """A handful of low-identity hits is noise, not a pattern worth flagging."""
    from openbiota.report import LOW_IDENTITY_MIN_FRAGMENTS, low_identity_capture

    lines = [
        hit_line(f"U{i}", 0, "URDA", "P00001", pident=67.0)
        for i in range(LOW_IDENTITY_MIN_FRAGMENTS - 1)
    ]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]
    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )

    assert low_identity_capture(result.panel("urda").targets[0]) is False


def test_low_identity_flag_appears_in_the_table(database, panel_set):
    lines = [hit_line(f"U{i}", 0, "URDA", "P00001", pident=67.0) for i in range(200)]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]
    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )
    text = render_summary(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=None,
        profile=None,
        diagnostics=[],
        run_meta={},
    )

    assert "LOW-ID" in text
    assert "distant-homolog capture rather than carriage" in _flat(text)


def test_panel_without_own_decoys_explains_cross_panel_competition(database, panel_set):
    """butyrate in the fixture declares no decoys."""
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    text = render_summary(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=None,
        profile=None,
        diagnostics=[],
        run_meta={},
    )

    assert "cross-panel competition" in _flat(text)
    assert "not that no competition took place" in _flat(text)


def test_diagnostics_without_qc_do_not_crash(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    diagnostics = run_diagnostics(tally=result, qc=None, selected=list(result.panels))
    assert isinstance(diagnostics, list)


# --------------------------------------------------------------------------- #
# rendering
# --------------------------------------------------------------------------- #


def test_summary_states_the_hard_constraints(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    text = render_summary(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=_qc(),
        profile=build_profile(
            phylum_counts=result.normalizer.phylum_counts,
            genus_counts=result.normalizer.genus_counts,
            organism_counts=result.normalizer.organism_counts,
        ),
        diagnostics=run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels)),
        run_meta={"threads": 8},
    )

    flat = _flat(text)
    # Part A's constraint: capacity rather than concentration.
    assert "genetic capacity, not metabolite concentrations" in flat
    # Part B's constraint: resemblance rather than diagnosis.
    assert "not a diagnosis and not a probability of disease" in flat
    assert "ARBITRARY HEURISTICS" in flat
    assert "NOT taxonomic identification" in flat
    assert "UPPER BOUND" in flat  # urda declares one
    # Limitations live in one document rather than being repeated throughout.
    assert "docs/VALIDATION.md" in flat
    # every panel appears
    for panel in result.panels:
        assert panel.panel.name in text


def test_caveats_are_consolidated_not_scattered():
    """One pointer, not a wall of disclaimers."""
    assert len(STANDING_CAVEATS) <= 3
    assert any("VALIDATION.md" in c for c in STANDING_CAVEATS)


def test_summary_reports_indeterminate_prominently(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=0)
    text = render_summary(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=None,
        profile=None,
        diagnostics=run_diagnostics(tally=result, qc=None, selected=list(result.panels)),
        run_meta={},
    )

    assert "INDETERMINATE" in text
    assert "indeterminate" in text


def test_summary_shows_mate_collapse_accounting(database, panel_set):
    lines = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(40)]
    result = tally(
        streams=[("R1", iter(lines)), ("R2", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )
    text = render_summary(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=None,
        profile=None,
        diagnostics=[],
        run_meta={},
    )

    assert "80" in text  # raw read hits
    assert "40 duplicate read hits removed" in _flat(text)


def test_compact_render_is_short_and_carries_the_caveat(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    text = render_compact(sample="TEST", tally=result, selected=list(result.panels))

    assert len(text.splitlines()) < 20
    assert "arbitrary heuristics" in _flat(text)
    assert "genetic capacity, not metabolites" in _flat(text)


def test_results_json_is_serialisable_and_complete(database, panel_set):
    result = _tally(database, panel_set, urda=60, butyrate=400, rpob=2000)
    payload = build_results_json(
        sample="TEST",
        tally=result,
        selected=list(result.panels),
        qc=_qc(),
        profile=build_profile(
            phylum_counts=result.normalizer.phylum_counts,
            genus_counts=result.normalizer.genus_counts,
            organism_counts=result.normalizer.organism_counts,
        ),
        diagnostics=run_diagnostics(tally=result, qc=_qc(), selected=list(result.panels)),
        run_meta={"threads": "8"},
        reference_entries={k: v.to_json() for k, v in database.entries.items()},
    )

    text = json.dumps(payload)  # must not raise
    assert len(text) > 500
    assert payload["validated_reference_range"] is False
    assert payload["qualitative_bands_are_arbitrary"] is True
    assert payload["measures"].startswith("genetic capacity")
    assert payload["standing_caveats"] == list(STANDING_CAVEATS)

    urda = next(p for p in payload["panels"] if p["name"] == "urda")
    assert urda["accepted_fragments"] == 60
    #   urdA: 60 / 600 aa   = 0.1
    #   rpoB: 2000 / 1200 aa = 1.6667
    #   0.1 / 1.6667 * 100   = 6.0
    assert urda["copies_per_100_genomes"] == pytest.approx(6.0)
    assert urda["upper_bound_reason"]
    assert urda["genes"][0]["residue_check"] is not None
    assert payload["normalisation"]["fragments"] == 2000
    assert "formula" in payload["normalisation"]

    # butyrate: 400 / 450 aa = 0.8889; 0.8889 / 1.6667 * 100 = 53.333
    butyrate = next(p for p in payload["panels"] if p["name"] == "butyrate")
    assert butyrate["copies_per_100_genomes"] == pytest.approx(53.3333, abs=1e-3)
