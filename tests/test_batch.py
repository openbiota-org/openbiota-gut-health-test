"""Batch discovery and the cross-sample comparison."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from openbiota.batch import (
    SampleOutcome,
    comparison_json,
    outcome_from_results,
    print_comparison,
    write_comparison,
)
from openbiota.errors import InputError
from openbiota.search import discover_samples


def _fastq(path: Path, *, records: int = 2, gzipped: bool = False) -> Path:
    body = "".join(f"@r{i}/1\nACGTACGT\n+\nIIIIIIII\n" for i in range(records))
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzipped:
        with gzip.open(path, "wt") as handle:
            handle.write(body)
    else:
        path.write_text(body, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# discovery
# --------------------------------------------------------------------------- #


def test_discovers_every_pair(tmp_path: Path) -> None:
    for name in ("SAMPLE2", "ALPHA", "ZULU"):
        _fastq(tmp_path / f"{name}_1.fastq.gz", gzipped=True)
        _fastq(tmp_path / f"{name}_2.fastq.gz", gzipped=True)

    samples = discover_samples(tmp_path)
    assert [s.sample for s in samples] == ["ALPHA", "SAMPLE2", "ZULU"]
    assert all(len(s.mates) == 2 for s in samples)


def test_sample_name_comes_from_the_filename(tmp_path: Path) -> None:
    """Renaming a file renames its report and its output directory."""
    _fastq(tmp_path / "Patient_A_R1.fastq.gz", gzipped=True)
    _fastq(tmp_path / "Patient_A_R2.fastq.gz", gzipped=True)
    (sample,) = discover_samples(tmp_path)
    assert sample.sample == "Patient_A"


@pytest.mark.parametrize(
    "names",
    [
        ("S_1.fastq.gz", "S_2.fastq.gz"),
        ("S_R1.fastq.gz", "S_R2.fastq.gz"),
        ("S_1.fq.gz", "S_2.fq.gz"),
        ("S_1.fastq", "S_2.fastq"),
    ],
)
def test_mate_naming_conventions(tmp_path: Path, names: tuple[str, str]) -> None:
    for name in names:
        _fastq(tmp_path / name, gzipped=name.endswith(".gz"))
    (sample,) = discover_samples(tmp_path)
    assert sample.sample == "S"
    assert len(sample.mates) == 2


def test_gzip_is_preferred_when_both_forms_exist(tmp_path: Path) -> None:
    _fastq(tmp_path / "S_1.fastq")
    _fastq(tmp_path / "S_1.fastq.gz", gzipped=True)
    _fastq(tmp_path / "S_2.fastq")
    _fastq(tmp_path / "S_2.fastq.gz", gzipped=True)
    (sample,) = discover_samples(tmp_path)
    assert all(m.path.suffix == ".gz" for m in sample.mates)


def test_unpaired_files_are_ignored_not_guessed(tmp_path: Path) -> None:
    """A stray FASTQ must not silently become a single-end sample."""
    _fastq(tmp_path / "S_1.fastq.gz", gzipped=True)
    _fastq(tmp_path / "S_2.fastq.gz", gzipped=True)
    _fastq(tmp_path / "stray.fastq.gz", gzipped=True)
    samples = discover_samples(tmp_path)
    assert [s.sample for s in samples] == ["S"]


def test_hidden_files_are_skipped(tmp_path: Path) -> None:
    _fastq(tmp_path / "S_1.fastq.gz", gzipped=True)
    _fastq(tmp_path / "S_2.fastq.gz", gzipped=True)
    _fastq(tmp_path / "._S_1.fastq.gz", gzipped=True)
    (sample,) = discover_samples(tmp_path)
    assert not any(m.path.name.startswith(".") for m in sample.mates)


def test_empty_directory_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="no mate pairs found"):
        discover_samples(tmp_path)


def test_missing_directory_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="not found"):
        discover_samples(tmp_path / "nope")


# --------------------------------------------------------------------------- #
# outcome extraction
# --------------------------------------------------------------------------- #


#: Shaped exactly as ``results.json`` is written — ``panels`` is a list keyed
#: by "name", and ``reference_comparison.panels`` is a *mapping* keyed by panel
#: name. Getting that second one wrong made a whole batch report every sample
#: as failed while the reports themselves were fine, so the fixture mirrors the
#: real file rather than a convenient approximation of it.
RESULTS = {
    "input": {"read_pairs": 8_029_909},
    "panels": [
        {"name": "butyrate", "copies_per_100_genomes": 82.8},
        {"name": "urda", "copies_per_100_genomes": 35.1},
    ],
    "reference_comparison": {
        "panels": {
            "butyrate": {"value": 82.8, "percentile": 71.0},
            "urda": {"value": 35.1, "percentile": 30.0},
        }
    },
    "profile_similarity": {
        "taxonomic_engine": {
            "n_species_detected": 92,
            "host_filter": {"host_fraction": 0.0},
        },
        "dysbiosis_anchor": {"score": -1.35, "band": "strongly dysbiotic"},
        "profiles": {
            "crc": {
                "combined": {"percentile": 23.0},
                "confidence": {"grade": "MODERATE"},
                "abstention": {"abstained": False},
                "modules": {
                    "ecological": {
                        "features": [{"name": "shannon_diversity", "raw_value": 3.17}]
                    }
                },
            },
            "mecfs": {
                "combined": {"percentile": 72.0},
                "confidence": {"grade": "INSUFFICIENT"},
                "abstention": {"abstained": True},
            },
        },
    },
}


def test_outcome_pulls_the_headline_numbers(tmp_path: Path) -> None:
    outcome = outcome_from_results("SAMPLE2", tmp_path, 12.5, RESULTS)
    assert outcome.ok
    assert outcome.read_pairs == 8_029_909
    assert outcome.panels["butyrate"] == pytest.approx(82.8)
    assert outcome.percentiles["urda"] == pytest.approx(30.0)
    assert outcome.n_species == 92
    assert outcome.gmwi2 == pytest.approx(-1.35)
    assert outcome.gmwi2_band == "strongly dysbiotic"
    assert outcome.host_fraction == pytest.approx(0.0)
    assert outcome.shannon == pytest.approx(3.17)


def test_abstained_profile_reports_no_percentile(tmp_path: Path) -> None:
    outcome = outcome_from_results("SAMPLE2", tmp_path, 1.0, RESULTS)
    assert outcome.profiles["crc"] == pytest.approx(23.0)
    assert outcome.profiles["mecfs"] is None
    assert outcome.profile_confidence["mecfs"] == "abstained"


def test_reference_comparison_may_also_be_a_list(tmp_path: Path) -> None:
    """Tolerate the row-list shape as well, so the reader is not brittle."""
    results = {
        **RESULTS,
        "reference_comparison": {
            "panels": [{"panel": "butyrate", "percentile": 64.0}]
        },
    }
    outcome = outcome_from_results("X", tmp_path, 1.0, results)
    assert outcome.percentiles["butyrate"] == pytest.approx(64.0)


def test_outcome_matches_a_real_results_file(tmp_path: Path) -> None:
    """Guard against schema drift using a committed run, when one exists."""
    from openbiota.samples import results_file

    real = results_file("SAMPLE2_A02")
    if not real.is_file():
        pytest.skip("no SAMPLE2_A02 run present; `make run` first")
    outcome = outcome_from_results(
        "SAMPLE2", tmp_path, 1.0, json.loads(real.read_text(encoding="utf-8"))
    )
    assert outcome.ok
    assert outcome.read_pairs and outcome.read_pairs > 0
    assert outcome.panels, "no panel values extracted"
    assert outcome.percentiles, "no percentiles extracted"
    assert all(v is not None for v in outcome.panels.values())


def test_outcome_survives_a_sparse_results_file(tmp_path: Path) -> None:
    """A run with profiles disabled still yields a usable comparison row."""
    outcome = outcome_from_results("X", tmp_path, 1.0, {"input": {"read_pairs": 10}})
    assert outcome.ok
    assert outcome.read_pairs == 10
    assert outcome.gmwi2 is None
    assert outcome.profiles == {}


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #


def _outcomes(tmp_path: Path) -> list[SampleOutcome]:
    good = outcome_from_results("SAMPLE2", tmp_path / "SAMPLE2", 10.0, RESULTS)
    other = outcome_from_results("A03", tmp_path / "A03", 11.0, RESULTS)
    other.gmwi2 = 0.8
    other.gmwi2_band = "not dysbiotic"
    bad = SampleOutcome(
        sample="BROKEN", ok=False, out_dir=tmp_path / "BROKEN",
        elapsed_s=0.4, error="InputError: truncated FASTQ",
    )
    return [good, other, bad]


def test_comparison_json_records_failures(tmp_path: Path) -> None:
    payload = comparison_json(_outcomes(tmp_path))
    assert payload["n_samples"] == 3
    assert payload["n_succeeded"] == 2
    assert "butyrate" in payload["panels_compared"]
    failed = [s for s in payload["samples"] if not s["ok"]]
    assert failed and "truncated" in failed[0]["error"]


def test_comparison_is_written_to_disk(tmp_path: Path) -> None:
    path = write_comparison(tmp_path / "comparison.json", _outcomes(tmp_path))
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["n_succeeded"] == 2


def test_comparison_table_prints_every_section(tmp_path: Path, capsys) -> None:
    print_comparison(_outcomes(tmp_path))
    out = capsys.readouterr().out
    assert "SAMPLE2" in out
    assert "A03" in out
    assert "GENERAL GUT HEALTH" in out
    assert "METABOLITE PATHWAYS" in out
    assert "DISEASE-PATTERN SIMILARITY" in out
    assert "abstained" in out
    assert "FAILED  BROKEN" in out


def test_comparison_handles_no_successes(capsys) -> None:
    print_comparison(
        [SampleOutcome(sample="X", ok=False, out_dir=Path(), elapsed_s=0.0, error="boom")]
    )
    assert "no samples completed successfully" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# the run log carries the total, and nobody else writes it
# --------------------------------------------------------------------------- #

ROOT = Path(__file__).resolve().parents[1]
RESULT_LOGS = sorted((ROOT / "results").glob("*/run.log")) if (ROOT / "results").is_dir() else []


def test_the_parent_never_opens_the_childs_run_log() -> None:
    """Two writers on one file left every batch run.log without a single timing.

    The child writes ``run.log``; the parent captures console output beside
    it. This pins the parent's file name so the clobber cannot come back.
    """
    import inspect

    from openbiota import cli

    source = inspect.getsource(cli._run_all_parallel)
    assert 'out_dir / "run.console.log"' in source
    assert 'out_dir / "run.log"' not in source, "the batch parent must not write the child's run.log"


def test_the_total_is_emitted_before_the_log_is_flushed() -> None:
    """The final 'done in' line must be buffered before write_log, or it never
    reaches the file — which is exactly what used to happen."""
    import inspect

    from openbiota import cli

    source = inspect.getsource(cli.cmd_run)
    done_at = source.index('f"done in {human_duration(total)}')
    flush_at = source.index("reporter.write_log(out_dir")
    assert done_at < flush_at, "the total must be emitted before the log buffer is written"
    # And it must be the true total, measured after the PDF, not the pre-PDF `elapsed`.
    assert 'human_duration(elapsed)} — ' not in source


@pytest.mark.skipif(not RESULT_LOGS, reason="no generated results to inspect")
def test_every_generated_run_log_ends_with_its_total() -> None:
    """On real output: last line is the total, header carries it, JSON carries it."""
    import json
    import re

    for log in RESULT_LOGS:
        text = log.read_text(encoding="utf-8")
        lines = [ln for ln in text.splitlines() if ln.strip()]
        assert re.search(r"^total wall clock: \S+", text, re.M), f"{log}: no total in header"
        assert re.match(r"^\[\s*[\d.]+s\] .*done in \S+ — ", lines[-1]), (
            f"{log}: last line is not the total: {lines[-1][:80]!r}"
        )
        run = json.loads((log.parent / "results.json").read_text(encoding="utf-8"))["run"]
        assert isinstance(run.get("total_s"), (int, float)) and run["total_s"] > 0
        assert run["total wall clock"]
