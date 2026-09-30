"""Input validation and read statistics."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from openbiota.errors import InputError
from openbiota.qc import count_records, profile_mate, profile_sample


def _fastq(records: list[tuple[str, str, str]]) -> str:
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)


def _write(path: Path, records: list[tuple[str, str, str]], *, gzipped: bool = False) -> Path:
    body = _fastq(records)
    if gzipped:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(body)
    else:
        path.write_text(body, encoding="utf-8")
    return path


def _records(n: int, *, length: int = 100, prefix: str = "READ") -> list[tuple[str, str, str]]:
    return [(f"{prefix}{i}", "ACGT" * (length // 4), "J" * (length // 4 * 4)) for i in range(n)]


# --------------------------------------------------------------------------- #
# counting
# --------------------------------------------------------------------------- #


def test_count_records_plain(tmp_path: Path):
    path = _write(tmp_path / "a.fastq", _records(37))
    assert count_records(path) == 37


def test_count_records_gzip(tmp_path: Path):
    path = _write(tmp_path / "a.fastq.gz", _records(41), gzipped=True)
    assert count_records(path) == 41


def test_count_records_rejects_truncated(tmp_path: Path):
    path = tmp_path / "bad.fastq"
    path.write_text("@a\nACGT\n+\nJJJJ\n@b\nACGT\n", encoding="utf-8")
    with pytest.raises(InputError, match="not a multiple of 4"):
        count_records(path)


# --------------------------------------------------------------------------- #
# per-mate profiling
# --------------------------------------------------------------------------- #


def test_profile_reports_length_gc_and_quality(tmp_path: Path):
    records = [
        ("r0", "GGCC" * 25, "J" * 100),
        ("r1", "ATAT" * 25, "J" * 100),
    ]
    path = _write(tmp_path / "a.fastq", records)

    profile, ids = profile_mate(path, "R1", count_all=True)

    assert profile.records == 2
    assert profile.sampled_reads == 2
    assert profile.min_length == profile.max_length == 100
    assert profile.mean_length == pytest.approx(100.0)
    assert profile.gc_percent == pytest.approx(50.0)  # one all-GC, one all-AT
    assert profile.mean_quality == pytest.approx(ord("J") - 33)
    assert profile.quality_alphabet == "J"
    assert ids == ["r0", "r1"]


def test_profile_detects_variable_read_length(tmp_path: Path):
    records = [("r0", "ACGT" * 25, "J" * 100), ("r1", "ACGT" * 37, "J" * 148)]
    path = _write(tmp_path / "a.fastq", records)

    profile, _ = profile_mate(path, "R1", count_all=False)

    assert profile.min_length == 100
    assert profile.max_length == 148
    assert profile.records is None  # counting was skipped


def test_profile_counts_exact_duplicates(tmp_path: Path):
    records = [("r0", "ACGT" * 25, "J" * 100)] * 4
    path = _write(tmp_path / "a.fastq", records)

    profile, _ = profile_mate(path, "R1", count_all=False)

    assert profile.duplicate_fraction == pytest.approx(0.75)


def test_profile_honours_the_sample_cap(tmp_path: Path):
    path = _write(tmp_path / "a.fastq", _records(500))

    profile, ids = profile_mate(path, "R1", sample_reads=10, count_all=True)

    assert profile.sampled_reads == 10
    assert len(ids) == 10
    assert profile.records == 500  # full count is independent of the cap


def test_profile_rejects_length_mismatch(tmp_path: Path):
    path = tmp_path / "bad.fastq"
    path.write_text("@r\nACGTACGT\n+\nJJJ\n", encoding="utf-8")
    with pytest.raises(InputError, match="quality values"):
        profile_mate(path, "R1", count_all=False)


def test_profile_rejects_bad_header(tmp_path: Path):
    path = tmp_path / "bad.fastq"
    path.write_text("r\nACGT\n+\nJJJJ\n", encoding="utf-8")
    with pytest.raises(InputError, match="does not start with '@'"):
        profile_mate(path, "R1", count_all=False)


def test_profile_rejects_empty_file(tmp_path: Path):
    path = tmp_path / "empty.fastq"
    path.write_text("", encoding="utf-8")
    with pytest.raises(InputError, match="no FASTQ records"):
        profile_mate(path, "R1", count_all=False)


# --------------------------------------------------------------------------- #
# sample-level validation
# --------------------------------------------------------------------------- #


def test_identical_mate_ids_produce_the_mate_collapse_note(tmp_path: Path, reporter):
    """This dataset's ids are byte-identical between R1 and R2 with no suffix."""
    records = _records(20)
    r1 = _write(tmp_path / "s_1.fastq", records)
    r2 = _write(tmp_path / "s_2.fastq", records)

    qc = profile_sample([("R1", r1), ("R2", r2)], reporter=reporter)

    assert qc.paired is True
    assert qc.ids_identical_between_mates is True
    assert qc.mate_suffixes_present is False
    assert qc.record_counts_match is True
    assert qc.read_pairs == 20
    assert any("NOT independent observations" in n for n in qc.notes)
    assert any("inflated by up to 2x" in n for n in qc.notes)


def test_slash_suffixed_ids_are_reported_as_such(tmp_path: Path, reporter):
    r1 = _write(tmp_path / "s_1.fastq", [(f"r{i}/1", "ACGT" * 25, "J" * 100) for i in range(5)])
    r2 = _write(tmp_path / "s_2.fastq", [(f"r{i}/2", "ACGT" * 25, "J" * 100) for i in range(5)])

    qc = profile_sample([("R1", r1), ("R2", r2)], reporter=reporter)

    assert qc.mate_suffixes_present is True
    assert any("stripped during mate collapse" in n for n in qc.notes)


def test_mismatched_record_counts_are_reported(tmp_path: Path, reporter):
    r1 = _write(tmp_path / "s_1.fastq", _records(10))
    r2 = _write(tmp_path / "s_2.fastq", _records(8))

    qc = profile_sample([("R1", r1), ("R2", r2)], reporter=reporter)

    assert qc.record_counts_match is False
    assert qc.read_pairs == 8
    assert any("record counts differ" in n for n in qc.notes)


def test_binned_quality_is_reported_as_normal_not_corruption(tmp_path: Path, reporter):
    """NovaSeq emits 3-4 distinct quality characters; that is binning."""
    records = [("r0", "ACGT" * 25, ("J" * 60) + ("9" * 30) + ("#" * 10))]
    path = _write(tmp_path / "s_1.fastq", records)

    qc = profile_sample([("R1", path)], reporter=reporter)

    assert len(qc.mates[0].quality_alphabet) == 3
    assert any("quality binning" in n and "not corruption" in n for n in qc.notes)


def test_high_duplication_is_noted(tmp_path: Path, reporter):
    records = [("r0", "ACGT" * 25, "J" * 100)] * 10
    path = _write(tmp_path / "s_1.fastq", records)

    qc = profile_sample([("R1", path)], reporter=reporter)

    assert any("exact sequence duplicates" in n for n in qc.notes)
    assert any("rpoB is affected equally" in n for n in qc.notes)


def test_single_mate_sample_is_not_marked_paired(tmp_path: Path, reporter):
    path = _write(tmp_path / "solo.fastq", _records(6))

    qc = profile_sample([("R1", path)], reporter=reporter)

    assert qc.paired is False
    assert qc.ids_identical_between_mates is None
    assert qc.total_bases == 600


def test_counting_works_on_both_forms(tmp_path: Path):
    """Gzipped input is streamed through a pipe; plain input is read directly.

    Inputs are kept compressed and never decompressed to a file — see
    `GZIP_STREAMER` for the measurements behind that choice.
    """
    plain = _write(tmp_path / "s_1.fastq", _records(12))
    gz = _write(tmp_path / "t_1.fastq.gz", _records(12), gzipped=True)

    assert count_records(plain) == 12
    assert count_records(gz) == 12


def test_streaming_and_zlib_agree(tmp_path: Path):
    """The fast path must return exactly what the library path would."""
    import gzip as _gzip

    from openbiota.qc import open_bytes

    gz = _write(tmp_path / "u_1.fastq.gz", _records(37), gzipped=True)
    with open_bytes(gz) as handle:
        streamed = handle.read()
    with _gzip.open(gz, "rb") as handle:
        via_zlib = handle.read()
    assert streamed == via_zlib


def test_early_close_is_not_reported_as_corruption(tmp_path: Path):
    """Stopping mid-stream sends SIGPIPE to gzip; that is not a bad file.

    Regression: taking a subsample closes the pipe after the requested
    records, the decompressor exits non-zero, and treating that as corruption
    rejected a perfectly good input.
    """
    from openbiota.qc import open_bytes

    gz = _write(tmp_path / "w_1.fastq.gz", _records(50_000), gzipped=True)
    with open_bytes(gz, partial_ok=True) as handle:
        handle.readline()  # read one line out of 200,000 and walk away
    # no exception


def test_full_read_of_a_good_file_succeeds(tmp_path: Path):
    gz = _write(tmp_path / "x_1.fastq.gz", _records(100), gzipped=True)
    from openbiota.qc import open_bytes

    with open_bytes(gz) as handle:
        assert handle.read().count(b"\n") == 400


def test_truncated_gzip_is_rejected(tmp_path: Path):
    """A record count that is not a multiple of four means a damaged file."""
    from openbiota.errors import InputError

    bad = tmp_path / "v_1.fastq"
    bad.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n", encoding="utf-8")
    with pytest.raises(InputError, match="not a multiple of 4"):
        count_records(bad)


# --------------------------------------------------------------------------- #
# caching
# --------------------------------------------------------------------------- #


def test_profile_is_cached_and_reused(tmp_path: Path, reporter):
    from openbiota.qc import profile_sample_cached

    r1 = _write(tmp_path / "s_1.fastq", _records(30))
    cache = tmp_path / "cache"

    first = profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)
    assert first.cached is False
    assert first.read_pairs == 30
    assert list(cache.glob("input_profile_*.json"))

    second = profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)
    assert second.cached is True
    assert second.read_pairs == 30
    assert second.mates[0].mean_length == first.mates[0].mean_length
    assert second.notes == first.notes


def test_cache_is_invalidated_when_the_input_changes(tmp_path: Path, reporter):
    from openbiota.qc import profile_sample_cached

    r1 = _write(tmp_path / "s_1.fastq", _records(30))
    cache = tmp_path / "cache"
    profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)

    # replace the input with a different number of records
    _write(r1, _records(40))
    refreshed = profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)

    assert refreshed.cached is False
    assert refreshed.read_pairs == 40


def test_cache_can_be_forced_to_refresh(tmp_path: Path, reporter):
    from openbiota.qc import profile_sample_cached

    r1 = _write(tmp_path / "s_1.fastq", _records(30))
    cache = tmp_path / "cache"
    profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)

    forced = profile_sample_cached(
        [("R1", r1)], cache_dir=cache, reporter=reporter, refresh=True
    )
    assert forced.cached is False


def test_corrupt_cache_falls_back_to_recomputing(tmp_path: Path, reporter):
    from openbiota.qc import profile_sample_cached

    r1 = _write(tmp_path / "s_1.fastq", _records(30))
    cache = tmp_path / "cache"
    profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)
    next(cache.glob("input_profile_*.json")).write_text("{not json", encoding="utf-8")

    recovered = profile_sample_cached([("R1", r1)], cache_dir=cache, reporter=reporter)
    assert recovered.cached is False
    assert recovered.read_pairs == 30


def test_qc_json_is_complete(tmp_path: Path, reporter):
    r1 = _write(tmp_path / "s_1.fastq", _records(10))
    qc = profile_sample([("R1", r1)], reporter=reporter)

    payload = qc.to_json()

    assert payload["read_pairs"] == 10
    assert payload["mates"][0]["read_length_mean"] == 100.0
    assert "gc_percent" in payload["mates"][0]
