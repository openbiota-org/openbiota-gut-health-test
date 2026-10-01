"""Input discovery, subsampling, command construction and caching."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from openbiota.errors import InputError, SearchError
from openbiota.net import Response, next_page_url
from openbiota.search import (
    BASE_FIELDS,
    RESIDUE_FIELDS,
    SearchConfig,
    _manifest_payload,
    build_command,
    count_lines,
    discover_sample,
    discover_samples,
    output_fields,
    sample_name_from,
    write_subsample,
)

READ = "@{name}\nACGTACGTAC\n+\nJJJJJJJJJJ\n"


def _write_fastq(path: Path, n: int, *, gzipped: bool = False, prefix: str = "READ") -> Path:
    body = "".join(READ.format(name=f"{prefix}{i}") for i in range(n))
    if gzipped:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(body)
    else:
        path.write_text(body, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# discovery
# --------------------------------------------------------------------------- #


def test_discovers_underscore_numbered_pair(tmp_path: Path):
    _write_fastq(tmp_path / "A02_1.fastq", 4)
    _write_fastq(tmp_path / "A02_2.fastq", 4)

    sample = discover_sample(fastq_dir=tmp_path)

    assert sample.sample == "A02"
    assert [m.label for m in sample.mates] == ["R1", "R2"]
    assert sample.mates[0].path.name == "A02_1.fastq"


def test_discovers_r1_r2_naming(tmp_path: Path):
    _write_fastq(tmp_path / "donor_R1.fq", 4)
    _write_fastq(tmp_path / "donor_R2.fq", 4)

    sample = discover_sample(fastq_dir=tmp_path)

    assert sample.sample == "donor"
    assert len(sample.mates) == 2


def test_prefers_gzip_when_both_forms_exist(tmp_path: Path):
    _write_fastq(tmp_path / "S_1.fastq", 4)
    _write_fastq(tmp_path / "S_1.fastq.gz", 4, gzipped=True)
    _write_fastq(tmp_path / "S_2.fastq", 4)
    _write_fastq(tmp_path / "S_2.fastq.gz", 4, gzipped=True)

    sample = discover_sample(fastq_dir=tmp_path)

    assert all(m.path.suffix == ".gz" for m in sample.mates)
    assert all(m.gzipped for m in sample.mates)


def test_uncompressed_can_be_forced(tmp_path: Path):
    _write_fastq(tmp_path / "S_1.fastq", 4)
    _write_fastq(tmp_path / "S_1.fastq.gz", 4, gzipped=True)

    sample = discover_sample(fastq_dir=tmp_path, prefer_gzip=False)

    assert sample.mates[0].path.suffix == ".fastq"


def test_single_unpaired_file_is_accepted(tmp_path: Path):
    _write_fastq(tmp_path / "solo.fastq", 4)

    sample = discover_sample(fastq_dir=tmp_path)

    assert len(sample.mates) == 1
    assert sample.sample == "solo"


def test_empty_directory_fails_clearly(tmp_path: Path):
    with pytest.raises(InputError, match="no FASTQ files found"):
        discover_sample(fastq_dir=tmp_path)


def test_missing_directory_fails_clearly(tmp_path: Path):
    with pytest.raises(InputError, match="directory not found"):
        discover_sample(fastq_dir=tmp_path / "absent")


def test_multiple_samples_require_disambiguation(tmp_path: Path):
    _write_fastq(tmp_path / "a_1.fastq", 4)
    _write_fastq(tmp_path / "a_2.fastq", 4)
    _write_fastq(tmp_path / "b_1.fastq", 4)
    _write_fastq(tmp_path / "b_2.fastq", 4)

    with pytest.raises(InputError, match="contains several samples"):
        discover_sample(fastq_dir=tmp_path)

    sample = discover_sample(fastq_dir=tmp_path, sample="b")
    assert sample.sample == "b"


def test_explicit_paths(tmp_path: Path):
    r1 = _write_fastq(tmp_path / "x_1.fastq", 4)
    r2 = _write_fastq(tmp_path / "x_2.fastq", 4)

    sample = discover_sample(r1=r1, r2=r2, sample="mine")

    assert sample.sample == "mine"
    assert sample.paths == (r1.resolve(), r2.resolve())


def test_explicit_missing_path_fails(tmp_path: Path):
    with pytest.raises(InputError, match="input file not found"):
        discover_sample(r1=tmp_path / "nope.fastq")


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("A02_1.fastq", "A02"),
        ("A02_2.fastq.gz", "A02"),
        ("donor_R1.fq.gz", "donor"),
        ("plain.fastq", "plain"),
        # the whole filename is the sample; prefixes are never parsed away
        ("SAMPLE6_A06_1.fastq.gz", "SAMPLE6_A06"),
        ("SAMPLE4_A04_R2.fastq.gz", "SAMPLE4_A04"),
        ("SAMPLE2_A02_1.fastq.gz", "SAMPLE2_A02"),
        ("SRR33789733_1.fastq.gz", "SRR33789733"),
        ("renamed.anything_goes-here_R1.fastq.gz", "renamed.anything_goes-here"),
    ],
)
def test_sample_name_from(name: str, expected: str):
    assert sample_name_from(Path(name)) == expected


def test_discovery_keeps_the_full_filename_as_the_sample(tmp_path: Path):
    for name in ("SAMPLE6_A06_1.fastq.gz", "SAMPLE6_A06_2.fastq.gz", "SAMPLE2_A02_1.fastq.gz", "SAMPLE2_A02_2.fastq.gz"):
        (tmp_path / name).write_bytes(b"")
    found = discover_samples(tmp_path)
    assert [s.sample for s in found] == ["SAMPLE2_A02", "SAMPLE6_A06"]
    assert all(len(s.mates) == 2 for s in found)


# --------------------------------------------------------------------------- #
# subsampling
# --------------------------------------------------------------------------- #


def test_subsample_writes_exactly_n_reads(tmp_path: Path, reporter):
    source = _write_fastq(tmp_path / "in.fastq", 1000)
    dest = tmp_path / "out.fastq"

    write_subsample(source, dest, 25, reporter)

    assert count_lines(dest) == 100
    assert dest.read_text(encoding="utf-8").startswith("@READ0\n")


def test_subsample_reads_gzip(tmp_path: Path, reporter):
    source = _write_fastq(tmp_path / "in.fastq.gz", 500, gzipped=True)
    dest = tmp_path / "out.fastq"

    write_subsample(source, dest, 10, reporter)

    assert count_lines(dest) == 40


def test_subsample_shorter_than_requested_is_fine(tmp_path: Path, reporter):
    source = _write_fastq(tmp_path / "in.fastq", 5)
    dest = tmp_path / "out.fastq"

    write_subsample(source, dest, 1000, reporter)

    assert count_lines(dest) == 20


def test_subsample_is_cached(tmp_path: Path, reporter):
    source = _write_fastq(tmp_path / "in.fastq", 100)
    dest = tmp_path / "out.fastq"

    write_subsample(source, dest, 10, reporter)
    dest.write_text("SENTINEL", encoding="utf-8")
    write_subsample(source, dest, 10, reporter)

    assert dest.read_text(encoding="utf-8") == "SENTINEL"


def test_subsample_rejects_truncated_fastq(tmp_path: Path, reporter):
    source = tmp_path / "bad.fastq"
    source.write_text("@READ0\nACGT\n+\n", encoding="utf-8")

    with pytest.raises(InputError, match="not a multiple of 4"):
        write_subsample(source, tmp_path / "out.fastq", 10, reporter)


def test_subsample_rejects_empty_input(tmp_path: Path, reporter):
    source = tmp_path / "empty.fastq"
    source.write_text("", encoding="utf-8")

    with pytest.raises(InputError, match="is empty"):
        write_subsample(source, tmp_path / "out.fastq", 10, reporter)


# --------------------------------------------------------------------------- #
# command construction
# --------------------------------------------------------------------------- #


#: `build_command` resolves the executable, so these need DIAMOND on PATH; a
#: clean checkout without it skips them rather than failing on the dependency.
needs_diamond = pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond is not installed")


def test_output_fields_toggle_residue_columns():
    assert output_fields(with_residues=False) == BASE_FIELDS
    assert output_fields(with_residues=True) == BASE_FIELDS + RESIDUE_FIELDS
    # "Do not request fields you will not parse."
    assert "qseq_gapped" not in output_fields(with_residues=False)


@needs_diamond
def test_command_carries_the_required_flags(tmp_path: Path):
    config = SearchConfig(threads=24, block_size=8.0, index_chunks=1, sensitivity="default")
    cmd = build_command(
        config=config,
        db_path=tmp_path / "db.dmnd",
        query_path=tmp_path / "q.fastq.gz",
        out_path=tmp_path / "out.tsv",
        fields=BASE_FIELDS,
    )
    text = " ".join(cmd)

    assert "blastx" in cmd
    assert "--max-target-seqs 1" in text  # best hit only
    assert "--index-chunks 1" in text  # single in-memory chunk
    assert "--block-size 8" in text
    assert "--threads 24" in text
    assert "--quiet" in cmd
    assert "--outfmt 6 qseqid sseqid pident length sstart send bitscore" in text
    # default sensitivity means no sensitivity flag at all
    assert not any(c.endswith("-sensitive") for c in cmd)


@needs_diamond
def test_sensitivity_flag_is_applied(tmp_path: Path):
    cmd = build_command(
        config=SearchConfig(sensitivity="very-sensitive"),
        db_path=tmp_path / "db.dmnd",
        query_path=tmp_path / "q.fastq",
        out_path=tmp_path / "o.tsv",
        fields=BASE_FIELDS,
    )
    assert "--very-sensitive" in cmd


def test_unknown_sensitivity_raises():
    with pytest.raises(SearchError, match="unknown sensitivity"):
        SearchConfig(sensitivity="hyper").sensitivity_args()


def test_count_lines(tmp_path: Path):
    path = tmp_path / "f.txt"
    path.write_text("a\nb\nc\n", encoding="utf-8")
    assert count_lines(path) == 3


@needs_diamond
def test_cache_manifest_ignores_thread_count_and_paths(tmp_path: Path):
    """A laptop and a workstation screening the same reads share one cache.

    Thread count changes how fast DIAMOND runs, not what it writes, so it must
    not be part of the cache key. Absolute paths vary per machine and per
    checkout, so neither must they. Anything that changes the output — here,
    the sensitivity mode — still must.
    """
    query = _write_fastq(tmp_path / "q.fastq", 3)

    def payload(threads: int, sensitivity: str = "default") -> dict[str, object]:
        cmd = build_command(
            config=SearchConfig(threads=threads, sensitivity=sensitivity),
            db_path=tmp_path / f"db_{threads}.dmnd",
            query_path=query,
            out_path=tmp_path / f"out_{threads}.tsv",
            fields=BASE_FIELDS,
        )
        return _manifest_payload(
            command=cmd, fields=BASE_FIELDS, db_fingerprint="abc", query_path=query, subsample=None
        )

    assert payload(threads=14) == payload(threads=32)
    assert "--threads" not in payload(threads=14)["command"]
    assert payload(threads=14) != payload(threads=14, sensitivity="very-sensitive")


# --------------------------------------------------------------------------- #
# HTTP helper
# --------------------------------------------------------------------------- #


def test_link_header_with_commas_inside_the_url():
    """UniProt's `fields=` parameter contains commas; splitting on ',' is wrong."""
    url = (
        "https://rest.uniprot.org/uniprotkb/search?format=tsv"
        "&fields=accession,protein_name,gene_primary,organism_name,lineage,length,sequence"
        "&query=x&cursor=abc&size=500"
    )
    response = Response(url="x", status=200, body=b"", headers={"link": f"<{url}>; rel=\"next\""})
    assert next_page_url(response) == url


def test_link_header_absent():
    assert next_page_url(Response(url="x", status=200, body=b"", headers={})) is None


def test_link_header_without_next_relation():
    response = Response(
        url="x", status=200, body=b"", headers={"Link": '<https://e.com/a>; rel="prev"'}
    )
    assert next_page_url(response) is None
