"""How the reference cohort chooses its reads, and why that is now measured.

The reference cohort used to be built from the first 600,000 read pairs of
each run. That was justified by "every figure is a ratio to rpoB, so depth
divides out". Measured on one sample screened three ways at the same
600k budget, against its full-depth answer:

    first 600k    median error 19.6%, 16 of 24 panels >10% off
    random 600k   median error 11.6%, 13 of 24 panels >10% off,
                  and every deviation inside Poisson counting error
                  (median |z| 0.87, none beyond 3σ, sign split 16/7)

So the head of a FASTQ is a biased sample of it — flowcell position changes
which fragments pass a fixed identity threshold — while depth, once sampling
is uniform, is only precision. These tests pin the consequences:

* sampling is uniform across the run and reproducible from a seed;
* mates stay paired;
* the ordering ENA returns is never trusted, so a cohort's membership does
  not depend on the minute it was built;
* every reference range records how it was made, and the active one was
  made the right way.
"""

from __future__ import annotations

import gzip
import http.server
import json
import random
import threading
from pathlib import Path

import pytest

from openbiota import fastcache
from openbiota.cohort import (
    DEFAULT_SEED,
    SAMPLE_WINDOWS,
    CohortRun,
    ReferenceRanges,
    _window_starts,
    build_ranges,
    sample_remote_fastq,
)

ROOT = Path(__file__).resolve().parents[1]


def _run(accession: str = "SRR000001", read_count: int = 32_000_000) -> CohortRun:
    return CohortRun(
        run=accession, sample="S", r1_url=f"https://x/{accession}_1.fastq.gz",
        r2_url=f"https://x/{accession}_2.fastq.gz", read_count=read_count,
    )


# --------------------------------------------------------------------------- #
# scattered windows
# --------------------------------------------------------------------------- #


def test_windows_are_spread_across_the_whole_run() -> None:
    run = _run(read_count=32_000_000)
    windows = _window_starts(run, 2_000_000, SAMPLE_WINDOWS, DEFAULT_SEED)
    assert len(windows) == SAMPLE_WINDOWS
    firsts = [a for a, _ in windows]
    lasts = [b for _, b in windows]
    # Not a prefix: the last window sits in the final few percent of the run.
    assert lasts[-1] > 0.9 * run.read_count
    assert firsts[0] < 0.1 * run.read_count
    # Evenly spaced, in order, non-overlapping.
    strides = [b - a for a, b in zip(firsts[:-1], firsts[1:], strict=True)]
    assert max(strides) - min(strides) <= 1
    assert all(lasts[i] < firsts[i + 1] for i in range(len(windows) - 1))
    # And they add up to the target.
    assert sum(b - a + 1 for a, b in windows) == pytest.approx(2_000_000, rel=0.01)


def test_windows_are_reproducible_and_run_specific() -> None:
    a = _window_starts(_run("SRR1"), 2_000_000, 30, DEFAULT_SEED)
    b = _window_starts(_run("SRR1"), 2_000_000, 30, DEFAULT_SEED)
    c = _window_starts(_run("SRR2"), 2_000_000, 30, DEFAULT_SEED)
    d = _window_starts(_run("SRR1"), 2_000_000, 30, DEFAULT_SEED + 1)
    assert a == b, "same run, same seed, same windows"
    assert a != c, "two runs must not sample the same spot positions"
    assert a != d, "a different seed must move the windows"


def test_a_shallow_run_is_taken_whole() -> None:
    run = _run(read_count=1_500_000)
    assert _window_starts(run, 2_000_000, 30, DEFAULT_SEED) == [(1, 1_500_000)]


# --------------------------------------------------------------------------- #
# streaming fallback: uniform, paired, and served from an HTTP server
# --------------------------------------------------------------------------- #


def _fastq(n: int, mate: int) -> bytes:
    rng = random.Random(7)
    out = []
    for i in range(n):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        out.append(f"@READ{i} {mate}\n{seq}\n+\n{'F' * 60}\n")
    return "".join(out).encode()


@pytest.fixture(scope="module")
def served_pair(tmp_path_factory: pytest.TempPathFactory) -> tuple[str, str, int]:
    """Two mates of a 4,000-record FASTQ, gzipped, behind a local HTTP server."""
    root = tmp_path_factory.mktemp("srv")
    n = 4_000
    (root / "r_1.fastq.gz").write_bytes(gzip.compress(_fastq(n, 1)))
    (root / "r_2.fastq.gz").write_bytes(gzip.compress(_fastq(n, 2)))

    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *_: object) -> None:
            return

    handler = lambda *a, **k: Quiet(*a, directory=str(root), **k)  # noqa: E731
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    yield f"{base}/r_1.fastq.gz", f"{base}/r_2.fastq.gz", n
    server.shutdown()


def test_streaming_sample_is_uniform_and_keeps_mates_paired(
    served_pair: tuple[str, str, int], tmp_path: Path
) -> None:
    url1, url2, n = served_pair
    out1, out2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"
    got1 = sample_remote_fastq(url1, out1, probability=0.25, seed=DEFAULT_SEED)
    got2 = sample_remote_fastq(url2, out2, probability=0.25, seed=DEFAULT_SEED)
    assert got1 == got2
    # Binomial(4000, 0.25): mean 1000, sd ~27. Five sd is a generous gate.
    assert abs(got1 - 1000) < 135, got1

    ids1 = [line.split()[0] for line in out1.read_text().splitlines()[::4]]
    ids2 = [line.split()[0] for line in out2.read_text().splitlines()[::4]]
    assert ids1 == ids2, "the same record indices must be kept in both mates"

    # Uniform, not a prefix: kept records span the file.
    index = [int(i[5:]) for i in ids1]
    assert index[0] < n * 0.05 and index[-1] > n * 0.95
    # Evenly spread: each quarter of the file contributes about a quarter.
    quarters = [sum(1 for i in index if q * n / 4 <= i < (q + 1) * n / 4) for q in range(4)]
    assert max(quarters) - min(quarters) < 0.4 * got1 / 4 + 60


def test_streaming_sample_is_cached_by_its_done_marker(
    served_pair: tuple[str, str, int], tmp_path: Path
) -> None:
    url1, _, _ = served_pair
    out = tmp_path / "c_1.fastq"
    first = sample_remote_fastq(url1, out, probability=0.1, seed=1)
    out.write_text("tampered")  # the marker, not the file, is the contract
    assert sample_remote_fastq(url1, out, probability=0.1, seed=1) == first


# --------------------------------------------------------------------------- #
# provenance
# --------------------------------------------------------------------------- #


def test_reference_ranges_record_how_they_were_built(tmp_path: Path) -> None:
    ranges = build_ranges(
        study="PRJNA0", description="d", reads_per_sample=2_000_000,
        panel_values={"butyrate": [1.0, 2.0, 3.0, 4.0, 5.0]},
        gene_values={"butyrate:BUT": [1.0, 2.0, 3.0, 4.0, 5.0]},
        sampling="random", seed=42, runs=["SRR3", "SRR1", "SRR2"], study_title="A study",
    )
    path = tmp_path / "r.json"
    ranges.save(path)
    back = ReferenceRanges.load(path)
    assert back.sampling == "random" and back.seed == 42
    assert back.runs == ("SRR3", "SRR1", "SRR2"), "order is the screening order"
    assert back.study_title == "A study"


def test_a_cohort_built_before_provenance_existed_is_a_prefix_cohort(tmp_path: Path) -> None:
    """Old files have no `sampling` key; they were all head-of-file cohorts."""
    legacy = {
        "study": "P", "description": "d", "n_samples": 5, "reads_per_sample": 600_000,
        "built_at": "t", "percentiles_reported": [5, 50, 95],
        "panels": {}, "genes": {},
    }
    path = tmp_path / "old.json"
    path.write_text(json.dumps(legacy))
    back = ReferenceRanges.load(path)
    assert back.sampling == "prefix" and back.seed is None and back.runs == ()


def test_ena_run_order_is_never_trusted(monkeypatch: pytest.MonkeyPatch) -> None:
    """ENA returns runs in no stable order; `list_runs` must sort them."""
    from openbiota import cohort as cohort_module

    tsv = (
        "run_accession\tsample_accession\tfastq_ftp\tread_count\n"
        "SRR9\tS\ta/1.gz;a/2.gz\t100\n"
        "SRR1\tS\tb/1.gz;b/2.gz\t100\n"
        "SRR5\tS\tc/1.gz;c/2.gz\t100\n"
    )

    class Response:
        def text(self) -> str:
            return tsv

    monkeypatch.setattr(cohort_module, "get", lambda *_a, **_k: Response())
    from openbiota.logging_util import Reporter

    runs = cohort_module.list_runs("P", Reporter(verbose=False))
    assert [r.run for r in runs] == ["SRR1", "SRR5", "SRR9"]


@pytest.mark.skipif(not (ROOT / "refs" / "reference_ranges.json").is_file(), reason="no built cohort")
def test_the_active_cohort_was_sampled_uniformly() -> None:
    """The regression this whole module exists to prevent."""
    active = ReferenceRanges.load(ROOT / "refs" / "reference_ranges.json")
    assert active.sampling == "random", (
        "refs/reference_ranges.json is a head-of-file cohort; rebuild with "
        "`openbiota cohort --sampling random`"
    )
    assert active.seed is not None
    assert len(active.runs) == active.n_samples >= 60
    assert active.reads_per_sample >= 2_000_000
    assert active.study_title, "the cohort should say who its population is"


# --------------------------------------------------------------------------- #
# the parsed-configuration cache
# --------------------------------------------------------------------------- #


def test_fastcache_self_test() -> None:
    assert fastcache.self_test() == 5


def test_cached_loaders_return_what_the_uncached_loaders_return(tmp_path: Path) -> None:
    """The cache must be invisible: same object either way."""
    from openbiota.panels import _load_panel_set_uncached, load_panel_set
    from openbiota.profiles import _load_profile_set_uncached, load_profile_set

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fastcache, "DEFAULT_CACHE_DIR", tmp_path)
        cached = load_panel_set(ROOT / "panels")
        cached_again = load_panel_set(ROOT / "panels")
    direct = _load_panel_set_uncached(ROOT / "panels")
    assert [p.name for p in cached.panels] == [p.name for p in direct.panels]
    assert [p.name for p in cached_again.panels] == [p.name for p in direct.panels]
    assert cached.panels[0].targets[0].key == direct.panels[0].targets[0].key

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fastcache, "DEFAULT_CACHE_DIR", tmp_path)
        profiles = load_profile_set(ROOT / "profiles")
    assert [p.name for p in profiles.profiles] == [
        p.name for p in _load_profile_set_uncached(ROOT / "profiles").profiles
    ]
    assert list(tmp_path.glob("*.pickle")), "the cache directory should have been populated"


def test_cache_invalidates_when_a_source_file_changes(tmp_path: Path) -> None:
    src = tmp_path / "x.yaml"
    src.write_text("a: 1")
    seen: list[int] = []

    def build() -> int:
        seen.append(1)
        return len(seen)

    assert fastcache.cached_load("t", [src], build, cache_dir=tmp_path / "c") == 1
    assert fastcache.cached_load("t", [src], build, cache_dir=tmp_path / "c") == 1
    src.write_text("a: 22")
    assert fastcache.cached_load("t", [src], build, cache_dir=tmp_path / "c") == 2
