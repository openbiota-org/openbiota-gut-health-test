"""Building the human-read filter index without being killed for it.

bowtie2-build's memory grows with its thread count. Asked for 32 threads on
a 62 GB machine it was killed by the kernel ten minutes in, and the only
word on stderr was "Building a LARGE index". These tests pin the three
things that stop that happening again: the thread cap from the machine's
memory, the retry with fewer threads when the build is killed, and the
temporary prefix that keeps a killed build from leaving six truncated files
that pass for an index.
"""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import pytest

from openbiota.engines import metaphlan as mp
from openbiota.errors import OpenBiotaError
from openbiota.logging_util import Reporter

BOWTIE2_BUILD = shutil.which("bowtie2-build")
PHIX = Path(__file__).resolve().parents[1] / "refs" / "host" / "phix.fna.gz"


def test_threads_are_capped_by_memory_and_at_eight():
    assert mp.host_build_threads(32, int(192e9)) == 8, "past eight threads the build barely speeds up"
    assert mp.host_build_threads(32, int(62e9)) == 8
    assert mp.host_build_threads(32, int(16e9)) == 4
    assert mp.host_build_threads(32, int(8e9)) == 1
    assert mp.host_build_threads(32, None) == 8, "memory unknown: the safe fixed cap"
    assert mp.host_build_threads(2, int(192e9)) == 2, "never more than asked for"
    assert mp.host_build_threads(0, int(192e9)) == 1


def _tiny_genome(host_dir: Path) -> None:
    """A host.fna small enough to index in a second: PhiX if the checkout has
    it, otherwise a synthetic 20 kb sequence."""
    host_dir.mkdir(parents=True, exist_ok=True)
    if PHIX.is_file():
        import gzip

        with gzip.open(PHIX) as handle:
            (host_dir / "host.fna").write_bytes(handle.read())
    else:
        import random

        rng = random.Random(7)
        seq = "".join(rng.choice("ACGT") for _ in range(20_000))
        (host_dir / "host.fna").write_text(">synthetic\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")


def _tools(bowtie2_build: str) -> mp.MetaphlanTools:
    return mp.MetaphlanTools(metaphlan="metaphlan", bowtie2="bowtie2", bowtie2_build=bowtie2_build, version="3.1.0")


def _fake_builder(tmp_path: Path, script: str) -> str:
    path = tmp_path / "fake-bowtie2-build"
    path.write_text("#!/usr/bin/env bash\n" + script)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


@pytest.mark.skipif(BOWTIE2_BUILD is None, reason="bowtie2-build not installed")
def test_a_real_build_lands_under_the_final_name_with_nothing_left_over(tmp_path: Path, monkeypatch):
    host_dir = tmp_path / "host"
    _tiny_genome(host_dir)
    monkeypatch.setattr(mp, "platform_info", lambda: type("P", (), {"memory_bytes": int(64e9)})(), raising=False)
    mp.build_host_index(_tools(BOWTIE2_BUILD), host_dir, Reporter(verbose=False), threads=4)
    assert mp.host_index_ready(host_dir)
    assert not list(host_dir.glob("building.*")), "the temporary prefix is gone"
    assert not (host_dir / "host.fna").exists(), "the concatenated genome is removed after a successful build"


@pytest.mark.skipif(BOWTIE2_BUILD is None, reason="bowtie2-build not installed")
def test_a_build_killed_once_is_retried_with_fewer_threads_and_finishes(tmp_path: Path, monkeypatch, capsys):
    """The out-of-memory shape: the process dies by signal, having already
    written part of the index. The retry must start clean and succeed."""
    host_dir = tmp_path / "host"
    _tiny_genome(host_dir)
    marker = tmp_path / "first-call-happened"
    fake = _fake_builder(tmp_path, f'''
if [ ! -e "{marker}" ]; then
    touch "{marker}"
    # leave a truncated shard behind, as a killed build does
    printf 'partial' > "${{@: -1}}.1.bt2l"
    kill -9 $$
fi
exec "{BOWTIE2_BUILD}" "$@"
''')
    monkeypatch.setattr(mp, "platform_info", lambda: type("P", (), {"memory_bytes": int(62e9)})(), raising=False)
    mp.build_host_index(_tools(fake), host_dir, Reporter(verbose=True), threads=32)
    err = capsys.readouterr().err
    assert "8 threads rather than 32" in err, "the cap is announced"
    assert "killed (signal 9)" in err and "out-of-memory" in err and "retrying once with 2" in err
    assert mp.host_index_ready(host_dir)
    assert not list(host_dir.glob("building.*"))
    assert (host_dir / "host.1.bt2l").stat().st_size > 7, "the final shard is not the truncated one"


def test_a_build_that_keeps_failing_explains_itself_and_leaves_no_partial_index(tmp_path: Path, monkeypatch):
    host_dir = tmp_path / "host"
    _tiny_genome(host_dir)
    fake = _fake_builder(tmp_path, '''
printf 'partial' > "${@: -1}.1.bt2l"
echo "Building a LARGE index" >&2
exit 1
''')
    monkeypatch.setattr(mp, "platform_info", lambda: type("P", (), {"memory_bytes": int(62e9)})(), raising=False)
    with pytest.raises(OpenBiotaError) as excinfo:
        mp.build_host_index(_tools(fake), host_dir, Reporter(verbose=False), threads=8)
    message = str(excinfo.value)
    assert "exit 1" in message and "--skip-host-index" in message and "LARGE index" in message
    assert not mp.host_index_ready(host_dir)
    assert not list(host_dir.glob("building.*")), "nothing half-built is left to be mistaken for an index"
    assert (host_dir / "host.fna").exists(), "the genome stays so the retry needs no new download"


def test_a_killed_build_that_dies_again_is_reported_as_killed(tmp_path: Path, monkeypatch):
    host_dir = tmp_path / "host"
    _tiny_genome(host_dir)
    fake = _fake_builder(tmp_path, "kill -9 $$\n")
    monkeypatch.setattr(mp, "platform_info", lambda: type("P", (), {"memory_bytes": int(62e9)})(), raising=False)
    with pytest.raises(OpenBiotaError, match="killed by signal 9"):
        mp.build_host_index(_tools(fake), host_dir, Reporter(verbose=False), threads=8)
    assert not mp.host_index_ready(host_dir)


def test_leftovers_from_an_interrupted_build_do_not_count_as_an_index(tmp_path: Path):
    host_dir = tmp_path / "host"
    host_dir.mkdir()
    for s in ("1", "2", "3", "4", "rev.1", "rev.2"):
        (host_dir / f"building.{s}.bt2l").write_bytes(b"partial")
    assert not mp.host_index_ready(host_dir)


def test_python_can_see_the_memory_the_cap_uses():
    from openbiota.search import platform_info

    assert platform_info().memory_bytes > 0 or os.environ.get("CI"), "the cap needs the machine's memory"
