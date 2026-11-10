"""`openbiota selftest`: a lane is installed when it has run.

The installer used to call a lane installed when its files were on disk.
SingleM's were all there on a fresh Ubuntu machine and its first real run
died on ``hmmsearch: command not found``. The self-test runs every lane on
made-up reads through the pipeline's own lane code; these tests cover the
reads it makes, how it judges a lane, what it prints and exits with, and
that the installer ends by running it. The lanes themselves need 190 GB of
references and are exercised by the installer, not here.
"""

from __future__ import annotations

import gzip
import io
import os
import random
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from openbiota import cli, selftest
from openbiota.selftest import LaneOutcome

REPO = Path(__file__).resolve().parents[1]
INSTALL = REPO / "install.sh"


def _genomes(folder: Path, length: int = 30_000) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(7)
    for filename, tag, _ in selftest.GENOMES:
        sequence = "".join(rng.choice("ACGT") for _ in range(length))
        (folder / filename).write_text(f">{tag}\n{sequence}\n")
    return folder


# --------------------------------------------------------------------------- #
# the reads
# --------------------------------------------------------------------------- #


def test_reads_are_a_gzipped_pair_with_a_distinct_name_per_read(tmp_path: Path):
    r1, r2 = selftest.simulate_reads(tmp_path / "work", n_pairs=300, genomes_dir=_genomes(tmp_path / "g"))
    assert r1.name == "selftest_1.fastq.gz" and r2.name == "selftest_2.fastq.gz"
    assert r1.read_bytes()[:2] == b"\x1f\x8b", "gzipped, as the reads a real run hands the lanes are"
    with gzip.open(r1, "rt") as fh:
        names1 = [line.split()[0] for i, line in enumerate(fh) if i % 4 == 0]
    with gzip.open(r2, "rt") as fh:
        names2 = [line.split()[0] for i, line in enumerate(fh) if i % 4 == 0]
    assert len(names1) == 300 and names1 == names2
    # A name is cut at its first space; cut names that collide are counted
    # once by the tools that count read pairs (mOTUs saw nothing at all).
    assert len(set(names1)) == 300
    assert not list((tmp_path / "work").glob("*.fastq")), "the uncompressed copies are not left behind"


def test_read_name_tags_have_no_spaces():
    for _, tag, _ in selftest.GENOMES:
        assert re.fullmatch(r"[A-Za-z0-9_]+", tag), tag


def test_the_reads_are_made_again_in_a_kept_folder(tmp_path: Path):
    genomes = _genomes(tmp_path / "g")
    first = selftest.simulate_reads(tmp_path / "work", n_pairs=50, genomes_dir=genomes)
    second = selftest.simulate_reads(tmp_path / "work", n_pairs=50, genomes_dir=genomes)
    assert first == second and all(p.is_file() for p in second)


def test_missing_genomes_point_at_the_installer(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match=r"run \./install\.sh"):
        selftest.simulate_reads(tmp_path / "work", genomes_dir=tmp_path / "nothing-here")


# --------------------------------------------------------------------------- #
# judging one lane
# --------------------------------------------------------------------------- #


def _lane_status(monkeypatch, status: dict, lineages: tuple[str, ...] = ()) -> None:
    from openbiota.expansion import lanes

    observations = [{"lineage": lineage, "rank": "species"} for lineage in lineages]
    monkeypatch.setattr(lanes, "run_all", lambda **kw: {
        "lane_status": {kw["only"][0]: status}, "lanes": {kw["only"][0]: {"observations": observations}}})


def _run_lane(tmp_path: Path, lane: str = "singlem_globdb") -> LaneOutcome:
    return selftest.run_lane(lane, r1=tmp_path / "r1", r2=tmp_path / "r2", refs_dir=tmp_path, work=tmp_path, threads=2)


@pytest.mark.parametrize("lineages", [
    # the same three genomes, as four catalogues spell them
    ("k__Bacteria|g__Bifidobacterium|s__Bifidobacterium_animalis", "k__Bacteria|g__Lacticaseibacillus|s__Lacticaseibacillus_rhamnosus"),
    ("d__Bacteria;g__Bifidobacterium;s__Bifidobacterium MOTU40_005954", "d__Bacteria;g__Lacticaseibacillus;s__Lacticaseibacillus rhamnosus"),
    ("d__Bacteria;g__Bifidobacterium;s__Unknown Bifidobacterium mOTUv4.0_005954",),
    ("d__Bacteria;g__Lacticaseibacillus;s__Lacticaseibacillus rhamnosus", "d__Bacteria;g__Blautia;s__MGYG000001580"),
])
def test_a_lane_that_found_the_test_organisms_passes(tmp_path: Path, monkeypatch, lineages: tuple[str, ...]):
    _lane_status(monkeypatch, {"state": "ran", "n_species": len(lineages)}, lineages)
    outcome = _run_lane(tmp_path)
    assert outcome.ok and outcome.reason == "" and outcome.found


def test_a_lane_that_ran_and_named_nothing_has_not_worked(tmp_path: Path, monkeypatch):
    _lane_status(monkeypatch, {"state": "ran", "n_species": 0})
    outcome = _run_lane(tmp_path)
    assert not outcome.ok and "did not find the test organisms" in outcome.reason


def test_a_lane_that_named_only_other_organisms_has_not_worked(tmp_path: Path, monkeypatch):
    """Reads from Bifidobacterium and Lacticaseibacillus, a lane answering Escherichia: its database is not its own."""
    _lane_status(monkeypatch, {"state": "ran", "n_species": 2},
                 ("d__Bacteria;g__Escherichia;s__Escherichia coli", "d__Bacteria;g__Blautia;s__Blautia wexlerae"))
    outcome = _run_lane(tmp_path)
    assert not outcome.ok and outcome.found == () and "Bifidobacterium" in outcome.reason


def test_a_failed_lane_gives_the_reason_a_person_can_act_on(tmp_path: Path, monkeypatch):
    _lane_status(monkeypatch, {
        "state": "failed", "reason": "RuntimeError: singlem pipe failed (exit 1): hmm - returned non-zero",
        "error_detail": "STDERR was: b'bash: line 1: hmmsearch: command not found\\n'",
    })
    outcome = _run_lane(tmp_path)
    assert not outcome.ok
    assert outcome.reason == "the program 'hmmsearch' is not installed - run ./install.sh to add it"


def test_a_lane_that_is_not_installed_says_so(tmp_path: Path, monkeypatch):
    _lane_status(monkeypatch, {"state": "not_assessed", "reason": "tool or reference not installed"})
    outcome = _run_lane(tmp_path)
    assert not outcome.ok and "not installed" in outcome.reason and "./install.sh" in outcome.reason


def test_a_lane_that_raises_is_reported_not_propagated(tmp_path: Path, monkeypatch):
    from openbiota.expansion import lanes

    def boom(**_kw):
        raise OSError("disk full")

    monkeypatch.setattr(lanes, "run_all", boom)
    outcome = _run_lane(tmp_path)
    assert not outcome.ok and "disk full" in outcome.reason


def test_each_lane_is_measured_in_a_process_of_its_own(tmp_path: Path):
    """The operating system reports the largest child a process ever had, so a lane's
    peak can only be read from a process that has run that lane and nothing else."""
    outcome = selftest.measure_lane("genome_profile", r1=tmp_path / "r1", r2=tmp_path / "r2",
                                    refs_dir=tmp_path, work=tmp_path, threads=1)
    # No sylph database in an empty folder: the lane reports that, from the child process.
    assert not outcome.ok and "./install.sh" in outcome.reason
    assert outcome.peak_gb is not None and outcome.peak_gb >= 0.0


def test_a_measuring_process_that_dies_is_a_failed_lane_not_a_hang(tmp_path: Path, monkeypatch):
    from concurrent.futures.process import BrokenProcessPool

    class Dead:
        def __init__(self, *_a, **_kw): pass
        def __enter__(self): return self
        def __exit__(self, *_a): return False
        def submit(self, *_a, **_kw): raise BrokenProcessPool("worker killed")

    import concurrent.futures

    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", Dead)
    outcome = selftest.measure_lane("motus4", r1=tmp_path / "r1", r2=tmp_path / "r2", refs_dir=tmp_path,
                                    work=tmp_path, threads=1)
    assert not outcome.ok and "ran out of memory" in outcome.reason


def test_the_whole_genome_lane_is_judged_by_its_confident_hits(tmp_path: Path, monkeypatch):
    from openbiota.engines import sylph

    class Hit:
        def __init__(self, lineage: str) -> None:
            self.lineage = lineage

    class Profile:
        confident_hits = (Hit("d__Bacteria;g__Bifidobacterium;s__Bifidobacterium longum"),
                          Hit("d__Bacteria;g__Lacticaseibacillus;s__Lacticaseibacillus rhamnosus"))

    monkeypatch.setattr(sylph, "available", lambda _db: True)
    monkeypatch.setattr(sylph, "run_sylph", lambda **_kw: Profile())
    assert _run_lane(tmp_path, "genome_profile").found == ("Bifidobacterium", "Lacticaseibacillus")
    Profile.confident_hits = ()
    assert not _run_lane(tmp_path, "genome_profile").ok
    monkeypatch.setattr(sylph, "available", lambda _db: False)
    missing = _run_lane(tmp_path, "genome_profile")
    assert not missing.ok and "./install.sh" in missing.reason


# --------------------------------------------------------------------------- #
# the command
# --------------------------------------------------------------------------- #


def _stub(monkeypatch, failing: dict[str, str] | None = None) -> list[str]:
    failing = failing or {}
    ran: list[str] = []
    monkeypatch.setattr(selftest, "simulate_reads", lambda work, **_kw: (work / "r1", work / "r2"))

    def measure_lane(lane, **_kw):
        ran.append(lane)
        if lane in failing:
            return LaneOutcome(lane, selftest.LANES[lane], False, 1.0, failing[lane])
        return LaneOutcome(lane, selftest.LANES[lane], True, 65.0, "", selftest.EXPECTED_GENERA, 33.6)

    monkeypatch.setattr(selftest, "measure_lane", measure_lane)
    return ran


def test_every_lane_passing_exits_zero_with_a_line_each(tmp_path: Path, monkeypatch):
    ran = _stub(monkeypatch)
    out = io.StringIO()
    assert selftest.run(out=out, keep=tmp_path / "keep") == 0
    assert ran == list(selftest.LANES), "every catalogue lane and the whole-genome profile"
    lines = out.getvalue().splitlines()
    assert lines[0] == ("Self-test: 6 detection lanes on 120,000 made-up read pairs from 3 gut genomes that ship "
                        "with the repository.")
    assert lines[1] == "A lane passes when it runs to the end and finds the organisms the reads were drawn from."
    assert sum(line.startswith("  ✓ ") for line in lines) == 6
    # Into a log, each lane is announced when it starts: the longest takes six minutes.
    assert lines[2] == "  · MetaPhlAn 4 (Jan26 catalogue) …"
    assert lines[3] == "  ✓ MetaPhlAn 4 (Jan26 catalogue)        1m05s    34 GB of memory"
    assert lines[-1].startswith("Every lane works (")


def test_on_a_terminal_the_running_lane_is_replaced_by_its_result(monkeypatch):
    _stub(monkeypatch)

    class Terminal(io.StringIO):
        def isatty(self) -> bool:
            return True

    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setenv("TERM", "xterm")
    out = Terminal()
    assert selftest.run(lanes=["motus4"], out=out) == 0
    shown = out.getvalue()
    assert "  · mOTUs 4 …\r\033[2K" in shown
    assert shown.count("\n") == 5, "two header lines, one result, a blank, the verdict"


def test_one_failing_lane_exits_one_and_says_which_and_why(tmp_path: Path, monkeypatch):
    why = "the program 'hmmsearch' is not installed - run ./install.sh to add it"
    _stub(monkeypatch, failing={"singlem_globdb": why})
    out = io.StringIO()
    assert selftest.run(out=out, keep=tmp_path / "keep") == 1
    text = out.getvalue()
    assert re.search(rf"  ✗ SingleM \(GlobDB\) +{re.escape(why)}", text)
    assert "1 of 6 lanes does not work on this machine" in text and "would be missing it." in text


def test_only_the_lanes_asked_for_are_run(tmp_path: Path, monkeypatch):
    ran = _stub(monkeypatch)
    assert selftest.run(lanes=["singlem_globdb", "motus4"], out=io.StringIO(), keep=tmp_path / "keep") == 0
    assert ran == ["singlem_globdb", "motus4"]


def test_an_unknown_lane_is_refused_with_the_known_ones(monkeypatch, capsys):
    ran = _stub(monkeypatch)
    assert selftest.run(lanes=["nonsense"], out=io.StringIO()) == 2
    assert ran == []
    assert "singlem_globdb" in capsys.readouterr().err


def test_reads_that_cannot_be_made_fail_the_test_plainly(tmp_path: Path, monkeypatch):
    def no_genomes(_work, **_kw):
        raise FileNotFoundError("the self-test genomes are not in refs/; run ./install.sh")

    monkeypatch.setattr(selftest, "simulate_reads", no_genomes)
    out = io.StringIO()
    assert selftest.run(out=out, keep=tmp_path / "keep") == 1
    assert "could not make the test reads" in out.getvalue() and "./install.sh" in out.getvalue()


def test_the_work_folder_is_removed_unless_kept(tmp_path: Path, monkeypatch):
    made: list[Path] = []
    _stub(monkeypatch)
    real = selftest.tempfile.mkdtemp

    def mkdtemp(**kw):
        made.append(Path(real(dir=tmp_path, **kw)))
        return str(made[-1])

    monkeypatch.setattr(selftest.tempfile, "mkdtemp", mkdtemp)
    assert selftest.run(out=io.StringIO()) == 0
    assert made and not made[0].exists()
    keep = tmp_path / "keep"
    assert selftest.run(out=io.StringIO(), keep=keep) == 0
    assert keep.is_dir()


def test_it_is_a_command(tmp_path: Path, monkeypatch):
    assert "selftest" in cli.SUBCOMMANDS
    called: dict = {}
    monkeypatch.setattr(selftest, "run", lambda **kw: called.update(kw) or 0)
    assert cli.main(["selftest", "--lanes", "motus4", "--threads", "3", "--keep", str(tmp_path)]) == 0
    assert called["lanes"] == ["motus4"] and called["threads"] == 3 and called["keep"] == tmp_path


# --------------------------------------------------------------------------- #
# the installer runs it
# --------------------------------------------------------------------------- #


BASH = shutil.which("bash")
needs_bash = pytest.mark.skipif(BASH is None, reason="bash is required to run the installer")


def _installer_with_a_fake_program(tmp_path: Path, script: str) -> Path:
    """The installer and launcher in a throwaway folder, with a stand-in for the program."""
    root = tmp_path / "repo"
    (root / "bin").mkdir(parents=True)
    shutil.copy(INSTALL, root / "install.sh")
    shutil.copy(REPO / "bin" / "openbiota", root / "bin" / "openbiota")
    fake = root / ".venv" / "bin" / "openbiota"
    fake.parent.mkdir(parents=True)
    fake.write_text("#!/usr/bin/env bash\n" + script)
    fake.chmod(0o755)
    return root


def _install(root: Path, *args: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "NO_COLOR": "1", "TERM": "dumb", "HOME": str(root.parent / "home")}
    return subprocess.run([BASH, "install.sh", *args], cwd=root, env=env, capture_output=True, text=True,
                          timeout=120, check=False)


def test_the_installer_ends_with_the_trial_run():
    text = INSTALL.read_text()
    steps = re.search(r'^STANDARD_STEPS="([^"]+)"', text, re.M).group(1).split()
    assert steps[-1] == "selftest"
    assert '"${STANDARD_STEPS% check selftest} $FULL_STEPS check selftest"' in text, "last with --full too"


@needs_bash
def test_the_installer_runs_the_program_itself_and_reports_success(tmp_path: Path):
    root = _installer_with_a_fake_program(tmp_path, 'echo "ran: $*" > "$(dirname "$0")/../../called"\nexit 0\n')
    proc = _install(root, "--only", "selftest", "--threads", "3")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (root / "called").read_text().strip() == "ran: selftest --threads 3"
    assert "every detection lane ran and found what was in the test reads" in proc.stdout
    assert "Trial run of every detection lane" in proc.stdout


@needs_bash
def test_a_lane_that_cannot_run_stops_the_install_and_is_named(tmp_path: Path):
    why = "the program 'hmmsearch' is not installed - run ./install.sh to add it"
    root = _installer_with_a_fake_program(tmp_path, f"""
echo "Self-test: 6 detection lanes"
echo "  ✓ MetaPhlAn 4 (Jan26 catalogue)        6m19s"
echo "  ✗ SingleM (GlobDB)                     {why}"
echo "1 of 6 lanes do not work on this machine (10m31s). A report made now would be missing them."
exit 1
""")
    proc = _install(root, "--only", "selftest")
    assert proc.returncode == 1
    assert "a detection lane does not work on this machine" in proc.stdout
    assert f"SingleM (GlobDB)                     {why}" in proc.stdout, "the failed lane and its reason, from the log"
    assert "MetaPhlAn 4" not in proc.stdout, "not the lanes that passed"
    assert "bin/openbiota selftest" in proc.stdout
    assert "Installed" not in proc.stdout


@needs_bash
def test_the_trial_run_is_skipped_with_the_lanes_and_never_run_dry(tmp_path: Path):
    root = _installer_with_a_fake_program(tmp_path, 'touch "$(dirname "$0")/../../called"\nexit 1\n')
    skipped = _install(root, "--only", "selftest", "--skip-catalogues")
    assert skipped.returncode == 0 and "the catalogue lanes are not installed" in skipped.stdout
    dry = _install(root, "--only", "selftest", "--dry-run")
    assert dry.returncode == 0 and "would run" in dry.stdout
    assert not (root / "called").exists()
