"""The installer, the launcher and the readiness checks.

Nothing here downloads anything or touches the real environment: the
installer is exercised in --dry-run / --list / --help modes and on throwaway
copies of the repository skeleton under tmp_path, and the launcher on a
skeleton with and without a .venv. What is tested is the contract a fresh
machine depends on: every advertised step exists, every option parses, bad
input fails with a clear message, the final instructions name the real
commands, and the script passes shellcheck.
"""

from __future__ import annotations

import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from openbiota import installcheck

REPO = Path(__file__).resolve().parents[1]
INSTALL = REPO / "install.sh"
LAUNCHER = REPO / "bin" / "openbiota"
BASH = shutil.which("bash")

pytestmark = pytest.mark.skipif(BASH is None, reason="bash is required to test the installer")


def _sh(*args: str, cwd: Path | None = None, env: dict | None = None, stdin: str = "") -> subprocess.CompletedProcess:
    full_env = {**os.environ, "NO_COLOR": "1", "TERM": "dumb"}
    if env:
        full_env.update(env)
    return subprocess.run(
        [BASH, *args], cwd=cwd, env=full_env, input=stdin,
        capture_output=True, text=True, timeout=120, check=False,
    )


def _steps(var: str) -> list[str]:
    text = INSTALL.read_text()
    match = re.search(rf'^{var}="([^"]+)"', text, re.M)
    assert match, f"{var} not found in install.sh"
    return match.group(1).split()


def _skeleton(tmp_path: Path, *, with_venv: bool = False) -> Path:
    """A fake repository: the installer, the launcher, the data folder."""
    root = tmp_path / "repo"
    (root / "bin").mkdir(parents=True)
    (root / "scripts").mkdir()
    (root / "data").mkdir()
    shutil.copy(INSTALL, root / "install.sh")
    shutil.copy(LAUNCHER, root / "bin" / "openbiota")
    shutil.copy(REPO / "scripts" / "derived_refs.py", root / "scripts" / "derived_refs.py")
    if with_venv:
        venv_bin = root / ".venv" / "bin"
        venv_bin.mkdir(parents=True)
        fake = venv_bin / "openbiota"
        fake.write_text('#!/usr/bin/env bash\necho "fake openbiota: $*"\necho "PATH=$PATH"\n')
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    return root


# --------------------------------------------------------------------------- #
# static properties of the script
# --------------------------------------------------------------------------- #


def test_scripts_are_executable_and_parse():
    for script in (INSTALL, LAUNCHER):
        assert os.access(script, os.X_OK), f"{script.name} must be executable"
        proc = _sh("-n", str(script))
        assert proc.returncode == 0, proc.stderr


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck not installed")
def test_shellcheck_is_clean():
    proc = subprocess.run(
        ["shellcheck", "-s", "bash", str(INSTALL), str(LAUNCHER)],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_every_step_has_a_done_and_a_run_function():
    text = INSTALL.read_text()
    for step in _steps("STANDARD_STEPS") + _steps("FULL_STEPS"):
        assert re.search(rf"^done_{step}\(\)", text, re.M), f"done_{step} missing"
        assert re.search(rf"^run_{step}\(\)", text, re.M), f"run_{step} missing"


def test_every_step_has_a_label_and_is_in_the_help():
    text = INSTALL.read_text()
    usage = text[text.index("usage()"):text.index("while [ $# -gt 0 ]")]
    for step in _steps("STANDARD_STEPS") + _steps("FULL_STEPS"):
        assert re.search(rf"^\s+{step}\) label=", text, re.M), f"no label for {step}"
        assert step in usage, f"{step} is not described in --help"


def test_standard_steps_begin_with_preflight_and_end_with_the_check_and_the_trial_run():
    steps = _steps("STANDARD_STEPS")
    assert steps[0] == "preflight"
    assert steps[-2:] == ["check", "selftest"], "in place first, then proven to run"
    # Readiness is asked of the pipeline, not guessed by the shell.
    for name in installcheck.CHECKS:
        assert name in steps or name == "python", name


def test_readiness_checks_cover_the_installed_pieces():
    assert set(installcheck.CHECKS) == {
        "python", "profiler", "data", "metaphlan_db", "diamond_db", "models", "host_index",
        "catalogue_tools", "catalogue_refs",
    }
    # Every standard step that produces something has a check with its name.
    steps = set(_steps("STANDARD_STEPS")) - {"preflight", "system", "launcher", "check", "selftest"}
    assert steps == set(installcheck.CHECKS)
    # A skipped host index is a supported choice and not a missing piece.
    assert "host_index" not in installcheck.standard_missing(Path("/nonexistent"))
    assert "catalogue_refs" in installcheck.standard_missing(Path("/nonexistent"))


def test_full_install_checks_are_exactly_the_installers_full_steps():
    """One list: the installer's --full steps and the checks `openbiota run
    --mode full` consults are the same names, so neither can drift."""
    assert set(installcheck.FULL_CHECKS) == set(_steps("FULL_STEPS"))
    text = INSTALL.read_text()
    for step in _steps("FULL_STEPS"):
        assert re.search(rf"^done_{step}\(\)\s+\{{ ready {step}; \}}", text, re.M), (
            f"done_{step} must ask openbiota.installcheck, not re-derive the answer")


def test_every_full_check_answers_false_on_an_empty_folder(tmp_path: Path, monkeypatch):
    # The tool checks look at the PATH too; an empty PATH makes the folder the whole world.
    monkeypatch.setenv("PATH", str(tmp_path))
    for name, check in installcheck.FULL_CHECKS.items():
        assert check(tmp_path) is False, name
    assert installcheck.full_missing(tmp_path) == list(installcheck.FULL_CHECKS)
    assert installcheck.catalogue_tools(tmp_path) is False
    assert installcheck.catalogue_refs(tmp_path) is False


def test_substrate_check_rejects_a_web_page_standing_in_for_the_file(tmp_path: Path):
    d = tmp_path / "refs" / "dbcan"
    d.mkdir(parents=True)
    (d / "snapshot.json").write_text("{}")
    (d / "fam-substrate-mapping.tsv").write_text("<!DOCTYPE html>...")
    assert installcheck.substrate_refs(tmp_path) is False
    (d / "fam-substrate-mapping.tsv").write_text("#comment\nSubstrate_high_level\tFamily\n")
    assert installcheck.substrate_refs(tmp_path) is True


def test_the_installer_is_written_for_the_bash_macos_ships():
    """macOS has bash 3.2: no associative arrays, mapfile or ${var,,}."""
    code = "\n".join(
        line for line in INSTALL.read_text().splitlines() if not line.lstrip().startswith("#")
    )
    for forbidden in ("declare -A", "mapfile", "readarray", "|&"):
        assert forbidden not in code, f"{forbidden!r} is not available in bash 3.2"
    assert not re.search(r"\$\{[A-Za-z_]+(,,|\^\^)\}", code), "case-modifying expansions need bash 4"


# --------------------------------------------------------------------------- #
# command-line behaviour
# --------------------------------------------------------------------------- #


def test_help_exits_zero_and_names_the_options(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--help", cwd=root)
    assert proc.returncode == 0
    for option in ("--full", "--skip-host-index", "--only", "--skip", "--threads",
                   "--force", "--yes", "--offline", "--dry-run", "--list", "--no-path"):
        assert option in proc.stdout, option


def test_unknown_option_is_refused(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--bogus", cwd=root)
    assert proc.returncode == 2
    assert "unknown option" in proc.stderr


@pytest.mark.parametrize("args", [("--threads", "lots"), ("--threads",), ("--only",), ("--skip",)])
def test_malformed_options_are_refused(tmp_path: Path, args: tuple[str, ...]):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", *args, cwd=root)
    assert proc.returncode == 2, proc.stdout + proc.stderr


def test_unknown_step_is_refused(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--only", "teleport", "--dry-run", cwd=root)
    assert proc.returncode == 1
    assert "unknown step: teleport" in proc.stdout + proc.stderr


def test_list_shows_every_standard_step_in_order(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--list", cwd=root)
    assert proc.returncode == 0, proc.stderr
    listed = [line.split()[0] for line in proc.stdout.splitlines() if line.startswith("  ")]
    assert listed == _steps("STANDARD_STEPS")
    # On a bare skeleton nothing that lives in the repository can be done.
    for step in ("python", "profiler", "data", "metaphlan_db", "diamond_db"):
        assert re.search(rf"^  {step}\s+to do$", proc.stdout, re.M), step


def test_list_with_full_appends_the_optional_lanes(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--list", "--full", cwd=root)
    listed = [line.split()[0] for line in proc.stdout.splitlines() if line.startswith("  ")]
    standard = _steps("STANDARD_STEPS")
    # The final check and the trial run move to the end, after the lanes they must also cover.
    assert listed == standard[:-2] + _steps("FULL_STEPS") + ["check", "selftest"]


def test_dry_run_launcher_step_changes_nothing(tmp_path: Path):
    root = _skeleton(tmp_path)
    home = tmp_path / "home"
    home.mkdir()
    proc = _sh("install.sh", "--dry-run", "--only", "launcher", cwd=root, env={"HOME": str(home)})
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "would link" in proc.stdout
    assert not (home / ".local").exists()
    assert "Installed" in proc.stdout


def test_launcher_step_links_into_local_bin_and_marks_the_rc(tmp_path: Path):
    root = _skeleton(tmp_path)
    home = tmp_path / "home"
    home.mkdir()
    env = {"HOME": str(home), "SHELL": "/bin/zsh", "PATH": os.environ["PATH"]}
    proc = _sh("install.sh", "--only", "launcher", cwd=root, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    link = home / ".local" / "bin" / "openbiota"
    assert link.is_symlink()
    assert link.readlink() == root / "bin" / "openbiota"
    rc = (home / ".zshrc").read_text()
    assert "# added by openbiota install.sh" in rc
    assert str(home / ".local" / "bin") in rc
    assert "new terminal" in proc.stdout
    # Running again is a no-op: the rc file is not appended to twice.
    again = _sh("install.sh", "--only", "launcher", cwd=root, env=env)
    assert again.returncode == 0
    assert "already done" in again.stdout
    assert rc == (home / ".zshrc").read_text()


def test_launcher_step_leaves_a_foreign_openbiota_alone(tmp_path: Path):
    root = _skeleton(tmp_path)
    home = tmp_path / "home"
    (home / ".local" / "bin").mkdir(parents=True)
    foreign = home / ".local" / "bin" / "openbiota"
    foreign.write_text("#!/bin/sh\necho someone else\n")
    env = {"HOME": str(home), "SHELL": "/bin/bash"}
    proc = _sh("install.sh", "--only", "launcher", cwd=root, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "not a link" in proc.stdout
    assert foreign.read_text().startswith("#!/bin/sh")


def test_no_path_skips_the_rc_file(tmp_path: Path):
    root = _skeleton(tmp_path)
    home = tmp_path / "home"
    home.mkdir()
    proc = _sh("install.sh", "--only", "launcher", "--no-path", cwd=root, env={"HOME": str(home)})
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not (home / ".local").exists()
    assert not list(home.glob(".*rc")) and not (home / ".bash_profile").exists()


def test_final_message_names_the_two_remaining_steps(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--only", "launcher", "--no-path", cwd=root,
               env={"HOME": str(tmp_path)})
    out = proc.stdout
    assert "fastq/my-sample/my-sample_1.fastq.gz" in out
    assert "bin/openbiota run" in out
    assert "results/my-sample/my-sample_report.pdf" in out


@pytest.mark.parametrize("fs,kind", [
    ("cifs", "network share"), ("smbfs", "network share"), ("nfs4", "network share"),
    ("sshfs", "network share"), ("exfat", "formatted for Windows"), ("ntfs", "formatted for Windows"),
    ("vfat", "formatted for Windows"), ("fuseblk", "formatted for Windows"),
])
def test_preflight_refuses_a_share_or_a_windows_drive(tmp_path: Path, fs: str, kind: str):
    """A Linux box cloning onto a Mac's SMB share was the first real failure
    a user hit: git itself fails there, and so would every later step."""
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--only", "preflight", "--dry-run", "--no-path", cwd=root,
               env={"INSTALL_FS_TYPE": fs})
    out = proc.stdout + proc.stderr
    assert proc.returncode == 1, out
    assert kind in out and f"({fs})" in out
    assert "cd ~ && git clone" in out, "the fix must be spelled out"
    assert not list(root.glob(".install-probe.*")), "the probe folder must not be left behind"


def test_preflight_accepts_a_local_filesystem(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--only", "preflight", "--dry-run", "--no-path", cwd=root,
               env={"INSTALL_FS_TYPE": "ext4"})
    assert "local filesystem (ext4)" in proc.stdout, proc.stdout + proc.stderr
    assert not list(root.glob(".install-probe.*"))


def test_preflight_only_warns_on_an_unknown_fuse_filesystem(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("install.sh", "--only", "preflight", "--dry-run", "--no-path", cwd=root,
               env={"INSTALL_FS_TYPE": "fuse.snapfuse"})
    assert "FUSE filesystem" in proc.stdout
    assert "local filesystem (fuse.snapfuse)" in proc.stdout


def _preflight(root: Path, memory_gb: int, *more: str) -> str:
    proc = _sh("install.sh", "--only", "preflight", "--dry-run", "--no-path", *more, cwd=root,
               env={"INSTALL_FS_TYPE": "ext4", "INSTALL_MEMORY_GB": str(memory_gb)})
    return proc.stdout + proc.stderr


def test_preflight_says_when_memory_is_short_for_the_metaphlan_lane(tmp_path: Path):
    """The Jan26 index is held in memory, 37 GB in one program (measured by the
    self-test). A machine without that has the lane killed on every run; it is
    told before the 130 GB download, not after."""
    root = _skeleton(tmp_path)
    text = INSTALL.read_text()
    assert re.search(r"^LANE_PEAK_GB=37$", text, re.M)
    small = _preflight(root, 31)
    assert "31 GB of memory: the MetaPhlAn 4 lane holds a 37 GB index in memory" in small
    assert "every report will be missing it" in small and "64 GB is comfortable" in small
    assert "MetaPhlAn 4 lane" not in _preflight(root, 62), "a 64 GB machine is not warned"
    assert "MetaPhlAn 4 lane" not in _preflight(root, 31, "--skip-catalogues"), "nor one that is not installing the lanes"


def test_the_filesystem_check_runs_before_anything_is_downloaded():
    text = INSTALL.read_text()
    preflight = text[text.index("run_preflight()"):text.index("# -- system")]
    assert preflight.index("check_filesystem") < preflight.index("free_gb")
    assert preflight.index("check_filesystem") < preflight.index("online")


def test_data_step_fails_clearly_without_the_archive(tmp_path: Path):
    root = _skeleton(tmp_path, with_venv=True)
    (root / ".venv" / "bin" / "python").write_text("#!/usr/bin/env bash\nexit 1\n")
    (root / ".venv" / "bin" / "python").chmod(0o755)
    proc = _sh("install.sh", "--only", "data", "--no-path", cwd=root, env={"HOME": str(tmp_path)})
    assert proc.returncode == 1
    assert "derived-refs.manifest.json is missing" in proc.stdout + proc.stderr
    assert "re-clone" in proc.stdout + proc.stderr


def test_a_second_installer_is_refused_while_one_runs(tmp_path: Path):
    root = _skeleton(tmp_path)
    (root / ".install.lock").write_text(str(os.getpid()))  # this test process is "running"
    proc = _sh("install.sh", "--only", "launcher", "--dry-run", cwd=root)
    assert proc.returncode == 1
    assert "another ./install.sh" in proc.stdout + proc.stderr


def test_a_stale_lock_is_cleared(tmp_path: Path):
    root = _skeleton(tmp_path)
    (root / ".install.lock").write_text("999999")
    proc = _sh("install.sh", "--only", "launcher", "--dry-run", cwd=root)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not (root / ".install.lock").exists()


def test_everything_the_installer_prints_is_also_in_the_log(tmp_path: Path):
    root = _skeleton(tmp_path)
    _sh("install.sh", "--only", "launcher", "--dry-run", cwd=root)
    log = (root / "install.log").read_text()
    assert "install.sh started" in log
    assert "step 1/1" in log


# --------------------------------------------------------------------------- #
# the launcher
# --------------------------------------------------------------------------- #


def test_launcher_without_a_venv_points_at_the_installer(tmp_path: Path):
    root = _skeleton(tmp_path)
    proc = _sh("bin/openbiota", "run", cwd=root)
    assert proc.returncode == 1
    assert "./install.sh" in proc.stderr
    assert str(root) in proc.stderr


def test_launcher_runs_the_venv_program_with_vendor_on_path(tmp_path: Path):
    root = _skeleton(tmp_path, with_venv=True)
    proc = _sh("bin/openbiota", "doctor", "--sample", "x", cwd=root)
    assert proc.returncode == 0, proc.stderr
    assert "fake openbiota: doctor --sample x" in proc.stdout
    assert proc.stdout.splitlines()[1].startswith(f"PATH={root}/vendor/bin:{root}/vendor/conda/bin:")


def test_launcher_passes_no_arguments_through_unchanged(tmp_path: Path):
    """A bare `openbiota` reaches the program with no arguments, which prints
    the help; the launcher no longer turns it into a run."""
    root = _skeleton(tmp_path, with_venv=True)
    proc = _sh("bin/openbiota", cwd=root)
    assert proc.stdout.splitlines()[0] == "fake openbiota: "


def test_launcher_works_through_a_symlink_elsewhere(tmp_path: Path):
    root = _skeleton(tmp_path, with_venv=True)
    elsewhere = tmp_path / "home" / ".local" / "bin"
    elsewhere.mkdir(parents=True)
    (elsewhere / "openbiota").symlink_to(root / "bin" / "openbiota")
    proc = _sh(str(elsewhere / "openbiota"), "--version", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert "fake openbiota: --version" in proc.stdout


# --------------------------------------------------------------------------- #
# the readiness checks
# --------------------------------------------------------------------------- #


def test_every_check_answers_false_on_an_empty_folder_without_raising(tmp_path: Path):
    for name, check in installcheck.CHECKS.items():
        assert check(tmp_path) is False, name


def test_check_cli_reports_and_exits_nonzero_when_missing(tmp_path: Path, capsys):
    assert installcheck.main(["data", "--repo", str(tmp_path)]) == 1
    assert installcheck.main(["--all", "--repo", str(tmp_path)]) == 1
    out = capsys.readouterr().out
    assert "data" in out and "MISSING" in out


def test_data_check_needs_every_restored_piece(tmp_path: Path):
    refs = tmp_path / "refs"
    for p in ("age", "crosswalk", "panels"):
        (refs / p).mkdir(parents=True)
    (refs / "reference_ranges.json").write_text("{}")
    assert installcheck.data(tmp_path) is False
    (refs / "taxonomic_cohort.json").write_text("{}")
    assert installcheck.data(tmp_path) is True


# --------------------------------------------------------------------------- #
# make install
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(shutil.which("make") is None, reason="make not installed")
def test_make_install_runs_the_script_with_its_arguments():
    proc = subprocess.run(
        ["make", "-n", "install", "INSTALL_ARGS=--skip-host-index"],
        cwd=REPO, capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "./install.sh --skip-host-index" in proc.stdout
    proc = subprocess.run(["make", "-n", "install-full"], cwd=REPO, capture_output=True, text=True, check=False)
    assert "./install.sh --full" in proc.stdout


def test_every_path_the_project_sets_up_agrees():
    """The installer, the Makefile and the launcher put the same vendored
    folders first on the PATH, so a tool found during the install is the
    tool found at run time: vendor/bin, the conda fallback, the cargo-built
    helpers (smafa, mfqe) and sylph."""
    for text, name in ((INSTALL.read_text(), "install.sh"), ((REPO / "Makefile").read_text(), "Makefile"),
                       (LAUNCHER.read_text(), "bin/openbiota")):
        path_line = next(line for line in text.splitlines() if line.lstrip().startswith("export PATH"))
        # install.sh spells the first two through its own variables.
        spelled = path_line.replace("$VENDOR_BIN", "vendor/bin").replace("$CONDA/bin", "vendor/conda/bin")
        for folder in ("vendor/bin", "vendor/conda/bin", "vendor/cargo/bin", "vendor/sylph/bin"):
            assert folder in spelled, f"{name} does not put {folder} on the PATH: {path_line}"
    mk = (REPO / "Makefile").read_text()
    assert "--root vendor/cargo smafa" in mk and "--root vendor/cargo mfqe" in mk, "cargo tools must be built into the project"
    assert "--root vendor/sylph sylph" in mk


def test_rust_from_a_distribution_is_treated_as_absent_when_too_old():
    text = INSTALL.read_text()
    assert "RUST_MIN_MINOR=85" in text
    assert re.search(r"^cargo_ok\(\)", text, re.M)
    assert "cc cmake pkg-config" in text, "the Rust builds need a compiler and cmake first"
