"""The install instructions describe the installer that actually exists.

A quickstart drifts: a flag is renamed, a step is added, a path changes, and
the page that a first-time user follows is the last thing anyone re-reads.
These tests pin the documentation to the script - every option the docs
mention is one the script accepts, every step they name is one it runs, the
three-step story is told in the right order, and the folder layout the docs
draw is the one the discovery code reads.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from openbiota.search import LAYOUT_HELP

REPO = Path(__file__).resolve().parents[1]
INSTALL = (REPO / "install.sh").read_text()
QUICKSTART = (REPO / "docs" / "QUICKSTART.md").read_text()
INSTALL_DOC = (REPO / "docs" / "INSTALL.md").read_text()
README = (REPO / "README.md").read_text()
FASTQ_README = (REPO / "fastq" / "README.md").read_text()
USAGE = (REPO / "docs" / "USAGE.md").read_text()
MAKEFILE = (REPO / "Makefile").read_text()


def _script_options() -> set[str]:
    """Every long option the argument loop accepts."""
    loop = INSTALL[INSTALL.index("while [ $# -gt 0 ]"):INSTALL.index("case \"${THREADS:-1}\"")]
    # `--only)`, `--only=*)` and `--yes|-y)` are all the same option.
    return set(re.findall(r"(--[a-z-]+)(?:=\*|\|-[a-z])?\)", loop))


def _script_steps() -> list[str]:
    steps: list[str] = []
    for var in ("STANDARD_STEPS", "FULL_STEPS"):
        steps += re.search(rf'^{var}="([^"]+)"', INSTALL, re.M).group(1).split()
    return steps


@pytest.mark.parametrize("doc_name,doc", [
    ("QUICKSTART", QUICKSTART), ("INSTALL", INSTALL_DOC), ("README", README), ("USAGE", USAGE),
], ids=["QUICKSTART", "INSTALL", "README", "USAGE"])
def test_every_installer_option_in_the_docs_exists(doc_name: str, doc: str):
    mentioned = set(re.findall(r"install\.sh (--[a-z-]+)", doc))
    mentioned |= set(re.findall(r"`(--[a-z-]+)`", doc)) & {
        "--full", "--skip-host-index", "--no-path", "--only", "--skip", "--threads",
        "--force", "--yes", "--offline", "--dry-run", "--list", "--no-color",
    }
    unknown = mentioned - _script_options()
    assert not unknown, f"{doc_name} mentions options install.sh does not have: {sorted(unknown)}"


def test_install_doc_describes_every_step_once():
    steps = _script_steps()
    for step in steps:
        assert f"`{step}`" in INSTALL_DOC, f"docs/INSTALL.md does not describe the {step} step"
    # The numbered table counts the standard steps, and so does the quickstart.
    n_standard = len(re.search(r'^STANDARD_STEPS="([^"]+)"', INSTALL, re.M).group(1).split())
    words = {11: "eleven", 10: "ten", 12: "twelve"}
    assert words[n_standard] in INSTALL_DOC
    assert words[n_standard] in QUICKSTART


def test_quickstart_tells_the_three_steps_in_order():
    order = [
        "open a terminal",
        "git clone https://github.com/openbiota-org/openbiota-gut-health-test",
        "cd openbiota-gut-health-test",
        "./install.sh",
        "fastq/my-sample/my-sample_1.fastq.gz",
        "bin/openbiota run",
        "results/my-sample/my-sample_report.pdf",
    ]
    # After the short version at the top, which is the same story in five lines.
    body = QUICKSTART[QUICKSTART.index("## Before you start"):].lower()
    positions = [body.find(item.lower()) for item in order]
    assert all(p >= 0 for p in positions), dict(zip(order, positions, strict=True))
    assert positions == sorted(positions), "the quickstart tells the steps out of order"


def test_quickstart_explains_the_terminal_before_asking_for_commands():
    first_command = QUICKSTART.index("## Step 1")
    explanation = QUICKSTART.index("## Before you start: open a terminal")
    assert explanation < first_command
    section = QUICKSTART[explanation:first_command]
    assert "git clone" in section and "downloads a copy" in section
    assert "images/terminal.svg" in section
    assert (REPO / "docs" / "images" / "terminal.svg").is_file()
    for platform in ("macOS", "Linux"):
        assert platform in section


def test_the_run_step_is_the_program_not_a_make_target():
    step3 = QUICKSTART[QUICKSTART.index("## Step 3"):QUICKSTART.index("## Optional")]
    assert "bin/openbiota run" in step3
    assert "make run" not in step3
    assert "bin/openbiota run" in README
    assert "make run" not in README.split("## Commands")[0]


def test_the_documented_layout_is_the_one_the_code_reads():
    # The same example everywhere: a folder named my-sample holding
    # my-sample_1 and my-sample_2, in the docs, the folder's own README and
    # the error message a user sees when the layout is wrong.
    for name, doc in (("QUICKSTART", QUICKSTART), ("INSTALL", INSTALL_DOC),
                      ("fastq/README", FASTQ_README), ("LAYOUT_HELP", LAYOUT_HELP)):
        assert "my-sample/" in doc, name
        assert "my-sample_1.fastq.gz" in doc and "my-sample_2.fastq.gz" in doc, name
    assert "fastq/my-sample/my-sample_1.fastq.gz" in LAYOUT_HELP


def test_installer_final_message_matches_the_quickstart():
    for text in ("fastq/my-sample/my-sample_1.fastq.gz", "bin/openbiota run",
                 "results/my-sample/my-sample_report.pdf"):
        assert text in INSTALL, f"install.sh's final message lost {text!r}"
        assert text in QUICKSTART


def test_make_install_exists_and_is_documented():
    assert re.search(r"^install:\n\t@\./install\.sh \$\(INSTALL_ARGS\)", MAKEFILE, re.M)
    assert re.search(r"^install-full:\n\t@\./install\.sh --full", MAKEFILE, re.M)
    assert "make install" in USAGE and "make install" in INSTALL_DOC
    assert "INSTALL_ARGS" in USAGE


def test_fastq_folder_is_ignored_but_its_readme_is_kept():
    gitignore = (REPO / ".gitignore").read_text().splitlines()
    assert "fastq/*" in gitignore
    assert "!fastq/README.md" in gitignore
    assert "install.log" in gitignore and ".install.lock" in gitignore
    assert (REPO / "fastq" / "README.md").is_file()


@pytest.mark.skipif(not (REPO / ".git").exists(), reason="not a git checkout")
def test_data_archive_is_tracked_not_ignored():
    import subprocess
    proc = subprocess.run(
        ["git", "check-ignore", "data/derived-refs.tar.gz", "data/derived-refs.manifest.json",
         "bin/openbiota", "install.sh", "fastq/README.md"],
        cwd=REPO, capture_output=True, text=True, check=False,
    )
    # check-ignore prints the ignored paths and exits 1 when there are none
    assert proc.returncode == 1 and proc.stdout == "", f"ignored: {proc.stdout}"


def test_disk_and_time_claims_agree_between_script_and_docs():
    standard = int(re.search(r"^DISK_STANDARD_GB=(\d+)", INSTALL, re.M).group(1))
    no_host = int(re.search(r"^DISK_NO_HOST_GB=(\d+)", INSTALL, re.M).group(1))
    assert f"{standard} GB" in QUICKSTART and f"~{standard} GB" in INSTALL_DOC
    assert f"{no_host} GB" in INSTALL_DOC


def test_ci_exercises_the_installer_on_both_platforms():
    workflow = (REPO / ".github" / "workflows" / "tests.yml").read_text()
    assert "ubuntu-latest" in workflow and "macos-latest" in workflow
    assert "./install.sh --only" in workflow
    assert "shellcheck -s bash install.sh bin/openbiota" in workflow
