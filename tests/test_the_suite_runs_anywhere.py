"""The suite gives the same answer wherever, and however many times at once, it is run.

Twice it did not, and neither time was the reason visible in the test that
failed:

* It depended on the terminal settings of whoever started it. Two tests
  compared printed output byte for byte and passed only where NO_COLOR was
  exported; GitHub's runner does not export it, and the suite went red there
  after passing on the machine it was written on. ``tests/conftest.py`` now
  clears those variables before every test.
* It depended on being the only run on the machine. Seven tests built a PDF
  at a fixed path in /tmp and read it back; two runs of the suite met in
  that file and one of them read it empty. They build in memory now.

And one way it could have, that was closed before it happened: the linter
was whichever ruff was newest on the day a machine was installed - and on
GitHub, newest at every run. ruff is pinned to one version now.
"""

from __future__ import annotations

import os
import re
import tomllib
from pathlib import Path

from tests import conftest

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent

#: A quoted path under the machine-wide temporary folder.
SHARED_PATH = re.compile(r"""["']/(?:var/)?tmp/""")
#: (file, text on the line): a path that is only a label in a fixture and is never opened.
LABELS = {("test_report.py", 'path="/tmp/x_1.fastq.gz"')}


def test_no_test_inherits_the_terminal_settings_of_whoever_runs_the_suite():
    for name in conftest.PRESENTATION_VARIABLES:
        assert name not in os.environ, name
    assert {"NO_COLOR", "FORCE_COLOR", "COLUMNS"} <= set(conftest.PRESENTATION_VARIABLES)
    assert os.environ["TERM"] == conftest.TEST_TERM


def test_no_test_writes_to_a_fixed_path_in_tmp():
    """A test's files go under ``tmp_path`` or stay in memory: a path every run shares is a
    path two runs can meet in - a second terminal, or a check that runs the suite under
    several environments at once."""
    found = []
    for path in sorted(TESTS.glob("*.py")):
        if path == Path(__file__).resolve():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if SHARED_PATH.search(line) and not any(path.name == name and text in line for name, text in LABELS):
                found.append(f"{path.name}:{number}: {line.strip()}")
    assert not found, "use tmp_path or io.BytesIO instead:\n" + "\n".join(found)


def test_the_check_above_would_see_one():
    assert SHARED_PATH.search('path = "/tmp/_render.pdf"')
    assert SHARED_PATH.search("report = Path('/var/tmp/report.tsv')")
    assert SHARED_PATH.search('path = f"/tmp/_layout_{fn}.pdf"')
    assert not SHARED_PATH.search("report = tmp_path / 'report.tsv'")
    assert not SHARED_PATH.search("# builds in memory, not under /tmp/")


# --------------------------------------------------------------------------- #
# the same tools everywhere
# --------------------------------------------------------------------------- #


def _dev_tools() -> dict[str, str]:
    """name -> requirement, from the `dev` extra everyone installs."""
    project = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    return {re.split(r"[<>=!~ ]", spec, maxsplit=1)[0].lower(): spec
            for spec in project["optional-dependencies"]["dev"]}


def test_the_linter_is_one_version_everywhere():
    """A ruff release can switch on more rules in the families pyproject.toml selects. Left
    floating, unchanged code fails the lint on GitHub, which installs the newest at every run,
    and passes on a machine installed a month earlier. Moving to a newer ruff is a line to
    change on purpose, not something that happens overnight."""
    spec = _dev_tools()["ruff"]
    assert re.fullmatch(r"ruff==\d+\.\d+\.\d+", spec), f"ruff is to be pinned to one version, not {spec!r}"


def test_the_test_tools_stop_before_their_next_major_version():
    tools = _dev_tools()
    for name in ("pytest", "pypdf", "pymupdf"):
        assert re.search(r"<\s*\d+", tools[name]), f"{name} has no upper bound: {tools[name]!r}"


def test_github_the_installer_and_make_all_install_those_tools():
    """A pin means something only if everyone installs through it."""
    workflow = (REPO / ".github" / "workflows" / "tests.yml").read_text(encoding="utf-8")
    assert "pip install -e '.[dev]'" in workflow
    assert not re.search(r"pip install[^\n]*\b(ruff|pytest|pypdf)\b", workflow), "a tool installed beside its pin"
    assert '-e ".[dev]"' in (REPO / "install.sh").read_text(encoding="utf-8")
    assert "-e '.[dev]'" in (REPO / "Makefile").read_text(encoding="utf-8")
