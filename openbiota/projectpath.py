"""The project's own tool folders, on the PATH whichever way the program starts.

The installer puts some tools inside the repository: micromamba's when the
system package manager had none or too old a one (``vendor/conda``), the
Rust builds (``vendor/cargo``, ``vendor/sylph``), single binaries
(``vendor/bin``). ``bin/openbiota`` and the Makefile put those folders on
the PATH; this does the same from inside Python, so ``.venv/bin/openbiota``,
``python -m openbiota`` and a test all find the same tools the launcher
does, in the same order.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Final

#: Relative to the repository, in the order they are searched. The launcher
#: (bin/openbiota) and the Makefile carry the same list; a test compares them.
TOOL_DIRS: Final = ("vendor/bin", "vendor/conda/bin", "vendor/cargo/bin", "vendor/sylph/bin")

REPO: Final = Path(__file__).resolve().parents[1]


def ensure(repo: Path | None = None) -> list[str]:
    """Put the folders that exist at the front of PATH; returns what was added."""
    root = repo or REPO
    current = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    add = [str(root / d) for d in TOOL_DIRS if (root / d).is_dir() and str(root / d) not in current]
    if add:
        os.environ["PATH"] = os.pathsep.join(add + current)
    return add
