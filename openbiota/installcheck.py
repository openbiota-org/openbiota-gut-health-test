"""Is each installed piece actually there? One question per name.

    python -m openbiota.installcheck diamond_db      # exit 0 when ready
    python -m openbiota.installcheck --all           # a table of everything

``install.sh`` asks these questions to decide what it can skip, and asks
them again at the end to prove the install is complete. They live here, not
in the shell script, so the same answer is available to ``openbiota doctor``
and to the tests, and so a check can use the pipeline's own notion of
"ready" - the DIAMOND database is keyed by a fingerprint of the panels and
the DIAMOND version, and only :mod:`openbiota.references` knows how to
compute it.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path

from openbiota.cli import REPO_ROOT

Check = Callable[[Path], bool]


def python_env(repo: Path) -> bool:
    return (repo / ".venv" / "bin" / "openbiota").exists()


def profiler(repo: Path) -> bool:
    return (repo / ".venv-mpa3" / "bin" / "metaphlan").exists()


def data(repo: Path) -> bool:
    """The archive in data/ has been unpacked into refs/."""
    refs = repo / "refs"
    needed = (
        refs / "reference_ranges.json",
        refs / "taxonomic_cohort.json",
        refs / "age",
        refs / "crosswalk",
        refs / "panels",
    )
    return all(p.exists() for p in needed)


def metaphlan_db(repo: Path) -> bool:
    from openbiota.engines.metaphlan import database_ready
    return database_ready(repo / "refs" / "metaphlan_db")


def diamond_db(repo: Path) -> bool:
    """Built for *this* machine's DIAMOND: the fingerprint includes its version."""
    from openbiota.cli import _entries_for_db
    from openbiota.errors import DependencyError
    from openbiota.panels import load_panel_set
    from openbiota.references import database_fingerprint, diamond_version
    try:
        version = diamond_version()
    except DependencyError:
        return False
    if not (repo / "panels").is_dir():
        return False
    panels = load_panel_set(repo / "panels")
    fingerprint = database_fingerprint(panels, _entries_for_db(panels, None), version)
    db_dir = repo / "refs" / "db" / fingerprint
    return (db_dir / "index.json").is_file() and (db_dir / "combined.dmnd").is_file()


def models(repo: Path) -> bool:
    try:
        import highspy  # noqa: F401
        import micom  # noqa: F401
    except ImportError:
        return False
    micom_dir = repo / "refs" / "micom"
    return (micom_dir / "agora201_refseq216_species_1.qza").is_file() and (
        micom_dir / "s10" / "S5_Data.xlsx").is_file()


def host_index(repo: Path) -> bool:
    from openbiota.engines.metaphlan import host_index_ready
    return host_index_ready(repo / "refs" / "host")


CHECKS: dict[str, Check] = {
    "python": python_env,
    "profiler": profiler,
    "data": data,
    "metaphlan_db": metaphlan_db,
    "diamond_db": diamond_db,
    "models": models,
    "host_index": host_index,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", nargs="?", choices=sorted(CHECKS))
    ap.add_argument("--all", action="store_true", help="report every check")
    ap.add_argument("--repo", type=Path, default=REPO_ROOT)
    args = ap.parse_args(argv)
    if not args.all and args.name is None:
        ap.error("give a check name or --all")
    if args.all:
        worst = 0
        for name, check in CHECKS.items():
            ready = check(args.repo)
            print(f"  {name:<14}{'ready' if ready else 'MISSING'}")
            worst |= 0 if ready else 1
        return worst
    return 0 if CHECKS[args.name](args.repo) else 1


if __name__ == "__main__":
    sys.exit(main())
