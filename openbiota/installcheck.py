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


#: The standard install: what a lite report needs.
CHECKS: dict[str, Check] = {
    "python": python_env,
    "profiler": profiler,
    "data": data,
    "metaphlan_db": metaphlan_db,
    "diamond_db": diamond_db,
    "models": models,
    "host_index": host_index,
}


# --------------------------------------------------------------------------- #
# The full install: the optional detection lanes. Each check is the question
# the installer's matching step asks before skipping itself, written once
# here so the installer, `openbiota run --mode full` and `openbiota doctor`
# cannot disagree about what "installed" means.
# --------------------------------------------------------------------------- #


def profiler4(repo: Path) -> bool:
    return (repo / ".venv-mpa4" / "bin" / "metaphlan").exists()


def vendor_tools(repo: Path) -> bool:
    return (repo / "vendor" / "bin" / "blastn").exists() or any(
        repo.glob("vendor/ncbi-blast-*/bin/blastn"))


def strain_tools(repo: Path) -> bool:
    return (repo / ".venv-strain" / "bin" / "ectyper").exists()


def refs_strain(repo: Path) -> bool:
    return (repo / "refs" / "strain_refs.lock.json").is_file() and (
        repo / "refs" / "panels" / "probiotics").is_dir()


def pathogens(repo: Path) -> bool:
    return (repo / "refs" / "pathogens" / "bundle" / "manifest.json").is_file()


def substrate_refs(repo: Path) -> bool:
    """The dbCAN mapping is present and is a table, not the web page the
    old host now serves in place of every file."""
    mapping = repo / "refs" / "dbcan" / "fam-substrate-mapping.tsv"
    if not mapping.is_file() or not (repo / "refs" / "dbcan" / "snapshot.json").is_file():
        return False
    try:
        with mapping.open("rb") as handle:
            return handle.read(1) not in (b"<", b"")
    except OSError:
        return False


def tools_expanded(repo: Path) -> bool:
    return all((repo / env / "bin" / tool).exists() for env, tool in (
        (".venv-mpa42", "metaphlan"), (".venv-motus", "motus"), (".venv-singlem", "singlem")))


def refs_expanded(repo: Path) -> bool:
    refs = repo / "refs"
    return (refs / "motus" / "db").is_dir() and (refs / "kraken2").is_dir() and (
        refs / "metaphlan4_db").is_dir()


def expansion(repo: Path) -> bool:
    return (repo / "refs" / "expanded" / "reconciliation" / "summary.json").is_file()


def genome_lane(repo: Path) -> bool:
    return (repo / "vendor" / "sylph" / "bin" / "sylph").exists() and any(
        (repo / "refs" / "sylph").glob("*.syl2db"))


def mycobiome_refs(repo: Path) -> bool:
    return (repo / "refs" / "mycobiome" / "reference_lock.json").is_file()


FULL_CHECKS: dict[str, Check] = {
    "profiler4": profiler4,
    "vendor_tools": vendor_tools,
    "strain_tools": strain_tools,
    "refs_strain": refs_strain,
    "pathogens": pathogens,
    "substrate_refs": substrate_refs,
    "tools_expanded": tools_expanded,
    "refs_expanded": refs_expanded,
    "expansion": expansion,
    "genome_lane": genome_lane,
    "mycobiome_refs": mycobiome_refs,
}

ALL_CHECKS: dict[str, Check] = {**CHECKS, **FULL_CHECKS}


def full_missing(repo: Path | None = None) -> list[str]:
    """The full-install pieces that are not present, in install order."""
    root = repo or REPO_ROOT
    return [name for name, check in FULL_CHECKS.items() if not check(root)]


def full_ready(repo: Path | None = None) -> bool:
    return not full_missing(repo)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", nargs="?", choices=sorted(ALL_CHECKS))
    ap.add_argument("--all", action="store_true", help="report every standard-install check")
    ap.add_argument("--full", action="store_true", help="report every full-install check")
    ap.add_argument("--repo", type=Path, default=REPO_ROOT)
    args = ap.parse_args(argv)
    if not (args.all or args.full) and args.name is None:
        ap.error("give a check name, --all or --full")
    if args.all or args.full:
        worst = 0
        table = {**(CHECKS if args.all else {}), **(FULL_CHECKS if args.full else {})}
        for name, check in table.items():
            ready = check(args.repo)
            print(f"  {name:<14}{'ready' if ready else 'MISSING'}")
            worst |= 0 if ready else 1
        return worst
    return 0 if ALL_CHECKS[args.name](args.repo) else 1


if __name__ == "__main__":
    sys.exit(main())
