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


def catalogue_tools_missing(repo: Path) -> list[str]:
    """Which of the catalogue lanes' tools are absent, by the name a person would install.

    Each engine looks in its own place (a venv under the repository, then
    the PATH); the same places are checked here, so a tool that came from
    conda or the system counts.
    """
    import shutil  # noqa: PLC0415

    from openbiota.engines.singlem import missing_path_tools  # noqa: PLC0415

    def on_path(tool: str) -> bool:
        return shutil.which(tool) is not None

    def venv_or_path(venv: str, tool: str) -> bool:
        return (repo / venv / "bin" / tool).exists() or on_path(tool)

    def vendored_or_path(sub: str, tool: str) -> bool:
        return (repo / "vendor" / sub / "bin" / tool).exists() or on_path(tool)

    present = {
        "MetaPhlAn 4.2": (repo / ".venv-mpa42" / "bin" / "metaphlan").exists(),   # the Jan26 lane looks only here
        "mOTUs": venv_or_path(".venv-motus", "motus"),
        "SingleM": venv_or_path(".venv-singlem", "singlem"),
        "kraken2": on_path("kraken2"), "bracken": on_path("bracken"),
        "bwa": on_path("bwa"), "bowtie2": on_path("bowtie2"),
        "sylph": vendored_or_path("sylph", "sylph"),
        "sylph-tax": venv_or_path(".venv", "sylph-tax"),
        "SPAdes": any(repo.glob("vendor/SPAdes-*/bin/metaspades.py")) or on_path("metaspades.py"),
        # SingleM runs these through the PATH: two Rust helpers, and the
        # programs it pipes reads through (OrfM into HMMER's hmmsearch, then DIAMOND).
        "smafa": vendored_or_path("cargo", "smafa"), "mfqe": vendored_or_path("cargo", "mfqe"),
    }
    absent = [name for name, there in present.items() if not there]
    return absent + [tool for tool in missing_path_tools() if tool not in absent]


def catalogue_tools(repo: Path) -> bool:
    """The tools behind the five 2025 catalogue lanes and targeted assembly."""
    return not catalogue_tools_missing(repo)


def catalogue_refs(repo: Path) -> bool:
    """Every catalogue lane's references, by each engine's own readiness test."""
    from openbiota.engines import (  # noqa: PLC0415
        kraken,
        metaphlan_jan26,
        motus,
        singlem,
        sylph,
        sylph_globdb,
    )

    refs = repo / "refs"
    try:
        return all((
            metaphlan_jan26.database_ready(refs / "metaphlan4_db"),
            motus.database_ready(refs),
            kraken.database_ready(refs),
            singlem.ensure_extracted(refs),
            sylph_globdb.database_ready(refs),
            (refs / "sylph" / sylph.DATABASE_FILE).is_file()
            and (refs / "sylph" / "tax" / "gtdb_r232_metadata.tsv.gz").is_file(),
        ))
    except Exception:  # noqa: BLE001 - an unreadable reference is a missing one
        return False


#: The standard install: what a lite report needs.
CHECKS: dict[str, Check] = {
    "python": python_env,
    "profiler": profiler,
    "data": data,
    "metaphlan_db": metaphlan_db,
    "diamond_db": diamond_db,
    "models": models,
    "host_index": host_index,
    "catalogue_tools": catalogue_tools,
    "catalogue_refs": catalogue_refs,
}


#: What each standard-install piece is, in words, for messages to a person.
LABELS: dict[str, str] = {
    "python": "the Python environment",
    "profiler": "MetaPhlAn 3",
    "data": "the reference data that ships in data/",
    "metaphlan_db": "the MetaPhlAn marker database",
    "diamond_db": "the protein search database",
    "models": "the metabolic model library",
    "host_index": "the human-read filter index",
    "catalogue_tools": "the catalogue lanes' tools",
    "catalogue_refs": "the catalogue lanes' references",
}


def describe_missing(names: list[str], repo: Path | None = None) -> str:
    """Missing pieces in words; for the tools, which ones: "the catalogue lanes' tools (hmmsearch)"."""
    root = repo or REPO_ROOT
    parts: list[str] = []
    for name in names:
        label = LABELS.get(name, name)
        if name == "catalogue_tools":
            try:
                which = catalogue_tools_missing(root)
            except Exception:  # noqa: BLE001 - a description must not fail the message it is for
                which = []
            if which:
                label += f" ({', '.join(which)})"
        parts.append(label)
    return ", ".join(parts)


def standard_missing(repo: Path | None = None) -> list[str]:
    """The standard-install pieces that are not present, in install order.

    ``host_index`` is left out: the installer's --skip-host-index is a
    supported choice and the run only loses the human-DNA fraction.
    """
    root = repo or REPO_ROOT
    return [name for name, check in CHECKS.items() if name != "host_index" and not check(root)]


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
    """The full install's further tools: skani, fastANI and the PhyloPhlAn
    toolchain StrainPhlAn needs (raxmlHPC, mafft, trimal, FastTree)."""
    import shutil  # noqa: PLC0415

    return all((
        shutil.which("skani") is not None,
        shutil.which("fastANI") is not None or shutil.which("fastani") is not None,
        (repo / "vendor" / "bin" / "raxmlHPC").exists() or shutil.which("raxmlHPC") is not None,
        shutil.which("mafft") is not None,
    ))


def globdb_genomes(repo: Path) -> bool:
    """The 272 GB GlobDB genome archive, unpacked: what confirms GlobDB-native clusters."""
    return (repo / "refs" / "globdb_r232" / "genomes" / "index.tsv").is_file()


def expansion(repo: Path) -> bool:
    return (repo / "refs" / "expanded" / "reconciliation" / "summary.json").is_file()


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
    "globdb_genomes": globdb_genomes,
    "expansion": expansion,
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
    from openbiota import projectpath  # noqa: PLC0415

    projectpath.ensure(REPO_ROOT)
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
            print(f"  {name:<16}{'ready' if ready else 'MISSING'}")
            worst |= 0 if ready else 1
        return worst
    return 0 if ALL_CHECKS[args.name](args.repo) else 1


if __name__ == "__main__":
    sys.exit(main())
