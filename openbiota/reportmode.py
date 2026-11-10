"""Full and lite: the two reports this software produces.

The standard install (``./install.sh``, about 190 GB) provides everything
for a **lite** report: the community, the microbial groups, the organisms
that need attention - named across the five 2025 catalogue lanes (MetaPhlAn
4 Jan26, mOTUs 4, Kraken 2 / UHGG, SingleM, sylph / GlobDB) - the metabolite
pathways, the disease-pattern resemblance, the microbiome age, biofilm, the
gut-skin axis and what the research says about each finding. The **full**
install (``./install.sh --full``, about 1.1 TB more) adds the pathogen
screen, the mycobiome, strain typing and placement, and the GlobDB genome
archive that confirms GlobDB-native clusters.

Before this module existed a run on the standard install attempted every
lane, printed "not installed" for the ones it lacked, rendered whole
sections whose only content was "this module did not run", and ended with
a red INCOMPLETE RUN. Now a run has a mode. In lite mode the full-only
lanes are not attempted, their sections are not rendered, the remaining
sections are numbered contiguously, and the completeness census asks only
about the stages lite promises. In full mode everything is required, and a
run refuses to start on a machine that does not have the full install
rather than spend hours producing an incomplete report.

One module defines the split so the pipeline, the census and the PDF
cannot disagree about it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

MODES: Final = ("full", "lite")

#: Census stage keys (see :mod:`openbiota.completeness`) that only the full
#: install can satisfy. Everything else is required in both modes.
FULL_ONLY_STAGES: Final = frozenset({
    # The Jun23 catalogue and the marker-consensus strain typing built on it,
    # StrainPhlAn placement, the pathogen screen and the mycobiome. The five
    # 2025 catalogue lanes, the whole-genome profile, confirmation and
    # targeted assembly are part of every report.
    "extended_catalogue", "strain_typing", "strain_analysis", "pathogens", "mycobiome",
})

#: Report sections (keys of ``pdfreport.SECTIONS``) that exist only in a
#: full report. A lite report leaves them out entirely - no heading, no
#: "did not run" page - and numbers the rest without gaps.
FULL_ONLY_SECTIONS: Final = frozenset({
    "pathogens", "pathogens_detail", "mycobiome", "mycobiome_detail", "strains",
})


@dataclass(frozen=True)
class Mode:
    name: str
    #: One line, for the interactive prompt and the help text.
    summary: str
    #: What a reader gets, in the words the report uses.
    sections: str
    #: Wall-clock on a desktop, first run of a sample / the same sample again.
    time_first: str
    time_again: str


FULL: Final = Mode(
    "full",
    "every section",
    "the pathogen screen (485 targets), the mycobiome, strain typing and strain placement, "
    "plus everything in lite",
    "two and a half to three hours", "about ten minutes",
)
LITE: Final = Mode(
    "lite",
    "the standard install's report",
    "your community, every organism the five 2025 catalogues can name, the organisms that need "
    "attention, 45 metabolite pathways, disease patterns, microbiome age and what the research "
    "says you can do",
    "one to two hours", "about five minutes",
)
BY_NAME: Final = {FULL.name: FULL, LITE.name: LITE}

#: What the full install costs beyond the standard one, for the prompt that offers it.
FULL_INSTALL_COST: Final = "about 1.1 TB more reference data and most of a day"


def excluded_stages(mode: str) -> frozenset[str]:
    """Census keys a run in this mode is not asked about."""
    return FULL_ONLY_STAGES if mode == "lite" else frozenset()


def omitted_sections(mode: str) -> frozenset[str]:
    """Section keys a report in this mode does not carry."""
    return FULL_ONLY_SECTIONS if mode == "lite" else frozenset()


def standard_install_missing(repo: Path | None = None) -> list[str]:
    """The standard-install pieces this machine lacks; empty when complete."""
    from openbiota.installcheck import standard_missing  # noqa: PLC0415

    return standard_missing(repo)


def full_install_missing(repo: Path | None = None) -> list[str]:
    """The full-install pieces this machine lacks; empty when it has them all."""
    from openbiota.installcheck import full_missing  # noqa: PLC0415

    return full_missing(repo)


def resolve(requested: str | None, repo: Path | None = None) -> tuple[str, str]:
    """Decide the mode for a run that did not ask interactively.

    Returns ``(mode, reason)``. ``None`` or ``"auto"`` picks full when the
    full install is present and lite otherwise, and says which it chose and
    why; an explicit mode is returned as given (full is checked by the
    caller, which can offer the installer).
    """
    if requested in MODES:
        return requested, f"--mode {requested}"
    missing = full_install_missing(repo)
    if not missing:
        return "full", "the full install is present"
    shown = ", ".join(missing[:3]) + (f" and {len(missing) - 3} more" if len(missing) > 3 else "")
    return "lite", f"the full install is not complete on this machine (missing: {shown})"


def mode_of(results: Mapping[str, Any]) -> str:
    """The mode a finished run was made in.

    Every run records it. A results file without the record predates modes,
    when the only run there was attempted everything, so it is read as full:
    the census stays strict about it and the report renders as it did then.
    """
    recorded = ((results.get("run") or {}).get("report_mode") or "").strip().lower()
    return recorded if recorded in MODES else "full"


def label(mode: str) -> str:
    """The word the report prints after its title."""
    return "full" if mode == "full" else "lite"
