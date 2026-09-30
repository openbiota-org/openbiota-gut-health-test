"""Pinned adapters for the external strain and comparison tools.

Each adapter declares the tool, its pinned version, what object it estimates,
and its installation state. An uninstalled adapter reports
`reference_unavailable` for its targets; it never degrades into a species-level
answer and never produces a negative (spec §11.3).

One correction here is mandatory rather than stylistic. **TRACS 1.1.4 ships
default clock and transmission parameters derived from SARS-CoV-2 examples**,
and computes transmission quantities when `--meta` is supplied. Handing it
specimen dates would make it emit bacterial "transmission" numbers from a viral
substitution rate. So the adapter below omits `--meta` unconditionally, keeps
dates in OpenBiota's own provenance, and consumes only SNP-distance and
callability outputs (§11.4, acceptance 14).

The other rule worth stating: these tools estimate **different objects**. A
marker phylogeny, a genome-wide population comparison and a reference-mixture
deconvolution are not three measurements of one quantity, so there is no
"two tools agree" bonus. Discordance is preserved (§4.2, acceptance 16).
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from .schema import (
    REFERENCE_UNAVAILABLE,
    SCHEDULED,
    UNRESOLVED,
    ResolutionCall,
    Validation,
)

#: Virtualenvs `make strain-tools` installs into. Separate venvs because the
#: tools have incompatible dependency pins (see the Makefile).
VENVS: Final[tuple[str, ...]] = (
    ".venv-strain", ".venv-instrain", ".venv-mpa4", ".venv"
)

@dataclass(frozen=True, slots=True)
class ToolAdapter:
    """One external tool, pinned, with the object it actually estimates."""

    adapter_id: str
    tool: str
    pinned_version: str
    #: The quantity this tool estimates. Distinct per tool on purpose.
    estimates: str
    #: Executable to probe on PATH, or None when it is a library.
    executable: str | None
    licence: str
    #: When set, resolve *only* inside this virtualenv and ignore PATH. Used
    #: where a wrong-version binary on PATH would be actively dangerous: a
    #: MetaPhlAn 4.2 `strainphlan` cannot read 4.1 consensus markers, and
    #: finding it first would silently produce incomparable results (§4.3).
    pinned_venv: str | None = None
    #: Flags that must never be passed, and why.
    forbidden_flags: Mapping[str, str] = field(default_factory=dict)
    #: Flags the default invocation uses.
    default_flags: tuple[str, ...] = ()
    limits: tuple[str, ...] = ()
    source: str = ""

    @property
    def path(self) -> str | None:
        """Where this tool actually is, PATH or a pinned venv.

        `make strain-tools` installs into per-tool virtualenvs rather than
        PATH, because several of these have mutually incompatible pins —
        inStrain needs Python 3.9 and biopython<=1.74 while everything else
        needs a modern interpreter. So detection has to look in those venvs,
        or the census would report installed tools as unavailable.
        """
        if self.executable is None:
            return None
        root = Path(__file__).resolve().parents[2]
        if self.pinned_venv is not None:
            candidate = root / self.pinned_venv / "bin" / self.executable
            return str(candidate) if candidate.exists() else None
        if found := shutil.which(self.executable):
            return found
        for venv in VENVS:
            candidate = root / venv / "bin" / self.executable
            if candidate.exists():
                return str(candidate)
        # Vendored tools sit outside the venvs (see the Makefile).
        for pattern in (f"vendor/*/{self.executable}", f"vendor/*/bin/{self.executable}"):
            for candidate in sorted(root.glob(pattern)):
                return str(candidate)
        return None

    @property
    def installed(self) -> bool:
        return self.path is not None

    def command(self, *args: str) -> list[str]:
        """Build an argv array, refusing any forbidden flag.

        Argv arrays rather than shell strings, so an identifier or path can
        never be interpolated into a shell (spec §11.1).
        """
        argv = [self.path or self.executable or self.tool, *self.default_flags, *args]
        for flag, reason in self.forbidden_flags.items():
            if flag in argv:
                raise ValueError(
                    f"{self.adapter_id}: refusing to pass {flag} \u2014 {reason}"
                )
        return argv

    def to_json(self) -> dict[str, Any]:
        return {
            "adapter_id": self.adapter_id,
            "tool": self.tool,
            "pinned_version": self.pinned_version,
            "estimates": self.estimates,
            "installed": self.installed,
            "path": self.path,
            "licence": self.licence,
            "forbidden_flags": dict(self.forbidden_flags),
            "default_flags": list(self.default_flags),
            "limits": list(self.limits),
            "source": self.source,
        }


ADAPTERS: Final[tuple[ToolAdapter, ...]] = (
    ToolAdapter(
        adapter_id="strainphlan",
        tool="StrainPhlAn",
        pinned_version="4.1.1",
        estimates="dominant marker-gene consensus and phylogenetic placement",
        # Lives in the pinned 4.1.1 venv, and must be resolved there rather
        # than on PATH: the repo's main venv carries MetaPhlAn 4.2, whose
        # strainphlan cannot read 4.1 consensus markers.
        executable="strainphlan",
        pinned_venv=".venv-mpa4",
        licence="MIT",
        default_flags=(),
        limits=(
            "Not exhaustive strain-mixture deconvolution: a dominant consensus can conceal "
            "a minority population.",
            "Marker genes are a conserved subset; agreement across them is not whole-genome "
            "identity.",
            "The 4.1 and 4.2 workflows must not be mixed; this adapter is pinned to 4.1.1 "
            "with the vJun23 database.",
        ),
        source="https://github.com/biobakery/MetaPhlAn",
    ),
    ToolAdapter(
        adapter_id="tracs",
        tool="TRACS",
        pinned_version="1.1.4",
        estimates="conservative lower bound on pairwise SNP distance with callability",
        executable="tracs",
        licence="MIT (verify pinned tag)",
        forbidden_flags={
            "--meta": (
                "TRACS 1.1.4 computes transmission quantities from default clock parameters "
                "derived from SARS-CoV-2 examples. Supplying specimen metadata would turn "
                "dates into assumed bacterial transmission rates. Dates stay in OpenBiota's "
                "provenance and only SNP-distance/callability outputs are consumed."
            ),
        },
        default_flags=(),
        limits=(
            "A conservative lower bound, not a haplotype reconstruction.",
            "Small distance with inadequate callable overlap is not identity.",
            "References below the configured depth or with excessive ambiguity are skipped; "
            "a skipped reference is insufficient resolution, not biological non-detection.",
        ),
        source="https://github.com/gtonkinhill/tracs",
    ),
    ToolAdapter(
        adapter_id="midas",
        tool="MIDAS v3",
        pinned_version="midasv3=1.0.0",
        estimates="pangenome gene presence and within-species SNVs (unphased)",
        executable="midas",
        licence="permissive (verify pinned commit)",
        limits=(
            "Gene and SNP profiles do not phase every gene into every strain; a gene union "
            "is not one strain's genome.",
            "Broad pangenome pruning can drop singleton families: curated targets must be "
            "retained by the independent targeted lane.",
        ),
        source="https://github.com/pollardlab/MIDAS",
    ),
    ToolAdapter(
        adapter_id="instrain",
        tool="inStrain",
        pinned_version="1.10.0",
        estimates="population microdiversity and pairwise popANI over compared bases",
        executable="inStrain",
        licence="MIT",
        limits=(
            "High population ANI over a small conserved portion of a genome is not a "
            "named-strain match; compared bases and masking must be reported.",
            "Runs on its own Python 3.9 interpreter: it requires biopython<=1.74 for "
            "Bio.codonalign.codonalphabet, which later releases removed.",
        ),
        source="https://github.com/MrOlm/inStrain",
    ),
    ToolAdapter(
        adapter_id="strainge",
        tool="StrainGE",
        pinned_version="1.3.9",
        estimates="closest-reference components and variation relative to them",
        executable="straingst",
        licence="BSD-3-Clause",
        limits=(
            "Reference-nearest is not exact-isolate identity.",
            "Low-depth performance from selected experiments is not a universal detection "
            "limit.",
        ),
        source="https://github.com/broadinstitute/StrainGE",
    ),
    ToolAdapter(
        adapter_id="stxtyper",
        tool="NCBI StxTyper",
        pinned_version="1.0.45",
        estimates="Shiga-toxin operon completeness and stx subtype",
        executable="stxtyper",
        licence="US government / public domain",
        limits=(
            "Operon completeness, frameshifts, truncations and novel states are distinct "
            "results; a one-gene positive/negative cannot replace them.",
            "Built from its own release tag: AMRFinderPlus bundles it as a git submodule "
            "that release tarballs omit.",
        ),
        source="https://github.com/ncbi/stxtyper",
    ),
    ToolAdapter(
        adapter_id="amrfinderplus",
        tool="NCBI AMRFinderPlus",
        pinned_version="4.2.7",
        estimates="curated acquired AMR genes, resistance substitutions and virulence loci",
        executable="amrfinder",
        licence="US government / public domain",
        limits=(
            "Genomic potential, not measured susceptibility: absence of catalogued genes "
            "does not establish that an organism is susceptible.",
            "Organism-specific mutation interpretation needs a justified host identity and "
            "discriminating coverage.",
            "A reported gene is not attributed to a carrier organism without linkage "
            "evidence; coordinates locate it on a contig, not in a strain.",
        ),
        source="https://github.com/ncbi/amr",
    ),
    ToolAdapter(
        adapter_id="samestr",
        tool="SameStr",
        pinned_version="pin from source",
        estimates="shared-strain calls from marker allele profiles",
        executable="samestr",
        licence="AGPL-3.0 (the MIT workflow wrapper does not relicense the engine)",
        limits=("Sharing is not disease transfer and not guaranteed engraftment.",),
        source="https://github.com/danielpodlesny/samestr",
    ),
)

BY_ADAPTER: Final[Mapping[str, ToolAdapter]] = {a.adapter_id: a for a in ADAPTERS}


def tracs_distance_command(
    *, msa: str, out: str, threads: int = 8, extra: Sequence[str] = ()
) -> list[str]:
    """The only sanctioned TRACS distance invocation.

    `--meta` is refused by construction, so no code path can produce
    transmission probabilities from viral-rate defaults.
    """
    adapter = BY_ADAPTER["tracs"]
    return adapter.command(
        "distance", "--msa", msa, "-o", out, "--filter", "-t", str(threads), *extra
    )


def calls_for_census(sample_id: str) -> list[ResolutionCall]:
    """One record per adapter, installed or not.

    Both states have to appear. An uninstalled adapter is a missing
    capability; an installed one that nothing invoked is a pending
    measurement. Only the first used to be emitted, so once the tools were
    installed they vanished from the census entirely - and the census
    exists precisely to stop "we have this and did not run it" from
    disappearing. Neither state is a negative result.
    """
    out: list[ResolutionCall] = []
    for adapter in ADAPTERS:
        if adapter.installed:
            status = SCHEDULED
            reason = "adapter_installed_not_run"
            plain = (
                f"{adapter.tool} ({adapter.pinned_version}) estimates "
                f"{adapter.estimates}. It is installed but was not run on this "
                "sample, so that object was not measured. Nothing here says it "
                "is absent."
            )
        else:
            status = REFERENCE_UNAVAILABLE
            reason = "adapter_not_installed"
            plain = (
                f"{adapter.tool} ({adapter.pinned_version}) estimates "
                f"{adapter.estimates}. It is not installed, so that object was not "
                "measured. This is a missing capability, not a negative finding, and "
                "it does not affect the analyses that did run."
            )
        out.append(
            ResolutionCall(
                sample_id=sample_id,
                target_id=f"adapter.{adapter.adapter_id}",
                identity_kind="strain_component",
                assay_status=status,
                analytical_call=UNRESOLVED,
                reason_codes=(reason,),
                validation=Validation(
                    analytical_status="not_validated",
                    reference_function_status="tool_defined",
                    phenotype_measured_in_sample=False,
                    clinical_predictive_status="not_established",
                ),
                plain=plain,
            )
        )
    return out


#: Retained name for callers that only wanted the uninstalled ones.
def calls_for_unavailable(sample_id: str) -> list[ResolutionCall]:
    """Only the uninstalled adapters. Prefer :func:`calls_for_census`."""
    return [
        c for c in calls_for_census(sample_id)
        if c.assay_status == REFERENCE_UNAVAILABLE
    ]


def registry_json() -> dict[str, Any]:
    return {
        "adapters": [a.to_json() for a in ADAPTERS],
        "n_adapters": len(ADAPTERS),
        "n_installed": sum(1 for a in ADAPTERS if a.installed),
        "no_consensus_bonus": (
            "These tools estimate different objects: a marker phylogeny, a genome-wide "
            "population comparison and a reference-mixture deconvolution are not three "
            "measurements of one quantity. Agreement between them is not multiplied into "
            "extra certainty, and discordance is preserved rather than averaged away."
        ),
    }


def self_test() -> int:
    """Guard acceptance tests 14-16. Returns a failure count."""
    failures = 0

    # Acceptance 14: TRACS must never receive --meta.
    tracs = BY_ADAPTER["tracs"]
    if "--meta" not in tracs.forbidden_flags:
        failures += 1
    try:
        tracs.command("distance", "--meta", "dates.csv")
        failures += 1
    except ValueError as exc:
        if "--meta" not in str(exc):
            failures += 1
    cmd = tracs_distance_command(msa="aln.fasta", out="d.csv")
    if "--meta" in cmd:
        failures += 1
    if "--filter" not in cmd:
        failures += 1
    # Argv array, never a shell string.
    if not isinstance(cmd, list) or any(" " in part for part in cmd if part.startswith("-")):
        failures += 1

    # Acceptance 15: a skipped reference is insufficient resolution.
    if not any("not biological non-detection" in limit for limit in tracs.limits):
        failures += 1

    # Acceptance 16: no agreement bonus, and each tool names a distinct object.
    reg = registry_json()
    if "not multiplied into" not in reg["no_consensus_bonus"]:
        failures += 1
    estimates = [a.estimates for a in ADAPTERS]
    if len(set(estimates)) != len(estimates):
        failures += 1

    # Uninstalled adapters produce non-negative records with reasons.
    for call in calls_for_unavailable("K1"):
        if call.assay_status != REFERENCE_UNAVAILABLE:
            failures += 1
        if call.is_negative:
            failures += 1
        if "not a negative finding" not in call.plain:
            failures += 1

    # Every adapter states its limits and licence.
    for adapter in ADAPTERS:
        if not adapter.limits or not adapter.licence:
            failures += 1

    # The version trap: a pinned-venv tool must never resolve to PATH.
    marker = BY_ADAPTER["strainphlan"]
    if marker.pinned_venv != ".venv-mpa4":
        failures += 1
    if marker.path is not None and ".venv-mpa4" not in marker.path:
        failures += 1
    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
