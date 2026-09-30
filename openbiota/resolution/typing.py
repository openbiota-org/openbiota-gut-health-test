"""Organism-specific typing: pathotypes, serotypes, sequence types, toxins.

Species names are the wrong unit for most of the organisms that matter here.
*E. coli* covers a harmless resident and O157:H7. *C. difficile* covers a
toxin-negative carriage strain and a fulminant one. So this module holds the
rules that turn locus evidence into a type — and, far more often, the rules
that refuse to.

The refusals are the substance, because every one of them is a real way that
mixed-specimen data fools isolate-derived tools (spec §6, acceptance 22-32):

* **Alleles from different strains cannot be combined.** An ST assembled from
  unphased alleles belonging to two populations is a sequence type that exists
  in no organism. Same for an O antigen from one population and an H antigen
  from another.
* **eae-only plus stx-only is not one strain.** A specimen containing an EPEC
  and a STEC must keep both possibilities, not invent a linked eae+stx hybrid.
* **A free phage carries stx too.** Shiga-toxin phage reads beside a commensal
  *E. coli* leave the carrier unresolved; they do not make that *E. coli*
  a STEC.
* **ipaH is shared.** It does not distinguish *Shigella* from
  enteroinvasive *E. coli*, and it certainly does not name a *Shigella*
  species.
* **mecA in a coagulase-negative staph is not MRSA**, and an unlinked *van*
  gene is not VRE.
* **Absence at low depth is not a negative criterion.** "stx-negative" is part
  of the EPEC definition, so it requires depth adequate to have seen stx.

Most tools in §6.1 were built for isolate assemblies. Accepting a FASTQ does
not make one validated on stool, so this module encodes the interpretation
layer and records which external typer would supply each call, with its status.
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from .schema import (
    AMBIGUOUS_HOMOLOGY,
    CANDIDATE_SEQUENCE,
    COMPLETED,
    INSUFFICIENT_DEPTH,
    MIXED_UNRESOLVED,
    NOT_REQUESTED,
    REFERENCE_UNAVAILABLE,
    SCHEDULED,
    UNRESOLVED,
    Coverage,
    Linkage,
    ResolutionCall,
    Validation,
)

#: Virtualenvs `make strain-tools` installs into.
_VENVS = (".venv-strain", ".venv-instrain", ".venv-mpa4", ".venv")

@dataclass(frozen=True, slots=True)
class TypingScheme:
    """An organism-specific typing capability and its real status."""

    scheme_id: str
    organism: str
    label: str
    #: The external tool that would supply this, per spec §6.1.
    tool: str
    #: Executable that proves the scheme's tool is present, resolved across
    #: the pinned virtualenvs. `None` means no single binary provides it, so
    #: the scheme stays unavailable until one is wired.
    executable: str | None = None
    #: Loci that must be present *and* linked for a positive type.
    required_loci: tuple[str, ...] = ()
    #: Loci that must be adequately assessed and negative.
    required_absent_loci: tuple[str, ...] = ()
    #: Why this scheme cannot be satisfied by species detection alone.
    species_insufficient_because: str = ""
    #: Notes on mixed-specimen behaviour.
    mixed_specimen_rule: str = ""
    validated_on_stool: bool = False

    @property
    def tool_installed(self) -> bool:
        """Is the scheme's tool actually present on this machine?

        Resolved across the pinned virtualenvs, because `make strain-tools`
        installs into those rather than onto PATH.
        """
        if self.executable is None:
            return False
        if shutil.which(self.executable):
            return True
        root = Path(__file__).resolve().parents[2]
        if any((root / venv / "bin" / self.executable).exists() for venv in _VENVS):
            return True
        # Vendored tools live outside any venv: BLAST+ and AMRFinderPlus are
        # installed under vendor/ because Homebrew cannot place them without
        # disturbing existing formulae (see the Makefile).
        return any(root.glob(f"vendor/*/{self.executable}")) or any(
            root.glob(f"vendor/*/bin/{self.executable}")
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "scheme_id": self.scheme_id,
            "organism": self.organism,
            "label": self.label,
            "tool": self.tool,
            "executable": self.executable,
            "tool_installed": self.tool_installed,
            "validated_on_mixed_stool": self.validated_on_stool,
            "required_loci": list(self.required_loci),
            "required_absent_loci": list(self.required_absent_loci),
            "species_insufficient_because": self.species_insufficient_because,
            "mixed_specimen_rule": self.mixed_specimen_rule,
        }


#: The mandatory initial families (spec §6.2). Availability is detected at
#: runtime; a scheme whose tool is absent drives `reference_unavailable` rather
#: than silence, and never a negative type.
SCHEMES: Final[tuple[TypingScheme, ...]] = (
    TypingScheme(
        scheme_id="ecoli.serotype",
        organism="Escherichia coli",
        label="O:H genoserotype",
        tool="ECTyper",
        executable="ectyper",
        required_loci=("O-antigen locus", "H-antigen (fliC) allele"),
        species_insufficient_because=(
            "E. coli spans commensals and every diarrhoeagenic pathotype; the species name "
            "carries no clinical information on its own."
        ),
        mixed_specimen_rule=(
            "O and H alleles from different populations must not be combined into one "
            "confident O:H type (acceptance 23)."
        ),
    ),
    TypingScheme(
        scheme_id="ecoli.pathotype",
        organism="Escherichia coli",
        label="Diarrhoeagenic pathotype (STEC/EPEC/ETEC/EAEC/EIEC)",
        tool="ECTyper pathotype rules + NCBI StxTyper",
        executable="stxtyper",
        required_loci=("pathotype-defining locus set with carrier linkage",),
        required_absent_loci=("stx1", "stx2"),
        species_insufficient_because=(
            "Pathotype is defined by virulence loci, not by species. Multiple pathotypes "
            "can coexist in one specimen."
        ),
        mixed_specimen_rule=(
            "eae-only and stx-only populations remain two possibilities; a linked eae+stx "
            "strain is never inferred (acceptance 24). Free Stx phage evidence stays "
            "carrier-unresolved (acceptance 25)."
        ),
    ),
    TypingScheme(
        scheme_id="shigella.identity",
        organism="Shigella / enteroinvasive E. coli",
        label="Shigella species versus EIEC",
        tool="ShigaTyper / ShigaPass",
        executable=None,  # neither is packaged for this platform yet
        required_loci=("species-discriminating loci beyond ipaH",),
        species_insufficient_because=(
            "Shigella and E. coli are one genomic species; ipaH is carried by both "
            "Shigella and EIEC."
        ),
        mixed_specimen_rule="ipaH alone never becomes a named Shigella species (acceptance 26).",
    ),
    TypingScheme(
        scheme_id="cdiff.toxinotype",
        organism="Clostridioides difficile",
        label="PaLoc/CdtLoc architecture and toxin subtype",
        tool="DiffBase",
        executable=None,  # reference set needs curation before a call
        required_loci=("tcdA", "tcdB", "PaLoc architecture"),
        species_insufficient_because=(
            "Carriage of toxin-negative C. difficile is common and is not disease; the "
            "toxin locus is the finding."
        ),
        mixed_specimen_rule=(
            "TcdA-negative/TcdB-positive strains must stay detectable, and toxin subtype is "
            "not ribotype (acceptance 32)."
        ),
    ),
    TypingScheme(
        scheme_id="klebsiella.loci",
        organism="Klebsiella pneumoniae complex",
        label="Species, ST, K/O loci, virulence and AMR",
        tool="Kleborate + Kaptive",
        executable="kleborate",
        required_loci=("coherent lineage with linked K/O and virulence loci",),
        species_insufficient_because=(
            "Capsule type, virulence plasmid content and AMR vary independently within the "
            "complex."
        ),
        mixed_specimen_rule=(
            "Capsule, AMR and virulence loci from different strains must not be merged into "
            "one artificial hypervirulent resistant strain (acceptance 28)."
        ),
    ),
    TypingScheme(
        scheme_id="staph.mrsa",
        organism="Staphylococcus spp.",
        label="mecA/mecC with SCCmec and host species",
        tool="AMRFinderPlus + species resolution",
        executable="amrfinder",
        required_loci=("mecA or mecC", "S. aureus host identity with linkage"),
        species_insufficient_because=(
            "mecA is common in coagulase-negative staphylococci, where it is not MRSA."
        ),
        mixed_specimen_rule=(
            "mecA in a non-aureus staphylococcus is never MRSA; an unlinked van gene is "
            "never VRE faecium (acceptance 29)."
        ),
    ),
    TypingScheme(
        scheme_id="salmonella.serovar",
        organism="Salmonella enterica",
        label="Serovar determinants and core lineage",
        tool="SeqSero2 / SISTR",
        executable=None,  # not packaged for this platform yet
        required_loci=("serovar determinant alleles with a coherent lineage",),
        species_insufficient_because="Serovars differ enormously in clinical meaning.",
        mixed_specimen_rule="Multiple serovars are separated or returned as a mixture.",
    ),
    TypingScheme(
        scheme_id="bfragilis.bft",
        organism="Bacteroides fragilis",
        label="bft subtype and BfPAI context",
        tool="targeted allele mapping",
        # The bft1/2/3 allele references are now on disk (refs/panels/
        # mechanisms), so this needs the mapping step wired, not a new tool.
        executable=None,
        required_loci=("bft1/bft2/bft3 discriminating CDS", "BfPAI context"),
        species_insufficient_because=(
            "Enterotoxigenic and nontoxigenic B. fragilis are the same species; "
            "metalloprotease homologs are decoys."
        ),
        mixed_specimen_rule="A gene-family fragment count does not establish a subtype.",
    ),
    TypingScheme(
        scheme_id="colibactin.island",
        organism="Colibactin-capable Enterobacteriaceae",
        label="clbA-S island architecture and carrier",
        tool="targeted locus mapping against AM229678.1",
        # The 55,140-bp pks island reference is on disk (refs/panels/
        # mechanisms); this needs the locus-architecture step wired.
        executable=None,
        required_loci=("clbA-S architecture", "carrier organism linkage"),
        species_insufficient_because=(
            "The carrier may be E. coli, Klebsiella, Citrobacter or another organism."
        ),
        mixed_specimen_rule="clbB or clbS alone is not an intact island (acceptance 31).",
    ),
)

BY_SCHEME: Final[Mapping[str, TypingScheme]] = {s.scheme_id: s for s in SCHEMES}


@dataclass(frozen=True, slots=True)
class LocusObservation:
    """One locus as actually observed in the specimen."""

    locus: str
    detected: bool
    #: Independent fragments supporting it, or None when unassessed.
    fragments: int | None = None
    #: Whether the locus was linked to a carrier organism.
    linked_to_carrier: bool = False
    #: Whether depth was adequate to read a non-detection as absence.
    depth_adequate_for_absence: bool = False
    #: Populations the locus could belong to, when more than one is present.
    possible_carriers: tuple[str, ...] = ()


def call_type(
    *,
    sample_id: str,
    scheme_id: str,
    observations: Sequence[LocusObservation],
    populations_of_species: int = 1,
    assay_ran: bool = False,
) -> ResolutionCall:
    """Apply a scheme's rules to observed loci and return a typed call.

    The return is almost always a refusal with a reason, and that is correct:
    short reads from a mixed specimen rarely support a confident type.

    ``assay_ran`` must be set by the caller when the scheme's tool was
    actually invoked on this sample. It defaults to False because the
    dangerous default is the other one. An empty ``observations`` sequence
    is ambiguous on its own - it means either "the tool ran and found no
    loci" or "no tool ever ran" - and those are opposite statements. Without
    this flag the registry census, which enumerates schemes rather than
    running them, produced calls reading ``completed`` and "required loci
    not detected" for tools that were never executed.
    """
    scheme = BY_SCHEME.get(scheme_id)
    if scheme is None:
        raise KeyError(f"unknown typing scheme {scheme_id!r}")

    by_locus = {o.locus: o for o in observations}
    reasons: list[str] = []
    identity_kind = (
        "serotype" if "serotype" in scheme_id
        else "sequence_type" if "st" in scheme_id.split(".")[-1]
        else "pathotype"
    )

    # The scheme's own tool is not installed: the call is unavailable, and
    # unavailable is not negative.
    if not scheme.tool_installed:
        return ResolutionCall(
            sample_id=sample_id,
            target_id=f"typing.{scheme_id}",
            identity_kind=identity_kind,
            assay_status=REFERENCE_UNAVAILABLE,
            analytical_call=UNRESOLVED,
            reason_codes=("typing_tool_not_installed",),
            plain=(
                f"{scheme.label}: the typing scheme is registered but {scheme.tool} is not "
                "installed, so the type was not determined. That is a missing capability, "
                "not a negative result. "
                + scheme.species_insufficient_because
            ),
            validation=Validation(
                analytical_status="not_validated",
                reference_function_status="scheme_defined",
                phenotype_measured_in_sample=False,
                clinical_predictive_status="not_established",
            ),
        )

    # The tool exists but has not been pointed at this sample. That is a
    # pending assay, not a negative one: reporting "required loci not
    # detected" here would be asserting the result of a test nobody ran.
    if not assay_ran:
        return ResolutionCall(
            sample_id=sample_id,
            target_id=f"typing.{scheme_id}",
            identity_kind=identity_kind,
            assay_status=SCHEDULED,
            analytical_call=UNRESOLVED,
            reason_codes=("typing_tool_installed_not_run",),
            plain=(
                f"{scheme.label}: {scheme.tool} is installed but was not run on "
                "this sample, so no type was determined. Nothing here says the "
                "loci are absent - they were not looked for. "
                + scheme.species_insufficient_because
            ),
            validation=Validation(
                analytical_status="not_validated",
                reference_function_status="scheme_defined",
                phenotype_measured_in_sample=False,
                clinical_predictive_status="not_established",
            ),
        )

    # More than one population of the species: alleles cannot be phased, so a
    # coherent type cannot be assembled from them.
    if populations_of_species > 1:
        reasons.append("multiple_conspecific_populations_unphased")

    missing = [
        locus for locus in scheme.required_loci
        if not (by_locus.get(locus) and by_locus[locus].detected)
    ]
    unlinked = [
        locus for locus in scheme.required_loci
        if by_locus.get(locus) and by_locus[locus].detected
        and not by_locus[locus].linked_to_carrier
    ]
    # A required-absent locus needs depth adequate to have seen it.
    shallow_negatives = [
        locus for locus in scheme.required_absent_loci
        if not (by_locus.get(locus) and by_locus[locus].depth_adequate_for_absence)
    ]
    ambiguous_carrier = [
        o.locus for o in observations if o.detected and len(o.possible_carriers) > 1
    ]

    if missing:
        reasons.append("required_loci_not_detected")
    if unlinked:
        reasons.append("required_loci_not_linked_to_a_carrier")
    if shallow_negatives:
        reasons.append("negative_criterion_depth_inadequate")
    if ambiguous_carrier:
        reasons.append("carrier_ambiguous_between_populations")

    if not reasons:
        call = CANDIDATE_SEQUENCE
        plain = (
            f"{scheme.label}: every required locus was detected and linked to one carrier, "
            "so a type is supported by the sequence. It remains a genotype prediction: "
            "expression, viability and clinical significance are not measured here."
        )
        status = COMPLETED
    else:
        if "multiple_conspecific_populations_unphased" in reasons:
            call = MIXED_UNRESOLVED
        elif "carrier_ambiguous_between_populations" in reasons:
            call = AMBIGUOUS_HOMOLOGY
        elif "negative_criterion_depth_inadequate" in reasons:
            call = INSUFFICIENT_DEPTH
        else:
            call = UNRESOLVED
        parts = [f"{scheme.label}: no type is asserted."]
        if missing:
            parts.append(f"Required loci not detected: {', '.join(missing)}.")
        if unlinked:
            parts.append(
                f"Detected but not linked to a carrier: {', '.join(unlinked)}. "
                "Two genes in one stool specimen are not automatically in one organism."
            )
        if shallow_negatives:
            parts.append(
                f"The definition requires {', '.join(shallow_negatives)} to be absent, but "
                "depth was not adequate to read non-detection as absence."
            )
        if ambiguous_carrier:
            parts.append(
                f"Carrier ambiguous for {', '.join(ambiguous_carrier)} between "
                f"{populations_of_species} populations."
            )
        if "multiple_conspecific_populations_unphased" in reasons:
            parts.append(
                "More than one population of this species is present and the alleles are "
                "unphased, so combining them would produce a type that exists in no "
                "organism."
            )
        parts.append(scheme.mixed_specimen_rule)
        plain = " ".join(p for p in parts if p)
        status = COMPLETED

    return ResolutionCall(
        sample_id=sample_id,
        target_id=f"typing.{scheme_id}",
        identity_kind=identity_kind,
        assay_status=status,
        call_id=f"{sample_id}:typing.{scheme_id}" if status == COMPLETED else None,
        analytical_call=call,
        reason_codes=tuple(reasons),
        coverage=Coverage(
            independent_fragments=sum(
                o.fragments or 0 for o in observations if o.detected
            ) or None,
        ),
        linkage=Linkage(
            state=(
                "read_pair_or_amplicon_linkage"
                if not unlinked and not missing else "unassigned"
            )
        ),
        validation=Validation(
            analytical_status=(
                "not_validated_on_mixed_specimen" if not scheme.validated_on_stool
                else "validated"
            ),
            reference_function_status="scheme_defined",
            phenotype_measured_in_sample=False,
            clinical_predictive_status="not_established",
        ),
        plain=plain,
    )


def registry_json() -> dict[str, Any]:
    installed = [s for s in SCHEMES if s.tool_installed]
    return {
        "schemes": [s.to_json() for s in SCHEMES],
        "n_schemes": len(SCHEMES),
        "n_installed": len(installed),
        "n_reference_unavailable": len(SCHEMES) - len(installed),
        "note": (
            "Registered typing schemes and whether their tool is installed. A scheme "
            "without its tool reports reference_unavailable for every sample; it never "
            "reports a negative type, and it never allows a species name to stand in for "
            "a type."
        ),
    }


def self_test() -> int:
    """Guard acceptance tests 22-32. Returns a failure count."""
    failures = 0

    # An uninstalled scheme is unavailable, never negative. `shigella.identity`
    # has no binary wired, so it exercises that branch whatever is installed.
    call = call_type(sample_id="K1", scheme_id="shigella.identity", observations=[])
    if call.assay_status != REFERENCE_UNAVAILABLE:
        failures += 1
    if call.is_negative:
        failures += 1
    if "not a negative result" not in call.plain:
        failures += 1

    # The rule branches need an installed scheme. ECTyper backs the E. coli
    # pathotype scheme, so use it when present and skip these checks when it
    # is genuinely absent rather than faking installation.
    if BY_SCHEME["ecoli.pathotype"].tool_installed:
        # Acceptance 27: a shallow stx non-detection cannot satisfy
        # "stx-negative".
        call = call_type(
            sample_id="K1", scheme_id="ecoli.pathotype",
            observations=[
                LocusObservation(
                    "pathotype-defining locus set with carrier linkage",
                    detected=True, fragments=40, linked_to_carrier=True,
                ),
                LocusObservation("stx1", detected=False, depth_adequate_for_absence=False),
                LocusObservation("stx2", detected=False, depth_adequate_for_absence=False),
            ],
        )
        if call.analytical_call != INSUFFICIENT_DEPTH:
            failures += 1
        if "negative_criterion_depth_inadequate" not in call.reason_codes:
            failures += 1

        # Acceptance 22/23: unphased alleles from two populations cannot type.
        call = call_type(
            sample_id="K1", scheme_id="ecoli.pathotype",
            observations=[
                LocusObservation(
                    "pathotype-defining locus set with carrier linkage",
                    detected=True, fragments=40, linked_to_carrier=True,
                ),
                LocusObservation("stx1", detected=False, depth_adequate_for_absence=True),
                LocusObservation("stx2", detected=False, depth_adequate_for_absence=True),
            ],
            populations_of_species=3,
        )
        if call.analytical_call != MIXED_UNRESOLVED:
            failures += 1
        if "exists in no" not in call.plain:
            failures += 1

        # Acceptance 25: an unlinked locus leaves the carrier unresolved.
        call = call_type(
            sample_id="K1", scheme_id="ecoli.pathotype",
            observations=[
                LocusObservation(
                    "pathotype-defining locus set with carrier linkage",
                    detected=True, fragments=5, linked_to_carrier=False,
                ),
                LocusObservation("stx1", detected=False, depth_adequate_for_absence=True),
                LocusObservation("stx2", detected=False, depth_adequate_for_absence=True),
            ],
        )
        if "required_loci_not_linked_to_a_carrier" not in call.reason_codes:
            failures += 1
        if call.linkage.state != "unassigned":
            failures += 1

        # All conditions met -> a candidate genotype, never a phenotype.
        call = call_type(
            sample_id="K1", scheme_id="ecoli.pathotype",
            observations=[
                LocusObservation(
                    "pathotype-defining locus set with carrier linkage",
                    detected=True, fragments=80, linked_to_carrier=True,
                ),
                LocusObservation("stx1", detected=False, depth_adequate_for_absence=True),
                LocusObservation("stx2", detected=False, depth_adequate_for_absence=True),
            ],
        )
        if call.analytical_call != CANDIDATE_SEQUENCE:
            failures += 1
        if call.validation.phenotype_measured_in_sample:
            failures += 1

    # Every scheme must explain why species is insufficient, and none may be
    # silently validated on stool.
    for s in SCHEMES:
        if not s.species_insufficient_because:
            failures += 1
        if s.validated_on_stool:
            failures += 1

    reg = registry_json()
    if reg["n_schemes"] != len(SCHEMES):
        failures += 1
    if reg["n_installed"] + reg["n_reference_unavailable"] != len(SCHEMES):
        failures += 1
    if NOT_REQUESTED != "not_requested":
        failures += 1
    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
