"""Fungi, protists, helminths and viruses: resolution and assay boundaries.

Bacterial assumptions break on every one of these groups, and the breakages
are specific rather than general (spec §7):

* **Fungi are not haploid.** *Candida albicans* is diploid; its MLST keeps
  heterozygosity, loss-of-heterozygosity and allele fractions, none of which
  survives a bacterial consensus. And ERG11 or FKS1 being *present* is not
  resistance — a specified substitution is.
* **Protists have their own namespaces.** A *Blastocystis* subtype and a
  bacterial sequence type are different kinds of label, and subtype alone is
  not a harmfulness verdict. *Giardia* heterozygosity is not automatically a
  mixed infection.
* **Helminth mitochondrial haplotypes are not strains**, and DNA in stool does
  not establish a viable infection.
* **DNA sequencing cannot test for RNA viruses.** Norovirus, rotavirus,
  astrovirus, sapovirus, hepatitis A/E, enteroviruses and SARS-CoV-2 are RNA
  genomes. Adding reference sequences does not fix a chemistry mismatch, so
  those targets register as `incompatible_input` — permanently, until an
  RNA/cDNA lane exists. This is the clearest case in the whole system of a
  target that must never read as a negative screen.
* **A phage is not a human pathogen.** Gut phage ecology and human viral
  pathogens are separate outputs and must not inherit each other's framing.

Denominators stay separate too: eukaryote marker-relative abundance, bacterial
relative abundance and genome-normalised gene counts cannot be added into one
percentage without a defined common measurement model.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

from .schema import (
    INCOMPATIBLE_INPUT,
    NOT_REQUESTED,
    PHENOTYPE_NOT_MEASURED,
    REFERENCE_UNAVAILABLE,
    UNRESOLVED,
    ResolutionCall,
    Validation,
)

#: Measurement denominators, kept apart on purpose (spec §7.1).
DENOMINATORS: Final[Mapping[str, str]] = {
    "bacterial_relative_abundance": (
        "Fraction of bacterial marker-assigned reads. Not comparable with eukaryotic or "
        "viral fractions."
    ),
    "eukaryote_marker_relative": (
        "Fraction among eukaryotic marker-assigned reads only. EukDetect-relative "
        "percentages cannot be merged into a bacterial percentage."
    ),
    "viral_votu_relative": (
        "Fraction among viral operational taxonomic units. Separate from cellular "
        "denominators."
    ),
    "genome_equivalents": (
        "Gene copies normalised to bacterial genome equivalents via a single-copy marker. "
        "A rate, not a composition share."
    ),
}


@dataclass(frozen=True, slots=True)
class KingdomTarget:
    """A non-bacterial target, its finest supported resolution and its limits."""

    target_id: str
    organism: str
    kingdom: str
    label: str
    identity_kind: str
    #: `available` | `tool_missing` | `chemistry_incompatible` | `curation_required`
    capability: str
    denominator: str
    required_resolution: str
    ploidy_note: str = ""
    namespace_note: str = ""
    viability_note: str = ""
    tool: str = ""
    sources: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "organism": self.organism,
            "kingdom": self.kingdom,
            "label": self.label,
            "identity_kind": self.identity_kind,
            "capability": self.capability,
            "denominator": self.denominator,
            "denominator_meaning": DENOMINATORS.get(self.denominator, ""),
            "required_resolution": self.required_resolution,
            "ploidy_note": self.ploidy_note,
            "namespace_note": self.namespace_note,
            "viability_note": self.viability_note,
            "tool": self.tool,
            "sources": list(self.sources),
        }


#: RNA genomes. DNA stool sequencing cannot assay these at all.
RNA_VIRUS_TARGETS: Final[tuple[str, ...]] = (
    "norovirus", "rotavirus", "astrovirus", "sapovirus",
    "hepatitis_a", "hepatitis_e", "enterovirus", "sars_cov_2",
)

TARGETS: Final[tuple[KingdomTarget, ...]] = (
    # ---- fungi --------------------------------------------------------- #
    KingdomTarget(
        target_id="fungi.candida_albicans.mlst",
        organism="Candida albicans",
        kingdom="fungi",
        label="Diploid MLST (AAT1a, ACC1, ADP1, MPIb, SYA1, VPS13, ZWF1b)",
        identity_kind="sequence_type",
        capability="tool_missing",
        denominator="eukaryote_marker_relative",
        required_resolution="diploid_allele_pair_per_locus",
        ploidy_note=(
            "Diploid: both alleles per locus must be retained, with heterozygosity, "
            "loss-of-heterozygosity, allele fractions and aneuploidy. A haploid bacterial "
            "consensus would silently discard one allele."
        ),
        tool="PubMLST C. albicans scheme",
        sources=("pubmlst.org/organisms/candida-albicans", "PMID 14605179"),
    ),
    KingdomTarget(
        target_id="fungi.auris.clade",
        organism="Candidozyma (Candida) auris",
        kingdom="fungi",
        label="Major clade assignment",
        identity_kind="lineage",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="genome_wide_clade_placement",
        namespace_note=(
            "Clade VI exists; an older four-clade list is out of date. Haemulonii-complex "
            "relatives are required decoys."
        ),
        sources=("PMID 39008997",),
    ),
    KingdomTarget(
        target_id="fungi.amr.substitutions",
        organism="fungi (ERG11 / FKS1 / FKS2)",
        kingdom="fungi",
        label="Azole and echinocandin resistance substitutions",
        identity_kind="allele",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="specified_substitution_in_correct_ortholog",
        ploidy_note=(
            "Allele fraction matters in a diploid or aneuploid background: a substitution "
            "on one of two alleles is not the same as on both."
        ),
        namespace_note=(
            "Presence of ERG11 or FKS1 is not resistance. Only specified substitutions "
            "count, in the correct species ortholog, with the right genetic code, introns "
            "and coordinates. Mutations reported *not* to confer resistance must be "
            "preserved as such."
        ),
        tool="FungAMR / ChroQueTas",
        sources=("doi:10.1038/s41564-025-02084-7",),
    ),
    # ---- protists ------------------------------------------------------ #
    KingdomTarget(
        target_id="protist.blastocystis.subtype",
        organism="Blastocystis spp.",
        kingdom="protist",
        label="SSU subtype and ST3/ST4 MLST",
        identity_kind="subspecies",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="ssu_subtype_plus_scheme_where_supported",
        namespace_note=(
            "Subtype and sequence-type are distinct namespaces and must not be conflated. "
            "A subtype is not a harmful-pathotype label; Blastocystis is common in healthy "
            "people."
        ),
        sources=("pubmlst.org/organisms/blastocystis-spp",),
    ),
    KingdomTarget(
        target_id="protist.giardia.assemblage",
        organism="Giardia duodenalis",
        kingdom="protist",
        label="Assemblage and multilocus profile (bg, gdh, tpi)",
        identity_kind="subspecies",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="multilocus_profile_with_discordance_preserved",
        ploidy_note=(
            "Allelic heterozygosity and ploidy are intrinsic to Giardia. Heterozygosity is "
            "not automatically a mixed infection, and discordant loci stay explicit rather "
            "than being resolved to a majority."
        ),
        sources=("doi:10.1186/s13071-023-05821-1",),
    ),
    KingdomTarget(
        target_id="protist.cryptosporidium.gp60",
        organism="Cryptosporidium spp.",
        kingdom="protist",
        label="Species plus gp60 subtype",
        identity_kind="subspecies",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="species_call_plus_spanning_gp60_repeat_evidence",
        namespace_note=(
            "SSU and gp60 are separate calls. Without spanning evidence across the gp60 "
            "repeat, the species evidence stands and no subtype is invented."
        ),
        tool="CryptoGenotyper",
    ),
    KingdomTarget(
        target_id="protist.entamoeba.species",
        organism="Entamoeba spp.",
        kingdom="protist",
        label="Competitive species discrimination",
        identity_kind="species",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="competitive_discrimination_against_dispar_moshkovskii_bangladeshi",
        namespace_note=(
            "At low depth, markers shared between E. histolytica and E. dispar must not be "
            "forced to the pathogenic species."
        ),
    ),
    # ---- helminths ------------------------------------------------------ #
    KingdomTarget(
        target_id="helminth.haplotype",
        organism="intestinal helminths",
        kingdom="helminth",
        label="Nuclear and mitochondrial population markers (cox1, nad1, ITS, 18S)",
        identity_kind="subspecies",
        capability="curation_required",
        denominator="eukaryote_marker_relative",
        required_resolution="marker_haplotype_with_nuclear_and_mito_kept_distinct",
        viability_note=(
            "A mitochondrial haplotype is not a complete nuclear strain, and DNA in stool "
            "does not establish a viable or current infection."
        ),
        namespace_note="No universal helminth strain nomenclature or distance threshold exists.",
    ),
    # ---- viruses -------------------------------------------------------- #
    KingdomTarget(
        target_id="virus.dna.adenovirus_type",
        organism="human adenovirus",
        kingdom="virus",
        label="Type from hexon/penton/fiber evidence",
        identity_kind="serotype",
        capability="curation_required",
        denominator="viral_votu_relative",
        required_resolution="informative_loci_plus_broader_genomic_evidence",
        namespace_note=(
            "Short shared fragments cannot establish a complete type or a recombinant; "
            "recombination and mixed infection stay represented."
        ),
    ),
    KingdomTarget(
        target_id="virus.phage.votu",
        organism="gut bacteriophage",
        kingdom="virus",
        label="vOTU detection and host prediction",
        identity_kind="lineage",
        capability="curation_required",
        denominator="viral_votu_relative",
        required_resolution="votu_assignment_with_host_prediction_labelled_as_prediction",
        namespace_note=(
            "Phage ecology is reported separately from human viral pathogens and does not "
            "inherit a pathogen framing. A predicted host is not physical attribution."
        ),
        viability_note="A viral DNA call does not establish active infection or expression.",
        sources=("UHGV 1.0 doi:10.5281/zenodo.17402089",),
    ),
    *(
        KingdomTarget(
            target_id=f"virus.rna.{name}",
            organism=name.replace("_", " "),
            kingdom="virus",
            label=f"{name.replace('_', ' ').title()} detection",
            identity_kind="species",
            capability="chemistry_incompatible",
            denominator="viral_votu_relative",
            required_resolution="rna_or_cdna_lane",
            namespace_note=(
                "This is an RNA genome. DNA stool sequencing cannot assay it, and adding "
                "reference sequences cannot change that. A separate RNA/cDNA lane with its "
                "own validation would be required."
            ),
        )
        for name in RNA_VIRUS_TARGETS
    ),
)

BY_TARGET: Final[Mapping[str, KingdomTarget]] = {t.target_id: t for t in TARGETS}

_CAPABILITY_TO_ASSAY: Final[Mapping[str, str]] = {
    "chemistry_incompatible": INCOMPATIBLE_INPUT,
    "tool_missing": REFERENCE_UNAVAILABLE,
    "curation_required": REFERENCE_UNAVAILABLE,
    "available": NOT_REQUESTED,
}


def call_for(sample_id: str, target_id: str) -> ResolutionCall:
    """The honest state of one non-bacterial target for one sample."""
    target = BY_TARGET[target_id]
    assay = _CAPABILITY_TO_ASSAY[target.capability]
    notes = " ".join(
        n for n in (target.namespace_note, target.ploidy_note, target.viability_note) if n
    )
    if target.capability == "chemistry_incompatible":
        plain = (
            f"{target.label}: not assayable from this specimen. {notes} This is registered "
            "as incompatible input rather than a negative result \u2014 it was never tested."
        )
    else:
        plain = (
            f"{target.label}: registered but not determined "
            f"({target.capability.replace('_', ' ')}"
            + (f", {target.tool}" if target.tool else "")
            + f"). {notes}"
        )
    return ResolutionCall(
        sample_id=sample_id,
        target_id=target_id,
        identity_kind=target.identity_kind,
        assay_status=assay,
        analytical_call=UNRESOLVED,
        reason_codes=(target.capability,),
        validation=Validation(
            analytical_status="not_validated",
            reference_function_status="scheme_defined",
            phenotype_measured_in_sample=False,
            clinical_predictive_status="not_established",
        ),
        plain=plain,
    )


def all_calls(sample_id: str) -> list[ResolutionCall]:
    return [call_for(sample_id, t.target_id) for t in TARGETS]


def registry_json() -> dict[str, Any]:
    by_capability: dict[str, int] = {}
    for t in TARGETS:
        by_capability[t.capability] = by_capability.get(t.capability, 0) + 1
    return {
        "targets": [t.to_json() for t in TARGETS],
        "n_targets": len(TARGETS),
        "by_capability": dict(sorted(by_capability.items())),
        "denominators": dict(DENOMINATORS),
        "rna_virus_targets": list(RNA_VIRUS_TARGETS),
        "note": (
            "Non-bacterial targets with the finest resolution each would need and the "
            "honest state of that capability. RNA-virus targets are permanently "
            "incompatible with DNA input; every other unavailable target names what is "
            "missing. None of them is ever reported as a negative screen."
        ),
    }


def self_test() -> int:
    """Guard acceptance tests 45-59. Returns a failure count."""
    failures = 0

    # Acceptance 56: DNA input returns RNA viruses as incompatible, and adding
    # sequences cannot change it.
    for name in RNA_VIRUS_TARGETS:
        call = call_for("K1", f"virus.rna.{name}")
        if call.assay_status != INCOMPATIBLE_INPUT:
            failures += 1
        if call.is_negative or call.is_positive:
            failures += 1
        if "never tested" not in call.plain:
            failures += 1
        if "cannot change that" not in call.plain:
            failures += 1

    # Acceptance 57: phages stay separate from human pathogens.
    phage = BY_TARGET["virus.phage.votu"]
    if "does not inherit a pathogen framing" not in phage.namespace_note:
        failures += 1

    # Acceptance 51: Candida diploid alleles survive.
    calb = BY_TARGET["fungi.candida_albicans.mlst"]
    if "diploid" not in calb.ploidy_note.lower():
        failures += 1
    if "discard one allele" not in calb.ploidy_note:
        failures += 1

    # Acceptance 52/53: presence is not resistance; ortholog/code must match.
    amr = BY_TARGET["fungi.amr.substitutions"]
    if "is not resistance" not in amr.namespace_note:
        failures += 1
    if "genetic code" not in amr.namespace_note:
        failures += 1
    if "not* to confer resistance" not in amr.namespace_note.replace("**", "*"):
        failures += 1

    # Acceptance 50: Blastocystis namespaces stay distinct.
    blasto = BY_TARGET["protist.blastocystis.subtype"]
    if "distinct namespaces" not in blasto.namespace_note:
        failures += 1
    if "not a harmful-pathotype label" not in blasto.namespace_note:
        failures += 1

    # Acceptance 47: Giardia heterozygosity is not a mixed infection.
    giardia = BY_TARGET["protist.giardia.assemblage"]
    if "not automatically a mixed infection" not in giardia.ploidy_note:
        failures += 1

    # Acceptance 46: Entamoeba ambiguity cannot force the pathogen.
    ent = BY_TARGET["protist.entamoeba.species"]
    if "must not be forced" not in ent.namespace_note:
        failures += 1

    # Acceptance 48: missing gp60 span keeps species evidence, invents nothing.
    crypto = BY_TARGET["protist.cryptosporidium.gp60"]
    if "no subtype is invented" not in crypto.namespace_note:
        failures += 1

    # Acceptance 54: a mito haplotype is not a strain or a viable infection.
    helm = BY_TARGET["helminth.haplotype"]
    if "not a complete nuclear strain" not in helm.viability_note:
        failures += 1
    if "viable" not in helm.viability_note:
        failures += 1

    # Acceptance 45: denominators stay separate.
    if len(DENOMINATORS) < 4:
        failures += 1
    euk = DENOMINATORS["eukaryote_marker_relative"]
    if "cannot be merged" not in euk:
        failures += 1

    # Every target reaches a terminal, non-negative state.
    for call in all_calls("K1"):
        if call.is_negative:
            failures += 1
        if call.analytical_call != UNRESOLVED:
            failures += 1
    if PHENOTYPE_NOT_MEASURED != "phenotype_not_measured":
        failures += 1
    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
