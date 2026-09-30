"""Assay eligibility: can this library answer this target's question at all?

Spec §2.2 table plus §12.3 `check_assay_eligibility`. This module exists so
that "we did not look" is structurally distinguishable from "we looked and
found nothing". It answers one narrow question per target — is this library
capable of carrying this target's molecule — and refuses to answer it
optimistically:

* An RNA-genome virus in a DNA library is **ineligible**, and its row reads
  "not assessed", never "negative". No amount of extra software can recover a
  molecule the extraction and library never sequenced.
* An **unknown** protocol yields `unknown`, which is not truthy anywhere.
  Findings stay explorable; assay-wide negative claims do not.
* An amplicon or capture library does not inherit shotgun coverage, because
  its primers determine what could possibly appear.
* The RNA route stays `unknown` until `research_rna_virus_v1` is installed and
  benchmarked. Shipping the plumbing is not the same as activating the claim.

Participant age deliberately plays no part in eligibility. Presence screening
for pathogen sequence is not an adult microbiome percentile, so the existing
adult-reference abstention rule must not suppress it in minors (spec §2.2,
acceptance P009). Age belongs to clinical interpretation, which happens later.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from openbiota.pathogens.schema import AssayManifest, TargetRecord

__all__ = ["EligibilityRecord", "check_assay_eligibility", "self_test"]


@dataclass(frozen=True)
class EligibilityRecord:
    """Why a target is or is not answerable by this library."""

    target_id: str
    eligibility: str
    reason_code: str
    statement: str
    #: True when evidence may still be *explored* and reported as an
    #: assay-inconsistent or candidate record, even though no supported or
    #: negative call is permitted.
    discovery_allowed: bool = False

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "eligibility": self.eligibility,
            "reason_code": self.reason_code,
            "statement": self.statement,
            "discovery_allowed": self.discovery_allowed,
        }


def check_assay_eligibility(
    assay: AssayManifest,
    target: TargetRecord,
    *,
    rna_profile_installed: bool = False,
) -> EligibilityRecord:
    """Decide whether `assay` can answer `target`.

    `rna_profile_installed` reflects whether `research_rna_virus_v1` has been
    installed *and* passed its own benchmark. Until then RNA libraries return
    `unknown` for RNA targets: the reads may exist, but no validated rule
    exists to adjudicate them, so neither a supported call nor a negative is
    available.
    """
    wants_dna = "DNA" in target.allowed_nucleic_acids
    wants_rna = "RNA" in target.allowed_nucleic_acids
    protocol = assay.nucleic_acid_protocol

    # A non-shotgun library's content is set by its primers or baits. It
    # cannot inherit the coverage of an unbiased library.
    if not assay.is_shotgun and assay.library_selection != "unknown":
        selection = assay.library_selection.replace("_", " ")
        if assay.library_selection in {"amplicon", "targeted_capture"}:
            return EligibilityRecord(
                target.target_id,
                "ineligible",
                "non_shotgun_library",
                f"Not assessed: this {selection} library does not carry the "
                "unbiased shotgun coverage this target's search rule assumes.",
                discovery_allowed=True,
            )

    # Unknown protocol: explore, but claim nothing either way.
    if protocol == "unknown":
        return EligibilityRecord(
            target.target_id,
            "unknown",
            "protocol_unresolved",
            "Not assessed: the laboratory nucleic-acid protocol for this "
            "sample is unresolved, so neither a supported finding nor a "
            "negative result can be claimed for this target.",
            discovery_allowed=True,
        )

    dna_available = assay.searches_dna
    rna_available = assay.searches_rna

    # The ordinary case: a DNA-genome target in a library that carries DNA.
    if wants_dna and dna_available:
        return EligibilityRecord(
            target.target_id,
            "eligible",
            "dna_route_available",
            "Assessed by the stool DNA sequence route.",
            discovery_allowed=True,
        )

    # RNA-genome target, RNA-capable library, but no benchmarked RNA rule.
    if wants_rna and rna_available:
        if rna_profile_installed:
            return EligibilityRecord(
                target.target_id,
                "eligible",
                "rna_route_available",
                "Assessed by the reverse-transcribed RNA sequence route.",
                discovery_allowed=True,
            )
        return EligibilityRecord(
            target.target_id,
            "unknown",
            "rna_profile_not_installed",
            "Not assessed: this library can carry RNA, but the separately "
            "benchmarked RNA-virus calling profile is not installed. "
            "Discovery candidates remain reviewable.",
            discovery_allowed=True,
        )

    # RNA-genome target in a DNA-only library. The defining not-assessed case.
    if wants_rna and not wants_dna and not rna_available:
        statement = target.dna_only_statement or (
            f"Not assessed: this library was prepared for DNA, while "
            f"{target.display_name} has an RNA genome."
        )
        return EligibilityRecord(
            target.target_id,
            "ineligible",
            "wrong_nucleic_acid",
            statement,
            discovery_allowed=True,
        )

    # DNA-genome target in an RNA-only library: transcripts are candidates,
    # not a genomic presence/absence screen.
    if wants_dna and not dna_available and rna_available:
        return EligibilityRecord(
            target.target_id,
            "ineligible",
            "rna_library_for_dna_target",
            "Not assessed: this RNA library can show transcript-associated "
            "candidates, but it is not a genomic-DNA presence or absence "
            "screen for this target.",
            discovery_allowed=True,
        )

    return EligibilityRecord(
        target.target_id,
        "unknown",
        "no_matching_route",
        "Not assessed: no installed route matches this target's molecule and "
        "this library's preparation.",
        discovery_allowed=True,
    )


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #


def _target(**kw: Any) -> TargetRecord:
    base: dict[str, Any] = {
        "schema_version": "pathogens.target.v1",
        "target_id": "t",
        "display_name": "Test organism",
        "aliases": (),
        "group": "bacteria",
        "interpretation_class": "established_enteric",
        "ncbi_taxids": (1,),
        "taxonomic_resolution": "species",
        "allowed_nucleic_acids": ("DNA",),
        "stool_role": "intestinal_shedding",
        "reference_status": "genome_supported",
        "reference_accessions": ("GCF_1.1",),
        "marker_accessions": (),
        "near_neighbor_target_ids": (),
        "sequence_call_profile_id": "research_dna_v1",
        "validation_id": None,
        "clinical_validation_status": "not_established_for_this_pipeline",
        "clinical_source_urls": (),
        "report_template_id": "enteric_bacterium",
        "confirmation_options": (),
        "negative_limitation_codes": (),
    }
    base.update(kw)
    return TargetRecord(**base)


def self_test() -> int:
    checks = 0
    dna_target = _target()
    rna_target = _target(
        target_id="rna",
        group="rna_viruses",
        allowed_nucleic_acids=("RNA",),
        host_category="human",
        host_evidence_level="established_human_pathogen",
        molecule_type="ssRNA_positive",
        rna_profile_required=True,
        dna_only_statement=(
            "Not assessed: this library was prepared for DNA, while norovirus "
            "has an RNA genome."
        ),
    )

    dna_assay = AssayManifest(sample_id="s", nucleic_acid_protocol="DNA")

    # P003: DNA-only library, RNA-only target => not assessed, not negative.
    rec = check_assay_eligibility(dna_assay, rna_target)
    assert rec.eligibility == "ineligible" and rec.reason_code == "wrong_nucleic_acid"
    assert "RNA genome" in rec.statement
    checks += 1

    # A DNA target in the same library is fine.
    assert check_assay_eligibility(dna_assay, dna_target).eligibility == "eligible"
    checks += 1

    # P004: unknown must not be truthy-eligible.
    unknown = AssayManifest(sample_id="s", nucleic_acid_protocol="unknown")
    rec = check_assay_eligibility(unknown, dna_target)
    assert rec.eligibility == "unknown" and rec.eligibility != "eligible"
    checks += 1

    # P005: RNA reads without the installed RNA policy.
    rna_assay = AssayManifest(
        sample_id="s",
        nucleic_acid_protocol="RNA_with_RT",
        reverse_transcription=True,
        library_selection="rrna_depleted_rna",
    )
    rec = check_assay_eligibility(rna_assay, rna_target)
    assert rec.eligibility == "unknown"
    assert rec.reason_code == "rna_profile_not_installed"
    assert rec.discovery_allowed
    checks += 1

    rec = check_assay_eligibility(rna_assay, rna_target, rna_profile_installed=True)
    assert rec.eligibility == "eligible"
    checks += 1

    # A DNA target in an RNA library is not a genomic screen.
    rec = check_assay_eligibility(rna_assay, dna_target)
    assert rec.eligibility == "ineligible"
    assert rec.reason_code == "rna_library_for_dna_target"
    checks += 1

    # P007: amplicon inherits nothing.
    amplicon = AssayManifest(
        sample_id="s", nucleic_acid_protocol="DNA", library_selection="amplicon"
    )
    rec = check_assay_eligibility(amplicon, dna_target)
    assert rec.eligibility == "ineligible" and rec.reason_code == "non_shotgun_library"
    checks += 1

    # Total nucleic acid with RT can serve both branches, RNA still gated.
    total = AssayManifest(
        sample_id="s",
        nucleic_acid_protocol="total_nucleic_acid_with_RT",
        reverse_transcription=True,
    )
    assert check_assay_eligibility(total, dna_target).eligibility == "eligible"
    assert check_assay_eligibility(total, rna_target).eligibility == "unknown"
    assert (
        check_assay_eligibility(
            total, rna_target, rna_profile_installed=True
        ).eligibility
        == "eligible"
    )
    checks += 3

    # P009: a minor's stool DNA is screened exactly like an adult's. Age is
    # not an input to this decision at all.
    minor = AssayManifest(
        sample_id="s", nucleic_acid_protocol="DNA", participant_age_years=8.0
    )
    assert check_assay_eligibility(minor, dna_target).eligibility == "eligible"
    checks += 1

    # Retroviruses declare both molecules: proviral DNA is searchable.
    retro = _target(
        target_id="retro",
        group="rna_viruses",
        allowed_nucleic_acids=("DNA", "RNA"),
        host_category="human",
        host_evidence_level="established_human_pathogen",
        molecule_type="ssRNA_RT",
    )
    assert check_assay_eligibility(dna_assay, retro).eligibility == "eligible"
    checks += 1

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} eligibility checks passed")
