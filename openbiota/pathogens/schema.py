"""Data contracts for the pathogen branch (BUILD_SPEC_v05.0 §2.2, §4).

Everything the pathogen branch stores or renders passes through the records
here. Three principles from the spec shape the design and are worth stating
because they explain choices that would otherwise look verbose:

* **Orthogonal fields, not one status.** `assay_eligibility`,
  `analysis_status`, `sequence_status`, `contamination_status`,
  `reference_status` and `validation_scope` are independent. Collapsing them
  is how a system ends up claiming "negative" for a target it never searched.
* **`null` means unknown, never false.** An unset metadata field cannot be
  read as an assurance. `AssayManifest` therefore uses `None` throughout and
  the eligibility rules treat `None` as blocking, not permissive.
* **Names change; identifiers do not.** `target_id` is a stable application
  key. Accepted names, taxids and aliases move underneath it via versioned
  mappings so a historical report keeps its identity.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.errors import OpenBiotaError

__all__ = [
    "AssayManifest",
    "CoverageLedger",
    "DeterminantResult",
    "GROUPS",
    "GROUP_ORDER",
    "GROUP_LABELS",
    "INTERPRETATION_CLASSES",
    "NUCLEIC_ACID_PROTOCOLS",
    "PathogenError",
    "PathogenResult",
    "REFERENCE_STATUSES",
    "ReferenceLock",
    "TargetRecord",
]


class PathogenError(OpenBiotaError):
    """A pathogen catalog, reference bundle or result record is invalid."""


# --------------------------------------------------------------------------- #
# Vocabularies
# --------------------------------------------------------------------------- #

#: Spec §4.1. Microsporidia appear once, under parasites, with fungal
#: relationship metadata — never counted twice (acceptance P082).
GROUPS: Final[frozenset[str]] = frozenset(
    {
        "bacteria",
        "protozoa",
        "helminths",
        "microsporidia",
        "fungi",
        "other_eukaryotes",
        "dna_viruses",
        "rna_viruses",
        "other_viruses",
        "virulence",
        "amr",
    }
)

#: Detailed section order, spec §13.2. Determinants and the coverage ledger
#: follow the organism groups.
GROUP_ORDER: Final[tuple[str, ...]] = (
    "bacteria",
    "protozoa",
    "helminths",
    "microsporidia",
    "fungi",
    "other_eukaryotes",
    "dna_viruses",
    "rna_viruses",
    "other_viruses",
    "virulence",
    "amr",
)

GROUP_LABELS: Final[dict[str, str]] = {
    "bacteria": "Bacterial pathogens",
    "protozoa": "Protozoa pathogens",
    "helminths": "Parasitic worms",
    "microsporidia": "Microsporidia pathogens",
    "fungi": "Fungal pathogens",
    "other_eukaryotes": "Other eukaryotes",
    "dna_viruses": "DNA viral pathogens",
    "rna_viruses": "RNA viral pathogens",
    "other_viruses": "Other and uncertain viral findings",
    "virulence": "Toxin and virulence determinants",
    "amr": "Antimicrobial-resistance determinants",
}

INTERPRETATION_CLASSES: Final[frozenset[str]] = frozenset(
    {
        "established_enteric",
        "toxin_or_pathotype_dependent",
        "conditional_opportunist",
        "rare_enteric",
        "uncertain_enteric_role",
        "extraintestinal_watch",
        "background_or_decoy",
    }
)

#: Spec §3.1/§4.2. Stored separately: a genus-level hit cannot satisfy
#: species-level coverage, and a marker cannot masquerade as a genome.
REFERENCE_STATUSES: Final[frozenset[str]] = frozenset(
    {
        "genome_supported",
        "marker_only",
        "organelle_only",
        "protein_only",
        "unresolved_taxonomy",
        "no_usable_reference",
        "license_unavailable",
        "not_in_this_bundle",
        # The exact assembly this row would install is already installed under
        # a more specific row, so it is searched through that row rather than
        # twice. Installing it twice is not a duplicate-storage problem: two
        # rows over one genome share every k-mer, the ownership pass masks
        # both, and each loses the sequence that identified it. Adding
        # genus-level rows this way silently cost *Giardia duodenalis* and
        # *Cryptosporidium parvum* their species-level calls.
        "covered_by_relative",
    }
)

#: Reference statuses from which a nucleotide sequence search can proceed.
USABLE_REFERENCE_STATUSES: Final[frozenset[str]] = frozenset(
    {"genome_supported", "marker_only", "organelle_only"}
)

#: Reference statuses that cap a result at `marker_signal` (spec §5.4 rule 5).
MARKER_REFERENCE_STATUSES: Final[frozenset[str]] = frozenset(
    {"marker_only", "organelle_only"}
)

NUCLEIC_ACID_PROTOCOLS: Final[frozenset[str]] = frozenset(
    {"DNA", "RNA_with_RT", "total_nucleic_acid_with_RT", "unknown"}
)

LIBRARY_SELECTIONS: Final[frozenset[str]] = frozenset(
    {
        "shotgun",
        "poly_a_selected_rna",
        "rrna_depleted_rna",
        "amplicon",
        "targeted_capture",
        "unknown",
    }
)

STOOL_ROLES: Final[frozenset[str]] = frozenset(
    {
        "intestinal_shedding",
        "variable_shedding",
        "carriage_common",
        "tissue_restricted",
        "environmental_or_dietary",
        "extraintestinal",
        "not_applicable",
    }
)

#: Stool roles for which a non-detection is not part of any "infections
#: excluded" denominator (spec §7 rule 6, §9.5).
NON_EXCLUDING_STOOL_ROLES: Final[frozenset[str]] = frozenset(
    {"tissue_restricted", "extraintestinal"}
)


# --------------------------------------------------------------------------- #
# Assay manifest (spec §2.2)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AssayManifest:
    """What the laboratory actually did to this specimen.

    Every optional field defaults to `None`, and `None` is read as *unknown*
    everywhere downstream. The FASTQ alphabet cannot distinguish genomic DNA
    from reverse-transcribed RNA, so `nucleic_acid_protocol` must come from
    the upstream protocol record rather than from inspecting reads.
    """

    schema_version: str = "pathogens.assay.v1"
    sample_id: str = ""
    specimen_type: str = "stool"
    collected_at: str | None = None
    participant_age_years: float | None = None
    nucleic_acid_protocol: str = "DNA"
    reverse_transcription: bool = False
    library_selection: str = "shotgun"
    extraction_protocol_id: str | None = None
    input_stool_mass_mg: float | None = None
    preservative: str | None = None
    mechanical_lysis: bool | None = None
    wetlab_batch_id: str | None = None
    library_batch_id: str | None = None
    sequencing_run_id: str | None = None
    read_layout: str = "paired"
    sequencer: str | None = None
    read_length_summary: str | None = None
    control_sample_ids: tuple[str, ...] = ()
    spike_in_id: str | None = None
    recent_antimicrobials: bool | None = None
    recent_probiotics: bool | None = None
    recent_live_oral_vaccine: bool | None = None
    symptoms: tuple[str, ...] = ()
    immunocompromise: bool | None = None
    raw_fastq_sha256: tuple[str, ...] = ()
    host_removal_bundle_id: str | None = None

    def __post_init__(self) -> None:
        if self.nucleic_acid_protocol not in NUCLEIC_ACID_PROTOCOLS:
            raise PathogenError(
                f"nucleic_acid_protocol {self.nucleic_acid_protocol!r} is not one of "
                + ", ".join(sorted(NUCLEIC_ACID_PROTOCOLS))
            )
        if self.library_selection not in LIBRARY_SELECTIONS:
            raise PathogenError(
                f"library_selection {self.library_selection!r} is not one of "
                + ", ".join(sorted(LIBRARY_SELECTIONS))
            )

    @property
    def is_shotgun(self) -> bool:
        """An amplicon or capture library cannot inherit shotgun coverage."""
        return self.library_selection == "shotgun"

    @property
    def searches_dna(self) -> bool:
        """Whether genomic-DNA targets can be searched at all."""
        return self.nucleic_acid_protocol in {"DNA", "total_nucleic_acid_with_RT"}

    @property
    def searches_rna(self) -> bool:
        """Whether RNA-genome targets could be present in this library."""
        return self.nucleic_acid_protocol in {
            "RNA_with_RT",
            "total_nucleic_acid_with_RT",
        }

    @property
    def capability_line(self) -> str:
        """The assay-capability sentence for the front of the report."""
        if self.nucleic_acid_protocol == "DNA":
            base = "Stool DNA sequence screen"
        elif self.nucleic_acid_protocol == "RNA_with_RT":
            base = "Stool RNA (reverse-transcribed) sequence screen"
        elif self.nucleic_acid_protocol == "total_nucleic_acid_with_RT":
            base = "Stool total-nucleic-acid sequence screen"
        else:
            base = "Stool sequence screen, laboratory protocol unresolved"
        if not self.is_shotgun:
            base += f" ({self.library_selection.replace('_', ' ')} library)"
        return base

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "sample_id": self.sample_id,
            "specimen_type": self.specimen_type,
            "collected_at": self.collected_at,
            "participant_age_years": self.participant_age_years,
            "nucleic_acid_protocol": self.nucleic_acid_protocol,
            "reverse_transcription": self.reverse_transcription,
            "library_selection": self.library_selection,
            "extraction_protocol_id": self.extraction_protocol_id,
            "input_stool_mass_mg": self.input_stool_mass_mg,
            "preservative": self.preservative,
            "mechanical_lysis": self.mechanical_lysis,
            "wetlab_batch_id": self.wetlab_batch_id,
            "library_batch_id": self.library_batch_id,
            "sequencing_run_id": self.sequencing_run_id,
            "read_layout": self.read_layout,
            "sequencer": self.sequencer,
            "read_length_summary": self.read_length_summary,
            "control_sample_ids": list(self.control_sample_ids),
            "spike_in_id": self.spike_in_id,
            "recent_antimicrobials": self.recent_antimicrobials,
            "recent_probiotics": self.recent_probiotics,
            "recent_live_oral_vaccine": self.recent_live_oral_vaccine,
            "symptoms": list(self.symptoms),
            "immunocompromise": self.immunocompromise,
            "raw_fastq_sha256": list(self.raw_fastq_sha256),
            "host_removal_bundle_id": self.host_removal_bundle_id,
            "capability_line": self.capability_line,
        }


# --------------------------------------------------------------------------- #
# Target record (spec §4.1)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TargetRecord:
    """One installed, resolvable screening target.

    A target with unresolved identifiers or no accessions can never receive
    `not_detected` or a species-supported status; `search_ready` is the gate
    that enforces it.
    """

    schema_version: str
    target_id: str
    display_name: str
    aliases: tuple[str, ...]
    group: str
    interpretation_class: str
    ncbi_taxids: tuple[int, ...]
    taxonomic_resolution: str
    allowed_nucleic_acids: tuple[str, ...]
    stool_role: str
    reference_status: str
    reference_accessions: tuple[str, ...]
    marker_accessions: tuple[str, ...]
    near_neighbor_target_ids: tuple[str, ...]
    sequence_call_profile_id: str
    validation_id: str | None
    clinical_validation_status: str
    clinical_source_urls: tuple[str, ...]
    report_template_id: str
    confirmation_options: tuple[str, ...]
    negative_limitation_codes: tuple[str, ...]
    # --- provenance and presentation ------------------------------------- #
    source_section: str = ""
    family: str = ""
    note: str = ""
    reference_gap_note: str | None = None
    parent_target_id: str | None = None
    member_target_ids: tuple[str, ...] = ()
    expansion_state: str = "installed"
    requires_determinants: bool = False
    # --- viral-only mandatory fields (spec §4.1) -------------------------- #
    host_category: str | None = None
    host_evidence_level: str | None = None
    molecule_type: str | None = None
    segmented: bool = False
    segment_manifest_id: str | None = None
    segment_count: int | None = None
    rna_profile_required: bool = False
    dna_only_statement: str | None = None
    vaccine_shedding_context: bool = False
    who_priority: bool = False
    not_all_alleles_are: str | None = None
    # --- how many usable bases the bundle actually installed -------------- #
    #: When this row's reference is held by a more specific row over the same
    #: assembly, the row that holds it. The reads were searched against that
    #: reference, so this target is answered at the rank the two share rather
    #: than not answered at all.
    covered_by_target_id: str | None = None
    reference_bases: int = 0
    informative_bases: int = 0

    def __post_init__(self) -> None:
        if self.group not in GROUPS:
            raise PathogenError(f"{self.target_id}: unknown group {self.group!r}")
        if self.interpretation_class not in INTERPRETATION_CLASSES:
            raise PathogenError(
                f"{self.target_id}: unknown interpretation_class "
                f"{self.interpretation_class!r}"
            )
        if self.reference_status not in REFERENCE_STATUSES:
            raise PathogenError(
                f"{self.target_id}: unknown reference_status {self.reference_status!r}"
            )
        if self.stool_role not in STOOL_ROLES:
            raise PathogenError(
                f"{self.target_id}: unknown stool_role {self.stool_role!r}"
            )
        bad = set(self.allowed_nucleic_acids) - {"DNA", "RNA"}
        if bad:
            raise PathogenError(
                f"{self.target_id}: allowed_nucleic_acids has {sorted(bad)}"
            )
        if not self.allowed_nucleic_acids:
            raise PathogenError(f"{self.target_id}: allowed_nucleic_acids is empty")
        if self.group in {"dna_viruses", "rna_viruses", "other_viruses"}:
            for name in ("host_category", "host_evidence_level", "molecule_type"):
                if getattr(self, name) is None:
                    raise PathogenError(
                        f"{self.target_id}: viral targets must declare {name}"
                    )
            if self.segmented and not self.segment_manifest_id:
                raise PathogenError(
                    f"{self.target_id}: segmented virus needs a segment_manifest_id"
                )

    @property
    def reference_usable(self) -> bool:
        """Whether the installed reference can support a nucleotide search."""
        return self.reference_status in USABLE_REFERENCE_STATUSES

    @property
    def marker_only(self) -> bool:
        """Whether the best installed route caps results at `marker_signal`."""
        return self.reference_status in MARKER_REFERENCE_STATUSES

    @property
    def search_ready(self) -> bool:
        """A target may only receive a negative or a species call when this holds.

        Requires resolved taxonomy *and* installed sequence. Spec §4.1
        compiler rule.
        """
        if not self.ncbi_taxids:
            return False
        if not self.reference_usable:
            return False
        return bool(self.reference_accessions or self.marker_accessions)

    @property
    def excludable_by_stool(self) -> bool:
        """Whether a non-detection here belongs in an "assessed" denominator.

        Tissue-restricted and extraintestinal organisms never do, however good
        the genomic reference is (spec §7 rule 6).
        """
        return self.stool_role not in NON_EXCLUDING_STOOL_ROLES

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "schema_version": self.schema_version,
            "target_id": self.target_id,
            "display_name": self.display_name,
            "aliases": list(self.aliases),
            "group": self.group,
            "family": self.family,
            "interpretation_class": self.interpretation_class,
            "ncbi_taxids": list(self.ncbi_taxids),
            "taxonomic_resolution": self.taxonomic_resolution,
            "allowed_nucleic_acids": list(self.allowed_nucleic_acids),
            "stool_role": self.stool_role,
            "reference_status": self.reference_status,
            "reference_accessions": list(self.reference_accessions),
            "marker_accessions": list(self.marker_accessions),
            "near_neighbor_target_ids": list(self.near_neighbor_target_ids),
            "sequence_call_profile_id": self.sequence_call_profile_id,
            "validation_id": self.validation_id,
            "clinical_validation_status": self.clinical_validation_status,
            "clinical_source_urls": list(self.clinical_source_urls),
            "report_template_id": self.report_template_id,
            "confirmation_options": list(self.confirmation_options),
            "negative_limitation_codes": list(self.negative_limitation_codes),
            "source_section": self.source_section,
            "note": self.note,
            "reference_gap_note": self.reference_gap_note,
            "parent_target_id": self.parent_target_id,
            "member_target_ids": list(self.member_target_ids),
            "expansion_state": self.expansion_state,
            "requires_determinants": self.requires_determinants,
            "reference_bases": self.reference_bases,
            "informative_bases": self.informative_bases,
            "search_ready": self.search_ready,
        }
        if self.group in {"dna_viruses", "rna_viruses", "other_viruses"}:
            out.update(
                {
                    "host_category": self.host_category,
                    "host_evidence_level": self.host_evidence_level,
                    "molecule_type": self.molecule_type,
                    "segmented": self.segmented,
                    "segment_manifest_id": self.segment_manifest_id,
                    "segment_count": self.segment_count,
                    "rna_profile_required": self.rna_profile_required,
                    "dna_only_statement": self.dna_only_statement,
                    "vaccine_shedding_context": self.vaccine_shedding_context,
                }
            )
        if self.who_priority:
            out["who_priority"] = True
        if self.not_all_alleles_are:
            out["not_all_alleles_are"] = self.not_all_alleles_are
        return out


# --------------------------------------------------------------------------- #
# Result objects (spec §4.3, §12.4)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PathogenResult:
    """One target's outcome for one sample.

    Quantitative fields that require independent calibration —
    `absolute_load`, `clinical_lod`, `infection_probability` — stay `None`.
    They are never filled from read fractions or literature values, because a
    normalised read count is not an organism count and a threshold from
    another assay is not this assay's limit of detection.
    """

    schema_version: str
    sample_id: str
    target_id: str
    group: str
    display_name: str
    interpretation_class: str
    assay_eligibility: str
    analysis_status: str
    reference_status: str
    sequence_status: str
    contamination_status: str
    resolution: str
    clinical_interpretation: str
    validation_scope: str
    calling_profile_id: str
    reference_bundle_id: str
    unique_supporting_fragments: int = 0
    ambiguous_fragments: int = 0
    informative_regions_supported: int = 0
    informative_bases_covered: int = 0
    reference_breadth_fraction: float | None = None
    informative_region_breadth_fraction: float | None = None
    median_alignment_identity: float | None = None
    normalized_fragments_per_million: float | None = None
    normalization_denominator: str | None = None
    absolute_load: None = None
    clinical_lod: None = None
    infection_probability: None = None
    linked_determinants: tuple[str, ...] = ()
    unlinked_determinants: tuple[str, ...] = ()
    #: What the reads were actually searched against for this row, in one
    #: word, so no part of the report has to infer it from a status code:
    #:
    #: ``own_rank``            searched and answerable at the row's own rank
    #: ``group_rank``          searched; its DNA is identical to a catalogue
    #:                         relative, so it answers only at the shared rank
    #:                         (the seven diarrhoeagenic *E. coli* pathotypes
    #:                         are one genomic species, separated by the toxin
    #:                         and virulence genes in the determinant screen)
    #: ``out_of_assay_scope``  an RNA genome, which a DNA library cannot contain
    #: ``reference_pending``   no public reference sequence installed yet
    #: For a `group_rank` row, the row whose reference carries its sequence.
    covered_by_target_id: str | None = None
    search_scope: str = "own_rank"
    reason_codes: tuple[str, ...] = ()
    evidence_artifact_ids: tuple[str, ...] = ()
    confirmation_options: tuple[str, ...] = ()
    # --- presentation ----------------------------------------------------- #
    display_status: str = ""
    display_qualifier: str | None = None
    plain_statement: str = ""
    #: Carried through from the target so the report can tell an ordinary
    #: resident or a food organism from a pathogen without reloading the
    #: catalogue. Presenting the two alike is a wrong result, not a style
    #: choice: it puts brewer's yeast in the same list as Shigella.
    stool_role: str = "intestinal_shedding"
    # --- species resolution and pathotype gating (openbiota.pathogens.resolve) #
    #: Whether the species name is supportable from this evidence.
    #: `resolved` | `resolved_by_marker` | `not_resolvable_within_complex` |
    #: `group_level_only`.
    species_resolution: str = "resolved"
    complex_id: str | None = None
    complex_label: str | None = None
    #: Species rows whose evidence this group row absorbed.
    absorbed_from: tuple[str, ...] = ()
    #: For an organism whose disease potential *is* a toxin or pathotype:
    #: `supported` | `not_detected` | `not_assessed` | `not_applicable`.
    pathotype_evidence: str = "not_applicable"
    pathotype_markers: tuple[str, ...] = ()
    #: The single gate the report's headline uses. False for carriage without
    #: the disease-causing gene, and for a species the data cannot name.
    counts_as_pathogen: bool = True
    #: Overrides the tier the report would infer from `interpretation_class`.
    report_tier: str | None = None
    #: Published healthy-carriage rate, so a row can answer "is this normal?".
    carriage_statement: str | None = None
    carriage_source: str | None = None
    #: Set when a stricter interpretation was applied to this row on request
    #: (``strict_cdiff``). Recorded so the reader knows which rule produced
    #: the call; never set by the default profile.
    strict_mode: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "sample_id": self.sample_id,
            "target_id": self.target_id,
            "group": self.group,
            "display_name": self.display_name,
            "interpretation_class": self.interpretation_class,
            "stool_role": self.stool_role,
            "assay_eligibility": self.assay_eligibility,
            "analysis_status": self.analysis_status,
            "reference_status": self.reference_status,
            "sequence_status": self.sequence_status,
            "contamination_status": self.contamination_status,
            "resolution": self.resolution,
            "clinical_interpretation": self.clinical_interpretation,
            "validation_scope": self.validation_scope,
            "calling_profile_id": self.calling_profile_id,
            "reference_bundle_id": self.reference_bundle_id,
            "unique_supporting_fragments": self.unique_supporting_fragments,
            "ambiguous_fragments": self.ambiguous_fragments,
            "informative_regions_supported": self.informative_regions_supported,
            "informative_bases_covered": self.informative_bases_covered,
            "reference_breadth_fraction": self.reference_breadth_fraction,
            "informative_region_breadth_fraction": (
                self.informative_region_breadth_fraction
            ),
            "median_alignment_identity": self.median_alignment_identity,
            "normalized_fragments_per_million": self.normalized_fragments_per_million,
            "normalization_denominator": self.normalization_denominator,
            "absolute_load": None,
            "clinical_lod": None,
            "infection_probability": None,
            "linked_determinants": list(self.linked_determinants),
            "unlinked_determinants": list(self.unlinked_determinants),
            "species_resolution": self.species_resolution,
            "complex_id": self.complex_id,
            "complex_label": self.complex_label,
            "absorbed_from": list(self.absorbed_from),
            "pathotype_evidence": self.pathotype_evidence,
            "pathotype_markers": list(self.pathotype_markers),
            "counts_as_pathogen": self.counts_as_pathogen,
            "report_tier": self.report_tier,
            "strict_mode": self.strict_mode,
            "carriage_statement": self.carriage_statement,
            "carriage_source": self.carriage_source,
            "reason_codes": list(self.reason_codes),
            "evidence_artifact_ids": list(self.evidence_artifact_ids),
            "confirmation_options": list(self.confirmation_options),
            "display_status": self.display_status,
            "display_qualifier": self.display_qualifier,
            "search_scope": self.search_scope,
            "covered_by_target_id": self.covered_by_target_id,
            "plain_statement": self.plain_statement,
        }


DETERMINANT_STATUSES: Final[frozenset[str]] = frozenset(
    {
        "supported_intact_sequence",
        "supported_partial_sequence",
        "candidate_homolog",
        "ambiguous_allele",
        "not_detected",
        "not_assessed",
    }
)

HOST_LINKAGE_STATES: Final[frozenset[str]] = frozenset(
    {
        "unlinked",
        "read_pair_supported",
        "contig_supported",
        "genome_supported",
        "ambiguous",
    }
)

#: Linkage strong enough to say "this organism carries this determinant"
#: (spec §6.1). Same-sample co-occurrence is explicitly not enough.
LINKED_STATES: Final[frozenset[str]] = frozenset(
    {"read_pair_supported", "contig_supported", "genome_supported"}
)


@dataclass(frozen=True)
class DeterminantResult:
    """One toxin, virulence or resistance determinant's outcome.

    The organism-linkage question is answered separately from the
    gene-presence question. An organism can be supported while its toxin
    genotype is `not_assessed`, and a toxin gene can be supported while its
    carrier is unresolved. Neither state may be rounded into the other.
    """

    schema_version: str
    sample_id: str
    determinant_id: str
    display_name: str
    gene_family: str
    kind: str
    determinant_status: str
    host_linkage: str
    assay_eligibility: str
    analysis_status: str
    reference_status: str
    amr_class: str | None = None
    aligned_identity: float | None = None
    covered_reference_fraction: float | None = None
    locus_completeness: str = "not_assessed"
    required_discriminating_positions_assessed: bool = False
    gene_function_evidence: str = "curated_reference_family"
    clinical_phenotype_inference: str = "not_inferred"
    supporting_fragments: int = 0
    associated_target_ids: tuple[str, ...] = ()
    linked_target_ids: tuple[str, ...] = ()
    contig_ids: tuple[str, ...] = ()
    reason_codes: tuple[str, ...] = ()
    plain_statement: str = ""
    validation_scope: str = "computational_research_rule"

    def __post_init__(self) -> None:
        if self.determinant_status not in DETERMINANT_STATUSES:
            raise PathogenError(
                f"{self.determinant_id}: unknown determinant_status "
                f"{self.determinant_status!r}"
            )
        if self.host_linkage not in HOST_LINKAGE_STATES:
            raise PathogenError(
                f"{self.determinant_id}: unknown host_linkage {self.host_linkage!r}"
            )

    @property
    def is_linked(self) -> bool:
        return self.host_linkage in LINKED_STATES

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "sample_id": self.sample_id,
            "determinant_id": self.determinant_id,
            "display_name": self.display_name,
            "gene_family": self.gene_family,
            "kind": self.kind,
            "amr_class": self.amr_class,
            "determinant_status": self.determinant_status,
            "host_linkage": self.host_linkage,
            "assay_eligibility": self.assay_eligibility,
            "analysis_status": self.analysis_status,
            "reference_status": self.reference_status,
            "aligned_identity": self.aligned_identity,
            "covered_reference_fraction": self.covered_reference_fraction,
            "locus_completeness": self.locus_completeness,
            "required_discriminating_positions_assessed": (
                self.required_discriminating_positions_assessed
            ),
            "gene_function_evidence": self.gene_function_evidence,
            "clinical_phenotype_inference": self.clinical_phenotype_inference,
            "supporting_fragments": self.supporting_fragments,
            "associated_target_ids": list(self.associated_target_ids),
            "linked_target_ids": list(self.linked_target_ids),
            "contig_ids": list(self.contig_ids),
            "reason_codes": list(self.reason_codes),
            "plain_statement": self.plain_statement,
            "validation_scope": self.validation_scope,
        }


# --------------------------------------------------------------------------- #
# Coverage ledger (spec §4.2, §13.1 item 5)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CoverageLedger:
    """What was assessed, what was not, and why — the counts shown up front.

    Organism targets, determinants and subtypes are counted separately and
    never summed into a single "organisms screened" figure. Every number the
    report prints comes from here rather than from a hardcoded total.
    """

    assessed: int
    not_assessed_assay: int
    not_assessed_reference: int
    not_assessed_route: int
    determinants_assessed: int
    determinants_not_assessed: int
    by_group: dict[str, dict[str, int]] = field(default_factory=dict)
    reference_gap_reasons: dict[str, int] = field(default_factory=dict)
    failed_routes: tuple[str, ...] = ()
    total_targets: int = 0
    #: Rows whose reference *is* installed and *was* aligned against, but whose
    #: every informative k-mer is shared with a catalogue relative, so the row
    #: can only be answered at the rank it shares. The seven diarrhoeagenic
    #: *E. coli* pathotypes are all one genomic species and are separated by
    #: the toxin and virulence genes in the determinant screen, not by species
    #: DNA; the complex-and-species pairs (Klebsiella, Enterobacter,
    #: Salmonella, Bacteroides) are two rows over one genome. Counting these as
    #: "not searched" understated the screen: the reads were searched, and for
    #: the pathotypes the answer comes from the gene screen.
    searched_group_level: int = 0
    #: Rows searched whose every k-mer is also in the host or decoy genomes.
    #: GRCh38 ships Epstein-Barr virus as a decoy contig, so an EBV read
    #: cannot be told from a read of the reader; the search happened and
    #: cannot return an answer either way.
    searched_masked_by_host: int = 0
    #: Rows a DNA library cannot contain. An RNA virus has no DNA genome to
    #: sequence, so no depth of DNA sequencing can answer for it. Naming this
    #: as a coverage failure implied a fixable gap; it is a property of the
    #: molecule.
    out_of_assay_scope: int = 0
    #: Rows with no public reference sequence yet, split by whether the gap is
    #: the organism's (nothing deposited) or the catalogue's (a curated set
    #: that needs its members listed before it can be built).
    reference_pending_no_public_sequence: int = 0
    reference_pending_needs_member_list: int = 0

    @property
    def not_assessed(self) -> int:
        return (
            self.not_assessed_assay
            + self.not_assessed_reference
            + self.not_assessed_route
        )

    @property
    def searched(self) -> int:
        """Rows the reads were actually searched against, at any rank."""
        return self.assessed + self.searched_group_level + self.searched_masked_by_host

    @property
    def answerable(self) -> int:
        """Rows this assay can answer at all: the honest denominator.

        A target excluded because an RNA genome cannot appear in a DNA library
        does not belong in the same fraction as one that was skipped.
        """
        return max(0, self.total_targets - self.out_of_assay_scope)

    def to_json(self) -> dict[str, Any]:
        return {
            "total_targets": self.total_targets,
            "assessed": self.assessed,
            "not_assessed": self.not_assessed,
            "not_assessed_assay": self.not_assessed_assay,
            "not_assessed_reference": self.not_assessed_reference,
            "not_assessed_route": self.not_assessed_route,
            "determinants_assessed": self.determinants_assessed,
            "determinants_not_assessed": self.determinants_not_assessed,
            "searched": self.searched,
            "answerable": self.answerable,
            "searched_group_level": self.searched_group_level,
            "searched_masked_by_host": self.searched_masked_by_host,
            "out_of_assay_scope": self.out_of_assay_scope,
            "reference_pending_no_public_sequence": self.reference_pending_no_public_sequence,
            "reference_pending_needs_member_list": self.reference_pending_needs_member_list,
            "by_group": self.by_group,
            "reference_gap_reasons": self.reference_gap_reasons,
            "failed_routes": list(self.failed_routes),
        }


@dataclass(frozen=True)
class ReferenceLock:
    """Immutable identity of the reference bundle a result was produced from.

    Old reports keep their old lock id. A database or threshold change
    produces a new lock and new results; it never silently relabels history.
    """

    bundle_id: str
    created_at: str
    catalog_version: str
    taxonomy_snapshot: str
    taxonomy_sha256: str
    software: tuple[dict[str, str], ...] = ()
    index_parameters: dict[str, Any] = field(default_factory=dict)
    validation_profile_ids: tuple[str, ...] = ()
    build_command_log_sha256: str = ""
    asset_count: int = 0
    total_reference_bases: int = 0
    target_coverage: dict[str, int] = field(default_factory=dict)

    @staticmethod
    def now() -> str:
        return dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat()

    def to_json(self) -> dict[str, Any]:
        return {
            "bundle_id": self.bundle_id,
            "created_at": self.created_at,
            "catalog_version": self.catalog_version,
            "taxonomy_snapshot": self.taxonomy_snapshot,
            "taxonomy_sha256": self.taxonomy_sha256,
            "software": [dict(s) for s in self.software],
            "index_parameters": dict(self.index_parameters),
            "validation_profile_ids": list(self.validation_profile_ids),
            "build_command_log_sha256": self.build_command_log_sha256,
            "asset_count": self.asset_count,
            "total_reference_bases": self.total_reference_bases,
            "target_coverage": dict(self.target_coverage),
        }
