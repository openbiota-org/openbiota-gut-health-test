"""The 14 public dataset records (BF-D01 .. BF-D14) and their correct use.

An accession is not a licence to train on something. Each record here
carries the assay that produced it, the specimen, whether the raw data are
actually obtainable, and whether the labels mean what a naive reader would
assume. Three of these datasets are routinely misdescribed in ways that
would silently corrupt a classifier:

* BF-D01 is 16S with an endoscopic phenotype. It is not shotgun training
  data, and its biopsy bloom rule is not a stool threshold.
* BF-D10's sequencing is of biopsies; the stool component is qPCR on a
  subset. Raw reads carrying human sequence are not public.
* BF-D12 compares mixed-species biofilms with single-species biofilms.
  Neither arm is planktonic, so it cannot train a biofilm-versus-planktonic
  model however convenient that would be.

``training_eligible`` is therefore computed, not declared: a dataset earns
it only by being shotgun, openly accessible, and carrying labels that have
actually been mapped. Everything else can inform the registry and appear in
``datasets audit`` output, but cannot enter supervised training.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

#: How the data were generated. Only ``shotgun_dna`` is assay-compatible
#: with this pipeline's own measurements.
ASSAYS: Final = (
    "16s",
    "shotgun_dna",
    "shotgun_dna_and_16s",
    "metatranscriptome",
    "qpcr",
    "imaging",
    "mixed",
    "isolate_genomes",
)

#: Whether the raw data can actually be downloaded.
ACCESS_STATES: Final = (
    "open",
    "request_only",
    "partially_open",
    "unconfirmed",
    "not_deposited",
)

#: Whether the dataset's labels have been mapped to this pipeline's
#: phenotype vocabulary. Unmapped labels cannot enter training.
LABEL_STATES: Final = ("mapped", "partially_mapped", "unmapped", "not_applicable")


@dataclass(frozen=True, slots=True)
class Dataset:
    """One public dataset, with the boundary of what it may be used for."""

    id: str
    accessions: tuple[str, ...]
    source_ids: tuple[str, ...]
    assay: str
    specimen: str
    host_species: str
    description: str
    correct_use: str
    exclusion: str
    access_status: str
    label_mapping_status: str
    #: Set when a specific miscount or mislabel is known and must be blocked.
    guard: str = ""

    def __post_init__(self) -> None:
        if self.assay not in ASSAYS:
            raise ValueError(f"{self.id}: unknown assay {self.assay!r}")
        if self.access_status not in ACCESS_STATES:
            raise ValueError(f"{self.id}: unknown access_status {self.access_status!r}")
        if self.label_mapping_status not in LABEL_STATES:
            raise ValueError(f"{self.id}: unknown label state {self.label_mapping_status!r}")

    @property
    def training_eligible(self) -> bool:
        """Whether this dataset may enter supervised training.

        Computed rather than declared, so a dataset cannot be waved through
        by editing a boolean. It must be shotgun (assay-compatible with our
        own measurements), openly downloadable, human, and carry labels
        that have actually been mapped.
        """
        return (
            self.assay in {"shotgun_dna", "shotgun_dna_and_16s"}
            and self.access_status == "open"
            and self.label_mapping_status == "mapped"
            and self.host_species == "human"
        )

    @property
    def blocking_reason(self) -> str:
        """Why this dataset is not training-eligible, in plain words."""
        if self.training_eligible:
            return ""
        reasons = []
        if self.assay not in {"shotgun_dna", "shotgun_dna_and_16s"}:
            reasons.append(f"assay is {self.assay}, not shotgun DNA")
        if self.access_status != "open":
            reasons.append(f"access is {self.access_status.replace('_', ' ')}")
        if self.label_mapping_status != "mapped":
            reasons.append(f"labels are {self.label_mapping_status.replace('_', ' ')}")
        if self.host_species != "human":
            reasons.append(f"host is {self.host_species}")
        return "; ".join(reasons)


DATASETS: Final[dict[str, Dataset]] = {
    d.id: d
    for d in (
        Dataset(
            id="BF-D01",
            accessions=("PRJNA644520",),
            source_ids=("BF-S01",),
            assay="16s",
            specimen="stool_and_biopsy",
            host_species="human",
            description=(
                "Stool and biopsy 16S with endoscopic biofilm phenotype from the "
                "Baumgartner cohort."
            ),
            correct_use=(
                "Phenotype and community context. Stool collection included the "
                "first stool after bowel preparation began, which must be preserved "
                "as a covariate."
            ),
            exclusion=(
                "Not shotgun training reads. The reported >3% bloom rule concerned "
                "biopsy analyses and is not a universal stool threshold."
            ),
            access_status="open",
            label_mapping_status="partially_mapped",
            guard="bloom_rule_is_biopsy_only",
        ),
        Dataset(
            id="BF-D02",
            accessions=("PRJEB66343",),
            source_ids=("BF-S06",),
            assay="shotgun_dna",
            specimen="stool",
            host_species="human",
            description=(
                "Original BAP study stool metagenomes. The 49-person analysis "
                "described 35 IBS and 14 controls."
            ),
            correct_use="Stool shotgun context for BAP-family features.",
            exclusion="Run-to-sample mapping must be verified before any use.",
            access_status="open",
            label_mapping_status="unmapped",
        ),
        Dataset(
            id="BF-D03",
            accessions=("PRJNA834801",),
            source_ids=("BF-S06",),
            assay="shotgun_dna",
            specimen="stool",
            host_species="human",
            description=(
                "Parkinson's cohort, 490 PD and 234 controls in the cited reanalysis."
            ),
            correct_use="Disease-labelled stool shotgun context.",
            exclusion=(
                "Disease labels are not biofilm imaging labels. A PD label does not "
                "make a sample biofilm-positive."
            ),
            access_status="open",
            label_mapping_status="unmapped",
            guard="disease_label_is_not_biofilm_label",
        ),
        Dataset(
            id="BF-D04",
            accessions=("PRJNA1279062",),
            source_ids=("BF-S13",),
            assay="mixed",
            specimen="laboratory_community",
            host_species="human",
            description="16S, metagenomic and metatranscriptomic data in one project.",
            correct_use="Separate specimen, treatment and assay before importing anything.",
            exclusion=(
                "Ex-vivo treated-community RNA is not a drop-in stool-DNA classifier."
            ),
            access_status="open",
            label_mapping_status="unmapped",
        ),
        Dataset(
            id="BF-D05",
            accessions=(
                "PRJNA1048869",
                "PRJNA473315",
                "PRJNA473316",
                "PRJNA857654",
            ),
            source_ids=("BF-S09",),
            assay="isolate_genomes",
            specimen="isolate",
            host_species="not_applicable",
            description="Klebsiella isolate, reference and outbreak contexts.",
            correct_use="Evaluating annotation and strain specificity of our own panels.",
            exclusion=(
                "Not a healthy-stool reference and not a validated gut-biofilm "
                "outcome dataset."
            ),
            access_status="open",
            label_mapping_status="not_applicable",
        ),
        Dataset(
            id="BF-D06",
            accessions=("PRJNA1228183",),
            source_ids=(),
            assay="16s",
            specimen="human_tissue",
            host_species="human",
            description="CRC tissue 16S context from the 2025 tissue study.",
            correct_use="Tissue community context.",
            exclusion="Not whole-gut shotgun truth.",
            access_status="open",
            label_mapping_status="unmapped",
        ),
        Dataset(
            id="BF-D07",
            accessions=(),
            source_ids=("BF-S05",),
            assay="shotgun_dna",
            specimen="stool",
            host_species="human",
            description=(
                "The bile-acid-diarrhoea study's own gene-family mappings and "
                "abundance tables."
            ),
            correct_use=(
                "Exact gene-family mappings, imported from supplements before any "
                "production feature mapping is selected."
            ),
            exclusion=(
                "A reusable public FASTQ accession was not confirmed in this audit. "
                "One must not be invented, and its absence must not block the "
                "independent gene module."
            ),
            access_status="unconfirmed",
            label_mapping_status="partially_mapped",
            guard="raw_access_unconfirmed",
        ),
        Dataset(
            id="BF-D08",
            accessions=("U_Net_bacteria_detection",),
            source_ids=("BF-S01",),
            assay="imaging",
            specimen="biopsy",
            host_species="human",
            description="Baumgartner bacterial image-segmentation code.",
            correct_use=(
                "An optional imaging research adapter, only after validation on "
                "compatible images."
            ),
            exclusion="Not a FASTQ biofilm predictor.",
            access_status="open",
            label_mapping_status="not_applicable",
        ),
        Dataset(
            id="BF-D09",
            accessions=("PRJNA1076354",),
            source_ids=("BF-S19",),
            assay="16s",
            specimen="laboratory_community",
            host_species="mixed",
            description=(
                "Public 16S from experimental biofilms and the accompanying mouse "
                "work."
            ),
            correct_use=(
                "Community and phenotype research. Runs must be audited against the "
                "article and supplements."
            ),
            exclusion=(
                "Human-derived communities, faecal inputs and mouse samples stay "
                "separate, and technical replicates are not independent donors. Not "
                "all original human faecal sequences or clinical covariates are "
                "deposited; INSPIRE clinical metadata require a governed request."
            ),
            access_status="partially_open",
            label_mapping_status="partially_mapped",
            guard="replicates_are_not_donors",
        ),
        Dataset(
            id="BF-D10",
            accessions=(),
            source_ids=("BF-S20",),
            assay="shotgun_dna_and_16s",
            specimen="biopsy",
            host_species="human",
            description=(
                "Paired tissue architecture and tissue metagenomes, with subset stool "
                "qPCR. Anonymised supplementary tables are public."
            ),
            correct_use=(
                "Openly available tables may be used immediately, with any "
                "unavailable sample mapping kept explicit."
            ),
            exclusion=(
                "Raw reads carrying human sequence are not public; processed data are "
                "available on request. There is no open stool-WGS training cohort "
                "here, and the sequencing was of biopsies."
            ),
            access_status="request_only",
            label_mapping_status="partially_mapped",
            guard="biopsy_wgs_and_stool_qpcr_are_distinct",
        ),
        Dataset(
            id="BF-D11",
            accessions=(),
            source_ids=("BF-S21",),
            assay="shotgun_dna",
            specimen="stool",
            host_species="human",
            description="Cirrhosis stool functional profiles and outcomes.",
            correct_use=(
                "Source-specific model and feature replication from available tables."
            ),
            exclusion=(
                "VA metadata are restricted and no specific raw accession was "
                "verified. Patient-level outcome validation cannot be invented."
            ),
            access_status="request_only",
            label_mapping_status="partially_mapped",
        ),
        Dataset(
            id="BF-D12",
            accessions=("PRJNA1188810",),
            source_ids=("BF-S22",),
            assay="metatranscriptome",
            specimen="laboratory_community",
            host_species="not_applicable",
            description=(
                "RNA from experimental mixed-species and single-species communities."
            ),
            correct_use="Preserve assay, species and experimental condition.",
            exclusion=(
                "It compares multi-species with single-species biofilms, NOT biofilm "
                "with planktonic state, and does not validate a protective "
                "human-stool score."
            ),
            access_status="open",
            label_mapping_status="not_applicable",
            guard="both_arms_are_biofilms",
        ),
        Dataset(
            id="BF-D13",
            accessions=("10.5281/zenodo.16886812", "PV739486"),
            source_ids=("BF-S25",),
            assay="16s",
            specimen="isolate",
            host_species="not_applicable",
            description=(
                "Figure and analysis data on Zenodo, plus the PV739486 16S sequence."
            ),
            correct_use="Validate actual files, version and licence before use.",
            exclusion=(
                "PV739486 is a 16S sequence. It cannot support strain-resolved "
                "protective-gene detection or stand in for a missing whole genome."
            ),
            access_status="open",
            label_mapping_status="not_applicable",
            guard="pv739486_is_16s_only",
        ),
        Dataset(
            id="BF-D14",
            accessions=(),
            source_ids=("BF-S32",),
            assay="imaging",
            specimen="human_tissue",
            host_species="human",
            description="Tissue imaging and supplementary data from the Indian CRC study.",
            correct_use="Tissue imaging and context validation, grouped by patient.",
            exclusion=(
                "Not a public stool-WGS calibration set. Paired samples remain "
                "grouped by patient."
            ),
            access_status="request_only",
            label_mapping_status="not_applicable",
            guard="paired_samples_stay_grouped",
        ),
    )
}

#: Metadata fields every dataset record must be able to answer before its
#: data may inform anything. Reported by ``biofilm datasets audit``.
REQUIRED_METADATA_FIELDS: Final = (
    "host_species",
    "specimen",
    "assay",
    "participant_id",
    "sample_id",
    "collection_timing",
    "phenotype_method",
    "phenotype_definition",
    "body_site",
    "disease_context",
    "age",
    "medications",
    "access_status",
    "label_mapping_status",
)


def get(dataset_id: str) -> Dataset:
    try:
        return DATASETS[dataset_id]
    except KeyError:
        raise KeyError(
            f"unknown dataset {dataset_id!r}; known IDs are BF-D01..BF-D14"
        ) from None


def audit() -> dict[str, object]:
    """Assay, label and access compatibility for every BF-D source.

    This is what ``openbiota biofilm datasets audit`` prints. The headline
    number is deliberately the count of training-eligible datasets, which
    is currently zero: every one of the fourteen is blocked by assay,
    access or label mapping. That is the honest state, and showing it
    stops a future change from quietly training on biopsy 16S.
    """
    rows = []
    for ds in DATASETS.values():
        rows.append(
            {
                "id": ds.id,
                "accessions": list(ds.accessions),
                "assay": ds.assay,
                "specimen": ds.specimen,
                "host_species": ds.host_species,
                "access_status": ds.access_status,
                "label_mapping_status": ds.label_mapping_status,
                "training_eligible": ds.training_eligible,
                "blocking_reason": ds.blocking_reason,
                "correct_use": ds.correct_use,
                "exclusion": ds.exclusion,
                "guard": ds.guard,
            }
        )
    eligible = [r for r in rows if r["training_eligible"]]
    return {
        "n_datasets": len(rows),
        "n_training_eligible": len(eligible),
        "training_eligible_ids": [r["id"] for r in eligible],
        "required_metadata_fields": list(REQUIRED_METADATA_FIELDS),
        "datasets": rows,
        "note": (
            "No dataset is currently training-eligible. Each is blocked for a "
            "stated reason: wrong assay, request-only access, or labels that have "
            "not been mapped. They inform the registry and the evidence cards; "
            "none of them silently enters supervised training."
        ),
    }
