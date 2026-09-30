"""Seed catalog compilation (BUILD_SPEC_v05.0 §7, §12.1).

The catalog YAML files under `pathogens/catalog/` are the machine-readable
form of spec sections 8–11. This module turns them into `TargetRecord`s by
joining each seed against the *frozen installed manifest* produced by
`pathogens refs build`.

The join is where the spec's central discipline lives. A seed is a request
("screen for Balantioides coli"); a target is an installed capability. When
the manifest has no usable sequence for a seed, the compiler does not drop the
row and it does not fabricate one — it emits a target whose
`reference_status` records exactly why, which makes `search_ready` false,
which makes every downstream status `not_assessed` with a reason. That chain
is deliberate: it is what stops an unavailable organism from silently
becoming a clean bill of health.

Compilation is deterministic and needs no network and no language model. Sample
runs consume the frozen manifest; only reference builds resolve anything.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.pathogens.schema import (
    GROUPS,
    INTERPRETATION_CLASSES,
    STOOL_ROLES,
    PathogenError,
    TargetRecord,
)

__all__ = [
    "CATALOG_DIR",
    "CallingProfile",
    "Catalog",
    "DeterminantSeed",
    "Seed",
    "load_calling_profile",
    "load_catalog",
]

CATALOG_DIR: Final[Path] = Path("pathogens/catalog")
PROFILE_DIR: Final[Path] = Path("pathogens/profiles")

_RESOLUTIONS: Final[frozenset[str]] = frozenset(
    {
        "species",
        "species_complex",
        "subspecies",
        "genus",
        "group",
        "genotype",
        "serotype",
        "serovar",
        "pathotype",
    }
)

_EXPANSION_STATES: Final[frozenset[str]] = frozenset(
    {"installed", "manual_review_required"}
)

_DETERMINANT_KINDS: Final[frozenset[str]] = frozenset(
    {"toxin", "pathotype_marker", "virulence_locus", "amr"}
)


# --------------------------------------------------------------------------- #
# Seeds
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Seed:
    """One requested screening target, before reference resolution.

    Mirrors the spec §12.1 compilation record. `resolution_state` and the
    accession lists are filled by reference builds, not by this loader, so a
    freshly parsed seed is honestly `pending_resolution`.
    """

    seed_id: str
    source_section: str
    requested_name: str
    display_name: str
    aliases: tuple[str, ...]
    group: str
    interpretation_class: str
    taxonomic_resolution: str
    allowed_nucleic_acids: tuple[str, ...]
    stool_role: str
    candidate_sources: tuple[str, ...]
    report_template_id: str
    confirmation_options: tuple[str, ...]
    negative_limitation_codes: tuple[str, ...]
    note: str
    clinical_source_urls: tuple[str, ...] = ()
    family: str = ""
    #: The taxon whose sequence represents this row, when the row itself is
    #: not a taxon. Pathotypes, complexes and viral groups need this: STEC is
    #: a trait of Escherichia coli, not a species, so its organism reference
    #: is E. coli's while its identity comes from determinants.
    resolution_name: str | None = None
    reference_gap_note: str | None = None
    parent_seed_id: str | None = None
    member_seed_ids: tuple[str, ...] = ()
    expansion_state: str = "installed"
    requires_determinants: bool = False
    #: Install this exact assembly instead of the taxon's best one. Needed when
    #: "best assembly for the taxon" is the wrong genome for the row's purpose:
    #: the commensal *E. coli* competitor must be a commensal strain, and the
    #: taxon's best assembly is the O157:H7 pathogen the pathotype rows use.
    pinned_accession: str | None = None
    # viral extras
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

    @property
    def negative_claim(self) -> str:
        """Spec §12.1: the only negative a seed can ever license."""
        return "not_detected_under_applied_rule_only"


@dataclass(frozen=True)
class DeterminantSeed:
    """One toxin, virulence or resistance determinant request (spec §8.2–8.3)."""

    determinant_id: str
    source_section: str
    display_name: str
    gene_family: str
    kind: str
    candidate_sources: tuple[str, ...]
    note: str
    amr_class: str | None = None
    associated_seed_ids: tuple[str, ...] = ()
    subunits: tuple[str, ...] = ()
    requires_host_linkage: bool = False
    linkage_absent_statement: str | None = None
    clinical_source_urls: tuple[str, ...] = ()
    not_all_alleles_are: str | None = None
    sufficient_alone: bool | None = None
    substitutes_for_major_toxins: bool | None = None
    organism_specific_model: bool = False


def _tup(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    return tuple(str(v) for v in value)


def _seed_from_mapping(raw: Mapping[str, Any], source: Path) -> Seed:
    where = f"{source.name}:{raw.get('seed_id', '<no seed_id>')}"
    for required in (
        "seed_id",
        "source_section",
        "requested_name",
        "display_name",
        "group",
        "interpretation_class",
        "taxonomic_resolution",
        "allowed_nucleic_acids",
        "stool_role",
        "candidate_sources",
        "report_template_id",
        "note",
    ):
        if required not in raw:
            raise PathogenError(f"{where}: missing required field {required!r}")

    group = str(raw["group"])
    if group not in GROUPS:
        raise PathogenError(f"{where}: unknown group {group!r}")
    if str(raw["interpretation_class"]) not in INTERPRETATION_CLASSES:
        raise PathogenError(
            f"{where}: unknown interpretation_class {raw['interpretation_class']!r}"
        )
    if str(raw["taxonomic_resolution"]) not in _RESOLUTIONS:
        raise PathogenError(
            f"{where}: unknown taxonomic_resolution {raw['taxonomic_resolution']!r}"
        )
    if str(raw["stool_role"]) not in STOOL_ROLES:
        raise PathogenError(f"{where}: unknown stool_role {raw['stool_role']!r}")
    nucleic = _tup(raw["allowed_nucleic_acids"])
    if not nucleic or set(nucleic) - {"DNA", "RNA"}:
        raise PathogenError(f"{where}: allowed_nucleic_acids is {list(nucleic)}")
    state = str(raw.get("expansion_state", "installed"))
    if state not in _EXPANSION_STATES:
        raise PathogenError(f"{where}: unknown expansion_state {state!r}")

    # A seed_id must be a namespaced, lowercase, ASCII identifier: it becomes
    # a stable application key and appears in historical reports.
    seed_id = str(raw["seed_id"])
    if seed_id != seed_id.lower() or not seed_id.isascii() or "." not in seed_id:
        raise PathogenError(
            f"{where}: seed_id must be lowercase ASCII of the form group.name"
        )

    # Viral rows carry mandatory host and molecule metadata (spec §4.1):
    # isolation from a human sample does not establish a human host.
    if group in {"dna_viruses", "rna_viruses", "other_viruses"}:
        for required in ("host_category", "host_evidence_level", "molecule_type"):
            if raw.get(required) is None:
                raise PathogenError(f"{where}: viral seed must declare {required!r}")

    # An RNA-genome seed must carry its own not-assessed wording and must
    # declare that it needs the separately benchmarked RNA profile.
    if list(nucleic) == ["RNA"]:
        if not raw.get("rna_profile_required"):
            raise PathogenError(f"{where}: RNA-only seed needs rna_profile_required")
        if not raw.get("dna_only_statement"):
            raise PathogenError(f"{where}: RNA-only seed needs dna_only_statement")
        if "wrong_nucleic_acid" not in _tup(raw.get("negative_limitation_codes")):
            raise PathogenError(
                f"{where}: RNA-only seed must list the wrong_nucleic_acid limitation"
            )

    segment_count = raw.get("segment_count")
    return Seed(
        seed_id=seed_id,
        source_section=str(raw["source_section"]),
        requested_name=str(raw["requested_name"]),
        display_name=str(raw["display_name"]),
        aliases=_tup(raw.get("aliases")),
        group=group,
        interpretation_class=str(raw["interpretation_class"]),
        taxonomic_resolution=str(raw["taxonomic_resolution"]),
        allowed_nucleic_acids=nucleic,
        stool_role=str(raw["stool_role"]),
        candidate_sources=_tup(raw["candidate_sources"]),
        report_template_id=str(raw["report_template_id"]),
        confirmation_options=_tup(raw.get("confirmation_options")),
        negative_limitation_codes=_tup(raw.get("negative_limitation_codes")),
        note=str(raw["note"]).strip(),
        clinical_source_urls=_tup(raw.get("clinical_source_urls")),
        family=str(raw.get("family", "")),
        resolution_name=(
            None if raw.get("resolution_name") is None
            else str(raw["resolution_name"]).strip()
        ),
        reference_gap_note=(
            None if raw.get("reference_gap_note") is None
            else str(raw["reference_gap_note"]).strip()
        ),
        parent_seed_id=(
            None if raw.get("parent_seed_id") is None else str(raw["parent_seed_id"])
        ),
        member_seed_ids=_tup(raw.get("member_seed_ids")),
        expansion_state=state,
        requires_determinants=bool(raw.get("requires_determinants", False)),
        pinned_accession=(
            str(raw["pinned_accession"]) if raw.get("pinned_accession") else None
        ),
        host_category=(
            None if raw.get("host_category") is None else str(raw["host_category"])
        ),
        host_evidence_level=(
            None if raw.get("host_evidence_level") is None
            else str(raw["host_evidence_level"])
        ),
        molecule_type=(
            None if raw.get("molecule_type") is None else str(raw["molecule_type"])
        ),
        segmented=bool(raw.get("segmented", False)),
        segment_manifest_id=(
            None if raw.get("segment_manifest_id") is None
            else str(raw["segment_manifest_id"])
        ),
        segment_count=None if segment_count is None else int(segment_count),
        rna_profile_required=bool(raw.get("rna_profile_required", False)),
        dna_only_statement=(
            None if raw.get("dna_only_statement") is None
            else str(raw["dna_only_statement"]).strip()
        ),
        vaccine_shedding_context=bool(raw.get("vaccine_shedding_context", False)),
        who_priority=bool(raw.get("who_priority", False)),
    )


def _determinant_from_mapping(
    raw: Mapping[str, Any], source: Path
) -> DeterminantSeed:
    where = f"{source.name}:{raw.get('determinant_id', '<none>')}"
    for required in (
        "determinant_id",
        "source_section",
        "display_name",
        "gene_family",
        "kind",
        "candidate_sources",
        "note",
    ):
        if required not in raw:
            raise PathogenError(f"{where}: missing required field {required!r}")
    kind = str(raw["kind"])
    if kind not in _DETERMINANT_KINDS:
        raise PathogenError(f"{where}: unknown kind {kind!r}")
    if kind == "amr" and not raw.get("amr_class"):
        raise PathogenError(f"{where}: an AMR determinant must declare amr_class")
    return DeterminantSeed(
        determinant_id=str(raw["determinant_id"]),
        source_section=str(raw["source_section"]),
        display_name=str(raw["display_name"]),
        gene_family=str(raw["gene_family"]),
        kind=kind,
        candidate_sources=_tup(raw["candidate_sources"]),
        note=str(raw["note"]).strip(),
        amr_class=None if raw.get("amr_class") is None else str(raw["amr_class"]),
        associated_seed_ids=_tup(raw.get("associated_seed_ids")),
        subunits=_tup(raw.get("subunits")),
        requires_host_linkage=bool(raw.get("requires_host_linkage", False)),
        linkage_absent_statement=(
            None if raw.get("linkage_absent_statement") is None
            else str(raw["linkage_absent_statement"]).strip()
        ),
        clinical_source_urls=_tup(raw.get("clinical_source_urls")),
        not_all_alleles_are=(
            None if raw.get("not_all_alleles_are") is None
            else str(raw["not_all_alleles_are"])
        ),
        sufficient_alone=raw.get("sufficient_alone"),
        substitutes_for_major_toxins=raw.get("substitutes_for_major_toxins"),
        organism_specific_model=bool(raw.get("organism_specific_model", False)),
    )


# --------------------------------------------------------------------------- #
# Catalog
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Catalog:
    """The whole installed seed catalog plus its content-addressed version."""

    seeds: tuple[Seed, ...]
    determinants: tuple[DeterminantSeed, ...]
    version: str
    sources: tuple[str, ...] = ()
    by_id: dict[str, Seed] = field(default_factory=dict, repr=False)

    def group_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for seed in self.seeds:
            counts[seed.group] = counts.get(seed.group, 0) + 1
        return counts

    def compile_targets(
        self,
        manifest: Mapping[str, Mapping[str, Any]] | None = None,
        *,
        default_profile_id: str = "research_dna_v1",
    ) -> tuple[TargetRecord, ...]:
        """Join seeds against a frozen installed manifest.

        `manifest` maps `seed_id` to the entry written by `refs build`:
        resolved taxids, accessions, reference status and usable base counts.
        A seed absent from the manifest becomes a target with
        `reference_status="not_in_this_bundle"` — visible in the coverage
        ledger, never eligible for a species call or a negative.
        """
        manifest = manifest or {}
        records: list[TargetRecord] = []
        for seed in self.seeds:
            entry = manifest.get(seed.seed_id, {})
            status = str(entry.get("reference_status", "not_in_this_bundle"))
            taxids = tuple(int(t) for t in entry.get("ncbi_taxids", ()) or ())
            gap = seed.reference_gap_note
            if not gap and status not in {"genome_supported", "marker_only",
                                          "organelle_only"}:
                gap = str(entry.get("reference_gap_reason", "") or "") or None
            records.append(
                TargetRecord(
                    schema_version="pathogens.target.v1",
                    target_id=seed.seed_id,
                    display_name=seed.display_name,
                    aliases=seed.aliases,
                    group=seed.group,
                    interpretation_class=seed.interpretation_class,
                    ncbi_taxids=taxids,
                    taxonomic_resolution=seed.taxonomic_resolution,
                    allowed_nucleic_acids=seed.allowed_nucleic_acids,
                    stool_role=seed.stool_role,
                    reference_status=status,
                    reference_accessions=_tup(entry.get("reference_accessions")),
                    marker_accessions=_tup(entry.get("marker_accessions")),
                    near_neighbor_target_ids=_tup(entry.get("near_neighbor_target_ids")),
                    covered_by_target_id=entry.get("covered_by_target_id") or None,
                    sequence_call_profile_id=str(
                        entry.get("sequence_call_profile_id", default_profile_id)
                    ),
                    validation_id=entry.get("validation_id"),
                    clinical_validation_status="not_established_for_this_pipeline",
                    clinical_source_urls=seed.clinical_source_urls,
                    report_template_id=seed.report_template_id,
                    confirmation_options=seed.confirmation_options,
                    negative_limitation_codes=seed.negative_limitation_codes,
                    source_section=seed.source_section,
                    family=seed.family,
                    note=seed.note,
                    reference_gap_note=gap,
                    parent_target_id=seed.parent_seed_id,
                    member_target_ids=seed.member_seed_ids,
                    expansion_state=seed.expansion_state,
                    requires_determinants=seed.requires_determinants,
                    host_category=seed.host_category,
                    host_evidence_level=seed.host_evidence_level,
                    molecule_type=seed.molecule_type,
                    segmented=seed.segmented,
                    segment_manifest_id=seed.segment_manifest_id,
                    segment_count=seed.segment_count,
                    rna_profile_required=seed.rna_profile_required,
                    dna_only_statement=seed.dna_only_statement,
                    vaccine_shedding_context=seed.vaccine_shedding_context,
                    who_priority=seed.who_priority,
                    reference_bases=int(entry.get("reference_bases", 0) or 0),
                    informative_bases=int(entry.get("informative_bases", 0) or 0),
                )
            )
        return tuple(records)


def _validate_catalog(
    seeds: Sequence[Seed], determinants: Sequence[DeterminantSeed]
) -> None:
    ids = [s.seed_id for s in seeds]
    duplicates = sorted({i for i in ids if ids.count(i) > 1})
    if duplicates:
        raise PathogenError(f"duplicate seed_ids: {', '.join(duplicates)}")
    known = set(ids)
    for seed in seeds:
        if seed.parent_seed_id and seed.parent_seed_id not in known:
            raise PathogenError(
                f"{seed.seed_id}: parent_seed_id {seed.parent_seed_id!r} is not a seed"
            )
        for member in seed.member_seed_ids:
            if member not in known:
                raise PathogenError(
                    f"{seed.seed_id}: member_seed_id {member!r} is not a seed"
                )

    dids = [d.determinant_id for d in determinants]
    dup = sorted({i for i in dids if dids.count(i) > 1})
    if dup:
        raise PathogenError(f"duplicate determinant_ids: {', '.join(dup)}")
    for det in determinants:
        for seed_id in det.associated_seed_ids:
            if seed_id not in known:
                raise PathogenError(
                    f"{det.determinant_id}: associated_seed_id {seed_id!r} "
                    "is not a seed"
                )

    # Entamoeba coli and Escherichia coli are unrelated names, not aliases
    # (spec §7 rule 7). A naive splitting or aliasing rule collapses them, so
    # assert the separation holds in the installed data.
    by_id = {s.seed_id: s for s in seeds}
    ent = by_id.get("protozoa.entamoeba_coli")
    if ent is not None:
        assert ent.group == "protozoa", "Entamoeba coli must be a protozoan"
        for alias in ent.aliases:
            if "escherichia" in alias.lower():
                raise PathogenError(
                    "protozoa.entamoeba_coli aliases Escherichia coli; these are "
                    "unrelated organisms"
                )


def load_catalog(directory: Path | str = CATALOG_DIR) -> Catalog:
    """Load and validate every catalog file in `directory`.

    The catalog version is a hash of the file bytes, so any edit to a seed
    produces a new version and therefore new, separately identified results.
    Served from a fingerprinted cache when no file has changed; the cache key
    covers the same bytes the version hash does, so the two cannot disagree.
    """
    from openbiota.fastcache import cached_load

    path = Path(directory)
    if not path.is_dir():
        raise PathogenError(f"pathogen catalog directory not found: {path}")
    return cached_load(
        "pathogen_catalog", sorted(path.glob("*.yaml")),
        lambda: _load_catalog_uncached(path),
        code_modules=("openbiota.pathogens.schema",),
    )


def _load_catalog_uncached(path: Path) -> Catalog:
    files = sorted(path.glob("*.yaml"))
    if not files:
        raise PathogenError(f"no catalog YAML files in {path}")

    seeds: list[Seed] = []
    determinants: list[DeterminantSeed] = []
    digest = hashlib.sha256()
    for file in files:
        raw_bytes = file.read_bytes()
        digest.update(file.name.encode())
        digest.update(raw_bytes)
        try:
            doc = yaml.safe_load(raw_bytes)
        except yaml.YAMLError as exc:
            raise PathogenError(f"{file.name}: unreadable YAML ({exc})") from exc
        if not isinstance(doc, Mapping):
            raise PathogenError(f"{file.name}: top level must be a mapping")
        for entry in doc.get("seeds") or ():
            if not isinstance(entry, Mapping):
                raise PathogenError(f"{file.name}: a seed entry is not a mapping")
            seeds.append(_seed_from_mapping(entry, file))
        for entry in doc.get("determinants") or ():
            if not isinstance(entry, Mapping):
                raise PathogenError(f"{file.name}: a determinant entry is not a mapping")
            determinants.append(_determinant_from_mapping(entry, file))

    if not seeds:
        raise PathogenError(f"{path}: catalog contains no seeds")
    _validate_catalog(seeds, determinants)
    return Catalog(
        seeds=tuple(seeds),
        determinants=tuple(determinants),
        version=digest.hexdigest()[:16],
        sources=tuple(f.name for f in files),
        by_id={s.seed_id: s for s in seeds},
    )


# --------------------------------------------------------------------------- #
# Calling profiles (spec §5.4)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CallingProfile:
    """A frozen, versioned detection policy.

    The clinical fields stay `None` by construction: these are engineering
    starting points to be benchmarked, and a profile file cannot assert a
    limit of detection or a sensitivity it has not measured.
    """

    profile_id: str
    kraken_confidence: float
    minimum_hit_groups: int
    quick_mode: bool
    preserve_unclassified: bool
    min_aligned_query_fraction: float
    min_identity: float
    min_aligned_bases: int
    min_target_level_mapq: int
    min_best_minus_other_taxon_alignment_score: int
    exclude_low_complexity_or_masked_only_support: bool
    genome_min_distinct_fragments: int
    genome_min_regions: int
    genome_min_informative_bases: int
    require_resolution_specific_regions: bool
    marker_min_distinct_fragments: int
    marker_min_windows: int
    marker_min_informative_bases: int
    absence_claim: str
    region_bin_bases: int = 10_000
    viral_region_policy: str = "three_nonoverlapping_thirds"
    clinical_lod: None = None
    clinical_sensitivity: None = None
    clinical_specificity: None = None

    def __post_init__(self) -> None:
        if self.quick_mode:
            raise PathogenError(
                f"{self.profile_id}: quick mode is forbidden in the production "
                "pathogen branch"
            )
        if not self.preserve_unclassified:
            raise PathogenError(
                f"{self.profile_id}: unclassified reads must be preserved"
            )

    def to_json(self) -> dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "candidate_generation": {
                "kraken_confidence": self.kraken_confidence,
                "minimum_hit_groups": self.minimum_hit_groups,
                "quick_mode": self.quick_mode,
                "preserve_unclassified": self.preserve_unclassified,
            },
            "informative_alignment": {
                "min_aligned_query_fraction": self.min_aligned_query_fraction,
                "min_identity": self.min_identity,
                "min_aligned_bases": self.min_aligned_bases,
                "min_target_level_mapq": self.min_target_level_mapq,
                "min_best_minus_other_taxon_alignment_score": (
                    self.min_best_minus_other_taxon_alignment_score
                ),
                "exclude_low_complexity_or_masked_only_support": (
                    self.exclude_low_complexity_or_masked_only_support
                ),
            },
            "genome_support": {
                "min_distinct_fragments": self.genome_min_distinct_fragments,
                "min_nonoverlapping_informative_regions": self.genome_min_regions,
                "min_informative_bases_covered": self.genome_min_informative_bases,
                "require_resolution_specific_regions": (
                    self.require_resolution_specific_regions
                ),
            },
            "marker_or_organelle_support": {
                "min_distinct_fragments": self.marker_min_distinct_fragments,
                "min_nonoverlapping_supported_windows": self.marker_min_windows,
                "min_informative_bases_covered": self.marker_min_informative_bases,
                "report_as": "marker_signal",
            },
            "region_policy": {
                "bin_bases": self.region_bin_bases,
                "viral": self.viral_region_policy,
            },
            "absence_claim": self.absence_claim,
            "clinical_lod": None,
            "clinical_sensitivity": None,
            "clinical_specificity": None,
        }


#: The spec §5.4 `research_dna_v1` values, in code so that the branch is
#: runnable without an external file, and loadable from YAML so that a new
#: policy version is a data change with its own identifier.
RESEARCH_DNA_V1: Final[CallingProfile] = CallingProfile(
    profile_id="research_dna_v1",
    kraken_confidence=0.0,
    minimum_hit_groups=2,
    quick_mode=False,
    preserve_unclassified=True,
    min_aligned_query_fraction=0.90,
    min_identity=0.95,
    min_aligned_bases=75,
    min_target_level_mapq=30,
    min_best_minus_other_taxon_alignment_score=10,
    exclude_low_complexity_or_masked_only_support=True,
    genome_min_distinct_fragments=5,
    genome_min_regions=3,
    genome_min_informative_bases=300,
    require_resolution_specific_regions=True,
    marker_min_distinct_fragments=5,
    marker_min_windows=3,
    marker_min_informative_bases=300,
    absence_claim="no_supported_sequence_under_this_rule",
)


def load_calling_profile(
    profile_id: str = "research_dna_v1", directory: Path | str = PROFILE_DIR
) -> CallingProfile:
    """Load a calling profile, falling back to the built-in research default.

    An RNA profile id is accepted only when its file exists: the RNA route may
    not inherit `research_dna_v1` by filename (spec §11.3).
    """
    path = Path(directory) / f"{profile_id}.yaml"
    if not path.is_file():
        if profile_id == "research_dna_v1":
            return RESEARCH_DNA_V1
        raise PathogenError(
            f"calling profile {profile_id!r} is not installed at {path}; an RNA or "
            "alternative profile must be installed and benchmarked separately and "
            "cannot inherit research_dna_v1"
        )
    doc = yaml.safe_load(path.read_bytes())
    if not isinstance(doc, Mapping):
        raise PathogenError(f"{path}: calling profile must be a mapping")
    body = doc.get(profile_id, doc)
    cand = body.get("candidate_generation", {})
    align = body.get("informative_alignment", {})
    genome = body.get("genome_support", {})
    marker = body.get("marker_or_organelle_support", {})
    region = body.get("region_policy", {})
    for forbidden in ("clinical_lod", "clinical_sensitivity", "clinical_specificity"):
        if body.get(forbidden) is not None:
            raise PathogenError(
                f"{path}: {forbidden} must be null — a calling profile cannot "
                "assert clinical performance it has not measured"
            )
    return CallingProfile(
        profile_id=profile_id,
        kraken_confidence=float(cand.get("kraken_confidence", 0.0)),
        minimum_hit_groups=int(cand.get("minimum_hit_groups", 2)),
        quick_mode=bool(cand.get("quick_mode", False)),
        preserve_unclassified=bool(cand.get("preserve_unclassified", True)),
        min_aligned_query_fraction=float(align.get("min_aligned_query_fraction", 0.90)),
        min_identity=float(align.get("min_identity", 0.95)),
        min_aligned_bases=int(align.get("min_aligned_bases", 75)),
        min_target_level_mapq=int(align.get("min_target_level_mapq", 30)),
        min_best_minus_other_taxon_alignment_score=int(
            align.get("min_best_minus_other_taxon_alignment_score", 10)
        ),
        exclude_low_complexity_or_masked_only_support=bool(
            align.get("exclude_low_complexity_or_masked_only_support", True)
        ),
        genome_min_distinct_fragments=int(genome.get("min_distinct_fragments", 5)),
        genome_min_regions=int(
            genome.get("min_nonoverlapping_informative_regions", 3)
        ),
        genome_min_informative_bases=int(
            genome.get("min_informative_bases_covered", 300)
        ),
        require_resolution_specific_regions=bool(
            genome.get("require_resolution_specific_regions", True)
        ),
        marker_min_distinct_fragments=int(marker.get("min_distinct_fragments", 5)),
        marker_min_windows=int(marker.get("min_nonoverlapping_supported_windows", 3)),
        marker_min_informative_bases=int(
            marker.get("min_informative_bases_covered", 300)
        ),
        absence_claim=str(
            body.get("absence_claim", "no_supported_sequence_under_this_rule")
        ),
        region_bin_bases=int(region.get("bin_bases", 10_000)),
        viral_region_policy=str(
            region.get("viral", "three_nonoverlapping_thirds")
        ),
    )


def iter_group(catalog: Catalog, group: str) -> Iterable[Seed]:
    return (s for s in catalog.seeds if s.group == group)
