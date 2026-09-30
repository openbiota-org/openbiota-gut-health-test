"""Input adapters: native `results.json` → typed observations (spec §4.3, §5.2).

Trust order is native typed measurements first, native summary second, PDF
observations last. This module implements the first two, which is what the
repository actually produces; a PDF-only path is declared unavailable in the
capability manifest rather than faked.

Two rules shape everything here:

* **Namespaces and denominators are preserved.** MetaPhlAn 3 species percent,
  MetaPhlAn 4 SGB percent, panel copies per 100 genomes and pathogen fragments
  per million are four different measurements. They are never added, averaged
  or silently converted (V7-010, V7-034).
* **Unknown stays unknown.** A missing measurement is null with a reason, not
  zero, not NaN, not −1 (V7-024).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, Literal

from openbiota.errors import OpenBiotaError

from .identity import MaterialUnit, artifact_id, content_id, make_material

FeatureKind = Literal["species", "guild", "gene_panel", "pathway", "ecology", "sgb", "profile", "determinant", "pathogen_target"]
MeasurementStatus = Literal[
    "quantified", "detected_unquantified", "not_detected", "unassessed", "failed_qc", "conflicting"
]

#: Reference lanes that must match across every compared material, or the
#: comparison is descriptive only (spec §4.3: never mix one upgraded donor
#: with older recipient abundances).
LANE_KEYS: Final = (
    "taxonomic engine",
    "taxonomic reference",
    "reference fingerprint",
    "panels reported",
)


class FmtInputError(OpenBiotaError):
    """No usable analytical input (CLI exit 3)."""


@dataclass(slots=True)
class Observation:
    """One measurement of one feature in one material (spec §5.2)."""

    observation_id: str
    material_id: str
    sample_id: str
    feature_kind: FeatureKind
    source_id: str
    source_namespace: str
    canonical_id: str | None
    mapping_status: str
    resolution: str
    value: float | None
    lower: float | None
    upper: float | None
    unit: str
    denominator_id: str
    status: MeasurementStatus
    assay_tool: str
    assay_version: str
    assay_database: str | None
    detection_calibration: str
    reference_bundle_id: str | None = None
    reference_comparison_type: str | None = None
    percentile: float | None = None
    reference_p25: float | None = None
    reference_median: float | None = None
    reference_p75: float | None = None
    reference_prevalence: float | None = None
    eligibility_status: str = "unresolved"
    evidence_level: str = "native_measurement"
    locator: str = "results.json"
    limitations: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "observation_id": self.observation_id,
            "material_id": self.material_id,
            "sample_id": self.sample_id,
            "feature": {
                "kind": self.feature_kind,
                "canonical_id": self.canonical_id,
                "source_id": self.source_id,
                "source_namespace": self.source_namespace,
                "mapping_status": self.mapping_status,
                "resolution": self.resolution,
            },
            "measurement": {
                "value": self.value,
                "lower": self.lower,
                "upper": self.upper,
                "unit": self.unit,
                "denominator_id": self.denominator_id,
                "status": self.status,
            },
            "assay": {
                "tool": self.assay_tool,
                "version": self.assay_version,
                "database_id": self.assay_database,
                "detection_calibration": self.detection_calibration,
            },
            "reference": {
                "bundle_id": self.reference_bundle_id,
                "comparison_type": self.reference_comparison_type,
                "percentile": self.percentile,
                "p25": self.reference_p25,
                "median": self.reference_median,
                "p75": self.reference_p75,
                "prevalence": self.reference_prevalence,
                "eligibility_status": self.eligibility_status,
            },
            "provenance": {
                "evidence_level": self.evidence_level,
                "locator": self.locator,
            },
            "limitations": list(self.limitations),
            "extra": dict(self.extra),
        }


@dataclass(slots=True)
class MaterialData:
    """Everything loaded for one material."""

    material: MaterialUnit
    results: dict[str, Any]
    observations: list[Observation]
    lane: dict[str, Any]
    by_key: dict[tuple[str, str], Observation] = field(default_factory=dict)
    #: Species detected by the MetaPhlAn 4 lane, as a second opinion on
    #: presence only. Kept separate from `observations`, which are the
    #: MetaPhlAn 3 lane the whole comparison is scored against.
    #:
    #: The two lanes are not interchangeable: MetaPhlAn 4 uses SGB species
    #: definitions and MetaPhlAn 3 uses NCBI taxonomy, and curatedMetagenomic-
    #: Data - the only reference cohort available - is MetaPhlAn 3. So every
    #: percentile has to stay on the older lane. But "is this organism here
    #: at all" needs no reference, and for that question two detectors are
    #: strictly better than one: MetaPhlAn 4 is four years newer and built
    #: from roughly a million genomes against MetaPhlAn 3's seventeen
    #: thousand.
    secondary_species: dict[str, float] = field(default_factory=dict)

    def observation(self, kind: str, source_id: str) -> Observation | None:
        return self.by_key.get((kind, source_id))

    def secondary_presence(self, source_id: str) -> float | None:
        """Abundance in the second lane, or None if it did not see it."""
        value = self.secondary_species.get(source_id)
        return value if value and value > 0 else None


def _results_path(root: Path) -> Path:
    if root.is_file():
        return root
    direct = root / "results.json"
    if direct.is_file():
        return direct
    nested = sorted(root.glob("*/results.json"))
    if len(nested) == 1:
        return nested[0]
    raise FmtInputError(
        f"no results.json under {root}. Point --recipient/--donor at a sample output "
        "directory produced by `openbiota run`."
    )


def load_material(
    *, person_id: str, role: str, root: Path
) -> MaterialData:
    """Load one material's native results and derive typed observations."""
    path = _results_path(Path(root))
    try:
        results = json.loads(path.read_text())
    except json.JSONDecodeError as exc:  # pragma: no cover - corrupt input
        raise FmtInputError(f"{path} is not valid JSON: {exc}") from exc
    sample_id = str(results.get("sample") or "")
    if not sample_id:
        raise FmtInputError(f"{path} has no sample ID")
    run = results.get("run") or {}
    material = make_material(
        person_id=person_id,
        sample_id=sample_id,
        role=role,  # type: ignore[arg-type]
        input_root=path.parent,
        run_facts=run,
        artifacts=[artifact_id(path)],
    )
    lane = {k: run.get(k) for k in LANE_KEYS}
    lane["functional_cohort"] = ((results.get("reference_comparison") or {}).get("cohort") or {}).get("study")
    lane["extended_database"] = (results.get("extended_catalogue") or {}).get("database")
    data = MaterialData(material=material, results=results, observations=[], lane=lane)
    _species(data)
    _guilds(data)
    _guild_member_species(data)
    _panels(data)
    _sgbs(data)
    _ecology(data)
    _profiles(data)
    data.by_key = {(o.feature_kind, o.source_id): o for o in data.observations}
    data.secondary_species = _secondary_species(results)
    return data


def _secondary_species(results: dict[str, Any]) -> dict[str, float]:
    """Named species from the MetaPhlAn 4 catalogue, summed per species.

    A species can appear as more than one SGB; presence is the total. Only
    named species are usable here, because a goal names an organism and an
    unnamed SGB cannot satisfy it.
    """
    out: dict[str, float] = {}
    for row in (results.get("extended_catalogue") or {}).get("sgbs") or []:
        name = row.get("species")
        if not name or row.get("unnamed"):
            continue
        out[str(name)] = out.get(str(name), 0.0) + float(row.get("percent") or 0.0)
    # The whole-genome lane answers the same presence question against a
    # far larger species universe. A named organism it found counts as
    # present under the name the goal uses; its GTDB name is tried both as
    # written and under the older name the goal may carry.
    from openbiota.inventory import former_name

    genome = results.get("genome_profile") or {}
    for hit in genome.get("hits") or []:
        if not hit.get("confident", True) or hit.get("placeholder"):
            continue
        name = str(hit.get("species") or "").replace(" ", "_")
        if not name:
            continue
        pct = float(hit.get("taxonomic_abundance") or 0.0)
        for key in (name, former_name(name) or ""):
            if key and key not in out:
                out[key] = pct
    # The merged inventory (spec 0.8.4): every organism the expanded lanes
    # *support* - two independent methods or competitive confirmation - is
    # present under its accepted name and its GTDB alias. Provisional calls
    # are not: one method at a marginal level does not establish presence
    # for a donor decision.
    for row in (results.get("organism_inventory") or {}).get("organisms") or []:
        if row.get("status") != "supported" or row.get("unnamed"):
            continue
        reading = float(row.get("percent") or 0.0) or float(row.get("secondary_percent") or 0.0)
        if reading <= 0:
            continue
        for key in (str(row.get("species") or ""), str(row.get("gtdb") or "").replace(" ", "_"),
                    str(row.get("formerly") or "")):
            if key and key not in out:
                out[key] = reading
    return out


def _add(data: MaterialData, **kw: Any) -> Observation:
    kw.setdefault("canonical_id", None)
    obs = Observation(
        observation_id=content_id(
            "obs", {"m": data.material.material_id, "k": kw["feature_kind"], "s": kw["source_id"]}
        ),
        material_id=data.material.material_id,
        sample_id=data.material.sample_id,
        **kw,
    )
    data.observations.append(obs)
    return obs


def _species(data: MaterialData) -> None:
    """MetaPhlAn 3 species lane: the lane the disease patterns score on."""
    sim = data.results.get("profile_similarity") or {}
    comm = sim.get("community") or {}
    engine = (sim.get("taxonomic_engine") or {})
    db = engine.get("database")
    version = str(engine.get("profiler_version") or "unknown")
    ref_n = comm.get("reference_n")
    ref_src = comm.get("reference_source")
    for row in comm.get("species") or []:
        pct = row.get("percent")
        trace = bool(row.get("trace"))
        detected = pct is not None and pct > 0
        limits = []
        if trace:
            limits.append("trace_call_near_detection_limit")
        if not detected:
            limits.append("uncalibrated_absence")
        _add(
            data,
            feature_kind="species",
            source_id=str(row["species"]),
            source_namespace=f"metaphlan3:{db}",
            mapping_status="exact",
            resolution="species",
            value=float(pct) if pct is not None else None,
            lower=None,
            upper=None,
            unit="relative_percent",
            denominator_id="classified_species_mp3",
            status="quantified" if detected else "not_detected",
            assay_tool="MetaPhlAn",
            assay_version=version,
            assay_database=db,
            detection_calibration="uncalibrated",
            reference_bundle_id=str(ref_src) if ref_src else None,
            reference_comparison_type="population",
            percentile=row.get("percentile"),
            reference_prevalence=row.get("cohort_prevalence"),
            eligibility_status="unmatched_population_reference",
            limitations=limits,
            extra={
                "genus": row.get("genus"),
                "in_catalogue": row.get("in_catalogue"),
                "trace": trace,
                "reference_n": ref_n,
                "groups": list(row.get("groups") or []),
            },
        )


def _guilds(data: MaterialData) -> None:
    """Curated microbial guilds: the only source of curated polarity we have."""
    comm = (data.results.get("profile_similarity") or {}).get("community") or {}
    for row in comm.get("groups") or []:
        _add(
            data,
            feature_kind="guild",
            source_id=str(row["group"]),
            source_namespace="openbiota:taxa_guild",
            mapping_status="aggregate",
            resolution="guild",
            value=row.get("percent"),
            lower=None,
            upper=None,
            unit="relative_percent",
            denominator_id="classified_species_mp3",
            status="quantified" if row.get("detected") else "not_detected",
            assay_tool="openbiota.taxongroups",
            assay_version="1",
            assay_database=None,
            detection_calibration="uncalibrated",
            reference_comparison_type="population",
            percentile=row.get("percentile"),
            reference_prevalence=row.get("cohort_prevalence"),
            limitations=["aggregate_of_curated_members"],
            extra={
                "label": row.get("label"),
                "higher_means": row.get("higher_means"),
                "category": row.get("category"),
                "n_members": row.get("n_members"),
                "n_detected": row.get("n_detected"),
                "members": [
                    {
                        "species": m.get("species"),
                        "percent": m.get("percent"),
                        "detected": m.get("detected"),
                        "percentile": m.get("percentile"),
                        "cohort_prevalence": m.get("cohort_prevalence"),
                        "in_catalogue": m.get("in_catalogue"),
                    }
                    for m in row.get("members") or []
                ],
            },
        )


def _guild_member_species(data: MaterialData) -> None:
    """Add the curated guild members the species table omits.

    The species table lists what was detected. A curated member missing from
    it was still looked for by the same profiler on the same database, and the
    guild row records that with `detected: false` — so it is an *assessed*
    non-detection, not an unassessed feature. Recording it keeps a donor from
    earning a wide uncertainty band for a species the assay did examine.
    """
    have = {o.source_id for o in data.observations if o.feature_kind == "species"}
    sim = data.results.get("profile_similarity") or {}
    engine = sim.get("taxonomic_engine") or {}
    comm = sim.get("community") or {}
    db = engine.get("database")
    for guild in [o for o in data.observations if o.feature_kind == "guild"]:
        for member in guild.extra.get("members") or []:
            name = str(member.get("species"))
            if not name or name in have:
                continue
            if not member.get("in_catalogue", True):
                continue
            have.add(name)
            # This observation is the MetaPhlAn 3 lane's. A member the merged
            # inventory marks detected only by the expanded lanes was not seen
            # by this lane: it is recorded here as this lane's non-detection
            # (its own reading lives in the inventory and reaches presence
            # goals through `secondary_presence`).
            mp3_detected = bool(member.get("detected")) and not member.get("expansion_only") \
                and float(member.get("percent") or 0.0) > 0
            _add(
                data,
                feature_kind="species",
                source_id=name,
                source_namespace=f"metaphlan3:{db}",
                mapping_status="exact",
                resolution="species",
                value=float(member.get("percent") or 0.0) if mp3_detected else 0.0,
                lower=None,
                upper=None,
                unit="relative_percent",
                denominator_id="classified_species_mp3",
                status="quantified" if mp3_detected else "not_detected",
                assay_tool="MetaPhlAn",
                assay_version=str(engine.get("profiler_version") or "unknown"),
                assay_database=db,
                detection_calibration="uncalibrated",
                reference_bundle_id=str(comm.get("reference_source") or None),
                reference_comparison_type="population",
                percentile=member.get("percentile"),
                reference_prevalence=member.get("cohort_prevalence"),
                eligibility_status="unmatched_population_reference",
                limitations=["uncalibrated_absence"] if not mp3_detected else [],
                extra={
                    "in_catalogue": member.get("in_catalogue"),
                    "trace": False,
                    "from_guild_member_row": guild.source_id,
                    "reference_n": comm.get("reference_n"),
                    "groups": [guild.source_id],
                },
            )


def _panels(data: MaterialData) -> None:
    """Gene-panel capacities with their own reference distribution.

    The panel lane is the one place with real quartiles, so it is the one place
    a graded supply threshold can be sourced rather than invented (§5.3).
    """
    ref = (data.results.get("reference_comparison") or {}).get("panels") or {}
    cohort = ((data.results.get("reference_comparison") or {}).get("cohort") or {})
    run = data.results.get("run") or {}
    panel_meta = {p["name"]: p for p in data.results.get("panels") or []}
    for name, row in ref.items():
        meta = panel_meta.get(name, {})
        _add(
            data,
            feature_kind="gene_panel",
            source_id=str(name),
            source_namespace="openbiota:panel/diamond",
            mapping_status="exact",
            resolution="gene_family_panel",
            value=row.get("value"),
            lower=None,
            upper=None,
            unit="copies_per_100_bacterial_genomes",
            denominator_id="rpob_genome_equivalents",
            status="quantified" if row.get("value") is not None else "unassessed",
            assay_tool="DIAMOND",
            assay_version=str(run.get("diamond") or "unknown"),
            assay_database=str(run.get("reference fingerprint") or None),
            detection_calibration="reference_distribution_only",
            reference_bundle_id=str(cohort.get("study") or None),
            reference_comparison_type="population",
            percentile=row.get("percentile"),
            reference_p25=row.get("cohort_p25"),
            reference_median=row.get("cohort_median"),
            reference_p75=row.get("cohort_p75"),
            eligibility_status="unmatched_population_reference",
            limitations=["genetic_capacity_not_metabolite_concentration"],
            extra={
                "higher_means": row.get("higher_means"),
                "status_label": row.get("status"),
                "confidence": row.get("confidence"),
                "metabolite": meta.get("metabolite"),
                "reference_n": cohort.get("n_samples"),
            },
        )


def _sgbs(data: MaterialData) -> None:
    """MetaPhlAn 4 SGB census: a separate namespace and denominator."""
    ext = data.results.get("extended_catalogue") or {}
    db = ext.get("database")
    for row in ext.get("sgbs") or []:
        _add(
            data,
            feature_kind="sgb",
            source_id=str(row.get("sgb")),
            source_namespace=f"metaphlan4:{db}",
            mapping_status="exact",
            resolution="species_level_genome_bin",
            value=row.get("percent"),
            lower=None,
            upper=None,
            unit="relative_percent",
            denominator_id="classified_sgb_mp4",
            status="quantified",
            assay_tool="MetaPhlAn",
            assay_version=str(ext.get("profiler_version") or "4"),
            assay_database=db,
            detection_calibration="uncalibrated",
            limitations=["separate_lane_not_comparable_to_mp3_percent"],
            extra={"species": row.get("species"), "genus": row.get("genus"), "unnamed": row.get("unnamed")},
        )


def _ecology(data: MaterialData) -> None:
    """Diversity and the community index: context, never a matching objective."""
    cp = data.results.get("community_profile") or {}
    ext = data.results.get("extended_catalogue") or {}
    sim = data.results.get("profile_similarity") or {}
    comm = sim.get("community") or {}
    for key, value, unit in (
        ("shannon_index", cp.get("shannon_index"), "index"),
        ("simpson_index", cp.get("simpson_index"), "index"),
        ("pielou_evenness", cp.get("pielou_evenness"), "index"),
        ("species_richness", comm.get("n_species_detected"), "count"),
        ("genus_richness", comm.get("n_genera_detected"), "count"),
        ("sgb_richness", ext.get("n_sgbs_detected"), "count"),
        ("unclassified_percent_mp4", ext.get("unclassified_percent"), "percent"),
    ):
        _add(
            data,
            feature_kind="ecology",
            source_id=key,
            source_namespace="openbiota:ecology",
            mapping_status="exact",
            resolution="community",
            value=float(value) if isinstance(value, (int, float)) else None,
            lower=None,
            upper=None,
            unit=unit,
            denominator_id="community",
            status="quantified" if isinstance(value, (int, float)) else "unassessed",
            assay_tool="openbiota",
            assay_version="1",
            assay_database=None,
            detection_calibration="not_applicable",
            limitations=["descriptor_no_universal_polarity"],
        )
    anchor = sim.get("dysbiosis_anchor") or {}
    if anchor:
        _add(
            data,
            feature_kind="ecology",
            source_id="gmwi2",
            source_namespace="openbiota:gmwi2",
            mapping_status="exact",
            resolution="community",
            value=anchor.get("score"),
            lower=None,
            upper=None,
            unit="index",
            denominator_id="community",
            status="quantified" if anchor.get("score") is not None else "unassessed",
            assay_tool="GMWI2",
            assay_version="published_coefficients",
            assay_database=anchor.get("database"),
            detection_calibration="not_applicable",
            limitations=["context_only_rule_C01", "shares_inputs_with_species_lane"],
            extra={"band": anchor.get("band"), "valid": anchor.get("valid")},
        )


def _profiles(data: MaterialData) -> None:
    """Disease-pattern percentiles: context only (rule C01)."""
    sim = data.results.get("profile_similarity") or {}
    for row in sim.get("ranked") or []:
        _add(
            data,
            feature_kind="profile",
            source_id=str(row.get("profile")),
            source_namespace="openbiota:disease_pattern",
            mapping_status="exact",
            resolution="community_pattern",
            value=row.get("percentile"),
            lower=None,
            upper=None,
            unit="percentile",
            denominator_id="matched_reference_group",
            status="quantified" if row.get("percentile") is not None else "unassessed",
            assay_tool="openbiota.similarity",
            assay_version="1",
            assay_database=None,
            detection_calibration="uncalibrated_pattern_concordance",
            reference_comparison_type="population",
            percentile=row.get("percentile"),
            limitations=["resemblance_not_diagnosis", "context_only_rule_C01"],
            extra={"label": row.get("label"), "status": row.get("status"),
                   "evidence_maturity": row.get("evidence_maturity")},
        )


def lane_compatibility(recipient: MaterialData, donor: MaterialData) -> dict[str, Any]:
    """Compare analytical lanes; incompatible layers become descriptive only."""
    mismatches = {
        key: {"recipient": recipient.lane.get(key), "donor": donor.lane.get(key)}
        for key in set(recipient.lane) | set(donor.lane)
        if recipient.lane.get(key) != donor.lane.get(key)
    }
    return {
        "compatible": not mismatches,
        "mismatched_keys": sorted(mismatches),
        "detail": mismatches,
        "scope": (
            "all compared lanes identical: species, guild, panel and SGB comparisons are numeric"
            if not mismatches
            else "lane mismatch: affected comparisons are descriptive only"
        ),
    }
