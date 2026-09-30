"""Orchestration of the profile-similarity stage.

Connects the two engines, the two reference cohorts, the profile definitions
and the dysbiosis anchor into one scored result per profile. This is the only
module that knows about all of them; `scoring` knows nothing about where its
inputs come from and `profiles` knows nothing about how they are measured.

The reference bundle is assembled per run because stratified matching depends
on the subject: a request for ``country=USA, age_band=45-60`` selects a
sub-cohort, or falls back to the full cohort and says so.
"""

from __future__ import annotations

import dataclasses
import json
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from openbiota.dysbiosis import DysbiosisAnchor, compute_anchor
from openbiota.engines.metaphlan import (
    PINNED_INDEX,
    HostFilterResult,
    MetaphlanTools,
    TaxonomicProfile,
    build_host_index,
    database_ready,
    filter_host,
    host_index_ready,
    install_database,
    run_metaphlan,
)
from openbiota.errors import OpenBiotaError
from openbiota.logging_util import Reporter
from openbiota.profiles import (
    Profile,
    ProfileSet,
)
from openbiota.refcohort import (
    CLR_FLOOR_PERCENT,
    ReferenceCohort,
    build_taxon_references,
    sample_clr,
)
from openbiota.scoring import (
    ProfileResult,
    ReferenceBundle,
    SampleMeasurements,
    ecological_metrics,
    score_profile,
    summarise_reference,
)
from openbiota.shape import (
    ConfounderLedger,
    ShapeAnalysis,
    annotate_competition,
    confounder_ledger,
    empirical_specificity,
    literature_uniqueness_prior,
    shape_analysis,
)
from openbiota.skin import ANTIBIOTIC_EXCLUSION_DAYS
from openbiota.tally import TallyResult


@dataclass(frozen=True, slots=True)
class Subject:
    """What is known about the person, for matching, strata and abstention."""

    country: str | None = None
    age: float | None = None
    sex: str | None = None
    illness_duration_years: float | None = None
    antibiotics_days_ago: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    #: Keys of ``metadata`` whose values were assumed, not supplied
    #: (see :mod:`openbiota.context`). Travels into the abstention rules.
    assumed_keys: frozenset[str] = frozenset()

    @property
    def age_band(self) -> str | None:
        from openbiota.refcohort import age_band

        return age_band(self.age)

    @property
    def all_metadata(self) -> dict[str, Any]:
        """Free-form metadata plus the structured fields, for matching rules."""
        return {
            **self.metadata,
            "country": self.country,
            "age": self.age,
            "sex": self.sex,
            "illness_duration_years": self.illness_duration_years,
            "antibiotics_days_ago": self.antibiotics_days_ago,
        }

    def to_json(self) -> dict[str, Any]:
        return {
            "country": self.country,
            "age": self.age,
            "age_band": self.age_band,
            "sex": self.sex,
            "illness_duration_years": self.illness_duration_years,
            "antibiotics_days_ago": self.antibiotics_days_ago,
            "provided_metadata": sorted(
                k for k, v in self.all_metadata.items()
                if v not in (None, "", False) and k not in self.assumed_keys
            ),
            "assumed_metadata": sorted(self.assumed_keys),
        }


# --------------------------------------------------------------------------- #
# taxonomic engine
# --------------------------------------------------------------------------- #


def run_taxonomic_engine(
    *,
    tools: MetaphlanTools,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    out_dir: Path,
    reporter: Reporter,
    threads: int,
    skip_host_filter: bool = False,
    known_total_pairs: int | None = None,
) -> TaxonomicProfile:
    """Host filter (unless skipped) then MetaPhlAn, both cached.

    ``known_total_pairs`` (from fastp or the input count) lets a probe-only
    host pass report the pair total it did not itself count."""
    db_dir = refs_dir / "metaphlan_db"
    host_dir = refs_dir / "host"
    work = out_dir / "taxonomy"

    if not database_ready(db_dir):
        install_database(tools, db_dir, reporter, threads=threads)

    host: HostFilterResult | None = None
    query_r1, query_r2 = r1, r2
    if skip_host_filter:
        reporter.record("    host filtering skipped by request")
    elif not host_index_ready(host_dir):
        reporter.warn(
            "host index not built; profiling without host filtering (build it with "
            "`openbiota build-host-index`). Stool host fractions are usually below 1%, so the "
            "effect on relative abundances is small, but the host fraction goes unreported."
        )
    else:
        host = filter_host(
            tools, r1=r1, r2=r2, host_dir=host_dir, work_dir=work,
            reporter=reporter, threads=threads, known_total_pairs=known_total_pairs,
        )
        query_r1, query_r2 = host.nonhost_r1, host.nonhost_r2

    profile = run_metaphlan(
        tools, r1=query_r1, r2=query_r2, db_dir=db_dir, work_dir=work,
        reporter=reporter, threads=threads, host=host,
    )
    if host is not None and host.upstream_removed and not host.total_pairs and profile.n_reads_processed:
        # Nothing upstream counted the pairs; MetaPhlAn did (reads, so halve
        # for a pair). Fills the depth gate rather than leaving it at zero.
        pairs = profile.n_reads_processed // (2 if r2 is not None else 1)
        profile.host = dataclasses.replace(host, total_pairs=pairs)
    return profile


def ensure_host_index(tools: MetaphlanTools, refs_dir: Path, reporter: Reporter, threads: int) -> None:
    build_host_index(tools, refs_dir / "host", reporter, threads=threads)


# --------------------------------------------------------------------------- #
# references
# --------------------------------------------------------------------------- #


def load_functional_samples(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return data if isinstance(data, list) else []


def genus_of(species: str) -> str:
    """MetaPhlAn species name -> genus (``Bacteroides_vulgatus`` -> ``Bacteroides``)."""
    return species.split("_", 1)[0]


def species_key(organism: str) -> str:
    """Panel organism name -> MetaPhlAn species key (``Roseburia intestinalis`` -> ``Roseburia_intestinalis``)."""
    return organism.strip().replace(" ", "_")


def _frame_mean_log(vector: Sequence[float], frame: Sequence[int]) -> float:
    import math
    import statistics

    logs = [math.log(v if v > 0 else CLR_FLOOR_PERCENT) for v in vector]
    denominator = [logs[i] for i in frame] if frame else logs
    return statistics.fmean(denominator) if denominator else 0.0


def aggregate_log_ratio(total_percent: float, frame_mean_log: float) -> float:
    """Log-ratio of a summed group against the sample's prevalent-core frame.

    The same transform the species CLR uses — fixed floor, same frame
    denominator — applied to a sum of species. Sample and reference are
    transformed identically so their values are comparable.
    """
    import math

    return math.log(total_percent if total_percent > 0 else CLR_FLOOR_PERCENT) - frame_mean_log


def cohort_group_series(
    cohort: ReferenceCohort, members: Sequence[str]
) -> tuple[list[float], list[bool]]:
    """Per-reference-sample aggregate log-ratio and detection for a species group."""
    rows = [cohort.index_of(m) for m in members]
    rows = [r for r in rows if r is not None]
    frame_means = cohort.frame_mean_log
    values: list[float] = []
    detected: list[bool] = []
    for j in range(cohort.n_samples):
        total = sum(cohort.abundance[r][j] for r in rows)
        values.append(aggregate_log_ratio(total, frame_means[j]))
        detected.append(total > 0)
    return values, detected


def cohort_genus_members(cohort: ReferenceCohort) -> dict[str, list[str]]:
    genera: dict[str, list[str]] = {}
    for taxon in cohort.taxa:
        genera.setdefault(genus_of(taxon), []).append(taxon)
    return genera


def carrier_map(panels: Sequence[Any]) -> dict[str, tuple[str, ...]]:
    """Panel name -> curated carrier species, from each panel's ``organisms``.

    Carriage is not production (spec 5.1): this is an independent line of
    evidence beside the gene count — organisms known to carry the pathway,
    counted one way, against the pathway's genes, counted another.
    """
    out: dict[str, tuple[str, ...]] = {}
    for panel in panels:
        organisms = getattr(panel, "organisms", ()) or ()
        if organisms:
            out[panel.name] = tuple(species_key(o) for o in organisms)
    return out


def build_reference_bundle(
    *,
    cohort: ReferenceCohort | None,
    functional_samples: Sequence[dict[str, Any]],
    subject: Subject,
    match_on: Sequence[str],
    profiles: Sequence[Profile],
    sample_profiler_family: str | None,
    carriers: Mapping[str, Sequence[str]] | None = None,
) -> tuple[ReferenceBundle, ReferenceCohort | None]:
    """Assemble per-engine reference distributions, matched where possible."""
    bundle = ReferenceBundle()
    matched: ReferenceCohort | None = None

    if cohort is not None:
        criteria = {
            "country": subject.country if "country" in match_on else None,
            "band": subject.age_band if "age_band" in match_on else None,
            "sex": subject.sex if "sex" in match_on else None,
        }
        matched, info = cohort.match(**criteria)
        bundle.match_info = info
        taxon_refs = build_taxon_references(matched)
        source_label = f"curatedMetagenomicData {cohort.snapshot} (n={matched.n_samples:,})"

        # Genus-level references for transported 16S features: the sum of
        # member species, log-ratio transformed like everything else.
        wanted_genera = {
            f.key for p in profiles for f in p.features_for_engine("metaphlan") if f.level == "genus"
        }
        genera = cohort_genus_members(matched)
        for genus in wanted_genera:
            members = genera.get(genus)
            if not members:
                continue
            series, detected = cohort_group_series(matched, members)
            bundle.genus[genus] = summarise_reference(
                genus, series, source_label,
                prevalence=sum(detected) / max(1, len(detected)),
                detected=detected,
            )

        # Carrier-abundance references: summed abundance of each panel's
        # curated carriers across the reference cohort.
        wanted_carriers = {f.key for p in profiles for f in p.features_for_engine("carrier")}
        for panel_name in wanted_carriers:
            members = list((carriers or {}).get(panel_name, ()))
            members = [m for m in members if matched.index_of(m) is not None]
            if not members:
                continue
            series, detected = cohort_group_series(matched, members)
            bundle.carrier[panel_name] = summarise_reference(
                panel_name, series, source_label,
                prevalence=sum(detected) / max(1, len(detected)),
                detected=detected,
            )

        wanted = {
            f.key for p in profiles for f in p.features_for_engine("metaphlan") if f.level != "genus"
        }
        for name in wanted:
            reference = taxon_refs.get(name)
            if reference is None:
                continue
            # Which reference samples actually carried the taxon, as opposed to
            # having its value imputed at the zero floor. This is what lets the
            # scoring engine treat absence as a category instead of a number.
            row = matched.index_of(name)
            detected = (
                [matched.abundance[row][j] > 0 for j in range(matched.n_samples)]
                if row is not None
                else None
            )
            bundle.taxonomic[name] = summarise_reference(
                name, reference.values,
                f"curatedMetagenomicData {cohort.snapshot} (n={matched.n_samples:,})",
                prevalence=reference.prevalence,
                detected=detected,
            )
        bundle.taxonomic_source = f"curatedMetagenomicData {cohort.snapshot}"
        bundle.taxonomic_n = matched.n_samples
        bundle.taxonomic_full_n = cohort.n_samples

        # Ecological metrics computed on the same reference samples.
        eco_series: dict[str, list[float]] = {}
        for j in range(matched.n_samples):
            column = {
                matched.taxa[i]: matched.abundance[i][j]
                for i in range(matched.n_taxa)
                if matched.abundance[i][j] > 0
            }
            for key, value in ecological_metrics(column).items():
                eco_series.setdefault(key, []).append(value)
        for key, series in eco_series.items():
            bundle.ecological[key] = summarise_reference(
                key, series, f"curatedMetagenomicData {cohort.snapshot} (n={matched.n_samples:,})"
            )
        bundle.ecological_source = bundle.taxonomic_source
        bundle.ecological_n = matched.n_samples

        if sample_profiler_family is not None:
            bundle.profiler_matches_reference = cohort.profiler.startswith(sample_profiler_family)

    if functional_samples:
        panels = {
            k for s in functional_samples for k, v in s.items()
            if isinstance(v, (int, float)) and k not in ("rpob_fragments",)
        }
        for panel in panels:
            series = [float(s[panel]) for s in functional_samples if panel in s]
            bundle.functional[panel] = summarise_reference(
                panel, series, f"DIAMOND reference cohort (n={len(series)})"
            )
        # Profile functional features may reference panels under other names.
        for profile in profiles:
            for feature in profile.features_for_engine("diamond"):
                alias = FUNCTIONAL_ALIASES.get(feature.key)
                if alias and alias in bundle.functional and feature.key not in bundle.functional:
                    bundle.functional[feature.key] = bundle.functional[alias]
        bundle.functional_source = "DIAMOND reference cohort"
        bundle.functional_n = len(functional_samples)

    return bundle, matched


#: Profile functional-feature names that map onto existing panel names.
FUNCTIONAL_ALIASES: dict[str, str] = {
    "pks_island": "pks",
    "bft_toxin": "bft",
}


# --------------------------------------------------------------------------- #
# measurements
# --------------------------------------------------------------------------- #


def rescale_to_classified(species: Mapping[str, float]) -> dict[str, float]:
    """Rescale species abundances to sum to 100%.

    MetaPhlAn run with unknown estimation reports species as a share of the
    *whole sample*, so they sum to 100 minus the unclassified fraction — about
    50% for a typical stool sample. The curatedMetagenomicData reference
    profiles carry no unknown row and sum to 100.

    That mismatch is not cosmetic. The CLR's zero-replacement floor is a fixed
    absolute constant, so halving every observed value moves the floor's
    position relative to the rest of the composition, and an *absent* taxon
    ends up scoring near the 80th percentile of a reference where it is also
    absent. Rescaling to the classified fraction puts both sides on one scale.
    It is the same normalisation the published GMWI2 implementation applies.
    """
    total = sum(v for v in species.values() if v > 0)
    if total <= 0:
        return dict(species)
    factor = 100.0 / total
    return {k: v * factor for k, v in species.items()}


def build_measurements(
    *,
    tally: TallyResult | None,
    taxonomy: TaxonomicProfile | None,
    cohort: ReferenceCohort | None,
    subject: Subject,
    carriers: Mapping[str, Sequence[str]] | None = None,
) -> SampleMeasurements:
    measurements = SampleMeasurements(
        antibiotics_days_ago=subject.antibiotics_days_ago,
        metadata=subject.all_metadata,
        assumed_keys=subject.assumed_keys,
    )
    if tally is not None:
        for panel_result in tally.panels:
            value = panel_result.copies_per_100_genomes
            if value is not None:
                measurements.functional[panel_result.panel.name] = float(value)
        for name, alias in FUNCTIONAL_ALIASES.items():
            if alias in measurements.functional:
                measurements.functional[name] = measurements.functional[alias]

    if taxonomy is not None:
        species = rescale_to_classified(taxonomy.species)
        measurements.taxonomic_abundance = dict(species)
        if cohort is not None:
            measurements.taxonomic_clr = sample_clr(species, cohort)
            # Genus sums and carrier sums use the cohort's frame so that the
            # sample and the reference are on one scale.
            vector = [max(0.0, species.get(taxon, 0.0)) for taxon in cohort.taxa]
            frame_mean = _frame_mean_log(vector, cohort.frame)
            genus_totals: dict[str, float] = {}
            for taxon, value in species.items():
                if value > 0:
                    genus_totals[genus_of(taxon)] = genus_totals.get(genus_of(taxon), 0.0) + value
            for genus in cohort_genus_members(cohort):
                total = genus_totals.get(genus, 0.0)
                measurements.genus_abundance[genus] = total
                measurements.genus_clr[genus] = aggregate_log_ratio(total, frame_mean)
            for panel_name, members in (carriers or {}).items():
                total = sum(species.get(m, 0.0) for m in members)
                measurements.carrier_abundance[panel_name] = total
                measurements.carrier_clr[panel_name] = aggregate_log_ratio(total, frame_mean)
        measurements.ecological = ecological_metrics(species)
        if taxonomy.host is not None:
            measurements.usable_nonhost_reads = taxonomy.host.nonhost_pairs * 2
        elif taxonomy.n_reads_processed:
            measurements.usable_nonhost_reads = taxonomy.n_reads_processed
    return measurements


# --------------------------------------------------------------------------- #
# the stage
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class SimilarityStage:
    """Everything the report needs from the profile-similarity stage."""

    results: list[ProfileResult]
    taxonomy: TaxonomicProfile | None
    anchor: DysbiosisAnchor | None
    references: ReferenceBundle
    matched_cohort: ReferenceCohort | None
    subject: Subject
    cohort_manifest: dict[str, Any]
    elapsed_s: float
    notes: list[str] = field(default_factory=list)
    shape: ShapeAnalysis | None = None
    ledger: ConfounderLedger | None = None
    measurements: SampleMeasurements | None = None
    #: Profiles present in the library but not scored (opt-in, not requested).
    withheld: list[Profile] = field(default_factory=list)
    #: Taxon-group readings and the full species table (openbiota.community).
    community: Any | None = None
    #: The urticaria (hives) directional pattern indices, keyed by panel id.
    #: A second summary of the same species measurements on the source
    #: studies' own fixed scale, with coverage and censoring bounds.
    urticaria: dict[str, Any] = field(default_factory=dict)
    #: Inflammatory-skin-disease research panels (spec v4.2), keyed by panel
    #: id, plus the registered tasks that deliberately produce no number.
    skin: dict[str, Any] = field(default_factory=dict)

    @property
    def ranked(self) -> list[ProfileResult]:
        """Every result, scored ones first by percentile, then abstained, then unscoreable."""
        def key(r: ProfileResult) -> tuple[int, float]:
            if r.reportable and r.combined_percentile is not None:
                return (0, -r.combined_percentile)
            if r.abstention.abstained:
                return (1, 0.0)
            return (2, 0.0)
        return sorted(self.results, key=key)

    def to_json(self) -> dict[str, Any]:
        return {
            "subject": self.subject.to_json(),
            "taxonomic_engine": None if self.taxonomy is None else self.taxonomy.to_json(),
            "dysbiosis_anchor": None if self.anchor is None else self.anchor.to_json(),
            "shape_analysis": None if self.shape is None else self.shape.to_json(),
            "confounder_ledger": None if self.ledger is None else self.ledger.to_json(),
            "community": None if self.community is None else self.community.to_json(),
            "urticaria_pattern_indices": {k: v.to_json() for k, v in self.urticaria.items()},
            "skin_pattern_indices": self.skin,
            "ranked": [
                {
                    "profile": r.profile.name,
                    "label": r.profile.label,
                    "status": r.status,
                    "percentile": (
                        None if r.combined_percentile is None else round(r.combined_percentile, 1)
                    ),
                    "pattern_concordance_percent": (
                        None if r.concordance is None else round(r.concordance, 1)
                    ),
                    "evidence_maturity": r.profile.evidence_maturity,
                }
                for r in self.ranked
            ],
            "withheld_opt_in_profiles": [p.name for p in self.withheld],
            "reference": {
                "taxonomic": {
                    "source": self.references.taxonomic_source,
                    "n": self.references.taxonomic_n,
                    "matching": self.references.match_info,
                    "profiler_matches_sample": self.references.profiler_matches_reference,
                    "manifest": self.cohort_manifest,
                },
                "functional": {
                    "source": self.references.functional_source,
                    "n": self.references.functional_n,
                },
            },
            "profiles": {r.profile.name: r.to_json() for r in self.results},
            "notes": self.notes,
            "elapsed_s": round(self.elapsed_s, 1),
            "what_this_is": (
                "Resemblance of this community to published, group-level disease-associated "
                "patterns. A percentile is the fraction of the matched reference set that is "
                "less concordant with the pattern than this sample. It is not a diagnosis and "
                "not a probability of disease, and the underlying associations largely do not "
                "distinguish one condition from general gut disturbance."
            ),
        }


def stratum_for(profile: Profile, subject: Subject) -> str | None:
    """Pick the stratum from the subject, or None when unknown.

    Criteria are of the form ``field < 4``, ``field > 10`` or ``field == female``.
    A numeric field that fits no declared stratum returns ``intermediate``;
    an unknown field value returns None so both strata are reported.
    """
    if not profile.strata or not profile.stratify_on:
        return None
    value = subject.all_metadata.get(profile.stratify_on)
    if value in (None, ""):
        return None
    numeric = isinstance(value, (int, float))
    for stratum in profile.strata:
        criterion = stratum.criterion.replace(" ", "")
        if "==" in criterion:
            if str(value).lower() == criterion.split("==")[1].lower():
                return stratum.name
        elif numeric and "<" in criterion and float(value) < float(criterion.split("<")[1]) or numeric and ">" in criterion and float(value) > float(criterion.split(">")[1]):
            return stratum.name
    return "intermediate" if numeric else None


def run_similarity_stage(
    *,
    profile_set: ProfileSet,
    profile_names: Sequence[str] | None,
    tally: TallyResult | None,
    taxonomy: TaxonomicProfile | None,
    cohort_path: Path,
    functional_samples_path: Path,
    subject: Subject,
    match_on: Sequence[str],
    reporter: Reporter,
    include_opt_in: bool = False,
    profile_validation: Mapping[str, Any] | None = None,
    taxon_groups: Any | None = None,
    usable_read_pairs: int | None = None,
    qc_passed: bool = True,
    mode: str = "research",
) -> SimilarityStage:
    started = time.monotonic()
    notes: list[str] = []
    profiles = profile_set.select(profile_names, include_opt_in=include_opt_in)
    withheld = [p for p in profile_set.profiles if p not in profiles and p.opt_in]
    if withheld:
        notes.append(
            f"{len(withheld)} opt-in profile(s) not scored ({', '.join(p.name for p in withheld)}); "
            "request with --include-opt-in."
        )
    # Carrier sets come from two places: each gene panel's curated organisms
    # (carriage beside the gene count) and the curated taxon groups (so a
    # profile can score "oral-origin organisms" as one quantity). Where both
    # name the same thing the taxon group wins: it is the one the report shows.
    carriers: dict[str, tuple[str, ...]] = (
        carrier_map([p.panel for p in tally.panels]) if tally is not None else {}
    )
    if taxon_groups is not None and getattr(taxon_groups, "groups", ()):
        from openbiota.community import group_carrier_map

        carriers.update(group_carrier_map(taxon_groups))

    cohort: ReferenceCohort | None = None
    manifest: dict[str, Any] = {}
    if cohort_path.is_file():
        cohort = ReferenceCohort.load(cohort_path)
        manifest = cohort.manifest()
        reporter.record(
            f"    taxonomic reference: {cohort.n_samples:,} samples, snapshot {cohort.snapshot}"
        )
    else:
        notes.append(
            f"No taxonomic reference cohort at {cohort_path}; taxonomic and ecological "
            "modules cannot be scored. Build one with `openbiota taxonomic-cohort`."
        )

    functional_samples = load_functional_samples(functional_samples_path)
    if not functional_samples:
        notes.append("No per-sample functional cohort values; functional modules unscored.")

    references, matched = build_reference_bundle(
        cohort=cohort,
        functional_samples=functional_samples,
        subject=subject,
        match_on=match_on,
        profiles=profiles,
        sample_profiler_family=None if taxonomy is None else taxonomy.family,
        carriers=carriers,
    )
    if references.match_info:
        reporter.record(f"    reference matching: {references.match_info.get('reason')}")

    measurements = build_measurements(
        tally=tally, taxonomy=taxonomy, cohort=cohort, subject=subject, carriers=carriers
    )

    anchor: DysbiosisAnchor | None = None
    if taxonomy is not None:
        anchor = compute_anchor(
            taxonomy.clades, unknown_percent=taxonomy.unknown_percent, database=taxonomy.index
        )
        reporter.record(f"    dysbiosis anchor GMWI2 {anchor.score:+.2f} ({anchor.band})")
    else:
        notes.append("Taxonomic engine did not run; the dysbiosis anchor is unavailable.")

    results: list[ProfileResult] = []
    for profile in profiles:
        stratum = stratum_for(profile, subject)
        specificity, _basis = empirical_specificity(profile.name, profile_validation)
        result = score_profile(
            profile=profile,
            measurements=measurements,
            references=references,
            dysbiosis=anchor,
            stratum=stratum,
            empirical_specificity=specificity,
            literature_uniqueness_prior=literature_uniqueness_prior(profile, profile_set),
            reference_manifest={
                "taxonomic_manifest_id": manifest.get("manifest_id"),
                "taxonomic_snapshot": manifest.get("snapshot"),
                "taxonomic_profiler": manifest.get("profiler"),
                "sample_profiler": None if taxonomy is None else f"{taxonomy.family} {taxonomy.profiler_version}",
                "sample_database": None if taxonomy is None else taxonomy.index,
                "profile_version": profile.version,
            },
        )
        results.append(result)
        tax = result.primary_module("metaphlan")
        summary = (
            f"ABSTAINED ({len(result.abstention.triggered)} reason(s))"
            if result.abstention.abstained
            else (
                f"{result.combined_percentile:.0f}th pct, {result.status}, "
                f"{profile.evidence_maturity}"
                if result.combined_percentile is not None
                else f"not computable ({', '.join(m.module.evidence_type for m in result.unbound_modules) or 'no measurable features'})"
            )
        )
        reporter.record(
            f"    {profile.name:<22} {summary}"
            + (f"  [taxa {tax.n_measured}/{len(tax.features)}]" if tax else "")
        )

    # The urticaria panels get a second summary on the source studies' fixed
    # 0-100 scale. Same species measurements, same reference participants; it
    # adds the coverage and censoring bounds that a cohort rank cannot carry.
    urticaria: dict[str, Any] = {}
    if taxonomy is not None:
        from openbiota.cuindex import urticaria_indices

        urticaria = urticaria_indices(measurements.taxonomic_abundance, matched)

    # The inflammatory-skin panels. Same species measurements again, but each
    # panel keeps its own source contrast and evidence tier, and the
    # population and exposure gates decide whether a number is issued at all.
    skin: dict[str, Any] = {}
    if taxonomy is not None:
        from openbiota.skin import skin_indices

        skin = skin_indices(
            measurements.taxonomic_abundance,
            matched,
            age_years=subject.age,
            body_site="stool",
            usable_read_pairs=usable_read_pairs,
            counts_are_pairs=True,
            qc_passed=qc_passed,
            antibiotics_within_90_days=(
                None if subject.antibiotics_days_ago is None
                else subject.antibiotics_days_ago <= ANTIBIOTIC_EXCLUSION_DAYS
            ),
            # The report's declared mode, not one inferred from what metadata
            # happens to be present: participant mode fails closed on unknown
            # age, and research mode says so in its reason codes.
            mode=mode,
        )
        scored = sum(1 for k, v in skin.items() if isinstance(v, dict) and v.get("index") is not None)
        reporter.record(
            f"    skin panels: {scored} of {sum(1 for v in skin.values() if isinstance(v, dict) and 'status' in v)} "
            f"scored; {len(skin.get('unavailable_tasks', {}))} registered without a number"
        )

    annotate_competition(results)
    shape = shape_analysis(results, anchor)
    reporter.record(f"    shape: {shape.verdict} — {shape.headline}")
    ledger = confounder_ledger(measurements.metadata, results, assumed_keys=measurements.assumed_keys)

    community = None
    if taxonomy is not None and taxon_groups is not None:
        from openbiota.community import build_community_overview

        community = build_community_overview(
            species_percent=taxonomy.species,
            cohort=matched,
            group_set=taxon_groups,
            source_label=references.taxonomic_source or "",
        )
        reporter.record(
            f"    community: {community.n_species_detected} species, "
            f"{community.n_genera_detected} genera, {len(community.groups)} groups"
        )

    return SimilarityStage(
        community=community,
        results=results,
        taxonomy=taxonomy,
        anchor=anchor,
        references=references,
        matched_cohort=matched,
        subject=subject,
        cohort_manifest=manifest,
        elapsed_s=time.monotonic() - started,
        notes=notes,
        shape=shape,
        ledger=ledger,
        measurements=measurements,
        withheld=withheld,
        urticaria=urticaria,
        skin=skin,
    )


def require_index_name(taxonomy: TaxonomicProfile | None) -> str:
    return PINNED_INDEX if taxonomy is None else taxonomy.index


__all__ = [
    "OpenBiotaError",
    "SimilarityStage",
    "Subject",
    "build_measurements",
    "build_reference_bundle",
    "ensure_host_index",
    "run_similarity_stage",
    "run_taxonomic_engine",
    "stratum_for",
]
