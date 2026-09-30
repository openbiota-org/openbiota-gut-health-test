"""Disease pattern profiles: declarative descriptions of published signatures.

A profile says "this published cohort reported these features moving in these
directions, measured this way, with this much evidence behind it." It contains
no scoring logic — that lives in `scoring.py` and has no knowledge of which
condition it is scoring.

Two schemas are accepted and produce the same in-memory objects:

* **Legacy (spec 2)** — three fixed modules ``taxonomic`` / ``functional`` /
  ``ecological`` plus an optional ``phenotype`` module, each feature carrying
  the four audit factors ``w q s c``.
* **Evidence-typed (spec 3)** — any number of named modules, each declaring an
  ``evidence_type`` (what was measured and by which engine), a ``claim_level``
  (what the measurement is allowed to mean), an ``assay_transport`` (whether
  the source measured something different from what this pipeline measures),
  and a ``study_group`` (which modules share participants and therefore do
  not count as independent replications).

Rules the loader enforces rather than documents:

* **Species, never genus, for conflicted taxa.** Aggregating a genus reverses
  direction in practice: genus *Faecalibacterium* has been reported higher in a
  cohort where species *F. prausnitzii* was lower. A genus on
  `CONFLICTED_TAXA` is refused unless the feature carries an explicit
  ``direction_conflict: material`` declaration.
* **Genus-level features only arrive by declared transport.** A 16S study
  reporting a genus can bind, but the module must say
  ``assay_transport: amplicon_taxon_to_shotgun`` so the report can label it.
* **Weights are audit trail, not multipliers.** Evidence quality, transport
  and confidence never multiply the score (that would pull weak evidence
  toward "typical" and make it read as biologically neutral rather than
  unknown). They travel alongside it in the confidence vector.
* **Same participants, one evidence group.** Two modules from one cohort's
  stool DNA share a ``study_group`` and fuse as one vote.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.errors import PanelError

# --------------------------------------------------------------------------- #
# vocabularies
# --------------------------------------------------------------------------- #

#: Genus-level taxa with no stable direction across cohorts. A profile may use
#: their constituent species, or the genus with an explicit conflict
#: declaration, but never the bare genus.
#:
#: Documented conflicts: *Streptococcus* and *Veillonella* are enriched in
#: Dutch and Hainan long-COVID cohorts but depleted in GI-predominant cases;
#: *Blautia* and *Lachnospiraceae* were enriched in a US long-COVID cohort
#: while treated as beneficial elsewhere; *Alistipes* direction varies by sex
#: and neurological phenotype; *Bacteroides* is case-enriched in severe
#: alopecia areata branches and lower in another; *Dorea* is higher at genus
#: level in one cohort while *D. longicatena* is control-associated in another.
CONFLICTED_TAXA: Final = frozenset(
    {
        "faecalibacterium",
        "streptococcus",
        "veillonella",
        "prevotella",
        "blautia",
        "lachnospiraceae",
        "alistipes",
        "bacteroides",
        "dorea",
    }
)

VALID_LEVELS: Final = frozenset({"species", "strain", "sgb", "genus"})
VALID_DIRECTIONS: Final = frozenset({"increased", "decreased"})
#: Legacy engine names. Evidence types map onto these.
VALID_ENGINES: Final = frozenset({
    "diamond", "metaphlan", "ecological", "carrier", "strainphlan", "unbound",
})
VALID_STATUS: Final = frozenset(
    {
        "research_only",
        "validation_control",
        "research_beta",
        "experimental",
        "production_candidate",
    }
)
VALID_MATURITY: Final = ("P0", "P1", "P2", "P3", "P4", "P5")
MATURITY_MEANING: Final = {
    "P0": "animal, case report, or mechanistic hypothesis only",
    "P1": "one small human cohort, amplicon-only, or inaccessible individual data",
    "P2": "usable human shotgun cohort, or multiple compatible human association cohorts",
    "P3": "independent human replication with coherent task and assay",
    "P4": "external population validation with a frozen model and proper participant split",
    "P5": "prospective intended-use validation with calibrated clinical performance",
}

VALID_TRANSPORTS: Final = frozenset(
    {"native", "amplicon_taxon_to_shotgun", "predicted_to_measured"}
)
TRANSPORT_LABEL: Final = {
    "native": "measured the same way as this sample",
    "amplicon_taxon_to_shotgun": "16S study; genus/species names carried to shotgun",
    "predicted_to_measured": "PICRUSt-inferred function tested against measured genes",
}

VALID_CLAIM_LEVELS: Final = frozenset(
    {
        "observed",
        "genomic_capacity",
        "carrier_proxy",
        "flux_predicted",
        "external_measured",
        "model_output",
        "association_only",
    }
)

#: Evidence type -> (engine in this pipeline, native claim level). An engine of
#: ``None`` means the type is registered but nothing here can measure it yet;
#: a module of that type loads, reports ``not_computable`` and says what would
#: bind it.
EVIDENCE_TYPES: Final[dict[str, tuple[str | None, str]]] = {
    "TAX_REL": ("metaphlan", "observed"),
    "TAX_PRES": ("metaphlan", "observed"),
    "GENE_ABUND": ("diamond", "genomic_capacity"),
    "CARRIER_ABUNDANCE": ("carrier", "carrier_proxy"),
    "ECO": ("ecological", "observed"),
    "STRAIN": (None, "observed"),
    "PANCNV": (None, "observed"),
    "SEQ_HMM": (None, "genomic_capacity"),
    "PATH_ABUND": (None, "genomic_capacity"),
    "PATH_COMPL": (None, "genomic_capacity"),
    "CAZY": (None, "genomic_capacity"),
    "AMR": (None, "genomic_capacity"),
    "FLUX": (None, "flux_predicted"),
    "METAB": (None, "external_measured"),
    "MODEL": (None, "model_output"),
}

#: Why a declared feature cannot be scored. These are read by the report, so
#: they must describe the real obstacle. STRAIN's entry changed once the marker
#: lane was built: a strain *engine* now exists and resolves population
#: fingerprints for every sufficiently covered organism (see the strain
#: section). What these particular features need is different — a validated,
#: target-specific reference panel for a named strain-level disease marker,
#: which for the entries in this library does not exist in the literature.
WOULD_BIND_IF: Final = {
    "STRAIN": (
        "a validated reference panel for this specific strain-level target. Strain "
        "resolution itself is available and reported separately: the marker lane "
        "fingerprints the dominant population of every organism with enough coverage. "
        "What is missing here is not the engine but a published, discriminating marker set "
        "for this named target"
    ),
    "PANCNV": "a pangenome copy-number engine (MIDAS2 / PanPhlAn) is added",
    "SEQ_HMM": "an HMMER stage with decoy models is added",
    "PATH_ABUND": "a HUMAnN 3 pathway/KO abundance stage is added",
    "PATH_COMPL": "assembly plus frozen reaction-completeness rules are added",
    "CAZY": "assembly plus dbCAN3 carbohydrate-active enzyme annotation is added",
    "AMR": "CARD / RGI resistome annotation is added",
    "FLUX": "an AGORA2 / MICOM flux stage with declared diet scenarios is added",
    "METAB": "a measured metabolomics file is supplied alongside the FASTQ",
    "MODEL": "the named compatibility container ships with a leakage-corrected reproduction",
}

MODULE_DIRECTION_PANEL: Final = "direction_panel"
MODULE_SENSITIVITY_PANEL: Final = "sensitivity_panel"
MODULE_MECHANISM_LANE: Final = "mechanism_lane"
MODULE_PHENOTYPE_PANEL: Final = "phenotype_panel"
VALID_MODULE_TYPES: Final = frozenset(
    {MODULE_DIRECTION_PANEL, MODULE_SENSITIVITY_PANEL, MODULE_MECHANISM_LANE, MODULE_PHENOTYPE_PANEL}
)

#: Legacy module names. Kept as constants because the report and the scoring
#: engine still look for them by name on legacy profiles.
MODULE_TAXONOMIC: Final = "taxonomic"
MODULE_FUNCTIONAL: Final = "functional"
MODULE_ECOLOGICAL: Final = "ecological"
MODULE_PHENOTYPE: Final = "phenotype"
COMBINED_MODULES: Final = (MODULE_TAXONOMIC, MODULE_FUNCTIONAL, MODULE_ECOLOGICAL)
VALID_MODULES: Final = frozenset({*COMBINED_MODULES, MODULE_PHENOTYPE})

LEGACY_EVIDENCE: Final = {
    MODULE_TAXONOMIC: "TAX_REL",
    MODULE_FUNCTIONAL: "GENE_ABUND",
    MODULE_ECOLOGICAL: "ECO",
    MODULE_PHENOTYPE: "ECO",
}

#: The engineering starting point for the legacy combined split.
DEFAULT_MODULE_WEIGHTS: Final = {
    MODULE_TAXONOMIC: 0.50,
    MODULE_FUNCTIONAL: 0.35,
    MODULE_ECOLOGICAL: 0.15,
}

COMBINED_SPLIT_NOTE: Final = (
    "The 50/35/15 module split is an engineering starting point, not an "
    "evidence-derived weighting."
)

FUSION_WEIGHTED_MODULES: Final = "weighted_modules"
FUSION_STUDY_AVERAGE: Final = "independent_study_module_average"
VALID_FUSION: Final = frozenset({FUSION_WEIGHTED_MODULES, FUSION_STUDY_AVERAGE})

STUDY_AVERAGE_NOTE: Final = (
    "Independent study modules are averaged with equal weight. Modules that share "
    "participants share one vote, so a study contributing twenty correlated taxa cannot "
    "outweigh a study contributing one replicated feature."
)

VALID_BODY_SITES: Final = frozenset({"stool", "scalp", "oral", "skin", "tumor", "vaginal"})


# --------------------------------------------------------------------------- #
# entities
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Feature:
    """One evidence assertion: what to measure, which way, and how it binds."""

    name: str
    module: str
    direction: str
    engine: str
    level: str | None
    #: evidence weight — effect precision, small cohort or single study lowers it
    w: float
    #: cohort independence — replication only within one recruitment family lowers it
    q: float
    #: disease specificity — also marking IBD, IBS, antibiotics lowers it
    s: float
    #: direction agreement — cohorts disagreeing on direction lowers it
    c: float
    sources: tuple[str, ...] = ()
    tier: str | None = None
    conflict_note: str | None = None
    note: str | None = None
    evidence_type: str = "TAX_REL"
    claim_level: str = "observed"
    #: Correlated features share a cluster and count once toward coverage.
    cluster: str | None = None
    #: "none" | "low" | "material" — declared by the profile author.
    direction_conflict: str = "none"
    #: For transported features: how confidently the source name maps here.
    mapping_confidence: str | None = None
    #: e.g. "not_below_0.05" for a feature the source reported at nominal p.
    nominal_status: str | None = None
    #: Which stratum this feature applies to (sex, subtype), or None for all.
    stratum: str | None = None
    #: The production catalogue name this source name is measured under, when
    #: the two differ (a reclassified species, an aggregated genus). ``name``
    #: stays the immutable identifier the source paper used; ``measured_as``
    #: is what this pipeline's profiler actually reports. Both are printed.
    measured_as: str | None = None
    #: Why the equivalence holds. Required whenever ``measured_as`` is set, so
    #: no mapping is ever silent.
    mapping_reason: str | None = None

    @property
    def d(self) -> int:
        """+1 for a reported increase, −1 for a decrease."""
        return 1 if self.direction == "increased" else -1

    @property
    def factor(self) -> float:
        """Scoring weight. The evidence weight alone — see the module docstring."""
        return self.w

    @property
    def audit_factor(self) -> float:
        """The legacy multiplied w·q·s·c, preserved for audit and never scored."""
        return self.w * self.q * self.s * self.c

    @property
    def key(self) -> str:
        """The name to look up in measurements and reference bundles."""
        return self.measured_as or self.name

    @property
    def display_name(self) -> str:
        """Source name, with the production name beside it when they differ."""
        return self.name if self.measured_as is None else f"{self.name} (measured as {self.measured_as})"

    @property
    def cluster_key(self) -> str:
        return self.cluster or self.name

    @property
    def bound(self) -> bool:
        return self.engine != "unbound"

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "module": self.module,
            "direction": self.direction,
            "engine": self.engine,
            "evidence_type": self.evidence_type,
            "claim_level": self.claim_level,
            "level": self.level,
            "w": self.w,
            "q": self.q,
            "s": self.s,
            "c": self.c,
            "scoring_weight": round(self.factor, 4),
            "audit_factor_wqsc": round(self.audit_factor, 4),
            "cluster": self.cluster_key,
            "direction_conflict": self.direction_conflict,
            "mapping_confidence": self.mapping_confidence,
            "nominal_status": self.nominal_status,
            "stratum": self.stratum,
            "measured_as": self.measured_as,
            "measurement_key": self.key,
            "mapping_reason": self.mapping_reason,
            "tier": self.tier,
            "sources": list(self.sources),
            "conflict_note": self.conflict_note,
            "note": self.note,
        }


@dataclass(frozen=True, slots=True)
class Module:
    """A set of assertions from one study or one coherent model, scored together."""

    name: str
    weight_in_combined: float
    features: tuple[Feature, ...]
    note: str = ""
    type: str = MODULE_DIRECTION_PANEL
    evidence_type: str = "TAX_REL"
    claim_level: str = "observed"
    assay_transport: str = "native"
    #: Modules sharing participants share a study group and one fused vote.
    study_group: str | None = None
    #: "validated" | "exploratory"
    status: str = "validated"
    maturity: str | None = None
    label: str = ""
    sources: tuple[str, ...] = ()
    #: Features naming the organism the module's finding is *about*. When a
    #: study's result is "X expanded in cases", the accessory taxa that moved
    #: with X are not the finding — X is. If every anchor is measured and
    #: absent, the pattern cannot be present whatever the accessory features
    #: do, and the module reports `anchor_absent` instead of a score.
    #:
    #: Without this, a module could report high resemblance to a
    #: Prevotella-expansion signature in a sample with no Prevotella, on the
    #: strength of one depleted Bacteroides species alone.
    anchor_features: tuple[str, ...] = ()

    @property
    def engine(self) -> str:
        engine, _ = EVIDENCE_TYPES.get(self.evidence_type, (None, ""))
        return engine or "unbound"

    @property
    def bound(self) -> bool:
        return self.engine != "unbound"

    @property
    def would_bind_if(self) -> str | None:
        if self.bound:
            return None
        return WOULD_BIND_IF.get(self.evidence_type, f"an engine for {self.evidence_type} is added")

    @property
    def fuses(self) -> bool:
        """Whether this module may contribute to the fused profile score."""
        return self.type == MODULE_DIRECTION_PANEL and self.status != "exploratory" and self.bound

    @property
    def group(self) -> str:
        return self.study_group or self.name

    @property
    def total_factor(self) -> float:
        return sum(f.factor for f in self.features)

    def to_json(self) -> dict[str, Any]:
        return {
            "module": self.name,
            "label": self.label or self.name,
            "type": self.type,
            "evidence_type": self.evidence_type,
            "claim_level": self.claim_level,
            "assay_transport": self.assay_transport,
            "transport_note": TRANSPORT_LABEL.get(self.assay_transport, ""),
            "study_group": self.group,
            "status": self.status,
            "maturity": self.maturity,
            "engine": self.engine,
            "bound": self.bound,
            "would_bind_if": self.would_bind_if,
            "fuses": self.fuses,
            "weight_in_combined": self.weight_in_combined,
            "note": " ".join(self.note.split()),
            "sources": list(self.sources),
            "n_features": len(self.features),
            "features": [f.to_json() for f in self.features],
        }


@dataclass(frozen=True, slots=True)
class CrossEngineCheck:
    """Two independent measurements of the same biology, compared."""

    taxonomic_feature: str
    functional_feature: str
    expect: str
    on_disagreement: str
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "taxonomic_feature": self.taxonomic_feature,
            "functional_feature": self.functional_feature,
            "expect": self.expect,
            "on_disagreement": self.on_disagreement,
            "note": " ".join(self.note.split()),
        }


@dataclass(frozen=True, slots=True)
class AbstentionRule:
    kind: str
    value: Any
    reason: str

    def to_json(self) -> dict[str, Any]:
        return {"kind": self.kind, "value": self.value, "reason": self.reason}


@dataclass(frozen=True, slots=True)
class ReferenceSpec:
    """Which contrast a module should be scored against."""

    primary_source: str | None
    primary_note: str
    secondary_source: str | None
    match_on: tuple[str, ...]

    @property
    def has_primary(self) -> bool:
        """A same-exposure, non-syndrome contrast, which most cohorts lack."""
        return bool(self.primary_source)

    def to_json(self) -> dict[str, Any]:
        return {
            "primary_source": self.primary_source,
            "primary_note": " ".join(self.primary_note.split()),
            "secondary_source": self.secondary_source,
            "match_on": list(self.match_on),
            "has_exposure_matched_primary": self.has_primary,
        }


@dataclass(frozen=True, slots=True)
class Stratum:
    """A duration, sex or subtype stratum with its own expectation."""

    name: str
    label: str
    criterion: str
    note: str = ""

    def to_json(self) -> dict[str, str]:
        return {
            "name": self.name,
            "label": self.label,
            "criterion": self.criterion,
            "note": " ".join(self.note.split()),
        }


@dataclass(frozen=True, slots=True)
class Confounder:
    name: str
    effect: str
    handling: str

    def to_json(self) -> dict[str, str]:
        return {"name": self.name, "effect": self.effect, "handling": self.handling}


@dataclass(frozen=True, slots=True)
class UnboundSource:
    """A published source the profile knows about but cannot yet measure."""

    source: str
    evidence_type: str
    note: str = ""

    @property
    def would_bind_if(self) -> str:
        return WOULD_BIND_IF.get(self.evidence_type, f"an engine for {self.evidence_type} is added")

    def to_json(self) -> dict[str, str]:
        return {
            "source": self.source,
            "evidence_type": self.evidence_type,
            "would_bind_if": self.would_bind_if,
            "note": " ".join(self.note.split()),
        }


@dataclass(frozen=True, slots=True)
class Profile:
    """A published disease-associated pattern, in declarative form."""

    name: str
    version: str
    label: str
    status: str
    summary: str
    modules: dict[str, Module]
    references: ReferenceSpec
    citation: str
    duration_dependent: bool = False
    strata: tuple[Stratum, ...] = ()
    cross_engine_checks: tuple[CrossEngineCheck, ...] = ()
    abstain_if: tuple[AbstentionRule, ...] = ()
    caveats: tuple[str, ...] = ()
    effect_size_note: str = ""
    combined_split_note: str = COMBINED_SPLIT_NOTE
    spec_version: int = 2
    disease: str = ""
    task: str = ""
    body_site: str = "stool"
    population: str = ""
    family: str = ""
    evidence_maturity: str = "P2"
    opt_in: bool = False
    #: Subject fields that must be supplied before this profile is scored.
    mandatory_match: tuple[str, ...] = ()
    #: Which subject field selects the stratum (e.g. "sex", "illness_duration_years").
    stratify_on: str | None = None
    confounders: tuple[Confounder, ...] = ()
    unbound_sources: tuple[UnboundSource, ...] = ()
    fusion: str = FUSION_WEIGHTED_MODULES
    #: Profiles this one is a specificity challenge for (spec 10.6).
    challenge_for: tuple[str, ...] = ()
    source_path: Path | None = field(default=None, compare=False)

    def features(self, module: str | None = None) -> tuple[Feature, ...]:
        if module is not None:
            return self.modules[module].features if module in self.modules else ()
        return tuple(f for m in self.modules.values() for f in m.features)

    def features_for_engine(self, engine: str) -> tuple[Feature, ...]:
        return tuple(f for f in self.features() if f.engine == engine)

    def modules_for_engine(self, engine: str) -> tuple[Module, ...]:
        return tuple(m for m in self.modules.values() if m.engine == engine)

    @property
    def fused_modules(self) -> tuple[Module, ...]:
        return tuple(m for m in self.modules.values() if m.fuses)

    @property
    def unbound_modules(self) -> tuple[Module, ...]:
        return tuple(m for m in self.modules.values() if not m.bound)

    @property
    def study_groups(self) -> dict[str, tuple[Module, ...]]:
        groups: dict[str, list[Module]] = {}
        for module in self.fused_modules:
            groups.setdefault(module.group, []).append(module)
        return {k: tuple(v) for k, v in groups.items()}

    @property
    def combined_weights(self) -> dict[str, float]:
        """Module -> weight in the fused score.

        Legacy profiles use their declared split. Evidence-typed profiles give
        every independent study group one equal vote, split evenly among the
        modules in that group.
        """
        if self.fusion == FUSION_WEIGHTED_MODULES:
            return {
                name: module.weight_in_combined
                for name, module in self.modules.items()
                if module.fuses and module.weight_in_combined > 0
            }
        groups = self.study_groups
        if not groups:
            return {}
        share = 1.0 / len(groups)
        weights: dict[str, float] = {}
        for members in groups.values():
            for module in members:
                weights[module.name] = share / len(members)
        return weights

    @property
    def fusion_note(self) -> str:
        return (
            self.combined_split_note
            if self.fusion == FUSION_WEIGHTED_MODULES
            else STUDY_AVERAGE_NOTE
        )

    @property
    def maturity_meaning(self) -> str:
        return MATURITY_MEANING.get(self.evidence_maturity, "")

    @property
    def is_transported(self) -> bool:
        return any(m.assay_transport != "native" for m in self.modules.values())

    @property
    def has_direction_conflict(self) -> bool:
        return any(f.direction_conflict == "material" for f in self.features())

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "spec_version": self.spec_version,
            "label": self.label,
            "status": self.status,
            "disease": self.disease,
            "task": self.task,
            "body_site": self.body_site,
            "population": self.population,
            "family": self.family,
            "evidence_maturity": self.evidence_maturity,
            "maturity_meaning": self.maturity_meaning,
            "opt_in": self.opt_in,
            "summary": " ".join(self.summary.split()),
            "duration_dependent": self.duration_dependent,
            "stratify_on": self.stratify_on,
            "mandatory_match": list(self.mandatory_match),
            "strata": [s.to_json() for s in self.strata],
            "references": self.references.to_json(),
            "fusion": {
                "method": self.fusion,
                "weights": self.combined_weights,
                "study_groups": {k: [m.name for m in v] for k, v in self.study_groups.items()},
                "note": " ".join(self.fusion_note.split()),
            },
            "modules": {name: module.to_json() for name, module in self.modules.items()},
            "unbound_sources": [u.to_json() for u in self.unbound_sources],
            "confounders": [c.to_json() for c in self.confounders],
            "challenge_for": list(self.challenge_for),
            "cross_engine_checks": [c.to_json() for c in self.cross_engine_checks],
            "abstain_if": [a.to_json() for a in self.abstain_if],
            "caveats": [" ".join(c.split()) for c in self.caveats],
            "effect_size_note": " ".join(self.effect_size_note.split()),
            "combined_split_note": " ".join(self.fusion_note.split()),
            "citation": " ".join(self.citation.split()),
        }


@dataclass(frozen=True, slots=True)
class ProfileSet:
    profiles: tuple[Profile, ...]

    def by_name(self, name: str) -> Profile:
        for profile in self.profiles:
            if profile.name == name:
                return profile
        raise PanelError(
            f"unknown profile {name!r}; available: "
            f"{', '.join(p.name for p in self.profiles)}"
        )

    def select(self, names: Sequence[str] | None, *, include_opt_in: bool = False) -> tuple[Profile, ...]:
        if names:
            return tuple(self.by_name(n) for n in names)
        return tuple(p for p in self.profiles if include_opt_in or not p.opt_in)

    @property
    def families(self) -> dict[str, tuple[Profile, ...]]:
        out: dict[str, list[Profile]] = {}
        for p in self.profiles:
            out.setdefault(p.family or p.disease or p.name, []).append(p)
        return {k: tuple(v) for k, v in out.items()}


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #


def _require(cond: bool, where: str, message: str) -> None:
    if not cond:
        raise PanelError(f"{where}: {message}")


def _unit(value: Any, where: str, key: str, default: float | None = None) -> float:
    if value is None and default is not None:
        return default
    _require(
        isinstance(value, (int, float)) and 0.0 <= float(value) <= 1.0,
        where,
        f"{key!r} must be a number in [0, 1], got {value!r}",
    )
    return float(value)


def _str_tuple(data: Mapping[str, Any], key: str, where: str) -> tuple[str, ...]:
    raw = data.get(key) or []
    if isinstance(raw, str):
        raw = [raw]
    _require(
        isinstance(raw, list) and all(isinstance(x, str) for x in raw),
        where,
        f"{key!r} must be a string or list of strings",
    )
    return tuple(raw)


def _parse_feature(
    raw: Any,
    *,
    module: str,
    profile: str,
    index: int,
    evidence_type: str,
    claim_level: str,
    transport: str,
) -> Feature:
    where = f"profile {profile!r} module {module!r} feature[{index}]"
    _require(isinstance(raw, dict), where, "each feature must be a mapping")

    name = str(raw.get("name") or raw.get("feature") or "").strip()
    _require(bool(name), where, "'name' must not be empty")

    direction = str(raw.get("direction", "")).strip()
    direction = {"higher_in_case": "increased", "lower_in_case": "decreased"}.get(direction, direction)
    _require(
        direction in VALID_DIRECTIONS,
        where,
        f"'direction' must be one of {sorted(VALID_DIRECTIONS)} (or higher_in_case / "
        f"lower_in_case), got {direction!r}",
    )

    default_engine, _ = EVIDENCE_TYPES.get(evidence_type, (None, ""))
    engine = str(raw.get("engine") or default_engine or "unbound").strip()
    _require(
        engine in VALID_ENGINES,
        where,
        f"'engine' must be one of {sorted(VALID_ENGINES)}, got {engine!r}",
    )

    conflict = str(raw.get("direction_conflict") or "none").strip()
    _require(
        conflict in ("none", "low", "material"),
        where,
        f"'direction_conflict' must be none, low or material, got {conflict!r}",
    )

    level = raw.get("level")
    level = str(level).strip() if level else None
    if engine == "metaphlan":
        _require(
            level in VALID_LEVELS,
            where,
            f"taxonomic features must declare 'level' as one of {sorted(VALID_LEVELS)}, "
            f"got {level!r}",
        )
        if name.strip().lower() in CONFLICTED_TAXA:
            _require(
                level == "genus" and conflict == "material" and transport == "amplicon_taxon_to_shotgun",
                where,
                f"{name!r} is on the conflicted-taxa list and has no stable direction when "
                "aggregated. Declare its constituent species, or — for a transported 16S "
                "finding only — declare it at 'level: genus' with 'direction_conflict: "
                "material' so the report can say so.",
            )
        if level == "genus":
            _require(
                transport == "amplicon_taxon_to_shotgun",
                where,
                f"genus-level feature {name!r} is only accepted in a module that declares "
                "'assay_transport: amplicon_taxon_to_shotgun'. Genus level is never accepted "
                "for shotgun-derived evidence: aggregating a genus reverses direction in practice.",
            )

    measured_as = raw.get("measured_as")
    mapping_reason = raw.get("mapping_reason")
    if measured_as is not None:
        measured_as = str(measured_as).strip()
        _require(bool(measured_as), where, "'measured_as' must not be empty when present")
        _require(
            measured_as != name,
            where,
            "'measured_as' equals 'name'; drop it rather than declaring an identity mapping",
        )
        _require(
            bool(mapping_reason and str(mapping_reason).strip()),
            where,
            "a feature that declares 'measured_as' must also declare 'mapping_reason': the "
            "equivalence has to be stated, never silent",
        )
        mapping_reason = str(mapping_reason).strip()
    else:
        _require(
            mapping_reason is None,
            where,
            "'mapping_reason' without 'measured_as' has nothing to explain",
        )

    w = _unit(raw.get("w"), where, "w", default=1.0)
    q = _unit(raw.get("q"), where, "q", default=1.0)
    s = _unit(raw.get("s"), where, "s", default=1.0)
    c = _unit(raw.get("c"), where, "c", default=1.0)
    _require(w > 0, where, "'w' must be greater than 0; drop the feature instead")

    return Feature(
        name=name,
        module=module,
        direction=direction,
        engine=engine,
        level=level,
        w=w,
        q=q,
        s=s,
        c=c,
        measured_as=measured_as,
        mapping_reason=mapping_reason,
        sources=_str_tuple(raw, "sources", where),
        tier=str(raw["tier"]).strip() if raw.get("tier") else None,
        conflict_note=str(raw["conflict_note"]).strip() if raw.get("conflict_note") else None,
        note=str(raw["note"]).strip() if raw.get("note") else None,
        evidence_type=evidence_type,
        claim_level=claim_level,
        cluster=str(raw["cluster"]).strip() if raw.get("cluster") else None,
        direction_conflict=conflict,
        mapping_confidence=(
            str(raw["mapping_confidence"]).strip() if raw.get("mapping_confidence") else None
        ),
        nominal_status=str(raw["nominal_status"]).strip() if raw.get("nominal_status") else None,
        stratum=str(raw["stratum"]).strip() if raw.get("stratum") else None,
    )


def _parse_module(raw: Any, *, name: str, profile: str, spec_version: int) -> Module:
    where = f"profile {profile!r} module {name!r}"
    _require(isinstance(raw, dict), where, "must be a mapping")

    legacy = spec_version == 2
    if legacy:
        _require(
            name in VALID_MODULES,
            where,
            f"unknown legacy module {name!r}; valid: {sorted(VALID_MODULES)}. "
            "Set 'spec_version: 3' to declare named evidence-typed modules.",
        )

    evidence_type = str(raw.get("evidence_type") or LEGACY_EVIDENCE.get(name, "")).strip()
    _require(
        evidence_type in EVIDENCE_TYPES,
        where,
        f"'evidence_type' must be one of {sorted(EVIDENCE_TYPES)}, got {evidence_type!r}",
    )
    _, native_claim = EVIDENCE_TYPES[evidence_type]
    claim_level = str(raw.get("claim_level") or native_claim).strip()
    _require(
        claim_level in VALID_CLAIM_LEVELS,
        where,
        f"'claim_level' must be one of {sorted(VALID_CLAIM_LEVELS)}, got {claim_level!r}",
    )
    # A module may only weaken a claim relative to its evidence type, never
    # strengthen it: gene abundance can be declared association_only but a
    # carrier proxy can never be declared observed.
    _require(
        not (native_claim == "carrier_proxy" and claim_level in ("observed", "external_measured")),
        where,
        "carrier abundance is 'carrier_proxy' and may not be declared observed or measured",
    )

    transport = str(raw.get("assay_transport") or "native").strip()
    _require(
        transport in VALID_TRANSPORTS,
        where,
        f"'assay_transport' must be one of {sorted(VALID_TRANSPORTS)}, got {transport!r}",
    )

    module_type = str(raw.get("type") or (MODULE_PHENOTYPE_PANEL if name == MODULE_PHENOTYPE else MODULE_DIRECTION_PANEL)).strip()
    _require(
        module_type in VALID_MODULE_TYPES,
        where,
        f"'type' must be one of {sorted(VALID_MODULE_TYPES)}, got {module_type!r}",
    )

    status = str(raw.get("status") or "validated").strip()
    _require(
        status in ("validated", "exploratory"),
        where,
        f"'status' must be validated or exploratory, got {status!r}",
    )
    if module_type == MODULE_SENSITIVITY_PANEL:
        status = "exploratory"

    maturity = raw.get("maturity")
    if maturity is not None:
        maturity = str(maturity).strip()
        _require(
            maturity in VALID_MATURITY,
            where,
            f"'maturity' must be one of {list(VALID_MATURITY)}, got {maturity!r}",
        )

    weight = raw.get("weight_in_combined", DEFAULT_MODULE_WEIGHTS.get(name, 0.0))
    if name == MODULE_PHENOTYPE or module_type != MODULE_DIRECTION_PANEL:
        weight = 0.0  # excluded from the fused figure by design
    weight = _unit(weight, where, "weight_in_combined", default=0.0)

    raw_features = raw.get("features") or []
    _require(isinstance(raw_features, list), where, "'features' must be a list")
    engine, _ = EVIDENCE_TYPES[evidence_type]
    if engine is not None:
        _require(bool(raw_features), where, "'features' must not be empty")
    # A module may declare one stratum for all its features (a female-only
    # cohort, say); a feature's own `stratum` still wins if it sets one.
    module_stratum = str(raw["stratum"]).strip() if raw.get("stratum") else None
    features = tuple(
        _parse_feature(
            {**f, "stratum": f.get("stratum") or module_stratum} if isinstance(f, dict) else f,
            module=name,
            profile=profile,
            index=i,
            evidence_type=evidence_type,
            claim_level=claim_level,
            transport=transport,
        )
        for i, f in enumerate(raw_features)
    )

    seen: set[str] = set()
    for feature in features:
        key = (feature.name, feature.stratum)
        _require(key not in seen, where, f"duplicate feature {feature.name!r}")
        seen.add(key)

    # An anchor must be one of the module's own features, or the gate could
    # never be evaluated and would silently do nothing.
    anchors = _str_tuple(raw, "anchor_features", where)
    names = {f.name for f in features}
    for anchor in anchors:
        _require(
            anchor in names,
            where,
            f"anchor_features names {anchor!r}, which is not a feature of this module",
        )

    return Module(
        name=name,
        weight_in_combined=weight,
        features=features,
        note=str(raw.get("note") or ""),
        type=module_type,
        evidence_type=evidence_type,
        claim_level=claim_level,
        assay_transport=transport,
        study_group=str(raw["study_group"]).strip() if raw.get("study_group") else None,
        status=status,
        maturity=maturity,
        label=str(raw.get("label") or ""),
        sources=_str_tuple(raw, "sources", where),
        anchor_features=anchors,
    )


def _parse_abstention(raw: Any, *, profile: str) -> tuple[AbstentionRule, ...]:
    where = f"profile {profile!r} abstain_if"
    if raw is None:
        return ()
    _require(isinstance(raw, list), where, "must be a list of single-key mappings")
    out: list[AbstentionRule] = []
    for index, item in enumerate(raw):
        _require(
            isinstance(item, dict) and len(item) >= 1,
            f"{where}[{index}]",
            "each rule must be a mapping",
        )
        kind = next(iter(k for k in item if k != "reason"))
        out.append(
            AbstentionRule(
                kind=kind,
                value=item[kind],
                reason=str(item.get("reason") or _default_abstention_reason(kind, item[kind])),
            )
        )
    return tuple(out)


def _default_abstention_reason(kind: str, value: Any) -> str:
    return {
        "recent_antibiotics_days_lt": (
            f"antibiotics within {value} days reshape the community enough that a profile "
            "score would measure the antibiotics"
        ),
        "usable_nonhost_reads_lt": (
            f"fewer than {value:,} usable non-host reads leaves the taxonomic profile too "
            "sparse to score" if isinstance(value, int) else "insufficient usable depth"
        ),
        "usable_nonhost_pairs_lt": (
            f"fewer than {value:,} usable non-host read pairs leaves the taxonomic profile "
            "too sparse to score" if isinstance(value, int) else "insufficient usable depth"
        ),
        "missing_metadata": (
            f"required metadata missing: {value}; without it the result cannot be placed "
            "in the right stratum"
        ),
        "outside_reference_envelope": (
            "the sample sits outside the multivariate reference distribution, so a "
            "percentile against that reference would be an extrapolation"
        ),
        "module_disagreement": (
            "the taxonomic and functional modules point opposite ways with no phenotype "
            "explanation, so at least one of them is wrong"
        ),
    }.get(kind, f"{kind} = {value}")


def _parse_confounders(raw: Any, *, where: str) -> tuple[Confounder, ...]:
    if not raw:
        return ()
    _require(isinstance(raw, list), where, "'confounders' must be a list")
    out: list[Confounder] = []
    for item in raw:
        _require(isinstance(item, dict) and item.get("name"), where, "each confounder needs a name")
        out.append(
            Confounder(
                name=str(item["name"]).strip(),
                effect=str(item.get("effect") or "").strip(),
                handling=str(item.get("handling") or "").strip(),
            )
        )
    return tuple(out)


def _parse_unbound(raw: Any, *, where: str) -> tuple[UnboundSource, ...]:
    if not raw:
        return ()
    _require(isinstance(raw, list), where, "'unbound_sources' must be a list")
    out: list[UnboundSource] = []
    for item in raw:
        _require(isinstance(item, dict) and item.get("source"), where, "each unbound source needs a source")
        evidence = str(item.get("evidence_type") or "").strip()
        _require(
            evidence in EVIDENCE_TYPES,
            where,
            f"unbound source {item['source']!r}: 'evidence_type' must be one of "
            f"{sorted(EVIDENCE_TYPES)}, got {evidence!r}",
        )
        out.append(
            UnboundSource(
                source=str(item["source"]).strip(),
                evidence_type=evidence,
                note=str(item.get("note") or ""),
            )
        )
    return tuple(out)


def parse_profile(data: Any, *, source_path: Path | None = None) -> Profile:
    origin = str(source_path) if source_path else "<inline>"
    _require(isinstance(data, dict), origin, "profile file must contain a YAML mapping")

    name = str(data.get("name", "")).strip()
    _require(bool(name), origin, "'name' must not be empty")
    where = f"profile {name!r}"
    _require(
        name.replace("_", "").replace("-", "").isalnum(),
        where,
        f"'name' must be alphanumeric plus '_' and '-', got {name!r}",
    )

    version = str(data.get("version", "")).strip()
    _require(
        bool(version) and version.count(".") >= 1,
        where,
        "'version' must be a semantic version; any weight change must bump it",
    )

    spec_version = data.get("spec_version", 2)
    _require(spec_version in (2, 3), where, f"'spec_version' must be 2 or 3, got {spec_version!r}")

    status = str(data.get("status", "research_only")).strip()
    _require(
        status in VALID_STATUS,
        where,
        f"'status' must be one of {sorted(VALID_STATUS)}, got {status!r}",
    )

    maturity = str(data.get("evidence_maturity") or "P2").strip()
    _require(
        maturity in VALID_MATURITY,
        where,
        f"'evidence_maturity' must be one of {list(VALID_MATURITY)}, got {maturity!r}",
    )

    body_site = str(data.get("body_site") or "stool").strip()
    _require(
        body_site in VALID_BODY_SITES,
        where,
        f"'body_site' must be one of {sorted(VALID_BODY_SITES)}, got {body_site!r}",
    )

    for key in ("label", "summary", "citation"):
        _require(
            isinstance(data.get(key), str) and bool(str(data[key]).strip()),
            where,
            f"{key!r} must be a non-empty string",
        )

    raw_modules = data.get("modules") or {}
    if isinstance(raw_modules, list):
        # Spec-3 files may write modules as a list with module_id keys.
        converted: dict[str, Any] = {}
        for item in raw_modules:
            _require(isinstance(item, dict) and item.get("module_id"), where, "list-form modules need 'module_id'")
            converted[str(item["module_id"])] = {k: v for k, v in item.items() if k != "module_id"}
        raw_modules = converted
    _require(isinstance(raw_modules, dict), where, "'modules' must be a mapping")
    _require(bool(raw_modules), where, "'modules' must not be empty")
    modules = {
        key: _parse_module(value, name=key, profile=name, spec_version=spec_version)
        for key, value in raw_modules.items()
    }

    # Same-participant modules declared independent is a compiler failure.
    for group, members in _group_modules(modules).items():
        declared = {m.study_group for m in members}
        _require(
            len(members) == 1 or None not in declared,
            where,
            f"modules {[m.name for m in members]} share study group {group!r} but at least one "
            "does not declare it; same-participant modules must all carry the study_group",
        )

    fusion_raw = data.get("fusion") or {}
    if isinstance(fusion_raw, str):
        fusion_raw = {"method": fusion_raw}
    _require(isinstance(fusion_raw, dict), where, "'fusion' must be a mapping")
    fusion = str(
        fusion_raw.get("method")
        or (FUSION_STUDY_AVERAGE if spec_version == 3 else FUSION_WEIGHTED_MODULES)
    ).strip()
    _require(fusion in VALID_FUSION, where, f"fusion method must be one of {sorted(VALID_FUSION)}")

    if fusion == FUSION_WEIGHTED_MODULES:
        combined = sum(m.weight_in_combined for m in modules.values() if m.fuses)
        _require(
            abs(combined - 1.0) < 1e-6 or combined == 0.0,
            where,
            f"combined module weights must sum to 1.0 (or all be 0), got {combined:.4f}",
        )

    raw_refs = data.get("references") or data.get("reference") or {}
    _require(isinstance(raw_refs, dict), where, "'references' must be a mapping")
    primary = raw_refs.get("primary") or {}
    secondary = raw_refs.get("secondary") or {}
    match_on = _str_tuple(secondary, "match_on", where) or _str_tuple(raw_refs, "match_on", where)
    references = ReferenceSpec(
        primary_source=str(primary.get("source")).strip() if primary.get("source") else None,
        primary_note=str(primary.get("note") or ""),
        secondary_source=(
            str(secondary.get("source")).strip()
            if secondary.get("source")
            else (str(raw_refs.get("cohort")).strip() if raw_refs.get("cohort") else None)
        ),
        match_on=match_on,
    )

    raw_checks = data.get("cross_engine_checks") or []
    _require(isinstance(raw_checks, list), where, "'cross_engine_checks' must be a list")
    checks: list[CrossEngineCheck] = []
    feature_names = {f.name for f in (f for m in modules.values() for f in m.features)}
    for index, item in enumerate(raw_checks):
        spot = f"{where} cross_engine_checks[{index}]"
        _require(isinstance(item, dict), spot, "must be a mapping")
        tax = str(item.get("taxonomic_feature", "")).strip()
        fun = str(item.get("functional_feature", "")).strip()
        for label, value in (("taxonomic_feature", tax), ("functional_feature", fun)):
            _require(bool(value), spot, f"{label!r} must not be empty")
            _require(
                value in feature_names,
                spot,
                f"{label} {value!r} is not a feature of this profile",
            )
        checks.append(
            CrossEngineCheck(
                taxonomic_feature=tax,
                functional_feature=fun,
                expect=str(item.get("expect", "concordant")).strip(),
                on_disagreement=str(item.get("on_disagreement", "downgrade_confidence")).strip(),
                note=str(item.get("note") or ""),
            )
        )

    raw_strata = data.get("strata") or []
    _require(isinstance(raw_strata, list), where, "'strata' must be a list")
    strata = tuple(
        Stratum(
            name=str(s.get("name", "")).strip(),
            label=str(s.get("label", "")).strip(),
            criterion=str(s.get("criterion", "")).strip(),
            note=str(s.get("note") or ""),
        )
        for s in raw_strata
        if isinstance(s, dict)
    )
    duration_dependent = bool(data.get("duration_dependent", False))
    if duration_dependent:
        _require(
            len(strata) >= 2,
            where,
            "'duration_dependent: true' requires at least two strata, since the point is "
            "that the strata differ and must not be averaged",
        )
    stratify_on = data.get("stratify_on")
    stratify_on = str(stratify_on).strip() if stratify_on else ("illness_duration_years" if duration_dependent else None)
    stratum_names = {s.name for s in strata}
    for feature in (f for m in modules.values() for f in m.features):
        _require(
            feature.stratum is None or feature.stratum in stratum_names,
            where,
            f"feature {feature.name!r} names stratum {feature.stratum!r} which is not declared",
        )

    return Profile(
        name=name,
        version=version,
        label=str(data["label"]).strip(),
        status=status,
        summary=str(data["summary"]).strip(),
        modules=modules,
        references=references,
        citation=str(data["citation"]).strip(),
        duration_dependent=duration_dependent,
        strata=strata,
        cross_engine_checks=tuple(checks),
        abstain_if=_parse_abstention(data.get("abstain_if"), profile=name),
        caveats=_str_tuple(data, "caveats", where),
        effect_size_note=str(data.get("effect_size_note") or ""),
        combined_split_note=str(data.get("combined_split_note") or COMBINED_SPLIT_NOTE),
        spec_version=int(spec_version),
        disease=str(data.get("disease") or "").strip(),
        task=str(data.get("task") or "").strip(),
        body_site=body_site,
        population=str(data.get("population") or "").strip(),
        family=str(data.get("family") or "").strip(),
        evidence_maturity=maturity,
        opt_in=bool(data.get("opt_in", False)),
        mandatory_match=_str_tuple(data, "mandatory_match", where),
        stratify_on=stratify_on,
        confounders=_parse_confounders(data.get("confounders"), where=where),
        unbound_sources=_parse_unbound(data.get("unbound_sources"), where=where),
        fusion=fusion,
        challenge_for=_str_tuple(data, "challenge_for", where),
        source_path=source_path,
    )


def _group_modules(modules: Mapping[str, Module]) -> dict[str, list[Module]]:
    groups: dict[str, list[Module]] = {}
    for module in modules.values():
        groups.setdefault(module.group, []).append(module)
    return groups


def load_profile_set(profiles_dir: Path) -> ProfileSet:
    """Every profile in `profiles_dir`, from a fingerprinted cache when unchanged."""
    from openbiota.fastcache import cached_load

    if not profiles_dir.is_dir():
        raise PanelError(f"profiles directory not found: {profiles_dir}")
    files = sorted(p for p in profiles_dir.glob("*.yaml") if not p.name.startswith("_"))
    return cached_load("profiles", files, lambda: _load_profile_set_uncached(profiles_dir))


def _load_profile_set_uncached(profiles_dir: Path) -> ProfileSet:
    files = sorted(p for p in profiles_dir.glob("*.yaml") if not p.name.startswith("_"))
    if not files:
        raise PanelError(f"no profile definitions found in {profiles_dir} (expected *.yaml)")

    profiles: list[Profile] = []
    seen: dict[str, Path] = {}
    for path in files:
        try:
            raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            raise PanelError(f"{path}: invalid YAML — {exc}") from exc
        except OSError as exc:
            raise PanelError(f"cannot read profile {path}: {exc}") from exc
        profile = parse_profile(raw, source_path=path)
        if profile.name in seen:
            raise PanelError(
                f"duplicate profile name {profile.name!r} in {path} "
                f"(already defined in {seen[profile.name]})"
            )
        seen[profile.name] = path
        profiles.append(profile)
    return ProfileSet(profiles=tuple(profiles))
