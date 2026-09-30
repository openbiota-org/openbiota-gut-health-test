"""The extension engine: assemble every new measurement for one sample.

BUILD_SPEC_v0.8.3 sections 3.2 and 13.1. This is what the CLI calls and what
writes `results.ext083.json`, `extension_manifest.json` and
`capability_coverage.json`.

Three properties it has to keep.

**It never touches a protected object.** Everything is written under one new
key. The preservation suite compares the protected keys before and after, and
this module gives it nothing to find.

**It degrades honestly.** A missing input makes one capability unavailable
with a stated reason; it does not fail the run or silently drop the feature.
Cached mode produces every view that existing structured data supports and
lists exactly what needs more analysis.

**It says what it is.** Every metric carries its dependency fingerprint, so a
cached value can be traced to the inputs, database releases and parameters
that produced it, and a changed reference invalidates only what depends on it.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.extension import SPEC_VERSION
from openbiota.extension import biotransform as BIO
from openbiota.extension import capacities as CAP
from openbiota.extension import communitytype as CT
from openbiota.extension import contexts as CX
from openbiota.extension import ecology as ECO
from openbiota.extension import explorer as EXP
from openbiota.extension import fermentation as FERM
from openbiota.extension import gusclasses as GUS
from openbiota.extension import inputs as INP
from openbiota.extension import longitudinal as LON
from openbiota.extension import nitrogen as NIT
from openbiota.extension import p9 as P9
from openbiota.extension import planner as PLAN
from openbiota.extension import readingscore as RSC
from openbiota.extension import registry as REG
from openbiota.extension import resistance as RES
from openbiota.extension import sources as SRC
from openbiota.extension import substrates as SUB
from openbiota.extension import vitamins as VIT
from openbiota.extension.schema import ExtensionMetric, fingerprint

#: The three execution modes of section 3.2.
MODES: Final[frozenset[str]] = frozenset({"cached", "incremental", "full-extension"})

RESULT_KEY: Final = "extension"


@dataclass(slots=True)
class Unavailable:
    """One capability that produced no value, and precisely why.

    Section 15.7 requires these reasons to stay apart: an unbuilt module and
    a sample that cannot answer are different facts, and only one of them is
    a limit of the science.
    """

    capability_id: str
    label: str
    reason: str
    detail: str
    what_would_unlock_it: str | None = None

    def to_json(self) -> dict[str, Any]:
        if self.reason not in REG.UNAVAILABLE_REASONS:
            raise ValueError(f"{self.capability_id}: unknown reason {self.reason!r}")
        return {
            "capability_id": self.capability_id,
            "label": self.label,
            "reason": self.reason,
            "is_engineering_debt": self.reason in REG.ENGINEERING_DEBT,
            "detail": self.detail,
            "what_would_unlock_it": self.what_would_unlock_it,
        }


@dataclass(slots=True)
class ExtensionResult:
    """Everything this release adds for one sample."""

    sample: str
    mode: str
    metrics: list[ExtensionMetric] = field(default_factory=list)
    views: dict[str, Any] = field(default_factory=dict)
    unavailable: list[Unavailable] = field(default_factory=list)
    manifest: dict[str, Any] = field(default_factory=dict)
    elapsed_s: float = 0.0

    def metric_ids(self) -> list[str]:
        return [m.metric_id for m in self.metrics]

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": "openbiota.extension-result/1.0",
            "spec_version": SPEC_VERSION,
            "sample": self.sample,
            "mode": self.mode,
            "n_metrics": len(self.metrics),
            "n_unavailable": len(self.unavailable),
            "metrics": [m.to_json() for m in self.metrics],
            "views": self.views,
            "unavailable": [u.to_json() for u in self.unavailable],
            "manifest": self.manifest,
            "elapsed_s": round(self.elapsed_s, 3),
        }


def _software_manifest() -> dict[str, Any]:
    """Program versions that can change a number."""
    versions: dict[str, str] = {
        "python": platform.python_version(),
        "openbiota_spec": SPEC_VERSION,
    }
    for module in ("numpy", "yaml", "micom", "cobra", "optlang"):
        try:
            imported = __import__(module)
            versions[module] = str(getattr(imported, "__version__", "unknown"))
        except ImportError:
            versions[module] = "not installed"
    return versions


def _input_manifest(results: Mapping[str, Any]) -> dict[str, Any]:
    """What the extension read, so a cached value can be traced to it."""
    run = results.get("run") or {}
    inventory = results.get("organism_inventory") or {}
    return {
        "sample": results.get("sample"),
        "baseline_version": results.get("version"),
        "reference_fingerprint": run.get("reference fingerprint"),
        "n_inventory_organisms": inventory.get("n_organisms"),
        "primary_lane": inventory.get("primary_lane"),
        "pathogen_bundle_id": (results.get("pathogens") or {}).get("bundle_id"),
        "mycobiome_schema": (results.get("mycobiome") or {}).get("schema_version"),
    }


# --------------------------------------------------------------------------- #
# the capability runners
# --------------------------------------------------------------------------- #


def _taxon_key(name: str) -> str:
    """The shared taxon normaliser, so this join cannot drift from the others."""
    from openbiota.extension.synbiotic import taxon_key  # noqa: PLC0415

    return taxon_key(name)


def _functional_redundancy(results: Mapping[str, Any]) -> list[dict[str, Any]]:
    """How many organisms hold each measured function — spec §6.4.

    The carriers are the organisms the panel already attributed its
    fragments to, so this asks a question of evidence that has already been
    gathered: is this function held by one dominant organism or spread
    across several? Two carriers at 50/50 give an effective count of two;
    two at 99/1 give barely more than one, which is the distinction that
    matters and the reason a raw carrier count will not do.

    Only the aggregated targets count. A decoy's organisms are the
    background the panel is discriminating against, not carriers of the
    function.
    """
    # A read's best match being some organism's copy of a gene does not put
    # that organism in this sample: the reference set may simply not hold the
    # one that is here. So carriage requires the organism to be in the
    # sample's own inventory as well, and anything else is recorded as
    # context. That join is why the counts below are lower than the number of
    # reference organisms a panel touches.
    inventory = {
        _taxon_key(str(o.get("species") or ""))
        for o in (results.get("organism_inventory") or {}).get("organisms") or []
    }
    out: list[dict[str, Any]] = []
    for panel in results.get("panels") or []:
        if not isinstance(panel, dict):
            continue
        aggregated = set(panel.get("aggregate_from") or ())
        carriers: dict[str, float] = {}
        for gene in panel.get("genes") or []:
            if aggregated and gene.get("entry_id") not in aggregated:
                continue
            for organism in gene.get("organisms") or []:
                name = str(organism.get("organism") or "").strip()
                if name:
                    carriers[name] = carriers.get(name, 0.0) + float(
                        organism.get("fragments") or 0.0)
        if not carriers:
            continue
        rows = [
            {
                "organism": name,
                "value": fragments,
                "is_carrier_evidence": _taxon_key(name) in inventory,
            }
            for name, fragments in carriers.items()
        ]
        attributed = sum(carriers.values())
        accepted = float(panel.get("accepted_fragments") or 0.0)
        reading = ECO.redundancy(
            f"panels.{panel['name']}",
            str(panel.get("metabolite") or panel["name"]),
            rows,
            unassigned_fraction=(
                max(0.0, 1.0 - attributed / accepted) if accepted else None
            ),
        )
        out.append(reading.to_json())
    return out


def _redundancy_metrics(readings: Sequence[Mapping[str, Any]]) -> list[ExtensionMetric]:
    """§6.4's effective carrier count, one metric per function.

    These are measurements of this sample and belong in the index with the
    rest, not only in the section that draws them: a reader looking up
    whether a function is held by one organism should find it where every
    other measurement is listed.

    The count is a number of organisms, not a score, so it carries no
    `score_0_100` - the bar the section draws is a presentation choice
    about where to stop counting, and that choice does not belong here.
    """
    out: list[ExtensionMetric] = []
    for reading in readings:
        effective = reading.get("effective_carriers")
        if effective is None:
            continue
        function_id = str(reading.get("function_id") or "")
        coverage = reading.get("carrier_coverage")
        out.append(ExtensionMetric(
            metric_id=f"A09.redundancy.{function_id}",
            label=f"{reading.get('label') or function_id}: organisms holding it",
            kind="community_composition",
            state="measured",
            method_id=ECO.METHOD_REDUNDANCY,
            feature_id="A09",
            view_ids=("ecology",),
            value=float(effective),
            unit="effective organisms",
            denominator="carriers resolved in this sample's own inventory",
            direction="higher_favourable_in_context",
            direction_context=(
                "a function held by several organisms survives losing one of them"
            ),
            # Under half the signal placed makes the count a floor, and a
            # floor is a weak reading of the number it stands in for.
            analytical_confidence="weak" if (coverage or 0) < 0.5 else "supported",
            evidence_maturity="mechanistic_hypothesis",
            limitations=(
                "Redundancy among the carriers this assay resolved, not a guarantee of "
                "ecological resilience.",
                *(
                    ("Most of this function's signal sits on organisms that could not be "
                     "placed in this sample, so the count is a floor.",)
                    if (coverage or 0) < 0.5 else ()
                ),
            ),
            extra={
                "n_resolved_carriers": reading.get("n_resolved_carriers"),
                "carrier_coverage": coverage,
            },
            input_fingerprint=fingerprint("A09", "redundancy", function_id, effective),
        ))
    return out


def _run_ecology(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A09 — companion diversity, dominance, ratios and aerotolerance."""
    inventory = results.get("organism_inventory") or {}
    organisms = inventory.get("organisms") or []
    if not organisms:
        out.unavailable.append(Unavailable(
            "A09", "Ecology companion dashboard", "missing_input",
            "the organism inventory is empty for this sample",
            "a completed taxonomic profile",
        ))
        return

    community = results.get("community_profile") or {}
    vector = ECO.vector_from_inventory(
        organisms,
        unresolved_percent=inventory.get("unclassified_percent"),
        depth_fragments=community.get("total_rpob_fragments"),
    )
    diversity = ECO.diversity(vector)

    by_level: dict[str, ECO.AbundanceVector] = {}
    phyla = community.get("phyla") or []
    if phyla:
        by_level["phylum"] = ECO.phylum_vector(
            {str(p.get("label")): float(p.get("percent") or 0.0) for p in phyla},
            lane="community_profile.rpoB",
            catalogue="rpoB single-copy marker attribution",
        )
    ratios = ECO.ratios(vector, by_level=by_level)

    aerotolerance = None
    try:
        aerotolerance = ECO.aerotolerance(vector, ECO.load_oxygen_traits())
    except OSError as exc:
        out.unavailable.append(Unavailable(
            "A09.aerotolerance", "Aerotolerance balance", "missing_input",
            f"the oxygen trait table could not be read: {exc}",
            "extension/oxygen_traits.yaml",
        ))

    # The composition view needs no model and is always available; the
    # classifier needs the frozen artifact and says so when it is absent.
    abundances = {
        str(o.get("species") or ""): float(o.get("percent") or 0.0)
        for o in organisms if o.get("percent")
    }
    model: CT.CommunityTypeModel | None = None
    try:
        model = CT.CommunityTypeModel.load()
    except CT.CommunityTypeError as exc:
        out.unavailable.append(Unavailable(
            "A09.community_type", "Gut community type", "missing_input", str(exc)[:240],
            "scripts/train_community_type.py",
        ))
    # The composition view names every genus it finds; only the classifier
    # works over the model's fixed universe. Restricting the chart to the
    # model would print a 60% "unresolved" bar that is an artefact of the
    # reference, not a property of the sample.
    comp = CT.composition(abundances, denominator=f"{vector.catalogue} ({vector.lane})")
    assignment = CT.assign(comp, model) if model else None

    redundancy = _functional_redundancy(results)
    out.metrics.extend(ECO.metrics(vector, diversity, ratios, aerotolerance))
    out.metrics.extend(_community_type_metrics(comp, assignment, model))
    out.metrics.extend(_redundancy_metrics(redundancy))
    out.views["ecology"] = {
        "feature_id": "A09",
        "vector": vector.declaration(),
        "diversity": diversity.to_json(),
        "ratios": [r.to_json() for r in ratios],
        "aerotolerance": aerotolerance.to_json() if aerotolerance else None,
        "composition": comp.to_json(),
        # §6.4: how many organisms hold each function, which is a different
        # question from how much of the function there is. The module was
        # written and exported and nothing ever called it, so the reading
        # existed in the code and in no report.
        "redundancy": redundancy,
        "community_type": assignment.to_json() if assignment else None,
        "community_type_model": {
            "model_id": model.model_id, "k": model.k, "stable": model.stable,
            "groups": list(model.labels), "evaluation": dict(model.evaluation),
            "cohort": dict(model.cohort),
        } if model else None,
        "legacy_values_unchanged": {
            "shannon_index": community.get("shannon_index"),
            "simpson_index": community.get("simpson_index"),
            "pielou_evenness": community.get("pielou_evenness"),
            "firmicutes_bacteroidetes_ratio": community.get("firmicutes_bacteroidetes_ratio"),
        },
    }
    if assignment is not None and assignment.state in {"mixed", "continuous_only", "unavailable"}:
        out.unavailable.append(Unavailable(
            "A09.community_type.label", "A named gut community type", "not_applicable",
            assignment.note, None,
        ))


def _community_type_metrics(
    comp: CT.Composition,
    assignment: CT.Assignment | None,
    model: CT.CommunityTypeModel | None,
) -> list[ExtensionMetric]:
    """The composition's headline shares and the distance to each group.

    A distance is a measurement; a group name is not a number and stays in
    the view. The unresolved share is reported as its own metric because a
    composition that names 55% of a sample is a different statement from one
    that names 95%.
    """
    metrics: list[ExtensionMetric] = []
    if comp.total_percent <= 0:
        return metrics
    metrics.append(ExtensionMetric(
        metric_id="ext083.ecology.composition_unresolved_share",
        label="Share of the community outside the composition's genus universe",
        kind="community_composition", state="measured",
        method_id=CT.METHOD_COMPOSITION,
        value=comp.unresolved_share, unit="fraction of 1",
        denominator=comp.denominator, feature_id="A09", view_ids=("ecology",),
        analytical_confidence="supported", evidence_maturity="not_applicable",
        limitations=(
            "A partition of the same total, not a residual that was dropped.",
            "Organisms outside the genus universe are not absent from the sample; they are "
            "outside the reference the composition is drawn against.",
        ),
        input_fingerprint=fingerprint("A09.composition", sorted(comp.shares.items())),
    ))
    for genus, share in comp.top(6):
        metrics.append(ExtensionMetric(
            metric_id=f"ext083.ecology.composition.{genus.lower()}",
            label=f"{genus} share of the community",
            kind="community_composition", state="measured",
            method_id=CT.METHOD_COMPOSITION, value=share, unit="fraction of 1",
            denominator=comp.denominator, feature_id="A09", view_ids=("ecology",),
            analytical_confidence="supported", evidence_maturity="not_applicable",
            limitations=(
                "A relative share of one lane's abundance, not an absolute count or load.",
                "A dominant genus is a composition, not a diet, a diagnosis or a cause.",
            ),
            input_fingerprint=fingerprint("A09.composition", genus, round(share, 12)),
        ))
    if assignment is None or model is None:
        return metrics
    for group, distance in assignment.distances:
        metrics.append(ExtensionMetric(
            metric_id=f"ext083.ecology.community_type.distance.{group}",
            label=f"Distance to the {group} reference group",
            kind="reference_percentile", state="measured",
            method_id=CT.METHOD_COMMUNITY_TYPE, value=distance,
            unit="sqrt(Jensen-Shannon divergence), nats^0.5",
            reference_id=model.model_id, reference_state="available",
            feature_id="A09", view_ids=("ecology",),
            analytical_confidence="supported", evidence_maturity="association_only",
            limitations=(
                "A resemblance to a reference composition, not a diagnosis.",
                "Relative composition cannot establish absolute microbial load.",
            ),
            input_fingerprint=fingerprint("A09.community_type", model.model_id, group,
                                          sorted(comp.shares.items())),
        ))
    return metrics


def _run_explorer(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A10 - the searchable explorer, rollups and taxonomy resolution."""
    organisms = (results.get("organism_inventory") or {}).get("organisms") or []
    if not organisms:
        out.unavailable.append(Unavailable(
            "A10", "Complete organism explorer", "missing_input",
            "the organism inventory is empty for this sample",
            "a completed taxonomic profile",
        ))
        return
    # Every named screening target from the capability manifest, so a
    # reader who looks one up gets an answer either way (spec 7.2).
    registry = REG.load()
    named = [
        {"id": c.capability_id, "concept": c.concept}
        for c in registry.of_kind("measurement_view")
        if "A10" in str(c.binding or "") and _looks_like_a_taxon(c.concept)
    ]
    view = EXP.build(results, named_targets=named)
    out.views["explorer"] = view
    out.metrics.extend(EXP.metrics(
        view, input_fingerprint=fingerprint(len(organisms), view["n_genera"]),
    ))


def _run_resistance(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A11 — class, mechanism and carrier tables over the existing calls."""
    pathogens = results.get("pathogens")
    if not pathogens or not (pathogens.get("determinants") or []):
        out.unavailable.append(Unavailable(
            "A11", "Resistance overview", "missing_input",
            "the pathogen screen did not run, so there are no determinant calls to group",
            "a completed pathogen screen",
        ))
        return
    view = RES.build(pathogens)
    classes = RES.summarise(RES.read_determinants(pathogens))
    view["carriage_note"] = RES.unlinked_note(classes)
    out.views["resistance"] = view
    out.metrics.extend(RES.metrics(
        classes,
        input_fingerprint=fingerprint(pathogens.get("bundle_id"), view["totals"]),
    ))


def _run_simulation(
    results: Mapping[str, Any],
    out: ExtensionResult,
    *,
    mode: str,
    simulate: str = "auto",
    lp_workers: int | None = None,
    progress: Any = None,
) -> None:
    """A16 — the model-assisted scenarios, solved, and their readiness.

    Every number the scenarios produce is the optimum of a linear programme,
    solved in parallel and cached per question. `simulate` is the feature
    flag: `never` reports readiness only; `auto` solves what the cache cannot
    serve (and in plain `cached` mode answers from the cache alone); `always`
    solves even when the cache would serve.
    """
    from openbiota.extension import simulation as SIM
    from openbiota.extension import synbiotic as SYN

    try:
        config = SIM.SimulationConfig.load()
    except (OSError, KeyError) as exc:
        out.unavailable.append(Unavailable(
            "A16.model", "Model-assisted scenarios", "missing_input",
            f"the simulation configuration could not be read: {exc}", None,
        ))
        return

    readiness = SIM.readiness(config)
    replay: dict[str, Any] | None = None
    if readiness["replay_fixture_installed"]:
        try:
            replay = SIM.replay_published_fixture(config).to_json()
        except SIM.SimulationError as exc:
            replay = {"error": str(exc)}

    view: dict[str, Any] = {
        "feature_id": "A16",
        "readiness": readiness,
        "protocol": config.protocol.to_json(),
        "substrates": [s.to_json() for s in config.substrates],
        "consortium": list(config.added_species),
        "media": [{k: v for k, v in m.items() if k != "sha256"} for m in config.media],
        "published_replay": replay,
        # §24's seventeen records are published trials, not model output, so
        # they belong in the view whether or not a solver is installed. Four
        # are negative and eight never ran a component-only arm; both facts
        # are findings, and dropping them would leave only the flattering ones.
        "evidence_seeds": SYN.seed_summary(),
        "feature_flag": simulate,
        "scenarios": None,
        "limitations": [
            "An experimental prediction from a metabolic model, not a measurement and not "
            "a guarantee of an individual response.",
            "Predicted flux is not a stool concentration, an absorbed exposure or a dose.",
        ],
    }
    out.views["simulation"] = view

    if not readiness["can_run"]:
        missing = [
            name for name, state in (
                ("MICOM", readiness["micom"]["installed"]),
                ("a linear-programming solver", readiness["solver"]["installed"]),
                ("the AGORA2 model library", readiness["model_library"]["installed"]),
            ) if not state
        ]
        out.unavailable.append(Unavailable(
            "A16.model", "Model-assisted scenarios", "missing_input",
            "the simulation layer is implemented but " + ", ".join(missing) + " is not installed",
            "make extension-refs",
        ))
        return

    if simulate == "never":
        out.unavailable.append(Unavailable(
            "A16.model", "Model-assisted scenarios", "missing_input",
            "the run was started with --simulate never; the solver, library and method are "
            "ready and nothing was solved",
            "--simulate auto",
        ))
        return

    # Cached mode may only answer from the cache; the other modes may solve.
    dry_run = (mode == "cached" and simulate != "always")
    try:
        scenarios = SYN.simulate(
            results, config, sample_id=str(out.sample), workers=lp_workers,
            progress=progress, dry_run=dry_run,
        )
    except SIM.SimulationError as exc:
        out.unavailable.append(Unavailable(
            "A16.model", "Model-assisted scenarios", "missing_input", str(exc)[:300], None,
        ))
        return
    view["scenarios"] = scenarios
    if scenarios["execution_state"] == "not_run":
        out.unavailable.append(Unavailable(
            "A16.model", "Model-assisted scenarios", "missing_input",
            scenarios["execution_state_meaning"], "--simulate auto with full-extension mode",
        ))
    elif scenarios["execution_state"] == "partial":
        out.unavailable.append(Unavailable(
            "A16.model.partial", "Some model-assisted scenarios", "missing_input",
            scenarios["execution_state_meaning"], None,
        ))
    for metric in _scenario_metrics(scenarios):
        out.metrics.append(metric)
    scores = _coverage_scores(scenarios)
    # The renderer reads the view, not the metric list, so the numbers have to
    # be here as well: a report built from the json must never need to
    # recompute one.
    view["coverage"] = scores
    for metric in _coverage_metrics(scores):
        out.metrics.append(metric)


def _coverage_scores(scenarios: Mapping[str, Any]) -> dict[str, Any]:
    """§11.4's C and P over the modelled candidates, as one record."""
    from openbiota.extension import synbiotic as SYN  # noqa: PLC0415

    baseline = scenarios.get("baseline_medium")
    goals = [
        {
            "goal": f"butyrate production on {c['substrate']['label']}",
            "weight": 1.0,
            "favourable": str(c.get("verdict")) in SYN.FAVOURABLE_VERDICTS,
            "pair_supported": str(c.get("verdict")) in SYN.PAIR_SUPPORTED_VERDICTS,
            "conflicted": str(c.get("verdict")) == "pair_worse",
            "evidence_setting": "in-silico community model, not a trial",
        }
        for c in scenarios.get("candidates") or []
        # A candidate whose model never solved is not an uncovered goal.
        # Counting it as one would read as "tried and found no benefit",
        # which is the opposite of what an unsolved arm means, and would
        # drag the percentage down by exactly the arithmetic §11.4 forbids.
        if c.get("medium_id") == baseline and str(c.get("verdict")) != "unavailable"
    ]
    try:
        return SYN.coverage(goals, n_components=2)
    except SYN.CoverageError as exc:
        return {
            "goal_coverage": None,
            "pair_supported_goal_coverage": None,
            "n_goals": len(goals),
            "unavailable_reason": str(exc),
            "goals": [],
        }


def _coverage_metrics(scores: Mapping[str, Any]) -> list[ExtensionMetric]:
    """§11.4's two coverage scores, over the modelled candidates.

    C asks what fraction of weighted goals the option has a favourable edge
    for; P asks the narrower question of whether the *pair* is what carries
    the endpoint. Both are percentages of a weighted goal set, so they are
    already 0-100 and become `score_0_100` directly.

    These describe the strength of an in-silico case, not a measurement of
    the sample, so they carry no reference percentile: there is no cohort of
    other people's coverage scores to stand next to.
    """
    common = {
        "kind": "model_prediction",
        "method_id": "A16.coverage/11.4",
        "feature_id": "A16",
        "view_ids": ("simulation",),
        "source_ids": ("S10", "S11"),
        "analytical_confidence": "provisional",
        "evidence_maturity": "mechanistic_hypothesis",
        "unit": "% of weighted goals",
        "score_definition_id": "openbiota.coverage/11.4",
    }
    limits = (
        "Coverage of a modelled goal set, not a clinical response rate.",
        "Built from in-silico predictions; no trial evidence is mixed in.",
    )
    out: list[ExtensionMetric] = []
    for key, label, note in (
        ("goal_coverage", "Goal coverage (C)",
         "share of weighted goals with a favourable modelled edge"),
        ("pair_supported_goal_coverage", "Pair-supported goal coverage (P)",
         "share where the combination, not either part alone, carries the endpoint"),
    ):
        value = scores.get(key)
        if value is None:
            reason = scores.get("pair_coverage_unavailable_reason") or scores.get(
                "unavailable_reason"
            ) or "not applicable"
            out.append(ExtensionMetric(
                metric_id=f"A16.{key}", label=label, state="not_applicable",
                value=None, limitations=(*limits, str(reason)),
                input_fingerprint=fingerprint("A16", key, "unavailable"),
                **common,
            ))
            continue
        out.append(ExtensionMetric(
            metric_id=f"A16.{key}", label=label, state="measured",
            value=float(value), score_0_100=float(value),
            direction="higher_favourable_in_context", direction_context=note,
            limitations=limits,
            input_fingerprint=fingerprint("A16", key, scores.get("n_goals")),
            extra={"n_goals": scores.get("n_goals"), "goals": scores.get("goals")},
            **common,
        ))
    return out


def _scenario_metrics(scenarios: Mapping[str, Any]) -> list[ExtensionMetric]:
    """One metric per arm on the baseline diet, plus each candidate's interaction.

    Model predictions carry no reference percentile and no direction: a
    higher butyrate potential is not "favourable" without a context the
    model does not have.
    """
    metrics: list[ExtensionMetric] = []
    baseline = scenarios.get("baseline_medium")
    protocol_id = str((scenarios.get("protocol") or {}).get("protocol_id", ""))
    common = {
        "kind": "model_prediction",
        "method_id": f"A16.lp-bounds/{protocol_id}",
        "feature_id": "A16",
        "view_ids": ("simulation",),
        "source_ids": ("S10", "S11"),
        "analytical_confidence": "provisional",
        "evidence_maturity": "mechanistic_hypothesis",
    }
    limits = (
        "An experimental prediction from a metabolic model, not a measurement.",
        "Not a stool concentration, an absorbed exposure or a dose.",
    )
    for candidate in scenarios.get("candidates") or []:
        if candidate.get("medium_id") != baseline:
            continue
        sid = candidate["substrate"]["substrate_id"]
        label = candidate["substrate"]["label"]
        for arm_name, arm in (candidate.get("arms") or {}).items():
            value = arm.get("raw_flux")
            metric_id = f"A16.{sid}.{arm_name}.butyrate_potential"
            arm_label = f"{label} — {arm_name.replace('_', ' ')}: butyrate potential"
            if value is None:
                metrics.append(ExtensionMetric(
                    metric_id=metric_id, label=arm_label, state="not_applicable",
                    value=None, unit="mmol/gDW/h",
                    limitations=(*limits, str(arm.get("unavailable_reason") or "not solved")),
                    input_fingerprint=fingerprint("A16", sid, arm_name, "unavailable"),
                    **common,
                ))
                continue
            metrics.append(ExtensionMetric(
                metric_id=metric_id, label=arm_label, state="measured",
                value=float(value), unit="mmol/gDW/h",
                limitations=limits,
                input_fingerprint=fingerprint(
                    "A16", sid, arm_name, arm.get("medium_id"), protocol_id,
                ),
                extra={
                    "lower_bound": arm.get("raw_flux_lower_bound"),
                    "community_growth": arm.get("community_growth"),
                    "flux_per_community_growth": arm.get("flux_per_community_growth"),
                    "reproducible": (arm.get("reproducibility") or {}).get("agrees"),
                },
                **common,
            ))
        interaction = (candidate.get("contrasts") or {}).get("interaction")
        if interaction is not None:
            metrics.append(ExtensionMetric(
                metric_id=f"A16.{sid}.interaction", state="measured",
                label=f"{label} + consortium: interaction contrast",
                value=float(interaction), unit="mmol/gDW/h",
                limitations=(*limits, "F11 - F10 - F01 + F00 in the model; not proof of "
                             "clinical synergy."),
                input_fingerprint=fingerprint("A16", sid, "interaction", protocol_id),
                extra={"verdict": candidate.get("verdict")},
                **common,
            ))
    return metrics


#: Which CAZy families each measured enzyme belongs to. The carbohydrate
#: panel searches by characterised enzyme name; the A01 specificity rules
#: are written in families, because that is how the substrate evidence is
#: published. One enzyme can sit in more than one family, and a family can
#: hold more than one enzyme; both are true here and neither is flattened.
CAZY_FAMILY: Final[Mapping[str, tuple[str, ...]]] = {
    # cellulose and mixed-linkage glucan
    "celA": ("GH5",), "celB": ("GH9",), "bglA": ("GH3",),
    "licB": ("GH16",), "glcA": ("GH55",),
    # starch and alpha-glucosides
    "amyA": ("GH13",), "pulA": ("GH13",), "gluA": ("GH15",),
    "aglA": ("GH31",), "malL": ("GH13_31", "GH13"),
    # chitin
    "chiA": ("GH18",), "chiD": ("CE4",), "nagZ": ("GH20",),
    # pectin
    "pgl": ("GH28",), "pelA": ("PL1",), "pme": ("CE8",), "rgl": ("GH28",),
    # fructans
    "sacA": ("GH32",), "levB": ("GH32",),
    # galactosides
    "lacZ": ("GH2",), "melA": ("GH27",),
    # xylan and arabinoxylan
    "xynA": ("GH10",), "xylB": ("GH43",), "abfA": ("GH51",),
    "axe": ("CE6",), "fae": ("CE1",),
    # mannans
    "manB": ("GH26",), "manA": ("GH2", "GH113"),
}

#: Where `make substrate-refs` puts the dbCAN snapshot.
REFS_ROOT: Final = Path(__file__).resolve().parents[2] / "refs"
DBCAN_DIR: Final = REFS_ROOT / "dbcan"


def _run_substrates(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A01 - the fourteen substrate views, and what each may claim.

    The discrimination rules are built and tested; the enzyme evidence they
    grade comes from a dbCAN snapshot that `make substrate-refs` installs.
    Until it is there, every card reads `not_assayed`, which is a different
    statement from "these enzymes were not found".
    """
    mapping_file = DBCAN_DIR / "fam-substrate-mapping.tsv"
    installed = mapping_file.is_file()
    reconciliation: dict[str, Any] | None = None
    if installed:
        try:
            reconciliation = SUB.validate_against_dbcan(
                SUB.load_dbcan_mapping(mapping_file),
            )
        except (OSError, SUB.SubstrateError) as exc:
            reconciliation = {"error": str(exc)[:300], "agrees": False}

    # The carbohydrate panel measures enzymes by gene symbol; the substrate
    # rules are written in CAZy families, because that is the vocabulary the
    # specificity evidence is published in. This is the join between them,
    # written out so a reviewer can check each assignment rather than trust
    # a lookup buried in a data file.
    hits_by_family: dict[str, list[SUB.GeneHit]] = {}
    panels = results.get("panels") or []
    for panel in panels if isinstance(panels, list) else []:
        if str(panel.get("name") or "").lower() not in {"carbohydrates", "cazy"}:
            continue
        for gene in panel.get("genes") or []:
            if str(gene.get("role")) != "target":
                continue
            fragments = int(gene.get("fragments") or 0)
            if fragments <= 0 or str(gene.get("confidence")) in {"not detected", "absent"}:
                continue
            for family in CAZY_FAMILY.get(str(gene.get("gene") or ""), ()):
                hits_by_family.setdefault(family, []).append(SUB.GeneHit(
                    gene_id=str(gene.get("gene")), family=family, fragments=fragments,
                ))
    assayed = bool(hits_by_family)
    findings = SUB.assess_all(hits_by_family, assayed=assayed)

    out.views["substrates"] = {
        "feature_id": "A01",
        "n_substrates": len(findings),
        "assayed": assayed,
        "reference_installed": installed,
        "reference_reconciliation": reconciliation,
        "substrates": [f.to_json() for f in findings],
        "shared_genes": SUB.shared_gene_report(findings),
        "specificity_states": list(SUB.SPECIFICITY),
        "limitations": [
            "Genetic capacity to act on a substrate, not a measurement of anything "
            "digested and not a food tolerance.",
            "A gene family shared between two substrates appears on both cards and is "
            "counted once in any total.",
        ],
    }
    if not assayed:
        out.unavailable.append(Unavailable(
            "A01", "Dietary substrate panel", "missing_input",
            "the fourteen substrate views and their specificity rules are implemented, and "
            "the enzyme reference set they grade is not installed"
            + ("" if installed else "; the dbCAN family-substrate mapping is absent"),
            "make substrate-refs",
        ))


def _run_fermentation(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A02 - the route breakdown and the cross-feeding network.

    The existing butyrate, propionate, methane, sulfide and TMA scores are
    read, not recomputed: this adds the routes behind them.
    """
    genes, assayed = _detected_genes(results)
    view = FERM.summarise(sorted(genes), assayed=assayed)
    out.views["fermentation"] = view
    if not assayed:
        out.unavailable.append(Unavailable(
            "A02", "Fermentation and cross-feeding", "missing_input",
            "the route definitions and the network are implemented; no gene panel ran "
            "for this sample",
            "a completed functional panel pass",
        ))
    elif view["n_supported"] == 0:
        out.unavailable.append(Unavailable(
            "A02.routes", "Supported fermentation routes", "missing_input",
            "no route had all of its required genes above the detection floor in this "
            "sample; each route says which genes were missing",
            "deeper sequencing, or the reference sets for the routes not yet covered",
        ))


#: Words that mark a capability's concept as a derived reading rather than
#: an organism the explorer could look up.
_NOT_A_TAXON: Final[frozenset[str]] = frozenset({
    "index", "abundance", "ratio", "score", "diversity", "richness", "count",
    "load", "profile", "balance", "fraction", "share", "coverage",
})


def _looks_like_a_taxon(concept: str) -> bool:
    """Whether this capability names an organism rather than a derived reading."""
    text = str(concept).strip()
    if not text or not text[:1].isupper():
        return False
    return not any(word in text.lower() for word in _NOT_A_TAXON)


def _searched_genes(results: Mapping[str, Any]) -> set[str]:
    """Every target gene this run searched for, detected or not.

    Needed to tell "looked for and not found" from "never looked for".
    A route whose genes no panel carries is `not_assayed`; calling it
    absent would report a gap in the assay as a property of the person.
    """
    out: set[str] = set()
    panels = results.get("panels") or []
    if not isinstance(panels, list):
        return out
    for panel in panels:
        for gene in panel.get("genes") or []:
            if str(gene.get("role")) == "target" and gene.get("gene"):
                out.add(str(gene["gene"]))
    return out


def _detected_genes(results: Mapping[str, Any]) -> tuple[set[str], bool]:
    """Target genes detected in this run, and whether any panel ran at all.

    Decoys are excluded: a panel searches for them in order to reject them,
    and treating one as present would build a route out of a negative
    control.
    """
    genes: set[str] = set()
    panels = results.get("panels") or []
    if not isinstance(panels, list):
        return genes, False
    for panel in panels:
        for gene in panel.get("genes") or []:
            if str(gene.get("role")) != "target":
                continue
            if str(gene.get("confidence")) in {"not detected", "absent"}:
                continue
            if int(gene.get("fragments") or 0) > 0 and gene.get("gene"):
                genes.add(str(gene["gene"]))
    return genes, bool(panels)


def _run_vitamins(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A03 - nine vitamins, with synthesis, salvage and uptake kept apart."""
    genes, assayed = _detected_genes(results)
    panels = {
        str(p.get("name")): p
        for p in (results.get("panels") or [])
        if isinstance(p, dict) and p.get("name")
    }
    view = VIT.dashboard(
        sorted(genes), existing_panels=panels, assayed=assayed,
        genes_searched=sorted(_searched_genes(results)),
    )
    out.views["vitamins"] = view
    unsearched = [
        v["vitamin_id"] for v in view["new"]
        if v["columns"]["synthesis"]["state"] == "not_assayed"
    ]
    incomplete = [
        v["vitamin_id"] for v in view["new"]
        if v["columns"]["synthesis"]["state"] in {"partial", "absent"}
    ]
    if unsearched:
        out.unavailable.append(Unavailable(
            "A03.new_vitamins", "Thiamine, niacin, pantothenate and B6 routes",
            "missing_input",
            "the routes for " + ", ".join(v.upper() for v in unsearched)
            + " are reconstructed and tested, and none of their genes is in the panels "
            "this run searched; that is not the same as their being absent from you",
            "the F06-F07 curated route reference sets",
        ))
    if incomplete:
        out.unavailable.append(Unavailable(
            "A03.synthesis", "Complete synthesis routes for the new vitamins",
            "missing_input",
            "no complete de-novo route was reconstructed for "
            + ", ".join(v.upper() for v in incomplete)
            + "; each route names which of its genes were missing",
            "deeper sequencing, or the remaining route references",
        ))


def _run_nitrogen(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A04 - protein, nitrogen and aromatic metabolism, kept apart."""
    detected, _assayed = _detected_genes(results)
    view = NIT.summarise(sorted(detected), searched=sorted(_searched_genes(results)))
    out.views["nitrogen"] = view
    unsearched = [
        m["module_id"] for m in view["modules"] if m["state"] == "not_assayed"
    ]
    if unsearched:
        out.unavailable.append(Unavailable(
            "A04.modules", "Protein and nitrogen modules", "missing_input",
            "the modules for " + ", ".join(unsearched)
            + " are defined and tested; none of their genes is in the panels this run "
            "searched, which is not the same as their being absent from you",
            "the F10-F14 curated reference sets",
        ))


def _run_biotransform(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A05, A06, A07 and A08 - named chemical steps and what they cannot claim."""
    detected, _ = _detected_genes(results)
    searched = sorted(_searched_genes(results))
    panels = {
        str(p.get("name")): p
        for p in (results.get("panels") or [])
        if isinstance(p, dict) and p.get("name")
    }
    view = BIO.summarise(sorted(detected), searched=searched, panels=panels)
    view["glp1"] = BIO.glp1_panel(sorted(detected), searched=searched)
    out.views["biotransformation"] = view

    # A matched negative control is a fact about the panel and has to surface.
    for group in view["by_feature"]:
        withheld = [
            step for step in group["steps"] if step["state"] == "negative_control_matched"
        ]
        if withheld:
            out.unavailable.append(Unavailable(
                f"{group['feature_id']}.control", group["feature"], "missing_input",
                "a declared negative control matched for "
                + ", ".join(s["label"] for s in withheld)
                + "; the step is withheld because the panel cannot currently tell the "
                "two enzymes apart in this sample",
                "residue-level discrimination over the diagnostic positions",
            ))
        unsearched = [s for s in group["steps"] if s["state"] == "not_assayed"]
        if unsearched and len(unsearched) == len(group["steps"]):
            out.unavailable.append(Unavailable(
                group["feature_id"], group["feature"], "missing_input",
                f"the {len(unsearched)} steps of this capability are defined and tested; "
                "none of their genes is in the panels this run searched",
                "the reference sets named in this capability's sources",
            ))


#: Which genes each functional reading requires, so it can be placed against
#: the cohort. Derived from the modules themselves rather than restated, so a
#: change to a route's requirement changes its score too.
def _reading_requirements(out: ExtensionResult) -> list[dict[str, Any]]:
    families: dict[str, list[str]] = {}
    for gene, fams in CAZY_FAMILY.items():
        for family in fams:
            families.setdefault(family, []).append(gene)

    readings: list[dict[str, Any]] = []
    for entry in (out.views.get("substrates") or {}).get("substrates") or []:
        genes = sorted({
            gene for family in entry.get("required_families") or ()
            for gene in families.get(str(family), ())
        })
        readings.append({"reading_id": entry.get("substrate_id"), "genes": genes})
    for route in (out.views.get("fermentation") or {}).get("routes") or []:
        readings.append({"reading_id": route.get("route_id"),
                         "genes": list(route.get("required_genes") or ())})
    for module in (out.views.get("nitrogen") or {}).get("modules") or []:
        readings.append({"reading_id": module.get("module_id"),
                         "genes": list(module.get("genes") or ())})
    for vitamin in (out.views.get("vitamins") or {}).get("new") or []:
        # The synthesis route is what the vitamin reading is about; salvage
        # and uptake are their own columns and are placed with the panel
        # entries that measure them, not folded in here.
        synthesis = [
            route for route in vitamin.get("routes") or []
            if route.get("column") == "synthesis"
        ]
        genes = sorted({
            gene
            for route in synthesis
            for branch in route.get("branches") or ()
            for gene in branch.get("genes") or ()
        } | {
            gene for route in synthesis for gene in route.get("coupling") or ()
        })
        readings.append({"reading_id": vitamin.get("vitamin_id"), "genes": genes})
    bio = out.views.get("biotransformation") or {}
    for feature in bio.get("by_feature") or []:
        for step in feature.get("steps") or []:
            readings.append({
                "reading_id": step.get("step_id"),
                "genes": list((step.get("requirement") or {}).get("genes") or ()),
            })
    for channel in (bio.get("glp1") or {}).get("channels") or []:
        readings.append({
            "reading_id": channel.get("step_id"),
            "genes": list((channel.get("requirement") or {}).get("genes") or ()),
        })
    for capacity in (out.views.get("additional_capacities") or {}).get("capacities") or []:
        readings.append({
            "reading_id": capacity.get("capacity_id"),
            "genes": list((capacity.get("requirement") or {}).get("genes") or ()),
        })
    return [r for r in readings if r.get("reading_id")]


def _run_reading_scores(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """Place every functional reading against the cohort, so each gets a bar.

    These arrived in the report as words - "class-level only", "part of it" -
    beside forty-two rows that each carried a percentile and a slider, which
    read as missing data. The distributions to compare against already
    existed, one per panel entry over the same 91 samples; what was missing
    was the join from a requirement to the entries that measure it.
    """
    path = REFS_ROOT / "reference_ranges.json"
    if not path.is_file():
        return
    ranges = json.loads(path.read_text(encoding="utf-8"))
    searched = sorted(_searched_genes(results))
    requirements = _defer_to_panel_curation(_reading_requirements(out), results)
    scores = RSC.score_all(
        requirements, results=results, ranges=ranges, searched=searched)
    out.views["reading_scores"] = {rid: s.to_json() for rid, s in scores.items()}


def _defer_to_panel_curation(
    requirements: Sequence[Mapping[str, Any]], results: Mapping[str, Any]
) -> list[dict[str, Any]]:
    """Let a panel's curated gene choice govern a reading of the same thing.

    A panel that aggregates one gene out of an operon has had that choice
    made deliberately: the urolithin panel counts *ucdO* alone, because
    that subunit is what performs the 9-dehydroxylation, and a stray
    *ucdC* fragment without it makes no urolithin A. The reading module
    lists the whole operon, so summing its three genes put the reading at
    the 46th percentile in a sample whose panel correctly measured zero.

    Where a reading's genes cover everything its panel aggregates, the two
    are the same quantity and the panel's narrower set wins. Where they
    are a strict subset - one fibre inside the carbohydrate panel - the
    reading is asking a narrower question and keeps its own.
    """
    from openbiota.pdfsynbiotic import owner_of  # noqa: PLC0415

    by_gene: dict[str, dict[str, set[str]]] = {}
    for panel in results.get("panels") or []:
        if not isinstance(panel, dict):
            continue
        aggregate = set(panel.get("aggregate_from") or ())
        if not aggregate:
            continue
        by_gene[str(panel.get("name"))] = {
            "aggregate": aggregate,
            "genes": {
                str(g.get("gene")) for g in panel.get("genes") or ()
                if g.get("entry_id") in aggregate and g.get("gene")
            },
            "all": {
                str(g.get("gene")) for g in panel.get("genes") or () if g.get("gene")
            },
        }

    out: list[dict[str, Any]] = []
    for requirement in requirements:
        row = dict(requirement)
        spec = by_gene.get(str(owner_of(str(row.get("reading_id") or "")) or ""))
        genes = {str(g) for g in row.get("genes") or ()}
        if spec and genes and spec["aggregate"] and genes >= spec["genes"] and (
            genes <= spec["all"]
        ):
            row["genes"] = sorted(spec["genes"])
        out.append(row)
    return out


def _run_contexts(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A13 - the fifteen symptom and condition cards of section 16.3.

    Runs after the readings it links to, so that a card can say which of
    them this particular run actually produced.
    """
    out.views["contexts"] = CX.resolve(results, out.views)


def _run_capacities(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """Section 5.8 - the additional capacities, and the drug-metabolism cards."""
    detected, _ = _detected_genes(results)
    searched = sorted(_searched_genes(results))
    panels = {
        str(p.get("name")): p
        for p in (results.get("panels") or [])
        if isinstance(p, dict) and p.get("name")
    }
    out.views["additional_capacities"] = CAP.summarise(
        sorted(detected), searched=searched, panels=panels)

    # A08's beta-glucuronidase substrate detail, which reports itself as
    # withheld unless its class assignment agrees with the published
    # structures. That check lives in the module; this only records the
    # outcome where the coverage ledger can see it.
    # §9.1's required P9 sequence panel. The anchors are verified against
    # the hashes the specification pins before the reading is published, so
    # a changed record is an error rather than a quietly different answer.
    anchors: dict[str, Any] = {}
    try:
        anchors = P9.verify_anchors(cache_dir=REFS_ROOT / "cache")
    except Exception as exc:  # noqa: BLE001 - recorded, not fatal
        out.unavailable.append(Unavailable(
            "A07.p9_anchors", "P9 sequence anchors", "missing_input",
            f"the pinned P9 records could not be verified: {type(exc).__name__}: {exc}"[:300],
            "network access to UniProt and ENA, or a cached anchor record",
        ))
    out.views["p9"] = P9.reading(results, anchors=anchors)

    view = GUS.analyse(results, refs_root=REFS_ROOT)
    if view is not None:
        out.views["gus_substrate_classes"] = view
        if not view.get("available"):
            out.unavailable.append(Unavailable(
                "A08.gus_substrate", "Beta-glucuronidase substrate classes",
                "scientifically_unresolved",
                str(view.get("reason"))[:600],
                view.get("what_would_unlock_it"),
            ))


#: Which extension state a functional reading's own state maps to. The
#: distinctions the modules make are finer than the schema's, and the
#: mapping keeps the important one: searched-and-absent is a measurement,
#: never-searched is not.
_FUNCTIONAL_STATE: Final[Mapping[str, str]] = {
    "substrate_specific": "measured",
    "class_level_only": "partial",
    "supported": "measured",
    "complete": "measured",
    "present": "measured",
    "partial": "partial",
    "insufficient": "not_detected_above_assay_threshold",
    "absent": "not_detected_above_assay_threshold",
    "not_assayed": "not_assayed",
    "not_applicable": "not_applicable",
    "negative_control_matched": "unsupported_by_assay",
}


#: The panel each newly reconstructed vitamin is measured by, so its metric
#: can carry the panel's own copies-per-100-genomes rather than a bare state.
_VITAMIN_PANEL: Final[Mapping[str, str]] = {
    "b1": "thiamine", "b3": "niacin", "b5": "pantothenate", "b6": "b6",
}


def _functional_metrics(
    out: ExtensionResult, panels: Mapping[str, Mapping[str, Any]] | None = None
) -> list[ExtensionMetric]:
    """One metric per functional reading, so each reaches the index and the export.

    The fragment count is the value: it is what was measured. The
    specificity or route state travels in `extra`, because "we found GH13"
    and "that does not establish resistant-starch use" are two different
    facts and only one of them is a number.
    """
    metrics: list[ExtensionMetric] = []

    def emit(metric_id: str, label: str, feature: str, state: str, value: float | None,
             unit: str, limits: tuple[str, ...], extra: dict[str, Any]) -> None:
        mapped = _FUNCTIONAL_STATE.get(state, "not_assayed")
        if mapped in {"measured", "partial"} and value is None:
            # A state that promises a number, with no number behind it, is a
            # reading that was made and could not be quantified.
            mapped = "insufficient_coverage"
        numeric = mapped in {"measured", "partial"}
        metrics.append(ExtensionMetric(
            metric_id=metric_id, label=label, kind="genetic_capacity",
            state=mapped, method_id=f"{feature}.functional/1.0",
            value=float(value) if (numeric and value is not None) else None,
            unit=unit if (numeric and value is not None) else None,
            denominator="fragments accepted against this panel's references"
            if (numeric and value is not None) else None,
            feature_id=feature, group="functions",
            analytical_confidence="supported" if mapped == "measured" else "provisional",
            evidence_maturity="genomic_prediction",
            limitations=limits or ("Genetic capacity, not a measured amount.",),
            input_fingerprint=fingerprint(feature, metric_id, state, value),
            extra={"reading_state": state, **extra},
        ))

    for entry in (out.views.get("substrates") or {}).get("substrates") or []:
        emit(f"ext083.{entry['substrate_id']}", entry["label"], "A01",
             str(entry["specificity"]), float(entry.get("fragments") or 0), "fragments",
             tuple(entry.get("limitations") or ()),
             {"families_present": entry.get("required_present")})

    for route in (out.views.get("fermentation") or {}).get("routes") or []:
        emit(f"ext083.fermentation.{route['route_id']}", route["label"], "A02",
             str(route["state"]), float(len(route.get("required_present") or [])),
             "genes present", tuple(route.get("limitations") or ()),
             {"genes_missing": route.get("required_missing")})

    for vitamin in (out.views.get("vitamins") or {}).get("new") or []:
        columns = vitamin.get("columns") or {}
        panel = (panels or {}).get(_VITAMIN_PANEL.get(str(vitamin.get("vitamin_id")), ""), {})
        emit(vitamin["metric_id"], vitamin["label"], "A03",
             str(columns.get("synthesis", {}).get("state")),
             panel.get("copies_per_100_genomes"), "copies per 100 bacterial genomes",
             tuple(vitamin.get("limitations") or ()),
             {"columns": {k: v.get("state") for k, v in columns.items()}})

    for module in (out.views.get("nitrogen") or {}).get("modules") or []:
        emit(f"ext083.nitrogen.{module['module_id']}", module["label"], "A04",
             str(module["state"]), float(len(module.get("informative_genes_present") or [])),
             "informative genes present", tuple(module.get("limitations") or ()),
             {"genes_missing": module.get("genes_missing")})

    bio = out.views.get("biotransformation") or {}
    for group in bio.get("by_feature") or []:
        for step in group.get("steps") or []:
            requirement = step.get("requirement") or {}
            emit(f"ext083.{step['step_id']}", step["label"], str(group["feature_id"]),
                 str(step["state"]), float(len(requirement.get("present") or [])),
                 "genes present", tuple(step.get("limitations") or ()),
                 {"genes_missing": requirement.get("missing")})
    for channel in (bio.get("glp1") or {}).get("channels") or []:
        requirement = channel.get("requirement") or {}
        emit(f"ext083.{channel['step_id']}", channel["label"], "A07",
             str(channel["state"]), float(len(requirement.get("present") or [])),
             "genes present", tuple(channel.get("limitations") or ()),
             {"genes_missing": requirement.get("missing")})

    # Section 5.8. These carry their panel's copies-per-100-genomes rather
    # than a gene count, because each has a panel of its own behind it.
    for capacity in (out.views.get("additional_capacities") or {}).get("capacities") or []:
        requirement = capacity.get("requirement") or {}
        emit(f"ext083.{capacity['capacity_id']}", capacity["label"], "5.8",
             str(capacity["state"]), capacity.get("copies_per_100_genomes"),
             "copies per 100 bacterial genomes",
             (capacity.get("does_not_establish") or "",),
             {"genes_missing": requirement.get("missing"),
              "evidence_level": capacity.get("evidence_level"),
              "context": capacity.get("context")})
    return metrics


def _run_planner(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A12 - one ranked view across every per-finding evidence card.

    Consolidation, not replacement: each option points back at the cards it
    came from, and nothing appears here that the registry does not already
    hold.
    """
    evidence = results.get("findings_and_evidence") or {}
    cards = evidence.get("evidence_summaries") or []
    if not cards:
        out.unavailable.append(Unavailable(
            "A12", "Consolidated action planner", "missing_input",
            "no evidence cards were produced for this sample, so there is nothing to "
            "consolidate",
            "a completed intervention-evidence pass",
        ))
        return

    identities: dict[str, Any] = {}
    identities_file = Path(__file__).resolve().parents[2] / "interventions" / "identities.yaml"
    if identities_file.is_file():
        import yaml  # noqa: PLC0415

        loaded = yaml.safe_load(identities_file.read_text(encoding="utf-8")) or {}
        identities = {
            str(record["intervention_id"]): record
            for record in loaded.get("interventions") or []
            if record.get("intervention_id")
        }

    flagged = [
        str(v.get("display") or str(v.get("species") or "").replace("_", " "))
        for v in (results.get("organism_verdicts") or {}).get("verdicts") or []
        if v.get("flag")
    ]
    # A supplied allergy or intolerance excludes an option from *this*
    # person's start-here list and from nothing else. An unknown fact is
    # not an exclusion, which is why only supplied values are read.
    supplied = (results.get("subject_context") or {}).get("supplied") or {}
    exclusions: dict[str, PLAN.Exclusion] = {}
    for key, value in supplied.items():
        if not value or "allerg" not in str(key).lower():
            continue
        for intervention_id, record in identities.items():
            if str(value).lower() in str(record.get("display_name", "")).lower():
                exclusions[intervention_id] = PLAN.Exclusion(
                    reason="declared_allergy",
                    detail=f"you reported an allergy to {value}",
                    from_supplied_fact=True,
                )

    plan = PLAN.build_plan(
        cards, sample_id=out.sample, flagged_findings=flagged,
        excluded=exclusions, identities=identities,
    )
    out.views["planner"] = plan.to_json()
    coverage = plan.coverage()
    if coverage["n_flagged_findings_addressed"] == 0 and flagged:
        out.unavailable.append(Unavailable(
            "A12.flagged", "Options for the organisms flagged here", "missing_input",
            f"{len(flagged)} organism(s) are flagged in this sample and the intervention "
            "registry holds no studied option that names any of them",
            "evidence for these specific organisms",
        ))


def _run_longitudinal(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A13 - the view across this participant's samples, or why there is none.

    A run analyses one sample. Grouping it with earlier ones needs a
    participant identifier and a collection date, and inventing either would
    manufacture a timeline. Where they are absent this reports the
    single-timepoint state, which §6.4 requires to stay visible rather than
    be suppressed.
    """
    inventory = (results.get("organism_inventory") or {}).get("organisms") or []
    abundances = {
        str(o.get("species") or ""): float(o.get("percent") or 0.0)
        for o in inventory if o.get("percent")
    }
    meta = results.get("run") or results.get("meta") or {}
    # Recorded exposures pick the follow-up template. Only *supplied* facts
    # count: an assumed "no antibiotics" is not a recorded exposure, and
    # letting it choose a template would be an assumption steering advice.
    context = results.get("subject_context") or {}
    supplied = context.get("supplied") or {}
    events = tuple(
        str(key).replace("_", " ")
        for key, value in supplied.items()
        if value not in (False, None, "", [])
    ) + tuple(str(m) for m in context.get("unclassified_medications") or [])

    point = LON.TimePoint(
        fingerprint=LON.Fingerprint(
            sample_id=out.sample,
            participant_id=meta.get("participant_id"),
            specimen_id=meta.get("specimen_id"),
            fields={
                "profiler": (results.get("taxonomy") or {}).get("profiler"),
                "specimen_type": "stool",
                "host_filter": (results.get("sequencing_quality") or {}).get("host_filter"),
                "read_depth": (results.get("sequencing_quality") or {}).get("usable_nonhost_pairs"),
                "collection_date": meta.get("collection_date"),
            },
        ),
        collected_at=LON.as_date(meta.get("collection_date")),
        abundances=abundances,
        events=events,
    )
    series = LON.Series.build(str(meta.get("participant_id") or out.sample), [point])
    if series.n_timepoints > 1:
        view = LON.summarise(series)
    else:
        view = LON.not_enough_history(out.sample, n_available=series.n_timepoints)
        view["follow_up"] = LON.follow_up(events)
        view["recorded_exposures"] = list(events)
    view["fingerprint"] = point.fingerprint.to_json()
    view["comparability_fields"] = {
        "blocking": list(LON.BLOCKING_FIELDS), "context": list(LON.CONTEXT_FIELDS),
    }
    out.views["longitudinal"] = view
    if view.get("state") == "single_timepoint":
        out.unavailable.append(Unavailable(
            "A13.trend", "Change over time", "missing_input",
            "only one sample is on file for this participant, so nothing can be compared "
            "with it yet",
            "a second sample collected and processed the same way",
        ))
        out.unavailable.append(Unavailable(
            "A13.recovery", "Return toward a pre-event community", "missing_input",
            f"a recovery estimate needs a recorded perturbation and at least "
            f"{LON.MIN_BASELINE_SAMPLES} eligible samples before it",
            f"{LON.MIN_BASELINE_SAMPLES} baseline samples and a recorded event",
        ))


def _run_inputs(results: Mapping[str, Any], out: ExtensionResult) -> None:
    """A14 - the register of everything supplied, assumed and absent.

    Built last, because it reads the other views: the diet a scenario
    assumed and the cohort a percentile came from are inputs to this report
    exactly as much as anything a person typed.
    """
    attached = dict(results)
    attached["extension"] = {"views": out.views}
    register = INP.build_register(attached, sample_id=out.sample)
    out.views["input_register"] = register.to_json()
    for record in register.missing:
        if record.category != "laboratory":
            continue
        out.unavailable.append(Unavailable(
            "A14.laboratory", "External laboratory results", "external_assay",
            record.limitations[0] if record.limitations else "not supplied",
            "an imported laboratory result",
        ))


def _run_owed(out: ExtensionResult) -> None:
    """Name every capability this release still owes, one by one.

    Section 15.7: an engineering TODO may not be disguised as a scientific
    limitation, so each outstanding feature says `not_implemented` in its own
    row rather than being absent from the output.
    """
    registry = REG.load()
    delivered = {
        "A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08",
        "A09", "A10", "A11", "A12", "A13", "A14", "A15", "A16",
    }
    for capability in registry.features:
        if capability.capability_id in delivered or capability.implemented:
            continue
        out.unavailable.append(Unavailable(
            capability.capability_id, capability.concept, "not_implemented",
            capability.binding[:240],
            "implementation of this capability",
        ))


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #


def analyze(
    results: Mapping[str, Any],
    *,
    mode: str = "cached",
    sample: str | None = None,
    simulate: str = "auto",
    lp_workers: int | None = None,
    progress: Any = None,
) -> ExtensionResult:
    """Build every extension measurement this sample's data supports."""
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}; expected one of {sorted(MODES)}")
    started = time.monotonic()
    out = ExtensionResult(sample=str(sample or results.get("sample") or "unknown"), mode=mode)

    for runner, capability in (
        (_run_ecology, "A09"),
        (_run_explorer, "A10"),
        (_run_resistance, "A11"),
    ):
        try:
            runner(results, out)
        except Exception as exc:  # noqa: BLE001 - one capability cannot fail the run
            out.unavailable.append(Unavailable(
                capability, capability, "not_implemented",
                f"{type(exc).__name__}: {exc}"[:300],
                "a fix to this capability",
            ))
    try:
        _run_simulation(results, out, mode=mode, simulate=simulate,
                        lp_workers=lp_workers, progress=progress)
    except Exception as exc:  # noqa: BLE001
        out.unavailable.append(Unavailable(
            "A16", "Personalized synbiotic options", "not_implemented",
            f"{type(exc).__name__}: {exc}"[:300], None,
        ))
    for runner, capability in (
        (_run_substrates, "A01"), (_run_fermentation, "A02"), (_run_vitamins, "A03"),
        (_run_nitrogen, "A04"), (_run_biotransform, "A05/A06/A07/A08"),
        (_run_capacities, "5.8"),
        (_run_planner, "A12"), (_run_longitudinal, "A13"),
        (_run_reading_scores, "12.1.scores"),
        (_run_contexts, "A13.contexts"),
    ):
        try:
            runner(results, out)
        except Exception as exc:  # noqa: BLE001 - one capability cannot fail the run
            out.unavailable.append(Unavailable(
                capability, capability, "not_implemented",
                f"{type(exc).__name__}: {exc}"[:300], "a fix to this capability",
            ))
    out.metrics.extend(_functional_metrics(out, {
        str(p.get("name")): p for p in (results.get("panels") or [])
        if isinstance(p, dict) and p.get("name")
    }))
    _run_inputs(results, out)
    _run_owed(out)

    out.elapsed_s = time.monotonic() - started
    out.manifest = {
        "schema_version": "openbiota.extension-manifest/1.0",
        "spec_version": SPEC_VERSION,
        "mode": mode,
        "inputs": _input_manifest(results),
        "software": _software_manifest(),
        "platform": f"{platform.system()} {platform.machine()}",
        "python": sys.version.split()[0],
        "parameter_fingerprint": fingerprint(mode, SPEC_VERSION, _software_manifest()),
        "n_metrics": len(out.metrics),
        "n_unavailable": len(out.unavailable),
        "engineering_debt": sorted(
            u.capability_id for u in out.unavailable if u.reason in REG.ENGINEERING_DEBT
        ),
    }
    return out


def coverage_report(result: ExtensionResult) -> dict[str, Any]:
    """`capability_coverage.json` for this sample.

    Joins the normative manifest with what actually happened, so a reader can
    see both what the release owes and what this run produced.
    """
    registry = REG.load()
    payload = registry.coverage(sample=result.sample)
    payload["this_run"] = {
        "mode": result.mode,
        "n_metrics": len(result.metrics),
        "metric_ids": result.metric_ids(),
        "unavailable": [u.to_json() for u in result.unavailable],
        "by_reason": {
            reason: sum(1 for u in result.unavailable if u.reason == reason)
            for reason in sorted({u.reason for u in result.unavailable})
        },
    }
    payload["sources"] = SRC.load().audit()
    return payload


def write_outputs(
    result: ExtensionResult, out_dir: Path, *, all_metrics_tsv: bool = True
) -> dict[str, Path]:
    """Write the section 13.1 artefacts. No hidden truncation anywhere."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    ext = out_dir / "results.ext083.json"
    ext.write_text(json.dumps(result.to_json(), indent=1) + "\n", encoding="utf-8")
    written["results.ext083.json"] = ext

    manifest = out_dir / "extension_manifest.json"
    manifest.write_text(json.dumps(result.manifest, indent=1) + "\n", encoding="utf-8")
    written["extension_manifest.json"] = manifest

    coverage = out_dir / "capability_coverage.json"
    coverage.write_text(json.dumps(coverage_report(result), indent=1) + "\n", encoding="utf-8")
    written["capability_coverage.json"] = coverage

    if all_metrics_tsv:
        tsv = out_dir / "all_metrics.tsv"
        columns = (
            "metric_id", "label", "group", "feature_id", "kind", "state", "value", "unit",
            "denominator", "direction", "reference_percentile", "reference_state",
            "analytical_confidence", "method_id", "input_fingerprint",
        )
        lines = ["\t".join(columns)]
        for metric in result.metrics:
            row = metric.to_json()
            lines.append("\t".join(
                "" if row.get(c) is None else str(row.get(c)).replace("\t", " ")
                for c in columns
            ))
        tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
        written["all_metrics.tsv"] = tsv

    return written


def attach(results: dict[str, Any], result: ExtensionResult) -> None:
    """Attach the extension to a result object, under its own key only.

    The one line that keeps the preservation contract: everything new lives
    under `extension`, and no protected object is touched.
    """
    from openbiota.extension.preservation import PROTECTED_KEYS

    if RESULT_KEY in PROTECTED_KEYS:  # pragma: no cover - guards a future edit
        raise RuntimeError("the extension key must never be a protected object")
    results[RESULT_KEY] = result.to_json()


__all__ = [
    "MODES",
    "RESULT_KEY",
    "ExtensionResult",
    "Unavailable",
    "analyze",
    "attach",
    "coverage_report",
    "write_outputs",
]
