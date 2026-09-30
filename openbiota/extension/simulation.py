"""A16 model-assisted scenarios: MICOM + AGORA2 community simulation.

BUILD_SPEC_v0.8.3 sections 11.6 and 11.7. The question this answers is the one
a reader actually asks about a probiotic and a prebiotic: *is the combination
worth more than either part on its own?* The engine builds four communities -
nothing added, the probiotic, the substrate, and both - and reports what
changes between them.

What it produces is an **experimental prediction from a metabolic model**, and
the code is built so that cannot be overstated:

* **Predicted flux is not a concentration.** Not stool butyrate, not absorbed
  exposure, not a dose, not a symptom. Every emitted metric says so.
* **Raw flux and flux per unit growth are different metrics.** They rank
  candidates differently - in the published fixture psyllium wins on raw
  butyrate and maltodextrin wins per unit of growth - so they are never
  interchangeable and never share a unit.
* **An infeasible solve is unavailable, not zero.** A nonoptimal status
  returns no number with the reason attached. Zero production is a result; a
  failed optimisation is not.
* **Mass is checked before anything is normalised.** The published notebook
  builds a perturbation summing to 1.15; dividing that away would silently
  change every abundance, so it fails instead.
* **Coverage is reported before renormalisation.** How much of the community
  had a model at all, and what the abundance cutoff discarded, are stated
  first; a model covering 60% of the reads cannot present itself as the
  community.
* **A published-output replay is not a solved model.** The two have different
  execution states, and neither is a measured outcome in a person.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.extension import lpsolve as LP
from openbiota.extension.schema import ExtensionMetric, fingerprint
from openbiota.extension.strict_ooxml import StrictWorkbook

CONFIG_FILE: Final = Path(__file__).resolve().parents[2] / "extension" / "simulation_config.yaml"
METHOD_SCENARIO: Final = "ext083.micom_scenario/1.0"
METHOD_REPLAY: Final = "ext083.micom_published_replay/1.0"

#: How a scenario's numbers came to exist. `published_replay` reads released
#: predictions; `solved_model` ran an optimisation here. Section 11.7 requires
#: these to be distinguishable, and neither is a measurement in a person.
EXECUTION_STATES: Final[frozenset[str]] = frozenset({
    "solved_model", "published_replay", "not_run", "failed",
})

#: Solver outcomes that yield a number. Everything else is unavailable.
OPTIMAL_STATUSES: Final[frozenset[str]] = frozenset({"optimal"})

#: How far two solves of the *same* community may disagree before its numbers
#: stop being reportable, as a fraction of the larger value.
#:
#: This is now a guard rather than a limitation. The first implementation
#: reported no numbers at all, because `cooperative_tradeoff` returned
#: "optimal" twice with butyrate fluxes of 19.50 and 12.50: its L2
#: regularisation makes the *growth rates* unique but leaves the flux vector
#: degenerate, and a single solve picks an arbitrary point among many.
#:
#: The fix is to stop asking a degenerate question. `solve_arm` maximises
#: community growth as a linear programme, pins it, and then bounds the target
#: exchange with two more linear programmes. An LP's optimal *value* is unique
#: even when its solution vector is not, and optlang's hybrid interface runs
#: them on a deterministic simplex. Three identical runs of the same community
#: now return `gmax=0.053888799182` and a butyrate maximum of
#: `139.864541107` to every digit - and take two seconds instead of seven
#: minutes, because no quadratic programme is solved at all.
#:
#: The guard remains so that a future change of solver, medium or model cannot
#: quietly reintroduce the problem.
REPRODUCIBILITY_TOLERANCE: Final = 1e-6

#: Linear-programming method to use. Simplex returns an exact vertex optimum
#: and is reproducible; the interior-point default converges to a tolerance,
#: which is what made the first implementation's numbers move between runs.
LP_METHOD: Final = LP.LP_METHOD

#: Solvers whose LP optimum is exact and reproducible. The open hybrid
#: interface qualifies once its LP method is simplex, which is why no
#: commercial licence is needed to quote a number.
DETERMINISTIC_SOLVERS: Final[frozenset[str]] = frozenset({"gurobi", "cplex", "osqp", "hybrid"})


class SimulationError(ValueError):
    """A scenario that would misrepresent a model as a measurement."""


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #


@dataclass(slots=True, frozen=True)
class Substrate:
    """One modelled substrate and its total exchange cap."""

    substrate_id: str
    label: str
    exchange: str
    cap: float
    food_context: str | None = None
    proxy_for: str | None = None

    @property
    def is_proxy(self) -> bool:
        return bool(self.proxy_for)

    def to_json(self) -> dict[str, Any]:
        return {
            "substrate_id": self.substrate_id,
            "label": self.label,
            "exchange_reaction": self.exchange,
            "cap_mmol_per_gDW_per_h": self.cap,
            "cap_semantics": "total cap on this exchange, not an addition to the medium's bound",
            "food_context": self.food_context,
            "is_model_proxy": self.is_proxy,
            "proxy_note": self.proxy_for,
        }


@dataclass(slots=True, frozen=True)
class Protocol:
    """The pinned, corrected paper protocol."""

    protocol_id: str
    micom_version: str
    model_asset: str
    model_asset_md5: str
    cutoff: float
    tradeoff: float
    strategy: str
    resident_fraction: float
    added_fraction_per_species: float
    required_total_abundance: float
    mass_tolerance: float
    normalize_input_abundance: bool
    merge_existing_and_added_species: str
    growth_threshold_methods: float
    growth_threshold_figure_caption: float
    model_source_id: str

    def to_json(self) -> dict[str, Any]:
        return {
            "protocol_id": self.protocol_id,
            "micom_version": self.micom_version,
            "model_asset": self.model_asset,
            "model_asset_md5": self.model_asset_md5,
            "model_source_id": self.model_source_id,
            "cutoff": self.cutoff,
            "tradeoff": self.tradeoff,
            "strategy": self.strategy,
            "resident_fraction": self.resident_fraction,
            "added_fraction_per_species": self.added_fraction_per_species,
            "required_total_abundance": self.required_total_abundance,
            "normalize_input_abundance": self.normalize_input_abundance,
            "merge_existing_and_added_species": self.merge_existing_and_added_species,
            "growth_threshold_discrepancy": {
                "methods_code": self.growth_threshold_methods,
                "figure_caption": self.growth_threshold_figure_caption,
                "resolution": (
                    "Preserved, not repaired. Any threshold sensitivity is research "
                    "output, never an engraftment probability."
                ),
            },
        }


@dataclass(slots=True)
class SimulationConfig:
    """The whole pinned configuration, loaded once."""

    protocol: Protocol
    added_species: tuple[str, ...]
    substrates: tuple[Substrate, ...]
    media: tuple[dict[str, Any], ...]
    replay: dict[str, Any]
    contrast_fixture: dict[str, Any]

    @classmethod
    def load(cls, path: Path | None = None) -> SimulationConfig:
        raw = yaml.safe_load((path or CONFIG_FILE).read_text(encoding="utf-8")) or {}
        p = raw["protocol"]
        return cls(
            protocol=Protocol(
                protocol_id=str(p["protocol_id"]),
                micom_version=str(p["micom_version"]),
                model_asset=str(p["model_asset"]),
                model_asset_md5=str(p["model_asset_md5"]),
                cutoff=float(p["cutoff"]),
                tradeoff=float(p["tradeoff"]),
                strategy=str(p["strategy"]),
                resident_fraction=float(p["resident_fraction"]),
                added_fraction_per_species=float(p["added_fraction_per_species"]),
                required_total_abundance=float(p["required_total_abundance"]),
                mass_tolerance=float(p["mass_tolerance"]),
                normalize_input_abundance=bool(p["normalize_input_abundance"]),
                merge_existing_and_added_species=str(p["merge_existing_and_added_species"]),
                growth_threshold_methods=float(p["growth_threshold_methods"]),
                growth_threshold_figure_caption=float(p["growth_threshold_figure_caption"]),
                model_source_id=str(p.get("model_source_id") or "S11"),
            ),
            added_species=tuple(raw.get("added_species") or ()),
            substrates=tuple(
                Substrate(
                    substrate_id=str(s["id"]), label=str(s["label"]),
                    exchange=str(s["exchange"]), cap=float(s["cap"]),
                    food_context=s.get("food_context"), proxy_for=s.get("proxy_for"),
                )
                for s in raw.get("substrates") or ()
            ),
            media=tuple(raw.get("media") or ()),
            replay=dict(raw.get("replay_fixture") or {}),
            contrast_fixture=dict(raw.get("contrast_fixture") or {}),
        )

    def substrate(self, substrate_id: str) -> Substrate:
        for s in self.substrates:
            if s.substrate_id == substrate_id:
                return s
        raise SimulationError(f"unknown substrate {substrate_id!r}")


# --------------------------------------------------------------------------- #
# 11.6 — the perturbation recipe and its mass check
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class Perturbation:
    """A community abundance recipe whose mass has actually been checked."""

    abundances: dict[str, float]
    resident_fraction: float
    added: tuple[str, ...]
    added_fraction_each: float
    total_before_normalisation: float
    residual_normalised: bool
    note: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "n_taxa": len(self.abundances),
            "resident_fraction": self.resident_fraction,
            "added_species": list(self.added),
            "added_fraction_each": self.added_fraction_each,
            "total_before_normalisation": self.total_before_normalisation,
            "floating_point_residual_normalised": self.residual_normalised,
            "note": self.note,
        }


def build_perturbation(
    baseline: Mapping[str, float],
    *,
    protocol: Protocol,
    added: Sequence[str] = (),
    resident_fraction: float | None = None,
    added_fraction_each: float | None = None,
) -> Perturbation:
    """Build an added-species recipe, refusing one whose mass is wrong.

    The order matters and is the point of AT101. The recipe is assembled,
    then its total is asserted against 1.0 within an absolute tolerance of
    1e-9 **before** anything is normalised. Only a floating-point residual
    inside that tolerance is then divided out.

    The published notebook uses resident .95 with five .04 additions, which
    totals 1.15. Normalising first would quietly rescale every resident
    abundance by 1/1.15 and produce a plausible-looking run from an invalid
    recipe, so this raises instead.
    """
    resident = protocol.resident_fraction if resident_fraction is None else resident_fraction
    each = (
        protocol.added_fraction_per_species
        if added_fraction_each is None else added_fraction_each
    )
    if resident < 0 or each < 0:
        raise SimulationError("negative abundance fractions are not a recipe")

    base_total = sum(float(v) for v in baseline.values())
    if base_total <= 0:
        raise SimulationError("baseline community has zero total abundance")
    if any(float(v) < 0 for v in baseline.values()):
        raise SimulationError("a negative baseline abundance cannot be modelled")
    if any(not math.isfinite(float(v)) for v in baseline.values()):
        raise SimulationError("a non-finite baseline abundance cannot be modelled")

    # Residents scaled to the resident fraction; true aliases merged once.
    recipe: dict[str, float] = {}
    for name, value in baseline.items():
        scaled = resident * float(value) / base_total
        recipe[name] = recipe.get(name, 0.0) + scaled
    for name in added:
        recipe[name] = recipe.get(name, 0.0) + each

    total = sum(recipe.values())
    if not math.isclose(
        total, protocol.required_total_abundance, rel_tol=0.0, abs_tol=protocol.mass_tolerance
    ):
        raise SimulationError(
            f"perturbation mass {total!r} does not equal "
            f"{protocol.required_total_abundance} within {protocol.mass_tolerance}. "
            f"resident {resident} plus {len(added)} x {each} = "
            f"{resident + len(added) * each}. An invalid recipe cannot be repaired by "
            "dividing its total away: that rescales every resident abundance."
        )
    residual = abs(total - protocol.required_total_abundance) > 0.0
    if residual:
        recipe = {k: v / total for k, v in recipe.items()}
    return Perturbation(
        abundances=recipe,
        resident_fraction=resident,
        added=tuple(added),
        added_fraction_each=each,
        total_before_normalisation=total,
        residual_normalised=residual,
        note=("merged an added species that was already resident: "
              + ", ".join(sorted(set(added) & set(baseline)))) if set(added) & set(baseline) else None,
    )


# --------------------------------------------------------------------------- #
# 11.6 — model coverage, reported before renormalisation
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class Coverage:
    """How much of the community the model library could represent."""

    n_input_taxa: int
    n_mapped_taxa: int
    abundance_mapped: float
    abundance_unmapped: float
    abundance_below_cutoff: float
    denominator: str
    unmapped_examples: tuple[str, ...]
    cutoff: float

    @property
    def mapped_fraction(self) -> float | None:
        total = self.abundance_mapped + self.abundance_unmapped + self.abundance_below_cutoff
        return (self.abundance_mapped / total) if total > 0 else None

    def to_json(self) -> dict[str, Any]:
        return {
            "n_input_taxa": self.n_input_taxa,
            "n_mapped_taxa": self.n_mapped_taxa,
            "abundance_mapped": self.abundance_mapped,
            "abundance_unmapped": self.abundance_unmapped,
            "abundance_below_cutoff": self.abundance_below_cutoff,
            "mapped_fraction_of_input": self.mapped_fraction,
            "denominator": self.denominator,
            "cutoff": self.cutoff,
            "largest_unmapped": list(self.unmapped_examples),
            "limitations": [
                "Coverage is of the bacterial abundance offered to the model. Fungal, "
                "viral and all-read denominators are different numbers and are not mixed in.",
                "High abundance coverage does not prove complete functional coverage: a "
                "mapped species model may still lack the pathway in question.",
            ],
        }


def map_to_models(
    abundances: Mapping[str, float],
    available: Iterable[str],
    *,
    cutoff: float,
    denominator: str = "bacterial_abundance_offered_to_the_model",
) -> tuple[dict[str, float], Coverage]:
    """Keep the taxa the library models, and report the loss first.

    Returns the mapped subset **unnormalised** together with the coverage
    record, so a caller cannot report a renormalised community without having
    the loss in hand. An unmapped organism stays unmapped: it is never
    reassigned to a convenient neighbour's model (AT102).
    """
    modelled = {str(m) for m in available}
    mapped: dict[str, float] = {}
    unmapped: dict[str, float] = {}
    below: float = 0.0
    total = sum(float(v) for v in abundances.values()) or 1.0
    for name, value in abundances.items():
        share = float(value) / total
        if share < cutoff:
            below += float(value)
            continue
        if name in modelled:
            mapped[name] = float(value)
        else:
            unmapped[name] = float(value)
    coverage = Coverage(
        n_input_taxa=len(abundances),
        n_mapped_taxa=len(mapped),
        abundance_mapped=sum(mapped.values()),
        abundance_unmapped=sum(unmapped.values()),
        abundance_below_cutoff=below,
        denominator=denominator,
        unmapped_examples=tuple(
            k for k, _ in sorted(unmapped.items(), key=lambda kv: -kv[1])[:8]
        ),
        cutoff=cutoff,
    )
    return mapped, coverage


# --------------------------------------------------------------------------- #
# 11.6 — the four-arm contrasts
# --------------------------------------------------------------------------- #


@dataclass(slots=True, frozen=True)
class Contrasts:
    """The four-arm comparison, or nothing if an arm is missing."""

    f00: float | None
    f10: float | None
    f01: float | None
    f11: float | None
    delta_pair: float | None
    delta_over_substrate: float | None
    delta_over_probiotic: float | None
    interaction: float | None
    missing_arms: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "arms": {
                "neither": self.f00, "probiotic_only": self.f10,
                "substrate_only": self.f01, "both": self.f11,
            },
            "delta_pair": self.delta_pair,
            "delta_over_substrate": self.delta_over_substrate,
            "delta_over_probiotic": self.delta_over_probiotic,
            "interaction": self.interaction,
            "missing_arms": list(self.missing_arms),
            "interpretation": (
                "An interaction contrast in the model, not proof of clinical synergy. "
                "No contrast is computed when a required arm is missing."
            ),
        }


def contrasts(
    f00: float | None, f10: float | None, f01: float | None, f11: float | None
) -> Contrasts:
    """Δpair, Δ over each component, and the interaction contrast.

    ``Δpair = F11 - F00``, ``Δover_substrate = F11 - F01``,
    ``Δover_probiotic = F11 - F10``, ``I = F11 - F10 - F01 + F00``.

    Every contrast needs all of its arms. A missing arm yields None for the
    contrasts that use it rather than a zero, because an unrun condition is
    not a condition that produced no change (AT104).
    """
    missing = tuple(
        name for name, value in
        (("neither", f00), ("probiotic_only", f10), ("substrate_only", f01), ("both", f11))
        if value is None
    )
    both_and = f11 is not None
    return Contrasts(
        f00=f00, f10=f10, f01=f01, f11=f11,
        delta_pair=(f11 - f00) if both_and and f00 is not None else None,
        delta_over_substrate=(f11 - f01) if both_and and f01 is not None else None,
        delta_over_probiotic=(f11 - f10) if both_and and f10 is not None else None,
        interaction=(
            f11 - f10 - f01 + f00
            if None not in (f00, f10, f01, f11) else None
        ),
        missing_arms=missing,
    )


# --------------------------------------------------------------------------- #
# 11.6 — sensitivity across predeclared assumptions
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class Sensitivity:
    """A scenario range across declared media and parameters.

    Deliberately min/median/max plus rank stability. Section 11.6: this is not
    a calibrated confidence interval and not a probability.
    """

    values: tuple[float, ...]
    minimum: float | None
    median: float | None
    maximum: float | None
    scenarios: tuple[str, ...]
    rank_changed: bool | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "min": self.minimum,
            "median": self.median,
            "max": self.maximum,
            "n_scenarios": len(self.values),
            "scenarios": list(self.scenarios),
            "rank_changed_across_scenarios": self.rank_changed,
            "interpretation": (
                "A range across predeclared assumptions, not a confidence interval "
                "and not a probability of benefit."
            ),
        }


def sensitivity(
    values: Mapping[str, float | None], *, rank_changed: bool | None = None
) -> Sensitivity:
    """min/median/max over the scenarios that actually solved."""
    solved = {k: float(v) for k, v in values.items() if v is not None and math.isfinite(float(v))}
    if not solved:
        return Sensitivity((), None, None, None, tuple(values), rank_changed)
    numbers = tuple(solved.values())
    return Sensitivity(
        values=numbers,
        minimum=min(numbers),
        median=statistics.median(numbers),
        maximum=max(numbers),
        scenarios=tuple(solved),
        rank_changed=rank_changed,
    )


# --------------------------------------------------------------------------- #
# 11.7 — the published-output replay
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ReplayResult:
    """The released-prediction replay, kept apart from a solved model."""

    execution_state: str
    sample_id: str
    raw: dict[str, float]
    normalized: dict[str, float]
    checks: tuple[dict[str, Any], ...]
    tolerance: float
    ranking_raw: tuple[str, ...]
    ranking_normalized: tuple[str, ...]

    @property
    def passed(self) -> bool:
        return all(c["within_tolerance"] for c in self.checks)

    def to_json(self) -> dict[str, Any]:
        return {
            "execution_state": self.execution_state,
            "sample_id": self.sample_id,
            "tolerance": self.tolerance,
            "passed": self.passed,
            "checks": list(self.checks),
            "ranking_raw_flux": list(self.ranking_raw),
            "ranking_flux_per_growth": list(self.ranking_normalized),
            "provenance": (
                "Values read from the released predictions of S10, not simulated here. "
                "This is a parser and ranker fixture; it is not an outcome measured in "
                "any person."
            ),
        }


def replay_published_fixture(
    config: SimulationConfig, *, workbook: Path | None = None
) -> ReplayResult:
    """Reproduce the published S5 butyrate predictions and their two orderings.

    Two objectives, asserted separately and under their correct names: raw
    butyrate flux, where the psyllium proxy ranks first, and butyrate flux per
    unit of community growth, where maltodextrin does. Swapping the units
    would swap the winner, which is why AT108 exists.
    """
    spec = config.replay
    path = workbook or (CONFIG_FILE.parents[1] / spec["file"])
    if not Path(path).is_file():
        raise SimulationError(
            f"the replay workbook is not installed at {path}; run `make extension-refs`"
        )
    excluded = set(spec.get("exclude_treatments") or ())
    tolerance = float(spec.get("tolerance") or 1e-8)
    sample_id = str(spec["sample_id"])

    with StrictWorkbook.open(path) as wb:
        raw_rows = wb.records(str(spec["sheet_raw"]))
        norm_rows = wb.records(str(spec["sheet_normalized"]))

    def index(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str], float]:
        out: dict[tuple[str, str], float] = {}
        for row in rows:
            if str(row.get("sample_id")) != sample_id:
                continue
            treatment = str(row.get("treatment") or "")
            probiotic = str(row.get("probiotic") or "-")
            if treatment in excluded:
                continue
            flux = row.get("flux")
            if flux is None:
                continue
            out[(treatment, probiotic)] = float(flux)
        return out

    raw = index(raw_rows)
    norm = index(norm_rows)
    if not raw:
        raise SimulationError(
            f"the replay workbook has no rows for sample {sample_id}; the released "
            "input and output sets differ in size (S10 audit item 4) and the "
            "intersection must be checked rather than assumed"
        )

    checks: list[dict[str, Any]] = []
    for expected in spec.get("raw_expected") or ():
        key = (str(expected["treatment"]), str(expected["probiotic"]))
        got = raw.get(key)
        want = float(expected["flux"])
        checks.append({
            "candidate": str(expected["candidate"]),
            "objective": "raw_butyrate_flux",
            "unit": "mmol/gDW/h",
            "expected": want,
            "actual": got,
            "within_tolerance": got is not None and abs(got - want) <= tolerance,
        })
        delta_want = expected.get("delta_from_no_addition")
        base = raw.get(("No Prebiotic", "-"))
        if delta_want is not None and got is not None and base is not None:
            checks.append({
                "candidate": str(expected["candidate"]),
                "objective": "raw_butyrate_flux_delta",
                "unit": "mmol/gDW/h",
                "expected": float(delta_want),
                "actual": got - base,
                "within_tolerance": abs((got - base) - float(delta_want)) <= tolerance,
            })
    for label, key in (
        ("baseline", ("No Prebiotic", "-")),
        ("maltodextrin_alone", ("Maltodextrin", "-")),
        ("psyllium_alone", ("Psyllium Husk", "-")),
    ):
        want = (spec.get("normalized_expected") or {}).get(label)
        if want is None:
            continue
        got = norm.get(key)
        checks.append({
            "candidate": label,
            "objective": "butyrate_flux_per_community_growth",
            "unit": "mmol/gDW/h per 1/h",
            "expected": float(want),
            "actual": got,
            "within_tolerance": got is not None and abs(got - float(want)) <= tolerance,
        })

    def order(table: Mapping[tuple[str, str], float]) -> tuple[str, ...]:
        alone = {t: v for (t, p), v in table.items() if p == "-" and t != "No Prebiotic"}
        return tuple(k for k, _ in sorted(alone.items(), key=lambda kv: -kv[1]))

    return ReplayResult(
        execution_state="published_replay",
        sample_id=sample_id,
        raw=dict(sorted({f"{t}|{p}": v for (t, p), v in raw.items()}.items())),
        normalized=dict(sorted({f"{t}|{p}": v for (t, p), v in norm.items()}.items())),
        checks=tuple(checks),
        tolerance=tolerance,
        ranking_raw=order(raw),
        ranking_normalized=order(norm),
    )


# --------------------------------------------------------------------------- #
# the solved model
# --------------------------------------------------------------------------- #


#: Member ranges are asked with the community held at this fraction of its
#: butyrate maximum. 0.99 rather than 1.0 because an optimum's face is thin:
#: at exactly the maximum, almost every member range collapses to a point that
#: says nothing about the alternatives one step away.
MEMBER_RANGE_PIN_FRACTION: Final = 0.99


@dataclass(frozen=True)
class MemberRange:
    """One member's contribution as a range: what it must do, what it can do.

    Values are abundance-weighted fluxes (mmol/gDW/h of community) and shares
    of the relevant community total. `must` above zero means no optimal state
    exists without this member's contribution; `can` at zero means the member
    carries the exchange but the community never uses it.
    """

    must_flux: float | None
    can_flux: float | None
    total: float
    kind: str

    @classmethod
    def of_production(cls, lo: float | None, hi: float | None, total: float) -> MemberRange:
        def clamp(v: float | None) -> float | None:
            return None if v is None else max(0.0, v)
        return cls(clamp(lo), clamp(hi), max(0.0, total), "production")

    @classmethod
    def of_uptake(cls, lo: float | None, hi: float | None, supply: float) -> MemberRange:
        # Uptake is negative flux, so the minimum flux is the most uptake.
        most = None if lo is None else max(0.0, -lo)
        least = None if hi is None else max(0.0, -hi)
        return cls(least, most, max(0.0, supply), "uptake")

    def share(self, value: float | None) -> float | None:
        if value is None or self.total <= 0:
            return None
        return min(1.0, value / self.total)

    @property
    def must_share(self) -> float | None:
        return self.share(self.must_flux)

    @property
    def can_share(self) -> float | None:
        return self.share(self.can_flux)

    @property
    def required(self) -> bool:
        return (self.must_flux or 0.0) > 1e-9

    @property
    def possible(self) -> bool:
        return (self.can_flux or 0.0) > 1e-9

    def to_json(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "must_flux": self.must_flux, "can_flux": self.can_flux,
            "must_share": self.must_share, "can_share": self.can_share,
            "required": self.required, "possible": self.possible,
            "unit": "mmol/gDW/h of community; share of community total",
        }


@dataclass(slots=True)
class ScenarioResult:
    """One solved arm: its objective values, or why there are none."""

    scenario_id: str
    arm: str
    execution_state: str
    solver_status: str | None
    raw_flux: float | None
    community_growth: float | None
    flux_per_growth: float | None
    target_exchange: str
    medium_id: str
    coverage: Coverage | None
    perturbation: Perturbation | None
    unavailable_reason: str | None = None
    elapsed_s: float | None = None
    other_exchanges: dict[str, float | None] = field(default_factory=dict)
    medium_state: dict[str, Any] = field(default_factory=dict)
    from_cache: bool = False
    #: The lower bound of the exchange at the same pinned growth. When it
    #: equals the upper bound the community's output is determined; when it is
    #: far below, the model cannot pin it down and says so.
    flux_lower: float | None = None
    #: Each member's share of the community's butyrate potential, as a range
    #: over every state in which the community makes at least 99% of its
    #: maximum: what it *must* contribute and what it *can*.
    taxon_butyrate_share: dict[str, MemberRange] = field(default_factory=dict)
    #: Each member's share of the arm's substrate supply, as the same kind of
    #: range. Who this fibre can feed - including organisms the report has
    #: flagged, which is the trade-off a reader needs to see.
    taxon_substrate_share: dict[str, MemberRange] = field(default_factory=dict)
    #: The substrate exchange this arm capped, if any.
    substrate_exchange: str | None = None
    #: Growth maxima recomputed by every worker process that touched this arm.
    #: Their agreement is the reproducibility check; it is recorded, not assumed.
    gmax_by_worker: dict[str, float | None] = field(default_factory=dict)
    n_lp_solved: int = 0
    n_lp_cached: int = 0
    #: Medium exchanges binding at the growth optimum, with marginal growth per
    #: unit of extra supply. Names what caps growth in this diet.
    growth_limits: dict[str, float] = field(default_factory=dict)
    #: Every member of the modelled community. Needed to say that an organism
    #: has *no* route to a substrate, which is a different and stronger
    #: statement than its not appearing in a table of those that do.
    community_taxa: tuple[str, ...] = ()

    @property
    def reproducible(self) -> bool:
        values = [v for v in self.gmax_by_worker.values() if v is not None]
        if len(values) < 2:
            return True
        low, high = min(values), max(values)
        return abs(high - low) <= REPRODUCIBILITY_TOLERANCE * max(abs(high), 1e-12)

    def to_json(self) -> dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "other_exchanges": dict(self.other_exchanges),
            "other_exchanges_note": (
                "Each is its own exchange with its own units. They are never summed "
                "into one short-chain-fatty-acid number."
            ),
            "medium": dict(self.medium_state),
            "served_from_cache": self.from_cache,
            "arm": self.arm,
            "execution_state": self.execution_state,
            "solver_status": self.solver_status,
            "target_exchange": self.target_exchange,
            "medium_id": self.medium_id,
            "raw_flux": self.raw_flux,
            "raw_flux_meaning": (
                "the most of this metabolite the community could release while growing at "
                "its optimum: a production potential, reproducible to every digit"
            ),
            "raw_flux_lower_bound": self.flux_lower,
            "raw_flux_unit": "mmol/gDW/h",
            "member_butyrate_share": {k: v.to_json() for k, v in self.taxon_butyrate_share.items()},
            "member_substrate_share": {k: v.to_json() for k, v in self.taxon_substrate_share.items()},
            "member_share_meaning": (
                "abundance-weighted contribution of each member as a fraction of the community "
                "total, over every state in which the community makes at least "
                f"{MEMBER_RANGE_PIN_FRACTION:.0%} of its butyrate maximum; `must` is the least "
                "it contributes in any such state and `can` the most"
            ),
            "substrate_exchange": self.substrate_exchange,
            "reproducibility": {
                "growth_maximum_by_worker": dict(self.gmax_by_worker),
                "agrees": self.reproducible,
                "tolerance": REPRODUCIBILITY_TOLERANCE,
                "method": (
                    "every value is the optimum of a linear programme solved by interior "
                    "point with crossover on one thread; the community growth maximum is "
                    "recomputed independently by each worker process and compared"
                ),
            },
            "linear_programmes": {"solved": self.n_lp_solved, "cached": self.n_lp_cached},
            "growth_limits": dict(self.growth_limits),
            "n_community_members": len(self.community_taxa),
            "growth_limits_meaning": (
                "medium exchanges whose declared supply bounds community growth in this "
                "diet, with the marginal growth per unit of extra supply; growth is a "
                "property of the diet definition here and is not compared between arms"
            ),
            "community_growth": self.community_growth,
            "community_growth_unit": "1/h",
            "flux_per_community_growth": self.flux_per_growth,
            "flux_per_community_growth_unit": "mmol/gDW/h per 1/h",
            "unavailable_reason": self.unavailable_reason,
            "elapsed_s": self.elapsed_s,
            "coverage": self.coverage.to_json() if self.coverage else None,
            "perturbation": self.perturbation.to_json() if self.perturbation else None,
            "limitations": [
                "A predicted exchange flux from a metabolic model. Not a stool "
                "concentration, not absorbed exposure, not a dose and not a symptom change.",
                "Raw flux and flux per unit of community growth are separate metrics with "
                "different units and different rankings.",
            ],
        }


def flux_per_growth(raw: float | None, growth: float | None) -> float | None:
    """Divide only by a finite positive growth rate.

    Section 11.6 keeps these as separate metrics. Dividing by a zero or
    negative growth would manufacture an enormous or negative productivity
    from a community that is not growing.
    """
    if raw is None or growth is None:
        return None
    if not math.isfinite(growth) or growth <= 0.0:
        return None
    return raw / growth


def solver_available() -> tuple[bool, str]:
    """Whether a usable QP solver is installed, and which one.

    A commercial solver is preferred when present, not for licensing reasons
    but for numerical ones: see `REPRODUCIBILITY_TOLERANCE`.
    """
    try:
        import optlang  # noqa: PLC0415
    except ImportError:
        return False, "optlang is not installed"
    for name in ("gurobi", "cplex", "osqp", "hybrid"):
        if optlang.available_solvers.get(name.upper()):
            return True, name
    return False, "no quadratic-programming solver found (gurobi, cplex or osqp)"


def solver_is_stable(name: str) -> bool:
    """Whether this solver's LP optimum is exact enough to quote a number from."""
    return name.lower() in DETERMINISTIC_SOLVERS


def check_reproducibility(
    community: Any, config: SimulationConfig, *, exchange: str | None = None,
) -> dict[str, Any]:
    """Solve the same community twice and report whether the answer holds.

    Run once per candidate, on the baseline arm, rather than for every arm:
    if a community's optimum is not reproducible there, it will not be
    reproducible in the other three either, and three more four-hundred-second
    solves would only confirm it more slowly.
    """
    exchange = exchange or BUTYRATE_EXCHANGE
    values: list[float] = []
    growths: list[float] = []
    for _ in range(2):
        solution = community.cooperative_tradeoff(
            fraction=config.protocol.tradeoff, fluxes=True, pfba=False,
        )
        growths.append(float(solution.growth_rate or 0.0))
        if exchange in solution.fluxes.columns and "medium" in solution.fluxes.index:
            raw = solution.fluxes.loc["medium", exchange]
            values.append(float(raw) if raw is not None and math.isfinite(float(raw)) else 0.0)
        else:
            values.append(0.0)
    spread = abs(values[0] - values[1])
    scale = max(abs(values[0]), abs(values[1]), 1e-9)
    growth_spread = abs(growths[0] - growths[1]) / max(abs(growths[0]), abs(growths[1]), 1e-12)
    relative = spread / scale
    return {
        "exchange": exchange,
        "solves": values,
        "growth_rates": growths,
        "relative_flux_spread": relative,
        "relative_growth_spread": growth_spread,
        "tolerance": REPRODUCIBILITY_TOLERANCE,
        "reproducible": relative <= REPRODUCIBILITY_TOLERANCE
        and growth_spread <= REPRODUCIBILITY_TOLERANCE,
    }


def micom_available() -> tuple[bool, str]:
    try:
        import micom  # noqa: PLC0415
    except ImportError:
        return False, "micom is not installed; run `make extension-refs`"
    return True, str(micom.__version__)


def model_library_installed(config: SimulationConfig) -> tuple[bool, Path | None]:
    path = CONFIG_FILE.parents[1] / "refs" / "micom" / config.protocol.model_asset
    return path.is_file(), (path if path.is_file() else None)


def readiness(config: SimulationConfig) -> dict[str, Any]:
    """Whether the simulation layer can run here, component by component.

    A feature flag controls whether simulations run; it does not decide
    whether they are implemented. This reports which of the three pieces -
    solver, MICOM, model library - is actually present.
    """
    has_micom, micom_note = micom_available()
    has_solver, solver_note = solver_available()
    has_models, model_path = model_library_installed(config)
    stable = has_solver and solver_is_stable(solver_note)
    return {
        "can_run": bool(has_micom and has_solver and has_models),
        "can_report_individual_fluxes": bool(stable),
        "micom": {"installed": has_micom, "detail": micom_note,
                  "pinned": config.protocol.micom_version},
        "solver": {
            "installed": has_solver, "detail": solver_note, "numerically_stable": stable,
            "lp_method": LP_METHOD,
            "note": (
                "Bounds are computed as linear programmes on a deterministic simplex, so "
                "repeated runs of the same community return the same number to every digit."
                if stable else
                "This solver's linear-programming optimum is not exact, so a quoted bound "
                "could move between runs."
            ),
        },
        "model_library": {
            "installed": has_models,
            "asset": config.protocol.model_asset,
            "path": str(model_path) if model_path else None,
            "source_id": config.protocol.model_source_id,
        },
        "replay_fixture_installed": (
            CONFIG_FILE.parents[1] / str(config.replay.get("file", ""))
        ).is_file(),
    }


# --------------------------------------------------------------------------- #
# metric emission
# --------------------------------------------------------------------------- #


def scenario_metrics(
    results: Sequence[ScenarioResult],
    contrast: Contrasts,
    *,
    candidate_id: str,
    substrate: Substrate | None,
    probiotic_label: str | None,
    config: SimulationConfig,
    sens: Sensitivity | None = None,
) -> tuple[ExtensionMetric, ...]:
    """Emit one candidate's scenario outputs as extension metrics."""
    fp = fingerprint(
        candidate_id, config.protocol.to_json(),
        [r.scenario_id for r in results], substrate.to_json() if substrate else None,
    )
    coverage = next((r.coverage for r in results if r.coverage), None)
    detail: dict[str, Any] = {
        "candidate_id": candidate_id,
        "protocol": config.protocol.to_json(),
        "substrate": substrate.to_json() if substrate else None,
        "probiotic": probiotic_label,
        "arms": [r.to_json() for r in results],
        "contrasts": contrast.to_json(),
        "coverage": coverage.to_json() if coverage else None,
        "sensitivity": sens.to_json() if sens else None,
    }
    limits = (
        "A model prediction, not a measurement: predicted flux is not a stool "
        "concentration, an absorbed exposure, a dose or a symptom change.",
        "An interaction contrast in the model is not proof of clinical synergy.",
        "Model coverage is of bacterial abundance only, and covering a species does "
        "not mean its model carries the pathway in question.",
    )
    out: list[ExtensionMetric] = []
    both = next((r for r in results if r.arm == "both"), None)
    state = "measured" if (both and both.raw_flux is not None) else "insufficient_coverage"
    if state != "measured" and both is not None and both.execution_state == "not_run":
        state = "not_assayed"

    out.append(ExtensionMetric(
        metric_id=f"ext083.simulation.{candidate_id}.butyrate_flux",
        label="Predicted butyrate exchange flux, both added",
        kind="model_prediction",
        state=state,
        value=both.raw_flux if (both and state == "measured") else None,
        unit="mmol/gDW/h" if state == "measured" else None,
        denominator="community_model_biomass" if state == "measured" else None,
        direction="descriptive",
        method_id=METHOD_SCENARIO,
        input_fingerprint=fp,
        group="simulation",
        feature_id="A16",
        analytical_confidence="provisional",
        evidence_maturity="genomic_prediction",
        assessable_fraction=(coverage.mapped_fraction if coverage else None),
        source_ids=("R08", "S10", "S11"),
        limitations=limits,
        extra=detail,
    ))
    if contrast.interaction is not None:
        out.append(ExtensionMetric(
            metric_id=f"ext083.simulation.{candidate_id}.interaction",
            label="Interaction contrast between the probiotic and the substrate",
            kind="model_prediction",
            state="measured",
            value=contrast.interaction,
            unit="mmol/gDW/h",
            denominator="community_model_biomass",
            direction="descriptive",
            method_id=METHOD_SCENARIO,
            input_fingerprint=fp,
            group="simulation",
            feature_id="A16",
            analytical_confidence="provisional",
            evidence_maturity="genomic_prediction",
            assessable_fraction=(coverage.mapped_fraction if coverage else None),
            source_ids=("R08", "S10", "S11"),
            limitations=limits,
            extra=detail,
        ))
    return tuple(out)




# --------------------------------------------------------------------------- #
# the runner: build a community, apply the medium, solve
# --------------------------------------------------------------------------- #

#: Where built community models and solved scenarios are cached. Keyed by the
#: content fingerprint of everything that can change an answer, so a cached
#: result can never belong to a different protocol, medium or community.
CACHE_DIR: Final = CONFIG_FILE.parents[1] / "refs" / "micom" / "cache"

#: The community-level butyrate exchange. MICOM gives each metabolite two
#: reactions: `EX_but_m` is the medium exchange, which is what leaves the whole
#: community, and `EX_but(e)` is the per-organism exchange. Community
#: production is the medium row of the `_m` form; reading the `(e)` form there
#: returns NaN, which must not be mistaken for zero production. A missing
#: reaction is an unavailable scenario, never a zero (AT103).
BUTYRATE_EXCHANGE: Final = "EX_but_m"
BUTYRATE_TAXON_EXCHANGE: Final = "EX_but(e)"

#: Metabolites reported beside butyrate. Each is a separate exchange with its
#: own units; they are never summed into one "SCFA" number.
TARGET_EXCHANGES: Final[Mapping[str, str]] = {
    "butyrate": "EX_but_m",
    "propionate": "EX_ppa_m",
    "acetate": "EX_ac_m",
    "lactate": "EX_lac_L_m",
}


#: The extracted model library, and the species index over it. Both are
#: cached on disk and in this process. `load_qiime_model_db` unpacks the whole
#: 2 GB archive into a fresh temporary directory on every call, so calling it
#: once per scenario arm meant extracting 8 GB to run one comparison - the
#: reason a four-arm run took half an hour and produced nothing.
EXTRACTED_DIR: Final = CONFIG_FILE.parents[1] / "refs" / "micom" / "extracted"
_MODEL_DB_CACHE: dict[str, Any] = {}


def extracted_model_db(config: SimulationConfig) -> tuple[Path, Any]:
    """Unpack the model library once and hand back its manifest.

    Keyed on the asset's own checksum, so a replaced library cannot be served
    from a stale extraction.
    """
    from micom.qiime_formats import load_qiime_model_db  # noqa: PLC0415

    installed, path = model_library_installed(config)
    if not installed or path is None:
        raise SimulationError(
            "the AGORA2 model library is not installed; run `make extension-refs`"
        )
    key = config.protocol.model_asset_md5
    if key in _MODEL_DB_CACHE:
        return _MODEL_DB_CACHE[key]

    import pandas as pd  # noqa: PLC0415

    root = EXTRACTED_DIR / key[:16]
    marker = root / ".extracted"
    if not marker.is_file():
        root.mkdir(parents=True, exist_ok=True)
        load_qiime_model_db(str(path), str(root))
        marker.write_text(key, encoding="utf-8")
    # The QIIME artifact unpacks to `<uuid>/data/`, which is the directory
    # MICOM wants as a `model_db`: it holds one JSON per species beside the
    # manifest that indexes them. Handing MICOM the `.qza` instead makes it
    # unpack all 2 GB again for every single call.
    data_dirs = sorted(root.glob("*/data"))
    if not data_dirs:
        raise SimulationError(f"the extracted model library at {root} has no data directory")
    folder = data_dirs[0]
    manifest = pd.read_csv(folder / "manifest.csv")
    _MODEL_DB_CACHE[key] = (folder, manifest)
    return folder, manifest


def model_species_index(config: SimulationConfig) -> dict[str, str]:
    """Species name -> model id, from the installed library's own manifest.

    Built from the manifest rather than by guessing an id format, and keyed on
    a normalised species name so a sample's `Blautia_wexlerae` and the
    library's `Blautia wexlerae` resolve to the same model. An organism with
    no model stays unmapped.
    """
    _folder, manifest = extracted_model_db(config)
    # Keyed to the manifest's *rank* column, because that is the column MICOM
    # joins a taxonomy on. The library's `id` is `Bacteroides_fragilis` while
    # its `species` is `Bacteroides fragilis`; handing MICOM the id under a
    # `species` heading matches nothing, builds a community of zero taxa, and
    # then fails deep inside the solver with an array-shape error instead of
    # saying so.
    rank = str(manifest["summary_rank"].iloc[0]) if "summary_rank" in manifest else "species"
    index: dict[str, str] = {}
    for _, row in manifest.iterrows():
        rank_value = str(row.get(rank) or row.get("id") or "")
        if not rank_value:
            continue
        for key in (str(row.get("id") or ""), rank_value):
            norm = _normalise_species(key)
            if norm:
                index.setdefault(norm, rank_value)
    return index


def community_from_inventory(
    organisms: Sequence[Mapping[str, Any]],
    index: Mapping[str, str],
) -> tuple[dict[str, float], dict[str, str]]:
    """The modelled community from the organism inventory, resolving renames.

    AGORA2 predates several GTDB renames, so the library holds *Ruminococcus
    gnavus* where the inventory says *Mediterraneibacter gnavus*. Those are
    one organism under two names, and the inventory already records the older
    one, so the alias is tried rather than the abundance being dropped. Only
    names the inventory itself supplies are used; no name is invented to find
    a model.

    Returns the abundance keyed by the *accepted* name, plus the alias that
    actually resolved, so the report can say which name matched.
    """
    abundance: dict[str, float] = {}
    resolved_via: dict[str, str] = {}
    for row in organisms:
        percent = row.get("percent")
        if percent is None or float(percent) <= 0:
            continue
        accepted = str(row.get("species") or "").replace("_", " ").strip()
        if not accepted:
            continue
        candidates = [accepted]
        for key in ("gtdb", "formerly"):
            alias = row.get(key)
            if alias:
                candidates.append(str(alias).replace("_", " ").strip())
        key_used = None
        for candidate in candidates:
            norm = _normalise_species(candidate)
            if norm in index:
                key_used = norm
                break
        name = key_used or _normalise_species(accepted)
        abundance[name] = abundance.get(name, 0.0) + float(percent)
        if key_used and key_used != _normalise_species(accepted):
            resolved_via[name] = accepted
    return abundance, resolved_via


def _normalise_species(name: str) -> str:
    """One spelling for a species, so two lanes' names can be compared."""
    cleaned = str(name).replace("_", " ").strip().lower()
    # Drop a GTDB split suffix and any strain tail after the binomial.
    parts = [p for p in cleaned.split() if p]
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return cleaned


@dataclass(slots=True)
class CommunityBuild:
    """A built community and what it cost to build."""

    taxonomy: Any
    coverage: Coverage
    model_folder: Path
    manifest: Any
    n_models: int
    model_db: str = ""

    @property
    def pickle_path(self) -> Path:
        return self.model_folder / "community.pickle"

    def members(self) -> dict[str, Any]:
        """The community's member ids and exchange reaction ids.

        Written once beside the pickle so that planning an arm's questions
        never needs the two-gigabyte community loaded in the planning process.
        """
        sidecar = self.model_folder / "members.json"
        if sidecar.is_file():
            import json as _json  # noqa: PLC0415

            return _json.loads(sidecar.read_text(encoding="utf-8"))
        return write_members_sidecar(self.load(), sidecar)

    def load(self) -> Any:
        """The built community, from the pickle if one was cached.

        Loading the pickle is what makes a second arm on the same community
        nearly free, and it is why the cache key is the community rather than
        the sample.
        """
        from micom import load_pickle  # noqa: PLC0415

        pickle_path = self.model_folder / "community.pickle"
        if pickle_path.is_file():
            community = load_pickle(str(pickle_path))
            # A pickled community restores its solver with whatever defaults
            # the interface supplies, and a first-order solver's answer moves
            # with its tolerances: the same arm returned 4.57 and then 3.83
            # because one community was freshly built and the other loaded.
            # Pinning the tolerances here makes the two paths agree.
            try:
                community.solver.configuration.tolerances.feasibility = 1e-9
                community.solver.configuration.tolerances.optimality = 1e-9
            except (AttributeError, ValueError):
                pass
            return community
        from micom import Community  # noqa: PLC0415

        return Community(
            self.taxonomy, model_db=self.model_db or None, progress=False, solver="osqp",
        )


def write_members_sidecar(community: Any, path: Path) -> dict[str, Any]:
    """Record who is in a community and which exchanges each member carries."""
    import json as _json  # noqa: PLC0415

    taxa = sorted(str(t) for t in community.taxa)
    member_exchanges: dict[str, list[str]] = {t: [] for t in taxa}
    medium_exchanges: list[str] = []
    for reaction in community.reactions:
        rid = str(reaction.id)
        if not rid.startswith("EX_"):
            continue
        if rid.endswith("_m"):
            medium_exchanges.append(rid)
            continue
        if "__" in rid:
            base, _, taxon = rid.rpartition("__")
            if taxon in member_exchanges:
                member_exchanges[taxon].append(base)
    abundance = {}
    try:
        for taxon, value in zip(community.taxonomy.index, community.taxonomy["abundance"], strict=True):
            abundance[str(taxon)] = float(value)
    except (AttributeError, KeyError, ValueError):
        pass
    payload = {
        "taxa": taxa,
        "abundance": abundance,
        "medium_exchanges": sorted(medium_exchanges),
        "member_exchanges": {t: sorted(v) for t, v in member_exchanges.items()},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    return payload


def build_community(
    abundances: Mapping[str, float],
    config: SimulationConfig,
    *,
    sample_id: str,
    threads: int = 4,  # noqa: ARG001 - built in-process; see the note below
) -> CommunityBuild:
    """Map the community onto the library and build it, once, in this process.

    Two things this deliberately does not do.

    It does not call `micom.workflows.build`. That wrapper builds samples in a
    multiprocessing pool, and on macOS - where the start method is `spawn` -
    each child re-imports `__main__`. Run from a script, the children re-run
    the script, which spawns more children: the first attempt at this sat at
    100% CPU for forty-five minutes and wrote nothing. One community per call
    needs no pool.

    And it computes coverage *before* MICOM renormalises anything, returning it
    alongside the build, so a caller always has the mapping loss in hand
    (AT102).
    """
    import pandas as pd  # noqa: PLC0415

    index = model_species_index(config)
    normalised = {_normalise_species(k): float(v) for k, v in abundances.items()}
    mapped, coverage = map_to_models(
        normalised, index.keys(), cutoff=config.protocol.cutoff,
    )
    if not mapped:
        raise SimulationError(
            f"{sample_id}: no organism in this community has a model in the library; "
            "the scenario is unavailable rather than empty"
        )
    extracted, library_manifest = extracted_model_db(config)
    rank = (
        str(library_manifest["summary_rank"].iloc[0])
        if "summary_rank" in library_manifest else "species"
    )
    # `id` is this community's own label for the member; the rank column is
    # what MICOM joins against the library. Both carry the library's rank
    # value so the join cannot miss.
    taxonomy = pd.DataFrame({
        "id": [index[name] for name in mapped],
        rank: [index[name] for name in mapped],
        "abundance": list(mapped.values()),
        "sample_id": sample_id,
    })
    # Keyed on the community itself, not the sample name: the four arms share
    # two distinct communities, and an identical community must never be built
    # twice.
    key = fingerprint(
        config.protocol.to_json(),
        sorted((k, round(v, 12)) for k, v in mapped.items()),
    )
    folder = CACHE_DIR / "models" / key.split(":")[1][:16]
    folder.mkdir(parents=True, exist_ok=True)
    pickle_path = folder / "community.pickle"
    if not pickle_path.is_file():
        from micom import Community  # noqa: PLC0415

        community = Community(
            taxonomy, model_db=str(extracted), progress=False, solver="osqp",
        )
        try:
            community.solver.configuration.tolerances.feasibility = 1e-9
            community.solver.configuration.tolerances.optimality = 1e-9
        except (AttributeError, ValueError):
            pass
        community.to_pickle(str(pickle_path))
        write_members_sidecar(community, folder / "members.json")
    return CommunityBuild(
        taxonomy=taxonomy, coverage=coverage, model_folder=folder,
        manifest=library_manifest, n_models=len(mapped), model_db=str(extracted),
    )


def _medium_frame(config: SimulationConfig, medium_id: str, caps: Mapping[str, float]):
    """The declared medium with the scenario's total exchange caps applied.

    A cap replaces the medium's bound rather than adding to it, which is what
    "total caps, not additions to the old cap" means in section 11.7. An
    exchange the medium does not carry is added, with its identifiers filled,
    so a substrate cannot be silently dropped.
    """
    import pandas as pd  # noqa: PLC0415

    entry = next((m for m in config.media if str(m["id"]) == medium_id), None)
    if entry is None:
        raise SimulationError(f"unknown medium {medium_id!r}")
    path = CONFIG_FILE.parents[1] / str(entry["file"])
    if not path.is_file():
        raise SimulationError(f"medium {medium_id!r} is not installed at {path}")
    frame = pd.read_csv(path)
    # The published CSVs carry duplicate `reaction`/`reaction.1` columns; the
    # specification requires asserting they agree rather than trusting one.
    if "reaction.1" in frame.columns and "reaction" in frame.columns:
        disagree = frame["reaction"].astype(str) != frame["reaction.1"].astype(str)
        if bool(disagree.any()):
            raise SimulationError(
                f"medium {medium_id!r}: duplicate reaction columns disagree on "
                f"{int(disagree.sum())} row(s); the medium cannot be trusted"
            )
    if "flux" not in frame.columns:
        raise SimulationError(f"medium {medium_id!r} has no flux column")
    frame = frame.copy()
    added: list[dict[str, Any]] = []
    for exchange, cap in caps.items():
        hit = frame["reaction"].astype(str) == exchange
        if bool(hit.any()):
            frame.loc[hit, "flux"] = float(cap)
        else:
            row = dict.fromkeys(frame.columns)
            row["reaction"] = exchange
            row["flux"] = float(cap)
            if "metabolite" in frame.columns:
                row["metabolite"] = exchange.replace("EX_", "").removesuffix("_m")
            if "global_id" in frame.columns:
                row["global_id"] = exchange
            if "reaction.1" in frame.columns:
                row["reaction.1"] = exchange
            added.append(row)
    if added:
        # One concat with typed columns: appending rows one at a time
        # emitted a FutureWarning per row into every run log.
        extra = pd.DataFrame(added, columns=frame.columns).astype(frame.dtypes.to_dict(), errors="ignore")
        frame = pd.concat([frame, extra], ignore_index=True)
    return frame


def _apply_medium(community: Any, config: SimulationConfig, medium_id: str,
                  caps: Mapping[str, float]) -> dict[str, Any]:
    """Set the declared medium with this scenario's total caps applied.

    Returns what was actually applied, because a medium entry the community
    has no exchange for is silently unusable and that has to be visible: a
    scenario running on 110 of 142 declared nutrients is not running on the
    declared diet.
    """
    frame = _medium_frame(config, medium_id, caps)
    available = {r.id for r in community.exchanges}
    wanted = {
        str(k): float(v)
        for k, v in zip(frame["reaction"].astype(str), frame["flux"], strict=True)
        if v is not None and float(v) > 0
    }
    applied = {k: v for k, v in wanted.items() if k in available}
    community.medium = applied
    missing_caps = sorted(k for k in caps if k not in available)
    return {
        "medium_id": medium_id,
        "n_declared": len(wanted),
        "n_applied": len(applied),
        "n_not_in_community": len(wanted) - len(applied),
        "caps_applied": {k: v for k, v in caps.items() if k in available},
        "caps_missing_from_community": missing_caps,
    }


def member_exchange_id(medium_exchange: str) -> str:
    """`EX_inulin_m` (the community's exchange) -> `EX_inulin(e)` (a member's)."""
    base = medium_exchange[:-2] if medium_exchange.endswith("_m") else medium_exchange
    return f"{base}(e)"


def plan_arm(
    build: CommunityBuild,
    config: SimulationConfig,
    *,
    arm_id: str,
    medium_id: str,
    caps: Mapping[str, float],
) -> LP.Arm:
    """Every linear programme one arm needs, planned without loading the model.

    Community growth; the target exchange's upper and lower bound at pinned
    growth; each other reported exchange's upper bound; each member's own
    butyrate capacity; and, when the arm caps a substrate, each member's
    maximum uptake of it. The last two are what turn a community number into
    "who makes it" and "who eats it".
    """
    members = build.members()
    frame = _medium_frame(config, medium_id, caps)
    medium = {
        str(k): float(v)
        for k, v in zip(frame["reaction"].astype(str), frame["flux"], strict=True)
        if v is not None and float(v) > 0
    }
    questions: list[LP.Question] = [
        LP.Question("growth", "community_objective", "max", pin_growth=False),
        LP.Question("exchange", BUTYRATE_EXCHANGE, "max"),
        LP.Question("exchange", BUTYRATE_EXCHANGE, "min"),
    ]
    for reaction in TARGET_EXCHANGES.values():
        if reaction != BUTYRATE_EXCHANGE:
            questions.append(LP.Question("exchange", reaction, "max"))
    # Member ranges are asked in the state the headline describes: growth
    # pinned, and butyrate held at 99% of its own maximum. Both directions,
    # so each member gets a "must" and a "can".
    member_exchanges: Mapping[str, Sequence[str]] = members.get("member_exchanges", {})
    pin = {"pin_exchange": BUTYRATE_EXCHANGE, "pin_fraction": MEMBER_RANGE_PIN_FRACTION}
    for taxon in sorted(member_exchanges):
        if BUTYRATE_TAXON_EXCHANGE in member_exchanges[taxon]:
            target = f"{BUTYRATE_TAXON_EXCHANGE}__{taxon}"
            questions.append(LP.Question("taxon_exchange", target, "min", **pin))
            questions.append(LP.Question("taxon_exchange", target, "max", **pin))
    for cap_exchange in sorted(caps):
        base = member_exchange_id(cap_exchange)
        for taxon in sorted(member_exchanges):
            if base in member_exchanges[taxon]:
                target = f"{base}__{taxon}"
                questions.append(LP.Question("taxon_exchange", target, "min", **pin))
                questions.append(LP.Question("taxon_exchange", target, "max", **pin))
    return LP.Arm(
        arm_id=arm_id, pickle_path=build.pickle_path, medium=medium,
        tradeoff=config.protocol.tradeoff, questions=tuple(questions),
    )


def assemble_arm(
    plan: LP.Arm,
    answers: LP.ArmAnswers,
    build: CommunityBuild,
    *,
    scenario_id: str,
    arm: str,
    medium_id: str,
    caps: Mapping[str, float],
    perturbation: Perturbation | None,
) -> ScenarioResult:
    """Turn an arm's answers into a result, or into the reason there is none."""
    common = {
        "scenario_id": scenario_id, "arm": arm, "target_exchange": BUTYRATE_EXCHANGE,
        "medium_id": medium_id, "coverage": build.coverage, "perturbation": perturbation,
    }
    members = build.members()
    available = set(members.get("medium_exchanges", []))
    medium_state = {
        "medium_id": medium_id,
        "n_declared": len(plan.medium),
        "n_applied": sum(1 for k in plan.medium if k in available),
        "n_not_in_community": sum(1 for k in plan.medium if k not in available),
        "caps_applied": {k: v for k, v in caps.items() if k in available},
        "caps_missing_from_community": sorted(k for k in caps if k not in available),
    }
    growth = answers.gmax
    gmax_by_worker = {str(k): v for k, v in answers.gmax_by_worker.items()}
    base_kwargs = {
        "solver_status": answers.gmax_status, "elapsed_s": answers.elapsed_s, **common,
    }

    def _fail(reason: str, state: str = "solved_model") -> ScenarioResult:
        result = ScenarioResult(
            execution_state=state, raw_flux=None, community_growth=growth,
            flux_per_growth=None, unavailable_reason=reason, **base_kwargs,
        )
        result.medium_state = medium_state
        result.gmax_by_worker = gmax_by_worker
        result.n_lp_solved, result.n_lp_cached = answers.n_solved, answers.n_cached
        result.growth_limits = dict(answers.growth_limits)
        result.community_taxa = tuple(members.get("taxa") or ())
        return result

    if growth is None or not math.isfinite(growth) or growth <= 0:
        return _fail(
            "this community cannot grow on the declared medium, so no production bound "
            "exists: unavailable, not zero",
        )
    high = answers.get(LP.Question("exchange", BUTYRATE_EXCHANGE, "max"))
    low = answers.get(LP.Question("exchange", BUTYRATE_EXCHANGE, "min"))
    if not high.optimal:
        return _fail(
            f"the {BUTYRATE_EXCHANGE} bound could not be solved ({high.status}); "
            "unavailable, not zero production",
        )
    if not answers.reproducible:
        return _fail(
            "two worker processes disagreed on this community's growth maximum "
            f"({gmax_by_worker}); the arm is withheld rather than quoted", "failed",
        )

    result = ScenarioResult(
        execution_state="solved_model", raw_flux=float(high.value or 0.0),
        community_growth=float(growth),
        flux_per_growth=flux_per_growth(float(high.value or 0.0), float(growth)),
        unavailable_reason=None, **base_kwargs,
    )
    result.flux_lower = low.value if low.optimal else None
    result.other_exchanges = {
        name: (a.value if (a := answers.get(LP.Question("exchange", rid, "max"))).optimal else None)
        for name, rid in TARGET_EXCHANGES.items() if rid != BUTYRATE_EXCHANGE
    }
    result.medium_state = medium_state
    result.gmax_by_worker = gmax_by_worker
    result.n_lp_solved, result.n_lp_cached = answers.n_solved, answers.n_cached
    result.growth_limits = dict(answers.growth_limits)
    result.community_taxa = tuple(members.get("taxa") or ())
    # Member fluxes are per gram of that member; multiplied by the member's
    # abundance they become its contribution to the community exchange, and
    # divided by the community total they become a share. An unweighted member
    # flux is meaningless here - it simply runs to the model's default bound.
    abundance: Mapping[str, float] = members.get("abundance", {})
    butyrate_total = float(high.value or 0.0)
    substrate_supply = {k: float(v) for k, v in caps.items()}
    ranges: dict[tuple[str, str], dict[str, float | None]] = {}
    for question in plan.questions:
        if question.kind != "taxon_exchange":
            continue
        base, _, taxon = question.target.rpartition("__")
        answer = answers.get(question)
        weight = float(abundance.get(taxon, 0.0))
        value = (float(answer.value) * weight) if answer.optimal and answer.value is not None else None
        ranges.setdefault((base, taxon), {})[question.direction] = value
    for (base, taxon), pair in sorted(ranges.items()):
        lo, hi = pair.get("min"), pair.get("max")
        if base == BUTYRATE_TAXON_EXCHANGE:
            # Production is a positive exchange flux.
            result.taxon_butyrate_share[taxon] = MemberRange.of_production(lo, hi, butyrate_total)
        else:
            medium_exchange = next((k for k in caps if member_exchange_id(k) == base), None)
            supply = substrate_supply.get(medium_exchange or "", 0.0)
            result.substrate_exchange = result.substrate_exchange or medium_exchange or base
            # Uptake is a negative exchange flux; min flux is maximum uptake.
            result.taxon_substrate_share[taxon] = MemberRange.of_uptake(lo, hi, supply)
    return result


def run_arm(
    build: CommunityBuild,
    config: SimulationConfig,
    *,
    scenario_id: str,
    arm: str,
    medium_id: str,
    caps: Mapping[str, float],
    perturbation: Perturbation | None,
    workers: int | None = None,
    progress: Any = None,
) -> ScenarioResult:
    """Plan, solve and assemble one arm on its own."""
    plan = plan_arm(build, config, arm_id=scenario_id, medium_id=medium_id, caps=caps)
    answers = LP.solve_arms([plan], cache_dir=CACHE_DIR, workers=workers, progress=progress)
    return assemble_arm(
        plan, answers[scenario_id], build, scenario_id=scenario_id, arm=arm,
        medium_id=medium_id, caps=caps, perturbation=perturbation,
    )


ARM_ORDER: Final[tuple[tuple[str, bool, bool], ...]] = (
    ("neither", False, False),
    ("probiotic_only", True, False),
    ("substrate_only", False, True),
    ("both", True, True),
)


def run_four_arms(
    baseline: Mapping[str, float],
    config: SimulationConfig,
    *,
    sample_id: str,
    candidate_id: str,
    substrate: Substrate | None,
    probiotic_species: Sequence[str] = (),
    medium_id: str = "european",
    workers: int | None = None,
    progress: Any = None,
    dry_run: bool = False,
) -> tuple[tuple[ScenarioResult, ...], Contrasts]:
    """The whole comparison: neither, probiotic, substrate, both.

    Non-intervention constraints are held fixed across the four arms - the
    same medium, the same tradeoff, the same community apart from the
    intervention itself - because that is the only way the contrasts mean
    anything. All four arms' linear programmes go to one worker pool.
    """
    caps = {substrate.exchange: substrate.cap} if substrate else {}
    total = sum(float(v) for v in baseline.values()) or 1.0
    resident = {k: float(v) / total for k, v in baseline.items()}
    perturbed = build_perturbation(
        baseline, protocol=config.protocol, added=tuple(probiotic_species),
    ) if probiotic_species else None

    builds: dict[bool, CommunityBuild] = {
        False: build_community(resident, config, sample_id=f"{sample_id}.resident"),
    }
    if perturbed is not None:
        builds[True] = build_community(perturbed.abundances, config, sample_id=f"{sample_id}.perturbed")

    plans: dict[str, LP.Arm] = {}
    for arm, with_probiotic, with_substrate in ARM_ORDER:
        if with_probiotic and perturbed is None:
            continue
        plans[arm] = plan_arm(
            builds[with_probiotic], config, arm_id=f"{candidate_id}.{arm}",
            medium_id=medium_id, caps=caps if with_substrate else {},
        )
    answers = LP.solve_arms(
        list(plans.values()), cache_dir=CACHE_DIR, workers=workers, progress=progress,
        dry_run=dry_run,
    )

    arms: list[ScenarioResult] = []
    values: dict[str, float | None] = {}
    for arm, with_probiotic, with_substrate in ARM_ORDER:
        if arm not in plans:
            arms.append(ScenarioResult(
                scenario_id=f"{candidate_id}.{arm}", arm=arm, execution_state="not_run",
                solver_status=None, raw_flux=None, community_growth=None,
                flux_per_growth=None, target_exchange=BUTYRATE_EXCHANGE,
                medium_id=medium_id, coverage=builds[False].coverage, perturbation=None,
                unavailable_reason="no probiotic species were supplied for this candidate",
            ))
            values[arm] = None
            continue
        result = assemble_arm(
            plans[arm], answers[plans[arm].arm_id], builds[with_probiotic],
            scenario_id=f"{candidate_id}.{arm}", arm=arm, medium_id=medium_id,
            caps=caps if with_substrate else {},
            perturbation=perturbed if with_probiotic else None,
        )
        arms.append(result)
        values[arm] = result.raw_flux

    return tuple(arms), contrasts(
        values["neither"], values["probiotic_only"], values["substrate_only"], values["both"],
    )


__all__ = [
    "BUTYRATE_EXCHANGE",
    "BUTYRATE_TAXON_EXCHANGE",
    "CACHE_DIR",
    "CONFIG_FILE",
    "EXTRACTED_DIR",
    "TARGET_EXCHANGES",
    "CommunityBuild",
    "build_community",
    "community_from_inventory",
    "extracted_model_db",
    "model_species_index",
    "run_arm",
    "run_four_arms",
    "plan_arm",
    "assemble_arm",
    "member_exchange_id",
    "write_members_sidecar",
    "ARM_ORDER",
    "DETERMINISTIC_SOLVERS",
    "EXECUTION_STATES",
    "OPTIMAL_STATUSES",
    "REPRODUCIBILITY_TOLERANCE",
    "Contrasts",
    "Coverage",
    "Perturbation",
    "Protocol",
    "ReplayResult",
    "ScenarioResult",
    "MemberRange",
    "MEMBER_RANGE_PIN_FRACTION",
    "Sensitivity",
    "SimulationConfig",
    "SimulationError",
    "Substrate",
    "build_perturbation",
    "check_reproducibility",
    "contrasts",
    "flux_per_growth",
    "LP_METHOD",
    "map_to_models",
    "micom_available",
    "model_library_installed",
    "readiness",
    "replay_published_fixture",
    "scenario_metrics",
    "sensitivity",
    "solver_available",
    "solver_is_stable",
]
