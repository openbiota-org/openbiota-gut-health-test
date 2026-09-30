"""Personalised synbiotic scenarios: every substrate against the consortium.

This is the layer between the solved arms (`simulation.run_four_arms`) and
the report. For one sample it:

* rebuilds the sample's community once (mapping losses recorded);
* runs the four arms for each declared substrate with the declared
  consortium, on the baseline diet and on each sensitivity diet;
* turns the arms into contrasts, and the member ranges into two tables a
  reader can act on - who would eat the substrate, and who would make the
  butyrate - with the sample's own flagged organisms marked;
* ranks candidates under both declared objectives, separately, and says when
  the two disagree;
* reports the reproducibility evidence and everything that is *not* a claim.

Terminology follows ISAPP: every combination here is a *candidate
combination* with `design_intent = synergistic` (it was built to have the
substrate feed the added organisms) and `definition_status = candidate`. A
positive interaction contrast in the model is evidence about the model, not
a demonstrated host benefit.
"""

from __future__ import annotations

import math
import pathlib
import re
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

import yaml

from openbiota.extension import simulation as SIM

#: Below this fraction of the community's butyrate potential, a member's
#: possible contribution is not worth a row. Rounding noise, not a role.
MEMBER_SHARE_FLOOR: Final = 0.005

#: An interaction contrast smaller than this fraction of the larger component
#: arm is called additive, not synergistic or antagonistic.
INTERACTION_FLOOR: Final = 0.01

#: Verdicts a candidate can receive, and what each means in plain words.
VERDICTS: Final[Mapping[str, str]] = {
    "pair_gains": (
        "the combination is predicted to release more butyrate than either component "
        "alone, and more than the sum of their separate gains"
    ),
    "pair_additive": (
        "the combination is predicted to release more than either component alone, "
        "about equal to the sum of their separate gains"
    ),
    "substrate_suffices": (
        "the substrate alone is predicted to do the work; adding the consortium adds "
        "little or nothing"
    ),
    "probiotic_suffices": (
        "the consortium alone is predicted to do the work; adding this substrate adds "
        "little or nothing"
    ),
    "no_gain": "neither the combination nor either component is predicted to raise butyrate",
    "pair_worse": (
        "the combination is predicted to release less than at least one component "
        "alone: the parts interfere in the model"
    ),
    "unavailable": "at least one arm could not be solved, so no contrast exists",
}

FLAG_CLASSES: Final[frozenset[str]] = frozenset({"opportunist", "conditional"})


class ProgressFn:  # pragma: no cover - typing helper
    def __call__(self, message: str) -> None: ...


@dataclass(frozen=True)
class OrganismContext:
    """What the rest of the report already says about one modelled member."""

    taxon: str
    display: str
    sample_species: str
    verdict_class: str | None
    flag_level: str | None
    flag_reason: str | None
    percentile: float | None
    percent: float | None

    @property
    def flagged(self) -> bool:
        return bool(self.flag_level) and (self.verdict_class in FLAG_CLASSES)

    @property
    def beneficial(self) -> bool:
        return self.verdict_class == "beneficial"

    def to_json(self) -> dict[str, Any]:
        return {
            "taxon": self.taxon, "display": self.display, "sample_species": self.sample_species,
            "class": self.verdict_class, "flag_level": self.flag_level,
            "flag_reason": self.flag_reason, "percentile": self.percentile,
            "percent": self.percent, "flagged": self.flagged, "beneficial": self.beneficial,
        }


def taxon_key(name: str) -> str:
    """One key for a member, whatever punctuation its library label carries.

    AGORA2 holds *Ruminococcus gnavus* as `[Ruminococcus] gnavus`, and MICOM
    turns the brackets into underscores when it names the member, so the
    library label and the community's own id for the same organism differ.
    Joined on the raw strings, this sample's second most abundant organism -
    flagged at the 98th percentile - lost its status silently and appeared in
    no table. Both sides go through here instead.
    """
    return re.sub(r"[^a-z0-9]+", "_", str(name).lower()).strip("_")


def context_for(
    contexts: Mapping[str, OrganismContext], taxon: str
) -> OrganismContext | None:
    """This member's entry, matched on the punctuation-free key."""
    return contexts.get(taxon_key(taxon))


def organism_contexts(
    results: Mapping[str, Any], index: Mapping[str, str]
) -> dict[str, OrganismContext]:
    """Join the organism verdicts onto model member ids, via the same aliases.

    The verdicts are what the report has already told the reader; a
    simulation that names an organism must use the same judgement of it.
    Keyed by `taxon_key`, not by the raw member id.
    """
    verdicts: dict[str, Mapping[str, Any]] = {}
    for row in (results.get("organism_verdicts") or {}).get("verdicts") or []:
        species = str(row.get("species") or "")
        if species:
            verdicts[SIM._normalise_species(species)] = row
    inventory = (results.get("organism_inventory") or {}).get("organisms") or []
    out: dict[str, OrganismContext] = {}
    for row in inventory:
        accepted = str(row.get("species") or "").replace("_", " ").strip()
        if not accepted:
            continue
        candidates = [accepted] + [
            str(row.get(k)).replace("_", " ").strip() for k in ("gtdb", "formerly") if row.get(k)
        ]
        model_label = next(
            (index[SIM._normalise_species(c)] for c in candidates
             if SIM._normalise_species(c) in index), None,
        )
        if model_label is None:
            continue
        taxon = taxon_key(model_label)
        verdict = verdicts.get(SIM._normalise_species(accepted), {})
        display = str(row.get("display") or verdict.get("display") or accepted)
        context = OrganismContext(
            taxon=taxon, display=display, sample_species=accepted,
            verdict_class=verdict.get("class"),
            flag_level=(verdict.get("flag") or None),
            flag_reason=verdict.get("flag_reason"),
            percentile=_float_or_none(verdict.get("percentile", row.get("percentile"))),
            percent=_float_or_none(row.get("percent")),
        )
        # Two sample rows can land on one model; keep the more abundant one.
        if taxon not in out or (context.percent or 0) > (out[taxon].percent or 0):
            out[taxon] = context
    return out


def _float_or_none(value: Any) -> float | None:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def display_name(taxon: str, contexts: Mapping[str, OrganismContext]) -> str:
    """What to call this member: the sample's own name for it where there is one."""
    context = context_for(contexts, taxon)
    if context:
        return context.display
    # A library label with no sample row: strip the bracket punctuation the
    # id carries so it reads as a species rather than as a fragment.
    return str(taxon).replace("_", " ").strip()


# --------------------------------------------------------------------------- #
# member tables
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MemberRow:
    taxon: str
    display: str
    must_share: float | None
    can_share: float | None
    required: bool
    possible: bool
    context: OrganismContext | None
    #: Present in this arm because the consortium added it.
    added: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "taxon": self.taxon, "display": self.display,
            "must_share": self.must_share, "can_share": self.can_share,
            "required": self.required, "possible": self.possible,
            "added_by_consortium": self.added,
            "flagged": bool(self.context and self.context.flagged),
            "beneficial": bool(self.context and self.context.beneficial),
            "class": self.context.verdict_class if self.context else None,
            "flag_level": self.context.flag_level if self.context else None,
            "percentile": self.context.percentile if self.context else None,
        }


def member_rows(
    ranges: Mapping[str, SIM.MemberRange],
    contexts: Mapping[str, OrganismContext],
    added: Sequence[str],
) -> list[MemberRow]:
    """Rows worth showing, required first, then by how much they can do."""
    added_set = {taxon_key(a) for a in added}
    rows: list[MemberRow] = []
    for taxon, rng in ranges.items():
        can = rng.can_share or 0.0
        if not rng.required and can < MEMBER_SHARE_FLOOR:
            continue
        rows.append(MemberRow(
            taxon=taxon, display=display_name(taxon, contexts),
            must_share=rng.must_share, can_share=rng.can_share,
            required=rng.required, possible=rng.possible,
            context=context_for(contexts, taxon), added=taxon_key(taxon) in added_set,
        ))
    rows.sort(key=lambda r: (not r.required, -(r.must_share or 0.0), -(r.can_share or 0.0), r.taxon))
    return rows


# --------------------------------------------------------------------------- #
# one candidate
# --------------------------------------------------------------------------- #


def verdict(c: SIM.Contrasts) -> str:
    """Name what the four arms say, without hiding a null or a negative."""
    if None in (c.f00, c.f10, c.f01, c.f11):
        return "unavailable"
    f00, f10, f01, f11 = float(c.f00), float(c.f10), float(c.f01), float(c.f11)  # type: ignore[arg-type]
    scale = max(abs(f00), abs(f10), abs(f01), abs(f11), 1e-9)
    tiny = INTERACTION_FLOOR * scale
    gain_pair = f11 - f00
    if f11 < max(f10, f01) - tiny:
        return "pair_worse"
    if gain_pair <= tiny:
        return "no_gain"
    over_sub, over_pro = f11 - f01, f11 - f10
    if over_sub <= tiny and over_pro > tiny:
        return "substrate_suffices"
    if over_pro <= tiny and over_sub > tiny:
        return "probiotic_suffices"
    if over_sub <= tiny and over_pro <= tiny:
        # Both components alone reach the pair: whichever you pick suffices.
        return "substrate_suffices" if f01 >= f10 else "probiotic_suffices"
    interaction = f11 - f10 - f01 + f00
    return "pair_gains" if interaction > tiny else "pair_additive"


@dataclass
class Candidate:
    candidate_id: str
    substrate: SIM.Substrate
    medium_id: str
    arms: tuple[SIM.ScenarioResult, ...]
    contrasts: SIM.Contrasts
    verdict: str
    eaters: dict[str, list[MemberRow]] = field(default_factory=dict)
    makers: dict[str, list[MemberRow]] = field(default_factory=dict)
    #: Flagged organisms in the modelled community that have no route to this
    #: substrate at all. Naming them is the strongest statement the model can
    #: make, and it is the one a reader worried about feeding an overgrown
    #: organism actually needs.
    flagged_without_route: list[OrganismContext] = field(default_factory=list)
    elapsed_s: float = 0.0

    def arm(self, name: str) -> SIM.ScenarioResult | None:
        return next((a for a in self.arms if a.arm == name), None)

    @property
    def solved(self) -> bool:
        return all(a.execution_state == "solved_model" and a.raw_flux is not None for a in self.arms)

    @property
    def reproducible(self) -> bool:
        return all(a.reproducible for a in self.arms)

    def flagged_eaters(self, arm: str) -> list[MemberRow]:
        return [r for r in self.eaters.get(arm, []) if r.context and r.context.flagged and r.possible]

    def per_growth(self, arm: str) -> float | None:
        result = self.arm(arm)
        return result.flux_per_growth if result else None

    def to_json(self) -> dict[str, Any]:
        both = self.arm("both")
        return {
            "candidate_id": self.candidate_id,
            "substrate": self.substrate.to_json(),
            "medium_id": self.medium_id,
            "design_intent": "synergistic",
            "definition_status": "candidate",
            "terminology_note": (
                "A candidate combination in ISAPP terms. Design intent says what the pairing "
                "was built to do; a model contrast is not a demonstrated host benefit."
            ),
            "arms": {a.arm: a.to_json() for a in self.arms},
            "contrasts": self.contrasts.to_json(),
            "verdict": self.verdict,
            "verdict_meaning": VERDICTS[self.verdict],
            "solved": self.solved,
            "reproducible": self.reproducible,
            "eaters": {arm: [r.to_json() for r in rows] for arm, rows in self.eaters.items()},
            "makers": {arm: [r.to_json() for r in rows] for arm, rows in self.makers.items()},
            "flagged_eaters_both": [r.to_json() for r in self.flagged_eaters("both")],
            "flagged_eaters_substrate_only": [
                r.to_json() for r in self.flagged_eaters("substrate_only")
            ],
            "flagged_without_route": [c.to_json() for c in self.flagged_without_route],
            "growth_limits": dict(both.growth_limits) if both else {},
            "elapsed_s": self.elapsed_s,
        }


def run_candidate(
    baseline: Mapping[str, float],
    config: SIM.SimulationConfig,
    substrate: SIM.Substrate,
    *,
    sample_id: str,
    medium_id: str,
    contexts: Mapping[str, OrganismContext],
    workers: int | None,
    progress: Callable[[str], None] | None,
    dry_run: bool = False,
) -> Candidate:
    started = time.monotonic()
    candidate_id = f"{substrate.substrate_id}+consortium@{medium_id}"
    arms, c = SIM.run_four_arms(
        baseline, config, sample_id=sample_id, candidate_id=candidate_id,
        substrate=substrate, probiotic_species=config.added_species,
        medium_id=medium_id, workers=workers, progress=progress, dry_run=dry_run,
    )
    candidate = Candidate(
        candidate_id=candidate_id, substrate=substrate, medium_id=medium_id,
        arms=arms, contrasts=c, verdict=verdict(c),
    )
    for arm in arms:
        if arm.taxon_substrate_share:
            candidate.eaters[arm.arm] = member_rows(
                arm.taxon_substrate_share, contexts, config.added_species,
            )
        if arm.taxon_butyrate_share:
            candidate.makers[arm.arm] = member_rows(
                arm.taxon_butyrate_share, contexts, config.added_species,
            )
    both = candidate.arm("both") or candidate.arm("substrate_only")
    if both is not None:
        reachable = {t for t, r in both.taxon_substrate_share.items() if r.possible}
        reachable_keys = {taxon_key(t) for t in reachable}
        candidate.flagged_without_route = sorted(
            (c for t in both.community_taxa
             if (c := context_for(contexts, t)) is not None and c.flagged
             and taxon_key(t) not in reachable_keys),
            key=lambda c: -(c.percentile or 0.0),
        )
    candidate.elapsed_s = round(time.monotonic() - started, 3)
    return candidate


# --------------------------------------------------------------------------- #
# ranking
# --------------------------------------------------------------------------- #


#: §24's seventeen strain-substrate evidence records.
SEEDS_FILE: Final = pathlib.Path(__file__).resolve().parents[2] / "extension" / "synbiotic_seeds.yaml"


def seeds() -> list[dict[str, Any]]:
    """The §24 records, as data.

    Most of them are negative. Four tested a strain, a substrate and the two
    together and found the combination no better than the strain alone, and
    eight lack a component-only arm and therefore cannot show synergy in
    either direction. Both facts are recorded rather than left as an absence,
    because a registry holding only the positives would describe a field that
    does not exist.
    """
    if not SEEDS_FILE.is_file():
        return []
    blob = yaml.safe_load(SEEDS_FILE.read_text(encoding="utf-8")) or {}
    return list(blob.get("records") or [])


def seed_summary() -> dict[str, Any]:
    """What the seed evidence base actually contains."""
    records = seeds()
    return {
        "n_records": len(records),
        "n_demonstrated_synergy": sum(1 for r in records if r.get("demonstrated_synergy") is True),
        "n_no_synergy_demonstrated": sum(
            1 for r in records if r.get("demonstrated_synergy") is False),
        "n_undetermined": sum(
            1 for r in records if r.get("demonstrated_synergy") is None),
        "n_missing_a_component_arm": sum(1 for r in records if r.get("arms_missing")),
        "n_with_incomplete_fields": sum(1 for r in records if r.get("incomplete_fields")),
        "settings": sorted({str(r.get("setting")) for r in records}),
        "note": (
            "Four of these are negative results: a strain, a substrate and the two "
            "together, with the combination no better than the strain alone. Eight "
            "lack a component-only arm and cannot show synergy in either direction. "
            "Neither fact is an absence of data."
        ),
    }


#: §11.4's two descriptive summaries, named as the specification names them.
COVERAGE_SCORES: Final[tuple[str, ...]] = ("goal_coverage", "pair_supported_goal_coverage")

COVERAGE_NOT: Final = (
    "Goal coverage is the weighted share of your eligible goals for which this "
    "candidate has a source-supported favourable edge. Pair coverage is the share the "
    "exact combination supports, rather than either component alone. Neither is a "
    "probability, a magnitude of benefit, or a health score, and neither enters the "
    "ranking."
)


class CoverageError(ValueError):
    """A coverage score that would be computed from invalid weights."""


def coverage(
    goals: Sequence[Mapping[str, Any]],
    *,
    n_components: int,
) -> dict[str, Any]:
    """§11.4's C and P, over the candidate's eligible goals.

    ``C = 100 * sum(w_j c_j) / sum(w_j)`` where ``c_j`` is 1 when the
    candidate has a source-supported favourable edge for that goal, and
    ``P`` the same over ``p_j``, which is 1 only when the exact combination
    supports the endpoint rather than either component alone.

    P is null for a single-component option, because there is no pair to
    support anything - not zero, which would read as a combination that was
    tried and failed. A conflicted goal keeps both its edges and is counted
    as covered *and* marked conflicted, because averaging the conflict away
    is how a contrary result disappears.
    """
    weights = [float(g.get("weight", 1.0)) for g in goals]
    for w in weights:
        if w < 0 or not math.isfinite(w):
            raise CoverageError(
                f"goal weight {w!r} is negative or not finite; a coverage score built "
                "from it would be meaningless"
            )
    total = sum(weights)
    if not goals or total <= 0:
        return {
            "goal_coverage": None,
            "pair_supported_goal_coverage": None,
            "n_goals": len(goals),
            "unavailable_reason": (
                "no eligible goals with positive weight, so there is nothing to cover"
            ),
            "goals": [],
            "what_these_are_not": COVERAGE_NOT,
        }

    covered = sum(w for w, g in zip(weights, goals, strict=True) if g.get("favourable"))
    paired = sum(w for w, g in zip(weights, goals, strict=True) if g.get("pair_supported"))
    single = n_components < 2
    return {
        "goal_coverage": round(100.0 * covered / total, 1),
        "pair_supported_goal_coverage": (
            None if single else round(100.0 * paired / total, 1)
        ),
        "pair_coverage_unavailable_reason": (
            "a single-component option has no pair to support an endpoint, so this is "
            "not applicable rather than zero" if single else None
        ),
        "n_goals": len(goals),
        "n_components": n_components,
        "goals": [
            {
                "goal": str(g.get("goal")),
                "weight": w,
                "favourable": bool(g.get("favourable")),
                "pair_supported": bool(g.get("pair_supported")),
                "conflicted": bool(g.get("conflicted")),
                # §11.4: for each covered goal show the evidence setting, so a
                # mouse result never wears the same percentage as a trial
                # without saying which it was.
                "evidence_setting": str(g.get("evidence_setting") or "unstated"),
            }
            for w, g in zip(weights, goals, strict=True)
        ],
        "n_conflicted": sum(1 for g in goals if g.get("conflicted")),
        "what_these_are_not": COVERAGE_NOT,
        "enters_the_ranking": False,
    }


#: Verdicts in which the exact combination, rather than either component
#: alone, supports the endpoint. This is §11.4's p_j, and the verdict
#: vocabulary already draws exactly that line, so it is read from there
#: rather than restated and left to drift.
PAIR_SUPPORTED_VERDICTS: Final[frozenset[str]] = frozenset({
    "pair_gains", "pair_additive",
})

#: Verdicts in which the candidate raises the endpoint at all, by any arm.
FAVOURABLE_VERDICTS: Final[frozenset[str]] = frozenset({
    "pair_gains", "pair_additive", "substrate_suffices", "probiotic_suffices",
})


def coverage_for(candidate: Candidate, *, endpoint: str = "butyrate production",
                 weight: float = 1.0) -> dict[str, Any]:
    """§11.4's C and P for one modelled candidate.

    The goal is the endpoint the model was run against, and the candidate's
    own verdict answers both questions the specification asks: whether it
    has a favourable edge at all, and whether the *pair* is what supports
    the endpoint rather than either component alone.

    The evidence setting is recorded as an in-silico model and never as a
    trial, because §11.4 is explicit that a percentage must not let one kind
    of evidence wear another's clothes.
    """
    verdict = str(candidate.verdict)
    return coverage(
        [{
            "goal": endpoint,
            "weight": weight,
            "favourable": verdict in FAVOURABLE_VERDICTS,
            "pair_supported": verdict in PAIR_SUPPORTED_VERDICTS,
            "conflicted": verdict == "pair_worse",
            "evidence_setting": "in-silico community model, not a trial",
        }],
        n_components=2,
    )


def rank(candidates: Sequence[Candidate], *, objective: str) -> list[dict[str, Any]]:
    """Rank the `both` arm under one objective. Unsolved candidates sit last, unranked."""
    def value(cand: Candidate) -> float | None:
        both = cand.arm("both")
        if both is None or both.raw_flux is None:
            return None
        return both.raw_flux if objective == "raw_butyrate_flux" else both.flux_per_growth

    scored = [(value(c), c) for c in candidates]
    ranked = sorted(
        [(v, c) for v, c in scored if v is not None], key=lambda vc: (-vc[0], vc[1].candidate_id),
    )
    out = [
        {"rank": i, "candidate_id": c.candidate_id, "substrate_id": c.substrate.substrate_id,
         "value": v, "verdict": c.verdict, "coverage": coverage_for(c)}
        for i, (v, c) in enumerate(ranked, start=1)
    ]
    out += [
        {"rank": None, "candidate_id": c.candidate_id, "substrate_id": c.substrate.substrate_id,
         "value": None, "verdict": c.verdict, "coverage": coverage_for(c)}
        for v, c in scored if v is None
    ]
    return out


def rank_stability(rankings: Mapping[str, Sequence[Mapping[str, Any]]]) -> dict[str, Any]:
    """Does each candidate keep its place across declared media? (AT111)

    A range and a rank, not a confidence interval.
    """
    positions: dict[str, dict[str, int | None]] = {}
    for medium_id, ranking in rankings.items():
        for row in ranking:
            positions.setdefault(str(row["substrate_id"]), {})[medium_id] = row["rank"]
    stable = {
        sid: len({r for r in by_medium.values() if r is not None}) <= 1
        for sid, by_medium in positions.items()
    }
    return {
        "positions": positions,
        "stable": stable,
        "all_stable": all(stable.values()) if stable else None,
        "meaning": (
            "rank of each candidate under each declared diet; a candidate is stable when "
            "its rank does not change. This is a range across predeclared assumptions, "
            "not an empirical confidence interval or a probability of response"
        ),
    }


# --------------------------------------------------------------------------- #
# the whole run
# --------------------------------------------------------------------------- #


def simulate(
    results: Mapping[str, Any],
    config: SIM.SimulationConfig | None = None,
    *,
    sample_id: str | None = None,
    workers: int | None = None,
    include_sensitivity: bool = True,
    progress: Callable[[str], None] | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Every declared substrate against the consortium, for this sample.

    `dry_run` answers only from the cache. Arms the cache cannot serve come
    back unsolved, and the view's `execution_state` says so; nothing is
    started.
    """
    started = time.monotonic()
    config = config or SIM.SimulationConfig.load()
    sample_id = str(sample_id or results.get("sample") or "sample")
    index = SIM.model_species_index(config)
    organisms = (results.get("organism_inventory") or {}).get("organisms") or []
    baseline, aliases = SIM.community_from_inventory(organisms, index)
    contexts = organism_contexts(results, index)

    media = [m for m in config.media if m.get("role") == "baseline"] or config.media[:1]
    if include_sensitivity:
        media += [m for m in config.media if m.get("role") == "sensitivity"]

    candidates: list[Candidate] = []
    for medium in media:
        for substrate in config.substrates:
            if progress:
                progress(f"candidate {substrate.substrate_id} on {medium['id']}")
            candidates.append(run_candidate(
                baseline, config, substrate, sample_id=sample_id, medium_id=str(medium["id"]),
                contexts=contexts, workers=workers, progress=progress, dry_run=dry_run,
            ))

    baseline_id = str(media[0]["id"])
    primary = [c for c in candidates if c.medium_id == baseline_id]
    rankings_raw = {
        m["id"]: rank([c for c in candidates if c.medium_id == m["id"]], objective="raw_butyrate_flux")
        for m in media
    }
    rankings_norm = {
        m["id"]: rank([c for c in candidates if c.medium_id == m["id"]],
                      objective="butyrate_flux_per_community_growth")
        for m in media
    }
    first_raw = next((r["substrate_id"] for r in rankings_raw[baseline_id] if r["rank"] == 1), None)
    first_norm = next((r["substrate_id"] for r in rankings_norm[baseline_id] if r["rank"] == 1), None)

    resident = primary[0].arm("neither") if primary else None
    n_lp_solved = sum(a.n_lp_solved for c in candidates for a in c.arms)
    n_lp_cached = sum(a.n_lp_cached for c in candidates for a in c.arms)
    n_unsolved_arms = sum(1 for c in candidates for a in c.arms if a.raw_flux is None)
    if all(c.solved for c in primary):
        state = "solved_model"
    elif dry_run:
        state = "not_run"
    else:
        state = "partial"
    return {
        "feature_id": "A16",
        "execution_state": state,
        "execution_state_meaning": {
            "solved_model": "every baseline-diet arm was solved for this sample",
            "not_run": (
                f"{n_unsolved_arms} arm(s) are not in the solution cache and the run was asked "
                "not to solve; rerun with --simulate auto (the default) to solve and cache them"
            ),
            "partial": f"{n_unsolved_arms} arm(s) could not be solved; each says why",
        }[state],
        "sample_id": sample_id,
        "protocol": config.protocol.to_json(),
        "consortium": list(config.added_species),
        "consortium_note": (
            "A species-level model perturbation following the published protocol: not an "
            "exact-strain product, not a dose."
        ),
        "community": {
            "n_members": resident.coverage.n_mapped_taxa if resident and resident.coverage else None,
            "coverage": resident.coverage.to_json() if resident and resident.coverage else None,
            "aliases_resolved": dict(aliases),
            "growth_limits": dict(resident.growth_limits) if resident else {},
            "growth_note": (
                "Community growth in this diet is capped by the nutrient(s) named in "
                "growth_limits and is identical across arms; growth is therefore a property of "
                "the diet definition here, and only metabolite potentials are compared."
            ),
        },
        "baseline_medium": baseline_id,
        "media": [{k: v for k, v in m.items() if k != "sha256"} for m in media],
        "candidates": [c.to_json() for c in candidates],
        "ranking": {
            "by_raw_butyrate_flux": rankings_raw,
            "by_butyrate_flux_per_community_growth": rankings_norm,
            "first_by_raw_flux": first_raw,
            "first_by_flux_per_growth": first_norm,
            "objectives_agree": (first_raw == first_norm) if first_raw and first_norm else None,
            "note": (
                "Two declared objectives, ranked separately; their units differ and their "
                "winners may differ. Neither ranking is a recommendation on its own."
            ),
            "stability_across_media": rank_stability(rankings_raw) if include_sensitivity else None,
        },
        "reproducibility": {
            "all_arms_agree": all(c.reproducible for c in candidates),
            "method": (
                "Every number is the optimum of a linear programme (interior point with "
                "crossover, one thread). Each arm's growth maximum is recomputed independently "
                "by every worker process that touches it; the copies agreed on every arm."
                if all(c.reproducible for c in candidates) else
                "At least one arm's independent re-solves disagreed; that arm is withheld."
            ),
            "linear_programmes_solved": n_lp_solved,
            "linear_programmes_from_cache": n_lp_cached,
        },
        "limitations": [
            "An experimental prediction from a metabolic model, not a measurement and not a "
            "guarantee of an individual response.",
            "A predicted exchange flux is not a stool concentration, an absorbed exposure or a "
            "dose. Model uptake caps are not oral doses.",
            "Members the model library lacks are absent from the community, and their share "
            "of the sample is recorded as unmapped rather than reassigned.",
            "Modelled feeding of a flagged organism is a hypothesis with a specific model edge "
            "behind it - the organism carries the transporter and the community state permits "
            "the uptake - not a rule that fibre feeds pathogens.",
            "The consortium is a species-level perturbation; strains, doses and survival "
            "through the stomach are outside the model.",
        ],
        "elapsed_s": round(time.monotonic() - started, 3),
    }


__all__ = [
    "FLAG_CLASSES",
    "INTERACTION_FLOOR",
    "MEMBER_SHARE_FLOOR",
    "VERDICTS",
    "Candidate",
    "MemberRow",
    "OrganismContext",
    "context_for",
    "display_name",
    "member_rows",
    "organism_contexts",
    "rank",
    "rank_stability",
    "run_candidate",
    "simulate",
    "taxon_key",
    "verdict",
]
