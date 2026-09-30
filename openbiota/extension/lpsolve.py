"""Deterministic linear programmes over a MICOM community, answered in parallel.

Why this exists
---------------
MICOM's `cooperative_tradeoff` answers "what does the community do?" with a
quadratic programme whose growth rates are unique but whose flux vector is
not: the same community returned a butyrate flux of 19.50 and then 12.50,
both labelled optimal. A report cannot quote a number that moves.

Every quantity this module produces is instead the optimal *value* of a
linear programme, and an LP's optimal value is unique even where its
solution vector is not. Community growth is maximised first; then, with
growth pinned to at least the protocol's tradeoff fraction of that maximum,
each exchange is maximised or minimised in turn. HiGHS's interior-point
method with crossover returns the same value to twelve digits as its dual
simplex, three times faster, on one thread; independent processes agree to
better than ten significant digits, and that agreement is checked, not
assumed.

Why it is parallel
------------------
One LP on a fifty-five-member community takes about thirteen seconds and
optlang rebuilds the whole matrix for every solve, so there is nothing to
warm-start. A candidate's four arms ask about two hundred such questions.
They are independent given each arm's growth maximum, so they are split
across worker processes; each worker loads the pickled community once,
recomputes the growth maximum (which doubles as a cross-process
reproducibility check) and answers its share. Values do not depend on which
worker answered them.

Why every answer is cached
--------------------------
Each answer is keyed on the community, the applied medium, the LP method
and the question. A regenerated report re-solves nothing it has already
solved, and adding a question to an arm solves only that question.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import json
import math
import multiprocessing
import os
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from .schema import fingerprint

#: Interior point with crossover: an exact vertex optimum, reproducible on one
#: thread, and about three times faster than dual simplex on these LPs.
LP_METHOD: Final = "interior point"

#: MICOM's own tolerances. A first-order or interior method's answer moves
#: with its tolerances; these are pinned so a freshly built and a pickled
#: community agree.
FEASIBILITY_TOLERANCE: Final = 1e-9

#: Answers from two independent processes must agree to this relative
#: tolerance or the run is flagged. With an exact LP method they agree to
#: every digit; this catches a future change of solver or model.
CROSS_PROCESS_TOLERANCE: Final = 1e-9

QUESTION_KINDS: Final[frozenset[str]] = frozenset({
    "growth",          # maximum community growth; the anchor for pinning
    "exchange",        # a community-level exchange, `EX_x_m`
    "taxon_exchange",  # one member's exchange, `EX_x(e)__Taxon`
    "taxon_growth",    # one member's growth rate
})
DIRECTIONS: Final[frozenset[str]] = frozenset({"max", "min"})


class LPError(RuntimeError):
    """A question that cannot be posed or a worker that cannot answer."""


@dataclass(frozen=True)
class Question:
    """One linear programme.

    `pin_growth` fixes community growth to the interval
    `[tradeoff * gmax, gmax]` before optimising the target, which is what makes
    the answer a *potential at near-optimal growth* rather than a value the
    community could only reach by not growing.
    """

    kind: str
    target: str
    direction: str = "max"
    pin_growth: bool = True
    #: Also hold a community exchange at or above `pin_fraction` of its own
    #: maximum (itself found at pinned growth). This is how a member's range is
    #: asked *in the state the headline describes*: "when the community makes
    #: the most butyrate it can, how much of it must - and can - come from you?"
    pin_exchange: str | None = None
    pin_fraction: float = 0.99

    def __post_init__(self) -> None:
        if self.kind not in QUESTION_KINDS:
            raise LPError(f"unknown question kind {self.kind!r}")
        if self.direction not in DIRECTIONS:
            raise LPError(f"unknown direction {self.direction!r}")
        if self.kind == "growth" and self.pin_growth:
            object.__setattr__(self, "pin_growth", False)
        if self.pin_exchange is not None and not self.pin_growth:
            raise LPError("an exchange can only be pinned at pinned growth")
        if not 0.0 < self.pin_fraction <= 1.0:
            raise LPError(f"pin_fraction must be in (0, 1], got {self.pin_fraction}")

    @property
    def key(self) -> str:
        parts: list[Any] = ["lp-question/1", self.kind, self.target, self.direction, self.pin_growth]
        if self.pin_exchange is not None:
            parts += ["pin", self.pin_exchange, self.pin_fraction]
        return fingerprint(*parts)


@dataclass(frozen=True)
class Answer:
    """The optimal value of one question, or why there is none."""

    value: float | None
    status: str
    #: Which worker answered; informational, never part of the value.
    worker: int = -1

    @property
    def optimal(self) -> bool:
        return self.value is not None and self.status == "optimal"


@dataclass(frozen=True)
class Job:
    """One worker's share of one arm.

    `medium` is the declared medium with the arm's caps already applied; the
    worker drops entries the community has no exchange for, exactly as the
    driver's own medium application does, so answers cannot depend on the
    worker.
    """

    job_id: str
    pickle_path: str
    medium: tuple[tuple[str, float], ...]
    tradeoff: float
    questions: tuple[Question, ...]
    lp_method: str = LP_METHOD


@dataclass(frozen=True)
class JobResult:
    job_id: str
    gmax: float | None
    gmax_status: str
    answers: dict[str, Answer]
    worker: int
    elapsed_s: float
    #: Medium exchanges whose bound was binding at the growth optimum, with the
    #: marginal growth per unit of extra supply. This is what caps growth.
    growth_limits: dict[str, float] = field(default_factory=dict)


def binding_medium_exchanges(community: Any, *, threshold: float = 1e-9) -> dict[str, float]:
    """Medium exchanges whose bound limits the objective just solved.

    Read from the reduced costs the solver leaves behind. The European diet
    caps a fifty-member community at 0.0539/h through a single nutrient -
    cobalt at 0.000167 mmol/gDW/h - and that is worth stating in a report
    rather than letting a reader compare growth rates that cannot differ.
    """
    try:
        vduals = community.solver.problem.vduals
    except AttributeError:
        return {}
    limits: dict[str, float] = {}
    for name, cost in vduals.items():
        if abs(cost) <= threshold or not name.startswith("EX_") or "__" in name:
            continue
        base = name.split("_reverse")[0]
        if not base.endswith("_m"):
            continue
        limits[base] = max(limits.get(base, 0.0), abs(float(cost)))
    return dict(sorted(limits.items(), key=lambda kv: -kv[1]))


#: One linear programme may not run longer than this. A community model
#: normally solves in seconds; one that has not solved in ten minutes is
#: cycling or degenerate, and its answer is "unanswered", not a two-hour wait
#: for every other job in the batch.
LP_TIME_LIMIT_S: Final = 600.0
#: The batch waits this long after its last completed job for the stragglers
#: before it gives them up as failed and moves on.
BATCH_STRAGGLER_S: Final = 1800.0


def _configure(community: Any, lp_method: str) -> None:
    community.solver.configuration.lp_method = lp_method
    try:
        community.solver.configuration.tolerances.feasibility = FEASIBILITY_TOLERANCE
        community.solver.configuration.tolerances.optimality = FEASIBILITY_TOLERANCE
    except (AttributeError, ValueError):
        pass
    with contextlib.suppress(AttributeError, ValueError, TypeError):
        community.solver.configuration.timeout = LP_TIME_LIMIT_S
    # One thread: HiGHS's interior point is deterministic on one thread and the
    # parallelism lives across processes instead.
    with contextlib.suppress(AttributeError, KeyError, TypeError):
        community.solver.problem.settings["threads"] = 1
    with contextlib.suppress(AttributeError, KeyError, TypeError):
        community.solver.problem.settings["time_limit"] = LP_TIME_LIMIT_S


def apply_medium(community: Any, medium: Mapping[str, float]) -> dict[str, float]:
    """Set the medium, keeping only entries the community can exchange.

    Returns what was applied. A declared nutrient the community has no
    exchange for is silently unusable, and that has to be visible upstream.
    """
    available = {r.id for r in community.exchanges}
    applied = {str(k): float(v) for k, v in medium.items() if k in available and float(v) > 0}
    community.medium = applied
    return applied


def _optimum(community: Any) -> tuple[float | None, str]:
    value = community.slim_optimize()
    status = str(getattr(community.solver, "status", "") or "unknown")
    if value is None or not math.isfinite(value):
        return None, status or "infeasible"
    return float(value), "optimal" if status in {"optimal", "unknown", ""} else status


def _objective_for(community: Any, question: Question) -> Any:
    if question.kind == "growth":
        return community.variables.community_objective
    if question.kind in {"exchange", "taxon_exchange"}:
        ids = {r.id for r in community.reactions}
        if question.target not in ids:
            raise LPError(f"no reaction {question.target!r} in this community")
        return community.reactions.get_by_id(question.target)
    if question.kind == "taxon_growth":
        try:
            return community.constraints["objective_" + question.target].expression
        except KeyError as exc:
            raise LPError(f"no member {question.target!r} in this community") from exc
    raise LPError(question.kind)


def answer_job(job: Job) -> JobResult:
    """Load one community, apply the medium and answer every question.

    Runs in a worker process. It must import cleanly under macOS's `spawn`
    start method, so nothing here depends on `__main__`.
    """
    import time  # noqa: PLC0415

    from micom import load_pickle  # noqa: PLC0415

    started = time.monotonic()
    worker = os.getpid()
    community = load_pickle(job.pickle_path)
    _configure(community, job.lp_method)
    apply_medium(community, dict(job.medium))

    community.objective = community.variables.community_objective
    community.objective_direction = "max"
    gmax, gmax_status = _optimum(community)
    growth_limits = binding_medium_exchanges(community) if gmax is not None else {}

    def _pin_growth() -> None:
        if gmax is not None:
            community.variables.community_objective.lb = job.tradeoff * gmax
            community.variables.community_objective.ub = gmax

    # Each pinned exchange's own maximum, found once per job at pinned growth.
    exchange_maxima: dict[str, tuple[float | None, str]] = {}

    def _exchange_maximum(exchange: str) -> tuple[float | None, str]:
        if exchange not in exchange_maxima:
            try:
                objective = _objective_for(community, Question("exchange", exchange, "max"))
            except LPError as exc:
                exchange_maxima[exchange] = (None, f"absent: {exc}")
            else:
                with community:
                    _pin_growth()
                    community.objective = objective
                    community.objective_direction = "max"
                    exchange_maxima[exchange] = _optimum(community)
        return exchange_maxima[exchange]

    answers: dict[str, Answer] = {}
    for question in job.questions:
        if question.kind == "growth":
            answers[question.key] = Answer(gmax, gmax_status, worker)
            continue
        if question.pin_growth and (gmax is None or gmax <= 0):
            answers[question.key] = Answer(None, "no_growth", worker)
            continue
        try:
            objective = _objective_for(community, question)
        except LPError as exc:
            answers[question.key] = Answer(None, f"absent: {exc}", worker)
            continue
        pinned_floor: float | None = None
        if question.pin_exchange is not None:
            emax, estatus = _exchange_maximum(question.pin_exchange)
            if emax is None:
                answers[question.key] = Answer(None, f"pin_unsolved: {estatus}", worker)
                continue
            # A floor only makes sense for a produced metabolite. If the
            # community cannot make it at all, the state "when it makes the
            # most it can" is just the pinned-growth state.
            pinned_floor = question.pin_fraction * emax if emax > 0 else None
        with community:
            if question.pin_growth:
                _pin_growth()
            if pinned_floor is not None and question.pin_exchange is not None:
                community.reactions.get_by_id(question.pin_exchange).lower_bound = pinned_floor
            community.objective = objective
            community.objective_direction = question.direction
            value, status = _optimum(community)
        answers[question.key] = Answer(value, status, worker)
    return JobResult(
        job_id=job.job_id, gmax=gmax, gmax_status=gmax_status, answers=answers,
        worker=worker, elapsed_s=round(time.monotonic() - started, 3),
        growth_limits=growth_limits,
    )


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Arm:
    """Everything one arm needs answered, and where its cache lives."""

    arm_id: str
    pickle_path: Path
    medium: Mapping[str, float]
    tradeoff: float
    questions: tuple[Question, ...]

    @property
    def cache_key(self) -> str:
        return fingerprint(
            "lp-arm/1", self.pickle_path.parent.name, LP_METHOD, self.tradeoff,
            sorted((k, round(float(v), 12)) for k, v in self.medium.items()),
        )


@dataclass
class ArmAnswers:
    arm_id: str
    gmax: float | None
    gmax_status: str
    answers: dict[str, Answer]
    #: Growth maxima reported by each worker that touched this arm. They must
    #: agree; this is the cross-process reproducibility check.
    gmax_by_worker: dict[int, float | None]
    n_solved: int
    n_cached: int
    elapsed_s: float
    growth_limits: dict[str, float] = field(default_factory=dict)

    def get(self, question: Question) -> Answer:
        return self.answers.get(question.key, Answer(None, "unanswered"))

    @property
    def reproducible(self) -> bool:
        values = [v for v in self.gmax_by_worker.values() if v is not None]
        if len(values) < 2:
            return True
        low, high = min(values), max(values)
        return abs(high - low) <= CROSS_PROCESS_TOLERANCE * max(abs(high), 1e-12)


def default_workers() -> int:
    """How many worker processes to use.

    Each holds a community in memory (about two gigabytes for fifty members),
    so this is capped well below the core count on large machines.
    """
    env = os.environ.get("OPENBIOTA_LP_WORKERS")
    if env and env.isdigit() and int(env) > 0:
        return int(env)
    cpus = os.cpu_count() or 4
    return max(1, min(16, cpus - 2))


def _cache_path(cache_dir: Path, arm: Arm) -> Path:
    return cache_dir / "lp" / f"{arm.cache_key.split(':')[1][:24]}.json"


def _read_cache(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _write_cache(path: Path, arm: Arm, gmax: float | None, gmax_status: str,
                 answers: Mapping[str, Answer], growth_limits: Mapping[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "arm_cache_key": arm.cache_key,
        "lp_method": LP_METHOD,
        "tradeoff": arm.tradeoff,
        "gmax": gmax,
        "gmax_status": gmax_status,
        "growth_limits": dict(growth_limits),
        "answers": {
            key: {"value": a.value, "status": a.status} for key, a in sorted(answers.items())
        },
    }
    path.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")


def _chunks(items: Sequence[Question], size: int) -> list[tuple[Question, ...]]:
    return [tuple(items[i:i + size]) for i in range(0, len(items), size)]


def solve_arms(
    arms: Sequence[Arm],
    *,
    cache_dir: Path,
    workers: int | None = None,
    questions_per_job: int = 8,
    progress: Any = None,
    dry_run: bool = False,
) -> dict[str, ArmAnswers]:
    """Answer every unanswered question of every arm, in parallel, then cache.

    Each arm's headline is solved twice - its growth maximum is recomputed by
    every worker that touches the arm - and the copies are compared. Where an
    arm has only one job, a second one-question job is added so that the check
    always runs.

    With `dry_run` nothing is solved: cached answers are returned and every
    unanswered question is left `unanswered`, so a caller can tell whether a
    cache would serve without starting a pool.
    """
    import time  # noqa: PLC0415

    started = time.monotonic()
    workers = workers or default_workers()
    out: dict[str, ArmAnswers] = {}
    jobs: list[Job] = []
    pending: dict[str, list[str]] = {}
    cached_answers: dict[str, dict[str, Answer]] = {}
    cached_gmax: dict[str, tuple[float | None, str]] = {}
    cached_limits: dict[str, dict[str, float]] = {}

    for arm in arms:
        cache = _read_cache(_cache_path(cache_dir, arm))
        have = {
            key: Answer(rec.get("value"), str(rec.get("status", "cached")))
            for key, rec in (cache.get("answers") or {}).items()
        }
        cached_answers[arm.arm_id] = have
        cached_gmax[arm.arm_id] = (cache.get("gmax"), str(cache.get("gmax_status", "")))
        cached_limits[arm.arm_id] = dict(cache.get("growth_limits") or {})
        todo = [q for q in arm.questions if q.key not in have]
        pending[arm.arm_id] = [q.key for q in todo]
        if not todo:
            continue
        chunks = _chunks(todo, questions_per_job)
        if len(chunks) == 1 and len(chunks[0]) > 1:
            # Split so two processes independently recompute the growth
            # maximum; that agreement is the reproducibility check.
            first = chunks[0]
            mid = max(1, len(first) // 2)
            chunks = [first[:mid], first[mid:]]
        for i, chunk in enumerate(chunks):
            jobs.append(Job(
                job_id=f"{arm.arm_id}#{i}", pickle_path=str(arm.pickle_path),
                medium=tuple(sorted((k, float(v)) for k, v in arm.medium.items())),
                tradeoff=arm.tradeoff, questions=chunk,
            ))

    results: dict[str, list[JobResult]] = {arm.arm_id: [] for arm in arms}
    if dry_run:
        jobs = []
    if jobs:
        if progress:
            progress(f"solving {sum(len(j.questions) for j in jobs)} linear programmes "
                     f"in {len(jobs)} jobs on {min(workers, len(jobs))} workers")
        # Spawned workers read the hash seed from the environment. Fixing it
        # removes one way two workers could order the same matrix differently.
        os.environ.setdefault("PYTHONHASHSEED", "0")
        context = multiprocessing.get_context("spawn")
        pool = concurrent.futures.ProcessPoolExecutor(max_workers=min(workers, len(jobs)), mp_context=context)
        futures = {pool.submit(answer_job, job): job for job in jobs}
        pending_futures = set(futures)
        done_count = 0
        try:
            while pending_futures:
                # Each job has a time limit per LP; the batch as a whole waits a
                # bounded time for stragglers after the others have finished,
                # then records them as failed rather than holding the report.
                finished, pending_futures = concurrent.futures.wait(
                    pending_futures, timeout=BATCH_STRAGGLER_S, return_when=concurrent.futures.FIRST_COMPLETED)
                if not finished:
                    for future in pending_futures:
                        job = futures[future]
                        results[job.job_id.split("#")[0]].append(JobResult(
                            job_id=job.job_id, gmax=None, gmax_status="failed: batch time limit",
                            answers={q.key: Answer(None, "failed: time limit") for q in job.questions},
                            worker=-1, elapsed_s=BATCH_STRAGGLER_S))
                        if progress:
                            progress(f"  job {job.job_id} gave up after {BATCH_STRAGGLER_S:.0f}s: recorded as unanswered")
                    for proc in list(getattr(pool, "_processes", {}).values()):
                        with contextlib.suppress(Exception):
                            proc.kill()
                    pending_futures = set()
                    break
                for future in finished:
                    done_count += 1
                    job = futures[future]
                    try:
                        result = future.result()
                    except Exception as exc:  # noqa: BLE001 - a failed job is a state
                        result = JobResult(
                            job_id=job.job_id, gmax=None, gmax_status=f"failed: {exc}"[:200],
                            answers={q.key: Answer(None, f"failed: {type(exc).__name__}")
                                     for q in job.questions},
                            worker=-1, elapsed_s=0.0,
                        )
                    results[job.job_id.split("#")[0]].append(result)
                    if progress:
                        progress(f"  job {done_count}/{len(jobs)} {job.job_id} "
                                 f"({len(job.questions)} LPs, {result.elapsed_s:.0f}s)")
        finally:
            pool.shutdown(wait=False, cancel_futures=True)

    for arm in arms:
        answers = dict(cached_answers[arm.arm_id])
        gmax_by_worker: dict[int, float | None] = {}
        gmax, gmax_status = cached_gmax[arm.arm_id]
        limits = dict(cached_limits[arm.arm_id])
        for result in results[arm.arm_id]:
            answers.update(result.answers)
            gmax_by_worker[result.worker] = result.gmax
            if result.gmax is not None:
                gmax, gmax_status = result.gmax, result.gmax_status
                limits = dict(result.growth_limits) or limits
        n_solved = 0 if dry_run else len(pending[arm.arm_id])
        n_cached = len(arm.questions) - len(pending[arm.arm_id])
        if n_solved:
            _write_cache(_cache_path(cache_dir, arm), arm, gmax, gmax_status, answers, limits)
        out[arm.arm_id] = ArmAnswers(
            arm_id=arm.arm_id, gmax=gmax, gmax_status=gmax_status or "cached",
            answers=answers, gmax_by_worker=gmax_by_worker,
            n_solved=n_solved, n_cached=n_cached,
            elapsed_s=round(time.monotonic() - started, 3), growth_limits=limits,
        )
    return out


def question_table(questions: Sequence[Question]) -> list[dict[str, Any]]:
    """A serialisable view of a question list, for audit output."""
    return [dict(asdict(q), key=q.key) for q in questions]


__all__ = [
    "CROSS_PROCESS_TOLERANCE",
    "DIRECTIONS",
    "FEASIBILITY_TOLERANCE",
    "LP_METHOD",
    "QUESTION_KINDS",
    "Answer",
    "Arm",
    "ArmAnswers",
    "Job",
    "JobResult",
    "LPError",
    "Question",
    "answer_job",
    "apply_medium",
    "binding_medium_exchanges",
    "default_workers",
    "question_table",
    "solve_arms",
]
