"""The deterministic LP engine and the member-range arithmetic on top of it.

None of these tests needs MICOM, a model library or a solver: the engine's
contract - question keys, chunking, the two-process reproducibility check,
the cache - is exercised with a stand-in worker, and the arithmetic that
turns member fluxes into shares is exercised directly.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.extension import lpsolve as LP
from openbiota.extension import simulation as SIM
from openbiota.extension import synbiotic as SYN

# --------------------------------------------------------------------------- #
# questions
# --------------------------------------------------------------------------- #


def test_question_key_is_stable_and_distinguishes_every_field():
    a = LP.Question("exchange", "EX_but_m", "max")
    b = LP.Question("exchange", "EX_but_m", "max")
    assert a.key == b.key
    assert a.key != LP.Question("exchange", "EX_but_m", "min").key
    assert a.key != LP.Question("exchange", "EX_ac_m", "max").key
    assert a.key != LP.Question("exchange", "EX_but_m", "max", pin_growth=False).key
    pinned = LP.Question("taxon_exchange", "EX_but(e)__X", "max", pin_exchange="EX_but_m")
    assert pinned.key != LP.Question("taxon_exchange", "EX_but(e)__X", "max").key
    assert pinned.key != LP.Question(
        "taxon_exchange", "EX_but(e)__X", "max", pin_exchange="EX_but_m", pin_fraction=0.9,
    ).key


def test_growth_question_never_pins_itself():
    q = LP.Question("growth", "community_objective", "max", pin_growth=True)
    assert q.pin_growth is False


@pytest.mark.parametrize("bad", [
    {"kind": "nonsense", "target": "x"},
    {"kind": "exchange", "target": "x", "direction": "sideways"},
    {"kind": "exchange", "target": "x", "pin_growth": False, "pin_exchange": "EX_but_m"},
    {"kind": "exchange", "target": "x", "pin_exchange": "EX_but_m", "pin_fraction": 0.0},
    {"kind": "exchange", "target": "x", "pin_exchange": "EX_but_m", "pin_fraction": 1.5},
])
def test_malformed_questions_are_refused(bad):
    with pytest.raises(LP.LPError):
        LP.Question(**bad)


# --------------------------------------------------------------------------- #
# the driver, with a stand-in worker
# --------------------------------------------------------------------------- #


def _fake_answer(job: LP.Job, *, gmax: float = 0.05, worker: int = 1) -> LP.JobResult:
    answers = {}
    for q in job.questions:
        if q.kind == "growth":
            answers[q.key] = LP.Answer(gmax, "optimal", worker)
        else:
            base = 10.0 if q.direction == "max" else 1.0
            answers[q.key] = LP.Answer(base + len(q.target) * 0.01, "optimal", worker)
    return LP.JobResult(job.job_id, gmax, "optimal", answers, worker, 0.1,
                        growth_limits={"EX_cobalt2_m": 323.0})


class _SerialPool:
    """Runs jobs in-process so the driver's bookkeeping can be tested."""

    def __init__(self, fn):
        self.fn = fn
        self.calls: list[LP.Job] = []

    def __call__(self, max_workers, mp_context):  # noqa: ARG002
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def shutdown(self, wait=True, cancel_futures=False):  # noqa: ARG002 - nothing to stop
        return None

    def submit(self, fn, job):  # noqa: ARG002 - the real worker is replaced
        import concurrent.futures

        self.calls.append(job)
        fut: concurrent.futures.Future = concurrent.futures.Future()
        try:
            fut.set_result(self.fn(job))
        except Exception as exc:  # noqa: BLE001 - a real pool captures this too
            fut.set_exception(exc)
        return fut


def _arm(tmp_path: Path, n_questions: int = 5, arm_id: str = "cand.neither") -> LP.Arm:
    questions = [LP.Question("growth", "community_objective", "max", pin_growth=False)]
    questions += [LP.Question("exchange", f"EX_m{i}_m", "max") for i in range(n_questions - 1)]
    return LP.Arm(
        arm_id=arm_id, pickle_path=tmp_path / "models" / "abc" / "community.pickle",
        medium={"EX_glc_D_m": 1.0, "EX_cobalt2_m": 0.000167}, tradeoff=0.99,
        questions=tuple(questions),
    )


def test_single_arm_is_split_across_two_jobs_so_growth_is_solved_twice(tmp_path, monkeypatch):
    counter = {"workers": 0}

    def worker_fn(job):
        counter["workers"] += 1
        return _fake_answer(job, worker=counter["workers"])

    pool = _SerialPool(worker_fn)
    monkeypatch.setattr(LP.concurrent.futures, "ProcessPoolExecutor", pool)
    out = LP.solve_arms([_arm(tmp_path)], cache_dir=tmp_path / "cache", workers=4)
    answers = out["cand.neither"]
    assert len(pool.calls) == 2, "one arm must still be answered by two independent processes"
    assert len(answers.gmax_by_worker) == 2
    assert answers.reproducible
    assert answers.n_solved == 5 and answers.n_cached == 0
    assert answers.growth_limits == {"EX_cobalt2_m": 323.0}


def test_disagreeing_workers_are_detected(tmp_path, monkeypatch):
    seen = {"n": 0}

    def worker_fn(job):
        seen["n"] += 1
        return _fake_answer(job, gmax=0.05 if seen["n"] == 1 else 0.0501, worker=seen["n"])

    monkeypatch.setattr(LP.concurrent.futures, "ProcessPoolExecutor", _SerialPool(worker_fn))
    out = LP.solve_arms([_arm(tmp_path)], cache_dir=tmp_path / "cache")
    assert out["cand.neither"].reproducible is False


def test_cache_is_written_and_a_second_call_solves_nothing(tmp_path, monkeypatch):
    pool = _SerialPool(_fake_answer)
    monkeypatch.setattr(LP.concurrent.futures, "ProcessPoolExecutor", pool)
    arm = _arm(tmp_path)
    first = LP.solve_arms([arm], cache_dir=tmp_path / "cache")
    assert first["cand.neither"].n_solved == 5
    files = list((tmp_path / "cache" / "lp").glob("*.json"))
    assert len(files) == 1
    payload = json.loads(files[0].read_text())
    assert payload["lp_method"] == LP.LP_METHOD
    assert payload["growth_limits"] == {"EX_cobalt2_m": 323.0}
    assert len(payload["answers"]) == 5

    calls_before = len(pool.calls)
    second = LP.solve_arms([arm], cache_dir=tmp_path / "cache")
    assert len(pool.calls) == calls_before, "a fully cached arm must not reach the pool"
    assert second["cand.neither"].n_solved == 0 and second["cand.neither"].n_cached == 5
    for q in arm.questions:
        assert first["cand.neither"].get(q).value == second["cand.neither"].get(q).value


def test_adding_one_question_solves_only_that_question(tmp_path, monkeypatch):
    pool = _SerialPool(_fake_answer)
    monkeypatch.setattr(LP.concurrent.futures, "ProcessPoolExecutor", pool)
    arm = _arm(tmp_path)
    LP.solve_arms([arm], cache_dir=tmp_path / "cache")
    extra = LP.Question("exchange", "EX_new_m", "min")
    bigger = LP.Arm(arm.arm_id, arm.pickle_path, arm.medium, arm.tradeoff, (*arm.questions, extra))
    out = LP.solve_arms([bigger], cache_dir=tmp_path / "cache")
    assert out["cand.neither"].n_solved == 1 and out["cand.neither"].n_cached == 5
    assert out["cand.neither"].get(extra).optimal


def test_arm_cache_key_depends_on_community_medium_and_tradeoff(tmp_path):
    a = _arm(tmp_path)
    b = LP.Arm(a.arm_id, a.pickle_path, {**a.medium, "EX_inulin_m": 6.14}, a.tradeoff, a.questions)
    c = LP.Arm(a.arm_id, a.pickle_path, a.medium, 0.5, a.questions)
    d = LP.Arm(a.arm_id, tmp_path / "models" / "zzz" / "community.pickle", a.medium, a.tradeoff, a.questions)
    assert len({a.cache_key, b.cache_key, c.cache_key, d.cache_key}) == 4
    # ...but not on the arm's label: the same LP under two names is one LP.
    e = LP.Arm("other.name", a.pickle_path, a.medium, a.tradeoff, a.questions)
    assert e.cache_key == a.cache_key


def test_failed_job_becomes_a_state_not_an_exception(tmp_path, monkeypatch):
    def worker_fn(job):  # noqa: ARG001
        raise RuntimeError("solver exploded")

    monkeypatch.setattr(LP.concurrent.futures, "ProcessPoolExecutor", _SerialPool(worker_fn))
    out = LP.solve_arms([_arm(tmp_path)], cache_dir=tmp_path / "cache")
    answers = out["cand.neither"]
    assert answers.gmax is None
    for q in _arm(tmp_path).questions:
        assert answers.get(q).optimal is False
        assert answers.get(q).status.startswith("failed")


def test_default_workers_respects_environment(monkeypatch):
    monkeypatch.setenv("OPENBIOTA_LP_WORKERS", "3")
    assert LP.default_workers() == 3
    monkeypatch.delenv("OPENBIOTA_LP_WORKERS")
    assert 1 <= LP.default_workers() <= 16


# --------------------------------------------------------------------------- #
# member ranges
# --------------------------------------------------------------------------- #


def test_production_range_clamps_and_shares():
    r = SIM.MemberRange.of_production(lo=12.0, hi=30.0, total=100.0)
    assert (r.must_share, r.can_share) == (0.12, 0.30)
    assert r.required and r.possible
    # A negative production bound is net uptake; it counts as no contribution.
    r2 = SIM.MemberRange.of_production(lo=-3.0, hi=0.0, total=100.0)
    assert (r2.must_flux, r2.can_flux) == (0.0, 0.0)
    assert not r2.required and not r2.possible


def test_uptake_range_flips_sign_and_orders_must_before_can():
    # Exchange flux -5 (most uptake) .. -1 (least uptake) of a 10 supply.
    r = SIM.MemberRange.of_uptake(lo=-5.0, hi=-1.0, supply=10.0)
    assert (r.must_flux, r.can_flux) == (1.0, 5.0)
    assert (r.must_share, r.can_share) == (0.1, 0.5)
    assert r.required
    # A member that can secrete but never must take up: must = 0.
    r2 = SIM.MemberRange.of_uptake(lo=-2.0, hi=4.0, supply=10.0)
    assert r2.must_flux == 0.0 and r2.can_flux == 2.0 and not r2.required


def test_share_is_capped_at_one_and_none_without_a_total():
    r = SIM.MemberRange.of_production(lo=0.0, hi=250.0, total=100.0)
    assert r.can_share == 1.0
    assert SIM.MemberRange.of_production(lo=1.0, hi=2.0, total=0.0).can_share is None


def test_unsolved_member_range_is_none_not_zero():
    r = SIM.MemberRange.of_production(lo=None, hi=None, total=100.0)
    assert r.must_share is None and r.can_share is None
    assert not r.required and not r.possible
    payload = r.to_json()
    assert payload["must_flux"] is None and payload["required"] is False


def test_member_exchange_id_maps_medium_to_member_form():
    assert SIM.member_exchange_id("EX_inulin_m") == "EX_inulin(e)"
    assert SIM.member_exchange_id("EX_but_m") == "EX_but(e)"


# --------------------------------------------------------------------------- #
# planning without a loaded model
# --------------------------------------------------------------------------- #


def _fake_build(tmp_path: Path) -> SIM.CommunityBuild:
    folder = tmp_path / "models" / "abc"
    folder.mkdir(parents=True)
    (folder / "members.json").write_text(json.dumps({
        "taxa": ["Bacteroides_fragilis", "Faecalibacterium_prausnitzii", "Ruminococcus_gnavus"],
        "abundance": {"Bacteroides_fragilis": 0.5, "Faecalibacterium_prausnitzii": 0.3,
                      "Ruminococcus_gnavus": 0.2},
        "medium_exchanges": ["EX_glc_D_m", "EX_inulin_m", "EX_but_m", "EX_ac_m", "EX_ppa_m",
                             "EX_lac_L_m"],
        "member_exchanges": {
            "Bacteroides_fragilis": ["EX_glc_D(e)", "EX_inulin(e)", "EX_ac(e)"],
            "Faecalibacterium_prausnitzii": ["EX_glc_D(e)", "EX_but(e)", "EX_ac(e)"],
            "Ruminococcus_gnavus": ["EX_glc_D(e)", "EX_ac(e)"],
        },
    }))
    coverage = SIM.Coverage(3, 3, 1.0, 0.0, 0.0, "bacterial", (), 0.001)
    return SIM.CommunityBuild(taxonomy=None, coverage=coverage, model_folder=folder,
                              manifest=None, n_models=3)


def test_plan_arm_asks_member_questions_only_of_carriers(tmp_path, monkeypatch):
    pd = pytest.importorskip("pandas")

    build = _fake_build(tmp_path)
    config = SIM.SimulationConfig.load()
    frame = pd.DataFrame({"reaction": ["EX_glc_D_m", "EX_cobalt2_m", "EX_inulin_m"],
                          "flux": [10.0, 0.000167, 6.14]})
    monkeypatch.setattr(SIM, "_medium_frame", lambda *_a, **_k: frame)
    plan = SIM.plan_arm(build, config, arm_id="c.both", medium_id="european",
                        caps={"EX_inulin_m": 6.14})
    kinds = [(q.kind, q.target, q.direction) for q in plan.questions]
    assert ("growth", "community_objective", "max") in kinds
    assert ("exchange", "EX_but_m", "max") in kinds and ("exchange", "EX_but_m", "min") in kinds
    # F. prausnitzii is the only butyrate carrier; B. fragilis the only inulin eater.
    member = [(t, d) for k, t, d in kinds if k == "taxon_exchange"]
    assert ("EX_but(e)__Faecalibacterium_prausnitzii", "min") in member
    assert ("EX_but(e)__Faecalibacterium_prausnitzii", "max") in member
    assert ("EX_inulin(e)__Bacteroides_fragilis", "min") in member
    assert not any("Ruminococcus_gnavus" in t for t, _ in member), (
        "an organism without the transporter is not asked about the substrate"
    )
    # Every member question is asked in the butyrate-optimal state.
    for q in plan.questions:
        if q.kind == "taxon_exchange":
            assert q.pin_exchange == SIM.BUTYRATE_EXCHANGE
            assert q.pin_fraction == SIM.MEMBER_RANGE_PIN_FRACTION
    assert plan.medium == {"EX_glc_D_m": 10.0, "EX_cobalt2_m": 0.000167, "EX_inulin_m": 6.14}


def test_assemble_arm_weights_member_fluxes_by_abundance(tmp_path, monkeypatch):
    pd = pytest.importorskip("pandas")

    build = _fake_build(tmp_path)
    config = SIM.SimulationConfig.load()
    frame = pd.DataFrame({"reaction": ["EX_glc_D_m", "EX_inulin_m"], "flux": [10.0, 6.14]})
    monkeypatch.setattr(SIM, "_medium_frame", lambda *_a, **_k: frame)
    caps = {"EX_inulin_m": 6.14}
    plan = SIM.plan_arm(build, config, arm_id="c.both", medium_id="european", caps=caps)

    def fake_answer(q: LP.Question) -> LP.Answer:
        if q.kind == "growth":
            return LP.Answer(0.05, "optimal")
        if q.target == "EX_but_m":
            return LP.Answer(40.0 if q.direction == "max" else 0.0, "optimal")
        if q.kind == "exchange":
            return LP.Answer(100.0, "optimal")
        if q.target.startswith("EX_but(e)__"):
            # per-gram flux 100 (must) .. 120 (can) for F. prausnitzii at abundance 0.3
            return LP.Answer(100.0 if q.direction == "min" else 120.0, "optimal")
        # inulin uptake by B. fragilis: flux -12 (most) .. -2 (least) at abundance 0.5
        return LP.Answer(-12.0 if q.direction == "min" else -2.0, "optimal")

    answers = LP.ArmAnswers(
        arm_id="c.both", gmax=0.05, gmax_status="optimal",
        answers={q.key: fake_answer(q) for q in plan.questions},
        gmax_by_worker={1: 0.05, 2: 0.05}, n_solved=len(plan.questions), n_cached=0,
        elapsed_s=1.0, growth_limits={"EX_cobalt2_m": 323.0},
    )
    result = SIM.assemble_arm(plan, answers, build, scenario_id="c.both", arm="both",
                              medium_id="european", caps=caps, perturbation=None)
    assert result.execution_state == "solved_model"
    assert result.raw_flux == 40.0 and result.flux_lower == 0.0
    assert result.community_growth == 0.05
    assert result.growth_limits == {"EX_cobalt2_m": 323.0}
    maker = result.taxon_butyrate_share["Faecalibacterium_prausnitzii"]
    assert maker.must_flux == pytest.approx(30.0) and maker.can_flux == pytest.approx(36.0)
    assert maker.must_share == pytest.approx(0.75) and maker.can_share == pytest.approx(0.90)
    eater = result.taxon_substrate_share["Bacteroides_fragilis"]
    assert eater.must_flux == pytest.approx(1.0) and eater.can_flux == pytest.approx(6.0)
    assert eater.can_share == pytest.approx(6.0 / 6.14)
    assert result.substrate_exchange == "EX_inulin_m"
    assert "Ruminococcus_gnavus" not in result.taxon_substrate_share
    payload = result.to_json()
    assert payload["reproducibility"]["agrees"] is True
    assert payload["member_butyrate_share"]["Faecalibacterium_prausnitzii"]["required"] is True


def test_assemble_arm_withholds_an_arm_whose_workers_disagree(tmp_path, monkeypatch):
    pd = pytest.importorskip("pandas")

    build = _fake_build(tmp_path)
    config = SIM.SimulationConfig.load()
    monkeypatch.setattr(SIM, "_medium_frame", lambda *_a, **_k: pd.DataFrame(
        {"reaction": ["EX_glc_D_m"], "flux": [10.0]}))
    plan = SIM.plan_arm(build, config, arm_id="c.neither", medium_id="european", caps={})
    answers = LP.ArmAnswers(
        arm_id="c.neither", gmax=0.05, gmax_status="optimal",
        answers={q.key: LP.Answer(1.0, "optimal") for q in plan.questions},
        gmax_by_worker={1: 0.05, 2: 0.06}, n_solved=1, n_cached=0, elapsed_s=1.0,
    )
    result = SIM.assemble_arm(plan, answers, build, scenario_id="c.neither", arm="neither",
                              medium_id="european", caps={}, perturbation=None)
    assert result.execution_state == "failed"
    assert result.raw_flux is None
    assert "disagreed" in (result.unavailable_reason or "")


def test_assemble_arm_reports_no_growth_as_unavailable_not_zero(tmp_path, monkeypatch):
    pd = pytest.importorskip("pandas")

    build = _fake_build(tmp_path)
    config = SIM.SimulationConfig.load()
    monkeypatch.setattr(SIM, "_medium_frame", lambda *_a, **_k: pd.DataFrame(
        {"reaction": ["EX_glc_D_m"], "flux": [10.0]}))
    plan = SIM.plan_arm(build, config, arm_id="c.neither", medium_id="european", caps={})
    answers = LP.ArmAnswers(
        arm_id="c.neither", gmax=0.0, gmax_status="optimal",
        answers={q.key: LP.Answer(None, "no_growth") for q in plan.questions},
        gmax_by_worker={1: 0.0}, n_solved=1, n_cached=0, elapsed_s=1.0,
    )
    result = SIM.assemble_arm(plan, answers, build, scenario_id="c.neither", arm="neither",
                              medium_id="european", caps={}, perturbation=None)
    assert result.raw_flux is None
    assert "not zero" in (result.unavailable_reason or "")


# --------------------------------------------------------------------------- #
# verdicts
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("arms", "expected"), [
    ((1.0, 3.0, 4.0, 8.0), "pair_gains"),          # the AT104 fixture: I = 2 > 0
    ((1.0, 3.0, 4.0, 6.0), "pair_additive"),       # I = 0
    ((54.8, 169.7, 54.8, 202.6), "pair_gains"),    # this sample's inulin result
    ((10.0, 12.0, 30.0, 30.2), "substrate_suffices"),
    ((10.0, 30.0, 12.0, 30.1), "probiotic_suffices"),
    ((10.0, 10.0, 10.0, 10.05), "no_gain"),
    ((10.0, 30.0, 25.0, 20.0), "pair_worse"),
    ((10.0, None, 25.0, 20.0), "unavailable"),
])
def test_verdict_names_what_the_arms_say(arms, expected):
    c = SIM.contrasts(*arms)
    assert SYN.verdict(c) == expected
    assert expected in SYN.VERDICTS


def test_verdict_keeps_a_negative_interaction_visible():
    # Component gains that do not add up: the pair is below the sum of parts
    # but still above each part. That is additive-or-less, never "gains".
    c = SIM.contrasts(10.0, 20.0, 20.0, 25.0)
    assert SYN.verdict(c) == "pair_additive"
    assert c.interaction == pytest.approx(-5.0)


# --------------------------------------------------------------------------- #
# organism context and member rows
# --------------------------------------------------------------------------- #


def _results_with_gnavus() -> dict:
    return {
        "organism_inventory": {"organisms": [
            {"species": "Mediterraneibacter_gnavus", "formerly": "Ruminococcus gnavus",
             "percent": 4.2, "percentile": 98.4},
            {"species": "Phocaeicola_vulgatus", "percent": 12.0, "percentile": 93.4},
            {"species": "Faecalibacterium_prausnitzii", "percent": 0.5, "percentile": 3.0},
            {"species": "GGB3000_SGB3991", "percent": 1.0},
        ]},
        "organism_verdicts": {"verdicts": [
            {"species": "Mediterraneibacter_gnavus", "class": "opportunist", "flag": "high",
             "flag_reason": "above the 95th percentile", "percentile": 98.4,
             "display": "Mediterraneibacter gnavus"},
            {"species": "Phocaeicola_vulgatus", "class": "conditional", "flag": "high",
             "percentile": 93.4, "display": "Phocaeicola vulgatus"},
            {"species": "Faecalibacterium_prausnitzii", "class": "beneficial", "flag": "low",
             "percentile": 3.0, "display": "Faecalibacterium prausnitzii"},
        ]},
    }


def test_organism_context_joins_verdicts_through_the_same_alias():
    index = {"ruminococcus gnavus": "Ruminococcus gnavus",
             "phocaeicola vulgatus": "Phocaeicola vulgatus",
             "faecalibacterium prausnitzii": "Faecalibacterium prausnitzii"}
    contexts = SYN.organism_contexts(_results_with_gnavus(), index)
    gnavus = SYN.context_for(contexts, "Ruminococcus_gnavus")
    assert gnavus is not None
    assert gnavus.sample_species == "Mediterraneibacter gnavus"
    assert gnavus.display == "Mediterraneibacter gnavus"
    assert gnavus.flagged and gnavus.verdict_class == "opportunist"
    assert gnavus.percentile == pytest.approx(98.4)
    vulgatus = SYN.context_for(contexts, "Phocaeicola_vulgatus")
    assert vulgatus is not None and vulgatus.flagged
    fp = SYN.context_for(contexts, "Faecalibacterium_prausnitzii")
    assert fp is not None
    assert fp.beneficial and not fp.flagged, "a low flag on a beneficial organism is not a warning"
    assert "GGB3000_SGB3991" not in {c.sample_species for c in contexts.values()}


def test_member_rows_put_required_first_and_drop_noise():
    contexts = SYN.organism_contexts(_results_with_gnavus(), {
        "ruminococcus gnavus": "Ruminococcus gnavus",
        "phocaeicola vulgatus": "Phocaeicola vulgatus",
    })
    ranges = {
        "Phocaeicola_vulgatus": SIM.MemberRange.of_uptake(-6.14, 0.0, 6.14),   # can 100%
        "Ruminococcus_gnavus": SIM.MemberRange.of_uptake(-0.01, 0.0, 6.14),    # can 0.16%
        "Anaerobutyricum_hallii": SIM.MemberRange.of_uptake(-2.0, -1.0, 6.14),  # must 16%
    }
    rows = SYN.member_rows(ranges, contexts, added=["Anaerobutyricum hallii"])
    assert [r.taxon for r in rows] == ["Anaerobutyricum_hallii", "Phocaeicola_vulgatus"]
    assert rows[0].required and rows[0].added
    assert rows[1].context is not None and rows[1].context.flagged
    payload = rows[1].to_json()
    assert payload["flagged"] is True and payload["class"] == "conditional"


def test_rank_orders_by_the_named_objective_and_leaves_unsolved_unranked():
    sub = SIM.Substrate("inulin", "Inulin", "EX_inulin_m", 6.14)
    sub2 = SIM.Substrate("pectin", "Pectin", "EX_pect_m", 0.4)
    sub3 = SIM.Substrate("starch", "Starch", "EX_strch1_m", 16.65)

    def cand(substrate, both_flux, growth):
        arms = tuple(
            SIM.ScenarioResult(
                scenario_id=f"{substrate.substrate_id}.{arm}", arm=arm,
                execution_state="solved_model" if both_flux is not None else "failed",
                solver_status="optimal", raw_flux=both_flux, community_growth=growth,
                flux_per_growth=SIM.flux_per_growth(both_flux, growth) if both_flux else None,
                unavailable_reason=None, target_exchange=SIM.BUTYRATE_EXCHANGE,
                medium_id="european", coverage=None, perturbation=None,
            ) for arm, *_ in SIM.ARM_ORDER
        )
        c = SIM.contrasts(both_flux, both_flux, both_flux, both_flux)
        return SYN.Candidate(substrate.substrate_id, substrate, "european", arms, c, SYN.verdict(c))

    cands = [cand(sub, 200.0, 0.05), cand(sub2, 300.0, 0.10), cand(sub3, None, None)]
    raw = SYN.rank(cands, objective="raw_butyrate_flux")
    assert [r["substrate_id"] for r in raw] == ["pectin", "inulin", "starch"]
    assert raw[2]["rank"] is None and raw[2]["value"] is None
    norm = SYN.rank(cands, objective="butyrate_flux_per_community_growth")
    # 200/0.05 = 4000 beats 300/0.10 = 3000: the two objectives pick different winners.
    assert [r["substrate_id"] for r in norm][:2] == ["inulin", "pectin"]


def test_rank_stability_is_a_rank_range_not_an_interval():
    rankings = {
        "european": [{"substrate_id": "inulin", "rank": 1}, {"substrate_id": "pectin", "rank": 2}],
        "high_fibre": [{"substrate_id": "inulin", "rank": 1}, {"substrate_id": "pectin", "rank": 3}],
    }
    out = SYN.rank_stability(rankings)
    assert out["stable"] == {"inulin": True, "pectin": False}
    assert out["all_stable"] is False
    assert "not an empirical confidence interval" in out["meaning"]


def test_bracketed_library_labels_still_join_to_their_verdict():
    """AGORA2's `[Ruminococcus] gnavus` and MICOM's `_Ruminococcus_gnavus`.

    These are one organism. Joined on the raw strings they are two, and the
    sample's second most abundant organism - flagged at the 98th percentile -
    lost its status and appeared in no table until v0.8.3.
    """
    index = {"ruminococcus gnavus": "[Ruminococcus] gnavus"}
    contexts = SYN.organism_contexts(_results_with_gnavus(), index)
    for member_id in ("_Ruminococcus_gnavus", "[Ruminococcus]_gnavus",
                      "[Ruminococcus] gnavus", "Ruminococcus_gnavus"):
        found = SYN.context_for(contexts, member_id)
        assert found is not None, f"{member_id} must find its verdict"
        assert found.flagged and found.verdict_class == "opportunist"
        assert SYN.display_name(member_id, contexts) == "Mediterraneibacter gnavus"


def test_taxon_key_ignores_punctuation_but_not_identity():
    assert SYN.taxon_key("[Ruminococcus] gnavus") == SYN.taxon_key("_Ruminococcus_gnavus")
    assert SYN.taxon_key("Clostridium sp. 7_2_43FAA") == SYN.taxon_key("Clostridium_sp_7_2_43FAA")
    assert SYN.taxon_key("Bacteroides fragilis") != SYN.taxon_key("Bacteroides ovatus")


def test_a_flagged_organism_with_no_transporter_is_named_not_omitted():
    """The reassuring case has to be stated, not left to absence."""
    index = {"ruminococcus gnavus": "[Ruminococcus] gnavus",
             "phocaeicola vulgatus": "Phocaeicola vulgatus"}
    contexts = SYN.organism_contexts(_results_with_gnavus(), index)
    arm = SIM.ScenarioResult(
        scenario_id="inulin.both", arm="both", execution_state="solved_model",
        solver_status="optimal", raw_flux=200.0, community_growth=0.05,
        flux_per_growth=4000.0, unavailable_reason=None,
        target_exchange=SIM.BUTYRATE_EXCHANGE, medium_id="european",
        coverage=None, perturbation=None,
    )
    arm.community_taxa = ("_Ruminococcus_gnavus", "Phocaeicola_vulgatus")
    arm.taxon_substrate_share = {
        "Phocaeicola_vulgatus": SIM.MemberRange.of_uptake(-6.0, 0.0, 6.14),
    }
    substrate = SIM.Substrate("inulin", "Inulin", "EX_inulin_m", 6.14)
    cand = SYN.Candidate("inulin", substrate, "european", (arm,),
                         SIM.contrasts(1.0, 2.0, 3.0, 4.0), "pair_gains")
    reachable = {t for t, r in arm.taxon_substrate_share.items() if r.possible}
    keys = {SYN.taxon_key(t) for t in reachable}
    cand.flagged_without_route = [
        c for t in arm.community_taxa
        if (c := SYN.context_for(contexts, t)) is not None and c.flagged
        and SYN.taxon_key(t) not in keys
    ]
    assert [c.display for c in cand.flagged_without_route] == ["Mediterraneibacter gnavus"]
    assert cand.to_json()["flagged_without_route"][0]["display"] == "Mediterraneibacter gnavus"
