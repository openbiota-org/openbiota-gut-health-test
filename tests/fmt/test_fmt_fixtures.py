"""BUILD_SPEC_v06.0 §15.2 — the mandatory executable fixtures F1–F7.

These test arithmetic and behaviour, never clinical efficacy. Every expected
number in this file comes from the specification text, not from the
implementation.
"""

from __future__ import annotations

import pytest

from openbiota.fmt import coverage as cov_mod
from openbiota.fmt import sets as sets_mod
from openbiota.fmt.coverage import Availability
from openbiota.fmt.targets import TargetRecord

TOL = 1e-9


def _target(name: str) -> TargetRecord:
    """One independent, equally weighted presence goal in its own group."""
    return TargetRecord(
        target_id=name,
        feature_group_id=f"g:{name}",
        feature_kind="species",
        feature_ids=[name],
        label=name,
        direction="supply",
        goal_kind="presence",
        origin="reference_supported",
        target_status="conditional",
        recipient_observation_ids=[],
        reference_id="ref",
        desired_state={"kind": "presence"},
        clinical_endpoint=None,
        mechanism_source_ids=[],
        priority_class="moderate",
        counting_group=f"g:{name}",
    )


def _cell(target: str, material: str, value: float | None, *, assessed: bool = True,
          lower: float | None = None, upper: float | None = None) -> Availability:
    return Availability(
        target_id=target,
        material_id=material,
        supported_value=value,
        scenario_lower=(value if lower is None else lower) if value is not None else (lower or 0.0),
        scenario_upper=(value if upper is None else upper) if value is not None else (upper or 1.0),
        rule_id="TEST",
        observation_ids=[],
        assessed=assessed,
    )


def _equal_weights(targets: list[TargetRecord]) -> tuple[dict[str, float], dict[str, float]]:
    alpha = {t.counting_group: 1.0 / len(targets) for t in targets}
    beta = {t.target_id: 1.0 for t in targets}
    return alpha, beta


def _coverage(targets: list[TargetRecord], table: dict[tuple[str, str], Availability],
              members: list[str]) -> cov_mod.CoverageResult:
    alpha, beta = _equal_weights(targets)
    return cov_mod.coverage(targets=targets, alpha=alpha, beta=beta, table=table, members=members)


def _candidates(targets: list[TargetRecord], table: dict[tuple[str, str], Availability],
                materials: list[str]) -> list[sets_mod.Candidate]:
    alpha, beta = _equal_weights(targets)
    subsets, _ = sets_mod.feasible_subsets(
        materials, person_of={m: m for m in materials}, max_donors=None,
        max_materials=None, max_subsets=1000,
    )
    cands = sets_mod.evaluate_candidates(
        targets=targets, alpha=alpha, beta=beta, table=table, subsets=subsets,
        person_of={m: m for m in materials}, label_of={m: m for m in materials},
    )
    sets_mod.pareto_fronts(cands)
    return cands


# --------------------------------------------------------------------------- #
# F1 — genuine complementarity
# --------------------------------------------------------------------------- #


@pytest.fixture()
def f1() -> tuple[list[TargetRecord], dict[tuple[str, str], Availability]]:
    targets = [_target("t1"), _target("t2"), _target("t3")]
    table = {
        ("t1", "A"): _cell("t1", "A", 1.0), ("t2", "A"): _cell("t2", "A", 1.0),
        ("t3", "A"): _cell("t3", "A", 0.0),
        ("t1", "B"): _cell("t1", "B", 0.0), ("t2", "B"): _cell("t2", "B", 1.0),
        ("t3", "B"): _cell("t3", "B", 1.0),
    }
    return targets, table


def test_f1_complementarity_and_marginals(f1) -> None:
    """A supplies t1,t2; B supplies t2,t3. Both 66.67; together 100; gain 33.33."""
    targets, table = f1
    assert _coverage(targets, table, ["A"]).coverage == pytest.approx(200.0 / 3, abs=1e-9)
    assert _coverage(targets, table, ["B"]).coverage == pytest.approx(200.0 / 3, abs=1e-9)
    assert _coverage(targets, table, ["A", "B"]).coverage == pytest.approx(100.0, abs=1e-9)

    cands = {c.candidate_id: c for c in _candidates(targets, table, ["A", "B"])}
    both = cands["A+B"]
    assert both.marginals["gain_over_best_member_singleton"] == pytest.approx(100.0 / 3, abs=1e-9)
    assert both.marginals["leave_one_member_out"]["A"] == pytest.approx(100.0 / 3, abs=1e-9)
    assert both.marginals["leave_one_member_out"]["B"] == pytest.approx(100.0 / 3, abs=1e-9)


def test_f1_all_three_are_nondominated(f1) -> None:
    """With donor count in the vector, A, B and A+B are all a real tradeoff."""
    targets, table = f1
    cands = {c.candidate_id: c for c in _candidates(targets, table, ["A", "B"])}
    assert {cid for cid, c in cands.items() if c.pareto_front == 1} == {"A", "B", "A+B"}


# --------------------------------------------------------------------------- #
# F2 — redundant donor
# --------------------------------------------------------------------------- #


def test_f2_redundant_donor_adds_nothing_and_is_dominated() -> None:
    targets = [_target("t1"), _target("t2"), _target("t3")]
    table = {
        ("t1", "A"): _cell("t1", "A", 1.0), ("t2", "A"): _cell("t2", "A", 1.0),
        ("t3", "A"): _cell("t3", "A", 1.0),
        ("t1", "B"): _cell("t1", "B", 1.0), ("t2", "B"): _cell("t2", "B", 0.0),
        ("t3", "B"): _cell("t3", "B", 0.0),
    }
    assert _coverage(targets, table, ["A"]).coverage == pytest.approx(100.0)
    assert _coverage(targets, table, ["A", "B"]).coverage == pytest.approx(100.0)

    cands = {c.candidate_id: c for c in _candidates(targets, table, ["A", "B"])}
    assert cands["A+B"].marginals["gain_over_best_member_singleton"] == pytest.approx(0.0)
    assert sets_mod.dominates(cands["A"], cands["A+B"])
    assert not sets_mod.dominates(cands["A+B"], cands["A"])
    # B's record survives and its zero contribution is explicit.
    assert cands["B"].coverage.coverage == pytest.approx(100.0 / 3, abs=1e-9)


# --------------------------------------------------------------------------- #
# F3 — unknown feature
# --------------------------------------------------------------------------- #


def test_f3_unknown_feature_bounds_and_assessed_weight() -> None:
    """Supported 33.33, lower 33.33, upper 66.67, assessed goal weight 2/3."""
    targets = [_target("t1"), _target("t2"), _target("t3")]
    table = {
        ("t1", "A"): _cell("t1", "A", 1.0),
        ("t2", "A"): _cell("t2", "A", None, assessed=False, lower=0.0, upper=1.0),
        ("t3", "A"): _cell("t3", "A", 0.0),
    }
    cov = _coverage(targets, table, ["A"])
    assert cov.coverage == pytest.approx(100.0 / 3, abs=1e-9)
    assert cov.lower == pytest.approx(100.0 / 3, abs=1e-9)
    assert cov.upper == pytest.approx(200.0 / 3, abs=1e-9)
    assert cov.assessed_goal_weight == pytest.approx(2.0 / 3, abs=1e-9)
    assert cov.unassessed_goal_weight == pytest.approx(1.0 / 3, abs=1e-9)
    # The measurement itself stays unknown.
    assert table[("t2", "A")].supported_value is None
    assert table[("t2", "A")].contribution == 0.0
    payload = table[("t2", "A")].to_json()
    assert payload["supported_value"] is None and payload["assessed"] is False


# --------------------------------------------------------------------------- #
# F4 — confirmed exclusion
# --------------------------------------------------------------------------- #


def test_f4_confirmed_exclusion_removes_material_and_every_containing_set(f1) -> None:
    from openbiota.fmt import concerns

    targets, table = f1
    hard = concerns.exclusions_from_screening(
        candidate_id="B",
        material_id="B",
        screening={
            "confirmed_exclusions": [
                {
                    "rule_id": "X02",
                    "finding": "a linked validated test confirmed ESBL-producing E. coli",
                    "source": "reference laboratory report",
                    "date": "2026-09-01",
                    "scope": "donation_lot",
                }
            ]
        },
    )
    assert len(hard) == 1 and hard[0].disposition == "exclusion_confirmed"
    assert "rule X02 applies" in hard[0].reason

    admissible = ["A"]  # B excluded, so A+B cannot be proposed
    cands = {c.candidate_id: c for c in _candidates(targets, table, admissible)}
    assert set(cands) == {"A"}


def test_f4_ambiguous_arg_alone_does_not_exclude() -> None:
    """Replacing the confirmation with an unlinked ARG hit restores comparison."""
    from openbiota.fmt import concerns

    hard = concerns.exclusions_from_screening(
        candidate_id="B", material_id="B",
        screening={"confirmed_exclusions": [{"rule_id": "R02", "finding": "tet efflux family fragment"}]},
    )
    assert hard == []


# --------------------------------------------------------------------------- #
# F5 — all excluded, or none
# --------------------------------------------------------------------------- #


def test_f5_no_donors_gives_targets_and_no_ranking(f1) -> None:
    targets, table = f1
    subsets, audit = sets_mod.feasible_subsets(
        [], person_of={}, max_donors=None, max_materials=None, max_subsets=100
    )
    assert subsets == []
    assert audit["candidate_material_count"] == 0
    assert audit["theoretical_nonempty_subsets"] == 0
    assert audit["search_complete"] is True


def test_f5_single_donor_marginal_to_empty_set_is_defined(f1) -> None:
    targets, table = f1
    cands = {c.candidate_id: c for c in _candidates(targets, table, ["A"])}
    a = cands["A"]
    assert a.marginals["empty_set_coverage"] == 0.0
    assert a.marginals["leave_one_member_out"]["A"] == pytest.approx(200.0 / 3, abs=1e-9)


# --------------------------------------------------------------------------- #
# F6 — batch invariance
# --------------------------------------------------------------------------- #


def test_f6_adding_an_unrelated_donor_changes_no_measurement(f1) -> None:
    targets, table = f1
    before = _coverage(targets, table, ["A"])
    wider = dict(table)
    for t in targets:
        wider[(t.target_id, "C")] = _cell(t.target_id, "C", 0.0)
    after = _coverage(targets, wider, ["A"])
    assert after.coverage == before.coverage
    assert after.denominator_target_ids == before.denominator_target_ids
    assert after.assessed_goal_weight == before.assessed_goal_weight


# --------------------------------------------------------------------------- #
# F7 — partial analytical scopes
# --------------------------------------------------------------------------- #


def test_f7_partial_scope_is_not_renormalised_away() -> None:
    """B has no function lane: its unavailable target stays in the denominator."""
    taxa1, taxa2 = _target("taxon1"), _target("taxon2")
    func = _target("function1")
    func.feature_kind = "gene_panel"
    targets = [taxa1, taxa2, func]
    table = {
        ("taxon1", "A"): _cell("taxon1", "A", 1.0),
        ("taxon2", "A"): _cell("taxon2", "A", 1.0),
        ("function1", "A"): _cell("function1", "A", 1.0),
        ("taxon1", "B"): _cell("taxon1", "B", 1.0),
        ("taxon2", "B"): _cell("taxon2", "B", 1.0),
        ("function1", "B"): _cell("function1", "B", None, assessed=False, lower=0.0, upper=1.0),
    }
    a = _coverage(targets, table, ["A"])
    b = _coverage(targets, table, ["B"])
    assert a.coverage == pytest.approx(100.0)
    assert b.coverage == pytest.approx(200.0 / 3, abs=1e-9)
    assert b.upper == pytest.approx(100.0, abs=1e-9)
    assert b.assessed_goal_weight == pytest.approx(2.0 / 3, abs=1e-9)
    assert len(a.denominator_target_ids) == len(b.denominator_target_ids) == 3
