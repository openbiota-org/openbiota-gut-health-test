"""Longitudinal views (A13) and the recovery envelope (§6.4).

The distances are the easy part. These tests concentrate on the three ways
a timeline lies: comparing across a change of method, counting a
reprocessing run as a new timepoint, and quoting a recovery percentage
from a displacement no larger than the noise.
"""

from __future__ import annotations

from datetime import date, timedelta

import pytest

from openbiota.extension import longitudinal as L

D0 = date(2026, 1, 1)
LANE = {"profiler": "metaphlan4", "database_version": "1", "specimen_type": "stool"}


def _fp(sample_id: str, *, specimen: str | None = None, **overrides) -> L.Fingerprint:
    return L.Fingerprint(
        sample_id, "P1", specimen_id=specimen or sample_id, fields={**LANE, **overrides},
    )


def _tp(sample_id: str, day: int, abundances: dict[str, float], *,
        specimen: str | None = None, events: tuple[str, ...] = ()) -> L.TimePoint:
    return L.TimePoint(
        _fp(sample_id, specimen=specimen), D0 + timedelta(days=day), abundances, events=events,
    )


# --------------------------------------------------------------------------- #
# distances
# --------------------------------------------------------------------------- #


def test_bray_curtis_matches_its_definition():
    p, q = {"a": 3.0, "b": 1.0}, {"a": 1.0, "c": 2.0}
    # |3-1| + |1-0| + |0-2| = 5 over (3+1) + (1+0) + (0+2) = 7
    assert L.bray_curtis(p, q) == pytest.approx(5 / 7)
    assert L.bray_curtis(p, p) == 0.0
    assert L.bray_curtis({"a": 1.0}, {"b": 1.0}) == 1.0


def test_two_absences_have_no_distance_rather_than_a_distance_of_zero():
    assert L.bray_curtis({}, {}) is None
    assert L.bray_curtis({"a": 0.0}, {"b": 0.0}) is None


def test_jaccard_uses_a_fixed_floor_so_depth_does_not_become_richness():
    shallow = {"a": 50.0, "b": 49.0, "c": 1.0}
    deep = {"a": 50.0, "b": 49.0, "c": 1.0, "trace": 0.02}
    # `trace` is 0.02% of the sample, under the 0.1% floor: not a new organism.
    assert L.presence_set(deep) == L.presence_set(shallow)
    assert L.jaccard_distance(shallow, deep) == 0.0


def test_an_empty_union_is_undefined_not_perfect_agreement():
    assert L.jaccard_distance({}, {}) is None
    assert L.jaccard_distance({"a": 0.0}, {}) is None


def test_jaccard_counts_membership_not_amount():
    a = {"x": 90.0, "y": 10.0}
    b = {"x": 10.0, "y": 90.0}
    assert L.jaccard_distance(a, b) == 0.0, "the same organisms, in different amounts"
    assert L.bray_curtis(a, b) == pytest.approx(0.8), "but the abundances moved a long way"


def test_closure_and_mean_composition_sum_to_one():
    closed = L.closure({"a": 2.0, "b": 2.0, "c": 0.0})
    assert closed == {"a": 0.5, "b": 0.5}
    mean = L.mean_composition([{"a": 1.0}, {"b": 1.0}])
    assert sum(mean.values()) == pytest.approx(1.0)
    assert mean["a"] == pytest.approx(0.5)


# --------------------------------------------------------------------------- #
# comparability
# --------------------------------------------------------------------------- #


def test_the_same_lane_is_natively_comparable():
    state, why = L.comparability(_fp("s1"), _fp("s2"))
    assert state == "native" and "same lane" in why


def test_a_database_upgrade_is_not_a_change_in_the_person():
    state, why = L.comparability(_fp("s1"), _fp("s2", database_version="2"))
    assert state == "descriptive_only"
    assert "change of method as much as any change in you" in why


@pytest.mark.parametrize("field", L.BLOCKING_FIELDS)
def test_every_blocking_field_actually_blocks(field):
    state, _ = L.comparability(_fp("s1"), _fp("s2", **{field: "something else"}))
    assert state == "descriptive_only", f"{field} should break native comparability"


@pytest.mark.parametrize("field", ("stool_form", "storage", "read_depth"))
def test_context_fields_are_recorded_without_blocking(field):
    state, _ = L.comparability(_fp("s1"), _fp("s2", **{field: "different"}))
    assert state == "native"


def test_reprocessing_together_or_a_bridge_restores_a_comparison():
    a, b = _fp("s1"), _fp("s2", database_version="2")
    assert L.comparability(a, b, reprocessed_together=True)[0] == "uniformly_reprocessed"
    assert L.comparability(a, b, bridge_validated=True)[0] == "bridge_validated"


def test_every_comparability_state_is_one_of_the_declared_four():
    a, b = _fp("s1"), _fp("s2", profiler="other")
    for kwargs in ({}, {"reprocessed_together": True}, {"bridge_validated": True}):
        assert L.comparability(a, b, **kwargs)[0] in L.COMPARABILITY
    assert L.comparability(a, _fp("s3"))[0] in L.COMPARABILITY


def test_a_reprocessed_sample_keeps_its_specimen_identity():
    original = _fp("run1", specimen="STOOL-A")
    rerun = _fp("run2", specimen="STOOL-A")
    assert L.is_repeat_processing(original, rerun)
    assert not L.is_repeat_processing(original, _fp("run3", specimen="STOOL-B"))


# --------------------------------------------------------------------------- #
# the series
# --------------------------------------------------------------------------- #


def test_reprocessing_does_not_become_a_second_timepoint():
    series = L.Series.build("P1", [
        _tp("run1", 0, {"a": 1.0}, specimen="STOOL-A"),
        _tp("run2", 30, {"a": 1.0}, specimen="STOOL-A"),
        _tp("run3", 60, {"a": 1.0}, specimen="STOOL-B"),
    ])
    assert series.n_timepoints == 2
    assert series.repeat_processing == ["run2"]


def test_points_are_ordered_by_collection_date_not_by_arrival():
    series = L.Series.build("P1", [
        _tp("late", 60, {"a": 1.0}), _tp("early", 0, {"a": 1.0}), _tp("mid", 30, {"a": 1.0}),
    ])
    assert [p.sample_id for p in series.points] == ["early", "mid", "late"]
    assert [p.sample_id for p in series.last(2)] == ["mid", "late"]


def test_a_version_change_starts_a_new_segment():
    series = L.Series.build("P1", [
        _tp("s1", 0, {"a": 1.0}), _tp("s2", 30, {"a": 1.0}),
        L.TimePoint(_fp("s3", database_version="2"), D0 + timedelta(days=60), {"a": 1.0}),
    ])
    segments = series.segments()
    assert len(segments) == 2
    assert segments[0]["sample_ids"] == ["s1", "s2"]
    assert segments[1]["sample_ids"] == ["s3"]


def test_no_distance_is_quoted_across_a_segment_boundary():
    series = L.Series.build("P1", [
        _tp("s1", 0, {"a": 60.0, "b": 40.0}),
        L.TimePoint(_fp("s2", profiler="other"), D0 + timedelta(days=30), {"a": 20.0, "b": 80.0}),
    ])
    comparison = series.pairwise()[0]
    assert comparison["comparability"] == "descriptive_only"
    assert comparison["bray_curtis"] is None
    assert comparison["jaccard_distance"] is None


def test_a_native_pair_carries_its_distances_and_metric_changes():
    a = L.TimePoint(_fp("s1"), D0, {"a": 60.0, "b": 40.0}, metrics={"shannon": 2.0})
    b = L.TimePoint(_fp("s2"), D0 + timedelta(days=30), {"a": 40.0, "b": 60.0},
                    metrics={"shannon": 2.5})
    comparison = L.Series.build("P1", [a, b]).pairwise()[0]
    assert comparison["comparability"] == "native"
    assert comparison["bray_curtis"] == pytest.approx(0.2)
    assert comparison["days_apart"] == 30
    assert comparison["metric_changes"]["shannon"] == pytest.approx(0.5)


# --------------------------------------------------------------------------- #
# recovery
# --------------------------------------------------------------------------- #


def _antibiotic_series() -> L.Series:
    return L.Series.build("P1", [
        _tp("b1", 0, {"Bacteroides": 40, "Faecalibacterium": 30, "Roseburia": 20, "Escherichia": 1}),
        _tp("b2", 7, {"Bacteroides": 42, "Faecalibacterium": 28, "Roseburia": 21, "Escherichia": 1}),
        _tp("b3", 14, {"Bacteroides": 39, "Faecalibacterium": 31, "Roseburia": 19, "Escherichia": 2}),
        _tp("e1", 21, {"Bacteroides": 20, "Faecalibacterium": 3, "Roseburia": 1, "Escherichia": 40},
            events=("amoxicillin course",)),
        _tp("a1", 40, {"Bacteroides": 30, "Faecalibacterium": 12, "Roseburia": 8, "Escherichia": 18}),
        _tp("a2", 70, {"Bacteroides": 37, "Faecalibacterium": 24, "Roseburia": 16, "Escherichia": 5}),
    ])


def test_recovery_measures_return_from_the_observed_peak():
    out = L.recovery(_antibiotic_series(), event_sample_id="e1")
    assert out.state == "measured"
    assert out.peak_sample == "e1"
    assert out.peak_distance > out.envelope, "displacement must clear the baseline's own spread"
    returns = {sample: ret for sample, _distance, ret in out.trajectory}
    assert returns["e1"] == pytest.approx(0.0)
    assert 50 < returns["a1"] < 70
    assert returns["a2"] > 85
    assert all(0.0 <= r <= 100.0 for r in returns.values())


def test_recovery_needs_three_baselines_and_says_so_without_hiding_the_rest():
    series = L.Series.build("P1", [
        _tp("b1", 0, {"a": 60.0, "b": 40.0}),
        _tp("b2", 7, {"a": 61.0, "b": 39.0}),
        _tp("e1", 14, {"a": 10.0, "b": 90.0}, events=("antibiotics",)),
    ])
    out = L.recovery(series, event_sample_id="e1")
    assert out.state == "insufficient_baseline"
    assert "at least 3" in out.reason
    assert "Changes and similarities are still shown" in out.reason
    assert out.to_json()["trajectory"] == []


def test_a_wobble_inside_the_baseline_spread_is_not_a_displacement():
    series = L.Series.build("P1", [
        _tp("b1", 0, {"a": 60.0, "b": 40.0}),
        _tp("b2", 7, {"a": 50.0, "b": 50.0}),
        _tp("b3", 14, {"a": 70.0, "b": 30.0}),
        _tp("e1", 21, {"a": 59.0, "b": 41.0}, events=("probiotic",)),
    ])
    out = L.recovery(series, event_sample_id="e1")
    assert out.state == "no_displacement"
    assert "within the baseline's own spread" in out.reason
    assert all(r is None for _s, _d, r in out.trajectory)


def test_a_technical_noise_bound_can_only_widen_the_envelope():
    series = _antibiotic_series()
    tight = L.recovery(series, event_sample_id="e1")
    wide = L.recovery(series, event_sample_id="e1", technical_noise_bound=0.9)
    assert wide.envelope == pytest.approx(0.9)
    assert wide.envelope > tight.envelope
    assert wide.state == "no_displacement", "a wide enough noise bound explains the move away"


def test_points_before_the_peak_show_displacement_not_return():
    series = L.Series.build("P1", [
        _tp("b1", 0, {"a": 60.0, "b": 40.0}),
        _tp("b2", 7, {"a": 61.0, "b": 39.0}),
        _tp("b3", 14, {"a": 59.0, "b": 41.0}),
        _tp("e1", 21, {"a": 40.0, "b": 60.0}, events=("antibiotics",)),
        _tp("p1", 28, {"a": 5.0, "b": 95.0}),   # still getting worse: this is the peak
        _tp("p2", 60, {"a": 55.0, "b": 45.0}),
    ])
    out = L.recovery(series, event_sample_id="e1")
    assert out.state == "measured" and out.peak_sample == "p1"
    by_sample = {s: r for s, _d, r in out.trajectory}
    assert by_sample["e1"] is None, "before the peak there is no return to measure"
    assert by_sample["p1"] == pytest.approx(0.0)
    assert by_sample["p2"] > 80


def test_recovery_never_claims_the_old_community_was_healthy():
    text = " ".join(L.recovery(_antibiotic_series(), event_sample_id="e1").to_json()["limitations"])
    assert "does not establish that the earlier community was healthy" in text
    assert "not a 95% confidence interval" in text


def test_an_unknown_event_sample_is_refused_rather_than_guessed():
    assert L.recovery(_antibiotic_series(), event_sample_id="nope").state == "unavailable"


# --------------------------------------------------------------------------- #
# follow-up and the single-sample view
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("events", "expected"), [
    ((), "routine_review"),
    (("amoxicillin 500mg course",), "recent_antibiotics"),
    (("started a probiotic",), "new_probiotic_or_prebiotic"),
    (("switched to a high-fibre diet",), "major_dietary_change"),
])
def test_each_exposure_picks_its_own_follow_up_question(events, expected):
    chosen = [t["template_id"] for t in L.follow_up(list(events))]
    assert expected in chosen


def test_no_single_universal_retesting_interval_is_hardcoded():
    intervals = {t["interval"] for t in L.FOLLOW_UP_TEMPLATES}
    assert len(intervals) == len(L.FOLLOW_UP_TEMPLATES), "each template needs its own interval"
    for template in L.FOLLOW_UP_TEMPLATES:
        assert template["question"] and template["rationale"]


def test_the_single_sample_view_is_shown_rather_than_suppressed():
    view = L.not_enough_history("S1")
    assert view["state"] == "single_timepoint"
    assert view["what_a_second_sample_adds"]
    assert view["what_makes_them_comparable"]
    assert view["follow_up"]
    assert "one day" in " ".join(view["limitations"])


def test_percentile_change_is_kept_apart_from_a_raw_change():
    out = L.percentile_change(40.0, 60.0)
    assert out["delta_percentile"] == 20.0
    assert out["unit"] == "percentile points"
    assert "not a change in the underlying amount" in out["note"]
    assert L.percentile_change(None, 60.0)["state"] == "unavailable"


def test_stability_needs_more_than_one_of_your_own_results():
    assert L.stability([1.0])["state"] == "insufficient_history"
    out = L.stability([1.0, 2.0, 3.0])
    assert out["state"] == "measured" and out["spread"] == 2.0
    assert "unusual for you" in out["note"]


def test_summarise_finds_the_recorded_event_without_being_told():
    view = L.summarise(_antibiotic_series())
    assert view["recovery"]["state"] == "measured"
    assert view["recovery"]["peak_sample"] == "e1"
    assert "recent_antibiotics" in [t["template_id"] for t in view["follow_up"]]
    assert view["presence_floor"] == L.PRESENCE_FLOOR
    assert len(view["comparisons"]) == view["n_timepoints"] - 1


def test_summarise_without_an_event_reports_no_recovery_rather_than_zero():
    series = L.Series.build("P1", [_tp("s1", 0, {"a": 1.0}), _tp("s2", 30, {"a": 1.0})])
    view = L.summarise(series)
    assert view["recovery"] is None
    assert view["n_timepoints"] == 2


# --------------------------------------------------------------------------- #
# the consolidated action planner (A12)
# --------------------------------------------------------------------------- #


def _option(**kw):
    from openbiota.extension import planner as P

    base = {"option_id": "x", "label": "X", "category": "food"}
    return P.Option(**{**base, **kw})


def test_options_are_ordered_by_findings_then_evidence_then_identity():
    from openbiota.extension import planner as P

    broad = _option(option_id="broad", addresses=("a", "b", "c"), best_lane="C")
    flagged = _option(option_id="flagged", addresses=("a",), addresses_flagged=("a",),
                      best_lane="E")
    strong = _option(option_id="strong", addresses=("b",), best_lane="A")
    plan = P.Plan("S1", [strong, broad, flagged])
    assert [o.option_id for o in plan.ordered] == ["flagged", "broad", "strong"]
    # ...and the same list, whatever order they arrived in.
    assert [o.option_id for o in P.Plan("S1", [broad, flagged, strong]).ordered] == \
           [o.option_id for o in plan.ordered]


def test_two_identical_options_never_swap_between_runs():
    from openbiota.extension import planner as P

    a = _option(option_id="aaa", addresses=("x",), best_lane="C")
    b = _option(option_id="bbb", addresses=("x",), best_lane="C")
    assert [o.option_id for o in P.Plan("S1", [b, a]).ordered] == ["aaa", "bbb"]


def test_start_here_holds_at_most_three_and_prefers_different_work():
    from openbiota.extension import planner as P

    options = [
        _option(option_id=f"p{i}", category="probiotic", addresses=("same",), best_lane="A")
        for i in range(5)
    ]
    options.append(_option(option_id="food1", category="food", addresses=("other",),
                           best_lane="C"))
    plan = P.Plan("S1", options)
    picked = plan.start_here
    assert len(picked) <= P.START_HERE_MAX
    assert "food1" in {o.option_id for o in picked}, (
        "a second idea beats a fifth probiotic addressing the same finding"
    )


def test_a_contraindication_leaves_the_catalogue_intact():
    from openbiota.extension import planner as P

    excluded = _option(
        option_id="nope", addresses=("a",), best_lane="A",
        exclusions=(P.Exclusion("declared_allergy", "you reported a dairy allergy", True),),
    )
    fine = _option(option_id="fine", addresses=("a",), best_lane="C")
    plan = P.Plan("S1", [excluded, fine])
    assert [o.option_id for o in plan.start_here] == ["fine"]
    assert "nope" in {o.option_id for o in plan.ordered}, "it stays in the catalogue"
    assert plan.coverage()["n_excluded_for_this_person"] == 1
    payload = excluded.to_json()["exclusions"][0]
    assert payload["from_supplied_fact"] is True
    assert "evidence intact" in payload["note"]


def test_a_missing_trial_does_not_hide_an_option():
    from openbiota.extension import planner as P

    lab_only = _option(option_id="lab", addresses=("a",), best_lane="E")
    plan = P.Plan("S1", [lab_only])
    assert lab_only.eligible_for_start_here
    assert not lab_only.has_human_evidence
    assert plan.coverage()["n_without_human_evidence"] == 1
    assert "not absence of anything to try" in plan.coverage()["coverage_note"]


def test_every_option_carries_an_answerable_monitoring_question():
    from openbiota.extension import planner as P

    for category, _label in P.CATEGORIES:
        question = P.monitoring_question(category)
        assert question and "?" in question or "clinician" in question
    assert "distinguishes passing through from settling in" in P.monitoring_question("probiotic")
    assert "colonis" not in P.monitoring_question("probiotic").lower().replace(
        "colonisation", ""), "no promise of colonisation"


def test_a_study_exposure_is_never_presented_as_a_dose():
    option = _option(studied_exposure="10 mg/kg in mice for 14 days")
    note = option.to_json()["exposure_note"]
    assert "not a dose to take" in note
    assert "not converted into a human regimen" in note


def test_the_checklist_merges_aliases_without_merging_different_products():
    from openbiota.extension import planner as P

    assert P.normalise_food("cacao") == P.normalise_food("cocoa powder") == "cocoa"
    assert P.normalise_food("cocoa husk") == "cocoa husk", "a feed ingredient is not cocoa powder"
    assert P.normalise_food("beans") == P.normalise_food("pulses") == "legumes"


def test_the_checklist_deduplicates_and_records_what_each_item_is_for():
    from openbiota.extension import planner as P

    a = _option(option_id="a", category="prebiotic", addresses=("low butyrate",),
                foods=("inulin",))
    b = _option(option_id="b", category="prebiotic", addresses=("low bifidobacteria",),
                foods=("inulin", "galacto-oligosaccharides"))
    groups = P.Plan("S1", [a, b]).checklist()
    items = {i["item"]: i for g in groups for i in g["items"]}
    assert set(items) == {"inulin", "galacto-oligosaccharides"}
    assert set(items["inulin"]["from_options"]) == {"a", "b"}
    assert set(items["inulin"]["addresses"]) == {"low butyrate", "low bifidobacteria"}
    assert items["inulin"]["group"] == "fibre_preparations"


def test_an_excluded_option_contributes_nothing_to_the_shopping_list():
    from openbiota.extension import planner as P

    excluded = _option(option_id="x", category="prebiotic", foods=("inulin",),
                       exclusions=(P.Exclusion("declared_intolerance", "reported", True),))
    assert P.Plan("S1", [excluded]).checklist() == []


def test_shopping_items_come_from_components_not_from_study_titles():
    from openbiota.extension import planner as P

    identity = {"components": [
        {"entity_type": "prebiotic", "molecule": "galacto-oligosaccharides"},
        {"entity_type": "diet", "name": "PREDIMED Mediterranean diet"},
        {"entity_type": "guideline_route", "name": "gastroenterology referral"},
    ]}
    assert P.shopping_items(identity) == (("galacto-oligosaccharides", "prebiotic"),)
    assert P.dietary_pattern(identity) == "PREDIMED Mediterranean diet"


def test_flagged_findings_match_on_the_organism_not_the_whole_display():
    from openbiota.extension import planner as P

    displays = ["Faecalibacterium prausnitzii: low", "Akkermansia muciniphila: typical"]
    assert P.addresses_flagged(displays, ["Faecalibacterium prausnitzii"]) == \
           ("Faecalibacterium prausnitzii: low",)
    assert P.addresses_flagged(displays, ["Bacteroides fragilis"]) == ()


def test_unknown_categories_and_lanes_are_refused():
    from openbiota.extension import planner as P

    with pytest.raises(P.PlannerError, match="unknown category"):
        _option(category="magic")
    with pytest.raises(P.PlannerError, match="unknown evidence lane"):
        _option(best_lane="Z")


def test_the_ranking_policy_disclaims_being_a_prediction():
    from openbiota.extension import planner as P

    payload = P.Plan("S1", [_option()]).to_json()
    assert "not a predicted probability of benefit" in payload["ranking_policy"]["note"]
    assert "not a prediction" in payload["start_here_note"]
    limits = " ".join(payload["limitations"])
    assert "not a dose" in limits
    assert "not automatically an effect of either component alone" in limits
    assert "Conflicting and null evidence stays" in limits
