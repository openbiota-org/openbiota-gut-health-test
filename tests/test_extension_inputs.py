"""The input register (A14): what was supplied, assumed, imported and absent.

The register's whole purpose is keeping four things apart that are easy to
merge: a fact, an assumption, an absence and an import. Most of these tests
check that a record which would blur two of them is refused.
"""

from __future__ import annotations

import json

import pytest

from openbiota.extension import inputs as IN

# --------------------------------------------------------------------------- #
# the states, and what each may carry
# --------------------------------------------------------------------------- #


def test_a_default_cannot_carry_a_reported_value():
    """The one failure this register exists to prevent."""
    with pytest.raises(IN.InputRegisterError, match="not something the person said"):
        IN.InputRecord(
            input_id="x", category="diet_tolerance", label="x",
            input_state="assumed_default", reported_value="high fibre",
            effective_value="high fibre", default_rationale="because",
        )


def test_a_default_must_say_why_it_is_the_default():
    with pytest.raises(IN.InputRegisterError, match="must say why"):
        IN.InputRecord(
            input_id="x", category="diet_tolerance", label="x",
            input_state="assumed_default", effective_value="european",
        )


def test_supplied_without_a_value_is_not_supplied():
    with pytest.raises(IN.InputRegisterError, match="is 'not_supplied'"):
        IN.InputRecord(input_id="x", category="laboratory", label="x", input_state="supplied")


def test_not_supplied_cannot_carry_an_effective_value():
    """An absent input with a value in it is a default wearing a disguise."""
    with pytest.raises(IN.InputRegisterError, match="that is a default"):
        IN.InputRecord(
            input_id="x", category="laboratory", label="x",
            input_state="not_supplied", effective_value=0.0,
        )


def test_a_derived_value_must_name_its_derivation():
    with pytest.raises(IN.InputRegisterError, match="must name its derivation"):
        IN.InputRecord(
            input_id="x", category="participant_specimen", label="x",
            input_state="derived", reported_value=None, effective_value="bmi 24",
        )


def test_a_conflict_records_both_sides_rather_than_picking_one():
    with pytest.raises(IN.InputRegisterError, match="records that disagree"):
        IN.InputRecord(
            input_id="x", category="clinical_history", label="x",
            input_state="conflicting", conflict=("only one",),
        )
    record = IN.InputRecord(
        input_id="x", category="clinical_history", label="Age",
        input_state="conflicting", conflict=("manifest says 41", "questionnaire says 43"),
    )
    assert len(record.to_json()["conflict"]) == 2
    assert not record.is_fact


def test_unknown_states_and_categories_are_refused():
    with pytest.raises(IN.InputRegisterError, match="unknown input state"):
        IN.InputRecord(input_id="x", category="laboratory", label="x", input_state="probably")
    with pytest.raises(IN.InputRegisterError, match="unknown category"):
        IN.InputRecord(input_id="x", category="astrology", label="x", input_state="not_supplied")


def test_only_supplied_and_derived_read_as_facts():
    def record(state: str, **kw):
        return IN.InputRecord(input_id="x", category="clinical_history", label="x",
                              input_state=state, **kw)
    assert record("supplied", reported_value=True).is_fact
    assert record("derived", effective_value=1, derivation="from height and weight").is_fact
    assert not record("assumed_default", effective_value=False,
                      default_rationale="not supplied").is_fact
    assert not record("inherited_default", effective_value=False,
                      default_rationale="not supplied").is_fact
    assert not record("not_supplied").is_fact


def test_a_value_with_no_consumer_says_so():
    record = IN.InputRecord(
        input_id="x", category="laboratory", label="Stool pH",
        input_state="supplied", reported_value=6.4, effective_value=6.4, unit="pH",
    )
    assert not record.is_used
    assert record.to_json()["usage_note"] == "Recorded; not used in this report"
    used = IN.InputRecord(
        input_id="y", category="laboratory", label="Calprotectin",
        input_state="supplied", reported_value=30, effective_value=30, unit="ug/g",
        affects=(IN.Influence("report.patterns", "lab import", "adds_context", "shown beside"),),
    )
    assert used.is_used and used.to_json()["usage_note"] is None


def test_an_influence_that_changes_nothing_does_not_count_as_use():
    record = IN.InputRecord(
        input_id="x", category="laboratory", label="x", input_state="supplied",
        reported_value=1, effective_value=1,
        affects=(IN.Influence("nothing", "none", "none", "stored only"),),
    )
    assert not record.is_used


def test_unknown_influence_kinds_are_refused():
    with pytest.raises(IN.InputRegisterError, match="unknown influence kind"):
        IN.Influence("t", "r", "magically_improves", "e")


# --------------------------------------------------------------------------- #
# laboratory results
# --------------------------------------------------------------------------- #


def _lab(**kw) -> IN.LabResult:
    base = {
        "analyte_id": "calprotectin", "panel_id": "inflammation",
        "original_name": "Faecal calprotectin", "specimen": "stool",
        "value": 42.0, "original_unit": "ug/g",
    }
    return IN.LabResult(**{**base, **kw})


def test_a_censored_result_keeps_its_censoring():
    """`<50` is neither the number 50 nor the number 0."""
    result = _lab(value=50.0, censoring="<")
    assert result.censored
    assert result.display_value == "<50 ug/g"
    assert result.to_json()["censoring"] == "<"


def test_a_censored_result_is_classified_only_where_the_bound_settles_it():
    """A censored result decides the question sometimes, and often does not."""
    # "<50" with an upper bound of 50 or more: the true value is under it either way.
    assert _lab(value=50.0, censoring="<", reference_high=50.0).in_reference_interval is True
    assert _lab(value=50.0, censoring="<", reference_high=200.0).in_reference_interval is True
    # "<50" with a lower bound too: could be 5 (below) or 30 (inside). Unanswerable.
    assert _lab(value=50.0, censoring="<", reference_low=10.0,
                reference_high=200.0).in_reference_interval is None
    # "<50" where even the ceiling sits under the floor: definitely below.
    assert _lab(value=50.0, censoring="<", reference_low=100.0).in_reference_interval is False
    # ">1000" against a top of 200: definitely outside.
    assert _lab(value=1000.0, censoring=">", reference_high=200.0).in_reference_interval is False
    # ">1000" with only a floor of 200: definitely above it, and nothing caps it.
    assert _lab(value=1000.0, censoring=">", reference_low=200.0).in_reference_interval is True


def test_no_interval_means_no_classification_not_normal():
    assert _lab().in_reference_interval is None
    assert _lab(reference_low=0.0, reference_high=50.0).in_reference_interval is True
    assert _lab(value=90.0, reference_high=50.0).in_reference_interval is False


def test_a_qualitative_result_needs_no_value_but_a_value_needs_a_unit():
    qualitative = IN.LabResult(
        analyte_id="pcr", panel_id="molecular_pathogen", original_name="C. difficile PCR",
        specimen="stool", qualitative="not_detected", target="tcdB", nucleic_acid="DNA",
    )
    assert qualitative.display_value == "not detected"
    with pytest.raises(IN.InputRegisterError, match="needs its original unit"):
        _lab(original_unit=None)
    with pytest.raises(IN.InputRegisterError, match="either a value or a qualitative"):
        IN.LabResult(analyte_id="x", panel_id="digestion", original_name="x", specimen="stool")


def test_a_result_must_name_its_specimen():
    with pytest.raises(IN.InputRegisterError, match="must name its specimen"):
        _lab(specimen="")


def test_another_body_site_is_never_read_as_a_stool_measurement():
    vaginal = _lab(specimen="vaginal swab", panel_id="molecular_pathogen",
                   value=None, qualitative="detected", original_unit=None)
    allowed, why = IN.link_is_allowed(vaginal, sample_specimen="stool")
    assert not allowed
    assert "vaginal swab" in why and "not compared" in why


def test_a_distant_collection_date_is_labelled_historical():
    old = _lab(collected_at="2024-01-01")
    allowed, why = IN.link_is_allowed(old, sample_collected_at="2026-09-01")
    assert not allowed and "historical" in why
    recent = _lab(collected_at="2026-08-20")
    allowed, why = IN.link_is_allowed(recent, sample_collected_at="2026-09-01")
    assert allowed and "comparable" in why


def test_unknown_panels_censoring_and_qualitative_values_are_refused():
    with pytest.raises(IN.InputRegisterError, match="unknown lab panel"):
        _lab(panel_id="astrology")
    with pytest.raises(IN.InputRegisterError, match="unknown censoring"):
        _lab(censoring="~")
    with pytest.raises(IN.InputRegisterError, match="unknown qualitative"):
        _lab(value=None, original_unit=None, qualitative="probably")


def test_document_hash_is_stable_and_prefixed():
    a, b = IN.document_hash("report text"), IN.document_hash(b"report text")
    assert a == b and a.startswith("sha256:")
    assert a != IN.document_hash("other text")


# --------------------------------------------------------------------------- #
# the register
# --------------------------------------------------------------------------- #


def test_an_empty_register_says_nothing_was_supplied_rather_than_nothing_is_wrong():
    register = IN.InputRegister(sample_id="S1")
    assert "sequencing data alone" in register.headline()
    assert register.summary()["external_laboratory_results"] == "Not supplied"


def test_every_supported_panel_is_listed_even_when_none_was_supplied():
    register = IN.InputRegister(sample_id="S1")
    panels = register.lab_panel_status()
    assert len(panels) == len(IN.LAB_PANELS)
    assert all(not p["supplied"] and p["n_results"] == 0 for p in panels)
    assert {p["panel_id"] for p in panels} >= {"scfa", "inflammation", "microbial_load"}


def test_a_supplied_panel_reports_its_results_and_destination():
    register = IN.InputRegister(sample_id="S1", lab_results=[_lab()])
    supplied = [p for p in register.lab_panel_status() if p["supplied"]]
    assert len(supplied) == 1
    assert supplied[0]["panel_id"] == "inflammation"
    assert supplied[0]["destination"] == "report.patterns_detail"
    assert "never an automatic diagnosis" in supplied[0]["relation"]
    assert "1 external laboratory result was imported" in register.headline()


def test_missing_optional_inputs_are_grouped_not_listed_one_per_page():
    register = IN.InputRegister(sample_id="S1")
    for i in range(6):
        register.add(IN.InputRecord(
            input_id=f"diet.{i}", category="diet_tolerance", label=f"thing {i}",
            input_state="not_supplied", limitations=("would sharpen the fibre ranking",),
        ))
    grouped = register.missing_categories()
    assert len(grouped) == 1
    assert grouped[0]["n_not_supplied"] == 6
    assert len(grouped[0]["examples"]) == 4, "a concise category, not every row"
    assert grouped[0]["what_it_would_add"]


def test_influential_assumptions_exclude_the_ones_nothing_reads():
    register = IN.InputRegister(sample_id="S1")
    register.add(IN.InputRecord(
        input_id="a", category="diet_tolerance", label="used",
        input_state="assumed_default", effective_value="x", default_rationale="r",
        affects=(IN.Influence("t", "r", "changes_calculation", "e"),),
    ))
    register.add(IN.InputRecord(
        input_id="b", category="diet_tolerance", label="unused",
        input_state="assumed_default", effective_value="y", default_rationale="r",
    ))
    assert len(register.assumed) == 2
    assert [r.input_id for r in register.influential_assumptions] == ["a"]


def test_the_contract_is_stated_in_the_output():
    text = " ".join(IN.InputRegister(sample_id="S1").to_json()["contract"]).lower()
    assert "never becomes a supplied fact" in text
    assert "no imputed normal value" in text
    assert "not used in this report" in text
    assert "neither is silently overwritten" in text


# --------------------------------------------------------------------------- #
# building from a real results object
# --------------------------------------------------------------------------- #


def _results() -> dict:
    return {
        "sample": "S1",
        "subject_context": {
            "mode": "research",
            "supplied": {},
            "assumed": [
                {"field": "metformin", "label": "metformin", "kind": "medication",
                 "assumed_value": False, "used_by": "the diabetes pattern",
                 "if_present": "Metformin raises Escherichia."},
            ],
            "unknown": [{"label": "sex", "if_present": "Sex-specific patterns could be scored."}],
        },
        "extension": {"views": {
            "simulation": {
                "readiness": {"can_run": True},
                "protocol": {"protocol_id": "p/1.0"},
                "media": [{"id": "european", "label": "Average European diet"}],
                "scenarios": {"baseline_medium": "european"},
            },
            "ecology": {"community_type_model": {
                "model_id": "ct/1.0", "cohort": {"n_reference_samples": 3027, "n_studies": 22},
            }},
        }},
    }


def test_build_register_records_the_assumed_diet_as_an_assumption():
    register = IN.build_register(_results(), sample_id="S1")
    diet = next(r for r in register.records if r.input_id == "diet.scenario_medium")
    assert diet.input_state == "assumed_default"
    assert diet.reported_value is None, "the person said nothing about their diet"
    assert diet.effective_value == "Average European diet"
    assert "not your actual diet" in " ".join(diet.limitations)
    assert diet.affects[0].kind == "changes_calculation"


def test_build_register_imports_the_existing_context_ledger_unchanged():
    register = IN.build_register(_results(), sample_id="S1")
    metformin = next(r for r in register.records if r.input_id == "context.metformin")
    assert metformin.input_state == "inherited_default"
    assert metformin.effective_value is False
    assert "Metformin raises Escherichia." in metformin.affects[0].explanation
    unknown = next(r for r in register.records if r.label == "sex")
    assert unknown.input_state == "not_supplied"


def test_build_register_names_the_frozen_cohort_as_an_assumption():
    register = IN.build_register(_results(), sample_id="S1")
    cohort = next(r for r in register.records if r.input_id == "reference.community_type_cohort")
    assert cohort.input_state == "assumed_default"
    assert "3027 samples from 22 studies" in (cohort.default_rationale or "")
    assert "never against the other samples" in cohort.affects[0].explanation


def test_build_register_says_no_laboratory_results_rather_than_omitting_the_topic():
    register = IN.build_register(_results(), sample_id="S1")
    lab = next(r for r in register.records if r.category == "laboratory")
    assert lab.input_state == "not_supplied"
    assert "none are needed" in " ".join(lab.limitations).lower()


def test_lab_results_become_records_that_do_not_overrule_the_genetic_reading():
    records = IN.records_from_lab_results([_lab()])
    assert len(records) == 1
    assert records[0].input_state == "supplied"
    assert records[0].category == "laboratory"
    text = " ".join(records[0].limitations)
    assert "does not replace, confirm or invalidate" in text


def test_the_whole_register_serialises_and_round_trips():
    register = IN.build_register(_results(), sample_id="S1", lab_results=[_lab()])
    payload = register.to_json()
    assert json.loads(json.dumps(payload)) == payload
    assert payload["feature_id"] == "A14"
    assert payload["summary"]["n_laboratory_results"] == 1
    categories = {g["category"] for g in payload["by_category"]}
    assert {"diet_tolerance", "medications_supplements", "analysis_configuration"} <= categories
