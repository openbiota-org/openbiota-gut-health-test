"""A16 model-assisted scenarios — BUILD_SPEC_v0.8.3 sections 11.6/11.7.

AT101-AT112. The engine answers "is the combination worth more than either
part?", and almost every test here is about the ways that answer could be
overstated: an invalid recipe quietly rescaled, an infeasible solve read as
zero production, two objectives with different winners reported in the same
units, or a published table presented as a measurement in a person.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.extension import simulation as SIM
from openbiota.extension.strict_ooxml import StrictOoxmlError, StrictWorkbook

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def config() -> SIM.SimulationConfig:
    return SIM.SimulationConfig.load()


# --------------------------------------------------------------------------- #
# AT101 — the mass check happens before anything is normalised
# --------------------------------------------------------------------------- #


def test_at101_the_published_notebook_recipe_fails_the_mass_fixture(config) -> None:
    """.95 plus five .04 additions totals 1.15 and must not enter production.

    Normalising first would rescale every resident abundance by 1/1.15 and
    produce a plausible-looking run from an invalid recipe, which is exactly
    the repair the specification forbids.
    """
    baseline = {"Bacteroides fragilis": 60.0, "Faecalibacterium prausnitzii": 40.0}
    with pytest.raises(SIM.SimulationError, match="does not equal 1.0"):
        SIM.build_perturbation(
            baseline, protocol=config.protocol,
            added=config.added_species, resident_fraction=0.95,
        )
    # The corrected protocol totals exactly one.
    ok = SIM.build_perturbation(baseline, protocol=config.protocol, added=config.added_species)
    assert ok.total_before_normalisation == pytest.approx(1.0, abs=1e-12)
    assert sum(ok.abundances.values()) == pytest.approx(1.0, abs=1e-12)
    assert len(ok.abundances) == len(baseline) + len(config.added_species)
    assert not ok.residual_normalised


def test_at101_the_error_names_the_arithmetic_rather_than_just_failing(config) -> None:
    baseline = {"a": 1.0}
    with pytest.raises(SIM.SimulationError) as excinfo:
        SIM.build_perturbation(
            baseline, protocol=config.protocol, added=("x", "y"), resident_fraction=0.95,
        )
    message = str(excinfo.value)
    assert "0.95" in message and "0.04" in message
    assert "dividing its total away" in message


def test_a_negative_or_empty_community_is_refused(config) -> None:
    for bad, match in (
        ({"a": -1.0, "b": 2.0}, "negative baseline"),
        ({"a": 0.0}, "zero total abundance"),
        ({"a": float("inf")}, "non-finite"),
    ):
        with pytest.raises(SIM.SimulationError, match=match):
            SIM.build_perturbation(bad, protocol=config.protocol, added=())


def test_an_added_species_already_resident_is_merged_once(config) -> None:
    """Section 11.6: merge true aliases once, and say that it happened."""
    baseline = {"Bifidobacterium longum": 50.0, "Bacteroides fragilis": 50.0}
    p = SIM.build_perturbation(
        baseline, protocol=config.protocol, added=("Bifidobacterium longum",),
        added_fraction_each=0.2, resident_fraction=0.8,
    )
    assert len(p.abundances) == 2
    assert p.abundances["Bifidobacterium longum"] == pytest.approx(0.4 + 0.2)
    assert sum(p.abundances.values()) == pytest.approx(1.0, abs=1e-12)
    assert p.note and "already resident" in p.note


# --------------------------------------------------------------------------- #
# AT102 — coverage is reported before renormalisation
# --------------------------------------------------------------------------- #


def test_at102_mapping_loss_is_reported_and_nothing_is_reassigned() -> None:
    abundances = {"modelled one": 0.5, "modelled two": 0.3, "no model here": 0.19, "tiny": 0.01}
    mapped, coverage = SIM.map_to_models(
        abundances, ["modelled one", "modelled two"], cutoff=0.05,
    )
    assert set(mapped) == {"modelled one", "modelled two"}
    assert coverage.abundance_unmapped == pytest.approx(0.19)
    assert coverage.abundance_below_cutoff == pytest.approx(0.01)
    assert coverage.mapped_fraction == pytest.approx(0.8)
    # The unmapped organism is named, not silently folded into a neighbour.
    assert "no model here" in coverage.unmapped_examples
    assert "no model here" not in mapped
    payload = coverage.to_json()
    assert payload["denominator"] == "bacterial_abundance_offered_to_the_model"
    assert any("Fungal, viral" in limit for limit in payload["limitations"])


def test_coverage_denominator_is_bacterial_and_says_so() -> None:
    """Section 11.6: do not mix bacterial coverage with other denominators."""
    _mapped, coverage = SIM.map_to_models({"a": 1.0}, ["a"], cutoff=0.0)
    assert "bacterial" in coverage.denominator
    assert any("does not prove complete functional coverage" in x
               for x in coverage.to_json()["limitations"])


# --------------------------------------------------------------------------- #
# AT104 — the four-arm contrasts
# --------------------------------------------------------------------------- #


def test_at104_the_worked_four_arm_fixture(config) -> None:
    """Values 1, 3, 4, 8 give pair 7, over-substrate 4, over-probiotic 5, I 2."""
    fixture = config.contrast_fixture
    c = SIM.contrasts(fixture["f00"], fixture["f10"], fixture["f01"], fixture["f11"])
    expected = fixture["expected"]
    assert c.delta_pair == pytest.approx(expected["delta_pair"])
    assert c.delta_over_substrate == pytest.approx(expected["delta_over_substrate"])
    assert c.delta_over_probiotic == pytest.approx(expected["delta_over_probiotic"])
    assert c.interaction == pytest.approx(expected["interaction"])
    assert c.missing_arms == ()


def test_at104_a_missing_arm_prevents_its_contrasts() -> None:
    c = SIM.contrasts(1.0, None, 4.0, 8.0)
    assert c.delta_over_probiotic is None
    assert c.interaction is None
    # The contrasts that do not need the missing arm still compute.
    assert c.delta_pair == pytest.approx(7.0)
    assert c.delta_over_substrate == pytest.approx(4.0)
    assert c.missing_arms == ("probiotic_only",)
    assert SIM.contrasts(None, None, None, None).delta_pair is None


def test_an_interaction_contrast_is_labelled_a_model_result(config: SIM.SimulationConfig) -> None:
    fixture = config.contrast_fixture
    payload = SIM.contrasts(
        fixture["f00"], fixture["f10"], fixture["f01"], fixture["f11"]
    ).to_json()
    assert "not proof of clinical synergy" in payload["interpretation"]


def test_at094_fibre_alone_is_allowed_to_beat_the_combination() -> None:
    """A negative interaction is a real answer, not something to suppress."""
    c = SIM.contrasts(4.57, 5.64, 25.47, 11.53)
    assert c.delta_over_substrate is not None and c.delta_over_substrate < 0
    assert c.interaction is not None and c.interaction < 0
    assert c.delta_pair is not None and c.delta_pair > 0


# --------------------------------------------------------------------------- #
# AT103 — an unavailable answer is not a zero
# --------------------------------------------------------------------------- #


def test_at103_flux_per_growth_refuses_a_nongrowing_community() -> None:
    assert SIM.flux_per_growth(10.0, 0.05) == pytest.approx(200.0)
    assert SIM.flux_per_growth(10.0, 0.0) is None
    assert SIM.flux_per_growth(10.0, -0.1) is None
    assert SIM.flux_per_growth(None, 0.05) is None
    assert SIM.flux_per_growth(10.0, float("nan")) is None


def test_at103_a_result_with_no_flux_carries_a_reason_not_a_zero() -> None:
    result = SIM.ScenarioResult(
        scenario_id="s", arm="both", execution_state="solved_model",
        solver_status="infeasible", raw_flux=None, community_growth=None,
        flux_per_growth=None, target_exchange="EX_but_m", medium_id="european",
        coverage=None, perturbation=None,
        unavailable_reason="solver status 'infeasible': unavailable, not zero production",
    )
    payload = result.to_json()
    assert payload["raw_flux"] is None
    assert "not zero production" in payload["unavailable_reason"]
    assert payload["raw_flux_unit"] == "mmol/gDW/h"


# --------------------------------------------------------------------------- #
# AT105/AT108 — units, caps and the two objectives
# --------------------------------------------------------------------------- #


def test_at105_exchange_caps_are_total_caps_with_their_units(config) -> None:
    expected = {
        "inulin": ("EX_inulin_m", 6.14), "pectin": ("EX_pect_m", 0.4),
        "starch": ("EX_strch1_m", 16.65), "maltodextrin": ("EX_dextrin_m", 30.28),
        "cellulose": ("EX_cellul_m", 0.37), "arabinoxylan": ("EX_arabinoxyl_m", 4.09),
    }
    for sid, (exchange, cap) in expected.items():
        s = config.substrate(sid)
        assert s.exchange == exchange
        assert s.cap == pytest.approx(cap)
        payload = s.to_json()
        assert payload["cap_semantics"].startswith("total cap")
        assert payload["cap_mmol_per_gDW_per_h"] == pytest.approx(cap)


def test_at105_a_model_proxy_says_it_is_not_the_food(config) -> None:
    """Cellulose is not hemp seed; arabinoxylan is not psyllium."""
    for sid in ("cellulose", "arabinoxylan"):
        s = config.substrate(sid)
        assert s.is_proxy
        assert s.proxy_for and "not" in s.proxy_for.lower()
    for sid in ("inulin", "pectin", "starch", "maltodextrin"):
        assert not config.substrate(sid).is_proxy


def test_at108_the_two_objectives_have_different_units_and_different_winners() -> None:
    """Swapping them would swap the winner, so they can never share a field."""
    assert SIM.TARGET_EXCHANGES["butyrate"] == SIM.BUTYRATE_EXCHANGE
    result = SIM.ScenarioResult(
        scenario_id="s", arm="both", execution_state="solved_model", solver_status="optimal",
        raw_flux=19.6188886985, community_growth=0.0528202764673425,
        flux_per_growth=SIM.flux_per_growth(19.6188886985, 0.0528202764673425),
        target_exchange="EX_but_m", medium_id="european", coverage=None, perturbation=None,
    ).to_json()
    assert result["raw_flux_unit"] == "mmol/gDW/h"
    assert result["flux_per_community_growth_unit"] == "mmol/gDW/h per 1/h"
    assert result["raw_flux"] != result["flux_per_community_growth"]


def test_metabolites_are_never_summed_into_one_scfa_number() -> None:
    result = SIM.ScenarioResult(
        scenario_id="s", arm="both", execution_state="solved_model", solver_status="optimal",
        raw_flux=11.5, community_growth=0.029, flux_per_growth=396.0,
        target_exchange="EX_but_m", medium_id="european", coverage=None, perturbation=None,
    )
    result.other_exchanges = {"propionate": 227.2, "acetate": 228.0, "lactate": 0.0}
    payload = result.to_json()
    assert set(payload["other_exchanges"]) == {"propionate", "acetate", "lactate"}
    assert "never summed" in payload["other_exchanges_note"]
    assert "total_scfa" not in payload


# --------------------------------------------------------------------------- #
# AT106/AT107/AT108/AT109 — the published replay
# --------------------------------------------------------------------------- #


def _fixture_installed(config) -> bool:
    return (REPO / str(config.replay.get("file", ""))).is_file()


def test_at106_a_parser_returning_zero_sheets_fails_loudly(tmp_path: Path) -> None:
    """The exact Strict OOXML failure the specification calls out.

    openpyxl opens this workbook and reports no sheets at all. An empty
    dataset silently passes every assertion about values it no longer holds,
    so the reader raises instead.
    """
    import zipfile

    broken = tmp_path / "no_sheets.xlsx"
    with zipfile.ZipFile(broken, "w") as z:
        z.writestr("xl/workbook.xml",
                   '<workbook xmlns="http://purl.oclc.org/ooxml/spreadsheetml/main"><sheets/></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", "<Relationships/>")
    with pytest.raises(StrictOoxmlError, match="no readable sheets"):
        StrictWorkbook.open(broken)

    not_a_book = tmp_path / "empty.xlsx"
    with zipfile.ZipFile(not_a_book, "w") as z:
        z.writestr("hello.txt", "nothing")
    with pytest.raises(StrictOoxmlError, match="not an xlsx package"):
        StrictWorkbook.open(not_a_book)


def test_at106_the_real_strict_workbook_reads(config) -> None:
    if not _fixture_installed(config):
        pytest.skip("replay fixture not installed; run `make extension-refs`")
    with StrictWorkbook.open(REPO / str(config.replay["file"])) as wb:
        assert wb.strict, "this workbook is the Strict OOXML variant"
        assert "production_but" in wb.sheet_names
        assert "production_but_normalized" in wb.sheet_names
        records = wb.records("production_but")
        assert len(records) > 1000
        assert {"sample_id", "flux", "treatment", "probiotic"} <= set(records[0])
        with pytest.raises(StrictOoxmlError, match="no sheet named"):
            wb.records("not_a_sheet")


def test_at107_at108_the_published_values_and_both_orderings(config) -> None:
    """Every released value to 1e-8, and the two objectives' different winners."""
    if not _fixture_installed(config):
        pytest.skip("replay fixture not installed; run `make extension-refs`")
    replay = SIM.replay_published_fixture(config)
    assert replay.passed, [c for c in replay.checks if not c["within_tolerance"]]
    assert replay.tolerance == pytest.approx(1e-8)
    assert len(replay.checks) >= 12

    # AT107: fibre alone wins on raw flux.
    assert replay.ranking_raw[0] == "Psyllium Husk"
    # AT108: maltodextrin wins per unit of community growth.
    assert replay.ranking_normalized[0] == "Maltodextrin"
    assert replay.ranking_raw[0] != replay.ranking_normalized[0]

    # The excluded non-standard diet is absent from the candidate set.
    assert all("High Fiber" not in c for c in replay.ranking_raw)


def test_at109_a_replay_is_not_a_solved_model_and_not_a_measurement(config) -> None:
    if not _fixture_installed(config):
        pytest.skip("replay fixture not installed")
    replay = SIM.replay_published_fixture(config)
    assert replay.execution_state == "published_replay"
    assert replay.execution_state in SIM.EXECUTION_STATES
    assert "solved_model" in SIM.EXECUTION_STATES
    provenance = replay.to_json()["provenance"]
    assert "not simulated here" in provenance
    assert "not an outcome measured in" in provenance


def test_at110_a_missing_sample_is_audited_not_invented(config) -> None:
    """The 154-input/156-output discrepancy must surface, not be patched."""
    if not _fixture_installed(config):
        pytest.skip("replay fixture not installed")
    import copy

    altered = copy.deepcopy(config)
    altered.replay = dict(config.replay) | {"sample_id": "00000000000000"}
    with pytest.raises(SIM.SimulationError, match="intersection"):
        SIM.replay_published_fixture(altered)


# --------------------------------------------------------------------------- #
# AT111 — sensitivity is a range, not an interval
# --------------------------------------------------------------------------- #


def test_at111_sensitivity_reports_a_range_and_says_it_is_not_a_confidence_interval() -> None:
    s = SIM.sensitivity({"european": 10.0, "high_fibre": 25.0, "failed": None},
                        rank_changed=True)
    assert s.minimum == pytest.approx(10.0)
    assert s.maximum == pytest.approx(25.0)
    assert s.median == pytest.approx(17.5)
    assert set(s.scenarios) == {"european", "high_fibre"}
    payload = s.to_json()
    assert payload["rank_changed_across_scenarios"] is True
    assert "not a confidence interval" in payload["interpretation"]
    # Nothing solved: a range over nothing is nothing.
    empty = SIM.sensitivity({"european": None})
    assert empty.minimum is None and empty.median is None


# --------------------------------------------------------------------------- #
# the pinned protocol and its recorded discrepancies
# --------------------------------------------------------------------------- #


def test_the_protocol_is_the_corrected_one_and_records_the_repairs(config) -> None:
    p = config.protocol
    assert p.micom_version == "0.37.0"
    assert p.cutoff == pytest.approx(0.001)
    assert p.tradeoff == pytest.approx(0.99)
    assert p.strategy == "none"
    assert p.resident_fraction == pytest.approx(0.8)
    assert p.added_fraction_per_species == pytest.approx(0.04)
    assert p.resident_fraction + 5 * p.added_fraction_per_species == pytest.approx(1.0)
    assert len(config.added_species) == 5

    payload = p.to_json()
    discrepancy = payload["growth_threshold_discrepancy"]
    assert discrepancy["methods_code"] == pytest.approx(0.01)
    assert discrepancy["figure_caption"] == pytest.approx(0.001)
    assert "Preserved, not repaired" in discrepancy["resolution"]


def test_the_emitted_metrics_are_labelled_predictions(config) -> None:
    arms = (
        SIM.ScenarioResult(
            scenario_id="c.both", arm="both", execution_state="solved_model",
            solver_status="optimal", raw_flux=11.53, community_growth=0.029,
            flux_per_growth=396.0, target_exchange="EX_but_m", medium_id="european",
            coverage=SIM.Coverage(157, 55, 73.1, 27.0, 1.7, "bacterial", (), 0.001),
            perturbation=None,
        ),
    )
    contrast = SIM.contrasts(4.57, 5.64, 25.47, 11.53)
    metrics = SIM.scenario_metrics(
        arms, contrast, candidate_id="inulin_consortium",
        substrate=config.substrate("inulin"), probiotic_label="5-species consortium",
        config=config,
    )
    assert metrics
    for m in metrics:
        assert m.kind == "model_prediction"
        assert m.feature_id == "A16"
        assert m.evidence_maturity == "genomic_prediction"
        assert m.reference_percentile is None
        joined = " ".join(m.limitations)
        assert "not a measurement" in joined
        assert "not proof of clinical synergy" in joined
        m.to_json()


def test_readiness_separates_the_three_components(config) -> None:
    """A feature flag controls whether simulations run, not whether they exist."""
    state = SIM.readiness(config)
    assert set(state) >= {"can_run", "micom", "solver", "model_library"}
    assert set(state["micom"]) >= {"installed", "pinned"}
    assert state["micom"]["pinned"] == "0.37.0"
    # Whatever is installed here, the three are reported separately.
    assert isinstance(state["solver"]["installed"], bool)
    assert isinstance(state["model_library"]["installed"], bool)


def test_the_model_library_checksum_is_pinned(config) -> None:
    assert config.protocol.model_asset_md5 == "6a13da2a3fd9b1bc059987d6344f32b6"
    assert config.protocol.model_asset == "agora201_refseq216_species_1.qza"
    assert config.protocol.model_source_id == "S11"


@pytest.mark.parametrize("medium_id", ["european", "high_fibre"])
def test_both_declared_media_are_pinned_and_present(config, medium_id: str) -> None:
    entry = next(m for m in config.media if m["id"] == medium_id)
    assert entry["sha256"]
    path = REPO / str(entry["file"])
    if not path.is_file():
        pytest.skip("media not installed; run `make extension-refs`")
    import hashlib

    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]


def test_a_cached_scenario_is_returned_unchanged() -> None:
    """A regenerated report shows the value that was solved, not a new one."""
    cache = SIM.CACHE_DIR / "solutions"
    if not cache.is_dir() or not any(cache.glob("*.json")):
        pytest.skip("no solved scenarios cached in this checkout")
    for path in sorted(cache.glob("*.json"))[:4]:
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["execution_state"] in SIM.EXECUTION_STATES
        assert "cache_key" in payload
        if payload["raw_flux"] is not None:
            assert payload["solver_status"] in SIM.OPTIMAL_STATUSES
