"""Essential checks for the mycobiome module (BUILD_SPEC_v08.2 §12.3).

Numbered MYC-xx where a test enforces one of the spec's listed checks.
Tests that need the reference lock or a run are skipped when neither exists;
the arithmetic, registry and rule tests always run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.mycobiome import registry as R
from openbiota.mycobiome import score as S
from openbiota.mycobiome.genome_lane import poisson_interval

RESULTS = sorted(Path("results").glob("*/results.json"))


def _c(rid, d, cap, act, alts=None, status="active", claim=None, fids=()):
    return S.Contribution(rule_id=rid, direction=d, biological_claim_id=claim or rid, source_ids=("s",),
                          required_identity="species", trigger_type="presence", base_cap=cap, activation=act,
                          alternatives=tuple(alts or (act,)), activation_status=status, finding_ids=tuple(fids))


# ---- MYC-77: worked examples ----------------------------------------------- #

@pytest.mark.parametrize("b,c,expected", [(0.5, 0, 75), (0, 0.5, 25), (0.5, 0.5, 37.5), (0.25, 0, 62.5), (0, 0.25, 37.5)])
def test_myc77_worked_examples(b, c, expected):
    assert 100 * S.score_01(b, c) == pytest.approx(expected)


# ---- MYC-75: bounded, monotone, deterministic ----------------------------- #

@pytest.mark.parametrize("b", [0.0, 0.125, 0.25, 0.5, 0.75, 1.0])
@pytest.mark.parametrize("c", [0.0, 0.125, 0.25, 0.5, 0.75, 1.0])
def test_myc75_bounds_and_monotonicity(b, c):
    v = 100 * S.score_01(b, c)
    assert 0.0 <= v <= 100.0
    if b < 1.0:
        assert 100 * S.score_01(min(1.0, b + 0.1), c) >= v
    if c < 1.0:
        assert 100 * S.score_01(b, min(1.0, c + 0.1)) <= v
    assert S.compute([_c("x", "B", b, 1.0)] if b else []).score_100 == S.compute([_c("x", "B", b, 1.0)] if b else []).score_100


# ---- MYC-76 / MYC-28: null with a reason, never 0/50/100 -------------------- #

def test_completed_search_with_nothing_to_weigh_is_neutral_50():
    """Product decision: a finished analysis that found no concerning and no
    beneficial fungal evidence reads as neutral (50), not as a failure."""
    s = S.compute([])
    assert s.score_100 == 50.0 and s.status == "computed" and s.direction == "neutral"
    assert "no_directional_evidence" in s.reason_codes and s.mixed_signals is False


def test_myc76_failures_are_null_with_reason_not_50():
    s2 = S.compute([_c("x", "B", 0.5, 0.0, status="inactive")], reason_if_empty="insufficient_fungal_information")
    assert s2.score_100 is None and "insufficient_fungal_information" in s2.reason_codes
    assert s2.direction == "not_assessed"
    s3 = S.compute([], reason_if_empty="incomplete_analysis", incomplete=True)
    assert s3.score_100 is None and s3.status == "not_computable"
    s4 = S.compute([], reason_if_empty="failed_qc")
    assert s4.score_100 is None


# ---- MYC-80: a dominant concern cannot be erased by benefit ---------------- #

def test_myc80_concern_caps_the_score():
    s = S.compute([_c("b", "B", 1.0, 1.0), _c("c", "C", 0.5, 1.0)])
    assert s.score_100 <= 50.0
    assert s.mixed_signals is True


# ---- MYC-78 / MYC-79: duplicates and unknowns do not move the score -------- #

def test_myc78_duplicate_claims_do_not_add():
    one = S.compute([_c("a", "B", 0.25, 1.0, claim="fam")]).score_100
    two = S.compute([_c("a", "B", 0.25, 1.0, claim="fam"), _c("a2", "B", 0.25, 1.0, claim="fam")]).score_100
    assert one == two == pytest.approx(62.5)


def test_myc79_unscored_findings_leave_score_unchanged_but_visible():
    base = S.compute([_c("a", "B", 0.25, 1.0)])
    more = S.compute([_c("a", "B", 0.25, 1.0)], unscored_findings=["Unknown fungus 1", "Unknown fungus 2"])
    assert base.score_100 == more.score_100 and len(more.unscored_findings) == 2


# ---- MYC-87 / MYC-88: sensitivity envelopes ------------------------------- #

def test_myc87_policy_and_combined_ranges():
    s = S.compute([_c("P4013B_B", "B", 0.5, 1.0)])
    assert (s.policy_range.low, s.policy_range.high) == (pytest.approx(62.5), pytest.approx(87.5))
    s2 = S.compute([_c("P4013B_B", "B", 0.5, 1.0), _c("CA_HD_C", "C", 0.5, 0.0, alts=(0.0, 1.0), status="unresolved")])
    assert s2.score_100 == pytest.approx(75.0)
    assert (s2.resolution_range.low, s2.resolution_range.high) == (pytest.approx(37.5), pytest.approx(75.0))
    assert (s2.display_range.low, s2.display_range.high) == (pytest.approx(15.625), pytest.approx(87.5))
    assert s2.display_range.low <= s2.score_100 <= s2.display_range.high


def test_myc88_range_crossing_midpoint_is_uncertain():
    s = S.compute([_c("P4013B_B", "B", 0.5, 1.0), _c("CA_HD_C", "C", 0.5, 0.0, alts=(0.0, 1.0), status="unresolved")])
    assert s.direction == "uncertain" and s.status == "provisional"
    assert S.compute([_c("b", "B", 0.5, 1.0)]).direction == "favorable_evidence"
    assert S.compute([_c("c", "C", 0.5, 1.0)]).direction == "concerning_evidence"


def test_myc93_interpretation_confidence_is_always_experimental():
    for s in (S.compute([]), S.compute([_c("b", "B", 0.5, 1.0)], analytical_confidence="supported")):
        assert s.interpretation_confidence == "experimental_limited"


# ---- MYC-86: abundance activation ramp ------------------------------------ #

def test_myc86_abundance_activation():
    assert S.abundance_activation(50.0) == 0.0
    assert S.abundance_activation(90.0) == pytest.approx(0.8)
    assert S.abundance_activation(None) is None
    assert S.abundance_activation(10.0) == 0.0


# ---- registry and rules ---------------------------------------------------- #

def test_myc12_synonyms_are_one_entry():
    reg = R.load()
    assert reg.lookup("Candida glabrata") is reg.lookup("Nakaseomyces glabratus")
    assert reg.lookup("Candida_krusei").name == "Pichia kudriavzevii"
    assert reg.lookup("Aspergillus fumigatus").name == "Aspergillus"   # genus-level entry
    assert reg.lookup("Not a fungus at all") is None


def test_myc26_evidence_labels_are_distinct_and_from_the_allowed_set():
    reg = R.load()
    allowed = set(R.EVIDENCE_LABELS)
    labels = {c.label for t in reg.taxa for c in t.evidence}
    assert labels <= allowed
    assert {"human_trial", "animal_experiment", "in_vitro"} <= labels
    assert len({R.EVIDENCE_LABELS[k] for k in labels}) == len(labels)


def test_myc43_every_claim_states_its_subject_and_limit():
    for t in R.load().taxa:
        for c in t.evidence:
            assert c.source.startswith("http"), (t.name, c.claim_id)
            assert c.endpoint and c.host and c.direction in ("favorable", "concerning", "context"), (t.name, c.claim_id)
            assert c.limit, (t.name, c.claim_id)


def _sup(name, state="supported", res="unresolved", **kw):
    return R.SupportedTaxon(finding_id=f"f:{name}", name=name, rank=kw.pop("rank", "species"), detection_state=state,
                            strain_resolution=res, **kw)


def test_myc83_p4013b_species_versus_genotype():
    reg = R.load()
    by = {c.rule_id: c for c in R.evaluate(reg, [_sup("Clavispora lusitaniae")])}
    assert by["CL_LUSITANIAE_B"].activation_status == "active" and by["CL_LUSITANIAE_B"].value == 0.25
    assert by["P4013B_B"].activation_status == "inactive"          # not even compatible without a strain attempt
    by2 = {c.rule_id: c for c in R.evaluate(reg, [_sup("Clavispora lusitaniae", res="reference_genotype",
                                                        reference_accessions=("GCA_023627835.1",))])}
    assert by2["P4013B_B"].activation_status == "active" and by2["P4013B_B"].value == 0.5
    by3 = {c.rule_id: c for c in R.evaluate(reg, [_sup("Clavispora lusitaniae", compatible_accessions=("GCA_023627835.1", "X"),
                                                        strain_analysis_status="complete")])}
    assert by3["P4013B_B"].activation_status == "unresolved" and by3["P4013B_B"].alternatives == (0.0, 1.0)


def test_myc81_myc23_generic_albicans_does_not_become_hd():
    reg = R.load()
    by = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida albicans", res="species_only")])}
    assert by["CA_HD_C"].activation_status == "inactive"
    hd = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida albicans", res="reference_genotype", reference_genotypes=("IDB311",))])}
    assert hd["CA_HD_C"].activation_status == "active" and hd["CA_HD_C"].value == 0.5
    ld = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida albicans", res="reference_genotype", reference_genotypes=("IDC561",))])}
    assert ld["CA_LD_CONTEXT"].value == 0.0 and ld["CA_HD_C"].activation_status == "inactive"


def test_myc82_hd_ld_equivalence_widens_the_range():
    reg = R.load()
    contribs = R.evaluate(reg, [_sup("Candida albicans", res="reference_equivalence_group",
                                     compatible_genotypes=("IDB311", "IDC561"), strain_analysis_status="complete"),
                                _sup("Clavispora lusitaniae")])
    s = S.compute(contribs)
    assert s.score_100 == pytest.approx(62.5) and s.direction == "uncertain"
    assert s.resolution_range.low == pytest.approx(31.25)


def test_myc84_myc22_boulardii_rule_is_unbound():
    reg = R.load()
    by = {c.rule_id: c for c in R.evaluate(reg, [_sup("Saccharomyces cerevisiae", res="lineage")])}
    assert by["SB_STUDIED_B"].activation_status == "inactive" and "unbound" in by["SB_STUDIED_B"].reason


def test_myc85_presence_of_food_yeast_awards_nothing():
    reg = R.load()
    s = S.compute(R.evaluate(reg, [_sup("Saccharomyces cerevisiae"), _sup("Pichia kudriavzevii"), _sup("Agaricus bisporus")]))
    assert s.score_100 == 50.0 and s.direction == "neutral" and s.benefit_contribution == 0.0


def test_co_presence_rule_needs_all_bacteria():
    reg = R.load()
    only = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida tropicalis")], bacteria_present=["Escherichia coli"])}
    assert only["CT_EC_SM_INTERACTION_C"].activation_status == "inactive"
    both = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida tropicalis")], bacteria_present=["Escherichia coli", "Serratia marcescens"])}
    assert both["CT_EC_SM_INTERACTION_C"].value == 0.125


def test_myc13_complex_only_finding_does_not_activate_species_rule():
    reg = R.load()
    by = {c.rule_id: c for c in R.evaluate(reg, [_sup("Candida parapsilosis", rank="complex")])}
    assert by["CPAR_GUT_C"].activation_status == "inactive"


# ---- denominators ----------------------------------------------------------- #

def test_myc06_fraction_and_per_million_agree():
    from openbiota.mycobiome.genome_lane import Ledger
    from openbiota.mycobiome.quantify import measurement
    m = measurement(Ledger(eligible_fragments=8_029_909, aligned_fragments=3_000_000, fungal_confident=244),
                    eligible=8_029_909, denominator_valid=True, lane_status="resolved")
    assert m.fungal_fragments_per_million == pytest.approx(1e6 * m.fungal_fragment_fraction)
    assert 0 <= m.fungal_fragment_fraction <= 1 and m.unassigned_fragments == 5_029_909


def test_myc20_myc02_no_denominator_is_not_zero_fungi():
    from openbiota.mycobiome.genome_lane import Ledger
    from openbiota.mycobiome.quantify import measurement
    m = measurement(Ledger(eligible_fragments=0, fungal_confident=12), eligible=None, denominator_valid=False, lane_status="resolved")
    assert m.fungal_fragment_fraction is None and m.quantification_status == "not_quantifiable_from_available_input"
    assert m.fungal_supported_fragments == 12
    n = measurement(None, eligible=100, denominator_valid=True, lane_status="not_assessed")
    assert n.quantification_status == "not_assessed" and n.fungal_supported_fragments is None


def test_myc21_wilson_interval_in_unit_range():
    lo, hi = poisson_interval(244, 8_029_909)
    assert 0 <= lo <= 244 / 8_029_909 <= hi <= 1


# ---- composition: one estimator, unresolved kept out ----------------------- #

def test_myc47_myc03_composition_excludes_unresolved_and_is_not_sample_percent():
    from openbiota.mycobiome.genome_lane import Ledger, TaxonSupport
    from openbiota.mycobiome.quantify import composition
    a = TaxonSupport(key="sp:1", rank="species", name="A", species_taxid=1, genus="G", accessions=["x"], fragments=200,
                     genome_bases=12_000_000, detection_state="supported")
    b = TaxonSupport(key="sp:2", rank="species", name="B", species_taxid=2, genus="G", accessions=["y"], fragments=100,
                     genome_bases=30_000_000, detection_state="supported")
    g = TaxonSupport(key="genus:G", rank="genus", name="G", species_taxid=None, genus="G", accessions=[], fragments=50,
                     detection_state="ambiguous_complex")
    comp = composition([a, b, g], fragment_bases=300.0, ledger=Ledger(eligible_fragments=1_000_000, fungal_confident=350))
    shares = {r.key: r.share for r in comp.rows}
    assert set(shares) == {"sp:1", "sp:2"} and sum(shares.values()) == pytest.approx(1.0)
    assert shares["sp:1"] > shares["sp:2"]                       # coverage-normalised: 200/12Mb > 100/30Mb
    assert comp.unresolved_fragments == 50 and comp.unresolved_fragment_fraction == pytest.approx(50 / 350)
    assert comp.quantity_type == "genome_coverage_normalised_share"


# ---- MYC-27: never a universal green interval ----------------------------- #

def test_myc27_no_green_interval_text():
    import inspect

    import openbiota.pdfmycobiome as P
    src = inspect.getsource(P)
    assert "0–1%" not in src and "0-1%" not in src and "excellent health" not in src.lower()
    assert "normal" not in inspect.getsource(P.MhsGauge.draw).lower()


# ---- MYC-73: the report renders the gauge and the JSON carries score_01 ----- #

@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_myc73_results_json_score_fields_are_consistent(path: Path):
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or myco.get("analysis_status") in (None, "not_assessed"):
        pytest.skip("mycobiome not run for this sample")
    hs = myco["health_score"]
    if hs["score_100"] is None:
        assert hs["score_01"] is None and hs["status"] in ("not_computable", "not_assessed")
        assert hs["reason_codes"]
        assert myco["analysis_status"] != "complete" or "failed_qc" in hs["reason_codes"]
    else:
        assert hs["score_100"] == pytest.approx(100 * hs["score_01"])
        assert 0 <= hs["score_100"] <= 100
        rng = hs["display_sensitivity_range"]
        assert rng["low"] <= hs["score_100"] <= rng["high"]
    assert hs["interpretation_confidence"] == "experimental_limited"
    assert hs["range_type"] == "scenario_sensitivity_not_confidence_interval"


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_myc07_myc36_each_pair_once_and_mass_conserved(path: Path):
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or not (myco.get("genome_lane") or {}).get("ledger"):
        pytest.skip("no genome lane ledger")
    led = myco["genome_lane"]["ledger"]
    parts = led["fungal_supported_fragments"] + led["fungal_below_identity_fragments"] + led["cross_kingdom_ambiguous_fragments"] + led["decoy_fragments"]
    assert parts <= led["aligned_to_index_fragments"] <= led["eligible_fragments"]
    assert led["unassigned_fragments"] == led["eligible_fragments"] - led["aligned_to_index_fragments"]
    # Fungal taxa never enter the bacterial inventory's composition.
    inv = d.get("organism_inventory") or {}
    fungal_names = {t["accepted_name"] for t in myco.get("taxa") or []}
    assert not any(o["species"].replace("_", " ") in fungal_names for o in inv.get("organisms") or [])


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_pathogen_screen_fungal_signals_all_appear_on_the_mycobiome_page(path: Path):
    """The two modules tell one story: a fungal target the pathogen screen
    reports with any fragments is on the mycobiome page, at the same or a
    more cautious tier, never missing."""
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or myco.get("analysis_status") in (None, "not_assessed") or "pathogen_screen_fungal_signals" not in myco:
        pytest.skip("mycobiome not run with cross-reference")
    listed = {t["accepted_name"] for t in myco["taxa"]} | {t["name"] for t in myco["taxa"]}
    crossed = {(t.get("pathogen_screen") or {}).get("target_id") for t in myco["taxa"] if t.get("pathogen_screen")}
    for sig in myco["pathogen_screen_fungal_signals"]:
        assert sig["target_id"] in crossed or sig["name"] in listed, f"{sig['name']} seen by the pathogen screen but absent here"


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_neutral_score_only_when_the_search_completed(path: Path):
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or myco.get("analysis_status") in (None, "not_assessed"):
        pytest.skip("mycobiome not run")
    hs = myco["health_score"]
    if hs.get("direction") == "neutral":
        assert hs["score_100"] == 50.0 and myco["analysis_status"] == "complete"
        assert not any(c["activation_status"] == "active" and c["value"] > 0 for c in hs["contributions"])


def _ms(**kw):
    base = {"opportunist_share": 0.0, "opportunist_fraction_all": 0.0, "opportunist_supported": 0,
            "alert_supported": 0, "saccharomyces_supported": False, "benefit": 0.0,
            "concern": 0.0, "complete": True}
    return S.myco_score(**(base | kw))["display"]


def test_myco_score_moves_the_way_the_reader_expects():
    """No or very little opportunistic fungi lifts the score; opportunists
    past a tenth of the fungal DNA (the dominance pattern the literature ties
    to disease) pull it down, a confirmed one further; presence at ordinary
    low share is neutral, as the review describes Candida carriage."""
    clean = _ms()
    trace = _ms(opportunist_share=0.0065, opportunist_fraction_all=6e-7)
    ordinary = _ms(opportunist_share=0.05, opportunist_fraction_all=2e-5)
    dominant = _ms(opportunist_share=0.6, opportunist_fraction_all=6e-5)
    dominant_confirmed = _ms(opportunist_share=0.6, opportunist_fraction_all=6e-5, opportunist_supported=1)
    assert clean == trace == 70
    assert ordinary == 50
    assert clean > ordinary > dominant > dominant_confirmed
    assert _ms(saccharomyces_supported=True) == clean + 5
    assert _ms(alert_supported=1) <= S.MYCO_ALERT_CEILING
    assert S.myco_score(opportunist_share=0, opportunist_fraction_all=0, opportunist_supported=0, alert_supported=0,
                            saccharomyces_supported=False, benefit=0, concern=0, complete=False)["value"] is None


def test_myco_score_is_monotone_in_opportunist_share():
    prev = 101
    for share in (0.1, 0.2, 0.4, 0.6, 0.8, 1.0):
        cur = _ms(opportunist_share=share, opportunist_fraction_all=1e-4)
        assert cur <= prev
        prev = cur
    assert prev == 20


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_opportunist_fraction_is_on_the_same_scale_as_fungal_dna(path: Path):
    """The page shows opportunistic fungi as a fraction of all read pairs,
    the same denominator as the fungal-DNA figure, never only as a share of
    the (tiny) fungal total."""
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or myco.get("analysis_status") in (None, "not_assessed") or "colonisation_resistance" not in myco:
        pytest.skip("mycobiome not run with the revised context")
    ctx, meas = myco["context"], myco["measurement"]
    frac_all = ctx["opportunist_fraction_of_all_fragments"]
    if meas.get("eligible_fragments"):
        assert frac_all == pytest.approx(ctx["opportunist_fragments"] / meas["eligible_fragments"])
        assert frac_all <= (meas.get("fungal_fragment_fraction") or 0) + 1e-12
    assert myco["myco_score"]["model_id"] == S.MYCO_SCORE_ID
    assert myco["myco_score"]["inputs"]["opportunist_fraction_of_all_fragments"] == frac_all
    names = {b["bacterium"] for b in myco["colonisation_resistance"]["bacteria"]}
    assert {"Bacteroides thetaiotaomicron", "Blautia producta", "Lacticaseibacillus rhamnosus"} <= names
    inv = {o["species"].replace("_", " ") for o in (d.get("organism_inventory") or {}).get("organisms") or []}
    for b in myco["colonisation_resistance"]["bacteria"]:
        if b["bacterium"] in inv:
            assert b["present"] and b["percent"] is not None


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_myc50_every_supported_fungus_has_a_strain_attempt(path: Path):
    d = json.loads(path.read_text())
    myco = d.get("mycobiome")
    if not myco or myco.get("analysis_status") in (None, "not_assessed"):
        pytest.skip("mycobiome not run")
    for t in myco.get("taxa") or []:
        if t["detection_state"] == "supported" and t["rank"] == "species" and t.get("genome_lane"):
            assert t["strain"]["analysis_status"] in ("complete", "partial", "failed"), t["name"]
            assert t["strain"]["resolution"] != "not_assessed", t["name"]
            if t["strain"]["resolution"] == "unresolved":
                assert t["strain"]["reason_codes"], t["name"]
