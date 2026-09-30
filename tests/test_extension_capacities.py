"""Section 5.8 capacities, and the beta-glucuronidase class assignment.

The second half of this file is the more important half. A structural class
that disagrees with the solved structures is worse than no class at all,
because it reads as substrate-specific information and is not, so the
withholding rule is pinned here rather than left to judgement.
"""

from __future__ import annotations

import pytest

from openbiota.extension import capacities as CAP
from openbiota.extension import gusclasses as GUS

# --------------------------------------------------------------------------- #
# the capacities themselves
# --------------------------------------------------------------------------- #


def test_every_capacity_states_what_it_does_not_establish():
    for capacity in CAP.CAPACITIES:
        assert capacity.does_not_establish.strip(), capacity.capacity_id
        assert len(capacity.does_not_establish.split()) >= 12, capacity.capacity_id


def test_every_capacity_explains_why_its_direction_is_not_simply_good_or_bad():
    """Polyamines are the reason this is enforced rather than encouraged.

    They are required for the gut lining to renew and are elevated in tumour
    tissue. Any report that picked a direction for them would be inventing
    one.
    """
    for capacity in CAP.CAPACITIES:
        assert capacity.context.strip(), capacity.capacity_id
        assert len(capacity.context.split()) >= 10, capacity.capacity_id


def test_a_capacity_cannot_both_require_and_exclude_the_same_gene():
    with pytest.raises(CAP.CapacityError, match="both required and excluded"):
        CAP.Capacity(
            capacity_id="x", label="x", requirement=CAP.Requirement(all_of=("gshA",)),
            panel="glutathione", establishes="x",
            does_not_establish="x" * 5, context="x" * 5,
            reported_separately={"gshA": "no"},
        )


def test_a_capacity_must_say_what_it_does_not_establish():
    with pytest.raises(CAP.CapacityError, match="does not establish"):
        CAP.Capacity(
            capacity_id="x", label="x", requirement=CAP.Requirement(all_of=("gshA",)),
            panel="glutathione", establishes="x", does_not_establish="", context="x",
        )


def test_glutathione_accepts_the_fused_enzyme_on_its_own():
    """gshF does both steps, and a two-enzyme search would miss every
    organism that carries it."""
    capacity = CAP.BY_ID["capacity.glutathione"]
    searched = frozenset({"gshA", "gshB", "gshF"})
    assert capacity.assess(frozenset({"gshF"}), searched) == "complete"
    assert capacity.assess(frozenset({"gshA", "gshB"}), searched) == "complete"
    assert capacity.assess(frozenset({"gshA"}), searched) == "partial"


def test_spermidine_accepts_the_route_the_bacteroidetes_actually_use():
    """Counting only spermidine synthase reports an absence in most
    Bacteroidetes, which is a bug in the pipeline and not a finding."""
    capacity = CAP.BY_ID["capacity.spermidine"]
    searched = frozenset({"speE", "casDC"})
    assert capacity.assess(frozenset({"casDC"}), searched) == "complete"
    assert capacity.assess(frozenset({"speE"}), searched) == "complete"


def test_urate_does_not_count_uricase_or_salvage_towards_the_anaerobic_route():
    capacity = CAP.BY_ID["capacity.urate_anaerobic"]
    assert "uox" in capacity.reported_separately
    assert "hpt" in capacity.reported_separately
    assert "uox" not in capacity.requirement.genes
    assert "hpt" not in capacity.requirement.genes
    # Carrying only the oxygen-requiring enzyme is not the anaerobic route.
    assert capacity.assess(
        frozenset({"uox", "hpt"}), frozenset({"xdhA", "xdhB", "uox", "hpt"})
    ) == "absent"


def test_a_gene_that_was_never_searched_for_is_not_reported_as_absent():
    capacity = CAP.BY_ID["capacity.ethanol"]
    assert capacity.assess(frozenset(), frozenset()) == "not_assayed"
    assert capacity.assess(frozenset(), frozenset({"adhE"})) == "absent"


def test_no_drug_card_claims_a_dose_or_an_exposure():
    banned = ("increase your dose", "reduce your dose", "stop taking", "mg ")
    for card in CAP.DRUG_REACTIONS:
        text = f"{card.mechanism} {card.does_not_mean}".lower()
        for phrase in banned:
            assert phrase not in text, f"{card.reaction_id}: {phrase!r}"
        assert card.does_not_mean.strip(), card.reaction_id
        assert card.source.strip(), card.reaction_id


def test_digoxin_card_records_that_arginine_inhibits_rather_than_increases():
    """The specification calls this out because the direction is commonly
    reported backwards."""
    card = next(c for c in CAP.DRUG_REACTIONS if c.reaction_id == "drug.digoxin")
    assert "inhibited by arginine" in card.mechanism
    assert "increased by arginine" not in card.mechanism


def test_an_unmeasured_drug_reaction_says_so_rather_than_implying_a_reading():
    for card in CAP.DRUG_REACTIONS:
        if card.gene is None:
            assert card.measured_by is None, card.reaction_id
            assert "not measure" in card.does_not_mean.lower(), card.reaction_id


def test_summarise_keeps_unsearched_capacities_out_of_the_assayed_count():
    view = CAP.summarise([], searched=[])
    assert view["n_assayed"] == 0
    assert view["n_complete"] == 0
    assert all(c["state"] == "not_assayed" for c in view["capacities"])


# --------------------------------------------------------------------------- #
# beta-glucuronidase structural classes
# --------------------------------------------------------------------------- #


def _alignment(query: str, subject: str) -> tuple[str, str]:
    assert len(query) == len(subject)
    return query, subject


def test_occupancy_counts_residues_sitting_in_the_anchor_window():
    # Anchor positions 1..10, window 3..6, query fully aligned: four residues.
    q, s = _alignment("ABCDEFGHIJ", "abcdefghij")
    assert GUS.occupancy_length(q, s, 1, (3, 6)) == 4


def test_occupancy_is_near_zero_for_a_query_that_lacks_the_loop():
    q, s = _alignment("AB----GHIJ", "abcdefghij")
    assert GUS.occupancy_length(q, s, 1, (3, 6)) == 0


def test_insertion_counts_only_what_the_query_adds():
    """The second loop is an insertion relative to the anchor, so counting
    occupancy there would score every reference that merely aligns."""
    q, s = _alignment("ABCDEFGHIJ", "abcdefghij")
    assert GUS.insertion_length(q, s, 1, (4, 5), slack=0) == 0
    q, s = _alignment("ABCDXXXEFGHIJ", "abcd---efghij")
    assert GUS.insertion_length(q, s, 1, (4, 5), slack=0) == 3


def test_an_alignment_that_stops_short_is_no_coverage_not_a_missing_loop():
    q, s = _alignment("ABCD", "abcd")
    assert GUS.occupancy_length(q, s, 1, (3, 20)) is None
    assert GUS.classify(None, 0) == "NC"
    assert GUS.classify(0, None) == "NC"


@pytest.mark.parametrize(
    ("loop1", "loop2", "expected"),
    [
        (25, 0, "L1"),    # E. coli, the type example
        (16, 0, "L1"),
        (15, 0, "mL1"),
        (12, 0, "mL1"),
        (9, 0, "NL"),
        (0, 19, "L2"),    # Bacteroides uniformis, the type example
        (0, 12, "L2"),
        (0, 10, "mL2"),
        (0, 0, "NL"),
        (12, 10, "mL1,2"),
    ],
)
def test_classification_follows_the_published_thresholds(loop1, loop2, expected):
    assert GUS.classify(loop1, loop2) == expected


def test_validation_requires_every_published_type_example_to_agree():
    organisms = {f"A{i}": needle for i, (needle, _) in enumerate(GUS.TYPE_EXAMPLES)}
    perfect = {
        f"A{i}": cls for i, (_, cls) in enumerate(GUS.TYPE_EXAMPLES)
    }
    assert GUS.validate({"assignment": perfect}, organisms)["passed"]

    wrong = dict(perfect)
    wrong["A3"] = "L1"  # Bacteroides fragilis called Loop 1; it is mini-Loop 1
    checked = GUS.validate({"assignment": wrong}, organisms)
    assert not checked["passed"]
    assert checked["n_agree"] == len(GUS.TYPE_EXAMPLES) - 1


def test_a_species_with_several_enzymes_agrees_if_any_carries_the_class():
    """Faecalibacterium prausnitzii encodes both a Loop 1 enzyme and the
    FMN-binding No Loop one, so the test has to be generous in this exact
    way or it would fail on correct data."""
    organisms = {"A": "Bacteroides fragilis", "B": "Bacteroides fragilis"}
    checked = GUS.validate({"assignment": {"A": "L1", "B": "mL1"}}, organisms)
    fragilis = next(c for c in checked["checks"] if c["organism"] == "Bacteroides fragilis")
    assert fragilis["agrees"]


def test_the_breakdown_is_withheld_when_validation_fails():
    results = {"panels": [{
        "name": "bglucuronidase",
        "genes": [{"entry_id": "GUS", "accession_fragments": [["P05804", 10]]}],
    }]}
    table = GUS.substrate_table(
        results, {"assignment": {"P05804": "L1"}},
        validation={"passed": False, "checks": []},
    )
    assert table is not None
    assert table["available"] is False
    assert table["classes"] == []
    assert "mini-Loop 1" in table["reason"]
    assert table["what_would_unlock_it"]


def test_the_breakdown_is_shown_when_validation_passes():
    results = {"panels": [{
        "name": "bglucuronidase",
        "genes": [{"entry_id": "GUS",
                   "accession_fragments": [["P05804", 10], ["Q1", 30]]}],
    }]}
    table = GUS.substrate_table(
        results,
        {"assignment": {"P05804": "L1", "Q1": "NL"},
         "reference_counts": {"L1": 1, "NL": 1}, "n_references": 2},
        validation={"passed": True, "checks": []},
    )
    assert table["available"] is True
    assert table["total_fragments"] == 40
    assert table["dominant"] == "NL"
    shown = {row["code"]: row["percent_of_gus"] for row in table["classes"]}
    assert shown["NL"] == 75.0
    assert shown["L1"] == 25.0


def test_no_class_description_promises_a_hormone_or_drug_exposure():
    for spec in GUS.CLASSES:
        text = f"{spec.substrates} {spec.does_not_mean}".lower()
        assert "your oestrogen" not in text
        assert "your estrogen" not in text
        assert spec.does_not_mean.strip(), spec.code


# --------------------------------------------------------------------------- #
# §6.4 functional redundancy — how many organisms hold each function
# --------------------------------------------------------------------------- #


def test_effective_carriers_distinguishes_spread_from_dominance():
    """A raw carrier count cannot. Two organisms at 50/50 is redundancy;
    two at 99/1 is one organism and a trace, and the effective count says
    so where a count of "2" would not."""
    from openbiota.extension import ecology as ECO

    assert ECO.effective_carriers([0.5, 0.5]) == pytest.approx(2.0)
    assert ECO.effective_carriers([0.99, 0.01]) < 1.1
    assert ECO.effective_carriers([1.0]) == pytest.approx(1.0)
    assert ECO.effective_carriers([]) is None, (
        "no resolved carrier is not the same as no redundancy"
    )


def test_a_reference_organism_is_not_a_carrier_unless_it_is_in_the_sample():
    """A read's best match being some organism's copy of a gene does not put
    that organism in the sample; the reference set may not hold the one that
    is. Counting those would inflate every redundancy reading."""
    from openbiota.extension import engine

    results = {
        "organism_inventory": {"organisms": [{"species": "Bacteroides_uniformis"}]},
        "panels": [{
            "name": "p", "metabolite": "Thing", "aggregate_from": ["A"],
            "accepted_fragments": 100,
            "genes": [{"entry_id": "A", "organisms": [
                {"organism": "Bacteroides uniformis", "fragments": 60},
                {"organism": "Some unrelated reference organism", "fragments": 40},
            ]}],
        }],
    }
    rows = engine._functional_redundancy(results)
    assert len(rows) == 1
    row = rows[0]
    assert row["n_resolved_carriers"] == 1, "only the organism in the sample carries it"
    assert row["effective_carriers"] == pytest.approx(1.0)
    assert row["carrier_coverage"] == pytest.approx(0.6), (
        "coverage must say how much of the signal the resolved carriers hold"
    )


def test_redundancy_uses_the_shared_taxon_normaliser():
    """AGORA2 writes '[Ruminococcus] gnavus' and the inventory writes it
    without brackets. That mismatch has already cost this project a silent
    join once."""
    from openbiota.extension import engine

    results = {
        "organism_inventory": {"organisms": [{"species": "Ruminococcus_gnavus"}]},
        "panels": [{
            "name": "p", "metabolite": "Thing", "aggregate_from": ["A"],
            "accepted_fragments": 10,
            "genes": [{"entry_id": "A", "organisms": [
                {"organism": "[Ruminococcus] gnavus", "fragments": 10},
            ]}],
        }],
    }
    row = engine._functional_redundancy(results)[0]
    assert row["n_resolved_carriers"] == 1, (
        "the bracketed AGORA2 name must join to the unbracketed inventory name"
    )


def test_only_aggregated_targets_contribute_carriers():
    """A decoy's organisms are the background a panel discriminates against,
    not carriers of its function."""
    from openbiota.extension import engine

    results = {
        "organism_inventory": {"organisms": [
            {"species": "A_one"}, {"species": "B_two"},
        ]},
        "panels": [{
            "name": "p", "metabolite": "Thing", "aggregate_from": ["KEEP"],
            "accepted_fragments": 20,
            "genes": [
                {"entry_id": "KEEP", "organisms": [{"organism": "A one", "fragments": 10}]},
                {"entry_id": "OTHER", "organisms": [{"organism": "B two", "fragments": 10}]},
            ],
        }],
    }
    row = engine._functional_redundancy(results)[0]
    assert row["n_resolved_carriers"] == 1
    assert [c["organism"] for c in row["carriers"]] == ["A one"]


# --------------------------------------------------------------------------- #
# §5.5's five neuroactive cards, and the boundaries they must keep
# --------------------------------------------------------------------------- #


def test_all_five_neuroactive_cards_of_section_5_5_exist():
    labels = " | ".join(c.label.lower() for c in CAP.NEUROACTIVE_CARDS)
    for required in ("dopamine", "serotonin", "acetylcholine", "norepinephrine",
                     "histamine"):
        assert required in labels, f"§5.5 card missing: {required}"
    assert len(CAP.NEUROACTIVE_CARDS) == 5


def test_every_card_states_where_its_evidence_stops():
    for card in CAP.NEUROACTIVE_CARDS:
        assert len(card.endpoint_boundary.split()) >= 15, card.card_id


def test_a_card_without_a_panel_must_say_why():
    """Otherwise a reader assumes the chemistry was searched for and absent,
    which is a different and much stronger claim."""
    with pytest.raises(CAP.CapacityError, match="must say why"):
        CAP.NeuroactiveCard(
            card_id="x", label="x", links_to=(), evidence="e",
            endpoint_boundary="b" * 90, sequence_panel=None,
        )
    for card in CAP.NEUROACTIVE_CARDS:
        if card.sequence_panel is None:
            assert card.no_panel_because, card.card_id


def test_serotonin_is_never_backed_by_the_tryptamine_enzyme():
    """The specification's sharpest instruction here: do not upgrade a
    tryptamine-producing enzyme to a serotonin detector merely because a
    homolog is retrievable."""
    card = next(c for c in CAP.NEUROACTIVE_CARDS if c.card_id == "neuro.serotonin")
    assert card.sequence_panel is None
    assert "tryptamine" in card.no_panel_because.lower()
    assert "different molecule" in card.no_panel_because
    # It may *link* the tryptamine reading, labelled as tryptamine.
    assert "tryptamine" in card.links_to


def test_acetylcholine_refuses_the_generic_acetyltransferase():
    card = next(c for c in CAP.NEUROACTIVE_CARDS if c.card_id == "neuro.acetylcholine")
    assert card.sequence_panel is None
    assert "acetyltransferase" in card.endpoint_boundary.lower()
    assert "not an acetylcholine assay" in card.endpoint_boundary


def test_norepinephrine_separates_microbial_handling_from_host_response():
    card = next(c for c in CAP.NEUROACTIVE_CARDS if c.card_id == "neuro.norepinephrine")
    assert "deconjugat" in card.evidence.lower()
    assert "no bacterial abundance converts" in card.endpoint_boundary.lower()


def test_histamine_reports_production_and_not_a_diagnosis():
    card = next(c for c in CAP.NEUROACTIVE_CARDS if c.card_id == "neuro.histamine")
    assert card.sequence_panel == "hdcA"
    assert "intolerance" in card.endpoint_boundary.lower()
    # Degradation is deliberately absent, with the reason recorded.
    assert "degradation" in card.no_panel_because.lower()


def test_the_cards_reach_the_view():
    view = CAP.summarise(["tyrDC"], searched=["tyrDC", "hdcA"])
    assert view["n_neuroactive_cards"] == 5
    states = {c["card_id"]: c["state"] for c in view["neuroactive_cards"]}
    assert states["neuro.dopamine"] == "gene_detected"
    assert states["neuro.histamine"] == "gene_not_detected"
    assert states["neuro.serotonin"] is None, (
        "a card with no panel has no detection state to report"
    )


# --------------------------------------------------------------------------- #
# §11.4's two coverage scores
# --------------------------------------------------------------------------- #


def test_the_coverage_formula_matches_the_specification():
    """C = 100·Σ(w·c)/Σw and P = 100·Σ(w·p)/Σw, worked by hand."""
    from openbiota.extension import synbiotic as SYN

    goals = [
        {"goal": "a", "weight": 2.0, "favourable": True, "pair_supported": True},
        {"goal": "b", "weight": 1.0, "favourable": True, "pair_supported": False},
        {"goal": "c", "weight": 1.0, "favourable": False, "pair_supported": False},
    ]
    out = SYN.coverage(goals, n_components=2)
    assert out["goal_coverage"] == pytest.approx(75.0)          # (2+1)/4
    assert out["pair_supported_goal_coverage"] == pytest.approx(50.0)  # 2/4


def test_pair_coverage_is_null_for_one_component_and_not_zero():
    """Zero would read as a combination that was tried and did nothing."""
    from openbiota.extension import synbiotic as SYN

    out = SYN.coverage(
        [{"goal": "a", "weight": 1.0, "favourable": True, "pair_supported": False}],
        n_components=1)
    assert out["goal_coverage"] == pytest.approx(100.0)
    assert out["pair_supported_goal_coverage"] is None
    assert "not applicable rather than zero" in out["pair_coverage_unavailable_reason"]


def test_both_scores_are_null_when_there_is_nothing_to_cover():
    from openbiota.extension import synbiotic as SYN

    out = SYN.coverage([], n_components=2)
    assert out["goal_coverage"] is None
    assert out["pair_supported_goal_coverage"] is None
    zero_weight = SYN.coverage(
        [{"goal": "a", "weight": 0.0, "favourable": True}], n_components=2)
    assert zero_weight["goal_coverage"] is None


def test_an_invalid_weight_is_refused_rather_than_coerced():
    from openbiota.extension import synbiotic as SYN

    for bad in (-1.0, float("inf"), float("nan")):
        with pytest.raises(SYN.CoverageError):
            SYN.coverage([{"goal": "a", "weight": bad, "favourable": True}],
                         n_components=2)


def test_a_conflicted_goal_is_counted_and_marked_not_averaged_away():
    """§11.4: the coverage display must not hide the distinction."""
    from openbiota.extension import synbiotic as SYN

    out = SYN.coverage([{
        "goal": "a", "weight": 1.0, "favourable": True, "pair_supported": True,
        "conflicted": True,
    }], n_components=2)
    assert out["n_conflicted"] == 1
    assert out["goals"][0]["conflicted"] is True


def test_every_covered_goal_shows_its_evidence_setting():
    """So a mouse result never wears the same percentage as a trial without
    saying which it was."""
    from openbiota.extension import synbiotic as SYN

    out = SYN.coverage([{"goal": "a", "weight": 1.0, "favourable": True}],
                       n_components=2)
    assert out["goals"][0]["evidence_setting"] == "unstated"


def test_neither_score_enters_the_ranking():
    from openbiota.extension import synbiotic as SYN

    out = SYN.coverage([{"goal": "a", "weight": 1.0, "favourable": True}],
                       n_components=2)
    assert out["enters_the_ranking"] is False
    refusal = out["what_these_are_not"]
    assert "Neither is a probability" in refusal
    assert "health score" in refusal
    assert "neither enters the ranking" in refusal


def test_a_substrate_that_suffices_covers_the_goal_but_not_as_a_pair():
    """The distinction the two scores exist to make."""
    from openbiota.extension import synbiotic as SYN

    class _Sub:
        substrate_id = "inulin"

    class _Cand:
        candidate_id, substrate = "c1", _Sub()
        verdict = "substrate_suffices"

    out = SYN.coverage_for(_Cand())
    assert out["goal_coverage"] == pytest.approx(100.0)
    assert out["pair_supported_goal_coverage"] == pytest.approx(0.0), (
        "the fibre does the work, so the pair supports nothing extra"
    )
    _Cand.verdict = "pair_gains"
    assert SYN.coverage_for(_Cand())["pair_supported_goal_coverage"] == pytest.approx(100.0)


def test_a_model_result_is_never_labelled_as_a_trial():
    from openbiota.extension import synbiotic as SYN

    class _Sub:
        substrate_id = "inulin"

    class _Cand:
        candidate_id, substrate, verdict = "c1", _Sub(), "pair_gains"

    setting = SYN.coverage_for(_Cand())["goals"][0]["evidence_setting"]
    assert "not a trial" in setting
    assert "in-silico" in setting
