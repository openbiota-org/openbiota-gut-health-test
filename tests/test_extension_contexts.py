"""A13 symptom and condition navigation — spec 0.8.3 sections 9.3 and 16.3.

The rules worth pinning are the ones a report drifts away from: that a
context never becomes a score, that a count of links never becomes a
likelihood, and that no contributor is hidden behind a truncation.
"""

from __future__ import annotations

import pytest

from openbiota.extension import contexts as CX

#: The fifteen contexts named in section 16.3, by the words the spec uses.
REQUIRED_TOPICS = (
    "anxiety", "autism", "serotonin", "atheroscler", "blood pressure",
    "ethanol", "constipation", "diarrhoea", "abdominal pain", "gluten",
    "irritable bowel", "thyroid", "allergy", "glp-1", "eczema",
)


def test_all_fifteen_contexts_of_section_16_3_are_present():
    labels = " | ".join(c.label.lower() for c in CX.CONTEXTS)
    missing = [topic for topic in REQUIRED_TOPICS if topic not in labels]
    assert not missing, f"contexts required by 16.3 and not implemented: {missing}"
    assert len(CX.CONTEXTS) == 15


def test_context_ids_are_unique():
    ids = [c.context_id for c in CX.CONTEXTS]
    assert len(ids) == len(set(ids))


def test_every_context_states_the_distinction_it_keeps():
    for context in CX.CONTEXTS:
        assert len(context.distinction.split()) >= 15, context.context_id


def test_every_context_says_what_it_does_not_know():
    for context in CX.CONTEXTS:
        assert context.missing_information, context.context_id
        for item in context.missing_information:
            assert len(item.split()) >= 6, f"{context.context_id}: {item!r}"


def test_a_context_cannot_be_declared_without_links():
    with pytest.raises(CX.ContextError, match="no links is a heading"):
        CX.Context(
            context_id="x", label="x", links=(), distinction="d" * 90,
            mechanisms=(), missing_information=("a b c d e f g",), follow_up=(),
        )


def test_a_link_must_say_why_it_is_relevant():
    with pytest.raises(CX.ContextError, match="must say why"):
        CX.Link("panel", "butyrate", "  ")


def test_a_link_kind_must_be_one_the_resolver_understands():
    """A typo here would silently produce a context that links to nothing."""
    with pytest.raises(CX.ContextError, match="unknown link kind"):
        CX.Link("pannel", "butyrate", "typo")


def test_no_context_invents_a_probability_or_a_diagnosis():
    # Phrases that assert something about the reader, as opposed to ones
    # that merely contain the same words: "if you have a measured TMAO" is
    # an invitation to supply a lab result, not a claim about anybody.
    banned = (
        "your risk of", "suggests you have", "indicates you have",
        "means you have", "likely that you", "probability of", "% chance",
        "you are at risk", "consistent with a diagnosis",
    )
    for context in CX.CONTEXTS:
        text = " ".join([
            context.distinction, *context.mechanisms,
            *context.missing_information, *context.follow_up,
        ]).lower()
        for phrase in banned:
            # "diagnosed", "diagnosis" appear only when refusing one, which
            # is the opposite of claiming it.
            if phrase == "diagnos":
                continue
            assert phrase not in text, f"{context.context_id}: {phrase!r}"


def test_contexts_that_reuse_an_existing_score_say_which_one():
    """Section 9.3 requires existing profiles to be linked unchanged rather
    than recomputed, and 16.3 requires IBS specifically to be preserved."""
    ibs = CX.BY_ID["context.ibs"]
    assert set(ibs.preserves) == {"ibs", "ibs_c", "ibs_d", "ibs_m"}
    for context in CX.CONTEXTS:
        for preserved in context.preserves:
            assert any(
                link.kind == "profile" and link.key == preserved
                for link in context.links
            ), f"{context.context_id} claims to preserve {preserved} but does not link it"


def test_hypertension_reuses_the_existing_result_rather_than_duplicating_it():
    context = CX.BY_ID["context.hypertension"]
    assert "hypertension" in context.preserves
    assert "does not produce a second number" in context.distinction


def test_the_ethanol_context_refuses_the_diagnosis_the_panel_invites():
    context = CX.BY_ID["context.auto_brewery_ethanol"]
    assert "not measured ethanol" in context.distinction
    assert any("blood alcohol" in item for item in context.missing_information)


def test_the_glp1_context_records_that_indole_has_no_single_direction():
    context = CX.BY_ID["context.glp1"]
    assert "one universal direction" in context.distinction


def test_resolve_lists_every_contributor_and_truncates_nothing():
    results = {
        "panels": [
            {"name": "butyrate", "copies_per_100_genomes": 12.0, "accepted_fragments": 40},
            {"name": "propionate", "copies_per_100_genomes": 3.0},
        ],
        "profile_similarity": {"profiles": {"hypertension": {"score": 1}}},
    }
    view = CX.resolve(results, {})
    card = next(c for c in view["contexts"] if c["context_id"] == "context.hypertension")
    assert card["n_contributing"] == 3
    assert len(card["contributing"]) == 3
    assert card["not_in_this_run"] == []
    keys = {row["key"] for row in card["contributing"]}
    assert keys == {"hypertension", "butyrate", "propionate"}


def test_a_link_that_this_run_did_not_produce_is_reported_not_dropped():
    view = CX.resolve({"panels": []}, {})
    card = next(c for c in view["contexts"] if c["context_id"] == "context.hypertension")
    assert card["n_contributing"] == 0
    assert {row["key"] for row in card["not_in_this_run"]} == {
        "hypertension", "butyrate", "propionate"
    }


def test_the_count_of_links_is_labelled_as_navigation_not_likelihood():
    view = CX.resolve({"panels": []}, {})
    for card in view["contexts"]:
        note = card["navigation_count_note"]
        assert "not a measure of how likely" in note
    assert "not a test for the condition" in view["standing_note"]


def test_reverse_navigation_finds_the_contexts_a_reading_feeds():
    topics = CX.contexts_for({}, "butyrate", "panel")
    assert "Blood pressure context" in topics
    assert "Eczema" in topics
    assert CX.contexts_for({}, "no_such_panel", "panel") == []
