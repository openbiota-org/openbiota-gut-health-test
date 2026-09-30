"""Spec v04 report behaviour: coherence notes, assumed context, polarity colour, ordering.

Each test pins a promise the report makes to its reader:

* two measurements of one function never contradict each other silently —
  a divergence is graded and explained from both sides with identical facts;
* subject context that was not supplied is assumed absent (research and
  clinician modes), the profile is scored, and the assumption is recorded —
  never a silent skip;
* an organism's colour is a judgement of level × curated polarity, so a
  beneficial species that is unusual to carry is green, not a warning;
* the document order is community → age → functions → patterns, in both
  Part A and Part B, and every cross-reference comes from one table.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from openbiota.coherence import (
    VERDICT_AGREE,
    VERDICT_DIVERGENT,
    VERDICT_OPPOSED,
    check_coherence,
)
from openbiota.context import resolve_context
from openbiota.pdfatlas import polarity_judgement, polarity_word
from openbiota.pdfreport import AMBER, CORAL, GREEN, SECTIONS, SLATE, Tile
from openbiota.profiles import parse_profile
from openbiota.scoring import (
    ReferenceBundle,
    SampleMeasurements,
    score_profile,
    summarise_reference,
)

# --------------------------------------------------------------------------- #
# coherence: gene panel vs taxon group for the same function
# --------------------------------------------------------------------------- #


def test_butyrate_genes_and_producers_agree_when_both_low() -> None:
    report = check_coherence(panel_percentiles={"butyrate": 5.0}, group_percentiles={"butyrate": 4.0})
    c = next(x for x in report.results if x.pair.panel == "butyrate")
    assert c.verdict == VERDICT_AGREE
    assert not report.divergent


def test_butyrate_split_is_graded_and_explained_from_both_sides() -> None:
    # The SAMPLE2_A02 case before the scope fix: genes 89th, producers 4th.
    report = check_coherence(
        panel_percentiles={"butyrate": 89.0}, group_percentiles={"butyrate": 4.0},
        out_of_scope_fractions={"butyrate": 0.66},
    )
    c = next(x for x in report.results if x.pair.panel == "butyrate")
    assert c.verdict in (VERDICT_DIVERGENT, VERDICT_OPPOSED)
    assert c in report.divergent
    cap, car = c.explanation(side="capacity"), c.explanation(side="carriers")
    # Same facts under both readings: both percentiles named on both sides.
    for text in (cap, car):
        assert "89th" in text and "4th" in text
        assert "\u2026" not in text and "..." not in text
    # A mechanism is offered, not just a shrug.
    assert c.mechanisms
    assert c.out_of_scope_fraction == pytest.approx(0.66)


def test_coherence_skips_pairs_that_were_not_measured() -> None:
    report = check_coherence(panel_percentiles={"butyrate": None}, group_percentiles={"butyrate": 4.0})
    assert all(x.pair.panel != "butyrate" for x in report.results)


# --------------------------------------------------------------------------- #
# assumed context: metformin and friends
# --------------------------------------------------------------------------- #


def _metformin_gated_profile():
    return parse_profile(
        {
            "name": "t2dlike",
            "version": "0.1.0",
            "label": "T2D-like",
            "status": "research_only",
            "summary": "s",
            "citation": "c",
            "mandatory_match": ["metformin"],
            "modules": {
                "taxonomic": {
                    "weight_in_combined": 1.0,
                    "features": [{"name": "Test_species", "level": "species", "direction": "decreased"}],
                }
            },
            "abstain_if": [{"missing_metadata": ["metformin"]}],
        }
    )


def _bundle(n: int = 500) -> ReferenceBundle:
    values = [float(v) for v in range(n)]
    bundle = ReferenceBundle(
        taxonomic={"Test_species": summarise_reference("Test_species", values, "t", prevalence=1.0)},
        taxonomic_source="test",
        taxonomic_n=n,
    )
    bundle.match_info = {"matched": True, "reason": "matched"}
    return bundle


def test_research_mode_assumes_no_metformin_and_records_it() -> None:
    ledger = resolve_context(mode="research")
    assert "metformin" in ledger.assumed_keys
    assert ledger.value("metformin") is False
    a = ledger.assumption_for("metformin")
    assert a is not None and "metformin" in a.sentence().lower()
    # The consequence is spelled out for the reader.
    assert a.field.if_present


def test_participant_mode_does_not_assume() -> None:
    ledger = resolve_context(mode="participant")
    assert "metformin" not in ledger.assumed_keys
    assert ledger.value("metformin") is None


def test_supplied_medication_overrides_the_assumption() -> None:
    ledger = resolve_context(mode="research", medications=["Glucophage 500 mg"])
    assert ledger.value("metformin") is True
    assert "metformin" not in ledger.assumed_keys


def test_metformin_gated_profile_scores_on_the_assumption_in_research_mode() -> None:
    ledger = resolve_context(mode="research")
    result = score_profile(
        profile=_metformin_gated_profile(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=10_000_000,
            metadata=ledger.metadata(),
            assumed_keys=ledger.assumed_keys,
        ),
        references=_bundle(),
    )
    assert not result.abstention.abstained
    assert "metformin" in result.abstention.assumptions
    assert result.combined_percentile is not None


def test_metformin_gated_profile_still_abstains_when_nothing_is_assumed() -> None:
    result = score_profile(
        profile=_metformin_gated_profile(),
        measurements=SampleMeasurements(
            taxonomic_clr={"Test_species": 10.0},
            taxonomic_abundance={"Test_species": 1.0},
            usable_nonhost_reads=10_000_000,
            metadata={},
        ),
        references=_bundle(),
    )
    assert result.abstention.abstained
    assert any("metformin" in r for r in result.abstention.triggered)


# --------------------------------------------------------------------------- #
# organism colour = level × polarity
# --------------------------------------------------------------------------- #


def _finding(bucket: int, polarity: str, guilds: tuple[str, ...] = ()) -> SimpleNamespace:
    return SimpleNamespace(bucket=bucket, level=None, guilds=guilds, interpretation=SimpleNamespace(polarity=polarity))


@pytest.mark.parametrize(
    ("bucket", "polarity", "guilds", "expected"),
    [
        (4, "health_associated", (), CORAL),            # beneficial species missing: red
        (3, "health_associated", (), AMBER),            # beneficial species low: amber
        (2, "health_associated", (), GREEN),            # beneficial species high: green
        (5, "health_associated", (), GREEN),            # beneficial species unusual: green
        (5, "context_dependent", ("probiotics",), GREEN),  # L. reuteri-style probiotic: green
        (2, "no_curated_interpretation", ("butyrate",), GREEN),  # butyrate producer high: green
        (5, "no_curated_interpretation", (), SLATE),    # unknown unusual: grey, never orange
        (5, "context_dependent", (), SLATE),            # neutral unusual: grey
        (2, "disease_associated", (), CORAL),           # adverse species high: red
        (5, "potential_pathobiont", (), CORAL),         # adverse species unusual: red
        (4, "disease_associated", (), GREEN),           # adverse species missing: green
        (6, "potential_pathobiont", (), CORAL),         # needs qualification: always red
        (1, "health_associated", (), SLATE),            # typical: no verdict
    ],
)
def test_polarity_judgement_colours(bucket, polarity, guilds, expected) -> None:
    colour, _bg = polarity_judgement(_finding(bucket, polarity, guilds))
    assert colour is expected


def test_polarity_word_prefers_curation_then_guild() -> None:
    assert polarity_word(_finding(5, "health_associated", ("opportunistic_pathogens",))) == "beneficial"
    assert polarity_word(_finding(5, "no_curated_interpretation", ("opportunistic_pathogens",))) == "adverse"
    assert polarity_word(_finding(5, "no_curated_interpretation", ())) == "unknown"
    assert polarity_word(_finding(5, "conflicting_evidence", ())) == "neutral"


# --------------------------------------------------------------------------- #
# document order and layout
# --------------------------------------------------------------------------- #


def test_part_a_and_part_b_follow_the_same_order() -> None:
    a = [SECTIONS[k] for k in ("community", "groups", "organisms", "age", "functions", "patterns", "actions")]
    b = [SECTIONS[k] for k in ("groups_detail", "organisms_detail", "functions_detail", "patterns_detail")]
    assert a == sorted(a) and b == sorted(b)
    assert SECTIONS["guide"] < SECTIONS["community"] < SECTIONS["age"] < SECTIONS["functions"] < SECTIONS["patterns"]
    assert max(a) < min(b)
    assert len(set(SECTIONS.values())) == len(SECTIONS)


def test_tile_value_never_collides_with_its_label(tmp_path) -> None:
    """The value's baseline must sit above the label's cap height at every tile height."""
    from reportlab.lib.units import mm
    from reportlab.pdfgen.canvas import Canvas

    drawn: dict[str, list[tuple[float, float, str]]] = {"strings": []}

    class Spy(Canvas):
        def drawString(self, x, y, text, *a, **k):  # noqa: N802 — reportlab API
            drawn["strings"].append((x, y, text))
            return super().drawString(x, y, text, *a, **k)

    for height in (12 * mm, 13 * mm, 15 * mm, 20 * mm):
        drawn["strings"].clear()
        c = Spy(str(tmp_path / "t.pdf"))
        t = Tile(value="3.17", label="Diversity", note="Shannon index", width=40 * mm, height=height)
        t.canv = c
        t.draw()
        ys = {text: y for _x, y, text in drawn["strings"]}
        assert ys["3.17"] > ys["DIVERSITY"] + 4.5, f"value overlaps label at height {height / mm:.0f} mm"
        assert ys["DIVERSITY"] > ys["Shannon index"]
