"""What may be offered as an option, and how a finding reaches one.

These run on synthetic records and on the shipped registry. None of them
needs a sample, private reads or a finished run.
"""

from __future__ import annotations

from typing import Any

from openbiota import delivery as DEL
from openbiota import solutions as sol

# --------------------------------------------------------------------------- #
# the evidence rule
# --------------------------------------------------------------------------- #


def test_no_study_design_is_excluded() -> None:
    """A dish, a mouse and a trial are all real findings about the organism.

    The report used to lean on human trials and let the rest fall off the
    page, which threw away most of what is known about botanicals, enzymes
    and antimicrobial food compounds. Every setting is reportable; what
    changes between them is the word beside it and a small marker.
    """
    for setting in ("human_trial", "human_observational", "animal", "ex_vivo", "in_vitro"):
        reach = DEL.assess(setting, "oral_tested")
        assert reach.reportable, f"{setting} must be reportable"
        assert reach.setting_words


def test_a_laboratory_result_is_marked_but_not_hidden() -> None:
    """The caveat is a superscript and a footnote, never a reason to drop it."""
    lab = DEL.assess("in_vitro", "inherently_stable")
    assert lab.reportable
    assert lab.markers == ("\u2021",)
    assert DEL.footnotes_for(set(lab.markers))[0].startswith("\u2021")

    trial = DEL.assess("human_trial", "oral_tested")
    assert trial.reportable
    assert trial.markers == (), "a human finding carries no marker"


def test_delivery_is_what_gates_a_culture_result() -> None:
    """A dish result counts when the material can arrive in the gut intact.

    By its own chemistry, or in an enteric-protected preparation - the only
    thing that disqualifies it is not being swallowed at all.
    """
    assert DEL.assess("in_vitro", "inherently_stable").reportable
    assert DEL.assess("in_vitro", "oral_tested").reportable

    enteric = DEL.assess("in_vitro", "enteric_required")
    assert enteric.reportable, "an enteric capsule is an ordinary thing to buy"
    assert "\u00a7" in enteric.markers

    assert not DEL.assess("in_vitro", "not_ingestible").reportable


def test_a_study_in_a_living_body_answers_the_delivery_question_itself() -> None:
    """A mouse that ate it has shown the material got where it acts."""
    fed = DEL.assess("animal", "not_ingestible")
    assert fed.reportable, "the animal was dosed; the route is not in doubt"
    assert "\u00a7" not in fed.markers


def test_footnotes_are_only_the_ones_earned() -> None:
    """A page with no animal evidence carries no line about animal evidence."""
    lab_marker = "\u2021"
    assert DEL.footnotes_for({lab_marker}) == [lab_marker + " " + DEL.FOOTNOTES[lab_marker]]
    assert DEL.footnotes_for(set()) == []


# --------------------------------------------------------------------------- #
# matching a finding to the research
# --------------------------------------------------------------------------- #


def test_the_crosswalk_reaches_research_written_under_the_old_name() -> None:
    """GTDB and the literature disagree, and the report has to bridge it.

    The organism this report calls *Mediterraneibacter gnavus* is the one the
    papers call *Ruminococcus gnavus*. Without the inventory's own `formerly`
    field none of that research would be found at all.
    """
    inventory = {
        "organisms": [
            {
                "species": "Mediterraneibacter_gnavus",
                "gtdb": "Mediterraneibacter gnavus",
                "formerly": "Ruminococcus_gnavus",
                "native_ids": {"scoring": "Ruminococcus_gnavus"},
            }
        ]
    }
    keys = sol.alias_map(inventory)[sol.norm("Mediterraneibacter_gnavus")]
    assert sol.norm("Ruminococcus gnavus") in keys


def test_gtdb_suffixes_and_underscores_do_not_block_a_match() -> None:
    """`Blautia_A luti` and `Blautia luti` are the same organism to a reader."""
    assert sol.norm("Blautia_A luti") == sol.norm("Blautia luti") == "blautia luti"
    assert sol.norm("Mediterraneibacter_gnavus") == "mediterraneibacter gnavus"


# --------------------------------------------------------------------------- #
# what may be offered as a solution
# --------------------------------------------------------------------------- #


def _option(**kw: Any) -> dict[str, Any]:
    base = {
        "option_id": "T-1", "material": "Test material", "kind": "herb",
        "setting": "in_vitro", "delivery": "inherently_stable",
        "direction": "decrease", "effect": "inhibited it", "boundary": "",
        "taxa": ["Testus organismus"], "genera": [], "panels": [], "concerns": [],
    }
    base.update(kw)
    return base


def test_an_option_that_would_worsen_a_finding_is_never_offered_as_a_fix(monkeypatch: Any) -> None:
    """A fibre shown to *feed* an overgrown organism is not a way to reduce it.

    This is the failure the direction rule exists to prevent: the row is real
    and belongs in the report, on the other side of the problem.
    """
    rows = (
        _option(option_id="DOWN", material="Lowers it", direction="decrease"),
        _option(option_id="UP", material="Feeds it", direction="increase"),
        _option(option_id="UNK", material="Unstated", direction="unclear"),
    )
    monkeypatch.setattr(sol, "_atlas", lambda: rows)
    sol._index.cache_clear()

    down = [o.option_id for o in sol.atlas_options(species=["Testus organismus"], want="decrease")]
    up = [o.option_id for o in sol.atlas_options(species=["Testus organismus"], want="increase")]
    sol._index.cache_clear()

    assert down == ["DOWN"]
    assert up == ["UP"]
    assert "UNK" not in down and "UNK" not in up, (
        "a row that never says which way it went is not a way to fix anything"
    )


def test_a_row_naming_several_organisms_uses_the_direction_for_each(monkeypatch: Any) -> None:
    """One direction for a whole cascade row is wrong for at least one of them."""
    rows = (
        _option(
            option_id="CASCADE", material="Inulin in a six-species community",
            direction="decrease",
            taxa=["Ruminococcus gnavus", "Bacteroides uniformis"],
            taxon_direction={"Ruminococcus gnavus": "increase",
                             "Bacteroides uniformis": "decrease"},
        ),
    )
    monkeypatch.setattr(sol, "_atlas", lambda: rows)
    sol._index.cache_clear()

    gnavus_up = sol.atlas_options(species=["Ruminococcus gnavus"], want="increase")
    gnavus_down = sol.atlas_options(species=["Ruminococcus gnavus"], want="decrease")
    uniformis_down = sol.atlas_options(species=["Bacteroides uniformis"], want="decrease")
    sol._index.cache_clear()

    assert [o.option_id for o in gnavus_up] == ["CASCADE"]
    assert gnavus_down == []
    assert [o.option_id for o in uniformis_down] == ["CASCADE"]


def test_a_problem_with_nothing_on_record_says_so(monkeypatch: Any) -> None:
    """Silence is a finding too, and is never left as an empty table."""
    monkeypatch.setattr(sol, "_atlas", tuple)
    sol._index.cache_clear()
    problems = sol.panel_problems(
        [{"panel": "butyrate", "percentile": 2.0, "metabolite": "Butyrate"}], section=19,
    )
    sol._index.cache_clear()
    assert len(problems) == 1
    assert problems[0].nothing_known
    assert problems[0].options == []


# --------------------------------------------------------------------------- #
# which findings become problems at all
# --------------------------------------------------------------------------- #


def test_a_pathway_is_a_problem_only_in_the_direction_that_matters() -> None:
    """Low butyrate is a finding; high butyrate is not. The reverse for H2S."""
    low_butyrate = sol.panel_problems(
        [{"panel": "butyrate", "percentile": 3.0, "metabolite": "Butyrate"}], section=19)
    high_butyrate = sol.panel_problems(
        [{"panel": "butyrate", "percentile": 97.0, "metabolite": "Butyrate"}], section=19)
    high_h2s = sol.panel_problems(
        [{"panel": "h2s", "percentile": 97.0, "metabolite": "Hydrogen sulfide"}], section=19)
    low_h2s = sol.panel_problems(
        [{"panel": "h2s", "percentile": 3.0, "metabolite": "Hydrogen sulfide"}], section=19)

    assert len(low_butyrate) == 1 and low_butyrate[0].want == sol.WANT_UP
    assert high_butyrate == []
    assert len(high_h2s) == 1 and high_h2s[0].want == sol.WANT_DOWN
    assert low_h2s == []


def test_a_typical_pathway_is_not_a_problem() -> None:
    """Sitting in the middle of the reference is not something to fix."""
    assert sol.panel_problems(
        [{"panel": "butyrate", "percentile": 50.0}], section=19) == []


def test_the_band_follows_the_reports_own_seven_band_scale() -> None:
    assert sol.band_of(3.0) == "notably_low"
    assert sol.band_of(10.0) == "low"
    assert sol.band_of(50.0) == "typical"
    assert sol.band_of(90.0) == "high"
    assert sol.band_of(99.0) == "notably_high"
    assert sol.band_of(None) == "unknown"


def test_problems_are_ordered_worst_first_and_answerable_first() -> None:
    """A pathogen outranks a pathway; within a tier, what can be acted on leads."""
    def problem(kind: str, *, known: bool, severity: float = 1.0) -> sol.Problem:
        p = sol.Problem(key=kind + str(severity), kind=kind, title=kind, detail="",
                        want=sol.WANT_DOWN, severity=severity, section=1)
        p.nothing_known = not known
        return p

    order = sol.ordered([
        problem("panel", known=True),
        problem("organism", known=False),
        problem("organism", known=True),
        problem("pathogen", known=True),
    ])
    assert [p.kind for p in order] == ["pathogen", "organism", "organism", "panel"]
    assert order[1].nothing_known is False


# --------------------------------------------------------------------------- #
# the shipped registry
# --------------------------------------------------------------------------- #


def test_the_imported_atlas_is_usable_if_it_is_present() -> None:
    """When the generated registry is in the checkout, it has to be coherent."""
    rows = sol._atlas()
    if not rows:  # not generated in this checkout
        return
    assert len(rows) > 200
    for row in rows:
        assert row.get("material"), "every option names its material"
        assert row.get("setting") in DEL.SETTINGS
        assert row.get("delivery") in DEL.DELIVERIES
        assert row.get("direction") in {"increase", "decrease", "null", "unclear"}
    bound = sum(1 for r in rows
                if r.get("taxa") or r.get("genera") or r.get("panels") or r.get("concerns"))
    assert bound > len(rows) // 2, "most options must be attached to something"


def test_laboratory_and_animal_evidence_actually_survives_into_the_registry() -> None:
    """The whole point of the rule: this is where most of the research is."""
    rows = sol._atlas()
    if not rows:
        return
    settings = {r["setting"] for r in rows}
    assert "in_vitro" in settings
    assert "animal" in settings
    assert sum(1 for r in rows if r["setting"] == "in_vitro") > 50


# --------------------------------------------------------------------------- #
# finding your way between a reading and what to do about it
# --------------------------------------------------------------------------- #


def test_a_reading_and_the_actions_section_agree_on_one_identity() -> None:
    """Both ends build the key from the organism, so neither has to be told."""
    from openbiota import pdflinks

    key = sol.problem_key("organism", "Mediterraneibacter_gnavus")
    assert key == "organism:mediterraneibacter gnavus"
    # A GTDB suffix must not produce a second, unreachable identity.
    assert sol.problem_key("organism", "Blautia_A luti") == sol.problem_key("organism", "Blautia luti")
    assert pdflinks.solution_dest(key) == "sol-organism-mediterraneibacter-gnavus"


def test_a_reading_links_on_only_when_there_is_something_at_the_far_end() -> None:
    """A link to a page the report never drew would crash the build."""
    from openbiota import pdflinks

    pdflinks.reset_solutions()
    key = sol.problem_key("panel", "butyrate")
    assert pdflinks.solutions_line(key) is None
    pdflinks.register_solutions([key])
    assert pdflinks.solutions_line(key) is not None
    pdflinks.reset_solutions()
    assert pdflinks.solutions_line(key) is None


def test_a_problem_points_at_the_reading_rather_than_the_section() -> None:
    """"See section 16" is no answer when section 16 is twenty pages long."""
    from openbiota import pdflinks

    problems = sol.panel_problems(
        [{"panel": "butyrate", "percentile": 3.0, "metabolite": "Butyrate"}], section=19)
    assert problems[0].detail_dest == pdflinks.function_dest("butyrate", detail=True)

    pathogens = sol.pathogen_problems(
        [{"species": "Clostridioides difficile", "target_id": "bacteria.clostridioides_difficile"}],
        section=18)
    assert pathogens[0].detail_dest == pdflinks.pathogen_dest(
        "bacteria.clostridioides_difficile", detail=True)


# --------------------------------------------------------------------------- #
# what the report may say about putting organisms back
# --------------------------------------------------------------------------- #


def test_the_report_never_claims_nothing_can_restore_an_organism() -> None:
    """Donor strains demonstrably establish and persist, so saying otherwise is wrong.

    The card for a missing commensal used to be headed "no supported way to
    'add it back'". Strain-resolved tracking puts 71% of donor strains still
    present five years after transplant, so a route plainly exists; what does
    not exist is a way to put back one chosen organism. The report must not
    tell someone there is no option.
    """
    import pathlib

    from openbiota.interventions import load_registry

    reg = load_registry(pathlib.Path("interventions"))
    maps = {m.action_mapping_id: m for m in reg.mappings}
    for mapping in maps.values():
        prose = f"{mapping.title} {mapping.context_note}".lower()
        assert "no supported way to 'add it back'" not in prose
        assert "no supported way to add it back" not in prose

    card = maps["MAP_TAXON_ABSENT_FMT_ENGRAFTMENT_V1"]
    note = card.context_note.lower()
    # It says the route exists, and cites the study that shows it.
    assert "persist" in note or "establish" in note
    assert "71%" in note and "34580445" in note, "the engraftment finding must be cited inline"
    # It must not read as though transfer only works for one disease: the
    # protocol it hangs off is the engraftment evidence, not the one approved
    # indication, and the note says so.
    assert card.intervention_protocol_ref == "PROTOCOL_FMT_ENGRAFTMENT_V1"
    assert "ten different conditions" in note or "not something that happens only in one" in note
    # ...and every limit still travels with it.
    assert "selective" in note
    assert "multidrug-resistant" in note
    assert card.fmt_status, "an FMT card must declare its regulatory status"
    assert any("DIY" in s or "home fecal" in s for s in card.never_say)


def test_a_studied_material_is_never_silently_dropped() -> None:
    """A source that did not state a direction is not a reason to hide it.

    A third of everything on record for an organism says what was given and
    what was measured but never which way the reading moved. Those rows are
    not claimed to help - they go in their own group - but they are shown.
    """
    rows = (
        _option(option_id="DOWN", material="Lowers it", direction="decrease"),
        _option(option_id="UNK", material="Studied, direction unstated", direction="unclear"),
    )
    import pytest as _pytest  # noqa: PLC0415

    mp = _pytest.MonkeyPatch()
    mp.setattr(sol, "_atlas", lambda: rows)
    sol._index.cache_clear()
    unclear = sol.atlas_options(species=["Testus organismus"], want="unclear")
    sol._index.cache_clear()
    mp.undo()
    assert [o.option_id for o in unclear] == ["UNK"]


def test_every_option_carries_its_paper_when_the_source_has_one() -> None:
    """800 citations were in the data and none reached the page."""
    rows = sol._atlas()
    if not rows:
        return
    with_ref = [r for r in rows if r.get("pmid") or r.get("doi") or r.get("url")]
    assert len(with_ref) > len(rows) // 2
    sample = sol._as_option(with_ref[0], binding="species")
    assert sample.source or sample.citation


def test_a_citation_is_never_mistaken_for_a_link() -> None:
    """"Ha CWY et al. Cell 2020" is not a URL, and linking it crashed the build."""
    summary = {"supporting": [{"url": "Ha CWY et al. Cell 2020;183:666", "citation": "Ha CWY et al."}]}
    assert sol._first_source(summary) == ""
    assert sol._citation_of(summary) == "Ha CWY et al."
    ok = {"supporting": [{"url": "https://pubmed.ncbi.nlm.nih.gov/1/", "citation": "A et al."}]}
    assert sol._first_source(ok).startswith("https://")


def test_a_card_that_proposes_nothing_is_not_counted_as_an_option() -> None:
    """"No supported targeted intervention" is not a thing you can do.

    Rendering those cards as options put that sentence in the column of
    things that might help and inflated the count of studied options on
    page one.
    """
    results = {
        "findings_and_evidence": {
            "triggers": [
                {"trigger_id": "T1", "trigger_type": "functional_capacity_finding",
                 "target": "butyrate"},
            ],
            "evidence_summaries": [
                {"summary_id": "S1", "render": "no_supported_targeted_intervention",
                 "title": "Something — no supported targeted intervention",
                 "intervention_display": None, "trigger_ids": ["T1"],
                 "context_note": "Why there is nothing."},
                {"summary_id": "S2", "render": "evidence_only", "title": "A real card",
                 "intervention_display": "An exact studied product", "trigger_ids": ["T1"],
                 "supporting_evidence_lane": "C", "intervention_class": "probiotic"},
            ],
        }
    }
    options = sol.registry_options(results)
    labels = [o.label for v in options.values() for o in v]
    assert labels == ["An exact studied product"]
    # The reasoning is kept, and lands on the problem instead.
    assert sol.registry_notes(results)["T1"] == ["Why there is nothing."]
