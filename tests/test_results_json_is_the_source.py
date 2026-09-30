"""results.json is the source of truth; the report is drawn from it.

The rule this file enforces: the analysis writes every displayed value into
`results.json`, and a renderer derives nothing further. Two reasons, and
the second is the one that bites.

A renderer that recomputes can disagree with the file beside it, and then
two artefacts make different claims about the same sample with no way to
tell which is wrong. That already happened here: the percentile shown on a
slider came from one code path and the "not detected" chip from another,
and they contradicted each other for methane in all five samples.

And the file is the interface. A web client will read the same JSON and
render its own view, so anything the PDF can see and the JSON cannot is a
value that client will be missing.

These tests read the JSON and check it *before* anything is drawn. The PDF
is validated separately, afterwards, against the same file.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

SAMPLES = sorted(
    p.name for p in RESULTS.iterdir()
    if p.is_dir() and (p / "results.json").is_file()
) if RESULTS.is_dir() else []

#: Every field a consumer needs to draw one reading without deriving anything.
REQUIRED_ROW_FIELDS = (
    "panel", "metabolite", "value", "percentile", "detected", "fragments",
    "status", "what_it_is", "made_from", "higher_means", "evidence_strength",
    "summary", "evidence_detail", "citation", "genes", "confidence",
    "implication", "category",
)

#: And the status has to carry its meaning, not a colour.
REQUIRED_STATUS_FIELDS = ("label", "level", "assessed", "note")


def _results(sample: str) -> dict[str, Any]:
    return json.loads((RESULTS / sample / "results.json").read_text(encoding="utf-8"))


def _rows(sample: str) -> list[dict[str, Any]]:
    rows = _results(sample).get("report_rows")
    if not rows:
        pytest.skip(f"{sample} predates report_rows in results.json")
    return rows


# --------------------------------------------------------------------------- #
# the file holds everything the report draws
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_results_json_carries_a_row_for_every_panel(sample: str) -> None:
    results = _results(sample)
    rows = _rows(sample)
    panels = {p["name"] for p in results["panels"]}
    in_rows = {r["panel"] for r in rows}
    assert in_rows == panels, (
        f"panels and report rows disagree: only in panels {sorted(panels - in_rows)}, "
        f"only in rows {sorted(in_rows - panels)}"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_row_carries_every_field_a_client_needs(sample: str) -> None:
    for row in _rows(sample):
        missing = [f for f in REQUIRED_ROW_FIELDS if f not in row]
        assert not missing, f"{row.get('panel')}: missing {missing}"
        status = row["status"]
        missing = [f for f in REQUIRED_STATUS_FIELDS if f not in status]
        assert not missing, f"{row['panel']}: status missing {missing}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_status_carries_meaning_and_not_presentation(sample: str) -> None:
    """A colour is a choice the client makes. The label, the level and
    whether the reading was placed at all are facts, and those travel."""
    for row in _rows(sample):
        status = row["status"]
        assert isinstance(status["assessed"], bool), row["panel"]
        assert isinstance(status["level"], int), row["panel"]
        for banned in ("colour", "color", "background", "hex"):
            assert banned not in status, (
                f"{row['panel']}: {banned!r} is presentation and does not belong in the "
                "canonical result"
            )


# --------------------------------------------------------------------------- #
# and what it holds is self-consistent
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_row_that_was_not_placed_says_why_in_the_file(sample: str) -> None:
    for row in _rows(sample):
        if row["status"]["assessed"]:
            continue
        assert row["status"].get("unassessed_reason"), (
            f"{row['panel']}: not assessed, and the file does not say why"
        )


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_zero_value_is_written_as_a_zero_position(sample: str) -> None:
    """The contradiction that started this: a slider at the 28th percentile
    beside a chip reading 'not detected'. One file, one answer."""
    for row in _rows(sample):
        if row["value"] in (0, 0.0) and row["percentile"] is not None:
            assert row["percentile"] in (0, 0.0), (
                f"{row['panel']}: value {row['value']} but percentile "
                f"{row['percentile']} in results.json"
            )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_row_agrees_with_the_panel_it_describes(sample: str) -> None:
    results = _results(sample)
    panels = {p["name"]: p for p in results["panels"]}
    for row in _rows(sample):
        panel = panels[row["panel"]]
        assert row["value"] == panel["copies_per_100_genomes"], (
            f"{row['panel']}: the row and the panel disagree about the value"
        )
        assert row["detected"] == bool(panel["copies_per_100_genomes"]), row["panel"]


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_row_agrees_with_the_reference_comparison(sample: str) -> None:
    results = _results(sample)
    comparison = results["reference_comparison"]["panels"]
    for row in _rows(sample):
        recorded = comparison.get(row["panel"]) or {}
        assert row["percentile"] == recorded.get("percentile"), (
            f"{row['panel']}: report_rows and reference_comparison disagree about the "
            "percentile, so the report and the JSON would show different numbers"
        )


# --------------------------------------------------------------------------- #
# round trip: reading the file back changes nothing
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_row_survives_a_round_trip_through_the_file(sample: str) -> None:
    from openbiota.pdfreport import MetaboliteRow

    for payload in _rows(sample):
        row = MetaboliteRow.from_json(payload)
        again = row.to_json()
        for field in REQUIRED_ROW_FIELDS:
            if field == "genes":
                assert [list(g) for g in again[field]] == [
                    list(g) for g in payload[field]
                ], f"{payload['panel']}.{field}"
            elif field == "status":
                assert again[field]["label"] == payload[field]["label"]
                assert again[field]["assessed"] == payload[field]["assessed"]
            else:
                assert again[field] == payload[field], f"{payload['panel']}.{field}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_rows_from_json_returns_one_row_per_panel(sample: str) -> None:
    from openbiota.pdfreport import rows_from_json

    results = _results(sample)
    if not results.get("report_rows"):
        pytest.skip(f"{sample} predates report_rows")
    rows = rows_from_json(results)
    assert len(rows) == len(results["panels"])
    assert all(r.status.label for r in rows)


def test_rows_from_json_refuses_a_file_that_predates_the_contract():
    from openbiota.pdfreport import rows_from_json

    with pytest.raises(ValueError, match="no 'report_rows'"):
        rows_from_json({"panels": []})


# --------------------------------------------------------------------------- #
# the extension's own readings travel in the file too
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_functional_reading_and_its_position_are_in_the_file(sample: str) -> None:
    views = (_results(sample).get("extension") or {}).get("views") or {}
    scores = views.get("reading_scores")
    if not scores:
        pytest.skip(f"{sample} predates the reading scores")
    for rid, score in scores.items():
        assert "percentile" in score, rid
        if score["percentile"] is None:
            assert score.get("unscored_reason"), f"{rid}: unplaced with no reason given"
        else:
            assert score.get("limiting_gene"), f"{rid}: placed without naming why"
            assert 0.0 <= score["percentile"] <= 100.0, rid


# --------------------------------------------------------------------------- #
# the file is the interface: every section's data is in it
# --------------------------------------------------------------------------- #

#: Each block the report draws, and the fields inside it that carry the
#: numbers a reader sees. A web client reading this file must be able to
#: render the same sections, so anything the PDF shows and the file lacks is
#: a gap in the interface rather than a detail of the renderer.
SECTION_CONTRACT: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("report_rows", ()),
    ("reference_comparison", ("panels",)),
    ("panels", ()),
    ("panel_bands", ()),
    ("community_profile", ("phyla",)),
    ("organism_inventory", ("organisms",)),
    ("organism_verdicts", ()),
    ("profile_similarity", ("profiles", "ranked", "dysbiosis_anchor", "subject")),
    ("findings_and_evidence", ("findings", "triggers", "evidence_summaries", "registry")),
    ("microbiome_age", ("age_residual_years", "calibration_population")),
    ("mycobiome", ()),
    ("pathogens", ()),
    ("biofilm", ()),
    ("strain_resolution", ()),
    ("sequencing_quality", ()),
    ("extension", ("views", "metrics", "manifest")),
    ("metabolite_drivers", ()),
    ("standing_caveats", ()),
    ("subject_context", ()),
    ("normalisation", ()),
)


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_section_the_report_draws_is_present_in_the_file(sample: str) -> None:
    _rows(sample)  # skips a file written before the contract existed
    results = _results(sample)
    missing: list[str] = []
    for block, fields in SECTION_CONTRACT:
        if block not in results:
            missing.append(block)
            continue
        payload = results[block]
        for field in fields:
            if not isinstance(payload, dict) or field not in payload:
                missing.append(f"{block}.{field}")
    assert not missing, (
        "a web client reading this file could not render these, because they are "
        f"not in it: {missing}"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_disease_patterns_carry_their_own_scores(sample: str) -> None:
    """The overview prints a percentile per pattern; it has to be in the file
    and not computed on the way to the page."""
    similarity = _results(sample)["profile_similarity"]
    for name, record in (similarity["profiles"] or {}).items():
        if str(record.get("status")) != "scored":
            continue
        combined = record.get("combined") or {}
        for field in ("score", "percentile"):
            assert field in combined, f"{name}: combined.{field} is not in the file"
    ranked = similarity.get("ranked") or []
    assert ranked, "the ranked list the overview reads is absent"
    for row in ranked:
        assert {"profile", "label", "status"} <= set(row), row


#: The fields inside each view that the report draws a whole section from.
#:
#: Naming the view is not enough: a client that has `ecology` but not
#: `ecology.redundancy` still cannot draw the page the PDF draws. Each entry
#: here is a section that exists on paper, so it has to exist in the file.
VIEW_FIELD_CONTRACT: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("ecology", ("redundancy",)),
    ("explorer", ("named_targets", "n_named_targets", "n_named_targets_detected")),
    ("simulation", ("coverage", "evidence_seeds")),
    ("fermentation", ("network",)),
    ("substrates", ("shared_genes",)),
    ("biotransformation", ("gaba_balance", "negative_controls")),
    ("nitrogen", ("three_different_concepts", "host_steps")),
)


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_view_field_a_section_is_drawn_from_is_in_the_file(sample: str) -> None:
    """A section on paper with nothing behind it in the file is a section the
    web client cannot rebuild."""
    _rows(sample)  # skips a file written before the contract existed
    views = (_results(sample)["extension"] or {}).get("views") or {}
    missing: list[str] = []
    for view, fields in VIEW_FIELD_CONTRACT:
        payload = views.get(view)
        if not isinstance(payload, dict):
            missing.append(view)
            continue
        missing.extend(f"{view}.{f}" for f in fields if f not in payload)
    assert not missing, (
        f"the report draws these and the file does not carry them: {missing}"
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_searched_target_carries_its_zero_in_the_file(sample: str) -> None:
    """The renderer must not be the place a measured absence becomes a number."""
    _rows(sample)
    views = (_results(sample)["extension"] or {}).get("views") or {}
    for target in (views.get("explorer") or {}).get("named_targets") or []:
        state = target.get("detection_state")
        if state == "no_supported_detection":
            assert target.get("percent") == 0, (
                f"{target.get('query')} was searched for and its absence is not "
                "written as a value, so the page would have to invent one"
            )
        elif state == "not_assessable":
            assert target.get("percent") is None, (
                f"{target.get('query')} could not be assessed yet carries a number"
            )


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_extension_views_are_all_in_the_file(sample: str) -> None:
    views = (_results(sample)["extension"] or {}).get("views") or {}
    expected = {
        "substrates", "fermentation", "vitamins", "nitrogen", "biotransformation",
        "additional_capacities", "contexts", "ecology", "explorer", "resistance",
        "planner", "longitudinal", "simulation", "input_register",
    }
    missing = sorted(expected - set(views))
    assert not missing, f"views the report renders and the file lacks: {missing}"


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_block_smuggles_presentation_into_the_canonical_result(sample: str) -> None:
    """Colours, fonts and page numbers are the renderer's business. A value
    like that in the file would mean the file had been written for one
    client, and the next one would inherit its choices."""
    blob = json.dumps(_results(sample))
    for banned in ('"colour":', '"color":', '"font"', '"page_number"', '"hex"'):
        assert banned not in blob, (
            f"{banned} appears in results.json; presentation belongs to whatever is "
            "drawing, not to the result"
        )
