"""AT157 — the embedded seed records, as records rather than as prose.

Every row has to arrive with a local identifier nothing else shares, a
route back to where it came from, and its gaps named rather than left
blank. The failure this guards against is the comfortable one: a record
that reads as complete because the fields it lacks were never declared.
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest
import yaml

from openbiota.extension import synbiotic as SYN

REPO = Path(__file__).resolve().parent.parent
SEEDS = REPO / "extension" / "synbiotic_seeds.yaml"
REGISTER = REPO / "extension" / "source_register.yaml"


@pytest.fixture(scope="module")
def records() -> list[dict]:
    return SYN.seeds()


def test_every_record_has_a_local_identifier_nothing_else_shares(records) -> None:
    ids = [r["syn_id"] for r in records]
    duplicated = [i for i, n in collections.Counter(ids).items() if n > 1]
    assert not duplicated, f"duplicate seed ids: {duplicated}"
    assert all(i.startswith("SYN") for i in ids)
    assert len(ids) == 17


@pytest.fixture(scope="module")
def register() -> dict[str, dict]:
    """The synbiotic sources, keyed by the id the seed rows cite."""
    raw = yaml.safe_load(REGISTER.read_text(encoding="utf-8"))
    return {str(entry["id"]): entry for entry in raw["synbiotic"]}


def test_every_record_names_where_it_came_from(records, register) -> None:
    """A claim whose source cannot be reached is not evidence."""
    for record in records:
        ref = record.get("source_ref")
        assert ref, f"{record['syn_id']} has no source_ref"
        assert ref in register, f"{record['syn_id']} cites unknown source {ref}"


def test_every_cited_source_carries_a_locator(records, register) -> None:
    """Provenance lives in the register, so the register must hold a way in.

    The seed row cites an id; that id has to resolve to a URL, a DOI or an
    accession, or the citation is a label rather than a route back.
    """
    for record in records:
        entry = register[str(record["source_ref"])]
        locators = (entry.get("urls") or []) + (entry.get("dois") or []) + (
            entry.get("accessions") or []
        )
        assert locators, (
            f"{record['syn_id']} cites {record['source_ref']}, which has nothing to follow"
        )


def test_an_incomplete_record_says_which_fields_are_open(records) -> None:
    """Silence about a missing field reads as the field not mattering."""
    for record in records:
        gaps = record.get("incomplete_fields")
        if gaps:
            assert all(isinstance(g, str) and g for g in gaps), record["syn_id"]


def test_a_record_that_cannot_show_synergy_says_which_arm_is_missing(records) -> None:
    """Eight of these never ran a component-only arm.

    Without that arm the trial cannot show synergy in either direction, and
    the reason has to be attached to the record rather than inferred.
    """
    missing_arm = [r for r in records if r.get("arms_missing")]
    assert len(missing_arm) == 8
    for record in missing_arm:
        assert record.get("demonstrated_synergy") in (False, None), (
            f"{record['syn_id']} claims synergy without the arm to show it"
        )


def test_the_negative_results_are_kept(records) -> None:
    """Four of these found the combination no better than the strain alone."""
    negative = [r for r in records if r.get("demonstrated_synergy") is False]
    assert len(negative) >= 4
    for record in negative:
        assert record.get("outcome"), f"{record['syn_id']} is negative with no outcome recorded"


def test_every_record_declares_its_evidence_setting(records) -> None:
    """A cell experiment and a randomised trial must not be read alike."""
    for record in records:
        assert record.get("setting"), f"{record['syn_id']} has no evidence setting"


def test_every_record_declares_what_it_may_be_lent_to(records) -> None:
    """The transfer key is what stops one strain's result becoming another's."""
    for record in records:
        assert record.get("evidence_transfer_key"), record["syn_id"]


def test_the_records_load_from_the_shipped_file_alone() -> None:
    """No research-note path, no network, no side file: the yaml is the record."""
    assert SEEDS.is_file()
    raw = yaml.safe_load(SEEDS.read_text(encoding="utf-8"))
    parsed = raw["records"] if isinstance(raw, dict) and "records" in raw else raw
    assert len(parsed) == len(SYN.seeds())


def test_the_summary_counts_agree_with_the_records(records) -> None:
    """A summary that drifts from its rows is worse than no summary."""
    summary = SYN.seed_summary()
    assert summary["n_records"] == len(records)
    assert summary["n_demonstrated_synergy"] == sum(
        1 for r in records if r.get("demonstrated_synergy") is True
    )
    assert summary["n_no_synergy_demonstrated"] == sum(
        1 for r in records if r.get("demonstrated_synergy") is False
    )
    assert summary["n_missing_a_component_arm"] == sum(
        1 for r in records if r.get("arms_missing")
    )
