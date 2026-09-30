"""A10 organism explorer — BUILD_SPEC_v0.8.3 sections 7.1/16.4, AT049/AT050/AT053.

The rule every test here circles: looking a name up is not the same as
finding the organism. A reader searching "Akkermansia muciniphila" must get a
straight answer whether it is present, absent or not something this assay can
speak to, and those must never be confused with each other.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import explorer as EXP
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parents[1]
SAMPLES = ["SAMPLE2_A02", "SAMPLE4_A04", "SAMPLE6_A06", "SAMPLE1_A01", "SAMPLE3_A03"]


def _results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


def _catalogue(sample: str = "SAMPLE2_A02") -> EXP.Catalogue:
    return EXP.Catalogue.from_inventory(_results(sample)["organism_inventory"]["organisms"])


def _tiny() -> EXP.Catalogue:
    return EXP.Catalogue.from_inventory([
        {"species": "Phocaeicola_vulgatus", "percent": 12.0, "genus": "Phocaeicola",
         "gtdb_genus": "Phocaeicola", "formerly": "Bacteroides_vulgatus", "percentile": 93.4},
        {"species": "Mediterraneibacter_gnavus", "percent": 0.5, "genus": "Mediterraneibacter",
         "gtdb_genus": "Mediterraneibacter", "formerly": "Ruminococcus_gnavus"},
        {"species": "Faecalibacterium_prausnitzii", "percent": 8.0,
         "genus": "Faecalibacterium", "gtdb_genus": "Faecalibacterium"},
        {"species": "GGB9480_SGB14874", "percent": 0.6},
    ])


# --------------------------------------------------------------------------- #
# AT049 — every fixture label resolves to an explicit state
# --------------------------------------------------------------------------- #


def test_at049_all_332_fixture_labels_resolve_and_none_is_dropped() -> None:
    fixture = EXP.load_fixture()
    assert len(fixture["resolution_fixture"]) == 332
    audit = EXP.resolution_audit(_catalogue(), fixture["resolution_fixture"])
    assert audit["n_queried"] == 332
    assert audit["unresolvable"] == [], audit["unresolvable"][:10]
    assert sum(audit["by_resolution_state"].values()) == 332
    for result in audit["results"]:
        assert result["resolution_state"] in EXP.RESOLUTION_STATES
        assert result["detection_state"] in EXP.DETECTION_STATES


def test_at049_the_fixture_is_not_a_list_of_required_detections() -> None:
    """Most of these names are legitimately absent from any one sample."""
    audit = EXP.resolution_audit(_catalogue(), EXP.load_fixture()["resolution_fixture"])
    assert audit["n_detected"] < audit["n_queried"], (
        "a fixture where everything is detected is not testing resolution"
    )
    absent = [r for r in audit["results"] if r["detection_state"] == "no_supported_detection"]
    assert absent, "some fixture labels must be absent for the negative path to be tested"
    for row in absent[:20]:
        assert row["resolution_state"] != "unresolvable"


def test_gtdb_split_suffixes_and_provisional_ids_are_preserved() -> None:
    """`Ruminococcus_B gnavus` and `sp900756035` are identifiers, not errors."""
    catalogue = _catalogue()
    split = catalogue.resolve("Ruminococcus_B gnavus")
    assert split.resolution_state in EXP.RESOLUTION_STATES
    assert split.query == "Ruminococcus_B gnavus"

    # A provisional identifier the catalogue does not carry resolves as a
    # cluster, not as a failure. Checked against a fixed catalogue: in a real
    # sample such a name may legitimately be an inventory alias.
    tiny = _tiny()
    provisional = tiny.resolve("Anaerostipes sp900756035")
    assert provisional.resolution_state == "cluster_correspondence"
    assert provisional.rank == "cluster"
    assert "not a name awaiting resolution" in (provisional.note or "")
    assert tiny.resolve("UBA9502 sp900538475").resolution_state == "cluster_correspondence"


# --------------------------------------------------------------------------- #
# resolving is not detecting
# --------------------------------------------------------------------------- #


def test_a_resolution_never_implies_a_detection() -> None:
    catalogue = _tiny()
    # A recognisable name this sample does not have.
    absent = catalogue.resolve("Bacteroides fragilis")
    assert absent.resolution_state == "not_in_catalogue"
    assert absent.detection_state == "no_supported_detection"
    assert absent.percent is None
    assert absent.note and "detected in this sample" in absent.note
    # Every record carries the reminder in its payload.
    assert "not a detection" in absent.to_json()["reminder"]


def test_an_older_name_resolves_without_becoming_a_new_finding() -> None:
    """*Ruminococcus gnavus* is *Mediterraneibacter gnavus*: one organism."""
    catalogue = _tiny()
    old = catalogue.resolve("Ruminococcus gnavus")
    assert old.resolution_state == "accepted_synonym"
    assert old.accepted_name == "Mediterraneibacter gnavus"
    assert old.detection_state == "supported_detection"
    new = catalogue.resolve("Mediterraneibacter gnavus")
    # The same organism and the same abundance under either name; not two.
    assert new.percent == old.percent
    assert new.accepted_name == old.accepted_name
    assert old.matched_via and "earlier name" in old.matched_via


def test_at050_an_alias_is_not_counted_as_a_second_organism() -> None:
    catalogue = _tiny()
    total = sum(float(r.get("percent") or 0) for r in catalogue.rows)
    both = catalogue.resolve("Bacteroides vulgatus"), catalogue.resolve("Phocaeicola vulgatus")
    assert both[0].accepted_name == both[1].accepted_name
    # The catalogue holds four rows however many names point at them.
    assert len(catalogue.rows) == 4
    assert total == pytest.approx(21.1)


def test_an_unresolvable_string_says_so() -> None:
    catalogue = _tiny()
    junk = catalogue.resolve("!!!")
    assert junk.resolution_state == "unresolvable"
    assert junk.detection_state == "not_assessable"


# --------------------------------------------------------------------------- #
# AT053 — a species label is not a strain, a toxin or a product
# --------------------------------------------------------------------------- #


def test_at053_a_species_detection_claims_nothing_about_strain_or_toxin() -> None:
    view = EXP.build(_results("SAMPLE2_A02"))
    joined = " ".join(view["limitations"])
    assert "does not establish a commercial probiotic strain" in joined
    assert "toxin" in joined
    for sentinel in view["sentinels"]:
        assert "strain" not in sentinel["resolution_state"]
        assert sentinel["rank"] in {"species", "genus", "cluster"}


# --------------------------------------------------------------------------- #
# rollups disclose their membership
# --------------------------------------------------------------------------- #


def test_a_genus_rollup_lists_its_members_and_says_it_is_a_sum() -> None:
    rows = EXP.rollup(_tiny().rows)
    by_genus = {r["genus"]: r for r in rows}
    phoc = by_genus["Phocaeicola"]
    assert phoc["percent"] == pytest.approx(12.0)
    assert phoc["members"] == ["Phocaeicola vulgatus"]
    assert "not an additional organism" in phoc["note"]
    assert "must not be charted beside them" in phoc["note"]


def test_organisms_with_no_genus_are_shown_not_dropped() -> None:
    rows = EXP.rollup(_tiny().rows)
    unassigned = [r for r in rows if r["genus"] == "(no genus assigned)"]
    assert unassigned, "an organism with no genus must still appear"
    assert unassigned[0]["percent"] == pytest.approx(0.6)
    # And the rollup accounts for everything.
    assert sum(r["percent"] for r in rows) == pytest.approx(
        sum(float(x.get("percent") or 0) for x in _tiny().rows)
    )


def test_a_genus_lookup_returns_a_rollup_with_its_species() -> None:
    catalogue = _catalogue()
    genus = catalogue.resolve("Bifidobacterium")
    assert genus.resolution_state == "rank_rollup"
    assert genus.rank == "genus"
    assert genus.candidates, "a rollup must disclose the species it summed"
    assert genus.percent is not None and genus.percent > 0
    assert "Genus total over" in (genus.note or "")


# --------------------------------------------------------------------------- #
# role composition keeps unknown visible
# --------------------------------------------------------------------------- #


def test_unassigned_mass_is_not_removed_from_the_denominator() -> None:
    organisms = [
        {"species": "A", "percent": 40.0},
        {"species": "B", "percent": 30.0},
        {"species": "C", "percent": 30.0},
    ]
    verdicts = [{"species": "A", "class": "beneficial"}]
    comp = EXP.role_composition(organisms, verdicts)
    assert comp["unassigned_percent"] == pytest.approx(60.0)
    assert comp["unassigned_share"] == pytest.approx(0.6)
    assert comp["total_percent"] == pytest.approx(100.0)
    # The one classified organism is not rescaled to look like the whole.
    assert comp["classes"][0]["share_of_named"] == pytest.approx(0.4)
    assert "rather than removed from the denominator" in comp["note"]


def test_unknown_is_an_evidence_statement_not_a_harm() -> None:
    comp = EXP.role_composition(
        [{"species": "A", "percent": 10.0}], [{"species": "A", "class": "unknown"}]
    )
    assert comp["classes"][0]["role"] == "unknown"
    assert "not about harm" in comp["note"]


# --------------------------------------------------------------------------- #
# the sentinels and the whole view
# --------------------------------------------------------------------------- #


def test_the_34_sentinels_are_a_closed_expanded_list() -> None:
    sentinels = EXP.load_fixture()["sentinels"]
    assert len(sentinels) == 34
    # Abbreviated genus names are expanded, as section 7.1 requires.
    assert not any(s.startswith(("B. ", "F. ", "A. ", "R. ")) for s in sentinels), sentinels
    for required in ("Akkermansia muciniphila", "Faecalibacterium prausnitzii",
                     "Methanobrevibacter smithii", "Oxalobacter formigenes",
                     "Bifidobacterium adolescentis", "Blautia wexlerae"):
        assert required in sentinels, required


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_explorer_builds_for_every_sample(sample: str) -> None:
    view = EXP.build(_results(sample))
    assert view["feature_id"] == "A10"
    assert view["n_organisms"] > 50
    assert len(view["sentinels"]) == 34
    # Every sentinel answers, including the ones that are absent.
    for sentinel in view["sentinels"]:
        assert sentinel["resolution_state"] in EXP.RESOLUTION_STATES
        assert sentinel["detection_state"] in EXP.DETECTION_STATES
    # Some are present and some are not: a view where everything is found is
    # not a search.
    states = {s["detection_state"] for s in view["sentinels"]}
    assert "supported_detection" in states
    assert "no_supported_detection" in states

    for metric in EXP.metrics(view):
        assert metric.feature_id == "A10"
        assert metric.direction == "descriptive"
        metric.to_json()


def test_a_second_lane_detection_is_not_reported_as_absent() -> None:
    """The lookup must agree with the catalogue page a reader just read.

    An organism the primary lane scores at zero can be quantified by the
    second lane, and the catalogue lists it with that abundance. Calling the
    same organism "not found" two pages later is the kind of contradiction
    that makes a reader distrust every other number.
    """
    catalogue = EXP.Catalogue.from_inventory([
        {"species": "Faecalibacterium_prausnitzii", "percent": 0.0,
         "secondary_percent": 0.0049, "genus": "Faecalibacterium"},
    ])
    found = catalogue.resolve("Faecalibacterium prausnitzii")
    assert found.detection_state == "supported_detection"
    assert found.percent == pytest.approx(0.0049)
    assert found.note and "second detection lane" in found.note


def test_a_dual_name_resolves_on_either_half() -> None:
    """"Agathobacter rectalis/Eubacterium rectale" is one organism."""
    catalogue = EXP.Catalogue.from_inventory([
        {"species": "Agathobacter_rectalis", "percent": 0.1, "genus": "Agathobacter"},
    ])
    for query in ("Agathobacter rectalis/Eubacterium rectale",
                  "Eubacterium rectale/Agathobacter rectalis"):
        hit = catalogue.resolve(query)
        assert hit.detection_state == "supported_detection", query
        assert hit.query == query, "the reader's own wording is preserved"
        assert hit.matched_via and "matched on" in hit.matched_via
    # A dual name of two organisms neither of which is here still says so.
    miss = catalogue.resolve("Nothing here/Also absent")
    assert miss.detection_state == "no_supported_detection"


def test_the_organism_count_is_not_presented_as_quality() -> None:
    """Section 12.4: a large count is not the definition of a good result."""
    view = EXP.build(_results("SAMPLE2_A02"))
    count = next(
        m for m in EXP.metrics(view) if m.metric_id == "ext083.explorer.organisms_listed"
    )
    assert "not a measure of quality" in " ".join(count.limitations)
    assert count.direction == "descriptive"
