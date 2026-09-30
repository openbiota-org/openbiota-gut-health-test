"""No section of the report may contradict another.

Three real contradictions motivated this file, all in shipped reports:

* the pathogen page called *Mediterraneibacter gnavus* "a large signal, 0.34%
  of analysed DNA" while every other section of the same report said the
  species was absent;
* a biofilm card printed "Mediterraneibacter gnavus 0%, not detected" for a
  sample whose organism list gave the same organism 6.9% of the community;
* *Clostridium innocuum*, a vancomycin-resistant opportunist 372% above the
  level of the adults who carry it, was printed in the full organism table
  and left out of the page headed "organisms that need attention".

Each check below reads only ``results.json`` - the single source of truth the
PDF is drawn from - so a disagreement is caught at its source.
"""

from __future__ import annotations

import pytest

from openbiota import inventory as inventory_mod
from openbiota import organisms as org

from ._data import results as _results

SAMPLES = ("SAMPLE1_A01", "SAMPLE2_A02", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE5_A05", "SAMPLE6_A06")

#: Pathogen statuses that name an organism as present.
NAMED = frozenset({"supported_sequence", "marker_signal"})


def _load(sample: str):
    blob = _results(sample)
    inv = inventory_mod.from_json(blob.get("organism_inventory"))
    if inv is None or not len(inv):
        pytest.skip(f"{sample}: no organism inventory")
    return blob, inv


def _found(o) -> bool:
    return bool(o) and ((o.in_primary and o.percent > 0)
                        or (o.status == "supported" and (o.best_percent or 0.0) > 0))


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_pathogen_is_named_that_the_inventory_says_is_absent(sample: str) -> None:
    """A bacterial call the pooled inventory did not find must either be
    reported as relatives' shared sequence, or say why it stands alone."""
    blob, inv = _load(sample)
    pathogens = blob.get("pathogens")
    if not isinstance(pathogens, dict) or not pathogens.get("results"):
        pytest.skip(f"{sample}: no pathogen screen")
    rec_block = pathogens.get("inventory_reconciliation") or {}
    if rec_block.get("status") != "completed":
        pytest.skip(f"{sample}: results predate the pathogen reconciliation")
    from openbiota.pathogens import reconcile

    for rec in pathogens["results"]:
        if rec.get("group") != "bacteria" or rec.get("sequence_status") not in NAMED:
            continue
        if rec.get("species_resolution") not in (None, "resolved", "resolved_by_marker"):
            continue
        agreement = rec.get("inventory_agreement")
        name = rec.get("inventory_species") or rec.get("display_name")
        assert agreement in ("found", "not_in_inventory"), (sample, name, agreement)
        if agreement == "found":
            # the organism the screen names is in the composition, and the card
            # can quote its share
            o = inv.get(str(name).replace(" ", "_")) or inv.get(str(name))
            if o is None:
                for cand in inv.organisms:
                    if reconcile._norm(cand.display) == reconcile._norm(str(name)):
                        o = cand
                        break
            assert _found(o), (sample, name)
        else:
            # it stands on the screen's own alignment: coverage shaped like a
            # population, and the row says so
            assert "not_in_organism_inventory" in (rec.get("reason_codes") or []), (sample, name)
            even = rec.get("inventory_evenness")
            assert even is None or even >= 0.25, (sample, name, even)


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_withdrawn_pathogen_call_states_no_amount_of_the_organism(sample: str) -> None:
    blob, _inv = _load(sample)
    pathogens = blob.get("pathogens")
    if not isinstance(pathogens, dict) or (pathogens.get("inventory_reconciliation") or {}).get("status") != "completed":
        pytest.skip(f"{sample}: no reconciled pathogen screen")
    from openbiota.pdfpathogens import amount_of, tier_of

    for rec in pathogens["results"]:
        if rec.get("inventory_agreement") != "shared_sequence_from_relatives":
            continue
        assert rec["sequence_status"] == "ambiguous_signal" and rec["counts_as_pathogen"] is False
        assert tier_of(rec) == "shared_sequence"
        amount = amount_of(rec)
        # no percentage of the sample is claimed for an organism that is not there
        assert amount["percent"] is None and amount.get("shared") is True
        assert "of analysed DNA" not in (amount["share"] or "")
        # the fragments are still shown, and the relatives named
        assert amount["fragments"] == rec["unique_supporting_fragments"]
        assert "did not find" in rec["plain_statement"]


@pytest.mark.parametrize("sample", SAMPLES)
def test_the_headline_pathogen_count_matches_the_records(sample: str) -> None:
    blob, _inv = _load(sample)
    pathogens = blob.get("pathogens")
    if not isinstance(pathogens, dict) or not pathogens.get("results"):
        pytest.skip(f"{sample}: no pathogen screen")
    named = [r for r in pathogens["results"]
             if r.get("sequence_status") in NAMED and r.get("display_status") != "technical_artifact"]
    expected = sum(1 for r in named if r.get("interpretation_class") != "background_or_decoy")
    assert pathogens["counts"]["pathogen_count"] == expected
    assert pathogens["counts"]["supported"] == len(named)


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_named_fungal_pathogen_is_one_the_mycobiome_also_found(sample: str) -> None:
    """The fungal half of the screen has a second opinion of its own: the
    mycobiome section maps its own competitive index over 961 representative
    genomes. A fungus named on the pathogen page that section did not see
    would be the same contradiction the bacterial reconciliation removes."""
    blob, _inv = _load(sample)
    pathogens = blob.get("pathogens")
    myco = blob.get("mycobiome")
    if not isinstance(pathogens, dict) or not isinstance(myco, dict):
        pytest.skip(f"{sample}: no pathogen screen or no mycobiome section")
    names: set[str] = set()
    for taxon in ((myco.get("genome_lane") or {}).get("taxa") or []):
        if isinstance(taxon, dict) and taxon.get("name"):
            names.add(str(taxon["name"]).lower())
    for taxon in ((myco.get("marker_lane") or {}).get("fungal") or []):
        if isinstance(taxon, dict):
            names.add(str(taxon.get("name") or taxon.get("species") or "").lower())
        elif isinstance(taxon, str):
            names.add(taxon.lower())
    named = [r for r in pathogens["results"]
             if r.get("group") == "fungi" and r.get("sequence_status") in NAMED]
    if not named:
        pytest.skip(f"{sample}: no fungus named on the pathogen page")
    if not names:
        pytest.skip(f"{sample}: the mycobiome lanes produced no taxon list")
    for rec in named:
        shown = str(rec.get("display_name") or "").lower()
        assert any(shown in n or n in shown for n in names if n), (sample, rec.get("display_name"))


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_biofilm_card_says_when_the_pooled_methods_disagree(sample: str) -> None:
    """The card's value must be the reference catalogue's, and where the pooled
    inventory found the organism anyway, the card has to say so."""
    blob, _inv = _load(sample)
    biofilm = blob.get("biofilm")
    if not isinstance(biofilm, dict) or not biofilm.get("cards"):
        pytest.skip(f"{sample}: no biofilm section")
    from openbiota.biofilm import engine, reference

    checked = any("inventory_note" in info
                  for card in biofilm["cards"]
                  for info in (card.get("feature_detail") or {}).values()
                  if isinstance(info, dict))
    if not checked:
        pytest.skip(f"{sample}: results predate the biofilm/inventory cross-check")
    seen = False
    for card in biofilm["cards"]:
        for name, info in (card.get("feature_detail") or {}).items():
            if not isinstance(info, dict) or info.get("detected"):
                continue
            leaves = reference.FEATURE_LEAVES.get(name)
            if not leaves:
                continue
            note = engine.inventory_note(blob, leaves)
            if note is None:
                continue
            seen = True
            assert info.get("inventory_note"), (sample, name)
            assert "section 6" in info["inventory_note"]
    if not seen:
        pytest.skip(f"{sample}: reference catalogue and pooled methods agree on every biofilm feature")


@pytest.mark.parametrize("sample", SAMPLES)
def test_every_opportunist_well_above_its_carriers_is_on_the_attention_page(sample: str) -> None:
    """An opportunist several times above the people who carry it, above the
    95th percentile of all reference adults and seen by more than one method,
    is a reading to act on however few adults carry it."""
    _blob, inv = _load(sample)
    from openbiota import pdfattention

    flagged = {v.organism.species for v in org.verdicts(list(inv.organisms)) if v.flagged}
    for o in inv.organisms:
        v = org.verdict(o)
        if v.cls != org.OPPORTUNIST or not (o.in_primary and o.percent > 0):
            continue
        if o.percentile is None or o.deviation_percent is None or o.level_percentile is None:
            continue
        if len(o.methods or ()) < org.UNCOMMON_MIN_METHODS:
            continue
        if (o.percentile >= org.UNCOMMON_POPULATION_AT
                and o.deviation_percent >= org.UNCOMMON_DEVIATION_AT
                and o.level_percentile >= org.UNCOMMON_CARRIER_AT):
            assert o.species in flagged, (sample, o.species, o.percent, o.deviation_percent)
            assert v.is_issue, (sample, o.species)
    # and the attention list contains every flagged organism, under one name
    items = pdfattention.build(inv, None)
    listed = {a.species for a in items}
    assert flagged <= listed, sorted(flagged - listed)
    assert len(listed) == len(items), "same organism listed twice"


@pytest.mark.parametrize("sample", SAMPLES)
def test_a_level_and_its_deviation_never_point_opposite_ways(sample: str) -> None:
    _blob, inv = _load(sample)
    for o in inv.organisms:
        lp, dev = o.level_percentile, o.deviation_percent
        if lp is None or dev is None or abs(dev) < 5:
            continue
        if lp == 50.0:
            # a midrank of exactly 50 is the tie case: the reading sits on the
            # median itself, and the printed deviation is inside the width of
            # that tie. Not a contradiction, and the renderer prints no arrow.
            continue
        assert (dev > 0) == (lp > 50.0), (sample, o.species, lp, dev)


@pytest.mark.parametrize("sample", SAMPLES)
def test_nothing_listed_as_missing_was_found_by_any_method(sample: str) -> None:
    _blob, inv = _load(sample)
    from openbiota import pdfattention

    for a in pdfattention.build(inv, None):
        if a.group != "missing":
            continue
        o = inv.get(a.species)
        assert not _found(o), (sample, a.species)


def test_a_catalogues_non_detection_carries_the_pooled_reading() -> None:
    """A disease-pattern row may say "not detected" only in its own catalogue's
    voice: where every method together found the organism, the row prints that
    share, so the table cannot contradict the organism list."""
    from openbiota import pdfcontext

    blob = {"primary_lane": "jan26", "unclassified_percent": 1.0, "organisms": [
        {"species": "Blautia_A_wexlerae", "gtdb": "Blautia_A wexlerae", "percent": 0.5875,
         "in_primary": True, "genus": "Blautia", "status": "supported", "detected_by": ["jan26", "globdb"],
         "aliases": ["Blautia_wexlerae"]},
        {"species": "Dorea_hominis", "gtdb": "Dorea_D hominis", "percent": 0.102, "in_primary": True,
         "genus": "Dorea", "status": "supported", "detected_by": ["jan26"],
         "formerly_listed_as": ["Mediterraneibacter_gnavus"]},
        {"species": "Roseburia_inulinivorans", "percent": 0.0, "in_primary": False, "genus": "Roseburia",
         "status": "provisional", "detected_by": ["rescue"], "share_basis": "member",
         "counted_within": "Agathobacter rectalis", "secondary_percent": 0.06},
    ]}
    inv = inventory_mod.from_json(blob)
    assert inv is not None
    pdfcontext.set_inventory(inv)
    try:
        # the reference catalogue's own name for the organism resolves to the current record
        assert pdfcontext.pooled_note("Blautia_wexlerae", catalogue_value=None) == "0.588% pooled"
        # an older catalogue's label for a population that turned out to be a
        # different species answers nothing: the bin the 2023 catalogue called
        # Ruminococcus gnavus is Dorea hominis, and quoting its share here
        # would assert that gnavus is present
        assert pdfcontext.pooled_note("Ruminococcus_gnavus", catalogue_value=None) is None
        # where that catalogue did measure it, nothing is added
        assert pdfcontext.pooled_note("Blautia_wexlerae", catalogue_value=0.3) is None
        # a population counted inside a relative's share has no share to quote
        assert pdfcontext.pooled_note("Roseburia_inulinivorans", catalogue_value=None) is None
        # and an organism no method found says nothing
        assert pdfcontext.pooled_note("Escherichia_coli", catalogue_value=None) is None
    finally:
        pdfcontext.set_inventory(None)


def test_the_pooled_reading_helper_is_safe_without_an_inventory() -> None:
    from openbiota import pdfcontext

    pdfcontext.set_inventory(None)
    assert pdfcontext.pooled_note("Blautia_wexlerae", catalogue_value=None) is None
    assert pdfcontext.pooled_share("Blautia_wexlerae") is None
    assert pdfcontext.organism("") is None
