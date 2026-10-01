"""The organism verdicts: class, flag, weight, and what reaches page one.

Every test here exists because the opposite once happened in a real report:
a probiotic flagged red for being present, a 13% opportunist bloom outranked
by a trace organism, an overgrowth cut from page one by a fixed category
order, and a real organism printed as an error code.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from openbiota import gtdb, inventory
from openbiota import organisms as org

RESULTS = sorted(Path("results").glob("*/results.json"))


def _o(species: str, percent: float, percentile: float | None = None, **kw) -> inventory.Organism:
    return inventory.Organism(species=species, percent=percent, percentile=percentile, **kw)


# --------------------------------------------------------------------------- #
# class
# --------------------------------------------------------------------------- #

def test_every_class_is_one_of_four() -> None:
    for name in ("Faecalibacterium_prausnitzii", "Ruminococcus_gnavus", "Blautia_wexlerae",
                 "Zzz_nothing_knows_this"):
        assert org.verdict(_o(name, 1.0)).cls in org.CLASSES


def test_the_classic_beneficial_and_opportunist_land_where_expected() -> None:
    assert org.verdict(_o("Faecalibacterium_prausnitzii", 1.0)).cls == org.BENEFICIAL
    assert org.verdict(_o("Ruminococcus_gnavus", 1.0)).cls == org.OPPORTUNIST
    assert org.verdict(_o("Escherichia_coli", 1.0)).cls == org.OPPORTUNIST


def test_a_probiotic_is_beneficial_not_conditional() -> None:
    """The registry files probiotics as context-dependent because they are
    transient. As a class that misleads; the genus table's verdict stands."""
    for name in ("Lactobacillus_plantarum", "Bifidobacterium_animalis", "Lactobacillus_rhamnosus"):
        assert org.verdict(_o(name, 1.0)).cls == org.BENEFICIAL, name


def test_an_unknown_genus_falls_back_to_its_gtdb_family() -> None:
    if not gtdb.available():
        pytest.skip("no GTDB mapping installed")
    # SGB63306 is Gemmiger in GTDB; Gemmiger is in the genus table.
    v = org.verdict(_o("GGB45596_SGB63306", 8.0, sgb="SGB63306", unnamed=True,
                       gtdb="Gemmiger sp937890665", gtdb_genus="Gemmiger"))
    assert v.cls == org.BENEFICIAL
    assert v.basis == "genus"


# --------------------------------------------------------------------------- #
# flag
# --------------------------------------------------------------------------- #

def test_a_probiotic_at_the_ceiling_is_never_flagged() -> None:
    """A supplement being taken is not an overgrowth."""
    for name in ("Lactobacillus_plantarum", "Bifidobacterium_animalis"):
        v = org.verdict(_o(name, 2.0, percentile=99.9))
        assert not v.flagged, (name, v.flag_reason)


def test_an_opportunist_high_is_an_issue() -> None:
    v = org.verdict(_o("Ruminococcus_gnavus", 13.0, percentile=98.0))
    assert v.flag == "high" and v.is_issue


def test_a_beneficial_organism_low_is_an_issue() -> None:
    v = org.verdict(_o("Faecalibacterium_prausnitzii", 0.01, percentile=2.0))
    assert v.flag == "low" and v.is_issue


def test_a_beneficial_organism_high_is_at_most_a_watch() -> None:
    v = org.verdict(_o("Akkermansia_muciniphila", 5.0, percentile=97.0))
    assert v.flag_level <= org.FLAG_WATCH, "abundant beneficial must never sit with opportunist blooms"


def test_no_percentile_means_no_flag() -> None:
    v = org.verdict(_o("Ruminococcus_gnavus", 13.0, percentile=None))
    assert not v.flagged


def test_no_flag_reason_ever_says_above_100_percent() -> None:
    v = org.verdict(_o("Ruminococcus_gnavus", 1.0, percentile=100.0))
    assert "100%" not in v.flag_reason
    assert ">99th" in v.flag_reason


# --------------------------------------------------------------------------- #
# weight
# --------------------------------------------------------------------------- #

def test_a_bloom_outranks_a_trace_organism_further_out() -> None:
    """13% at the 98th percentile is the finding; 0.3% at the 99.6th is not."""
    bloom = org.verdict(_o("Ruminococcus_gnavus", 13.3, percentile=98.0))
    trace = org.verdict(_o("Blautia_sp_N6H1_15", 0.27, percentile=99.6))
    assert org.importance(bloom) > org.importance(trace)


def test_issues_are_ordered_by_weight() -> None:
    vs = [
        org.verdict(_o("Blautia_sp_N6H1_15", 0.27, percentile=99.6)),
        org.verdict(_o("Ruminococcus_gnavus", 13.3, percentile=98.0)),
        org.verdict(_o("Faecalibacterium_prausnitzii", 0.01, percentile=2.0)),
    ]
    ranked = org.issues(vs)
    assert ranked[0].organism.species == "Ruminococcus_gnavus"


def test_composition_sums_to_one_hundred() -> None:
    vs = [org.verdict(_o("Faecalibacterium_prausnitzii", 30.0)),
          org.verdict(_o("Ruminococcus_gnavus", 10.0)),
          org.verdict(_o("Zzz_nothing", 60.0))]
    assert sum(org.composition(vs).values()) == pytest.approx(100.0)


# --------------------------------------------------------------------------- #
# naming
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_unnamed_bins_print_a_real_name_where_gtdb_has_one() -> None:
    if not gtdb.available():
        pytest.skip("no GTDB mapping installed")
    for path in RESULTS:
        inv = inventory.from_json(json.loads(path.read_text()).get("organism_inventory"))
        if inv is None:
            continue
        unnamed = inv.unnamed
        if not unnamed:
            continue
        resolved = [o for o in unnamed if o.gtdb]
        # Nearly every bin maps; three of 164 did not on the test sample.
        assert len(resolved) >= 0.9 * len(unnamed), path
        for o in resolved:
            assert not re.match(r"^GGB\d+", o.display), (
                f"{path}: {o.species} still prints as a bin identifier")


def test_a_named_species_keeps_its_name_and_carries_gtdb_as_alias() -> None:
    if not gtdb.available():
        pytest.skip("no GTDB mapping installed")
    hit = gtdb.lookup("SGB4584")  # Ruminococcus_B gnavus in GTDB
    assert hit is not None
    o = _o("Mediterraneibacter_gnavus", 13.0, sgb="SGB4584", gtdb=hit.display)
    assert o.display == "Mediterraneibacter gnavus"
    assert "gnavus" in (o.gtdb or "")


# --------------------------------------------------------------------------- #
# every bar is the reader's own position
# --------------------------------------------------------------------------- #

def test_no_bar_is_ever_fed_the_reference_carriage_rate() -> None:
    """A percentile bar shows where *this sample* sits. It once showed the
    share of reference adults carrying an organism the sample did not have,
    so a missing organism drew its marker at 88. The carriage rate belongs
    in words, never on the bar."""
    offenders = []
    for path in sorted(Path("openbiota").glob("pdf*.py")):
        text = path.read_text()
        for m in re.finditer(r"PercentileBar\((?:[^()]|\([^()]*\))*\)", text, re.S):
            call = m.group(0)
            if re.search(r"percentile\s*=\s*[^,)]*prevalence", call):
                line = text[:m.start()].count("\n") + 1
                offenders.append(f"{path.name}:{line}")
        # Also the pattern of pre-scaling prevalence into a percentile.
        for m in re.finditer(r"100(?:\.0)?\s*\*\s*\w+\.prevalence", text):
            line = text[:m.start()].count("\n") + 1
            nearby = text[m.start():m.start() + 400]
            if "PercentileBar" in nearby:
                offenders.append(f"{path.name}:{line} (prevalence scaled toward a bar)")
    assert not offenders, offenders


# --------------------------------------------------------------------------- #
# page one
# --------------------------------------------------------------------------- #

def test_an_acute_overgrowth_leads_the_organism_headlines() -> None:
    """A fixed Missing/Depleted/Expanded order once cut every overgrowth
    from page one. Ranked by weight, the bloom leads."""
    from openbiota.pdfsummary import organism_headlines

    inv = inventory.Inventory(organisms=[
        _o("Ruminococcus_gnavus", 13.3, percentile=98.0),
        _o("Faecalibacterium_prausnitzii", 0.01, percentile=2.0),
        _o("Blautia_wexlerae", 3.0, percentile=50.0),
    ])
    items = organism_headlines(None, inv, limit=3)
    assert items, "an overgrowth and a depletion should each produce a headline"
    first = re.sub(r"<[^>]+>", "", items[0][1])
    assert first.startswith("Overgrown"), first
    assert "gnavus" in first


def test_every_organism_headline_lands_on_the_attention_section() -> None:
    """The overgrown headline once linked to the full catalogue - every
    organism detected - instead of the section that lists what is overgrown,
    missing and low. The reader clicked on a problem; land them on it."""
    from openbiota.pdflinks import section_dest
    from openbiota.pdfsummary import organism_headlines

    inv = inventory.Inventory(organisms=[
        _o("Ruminococcus_gnavus", 13.3, percentile=98.0),
        _o("Faecalibacterium_prausnitzii", 0.01, percentile=2.0),
    ])
    for _, text, dest in organism_headlines(None, inv, limit=3):
        assert dest == section_dest("organisms"), (re.sub(r"<[^>]+>", "", text), dest)


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_page_one_headline_counts_only_issues_not_watches() -> None:
    from openbiota.pdfsummary import organism_headlines

    for path in RESULTS:
        inv = inventory.from_json(json.loads(path.read_text()).get("organism_inventory"))
        if inv is None:
            continue
        for _, text, _ in organism_headlines(None, inv, limit=3):
            plain = re.sub(r"<[^>]+>", "", text)
            m = re.search(r"and (\d+) more at a level of concern", plain)
            if m:
                vs = org.verdicts(list(inv.organisms))
                serious = [v for v in org.issues(vs)
                           if v.flag == "high" and v.cls in (org.OPPORTUNIST, org.CONDITIONAL)
                           and v.is_issue]
                assert int(m.group(1)) == len(serious) - 2, path


# --------------------------------------------------------------------------- #
# the attention rule, every branch, on synthetic organisms
# --------------------------------------------------------------------------- #

def _ranked(species: str, percent: float, carrier_pct: float, **kw) -> inventory.Organism:
    """An organism the scoring cohort ranked: share, rank among carriers, and
    the fields the rule reads, with sensible defaults for the rest."""
    base = {
        "species": species, "percent": percent, "in_primary": True, "status": "supported",
        "percentile": kw.pop("population_pct", carrier_pct), "carrier_percentile": carrier_pct,
        "reference_carriers": kw.pop("carriers", 500), "reference_percent": kw.pop("typical", 0.1),
        "reference_reading": kw.pop("reading", percent), "scoring_percent": percent,
        "percentile_source": "scoring cohort", "prevalence": kw.pop("prevalence", 0.6),
        "methods": kw.pop("methods", ("marker", "genome_sketch")),
    }
    base.update(kw)
    return inventory.Organism(**base)


def test_the_attention_rule_branch_by_branch() -> None:
    """`organisms._flag` plus `_expanded_uncommon` is the whole rule; the page is
    its image. Each row: what the organism is, where it sits, what the page says."""
    cases = [
        # opportunist at or above the 90th among carriers: overgrown
        (_ranked("Ruminococcus_gnavus", 2.0, 90.0), "high", org.FLAG_ISSUE),
        (_ranked("Ruminococcus_gnavus", 2.0, 89.9), "", org.FLAG_NONE),
        # a conditional resident: watch between the 90th and 97th, overgrown from the 97th
        (_ranked("Phocaeicola_vulgatus", 2.0, 92.0), "high", org.FLAG_WATCH),
        (_ranked("Phocaeicola_vulgatus", 2.0, 97.0), "high", org.FLAG_ISSUE),
        # a beneficial organism high is not a finding, however high
        (_ranked("Faecalibacterium_prausnitzii", 20.0, 99.9), "", org.FLAG_NONE),
        (_ranked("Lactobacillus_plantarum", 2.0, 100.0), "", org.FLAG_NONE),
        # a beneficial organism low: depleted, an issue from the 3rd percentile down
        (_ranked("Faecalibacterium_prausnitzii", 0.05, 10.0), "low", org.FLAG_WATCH),
        (_ranked("Faecalibacterium_prausnitzii", 0.01, 3.0), "low", org.FLAG_ISSUE),
        (_ranked("Faecalibacterium_prausnitzii", 0.2, 10.1), "", org.FLAG_NONE),
        # an opportunist low is nothing
        (_ranked("Ruminococcus_gnavus", 0.001, 1.0), "", org.FLAG_NONE),
        # the expanded-uncommon opportunist: top 5% of all adults, >=3x the typical carrier, upper quarter of carriers, >1 method
        (_ranked("Clostridium_innocuum", 0.04, 80.7, population_pct=95.5, typical=0.0183, reading=0.0863, prevalence=0.29),
         "high", org.FLAG_ISSUE),
        (_ranked("Clostridium_innocuum", 0.04, 80.7, population_pct=94.0, typical=0.0183, reading=0.0863, prevalence=0.29),
         "", org.FLAG_NONE),
        (_ranked("Clostridium_innocuum", 0.04, 80.7, population_pct=95.5, typical=0.0183, reading=0.0863, prevalence=0.29,
                 methods=("marker",)), "", org.FLAG_NONE),
        # an uncommon opportunist merely present: watch
        (_ranked("Enterococcus_faecium", 0.02, 66.0, prevalence=0.05), "uncommon", org.FLAG_WATCH),
        # no rank, no flag: not in the reference set, too few carriers, or the catalogue reads it differently
        (_ranked("Ruminococcus_gnavus", 5.0, 99.0, carriers=4), "", org.FLAG_NONE),
        (_ranked("Ruminococcus_gnavus", 5.0, 99.0, reference_conflict=True), "", org.FLAG_NONE),
        (inventory.Organism(species="Ruminococcus_gnavus", percent=5.0), "", org.FLAG_NONE),
    ]
    for o, flag, level in cases:
        v = org.verdict(o)
        assert (v.flag, v.flag_level) == (flag, level), (o.species, o.percent, o.carrier_percentile, v.flag, v.flag_level, v.flag_reason)


def test_a_population_counted_within_a_relative_is_judged_once() -> None:
    """The same population under two catalogues' names is judged as the
    relative; a distinct population the competition told apart, or one whose
    relative has no rank, is judged on its own lane rank."""
    same = _ranked("Phocaeicola_SPECIV4_34405", 0.0, 92.0, in_primary=False, share_basis="member",
                   counted_within="Phocaeicola vulgatus", judged_on_own_rank=False, secondary_percent=11.3,
                   gtdb="Phocaeicola SPECIV4_34405", unnamed=True, genus="Phocaeicola")
    distinct = _ranked("Blautia_A_luti", 0.0, 100.0, in_primary=False, share_basis="member",
                       counted_within="Blautia A wexlerae", judged_on_own_rank=True, secondary_percent=0.39,
                       gtdb="Blautia_A luti", genus="Blautia")
    assert not org.verdict(same).flagged
    assert org.verdict(distinct).flagged and org.verdict(distinct).is_issue
    # a record from before the field existed is treated as the same population: silent, never a duplicate
    legacy = _ranked("Blautia_A_luti", 0.0, 100.0, in_primary=False, share_basis="member",
                     counted_within="Blautia A wexlerae", secondary_percent=0.39, gtdb="Blautia_A luti", genus="Blautia")
    assert not org.verdict(legacy).flagged
