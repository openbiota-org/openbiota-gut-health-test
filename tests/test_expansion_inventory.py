"""The merged inventory under the expanded lanes (spec 0.8.4 §6, §7).

Pure fixtures: no reference files, no sample. What is pinned:
  * one organism under two names across lanes renders once
  * abundances are never summed across lanes; the primary lane's value
    is the share and a detection-only lane never supplies a number
  * two independent methods make a call supported; one marginal method
    leaves it provisional; a split is ambiguous and a complex
  * count categories and the baseline/expansion split are right
  * a rejected confirmation moves the record out of the organism list and
    into the rejected list with its reason, and the counts follow
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from openbiota import inventory
from openbiota.expansion import confirm


def _lanes():
    scoring = inventory.Lane(name="scoring", species={"Bacteroides_vulgatus": 5.0, "Blautia_wexlerae": 7.0},
                             rankable=True, method="marker")
    jan26 = inventory.Lane(name="jan26", species={"Phocaeicola_vulgatus": 4.8, "Blautia_wexlerae": 7.2}, rankable=False,
                           primary=True, method="marker", rows=[
        {"species": "Phocaeicola_vulgatus", "genus": "Phocaeicola", "percent": 4.8, "sgb": "SGB1814",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB1814", "gtdb": "Phocaeicola vulgatus", "unnamed": False},
        {"species": "Blautia_wexlerae", "genus": "Blautia", "percent": 7.2, "sgb": "SGB4837",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB4837", "gtdb": "Blautia_A wexlerae", "unnamed": False},
        {"species": "GGB1_SGB999", "genus": "Ruminococcus", "percent": 0.02, "sgb": "SGB999",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB999", "gtdb": "", "unnamed": True},
    ])
    globdb = inventory.Lane(name="globdb", species={}, rankable=False, method="genome_sketch", rows=[
        {"species": "Phocaeicola_vulgatus", "gtdb": "Phocaeicola vulgatus", "genus": "Phocaeicola", "percent": 4.1,
         "unnamed": False, "native_id": "globdb:GCF_000012825", "status": "supported"},
        {"species": "Blautia_A_MGYG000001338", "gtdb": "Blautia_A MGYG000001338", "genus": "Blautia_A", "percent": 6.35,
         "unnamed": True, "native_id": "globdb:MGYG000001338", "status": "supported"},
        {"species": "Gemmiger_MOTU40_023295", "gtdb": "Gemmiger MOTU40_023295", "genus": "Gemmiger", "percent": 0.03,
         "unnamed": True, "native_id": "globdb:MOTU40_023295", "status": "provisional"},
    ])
    motus = inventory.Lane(name="motus", species={}, rankable=False, method="universal_marker", rows=[
        {"species": "Phocaeicola_vulgatus", "gtdb": "Phocaeicola vulgatus", "genus": "Phocaeicola", "percent": 3.9,
         "unnamed": False, "native_id": "mOTUv4.0_000123", "status": "supported"},
    ])
    singlem = inventory.Lane(name="singlem", species={}, rankable=False, method="marker_window", detection_only=True,
                             unit="coverage", rows=[
        {"species": "Blautia_A_MGYG000001338", "gtdb": "Blautia_A MGYG000001338", "genus": "Blautia_A", "percent": None,
         "unnamed": True, "native_id": "singlem:s__Blautia_A MGYG000001338", "status": "supported"},
    ])
    return [scoring, jan26, globdb, motus, singlem]


def test_one_organism_under_two_names_renders_once() -> None:
    inv = inventory.build(_lanes())
    vulgatus = [o for o in inv.organisms if "vulgatus" in o.species]
    assert len(vulgatus) == 1, [o.species for o in vulgatus]
    o = vulgatus[0]
    assert set(o.lanes) == {"scoring", "jan26", "globdb", "motus"}
    assert o.native_ids["jan26"] == "SGB1814" and o.native_ids["motus"] == "mOTUv4.0_000123"


def test_abundance_is_the_primary_lane_and_never_a_sum() -> None:
    inv = inventory.build(_lanes())
    o = next(o for o in inv.organisms if "vulgatus" in o.species)
    assert o.percent == 4.8, "the share is the primary lane's own value"
    assert o.secondary_percent == 5.0, "the secondary reading is the largest other percent-like value, not a sum"
    assert abs(o.percent + (o.secondary_percent or 0) - 9.8) < 1e-9  # nothing was added into `percent`


def test_a_detection_only_lane_never_supplies_a_number() -> None:
    inv = inventory.build(_lanes())
    cluster = next(o for o in inv.organisms if "MGYG000001338" in (o.gtdb or ""))
    assert "singlem" in cluster.lanes
    assert cluster.secondary_percent == 6.35, "coverage from SingleM must not become a percent"
    assert not cluster.in_primary and cluster.percent == 0.0


def test_statuses_follow_the_methods() -> None:
    inv = inventory.build(_lanes())
    by = {o.gtdb or o.species: o for o in inv.organisms}
    assert by["Phocaeicola vulgatus"].status == "supported"          # four lanes, three methods
    # Two methods, but both read GlobDB's genomes: a cluster that duplicates a
    # GTDB species would absorb its reads in both. Provisional until confirmed.
    cluster = by["Blautia_A MGYG000001338"]
    assert cluster.status == "provisional" and "GlobDB" in cluster.confidence_basis
    weak = by["Gemmiger MOTU40_023295"]
    assert weak.status == "provisional" and "confirmation" in weak.confidence_basis
    tiny = by["GGB1_SGB999"]
    assert tiny.status == "provisional", "0.02% on the primary lane alone is marginal"


def test_count_categories_and_incremental_gain() -> None:
    inv = inventory.build(_lanes())
    by = {o.gtdb or o.species: o for o in inv.organisms}
    assert by["Phocaeicola vulgatus"].count_category == "named_species"
    assert by["Blautia_A MGYG000001338"].count_category == "unnamed_species_cluster"
    assert by["Blautia_A MGYG000001338"].incremental_gain == "expansion"
    assert by["Phocaeicola vulgatus"].incremental_gain == "baseline"
    c = inv.counts()
    assert c["category_named_species"] + c["category_unnamed_species_cluster"] + c["category_unresolved_complex"] + c["category_higher_rank"] == c["organisms"]
    # Blautia_A MGYG..., Gemmiger MOTU40... and the Jan26-only SGB999: all three are
    # gains over the installed baseline (Jan26 is itself a new lane); one is supported.
    assert c["expansion"] == 3 and c["expansion_supported"] == 0, "nothing new is supported before confirmation"


def test_json_round_trip_keeps_the_new_fields() -> None:
    inv = inventory.build(_lanes())
    back = inventory.from_json(inv.to_json())
    assert back is not None
    a = {o.species: o for o in inv.organisms}
    b = {o.species: o for o in back.organisms}
    for k, o in a.items():
        assert b[k].status == o.status and b[k].count_category == o.count_category
        assert b[k].incremental_gain == o.incremental_gain and dict(b[k].native_ids) == dict(o.native_ids)
    blob = inv.to_json()
    assert "comparison" in blob and "newly_supported" in blob["comparison"]
    assert blob["comparison"]["newly_supported"] == [], "newly supported is earned by confirmation"
    confirmed = confirm.apply_verdicts(blob, {"verdicts": [
        {"organism": "Blautia_A MGYG000001338", "status": "supported", "reason": "4,000 fragments over 60% of the genome",
         "unique_fragments": 4000, "breadth": 0.6, "identity": 0.99, "evenness": 0.8, "competitors": {"Blautia wexlerae": 300}}]})
    rec = next(r for r in confirmed["organisms"] if r.get("gtdb") == "Blautia_A MGYG000001338")
    assert rec["status"] == "supported"
    assert confirmed["counts"]["supported"] == sum(1 for r in confirmed["organisms"] if r["status"] == "supported")


def test_a_rejected_confirmation_leaves_the_list_but_not_the_file() -> None:
    inv = inventory.build(_lanes())
    blob = inv.to_json()
    verdict = {"verdicts": [{"organism": "Gemmiger MOTU40_023295", "status": "not_detected",
                             "reason": "12 unique fragments; reads in this genus map to Gemmiger formicilis",
                             "unique_fragments": 12, "breadth": 0.001, "identity": 0.95, "evenness": 0.1, "competitors": {}},
                            {"organism": "Blautia_A MGYG000001338", "status": "supported", "reason": "ok",
                             "unique_fragments": 4000, "breadth": 0.6, "identity": 0.99, "evenness": 0.8, "competitors": {}}]}
    out = confirm.apply_verdicts(blob, verdict)
    names = [(r.get("gtdb") or r["species"]) for r in out["organisms"]]
    assert "Gemmiger MOTU40_023295" not in names
    assert out["rejected"][0]["status"] == "rejected" and "Gemmiger formicilis" in out["rejected"][0]["confidence_basis"]
    assert out["n_organisms"] == len(out["organisms"]) == len(blob["organisms"])
    assert out["counts"]["rejected"] == 1
    assert out["comparison"]["rejected"][0]["organism"] == "Gemmiger MOTU40_023295"
    blautia = next(r for r in out["organisms"] if r.get("gtdb") == "Blautia_A MGYG000001338")
    assert blautia["status"] == "supported" and blautia["confirmation"]["unique_fragments"] == 4000


def test_strain_placements_attach_to_their_organism_only() -> None:
    from openbiota.expansion import strainphlan

    blob = inventory.build(_lanes()).to_json()
    analysis = {"status": "completed", "clades": [
        {"sgb": "SGB1814", "status": "placed", "nearest_reference": "GCA_000012825.1.fna", "distance_to_nearest": 0.02,
         "n_markers": 180, "polymorphic_rate": 0.001, "references": ["GCA_000012825.1.fna"], "note": "a placement"},
        {"sgb": "SGB4837", "status": "unresolved", "reason": "only 2 reference genomes"},
        {"sgb": "SGB999999", "status": "placed", "nearest_reference": "x", "n_markers": 30},
    ]}
    n = strainphlan.attach_to_inventory(blob, analysis)
    assert n == 1
    recs = {r.get("sgb"): r for r in blob["organisms"]}
    assert recs["SGB1814"]["strain"]["nearest_reference"] == "GCA_000012825.1.fna"
    assert "strain" not in recs["SGB4837"] or not recs["SGB4837"].get("strain")
    assert blob["counts"]["strain_resolved"] == 1
    back = inventory.from_json(blob)
    assert back is not None and len(back.strain_resolved) == 1
    # a placement never changes how many organisms there are
    assert len(back.organisms) == len(inventory.build(_lanes()).organisms)


def test_placement_joins_an_existing_marker_typing_record() -> None:
    """An organism the marker-typing lane resolved keeps that record and gains the placement."""
    from openbiota.expansion import strainphlan

    blob = {"organisms": [{"sgb": "SGB1814", "species": "X_y",
                           "strain": {"sgb": "SGB1814", "markers_resolved": 120, "callable_bases": 90000}}],
            "counts": {"strain_resolved": 1}}
    analysis = {"status": "completed", "clades": [
        {"sgb": "SGB1814", "status": "placed", "nearest_reference": "GCA_000012825.1.fna",
         "distance_to_nearest": 0.01, "n_markers": 100, "polymorphic_rate": 0.002}]}
    assert strainphlan.attach_to_inventory(blob, analysis) == 1
    s = blob["organisms"][0]["strain"]
    assert s["markers_resolved"] == 120 and s["nearest_reference"] == "GCA_000012825.1.fna"
    assert s["placement_method"].startswith("StrainPhlAn 4")
    assert blob["counts"]["strain_resolved"] == 1
    # a second attach is idempotent
    assert strainphlan.attach_to_inventory(blob, analysis) == 0


def test_sibling_species_from_disjoint_lanes_fold_into_one_complex() -> None:
    """A marker bin labelled X_B and a genome sketch placing the reads in X are one organism."""
    jan26 = inventory.Lane(name="jan26", species={"Blautia_wexlerae": 7.0}, rankable=False, primary=True, method="marker", rows=[
        {"species": "Blautia_wexlerae", "genus": "Blautia", "percent": 7.0, "sgb": "SGB4837",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB4837", "gtdb": "Blautia_A wexlerae_B", "unnamed": False}])
    genome = inventory.Lane(name="genome", species={"Blautia_A_wexlerae": 7.2}, rankable=False, method="genome_sketch", rows=[
        {"species": "Blautia_A_wexlerae", "gtdb": "Blautia_A wexlerae", "genus": "Blautia_A", "percent": 7.2, "unnamed": False,
         "native_id": "GCF_000000000.1"}])
    inv = inventory.build([jan26, genome], lane_cohorts={})
    hits = [o for o in inv.organisms if "wexlerae" in (o.gtdb or "")]
    assert len(hits) == 1, [(o.species, o.gtdb) for o in hits]
    o = hits[0]
    assert o.status == "ambiguous" and o.count_category == "unresolved_complex"
    assert o.gtdb == "Blautia_A wexlerae / wexlerae_B"
    assert set(o.lanes) == {"jan26", "genome"}
    assert o.percent == 7.0  # the primary lane's share, never a sum
    # a lane that saw both siblings keeps them apart
    genome2 = inventory.Lane(name="genome", species={"Blautia_A_wexlerae": 5.0, "Blautia_A_wexlerae_B": 2.0}, rankable=False,
                             method="genome_sketch", rows=[
        {"species": "Blautia_A_wexlerae", "gtdb": "Blautia_A wexlerae", "genus": "Blautia_A", "percent": 5.0, "unnamed": False},
        {"species": "Blautia_A_wexlerae_B", "gtdb": "Blautia_A wexlerae_B", "genus": "Blautia_A", "percent": 2.0, "unnamed": False}])
    inv2 = inventory.build([jan26, genome2], lane_cohorts={})
    assert sum(1 for o in inv2.organisms if "wexlerae" in (o.gtdb or "")) == 2


# --- one organism under two labels -------------------------------------------


class _Row:
    """A scoring-cohort row, as community.SpeciesRow presents it."""

    def __init__(self, species: str, percentile: float) -> None:
        self.species, self.percentile, self.cohort_prevalence, self.groups = species, percentile, 0.5, ()
        self.trace = self.rare_in_cohort = False


def test_placeholder_recognises_catalogue_identifiers() -> None:
    ph = inventory._placeholder
    for label in ("Gemmiger sp937890665", "Blautia_A MGYG000001338", "Lachnospira MOTU40_010844",
                  "Phocaeicola SPECIV4_34405", "Allofournierella GCMETA_00520113", "Alistipes_A HRGMV2_0140",
                  "MOTU40_058831", "MGYG000001032", "Unknown Eggerthella mOTUv4.0_078089", "UBA1417 MGYG000004815",
                  "Bacteroides GWHFDIJ00000000", "CAG-217 MGYG000001029"):
        assert ph(label), label
    for label in ("Blautia_A wexlerae_B", "Phocaeicola vulgatus", "Faecalibacterium prausnitzii_C",
                  "Collinsella aerofaciens G", "Escherichia coli"):
        assert not ph(label), label


def test_scoring_row_keeps_its_percentile_when_the_record_is_renamed_by_a_merge() -> None:
    """MP3 says Blautia_wexlerae (ranked); the genome lane keys it Blautia_A_wexlerae. One record, ranked."""
    scoring = inventory.Lane(name="scoring", species={"Blautia_wexlerae": 0.9}, rankable=True, method="marker",
                             rows=[{"species": "Blautia_wexlerae", "percent": 0.9, "gtdb": "Blautia_A wexlerae",
                                    "genus": "Blautia", "unnamed": False}])
    genome = inventory.Lane(name="genome", species={"Blautia_A_wexlerae": 7.2}, rankable=False, method="genome_sketch",
                            primary=True, rows=[{"species": "Blautia_A_wexlerae", "gtdb": "Blautia_A wexlerae",
                                                 "genus": "Blautia_A", "percent": 7.2, "unnamed": False,
                                                 "accession": "GCF_000000001.1"}])
    inv = inventory.build([scoring, genome], ranked_rows=[_Row("Blautia_wexlerae", 30.6)], lane_cohorts={})
    assert len(inv.organisms) == 1
    o = inv.organisms[0]
    assert o.percentile == 30.6 and o.percentile_source == "scoring cohort"
    assert "Blautia_wexlerae" in o.aliases or o.species == "Blautia_wexlerae"
    assert o.scoring_percent == 0.9


def test_records_naming_the_same_genome_are_one_organism() -> None:
    """A mOTUs cluster and the GlobDB genome that represents it: two lanes, one record, supported."""
    globdb = inventory.Lane(name="globdb", species={}, rankable=False, method="genome_sketch", rows=[
        {"species": "Eggerthella_MOTU40_078089", "gtdb": "Eggerthella MOTU40_078089", "genus": "Eggerthella",
         "percent": 0.06, "unnamed": True, "native_id": "globdb:MOTU40_078089", "genome_id": "MOTU40_078089"}])
    motus = inventory.Lane(name="motus", species={}, rankable=False, method="universal_marker", rows=[
        {"species": "Eggerthella_MOTU40_078089", "gtdb": "Eggerthella MOTU40_078089", "genus": "Eggerthella",
         "percent": 0.05, "unnamed": True, "native_id": "mOTUv4.0_078089", "genome_id": "MOTU40_078089"}])
    # and the genus-less form GlobDB writes, matched to the mOTUs genus by the shared genome
    globdb2 = inventory.Lane(name="globdb", species={}, rankable=False, method="genome_sketch", rows=[
        {"species": "MOTU40_058831", "gtdb": "MOTU40_058831", "genus": "", "percent": 0.01, "unnamed": True,
         "native_id": "globdb:MOTU40_058831", "genome_id": "MOTU40_058831"}])
    motus2 = inventory.Lane(name="motus", species={}, rankable=False, method="universal_marker", rows=[
        {"species": "JAAWCD01_MOTU40_058831", "gtdb": "JAAWCD01 MOTU40_058831", "genus": "JAAWCD01", "percent": 0.03,
         "unnamed": True, "native_id": "mOTUv4.0_058831", "genome_id": "MOTU40_058831"}])
    inv = inventory.build([globdb, motus], lane_cohorts={})
    assert len(inv.organisms) == 1
    o = inv.organisms[0]
    assert set(o.lanes) == {"globdb", "motus"} and o.status == "supported"
    assert o.genome_ids == ("MOTU40_078089",) and o.count_category == "unnamed_species_cluster"
    inv2 = inventory.build([globdb2, motus2], lane_cohorts={})
    assert len(inv2.organisms) == 1
    o2 = inv2.organisms[0]
    assert o2.genus == "JAAWCD01" and o2.display in ("JAAWCD01 MOTU40_058831", "MOTU40_058831")
    assert o2.native_ids == {"globdb": "globdb:MOTU40_058831", "motus": "mOTUv4.0_058831"}


def test_confirmation_genome_comes_from_the_lanes_before_the_species_representative() -> None:
    o = inventory.Organism(species="Phocaeicola_SPECIV4_34405", percent=0.0, genus="Phocaeicola", unnamed=True,
                           gtdb="Phocaeicola SPECIV4_34405", native_ids={"globdb": "globdb:SPECIV4_34405"},
                           genome_ids=("SPECIV4_34405",), lanes=("globdb",))
    from unittest import mock

    with mock.patch("openbiota.expansion.genomes.in_globdb", return_value=True):
        assert confirm._genome_id_for(o) == "SPECIV4_34405"
    o2 = inventory.Organism(species="Escherichia_coli", percent=0.0, genus="Escherichia", gtdb="Escherichia coli",
                            native_ids={"rescue": "panel:BackhedF_2015_ERR525735_bin.17.fna.gz"},
                            genome_ids=("file:/nonexistent/genome.fna.gz",), lanes=("rescue",))
    with mock.patch("openbiota.expansion.confirm._r232_representative", return_value="GCF_000005845.2"), \
            mock.patch("openbiota.expansion.genomes.in_globdb", return_value=False):
        assert confirm._genome_id_for(o2) == "GCF_000005845.2"


def test_evenness_is_breadth_against_the_random_expectation() -> None:
    """46,000 fragments over 64% of a genome at 7x is even coverage of a divergent strain, not a hotspot."""
    assert confirm._evenness(0.636, 7.18) > 0.6
    # a shadow: 0.5% of the genome at a depth that would have covered a quarter of it
    assert confirm._evenness(0.005, 0.3) < 0.05
    assert confirm._evenness(0.0, 0.0) == 0.0
    assert confirm._evenness(1.0, 30.0) == 1.0 and confirm._evenness(0.99, 30.0) > 0.98


# --- one share of the whole ------------------------------------------------------


def test_unify_shares_splits_a_genus_total_by_unique_mapping_and_estimates_the_rest() -> None:
    blob = {"primary_lane": "jan26", "unclassified_percent": 10.0, "organisms": [
        {"species": "Lachnospira_eligens", "percent": 5.0, "in_primary": True, "genus": "Lachnospira",
         "gtdb": "Lachnospira eligens", "status": "supported", "detected_by": ["jan26", "globdb"], "secondary_percent": 5.2},
        {"species": "Lachnospira_eligens_A", "percent": 0.0, "in_primary": False, "genus": "Lachnospira",
         "gtdb": "Lachnospira eligens_A", "status": "supported", "detected_by": ["globdb", "motus"], "secondary_percent": 7.0},
        {"species": "Lachnospira_sp000436475", "percent": 0.0, "in_primary": False, "genus": "Lachnospira",
         "gtdb": "Lachnospira sp000436475", "unnamed": True, "status": "provisional", "detected_by": ["rescue"],
         "secondary_percent": 0.04},
        {"species": "Zag1_MGYG000001838", "percent": 0.0, "in_primary": False, "genus": "Zag1", "unnamed": True,
         "gtdb": "Zag1 MGYG000001838", "status": "supported", "detected_by": ["globdb", "kraken"], "secondary_percent": 0.5},
        {"species": "Phocaeicola_vulgatus", "percent": 20.0, "in_primary": True, "genus": "Phocaeicola",
         "gtdb": "Phocaeicola vulgatus", "status": "supported", "detected_by": ["jan26", "globdb"], "secondary_percent": 25.0},
    ]}
    confirmation = {"genomes": {
        "g1": {"organism": "Lachnospira eligens", "unique_fragments": 3000, "genome_length": 3_000_000},
        "g2": {"organism": "Lachnospira eligens_A", "unique_fragments": 6000, "genome_length": 3_000_000},
    }}
    inventory.unify_shares(blob, confirmation)
    by = {o["species"]: o for o in blob["organisms"]}
    # the genus total (5.0) is conserved and divided 1:2 by unique reads per megabase
    assert abs(by["Lachnospira_eligens"]["percent"] - 5.0 / 3) < 1e-6
    assert abs(by["Lachnospira_eligens_A"]["percent"] - 10.0 / 3) < 1e-6
    assert by["Lachnospira_eligens"]["share_basis"] == by["Lachnospira_eligens_A"]["share_basis"] == "split"
    assert by["Lachnospira_eligens"]["marker_percent"] == 5.0 and by["Lachnospira_eligens_A"]["in_primary"]
    # a single-method relative with no mapping is counted within the genus, without a share
    weak = by["Lachnospira_sp000436475"]
    assert weak["share_basis"] == "member" and weak["percent"] == 0.0 and not weak["in_primary"]
    assert weak["counted_within"] in ("Lachnospira eligens", "Lachnospira eligens_A")
    # a genus the marker lane never placed is estimated (too few pairs for a scale: factor 1)
    est = by["Zag1_MGYG000001838"]
    assert est["share_basis"] == "estimated" and abs(est["percent"] - 0.5) < 1e-6 and est["in_primary"]
    # untouched marker readings stay as they were
    assert by["Phocaeicola_vulgatus"]["share_basis"] == "marker" and by["Phocaeicola_vulgatus"]["percent"] == 20.0
    inv = inventory.from_json(blob)
    assert inv is not None and abs(inv.unplaced_percent - 9.5) < 1e-6
    assert blob["share_unification"]["n_split"] == 2 and blob["share_unification"]["n_estimated"] == 1
    # running it again starts from the marker readings, not from the split ones
    inventory.unify_shares(blob, confirmation)
    by = {o["species"]: o for o in blob["organisms"]}
    assert abs(by["Lachnospira_eligens"]["percent"] - 5.0 / 3) < 1e-6


def test_a_rejected_calls_marker_share_goes_to_the_relative_the_reads_belong_to() -> None:
    """Blautia massiliensis at 5.5% was rejected because its reads mapped to
    B. caecimuris; the whole must still add up, so the share goes with the reads."""
    blob = {"primary_lane": "jan26", "unclassified_percent": 10.0, "organisms": [
        {"species": "Blautia_caecimuris", "percent": 0.9, "in_primary": True, "genus": "Blautia",
         "gtdb": "Blautia_A caecimuris", "status": "supported", "detected_by": ["jan26", "globdb"],
         "reference_reading": 0.9, "scoring_percent": 0.9, "percentile": 60.0, "secondary_percent": 1.0},
        {"species": "Hominimerdicola_sp900066445", "percent": 0.0, "in_primary": False, "genus": "Hominimerdicola",
         "gtdb": "Hominimerdicola sp900066445", "unnamed": True, "status": "supported", "detected_by": ["globdb", "motus"],
         "secondary_percent": 0.7},
    ], "rejected": [
        {"species": "Blautia_massiliensis", "percent": 5.5, "in_primary": True, "genus": "Blautia", "status": "rejected",
         "detected_by": ["jan26"], "reads_belong_to": "Blautia_A caecimuris"},
        {"species": "Ruminococcus_bicirculans", "percent": 0.8, "in_primary": True, "genus": "Ruminococcus",
         "status": "rejected", "detected_by": ["jan26"], "reads_belong_to": "Hominimerdicola sp900066445"},
        {"species": "Streptococcus_anginosus", "percent": 0.004, "in_primary": True, "status": "rejected",
         "detected_by": ["jan26"], "reads_belong_to": "Streptococcus thermophilus"},   # no such record: nothing to do
    ]}
    inventory.unify_shares(blob, {"genomes": {}})
    by = {o["species"]: o for o in blob["organisms"]}
    target = by["Blautia_caecimuris"]
    assert abs(target["percent"] - 6.4) < 1e-9 and target["marker_percent"] == 0.9 and target["share_basis"] == "marker"
    assert target["absorbed_percent"] == 5.5 and target["absorbed_from"] == ["Blautia massiliensis"]
    # the share now differs sevenfold from what the reference catalogue read under this name: not comparable
    assert target["reference_conflict"] is True
    # a relative the marker catalogue had no entry for was estimated from its whole-genome reading; the
    # marker lane's reading of the same population (under the rejected name) replaces that estimate,
    # which is no longer drawn from the unclassified band
    newly = by["Hominimerdicola_sp900066445"]
    assert newly["in_primary"] and newly["share_basis"] == "absorbed" and abs(newly["percent"] - 0.8) < 1e-9
    assert newly["absorbed_from"] == ["Ruminococcus bicirculans"]
    su = blob["share_unification"]
    assert su["n_absorbed"] == 2 and abs(su["absorbed_percent_total"] - 6.3) < 1e-9
    assert su["n_estimated"] == 0 and su["estimated_percent_total"] == 0
    inv = inventory.from_json(blob)
    assert inv is not None and abs(inv.composition_total - 7.2) < 1e-9 and abs(inv.unplaced_percent - 10.0) < 1e-9
    # the fields round-trip, and the estimated band is never charged for an absorbed share
    again = inventory.from_json(inv.to_json())
    assert again is not None and again.get("Blautia_caecimuris").absorbed_percent == 5.5
    assert again.get("Blautia_caecimuris").absorbed_from == ("Blautia massiliensis",)
    assert again.get("Blautia_caecimuris").level_percentile is None
    # running it again starts over: the absorbed share is not added twice
    inventory.unify_shares(blob, {"genomes": {}})
    by = {o["species"]: o for o in blob["organisms"]}
    assert abs(by["Blautia_caecimuris"]["percent"] - 6.4) < 1e-9 and by["Blautia_caecimuris"]["absorbed_from"] == ["Blautia massiliensis"]
    assert abs(by["Hominimerdicola_sp900066445"]["percent"] - 0.8) < 1e-9


def test_status_does_not_depend_on_which_lane_a_set_yields_first() -> None:
    """Three lanes of one method with the primary at 1.5%: supported, whatever
    order the set of lane names comes out in (it varied by hash seed, and with
    it the confirmation cache key)."""
    lanes = {n: inventory.Lane(name=n, species={}, rankable=False, method="marker", primary=(n == "jan26"))
             for n in ("scoring", "extended", "jan26")}
    for order in (("scoring", "extended", "jan26"), ("extended", "jan26", "scoring"), ("jan26", "scoring", "extended")):
        rec = {"lanes": set(order), "by_lane": {"scoring": 1.1, "extended": 1.4, "jan26": 1.5}, "native_ids": {}}
        status, basis, methods = inventory._status_of(rec, lanes, "jan26", 1.5)
        assert status == "supported" and basis.startswith("primary marker lane at 1.50%") and methods == ("marker",)
    # without the primary lane, the provisional note names the lane with the largest reading
    rec = {"lanes": {"scoring", "extended"}, "by_lane": {"scoring": 0.2, "extended": 0.9}, "native_ids": {}}
    status, basis, _ = inventory._status_of(rec, lanes, "jan26", 0.0)
    assert status == "provisional" and "one method (extended) at 0.90" in basis


def test_the_relative_is_read_back_out_of_a_cached_verdict() -> None:
    from openbiota.expansion.confirm import relative_of
    assert relative_of({"status": "not_detected", "relative": "Blautia caecimuris", "reason": "x"}) == "Blautia caecimuris"
    assert relative_of({"status": "not_detected", "reason": (
        "11,674 unique fragments at 98.6% identity, while Blautia caecimuris took 107,278 at 99.0% in the same "
        "competition; the reads belong to the relative")}) == "Blautia caecimuris"
    assert relative_of({"status": "not_detected", "reason": "7 unique fragments; reads in this genus map to "
                        "CAG-433 sp000433675 (1,234 fragments)"}) == "CAG-433 sp000433675"
    assert relative_of({"status": "not_detected", "reason": "3 unique fragments after competition"}) == ""
    assert relative_of({"status": "supported", "reason": "while X took 1 at"}) == ""


def test_deviation_is_against_the_typical_carrier_on_the_reference_lane() -> None:
    o = inventory.Organism(species="Blautia_wexlerae", percent=7.0, percentile=97.0, scoring_percent=6.0,
                           reference_percent=1.2, reference_reading=6.0)
    assert abs(o.deviation_percent - 400.0) < 1e-9
    low = inventory.Organism(species="X_y", percent=0.2, percentile=3.0, scoring_percent=0.2,
                             reference_percent=2.0, reference_reading=0.2)
    assert abs(low.deviation_percent - (-90.0)) < 1e-9
    assert inventory.Organism(species="Z_z", percent=1.0).deviation_percent is None
    from openbiota.pdforganisms import _deviation_text
    assert _deviation_text(400.0) == "+400%" and _deviation_text(-90.0) == "\u221290%" and _deviation_text(2.0) == "\u2248 typical"
    assert _deviation_text(10883.0) == "\u00d7110"


def test_a_rank_among_a_handful_of_carriers_is_not_stated() -> None:
    """Four reference carriers cannot place a level; the row says so and carries no flag."""
    o = inventory.Organism(species="Mediterraneibacter_gnavus", percent=0.45, percentile=97.0, carrier_percentile=100.0,
                           reference_percent=0.108, reference_reading=0.1135, reference_carriers=4, prevalence=0.04,
                           percentile_source="globdb cohort, n=100")
    assert o.few_reference_carriers and o.level_percentile is None and o.deviation_percent is None
    from openbiota import organisms as org
    assert not org.verdict(o).flagged
    ok = inventory.replace(o, reference_carriers=40)
    assert ok.level_percentile == 100.0 and ok.deviation_percent is not None and ok.deviation_percent > 0


def test_a_label_naming_another_species_is_not_a_lookup_alias() -> None:
    """Jun23 called SGB4571 Ruminococcus gnavus; Jan26 and GTDB place it in Dorea hominis. The old label is
    kept as a record, and a lookup for R. gnavus does not return D. hominis."""
    jan26 = inventory.Lane(name="jan26", species={"Dorea_hominis": 0.1}, rankable=False, primary=True, method="marker", rows=[
        {"species": "Dorea_hominis", "genus": "Dorea", "percent": 0.1, "sgb": "SGB4571",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB4571", "gtdb": "Dorea_D hominis", "unnamed": False}])
    extended = inventory.Lane(name="extended", species={"Ruminococcus_gnavus": 0.1}, rankable=False, method="marker", rows=[
        {"species": "Ruminococcus_gnavus", "genus": "Ruminococcus", "percent": 0.1, "sgb": "SGB4571",
         "sgb_db": "mpa_vJun23_CHOCOPhlAnSGB_202403", "native_id": "SGB4571", "gtdb": "Dorea_D hominis", "unnamed": False}])
    from openbiota.expansion import names

    if not names.ncbi_species_to_r232():
        pytest.skip("the NCBI-name bridge needs the GTDB R232 metadata on disk (make refs-expanded)")
    inv = inventory.build([jan26, extended], lane_cohorts={})
    assert len(inv.organisms) == 1
    o = inv.organisms[0]
    assert o.species == "Dorea_hominis", "the newer catalogue names the bin"
    assert "Mediterraneibacter_gnavus" in o.formerly_listed_as
    assert "Mediterraneibacter_gnavus" not in o.aliases
    assert inv.get("Ruminococcus_gnavus") is None and inv.get("Mediterraneibacter_gnavus") is None
    assert inv.get("Dorea_hominis") is o


def test_the_reference_catalogues_own_name_answers_a_lookup() -> None:
    """MetaPhlAn 3 called a bin Blautia producta; the merged record is keyed by the newer name, and a lookup
    by the MetaPhlAn 3 name finds it because that catalogue itself reported it under that name."""
    scoring = inventory.Lane(name="scoring", species={"Blautia_producta": 0.3}, rankable=True, method="marker",
                             rows=[{"species": "Blautia_producta", "percent": 0.3, "gtdb": "Blautia producta",
                                    "genus": "Blautia", "unnamed": False, "native_id": "Blautia_producta"}])
    jan26 = inventory.Lane(name="jan26", species={"Blautia_celeris": 0.3}, rankable=False, primary=True, method="marker", rows=[
        {"species": "Blautia_celeris", "genus": "Blautia", "percent": 0.3, "sgb": "SGB4794",
         "sgb_db": "mpa_vJan26_CHOCOPhlAnSGB_202605", "native_id": "SGB4794", "gtdb": "Blautia celeris", "unnamed": False}])
    inv = inventory.build([scoring, jan26], lane_cohorts={})
    # the two records may or may not fold depending on the crosswalk on disk; either way the MP3 name answers
    hit = inv.get("Blautia_producta")
    assert hit is not None and (hit.native_ids or {}).get("scoring") == "Blautia_producta"


def test_an_uncommon_opportunist_far_above_its_carriers_is_overgrown() -> None:
    """Clostridium innocuum: 0.041%, +372% of the typical carrier, above 95% of all
    reference adults, 81st among the 29% who carry it. Ranked among carriers alone it
    escaped the attention page entirely; it is a vancomycin-resistant opportunist."""
    from openbiota import organisms as org

    innocuum = inventory.Organism(
        species="Clostridium_innocuum", percent=0.0405, in_primary=True, genus="Clostridium",
        percentile=95.5, carrier_percentile=80.7, prevalence=0.2904, reference_percent=0.0183,
        reference_reading=0.0863, scoring_percent=0.0433, percentile_source="scoring cohort",
        reference_carriers=879, methods=("genome_sketch", "marker", "read_classification", "universal_marker"),
        status="supported")
    v = org.verdict(innocuum)
    assert v.cls == org.OPPORTUNIST
    assert v.is_issue and v.flag == "high"
    assert "expanded and uncommon" in v.flag_reason and "372%" in v.flag_reason and "29% carry it" in v.flag_reason
    # one method is not enough, a trace is not enough, and a common organism is judged among its carriers
    single = replace(innocuum, methods=("marker",))
    assert not org.verdict(single).flagged
    absent = replace(innocuum, in_primary=False, percent=0.0)
    assert not org.verdict(absent).flagged
    mild = replace(innocuum, reference_reading=0.02)        # +9% of the typical carrier
    assert not org.verdict(mild).flagged
    low_in_carriers = replace(innocuum, carrier_percentile=40.0)
    assert not org.verdict(low_in_carriers).flagged
    # a conditional resident that is simply present in a population that often lacks it stays off the page
    plebeius = replace(innocuum, species="Phocaeicola_plebeius", genus="Phocaeicola", percent=3.747)
    assert org.verdict(plebeius).cls != org.OPPORTUNIST and not org.verdict(plebeius).is_issue


def test_a_member_is_judged_on_its_own_rank_only_when_it_is_a_distinct_population() -> None:
    """GlobDB's "Phocaeicola SPECIV4_34405" at 11.3% inside Phocaeicola vulgatus's
    11.7% is the same population under another catalogue's name and is judged
    once, as vulgatus. Blautia luti at 0.39% inside Blautia wexlerae's 7.2% is a
    distinct population the competition told apart, with a level of its own.
    And where the relative has no rank at all, the member's lane rank is the
    only reading of that population the report has."""
    from openbiota import organisms as org

    blob = {"primary_lane": "jan26", "unclassified_percent": 5.0, "organisms": [
        {"species": "Phocaeicola_vulgatus", "gtdb": "Phocaeicola vulgatus", "percent": 11.75, "in_primary": True,
         "genus": "Phocaeicola", "status": "supported", "detected_by": ["jan26", "globdb"], "percentile": 93.4,
         "carrier_percentile": 92.9, "reference_carriers": 2687, "reference_percent": 4.29, "reference_reading": 22.3},
        {"species": "Phocaeicola_SPECIV4_34405", "gtdb": "Phocaeicola SPECIV4_34405", "percent": 0.0, "in_primary": False,
         "genus": "Phocaeicola", "unnamed": True, "status": "supported", "detected_by": ["globdb", "kraken"],
         "secondary_percent": 11.27, "percentile": 92.0, "carrier_percentile": 91.9, "reference_carriers": 60,
         "percentile_source": "globdb cohort, n=100", "reference_percent": 3.0, "reference_reading": 11.27},
        {"species": "Blautia_A_wexlerae", "gtdb": "Blautia_A wexlerae", "percent": 7.21, "in_primary": True,
         "genus": "Blautia", "status": "supported", "detected_by": ["jan26", "globdb"], "reference_conflict": True,
         "percentile": 13.9, "carrier_percentile": 13.9, "reference_carriers": 2644},
        {"species": "Blautia_A_luti", "gtdb": "Blautia_A luti", "percent": 0.0, "in_primary": False, "genus": "Blautia",
         "status": "supported", "detected_by": ["globdb", "motus"], "secondary_percent": 0.391,
         "percentile": 100.0, "carrier_percentile": 100.0, "reference_carriers": 40,
         "percentile_source": "globdb cohort, n=100", "reference_percent": 0.05, "reference_reading": 0.391},
        {"species": "Blautia_A_MGYG000001338", "gtdb": "Blautia_A MGYG000001338", "percent": 0.0, "in_primary": False,
         "genus": "Blautia", "unnamed": True, "status": "supported", "detected_by": ["globdb", "kraken"],
         "secondary_percent": 6.35, "percentile": 90.3, "carrier_percentile": 90.3, "reference_carriers": 70,
         "percentile_source": "globdb cohort, n=100", "reference_percent": 2.0, "reference_reading": 6.35},
    ]}
    inventory.unify_shares(blob, {"genomes": {}})
    by = {o["species"]: o for o in blob["organisms"]}
    assert by["Phocaeicola_SPECIV4_34405"]["share_basis"] == "member"
    assert by["Phocaeicola_SPECIV4_34405"]["judged_on_own_rank"] is False      # same population, vulgatus has a rank
    assert by["Blautia_A_luti"]["judged_on_own_rank"] is True                  # distinct population
    assert by["Blautia_A_MGYG000001338"]["judged_on_own_rank"] is True         # same population, but wexlerae has no rank
    inv = inventory.from_json(blob)
    assert inv is not None
    flags = {v.organism.species: v for v in org.verdicts(list(inv.organisms))}
    assert flags["Phocaeicola_vulgatus"].flagged and not flags["Phocaeicola_SPECIV4_34405"].flagged
    assert flags["Blautia_A_luti"].flagged and flags["Blautia_A_luti"].is_issue
    assert flags["Blautia_A_MGYG000001338"].flagged
    # the flag survives a round trip through results.json
    again = inventory.from_json(inv.to_json())
    assert again is not None and again.get("Blautia_A_luti").judged_on_own_rank is True
    assert inventory._distinct_from(0.39, 7.2) and inventory._distinct_from(0.09, 0.008)
    assert not inventory._distinct_from(11.27, 11.75) and not inventory._distinct_from(10.2, 6.4)
