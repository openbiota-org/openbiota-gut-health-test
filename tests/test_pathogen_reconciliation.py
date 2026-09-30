"""The pathogen screen must not contradict the organism inventory.

The screen aligns to its own 485 reference genomes; a relative it carries no
reference for has nowhere to put its conserved reads, so they land on the
nearest target. One sample read "Mediterraneibacter gnavus - large signal,
0.34% of analysed DNA, 19,156 fragments" on the pathogen page while every
other section of the same report said the species was absent: the reads were
Dorea hominis's and Mediterraneibacter lactaris's.
"""

from __future__ import annotations

import pytest

from openbiota import inventory
from openbiota.pathogens import reconcile

pytest.importorskip("openbiota.pathogens.reconcile")


def _record(**over):
    rec = {
        "target_id": "bacteria.mediterraneibacter_gnavus",
        "display_name": "Mediterraneibacter gnavus",
        "group": "bacteria",
        "interpretation_class": "background_or_decoy",
        "sequence_status": "supported_sequence",
        "display_status": "supported_sequence",
        "species_resolution": "resolved",
        "unique_supporting_fragments": 19156,
        "informative_regions_supported": 106,
        "informative_bases_covered": 286566,
        "informative_region_breadth_fraction": 0.08303296,
        "reference_breadth_fraction": 0.08074122,
        "median_alignment_identity": 0.98,
        "normalized_fragments_per_million": 3384.37,
        "counts_as_pathogen": True,
        "plain_statement": "Dna sequences consistent with Mediterraneibacter gnavus detected.",
        "reason_codes": [],
    }
    rec.update(over)
    return rec


def _inventory(records):
    blob = {"primary_lane": "jan26", "unclassified_percent": 1.0, "organisms": records}
    inv = inventory.from_json(blob)
    assert inv is not None
    return inv


def _relatives_only():
    return _inventory([
        {"species": "Dorea_hominis", "gtdb": "Dorea_D hominis", "percent": 0.102, "in_primary": True,
         "genus": "Dorea", "status": "supported", "detected_by": ["jan26", "globdb"],
         "formerly_listed_as": ["Mediterraneibacter_gnavus"], "methods": ["marker", "genome_sketch"]},
        {"species": "Mediterraneibacter_lactaris", "gtdb": "Mediterraneibacter lactaris", "percent": 1.125,
         "in_primary": True, "genus": "Mediterraneibacter", "status": "supported",
         "detected_by": ["jan26", "globdb"], "methods": ["marker", "genome_sketch"]},
    ])


def test_shared_sequence_is_not_a_finding_and_names_the_relatives() -> None:
    pathogens = {"results": [_record()], "counts": {"pathogen_count": 1, "supported": 1}}
    summary = reconcile.reconcile(pathogens, _relatives_only(), read_length_bp=147.0, species_of={})
    rec = pathogens["results"][0]
    assert summary["n_shared_sequence"] == 1 and summary["n_found"] == 0
    assert rec["inventory_agreement"] == "shared_sequence_from_relatives"
    assert rec["sequence_status"] == "ambiguous_signal" and rec["counts_as_pathogen"] is False
    assert rec["report_tier"] == "shared_sequence"
    assert "shared_sequence_from_relatives" in rec["reason_codes"]
    # the relatives that carry the reads are named, the older catalogue's label explained
    assert "Mediterraneibacter lactaris (1.12%)" in rec["inventory_relatives"]
    assert "Dorea hominis" in rec["plain_statement"] and "did not find" in rec["plain_statement"]
    # the evidence survives: fragments, regions and breadth are untouched
    assert rec["unique_supporting_fragments"] == 19156 and rec["reference_breadth_fraction"] == 0.08074122
    assert rec["inventory_evenness"] is not None and rec["inventory_evenness"] < 0.25
    # and the headline count no longer includes it
    assert pathogens["counts"]["pathogen_count"] == 0 and pathogens["counts"]["supported"] == 0


def test_a_species_the_inventory_found_keeps_its_call_and_quotes_the_share() -> None:
    inv = _inventory([
        {"species": "Mediterraneibacter_gnavus", "gtdb": "Mediterraneibacter gnavus", "percent": 6.9336,
         "in_primary": True, "genus": "Mediterraneibacter", "status": "supported",
         "detected_by": ["jan26", "globdb", "motus"], "methods": ["marker", "genome_sketch", "universal_marker"]},
    ])
    rec = _record(reference_breadth_fraction=0.851, informative_region_breadth_fraction=0.875,
                  unique_supporting_fragments=326669, informative_bases_covered=2_900_000,
                  median_alignment_identity=1.0, normalized_fragments_per_million=40681.5)
    pathogens = {"results": [rec], "counts": {}}
    summary = reconcile.reconcile(pathogens, inv, read_length_bp=147.0, species_of={})
    assert summary["n_found"] == 1 and summary["n_shared_sequence"] == 0
    assert rec["inventory_agreement"] == "found" and rec["sequence_status"] == "supported_sequence"
    assert rec["inventory_share_percent"] == 6.9336
    assert "6.93%" in rec["plain_statement"]


def test_the_screens_count_is_explained_when_it_dwarfs_the_share() -> None:
    inv = _inventory([
        {"species": "Mediterraneibacter_gnavus", "gtdb": "Mediterraneibacter gnavus", "percent": 0.0389,
         "in_primary": True, "genus": "Mediterraneibacter", "status": "supported",
         "detected_by": ["jan26", "globdb"], "methods": ["marker", "genome_sketch"]},
    ])
    pathogens = {"results": [_record(normalized_fragments_per_million=3507.6, reference_breadth_fraction=0.19,
                                     informative_region_breadth_fraction=0.198,
                                     informative_bases_covered=700_000)], "counts": {}}
    reconcile.reconcile(pathogens, inv, read_length_bp=147.0, species_of={})
    rec = pathogens["results"][0]
    assert rec["inventory_agreement"] == "found"
    assert "count_includes_shared_sequence" in rec["reason_codes"]
    assert "not competitive" in rec["plain_statement"]


def test_a_genuine_call_the_inventory_cannot_see_is_kept() -> None:
    """The inventory's catalogues have a detection floor this screen does not:
    a pathogen at 0.001% with evenly spread coverage must still be reported."""
    rec = _record(target_id="bacteria.listeria_monocytogenes", display_name="Listeria monocytogenes",
                  interpretation_class="conditional_opportunist", unique_supporting_fragments=900,
                  informative_regions_supported=400, informative_bases_covered=120_000,
                  informative_region_breadth_fraction=0.04, reference_breadth_fraction=0.04,
                  median_alignment_identity=0.972, normalized_fragments_per_million=159.0)
    pathogens = {"results": [rec], "counts": {}}
    summary = reconcile.reconcile(pathogens, _relatives_only(), read_length_bp=147.0, species_of={})
    assert summary["n_not_in_inventory"] == 1 and summary["n_shared_sequence"] == 0
    assert rec["sequence_status"] == "supported_sequence" and rec["counts_as_pathogen"] is True
    assert "detection floor" in rec["plain_statement"]


def test_a_pathotype_target_is_checked_under_the_species_it_resolves_to() -> None:
    """The ETBF target's genome is Bacteroides fragilis; the inventory holds no
    organism called "Enterotoxigenic Bacteroides fragilis (ETBF)"."""
    inv = _inventory([
        {"species": "Bacteroides_fragilis", "gtdb": "Bacteroides fragilis", "percent": 0.9893,
         "in_primary": True, "genus": "Bacteroides", "status": "supported",
         "detected_by": ["jan26", "globdb"], "methods": ["marker", "genome_sketch"]},
    ])
    rec = _record(target_id="bacteria.bacteroides_fragilis_etbf",
                  display_name="Enterotoxigenic Bacteroides fragilis (ETBF)",
                  interpretation_class="toxin_or_pathotype_dependent",
                  reference_breadth_fraction=0.785, informative_region_breadth_fraction=0.79,
                  informative_bases_covered=4_000_000, median_alignment_identity=1.0)
    pathogens = {"results": [rec], "counts": {}}
    reconcile.reconcile(pathogens, inv, read_length_bp=147.0,
                        species_of={"bacteria.bacteroides_fragilis_etbf": "Bacteroides fragilis"})
    assert rec["inventory_species"] == "Bacteroides fragilis"
    assert rec["inventory_agreement"] == "found" and rec["inventory_share_percent"] == 0.9893


def test_a_call_the_competitive_confirmation_rejected_is_shared_sequence() -> None:
    inv = _relatives_only()
    rejected = [{"species": "Mediterraneibacter_gnavus", "gtdb": "Mediterraneibacter gnavus",
                 "reads_belong_to": "Dorea hominis", "percent": 0.0, "in_primary": False}]
    rec = _record(reference_breadth_fraction=0.9, informative_region_breadth_fraction=0.9,
                  informative_bases_covered=3_000_000)   # even coverage: only the verdict withdraws it
    pathogens = {"results": [rec], "counts": {}}
    reconcile.reconcile(pathogens, inv, read_length_bp=147.0, rejected=rejected, species_of={})
    assert rec["inventory_agreement"] == "shared_sequence_from_relatives"
    assert "Dorea hominis" in rec["inventory_relatives"]
    assert "competitive whole-genome confirmation rejected" in rec["plain_statement"]


def test_only_named_bacterial_calls_are_checked() -> None:
    fungal = _record(target_id="fungi.candida_albicans", display_name="Candida albicans", group="fungi")
    trace = _record(target_id="bacteria.listeria_monocytogenes", display_name="Listeria monocytogenes",
                    sequence_status="candidate_signal", display_status="candidate_signal")
    pathogens = {"results": [fungal, trace], "counts": {}}
    summary = reconcile.reconcile(pathogens, _relatives_only(), read_length_bp=147.0, species_of={})
    assert summary["n_checked"] == 0
    assert "inventory_agreement" not in fungal and "inventory_agreement" not in trace


def test_no_inventory_leaves_every_record_untouched() -> None:
    rec = _record()
    pathogens = {"results": [rec], "counts": {"pathogen_count": 1}}
    summary = reconcile.reconcile(pathogens, None)
    assert summary["status"] == "inventory unavailable"
    assert rec["sequence_status"] == "supported_sequence" and pathogens["counts"]["pathogen_count"] == 1
