"""The whole-genome lane and how it joins the inventory.

Three things must hold. A genome-lane hit for an organism a marker lane
already named must fold into that record, under either the current name or
the GTDB alias, never a second row. A genome-lane hit for an organism no
marker lane saw is a new organism with no share of the composition. And
the lane's own quality flag is respected: a hit whose coverage model did
not fit is kept but not confident.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from openbiota import inventory
from openbiota.engines import sylph

PROFILE = sorted(Path("results/_sylph").glob("*.sylph.tsv")) if Path("results/_sylph").is_dir() else []


def _lane(name: str, species: dict[str, float], *, rankable: bool = False, primary: bool = False,
          rows=()) -> inventory.Lane:
    return inventory.Lane(name=name, species=species, rankable=rankable, primary=primary, rows=rows)


def test_accession_is_read_from_the_genome_path() -> None:
    assert sylph._accession(
        "gtdb_genomes_reps_r232/database/GCF/964/248/265/GCF_964248265.1_genomic.fna.gz"
    ) == "GCF_964248265.1"
    assert sylph._accession("GCA_000432195.1_genomic.fna") == "GCA_000432195.1"


def test_a_low_correction_hit_is_kept_but_not_confident() -> None:
    hit = sylph.GenomeHit(
        lineage="d__Bacteria;s__X y", species="X y", genus="X", family="", domain="Bacteria",
        taxonomic_abundance=0.01, sequence_abundance=0.01, ani=96.0, coverage=0.1,
        correction="LOW", accession="GCA_1.1",
    )
    assert not hit.confident
    from dataclasses import replace

    assert replace(hit, correction="HIGH").confident
    assert replace(hit, correction="1.42").confident


def test_placeholder_names_are_recognised() -> None:
    mk = lambda sp: sylph.GenomeHit(  # noqa: E731
        lineage="", species=sp, genus=sp.split()[0], family="", domain="Bacteria",
        taxonomic_abundance=1, sequence_abundance=1, ani=99, coverage=1,
        correction="1.0", accession="")
    assert mk("Gemmiger sp963557315").placeholder
    assert not mk("Phocaeicola vulgatus").placeholder
    assert not mk("Blautia_A fusiformis").placeholder


def test_a_genome_hit_folds_into_the_marker_record_under_the_current_name() -> None:
    """Marker catalogue says Mediterraneibacter_gnavus; GTDB r232 agrees."""
    inv = inventory.build([
        _lane("extended", {"Mediterraneibacter_gnavus": 13.0}, primary=True,
              rows=[{"species": "Mediterraneibacter_gnavus", "percent": 13.0, "genus": "Mediterraneibacter",
                     "sgb": "SGB4584"}]),
        _lane("genome", {"Mediterraneibacter_gnavus": 4.7},
              rows=[{"species": "Mediterraneibacter_gnavus", "gtdb": "Mediterraneibacter gnavus",
                     "genus": "Mediterraneibacter", "percent": 4.7}]),
    ])
    assert len(inv) == 1
    o = inv.organisms[0]
    assert set(o.lanes) == {"extended", "genome"}
    assert o.in_primary and o.percent == pytest.approx(13.0)
    assert o.secondary_percent == pytest.approx(4.7)


def test_a_genome_hit_folds_into_an_unnamed_bin_through_the_gtdb_alias() -> None:
    """The marker lane's bare bin and the genome lane's GTDB name are one organism."""
    from openbiota import gtdb

    if not gtdb.available() or gtdb.lookup("SGB5809") is None:
        pytest.skip("no GTDB mapping installed")
    alias = gtdb.lookup("SGB5809").display  # Dialister sp000434475 in r220
    inv = inventory.build([
        _lane("extended", {"GGB4266_SGB5809": 0.5}, primary=True,
              rows=[{"species": "GGB4266_SGB5809", "percent": 0.5, "genus": "GGB4266",
                     "sgb": "SGB5809", "unnamed": True}]),
        _lane("genome", {alias.replace(" ", "_"): 0.4},
              rows=[{"species": alias.replace(" ", "_"), "gtdb": alias,
                     "genus": alias.split()[0], "percent": 0.4, "unnamed": True}]),
    ])
    assert len(inv) == 1, [o.species for o in inv]
    assert set(inv.organisms[0].lanes) == {"extended", "genome"}


def test_one_placeholder_each_side_in_a_genus_is_one_organism() -> None:
    """GTDB's sp codes move with the representative between releases."""
    inv = inventory.build([
        _lane("extended", {"GGB1_SGB999999": 8.0}, primary=True,
              rows=[{"species": "GGB1_SGB999999", "percent": 8.0, "genus": "Gemmiger",
                     "sgb": "SGB999999", "unnamed": True}]),
        _lane("genome", {"Gemmiger_sp963557315": 8.7},
              rows=[{"species": "Gemmiger_sp963557315", "gtdb": "Gemmiger sp963557315",
                     "genus": "Gemmiger", "percent": 8.7, "unnamed": True}]),
    ])
    assert len(inv) == 1


def test_two_placeholders_one_side_stay_apart() -> None:
    """Ambiguity is preserved, not resolved by guessing."""
    inv = inventory.build([
        _lane("extended", {"GGB1_SGB999999": 1.0}, primary=True,
              rows=[{"species": "GGB1_SGB999999", "percent": 1.0, "genus": "Blautia",
                     "sgb": "SGB999999", "unnamed": True}]),
        _lane("genome", {"Blautia_sp900555025": 0.5, "Blautia_sp958352855": 0.3},
              rows=[{"species": "Blautia_sp900555025", "gtdb": "Blautia sp900555025",
                     "genus": "Blautia", "percent": 0.5, "unnamed": True},
                    {"species": "Blautia_sp958352855", "gtdb": "Blautia sp958352855",
                     "genus": "Blautia", "percent": 0.3, "unnamed": True}]),
    ])
    assert len(inv) == 3


def test_a_genome_only_organism_has_no_share_of_the_composition() -> None:
    inv = inventory.build([
        _lane("extended", {"A_b": 50.0}, primary=True,
              rows=[{"species": "A_b", "percent": 50.0, "genus": "A"}]),
        _lane("genome", {"C_d": 2.0}, rows=[{"species": "C_d", "gtdb": "C d", "genus": "C", "percent": 2.0}]),
    ])
    c = inv.get("C_d")
    assert c is not None and not c.in_primary
    assert c.percent == 0.0 and c.secondary_percent == pytest.approx(2.0)
    assert inv.composition_total == pytest.approx(50.0)


@pytest.mark.skipif(not PROFILE, reason="no Sylph profile on disk")
def test_real_profile_parses_with_a_lineage_for_every_hit() -> None:
    hits = sylph._parse(PROFILE[0])
    assert hits, PROFILE[0]
    assert all(h.lineage and h.species for h in hits)
    # Every hit clears GTDB's species boundary; that is the profiling floor.
    assert all(h.ani >= sylph.MINIMUM_ANI for h in hits)


def test_parser_reads_the_documented_columns() -> None:
    header = ("Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t"
              "True_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\t"
              "Mean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name")
    row = ("S\tx/GCF/964/248/265/GCF_964248265.1_genomic.fna.gz\t11.70\t17.27\t99.21\t67.07\tNA-NA\tHIGH\t"
           "NA-NA\t47\t45.97\t17830/22779\t99.21\t136\tNZ_X Phocaeicola vulgatus")
    text = header + "\n" + row + "\n"
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    assert sylph._accession(rows[0]["Genome_file"]) == "GCF_964248265.1"
    assert rows[0]["Eff_lambda"] == "HIGH"
