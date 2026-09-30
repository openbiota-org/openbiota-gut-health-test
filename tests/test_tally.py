"""Tally, mate collapse, decoy rejection and normalisation.

Covers spec test cases 1-5 and 7.
"""

from __future__ import annotations

import pytest

from openbiota.panels import PanelSet
from openbiota.references import ReferenceDatabase
from openbiota.tally import HitParser, normalise_query_id, tally

from .conftest import FIELDS, hit_line


def _run(
    database: ReferenceDatabase,
    panel_set: PanelSet,
    r1: list[str],
    r2: list[str] | None = None,
):
    streams = [("R1", iter(r1))]
    if r2 is not None:
        streams.append(("R2", iter(r2)))
    return tally(streams=streams, database=database, panel_set=panel_set, fields=FIELDS)


# --------------------------------------------------------------------------- #
# 1. mate collapse
# --------------------------------------------------------------------------- #


def test_mate_collapse_counts_fragments_not_reads(database, panel_set):
    """40 fragments appearing as 80 read hits must count as 40."""
    r1 = [hit_line(f"READ{i:04d}", 0, "URDA", "P00001") for i in range(40)]
    r2 = [hit_line(f"READ{i:04d}", 0, "URDA", "P00001") for i in range(40)]

    result = _run(database, panel_set, r1, r2)

    assert result.total_hit_lines == 80
    assert result.fragments_with_hit == 40
    assert result.panel("urda").accepted_fragments == 40


def test_mate_collapse_strips_slash_suffixes(database, panel_set):
    """Ids that do carry /1 and /2 must still collapse to one fragment."""
    r1 = [hit_line(f"READ{i:04d}/1", 0, "URDA", "P00001") for i in range(10)]
    r2 = [hit_line(f"READ{i:04d}/2", 0, "URDA", "P00001") for i in range(10)]

    result = _run(database, panel_set, r1, r2)

    assert result.fragments_with_hit == 10
    assert result.panel("urda").accepted_fragments == 10


def test_normalise_query_id():
    assert normalise_query_id("LH00469:636:23VMTGLT4:8:1101:1003:20056") == (
        "LH00469:636:23VMTGLT4:8:1101:1003:20056"
    )
    assert normalise_query_id("READ1/1") == "READ1"
    assert normalise_query_id("READ1/2") == "READ1"
    # A trailing digit that is not a mate marker must survive untouched.
    assert normalise_query_id("SRR12345.1") == "SRR12345.1"


def test_higher_scoring_mate_wins(database, panel_set):
    """When the two mates disagree, the higher bit score decides the fragment."""
    r1 = [hit_line("READ1", 0, "URDA", "P00001", bitscore=50.0)]
    r2 = [hit_line("READ1", 2, "FRDA", "P00003", bitscore=120.0)]

    result = _run(database, panel_set, r1, r2)
    urda = result.panel("urda")

    assert result.fragments_with_hit == 1
    assert urda.accepted_fragments == 0
    assert urda.decoy_fragments == 1


# --------------------------------------------------------------------------- #
# 2. decoy rejection
# --------------------------------------------------------------------------- #


def test_decoy_best_hit_counts_as_rejected_never_as_target(database, panel_set):
    r1 = [hit_line(f"T{i}", 0, "URDA", "P00001") for i in range(12)]
    r1 += [hit_line(f"D{i}", 2, "FRDA", "P00003") for i in range(30)]

    result = _run(database, panel_set, r1)
    urda = result.panel("urda")

    assert urda.accepted_fragments == 12
    assert urda.decoy_fragments == 30
    assert sum(t.fragments for t in urda.targets) == 12
    assert sum(d.fragments for d in urda.decoys) == 30


def test_decoy_below_its_own_identity_floor_is_not_counted_as_decoy_match(database, panel_set):
    """A weak decoy hit is noise, not evidence of the decoy family."""
    r1 = [hit_line("READ1", 2, "FRDA", "P00003", pident=30.0)]

    result = _run(database, panel_set, r1)
    urda = result.panel("urda")

    assert urda.decoy_fragments == 0
    assert urda.decoys[0].rejected_low_identity == 1


# --------------------------------------------------------------------------- #
# 3. identity threshold
# --------------------------------------------------------------------------- #


def test_target_hit_below_min_identity_is_dropped(database, panel_set):
    """urda:URDA has min_identity 60."""
    r1 = [
        hit_line("KEEP1", 0, "URDA", "P00001", pident=60.0),
        hit_line("KEEP2", 0, "URDA", "P00001", pident=99.9),
        hit_line("DROP1", 0, "URDA", "P00001", pident=59.9),
        hit_line("DROP2", 0, "URDA", "P00001", pident=12.0),
    ]

    result = _run(database, panel_set, r1)
    urda = result.panel("urda")

    assert urda.accepted_fragments == 2
    assert urda.targets[0].rejected_low_identity == 2


def test_short_alignment_is_dropped(database, panel_set):
    """min_alignment_aa is 25 for this panel."""
    r1 = [
        hit_line("KEEP", 0, "URDA", "P00001", aln_len=25),
        hit_line("DROP", 0, "URDA", "P00001", aln_len=24),
    ]

    result = _run(database, panel_set, r1)
    urda = result.panel("urda")

    assert urda.accepted_fragments == 1
    assert urda.targets[0].rejected_short_alignment == 1


def test_per_panel_identity_floors_are_independent(database, panel_set):
    """butyrate:BUT has min_identity 55, urda:URDA has 60."""
    r1 = [
        hit_line("A", 3, "BUT", "P00004", pident=57.0),
        hit_line("B", 0, "URDA", "P00001", pident=57.0),
    ]

    result = _run(database, panel_set, r1)

    assert result.panel("butyrate").accepted_fragments == 1
    assert result.panel("urda").accepted_fragments == 0


# --------------------------------------------------------------------------- #
# 4. normalisation arithmetic
# --------------------------------------------------------------------------- #


def test_normalisation_matches_hand_computed_value(database, panel_set):
    """Hand computation.

        urdA:  120 fragments / 600 aa mean reference length = 0.2
        rpoB:  600 fragments / 1200 aa                      = 0.5
        copies per genome      = 0.2 / 0.5 = 0.4
        copies per 100 genomes = 40.0
    """
    r1 = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(120)]
    r1 += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(600)]

    result = _run(database, panel_set, r1)

    assert result.normalizer.fragments == 600
    assert result.normalizer.reference_mean_length_aa == pytest.approx(1200.0)
    assert result.normalizer.rate == pytest.approx(0.5)

    urda = result.panel("urda")
    assert urda.accepted_fragments == 120
    assert urda.copies_per_100_genomes == pytest.approx(40.0)
    assert urda.targets[0].copies_per_100_genomes == pytest.approx(40.0)


def test_length_correction_makes_different_gene_lengths_comparable(database, panel_set):
    """A 600 aa and a 450 aa gene at equal per-genome copy number must agree.

    Equal copy number means fragment counts proportional to reference length:
    600 aa -> 120 fragments, 450 aa -> 90 fragments.
    """
    r1 = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(120)]
    r1 += [hit_line(f"B{i}", 3, "BUT", "P00004") for i in range(90)]
    r1 += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(600)]

    result = _run(database, panel_set, r1)

    urda = result.panel("urda").copies_per_100_genomes
    butyrate = result.panel("butyrate").copies_per_100_genomes
    assert urda == pytest.approx(40.0)
    assert butyrate == pytest.approx(40.0)


def test_depth_independence(database, panel_set):
    """Doubling depth must not change the normalised figure."""

    def build(scale: int) -> list[str]:
        lines = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(60 * scale)]
        lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(300 * scale)]
        return lines

    shallow = _run(database, panel_set, build(1)).panel("urda").copies_per_100_genomes
    deep = _run(database, panel_set, build(2)).panel("urda").copies_per_100_genomes

    assert shallow == pytest.approx(deep)


def test_panel_aggregate_sum_over_genes(database, panel_set):
    """urda aggregates by sum; with one gene the aggregate is that gene."""
    r1 = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(120)]
    r1 += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(600)]

    urda = _run(database, panel_set, r1).panel("urda")

    assert urda.aggregate_method == "sum"
    assert urda.copies_per_100_genomes == pytest.approx(
        sum(t.copies_per_100_genomes or 0.0 for t in urda.targets)
    )


# --------------------------------------------------------------------------- #
# 5. zero rpoB
# --------------------------------------------------------------------------- #


def test_zero_rpob_is_indeterminate_not_a_division_by_zero(database, panel_set):
    r1 = [hit_line(f"U{i}", 0, "URDA", "P00001") for i in range(50)]

    result = _run(database, panel_set, r1)

    assert result.normalizer.fragments == 0
    assert result.normalizer.indeterminate is True
    assert result.normalizer.rate is None

    urda = result.panel("urda")
    assert urda.indeterminate is True
    assert urda.copies_per_100_genomes is None
    assert urda.targets[0].copies_per_100_genomes is None
    # The raw fragment count is still reported.
    assert urda.accepted_fragments == 50


def test_zero_rpob_reports_indeterminate_band():
    from openbiota.report import band_for

    assert band_for(None) == "indeterminate"
    assert band_for(0.0) == "not detected"
    assert band_for(0.5) == "trace"
    assert band_for(5.0) == "low"
    assert band_for(25.0) == "moderate"
    assert band_for(75.0) == "high"
    assert band_for(400.0) == "very high"


def test_rpob_below_identity_floor_does_not_count(database, panel_set):
    r1 = [hit_line(f"R{i}", 4, "RPOB", "P00005", pident=40.0) for i in range(100)]

    result = _run(database, panel_set, r1)

    assert result.normalizer.fragments == 0
    assert result.normalizer.rejected_low_identity == 100
    assert result.normalizer.indeterminate is True


# --------------------------------------------------------------------------- #
# 7. organism attribution reconciles with fragments, not reads
# --------------------------------------------------------------------------- #


def test_organism_counts_sum_to_fragment_total_not_read_total(database, panel_set):
    """This was a real bug in the prototype: attribution counted reads."""
    r1 = [hit_line(f"A{i}", 0, "URDA", "P00001") for i in range(30)]
    r1 += [hit_line(f"B{i}", 1, "URDA", "P00002") for i in range(20)]
    # Every fragment appears again in R2 — 100 read hits, 50 fragments.
    r2 = list(r1)

    result = _run(database, panel_set, r1, r2)
    gene = result.panel("urda").targets[0]

    assert result.total_hit_lines == 100
    assert gene.fragments == 50
    assert sum(o.fragments for o in gene.organisms) == 50
    assert {o.organism for o in gene.organisms} == {
        "Streptococcus pasteurianus",
        "Eggerthella lenta",
    }


def test_organism_attribution_survives_cross_reference_mate_disagreement(database, panel_set):
    """A fragment must be attributed exactly once, to its winning reference."""
    r1 = [hit_line("READ1", 0, "URDA", "P00001", bitscore=80.0)]
    r2 = [hit_line("READ1", 1, "URDA", "P00002", bitscore=140.0)]

    gene = _run(database, panel_set, r1, r2).panel("urda").targets[0]

    assert gene.fragments == 1
    assert sum(o.fragments for o in gene.organisms) == 1
    assert gene.organisms[0].organism == "Eggerthella lenta"


def test_normalizer_organism_attribution_reconciles(database, panel_set):
    r1 = [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(60)]
    r1 += [hit_line(f"S{i}", 5, "RPOB", "P00006") for i in range(40)]

    normalizer = _run(database, panel_set, r1).normalizer

    assert normalizer.fragments == 100
    assert sum(normalizer.organism_counts.values()) == 100
    assert sum(normalizer.phylum_counts.values()) == 100
    assert normalizer.phylum_counts == {"Bacteroidota": 60, "Bacillota": 40}


# --------------------------------------------------------------------------- #
# parser and robustness
# --------------------------------------------------------------------------- #


def test_parser_rejects_field_list_missing_required_columns():
    with pytest.raises(ValueError, match="missing"):
        HitParser(("qseqid", "sseqid", "pident"))


def test_parser_tolerates_malformed_lines(database, panel_set):
    r1 = [
        "not a tsv line at all",
        "",
        "\t".join(["READ1", "0~URDA~P00001", "notanumber", "45", "1", "45", "90", "", ""]),
        hit_line("READ2", 0, "URDA", "P00001"),
    ]

    result = _run(database, panel_set, r1)

    assert result.panel("urda").accepted_fragments == 1


def test_unknown_subject_id_is_counted_not_crashed(database, panel_set):
    r1 = [
        hit_line("READ1", 999, "GHOST", "P99999"),
        "\t".join(["READ2", "notanumber~X~Y", "90", "45", "1", "45", "90", "", ""]),
        hit_line("READ3", 0, "URDA", "P00001"),
    ]

    result = _run(database, panel_set, r1)

    assert result.unresolved_sseqids == 2
    assert result.panel("urda").accepted_fragments == 1


def test_hits_per_mate_recorded(database, panel_set):
    r1 = [hit_line(f"A{i}", 0, "URDA", "P00001") for i in range(7)]
    r2 = [hit_line(f"A{i}", 0, "URDA", "P00001") for i in range(3)]

    result = _run(database, panel_set, r1, r2)

    assert result.hits_per_mate == {"R1": 7, "R2": 3}
