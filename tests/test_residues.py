"""Active-site residue mapping and read-level checks."""

from __future__ import annotations

import pytest

from openbiota.residues import (
    NO_MAPPING,
    NOT_SPANNED,
    SPAN_FAIL,
    SPAN_GAP,
    SPAN_PASS,
    AnchorMapping,
    classify_read_residue,
    map_position_through_alignment,
    residue_at_subject_position,
)

ACCEPTED = frozenset({"Y", "M"})


# --------------------------------------------------------------------------- #
# position mapping through a gapped alignment
# --------------------------------------------------------------------------- #


def test_ungapped_alignment_maps_position_directly():
    #             positions 10..14 in both
    query = "ACDEF"
    subject = "ACDEF"
    assert map_position_through_alignment(query, subject, 10, 10, 12) == 12


def test_offset_alignment_maps_correctly():
    """Subject starts at 100, query at 1; subject 103 is the fourth column."""
    assert map_position_through_alignment("ACDEF", "ACDEF", 100, 1, 103) == 4


def test_gap_in_query_shifts_the_subject_frame():
    #  subject: A C D E F   -> positions 100 101 102 103 104
    #  query:   A - D E F   -> positions   1       2   3   4
    assert map_position_through_alignment("A-DEF", "ACDEF", 100, 1, 100) == 1
    # subject 101 aligns to a query gap -> unusable
    assert map_position_through_alignment("A-DEF", "ACDEF", 100, 1, 101) is None
    assert map_position_through_alignment("A-DEF", "ACDEF", 100, 1, 102) == 2


def test_gap_in_subject_does_not_advance_subject_position():
    #  subject: A - C D E   -> positions 100     101 102 103
    #  query:   A B C D E   -> positions   1   2   3   4   5
    assert map_position_through_alignment("ABCDE", "A-CDE", 100, 1, 101) == 3
    assert map_position_through_alignment("ABCDE", "A-CDE", 100, 1, 103) == 5


def test_position_outside_alignment_returns_none():
    assert map_position_through_alignment("ACDEF", "ACDEF", 10, 10, 99) is None


def test_mismatched_alignment_lengths_return_none():
    assert map_position_through_alignment("ACDEF", "ACD", 1, 1, 2) is None


# --------------------------------------------------------------------------- #
# reading the residue off a read alignment
# --------------------------------------------------------------------------- #


def test_residue_at_subject_position():
    #  subject positions 370..375  = A C Y D E F ; 373 -> 'Y' column -> query 'M'
    assert residue_at_subject_position("ACMDEF", "ACYDEF", 370, 372) == "M"
    assert residue_at_subject_position("ACMDEF", "ACYDEF", 370, 370) == "A"
    assert residue_at_subject_position("ACMDEF", "ACYDEF", 370, 999) is None


def test_residue_lookup_skips_subject_gaps():
    #  subject: A - Y D  -> positions 370     371 372
    assert residue_at_subject_position("AQMD", "A-YD", 370, 371) == "M"


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #


def test_classify_pass():
    outcome, residue = classify_read_residue(
        anchor_position=373,
        subject_start=350,
        subject_end=400,
        query_gapped="Y" * 51,
        subject_gapped="Y" * 51,
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (SPAN_PASS, "Y")


def test_classify_accepts_methionine_too():
    outcome, residue = classify_read_residue(
        anchor_position=373,
        subject_start=350,
        subject_end=400,
        query_gapped="M" * 51,
        subject_gapped="Y" * 51,
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (SPAN_PASS, "M")


def test_classify_fail_on_wrong_residue():
    outcome, residue = classify_read_residue(
        anchor_position=373,
        subject_start=350,
        subject_end=400,
        query_gapped="A" * 51,
        subject_gapped="Y" * 51,
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (SPAN_FAIL, "A")


def test_classify_gap_when_read_has_a_deletion():
    length = 51
    query = ["A"] * length
    subject = ["Y"] * length
    query[373 - 350] = "-"
    outcome, residue = classify_read_residue(
        anchor_position=373,
        subject_start=350,
        subject_end=400,
        query_gapped="".join(query),
        subject_gapped="".join(subject),
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (SPAN_GAP, "-")


def test_classify_not_spanned_when_alignment_misses_the_position():
    outcome, residue = classify_read_residue(
        anchor_position=373,
        subject_start=10,
        subject_end=60,
        query_gapped="A" * 51,
        subject_gapped="A" * 51,
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (NOT_SPANNED, None)


def test_classify_no_mapping_when_reference_has_no_anchor_position():
    outcome, residue = classify_read_residue(
        anchor_position=None,
        subject_start=350,
        subject_end=400,
        query_gapped="Y" * 51,
        subject_gapped="Y" * 51,
        accepted=ACCEPTED,
    )
    assert (outcome, residue) == (NO_MAPPING, None)


def test_classify_handles_reversed_subject_coordinates():
    outcome, _ = classify_read_residue(
        anchor_position=373,
        subject_start=400,
        subject_end=350,
        query_gapped="Y" * 51,
        subject_gapped="Y" * 51,
        accepted=ACCEPTED,
    )
    assert outcome == SPAN_PASS


def test_classify_not_spanned_without_alignment_strings():
    """Residue columns are only requested when a panel needs them."""
    outcome, _ = classify_read_residue(
        anchor_position=373,
        subject_start=350,
        subject_end=400,
        query_gapped="",
        subject_gapped="",
        accepted=ACCEPTED,
    )
    assert outcome == NOT_SPANNED


# --------------------------------------------------------------------------- #
# residue results in the tally
# --------------------------------------------------------------------------- #


def test_residue_summary_partitions_every_accepted_fragment(database, panel_set):
    """The five outcome codes must account for all accepted fragments."""
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    lines = []
    # spans the mapped position (300) and passes
    for i in range(5):
        lines.append(
            hit_line(
                f"P{i}", 0, "URDA", "P00001",
                sstart=280, send=330,
                qseq_gapped="Y" * 51, sseq_gapped="Y" * 51,
            )
        )
    # spans and fails
    for i in range(3):
        lines.append(
            hit_line(
                f"F{i}", 0, "URDA", "P00001",
                sstart=280, send=330,
                qseq_gapped="A" * 51, sseq_gapped="Y" * 51,
            )
        )
    # does not span
    for i in range(7):
        lines.append(
            hit_line(
                f"N{i}", 0, "URDA", "P00001",
                sstart=10, send=60,
                qseq_gapped="A" * 51, sseq_gapped="A" * 51,
            )
        )
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(200)]

    result = tally(
        streams=[("R1", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )
    gene = result.panel("urda").targets[0]
    residue = gene.residue

    assert residue is not None
    assert residue.available is True
    assert gene.fragments == 15
    assert residue.fragments_spanning == 8
    assert residue.passed == 5
    assert residue.failed == 3
    assert residue.not_spanned == 7
    assert residue.pass_rate == pytest.approx(5 / 8)
    # the headline count is unchanged by the residue check
    assert result.panel("urda").accepted_fragments == 15


def test_reference_side_filter_splits_fragments_by_reference_residue(database, panel_set):
    """Reference idx 0 carries an accepted residue; idx 1 does not.

    The reference-side half of the filter applies to every accepted fragment,
    unlike the read-side half which needs the read to reach the site.
    """
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    lines = [hit_line(f"OK{i}", 0, "URDA", "P00001") for i in range(30)]
    lines += [hit_line(f"BAD{i}", 1, "URDA", "P00002") for i in range(90)]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(600)]

    result = tally(
        streams=[("R1", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )
    panel = result.panel("urda")
    gene = panel.targets[0]
    residue = gene.residue
    assert residue is not None

    assert gene.fragments == 120
    assert residue.fragments_on_residue_ok_reference == 30
    assert residue.fragments_on_residue_bad_reference == 90
    assert residue.reference_consistent_fraction == pytest.approx(0.25)
    assert residue.references_residue_ok == 1

    #   headline:   120 / 600 aa  vs 600 / 1200 aa -> 40.0 per 100 genomes
    #   consistent:  30 / 600 aa  vs 600 / 1200 aa -> 10.0 per 100 genomes
    assert gene.copies_per_100_genomes == pytest.approx(40.0)
    assert gene.copies_per_100_genomes_residue_consistent == pytest.approx(10.0)
    assert panel.copies_per_100_genomes_residue_consistent == pytest.approx(10.0)
    assert panel.residue_consistent_fragments == 30
    # the tighter estimate never exceeds the upper bound
    assert (
        gene.copies_per_100_genomes_residue_consistent <= gene.copies_per_100_genomes
    )


def test_reference_side_filter_is_indeterminate_without_rpob(database, panel_set):
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    lines = [hit_line(f"OK{i}", 0, "URDA", "P00001") for i in range(30)]
    result = tally(
        streams=[("R1", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )
    gene = result.panel("urda").targets[0]

    assert gene.copies_per_100_genomes is None
    assert gene.copies_per_100_genomes_residue_consistent is None
    # counts are still reported
    assert gene.residue is not None
    assert gene.residue.fragments_on_residue_ok_reference == 30


def test_concordant_residue_checks_are_reported_as_corroboration(database, panel_set):
    """The read-side and reference-side checks share no evidence.

    Reference 0 is residue-consistent, reference 1 is not. Sending 40 fragments
    to reference 0 (all spanning and passing) and 60 to reference 1 (all
    spanning and failing) makes both estimates land on 40.
    """
    from openbiota.report import run_diagnostics
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    lines = [
        hit_line(f"OK{i}", 0, "URDA", "P00001", sstart=280, send=330,
                 qseq_gapped="Y" * 51, sseq_gapped="Y" * 51)
        for i in range(40)
    ]
    lines += [
        hit_line(f"BAD{i}", 1, "URDA", "P00002", sstart=280, send=330,
                 qseq_gapped="A" * 51, sseq_gapped="Y" * 51)
        for i in range(60)
    ]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]
    lines += [hit_line(f"B{i}", 3, "BUT", "P00004") for i in range(400)]

    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )
    residue = result.panel("urda").targets[0].residue
    assert residue is not None
    assert residue.passed == 40
    assert residue.failed == 60
    assert residue.pass_rate == pytest.approx(0.40)
    assert residue.reference_consistent_fraction == pytest.approx(0.40)

    diagnostics = run_diagnostics(tally=result, qc=None, selected=list(result.panels))
    hit = next(d for d in diagnostics if d.code == "residue-concordant-urda:URDA")
    assert hit.level == "info"
    assert "agree" in hit.message
    assert "share no evidence" in hit.message


def test_discordant_residue_checks_warn(database, panel_set):
    """Reference 0 is consistent, so read-side passes but reference-side is 100%."""
    from openbiota.report import run_diagnostics
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    # all 100 fragments land on the residue-consistent reference, but every
    # read that reaches the site carries a rejected residue
    lines = [
        hit_line(f"X{i}", 0, "URDA", "P00001", sstart=280, send=330,
                 qseq_gapped="A" * 51, sseq_gapped="Y" * 51)
        for i in range(100)
    ]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]

    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )
    residue = result.panel("urda").targets[0].residue
    assert residue is not None
    assert residue.pass_rate == pytest.approx(0.0)
    assert residue.reference_consistent_fraction == pytest.approx(1.0)

    diagnostics = run_diagnostics(tally=result, qc=None, selected=list(result.panels))
    hit = next(d for d in diagnostics if d.code == "residue-discordant-urda:URDA")
    assert hit.level == "warn"
    assert "disagree" in hit.message


def test_thin_read_side_subset_is_not_cross_checked(database, panel_set):
    from openbiota.report import run_diagnostics
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    # only 5 fragments reach the site — too few to corroborate anything
    lines = [
        hit_line(f"S{i}", 0, "URDA", "P00001", sstart=280, send=330,
                 qseq_gapped="Y" * 51, sseq_gapped="Y" * 51)
        for i in range(5)
    ]
    lines += [
        hit_line(f"N{i}", 0, "URDA", "P00001", sstart=10, send=60,
                 qseq_gapped="Y" * 51, sseq_gapped="Y" * 51)
        for i in range(95)
    ]
    lines += [hit_line(f"R{i}", 4, "RPOB", "P00005") for i in range(2000)]

    result = tally(
        streams=[("R1", iter(lines))], database=database, panel_set=panel_set, fields=FIELDS
    )
    diagnostics = run_diagnostics(tally=result, qc=None, selected=list(result.panels))
    codes = {d.code for d in diagnostics}

    assert "residue-read-side-thin-urda:URDA" in codes
    assert "residue-concordant-urda:URDA" not in codes
    assert "residue-discordant-urda:URDA" not in codes


def test_residue_check_absent_for_panels_without_one(database, panel_set):
    from openbiota.tally import tally

    from .conftest import FIELDS, hit_line

    lines = [hit_line("A", 3, "BUT", "P00004")]
    result = tally(
        streams=[("R1", iter(lines))],
        database=database,
        panel_set=panel_set,
        fields=FIELDS,
    )
    assert result.panel("butyrate").targets[0].residue is None


# --------------------------------------------------------------------------- #
# AnchorMapping serialisation
# --------------------------------------------------------------------------- #


def test_anchor_mapping_round_trip():
    mapping = AnchorMapping(
        panel="urda",
        anchor_accession="Q8CVD0",
        anchor_length=582,
        canonical_position=373,
        canonical_residue="Y",
        accepted_residues=("M", "Y"),
        positions={0: 300, 7: 412},
        residue_ok={0: True, 7: False},
        n_candidates=10,
        n_aligned=9,
        n_position_spanned=2,
        n_reference_residue_ok=1,
    )
    restored = AnchorMapping.from_json(mapping.to_json())
    assert restored == mapping
    assert restored.positions[7] == 412
    assert restored.residue_ok == {0: True, 7: False}
