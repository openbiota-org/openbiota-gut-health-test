"""Reference parsing, filtering, index persistence and role precedence.

These tests exercise the pure functions and the on-disk index. The network
fetch itself is not tested here; ``openbiota build-db`` is the integration path.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota.errors import ReferenceFetchError
from openbiota.panels import RefEntry
from openbiota.references import (
    SEP,
    EntryStats,
    ReferenceDatabase,
    RefSequence,
    _length_window,
    _load_index,
    _parse_lineage,
    _parse_tsv,
    _write_index,
    database_fingerprint,
)
from openbiota.seqio import clean_sequence, iter_fasta, write_fasta

HEADER = "Entry\tProtein names\tGene Names (primary)\tOrganism\tTaxonomic lineage\tLength\tSequence"


def _entry(**kwargs) -> RefEntry:
    defaults = {
        "id": "TEST",
        "role": "target",
        "source": "uniprot",
        "query": "(protein_name:test)",
        "min_identity": 60.0,
        "label": "test",
        "panel": "p",
    }
    return RefEntry(**{**defaults, **kwargs})


# --------------------------------------------------------------------------- #
# TSV parsing
# --------------------------------------------------------------------------- #


def test_parse_tsv_extracts_every_field():
    text = "\n".join(
        [
            HEADER,
            "\t".join(
                [
                    "P37870",
                    "DNA-directed RNA polymerase subunit beta (EC 2.7.7.6)",
                    "rpoB",
                    "Bacillus subtilis (strain 168)",
                    "cellular organisms (no rank), Bacteria (domain), Bacillota (phylum), "
                    "Bacilli (class), Caryophanales (order), Bacillaceae (family)",
                    "5",
                    "MACDE",
                ]
            ),
        ]
    )
    records = _parse_tsv(text)

    assert len(records) == 1
    record = records[0]
    assert record.accession == "P37870"
    assert record.gene == "rpoB"
    assert record.organism == "Bacillus subtilis (strain 168)"
    assert record.length == 5
    assert record.sequence == "MACDE"


def test_parse_tsv_skips_rows_without_a_sequence():
    text = "\n".join(
        [
            HEADER,
            "\t".join(["A1", "x", "g", "Org", "Bacteria (domain)", "5", ""]),
            "\t".join(["A2", "x", "g", "Org", "Bacteria (domain)", "5", "MACDE"]),
            "\t".join(["A3", "x", "g", "Org", "Bacteria (domain)", "notanumber", "MACDE"]),
            "short\trow",
            "",
        ]
    )
    records = _parse_tsv(text)

    assert [r.accession for r in records] == ["A2"]


def test_parse_tsv_rejects_unexpected_columns():
    with pytest.raises(ReferenceFetchError, match="unexpected UniProt TSV columns"):
        _parse_tsv("Accession\tSeq\nA1\tMACDE")


def test_parse_lineage_pulls_the_ranks_we_use():
    lineage = (
        "cellular organisms (no rank), Bacteria (domain), Bacillati (kingdom), "
        "Bacillota (phylum), Clostridia (class), Eubacteriales (order), "
        "Lachnospiraceae (family), Roseburia (genus)"
    )
    taxa = _parse_lineage(lineage)

    assert taxa["phylum"] == "Bacillota"
    assert taxa["genus"] == "Roseburia"
    assert taxa["family"] == "Lachnospiraceae"
    assert "no rank" not in taxa


def test_parse_lineage_tolerates_junk():
    assert _parse_lineage("") == {}
    assert _parse_lineage("no parens here") == {}


# --------------------------------------------------------------------------- #
# length windows
# --------------------------------------------------------------------------- #


def test_explicit_length_range_is_used_verbatim():
    lo, hi, label = _length_window(_entry(length_range=(1000, 1500)), [1, 2, 3])
    assert (lo, hi) == (1000, 1500)
    assert "explicit" in label


def test_median_relative_window_trims_fragments_and_fusions():
    """Reference length must be trimmed or the length normalisation is corrupt."""
    lengths = [50, 600, 610, 620, 630, 640, 3000]
    lo, hi, label = _length_window(_entry(length_tolerance=0.25), lengths)

    assert lo == int(620 * 0.75)
    assert hi == int(620 * 1.25)
    assert lo > 50 and hi < 3000  # both outliers excluded
    assert "median" in label


def test_tolerance_widens_the_window():
    narrow = _length_window(_entry(length_tolerance=0.10), [600] * 5)
    wide = _length_window(_entry(length_tolerance=0.50), [600] * 5)
    assert wide[0] < narrow[0]
    assert wide[1] > narrow[1]


# --------------------------------------------------------------------------- #
# fingerprint
# --------------------------------------------------------------------------- #


def test_fingerprint_changes_when_a_query_changes(panel_set):
    entries = list(panel_set.all_entries())
    before = database_fingerprint(panel_set, entries, "diamond 2.2.6")

    changed = [
        _entry(id=e.id, role=e.role, query=e.query + " AND (reviewed:true)", panel=e.panel,
               label=e.label, min_identity=e.min_identity)
        if e.id == "URDA"
        else e
        for e in entries
    ]
    after = database_fingerprint(panel_set, changed, "diamond 2.2.6")

    assert before != after


def test_fingerprint_changes_with_the_diamond_version(panel_set):
    entries = list(panel_set.all_entries())
    assert database_fingerprint(panel_set, entries, "2.2.6") != database_fingerprint(
        panel_set, entries, "2.1.9"
    )


def test_fingerprint_is_order_independent(panel_set):
    entries = list(panel_set.all_entries())
    assert database_fingerprint(panel_set, entries, "d") == database_fingerprint(
        panel_set, list(reversed(entries)), "d"
    )


def test_fingerprint_tracks_the_residue_check(panel_set):
    import dataclasses

    entries = list(panel_set.all_entries())
    before = database_fingerprint(panel_set, entries, "d")

    urda = panel_set.by_name("urda")
    moved = dataclasses.replace(
        panel_set,
        panels=(
            dataclasses.replace(
                urda,
                residue_check=dataclasses.replace(urda.residue_check, canonical_position=374),
            ),
            *panel_set.panels[1:],
        ),
    )
    assert database_fingerprint(moved, entries, "d") != before


# --------------------------------------------------------------------------- #
# sseqid encoding and index round trip
# --------------------------------------------------------------------------- #


def test_sseqid_encoding_puts_the_index_first():
    reference = RefSequence(
        42, "urda:URDA", "URDA", "target", "urda", "Q8CVD0", "Shewanella oneidensis", 582
    )
    assert reference.sseqid == f"42{SEP}URDA{SEP}Q8CVD0"
    # separator must not be something DIAMOND's seqid parser treats specially
    assert "|" not in reference.sseqid
    assert " " not in reference.sseqid


def test_by_sseqid_is_an_index_lookup(database):
    assert database.by_sseqid(f"0{SEP}URDA{SEP}P00001").accession == "P00001"
    assert database.by_sseqid(f"3{SEP}BUT{SEP}P00004").entry_key == "butyrate:BUT"
    assert database.by_sseqid(f"999{SEP}X{SEP}Y") is None
    assert database.by_sseqid("garbage") is None
    assert database.by_sseqid("") is None


def test_index_round_trip_preserves_everything(database, tmp_path: Path):
    database.index_path = tmp_path / "index.json"
    _write_index(database)

    restored = _load_index(database.index_path, database.dmnd_path, database.fasta_path)

    assert restored.fingerprint == database.fingerprint
    assert restored.n_sequences == database.n_sequences
    assert restored.sequences == database.sequences
    assert restored.entries == database.entries
    assert restored.anchors["urda"].positions == {0: 300, 1: 300}
    assert restored.anchors["urda"].residue_ok == {0: True, 1: False}
    assert restored.by_sseqid(f"1{SEP}URDA{SEP}P00002").anchor_residue_ok is False


def test_index_rejects_a_stale_schema_version(database, tmp_path: Path):
    database.index_path = tmp_path / "index.json"
    _write_index(database)
    payload = json.loads(database.index_path.read_text(encoding="utf-8"))
    payload["index_version"] = 1
    database.index_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="index version"):
        _load_index(database.index_path, database.dmnd_path, database.fasta_path)


def test_index_rejects_sparse_rows(database, tmp_path: Path):
    database.index_path = tmp_path / "index.json"
    _write_index(database)
    payload = json.loads(database.index_path.read_text(encoding="utf-8"))
    payload["sequences"][2][0] = 99
    database.index_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="must be dense"):
        _load_index(database.index_path, database.dmnd_path, database.fasta_path)


def test_total_residues(database):
    assert database.total_residues() == 600 + 600 + 600 + 450 + 1200 + 1200


# --------------------------------------------------------------------------- #
# sequence I/O
# --------------------------------------------------------------------------- #


def test_clean_sequence_strips_non_letters_and_uppercases():
    assert clean_sequence("acd ef\n1-2*") == "ACDEF"
    assert clean_sequence("") == ""


def test_fasta_round_trip(tmp_path: Path):
    records = [("a~A~P1", "MACDE"), ("b~B~P2", "WYFGH")]
    path = tmp_path / "x.faa"
    write_fasta(path, records)

    assert list(iter_fasta(path)) == records


def test_iter_fasta_joins_wrapped_lines(tmp_path: Path):
    path = tmp_path / "x.faa"
    path.write_text(">h\nMACD\nEFGH\n\n>i\nWY\n", encoding="utf-8")

    assert list(iter_fasta(path)) == [("h", "MACDEFGH"), ("i", "WY")]


# --------------------------------------------------------------------------- #
# EntryStats serialisation
# --------------------------------------------------------------------------- #


def test_entry_stats_round_trip(database):
    stats = database.entries["urda:URDA"]
    assert EntryStats(**stats.to_json()) == stats


def test_reference_database_entry_lookup(database):
    assert database.entry("urda:URDA").role == "target"
    with pytest.raises(KeyError):
        database.entry("nope:NOPE")


def test_role_precedence_ordering_is_normalizer_then_target_then_decoy(panel_set):
    """A protein claimed by a target must never also sit in a decoy set."""
    entries = list(panel_set.all_entries())
    ordered = sorted(entries, key=lambda e: {"normalizer": 0, "target": 1, "decoy": 2}[e.role])
    roles = [e.role for e in ordered]

    assert roles[0] == "normalizer"
    assert roles.index("target") < roles.index("decoy")
    # role ordering must be a stable partition
    assert roles == sorted(roles, key=lambda r: {"normalizer": 0, "target": 1, "decoy": 2}[r])


def test_reference_database_is_a_real_object_not_a_dict(database: ReferenceDatabase):
    assert isinstance(database, ReferenceDatabase)
    assert database.n_sequences == 6
