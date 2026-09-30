"""MetaPhlAn 4 extended catalogue: parsing, SGB rows, and the bridge to the scoring lane."""

from __future__ import annotations

import gzip

from openbiota.engines.metaphlan import _gzip_intact
from openbiota.engines.metaphlan4 import RENAMED, ExtendedCatalogue, parse_profile4

PROFILE = """\
#mpa_vJun23_CHOCOPhlAnSGB_202403
#/path/metaphlan reads.fastq --nproc 8
#7949610 reads processed
#SampleID\tMetaphlan_Analysis
#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species
UNCLASSIFIED\t-1\t38.5\t
k__Bacteria\t2\t60.0\t
k__Archaea\t2157\t1.5\t
k__Bacteria|p__Firmicutes\t2|1239\t45.0\t
k__Bacteria|p__Bacteroidota\t2|976\t15.0\t
k__Archaea|p__Euryarchaeota\t2157|28890\t1.5\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium\t2|1239|186801|186802|216572|216851\t20.0\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii\t2|1239|186801|186802|216572|216851|853\t20.0\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__SGB15318\t2|1239|186801|186802|216572|216851|853|\t20.0\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__GGB9758\t2|1239|186801|186802|186803|\t12.0\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__GGB9758|s__GGB9758_SGB15368\t2|1239|186801|186802|186803||\t12.0\t
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__GGB9758|s__GGB9758_SGB15368|t__SGB15368\t2|1239|186801|186802|186803|||\t12.0\t
k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola\t2|976|200643|171549|815|909656\t15.0\t
k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus\t2|976|200643|171549|815|909656|821\t15.0\t
k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Phocaeicola|s__Phocaeicola_vulgatus|t__SGB1814\t2|976|200643|171549|815|909656|821|\t15.0\t
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter\t2157|28890|183925|2158|2159|2172\t1.5\t
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii\t2157|28890|183925|2158|2159|2172|2173\t1.5\t
k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii|t__SGB714\t2157|28890|183925|2158|2159|2172|2173|\t1.5\t
"""


def _catalogue(tmp_path) -> ExtendedCatalogue:
    path = tmp_path / "mpa4.tsv"
    path.write_text(PROFILE)
    clades, unclassified, n_reads = parse_profile4(path)
    return ExtendedCatalogue(
        index="mpa_vJun23_CHOCOPhlAnSGB_202403", profiler_version="4.1.1", clades=clades,
        unclassified_percent=unclassified, n_reads_processed=n_reads, elapsed_s=1.0, cached=False,
    )


def test_parse_reads_unclassified_and_read_count(tmp_path):
    cat = _catalogue(tmp_path)
    assert cat.unclassified_percent == 38.5
    assert cat.n_reads_processed == 7_949_610
    assert cat.species["Faecalibacterium_prausnitzii"] == 20.0
    assert cat.kingdoms == {"Bacteria": 60.0, "Archaea": 1.5}
    assert cat.phyla["Firmicutes"] == 45.0


def test_sgb_rows_carry_lineage_and_flag_unnamed_bins(tmp_path):
    cat = _catalogue(tmp_path)
    rows = cat.rows()
    assert [r.sgb for r in rows] == ["SGB15318", "SGB1814", "SGB15368", "SGB714"]
    named = {r.sgb: r for r in rows}
    assert named["SGB15368"].unnamed and named["SGB15368"].genus == "GGB9758"
    assert not named["SGB15318"].unnamed
    assert named["SGB714"].kingdom == "Archaea"
    assert cat.n_sgbs == 4 and cat.n_species == 4 and cat.n_genera == 4
    assert cat.n_unnamed == 1 and cat.unnamed_percent == 12.0
    assert cat.kingdom_counts() == {"Bacteria": 3, "Archaea": 1}


def test_bridge_translates_renamed_species_before_comparing(tmp_path):
    cat = _catalogue(tmp_path)
    # MetaPhlAn 3 still calls P. vulgatus "Bacteroides vulgatus".
    assert RENAMED["Bacteroides_vulgatus"] == "Phocaeicola_vulgatus"
    mpa3 = {"Faecalibacterium_prausnitzii": 18.0, "Bacteroides_vulgatus": 14.0, "Roseburia_hominis": 2.0}
    bridge = cat.bridge(mpa3)
    assert bridge["named_by_both"] == 2
    assert bridge["named_only_by_metaphlan3"] == 1
    assert bridge["metaphlan3_only_examples"] == ["Roseburia_hominis"]
    # The uSGB and the archaeon are only nameable by the extended catalogue.
    assert bridge["named_only_by_metaphlan4"] == 2
    assert bridge["metaphlan4_only_examples"][0] == "GGB9758_SGB15368"


def test_to_json_is_explicit_that_this_lane_is_not_scored(tmp_path):
    payload = _catalogue(tmp_path).to_json(mpa3_species={"Faecalibacterium_prausnitzii": 18.0})
    assert payload["engine"] == "metaphlan4"
    assert payload["n_sgbs_detected"] == 4
    assert payload["n_unnamed_sgbs"] == 1
    assert payload["bridge_to_scoring_lane"]["named_by_both"] == 1
    assert "MetaPhlAn 3 lane" in payload["what_this_is"]


# --------------------------------------------------------------------------- #
# host-filter output integrity (shared by both taxonomic lanes)
# --------------------------------------------------------------------------- #


def test_gzip_integrity_rejects_trailing_garbage_and_empty_files(tmp_path):
    good = tmp_path / "good.fastq.gz"
    with gzip.open(good, "wb") as fh:
        fh.write(b"@r1\nACGT\n+\nIIII\n" * 1000)
    assert _gzip_intact(good)

    # Bowtie2's writer has been seen to leave non-gzip bytes after the last
    # member; Python's reader (and MetaPhlAn's) reject exactly this.
    bad = tmp_path / "bad.fastq.gz"
    bad.write_bytes(good.read_bytes() + b"Y\xc6garbage")
    assert not _gzip_intact(bad)

    truncated = tmp_path / "trunc.fastq.gz"
    truncated.write_bytes(good.read_bytes()[:-20])
    assert not _gzip_intact(truncated)

    empty = tmp_path / "empty.fastq.gz"
    empty.write_bytes(b"")
    assert not _gzip_intact(empty)
    assert not _gzip_intact(tmp_path / "missing.fastq.gz")
