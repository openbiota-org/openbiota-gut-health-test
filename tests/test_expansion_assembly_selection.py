"""Targeted assembly assembles the reads no gut catalogue placed at species.

Pinned: Kraken's unclassified reads and reads whose LCA is above species are
selected; reads at species or below (including strains under a species) are
not; the counts add up to the total so the operating point is auditable.
"""

from __future__ import annotations

import gzip
from pathlib import Path

from openbiota.engines import kraken


def test_unplaced_read_ids(tmp_path: Path) -> None:
    tax = tmp_path / "taxonomy"
    tax.mkdir()
    (tax / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tgenus\t|\n"
        "11\t|\t10\t|\tspecies\t|\n12\t|\t11\t|\tstrain\t|\n")
    (tax / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
    assignments = tmp_path / "a.tsv.gz"
    with gzip.open(assignments, "wt") as fh:
        fh.write("C\tr1\t11\t151|151\tx\nC\tr2\t12\t151|151\tx\nC\tr3\t10\t151|151\tx\n"
                 "U\tr4\t0\t151|151\tx\nC\tr5\t2\t151|151\tx\n")
    ids, counts = kraken.unplaced_read_ids(assignments, tmp_path)
    assert sorted(ids) == ["r3", "r4", "r5"]
    assert counts == {"total": 5, "unclassified": 1, "above_species": 2, "at_species": 2}
    assert counts["unclassified"] + counts["above_species"] + counts["at_species"] == counts["total"]
