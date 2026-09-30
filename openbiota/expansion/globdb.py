"""GlobDB r232: representatives, lineages, source identities, quality.

GlobDB is one clustering over 26 source datasets under its own rules (96%
ANI / 50% aligned fraction, GTDB priority). Its identities are preserved
here exactly as published: a `MGYG000001338` representative is a GlobDB
unit whose name is `Blautia_A MGYG000001338`, not a GTDB species and not
"a new species" until competitive confirmation and the crosswalk say what
it is relative to GTDB R232.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

GLOBDB_DIR: Final = Path("refs/globdb_r232")
TAXONOMY: Final = GLOBDB_DIR / "globdb_r232_taxonomy.tsv.gz"
CHECKM2: Final = GLOBDB_DIR / "globdb_r232_checkm2.tsv.gz"
DICTS: Final = GLOBDB_DIR / "globdb_r232_dictionaries"
DATASETS: Final = GLOBDB_DIR / "globdb_r232_dataset_list.tsv"

#: Identifier prefix -> source dataset, from the dataset list.
SOURCE_OF_PREFIX: Final[dict[str, str]] = {
    "GCA": "GTDB", "GCF": "GTDB", "MOTU40": "mOTU", "SPIREOTU": "SPIRE", "GCMETA": "GCMETA", "GWH": "NGDC",
    "BCRBG": "RBG", "TPMCOTU": "TPMC", "GEMOTU": "GEM", "MGYG": "MGnify", "TPMCS": "TPMCS", "GOMCOTU": "GOMC",
    "SMAGOTU": "SMAG", "QXLSG": "QXLSG", "CRBC": "CRBC", "HOGU": "HOGU", "TPLM": "TPLM", "HRGMV2": "HRGM2",
    "AMXMAG": "AMXMAG", "TG2G": "TG2G", "PREC": "PREC", "FMDMAG": "cFMD", "DAWW": "DAWW", "CRLG": "CRLG",
    "MRGM": "MRGM", "SHGOMAG": "SHGO", "SCSSF": "SCSSF",
}
#: Sources whose genomes are food or non-human-gut habitats: competitors and
#: background for the rescue panel, never gut organisms by default.
BACKGROUND_SOURCES: Final = frozenset({"cFMD", "DAWW", "CRLG", "GOMC", "SMAG", "TPMC", "TPMCS", "TPLM", "QXLSG",
                                       "AMXMAG", "TG2G", "PREC", "SCSSF", "CRBC", "HOGU", "MRGM", "SHGO"})


@dataclass(frozen=True)
class Representative:
    globdb_id: str
    lineage: str
    source: str

    @property
    def species(self) -> str:
        for seg in self.lineage.split(";"):
            if seg.startswith("s__"):
                return seg[3:].strip()
        return ""

    @property
    def genus(self) -> str:
        for seg in self.lineage.split(";"):
            if seg.startswith("g__"):
                return seg[3:].strip()
        return ""

    @property
    def is_gtdb(self) -> bool:
        return self.source == "GTDB"


def source_of(globdb_id: str) -> str:
    head = globdb_id.split("_", 1)[0]
    if head in ("GCA", "GCF"):
        return "GTDB"
    for prefix in sorted(SOURCE_OF_PREFIX, key=len, reverse=True):
        if globdb_id.startswith(prefix):
            return SOURCE_OF_PREFIX[prefix]
    return "unknown"


def representatives() -> dict[str, Representative]:
    out: dict[str, Representative] = {}
    with gzip.open(TAXONOMY, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                gid = parts[0].strip()
                out[gid] = Representative(gid, parts[1].strip(), source_of(gid))
    return out


def quality() -> dict[str, tuple[float | None, float | None]]:
    """globdb id -> (completeness, contamination) from the publisher's CheckM2 table."""
    out: dict[str, tuple[float | None, float | None]] = {}
    if not CHECKM2.is_file():
        return out
    with gzip.open(CHECKM2, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            name = (row.get("Name") or row.get("name") or row.get("genome") or "").strip()
            name = name.rsplit("/", 1)[-1]
            for suf in (".fa.gz", ".fna.gz", ".fa", ".fna"):
                if name.endswith(suf):
                    name = name[: -len(suf)]

            def _f(k: str, row: Mapping[str, Any] = row) -> float | None:
                try:
                    return float(row.get(k) or row.get(k.lower()) or "")
                except ValueError:
                    return None
            out[name] = (_f("Completeness"), _f("Contamination"))
    return out


def dictionaries() -> dict[str, dict[str, str]]:
    """source dataset -> {source genome id -> GlobDB id}, where a dictionary exists."""
    out: dict[str, dict[str, str]] = {}
    if not DICTS.is_dir():
        return out
    for path in sorted(DICTS.glob("*_dictionary.tsv")):
        source = path.name.replace("_dictionary.tsv", "")
        table: dict[str, str] = {}
        with path.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and not parts[0].lower().startswith(("genome", "cluster", "species")):
                    table[parts[0]] = parts[-1]
        out[source] = table
    return out


def summary() -> Mapping[str, object]:
    reps = representatives()
    by_source: dict[str, int] = {}
    for r in reps.values():
        by_source[r.source] = by_source.get(r.source, 0) + 1
    return {"n_representatives": len(reps), "by_source": dict(sorted(by_source.items(), key=lambda kv: -kv[1])),
            "n_non_gtdb": sum(v for k, v in by_source.items() if k != "GTDB")}


if __name__ == "__main__":  # pragma: no cover
    print(json.dumps(summary(), indent=1))
