"""Genome-wide ANI with alignment fraction, via skani.

Spec 0.8.4 §3.2: unmatched representatives are compared by ANI *and*
alignment fraction, both directions kept. The project rule for a
provisional non-GTDB cluster is 95% ANI and at least 65% aligned fraction
of the shorter genome; boundary and incomplete-genome cases are flagged,
not forced. No transitive chaining: a query is compared to
representatives, and matched or not, one pair at a time.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

ANI_SAME_SPECIES: Final = 95.0
AF_SAME_SPECIES: Final = 65.0
ANI_BOUNDARY: Final = (94.0, 96.0)
AF_BOUNDARY: Final = (55.0, 65.0)


@dataclass(frozen=True)
class Pair:
    query: str
    reference: str
    ani: float
    af_query: float
    af_reference: float

    @property
    def af_shorter(self) -> float:
        """The aligned fraction of the shorter genome is the larger of the two."""
        return max(self.af_query, self.af_reference)

    @property
    def verdict(self) -> str:
        if self.ani >= ANI_SAME_SPECIES and self.af_shorter >= AF_SAME_SPECIES:
            return "same_species"
        if ANI_BOUNDARY[0] <= self.ani <= ANI_BOUNDARY[1] or (
            self.ani >= ANI_SAME_SPECIES and AF_BOUNDARY[0] <= self.af_shorter < AF_BOUNDARY[1]
        ):
            return "boundary"
        return "different"


def skani_version() -> str:
    exe = shutil.which("skani")
    if not exe:
        return ""
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, check=False)
    return (out.stdout or out.stderr).strip()


def available() -> bool:
    return shutil.which("skani") is not None


def _list_file(paths: Iterable[Path], where: Path) -> Path:
    where.parent.mkdir(parents=True, exist_ok=True)
    where.write_text("\n".join(str(p) for p in paths) + "\n")
    return where


def dist(
    queries: Sequence[Path], references: Sequence[Path], *, work_dir: Path, threads: int = 16,
    min_af: float = 15.0,
) -> list[Pair]:
    """skani dist, queries against references; every pair above skani's floor."""
    exe = shutil.which("skani")
    if exe is None:
        raise RuntimeError("skani is not installed")
    work_dir.mkdir(parents=True, exist_ok=True)
    ql = _list_file(queries, work_dir / "queries.txt")
    rl = _list_file(references, work_dir / "references.txt")
    out = work_dir / "skani_dist.tsv"
    cmd = [exe, "dist", "--ql", str(ql), "--rl", str(rl), "-o", str(out), "-t", str(threads),
           "--min-af", str(min_af)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"skani dist failed: {proc.stderr[-1500:]}")
    return parse(out)


def parse(path: Path) -> list[Pair]:
    pairs: list[Pair] = []
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                pairs.append(Pair(
                    query=Path(row["Query_file"]).name, reference=Path(row["Ref_file"]).name,
                    ani=float(row["ANI"]), af_query=float(row["Align_fraction_query"]),
                    af_reference=float(row["Align_fraction_ref"]),
                ))
            except (KeyError, ValueError):
                continue
    return pairs


def best_by_query(pairs: Iterable[Pair]) -> dict[str, Pair]:
    best: dict[str, Pair] = {}
    for p in pairs:
        cur = best.get(p.query)
        if cur is None or (p.ani, p.af_shorter) > (cur.ani, cur.af_shorter):
            best[p.query] = p
    return best
