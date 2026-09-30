#!/usr/bin/env python3
"""Delimit the CTnPc accessory region by presence/absence, not by annotation.

The published supplement identifies CTnPc by 42 RAST CDS identifiers
(`fig|165179.43.peg.373-414`), which are annotation-server IDs and do not map
onto NCBI's annotation of the public assembly. Picking the region from generic
mobile-element annotation instead is unsafe: *S. copri* devotes ~20% of its
CDSs to conjugative machinery, and how many "~100-kb conjugative regions" you
count depends entirely on an arbitrary clustering distance.

So this takes the other route, which is also the better one scientifically: the
element is defined by the property the study actually associated with
arthritis-promoting activity — **presence in the CTnPc-positive isolates and
absence from the negatives**.

Method:

1. Build a DIAMOND protein database from the two CTnPc-negative isolates
   (N115-17, an RA-origin isolate, and H012_6, a healthy control).
2. Query each positive isolate's proteins against it.
3. Proteins with no acceptable hit are candidate accessory content.
4. Map those back to genomic coordinates and report contiguous runs.
5. Intersect the two positives: content specific to one isolate is that
   isolate's, not the shared element.

What this does and does not establish is stated in the output. A region
delimited this way is a **candidate** consistent with the study's definition;
it is not a validated assay, and a marker-positive sample is not "CTnPc
present" without breadth, depth and carrier evidence.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

PANELS = Path("refs/panels/ctnpc")
OUT = Path("refs/ctnpc/CTNPC_REGION.json")

POSITIVES = ("ctnpc.ra_n001_13", "ctnpc.rap9_13")
NEGATIVES = ("ctnpc.n115_17", "ctnpc.h012_6")

#: A protein counts as "present in the negatives" at or above this identity
#: over this much of its length. Deliberately permissive: the question is
#: whether a homolog exists at all, so a loose threshold makes the
#: accessory call conservative.
MIN_IDENTITY = 70.0
MIN_QUERY_COVERAGE = 70.0

#: Accessory proteins this far apart (bp) still count as one run.
MAX_GAP = 12_000

#: A run needs at least this many accessory CDSs to be worth reporting.
MIN_RUN_CDS = 8


@dataclass(frozen=True, slots=True)
class Cds:
    protein_id: str
    contig: str
    start: int
    end: int
    product: str


def _decompress(path: Path, suffix: str, work: Path) -> Path:
    out = work / (path.stem + suffix)
    out.write_bytes(gzip.decompress(path.read_bytes()))
    return out


def load_cds(gff_gz: Path) -> dict[str, Cds]:
    """protein_id -> its coordinates and product, from a GFF."""
    out: dict[str, Cds] = {}
    with gzip.open(gff_gz, "rt") as handle:
        lines = handle.readlines()
    for raw in lines:
        if raw.startswith("#"):
            continue
        parts = raw.rstrip().split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        attrs = parts[8]
        pid = re.search(r"protein_id=([^;]+)", attrs)
        if not pid:
            continue
        product = re.search(r"product=([^;]+)", attrs)
        out[pid.group(1)] = Cds(
            protein_id=pid.group(1),
            contig=parts[0],
            start=int(parts[3]),
            end=int(parts[4]),
            product=product.group(1) if product else "",
        )
    return out


def diamond_hits(query: Path, db: Path, threads: int) -> set[str]:
    """Query protein IDs with an acceptable homolog in the database."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as handle:
        out_path = Path(handle.name)
    subprocess.run(
        ["diamond", "blastp", "--quiet", "--threads", str(threads),
         "--db", str(db), "--query", str(query), "--out", str(out_path),
         "--outfmt", "6", "qseqid", "pident", "qcovhsp",
         "--max-target-seqs", "1", "--evalue", "1e-5"],
        check=True, capture_output=True,
    )
    found: set[str] = set()
    for line in out_path.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        if float(cols[1]) >= MIN_IDENTITY and float(cols[2]) >= MIN_QUERY_COVERAGE:
            found.add(cols[0])
    out_path.unlink(missing_ok=True)
    return found


def runs_of(accessory: list[Cds]) -> list[dict[str, Any]]:
    """Contiguous genomic runs of accessory CDSs."""
    by_contig: dict[str, list[Cds]] = defaultdict(list)
    for cds in accessory:
        by_contig[cds.contig].append(cds)
    out: list[dict[str, Any]] = []
    for contig, items in by_contig.items():
        items.sort(key=lambda c: c.start)
        current: list[Cds] = []
        for cds in items:
            if current and cds.start - current[-1].end > MAX_GAP:
                out.append(_run(contig, current))
                current = []
            current.append(cds)
        if current:
            out.append(_run(contig, current))
    return sorted(
        (r for r in out if r["accessory_cds"] >= MIN_RUN_CDS),
        key=lambda r: -r["span_bp"],
    )


def _run(contig: str, items: list[Cds]) -> dict[str, Any]:
    mobile = re.compile(
        r"conjugat|relaxase|integrase|mobilization|TraG|TraM|TraN|TraK|type IV secret",
        re.I,
    )
    return {
        "contig": contig,
        "start": items[0].start,
        "end": items[-1].end,
        "span_bp": items[-1].end - items[0].start + 1,
        "accessory_cds": len(items),
        "mobile_element_cds": sum(1 for c in items if mobile.search(c.product)),
        "protein_ids": [c.protein_id for c in items],
        "products": sorted({c.product for c in items if c.product})[:24],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args(argv)

    for key in (*POSITIVES, *NEGATIVES):
        for suffix in (".protein.faa.gz", ".gff.gz"):
            if not (PANELS / f"{key}{suffix}").exists():
                print(f"missing {key}{suffix}; run tools/fetch_strain_refs.py --group ctnpc")
                return 1

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        # One database from both negatives: a protein must be absent from
        # *both* to count as accessory.
        negative_faa = work / "negatives.faa"
        with negative_faa.open("wb") as handle:
            for key in NEGATIVES:
                handle.write(gzip.decompress((PANELS / f"{key}.protein.faa.gz").read_bytes()))
        db = work / "negatives.dmnd"
        subprocess.run(
            ["diamond", "makedb", "--quiet", "--in", str(negative_faa), "--db", str(db)],
            check=True, capture_output=True,
        )

        per_positive: dict[str, Any] = {}
        accessory_products: dict[str, set[str]] = {}
        for key in POSITIVES:
            faa = _decompress(PANELS / f"{key}.protein.faa.gz", ".faa", work)
            cds = load_cds(PANELS / f"{key}.gff.gz")
            shared = diamond_hits(faa, db, args.threads)
            accessory = [c for pid, c in cds.items() if pid not in shared]
            found_runs = runs_of(accessory)
            per_positive[key] = {
                "total_cds": len(cds),
                "shared_with_negatives": len(shared),
                "accessory_cds": len(accessory),
                "accessory_fraction": round(len(accessory) / max(1, len(cds)), 4),
                "runs": found_runs,
            }
            accessory_products[key] = {
                c.product for c in accessory if c.product and "hypothetical" not in c.product
            }
            print(f"{key}: {len(cds)} CDS, {len(accessory)} accessory "
                  f"({100 * len(accessory) / max(1, len(cds)):.1f}%), "
                  f"{len(found_runs)} run(s) >= {MIN_RUN_CDS} CDS")
            for run in found_runs[:6]:
                print(f"    {run['contig']} {run['start']:>9,}-{run['end']:<9,} "
                      f"{run['span_bp']:>7,} bp  {run['accessory_cds']:3d} accessory CDS "
                      f"({run['mobile_element_cds']} mobile)")

    shared_products = set.intersection(*accessory_products.values())
    result = {
        "target_id": "ra.ctnpc",
        "method": "protein presence/absence against CTnPc-negative isolates (DIAMOND blastp)",
        "positives": list(POSITIVES),
        "negatives": list(NEGATIVES),
        "thresholds": {
            "min_identity_percent": MIN_IDENTITY,
            "min_query_coverage_percent": MIN_QUERY_COVERAGE,
            "max_gap_bp": MAX_GAP,
            "min_run_cds": MIN_RUN_CDS,
            "rationale": (
                "Permissive homology thresholds on purpose: the question is whether any "
                "homolog exists in the negatives, so a loose threshold makes the accessory "
                "call conservative rather than inflated."
            ),
        },
        "per_positive": per_positive,
        "shared_accessory_products": sorted(shared_products)[:80],
        "n_shared_accessory_products": len(shared_products),
        "status": "candidate_region_delimited",
        "what_this_establishes": (
            "Genomic runs present in both CTnPc-positive isolates and absent from both "
            "negatives, one of which (N115-17) is itself RA-origin. This is the study's own "
            "definition of the element applied to public sequence, rather than a guess from "
            "generic mobile-element annotation."
        ),
        "what_this_does_not_establish": (
            "It is not a validated assay and not a human risk marker. Two positive isolates "
            "cannot separate the element from other content those two happen to share, the "
            "runs are delimited by annotation boundaries rather than by the element's real "
            "termini, and a read-level hit against these regions still requires breadth, "
            "depth and carrier evidence before it is reported as anything."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"\nshared accessory products (in both positives, neither negative): "
          f"{len(shared_products)}")
    print(f"written: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
