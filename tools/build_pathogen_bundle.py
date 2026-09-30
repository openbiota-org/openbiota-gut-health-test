#!/usr/bin/env python3
"""Rebuild the pathogen reference bundle from the installed catalog.

    .venv/bin/python tools/build_pathogen_bundle.py --threads 32

Every stage is incremental where it can be:

* **Resolution** reuses `refs/pathogens/resolutions.cache.json`, so only a seed
  that is new or has changed hits NCBI.
* **Acquisition** skips any asset that already has an `asset.json` stamp, so
  only new genomes are downloaded.
* **The k-mer ownership pass is not incremental, and must not be.** Ownership
  is decided across the whole pooled k-mer set: which sequence is specific to
  one organism depends on every other organism installed. Adding a reference
  therefore has to re-decide ownership for all of them — that is the point of
  adding it.
* **The competitive bowtie2 index** is keyed by its member set, so a new member
  means a new index. The old one stays on disk until `openbiota prune` removes
  it, which is what makes the switch atomic.

Why this exists: the bundle had no commensal *Escherichia coli*. Ordinary gut
E. coli reads had nowhere correct to land, so competitive alignment placed them
on the nearest reference — a *Shigella* genome, which is the same genomic
species — and the ownership pass never masked the sequence the two share. The
result was a report naming three *Shigella* species in a sample whose only
Enterobacteriaceae was E. coli.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # pragma: no cover - script entry
    sys.path.insert(0, str(ROOT))

from openbiota.pathogens.catalog import load_catalog  # noqa: E402
from openbiota.pathogens.refs import (  # noqa: E402
    acquire,
    build_bundle,
    deduplicate_references,
    resolve_catalog,
)

#: The background genomes whose k-mers mask target-specific sequence. Human
#: plus the food animals a stool sample plausibly carries.
DECOY_NAMES = ("bos_taurus", "gallus_gallus", "sus_scrofa")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs", default="refs", help="reference root (default: refs)")
    ap.add_argument("--catalog", default="pathogens/catalog")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument(
        "--skip-aligner-index",
        action="store_true",
        help="write the k-mer bundle but leave the bowtie2 index to the next run",
    )
    args = ap.parse_args(argv)

    refs = Path(args.refs)
    pathogens = refs / "pathogens"
    started = time.monotonic()

    def say(message: str) -> None:
        print(f"[{time.monotonic() - started:8.1f}s] {message}", flush=True)

    catalog = load_catalog(Path(args.catalog))
    say(f"catalog {catalog.version} — {len(catalog.seeds)} seeds")

    resolutions = resolve_catalog(
        catalog,
        cache=pathogens / "resolutions.cache.json",
        progress=say,
    )
    resolved = sum(1 for r in resolutions.values() if r.resolution_state == "resolved")
    say(f"resolved {resolved}/{len(resolutions)}")

    folded = deduplicate_references(catalog, resolutions, progress=say)
    if folded:
        say(f"folded {folded} row(s) onto a more specific row's reference")

    assets = acquire(resolutions, pathogens / "assets", progress=say)
    say(f"{len(assets)} assets on disk")

    decoys = [refs / "decoys" / name / "genome.fna" for name in DECOY_NAMES]
    host = refs / "host" / "grch38.fna.gz"
    decoy_fastas = [p for p in [host, *decoys] if p.is_file()]
    missing = [p for p in [host, *decoys] if not p.is_file()]
    for p in missing:
        say(f"WARNING: decoy missing, shared sequence will not be masked against it: {p}")

    manifest = build_bundle(
        catalog,
        resolutions,
        assets,
        out_dir=pathogens / "bundle",
        decoy_fastas=decoy_fastas,
        threads=args.threads,
        progress=say,
        build_aligner_index=not args.skip_aligner_index,
        decoy_cache=pathogens / "decoy_kmers.npy",
    )
    say(
        f"bundle {manifest.bundle_id} — catalog {manifest.catalog_version}, "
        f"{len(manifest.entries)} targets"
    )
    say("done. Re-run the screen for every sample: the bundle ID has changed.")
    return 0


if __name__ == "__main__":  # pragma: no cover - script entry
    raise SystemExit(main())
