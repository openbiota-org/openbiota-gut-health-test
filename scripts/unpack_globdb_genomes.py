#!/usr/bin/env python3
"""Unpack the GlobDB r232 genome archive into refs/globdb_r232/genomes/.

GlobDB publishes its 346,233 representative genomes as one 272 GB tar of
per-genome `.fa.gz` files in chunk directories. Competitive confirmation of
a GlobDB-native cluster needs its genome, and `openbiota.expansion.genomes`
looks for `<id>.fa.gz` anywhere under refs/globdb_r232/genomes/, so the
archive's own layout is kept.

Resumable: tar's -k keeps files already on disk, so a second run after an
interruption only writes what is missing. A lock (`genomes.unpack.json`)
records the archive's SHA-256 from the fetch lock, the file count and the
time, and `--verify` recounts.

    scripts/unpack_globdb_genomes.py            # unpack everything
    scripts/unpack_globdb_genomes.py --verify   # count what is there
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ARCHIVE = REPO / "refs" / "globdb_r232" / "globdb_r232_genome_fasta.tar.gz"
DEST = REPO / "refs" / "globdb_r232" / "genomes"
LOCK = REPO / "refs" / "expanded" / "locks" / "globdb_genomes.lock.json"
UNPACK_LOCK = DEST / "genomes.unpack.json"
EXPECTED = 346_233


def _count() -> int:
    return sum(1 for _ in DEST.rglob("*.fa.gz"))


def _fetched_and_verified() -> tuple[bool, str]:
    if not LOCK.is_file():
        return False, "no fetch lock; run scripts/fetch_expanded_refs.py --only globdb_genomes --include-deferred"
    lock = json.loads(LOCK.read_text())
    files = lock.get("files") or {}
    rec = next(iter(files.values()), {}) if files else {}
    if not rec.get("verified"):
        return False, "archive not yet verified against the publisher checksum (fetch still running or failed)"
    return True, str(rec.get("sha256") or "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify", action="store_true", help="only count genomes on disk")
    ap.add_argument("--force", action="store_true", help="unpack even if the fetch lock does not say verified")
    ap.add_argument("--rate-limit", type=int, default=0, metavar="MB_PER_S",
                    help="cap the archive read rate (needs `pv`); keeps a 290 GB extraction from monopolising the disk")
    args = ap.parse_args()
    if args.verify:
        n = _count()
        print(f"{n:,} of {EXPECTED:,} GlobDB genomes under {DEST}")
        return 0 if n >= EXPECTED else 1
    ok, detail = _fetched_and_verified()
    if not ok and not args.force:
        print(detail)
        return 2
    if not ARCHIVE.is_file():
        print(f"archive missing: {ARCHIVE}")
        return 2
    DEST.mkdir(parents=True, exist_ok=True)
    before = _count()
    print(f"{before:,} genomes already on disk; unpacking {ARCHIVE.name} ({ARCHIVE.stat().st_size / 1e9:.0f} GB)")
    t0 = time.monotonic()
    # --strip-components drops the top directory; -k never overwrites, so a
    # rerun only fills gaps. Errors for existing files are expected with -k.
    tar_cmd = ["tar", "-xzkf", "-", "-C", str(DEST), "--strip-components", "1"]
    import shutil

    pv = shutil.which("pv") if args.rate_limit else None
    if args.rate_limit and pv is None:
        print("pv not installed; extracting without a rate limit")
    if pv:
        reader = subprocess.Popen([pv, "-q", "-L", f"{args.rate_limit}m", str(ARCHIVE)], stdout=subprocess.PIPE)
        proc = subprocess.run(tar_cmd, stdin=reader.stdout, capture_output=True, text=True, check=False)
        reader.stdout.close()
        reader.wait()
    else:
        with ARCHIVE.open("rb") as fh:
            proc = subprocess.run(tar_cmd, stdin=fh, capture_output=True, text=True, check=False)
    after = _count()
    took = time.monotonic() - t0
    if proc.returncode != 0 and after <= before:
        print(proc.stderr[-2000:])
        return 1
    lock = {
        "unpacked_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "archive": str(ARCHIVE), "archive_sha256": detail if ok else "unverified (--force)",
        "n_genomes": after, "expected": EXPECTED, "complete": after >= EXPECTED,
        "seconds": round(took, 1), "layout": "GlobDB's own chunk directories; <id>.fa.gz",
    }
    UNPACK_LOCK.write_text(json.dumps(lock, indent=2) + "\n")
    print(f"{after:,} genomes on disk after {took / 60:.0f} min "
          f"({'complete' if lock['complete'] else f'{EXPECTED - after:,} short of the published count'})")
    return 0 if lock["complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
