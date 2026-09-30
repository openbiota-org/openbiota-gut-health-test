#!/usr/bin/env python3
"""Fetch the expanded references once: resumable, verified, locked.

    scripts/fetch_expanded_refs.py [--only SOURCE_ID ...] [--jobs N] [--status]

Every file in `openbiota.refsources.SOURCES` is downloaded with resume
(`curl -C -`), verified against the publisher's checksum where one is
published, hashed locally (SHA-256), and recorded in
`refs/expanded/locks/<source>.lock.json` with URL, size, retrieval time,
release, taxonomy release and licence. A file whose lock already verifies
is not downloaded again - there are no unbounded repeated whole-reference
downloads (spec §7).

Progress goes to `refs/expanded/fetch.progress`, one line per file, so a
watcher can see bytes done against bytes expected without touching the
downloads. `--status` prints that view and exits.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from openbiota.refsources import LOCK_DIR, PROGRESS_FILE, SOURCES, RefFile  # noqa: E402

CONNECTIONS = 8

_lock = threading.Lock()
_progress: dict[str, tuple[int, int, str]] = {}   # dest -> (done, total, state)
_md5_lists: dict[str, dict[str, str]] = {}


def _write_progress() -> None:
    PROGRESS_FILE.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{dt.datetime.now().strftime('%H:%M:%S')}"]
    done_total = sum(d for d, _, _ in _progress.values())
    all_total = sum(t for _, t, _ in _progress.values())
    lines.append(f"TOTAL {done_total/1e9:.2f} / {all_total/1e9:.2f} GB")
    for dest, (d, t, state) in sorted(_progress.items()):
        pct = (100.0 * d / t) if t else 0.0
        lines.append(f"{state:9} {pct:5.1f}%  {d/1e9:7.2f}/{t/1e9:7.2f} GB  {dest}")
    PROGRESS_FILE.write_text("\n".join(lines) + "\n")


def _set(dest: str, done: int, total: int, state: str) -> None:
    with _lock:
        _progress[dest] = (done, total, state)
        _write_progress()


def _remote_size(url: str) -> int:
    try:
        out = subprocess.run(["curl", "-sIL", "--max-time", "60", url], capture_output=True, text=True, check=False).stdout
        sizes = re.findall(r"(?i)content-length:\s*(\d+)", out)
        return int(sizes[-1]) if sizes else 0
    except Exception:  # noqa: BLE001
        return 0


def _digests(path: Path) -> tuple[str, str]:
    """(sha256, md5) in one pass: a 290 GB archive is read once, not twice."""
    sha = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - publisher checksums are md5
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest()


def _sha256(path: Path) -> str:
    return _digests(path)[0]


def _md5(path: Path) -> str:
    return _digests(path)[1]


def _aria2c() -> str | None:
    """aria2c when installed: parallel range requests finish a large file in
    minutes where one curl stream takes hours. Resumes a partial file."""
    import shutil

    return shutil.which("aria2c")


def _publisher_checksum(f: RefFile) -> tuple[str, str] | None:
    """(algorithm, hex) from the declaration, resolving md5 lists by basename."""
    spec = f.publisher_checksum
    if not spec:
        return None
    kind, _, value = spec.partition(":")
    if kind in ("md5", "sha256"):
        return kind, value.lower()
    if kind == "md5-list":
        with _lock:
            table = _md5_lists.get(value)
        if table is None:
            table = {}
            try:
                text = urllib.request.urlopen(value, timeout=60).read().decode("utf-8", "replace")
                for line in text.splitlines():
                    parts = line.split()
                    if len(parts) >= 2:
                        table[Path(parts[-1]).name] = parts[0].lower()
            except Exception:  # noqa: BLE001
                pass
            with _lock:
                _md5_lists[value] = table
        want = table.get(f.path.name)
        return ("md5", want) if want else None
    return None


def _lock_path(f: RefFile) -> Path:
    return LOCK_DIR / f"{f.source_id}.lock.json"


def _read_lock(f: RefFile) -> dict:
    p = _lock_path(f)
    return json.loads(p.read_text()) if p.is_file() else {"source_id": f.source_id, "files": {}}


def _write_lock(f: RefFile, record: dict) -> None:
    with _lock:
        lock = _read_lock(f)
        lock.update({
            "source_id": f.source_id, "release_id": f.release_id,
            "taxonomy_release": f.taxonomy_release, "license": f.license,
            "count_unit": f.count_unit, "source_count": f.source_count,
        })
        lock["files"][f.dest] = record
        LOCK_DIR.mkdir(parents=True, exist_ok=True)
        _lock_path(f).write_text(json.dumps(lock, indent=2) + "\n")


def _verified(f: RefFile) -> bool:
    """True when the lock records this file and the file on disk matches it."""
    rec = _read_lock(f)["files"].get(f.dest)
    if not rec or not f.path.is_file():
        return False
    return f.path.stat().st_size == rec.get("size") and rec.get("verified") is True


def fetch(f: RefFile) -> tuple[RefFile, str]:
    if _verified(f):
        _set(f.dest, f.path.stat().st_size, f.path.stat().st_size, "ok")
        return f, "already verified"
    f.path.parent.mkdir(parents=True, exist_ok=True)
    total = _remote_size(f.url)
    _set(f.dest, f.path.stat().st_size if f.path.exists() else 0, total, "fetching")

    # Resume with aria2c (parallel connections) when available, else curl.
    # Both continue a partial file; the poll below reports the size.
    aria = _aria2c()
    if aria and (total or 0) >= 1 << 30:
        cmd = [aria, "-c", "-x", str(CONNECTIONS), "-s", str(CONNECTIONS), "-k", "64M", "--file-allocation=none",
               "--max-tries=8", "--retry-wait=5", "--summary-interval=0", "--console-log-level=warn",
               "-d", str(f.path.parent), "-o", f.path.name, f.url]
    else:
        cmd = ["curl", "-sSL", "--fail", "--retry", "8", "--retry-delay", "5", "--retry-all-errors",
               "-C", "-", "-o", str(f.path), f.url]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    while proc.poll() is None:
        time.sleep(3)
        if f.path.exists():
            _set(f.dest, f.path.stat().st_size, total, "fetching")
    err = (proc.stderr.read() if proc.stderr else "").strip()
    if proc.returncode not in (0,):
        # curl 33 = server does not support resume on a complete file; 22 with 416 similar
        if f.path.exists() and total and f.path.stat().st_size == total:
            pass
        else:
            _set(f.dest, f.path.stat().st_size if f.path.exists() else 0, total, "FAILED")
            return f, f"{Path(cmd[0]).name} exit {proc.returncode}: {err[-200:]}"
    size = f.path.stat().st_size
    if total and size != total:
        _set(f.dest, size, total, "FAILED")
        return f, f"size {size} != expected {total}"

    _set(f.dest, size, total or size, "verifying")
    sha256_hex, md5_hex = _digests(f.path)
    record = {
        "url": f.url, "size": size, "retrieved_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "sha256": sha256_hex, "verified": False, "publisher_checksum": None, "note": f.note,
        "downloader": Path(cmd[0]).name,
    }
    pub = _publisher_checksum(f)
    if pub:
        algo, want = pub
        got = md5_hex if algo == "md5" else record["sha256"]
        record["publisher_checksum"] = f"{algo}:{want}"
        if got != want:
            record["verified"] = False
            record["mismatch"] = f"{algo} {got} != publisher {want}"
            _write_lock(f, record)
            _set(f.dest, size, total or size, "MISMATCH")
            return f, f"checksum mismatch ({algo})"
        record["verified"] = True
        record["verification"] = f"publisher {algo}"
    else:
        record["verified"] = True
        record["verification"] = "size matched Content-Length; no publisher checksum published"
    _write_lock(f, record)
    _set(f.dest, size, total or size, "ok")
    return f, "fetched and verified" if pub else "fetched (size-verified)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None, help="source ids to fetch")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--connections", type=int, default=8, help="parallel connections per large file (aria2c)")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--include-deferred", action="store_true")
    args = ap.parse_args()
    global CONNECTIONS
    CONNECTIONS = max(1, args.connections)

    if args.status:
        print(PROGRESS_FILE.read_text() if PROGRESS_FILE.is_file() else "no fetch has run")
        return 0

    wanted = [f for f in SOURCES if (args.only is None or f.source_id in args.only)
              and (args.include_deferred or not f.deferred)]
    print(f"{len(wanted)} files across {len({f.source_id for f in wanted})} sources; {args.jobs} parallel")
    for f in wanted:
        _progress[f.dest] = (f.path.stat().st_size if f.path.exists() else 0, 0, "queued")
    _write_progress()

    failures = 0
    # Big files first so the long tail starts early.
    wanted.sort(key=lambda f: -(f.path.stat().st_size if f.path.exists() else 0))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(fetch, f): f for f in wanted}
        for fut in as_completed(futures):
            f, msg = fut.result()
            ok = not (msg.startswith("curl") or "mismatch" in msg or "!=" in msg)
            failures += 0 if ok else 1
            print(f"  {'ok ' if ok else 'ERR'} {f.dest:70} {msg}")
    print(f"\n{'all verified' if not failures else f'{failures} failure(s)'}; locks in {LOCK_DIR}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
