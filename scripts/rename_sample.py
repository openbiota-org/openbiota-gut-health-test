#!/usr/bin/env python3
"""Rename a sample - its FASTQ files already renamed - and carry every cache with it.

    scripts/rename_sample.py OLD_ID NEW_ID [--old-root /path/of/the/repo/before/it/moved] [--dry-run]

Every stage of a run caches its work under results/<sample>/ and keys that
cache on the input reads' *name* (and sometimes absolute path), size and
mtime - a cheap identity that tells two files apart without hashing gigabytes.
Renaming the reads, or moving the repository, therefore invalidates hours of
cached work that is still perfectly valid: the bytes did not change. This tool
moves the results directory, renames every file that carries the old id, and
re-keys each cache by computing its key twice with the pipeline's own
functions - once as the file was named, once as it is now - so the next run
finds everything "cached".

Caches handled (each by the function the stage itself uses):
  input validation      refs/cache/input_profile_<key>.json         qc._cache_key (path)
  fastp                 preprocess/fastp.<mode>.stamp.json          preprocess._fingerprint (resolved path)
  DIAMOND               alignments/*.done.json                      search._manifest_payload (path string)
  MetaPhlAn 3 / Jun23   taxonomy/metaphlan*.json "inputs"           engines.metaphlan._fingerprint (name)
  host filter           taxonomy/host_filter.json "inputs"          same
  Jan26, mOTUs, SingleM, sylph, sylph-GlobDB, Kraken x2
                        <lane>/<lane>.<key>.json                    each lane's cache_file()
  pathogen candidates   pathogens/candidates.json                   relative paths rewritten
  competitive confirm.  confirmation/<sample>.confirmation.*.json   renamed; confirm.py adopts a match on content
The LP caches (refs/micom/cache) are keyed on the community, not the sample,
and need nothing. results.json and the PDF are rewritten by the next run.

Sizes and mtimes must be unchanged (mv and directory renames preserve them);
the tool refuses to proceed if the new reads are missing.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota import preprocess, qc  # noqa: E402
from openbiota.engines import (  # noqa: E402
    kraken,
    metaphlan,
    metaphlan_jan26,
    motus,
    singlem,
    sylph,
    sylph_globdb,
)


class AsNamed:
    """A path as it *was*: the old name and location, the real file's stat.

    The fingerprint functions read `.name`, `.stat()`, `.resolve()` and
    `str()`; nothing else. This stands in for the old file without needing
    it to exist.
    """

    def __init__(self, old: Path, real: Path) -> None:
        self._old, self._real = old, real

    @property
    def name(self) -> str:
        return self._old.name

    def stat(self) -> os.stat_result:
        return self._real.stat()

    def resolve(self) -> Path:
        return self._old

    def __str__(self) -> str:
        return str(self._old)

    def __fspath__(self) -> str:
        return str(self._old)


def _old_root_from_manifests(results_dir: Path) -> Path | None:
    """The repository path recorded by the last run, from any DIAMOND manifest."""
    for m in sorted((results_dir / "alignments").glob("*.done.json")):
        try:
            q = json.loads(m.read_text()).get("query", "")
        except (OSError, ValueError):
            continue
        if q:
            p = Path(q)
            if p.parent.name == "fastq":
                return p.parent.parent
    return None


def _rename_tree(results_dir: Path, old_id: str, new_id: str, *, dry: bool) -> int:
    """Rename every file and directory under results_dir whose name carries old_id (deepest first)."""
    n = 0
    for path in sorted(results_dir.rglob(f"*{old_id}*"), key=lambda p: -len(p.parts)):
        target = path.with_name(path.name.replace(old_id, new_id))
        print(f"  mv {path.relative_to(REPO)} -> {target.name}")
        if not dry:
            path.rename(target)
        n += 1
    return n


def _rewrite_text(path: Path, pairs: list[tuple[str, str]], *, dry: bool) -> bool:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False
    new = text
    for old, repl in pairs:
        new = new.replace(old, repl)
    if new != text:
        print(f"  rewrite {path.relative_to(REPO)}")
        if not dry:
            path.write_text(new, encoding="utf-8")
        return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old_id")
    ap.add_argument("new_id")
    ap.add_argument("--old-root", type=Path, default=None,
                    help="repository path at the time of the last run, if it moved (default: read from a DIAMOND manifest)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    old_id, new_id, dry = args.old_id, args.new_id, args.dry_run
    moved_only = old_id == new_id  # the reads kept their names; only the repository moved

    new_r1, new_r2 = REPO / "fastq" / f"{new_id}_1.fastq.gz", REPO / "fastq" / f"{new_id}_2.fastq.gz"
    if not new_r1.is_file():
        sys.exit(f"{new_r1} does not exist - rename the FASTQ files first")
    r2_exists = new_r2.is_file()
    results_dir = REPO / "results" / old_id
    if not results_dir.is_dir():
        results_dir = REPO / "results" / new_id
        if not results_dir.is_dir():
            sys.exit(f"no results directory for {old_id} or {new_id}")
    old_root = args.old_root or _old_root_from_manifests(results_dir) or REPO
    print(f"old id {old_id} -> new id {new_id}; repository then {old_root}, now {REPO}")

    old_r1 = AsNamed(old_root / "fastq" / f"{old_id}_1.fastq.gz", new_r1)
    old_r2 = AsNamed(old_root / "fastq" / f"{old_id}_2.fastq.gz", new_r2) if r2_exists else None
    cur_r2 = new_r2 if r2_exists else None
    pairs = [(str(old_root), str(REPO)), (old_id, new_id)]
    if str(old_root) == str(REPO):
        pairs = [(old_id, new_id)]
    if moved_only and str(old_root) == str(REPO):
        sys.exit("nothing to do: same id and the repository has not moved")

    # 1. input validation cache in refs/cache (keyed on the full path)
    print("\ninput validation")
    cache_dir = REPO / "refs" / "cache"
    for sample_reads, count_all in ((qc.DEFAULT_SAMPLE_READS, True), (qc.DEFAULT_SAMPLE_READS, False)):
        mates_old = [("R1", old_r1)] + ([("R2", old_r2)] if old_r2 else [])
        mates_new = [("R1", new_r1)] + ([("R2", cur_r2)] if cur_r2 else [])
        k_old = qc._cache_key(mates_old, sample_reads, count_all)
        k_new = qc._cache_key(mates_new, sample_reads, count_all)
        src, dst = cache_dir / f"input_profile_{k_old}.json", cache_dir / f"input_profile_{k_new}.json"
        if src.is_file():
            print(f"  mv refs/cache/{src.name} -> {dst.name}")
            if not dry:
                src.rename(dst)

    # 2. the results tree: directory, then every file named with the old id
    print("\nresults tree")
    if moved_only:
        print("  same id - nothing to rename")
    elif results_dir.name == old_id:
        print(f"  mv results/{old_id} -> results/{new_id}")
        if not dry:
            results_dir = results_dir.rename(REPO / "results" / new_id)
        else:
            results_dir = REPO / "results" / new_id
    n = 0 if moved_only else _rename_tree(results_dir if not dry else REPO / "results" / old_id, old_id, new_id, dry=dry)
    print(f"  {n} path(s) renamed")
    for extra in () if moved_only else sorted((REPO / "results").glob(f"*{old_id}*")):
        target = extra.with_name(extra.name.replace(old_id, new_id))
        print(f"  mv results/{extra.name} -> {target.name}")
        if not dry:
            extra.rename(target)
    for extra in () if moved_only else sorted((REPO / "results" / "_sylph").rglob(f"*{old_id}*")):
        target = extra.with_name(extra.name.replace(old_id, new_id))
        print(f"  mv {extra.relative_to(REPO)} -> {target.name}")
        if not dry:
            extra.rename(target)
    if dry:
        results_dir = REPO / "results" / old_id

    # 3. fastp stamp (resolved path)
    print("\nfastp")
    for stamp in sorted((results_dir / "preprocess").glob("fastp.*.stamp.json")):
        meta = json.loads(stamp.read_text())
        mode = stamp.name.split(".")[1]
        cmd = meta.get("command") or []
        # the arguments that formed the key, as the stage builds them
        detect = "--detect_adapter_for_pe" in cmd
        min_len = next((cmd[i + 1] for i, a in enumerate(cmd) if a in ("-l", "--length_required")), None)
        # the installed fastp's version is part of the key; the stamp does not record it
        try:
            version = preprocess.locate_fastp().version
        except Exception:  # noqa: BLE001 - without fastp there is nothing to re-key
            version = ""
        extra_old = f"{mode}|{detect}|{min_len}|{version}"
        inputs_old = [old_r1] + ([old_r2] if old_r2 else [])
        inputs_new = [new_r1] + ([cur_r2] if cur_r2 else [])
        if preprocess._fingerprint(inputs_old, extra_old) != meta.get("key"):
            # the stage's extras are not fully recoverable from the stamp; say so and let fastp re-run (about a minute)
            print(f"  {stamp.name}: key does not reproduce from the stamp; fastp will re-run")
            continue
        meta["key"] = preprocess._fingerprint(inputs_new, extra_old)
        meta["command"] = [str(c).replace(str(old_root), str(REPO)).replace(old_id, new_id) for c in cmd]
        print(f"  {stamp.name}: key re-computed")
        if not dry:
            stamp.write_text(json.dumps(meta, indent=1))

    # 4. DIAMOND manifests (path string)
    print("\nDIAMOND")
    for m in sorted((results_dir / "alignments").glob("*.done.json")):
        _rewrite_text(m, pairs, dry=dry)

    # 5. MetaPhlAn 3 / Jun23 / host filter stamps (name)
    print("\nMetaPhlAn 3, MetaPhlAn 4 (Jun23), host filter")
    fp_old = metaphlan._fingerprint(old_r1, old_r2)
    fp_new = metaphlan._fingerprint(new_r1, cur_r2)
    for stamp in sorted((results_dir / "taxonomy").glob("*.json")):
        try:
            meta = json.loads(stamp.read_text())
        except ValueError:
            continue
        if isinstance(meta, dict) and meta.get("inputs") == fp_old:
            meta["inputs"] = fp_new
            if "command" in meta:
                meta["command"] = [str(c).replace(str(old_root), str(REPO)).replace(old_id, new_id) for c in meta["command"]]
            print(f"  {stamp.name}: inputs {fp_old} -> {fp_new}")
            if not dry:
                stamp.write_text(json.dumps(meta, indent=1))

    # 6. the expansion lanes: <lane>.<key>.json keyed through each lane's cache_file()
    print("\nexpansion lanes")
    lanes = [
        (results_dir / "taxonomy", lambda r1, r2: metaphlan_jan26.cache_file(results_dir / "taxonomy", r1, r2)),
        (results_dir / "motus", lambda r1, r2: motus.cache_file(results_dir / "motus", r1, r2)),
        (results_dir / "singlem", lambda r1, r2: singlem.cache_file(results_dir / "singlem", r1, r2)),
        (results_dir / "globdb", lambda r1, r2: sylph_globdb.cache_file(results_dir / "globdb", r1, r2)),
        (results_dir / "genome", lambda r1, r2: sylph.cache_file(results_dir / "genome", r1, r2)),
        (results_dir / "kraken", lambda r1, r2: kraken.cache_file(results_dir / "kraken", r1, r2, refs_dir=REPO / "refs", panel="uhgg")),
        (results_dir / "kraken_rescue", lambda r1, r2: kraken.cache_file(results_dir / "kraken_rescue", r1, r2, refs_dir=REPO / "refs", panel="rescue")),
    ]
    for work_dir, keyed in lanes:
        if not work_dir.is_dir():
            continue
        try:
            src, dst = keyed(old_r1, old_r2), keyed(new_r1, cur_r2)
        except Exception as exc:  # noqa: BLE001 - a lane whose tool is absent cannot be re-keyed here
            print(f"  {work_dir.name}: cannot compute key ({type(exc).__name__}: {exc}); the lane will re-run")
            continue
        if src.is_file() and src != dst:
            print(f"  mv {src.relative_to(REPO)} -> {dst.name}")
            if not dry:
                src.rename(dst)
        elif dst.is_file():
            print(f"  {dst.name}: already keyed for the new name")
        else:
            print(f"  {work_dir.name}: no cache file matched the old key; the lane will re-run")
        # paths inside the lane records
        for rec in sorted(work_dir.glob("*.json")):
            _rewrite_text(rec, pairs, dry=dry)

    # 7. everything else that records a path or the id as text
    print("\nother records")
    for rel in ("pathogens/candidates.json", "strain/consensus_markers.json", "extension_manifest.json",
                "depth_check.json", "capability_coverage.json"):
        path = results_dir / rel
        if path.is_file():
            _rewrite_text(path, pairs, dry=dry)
    for sub in ("confirmation", "strain", "assembly", "mycobiome", "genome", "globdb"):
        d = results_dir / sub
        if d.is_dir():
            for path in sorted(d.rglob("*.json")):
                _rewrite_text(path, pairs, dry=dry)

    # 8. the alias file
    alias_file = REPO / "results" / ".sample_aliases.json"
    if alias_file.is_file():
        aliases = json.loads(alias_file.read_text())
        changed = {k: (new_id if v == old_id else v) for k, v in aliases.items()}
        if changed != aliases:
            print(f"\nresults/.sample_aliases.json: {old_id} -> {new_id}")
            if not dry:
                alias_file.write_text(json.dumps(changed, indent=2) + "\n")

    # sanity: nothing left under the old name
    left = [] if (dry or moved_only) else [str(p.relative_to(REPO)) for p in (REPO / "results").rglob(f"*{old_id}*")]
    if left:
        print(f"\nstill named with {old_id}: {left[:5]}")
        return 1
    print("\ndone" if not dry else "\ndry run - nothing changed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
