#!/usr/bin/env python3
"""Pack and unpack the reference data that travels with the repository.

    scripts/derived_refs.py pack [BUNDLE]      # rebuild data/<bundle>.tar.gz* from refs/
    scripts/derived_refs.py restore [BUNDLE]   # unpack into refs/
    scripts/derived_refs.py verify [BUNDLE]    # check the archives against their manifests
    scripts/derived_refs.py list               # what is in them, and why

``refs/`` is about 1.3 TB, and almost all of it is somebody else's: GlobDB,
GTDB, the MetaPhlAn indexes, Kraken2, the AGORA2 models. Those have homes of
their own and the installer fetches them from those homes, so the project
never has to host them and nobody has to pay to serve them.

Two bundles travel with the repository instead:

``derived``
    What this project computed itself: the reference ranges the whole report
    is scored against, the scoring cohort matrix, the trained age model, the
    GTDB crosswalks, the reconciliation locks, the resolution lanes'
    reference sequences and the BacDive cache. About 50 MB compressed.
    Rebuilding them is not a download - it is hours of compute and, for
    BacDive, an API account.

``uniprot``
    The UniProt REST responses behind the protein search database: one
    gzipped table per reference entry, 244 of them, about 200 MB. UniProt
    answers these queries at roughly ninety seconds each, which makes the
    database build a six-hour step on a fresh machine; with the responses
    present it is a two-minute local DIAMOND build. UniProt's data is
    CC BY 4.0; every file's sidecar records the query that produced it.

Each bundle is a plain ``tar.gz`` - every machine that can clone the
repository can read it - cut into parts under GitHub's 50 MB per-file
warning, and ``restore`` joins the parts back. A manifest beside the parts
holds their sizes and SHA-256s, and ``restore`` refuses an archive that does
not match it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tarfile
import tempfile
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ARCHIVE_DIR = REPO / "data"

#: Parts stay under GitHub's 50 MB per-file warning (its hard limit is 100 MB).
PART_BYTES: int = 45 * 1024 * 1024


@dataclass(frozen=True)
class Item:
    """One thing in a bundle, and the reason it is there."""

    path: str
    why: str
    #: False when the installer can do without it (and so can `pack`).
    required: bool = True


#: Everything this project computed rather than downloaded. The paths are
#: relative to the repository root and may be files or directories.
DERIVED: tuple[Item, ...] = (
    Item("refs/reference_ranges.json",
         "the percentile ranges every pathway in the report is scored against, "
         "computed from the reference cohort"),
    Item("refs/reference_cohort_samples.json",
         "which public samples the pathway cohort was built from"),
    Item("refs/taxonomic_cohort.json",
         "the scoring cohort matrix: species abundances across the reference "
         "adults, from curatedMetagenomicData"),
    Item("refs/taxonomic_cohort_manifest.json", "what that matrix was built from"),
    Item("refs/cmd",
         "the curatedMetagenomicData extract behind the scoring cohort; "
         "rebuilding it needs R and ExperimentHub"),
    Item("refs/age", "the trained microbiome-age models and their model cards"),
    Item("refs/crosswalk",
         "species-name crosswalks between GTDB releases and the MetaPhlAn "
         "catalogues, derived from the bridges"),
    Item("refs/expanded",
         "the reference locks and the reconciliation summary: what every "
         "expansion lane was pinned to"),
    Item("refs/expansion_cohort",
         "expansion lanes profiled locally from public reads - days of compute"),
    Item("refs/ctnpc", "the curated toxin and non-pathogen control references"),
    Item("refs/panels",
         "the toxin-locus, probiotic-strain and control reference sequences "
         "the resolution lanes align against; small, lock-verified NCBI "
         "fetches that are flaky to repeat"),
    Item("refs/bacdive",
         "the BacDive API cache, so organism enrichment works without an "
         "account of your own", required=False),
    Item("refs/strain_refs.lock.json", "what the strain lane was pinned to"),
    Item("refs/pathogens/resolutions.cache.json",
         "how every pathogen catalogue name resolved at NCBI, so the full "
         "install's pathogen build does not repeat 377 taxonomy lookups",
         required=False),
    Item("refs/MANIFEST.json", "the reference inventory"),
    Item("refs/MANIFEST.md", "the reference inventory, readable"),
)


def uniprot_items() -> tuple[Item, ...]:
    """The cached UniProt responses the *current* panels need - no more.

    Computed from the panel files through the pipeline's own cache-key
    function, so the bundle always matches what ``openbiota build-db`` will
    look for; the historical responses of panels since changed stay out.
    Needs the package importable, which ``pack`` has and ``restore`` does
    not need.
    """
    sys.path.insert(0, str(REPO))
    from openbiota.cli import _entries_for_db  # noqa: PLC0415
    from openbiota.panels import load_panel_set  # noqa: PLC0415
    from openbiota.references import _entry_cache_key  # noqa: PLC0415

    panels = load_panel_set(REPO / "panels")
    out: list[Item] = []
    for entry in _entries_for_db(panels, None):
        if entry.source != "uniprot":
            continue
        limit = max(entry.max_sequences * 2, 1000)   # as fetch_entry_records does
        stem = f"refs/cache/{entry.id.lower()}_{_entry_cache_key(entry, limit)}"
        why = f"UniProt response for {entry.key} ({entry.role})"
        out.append(Item(f"{stem}.tsv.gz", why))
        out.append(Item(f"{stem}.tsv.meta.json", f"query and counts behind {entry.key}"))
    return tuple(out)


@dataclass(frozen=True)
class Bundle:
    name: str
    items: Callable[[], Sequence[Item]]
    summary: str

    @property
    def archive(self) -> Path:
        return ARCHIVE_DIR / f"{self.name}.tar.gz"

    @property
    def manifest(self) -> Path:
        return ARCHIVE_DIR / f"{self.name}.manifest.json"

    def parts(self) -> list[Path]:
        """The archive as one file, or the parts it was split into."""
        if self.archive.is_file():
            return [self.archive]
        return sorted(ARCHIVE_DIR.glob(self.archive.name + ".part*"))


BUNDLES: dict[str, Bundle] = {
    "derived-refs": Bundle(
        "derived-refs", lambda: DERIVED,
        "reference data this project computed: ranges, cohorts, age model, "
        "crosswalks, locks, resolution references, BacDive cache"),
    "uniprot-cache": Bundle(
        "uniprot-cache", uniprot_items,
        "the UniProt responses behind the protein search database "
        "(CC BY 4.0), so build-db is a local two-minute step"),
}
# What `restore` does when no bundle is named: the derived data first, the
# UniProt cache second, so a machine without the second still has the first.
DEFAULT_ORDER = ("derived-refs", "uniprot-cache")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _reproducible(info: tarfile.TarInfo) -> tarfile.TarInfo:
    """Drop the things that make two packs of the same data differ."""
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.mtime = 0
    return info


# --------------------------------------------------------------------------- #
# pack
# --------------------------------------------------------------------------- #


def pack(bundle: Bundle, *, quiet: bool = False) -> int:
    items = list(bundle.items())
    missing = [i for i in items if i.required and not (REPO / i.path).exists()]
    if missing:
        print(f"cannot pack {bundle.name}: these are required and not present:", file=sys.stderr)
        for item in missing:
            print(f"  {item.path}", file=sys.stderr)
        return 1

    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    for stale in (*bundle.parts(), bundle.manifest):
        stale.unlink(missing_ok=True)

    included = [i for i in items if (REPO / i.path).exists()]
    started = time.monotonic()
    n_files = 0
    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp) / bundle.archive.name
        # Already-compressed members (the UniProt tables) do not shrink again;
        # level 6 keeps the pack quick without costing the others anything.
        with tarfile.open(staged, "w:gz", compresslevel=6) as tar:
            for item in included:
                source = REPO / item.path
                if not quiet and len(included) <= 40:
                    print(f"  packing {item.path}")
                tar.add(source, arcname=item.path, filter=_reproducible)
                n_files += 1 if source.is_file() else sum(
                    1 for _ in source.rglob("*") if _.is_file())
        if staged.stat().st_size <= PART_BYTES:
            shutil.move(str(staged), bundle.archive)
            written = [bundle.archive]
        else:
            written = _split(bundle, staged)

    manifest = {
        "schema": "openbiota.bundle/1",
        "bundle": bundle.name,
        "summary": bundle.summary,
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_bytes": sum(p.stat().st_size for p in written),
        "n_files": n_files,
        "parts": [
            {"name": p.name, "bytes": p.stat().st_size, "sha256": _sha256(p)}
            for p in written
        ],
        "contents": [
            {"path": i.path, "why": i.why, "required": i.required} for i in included
        ],
    }
    bundle.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    if not quiet:
        print(f"{bundle.name}: {len(written)} part(s), {_human(manifest['total_bytes'])}, "
              f"{n_files:,} files in {time.monotonic() - started:.0f}s")
        for part in written:
            print(f"  {part.relative_to(REPO)}  {_human(part.stat().st_size)}")
    return 0


def _split(bundle: Bundle, source: Path) -> list[Path]:
    out: list[Path] = []
    with source.open("rb") as handle:
        index = 0
        while True:
            chunk = handle.read(PART_BYTES)
            if not chunk:
                break
            part = ARCHIVE_DIR / f"{bundle.archive.name}.part{index:02d}"
            part.write_bytes(chunk)
            out.append(part)
            index += 1
    return out


# --------------------------------------------------------------------------- #
# restore
# --------------------------------------------------------------------------- #


def restore(bundle: Bundle, *, force: bool = False, quiet: bool = False) -> int:
    parts = bundle.parts()
    if not parts:
        print(f"no archive for {bundle.name} under {ARCHIVE_DIR.relative_to(REPO)}/", file=sys.stderr)
        return 1
    if not verify(bundle, quiet=True):
        print(f"{bundle.name} does not match its manifest; refusing to restore", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory() as tmp:
        joined = parts[0]
        if len(parts) > 1:
            joined = Path(tmp) / bundle.archive.name
            with joined.open("wb") as out:
                for part in parts:
                    out.write(part.read_bytes())
        with tarfile.open(joined, "r:gz") as tar:
            members = tar.getmembers()
            skipped = 0
            root = REPO.resolve()
            for member in members:
                target = (REPO / member.name).resolve()
                # A path that escapes the repository is a broken or hostile
                # archive, not something to extract and find out.
                if root != target and root not in target.parents:
                    print(f"refusing path outside the repository: {member.name}", file=sys.stderr)
                    return 1
                if target.exists() and not force:
                    skipped += 1
                    continue
                try:
                    tar.extract(member, REPO, filter="data")
                except TypeError:
                    # Python before 3.11.4 has no extraction filter; the path
                    # check above is the protection it would have given.
                    tar.extract(member, REPO)
    if not quiet:
        kept = f", {skipped} already present and left alone" if skipped else ""
        print(f"{bundle.name}: restored {len(members) - skipped} entries into refs/{kept}")
    return 0


# --------------------------------------------------------------------------- #
# check / verify / list
# --------------------------------------------------------------------------- #


def check(bundle: Bundle, *, quiet: bool = False) -> bool:
    """Is everything the bundle carries present under refs/?

    The installer's question before it unpacks: a yes means the step can be
    skipped. Answered from the manifest, not the archive, so it costs
    nothing; a path listed as required that is absent is a no.
    """
    if not bundle.manifest.is_file():
        if not quiet:
            print(f"{bundle.name}: no manifest, nothing to check against", file=sys.stderr)
        return False
    manifest = json.loads(bundle.manifest.read_text())
    missing = [
        entry["path"] for entry in manifest.get("contents", [])
        if entry.get("required", True) and not (REPO / entry["path"]).exists()
    ]
    if not quiet:
        if missing:
            print(f"{bundle.name}: {len(missing)} of {len(manifest.get('contents', []))} "
                  f"items missing from refs/ (first: {missing[0]})")
        else:
            print(f"{bundle.name}: present")
    return not missing


def verify(bundle: Bundle, *, quiet: bool = False) -> bool:
    if not bundle.manifest.is_file():
        if not quiet:
            print(f"no manifest for {bundle.name}", file=sys.stderr)
        return False
    try:
        manifest = json.loads(bundle.manifest.read_text())
    except (OSError, ValueError) as exc:
        if not quiet:
            print(f"{bundle.name}: manifest unreadable: {exc}", file=sys.stderr)
        return False
    ok = True
    for part in manifest.get("parts", []):
        path = ARCHIVE_DIR / part["name"]
        if not path.is_file():
            ok = False
            if not quiet:
                print(f"  MISSING  {part['name']}", file=sys.stderr)
            continue
        if path.stat().st_size != part["bytes"]:
            ok = False
            if not quiet:
                print(f"  WRONG SIZE  {part['name']}", file=sys.stderr)
            continue
        if _sha256(path) != part["sha256"]:
            ok = False
            if not quiet:
                print(f"  CHECKSUM MISMATCH  {part['name']}", file=sys.stderr)
            continue
        if not quiet:
            print(f"  ok  {part['name']}  {_human(path.stat().st_size)}")
    if not quiet and ok:
        print(f"{bundle.name} matches its manifest ({manifest.get('n_files', 0):,} files)")
    return ok


def show() -> int:
    print("What travels with the repository, and why:\n")
    for name in DEFAULT_ORDER:
        bundle = BUNDLES[name]
        parts = bundle.parts()
        size = _human(sum(p.stat().st_size for p in parts)) if parts else "not packed"
        print(f"{bundle.name}  ({len(parts)} part(s), {size})\n  {bundle.summary}\n")
        if name == "derived-refs":
            width = max(len(i.path) for i in DERIVED)
            for item in DERIVED:
                mark = " " if (REPO / item.path).exists() else "?"
                need = "" if item.required else "  (optional)"
                print(f"  {mark} {item.path:<{width}}  {item.why}{need}")
        elif bundle.manifest.is_file():
            n = len(json.loads(bundle.manifest.read_text()).get("contents", []))
            print(f"    {n} files; `verify` checks them, the manifest lists each with its query")
        print()
    print("Everything else under refs/ is downloaded from its own home by the "
          "installer, and is not in here.")
    return 0


def _selected(name: str | None) -> list[Bundle]:
    if name is None or name == "all":
        return [BUNDLES[n] for n in DEFAULT_ORDER]
    if name not in BUNDLES:
        raise SystemExit(f"unknown bundle {name!r}; one of: {', '.join(DEFAULT_ORDER)}, all")
    return [BUNDLES[name]]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=("pack", "restore", "verify", "check", "list"))
    ap.add_argument("bundle", nargs="?", help="derived-refs, uniprot-cache or all (default)")
    ap.add_argument("--force", action="store_true",
                    help="restore over files that are already there")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    if args.action == "list":
        return show()
    worst = 0
    for bundle in _selected(args.bundle):
        if args.action == "pack":
            worst |= pack(bundle, quiet=args.quiet)
        elif args.action == "restore":
            worst |= restore(bundle, force=args.force, quiet=args.quiet)
        elif args.action == "check":
            worst |= 0 if check(bundle, quiet=args.quiet) else 1
        else:
            worst |= 0 if verify(bundle, quiet=args.quiet) else 1
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
