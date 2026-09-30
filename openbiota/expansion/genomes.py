"""Fetch one genome by identifier, from whichever publisher holds it.

Accession-versioned where the source has accessions (NCBI), publisher-local
identifiers where it does not (UHGG species representatives, HRGM2, HROM,
GlobDB units). Every fetch is cached under refs/genomes/<id>.fna.gz with a
sidecar recording URL, size and SHA-256; nothing is fetched twice.

GlobDB publishes no per-genome endpoint. Its units become fetchable once
`refs/globdb_r232/genomes/` has been populated from the one-time archive
(`scripts/fetch_expanded_refs.py --only globdb_genomes --include-deferred`,
then `scripts/unpack_globdb_genomes.py`). Until then a GlobDB-only unit
reports `pending_archive`, and a candidate that needs it stays provisional
with that reason - never confirmed by default.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Final

CACHE: Final = Path("refs/genomes")
GLOBDB_GENOMES: Final = Path("refs/globdb_r232/genomes")
HRGM2_REPS: Final = Path("refs/supplements/hrgm2/HRGMv2_Rep_Genome")
UHGG_SPECIES: Final = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/species_catalogue/{prefix}/{gid}/genome/{gid}.fna"
HROM_GENOME: Final = "https://www.decodebiome.org/HROM/data/genome_catalog/HROM_nonredundant_genomes/{gid}/{gid}_1.fna"
NCBI_FTP: Final = "https://ftp.ncbi.nlm.nih.gov/genomes/all/{kind}/{a}/{b}/{c}/"

_ACC = re.compile(r"^(GC[AF])_(\d{9})(?:\.(\d+))?$")


class Pending(Exception):
    """The genome exists at a source we have not acquired yet."""


def _sidecar(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".json")


def _record(path: Path, url: str) -> None:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    _sidecar(path).write_text(json.dumps({
        "url": url, "size": path.stat().st_size, "sha256": h.hexdigest(),
        "retrieved_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
    }, indent=1))


RETRY_DELAYS_S: Final = (2.0, 6.0, 18.0)


def _open(url: str, *, timeout: int):
    """urlopen with backoff: NCBI answers 429/503 under concurrent requests,
    and a rate-limited attempt is not a missing genome. 404 is returned at
    once - that is an answer."""
    import time as _time
    import urllib.error

    req = urllib.request.Request(url, headers={"User-Agent": "openbiota/0.8.4"})
    for i, delay in enumerate((*RETRY_DELAYS_S, None)):
        try:
            return urllib.request.urlopen(req, timeout=timeout)
        except urllib.error.HTTPError as exc:
            if exc.code == 404 or delay is None:
                raise
        except (urllib.error.URLError, TimeoutError, OSError):
            if delay is None:
                raise
        _time.sleep(delay * (1 + 0.5 * i))
    raise RuntimeError("unreachable")


def _download(url: str, dest: Path, *, gzip_it: bool) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    import os

    tmp = dest.with_suffix(dest.suffix + f".part{os.getpid()}")  # concurrent fetchers never share a temp file
    with _open(url, timeout=120) as resp:
        if gzip_it:
            with gzip.open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out, 1 << 20)
        else:
            with tmp.open("wb") as out:
                shutil.copyfileobj(resp, out, 1 << 20)
    tmp.replace(dest)
    _record(dest, url)
    return dest


def _ncbi(accession: str) -> str | None:
    m = _ACC.match(accession)
    if not m:
        return None
    kind, digits, _ = m.groups()
    listing = NCBI_FTP.format(kind=kind, a=digits[0:3], b=digits[3:6], c=digits[6:9])
    try:
        with _open(listing, timeout=60) as resp:
            html = resp.read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        import urllib.error

        if isinstance(exc, urllib.error.HTTPError) and exc.code == 404:
            return None
        raise
    dirs = re.findall(rf'href="({kind}_{digits}\.\d+_[^"/]+)/"', html)
    if not dirs:
        return None
    # exact version if given, else the highest
    want = accession if "." in accession else None
    chosen = next((d for d in dirs if want and d.startswith(want + "_")), None) or sorted(dirs)[-1]
    return f"{listing}{chosen}/{chosen}_genomic.fna.gz"


def normalise(genome_id: str) -> str:
    """The bare identifier: ``globdb:MOTU40_1.fa.gz_`` -> ``MOTU40_1``.

    Lane records carry ids with a source prefix, sometimes with the archive
    file's extension and a stray trailing underscore from the panel build;
    every source lookup starts from the bare id.
    """
    gid = genome_id.strip()
    if gid.startswith("file:"):
        return gid
    gid = gid.split(":", 1)[-1].strip().rstrip("_")
    for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if gid.endswith(suffix):
            gid = gid[: -len(suffix)]
            break
    return gid.rstrip("_")


def in_globdb(genome_id: str) -> bool:
    """Whether the unpacked GlobDB archive holds this genome."""
    gid = normalise(genome_id)
    return bool(gid) and not gid.startswith("file:") and _globdb_path(gid) is not None


def _place(src: Path, dest: Path, *, gzip_it: bool, url: str) -> Path:
    """Copy a local genome into the cache atomically.

    Two runs fetching the same genome at the same moment used to write the
    same file together and leave a gzip stream neither could read; each
    writer now fills its own temporary file and the last rename wins.
    """
    import os

    tmp = dest.with_suffix(dest.suffix + f".part{os.getpid()}")
    if gzip_it:
        with src.open("rb") as fh, gzip.open(tmp, "wb") as out:
            shutil.copyfileobj(fh, out, 1 << 20)
    else:
        shutil.copyfile(src, tmp)
    tmp.replace(dest)
    _record(dest, url)
    return dest


def discard(genome_id: str) -> None:
    """Drop a cached genome that turned out unreadable, so the next fetch refills it."""
    gid = normalise(genome_id)
    if gid.startswith("file:"):
        return
    dest = CACHE / f"{gid}.fna.gz"
    dest.unlink(missing_ok=True)
    _sidecar(dest).unlink(missing_ok=True)


def _fetch_local_file(spec: str, dest_dir: Path) -> Path:
    """``file:refs/...genome.fna(.gz)`` -> a gzipped copy in the cache."""
    src = Path(spec[len("file:"):])
    if not src.is_file():
        raise FileNotFoundError(spec)
    name = src.name
    for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    dest = dest_dir / f"{name}.fna.gz"
    if dest.is_file() and _sidecar(dest).is_file():
        return dest
    return _place(src, dest, gzip_it=src.suffix != ".gz", url=f"file://{src}")


def fetch(genome_id: str) -> Path:
    """Path to a gzipped FASTA for `genome_id`, fetching it if needed.

    Raises `Pending` when the only source is an archive not yet unpacked,
    and `FileNotFoundError` when no source knows the identifier.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    gid = normalise(genome_id)
    if gid.startswith("file:"):
        return _fetch_local_file(gid, CACHE)
    dest = CACHE / f"{gid}.fna.gz"
    if dest.is_file() and _sidecar(dest).is_file():
        return dest

    if _ACC.match(gid):
        # The unpacked GlobDB archive holds every GTDB representative under
        # its unversioned accession: a local copy before a network fetch.
        local = _globdb_path(gid) or _globdb_path(gid.split(".")[0])
        if local is not None:
            return _place(local, dest, gzip_it=False, url=f"file://{local}")
        url = _ncbi(gid)
        if url is None:
            raise FileNotFoundError(f"NCBI has no assembly directory for {gid}")
        return _download(url, dest, gzip_it=False)

    if gid.startswith("MGYG"):
        # GlobDB holds the UHGG representatives it kept: a local copy first,
        # then MGnify's FTP for the rest.
        local = _globdb_path(gid)
        if local is not None:
            return _place(local, dest, gzip_it=False, url=f"file://{local}")
        return _download(UHGG_SPECIES.format(prefix=gid[:-2], gid=gid), dest, gzip_it=True)

    if gid.startswith("HRGMV2_") or gid.startswith("HRGMv2_"):
        local = HRGM2_REPS / gid.replace("HRGMV2_", "HRGMv2_")
        hits = list(HRGM2_REPS.rglob(f"{gid.replace('HRGMV2_', 'HRGMv2_')}.fna")) if HRGM2_REPS.is_dir() else []
        if hits:
            return _place(hits[0], dest, gzip_it=True, url=f"file://{hits[0]}")
        raise Pending(f"{gid}: HRGM2 representative archive not unpacked")

    if gid.startswith("HROM_Genome_"):
        return _download(HROM_GENOME.format(gid=gid), dest, gzip_it=True)

    # GlobDB-native units: only from the unpacked archive.
    src = _globdb_path(gid)
    if src is not None:
        return _place(src, dest, gzip_it=False, url=f"file://{src}")
    raise Pending(f"{gid}: GlobDB genome archive not yet acquired or unpacked")


_GLOBDB_INDEX: dict[str, str] | None = None
GLOBDB_INDEX_FILE = GLOBDB_GENOMES / "index.tsv"


def _globdb_path(gid: str) -> Path | None:
    """Where the unpacked archive holds `gid`, via a one-time index.

    The archive unpacks into chunk directories (346,233 files); walking
    them for every lookup cost seconds per genome. The index maps id ->
    relative path, is written after the first walk, and is rebuilt when
    the genome directory has grown since (an interrupted unpack resumed).
    """
    global _GLOBDB_INDEX
    if not GLOBDB_GENOMES.is_dir():
        return None
    direct = GLOBDB_GENOMES / f"{gid}.fa.gz"
    if direct.is_file():
        return direct
    if _GLOBDB_INDEX is None:
        _GLOBDB_INDEX = _load_or_build_globdb_index()
    rel = _GLOBDB_INDEX.get(gid)
    if rel is None:
        # not indexed: the unpack may have progressed since; one cheap retry
        # by rebuilding only if the directory count changed
        n_dirs = str(sum(1 for p in GLOBDB_GENOMES.iterdir() if p.is_dir()))
        if n_dirs != _GLOBDB_INDEX.get("__n_dirs__"):
            _GLOBDB_INDEX = _load_or_build_globdb_index(force=True)
            rel = _GLOBDB_INDEX.get(gid)
    if rel is None:
        return None
    p = GLOBDB_GENOMES / rel
    return p if p.is_file() else None


def _load_or_build_globdb_index(force: bool = False) -> dict[str, str]:
    index: dict[str, str] = {}
    if GLOBDB_INDEX_FILE.is_file() and not force:
        with GLOBDB_INDEX_FILE.open() as fh:
            for line in fh:
                gid, _, rel = line.rstrip("\n").partition("\t")
                if gid:
                    index[gid] = rel
        if index:
            index["__n_dirs__"] = str(sum(1 for p in GLOBDB_GENOMES.iterdir() if p.is_dir()))  # type: ignore[assignment]
            return index
    n_dirs = 0
    for chunk in sorted(GLOBDB_GENOMES.iterdir()):
        if not chunk.is_dir():
            continue
        n_dirs += 1
        for f in chunk.iterdir():
            name = f.name
            if name.endswith(".fa.gz"):
                index[name[:-6]] = f"{chunk.name}/{name}"
    tmp = GLOBDB_INDEX_FILE.with_suffix(".tmp")
    with tmp.open("w") as out:
        for gid, rel in index.items():
            out.write(f"{gid}\t{rel}\n")
    tmp.replace(GLOBDB_INDEX_FILE)
    index["__n_dirs__"] = str(n_dirs)  # type: ignore[assignment]
    return index


def available(genome_id: str) -> str:
    """'cached' | 'fetchable' | 'pending_archive' | 'unknown' without fetching."""
    gid = normalise(genome_id)
    if gid.startswith("file:"):
        return "fetchable" if Path(gid[len("file:"):]).is_file() else "unknown"
    if (CACHE / f"{gid}.fna.gz").is_file():
        return "cached"
    if _ACC.match(gid) or gid.startswith(("MGYG", "HROM_Genome_")):
        return "fetchable"
    if gid.startswith(("HRGMV2_", "HRGMv2_")):
        return "fetchable" if HRGM2_REPS.is_dir() else "pending_archive"
    if _globdb_path(gid) is not None:
        return "fetchable"
    return "pending_archive"
