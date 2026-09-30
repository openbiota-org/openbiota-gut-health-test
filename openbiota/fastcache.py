"""A disk cache for parsed configuration, keyed on what could change it.

Every report run re-parses the same YAML: 485 pathogen seeds, 33 disease
profiles, 97 organisms of agent advice, 25 panels, and a 14 MB taxonomic
cohort. None of it depends on the sample. Measured, it costs about 2.9 s a
run — small against a cold run's twenty minutes, but paid on every one of
them, and pure waste when the files have not changed.

The cache is a pickle of the parsed object, valid only while a fingerprint
matches. The fingerprint covers everything that could make the cached object
wrong:

* the source files — path, size and modification time of every one;
* the code that parsed them — a hash of the parsing module's source, so a
  change to a dataclass or loader invalidates the pickle instead of
  returning an object of the old shape;
* the package version.

A stale or unreadable cache is never an error: the loader simply runs and
the cache is rewritten. Set ``OPENBIOTA_NO_FASTCACHE=1`` to bypass it.
"""

from __future__ import annotations

import hashlib
import inspect
import os
import pickle
import sys
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Final, TypeVar

from openbiota import __version__

T = TypeVar("T")

#: Where pickles live. `refs/` is the repository's generated-state directory
#: and is ignored by git, so a cache here can never be committed.
DEFAULT_CACHE_DIR: Final = Path(__file__).resolve().parents[1] / "refs" / "cache" / "parsed"

_DISABLE_ENV: Final = "OPENBIOTA_NO_FASTCACHE"


def _module_source_hash(module_name: str) -> str:
    """Hash of a module's source.

    Parsed objects are dataclasses; if the module defining them changes, an
    old pickle may unpickle into an object with the wrong fields or
    defaults. Hashing the whole module is coarse but safe.
    """
    module = sys.modules.get(module_name)
    if module is None:
        try:
            module = __import__(module_name, fromlist=["_"])
        except ImportError:
            return f"{module_name}:no-module"
    try:
        return hashlib.sha256(inspect.getsource(module).encode("utf-8")).hexdigest()[:16]
    except (OSError, TypeError):
        return f"{module_name}:no-source"


def fingerprint(
    sources: Iterable[Path],
    build: Callable[..., object],
    code_modules: Iterable[str] = (),
) -> str:
    """One hash over the source files, the parsing code and the version.

    `code_modules` names every module whose classes end up inside the cached
    object; the builder's own module is always included.
    """
    digest = hashlib.sha256()
    digest.update(__version__.encode("utf-8"))
    modules = {getattr(build, "__module__", "") or "", *code_modules}
    for module_name in sorted(m for m in modules if m):
        digest.update(_module_source_hash(module_name).encode("utf-8"))
    for path in sorted(Path(p) for p in sources):
        try:
            stat = path.stat()
        except OSError:
            digest.update(f"{path}:missing".encode())
            continue
        digest.update(f"{path}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return digest.hexdigest()[:24]


def cached_load(
    name: str,
    sources: Iterable[Path],
    build: Callable[[], T],
    *,
    code_modules: Iterable[str] = (),
    cache_dir: Path | None = None,
) -> T:
    """Return `build()`, from disk when nothing it depends on has changed.

    `name` distinguishes caches for different loaders; `sources` are the
    files whose change must invalidate the result; `code_modules` are the
    modules whose classes the result contains.
    """
    if os.environ.get(_DISABLE_ENV):
        return build()

    sources = list(sources)
    directory = cache_dir or DEFAULT_CACHE_DIR
    key = fingerprint(sources, build, code_modules)
    path = directory / f"{name}.{key}.pickle"

    if path.is_file():
        try:
            with path.open("rb") as fh:
                return pickle.load(fh)  # noqa: S301 - our own pickle, in our own refs/ dir
        except Exception:  # noqa: BLE001 - any unreadable cache is just rebuilt
            with _suppress_oserror():
                path.unlink()

    result = build()
    try:
        directory.mkdir(parents=True, exist_ok=True)
        # Write beside, then rename: a reader never sees a half-written file.
        fd, tmp = tempfile.mkstemp(dir=directory, prefix=f".{name}.", suffix=".tmp")
        with os.fdopen(fd, "wb") as fh:
            pickle.dump(result, fh, protocol=pickle.HIGHEST_PROTOCOL)
        Path(tmp).replace(path)
        _sweep_stale(directory, name, keep=path)
    except Exception:  # noqa: BLE001 - a failed cache write must never fail a run
        pass
    return result


def _sweep_stale(directory: Path, name: str, *, keep: Path) -> None:
    """Remove older fingerprints of the same cache so the directory stays small."""
    for old in directory.glob(f"{name}.*.pickle"):
        if old != keep:
            with _suppress_oserror():
                old.unlink()


class _suppress_oserror:
    def __enter__(self) -> None:
        return None

    def __exit__(self, exc_type: object, *_: object) -> bool:
        return exc_type is not None and issubclass(exc_type, OSError)  # type: ignore[arg-type]


def self_test() -> int:
    """The contract: same inputs hit, any change misses, corruption heals."""
    import time

    checks = 0
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "a.yaml"
        source.write_text("x: 1", encoding="utf-8")
        calls = {"n": 0}

        def build() -> dict[str, int]:
            calls["n"] += 1
            return {"x": calls["n"]}

        first = cached_load("t", [source], build, cache_dir=root / "c")
        second = cached_load("t", [source], build, cache_dir=root / "c")
        assert first == second == {"x": 1} and calls["n"] == 1, "second call should hit"
        checks += 1

        time.sleep(0.01)
        source.write_text("x: 2", encoding="utf-8")  # new mtime and size
        third = cached_load("t", [source], build, cache_dir=root / "c")
        assert third == {"x": 2} and calls["n"] == 2, "a changed source must miss"
        checks += 1

        # Corruption heals rather than raises.
        for pickle_file in (root / "c").glob("t.*.pickle"):
            pickle_file.write_bytes(b"not a pickle")
        fourth = cached_load("t", [source], build, cache_dir=root / "c")
        assert fourth == {"x": 3} and calls["n"] == 3
        checks += 1

        # Only one fingerprint kept per name.
        assert len(list((root / "c").glob("t.*.pickle"))) == 1
        checks += 1

        os.environ[_DISABLE_ENV] = "1"
        try:
            fifth = cached_load("t", [source], build, cache_dir=root / "c")
            assert fifth == {"x": 4} and calls["n"] == 4, "bypass must call build"
        finally:
            del os.environ[_DISABLE_ENV]
        checks += 1
    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} fastcache checks passed")
