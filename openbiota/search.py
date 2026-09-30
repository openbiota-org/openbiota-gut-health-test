"""DIAMOND blastx wrapper: input discovery, flags, caching, resume, progress.

Two design points worth keeping:

* The ``.gz`` inputs are handed straight to DIAMOND, which reads them natively.
  That halves disk I/O and lets the uncompressed copies be deleted.
* Hit output is cached per mate file with a manifest covering the database
  fingerprint, the exact DIAMOND flags, the output field list and the input
  file's size and mtime. An interrupted run resumes at the mate boundary
  instead of restarting.
"""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import subprocess
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from openbiota.errors import DependencyError, InputError, SearchError
from openbiota.logging_util import Reporter, human_bytes, human_duration
from openbiota.qc import open_bytes
from openbiota.references import ReferenceDatabase

#: Base output fields. Nothing here is unparsed — the build spec forbids
#: requesting fields that are not used.
BASE_FIELDS: Final = ("qseqid", "sseqid", "pident", "length", "sstart", "send", "bitscore")

#: Added only when at least one selected panel runs a residue check; the gapped
#: alignment strings are what make the position walk possible.
RESIDUE_FIELDS: Final = ("qseq_gapped", "sseq_gapped")

SENSITIVITY_FLAGS: Final = {
    "default": (),
    "faster": ("--faster",),
    "fast": ("--fast",),
    "mid-sensitive": ("--mid-sensitive",),
    "sensitive": ("--sensitive",),
    "more-sensitive": ("--more-sensitive",),
    "very-sensitive": ("--very-sensitive",),
    "ultra-sensitive": ("--ultra-sensitive",),
}

_R1_PATTERNS: Final = (r"_R?1(?=[._])", r"_1$")
_SUFFIX_RE: Final = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)

@dataclass(frozen=True, slots=True)
class SearchConfig:
    """Everything that affects the DIAMOND command line."""

    diamond: str = "diamond"
    threads: int = 1
    block_size: float = 8.0
    index_chunks: int = 1
    sensitivity: str = "default"
    evalue: float = 1e-5
    max_target_seqs: int = 1
    quiet: bool = True
    extra_args: tuple[str, ...] = ()

    def sensitivity_args(self) -> tuple[str, ...]:
        try:
            return SENSITIVITY_FLAGS[self.sensitivity]
        except KeyError as exc:
            raise SearchError(
                f"unknown sensitivity {self.sensitivity!r}; "
                f"choose from {', '.join(SENSITIVITY_FLAGS)}"
            ) from exc


@dataclass(frozen=True, slots=True)
class MateInput:
    """One FASTQ file of a pair."""

    label: str
    path: Path

    @property
    def gzipped(self) -> bool:
        return self.path.suffix.lower() == ".gz"


@dataclass(frozen=True, slots=True)
class SampleInput:
    sample: str
    mates: tuple[MateInput, ...]

    @property
    def paths(self) -> tuple[Path, ...]:
        return tuple(m.path for m in self.mates)


@dataclass(frozen=True, slots=True)
class SearchOutput:
    mate: str
    query_path: Path
    hits_path: Path
    n_hit_lines: int
    elapsed_s: float
    cached: bool
    command: tuple[str, ...] = field(default=())


# --------------------------------------------------------------------------- #
# input discovery
# --------------------------------------------------------------------------- #


def strip_fastq_suffix(name: str) -> str:
    return _SUFFIX_RE.sub("", name)


def sample_name_from(path: Path) -> str:
    """The sample is named by its file: the FASTQ name minus the extension and
    the mate number, nothing else. ``SAMPLE2_A02_R1.fastq.gz`` and
    ``SAMPLE3_A03_1.fq.gz`` are the samples ``SAMPLE2_A02`` and ``SAMPLE3_A03``.
    Renaming the file renames the sample; there is no other parsing rule."""
    stem = strip_fastq_suffix(path.name)
    for pattern in (r"_R?[12]$", r"\.R?[12]$"):
        stem = re.sub(pattern, "", stem)
    return stem or "sample"


def _mate_key(path: Path) -> tuple[str, str] | None:
    """Return ``(sample, mate_label)`` if the filename encodes a mate number."""
    stem = strip_fastq_suffix(path.name)
    match = re.search(r"(_R?|\.R?)([12])$", stem)
    if not match:
        return None
    return stem[: match.start()], f"R{match.group(2)}"


def discover_samples(fastq_dir: Path, *, prefer_gzip: bool = True) -> list[SampleInput]:
    """Every mate pair in a directory, for batch runs.

    Unpaired files are ignored rather than guessed at: a lone FASTQ in a
    directory of pairs is far more likely to be a stray than a single-end
    sample, and silently screening it would produce a report nobody asked for.
    """
    if not fastq_dir.is_dir():
        raise InputError(f"FASTQ directory not found: {fastq_dir}")

    candidates = [
        p
        for p in sorted(fastq_dir.iterdir())
        if p.is_file() and _SUFFIX_RE.search(p.name) and not p.name.startswith(".")
    ]
    groups: dict[str, dict[str, dict[str, Path]]] = {}
    for path in candidates:
        key = _mate_key(path)
        if key is None:
            continue
        sample_id, mate_label = key
        slot = groups.setdefault(sample_id, {}).setdefault(mate_label, {})
        slot["gz" if path.suffix.lower() == ".gz" else "plain"] = path

    out: list[SampleInput] = []
    for sample_id in sorted(groups):
        mates: list[MateInput] = []
        for label in sorted(groups[sample_id]):
            slot = groups[sample_id][label]
            path = slot["gz"] if (prefer_gzip and "gz" in slot) else slot.get("plain") or slot["gz"]
            mates.append(MateInput(label, path.resolve()))
        if mates:
            out.append(SampleInput(sample_id, tuple(mates)))
    if not out:
        raise InputError(
            f"no mate pairs found in {fastq_dir}. Expected names like "
            f"NAME_1.fastq.gz / NAME_2.fastq.gz or NAME_R1.fastq.gz / NAME_R2.fastq.gz"
        )
    return out


def discover_sample(
    *,
    fastq_dir: Path | None = None,
    r1: Path | None = None,
    r2: Path | None = None,
    sample: str | None = None,
    prefer_gzip: bool = True,
    reporter: Reporter | None = None,
) -> SampleInput:
    """Resolve the FASTQ inputs, preferring ``.gz`` when both forms exist."""
    if r1 is not None:
        mates = [MateInput("R1", r1.resolve())]
        if r2 is not None:
            mates.append(MateInput("R2", r2.resolve()))
        for mate in mates:
            if not mate.path.is_file():
                raise InputError(f"input file not found: {mate.path}")
        return SampleInput(sample or sample_name_from(mates[0].path), tuple(mates))

    if fastq_dir is None:
        raise InputError("either --fastq-dir or --r1/--r2 must be given")
    if not fastq_dir.is_dir():
        raise InputError(f"FASTQ directory not found: {fastq_dir}")

    candidates = [
        p
        for p in sorted(fastq_dir.iterdir())
        if p.is_file() and _SUFFIX_RE.search(p.name) and not p.name.startswith(".")
    ]
    if not candidates:
        raise InputError(
            f"no FASTQ files found in {fastq_dir} (looking for *.fastq, *.fq, *.fastq.gz, *.fq.gz)"
        )

    # sample -> mate -> {gz: path, plain: path}
    groups: dict[str, dict[str, dict[str, Path]]] = {}
    unpaired: list[Path] = []
    for path in candidates:
        key = _mate_key(path)
        if key is None:
            unpaired.append(path)
            continue
        sample_id, mate_label = key
        slot = groups.setdefault(sample_id, {}).setdefault(mate_label, {})
        slot["gz" if path.suffix.lower() == ".gz" else "plain"] = path

    if not groups:
        if len(unpaired) == 1:
            only = unpaired[0]
            return SampleInput(
                sample or sample_name_from(only), (MateInput("R1", only.resolve()),)
            )
        raise InputError(
            f"could not identify mate pairs in {fastq_dir}. Expected names like "
            f"NAME_1.fastq.gz / NAME_2.fastq.gz or NAME_R1.fastq / NAME_R2.fastq. "
            f"Found: {', '.join(p.name for p in candidates[:8])}"
        )

    if sample is not None and sample in groups:
        chosen = sample
    elif len(groups) == 1:
        chosen = next(iter(groups))
    else:
        raise InputError(
            f"{fastq_dir} contains several samples ({', '.join(sorted(groups))}); "
            "pass --sample to pick one, or --r1/--r2 explicitly"
        )

    mates: list[MateInput] = []
    for label in sorted(groups[chosen]):
        slot = groups[chosen][label]
        if prefer_gzip and "gz" in slot:
            path = slot["gz"]
            if "plain" in slot and reporter is not None:
                reporter.record(
                    f"    {label}: using {path.name} in preference to the uncompressed "
                    f"{slot['plain'].name} (DIAMOND reads gzip natively; halves disk I/O)"
                )
        else:
            path = slot.get("plain") or slot["gz"]
        mates.append(MateInput(label, path.resolve()))

    return SampleInput(chosen, tuple(mates))


# --------------------------------------------------------------------------- #
# environment checks
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class PlatformInfo:
    machine: str
    system: str
    translated: bool
    cpu_brand: str
    logical_cpus: int
    memory_bytes: int

    @property
    def under_rosetta(self) -> bool:
        return self.translated


def _sysctl(name: str) -> str:
    exe = shutil.which("sysctl")
    if exe is None:
        return ""
    try:
        proc = subprocess.run(  # noqa: S603
            [exe, "-n", name], capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return proc.stdout.strip() if proc.returncode == 0 else ""


def platform_info() -> PlatformInfo:
    memory = _sysctl("hw.memsize")
    return PlatformInfo(
        machine=platform.machine(),
        system=platform.system(),
        translated=_sysctl("sysctl.proc_translated") == "1",
        cpu_brand=_sysctl("machdep.cpu.brand_string") or platform.processor() or "unknown",
        logical_cpus=os.cpu_count() or 1,
        memory_bytes=int(memory) if memory.isdigit() else 0,
    )


def diamond_binary_arch(diamond: str = "diamond") -> str:
    """Best-effort architecture readout for the DIAMOND executable."""
    exe = shutil.which(diamond)
    file_tool = shutil.which("file")
    if exe is None or file_tool is None:
        return "unknown"
    try:
        proc = subprocess.run(  # noqa: S603
            [file_tool, "-b", exe], capture_output=True, text=True, timeout=15, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    text = proc.stdout.strip()
    for arch in ("arm64", "x86_64", "arm64e", "i386"):
        if arch in text:
            return arch
    return text[:60] or "unknown"


def check_platform(reporter: Reporter, diamond: str = "diamond") -> PlatformInfo:
    info = platform_info()
    reporter.record(
        f"    platform: {info.system}/{info.machine}, {info.logical_cpus} logical CPUs, "
        f"{human_bytes(info.memory_bytes)} RAM, {info.cpu_brand}"
    )
    if info.under_rosetta:
        reporter.warn(
            "this Python process is running translated under Rosetta — expect roughly 2x "
            "slowdown. Reinstall a native arm64 Python and rebuild the environment."
        )
    arch = diamond_binary_arch(diamond)
    if info.machine == "arm64" and arch == "x86_64":
        reporter.warn(
            f"the DIAMOND binary is {arch} on an {info.machine} host, so it runs under Rosetta "
            "and costs roughly 2x. Install a native build: `brew install diamond`, or compile "
            "from source. Bioconda's osx-64 package is the usual cause."
        )
    else:
        reporter.record(f"    diamond binary architecture: {arch}")
    return info


# --------------------------------------------------------------------------- #
# subsampling
# --------------------------------------------------------------------------- #


def write_subsample(source: Path, dest: Path, n_reads: int, reporter: Reporter) -> Path:
    """Materialise the first ``n_reads`` records of ``source`` as plain FASTQ."""
    if dest.is_file():
        marker = dest.with_suffix(dest.suffix + ".done")
        if marker.is_file() and marker.read_text(encoding="utf-8").strip() == str(n_reads):
            reporter.record(f"    reusing cached subsample {dest.name} ({n_reads:,} reads)")
            return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    wanted = n_reads * 4
    # Streamed through the system decompressor rather than zlib, and only the
    # requested prefix is inflated — the pipe is closed as soon as enough
    # records have been taken, so the rest of the file is never touched.
    written = 0
    with open_bytes(source, partial_ok=True) as src, dest.open("wb") as out:
        for line in src:
            out.write(line)
            written += 1
            if written >= wanted:
                break
    if written == 0:
        raise InputError(f"{source} is empty")
    if written % 4 != 0:
        raise InputError(
            f"{source}: truncated FASTQ — read {written} lines, which is not a multiple of 4"
        )
    dest.with_suffix(dest.suffix + ".done").write_text(str(n_reads), encoding="utf-8")
    reporter.record(f"    subsample written: {dest.name} ({written // 4:,} reads)")
    return dest


# --------------------------------------------------------------------------- #
# search
# --------------------------------------------------------------------------- #


def output_fields(*, with_residues: bool) -> tuple[str, ...]:
    return BASE_FIELDS + (RESIDUE_FIELDS if with_residues else ())


def build_command(
    *,
    config: SearchConfig,
    db_path: Path,
    query_path: Path,
    out_path: Path,
    fields: Sequence[str],
) -> tuple[str, ...]:
    exe = shutil.which(config.diamond)
    if exe is None:
        raise DependencyError(
            f"{config.diamond!r} not found on PATH. Run `openbiota doctor` for install guidance."
        )
    cmd: list[str] = [
        exe, "blastx",
        "--db", str(db_path),
        "--query", str(query_path),
        "--out", str(out_path),
        # Best hit only. This is both the specificity mechanism (a read is
        # credited to a target only when that target beats every decoy) and a
        # large speedup.
        "--max-target-seqs", str(config.max_target_seqs),
        # Single in-memory index chunk: the reference database is a few tens of
        # thousands of proteins and there is ample RAM.
        "--index-chunks", str(config.index_chunks),
        "--block-size", f"{config.block_size:g}",
        "--evalue", f"{config.evalue:g}",
        "--threads", str(config.threads),
        "--outfmt", "6", *fields,
    ]
    cmd.extend(config.sensitivity_args())
    if config.quiet:
        cmd.append("--quiet")
    cmd.extend(config.extra_args)
    return tuple(cmd)


_PERFORMANCE_ONLY_FLAGS = frozenset({"--threads", "-p"})
#: Path-valued flags: the path is recorded separately (query) or is an output
#: location (db, out). Dropped by value, not by shape, so a relative path from
#: a batch run and an absolute path from a direct run key the same cache.
_PATH_FLAGS = frozenset({"--db", "-d", "--query", "-q", "--out", "-o"})


def _significant_args(command: Sequence[str]) -> list[str]:
    """Drop paths and performance-only flags (with their values)."""
    kept: list[str] = []
    skip_value = False
    drop_path = False
    for arg in command:
        if skip_value:
            skip_value = False
            continue
        if drop_path:
            drop_path = False
            if not arg.startswith("-"):
                # the path itself; an already-normalised list has none here
                continue
        if arg in _PERFORMANCE_ONLY_FLAGS:
            skip_value = True
            continue
        if arg in _PATH_FLAGS:
            kept.append(arg)
            drop_path = True
            continue
        if arg.startswith("/"):
            continue
        kept.append(arg)
    return kept


def _manifest_payload(
    *,
    command: Sequence[str],
    fields: Sequence[str],
    db_fingerprint: str,
    query_path: Path,
    subsample: int | None,
) -> dict[str, object]:
    stat = query_path.stat()
    return {
        "db_fingerprint": db_fingerprint,
        "fields": list(fields),
        # The db path and out path are absolute and vary between runs, so
        # exclude them. So is the thread count: it changes how fast DIAMOND
        # runs, not what it writes, and a laptop and a workstation screening
        # the same reads must be able to share a cache. Everything else about
        # the command is significant.
        "command": _significant_args(command),
        "query": str(query_path),
        "query_size": stat.st_size,
        "query_mtime": int(stat.st_mtime),
        "subsample": subsample,
    }


def count_lines(path: Path, chunk_size: int = 1 << 22) -> int:
    total = 0
    with path.open("rb") as fh:
        while chunk := fh.read(chunk_size):
            total += chunk.count(b"\n")
    return total


def run_search(
    *,
    database: ReferenceDatabase,
    sample_input: SampleInput,
    out_dir: Path,
    config: SearchConfig,
    fields: Sequence[str],
    reporter: Reporter,
    subsample: int | None = None,
    force: bool = False,
) -> list[SearchOutput]:
    """Run (or reuse) one DIAMOND blastx pass per mate file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[SearchOutput] = []

    for mate in sample_input.mates:
        query_path = mate.path
        if subsample is not None:
            query_path = write_subsample(
                mate.path,
                out_dir / "alignments" / f"subsample_{subsample}_{mate.label}.fastq",
                subsample,
                reporter,
            )

        # Raw alignments go in their own subfolder: there are two of them plus
        # eighteen per-panel views, and burying the four files a human actually
        # opens among twenty machine artefacts makes the directory unusable.
        hits_dir = out_dir / "alignments"
        hits_dir.mkdir(parents=True, exist_ok=True)
        hits_path = hits_dir / f"hits_all_{mate.label}.tsv"
        manifest_path = hits_path.with_suffix(".tsv.done.json")
        command = build_command(
            config=config,
            db_path=database.dmnd_path,
            query_path=query_path,
            out_path=hits_path,
            fields=fields,
        )
        expected = _manifest_payload(
            command=command,
            fields=fields,
            db_fingerprint=database.fingerprint,
            query_path=query_path,
            subsample=subsample,
        )

        if not force and hits_path.is_file() and manifest_path.is_file():
            try:
                recorded = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                recorded = None
            if isinstance(recorded, dict) and isinstance(recorded.get("command"), list):
                # Manifests written before performance flags were excluded
                # still carry them; normalise so an old cache stays valid.
                recorded["command"] = _significant_args(recorded["command"])
            if isinstance(recorded, dict) and {
                k: recorded.get(k) for k in expected
            } == expected:
                n_lines = int(recorded.get("n_hit_lines") or count_lines(hits_path))
                reporter.ok(
                    f"{mate.label}: reusing cached DIAMOND output "
                    f"({n_lines:,} hit lines, {human_bytes(hits_path.stat().st_size)})"
                )
                results.append(
                    SearchOutput(
                        mate=mate.label,
                        query_path=query_path,
                        hits_path=hits_path,
                        n_hit_lines=n_lines,
                        elapsed_s=0.0,
                        cached=True,
                        command=command,
                    )
                )
                continue

        # Write to a temporary path so an interrupted run never leaves a
        # partial hits file that a later run could mistake for complete.
        partial_path = hits_path.with_suffix(".tsv.partial")
        partial_command = build_command(
            config=config,
            db_path=database.dmnd_path,
            query_path=query_path,
            out_path=partial_path,
            fields=fields,
        )
        manifest_path.unlink(missing_ok=True)

        size_label = human_bytes(query_path.stat().st_size)
        reporter.step(
            f"{mate.label}: diamond blastx on {query_path.name} ({size_label}), "
            f"{config.threads} threads, sensitivity={config.sensitivity}"
        )
        reporter.record(f"    $ {' '.join(partial_command)}")

        def probe(p: Path = partial_path) -> str:
            try:
                return f"{human_bytes(p.stat().st_size)} of hits written"
            except FileNotFoundError:
                return "no hits written yet"

        started = time.monotonic()
        with reporter.heartbeat(f"{mate.label} diamond blastx", probe):
            try:
                proc = subprocess.run(  # noqa: S603 — fixed argv, no shell
                    partial_command, capture_output=True, text=True, check=True
                )
            except subprocess.CalledProcessError as exc:
                partial_path.unlink(missing_ok=True)
                raise SearchError(
                    "diamond blastx failed.\n"
                    f"  command: {' '.join(partial_command)}\n"
                    f"  exit code: {exc.returncode}\n"
                    f"  stderr: {(exc.stderr or '').strip()[-2000:]}"
                ) from exc
        elapsed = time.monotonic() - started

        stderr = (proc.stderr or "").strip()
        if stderr:
            for line in stderr.splitlines():
                reporter.record(f"    diamond: {line}")

        if not partial_path.is_file():
            raise SearchError(
                f"diamond blastx reported success but produced no output file at {partial_path}"
            )
        partial_path.replace(hits_path)
        n_lines = count_lines(hits_path)
        manifest_path.write_text(
            json.dumps({**expected, "n_hit_lines": n_lines, "elapsed_s": elapsed}, indent=2),
            encoding="utf-8",
        )
        reporter.ok(
            f"{mate.label}: {n_lines:,} hit lines in {human_duration(elapsed)} "
            f"({human_bytes(hits_path.stat().st_size)})"
        )
        results.append(
            SearchOutput(
                mate=mate.label,
                query_path=query_path,
                hits_path=hits_path,
                n_hit_lines=n_lines,
                elapsed_s=elapsed,
                cached=False,
                command=partial_command,
            )
        )

    if not results:
        raise SearchError("no mate files were searched")
    return results
