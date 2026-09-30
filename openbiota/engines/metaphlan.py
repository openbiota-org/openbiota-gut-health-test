"""Taxonomic engine: MetaPhlAn with host-read removal and version pinning.

Why a second engine
-------------------
The DIAMOND engine answers "how much of gene X"; disease profiles need "how
much of species Y", which requires marker-gene profiling. The `rpoB`-derived
composition already in the report is a useful by-product but is marker
limited — reliable at phylum, indicative at genus. Profile scoring needs
species.

Why MetaPhlAn 3 rather than 4
-----------------------------
The reference cohort comes from curatedMetagenomicData, whose profiles were
computed with MetaPhlAn 3 against the CHOCOPhlAn 201901 markers. Species
abundances from different profiler versions are not comparable: MetaPhlAn 4
reorganised the taxonomy into SGBs, renamed genera (*Bacteroides vulgatus*
became *Phocaeicola vulgatus*), and changed detection behaviour. A sample
profiled with one version and scored against a reference built with another
would measure the version difference. So the sample is profiled with the same
version as the reference, and both are recorded. The GMWI2 dysbiosis anchor
was also trained on MetaPhlAn 3 output, which settles it.

Host reads
----------
The gene-capacity pipeline needs no host removal because human reads carry
neither pathway genes nor `rpoB` and cancel in the ratio. MetaPhlAn has no
such property: host reads inflate the denominator and, in small numbers, can
hit markers. Reads are filtered against GRCh38 with Bowtie2 before profiling,
the host fraction is reported as a QC field, and the host-aligned reads are
never written to disk — only the non-host FASTQ is kept.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import DependencyError, OpenBiotaError
from openbiota.logging_util import Reporter, human_duration

#: Database index this pipeline is pinned to. Changing it invalidates every
#: cached profile and every reference comparison; see the module docstring.
PINNED_INDEX: Final = "mpa_v31_CHOCOPhlAn_201901"

#: Profiler family the reference cohort was built with. Compared at run time.
PINNED_PROFILER_FAMILY: Final = "MetaPhlAn 3"

#: Human reference for host filtering. GRCh38 primary assembly; the Bowtie2
#: index is built once and cached. CHM13 and PhiX can be appended by placing
#: additional FASTA files in the host directory before the index is built.
GRCH38_URL: Final = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
PHIX_URL: Final = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/"
    "GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
)

#: Cache schema. Bump when the output layout changes.
CACHE_VERSION: Final = 2


@dataclass(frozen=True, slots=True)
class MetaphlanTools:
    metaphlan: str
    bowtie2: str
    bowtie2_build: str
    version: str

    @property
    def family(self) -> str:
        major = self.version.split(".")[0] if self.version else "?"
        return f"MetaPhlAn {major}"


def locate_tools(*, explicit: str | None = None) -> MetaphlanTools:
    """Find a MetaPhlAn 3 executable and Bowtie2.

    Search order: an explicit path, ``$OPENBIOTA_METAPHLAN``, the project's
    dedicated ``.venv-mpa3`` environment, then ``PATH``. MetaPhlAn 3 and 4
    share a package name and cannot coexist in one environment, which is why
    the dedicated environment exists.
    """
    candidates: list[str] = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get("OPENBIOTA_METAPHLAN"):
        candidates.append(os.environ["OPENBIOTA_METAPHLAN"])
    repo_root = Path(__file__).resolve().parents[2]
    candidates.append(str(repo_root / ".venv-mpa3" / "bin" / "metaphlan"))
    on_path = shutil.which("metaphlan")
    if on_path:
        candidates.append(on_path)

    metaphlan = next((c for c in candidates if Path(c).is_file()), None)
    if metaphlan is None:
        raise DependencyError(
            "MetaPhlAn 3 not found. Create the dedicated environment with\n"
            "    python3 -m venv .venv-mpa3 && .venv-mpa3/bin/pip install 'metaphlan==3.1.0'\n"
            "or set OPENBIOTA_METAPHLAN to a metaphlan executable."
        )
    version = _metaphlan_version(metaphlan)
    if not version.startswith("3."):
        raise DependencyError(
            f"MetaPhlAn {version} found at {metaphlan}, but this pipeline is pinned to "
            f"{PINNED_PROFILER_FAMILY} so that the sample is profiled with the same version "
            "as the reference cohort. Install metaphlan==3.1.0 in .venv-mpa3."
        )

    bowtie2 = shutil.which("bowtie2")
    bowtie2_build = shutil.which("bowtie2-build")
    if not bowtie2 or not bowtie2_build:
        raise DependencyError(
            "bowtie2 and bowtie2-build are required for host filtering and MetaPhlAn; "
            "install with `brew install bowtie2` or `apt install bowtie2`."
        )
    return MetaphlanTools(
        metaphlan=metaphlan, bowtie2=bowtie2, bowtie2_build=bowtie2_build, version=version
    )


def _env_for(tools_or_exe: MetaphlanTools | str) -> dict[str, str]:
    """Subprocess environment with the MetaPhlAn bin directory on PATH.

    MetaPhlAn 3 launches its own helper scripts (``read_fastx.py``) by name
    through the shell, so running the venv's executable directly — without
    activating the venv — fails unless its ``bin/`` is on PATH.
    """
    executable = tools_or_exe.metaphlan if isinstance(tools_or_exe, MetaphlanTools) else tools_or_exe
    env = dict(os.environ)
    bin_dir = str(Path(executable).resolve().parent)
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def _metaphlan_version(executable: str) -> str:
    try:
        out = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=120, check=False,
            env=_env_for(executable),
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise DependencyError(f"cannot run {executable}: {exc}") from exc
    text = (out.stdout + out.stderr).strip()
    for token in text.replace("(", " ").split():
        if token[:1].isdigit() and token.count(".") >= 1:
            return token
    return text.split()[-1] if text else ""


# --------------------------------------------------------------------------- #
# database
# --------------------------------------------------------------------------- #


def database_ready(db_dir: Path, index: str = PINNED_INDEX) -> bool:
    """True when the Bowtie2 index for the pinned database exists."""
    needed = [db_dir / f"{index}.{n}.bt2" for n in ("1", "2", "3", "4", "rev.1", "rev.2")]
    large = [db_dir / f"{index}.{n}.bt2l" for n in ("1", "2", "3", "4", "rev.1", "rev.2")]
    return (all(p.is_file() for p in needed) or all(p.is_file() for p in large)) and (
        db_dir / f"{index}.pkl"
    ).is_file()


def install_database(
    tools: MetaphlanTools, db_dir: Path, reporter: Reporter, *, index: str = PINNED_INDEX,
    threads: int = 8,
) -> None:
    """Download and index the pinned MetaPhlAn database."""
    db_dir.mkdir(parents=True, exist_ok=True)
    if database_ready(db_dir, index):
        return
    reporter.info(f"  installing MetaPhlAn database {index} into {db_dir} (one-time)")
    command = [
        tools.metaphlan, "--install", "--index", index,
        "--bowtie2db", str(db_dir), "--nproc", str(threads),
    ]
    started = time.monotonic()
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, env=_env_for(tools)
    )
    if result.returncode != 0 or not database_ready(db_dir, index):
        raise OpenBiotaError(
            f"MetaPhlAn database install failed (exit {result.returncode}):\n"
            f"{result.stderr[-2000:]}"
        )
    reporter.ok(f"MetaPhlAn database ready in {human_duration(time.monotonic() - started)}")


# --------------------------------------------------------------------------- #
# host filtering
# --------------------------------------------------------------------------- #


def host_index_ready(host_dir: Path) -> bool:
    return all((host_dir / f"host.{n}.bt2l").is_file() for n in ("1", "2", "3", "4", "rev.1", "rev.2")) or all(
        (host_dir / f"host.{n}.bt2").is_file() for n in ("1", "2", "3", "4", "rev.1", "rev.2")
    )


def build_host_index(
    tools: MetaphlanTools, host_dir: Path, reporter: Reporter, *, threads: int = 8
) -> None:
    """Download GRCh38 + PhiX and build the Bowtie2 host index (one-time, ~1 h)."""
    from openbiota.net import download_to

    host_dir.mkdir(parents=True, exist_ok=True)
    if host_index_ready(host_dir):
        return
    fasta = host_dir / "host.fna"
    if not fasta.is_file():
        reporter.info("  downloading GRCh38 (~900 MB) and PhiX for host filtering (one-time)")
        parts = []
        for url, name in ((GRCH38_URL, "grch38.fna.gz"), (PHIX_URL, "phix.fna.gz")):
            dest = host_dir / name
            if not dest.is_file():
                download_to(url, dest)
            parts.append(dest)
        with fasta.open("wb") as out:
            for part in parts:
                with gzip.open(part, "rb") as handle:
                    shutil.copyfileobj(handle, out, 1 << 22)
    reporter.info("  building Bowtie2 host index (one-time; this takes a while)")
    started = time.monotonic()
    result = subprocess.run(
        [tools.bowtie2_build, "--threads", str(threads), "--large-index",
         str(fasta), str(host_dir / "host")],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0 or not host_index_ready(host_dir):
        raise OpenBiotaError(f"bowtie2-build failed:\n{result.stderr[-2000:]}")
    fasta.unlink(missing_ok=True)
    reporter.ok(f"host index built in {human_duration(time.monotonic() - started)}")


#: Pairs aligned in the host probe before the full pass is attempted. At this
#: size a sample with even 0.002% human DNA has a >99% chance of showing at
#: least one host pair, so zero hits is strong evidence the provider stripped
#: human reads upstream — the routine case for consumer stool kits.
HOST_PROBE_PAIRS: Final = 250_000

#: Host filter outcomes. ``filtered`` is the ordinary full pass;
#: ``not_assessable_upstream_removed`` records that the probe found no human
#: reads, so the host fraction of the original specimen cannot be measured
#: here and the full pass was skipped as redundant (spec v04 §1.3, §4.1).
HOST_STATUS_FILTERED: Final = "filtered"
HOST_STATUS_UPSTREAM_REMOVED: Final = "not_assessable_upstream_removed"


@dataclass(frozen=True, slots=True)
class HostFilterResult:
    nonhost_r1: Path
    nonhost_r2: Path | None
    total_pairs: int
    host_pairs: int
    cached: bool
    status: str = HOST_STATUS_FILTERED
    #: Pairs examined by the probe (0 when the full pass ran without one).
    probe_pairs: int = 0
    probe_host_pairs: int = 0

    @property
    def upstream_removed(self) -> bool:
        return self.status == HOST_STATUS_UPSTREAM_REMOVED

    @property
    def host_fraction(self) -> float | None:
        """Measured host fraction, or ``None`` when it is not assessable."""
        if self.upstream_removed:
            return None
        return self.host_pairs / self.total_pairs if self.total_pairs else 0.0

    @property
    def host_fraction_upper_bound(self) -> float | None:
        """One-sided 95% upper bound on the host fraction after a zero-hit
        probe (rule of three: 3/n)."""
        if self.probe_pairs and self.probe_host_pairs == 0:
            return 3.0 / self.probe_pairs
        return None

    @property
    def nonhost_pairs(self) -> int:
        return self.total_pairs - self.host_pairs

    def describe(self) -> str:
        """One clause for prose: what the host filter found."""
        if self.upstream_removed:
            return (
                f"none of the first {self.probe_pairs:,} pairs were human — human reads were "
                "removed by the sequencing provider before delivery, so the host fraction of "
                "the original specimen is not assessable here"
            )
        return f"{self.host_pairs:,} of {self.total_pairs:,} pairs ({self.host_fraction:.2%})"

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "total_pairs": self.total_pairs,
            "host_pairs": self.host_pairs,
            "nonhost_pairs": self.nonhost_pairs,
            "host_fraction": None if self.host_fraction is None else round(self.host_fraction, 5),
            "host_fraction_upper_bound_95": (
                None if self.host_fraction_upper_bound is None
                else round(self.host_fraction_upper_bound, 7)
            ),
            "probe_pairs": self.probe_pairs,
            "probe_host_pairs": self.probe_host_pairs,
            "reference": "GRCh38 no-alt analysis set + PhiX174",
        }


def _run_host_bowtie2(
    tools: MetaphlanTools,
    *,
    r1: Path,
    r2: Path | None,
    host_dir: Path,
    work_dir: Path,
    threads: int,
    compressed: bool,
    probe: int | None = None,
) -> str:
    """One Bowtie2 pass writing only the unaligned (non-host) reads; returns stderr.

    With ``probe`` set, only the first ``probe`` pairs are aligned and nothing
    is written: the summary alone is wanted."""
    command = [
        tools.bowtie2, "-p", str(threads), "--very-fast", "-x", str(host_dir / "host"),
    ]
    suffix = ".fastq.gz" if compressed else ".fastq"
    if probe is not None:
        command += ["-u", str(probe)]
        if r2 is not None:
            command += ["-1", str(r1), "-2", str(r2)]
        else:
            command += ["-U", str(r1)]
    elif r2 is not None:
        flag = "--un-conc-gz" if compressed else "--un-conc"
        command += ["-1", str(r1), "-2", str(r2), flag, str(work_dir / f"nonhost.%{suffix}")]
    else:
        flag = "--un-gz" if compressed else "--un"
        command += ["-U", str(r1), flag, str(work_dir / f"nonhost.1{suffix}")]
    command += ["-S", os.devnull]
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, env=_env_for(tools)
    )
    if result.returncode != 0:
        raise OpenBiotaError(f"bowtie2 host filtering failed:\n{result.stderr[-2000:]}")
    return result.stderr


def _gzip_intact(path: Path) -> bool:
    """True when ``path`` decompresses end-to-end without error.

    Python's reader rejects trailing garbage after the last member, which is
    exactly the failure MetaPhlAn's FASTQ reader trips over.
    """
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1 << 22):
                pass
    except (OSError, EOFError, gzip.BadGzipFile):
        return False
    return True


def _gzip_file(src: Path, dest: Path) -> None:
    """Compress ``src`` to ``dest`` (level 1: these are scratch reads)."""
    with src.open("rb") as inp, gzip.open(dest, "wb", compresslevel=1) as out:
        shutil.copyfileobj(inp, out, 1 << 22)


def filter_host(
    tools: MetaphlanTools,
    *,
    r1: Path,
    r2: Path | None,
    host_dir: Path,
    work_dir: Path,
    reporter: Reporter,
    threads: int = 8,
    probe_pairs: int = HOST_PROBE_PAIRS,
    known_total_pairs: int | None = None,
) -> HostFilterResult:
    """Remove read pairs that align to the host reference.

    Only the non-host reads are written. Bowtie2's ``--un-conc-gz`` writes
    pairs where neither mate aligned; the SAM stream itself goes to /dev/null,
    so host-aligned sequence never touches disk.

    A probe of the first ``probe_pairs`` pairs runs first. When it finds no
    host read at all, the provider removed human reads before delivery: the
    full pass would only rewrite every read unchanged (five minutes on a
    30M-pair sample for nothing), so it is skipped, the original files are
    used as the non-host reads and the result says the host fraction is
    ``not_assessable_upstream_removed``. ``probe_pairs=0`` disables the probe.
    ``known_total_pairs`` (from fastp or the input count) fills in the pair
    total the skipped pass would have reported.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    stamp = work_dir / "host_filter.json"
    out1 = work_dir / "nonhost.1.fastq.gz"
    out2 = work_dir / "nonhost.2.fastq.gz"
    fingerprint = _fingerprint(r1, r2)
    if stamp.is_file():
        payload = json.loads(stamp.read_text())
        fresh = payload.get("cache_version") == CACHE_VERSION and payload.get("inputs") == fingerprint
        if fresh and payload.get("status") == HOST_STATUS_UPSTREAM_REMOVED:
            return HostFilterResult(
                nonhost_r1=r1, nonhost_r2=r2,
                total_pairs=known_total_pairs or payload.get("total_pairs", 0),
                host_pairs=0, cached=True, status=HOST_STATUS_UPSTREAM_REMOVED,
                probe_pairs=payload.get("probe_pairs", 0), probe_host_pairs=0,
            )
        if (
            fresh
            and out1.is_file() and (r2 is None or out2.is_file())
            and _gzip_intact(out1)
            and (r2 is None or _gzip_intact(out2))
        ):
            return HostFilterResult(
                nonhost_r1=out1,
                nonhost_r2=out2 if r2 else None,
                total_pairs=payload["total_pairs"],
                host_pairs=payload["host_pairs"],
                cached=True,
                probe_pairs=payload.get("probe_pairs", 0),
                probe_host_pairs=payload.get("probe_host_pairs", 0),
            )

    started = time.monotonic()
    probe_host = 0
    if probe_pairs > 0:
        reporter.info(f"  probing the first {probe_pairs:,} pairs for host reads")
        summary = _run_host_bowtie2(tools, r1=r1, r2=r2, host_dir=host_dir, work_dir=work_dir,
                                    threads=threads, compressed=True, probe=probe_pairs)
        probe_total, probe_host = _parse_bowtie2_summary(summary)
        if probe_total and probe_host == 0:
            elapsed = time.monotonic() - started
            stamp.write_text(json.dumps({
                "cache_version": CACHE_VERSION, "inputs": fingerprint,
                "status": HOST_STATUS_UPSTREAM_REMOVED,
                "probe_pairs": probe_total, "probe_host_pairs": 0,
                "total_pairs": known_total_pairs or 0, "elapsed_s": round(elapsed, 1),
            }))
            reporter.ok(
                f"host filtering: 0 of the first {probe_total:,} pairs were human — reads were "
                f"host-depleted before delivery; full pass skipped ({human_duration(elapsed)})"
            )
            return HostFilterResult(
                nonhost_r1=r1, nonhost_r2=r2, total_pairs=known_total_pairs or 0,
                host_pairs=0, cached=False, status=HOST_STATUS_UPSTREAM_REMOVED,
                probe_pairs=probe_total, probe_host_pairs=0,
            )
        reporter.record(f"    probe: {probe_host:,} of {probe_total:,} pairs were host; running the full pass")
    reporter.info("  removing host reads (Bowtie2 vs GRCh38 + PhiX)")
    stderr = _run_host_bowtie2(tools, r1=r1, r2=r2, host_dir=host_dir, work_dir=work_dir,
                               threads=threads, compressed=True)
    outputs = [out1] + ([out2] if r2 is not None else [])
    bad = [p for p in outputs if not _gzip_intact(p)]
    if bad:
        # Bowtie2's built-in gzip writer has produced files with trailing
        # garbage under multi-threaded runs (seen intermittently with 2.5.5).
        # Rather than hand MetaPhlAn a broken stream, redo the pass writing
        # plain FASTQ and compress it here, where the result can be trusted.
        reporter.warn(
            f"  bowtie2 wrote a corrupt gzip stream ({', '.join(p.name for p in bad)}); "
            "re-running with plain output and compressing in-process"
        )
        for p in outputs:
            p.unlink(missing_ok=True)
        stderr = _run_host_bowtie2(tools, r1=r1, r2=r2, host_dir=host_dir, work_dir=work_dir,
                                   threads=threads, compressed=False)
        for p in outputs:
            plain = p.with_suffix("")  # nonhost.1.fastq
            _gzip_file(plain, p)
            plain.unlink(missing_ok=True)
        still_bad = [p for p in outputs if not _gzip_intact(p)]
        if still_bad:
            raise OpenBiotaError(
                f"host-filtered reads could not be written intact: {[str(p) for p in still_bad]}"
            )

    total, host = _parse_bowtie2_summary(stderr)
    stamp.write_text(
        json.dumps(
            {
                "cache_version": CACHE_VERSION,
                "inputs": _fingerprint(r1, r2),
                "status": HOST_STATUS_FILTERED,
                "total_pairs": total,
                "host_pairs": host,
                "probe_pairs": probe_pairs if probe_pairs > 0 else 0,
                "probe_host_pairs": probe_host,
                "elapsed_s": round(time.monotonic() - started, 1),
            }
        )
    )
    reporter.ok(
        f"host filtering: {host:,} of {total:,} pairs ({host / max(total, 1):.2%}) were host; "
        f"{human_duration(time.monotonic() - started)}"
    )
    return HostFilterResult(
        nonhost_r1=out1, nonhost_r2=out2 if r2 else None,
        total_pairs=total, host_pairs=host, cached=False,
        probe_pairs=probe_pairs if probe_pairs > 0 else 0, probe_host_pairs=probe_host,
    )


def _parse_bowtie2_summary(stderr: str) -> tuple[int, int]:
    """Total pairs and pairs aligning (any mode) from Bowtie2's summary."""
    total = 0
    unaligned = 0
    for line in stderr.splitlines():
        text = line.strip()
        if text.endswith("reads; of these:"):
            total = int(text.split()[0])
        elif "aligned concordantly 0 times" in text and "(" in text and "of these" not in text:
            unaligned = int(text.split()[0])
    # Unaligned-concordant pairs are the ones written to --un-conc. Anything
    # else aligned in some fashion and is treated as host.
    return total, max(0, total - unaligned) if total else (0, 0)


def _fingerprint(r1: Path, r2: Path | None) -> str:
    digest = hashlib.sha256()
    for path in (r1, r2):
        if path is None:
            continue
        stat = path.stat()
        digest.update(f"{path.name}|{stat.st_size}|{int(stat.st_mtime)}".encode())
    return digest.hexdigest()[:16]


# --------------------------------------------------------------------------- #
# profiling
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class TaxonomicProfile:
    """One sample's MetaPhlAn output, parsed."""

    index: str
    profiler_version: str
    #: full clade lineage -> relative abundance percent (all levels)
    clades: dict[str, float]
    unknown_percent: float
    n_reads_processed: int
    elapsed_s: float
    cached: bool
    host: HostFilterResult | None = None
    command: tuple[str, ...] = field(default=())

    @property
    def family(self) -> str:
        return f"MetaPhlAn {self.profiler_version.split('.')[0]}"

    def level(self, prefix: str) -> dict[str, float]:
        """Abundances at one rank, keyed by bare name (e.g. ``s__`` -> species)."""
        out: dict[str, float] = {}
        for lineage, value in self.clades.items():
            last = lineage.rsplit("|", 1)[-1]
            if last.startswith(prefix) and "|t__" not in lineage:
                out[last[len(prefix):]] = value
        return out

    @property
    def species(self) -> dict[str, float]:
        return self.level("s__")

    @property
    def genera(self) -> dict[str, float]:
        return self.level("g__")

    @property
    def phyla(self) -> dict[str, float]:
        return self.level("p__")

    def to_json(self) -> dict[str, Any]:
        species = self.species
        return {
            "engine": "metaphlan",
            "profiler_version": self.profiler_version,
            "database": self.index,
            "unknown_percent": round(self.unknown_percent, 3),
            "n_reads_processed": self.n_reads_processed,
            "n_species_detected": sum(1 for v in species.values() if v > 0),
            "host_filter": None if self.host is None else self.host.to_json(),
            "cached": self.cached,
            "elapsed_s": round(self.elapsed_s, 1),
        }


def run_metaphlan(
    tools: MetaphlanTools,
    *,
    r1: Path,
    r2: Path | None,
    db_dir: Path,
    work_dir: Path,
    reporter: Reporter,
    threads: int = 8,
    index: str = PINNED_INDEX,
    host: HostFilterResult | None = None,
) -> TaxonomicProfile:
    """Profile one sample. Caches the Bowtie2 alignment so re-profiling is instant."""
    work_dir.mkdir(parents=True, exist_ok=True)
    bowtie2out = work_dir / f"metaphlan.{index}.bowtie2.bz2"
    profile_path = work_dir / f"metaphlan.{index}.tsv"
    stamp = work_dir / f"metaphlan.{index}.json"

    fingerprint = _fingerprint(r1, r2)
    if stamp.is_file() and profile_path.is_file():
        meta = json.loads(stamp.read_text())
        if meta.get("cache_version") == CACHE_VERSION and meta.get("inputs") == fingerprint:
            clades, unknown, n_reads = parse_profile(profile_path)
            return TaxonomicProfile(
                index=index, profiler_version=meta.get("version", tools.version),
                clades=clades, unknown_percent=unknown, n_reads_processed=n_reads,
                elapsed_s=meta.get("elapsed_s", 0.0), cached=True, host=host,
                command=tuple(meta.get("command", ())),
            )

    inputs = str(r1) if r2 is None else f"{r1},{r2}"
    command = [
        tools.metaphlan, inputs,
        "--input_type", "fastq",
        "--bowtie2db", str(db_dir),
        "--index", index,
        "--nproc", str(threads),
        "--unknown_estimation",
        "-o", str(profile_path),
    ]
    if bowtie2out.is_file():
        # Re-use the cached alignment; MetaPhlAn skips Bowtie2 entirely.
        command = [
            tools.metaphlan, str(bowtie2out), "--input_type", "bowtie2out",
            "--bowtie2db", str(db_dir), "--index", index, "--nproc", str(threads),
            "--unknown_estimation", "-o", str(profile_path),
        ]
    else:
        command += ["--bowtie2out", str(bowtie2out)]

    reporter.info(f"  MetaPhlAn {tools.version} against {index}")
    started = time.monotonic()
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, env=_env_for(tools)
    )
    elapsed = time.monotonic() - started
    if result.returncode != 0 or not profile_path.is_file():
        raise OpenBiotaError(
            f"MetaPhlAn failed (exit {result.returncode}):\n{result.stderr[-3000:]}"
        )
    clades, unknown, n_reads = parse_profile(profile_path)
    stamp.write_text(
        json.dumps(
            {
                "cache_version": CACHE_VERSION,
                "inputs": fingerprint,
                "version": tools.version,
                "index": index,
                "elapsed_s": round(elapsed, 1),
                "command": [Path(command[0]).name, *command[1:]],
            }
        )
    )
    reporter.ok(
        f"taxonomic profile: {sum(1 for k, v in clades.items() if '|s__' in k and '|t__' not in k and v > 0)} "
        f"species, {unknown:.1f}% unclassified, {human_duration(elapsed)}"
    )
    return TaxonomicProfile(
        index=index, profiler_version=tools.version, clades=clades,
        unknown_percent=unknown, n_reads_processed=n_reads, elapsed_s=elapsed,
        cached=False, host=host, command=tuple(command),
    )


def parse_profile(path: Path) -> tuple[dict[str, float], float, int]:
    """Parse a MetaPhlAn 3 profile into (clades, unknown_percent, n_reads)."""
    clades: dict[str, float] = {}
    unknown = 0.0
    n_reads = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "reads processed" in line:
                digits = "".join(ch for ch in line if ch.isdigit())
                n_reads = int(digits) if digits else 0
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name = parts[0]
        try:
            value = float(parts[2])
        except ValueError:
            continue
        if name == "UNKNOWN":
            unknown = value
            continue
        clades[name] = value
    return clades, unknown, n_reads
