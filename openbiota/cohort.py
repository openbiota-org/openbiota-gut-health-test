"""Reference cohort: turn public metagenomes into empirical expected ranges.

A single sample's "35 copies per 100 genomes" means nothing on its own. It
means something once you know where it falls in a population. This module
builds that population from public data and reduces it to percentiles.

The cohort used by default is ENA study ``PRJNA1271016`` — the 294 adult stool
shotgun metagenomes from the urdA source study. It is the best available match
for this pipeline's inputs: same body site, same library strategy, and the same
NovaSeq X platform family, which keeps read length and error profile
comparable.

Downloads are HTTP range requests for the first N bytes of each gzipped FASTQ,
decompressed as a truncated stream. Fetching a 600,000-read prefix costs about
36 MB per mate instead of the 2 GB full file. That is sound here because every
figure this tool reports is depth-independent by construction — the rpoB
denominator scales with the numerator — so a shallow prefix gives the same
expected value with wider Poisson error.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import gzip
import io
import json
import os
import random
import shutil
import statistics
import subprocess
import time
import urllib.error
import urllib.request
import zlib
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import InputError, OpenBiotaError
from openbiota.logging_util import Reporter, human_bytes, human_duration
from openbiota.net import USER_AGENT, get

ENA_PORTAL: Final = "https://www.ebi.ac.uk/ena/portal/api/filereport"

#: Default cohort: the urdA source study's own 294 stool metagenomes.
DEFAULT_STUDY: Final = "PRJNA1271016"

#: Compressed bytes needed per 1,000 read pairs, measured on this study
#: (8 MB of gzip yielded 132,511 records). Used to size the range request.
BYTES_PER_KILOREAD: Final = 63_000

#: Percentiles reported for each panel.
PERCENTILES: Final = (5, 10, 25, 50, 75, 90, 95)


@dataclass(frozen=True, slots=True)
class CohortRun:
    run: str
    sample: str
    r1_url: str
    r2_url: str
    read_count: int

    @property
    def urls(self) -> tuple[str, str]:
        return self.r1_url, self.r2_url


@dataclass(frozen=True, slots=True)
class PanelRange:
    """Empirical distribution of one panel across the cohort."""

    panel: str
    n: int
    percentiles: dict[int, float]
    mean: float
    stdev: float
    n_detected: int
    min_value: float
    max_value: float

    @property
    def median(self) -> float:
        return self.percentiles[50]

    @property
    def detection_rate(self) -> float:
        return self.n_detected / self.n if self.n else 0.0

    def percentile_of(self, value: float) -> float:
        """Where ``value`` falls in this distribution, 0-100.

        A value that ties with a run of equal stored quantiles gets the
        **midpoint** of that run, not its floor. This matters most on
        zero-inflated panels: methane has p5 = p25 = p50 = 0 because 57% of
        the cohort has none, and returning the floor put a zero at the 5th
        percentile - "below the usual range" for the commonest value there
        is. The midpoint is the same convention as the midrank used for
        the taxonomic and biofilm references, so the three agree.
        """
        points = sorted(self.percentiles.items())

        def tie_span(target: float) -> float | None:
            """Midpoint of the stored percentiles whose value equals target."""
            tied = [p for p, v in points if v == target]
            return (tied[0] + tied[-1]) / 2.0 if tied else None

        lowest_p, lowest_v = points[0]
        highest_p, highest_v = points[-1]
        if value <= lowest_v:
            # Strictly below the lowest stored quantile: it is at the floor.
            # Equal to it: share the tie block.
            return float(lowest_p) if value < lowest_v else (
                tie_span(value) or float(lowest_p)
            )
        if value >= highest_v:
            return float(highest_p) if value > highest_v else (
                tie_span(value) or float(highest_p)
            )
        for (p_lo, v_lo), (p_hi, v_hi) in zip(points, points[1:], strict=False):
            if v_lo <= value <= v_hi:
                if v_hi == v_lo:
                    return tie_span(value) or float(p_lo)
                frac = (value - v_lo) / (v_hi - v_lo)
                return p_lo + frac * (p_hi - p_lo)
        return 50.0

    def classify(self, value: float | None) -> str:
        """Plain-language position: below / within / above the usual range."""
        if value is None:
            return "not measurable"
        pct = self.percentile_of(value)
        if pct < 10:
            return "below the usual range"
        if pct < 25:
            return "lower end of usual"
        if pct <= 75:
            return "within the usual range"
        if pct <= 90:
            return "upper end of usual"
        return "above the usual range"

    def to_json(self) -> dict[str, Any]:
        return {
            "panel": self.panel,
            "cohort_n": self.n,
            "percentiles": {str(k): round(v, 4) for k, v in sorted(self.percentiles.items())},
            "mean": round(self.mean, 4),
            "stdev": round(self.stdev, 4),
            "detected_in": self.n_detected,
            "detection_rate": round(self.detection_rate, 4),
            "min": round(self.min_value, 4),
            "max": round(self.max_value, 4),
        }

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> PanelRange:
        return cls(
            panel=payload["panel"],
            n=int(payload["cohort_n"]),
            percentiles={int(k): float(v) for k, v in payload["percentiles"].items()},
            mean=float(payload["mean"]),
            stdev=float(payload["stdev"]),
            n_detected=int(payload["detected_in"]),
            min_value=float(payload["min"]),
            max_value=float(payload["max"]),
        )


@dataclass(frozen=True, slots=True)
class ReferenceRanges:
    """Percentile ranges for every panel, plus provenance."""

    study: str
    description: str
    n_samples: int
    reads_per_sample: int
    built_at: str
    panels: dict[str, PanelRange] = field(default_factory=dict)
    genes: dict[str, PanelRange] = field(default_factory=dict)
    #: How each cohort sample's reads were chosen. "prefix" means the head of
    #: the file, which is measurably biased; "random" means uniform across the
    #: whole file. Recorded because it changes every percentile in the report.
    sampling: str = "prefix"
    #: Seed for `sampling="random"`, so a cohort can be rebuilt exactly.
    seed: int | None = None
    #: ENA's title for the study, so the report can name the population.
    study_title: str = ""
    #: The run accessions that produced these ranges, in the order screened.
    #: With `sampling` and `seed`, this is everything needed to rebuild the
    #: cohort byte-for-byte; `openbiota cohort --runs` accepts it directly.
    runs: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "study": self.study,
            "description": self.description,
            "n_samples": self.n_samples,
            "reads_per_sample": self.reads_per_sample,
            "built_at": self.built_at,
            "percentiles_reported": list(PERCENTILES),
            "sampling": self.sampling,
            "seed": self.seed,
            "runs": list(self.runs),
            "study_title": self.study_title,
            "panels": {k: v.to_json() for k, v in sorted(self.panels.items())},
            "genes": {k: v.to_json() for k, v in sorted(self.genes.items())},
        }

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> ReferenceRanges:
        return cls(
            study=payload["study"],
            description=payload["description"],
            n_samples=int(payload["n_samples"]),
            reads_per_sample=int(payload["reads_per_sample"]),
            built_at=payload["built_at"],
            panels={k: PanelRange.from_json(v) for k, v in payload.get("panels", {}).items()},
            genes={k: PanelRange.from_json(v) for k, v in payload.get("genes", {}).items()},
            # A cohort built before sampling was recorded is a prefix cohort.
            sampling=str(payload.get("sampling", "prefix")),
            seed=payload.get("seed"),
            runs=tuple(payload.get("runs", ())),
            study_title=str(payload.get("study_title", "")),
        )

    @classmethod
    def load(cls, path: Path) -> ReferenceRanges:
        return cls.from_json(json.loads(path.read_text(encoding="utf-8")))

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- #
# ENA discovery and prefix download
# --------------------------------------------------------------------------- #


def study_title(study: str) -> str:
    """The study's own title from ENA, so a cohort can say who it is.

    "294 adult stool metagenomes" says nothing about the people. The title
    does — for `PRJNA1271016` it is "Gut microbiome and AD study in MARS
    cohort", which tells a reader the reference population is an older-adult
    Alzheimer's-disease cohort before they compare a child against it.
    """
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=study"
        f"&query=study_accession%3D{study}&fields=study_title,study_description&format=tsv"
    )
    try:
        rows = get(url, timeout=60).text().splitlines()
    except Exception:  # noqa: BLE001 - a missing title never blocks a build
        return ""
    if len(rows) < 2:
        return ""
    header = rows[0].split("\t")
    values = dict(zip(header, rows[1].split("\t"), strict=False))
    title = values.get("study_title", "").strip()
    description = values.get("study_description", "").strip()
    # The title is often a project codename; the description says what it is.
    if description and description.lower() != title.lower():
        return f"{title} — {description}" if title else description
    return title


def list_runs(study: str, reporter: Reporter) -> list[CohortRun]:
    """Query the ENA portal for every paired run in ``study``."""
    url = (
        f"{ENA_PORTAL}?accession={study}&result=read_run"
        "&fields=run_accession,sample_accession,fastq_ftp,read_count&format=tsv"
    )
    rows = get(url, timeout=180).text().splitlines()
    if len(rows) < 2:
        raise OpenBiotaError(f"ENA returned no runs for study {study!r}")

    runs: list[CohortRun] = []
    for row in rows[1:]:
        parts = row.split("\t")
        if len(parts) < 4:
            continue
        ftp = parts[2].split(";")
        if len(ftp) < 2 or not ftp[0] or not ftp[1]:
            continue
        try:
            count = int(parts[3])
        except ValueError:
            continue
        runs.append(
            CohortRun(
                run=parts[0],
                sample=parts[1],
                r1_url=f"https://{ftp[0]}",
                r2_url=f"https://{ftp[1]}",
                read_count=count,
            )
        )
    # ENA returns runs in no stable order — two calls minutes apart can differ.
    # Everything downstream slices this list, so an unsorted list makes the
    # cohort's membership depend on when it was built. Sort by accession.
    runs.sort(key=lambda r: r.run)
    reporter.info(f"  ENA study {study}: {len(runs)} paired runs available")
    return runs


#: Minimum read pairs a sample must yield to be worth screening.
#:
#: Shallow samples cost precision, not correctness: every reported figure is a
#: ratio to `rpoB`, and the denominator scales with the numerator. That was
#: asserted here for a long time; it has now been measured, and the measurement
#: is worth recording because it decided how this module samples.
#:
#: Screening one sample three ways at a 600k budget — the first 600k read
#: pairs, a uniform random 600k, and full depth at 8.03M — gives, against the
#: full-depth answer:
#:
#:     first 600k     median error 19.6%, 16 of 24 panels more than 10% off
#:     random 600k    median error 11.6%, 13 of 24 panels more than 10% off
#:
#: and the random-600k deviations sit inside Poisson counting error (median
#: |z| = 0.87 against an expectation of 1.0, no panel beyond 3σ, sign split
#: 16 high / 7 low). So:
#:
#: * **Position in the file is a real bias.** The head of a FASTQ is not a
#:   random sample of it — the first tiles of a flowcell carry a different
#:   quality profile, and at a fixed identity threshold that changes how many
#:   fragments each panel accepts. Sampling uniformly removes it.
#: * **Depth is only precision.** Once sampling is uniform there is no
#:   directional drift with depth, so the cohort does not need to be
#:   depth-matched to the query — it needs to be deep enough that its counting
#:   noise is small against the biological spread it is trying to describe.
#:
#: Hence `sampling="random"` is the default, and `--reads` is chosen so that
#: per-panel counting noise (~2.5% at 12M pairs) is well under the cohort's
#: biological coefficient of variation (10–25% per panel).
MIN_USABLE_PAIRS: Final = 120_000

#: Sampling strategies for cohort reads.
SAMPLING_MODES: Final = ("random", "prefix")

#: Where the reads come from. "auto" uses the SRA toolkit when it is installed
#: and falls back to ENA per run. "ena" refuses the toolkit outright, which is
#: what you want when NCBI's resolver is throttling: a blocked fastq-dump costs
#: a two-minute timeout per run before the fallback even starts, and a hundred
#: of those is three wasted hours. "sra" is the mirror image, for when ENA is
#: the sick one. "odp" hands fastq-dump the AWS Open Data object directly.
TRANSPORTS: Final = ("auto", "ena", "sra", "odp")

#: The AWS Open Data mirror of the SRA, addressed as a plain object.
#:
#: What fails under load is NCBI's *name resolver*, not its data: the toolkit
#: asks trace.ncbi.nlm.nih.gov where a run lives before reading a byte of it,
#: and that service is the one that starts refusing. Handing fastq-dump this
#: URL skips the question entirely. Spot-range access still works — the SRA
#: format seeks over HTTP range requests — so `-N`/`-X` windows cost what they
#: always did, and random sampling survives intact. Measured here at 31 MB/s
#: against ENA's 1.1, and it scales with concurrency where ENA does not.
ODP_URL: Final = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"


def _dump_target(run: str, transport: str) -> str:
    """What to hand fastq-dump: a bare accession, or the Open Data object."""
    return ODP_URL.format(run=run) if transport == "odp" else run

#: Default seed for cohort subsampling. Fixed so a rebuild reproduces the
#: same cohort exactly; change it only to measure sampling variance.
DEFAULT_SEED: Final = 20260909

#: Bytes pulled per socket read.
_CHUNK: Final = 1 << 20

#: SRA toolkit config, written into the work directory so `fastq-dump` never
#: drops into its interactive first-run prompt.
_VDB_CONFIG: Final = (
    '/LIBS/GUID = "00000000-0000-0000-0000-000000000000"\n'
    '/libs/cloud/report_instance_identity = "false"\n'
)


def sra_toolkit_available() -> bool:
    return shutil.which("fastq-dump") is not None


def fetch_with_fastq_dump(
    run: str,
    dest_dir: Path,
    read_pairs: int,
    *,
    config_dir: Path,
    timeout: float = 1800.0,
    transport: str = "auto",
) -> tuple[Path, Path, int]:
    """Fetch the first ``read_pairs`` spots of an SRA run with `fastq-dump -X`.

    Preferred over range-requesting the ENA gzip: it is the purpose-built path,
    it returns both mates consistently, and `--origfmt` preserves the original
    instrument read ids so mate collapse is exercised exactly as it is on real
    input.
    """
    r1 = dest_dir / f"{run}_1.fastq"
    r2 = dest_dir / f"{run}_2.fastq"
    marker = dest_dir / f"{run}.done"
    if r1.is_file() and r2.is_file() and marker.is_file():
        try:
            return r1, r2, int(marker.read_text(encoding="utf-8").strip())
        except ValueError:
            pass

    exe = shutil.which("fastq-dump")
    if exe is None:
        raise OpenBiotaError("fastq-dump not found on PATH (brew install sratoolkit)")

    dest_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "vdb.kfg"
    if not config_path.is_file():
        config_path.write_text(_VDB_CONFIG, encoding="utf-8")

    cmd = [
        exe,
        "-X", str(read_pairs),
        "--split-files",
        "--skip-technical",
        # keep the instrument read id, identical across mates and with no
        # /1 /2 suffix, matching the screened data
        "--origfmt",
        "--outdir", str(dest_dir),
        _dump_target(run, transport),
    ]
    try:
        subprocess.run(  # noqa: S603 — fixed argv, no shell
            cmd,
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout,
            env={**os.environ, "NCBI_SETTINGS": str(config_path), "HOME": str(config_dir)},
        )
    except subprocess.CalledProcessError as exc:
        raise OpenBiotaError(
            f"{run}: fastq-dump failed — {(exc.stderr or '').strip()[-300:]}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise OpenBiotaError(f"{run}: fastq-dump timed out after {timeout:.0f}s") from exc

    if not (r1.is_file() and r2.is_file()):
        raise OpenBiotaError(f"{run}: fastq-dump produced no paired output")

    pairs = _count_records(r1)
    other = _count_records(r2)
    if pairs != other:
        pairs = min(pairs, other)
        _truncate_to(r1, pairs)
        _truncate_to(r2, pairs)
    if pairs < MIN_USABLE_PAIRS:
        raise OpenBiotaError(f"{run}: only {pairs:,} read pairs retrieved")

    marker.write_text(str(pairs), encoding="utf-8")
    return r1, r2, pairs


def _count_records(path: Path) -> int:
    lines = 0
    with path.open("rb") as fh:
        while chunk := fh.read(_CHUNK):
            lines += chunk.count(b"\n")
    return lines // 4


def _read_prefix_bytes(url: str, want_bytes: int, timeout: float) -> bytes:
    """Range-request ``want_bytes``, keeping whatever actually arrives.

    ENA frequently drops a range transfer partway through. That is not an
    error for this purpose — the stream is being truncated deliberately — so
    a short read is accepted rather than retried into the ground.
    """
    buffer = bytearray()
    # A dropped transfer is resumed from the byte it stopped at: the gzip
    # stream is byte-addressable, so appending the next range continues it
    # exactly. Bounded so a server that keeps dropping cannot loop forever.
    for _attempt in range(8):
        start = len(buffer)
        if start >= want_bytes:
            break
        request = urllib.request.Request(
            url, headers={"Range": f"bytes={start}-{want_bytes - 1}", "User-Agent": USER_AGENT}
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if start and response.status != 206:
                    break  # server ignored the range; cannot resume safely
                while len(buffer) < want_bytes:
                    try:
                        chunk = response.read(_CHUNK)
                    except Exception:  # noqa: BLE001 — a dropped transfer keeps its prefix
                        break
                    if not chunk:
                        break
                    buffer.extend(chunk)
        except (urllib.error.URLError, TimeoutError, OSError):
            if not buffer:
                raise
        if len(buffer) == start:
            break  # no progress; stop rather than hammer the mirror
        if len(buffer) < want_bytes:
            time.sleep(1.5)
    return bytes(buffer)


class _ResumableStream(io.RawIOBase):
    """A read-only byte stream over a remote file that survives dropped transfers.

    ENA drops long transfers routinely — reliably enough that the prefix
    fetcher above is written around it. A drop is recoverable because the
    object is a plain byte range: reconnect with a `Range` header at the byte
    already consumed and carry on. The consumer sees one continuous stream and
    never learns the socket was replaced, which is what lets a gzip
    decompressor run over the whole of a multi-gigabyte file.
    """

    def __init__(self, url: str, *, timeout: float = 120.0, attempts: int = 40) -> None:
        super().__init__()
        self._url = url
        self._timeout = timeout
        self._attempts = attempts
        self._pos = 0
        self._response: Any | None = None
        self._exhausted = False
        self._reconnects = 0

    @property
    def reconnects(self) -> int:
        return self._reconnects

    def readable(self) -> bool:
        return True

    def _connect(self) -> None:
        request = urllib.request.Request(
            self._url,
            headers={"Range": f"bytes={self._pos}-", "User-Agent": USER_AGENT},
        )
        response = urllib.request.urlopen(request, timeout=self._timeout)  # noqa: S310 - ENA
        if self._pos and response.status != 206:
            # The server ignored the range. Continuing would silently splice
            # the file's beginning into the middle of the stream.
            response.close()
            raise OpenBiotaError(f"{self._url}: server will not resume from byte {self._pos}")
        self._response = response

    def readinto(self, buffer: Any) -> int:  # noqa: ANN401 - buffer protocol
        if self._exhausted:
            return 0
        last: Exception | None = None
        for _ in range(self._attempts):
            if self._response is None:
                try:
                    self._connect()
                except (urllib.error.URLError, TimeoutError, OSError, OpenBiotaError) as exc:
                    last = exc
                    time.sleep(2.0)
                    continue
            try:
                chunk = self._response.read(len(buffer))
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last = exc
                self._drop()
                continue
            if not chunk:
                # A clean empty read is end of file. A transfer the server
                # cut short raises instead, and is retried above; a truncated
                # gzip tail is caught by the decompressor downstream.
                self._exhausted = True
                return 0
            buffer[: len(chunk)] = chunk
            self._pos += len(chunk)
            return len(chunk)
        raise OpenBiotaError(f"{self._url}: gave up at byte {self._pos} ({last})")

    def _drop(self) -> None:
        if self._response is not None:
            with contextlib.suppress(Exception):
                self._response.close()
        self._response = None
        self._reconnects += 1
        time.sleep(1.0)

    def close(self) -> None:
        self._drop()
        super().close()


def sample_remote_fastq(
    url: str,
    dest: Path,
    *,
    probability: float,
    seed: int,
    limit: int | None = None,
) -> int:
    """Write a uniform random subsample of a remote gzipped FASTQ.

    Every record is kept with probability `probability`, decided by a
    generator seeded with `seed`. Two mates given the same seed keep the same
    record indices, because both files hold their records in the same order
    and one draw is consumed per record — so pairing survives without either
    pass knowing about the other.

    The whole remote file is streamed and decompressed; only the retained
    records are written, so the disk cost is the sample rather than the
    source. Returns the number of records written.
    """
    if not 0.0 < probability <= 1.0:
        raise OpenBiotaError(f"sampling probability out of range: {probability}")

    marker = dest.with_suffix(dest.suffix + ".done")
    if dest.is_file() and marker.is_file():
        with contextlib.suppress(ValueError):
            return int(marker.read_text(encoding="utf-8").strip())

    rng = random.Random(seed)
    kept = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    stream = _ResumableStream(url)
    try:
        with (
            gzip.open(io.BufferedReader(stream, buffer_size=1 << 22), "rb") as source,
            partial.open("wb") as out,
        ):
            while True:
                header = source.readline()
                if not header:
                    break
                record = header + source.readline() + source.readline() + source.readline()
                if rng.random() < probability:
                    out.write(record)
                    kept += 1
                    if limit is not None and kept >= limit:
                        break
    except EOFError:
        # A truncated gzip tail. Everything already written is intact, and a
        # short sample is reported rather than thrown away.
        pass
    finally:
        stream.close()

    partial.replace(dest)
    marker.write_text(str(kept), encoding="utf-8")
    return kept


#: How many scattered windows make up one random sample. A window costs a
#: fixed ~6.4 s of `fastq-dump` start-up plus ~44,000 spots/second, so window
#: count trades run time against how widely the sample is spread across the
#: flowcell. Thirty windows over a 32M-spot run covers roughly 450 of its
#: ~780 tiles and costs about four minutes per run.
SAMPLE_WINDOWS: Final = 30


def _window_starts(run: CohortRun, target_pairs: int, windows: int, seed: int) -> list[tuple[int, int]]:
    """Evenly spaced spot ranges covering the whole run, with a seeded phase.

    Systematic sampling with a random offset: the windows are spread across
    the entire run rather than taken from its head, which is the property the
    measurement in `MIN_USABLE_PAIRS` showed actually matters. The phase is
    derived from the run accession so two runs do not sample the same
    positions, and from `seed` so a rebuild reproduces them exactly.
    """
    total = max(1, run.read_count)
    if target_pairs >= total:
        return [(1, total)]
    per_window = max(1, target_pairs // windows)
    stride = total // windows
    if stride <= per_window:  # the sample is most of the run; just take a prefix
        return [(1, min(total, target_pairs))]
    rng = random.Random(seed ^ zlib.crc32(run.run.encode()))
    phase = rng.randrange(stride - per_window)
    ranges: list[tuple[int, int]] = []
    for index in range(windows):
        first = 1 + index * stride + phase
        last = min(total, first + per_window - 1)
        if first <= total:
            ranges.append((first, last))
    return ranges


def fetch_random_sample(
    run: CohortRun,
    dest_dir: Path,
    target_pairs: int,
    *,
    seed: int,
    windows: int = SAMPLE_WINDOWS,
    config_dir: Path | None = None,
    timeout: float = 2400.0,
    transport: str = "auto",
) -> tuple[Path, Path, int]:
    """Sample `target_pairs` read pairs spread across the whole of a run.

    Uses `fastq-dump`'s spot-range access, which seeks into the SRA archive
    rather than transferring the file: a window at spot 18,000,000 costs the
    same as one at spot 1. That is what makes an unbiased sample affordable —
    streaming the whole of a 32M-spot run to keep 3M of it takes fifteen
    minutes a run, while thirty scattered windows take four.

    Falls back to streaming the ENA gzip when the SRA toolkit is missing,
    which is correct but much slower.
    """
    r1 = dest_dir / f"{run.run}_1.fastq"
    r2 = dest_dir / f"{run.run}_2.fastq"
    marker = dest_dir / f"{run.run}.done"
    if r1.is_file() and r2.is_file() and marker.is_file():
        with contextlib.suppress(ValueError):
            return r1, r2, int(marker.read_text(encoding="utf-8").strip())

    if not sra_toolkit_available():
        total = max(1, run.read_count)
        probability = min(1.0, target_pairs / total)
        got1 = sample_remote_fastq(run.r1_url, r1, probability=probability, seed=seed)
        got2 = sample_remote_fastq(run.r2_url, r2, probability=probability, seed=seed)
        pairs = min(got1, got2)
        if got1 != got2:
            _truncate_to(r1, pairs)
            _truncate_to(r2, pairs)
        if pairs < MIN_USABLE_PAIRS:
            raise OpenBiotaError(f"{run.run}: only {pairs:,} read pairs sampled")
        marker.write_text(str(pairs), encoding="utf-8")
        return r1, r2, pairs

    exe = shutil.which("fastq-dump")
    if exe is None:  # pragma: no cover - guarded by sra_toolkit_available
        raise OpenBiotaError("fastq-dump disappeared between check and use")

    dest_dir.mkdir(parents=True, exist_ok=True)
    config_dir = config_dir or dest_dir.parent / "sra"
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "vdb.kfg"
    if not config_path.is_file():
        config_path.write_text(_VDB_CONFIG, encoding="utf-8")
    env = {**os.environ, "NCBI_SETTINGS": str(config_path), "HOME": str(config_dir)}

    scratch = dest_dir / f".{run.run}.windows"
    if scratch.exists():
        shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True)
    pairs = 0
    try:
        with r1.open("wb") as out1, r2.open("wb") as out2:
            for index, (first, last) in enumerate(
                _window_starts(run, target_pairs, windows, seed)
            ):
                window_dir = scratch / str(index)
                cmd = [
                    exe,
                    "-N", str(first),
                    "-X", str(last),
                    "--split-files",
                    "--skip-technical",
                    # keep the instrument read id, identical across mates and
                    # with no /1 /2 suffix, matching the screened data
                    "--origfmt",
                    "--outdir", str(window_dir),
                    _dump_target(run.run, transport),
                ]
                try:
                    subprocess.run(  # noqa: S603 - fixed argv, no shell
                        cmd, check=True, capture_output=True, text=True,
                        timeout=timeout, env=env,
                    )
                except subprocess.CalledProcessError as exc:
                    raise OpenBiotaError(
                        f"{run.run}: fastq-dump failed on spots {first}-{last} — "
                        f"{(exc.stderr or '').strip()[-200:]}"
                    ) from exc
                except subprocess.TimeoutExpired as exc:
                    raise OpenBiotaError(
                        f"{run.run}: fastq-dump timed out on spots {first}-{last}"
                    ) from exc
                w1 = window_dir / f"{run.run}_1.fastq"
                w2 = window_dir / f"{run.run}_2.fastq"
                if not (w1.is_file() and w2.is_file()):
                    raise OpenBiotaError(f"{run.run}: window {first}-{last} produced no pair")
                got = min(_count_records(w1), _count_records(w2))
                # Append only whole pairs, so the concatenation stays in sync.
                for source, sink in ((w1, out1), (w2, out2)):
                    with source.open("rb") as fh:
                        for line_no, line in enumerate(fh):
                            if line_no >= got * 4:
                                break
                            sink.write(line)
                pairs += got
                shutil.rmtree(window_dir, ignore_errors=True)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if pairs < MIN_USABLE_PAIRS:
        raise OpenBiotaError(f"{run.run}: only {pairs:,} read pairs sampled")
    marker.write_text(str(pairs), encoding="utf-8")
    return r1, r2, pairs


def fetch_prefix(
    url: str,
    dest: Path,
    read_pairs: int,
    *,
    attempts: int = 4,
    min_pairs: int = MIN_USABLE_PAIRS,
) -> int:
    """Download the gzip prefix of ``url`` and write up to ``read_pairs`` records.

    A gzip stream decodes from the start, so a range request for the first N
    bytes yields a usable — if truncated — FASTQ. Output is trimmed to the last
    complete 4-line record.
    """
    marker = dest.with_suffix(dest.suffix + ".done")
    if dest.is_file() and marker.is_file():
        try:
            return int(marker.read_text(encoding="utf-8").strip())
        except ValueError:
            pass

    want_bytes = int(read_pairs / 1000 * BYTES_PER_KILOREAD * 1.3)
    best_lines: list[str] = []
    last_error = "no attempt"

    for attempt in range(1, attempts + 1):
        try:
            raw = _read_prefix_bytes(url, want_bytes, timeout=300.0)
        except Exception as exc:  # noqa: BLE001 — any transport failure is retryable
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < attempts:
                time.sleep(1.5 * attempt)
            continue

        if not raw:
            last_error = "server returned no data"
            if attempt < attempts:
                time.sleep(1.5 * attempt)
            continue

        decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            text = decompressor.decompress(raw).decode("utf-8", "replace")
        except zlib.error as exc:
            last_error = f"gzip prefix undecodable: {exc}"
            if attempt < attempts:
                time.sleep(1.5 * attempt)
            continue

        lines = text.split("\n")
        if len(lines) > len(best_lines):
            best_lines = lines
        if (len(best_lines) - 1) // 4 >= min_pairs:
            break
        last_error = f"prefix yielded only {(len(best_lines) - 1) // 4:,} pairs"
        if attempt < attempts:
            time.sleep(1.5 * attempt)

    complete = min(len(best_lines) - 1, read_pairs * 4) // 4 * 4
    pairs = complete // 4
    if pairs < min_pairs:
        raise OpenBiotaError(
            f"{url}: could only retrieve {pairs:,} read pairs, below the {min_pairs:,} "
            f"minimum — {last_error}"
        )

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(best_lines[:complete]) + "\n", encoding="utf-8")
    marker.write_text(str(pairs), encoding="utf-8")
    return pairs


def _truncate_to(path: Path, pairs: int) -> None:
    """Trim a FASTQ in place to its first ``pairs`` records."""
    lines = path.read_text(encoding="utf-8").split("\n")
    keep = pairs * 4
    if len(lines) - 1 <= keep:
        return
    path.write_text("\n".join(lines[:keep]) + "\n", encoding="utf-8")
    path.with_suffix(path.suffix + ".done").write_text(str(pairs), encoding="utf-8")


def download_cohort(
    runs: Sequence[CohortRun],
    *,
    dest_dir: Path,
    read_pairs: int,
    workers: int,
    reporter: Reporter,
    sampling: str = "random",
    seed: int = DEFAULT_SEED,
    transport: str = "auto",
) -> list[tuple[CohortRun, Path, Path]]:
    """Fetch reads for every run, in parallel.

    With `sampling="random"` (the default) each run's whole FASTQ pair is
    streamed and a uniform random `read_pairs` is retained; see
    `MIN_USABLE_PAIRS` for the measurement that made this the default. With
    `sampling="prefix"` only the head of each file is fetched, which is much
    cheaper and measurably biased.
    """
    if sampling not in SAMPLING_MODES:
        raise OpenBiotaError(f"unknown sampling mode {sampling!r}; expected one of {SAMPLING_MODES}")
    if transport not in TRANSPORTS:
        raise OpenBiotaError(f"unknown transport {transport!r}; expected one of {TRANSPORTS}")
    if sampling == "random" and transport == "ena":
        # Random sampling reads the whole file to scatter its picks. ENA serves
        # at about a megabyte a second and the files are gigabytes, so this
        # combination is half an hour per mate; say so now rather than after.
        raise OpenBiotaError(
            "random sampling streams the entire FASTQ, which over ENA is tens of minutes "
            "per run; use --sampling prefix with --transport ena, or --transport auto"
        )
    dest_dir.mkdir(parents=True, exist_ok=True)
    results: list[tuple[CohortRun, Path, Path]] = []
    failures: list[str] = []
    started = time.monotonic()
    done = 0

    use_sra = sampling == "prefix" and transport != "ena" and sra_toolkit_available()
    if sampling == "random":
        reporter.record(
            f"    sampling: uniform random {read_pairs:,} read pairs per run, seed {seed}, "
            "spot windows across the whole run"
        )
        reporter.record(
            f"    transport: {'AWS Open Data (no NCBI resolver)' if transport == 'odp' else transport}"
        )
    else:
        reporter.record(
            f"    transport: {'fastq-dump (SRA toolkit)' if use_sra else 'ENA gzip range requests'}"
        )

    def job(run: CohortRun) -> tuple[CohortRun, Path, Path]:
        if sampling == "random":
            r1, r2, _ = fetch_random_sample(
                run, dest_dir, read_pairs, seed=seed, transport=transport
            )
            return run, r1, r2

        if use_sra:
            try:
                r1, r2, _ = fetch_with_fastq_dump(
                    run.run, dest_dir, read_pairs,
                    config_dir=dest_dir.parent / "sra", transport=transport,
                )
                return run, r1, r2
            except OpenBiotaError as exc:
                # NCBI's resolver is down more often than ENA's mirror; fall
                # through to ENA for this run rather than lose the sample.
                reporter.record(f"    {run.run}: fastq-dump unavailable ({exc}); using ENA")

        # Fallback: range-request the ENA gzip prefixes. ENA truncates these
        # transfers unpredictably, so the two mates can land at different
        # depths; both files are in the same order, so trimming to the shorter
        # one restores a clean pair.
        r1 = dest_dir / f"{run.run}_1.fastq"
        r2 = dest_dir / f"{run.run}_2.fastq"
        got1 = fetch_prefix(run.r1_url, r1, read_pairs)
        got2 = fetch_prefix(run.r2_url, r2, read_pairs)
        if got1 != got2:
            _truncate_to(r1, min(got1, got2))
            _truncate_to(r2, min(got1, got2))
        return run, r1, r2

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(job, run): run for run in runs}
        for future in concurrent.futures.as_completed(futures):
            run = futures[future]
            done += 1
            try:
                results.append(future.result())
            except OpenBiotaError as exc:
                failures.append(f"{run.run}: {exc}")
                reporter.warn(f"  {run.run} skipped — {exc}")
                continue
            if done % 5 == 0 or done == len(runs):
                downloaded = sum(
                    p.stat().st_size for _, a, b in results for p in (a, b) if p.is_file()
                )
                reporter.info(
                    f"  cohort download {done}/{len(runs)} runs, "
                    f"{human_bytes(downloaded)} written, "
                    f"{human_duration(time.monotonic() - started)} elapsed"
                )

    if not results:
        raise InputError(
            "no cohort samples could be downloaded. Check network access to "
            "ftp.sra.ebi.ac.uk. Failures:\n  " + "\n  ".join(failures[:5])
        )
    if failures:
        reporter.warn(f"  {len(failures)} of {len(runs)} runs failed and were skipped")
    return results


# --------------------------------------------------------------------------- #
# range construction
# --------------------------------------------------------------------------- #


def _percentiles(values: Sequence[float]) -> dict[int, float]:
    if not values:
        return dict.fromkeys(PERCENTILES, 0.0)
    ordered = sorted(values)
    out: dict[int, float] = {}
    for p in PERCENTILES:
        # linear interpolation between order statistics
        position = (len(ordered) - 1) * (p / 100.0)
        lower = int(position)
        upper = min(lower + 1, len(ordered) - 1)
        weight = position - lower
        out[p] = ordered[lower] * (1 - weight) + ordered[upper] * weight
    return out


def build_ranges(
    *,
    study: str,
    description: str,
    reads_per_sample: int,
    panel_values: dict[str, list[float]],
    gene_values: dict[str, list[float]],
    sampling: str = "random",
    seed: int | None = DEFAULT_SEED,
    runs: Sequence[str] = (),
    study_title: str = "",
) -> ReferenceRanges:
    """Reduce per-sample values to percentile ranges."""

    def summarise(name: str, values: list[float]) -> PanelRange:
        return PanelRange(
            panel=name,
            n=len(values),
            percentiles=_percentiles(values),
            mean=statistics.fmean(values) if values else 0.0,
            stdev=statistics.stdev(values) if len(values) > 1 else 0.0,
            n_detected=sum(1 for v in values if v > 0),
            min_value=min(values) if values else 0.0,
            max_value=max(values) if values else 0.0,
        )

    n = max((len(v) for v in panel_values.values()), default=0)
    return ReferenceRanges(
        study=study,
        description=description,
        n_samples=n,
        reads_per_sample=reads_per_sample,
        built_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        panels={k: summarise(k, v) for k, v in panel_values.items()},
        genes={k: summarise(k, v) for k, v in gene_values.items()},
        sampling=sampling,
        seed=seed,
        runs=tuple(runs),
        study_title=study_title,
    )
