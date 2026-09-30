"""Input validation and read statistics — an extension, not spec.

The build spec says not to re-derive the documented properties of the input but
to validate them at runtime. That is what this module does, and it reports a
few extra sequencing metrics while it is already streaming the data.

Cost control: the exact record count needs a full pass, which is done with
chunked byte counting (no per-line Python work). Everything more expensive —
lengths, GC, quality, duplicate rate — is computed on a bounded prefix.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import json
import shutil
import statistics
import subprocess
import time
from collections import Counter
from collections.abc import Iterator
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import IO, Any, Final

from openbiota.errors import InputError
from openbiota.logging_util import Reporter, human_bytes, human_duration

CHUNK: Final = 1 << 23  # 8 MiB

#: Reads examined per mate for the sampled metrics.
DEFAULT_SAMPLE_READS: Final = 400_000

#: Phred+33 offset. Phred+64 encodings have not been produced by Illumina
#: since the HiSeq era and are detected rather than supported.
PHRED33_OFFSET: Final = 33


#: External decompressor, preferred over Python's zlib for whole-file passes.
#:
#: Inputs stay gzipped on disk — they are 600 MB compressed against 2.6 GB
#: plain, and every consumer in the pipeline (DIAMOND, Bowtie2, MetaPhlAn)
#: reads gzip natively. Nothing is ever decompressed to a file.
#:
#: Where a full pass *is* needed, streaming through the system `gzip -dc` beats
#: Python's `gzip` module by about 5x on this hardware — 3.9 s against 19.8 s
#: for a 594 MB mate — because zlib in Python pays interpreter overhead per
#: chunk. `pigz` was measured too and is *slower* here (7.2 s): a single gzip
#: stream cannot be inflated in parallel, so its extra threads only add
#: coordination cost.
GZIP_STREAMER: Final = "gzip"


def _stream_decompressed(path: Path) -> subprocess.Popen[bytes] | None:
    """Open a streaming decompressor for ``path``, or None to read directly."""
    if path.suffix.lower() != ".gz" or shutil.which(GZIP_STREAMER) is None:
        return None
    return subprocess.Popen(
        [GZIP_STREAMER, "-dc", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        bufsize=CHUNK,
    )


@contextlib.contextmanager
def open_bytes(path: Path, *, partial_ok: bool = False) -> Iterator[IO[bytes]]:
    """Binary stream over a FASTQ, transparently decompressing in a pipe.

    Set ``partial_ok`` when the caller intends to stop reading before the end
    of the file — taking a subsample, for instance. Closing the pipe early
    sends SIGPIPE to the decompressor, which then exits non-zero through no
    fault of the data, and treating that as corruption would reject a
    perfectly good input.

    Left at the default, a non-zero exit *is* reported, because a stream that
    was read to EOF and still failed means the file is genuinely damaged.
    """
    process = _stream_decompressed(path)
    if process is None:
        with path.open("rb") as handle:
            yield handle
        return
    assert process.stdout is not None
    try:
        yield process.stdout
    finally:
        process.stdout.close()
        status = process.wait()
        if not partial_ok and status not in (0, None):
            raise InputError(
                f"{path}: {GZIP_STREAMER} -dc exited with status {status} after reading the "
                "whole stream — the file is truncated or corrupt"
            )


def _open_text(path: Path, *, prefer_plain: bool = False) -> IO[str]:
    """Text stream over a FASTQ.

    ``prefer_plain`` is retained for call compatibility and ignored: inputs
    are read in whichever form they are on disk, and no plain-text twin is
    ever produced.
    """
    del prefer_plain
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def count_records(path: Path, reporter: Reporter | None = None) -> int:
    """Exact FASTQ record count via chunked newline counting."""
    lines = 0
    with open_bytes(path) as handle:
        while chunk := handle.read(CHUNK):
            lines += chunk.count(b"\n")
    if lines % 4 != 0:
        raise InputError(
            f"{path}: {lines} lines is not a multiple of 4 — the FASTQ is truncated or malformed"
        )
    if reporter is not None:
        reporter.record(f"    {path.name}: {lines // 4:,} records")
    return lines // 4


@dataclass(frozen=True, slots=True)
class MateQC:
    label: str
    path: str
    file_bytes: int
    gzipped: bool
    records: int | None
    sampled_reads: int
    min_length: int
    max_length: int
    mean_length: float
    total_bases_sampled: int
    gc_percent: float
    mean_quality: float
    quality_alphabet: str
    n_percent: float
    duplicate_fraction: float
    first_read_id: str
    last_sampled_read_id: str
    estimated_bases: int | None

    def to_json(self) -> dict[str, Any]:
        return {
            "mate": self.label,
            "path": self.path,
            "file_bytes": self.file_bytes,
            "gzipped": self.gzipped,
            "records": self.records,
            "sampled_reads": self.sampled_reads,
            "read_length_min": self.min_length,
            "read_length_max": self.max_length,
            "read_length_mean": round(self.mean_length, 1),
            "gc_percent": round(self.gc_percent, 2),
            "mean_phred": round(self.mean_quality, 1),
            "quality_alphabet": self.quality_alphabet,
            "n_base_percent": round(self.n_percent, 4),
            "sampled_duplicate_fraction": round(self.duplicate_fraction, 4),
            "first_read_id": self.first_read_id,
            "estimated_total_bases": self.estimated_bases,
        }


@dataclass(frozen=True, slots=True)
class SampleQC:
    mates: tuple[MateQC, ...]
    paired: bool
    ids_identical_between_mates: bool | None
    mate_suffixes_present: bool
    record_counts_match: bool | None
    notes: tuple[str, ...] = field(default=())
    cached: bool = False

    @property
    def read_pairs(self) -> int | None:
        counts = [m.records for m in self.mates if m.records is not None]
        if not counts:
            return None
        return min(counts)

    @property
    def failed(self) -> bool:
        """Whether input validation found a defect that invalidates depth.

        Only structural defects count. A mate-count mismatch or
        non-corresponding read IDs mean the pairing is broken, so any
        depth-dependent claim downstream is unsound. ``None`` means the
        check was not run, which is not a failure - callers must not read
        "not assessed" as "failed", or a fast run would silently suppress
        panels.

        This property exists because callers were reaching for
        ``getattr(qc, "failed", False)`` on a class that had no such
        attribute, so the default fired every time and every sample was
        treated as having passed.
        """
        return self.record_counts_match is False or (
            self.ids_identical_between_mates is False
        )

    @property
    def failure_reasons(self) -> tuple[str, ...]:
        """Why :attr:`failed` is true, for the report and the log."""
        reasons: list[str] = []
        if self.record_counts_match is False:
            reasons.append("mate record counts differ")
        if self.ids_identical_between_mates is False:
            reasons.append("read IDs do not correspond between mates")
        return tuple(reasons)

    @property
    def total_bases(self) -> int | None:
        values = [m.estimated_bases for m in self.mates if m.estimated_bases is not None]
        return sum(values) if values else None

    def to_json(self) -> dict[str, Any]:
        return {
            "paired": self.paired,
            "read_pairs": self.read_pairs,
            "estimated_total_bases": self.total_bases,
            "read_ids_identical_between_mates": self.ids_identical_between_mates,
            "mate_suffixes_present": self.mate_suffixes_present,
            "record_counts_match": self.record_counts_match,
            "notes": list(self.notes),
            "mates": [m.to_json() for m in self.mates],
        }


def profile_mate(
    path: Path,
    label: str,
    *,
    sample_reads: int = DEFAULT_SAMPLE_READS,
    count_all: bool = True,
    reporter: Reporter | None = None,
) -> tuple[MateQC, list[str]]:
    """Sampled metrics for one mate file, plus the sampled read ids."""
    lengths: list[int] = []
    gc = n_bases = 0
    quality_sum = 0
    quality_chars: Counter[str] = Counter()
    seen: set[bytes] = set()
    duplicates = 0
    ids: list[str] = []

    with _open_text(path, prefer_plain=True) as fh:
        while len(lengths) < sample_reads:
            header = fh.readline()
            if not header:
                break
            sequence = fh.readline().rstrip("\n")
            plus = fh.readline()
            quality = fh.readline().rstrip("\n")
            if not plus or not quality:
                raise InputError(
                    f"{path}: FASTQ record {len(lengths) + 1} is truncated "
                    "(missing '+' or quality line)"
                )
            if not header.startswith("@"):
                raise InputError(
                    f"{path}: record {len(lengths) + 1} header does not start with '@': "
                    f"{header[:60]!r}"
                )
            if len(sequence) != len(quality):
                raise InputError(
                    f"{path}: record {len(lengths) + 1} has {len(sequence)} bases but "
                    f"{len(quality)} quality values"
                )
            ids.append(header[1:].split()[0] if len(header) > 1 else "")
            lengths.append(len(sequence))
            upper = sequence.upper()
            gc += upper.count("G") + upper.count("C")
            n_bases += upper.count("N")
            quality_chars.update(quality)
            quality_sum += sum(ord(c) for c in quality) - PHRED33_OFFSET * len(quality)
            digest = hashlib.blake2b(upper.encode("ascii", "replace"), digest_size=12).digest()
            if digest in seen:
                duplicates += 1
            else:
                seen.add(digest)

    if not lengths:
        raise InputError(f"{path}: no FASTQ records found")

    total_bases = sum(lengths)
    records = count_records(path, reporter) if count_all else None
    mean_length = statistics.fmean(lengths)
    estimated = int(records * mean_length) if records is not None else None

    return (
        MateQC(
            label=label,
            path=str(path),
            file_bytes=path.stat().st_size,
            gzipped=path.suffix.lower() == ".gz",
            records=records,
            sampled_reads=len(lengths),
            min_length=min(lengths),
            max_length=max(lengths),
            mean_length=mean_length,
            total_bases_sampled=total_bases,
            gc_percent=gc / total_bases * 100.0 if total_bases else 0.0,
            mean_quality=quality_sum / total_bases if total_bases else 0.0,
            quality_alphabet="".join(sorted(quality_chars)),
            n_percent=n_bases / total_bases * 100.0 if total_bases else 0.0,
            duplicate_fraction=duplicates / len(lengths),
            first_read_id=ids[0] if ids else "",
            last_sampled_read_id=ids[-1] if ids else "",
            estimated_bases=estimated,
        ),
        ids,
    )


def profile_sample(
    mates: list[tuple[str, Path]],
    *,
    sample_reads: int = DEFAULT_SAMPLE_READS,
    count_all: bool = True,
    reporter: Reporter,
) -> SampleQC:
    """Validate and profile every mate file of one sample."""
    profiles: list[MateQC] = []
    id_lists: list[list[str]] = []
    for label, path in mates:
        started = time.monotonic()
        reporter.info(
            f"  {label}: profiling {path.name} ({human_bytes(path.stat().st_size)})"
            + ("" if count_all else " (record count skipped)")
        )
        profile, ids = profile_mate(
            path, label, sample_reads=sample_reads, count_all=count_all, reporter=reporter
        )
        profiles.append(profile)
        id_lists.append(ids)
        reporter.record(
            f"    {label}: {profile.sampled_reads:,} reads sampled, "
            f"mean length {profile.mean_length:.1f} bp, GC {profile.gc_percent:.1f}%, "
            f"mean Phred {profile.mean_quality:.1f}, "
            f"quality alphabet {profile.quality_alphabet!r} "
            f"in {human_duration(time.monotonic() - started)}"
        )

    notes: list[str] = []
    paired = len(profiles) >= 2
    identical: bool | None = None
    suffixes = False
    counts_match: bool | None = None

    if paired:
        head = min(len(id_lists[0]), len(id_lists[1]), 100_000)
        identical = id_lists[0][:head] == id_lists[1][:head]
        suffixes = any(i.endswith(("/1", "/2")) for i in id_lists[0][:1000] + id_lists[1][:1000])
        r1, r2 = profiles[0].records, profiles[1].records
        if r1 is not None and r2 is not None:
            counts_match = r1 == r2
            if not counts_match:
                notes.append(
                    f"mate record counts differ: {profiles[0].label}={r1:,} vs "
                    f"{profiles[1].label}={r2:,}. Mate collapse still works — it keys on read id "
                    "— but the files are not a clean pair."
                )
        if identical and not suffixes:
            notes.append(
                "Read ids are identical between mates and carry no /1 or /2 suffix, so the two "
                "mates of a fragment are NOT independent observations. Mate collapse by read id "
                "is required; without it every count would be inflated by up to 2x."
            )
        elif suffixes:
            notes.append("Mate suffixes (/1, /2) detected and stripped during mate collapse.")

    for profile in profiles:
        if profile.min_length != profile.max_length:
            notes.append(
                f"{profile.label}: variable read length "
                f"({profile.min_length}-{profile.max_length} bp) — consistent with reads that "
                "have already been quality- and adapter-trimmed. No trimming is applied."
            )
        if len(profile.quality_alphabet) <= 8:
            notes.append(
                f"{profile.label}: only {len(profile.quality_alphabet)} distinct quality "
                "characters. This is quality binning as used by NovaSeq, not corruption."
            )
        if any(ord(c) < PHRED33_OFFSET for c in profile.quality_alphabet):
            notes.append(
                f"{profile.label}: quality characters below Phred+33 range detected — the "
                "encoding may not be Phred+33."
            )
        if profile.duplicate_fraction > 0.25:
            notes.append(
                f"{profile.label}: {profile.duplicate_fraction:.1%} of sampled reads are exact "
                "sequence duplicates. High duplication reduces effective depth; normalised "
                "figures stay valid because rpoB is affected equally."
            )

    return SampleQC(
        mates=tuple(profiles),
        paired=paired,
        ids_identical_between_mates=identical,
        mate_suffixes_present=suffixes,
        record_counts_match=counts_match,
        notes=tuple(dict.fromkeys(notes)),
    )


# --------------------------------------------------------------------------- #
# caching
# --------------------------------------------------------------------------- #


def _cache_key(mates: list[tuple[str, Path]], sample_reads: int, count_all: bool) -> str:
    parts = [f"v1|{sample_reads}|{count_all}"]
    for label, path in mates:
        stat = path.stat()
        parts.append(f"{label}|{path}|{stat.st_size}|{int(stat.st_mtime)}")
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()[:24]


def profile_sample_cached(
    mates: list[tuple[str, Path]],
    *,
    cache_dir: Path,
    sample_reads: int = DEFAULT_SAMPLE_READS,
    count_all: bool = True,
    reporter: Reporter,
    refresh: bool = False,
) -> SampleQC:
    """``profile_sample`` with an on-disk cache.

    The exact record count is a full pass over the input, so re-running the
    report — which is otherwise seconds once DIAMOND output is cached — would
    otherwise pay for it every time. The key covers each input's size and
    mtime, so editing or replacing a FASTQ invalidates it.
    """
    key = _cache_key(mates, sample_reads, count_all)
    path = cache_dir / f"input_profile_{key}.json"

    if path.is_file() and not refresh:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            qc = _qc_from_json(payload)
        except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
            reporter.record(f"    input profile cache unusable ({exc}); recomputing")
        else:
            reporter.ok(
                f"input validation reused from cache "
                f"({qc.read_pairs:,} read pairs)" if qc.read_pairs else "input validation reused"
            )
            return qc

    qc = profile_sample(
        mates, sample_reads=sample_reads, count_all=count_all, reporter=reporter
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_qc_to_cache(qc), indent=2), encoding="utf-8")
    return qc


def _qc_to_cache(qc: SampleQC) -> dict[str, Any]:
    return {
        "paired": qc.paired,
        "ids_identical_between_mates": qc.ids_identical_between_mates,
        "mate_suffixes_present": qc.mate_suffixes_present,
        "record_counts_match": qc.record_counts_match,
        "notes": list(qc.notes),
        "mates": [asdict(m) for m in qc.mates],
    }


def _qc_from_json(payload: dict[str, Any]) -> SampleQC:
    return SampleQC(
        mates=tuple(MateQC(**m) for m in payload["mates"]),
        paired=bool(payload["paired"]),
        ids_identical_between_mates=payload["ids_identical_between_mates"],
        mate_suffixes_present=bool(payload["mate_suffixes_present"]),
        record_counts_match=payload["record_counts_match"],
        notes=tuple(payload.get("notes", ())),
        cached=True,
    )
