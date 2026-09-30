"""Read preprocessing and sequencing QC via fastp, plus the spec 4.2 gates.

``fastp`` (Chen et al., Bioinformatics 2018) is a C++ FASTQ preprocessor. It
is used here as an external binary, not as a library: its value is a single
multithreaded pass that produces adapter content, duplication rate, Q20/Q30,
insert-size estimate, per-cycle quality and GC — everything the sequencing
quality tile needs — in about the time it takes to decompress the input.

Two modes:

* **assess** (default) — read both mates once, write the JSON report, write no
  reads. Downstream engines see the original files, so every DIAMOND and
  MetaPhlAn cache stays valid.
* **trim** (``--trim``) — additionally write adapter- and polyG-trimmed
  FASTQs and hand those downstream. This changes every downstream result, so
  it is opt-in and the trimmed files carry their own cache fingerprint.

Also here: run/flowcell/lane provenance from Illumina read headers (spec 4.2:
"batch effects routinely exceed the biological differences being measured"),
and the QC-gate table that the report's sequencing-quality section shows.

A note on adapter detection. ``--detect_adapter_for_pe`` in fastp 1.x is a
single-threaded pre-pass that took 49 s on a 200 k-pair prefix here, against
1 s for the whole run without it. Paired-end mode already trims adapters from
overlap analysis, which needs no adapter sequence, so detection is off by
default and can be requested with ``adapter_detect=True``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import DependencyError, OpenBiotaError
from openbiota.logging_util import Reporter, human_duration

CACHE_VERSION: Final = 1

#: Spec 4.2 prespecified minimum usable non-host read pairs.
MIN_USABLE_NONHOST_PAIRS: Final = 500_000
#: Below this, results are shown but every taxonomic non-detection is
#: ``missing`` rather than absence (spec 7.2).
LOW_DEPTH_WARN_PAIRS: Final = 1_000_000
#: Duplication above this suggests a low-complexity library or over-sequencing.
DUP_RATE_WARN: Final = 0.20
DUP_RATE_FAIL: Final = 0.50
Q30_WARN: Final = 0.80
Q30_FAIL: Final = 0.70
HOST_FRACTION_WARN: Final = 0.05
HOST_FRACTION_FAIL: Final = 0.50
ADAPTER_READS_WARN: Final = 0.05
#: Stool metagenomes sit near 45-50% GC; far outside that suggests contamination
#: or a wrong sample type.
GC_LOW_WARN: Final = 0.38
GC_HIGH_WARN: Final = 0.58
#: rpoB fragments needed for a stable copies-per-100-genomes denominator.
RPOB_WARN: Final = 2_000
RPOB_FAIL: Final = 500


# --------------------------------------------------------------------------- #
# locating fastp
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class FastpTools:
    fastp: str
    version: str


def locate_fastp(explicit: str | None = None) -> FastpTools:
    """Find a fastp executable and read its version."""
    candidate = explicit or shutil.which("fastp")
    if candidate is None or not Path(candidate).exists() and shutil.which(candidate) is None:
        raise DependencyError(
            "fastp not found. Install with `brew install fastp` or `conda install -c "
            "bioconda fastp`, or pass --fastp /path/to/fastp. Run with --no-fastp to skip "
            "read preprocessing (sequencing-quality metrics will come from the built-in "
            "sampler only)."
        )
    proc = subprocess.run([candidate, "--version"], capture_output=True, text=True, check=False)
    text = (proc.stdout or proc.stderr).strip()
    m = re.search(r"fastp\s+v?(\d+(?:\.\d+)+)", text)
    if not m:
        raise DependencyError(f"could not read fastp version from {candidate!r}: {text!r}")
    return FastpTools(fastp=candidate, version=m.group(1))


# --------------------------------------------------------------------------- #
# header provenance
# --------------------------------------------------------------------------- #

#: Illumina CASAVA 1.8+ header: @instrument:run:flowcell:lane:tile:x:y [read:filter:control:index]
_ILLUMINA = re.compile(
    r"^@(?P<instrument>[^:\s]+):(?P<run>\d+):(?P<flowcell>[^:\s]+):(?P<lane>\d+):(?P<tile>\d+):"
    r"(?P<x>\d+):(?P<y>\d+)(?::[ACGTN+]+)?(?:\s+(?P<read>\d):(?P<filter>[YN]):(?P<control>\d+):(?P<index>\S*))?"
)
#: SRA-normalised header: @SRR1234567.1 1 length=151
_SRA = re.compile(r"^@(?P<accession>[SED]RR\d+)\.(?P<spot>\d+)")


@dataclass(frozen=True, slots=True)
class HeaderProvenance:
    """Where the reads came from, as far as the FASTQ header says."""

    style: str  # "illumina" | "sra" | "unknown"
    instrument: str | None = None
    run: str | None = None
    flowcells: tuple[str, ...] = ()
    lanes: tuple[int, ...] = ()
    #: Index (barcode) sequences seen in the header, when present.
    indexes: tuple[str, ...] = ()
    sra_accession: str | None = None
    n_examined: int = 0
    example: str = ""

    @property
    def batch_key(self) -> str | None:
        """One string per sequencing batch: flowcell + lane set. Two samples
        sharing it were sequenced together; two that differ carry a batch
        difference that can exceed the biology (spec 4.2)."""
        if self.style != "illumina" or not self.flowcells:
            return None
        return "+".join(self.flowcells) + ":" + ",".join(str(lane) for lane in self.lanes)

    def to_json(self) -> dict[str, Any]:
        return {
            "header_style": self.style,
            "instrument": self.instrument,
            "run": self.run,
            "flowcells": list(self.flowcells),
            "lanes": list(self.lanes),
            "indexes": list(self.indexes),
            "sra_accession": self.sra_accession,
            "batch_key": self.batch_key,
            "headers_examined": self.n_examined,
            "example_header": self.example,
        }


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_header_provenance(path: Path, *, n_reads: int = 2_000) -> HeaderProvenance:
    """Parse the first ``n_reads`` headers for instrument, run, flowcell, lane."""
    instruments: Counter[str] = Counter()
    runs: Counter[str] = Counter()
    flowcells: Counter[str] = Counter()
    lanes: Counter[int] = Counter()
    indexes: Counter[str] = Counter()
    sra: str | None = None
    example = ""
    seen = 0
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            if seen >= n_reads:
                break
            seen += 1
            line = line.rstrip("\n")
            if not example:
                example = line[:120]
            m = _ILLUMINA.match(line)
            if m:
                instruments[m.group("instrument")] += 1
                runs[m.group("run")] += 1
                flowcells[m.group("flowcell")] += 1
                lanes[int(m.group("lane"))] += 1
                if m.group("index"):
                    indexes[m.group("index")] += 1
                continue
            s = _SRA.match(line)
            if s and sra is None:
                sra = s.group("accession")
    if flowcells:
        return HeaderProvenance(
            style="illumina",
            instrument=instruments.most_common(1)[0][0],
            run=runs.most_common(1)[0][0],
            flowcells=tuple(sorted(flowcells)),
            lanes=tuple(sorted(lanes)),
            indexes=tuple(k for k, _ in indexes.most_common(4)),
            n_examined=seen,
            example=example,
        )
    if sra:
        return HeaderProvenance(style="sra", sra_accession=sra, n_examined=seen, example=example)
    return HeaderProvenance(style="unknown", n_examined=seen, example=example)


# --------------------------------------------------------------------------- #
# fastp
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class FastpResult:
    """The metrics from one fastp pass, parsed from its JSON report."""

    fastp_version: str
    mode: str  # "assess" | "trim"
    sequencing: str
    #: Read *pairs* (fastp reports reads; halved here for paired input).
    pairs_in: int
    pairs_out: int
    bases_in: int
    bases_out: int
    q20_rate: float
    q30_rate: float
    gc_content: float
    read1_mean_length: float
    read2_mean_length: float | None
    duplication_rate: float
    #: Fraction of reads where an adapter was trimmed (overlap analysis).
    adapter_trimmed_fraction: float
    adapter_trimmed_bases: int
    #: Insert size peak among pairs that overlap; None for single-end.
    insert_peak: int | None
    #: Fraction of pairs whose insert could be measured (i.e. that overlap).
    insert_measurable_fraction: float | None
    low_quality_reads: int
    too_many_n_reads: int
    too_short_reads: int
    polyg_trimmed_reads: int
    total_cycles: int
    elapsed_s: float
    cached: bool
    command: tuple[str, ...]
    json_path: Path
    #: Trimmed outputs when mode == "trim".
    trimmed: tuple[Path, ...] = ()

    @property
    def paired(self) -> bool:
        return self.read2_mean_length is not None

    @property
    def pass_fraction(self) -> float:
        return self.pairs_out / self.pairs_in if self.pairs_in else 0.0

    @property
    def insert_note(self) -> str:
        if self.insert_peak is None or self.insert_measurable_fraction is None:
            return "single-end input; insert size not estimable"
        if self.insert_measurable_fraction < 0.25:
            return (
                f"{self.insert_measurable_fraction:.0%} of pairs overlap; most inserts exceed "
                f"{2 * int(self.read1_mean_length)} bp (2× read length), which is the normal "
                "state for a stool library and means very little adapter read-through"
            )
        head = (
            f"{self.insert_measurable_fraction:.0%} of pairs overlap, peak insert {self.insert_peak} bp"
        )
        if self.pre_trimmed:
            return (
                f"{head}; reads average {self.read1_mean_length:.0f} bp against "
                f"{self.total_cycles} cycles with almost no adapter found, so the provider "
                "had already trimmed adapters to the insert"
            )
        return f"{head}; a short-insert library, so adapter trimming matters more here"

    @property
    def pre_trimmed(self) -> bool:
        """Reads shorter than the cycle count with no adapter left to find:
        the sequencing provider trimmed before delivery."""
        return (
            self.total_cycles > 0
            and self.read1_mean_length < self.total_cycles - 1
            and self.adapter_trimmed_fraction < 0.005
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "tool": f"fastp {self.fastp_version}",
            "mode": self.mode,
            "sequencing": self.sequencing,
            "pairs_in": self.pairs_in,
            "pairs_out": self.pairs_out,
            "pass_fraction": round(self.pass_fraction, 5),
            "bases_in": self.bases_in,
            "bases_out": self.bases_out,
            "q20_rate": round(self.q20_rate, 5),
            "q30_rate": round(self.q30_rate, 5),
            "gc_content": round(self.gc_content, 4),
            "read1_mean_length": self.read1_mean_length,
            "read2_mean_length": self.read2_mean_length,
            "total_cycles": self.total_cycles,
            "duplication_rate": round(self.duplication_rate, 5),
            "adapter_trimmed_fraction": round(self.adapter_trimmed_fraction, 5),
            "adapter_trimmed_bases": self.adapter_trimmed_bases,
            "insert_peak": self.insert_peak,
            "insert_measurable_fraction": (
                None if self.insert_measurable_fraction is None
                else round(self.insert_measurable_fraction, 4)
            ),
            "insert_note": self.insert_note,
            "low_quality_reads": self.low_quality_reads,
            "too_many_n_reads": self.too_many_n_reads,
            "too_short_reads": self.too_short_reads,
            "polyg_trimmed_reads": self.polyg_trimmed_reads,
            "elapsed_s": round(self.elapsed_s, 1),
            "cached": self.cached,
            "command": list(self.command),
            "trimmed_outputs": [str(p) for p in self.trimmed],
        }


def _fingerprint(paths: list[Path], extra: str) -> str:
    h = hashlib.sha1()
    for p in paths:
        st = p.stat()
        h.update(f"{p.resolve()}|{st.st_size}|{int(st.st_mtime)}|".encode())
    h.update(extra.encode())
    h.update(f"|v{CACHE_VERSION}".encode())
    return h.hexdigest()[:16]


def parse_fastp_json(payload: dict[str, Any], *, mode: str, json_path: Path,
                     elapsed_s: float, cached: bool, command: tuple[str, ...],
                     trimmed: tuple[Path, ...] = ()) -> FastpResult:
    summary = payload["summary"]
    before = summary["before_filtering"]
    after = summary["after_filtering"]
    paired = "read2_before_filtering" in payload
    divisor = 2 if paired else 1
    filt = payload.get("filtering_result", {})
    dup = payload.get("duplication", {}).get("rate", 0.0)
    ins = payload.get("insert_size", {})
    adapters = payload.get("adapter_cutting", {})
    r1_before = payload.get("read1_before_filtering", {})
    total_reads = int(before["total_reads"])
    pairs_in = total_reads // divisor
    hist = ins.get("histogram") or []
    measurable = sum(hist)
    unknown = int(ins.get("unknown", 0))
    insert_frac = (measurable / (measurable + unknown)) if paired and (measurable + unknown) else None
    polyg = int(payload.get("polyx_trimming", {}).get("read1_polyx_trimmed_reads", 0) or 0)
    polyg += int(r1_before.get("polyg_trimmed_reads", 0) or 0)
    return FastpResult(
        fastp_version=str(summary.get("fastp_version", "?")),
        mode=mode,
        sequencing=str(summary.get("sequencing", "")),
        pairs_in=pairs_in,
        pairs_out=int(after["total_reads"]) // divisor,
        bases_in=int(before["total_bases"]),
        bases_out=int(after["total_bases"]),
        q20_rate=float(before["q20_rate"]),
        q30_rate=float(before["q30_rate"]),
        gc_content=float(before["gc_content"]),
        read1_mean_length=float(before.get("read1_mean_length", 0)),
        read2_mean_length=float(before["read2_mean_length"]) if paired else None,
        duplication_rate=float(dup),
        adapter_trimmed_fraction=(
            int(adapters.get("adapter_trimmed_reads", 0)) / total_reads if total_reads else 0.0
        ),
        adapter_trimmed_bases=int(adapters.get("adapter_trimmed_bases", 0)),
        insert_peak=int(ins["peak"]) if paired and "peak" in ins else None,
        insert_measurable_fraction=insert_frac,
        low_quality_reads=int(filt.get("low_quality_reads", 0)),
        too_many_n_reads=int(filt.get("too_many_N_reads", 0)),
        too_short_reads=int(filt.get("too_short_reads", 0)),
        polyg_trimmed_reads=polyg,
        total_cycles=int(r1_before.get("total_cycles", 0)),
        elapsed_s=elapsed_s,
        cached=cached,
        command=command,
        json_path=json_path,
        trimmed=trimmed,
    )


def run_fastp(
    tools: FastpTools,
    *,
    r1: Path,
    r2: Path | None,
    work_dir: Path,
    reporter: Reporter,
    threads: int = 8,
    trim: bool = False,
    adapter_detect: bool = False,
    min_length: int = 60,
    force: bool = False,
) -> FastpResult:
    """Run fastp once over the sample; cache the JSON (and trimmed reads) by input.

    ``threads`` is capped at 16, fastp's own maximum.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    mode = "trim" if trim else "assess"
    inputs = [r1] + ([r2] if r2 else [])
    key = _fingerprint(inputs, f"{mode}|{adapter_detect}|{min_length}|{tools.version}")
    json_path = work_dir / f"fastp.{mode}.json"
    stamp = work_dir / f"fastp.{mode}.stamp.json"
    out1 = work_dir / "trimmed.1.fastq.gz"
    out2 = work_dir / "trimmed.2.fastq.gz"
    trimmed: tuple[Path, ...] = ()
    if trim:
        trimmed = (out1, out2) if r2 else (out1,)

    if not force and json_path.is_file() and stamp.is_file():
        meta = json.loads(stamp.read_text())
        if meta.get("key") == key and all(p.is_file() for p in trimmed):
            payload = json.loads(json_path.read_text())
            result = parse_fastp_json(
                payload, mode=mode, json_path=json_path, elapsed_s=meta.get("elapsed_s", 0.0),
                cached=True, command=tuple(meta.get("command", ())), trimmed=trimmed,
            )
            reporter.record(f"    fastp {mode}: cached ({json_path.name})")
            return result

    command = [
        tools.fastp,
        "-i", str(r1),
        "-w", str(max(1, min(16, threads))),
        "-j", str(json_path),
        "-h", str(work_dir / f"fastp.{mode}.html"),
        "--trim_poly_g",
        "--length_required", str(min_length),
        "--report_title", f"openbiota {mode}",
    ]
    if r2 is not None:
        command += ["-I", str(r2)]
        if adapter_detect:
            command.append("--detect_adapter_for_pe")
    else:
        # fastp auto-detects adapters for single-end input unless told not to;
        # that pre-pass is the slow part, so honour the same switch.
        if not adapter_detect:
            command.append("--disable_adapter_trimming")
    if trim:
        command += ["-o", str(out1)]
        if r2 is not None:
            command += ["-O", str(out2)]
        command += ["--compression", "4"]

    reporter.info(f"  fastp {tools.version} ({mode})")
    started = time.monotonic()
    proc = subprocess.run(command, capture_output=True, text=True, check=False)
    elapsed = time.monotonic() - started
    if proc.returncode != 0 or not json_path.is_file():
        raise OpenBiotaError(f"fastp failed (exit {proc.returncode}):\n{proc.stderr[-3000:]}")
    payload = json.loads(json_path.read_text())
    stamp.write_text(
        json.dumps({"key": key, "elapsed_s": round(elapsed, 1),
                    "command": [Path(command[0]).name, *command[1:]]})
    )
    result = parse_fastp_json(
        payload, mode=mode, json_path=json_path, elapsed_s=elapsed, cached=False,
        command=(Path(command[0]).name, *command[1:]), trimmed=trimmed,
    )
    reporter.ok(
        f"fastp: {result.pairs_in:,} pairs, Q30 {result.q30_rate:.1%}, "
        f"dup {result.duplication_rate:.2%}, adapters in {result.adapter_trimmed_fraction:.2%} "
        f"of reads, {human_duration(elapsed)}"
    )
    return result


# --------------------------------------------------------------------------- #
# QC gates (spec 4.2)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class QCGate:
    """One row of the sequencing-quality table."""

    name: str
    label: str
    value: float | int | str | None
    display: str
    status: str  # "pass" | "warn" | "fail" | "unknown" | "not_assessable"
    note: str

    def to_json(self) -> dict[str, Any]:
        return {
            "gate": self.name,
            "label": self.label,
            "value": self.value,
            "display": self.display,
            "status": self.status,
            "note": self.note,
        }


@dataclass(slots=True)
class QCGates:
    gates: list[QCGate]
    provenance: HeaderProvenance | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def overall(self) -> str:
        statuses = {g.status for g in self.gates}
        if "fail" in statuses:
            return "fail"
        if "warn" in statuses:
            return "warn"
        return "pass"

    @property
    def usable_pairs(self) -> int | None:
        g = next((g for g in self.gates if g.name == "nonhost_pairs"), None)
        return None if g is None or not isinstance(g.value, int) else g.value

    @property
    def depth_adequate(self) -> bool | None:
        """Spec 4.2: at least 500,000 usable non-host pairs."""
        if self.usable_pairs is None:
            return None
        return self.usable_pairs >= MIN_USABLE_NONHOST_PAIRS

    def to_json(self) -> dict[str, Any]:
        return {
            "overall": self.overall,
            "usable_nonhost_pairs": self.usable_pairs,
            "depth_adequate": self.depth_adequate,
            "minimum_usable_pairs": MIN_USABLE_NONHOST_PAIRS,
            "gates": [g.to_json() for g in self.gates],
            "provenance": None if self.provenance is None else self.provenance.to_json(),
            "notes": list(self.notes),
        }


def _gate(name: str, label: str, value: Any, display: str, status: str, note: str) -> QCGate:
    return QCGate(name, label, value, display, status, note)


def assemble_gates(
    *,
    input_pairs: int | None,
    fastp: FastpResult | None,
    host_fraction: float | None,
    nonhost_pairs: int | None,
    classified_fraction: float | None,
    host_status: str | None = None,
    rpob_fragments: int | None,
    sampled_duplicate_fraction: float | None = None,
    sampled_gc_percent: float | None = None,
    sampled_mean_length: float | None = None,
    provenance: HeaderProvenance | None = None,
) -> QCGates:
    """Build the spec 4.2 gate table from whatever stages ran.

    Every value is sourced from the stage that measured it; where fastp did
    not run, the built-in sampler's estimates stand in and the note says so.
    """
    gates: list[QCGate] = []
    notes: list[str] = []

    # Input pairs
    if input_pairs is not None:
        status = "pass" if input_pairs >= LOW_DEPTH_WARN_PAIRS else "warn" if input_pairs >= MIN_USABLE_NONHOST_PAIRS else "fail"
        gates.append(_gate("input_pairs", "Input read pairs", input_pairs, f"{input_pairs:,}", status,
                           "read pairs in the FASTQ files as delivered"))
    else:
        gates.append(_gate("input_pairs", "Input read pairs", None, "—", "unknown", "not counted"))

    # Post-trim pairs
    if fastp is not None:
        frac = fastp.pass_fraction
        status = "pass" if frac >= 0.95 else "warn" if frac >= 0.85 else "fail"
        gates.append(_gate("post_trim_pairs", "Pairs passing quality filters", fastp.pairs_out,
                           f"{fastp.pairs_out:,} ({frac:.1%})", status,
                           "fastp: length ≥ 60 after polyG and adapter trimming, ≤ 5 N, mean Q ≥ 15"
                           + ("; trimmed reads were used downstream" if fastp.mode == "trim"
                              else "; assessment only — original reads went downstream")))
    else:
        gates.append(_gate("post_trim_pairs", "Pairs passing quality filters", None, "—", "unknown",
                           "fastp did not run"))

    # Host fraction
    if host_status == "not_assessable_upstream_removed" or (
        host_fraction == 0.0 and (input_pairs or 0) >= 1_000_000
    ):
        # Real stool always carries some human DNA. Not one pair means the
        # provider removed human reads before delivery, which is standard for
        # consumer kits; the figure then says nothing about the specimen and
        # is recorded as not assessable rather than as a 0.00% pass.
        gates.append(_gate("host_fraction", "Human DNA fraction", None, "not assessable", "not_assessable",
                           "no human reads at all — the sequencing provider removed them before "
                           "delivery (routine for consumer kits), so the host fraction of the original "
                           "specimen is not measurable here (not_assessable_upstream_removed)"))
    elif host_fraction is not None:
        status = "pass" if host_fraction < HOST_FRACTION_WARN else "warn" if host_fraction < HOST_FRACTION_FAIL else "fail"
        note = ("pairs aligning to GRCh38 (Bowtie2); healthy stool is typically < 1%, "
                "inflamed or bloody stool higher")
        gates.append(_gate("host_fraction", "Human DNA fraction", host_fraction, f"{host_fraction:.2%}", status, note))
    else:
        gates.append(_gate("host_fraction", "Human DNA fraction", None, "—", "unknown",
                           "host filter did not run"))

    # Non-host pairs — the prespecified gate
    if nonhost_pairs is not None:
        status = "pass" if nonhost_pairs >= LOW_DEPTH_WARN_PAIRS else "warn" if nonhost_pairs >= MIN_USABLE_NONHOST_PAIRS else "fail"
        note = f"prespecified minimum {MIN_USABLE_NONHOST_PAIRS:,}"
        if nonhost_pairs < MIN_USABLE_NONHOST_PAIRS:
            note += "; below it, a species not detected is 'missing', not absent"
        gates.append(_gate("nonhost_pairs", "Usable non-host pairs", nonhost_pairs, f"{nonhost_pairs:,}", status, note))
    else:
        gates.append(_gate("nonhost_pairs", "Usable non-host pairs", None, "—", "unknown", "host filter did not run"))

    # Classified fraction
    if classified_fraction is not None:
        # Catalogue coverage, not sequencing quality: MetaPhlAn 3 leaves
        # 30-55% of healthy stool unclassified. Only a very low fraction —
        # the shape of a non-stool or heavily contaminated sample — fails.
        status = "pass" if classified_fraction >= 0.50 else "warn" if classified_fraction >= 0.30 else "fail"
        gates.append(_gate("classified_fraction", "Reads assigned to a known species", classified_fraction,
                           f"{classified_fraction:.1%}", status,
                           "the remainder is organisms without a reference genome, "
                           "not a quality problem — 30-55% unclassified is normal for stool"))
    else:
        gates.append(_gate("classified_fraction", "Reads assigned to a known species", None, "—", "unknown",
                           "taxonomic engine did not run"))

    # Duplicate rate
    dup = fastp.duplication_rate if fastp is not None else sampled_duplicate_fraction
    if dup is not None:
        status = "pass" if dup < DUP_RATE_WARN else "warn" if dup < DUP_RATE_FAIL else "fail"
        src = "fastp, whole file" if fastp is not None else "built-in sampler, prefix estimate"
        gates.append(_gate("duplicate_rate", "Duplicate read rate", dup, f"{dup:.2%}", status,
                           f"{src}; high duplication means less unique information than the read count suggests"))
    else:
        gates.append(_gate("duplicate_rate", "Duplicate read rate", None, "—", "unknown", "not measured"))

    # Base quality
    if fastp is not None:
        q30 = fastp.q30_rate
        status = "pass" if q30 >= Q30_WARN else "warn" if q30 >= Q30_FAIL else "fail"
        gates.append(_gate("q30", "Bases at Q30 or better", q30, f"{q30:.1%}", status,
                           "fraction of bases with ≤ 0.1% error probability"))
        adapt = fastp.adapter_trimmed_fraction
        status = "pass" if adapt < ADAPTER_READS_WARN else "warn"
        gates.append(_gate("adapter_content", "Reads with adapter sequence", adapt, f"{adapt:.2%}", status,
                           "detected by read-pair overlap; " + fastp.insert_note))

    # Read length
    length = fastp.read1_mean_length if fastp is not None else sampled_mean_length
    if length is not None:
        status = "pass" if length >= 100 else "warn" if length >= 75 else "fail"
        gates.append(_gate("read_length", "Mean read length", length, f"{length:.0f} bp", status,
                           "translated search needs ≥ 75 bp to place a read on a protein; ≥ 100 bp preferred"))

    # GC
    gc = fastp.gc_content if fastp is not None else (sampled_gc_percent / 100.0 if sampled_gc_percent is not None else None)
    if gc is not None:
        status = "pass" if GC_LOW_WARN <= gc <= GC_HIGH_WARN else "warn"
        gates.append(_gate("gc_content", "GC content", gc, f"{gc:.1%}", status,
                           "stool metagenomes sit near 45-50%; far outside suggests contamination or a different sample type"))

    # Contamination flags
    flags: list[str] = []
    if host_fraction is not None and host_fraction >= HOST_FRACTION_FAIL:
        flags.append("host DNA dominates")
    if gc is not None and not (GC_LOW_WARN <= gc <= GC_HIGH_WARN):
        flags.append("GC outside stool range")
    if classified_fraction is not None and classified_fraction < 0.40:
        flags.append("low classified fraction")
    gates.append(_gate("contamination", "Contamination flags", len(flags),
                       "; ".join(flags) if flags else "none raised",
                       "warn" if flags else "pass",
                       "composite of host fraction, GC and classified fraction"))

    # rpoB
    if rpob_fragments is not None:
        status = "pass" if rpob_fragments >= RPOB_WARN else "warn" if rpob_fragments >= RPOB_FAIL else "fail"
        gates.append(_gate("rpob_fragments", "rpoB fragments (genome denominator)", rpob_fragments,
                           f"{rpob_fragments:,}", status,
                           "single-copy marker used to express gene counts per 100 genomes; "
                           f"≥ {RPOB_WARN:,} for a stable denominator"))

    # Batch
    if provenance is not None:
        if provenance.style == "illumina":
            gates.append(_gate("batch", "Sequencing batch", provenance.batch_key,
                               f"{provenance.instrument} run {provenance.run}, flowcell "
                               f"{'+'.join(provenance.flowcells)}, lane{'s' if len(provenance.lanes) > 1 else ''} "
                               f"{','.join(str(lane) for lane in provenance.lanes)}",
                               "pass", "from read headers; samples from different flowcells carry a batch difference"))
        else:
            gates.append(_gate("batch", "Sequencing batch", None,
                               "not recoverable from headers" if provenance.style == "unknown"
                               else f"SRA {provenance.sra_accession} (headers normalised)",
                               "unknown", "run and flowcell could not be read from the FASTQ headers"))

    return QCGates(gates=gates, provenance=provenance, notes=notes)


__all__ = [
    "FastpResult",
    "FastpTools",
    "HeaderProvenance",
    "MIN_USABLE_NONHOST_PAIRS",
    "QCGate",
    "QCGates",
    "assemble_gates",
    "locate_fastp",
    "parse_fastp_json",
    "read_header_provenance",
    "run_fastp",
]
