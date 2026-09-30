"""The pathogen branch: reads in, adjudicated findings out (spec §5.2, §6, §13.3).

Two stages, for a reason. A single pass that both finds and confirms would
have to choose between speed and rigour; splitting them lets the first stage
be permissive and cheap over the full library, and the second stage be
expensive and strict over the handful of reads that matter.

**Stage 1 — candidate generation.** Every QC-passed read is hashed against the
frozen bundle index. No abundance filter, no downsampling, unclassified reads
preserved. Output is a list of targets worth aligning. This is explicitly not
a call: k-mer classification is candidate generation, and treating its output
as positive is the documented cause of pathogen false positives.

**Stage 2 — competitive confirmation.** Candidate reads *and their mates* are
realigned against the whole bundle, so every target competes with its
relatives, the host and the food decoys simultaneously. Only here can a
finding become supported.

Then adjudication, which is deliberately dull: `openbiota.pathogens.kernel` decides
the status, this module only supplies evidence and wording. The wording
functions are as load-bearing as the arithmetic — "Candida DNA detected" and
"you have a Candida infection" are different claims, and the whole point of
the branch is to make only the first one.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import numpy as np

from openbiota.errors import DependencyError, OpenBiotaError
from openbiota.pathogens.align import (
    AlignmentEvidence,
    TargetEvidence,
    iter_sam,
    run_bowtie2,
    summarise_alignments,
)
from openbiota.pathogens.catalog import CallingProfile, Catalog, load_calling_profile
from openbiota.pathogens.eligibility import EligibilityRecord, check_assay_eligibility
from openbiota.pathogens.kernel import (
    SUPPORTED_STATUSES,
    NormalizedEvidence,
    SupportPolicy,
    display_precedence,
    sequence_status,
)
from openbiota.pathogens.kmers import SHARED, K, KmerIndex, select_kmers
from openbiota.pathogens.refs import BundleManifest, build_competitive_index
from openbiota.pathogens.schema import (
    AssayManifest,
    CoverageLedger,
    DeterminantResult,
    PathogenResult,
    TargetRecord,
)

__all__ = [
    "PathogenBranchResult",
    "compile_coverage",
    "discover_candidates",
    "interpret_finding",
    "run_pathogen_branch",
    "self_test",
]

#: A read must share at least this many distinct sampled k-mers with one
#: target before that target is worth aligning. Two is the spec's
#: `minimum_hit_groups`; a single shared 31-mer is noise at this scale.
MIN_READ_KMERS: Final[int] = 2

#: A target enters the competitive stage once this many reads nominate it.
#: Set to one so nothing is filtered out by abundance — the gate exists only
#: to avoid aligning against targets no read mentioned at all.
MIN_CANDIDATE_READS: Final[int] = 1

#: Reads are hashed in blocks of roughly this many bases. One vectorised call
#: per block instead of one per read is what makes a full-depth pass cost
#: about a minute rather than an hour.
BLOCK_BASES: Final[int] = 24_000_000

#: Separator inserted between concatenated reads. Longer than k, so no k-mer
#: can span two reads, and non-ACGT so the window is dropped outright.
_SEP: Final[bytes] = b"N" * (K + 1)

#: Normalisation denominator label (spec §5.5).
_DENOM: Final = "qc_pass_nonhost_fragments"


# --------------------------------------------------------------------------- #
# FASTQ streaming
# --------------------------------------------------------------------------- #


def _open_reads(path: Path) -> Iterator[bytes]:
    """Yield raw lines from a FASTQ, decompressing out of process when possible.

    `gzip -dc` (or `pigz`) in a pipe moves decompression onto another core,
    which roughly halves wall time on the full-depth pass. The pure-Python
    fallback keeps the branch runnable without either binary.
    """
    path = Path(path)
    if path.suffix != ".gz":
        with path.open("rb") as handle:
            yield from handle
        return
    exe = shutil.which("pigz") or shutil.which("gzip")
    if exe is None:
        with gzip.open(path, "rb") as handle:
            yield from handle
        return
    argv = [exe, "-dc", str(path)]
    proc = subprocess.Popen(  # noqa: S603 - fixed argv, no shell
        argv, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, bufsize=1 << 22
    )
    try:
        assert proc.stdout is not None
        yield from proc.stdout
    finally:
        if proc.poll() is None:
            proc.terminate()
        proc.wait(timeout=30)


@dataclass
class _Block:
    """One batch of reads, ready for a single vectorised hash pass."""

    names: list[str]
    seqs: list[bytes]
    quals: list[bytes]
    mate: list[int]


def _iter_blocks(
    r1: Path, r2: Path | None, *, block_bases: int = BLOCK_BASES
) -> Iterator[_Block]:
    """Stream mate-synchronised read blocks.

    Mates travel together so that a hit on either mate retains both. Losing
    the mate of a candidate read would discard exactly the evidence the
    competitive stage needs.
    """
    it1 = _open_reads(r1)
    it2 = _open_reads(r2) if r2 is not None else None
    block = _Block([], [], [], [])
    bases = 0
    while True:
        try:
            name1 = next(it1)
        except StopIteration:
            break
        try:
            seq1 = next(it1).rstrip()
            next(it1)
            qual1 = next(it1).rstrip()
        except StopIteration:
            break
        stem = name1[1:].split()[0].decode("ascii", "replace")
        block.names.append(stem)
        block.seqs.append(seq1)
        block.quals.append(qual1)
        block.mate.append(1)
        bases += len(seq1)
        if it2 is not None:
            try:
                next(it2)
                seq2 = next(it2).rstrip()
                next(it2)
                qual2 = next(it2).rstrip()
            except StopIteration:
                it2 = None
            else:
                block.names.append(stem)
                block.seqs.append(seq2)
                block.quals.append(qual2)
                block.mate.append(2)
                bases += len(seq2)
        if bases >= block_bases:
            yield block
            block = _Block([], [], [], [])
            bases = 0
    if block.names:
        yield block


# --------------------------------------------------------------------------- #
# Stage 1: candidate generation
# --------------------------------------------------------------------------- #


@dataclass
class CandidateSet:
    """Stage-1 output: what to align, and the raw k-mer support behind it."""

    reads_examined: int = 0
    fragments_examined: int = 0
    bases_examined: int = 0
    unclassified_reads: int = 0
    shared_only_reads: int = 0
    per_target_reads: dict[str, int] = field(default_factory=dict)
    candidate_fastq: tuple[Path, ...] = ()
    candidate_fragments: int = 0
    elapsed_s: float = 0.0

    def to_json(self) -> dict[str, Any]:
        return {
            "reads_examined": self.reads_examined,
            "fragments_examined": self.fragments_examined,
            "bases_examined": self.bases_examined,
            "unclassified_reads": self.unclassified_reads,
            "shared_only_reads": self.shared_only_reads,
            "candidate_fragments": self.candidate_fragments,
            "candidate_targets": len(self.per_target_reads),
            "per_target_reads": dict(sorted(self.per_target_reads.items())),
            "elapsed_s": round(self.elapsed_s, 1),
        }


def discover_candidates(
    r1: Path,
    r2: Path | None,
    index: KmerIndex,
    *,
    work_dir: Path,
    progress: Callable[[str], None] | None = None,
) -> CandidateSet:
    """Hash every read against the bundle and collect candidate fragments.

    Reads whose only support is `SHARED` sequence are counted and retained as
    candidate material rather than discarded: they are exactly the evidence
    that has to appear as ambiguous at a higher rank instead of vanishing.
    """
    say = progress or (lambda _m: None)
    started = time.monotonic()
    out = CandidateSet()
    work_dir.mkdir(parents=True, exist_ok=True)
    cand1 = work_dir / "candidates_R1.fastq"
    cand2 = work_dir / "candidates_R2.fastq" if r2 is not None else None
    handles = [cand1.open("wb")]
    if cand2 is not None:
        handles.append(cand2.open("wb"))

    n_targets = len(index.target_ids)
    read_hits = np.zeros(n_targets, dtype=np.int64)

    try:
        for block_no, block in enumerate(_iter_blocks(r1, r2), 1):
            n = len(block.seqs)
            if not n:
                continue
            lengths = np.fromiter((len(s) for s in block.seqs), dtype=np.int64, count=n)
            # Offsets of each read inside the concatenated block.
            starts = np.zeros(n, dtype=np.int64)
            np.cumsum(lengths[:-1] + len(_SEP), out=starts[1:])
            joined = _SEP.join(block.seqs)

            kmers, positions = select_kmers(joined, positions=True)
            out.reads_examined += n
            out.bases_examined += int(lengths.sum())
            if kmers.size == 0:
                out.unclassified_reads += n
                continue

            assert positions is not None
            owners = index.lookup(kmers)
            # Map each retained k-mer back to the read it came from.
            read_of = np.searchsorted(starts, positions, side="right") - 1

            resolved = owners >= 0
            specific = resolved & (owners != SHARED)

            # Distinct (read, target) support. `unique` over the pair is the
            # per-read hit-group count the profile's minimum applies to.
            if specific.any():
                pairs = np.stack(
                    (read_of[specific], owners[specific].astype(np.int64)), axis=1
                )
                uniq, counts = np.unique(pairs, axis=0, return_counts=True)
                keep = counts >= MIN_READ_KMERS
                nominated_reads = uniq[keep, 0]
                nominated_targets = uniq[keep, 1]
                if nominated_targets.size:
                    np.add.at(read_hits, nominated_targets, 1)
            else:
                nominated_reads = np.empty(0, dtype=np.int64)

            shared_reads = np.unique(read_of[resolved & (owners == SHARED)])
            interesting = np.unique(nominated_reads)
            out.shared_only_reads += int(
                np.setdiff1d(shared_reads, interesting, assume_unique=False).size
            )
            out.unclassified_reads += n - int(np.unique(read_of[resolved]).size)

            if interesting.size:
                # Retain the nominated read and its mate: fragments, not reads.
                names = np.array(block.names, dtype=object)
                wanted_names = set(names[interesting].tolist())
                for i in range(n):
                    if block.names[i] in wanted_names:
                        handle = handles[block.mate[i] - 1] if len(handles) > 1 else handles[0]
                        handle.write(
                            b"@"
                            + block.names[i].encode()
                            + b"\n"
                            + block.seqs[i]
                            + b"\n+\n"
                            + block.quals[i]
                            + b"\n"
                        )
                out.candidate_fragments += len(wanted_names)
            if block_no % 10 == 0:
                say(
                    f"candidate scan: {out.reads_examined:,} reads, "
                    f"{out.candidate_fragments:,} candidate fragments"
                )
    finally:
        for handle in handles:
            handle.close()

    out.fragments_examined = (
        out.reads_examined // 2 if r2 is not None else out.reads_examined
    )
    for i, target_id in enumerate(index.target_ids):
        if read_hits[i] >= MIN_CANDIDATE_READS:
            out.per_target_reads[target_id] = int(read_hits[i])
    out.candidate_fastq = tuple(
        p for p in (cand1, cand2) if p is not None and p.stat().st_size
    )
    out.elapsed_s = time.monotonic() - started
    say(
        f"candidate scan complete: {out.reads_examined:,} reads in "
        f"{out.elapsed_s:.0f}s, {len(out.per_target_reads)} candidate targets"
    )
    return out


# --------------------------------------------------------------------------- #
# Interpretation wording (spec §13.3)
# --------------------------------------------------------------------------- #

_CLASS_MEANING: Final[dict[str, str]] = {
    "established_enteric": (
        "an organism capable of causing gastrointestinal infection"
    ),
    "toxin_or_pathotype_dependent": (
        "an organism whose disease potential depends on strain traits that are "
        "assessed separately"
    ),
    "conditional_opportunist": (
        "an organism that commonly lives in or passes through the gut and "
        "causes disease only in particular hosts or sites"
    ),
    "rare_enteric": (
        "an uncommon but documented cause of human intestinal or hepatobiliary "
        "infection"
    ),
    "uncertain_enteric_role": (
        "an organism whose role in disease is uncertain or disputed"
    ),
    "extraintestinal_watch": (
        "an organism for which stool is not a diagnostic specimen"
    ),
    "background_or_decoy": (
        "a commensal, dietary, environmental or harmless related organism"
    ),
}


def interpret_finding(
    target: TargetRecord,
    status: str,
    eligibility: EligibilityRecord,
    *,
    marker_only: bool = False,
    ambiguous_with: str | None = None,
    contamination: str = "controls_unavailable",
) -> str:
    """The participant-facing sentence for one result.

    Every branch here is constrained by the spec's approved wording table. The
    rules that generate the awkward-sounding sentences are the important ones:
    a non-detection states its limitations rather than reassuring, and a
    detection names what it does *not* establish.
    """
    name = target.display_name

    if status == "not_assessed":
        return eligibility.statement if eligibility.eligibility != "eligible" else (
            f"Not assessed for {name}: "
            + (target.reference_gap_note or "no usable installed reference")
            + "."
        )

    if status == "not_detected":
        limits = _limitation_sentence(target)
        if not target.excludable_by_stool:
            return (
                f"No supported {name} sequence was found. Stool is not a "
                f"diagnostic specimen for this organism, so this result does "
                f"not screen out or exclude its infection. {limits}"
            )
        return (
            f"No supported {name} sequence detected in this sample. {limits}"
        )

    if status == "ambiguous_signal":
        rank = target.taxonomic_resolution.replace("_", " ")
        extra = (
            f" The evidence does not separate it from {ambiguous_with.split('.')[-1].replace('_', ' ')}."
            if ambiguous_with
            else ""
        )
        return (
            f"Low-level sequence consistent with {name} was found, but the "
            f"assignment is unresolved at {rank} level and is reported at the "
            f"rank the evidence supports.{extra}"
        )

    if status == "candidate_signal":
        return (
            f"A weak signal consistent with {name} did not meet the documented "
            f"research support rule. It is retained here as a candidate, not as "
            f"a finding, and is not evidence of absence either."
        )

    # Supported or marker-level.
    route = (
        "single-locus or organelle sequence"
        if marker_only or status == "marker_signal"
        else "DNA sequences"
    )
    lead = f"{route.capitalize()} consistent with {name} detected under the research rule."
    meaning = _CLASS_MEANING.get(target.interpretation_class, "")
    body = f" This is {meaning}." if meaning else ""

    if target.interpretation_class == "background_or_decoy":
        tail = (
            " Detection of this organism is not a pathogen alarm and does not "
            "by itself indicate a parasitic or infectious problem."
        )
    elif target.interpretation_class == "extraintestinal_watch":
        tail = (
            " An unexpected finding of this kind warrants appropriate "
            "confirmation with the correct specimen type; stool sequencing "
            "neither diagnoses nor excludes its infection."
        )
    elif target.interpretation_class == "uncertain_enteric_role":
        tail = (
            " Its clinical significance is uncertain, and this result is not by "
            "itself a reason for treatment."
        )
    elif target.interpretation_class == "toxin_or_pathotype_dependent":
        tail = (
            " Species-level identification alone cannot establish a pathogenic "
            "strain; any toxin or pathotype evidence is reported separately."
        )
    elif target.interpretation_class == "conditional_opportunist":
        tail = (
            " Carriage is common, and this result does not diagnose an "
            "infection or establish a need for treatment."
        )
    else:
        tail = (
            " This sequence finding may warrant confirmation with an "
            "appropriate clinical test; it does not by itself establish the "
            "cause of any symptoms."
        )

    if status == "marker_signal":
        tail += (
            " Support comes from a single locus or organelle reference, which is "
            "reported as marker-level evidence and not as whole-genome support."
        )
    if contamination == "controls_unavailable":
        tail += (
            " This sample has no extraction or library blanks, so contamination "
            "checks could not be completed."
        )
    return lead + body + tail


_LIMITATION_TEXT: Final[dict[str, str]] = {
    "sampling": "a single stool sample may not contain the organism",
    "shedding": "shedding can be low or intermittent",
    "extraction": "extraction may not release the organism's DNA efficiently",
    "depth": "sequencing depth limits sensitivity for low-abundance targets",
    "reference_diversity": "divergent strains may be missed by installed references",
    "reference_absent": "no discriminating reference is installed",
    "no_patent_stool_stage": "human infection may not produce a stool-detectable stage",
    "tissue_restricted": "the infection is tissue-restricted",
    "wrong_nucleic_acid": "this library did not sequence the target's molecule",
    "wrong_specimen_type": "stool is not the appropriate specimen",
    "hard_to_lyse_cell_wall": "the cell wall resists standard lysis",
    "low_abundance": "very low abundance may fall below detection",
}


def _limitation_sentence(target: TargetRecord) -> str:
    """Say which specific limitations bound this target's non-detection.

    Generic caveats are easy to ignore. Naming the mechanisms — and never
    implying that more stool or more reads would fix them — is the point.
    """
    parts = [
        _LIMITATION_TEXT[code]
        for code in target.negative_limitation_codes
        if code in _LIMITATION_TEXT
    ]
    if not parts:
        parts = [_LIMITATION_TEXT["sampling"], _LIMITATION_TEXT["depth"]]
    joined = "; ".join(parts[:4])
    return (
        f"This means no signal met the documented analytical rule, not that the "
        f"organism is absent: {joined}."
    )


# --------------------------------------------------------------------------- #
# Coverage ledger
# --------------------------------------------------------------------------- #


def compile_coverage(
    targets: Sequence[TargetRecord],
    results: Mapping[str, PathogenResult],
    determinants: Sequence[DeterminantResult],
    *,
    failed_routes: Sequence[str] = (),
) -> CoverageLedger:
    """Count what was and was not assessed, by reason.

    Organism targets, determinants and subtypes are tallied separately. Adding
    them together would inflate a coverage claim, which is how a screen ends up
    advertising more than it can do.
    """
    assessed = 0
    by_assay = 0
    by_reference = 0
    by_route = 0
    group_level = 0
    masked = 0
    out_of_scope = 0
    pending_public = 0
    pending_members = 0
    gaps: dict[str, int] = {}
    by_group: dict[str, dict[str, int]] = {}

    for target in targets:
        result = results.get(target.target_id)
        group = by_group.setdefault(
            target.group,
            {"assessed": 0, "not_assessed": 0, "supported": 0, "candidate": 0},
        )
        if result is None:
            by_route += 1
            group["not_assessed"] += 1
            gaps["route_not_run"] = gaps.get("route_not_run", 0) + 1
            continue
        if result.sequence_status == "not_assessed":
            group["not_assessed"] += 1
            if result.assay_eligibility != "eligible":
                by_assay += 1
                # An RNA genome cannot appear in a DNA library. That is the
                # molecule, not a gap in this screen, and it is counted apart
                # from anything a deeper run or a better reference would fix.
                if "wrong_nucleic_acid" in (result.reason_codes or ()):
                    out_of_scope += 1
            elif not target.reference_usable:
                by_reference += 1
                reason = target.reference_status
                gaps[reason] = gaps.get(reason, 0) + 1
                # Installed sequence, aligned against, but every informative
                # k-mer is shared with a catalogue relative: searched, and
                # answerable only at the rank it shares.
                if target.reference_status == "covered_by_relative":
                    group_level += 1
                    group["not_assessed"] -= 1
                    group["searched_group_level"] = group.get("searched_group_level", 0) + 1
                elif target.reference_bases and not target.informative_bases:
                    if target.near_neighbor_target_ids:
                        group_level += 1
                    else:
                        masked += 1
                    group["not_assessed"] -= 1
                    group["searched_group_level"] = group.get("searched_group_level", 0) + 1
                elif target.taxonomic_resolution in {"genus", "group", "species_complex"} and not target.ncbi_taxids:
                    pending_members += 1
                else:
                    pending_public += 1
            else:
                by_route += 1
                gaps["route_not_run"] = gaps.get("route_not_run", 0) + 1
            continue
        assessed += 1
        group["assessed"] += 1
        if result.sequence_status in SUPPORTED_STATUSES:
            group["supported"] += 1
        elif result.sequence_status in {"candidate_signal", "ambiguous_signal"}:
            group["candidate"] += 1

    det_assessed = sum(
        1 for d in determinants if d.determinant_status != "not_assessed"
    )
    return CoverageLedger(
        assessed=assessed,
        searched_group_level=group_level,
        searched_masked_by_host=masked,
        out_of_assay_scope=out_of_scope,
        reference_pending_no_public_sequence=pending_public,
        reference_pending_needs_member_list=pending_members,
        not_assessed_assay=by_assay,
        not_assessed_reference=by_reference,
        not_assessed_route=by_route,
        determinants_assessed=det_assessed,
        determinants_not_assessed=len(determinants) - det_assessed,
        by_group=by_group,
        reference_gap_reasons=dict(sorted(gaps.items())),
        failed_routes=tuple(failed_routes),
        total_targets=len(targets),
    )


# --------------------------------------------------------------------------- #
# Branch result
# --------------------------------------------------------------------------- #


@dataclass
class PathogenBranchResult:
    """Everything the report and the JSON need from one sample's screen."""

    sample_id: str
    assay: AssayManifest
    bundle_id: str
    catalog_version: str
    calling_profile: CallingProfile
    results: dict[str, PathogenResult] = field(default_factory=dict)
    determinants: tuple[DeterminantResult, ...] = ()
    coverage: CoverageLedger | None = None
    candidates: CandidateSet | None = None
    alignment: AlignmentEvidence | None = None
    analysis_status: str = "completed"
    not_run_reason: str | None = None
    failed_routes: tuple[str, ...] = ()
    normalization_fragments: int = 0
    timings: dict[str, float] = field(default_factory=dict)
    lock: dict[str, Any] = field(default_factory=dict)
    #: Interpretation modes applied on request. Empty under the default
    #: profile. Recorded so a report states which rule produced its calls.
    strict_modes: tuple[str, ...] = ()

    # -- summary views used by the report ---------------------------------- #

    @property
    def supported(self) -> list[PathogenResult]:
        """Findings the front page counts, in report order.

        Confirmed technical artifacts are excluded here by construction: the
        display precedence rule rewrites their status, so they never reach
        this list while remaining visible in their own section.
        """
        return [
            r
            for r in self.results.values()
            if r.sequence_status in SUPPORTED_STATUSES
            and r.display_status != "technical_artifact"
        ]

    @property
    def attention(self) -> list[PathogenResult]:
        """Supported findings of established, rare or unexpected pathogens."""
        return [
            r
            for r in self.supported
            if r.interpretation_class
            in {
                "established_enteric",
                "rare_enteric",
                "toxin_or_pathotype_dependent",
                "extraintestinal_watch",
            }
        ]

    @property
    def opportunists(self) -> list[PathogenResult]:
        return [
            r for r in self.supported if r.interpretation_class == "conditional_opportunist"
        ]

    @property
    def uncertain(self) -> list[PathogenResult]:
        """Candidate, marker-only, ambiguous and contamination-qualified rows."""
        return [
            r
            for r in self.results.values()
            if r.sequence_status in {"candidate_signal", "ambiguous_signal"}
            or r.display_qualifier is not None
            and r.sequence_status in SUPPORTED_STATUSES
        ]

    @property
    def pathogen_count(self) -> int:
        """The number shown beside the bug on the front page.

        Supported findings of organisms that can cause disease. Background,
        dietary and decoy organisms are deliberately excluded: counting
        harmless relatives would make the number meaningless.
        """
        return len(
            [
                r
                for r in self.supported
                if r.interpretation_class != "background_or_decoy"
            ]
        )

    def by_group(self, group: str) -> list[PathogenResult]:
        return sorted(
            (r for r in self.results.values() if r.group == group),
            key=lambda r: (
                _STATUS_RANK.get(r.sequence_status, 9),
                r.display_name,
            ),
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": "pathogens.branch.v1",
            "sample_id": self.sample_id,
            "analysis_status": self.analysis_status,
            "not_run_reason": self.not_run_reason,
            "assay": self.assay.to_json(),
            "bundle_id": self.bundle_id,
            "catalog_version": self.catalog_version,
            "calling_profile": self.calling_profile.to_json(),
            "strict_modes": list(self.strict_modes),
            "reference_lock": self.lock,
            "normalization": {
                "denominator": _DENOM,
                "fragments": self.normalization_fragments,
            },
            "coverage": None if self.coverage is None else self.coverage.to_json(),
            "candidate_stage": (
                None if self.candidates is None else self.candidates.to_json()
            ),
            "counts": {
                "pathogen_count": self.pathogen_count,
                "supported": len(self.supported),
                "attention": len(self.attention),
                "opportunists": len(self.opportunists),
                "uncertain": len(self.uncertain),
            },
            "failed_routes": list(self.failed_routes),
            "results": [
                r.to_json()
                for r in sorted(
                    self.results.values(),
                    key=lambda r: (r.group, _STATUS_RANK.get(r.sequence_status, 9),
                                   r.display_name),
                )
            ],
            "determinants": [d.to_json() for d in self.determinants],
            "timings_s": {k: round(v, 1) for k, v in self.timings.items()},
        }


_STATUS_RANK: Final[dict[str, int]] = {
    "supported_sequence": 0,
    "marker_signal": 1,
    "ambiguous_signal": 2,
    "candidate_signal": 3,
    "not_detected": 4,
    "not_assessed": 5,
}


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


def run_pathogen_branch(
    *,
    sample_id: str,
    r1: Path | None,
    r2: Path | None,
    catalog: Catalog,
    bundle_dir: Path | None,
    out_dir: Path,
    assay: AssayManifest | None = None,
    qc_passed: bool = True,
    normalization_fragments: int = 0,
    threads: int = 8,
    profile_id: str = "research_dna_v1",
    rna_profile_installed: bool = False,
    host_fasta: Path | None = None,
    determinant_dir: Path | None = None,
    progress: Callable[[str], None] | None = None,
    reuse: bool = True,
    strict_cdiff: bool = False,
) -> PathogenBranchResult:
    """Run the whole branch for one sample.

    Three ways this returns without findings, all of them explicit:
    no raw reads, no installed bundle, or a failed route. None of them
    produces a negative panel — an abundance table cannot support this
    analysis, and saying so is the correct output (acceptance P002).
    """
    say = progress or (lambda _m: None)
    out_dir.mkdir(parents=True, exist_ok=True)
    profile = load_calling_profile(profile_id)
    assay = assay or AssayManifest(sample_id=sample_id)
    timings: dict[str, float] = {}

    manifest = None
    if bundle_dir is not None and (Path(bundle_dir) / "manifest.json").is_file():
        manifest = BundleManifest.load(Path(bundle_dir) / "manifest.json")
    targets = catalog.compile_targets(
        manifest.entries if manifest is not None else None,
        default_profile_id=profile_id,
    )

    def barren(status: str, reason: str) -> PathogenBranchResult:
        """Emit a fully-populated not-assessed panel with an explicit reason."""
        results: dict[str, PathogenResult] = {}
        for target in targets:
            elig = check_assay_eligibility(
                assay, target, rna_profile_installed=rna_profile_installed
            )
            results[target.target_id] = _result_for(
                sample_id=sample_id,
                target=target,
                eligibility=elig,
                analysis_status=status,
                sequence="not_assessed",
                profile=profile,
                bundle_id=manifest.bundle_id if manifest else "none",
                reason_codes=(reason,),
                statement=reason_text,
            )
        # Determinants get gap records too. Omitting them would leave the
        # report with an empty toxin-and-resistance section, which reads as
        # "none found" — the exact inference this whole path exists to avoid.
        from openbiota.pathogens.determinants import screen_determinants

        gap_determinants = screen_determinants(
            catalog.determinants,
            sample_id=sample_id,
            ndaro=None,
            analysis_status=status,
        )
        branch = PathogenBranchResult(
            sample_id=sample_id,
            assay=assay,
            bundle_id=manifest.bundle_id if manifest else "none",
            catalog_version=catalog.version,
            calling_profile=profile,
            results=results,
            determinants=gap_determinants,
            analysis_status=status,
            not_run_reason=reason_text,
            failed_routes=(reason,),
            normalization_fragments=normalization_fragments,
            lock=manifest.lock if manifest else {},
        )
        branch.coverage = compile_coverage(
            targets, results, gap_determinants, failed_routes=(reason,)
        )
        return branch

    if r1 is None or not Path(r1).is_file():
        reason_text = (
            "Pathogen screening not run: raw reads required. An abundance table "
            "cannot support this analysis, so no target was assessed and no "
            "negative result is implied."
        )
        say(reason_text)
        return barren("not_run", "raw_reads_unavailable")

    if manifest is None:
        reason_text = (
            "Pathogen screening not run: no installed reference bundle. Every "
            "target is reported as not assessed with its reference gap; none is "
            "reported as negative."
        )
        say(reason_text)
        return barren("not_run", "reference_bundle_unavailable")

    index_path = Path(bundle_dir) / "kmer_index.npz"  # type: ignore[arg-type]
    if not index_path.is_file():
        reason_text = (
            "Pathogen screening not run: the reference bundle has no k-mer index."
        )
        say(reason_text)
        return barren("failed", "reference_index_missing")

    # ---- stage 1 ---------------------------------------------------------- #
    t0 = time.monotonic()
    stage_dir = out_dir / "pathogens"
    cache = stage_dir / "candidates.json"
    say("loading reference bundle index")
    index = KmerIndex.load(index_path)
    candidates = None
    if reuse and cache.is_file():
        try:
            prior = json.loads(cache.read_text())
            if prior.get("bundle_id") == manifest.bundle_id:
                candidates = CandidateSet(
                    reads_examined=prior["reads_examined"],
                    fragments_examined=prior["fragments_examined"],
                    bases_examined=prior["bases_examined"],
                    unclassified_reads=prior["unclassified_reads"],
                    shared_only_reads=prior["shared_only_reads"],
                    per_target_reads=dict(prior["per_target_reads"]),
                    candidate_fragments=prior["candidate_fragments"],
                    candidate_fastq=tuple(
                        Path(p) for p in prior["candidate_fastq"] if Path(p).is_file()
                    ),
                    elapsed_s=prior.get("elapsed_s", 0.0),
                )
                say("reusing cached candidate scan")
        except (OSError, KeyError, json.JSONDecodeError):
            candidates = None
    if candidates is None:
        candidates = discover_candidates(
            Path(r1), Path(r2) if r2 else None, index,
            work_dir=stage_dir, progress=say,
        )
        payload = candidates.to_json()
        payload["bundle_id"] = manifest.bundle_id
        payload["candidate_fastq"] = [str(p) for p in candidates.candidate_fastq]
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=1))
    timings["candidate scan"] = time.monotonic() - t0

    if not normalization_fragments:
        normalization_fragments = candidates.fragments_examined

    # ---- stage 2 ---------------------------------------------------------- #
    t0 = time.monotonic()
    failed_routes: list[str] = []
    evidence = AlignmentEvidence()
    if candidates.candidate_fastq:
        try:
            say(
                f"competitive alignment of {candidates.candidate_fragments:,} "
                "candidate fragments"
            )
            paths = candidates.candidate_fastq
            # Align into one index over the whole bundle, built once and
            # then reused by every sample. Restricting it to this sample's
            # candidates and their measured relatives saves little — on a
            # real library that set was already 197 of 377 references — and
            # costs a full rebuild per sample, because the member set (and
            # so the cache key) changes with every sample. The shared index
            # is also the more competitive of the two.
            prefix = build_competitive_index(
                Path(bundle_dir),  # type: ignore[arg-type]
                None,
                manifest=manifest,
                include_host=host_fasta,
                threads=threads,
                progress=say,
            )
            if prefix is None:
                raise DependencyError(
                    "bowtie2-build is required to confirm candidate findings"
                )
            lines = run_bowtie2(
                prefix,
                r1=paths[0],
                r2=paths[1] if len(paths) > 1 else None,
                threads=threads,
            )
            owner_of = {
                t.target_id: int(manifest.entries[t.target_id]["index_owner"])
                for t in targets
                if "index_owner" in manifest.entries.get(t.target_id, {})
            }
            evidence = summarise_alignments(
                list(iter_sam(lines)),
                profile=profile,
                index=index,
                owner_of=owner_of,
                genome_lengths={
                    t.target_id: t.reference_bases for t in targets
                },
                viral_targets=frozenset(
                    t.target_id
                    for t in targets
                    if t.group in {"dna_viruses", "rna_viruses", "other_viruses"}
                ),
                marker_targets=frozenset(
                    t.target_id for t in targets if t.marker_only
                ),
            )
            say(
                f"confirmed evidence for {len(evidence.by_target)} targets from "
                f"{evidence.fragments_examined:,} aligned fragments"
            )
        except (DependencyError, OpenBiotaError) as exc:
            failed_routes.append("competitive_alignment")
            say(f"competitive confirmation unavailable: {exc}")
    timings["competitive alignment"] = time.monotonic() - t0

    # ---- adjudication ----------------------------------------------------- #
    alignment_failed = "competitive_alignment" in failed_routes
    policy = SupportPolicy(
        min_distinct_fragments=profile.genome_min_distinct_fragments,
        min_regions=profile.genome_min_regions,
        min_informative_bases=profile.genome_min_informative_bases,
    )
    marker_policy = SupportPolicy(
        min_distinct_fragments=profile.marker_min_distinct_fragments,
        min_regions=profile.marker_min_windows,
        min_informative_bases=profile.marker_min_informative_bases,
    )

    results: dict[str, PathogenResult] = {}
    for target in targets:
        elig = check_assay_eligibility(
            assay, target, rna_profile_installed=rna_profile_installed
        )
        ev = evidence.by_target.get(target.target_id)
        nominated = candidates.per_target_reads.get(target.target_id, 0)

        # A route failure is not a negative. Everything the aligner could not
        # confirm becomes not-assessed with the failure recorded.
        analysis = "completed"
        if alignment_failed and nominated:
            analysis = "failed"

        normalized = NormalizedEvidence(
            assay_eligibility=elig.eligibility,
            analysis_status=analysis,
            reference_usable=target.search_ready,
            qc_pass=qc_passed,
            qualifying_fragments=ev.qualifying_fragments if ev else 0,
            nonoverlapping_informative_regions=ev.informative_regions if ev else 0,
            informative_bases_covered=ev.informative_bases_covered if ev else 0,
            resolution_specific_support=ev.resolution_specific_support if ev else False,
            only_unresolved_taxonomic_support=(
                ev.only_unresolved_taxonomic_support if ev else False
            ),
            # Post-competition only. A k-mer nomination is candidate
            # *generation* (spec §2, Bradford 2024); competitive alignment is
            # what adjudicates it. Feeding the nomination back in here would
            # let every target a stray 31-mer touched resurface as a
            # candidate, which is the competitive stage overruled by the
            # stage it exists to check — in practice ~70 spurious candidates
            # per stool library, including systemic fungi and hookworms that
            # no read actually aligned to. The nomination is not lost: it is
            # the `nominated_but_not_confirmed_competitively` reason code and
            # the candidate-stage ledger.
            any_candidate_signal=bool(ev.any_candidate_signal if ev else False),
            marker_or_organelle_only=target.marker_only,
        )
        status = sequence_status(
            normalized, marker_policy if target.marker_only else policy
        )
        results[target.target_id] = _result_for(
            sample_id=sample_id,
            target=target,
            eligibility=elig,
            analysis_status=analysis,
            sequence=status,
            profile=profile,
            bundle_id=manifest.bundle_id,
            evidence=ev,
            nominated_reads=nominated,
            normalization_fragments=normalization_fragments,
            reason_codes=tuple(ev.reason_codes) if ev else (),
        )

    # ---- determinants ----------------------------------------------------- #
    # Screened against their own small reference bundle rather than the
    # organism bundle: these are short genes, and the question they answer
    # (is this gene here, and what carries it) is not the question the
    # organism k-mer index is built for.
    t0 = time.monotonic()
    determinant_results = _screen_determinants(
        catalog,
        sample_id=sample_id,
        r1=Path(r1),
        r2=Path(r2) if r2 else None,
        determinant_dir=determinant_dir,
        qc_passed=qc_passed,
        threads=threads,
        say=say,
    )
    timings["determinants"] = time.monotonic() - t0

    # ---- species resolution and pathotype gating -------------------------- #
    # Runs with both screens in hand: a species inside an indistinguishable
    # complex needs its discriminating marker, and an organism whose disease
    # potential is a toxin needs the toxin gene, before either is called a
    # pathogen. Nothing is deleted; rows are re-stated and re-tiered.
    from openbiota.pathogens.resolve import resolve as _resolve_calls

    results = _resolve_calls(
        results,
        determinant_results,
        families={t.target_id: getattr(t, "family", "") for t in targets},
        requires={
            t.target_id: bool(getattr(t, "requires_determinants", False)) for t in targets
        },
        strict_cdiff=strict_cdiff,
    )

    branch = PathogenBranchResult(
        strict_modes=("strict_cdiff",) if strict_cdiff else (),
        sample_id=sample_id,
        assay=assay,
        bundle_id=manifest.bundle_id,
        catalog_version=catalog.version,
        calling_profile=profile,
        results=results,
        determinants=determinant_results,
        alignment=evidence,
        candidates=candidates,
        analysis_status="completed",
        failed_routes=tuple(failed_routes),
        normalization_fragments=normalization_fragments,
        timings=timings,
        lock=manifest.lock,
    )
    branch.coverage = compile_coverage(
        targets, results, branch.determinants, failed_routes=failed_routes
    )
    say(
        f"pathogen branch: {branch.pathogen_count} supported finding(s), "
        f"{branch.coverage.assessed} targets assessed, "
        f"{branch.coverage.not_assessed} not assessed"
    )
    return branch


def _screen_determinants(
    catalog: Catalog,
    *,
    sample_id: str,
    r1: Path,
    r2: Path | None,
    determinant_dir: Path | None,
    qc_passed: bool,
    threads: int,
    say: Callable[[str], None],
) -> tuple[DeterminantResult, ...]:
    """Run the determinant screen, degrading to explicit gaps on any failure.

    Every exit path returns one record per seed. There is no path that
    returns nothing, because an empty determinant list in the report would
    read as "no resistance or toxin genes found" when the truth may be
    "the gene reference is not installed".
    """
    from openbiota.pathogens.determinants import (
        DeterminantBundle,
        load_ndaro_map,
        screen_determinants,
        screen_determinants_from_reads,
    )

    seeds = catalog.determinants
    if not seeds:
        return ()

    map_path = Path("pathogens/catalog/ndaro_map.yaml")
    manifest_path = (
        Path(determinant_dir) / "manifest.json" if determinant_dir else None
    )
    if not map_path.is_file() or manifest_path is None or not manifest_path.is_file():
        say("determinant reference not installed; recording gaps")
        return screen_determinants(
            seeds, sample_id=sample_id, ndaro=None,
            analysis_status="not_run", progress=say,
        )

    try:
        ndaro = load_ndaro_map(map_path)
        doc = json.loads(manifest_path.read_text())
        prefix = doc.get("index_prefix")
        bundle = DeterminantBundle(
            bundle_id=str(doc["bundle_id"]),
            catalog_version=str(doc.get("catalog_version", "unknown")),
            fasta=Path(doc["fasta"]),
            index_prefix=Path(prefix) if prefix else None,
            families_by_determinant={
                k: tuple(v) for k, v in (doc.get("families_by_determinant") or {}).items()
            },
            family_bases={
                k: int(v) for k, v in (doc.get("family_bases") or {}).items()
            },
            determinant_of_family=dict(doc.get("determinant_of_family") or {}),
            gaps=dict(doc.get("gaps") or {}),
        )
    except (OSError, KeyError, json.JSONDecodeError, OpenBiotaError) as exc:
        say(f"determinant reference unusable ({exc}); recording gaps")
        return screen_determinants(
            seeds, sample_id=sample_id, ndaro=None,
            analysis_status="failed", progress=say,
        )

    return screen_determinants_from_reads(
        seeds,
        sample_id=sample_id,
        r1=r1,
        r2=r2,
        ndaro=ndaro,
        bundle=bundle,
        qc_pass=qc_passed,
        threads=threads,
        progress=say,
    )


def _result_for(
    *,
    sample_id: str,
    target: TargetRecord,
    eligibility: EligibilityRecord,
    analysis_status: str,
    sequence: str,
    profile: CallingProfile,
    bundle_id: str,
    evidence: TargetEvidence | None = None,
    nominated_reads: int = 0,
    normalization_fragments: int = 0,
    reason_codes: Sequence[str] = (),
    statement: str | None = None,
) -> PathogenResult:
    """Assemble one result record from adjudicated evidence."""
    # Historical samples have no blanks. That is `controls_unavailable`, which
    # is emphatically not "contamination checks passed".
    contamination = "controls_unavailable"
    display, qualifier = display_precedence(
        sequence=sequence,
        contamination=contamination if sequence in SUPPORTED_STATUSES else "not_evaluated",
        reference_gap_reason=(
            target.reference_gap_note
            or (None if target.reference_usable else target.reference_status)
        ),
    )
    codes = list(dict.fromkeys(reason_codes))
    # How many reads the k-mer stage nominated is the audit trail for stage 2:
    # a target that was nominated and then rejected by competitive alignment
    # is a different story from one no read ever mentioned, and only this
    # code distinguishes them after the fact.
    if nominated_reads and (evidence is None or not evidence.qualifying_fragments):
        codes.append("nominated_but_not_confirmed_competitively")
    if not target.reference_usable:
        codes.append(f"reference_{target.reference_status}")
    if eligibility.eligibility != "eligible":
        codes.append(eligibility.reason_code)
    if sequence in SUPPORTED_STATUSES:
        codes.append("controls_unavailable")
    if target.taxonomic_resolution in {"species_complex", "group", "genus"}:
        codes.append("species_resolution_limited")

    breadth = None
    informative_breadth = None
    if evidence is not None and target.reference_bases:
        breadth = round(
            evidence.informative_bases_covered / target.reference_bases, 8
        )
    if evidence is not None and target.informative_bases:
        informative_breadth = round(
            min(1.0, evidence.informative_bases_covered / target.informative_bases), 8
        )
    normalized_fpm = None
    if evidence is not None and normalization_fragments:
        normalized_fpm = round(
            evidence.qualifying_fragments / normalization_fragments * 1e6, 4
        )

    # One word for what the reads were searched against, decided here so the
    # report never has to infer coverage from a status code. A row whose
    # sequence is installed and aligned but shares every informative k-mer
    # with a catalogue relative was searched: it answers at the shared rank.
    if eligibility.eligibility != "eligible" and "wrong_nucleic_acid" in codes:
        search_scope = "out_of_assay_scope"
    elif target.reference_usable:
        search_scope = "own_rank"
    elif target.reference_status == "covered_by_relative":
        search_scope = "group_rank"
    elif target.reference_bases and not target.informative_bases:
        # Sequence installed and aligned against, but nothing in it is this
        # target's alone. Either a catalogue relative has the same genome, or
        # the host does: GRCh38 carries Epstein-Barr virus as a decoy contig,
        # so an EBV read cannot be told from a read of the reader.
        search_scope = "group_rank" if target.near_neighbor_target_ids else "masked_by_host"
    else:
        search_scope = "reference_pending"

    return PathogenResult(
        schema_version="pathogens.result.v1",
        sample_id=sample_id,
        target_id=target.target_id,
        search_scope=search_scope,
        group=target.group,
        display_name=target.display_name,
        interpretation_class=target.interpretation_class,
        stool_role=target.stool_role,
        assay_eligibility=eligibility.eligibility,
        analysis_status=analysis_status,
        reference_status=target.reference_status,
        sequence_status=sequence,
        contamination_status=contamination,
        resolution=target.taxonomic_resolution,
        covered_by_target_id=target.covered_by_target_id,
        clinical_interpretation="organism_sequence_not_infection_diagnosis",
        validation_scope="computational_research_rule",
        calling_profile_id=profile.profile_id,
        reference_bundle_id=bundle_id,
        unique_supporting_fragments=evidence.qualifying_fragments if evidence else 0,
        ambiguous_fragments=evidence.ambiguous_fragments if evidence else 0,
        informative_regions_supported=evidence.informative_regions if evidence else 0,
        informative_bases_covered=(
            evidence.informative_bases_covered if evidence else 0
        ),
        reference_breadth_fraction=breadth,
        informative_region_breadth_fraction=informative_breadth,
        median_alignment_identity=evidence.median_identity if evidence else None,
        normalized_fragments_per_million=normalized_fpm,
        normalization_denominator=_DENOM if normalized_fpm is not None else None,
        reason_codes=tuple(dict.fromkeys(codes)),
        confirmation_options=target.confirmation_options,
        display_status=display,
        display_qualifier=qualifier,
        plain_statement=statement
        or interpret_finding(
            target,
            sequence,
            eligibility,
            marker_only=target.marker_only,
            ambiguous_with=evidence.best_competing_target if evidence else None,
            contamination=contamination,
        ),
    )


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #


def self_test() -> int:
    from openbiota.pathogens.eligibility import _target  # noqa: PLC2701 - test helper

    checks = 0
    dna = AssayManifest(sample_id="s", nucleic_acid_protocol="DNA")

    # A not-detected established pathogen must state its limitations and must
    # not reassure.
    strong = _target(
        target_id="protozoa.giardia_duodenalis",
        display_name="Giardia duodenalis",
        group="protozoa",
        negative_limitation_codes=("sampling", "shedding", "extraction", "depth"),
    )
    elig = check_assay_eligibility(dna, strong)
    text = interpret_finding(strong, "not_detected", elig)
    assert "not that the organism is absent" in text
    assert "intermittent" in text or "low or intermittent" in text or "shedding" in text
    for banned in ("all clear", "clear of", "infection-free", "you are healthy"):
        assert banned not in text.lower()
    checks += 1

    # P077: a Strongyloides non-detection carries sampling and shedding limits.
    strongy = _target(
        target_id="helminths.strongyloides_stercoralis",
        display_name="Strongyloides stercoralis",
        group="helminths",
        stool_role="variable_shedding",
        negative_limitation_codes=("sampling", "shedding", "depth"),
        confirmation_options=("strongyloides_specialist_stool_and_serology",),
    )
    text = interpret_finding(strongy, "not_detected", check_assay_eligibility(dna, strongy))
    assert "shedding can be low or intermittent" in text
    checks += 1

    # P080: a tissue-restricted target's non-detection is not an exclusion.
    tissue = _target(
        target_id="helminths.echinococcus_granulosus",
        display_name="Echinococcus granulosus",
        group="helminths",
        stool_role="tissue_restricted",
        interpretation_class="extraintestinal_watch",
        negative_limitation_codes=("tissue_restricted", "no_patent_stool_stage"),
    )
    text = interpret_finding(tissue, "not_detected", check_assay_eligibility(dna, tissue))
    assert "does not screen out or exclude" in text
    checks += 1

    # A supported Candida finding names what it does not establish.
    candida = _target(
        target_id="fungi.candida_albicans",
        display_name="Candida albicans",
        group="fungi",
        interpretation_class="conditional_opportunist",
        stool_role="carriage_common",
    )
    text = interpret_finding(
        candida, "supported_sequence", check_assay_eligibility(dna, candida)
    )
    assert "does not diagnose an infection" in text
    assert "need for treatment" in text
    checks += 1

    # A marker-only finding says so explicitly.
    text = interpret_finding(
        strongy,
        "marker_signal",
        check_assay_eligibility(dna, strongy),
        marker_only=True,
    )
    assert "marker-level evidence" in text and "not as whole-genome support" in text
    checks += 1

    # Ambiguity reports the competing relative rather than picking one.
    text = interpret_finding(
        strong,
        "ambiguous_signal",
        check_assay_eligibility(dna, strong),
        ambiguous_with="protozoa.entamoeba_dispar",
    )
    assert "unresolved" in text and "entamoeba dispar" in text.lower()
    checks += 1

    # A candidate signal is neither a finding nor an absence.
    text = interpret_finding(
        strong, "candidate_signal", check_assay_eligibility(dna, strong)
    )
    assert "not as a finding" in text and "not evidence of absence" in text
    checks += 1

    # An RNA target in a DNA library reproduces the eligibility wording.
    rna = _target(
        target_id="rna_viruses.norovirus_group",
        display_name="Human noroviruses",
        group="rna_viruses",
        allowed_nucleic_acids=("RNA",),
        host_category="human",
        host_evidence_level="established_human_pathogen",
        molecule_type="ssRNA_positive",
        dna_only_statement=(
            "Not assessed: this library was prepared for DNA, while norovirus "
            "has an RNA genome."
        ),
    )
    text = interpret_finding(rna, "not_assessed", check_assay_eligibility(dna, rna))
    assert text.startswith("Not assessed:") and "RNA genome" in text
    checks += 1

    # Coverage arithmetic keeps organisms and determinants apart.
    ledger = compile_coverage([strong, rna], {}, [])
    assert ledger.total_targets == 2 and ledger.not_assessed == 2
    assert ledger.determinants_assessed == 0
    checks += 1

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} detection checks passed")
