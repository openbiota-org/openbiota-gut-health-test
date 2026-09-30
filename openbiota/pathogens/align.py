"""Competitive alignment and evidence extraction (BUILD_SPEC_v05.0 §5.3).

Candidate generation asks "is there anything here?". This module asks the far
harder question: *is this sequence really this target, rather than a relative,
the host, dinner, a repeat, or a duplicate of one molecule counted five
times?* Everything here exists to answer that honestly.

The evidence rules, and why each one is needed:

**Fragments, not reads.** Overlapping mates are one fragment. Different read
names do not establish independent support, so non-UMI libraries collapse
conservative duplicate families — same canonical template start, end and
orientation — to one fragment for the support gate (acceptance P031, P032).
Raw and deduplicated counts are both retained, because coordinate
deduplication is not exact molecule counting and can merge genuine molecules.

**Target-level uniqueness, not raw MAPQ.** A read that maps equally well to
five strains of the same target is unique *at target level*; bowtie2 reports
MAPQ 0 or 1 for it. Rejecting it would penalise well-sampled targets for
being well-sampled (acceptance P027). Conversely, a read that maps equally
well to another target is ambiguous no matter how high its MAPQ, and yields
group-level evidence instead of a species call (acceptance P028).

**Independent regions, not adjacent windows.** Regions are fixed canonical
10-kb bins of the reference. Three adjacent windows inside one locus collapse
to one region (acceptance P030). Small viral genomes use nonoverlapping
thirds instead, so a 7-kb genome can still demonstrate distributed support.

**Region specificity.** A qualifying genomic fragment must land in a bin the
reference build proved carries sequence absent from every relative, the host
and the decoys. Five reads in a conserved rDNA locus therefore cannot
establish a species; they stay marker-like ambiguous evidence (P029).
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import DependencyError
from openbiota.pathogens.catalog import CallingProfile
from openbiota.pathogens.kmers import KmerIndex

__all__ = [
    "AlignmentEvidence",
    "SamRecord",
    "TargetEvidence",
    "iter_sam",
    "run_bowtie2",
    "summarise_alignments",
    "self_test",
]

_CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")

#: Reads whose best alignment beats the best *other-target* alignment by less
#: than the profile's margin are ambiguous. Reported at the rank the data
#: supports rather than forced onto one species.
_LOW_COMPLEXITY_MAX_ENTROPY: Final[float] = 1.2


@dataclass(frozen=True)
class SamRecord:
    """One alignment line, reduced to the fields the evidence rules use."""

    query: str
    flag: int
    target_id: str
    contig: str
    pos: int
    mapq: int
    query_length: int
    aligned_bases: int
    edit_distance: int
    score: int
    sequence: bytes

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def is_read2(self) -> bool:
        return bool(self.flag & 0x80)

    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & 0x100)

    @property
    def identity(self) -> float:
        """Aligned identity over aligned bases."""
        if self.aligned_bases <= 0:
            return 0.0
        return max(0.0, 1.0 - self.edit_distance / self.aligned_bases)

    @property
    def aligned_fraction(self) -> float:
        if self.query_length <= 0:
            return 0.0
        return min(1.0, self.aligned_bases / self.query_length)

    @property
    def end(self) -> int:
        return self.pos + self.aligned_bases


def _cigar_aligned(cigar: bytes) -> tuple[int, int]:
    """Return ``(reference_span, aligned_query_bases)`` for a CIGAR string."""
    span = 0
    aligned = 0
    for count, op in _CIGAR.findall(cigar):
        n = int(count)
        if op in b"M=X":
            span += n
            aligned += n
        elif op in b"DN":
            span += n
        elif op == b"I"[0:1] or op == b"I":
            aligned += n
    return span, aligned


def iter_sam(lines: Iterable[bytes]) -> Iterator[SamRecord]:
    """Parse SAM text into `SamRecord`s.

    We parse SAM directly rather than requiring samtools: the competitive
    stage aligns only candidate reads, so the volume is small and one fewer
    external binary is one fewer thing to pin and version.

    Reference names are `target_id|contig` as written by the bundle build, so
    the taxonomic identity of every alignment comes from the frozen bundle
    rather than from a parsed organism name.
    """
    for line in lines:
        if not line or line.startswith(b"@"):
            continue
        parts = line.rstrip(b"\n").split(b"\t")
        if len(parts) < 11:
            continue
        flag = int(parts[1])
        if flag & 0x4:  # unmapped
            continue
        rname = parts[2].decode("utf-8", "replace")
        target_id, _, contig = rname.partition("|")
        cigar = parts[5]
        seq = parts[9]
        _span, aligned = _cigar_aligned(cigar)
        nm = 0
        score = 0
        for field_ in parts[11:]:
            if field_.startswith(b"NM:i:"):
                nm = int(field_[5:])
            elif field_.startswith(b"AS:i:"):
                score = int(field_[5:])
        query_length = len(seq) if seq != b"*" else aligned
        yield SamRecord(
            query=parts[0].decode("utf-8", "replace"),
            flag=flag,
            target_id=target_id,
            contig=contig or rname,
            pos=int(parts[3]) - 1,
            mapq=int(parts[4]),
            query_length=query_length,
            aligned_bases=aligned,
            edit_distance=nm,
            score=score,
            sequence=b"" if seq == b"*" else seq,
        )


def run_bowtie2(
    index_prefix: Path,
    *,
    r1: Path,
    r2: Path | None = None,
    threads: int = 8,
    binary: str = "bowtie2",
    max_alignments: int = 20,
    timeout: float = 7200.0,
) -> list[bytes]:
    """Align candidate reads against the frozen competitive bundle.

    `-k` reporting is essential: we need the *alternatives* to a read's best
    placement in order to judge ambiguity. A single best-hit alignment against
    a pathogen-only reference is exactly the setup that produces confident
    false assignments.
    """
    exe = shutil.which(binary)
    if exe is None:
        raise DependencyError(
            f"{binary} is not installed; competitive confirmation cannot run and "
            "candidate findings therefore cannot be promoted to supported findings"
        )
    argv = [
        exe,
        "--threads", str(threads),
        "-x", str(index_prefix),
        "-k", str(max_alignments),
        "--local",
        "--no-unal",
        "--reorder",
        "--seed", "20260907",
    ]
    if r2 is not None:
        argv += ["-1", str(r1), "-2", str(r2)]
    else:
        argv += ["-U", str(r1)]
    proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
        argv, capture_output=True, timeout=timeout, check=False
    )
    if proc.returncode != 0:
        tail = proc.stderr.decode("utf-8", "replace").strip().splitlines()[-3:]
        raise DependencyError(
            f"bowtie2 failed (rc={proc.returncode}): " + " | ".join(tail)
        )
    return proc.stdout.splitlines()


# --------------------------------------------------------------------------- #
# Evidence
# --------------------------------------------------------------------------- #


@dataclass
class TargetEvidence:
    """Normalised alignment evidence for one target in one sample.

    Field names deliberately mirror `NormalizedEvidence` in the kernel: this
    is the adapter boundary, and keeping the vocabulary identical makes it
    obvious that no extra massaging happens in between.
    """

    target_id: str
    qualifying_fragments: int = 0
    raw_fragments: int = 0
    duplicate_like_fragments: int = 0
    ambiguous_fragments: int = 0
    informative_regions: int = 0
    informative_bases_covered: int = 0
    resolution_specific_support: bool = False
    only_unresolved_taxonomic_support: bool = False
    any_candidate_signal: bool = False
    median_identity: float | None = None
    max_score: int = 0
    shared_support_fraction: float | None = None
    best_competing_target: str | None = None
    support_shape: tuple[str, ...] = ()
    regions: tuple[int, ...] = ()
    reason_codes: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "qualifying_fragments": self.qualifying_fragments,
            "raw_fragments": self.raw_fragments,
            "duplicate_like_fragments": self.duplicate_like_fragments,
            "ambiguous_fragments": self.ambiguous_fragments,
            "informative_regions": self.informative_regions,
            "informative_bases_covered": self.informative_bases_covered,
            "resolution_specific_support": self.resolution_specific_support,
            "only_unresolved_taxonomic_support": self.only_unresolved_taxonomic_support,
            "any_candidate_signal": self.any_candidate_signal,
            "median_identity": self.median_identity,
            "shared_support_fraction": self.shared_support_fraction,
            "best_competing_target": self.best_competing_target,
            "support_shape": list(self.support_shape),
            "reason_codes": list(self.reason_codes),
        }


@dataclass
class AlignmentEvidence:
    """All targets' evidence from one competitive alignment run."""

    by_target: dict[str, TargetEvidence] = field(default_factory=dict)
    fragments_examined: int = 0
    fragments_aligned: int = 0

    def to_json(self) -> dict[str, Any]:
        return {
            "fragments_examined": self.fragments_examined,
            "fragments_aligned": self.fragments_aligned,
            "by_target": {k: v.to_json() for k, v in sorted(self.by_target.items())},
        }



def _add_reason(existing: Sequence[str], code: str) -> tuple[str, ...]:
    """Append a reason code once, preserving order.

    Built with a set literal until v0.8.3, which meant the same codes
    serialised in a different order on every run: Python randomises string
    hashing per process, so `tuple({*codes, "x"})` is not reproducible. The
    content never changed, but a preservation diff of two identical runs
    reported 129 altered values, and a reader comparing two reports would
    have seen the reasons shuffle for no reason.
    """
    out = list(dict.fromkeys(existing))
    if code not in out:
        out.append(code)
    return tuple(out)


def _entropy(sequence: bytes) -> float:
    """Shannon entropy over 2-mers, used to reject low-complexity support."""
    if len(sequence) < 8:
        return 4.0
    counts: dict[bytes, int] = {}
    for i in range(len(sequence) - 1):
        pair = sequence[i : i + 2]
        counts[pair] = counts.get(pair, 0) + 1
    total = sum(counts.values())
    from math import log2

    return -sum((n / total) * log2(n / total) for n in counts.values())


def _region_of(
    record: SamRecord,
    *,
    bin_bases: int,
    genome_length: int,
    viral: bool,
) -> int:
    """Canonical independent-region id for one alignment.

    Fixed bins, not alignment clusters: this is what makes three adjacent
    windows of one locus count once. Small viral genomes are divided into
    thirds instead so that distributed support remains demonstrable.
    """
    if viral and genome_length > 0:
        third = max(1, genome_length // 3)
        return min(2, record.pos // third)
    return record.pos // bin_bases


def summarise_alignments(
    records: Sequence[SamRecord],
    *,
    profile: CallingProfile,
    index: KmerIndex | None = None,
    owner_of: Mapping[str, int] | None = None,
    genome_lengths: Mapping[str, int] | None = None,
    viral_targets: frozenset[str] = frozenset(),
    marker_targets: frozenset[str] = frozenset(),
) -> AlignmentEvidence:
    """Turn raw alignments into per-target normalised evidence.

    Processing is per fragment, because every rule that prevents
    double-counting — mate overlap, duplicate families, ambiguity between
    targets — is a property of the fragment rather than of one alignment line.
    """
    owner_of = owner_of or {}
    genome_lengths = genome_lengths or {}

    # ---- group alignments by fragment ------------------------------------- #
    by_fragment: dict[str, list[SamRecord]] = {}
    for record in records:
        by_fragment.setdefault(record.query, []).append(record)

    evidence: dict[str, TargetEvidence] = {}

    def slot(target_id: str) -> TargetEvidence:
        if target_id not in evidence:
            evidence[target_id] = TargetEvidence(target_id=target_id)
        return evidence[target_id]

    # Per target: duplicate-family keys already counted, region ids, covered
    # intervals inside informative regions, and identity samples.
    families: dict[str, set[tuple[int, int, int, bool]]] = {}
    regions: dict[str, set[int]] = {}
    intervals: dict[str, list[tuple[int, int]]] = {}
    identities: dict[str, list[float]] = {}
    shared_counts: dict[str, int] = {}

    for _query, group in by_fragment.items():
        # Best alignment per target within this fragment. Mates and secondary
        # placements of the same target collapse here, which is the mate-overlap
        # and conspecific-multimapping rule in one step.
        best_by_target: dict[str, SamRecord] = {}
        for record in group:
            prior = best_by_target.get(record.target_id)
            if prior is None or record.score > prior.score:
                best_by_target[record.target_id] = record
        if not best_by_target:
            continue

        ranked = sorted(best_by_target.values(), key=lambda r: -r.score)
        winner = ranked[0]
        runner_up = ranked[1] if len(ranked) > 1 else None
        margin = (
            winner.score - runner_up.score
            if runner_up is not None
            else profile.min_best_minus_other_taxon_alignment_score
        )

        # Which targets this fragment is evidence *for*. Only the winner and
        # anything within the ambiguity margin of it: a target the fragment
        # aligned to but lost to by more than that margin has been decided
        # against, and the read belongs to the winner. Crediting every target
        # a fragment merely touched is how conserved and repetitive sequence
        # turns into a page of systemic fungi and hookworms in an ordinary
        # stool library — the competitive stage overruled by its own inputs.
        contenders = [
            r.target_id
            for r in ranked
            if winner.score - r.score
            < profile.min_best_minus_other_taxon_alignment_score
        ] or [winner.target_id]
        for target_id in contenders:
            slot(target_id).any_candidate_signal = True
            slot(target_id).raw_fragments += 1

        # Ambiguity between *different* targets: report at the supported rank,
        # never forced onto the better-scoring species.
        if runner_up is not None and margin < profile.min_best_minus_other_taxon_alignment_score:
            tied = list(contenders)
            for target_id in tied:
                ev = slot(target_id)
                ev.ambiguous_fragments += 1
                ev.only_unresolved_taxonomic_support = True
                shared_counts[target_id] = shared_counts.get(target_id, 0) + 1
                others = [t for t in tied if t != target_id]
                if others and ev.best_competing_target is None:
                    ev.best_competing_target = others[0]
            continue

        ev = slot(winner.target_id)
        # Unambiguous at target level: this fragment can resolve the target,
        # so it no longer forces group-level reporting.
        ev.only_unresolved_taxonomic_support = False
        if runner_up is not None and ev.best_competing_target is None:
            ev.best_competing_target = runner_up.target_id

        # ---- alignment quality gates ------------------------------------- #
        if winner.identity < profile.min_identity:
            ev.reason_codes = _add_reason(ev.reason_codes, "below_identity_gate")
            continue
        if winner.aligned_fraction < profile.min_aligned_query_fraction:
            ev.reason_codes = _add_reason(ev.reason_codes, "partial_query_alignment")
            continue
        if winner.aligned_bases < profile.min_aligned_bases:
            ev.reason_codes = _add_reason(ev.reason_codes, "short_alignment")
            continue
        if (
            profile.exclude_low_complexity_or_masked_only_support
            and winner.sequence
            and _entropy(winner.sequence) < _LOW_COMPLEXITY_MAX_ENTROPY
        ):
            ev.reason_codes = _add_reason(ev.reason_codes, "low_complexity_support")
            continue

        # ---- duplicate families ------------------------------------------- #
        # Conservative: identical canonical template coordinates and
        # orientation collapse to one supporting fragment regardless of read
        # name. Raw counts are kept alongside.
        key = (
            hash(winner.contig) & 0xFFFFFFFF,
            winner.pos,
            winner.end,
            winner.is_reverse,
        )
        seen = families.setdefault(winner.target_id, set())
        if key in seen:
            ev.duplicate_like_fragments += 1
            ev.reason_codes = _add_reason(ev.reason_codes, "duplicate_family_collapsed")
            continue
        seen.add(key)

        # ---- region accounting -------------------------------------------- #
        region = _region_of(
            winner,
            bin_bases=profile.region_bin_bases,
            genome_length=genome_lengths.get(winner.target_id, 0),
            viral=winner.target_id in viral_targets,
        )
        specific = True
        if index is not None and winner.target_id in owner_of:
            specific = index.is_informative_bin(owner_of[winner.target_id], region)
        if not specific:
            # Support in a region indistinguishable from a relative: real
            # sequence, but not species-resolving. Keep it as candidate-level
            # evidence rather than counting it toward the species gate.
            ev.ambiguous_fragments += 1
            ev.reason_codes = _add_reason(ev.reason_codes, "support_in_shared_region_only")
            continue

        ev.qualifying_fragments += 1
        ev.resolution_specific_support = True
        regions.setdefault(winner.target_id, set()).add(region)
        intervals.setdefault(winner.target_id, []).append((winner.pos, winner.end))
        identities.setdefault(winner.target_id, []).append(winner.identity)
        ev.max_score = max(ev.max_score, winner.score)

    # ---- finalise per target ---------------------------------------------- #
    for target_id, ev in evidence.items():
        ev.informative_regions = len(regions.get(target_id, ()))
        ev.regions = tuple(sorted(regions.get(target_id, ())))
        ev.informative_bases_covered = _merged_length(intervals.get(target_id, []))
        vals = sorted(identities.get(target_id, []))
        if vals:
            mid = len(vals) // 2
            ev.median_identity = (
                vals[mid] if len(vals) % 2 else (vals[mid - 1] + vals[mid]) / 2
            )
        total = ev.qualifying_fragments + ev.ambiguous_fragments
        if total:
            ev.shared_support_fraction = ev.ambiguous_fragments / total
        shapes: list[str] = []
        if target_id in marker_targets:
            shapes.append("marker_or_organelle_reference_only")
        if ev.informative_regions == 1 and ev.qualifying_fragments >= 3:
            shapes.append("single_locus_support")
        if ev.duplicate_like_fragments > ev.qualifying_fragments:
            shapes.append("duplicate_dominated")
        if ev.ambiguous_fragments and not ev.qualifying_fragments:
            shapes.append("shared_regions_only")
        ev.support_shape = tuple(shapes)

    return AlignmentEvidence(
        by_target=evidence,
        fragments_examined=len(by_fragment),
        fragments_aligned=sum(1 for g in by_fragment.values() if g),
    )


def _merged_length(intervals: Sequence[tuple[int, int]]) -> int:
    """Total length of the union of aligned intervals.

    Overlapping alignments must not each contribute their full length to
    "informative bases covered": twenty reads stacked on one 150 bp site cover
    150 bases, not 3,000.
    """
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    start, end = ordered[0]
    for lo, hi in ordered[1:]:
        if lo > end:
            total += end - start
            start, end = lo, hi
        else:
            end = max(end, hi)
    total += end - start
    return total


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #


def _rec(
    query: str,
    target: str,
    pos: int,
    *,
    score: int = 300,
    nm: int = 0,
    aligned: int = 150,
    qlen: int = 150,
    reverse: bool = False,
    seq: bytes | None = None,
    read2: bool = False,
) -> SamRecord:
    flag = (0x10 if reverse else 0) | (0x80 if read2 else 0x40)
    return SamRecord(
        query=query,
        flag=flag,
        target_id=target,
        contig="c1",
        pos=pos,
        mapq=42,
        query_length=qlen,
        aligned_bases=aligned,
        edit_distance=nm,
        score=score,
        sequence=seq if seq is not None else b"ACGT" * 40,
    )


def self_test() -> int:
    import numpy as np

    from openbiota.pathogens.catalog import RESEARCH_DNA_V1 as P

    checks = 0
    # An index where target 0's bins 0..9 are all informative.
    index = KmerIndex(
        kmers=np.array([1], dtype=np.uint64),
        owners=np.array([0], dtype=np.uint16),
        target_ids=("t",),
        informative_bins={0: np.arange(10, dtype=np.int64)},
    )
    owner_of = {"t": 0}

    # P034: five fragments in three distinct bins with 300+ informative bases.
    records = [
        _rec(f"f{i}", "t", i * 10_000 + 100) for i in range(5)
    ]
    ev = summarise_alignments(
        records, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.qualifying_fragments == 5, ev
    assert ev.informative_regions == 5
    assert ev.informative_bases_covered == 750
    assert ev.resolution_specific_support
    checks += 1

    # P032: overlapping mates of one fragment count once.
    mates = [
        _rec("frag", "t", 100),
        _rec("frag", "t", 120, read2=True, score=290),
    ]
    ev = summarise_alignments(
        mates, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.qualifying_fragments == 1, ev
    checks += 1

    # P031: PCR duplicates with different names cannot satisfy the gate.
    dups = [_rec(f"dup{i}", "t", 500) for i in range(8)]
    ev = summarise_alignments(
        dups, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.qualifying_fragments == 1 and ev.duplicate_like_fragments == 7, ev
    assert "duplicate_dominated" in ev.support_shape
    checks += 1

    # P030: three adjacent windows of one locus are one region.
    adjacent = [_rec(f"a{i}", "t", 1000 + i * 200) for i in range(3)]
    ev = summarise_alignments(
        adjacent, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.qualifying_fragments == 3 and ev.informative_regions == 1, ev
    assert "single_locus_support" in ev.support_shape
    checks += 1

    # P029: support only in a shared (non-informative) bin is not species-level.
    shared_only = [_rec(f"s{i}", "t", 200_000 + i * 300) for i in range(20)]
    ev = summarise_alignments(
        shared_only, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.qualifying_fragments == 0 and ev.ambiguous_fragments == 20, ev
    assert not ev.resolution_specific_support
    assert "shared_regions_only" in ev.support_shape
    checks += 1

    # P028: equally good placements on another target give ambiguity, not a
    # forced species call.
    two_index = KmerIndex(
        kmers=np.array([1], dtype=np.uint64),
        owners=np.array([0], dtype=np.uint16),
        target_ids=("t", "u"),
        informative_bins={0: np.arange(10, dtype=np.int64),
                          1: np.arange(10, dtype=np.int64)},
    )
    tie = []
    for i in range(6):
        tie.append(_rec(f"x{i}", "t", i * 10_000, score=300))
        tie.append(_rec(f"x{i}", "u", i * 10_000, score=300))
    out = summarise_alignments(
        tie, profile=P, index=two_index, owner_of={"t": 0, "u": 1}
    )
    for name in ("t", "u"):
        assert out.by_target[name].qualifying_fragments == 0
        assert out.by_target[name].ambiguous_fragments == 6
        assert out.by_target[name].only_unresolved_taxonomic_support
    checks += 1

    # A clear winner over a relative resolves normally.
    clear = []
    for i in range(6):
        clear.append(_rec(f"y{i}", "t", i * 10_000, score=300))
        clear.append(_rec(f"y{i}", "u", i * 10_000, score=200))
    out = summarise_alignments(
        clear, profile=P, index=two_index, owner_of={"t": 0, "u": 1}
    )
    assert out.by_target["t"].qualifying_fragments == 6
    assert out.by_target["t"].any_candidate_signal
    # The relative lost by more than the ambiguity margin, so the fragments
    # are the winner's and the loser is not a candidate on their strength.
    # Anything else and every conserved region in the bundle would nominate
    # its whole neighbourhood.
    assert "u" not in out.by_target
    checks += 1

    # Identity, aligned-fraction and length gates each reject on their own.
    assert summarise_alignments(
        [_rec("z", "t", 10, nm=30)], profile=P, index=index, owner_of=owner_of
    ).by_target["t"].qualifying_fragments == 0
    assert summarise_alignments(
        [_rec("z", "t", 10, aligned=100, qlen=150)],
        profile=P, index=index, owner_of=owner_of,
    ).by_target["t"].qualifying_fragments == 0
    assert summarise_alignments(
        [_rec("z", "t", 10, aligned=50, qlen=50)],
        profile=P, index=index, owner_of=owner_of,
    ).by_target["t"].qualifying_fragments == 0
    checks += 3

    # Low-complexity support is excluded.
    assert summarise_alignments(
        [_rec("z", "t", 10, seq=b"A" * 150)],
        profile=P, index=index, owner_of=owner_of,
    ).by_target["t"].qualifying_fragments == 0
    checks += 1

    # Viral thirds: a 9-kb genome gives three regions, not one bin.
    viral = [_rec(f"v{i}", "t", i * 3_000, aligned=150) for i in range(3)]
    ev = summarise_alignments(
        viral,
        profile=P,
        index=index,
        owner_of=owner_of,
        genome_lengths={"t": 9_000},
        viral_targets=frozenset({"t"}),
    ).by_target["t"]
    assert ev.informative_regions == 3, ev
    checks += 1

    # Overlapping alignments do not inflate covered bases.
    stacked = [_rec(f"o{i}", "t", 100 + i, aligned=150) for i in range(5)]
    ev = summarise_alignments(
        stacked, profile=P, index=index, owner_of=owner_of
    ).by_target["t"]
    assert ev.informative_bases_covered == 154, ev.informative_bases_covered
    checks += 1

    # CIGAR arithmetic and SAM parsing.
    line = (
        b"r1\t99\tbacteria.x|c1\t101\t42\t10S140M\t=\t300\t199\t"
        + b"A" * 150
        + b"\t"
        + b"I" * 150
        + b"\tAS:i:280\tNM:i:2"
    )
    parsed = list(iter_sam([line]))
    assert len(parsed) == 1
    rec = parsed[0]
    assert rec.target_id == "bacteria.x" and rec.contig == "c1"
    assert rec.pos == 100 and rec.aligned_bases == 140 and rec.edit_distance == 2
    assert rec.score == 280 and abs(rec.identity - (1 - 2 / 140)) < 1e-9
    checks += 1

    # Headers and unmapped records are skipped.
    assert list(iter_sam([b"@HD\tVN:1.0", b"r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*"])) == []
    checks += 1

    assert _merged_length([]) == 0
    assert _merged_length([(0, 10), (5, 20), (30, 40)]) == 30
    checks += 2

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} alignment evidence checks passed")
