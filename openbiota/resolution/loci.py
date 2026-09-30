"""Targeted locus mapping: subtype discrimination and locus architecture.

A gene-family fragment count from a translated protein search is the weakest
evidence there is. It cannot say *which* allele of a toxin is present, and it
cannot say whether a multi-gene biosynthetic island is intact or represented by
one stray gene. Those are the two questions this stage answers, by mapping the
sample's reads onto exact accession-versioned nucleotide references.

Two targets, two different questions:

**bft subtype.** The three alleles (bft1/bft2/bft3) are ~1.5 kb and highly
similar, so reads map to all of them. A subtype is only called when one allele
is better supported than the runner-up by a margin — otherwise the honest
answer is "bft present, subtype not discriminated", which is what a fragment
count already said.

**pks island architecture.** The island is 55,140 bp across clbA-S. Colibactin
capacity needs the locus, not a gene, so coverage is measured in windows: a
contiguous run across most of the island means something quite different from
the same number of reads scattered over 5% of it. The window profile is
reported, never collapsed to a yes/no.

Both refuse to over-read in the same way the rest of the resolution layer does:
breadth and depth are reported with every call, a non-detection carries no
calibrated limit, and nothing here attributes a locus to a carrier organism —
that still needs linkage evidence.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from .schema import (
    CANDIDATE_SEQUENCE,
    COMPLETED,
    NOT_DETECTED_UNVALIDATED,
    REFERENCE_UNAVAILABLE,
    SUPPORTED_DETECTION,
    UNRESOLVED,
    Coverage,
    Linkage,
    Provenance,
    ResolutionCall,
    Validation,
)

PANELS: Final = Path("refs/panels/mechanisms")

#: Window size for the island architecture profile.
WINDOW_BP: Final = 1_000

#: A window counts as covered at or above this depth.
MIN_WINDOW_DEPTH: Final = 1

#: Breadth below this and the locus is a scattered signal, not a locus.
ARCHITECTURE_BREADTH: Final = 0.60

#: One allele must beat the runner-up by this much breadth to be named.
SUBTYPE_MARGIN: Final = 0.15

#: Minimum breadth before any subtype call is attempted at all.
SUBTYPE_MIN_BREADTH: Final = 0.50


@dataclass(frozen=True, slots=True)
class LocusTarget:
    """One nucleotide locus to map against."""

    target_id: str
    label: str
    #: Reference stems under `refs/panels/mechanisms`, without `.fasta`.
    references: tuple[str, ...]
    question: str
    #: `subtype` picks between alleles; `architecture` profiles one locus.
    mode: str
    accessions: tuple[str, ...]
    #: What a confirmed call would still require beyond this stage.
    still_required: str


TARGETS: Final[tuple[LocusTarget, ...]] = (
    LocusTarget(
        target_id="locus.bft.subtype",
        label="B. fragilis toxin subtype (bft1/bft2/bft3)",
        references=("mech.bft1", "mech.bft2", "mech.bft3"),
        question="which bft allele is present, if one can be distinguished",
        mode="subtype",
        accessions=("AB026625.1", "AB026626.1", "AB026624.1"),
        still_required=(
            "BfPAI context and carrier linkage. An allele call names the toxin variant; it "
            "does not establish which organism carries it or that the gene is intact."
        ),
    ),
    LocusTarget(
        target_id="locus.pks_island.architecture",
        label="Colibactin pks island architecture (clbA-S, 55,140 bp)",
        references=("mech.pks_island",),
        question="how much of the island is present, and contiguously",
        mode="architecture",
        accessions=("AM229678.1",),
        still_required=(
            "Carrier linkage and per-gene integrity. Breadth across the island is necessary "
            "for colibactin capacity but not sufficient: a disrupted clb gene inside an "
            "otherwise complete locus would still abolish the product."
        ),
    ),
)

BY_TARGET: Final[Mapping[str, LocusTarget]] = {t.target_id: t for t in TARGETS}


def _tools_present() -> bool:
    return all(shutil.which(t) for t in ("bowtie2", "bowtie2-build", "samtools"))


@dataclass(frozen=True, slots=True)
class RefCoverage:
    """Coverage of one reference by the sample's reads."""

    reference: str
    length: int
    covered_bases: int
    total_depth: int
    windows: tuple[int, ...]
    reads: int

    @property
    def breadth(self) -> float:
        return self.covered_bases / self.length if self.length else 0.0

    @property
    def mean_depth(self) -> float:
        return self.total_depth / self.length if self.length else 0.0

    @property
    def longest_run_windows(self) -> int:
        best = current = 0
        for depth in self.windows:
            current = current + 1 if depth >= MIN_WINDOW_DEPTH else 0
            best = max(best, current)
        return best


def _depth_profile(bam: Path, lengths: Mapping[str, int]) -> dict[str, RefCoverage]:
    """Per-reference coverage from `samtools depth`."""
    covered: dict[str, int] = dict.fromkeys(lengths, 0)
    total: dict[str, int] = dict.fromkeys(lengths, 0)
    windows: dict[str, list[int]] = {
        name: [0] * (length // WINDOW_BP + 1) for name, length in lengths.items()
    }
    result = subprocess.run(
        ["samtools", "depth", "-a", str(bam)], capture_output=True, text=True, check=True
    )
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name, position, depth_str = parts[0], int(parts[1]), int(parts[2])
        if name not in lengths:
            continue
        if depth_str > 0:
            covered[name] += 1
            total[name] += depth_str
            index = min((position - 1) // WINDOW_BP, len(windows[name]) - 1)
            windows[name][index] += depth_str

    counts = subprocess.run(
        ["samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True
    )
    reads: dict[str, int] = {}
    for line in counts.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] in lengths:
            reads[parts[0]] = int(parts[2])

    return {
        name: RefCoverage(
            reference=name,
            length=length,
            covered_bases=covered[name],
            total_depth=total[name],
            windows=tuple(windows[name]),
            reads=reads.get(name, 0),
        )
        for name, length in lengths.items()
    }


def map_locus(
    target: LocusTarget,
    *,
    reads: Sequence[Path],
    threads: int = 8,
    panels: Path = PANELS,
) -> dict[str, RefCoverage] | None:
    """Align the sample's reads to a target's references.

    Returns per-reference coverage, or None when the stage cannot run. A
    competitive index over all of a target's references at once is deliberate:
    for the bft alleles it is what makes the best-supported one meaningful.
    """
    fastas = [panels / f"{name}.fasta" for name in target.references]
    if not _tools_present() or not all(f.exists() for f in fastas):
        return None

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        combined = work / "refs.fasta"
        lengths: dict[str, int] = {}
        with combined.open("w") as handle:
            for name, fasta in zip(target.references, fastas, strict=True):
                seq = "".join(
                    line.strip()
                    for line in fasta.read_text().splitlines()
                    if not line.startswith(">")
                )
                lengths[name] = len(seq)
                handle.write(f">{name}\n")
                for i in range(0, len(seq), 70):
                    handle.write(seq[i : i + 70] + "\n")

        index = work / "idx"
        subprocess.run(
            ["bowtie2-build", "--quiet", "--threads", str(threads),
             str(combined), str(index)],
            check=True, capture_output=True,
        )
        sam = work / "aln.sam"
        argv = [
            "bowtie2", "--quiet", "--threads", str(threads),
            "--very-sensitive", "--no-unal", "-x", str(index), "-S", str(sam),
        ]
        if len(reads) >= 2:
            argv += ["-1", str(reads[0]), "-2", str(reads[1])]
        else:
            argv += ["-U", str(reads[0])]
        subprocess.run(argv, check=True, capture_output=True)

        bam = work / "aln.bam"
        subprocess.run(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)],
            check=True, capture_output=True,
        )
        subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
        return _depth_profile(bam, lengths)


def call_for(
    target: LocusTarget,
    coverage: Mapping[str, RefCoverage] | None,
    *,
    sample_id: str,
) -> ResolutionCall:
    """Turn locus coverage into a typed call, refusing to over-read it."""
    base = {
        "sample_id": sample_id,
        "target_id": target.target_id,
        "identity_kind": "allele" if target.mode == "subtype" else "operon",
    }
    if coverage is None:
        return ResolutionCall(
            **base,
            assay_status=REFERENCE_UNAVAILABLE,
            analytical_call=UNRESOLVED,
            reason_codes=("mapping_stage_unavailable",),
            plain=(
                f"{target.label}: the locus-mapping stage could not run (bowtie2, samtools "
                "or the reference is missing), so the question — "
                f"{target.question} — was not answered. Not a negative result."
            ),
        )

    provenance = Provenance(
        reference_accession_versions=target.accessions,
        tool_version="bowtie2 + samtools",
        parameters={
            "preset": "--very-sensitive",
            "window_bp": WINDOW_BP,
            "competitive_index": list(target.references),
        },
        threshold_version="locus/1",
    )
    validation = Validation(
        analytical_status="not_locally_calibrated",
        reference_function_status="experimentally_characterized_reference",
        phenotype_measured_in_sample=False,
        clinical_predictive_status="not_established",
    )
    total_reads = sum(c.reads for c in coverage.values())

    if total_reads == 0:
        return ResolutionCall(
            **base,
            assay_status=COMPLETED,
            call_id=f"{sample_id}:{target.target_id}",
            analytical_call=NOT_DETECTED_UNVALIDATED,
            reason_codes=("no_reads_mapped", "detection_limit_not_calibrated"),
            coverage=Coverage(independent_fragments=0),
            provenance=provenance,
            validation=validation,
            plain=(
                f"{target.label}: no reads mapped to "
                f"{', '.join(target.accessions)}. This stage has no calibrated detection "
                "limit, so that is a non-detection rather than a demonstrated absence."
            ),
        )

    if target.mode == "subtype":
        ranked = sorted(coverage.values(), key=lambda c: -c.breadth)
        best, runner_up = ranked[0], (ranked[1] if len(ranked) > 1 else None)
        margin = best.breadth - (runner_up.breadth if runner_up else 0.0)
        discriminated = (
            best.breadth >= SUBTYPE_MIN_BREADTH and margin >= SUBTYPE_MARGIN
        )
        allele = best.reference.replace("mech.", "")
        detail = "; ".join(
            f"{c.reference.replace('mech.', '')} {100 * c.breadth:.0f}% breadth, "
            f"{c.mean_depth:.1f}x"
            for c in ranked
        )
        if discriminated:
            plain = (
                f"{target.label}: best supported by <b>{allele}</b> "
                f"({100 * best.breadth:.0f}% of the allele covered at "
                f"{best.mean_depth:.1f}x), ahead of the runner-up by "
                f"{100 * margin:.0f} points of breadth. All alleles: {detail}. "
                + target.still_required
            )
            reasons: tuple[str, ...] = ("subtype_discriminated", "carrier_unresolved")
            call = CANDIDATE_SEQUENCE
        else:
            plain = (
                f"{target.label}: reads map to the alleles but no subtype can be "
                f"distinguished — the best is {allele} at {100 * best.breadth:.0f}% breadth, "
                f"only {100 * margin:.0f} points ahead of the next. These alleles are ~1.5 kb "
                f"and highly similar, so that margin is not enough to name one. All alleles: "
                f"{detail}. " + target.still_required
            )
            reasons = (
                "subtype_not_discriminated", "insufficient_discriminating_margin",
                "carrier_unresolved",
            )
            call = CANDIDATE_SEQUENCE
        return ResolutionCall(
            **base,
            assay_status=COMPLETED,
            call_id=f"{sample_id}:{target.target_id}",
            analytical_call=call,
            reason_codes=reasons,
            coverage=Coverage(
                independent_fragments=total_reads,
                target_breadth=round(best.breadth, 4),
                median_depth=round(best.mean_depth, 3),
                callable_bases=best.covered_bases,
                discriminatory_bases=best.covered_bases if discriminated else 0,
            ),
            linkage=Linkage(state="sample_co_detection"),
            provenance=provenance,
            validation=validation,
            plain=plain,
        )

    # architecture
    only = next(iter(coverage.values()))
    windows_covered = sum(1 for d in only.windows if d >= MIN_WINDOW_DEPTH)
    n_windows = len(only.windows)
    intact = only.breadth >= ARCHITECTURE_BREADTH
    plain = (
        f"{target.label}: {100 * only.breadth:.1f}% of the island covered "
        f"({only.covered_bases:,} of {only.length:,} bp) at {only.mean_depth:.2f}x mean "
        f"depth, from {only.reads:,} mapped reads. Contiguity: {windows_covered} of "
        f"{n_windows} one-kilobase windows have coverage, longest unbroken run "
        f"{only.longest_run_windows} windows. "
        + (
            "That is broad enough to discuss the locus rather than a gene, though "
            "per-gene integrity is still unassessed. "
            if intact else
            "That is a scattered signal, not a locus: fragments spread thinly over an "
            "island of this size do not establish that the island is present. "
        )
        + target.still_required
    )
    return ResolutionCall(
        **base,
        assay_status=COMPLETED,
        call_id=f"{sample_id}:{target.target_id}",
        analytical_call=SUPPORTED_DETECTION if intact else CANDIDATE_SEQUENCE,
        reason_codes=(
            ("locus_breadth_supported", "per_gene_integrity_unassessed", "carrier_unresolved")
            if intact else
            ("scattered_fragments_only", "locus_architecture_not_supported",
             "carrier_unresolved")
        ),
        coverage=Coverage(
            independent_fragments=only.reads,
            target_breadth=round(only.breadth, 4),
            median_depth=round(only.mean_depth, 3),
            callable_bases=only.covered_bases,
        ),
        linkage=Linkage(state="sample_co_detection"),
        provenance=provenance,
        validation=validation,
        plain=plain,
    )


def run_all(
    *, sample_id: str, reads: Sequence[Path], threads: int = 8, panels: Path = PANELS
) -> list[ResolutionCall]:
    """Map every locus target for one sample."""
    out: list[ResolutionCall] = []
    for target in TARGETS:
        try:
            coverage = map_locus(target, reads=reads, threads=threads, panels=panels)
        except (subprocess.CalledProcessError, OSError):
            coverage = None
        out.append(call_for(target, coverage, sample_id=sample_id))
    return out


def self_test() -> int:
    """Check the refusals without invoking an aligner. Returns failures."""
    failures = 0

    bft = BY_TARGET["locus.bft.subtype"]
    pks = BY_TARGET["locus.pks_island.architecture"]

    # An unavailable stage is never a negative.
    call = call_for(bft, None, sample_id="S")
    if call.assay_status != REFERENCE_UNAVAILABLE or call.is_negative:
        failures += 1
    if "Not a negative result" not in call.plain:
        failures += 1

    def cov(name: str, length: int, covered: int, depth: int, reads: int) -> RefCoverage:
        n = length // WINDOW_BP + 1
        filled = max(0, min(n, covered // WINDOW_BP))
        return RefCoverage(name, length, covered, depth * length, (1,) * filled + (0,) * (n - filled), reads)

    # Zero reads: non-detection with no calibrated limit.
    empty = {r: cov(r, 1546, 0, 0, 0) for r in bft.references}
    call = call_for(bft, empty, sample_id="S")
    if call.analytical_call != NOT_DETECTED_UNVALIDATED:
        failures += 1
    if call.coverage.limit_of_detection is not None:
        failures += 1

    # Similar alleles: no subtype may be named.
    tie = {
        "mech.bft1": cov("mech.bft1", 1546, 1400, 12, 300),
        "mech.bft2": cov("mech.bft2", 1546, 1380, 11, 290),
        "mech.bft3": cov("mech.bft3", 1547, 1370, 11, 285),
    }
    call = call_for(bft, tie, sample_id="S")
    if "subtype_not_discriminated" not in call.reason_codes:
        failures += 1
    if call.coverage.discriminatory_bases != 0:
        failures += 1
    if "not enough to name one" not in call.plain:
        failures += 1

    # A clear winner may be named.
    clear = {
        "mech.bft1": cov("mech.bft1", 1546, 1500, 20, 400),
        "mech.bft2": cov("mech.bft2", 1546, 400, 3, 60),
        "mech.bft3": cov("mech.bft3", 1547, 380, 3, 55),
    }
    call = call_for(bft, clear, sample_id="S")
    if "subtype_discriminated" not in call.reason_codes:
        failures += 1
    if "bft1" not in call.plain:
        failures += 1
    # Even a named subtype leaves the carrier open.
    if call.linkage.state != "sample_co_detection":
        failures += 1
    if "carrier_unresolved" not in call.reason_codes:
        failures += 1

    # Acceptance 31, at nucleotide level: scattered reads are not an island.
    scattered = {"mech.pks_island": cov("mech.pks_island", 55140, 1200, 1, 20)}
    call = call_for(pks, scattered, sample_id="S")
    if call.analytical_call != CANDIDATE_SEQUENCE:
        failures += 1
    if "scattered_fragments_only" not in call.reason_codes:
        failures += 1
    if "not a locus" not in call.plain:
        failures += 1

    # Broad coverage supports the locus but not per-gene integrity.
    broad = {"mech.pks_island": cov("mech.pks_island", 55140, 50000, 8, 3000)}
    call = call_for(pks, broad, sample_id="S")
    if call.analytical_call != SUPPORTED_DETECTION:
        failures += 1
    if "per_gene_integrity_unassessed" not in call.reason_codes:
        failures += 1
    if "integrity is still unassessed" not in call.plain:
        failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
