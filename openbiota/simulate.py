"""Paired-end read simulator for synthetic validation communities.

Deliberately simple and deterministic. The goal is not to model an Illumina
run in detail — it is to produce reads whose *provenance is known exactly*, so
that what the pipeline reports can be compared against ground truth.

Matched to the screened data: 100-151 bp reads with the same length
distribution, Phred+33 with NovaSeq-style binned qualities, and a uniform
substitution error rate.
"""

from __future__ import annotations

import gzip
import math
import random
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from openbiota.seqio import iter_fasta

COMPLEMENT: Final = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES: Final = "ACGT"

#: Read-length distribution measured on the screened sample: mean 148.3 bp
#: over a 100-151 bp range, i.e. mostly full-length with a trimmed tail.
LENGTH_CHOICES: Final = (151,) * 88 + (150,) * 3 + (140,) * 3 + (130,) * 2 + (120,) * 2 + (110,) * 1 + (100,) * 1

#: NovaSeq emits a handful of binned quality values rather than a continuum.
QUALITY_BINS: Final = (("I", 0.93), ("9", 0.05), ("-", 0.015), ("#", 0.005))

#: Typical insert size for a shotgun library.
MEAN_INSERT: Final = 350
INSERT_SD: Final = 80


@dataclass(frozen=True, slots=True)
class GenomeSpec:
    """One community member: a genome and its relative cell abundance."""

    name: str
    path: Path
    abundance: float
    genome_size: int = 0


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def load_contigs(path: Path) -> list[str]:
    """Read a (possibly gzipped) genome FASTA into a list of contig sequences."""
    if path.suffix.lower() == ".gz":
        text = gzip.decompress(path.read_bytes()).decode("utf-8", "replace")
        contigs: list[str] = []
        current: list[str] = []
        for line in text.splitlines():
            if line.startswith(">"):
                if current:
                    contigs.append("".join(current))
                current = []
            elif line:
                current.append(line.strip())
        if current:
            contigs.append("".join(current))
        return [c.upper() for c in contigs if c]
    return [seq.upper() for _, seq in iter_fasta(path) if seq]


#: Substitution lookup, so mutating a base is a dict hit rather than a list
#: comprehension plus a choice.
_SUBSTITUTIONS: Final = {b: tuple(x for x in BASES if x != b) for b in BASES}

#: Pre-generated pool of quality strings. Drawing from a pool rather than
#: generating per read matters: at 3 million pairs, per-character RNG calls
#: dominate the whole run.
_QUALITY_POOL_SIZE: Final = 512


def _build_quality_pool(rng: random.Random) -> dict[int, list[str]]:
    symbols = [c for c, _ in QUALITY_BINS]
    weights = [w for _, w in QUALITY_BINS]
    pool: dict[int, list[str]] = {}
    for length in set(LENGTH_CHOICES):
        pool[length] = [
            "".join(rng.choices(symbols, weights=weights, k=length))
            for _ in range(_QUALITY_POOL_SIZE)
        ]
    return pool


def _mutate(seq: str, n_errors: int, rng: random.Random) -> str:
    """Apply exactly ``n_errors`` substitutions at random positions."""
    if n_errors <= 0:
        return seq
    out = list(seq)
    length = len(out)
    for _ in range(n_errors):
        position = rng.randrange(length)
        alternatives = _SUBSTITUTIONS.get(out[position])
        if alternatives is not None:
            out[position] = alternatives[rng.randrange(3)]
    return "".join(out)


def simulate_pairs(
    genomes: Sequence[GenomeSpec],
    *,
    n_pairs: int,
    seed: int,
    error_rate: float,
) -> Iterator[tuple[str, str, str, str, str]]:
    """Yield ``(read_id, seq1, qual1, seq2, qual2)`` for a synthetic community.

    Reads are drawn in proportion to each genome's *DNA mass*, which is its
    cell abundance times its genome size. That is what a real shotgun library
    samples, and getting it wrong would bias the calibration this exists to
    measure.
    """
    rng = random.Random(seed)

    loaded: list[tuple[GenomeSpec, list[str], int]] = []
    for spec in genomes:
        contigs = load_contigs(spec.path)
        if not contigs:
            continue
        size = sum(len(c) for c in contigs)
        loaded.append((spec, contigs, size))
    if not loaded:
        raise ValueError("no usable genomes supplied to the simulator")

    mass = [spec.abundance * size for spec, _, size in loaded]
    total_mass = sum(mass)
    if total_mass <= 0:
        raise ValueError("community has zero total DNA mass")

    quality_pool = _build_quality_pool(rng)
    pool_last = _QUALITY_POOL_SIZE - 1

    # Draw genome and contig assignments in bulk; per-item rng.choices calls
    # are the single biggest cost at millions of reads.
    genome_picks = rng.choices(range(len(loaded)), weights=mass, k=n_pairs)
    contig_pickers = [
        (contigs, [len(c) for c in contigs]) for _, contigs, _ in loaded
    ]
    lengths1 = rng.choices(LENGTH_CHOICES, k=n_pairs)
    lengths2 = rng.choices(LENGTH_CHOICES, k=n_pairs)
    # Errors per read: mean = read_length * error_rate, which at these rates is
    # well approximated by a Poisson draw and costs one call instead of 151.
    mean_errors = error_rate * 151

    for index in range(n_pairs):
        which = genome_picks[index]
        spec, _, _ = loaded[which]
        contigs, contig_weights = contig_pickers[which]
        contig = (
            contigs[0]
            if len(contigs) == 1
            else rng.choices(contigs, weights=contig_weights, k=1)[0]
        )

        len1 = lengths1[index]
        len2 = lengths2[index]
        insert = max(int(rng.gauss(MEAN_INSERT, INSERT_SD)), len1, len2, 200)
        if len(contig) <= insert:
            continue
        start = rng.randrange(0, len(contig) - insert)
        fragment = contig[start : start + insert]
        if rng.random() < 0.5:
            fragment = reverse_complement(fragment)

        seq1 = fragment[:len1]
        seq2 = reverse_complement(fragment)[:len2]
        if mean_errors > 0:
            seq1 = _mutate(seq1, _poisson(mean_errors, rng), rng)
            seq2 = _mutate(seq2, _poisson(mean_errors, rng), rng)

        # Mate ids are identical with no /1 /2 suffix, matching real NovaSeq
        # output, so mate collapse is exercised the same way.
        yield (
            f"SIM:{spec.name}:{index}",
            seq1,
            quality_pool[len1][rng.randint(0, pool_last)][: len(seq1)],
            seq2,
            quality_pool[len2][rng.randint(0, pool_last)][: len(seq2)],
        )


def _poisson(mean: float, rng: random.Random) -> int:
    """Knuth's Poisson sampler. Cheap at the small means used here."""
    if mean <= 0:
        return 0
    limit = math.exp(-mean)
    product = rng.random()
    count = 0
    while product > limit:
        count += 1
        product *= rng.random()
    return count


def write_community(
    genomes: Sequence[GenomeSpec],
    *,
    out_prefix: Path,
    n_pairs: int,
    seed: int,
    error_rate: float = 0.002,
) -> tuple[Path, Path, int]:
    """Write a simulated paired FASTQ. Returns ``(r1, r2, pairs_written)``."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    r1 = out_prefix.with_name(out_prefix.name + "_1.fastq")
    r2 = out_prefix.with_name(out_prefix.name + "_2.fastq")
    marker = out_prefix.with_name(out_prefix.name + ".done")

    signature = f"{n_pairs}|{seed}|{error_rate}|" + ";".join(
        f"{g.name}:{g.abundance}" for g in genomes
    )
    if marker.is_file() and marker.read_text(encoding="utf-8").strip() == signature:
        return r1, r2, n_pairs

    written = 0
    with r1.open("w", encoding="utf-8") as f1, r2.open("w", encoding="utf-8") as f2:
        for read_id, seq1, qual1, seq2, qual2 in simulate_pairs(
            genomes, n_pairs=n_pairs, seed=seed, error_rate=error_rate
        ):
            f1.write(f"@{read_id}\n{seq1}\n+\n{qual1}\n")
            f2.write(f"@{read_id}\n{seq2}\n+\n{qual2}\n")
            written += 1

    marker.write_text(signature, encoding="utf-8")
    return r1, r2, written
