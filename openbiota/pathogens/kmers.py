"""Frozen k-mer engine for candidate generation (BUILD_SPEC_v05.0 §5.2 step 4).

Candidate generation needs to touch every read of every eligible library
without abundance-filtering or downsampling, so it has to be fast; but it is
*only* candidate generation. Nothing here is a positive call. A read that
matches is a read worth aligning competitively, and that is all.

**Scheme.** Canonical 31-mers, selected by FracMinHash: a k-mer is kept when
`splitmix64(kmer) < 2**64 / SCALED`. The same hash and the same threshold are
applied to references and to reads, so selection is consistent by
construction and the whole scheme is one frozen integer (`SCALED`) rather than
a tuning knob. With `SCALED = 8` a 150 bp read retains ~15 k-mers on average
and the chance of retaining none is about one in ten thousand.

**Why not minimizers.** A true minimizer scheme needs a sliding-window
minimum per position; FracMinHash needs one comparison. On 2.4 Gbp of reads
that difference is the difference between minutes and an hour, and the
sensitivity cost at this density is negligible for a stage whose output is
"align this read properly".

**Ownership, not classification.** Each retained reference k-mer records which
target it came from. A k-mer seen in more than one target — or in the host or
a food or near-neighbour decoy — is marked `SHARED` rather than assigned to
whichever genome happened to be read first. Reads supported only by shared
k-mers become ambiguous evidence at the rank the data supports; they are never
promoted to a species.

Bit-field extraction: bases are packed 2 bits each into a uint64 stream, and
each k-mer is read as an unaligned 62-bit field spanning at most two words.
That costs a handful of vector operations per position instead of 31
shift-or steps, which is what makes a full-depth pass affordable.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np

__all__ = [
    "K",
    "MASK",
    "SCALED",
    "SHARED",
    "KmerIndex",
    "hash64",
    "iter_fasta",
    "select_kmers",
    "self_test",
]

#: k-mer length. 31 fits a uint64 with two bits to spare and is the de-facto
#: standard for species-level nucleotide identity.
K: Final[int] = 31

#: FracMinHash denominator. Frozen: it is part of the reference lock, and
#: changing it changes every count, so it gets a new bundle id.
SCALED: Final[int] = 8

#: Sentinel owner meaning "this k-mer occurs in more than one target, or in a
#: host or decoy reference". Never a taxonomic assignment.
SHARED: Final[int] = 0xFFFF

MASK: Final[int] = (1 << (2 * K)) - 1
_THRESHOLD: Final[int] = (1 << 64) // SCALED

_U64: Final = np.uint64
_SHIFT2: Final = _U64(2)


# A→0, C→1, G→2, T→3; everything else 255 and therefore rejected.
_CODES: Final[np.ndarray] = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODES[_b] = _i
    _CODES[_b + 32] = _i  # lower case


def hash64(values: np.ndarray) -> np.ndarray:
    """SplitMix64 finalizer, vectorised.

    A fixed, well-mixed integer hash. It must never change: it defines which
    k-mers exist in the index, so a different hash is a different reference
    bundle, not a refactor.
    """
    x = values.astype(np.uint64, copy=True)
    x ^= x >> _U64(30)
    x *= _U64(0xBF58476D1CE4E5B9)
    x ^= x >> _U64(27)
    x *= _U64(0x94D049BB133111EB)
    x ^= x >> _U64(31)
    return x


def _pack(codes: np.ndarray) -> np.ndarray:
    """Pack 2-bit base codes into a big-endian uint64 bit stream."""
    n = codes.size
    words = (n + 31) // 32 + 1  # one spare word so the last field can span
    padded = np.zeros(words * 32, dtype=np.uint64)
    padded[:n] = codes
    grouped = padded.reshape(words, 32)
    out = np.zeros(words, dtype=np.uint64)
    for i in range(32):
        out = (out << _SHIFT2) | grouped[:, i]
    return out


def _extract(packed: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Read the 2K-bit field beginning at each base position in `starts`."""
    bit = starts.astype(np.uint64) * _SHIFT2
    word = bit >> _U64(6)  # // 64
    off = bit & _U64(63)
    hi = packed[word.astype(np.int64)]
    lo = packed[(word + _U64(1)).astype(np.int64)]
    # `hi << off` drops the bits before the field; the complementary shift of
    # `lo` supplies the tail when the field straddles a word boundary. numpy
    # shifts by 64 are undefined, so the off==0 case is masked out explicitly.
    shifted = hi << off
    tail = np.where(off == _U64(0), _U64(0), lo >> (_U64(64) - off))
    return (shifted | tail) >> _U64(64 - 2 * K)


def _revcomp(kmers: np.ndarray) -> np.ndarray:
    """Reverse-complement packed k-mers by reversing 2-bit groups.

    Complementing is `3 - base`, i.e. a bitwise NOT of the 2K used bits.
    Reversal is done with the standard pairwise/quad/byte swap ladder, then a
    shift to discard the unused high bits.
    """
    x = (~kmers) & _U64(MASK)
    # swap adjacent 2-bit groups into reversed order
    x = ((x & _U64(0x3333333333333333)) << _U64(2)) | (
        (x >> _U64(2)) & _U64(0x3333333333333333)
    )
    x = ((x & _U64(0x0F0F0F0F0F0F0F0F)) << _U64(4)) | (
        (x >> _U64(4)) & _U64(0x0F0F0F0F0F0F0F0F)
    )
    x = ((x & _U64(0x00FF00FF00FF00FF)) << _U64(8)) | (
        (x >> _U64(8)) & _U64(0x00FF00FF00FF00FF)
    )
    x = ((x & _U64(0x0000FFFF0000FFFF)) << _U64(16)) | (
        (x >> _U64(16)) & _U64(0x0000FFFF0000FFFF)
    )
    x = (x << _U64(32)) | (x >> _U64(32))
    return x >> _U64(64 - 2 * K)


def select_kmers(
    sequence: bytes | bytearray, *, positions: bool = False
) -> tuple[np.ndarray, np.ndarray | None]:
    """Canonical FracMinHash k-mers of one sequence.

    Returns ``(kmers, positions)``. Positions are the 0-based start offsets of
    the retained k-mers and are only computed when asked for, because the
    read-side pass does not need them.

    Any window containing a non-ACGT base is dropped rather than substituted.
    Imputing an ambiguity code as a real base is how low-complexity and
    masked regions leak into "support".
    """
    codes = _CODES[np.frombuffer(bytes(sequence), dtype=np.uint8)]
    n = codes.size
    if n < K:
        empty = np.empty(0, dtype=np.uint64)
        return empty, (np.empty(0, dtype=np.int64) if positions else None)

    valid = codes != 255
    clean = np.where(valid, codes, 0).astype(np.uint64)
    packed = _pack(clean)

    starts = np.arange(n - K + 1, dtype=np.int64)
    fwd = _extract(packed, starts)

    # A k-mer is usable only when all K of its bases were ACGT. A sliding
    # all() is a difference of prefix sums over the invalid mask.
    bad = np.cumsum(~valid, dtype=np.int64)
    bad_in_window = bad[K - 1 :] - np.concatenate(([0], bad[: n - K]))
    usable = bad_in_window == 0
    if not usable.any():
        empty = np.empty(0, dtype=np.uint64)
        return empty, (np.empty(0, dtype=np.int64) if positions else None)

    fwd = fwd[usable]
    starts = starts[usable]
    canonical = np.minimum(fwd, _revcomp(fwd))
    keep = hash64(canonical) < _U64(_THRESHOLD)
    return canonical[keep], (starts[keep] if positions else None)


# --------------------------------------------------------------------------- #
# FASTA reading
# --------------------------------------------------------------------------- #


def iter_fasta(path: Path | str) -> Iterator[tuple[str, bytes]]:
    """Yield ``(header, sequence)`` from a plain or gzipped FASTA.

    Contigs are yielded separately so that no k-mer ever spans a contig
    boundary and invents sequence that does not exist.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    chunks: list[bytes] = []
    with opener(path, "rb") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(chunks)
                header = line[1:].strip().decode("utf-8", "replace")
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, b"".join(chunks)


# --------------------------------------------------------------------------- #
# Index
# --------------------------------------------------------------------------- #


@dataclass
class KmerIndex:
    """A sorted, memory-mappable k-mer → owner index.

    `kmers` is sorted ascending and unique; `owners[i]` is either a target
    index or `SHARED`. Lookup is a single `searchsorted`, which is why a
    full-depth read pass costs about a minute rather than an hour.
    """

    kmers: np.ndarray
    owners: np.ndarray
    target_ids: tuple[str, ...]
    informative_bins: dict[int, np.ndarray]
    scaled: int = SCALED
    k: int = K

    def __post_init__(self) -> None:
        if self.kmers.dtype != np.uint64 or self.owners.dtype != np.uint16:
            raise ValueError("kmer index dtypes must be uint64 / uint16")
        if self.kmers.size != self.owners.size:
            raise ValueError("kmer and owner arrays disagree in length")

    def lookup(self, query: np.ndarray) -> np.ndarray:
        """Map query k-mers to owner ids; absent k-mers map to -1."""
        if query.size == 0:
            return np.empty(0, dtype=np.int32)
        idx = np.searchsorted(self.kmers, query)
        idx_clipped = np.minimum(idx, self.kmers.size - 1)
        hit = self.kmers[idx_clipped] == query
        out = np.full(query.size, -1, dtype=np.int32)
        out[hit] = self.owners[idx_clipped[hit]].astype(np.int32)
        return out

    def is_informative_bin(self, target_index: int, bin_index: int) -> bool:
        """Whether a canonical 10-kb bin carries target-specific sequence.

        Used to enforce the spec's requirement that qualifying genomic support
        includes resolution-specific regions: five reads in a bin that is
        indistinguishable from a relative cannot establish a species.
        """
        bins = self.informative_bins.get(target_index)
        if bins is None or bins.size == 0:
            return False
        pos = np.searchsorted(bins, bin_index)
        return bool(pos < bins.size and bins[pos] == bin_index)

    # -- persistence -------------------------------------------------------- #

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload: dict[str, np.ndarray] = {
            "kmers": self.kmers,
            "owners": self.owners,
            "target_ids": np.array(self.target_ids, dtype=object),
            "meta": np.array([self.k, self.scaled], dtype=np.int64),
        }
        for index, bins in self.informative_bins.items():
            payload[f"bins_{index}"] = bins
        np.savez(path, **payload)

    @classmethod
    def load(cls, path: Path, *, mmap: bool = True) -> KmerIndex:
        data = np.load(path, allow_pickle=True, mmap_mode="r" if mmap else None)
        bins = {
            int(name.split("_", 1)[1]): np.asarray(data[name])
            for name in data.files
            if name.startswith("bins_")
        }
        meta = np.asarray(data["meta"])
        return cls(
            kmers=np.asarray(data["kmers"]),
            owners=np.asarray(data["owners"]),
            target_ids=tuple(str(t) for t in np.asarray(data["target_ids"])),
            informative_bins=bins,
            k=int(meta[0]),
            scaled=int(meta[1]),
        )


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #


def _naive_canonical(seq: str) -> list[int]:
    table = str.maketrans("ACGT", "TGCA")
    out: list[int] = []
    for i in range(len(seq) - K + 1):
        window = seq[i : i + K]
        if any(c not in "ACGT" for c in window):
            continue
        rev = window.translate(table)[::-1]
        enc = lambda s: int("".join("ACGT".index(c).__format__("02b") for c in s), 2)  # noqa: E731
        out.append(min(enc(window), enc(rev)))
    return out


def self_test() -> int:
    rng = np.random.default_rng(20260907)
    checks = 0

    # Bit-field extraction and canonicalisation must match a naive encoder
    # exactly. This is the load-bearing correctness claim of the module.
    for _ in range(20):
        seq = "".join(rng.choice(list("ACGT"), size=rng.integers(K, 400)))
        expected = _naive_canonical(seq)
        kept = [c for c in expected if int(hash64(np.array([c], np.uint64))[0]) < _THRESHOLD]
        got, _ = select_kmers(seq.encode())
        assert sorted(int(x) for x in got) == sorted(kept), seq[:60]
        checks += 1

    # Reverse complement of a sequence yields exactly the same canonical set.
    seq = "".join(rng.choice(list("ACGT"), size=5_000))
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    a, _ = select_kmers(seq.encode())
    b, _ = select_kmers(rc.encode())
    assert sorted(a.tolist()) == sorted(b.tolist())
    checks += 1

    # Ambiguity codes are dropped, not imputed: a run of Ns must remove every
    # window that touches it rather than encoding as poly-A.
    with_n = seq[:100] + "N" * 5 + seq[105:]
    c, _ = select_kmers(with_n.encode())
    assert c.size < a.size
    poly_a, _ = select_kmers(("A" * 200).encode())
    assert not set(c.tolist()) & set(poly_a.tolist())
    checks += 1

    # Density is close to the frozen 1/SCALED.
    long_seq = "".join(rng.choice(list("ACGT"), size=2_000_000))
    kmers, positions = select_kmers(long_seq.encode(), positions=True)
    density = kmers.size / (len(long_seq) - K + 1)
    assert 0.8 / SCALED < density < 1.25 / SCALED, density
    assert positions is not None and positions.size == kmers.size
    checks += 2

    # Positions really point at the k-mers they claim to.
    sample = rng.choice(positions.size, size=25, replace=False)
    for i in sample:
        start = int(positions[i])
        window = long_seq[start : start + K]
        again, _ = select_kmers(window.encode())
        assert again.size == 1 and int(again[0]) == int(kmers[i])
    checks += 1

    # Index lookup: hits resolve, misses return -1, shared stays shared.
    uniq = np.unique(kmers[:1000])
    owners = np.zeros(uniq.size, dtype=np.uint16)
    owners[::3] = SHARED
    index = KmerIndex(uniq, owners, ("t0",), {0: np.array([0, 5], dtype=np.int64)})
    got = index.lookup(uniq)
    assert (got[::3] == SHARED).all()
    absent = np.array([np.uint64(1)], dtype=np.uint64)
    while absent[0] in set(uniq.tolist()):
        absent[0] += np.uint64(1)
    assert index.lookup(absent)[0] == -1
    assert index.lookup(np.empty(0, dtype=np.uint64)).size == 0
    checks += 3

    assert index.is_informative_bin(0, 5)
    assert not index.is_informative_bin(0, 6)
    assert not index.is_informative_bin(7, 5)
    checks += 3

    # Short sequences are handled without raising.
    short, pos = select_kmers(b"ACGT", positions=True)
    assert short.size == 0 and pos is not None and pos.size == 0
    checks += 1

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} k-mer engine checks passed")
