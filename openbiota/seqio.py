"""Tiny FASTA helpers shared by the reference and residue modules."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path


def clean_sequence(sequence: str) -> str:
    """Keep only alphabetic residues; DIAMOND rejects anything else."""
    return "".join(c for c in sequence if c.isalpha()).upper()


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8"
    )


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs, streaming rather than slurping."""
    header: str | None = None
    buffer: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buffer)
                header = line[1:]
                buffer = []
            elif line:
                buffer.append(line)
    if header is not None:
        yield header, "".join(buffer)
