"""Active-site residue checks — partial recovery of a full-length filter.

The published urdA method separates urocanate reductase from fumarate
reductase flavoprotein by requiring tyrosine or methionine at FAD-binding-site
residue 373 of the *full-length* predicted protein. A 151 bp read translates to
about 50 residues, so that filter is mostly unavailable at read level.

What is recoverable:

1.  At database build time, align every target reference to a designated anchor
    protein and record, for each reference, which of *its own* positions
    corresponds to the anchor's canonical position. The alignment is done with
    ``diamond blastp`` — no extra dependency, and it returns the gapped
    alignment strings directly.

2.  At search time, DIAMOND reports subject alignment coordinates. For the
    subset of reads whose alignment happens to span the mapped position, read
    off the aligned query residue and test it.

The result is reported as its own subset ("N accepted, of which M spanned the
mapped position; K passed"), never folded into the headline count.
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import DependencyError, ReferenceFetchError
from openbiota.logging_util import Reporter
from openbiota.net import get
from openbiota.panels import ResidueCheck
from openbiota.seqio import clean_sequence, write_fasta

UNIPROT_ENTRY: Final = "https://rest.uniprot.org/uniprotkb"

GAP: Final = "-"

#: Outcome codes for a single residue check.
#: Sentinel query index used to align the anchor against itself as a self-test.
#: Negative so it can never collide with a real database index.
_SELF_CHECK_IDX: Final = -1

SPAN_PASS: Final = "pass"
SPAN_FAIL: Final = "fail"
SPAN_GAP: Final = "gap"
NOT_SPANNED: Final = "not_spanned"
NO_MAPPING: Final = "no_mapping"


@dataclass(frozen=True, slots=True)
class AnchorMapping:
    """Mapping from database index to that reference's own anchor-equivalent position."""

    panel: str
    anchor_accession: str
    anchor_length: int
    canonical_position: int
    canonical_residue: str
    accepted_residues: tuple[str, ...]
    #: database index -> 1-based position in that reference sequence
    positions: dict[int, int] = field(default_factory=dict)
    #: database index -> whether that reference's own residue at the mapped
    #: position is one of ``accepted_residues``. This is the reference-side
    #: half of the filter and, unlike the read-side check, it applies to every
    #: fragment rather than only to those that happen to span the site.
    residue_ok: dict[int, bool] = field(default_factory=dict)
    n_candidates: int = 0
    n_aligned: int = 0
    n_position_spanned: int = 0
    n_reference_residue_ok: int = 0

    def to_json(self) -> dict[str, Any]:
        return {
            "panel": self.panel,
            "anchor_accession": self.anchor_accession,
            "anchor_length": self.anchor_length,
            "canonical_position": self.canonical_position,
            "canonical_residue": self.canonical_residue,
            "accepted_residues": list(self.accepted_residues),
            "positions": {str(k): v for k, v in self.positions.items()},
            "residue_ok": sorted(k for k, v in self.residue_ok.items() if v),
            "n_candidates": self.n_candidates,
            "n_aligned": self.n_aligned,
            "n_position_spanned": self.n_position_spanned,
            "n_reference_residue_ok": self.n_reference_residue_ok,
        }

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> AnchorMapping:
        positions = {int(k): int(v) for k, v in payload["positions"].items()}
        ok = {int(k) for k in payload.get("residue_ok", [])}
        return cls(
            panel=payload["panel"],
            anchor_accession=payload["anchor_accession"],
            anchor_length=int(payload["anchor_length"]),
            canonical_position=int(payload["canonical_position"]),
            canonical_residue=payload.get("canonical_residue", "?"),
            accepted_residues=tuple(payload["accepted_residues"]),
            positions=positions,
            residue_ok={idx: idx in ok for idx in positions},
            n_candidates=int(payload.get("n_candidates", 0)),
            n_aligned=int(payload.get("n_aligned", 0)),
            n_position_spanned=int(payload.get("n_position_spanned", 0)),
            n_reference_residue_ok=int(payload.get("n_reference_residue_ok", 0)),
        )


def fetch_anchor_sequence(accession: str, cache_dir: Path, *, refresh: bool = False) -> str:
    """Download (and cache) one UniProt sequence by accession."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"anchor_{accession}.fasta"
    if path.is_file() and not refresh:
        text = path.read_text(encoding="utf-8")
    else:
        text = get(f"{UNIPROT_ENTRY}/{accession}.fasta").text()
        if not text.startswith(">"):
            raise ReferenceFetchError(
                f"anchor accession {accession!r}: UniProt did not return a FASTA record. "
                "Check 'anchor_accession' in the panel's residue_check block."
            )
        path.write_text(text, encoding="utf-8")
    sequence = clean_sequence("".join(text.splitlines()[1:]))
    if not sequence:
        raise ReferenceFetchError(f"anchor accession {accession!r}: empty sequence")
    return sequence


def map_position_through_alignment(
    query_gapped: str,
    subject_gapped: str,
    subject_start: int,
    query_start: int,
    target_subject_position: int,
) -> int | None:
    """Return the query position aligned to ``target_subject_position``.

    ``*_gapped`` are the two alignment strings of equal length; ``*_start`` are
    1-based coordinates of their first non-gap character. Returns ``None`` when
    the alignment does not span the position or aligns a gap to it.
    """
    if len(query_gapped) != len(subject_gapped):
        return None
    q_pos = query_start
    s_pos = subject_start
    for q_char, s_char in zip(query_gapped, subject_gapped, strict=True):
        if s_char != GAP and s_pos == target_subject_position:
            return None if q_char == GAP else q_pos
        if s_char != GAP:
            s_pos += 1
        if q_char != GAP:
            q_pos += 1
    return None


def residue_at_subject_position(
    query_gapped: str,
    subject_gapped: str,
    subject_start: int,
    target_subject_position: int,
) -> str | None:
    """Return the query *residue* aligned to ``target_subject_position``.

    Used at search time, where the "query" is the translated read. ``None``
    means the alignment does not span the position; ``GAP`` means it does but
    the read has a deletion there.
    """
    if len(query_gapped) != len(subject_gapped):
        return None
    s_pos = subject_start
    for q_char, s_char in zip(query_gapped, subject_gapped, strict=True):
        if s_char != GAP:
            if s_pos == target_subject_position:
                return q_char
            s_pos += 1
    return None


def classify_read_residue(
    *,
    anchor_position: int | None,
    subject_start: int,
    subject_end: int,
    query_gapped: str,
    subject_gapped: str,
    accepted: frozenset[str] | set[str] | Sequence[str],
) -> tuple[str, str | None]:
    """Classify one read hit against the mapped active-site position."""
    if anchor_position is None:
        return NO_MAPPING, None
    lo, hi = (subject_start, subject_end) if subject_start <= subject_end else (subject_end, subject_start)
    if not lo <= anchor_position <= hi:
        return NOT_SPANNED, None
    if not query_gapped or not subject_gapped:
        return NOT_SPANNED, None
    residue = residue_at_subject_position(query_gapped, subject_gapped, lo, anchor_position)
    if residue is None:
        return NOT_SPANNED, None
    if residue == GAP:
        return SPAN_GAP, GAP
    return (SPAN_PASS if residue.upper() in set(accepted) else SPAN_FAIL), residue.upper()


def map_anchor_positions(
    *,
    panel_name: str,
    residue_check: ResidueCheck,
    refs: Sequence[tuple[int, str, str]],
    work_dir: Path,
    reporter: Reporter,
    cache_dir: Path,
    diamond: str = "diamond",
    threads: int = 1,
    refresh: bool = False,
) -> AnchorMapping:
    """Align each reference to the anchor and record its equivalent position.

    ``refs`` is ``(database_index, accession, sequence)``.
    """
    anchor = fetch_anchor_sequence(residue_check.anchor_accession, cache_dir, refresh=refresh)
    position = residue_check.canonical_position
    if position > len(anchor):
        raise ReferenceFetchError(
            f"panel {panel_name!r}: residue_check.canonical_position {position} exceeds the "
            f"length of anchor {residue_check.anchor_accession} ({len(anchor)} aa). "
            "The anchor accession is almost certainly wrong for this position."
        )
    canonical_residue = anchor[position - 1]
    if canonical_residue not in residue_check.accepted_residues:
        reporter.warn(
            f"panel {panel_name!r}: anchor {residue_check.anchor_accession} carries "
            f"{canonical_residue!r} at position {position}, which is not in the accepted set "
            f"{sorted(residue_check.accepted_residues)}. Verify the anchor and the position."
        )
    else:
        reporter.record(
            f"    anchor {residue_check.anchor_accession} ({len(anchor)} aa) position {position} "
            f"= {canonical_residue!r} — consistent with accepted residues"
        )

    work_dir.mkdir(parents=True, exist_ok=True)
    anchor_fasta = work_dir / f"{panel_name}_anchor.faa"
    refs_fasta = work_dir / f"{panel_name}_refs.faa"
    anchor_db = work_dir / f"{panel_name}_anchor"
    out_tsv = work_dir / f"{panel_name}_anchor_hits.tsv"

    write_fasta(anchor_fasta, [(f"anchor{residue_check.anchor_accession}", anchor)])
    # The anchor is included in its own query set under a sentinel index. Its
    # mapped position must come back equal to the canonical position; if it
    # does not, the alignment walk is wrong and every other mapping is
    # untrustworthy. Cheap, and it catches off-by-one errors that would
    # otherwise silently shift the whole check by one residue.
    write_fasta(
        refs_fasta,
        [(str(idx), sequence) for idx, _, sequence in refs] + [(str(_SELF_CHECK_IDX), anchor)],
    )

    exe = shutil.which(diamond)
    if exe is None:  # pragma: no cover — checked upstream
        raise DependencyError(f"{diamond!r} not found on PATH")

    _run(
        [exe, "makedb", "--in", str(anchor_fasta), "--db", str(anchor_db),
         "--threads", str(threads), "--quiet"],
        reporter,
    )
    _run(
        [
            exe, "blastp",
            "--db", str(anchor_db),
            "--query", str(refs_fasta),
            "--out", str(out_tsv),
            "--outfmt", "6", "qseqid", "qstart", "qend", "sstart", "send",
            "qseq_gapped", "sseq_gapped",
            "--max-target-seqs", "1",
            "--index-chunks", "1",
            # These are within-family alignments over a handful of thousand
            # short sequences, so maximum sensitivity costs nothing here and
            # maximises how many references get a usable mapping.
            "--very-sensitive",
            "--evalue", "1e-3",
            "--threads", str(threads),
            "--quiet",
        ],
        reporter,
    )

    positions: dict[int, int] = {}
    residue_ok: dict[int, bool] = {}
    n_aligned = 0
    n_ok = 0
    sequences = {idx: sequence for idx, _, sequence in refs}
    sequences[_SELF_CHECK_IDX] = anchor

    with out_tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            try:
                idx = int(parts[0])
                q_start, q_end = int(parts[1]), int(parts[2])
                s_start, s_end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            del q_end
            n_aligned += 1
            if not (min(s_start, s_end) <= position <= max(s_start, s_end)):
                continue
            mapped = map_position_through_alignment(
                query_gapped=parts[5],
                subject_gapped=parts[6],
                subject_start=min(s_start, s_end),
                query_start=q_start,
                target_subject_position=position,
            )
            if mapped is None:
                continue
            sequence = sequences.get(idx, "")
            if not 1 <= mapped <= len(sequence):
                continue
            positions[idx] = mapped
            ok = sequence[mapped - 1].upper() in residue_check.accepted_residues
            residue_ok[idx] = ok
            if ok:
                n_ok += 1

    # ---- self-test ---------------------------------------------------------- #
    self_mapped = positions.pop(_SELF_CHECK_IDX, None)
    residue_ok.pop(_SELF_CHECK_IDX, None)
    n_aligned = max(0, n_aligned - 1)
    if self_mapped is None:
        raise ReferenceFetchError(
            f"panel {panel_name!r}: the anchor {residue_check.anchor_accession} failed to map "
            f"position {position} onto itself. The alignment walk is broken, so no residue "
            "check can be trusted."
        )
    if self_mapped != position:
        raise ReferenceFetchError(
            f"panel {panel_name!r}: anchor self-test failed — position {position} mapped onto "
            f"itself as {self_mapped}. This is an off-by-one in the alignment walk; every "
            "residue check would be shifted by the same amount."
        )
    if n_ok:
        n_ok -= 1 if anchor[position - 1].upper() in residue_check.accepted_residues else 0
    reporter.record(
        f"    anchor self-test passed: position {position} maps to itself in "
        f"{residue_check.anchor_accession}"
    )

    reporter.info(
        f"  residue anchor for panel {panel_name!r}: {len(positions)}/{len(refs)} references "
        f"mapped to position {position}; {n_ok} carry an accepted residue themselves "
        f"({n_ok / len(positions):.0%} of mapped)" if positions else
        f"  residue anchor for panel {panel_name!r}: no reference mapped"
    )
    if not positions:
        reporter.warn(
            f"panel {panel_name!r}: no reference could be aligned to anchor "
            f"{residue_check.anchor_accession}; the residue check will report as unavailable"
        )

    return AnchorMapping(
        panel=panel_name,
        anchor_accession=residue_check.anchor_accession,
        anchor_length=len(anchor),
        canonical_position=position,
        canonical_residue=canonical_residue,
        accepted_residues=tuple(sorted(residue_check.accepted_residues)),
        positions=positions,
        residue_ok=residue_ok,
        n_candidates=len(refs),
        n_aligned=n_aligned,
        n_position_spanned=len(positions),
        n_reference_residue_ok=n_ok,
    )


def _run(cmd: list[str], reporter: Reporter) -> None:
    reporter.record(f"    $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        raise DependencyError(
            f"command failed: {' '.join(cmd)}\n  stderr: {(exc.stderr or '').strip()}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise DependencyError(f"command timed out: {' '.join(cmd)} ({exc})") from exc
