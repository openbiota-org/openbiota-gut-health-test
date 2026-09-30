"""Hit parsing, mate collapse, decoy competition, and normalisation.

Three things happen here and the order matters.

1.  **Mate collapse.** The two mates of one DNA fragment are not independent
    observations; counting both inflates every result by up to 2x. Hits are
    keyed by query id across both files and the higher-scoring hit wins. On
    this dataset the read ids are byte-identical between R1 and R2, so a plain
    set union on the id suffices — but the suffix-stripping is implemented
    anyway so the tool is correct on data that does use ``/1`` and ``/2``.

2.  **Classification.** Every fragment carries exactly one best hit across the
    whole combined database. If that hit is a decoy, the fragment is a decoy
    rejection and can never also be a target acceptance. This is the primary
    specificity mechanism.

3.  **Normalisation.**

        copies_per_genome = (target_fragments / mean_target_ref_length_aa)
                          / (rpoB_fragments   / mean_rpoB_ref_length_aa)

    reported as copies per 100 bacterial genomes. Depth-independent,
    gene-length-corrected, and immune to human read contamination because
    human reads carry no bacterial rpoB and so inflate neither numerator nor
    denominator. If rpoB fragments are zero the result is *indeterminate* and
    is reported as such — never as a per-million-reads figure dressed up as an
    equivalent.
"""

from __future__ import annotations

import re
import statistics
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.confidence import classify_confidence
from openbiota.panels import ROLE_DECOY, ROLE_NORMALIZER, ROLE_TARGET, Panel, PanelSet
from openbiota.references import OUT_OF_SCOPE_SUFFIX, ReferenceDatabase
from openbiota.residues import (
    NO_MAPPING,
    NOT_SPANNED,
    SPAN_FAIL,
    SPAN_GAP,
    SPAN_PASS,
    classify_read_residue,
)

#: Trailing mate markers seen in the wild. This dataset has none — the ids are
#: identical in R1 and R2 — but stripping them is required for correctness on
#: any input that does carry them.
_MATE_SUFFIX_RE: Final = re.compile(r"/[12]$")

#: Minimum aligned residues for the rpoB denominator. Kept equal to the panel
#: default so numerator and denominator are filtered on the same footing.
NORMALIZER_MIN_ALIGNMENT_AA: Final = 25

#: Entries whose fragments are kept per individual reference accession, not
#: only per organism. One species can encode several enzymes of the same
#: family in different structural classes - Faecalibacterium prausnitzii
#: encodes both a Loop 1 beta-glucuronidase and the FMN-binding No Loop one -
#: so attributing a class by organism would be wrong for exactly the largest
#: contributors. The class belongs to the protein.
TRACK_ACCESSIONS: Final[frozenset[str]] = frozenset({"bglucuronidase:GUS"})

#: How many nearest reference organisms to retain per entry.
TOP_ORGANISMS: Final = 8


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #


def normalise_query_id(qseqid: str) -> str:
    """Collapse a read id to its fragment id.

    DIAMOND's ``qseqid`` is already the first whitespace-delimited token, so
    the Illumina ``1:N:0:INDEX`` comment field is gone. What remains is a
    possible ``/1`` or ``/2`` mate marker.
    """
    return _MATE_SUFFIX_RE.sub("", qseqid)


@dataclass(frozen=True, slots=True)
class HitRecord:
    qseqid: str
    sseqid: str
    pident: float
    aln_len: int
    sstart: int
    send: int
    bitscore: float
    qseq_gapped: str = ""
    sseq_gapped: str = ""


class HitParser:
    """Column-index-driven parser for ``diamond -f 6`` output."""

    __slots__ = ("_idx", "_n_min", "fields")

    def __init__(self, fields: Sequence[str]) -> None:
        self.fields = tuple(fields)
        required = ("qseqid", "sseqid", "pident", "length", "sstart", "send", "bitscore")
        missing = [f for f in required if f not in self.fields]
        if missing:
            raise ValueError(
                f"DIAMOND output field list is missing {missing}; got {list(self.fields)}"
            )
        self._idx = {name: i for i, name in enumerate(self.fields)}
        self._n_min = max(self._idx[f] for f in required) + 1

    @property
    def has_residue_fields(self) -> bool:
        return "qseq_gapped" in self._idx and "sseq_gapped" in self._idx

    def parse(self, line: str) -> HitRecord | None:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < self._n_min:
            return None
        i = self._idx
        try:
            return HitRecord(
                qseqid=parts[i["qseqid"]],
                sseqid=parts[i["sseqid"]],
                pident=float(parts[i["pident"]]),
                aln_len=int(parts[i["length"]]),
                sstart=int(parts[i["sstart"]]),
                send=int(parts[i["send"]]),
                bitscore=float(parts[i["bitscore"]]),
                qseq_gapped=parts[i["qseq_gapped"]] if "qseq_gapped" in i and len(parts) > i["qseq_gapped"] else "",
                sseq_gapped=parts[i["sseq_gapped"]] if "sseq_gapped" in i and len(parts) > i["sseq_gapped"] else "",
            )
        except (ValueError, IndexError):
            return None


def iter_lines(path: Path) -> Iterator[str]:
    """Stream a hits TSV line by line; never load it fully into memory."""
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        yield from fh


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class OrganismCount:
    organism: str
    genus: str
    phylum: str
    fragments: int
    mean_identity: float

    def to_json(self) -> dict[str, Any]:
        return {
            "organism": self.organism,
            "genus": self.genus,
            "phylum": self.phylum,
            "fragments": self.fragments,
            "mean_identity": round(self.mean_identity, 1),
        }


@dataclass(frozen=True, slots=True)
class ResidueSummary:
    """Outcome of the active-site residue check for one entry.

    Two independent halves of the same filter are reported.

    **Reference side** — every target reference was aligned to the anchor at
    database build time and its own residue at the canonical position recorded.
    A reference that does not carry an accepted residue is very likely a
    misannotated member of the homologous decoy family, so fragments whose best
    hit lands on one are separated out. This applies to *every* accepted
    fragment, which is what makes it useful at read level.

    **Read side** — the published filter, applied to those fragments whose
    alignment happens to span the mapped position. On 100-151 bp reads that is
    a small minority, which is precisely why the reference-side half exists.
    """

    available: bool
    anchor_accession: str
    canonical_position: int
    canonical_residue: str
    accepted_residues: tuple[str, ...]
    references_total: int
    references_mapped: int
    references_residue_ok: int
    fragments_checked: int
    fragments_on_mapped_reference: int
    fragments_on_residue_ok_reference: int
    fragments_on_residue_bad_reference: int
    fragments_spanning: int
    passed: int
    failed: int
    gapped: int
    not_spanned: int
    no_mapping: int
    observed: dict[str, int] = field(default_factory=dict)

    @property
    def pass_rate(self) -> float | None:
        """Read-side pass rate over the fragments that could be evaluated."""
        evaluable = self.passed + self.failed
        return self.passed / evaluable if evaluable else None

    @property
    def reference_consistent_fraction(self) -> float | None:
        """Fraction of accepted fragments landing on a residue-consistent reference."""
        evaluable = self.fragments_on_residue_ok_reference + self.fragments_on_residue_bad_reference
        return self.fragments_on_residue_ok_reference / evaluable if evaluable else None

    def to_json(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "anchor_accession": self.anchor_accession,
            "canonical_position": self.canonical_position,
            "canonical_residue": self.canonical_residue,
            "accepted_residues": list(self.accepted_residues),
            "reference_side": {
                "references_total": self.references_total,
                "references_mapped": self.references_mapped,
                "references_carrying_accepted_residue": self.references_residue_ok,
                "fragments_on_residue_consistent_reference": (
                    self.fragments_on_residue_ok_reference
                ),
                "fragments_on_residue_inconsistent_reference": (
                    self.fragments_on_residue_bad_reference
                ),
                "fragments_on_unmapped_reference": (
                    self.fragments_checked
                    - self.fragments_on_residue_ok_reference
                    - self.fragments_on_residue_bad_reference
                ),
                "consistent_fraction": (
                    None
                    if self.reference_consistent_fraction is None
                    else round(self.reference_consistent_fraction, 4)
                ),
            },
            "read_side": {
                "fragments_checked": self.fragments_checked,
                "fragments_spanning_position": self.fragments_spanning,
                "passed": self.passed,
                "failed": self.failed,
                "gapped": self.gapped,
                "not_spanned": self.not_spanned,
                "no_mapping": self.no_mapping,
                "observed_residues": dict(sorted(self.observed.items())),
                "pass_rate": None if self.pass_rate is None else round(self.pass_rate, 4),
            },
        }


@dataclass(frozen=True, slots=True)
class EntryResult:
    """Per-gene (or per-decoy-family) counts and normalised figures."""

    key: str
    panel: str
    entry_id: str
    label: str
    gene: str | None
    role: str
    fragments: int
    rejected_low_identity: int
    rejected_short_alignment: int
    mean_identity: float | None
    median_identity: float | None
    reference_count: int
    reference_mean_length_aa: float
    reference_truncated: bool
    min_identity: float
    copies_per_100_genomes: float | None
    #: Same figure restricted to fragments whose best-hit reference itself
    #: carries an accepted active-site residue. Present only for entries with
    #: a residue check; a tighter estimate than the headline upper bound.
    copies_per_100_genomes_residue_consistent: float | None = None
    organisms: tuple[OrganismCount, ...] = ()
    residue: ResidueSummary | None = None
    #: Fragments per individual reference accession. Retained only for the
    #: entries in ``TRACK_ACCESSIONS``, where the reading is a mixture of
    #: structurally distinct enzymes and the mixture is the point: a
    #: beta-glucuronidase total says nothing about which substrates the
    #: enzymes act on, and the structural class is a property of the
    #: reference protein rather than of the read that hit it.
    accession_fragments: tuple[tuple[str, int], ...] = ()

    @property
    def stable(self) -> bool:
        return self.fragments > 0

    @property
    def confidence(self) -> str:
        return classify_confidence(
            fragments=self.fragments,
            mean_identity=self.mean_identity,
            indeterminate=self.copies_per_100_genomes is None and self.fragments > 0,
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "confidence": self.confidence,
            "panel": self.panel,
            "entry_id": self.entry_id,
            "label": self.label,
            "gene": self.gene,
            "role": self.role,
            "fragments": self.fragments,
            "rejected_low_identity": self.rejected_low_identity,
            "rejected_short_alignment": self.rejected_short_alignment,
            "mean_identity": None if self.mean_identity is None else round(self.mean_identity, 1),
            "median_identity": None if self.median_identity is None else round(self.median_identity, 1),
            "reference_count": self.reference_count,
            "reference_mean_length_aa": round(self.reference_mean_length_aa, 1),
            "reference_truncated": self.reference_truncated,
            "min_identity": self.min_identity,
            "copies_per_100_genomes": (
                None if self.copies_per_100_genomes is None else round(self.copies_per_100_genomes, 4)
            ),
            "copies_per_100_genomes_residue_consistent": (
                None
                if self.copies_per_100_genomes_residue_consistent is None
                else round(self.copies_per_100_genomes_residue_consistent, 4)
            ),
            "organisms": [o.to_json() for o in self.organisms],
            "residue_check": None if self.residue is None else self.residue.to_json(),
            "accession_fragments": (
                [list(pair) for pair in self.accession_fragments]
                if self.accession_fragments else None
            ),
        }


@dataclass(frozen=True, slots=True)
class NormalizerResult:
    fragments: int
    rejected_low_identity: int
    rejected_short_alignment: int
    reference_count: int
    reference_mean_length_aa: float
    mean_identity: float | None
    phylum_counts: dict[str, int] = field(default_factory=dict)
    genus_counts: dict[str, int] = field(default_factory=dict)
    organism_counts: dict[str, int] = field(default_factory=dict)

    @property
    def indeterminate(self) -> bool:
        return self.fragments == 0

    @property
    def rate(self) -> float | None:
        """Fragments per amino acid of reference length."""
        if self.fragments == 0 or self.reference_mean_length_aa <= 0:
            return None
        return self.fragments / self.reference_mean_length_aa

    def to_json(self) -> dict[str, Any]:
        return {
            "gene": "rpoB",
            "fragments": self.fragments,
            "rejected_low_identity": self.rejected_low_identity,
            "rejected_short_alignment": self.rejected_short_alignment,
            "reference_count": self.reference_count,
            "reference_mean_length_aa": round(self.reference_mean_length_aa, 1),
            "mean_identity": None if self.mean_identity is None else round(self.mean_identity, 1),
            "indeterminate": self.indeterminate,
            "rate_fragments_per_aa": None if self.rate is None else round(self.rate, 6),
        }


@dataclass(frozen=True, slots=True)
class PanelResult:
    panel: Panel
    targets: tuple[EntryResult, ...]
    decoys: tuple[EntryResult, ...]
    accepted_fragments: int
    decoy_fragments: int
    rejected_low_identity: int
    rejected_short_alignment: int
    copies_per_100_genomes: float | None
    aggregate_method: str
    indeterminate: bool
    #: Same aggregate restricted to residue-consistent references, where the
    #: panel has a residue check. ``None`` when the panel has none.
    copies_per_100_genomes_residue_consistent: float | None = None
    residue_consistent_fragments: int | None = None

    @property
    def aggregate_fragments(self) -> int:
        """Fragments on the genes that actually set the headline figure.

        ``accepted_fragments`` is the whole panel, context genes included.
        For a panel whose aggregate is one gene, that total can be large
        while the gene itself has nothing: the histamine panel aggregates
        over the decarboxylase HDCA, and a sample can carry thousands of
        fragments of the HDCTRANS transporter with zero HDCA. The total is
        the right number to print as "fragments matched"; it is the wrong
        number to decide whether the headline measurement exists.
        """
        aggregate_ids = set(self.panel.aggregate_target_ids())
        return sum(
            t.fragments for t in self.targets if t.entry_id in aggregate_ids
        )

    @property
    def stable(self) -> bool:
        return self.aggregate_fragments >= self.panel.min_fragments_for_stability

    @property
    def detected(self) -> bool:
        """Whether the headline measurement exists at all.

        Restricted to the aggregate genes so this agrees with
        :attr:`confidence`, which is already restricted to them. They
        disagreed: a panel could report ``confidence="not detected"`` in
        results.json while the report, reading ``detected``, placed it on
        the reference scale and called it "notably low" - a judged reading
        for a gene with zero fragments.
        """
        return self.aggregate_fragments > 0

    @property
    def confidence(self) -> str:
        """Tier of the genes that actually set this panel's headline figure."""
        aggregate_ids = set(self.panel.aggregate_target_ids())
        tiers = [t.confidence for t in self.targets if t.entry_id in aggregate_ids]
        if self.indeterminate:
            return "indeterminate"
        if any(t == "confirmed" for t in tiers):
            return "confirmed"
        if any(t == "provisional" for t in tiers):
            return "provisional"
        return "not detected"

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.panel.name,
            "confidence": self.confidence,
            "description": " ".join(self.panel.description.split()),
            "metabolite": self.panel.metabolite,
            "pathway": self.panel.pathway,
            "extension": self.panel.extension,
            "citation": " ".join(self.panel.citation.split()),
            "accepted_fragments": self.accepted_fragments,
            "decoy_fragments": self.decoy_fragments,
            "rejected_low_identity": self.rejected_low_identity,
            "rejected_short_alignment": self.rejected_short_alignment,
            "copies_per_100_genomes": (
                None if self.copies_per_100_genomes is None else round(self.copies_per_100_genomes, 4)
            ),
            "copies_per_100_genomes_residue_consistent": (
                None
                if self.copies_per_100_genomes_residue_consistent is None
                else round(self.copies_per_100_genomes_residue_consistent, 4)
            ),
            "residue_consistent_fragments": self.residue_consistent_fragments,
            "aggregate_method": self.aggregate_method,
            "aggregate_from": list(self.panel.aggregate_target_ids()),
            "aggregate_reason": " ".join(self.panel.aggregate_reason.split()),
            "indeterminate": self.indeterminate,
            "stable": self.stable,
            "min_fragments_for_stability": self.panel.min_fragments_for_stability,
            "upper_bound_reason": (
                " ".join(self.panel.upper_bound_reason.split())
                if self.panel.upper_bound_reason
                else None
            ),
            "caveats": [" ".join(c.split()) for c in self.panel.caveats],
            "known_producers": list(self.panel.organisms),
            "genes": [t.to_json() for t in self.targets],
            "decoys": [d.to_json() for d in self.decoys],
        }


@dataclass(frozen=True, slots=True)
class TallyResult:
    total_hit_lines: int
    fragments_with_hit: int
    unresolved_sseqids: int
    normalizer: NormalizerResult
    panels: tuple[PanelResult, ...]
    hits_per_mate: dict[str, int] = field(default_factory=dict)

    def panel(self, name: str) -> PanelResult:
        for p in self.panels:
            if p.panel.name == name:
                return p
        raise KeyError(name)

    def to_json(self) -> dict[str, Any]:
        return {
            "total_hit_lines": self.total_hit_lines,
            "fragments_with_hit": self.fragments_with_hit,
            "unresolved_sseqids": self.unresolved_sseqids,
            "hits_per_mate": dict(self.hits_per_mate),
            "normalizer": self.normalizer.to_json(),
            "panels": [p.to_json() for p in self.panels],
        }


# --------------------------------------------------------------------------- #
# best-hit collapse
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class _BestHit:
    idx: int
    pident: float
    aln_len: int
    sstart: int
    send: int
    bitscore: float
    qseq_gapped: str
    sseq_gapped: str


def collapse_mates(
    streams: Iterable[tuple[str, Iterable[str]]],
    database: ReferenceDatabase,
    parser: HitParser,
) -> tuple[dict[str, _BestHit], int, int, dict[str, int]]:
    """Reduce per-read hits to one best hit per DNA fragment.

    Returns ``(best_by_fragment, total_lines, unresolved_sseqids, lines_per_mate)``.

    Memory is bounded by the number of *hits*, not by the number of reads: only
    reads that aligned to something in the reference database are held, which
    on a 16 M read pair sample is order 10^5 entries.
    """
    best: dict[str, _BestHit] = {}
    total = 0
    unresolved = 0
    per_mate: dict[str, int] = {}

    for mate_label, lines in streams:
        mate_total = 0
        for line in lines:
            record = parser.parse(line)
            if record is None:
                continue
            mate_total += 1
            reference = database.by_sseqid(record.sseqid)
            if reference is None:
                unresolved += 1
                continue
            fragment = normalise_query_id(record.qseqid)
            current = best.get(fragment)
            if current is not None and current.bitscore >= record.bitscore:
                continue
            # The gapped alignment strings are only retained where a residue
            # check can actually use them, which keeps this dictionary small.
            keep_alignment = reference.anchor_pos is not None
            best[fragment] = _BestHit(
                idx=reference.idx,
                pident=record.pident,
                aln_len=record.aln_len,
                sstart=record.sstart,
                send=record.send,
                bitscore=record.bitscore,
                qseq_gapped=record.qseq_gapped if keep_alignment else "",
                sseq_gapped=record.sseq_gapped if keep_alignment else "",
            )
        per_mate[mate_label] = mate_total
        total += mate_total

    return best, total, unresolved, per_mate


# --------------------------------------------------------------------------- #
# accumulation
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class _EntryAcc:
    fragments: int = 0
    rejected_low_identity: int = 0
    rejected_short_alignment: int = 0
    identities: list[float] = field(default_factory=list)
    organisms: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    organism_identity: dict[tuple[str, str, str], list[float]] = field(
        default_factory=lambda: defaultdict(list)
    )
    accessions: Counter[str] = field(default_factory=Counter)
    residue_outcomes: Counter[str] = field(default_factory=Counter)
    residue_observed: Counter[str] = field(default_factory=Counter)
    fragments_on_mapped_reference: int = 0
    fragments_on_residue_ok_reference: int = 0
    fragments_on_residue_bad_reference: int = 0


def _aggregate(values: Sequence[float], method: str) -> float | None:
    if not values:
        return None
    if method == "sum":
        return float(sum(values))
    if method == "median":
        return float(statistics.median(values))
    if method == "max":
        return float(max(values))
    if method == "mean":
        return float(statistics.fmean(values))
    raise ValueError(f"unknown aggregate method {method!r}")


def tally(
    *,
    streams: Iterable[tuple[str, Iterable[str]]],
    database: ReferenceDatabase,
    panel_set: PanelSet,
    fields: Sequence[str],
) -> TallyResult:
    """Turn streamed DIAMOND output into normalised per-panel results."""
    parser = HitParser(fields)
    best, total_lines, unresolved, per_mate = collapse_mates(streams, database, parser)

    panels_by_name = {p.name: p for p in panel_set.panels}
    entry_min_identity = {e.key: e.min_identity for e in panel_set.all_entries()}
    # A scoped target's shadow decoy inherits its identity floor: the demoted
    # references are the same gene family, so the same threshold applies.
    entry_min_identity.update({
        f"{e.key}{OUT_OF_SCOPE_SUFFIX}": e.min_identity for e in panel_set.all_entries()
    })
    accumulators: dict[str, _EntryAcc] = defaultdict(_EntryAcc)

    normalizer_key = panel_set.normalizer.key
    normalizer_min_identity = panel_set.normalizer.min_identity

    for hit in best.values():
        reference = database.sequences[hit.idx]
        key = reference.entry_key
        acc = accumulators[key]
        taxon = (reference.organism, reference.genus, reference.phylum)

        if reference.role == ROLE_NORMALIZER:
            if hit.aln_len < NORMALIZER_MIN_ALIGNMENT_AA:
                acc.rejected_short_alignment += 1
                continue
            if hit.pident < normalizer_min_identity:
                acc.rejected_low_identity += 1
                continue
            acc.fragments += 1
            acc.identities.append(hit.pident)
            acc.organisms[taxon] += 1
            acc.organism_identity[taxon].append(hit.pident)
            continue

        panel = panels_by_name.get(reference.panel)
        min_alignment = panel.min_alignment_aa if panel else 25
        if hit.aln_len < min_alignment:
            acc.rejected_short_alignment += 1
            continue
        if hit.pident < entry_min_identity.get(key, 50.0):
            acc.rejected_low_identity += 1
            continue

        acc.fragments += 1
        acc.identities.append(hit.pident)
        acc.organisms[taxon] += 1
        acc.organism_identity[taxon].append(hit.pident)
        if key in TRACK_ACCESSIONS:
            acc.accessions[reference.accession] += 1

        if (
            reference.role == ROLE_TARGET
            and panel is not None
            and panel.residue_check is not None
            and panel.residue_check.enabled
            and reference.entry_id in panel.residue_check.applies_to
        ):
            if reference.anchor_pos is not None:
                acc.fragments_on_mapped_reference += 1
                if reference.anchor_residue_ok:
                    acc.fragments_on_residue_ok_reference += 1
                else:
                    acc.fragments_on_residue_bad_reference += 1
            outcome, residue = classify_read_residue(
                anchor_position=reference.anchor_pos,
                subject_start=hit.sstart,
                subject_end=hit.send,
                query_gapped=hit.qseq_gapped,
                subject_gapped=hit.sseq_gapped,
                accepted=panel.residue_check.accepted_residues,
            )
            acc.residue_outcomes[outcome] += 1
            if residue:
                acc.residue_observed[residue] += 1

    # ---- normaliser -------------------------------------------------------- #
    norm_acc = accumulators.get(normalizer_key, _EntryAcc())
    norm_stats = database.entries.get(normalizer_key)
    if norm_stats is None:
        raise ValueError(
            f"reference database has no entry {normalizer_key!r}; rebuild it with --rebuild-db"
        )
    normalizer = NormalizerResult(
        fragments=norm_acc.fragments,
        rejected_low_identity=norm_acc.rejected_low_identity,
        rejected_short_alignment=norm_acc.rejected_short_alignment,
        reference_count=norm_stats.n_sequences,
        reference_mean_length_aa=norm_stats.mean_length_aa,
        mean_identity=statistics.fmean(norm_acc.identities) if norm_acc.identities else None,
        phylum_counts=_counter_by(norm_acc.organisms, 2),
        genus_counts=_counter_by(norm_acc.organisms, 1),
        organism_counts=_counter_by(norm_acc.organisms, 0),
    )
    norm_rate = normalizer.rate

    # invariant: organism attribution must reconcile with the fragment total,
    # not with the read total. This was a real bug in the prototype.
    _assert_reconciles(norm_acc, normalizer_key)

    # ---- panels ------------------------------------------------------------ #
    panel_results: list[PanelResult] = []
    for panel in panel_set.panels:
        targets = tuple(
            _entry_result(
                entry_key=entry.key,
                panel=panel,
                database=database,
                acc=accumulators.get(entry.key, _EntryAcc()),
                norm_rate=norm_rate,
                role=ROLE_TARGET,
            )
            for entry in panel.targets
        )
        # Declared decoys, plus the shadow decoy a scoped target spawns for the
        # references it may not claim. Those fragments are real and are reported
        # as background rather than silently dropped: they are the measure of
        # how much of this gene family sits in organisms that carry the fold
        # without running the pathway.
        decoy_keys = [entry.key for entry in panel.decoys]
        decoy_keys += [
            f"{entry.key}{OUT_OF_SCOPE_SUFFIX}"
            for entry in panel.targets
            if f"{entry.key}{OUT_OF_SCOPE_SUFFIX}" in database.entries
        ]
        decoys = tuple(
            _entry_result(
                entry_key=key,
                panel=panel,
                database=database,
                acc=accumulators.get(key, _EntryAcc()),
                norm_rate=norm_rate,
                role=ROLE_DECOY,
            )
            for key in decoy_keys
        )
        # Every target is reported, but the headline aggregate may be narrowed
        # to the genes whose reference sets are specific enough to carry it.
        aggregate_ids = set(panel.aggregate_target_ids())
        aggregated = [t for t in targets if t.entry_id in aggregate_ids]
        values = [
            t.copies_per_100_genomes for t in aggregated if t.copies_per_100_genomes is not None
        ]
        consistent_values = [
            t.copies_per_100_genomes_residue_consistent
            for t in aggregated
            if t.copies_per_100_genomes_residue_consistent is not None
        ]
        consistent_fragments = [
            t.residue.fragments_on_residue_ok_reference for t in targets if t.residue is not None
        ]
        panel_results.append(
            PanelResult(
                panel=panel,
                targets=targets,
                decoys=decoys,
                accepted_fragments=sum(t.fragments for t in targets),
                decoy_fragments=sum(d.fragments for d in decoys),
                rejected_low_identity=sum(t.rejected_low_identity for t in targets),
                rejected_short_alignment=sum(t.rejected_short_alignment for t in targets),
                copies_per_100_genomes=_aggregate(values, panel.aggregate),
                aggregate_method=panel.aggregate,
                indeterminate=norm_rate is None,
                copies_per_100_genomes_residue_consistent=(
                    _aggregate(consistent_values, panel.aggregate) if consistent_values else None
                ),
                residue_consistent_fragments=(
                    sum(consistent_fragments) if consistent_fragments else None
                ),
            )
        )

    return TallyResult(
        total_hit_lines=total_lines,
        fragments_with_hit=len(best),
        unresolved_sseqids=unresolved,
        normalizer=normalizer,
        panels=tuple(panel_results),
        hits_per_mate=per_mate,
    )


def _counter_by(counter: Counter[tuple[str, str, str]], position: int) -> dict[str, int]:
    out: Counter[str] = Counter()
    for taxon, count in counter.items():
        label = taxon[position] or "unassigned"
        out[label] += count
    return dict(out.most_common())


def _assert_reconciles(acc: _EntryAcc, key: str) -> None:
    attributed = sum(acc.organisms.values())
    if attributed != acc.fragments:
        raise AssertionError(
            f"{key}: organism attribution sums to {attributed} but {acc.fragments} fragments "
            "were accepted. Attribution must be derived from mate-collapsed fragments."
        )


def _entry_result(
    *,
    entry_key: str,
    panel: Panel,
    database: ReferenceDatabase,
    acc: _EntryAcc,
    norm_rate: float | None,
    role: str,
) -> EntryResult:
    stats = database.entries.get(entry_key)
    if stats is None:
        raise ValueError(
            f"reference database has no entry {entry_key!r}; the panel definition changed since "
            "the database was built — re-run with --rebuild-db"
        )

    copies: float | None = None
    if norm_rate is not None and stats.mean_length_aa > 0:
        rate = acc.fragments / stats.mean_length_aa
        copies = (rate / norm_rate) * 100.0

    organisms = tuple(
        OrganismCount(
            organism=taxon[0],
            genus=taxon[1] or "unassigned",
            phylum=taxon[2] or "unassigned",
            fragments=count,
            mean_identity=statistics.fmean(acc.organism_identity[taxon]),
        )
        for taxon, count in acc.organisms.most_common(TOP_ORGANISMS)
    )

    residue: ResidueSummary | None = None
    copies_consistent: float | None = None
    rc = panel.residue_check
    if role == ROLE_TARGET and rc is not None and rc.enabled and stats.entry_id in rc.applies_to:
        mapping = database.anchors.get(panel.name)
        # Normalise the residue-consistent subset against the mean length of
        # the residue-consistent references, falling back to the entry mean.
        consistent_length = stats.mean_length_residue_ok_aa or stats.mean_length_aa
        if norm_rate is not None and consistent_length > 0:
            consistent_rate = acc.fragments_on_residue_ok_reference / consistent_length
            copies_consistent = (consistent_rate / norm_rate) * 100.0
        residue = ResidueSummary(
            available=mapping is not None and bool(mapping.positions),
            anchor_accession=rc.anchor_accession,
            canonical_position=rc.canonical_position,
            canonical_residue=mapping.canonical_residue if mapping else "?",
            accepted_residues=tuple(sorted(rc.accepted_residues)),
            references_total=stats.n_sequences,
            references_mapped=stats.n_anchor_mapped,
            references_residue_ok=stats.n_anchor_residue_ok,
            fragments_checked=acc.fragments,
            fragments_on_mapped_reference=acc.fragments_on_mapped_reference,
            fragments_on_residue_ok_reference=acc.fragments_on_residue_ok_reference,
            fragments_on_residue_bad_reference=acc.fragments_on_residue_bad_reference,
            fragments_spanning=(
                acc.residue_outcomes[SPAN_PASS]
                + acc.residue_outcomes[SPAN_FAIL]
                + acc.residue_outcomes[SPAN_GAP]
            ),
            passed=acc.residue_outcomes[SPAN_PASS],
            failed=acc.residue_outcomes[SPAN_FAIL],
            gapped=acc.residue_outcomes[SPAN_GAP],
            not_spanned=acc.residue_outcomes[NOT_SPANNED],
            no_mapping=acc.residue_outcomes[NO_MAPPING],
            observed=dict(acc.residue_observed),
        )
        # Sanity: the five outcome codes must partition every accepted fragment.
        accounted = sum(acc.residue_outcomes.values())
        if accounted != acc.fragments:
            raise AssertionError(
                f"{entry_key}: residue outcomes account for {accounted} of {acc.fragments} "
                "accepted fragments"
            )

    return EntryResult(
        key=entry_key,
        panel=panel.name,
        entry_id=stats.entry_id,
        label=stats.label,
        gene=stats.gene,
        role=role,
        fragments=acc.fragments,
        rejected_low_identity=acc.rejected_low_identity,
        rejected_short_alignment=acc.rejected_short_alignment,
        mean_identity=statistics.fmean(acc.identities) if acc.identities else None,
        median_identity=float(statistics.median(acc.identities)) if acc.identities else None,
        reference_count=stats.n_sequences,
        reference_mean_length_aa=stats.mean_length_aa,
        reference_truncated=stats.truncated,
        min_identity=stats.min_identity,
        copies_per_100_genomes=copies,
        copies_per_100_genomes_residue_consistent=copies_consistent,
        organisms=organisms,
        residue=residue,
        accession_fragments=tuple(sorted(acc.accessions.items())),
    )


# --------------------------------------------------------------------------- #
# per-panel hit views
# --------------------------------------------------------------------------- #


def write_panel_hit_views(
    *,
    hit_files: Sequence[tuple[str, Path]],
    database: ReferenceDatabase,
    panel_names: Sequence[str],
    fields: Sequence[str],
    out_dir: Path,
) -> list[Path]:
    """Split the combined hits TSV into ``alignments/hits_<panel>_<mate>.tsv`` views.

    The combined file is the raw DIAMOND output; these are filtered views of it,
    written with a header so they are readable on their own.
    """
    parser = HitParser(fields)
    wanted = set(panel_names)
    written: list[Path] = []

    for mate_label, path in hit_files:
        handles: dict[str, Any] = {}
        try:
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    record = parser.parse(line)
                    if record is None:
                        continue
                    reference = database.by_sseqid(record.sseqid)
                    if reference is None or reference.panel not in wanted:
                        continue
                    handle = handles.get(reference.panel)
                    if handle is None:
                        target = out_dir / "alignments" / f"hits_{reference.panel}_{mate_label}.tsv"
                        target.parent.mkdir(parents=True, exist_ok=True)
                        handle = target.open("w", encoding="utf-8")
                        handle.write("#" + "\t".join([*fields, "role", "entry", "organism"]) + "\n")
                        handles[reference.panel] = handle
                        written.append(target)
                    handle.write(
                        line.rstrip("\n")
                        + f"\t{reference.role}\t{reference.entry_id}\t{reference.organism}\n"
                    )
        finally:
            for handle in handles.values():
                handle.close()
    return written
