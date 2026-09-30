"""Reference sequence acquisition and DIAMOND database construction.

One combined protein database holds every panel's targets, every panel's
decoys, and the rpoB normaliser. That is deliberate:

* the whole FASTQ is translated once instead of once per panel, which is where
  essentially all the runtime goes;
* best-hit competition becomes global, so a read is only credited to a target
  when that target beats every decoy *and* every other panel's target. cutC
  versus hpdB — both glycyl-radical enzymes, and both targets of different
  panels — is resolved for free by this arrangement.

Everything is cached. A re-run with unchanged panels re-downloads nothing and
rebuilds nothing.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import statistics
import subprocess
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota import __version__
from openbiota.errors import DependencyError, ReferenceFetchError
from openbiota.logging_util import Reporter, human_bytes
from openbiota.net import build_url, get, next_page_url
from openbiota.panels import PanelSet, RefEntry
from openbiota.residues import AnchorMapping, map_anchor_positions
from openbiota.seqio import clean_sequence

UNIPROT_SEARCH: Final = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS: Final = "accession,protein_name,gene_primary,organism_name,lineage,length,sequence"
UNIPROT_PAGE_SIZE: Final = 500

#: FASTA header field separator. Must not appear in UniProt accessions, in
#: panel entry ids (validated in panels.py), or be special to DIAMOND's seqid
#: parser — verified against DIAMOND 2.2.6.
SEP: Final = "~"

#: Sequences shorter than this are unusable as blastx subjects for 33-50
#: residue query translations.
MIN_REF_LENGTH_AA: Final = 60

#: On-disk index schema. Bumping this rebuilds the DIAMOND database.
INDEX_VERSION: Final = 5

#: Appended to a scoped target's entry id for the references it may not claim.
#: Must avoid ``SEP`` and whitespace (the FASTA header encoding) and stay
#: inside the panel id charset (alphanumeric plus ``-`` and ``.``).
OUT_OF_SCOPE_SUFFIX: Final = ".OOS"

#: UniProt fetch cache schema, kept separate from INDEX_VERSION so that
#: changing how the index is *stored* does not throw away a slow download.
FETCH_CACHE_VERSION: Final = 3

#: Ranks pulled out of the UniProt "Taxonomic lineage" string.
_WANTED_RANKS: Final = ("phylum", "class", "order", "family", "genus")


@dataclass(frozen=True, slots=True)
class RefSequence:
    """One reference protein in the combined database."""

    idx: int
    entry_key: str
    entry_id: str
    role: str
    panel: str
    accession: str
    organism: str
    length: int
    phylum: str = ""
    genus: str = ""
    family: str = ""
    #: 1-based position in this sequence equivalent to the panel's canonical
    #: active-site position, or ``None`` if it could not be mapped.
    anchor_pos: int | None = None
    #: Whether this reference's own residue at ``anchor_pos`` is one the panel
    #: accepts. ``None`` means unmapped. A reference that fails this is very
    #: likely a misannotated member of the homologous decoy family, so reads
    #: whose best hit lands on one are reported as a separate, weaker subset.
    anchor_residue_ok: bool | None = None

    @property
    def sseqid(self) -> str:
        return f"{self.idx}{SEP}{self.entry_id}{SEP}{self.accession}"


@dataclass(frozen=True, slots=True)
class EntryStats:
    """Per reference-entry provenance, reported verbatim in run.log."""

    key: str
    panel: str
    entry_id: str
    role: str
    label: str
    gene: str | None
    source: str
    query: str
    min_identity: float
    total_available: int
    n_fetched: int
    n_dropped_length: int
    n_dropped_identical: int
    n_dropped_claimed: int
    n_sequences: int
    mean_length_aa: float
    median_length_aa: float
    min_length_aa: int
    max_length_aa: int
    length_window: str
    n_anchor_mapped: int = 0
    n_anchor_residue_ok: int = 0
    mean_length_residue_ok_aa: float = 0.0
    truncated: bool = False
    #: set on the shadow decoy entry a scoped target spawns (see Scope)
    scope_reason: str = ""
    #: set on the target itself: how many of its references were demoted
    n_out_of_scope: int = 0

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ReferenceDatabase:
    """A built, on-disk DIAMOND database plus the metadata needed to read hits."""

    fingerprint: str
    dmnd_path: Path
    fasta_path: Path
    index_path: Path
    sequences: list[RefSequence]
    entries: dict[str, EntryStats]
    diamond_version: str
    built_at: str
    anchors: dict[str, AnchorMapping] = field(default_factory=dict)

    # -- lookups ---------------------------------------------------------- #

    def by_sseqid(self, sseqid: str) -> RefSequence | None:
        """Resolve a DIAMOND ``sseqid`` back to its reference metadata.

        The leading field is the database index, so this is an O(1) list
        lookup rather than a dictionary of tens of thousands of strings.
        """
        head, _, _ = sseqid.partition(SEP)
        try:
            idx = int(head)
        except ValueError:
            return None
        if 0 <= idx < len(self.sequences):
            return self.sequences[idx]
        return None

    def entry(self, key: str) -> EntryStats:
        return self.entries[key]

    @property
    def n_sequences(self) -> int:
        return len(self.sequences)

    def total_residues(self) -> int:
        return sum(s.length for s in self.sequences)


# --------------------------------------------------------------------------- #
# UniProt fetch
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class FetchedRecord:
    accession: str
    protein_name: str
    gene: str
    organism: str
    lineage: str
    length: int
    sequence: str


#: Where `literal` entries are read from, relative to the repository root.
LITERAL_ROOT: Final = Path(__file__).resolve().parents[1]


def _read_literal(entry: Any) -> list[FetchedRecord]:
    """Sequences shipped with this repository, addressed as `file.json#key`.

    For a protein with no UniProt entry the alternative is a gene-symbol
    query, and for this operon that returns a fungal protein which happens
    to share the symbol. A named accession and the sequence beside it is
    both more specific and auditable.
    """
    spec = str(entry.query)
    if "#" not in spec:
        raise ReferenceFetchError(
            f"{entry.key}: a literal source must be addressed as 'file.json#key', "
            f"got {spec!r}"
        )
    rel, _, key = spec.partition("#")
    path = LITERAL_ROOT / rel
    if not path.is_file():
        raise ReferenceFetchError(f"{entry.key}: {path} does not exist")
    blob = json.loads(path.read_text(encoding="utf-8"))
    proteins = blob.get("proteins") or {}
    wanted = [key] if key != "*" else sorted(proteins)
    records: list[FetchedRecord] = []
    for name in wanted:
        record = proteins.get(name)
        if record is None:
            raise ReferenceFetchError(
                f"{entry.key}: {rel} has no protein {name!r}; it holds "
                f"{sorted(proteins)}"
            )
        sequence = clean_sequence(str(record["sequence"]))
        if not sequence:
            raise ReferenceFetchError(f"{entry.key}: {name} has an empty sequence")
        records.append(FetchedRecord(
            accession=str(record.get("protein_id") or name),
            protein_name=str(blob.get("source") or name),
            gene=name,
            organism=str(blob.get("organism") or blob.get("source") or "unknown organism"),
            lineage="", length=len(sequence), sequence=sequence,
        ))
    if not records:
        raise ReferenceFetchError(
            f"{entry.key}: {spec} resolved to no sequences, which would silently become "
            "a count of zero"
        )
    return records


def _parse_lineage(lineage: str) -> dict[str, str]:
    """Parse UniProt's ``name (rank), name (rank), ...`` lineage string."""
    out: dict[str, str] = {}
    for chunk in lineage.split(","):
        chunk = chunk.strip()
        if not chunk.endswith(")") or "(" not in chunk:
            continue
        name, _, rank = chunk[:-1].rpartition("(")
        rank = rank.strip().lower()
        if rank in _WANTED_RANKS:
            out.setdefault(rank, name.strip())
    return out


def _entry_cache_key(entry: RefEntry, fetch_limit: int) -> str:
    payload = f"{entry.source}|{entry.query}|{UNIPROT_FIELDS}|{fetch_limit}|v{FETCH_CACHE_VERSION}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _fetch_uniprot_tsv(query: str, fetch_limit: int, reporter: Reporter) -> tuple[str, int]:
    """Return the concatenated TSV body and the server's total result count."""
    url: str | None = build_url(
        UNIPROT_SEARCH,
        {"query": query, "format": "tsv", "fields": UNIPROT_FIELDS, "size": UNIPROT_PAGE_SIZE},
    )
    chunks: list[str] = []
    rows = 0
    total = -1
    header: str | None = None

    while url and rows < fetch_limit:
        response = get(url)
        if total < 0:
            raw_total = response.header("x-total-results")
            total = int(raw_total) if raw_total and raw_total.isdigit() else -1
        text = response.text()
        lines = text.split("\n")
        if lines and lines[0].startswith("Entry\t"):
            if header is None:
                header = lines[0]
            lines = lines[1:]
        body = [ln for ln in lines if ln.strip()]
        rows += len(body)
        chunks.extend(body)
        url = next_page_url(response)
        if url:
            reporter.record(f"    paging: {rows} rows so far (total available {total})")

    if header is None:
        raise ReferenceFetchError(
            f"UniProt returned no TSV header for query {query!r}; the query is probably malformed"
        )
    return "\n".join([header, *chunks[:fetch_limit]]), total


def _parse_tsv(text: str) -> list[FetchedRecord]:
    lines = text.split("\n")
    if not lines:
        return []
    columns = lines[0].split("\t")
    try:
        i_acc = columns.index("Entry")
        i_name = columns.index("Protein names")
        i_gene = columns.index("Gene Names (primary)")
        i_org = columns.index("Organism")
        i_lin = columns.index("Taxonomic lineage")
        i_len = columns.index("Length")
        i_seq = columns.index("Sequence")
    except ValueError as exc:
        raise ReferenceFetchError(
            f"unexpected UniProt TSV columns {columns!r}; the REST field names may have changed"
        ) from exc

    out: list[FetchedRecord] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= i_seq:
            continue
        try:
            length = int(parts[i_len])
        except ValueError:
            continue
        sequence = parts[i_seq].strip().upper()
        if not sequence:
            continue
        out.append(
            FetchedRecord(
                accession=parts[i_acc].strip(),
                protein_name=parts[i_name].strip(),
                gene=parts[i_gene].strip(),
                organism=parts[i_org].strip(),
                lineage=parts[i_lin].strip(),
                length=length,
                sequence=sequence,
            )
        )
    return out


def fetch_entry_records(
    entry: RefEntry,
    cache_dir: Path,
    reporter: Reporter,
    *,
    refresh: bool = False,
) -> tuple[list[FetchedRecord], int]:
    """Fetch (or read from cache) the raw records for one reference entry."""
    if entry.source == "literal":
        records = _read_literal(entry)
        reporter.info(
            f"  {entry.key} ({entry.role}): {len(records)} sequence(s) from "
            f"{entry.query} in this repository"
        )
        return records, len(records)
    if entry.source != "uniprot":  # pragma: no cover — validated in panels.py
        raise ReferenceFetchError(f"{entry.key}: unsupported source {entry.source!r}")

    # Over-fetch so that the length window is computed on a representative
    # sample rather than on a truncated head of the result list.
    fetch_limit = max(entry.max_sequences * 2, 1000)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{entry.id.lower()}_{_entry_cache_key(entry, fetch_limit)}.tsv.gz"
    meta_path = cache_path.with_suffix(".meta.json")

    if cache_path.is_file() and meta_path.is_file() and not refresh:
        with gzip.open(cache_path, "rt", encoding="utf-8") as fh:
            text = fh.read()
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        records = _parse_tsv(text)
        reporter.record(f"    {entry.key}: {len(records)} records from cache {cache_path.name}")
        return records, int(meta.get("total_available", -1))

    reporter.info(f"  fetching {entry.key} ({entry.role}) from UniProt")
    text, total = _fetch_uniprot_tsv(entry.query, fetch_limit, reporter)
    records = _parse_tsv(text)
    if not records:
        raise ReferenceFetchError(
            f"{entry.key}: UniProt query returned zero sequences.\n"
            f"  query: {entry.query}\n"
            "A zero-sequence reference fetch is a hard error: it would silently turn into a "
            "count of zero. Fix the query in the panel YAML, or drop the entry."
        )
    with gzip.open(cache_path, "wt", encoding="utf-8") as fh:
        fh.write(text)
    meta_path.write_text(
        json.dumps({"query": entry.query, "total_available": total, "fetched": len(records)}, indent=2),
        encoding="utf-8",
    )
    return records, total


# --------------------------------------------------------------------------- #
# filtering
# --------------------------------------------------------------------------- #


def _length_window(entry: RefEntry, lengths: Sequence[int]) -> tuple[int, int, str]:
    if entry.length_range is not None:
        lo, hi = entry.length_range
        return lo, hi, f"explicit {lo}-{hi} aa"
    median = statistics.median(lengths)
    lo = int(median * (1.0 - entry.length_tolerance))
    hi = int(median * (1.0 + entry.length_tolerance))
    return lo, hi, f"median {median:.0f} aa +/-{entry.length_tolerance:.0%} -> {lo}-{hi} aa"


# --------------------------------------------------------------------------- #
# database build
# --------------------------------------------------------------------------- #


def database_fingerprint(panel_set: PanelSet, entries: Sequence[RefEntry], diamond_version: str) -> str:
    # Deliberately NOT salted with ``__version__``. The fingerprint exists to
    # invalidate the database and the cached searches when something that
    # changes their content changes: the index schema (``INDEX_VERSION``), the
    # DIAMOND build, the reference entries and the residue checks. Releasing a
    # new report version does none of that, and salting with the version threw
    # away an eight-minute translated search per sample on every bump. Bump
    # ``INDEX_VERSION`` when the build itself changes.
    parts = ["openbiota", f"index{INDEX_VERSION}", diamond_version.strip()]
    parts.extend(sorted(e.fingerprint() for e in entries))
    for panel in panel_set.panels:
        rc = panel.residue_check
        if rc and rc.enabled:
            parts.append(
                f"rc:{panel.name}:{rc.anchor_accession}:{rc.canonical_position}:"
                f"{''.join(sorted(rc.accepted_residues))}:{','.join(rc.applies_to)}"
            )
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()[:16]


def diamond_version(diamond: str = "diamond") -> str:
    exe = shutil.which(diamond)
    if exe is None:
        raise DependencyError(
            f"{diamond!r} not found on PATH.\n"
            "Install it natively — on macOS: `brew install diamond`. Bioconda's osx-64 "
            "builds run under Rosetta on Apple Silicon and cost roughly 2x.\n"
            "Run `openbiota doctor` for a full dependency check."
        )
    try:
        proc = subprocess.run(  # noqa: S603 — fixed argv, no shell
            [exe, "version"], capture_output=True, text=True, check=True, timeout=60
        )
    except subprocess.CalledProcessError as exc:
        raise DependencyError(f"`{exe} version` failed: {exc.stderr.strip() or exc}") from exc
    except subprocess.TimeoutExpired as exc:
        raise DependencyError(f"`{exe} version` timed out after 60s") from exc
    return proc.stdout.strip() or "unknown"


def build_reference_database(
    panel_set: PanelSet,
    entries: Sequence[RefEntry],
    *,
    refs_dir: Path,
    reporter: Reporter,
    refresh: bool = False,
    rebuild: bool = False,
    diamond: str = "diamond",
    threads: int = 1,
) -> ReferenceDatabase:
    """Fetch references, filter them, and build (or reuse) the DIAMOND database."""
    version = diamond_version(diamond)
    fingerprint = database_fingerprint(panel_set, entries, version)
    db_dir = refs_dir / "db" / fingerprint
    index_path = db_dir / "index.json"
    dmnd_path = db_dir / "combined.dmnd"
    fasta_path = db_dir / "combined.faa"

    if not rebuild and index_path.is_file() and dmnd_path.is_file():
        try:
            database = _load_index(index_path, dmnd_path, fasta_path)
        except (OSError, ValueError, KeyError) as exc:
            reporter.warn(f"cached reference index unusable ({exc}); rebuilding")
        else:
            reporter.ok(
                f"reference database reused: {database.n_sequences:,} proteins, "
                f"fingerprint {fingerprint} ({human_bytes(dmnd_path.stat().st_size)})"
            )
            return database

    db_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = refs_dir / "cache"

    sequences: list[RefSequence] = []
    stats: dict[str, EntryStats] = {}
    claimed: dict[str, str] = {}
    fasta_lines: list[str] = []
    # Sequence text keyed by database index, needed for anchor mapping below.
    seq_text: dict[int, str] = {}

    # Role precedence matters: a protein claimed by the normaliser or by a
    # target must never also sit in a decoy set, or best-hit competition
    # between the two becomes a coin flip on identical sequences.
    ordered = sorted(entries, key=lambda e: {"normalizer": 0, "target": 1, "decoy": 2}[e.role])

    for entry in ordered:
        records, total_available = fetch_entry_records(entry, cache_dir, reporter, refresh=refresh)
        n_fetched = len(records)

        lengths = [r.length for r in records]
        lo, hi, window = _length_window(entry, lengths)

        kept: list[FetchedRecord] = []
        dropped_length = dropped_identical = dropped_claimed = 0
        seen_sequences: set[bytes] = set()

        for record in records:
            if len(kept) >= entry.max_sequences:
                break
            sequence = clean_sequence(record.sequence)
            if len(sequence) < MIN_REF_LENGTH_AA or not (lo <= len(sequence) <= hi):
                dropped_length += 1
                continue
            if record.accession in claimed:
                dropped_claimed += 1
                continue
            digest = hashlib.blake2b(sequence.encode("ascii"), digest_size=16).digest()
            if digest in seen_sequences:
                dropped_identical += 1
                continue
            seen_sequences.add(digest)
            claimed[record.accession] = entry.key
            kept.append(record)

        if not kept:
            raise ReferenceFetchError(
                f"{entry.key}: every one of {n_fetched} fetched sequences was filtered out.\n"
                f"  length window: {window}\n"
                f"  dropped: {dropped_length} on length, {dropped_identical} identical, "
                f"{dropped_claimed} already claimed by another entry\n"
                "Widen 'length_tolerance' / 'length_range', or fix the query, in the panel YAML."
            )

        kept_lengths: list[int] = []
        n_out_of_scope = 0
        for record in kept:
            sequence = clean_sequence(record.sequence)
            taxa = _parse_lineage(record.lineage)
            phylum = taxa.get("phylum", "")
            # A scoped target keeps its out-of-scope references in the database
            # but demotes them to background, so a read from an organism that
            # carries the fold without the pathway best-hits its own protein and
            # is counted as a decoy rather than as capacity.
            in_scope = entry.scope is None or entry.scope.admits(
                phylum=phylum, protein_name=record.protein_name
            )
            if not in_scope:
                n_out_of_scope += 1
            idx = len(sequences)
            ref = RefSequence(
                idx=idx,
                entry_key=entry.key if in_scope else f"{entry.key}{OUT_OF_SCOPE_SUFFIX}",
                entry_id=entry.id if in_scope else f"{entry.id}{OUT_OF_SCOPE_SUFFIX}",
                role=entry.role if in_scope else "decoy",
                panel=entry.panel,
                accession=record.accession,
                organism=record.organism or "unknown organism",
                length=len(sequence),
                phylum=phylum,
                genus=taxa.get("genus", "") or _genus_from_organism(record.organism),
                family=taxa.get("family", ""),
            )
            sequences.append(ref)
            seq_text[idx] = sequence
            fasta_lines.append(f">{ref.sseqid}\n{sequence}")
            kept_lengths.append(len(sequence))

        if n_out_of_scope:
            scope = entry.scope
            assert scope is not None
            stats[f"{entry.key}{OUT_OF_SCOPE_SUFFIX}"] = EntryStats(
                key=f"{entry.key}{OUT_OF_SCOPE_SUFFIX}",
                panel=entry.panel,
                entry_id=f"{entry.id}{OUT_OF_SCOPE_SUFFIX}",
                role="decoy",
                label=f"{entry.label} — out of scope (counted as background)",
                gene=entry.gene,
                source=entry.source,
                query=entry.query,
                min_identity=entry.min_identity,
                total_available=total_available,
                n_fetched=n_fetched,
                n_dropped_length=0,
                n_dropped_identical=0,
                n_dropped_claimed=0,
                n_sequences=n_out_of_scope,
                mean_length_aa=statistics.fmean(kept_lengths),
                median_length_aa=float(statistics.median(kept_lengths)),
                min_length_aa=min(kept_lengths),
                max_length_aa=max(kept_lengths),
                length_window=window,
                scope_reason=scope.reason,
            )
            reporter.record(
                f"    {entry.key}: {n_out_of_scope}/{len(kept)} references are out of scope "
                f"and count as background — {scope.reason}"
            )

        stats[entry.key] = EntryStats(
            key=entry.key,
            panel=entry.panel,
            entry_id=entry.id,
            role=entry.role,
            label=entry.label,
            gene=entry.gene,
            source=entry.source,
            query=entry.query,
            min_identity=entry.min_identity,
            total_available=total_available,
            n_fetched=n_fetched,
            n_dropped_length=dropped_length,
            n_dropped_identical=dropped_identical,
            n_dropped_claimed=dropped_claimed,
            n_sequences=len(kept),
            mean_length_aa=statistics.fmean(kept_lengths),
            median_length_aa=float(statistics.median(kept_lengths)),
            min_length_aa=min(kept_lengths),
            max_length_aa=max(kept_lengths),
            length_window=window,
            truncated=total_available > len(kept) and len(kept) >= entry.max_sequences,
            scope_reason="" if entry.scope is None else entry.scope.reason,
            n_out_of_scope=n_out_of_scope,
        )
        reporter.record(
            f"    {entry.key}: kept {len(kept)}/{n_fetched} "
            f"(mean {statistics.fmean(kept_lengths):.0f} aa; {window})"
        )

    # ---- active-site residue anchors --------------------------------------- #
    anchors: dict[str, AnchorMapping] = {}
    for panel in panel_set.panels:
        rc = panel.residue_check
        if rc is None or not rc.enabled:
            continue
        applicable = [
            (s.idx, s.accession, seq_text[s.idx])
            for s in sequences
            if s.panel == panel.name and s.entry_id in rc.applies_to
        ]
        if not applicable:
            continue
        mapping = map_anchor_positions(
            panel_name=panel.name,
            residue_check=rc,
            refs=applicable,
            work_dir=db_dir / "anchors",
            reporter=reporter,
            cache_dir=cache_dir,
            diamond=diamond,
            threads=threads,
            refresh=refresh,
        )
        anchors[panel.name] = mapping
        for idx, position in mapping.positions.items():
            current = sequences[idx]
            sequences[idx] = RefSequence(
                **{
                    **asdict(current),
                    "anchor_pos": position,
                    "anchor_residue_ok": mapping.residue_ok.get(idx),
                }
            )
        for entry_id in rc.applies_to:
            key = f"{panel.name}:{entry_id}"
            if key not in stats:
                continue
            mapped = [s for s in sequences if s.entry_key == key and s.anchor_pos is not None]
            consistent = [s for s in mapped if s.anchor_residue_ok]
            stats[key] = EntryStats(
                **{
                    **stats[key].to_json(),
                    "n_anchor_mapped": len(mapped),
                    "n_anchor_residue_ok": len(consistent),
                    "mean_length_residue_ok_aa": (
                        statistics.fmean([s.length for s in consistent]) if consistent else 0.0
                    ),
                }
            )
            reporter.record(
                f"    {key}: {len(consistent)}/{len(mapped)} mapped references carry an "
                f"accepted residue at position {rc.canonical_position}"
            )

    # ---- write FASTA + build DIAMOND index --------------------------------- #
    fasta_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    total_aa = sum(s.length for s in sequences)
    reporter.info(
        f"  combined reference set: {len(sequences):,} proteins, {total_aa:,} aa "
        f"({human_bytes(fasta_path.stat().st_size)})"
    )

    _run_makedb(fasta_path, dmnd_path, diamond=diamond, threads=threads, reporter=reporter)

    built_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    database = ReferenceDatabase(
        fingerprint=fingerprint,
        dmnd_path=dmnd_path,
        fasta_path=fasta_path,
        index_path=index_path,
        sequences=sequences,
        entries=stats,
        diamond_version=version,
        built_at=built_at,
        anchors=anchors,
    )
    _write_index(database)
    reporter.ok(
        f"reference database built: {len(sequences):,} proteins across {len(stats)} entries, "
        f"fingerprint {fingerprint}"
    )
    return database


def _genus_from_organism(organism: str) -> str:
    """First whitespace-delimited token of a binomial, minus bracket notation."""
    token = organism.strip().lstrip("[").split(" ")[0].rstrip("]")
    return token if token[:1].isupper() else ""


def _run_makedb(
    fasta_path: Path,
    dmnd_path: Path,
    *,
    diamond: str,
    threads: int,
    reporter: Reporter,
) -> None:
    exe = shutil.which(diamond)
    if exe is None:  # pragma: no cover — checked earlier
        raise DependencyError(f"{diamond!r} not found on PATH")
    cmd = [
        exe, "makedb",
        "--in", str(fasta_path),
        "--db", str(dmnd_path),
        "--threads", str(threads),
        "--quiet",
    ]
    reporter.record(f"    $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1800)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        raise DependencyError(
            "diamond makedb failed.\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {(exc.stderr or '').strip()}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise DependencyError(f"diamond makedb timed out after 1800s: {exc}") from exc


# --------------------------------------------------------------------------- #
# index persistence
# --------------------------------------------------------------------------- #


def _write_index(database: ReferenceDatabase) -> None:
    payload = {
        "index_version": INDEX_VERSION,
        "openbiota_version": __version__,
        "fingerprint": database.fingerprint,
        "diamond_version": database.diamond_version,
        "built_at": database.built_at,
        "entries": {k: v.to_json() for k, v in database.entries.items()},
        "anchors": {k: v.to_json() for k, v in database.anchors.items()},
        # Compact positional rows keep this file small enough to parse in
        # well under a second even with tens of thousands of references.
        "sequence_columns": [
            "idx", "entry_key", "entry_id", "role", "panel", "accession", "organism",
            "length", "phylum", "genus", "family", "anchor_pos", "anchor_residue_ok",
        ],
        "sequences": [
            [
                s.idx, s.entry_key, s.entry_id, s.role, s.panel, s.accession, s.organism,
                s.length, s.phylum, s.genus, s.family, s.anchor_pos, s.anchor_residue_ok,
            ]
            for s in database.sequences
        ],
    }
    database.index_path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")


def _load_index(index_path: Path, dmnd_path: Path, fasta_path: Path) -> ReferenceDatabase:
    payload = json.loads(index_path.read_text(encoding="utf-8"))
    if int(payload.get("index_version", 0)) != INDEX_VERSION:
        raise ValueError(
            f"index version {payload.get('index_version')} != expected {INDEX_VERSION}"
        )
    sequences = [
        RefSequence(
            idx=row[0], entry_key=row[1], entry_id=row[2], role=row[3], panel=row[4],
            accession=row[5], organism=row[6], length=row[7], phylum=row[8],
            genus=row[9], family=row[10], anchor_pos=row[11], anchor_residue_ok=row[12],
        )
        for row in payload["sequences"]
    ]
    for position, sequence in enumerate(sequences):
        if sequence.idx != position:
            raise ValueError(f"index row {position} has idx {sequence.idx}; rows must be dense")
    entries = {k: EntryStats(**v) for k, v in payload["entries"].items()}
    anchors = {k: AnchorMapping.from_json(v) for k, v in payload.get("anchors", {}).items()}
    return ReferenceDatabase(
        fingerprint=payload["fingerprint"],
        dmnd_path=dmnd_path,
        fasta_path=fasta_path,
        index_path=index_path,
        sequences=sequences,
        entries=entries,
        diamond_version=payload.get("diamond_version", "unknown"),
        built_at=payload.get("built_at", "unknown"),
        anchors=anchors,
    )


def entry_length_stats(database: ReferenceDatabase, key: str) -> tuple[int, float]:
    """Return ``(n_sequences, mean_length_aa)`` for one reference entry."""
    stats = database.entries[key]
    return stats.n_sequences, stats.mean_length_aa
