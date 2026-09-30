"""Declarative gene panels: load and validate the YAML in ``panels/``.

Adding a gene family is a config edit, never a code change. Every field is
validated here with an actionable error message, because a typo in a UniProt
query silently changes what the tool measures.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.errors import PanelError

#: The universal single-copy normaliser is not a per-panel entry; it lives in
#: this file and is included in every run.
NORMALIZER_FILE: Final = "_normalizer.yaml"

ROLE_TARGET: Final = "target"
ROLE_DECOY: Final = "decoy"
ROLE_NORMALIZER: Final = "normalizer"

#: Where an entry's references come from. `literal` reads sequences shipped
#: in this repository, which is the only way to hold a protein that has no
#: UniProt entry: the urolithin dehydroxylase operon was published to
#: GenBank in 2025 and searching for its gene symbols instead returns an
#: *Aspergillus* protein of the same name.
VALID_SOURCES: Final = frozenset({"uniprot", "literal"})
VALID_AGGREGATES: Final = frozenset({"sum", "median", "max", "mean"})
AMINO_ACIDS: Final = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Default cap on sequences fetched per reference entry. Reference sets are
#: small by design; an unbounded fetch would be a silent correctness hazard.
DEFAULT_MAX_SEQUENCES: Final = 4000

#: Reference sequences whose length falls outside
#: ``[1-tol, 1+tol] x median`` are dropped. This removes database fragments and
#: multi-domain fusions, both of which would corrupt the length normalisation.
DEFAULT_LENGTH_TOLERANCE: Final = 0.25

#: Minimum aligned amino acids for a read hit to be considered at all.
#: Reads are 100-151 bp, i.e. 33-50 translated residues.
DEFAULT_MIN_ALIGNMENT_AA: Final = 25


#: Protein-name markers that UniProt uses when an annotation is inferred by
#: rule rather than curated from experiment. A name carrying one of these is
#: an assertion about sequence similarity, not about measured activity, which
#: is exactly the case where a promiscuous enzyme family produces a reference
#: set full of organisms that do not run the pathway.
INFERRED_NAME_MARKERS: Final = ("probable ", "putative ", "predicted ", "uncharacterized")


@dataclass(frozen=True, slots=True)
class Scope:
    """Which reference proteins a target may legitimately claim.

    A DIAMOND target set built from a protein-name query is a similarity set,
    not a function set. For a promiscuous family (butyrate kinase and the
    branched-chain carboxylic acid kinases; GABA transaminase; the CoA
    transferases) that set fills with organisms that carry the fold but not
    the pathway, and a sample dominated by those organisms then scores high
    for a capacity it does not have.

    A scope names the clades where the function is actually characterised.
    Out-of-scope references are not discarded — they are moved into the
    panel's decoy pool, so a read from such an organism still best-hits its
    own protein and is counted as background instead of as capacity. Removing
    them from the database entirely would be worse: the read would fall
    through to a distant in-scope reference and be miscounted.
    """

    #: Phyla (UniProt lineage names) where the function is characterised.
    include_phyla: tuple[str, ...] = ()
    #: Phyla explicitly known to carry the fold without the function.
    exclude_phyla: tuple[str, ...] = ()
    #: Treat rule-inferred annotations as out of scope (see markers above).
    exclude_inferred: bool = False
    #: Substrings in the protein name that mark the wrong enzyme outright.
    exclude_name_contains: tuple[str, ...] = ()
    #: Why, for the report and the audit trail.
    reason: str = ""

    @property
    def active(self) -> bool:
        return bool(
            self.include_phyla or self.exclude_phyla
            or self.exclude_inferred or self.exclude_name_contains
        )

    def admits(self, *, phylum: str, protein_name: str) -> bool:
        """True when this reference may count as target capacity."""
        name = protein_name.lower()
        if self.exclude_inferred and any(m in name for m in INFERRED_NAME_MARKERS):
            return False
        if any(s.lower() in name for s in self.exclude_name_contains):
            return False
        if self.exclude_phyla and phylum in self.exclude_phyla:
            return False
        return not (self.include_phyla and phylum not in self.include_phyla)

    def fingerprint(self) -> str:
        return (
            f"{','.join(self.include_phyla)}|{','.join(self.exclude_phyla)}|"
            f"{int(self.exclude_inferred)}|{','.join(self.exclude_name_contains)}"
        )


@dataclass(frozen=True, slots=True)
class RefEntry:
    """One reference sequence set: a target gene, a decoy family, or rpoB."""

    id: str
    role: str
    source: str
    query: str
    min_identity: float
    label: str
    panel: str
    gene: str | None = None
    max_sequences: int = DEFAULT_MAX_SEQUENCES
    length_range: tuple[int, int] | None = None
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE
    note: str | None = None
    scope: Scope | None = None

    @property
    def key(self) -> str:
        """Globally unique identifier: ``panel:ENTRY``."""
        return f"{self.panel}:{self.id}"

    def fingerprint(self) -> str:
        """Everything that changes which sequences end up in the database."""
        lr = f"{self.length_range[0]}-{self.length_range[1]}" if self.length_range else "median"
        sc = self.scope.fingerprint() if self.scope is not None and self.scope.active else "-"
        return (
            f"{self.key}|{self.role}|{self.source}|{self.query}|{self.max_sequences}|"
            f"{lr}|{self.length_tolerance}|{sc}"
        )


@dataclass(frozen=True, slots=True)
class ResidueCheck:
    """Partial recovery of a full-length active-site residue filter.

    At database build time every target reference is aligned to ``anchor_accession``
    and the reference's own position corresponding to ``canonical_position`` is
    recorded. At search time, reads whose alignment spans that mapped position
    are checked for one of ``accepted_residues``.
    """

    enabled: bool
    anchor_accession: str
    canonical_position: int
    accepted_residues: frozenset[str]
    applies_to: tuple[str, ...]
    description: str = ""


VALID_DIRECTIONS: Final = frozenset({"adverse", "favourable", "context-dependent", "unclear"})
VALID_STRENGTHS: Final = frozenset({"strong", "moderate", "limited", "preliminary"})

#: Report groupings for panels, in the order the report presents them.
PANEL_CATEGORIES: Final = (
    ("metabolite", "Metabolite production capacity"),
    ("neuroactive", "Neuroactive compounds"),
    ("vitamin", "Vitamin synthesis"),
    # Spec 0.8.3 A01: what the community can break down, as distinct from
    # what it can make. A substrate panel answers "can these microbes open
    # this fibre", which is a different question from every other category.
    ("substrate", "Dietary substrate breakdown"),
    ("gas", "Gas production"),
    ("detox", "Detoxification and breakdown"),
    ("virulence", "Toxins and virulence"),
)
VALID_CATEGORIES: Final = frozenset(k for k, _ in PANEL_CATEGORIES)


@dataclass(frozen=True, slots=True)
class Interpretation:
    """What the literature reports about this metabolite.

    Attached to the panel so the consumer report can say something meaningful
    about a number instead of just printing it. ``higher_means`` is the
    direction of the *reported association for the metabolite* — never a
    statement about the person being screened.
    """

    what_it_is: str
    made_from: str
    higher_means: str
    evidence_strength: str
    summary: str
    #: The specific findings with study context. Named ``cognitive_evidence``
    #: when this report only covered Alzheimer's disease; panels may now
    #: write it as ``evidence_detail``, which is the accurate name for a
    #: metabolite whose literature is cardiovascular or gastrointestinal.
    evidence_detail: str
    citation: str
    #: What a reading in each direction means for this metabolite. Keyed
    #: "higher" / "typical" / "lower" so a status is never shown without an
    #: explanation of what that status implies.
    implications: dict[str, str] = field(default_factory=dict)
    #: The direction a reader should treat the reading as having, where the
    #: literature's own answer is "it depends" but the dependence has a
    #: usual case. Falls back to `higher_means`. `direction_note` says why.
    reader_direction: str | None = None
    direction_note: str = ""
    #: Organism-name prefix -> direction, for pathways whose consequence
    #: depends on the producer. Resolved against the organisms actually
    #: behind the reading in each sample.
    driver_exceptions: dict[str, str] = field(default_factory=dict)

    @property
    def effective_direction(self) -> str:
        """The direction to colour and rank by, absent per-sample evidence."""
        return self.reader_direction or self.higher_means

    def implication_for(self, percentile: float | None) -> str:
        """Pick the explanation matching where the reading fell."""
        if percentile is None:
            return ""
        if percentile >= 75:
            key = "higher"
        elif percentile <= 25:
            key = "lower"
        else:
            key = "typical"
        return " ".join(self.implications.get(key, "").split())

    def to_json(self) -> dict[str, Any]:
        return {
            "what_it_is": " ".join(self.what_it_is.split()),
            "made_from": " ".join(self.made_from.split()),
            "higher_means": self.higher_means,
            "evidence_strength": self.evidence_strength,
            "summary": " ".join(self.summary.split()),
            "evidence_detail": " ".join(self.evidence_detail.split()),
            "citation": " ".join(self.citation.split()),
            "implications": {
                k: " ".join(v.split()) for k, v in sorted(self.implications.items())
            },
        }


@dataclass(frozen=True, slots=True)
class Panel:
    """A gene family screened as one unit."""

    name: str
    description: str
    metabolite: str
    targets: tuple[RefEntry, ...]
    decoys: tuple[RefEntry, ...]
    citation: str
    pathway: str | None = None
    organisms: tuple[str, ...] = ()
    residue_check: ResidueCheck | None = None
    min_fragments_for_stability: int = 20
    min_alignment_aa: int = DEFAULT_MIN_ALIGNMENT_AA
    #: How many copies of one headline gene a carrier genome typically holds.
    #: Most marker genes are single-copy (1); a paralog-rich family such as the
    #: beta-glucuronidases, where one Bacteroides genome can carry several, is
    #: declared higher so the plausibility ceiling is not tripped by biology.
    expected_copies_per_genome: float = 1.0
    aggregate: str = "sum"
    aggregate_reason: str = ""
    #: Target ids the panel-level aggregate is computed from. Every target is
    #: always reported; this only narrows the headline figure to the genes
    #: whose reference sets are specific enough to carry it. Empty means all.
    aggregate_from: tuple[str, ...] = ()
    upper_bound_reason: str | None = None
    caveats: tuple[str, ...] = ()
    extension: bool = False
    interpretation: Interpretation | None = None
    #: Report grouping. The original nine panels are all metabolite pathways;
    #: later panels cover neuroactive compounds, vitamins, gases, detoxification
    #: and virulence, and the report presents each group under its own heading.
    category: str = "metabolite"
    source_path: Path | None = field(default=None, compare=False)

    def aggregate_target_ids(self) -> tuple[str, ...]:
        return self.aggregate_from or tuple(t.id for t in self.targets)

    def plausible_ceiling(self, per_gene_ceiling: float) -> float:
        """Copies per 100 genomes above which the figure is implausible.

        ``per_gene_ceiling`` is the ceiling for one single-copy gene. A ``sum``
        panel adds one such figure per headline gene, and a multi-copy family
        scales by its declared copies per genome.
        """
        n_headline = len(self.aggregate_target_ids()) if self.aggregate == "sum" else 1
        return per_gene_ceiling * max(1, n_headline) * self.expected_copies_per_genome

    @property
    def entries(self) -> tuple[RefEntry, ...]:
        return self.targets + self.decoys

    def target_ids(self) -> tuple[str, ...]:
        return tuple(t.id for t in self.targets)


@dataclass(frozen=True, slots=True)
class PanelSet:
    """All panels plus the built-in rpoB normaliser."""

    panels: tuple[Panel, ...]
    normalizer: RefEntry
    normalizer_citation: str
    normalizer_description: str

    def by_name(self, name: str) -> Panel:
        for p in self.panels:
            if p.name == name:
                return p
        raise PanelError(
            f"unknown panel {name!r}; available: {', '.join(p.name for p in self.panels)}"
        )

    def select(self, names: Sequence[str] | None) -> tuple[Panel, ...]:
        if not names:
            return self.panels
        return tuple(self.by_name(n) for n in names)

    def restrict(self, names: Sequence[str] | None) -> PanelSet:
        """The same set narrowed to ``names`` (all panels when ``names`` is empty)."""
        if not names:
            return self
        return PanelSet(
            panels=self.select(names),
            normalizer=self.normalizer,
            normalizer_citation=self.normalizer_citation,
            normalizer_description=self.normalizer_description,
        )

    def all_entries(self) -> tuple[RefEntry, ...]:
        out: list[RefEntry] = [self.normalizer]
        for p in self.panels:
            out.extend(p.entries)
        return tuple(out)


# --------------------------------------------------------------------------- #
# validation helpers
# --------------------------------------------------------------------------- #


def _require(cond: bool, where: str, message: str) -> None:
    if not cond:
        raise PanelError(f"{where}: {message}")


def _get(data: Mapping[str, Any], key: str, where: str, *, expected: type | tuple[type, ...]) -> Any:
    if key not in data:
        raise PanelError(f"{where}: missing required key {key!r}")
    value = data[key]
    if not isinstance(value, expected):
        names = expected.__name__ if isinstance(expected, type) else "/".join(t.__name__ for t in expected)
        raise PanelError(f"{where}: key {key!r} must be {names}, got {type(value).__name__}")
    return value


def _opt_str_tuple(data: Mapping[str, Any], key: str, where: str) -> tuple[str, ...]:
    raw = data.get(key) or []
    if isinstance(raw, str):
        raw = [raw]
    if not isinstance(raw, list) or not all(isinstance(x, str) for x in raw):
        raise PanelError(f"{where}: key {key!r} must be a string or a list of strings")
    return tuple(raw)


def _parse_length_range(raw: Any, where: str) -> tuple[int, int] | None:
    if raw is None:
        return None
    if (
        not isinstance(raw, list)
        or len(raw) != 2
        or not all(isinstance(x, int) for x in raw)
    ):
        raise PanelError(f"{where}: 'length_range' must be a two-item list of integers [min, max]")
    lo, hi = int(raw[0]), int(raw[1])
    _require(0 < lo < hi, where, f"'length_range' must satisfy 0 < min < max, got [{lo}, {hi}]")
    return lo, hi


def _parse_entry(raw: Any, *, role: str, panel: str, index: int) -> RefEntry:
    where = f"panel {panel!r} {role}[{index}]"
    if not isinstance(raw, dict):
        raise PanelError(f"{where}: each entry must be a mapping, got {type(raw).__name__}")

    entry_id = _get(raw, "id", where, expected=str).strip()
    _require(bool(entry_id), where, "'id' must not be empty")
    _require(
        entry_id.replace("-", "").replace(".", "").isalnum(),
        where,
        f"'id' must be alphanumeric (plus '-' and '.'), got {entry_id!r}; "
        "'~' and whitespace are reserved by the FASTA header encoding",
    )

    source = _get(raw, "source", where, expected=str).strip()
    _require(
        source in VALID_SOURCES,
        where,
        f"unsupported source {source!r}; supported: {', '.join(sorted(VALID_SOURCES))}",
    )

    query = _get(raw, "query", where, expected=str).strip()
    _require(bool(query), where, "'query' must not be empty")

    min_identity = raw.get("min_identity", 50.0)
    _require(
        isinstance(min_identity, (int, float)) and 0.0 <= float(min_identity) <= 100.0,
        where,
        f"'min_identity' must be a number in [0, 100], got {min_identity!r}",
    )

    max_sequences = raw.get("max_sequences", DEFAULT_MAX_SEQUENCES)
    _require(
        isinstance(max_sequences, int) and max_sequences > 0,
        where,
        f"'max_sequences' must be a positive integer, got {max_sequences!r}",
    )

    tolerance = raw.get("length_tolerance", DEFAULT_LENGTH_TOLERANCE)
    _require(
        isinstance(tolerance, (int, float)) and 0.0 < float(tolerance) <= 1.0,
        where,
        f"'length_tolerance' must be in (0, 1], got {tolerance!r}",
    )

    return RefEntry(
        id=entry_id,
        role=role,
        source=source,
        query=query,
        min_identity=float(min_identity),
        label=str(raw.get("label") or entry_id),
        panel=panel,
        gene=str(raw["gene"]) if raw.get("gene") else None,
        max_sequences=int(max_sequences),
        length_range=_parse_length_range(raw.get("length_range"), where),
        length_tolerance=float(tolerance),
        note=str(raw["note"]) if raw.get("note") else None,
        scope=_parse_scope(raw.get("scope"), where=where, role=role),
    )


def _parse_scope(raw: Any, *, where: str, role: str) -> Scope | None:
    if raw is None:
        return None
    _require(isinstance(raw, dict), where, f"'scope' must be a mapping, got {type(raw).__name__}")
    _require(
        role == ROLE_TARGET, where,
        "'scope' only applies to targets: it decides what counts as capacity, and a decoy "
        "or the normaliser has no capacity to count",
    )
    unknown = set(raw) - {
        "include_phyla", "exclude_phyla", "exclude_inferred", "exclude_name_contains", "reason",
    }
    _require(not unknown, where, f"unknown 'scope' keys: {', '.join(sorted(unknown))}")

    def strs(key: str) -> tuple[str, ...]:
        if key not in raw:
            return ()
        value = raw[key]
        _require(
            isinstance(value, list) and bool(value)
            and all(isinstance(v, str) and v.strip() for v in value),
            where, f"'scope.{key}' must be a non-empty list of non-empty strings",
        )
        return tuple(v.strip() for v in value)

    inferred = raw.get("exclude_inferred", False)
    _require(isinstance(inferred, bool), where, "'scope.exclude_inferred' must be a boolean")
    reason = str(raw.get("reason") or "").strip()
    scope = Scope(
        include_phyla=strs("include_phyla"),
        exclude_phyla=strs("exclude_phyla"),
        exclude_inferred=bool(inferred),
        exclude_name_contains=strs("exclude_name_contains"),
        reason=" ".join(reason.split()),
    )
    _require(scope.active, where, "'scope' is present but empty; remove it or give it a rule")
    _require(
        bool(scope.reason), where,
        "'scope' must carry a 'reason' citing why these references cannot count as capacity: "
        "narrowing a target changes every percentile downstream and must be justifiable",
    )
    _require(
        not (scope.include_phyla and scope.exclude_phyla), where,
        "'scope' takes include_phyla or exclude_phyla, not both",
    )
    return scope


def _parse_residue_check(raw: Any, *, panel: str, target_ids: Sequence[str]) -> ResidueCheck | None:
    if raw is None:
        return None
    where = f"panel {panel!r} residue_check"
    if not isinstance(raw, dict):
        raise PanelError(f"{where}: must be a mapping, got {type(raw).__name__}")

    enabled = raw.get("enabled", True)
    _require(isinstance(enabled, bool), where, "'enabled' must be a boolean")

    anchor = _get(raw, "anchor_accession", where, expected=str).strip()
    _require(bool(anchor), where, "'anchor_accession' must not be empty")

    position = _get(raw, "canonical_position", where, expected=int)
    _require(position > 0, where, f"'canonical_position' must be >= 1, got {position}")

    residues_raw = _get(raw, "accepted_residues", where, expected=list)
    residues = frozenset(str(r).strip().upper() for r in residues_raw)
    _require(bool(residues), where, "'accepted_residues' must not be empty")
    bad = sorted(r for r in residues if r not in AMINO_ACIDS)
    _require(not bad, where, f"'accepted_residues' contains non-amino-acid codes: {bad}")

    applies_to = _opt_str_tuple(raw, "applies_to", where) or tuple(target_ids)
    unknown = sorted(set(applies_to) - set(target_ids))
    _require(
        not unknown,
        where,
        f"'applies_to' names targets that do not exist in this panel: {unknown}; "
        f"available targets: {list(target_ids)}",
    )

    return ResidueCheck(
        enabled=bool(enabled),
        anchor_accession=anchor,
        canonical_position=int(position),
        accepted_residues=residues,
        applies_to=applies_to,
        description=str(raw.get("description") or ""),
    )


def _parse_interpretation(raw: Any, *, panel: str) -> Interpretation | None:
    if raw is None:
        return None
    where = f"panel {panel!r} interpretation"
    if not isinstance(raw, dict):
        raise PanelError(f"{where}: must be a mapping, got {type(raw).__name__}")

    direction = str(raw.get("higher_means", "unclear")).strip()
    _require(
        direction in VALID_DIRECTIONS,
        where,
        f"'higher_means' must be one of {sorted(VALID_DIRECTIONS)}, got {direction!r}",
    )
    strength = str(raw.get("evidence_strength", "limited")).strip()
    _require(
        strength in VALID_STRENGTHS,
        where,
        f"'evidence_strength' must be one of {sorted(VALID_STRENGTHS)}, got {strength!r}",
    )
    # `evidence_detail` is the current name; `cognitive_evidence` is what the
    # Alzheimer's-only version of this report called it, and the nine original
    # panels still use it.
    detail = raw.get("evidence_detail") or raw.get("cognitive_evidence")
    _require(
        isinstance(detail, str) and bool(detail.strip()),
        where,
        "'evidence_detail' must be a non-empty string",
    )
    for key in ("what_it_is", "made_from", "summary", "citation"):
        value = raw.get(key)
        _require(
            isinstance(value, str) and bool(value.strip()),
            where,
            f"{key!r} must be a non-empty string",
        )
    implications_raw = raw.get("implications") or {}
    if not isinstance(implications_raw, dict):
        raise PanelError(f"{where}: 'implications' must be a mapping")
    missing = sorted({"higher", "typical", "lower"} - set(implications_raw))
    _require(
        not missing,
        where,
        f"'implications' must define {sorted(('higher', 'typical', 'lower'))}; missing {missing}",
    )
    for key, value in implications_raw.items():
        _require(
            isinstance(value, str) and bool(value.strip()),
            where,
            f"implications[{key!r}] must be a non-empty string",
        )

    reader_direction = raw.get("reader_direction")
    if reader_direction is not None:
        reader_direction = str(reader_direction).strip()
        _require(
            reader_direction in {"adverse", "favourable"},
            where,
            f"reader_direction must be adverse or favourable, got {reader_direction!r}",
        )
    exceptions_raw = raw.get("driver_exceptions") or {}
    _require(isinstance(exceptions_raw, dict), where, "driver_exceptions must be a mapping")
    for k, v in exceptions_raw.items():
        _require(str(v) in {"adverse", "favourable"}, where,
                 f"driver_exceptions[{k!r}] must be adverse or favourable")
    return Interpretation(
        what_it_is=str(raw["what_it_is"]).strip(),
        made_from=str(raw["made_from"]).strip(),
        higher_means=direction,
        evidence_strength=strength,
        summary=str(raw["summary"]).strip(),
        evidence_detail=detail.strip(),
        citation=str(raw["citation"]).strip(),
        implications={k: str(v).strip() for k, v in implications_raw.items()},
        reader_direction=reader_direction,
        direction_note=" ".join(str(raw.get("direction_note") or "").split()),
        driver_exceptions={str(k): str(v) for k, v in exceptions_raw.items()},
    )


def parse_panel(data: Any, *, source_path: Path | None = None) -> Panel:
    """Validate one panel mapping. Raises :class:`PanelError` with context."""
    origin = str(source_path) if source_path else "<inline>"
    if not isinstance(data, dict):
        raise PanelError(f"{origin}: panel file must contain a YAML mapping, got {type(data).__name__}")

    name = _get(data, "name", origin, expected=str).strip()
    _require(bool(name), origin, "'name' must not be empty")
    _require(
        name.replace("_", "").replace("-", "").isalnum(),
        origin,
        f"'name' must be alphanumeric (plus '_' and '-'), got {name!r}",
    )

    where = f"panel {name!r}"
    description = _get(data, "description", where, expected=str).strip()
    _require(bool(description), where, "'description' must not be empty")

    metabolite = _get(data, "metabolite", where, expected=str).strip()
    _require(bool(metabolite), where, "'metabolite' must not be empty")

    raw_targets = _get(data, "targets", where, expected=list)
    _require(bool(raw_targets), where, "'targets' must contain at least one entry")
    targets = tuple(
        _parse_entry(t, role=ROLE_TARGET, panel=name, index=i) for i, t in enumerate(raw_targets)
    )

    raw_decoys = data.get("decoys") or []
    _require(isinstance(raw_decoys, list), where, "'decoys' must be a list when present")
    decoys = tuple(
        _parse_entry(d, role=ROLE_DECOY, panel=name, index=i) for i, d in enumerate(raw_decoys)
    )

    seen: set[str] = set()
    for entry in targets + decoys:
        _require(entry.id not in seen, where, f"duplicate entry id {entry.id!r}")
        seen.add(entry.id)

    citation = _get(data, "citation", where, expected=str).strip()
    _require(bool(citation), where, "'citation' must not be empty")

    min_fragments = data.get("min_fragments_for_stability", 20)
    _require(
        isinstance(min_fragments, int) and min_fragments >= 0,
        where,
        f"'min_fragments_for_stability' must be a non-negative integer, got {min_fragments!r}",
    )

    expected_copies = data.get("expected_copies_per_genome", 1.0)
    _require(
        isinstance(expected_copies, (int, float)) and not isinstance(expected_copies, bool)
        and 1.0 <= float(expected_copies) <= 10.0,
        where,
        f"'expected_copies_per_genome' must be a number in [1, 10], got {expected_copies!r}",
    )

    min_alignment = data.get("min_alignment_aa", DEFAULT_MIN_ALIGNMENT_AA)
    _require(
        isinstance(min_alignment, int) and min_alignment > 0,
        where,
        f"'min_alignment_aa' must be a positive integer, got {min_alignment!r}",
    )

    aggregate = str(data.get("aggregate", "sum")).strip()
    _require(
        aggregate in VALID_AGGREGATES,
        where,
        f"'aggregate' must be one of {sorted(VALID_AGGREGATES)}, got {aggregate!r}",
    )

    extension = data.get("extension", False)
    _require(isinstance(extension, bool), where, "'extension' must be a boolean")

    category = str(data.get("category") or "metabolite").strip()
    _require(
        category in VALID_CATEGORIES,
        where,
        f"'category' must be one of {sorted(VALID_CATEGORIES)}, got {category!r}",
    )

    aggregate_from = _opt_str_tuple(data, "aggregate_from", where)
    target_ids = [t.id for t in targets]
    unknown_aggregate = sorted(set(aggregate_from) - set(target_ids))
    _require(
        not unknown_aggregate,
        where,
        f"'aggregate_from' names targets that do not exist in this panel: {unknown_aggregate}; "
        f"available targets: {target_ids}",
    )

    return Panel(
        name=name,
        description=description,
        metabolite=metabolite,
        pathway=str(data["pathway"]).strip() if data.get("pathway") else None,
        organisms=_opt_str_tuple(data, "organisms", where),
        targets=targets,
        decoys=decoys,
        residue_check=_parse_residue_check(
            data.get("residue_check"), panel=name, target_ids=[t.id for t in targets]
        ),
        min_fragments_for_stability=int(min_fragments),
        min_alignment_aa=int(min_alignment),
        expected_copies_per_genome=float(expected_copies),
        aggregate=aggregate,
        aggregate_reason=str(data.get("aggregate_reason") or ""),
        aggregate_from=aggregate_from,
        upper_bound_reason=str(data["upper_bound_reason"]) if data.get("upper_bound_reason") else None,
        caveats=_opt_str_tuple(data, "caveats", where),
        citation=citation,
        extension=bool(extension),
        interpretation=_parse_interpretation(data.get("interpretation"), panel=name),
        category=category,
        source_path=source_path,
    )


def _load_yaml(path: Path) -> Any:
    if not path.is_file():
        raise PanelError(f"panel file not found: {path}")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise PanelError(f"cannot read panel file {path}: {exc}") from exc
    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise PanelError(f"{path}: invalid YAML — {exc}") from exc


def parse_normalizer(data: Any, *, source_path: Path | None = None) -> tuple[RefEntry, str, str]:
    origin = str(source_path) if source_path else "<inline>"
    if not isinstance(data, dict):
        raise PanelError(f"{origin}: normaliser file must contain a YAML mapping")
    raw = _get(data, "normalizer", origin, expected=dict)
    entry = _parse_entry(raw, role=ROLE_NORMALIZER, panel="_normalizer", index=0)
    citation = str(data.get("citation") or "").strip()
    _require(bool(citation), origin, "'citation' must not be empty")
    description = str(data.get("description") or "").strip()
    _require(bool(description), origin, "'description' must not be empty")
    return entry, citation, description


def load_panel_set(panels_dir: Path) -> PanelSet:
    """Load ``panels/_normalizer.yaml`` plus every other ``*.yaml`` in the directory.

    Served from a fingerprinted cache when none of the files has changed.
    """
    from openbiota.fastcache import cached_load

    if not panels_dir.is_dir():
        raise PanelError(f"panels directory not found: {panels_dir}")
    return cached_load(
        "panels", sorted(panels_dir.glob("*.yaml")), lambda: _load_panel_set_uncached(panels_dir)
    )


def _load_panel_set_uncached(panels_dir: Path) -> PanelSet:

    normalizer_path = panels_dir / NORMALIZER_FILE
    entry, citation, description = parse_normalizer(
        _load_yaml(normalizer_path), source_path=normalizer_path
    )

    files = sorted(
        p for p in panels_dir.glob("*.yaml") if p.name != NORMALIZER_FILE and not p.name.startswith(".")
    )
    if not files:
        raise PanelError(f"no panel definitions found in {panels_dir} (expected *.yaml)")

    panels: list[Panel] = []
    names: dict[str, Path] = {}
    for path in files:
        panel = parse_panel(_load_yaml(path), source_path=path)
        if panel.name in names:
            raise PanelError(
                f"duplicate panel name {panel.name!r} in {path} (already defined in {names[panel.name]})"
            )
        names[panel.name] = path
        panels.append(panel)

    _check_global_entry_ids(panels, entry)
    return PanelSet(
        panels=tuple(panels),
        normalizer=entry,
        normalizer_citation=citation,
        normalizer_description=description,
    )


def _check_global_entry_ids(panels: Iterable[Panel], normalizer: RefEntry) -> None:
    """Entry keys must be globally unique — they index the shared DIAMOND database."""
    seen: dict[str, str] = {normalizer.key: "built-in normaliser"}
    for panel in panels:
        for e in panel.entries:
            if e.key in seen:
                raise PanelError(f"duplicate reference entry key {e.key!r} (also in {seen[e.key]})")
            seen[e.key] = f"panel {panel.name!r}"
