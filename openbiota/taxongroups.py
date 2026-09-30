"""Curated groups of species, counted together.

A taxon group is a named list of species that share a function or an origin —
butyrate producers, oral-origin translocators, sulfate reducers. Two things
use them.

The report shows each group's summed abundance against the reference cohort,
which is the "microbial groups" section: a reader who cannot interpret
*Roseburia intestinalis* can interpret "butyrate producers, notably low".

Profiles use them as ``CARRIER_ABUNDANCE`` features. That is a deliberately
weaker claim than a gene count and is labelled as one: carriage is not
production. Its value is as a second, independent line on the same biology —
genes counted one way, organisms carrying them counted another. Agreement is
a cross-evidence check; disagreement means a broken map or an unusual
community, and either is worth saying out loud.

Membership is curated from the literature, not inferred, and every group
records where its list came from. A species belongs to a group when
published work puts it there, not when a genome annotation suggests it might.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from .errors import PanelError

#: Report groupings, in the order the report presents them.
GROUP_CATEGORIES: Final = (
    ("function", "Functional groups"),
    ("origin", "Organisms from elsewhere in the body"),
    ("beneficial", "Commonly beneficial organisms"),
    ("opportunist", "Opportunistic organisms"),
    ("surface", "Cell-surface characteristics"),
)
VALID_GROUP_CATEGORIES: Final = frozenset(k for k, _ in GROUP_CATEGORIES)

#: Which direction a raised group abundance points, for report wording.
VALID_GROUP_DIRECTIONS: Final = frozenset(
    {"adverse", "favourable", "context-dependent", "unclear"}
)


def _reader_direction(data: Any, where: str) -> str | None:
    """The direction to colour by, when the curated one is "it depends".

    A group whose literature genuinely cuts both ways stays neutral. One
    with a usual case - mucus degraders, sulfate reducers - says so, and
    grey stops standing in for "no opinion" when an opinion exists.
    """
    value = data.get("reader_direction")
    if value is None:
        return None
    value = str(value).strip()
    _require(
        value in {"adverse", "favourable"},
        where,
        f"'reader_direction' must be adverse or favourable, got {value!r}",
    )
    return value


def _require(condition: bool, where: str, message: str) -> None:
    if not condition:
        raise PanelError(f"{where}: {message}")


def species_key(name: str) -> str:
    """Normalise a species name to the catalogue's underscore form."""
    return name.strip().replace(" ", "_")


@dataclass(frozen=True, slots=True)
class TaxonGroup:
    """A curated set of species reported and scored as one quantity."""

    name: str
    label: str
    category: str
    description: str
    members: tuple[str, ...]
    higher_means: str
    summary: str
    citation: str
    #: Species named in the source literature that this pipeline's catalogue
    #: does not resolve. Recorded rather than dropped, so the report can say
    #: how complete the group is.
    unresolved: tuple[str, ...] = ()
    caveats: tuple[str, ...] = ()
    source_path: Path | None = field(default=None, compare=False)
    #: Direction to colour and rank by where the curated direction is "it
    #: depends" but the dependence has a usual case. Falls back to
    #: `higher_means`; `direction_note` says why.
    reader_direction: str | None = None
    direction_note: str = ""

    @property
    def effective_direction(self) -> str:
        """The direction to colour and rank by."""
        return self.reader_direction or self.higher_means

    def to_json(self) -> dict[str, Any]:
        return {
            "group": self.name,
            "label": self.label,
            "category": self.category,
            "description": " ".join(self.description.split()),
            "members": list(self.members),
            "n_members": len(self.members),
            "unresolved": list(self.unresolved),
            "higher_means": self.higher_means,
            "reader_direction": self.reader_direction,
            "direction_note": self.direction_note,
            "summary": " ".join(self.summary.split()),
            "caveats": [" ".join(c.split()) for c in self.caveats],
            "citation": " ".join(self.citation.split()),
        }


@dataclass(frozen=True, slots=True)
class TaxonGroupSet:
    groups: tuple[TaxonGroup, ...]

    def by_name(self, name: str) -> TaxonGroup | None:
        for group in self.groups:
            if group.name == name:
                return group
        return None

    def carrier_map(self) -> dict[str, tuple[str, ...]]:
        """Group name -> member species, for CARRIER_ABUNDANCE features."""
        return {g.name: g.members for g in self.groups}

    def in_category(self, category: str) -> tuple[TaxonGroup, ...]:
        return tuple(g for g in self.groups if g.category == category)

    def to_json(self) -> list[dict[str, Any]]:
        return [g.to_json() for g in self.groups]


def _str_tuple(data: Any, key: str, where: str) -> tuple[str, ...]:
    raw = data.get(key) or []
    _require(isinstance(raw, list), where, f"{key!r} must be a list")
    return tuple(str(x).strip() for x in raw if str(x).strip())


def parse_taxon_group(data: Any, *, source_path: Path | None = None) -> TaxonGroup:
    where = f"taxon group {source_path.name if source_path else '<inline>'}"
    _require(isinstance(data, dict), where, "must be a mapping")

    name = str(data.get("name", "")).strip()
    _require(bool(name), where, "'name' must not be empty")
    where = f"taxon group {name!r}"

    category = str(data.get("category", "")).strip()
    _require(
        category in VALID_GROUP_CATEGORIES,
        where,
        f"'category' must be one of {sorted(VALID_GROUP_CATEGORIES)}, got {category!r}",
    )

    direction = str(data.get("higher_means", "unclear")).strip()
    _require(
        direction in VALID_GROUP_DIRECTIONS,
        where,
        f"'higher_means' must be one of {sorted(VALID_GROUP_DIRECTIONS)}, got {direction!r}",
    )

    members = tuple(species_key(m) for m in _str_tuple(data, "members", where))
    _require(bool(members), where, "'members' must not be empty")
    _require(
        len(set(members)) == len(members),
        where,
        "'members' contains a duplicate species",
    )

    for key in ("label", "description", "summary", "citation"):
        value = data.get(key)
        _require(
            isinstance(value, str) and bool(value.strip()),
            where,
            f"{key!r} must be a non-empty string",
        )

    return TaxonGroup(
        name=name,
        label=str(data["label"]).strip(),
        category=category,
        description=str(data["description"]).strip(),
        members=members,
        higher_means=direction,
        reader_direction=_reader_direction(data, where),
        direction_note=" ".join(str(data.get("direction_note") or "").split()),
        summary=str(data["summary"]).strip(),
        citation=str(data["citation"]).strip(),
        unresolved=tuple(species_key(u) for u in _str_tuple(data, "unresolved", where)),
        caveats=_str_tuple(data, "caveats", where),
        source_path=source_path,
    )


def load_taxon_group_set(groups_dir: Path) -> TaxonGroupSet:
    """Load every ``*.yaml`` in ``groups_dir``, sorted by name."""
    if not groups_dir.is_dir():
        return TaxonGroupSet(groups=())
    groups: list[TaxonGroup] = []
    for path in sorted(groups_dir.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        groups.append(parse_taxon_group(data, source_path=path))
    names = [g.name for g in groups]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise PanelError(f"duplicate taxon group names: {duplicates}")
    return TaxonGroupSet(groups=tuple(sorted(groups, key=lambda g: g.name)))


def resolve_against_catalogue(
    group_set: TaxonGroupSet, catalogue: Sequence[str]
) -> dict[str, tuple[tuple[str, ...], tuple[str, ...]]]:
    """Split each group's members into (present, absent) against a catalogue.

    A member absent from the taxonomic catalogue is not evidence of absence
    from the sample: the profiler cannot report what it cannot name. The
    report shows the resolved fraction so a group built from twelve species
    of which four are nameable does not read as a complete measurement.
    """
    known = set(catalogue)
    out: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {}
    for group in group_set.groups:
        present = tuple(m for m in group.members if m in known)
        absent = tuple(m for m in group.members if m not in known)
        out[group.name] = (present, absent)
    return out
