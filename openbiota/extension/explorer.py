"""A10 — the complete organism explorer and taxonomy resolution.

BUILD_SPEC_v0.8.3 section 7.1. The existing organism pages keep their values;
this adds the view that lets a reader look an organism up by any name they
have for it and get a straight answer, including "not detected here".

The distinctions this module exists to keep:

* **Resolving a name is not a detection.** Every lookup returns a resolution
  state *and*, separately, a detection state. Finding that *Ruminococcus
  gnavus* is now *Mediterraneibacter gnavus* tells you nothing about whether
  it is in the sample.
* **A rollup discloses its members.** A genus total says which species it
  summed and which it excluded, and a genus total is never placed beside its
  own species as though they were independent slices of one pie.
* **Unknown is an interpretation, not a harm.** The role composition keeps
  unassigned mass visible rather than shrinking the denominator until the
  categories look complete.
* **A species label is not a strain.** Detecting *Lacticaseibacillus
  rhamnosus* does not identify the commercial strain, its toxin carriage or
  a disease-causing subtype.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.extension.schema import ExtensionMetric, fingerprint

METHOD: Final = "ext083.organism_explorer/1.0"
FIXTURE_FILE: Final = Path(__file__).resolve().parents[2] / "extension" / "taxonomy_fixture.yaml"

#: How a searched name resolved against this run's catalogue. Section 7.1
#: requires each of these to be reachable and distinguishable.
RESOLUTION_STATES: Final[frozenset[str]] = frozenset({
    "exact_identity",           # the catalogue holds this exact name
    "accepted_synonym",         # an older or alternative name for one organism
    "cluster_correspondence",   # resolves to a genome cluster, not a binomial
    "split_or_merge_ambiguous", # one name, several current concepts
    "rank_rollup",              # a genus or phylum, answered as a rollup
    "not_in_catalogue",         # a valid name this assay does not carry
    "unresolvable",             # cannot be resolved to any concept
})

#: What the sample says about an organism, independently of how its name
#: resolved. Kept apart so a resolution can never read as a finding.
DETECTION_STATES: Final[frozenset[str]] = frozenset({
    "supported_detection",
    "below_quantification_limit",
    "no_supported_detection",
    "not_assessable",
})

#: GTDB split suffixes: `Ruminococcus_B`, `Blautia_A`. The suffix marks a
#: genus boundary in that taxonomy, not a different organism concept.
_SPLIT = re.compile(r"^(?P<base>[A-Z][a-z]+)_(?P<suffix>[A-Z]{1,2})$")
#: Provisional identifiers with no binomial: `sp900756035`, `UBA9502`, `CAG-269`.
_PROVISIONAL = re.compile(r"^(?:sp\d{6,}|UBA\d+|CAG-\d+|GGB\d+|SGB\d+|[A-Z]{2,}\d{3,})$")


#: Genera this assay cannot speak to, because they are not bacteria. Read by
#: name rather than by rank: *Candida* looks exactly like a bacterial genus
#: and answering it as one would report a fungus as absent on bacterial
#: evidence. The report covers these in its own sections.
NON_BACTERIAL_GENERA: Final[frozenset[str]] = frozenset({
    # fungi
    "aspergillus", "candida", "cryptococcus", "malassezia", "microsporum",
    "rhodotorula", "saccharomyces", "saprochaete", "trichophyton", "penicillium",
    "mucor", "rhizopus", "geotrichum", "pichia", "debaryomyces", "trichosporon",
    # single-celled parasites
    "blastocystis", "cryptosporidium", "giardia", "entamoeba", "cyclospora",
    "dientamoeba", "toxoplasma",
})

#: Suffixes that mark a rank above genus. A phylum is not an organism this
#: table can count; its members are listed individually.
_ABOVE_GENUS_SUFFIXES: Final[tuple[str, ...]] = (
    "ota", "aceae", "bacteria", "mycota", "ales", "idae",
)


def _above_genus(name: str) -> bool:
    return _norm(name).endswith(_ABOVE_GENUS_SUFFIXES)


def _is_genus_like(name: str) -> bool:
    """A single capitalised word of letters, as a genus is written."""
    token = name.strip()
    return bool(token) and token[:1].isupper() and token.isalpha()


def _norm(name: str) -> str:
    """One spelling of a name for comparison: spaces, case, no suffix noise."""
    return re.sub(r"\s+", " ", str(name).replace("_", " ")).strip().lower()


def _strip_split(token: str) -> str:
    match = _SPLIT.match(token)
    return match.group("base") if match else token


@dataclass(slots=True)
class Resolution:
    """One searched name, and what this catalogue can say about it."""

    query: str
    resolution_state: str
    detection_state: str
    accepted_name: str | None = None
    matched_via: str | None = None
    rank: str = "species"
    genus: str | None = None
    percent: float | None = None
    percentile: float | None = None
    candidates: tuple[str, ...] = ()
    note: str | None = None

    def __post_init__(self) -> None:
        if self.resolution_state not in RESOLUTION_STATES:
            raise ValueError(f"{self.query}: unknown resolution state {self.resolution_state!r}")
        if self.detection_state not in DETECTION_STATES:
            raise ValueError(f"{self.query}: unknown detection state {self.detection_state!r}")

    def to_json(self) -> dict[str, Any]:
        return {
            "query": self.query,
            "resolution_state": self.resolution_state,
            "detection_state": self.detection_state,
            "accepted_name": self.accepted_name,
            "matched_via": self.matched_via,
            "rank": self.rank,
            "genus": self.genus,
            # An organism searched for and not found was measured, and the
            # measurement is nought. Writing null here would leave the
            # renderer to invent the zero it draws, which puts a number on
            # the page that is not in the file - and a client reading this
            # file could not tell the absence apart from an unasked question,
            # which is the whole distinction the detection state carries.
            "percent": (
                0.0
                if self.percent is None
                and self.detection_state == "no_supported_detection"
                else self.percent
            ),
            "percentile": self.percentile,
            "candidates": list(self.candidates),
            "note": self.note,
            "reminder": (
                "Resolving a name is not a detection: the resolution state says what this "
                "catalogue calls the organism, the detection state says whether it is here."
            ),
        }


@dataclass(slots=True)
class Catalogue:
    """The run's organisms, indexed by every name they are known under."""

    rows: tuple[Mapping[str, Any], ...]
    by_name: dict[str, Mapping[str, Any]] = field(default_factory=dict)
    aliases: dict[str, str] = field(default_factory=dict)
    genera: dict[str, list[Mapping[str, Any]]] = field(default_factory=dict)
    #: Candidate calls competitive confirmation rejected, by name: a lookup
    #: for one of these says so rather than "not detected".
    rejected: dict[str, Mapping[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_inventory(cls, organisms: Sequence[Mapping[str, Any]],
                       rejected: Sequence[Mapping[str, Any]] = ()) -> Catalogue:
        cat = cls(rows=tuple(organisms))
        for rec in rejected:
            for name in (rec.get("species"), rec.get("gtdb")):
                if name:
                    cat.rejected.setdefault(_norm(str(name).replace("_", " ")), rec)
        for row in organisms:
            accepted = str(row.get("species") or "").replace("_", " ").strip()
            if not accepted:
                continue
            cat.by_name[_norm(accepted)] = row
            # Every other name this organism is known under, from the
            # inventory itself. No alias is invented: the GTDB name, the
            # former name, each member of a sibling complex ("Blautia_A
            # wexlerae / wexlerae_B" answers to both), and the accepted name
            # without GTDB's genus/epithet split suffixes ("Blautia_A
            # wexlerae" answers to "Blautia wexlerae").
            for key in ("gtdb", "formerly"):
                alias = row.get(key)
                if not alias:
                    continue
                alias = str(alias)
                if " / " in alias:
                    genus, _, rest = alias.partition(" ")
                    for epithet in rest.split(" / "):
                        cat.aliases.setdefault(_norm(f"{genus} {epithet}"), accepted)
                        cat.aliases.setdefault(_norm(f"{_strip_split(genus)} {_strip_split(epithet)}"), accepted)
                else:
                    cat.aliases.setdefault(_norm(alias), accepted)
                    cat.aliases.setdefault(_norm(" ".join(_strip_split(t) for t in alias.split(" "))), accepted)
            for sib in row.get("sibling_species") or ():
                cat.aliases.setdefault(_norm(str(sib)), accepted)
            stripped = " ".join(_strip_split(t) for t in accepted.split(" "))
            if stripped != accepted:
                cat.aliases.setdefault(_norm(stripped), accepted)
            genus = str(row.get("gtdb_genus") or row.get("genus") or "").replace("_", " ").strip()
            if genus:
                cat.genera.setdefault(_norm(genus), []).append(row)
                cat.genera.setdefault(_norm(_strip_split(genus)), []).append(row)
        return cat

    def resolve(self, query: str) -> Resolution:
        """Resolve one name to a state. Never guesses, never invents a synonym."""
        raw = str(query).strip()
        # A dual name such as "Agathobacter rectalis/Eubacterium rectale" is
        # two names for one organism. Each half is tried in turn; the first
        # that resolves answers, and the query string is preserved.
        if "/" in raw:
            for half in (h.strip() for h in raw.split("/")):
                if not half:
                    continue
                attempt = self.resolve(half)
                if attempt.resolution_state not in {"not_in_catalogue", "unresolvable"}:
                    return Resolution(
                        query=raw, resolution_state=attempt.resolution_state,
                        detection_state=attempt.detection_state,
                        accepted_name=attempt.accepted_name,
                        matched_via=f"matched on {half!r}", rank=attempt.rank,
                        genus=attempt.genus, percent=attempt.percent,
                        percentile=attempt.percentile, candidates=attempt.candidates,
                        note=attempt.note,
                    )
        key = _norm(raw)
        tokens = raw.replace("_", " ").split()

        # 1. Exact identity in this run's catalogue.
        row = self.by_name.get(key)
        if row is not None:
            return self._detected(raw, row, "exact_identity", matched_via="accepted name")

        # 2. A name the inventory itself records as an alias.
        accepted = self.aliases.get(key)
        if accepted is not None:
            row = self.by_name.get(_norm(accepted))
            if row is not None:
                return self._detected(
                    raw, row, "accepted_synonym",
                    matched_via=f"recorded as an earlier name for {accepted}",
                )

        # 3. A GTDB split of a name we hold, or the reverse.
        stripped = " ".join(_strip_split(t) for t in tokens)
        row = self.by_name.get(_norm(stripped))
        if row is not None and _norm(stripped) != key:
            return self._detected(
                raw, row, "accepted_synonym",
                matched_via="the same organism without its GTDB split suffix",
            )

        # 3b. A candidate that competitive confirmation rejected: the reads
        # belong to a relative. Said plainly, not "not detected".
        rej = self.rejected.get(key) or self.rejected.get(_norm(stripped))
        if rej is not None:
            reason = str(rej.get("confidence_basis") or "").replace("competitive confirmation: ", "")
            return Resolution(
                query=raw, resolution_state="exact_identity", detection_state="no_supported_detection",
                accepted_name=str(rej.get("species") or raw).replace("_", " "), matched_via="rejected candidate",
                rank="species", genus=str(rej.get("gtdb_genus") or rej.get("genus") or "").replace("_", " ") or None,
                note=f"Called by {', '.join(rej.get('detected_by') or [])} and rejected by competitive confirmation: {reason}",
            )

        # 4. A provisional cluster identifier.
        if any(_PROVISIONAL.match(t) for t in tokens):
            row = self.by_name.get(key)
            return Resolution(
                query=raw,
                resolution_state="cluster_correspondence",
                detection_state="no_supported_detection",
                accepted_name=raw,
                rank="cluster",
                note=(
                    "A provisional genome-cluster identifier with no binomial. It is a stable "
                    "identifier, not a name awaiting resolution."
                ),
            )

        # 5. A genus or higher rank: answer as a rollup.
        members = self.genera.get(key) or self.genera.get(_norm(_strip_split(raw)))
        if members and len(tokens) == 1:
            unique = {id(m): m for m in members}.values()
            total = sum(float(m.get("percent") or 0.0) for m in unique)
            return Resolution(
                query=raw,
                resolution_state="rank_rollup",
                detection_state=("supported_detection" if total > 0 else "no_supported_detection"),
                accepted_name=raw,
                rank="genus",
                genus=raw,
                percent=total or None,
                candidates=tuple(sorted(
                    str(m.get("species") or "").replace("_", " ") for m in unique
                )),
                note=f"Genus total over {len(unique)} species detected here.",
            )

        # 6. A valid-looking binomial this assay simply did not detect.
        if len(tokens) >= 2 and tokens[0][:1].isupper():
            genus_members = self.genera.get(_norm(tokens[0])) or self.genera.get(
                _norm(_strip_split(tokens[0]))
            )
            return Resolution(
                query=raw,
                resolution_state="not_in_catalogue",
                detection_state="no_supported_detection",
                accepted_name=raw,
                rank="species",
                genus=tokens[0],
                candidates=tuple(sorted(
                    {str(m.get("species") or "").replace("_", " ") for m in (genus_members or ())}
                )),
                note=(
                    "A recognisable species name with no supported detection in this sample. "
                    "That is a negative result for this assay, not an absence from the catalogue."
                    if genus_members else
                    "No organism of this name or genus was detected in this sample."
                ),
            )

        # 7. A bare genus name. Nothing of it is in this sample, which is the
        # same answer a binomial gets at step 6 and deserves the same words:
        # searched for, not found. Reporting it as "not assessable" made the
        # reader's question depend on whether they typed one word or two -
        # "Klebsiella pneumoniae" read as absent while "Klebsiella" read as
        # unanswered, on identical evidence.
        if len(tokens) == 1 and _is_genus_like(raw) and not _above_genus(raw):
            if _norm(raw) in NON_BACTERIAL_GENERA:
                return Resolution(
                    query=raw,
                    resolution_state="unresolvable",
                    detection_state="not_assessable",
                    note=(
                        "Not a bacterium, so a bacterial catalogue cannot answer it. This "
                        "report covers these elsewhere; absence here is not absence from you."
                    ),
                )
            return Resolution(
                query=raw,
                resolution_state="not_in_catalogue",
                detection_state="no_supported_detection",
                accepted_name=raw,
                rank="genus",
                genus=raw,
                note=(
                    "A recognisable genus with no member detected in this sample. That is a "
                    "negative result for this assay, not a question left unasked."
                ),
            )

        return Resolution(
            query=raw,
            resolution_state="unresolvable",
            detection_state="not_assessable",
            note=(
                "A rank above the species level this counts; the organisms beneath it are "
                "listed individually."
                if _above_genus(raw)
                else "This string could not be resolved to any organism concept."
            ),
        )

    def _detected(
        self, query: str, row: Mapping[str, Any], state: str, *, matched_via: str
    ) -> Resolution:
        # The primary lane's percent can be zero for an organism the second
        # lane quantified: the catalogue page shows it at its secondary
        # abundance, and a lookup that called the same organism "not found"
        # would contradict the page a reader just read.
        primary = row.get("percent")
        secondary = row.get("secondary_percent")
        percent = float(primary) if primary else None
        if not percent and secondary:
            percent = float(secondary)
        return Resolution(
            query=query,
            resolution_state=state,
            detection_state=(
                "supported_detection" if percent and percent > 0
                else "below_quantification_limit"
            ),
            accepted_name=str(row.get("species") or "").replace("_", " "),
            matched_via=matched_via,
            rank="species",
            genus=str(row.get("gtdb_genus") or row.get("genus") or "").replace("_", " ") or None,
            percent=percent,
            percentile=row.get("percentile"),
            note=(None if primary else
                  "Quantified by the second detection lane; the primary lane reports zero here."),
        )


def load_fixture(path: Path | None = None) -> dict[str, list[str]]:
    """The 332 resolution labels and the 34 sentinels."""
    raw = yaml.safe_load((path or FIXTURE_FILE).read_text(encoding="utf-8")) or {}
    return {
        "resolution_fixture": list(raw.get("resolution_fixture") or ()),
        "sentinels": list(raw.get("sentinels") or ()),
    }


def rollup(
    organisms: Sequence[Mapping[str, Any]], level: str = "genus"
) -> list[dict[str, Any]]:
    """Genus or phylum totals that disclose their membership.

    Section 7.1: a rollup must say what it summed and what it excluded, and a
    genus total may not sit beside its own species as an independent slice.
    """
    if level != "genus":
        raise ValueError(f"unsupported rollup level {level!r}")
    groups: dict[str, list[Mapping[str, Any]]] = {}
    unassigned: list[Mapping[str, Any]] = []
    for row in organisms:
        genus = str(row.get("gtdb_genus") or row.get("genus") or "").replace("_", " ").strip()
        if genus:
            groups.setdefault(genus, []).append(row)
        else:
            unassigned.append(row)
    out: list[dict[str, Any]] = []
    for genus, members in sorted(groups.items()):
        total = sum(float(m.get("percent") or 0.0) for m in members)
        out.append({
            "genus": genus,
            "percent": total,
            "n_members": len(members),
            "members": sorted(str(m.get("species") or "").replace("_", " ") for m in members),
            "excluded": [],
            "note": (
                "A genus total over the species listed. It is a sum of those rows, not an "
                "additional organism, and must not be charted beside them."
            ),
        })
    out.sort(key=lambda r: -r["percent"])
    if unassigned:
        out.append({
            "genus": "(no genus assigned)",
            "percent": sum(float(m.get("percent") or 0.0) for m in unassigned),
            "n_members": len(unassigned),
            "members": sorted(str(m.get("species") or "").replace("_", " ") for m in unassigned),
            "excluded": [],
            "note": "Organisms whose lane gives no genus. Shown rather than dropped.",
        })
    return out


def role_composition(
    organisms: Sequence[Mapping[str, Any]], verdicts: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Share of the community by the existing interpretation classes.

    Uses the registry's own classes with their existing meanings. Unassigned
    mass stays visible: shrinking the denominator until the categories look
    complete is how an "unknown" slice disappears.
    """
    by_species = {
        _norm(str(v.get("species") or "")): str(v.get("class") or "unknown")
        for v in verdicts
    }
    buckets: dict[str, float] = {}
    unassigned = 0.0
    for row in organisms:
        percent = float(row.get("percent") or 0.0)
        if percent <= 0:
            continue
        name = _norm(str(row.get("species") or ""))
        role = by_species.get(name)
        if role is None:
            unassigned += percent
        else:
            buckets[role] = buckets.get(role, 0.0) + percent
    total = sum(buckets.values()) + unassigned
    return {
        "classes": [
            {"role": role, "percent": value,
             "share_of_named": (value / total if total else None)}
            for role, value in sorted(buckets.items(), key=lambda kv: -kv[1])
        ],
        "unassigned_percent": unassigned,
        "unassigned_share": (unassigned / total if total else None),
        "total_percent": total,
        "note": (
            "Uses the existing organism classes unchanged. Mass with no assigned class is "
            "shown as unassigned rather than removed from the denominator, and 'unknown' "
            "is a statement about the evidence, not about harm."
        ),
    }


def named_target_coverage(
    catalogue: Catalogue, targets: Sequence[Mapping[str, str]]
) -> list[dict[str, Any]]:
    """Every named screening target, resolved, whether or not it was found.

    §7.2. A reader who looks up *Klebsiella oxytoca* needs an answer
    whether or not it is in their sample, and "searched for and not found"
    is an answer. Leaving the row out would make an absent organism
    indistinguishable from one nobody looked for.
    """
    out: list[dict[str, Any]] = []
    for target in targets:
        name = str(target.get("concept") or "").strip()
        if not name:
            continue
        resolution = catalogue.resolve(name).to_json()
        out.append({
            "capability_id": target.get("id"),
            "query": name,
            **{k: v for k, v in resolution.items() if k != "query"},
        })
    return out


def build(
    results: Mapping[str, Any], *, fixture: Mapping[str, list[str]] | None = None,
    named_targets: Sequence[Mapping[str, str]] = (),
) -> dict[str, Any]:
    """The whole A10 view."""
    inventory = results.get("organism_inventory") or {}
    organisms = list(inventory.get("organisms") or ())
    catalogue = Catalogue.from_inventory(organisms, rejected=list(inventory.get("rejected") or ()))
    fixture = fixture or load_fixture()

    sentinels = [catalogue.resolve(name).to_json() for name in fixture["sentinels"]]
    verdicts = (results.get("organism_verdicts") or {}).get("verdicts") or []

    detected = sum(1 for s in sentinels if s["detection_state"] == "supported_detection")
    coverage = named_target_coverage(catalogue, named_targets)
    found = sum(1 for t in coverage if t["detection_state"] == "supported_detection")
    return {
        "schema_version": "openbiota.organism-explorer/1.0",
        "feature_id": "A10",
        "method_id": METHOD,
        "n_organisms": len(organisms),
        "n_genera": len({
            str(r.get("gtdb_genus") or r.get("genus") or "") for r in organisms
        } - {""}),
        "sentinels": sentinels,
        "n_sentinels_detected": detected,
        "named_targets": coverage,
        "n_named_targets": len(coverage),
        "n_named_targets_detected": found,
        "named_target_note": (
            "Every organism this report screens for by name, with what was found. A "
            "target that is listed and not detected was searched for and not found, "
            "which is different from one nobody looked for."
        ),
        "genus_rollup": rollup(organisms),  # every genus; no hidden top-N (spec 0.8.4 §6)
        "role_composition": role_composition(organisms, verdicts),
        "limitations": [
            "Resolving a name is not a detection, and a name absent from this list is not "
            "absent from the catalogue.",
            "A species label does not establish a commercial probiotic strain, toxin "
            "carriage or a disease-causing subtype: those are strain-level questions.",
            "A genus total is a sum of the species listed under it, not an extra organism.",
        ],
    }


def metrics(view: Mapping[str, Any], *, input_fingerprint: str | None = None) -> tuple[ExtensionMetric, ...]:
    """The explorer's countable outputs."""
    fp = input_fingerprint or fingerprint(view.get("n_organisms"), view.get("n_genera"))
    role = view.get("role_composition") or {}
    out = [
        ExtensionMetric(
            metric_id="ext083.explorer.organisms_listed",
            label="Organisms in the complete explorer",
            kind="taxon_abundance",
            state="measured",
            value=float(view.get("n_organisms") or 0),
            unit="organisms",
            denominator="all supported detections, unfiltered",
            direction="descriptive",
            method_id=METHOD,
            input_fingerprint=fp,
            group="explorer",
            feature_id="A10",
            analytical_confidence="supported",
            source_ids=("E08", "E09"),
            limitations=(
                "A count of what was detected, not a measure of quality: a larger number is "
                "not a better result.",
            ),
            extra={"n_genera": view.get("n_genera")},
        ),
    ]
    if role.get("unassigned_share") is not None:
        out.append(ExtensionMetric(
            metric_id="ext083.explorer.interpretation_coverage",
            label="Share of the community with an interpretation",
            kind="community_composition",
            state="measured",
            value=1.0 - float(role["unassigned_share"]),
            unit="fraction",
            denominator="named abundance",
            direction="descriptive",
            method_id=METHOD,
            input_fingerprint=fp,
            group="explorer",
            feature_id="A10",
            analytical_confidence="supported",
            source_ids=("E09",),
            limitations=(
                "How much of the community the registry can interpret. Unassigned mass is "
                "shown, not removed from the denominator.",
            ),
            extra=role,
        ))
    return tuple(out)


def resolution_audit(
    catalogue: Catalogue, names: Iterable[str]
) -> dict[str, Any]:
    """Resolve every fixture label and report the states reached.

    Section 16.4: every label must produce an explicit state and none may be
    silently dropped.
    """
    results = [catalogue.resolve(name) for name in names]
    by_state: dict[str, int] = {}
    for r in results:
        by_state[r.resolution_state] = by_state.get(r.resolution_state, 0) + 1
    return {
        "n_queried": len(results),
        "by_resolution_state": dict(sorted(by_state.items())),
        "n_detected": sum(1 for r in results if r.detection_state == "supported_detection"),
        "unresolvable": [r.query for r in results if r.resolution_state == "unresolvable"],
        "results": [r.to_json() for r in results],
    }


__all__ = [
    "DETECTION_STATES",
    "FIXTURE_FILE",
    "RESOLUTION_STATES",
    "Catalogue",
    "Resolution",
    "build",
    "load_fixture",
    "named_target_coverage",
    "metrics",
    "resolution_audit",
    "role_composition",
    "rollup",
]
