"""The embedded research register of BUILD_SPEC_v0.8.3 section 15.

Every new claim, metric and action in the extension resolves to a source ID
here. Two rules shape the module.

**A reachable URL is not an installed database.** Readiness is six independent
flags (section 15.6), and only the asset verifier may set the later ones,
after it has fetched bytes and matched a recorded hash. Nothing in this module
sets `bytes_verified` from a successful HTTP request.

**Terms are not one boolean.** Several of these sources carry a different
licence for their code than for their data, and one bundles an
academic-noncommercial binary inside an otherwise open wrapper. Flattening
that into "open: true" is how a noncommercial artefact becomes a production
dependency, so licences are kept as a list of the specification's own strings
and restrictions are carried verbatim.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.extension import SPEC_VERSION
from openbiota.extension.registry import READINESS_AXES, Readiness

DATA_DIR: Final = Path(__file__).resolve().parents[2] / "extension"
REGISTER_FILE: Final = DATA_DIR / "source_register.yaml"

#: The register's groups, and what each is for.
GROUPS: Final[Mapping[str, str]] = {
    "functional": "F — sequence, pathway and enzyme sources for the functional lane",
    "ecology": "E — ecology, reference-cohort and interpretation resources",
    "context": "C — host-signalling and symptom/condition evidence",
    "benchmark": "D — profiler, AMR and validation benchmarks",
    "intervention": "R — intervention, modelling and dietary-response research",
    "synbiotic": "S — strain–substrate, product and simulation assets",
}

#: Sources whose terms forbid, or condition, bundling them in production. Kept
#: as an explicit set so a core dependency can be tested against it: section
#: 13.2 requires the report to work without any of these.
NON_CORE_TERMS: Final[frozenset[str]] = frozenset({
    "noncommercial", "commercial use paid", "academic-login", "requires verify",
    "request-only", "request/approval", "approved access", "controlled", "restricted",
    "require approved access", "access-restricted",
})


class SourceError(ValueError):
    """A source record that would misrepresent its terms or its readiness."""


@dataclass(slots=True)
class Source:
    """One primary source, with its assets, terms and actual readiness."""

    source_id: str
    group: str
    note: str
    urls: tuple[str, ...] = ()
    dois: tuple[str, ...] = ()
    accessions: tuple[str, ...] = ()
    licence: tuple[str, ...] = ()
    access_restrictions: tuple[str, ...] = ()
    sha256: tuple[str, ...] = ()
    md5: tuple[str, ...] = ()
    code_commit: tuple[str, ...] = ()
    caution: str | None = None
    correction: str | None = None
    readiness: Readiness = field(default_factory=Readiness)

    def __post_init__(self) -> None:
        if self.group not in GROUPS:
            raise SourceError(f"{self.source_id}: unknown group {self.group!r}")
        if not self.note:
            raise SourceError(f"{self.source_id}: a source must say what it is for")
        if not (self.urls or self.dois or self.accessions):
            raise SourceError(
                f"{self.source_id}: a source with no URL, DOI or accession cannot be located"
            )
        # Readiness beyond retrieval requires something to have been checked.
        if self.readiness.bytes_verified and not (self.sha256 or self.md5):
            raise SourceError(
                f"{self.source_id}: bytes_verified with no recorded checksum; an HTTP 200 "
                "is not a verified download"
            )

    @property
    def core_eligible(self) -> bool:
        """Whether this source may be a required production dependency.

        Section 13.2: core dependencies must remain usable without a
        restricted commercial knowledge base, so anything carrying a
        conditional term is adapter-only.
        """
        return not any(
            r.lower() in NON_CORE_TERMS for r in (*self.access_restrictions, *self.licence)
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.source_id,
            "group": self.group,
            "note": self.note,
            "caution": self.caution,
            "correction": self.correction,
            "urls": list(self.urls),
            "dois": list(self.dois),
            "accessions": list(self.accessions),
            "licence": list(self.licence),
            "access_restrictions": list(self.access_restrictions),
            "core_eligible": self.core_eligible,
            "sha256": list(self.sha256),
            "md5": list(self.md5),
            "code_commit": list(self.code_commit),
            "readiness": self.readiness.to_json(),
        }


@dataclass(slots=True)
class SourceRegister:
    """Every source the extension may cite, keyed by ID."""

    sources: dict[str, Source] = field(default_factory=dict)
    spec_version: str = SPEC_VERSION

    @classmethod
    def load(cls, path: Path | None = None) -> SourceRegister:
        path = path or REGISTER_FILE
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        reg = cls(spec_version=str(raw.get("spec_version") or SPEC_VERSION))
        for group in GROUPS:
            for row in raw.get(group) or ():
                sid = str(row["id"])
                if sid in reg.sources:
                    raise SourceError(f"duplicate source {sid}")
                reg.sources[sid] = Source(
                    source_id=sid,
                    group=group,
                    note=str(row.get("note") or ""),
                    urls=tuple(row.get("urls") or ()),
                    dois=tuple(row.get("dois") or ()),
                    accessions=tuple(row.get("accessions") or ()),
                    licence=tuple(row.get("licence") or ()),
                    access_restrictions=tuple(row.get("access_restrictions") or ()),
                    sha256=tuple(row.get("sha256") or ()),
                    md5=tuple(row.get("md5") or ()),
                    code_commit=tuple(row.get("code_commit") or ()),
                    caution=row.get("caution"),
                    correction=row.get("correction"),
                    readiness=Readiness(**(row.get("readiness") or {})),
                )
        if not reg.sources:
            raise SourceError(f"{path} contains no sources")
        return reg

    def require(self, ids: Iterable[str]) -> None:
        """Assert every cited ID exists. Called wherever a metric is built."""
        missing = sorted(set(map(str, ids)) - set(self.sources))
        if missing:
            raise SourceError(f"cited source ID(s) not in the register: {missing}")

    def of_group(self, group: str) -> tuple[Source, ...]:
        return tuple(s for s in self.sources.values() if s.group == group)

    def restricted(self) -> tuple[Source, ...]:
        """Sources that may not be a required production dependency."""
        return tuple(s for s in self.sources.values() if not s.core_eligible)

    def with_checksums(self) -> tuple[Source, ...]:
        return tuple(s for s in self.sources.values() if s.sha256 or s.md5)

    def audit(self) -> dict[str, Any]:
        """Asset-level provenance state, for the extension manifest."""
        by_group = {g: len(self.of_group(g)) for g in GROUPS}
        axis_counts = {
            axis: sum(1 for s in self.sources.values() if getattr(s.readiness, axis))
            for axis in READINESS_AXES
        }
        return {
            "n_sources": len(self.sources),
            "by_group": by_group,
            "readiness_axis_counts": axis_counts,
            "n_restricted": len(self.restricted()),
            "restricted_ids": sorted(s.source_id for s in self.restricted()),
            "n_with_recorded_checksum": len(self.with_checksums()),
            "n_with_accessions": sum(1 for s in self.sources.values() if s.accessions),
        }

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": "openbiota.source-register/1.0",
            "spec_version": self.spec_version,
            "audit": self.audit(),
            "sources": {k: v.to_json() for k, v in sorted(self.sources.items())},
        }


def load(path: Path | None = None) -> SourceRegister:
    return SourceRegister.load(path)


__all__ = [
    "GROUPS",
    "NON_CORE_TERMS",
    "REGISTER_FILE",
    "Source",
    "SourceError",
    "SourceRegister",
    "load",
]
