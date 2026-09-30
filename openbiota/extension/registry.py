"""The capability registry: what this release owes, and what it actually does.

BUILD_SPEC_v0.8.3 section 16 is normative for completeness, and section 15.7
is blunt about why this module exists: "An engineering TODO cannot be
disguised as a scientific limitation." So the states here keep apart the
reasons a feature has no number:

* ``not_implemented``   — owed, and nobody has built it. An honest debt.
* ``missing_input``     — built, and this sample lacks the data to run it.
* ``access_restricted`` — built, and the reference data needs an agreement.
* ``external_assay``    — no sequence can answer; a laboratory result can.
* ``unsupported_by_assay`` — this assay cannot answer at any depth.
* ``scientifically_unresolved`` — the biology is genuinely not established.

Only the last is a statement about science. Collapsing these into one
"unavailable" flag is what lets an unfinished build read as a limit of
knowledge, and the coverage report refuses to do it.

Readiness is tracked on six independent axes (section 15.7), because a URL
that returns 200 is not an installed, validated database.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.extension import SPEC_VERSION

#: How a required capability is implemented. `legacy_surfaced` is a promise
#: that the underlying value is untouched.
IMPLEMENTATION_STATES: Final[frozenset[str]] = frozenset({
    "legacy_surfaced",
    "extension_computed",
    "evidence_card",
    "external_assay",
    "not_implemented",
})

#: Why a capability has no value for a given sample. Distinct from *how* it is
#: implemented: a fully built module still has no number when the input is
#: absent, and that is not the same as never having been written.
UNAVAILABLE_REASONS: Final[frozenset[str]] = frozenset({
    "computed",
    "not_implemented",
    "missing_input",
    "access_restricted",
    "external_assay",
    "unsupported_by_assay",
    "scientifically_unresolved",
    "not_applicable",
})

#: Reasons that are an engineering debt rather than a property of the world.
ENGINEERING_DEBT: Final[frozenset[str]] = frozenset({"not_implemented"})

#: The six readiness axes of section 15.7, in order of increasing strength.
#: Each is independent: `retrievable` says a URL answered, `bytes_verified`
#: says the bytes matched a recorded hash, and `application_validated` says
#: this application produces a tested result from it.
READINESS_AXES: Final[tuple[str, ...]] = (
    "source_identified",
    "retrievable",
    "bytes_verified",
    "registry_validated",
    "example_reproduced",
    "application_validated",
)

DATA_DIR: Final = Path(__file__).resolve().parents[2] / "extension"
MANIFEST_FILE: Final = DATA_DIR / "capability_manifest.yaml"


class RegistryError(ValueError):
    """A capability manifest that does not describe a possible state."""


@dataclass(slots=True)
class Readiness:
    """Where one capability or asset actually stands, on six axes."""

    source_identified: bool = False
    retrievable: bool = False
    bytes_verified: bool = False
    registry_validated: bool = False
    example_reproduced: bool = False
    application_validated: bool = False
    note: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {axis: getattr(self, axis) for axis in READINESS_AXES} | {"note": self.note}

    @property
    def highest_axis(self) -> str | None:
        reached = [a for a in READINESS_AXES if getattr(self, a)]
        return reached[-1] if reached else None


@dataclass(slots=True)
class Capability:
    """One row of the normative manifest, with what implements it."""

    capability_id: str
    kind: str                      # feature | measurement_view | function_view
    concept: str
    binding: str
    implementation: str = "not_implemented"
    module: str | None = None
    metric_ids: tuple[str, ...] = ()
    legacy_metric_ids: tuple[str, ...] = ()
    shares_measurement_with: tuple[str, ...] = ()
    report_section: str | None = None
    source_ids: tuple[str, ...] = ()
    tests: tuple[str, ...] = ()
    readiness: Readiness = field(default_factory=Readiness)
    note: str | None = None

    def __post_init__(self) -> None:
        if self.implementation not in IMPLEMENTATION_STATES:
            raise RegistryError(
                f"{self.capability_id}: unknown implementation state {self.implementation!r}"
            )
        # A legacy view must say which legacy metric it surfaces, or the
        # promise that the value is unchanged is unverifiable.
        if self.implementation == "legacy_surfaced" and not self.legacy_metric_ids:
            raise RegistryError(
                f"{self.capability_id}: 'legacy_surfaced' must name the legacy metric it shows"
            )
        # A computed capability must name what it produces.
        if self.implementation == "extension_computed" and not self.metric_ids:
            raise RegistryError(
                f"{self.capability_id}: 'extension_computed' must name its metric IDs"
            )
        if self.implementation != "not_implemented" and not self.module:
            raise RegistryError(
                f"{self.capability_id}: an implemented capability must name its module"
            )

    @property
    def implemented(self) -> bool:
        return self.implementation != "not_implemented"

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.capability_id,
            "kind": self.kind,
            "concept": self.concept,
            "binding": self.binding,
            "implementation": self.implementation,
            "implemented": self.implemented,
            "module": self.module,
            "metric_ids": list(self.metric_ids),
            "legacy_metric_ids": list(self.legacy_metric_ids),
            "shares_measurement_with": list(self.shares_measurement_with),
            "report_section": self.report_section,
            "source_ids": list(self.source_ids),
            "tests": list(self.tests),
            "readiness": self.readiness.to_json(),
            "note": self.note,
        }


@dataclass(slots=True)
class CapabilityRegistry:
    """Every capability this release owes, keyed by ID."""

    capabilities: dict[str, Capability] = field(default_factory=dict)
    spec_version: str = SPEC_VERSION

    # -- construction ------------------------------------------------------ #

    @classmethod
    def load(cls, path: Path | None = None) -> CapabilityRegistry:
        """Read the normative manifest."""
        path = path or MANIFEST_FILE
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        reg = cls(spec_version=str(raw.get("spec_version") or SPEC_VERSION))
        for kind, key in (
            ("feature", "features"),
            ("measurement_view", "measurement_views"),
            ("function_view", "function_views"),
        ):
            for row in raw.get(key) or ():
                cid = str(row["id"])
                if cid in reg.capabilities:
                    raise RegistryError(f"duplicate capability {cid}")
                concept = str(
                    row.get("concept") or row.get("addition") or ""
                )
                binding = str(
                    row.get("binding") or row.get("required_output") or ""
                )
                reg.capabilities[cid] = Capability(
                    capability_id=cid,
                    kind=kind,
                    concept=concept,
                    binding=binding,
                    implementation=str(row.get("implementation") or "not_implemented"),
                    module=row.get("module"),
                    metric_ids=tuple(row.get("metric_ids") or ()),
                    legacy_metric_ids=tuple(row.get("legacy_metric_ids") or ()),
                    shares_measurement_with=tuple(row.get("shares_measurement_with") or ()),
                    report_section=row.get("report_section"),
                    source_ids=tuple(row.get("source_ids") or ()),
                    tests=tuple(row.get("tests") or ()),
                    readiness=Readiness(**(row.get("readiness") or {})),
                    note=row.get("note"),
                )
        if not reg.capabilities:
            raise RegistryError(f"{path} contains no capabilities")
        return reg

    # -- views ------------------------------------------------------------- #

    def of_kind(self, kind: str) -> tuple[Capability, ...]:
        return tuple(c for c in self.capabilities.values() if c.kind == kind)

    @property
    def features(self) -> tuple[Capability, ...]:
        return self.of_kind("feature")

    def missing(self) -> tuple[Capability, ...]:
        """Capabilities still owed. Named, never aggregated into a percentage."""
        return tuple(c for c in self.capabilities.values() if not c.implemented)

    def require_ids(self, ids: Iterable[str]) -> None:
        """Assert the manifest carries every ID the specification lists."""
        missing = sorted(set(map(str, ids)) - set(self.capabilities))
        if missing:
            raise RegistryError(
                f"the capability manifest is missing {len(missing)} required ID(s): {missing[:12]}"
            )

    def coverage(self, *, sample: str | None = None) -> dict[str, Any]:
        """The `capability_coverage.json` payload of section 13.1.

        Every required capability, whether it is computed, data-limited,
        externally measured or not applicable - and, separately, whether the
        reason is an engineering debt. A reader of this file can tell the
        difference between "we have not built it" and "the sample cannot
        answer it", which is the whole point.
        """
        by_state: dict[str, int] = {}
        by_kind: dict[str, dict[str, int]] = {}
        for cap in self.capabilities.values():
            by_state[cap.implementation] = by_state.get(cap.implementation, 0) + 1
            per = by_kind.setdefault(cap.kind, {})
            per[cap.implementation] = per.get(cap.implementation, 0) + 1
        debt = sorted(c.capability_id for c in self.missing())
        readiness_counts = {
            axis: sum(1 for c in self.capabilities.values() if getattr(c.readiness, axis))
            for axis in READINESS_AXES
        }
        return {
            "schema_version": "openbiota.capability-coverage/1.0",
            "spec_version": self.spec_version,
            "sample": sample,
            "n_capabilities": len(self.capabilities),
            "by_implementation_state": dict(sorted(by_state.items())),
            "by_kind": {k: dict(sorted(v.items())) for k, v in sorted(by_kind.items())},
            "readiness_axis_counts": readiness_counts,
            "engineering_debt": debt,
            "n_engineering_debt": len(debt),
            "capabilities": [c.to_json() for c in sorted(
                self.capabilities.values(), key=lambda c: (c.kind, c.capability_id)
            )],
        }

    def repeats(self) -> dict[str, tuple[str, ...]]:
        """Concepts that deliberately share one canonical measurement.

        Section 16.1: M064/M090 are the same organism in two contexts and
        M068/M073 the same genus total. They must not be counted as four
        independent assays, and this is where that is recorded.
        """
        out: dict[str, tuple[str, ...]] = {}
        for cap in self.capabilities.values():
            if cap.shares_measurement_with:
                out[cap.capability_id] = cap.shares_measurement_with
        return dict(sorted(out.items()))

    def unique_measurement_count(self) -> int:
        """Capabilities, minus those that are another's second view.

        Section 12.3 requires per-section counts to distinguish unique
        measurements from repeated contextual views.
        """
        seen: set[frozenset[str]] = set()
        for cap in self.capabilities.values():
            if cap.kind != "measurement_view":
                continue
            group = frozenset({cap.capability_id, *cap.shares_measurement_with})
            seen.add(group)
        return len(seen)


def load(path: Path | None = None) -> CapabilityRegistry:
    """Convenience loader."""
    return CapabilityRegistry.load(path)


__all__ = [
    "ENGINEERING_DEBT",
    "IMPLEMENTATION_STATES",
    "MANIFEST_FILE",
    "READINESS_AXES",
    "UNAVAILABLE_REASONS",
    "Capability",
    "CapabilityRegistry",
    "Readiness",
    "RegistryError",
    "load",
]
