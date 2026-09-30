"""Person, sample and material lineage (spec §4.3, §5.2, §9.1).

A donor candidate is a *person × material unit*, not a name. Two samples from
one person are replicates or longitudinal evidence, never two donors and never
extra diversity; the same sample handed in twice is one material. Every
identity conclusion here is a computational condition, not a disease
hypothesis: a mismatch excludes the *record*, never the person.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, Literal

Role = Literal["recipient", "donor"]

#: Facts an older single-sample report assumed rather than measured. They are
#: carried as `assumed` and can never be promoted to verified metadata (§5.2).
ASSUMED_FACT_NOTE: Final = (
    "assumed by the source report because the fact was not supplied; not participant metadata"
)


def content_id(prefix: str, payload: Any) -> str:
    """A stable, content-addressed ID: same inputs always give the same ID."""
    blob = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return f"{prefix}-{hashlib.sha256(blob.encode('utf-8')).hexdigest()[:16]}"


@dataclass(slots=True)
class MetadataFact:
    """One fact about a person or material, with where it came from."""

    field_name: str
    value: Any
    unit: str | None
    source: str
    date: str | None
    verification: Literal["verified", "self_reported", "assumed", "unknown"]
    assumed_by_report: bool = False
    note: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "field": self.field_name,
            "value": self.value,
            "unit": self.unit,
            "source": self.source,
            "date": self.date,
            "verification": self.verification,
            "assumed_by_source_report": self.assumed_by_report,
            "note": self.note,
        }


@dataclass(slots=True)
class MaterialUnit:
    """One sequenced material from one person (spec §5.2)."""

    material_id: str
    person_id: str
    sample_id: str
    role: Role
    input_root: str
    material_type: str = "stool_dna_sequencing"
    body_site: str = "stool"
    assay_type: str = "shotgun_metagenome_dna"
    collection_time: str | None = None
    donation_id: str | None = None
    lot_id: str | None = None
    related_person_group: str | None = None
    replicate_group: str | None = None
    processing_batch: str | None = None
    identity_status: str = "declared_unverified"
    identity_notes: list[str] = field(default_factory=list)
    metadata_facts: list[MetadataFact] = field(default_factory=list)
    storage_processing_facts: dict[str, Any] = field(default_factory=dict)
    input_artifact_ids: list[str] = field(default_factory=list)

    def to_json(self) -> dict[str, Any]:
        return {
            "material_id": self.material_id,
            "person_id": self.person_id,
            "sample_id": self.sample_id,
            "role": self.role,
            "material_type": self.material_type,
            "body_site": self.body_site,
            "assay_type": self.assay_type,
            "collection_time": self.collection_time,
            "donation_id": self.donation_id,
            "lot_id": self.lot_id,
            "related_person_group": self.related_person_group,
            "replicate_group": self.replicate_group,
            "processing_batch": self.processing_batch,
            "identity_status": self.identity_status,
            "identity_notes": list(self.identity_notes),
            "metadata_facts": [f.to_json() for f in self.metadata_facts],
            "storage_processing_facts": dict(self.storage_processing_facts),
            "input_artifact_ids": list(self.input_artifact_ids),
        }


def artifact_id(path: Path) -> str:
    """Digest-addressed ID for an input file, so a run can be reproduced."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return f"artifact-{h.hexdigest()[:16]}"


def make_material(
    *,
    person_id: str,
    sample_id: str,
    role: Role,
    input_root: Path,
    run_facts: dict[str, Any],
    artifacts: list[str],
) -> MaterialUnit:
    """Build a material unit from a loaded results bundle.

    `collection_time`, `donation_id`, `lot_id` and viability are genuinely
    unknown from sequencing outputs, so they stay null rather than being
    back-filled from the report date (spec §3, §6.3).
    """
    material_id = content_id("material", {"person": person_id, "sample": sample_id, "artifacts": artifacts})
    facts: list[MetadataFact] = []
    started = run_facts.get("started")
    if started:
        facts.append(
            MetadataFact(
                field_name="sequencing_run_started",
                value=started,
                unit=None,
                source="pipeline_run_log",
                date=str(started)[:10],
                verification="verified",
                note="when the analysis ran, not when the material was collected",
            )
        )
    return MaterialUnit(
        material_id=material_id,
        person_id=person_id,
        sample_id=sample_id,
        role=role,
        input_root=str(input_root),
        processing_batch=str(run_facts.get("reference fingerprint") or "") or None,
        metadata_facts=facts,
        input_artifact_ids=list(artifacts),
        identity_status="sample_id_matches_declared_input",
    )


@dataclass(slots=True)
class IdentityConflict:
    """A computational blocker: wrong person, wrong material or bad integrity."""

    material_id: str | None
    sample_id: str | None
    person_id: str | None
    kind: str
    detail: str

    def to_json(self) -> dict[str, Any]:
        return {
            "material_id": self.material_id,
            "sample_id": self.sample_id,
            "person_id": self.person_id,
            "kind": self.kind,
            "detail": self.detail,
            "scope": "analytical_record_only",
            "note": "an identity conflict excludes this record, not the person (rule X03)",
        }


def dedupe_materials(materials: list[MaterialUnit]) -> tuple[list[MaterialUnit], list[str]]:
    """Collapse identical inputs; keep separate donations from one person apart.

    Duplicate submission of the same sample must not create a second donor or
    inflate coverage or diversity (V7-004). Distinct samples from one person
    stay as distinct materials sharing a person ID (V7-005).
    """
    seen: dict[tuple[str, str], MaterialUnit] = {}
    notes: list[str] = []
    for m in materials:
        key = (m.person_id, m.sample_id)
        if key in seen:
            notes.append(
                f"duplicate input for {m.person_id}/{m.sample_id} collapsed into one material unit"
            )
            continue
        seen[key] = m
    kept = list(seen.values())
    by_person: dict[str, list[MaterialUnit]] = {}
    for m in kept:
        by_person.setdefault(m.person_id, []).append(m)
    for person, group in by_person.items():
        if len(group) > 1:
            for m in group:
                m.replicate_group = f"person:{person}"
                m.identity_notes.append(
                    "several materials from this person: replicate or longitudinal evidence, "
                    "not several donors"
                )
    return kept, notes
