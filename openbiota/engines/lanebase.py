"""What every detection lane hands back, in one shape.

Spec 0.8.4 §6 asks for a `DetectionObservation` per lane per organism:
native identifier, rank, status, support metrics, alternatives, the
lane's own abundance value *with its own unit and denominator*, and where
the raw result lives. Lanes normalise differently - MetaPhlAn to percent
of classified, sylph to taxonomic abundance, mOTUs to relative abundance
of mOTUs, Kraken/Bracken to read fractions - so the unit travels with the
number and nothing downstream may add two lanes' numbers together.

`LaneResult` is the JSON record a lane writes into results.json; the
inventory turns its observations into `Lane` objects. A lane that did not
run returns `None` from its runner and is recorded as `not_assessed`, never
as a set of absences.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

STATUSES = ("supported", "provisional", "ambiguous", "not_detected", "not_assessed")


@dataclass(frozen=True)
class Observation:
    native_id: str                 # SGB id, GTDB accession, mOTU id, MGYG id, taxid...
    lineage: str                   # the lane's own lineage string
    rank: str                      # "species" | "genus" | ... (rank of the call)
    status: str = "supported"      # one of STATUSES; the lane's own opinion
    abundance_value: float | None = None
    abundance_unit: str = ""       # e.g. "percent of classified reads"
    denominator: str = ""          # e.g. "classified reads", "mOTU-assigned reads"
    support_metrics: dict[str, Any] = field(default_factory=dict)
    alternatives: tuple[str, ...] = ()

    @property
    def species(self) -> str:
        for seg in self.lineage.replace("|", ";").split(";"):
            seg = seg.strip()
            if seg.startswith("s__"):
                return seg[3:].strip().replace("_", " ")
        return ""

    @property
    def genus(self) -> str:
        for seg in self.lineage.replace("|", ";").split(";"):
            seg = seg.strip()
            if seg.startswith("g__"):
                return seg[3:].strip()
        return ""

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class LaneResult:
    lane_id: str
    tool: str
    tool_version: str
    reference_release_id: str
    taxonomy_release: str
    observations: list[Observation]
    elapsed_s: float
    cached: bool
    command: tuple[str, ...]
    input_sha256: str
    raw_result_uri: str
    #: What the lane could not do, in words a reader can use.
    notes: list[str] = field(default_factory=list)
    #: Anything the lane measured about the whole run (reads assigned etc).
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def species_percent(self) -> dict[str, float]:
        """Species -> the lane's abundance, in the lane's own unit.

        Kept for `inventory.Lane.species`; the unit is on the LaneResult.
        """
        out: dict[str, float] = {}
        for o in self.observations:
            if o.rank == "species" and o.species and o.abundance_value is not None:
                out[o.species] = out.get(o.species, 0.0) + float(o.abundance_value)
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "lane_id": self.lane_id, "tool": self.tool, "tool_version": self.tool_version,
            "reference_release_id": self.reference_release_id, "taxonomy_release": self.taxonomy_release,
            "n_observations": len(self.observations),
            "n_species": sum(1 for o in self.observations if o.rank == "species"),
            "elapsed_s": round(self.elapsed_s, 1), "cached": self.cached,
            "command": list(self.command), "input_sha256": self.input_sha256,
            "raw_result_uri": self.raw_result_uri, "notes": list(self.notes), "summary": dict(self.summary),
            "observations": [o.to_json() for o in self.observations],
        }

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> LaneResult:
        obs = [Observation(**{k: v for k, v in o.items() if k in Observation.__dataclass_fields__})
               for o in data.get("observations") or []]
        for i, o in enumerate(obs):
            if isinstance(o.alternatives, list):
                obs[i] = Observation(**{**o.to_json(), "alternatives": tuple(o.alternatives)})
        return cls(
            lane_id=str(data["lane_id"]), tool=str(data.get("tool", "")), tool_version=str(data.get("tool_version", "")),
            reference_release_id=str(data.get("reference_release_id", "")), taxonomy_release=str(data.get("taxonomy_release", "")),
            observations=obs, elapsed_s=float(data.get("elapsed_s") or 0.0), cached=bool(data.get("cached")),
            command=tuple(data.get("command") or ()), input_sha256=str(data.get("input_sha256", "")),
            raw_result_uri=str(data.get("raw_result_uri", "")), notes=list(data.get("notes") or []),
            summary=dict(data.get("summary") or {}),
        )


def input_sha256(paths: Iterable[Path | None]) -> str:
    """A cheap, stable identity for the input reads: name, size, mtime.

    Hashing 6 GB of FASTQ per lane per run is not worth it; the sequencing
    provider's file is what it is, and the pipeline records the full
    content hash once in the input register.
    """
    h = hashlib.sha256()
    for p in paths:
        if p is None:
            continue
        st = p.stat()
        h.update(f"{p.name}:{st.st_size}:{int(st.st_mtime)}".encode())
    return h.hexdigest()[:24]


def cache_key(*parts: str) -> str:
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()[:16]


def read_cached(path: Path) -> LaneResult | None:
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    result = LaneResult.from_json(data)
    result.cached = True
    return result
