"""The protected-output manifest, and the check that nothing existing moved.

BUILD_SPEC_v0.8.3 section 1 makes one promise about every number already in
the report: for identical resolved inputs it is identical with the extension
on or off. That promise is only worth something if it is machine-checked, so
this module walks a baseline `results.json` and a candidate, compares every
leaf under the protected keys, and reports each difference with its path.

Two things make this harder than a dictionary comparison.

*Volatile run metadata.* A second run of the same sample legitimately differs
in wall-clock timings, log paths and a report date. Those are enumerated
explicitly rather than pattern-matched, so a genuinely changed value can never
hide behind a loose rule.

*Floating point.* A recomputed mean can differ in its last bits without any
change in method. Numbers compare with a relative tolerance, and the tolerance
is part of the recorded result rather than a silent constant.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.extension import SPEC_VERSION

#: Result keys whose contents are protected. Every existing score, profile,
#: organism call, pathogen result and functional panel lives under one of
#: these, and the extension may not change any leaf inside them.
PROTECTED_KEYS: Final[tuple[str, ...]] = (
    "biofilm",
    "community_profile",
    "copri_complex",
    "does_not_measure",
    "extended_catalogue",
    "findings_and_evidence",
    "genome_profile",
    "locus_resolution",
    "measures",
    "metabolite_drivers",
    "microbiome_age",
    "mycobiome",
    "normalisation",
    "organism_inventory",
    "organism_verdicts",
    "panel_bands",
    "panels",
    "pathogens",
    "profile_similarity",
    "qualitative_bands_are_arbitrary",
    "reference_comparison",
    "reference_entries",
    "resolution_census",
    "sequencing_quality",
    "standing_caveats",
    "strain_resolution",
    "subject_context",
    "validated_reference_range",
)

#: Paths that legitimately differ between two runs of the same input. Listed
#: one by one: a pattern such as "anything containing time" would also excuse a
#: changed measurement whose label happens to mention timing.
VOLATILE_PATHS: Final[frozenset[str]] = frozenset({
    "run.total wall clock",
    "run.total_s",
    "run.wall clock",
    "run.elapsed_s",
    "run.started",
    "run.finished",
    "run.report date",
    "run.log",
    "run.host",
    "run.threads",
    "run.command",
    "run.output directory",
    "run.stage timings",
    "run.timings_s",
    "input.path",
    "input.fastq",
    "diagnostics",
    "search.elapsed_s",
    "pathogens.timings_s",
    "mycobiome.timings_s",
    "biofilm.timings_s",
    # Stage timings that live *inside* a protected object. How long a stage
    # took is run metadata wherever it is stored, and it moves by tenths of a
    # second between identical runs. Each is listed by its exact path: a rule
    # like "skip anything ending _s" would also excuse a measurement.
    "extended_catalogue.elapsed_s",
    "genome_profile.elapsed_s",
    "pathogens.candidate_stage.elapsed_s",
    "profile_similarity.elapsed_s",
    "profile_similarity.taxonomic_engine.elapsed_s",
    "sequencing_quality.fastp.elapsed_s",
})

#: Relative tolerance for numeric comparison. A recomputed float may differ in
#: its last bits; a changed calculation will not.
DEFAULT_RTOL: Final = 1e-9


@dataclass(slots=True)
class Difference:
    """One protected value that changed, with enough detail to act on."""

    path: str
    baseline: Any
    candidate: Any
    reason: str

    def to_json(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "baseline": _safe(self.baseline),
            "candidate": _safe(self.candidate),
            "reason": self.reason,
        }


@dataclass(slots=True)
class PreservationReport:
    """The verdict, the differences, and what was actually compared."""

    baseline_path: str
    candidate_path: str
    compared_leaves: int = 0
    protected_keys_present: tuple[str, ...] = ()
    protected_keys_missing: tuple[str, ...] = ()
    differences: list[Difference] = field(default_factory=list)
    rtol: float = DEFAULT_RTOL
    spec_version: str = SPEC_VERSION

    @property
    def preserved(self) -> bool:
        """True when every protected leaf matched and none went missing."""
        return not self.differences and not self.protected_keys_missing

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": "openbiota.preservation-report/1.0",
            "spec_version": self.spec_version,
            "baseline": self.baseline_path,
            "candidate": self.candidate_path,
            "preserved": self.preserved,
            "compared_leaves": self.compared_leaves,
            "rtol": self.rtol,
            "protected_keys_present": list(self.protected_keys_present),
            "protected_keys_missing": list(self.protected_keys_missing),
            "n_differences": len(self.differences),
            "differences": [d.to_json() for d in self.differences],
        }

    def summary(self) -> str:
        if self.preserved:
            return (
                f"preserved: {self.compared_leaves:,} protected values identical across "
                f"{len(self.protected_keys_present)} result objects"
            )
        parts = [f"{len(self.differences)} protected value(s) changed"]
        if self.protected_keys_missing:
            parts.append(f"{len(self.protected_keys_missing)} result object(s) missing")
        return "NOT preserved: " + ", ".join(parts)


def _safe(value: Any) -> Any:
    """A JSON-safe, size-bounded rendering of one value for the report."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        if isinstance(value, float) and not math.isfinite(value):
            return repr(value)
        if isinstance(value, str) and len(value) > 300:
            return value[:300] + "…"
        return value
    if isinstance(value, Mapping):
        return f"<{len(value)} keys>"
    if isinstance(value, Sequence):
        return f"<{len(value)} items>"
    return repr(value)[:300]


def _is_volatile(path: str) -> bool:
    """Whether this path, or a parent of it, is enumerated as volatile."""
    if path in VOLATILE_PATHS:
        return True
    return any(path.startswith(v + ".") or path.startswith(v + "[") for v in VOLATILE_PATHS)


def _numbers_match(a: float, b: float, rtol: float) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isinf(a) or math.isinf(b):
        return a == b
    return math.isclose(a, b, rel_tol=rtol, abs_tol=0.0) or a == b


def _walk(
    baseline: Any,
    candidate: Any,
    path: str,
    out: list[Difference],
    rtol: float,
    counter: list[int],
) -> None:
    """Compare two values leaf by leaf, recording each difference."""
    if _is_volatile(path):
        return

    if isinstance(baseline, Mapping):
        if not isinstance(candidate, Mapping):
            out.append(Difference(path, baseline, candidate, "type changed"))
            return
        for key in baseline:
            child = f"{path}.{key}" if path else str(key)
            if key not in candidate:
                if not _is_volatile(child):
                    out.append(Difference(child, baseline[key], None, "key removed"))
                continue
            _walk(baseline[key], candidate[key], child, out, rtol, counter)
        # An added key inside a protected object is a change to that object.
        for key in candidate:
            child = f"{path}.{key}" if path else str(key)
            if key not in baseline and not _is_volatile(child):
                out.append(Difference(child, None, candidate[key], "key added"))
        return

    if isinstance(baseline, (list, tuple)):
        if not isinstance(candidate, (list, tuple)):
            out.append(Difference(path, baseline, candidate, "type changed"))
            return
        if len(baseline) != len(candidate):
            out.append(Difference(
                path, f"<{len(baseline)} items>", f"<{len(candidate)} items>", "length changed",
            ))
            return
        for i, (b, c) in enumerate(zip(baseline, candidate, strict=True)):
            _walk(b, c, f"{path}[{i}]", out, rtol, counter)
        return

    counter[0] += 1
    if isinstance(baseline, bool) or isinstance(candidate, bool):
        # bool is an int in Python; compare it as itself so True != 1 here.
        if baseline is not candidate and baseline != candidate:
            out.append(Difference(path, baseline, candidate, "value changed"))
        return
    if isinstance(baseline, (int, float)) and isinstance(candidate, (int, float)):
        if not _numbers_match(float(baseline), float(candidate), rtol):
            out.append(Difference(path, baseline, candidate, "value changed"))
        return
    if baseline != candidate:
        out.append(Difference(path, baseline, candidate, "value changed"))


def compare(
    baseline: Mapping[str, Any],
    candidate: Mapping[str, Any],
    *,
    keys: Iterable[str] = PROTECTED_KEYS,
    rtol: float = DEFAULT_RTOL,
    baseline_path: str = "baseline",
    candidate_path: str = "candidate",
) -> PreservationReport:
    """Compare every protected leaf of two result objects.

    An object absent from the baseline is not a failure: the baseline may
    predate a module. An object present in the baseline and absent from the
    candidate is a failure, because the extension has removed a result.
    """
    report = PreservationReport(
        baseline_path=baseline_path, candidate_path=candidate_path, rtol=rtol,
    )
    present: list[str] = []
    missing: list[str] = []
    counter = [0]
    for key in keys:
        if key not in baseline:
            continue
        if key not in candidate:
            missing.append(key)
            continue
        present.append(key)
        _walk(baseline[key], candidate[key], key, report.differences, rtol, counter)
    report.compared_leaves = counter[0]
    report.protected_keys_present = tuple(present)
    report.protected_keys_missing = tuple(missing)
    return report


def anonymise(results: Mapping[str, Any], *, local: str, published: str, root: Path | None = None) -> Any:
    """A copy of a result object with the local sample name replaced by its
    published identifier and the checkout's absolute path made relative.

    Frozen baselines are tracked by git and read by anyone who clones the
    repository; they carry the sample's identifier in call ids, trigger ids
    and file paths. Both the frozen copy and the live result are passed
    through this before they are compared, so the comparison sees the same
    strings on both sides and the repository never holds a local name.
    """
    root_str = str((root or Path(__file__).resolve().parents[2]).resolve())

    def fix(value: Any) -> Any:
        if isinstance(value, str):
            out = value
            if root_str and root_str in out:
                out = out.replace(root_str + "/", "").replace(root_str, ".")
            if local and local != published and local in out:
                out = out.replace(local, published)
            return out
        if isinstance(value, Mapping):
            return {fix(k) if isinstance(k, str) else k: fix(v) for k, v in value.items()}
        if isinstance(value, list):
            return [fix(v) for v in value]
        if isinstance(value, tuple):
            return tuple(fix(v) for v in value)
        return value

    return fix(results)


def compare_files(baseline: Path, candidate: Path, *, rtol: float = DEFAULT_RTOL) -> PreservationReport:
    """Compare two `results.json` files on disk."""
    return compare(
        json.loads(baseline.read_text(encoding="utf-8")),
        json.loads(candidate.read_text(encoding="utf-8")),
        rtol=rtol,
        baseline_path=str(baseline),
        candidate_path=str(candidate),
    )


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #


def metric_ids(results: Mapping[str, Any]) -> dict[str, list[str]]:
    """Every existing metric identifier, grouped by the object that owns it.

    Read from a real result object rather than written down, so the manifest
    cannot drift from the application. This is the list section 1 requires to
    be captured before any addition is implemented.
    """
    out: dict[str, list[str]] = {}

    def ids_from(rows: Any, *fields: str) -> list[str]:
        """Stable IDs from a list of records, trying each field in order."""
        found: set[str] = set()
        for row in rows or ():
            if not isinstance(row, Mapping):
                continue
            for f in fields:
                if row.get(f):
                    found.add(str(row[f]))
                    break
        return sorted(found)

    out["panels"] = ids_from(results.get("panels"), "name", "metabolite")
    out["panel_bands"] = sorted(str(k) for k in (results.get("panel_bands") or {}))
    out["metabolite_drivers"] = sorted(str(k) for k in (results.get("metabolite_drivers") or {}))

    sim = results.get("profile_similarity") or {}
    out["disease_profiles"] = ids_from(sim.get("profiles"), "profile_id", "id", "profile", "name")
    out["skin_pattern_indices"] = sorted(str(k) for k in (sim.get("skin_pattern_indices") or {}))
    out["urticaria_pattern_indices"] = sorted(str(k) for k in (sim.get("urticaria_pattern_indices") or {}))

    out["organism_verdicts"] = ids_from(
        (results.get("organism_verdicts") or {}).get("verdicts"), "species", "display",
    )

    bio = results.get("biofilm") or {}
    out["biofilm_cards"] = ids_from(bio.get("cards"), "card", "proxy_id")
    out["biofilm_modules"] = ids_from(bio.get("module_results"), "module_id")
    out["biofilm_findings"] = ids_from(bio.get("ranked_findings"), "id")

    myco = results.get("mycobiome") or {}
    out["mycobiome_scores"] = sorted(
        k for k in ("health_score", "myco_score") if k in myco
    )

    paths = results.get("pathogens") or {}
    out["pathogen_targets"] = sorted(
        str(r.get("target_id")) for r in (paths.get("results") or [])
        if isinstance(r, Mapping) and r.get("target_id")
    )
    out["pathogen_determinants"] = sorted(
        str(d.get("determinant_id")) for d in (paths.get("determinants") or [])
        if isinstance(d, Mapping) and d.get("determinant_id")
    )

    age = results.get("microbiome_age") or {}
    if age.get("model_bundle_id"):
        out["age_models"] = [str(age["model_bundle_id"])]
    elif age:
        out["age_models"] = ["microbiome_age"]

    inv = results.get("organism_inventory") or {}
    # An organism id is still present when the record it named was merged
    # into another (kept as that record's alias) or was rejected by
    # competitive confirmation (kept in the rejected list with its reason);
    # only an id that vanished from all three was silently dropped.
    organism_ids = set(ids_from(inv.get("organisms"), "species"))
    for row in inv.get("organisms") or ():
        if isinstance(row, Mapping):
            organism_ids.update(str(a) for a in (row.get("aliases") or ()))
    organism_ids.update(ids_from(inv.get("rejected"), "species"))
    out["inventory_organisms"] = sorted(organism_ids)

    myco_taxa = (results.get("mycobiome") or {}).get("taxa") or []
    out["mycobiome_taxa"] = ids_from(myco_taxa, "accepted_name", "name")

    strain = results.get("strain_resolution") or {}
    out["strain_records"] = ids_from(strain.get("records") or strain.get("strains"), "species", "id")

    return {k: v for k, v in out.items() if v}


def build_manifest(results: Mapping[str, Any], *, sample: str | None = None) -> dict[str, Any]:
    """The protected-output manifest for one representative result object."""
    ids = metric_ids(results)
    return {
        "schema_version": "openbiota.protected-manifest/1.0",
        "spec_version": SPEC_VERSION,
        "sample": sample or str(results.get("sample") or "unknown"),
        "baseline_version": str(results.get("version") or "unknown"),
        "protected_keys": [k for k in PROTECTED_KEYS if k in results],
        "volatile_paths": sorted(VOLATILE_PATHS),
        "metric_ids": ids,
        "metric_id_count": sum(len(v) for v in ids.values()),
    }


__all__ = [
    "DEFAULT_RTOL",
    "PROTECTED_KEYS",
    "VOLATILE_PATHS",
    "Difference",
    "PreservationReport",
    "build_manifest",
    "compare",
    "compare_files",
    "metric_ids",
]
