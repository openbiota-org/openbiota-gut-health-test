"""Batch runs: every sample in a directory, plus a cross-sample comparison.

Each sample gets its own directory and its own named report, exactly as a
single run would. What batch mode adds is the comparison: once several samples
have been through the identical pipeline, the interesting question stops being
"where does this sample sit against a public cohort" and becomes "how do these
samples differ from each other".

The comparison is deliberately thin. It reports what was measured side by
side and does not attempt statistics across a handful of samples.
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from openbiota.logging_util import Reporter, human_duration


@dataclass(slots=True)
class SampleOutcome:
    """One sample's headline numbers, for the comparison table."""

    sample: str
    ok: bool
    out_dir: Path
    elapsed_s: float
    read_pairs: int | None = None
    #: panel name -> copies per 100 genomes
    panels: dict[str, float | None] = field(default_factory=dict)
    #: panel name -> percentile against the reference cohort
    percentiles: dict[str, float | None] = field(default_factory=dict)
    n_species: int | None = None
    shannon: float | None = None
    gmwi2: float | None = None
    gmwi2_band: str = ""
    #: profile name -> combined percentile (None when abstained)
    profiles: dict[str, float | None] = field(default_factory=dict)
    profile_confidence: dict[str, str] = field(default_factory=dict)
    host_fraction: float | None = None
    error: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "sample": self.sample,
            "ok": self.ok,
            "output_directory": str(self.out_dir),
            "elapsed_s": round(self.elapsed_s, 1),
            "read_pairs": self.read_pairs,
            "panels_copies_per_100_genomes": self.panels,
            "panel_percentiles": self.percentiles,
            "n_species": self.n_species,
            "shannon_diversity": self.shannon,
            "gmwi2": self.gmwi2,
            "gmwi2_band": self.gmwi2_band,
            "profile_percentiles": self.profiles,
            "profile_confidence": self.profile_confidence,
            "host_fraction": self.host_fraction,
            "error": self.error,
        }


def outcome_from_results(
    sample: str, out_dir: Path, elapsed_s: float, results: dict[str, Any]
) -> SampleOutcome:
    """Extract the comparison fields from a finished results.json."""
    outcome = SampleOutcome(sample=sample, ok=True, out_dir=out_dir, elapsed_s=elapsed_s)

    qc = results.get("input") or {}
    outcome.read_pairs = qc.get("read_pairs")

    for panel in results.get("panels") or []:
        if isinstance(panel, dict) and panel.get("name"):
            outcome.panels[panel["name"]] = panel.get("copies_per_100_genomes")

    # reference_comparison.panels is keyed by panel name rather than being a
    # list of rows, so it has to be iterated as a mapping.
    comparison = (results.get("reference_comparison") or {}).get("panels") or {}
    if isinstance(comparison, dict):
        for name, row in comparison.items():
            if isinstance(row, dict):
                outcome.percentiles[name] = row.get("percentile")
    else:
        for row in comparison:
            if isinstance(row, dict) and row.get("panel"):
                outcome.percentiles[row["panel"]] = row.get("percentile")

    similarity = results.get("profile_similarity") or {}
    engine = similarity.get("taxonomic_engine") or {}
    outcome.n_species = engine.get("n_species_detected")
    host = engine.get("host_filter") or {}
    outcome.host_fraction = host.get("host_fraction")

    anchor = similarity.get("dysbiosis_anchor") or {}
    outcome.gmwi2 = anchor.get("score")
    outcome.gmwi2_band = anchor.get("band", "")

    for name, payload in (similarity.get("profiles") or {}).items():
        abstained = (payload.get("abstention") or {}).get("abstained")
        combined = (payload.get("combined") or {}).get("percentile")
        outcome.profiles[name] = None if abstained else combined
        outcome.profile_confidence[name] = (
            "abstained" if abstained else (payload.get("confidence") or {}).get("grade", "")
        )

    # Shannon diversity is measured once and appears in whichever profiles
    # declare it as an ecological feature; take it from the first that has it.
    for payload in (similarity.get("profiles") or {}).values():
        features = (
            (payload.get("modules") or {}).get("ecological") or {}
        ).get("features") or []
        for feature in features:
            if isinstance(feature, dict) and feature.get("name") == "shannon_diversity":
                outcome.shannon = feature.get("raw_value")
                break
        if outcome.shannon is not None:
            break

    return outcome


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #


def comparison_json(outcomes: Sequence[SampleOutcome]) -> dict[str, Any]:
    ok = [o for o in outcomes if o.ok]
    panels = sorted({k for o in ok for k in o.panels})
    profiles = sorted({k for o in ok for k in o.profiles})
    return {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "n_samples": len(outcomes),
        "n_succeeded": len(ok),
        "samples": [o.to_json() for o in outcomes],
        "panels_compared": panels,
        "profiles_compared": profiles,
        "note": (
            "Every sample went through the identical pipeline, so the columns are directly "
            "comparable. With a handful of samples no statistics are computed across them; "
            "this is a side-by-side view, not a study."
        ),
    }


def _fmt(value: float | None, digits: int = 1) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


#: Column width is 13 characters, so the full band names have to be shortened
#: deliberately rather than truncated — "strongly dys" reads as an error.
_BAND_SHORT: dict[str, str] = {
    "strongly dysbiotic": "dysbiotic!!",
    "mildly dysbiotic": "dysbiotic",
    "indeterminate": "indeterm.",
    "not dysbiotic": "healthy",
    "not computed": "—",
}


def _short_band(band: str) -> str:
    if not band:
        return "—"
    return _BAND_SHORT.get(band, band[:11])


def print_comparison(outcomes: Sequence[SampleOutcome]) -> None:
    """Side-by-side table across samples, printed to stdout."""
    ok = [o for o in outcomes if o.ok]
    failed = [o for o in outcomes if not o.ok]
    if not ok:
        print("  no samples completed successfully")
        return

    names = [o.sample for o in ok]
    width = max(22, max(len(n) for n in names) + 2)

    def header(title: str) -> None:
        print()
        print(f"  {title}")
        print(f"  {'':<{width}}" + "".join(f"{n[:11]:>13}" for n in names))
        print("  " + "-" * (width + 13 * len(names)))

    def line(label: str, values: Sequence[str]) -> None:
        print(f"  {label:<{width}}" + "".join(f"{v:>13}" for v in values))

    header("SEQUENCING")
    line("read pairs (millions)", [_fmt((o.read_pairs or 0) / 1e6, 1) for o in ok])
    line("species detected", [str(o.n_species or "—") for o in ok])
    line(
        "human DNA",
        ["—" if o.host_fraction is None else f"{o.host_fraction:.2%}" for o in ok],
    )

    header("GENERAL GUT HEALTH (GMWI2)")
    line("index", [_fmt(o.gmwi2, 2) for o in ok])
    line("band", [_short_band(o.gmwi2_band) for o in ok])
    line("diversity (Shannon)", [_fmt(o.shannon, 2) for o in ok])

    panels = sorted({k for o in ok for k in o.panels})
    if panels:
        header("METABOLITE PATHWAYS — copies per 100 genomes")
        for panel in panels:
            line(panel, [_fmt(o.panels.get(panel), 1) for o in ok])
        header("METABOLITE PATHWAYS — percentile vs reference cohort")
        for panel in panels:
            line(panel, [_fmt(o.percentiles.get(panel), 0) for o in ok])

    profiles = sorted({k for o in ok for k in o.profiles})
    if profiles:
        header("DISEASE-PATTERN SIMILARITY — percentile (resemblance, not diagnosis)")
        for profile in profiles:
            values = []
            for outcome in ok:
                if outcome.profile_confidence.get(profile) == "abstained":
                    values.append("abstained")
                else:
                    values.append(_fmt(outcome.profiles.get(profile), 0))
            line(profile, values)

    print()
    for outcome in ok:
        print(
            f"  {outcome.sample:<{width}}{outcome.out_dir}/{outcome.sample}_report.pdf"
        )
    if failed:
        print()
        for outcome in failed:
            print(f"  FAILED  {outcome.sample}: {outcome.error}")
    print()


def write_comparison(path: Path, outcomes: Sequence[SampleOutcome]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(comparison_json(outcomes), indent=2) + "\n", encoding="utf-8")
    return path


def report_progress(
    reporter: Reporter, index: int, total: int, sample: str, elapsed: float
) -> None:
    reporter.ok(
        f"[{index}/{total}] {sample} done in {human_duration(elapsed)}"
    )
