"""Marker-level population fingerprints, and comparison between specimens.

This is the layer that answers "which one" for the organisms the data can
support it for. It reads the consensus marker sequences produced by the pinned
MetaPhlAn 4.1.1 marker lane and turns them into typed evidence records.

What a marker consensus is, precisely — because the limits matter as much as
the capability:

* It is a **dominant** consensus. Where two populations of one species are
  present, the consensus reports the majority allele and the minority
  population can be concealed entirely (spec §3.1, acceptance 10). This module
  records that limitation on every call rather than in a footnote.
* It is a **marker** fingerprint, not a genome. Agreement across a few hundred
  conserved marker genes cannot establish whole-genome exact-strain identity
  (acceptance 12), and it cannot assign an unlinked resistance or toxin gene to
  the organism it resolved (acceptance 80).
* Agreement over a **small callable span is not identity**. Zero differences
  across 200 compared positions means something quite different from zero
  across 200,000, so compared positions are reported with every distance and
  a comparison below the floor returns "insufficient resolution" rather than
  a reassuring number (acceptance 13).

What it is good for is real: a per-organism population identity that can be
compared between a donor and a recipient at baseline, which is the foundation
of the FMT strain comparison, and which species abundance cannot approach.
"""

from __future__ import annotations

import bz2
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from .schema import (
    COMPLETED,
    INSUFFICIENT_DEPTH,
    SUPPORTED_DETECTION,
    Coverage,
    Linkage,
    Provenance,
    ResolutionCall,
    Validation,
)

#: Marker names look like `UniRef90_XXXX|1__6|SGB4993`.
_SGB_RE: Final = re.compile(r"\|(SGB\d+)")

#: A marker is only counted when the consensus covered enough of it to be
#: worth comparing. This is the MetaPhlAn workflow's own breadth filter
#: restated as a method parameter, not a clinical detection limit (§11.2).
MIN_MARKER_BREADTH: Final = 80.0

#: Fewer resolved markers than this and the organism's population is reported
#: as detected-but-unresolved. A handful of markers cannot characterise a
#: population.
MIN_MARKERS_FOR_FINGERPRINT: Final = 20

#: Comparisons below this many shared callable bases return insufficient
#: resolution instead of a distance (acceptance 13).
MIN_COMPARED_BASES: Final = 5_000

#: Nucleotide codes that are not a callable base.
_AMBIGUOUS: Final = frozenset("nN-.*?")

DOMINANT_CONSENSUS_LIMIT: Final = (
    "This is the dominant population's consensus across marker genes. A second, less "
    "abundant population of the same species can be concealed by it, so 'one population' "
    "here means 'one detectable dominant population', not proof of a single strain."
)

MARKER_SCOPE_LIMIT: Final = (
    "Marker genes are a conserved subset of the genome. Agreement across them supports a "
    "related population; it does not establish whole-genome identity to a named isolate, "
    "and it cannot attribute an unlinked toxin or resistance gene to this organism."
)


@dataclass(frozen=True, slots=True)
class MarkerFingerprint:
    """One organism's dominant population fingerprint in one specimen."""

    sample_id: str
    sgb: str
    n_markers: int
    callable_bases: int
    median_breadth: float
    median_depth: float
    #: marker name -> consensus sequence, for pairwise comparison.
    sequences: Mapping[str, str]

    @property
    def resolved(self) -> bool:
        return self.n_markers >= MIN_MARKERS_FOR_FINGERPRINT


def _callable_bases(sequence: str) -> int:
    return sum(1 for base in sequence if base not in _AMBIGUOUS)


def load_fingerprints(path: Path, *, sample_id: str) -> tuple[dict[str, MarkerFingerprint], str]:
    """Read a `sample2markers` output into per-SGB fingerprints.

    Returns the fingerprints and the database release the markers came from,
    so a later comparison can refuse to mix database versions.
    """
    with bz2.open(path, "rt") as handle:
        payload = json.load(handle)
    database = str(payload.get("database_name") or "")

    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for marker in payload.get("consensus_markers") or []:
        name = str(marker.get("marker") or "")
        breadth = float(marker.get("breath") or 0.0)
        if breadth < MIN_MARKER_BREADTH:
            continue
        found = _SGB_RE.search(name)
        if not found:
            continue
        grouped.setdefault(found.group(1), []).append(marker)

    out: dict[str, MarkerFingerprint] = {}
    for sgb, markers in grouped.items():
        breadths = sorted(float(m.get("breath") or 0.0) for m in markers)
        depths = sorted(float(m.get("avg_depth") or 0.0) for m in markers)
        sequences = {
            str(m.get("marker")): str(m.get("sequence") or "") for m in markers
        }
        out[sgb] = MarkerFingerprint(
            sample_id=sample_id,
            sgb=sgb,
            n_markers=len(markers),
            callable_bases=sum(_callable_bases(s) for s in sequences.values()),
            median_breadth=breadths[len(breadths) // 2],
            median_depth=depths[len(depths) // 2],
            sequences=sequences,
        )
    return out, database


def calls_for_sample(
    fingerprints: Mapping[str, MarkerFingerprint],
    *,
    sample_id: str,
    database: str,
    species_of: Mapping[str, str] | None = None,
) -> list[ResolutionCall]:
    """Typed population-fingerprint calls, one per organism."""
    names = species_of or {}
    out: list[ResolutionCall] = []
    for sgb, fp in sorted(fingerprints.items()):
        species = names.get(sgb, "")
        resolved = fp.resolved
        call = ResolutionCall(
            sample_id=sample_id,
            target_id=f"strain.population.{sgb}",
            identity_kind="population_fingerprint",
            assay_status=COMPLETED,
            call_id=f"{sample_id}:strain.population.{sgb}",
            # A fingerprint is a supported detection of a population, or an
            # organism whose population could not be characterised. It is
            # never a named strain.
            analytical_call=SUPPORTED_DETECTION if resolved else INSUFFICIENT_DEPTH,
            reason_codes=(
                ("dominant_consensus_only",)
                if resolved
                else ("dominant_consensus_only", "too_few_resolved_markers")
            ),
            taxonomic_assertions=(
                {
                    "rank": "sgb",
                    "identifier": sgb,
                    "species_label": species,
                    "database": database,
                    # Spec §3.1 and acceptance 6: an SGB is not a named strain.
                    "is_named_strain": False,
                },
            ),
            population_components=(
                {
                    "component_id": f"{sample_id}:{sgb}:dominant",
                    "method": "metaphlan_marker_consensus",
                    "role": "dominant",
                    "n_markers": fp.n_markers,
                    "median_breadth_percent": round(fp.median_breadth, 2),
                    "median_depth": round(fp.median_depth, 3),
                    "minority_populations_resolved": False,
                },
            ),
            coverage=Coverage(
                callable_bases=fp.callable_bases,
                target_breadth=round(fp.median_breadth / 100.0, 4),
                median_depth=round(fp.median_depth, 3),
            ),
            # Resolving a population does not link anything else to it.
            linkage=Linkage(state="unassigned"),
            validation=Validation(
                analytical_status="method_defaults_not_locally_calibrated",
                reference_function_status="not_applicable",
                phenotype_measured_in_sample=False,
                clinical_predictive_status="not_established",
            ),
            provenance=Provenance(
                database_release=database,
                tool_version="metaphlan 4.1.1 sample2markers",
                parameters={
                    "min_marker_breadth": MIN_MARKER_BREADTH,
                    "min_markers_for_fingerprint": MIN_MARKERS_FOR_FINGERPRINT,
                },
                threshold_version="strain.marker/1",
            ),
            plain=(
                (
                    f"{species or sgb}: dominant population characterised across "
                    f"{fp.n_markers} marker genes ({fp.callable_bases:,} callable bases, "
                    f"median breadth {fp.median_breadth:.0f}%). "
                )
                if resolved
                else (
                    f"{species or sgb}: detected, but only {fp.n_markers} marker gene(s) "
                    "were resolved \u2014 too few to characterise the population. "
                )
            )
            + DOMINANT_CONSENSUS_LIMIT,
        )
        out.append(call)
    return out


# --------------------------------------------------------------------------- #
# Pairwise comparison (the FMT baseline capability, spec §10.1)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Comparison:
    """Relatedness of two specimens' populations of one organism."""

    sgb: str
    left_sample: str
    right_sample: str
    shared_markers: int
    compared_bases: int
    differing_bases: int
    #: None when the comparison did not clear the callable-span floor.
    difference_per_kb: float | None
    status: str
    plain: str

    def to_json(self) -> dict[str, Any]:
        return {
            "sgb": self.sgb,
            "left_sample": self.left_sample,
            "right_sample": self.right_sample,
            "shared_markers": self.shared_markers,
            "compared_bases": self.compared_bases,
            "differing_bases": self.differing_bases,
            "difference_per_kb": self.difference_per_kb,
            "status": self.status,
            "plain": self.plain,
            "limits": [DOMINANT_CONSENSUS_LIMIT, MARKER_SCOPE_LIMIT],
            # Relatedness is not a forecast (spec §10.1).
            "predicts_establishment": False,
            "disease_transmission_probability": None,
        }


def compare(
    left: MarkerFingerprint, right: MarkerFingerprint, *, sgb: str
) -> Comparison:
    """Compare two dominant consensuses over their shared callable positions.

    Only positions callable in *both* specimens are compared, and the count is
    reported: a distance without its denominator is uninterpretable.
    """
    shared = set(left.sequences) & set(right.sequences)
    compared = 0
    differing = 0
    for marker in shared:
        a, b = left.sequences[marker], right.sequences[marker]
        for base_a, base_b in zip(a, b, strict=False):
            if base_a in _AMBIGUOUS or base_b in _AMBIGUOUS:
                continue
            compared += 1
            if base_a != base_b:
                differing += 1

    if compared < MIN_COMPARED_BASES:
        return Comparison(
            sgb=sgb,
            left_sample=left.sample_id,
            right_sample=right.sample_id,
            shared_markers=len(shared),
            compared_bases=compared,
            differing_bases=differing,
            difference_per_kb=None,
            status="insufficient_callable_overlap",
            plain=(
                f"Only {compared:,} positions were callable in both specimens, below the "
                f"{MIN_COMPARED_BASES:,} needed to say anything about relatedness. No "
                "distance is reported: agreement over a short span is not evidence of "
                "sharing."
            ),
        )

    per_kb = 1000.0 * differing / compared
    if differing == 0:
        status = "indistinguishable_over_compared_sites"
        plain = (
            f"No differences across {compared:,} positions callable in both specimens "
            f"({len(shared)} shared markers). These populations are indistinguishable at "
            "marker resolution, which is consistent with a shared population but is not "
            "proof of whole-genome identity."
        )
    elif per_kb < 1.0:
        status = "closely_related"
        plain = (
            f"{differing} difference(s) across {compared:,} compared positions "
            f"({per_kb:.2f} per kb) \u2014 closely related populations."
        )
    else:
        status = "distinct_populations"
        plain = (
            f"{differing:,} differences across {compared:,} compared positions "
            f"({per_kb:.2f} per kb) \u2014 distinct populations of this organism."
        )
    return Comparison(
        sgb=sgb,
        left_sample=left.sample_id,
        right_sample=right.sample_id,
        shared_markers=len(shared),
        compared_bases=compared,
        differing_bases=differing,
        difference_per_kb=round(per_kb, 4),
        status=status,
        plain=plain,
    )


def compare_all(
    left: Mapping[str, MarkerFingerprint],
    right: Mapping[str, MarkerFingerprint],
) -> list[Comparison]:
    """Every organism both specimens resolved, most-compared first."""
    out = [
        compare(left[sgb], right[sgb], sgb=sgb)
        for sgb in sorted(set(left) & set(right))
        if left[sgb].resolved and right[sgb].resolved
    ]
    return sorted(out, key=lambda c: -c.compared_bases)


def self_test() -> int:
    """Guard the comparison limits. Returns a failure count."""
    failures = 0

    def fp(sample: str, seqs: dict[str, str]) -> MarkerFingerprint:
        return MarkerFingerprint(
            sample_id=sample, sgb="SGB1", n_markers=len(seqs),
            callable_bases=sum(_callable_bases(s) for s in seqs.values()),
            median_breadth=95.0, median_depth=10.0, sequences=seqs,
        )

    # Acceptance 13: a short span cannot establish sharing, even at zero SNPs.
    short = {"m1|1__2|SGB1": "ACGT" * 100}
    result = compare(fp("A", short), fp("B", short), sgb="SGB1")
    if result.status != "insufficient_callable_overlap":
        failures += 1
    if result.difference_per_kb is not None:
        failures += 1
    if "not evidence of" not in result.plain:
        failures += 1

    # A long identical span is reported as indistinguishable, not identical.
    long_seq = {"m1|1__2|SGB1": "ACGT" * 3000}
    result = compare(fp("A", long_seq), fp("B", long_seq), sgb="SGB1")
    if result.status != "indistinguishable_over_compared_sites":
        failures += 1
    if "not proof of whole-genome identity" not in result.plain:
        failures += 1
    if result.compared_bases != 12000:
        failures += 1

    # Ambiguous positions are excluded from the denominator.
    masked = {"m1|1__2|SGB1": ("ACGT" * 2999) + "NNNN"}
    result = compare(fp("A", long_seq), fp("B", masked), sgb="SGB1")
    if result.compared_bases != 11996:
        failures += 1

    # Differences are counted and normalised.
    changed = {"m1|1__2|SGB1": ("TCGT" + "ACGT" * 2999)}
    result = compare(fp("A", long_seq), fp("B", changed), sgb="SGB1")
    if result.differing_bases != 1 or result.status != "closely_related":
        failures += 1

    # Every comparison must carry its limits and no forecast.
    blob = result.to_json()
    if blob["predicts_establishment"] is not False:
        failures += 1
    if blob["disease_transmission_probability"] is not None:
        failures += 1
    if len(blob["limits"]) != 2:
        failures += 1

    # A thin fingerprint is detected-but-unresolved, never a named strain.
    thin = fp("A", {f"m{i}|1__2|SGB1": "ACGT" * 50 for i in range(3)})
    calls = calls_for_sample({"SGB1": thin}, sample_id="A", database="db")
    if calls[0].analytical_call != INSUFFICIENT_DEPTH:
        failures += 1
    if calls[0].taxonomic_assertions[0]["is_named_strain"] is not False:
        failures += 1
    if calls[0].linkage.state != "unassigned":
        failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
