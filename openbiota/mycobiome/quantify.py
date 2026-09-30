"""Fungal signal and within-fungi composition (spec §7.1-7.3).

Two denominators, never mixed:

* **Fungal sequence signal** - ``F / N``: confidently fungal fragments over
  *all* eligible QC non-host fragments, counted before any bacterial
  filtering. Also given per million; the two agree exactly before display
  rounding. The sampling interval is a Wilson interval on the fraction and
  covers counting noise only, not extraction or classification bias.
* **Within-fungi composition** - one estimator, coverage-normalised:
  ``c_i = fragments_i × fragment_bases / backbone_bases_i``, ``q_i = c_i /
  Σc``. Species not quantifiable (genus-level or provisional support) are
  reported as fragment counts against the fungal-fragment denominator, not
  inserted into the composition bar.

A quantity that cannot be measured is null with a state, never zero.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Final

from openbiota.mycobiome.genome_lane import Ledger, TaxonSupport, poisson_interval

#: Starting limit of quantification for the total-sample fungal fraction:
#: below this many supported fragments the count is reported but the
#: percentage is labelled below the quantification limit (Poisson CV > 14%).
LOQ_FRAGMENTS: Final = 50
ESTIMATOR_ID: Final = "coverage_normalised_fragments_v1"


@dataclass(slots=True)
class Measurement:
    denominator_stage: str
    eligible_fragments: int | None
    fungal_supported_fragments: int | None
    cross_kingdom_ambiguous_fragments: int | None
    unassigned_fragments: int | None
    fungal_fragment_fraction: float | None
    fungal_fragments_per_million: float | None
    sampling_interval: tuple[float, float] | None
    quantification_status: str
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "denominator_stage": self.denominator_stage,
            "eligible_fragments": self.eligible_fragments,
            "fungal_supported_fragments": self.fungal_supported_fragments,
            "cross_kingdom_ambiguous_fragments": self.cross_kingdom_ambiguous_fragments,
            "unassigned_fragments": self.unassigned_fragments,
            "fungal_fragment_fraction": self.fungal_fragment_fraction,
            "fungal_fragments_per_million": self.fungal_fragments_per_million,
            "sampling_interval": None if self.sampling_interval is None else
                {"low": self.sampling_interval[0], "high": self.sampling_interval[1], "kind": "wilson_95_counting_only"},
            "quantification_status": self.quantification_status,
            "label": "fungal sequence signal (% of analyzed nonhost fragments)",
            "note": self.note,
        }


@dataclass(slots=True)
class CompositionRow:
    key: str
    name: str
    rank: str
    fragments: int
    backbone_bases: int
    coverage_estimate: float
    share: float
    detection_state: str

    def to_json(self) -> dict[str, Any]:
        return {"key": self.key, "name": self.name, "rank": self.rank, "fragments": self.fragments,
                "backbone_bases": self.backbone_bases, "coverage_estimate": round(self.coverage_estimate, 6),
                "share_of_quantified_fungi": round(self.share, 6), "detection_state": self.detection_state}


@dataclass(slots=True)
class Composition:
    estimator_id: str
    quantity_type: str
    denominator_id: str
    rows: list[CompositionRow]
    unresolved_fragment_fraction: float | None
    unresolved_fragment_denominator_id: str
    unresolved_fragments: int
    quantifiable: bool
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "estimator_id": self.estimator_id, "quantity_type": self.quantity_type, "denominator_id": self.denominator_id,
            "unresolved_composition_fraction": None,
            "unresolved_fragment_fraction": self.unresolved_fragment_fraction,
            "unresolved_fragment_denominator_id": self.unresolved_fragment_denominator_id,
            "unresolved_fragments": self.unresolved_fragments, "quantifiable": self.quantifiable,
            "composition": [r.to_json() for r in self.rows], "note": self.note,
        }


@dataclass(slots=True)
class Diversity:
    supported_species: int
    provisional_taxa: int
    unresolved_complexes: int
    note: str = ("Observed richness of supported species only. Not a directional health measure; zero detections "
                 "is not a reassuring diversity of zero.")

    def to_json(self) -> dict[str, Any]:
        return {"supported_species": self.supported_species, "provisional_taxa": self.provisional_taxa,
                "unresolved_complexes": self.unresolved_complexes, "shannon": None,
                "shannon_note": "not reported: rarefaction stability not demonstrated at this fungal depth", "note": self.note}


def measurement(ledger: Ledger | None, *, eligible: int | None, denominator_valid: bool, lane_status: str) -> Measurement:
    if ledger is None or lane_status != "resolved":
        return Measurement("qc_nonhost_before_bacterial_filter", eligible, None, None, None, None, None, None,
                           "not_assessed", note="whole-genome lane did not complete")
    if not denominator_valid or not eligible:
        return Measurement("qc_nonhost_before_bacterial_filter", eligible, ledger.fungal_confident,
                           ledger.cross_kingdom_ambiguous, None, None, None, None,
                           "not_quantifiable_from_available_input",
                           note="no valid original denominator: fungal detections stand, the total-sample fraction does not")
    f, n = ledger.fungal_confident, eligible
    frac = f / n
    ppm = 1_000_000.0 * f / n
    lo, hi = poisson_interval(f, n)
    status = "quantifiable" if f >= LOQ_FRAGMENTS else ("detected_below_quantification_limit" if f > 0 else "quantifiable")
    return Measurement(
        "qc_nonhost_before_bacterial_filter", n, f, ledger.cross_kingdom_ambiguous, max(0, n - ledger.aligned_fragments),
        frac, ppm, (lo, hi), status,
        note=(f"{f} fungal fragments among {n:,} eligible; percentage precision limited below {LOQ_FRAGMENTS} fragments"
              if status == "detected_below_quantification_limit" else ""),
    )


def composition(taxa: list[TaxonSupport], *, fragment_bases: float, ledger: Ledger | None) -> Composition:
    """Coverage-normalised shares among supported species with a backbone."""
    quant = [t for t in taxa if t.detection_state == "supported" and t.rank == "species" and t.genome_bases > 0]
    covs = {t.key: (t.fragments * fragment_bases) / t.genome_bases for t in quant}
    total = sum(covs.values())
    rows = [
        CompositionRow(key=t.key, name=t.name, rank=t.rank, fragments=t.fragments, backbone_bases=t.genome_bases,
                       coverage_estimate=covs[t.key], share=(covs[t.key] / total) if total > 0 else 0.0,
                       detection_state=t.detection_state)
        for t in sorted(quant, key=lambda t: -covs[t.key])
    ]
    unresolved = sum(t.fragments for t in taxa if t not in quant)
    f_total = ledger.fungal_confident if ledger else None
    return Composition(
        estimator_id=ESTIMATOR_ID, quantity_type="genome_coverage_normalised_share",
        denominator_id="sum_of_coverage_estimates_over_supported_species", rows=rows,
        unresolved_fragment_fraction=(unresolved / f_total) if f_total else None,
        unresolved_fragment_denominator_id="fungal_supported_fragments", unresolved_fragments=unresolved,
        quantifiable=bool(rows),
        note=("Shares are genome-coverage estimates among supported species; they are not cell counts and they do not "
              "include genus-level or provisional support, which is reported as fragments against the fungal total."),
    )


def diversity(taxa: list[TaxonSupport]) -> Diversity:
    return Diversity(
        supported_species=sum(1 for t in taxa if t.detection_state == "supported" and t.rank == "species"),
        provisional_taxa=sum(1 for t in taxa if t.detection_state == "provisional"),
        unresolved_complexes=sum(1 for t in taxa if t.detection_state == "ambiguous_complex"),
    )


__all__ = ["Composition", "CompositionRow", "Diversity", "ESTIMATOR_ID", "LOQ_FRAGMENTS", "Measurement", "composition",
           "diversity", "measurement"]
