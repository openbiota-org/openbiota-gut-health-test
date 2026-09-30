"""Chronic urticaria (chronic hives) directional pattern indices — spec v04.1.

Two registered tasks: chronic spontaneous urticaria (chronic hives) from Zhu
2024, and symptomatic dermographism (the hives raised by stroking the skin)
from Liu 2021. The arithmetic, the missingness policy and the reference
fitting all live in :mod:`openbiota.dirindex`, which the inflammatory-skin panels
share; this module holds only what is specific to hives — the two feature
manifests, the invariants that pin them to their source tables, and the
top-level entry point the pipeline calls.

Everything the shared kernel guarantees applies here: equal weights, source
coefficients recorded as provenance only, unmeasurable slots lowering
coverage instead of scoring zero, and coverage-and-censoring bounds that are
not confidence intervals.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Final

from openbiota.dirindex import (
    BOUNDS_TYPE,
    EPS,
    MAD_TO_SIGMA,
    MEASUREMENT_STATES,
    MIN_CONTROLS,
    MIN_REFERENCE_PREVALENCE,
    SCALE_FLOOR,
    SCORE_SEMANTICS,
    SCORE_SEMANTICS_TRANSPORTED,
    TEMPERATURE,
    UNAVAILABLE_STATES,
    Z_CLIP,
    DirIndexError,
    FeatureRow,
    FrozenReference,
    Panel,
    PanelFeature,
    PanelResult,
    comparable,
    evaluate,
    fit_reference,
    fit_reference_detected,
    fraction,
    freeze_reference,
    score_panel,
    self_test,
    summarize,
    support,
    weighted_midrank,
    z_value,
)

#: Historical alias: this module raised its own error type before the kernel
#: was shared, and callers still catch that name.
CuIndexError = DirIndexError
#: Historical alias for the shared result type.
CuResult = PanelResult

__all__ = [
    "BOUNDS_TYPE", "CSU_PANEL", "CSU_STRICT_PANEL", "EPS", "MAD_TO_SIGMA",
    "MEASUREMENT_STATES", "MIN_CONTROLS", "MIN_REFERENCE_PREVALENCE", "PANELS",
    "SCALE_FLOOR", "SCORE_SEMANTICS", "SCORE_SEMANTICS_TRANSPORTED", "SD_PANEL",
    "TEMPERATURE", "UNAVAILABLE_STATES", "Z_CLIP", "CuIndexError", "CuResult",
    "DirIndexError", "FeatureRow", "FrozenReference", "Panel", "PanelFeature",
    "PanelResult", "assert_panel_invariants", "comparable", "evaluate",
    "fit_reference", "fit_reference_detected", "fraction", "freeze_reference",
    "score_panel", "self_test", "summarize", "support", "urticaria_indices",
    "weighted_midrank", "z_value",
]


#: Zhu 2024, source workbook sheet ``Fig 1b``: fifteen species-level
#: directions from 26 adults with chronic spontaneous urticaria (chronic
#: hives) and 26 healthy controls. Ten decreases, five increases. Coefficients
#: and q-values are provenance; all computational weights are equal.
CSU_PANEL: Final = Panel(
    panel_id="CSU_STOOL_WGS_V1",
    revision="4.1.1",
    label="Chronic spontaneous urticaria (chronic hives) — Zhu 2024 fifteen-species panel",
    source_id="ZHU_2024",
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    score_semantics=SCORE_SEMANTICS,
    discovery_sample_count=52,
    features=(
        PanelFeature("Alistipes_onderdonkii", -1, coef=-2.469669346, q_value=0.015384911),
        PanelFeature("Alistipes_putredinis", -1, coef=-2.459051270, q_value=0.027690473),
        PanelFeature("Alistipes_shahii", -1, coef=-2.039284071, q_value=0.027690473),
        PanelFeature("Bacteroides_cellulosilyticus", -1, coef=-1.642682855, q_value=0.100549642),
        PanelFeature("Bacteroides_intestinalis", -1, coef=-1.888324131, q_value=0.042168812),
        PanelFeature("Coprococcus_catus", -1, coef=-2.167381929, q_value=0.002208388),
        PanelFeature("Odoribacter_splanchnicus", -1, coef=-1.955020950, q_value=0.098787634),
        PanelFeature("Roseburia_hominis", -1, coef=-1.359086378, q_value=0.128357101),
        PanelFeature(
            "Ruminococcus_obeum",
            -1,
            measured_as="Blautia_obeum",
            mapping_reason=(
                "Ruminococcus obeum was transferred to Blautia obeum (Liu 2008, Int J Syst Evol "
                "Microbiol 58:1896). Same organism and same biological scope, so the source name "
                "identifies the slot and the current name is what the profiler reports."
            ),
            coef=-0.831660425,
            q_value=0.022195148,
        ),
        PanelFeature("Veillonella_parvula", 1, coef=1.350313455, q_value=0.047789302),
        PanelFeature("Bacteroides_stercoris", 1, coef=1.827477087, q_value=0.119496077),
        PanelFeature("Bilophila_wadsworthia", -1, coef=-1.083777733, q_value=0.119496077),
        PanelFeature("Escherichia_coli", 1, coef=1.389339140, q_value=0.107643794),
        PanelFeature("Klebsiella_pneumoniae", 1, coef=1.713724620, q_value=0.124562782),
        PanelFeature("Ruminococcus_gnavus", 1, coef=3.057857054, q_value=0.002208388),
    ),
)

#: The eight of fifteen that pass Benjamini-Hochberg correction. A sensitivity
#: analysis on the same measurements, not a second independent test.
CSU_STRICT_PANEL: Final = CSU_PANEL.subset(
    "CSU_STOOL_WGS_V1_STRICT8",
    "Chronic hives — strict subset (BH q < 0.05), sensitivity analysis",
    max_q=0.05,
)

#: Liu 2021, symptomatic dermographism (hives raised by stroking the skin):
#: two depletion-oriented candidates from the healthy-control core community
#: of a 16S cohort, carried across to shotgun measurement.
SD_PANEL: Final = Panel(
    panel_id="SD_STOOL_TRANSLATED_V1",
    revision="4.1.1",
    label="Symptomatic dermographism (hives from stroking the skin) — Liu 2021, translated",
    source_id="LIU_SD_2021",
    source_assay="16S_amplicon",
    scoring_assay="shotgun_metagenomics",
    score_semantics=SCORE_SEMANTICS_TRANSPORTED,
    assay_transport="amplicon_taxon_to_shotgun",
    features=(
        PanelFeature(
            "Subdoligranulum",
            -1,
            rank="genus",  # summed over the catalogue's children at fit time
            mapping_reason=(
                "Summed over this catalogue's Subdoligranulum species: a non-overlapping set of "
                "children with no parent total included, and no overlap with Ruminococcus bromii."
            ),
        ),
        PanelFeature("Ruminococcus_bromii", -1),
    ),
)

PANELS: Final = {p.panel_id: p for p in (CSU_PANEL, CSU_STRICT_PANEL, SD_PANEL)}


def assert_panel_invariants() -> None:
    """The manifest facts the specification pins down, checked at import.

    Cheap, and it means a typo in a direction or a q-value cannot reach a
    report: fifteen unique species, five up and ten down, exactly eight
    surviving correction.
    """
    ids = CSU_PANEL.ids
    if len(ids) != 15 or len(set(ids)) != 15:
        raise CuIndexError(f"CSU panel must hold 15 unique species, got {len(set(ids))}")
    up = sum(1 for f in CSU_PANEL.features if f.direction == 1)
    down = sum(1 for f in CSU_PANEL.features if f.direction == -1)
    if (up, down) != (5, 10):
        raise CuIndexError(f"CSU panel must be 5 increases and 10 decreases, got {up} and {down}")
    strict = set(CSU_STRICT_PANEL.ids)
    expected = {
        "Alistipes_onderdonkii", "Alistipes_putredinis", "Alistipes_shahii",
        "Bacteroides_intestinalis", "Coprococcus_catus", "Ruminococcus_obeum",
        "Veillonella_parvula", "Ruminococcus_gnavus",
    }
    if strict != expected:
        raise CuIndexError(f"strict subset must be exactly the eight q<0.05 species, got {sorted(strict)}")
    if len(SD_PANEL.features) != 2:
        raise CuIndexError("SD panel holds exactly two curated candidates")
    keys = {f.key for f in SD_PANEL.features}
    if len(keys) != 2:
        raise CuIndexError("SD panel slots must not overlap")


assert_panel_invariants()


def urticaria_indices(
    species_percent: Mapping[str, float],
    cohort: Any | None,
    *,
    bundle_id: str | None = None,
    reporting_limit: float | None = None,
) -> dict[str, CuResult]:
    """Score both urticaria (hives) panels for one sample.

    ``species_percent`` is this sample's species table on the percent scale
    the profiler reports, and ``cohort`` is the matched reference cohort — the
    same object the rest of the report ranks against, so the two summaries
    rest on one set of reference participants rather than two.

    No reporting limit is passed by the pipeline. This profiler's detection
    limit depends on sequencing depth and on the organism, and no validated
    per-feature limit exists for these fifteen species, so a nondetection is
    left unquantifiable: it keeps its bounds and drops its point rather than
    being scored as a zero that was never measured.
    """
    if cohort is None:
        return {
            CSU_PANEL.panel_id: score_panel(CSU_PANEL, {}, None),
            SD_PANEL.panel_id: score_panel(SD_PANEL, {}, None),
        }

    abundances = {k: min(max(float(v) / 100.0, 0.0), 1.0) for k, v in species_percent.items()}
    label = bundle_id or f"{getattr(cohort, 'snapshot', 'reference')}:n={getattr(cohort, 'n_samples', 0)}"
    out: dict[str, CuResult] = {}
    csu: CuResult | None = None
    for panel in (CSU_PANEL, CSU_STRICT_PANEL, SD_PANEL):
        reference = freeze_reference(
            panel,
            catalogue=cohort.taxa,
            matrix=cohort.abundance,
            bundle_id=label,
        )
        result = score_panel(panel, abundances, reference, reporting_limit=reporting_limit)
        if panel is CSU_PANEL:
            csu = result
        elif panel is CSU_STRICT_PANEL and csu is not None:
            csu.sensitivity = result
            continue
        out[panel.panel_id] = result
    return out

