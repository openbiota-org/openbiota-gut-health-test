"""Capability manifest (spec §4.1, §10, §15.1).

Every optional extension declares an honest status and, when it is blocked,
the exact artifacts that would unblock it. The spec is explicit that an
extension is not "implemented" because its paper is linked, and that a release
must not be marked complete with placeholder adapter classes — so the
adapters below report `artifact_blocked` rather than pretending.

An unavailable optional capability never fails a complete core comparison.
"""

from __future__ import annotations

import importlib.util
import shutil
from typing import Any


def _module_present(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return False


def manifest() -> list[dict[str, Any]]:
    """The capability list for this environment, computed at run time."""
    have_micom = _module_present("micom")
    have_samestr = shutil.which("samestr") is not None
    have_instrain = shutil.which("inStrain") is not None
    have_r = shutil.which("Rscript") is not None

    return [
        {
            "capability_id": "core_measured_coverage",
            "kind": "core",
            "status": "implemented",
            "endpoint": "feature_available_in_donor",
            "claim_class": "measured_comparison",
            "note": (
                "recipient target ledger, per-material supported availability, frozen-denominator "
                "coverage, exhaustive set search, Pareto fronts, marginal value and sensitivity"
            ),
        },
        {
            "capability_id": "qualified_concern_ledger",
            "kind": "core",
            "status": "implemented",
            "endpoint": "donor_derived_adverse_event",
            "claim_class": "qualified_evidence",
            "note": "v05 pathogen and determinant evidence mapped through the §8.2 rule registry",
        },
        {
            "capability_id": "strain_resolution",
            "kind": "research_extension",
            "status": "unavailable",
            "endpoint": "donor_strain_engraftment",
            "claim_class": None,
            "missing_artifacts": [
                "SameStr v1.2025.111 or inStrain 1.10.0 on PATH"
                + ("" if (have_samestr or have_instrain) else " (not installed)"),
                "retained per-sample marker alignments or BAM/CRAM mappings for both materials",
            ],
            "fallback": "shared species are reported as not_resolved, never as same or different strain",
        },
        {
            "capability_id": "public_strain_baseline_v1",
            "kind": "predictive_model",
            "status": "not_trained",
            "endpoint": "donor_strain_engraftment",
            "claim_class": None,
            "missing_artifacts": [
                "Schmidt 2022 linked strain objects (doi:10.5281/zenodo.6611040) and MAGs (doi:10.5281/zenodo.5534163)",
                "donor/baseline/post-FMT triad mappings with prespecified follow-up windows",
                "grouped leave-one-study-out split recipe and a held-out calibration set",
            ],
            "fallback": "measured coverage only; no engraftment probability is produced",
            "note": "the training recipe is specified in §10.1; no cohort data is present in this repository",
        },
        {
            "capability_id": "ianiro2022_species",
            "kind": "published_adapter",
            "status": "artifact_blocked",
            "endpoint": "post_fmt_species_presence",
            "missing_artifacts": [
                "species-presence feature table rebuilt from the public mapped triads (PRJEB47909)",
                "the cohort/study grouping needed for leave-study-out evaluation",
            ],
            "fallback": "core measured coverage",
        },
        {
            "capability_id": "mozaic2026",
            "kind": "published_adapter",
            "status": "artifact_blocked",
            "endpoint": "community_convergence",
            "missing_artifacts": [
                "trained inference weights (the author driver expects local arrays that are not distributed)",
                "exact feature definitions, metadata and train/test splits",
                "resolution of the duplicate B010/B011 → SRR33557491 mapping in the exported metadata",
            ],
            "fallback": "core measured coverage; no convergence score is invented from the abstract",
        },
        {
            "capability_id": "zhang2026_replay",
            "kind": "published_adapter",
            "status": "artifact_blocked",
            "endpoint": "modeled_exchange_flux",
            "missing_artifacts": [
                "Rscript with ranger" + ("" if have_r else " (Rscript not installed)"),
                "rf_model_IBS.rds and the four committed predictors (average_distance, diversity_ratio, met_independence, log10_pre_abund)",
                "agora_species_103.qza and western_diet_gut_agora.qza, which the repository references but does not bundle",
                "an anvi'o workflow for met_independence",
            ],
            "fallback": "core functional comparison; no fold-change prediction is produced",
        },
        {
            "capability_id": "metabolic_complementarity_v1",
            "kind": "mechanistic_model",
            "status": "unavailable",
            "endpoint": "modeled_exchange_flux",
            "missing_artifacts": [
                "MICOM v0.39.1" + ("" if have_micom else " (not installed)"),
                "a pinned metabolic model collection, medium and solver",
            ],
            "fallback": "gene-capacity panels are compared directly and never called metabolite levels",
        },
        {
            "capability_id": "host_immune_compatibility",
            "kind": "external_modality",
            "status": "not_supplied",
            "endpoint": "clinical_response",
            "missing_artifacts": ["an actual GMLR or equivalent validated ex-vivo assay result"],
            "fallback": "absent; a host immune reaction is never reconstructed from microbial DNA",
        },
        {
            "capability_id": "pdf_only_mode",
            "kind": "core",
            "status": "not_implemented",
            "missing_artifacts": ["a PDF text extraction adapter with page locators"],
            "fallback": (
                "native results.json is required. This repository always produces it alongside the PDF, "
                "so the degraded path is unnecessary here"
            ),
        },
        {
            "capability_id": "followup_attribution",
            "kind": "core",
            "status": "not_implemented",
            "endpoint": "persistent_donor_strain",
            "missing_artifacts": [
                "post-exposure recipient samples",
                "strain resolution (see strain_resolution)",
                "unrelated-donor controls",
            ],
            "fallback": "no post-exposure claim of any kind is made",
        },
        {
            "capability_id": "clinical_matching_validation",
            "kind": "evidence_status",
            "status": "not_established_in_retrieved_evidence",
            "note": (
                "no prospectively validated donor-selection algorithm for post-acute COVID or ME/CFS was "
                "identified. Lau 2024 is a non-randomised insomnia study; Salonen 2023 is an 11-person "
                "autologous-controlled pilot with no significant fatigue benefit; COMEBACK has no posted "
                "results. This label does not suppress goal-relative matching."
            ),
        },
    ]


def core_ready(man: list[dict[str, Any]]) -> bool:
    return all(
        c["status"] == "implemented" for c in man if c["kind"] == "core" and c["capability_id"].startswith("core")
    )
