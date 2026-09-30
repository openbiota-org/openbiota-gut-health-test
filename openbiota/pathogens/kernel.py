"""Dependency-free sequence-status decision kernel (BUILD_SPEC_v05.0 §5.4).

This module is the adjudication boundary between alignment evidence and the
result schema. It is deliberately free of imports beyond the standard library
so that it can be executed, tested and audited on its own:

    python -m openbiota.pathogens.kernel

The spec supplies a reference implementation plus a fixed contract-check
suite. `sequence_status` here is *behaviourally identical* to that reference:
`reference_contract_checks()` runs the spec's own cases verbatim. Everything
else in this file is the surrounding machinery the spec requires to stay
orthogonal to the status decision — assay eligibility, contamination and
validation scope are separate mandatory fields and must never be folded into
one overloaded red/green verdict.

Three rules drive almost every surprising behaviour here, and all three exist
to stop the system inventing negatives:

1. `not_detected` is a *claim about an analysis*, not about an organism. It
   requires an eligible assay, a usable reference, a completed search and
   passed QC. Absent any of those the answer is `not_assessed`.
2. Support gates are conjunctive and include a specificity term. Five reads in
   a conserved gene never establish a species.
3. Protein-only homology can generate a candidate but can never populate
   nucleotide support counts, and so can never produce a species call or a
   negative.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Final

__all__ = [
    "ANALYSIS_STATUSES",
    "ASSAY_ELIGIBILITIES",
    "CONTAMINATION_STATUSES",
    "NormalizedEvidence",
    "SEQUENCE_STATUSES",
    "SUPPORTED_STATUSES",
    "SupportPolicy",
    "VALIDATION_SCOPES",
    "display_precedence",
    "reference_contract_checks",
    "self_test",
    "sequence_status",
]


# --------------------------------------------------------------------------- #
# Vocabularies (spec §4.4). Frozen: a value outside these sets is a bug, not
# a new category, and the kernel raises rather than guessing.
# --------------------------------------------------------------------------- #

ASSAY_ELIGIBILITIES: Final[frozenset[str]] = frozenset(
    {"eligible", "ineligible", "unknown"}
)
ANALYSIS_STATUSES: Final[frozenset[str]] = frozenset(
    {"completed", "not_run", "failed", "qc_failed"}
)
SEQUENCE_STATUSES: Final[frozenset[str]] = frozenset(
    {
        "supported_sequence",
        "marker_signal",
        "candidate_signal",
        "ambiguous_signal",
        "not_detected",
        "not_assessed",
    }
)
CONTAMINATION_STATUSES: Final[frozenset[str]] = frozenset(
    {
        "no_flag_under_applied_checks",
        "suspected",
        "confirmed_technical_artifact",
        "controls_unavailable",
        "not_evaluated",
    }
)
VALIDATION_SCOPES: Final[frozenset[str]] = frozenset(
    {
        "computational_research_rule",
        "analytical_assay_validated",
        "clinical_assay_validated",
    }
)

#: Statuses that represent positive sequence evidence at *some* strength. Used
#: for counting; `supported_sequence` and `marker_signal` are the only two
#: that the participant-facing summary promotes to a finding.
SUPPORTED_STATUSES: Final[frozenset[str]] = frozenset(
    {"supported_sequence", "marker_signal"}
)


@dataclass(frozen=True)
class SupportPolicy:
    """The three conjunctive support gates of a calling profile.

    Defaults are the spec's ``research_dna_v1`` values. They are engineering
    starting points to be benchmarked, not thresholds borrowed from a clinical
    guideline, and carry no sensitivity or limit-of-detection claim.
    """

    min_distinct_fragments: int = 5
    min_regions: int = 3
    min_informative_bases: int = 300


@dataclass(frozen=True)
class NormalizedEvidence:
    """Alignment evidence after the §5.3 normalisation rules have been applied.

    "Normalised" means the alignment identity/MAPQ filters, duplicate-family
    collapsing, independent-region merging and reference-specificity checks
    have *already* run. The kernel does no filtering of its own: it only
    adjudicates. Counts here are therefore qualifying counts, and inflating
    them upstream is the one way to make this kernel lie.
    """

    assay_eligibility: str = "unknown"
    analysis_status: str = "not_run"
    reference_usable: bool = False
    qc_pass: bool = False
    qualifying_fragments: int = 0
    nonoverlapping_informative_regions: int = 0
    informative_bases_covered: int = 0
    resolution_specific_support: bool = False
    only_unresolved_taxonomic_support: bool = False
    any_candidate_signal: bool = False
    marker_or_organelle_only: bool = False
    protein_only: bool = False


def sequence_status(e: NormalizedEvidence, policy: SupportPolicy) -> str:
    """Adjudicate one target's evidence into a `sequence_status` value.

    Behaviourally identical to the spec §5.4 reference kernel. Order matters:
    the assay/analysis/reference/QC gate runs before any evidence is
    considered, so a strong signal in an ineligible or failed analysis is
    `not_assessed` rather than a finding — the evidence survives elsewhere as
    an assay-inconsistent record, but it cannot unlock a negative panel.
    """
    if e.assay_eligibility not in ASSAY_ELIGIBILITIES:
        raise ValueError("Invalid assay eligibility")
    if e.analysis_status not in ANALYSIS_STATUSES:
        raise ValueError("Invalid analysis status")
    counts = (
        e.qualifying_fragments,
        e.nonoverlapping_informative_regions,
        e.informative_bases_covered,
    )
    if any(type(x) is not int or x < 0 for x in counts):
        raise ValueError("Support counts must be nonnegative integers")
    if (
        e.assay_eligibility != "eligible"
        or e.analysis_status != "completed"
        or not e.reference_usable
        or not e.qc_pass
    ):
        return "not_assessed"
    if e.protein_only:
        # Homology-only evidence: no nucleotide-negative and no species call.
        return "candidate_signal"
    if e.only_unresolved_taxonomic_support:
        return "ambiguous_signal"
    if e.qualifying_fragments == 0:
        return "candidate_signal" if e.any_candidate_signal else "not_detected"
    supported = (
        e.qualifying_fragments >= policy.min_distinct_fragments
        and e.nonoverlapping_informative_regions >= policy.min_regions
        and e.informative_bases_covered >= policy.min_informative_bases
        and e.resolution_specific_support
    )
    if not supported:
        return "candidate_signal"
    return "marker_signal" if e.marker_or_organelle_only else "supported_sequence"


def display_precedence(
    *,
    sequence: str,
    contamination: str,
    reference_gap_reason: str | None = None,
) -> tuple[str, str | None]:
    """Resolve what the participant-facing row shows (spec §4.4).

    Returns ``(display_status, qualifier)``.

    Two asymmetries are load-bearing. A confirmed technical artifact is
    *excluded* from the supported-organism count but its evidence stays
    visible as an artifact record. Merely *suspected* contamination is shown
    beside the finding and downgrades its prominence — it never silently
    becomes a negative, because "this might be lab carry-over" is not the
    same claim as "this organism is absent".
    """
    if sequence not in SEQUENCE_STATUSES:
        raise ValueError(f"Invalid sequence status: {sequence!r}")
    if contamination not in CONTAMINATION_STATUSES:
        raise ValueError(f"Invalid contamination status: {contamination!r}")
    if sequence == "not_assessed":
        return "not_assessed", reference_gap_reason
    if contamination == "confirmed_technical_artifact":
        return "technical_artifact", "excluded from supported findings"
    if contamination == "suspected" and sequence in SUPPORTED_STATUSES:
        return sequence, "suspected contamination; confirmation needed"
    if contamination == "controls_unavailable" and sequence in SUPPORTED_STATUSES:
        return sequence, "no extraction or library controls for this sample"
    return sequence, None


# --------------------------------------------------------------------------- #
# Contract checks
# --------------------------------------------------------------------------- #


def reference_contract_checks() -> int:
    """The spec §5.4 status-contract suite, verbatim.

    Kept byte-for-byte equivalent to the specification so that any drift in
    our kernel is caught here rather than discovered in a report.
    """
    p = SupportPolicy()
    ready = NormalizedEvidence("eligible", "completed", True, True)
    strong = replace(
        ready,
        qualifying_fragments=5,
        nonoverlapping_informative_regions=3,
        informative_bases_covered=300,
        resolution_specific_support=True,
    )
    cases = [
        (NormalizedEvidence(), "not_assessed"),
        (ready, "not_detected"),
        (strong, "supported_sequence"),
        (replace(strong, assay_eligibility="unknown"), "not_assessed"),
        (replace(strong, assay_eligibility="ineligible"), "not_assessed"),
        (replace(strong, analysis_status="failed"), "not_assessed"),
        (replace(strong, reference_usable=False), "not_assessed"),
        (replace(strong, qc_pass=False), "not_assessed"),
        (replace(strong, qualifying_fragments=4), "candidate_signal"),
        (replace(strong, nonoverlapping_informative_regions=1), "candidate_signal"),
        (replace(strong, informative_bases_covered=299), "candidate_signal"),
        (replace(strong, resolution_specific_support=False), "candidate_signal"),
        (replace(strong, marker_or_organelle_only=True), "marker_signal"),
        (replace(strong, only_unresolved_taxonomic_support=True), "ambiguous_signal"),
        (replace(strong, protein_only=True), "candidate_signal"),
        (replace(ready, any_candidate_signal=True), "candidate_signal"),
    ]
    for evidence, expected in cases:
        assert sequence_status(evidence, p) == expected
    try:
        sequence_status(replace(ready, qualifying_fragments=-1), p)
    except ValueError:
        pass
    else:
        raise AssertionError("Negative counts accepted")
    return len(cases) + 1


def self_test() -> int:
    """Run the spec suite plus this module's own precedence checks."""
    checks = reference_contract_checks()

    # A confirmed artifact never counts as a supported organism.
    status, note = display_precedence(
        sequence="supported_sequence", contamination="confirmed_technical_artifact"
    )
    assert status == "technical_artifact" and note is not None
    checks += 1

    # Suspicion qualifies a finding; it does not delete it.
    status, note = display_precedence(
        sequence="supported_sequence", contamination="suspected"
    )
    assert status == "supported_sequence" and note is not None
    checks += 1

    # Historical samples without blanks are qualified, never "checks passed".
    status, note = display_precedence(
        sequence="marker_signal", contamination="controls_unavailable"
    )
    assert status == "marker_signal" and note is not None
    checks += 1

    # A reference gap surfaces its reason on the not-assessed row.
    status, note = display_precedence(
        sequence="not_assessed",
        contamination="not_evaluated",
        reference_gap_reason="no_usable_reference",
    )
    assert status == "not_assessed" and note == "no_usable_reference"
    checks += 1

    # An artifact flag cannot manufacture a finding out of a non-detection.
    status, _ = display_precedence(
        sequence="not_detected", contamination="controls_unavailable"
    )
    assert status == "not_detected"
    checks += 1

    for bad in ("Detected", "positive", ""):
        try:
            display_precedence(sequence=bad, contamination="not_evaluated")
        except ValueError:
            checks += 1
        else:  # pragma: no cover - guard
            raise AssertionError(f"accepted invalid sequence status {bad!r}")

    # A marker-only reference caps the result at marker_signal even when the
    # evidence would otherwise be genome-supported.
    p = SupportPolicy()
    marker = NormalizedEvidence(
        "eligible", "completed", True, True,
        qualifying_fragments=500,
        nonoverlapping_informative_regions=40,
        informative_bases_covered=50_000,
        resolution_specific_support=True,
        marker_or_organelle_only=True,
    )
    assert sequence_status(marker, p) == "marker_signal"
    checks += 1

    # A stricter profile can demote what the default would support.
    strict = SupportPolicy(min_distinct_fragments=50, min_regions=10,
                           min_informative_bases=5_000)
    borderline = NormalizedEvidence(
        "eligible", "completed", True, True,
        qualifying_fragments=5,
        nonoverlapping_informative_regions=3,
        informative_bases_covered=300,
        resolution_specific_support=True,
    )
    assert sequence_status(borderline, p) == "supported_sequence"
    assert sequence_status(borderline, strict) == "candidate_signal"
    checks += 2

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{reference_contract_checks()} reference contract checks passed")
    print(f"{self_test()} kernel checks passed")
