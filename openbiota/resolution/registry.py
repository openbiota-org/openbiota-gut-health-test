"""Target definitions, compiled in one place rather than scattered.

Every disease, pathogen, functional and intervention module declares what it
needs resolved. The compiler expands those declarations into concrete assays,
checks that each target's references are actually ready, and builds a census so
an unrun target stays visible instead of quietly reading as a negative
(spec §9.2).

The rheumatoid-arthritis targets are defined here first because RA is the case
that exposed why this layer is needed. The old profile printed a single
"RA-associated *P. copri* clade" module, which is not a real biological entity:

* Nii et al. compared RA and healthy *S. copri* isolates and found both sat
  mostly in **clade A**. The difference associated with arthritis-promoting
  activity was a ~100-kb accessory element, **CTnPc** — not the clade.
* **Pc-p27**, the antigen from Pianta et al., also occurs in isolates from
  healthy people.
* Manghi et al. resolved the whole complex into 13 species-level clades and
  found none associated with any condition, RA included.

So "clade A present" is not an RA finding, and neither is the antigen alone.
The honest decomposition is four separate targets with separate evidence, which
is what this module registers. Three of them are genuinely measurable from the
data we hold; CTnPc needs coordinate curation against the public assembly
before it can be called complete, and says so.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

from .schema import (
    ANIMAL_MECHANISTIC,
    HUMAN_ASSOCIATION,
    IDENTITY_KINDS,
    IDENTITY_ONLY,
    BiologicalEvidence,
)


@dataclass(frozen=True, slots=True)
class Target:
    """One thing the system is asked to resolve, and on what terms."""

    target_id: str
    identity_kind: str
    label: str
    #: Which module/profile asked for it, for the migration matrix.
    requested_by: tuple[str, ...]
    required_resolution: str
    #: Run this screen even when the parent species looks absent. An
    #: incomplete marker catalogue must not create a permanent blind spot
    #: (spec §4.1).
    screen_independent_of_parent_profile: bool = False
    parent_complex: str | None = None
    reference_anchor: str | None = None
    #: `ready` | `coordinates_require_reconciliation` | `curation_required`
    #: | `patent_derived_candidate`
    reference_readiness: str = "ready"
    canonical_region_coordinates: str | None = None
    allow_generic_mobile_element_match: bool = False
    allow_species_fallback: bool = False
    positive_reference_groups: tuple[str, ...] = ()
    negative_reference_groups: tuple[str, ...] = ()
    result_when_incomplete: str = "preserve_partial_evidence_and_reason"
    biological_evidence: tuple[BiologicalEvidence, ...] = ()
    #: Interpretation guards, printed with the result.
    candidate_mechanism: bool = False
    broad_clade_is_pathogenicity: bool = False
    validated_human_risk_predictor: bool = False
    source_doi: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if self.identity_kind not in IDENTITY_KINDS:
            raise ValueError(f"{self.target_id}: unknown identity_kind {self.identity_kind!r}")

    @property
    def operational(self) -> bool:
        """Can this target produce a complete call, as opposed to partial?

        A target whose canonical coordinates are still null can produce
        explicitly partial candidate evidence while curation proceeds. It
        cannot pass as a complete-region assay (spec §9.2).
        """
        return self.reference_readiness == "ready"

    @property
    def region_delimited(self) -> bool:
        """Coordinates exist, but the assay is not validated.

        A delimited candidate region is real progress and must be visible, yet
        it is still not a `ready` reference: it can support partial candidate
        evidence with breadth and discriminatory loci, never a confirmed call.
        """
        return self.reference_readiness == "candidate_region_delimited"

    def to_json(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "identity_kind": self.identity_kind,
            "label": self.label,
            "requested_by": list(self.requested_by),
            "required_resolution": self.required_resolution,
            "screen_independent_of_parent_profile": self.screen_independent_of_parent_profile,
            "parent_complex": self.parent_complex,
            "reference_anchor": self.reference_anchor,
            "reference_readiness": self.reference_readiness,
            "canonical_region_coordinates": self.canonical_region_coordinates,
            "allow_generic_mobile_element_match": self.allow_generic_mobile_element_match,
            "allow_species_fallback": self.allow_species_fallback,
            "positive_reference_groups": list(self.positive_reference_groups),
            "negative_reference_groups": list(self.negative_reference_groups),
            "result_when_incomplete": self.result_when_incomplete,
            "biological_evidence": [e.to_json() for e in self.biological_evidence],
            "interpretation": {
                "candidate_mechanism": self.candidate_mechanism,
                "broad_clade_is_pathogenicity": self.broad_clade_is_pathogenicity,
                "validated_human_risk_predictor": self.validated_human_risk_predictor,
            },
            "operational": self.operational,
            "region_delimited": self.region_delimited,
            "source_doi": self.source_doi,
            "notes": self.notes,
        }


# --------------------------------------------------------------------------- #
# Rheumatoid arthritis (spec §8.1)
# --------------------------------------------------------------------------- #

RA_TARGETS: Final[tuple[Target, ...]] = (
    Target(
        target_id="ra.segatella_complex",
        identity_kind="species_complex",
        label="Segatella (Prevotella) copri complex, resolved to clade",
        requested_by=("profile:ra",),
        required_resolution="clade_level_composition",
        parent_complex="segatella_copri_complex",
        reference_anchor="mpa_vJun23_CHOCOPhlAnSGB_202403",
        reference_readiness="ready",
        source_doi="10.1016/j.chom.2023.09.013",
        biological_evidence=(
            BiologicalEvidence(
                kind=HUMAN_ASSOCIATION,
                statement=(
                    "Expansion of the species was reported in new-onset untreated RA, but "
                    "replication has been poor and, at clade resolution, none of the 13 "
                    "species of the complex was associated with any condition across 1,635 "
                    "cases and 1,854 controls."
                ),
                studied_entity="Segatella copri complex, clades A-M",
                endpoint="case-control community association",
                sources=(
                    "Scher 2013 eLife doi:10.7554/eLife.01202",
                    "Manghi 2023 Cell Host Microbe doi:10.1016/j.chom.2023.09.013",
                ),
                boundary=(
                    "Composition only. A clade is taxonomy, not a risk category, and clade A "
                    "carries the organisms found in both RA patients and healthy controls."
                ),
            ),
        ),
        # The guard that the old module violated by construction.
        broad_clade_is_pathogenicity=False,
        validated_human_risk_predictor=False,
        notes=(
            "Replaces the former 'RA-associated P. copri clade' module, which named an "
            "entity that does not exist in the literature."
        ),
    ),
    Target(
        target_id="ra.ctnpc",
        identity_kind="accessory_region",
        label="CTnPc accessory element (~100 kb) in Segatella copri",
        requested_by=("profile:ra",),
        required_resolution="region_with_carrier_context",
        parent_complex="segatella_copri_complex",
        screen_independent_of_parent_profile=True,
        reference_anchor="GCA_026015645.1",
        # Delimited by presence/absence against the study's own comparators
        # rather than by the RAST identifiers, which do not map onto the public
        # assembly. Still not `ready`: the region is a candidate supported by
        # two positive isolates, not a validated assay.
        reference_readiness="candidate_region_delimited",
        canonical_region_coordinates="JAPDUN010000001.1:399539-497219",
        # A conjugative element is made of integrases and transfer genes that
        # every such element shares. Matching those would call CTnPc in any
        # sample carrying any mobile element.
        allow_generic_mobile_element_match=False,
        allow_species_fallback=False,
        positive_reference_groups=("RA-N001-13", "RAP9-13"),
        negative_reference_groups=("N115-17", "H012_6"),
        source_doi="10.1136/ard-2022-222881",
        biological_evidence=(
            BiologicalEvidence(
                kind=ANIMAL_MECHANISTIC,
                statement=(
                    "RA and healthy S. copri isolates largely shared clade A; isolates "
                    "carrying the ~100-kb CTnPc accessory content showed greater "
                    "arthritis-promoting activity in the tested experimental system."
                ),
                studied_entity="S. copri isolates RA-N001-13, RAP9-13 (positive); N115-17, H012_6 (negative)",
                endpoint="arthritis-promoting activity in a susceptible animal model",
                sources=("Nii 2023 Ann Rheum Dis doi:10.1136/ard-2022-222881",),
                boundary=(
                    "CTnPc is not a proven sufficient causal determinant and not a validated "
                    "human risk marker. N115-17 is an RA-origin isolate that is CTnPc-negative, "
                    "so patient origin is not the call rule."
                ),
            ),
        ),
        candidate_mechanism=True,
        validated_human_risk_predictor=False,
        notes=(
            "The published supplement lists 42 annotated CDSs (fig|165179.43.peg.373-414) "
            "inside the region, which are RAST annotation-server identifiers and do not map "
            "onto NCBI's annotation of the public assembly. The region was therefore "
            "delimited the other way, by the study's own definition: content present in both "
            "CTnPc-positive isolates and absent from both negatives (one of which, N115-17, "
            "is RA-origin). That yields JAPDUN010000001.1:399539-497219 \u2014 97,681 bp "
            "carrying 83 shared accessory CDSs, against a published description of ~100 kb "
            "with 90-100 CDSs. "
            "A first attempt used each isolate's largest accessory run instead; those proved "
            "entirely unrelated to each other (zero orthologs), so the size match to the "
            "published description was coincidence. Requiring presence in both positives is "
            "what separated the element from isolate-specific content. "
            "It stays non-operational: two positives cannot fully separate the element from "
            "content they share for unrelated reasons, the boundaries follow annotated CDS "
            "edges rather than the element's termini, and a read hit still requires breadth, "
            "depth and carrier evidence. See refs/ctnpc/CTNPC_REGION.json."
        ),
    ),
    Target(
        target_id="ra.pc_p27_antigen",
        identity_kind="locus",
        label="Pc-p27 antigen locus (Segatella copri)",
        requested_by=("profile:ra",),
        required_resolution="locus_with_allele_evidence",
        parent_complex="segatella_copri_complex",
        screen_independent_of_parent_profile=True,
        reference_readiness="curation_required",
        allow_species_fallback=False,
        source_doi="10.1002/art.40003",
        biological_evidence=(
            BiologicalEvidence(
                kind=HUMAN_ASSOCIATION,
                statement=(
                    "T-cell and antibody reactivity to the Pc-p27 peptide was reported in a "
                    "subset of RA patients."
                ),
                studied_entity="Pc-p27 peptide of P. copri",
                endpoint="immune reactivity in a patient subset",
                sources=(
                    "Pianta 2017 Arthritis Rheumatol doi:10.1002/art.40003",
                ),
                boundary=(
                    "The antigen also occurs in isolates from healthy people, and attempts to "
                    "reproduce the antibody finding returned modest or null results with "
                    "reactivity varying by strain. Presence of the locus is not an RA finding."
                ),
            ),
        ),
        candidate_mechanism=True,
        validated_human_risk_predictor=False,
    ),
    Target(
        target_id="ra.subdoligranulum_d8_marker",
        identity_kind="allele",
        label="Candidate D8-associated marker (Subdoligranulum)",
        requested_by=("profile:ra",),
        required_resolution="candidate_marker_with_specificity_controls",
        reference_readiness="patent_derived_candidate",
        allow_species_fallback=False,
        negative_reference_groups=("H3",),
        source_doi="10.1002/art.42381",
        biological_evidence=(
            BiologicalEvidence(
                kind=ANIMAL_MECHANISTIC,
                statement=(
                    "A specific Subdoligranulum isolate, unlike a comparison isolate from the "
                    "same study, produced RA-relevant joint and immune findings in mice."
                ),
                studied_entity="Subdoligranulum isolate 7 (proposed S. didolesgii D8); isolate 1 (H3) as comparator",
                endpoint="joint swelling and RA-relevant immune findings in mice",
                sources=(
                    "Chriswell 2022 Sci Transl Med PMC9804515",
                    "Patent WO2024006983A1 (candidate marker regions)",
                ),
                boundary=(
                    "The marker regions are patent-derived and not independently validated. "
                    "Region A is AT-rich, so specificity against H3 and low-complexity decoys "
                    "must be established first. A positive marker with unresolved backbone is "
                    "not 'D8 present', and a marker-negative shallow sample is not 'no "
                    "arthritogenic strain'."
                ),
            ),
        ),
        candidate_mechanism=True,
        validated_human_risk_predictor=False,
        notes=(
            "Report only as 'candidate D8-associated marker evidence' until the official "
            "sequence listing and complete D8/H3 comparator genomes are reconciled."
        ),
    ),
)

#: The community-resemblance percentile that the legacy RA profile produces.
#: Kept deliberately, as a descriptive measure, and kept separate from every
#: target above (spec §8.1, acceptance 37).
RA_LEGACY_SCORE_MEANING: Final = (
    "The RA profile percentile is community resemblance to a published case-control "
    "pattern, on the profile's frozen reference scale. It is not a detected transferable "
    "agent and not a risk of developing or transmitting rheumatoid arthritis. In the "
    "supplied samples most of the difference between donors came from Bacteroides "
    "vulgatus being undetected rather than from any RA-associated organism: a depletion "
    "term cannot become a detected causal agent."
)

#: Every registered target, by id.
TARGETS: Final[Mapping[str, Target]] = {t.target_id: t for t in RA_TARGETS}


@dataclass(frozen=True, slots=True)
class MigrationRow:
    """One line of the additive migration matrix (spec §9.2)."""

    profile_id: str
    current_feature_type: str
    required_resolution: str
    supporting_source: str
    new_adapter: str
    reference_readiness: str
    downstream_consumers: tuple[str, ...]
    fixture: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "current_feature_type": self.current_feature_type,
            "required_resolution": self.required_resolution,
            "supporting_source": self.supporting_source,
            "new_adapter": self.new_adapter,
            "reference_readiness": self.reference_readiness,
            "downstream_consumers": list(self.downstream_consumers),
            "fixture": self.fixture,
        }


def targets_for(requester: str) -> list[Target]:
    """Targets a given profile or module declared."""
    return [t for t in TARGETS.values() if requester in t.requested_by]


def unresolved_placeholders() -> list[str]:
    """Targets that cannot yet produce a complete call, and must say so.

    Acceptance test 1 forbids a hidden strain-only placeholder outside the
    census; this is the list that keeps them visible.
    """
    return sorted(t.target_id for t in TARGETS.values() if not t.operational)


def registry_json() -> dict[str, Any]:
    return {
        "targets": [t.to_json() for t in TARGETS.values()],
        "n_targets": len(TARGETS),
        "non_operational": unresolved_placeholders(),
        "ra_legacy_score_meaning": RA_LEGACY_SCORE_MEANING,
    }


def self_test() -> int:
    """Guard the interpretation rules that RA got wrong. Returns failures."""
    failures = 0

    # Acceptance 33: clade A or the antigen alone cannot be a pathogenic call.
    for tid in ("ra.segatella_complex", "ra.pc_p27_antigen"):
        target = TARGETS[tid]
        if target.broad_clade_is_pathogenicity:
            failures += 1
        if target.validated_human_risk_predictor:
            failures += 1

    # Acceptance 35: generic mobile-element hits cannot establish CTnPc, and
    # missing coordinates cannot masquerade as a complete assay.
    ctnpc = TARGETS["ra.ctnpc"]
    if ctnpc.allow_generic_mobile_element_match:
        failures += 1
    # Coordinates may now exist (delimited by presence/absence), but having
    # them must never promote the target to a validated assay.
    if ctnpc.operational:
        failures += 1
    if ctnpc.canonical_region_coordinates and not ctnpc.region_delimited:
        failures += 1
    if ctnpc.region_delimited and ctnpc.operational:
        failures += 1
    # Acceptance 34: both study groups present, so origin is not the rule.
    if not (ctnpc.positive_reference_groups and ctnpc.negative_reference_groups):
        failures += 1

    # Acceptance 36: D8 stays candidate-status.
    d8 = TARGETS["ra.subdoligranulum_d8_marker"]
    if d8.operational or d8.validated_human_risk_predictor:
        failures += 1

    # Acceptance 17: the accessory screens run even if the parent looks absent.
    for tid in ("ra.ctnpc", "ra.pc_p27_antigen"):
        if not TARGETS[tid].screen_independent_of_parent_profile:
            failures += 1

    # No target may fall back to a species-level answer for a strain question.
    if any(t.allow_species_fallback for t in TARGETS.values()):
        failures += 1

    # Every target must carry its evidence with a stated boundary.
    for target in TARGETS.values():
        if not target.biological_evidence:
            failures += 1
        for evidence in target.biological_evidence:
            if evidence.kind != IDENTITY_ONLY and not evidence.boundary:
                failures += 1
            if not evidence.studied_entity:
                failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
