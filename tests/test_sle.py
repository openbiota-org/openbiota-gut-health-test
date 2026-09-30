"""The lupus profile, and the restraint it is built out of.

A twenty-five-study audit concluded that this literature does not support a
lupus-specific stool score. The profile that exists is therefore deliberately
small and deliberately misnamed away from "lupus": it scores a shared
inflammatory pattern, excludes the organism lupus is famous for, and declares
the strain measurement it cannot make.

These tests exist because every one of those restraints is the kind of thing
a later well-meaning edit removes. Scoring Ruminococcus gnavus here, or
raising its specificity, or dropping the serology caveat, would each turn an
honest non-specific descriptor back into a claim the evidence cannot carry.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from openbiota.profiles import load_profile_set

PROFILES = Path("profiles")


@pytest.fixture(scope="module")
def sle():
    ps = load_profile_set(PROFILES)
    return next((p for p in ps.profiles if p.name == "sle"), None)


@pytest.fixture(scope="module")
def profile_set():
    return load_profile_set(PROFILES)


def test_the_profile_exists_and_is_labelled_research_only(sle) -> None:
    assert sle is not None, "lupus profile is missing"
    assert sle.status == "research_only"
    # 16S-dominant evidence with no independent shotgun replication of the
    # scored directions caps maturity at one small-cohort tier.
    assert sle.evidence_maturity in {"P1", "P2"}


def test_the_label_does_not_claim_lupus_specificity(sle) -> None:
    """The name is the finding: stool cannot separate lupus from Sjögren's."""
    label = sle.label.lower()
    assert "not lupus-specific" in label or "shared inflammatory" in label
    # And it must not be dressed up as a diagnostic or a risk score.
    for word in ("diagnos", "risk score", "probability", "screen for"):
        assert word not in label


def test_ruminococcus_gnavus_is_not_scored(sle) -> None:
    """The flagship finding is strain-level and a transient bloom.

    Species abundance is measured accurately and answers a different
    question. It is also already a scored feature in fifteen other profiles,
    so adding it here would mostly import general dysbiosis.
    """
    scored = [
        f.name
        for m in sle.modules.values()
        if m.bound
        for f in m.features
    ]
    assert not any(
        "gnavus" in n.lower() and "strain" not in n.lower() for n in scored
    ), f"R. gnavus must not be a scored species feature; scored: {scored}"


def test_the_strain_finding_is_declared_but_not_computable(sle) -> None:
    """The one lupus-specific mechanism, stated as unavailable rather than dropped."""
    lane = next(
        (m for m in sle.modules.values() if m.evidence_type == "STRAIN"), None
    )
    assert lane is not None, "the strain lane must be declared"
    assert not lane.bound, "no strain engine exists, so this must not score"
    assert lane.weight_in_combined == 0.0
    names = [f.name for f in lane.features]
    assert any("lipoglycan" in n.lower() for n in names)


def test_every_scored_feature_carries_low_specificity(sle) -> None:
    """Both scored directions are pan-rheumatic markers.

    Faecalibacterium depletion and Streptococcus enrichment appear across
    rheumatic and inflammatory disease. A high specificity weight here would
    assert a discrimination the head-to-head study failed to find.
    """
    for module in sle.modules.values():
        if not module.bound:
            continue
        for feature in module.features:
            assert feature.s <= 0.3, (
                f"{feature.name}: specificity {feature.s} overstates a "
                "shared inflammatory marker"
            )


def test_competing_inflammatory_patterns_are_scored_beside_it(sle, profile_set) -> None:
    """Specificity is handled by scoring the alternatives, not by assertion."""
    challenged = set(sle.challenge_for)
    assert {"crohns", "uc", "ibd"} <= challenged
    known = {p.name for p in profile_set.profiles}
    assert challenged <= known, f"unknown challenge targets: {challenged - known}"


def test_the_conflicted_genus_declares_its_conflict(sle) -> None:
    """A Faecalibacterium aggregate spans species that do not move together."""
    feature = next(
        f
        for m in sle.modules.values()
        for f in m.features
        if f.name == "Faecalibacterium"
    )
    assert feature.level == "genus"
    assert feature.direction_conflict == "material"


def test_the_unscorable_evidence_is_recorded_not_dropped(sle) -> None:
    """Four things the audit found that this assay cannot turn into a score."""
    blob = " ".join(str(s) for s in sle.unbound_sources).lower()
    # The only Bonferroni-surviving taxonomic finding — and a family, which
    # this schema does not score.
    assert "ruminococcaceae" in blob
    # The founding R. gnavus observation, kept as provenance.
    assert "gnavus" in blob
    # The historical F/B ratio claim, at zero weight.
    assert "bacteroidetes" in blob
    # Temporal instability: the one design the evidence endorses, and one
    # this single-specimen report cannot compute at all.
    assert "instability" in blob or "temporal" in blob


def test_the_report_says_the_real_lupus_tests_are_blood_tests(sle) -> None:
    """The most useful sentence in the profile.

    Every validated biomarker in this literature is serological. A reader who
    takes nothing else away should take that.
    """
    text = " ".join(sle.caveats).lower() + " " + " ".join(
        str(s) for s in sle.unbound_sources
    ).lower()
    assert "blood test" in text or "serolog" in text
    assert "dsdna" in text or "autoantibody" in text


def test_medication_confounds_are_declared(sle) -> None:
    """PPIs erased the diversity deficit and out-explained every clinical variable."""
    names = " ".join(c.name for c in sle.confounders).lower()
    assert "proton pump" in names
    assert "hydroxychloroquine" in names
    # No diversity term may be scored, because PPI status governs it.
    for module in sle.modules.values():
        for feature in module.features:
            assert "diversity" not in feature.name.lower()
            assert "shannon" not in feature.name.lower()


def test_no_probiotic_recommendation_can_be_read_out_of_it(sle) -> None:
    """L. reuteri improved one lupus model and worsened another.

    A depletion scored anywhere near this condition must not read as an
    argument for supplementing the organism.
    """
    text = " ".join(sle.caveats).lower()
    assert "lactobacillus" in text
    assert "immunosuppress" in text or "bacteraemia" in text
    scored = [f.name.lower() for m in sle.modules.values() if m.bound for f in m.features]
    assert not any("lactobacillus" in n or "reuteri" in n for n in scored)
