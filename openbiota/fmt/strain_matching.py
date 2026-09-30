"""Strain matching for donor comparison: three separate capabilities.

Spec §10.1 insists these are kept apart, and the reason is that conflating
them is how a comparison becomes a promise:

1. **Baseline comparison** — how related are the donor's and the recipient's
   populations of each organism, right now? Needs only the specimens in hand.
   This is computable today and is computed here.

2. **Observed engraftment** — did the donor's population actually establish?
   Needs a dated recipient follow-up specimen. With no follow-up the answer is
   `followup_required`, never "no engraftment" and never zero.

3. **Prospective establishment prediction** — will it establish? Needs a
   frozen, evaluated model. Without one the answer is `model_unavailable`, and
   the baseline comparison stands on its own rather than being deleted.

The distinction that matters most for reading the output: a donor strain that
differs from the recipient's is a *candidate for changing something*, and one
that is indistinguishable adds nothing new. Neither statement is a forecast.
Phylogenetic distance is not converted into a probability anywhere in this
module, and `match_score` keeps its original meaning — target availability and
complementarity — rather than being quietly renamed (§10.1).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.resolution import strains as strain_mod

#: Where the marker lane leaves a sample's fingerprints, relative to the
#: material's input root.
CONSENSUS_GLOB: Final = "strain/consensus_markers/*.json.bz2"

PENDING: Final = "pending"
COMPLETED: Final = "completed"
UNAVAILABLE_INPUT: Final = "strain_inputs_unavailable"
FOLLOWUP_REQUIRED: Final = "followup_required"
MODEL_UNAVAILABLE: Final = "model_unavailable"

MATCH_SCORE_SEMANTICS: Final = "target_availability_and_complementarity"


@dataclass(slots=True)
class StrainMatching:
    """The nullable strain block attached to a match result."""

    baseline_comparison: dict[str, Any] = field(
        default_factory=lambda: {"status": PENDING, "evidence_ids": [], "results": None}
    )
    observed_engraftment: dict[str, Any] = field(
        default_factory=lambda: {
            "status": FOLLOWUP_REQUIRED,
            "followup_sample_ids": [],
            "results": None,
            "why": (
                "Observed engraftment compares a dated recipient follow-up specimen against "
                "the donor and the recipient's own baseline. No follow-up specimen was "
                "supplied, so this is unmeasured \u2014 which is not the same as no "
                "engraftment, and is not zero."
            ),
        }
    )
    prospective_establishment: dict[str, Any] = field(
        default_factory=lambda: {
            "status": MODEL_UNAVAILABLE,
            "model_id": None,
            "predictions": None,
            "why": (
                "No frozen, evaluated establishment model is installed. A phylogenetic "
                "distance is not a probability, and a model trained on individual donors "
                "could not be applied to a donor combination without combination-specific "
                "validation. The baseline comparison below is unaffected."
            ),
        }
    )
    disease_transmission_probability: None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "baseline_comparison": self.baseline_comparison,
            "observed_engraftment": self.observed_engraftment,
            "prospective_establishment": self.prospective_establishment,
            "disease_transmission_probability": None,
        }


def _consensus_path(root: Path) -> Path | None:
    hits = sorted(root.glob(CONSENSUS_GLOB))
    return hits[0] if hits else None


def load_for(materials: Mapping[str, Path]) -> tuple[dict[str, Any], list[str]]:
    """Load fingerprints for every person with a marker lane result.

    Returns `(fingerprints_by_person, missing_person_ids)`. A person without
    marker artifacts is named rather than silently skipped, so the reason a
    comparison is incomplete stays visible.
    """
    loaded: dict[str, Any] = {}
    missing: list[str] = []
    databases: set[str] = set()
    for person, root in materials.items():
        path = _consensus_path(Path(root))
        if path is None:
            missing.append(person)
            continue
        fingerprints, database = strain_mod.load_fingerprints(path, sample_id=person)
        databases.add(database)
        loaded[person] = fingerprints
    if len(databases) > 1:
        # Mixing marker databases would compare incomparable coordinates.
        raise ValueError(
            f"marker fingerprints span more than one database release: {sorted(databases)}"
        )
    return loaded, sorted(missing)


def baseline(
    *,
    recipient: str,
    donors: Sequence[str],
    materials: Mapping[str, Path],
) -> StrainMatching:
    """Per-donor baseline strain comparison against the recipient.

    Each donor is evaluated separately before any set is considered, because a
    set's strain relatedness is not the average of its members' (§10.1).
    """
    out = StrainMatching()
    fingerprints, missing = load_for(materials)

    if recipient not in fingerprints:
        out.baseline_comparison = {
            "status": UNAVAILABLE_INPUT,
            "evidence_ids": [],
            "results": None,
            "why": (
                f"No marker fingerprints for the recipient ({recipient}), so no donor "
                "comparison is possible. Run the marker lane for this specimen."
            ),
            "missing_marker_lane": missing,
        }
        return out

    per_donor: dict[str, Any] = {}
    evidence_ids: list[str] = []
    for donor in donors:
        if donor not in fingerprints:
            per_donor[donor] = {
                "status": UNAVAILABLE_INPUT,
                "why": "no marker fingerprints for this donor; other donors are unaffected",
                "organisms_compared": 0,
                "comparisons": [],
            }
            continue
        comparisons = strain_mod.compare_all(fingerprints[recipient], fingerprints[donor])
        distinct = [c for c in comparisons if c.status == "distinct_populations"]
        same = [c for c in comparisons if c.status == "indistinguishable_over_compared_sites"]
        close = [c for c in comparisons if c.status == "closely_related"]
        thin = [c for c in comparisons if c.status == "insufficient_callable_overlap"]
        donor_only = sorted(
            sgb for sgb, fp in fingerprints[donor].items()
            if fp.resolved and sgb not in fingerprints[recipient]
        )
        evidence_ids.extend(f"{donor}:{c.sgb}" for c in comparisons)
        per_donor[donor] = {
            "status": COMPLETED,
            "organisms_compared": len(comparisons),
            "distinct_populations": len(distinct),
            "closely_related": len(close),
            "indistinguishable": len(same),
            "insufficient_overlap": len(thin),
            "organisms_resolved_in_donor_only": donor_only,
            "n_organisms_resolved_in_donor_only": len(donor_only),
            "comparisons": [c.to_json() for c in comparisons],
            "plain": (
                f"Of {len(comparisons)} organisms resolved in both specimens, "
                f"{len(distinct)} are carried as distinct populations, {len(close)} are "
                f"closely related and {len(same)} are indistinguishable at marker "
                f"resolution. {len(donor_only)} further organism(s) were strain-resolved in "
                f"{donor} but not in the recipient. A distinct population is a candidate for "
                "changing something; an indistinguishable one would add nothing new. Neither "
                "is a prediction that it would establish."
            ),
        }

    out.baseline_comparison = {
        "status": COMPLETED,
        "recipient": recipient,
        "method": "metaphlan_4.1.1_marker_consensus",
        "per_donor": per_donor,
        "evidence_ids": evidence_ids,
        "missing_marker_lane": missing,
        "limits": [
            strain_mod.DOMINANT_CONSENSUS_LIMIT,
            strain_mod.MARKER_SCOPE_LIMIT,
        ],
        "note": (
            "Donors are compared to the recipient one at a time. A donor combination's "
            "strain relatedness is not the average of its members', and no combination "
            "figure is produced here."
        ),
    }
    return out


def self_test() -> int:
    """Check the capability states stay separate. Returns failures."""
    failures = 0
    empty = StrainMatching().to_json()

    # Absent follow-up must not read as zero engraftment.
    engraftment = empty["observed_engraftment"]
    if engraftment["status"] != FOLLOWUP_REQUIRED:
        failures += 1
    if engraftment["results"] is not None:
        failures += 1
    if "not zero" not in engraftment["why"]:
        failures += 1

    # Absent model must not zero the comparison.
    prospective = empty["prospective_establishment"]
    if prospective["status"] != MODEL_UNAVAILABLE:
        failures += 1
    if prospective["predictions"] is not None:
        failures += 1
    if "not a probability" not in prospective["why"]:
        failures += 1

    if empty["disease_transmission_probability"] is not None:
        failures += 1

    # A recipient without a marker lane fails loudly, not silently.
    out = baseline(recipient="R", donors=["D"], materials={"R": Path("/nonexistent")})
    if out.baseline_comparison["status"] != UNAVAILABLE_INPUT:
        failures += 1
    if out.baseline_comparison["results"] is not None:
        failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
