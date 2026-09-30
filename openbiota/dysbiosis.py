"""General dysbiosis anchor: GMWI2.

Every profile score is reported beside a general gut-health index, and the
reason is specificity. If the index says the community is generally
disturbed, a high score on any profile most likely reflects that disturbance
rather than anything condition-specific — and the report has to say so
instead of presenting the profile score alone. The cross-profile specificity
matrix (see `profilevalidate`) is what makes this anchor interpretable: it
shows that the profiles' features overlap heavily with general dysbiosis.

GMWI2 (Chang & Gupta et al., Nat Commun 2024) is a Lasso-penalised logistic
regression on the presence/absence of MetaPhlAn 3 clades, trained on 8,069
health-status-labelled stool metagenomes from 54 studies. Positive scores
lean healthy, negative lean unhealthy, and the magnitude is the log-odds.
The published model has 95 non-zero coefficients; they are stored verbatim in
`openbiota/data/gmwi2_model.json` and applied exactly as the reference
implementation does — normalise by the classified fraction, threshold
presence at 1e-5, dot with the coefficients. No sklearn dependency is needed
for a dot product.

The model was trained on MetaPhlAn 3 profiles against the CHOCOPhlAn 201901
markers, which is the database this pipeline runs. Applying it to output from
a different profiler version would not be valid, and the module refuses to.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from typing import Any, Final

#: MetaPhlAn database families the published model is valid for.
COMPATIBLE_DATABASES: Final = ("mpa_v30_CHOCOPhlAn_201901", "mpa_v31_CHOCOPhlAn_201901")

#: Interpretation bands. The authors describe |GMWI2| > 1 as a strong call
#: and values near 0 as indeterminate; the bands below follow that.
STRONGLY_NEGATIVE: Final = -1.0
MILDLY_NEGATIVE: Final = -0.25
MILDLY_POSITIVE: Final = 0.25


def load_model() -> dict[str, Any]:
    with resources.files("openbiota.data").joinpath("gmwi2_model.json").open("r") as handle:
        return json.load(handle)


@dataclass(frozen=True, slots=True)
class DysbiosisAnchor:
    score: float
    n_health_taxa_present: int
    n_disease_taxa_present: int
    contributions: tuple[tuple[str, float], ...]
    database: str
    valid: bool
    note: str = ""

    @property
    def band(self) -> str:
        if self.score <= STRONGLY_NEGATIVE:
            return "strongly dysbiotic"
        if self.score <= MILDLY_NEGATIVE:
            return "mildly dysbiotic"
        if self.score < MILDLY_POSITIVE:
            return "indeterminate"
        return "not dysbiotic"

    @property
    def strongly_negative(self) -> bool:
        return self.score <= STRONGLY_NEGATIVE

    @property
    def caution(self) -> str:
        """The sentence the report must print beside every profile score."""
        if not self.valid:
            return (
                "The dysbiosis anchor could not be computed, so profile scores cannot be "
                "checked against general gut disturbance."
            )
        if self.strongly_negative:
            return (
                f"The general dysbiosis index is strongly negative ({self.score:+.2f}). This "
                "community is broadly disturbed, so a high score on any profile most likely "
                "reflects that general disturbance rather than anything condition-specific."
            )
        if self.score <= MILDLY_NEGATIVE:
            return (
                f"The general dysbiosis index is mildly negative ({self.score:+.2f}). Read "
                "profile scores with that in mind: part of any elevation may be general "
                "disturbance."
            )
        return (
            f"The general dysbiosis index is {self.score:+.2f} ({self.band}). Profile scores "
            "are not being driven by broad gut disturbance."
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "index": "GMWI2",
            "score": round(self.score, 4),
            "band": self.band,
            "valid": self.valid,
            "database": self.database,
            "health_associated_taxa_present": self.n_health_taxa_present,
            "disease_associated_taxa_present": self.n_disease_taxa_present,
            "top_contributions": [
                {"clade": clade, "coefficient": round(coef, 4)}
                for clade, coef in self.contributions[:12]
            ],
            "caution": self.caution,
            "note": self.note,
            "citation": (
                "Chang D., Gupta V.K. et al., Nat Commun 15:7447 (2024), "
                "doi:10.1038/s41467-024-51651-9"
            ),
        }


def compute_anchor(
    clade_abundance: dict[str, float], *, unknown_percent: float, database: str
) -> DysbiosisAnchor:
    """Apply the published GMWI2 model to a MetaPhlAn 3 profile.

    ``clade_abundance`` maps full clade lineage strings (as MetaPhlAn prints
    them) to relative abundance in percent. ``unknown_percent`` is the
    UNKNOWN row. ``database`` is the MetaPhlAn index the profile came from.
    """
    model = load_model()
    valid = any(database.startswith(d) for d in COMPATIBLE_DATABASES)
    if not valid:
        return DysbiosisAnchor(
            score=0.0,
            n_health_taxa_present=0,
            n_disease_taxa_present=0,
            contributions=(),
            database=database,
            valid=False,
            note=(
                f"GMWI2 was trained on {COMPATIBLE_DATABASES[0]}; profiles from {database} "
                "are not comparable and the index was not computed."
            ),
        )

    classified = max(1e-9, 100.0 - unknown_percent)
    cutoff = float(model["presence_cutoff"])
    coefficients: dict[str, float] = model["coefficients"]

    score = float(model["intercept"])
    contributions: list[tuple[str, float]] = []
    n_health = n_disease = 0
    for clade, coefficient in coefficients.items():
        value = clade_abundance.get(clade, 0.0) / classified
        if value > cutoff:
            score += coefficient
            contributions.append((clade, coefficient))
            if coefficient > 0:
                n_health += 1
            else:
                n_disease += 1

    contributions.sort(key=lambda kv: -abs(kv[1]))
    return DysbiosisAnchor(
        score=score,
        n_health_taxa_present=n_health,
        n_disease_taxa_present=n_disease,
        contributions=tuple(contributions),
        database=database,
        valid=True,
        note=(
            f"{n_health} health-associated and {n_disease} disease-associated model taxa "
            f"present out of {len(coefficients)} in the published model."
        ),
    )
