"""Donor suitability: is this material healthy, and would it bring disease?

The matcher used to rank on one thing — the share of the recipient's gaps a
material could supply — and the report said so plainly: "What is
deliberately not in the score: GMWI2, diversity, microbiome age and
disease-pattern percentiles."

That is the wrong objective for a transplant. Restoration measures what a
donor *could give*. It says nothing about whether the donor's own community
is healthy, or whether it would bring a disease pattern the recipient does
not have. On the live data it ranked a donor with a **negative** gut-health
index and a 90th-percentile colorectal-cancer pattern (against a recipient
at the 15th) into the top three, because that donor happened to cover gaps.

This module supplies the two missing axes and the gates that go with them.

Three ideas do the work.

**A donor's own gut health is a precondition, not a tie-break.** GMWI2 is a
published index trained on 8,069 labelled stool samples; negative means the
community reads as dysbiotic. Transplanting a dysbiotic community to treat
dysbiosis is self-defeating, so a negative index disqualifies a material
from being recommended rather than costing it a few rank positions.

**"New" means new to this recipient.** A donor at the 97th percentile for a
pattern the recipient is already at the 100th for is not introducing
anything. Only the gap above the recipient counts.

**Presence is not a dose.** The previous organism-level rule flagged any
disease-associated species the donor carried and the recipient lacked. That
ranked the healthiest donor in the pool last, over *Alistipes timonensis* at
0.002% and *Escherichia coli* at 0.112% — normal carriage in a majority of
healthy adults. An organism only counts when the donor carries **more of it
than healthy adults normally do**, judged against the same reference cohort
the rest of the report uses.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

#: GMWI2 below this reads as a dysbiotic community on the published index.
#: Zero is the index's own neutral point, not a threshold invented here.
HEALTH_FLOOR: Final = 0.0

#: Observed GMWI2 range used to put the index on a 0-100 scale for display
#: and weighting. The index is unbounded in principle; this window covers
#: the range seen across the reference literature and the local samples.
HEALTH_SCALE_LOW: Final = -3.0
HEALTH_SCALE_HIGH: Final = 3.0

#: A donor must exceed this percentile of the healthy-adult distribution for
#: an organism before carrying it counts as introducing a burden. Below it
#: the donor carries no more than a normal person does, and a transplant
#: that moves a normal amount of a normal organism has introduced nothing
#: abnormal. p90 rather than the median because the question is "is this
#: unusual", not "is this present".
ABNORMAL_ABUNDANCE_PERCENTILE: Final = 0.90

#: Weights for the suitability composite, within an eligibility tier. They
#: are transparent, they sum to 1, and they are printed in the report.
#: Restoration and donor health are weighted near-equally because the brief
#: is to maximise healthy transfer, which needs both: a healthy donor who
#: supplies nothing, and a donor who supplies everything but is unwell, are
#: both poor choices. Disease burden carries the smallest weight here only
#: because its severe cases are handled by the gates above, not by weight.
WEIGHT_RESTORATION: Final = 0.45
WEIGHT_HEALTH: Final = 0.40
WEIGHT_LOW_BURDEN: Final = 0.15


@dataclass(frozen=True, slots=True)
class OrganismVerdict:
    """One organism a donor would add, with the evidence for judging it."""

    species: str
    donor_percent: float
    recipient_percent: float
    #: Healthy-adult reference quantiles for this species, or None when the
    #: species is absent from the reference cohort.
    reference_p50: float | None
    reference_p90: float | None
    carried_by_percent_of_adults: float | None
    #: Profiles that mark this organism as raised in their condition.
    raised_in: tuple[str, ...]
    is_known_pathogen: bool
    pathogen_note: str

    @property
    def abnormal(self) -> bool:
        """Whether the donor carries more than healthy adults normally do."""
        if self.reference_p90 is None:
            # Not in the reference: cannot say it is abnormal, so it is not
            # treated as such. Reported with that stated.
            return False
        return self.donor_percent > self.reference_p90

    @property
    def verdict(self) -> str:
        if self.is_known_pathogen:
            return "known pathogen"
        if self.reference_p90 is None:
            return "no healthy-adult reference"
        if self.abnormal:
            return "above the healthy range"
        if self.reference_p50 is not None and self.donor_percent <= self.reference_p50:
            return "below the healthy median"
        return "within the healthy range"

    @property
    def plain(self) -> str:
        """One sentence a reader can act on."""
        name = self.species.replace("_", " ")
        if self.reference_p90 is None:
            return (
                f"{name} at {self.donor_percent:.3f}%. This species is not in the "
                "healthy-adult reference, so there is no normal range to compare "
                "against and no judgement is made."
            )
        carriers_note = (
            " Quantiles are among adults who carry it, not across everyone: for an "
            "uncommon organism most people are at zero, and a trace detection would "
            "otherwise read as above the normal range."
        )
        share = (
            f"{self.carried_by_percent_of_adults:.0f}% of healthy adults carry it"
            if self.carried_by_percent_of_adults is not None else "carriage unknown"
        )
        if self.abnormal:
            return (
                f"{name} at {self.donor_percent:.3f}%, above the healthy-adult 90th "
                f"percentile of {self.reference_p90:.3f}% ({share}). This is more "
                "than a normal amount and is the kind of finding worth weighing."
                + carriers_note
            )
        return (
            f"{name} at {self.donor_percent:.3f}%, within normal carriage "
            f"(healthy-adult median {self.reference_p50:.3f}%, 90th percentile "
            f"{self.reference_p90:.3f}%; {share}). A transplant moving a normal "
            "amount of a common organism introduces nothing abnormal."
            + carriers_note
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "donor_percent": round(self.donor_percent, 4),
            "recipient_percent": round(self.recipient_percent, 4),
            "reference_median_percent": (
                None if self.reference_p50 is None else round(self.reference_p50, 4)
            ),
            "reference_p90_percent": (
                None if self.reference_p90 is None else round(self.reference_p90, 4)
            ),
            "carried_by_percent_of_healthy_adults": (
                None if self.carried_by_percent_of_adults is None
                else round(self.carried_by_percent_of_adults, 1)
            ),
            "raised_in_conditions": list(self.raised_in),
            "is_known_pathogen": self.is_known_pathogen,
            "pathogen_note": self.pathogen_note,
            "abnormal": self.abnormal,
            "verdict": self.verdict,
            "plain": self.plain,
        }


@dataclass(frozen=True, slots=True)
class DonorHealth:
    """A donor's own community health, independent of the recipient."""

    person_id: str
    gmwi2: float | None
    shannon: float | None
    #: Transfer-evidenced patterns this donor sits high on, with the margin
    #: over the recipient. Only positive margins are kept.
    new_patterns: tuple[tuple[str, float, float], ...] = ()
    organisms: tuple[OrganismVerdict, ...] = ()

    @property
    def health_score(self) -> float:
        """GMWI2 on a 0-100 display scale. 50 is the index's neutral point."""
        if self.gmwi2 is None:
            return 50.0
        span = HEALTH_SCALE_HIGH - HEALTH_SCALE_LOW
        return max(0.0, min(100.0, (self.gmwi2 - HEALTH_SCALE_LOW) / span * 100.0))

    @property
    def passes_health_floor(self) -> bool:
        return self.gmwi2 is None or self.gmwi2 >= HEALTH_FLOOR

    @property
    def disease_burden(self) -> float:
        """Total percentile margin over the recipient across new patterns."""
        return sum(margin for _, _, margin in self.new_patterns)

    @property
    def abnormal_organisms(self) -> tuple[OrganismVerdict, ...]:
        return tuple(o for o in self.organisms if o.abnormal or o.is_known_pathogen)

    @property
    def gate_failures(self) -> tuple[str, ...]:
        """Every reason this donor cannot be a primary recommendation."""
        out: list[str] = []
        if not self.passes_health_floor:
            out.append(
                f"gut-health index {self.gmwi2:+.2f} is below zero: this donor's own "
                "community reads as dysbiotic on an index trained on 8,069 labelled "
                "stool samples"
            )
        for profile, donor_pct, margin in self.new_patterns:
            out.append(
                f"{profile} pattern at the {donor_pct:.0f}th percentile, "
                f"{margin:.0f} points above the recipient, for a condition where "
                "patient stool reproduced disease features in recipient animals"
            )
        for organism in self.abnormal_organisms:
            out.append(
                f"{organism.species.replace('_', ' ')} at {organism.donor_percent:.3f}%, "
                + (
                    "a known pathogen"
                    if organism.is_known_pathogen
                    else f"above the healthy-adult 90th percentile of "
                         f"{organism.reference_p90:.3f}%"
                )
            )
        return tuple(out)

    @property
    def eligible(self) -> bool:
        return not self.gate_failures

    def to_json(self) -> dict[str, Any]:
        return {
            "person_id": self.person_id,
            "gmwi2": self.gmwi2,
            "health_score": round(self.health_score, 1),
            "shannon_index": self.shannon,
            "passes_health_floor": self.passes_health_floor,
            "health_floor": HEALTH_FLOOR,
            "new_transfer_evidenced_patterns": [
                {"profile": p, "donor_percentile": round(d, 1),
                 "margin_over_recipient": round(m, 1)}
                for p, d, m in self.new_patterns
            ],
            "disease_burden": round(self.disease_burden, 1),
            "organisms_it_would_add": [o.to_json() for o in self.organisms],
            "abnormal_organism_count": len(self.abnormal_organisms),
            "eligible": self.eligible,
            "gate_failures": list(self.gate_failures),
        }


@lru_cache(maxsize=1)
def _reference_quantiles(cohort_path: str) -> dict[str, tuple[float, float, float]]:
    """``species -> (p50, p90, carriage%)`` from the healthy-adult cohort.

    The same cohort the rest of the report compares against, so "normal for
    a healthy adult" means one thing throughout.
    """
    path = Path(cohort_path)
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, tuple[float, float, float]] = {}
    for index, species in enumerate(data.get("taxa", [])):
        column = data["abundance"][index]
        n = len(column)
        if not n:
            continue
        carriers = sorted(v for v in column if v > 0)
        carried = len(carriers) / n * 100.0
        if not carriers:
            continue
        # Quantiles among the adults who carry the organism at all, not
        # across everyone. Across everyone, a species carried by 7% of
        # people has a 90th percentile of exactly zero, so any trace
        # detection "exceeds the healthy range" - which is how Alistipes
        # timonensis at 0.002% came to disqualify the healthiest donor in
        # the pool. The question worth asking is: among people who have
        # this organism, does this donor have an unusual amount?
        m = len(carriers)
        out[species] = (
            carriers[min(m - 1, int(0.50 * m))],
            carriers[min(m - 1, int(ABNORMAL_ABUNDANCE_PERCENTILE * m))],
            carried,
        )
    return out


def _species_abundance(material: Any) -> dict[str, float]:
    """Species abundances from the lane the matching actually runs on.

    The typed species observations are built from the MetaPhlAn 3 profile
    (`inputs._species`, namespace ``metaphlan3:``), and every coverage
    decision and every goal is expressed against that lane — the targets
    say so: "detected at any abundance in the same species lane".

    Reading `extended_catalogue` instead would read MetaPhlAn 4, and the
    two lanes disagree about real organisms. SAMPLE4 carries *Adlercreutzia
    equolifaciens* in the MetaPhlAn 4 lane at 0.07% and not at all in
    MetaPhlAn 3, so a page built on the wrong lane listed it as something
    SAMPLE4 would add while the coverage table on the same page listed it as
    the one gap no candidate could fill. Both cannot be true.
    """
    out: dict[str, float] = {}
    for obs in getattr(material, "observations", ()):
        if obs.feature_kind != "species" or obs.value is None:
            continue
        value = float(obs.value)
        if value > 0:
            # A species can appear from more than one source in the lane
            # (direct profile row, guild member list). Keep the larger.
            out[str(obs.source_id)] = max(out.get(str(obs.source_id), 0.0), value)
    return out


def _percentiles(material: Any) -> dict[str, float]:
    ranked = (material.results.get("profile_similarity") or {}).get("ranked") or []
    return {
        str(row.get("profile")): float(row.get("percentile") or 0.0)
        for row in ranked
        if row.get("percentile") is not None
    }


def _gmwi2(material: Any) -> float | None:
    anchor = (material.results.get("profile_similarity") or {}).get("dysbiosis_anchor")
    if not isinstance(anchor, dict):
        return None
    score = anchor.get("score")
    return None if score is None else float(score)


def _pathogen_names(material: Any) -> dict[str, str]:
    """Species this sample's pathogen screen actually called, with the call."""
    records = (material.results.get("pathogens") or {}).get("results") or []
    out: dict[str, str] = {}
    for row in records:
        name = str(row.get("species") or row.get("display_name") or "")
        tier = str(row.get("report_tier") or "")
        if name and tier in {"pathogen", "pathogen_candidate"}:
            out[name.replace(" ", "_")] = tier
    return out


def assess_donor(
    *,
    donor: Any,
    recipient: Any,
    person_id: str,
    evidence: Any,
    actionable_tiers: frozenset[str],
    adverse_organisms: Any,
    profile_set: Any,
    high_percentile: float,
    new_margin: float,
    cohort_path: str = "refs/taxonomic_cohort.json",
) -> DonorHealth:
    """Build one donor's health record against this recipient.

    Everything the gates need, computed once and carried as data so the
    report can show the evidence for each decision rather than restating a
    conclusion.
    """
    reference = _reference_quantiles(cohort_path)
    donor_abundance = _species_abundance(donor)
    recipient_abundance = _species_abundance(recipient)
    donor_pct = _percentiles(donor)
    recipient_pct = _percentiles(recipient)
    pathogens = _pathogen_names(donor)

    new_patterns: list[tuple[str, float, float]] = []
    raised_by_species: dict[str, set[str]] = {}
    for profile, ev in evidence.items():
        if ev.tier not in actionable_tiers:
            continue
        theirs = donor_pct.get(profile, 0.0)
        mine = recipient_pct.get(profile, 0.0)
        if theirs < high_percentile or theirs - mine < new_margin:
            continue
        new_patterns.append((profile, theirs, theirs - mine))
        if profile_set is not None:
            for species in adverse_organisms(profile, profile_set):
                raised_by_species.setdefault(species, set()).add(profile)

    organisms: list[OrganismVerdict] = []
    for species, profiles in sorted(raised_by_species.items()):
        theirs = donor_abundance.get(species, 0.0)
        mine = recipient_abundance.get(species, 0.0)
        # Only organisms the donor actually has and the recipient does not.
        if theirs <= 0.0 or mine > 0.0:
            continue
        p50, p90, carriage = reference.get(species, (None, None, None))
        organisms.append(
            OrganismVerdict(
                species=species,
                donor_percent=theirs,
                recipient_percent=mine,
                reference_p50=p50,
                reference_p90=p90,
                carried_by_percent_of_adults=carriage,
                raised_in=tuple(sorted(profiles)),
                is_known_pathogen=species in pathogens,
                pathogen_note=(
                    f"called by the pathogen screen as {pathogens[species]}"
                    if species in pathogens else ""
                ),
            )
        )

    community = donor.results.get("community_profile") or {}
    return DonorHealth(
        person_id=person_id,
        gmwi2=_gmwi2(donor),
        shannon=community.get("shannon_index"),
        new_patterns=tuple(sorted(new_patterns, key=lambda t: -t[2])),
        organisms=tuple(organisms),
    )


def combine(person_health: dict[str, DonorHealth], members: list[str]) -> DonorHealth:
    """A set is only as suitable as its least suitable member.

    Health is the minimum, not the mean: pooling a dysbiotic material with
    a healthy one does not produce a healthy material, and averaging would
    let a good donor launder a bad one. Patterns and organisms are unions,
    for the same reason.
    """
    members_health = [person_health[m] for m in members if m in person_health]
    if not members_health:
        return DonorHealth(person_id="+".join(members), gmwi2=None, shannon=None)
    gmwis = [h.gmwi2 for h in members_health if h.gmwi2 is not None]
    seen: set[str] = set()
    organisms: list[OrganismVerdict] = []
    for h in members_health:
        for o in h.organisms:
            if o.species not in seen:
                seen.add(o.species)
                organisms.append(o)
    patterns: dict[str, tuple[str, float, float]] = {}
    for h in members_health:
        for profile, pct, margin in h.new_patterns:
            best = patterns.get(profile)
            if best is None or margin > best[2]:
                patterns[profile] = (profile, pct, margin)
    return DonorHealth(
        person_id="+".join(members),
        gmwi2=min(gmwis) if gmwis else None,
        shannon=min((h.shannon for h in members_health if h.shannon is not None),
                    default=None),
        new_patterns=tuple(sorted(patterns.values(), key=lambda t: -t[2])),
        organisms=tuple(sorted(organisms, key=lambda o: o.species)),
    )


def suitability_score(
    *, restoration: float, health: DonorHealth
) -> float:
    """The composite used to order candidates inside an eligibility tier.

    Deliberately a weighted sum of three quantities a reader can check, not
    a fitted model: there is no labelled outcome data to fit one against,
    and an opaque number here would be worse than an explicit arithmetic
    one. The weights are module constants and are printed in the report.
    """
    burden_term = max(0.0, 100.0 - min(100.0, health.disease_burden))
    return (
        WEIGHT_RESTORATION * restoration
        + WEIGHT_HEALTH * health.health_score
        + WEIGHT_LOW_BURDEN * burden_term
    )


def score_explanation(*, restoration: float, health: DonorHealth) -> str:
    """The arithmetic, spelled out, for the report."""
    burden_term = max(0.0, 100.0 - min(100.0, health.disease_burden))
    total = suitability_score(restoration=restoration, health=health)
    return (
        f"{WEIGHT_RESTORATION:.2f} x {restoration:.1f} restoration "
        f"+ {WEIGHT_HEALTH:.2f} x {health.health_score:.1f} donor gut health "
        f"+ {WEIGHT_LOW_BURDEN:.2f} x {burden_term:.1f} freedom from new disease "
        f"patterns = {total:.1f}"
    )
