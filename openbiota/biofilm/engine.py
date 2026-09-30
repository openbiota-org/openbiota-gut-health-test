"""Assembling one sample's biofilm result.

The engine's responsibility is to produce the typed result contract from
section 9 without ever inventing a value. Everything it emits is either a
measurement, a rank against a named reference, a registry fact, or an
explicit unavailability with a reason code.

Two rules shape most of the code here.

**The axes never interact.** ``harmful_associated`` and
``protective_associated`` are computed independently, from disjoint
features, against separately hashed calibrations. There is no subtraction,
no ratio, no net score and no combined "biofilm health" number anywhere in
this module, because the two things genuinely can be high at once and the
arithmetic that would hide that is the arithmetic the spec forbids.

**Unavailability is typed.** A missing module produces a reason code, not
a zero. :func:`_card` returns ``status`` alongside the percentile, and a
null percentile always carries the reason it is null - out of domain, too
few reference participants, a constant reference, or a feature that was
never assayed. The report renders those states as text, never as a green
zero.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from . import datasets, interventions, reference, registry, scoring
from . import sources as src

SCHEMA_VERSION: Final = "openbiota.biofilm/1.1"
SPEC_VERSION: Final = "08.0.1"
REGISTRY_RELEASE: Final = "biofilm-v08"

#: Headline label for a DNA-only result. "Biofilm state" is reserved for a
#: described multimodal observation and is never produced from sequence.
MEASUREMENT_LABEL: Final = "Biofilm-related genetic potential"

#: Percentile band names. Deliberately descriptive quartiles, not clinical
#: normal/abnormal thresholds.
def band(percentile: float | None) -> str:
    if percentile is None:
        return "not computed"
    if percentile < 25:
        return "lower reference range"
    if percentile <= 75:
        return "middle reference range"
    return "upper reference range"


@dataclass(frozen=True, slots=True)
class Card:
    """One of the four quantitative cards from section 4."""

    card: str
    heading: str
    label: str
    proxy_id: str
    status: str
    reference_percentile: float | None
    raw_index: float | None
    feature_percentiles: dict[str, float]
    feature_values: dict[str, float]
    feature_detail: dict[str, Any]
    reference_n: int | None
    calibration_id: str | None
    intervals: dict[str, Any]
    scope: str
    scope_names: tuple[str, ...]
    module_ids: tuple[str, ...]
    source_ids: tuple[str, ...]
    reason_codes: tuple[str, ...]
    omissions: tuple[str, ...]
    what_it_is_not: str
    direction_note: str
    transport_change: str
    coverage_supported: int
    coverage_assayed: int
    largest_driver: dict[str, Any] | None = None

    @property
    def is_proxy(self) -> bool:
        """Both current cards are proxies and must say so on their face."""
        return self.proxy_id.startswith("BF-PROXY")

    def to_json(self) -> dict[str, Any]:
        return {
            "card": self.card,
            "heading": self.heading,
            "label": self.label,
            "proxy_id": self.proxy_id,
            "status": self.status,
            "reference_percentile": (
                None
                if self.reference_percentile is None
                else round(self.reference_percentile, 2)
            ),
            "band": band(self.reference_percentile),
            "raw_index": None if self.raw_index is None else round(self.raw_index, 4),
            "feature_percentiles": {
                k: round(v, 2) for k, v in self.feature_percentiles.items()
            },
            "feature_values": {k: round(v, 6) for k, v in self.feature_values.items()},
            "feature_detail": self.feature_detail,
            "reference_n": self.reference_n,
            "calibration_id": self.calibration_id,
            "intervals": self.intervals,
            "scope": self.scope,
            "scope_names": list(self.scope_names),
            "module_ids": list(self.module_ids),
            "source_ids": list(self.source_ids),
            "reason_codes": list(self.reason_codes),
            "omissions": list(self.omissions),
            "what_it_is_not": self.what_it_is_not,
            "direction_note": self.direction_note,
            "transport_change": self.transport_change,
            "coverage": {
                "supported": self.coverage_supported,
                "assayed": self.coverage_assayed,
                "unit": "modules, not a health score",
            },
            "largest_driver": self.largest_driver,
            "research_index_note": (
                "research index; not measured biofilm quantity"
                if self.is_proxy
                else ""
            ),
        }


def _regroup(
    groups: tuple[tuple[str, ...], ...], keep: tuple[str, ...]
) -> tuple[tuple[str, ...], ...]:
    """Restrict a frozen dependence structure to the features still present.

    Groups keep their identity: dropping one genus from a four-genus group
    leaves one group of three, not three groups of one. That matters,
    because splitting a correlated group into singletons would give each
    surviving feature its own group weight and change the arithmetic in
    the sample's favour or against it arbitrarily. A group that loses all
    its members disappears entirely.
    """
    kept = set(keep)
    return tuple(
        tuple(name for name in group if name in kept)
        for group in groups
        if any(name in kept for name in group)
    )


def _card_for_proxy(
    proxy: scoring.Proxy,
    features: dict[str, float],
    detail: dict[str, Any],
    cohort: reference.ReferenceCohort | None,
    *,
    out_of_domain: str,
    input_hashes: tuple[str, ...],
) -> Card:
    """Rank one proxy, or explain precisely why it could not be ranked.

    When the reference cannot measure one of the proxy's features, the
    named full proxy is not computed. Section 8.7 is explicit that all the
    named features are required, so the card falls back to a partial mask
    with its own calibration, its own label, and the dropped feature
    recorded as an omission - it does not quietly average over whatever
    happens to be available.
    """
    usable = (
        proxy.features
        if cohort is None
        else cohort.usable_features(proxy.features)
    )
    dropped = tuple(f for f in proxy.features if f not in usable)
    partial = bool(dropped)
    # The whole proxy rests on features the reference cannot rank.
    if not usable:
        partial = False
        usable = proxy.features

    groups = _regroup(proxy.groups, usable)
    label = proxy.label if not partial else f"{proxy.label} (partial panel)"
    proxy_id = proxy.id if not partial else f"{proxy.id}-partial{len(usable)}"

    sub_features = {name: features[name] for name in usable}
    sub_detail = {name: detail[name] for name in usable}
    base = {
        "card": proxy.card,
        "heading": proxy.heading,
        "label": label,
        "proxy_id": proxy_id,
        "feature_values": sub_features,
        "feature_detail": sub_detail,
        "scope": (
            "single_mechanism" if len(groups) == 1 else "multiple_mechanisms"
        ),
        "scope_names": (label,),
        "module_ids": (proxy.module_id,),
        "source_ids": proxy.source_ids,
        "what_it_is_not": (
            proxy.what_it_is_not
            + (
                ""
                if not partial
                else " This is a partial panel: "
                + "; ".join(
                    f"{name} could not be measured in this reference"
                    for name in dropped
                )
                + ". It is calibrated separately and is not the named "
                f"{len(proxy.features)}-feature index."
            )
        ),
        "direction_note": proxy.direction_note,
        "transport_change": proxy.transport_change,
        "coverage_supported": sum(
            1 for name in usable if detail[name]["detected"]
        ),
        "coverage_assayed": len(usable),
        "raw_index": None,
        "feature_percentiles": {},
        "largest_driver": None,
    }

    if cohort is None:
        return Card(
            **base,
            status="reference_unavailable",
            reference_percentile=None,
            reference_n=None,
            calibration_id=None,
            intervals={"analytical": None, "reference": None},
            reason_codes=("reference_cohort_unavailable",),
            omissions=(),
        )

    if out_of_domain:
        # An adult-only calibration cannot rank a child. The measurement is
        # still shown; only the percentile is withheld.
        return Card(
            **base,
            status="out_of_domain",
            reference_percentile=None,
            reference_n=cohort.n,
            calibration_id=None,
            intervals={"analytical": None, "reference": None},
            reason_codes=("out_of_domain", out_of_domain),
            omissions=(),
        )

    if cohort.n < scoring.MIN_REFERENCE_N:
        return Card(
            **base,
            status="insufficient_reference",
            reference_percentile=None,
            reference_n=cohort.n,
            calibration_id=None,
            intervals={"analytical": None, "reference": None},
            reason_codes=(f"reference_n_below_{scoring.MIN_REFERENCE_N}",),
            omissions=(),
        )

    rows = [{name: row[name] for name in usable} for row in cohort.rows]
    calibration = scoring.calibration_hash(
        proxy, rows, f"{cohort.id}|{proxy_id}"
    )
    try:
        ranked = scoring.rank_panel(sub_features, rows, groups)
    except scoring.ReferenceNonDiscriminating:
        return Card(
            **base,
            status="reference_non_discriminating",
            reference_percentile=None,
            reference_n=cohort.n,
            calibration_id=calibration,
            intervals={"analytical": None, "reference": None},
            reason_codes=("reference_non_discriminating",),
            omissions=(),
        )
    except (ValueError, scoring.MaskMismatch) as exc:
        return Card(
            **base,
            status="not_computable",
            reference_percentile=None,
            reference_n=cohort.n,
            calibration_id=calibration,
            intervals={"analytical": None, "reference": None},
            reason_codes=(type(exc).__name__,),
            omissions=(str(exc),),
        )

    feature_percentiles = {
        k: float(v) for k, v in ranked["feature_percentiles"].items()
    }
    # The single largest driver, named plainly, per section 8.8.
    driver_name = max(feature_percentiles, key=lambda k: feature_percentiles[k])
    driver = {
        "feature": driver_name,
        "percentile": round(feature_percentiles[driver_name], 2),
        "value": round(sub_features[driver_name], 6),
        "unit": "relative abundance (% of classified)",
        "detected": sub_detail[driver_name]["detected"],
        "note": (
            f"{driver_name} ranks highest of this card's "
            f"{len(usable)} features against the reference."
        ),
    }

    ref_interval = scoring.reference_interval(
        sub_features,
        rows,
        groups,
        input_hashes=input_hashes,
        registry_release=REGISTRY_RELEASE,
        calibration_hash=calibration,
    )
    # Fragment-level resampling needs the reads. These proxies are computed
    # from a summary taxonomic table, so the analytical interval is
    # unavailable rather than invented.
    ana_interval = scoring.analytical_interval(
        available=False, reason="taxonomic_summary_table_only"
    )

    undetected = tuple(
        [
            f"{name}: not detected above this assay's limit"
            for name in usable
            if not sub_detail[name]["detected"]
        ]
        + [
            f"{name}: dropped from the panel - "
            + (cohort.unmeasurable.get(name, "not measurable in this reference"))
            for name in dropped
        ]
    )

    return Card(
        **{**base, "raw_index": float(ranked["raw_index"]),
           "feature_percentiles": feature_percentiles,
           "largest_driver": driver},
        status="available",
        reference_percentile=float(ranked["reference_percentile"]),
        reference_n=cohort.n,
        calibration_id=calibration,
        intervals={
            "analytical": ana_interval.to_json(),
            "reference": ref_interval.to_json(),
        },
        reason_codes=(),
        omissions=undetected,
    )


def _mechanism_card(
    *,
    card: str,
    heading: str,
    label: str,
    axis: str,
) -> Card:
    """H-M and P-M, which have no validated sequence panel yet.

    Both are reported as unavailable with the specific reason, because the
    alternative - filling them from the proxies - is exactly what BF-T098
    forbids. A protective ecological signal does not become a protective
    mechanism measurement by being the only number available.
    """
    axis_modules = tuple(
        m.id
        for m in registry.MODULES.values()
        if m.candidate_axis == axis and m.id not in {"BF-M18", "BF-M19"}
    )
    blocked = tuple(
        f"{m}: {registry.get(m).ineligibility_reason}" for m in axis_modules
    )
    return Card(
        card=card,
        heading=heading,
        label=label,
        proxy_id="",
        status="no_validated_panel",
        reference_percentile=None,
        raw_index=None,
        feature_percentiles={},
        feature_values={},
        feature_detail={},
        reference_n=None,
        calibration_id=None,
        intervals={"analytical": None, "reference": None},
        scope="none",
        scope_names=(),
        module_ids=axis_modules,
        source_ids=(),
        reason_codes=("no_module_eligible_for_axis",),
        omissions=blocked,
        what_it_is_not=(
            "This card is blank because no mechanism module has a validated "
            "sequence panel yet, not because nothing was found. The community "
            "card beside it is a different measurement and is not a substitute."
        ),
        direction_note="",
        transport_change="",
        coverage_supported=0,
        coverage_assayed=len(axis_modules),
    )


def _domain_state(results_json: dict[str, Any]) -> tuple[str, str]:
    """Whether the adult reference applies to this subject.

    Returns ``(out_of_domain_reason, assumption_caveat)``. Exactly one is
    ever non-empty.

    Age is read from the existing subject-context record, which already
    distinguishes supplied facts from assumed ones, and is never inferred
    from a file name. Three cases:

    * Age supplied and under 18, or ``adult`` supplied as false: out of
      domain. The measurement is shown and the percentile is withheld,
      because an adult reference tail is not pediatric abnormality.
    * ``adult`` merely assumed: the percentile is computed, but the
      assumption travels with it as a caveat rather than disappearing.
    * Age supplied and adult: no caveat.
    """
    ctx = results_json.get("subject_context") or {}
    supplied = ctx.get("supplied") or {}

    age = supplied.get("age", supplied.get("age_years"))
    if age is not None:
        try:
            if float(age) < reference.ADULT_MIN_AGE:
                return (
                    f"subject age {float(age):g} is below the adult reference "
                    f"domain ({reference.ADULT_MIN_AGE}+). An adult reference "
                    "tail is not pediatric abnormality.",
                    "",
                )
            return "", ""
        except (TypeError, ValueError):
            pass

    if supplied.get("adult") is False:
        return (
            "the subject is recorded as under 18. The adult reference cannot "
            "rank a child, and an adult tail is not pediatric abnormality.",
            "",
        )

    assumed_fields = {a.get("field") for a in (ctx.get("assumed") or [])}
    unknown_fields = {u.get("field") for u in (ctx.get("unknown") or [])}
    if "adult" in assumed_fields or "age" in unknown_fields:
        return (
            "",
            "Age was not supplied. This percentile assumes the subject is an "
            "adult, because the reference is adults aged "
            f"{reference.ADULT_MIN_AGE} and over. If the subject is a child "
            "the percentile does not apply and should be disregarded.",
        )
    return "", ""


def _rank_findings(
    cards: tuple[Card, ...]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Two ranked lists, per section 8.8.

    Returns ``(needing_review, supportive)``. They are separate because the
    two axes read in opposite directions: a high concerning percentile is
    what draws attention, and a high protective percentile is not. Merging
    them into one list sorted by percentile would rank a favourable 74
    above a concerning 8 and imply the wrong thing about both.

    Within each list the key is evidence tier, then human-gut
    applicability, then the measured rank in that list's own direction,
    then the stable ID. It is an order for navigating evidence and says so
    on every row.
    """
    review: list[dict[str, Any]] = []
    supportive: list[dict[str, Any]] = []

    for card in cards:
        if card.status != "available" or card.reference_percentile is None:
            continue
        module = registry.get(card.module_ids[0])
        sources = src.resolve_all(card.source_ids)
        human_gut = any(s.human_gut_applicable for s in sources)
        concerning = card.heading.startswith("Concerning")
        row = {
            "id": card.proxy_id,
            "card": card.card,
            "label": card.label,
            "heading": card.heading,
            "list": "needing_review" if concerning else "supportive",
            "percentile": round(card.reference_percentile, 2),
            "band": band(card.reference_percentile),
            "evidence_class": module.biological_evidence_status,
            "human_gut_applicable": human_gut,
            "largest_driver": card.largest_driver,
            "direction_note": (
                "Higher means more of this concerning pattern."
                if concerning
                else "Higher means more of this supportive pattern."
            ),
            "why_ranked_here": (
                "Placed by evidence tier, then by rank within this list. It is "
                "an order for reviewing evidence, not a severity score and not "
                "a risk probability."
            ),
            "source_ids": list(card.source_ids),
            "what_would_resolve": (
                "A carrier-resolved mechanism measurement, or an independently "
                "measured phenotype, would say more than this community proxy "
                "can."
            ),
            # Community-proxy tier is item 4 of the section 8.8 key; nothing
            # here reaches the measured-phenotype or same-carrier tiers.
            "_sort": (
                3,
                0 if human_gut else 1,
                # Each list sorts by what matters in its own direction: the
                # most concerning first, and the strongest support first.
                -card.reference_percentile,
                card.proxy_id,
            ),
        }
        (review if concerning else supportive).append(row)

    for group in (review, supportive):
        group.sort(key=lambda f: f["_sort"])
        for i, f in enumerate(group, 1):
            f.pop("_sort")
            f["rank"] = i
    return review, supportive


#: The observation types that would describe an actual biofilm rather than
#: a genetic potential. None is derivable from sequence, so each row is
#: either an externally supplied measurement or an explicit blank.
ACTUAL_STATE_ROWS: Final = (
    (
        "Attached biofilm seen directly",
        "endoscopy or histology",
        "Would show whether a biofilm is physically present, and where.",
    ),
    (
        "Community behaviour",
        "ex-vivo dispersal, viable burden, epithelial effect",
        "Would show what the community does, not just what it could do.",
    ),
    (
        "Expression or activity",
        "RNA or protein from a matched specimen",
        "Would show which of these capabilities are switched on.",
    ),
    (
        "Host response",
        "measured host protein, permeability or inflammation",
        "Would show whether the host is actually affected.",
    ),
)


def _actual_state(external: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The actual-state strip. Never auto-filled from DNA.

    Every row is ``not_supplied`` unless an external assay result was
    actually provided, because none of these can be inferred from sequence.
    The strip exists so that the absence is visible rather than implied.
    """
    supplied = {str(a.get("analyte", "")).lower() for a in external}
    return [
        {
            "observation": label,
            "would_come_from": source,
            "why_it_matters": why,
            "state": "supplied" if label.lower() in supplied else "not_supplied",
        }
        for label, source, why in ACTUAL_STATE_ROWS
    ]


def _agreement_table(cards: tuple[Card, ...]) -> list[dict[str, str]]:
    """The agreement table from section 8.8.

    Deliberately not a combined confidence number. Each row states what it
    measured about the named concern, and rows that were not measured say
    so rather than being filled from a correlated proxy.
    """
    concern = "biofilm-associated community pattern (BF-S01 association)"
    hc = next((c for c in cards if c.card == "H-C"), None)
    hc_state = "not_measured"
    if hc is not None and hc.status == "available" and hc.reference_percentile is not None:
        hc_state = (
            "supports_named_concern"
            if hc.reference_percentile > 75
            else "opposes_named_concern"
            if hc.reference_percentile < 25
            else "mixed"
        )
    return [
        {
            "row": "Molecular mechanism",
            "state": "not_measured",
            "detail": (
                "No mechanism module has a validated sequence panel, so no "
                "molecular mechanism measurement exists to agree or disagree."
            ),
        },
        {
            "row": "Community proxy",
            "state": hc_state,
            "detail": (
                "The two-feature H-C pattern, ranked against the adult control "
                "reference."
            ),
        },
        {
            "row": "Expression or activity",
            "state": "not_measured",
            "detail": "No RNA or protein assay was supplied.",
        },
        {
            "row": "Directly measured architecture",
            "state": "not_measured",
            "detail": (
                "No endoscopy or imaging was supplied. Generic inflammation "
                "cannot fill this row."
            ),
        },
        {
            "row": "Measured host consequence",
            "state": "not_measured",
            "detail": "No host protein, permeability or injury assay was supplied.",
        },
        {
            "row": "_concern",
            "state": concern,
            "detail": (
                "Two scores computed from the same reads or taxa are not "
                "independent confirmations, and no majority vote across "
                "correlated proxies becomes proof."
            ),
        },
    ]


def analyze(
    sample_dir: Path,
    results_json: dict[str, Any],
    *,
    cohort_path: str = "refs/taxonomic_cohort.json",
) -> dict[str, Any]:
    """Produce the full typed biofilm result for one sample."""
    sample_id = str(results_json.get("sample") or sample_dir.name)

    # Inputs. The hash covers the profile actually read, so a changed
    # profile changes the bootstrap seed and the cache key.
    profile_path = sample_dir / "taxonomy" / reference.MPA3_PROFILE
    input_hashes: tuple[str, ...] = ()
    if profile_path.exists():
        input_hashes = (
            hashlib.sha256(profile_path.read_bytes()).hexdigest()[:16],
        )

    try:
        cohort: reference.ReferenceCohort | None = reference.load_reference(
            cohort_path
        )
    except reference.ReferenceUnavailable:
        cohort = None

    profile_status = "available"
    try:
        features, detail = reference.sample_features(sample_dir)
    except reference.IncompatibleProfile as exc:
        profile_status = str(exc)
        features = dict.fromkeys(reference.FEATURE_LEAVES, 0.0)
        detail = {
            name: {
                "value": None,
                "unit": "relative abundance (%), MetaPhlAn 3 species level",
                "leaves_summed": list(leaves),
                "leaves_detected": {},
                "n_leaves_detected": 0,
                "detected": False,
                "censoring": "not_assayed",
                "transport_note": reference.TRANSPORT_NOTES[name],
            }
            for name, leaves in reference.FEATURE_LEAVES.items()
        }
        cohort = None

    out_of_domain, domain_caveat = _domain_state(results_json)

    hc = _card_for_proxy(
        scoring.COMMUNITY_PROXY,
        features,
        detail,
        cohort,
        out_of_domain=out_of_domain,
        input_hashes=input_hashes,
    )
    pe = _card_for_proxy(
        scoring.ECOLOGY_PROXY,
        features,
        detail,
        cohort,
        out_of_domain=out_of_domain,
        input_hashes=input_hashes,
    )
    hm = _mechanism_card(
        card="H-M",
        heading="Concerning biofilm potential",
        label="Harm-associated mechanisms",
        axis="harmful_associated",
    )
    pm = _mechanism_card(
        card="P-M",
        heading="Protective biofilm/ecosystem support",
        label="Protective-associated mechanisms",
        axis="protective_associated",
    )
    cards = (hm, hc, pm, pe)

    needing_review, supportive = _rank_findings(cards)

    # Action candidates, routed from the findings that actually exist.
    tags: list[str] = []
    if hc.status == "available" and (hc.reference_percentile or 0) > 75:
        tags.append("high_hc_proxy")
    if pe.status == "available" and (pe.reference_percentile or 100) < 25:
        tags.append("low_pe_proxy")
    tags.extend(["protective_support"])
    matched = interventions.retrieve(tuple(tags))

    module_results = [m.manifest() for m in registry.MODULES.values()]

    return {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "spec_version": SPEC_VERSION,
        "input_assays": ["stool_shotgun_dna"],
        "registry_release": REGISTRY_RELEASE,
        "measurement_label": MEASUREMENT_LABEL,
        "headline_note": (
            "DNA measures potential. Activity, attachment, thickness and "
            "location are not measured by this test."
        ),
        "profile_status": profile_status,
        "summary": {
            "measurement_label": MEASUREMENT_LABEL,
            "harmful_associated": {
                "status": hm.status,
                "reference_percentile": None,
                "intervals": {"analytical": None, "reference": None},
                "calibration": None,
                "omissions": list(hm.omissions),
                "module_ids": list(hm.module_ids),
                "reason_codes": list(hm.reason_codes),
            },
            "protective_associated": {
                "status": pm.status,
                "reference_percentile": None,
                "intervals": {"analytical": None, "reference": None},
                "calibration": None,
                "omissions": list(pm.omissions),
                "module_ids": list(pm.module_ids),
                "reason_codes": list(pm.reason_codes),
            },
            "community_pattern_proxy": {
                "status": hc.status,
                "reference_percentile": (
                    None
                    if hc.reference_percentile is None
                    else round(hc.reference_percentile, 2)
                ),
            },
            "protective_ecology_proxy": {
                "status": pe.status,
                "reference_percentile": (
                    None
                    if pe.reference_percentile is None
                    else round(pe.reference_percentile, 2)
                ),
            },
            "active_biofilm_burden": None,
            "anatomical_location": None,
            "clinical_disease_probability": None,
            "no_cancellation_note": (
                "The concerning and protective readings are independent. They "
                "are never subtracted, averaged or offset against each other, "
                "because both can genuinely be high at once."
            ),
        },
        "cards": [c.to_json() for c in cards],
        # Two lists, per section 8.8. ``ranked_findings`` is retained as the
        # concatenation for consumers that want one sequence, but the two
        # directions are never sorted against each other.
        "findings_needing_review": needing_review,
        "supportive_observations": supportive,
        "ranked_findings": [*needing_review, *supportive],
        "actual_state": _actual_state([]),
        "agreement_table": _agreement_table(cards),
        "module_results": module_results,
        "module_census": registry.census(),
        "intervention_candidates": [
            {
                "id": i.id,
                "compound": i.compound,
                "section": i.section,
                "endpoint": i.endpoint,
                "endpoint_plain": i.endpoint_plain,
                "direction": i.direction,
                "review_label": i.review_label,
                "model": i.model,
                "organisms": list(i.organisms),
                "strain": i.strain,
                "result": i.result,
                "lab_exposure": i.lab_exposure,
                "delivery": i.delivery,
                "delivery_note": i.delivery_note,
                "opposes": list(i.opposes),
                "caution": i.caution,
                "collateral": i.collateral,
                "funding_conflict": i.funding_conflict,
                "safety": i.safety,
                "doi": i.doi,
                "url": i.url,
            }
            for i in matched
        ],
        "intervention_census": interventions.census(),
        "dataset_audit": datasets.audit(),
        "source_count": len(src.SOURCES),
        "reference": (cohort.provenance() if cohort is not None else None),
        "out_of_domain": out_of_domain,
        "domain_caveat": domain_caveat,
        "external_assay_results": [],
        "context_results": [],
        "provenance": {
            "input_sha256": list(input_hashes),
            "reference_manifest_sha256": (cohort.manifest_id if cohort else None),
            "tool_versions": {"profiler": "MetaPhlAn 3.1 (reference-compatible lane)"},
            "parameters": {
                "adult_min_age": reference.ADULT_MIN_AGE,
                "bootstrap_replicates": scoring.BOOTSTRAP_REPLICATES,
                "min_reference_n": scoring.MIN_REFERENCE_N,
            },
            "calibration_id": hc.calibration_id,
            "spec_version": SPEC_VERSION,
        },
        "limits": [
            *( [domain_caveat] if domain_caveat else [] ),
            "DNA shows genetic potential, not activity, attachment or location.",
            "Both computable cards are experimental research proxies, not "
            "measured biofilm quantity.",
            "No validated stool-DNA classifier of protective versus harmful "
            "biofilm exists; this module does not claim one.",
            "A negative or low result never means biofilm-free, infection-free "
            "or cleared for donation.",
        ],
    }
