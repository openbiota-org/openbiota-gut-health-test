"""Transmission concerns, screening status and evidence-based exclusions.

Implements the rule registry of spec §8.2 over this pipeline's own qualified
pathogen and determinant evidence (v05), plus an optional screening manifest
of validated laboratory results.

The hard separations, which this module exists to enforce:

* A **confirmed, material-linked** hazard or a documented program criterion can
  exclude material (X01–X04). Nothing else can.
* An unresolved high-consequence signal stays visible as review pending and
  never becomes a negative screen (R01, V7-051).
* A family-level resistance gene with no resolved carrier is a review item, not
  an MDRO (R02, V7-038).
* A disease-pattern percentile is context and can never exclude anyone (C01,
  V7-041).
* Sets inherit the union of member concerns; nothing dilutes (§8.4, V7-048).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from .identity import content_id
from .inputs import MaterialData

#: Compiled from spec §8.2. `disposition` is what the rule may do, never more.
RULES: Final = {
    "X01": {
        "condition": "validated, material-linked confirmation of a qualifying transmissible enteric pathogen",
        "requires": "external validated assay result linked to this donation or lot",
        "disposition": "exclusion_confirmed",
        "sources": ["S01", "S02", "S03"],
    },
    "X02": {
        "condition": "confirmed qualifying ESBL-producing Enterobacteriaceae, CRE, VRE or MRSA",
        "requires": "organism-linked confirmation; a taxon plus an unlinked resistance fragment is insufficient",
        "disposition": "exclusion_confirmed",
        "sources": ["S01"],
    },
    "X03": {
        "condition": "confirmed specimen identity mix-up, wrong person or failed analytical integrity",
        "requires": "identity or integrity failure in this analytical record",
        "disposition": "exclusion_confirmed",
        "sources": [],
    },
    "X04": {
        "condition": "an explicit documented clinical-program criterion is met",
        "requires": "the criterion text, version and the evidence it was applied to",
        "disposition": "exclusion_confirmed",
        "sources": ["S04"],
    },
    "R01": {
        "condition": "credible but unconfirmed pathogen, pathotype or toxin signal",
        "requires": "sequence evidence in this material",
        "disposition": "review_pending",
        "sources": ["S01", "S02", "S03"],
    },
    "R02": {
        "condition": "resistance gene or mobile element without resolved carrier, intact allele or phenotype",
        "requires": "determinant evidence without host linkage",
        "disposition": "review_pending",
        "sources": ["S01"],
    },
    "R03": {
        "condition": "credible pks/clb island or bft/ETBF-related evidence without confirmed pathogen classification",
        "requires": "toxin or genotoxicity determinant evidence",
        "disposition": "review_pending",
        "sources": ["S07", "S08"],
    },
    "R04": {
        "condition": "screening records absent, stale under their policy, or not linked to the analysed donation",
        "requires": "absence of a linked screening record",
        "disposition": "review_pending",
        "sources": [],
    },
    "R05": {
        "condition": "shared-locus hit, cross-species ambiguity, low breadth, or absent controls",
        "requires": "analytical qualification already recorded by the pathogen lane",
        "disposition": "review_pending",
        "sources": [],
    },
    "C01": {
        "condition": "disease resemblance, community index, microbiome age or broad dysbiosis pattern",
        "requires": "nothing beyond the score itself",
        "disposition": "context_only",
        "sources": [],
    },
    "C02": {
        "condition": "confirmed common coloniser without a relevant pathogenic determinant or context",
        "requires": "organism evidence without determinant linkage",
        "disposition": "context_only",
        "sources": [],
    },
    "C03": {
        "condition": "noninfectious animal-transfer or mechanistic evidence",
        "requires": "the specific model and its scope",
        "disposition": "context_only",
        "sources": ["S06", "S09", "S10"],
    },
}

#: Human transmission and causal-mechanism evidence (spec §8.3). Kept as data
#: so a report can cite the exact basis of a concern's severity.
SOURCES: Final = {
    "S01": {
        "statement": "Two immunocompromised FMT recipients developed invasive ESBL-producing E. coli infection; one died. Retained donor material contained matching organisms.",
        "tier": "human_fmt_transmission",
        "url": "https://www.fda.gov/safety/medical-product-safety-information/fecal-microbiota-transplantation-safety-communication-risk-serious-adverse-reactions-due",
        "scope": "demonstrated donor-transmission risk; not a reason to treat all E. coli as pathogenic",
    },
    "S02": {
        "statement": "FDA described six EPEC/STEC infections and four hospitalisations after stool-bank FMT; two further deaths were investigated with unknown attribution.",
        "tier": "human_fmt_transmission",
        "url": "https://www.fda.gov/safety/medical-product-safety-information/fecal-microbiota-transplantation-safety-alert-risk-serious-adverse-events-likely-due-transmission",
        "scope": "suspected transmission and uncertain mortality attribution are not proven claims",
    },
    "S03": {
        "statement": "FDA SARS-CoV-2 and mpox alerts require additional protections for potential transmission.",
        "tier": "precautionary_policy",
        "url": "https://www.fda.gov/vaccines-blood-biologics/safety-availability-biologics/information-pertaining-additional-safety-protections-regarding-use-fecal-microbiota-transplantation-0",
        "scope": "DNA stool sequencing does not clear RNA-virus risk",
    },
    "S04": {
        "statement": "The 2024 BSG/HIS second-edition guideline updates donor screening and pathogen-transmission requirements.",
        "tier": "clinical_program_policy",
        "url": "https://pubmed.ncbi.nlm.nih.gov/38609760/",
        "scope": "use a current fully sourced policy; detailed thresholds were not verified in the spec audit",
    },
    "S07": {
        "statement": "pks-positive E. coli produced a mutational signature in human intestinal organoids that is also found in human cancers.",
        "tier": "mechanistic",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8142898/",
        "scope": "an isolated low-identity clbB hit does not establish an intact active island",
    },
    "S08": {
        "statement": "Enterotoxigenic B. fragilis caused colitis and tumours in susceptible mice; nontoxigenic B. fragilis did not.",
        "tier": "animal_transfer",
        "url": "https://www.nature.com/articles/nm.2015",
        "scope": "species-level B. fragilis is not ETBF",
    },
    "S06": {
        "statement": "Obesity-discordant human twin microbiota transferred adiposity to germ-free mice; diet modified the effect.",
        "tier": "animal_transfer",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3829625/",
        "scope": "causal in that model; not a human donor-risk probability",
    },
    "S09": {
        "statement": "Parkinson's-patient microbiota aggravated impairment in alpha-synuclein-overexpressing mice.",
        "tier": "animal_transfer",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5718049/",
        "scope": "host-dependent experimental transfer, not human FMT-transmitted Parkinson's",
    },
    "S10": {
        "statement": "MS-derived microbiota worsened induced autoimmune encephalomyelitis in germ-free mice.",
        "tier": "animal_transfer",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5635915/",
        "scope": "supports mechanistic research, not exclusion by an MS-like microbiome score",
    },
}

#: Plain-language names for the catalogue's confirmation routes, and what each
#: would actually be run on. A donor candidate is usually a banked, frozen
#: aliquot, so "get a lab test" has to say *which* test and *on what*: a
#: culture needs viable organisms and may fail on a frozen aliquot, while a
#: nucleic-acid test works on frozen or fixed material.
CONFIRMATION_TESTS: Final[dict[str, dict[str, str]]] = {
    "clinical_stool_NAAT": {
        "name": "stool nucleic-acid test (PCR panel)",
        "material": "works on a thawed aliquot of the banked sample, or on a fresh sample from "
        "the donor",
    },
    "clinical_stool_culture": {
        "name": "stool culture",
        "material": "needs viable organisms, so it should be run on a fresh sample from the "
        "donor; a frozen aliquot may not grow even when the organism was there",
    },
    "clinical_toxin_immunoassay": {
        "name": "toxin immunoassay (EIA)",
        "material": "detects the toxin protein in a thawed aliquot or a fresh sample",
    },
    "clinical_stool_antigen": {
        "name": "stool antigen test",
        "material": "runs on a thawed aliquot or a fresh sample",
    },
    "modified_acid_fast_stain": {
        "name": "modified acid-fast microscopy",
        "material": "runs on a preserved or thawed aliquot",
    },
    "clinician_selected_stool_NAAT_or_culture": {
        "name": "a stool NAAT or culture chosen by the clinician",
        "material": "NAAT on a thawed aliquot; culture needs a fresh sample",
    },
}


def confirmation_sentence(options: Sequence[str]) -> str:
    """"Needs a lab test" is useless without naming the test and the material."""
    seen = [CONFIRMATION_TESTS[o] for o in options if o in CONFIRMATION_TESTS]
    if not seen:
        return (
            "No validated confirmation route is recorded for this target, so nothing here can be "
            "settled by a test we can name."
        )
    first = seen[0]
    extra = (
        " Alternatives: " + ", ".join(x["name"] for x in seen[1:]) + "."
        if len(seen) > 1 else ""
    )
    return f"What would settle it: a {first['name']} — {first['material']}.{extra}"


#: Screening that stool DNA sequencing cannot cover at all (spec §8.4).
UNCOVERED_SCREENING: Final = (
    "blood-borne viruses (HIV, HBV, HCV, HTLV) and syphilis serology",
    "RNA enteric viruses (norovirus, rotavirus, sapovirus, astrovirus, SARS-CoV-2)",
    "culture-confirmed enteric bacterial pathogens and antimicrobial phenotype",
    "ova and parasite microscopy, and protozoal antigen or PCR panels",
    "Helicobacter pylori stool antigen",
    "medical history, examination, travel and exposure screening",
    "viability of any organism in the material actually administered",
)


@dataclass(slots=True)
class ConcernRecord:
    """One qualified hazard record (spec §5.4)."""

    concern_id: str
    candidate_id: str
    material_ids: list[str]
    observation_refs: list[str]
    concern_type: str
    feature_label: str
    identity_resolution: str
    analytical_status: str
    causal_evidence: str
    severity: str
    rule_id: str
    rule_source_ids: list[str]
    disposition: str
    scope: str
    reason: str
    tier: str = "context"
    next_evidence_needed: list[str] = field(default_factory=list)
    measurements: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "concern_id": self.concern_id,
            "candidate_id": self.candidate_id,
            "material_ids": list(self.material_ids),
            "observation_refs": list(self.observation_refs),
            "concern_type": self.concern_type,
            "feature_label": self.feature_label,
            "identity_resolution": self.identity_resolution,
            "analytical_status": self.analytical_status,
            "causal_evidence": self.causal_evidence,
            "severity": self.severity,
            "tier": self.tier,
            "rule_id": self.rule_id,
            "rule": RULES.get(self.rule_id, {}).get("condition"),
            "rule_source_ids": list(self.rule_source_ids),
            "disposition": self.disposition,
            "scope": self.scope,
            "reason": self.reason,
            "next_evidence_needed": list(self.next_evidence_needed),
            "measurements": dict(self.measurements),
        }


def _severity_for(target: dict[str, Any]) -> str:
    cls = str(target.get("interpretation_class") or "")
    if cls == "established_enteric":
        return "high_consequence"
    if cls in ("conditional_opportunist", "uncertain_enteric_role"):
        return "context_dependent"
    return "uncertain"


def evaluate(
    data: MaterialData,
    *,
    candidate_id: str,
    screening: dict[str, Any] | None = None,
) -> tuple[list[ConcernRecord], dict[str, Any]]:
    """Build every concern for one material, plus its screening record."""
    concerns: list[ConcernRecord] = []
    pth = data.results.get("pathogens") or {}
    mid = data.material.material_id

    def add(**kw: Any) -> None:
        cid = content_id("concern", {"c": candidate_id, "f": kw["feature_label"], "r": kw["rule_id"]})
        concerns.append(ConcernRecord(concern_id=cid, candidate_id=candidate_id, material_ids=[mid], **kw))

    for row in sorted(pth.get("results") or [], key=lambda r: str(r.get("target_id"))):
        status = str(row.get("display_status") or "")
        if status in ("not_detected", "not_assessed"):
            continue
        frags = row.get("unique_supporting_fragments") or 0
        ambig = row.get("ambiguous_fragments") or 0
        regions = row.get("informative_regions_supported") or 0
        name = str(row.get("display_name") or row.get("target_id"))
        single_region = regions <= 1 and frags > 0
        no_specific = frags == 0
        resolution = str(row.get("resolution") or "unresolved")
        # Severity follows the analytical strength, not the organism's name. A
        # row with no organism-specific fragment has nothing to be severe
        # about: the pathogen lane already calls it ambiguous because the
        # fragments are shared with other taxa. Promoting those to
        # high-consequence review would bury the rows that do have support.
        rule = "R05" if (single_region or no_specific or status == "ambiguous_signal") else "R01"
        if str(row.get("interpretation_class")) == "background_or_decoy":
            rule = "C02"
        measurements = {
            "unique_supporting_fragments": frags,
            "ambiguous_fragments": row.get("ambiguous_fragments"),
            "informative_regions_supported": regions,
            "reference_breadth_fraction": row.get("reference_breadth_fraction"),
            "median_alignment_identity": row.get("median_alignment_identity"),
            "normalized_fragments_per_million": row.get("normalized_fragments_per_million"),
            "display_status": status,
            "contamination_status": row.get("contamination_status"),
        }
        # Where the pathogen lane's resolution pass has ruled, its verdict wins:
        # the two reports must not disagree about the same organism. A species
        # the data cannot name, and carriage without the toxin gene, are not
        # high-consequence findings however alarming the organism's name is.
        report_tier = str(row.get("report_tier") or "")
        pathotype = str(row.get("pathotype_evidence") or "not_applicable")
        if report_tier in {"unresolved_complex", "pathotype_negative"} and frags == 0:
            # Nothing organism-specific was found. "DNA present (0 fragments
            # across 0 regions)" is a sentence that contradicts itself, and it
            # is the single most alarming-looking line in the report for a row
            # that carries no evidence at all. Fall through to the generic
            # no-specific-sequence branch, which says that plainly.
            report_tier = ""
        if report_tier == "unresolved_complex":
            complex_label = str(row.get("complex_label") or "an indistinguishable complex")
            add(
                observation_refs=[str(row.get("target_id"))],
                concern_type="technical",
                feature_label=name,
                identity_resolution="not_resolvable_within_complex",
                analytical_status="ambiguous",
                causal_evidence="association",
                severity="uncertain",
                tier="species_unresolved",
                rule_id="R05",
                rule_source_ids=[],
                disposition="review_pending",
                scope="sequenced_sample",
                reason=(
                    f"{name}: {frags} fragment(s) aligned, but this organism is part of "
                    f"{complex_label} and cannot be named from sequence. In stool such evidence "
                    "usually comes from the harmless relative that most healthy guts carry. A "
                    "clinical culture or NAAT is what resolves it."
                ),
                next_evidence_needed=[
                    "the discriminating virulence marker for a species-level call",
                    confirmation_sentence(list(row.get("confirmation_options") or [])),
                ],
                measurements={
                    "unique_supporting_fragments": frags,
                    "informative_regions_supported": regions,
                    "reference_breadth_fraction": row.get("reference_breadth_fraction"),
                    "species_resolution": row.get("species_resolution"),
                    "display_status": status,
                },
            )
            continue
        if report_tier == "pathotype_negative":
            add(
                observation_refs=[str(row.get("target_id"))],
                concern_type="toxin",
                feature_label=name,
                identity_resolution=resolution,
                analytical_status="supported_genomic_evidence",
                causal_evidence="association",
                severity="context_dependent",
                tier="carriage_toxin_negative",
                rule_id="C02",
                rule_source_ids=[],
                disposition="review_pending",
                scope="sequenced_sample",
                reason=(
                    f"{name}: organism DNA present ({frags} fragment(s) across {regions} "
                    "region(s)), and the genes that make it capable of causing disease were "
                    "looked for and not found. That is carriage, not a disease finding. "
                    + (str(row.get("carriage_statement") or ""))
                ).strip(),
                next_evidence_needed=[
                    confirmation_sentence(list(row.get("confirmation_options") or [])),
                ],
                measurements={
                    "unique_supporting_fragments": frags,
                    "informative_regions_supported": regions,
                    "pathotype_evidence": pathotype,
                    "display_status": status,
                },
            )
            continue
        if no_specific:
            reason = (
                f"{name}: no organism-specific fragment in this material. {ambig} fragment(s) are shared "
                f"with other taxa, which is why the pathogen lane reports '{status}'. There is nothing "
                "here to identify this organism, and nothing here that clears it either."
            )
            tier = "technical_ambiguity"
            severity = "uncertain"
        elif single_region:
            reason = (
                f"{name}: {frags} organism-specific fragments but only {regions} informative region. A "
                "single region cannot be repaired by fragment count — a shared or conserved locus, a "
                "database gap and a real organism look the same — so identity stays unresolved."
            )
            tier = (
                "high_consequence_unresolved"
                if _severity_for(row) == "high_consequence"
                else "technical_ambiguity"
            )
            severity = _severity_for(row)
        else:
            reason = (
                f"{name}: {frags} organism-specific fragments across {regions} informative regions, "
                f"reported as '{status}'. Sequence evidence is not an infection diagnosis, no validated "
                "assay confirmed it, and this run had no extraction or library controls."
            )
            counts = row.get("counts_as_pathogen")
            tier = (
                "high_consequence_supported"
                if _severity_for(row) == "high_consequence"
                and status == "supported_sequence"
                and counts is not False
                else ("high_consequence_unresolved" if _severity_for(row) == "high_consequence" else "context")
            )
            severity = _severity_for(row)
        add(
            observation_refs=[str(row.get("target_id"))],
            concern_type="infectious_agent",
            feature_label=name,
            identity_resolution=resolution,
            analytical_status=(
                "ambiguous" if no_specific
                else ("supported_genomic_evidence" if status == "supported_sequence" else "sequence_signal")
            ),
            causal_evidence="potential_transmission" if rule == "R01" else "association",
            severity=severity,
            tier=tier,
            rule_id=rule,
            rule_source_ids=list(RULES[rule]["sources"]),
            disposition=RULES[rule]["disposition"],
            scope="sequenced_sample",
            reason=reason,
            next_evidence_needed=[
                confirmation_sentence(list(row.get("confirmation_options") or [])),
                "extraction and library negative controls for this batch, which this run did not have",
            ]
            + (["reference breadth across more than one informative region"] if single_region else []),
            measurements=measurements,
        )

    for det in sorted(pth.get("determinants") or [], key=lambda d: str(d.get("determinant_id"))):
        status = str(det.get("determinant_status") or "")
        if status in ("not_detected", "not_assessed"):
            continue
        kind = str(det.get("kind") or "")
        name = str(det.get("display_name") or det.get("determinant_id"))
        is_toxin = kind in ("toxin", "virulence", "virulence_locus") or str(
            det.get("determinant_id")
        ).startswith(("toxin.", "virulence"))
        rule = "R03" if is_toxin else "R02"
        add(
            observation_refs=[str(det.get("determinant_id"))],
            concern_type="toxin" if is_toxin else "mdro",
            feature_label=name,
            identity_resolution=str(det.get("host_linkage") or "unlinked"),
            analytical_status=(
                "supported_genomic_evidence" if status == "supported_intact_sequence" else "sequence_signal"
            ),
            causal_evidence="mechanistic" if is_toxin else "association",
            severity="context_dependent",
            tier="toxin_or_resistance_review" if is_toxin else "technical_ambiguity",
            rule_id=rule,
            rule_source_ids=list(RULES[rule]["sources"]),
            disposition=RULES[rule]["disposition"],
            scope="sequenced_sample",
            reason=(
                f"{name}: {det.get('supporting_fragments')} fragments, host linkage "
                f"'{det.get('host_linkage')}', locus completeness '{det.get('locus_completeness')}'. "
                "Gene-family sequence evidence without a resolved carrier organism, an intact allele "
                "and a phenotype is a review item, not a resistant-organism or toxin-producer call."
            ),
            next_evidence_needed=[
                "carrier organism attribution for this determinant",
                "allele completeness and discriminating positions",
                "phenotypic susceptibility testing where clinically relevant",
            ],
            measurements={
                "supporting_fragments": det.get("supporting_fragments"),
                "determinant_status": status,
                "host_linkage": det.get("host_linkage"),
                "locus_completeness": det.get("locus_completeness"),
                "amr_class": det.get("amr_class"),
                "gene_family": det.get("gene_family"),
            },
        )

    # Disease-pattern context (C01): never an exclusion, always visible.
    for obs in data.observations:
        if obs.feature_kind == "profile" and (obs.percentile or 0) >= 90:
            add(
                observation_refs=[obs.observation_id],
                concern_type="disease_association",
                feature_label=str(obs.extra.get("label") or obs.source_id),
                identity_resolution="community_pattern",
                analytical_status="not_assessed",
                causal_evidence="association",
                severity="uncertain",
                tier="context",
                rule_id="C01",
                rule_source_ids=["S09", "S10"],
                disposition="context_only",
                scope="sequenced_sample",
                reason=(
                    f"resemblance to the published {obs.extra.get('label') or obs.source_id} pattern at the "
                    f"{obs.percentile:.0f}th percentile. Shared-feature resemblance is not a diagnosis, not a "
                    "transmissible entity and cannot exclude this material."
                ),
                measurements={"percentile": obs.percentile, "status": obs.extra.get("status")},
            )

    coverage = (pth.get("coverage") or {})
    record = {
        "material_id": mid,
        "candidate_id": candidate_id,
        "screening_status": "not_provided" if not screening else str(screening.get("status") or "incomplete"),
        "policy": (screening or {}).get("policy"),
        "records": list((screening or {}).get("records") or []),
        "assay_coverage": {
            "pathogen_targets_total": coverage.get("total_targets"),
            "pathogen_targets_assessed": coverage.get("assessed"),
            "pathogen_targets_not_assessed": coverage.get("not_assessed"),
            "determinants_assessed": coverage.get("determinants_assessed"),
            "determinants_not_assessed": coverage.get("determinants_not_assessed"),
            "laboratory_controls": "not_available_in_this_run",
        },
        "not_covered_by_this_assay": list(UNCOVERED_SCREENING),
        "statement": (
            "No confirmed exclusion was identified in the available evidence. That is not a negative "
            "comprehensive pathogen screen and not a clinical clearance."
        ),
    }
    if record["screening_status"] != "records_complete_under_named_policy":
        concerns.append(
            ConcernRecord(
                concern_id=content_id("concern", {"c": candidate_id, "r": "R04"}),
                candidate_id=candidate_id,
                material_ids=[mid],
                observation_refs=[],
                concern_type="clinical_criterion",
                feature_label="clinical screening records",
                identity_resolution="not_applicable",
                analytical_status="not_assessed",
                causal_evidence="precautionary_policy",
                severity="context_dependent",
                tier="screening_gap",
                rule_id="R04",
                rule_source_ids=["S04"],
                disposition="review_pending",
                scope="person_with_time_window",
                reason=(
                    "No validated screening record linked to this donation was supplied, so screening is "
                    "incomplete. Missing evidence is neither a positive hazard nor a clearance."
                ),
                next_evidence_needed=[
                    "screening results under a named, current policy, linked to the analysed donation",
                    *UNCOVERED_SCREENING,
                ],
            )
        )
    return concerns, record


def exclusions_from_screening(
    *, candidate_id: str, material_id: str, screening: dict[str, Any] | None
) -> list[ConcernRecord]:
    """Confirmed, material-linked exclusions — the only kind that can exclude.

    Each entry must name its rule, the validated finding and the material it
    is linked to. Nothing here can be triggered by a metagenomic flag.
    """
    out: list[ConcernRecord] = []
    for entry in (screening or {}).get("confirmed_exclusions") or []:
        rule = str(entry.get("rule_id") or "X01")
        if rule not in ("X01", "X02", "X03", "X04"):
            continue
        finding = str(entry.get("finding") or "a validated confirmed finding")
        out.append(
            ConcernRecord(
                concern_id=content_id("concern", {"c": candidate_id, "r": rule, "f": finding}),
                candidate_id=candidate_id,
                material_ids=[material_id],
                observation_refs=[str(x) for x in entry.get("evidence_ids") or []],
                concern_type=str(entry.get("concern_type") or "infectious_agent"),
                feature_label=finding,
                identity_resolution="confirmed",
                analytical_status="confirmed_validated_assay",
                causal_evidence=str(entry.get("causal_evidence") or "human_fmt_transmission"),
                severity="high_consequence",
                tier="confirmed_exclusion",
                rule_id=rule,
                rule_source_ids=list(RULES[rule]["sources"]),
                disposition="exclusion_confirmed",
                scope=str(entry.get("scope") or "donation_lot"),
                reason=(
                    f"This material is excluded from candidate sets because {finding}; rule {rule} applies "
                    f"({RULES[rule]['condition']}). Source: {entry.get('source') or 'supplied screening manifest'}"
                    f", dated {entry.get('date') or 'date not supplied'}."
                ),
                measurements={"laboratory": entry.get("laboratory"), "assay": entry.get("assay")},
            )
        )
    return out
