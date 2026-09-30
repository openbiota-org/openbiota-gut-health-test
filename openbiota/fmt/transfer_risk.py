"""Disease-pattern resemblance: descriptive context, never a donor exclusion.

This module used to implement a causal bridge: a donor at or above the 75th
percentile of a disease profile, at least 20 points above the recipient, plus
any animal-transfer paper for that disease, was treated as introducing a
transferable causal pattern and was blocked from the recommendation.

**That rule is retired (spec §8.4).** It was wrong in a way worth recording,
because the failure is easy to repeat:

* A percentile is resemblance to a case-control *community pattern*. It is not
  a feature of the donor, so there is nothing there to transfer.
* A paper showing that patient stool produced disease features in animals says
  something about the organisms those patients carried. It says nothing about
  whether this donor carries them.
* In the supplied samples the rule blocked SAMPLE1 on a colorectal-cancer
  resemblance score, while SAMPLE1's actual measured colibactin evidence was a
  single `clbB` fragment — not an intact island by any standard. The rule was
  simultaneously alarming about the wrong thing and silent about the measured
  evidence.

What replaces it lives in `openbiota.resolution.mechanisms`: candidate
mechanisms counted only where the sequence evidence was **measured in that
donor**, reported per donor, never pooled or averaged.

What remains here is the resemblance itself, kept because it is informative
context and because the spec requires the historical score to stay visible
(§8.1, acceptance 37). It is reported with the recipient's own reading beside
it, and it no longer decides anything. Two guards apply:

* A resemblance percentile can never exclude a donor by itself (§10).
* A recipient's own high percentile can never suppress a donor's measured
  finding (acceptance 62) — that comparison belongs to the mechanism ledger,
  not to this one.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

#: A donor pattern at or above this percentile is "high" enough to discuss.
HIGH_PERCENTILE: Final = 75.0

#: And it is only *new* to this recipient if it exceeds the recipient's own
#: reading by this margin. A donor at the 80th where the recipient is at the
#: 78th introduces nothing.
NEW_MARGIN: Final = 20.0


@dataclass(frozen=True)
class Mechanism:
    """An organism the transfer evidence runs through.

    Some transfer experiments are not about a diagnosis at all — they are about
    an organism that happened to dominate the donors' guts. Maeda's arthritis
    result came from RA patients whose microbiota was *Prevotella
    copri*-dominated, and monocolonising mice with *P. copri* alone reproduced
    it. Where the evidence is organism-mediated like this and the candidate
    does not carry the organism, the published pathway is not present in that
    candidate, whatever a pattern-resemblance percentile says.

    So a mechanism is checked against the candidate's own data. If it is
    absent, the exposure is still reported — with the reason — but it stops
    counting against the candidate, because the thing the evidence is about is
    not there to transfer.
    """

    label: str
    #: Which check to run. `segatella_copri_complex` resolves all 13 clades
    #: from the extended SGB lane.
    check: str
    note: str


@dataclass(frozen=True)
class TransferEvidence:
    """What is published about transferring this condition's community."""

    #: `human_donor_transfer` — patient stool into an animal reproduced
    #: disease features against healthy-donor controls.
    #: `model_donor_transfer` — only animal-model donors were used.
    #: `association_only` — no transfer evidence; association or rescue only.
    tier: str
    statement: str
    boundary: str
    sources: tuple[str, ...]
    mechanism: Mechanism | None = None


_HUMAN = "human_donor_transfer"
_MODEL = "model_donor_transfer"
_ASSOC = "association_only"

#: Public aliases. The tier a caller most often needs to name is the
#: strongest one, and reaching for a leading-underscore module private to
#: do it invites a copy of the literal instead.
HUMAN_DONOR_TRANSFER: Final = _HUMAN
MODEL_DONOR_TRANSFER: Final = _MODEL
ASSOCIATION_ONLY: Final = _ASSOC

#: Per-condition transfer evidence. A profile absent from this table is
#: reported as `not_assessed` rather than guessed at.
EVIDENCE: Final[Mapping[str, TransferEvidence]] = {
    "crc": TransferEvidence(
        _HUMAN,
        "Stool from colorectal-cancer patients, gavaged into mice, promoted intestinal "
        "carcinogenesis, accelerated adenoma progression and produced oncogenic epigenetic "
        "signatures, each against healthy-donor controls.",
        "Every experiment used a predisposed, germ-free or carcinogen-sensitised animal. None "
        "shows spontaneous cancer from stool alone, and none is a human FMT outcome.",
        (
            "Wong 2017 Gastroenterology doi:10.1053/j.gastro.2017.08.022",
            "Li 2019 EBioMedicine doi:10.1016/j.ebiom.2019.09.021",
            "Sobhani 2019 PNAS doi:10.1073/pnas.1912129116",
            "Zackular 2013 mBio doi:10.1128/mbio.00692-13",
        ),
    ),
    "adenoma": TransferEvidence(
        _HUMAN,
        "Colorectal-cancer stool enhanced adenoma progression in Apc-mutant mice against "
        "healthy-donor controls.",
        "Genetically predisposed recipient; adenoma progression is not adenoma induction.",
        ("Li 2019 EBioMedicine doi:10.1016/j.ebiom.2019.09.021",),
    ),
    "mdd": TransferEvidence(
        _HUMAN,
        "Stool from people with major depression induced depression-like behaviour in germ-free "
        "mice and in microbiota-depleted rats, with the donors' metabolic shifts.",
        "Behavioural features in animals, not a human psychiatric diagnosis. One human-donor "
        "experiment found no worsening at all.",
        ("Zheng 2016 Mol Psychiatry doi:10.1038/mp.2016.44",
         "Kelly 2016 J Psychiatr Res doi:10.1016/j.jpsychires.2016.07.019",
         "Knudsen 2021 Sci Rep doi:10.1038/s41598-021-01248-9 (null result)"),
    ),
    "parkinsons": TransferEvidence(
        _HUMAN,
        "Parkinson's-patient microbiota aggravated motor impairment in alpha-synuclein-"
        "overexpressing mice.",
        "Host-dependent: the recipient was already genetically predisposed.",
        ("Sampson 2016 Cell PMC5718049",),
    ),
    "ms": TransferEvidence(
        _HUMAN,
        "MS-patient microbiota increased the incidence or severity of experimental autoimmune "
        "encephalomyelitis in germ-free and transgenic mice against healthy-twin controls.",
        "Strongly primed hosts, and both healthy and disease communities can trigger EAE in "
        "these models.",
        ("Cekanaviciute 2017 PNAS PMC5635915", "Berer 2017 PNAS (discordant twins)"),
    ),
    "ad_clinical": TransferEvidence(
        _HUMAN,
        "Alzheimer's-patient microbiota impaired memory and hippocampal neurogenesis in "
        "microbiota-depleted rodents.",
        "Cognition and neurogenesis, not amyloid or tau pathology, and not a human diagnosis.",
        ("Grabrucker 2023 Brain PMC; Fujii 2019 gnotobiotic AD transfer",),
    ),
    "ad_mci": TransferEvidence(
        _HUMAN,
        "Mild-cognitive-impairment donor stool impaired learning and cerebral glucose uptake in "
        "wild-type mice against healthy-donor controls.",
        "Cognitive and metabolic features only.",
        ("MCI FMT study, Disease_Evidence.md N06",),
    ),
    "ad_preclinical_amyloid": TransferEvidence(
        _MODEL,
        "Amyloid pathology was amplified by transgenic-donor flora in transgenic recipients; the "
        "human preclinical-amyloid community itself has not been transferred.",
        "Donor and recipient were both AD-transgenic. This is amplification, not induction.",
        ("APPPS1 germ-free transfer, Disease_Evidence.md N01",),
    ),
    "ibs": TransferEvidence(
        _HUMAN,
        "IBS-patient stool produced visceral hypersensitivity and altered transit in "
        "germ-free rodents against healthy-donor controls.",
        "Gut-function features, not the full syndrome.",
        ("Crouzet 2013; De Palma 2017 Sci Transl Med",),
    ),
    "ibs_d": TransferEvidence(
        _HUMAN,
        "IBS-D donor stool transferred faster transit, barrier defects and anxiety-like "
        "behaviour to germ-free mice.",
        "Features, not diagnosis; donor-dependent.",
        ("De Palma 2017 Sci Transl Med",),
    ),
    "ibs_c": TransferEvidence(
        _HUMAN,
        "Constipation-predominant donor stool slowed transit in recipient rodents.",
        "Motility feature only.",
        ("Disease_Evidence.md, intestinal disease section",),
    ),
    "ibs_m": TransferEvidence(
        _ASSOC,
        "No mixed-type-specific transfer experiment was identified; the subtype rests on "
        "association.",
        "Treat as association only.",
        ("Disease_Evidence.md, intestinal disease section",),
    ),
    "ibd": TransferEvidence(
        _HUMAN,
        "IBD-patient communities worsened colitis in susceptible or chemically challenged "
        "recipients.",
        "Full inflammatory disease generally required a susceptible or challenged recipient; one "
        "discordant-twin experiment produced the opposite result.",
        ("Disease_Evidence.md, intestinal disease section; Knudsen 2024 UC twin counterexample",),
    ),
    "crohns": TransferEvidence(
        _HUMAN,
        "Crohn's-associated communities increased inflammatory susceptibility in challenged "
        "recipients.",
        "Susceptible or challenged recipient required.",
        ("Disease_Evidence.md, intestinal disease section",),
    ),
    "uc": TransferEvidence(
        _HUMAN,
        "UC-associated communities increased colitis severity in challenged recipients.",
        "One monozygotic-twin experiment found the healthy co-twin's community produced *more* "
        "disease activity, which is a direct warning against reading these labels simply.",
        ("Knudsen 2024 Comp Med doi:10.30802/aalas-cm-23-000065",),
    ),
    "t2d": TransferEvidence(
        _HUMAN,
        "Type-2-diabetes donor stool impaired glucose and insulin tolerance in recipient mice.",
        "Glucose dysregulation is narrower than T2D, and one key study had no healthy-donor "
        "comparator.",
        ("Wang 2022 BMC Endocr Disord doi:10.1186/s12902-022-01155-8",),
    ),
    "t1d": TransferEvidence(
        _MODEL,
        "Glucose impairment transferred in model systems; patient-community transfer of "
        "autoimmune T1D was not established.",
        "Glucose impairment is not autoimmune diabetes.",
        ("Disease_Evidence.md, nutrition/metabolism section",),
    ),
    "obesity": TransferEvidence(
        _HUMAN,
        "Stool from the obese twin of discordant human twin pairs produced greater fat mass in "
        "germ-free mice than the lean twin's.",
        "Diet-dependent, and cohousing abolished the effect.",
        ("Ridaura 2013 Science PMC3829625",),
    ),
    "masld": TransferEvidence(
        _HUMAN,
        "Fatty-liver-associated communities transferred hepatic steatosis features.",
        "Liver injury, not the human diagnosis.",
        ("Disease_Evidence.md, liver and biliary section",),
    ),
    "cirrhosis": TransferEvidence(
        _HUMAN,
        "Cirrhosis-associated communities transferred injury and encephalopathy-related features.",
        "Distinguish liver injury from its neurological complications.",
        ("Disease_Evidence.md, liver and biliary section",),
    ),
    "hypertension": TransferEvidence(
        _HUMAN,
        "Hypertensive-donor stool raised blood pressure in recipient animals.",
        "Often required an independent vascular insult.",
        ("Disease_Evidence.md, cardiovascular/renal section",),
    ),
    "acvd": TransferEvidence(
        _HUMAN,
        "Atherosclerosis-associated communities worsened plaque and thrombosis potential in "
        "predisposed recipients.",
        "Predisposed recipients; worsening of an induced lesion.",
        ("Disease_Evidence.md, cardiovascular/renal section",),
    ),
    "ckd": TransferEvidence(
        _HUMAN,
        "CKD-associated communities worsened kidney injury markers in recipients.",
        "Worsened response to an independent kidney insult.",
        ("Disease_Evidence.md, cardiovascular/renal section",),
    ),
    "ra": TransferEvidence(
        _HUMAN,
        "Stool from early rheumatoid-arthritis patients whose microbiota was Prevotella "
        "copri-dominated raised intestinal Th17 counts in germ-free SKG mice, which then "
        "developed severe arthritis. P. copri monocolonisation alone reproduced it.",
        "Three conditions had to hold at once, and none of them describes a human recipient: "
        "the host carried the SKG ZAP-70 mutation that makes it arthritis-prone, it was "
        "germ-free, and arthritis appeared only after an injected fungal challenge "
        "(zymosan/curdlan) \u2014 germ-free SKG mice given the same stool and no challenge did not "
        "develop arthritis. The donors were also a P. copri-dominated subset of RA patients, so "
        "the effect is tied to that organism rather than to an RA diagnosis. The species-level "
        "association itself replicates poorly, and at clade resolution none of the 13 species of "
        "the S. copri complex was associated with any condition across 1,635 cases including RA. "
        "In humans, FMT has been used to treat RA rather than cause it, and long-term recipient "
        "follow-up has not attributed new-onset RA to a transplant.",
        (
            "Maeda 2016 Arthritis Rheumatol doi:10.1002/art.39783",
            "Scher 2013 eLife doi:10.7554/eLife.01202",
            "Manghi 2023 Cell Host Microbe doi:10.1016/j.chom.2023.09.013 (clade-level null)",
            "Ten-year FMT follow-up, Microorganisms 2021 doi:10.3390/microorganisms9030548",
        ),
        mechanism=Mechanism(
            label="the Segatella (Prevotella) copri complex",
            check="segatella_copri_complex",
            note=(
                "The arthritis transfer evidence runs through P. copri: the donors were the "
                "P. copri-dominated subset of RA patients, and monocolonising germ-free SKG mice "
                "with P. copri alone reproduced the arthritis. Where no member of the complex is "
                "detected \u2014 checked across all 13 clades, not the species label \u2014 the "
                "published pathway is absent from this candidate, and any remaining resemblance "
                "rests on accessory taxa from one cohort rather than on the organism the "
                "experiments were about."
            ),
        ),
    ),
    "ankylosing_spondylitis": TransferEvidence(
        _HUMAN,
        "Ankylosing-spondylitis-associated communities transferred inflammatory and bone "
        "components in susceptible hosts.",
        "Components, in genetically susceptible recipients.",
        ("Disease_Evidence.md, immune/eye section",),
    ),
    "sle": TransferEvidence(
        _HUMAN,
        "Lupus-associated communities transferred autoimmune components.",
        "Components, not the systemic disease.",
        ("Disease_Evidence.md, immune/eye section",),
    ),
    "csu": TransferEvidence(
        _HUMAN,
        "Chronic-spontaneous-urticaria-associated communities produced mast-cell responses in "
        "challenge models.",
        "Challenge model; CSU-like responses rather than the condition.",
        ("Disease_Evidence.md, skin/allergy/lung section",),
    ),
    "longcovid": TransferEvidence(
        _HUMAN,
        "Stool from people with Long COVID produced lung and intestinal inflammation and "
        "anxiety-like behaviour in recipient mice.",
        "The comparator was unexposed controls rather than recovered people, so Long-COVID-"
        "specific attribution remains incomplete.",
        ("Disease_Evidence.md LC02 doi:10.1007/s12602-026-11197-2",),
    ),
    "celiac": TransferEvidence(
        _ASSOC,
        "No coeliac transfer experiment was identified in the catalogue.",
        "Association only; coeliac disease requires gluten and HLA susceptibility.",
        ("Disease_Evidence.md (absent from the transfer inventory)",),
    ),
    "mecfs": TransferEvidence(
        _ASSOC,
        "Listed in the catalogue's explicit boundary section: association must not be relabelled "
        "transfer.",
        "No disease-community transfer experiment supports it.",
        ("Disease_Evidence.md, boundaries section",),
    ),
    "alopecia_areata": TransferEvidence(
        _ASSOC,
        "Boundary case: hair-regrowth case reports are rescue observations, not disease transfer.",
        "Association and rescue only.",
        ("Disease_Evidence.md, boundaries section",),
    ),
    "androgenetic_alopecia": TransferEvidence(
        _ASSOC,
        "Boundary case in the catalogue; no transfer evidence.",
        "Association only.",
        ("Disease_Evidence.md, boundaries section",),
    ),
    "symptomatic_dermographism": TransferEvidence(
        _ASSOC,
        "No transfer experiment identified.",
        "Association only.",
        ("Disease_Evidence.md (absent from the transfer inventory)",),
    ),
}

#: The evidence is organism-mediated and the organism is not in this candidate.
MECHANISM_ABSENT: Final = "mechanism_absent"

#: Tiers that disqualify a candidate from being recommended.
#:
#: History, because this line has been changed twice and the reasoning
#: matters. BUILD_SPEC_v07.0 §8.4 asked for the percentile bridge to be
#: retired on the grounds that a resemblance percentile is not itself a
#: transferable feature. That reasoning is sound as far as it goes, and it
#: was implemented by emptying this set. The consequence was that the
#: highest-scoring candidate could carry a colorectal-cancer pattern at the
#: 90th percentile against a recipient at the 15th and still be recommended,
#: with the tie against an equally-scoring candidate that carried no such
#: pattern decided by a count of unconfirmed trace findings.
#:
#: The recipient's standing instruction is the opposite of the spec's, and on
#: a question of what material goes into their own body their instruction
#: governs. The bridge is re-armed for the one tier where patient stool has
#: actually reproduced disease features in recipient animals against
#: healthy-donor controls.
#:
#: What this does and does not claim: it does not claim the percentile
#: predicts transfer. It applies the precautionary rule that where a
#: condition has that grade of evidence and a candidate sits far above the
#: recipient on its pattern, an equally-scoring candidate without it is
#: preferred. Association-only tiers do not disqualify; they are reported and
#: they cost rank, but they do not veto.
ACTIONABLE_TIERS: Final[frozenset[str]] = frozenset({_HUMAN})

#: Rank penalty per carried pattern, by evidence tier. Applied to the
#: ordering, never to `match_score` itself: the score keeps meaning "share of
#: reachable gaps supplied" and stays comparable across runs.
TIER_RANK_PENALTY: Final[dict[str, int]] = {
    _HUMAN: 1000,
    _MODEL: 10,
    _ASSOC: 1,
}

#: Organisms that drive a profile in the adverse direction, per profile.
#: Built from the profile registry at call time so it cannot drift from the
#: definitions the resemblance score itself uses.
def adverse_organisms(profile_name: str, profile_set: Any) -> frozenset[str]:
    """Species a profile marks as ENRICHED in its condition.

    Direction is the whole point. A disease profile is a pattern with two
    halves: organisms raised in the disease (``d > 0``) and organisms
    depleted by it (``d < 0``). Faecalibacterium prausnitzii sits in the
    Crohn's, UC and IBD profiles at ``d = -1`` because it is *lost* in
    those conditions.

    A rule that vetoed any donor bringing an organism from a disease
    profile would therefore veto every healthy donor for every condition,
    because restoring the depleted half is precisely what the transfer is
    for. Only the enriched half can be "introducing the disease pattern".
    """
    try:
        profile = profile_set.by_name(profile_name)
    except Exception:  # noqa: BLE001 - an absent profile contributes nothing
        return frozenset()
    return frozenset(
        feature.key
        for module in profile.modules.values()
        for feature in module.features
        if getattr(feature, "engine", None) == "metaphlan"
        and getattr(feature, "d", 0) > 0
    )


#: How hard association-only patterns push a candidate down the ranking.
#:
#: ``strict``  - any carried pattern outranks restoration score. The
#:               candidate carrying nothing comes first even if it restores
#:               less. This is the literal reading of "downscore any and all
#:               disease profiles".
#: ``weighted`` - a transfer-evidenced pattern still vetoes absolutely, but
#:               association-only patterns cost a tie-break rather than
#:               leading the sort, so a materially better restoration score
#:               is not given up to avoid a weakly-evidenced association.
#:
#: This is a real trade and the two settings disagree on this data:
#: strict recommends SAMPLE6 alone (81.2, 13 of 16 gaps, carries nothing);
#: weighted recommends SAMPLE6+SAMPLE4 (100.0, 16 of 16 gaps, carries coeliac and
#: alopecia, both association-only). Neither recommends a candidate carrying
#: colorectal cancer - that is the veto, and it is not configurable.
ASSOCIATION_WEIGHTING: Final[str] = "strict"

#: Why an exposure no longer blocks, printed with every one of them.
RETIREMENT_NOTE: Final = (
    "Resemblance to a disease pattern is descriptive context. It is not a feature that "
    "can be transferred, and it no longer excludes a candidate. Candidate mechanisms are "
    "counted only where the sequence evidence was measured in that donor \u2014 see the "
    "measured-mechanism ledger."
)


def mechanism_present(mechanism: Mechanism, material: Any) -> tuple[bool, str]:
    """Is the organism the transfer evidence runs through in this material?

    Returns `(present, statement)`. An unrunnable check counts as present, so
    a missing lane can never quietly clear a candidate.
    """
    if mechanism.check == "segatella_copri_complex":
        from openbiota import copri

        extended = material.results.get("extended_catalogue")
        if not extended:
            return True, (
                "The extended SGB lane did not run for this material, so the complex could not "
                "be checked; the exposure is left counting against the candidate."
            )
        resolved = copri.resolve(extended.get("sgbs"))
        if resolved["status"] == "absent":
            return False, (
                f"No member of {mechanism.label} was detected in this material. "
                f"{resolved['clades_assessed']} clades were assessed, including all "
                f"{resolved['human_clades_assessed']} found in humans."
            )
        if resolved["status"] == "resolved":
            return True, (
                f"{mechanism.label} is present at {resolved['total_percent']:.3f}% "
                f"({len(resolved['clades_present'])} clade(s) resolved). No clade has an "
                "established disease direction, so this is presence of the organism the "
                "experiments used, not a clade-specific finding."
            )
        return True, str(resolved.get("reason") or "the complex could not be checked")
    return True, f"No check is implemented for {mechanism.label}."

#: Short names, because a table cell has no room for "Colorectal
#: cancer-associated community pattern".
SHORT_LABELS: Final[Mapping[str, str]] = {
    "crc": "Colorectal cancer",
    "adenoma": "Colorectal adenoma",
    "mdd": "Major depression",
    "csu": "Chronic hives (CSU)",
    "ra": "Rheumatoid arthritis",
    "ms": "Multiple sclerosis",
    "ankylosing_spondylitis": "Ankylosing spondylitis",
    "sle": "Lupus",
    "ad_clinical": "Alzheimer's disease",
    "ad_mci": "Mild cognitive impairment",
    "ad_preclinical_amyloid": "Preclinical amyloid",
    "longcovid": "Long COVID",
    "acvd": "Atherosclerosis",
    "ckd": "Chronic kidney disease",
    "t2d": "Type 2 diabetes",
    "t1d": "Type 1 diabetes",
    "masld": "Fatty liver (MASLD)",
    "ibd": "Inflammatory bowel disease",
    "ibs_d": "IBS with diarrhoea",
    "ibs_c": "IBS with constipation",
    "ibs_m": "IBS, mixed type",
    "celiac": "Coeliac disease",
    "mecfs": "ME/CFS",
    "alopecia_areata": "Alopecia areata",
    "androgenetic_alopecia": "Androgenetic alopecia",
    "symptomatic_dermographism": "Symptomatic dermographism",
}


def short_label(profile: str, label: str) -> str:
    if profile in SHORT_LABELS:
        return SHORT_LABELS[profile]
    return label.split(" — ")[0].split(" (")[0].strip()


@dataclass(frozen=True)
class Exposure:
    """One pattern a candidate would introduce that the recipient lacks."""

    profile: str
    label: str
    donor_percentile: float
    recipient_percentile: float
    donor_person_ids: tuple[str, ...]
    tier: str
    statement: str
    boundary: str
    sources: tuple[str, ...]
    mechanism_label: str = ""
    mechanism_status: str = ""
    #: Species this profile marks as enriched in its condition that the
    #: donor carries and the recipient does not. Empty when the percentile
    #: flagged the pattern but no such organism would actually be brought.
    introduced_organisms: tuple[str, ...] = ()

    @property
    def actionable(self) -> bool:
        """Whether this exposure disqualifies a candidate.

        Two conditions, both required. The evidence tier must be one where
        patient stool reproduced disease features in recipient animals
        against healthy-donor controls, and the donor must actually be
        bringing an organism the disease is enriched for that the
        recipient does not already carry. A high percentile on its own is
        not enough: it is a composite of many species in both directions,
        and the recipient may already have every organism behind it.
        """
        return self.tier in ACTIONABLE_TIERS and bool(self.introduced_organisms)

    def to_json(self) -> dict[str, Any]:
        return {
            "profile": self.profile,
            "label": self.label,
            "short_label": short_label(self.profile, self.label),
            "donor_percentile": self.donor_percentile,
            "recipient_percentile": self.recipient_percentile,
            "donor_person_ids": list(self.donor_person_ids),
            "transfer_evidence_tier": self.tier,
            "evidence": self.statement,
            "evidence_boundary": self.boundary,
            "sources": list(self.sources),
            "introduced_organisms": list(self.introduced_organisms),
            "mechanism": self.mechanism_label,
            "mechanism_status": self.mechanism_status,
            "counts_against_selection": self.actionable,
        }


def recipient_patterns(recipient: Any) -> dict[str, tuple[str, float]]:
    """profile → (label, percentile) for the recipient."""
    out: dict[str, tuple[str, float]] = {}
    sim = (recipient.results.get("profile_similarity") or {})
    for row in sim.get("ranked") or []:
        pct = row.get("percentile")
        if pct is None:
            continue
        out[str(row.get("profile"))] = (str(row.get("label") or row.get("profile")), float(pct))
    return out


def _detected_species(material: Any) -> frozenset[str]:
    """Named species detected in the lane the matching runs on."""
    from .suitability import _species_abundance

    return frozenset(_species_abundance(material))


def _species_percent(material: Any) -> dict[str, float]:
    """Species abundance in the lane the matching runs on.

    Delegates rather than re-reading the catalogue: this module and
    `suitability` were reading MetaPhlAn 4 and MetaPhlAn 3 respectively,
    and the two lanes disagree about which organisms a donor has.
    """
    from .suitability import _species_abundance

    return _species_abundance(material)


def exposures_for(
    *,
    recipient: Any,
    donors: Sequence[Any],
    person_of: Mapping[str, str],
    profile_set: Any = None,
) -> list[Exposure]:
    """Patterns these donors would introduce that the recipient does not have.

    "Does not have" is decided at organism level when ``profile_set`` is
    supplied: the donor must carry a species the profile marks as enriched
    in its condition, and the recipient must not carry it. Without the
    profile set the function falls back to the percentile comparison, which
    is weaker and is documented as such at the call site.
    """
    mine = recipient_patterns(recipient)
    recipient_organisms = _detected_species(recipient)
    found: dict[str, Exposure] = {}
    for donor in donors:
        donor_organisms = _detected_species(donor)
        donor_abundance = _species_percent(donor)
        label_person = person_of.get(donor.material.material_id, donor.material.person_id)
        for row in ((donor.results.get("profile_similarity") or {}).get("ranked") or []):
            profile = str(row.get("profile"))
            pct = row.get("percentile")
            if pct is None:
                continue
            donor_pct = float(pct)
            rec_label, rec_pct = mine.get(profile, (profile, 0.0))

            # Which disease-enriched organisms this donor would actually
            # bring that the recipient does not already carry.
            #
            # This replaces a percentile comparison as the test for "the
            # recipient already has this". A percentile is a composite of
            # many organisms in two directions, so two people can sit at
            # the same rank on entirely different species, and a donor
            # could bring organisms the recipient lacks while the
            # percentile gap said there was nothing to introduce. On this
            # data that is not hypothetical: the recipient carries no
            # E. coli, one donor does, E. coli is enriched in ten separate
            # transfer-evidenced profiles, and every one of them was
            # silenced because the recipient's own composite percentile was
            # already high.
            # Abundance, not presence. An organism counts only where this
            # donor carries more of it than healthy adults normally do;
            # `suitability` applies the same test, and the two must not
            # drift apart or the front page and the detail table disagree
            # about the same donor. See suitability.ABNORMAL_ABUNDANCE_PERCENTILE.
            introduced = ()
            if profile_set is not None:
                from .suitability import _reference_quantiles

                reference = _reference_quantiles("refs/taxonomic_cohort.json")
                candidates = (
                    adverse_organisms(profile, profile_set)
                    & donor_organisms - recipient_organisms
                )
                introduced = tuple(sorted(
                    species for species in candidates
                    if (reference.get(species) or (None, None, None))[1] is not None
                    and donor_abundance.get(species, 0.0) > reference[species][1]
                ))

            # Report when the donor brings disease-enriched organisms the
            # recipient lacks, or when the old percentile rule fires. The
            # organism test is the one that governs selection; the
            # percentile is kept as context because it is what the
            # resemblance score actually measured.
            percentile_flag = (
                donor_pct >= HIGH_PERCENTILE and donor_pct - rec_pct >= NEW_MARGIN
            )
            if not introduced and not percentile_flag:
                continue
            ev = EVIDENCE.get(profile)
            existing = found.get(profile)
            people = (label_person,) if existing is None else tuple(
                sorted({*existing.donor_person_ids, label_person})
            )
            # Where the evidence runs through a named organism, check this
            # donor for it. Absent organism, absent pathway.
            tier = ev.tier if ev else "not_assessed"
            mech_label = mech_status = ""
            if ev is not None and ev.mechanism is not None:
                mech_label = ev.mechanism.label
                present, mech_status = mechanism_present(ev.mechanism, donor)
                mech_status = f"{ev.mechanism.note} {mech_status}"
                if not present:
                    tier = MECHANISM_ABSENT
            found[profile] = Exposure(
                profile=profile,
                label=str(row.get("label") or rec_label),
                donor_percentile=max(donor_pct, existing.donor_percentile if existing else 0.0),
                recipient_percentile=rec_pct,
                donor_person_ids=people,
                tier=tier,
                mechanism_label=mech_label,
                mechanism_status=mech_status,
                statement=(
                    ev.statement if ev else
                    "This condition is not in the transfer-evidence catalogue, so no transfer "
                    "claim is made either way."
                ),
                boundary=(
                    ev.boundary if ev else
                    "Absence from the catalogue is a gap in our review, not evidence of safety."
                ),
                sources=ev.sources if ev else (),
                introduced_organisms=introduced,
            )
    return sorted(
        found.values(),
        key=lambda e: (not e.actionable, -e.donor_percentile, e.profile),
    )


def candidate_exposures(
    all_exposures: Sequence[Exposure], member_person_ids: Sequence[str]
) -> list[Exposure]:
    """The subset introduced by this candidate's members."""
    members = set(member_person_ids)
    return [e for e in all_exposures if members & set(e.donor_person_ids)]
