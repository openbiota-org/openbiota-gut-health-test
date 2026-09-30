"""What the pipeline knows about the person, what it assumed, and why it matters.

A stool metagenome is anonymous. Everything else about the person —
medications, recent antibiotics, pregnancy, immune status, age — arrives only
if someone supplies it, on the command line or in a sample manifest. Several
readings depend on those facts: metformin produces much of the published type
2 diabetes signature, proton-pump inhibitors push mouth bacteria into stool,
a recent antibiotic course flattens everything, and every intervention card
has safety rules keyed on pregnancy, immune status and age.

The question is what to do when a fact is *not* supplied. There are two
honest answers, and the report mode chooses between them:

* **participant** mode fails closed. An unknown is an unknown; anything that
  depends on it is withheld. This is the regulatory posture for a report a
  person reads about themselves without a clinician (spec v04 §8.7).
* **research** and **clinician** modes assume the ordinary case — no
  medication, no recent antibiotics, not pregnant, immune-competent, adult —
  score everything, and *say so*. Each assumption is recorded here with the
  concrete way the result would change if it were wrong, and that ledger is
  printed on the reading-guide page, on every profile that leaned on an
  assumption, and in ``results.json``.

The second posture is the one a researcher or clinician actually wants: a
score they can re-read against the stated assumption, rather than a blank.
What is never acceptable is a silent default — assuming "no metformin" and
not telling anyone is how a drug effect gets read as disease.

Fields the pipeline cannot reasonably default (sex, age in years, country)
are never assumed a value. Age is the one partial exception: for the safety
rules only, an adult is assumed, because a stool sample submitted for adult
screening is overwhelmingly from an adult and the alternative is to show a
clinician no intervention content at all.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

MODES_THAT_ASSUME: Final = frozenset({"research", "clinician"})

KIND_MEDICATION: Final = "medication"
KIND_EXPOSURE: Final = "exposure"
KIND_SAFETY: Final = "safety"
KIND_DEMOGRAPHIC: Final = "demographic"


@dataclass(frozen=True, slots=True)
class ContextField:
    """One fact about the person that some reading depends on."""

    key: str
    label: str
    kind: str
    #: Value assumed when not supplied, in the modes that assume. ``None``
    #: means the field is never assumed and stays unknown.
    assumed: Any
    #: Which readings lean on it, for the ledger.
    used_by: str
    #: How the results would differ if the fact were present after all.
    if_present: str
    #: Names that identify it in a free-text medication list.
    aliases: tuple[str, ...] = ()
    #: Dotted path into the safety-rule context, when a rule reads it.
    safety_path: str | None = None

    @property
    def assumable(self) -> bool:
        return self.assumed is not None


# Order is the order the ledger prints in: the ones that move scores first.
CONTEXT_FIELDS: Final[tuple[ContextField, ...]] = (
    ContextField(
        key="metformin", label="metformin", kind=KIND_MEDICATION, assumed=False,
        aliases=("metformin", "glucophage", "fortamet", "glumetza", "riomet"),
        used_by="the type 2 diabetes pattern, and the Escherichia and bile-salt-hydrolase readings it contains",
        if_present=(
            "Metformin raises Escherichia, lowers Intestinibacter and shifts bile-salt hydrolase carriage "
            "(Forslund 2015). The type 2 diabetes pattern would then be partly measuring the drug, and its "
            "resemblance percentile would read higher than the disease alone warrants; it would need to be "
            "re-read against the metformin-treated stratum."
        ),
    ),
    ContextField(
        key="ppi", label="proton-pump inhibitor", kind=KIND_MEDICATION, assumed=False,
        aliases=("omeprazole", "esomeprazole", "lansoprazole", "pantoprazole", "rabeprazole",
                 "dexlansoprazole", "nexium", "prilosec", "protonix", "prevacid", "ppi", "proton pump"),
        used_by="the oral-origin organism group, several disease patterns weighted on oral taxa, and the community diversity strips",
        if_present=(
            "Proton-pump inhibitors are the medication with the largest known effect on stool composition: "
            "they let Streptococcus, Veillonella and other mouth bacteria reach the colon and lower diversity "
            "(Imhann 2016, Jackson 2016). The oral-origin group and any pattern that leans on oral taxa would "
            "read higher than the underlying gut state."
        ),
    ),
    ContextField(
        key="recent_antibiotics", label="antibiotics in the last 90 days", kind=KIND_EXPOSURE, assumed=False,
        used_by="every disease pattern (all abstain within 90 days of a course), the diversity strips and the butyrate-producer group",
        if_present=(
            "A recent antibiotic course depresses diversity, butyrate producers and most of the species the "
            "disease patterns are built on. Every pattern would abstain, and the low readings in the organism "
            "section would be attributed to the drug rather than to the person's baseline community."
        ),
    ),
    ContextField(
        key="statin", label="statin", kind=KIND_MEDICATION, assumed=False,
        aliases=("atorvastatin", "rosuvastatin", "simvastatin", "pravastatin", "lovastatin", "pitavastatin",
                 "fluvastatin", "lipitor", "crestor", "zocor", "statin"),
        used_by="the general gut-health index and the cardiometabolic patterns",
        if_present=(
            "Statin use is associated with a lower prevalence of the Bacteroides 2 enterotype and a less "
            "disturbed community in obese adults (Vieira-Silva 2020); the general index and the "
            "cardiometabolic patterns would need to be read with the drug in view."
        ),
    ),
    ContextField(
        key="glp1_agonist", label="GLP-1 receptor agonist", kind=KIND_MEDICATION, assumed=False,
        aliases=("semaglutide", "liraglutide", "dulaglutide", "tirzepatide", "exenatide", "ozempic",
                 "wegovy", "mounjaro", "zepbound", "victoza", "trulicity", "glp-1", "glp1"),
        used_by="the cardiometabolic patterns and the bile-acid readings",
        if_present=(
            "GLP-1 agonists slow gastric emptying and change bile-acid handling; the type 2 diabetes, obesity "
            "and MASLD patterns would be read against a treated rather than an untreated community."
        ),
    ),
    ContextField(
        key="laxative_or_bowel_prep", label="laxative or bowel preparation", kind=KIND_EXPOSURE, assumed=False,
        aliases=("laxative", "polyethylene glycol", "miralax", "bisacodyl", "senna", "lactulose",
                 "bowel prep", "colonoscopy prep", "picolax", "moviprep"),
        used_by="the diversity strips, the mucin-degrader group and every pattern",
        if_present=(
            "Osmotic laxatives and colonoscopy preparation transiently collapse diversity and enrich "
            "Proteobacteria and mucin degraders; readings taken within about two weeks would reflect the "
            "preparation rather than the baseline community."
        ),
    ),
    ContextField(
        key="immunosuppressant", label="immunosuppressant", kind=KIND_MEDICATION, assumed=False,
        aliases=("prednisone", "prednisolone", "methylprednisolone", "tacrolimus", "cyclosporine",
                 "mycophenolate", "azathioprine", "methotrexate", "infliximab", "adalimumab", "vedolizumab",
                 "ustekinumab", "rituximab", "chemotherapy"),
        used_by="the safety rules on live-microbe products (probiotics, fermented foods, FMT)",
        safety_path="context.immunocompromise.status",
        if_present=(
            "Live-microbe interventions carry case reports of bacteraemia and fungaemia in immunocompromised "
            "adults; the safety rules would move every probiotic and fermented-food card to clinician review."
        ),
    ),
    ContextField(
        key="anticoagulant", label="anticoagulant or antiplatelet", kind=KIND_MEDICATION, assumed=False,
        aliases=("warfarin", "coumadin", "apixaban", "eliquis", "rivaroxaban", "xarelto", "dabigatran",
                 "pradaxa", "edoxaban", "clopidogrel", "plavix", "heparin", "enoxaparin"),
        used_by="the safety rules on vitamin K2 and on botanical supplements with antiplatelet activity",
        safety_path="context.bleeding_disorder",
        if_present=(
            "Vitamin K2 and several botanicals (garlic, ginkgo, high-dose fish oil) interact with "
            "anticoagulation; those cards would carry a caution and move to clinician review."
        ),
    ),
    ContextField(
        key="pregnant", label="pregnancy", kind=KIND_SAFETY, assumed=False,
        safety_path="context.pregnant",
        used_by="the safety rules on every supplement, botanical and drug card",
        if_present=(
            "Most botanical and high-dose supplement evidence excludes pregnancy; those cards would be "
            "restricted to clinician review and berberine in particular would be flagged as contraindicated."
        ),
    ),
    ContextField(
        key="breastfeeding", label="breastfeeding", kind=KIND_SAFETY, assumed=False,
        safety_path="context.breastfeeding",
        used_by="the safety rules on botanical and drug cards",
        if_present="Botanical and drug cards would be restricted to clinician review.",
    ),
    ContextField(
        key="immunocompromised", label="immune compromise", kind=KIND_SAFETY, assumed=False,
        safety_path="context.immunocompromise.status",
        used_by="the safety rules on live-microbe products",
        if_present=(
            "Probiotic, fermented-food and FMT cards would move to clinician review because of the "
            "documented risk of invasive infection from live organisms."
        ),
    ),
    ContextField(
        key="central_venous_catheter", label="central venous catheter", kind=KIND_SAFETY, assumed=False,
        safety_path="context.central_venous_catheter",
        used_by="the safety rules on live-microbe products",
        if_present="Live-microbe cards would be suppressed: catheter-associated probiotic bacteraemia is documented.",
    ),
    ContextField(
        key="critical_illness", label="critical illness or ICU care", kind=KIND_SAFETY, assumed=False,
        safety_path="context.critical_illness",
        used_by="the safety rules on live-microbe products",
        if_present="Live-microbe cards would be suppressed.",
    ),
    ContextField(
        key="short_bowel_or_structural_gi", label="short bowel or structural GI disease", kind=KIND_SAFETY,
        assumed=False, safety_path="context.short_bowel_or_structural_gi",
        used_by="the safety rules on fibre, prebiotic and live-microbe cards",
        if_present="Fibre, prebiotic and live-microbe cards would move to clinician review (obstruction and D-lactic acidosis risk).",
    ),
    ContextField(
        key="active_ibd_flare", label="active IBD flare", kind=KIND_SAFETY, assumed=False,
        safety_path="context.active_ibd_flare",
        used_by="the safety rules on fibre and diet cards",
        if_present="High-fibre and fermentable-substrate cards would move to clinician review.",
    ),
    ContextField(
        key="kidney_disease", label="kidney disease", kind=KIND_SAFETY, assumed=False,
        safety_path="context.kidney_disease",
        used_by="the safety rules on mineral, oxalate and protein-related cards",
        if_present="Magnesium, potassium-rich diet and high-oxalate cards would carry a caution.",
    ),
    ContextField(
        key="liver_disease", label="liver disease", kind=KIND_SAFETY, assumed=False,
        safety_path="context.liver_disease",
        used_by="the safety rules on botanical cards",
        if_present="Botanical cards with hepatotoxicity signals would move to clinician review.",
    ),
    ContextField(
        key="kidney_stone_history", label="kidney stone history", kind=KIND_SAFETY, assumed=False,
        safety_path="context.kidney_stone_history",
        used_by="the safety rules on the oxalate reading and vitamin C cards",
        if_present="Oxalate-related cards would carry a caution.",
    ),
    ContextField(
        key="gi_bleeding_history", label="GI bleeding history", kind=KIND_SAFETY, assumed=False,
        safety_path="context.gi_bleeding_history",
        used_by="the safety rules on antiplatelet botanicals",
        if_present="Garlic, ginkgo and fish-oil cards would carry a caution.",
    ),
    ContextField(
        key="eating_disorder_or_malnutrition", label="eating disorder or malnutrition", kind=KIND_SAFETY,
        assumed=False, safety_path="context.eating_disorder_or_malnutrition",
        used_by="the safety rules on every diet card",
        if_present="Restrictive-diet cards would be suppressed.",
    ),
    ContextField(
        key="surgery_planned_within_14_days", label="surgery in the next 14 days", kind=KIND_SAFETY,
        assumed=False, safety_path="context.surgery_planned_within_14_days",
        used_by="the safety rules on antiplatelet botanicals and live-microbe products",
        if_present="Those cards would carry a pre-operative caution.",
    ),
    ContextField(
        key="adult", label="age 18 or over", kind=KIND_DEMOGRAPHIC, assumed=True,
        used_by="the safety rules that withhold adult-derived intervention content from minors",
        if_present=(
            "If the person is under 18, every supplement, probiotic, diet and medication card would be "
            "withheld: none of the underlying evidence was generated in children."
        ),
    ),
    # Never assumed: there is no ordinary case to assume.
    ContextField(
        key="sex", label="sex", kind=KIND_DEMOGRAPHIC, assumed=None,
        used_by="reference-cohort matching and any sex-specific pattern",
        if_present="Sex-specific patterns could be scored and the reference cohort matched more tightly.",
    ),
    ContextField(
        key="age", label="age in years", kind=KIND_DEMOGRAPHIC, assumed=None,
        used_by="reference-cohort matching, the age-band strata, and the microbiome-age comparison",
        if_present=(
            "The reference cohort would be matched by age band and the microbiome-age page would show the "
            "difference between the estimate and the actual age."
        ),
    ),
    ContextField(
        key="country", label="country", kind=KIND_DEMOGRAPHIC, assumed=None,
        used_by="reference-cohort matching",
        if_present="The reference cohort would be matched by country, which is the strongest single predictor of stool composition.",
    ),
)

FIELD_BY_KEY: Final = {f.key: f for f in CONTEXT_FIELDS}


@dataclass(frozen=True, slots=True)
class Assumption:
    """One fact the run proceeded without, and what it assumed instead."""

    field: ContextField
    value: Any

    def sentence(self) -> str:
        assumed = {True: "yes", False: "none"}.get(self.value, str(self.value))
        if self.field.kind == KIND_MEDICATION:
            assumed = "not taken"
        elif self.field.kind == KIND_EXPOSURE:
            assumed = "none"
        elif self.field.kind == KIND_SAFETY:
            assumed = "no"
        elif self.field.key == "adult":
            assumed = "yes"
        return f"{self.field.label}: assumed {assumed}"

    def to_json(self) -> dict[str, Any]:
        return {
            "field": self.field.key, "label": self.field.label, "kind": self.field.kind,
            "assumed_value": self.value, "used_by": self.field.used_by,
            "if_present": self.field.if_present,
        }


@dataclass(slots=True)
class ContextLedger:
    """Everything supplied, everything assumed, everything left unknown."""

    mode: str
    supplied: dict[str, Any] = field(default_factory=dict)
    assumptions: list[Assumption] = field(default_factory=list)
    unknown: list[ContextField] = field(default_factory=list)
    #: Free-text medication entries that matched no known class. Kept so a
    #: reader can see that the list was read, not just the classes.
    unclassified_medications: list[str] = field(default_factory=list)

    @property
    def assumes(self) -> bool:
        return self.mode in MODES_THAT_ASSUME

    @property
    def assumed_keys(self) -> frozenset[str]:
        return frozenset(a.field.key for a in self.assumptions)

    def value(self, key: str) -> Any:
        """Supplied value, else assumed value, else ``None``."""
        if key in self.supplied:
            return self.supplied[key]
        for a in self.assumptions:
            if a.field.key == key:
                return a.value
        return None

    def assumption_for(self, key: str) -> Assumption | None:
        for a in self.assumptions:
            if a.field.key == key:
                return a
        return None

    # ---- what the engines consume ------------------------------------------ #

    def metadata(self) -> dict[str, Any]:
        """Flat mapping for the profile abstention rules and matching."""
        out: dict[str, Any] = {}
        for f in CONTEXT_FIELDS:
            v = self.value(f.key)
            if v is not None:
                out[f.key] = v
        # the historical key the antibiotic rules read
        out["antibiotics"] = bool(self.value("recent_antibiotics"))
        return out

    def safety_context(
        self, *, subject_age: float | None, subject_sex: str | None, subject_country: str | None,
    ) -> dict[str, Any]:
        """The nested context the safety-rule DSL evaluates.

        In the assuming modes every safety field the rules can read is set to
        its assumed value, so no rule falls through to its fail-closed unknown
        branch for want of a manifest; the assumption is then printed instead.
        ``adult_status`` carries the age posture explicitly: ``minor``,
        ``adult``, ``assumed_adult`` or ``unknown``.
        """
        meds = list(self.supplied.get("medication_list") or [])
        if subject_age is not None:
            adult_status = "minor" if subject_age < 18 else "adult"
        elif self.assumes:
            adult_status = "assumed_adult"
        else:
            adult_status = "unknown"
        ctx: dict[str, Any] = {
            "medications": meds,
            "conditions": [], "self_reported_conditions": [], "symptoms": [], "allergies": [],
            "external_labs": [],
            "immunocompromise": {},
        }
        for f in CONTEXT_FIELDS:
            if f.safety_path is None:
                continue
            v = self.value(f.key)
            if v is None:
                continue
            _, _, path = f.safety_path.partition(".")
            head, _, tail = path.partition(".")
            if tail:
                ctx.setdefault(head, {})
                # several fields may feed one flag (immunosuppressant → immunocompromise); any yes wins
                ctx[head][tail] = bool(ctx[head].get(tail)) or bool(v)
            else:
                ctx[head] = bool(ctx.get(head)) or bool(v)
        if ctx["immunocompromise"] and "severe" not in ctx["immunocompromise"]:
            ctx["immunocompromise"]["severe"] = False
        if self.assumes:
            # a manifest-less run in an assuming mode: everything the rules can
            # read is either supplied or assumed absent
            for key in ("recent_hospitalization", "preterm_infant", "diabetes_on_glucose_lowering",
                        "bleeding_disorder"):
                ctx.setdefault(key, False)
            # eGFR is only read by the kidney rules; "no kidney disease" assumed
            # means normal function for their purpose, so the numeric branch
            # does not fall through to unknown while the boolean one is assumed.
            if not ctx.get("kidney_disease"):
                ctx.setdefault("egfr", 90.0)
            ctx["immunocompromise"].setdefault("status", False)
            ctx["immunocompromise"].setdefault("severe", False)
            if self.value("metformin") or self.value("glp1_agonist"):
                ctx["diabetes_on_glucose_lowering"] = True
        return {
            "subject": {
                "age_years": subject_age, "sex": subject_sex, "country": subject_country,
                "mode": self.mode, "age_known": subject_age is not None,
                "is_minor": (
                    None if subject_age is None and not self.assumes
                    else (subject_age is not None and subject_age < 18)
                ),
                "adult_status": adult_status,
            },
            "context": ctx,
        }

    # ---- what the report prints -------------------------------------------- #

    def headline(self) -> str:
        """One sentence for the reading-guide page."""
        n_sup = len([k for k in self.supplied if k != "medication_list"])
        if not self.assumes:
            return (
                "This report was generated in participant mode: any fact about you that was not supplied "
                "is treated as unknown, and every reading that depends on it is withheld rather than guessed."
            )
        if not self.assumptions:
            return "Every fact the readings depend on was supplied; nothing about you was assumed."
        return (
            f"{'No' if n_sup == 0 else str(n_sup)} personal fact{'s were' if n_sup != 1 else ' was'} supplied "
            f"with this sample, so the report proceeds on {len(self.assumptions)} stated assumptions rather "
            "than withholding readings. Each one is listed below with how the result would change if it "
            "were wrong. Supplying the missing facts (--medications, --subject-age, --subject-sex, "
            "--subject-country, --antibiotics-days-ago, or a sample manifest) replaces the assumption with "
            "the fact and re-scores."
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "what_this_is": (
                "Facts about the person that readings depend on: which were supplied, which were assumed "
                "because they were not supplied, and which stay unknown because no default is honest."
            ),
            "mode": self.mode,
            "assumes_when_missing": self.assumes,
            "supplied": dict(self.supplied.items()),
            "assumed": [a.to_json() for a in self.assumptions],
            "unknown": [
                {"field": f.key, "label": f.label, "used_by": f.used_by, "if_present": f.if_present}
                for f in self.unknown
            ],
            "unclassified_medications": list(self.unclassified_medications),
        }


# --------------------------------------------------------------------------- #
# resolution
# --------------------------------------------------------------------------- #

_SPLIT: Final = re.compile(r"[;,/]+|\s+and\s+|\s*\+\s*")


def classify_medications(entries: Sequence[str]) -> tuple[dict[str, bool], list[str]]:
    """Map a free-text medication list onto the medication classes above.

    Returns ``(flags, unclassified)``: every medication-kind field that matched
    is ``True``; entries matching nothing are returned so they can be shown.
    Matching is case-insensitive substring on the alias list, which is the
    right tolerance for "Metformin 500mg bd" or "omeprazole (20 mg)".
    """
    flags: dict[str, bool] = {}
    unclassified: list[str] = []
    for raw in entries:
        for part in (p.strip() for p in _SPLIT.split(raw) if p.strip()):
            low = part.lower()
            hit = False
            for f in CONTEXT_FIELDS:
                if f.aliases and any(a in low for a in f.aliases):
                    flags[f.key] = True
                    hit = True
            if not hit:
                unclassified.append(part)
    return flags, unclassified


def resolve_context(
    *,
    mode: str,
    medications: Sequence[str] | None = None,
    antibiotics_days_ago: int | None = None,
    subject_age: float | None = None,
    subject_sex: str | None = None,
    subject_country: str | None = None,
    manifest: Mapping[str, Any] | None = None,
) -> ContextLedger:
    """Build the ledger from whatever was supplied.

    ``manifest`` is an optional flat mapping of ``field key → value`` (the
    sample-manifest namespace); explicit command-line values win over it.
    """
    ledger = ContextLedger(mode=mode)
    manifest = dict(manifest or {})

    # ---- supplied --------------------------------------------------------- #
    med_list = [m for m in (medications or []) if str(m).strip()]
    if med_list:
        flags, unclassified = classify_medications([str(m) for m in med_list])
        ledger.supplied["medication_list"] = [str(m) for m in med_list]
        ledger.supplied.update(flags)
        ledger.unclassified_medications = unclassified
        # a supplied list is evidence of absence for the classes it lacks
        for f in CONTEXT_FIELDS:
            if f.kind == KIND_MEDICATION and f.key not in flags:
                ledger.supplied[f.key] = False
    if antibiotics_days_ago is not None:
        ledger.supplied["recent_antibiotics"] = antibiotics_days_ago < 90
        ledger.supplied["antibiotics_days_ago"] = antibiotics_days_ago
    if subject_age is not None:
        ledger.supplied["age"] = subject_age
        ledger.supplied["adult"] = subject_age >= 18
    if subject_sex:
        ledger.supplied["sex"] = subject_sex
    if subject_country:
        ledger.supplied["country"] = subject_country
    for f in CONTEXT_FIELDS:
        if f.key in manifest and manifest[f.key] is not None and f.key not in ledger.supplied:
            ledger.supplied[f.key] = manifest[f.key]

    # ---- assumed or unknown ----------------------------------------------- #
    for f in CONTEXT_FIELDS:
        if f.key in ledger.supplied:
            continue
        if f.assumable and ledger.assumes:
            ledger.assumptions.append(Assumption(field=f, value=f.assumed))
        else:
            ledger.unknown.append(f)
    return ledger
