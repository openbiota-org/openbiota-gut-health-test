"""Symptom and condition navigation — spec 0.8.3 section 9.3 and 16.3.

Fifteen contexts, named in the specification, each linking what this report
already measured to the literature about a condition. Three rules shape the
whole module and are enforced by its tests.

**Nothing here is a new score.** Where a disease profile already exists it is
linked exactly as it stands. This module never recomputes one, never
combines several into a verdict, and never introduces a probability for a
condition that has no profile. A context with no profile behind it is an
evidence card, and says so.

**A marker count is a navigation count.** That a context lists nine
contributing readings does not make it nine times more likely to apply to
you. The count is how many things are worth reading, nothing else.

**Every contributing feature is listed.** No truncation, no "and 4 more".
If a context draws on twelve readings, all twelve are named, because a
hidden contributor is one the reader cannot check.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

METHOD_CONTEXTS: Final = "ext083.context_navigation/1.0"

#: The kinds of thing a context can point at.
LINK_KINDS: Final[frozenset[str]] = frozenset({
    "profile",    # an existing disease-profile similarity, linked unchanged
    "panel",      # a functional panel in this report
    "capacity",   # a section 5.8 capacity
    "step",       # a named biotransformation step
    "view",       # an extension view
    "index",      # a pattern index, e.g. the skin ones
})


class ContextError(ValueError):
    """A context that would invent a score or hide a contributor."""


@dataclass(frozen=True)
class Link:
    """One thing this report measured, and why it is relevant here."""

    kind: str
    key: str
    why: str

    def __post_init__(self) -> None:
        if self.kind not in LINK_KINDS:
            raise ContextError(f"{self.key}: unknown link kind {self.kind!r}")
        if not self.why.strip():
            raise ContextError(f"{self.key}: a link must say why it is relevant")

    def to_json(self) -> dict[str, Any]:
        return {"kind": self.kind, "key": self.key, "why": self.why}


@dataclass(frozen=True)
class Context:
    """One symptom or condition, and how to navigate to what was measured."""

    context_id: str
    label: str
    links: tuple[Link, ...]
    #: The distinction the specification requires this card to preserve.
    distinction: str
    #: What the literature establishes, as mechanism rather than as verdict.
    mechanisms: tuple[str, ...]
    #: What this report does not know, stated before the reader asks.
    missing_information: tuple[str, ...]
    #: What could actually be done next, including doing nothing.
    follow_up: tuple[str, ...]
    sources: tuple[str, ...] = ()
    #: Existing scores that are linked as they are and never recomputed.
    preserves: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.links:
            raise ContextError(f"{self.context_id}: a context with no links is a heading")
        if not self.distinction.strip():
            raise ContextError(
                f"{self.context_id}: every context must state the distinction it keeps"
            )
        if not self.missing_information:
            raise ContextError(
                f"{self.context_id}: every context must say what it does not know"
            )

    def to_json(self) -> dict[str, Any]:
        return {
            "context_id": self.context_id, "label": self.label,
            "links": [link.to_json() for link in self.links],
            "distinction": self.distinction,
            "mechanisms": list(self.mechanisms),
            "missing_information": list(self.missing_information),
            "follow_up": list(self.follow_up),
            "sources": list(self.sources),
            "preserves": list(self.preserves),
        }


NOT_A_DIAGNOSIS: Final = (
    "This is a way of navigating what was measured, not a test for the condition. "
    "Nothing in this card raises or lowers the chance that you have it."
)


CONTEXTS: Final[tuple[Context, ...]] = (
    Context(
        context_id="context.anxiety_depression",
        label="Anxiety and depression",
        links=(
            Link("profile", "mdd", "the existing depression gut-pattern similarity, linked as it stands"),
            Link("panel", "gaba", "GABA production capacity, one of the routes this literature discusses"),
            Link("panel", "gababreak", "GABA degradation, shown beside production and never subtracted from it"),
            Link("panel", "indole", "indole formation from tryptophan, the other much-discussed route"),
            Link("panel", "tryptamine", "tryptamine formation, which competes for the same tryptophan"),
            Link("step", "neuro.gaba_synthesis", "the specific decarboxylation step, with its own limits"),
        ),
        distinction=(
            "Human microbiome associations with depression are observational and mostly "
            "cross-sectional; the behavioural findings are from animals. Neither "
            "establishes cause, and no stool measurement diagnoses anxiety or depression."
        ),
        mechanisms=(
            "Several gut bacteria decarboxylate glutamate to GABA, and GABA-modulating "
            "strains were isolated from human stool. The blood-brain barrier is not "
            "permeable to GABA, so any effect would have to run through the vagus nerve "
            "or the enteric nervous system rather than through GABA reaching the brain.",
            "Tryptophan is the shared precursor of serotonin, indoles and kynurenines. "
            "Bacteria competing for it is a plausible mechanism and is not a measurement "
            "of anybody's serotonin.",
        ),
        missing_information=(
            "Whether any of these genes are being expressed, which sequencing DNA cannot show.",
            "Your symptoms, which are not in this data at all and matter more than any of it.",
            "Whether you take medication affecting these systems, which changes the picture entirely.",
        ),
        follow_up=(
            "If you are looking for help with mood, that is a conversation with a clinician, "
            "not a microbiome result.",
            "A symptom diary alongside a repeat sample would at least make a within-person "
            "comparison possible, which a single timepoint cannot support.",
        ),
        sources=("F15", "C06"),
        preserves=("mdd",),
    ),
    Context(
        context_id="context.autism_research",
        label="Autism research",
        links=(
            Link("panel", "thiamine", "B1 synthesis, one of the pathways the multikingdom work highlighted"),
            Link("panel", "k2", "quinone synthesis, highlighted in the same work"),
            Link("view", "mycobiome", "the fungal and multikingdom component those studies used"),
            Link("view", "ecology", "community composition and type, which those cohorts stratified on"),
        ),
        distinction=(
            "The published classifiers were built in paediatric cohorts, and diet differs "
            "systematically between the groups they compared. This report has no autism "
            "classifier and produces no score; these are the pathways the research names."
        ),
        mechanisms=(
            "A 2024 multikingdom study reported a panel spanning bacteria, fungi, archaea "
            "and viruses, and highlighted B1 and quinone pathways among its features.",
            "Whether these differences are a cause, a consequence of diet and behaviour, "
            "or an artefact of cohort design is unresolved in the source literature itself.",
        ),
        missing_information=(
            "Diet, which differs between the compared groups and is not recorded here.",
            "Age: these cohorts were children, and the reference set behind this report is adult.",
            "Everything a diagnosis actually rests on, none of which is microbial.",
        ),
        follow_up=(
            "There is no action to take from this card. It exists so the pathways named in "
            "the research can be found in your own results.",
        ),
        sources=("C07", "C08"),
    ),
    Context(
        context_id="context.serotonin_melatonin_sleep",
        label="Serotonin, melatonin and sleep",
        links=(
            Link("panel", "tryptamine", "tryptamine, a distinct molecule from serotonin, made from the same precursor"),
            Link("panel", "indole", "the bacterial route that competes for that precursor"),
            Link("panel", "dopamine", "the other monoamine chemistry measured here, for contrast"),
            Link("panel", "histamine", "another amine from the same decarboxylase chemistry"),
        ),
        distinction=(
            "Serotonin, melatonin and tryptamine are three different molecules made by "
            "three different reactions. Gut serotonin and brain serotonin are separate "
            "pools that do not exchange. No number here is a concentration of any of them."
        ),
        mechanisms=(
            "Most of the body's serotonin is made in the gut by host enterochromaffin "
            "cells, and specific bacteria can increase that host production in gnotobiotic "
            "mice. That is host synthesis being modulated, not bacterial synthesis.",
            "A traceable candidate decarboxylase produces tryptamine from tryptophan. "
            "Tryptamine is not serotonin, and an enzyme that makes one is not a detector "
            "for the other.",
        ),
        missing_information=(
            "Your sleep, which is not measurable from stool.",
            "Whether the candidate enzymes here act on the substrates the papers used, "
            "which is unresolved for some of them in the source.",
        ),
        follow_up=(
            "A sleep or symptom diary is the only thing that would make this comparable "
            "over time, and it is worth more than the gene counts.",
        ),
        sources=("F17", "F18", "C09", "C10", "C11"),
    ),
    Context(
        context_id="context.atherosclerosis",
        label="Atherosclerosis and cardiovascular context",
        links=(
            Link("profile", "acvd", "the existing cardiovascular gut-pattern similarity, linked unchanged"),
            Link("panel", "cutc", "choline TMA-lyase, the most studied microbial link"),
            Link("panel", "bcaa", "amino-acid fermentation, the other route to circulating microbial metabolites"),
            Link("panel", "pcresol", "p-cresol, another microbial metabolite the liver conjugates before it circulates"),
        ),
        distinction=(
            "Microbes make TMA. Your liver makes TMAO from it. These are different steps "
            "in different organs, and the existing cardiovascular score is linked exactly "
            "as it is rather than recomputed from these genes."
        ),
        mechanisms=(
            "Gut bacteria cleave choline and carnitine to trimethylamine, which the liver "
            "oxidises to TMAO. Higher circulating TMAO is associated with cardiovascular "
            "events in cohort studies, and the causal question remains open.",
            "Host flavin-containing monooxygenase 3 activity, which varies genetically and "
            "by sex, sits between the bacterial step and the circulating molecule.",
        ),
        missing_information=(
            "Your blood TMAO, which this cannot estimate and which a laboratory can measure.",
            "Your lipids, blood pressure and family history, which dominate cardiovascular risk.",
        ),
        follow_up=(
            "Cardiovascular risk is assessed with blood pressure, lipids and history. If that "
            "is the question, those are the measurements.",
            "A measured TMAO, if you have one, can be entered in the lab results section and "
            "will be shown beside this rather than replaced by it.",
        ),
        preserves=("acvd",),
    ),
    Context(
        context_id="context.hypertension",
        label="Blood pressure context",
        links=(
            Link("profile", "hypertension", "the existing gut-pattern similarity, reused rather than duplicated"),
            Link("panel", "butyrate", "short-chain fatty acid production, the mechanism most often proposed"),
            Link("panel", "propionate", "propionate specifically, which has its own receptor literature"),
        ),
        distinction=(
            "The existing profile is the result. This card adds navigation and the "
            "functional context, and deliberately does not produce a second number for "
            "the same thing."
        ),
        mechanisms=(
            "Short-chain fatty acids act on host receptors including GPR41 and Olfr78 that "
            "have opposing effects on vascular tone in rodent experiments, which is why no "
            "single direction can be assigned to a butyrate reading.",
        ),
        missing_information=(
            "Your blood pressure, which takes thirty seconds to measure properly and is "
            "not inferable from stool.",
            "Salt intake, weight and medication, which are the things that actually move it.",
        ),
        follow_up=(
            "If blood pressure is the question, measure blood pressure. Repeatedly, at home, "
            "is better than once in a clinic.",
        ),
        preserves=("hypertension",),
    ),
    Context(
        context_id="context.auto_brewery_ethanol",
        label="Endogenous ethanol and auto-brewery evaluation",
        links=(
            Link("capacity", "capacity.ethanol", "the ethanol-forming route, newly measured"),
            Link("panel", "ethanol", "the underlying gene readings, including the ethanolamine route"),
            Link("profile", "masld", "the existing fatty-liver gut-pattern similarity, linked unchanged"),
            Link("view", "mycobiome", "yeasts, which ferment ethanol and are measured separately"),
            Link("view", "strain_resolution", "strain-level evidence, since this is a strain-level question"),
        ),
        distinction=(
            "Carrying adhE is ordinary; the liver findings concern particular "
            "high-alcohol-producing strains. Genetic potential is not measured ethanol, "
            "and auto-brewery syndrome is diagnosed with a carbohydrate challenge and "
            "serial blood alcohol measurements, not with a stool gene count."
        ),
        mechanisms=(
            "Specific high-alcohol-producing Klebsiella pneumoniae strains isolated from "
            "people with fatty liver disease produced liver changes when colonised into "
            "mice. Strains of the same species differ enormously in how much they produce.",
            "Acetaldehyde, the intermediate, is the genotoxic part of this chemistry, and "
            "is handled by an enzyme family too broad to measure specifically here.",
        ),
        missing_information=(
            "Whether the organisms carrying these genes here are high producers, which "
            "requires culture and measurement rather than sequence.",
            "Any actual blood alcohol measurement, which is what the diagnosis rests on.",
        ),
        follow_up=(
            "If episodes of unexplained intoxication are the reason for asking, that is a "
            "medical evaluation with a supervised carbohydrate challenge, and this card is "
            "not a substitute for it.",
        ),
        sources=("C03", "C04", "C05"),
        preserves=("masld",),
    ),
    Context(
        context_id="context.constipation",
        label="Chronic constipation",
        links=(
            Link("profile", "ibs_c", "the existing constipation-predominant IBS similarity, linked unchanged"),
            Link("panel", "methane", "methanogen capacity, associated with slower transit"),
            Link("view", "substrates", "which fibres this community can actually open, newly measured"),
            Link("panel", "butyrate", "the fermentation end-product most linked to colonic function"),
        ),
        distinction=(
            "Methanogen genes in stool are not a breath test. The association between "
            "methane and slow transit comes from breath measurements, and a gene count "
            "cannot substitute for one or diagnose intestinal methanogen overgrowth."
        ),
        mechanisms=(
            "Methane slows intestinal transit in animal models, and breath methane is "
            "associated with constipation in humans. Methanogens are archaea and are "
            "measured here by their own marker genes.",
            "Fermentable substrate reaching the colon affects stool water and bulk, which "
            "is why the substrate readings are linked here.",
        ),
        missing_information=(
            "Your stool form and frequency, which can be recorded and are more informative "
            "than any of this.",
            "Transit time, which is directly measurable and is not inferable from stool DNA.",
        ),
        follow_up=(
            "Recording stool form on the Bristol scale for two weeks costs nothing and "
            "gives a better picture than a single sequencing run.",
            "A breath test is the measurement that corresponds to the methane literature.",
        ),
        preserves=("ibs_c",),
    ),
    Context(
        context_id="context.diarrhea",
        label="Chronic diarrhoea",
        links=(
            Link("profile", "ibs_d", "the existing diarrhoea-predominant IBS similarity, linked unchanged"),
            Link("panel", "bsh", "bile salt hydrolase, the first step of bile-acid transformation"),
            Link("panel", "bai", "the 7-alpha-dehydroxylation route to secondary bile acids"),
            Link("view", "pathogens", "the pathogen screen, which is a different question with a different answer"),
        ),
        distinction=(
            "Bile-acid diarrhoea is one cause among many, and it is diagnosed with a SeHCAT "
            "scan or a measured faecal bile acid, not with these genes. The pathogen screen "
            "is linked here because it answers a separate question that matters more when "
            "diarrhoea is acute."
        ),
        mechanisms=(
            "Bacteria deconjugate and then dehydroxylate bile acids. Reduced conversion "
            "leaves more primary bile acid reaching the colon, which draws in water. The "
            "bacterial step is real; the amount arriving is set by the liver and the ileum.",
        ),
        missing_information=(
            "Whether your ileum is reabsorbing bile acids normally, which is the actual "
            "question in bile-acid diarrhoea and needs a scan or a stool measurement.",
            "Recent antibiotics and travel, which change the interpretation completely.",
        ),
        follow_up=(
            "Persistent diarrhoea, and particularly diarrhoea with weight loss, blood or "
            "fever, needs a clinician rather than a microbiome report.",
        ),
        sources=("C12",),
        preserves=("ibs_d",),
    ),
    Context(
        context_id="context.abdominal_pain",
        label="Recurring abdominal pain",
        links=(
            Link("profile", "ibs", "the existing IBS similarity, linked as it stands"),
            Link("panel", "histamine", "bacterial histamine production, a strain-specific mechanism"),
            Link("view", "fermentation", "gas-forming fermentation routes, newly measured"),
            Link("panel", "h2s", "hydrogen sulfide, which has its own visceral sensitivity literature"),
        ),
        distinction=(
            "The histamine mechanism is strain-dependent. A species that contains "
            "histamine-producing strains is not evidence that yours produces histamine, "
            "and species-level abundance cannot establish it."
        ),
        mechanisms=(
            "A specific Klebsiella aerogenes strain producing histamine was linked to "
            "visceral hypersensitivity in a mouse model and to a subgroup of patients. The "
            "work is explicitly strain-level.",
            "Gas production from fermentation distends the colon, and distension is painful "
            "at lower volumes in people with visceral hypersensitivity.",
        ),
        missing_information=(
            "Where the pain is, when it happens and what changes it, none of which is here "
            "and all of which is more diagnostic than this report.",
            "Whether the strains present are the producing ones, which needs strain-level "
            "resolution this depth cannot always reach.",
        ),
        follow_up=(
            "Pain with weight loss, bleeding, fever or waking you at night is a reason to "
            "see a clinician promptly, whatever a microbiome report says.",
        ),
        sources=("C13",),
        preserves=("ibs",),
    ),
    Context(
        context_id="context.gluten_fructan",
        label="Wheat, gluten and fructan sensitivity",
        links=(
            Link("profile", "celiac", "the existing coeliac gut-pattern similarity, linked unchanged"),
            Link("view", "substrates", "fructan utilisation specifically, newly measured"),
            Link("panel", "carbohydrates", "the wider carbohydrate-active enzyme readings"),
        ),
        distinction=(
            "Coeliac disease, wheat allergy and non-coeliac wheat sensitivity are three "
            "different conditions with three different mechanisms. Coeliac disease is "
            "diagnosed serologically and by biopsy, while still eating gluten. Nothing "
            "here tests for any of them."
        ),
        mechanisms=(
            "Trials that separated gluten from fructans found that fructans, not gluten, "
            "reproduced symptoms in many people who believed they were gluten-sensitive. "
            "Wheat carries both, which is why the substrate reading is linked here.",
        ),
        missing_information=(
            "Whether you are currently eating gluten, which determines whether coeliac "
            "testing would even be valid.",
            "Your symptoms and their timing relative to meals.",
        ),
        follow_up=(
            "If coeliac disease has never been excluded, that test is done before removing "
            "gluten, not after: stopping first makes the test unreliable.",
            "A supervised fructan reduction is the way to separate the two, and it is "
            "temporary rather than permanent.",
        ),
        sources=("C14",),
        preserves=("celiac",),
    ),
    Context(
        context_id="context.ibs",
        label="Irritable bowel syndrome",
        links=(
            Link("profile", "ibs", "the existing overall IBS similarity"),
            Link("profile", "ibs_c", "the constipation-predominant variant"),
            Link("profile", "ibs_d", "the diarrhoea-predominant variant"),
            Link("profile", "ibs_m", "the mixed variant"),
            Link("view", "substrates", "which fermentable substrates this community can open"),
            Link("view", "fermentation", "the gas-forming routes downstream of them"),
        ),
        distinction=(
            "All four IBS scores already exist in this report and are shown as they are. "
            "This card adds navigation and substrate context only, and produces no new "
            "score. IBS is a clinical diagnosis made on symptom criteria."
        ),
        mechanisms=(
            "The microbiome literature in IBS is large, mostly 16S, and inconsistent "
            "between cohorts, which is what a symptom-defined diagnosis tends to produce.",
            "Fermentable substrate reaching the colon is the mechanism behind the low-FODMAP "
            "response, and substrate capacity is the part of that this report can measure.",
        ),
        missing_information=(
            "Your symptom pattern against the Rome criteria, which is what the diagnosis is.",
            "What you eat, which determines how much fermentable substrate arrives.",
        ),
        follow_up=(
            "A low-FODMAP trial is a short diagnostic exercise with a structured "
            "reintroduction, not a long-term diet, and is best done with a dietitian.",
        ),
        preserves=("ibs", "ibs_c", "ibs_d", "ibs_m"),
    ),
    Context(
        context_id="context.thyroid",
        label="Thyroid context",
        links=(
            Link("panel", "bglucuronidase", "beta-glucuronidase, often invoked in this context"),
            Link("view", "ecology", "overall community composition, which the observational work compared"),
            Link("view", "input_register", "where a measured TSH or free T4 can be entered and shown"),
        ),
        distinction=(
            "Hypothyroidism, Hashimoto's thyroiditis and Graves' disease are separate "
            "conditions. Host blood markers and host gene expression cannot be inferred "
            "from stool DNA, and this report has no thyroid score."
        ),
        mechanisms=(
            "Observational studies report compositional differences in thyroid conditions. "
            "They are small, cross-sectional and confounded by treatment, and no mechanism "
            "has been established in humans.",
            "Selenium and iodine handling are frequently discussed, but the bacterial "
            "contribution to either has not been quantified in people.",
        ),
        missing_information=(
            "Your TSH, free T4 and antibody status, which are blood tests and are the "
            "actual measurements for this.",
            "Whether you take thyroid hormone, which changes everything about interpretation.",
        ),
        follow_up=(
            "If thyroid function is the question, a blood test answers it directly. Results "
            "can be entered in the lab section and will be shown alongside, not merged in.",
        ),
        sources=("C15", "C16"),
    ),
    Context(
        context_id="context.food_allergy",
        label="Food allergy and immune tolerance",
        links=(
            Link("profile", "csu", "the existing chronic urticaria pattern, linked unchanged"),
            Link("index", "ATD_ADULT_WANG2023_TRANSLATED_V1", "the existing atopic dermatitis pattern index"),
            Link("panel", "butyrate", "butyrate, the metabolite in the regulatory T cell literature"),
            Link("view", "substrates", "the fibre substrates upstream of that butyrate"),
        ),
        distinction=(
            "Mouse transfer and consortium experiments show that particular bacteria affect "
            "tolerance in mice. They do not establish that any food is safe or unsafe for "
            "you, and nothing here identifies an allergen."
        ),
        mechanisms=(
            "Consortia from healthy infants protected mice against cow's milk allergy, and "
            "consortia from allergic infants did not. The effect travelled with specific "
            "organisms rather than with diversity.",
            "Butyrate induces regulatory T cells in the colon in mouse work, which is the "
            "usual proposed bridge between fibre and tolerance.",
        ),
        missing_information=(
            "Which foods, if any, cause you problems, which is established by history and "
            "by supervised challenge rather than by sequencing.",
            "Your IgE status, which is a blood or skin test.",
        ),
        follow_up=(
            "Never use a microbiome result to decide whether a food is safe. Suspected "
            "allergy is assessed by an allergist, and a genuine anaphylaxis risk is not "
            "something to experiment with.",
        ),
        sources=("C17", "C18", "C19"),
        preserves=("csu",),
    ),
    Context(
        context_id="context.glp1",
        label="GLP-1-related microbial mechanisms",
        links=(
            Link("view", "biotransformation", "the GLP-1 channels, kept as separate mechanisms"),
            Link("panel", "butyrate", "the short-chain fatty acid channel"),
            Link("panel", "propionate", "the colonic propionate experiments specifically"),
            Link("panel", "bsh", "the bile-acid channel, which acts through a different receptor"),
            Link("panel", "indole", "the indole channel, whose direction depends on exposure duration"),
        ),
        distinction=(
            "Microbial genetic potential is not a hormone concentration and is not "
            "comparable to a GLP-1 medication. Host receptors are host genes and are not "
            "targets in stool sequencing. Indole does not have one universal direction: "
            "prolonged exposure reduced secretion in the source experiments."
        ),
        mechanisms=(
            "Short-chain fatty acids trigger GLP-1 secretion through FFAR2 in experimental "
            "systems, and a targeted colonic propionate delivery study in humans increased "
            "secretion.",
            "Bile acids act through TGR5 on the same cells, which is a separate mechanism "
            "that happens to converge on the same hormone.",
            "A specific Akkermansia muciniphila protein has experimental support, and "
            "detecting it requires the tested sequence rather than the species.",
        ),
        missing_information=(
            "Your actual GLP-1, which requires a blood assay with the active-versus-total "
            "distinction, fasting or postprandial timing, and correct sample handling.",
            "Whether any of these genes are expressed or their products secreted.",
        ),
        follow_up=(
            "A measured GLP-1 can be entered in the lab section with its assay and timing, "
            "and will be shown next to this rather than combined with it.",
        ),
        sources=("F43", "C01", "C02"),
    ),
    Context(
        context_id="context.eczema",
        label="Eczema",
        links=(
            Link("index", "ATD_ADULT_WANG2023_TRANSLATED_V1", "the existing atopic dermatitis pattern index, preserved"),
            Link("panel", "butyrate", "the short-chain fatty acid arm of the tolerance literature"),
            Link("view", "substrates", "the fibre substrates upstream of it"),
        ),
        distinction=(
            "The existing skin pattern index is preserved exactly and is not recomputed "
            "here. Gut composition is one factor among genetics, skin barrier function and "
            "the skin's own microbiome, which a stool sample does not sample."
        ),
        mechanisms=(
            "Filaggrin variants and skin barrier function dominate atopic dermatitis risk. "
            "Gut microbial associations are real in cohort studies and much weaker than "
            "those factors.",
        ),
        missing_information=(
            "Your skin, which is not what was sampled. The skin microbiome is a different "
            "community and would need a different swab.",
            "Whether the flares follow anything identifiable, which a diary answers.",
        ),
        follow_up=(
            "Emollients and the treatments a dermatologist recommends have far stronger "
            "evidence than anything derived from a stool sample.",
        ),
    ),
)

BY_ID: Final[Mapping[str, Context]] = {c.context_id: c for c in CONTEXTS}


# --------------------------------------------------------------------------- #
# resolution against one sample
# --------------------------------------------------------------------------- #


def _panel_names(results: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    return {
        str(p.get("name")): p
        for p in (results.get("panels") or [])
        if isinstance(p, dict) and p.get("name")
    }


def _profiles(results: Mapping[str, Any]) -> Mapping[str, Any]:
    payload = (results.get("profile_similarity") or {}).get("profiles")
    return payload if isinstance(payload, Mapping) else {}


def _indices(results: Mapping[str, Any]) -> Mapping[str, Any]:
    similarity = results.get("profile_similarity") or {}
    merged: dict[str, Any] = {}
    for key in ("skin_pattern_indices", "urticaria_pattern_indices"):
        block = similarity.get(key)
        if isinstance(block, Mapping):
            merged.update(block)
    return merged


def resolve(
    results: Mapping[str, Any], views: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Every context, with the readings behind it that this sample has.

    A link whose target is missing is reported as missing rather than
    dropped: that a context draws on something this run did not produce is
    exactly the sort of thing that should be visible.
    """
    views = views or {}
    panels = _panel_names(results)
    profiles = _profiles(results)
    indices = _indices(results)
    capacities = {
        str(c.get("capacity_id")): c
        for c in (views.get("additional_capacities") or {}).get("capacities") or []
    }
    steps = {
        str(step.get("step_id")): step
        for feature in (views.get("biotransformation") or {}).get("by_feature") or []
        for step in feature.get("steps") or []
    }

    rows: list[dict[str, Any]] = []
    for context in CONTEXTS:
        contributing: list[dict[str, Any]] = []
        missing: list[dict[str, Any]] = []
        for link in context.links:
            found: Any = None
            detail: dict[str, Any] = {}
            if link.kind == "profile":
                found = profiles.get(link.key)
            elif link.kind == "panel":
                panel = panels.get(link.key)
                found = panel
                if panel:
                    detail = {
                        "copies_per_100_genomes": panel.get("copies_per_100_genomes"),
                        "fragments": panel.get("accepted_fragments"),
                    }
            elif link.kind == "capacity":
                found = capacities.get(link.key)
                if found:
                    detail = {"state": found.get("state")}
            elif link.kind == "step":
                found = steps.get(link.key)
                if found:
                    detail = {"state": found.get("state")}
            elif link.kind == "view":
                found = views.get(link.key) or results.get(link.key)
            elif link.kind == "index":
                found = indices.get(link.key)
            row = {**link.to_json(), **detail}
            (contributing if found else missing).append(row)
        rows.append({
            **context.to_json(),
            "contributing": contributing,
            "not_in_this_run": missing,
            # Spec 16.3: all contributors listed, never a truncated count.
            "n_contributing": len(contributing),
            "navigation_count_note": (
                f"{len(contributing)} readings in this report are relevant to this topic. "
                "That is a count of things to read, not a measure of how likely the "
                "condition is."
            ),
        })

    return {
        "method": METHOD_CONTEXTS,
        "feature_id": "A13",
        "n_contexts": len(rows),
        "contexts": rows,
        "standing_note": NOT_A_DIAGNOSIS,
        "what_this_is": (
            "Fifteen topics people ask about, each pointing at what this report actually "
            "measured and at what it did not. Where a pattern score for a condition "
            "already exists it is linked exactly as calculated; no score here is new, and "
            "no topic without a score has been given one."
        ),
    }


def contexts_for(results: Mapping[str, Any], key: str, kind: str) -> list[str]:
    """Which contexts point at a given reading — the reverse link.

    Lets a panel or profile carry a list of the topics it feeds, so
    navigation runs in both directions rather than only outward from the
    context cards.
    """
    del results
    return [
        context.label
        for context in CONTEXTS
        if any(link.kind == kind and link.key == key for link in context.links)
    ]
