"""The 50 intervention seed records (BF-I01 .. BF-I50).

These are evidence cards, not recommendations. The module's job is to make
matched evidence findable and to make its limits impossible to miss - not
to assemble a supplement stack. Several design choices follow from that.

**Opposing evidence is structural, not editorial.** ``opposes`` links a
card to the records that contradict it, and :func:`retrieve` pulls those
companions in automatically. Berberine cannot be displayed for E. coli
without BF-I06's adaptation finding; L. plantarum CIRM653 cannot be
displayed without BF-I45's mouse-carriage result. The acceptance suite
asserts these pairs travel together.

**Endpoints are typed and ranked, because they are routinely conflated.**
A crystal-violet biomass decline is ``matrix_biomass``; it is not
``viable_burden`` and not ``eradication``. A study titled "disrupts
biofilm" that ran a formation assay is recorded as ``formation``. The
ranking key in section 10 sorts on the endpoint that was actually measured,
so a formation-only result cannot outrank a viable-burden result by having
a better title.

**Laboratory exposures stay laboratory exposures.** ``lab_exposure`` holds
the tested concentration verbatim. There is deliberately no field for a
suggested dose, and no code path converts one into the other.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from . import sources as src

#: What the study actually measured, ordered by how directly it speaks to
#: "did this help". Index order is used by the ranking key.
ENDPOINTS: Final = (
    "viable_burden",
    "harmful_behaviour",
    "protective_function",
    "eradication",
    "matrix_biomass",
    "formation",
    "adhesion",
    "organism_load",
    "clinical_outcome_no_biofilm_endpoint",
    "expression_only",
    "hypothesis",
)

#: Where the work was done, ordered by applicability to a human gut.
MODELS: Final = (
    "human_gut_controlled",
    "human_derived_gut_community",
    "animal_gut",
    "gut_relevant_in_vitro",
    "extraintestinal",
    "in_silico",
)

#: How well the evidence matches the finding it is retrieved for.
TARGET_MATCHES: Final = (
    "exact_tested_strain",
    "same_species",
    "supported_mechanism_other_carrier",
    "community_model_match",
    "context_only",
)

#: Whether the tested molecule plausibly reaches the colon.
DELIVERY_STATES: Final = (
    "target_exposure_measured",
    "indirect_delivery_demonstrated",
    "plausible_unverified",
    "route_site_mismatch",
    "unknown",
)

#: The four action sections from section 10.
SECTIONS: Final = (
    "support_protective_ecology",
    "reduce_matched_concerning_mechanism",
    "clarify_or_treat_established_condition",
    "investigational_precision",
)

#: Direction of the observed effect. ``unfavourable`` cards stay in the
#: matched list rather than being hidden.
DIRECTIONS: Final = ("favourable", "unfavourable", "mixed", "neutral")


@dataclass(frozen=True, slots=True)
class Intervention:
    """One intervention evidence record."""

    id: str
    compound: str
    source_ids: tuple[str, ...]
    doi: str
    study_type: str
    organisms: tuple[str, ...]
    strain: str
    model: str
    endpoint: str
    direction: str
    #: What the study found, in its own terms.
    result: str
    #: Tested concentration, verbatim. Never a dose.
    lab_exposure: str
    delivery: str
    delivery_note: str
    section: str
    #: IDs of records that contradict or qualify this one.
    opposes: tuple[str, ...] = ()
    #: Mechanisms or findings this card is retrievable for.
    addresses: tuple[str, ...] = ()
    #: Effects on protective commensals, where measured.
    collateral: str = ""
    funding_conflict: str = ""
    safety: str = ""
    #: Biofilm maturity actually tested.
    maturity: str = "not_stated"
    caution: str = ""

    def __post_init__(self) -> None:
        if self.endpoint not in ENDPOINTS:
            raise ValueError(f"{self.id}: bad endpoint {self.endpoint!r}")
        if self.model not in MODELS:
            raise ValueError(f"{self.id}: bad model {self.model!r}")
        if self.delivery not in DELIVERY_STATES:
            raise ValueError(f"{self.id}: bad delivery {self.delivery!r}")
        if self.section not in SECTIONS:
            raise ValueError(f"{self.id}: bad section {self.section!r}")
        if self.direction not in DIRECTIONS:
            raise ValueError(f"{self.id}: bad direction {self.direction!r}")
        for sid in self.source_ids:
            src.get(sid)

    @property
    def url(self) -> str:
        return f"https://doi.org/{self.doi}"

    @property
    def endpoint_rank(self) -> int:
        return ENDPOINTS.index(self.endpoint)

    @property
    def model_rank(self) -> int:
        return MODELS.index(self.model)

    @property
    def delivery_rank(self) -> int:
        return DELIVERY_STATES.index(self.delivery)

    @property
    def is_adverse(self) -> bool:
        return self.direction == "unfavourable"

    @property
    def review_label(self) -> str:
        """What the card says about itself before anything else."""
        if self.is_adverse:
            return "opposing evidence - review before considering"
        if self.direction == "mixed":
            return "mixed evidence - both directions observed"
        return "matched research evidence"

    @property
    def endpoint_plain(self) -> str:
        """The measured endpoint in words a reader can check against a claim."""
        return {
            "viable_burden": "viable organism burden",
            "harmful_behaviour": "harmful community behaviour",
            "protective_function": "protective function",
            "eradication": "eradication of established biofilm",
            "matrix_biomass": "matrix biomass (stain-based)",
            "formation": "biofilm formation (prevention assay)",
            "adhesion": "adhesion",
            "organism_load": "organism load only",
            "clinical_outcome_no_biofilm_endpoint": (
                "clinical outcome, no biofilm measured"
            ),
            "expression_only": "gene expression only",
            "hypothesis": "hypothesis",
        }[self.endpoint]


def _i(**kw: object) -> Intervention:
    return Intervention(**kw)  # type: ignore[arg-type]


_ALL: Final = (
    Intervention(
        id="BF-I01",
        compound="Allicin",
        source_ids=(),
        doi="10.3390/ijms17070979",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="UPEC CFT073 and J96",
        model="extraintestinal",
        endpoint="matrix_biomass",
        direction="favourable",
        result=(
            "At 50 ug/mL, pre-exposure reduced formation by about 33% and 17% in "
            "the two strains; 1 h post-treatment dispersed about 40% and 30%. "
            "Crystal violet biomass, microscopy, adhesion/motility and RT-qPCR. "
            "Covers both prevention and established-biofilm dispersal."
        ),
        lab_exposure="12-50 ug/mL, growth-sparing concentrations",
        delivery="unknown",
        delivery_note=(
            "Urinary isolates on plates. An oral dose is not equivalent to these "
            "well concentrations."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("ecoli_adhesion", "curli", "BF-M01", "BF-M18"),
        maturity="both_formation_and_established",
        caution=(
            "fimH presence does not predict allicin susceptibility. Dispersal is "
            "not sterilisation."
        ),
    ),
    Intervention(
        id="BF-I02",
        compound="Allicin oral delivery (garlic supplement forms)",
        source_ids=(),
        doi="10.3390/nu10070812",
        study_type="human_crossover",
        organisms=(),
        strain="",
        model="human_gut_controlled",
        endpoint="hypothesis",
        direction="neutral",
        result=(
            "13 subjects; garlic products compared by breath allyl methyl sulfide "
            "metabolite. Enteric-tablet bioequivalence 36-104%, reduced 22-57% "
            "with a high-protein meal, with marked product variability."
        ),
        lab_exposure="commercial product servings, not standardised",
        delivery="indirect_delivery_demonstrated",
        delivery_note=(
            "Supports oral generation of allicin-derived compounds. A systemic "
            "metabolite is not measured intact allicin at a colonic biofilm. "
            "Enteric formulation is plausible and is not ruled out; disintegration "
            "site, alliinase activity and meal effects must be recorded."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("ecoli_adhesion", "delivery"),
        caution="Metabolite detection is not intact colonic allicin exposure.",
    ),
    Intervention(
        id="BF-I03",
        compound="Allicin nanoemulsion with epsilon-polylysine",
        source_ids=(),
        doi="10.1016/j.foodchem.2025.142949",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="favourable",
        result=(
            "Planktonic and mature biofilm disruption in a food-preservation and "
            "raw-beef application."
        ),
        lab_exposure="formulation-specific, food matrix",
        delivery="route_site_mismatch",
        delivery_note="Food-preservation context, not an oral gut treatment.",
        section="investigational_precision",
        addresses=("ecoli_adhesion",),
        maturity="mature",
        caution="Formulation research only; does not validate oral use.",
    ),
    Intervention(
        id="BF-I04",
        compound="Berberine, chitosan, eugenol, linoleic acid, curcumin",
        source_ids=(),
        doi="10.0000/pubmed.24377137",
        study_type="in_vitro",
        organisms=("Klebsiella pneumoniae",),
        strain="7 strong-biofilm clinical isolates from 35 screened",
        model="extraintestinal",
        endpoint="formation",
        direction="favourable",
        result=(
            "Berberine MBIC 0.0635 mg/mL; chitosan and eugenol the same; linoleic "
            "acid 0.0312 mg/mL; curcumin 0.25 mg/mL. Formation endpoint, not "
            "mature eradication."
        ),
        lab_exposure="MBIC 0.0312-0.25 mg/mL depending on compound",
        delivery="unknown",
        delivery_note="Isolate context and culture concentration travel with the card.",
        section="reduce_matched_concerning_mechanism",
        addresses=("klebsiella_adhesion", "BF-M05"),
        maturity="formation_only",
    ),
    Intervention(
        id="BF-I05",
        compound="Berberine",
        source_ids=(),
        doi="10.3389/fmicb.2019.02584",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="antimicrobial-resistant isolates",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result=(
            "Sub-MIC berberine reduced formation and downregulated quorum-related "
            "genes including luxS, pfs and sdiA."
        ),
        lab_exposure="sub-MIC",
        delivery="plausible_unverified",
        delivery_note="Oral berberine is widely used; colonic active concentration unknown.",
        section="reduce_matched_concerning_mechanism",
        opposes=("BF-I06",),
        addresses=("ecoli_adhesion", "BF-M07", "BF-M18"),
        maturity="formation_only",
        caution="Not proof of mature intestinal removal.",
    ),
    Intervention(
        id="BF-I06",
        compound="Berberine (adaptation counterevidence)",
        source_ids=(),
        doi="10.3389/fcimb.2025.1565714",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="unfavourable",
        result=(
            "Repeated half-MIC exposure selected a >32-fold MIC increase with "
            "INCREASED biofilm and csgD expression. A csgD-overexpression "
            "experiment supports the mechanism."
        ),
        lab_exposure="repeated half-MIC passage",
        delivery="plausible_unverified",
        delivery_note="Prolonged low-level exposure is the condition that produced harm.",
        section="reduce_matched_concerning_mechanism",
        opposes=("BF-I05",),
        addresses=("ecoli_adhesion", "BF-M01"),
        caution=(
            "Mandatory companion to any E. coli berberine card. It does not mean "
            "berberine never works; it means prolonged sub-inhibitory exposure "
            "carries an adaptation risk."
        ),
    ),
    Intervention(
        id="BF-I07",
        compound="Berberine with vancomycin",
        source_ids=(),
        doi="10.1007/s10096-020-03857-0",
        study_type="in_vitro",
        organisms=("Clostridioides difficile",),
        strain="12 strains",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="mixed",
        result=(
            "Planktonic antibacterial and synergistic effects, but half-MIC "
            "berberine or the combination ENHANCED biofilm in particular strains."
        ),
        lab_exposure="half-MIC and combination",
        delivery="plausible_unverified",
        delivery_note="",
        section="clarify_or_treat_established_condition",
        addresses=("cdifficile", "BF-M08"),
        caution=(
            "Planktonic efficacy and biofilm direction are separate outcomes. An "
            "antibacterial hit is not an antibiofilm hit."
        ),
    ),
    Intervention(
        id="BF-I08",
        compound="Berberine",
        source_ids=(),
        doi="10.2147/DDDT.S230857",
        study_type="in_vitro",
        organisms=("Candida",),
        strain="clinical isolates",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result=(
            "Candida isolates assessed separately in planktonic and biofilm "
            "conditions, with antifungal and biofilm effects."
        ),
        lab_exposure="isolate-specific MIC range",
        delivery="plausible_unverified",
        delivery_note="",
        section="clarify_or_treat_established_condition",
        addresses=("candida", "BF-M09"),
        caution=(
            "Match the exact Candida species and assay. Intestinal colonisation is "
            "not candidiasis, and stool DNA does not justify treatment."
        ),
    ),
    Intervention(
        id="BF-I09",
        compound="Carvacrol, thymol and geraniol",
        source_ids=(),
        doi="10.3390/antibiotics11020147",
        study_type="in_vitro",
        organisms=("Klebsiella pneumoniae",),
        strain="uropathogenic NDM-1-producing",
        model="extraintestinal",
        endpoint="formation",
        direction="favourable",
        result=(
            "15 essential-oil components compared; thymol, carvacrol and geraniol "
            "performed best for antibacterial and biofilm activity."
        ),
        lab_exposure="component-specific MIC panel",
        delivery="plausible_unverified",
        delivery_note=(
            "A concentrated purified compound differs from whole oregano or thyme "
            "oil. No validated human gut MBEC or local exposure."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("klebsiella_adhesion", "BF-M05"),
        maturity="formation_only",
        safety=(
            "Essential oils can irritate mucosa. Products manufactured for "
            "fragrance or topical use are not for oral use."
        ),
    ),
    Intervention(
        id="BF-I10",
        compound="Carvacrol, alone and with cefixime",
        source_ids=(),
        doi="10.1186/s12866-023-02797-x",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result="Antibacterial and antibiofilm effects, alone and in combination.",
        lab_exposure="MIC and sub-MIC",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        addresses=("ecoli_adhesion",),
        caution="The antibiotic combination belongs in clinician-review context.",
    ),
    Intervention(
        id="BF-I11",
        compound="Thyme and oregano oils, polymicrobial",
        source_ids=(),
        doi="10.1007/s10482-026-02284-z",
        study_type="in_vitro",
        organisms=("Klebsiella pneumoniae", "Acinetobacter baumannii"),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result="Oils alone and in combination against a polymicrobial biofilm in vitro.",
        lab_exposure="oil dilutions, in vitro",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        addresses=("klebsiella_adhesion",),
        caution=(
            "Broadens from single-species evidence but is not intact human "
            "gut-community evidence. Import full methods before reusing effects."
        ),
    ),
    Intervention(
        id="BF-I12",
        compound="Thymoquinone",
        source_ids=(),
        doi="10.1007/s12088-024-01231-8",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="carbapenem-resistant uropathogenic",
        model="extraintestinal",
        endpoint="eradication",
        direction="favourable",
        result=(
            "Formation inhibition and eradication assays. 256 ug/mL antibacterial; "
            "128 ug/mL reduced motility; changes in blaKPC, efflux and motility "
            "gene expression."
        ),
        lab_exposure="128-256 ug/mL",
        delivery="plausible_unverified",
        delivery_note=(
            "Purified thymoquinone is not equivalent to an arbitrary black-seed-oil "
            "dose or content. Active local concentration unknown."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("ecoli_adhesion",),
        maturity="both_formation_and_established",
    ),
    Intervention(
        id="BF-I13",
        compound="Carvacrol-thymoquinone nanocarrier",
        source_ids=(),
        doi="10.1016/j.diagmicrobio.2024.116606",
        study_type="in_vitro",
        organisms=("Candida albicans", "Candida glabrata"),
        strain="vaginal isolates",
        model="extraintestinal",
        endpoint="formation",
        direction="favourable",
        result="Formation prevention at half-MIC and MIC using 50 nm carriers.",
        lab_exposure="half-MIC and MIC, 50 nm carriers",
        delivery="route_site_mismatch",
        delivery_note="Vaginal context, not human gut.",
        section="investigational_precision",
        addresses=("candida", "BF-M09"),
        maturity="formation_only",
    ),
    Intervention(
        id="BF-I14",
        compound="N-acetylcysteine, metformin, secnidazole",
        source_ids=(),
        doi="10.1186/s12866-023-02969-9",
        study_type="in_vitro_and_animal",
        organisms=("Klebsiella pneumoniae",),
        strain="",
        model="animal_gut",
        endpoint="formation",
        direction="favourable",
        result=(
            "At one-eighth MIC each reduced formation and virulence expression. "
            "Mouse protection was tested after experimental manipulation, not as a "
            "human gut biofilm trial."
        ),
        lab_exposure="one-eighth MIC",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        addresses=("klebsiella_adhesion", "BF-M05"),
        caution="Stool genes do not establish a need for any of these drugs.",
    ),
    Intervention(
        id="BF-I15",
        compound="N-acetylcysteine (pH-qualified)",
        source_ids=(),
        doi="10.3390/microorganisms14020512",
        study_type="in_vitro",
        organisms=("Klebsiella pneumoniae",),
        strain="34 strains for susceptibility, selected strains for biofilms",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="mixed",
        result=(
            "Acidic NAC with polymyxin B was additive; NEUTRALISED NAC was much "
            "weaker, with significant activity generally only at 32 mg/mL."
        ),
        lab_exposure="up to 32 mg/mL; activity depends on medium pH",
        delivery="route_site_mismatch",
        delivery_note=(
            "Colonic and host buffering can erase activity that depended on "
            "acidification. Acidifying a person's gut to laboratory pH is not a "
            "proposal."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("klebsiella_adhesion",),
        maturity="formation_only",
        caution=(
            "The pH qualifier is the finding. The study's word 'disrupts' does not "
            "override its formation-assay design."
        ),
    ),
    Intervention(
        id="BF-I16",
        compound="N-acetylcysteine before culture-guided antibiotics",
        source_ids=(),
        doi="10.1016/j.cgh.2010.05.006",
        study_type="randomised_controlled_trial",
        organisms=("Helicobacter pylori",),
        strain="",
        model="human_gut_controlled",
        endpoint="eradication",
        direction="favourable",
        result=(
            "40 refractory patients with at least four prior failures. Eradication "
            "13/20 with NAC versus 4/20 without. Gastric biofilm rationale with "
            "biopsy evaluation."
        ),
        lab_exposure="oral NAC pretreatment, clinical regimen",
        delivery="target_exposure_measured",
        delivery_note="Gastric site, where oral NAC does reach.",
        section="clarify_or_treat_established_condition",
        opposes=("BF-I17",),
        addresses=("hpylori",),
        caution=(
            "A gastric H. pylori salvage strategy in refractory patients. Not "
            "NAC-alone and not a colonic dysbiosis treatment."
        ),
    ),
    Intervention(
        id="BF-I17",
        compound="N-acetylcysteine as first-line adjunct",
        source_ids=(),
        doi="10.1177/1756284820927306",
        study_type="randomised_controlled_trial",
        organisms=("Helicobacter pylori",),
        strain="",
        model="human_gut_controlled",
        endpoint="eradication",
        direction="unfavourable",
        result=(
            "680 treatment-naive patients. Adjunct NAC gave ITT eradication 81.7% "
            "versus 84.3% - no improvement."
        ),
        lab_exposure="oral NAC adjunct, clinical regimen",
        delivery="target_exposure_measured",
        delivery_note="",
        section="clarify_or_treat_established_condition",
        opposes=("BF-I16",),
        addresses=("hpylori",),
        caution=(
            "Mandatory companion to BF-I16. The regimen and the prior-treatment "
            "context differ; both must be shown."
        ),
    ),
    Intervention(
        id="BF-I18",
        compound="Serrapeptase",
        source_ids=(),
        doi="10.3390/microorganisms13081875",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="ATCC 25922",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="favourable",
        result=(
            "Capsule-derived preparation tested for formation and for 24 h "
            "treatment of preformed biofilms, with biomass, viability and amyloid "
            "stains. Formation IC50 14.2 ng/mL, with lower curli-associated signal. "
            "Proposed direct CsgA/CsgB binding comes from docking, not biochemistry."
        ),
        lab_exposure="IC50 14.2 ng/mL",
        delivery="plausible_unverified",
        delivery_note=(
            "Oral enteric delivery is possible, but active enzyme reaching the "
            "human colon at an effective concentration was not established."
        ),
        section="reduce_matched_concerning_mechanism",
        addresses=("ecoli_adhesion", "curli", "BF-M01"),
        maturity="both_formation_and_established",
        caution=(
            "One strain. No dose conversion from ng/mL, and fibrinolytic activity "
            "is not antibiofilm action."
        ),
    ),
    Intervention(
        id="BF-I19",
        compound="Lactoferrin, lysozyme and dextranase",
        source_ids=(),
        doi="10.2436/20.1501.01.171",
        study_type="in_vitro",
        organisms=("Escherichia coli", "Klebsiella pneumoniae"),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="viable_burden",
        direction="mixed",
        result=(
            "Single-species MBEC-device biofilms with an ATP-based viable-biomass "
            "readout. Reductions varied and no agent destroyed both completely."
        ),
        lab_exposure="MBEC device, agent-specific",
        delivery="plausible_unverified",
        delivery_note="Protein digestion and intact local exposure both matter.",
        section="reduce_matched_concerning_mechanism",
        opposes=("BF-I20", "BF-I21"),
        addresses=("ecoli_adhesion", "klebsiella_adhesion"),
        caution="ATP reduction is not histologic colon clearance.",
    ),
    Intervention(
        id="BF-I20",
        compound="Lactoferrin (organism-opposite directions)",
        source_ids=(),
        doi="10.1007/s11259-026-11512-w",
        study_type="in_vitro",
        organisms=("Escherichia coli", "Klebsiella pneumoniae"),
        strain="24 E. coli and 20 K. pneumoniae bovine/environment isolates",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="mixed",
        result=(
            "At 200 and 1,000 ug/mL growth decreased. K. pneumoniae biofilm was "
            "INHIBITED but E. coli biofilm was ENHANCED."
        ),
        lab_exposure="200 and 1,000 ug/mL",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        opposes=("BF-I19",),
        addresses=("ecoli_adhesion", "klebsiella_adhesion"),
        caution=(
            "Mandatory organism-direction counterevidence. Opposite directions "
            "must not be averaged into universal benefit."
        ),
    ),
    Intervention(
        id="BF-I21",
        compound="Lactoferrin (commensal collateral effects)",
        source_ids=(),
        doi="10.1016/j.anaerobe.2020.102232",
        study_type="in_vitro",
        organisms=("Bacteroides fragilis", "Bacteroides thetaiotaomicron"),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="adhesion",
        direction="unfavourable",
        result=(
            "Formation and laminin binding inhibited at 12.5 ug/mL. Growth was not "
            "inhibited by the physiological concentration of 2 mg/mL. Antibiotic "
            "synergy was negative."
        ),
        lab_exposure="12.5 ug/mL for adhesion; 2 mg/mL physiological",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        opposes=("BF-I19",),
        addresses=("commensal_collateral",),
        collateral=(
            "Matrix and colonisation effects extend to protective commensals. This "
            "is the collateral-effect record for lactoferrin."
        ),
        caution="Antibiofilm activity is not selective for unwanted organisms.",
    ),
    Intervention(
        id="BF-I22",
        compound="BioDisrupt enzyme and botanical mixture",
        source_ids=(),
        doi="10.4014/jmb.2212.12010",
        study_type="in_vitro",
        organisms=(
            "Candida",
            "Staphylococcus",
            "Pseudomonas aeruginosa",
            "Borrelia burgdorferi",
        ),
        strain="Candida and two Staphylococcus strains",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="mixed",
        result=(
            "Established biofilms: Candida and two Staphylococcus strains improved. "
            "P. aeruginosa REVERSED at a high dose with increased biomass and "
            "activity. B. burgdorferi residual metabolic activity increased. "
            "Combination of enzymes, NAC, cranberry, berberine, rosemary and "
            "peppermint."
        ),
        lab_exposure="product at 1-40 mg/mL",
        delivery="plausible_unverified",
        delivery_note="No gut exposure demonstrated.",
        section="reduce_matched_concerning_mechanism",
        addresses=("candida", "mixed_product"),
        funding_conflict="Manufacturer-sponsored.",
        maturity="mature",
        caution=(
            "Evidence for this exact mixture only - not for all marketed 'biofilm "
            "breakers', not for the ingredients individually, and not for systemic "
            "eradication. The dose-direction reversal must display."
        ),
    ),
    Intervention(
        id="BF-I23",
        compound="Herbal protocol with NAC or enzyme-blend adjunct",
        source_ids=(),
        doi="10.7759/cureus.99116",
        study_type="retrospective_case_series",
        organisms=(),
        strain="",
        model="human_gut_controlled",
        endpoint="clinical_outcome_no_biofilm_endpoint",
        direction="unfavourable",
        result=(
            "13 patients: herbs alone (n=5) versus herbs plus adjunct (n=8; NAC "
            "n=2, enzyme blend n=6). SIBO eradication was 60% WITH adjunct versus "
            "100% control, not significant. IMO eradication was zero in both. "
            "Within-group gas changes only; no symptom or direct biofilm "
            "measurement."
        ),
        lab_exposure="clinical protocol, retrospective",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        addresses=("sibo", "mixed_product"),
        funding_conflict="Formulation-owner conflict disclosed.",
        caution=(
            "n=13. Within-group significance is not between-group superiority, and "
            "breath-gas change does not demonstrate biofilm removal. The adjunct "
            "arm did numerically worse."
        ),
    ),
    Intervention(
        id="BF-I24",
        compound="Saccharomyces boulardii CNCM I-745",
        source_ids=(),
        doi="10.3390/microorganisms10061082",
        study_type="in_vitro",
        organisms=("Clostridioides difficile",),
        strain="CNCM I-745 against three C. difficile strains",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result=(
            "Live yeast co-culture reduced formation, thickness and eDNA. "
            "Genetically close S. cerevisiae did NOT reproduce the effect, and "
            "direct contact appeared required."
        ),
        lab_exposure="live yeast co-culture",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("cdifficile", "BF-M08"),
        maturity="formation_only",
        safety=(
            "Yeast infection risk is relevant for severely immunocompromised, "
            "critically ill or central-line patients."
        ),
        caution=(
            "Exact strain. The S. cerevisiae non-result is why this evidence "
            "cannot be assigned to other strains."
        ),
    ),
    Intervention(
        id="BF-I25",
        compound="Polydextrose",
        source_ids=(),
        doi="10.1128/spectrum.01017-25",
        study_type="animal_and_in_vitro",
        organisms=("Klebsiella pneumoniae",),
        strain="",
        model="animal_gut",
        endpoint="organism_load",
        direction="favourable",
        result=(
            "Mouse prevention work associated polydextrose with reduced bacterial "
            "loads and changes in TamA-related adhesion and biofilm mechanisms."
        ),
        lab_exposure="dietary supplementation, mouse",
        delivery="plausible_unverified",
        delivery_note="A food-grade prebiotic; colonic delivery is inherent.",
        section="support_protective_ecology",
        addresses=("klebsiella_adhesion", "BF-M05"),
        caution="Prevention, not removal of established human biofilms.",
    ),
    Intervention(
        id="BF-I26",
        compound="Wheat fibre supporting a nine-species consortium",
        source_ids=("BF-S16",),
        doi="10.1039/D4FO01294A",
        study_type="in_vitro",
        organisms=(),
        strain="nine-species gut consortium",
        model="gut_relevant_in_vitro",
        endpoint="protective_function",
        direction="favourable",
        result=(
            "A stable biofilm formed on wheat fibre and withstood relevant "
            "environmental challenge."
        ),
        lab_exposure="in vitro consortium on fibre",
        delivery="plausible_unverified",
        delivery_note="Dietary fibre reaches the colon by definition.",
        section="support_protective_ecology",
        addresses=("low_pe_proxy", "BF-M12", "BF-M19"),
        caution=(
            "Supports the protective ecological hypothesis. More biofilm is not "
            "automatically bad, and this is not established human clinical benefit."
        ),
    ),
    Intervention(
        id="BF-I27",
        compound="EGCG and myricetin",
        source_ids=(),
        doi="10.1038/s41598-018-26748-z",
        study_type="in_vitro",
        organisms=("Escherichia coli",),
        strain="K-12",
        model="gut_relevant_in_vitro",
        endpoint="formation",
        direction="favourable",
        result="Curli-related formation and regulation assays; EGCG was active.",
        lab_exposure="in vitro, K-12",
        delivery="plausible_unverified",
        delivery_note="A beverage and a concentrated extract are different interventions.",
        section="reduce_matched_concerning_mechanism",
        addresses=("curli", "ecoli_adhesion", "BF-M01"),
        maturity="formation_only",
        safety="Concentrated green tea extract has documented hepatotoxicity reports.",
    ),
    Intervention(
        id="BF-I28",
        compound="EGCG",
        source_ids=(),
        doi="10.3389/fcimb.2024.1432111",
        study_type="in_vitro_and_animal",
        organisms=("Salmonella",),
        strain="",
        model="animal_gut",
        endpoint="organism_load",
        direction="favourable",
        result="Laboratory quorum and biofilm inhibition plus a mouse enteritis model.",
        lab_exposure="in vitro and mouse dietary",
        delivery="plausible_unverified",
        delivery_note="",
        section="reduce_matched_concerning_mechanism",
        addresses=("salmonella",),
        caution="Organism-specific. Not proof that green tea removes human intestinal biofilms.",
    ),
    Intervention(
        id="BF-I29",
        compound="D-amino acids",
        source_ids=(),
        doi="10.1126/science.1188628",
        study_type="in_vitro",
        organisms=("Bacillus subtilis",),
        strain="laboratory strains",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="favourable",
        result="Matrix-disassembly and formation experiments in defined organisms.",
        lab_exposure="millimolar D-amino acids",
        delivery="plausible_unverified",
        delivery_note="",
        section="investigational_precision",
        opposes=("BF-I30",),
        addresses=("matrix_generic",),
        caution="Link the reproduction findings in BF-I30.",
    ),
    Intervention(
        id="BF-I30",
        compound="D-amino acids (counterevidence and reproduction)",
        source_ids=(),
        doi="10.1128/JB.00975-13",
        study_type="in_vitro",
        organisms=("Bacillus subtilis", "Staphylococcus"),
        strain="",
        model="gut_relevant_in_vitro",
        endpoint="matrix_biomass",
        direction="unfavourable",
        result=(
            "A laboratory-strain defect explained much of the apparent "
            "susceptibility, and a separate 2015 study did not reproduce "
            "inhibition in the tested Bacillus and Staphylococcus contexts."
        ),
        lab_exposure="as in the original reports",
        delivery="plausible_unverified",
        delivery_note="",
        section="investigational_precision",
        opposes=("BF-I29",),
        addresses=("matrix_generic",),
        caution=(
            "Preserve isomer, organism, strain and exposure. No universal "
            "D-amino-acid efficacy claim survives this pair."
        ),
    ),
    Intervention(
        id="BF-I31",
        compound="Phage consortium",
        source_ids=(),
        doi="10.1016/j.cell.2022.07.003",
        study_type="animal_and_human_safety",
        organisms=("Klebsiella pneumoniae",),
        strain="IBD-associated",
        model="animal_gut",
        endpoint="organism_load",
        direction="favourable",
        result=(
            "Suppressed IBD-associated K. pneumoniae and reduced inflammation in "
            "mice. Artificial-gut and healthy-volunteer work supported delivery "
            "and safety."
        ),
        lab_exposure="phage cocktail, oral",
        delivery="indirect_delivery_demonstrated",
        delivery_note="Delivery and safety were studied in humans; efficacy was not.",
        section="investigational_precision",
        addresses=("klebsiella_adhesion", "BF-M05"),
        caution=(
            "Pathobiont suppression is not a measured biofilm-eradication endpoint. "
            "Host range requires separate evidence."
        ),
    ),
    Intervention(
        id="BF-I32",
        compound="Targeted Klebsiella phages",
        source_ids=(),
        doi="10.1038/s41467-023-39029-9",
        study_type="animal",
        organisms=("Klebsiella pneumoniae",),
        strain="PSC-associated and high-alcohol-producing models",
        model="animal_gut",
        endpoint="organism_load",
        direction="favourable",
        result=(
            "Targeted gut Klebsiella reduction improved experimental hepatobiliary "
            "and metabolic disease outcomes."
        ),
        lab_exposure="phage, mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="investigational_precision",
        addresses=("klebsiella_adhesion",),
        caution="No biofilm endpoint was assessed in these studies.",
    ),
    Intervention(
        id="BF-I33",
        compound="Garlic oil, peppermint oil, curcumin, cinnamaldehyde",
        source_ids=(),
        doi="10.3390/pharmaceutics16040478",
        study_type="in_vitro",
        organisms=("Clostridioides difficile",),
        strain="strain-dependent",
        model="gut_relevant_in_vitro",
        endpoint="adhesion",
        direction="mixed",
        result=(
            "Strain-dependent antibacterial and adhesion changes. Biofilm "
            "architecture could also become DENSER and THICKER."
        ),
        lab_exposure="compound-specific in vitro",
        delivery="plausible_unverified",
        delivery_note="",
        section="clarify_or_treat_established_condition",
        addresses=("cdifficile", "BF-M08"),
        caution=(
            "Growth inhibition cannot be converted into uniformly beneficial "
            "antibiofilm effects."
        ),
    ),
    Intervention(
        id="BF-I34",
        compound="L. reuteri ATCC 23272 in microsphere formulation",
        source_ids=("BF-S15",),
        doi="10.3389/fmicb.2017.00489",
        study_type="in_vitro",
        organisms=("Limosilactobacillus reuteri",),
        strain="ATCC 23272",
        model="gut_relevant_in_vitro",
        endpoint="protective_function",
        direction="favourable",
        result="Improved probiotic properties in the biofilm state, with GtfW implicated.",
        lab_exposure="microsphere formulation, in vitro",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("BF-M10", "protective_support"),
        caution=(
            "The ATCC 23272 and GtfW context is not a generic probiotic-species "
            "effect, and the microsphere formulation is not naturally present in a "
            "donor sample."
        ),
    ),
    Intervention(
        id="BF-I35",
        compound="SB-121 (L. reuteri with dextran microparticles and maltose)",
        source_ids=(),
        doi="10.1038/s41598-023-30909-0",
        study_type="randomised_controlled_trial",
        organisms=("Limosilactobacillus reuteri",),
        strain="SB-121 formulation",
        model="human_gut_controlled",
        endpoint="clinical_outcome_no_biofilm_endpoint",
        direction="neutral",
        result=(
            "Small randomised placebo-controlled crossover study in 15 participants "
            "with autism, reporting safety, tolerability and exploratory outcomes."
        ),
        lab_exposure="oral SB-121, clinical",
        delivery="target_exposure_measured",
        delivery_note="An oral formulation studied in humans.",
        section="support_protective_ecology",
        addresses=("BF-M10", "protective_support"),
        caution=(
            "Relevant evidence that biofilm-oriented supportive formulations have "
            "reached human research. Not validation of a stool protective score and "
            "not proof of broad clinical efficacy."
        ),
    ),
    Intervention(
        id="BF-I36",
        compound="Host-glycan support of F. prausnitzii",
        source_ids=("BF-S14",),
        doi="10.1080/19490976.2026.2665870",
        study_type="animal",
        organisms=("Faecalibacterium prausnitzii",),
        strain="",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result=(
            "Protective mechanism in experimental colitis with correlations in "
            "human tissue."
        ),
        lab_exposure="glycan intervention, mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("BF-M11", "low_pe_proxy", "protective_support"),
        caution=(
            "Not a validated oral treatment and not a unique DNA biomarker."
        ),
    ),
    Intervention(
        id="BF-I37",
        compound="Serrapeptase (staphylococcal and pseudomonal)",
        source_ids=(),
        doi="10.1007/s00253-022-12356-5",
        study_type="in_vitro",
        organisms=("Staphylococcus aureus", "Pseudomonas aeruginosa"),
        strain="strain-specific",
        model="extraintestinal",
        endpoint="matrix_biomass",
        direction="favourable",
        result="Strain-specific laboratory enzyme activity against these organisms.",
        lab_exposure="enzyme units, in vitro",
        delivery="plausible_unverified",
        delivery_note="",
        section="investigational_precision",
        addresses=("staph", "BF-M14"),
        caution=(
            "Separate from the E. coli study BF-I18. No automatic systemic or gut "
            "clearance claim."
        ),
    ),
    Intervention(
        id="BF-I38",
        compound="Nattokinase",
        source_ids=(),
        doi="10.3390/pathogens13040286",
        study_type="in_vitro",
        organisms=("Streptococcus mutans",),
        strain="46 dental-plaque isolates",
        model="extraintestinal",
        endpoint="matrix_biomass",
        direction="favourable",
        result="Organism-specific activity in 46 dental-plaque isolates.",
        lab_exposure="enzyme units, in vitro",
        delivery="route_site_mismatch",
        delivery_note="Oral cavity, not colon.",
        section="investigational_precision",
        addresses=("fibrinolytic_enzymes",),
        caution="Not evidence that oral nattokinase clears colonic biofilms.",
    ),
    Intervention(
        id="BF-I39",
        compound="Lumbrokinase with emodin in a biomaterial",
        source_ids=(),
        doi="10.3390/bioengineering10080906",
        study_type="in_vitro",
        organisms=(),
        strain="",
        model="extraintestinal",
        endpoint="matrix_biomass",
        direction="favourable",
        result="Local combination and material context.",
        lab_exposure="biomaterial, local",
        delivery="route_site_mismatch",
        delivery_note="A local material, not an oral enzyme.",
        section="investigational_precision",
        addresses=("fibrinolytic_enzymes",),
        caution="The result belongs to the combination material, not to oral lumbrokinase.",
    ),
    Intervention(
        id="BF-I40",
        compound="Anti-DNABII antibody and host HMGB1 approaches",
        source_ids=(),
        doi="10.3389/fmicb.2023.1202215",
        study_type="animal",
        organisms=(),
        strain="",
        model="extraintestinal",
        endpoint="matrix_biomass",
        direction="favourable",
        result="Investigational matrix-directed approaches with model clearance.",
        lab_exposure="antibody, animal models",
        delivery="route_site_mismatch",
        delivery_note="Clinical availability and route are separate from mechanistic relevance.",
        section="investigational_precision",
        addresses=("BF-M07", "matrix_generic"),
        caution=(
            "Increased susceptibility after matrix release is not proof that "
            "unsupervised dispersal is harmless."
        ),
    ),
    Intervention(
        id="BF-I41",
        compound="Biofilm-targeted antibody versus antibiotic, microbiome effects",
        source_ids=(),
        doi="10.1038/s41522-024-00481-0",
        study_type="animal",
        organisms=(),
        strain="",
        model="extraintestinal",
        endpoint="organism_load",
        direction="favourable",
        result=(
            "Duff et al.: oral and middle-ear antibiotic versus targeted-antibody "
            "microbiome effects in chinchillas."
        ),
        lab_exposure="antibody versus antibiotic, animal",
        delivery="route_site_mismatch",
        delivery_note="",
        section="investigational_precision",
        addresses=("commensal_collateral",),
        collateral="The comparison's point is differential effect on commensals.",
        caution=(
            "An animal model with local treatment. Not a human gut therapy and not "
            "proof that all antibody approaches spare every commensal."
        ),
    ),
    Intervention(
        id="BF-I42",
        compound="Grape-pomace polyphenols",
        source_ids=("BF-S19",),
        doi="10.1038/s41522-025-00716-8",
        study_type="human_derived_ex_vivo",
        organisms=(),
        strain="",
        model="human_derived_gut_community",
        endpoint="harmful_behaviour",
        direction="favourable",
        result=(
            "In human-derived prefrail communities the tested extract REDUCED "
            "dispersal and IL-6-inducing effects while biomass INCREASED. "
            "Epithelial adhesion, IL-8 and TNF did not significantly improve. "
            "Small donor subset."
        ),
        lab_exposure="extract concentration in ex-vivo community culture",
        delivery="plausible_unverified",
        delivery_note=(
            "Oral exposure at the tested community conditions is unestablished, "
            "and an ordinary grape food is not this extract."
        ),
        section="support_protective_ecology",
        addresses=("measured_dispersal", "BF-M21", "community_stabilisation"),
        caution=(
            "A community-stabilisation card: the favourable endpoint is behaviour, "
            "not biomass. Biomass went up and that is not a failure here, which is "
            "precisely why biomass alone is not the success criterion."
        ),
    ),
    Intervention(
        id="BF-I43",
        compound="Pomegranate extract",
        source_ids=("BF-S24",),
        doi="10.3390/nu15071771",
        study_type="animal",
        organisms=("Citrobacter freundii",),
        strain="",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result=(
            "Improved experimental-colitis recovery and barrier outcomes. Broad "
            "Enterobacterales biofilm inhibition was NOT demonstrated, with only a "
            "slight effect on C. freundii."
        ),
        lab_exposure="extract, mouse dietary",
        delivery="plausible_unverified",
        delivery_note="Product composition and gut delivery remain uncertain.",
        section="support_protective_ecology",
        addresses=("barrier_support", "low_pe_proxy"),
        caution=(
            "A barrier-support candidate. It must not be described as proven "
            "Klebsiella or E. coli biofilm clearance - the study does not show that."
        ),
    ),
    Intervention(
        id="BF-I44",
        compound="L. johnsonii studied isolate",
        source_ids=("BF-S25",),
        doi="10.3389/fimmu.2026.1749001",
        study_type="in_vitro_and_animal",
        organisms=("Lactobacillus johnsonii", "Escherichia coli"),
        strain="the studied isolate; PV739486 is a 16S sequence only",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result=(
            "Direct EPEC biofilm experiments plus improved pathogen and "
            "inflammation outcomes in a C. rodentium mouse model."
        ),
        lab_exposure="live isolate, in vitro and mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("protective_support", "BF-M24", "ecoli_adhesion"),
        caution=(
            "Do not substitute a commercial L. johnsonii strain without identity "
            "evidence. Public 16S data do not provide an exact strain genome, and "
            "putatively annotated metabolites are not proven active ingredients."
        ),
    ),
    Intervention(
        id="BF-I45",
        compound="L. plantarum CIRM653",
        source_ids=("BF-S23",),
        doi="10.3920/bm2017.0002",
        study_type="in_vitro_and_animal",
        organisms=("Lactiplantibacillus plantarum", "Klebsiella pneumoniae"),
        strain="CIRM653",
        model="animal_gut",
        endpoint="viable_burden",
        direction="unfavourable",
        result=(
            "Reduced K. pneumoniae biofilm in vitro, but treated mice MAINTAINED "
            "intestinal carriage while untreated controls declined."
        ),
        lab_exposure="live strain, in vitro and mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("klebsiella_adhesion", "protective_support", "BF-M24"),
        caution=(
            "The essential regression case: a favourable in-vitro biomass result "
            "with an unfavourable in-vivo carriage result. It must never be ranked "
            "as a selected decolonisation treatment."
        ),
    ),
    Intervention(
        id="BF-I46",
        compound="Resveratrol with L. paracasei ATCC334",
        source_ids=("BF-S26",),
        doi="10.3390/ijms21155423",
        study_type="in_vitro",
        organisms=("Lacticaseibacillus paracasei",),
        strain="ATCC334",
        model="gut_relevant_in_vitro",
        endpoint="adhesion",
        direction="favourable",
        result="Strain-dependent promotion of adhesion and biofilm formation in vitro.",
        lab_exposure="resveratrol in culture",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("protective_support", "BF-M13", "low_pe_proxy"),
        caution=(
            "A supportive formulation hypothesis. Not proven human engraftment, not "
            "a beneficial-biofilm quantity, and not a reason to take concentrated "
            "resveratrol."
        ),
    ),
    Intervention(
        id="BF-I47",
        compound="Biofilm-state L. plantarum Y42 and its matrix products",
        source_ids=("BF-S27",),
        doi="10.1021/acs.jafc.4c00460",
        study_type="animal",
        organisms=("Lactiplantibacillus plantarum",),
        strain="Y42",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result="Barrier and immune outcomes improved in a mouse challenge model.",
        lab_exposure="biofilm-state cells and components, mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("protective_support", "BF-M13", "barrier_support"),
        caution=(
            "Exact strain, state and components matter. Ordinary planktonic "
            "products and generic EPS genes are not equivalent."
        ),
    ),
    Intervention(
        id="BF-I48",
        compound="Biofilm-state Lactiplantibacillus LR-1",
        source_ids=("BF-S28",),
        doi="10.1039/d3fo02733c",
        study_type="animal",
        organisms=("Lactiplantibacillus plantarum",),
        strain="LR-1",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result=(
            "Biofilm-state cells outperformed planktonic cells on selected colitis "
            "outcomes."
        ),
        lab_exposure="biofilm-state cells, mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("protective_support", "BF-M13"),
        caution="Experimental protective support, with no claimed human treatment effect.",
    ),
    Intervention(
        id="BF-I49",
        compound="Multi-strain restorative consortium",
        source_ids=("BF-S29",),
        doi="10.1007/s11274-026-05071-0",
        study_type="animal",
        organisms=(),
        strain="the administered consortium",
        model="animal_gut",
        endpoint="protective_function",
        direction="favourable",
        result="Mouse recovery after antibiotics.",
        lab_exposure="administered consortium, mouse",
        delivery="plausible_unverified",
        delivery_note="",
        section="support_protective_ecology",
        addresses=("protective_support", "BF-M24"),
        caution=(
            "Full tested composition and methods must be imported before any exact "
            "product matching. Co-detection of these taxa in stool or donor "
            "material does not recreate the consortium."
        ),
    ),
    Intervention(
        id="BF-I50",
        compound="Colon-targeted anti-E. coli biofilm prodrug",
        source_ids=("BF-S30",),
        doi="10.1021/acsmedchemlett.5c00675",
        study_type="animal",
        organisms=("Escherichia coli",),
        strain="",
        model="animal_gut",
        endpoint="matrix_biomass",
        direction="favourable",
        result="An experimental prodrug improved colonic delivery and mouse-colitis outcomes.",
        lab_exposure="prodrug, mouse",
        delivery="indirect_delivery_demonstrated",
        delivery_note="Colonic delivery was the point of the design and was demonstrated in mice.",
        section="investigational_precision",
        addresses=("ecoli_adhesion", "BF-M18"),
        caution=(
            "Investigational. No validated human regimen and no generally available "
            "supplement. No synthesis or unsupervised experimental treatment."
        ),
    ),
)

INTERVENTIONS: Final[dict[str, Intervention]] = {i.id: i for i in _ALL}

#: Finding-to-action routing from section 10. Keys are the ``addresses``
#: tags a finding emits; values are the interventions that must be
#: retrievable for it.
ROUTING: Final[dict[str, str]] = {
    "ecoli_adhesion": (
        "Confirm carrier and pathotype where clinically relevant. Follow the raw "
        "signal, carrier abundance and symptoms separately. None of these options "
        "is promised to remove a human intestinal biofilm."
    ),
    "klebsiella_adhesion": (
        "Distinguish colonisation from infection. Actual viable burden and a "
        "relevant clinical assessment are more informative than a generic matrix "
        "reduction. BF-I45's counterexample belongs in any probiotic "
        "decolonisation discussion."
    ),
    "high_hc_proxy": (
        "Explain the two contributing taxa. Gene or carrier resolution, or an "
        "independently measured phenotype, would clarify relevance. Proxy "
        "elevation alone is not an eradication instruction."
    ),
    "low_pe_proxy": (
        "Track the measured ecological pattern alongside tolerance and clinical "
        "outcomes. This does not diagnose absent protective biofilm and does not "
        "promise stable engraftment."
    ),
    "measured_dispersal": (
        "Repeat the same measured endpoint under comparable conditions, not just "
        "total biomass."
    ),
    "protective_support": (
        "Improvement is not defined as the disappearance of every biofilm-related "
        "gene. Preservation-oriented options belong alongside target-specific ones."
    ),
    "matrix_generic": (
        "Explain the context and show the measured inventory. No organism-"
        "eradication plan follows from generic matrix or housekeeping genes."
    ),
}


def get(intervention_id: str) -> Intervention:
    try:
        return INTERVENTIONS[intervention_id]
    except KeyError:
        raise KeyError(
            f"unknown intervention {intervention_id!r}; known IDs are BF-I01..BF-I50"
        ) from None


def with_opposing(ids: tuple[str, ...]) -> tuple[str, ...]:
    """Close a set of intervention IDs under the ``opposes`` relation.

    Retrieving BF-I05 always also retrieves BF-I06. This is the mechanism
    that makes "strong contradictory results appear on the same card, not
    hidden in a bibliography" structural rather than a matter of whoever
    writes the report template remembering to do it.
    """
    out: list[str] = []
    seen: set[str] = set()
    for iid in ids:
        for candidate in (iid, *get(iid).opposes):
            if candidate not in seen:
                seen.add(candidate)
                out.append(candidate)
    return tuple(out)


def retrieve(tags: tuple[str, ...]) -> tuple[Intervention, ...]:
    """Every intervention matching any tag, closed under opposing evidence
    and sorted by the deterministic lexicographic key from section 10.

    The key is, in order: endpoint directness, model applicability,
    delivery evidence, whether the record is adverse (adverse records sort
    after favourable ones but are never dropped), then the stable ID.
    Nothing in this ordering is a claim of net benefit.
    """
    wanted = set(tags)
    matched = tuple(
        i.id for i in _ALL if wanted & set(i.addresses)
    )
    closed = with_opposing(matched)
    records = [get(i) for i in closed]
    records.sort(
        key=lambda i: (
            i.endpoint_rank,
            i.model_rank,
            i.delivery_rank,
            i.is_adverse,
            i.id,
        )
    )
    return tuple(records)


def census() -> dict[str, object]:
    """Counts by section, direction and endpoint for the report."""
    by_section: dict[str, int] = {}
    by_direction: dict[str, int] = {}
    for i in _ALL:
        by_section[i.section] = by_section.get(i.section, 0) + 1
        by_direction[i.direction] = by_direction.get(i.direction, 0) + 1
    linked = [(i.id, list(i.opposes)) for i in _ALL if i.opposes]
    return {
        "n_interventions": len(_ALL),
        "by_section": by_section,
        "by_direction": by_direction,
        "opposing_pairs": linked,
        "n_with_opposing": len(linked),
    }
