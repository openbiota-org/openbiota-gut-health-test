"""The 32 primary sources behind the biofilm module (BF-S01 .. BF-S32).

Every scoring rule, module direction and intervention card in this package
traces back to a record here. The point of storing them as data rather than
as prose in a report template is that a claim and its limits travel
together: a module cannot cite BF-S01 for a stool measurement, because
BF-S01's ``specimen`` says biopsy, and the code checks.

Three fields do the load-bearing work.

``result`` is what the study actually observed, with its real denominators.
Where a paper is routinely miscited by conflating two numbers, both numbers
are here and :mod:`openbiota.biofilm.selftest` asserts they stay distinct.

``supports`` is the narrow claim the source can carry. ``cannot_support``
is the claim it is most likely to be stretched into. The second field is
not decoration - the acceptance suite reads it, and a module that asserts
something listed in its source's ``cannot_support`` fails to build.

``evidence_class`` places the observation on the transport ladder from
human association down to in-silico hypothesis. It sets the ceiling on what
any module resting on this source may claim, independent of how good the
molecular detection was.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

#: Where the observation sits on the transport ladder. Ordered: earlier
#: entries transport to a human gut claim more readily than later ones.
EVIDENCE_CLASSES: Final = (
    "human_association",
    "human_derived_ex_vivo",
    "animal_mechanism",
    "in_vitro_mechanism",
    "in_silico_hypothesis",
    "mixed",
    "review",
    "resource",
)

#: Specimen the observation was actually made on. A stool-DNA module may not
#: cite a biopsy-only source as if it measured stool.
SPECIMENS: Final = (
    "stool",
    "biopsy",
    "stool_and_biopsy",
    "isolate",
    "laboratory_community",
    "animal_tissue",
    "human_tissue",
    "not_applicable",
)


@dataclass(frozen=True, slots=True)
class Source:
    """One primary source, with the boundary of what it can be cited for."""

    id: str
    citation: str
    doi: str
    year: int
    #: What the study actually found, with real denominators.
    result: str
    #: The narrow claim this source may be cited for.
    supports: str
    #: The claim it must never be stretched into.
    cannot_support: str
    evidence_class: str
    specimen: str
    #: Organisms actually studied, empty when not organism-specific.
    organisms: tuple[str, ...] = ()
    #: Denominators that are commonly conflated, kept apart on purpose.
    denominators: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.evidence_class not in EVIDENCE_CLASSES:
            raise ValueError(f"{self.id}: unknown evidence_class {self.evidence_class!r}")
        if self.specimen not in SPECIMENS:
            raise ValueError(f"{self.id}: unknown specimen {self.specimen!r}")

    @property
    def url(self) -> str:
        return f"https://doi.org/{self.doi}"

    @property
    def human_gut_applicable(self) -> bool:
        """True when the observation was made in a human gut context.

        Used by the finding-ranking key in section 8.8, which places human
        gut evidence above human-derived community evidence above animal
        gut evidence. An in-vitro isolate result is not human gut evidence
        however strong the effect was.
        """
        return self.evidence_class in {"human_association", "human_derived_ex_vivo"} and (
            self.specimen in {"stool", "biopsy", "stool_and_biopsy", "human_tissue"}
        )


def _s(*args: object, **kwargs: object) -> Source:
    return Source(*args, **kwargs)  # type: ignore[arg-type]


SOURCES: Final[dict[str, Source]] = {
    s.id: s
    for s in (
        Source(
            id="BF-S01",
            citation="Baumgartner et al., Gastroenterology",
            doi="10.1053/j.gastro.2021.06.024",
            year=2021,
            result=(
                "Visible ileocolonic biofilms in 57% of IBS, 34% of UC and 6% of "
                "controls. 1,426 participants were screened and 1,112 met inclusion "
                "criteria; 212 were biofilm-positive overall. Molecular analyses used "
                "smaller subsets, not the full included cohort."
            ),
            supports=(
                "That human mucosal biofilms are real, endoscopically visible, and "
                "more common in IBS and UC than in controls."
            ),
            cannot_support=(
                "A validated stool-shotgun biofilm classifier, or a claim that the "
                "association is universal or causal."
            ),
            evidence_class="human_association",
            specimen="stool_and_biopsy",
            denominators={"screened": 1426, "included": 1112, "biofilm_positive": 212},
        ),
        Source(
            id="BF-S02",
            citation="Buret and Allain, Journal of Experimental Medicine",
            doi="10.1084/jem.20221743",
            year=2023,
            result=(
                "Synthesis of host regulation of gut biofilms and of the disruption "
                "of protective communities."
            ),
            supports=(
                "That dissolving every biofilm is not a coherent treatment "
                "objective, because some matrix communities are protective."
            ),
            cannot_support="Independently validated diagnostic weights of any kind.",
            evidence_class="review",
            specimen="not_applicable",
        ),
        Source(
            id="BF-S03",
            citation="Jandl et al., Clinical Microbiology Reviews",
            doi="10.1128/cmr.00133-23",
            year=2024,
            result="Review of intestinal biofilms, host defence and therapeutic opportunities.",
            supports="An evidence map for which mechanisms are worth measuring.",
            cannot_support=(
                "Independently validated diagnostic weights. Review statements do "
                "not become coefficients."
            ),
            evidence_class="review",
            specimen="not_applicable",
        ),
        Source(
            id="BF-S04",
            citation="Biophysical determinants of gut biofilms, Current Opinion in Biomedical Engineering",
            doi="10.1016/j.cobme.2021.100275",
            year=2021,
            result="Flow, mucus and spatial conditions shape where and how gut biofilms form.",
            supports=(
                "That architecture depends on physical conditions DNA does not "
                "observe."
            ),
            cannot_support=(
                "Reconstruction of thickness, attachment or location from sequence "
                "data."
            ),
            evidence_class="review",
            specimen="not_applicable",
        ),
        Source(
            id="BF-S05",
            citation="Hillman et al., Gastro Hep Advances",
            doi="10.1016/j.gastha.2025.100712",
            year=2025,
            result=(
                "Shotgun stool study of 26 bile-acid-diarrhoea cases and 21 controls. "
                "Higher biofilm-associated gene abundance and altered bile-acid "
                "metabolism in cases. MetaPhlAn 3.1 and HUMAnN 3.0.3. The abstract "
                "emphasises bssS, pgaA and pgaB; the results discuss the five-gene "
                "set bssS and pgaA/B/C/D."
            ),
            supports=(
                "An exploratory stool-DNA gene pattern associated with bile-acid "
                "diarrhoea in one small cross-sectional cohort."
            ),
            cannot_support=(
                "A bile-acid-diarrhoea diagnosis, a measurement of mucosal biofilm, "
                "attribution of the community genes to any one organism, or "
                "reuse of its fold differences as another pipeline's reference."
            ),
            evidence_class="human_association",
            specimen="stool",
            organisms=("Mediterraneibacter gnavus", "Escherichia coli"),
            denominators={"cases": 26, "controls": 21},
        ),
        Source(
            id="BF-S06",
            citation="Biofilm-associated protein and amyloid study, Nature Communications",
            doi="10.1038/s41467-024-48309-x",
            year=2024,
            result=(
                "Accessory biofilm-associated proteins (BAP) with amyloid-related "
                "experiments, stool molecular measurements, and a reanalysis of a "
                "Parkinson's metagenome cohort. Table S2 carries the authors' BAP "
                "reference identifiers."
            ),
            supports=(
                "A BAP genetic-potential inventory, and separately, measured-amyloid "
                "evidence where an actual assay was run."
            ),
            cannot_support=(
                "Brain amyloid, Parkinson's transmission, or harmful activity in "
                "every carrier of a BAP homolog."
            ),
            evidence_class="mixed",
            specimen="stool",
        ),
        Source(
            id="BF-S07",
            citation="BBSdb, Frontiers in Cellular and Infection Microbiology",
            doi="10.3389/fcimb.2024.1428784",
            year=2024,
            result=(
                "Curated database of biofilm-associated proteins with transcriptomic "
                "and proteomic evidence across organisms. Prose and tables report "
                "inconsistent record totals, so imported counts must be counted, not "
                "quoted."
            ),
            supports="Candidate protein and function annotations for discovery.",
            cannot_support=(
                "A person-level biofilm classifier, a good/bad call, or patient "
                "diagnostic accuracy derived from protein-classification performance."
            ),
            evidence_class="resource",
            specimen="not_applicable",
        ),
        Source(
            id="BF-S08",
            citation="Klebsiella clinical isolate biofilm study, npj Biofilms and Microbiomes",
            doi="10.1038/s41522-024-00629-y",
            year=2024,
            result=(
                "Substantial biofilm-phenotype variation across 100 clinical "
                "Klebsiella isolates carrying broadly similar adhesion gene content."
            ),
            supports=(
                "That adhesion gene presence does not predict the biofilm phenotype "
                "of the strain carrying it."
            ),
            cannot_support=(
                "A strain's biofilm phenotype inferred from species abundance plus "
                "one adhesion-gene hit."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="isolate",
            organisms=("Klebsiella pneumoniae",),
            denominators={"isolates": 100},
        ),
        Source(
            id="BF-S09",
            citation="Klebsiella biofilm evolution, Nature Communications",
            doi="10.1038/s41467-026-71505-w",
            year=2026,
            result=(
                "Background- and environment-dependent genomic changes during biofilm "
                "adaptation, in non-gut infection models."
            ),
            supports="Research context for curated variants.",
            cannot_support="Universal unvalidated gut-risk weights for those variants.",
            evidence_class="animal_mechanism",
            specimen="isolate",
            organisms=("Klebsiella pneumoniae",),
        ),
        Source(
            id="BF-S10",
            citation="Tomkovich et al., Journal of Clinical Investigation",
            doi="10.1172/JCI124196",
            year=2019,
            result=(
                "Human biofilm-positive mucosal communities promoted tumours when "
                "transferred into susceptible mice."
            ),
            supports="That the mechanistic concern about biofilm communities is real.",
            cannot_support=(
                "A donor-specific human transmission probability, or any numeric "
                "risk of transferring disease by transplant."
            ),
            evidence_class="animal_mechanism",
            specimen="human_tissue",
        ),
        Source(
            id="BF-S11",
            citation="Dejea et al., PNAS",
            doi="10.1073/pnas.1406199111",
            year=2014,
            result=(
                "Spatially organised mucosal communities associated with colorectal "
                "cancer, concentrated in the proximal colon."
            ),
            supports="That tissue location and community organisation matter.",
            cannot_support=(
                "A risk call from a common-genus list, without spatial or functional "
                "evidence."
            ),
            evidence_class="human_association",
            specimen="biopsy",
        ),
        Source(
            id="BF-S12",
            citation="Motta et al., Nature Communications",
            doi="10.1038/s41467-019-11140-w",
            year=2019,
            result=(
                "Physiological epithelial thrombin contributes to control of mucosal "
                "biofilms."
            ),
            supports=(
                "That host protease activity participates in normal biofilm control "
                "and is not inherently adverse."
            ),
            cannot_support="Any inference of host thrombin activity from DNA.",
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
        ),
        Source(
            id="BF-S13",
            citation="Excess thrombin and Crohn's disease, Gut Microbes",
            doi="10.1080/19490976.2026.2687903",
            year=2026,
            result=(
                "Excessive thrombin induced harmful community behaviour with limited "
                "accompanying taxonomic change."
            ),
            supports=(
                "That important state changes can occur without taxonomic change, so "
                "static gene counts miss them."
            ),
            cannot_support=(
                "A thrombin activity value inferred from sequence, or fabricated "
                "protein and RNA inputs."
            ),
            evidence_class="mixed",
            specimen="laboratory_community",
        ),
        Source(
            id="BF-S14",
            citation="F. prausnitzii and syndecan-1, Gut Microbes",
            doi="10.1080/19490976.2026.2665870",
            year=2026,
            result=(
                "Host glycans supported protective bacterial biofilm-related behaviour "
                "in experimental colitis, with correlations in human tissue."
            ),
            supports=(
                "That protective biofilm biology exists and deserves its own axis."
            ),
            cannot_support=(
                "A validated protective-biofilm score built from F. prausnitzii "
                "abundance, generic vitamin genes, or residual host SDC1 DNA."
            ),
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
            organisms=("Faecalibacterium prausnitzii",),
        ),
        Source(
            id="BF-S15",
            citation="L. reuteri microsphere biofilm study, Frontiers in Microbiology",
            doi="10.3389/fmicb.2017.00489",
            year=2017,
            result=(
                "A specified L. reuteri strain in a defined microsphere formulation "
                "showed improved probiotic properties in the biofilm state, with GtfW "
                "implicated."
            ),
            supports=(
                "A strain- and formulation-qualified protective mechanism card."
            ),
            cannot_support=(
                "That all L. reuteri, or all biofilms of it, carry the same benefit, "
                "or that the microsphere formulation exists in a stool sample."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="isolate",
            organisms=("Limosilactobacillus reuteri ATCC 23272",),
        ),
        Source(
            id="BF-S16",
            citation="Wheat-fibre community biofilm study, Food & Function",
            doi="10.1039/D4FO01294A",
            year=2024,
            result=(
                "A nine-species consortium formed a stable biofilm on wheat fibre in "
                "vitro and withstood relevant environmental challenge."
            ),
            supports=(
                "That fibre-associated biofilm formation can be supportive, so more "
                "matrix is not automatically harmful."
            ),
            cannot_support=(
                "An exclusive DNA signature of protective biofilm, or established "
                "human clinical benefit."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="laboratory_community",
        ),
        Source(
            id="BF-S17",
            citation="Staphylococcal fibrin model, Biofilm",
            doi="10.1016/j.bioflm.2025.100261",
            year=2025,
            result=(
                "Host fibrin recruitment into staphylococcal biofilms was specific to "
                "organism, strain and substrate."
            ),
            supports=(
                "A microbial fibrin-interaction capability card, kept separate from "
                "measured host fibrin."
            ),
            cannot_support=(
                "Measured host fibrin, systemic infection, or an intestinal fibrin "
                "biofilm."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="isolate",
            organisms=("Staphylococcus aureus",),
        ),
        Source(
            id="BF-S18",
            citation="Staphylococcal context-dependence study, Biofilm",
            doi="10.1016/j.bioflm.2026.100389",
            year=2026,
            result=(
                "The same genetic changes had different effects in surface biofilms "
                "versus abscess communities."
            ),
            supports=(
                "That a regulatory gene has no universal positive or negative weight "
                "across sites."
            ),
            cannot_support="A single directional weight for a regulator.",
            evidence_class="in_vitro_mechanism",
            specimen="isolate",
            organisms=("Staphylococcus aureus",),
        ),
        Source(
            id="BF-S19",
            citation="Prefrail human gut community biofilms, npj Biofilms and Microbiomes",
            doi="10.1038/s41522-025-00716-8",
            year=2025,
            result=(
                "42 donors: 13 adults aged 30-70, 15 robust adults over 70, and 14 "
                "prefrail adults over 70. Faeces-derived laboratory biofilms differed "
                "in dispersal and inflammatory effects. Grape polyphenols reduced "
                "selected adverse effects while increasing biomass. The transfer "
                "experiment used mouse microbiota, not a human treatment trial."
            ),
            supports=(
                "Community-behaviour measurements and stabilisation-oriented action "
                "candidates, as human-derived ex-vivo evidence."
            ),
            cannot_support=(
                "Observed in-vivo intestinal architecture, calibration for minors, or "
                "a human treatment result."
            ),
            evidence_class="human_derived_ex_vivo",
            specimen="laboratory_community",
            denominators={"donors": 42, "adults_30_70": 13, "robust_over_70": 15, "prefrail_over_70": 14},
        ),
        Source(
            id="BF-S20",
            citation="UC biofilms and oncotraits, Journal of Crohn's and Colitis",
            doi="10.1093/ecco-jcc/jjad092",
            year=2023,
            result=(
                "80 UC patients and 35 controls; 873 longitudinal biopsies and 265 "
                "shotgun-sequenced specimens. The sequencing was of BIOPSIES; stool "
                "oncotraits were measured by qPCR in a subset. Biofilm presence was "
                "not significantly associated with dysplasia (aOR 1.45, 95% CI "
                "0.63-3.40), whereas colibactin-marker carriage was. FadA showed an "
                "inverse association in this UC cohort."
            ),
            supports=(
                "That specific damage mechanisms outrank generic biofilm presence, "
                "and that context reverses direction."
            ),
            cannot_support=(
                "A universal rule from all biofilm, all Fusobacterium or all FadA. "
                "265 specimens are not 265 people and not stool WGS."
            ),
            evidence_class="human_association",
            specimen="biopsy",
            denominators={
                "uc_patients": 80,
                "controls": 35,
                "biopsies": 873,
                "sequenced_specimens": 265,
            },
        ),
        Source(
            id="BF-S21",
            citation="Cirrhosis virulence-factor study, Gut Microbes",
            doi="10.1080/19490976.2021.1993584",
            year=2021,
            result=(
                "233 people including 40 controls. Stool adhesion and "
                "biofilm-associated factors tracked decompensation and outcomes; "
                "changes were also studied around faecal transplant."
            ),
            supports=(
                "That disease-specific stool functional patterns exist beyond "
                "taxonomy."
            ),
            cannot_support=(
                "The same risk weights in otherwise healthy donors, or that lowering "
                "a score changes an outcome."
            ),
            evidence_class="human_association",
            specimen="stool",
            denominators={"participants": 233, "controls": 40},
        ),
        Source(
            id="BF-S22",
            citation="Nine-species gut-community RNA and metabolite study, Microorganisms",
            doi="10.3390/microorganisms13020234",
            year=2025,
            result=(
                "Mixed-species versus single-species laboratory biofilms differed in "
                "expression and metabolites. BOTH comparison arms were biofilms."
            ),
            supports="An RNA research dataset and an interaction hypothesis resource.",
            cannot_support=(
                "An active-biofilm-versus-planktonic training set. It cannot be "
                "relabelled as one."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="laboratory_community",
        ),
        Source(
            id="BF-S23",
            citation="Opposing Lactobacillus and Klebsiella effects, Beneficial Microbes",
            doi="10.3920/bm2017.0002",
            year=2017,
            result=(
                "L. plantarum CIRM653 reduced K. pneumoniae biofilm in vitro, but "
                "treated mice MAINTAINED intestinal carriage while controls declined."
            ),
            supports=(
                "A mandatory counterexample: in-vitro disruption did not improve "
                "in-vivo clearance and went the wrong way."
            ),
            cannot_support=(
                "Ranking this strain as a proven helpful intestinal treatment on the "
                "strength of the biomass assay."
            ),
            evidence_class="animal_mechanism",
            specimen="isolate",
            organisms=("Lactiplantibacillus plantarum CIRM653", "Klebsiella pneumoniae"),
        ),
        Source(
            id="BF-S24",
            citation="Pomegranate extract and barrier repair, Nutrients",
            doi="10.3390/nu15071771",
            year=2023,
            result=(
                "Improved experimental-colitis recovery and barrier-related outcomes. "
                "Broad Enterobacterales biofilm inhibition was NOT demonstrated, with "
                "only a slight effect on C. freundii."
            ),
            supports="A barrier-support intervention candidate.",
            cannot_support=(
                "Strong Klebsiella or E. coli biofilm clearance. The study does not "
                "show it."
            ),
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
        ),
        Source(
            id="BF-S25",
            citation="L. johnsonii against attaching and effacing pathogens, Frontiers in Immunology",
            doi="10.3389/fimmu.2026.1749001",
            year=2026,
            result=(
                "EPEC biofilm effects in vitro plus improved pathogen and inflammation "
                "outcomes in a C. rodentium mouse model. The deposited PV739486 is a "
                "16S sequence, not a complete genome."
            ),
            supports="A strain-qualified protective-antagonism candidate.",
            cannot_support=(
                "Strain-resolved protective-gene detection from a 16S accession, or "
                "activity of putatively annotated metabolites."
            ),
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
            organisms=("Lactobacillus johnsonii",),
        ),
        Source(
            id="BF-S26",
            citation="Resveratrol and L. paracasei ATCC334, International Journal of Molecular Sciences",
            doi="10.3390/ijms21155423",
            year=2020,
            result=(
                "Enhanced adhesion and biofilm formation in vitro, with clear strain "
                "dependence."
            ),
            supports="A protective-support research card for that strain.",
            cannot_support=(
                "That every Lactobacillus, or every resveratrol exposure, benefits "
                "the gut."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="isolate",
            organisms=("Lacticaseibacillus paracasei ATCC334",),
        ),
        Source(
            id="BF-S27",
            citation="L. plantarum Y42 biofilm state, Journal of Agricultural and Food Chemistry",
            doi="10.1021/acs.jafc.4c00460",
            year=2024,
            result=(
                "Biofilm-state organisms and their EPS and surface proteins supported "
                "barrier and immune outcomes in a mouse challenge model."
            ),
            supports="Broader protective mechanisms with strain, formulation and model intact.",
            cannot_support="EPS genes alone as a unique protective signature.",
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
            organisms=("Lactiplantibacillus plantarum Y42",),
        ),
        Source(
            id="BF-S28",
            citation="Lactiplantibacillus LR-1 biofilm state, Food & Function",
            doi="10.1039/d3fo02733c",
            year=2023,
            result=(
                "Biofilm-state cells outperformed planktonic cells on selected "
                "outcomes in experimental colitis."
            ),
            supports="That functional state is worth testing, not just species identity.",
            cannot_support="Universal probiotic efficacy, or a DNA-based activity claim.",
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
            organisms=("Lactiplantibacillus plantarum LR-1",),
        ),
        Source(
            id="BF-S29",
            citation="Multi-strain restorative biofilm consortium, World Journal of Microbiology and Biotechnology",
            doi="10.1007/s11274-026-05071-0",
            year=2026,
            result=(
                "Mouse studies of biofilm formation and recovery after antibiotics. "
                "Abstract-level evidence in this audit; full methods are needed before "
                "importing quantitative effects."
            ),
            supports="That consortia may restore functions, as a hypothesis.",
            cannot_support=(
                "That taxonomic co-presence in a sample recreates the administered "
                "formulation or its effect."
            ),
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
        ),
        Source(
            id="BF-S30",
            citation="Colon-targeted E. coli biofilm inhibitor, ACS Medicinal Chemistry Letters",
            doi="10.1021/acsmedchemlett.5c00675",
            year=2026,
            result=(
                "An experimental prodrug improved colonic delivery and selected "
                "mouse-colitis outcomes."
            ),
            supports="That colonic delivery can be engineered; an investigational approach.",
            cannot_support="An available supplement or a proven human treatment.",
            evidence_class="animal_mechanism",
            specimen="animal_tissue",
            organisms=("Escherichia coli",),
        ),
        Source(
            id="BF-S31",
            citation="B. subtilis protective matrix and delivery, Gut Microbes",
            doi="10.1080/19490976.2026.2684066",
            year=2026,
            result=(
                "Experimental systems support matrix-mediated stress protection and "
                "delivery."
            ),
            supports="Context for beneficial biofilm formulations.",
            cannot_support=(
                "That ordinary species carriage or stool gene abundance establishes "
                "the studied formulation or any clinical benefit."
            ),
            evidence_class="in_vitro_mechanism",
            specimen="laboratory_community",
            organisms=("Bacillus subtilis",),
        ),
        Source(
            id="BF-S32",
            citation="Indian CRC tissue biofilms, Applied Microbiology and Biotechnology",
            doi="10.1007/s00253-025-13537-8",
            year=2025,
            result=(
                "Imaging of 15 tumours and 15 paired adjacent samples found "
                "heterogeneous species distributions."
            ),
            supports="Geographic and tissue validation context.",
            cannot_support=(
                "A validated faecal signature. Adjacent tissue is not an independent "
                "healthy person."
            ),
            evidence_class="human_association",
            specimen="biopsy",
            denominators={"tumours": 15, "paired_adjacent": 15},
        ),
    )
}

#: Sources whose denominators are most often conflated in secondary citation.
#: The self-test asserts each pair stays numerically distinct.
GUARDED_DENOMINATORS: Final = {
    "BF-S01": ("screened", "included"),
    "BF-S20": ("sequenced_specimens", "uc_patients"),
}


def get(source_id: str) -> Source:
    """Look up one source, with a useful error when the ID is a typo."""
    try:
        return SOURCES[source_id]
    except KeyError:
        raise KeyError(
            f"unknown source {source_id!r}; known IDs are BF-S01..BF-S32"
        ) from None


def resolve_all(source_ids: tuple[str, ...]) -> tuple[Source, ...]:
    """Resolve a tuple of IDs, raising on the first unknown one."""
    return tuple(get(sid) for sid in source_ids)


def best_evidence_class(source_ids: tuple[str, ...]) -> str:
    """The strongest transport class among these sources.

    Used to cap a module's claim. Ordering follows :data:`EVIDENCE_CLASSES`,
    so ``human_association`` wins over ``animal_mechanism``. ``review`` and
    ``resource`` never win: a review of primary work is not itself primary
    evidence, so if a module cites only reviews it is capped at review.
    """
    if not source_ids:
        return "in_silico_hypothesis"
    order = {name: i for i, name in enumerate(EVIDENCE_CLASSES)}
    return min((get(sid).evidence_class for sid in source_ids), key=lambda c: order[c])
