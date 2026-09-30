"""Additional metabolic capacities — spec 0.8.3 section 5.8.

Urate degradation, the three polyamines, glutathione, and microbial drug
metabolism. The specification calls these "research additions" and is
explicit that research is a visible evidence level rather than permission to
use untraceable markers, so each capacity declares the same things a
required panel declares: the reaction it establishes, the claim it refuses,
the genes that carry it, and the source.

Two of these are easy to overstate and are written defensively.

Polyamines have no good direction. They are required for the gut lining to
renew itself, and they are elevated in colorectal tumour tissue. A report
that called a high reading "good" or "bad" would be inventing a direction
the literature does not have, so all three are context-dependent and say so.

Drug metabolism is not a measurement here at all. Whether a drug is
metabolised by your gut bacteria depends on the enzyme, the strain, the dose
and the timing, and none of that is visible in a stool gene count. These
appear as mechanism cards attached to genes that are measured elsewhere, and
each one states plainly that no dose or medication decision follows from it.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.extension.biotransform import Requirement

METHOD_CAPACITIES: Final = "ext083.additional_capacities/1.0"

#: Every capacity here is a research addition, and says so in the report.
EVIDENCE_RESEARCH: Final = "research"

STATES: Final[tuple[str, ...]] = (
    "complete", "partial", "absent", "not_assayed",
)


class CapacityError(ValueError):
    """A capacity that would claim more than its genes support."""


@dataclass(frozen=True)
class Capacity:
    """One additional metabolic capacity."""

    capacity_id: str
    label: str
    requirement: Requirement
    #: The panel whose reading this capacity is reported against.
    panel: str
    establishes: str
    does_not_establish: str
    #: Why the direction of this reading is not simply good or bad. Required
    #: for anything context-dependent, which is most of this module.
    context: str
    sources: tuple[str, ...] = ()
    #: Genes measured and shown, but deliberately not counted towards the
    #: route, with the reason each is excluded.
    reported_separately: Mapping[str, str] = field(default_factory=dict)
    note: str = ""

    def __post_init__(self) -> None:
        if not self.does_not_establish:
            raise CapacityError(
                f"{self.capacity_id}: every capacity must say what it does not establish"
            )
        if not self.context:
            raise CapacityError(
                f"{self.capacity_id}: every capacity must say why its direction is not "
                "simply favourable or adverse"
            )
        overlap = set(self.requirement.genes) & set(self.reported_separately)
        if overlap:
            raise CapacityError(
                f"{self.capacity_id}: {sorted(overlap)} is both required and excluded"
            )

    def assess(self, detected: frozenset[str], searched: frozenset[str]) -> str:
        genes = set(self.requirement.genes)
        if not genes & searched:
            return "not_assayed"
        if self.requirement.satisfied_by(detected):
            return "complete"
        return "partial" if genes & detected else "absent"

    def to_json(
        self, detected: frozenset[str] | None = None,
        searched: frozenset[str] | None = None,
    ) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "capacity_id": self.capacity_id,
            "label": self.label,
            "panel": self.panel,
            "evidence_level": EVIDENCE_RESEARCH,
            "requirement": self.requirement.to_json(detected),
            "establishes": self.establishes,
            "does_not_establish": self.does_not_establish,
            "context": self.context,
            "reported_separately": dict(self.reported_separately),
            "sources": list(self.sources),
            "note": self.note,
        }
        if detected is not None and searched is not None:
            payload["state"] = self.assess(detected, searched)
        return payload


CAPACITIES: Final[tuple[Capacity, ...]] = (
    Capacity(
        capacity_id="capacity.urate_anaerobic",
        label="Anaerobic urate degradation",
        requirement=Requirement(all_of=("xdhA",), any_of=(("xdhB", "xdhC"),)),
        panel="urate",
        establishes=(
            "the xanthine dehydrogenase route by which gut bacteria degrade urate "
            "without oxygen"
        ),
        does_not_establish=(
            "your blood urate level, and nothing about gout. Urate in blood is set "
            "mostly by how much your body makes and how much your kidneys excrete; a "
            "gene count in stool cannot substitute for a blood test"
        ),
        context=(
            "Degrading urate is neither good nor bad in itself. The mouse work links "
            "this route to lower host urate, but the same enzymes participate in purine "
            "fermentation generally, and the human evidence is association."
        ),
        sources=("F28", "F29"),
        reported_separately={
            "uox": (
                "uricase needs molecular oxygen, which is scarce in the colon, so it is "
                "a different reaction and is not counted towards the anaerobic route"
            ),
            "hpt": (
                "purine salvage recycles purines back into nucleotides instead of "
                "degrading them, which is the opposite of what this capacity measures"
            ),
        },
        note=(
            "The terminal steps of purine fermentation run through glycine, whose "
            "reductase is measured by the Stickland fermentation panel rather than here."
        ),
    ),
    Capacity(
        capacity_id="capacity.putrescine",
        label="Putrescine formation",
        requirement=Requirement(any_of=(("speC", "speB", "aguB"),)),
        panel="putrescine",
        establishes="at least one of the three bacterial routes that release putrescine",
        does_not_establish=(
            "a putrescine concentration, and nothing about polyamines measured in blood "
            "or urine, which come overwhelmingly from your own cells and your diet"
        ),
        context=(
            "Polyamines are required for the gut lining to renew itself, and they are "
            "also elevated in colorectal tumour tissue. The same molecule is doing both, "
            "so more is not better and less is not safer."
        ),
        sources=("F30", "F31"),
        reported_separately={
            "aguA": (
                "the deiminase is the step before aguB in the same two-step route, so "
                "counting both would count one organism's single capacity twice"
            ),
        },
    ),
    Capacity(
        capacity_id="capacity.agmatine",
        label="Agmatine formation",
        requirement=Requirement(all_of=("adiA",)),
        panel="agmatine",
        establishes="arginine decarboxylation to agmatine",
        does_not_establish=(
            "any agmatine reaching your circulation or your brain, and nothing about "
            "the doses used in the animal work on pain and mood, which are "
            "pharmacological"
        ),
        context=(
            "The same enzyme is part of how several gut bacteria survive stomach acid, "
            "so a high reading can be about acid resistance rather than about polyamine "
            "or neuroactive chemistry."
        ),
        sources=("F30", "F31"),
    ),
    Capacity(
        capacity_id="capacity.spermidine",
        label="Spermidine formation",
        requirement=Requirement(any_of=(("speE", "casDC"),)),
        panel="spermidine",
        establishes=(
            "spermidine synthesis by either the aminopropyltransferase route or the "
            "carboxyspermidine route that most gut Bacteroidetes actually use"
        ),
        does_not_establish=(
            "a spermidine level, and nothing about the autophagy and longevity work, "
            "which is largely from model organisms and from dietary intake rather than "
            "from bacterial synthesis"
        ),
        context=(
            "Counting only spermidine synthase would report an absence in most "
            "Bacteroidetes, which make spermidine by a different route entirely. Both "
            "routes are counted here, and neither direction of the result is a health "
            "claim."
        ),
        sources=("F30", "F31"),
        reported_separately={
            "speD": (
                "this supplies the aminopropyl group rather than forming spermidine, so "
                "it is route detail and not a second spermidine-forming reaction"
            ),
            "casDH": (
                "the dehydrogenase is the step before the decarboxylase in the same "
                "route, and counting both would double one capacity"
            ),
        },
    ),
    Capacity(
        capacity_id="capacity.glutathione",
        label="Glutathione synthesis",
        requirement=Requirement(any_of=(("gshF", "gshA"), ("gshF", "gshB"))),
        panel="glutathione",
        establishes=(
            "bacterial glutathione synthesis, by the two-enzyme route or by the fused "
            "single-enzyme version"
        ),
        does_not_establish=(
            "your own glutathione status. Host glutathione is made by your own cells "
            "from cysteine, and there is no established route by which bacterial "
            "glutathione in the colon becomes glutathione in your tissues"
        ),
        context=(
            "Many gut anaerobes import glutathione instead of making it, so a low "
            "synthesis reading is not an absence of glutathione in the community. "
            "Recycling and breakdown are measured separately for the same reason."
        ),
        sources=("F32",),
        reported_separately={
            "gor": "the reductase recycles glutathione that already exists rather than making it",
            "ggt": "the transpeptidase takes glutathione apart, which is the opposite reaction",
        },
        note=(
            "Requiring gshA and gshB together, or gshF alone, is what makes this a route "
            "rather than a gene count: the fused enzyme does both steps, and a search "
            "for the two-enzyme route alone would miss every organism that carries it."
        ),
    ),
    Capacity(
        capacity_id="capacity.ethanol",
        label="Microbial ethanol formation",
        requirement=Requirement(all_of=("adhE",)),
        panel="ethanol",
        establishes=(
            "the bifunctional route from acetyl-CoA through acetaldehyde to ethanol"
        ),
        does_not_establish=(
            "auto-brewery syndrome, fatty liver disease, or any blood alcohol level. The "
            "liver findings concern particular high-producing strains, and carrying the "
            "gene is ordinary"
        ),
        context=(
            "Low-level endogenous ethanol production is normal. What the fatty liver "
            "work identified was strain-level over-production, which a gene count cannot "
            "distinguish from the ordinary case."
        ),
        sources=("C03", "C04", "C05"),
        reported_separately={
            "eutB": "the ethanolamine route produces acetaldehyde from a different substrate entirely",
            "eutC": "the second subunit of that same separate route",
            "pdc": "decarboxylating pyruvate is a different first step, common in yeast and rare in gut bacteria",
            "adhP": "a separate alcohol dehydrogenase that prefers propanol, reported for completeness",
        },
        note=(
            "No acetaldehyde disposal reading is given. The aldehyde dehydrogenase family "
            "is not substrate-specific, so a count over it would measure aldehyde "
            "chemistry in general rather than acetaldehyde, and would take references "
            "from the propanediol and GABA routes that use the same family."
        ),
    ),
)

BY_ID: Final[Mapping[str, Capacity]] = {c.capacity_id: c for c in CAPACITIES}


# --------------------------------------------------------------------------- #
# drug metabolism — cards, not measurements
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DrugReaction:
    """One experimentally characterised microbial drug transformation."""

    reaction_id: str
    drug: str
    enzyme: str
    organism: str
    mechanism: str
    #: The inference a reader will reach for, refused explicitly.
    does_not_mean: str
    #: The gene symbol, when this report measures it at all.
    gene: str | None
    measured_by: str | None
    source: str

    def to_json(self, detected: frozenset[str] | None = None) -> dict[str, Any]:
        state: str | None = None
        if detected is not None and self.gene:
            state = "gene_detected" if self.gene in detected else "gene_not_detected"
        return {
            "reaction_id": self.reaction_id, "drug": self.drug, "enzyme": self.enzyme,
            "organism": self.organism, "mechanism": self.mechanism,
            "does_not_mean": self.does_not_mean, "gene": self.gene,
            "measured_by": self.measured_by, "source": self.source,
            "evidence_level": EVIDENCE_RESEARCH, "state": state,
        }


@dataclass(frozen=True)
class NeuroactiveCard:
    """One neurotransmitter, and what microbial evidence exists for it.

    Spec §5.5 asks for five of these and gives each an endpoint boundary,
    because this is the part of the field where the gap between what can be
    measured and what people want to conclude is widest. Four of the five
    have no validated sequence panel, and the specification says so
    directly: a generic acetyltransferase is not an acetylcholine assay, and
    a tryptamine-producing enzyme must not be upgraded to a serotonin
    detector merely because a homolog is retrievable.

    So each card records what it links to, what the experiments actually
    showed, and where the evidence stops. ``measured_by`` is the existing
    reading it points at; ``sequence_panel`` is whether a specific detector
    exists at all, and for most of these the honest answer is no.
    """

    card_id: str
    label: str
    #: Existing readings in this report that bear on it.
    links_to: tuple[str, ...]
    #: What the experimental work established, at strain level where that is
    #: the level it was done at.
    evidence: str
    #: The boundary the specification sets for this card, verbatim in intent.
    endpoint_boundary: str
    #: Whether a validated substrate-specific sequence panel exists here.
    sequence_panel: str | None
    #: Why there is no panel, when there is none.
    no_panel_because: str = ""
    #: An imported laboratory result that would bear on it, if supplied.
    imported_measurement: str = ""
    sources: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.endpoint_boundary:
            raise CapacityError(f"{self.card_id}: every card must state its boundary")
        if self.sequence_panel is None and not self.no_panel_because:
            raise CapacityError(
                f"{self.card_id}: a card with no sequence panel must say why, or a "
                "reader will assume the chemistry was searched for and absent"
            )

    def to_json(self, detected: frozenset[str] | None = None) -> dict[str, Any]:
        return {
            "card_id": self.card_id, "label": self.label,
            "links_to": list(self.links_to),
            "evidence": self.evidence,
            "endpoint_boundary": self.endpoint_boundary,
            "sequence_panel": self.sequence_panel,
            "no_panel_because": self.no_panel_because,
            "imported_measurement": self.imported_measurement,
            "sources": list(self.sources),
            "evidence_level": EVIDENCE_RESEARCH,
            "state": (
                None if detected is None or not self.sequence_panel
                else ("gene_detected" if self.sequence_panel in detected
                      else "gene_not_detected")
            ),
        }


NEUROACTIVE_CARDS: Final[tuple[NeuroactiveCard, ...]] = (
    NeuroactiveCard(
        card_id="neuro.dopamine",
        label="Dopamine-related microbial metabolism",
        links_to=("dopamine", "tyrDC", "dadh"),
        evidence=(
            "Enterococcus faecalis tyrosine decarboxylase decarboxylates levodopa to "
            "dopamine in the small intestine, and Eggerthella lenta dopamine "
            "dehydroxylase then removes a hydroxyl to give m-tyramine. Both enzymes "
            "were purified and both reactions were demonstrated; both genes are "
            "measured in this report."
        ),
        endpoint_boundary=(
            "Not brain dopamine, and not a medication instruction. These reactions "
            "happen in the gut lumen on a drug or on dietary tyrosine; dopamine made "
            "there does not cross into the brain, and nothing here indicates whether "
            "your levodopa, if you take any, is being consumed."
        ),
        sequence_panel="tyrDC",
        sources=("F15", "F16"),
    ),
    NeuroactiveCard(
        card_id="neuro.serotonin",
        label="Serotonin-related microbial metabolism",
        links_to=("tryptamine", "indole", "ipa"),
        evidence=(
            "A 2025 result describes a Limosilactobacillus mucosae and Ligilactobacillus "
            "ruminis consortium whose effect depends on 5-hydroxytryptophan being "
            "available, with a traceable candidate decarboxylase. The chemistry is "
            "consortium-level and 5-HTP-dependent, and the construct's substrate "
            "specificity is recorded as unresolved in the source itself."
        ),
        endpoint_boundary=(
            "This is not a serotonin measurement and not an inference about serotonin "
            "in your blood or brain. Most of the body's serotonin is made by your own "
            "enterochromaffin cells, not by bacteria; the two pools do not exchange "
            "across the blood-brain barrier."
        ),
        sequence_panel=None,
        no_panel_because=(
            "There is no validated serotonin-specific gene panel. The available "
            "decarboxylase evidence is for tryptamine, which is a different molecule, "
            "and turning a tryptamine enzyme into a serotonin detector because a "
            "homolog can be retrieved is precisely the error this card exists to avoid. "
            "The tryptamine reading is linked instead, labelled as tryptamine."
        ),
        imported_measurement=(
            "A measured serotonin or 5-HIAA result can be entered in the lab section "
            "and will be shown beside this card rather than merged into it."
        ),
        sources=("F17", "F18"),
    ),
    NeuroactiveCard(
        card_id="neuro.acetylcholine",
        label="Acetylcholine",
        links_to=("cutc",),
        evidence=(
            "Some lactobacilli have been reported to produce acetylcholine in culture, "
            "and choline is the shared precursor of both acetylcholine and the "
            "trimethylamine route this report does measure. The enzymes responsible in "
            "bacteria are not characterised to the level a specific detector needs."
        ),
        endpoint_boundary=(
            "A generic acetyltransferase count is not an acetylcholine assay. Nothing "
            "here measures acetylcholine, and nothing here bears on cognition or "
            "neuromuscular function."
        ),
        sequence_panel=None,
        no_panel_because=(
            "No validated substrate-specific panel exists. Acetyltransferases are one "
            "of the largest and least substrate-specific enzyme families in bacteria, "
            "so a count over them would measure the family and be reported as "
            "acetylcholine, which is worse than reporting nothing."
        ),
        sources=("F15",),
    ),
    NeuroactiveCard(
        card_id="neuro.norepinephrine",
        label="Norepinephrine",
        links_to=("bglucuronidase",),
        evidence=(
            "The microbial contribution here is handling rather than synthesis: "
            "bacterial beta-glucuronidase can deconjugate catecholamine glucuronides in "
            "the gut lumen, freeing the parent compound. Host adrenergic signalling is "
            "host physiology and is a separate matter from anything in a stool sample."
        ),
        endpoint_boundary=(
            "No bacterial abundance converts to a norepinephrine level. Microbial "
            "deconjugation and your own adrenergic response are different things, and "
            "this card keeps them apart rather than joining them."
        ),
        sequence_panel=None,
        no_panel_because=(
            "Bacterial norepinephrine synthesis is not established, so there is nothing "
            "specific to search for. The deconjugation capacity is measured by the "
            "existing beta-glucuronidase reading, which is linked here and is not "
            "substrate-specific to catecholamines."
        ),
        imported_measurement=(
            "A measured plasma or urinary catecholamine result belongs with its own "
            "assay context and can be entered in the lab section."
        ),
        sources=("F39",),
    ),
    NeuroactiveCard(
        card_id="neuro.histamine",
        label="Histamine detail",
        links_to=("histamine", "hdcA"),
        evidence=(
            "Bacterial histidine decarboxylase converts histidine to histamine, and it "
            "is measured in this report whatever organism carries it. A specific "
            "Klebsiella aerogenes strain producing histamine was linked to visceral "
            "hypersensitivity in a mouse model and to a patient subgroup, which is a "
            "strain-level finding rather than a species-level one."
        ),
        endpoint_boundary=(
            "Potential production is not a diagnosis of histamine intolerance. That "
            "condition is defined clinically and is not established from stool DNA, and "
            "a high reading here is not a reason to start a low-histamine diet."
        ),
        sequence_panel="hdcA",
        no_panel_because=(
            "Production is measured; degradation is not. Bacterial histamine-degrading "
            "enzymes are not discriminable at the sequence level from the wider amine "
            "oxidase and dehydrogenase families, so a degradation route is not reported "
            "rather than being reported unreliably."
        ),
        sources=("F15", "C13"),
    ),
)


DRUG_REACTIONS: Final[tuple[DrugReaction, ...]] = (
    DrugReaction(
        reaction_id="drug.digoxin",
        drug="Digoxin",
        enzyme="cardiac glycoside reductase (Cgr2)",
        organism="Eggerthella lenta, and only some strains of it",
        mechanism=(
            "Cgr2 reduces the lactone ring of digoxin to dihydrodigoxin, which is "
            "inactive. Strains of E. lenta differ in whether they carry a functional "
            "cgr operon, and the reaction is inhibited by arginine rather than "
            "increased by it."
        ),
        does_not_mean=(
            "Nothing here estimates a digoxin level or supports a dose change. Carrying "
            "E. lenta is not the same as carrying a reducing strain, and this report "
            "does not measure the cgr operon: too few characterised sequences exist for "
            "a read-based search to be sensitive."
        ),
        gene=None,
        measured_by=None,
        source=(
            "Haiser H.J. et al., Science 341:295-298 (2013), doi:10.1126/science.1235872; "
            "Koppel N. et al., eLife 7:e33953 (2018), doi:10.7554/eLife.33953."
        ),
    ),
    DrugReaction(
        reaction_id="drug.levodopa",
        drug="Levodopa (L-DOPA)",
        enzyme="tyrosine decarboxylase (TyrDC), then dopamine dehydroxylase (Dadh)",
        organism="Enterococcus faecalis, then Eggerthella lenta",
        mechanism=(
            "E. faecalis TyrDC decarboxylates levodopa to dopamine in the small "
            "intestine, before the drug reaches the brain. E. lenta Dadh then "
            "dehydroxylates that dopamine to m-tyramine. The first step is not blocked "
            "by carbidopa, the inhibitor given alongside levodopa for the human enzyme."
        ),
        does_not_mean=(
            "Carrying these genes does not mean your levodopa is being consumed, and no "
            "dose change follows from this. Both genes are measured in this report "
            "because they have other roles; their presence is common."
        ),
        gene="tyrDC",
        measured_by="dopamine",
        source=(
            "Maini Rekdal V. et al., Science 364:eaau6323 (2019), "
            "doi:10.1126/science.aau6323."
        ),
    ),
    DrugReaction(
        reaction_id="drug.irinotecan",
        drug="Irinotecan (via its metabolite SN-38)",
        enzyme="beta-glucuronidase",
        organism="many gut bacteria, with large differences between enzymes",
        mechanism=(
            "The liver inactivates SN-38 by attaching glucuronic acid. Bacterial "
            "beta-glucuronidases remove it again in the gut, regenerating the active "
            "drug in the intestinal lining, which is a proposed mechanism of the "
            "diarrhoea that limits this chemotherapy."
        ),
        does_not_mean=(
            "Your total beta-glucuronidase reading does not estimate this. The enzymes "
            "that perform the reaction efficiently are one structural class among "
            "several, and this report cannot resolve which class yours belong to."
        ),
        gene="gus",
        measured_by="bglucuronidase",
        source=(
            "Wallace B.D. et al., Science 330:831-835 (2010), doi:10.1126/science.1191175."
        ),
    ),
    DrugReaction(
        reaction_id="drug.sulfasalazine",
        drug="Sulfasalazine",
        enzyme="azoreductase",
        organism="a broad range of gut bacteria",
        mechanism=(
            "Sulfasalazine is a prodrug: it passes the small intestine intact and gut "
            "bacterial azoreductases cleave its azo bond to release the active "
            "5-aminosalicylate in the colon, which is how the drug is designed to work."
        ),
        does_not_mean=(
            "This is the intended mechanism of the drug rather than a problem, and no "
            "reading here estimates how much is released. Azoreductase is not measured "
            "in this report: the protein name covers an enzyme family too broad for a "
            "count over it to mean this reaction."
        ),
        gene=None,
        measured_by=None,
        source=(
            "Zimmermann M. et al., Nature 570:462-467 (2019), "
            "doi:10.1038/s41586-019-1291-3."
        ),
    ),
)


# --------------------------------------------------------------------------- #
# assembly
# --------------------------------------------------------------------------- #


def summarise(
    detected: Sequence[str], *, searched: Sequence[str],
    panels: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Every additional capacity, with its state in this sample."""
    have = frozenset(detected)
    looked = frozenset(searched)
    panels = panels or {}
    rows: list[dict[str, Any]] = []
    for capacity in CAPACITIES:
        payload = capacity.to_json(have, looked)
        panel = panels.get(capacity.panel)
        if isinstance(panel, Mapping):
            payload["copies_per_100_genomes"] = panel.get("copies_per_100_genomes")
            payload["fragments"] = panel.get("accepted_fragments")
        rows.append(payload)
    counted = [r for r in rows if r["state"] != "not_assayed"]
    return {
        "method": METHOD_CAPACITIES,
        "evidence_level": EVIDENCE_RESEARCH,
        "n_capacities": len(rows),
        "n_assayed": len(counted),
        "n_complete": sum(1 for r in rows if r["state"] == "complete"),
        "capacities": rows,
        "drug_reactions": [d.to_json(have) for d in DRUG_REACTIONS],
        # §5.5's five cards. Four of them have no validated sequence panel
        # and say so, which is the point: the gap between what can be
        # measured and what a reader wants to conclude is widest here.
        "neuroactive_cards": [c.to_json(have) for c in NEUROACTIVE_CARDS],
        "n_neuroactive_cards": len(NEUROACTIVE_CARDS),
        "standing_note": (
            "These are research additions. Each is a count of genes in stool, which is a "
            "statement about what the community could do and never about what it is "
            "doing, or about a concentration in you. None of them supports changing a "
            "medication."
        ),
    }
