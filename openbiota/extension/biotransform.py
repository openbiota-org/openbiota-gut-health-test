"""Neuroactive, plant-compound and surface biotransformations.

A05 (§5.5), A06 (§5.6) and A08 (§5.7) share one shape: a named chemical
step, a specific set of characterised enzymes, and a long list of things
the enzymes do *not* establish. They live together because the refusals
are the same refusals.

The refusals, each enforced rather than described:

**A difference between two percentiles is not a net flux.** GABA
synthesis and GABA degradation are shown side by side, never subtracted.
Where a ratio is shown at all it uses compatible raw units, states its
formula, and prints both components.

**A domain is not a reaction.** Equol conversion needs the validated
DZNR, DHDR and THDR chemistry; a general reductase domain is not
evidence. Urolithin needs all three ucdCFO components, and even then
supports one step of the pathway rather than the whole of it.

**A negative control is part of the method.** UrdA acts on urocanate and
makes imidazole propionate. It is close enough to the urolithin
dehydroxylases to be mistaken for one, so it is declared as an explicit
negative control and a panel that matches it is failing, not finding.

**Capacity is not damage.** An organism that uses mucin glycans is doing
something physiological. Genetic capacity is not an image of an eroded
mucus layer, and total Gram-negative abundance is not a TLR4 stimulus.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, Final

METHOD_BIOTRANSFORM: Final = "ext083.biotransformation/1.0"

FEATURES: Final[Mapping[str, str]] = {
    "A05": "Neuroactive metabolism",
    "A06": "Plant-compound conversion",
    "A07": "GLP-1-related microbial mechanisms",
    "A08": "Mucin, lipid A and bile acids",
}


class BiotransformError(ValueError):
    """A step that would claim a reaction its genes do not establish."""


@dataclass(frozen=True)
class Requirement:
    """A boolean requirement over gene symbols.

    Written as a small expression rather than a flat list because the
    specification's requirements genuinely are boolean: the glucosinolate
    core is `BT2158 AND (BT2156 OR BT2157)`, and flattening that to five
    genes would either over- or under-call it.
    """

    all_of: tuple[str, ...] = ()
    any_of: tuple[tuple[str, ...], ...] = ()

    def __post_init__(self) -> None:
        if not self.all_of and not self.any_of:
            raise BiotransformError("a requirement must require something")

    @property
    def genes(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(
            [*self.all_of, *[g for group in self.any_of for g in group]]
        ))

    def satisfied_by(self, present: frozenset[str]) -> bool:
        if not all(g in present for g in self.all_of):
            return False
        return all(any(g in present for g in group) for group in self.any_of)

    def describe(self) -> str:
        parts = list(self.all_of)
        parts += [f"({' or '.join(group)})" for group in self.any_of]
        return " and ".join(parts)

    def to_json(self, present: frozenset[str] | None = None) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "expression": self.describe(), "genes": list(self.genes),
            "all_of": list(self.all_of),
            "any_of": [list(group) for group in self.any_of],
        }
        if present is not None:
            payload["satisfied"] = self.satisfied_by(present)
            payload["present"] = [g for g in self.genes if g in present]
            payload["missing"] = [g for g in self.genes if g not in present]
        return payload


@dataclass(frozen=True)
class Step:
    """One characterised chemical step."""

    step_id: str
    feature_id: str
    label: str
    requirement: Requirement
    #: The step this establishes, stated narrowly.
    establishes: str
    #: What it does *not* establish, stated at least as narrowly.
    does_not_establish: str
    #: Enzymes close enough to be mistaken for these, whose match is a
    #: failure of the panel rather than a finding.
    negative_controls: Mapping[str, str] = field(default_factory=dict)
    food_context: str = ""
    links_to: tuple[str, ...] = ()
    source_ids: tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.feature_id not in FEATURES:
            raise BiotransformError(f"{self.step_id}: unknown feature {self.feature_id!r}")
        if not self.does_not_establish:
            raise BiotransformError(
                f"{self.step_id}: every step must say what it does not establish"
            )
        overlap = set(self.requirement.genes) & set(self.negative_controls)
        if overlap:
            raise BiotransformError(
                f"{self.step_id}: {sorted(overlap)} is both required and a negative control"
            )

    def to_json(self, present: frozenset[str] | None = None) -> dict[str, Any]:
        return {
            "step_id": self.step_id, "feature_id": self.feature_id,
            "feature": FEATURES[self.feature_id], "label": self.label,
            "requirement": self.requirement.to_json(present),
            "establishes": self.establishes,
            "does_not_establish": self.does_not_establish,
            "negative_controls": dict(self.negative_controls),
            "food_context": self.food_context,
            "links_to": list(self.links_to),
            "source_ids": list(self.source_ids),
            "note": self.note,
        }


# --------------------------------------------------------------------------- #
# A05 — neuroactive
# --------------------------------------------------------------------------- #

NEUROACTIVE: Final[tuple[Step, ...]] = (
    Step(
        step_id="neuro.gaba_synthesis", feature_id="A05",
        label="GABA synthesis",
        requirement=Requirement(any_of=(("gadA", "gadB"),)),
        establishes="glutamate decarboxylation to GABA",
        does_not_establish=(
            "a GABA concentration in your gut, and nothing about GABA reaching your "
            "brain: the blood-brain barrier is not permeable to it"
        ),
        links_to=("gaba",),
        source_ids=("F15", "F16"),
        note="MGB020-022 in the Gut Brain Module definitions.",
    ),
    Step(
        step_id="neuro.gaba_degradation", feature_id="A05",
        label="GABA degradation",
        requirement=Requirement(all_of=("gabT", "gabD")),
        establishes="the GabT/GabD route from GABA to succinate semialdehyde and on",
        does_not_establish=(
            "a net GABA balance. Degradation and synthesis are separate measurements "
            "and subtracting one percentile from another is not a flux"
        ),
        links_to=("gaba",),
        source_ids=("F15", "F16"),
        note="MGB019. The required measured addition of section 5.5.",
    ),
    Step(
        step_id="neuro.gaba_to_butyrate", feature_id="A05",
        label="4-aminobutyrate to butyrate",
        requirement=Requirement(all_of=("gabT", "abfH", "abfD")),
        establishes="the supported route from 4-aminobutyrate through to butyrate",
        does_not_establish=(
            "that GABA is the main butyrate source here; the fibre routes are usually "
            "much larger and are measured separately"
        ),
        links_to=("gaba", "butyrate"),
        source_ids=("F15",),
    ),
)


# --------------------------------------------------------------------------- #
# A06 — plant compounds
# --------------------------------------------------------------------------- #

PLANT_COMPOUNDS: Final[tuple[Step, ...]] = (
    Step(
        step_id="polyphenol.urolithin_9_dehydroxylation", feature_id="A06",
        label="Urolithin C to urolithin A (9-dehydroxylation)",
        requirement=Requirement(all_of=("ucdC", "ucdF", "ucdO")),
        establishes=(
            "the characterised ucdCFO 9-dehydroxylation, one step of the ellagitannin "
            "pathway"
        ),
        does_not_establish=(
            "the whole ellagitannin-to-urolithin pathway, and not a urolithin producer "
            "metabotype: that is a measured phenotype from a substrate challenge, not a "
            "DNA finding"
        ),
        negative_controls={
            "urdA": (
                "UrdA acts on urocanate and makes imidazole propionate. It is close "
                "enough to be mistaken for a urolithin dehydroxylase, so a panel that "
                "matches it is failing rather than finding. Residue discrimination must "
                "use aligned reference numbering over covered diagnostic positions."
            ),
        },
        food_context="Pomegranate, walnuts, berries and other ellagitannin foods",
        source_ids=("F19", "F20", "F22"),
        note=(
            "Not measured in this release. The operon is described in GenBank "
            "PQ855390.1, but UniProt carries no bacterial entry under ucdC, ucdF or "
            "ucdO - those symbols there belong to an unrelated Aspergillus ustus gene "
            "cluster. Measuring against a mould would manufacture exactly the false "
            "positive this step's negative control exists to prevent, so the reading "
            "waits for the operon's sequences rather than approximating them."
        ),
    ),
    Step(
        step_id="polyphenol.equol_daidzein_conversion", feature_id="A06",
        label="Daidzein to equol",
        requirement=Requirement(all_of=("dznr", "ddr", "tdr")),
        establishes="the validated DZNR, DHDR and THDR chemistry of equol formation",
        does_not_establish=(
            "an equol producer metabotype, which is a measured phenotype; and a general "
            "reductase domain is not evidence for any of these three reactions"
        ),
        negative_controls={
            "generic_reductase": (
                "An unassigned oxidoreductase domain matches many things. The aliases "
                "eqlA, eqlB and eqlC need source- and strain-specific mapping before "
                "they can be read as these reactions."
            ),
        },
        food_context="Soy and other isoflavone foods",
        source_ids=("F21",),
    ),
    Step(
        step_id="diet.glucosinolate_isothiocyanate_conversion", feature_id="A06",
        label="Glucosinolate to isothiocyanate",
        # The experimentally supported in-vitro core, exactly as specified.
        requirement=Requirement(all_of=("BT2158",), any_of=(("BT2156", "BT2157"),)),
        establishes=(
            "the bacterial route characterised in B. thetaiotaomicron, whose in-vitro "
            "core is BT2158 with either BT2156 or BT2157"
        ),
        does_not_establish=(
            "plant myrosinase activity, which is a different enzyme from the vegetable "
            "itself, nor the other bacterial routes that have been described"
        ),
        food_context="Broccoli, rocket and other cruciferous vegetables",
        source_ids=("F23",),
        note="Full BT2159-BT2156 locus context is retained where assembly permits.",
    ),
)


# --------------------------------------------------------------------------- #
# A08 — mucin, lipid A, bile acids
# --------------------------------------------------------------------------- #

SURFACE: Final[tuple[Step, ...]] = (
    Step(
        step_id="mucin.glycan_foraging", feature_id="A08",
        label="Mucin glycan foraging",
        requirement=Requirement(all_of=("nanH",), any_of=(("fucA", "afcA"), ("sulf", "sgl"))),
        establishes=(
            "sialidase with fucosidase and sulfatase activity on mucin glycans, which is "
            "how a mucin forager makes a living"
        ),
        does_not_establish=(
            "damage to your mucus barrier. Mucin foraging is physiological and these "
            "organisms are normal residents; genetic capacity is not an image of "
            "mucosal erosion"
        ),
        links_to=("biofilm",),
        source_ids=("F41", "F42"),
        note=(
            "Gastric and colonic mucin carry different glycans; substrate context is "
            "preserved rather than pooled, and these enzymes are kept apart from the "
            "dietary-fibre CAZymes."
        ),
    ),
    Step(
        step_id="lps.lipid_a_modification", feature_id="A08",
        label="Lipid A acylation and modification",
        requirement=Requirement(all_of=("lpxA", "lpxC"), any_of=(("lpxL", "lpxM", "lpxE"),)),
        establishes="the presence of lipid A biosynthesis with at least one modification enzyme",
        does_not_establish=(
            "a TLR4 stimulus or an endotoxin concentration. The number of acyl chains "
            "decides the immune effect, one modification gene does not resolve the "
            "structure, and total Gram-negative abundance resolves nothing at all"
        ),
        links_to=("biofilm",),
        source_ids=("F41",),
        note=(
            "Where the structure cannot be resolved the gene evidence is shown with its "
            "uncertainty, rather than a fabricated endotoxin figure."
        ),
    ),
    Step(
        step_id="bile.deconjugation", feature_id="A08",
        label="Bile-salt deconjugation",
        requirement=Requirement(all_of=("bsh",)),
        establishes="bile-salt hydrolase presence",
        does_not_establish=(
            "which bile salts are preferred, nor the direction: BSH has both hydrolase "
            "and acyltransferase activity, so the same enzyme can also reconjugate"
        ),
        links_to=("bsh",),
        source_ids=("F33", "F36"),
    ),
    Step(
        step_id="bile.dehydroxylation", feature_id="A08",
        label="7-alpha-dehydroxylation (bai pathway)",
        requirement=Requirement(all_of=("baiB", "baiCD", "baiE", "baiA")),
        establishes="a complete bai operon for primary-to-secondary bile-acid conversion",
        does_not_establish=(
            "a secondary bile-acid concentration. The pathway needs all its steps, which "
            "is why completeness is reported rather than a count of bai genes"
        ),
        links_to=("bai",),
        source_ids=("F33", "F34"),
    ),
    Step(
        step_id="bile.hsdh_transformation", feature_id="A08",
        label="Hydroxysteroid dehydrogenase transformations",
        requirement=Requirement(any_of=(("hdhA", "hsdh3", "hsdh7", "hsdh12"),)),
        establishes="position- and stereo-specific bile-acid oxidation or epimerisation",
        does_not_establish=(
            "one linear sequence of bile-acid chemistry. Some HSDHs act on conjugated "
            "substrates, so these reactions form a graph rather than a chain"
        ),
        links_to=("bai", "bsh"),
        source_ids=("F35", "F36"),
        note=(
            "BileActome's HMM families are candidate annotations validated against "
            "human-gut biochemical examples; its rumen validation population supplies no "
            "human percentile, and its E-value rule is a starting point rather than a "
            "clinical threshold."
        ),
    ),
)

ALL_STEPS: Final[tuple[Step, ...]] = (*NEUROACTIVE, *PLANT_COMPOUNDS, *SURFACE)
BY_ID: Final[Mapping[str, Step]] = {s.step_id: s for s in ALL_STEPS}


# --------------------------------------------------------------------------- #
# grading
# --------------------------------------------------------------------------- #

STATES: Final[tuple[str, ...]] = (
    "supported", "partial", "absent", "not_assayed", "negative_control_matched",
)


@dataclass(frozen=True)
class StepFinding:
    step: Step
    state: str
    present: tuple[str, ...]
    missing: tuple[str, ...]
    controls_matched: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        payload = {
            "method_id": METHOD_BIOTRANSFORM,
            **self.step.to_json(frozenset(self.present)),
            "state": self.state,
            "limitations": [
                "Genetic capacity for a named reaction, not a measured amount.",
                self.step.does_not_establish,
            ],
        }
        if self.controls_matched:
            payload["negative_controls_matched"] = {
                gene: self.step.negative_controls[gene] for gene in self.controls_matched
            }
            payload["limitations"].append(
                "A declared negative control matched. That is a signal about the panel, "
                "not a finding about you, and this step is withheld."
            )
        return payload


def assess_step(
    step: Step, detected: Iterable[str], *, searched: Iterable[str] | None = None
) -> StepFinding:
    """Grade one step, withholding it when a negative control matches."""
    present = frozenset(str(g) for g in detected)
    looked = frozenset(str(g) for g in searched) if searched is not None else None
    genes = step.requirement.genes
    controls = tuple(g for g in step.negative_controls if g in present)

    if looked is not None and not (set(genes) & looked):
        return StepFinding(step, "not_assayed", (), genes, controls)
    found = tuple(g for g in genes if g in present)
    missing = tuple(g for g in genes if g not in present)
    if controls:
        # The panel matched something it declared it must not match. Reporting
        # the step anyway would publish the false positive.
        return StepFinding(step, "negative_control_matched", found, missing, controls)
    if step.requirement.satisfied_by(present):
        state = "supported"
    elif found:
        state = "partial"
    else:
        state = "absent"
    return StepFinding(step, state, found, missing, controls)


def gaba_balance(synthesis: float | None, degradation: float | None) -> dict[str, Any]:
    """Both components and, only with compatible raw units, their ratio.

    §5.5: a difference between two percentiles is not a net flux. This
    returns the two components always, and a ratio only from raw
    copies-per-100-genomes, with its formula printed beside it.
    """
    out: dict[str, Any] = {
        "synthesis": synthesis,
        "degradation": degradation,
        "unit": "copies per 100 bacterial genomes",
        "ratio": None,
        "formula": "degradation / synthesis, both in copies per 100 bacterial genomes",
        "note": (
            "Two separate measurements, shown side by side. Their difference is not a "
            "net flux and neither is their ratio: both describe how many genes are "
            "present, not how much GABA is made or broken down."
        ),
    }
    if synthesis and degradation is not None and synthesis > 0:
        out["ratio"] = degradation / synthesis
    elif synthesis == 0:
        out["ratio_withheld"] = (
            "no synthesis genes were found, so a ratio would divide by zero and would "
            "say nothing"
        )
    return out


def summarise(
    detected: Iterable[str], *, searched: Iterable[str] | None = None,
    panels: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """A05, A06 and A08 together."""
    genes = list(detected)
    findings = [assess_step(s, genes, searched=searched) for s in ALL_STEPS]
    by_feature: dict[str, list[dict[str, Any]]] = {}
    for finding in findings:
        by_feature.setdefault(finding.step.feature_id, []).append(finding.to_json())
    gaba = (panels or {}).get("gaba") or {}
    return {
        "feature_ids": list(FEATURES),
        "n_steps": len(findings),
        "by_feature": [
            {"feature_id": fid, "feature": FEATURES[fid], "steps": steps}
            for fid, steps in sorted(by_feature.items())
        ],
        "gaba_balance": gaba_balance(
            gaba.get("copies_per_100_genomes"), None,
        ),
        "negative_controls": {
            step.step_id: dict(step.negative_controls)
            for step in ALL_STEPS if step.negative_controls
        },
        "limitations": [
            "A DNA finding is not a metabotype. Whether you actually produce urolithin A "
            "or equol is settled by a substrate challenge with metabolomics, which can "
            "be imported beside this as a separate measurement.",
            "Mucin foraging is physiological; capacity is not an image of an eroded "
            "mucus layer.",
            "Lipid A gene evidence does not resolve a structure, and total Gram-negative "
            "abundance does not resolve a TLR4 stimulus.",
        ],
    }


__all__ = [
    "ALL_STEPS",
    "GLP1_CHANNELS",
    "NOT_MICROBIAL_TARGETS",
    "BY_ID",
    "FEATURES",
    "METHOD_BIOTRANSFORM",
    "NEUROACTIVE",
    "PLANT_COMPOUNDS",
    "STATES",
    "SURFACE",
    "BiotransformError",
    "Requirement",
    "Step",
    "StepFinding",
    "assess_step",
    "gaba_balance",
    "glp1_panel",
    "summarise",
]


# --------------------------------------------------------------------------- #
# A07 — GLP-1-related microbial mechanisms (§9.1)
# --------------------------------------------------------------------------- #

#: Host receptors and enzymes that appear in the GLP-1 literature and are
#: *not* microbial genes. Naming them is the point: a stool metagenome
#: cannot contain FFAR2, and a panel that claims to measure it is
#: measuring something else.
NOT_MICROBIAL_TARGETS: Final[Mapping[str, str]] = {
    "FFAR2": "a human receptor (GPR43). Your cells carry it; your microbes do not.",
    "FFAR3": "a human receptor (GPR41), likewise host rather than microbial.",
    "GCG": "the human proglucagon gene. GLP-1 is its product, made by your L-cells.",
    "TGR5": "a human bile-acid receptor (GPBAR1).",
    "DPP4": (
        "the human enzyme that degrades GLP-1. Some bacteria carry DPP4-like "
        "peptidases, but they are different proteins with different substrates and are "
        "not this."
    ),
}

#: Four channels, kept apart because they act through different chemistry
#: and their evidence is of different strengths.
GLP1_CHANNELS: Final[tuple[Step, ...]] = (
    Step(
        step_id="glp1.scfa_signalling", feature_id="A07",
        label="Short-chain fatty acid signalling",
        requirement=Requirement(any_of=(("but", "buk"), ("pta", "ackA"))),
        establishes=(
            "capacity to make the short-chain fatty acids that are the best-supported "
            "microbial input to L-cell GLP-1 release"
        ),
        does_not_establish=(
            "a GLP-1 level. The receptors that read these signals are yours, not your "
            "microbes', and the step from a fatty acid in the lumen to a hormone in "
            "your blood is host physiology this test does not see"
        ),
        links_to=("butyrate", "propionate"),
        source_ids=("F43",),
        note="The existing butyrate and propionate readings are unchanged.",
    ),
    Step(
        step_id="glp1.bile_acid_signalling", feature_id="A07",
        label="Bile-acid transformation",
        requirement=Requirement(all_of=("bsh",)),
        establishes="capacity to transform bile acids, which changes which ligands reach TGR5",
        does_not_establish=(
            "TGR5 activation. Which bile acids you end up with depends on your own "
            "synthesis and reabsorption as much as on the microbes"
        ),
        links_to=("bsh", "bai"),
        source_ids=("F43",),
    ),
    Step(
        step_id="glp1.indole_signalling", feature_id="A07",
        label="Indole and tryptophan routes",
        requirement=Requirement(any_of=(("tnaA", "fldH"),)),
        establishes="capacity for the indole routes studied in relation to L-cell signalling",
        does_not_establish=(
            "a direction. Indole effects depend on the exposure and the context, and "
            "the literature reports both stimulation and inhibition; no universal "
            "GLP-1-stimulating direction is assigned here"
        ),
        links_to=("indole", "ipa"),
        source_ids=("F43",),
        note="Section 9.1 is explicit that indole routes carry no universal polarity.",
    ),
    Step(
        step_id="glp1.microbial_protein", feature_id="A07",
        label="Experimentally supported microbial proteins",
        # §9.1 requires the detector to be the exact experimentally tested
        # sequence. Amuc_1631 is P9 and is measured by its own panel, which
        # competes it against the protease family it belongs to; clpB is a
        # separate protein with its own satiety literature.
        requirement=Requirement(all_of=("Amuc_1631",)),
        establishes=(
            "presence of a microbial protein with experimental evidence for an effect on "
            "satiety signalling"
        ),
        does_not_establish=(
            "an effect in you. These are single proteins with preclinical evidence, and "
            "presence in a stool metagenome is a long way from a hormone response - it "
            "is not even evidence that the protein is being made, let alone secreted"
        ),
        links_to=("p9",),
        source_ids=("F43", "C02"),
        note=(
            "P9 is detected from the exact sequence that was tested, not from "
            "Akkermansia abundance, and the S41A family it belongs to is counted "
            "against it rather than towards it."
        ),
    ),
)


def glp1_panel(
    detected: Iterable[str], *, searched: Iterable[str] | None = None
) -> dict[str, Any]:
    """The four GLP-1 channels, with the host boundary stated."""
    genes = list(detected)
    findings = [assess_step(step, genes, searched=searched) for step in GLP1_CHANNELS]
    return {
        "feature_id": "A07",
        "n_channels": len(findings),
        "channels": [f.to_json() for f in findings],
        "not_microbial_targets": dict(NOT_MICROBIAL_TARGETS),
        "host_boundary": (
            "Every channel here is a microbial capacity. GLP-1 itself is a human "
            "hormone made by your L-cells from the human proglucagon gene, and the "
            "receptors that read these microbial signals are human receptors. A stool "
            "metagenome contains none of them, so nothing in this panel is a GLP-1 "
            "measurement or a prediction of one."
        ),
        "limitations": [
            "Gene potential, labelled as such, never a host response.",
            "Indole routes carry no universal direction; the evidence goes both ways "
            "and depends on exposure and context.",
            "These channels are separate mechanisms and are not summed into a single "
            "GLP-1 score, because they are not measurements of the same thing.",
        ],
    }
