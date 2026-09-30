"""Dietary substrate panel — A01, BUILD_SPEC_v0.8.3 §5.1.

Fourteen substrate views over one body of enzyme evidence. The hard part
is not finding the enzymes; it is refusing to overclaim from them, and the
specification names the exact overclaims to refuse:

    GH13 alone supports general amylolysis, not every resistant starch.
    GH32 alone does not separate FOS from inulin.
    GH2/GH42 alone do not establish specific GOS utilization.
    GH18 alone does not demonstrate an intact chitin-utilization pathway.

So every substrate declares three things separately: the families that are
**necessary** to act on it at all, the families or context that make the
claim **specific** to this substrate rather than to its whole class, and
the families it **shares** with a neighbouring substrate. A card that has
the first but not the second says so, in those words, instead of quietly
reading as a positive.

The second rule is about arithmetic. A shared gene may appear as evidence
on two cards - that is what "shared" means - but its fragments are counted
once in any aggregate. `aggregate_fragments` is the only supported way to
add these up, and it deduplicates by gene.

Family assignments here are declarations to be checked, not assertions to
be trusted: `validate_against_dbcan` reconciles them with the installed
dbCAN family-substrate mapping and reports every disagreement rather than
silently preferring one side.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

METHOD_SUBSTRATE: Final = "ext083.substrate_capacity/1.0"

#: How specific the evidence on a card is. The middle state is the one the
#: specification exists to protect: enzymes that act on the class, without
#: anything showing they act on *this* member of it.
SPECIFICITY: Final[tuple[str, ...]] = (
    "substrate_specific",   # families or context specific to this substrate
    "class_level_only",     # acts on the class; cannot separate this member
    "insufficient",         # not even the class-level requirement is met
    "not_assayed",          # nothing was searched for
)

#: What a card may claim at each level of specificity.
CLAIM: Final[Mapping[str, str]] = {
    "substrate_specific": (
        "genes whose characterised activity is specific to this substrate were found"
    ),
    "class_level_only": (
        "genes were found that act on this class of carbohydrate, but nothing that "
        "distinguishes this substrate from the others in the class"
    ),
    "insufficient": "the genes required to act on this substrate at all were not found",
    "not_assayed": "this substrate was not searched for in this run",
}


class SubstrateError(ValueError):
    """A substrate definition that would let a class-level finding read as specific."""


@dataclass(frozen=True)
class Substrate:
    """One substrate view and the evidence that would make it specific."""

    substrate_id: str
    label: str
    #: Families without which the organism cannot act on this substrate at
    #: all. Necessary, never sufficient.
    required_families: tuple[str, ...]
    #: What makes the finding about *this* substrate rather than its class.
    #: Either a family only this substrate needs, or a named context such as
    #: a binding module or transporter.
    discriminating_families: tuple[str, ...] = ()
    discriminating_context: tuple[str, ...] = ()
    #: Substrates this one shares required families with, and why sharing is
    #: expected rather than a mistake.
    shares_with: Mapping[str, str] = field(default_factory=dict)
    #: The specific confusion this card must not create.
    must_not_conclude: str = ""
    food_context: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.required_families:
            raise SubstrateError(f"{self.substrate_id}: needs at least one required family")
        if not (self.discriminating_families or self.discriminating_context):
            raise SubstrateError(
                f"{self.substrate_id}: needs something that discriminates it from its class, "
                "or every class-level finding here reads as substrate-specific"
            )
        if self.shares_with and not self.must_not_conclude:
            raise SubstrateError(
                f"{self.substrate_id}: shares families with "
                f"{sorted(self.shares_with)} and must say what that sharing does not prove"
            )

    @property
    def metric_id(self) -> str:
        return f"ext083.{self.substrate_id}"

    def to_json(self) -> dict[str, Any]:
        return {
            "substrate_id": self.substrate_id,
            "metric_id": self.metric_id,
            "label": self.label,
            "required_families": list(self.required_families),
            "discriminating_families": list(self.discriminating_families),
            "discriminating_context": list(self.discriminating_context),
            "shares_with": dict(self.shares_with),
            "must_not_conclude": self.must_not_conclude,
            "food_context": self.food_context,
            "notes": self.notes,
        }


#: The fourteen views of §5.1. Family assignments follow the CAZy
#: classification and are reconciled with the installed dbCAN
#: family-substrate mapping by `validate_against_dbcan`; where the two
#: disagree, the disagreement is reported rather than resolved silently.
SUBSTRATES: Final[tuple[Substrate, ...]] = (
    Substrate(
        substrate_id="carb.cellulose", label="Cellulose",
        required_families=("GH5", "GH9", "GH44", "GH45", "GH48"),
        discriminating_families=("GH48",),
        discriminating_context=("CBM3 cellulose-binding module", "AA10 LPMO", "cellulosome scaffoldin"),
        shares_with={
            "carb.beta_glucan": (
                "GH5 and GH9 include beta-1,4-glucanases that also act on mixed-linkage "
                "beta-glucans"
            ),
        },
        must_not_conclude=(
            "A beta-glucanase is not a cellulase. Cellulose is crystalline beta-1,4-glucan "
            "and needs processive enzymes and a binding module; an enzyme that clears the "
            "soluble mixed-linkage glucan of oats does not touch it."
        ),
        food_context="Vegetables, bran and intact plant material",
    ),
    Substrate(
        substrate_id="carb.starch_general", label="General starch",
        required_families=("GH13", "GH14", "GH15", "GH31", "GH97"),
        discriminating_families=("GH13",),
        discriminating_context=("SusC/SusD-like starch uptake", "GH13 with CBM20/CBM48"),
        shares_with={
            "carb.resistant_starch": "resistant starch is starch; the same amylases open both",
            "carb.imo": "GH13 subfamilies and GH31 also act on alpha-1,6 glucosides",
        },
        must_not_conclude=(
            "Amylolytic capacity is the broad comparator, not evidence about any particular "
            "starch. It says the community can digest starch it can reach."
        ),
        food_context="Starchy foods; cooking and cooling change what reaches the colon",
    ),
    Substrate(
        substrate_id="carb.resistant_starch", label="Resistant starch",
        required_families=("GH13",),
        discriminating_families=("GH13_36", "GH13_39"),
        discriminating_context=(
            "CBM74 starch-binding module", "CBM26 with an amylase",
            "amylosome components (Sas6, Sas20, Amy4)",
        ),
        shares_with={"carb.starch_general": "the same amylase families act on both"},
        must_not_conclude=(
            "GH13 alone supports general amylolysis, not resistant-starch utilisation. "
            "Reaching granular or retrograded starch needs the binding modules and the "
            "surface machinery that a general amylase does not carry, and the physical "
            "type of the starch (RS1 to RS5) is a property of the food rather than of you."
        ),
        food_context="Legumes, unripe banana, cooked-and-cooled starches, named preparations",
        notes="Record the supported physical type where the evidence names one.",
    ),
    Substrate(
        substrate_id="carb.chitin", label="Chitin",
        required_families=("GH18", "GH19"),
        discriminating_families=("GH20",),
        discriminating_context=("CE4 chitin deacetylase", "CBM5/CBM12 chitin-binding module"),
        shares_with={},
        must_not_conclude=(
            "GH18 alone does not demonstrate an intact chitin-utilisation pathway. Many "
            "GH18 proteins act on peptidoglycan or on host glycans, and releasing "
            "chito-oligosaccharides without the downstream hexosaminidase leaves them unused."
        ),
        food_context="Mushrooms and other chitin-containing foods",
    ),
    Substrate(
        substrate_id="carb.pectin", label="Pectin",
        required_families=("GH28", "PL1", "PL9", "PL10", "PL11"),
        discriminating_families=("PL1", "PL9", "CE8"),
        discriminating_context=("CE8 pectin methylesterase", "CE12 rhamnogalacturonan acetylesterase"),
        shares_with={},
        must_not_conclude=(
            "A polygalacturonase acts on the homogalacturonan backbone. The branched "
            "rhamnogalacturonan regions of real pectin need their own enzymes, so backbone "
            "capacity is not whole-pectin capacity."
        ),
        food_context="Apples, citrus, vegetables",
    ),
    Substrate(
        substrate_id="carb.inulin", label="Inulin",
        required_families=("GH32",),
        discriminating_families=("GH91",),
        discriminating_context=(
            "extracellular endo-inulinase activity", "long-chain fructan ABC transporter",
        ),
        shares_with={"carb.fos": "GH32 hydrolyses both; chain length is the difference"},
        must_not_conclude=(
            "GH32 alone does not separate inulin from fructooligosaccharides. Long-chain "
            "inulin has to be cut outside the cell before anything can take it up, so a "
            "GH32 that handles short fructans is not evidence of inulin utilisation."
        ),
        food_context="Chicory root, alliums, artichoke",
    ),
    Substrate(
        substrate_id="carb.fos", label="Fructooligosaccharides",
        required_families=("GH32",),
        discriminating_families=(),
        discriminating_context=(
            "fructooligosaccharide ABC transporter (fosABCDE-like)",
            "intracellular beta-fructofuranosidase with short-chain preference",
        ),
        shares_with={"carb.inulin": "GH32 hydrolyses both; chain length is the difference"},
        must_not_conclude=(
            "Sharing GH32 with inulin is expected and is not double evidence. Short-chain "
            "fructans are taken up whole and cut inside, which is a different route from "
            "the one inulin needs."
        ),
        food_context="Onion, garlic, banana and defined FOS preparations",
    ),
    Substrate(
        substrate_id="carb.gos", label="Galactooligosaccharides",
        required_families=("GH2", "GH42"),
        discriminating_families=("GH53",),
        discriminating_context=(
            "galactooligosaccharide ABC transporter", "LacS-type permease with GH42",
        ),
        shares_with={
            "carb.lactose": "GH2 and GH42 beta-galactosidases hydrolyse both",
        },
        must_not_conclude=(
            "GH2 and GH42 alone do not establish specific GOS utilisation. Every "
            "lactose-using organism carries a beta-galactosidase; taking up an "
            "oligosaccharide of three to eight sugars needs a transporter that lactose "
            "does not."
        ),
        food_context="Defined GOS preparations, which are not the same as legume carbohydrates",
    ),
    Substrate(
        substrate_id="carb.xos", label="Xylooligosaccharides",
        required_families=("GH43", "GH120"),
        discriminating_families=("GH120",),
        discriminating_context=(
            "xylooligosaccharide ABC transporter", "intracellular beta-xylosidase",
        ),
        shares_with={
            "carb.arabinoxylan": "GH43 acts on the xylan backbone and on its oligosaccharides",
        },
        must_not_conclude=(
            "Using xylooligosaccharides does not mean using arabinoxylan. The polymer has "
            "to be cut by an endo-xylanase first, and an organism that lives on the "
            "released oligosaccharides may carry no xylanase at all."
        ),
        food_context="Defined XOS preparations; whole-grain substrates release them",
    ),
    Substrate(
        substrate_id="carb.imo", label="Isomaltooligosaccharides",
        required_families=("GH13", "GH31"),
        discriminating_families=("GH13_31",),
        discriminating_context=("oligo-1,6-glucosidase activity", "alpha-glucoside transporter"),
        shares_with={"carb.starch_general": "GH13 and GH31 also act on alpha-1,4 starch"},
        must_not_conclude=(
            "An alpha-glucosidase is not specific to the alpha-1,6 linkages that define "
            "IMO. Commercial IMO preparations also differ widely in how much of them "
            "actually resists digestion, which is a property of the product."
        ),
        food_context="Defined IMO preparations; formulation matters",
    ),
    Substrate(
        substrate_id="carb.lactose", label="Lactose",
        required_families=("GH2", "GH42", "GH35"),
        discriminating_families=("GH35",),
        discriminating_context=("LacS or LacY lactose permease", "lactose phosphotransferase system"),
        shares_with={"carb.gos": "the same beta-galactosidases hydrolyse both"},
        must_not_conclude=(
            "This is the community's capacity to use lactose that reaches the colon. It is "
            "not a lactase-persistence genotype, not a breath test, and not an explanation "
            "of dairy symptoms, which depend on your own small-intestinal lactase."
        ),
        food_context="Dairy lactose; fermented dairy carries less",
    ),
    Substrate(
        substrate_id="carb.beta_glucan", label="Beta-glucans",
        required_families=("GH16", "GH3"),
        discriminating_families=("GH55", "GH64", "GH81", "GH128"),
        discriminating_context=(
            "licheninase (beta-1,3-1,4) activity for cereal glucan",
            "beta-1,3-glucanase for fungal glucan",
        ),
        shares_with={
            "carb.cellulose": "some GH5 and GH9 enzymes act on both beta-glucan linkages",
        },
        must_not_conclude=(
            "Cereal and fungal beta-glucans are different molecules with different "
            "linkages, and the enzymes are not interchangeable. An oat beta-glucan finding "
            "says nothing about a mushroom one."
        ),
        food_context="Oats and barley, versus fungal preparations",
    ),
    Substrate(
        substrate_id="carb.arabinoxylan", label="Arabinoxylan",
        required_families=("GH10", "GH11", "GH43", "GH51"),
        discriminating_families=("GH11", "GH51"),
        discriminating_context=(
            "CE1 feruloyl esterase", "CE6 acetyl xylan esterase",
            "arabinofuranosidase acting on doubly substituted xylose",
        ),
        shares_with={"carb.xos": "GH43 acts on the backbone and on the oligosaccharides"},
        must_not_conclude=(
            "A xylanase opens the backbone; the arabinose side chains and their ferulate "
            "cross-links are what make cereal arabinoxylan hard to use, and they need "
            "their own enzymes."
        ),
        food_context="Wheat and rye bran, psyllium and specified extracts",
    ),
    Substrate(
        substrate_id="carb.galactomannan", label="Galactomannan and guar",
        required_families=("GH26", "GH5_8"),
        discriminating_families=("GH26", "GH113"),
        discriminating_context=("GH27 alpha-galactosidase for the side chains",
                                "beta-mannosidase for the released oligosaccharides"),
        shares_with={},
        must_not_conclude=(
            "The mannan backbone and its galactose side chains need different enzymes, and "
            "partially hydrolysed guar behaves differently from the intact gum, which is a "
            "property of the preparation rather than of the community."
        ),
        food_context="Guar gum and partially hydrolysed guar gum",
    ),
)
BY_ID: Final[Mapping[str, Substrate]] = {s.substrate_id: s for s in SUBSTRATES}


# --------------------------------------------------------------------------- #
# evidence on one card
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class GeneHit:
    """One measured gene, with the family it belongs to."""

    gene_id: str
    family: str
    fragments: int
    rpkm: float | None = None


def families_found(hits: Iterable[GeneHit]) -> frozenset[str]:
    """Every family present, including the parent of any subfamily.

    `GH13_36` is a GH13, so a subfamily hit satisfies a GH13 requirement.
    The reverse never holds, which is the whole point of the subfamily.
    """
    out: set[str] = set()
    for hit in hits:
        family = str(hit.family)
        out.add(family)
        if "_" in family:
            out.add(family.split("_", 1)[0])
    return frozenset(out)


def aggregate_fragments(hits: Iterable[GeneHit]) -> int:
    """Total fragments, counting each gene once.

    §5.1: a shared gene may appear as evidence on two substrate cards, but
    its fragments cannot be counted twice in an aggregate. Summing the
    cards would do exactly that, so this is the only supported way to add
    them up.
    """
    by_gene: dict[str, int] = {}
    for hit in hits:
        by_gene[hit.gene_id] = max(by_gene.get(hit.gene_id, 0), int(hit.fragments))
    return sum(by_gene.values())


@dataclass(frozen=True)
class SubstrateFinding:
    """What one substrate card may and may not say."""

    substrate: Substrate
    specificity: str
    families_present: tuple[str, ...]
    required_present: tuple[str, ...]
    required_missing: tuple[str, ...]
    discriminators_present: tuple[str, ...]
    fragments: int
    genes: tuple[str, ...]
    assayed: bool = True

    @property
    def claim(self) -> str:
        return CLAIM[self.specificity]

    def to_json(self) -> dict[str, Any]:
        return {
            "method_id": METHOD_SUBSTRATE,
            **self.substrate.to_json(),
            "specificity": self.specificity,
            "claim": self.claim,
            "families_present": list(self.families_present),
            "required_present": list(self.required_present),
            "required_missing": list(self.required_missing),
            "discriminators_present": list(self.discriminators_present),
            "fragments": self.fragments,
            "n_genes": len(self.genes),
            "genes": list(self.genes),
            "aggregate_rule": (
                "A gene shared with another substrate appears on both cards and is counted "
                "once in any total."
            ),
            "limitations": [
                "Genetic capacity, not a measured amount of anything digested.",
                self.substrate.must_not_conclude,
            ],
        }


def assess(
    substrate: Substrate, hits: Sequence[GeneHit], *, assayed: bool = True
) -> SubstrateFinding:
    """Grade one substrate's evidence, refusing to let class pass as specific."""
    present = families_found(hits)
    required_present = tuple(f for f in substrate.required_families if f in present)
    required_missing = tuple(f for f in substrate.required_families if f not in present)
    discriminators = tuple(f for f in substrate.discriminating_families if f in present)

    if not assayed:
        specificity = "not_assayed"
    elif not required_present:
        specificity = "insufficient"
    elif discriminators or (
        substrate.discriminating_context and not substrate.discriminating_families
        and len(required_present) == len(substrate.required_families)
    ):
        # Where a substrate has no discriminating family of its own, the
        # discrimination is contextual (a transporter, a binding module) and
        # cannot be settled from family membership alone. Requiring the full
        # complement is the strictest thing family evidence can support.
        specificity = "substrate_specific" if discriminators else "class_level_only"
    else:
        specificity = "class_level_only"

    return SubstrateFinding(
        substrate=substrate, specificity=specificity,
        families_present=tuple(sorted(present)),
        required_present=required_present, required_missing=required_missing,
        discriminators_present=discriminators,
        fragments=aggregate_fragments(hits),
        genes=tuple(sorted({h.gene_id for h in hits})),
        assayed=assayed,
    )


def assess_all(
    hits_by_family: Mapping[str, Sequence[GeneHit]], *, assayed: bool = True
) -> list[SubstrateFinding]:
    """Every substrate card from one body of gene evidence."""
    out: list[SubstrateFinding] = []
    for substrate in SUBSTRATES:
        wanted = set(substrate.required_families) | set(substrate.discriminating_families)
        hits = [
            hit for family, family_hits in hits_by_family.items()
            for hit in family_hits
            if family in wanted or family.split("_", 1)[0] in wanted
        ]
        out.append(assess(substrate, hits, assayed=assayed))
    return out


def shared_gene_report(findings: Sequence[SubstrateFinding]) -> dict[str, Any]:
    """Which genes support more than one card, and the arithmetic that follows."""
    by_gene: dict[str, list[str]] = {}
    for finding in findings:
        for gene in finding.genes:
            by_gene.setdefault(gene, []).append(finding.substrate.substrate_id)
    shared = {g: cards for g, cards in by_gene.items() if len(cards) > 1}
    naive = sum(f.fragments for f in findings)
    return {
        "n_genes": len(by_gene),
        "n_shared_genes": len(shared),
        "shared": {g: sorted(set(cards)) for g, cards in sorted(shared.items())},
        "sum_of_cards": naive,
        "note": (
            "Adding the cards together would count every shared gene once per card. "
            "Any total over these substrates must be taken with `aggregate_fragments`, "
            "which counts each gene once."
        ),
    }


# --------------------------------------------------------------------------- #
# reconciliation with the installed dbCAN mapping
# --------------------------------------------------------------------------- #

#: Words in the dbCAN substrate column that name each of our substrates.
DBCAN_SUBSTRATE_WORDS: Final[Mapping[str, tuple[str, ...]]] = {
    "carb.cellulose": ("cellulose",),
    "carb.starch_general": ("starch",),
    "carb.resistant_starch": ("resistant starch", "starch"),
    "carb.chitin": ("chitin",),
    "carb.pectin": ("pectin", "homogalacturonan", "rhamnogalacturonan"),
    "carb.inulin": ("inulin", "fructan"),
    "carb.fos": ("fructan", "fructooligosaccharide", "levan"),
    "carb.gos": ("galactan", "galactooligosaccharide", "lactose"),
    "carb.xos": ("xylan", "xylooligosaccharide"),
    "carb.imo": ("starch", "dextran", "isomaltose"),
    "carb.lactose": ("lactose", "galactan"),
    "carb.beta_glucan": ("beta-glucan", "glucan", "lichenan", "curdlan"),
    "carb.arabinoxylan": ("arabinoxylan", "xylan", "arabinan"),
    "carb.galactomannan": ("mannan", "galactomannan", "glucomannan"),
}

_SPLIT = re.compile(r"[|,;]+")


def load_dbcan_mapping(path: Path) -> dict[str, frozenset[str]]:
    """family -> substrate words, from dbCAN's `fam-substrate-mapping.tsv`."""
    out: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            family_col = next(
                i for i, name in enumerate(header) if "family" in name.lower()
            )
            substrate_col = next(
                i for i, name in enumerate(header) if "substrate" in name.lower()
            )
        except StopIteration as exc:
            raise SubstrateError(
                f"{path}: expected a family column and a substrate column, got {header}"
            ) from exc
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(family_col, substrate_col):
                continue
            family = parts[family_col].strip()
            if not family:
                continue
            words = {
                w.strip().lower() for w in _SPLIT.split(parts[substrate_col]) if w.strip()
            }
            out.setdefault(family, set()).update(words)
    return {k: frozenset(v) for k, v in out.items()}


def validate_against_dbcan(mapping: Mapping[str, frozenset[str]]) -> dict[str, Any]:
    """Reconcile the declared families with the installed dbCAN mapping.

    Disagreements are reported, not resolved. A family this module claims
    for a substrate that dbCAN does not associate with it is a question for
    a human, and silently preferring either source would hide it.
    """
    unknown: dict[str, list[str]] = {}
    unsupported: dict[str, list[str]] = {}
    for substrate in SUBSTRATES:
        words = DBCAN_SUBSTRATE_WORDS.get(substrate.substrate_id, ())
        for family in (*substrate.required_families, *substrate.discriminating_families):
            base = family.split("_", 1)[0]
            declared = mapping.get(family) or mapping.get(base)
            if declared is None:
                unknown.setdefault(substrate.substrate_id, []).append(family)
                continue
            if words and not any(
                any(word in entry for entry in declared) for word in words
            ):
                unsupported.setdefault(substrate.substrate_id, []).append(family)
    return {
        "n_families_in_mapping": len(mapping),
        "families_not_in_mapping": {k: sorted(v) for k, v in sorted(unknown.items())},
        "families_dbcan_does_not_link_to_this_substrate": {
            k: sorted(v) for k, v in sorted(unsupported.items())
        },
        "agrees": not unknown and not unsupported,
        "note": (
            "Disagreements are reported rather than resolved. A family declared here that "
            "the installed dbCAN release does not associate with this substrate is a "
            "question for a person, not something to silently drop or keep."
        ),
    }


__all__ = [
    "BY_ID",
    "CLAIM",
    "DBCAN_SUBSTRATE_WORDS",
    "METHOD_SUBSTRATE",
    "SPECIFICITY",
    "SUBSTRATES",
    "GeneHit",
    "Substrate",
    "SubstrateError",
    "SubstrateFinding",
    "aggregate_fragments",
    "assess",
    "assess_all",
    "families_found",
    "load_dbcan_mapping",
    "shared_gene_report",
    "validate_against_dbcan",
]
