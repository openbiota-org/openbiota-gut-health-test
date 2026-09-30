"""Protein, nitrogen and aromatic metabolism — A04, BUILD_SPEC_v0.8.3 §5.4.

Six modules over one question a reader asks badly: "is my gut breaking
down too much protein?" The honest answer needs several things kept
apart that are easy to blur, and the specification names each of them.

**Protein breakdown, BCAA synthesis and BCAA fermentation are three
different concepts.** The existing `bcaa` panel measures the second. This
module measures the third, and they move in opposite directions: a
community that *makes* branched-chain amino acids is not one that
*ferments* them into branched-chain fatty acids. The module refuses to
share a metric ID with the existing panel and says so on the card.

**Generic peptidase abundance is not dietary protein digestion.** Every
organism carries peptidases for its own turnover. A count of them
measures the census, not your dinner, so the peptidase module reports
localisation and substrate specificity where they are established and
declines to aggregate where they are not.

**A microbial precursor is not a host metabolite.** Microbes make
phenylacetate; the *host* conjugates it to phenylacetylglutamine. Microbes
make TMA; the host oxidises it to TMAO. Both relationships are carried as
explicit links with the host step named, never as one measurement.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, Final

METHOD_NITROGEN: Final = "ext083.nitrogen_metabolism/1.0"

MODULES: Final[tuple[tuple[str, str], ...]] = (
    ("proteolysis", "Protein and peptide breakdown"),
    ("peptide_transport", "Peptide uptake and amino-acid use"),
    ("amino_acid_fermentation", "Amino-acid fermentation, including Stickland pairs"),
    ("bcfa", "Branched-chain fatty acid formation"),
    ("ammonia", "Ammonia generation and assimilation"),
    ("aromatic", "Aromatic precursor routes"),
)
MODULE_LABELS: Final[Mapping[str, str]] = dict(MODULES)


class NitrogenError(ValueError):
    """A definition that would blur two of the concepts this module separates."""


@dataclass(frozen=True)
class HostStep:
    """A reaction the host performs on a microbial product.

    Named so that a microbial precursor is never printed as the
    circulating metabolite it is a precursor to.
    """

    microbial_product: str
    host_enzyme: str
    host_product: str
    note: str

    def to_json(self) -> dict[str, Any]:
        return {
            "microbial_product": self.microbial_product,
            "host_enzyme": self.host_enzyme,
            "host_product": self.host_product,
            "note": self.note,
            "measured_here": self.microbial_product,
            "not_measured_here": self.host_product,
        }


#: The two host steps §5.4 insists on carrying as links rather than merging.
HOST_STEPS: Final[tuple[HostStep, ...]] = (
    HostStep(
        microbial_product="phenylacetate",
        host_enzyme="hepatic glutamine N-acyltransferase",
        host_product="phenylacetylglutamine (PAGln)",
        note=(
            "Microbes make phenylacetate; your liver conjugates it to PAGln. The "
            "circulating metabolite associated with cardiovascular outcomes is the host "
            "product, and this measures only the microbial precursor."
        ),
    ),
    HostStep(
        microbial_product="trimethylamine (TMA)",
        host_enzyme="hepatic FMO3",
        host_product="trimethylamine N-oxide (TMAO)",
        note=(
            "Microbes make TMA; your liver oxidises it to TMAO. FMO3 activity varies "
            "widely between people, so microbial TMA capacity does not fix a TMAO level."
        ),
    ),
)


@dataclass(frozen=True)
class NitrogenModule:
    """One module, its genes, and the conclusion it must not support."""

    module_id: str
    label: str
    genes: tuple[str, ...]
    #: Genes whose presence is expected in every organism and therefore
    #: carries little information about diet or disease.
    housekeeping: tuple[str, ...] = ()
    #: What this module is *not*, named against the concept it is confused
    #: with most often.
    distinct_from: Mapping[str, str] = field(default_factory=dict)
    #: Existing report readings this one links to without replacing.
    links_to: tuple[str, ...] = ()
    host_steps: tuple[HostStep, ...] = ()
    aggregatable: bool = True
    note: str = ""

    def __post_init__(self) -> None:
        if not self.genes:
            raise NitrogenError(f"{self.module_id}: a module needs genes")
        if self.module_id not in MODULE_LABELS:
            raise NitrogenError(f"{self.module_id}: unknown module")
        unknown = set(self.housekeeping) - set(self.genes)
        if unknown:
            raise NitrogenError(
                f"{self.module_id}: {sorted(unknown)} marked housekeeping but not in the "
                "module's genes"
            )

    @property
    def informative_genes(self) -> tuple[str, ...]:
        """Genes whose presence says something beyond "a cell lives here"."""
        return tuple(g for g in self.genes if g not in self.housekeeping)

    def to_json(self) -> dict[str, Any]:
        return {
            "module_id": self.module_id, "label": self.label,
            "genes": list(self.genes),
            "housekeeping_genes": list(self.housekeeping),
            "informative_genes": list(self.informative_genes),
            "distinct_from": dict(self.distinct_from),
            "links_to": list(self.links_to),
            "host_steps": [h.to_json() for h in self.host_steps],
            "aggregatable": self.aggregatable,
            "note": self.note,
        }


MODULE_DEFINITIONS: Final[tuple[NitrogenModule, ...]] = (
    NitrogenModule(
        module_id="proteolysis", label="Protein and peptide breakdown",
        genes=("pepA", "pepD", "pepN", "pepT", "clpP", "lon", "dppA", "sspA"),
        housekeeping=("clpP", "lon", "pepA", "pepN"),
        aggregatable=False,
        distinct_from={
            "dietary protein digestion": (
                "Every organism carries peptidases for its own protein turnover. Counting "
                "them measures how many organisms are present, not how much of your "
                "dinner reaches the colon undigested."
            ),
        },
        note=(
            "Reported by MEROPS family with localisation and substrate specificity where "
            "they are established. Not summed into one proteolysis score, because the "
            "housekeeping enzymes would dominate any such total."
        ),
    ),
    NitrogenModule(
        module_id="peptide_transport", label="Peptide uptake and amino-acid use",
        genes=("oppA", "oppB", "oppD", "dtpT", "livK", "brnQ"),
        distinct_from={
            "proteolysis": (
                "Taking peptides in is a separate capability from cutting proteins up, "
                "and an organism can do either without the other."
            ),
        },
    ),
    NitrogenModule(
        module_id="amino_acid_fermentation",
        label="Amino-acid fermentation, including Stickland pairs",
        genes=("grdA", "grdB", "prdA", "prdB", "fldH", "acdA"),
        distinct_from={
            "amino acid synthesis": (
                "Fermenting an amino acid for energy is the opposite of making it. The "
                "two capacities can coexist in a community and say different things."
            ),
        },
        note=(
            "Stickland fermentation pairs a donor and an acceptor amino acid; the "
            "reductive (grd, prd) and oxidative branches are reported together because "
            "neither runs alone."
        ),
    ),
    NitrogenModule(
        module_id="bcfa", label="Branched-chain fatty acid formation",
        genes=("ilvE", "kivD", "bkdA", "bkdB"),
        links_to=("bcaa",),
        distinct_from={
            "BCAA biosynthesis": (
                "The existing BCAA panel measures the capacity to *make* branched-chain "
                "amino acids. This measures the capacity to *ferment* them into "
                "branched-chain fatty acids. They are different pathways running in "
                "opposite directions, and a community can be high in both."
            ),
        },
        note=(
            "Isobutyrate and isovalerate come from valine and leucine fermentation, which "
            "is a protein-fermentation signal rather than a fibre one."
        ),
    ),
    NitrogenModule(
        module_id="ammonia", label="Ammonia generation and assimilation",
        genes=("gdhA", "gdhB", "glnA", "gltB", "asnB", "ansA", "aspA"),
        links_to=("urease",),
        distinct_from={
            "urease": (
                "Urease is its own existing reading and is unchanged. These are the "
                "non-urease routes to ammonia, plus the routes that assimilate it back."
            ),
        },
        note=(
            "Generation and assimilation are reported separately: a community that makes "
            "ammonia and reuses it is not the same as one that only makes it."
        ),
    ),
    NitrogenModule(
        module_id="aromatic", label="Aromatic precursor routes",
        genes=("ppdC", "fldC", "hpdB", "padI", "porA"),
        links_to=("indole", "ipa", "pcresol"),
        host_steps=HOST_STEPS,
        distinct_from={
            "the circulating metabolite": (
                "These are microbial precursors. The host performs a further step on "
                "several of them, and the blood metabolite is the host's product."
            ),
        },
        note=(
            "Linked to the existing indole, IPA and p-cresol readings without merging "
            "them: they are chemically distinct endpoints and stay separate."
        ),
    ),
)
BY_ID: Final[Mapping[str, NitrogenModule]] = {m.module_id: m for m in MODULE_DEFINITIONS}


# --------------------------------------------------------------------------- #
# grading
# --------------------------------------------------------------------------- #

STATES: Final[tuple[str, ...]] = ("present", "partial", "absent", "not_assayed")


@dataclass(frozen=True)
class ModuleFinding:
    module: NitrogenModule
    state: str
    present: tuple[str, ...]
    informative_present: tuple[str, ...]
    missing: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        payload = {
            "method_id": METHOD_NITROGEN, **self.module.to_json(),
            "state": self.state,
            "genes_present": list(self.present),
            "informative_genes_present": list(self.informative_present),
            "genes_missing": list(self.missing),
            "limitations": [
                "Genetic capacity, not a measured amount and not a dietary intake.",
                *self.module.distinct_from.values(),
            ],
        }
        if not self.module.aggregatable:
            payload["no_aggregate_reason"] = (
                "This module is not summed into a single score; its housekeeping members "
                "would dominate the total and the number would track biomass rather than "
                "the capability."
            )
        return payload


def assess_module(
    module: NitrogenModule, detected: Iterable[str], *,
    searched: Iterable[str] | None = None,
) -> ModuleFinding:
    """Grade one module, counting only genes that carry information."""
    found = frozenset(str(g) for g in detected)
    looked = frozenset(str(g) for g in searched) if searched is not None else None
    if looked is not None and not (set(module.genes) & looked):
        return ModuleFinding(module, "not_assayed", (), (), module.genes)
    present = tuple(g for g in module.genes if g in found)
    informative = tuple(g for g in module.informative_genes if g in found)
    missing = tuple(g for g in module.genes if g not in found)
    if informative and len(present) == len(module.genes):
        state = "present"
    elif informative:
        state = "partial"
    elif present:
        # Only housekeeping genes found: that is the census, not the capability.
        state = "partial"
    else:
        state = "absent"
    return ModuleFinding(module, state, present, informative, missing)


def summarise(
    detected: Iterable[str], *, searched: Iterable[str] | None = None
) -> dict[str, Any]:
    """The whole A04 view."""
    genes = list(detected)
    findings = [assess_module(m, genes, searched=searched) for m in MODULE_DEFINITIONS]
    return {
        "feature_id": "A04",
        "n_modules": len(findings),
        "modules": [f.to_json() for f in findings],
        "three_different_concepts": {
            "protein_breakdown": "cutting proteins and peptides up (proteolysis)",
            "bcaa_synthesis": "making branched-chain amino acids (the existing bcaa panel)",
            "bcaa_fermentation": "fermenting them into branched-chain fatty acids (bcfa)",
            "note": (
                "These are three different things and this report keeps three separate "
                "readings for them. A community can be high in all three."
            ),
        },
        "host_steps": [h.to_json() for h in HOST_STEPS],
        "preserved_unchanged": ["urease", "indole", "ipa", "pcresol", "bcaa", "cutc"],
        "limitations": [
            "Generic peptidase abundance is not a measurement of dietary protein "
            "digestion; every organism carries peptidases for its own turnover.",
            "A microbial precursor is not the circulating metabolite the host makes "
            "from it.",
        ],
    }


__all__ = [
    "BY_ID",
    "HOST_STEPS",
    "METHOD_NITROGEN",
    "MODULES",
    "MODULE_DEFINITIONS",
    "MODULE_LABELS",
    "STATES",
    "HostStep",
    "ModuleFinding",
    "NitrogenError",
    "NitrogenModule",
    "assess_module",
    "summarise",
]
