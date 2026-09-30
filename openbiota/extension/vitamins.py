"""Nine-vitamin dashboard — A03, BUILD_SPEC_v0.8.3 §5.3.

Four new vitamins (B1, B3, B5, B6) reconstructed here; five existing ones
(B2, B7, B9, B12, K2) reused unchanged with a route-detail view only.

Three distinctions the specification insists on, each enforced rather
than described:

**A transport gene is not synthesis.** Uptake, salvage and de-novo
synthesis are three columns, never summed. An organism that imports
thiamine cannot make it, and a dashboard that adds the columns together
says it can.

**An alternative route is an alternative.** B6 has a DXP-dependent
(PdxA/PdxJ) and a DXP-independent (PdxS/PdxT) pathway. Requiring both
would report nearly every community as unable to make B6; scoring one
isolated homolog as a complete route would report nearly all of them as
able to. Each route is graded on its own genes and the vitamin takes the
best complete one.

**A branch is not a pathway.** B1's thiazole and pyrimidine branches
must both be present *and* be coupled; B5 needs pantoate and
beta-alanine *and* the ligase. Downstream use of the product - CoA from
pantothenate - is not a second way of producing it.

Nothing here is a serum level. It is what the organisms present could
make, which is a different quantity from what you absorb, and it is
never a reason to stop a prescribed supplement.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Final

METHOD_VITAMIN: Final = "ext083.vitamin_capability/1.0"

#: The three columns. They answer different questions and are never added.
COLUMNS: Final[tuple[tuple[str, str], ...]] = (
    ("synthesis", "Can make it from simpler precursors"),
    ("salvage", "Can rebuild it from a partial form it takes in"),
    ("uptake", "Can take the finished vitamin in, but not make it"),
)
COLUMN_LABELS: Final[Mapping[str, str]] = dict(COLUMNS)

#: How complete a route is.
ROUTE_STATES: Final[tuple[str, ...]] = (
    "complete", "partial", "absent", "not_assayed", "not_applicable",
)


class VitaminError(ValueError):
    """A definition that would let uptake read as synthesis."""


@dataclass(frozen=True)
class Branch:
    """One arm of a pathway, which is not a pathway on its own."""

    branch_id: str
    label: str
    genes: tuple[str, ...]
    note: str = ""

    def __post_init__(self) -> None:
        if not self.genes:
            raise VitaminError(f"{self.branch_id}: a branch needs genes")

    def state(
        self, present: frozenset[str], *, assayed: bool = True,
        searched: frozenset[str] | None = None,
    ) -> str:
        # "Absent" means looked for and not found. A gene no panel searches
        # for is `not_assayed`, and calling it absent would report a gap in
        # the assay as a property of the person.
        if not assayed or (searched is not None and not (set(self.genes) & searched)):
            return "not_assayed"
        found = [g for g in self.genes if g in present]
        if len(found) == len(self.genes):
            return "complete"
        return "partial" if found else "absent"

    def to_json(self, present: frozenset[str] | None = None) -> dict[str, Any]:
        payload = {
            "branch_id": self.branch_id, "label": self.label,
            "genes": list(self.genes), "note": self.note,
        }
        if present is not None:
            payload["genes_present"] = [g for g in self.genes if g in present]
            payload["genes_missing"] = [g for g in self.genes if g not in present]
            payload["state"] = self.state(present)
        return payload


@dataclass(frozen=True)
class VitaminRoute:
    """One complete way to arrive at a vitamin.

    A route is complete only when every branch is complete *and* the
    coupling step is present. Branches without coupling make the
    intermediates and stop.
    """

    route_id: str
    label: str
    column: str
    branches: tuple[Branch, ...]
    coupling: tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        if self.column not in COLUMN_LABELS:
            raise VitaminError(f"{self.route_id}: unknown column {self.column!r}")
        if not self.branches:
            raise VitaminError(f"{self.route_id}: a route needs at least one branch")
        if self.column == "synthesis" and len(self.branches) > 1 and not self.coupling:
            raise VitaminError(
                f"{self.route_id}: a multi-branch synthesis route must name its coupling "
                "step, or two half-pathways would read as a whole one"
            )

    @property
    def all_genes(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(
            [g for b in self.branches for g in b.genes] + list(self.coupling)
        ))

    def state(
        self, present: frozenset[str], *, assayed: bool = True,
        searched: frozenset[str] | None = None,
    ) -> str:
        if not assayed:
            return "not_assayed"
        if searched is not None and not (set(self.all_genes) & searched):
            return "not_assayed"
        branch_states = [b.state(present, searched=searched) for b in self.branches]
        coupled = all(g in present for g in self.coupling)
        if all(s == "complete" for s in branch_states) and coupled:
            return "complete"
        if any(s in {"complete", "partial"} for s in branch_states):
            return "partial"
        if all(s == "not_assayed" for s in branch_states):
            return "not_assayed"
        return "absent"

    def to_json(self, present: frozenset[str] | None = None) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "route_id": self.route_id, "label": self.label,
            "column": self.column, "column_meaning": COLUMN_LABELS[self.column],
            "branches": [b.to_json(present) for b in self.branches],
            "coupling_genes": list(self.coupling),
            "note": self.note,
        }
        if present is not None:
            payload["state"] = self.state(present)
            payload["coupling_present"] = [g for g in self.coupling if g in present]
            payload["coupling_missing"] = [g for g in self.coupling if g not in present]
        return payload


@dataclass(frozen=True)
class Vitamin:
    """One vitamin and every way this dashboard knows of arriving at it."""

    vitamin_id: str
    label: str
    common_name: str
    routes: tuple[VitaminRoute, ...]
    #: Reused from the existing panels rather than recomputed here.
    reused_panel: str | None = None
    must_not_conclude: str = ""

    def __post_init__(self) -> None:
        if not self.routes:
            raise VitaminError(f"{self.vitamin_id}: needs at least one route")
        if not any(r.column == "synthesis" for r in self.routes):
            raise VitaminError(
                f"{self.vitamin_id}: has no synthesis route, so its dashboard could only "
                "report uptake as though it were production"
            )

    @property
    def metric_id(self) -> str:
        return f"ext083.vitamin.{self.vitamin_id}"

    def routes_in(self, column: str) -> tuple[VitaminRoute, ...]:
        return tuple(r for r in self.routes if r.column == column)


#: B1. Two branches that have to be coupled, plus salvage kept separate.
_B1 = Vitamin(
    vitamin_id="b1", label="Vitamin B1", common_name="thiamine",
    must_not_conclude=(
        "Taking thiamine in is not making it. An organism with only the transporter "
        "depends on thiamine that is already there."
    ),
    routes=(
        VitaminRoute(
            route_id="b1.denovo", label="De-novo thiamine synthesis", column="synthesis",
            branches=(
                Branch("b1.thiazole", "Thiazole branch", ("thiG", "thiH", "thiS"),
                       note="Forms the thiazole ring."),
                Branch("b1.pyrimidine", "Hydroxymethylpyrimidine branch", ("thiC", "thiD"),
                       note="Forms the pyrimidine ring."),
            ),
            coupling=("thiE",),
            note=(
                "Both rings and the thiamine phosphate synthase that joins them. Either "
                "branch alone makes an intermediate and stops."
            ),
        ),
        VitaminRoute(
            route_id="b1.salvage", label="Thiamine and precursor salvage", column="salvage",
            branches=(Branch("b1.salvage.genes", "Salvage kinases", ("thiM", "thiD")),),
            note="Rebuilds thiamine from thiazole or pyrimidine taken in from outside.",
        ),
        VitaminRoute(
            route_id="b1.uptake", label="Thiamine uptake", column="uptake",
            branches=(Branch("b1.transport", "ThiBPQ / PnuT transport", ("thiB", "thiP", "thiQ")),),
            note="Transport only. This is the column that must never read as synthesis.",
        ),
    ),
)

#: B3. De-novo and salvage are genuinely different chemistries here.
_B3 = Vitamin(
    vitamin_id="b3", label="Vitamin B3", common_name="niacin and NAD cofactors",
    must_not_conclude=(
        "Salvaging nicotinamide is not making the ring. The de-novo aspartate route and "
        "the salvage route arrive at NAD by different chemistry and are reported apart."
    ),
    routes=(
        VitaminRoute(
            route_id="b3.denovo", label="De-novo (aspartate to quinolinate to NAD)",
            column="synthesis",
            branches=(
                Branch("b3.quinolinate", "Aspartate to quinolinate", ("nadA", "nadB")),
                Branch("b3.namn", "Quinolinate to NaMN and NAD", ("nadC", "nadD", "nadE")),
            ),
            coupling=("nadC",),
            note="The ring is built from aspartate rather than recovered.",
        ),
        VitaminRoute(
            route_id="b3.salvage", label="Nicotinamide and nicotinate salvage",
            column="salvage",
            branches=(Branch("b3.salvage.genes", "Salvage route", ("pncA", "pncB")),),
            note="Recovers the intact pyridine ring; it does not make one.",
        ),
        VitaminRoute(
            route_id="b3.uptake", label="Niacin uptake", column="uptake",
            branches=(Branch("b3.transport", "PnuC-type transport", ("pnuC",)),),
        ),
    ),
)

#: B5. Two branches, a ligase, and the CoA trap called out explicitly.
_B5 = Vitamin(
    vitamin_id="b5", label="Vitamin B5", common_name="pantothenate",
    must_not_conclude=(
        "Using pantothenate to make coenzyme A is downstream consumption, not a second "
        "way of producing pantothenate. Counting coaA as production would turn every "
        "organism into a producer."
    ),
    routes=(
        VitaminRoute(
            route_id="b5.denovo", label="De-novo pantothenate synthesis", column="synthesis",
            branches=(
                Branch("b5.pantoate", "Pantoate branch", ("panB", "panE")),
                Branch("b5.beta_alanine", "Beta-alanine branch", ("panD",)),
            ),
            coupling=("panC",),
            note="Pantoate and beta-alanine, joined by pantothenate synthetase.",
        ),
        VitaminRoute(
            route_id="b5.uptake", label="Pantothenate uptake", column="uptake",
            branches=(Branch("b5.transport", "PanF / PanT transport", ("panF",)),),
        ),
    ),
)

#: B6. Two documented alternatives; neither is required of the other.
_B6 = Vitamin(
    vitamin_id="b6", label="Vitamin B6", common_name="pyridoxine and pyridoxal phosphate",
    must_not_conclude=(
        "The DXP-dependent and DXP-independent routes are alternatives. Requiring both "
        "would call almost every community unable to make B6; accepting one isolated "
        "homolog would call almost all of them able to."
    ),
    routes=(
        VitaminRoute(
            route_id="b6.dxp_dependent", label="DXP-dependent route (PdxA/PdxJ)",
            column="synthesis",
            branches=(Branch("b6.dxp", "PdxA and PdxJ", ("pdxA", "pdxJ")),),
            note="The route found in E. coli and its relatives.",
        ),
        VitaminRoute(
            route_id="b6.dxp_independent", label="DXP-independent route (PdxS/PdxT)",
            column="synthesis",
            branches=(Branch("b6.pdxst", "PdxS and PdxT", ("pdxS", "pdxT")),),
            note=(
                "The route found in most other bacteria. Both subunits are needed: PdxS "
                "alone is a glutaminase-dependent synthase missing its glutaminase."
            ),
        ),
        VitaminRoute(
            route_id="b6.salvage", label="B6 salvage", column="salvage",
            branches=(Branch("b6.salvage.genes", "Pyridoxal kinase route", ("pdxK", "pdxH")),),
        ),
    ),
)

NEW_VITAMINS: Final[tuple[Vitamin, ...]] = (_B1, _B3, _B5, _B6)

#: Reused unchanged. The dashboard shows their existing value and adds a
#: route-detail view; §5.3 is explicit that these are not rescored.
REUSED: Final[tuple[dict[str, str], ...]] = (
    {"vitamin_id": "b2", "label": "Vitamin B2", "common_name": "riboflavin",
     "panel": "riboflavin",
     "route_detail": "The ribB/ribD/ribE route to the isoalloxazine ring."},
    {"vitamin_id": "b7", "label": "Vitamin B7", "common_name": "biotin",
     "panel": "biotin",
     "route_detail": "The pimelate to biotin route through bioF, bioD and bioB."},
    {"vitamin_id": "b9", "label": "Vitamin B9", "common_name": "folate",
     "panel": "folate",
     "route_detail": "The pterin and pABA branches joined by folC."},
    {"vitamin_id": "b12", "label": "Vitamin B12", "common_name": "cobalamin",
     "panel": "b12",
     "route_detail": (
         "Corrinoid genes support more than one corrinoid, and several support salvage "
         "rather than synthesis. A B12 gene is not always a B12 result."
     )},
    {"vitamin_id": "k2", "label": "Vitamin K2", "common_name": "menaquinone",
     "panel": "k2",
     "route_detail": (
         "Menaquinone genes do not identify a chain length. MK-4 and MK-10 are different "
         "molecules with different behaviour, and these genes do not tell them apart."
     )},
)


# --------------------------------------------------------------------------- #
# grading
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class VitaminFinding:
    """One vitamin's three columns, kept apart."""

    vitamin: Vitamin
    by_column: Mapping[str, str]
    routes: tuple[dict[str, Any], ...]
    assayed: bool = True

    @property
    def can_synthesise(self) -> bool:
        return self.by_column.get("synthesis") == "complete"

    @property
    def searched(self) -> bool:
        """Whether any modelled route for this vitamin was searched for."""
        modelled = [s for s in self.by_column.values() if s != "not_applicable"]
        return bool(modelled) and not all(s == "not_assayed" for s in modelled)

    @property
    def uptake_only(self) -> bool:
        """The case the columns exist to make visible."""
        return (
            self.by_column.get("synthesis") != "complete"
            and self.by_column.get("salvage") != "complete"
            and self.by_column.get("uptake") == "complete"
        )

    @property
    def headline(self) -> str:
        modelled = [
            state for state in self.by_column.values() if state != "not_applicable"
        ]
        if not self.assayed or (modelled and all(s == "not_assayed" for s in modelled)):
            return "the genes for this vitamin were not searched for in this run"
        if self.can_synthesise:
            return "organisms here carry a complete route to make it"
        if self.by_column.get("salvage") == "complete":
            return "organisms here can rebuild it from a partial form, but not make it outright"
        if self.uptake_only:
            return "organisms here can take it in, but none carries a route to make it"
        if "partial" in self.by_column.values():
            return "part of a route is here; no complete one is"
        return "no route to it was found"

    def to_json(self) -> dict[str, Any]:
        return {
            "method_id": METHOD_VITAMIN,
            "vitamin_id": self.vitamin.vitamin_id,
            "metric_id": self.vitamin.metric_id,
            "label": self.vitamin.label,
            "common_name": self.vitamin.common_name,
            "columns": {
                column: {
                    "state": self.by_column.get(column, "absent"),
                    "meaning": COLUMN_LABELS[column],
                }
                for column, _ in COLUMNS
            },
            "columns_note": (
                "Three separate questions. They are never added together: taking a "
                "vitamin in is not making it, and an organism that only imports it "
                "depends on the vitamin already being there."
            ),
            "headline": self.headline,
            "can_synthesise": self.can_synthesise,
            "uptake_only": self.uptake_only,
            "routes": list(self.routes),
            "must_not_conclude": self.vitamin.must_not_conclude,
            "limitations": [
                "What the organisms present could make, which is not a serum level, not "
                "how much you absorb, and not a reason to stop a prescribed supplement.",
                self.vitamin.must_not_conclude,
            ],
        }


def assess_vitamin(
    vitamin: Vitamin, genes_present: Iterable[str], *, assayed: bool = True,
    genes_searched: Iterable[str] | None = None,
) -> VitaminFinding:
    """Grade one vitamin, column by column.

    A column is `complete` when any one of its routes is complete: two
    alternative synthesis routes are alternatives, and a community needs
    only one of them.
    """
    present = frozenset(str(g) for g in genes_present)
    searched = frozenset(str(g) for g in genes_searched) if genes_searched is not None else None
    by_column: dict[str, str] = {}
    payloads: list[dict[str, Any]] = []
    for column, _ in COLUMNS:
        states = []
        for route in vitamin.routes_in(column):
            state = route.state(present, assayed=assayed, searched=searched)
            states.append(state)
            payloads.append(route.to_json(present if assayed else None) | {"state": state})
        if not states:
            # This dashboard models no such route for this vitamin, which is
            # not the same as the route being absent from the sample.
            by_column[column] = "not_applicable"
        elif not assayed or all(x == "not_assayed" for x in states):
            by_column[column] = "not_assayed"
        elif "complete" in states:
            by_column[column] = "complete"
        elif "partial" in states:
            by_column[column] = "partial"
        else:
            by_column[column] = "absent"
    return VitaminFinding(vitamin, by_column, tuple(payloads), assayed=assayed)


def dashboard(
    genes_present: Iterable[str],
    *,
    existing_panels: Mapping[str, Any] | None = None,
    assayed: bool = True,
    genes_searched: Iterable[str] | None = None,
) -> dict[str, Any]:
    """All nine vitamins: four reconstructed here, five reused unchanged."""
    genes = list(genes_present)
    searched = list(genes_searched) if genes_searched is not None else None
    findings = [
        assess_vitamin(v, genes, assayed=assayed, genes_searched=searched)
        for v in NEW_VITAMINS
    ]
    panels = dict(existing_panels or {})
    reused = []
    for record in REUSED:
        panel = panels.get(record["panel"])
        reused.append({
            **record,
            "reused_unchanged": True,
            "existing_value": (panel or {}).get("copies_per_100_genomes"),
            "existing_confidence": (panel or {}).get("confidence"),
            "note": (
                "Reused from the existing panel without rescoring. The route detail "
                "explains what its genes do and do not establish."
            ),
        })
    return {
        "feature_id": "A03",
        "n_vitamins": len(findings) + len(reused),
        "assayed": assayed,
        "new": [f.to_json() for f in findings],
        "reused": reused,
        "columns": [{"column": c, "meaning": m} for c, m in COLUMNS],
        "limitations": [
            "Microbial capability, not serum vitamin status and not host absorption.",
            "A transport gene is never counted as synthesis.",
            "Never a reason to stop a prescribed supplement.",
        ],
        "provenance": (
            "Route definitions follow the curated reconstructions of Magnusdottir 2015 "
            "and Rodionov 2019 [F06-F07]. Those supply the routes and their test "
            "genomes; they are not evidence about this sample."
        ),
    }


def imported_gene_inventory() -> list[dict[str, Any]]:
    """Every gene this module uses, with the route and branch it belongs to.

    §5.3 asks for every imported gene, alternative and source table to be
    recorded; this is that record, generated from the definitions rather
    than maintained beside them.
    """
    out: list[dict[str, Any]] = []
    for vitamin in NEW_VITAMINS:
        for route in vitamin.routes:
            for branch in route.branches:
                for gene in branch.genes:
                    out.append({
                        "gene": gene, "vitamin_id": vitamin.vitamin_id,
                        "route_id": route.route_id, "branch_id": branch.branch_id,
                        "column": route.column, "role": "branch",
                    })
            for gene in route.coupling:
                out.append({
                    "gene": gene, "vitamin_id": vitamin.vitamin_id,
                    "route_id": route.route_id, "branch_id": None,
                    "column": route.column, "role": "coupling",
                })
    return out


__all__ = [
    "COLUMNS",
    "COLUMN_LABELS",
    "METHOD_VITAMIN",
    "NEW_VITAMINS",
    "REUSED",
    "ROUTE_STATES",
    "Branch",
    "Vitamin",
    "VitaminError",
    "VitaminFinding",
    "VitaminRoute",
    "assess_vitamin",
    "dashboard",
    "imported_gene_inventory",
]
