"""Fermentation and cross-feeding — A02, BUILD_SPEC_v0.8.3 §5.2.

The substrate → intermediate → terminal-product network, as a graph of
*supported biochemical opportunities*. What it is not, and the code
enforces each of these:

**Not a flow rate.** An edge says a reaction has support in this
community, not that anything is flowing along it. Nothing here is a rate,
a concentration or an amount.

**Not an assumption that the substrate is present.** Lactate utilisation
is capacity to consume lactate if there is lactate. It is not a
subtraction from lactate production, and the two are reported separately
because they are different measurements of different genes.

**Not a polarity where none exists.** Hydrogen production has no
universal adverse direction; a healthy community produces a great deal of
it. Production and consumption are shown apart for that reason, and the
direction word on any node is `descriptive` unless a cited context says
otherwise.

Two discriminations the specification calls out by name are enforced in
the route definitions: an AMP-forming `acs` supports acetate
*assimilation* and does not by itself support formation, and a glycyl
radical enzyme is not `IslA` without its activase.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

METHOD_ROUTE: Final = "ext083.fermentation_route/1.0"
METHOD_NETWORK: Final = "ext083.cross_feeding_network/1.0"

#: How much support an edge has. Only the first draws a solid line.
EDGE_SUPPORT: Final[tuple[str, ...]] = ("supported", "hypothesis", "unsupported", "not_assayed")

#: Node roles in the substrate → intermediate → product picture.
NODE_ROLES: Final[tuple[str, ...]] = ("substrate", "intermediate", "terminal_product", "sink")


class FermentationError(ValueError):
    """A route that would claim more than its genes support."""


@dataclass(frozen=True)
class Route:
    """One reconstructable route, and what its genes do and do not establish."""

    route_id: str
    label: str
    #: Genes without which the route cannot run. All of them are needed.
    required_genes: tuple[str, ...]
    #: Genes that would make the route unambiguous if present, but whose
    #: absence does not rule it out.
    corroborating_genes: tuple[str, ...] = ()
    #: Genes that are *not* evidence for this route however suggestive.
    #: Named so that a reviewer can see the trap that was avoided.
    insufficient_alone: Mapping[str, str] = field(default_factory=dict)
    produces: tuple[str, ...] = ()
    consumes: tuple[str, ...] = ()
    direction: str = "descriptive"
    note: str = ""

    def __post_init__(self) -> None:
        if not self.required_genes:
            raise FermentationError(f"{self.route_id}: a route needs required genes")
        if not (self.produces or self.consumes):
            raise FermentationError(
                f"{self.route_id}: a route must say what it produces or consumes"
            )
        overlap = set(self.required_genes) & set(self.insufficient_alone)
        if overlap:
            raise FermentationError(
                f"{self.route_id}: {sorted(overlap)} cannot be both required and "
                "insufficient-alone"
            )

    def to_json(self) -> dict[str, Any]:
        return {
            "route_id": self.route_id, "label": self.label,
            "required_genes": list(self.required_genes),
            "corroborating_genes": list(self.corroborating_genes),
            "insufficient_alone": dict(self.insufficient_alone),
            "produces": list(self.produces), "consumes": list(self.consumes),
            "direction": self.direction, "note": self.note,
        }


#: The seven measurements of §5.2, as routes. Gene symbols are the
#: conventional ones; the reference sequences behind them are installed
#: separately and each route reports `not_assayed` until they are.
ROUTES: Final[tuple[Route, ...]] = (
    Route(
        route_id="acetate.pta_ack", label="Acetate production",
        required_genes=("pta",), corroborating_genes=("ackA",), produces=("acetate",),
        insufficient_alone={
            "acs": (
                "AMP-forming acetyl-CoA synthetase runs towards acetate *assimilation*. "
                "Finding it is evidence that acetate can be consumed, not made."
            ),
        },
        note=(
            "Phosphotransacetylase is the committed step. Acetate kinase completes the "
            "route and is measured by the butyrate panel as a decoy, so it corroborates "
            "here rather than being required: duplicating it as a target would take "
            "reads from that panel rather than add a reading."
        ),
    ),
    Route(
        route_id="acetate.wood_ljungdahl", label="Acetate made from hydrogen and CO2",
        required_genes=("acsB", "cooS"), corroborating_genes=("fhs", "metF", "acsE"),
        produces=("acetate",), consumes=("hydrogen", "co2"),
        insufficient_alone={
            "fhs": (
                "Formate-tetrahydrofolate ligase is shared with other one-carbon "
                "metabolism and does not establish the acetogenic pathway."
            ),
        },
        note=(
            "Reconstructed separately from Pta-AckA because it is a different chemistry "
            "with a different substrate: it consumes hydrogen rather than sugar."
        ),
    ),
    Route(
        route_id="lactate.d_formation", label="D-lactate production",
        required_genes=("ldhA",), produces=("d_lactate",),
        note="Stereoisomer reported where the enzyme resolves it.",
    ),
    Route(
        route_id="lactate.l_formation", label="L-lactate production",
        required_genes=("ldh",), produces=("l_lactate",),
        note="Stereoisomer reported where the enzyme resolves it.",
    ),
    Route(
        route_id="lactate.utilisation_butyrate", label="Turning lactate into butyrate",
        required_genes=("lctA", "lctB", "lctC", "but"), corroborating_genes=("bcd", "etfAB"),
        consumes=("d_lactate", "l_lactate"), produces=("butyrate",),
        note=(
            "A cross-feeding sink for lactate, not a subtraction from the lactate a "
            "community can make. Both are capacities and neither measures an amount."
        ),
    ),
    Route(
        route_id="lactate.utilisation_propionate", label="Turning lactate into propionate",
        required_genes=("lcdA", "lcdB"), consumes=("d_lactate", "l_lactate"),
        produces=("propionate",),
        insufficient_alone={
            "pct": (
                "Propionate CoA-transferase appears in more than one propionate route "
                "and does not identify the acrylate pathway on its own."
            ),
        },
        note="The acrylate route, distinct from the succinate and propanediol routes.",
    ),
    Route(
        route_id="succinate.formation", label="Succinate production",
        required_genes=("mdh",), corroborating_genes=("frdA", "fumB", "pckA"),
        produces=("succinate",),
        note="A fermentative intermediate that other organisms consume.",
    ),
    Route(
        route_id="succinate.to_propionate", label="Turning succinate into propionate",
        required_genes=("mmdA", "scpA"), corroborating_genes=("mutA", "mutB"),
        consumes=("succinate",), produces=("propionate",),
        note=(
            "The succinate route to propionate, which is why succinate rarely "
            "accumulates in a community that carries it."
        ),
    ),
    Route(
        route_id="propionate.propanediol", label="Turning propanediol into propionate",
        required_genes=("pduC", "pduP"), corroborating_genes=("pduL", "pduW"),
        consumes=("fucose", "rhamnose"), produces=("propionate",),
        note="The deoxy-sugar route; a different chemistry from succinate or acrylate.",
    ),
    Route(
        route_id="hydrogen.production_fefe", label="Hydrogen production",
        required_genes=("hydA",), corroborating_genes=("hydE", "hydF", "hydG"),
        produces=("hydrogen",), direction="descriptive",
        note=(
            "Direction-resolved where the HydDB class allows it. A healthy community "
            "produces a great deal of hydrogen; there is no adverse direction here."
        ),
    ),
    Route(
        route_id="hydrogen.consumption_methanogenesis",
        label="Hydrogen used to make methane",
        required_genes=("mcrA",), consumes=("hydrogen",), produces=("methane",),
        note="Links to the existing methane reading without changing it.",
    ),
    Route(
        route_id="hydrogen.consumption_sulfate",
        label="Hydrogen used to make sulfide",
        required_genes=("dsrA",), corroborating_genes=("dsrB", "aprA", "sat"),
        consumes=("hydrogen", "sulfate"), produces=("hydrogen_sulfide",),
        insufficient_alone={
            "asrA": (
                "Anaerobic sulfite reductase is a different endpoint from dissimilatory "
                "sulfate reduction and is not interchangeable with it."
            ),
        },
        note=(
            "One of several sulfur routes, kept apart from the taurine, isethionate, "
            "sulfoacetate and cysteine routes because their endpoints differ. dsrA is "
            "measured by the hydrogen-sulfide panel, which owns it; this route reads "
            "that measurement rather than duplicating the target."
        ),
    ),
    Route(
        route_id="sulfur.taurine", label="Turning taurine into sulfide",
        required_genes=("tpa", "isl"), corroborating_genes=("xsc",),
        consumes=("taurine",), produces=("hydrogen_sulfide",),
        insufficient_alone={
            "grdB": (
                "A glycyl-radical enzyme is not IslA. Without its activase and the "
                "reaction context, the family alone does not establish this route."
            ),
        },
        note=(
            "Bilophila abundance is not this pathway, and this pathway is not every "
            "sulfur pathway."
        ),
    ),
)
BY_ID: Final[Mapping[str, Route]] = {r.route_id: r for r in ROUTES}


# --------------------------------------------------------------------------- #
# grading a route
# --------------------------------------------------------------------------- #

ROUTE_STATES: Final[tuple[str, ...]] = (
    "supported",        # every required gene present
    "partial",          # some required genes present, not all
    "insufficient",     # only genes that do not establish the route
    "absent",           # searched, nothing found
    "not_assayed",      # not searched in this run
)


@dataclass(frozen=True)
class RouteFinding:
    route: Route
    state: str
    present: tuple[str, ...]
    missing: tuple[str, ...]
    corroborating_present: tuple[str, ...]
    insufficient_present: tuple[str, ...]

    @property
    def supported(self) -> bool:
        return self.state == "supported"

    def to_json(self) -> dict[str, Any]:
        payload = {
            "method_id": METHOD_ROUTE, **self.route.to_json(),
            "state": self.state,
            "required_present": list(self.present),
            "required_missing": list(self.missing),
            "corroborating_present": list(self.corroborating_present),
            "limitations": [
                "Genetic capacity for a reaction, not a rate and not an amount.",
            ],
        }
        if self.insufficient_present:
            payload["insufficient_alone_present"] = {
                gene: self.route.insufficient_alone[gene]
                for gene in self.insufficient_present
            }
            payload["limitations"].append(
                "Genes were found that are suggestive but do not establish this route on "
                "their own; each says why."
            )
        return payload


def assess_route(
    route: Route, genes_present: Iterable[str], *, assayed: bool = True
) -> RouteFinding:
    """Grade one route. Partial is its own state, not a quiet negative."""
    present_set = {str(g) for g in genes_present}
    present = tuple(g for g in route.required_genes if g in present_set)
    missing = tuple(g for g in route.required_genes if g not in present_set)
    corroborating = tuple(g for g in route.corroborating_genes if g in present_set)
    insufficient = tuple(g for g in route.insufficient_alone if g in present_set)

    if not assayed:
        state = "not_assayed"
    elif not missing:
        state = "supported"
    elif present:
        state = "partial"
    elif insufficient:
        state = "insufficient"
    else:
        state = "absent"
    return RouteFinding(route, state, present, missing, corroborating, insufficient)


def assess_all(
    genes_present: Iterable[str], *, assayed: bool = True
) -> list[RouteFinding]:
    genes = list(genes_present)
    return [assess_route(route, genes, assayed=assayed) for route in ROUTES]


# --------------------------------------------------------------------------- #
# the network
# --------------------------------------------------------------------------- #

#: Where each metabolite sits in the picture.
NODE_ROLE: Final[Mapping[str, str]] = {
    "fucose": "substrate", "rhamnose": "substrate", "taurine": "substrate",
    "sulfate": "substrate", "co2": "substrate",
    "d_lactate": "intermediate", "l_lactate": "intermediate",
    "succinate": "intermediate", "hydrogen": "intermediate",
    "acetate": "terminal_product", "butyrate": "terminal_product",
    "propionate": "terminal_product", "methane": "terminal_product",
    "hydrogen_sulfide": "terminal_product",
}

NODE_LABELS: Final[Mapping[str, str]] = {
    "d_lactate": "D-lactate", "l_lactate": "L-lactate",
    "hydrogen_sulfide": "hydrogen sulfide", "co2": "carbon dioxide",
}


def node_label(node: str) -> str:
    return NODE_LABELS.get(node, node.replace("_", " "))


@dataclass(frozen=True)
class Edge:
    """One arrow: a reaction this community has support for."""

    source: str
    target: str
    route_id: str
    support: str

    @property
    def style(self) -> str:
        """Solid only where the reaction evidence is there."""
        return "solid" if self.support == "supported" else "dashed"

    def to_json(self) -> dict[str, Any]:
        return {
            "from": self.source, "from_label": node_label(self.source),
            "to": self.target, "to_label": node_label(self.target),
            "route_id": self.route_id, "support": self.support, "style": self.style,
        }


def network(findings: Sequence[RouteFinding]) -> dict[str, Any]:
    """The cross-feeding graph: supported opportunities, not flows."""
    edges: list[Edge] = []
    for finding in findings:
        support = {
            "supported": "supported", "partial": "hypothesis",
            "insufficient": "hypothesis", "absent": "unsupported",
            "not_assayed": "not_assayed",
        }[finding.state]
        if support == "not_assayed":
            continue
        sources = finding.route.consumes or ("substrate",)
        for source in sources:
            for target in finding.route.produces:
                edges.append(Edge(source, target, finding.route.route_id, support))
    nodes = sorted({e.source for e in edges} | {e.target for e in edges})
    return {
        "method_id": METHOD_NETWORK,
        "nodes": [
            {"id": n, "label": node_label(n), "role": NODE_ROLE.get(n, "intermediate")}
            for n in nodes
        ],
        "edges": [e.to_json() for e in edges if e.support != "unsupported"],
        "unsupported_edges": [e.to_json() for e in edges if e.support == "unsupported"],
        "legend": {
            "solid": "a reaction this community has the genes for",
            "dashed": "compatible with the literature, not established by these genes",
        },
        "limitations": [
            "Supported biochemical opportunities, not flows. No edge carries a rate.",
            "An edge does not assume its substrate is present: a lactate sink is "
            "capacity to use lactate, not evidence that there is lactate to use.",
            "Production and consumption of the same metabolite are separate "
            "measurements; neither is subtracted from the other.",
            "Filling a node does not mean an added organism would establish there.",
        ],
    }


def production_and_consumption(findings: Sequence[RouteFinding]) -> dict[str, Any]:
    """Each metabolite's routes in, and its routes out, side by side.

    Never netted. §5.2 is explicit that a consumption route is a
    complementary sink and not a measured subtraction, and hydrogen is the
    case that makes the reason obvious: a community that both makes and
    uses a great deal of hydrogen is not a community with no hydrogen.
    """
    out: dict[str, dict[str, list[str]]] = {}
    for finding in findings:
        if finding.state in {"absent", "not_assayed"}:
            continue
        for node in finding.route.produces:
            out.setdefault(node, {"produced_by": [], "consumed_by": []})
            out[node]["produced_by"].append(finding.route.route_id)
        for node in finding.route.consumes:
            out.setdefault(node, {"produced_by": [], "consumed_by": []})
            out[node]["consumed_by"].append(finding.route.route_id)
    return {
        node: {
            "label": node_label(node),
            "role": NODE_ROLE.get(node, "intermediate"),
            "produced_by": sorted(routes["produced_by"]),
            "consumed_by": sorted(routes["consumed_by"]),
            "note": (
                "Routes in and routes out, side by side. They are not netted against "
                "each other: both are capacities, and neither is an amount."
            ),
        }
        for node, routes in sorted(out.items())
    }


def summarise(
    genes_present: Iterable[str], *, assayed: bool = True
) -> dict[str, Any]:
    """The whole A02 view."""
    findings = assess_all(genes_present, assayed=assayed)
    return {
        "feature_id": "A02",
        "assayed": assayed,
        "n_routes": len(findings),
        "n_supported": sum(1 for f in findings if f.supported),
        "routes": [f.to_json() for f in findings],
        "network": network(findings),
        "metabolites": production_and_consumption(findings),
        "preserved_unchanged": [
            "butyrate", "propionate", "methane", "hydrogen_sulfide", "tma",
        ],
        "preserved_note": (
            "The existing butyrate, propionate, methane, sulfide and TMA readings are "
            "unchanged. This adds a breakdown of the routes behind them; it does not "
            "rescore them."
        ),
        "limitations": [
            "Hydrogen production has no universal adverse direction, which is why "
            "production and consumption are shown separately rather than netted.",
            "Bilophila abundance is not the taurine pathway, and the taurine pathway is "
            "not every sulfur pathway.",
        ],
    }


__all__ = [
    "BY_ID",
    "EDGE_SUPPORT",
    "METHOD_NETWORK",
    "METHOD_ROUTE",
    "NODE_LABELS",
    "NODE_ROLE",
    "NODE_ROLES",
    "ROUTES",
    "ROUTE_STATES",
    "Edge",
    "FermentationError",
    "Route",
    "RouteFinding",
    "assess_all",
    "assess_route",
    "network",
    "node_label",
    "production_and_consumption",
    "summarise",
]
