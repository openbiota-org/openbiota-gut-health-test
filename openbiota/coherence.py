"""Do two measurements of the same biology agree, and if not, why not?

The pipeline measures some functions twice, by two methods that do not share
evidence:

* a **gene panel** — translated search of every read against curated protein
  reference sets, normalised to ``rpoB``. It counts a pathway's genes wherever
  they sit in the community, in any organism, named or not.
* a **taxon group** — the summed relative abundance of the species that are
  characterised as doing that job. It counts organisms, not genes, and only
  the organisms someone has characterised.

These answer different questions, so they can legitimately differ, and when
they do the difference is information rather than noise. What is *not*
acceptable is for the report to print "more butyrate capacity than most
people" beside "butyrate producers notably low" with nothing in between: a
reader cannot reconcile that, and it destroys trust in both numbers.

This module pairs every panel that has a curated carrier group, grades the
agreement, and produces one sentence that names the mechanism. The report
prints that sentence under *both* readings, so whichever one the reader meets
first, they meet the explanation with it.

Three mechanisms account for essentially every real divergence, and the pair
table records which ones apply to each function:

``carriers_outside_set``
    The genes are real but sit in organisms outside the curated group — either
    uncharacterised species or ones nobody has assigned to the guild. Gene
    capacity is then the broader measurement and the group is the narrower.
``homolog_inflation``
    The gene family is promiscuous, so a reference set built from protein
    names admits organisms that carry the fold without running the pathway.
    The group is then the more trustworthy of the two. Scope rules in the
    panel YAML suppress the known cases (see :class:`openbiota.panels.Scope`); this
    mechanism stays listed because suppression is never complete.
``below_detection``
    The characterised species are present but under the profiler's limit of
    detection, so the group reads low while their genes are still counted.
``incomplete_pathway``
    The carriers are present but the panel found fewer terminal-step genes
    than expected — a strain-level gene-variant miss, or organisms that are
    assigned to the guild on 16S grounds without the terminal gene.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Final

# --------------------------------------------------------------------------- #
# what pairs with what
# --------------------------------------------------------------------------- #

#: Percentile gap at which the two readings stop being a plausible pair. Below
#: this they are treated as agreeing to within the noise of two different
#: measurements of a skewed quantity; above it the report explains the split.
DIVERGENCE_GAP: Final = 30.0
#: Gap at which the readings land in opposite halves of the reference and the
#: split is the headline rather than a footnote.
CONTRADICTION_GAP: Final = 50.0


@dataclass(frozen=True, slots=True)
class Pair:
    """One function measured both ways."""

    panel: str
    group: str
    #: What the gene panel measures, in the report's own words.
    capacity_label: str
    #: What the taxon group measures.
    carrier_label: str
    #: Mechanisms that can separate the two, most likely first.
    mechanisms: tuple[str, ...]
    #: Which reading is the better guide when they disagree, and for what.
    #: One of ``capacity``, ``carriers`` or ``neither``.
    prefer: str
    #: Function-specific sentence, appended to the generic explanation.
    note: str = ""


#: Every panel with a curated carrier group. A panel absent from this table is
#: measured once and cannot contradict a group; a group absent from it (the
#: ecological ones — oral-origin, opportunistic pathogens, hexa-LPS) has no
#: pathway panel to pair with.
PAIRS: Final = (
    Pair(
        panel="butyrate", group="butyrate",
        capacity_label="butyrate production capacity",
        carrier_label="butyrate-producing species",
        mechanisms=("carriers_outside_set", "below_detection", "homolog_inflation"),
        prefer="carriers",
        note=(
            "The terminal butyrate gene sits in many Lachnospiraceae and Ruminococcaceae that no "
            "curated list names, so gene capacity routinely reads higher than the named producers. "
            "Where the two split, the named producers are the better guide to how much butyrate this "
            "community actually makes, because butyrate output tracks the abundance of the organisms "
            "that specialise in it. Neither number is a butyrate measurement: production also needs "
            "fermentable fibre reaching the colon."
        ),
    ),
    Pair(
        panel="h2s", group="sulfate_reducers",
        capacity_label="hydrogen-sulfide capacity",
        carrier_label="sulfate-reducing species",
        mechanisms=("carriers_outside_set", "below_detection"),
        prefer="capacity",
        note=(
            "Hydrogen sulfide has two independent sources: dissimilatory sulfate reduction by "
            "Desulfovibrio and relatives, and cysteine degradation, which is widespread in ordinary "
            "gut bacteria. The gene panel sees both; the group sees only the sulfate reducers, which "
            "are a minority of stool even when sulfide capacity is high."
        ),
    ),
    Pair(
        panel="methane", group="methanogens",
        capacity_label="methane capacity (mcrA)",
        carrier_label="methanogen abundance",
        mechanisms=("below_detection", "incomplete_pathway"),
        prefer="neither",
        note=(
            "mcrA is exclusive to methanogenic archaea, so these two should track each other tightly. "
            "A split means one of the two engines is at its limit: methanogens are archaea, they are "
            "a small fraction of stool DNA, and the bacterial profiler under-reports them."
        ),
    ),
    Pair(
        panel="bsh", group="probiotics",
        capacity_label="bile-salt deconjugation capacity",
        carrier_label="probiotic-genus species",
        mechanisms=("carriers_outside_set",),
        prefer="capacity",
        note=(
            "Bile salt hydrolase is not a probiotic trait: it is carried across Lactobacillaceae, "
            "Bifidobacterium, Clostridium, Bacteroides and Enterococcus alike. A high capacity with "
            "few probiotic-genus species simply means the ordinary residents are doing it."
        ),
    ),
    Pair(
        panel="urda", group="polyphenol_metabolisers",
        capacity_label="urolithin and equol capacity",
        carrier_label="polyphenol-metabolising species",
        mechanisms=("carriers_outside_set", "homolog_inflation"),
        prefer="carriers",
        note=(
            "Urolithin and equol production is strain-level: two isolates of the same species differ, "
            "and only a handful of characterised strains do it. The curated species are the better "
            "guide; a high gene reading without them is most likely a related enzyme."
        ),
    ),
)

PAIR_BY_PANEL: Final = {p.panel: p for p in PAIRS}
PAIR_BY_GROUP: Final = {p.group: p for p in PAIRS}

MECHANISM_TEXT: Final = {
    "carriers_outside_set": (
        "the genes are carried by organisms outside the curated list, which the gene search counts "
        "and the species list does not"
    ),
    "homolog_inflation": (
        "the gene family is shared with enzymes that do a different job, so some matched fragments "
        "are not this pathway"
    ),
    "below_detection": (
        "the characterised species are present below the profiler's limit of detection, so their "
        "genes are counted while the organisms are not"
    ),
    "incomplete_pathway": (
        "the organisms are present but fewer of the terminal-step genes were found than their "
        "genomes would predict"
    ),
}

VERDICT_AGREE: Final = "agree"
VERDICT_DIVERGENT: Final = "divergent"
VERDICT_OPPOSED: Final = "opposed"


# --------------------------------------------------------------------------- #
# the check
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class Coherence:
    """The result of measuring one function two ways."""

    pair: Pair
    panel_percentile: float
    group_percentile: float
    verdict: str
    #: Mechanisms that could produce a split this large, most likely first.
    mechanisms: tuple[str, ...] = ()
    #: Fraction of this panel's matched fragments that landed on references
    #: demoted to background by a scope rule, when the panel declares one.
    out_of_scope_fraction: float | None = None

    @property
    def gap(self) -> float:
        return abs(self.panel_percentile - self.group_percentile)

    @property
    def agrees(self) -> bool:
        return self.verdict == VERDICT_AGREE

    @property
    def capacity_is_higher(self) -> bool:
        return self.panel_percentile > self.group_percentile

    def headline(self) -> str:
        """One line naming the split, for a heading or a chip."""
        hi, lo = (
            (self.pair.capacity_label, self.pair.carrier_label)
            if self.capacity_is_higher
            else (self.pair.carrier_label, self.pair.capacity_label)
        )
        return f"Higher {hi} than {lo}"

    def explanation(self, side: str = "capacity") -> str:
        """The sentence printed under a reading.

        ``side`` is where the reader is standing: ``capacity`` under the gene
        panel, ``carriers`` under the taxon group. The facts, the mechanism and
        the verdict on which number to weight are identical from both sides;
        only the pointer to the other reading changes. That way the two halves
        of the report never say different things about the same split.
        """
        cap, car = _ord(self.panel_percentile), _ord(self.group_percentile)
        this_label, other_label, this_pct, other_pct = (
            (self.pair.capacity_label, self.pair.carrier_label, cap, car)
            if side == "capacity"
            else (self.pair.carrier_label, self.pair.capacity_label, car, cap)
        )
        how_other = (
            "counting the organisms characterised for the job instead of the genes"
            if side == "capacity"
            else "counting the pathway's genes across the whole community instead of these species"
        )
        if self.agrees:
            return (
                f"This agrees with the {other_label} reading ({other_pct} percentile), which measures "
                f"the same function by {how_other}."
            )
        why = "; ".join(MECHANISM_TEXT[m] for m in self.mechanisms[:2])
        prefer = {
            "capacity": (
                "Of the two, the gene reading is the better guide here, because this capacity does "
                "not depend on any one set of species."
            ),
            "carriers": "Of the two, the species reading is the better guide here.",
            "neither": (
                "Neither reading is preferred: a split this size means one of the two engines is at "
                "its limit for this function."
            ),
        }[self.pair.prefer]
        return (
            f"<b>Read this beside the {other_label} reading.</b> The {this_label} reading sits at the "
            f"{this_pct} percentile; the {other_label} reading is at the {other_pct}. The two are measured "
            "differently — genes in the whole community against the abundance of the species "
            f"characterised for the job — and they can separate when {why}. {prefer} {self.pair.note}"
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "panel": self.pair.panel,
            "group": self.pair.group,
            "panel_percentile": round(self.panel_percentile, 1),
            "group_percentile": round(self.group_percentile, 1),
            "gap": round(self.gap, 1),
            "verdict": self.verdict,
            "mechanisms": list(self.mechanisms),
            "prefer": self.pair.prefer,
            "out_of_scope_fraction": (
                None if self.out_of_scope_fraction is None else round(self.out_of_scope_fraction, 4)
            ),
            "explanation": _plain(self.explanation()),
        }


@dataclass(slots=True)
class CoherenceReport:
    """Every double-measured function, indexed for the report."""

    results: list[Coherence] = field(default_factory=list)

    def by_panel(self, panel: str) -> Coherence | None:
        for c in self.results:
            if c.pair.panel == panel:
                return c
        return None

    def by_group(self, group: str) -> Coherence | None:
        for c in self.results:
            if c.pair.group == group:
                return c
        return None

    @property
    def divergent(self) -> list[Coherence]:
        return [c for c in self.results if not c.agrees]

    def to_json(self) -> dict[str, Any]:
        return {
            "what_this_is": (
                "Functions this pipeline measures twice, by gene search and by carrier abundance. "
                "Agreement is corroboration; divergence is reported with the mechanism that "
                "explains it and is never averaged away."
            ),
            "n_double_measured": len(self.results),
            "n_divergent": len(self.divergent),
            "divergence_gap_percentile_points": DIVERGENCE_GAP,
            "results": [c.to_json() for c in self.results],
        }


def check_coherence(
    *,
    panel_percentiles: Mapping[str, float | None],
    group_percentiles: Mapping[str, float | None],
    out_of_scope_fractions: Mapping[str, float] | None = None,
) -> CoherenceReport:
    """Grade every double-measured function.

    Both inputs are percentiles against their own reference cohorts, which is
    what makes them comparable at all: the gene figure is copies per 100
    genomes against 34 adult stool metagenomes screened identically, the group
    figure is a summed relative abundance against 3,027 adult stool profiles.
    Neither raw scale means anything next to the other; their positions do.
    """
    oos = out_of_scope_fractions or {}
    report = CoherenceReport()
    for pair in PAIRS:
        pp = panel_percentiles.get(pair.panel)
        gp = group_percentiles.get(pair.group)
        if pp is None or gp is None:
            continue
        gap = abs(pp - gp)
        if gap < DIVERGENCE_GAP:
            verdict = VERDICT_AGREE
        elif gap < CONTRADICTION_GAP:
            verdict = VERDICT_DIVERGENT
        else:
            verdict = VERDICT_OPPOSED

        # Order the mechanisms by which direction the split runs. Capacity
        # above carriers is explained by genes outside the set or by homolog
        # inflation; carriers above capacity is explained by a pathway the
        # panel under-counts.
        if verdict == VERDICT_AGREE:
            mechanisms: tuple[str, ...] = ()
        elif pp > gp:
            mechanisms = tuple(
                m for m in pair.mechanisms
                if m in ("carriers_outside_set", "homolog_inflation", "below_detection")
            ) or pair.mechanisms
        else:
            mechanisms = tuple(
                m for m in pair.mechanisms if m in ("incomplete_pathway", "below_detection")
            ) or ("incomplete_pathway",)

        report.results.append(Coherence(
            pair=pair, panel_percentile=float(pp), group_percentile=float(gp),
            verdict=verdict, mechanisms=mechanisms,
            out_of_scope_fraction=oos.get(pair.panel),
        ))
    return report


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _ord(percentile: float) -> str:
    """``63rd``, with the extremes named as positions rather than as 0 or 100."""
    if percentile >= 99.5:
        return ">99th"
    if percentile < 1.0:
        return "<1st"
    n = int(round(percentile))
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def _plain(markup: str) -> str:
    """Strip the report's inline bold for the JSON copy."""
    return markup.replace("<b>", "").replace("</b>", "")
