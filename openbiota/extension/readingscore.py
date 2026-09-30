"""A percentile for every functional reading — so each one gets a slider.

The substrates, routes, modules, steps and capacities arrived in the report
as words: "class-level only", "part of it", "supported". Beside forty-two
panel rows that each carry a percentile and a bar, a column of adjectives
reads as missing data, and a reader cannot tell high from low in it at all.

They had no percentile because the cohort stores reference distributions
per panel and per *panel entry*, and none of these readings is a panel
entry. What each one is, though, is a requirement over entries: cellulose
needs the GH5 and GH9 families, the acrylate route needs lcdA, lcdB and
lcdC. So the distributions to compare against already exist - 162 of them,
each over the same 91 samples - and what was missing was the join.

A composite of them, however, is not something this module may invent, and
an earlier version of it did.

It combined a reading's genes with a limiting-step rule - a reading is no
more present than its scarcest gene - on the reasoning that a route needs
every one of its parts. That reasoning is wrong for this data. Every panel
in the reference set declares ``aggregate: sum`` and none declares a
chain, because the gene sets here are alternatives rather than sequences:
*gshF* is the bifunctional enzyme that replaces *gshA* plus *gshB*,
putrescine has three independent routes to it, and cellulose is degraded
by either the GH5 or the GH9 family. Taking a minimum over alternatives
answers no question at all.

What it did instead was contradict the panels. The two agreed exactly
whenever a reading had one gene and disagreed whenever it had more: the
glutathione panel read 81st on the sum of its routes while this module
read zero, having found one absent alternative and called the route dead.
Both numbers were on the same page under the same name.

So the composite is gone. A reading carries a position when, and only
when, one exists to carry: a single-gene reading *is* its gene and takes
that gene's place in the cohort, which is a real distribution over real
samples. A multi-gene reading has no distribution of its own - the cohort
stores panels and entries, not arbitrary subsets - and reports its genes
with their own positions rather than a number nobody measured.

The scored panel remains the one headline for its quantity. These are its
components, and they are shown inside it.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

METHOD_READING_SCORE: Final = "ext083.reading_percentile/1.0"

#: How a reading's value is formed from its genes.
#:
#: The sum, which is what every panel in this reference declares and what
#: these gene sets are: alternatives and parallel routes, not the steps of
#: a chain. Summing them measures the capacity present; the minimum this
#: module used to take measured which alternative was rarest, which is not
#: a question anybody asked and which contradicted the panels.
RULE: Final = "sum"

#: Why there is a position, in one sentence the report can print.
RULE_NOTE: Final = (
    "This reading is the sum of its genes, placed against the same sum measured across "
    "the reference cohort. Each contributing gene is listed with its own position as "
    "well, so the total can be taken apart."
)

#: Said when the cohort has no distribution for this reading's gene set.
NO_COMPOSITE_NOTE: Final = (
    "No reference distribution has been built for this reading's combination of genes, "
    "so it is reported by its parts rather than with a position nothing measured."
)


class ReadingScoreError(ValueError):
    """A score that would be derived from nothing."""


@dataclass(frozen=True)
class GenePosition:
    """One contributing gene, and where it sits in the reference."""

    key: str
    gene: str
    #: What the enzyme does, from the panel entry's own label. A reader
    #: cannot act on "pta"; they can act on "phosphate acetyltransferase".
    label: str
    #: The organisms this gene's reads were attributed to in this sample,
    #: largest first. A reading is caused by organisms, and a card that
    #: lists only genes cannot tell anybody what to act on.
    organisms: tuple[dict[str, Any], ...]
    value: float | None
    percentile: float | None
    detected_in: int | None
    cohort_n: int | None

    def to_json(self) -> dict[str, Any]:
        return {
            "key": self.key, "gene": self.gene, "label": self.label,
            "organisms": list(self.organisms),
            "value": self.value,
            "percentile": self.percentile, "detected_in": self.detected_in,
            "cohort_n": self.cohort_n,
        }


@dataclass(frozen=True)
class ReadingScore:
    """A functional reading's position against the reference cohort."""

    reading_id: str
    percentile: float | None
    value: float | None
    rule: str
    contributors: tuple[GenePosition, ...]
    limiting: str | None
    unscored_reason: str | None = None

    @property
    def scored(self) -> bool:
        return self.percentile is not None

    def to_json(self) -> dict[str, Any]:
        return {
            "method": METHOD_READING_SCORE,
            "reading_id": self.reading_id,
            "percentile": self.percentile,
            "value": self.value,
            "rule": self.rule,
            "rule_note": RULE_NOTE if self.scored else None,
            "limiting_gene": self.limiting,
            "contributors": [c.to_json() for c in self.contributors],
            "unscored_reason": self.unscored_reason,
        }


def percentile_of(percentiles: Mapping[str, float], value: float) -> float:
    """Where a value falls in a stored percentile ladder, 0-100.

    A measured zero is the bottom of the scale and returns zero, for the
    same reason the panel readings do: a sample with none of something is
    not above a quarter of the reference merely because the reference is
    mostly zeros too.
    """
    if not value:
        return 0.0
    points = sorted((int(p), v) for p, v in percentiles.items())
    if value <= points[0][1]:
        return float(points[0][0])
    if value >= points[-1][1]:
        return float(points[-1][0])
    for (p_lo, v_lo), (p_hi, v_hi) in zip(points, points[1:], strict=False):
        if v_lo <= value <= v_hi:
            if v_hi == v_lo:
                return float(p_hi)
            span = (value - v_lo) / (v_hi - v_lo)
            return float(p_lo) + span * (p_hi - p_lo)
    return float(points[-1][0])


def _organisms_of(entry: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    """The organisms a gene's reads came closest to, largest share first.

    This is what makes a reading actionable: the genes say what chemistry
    is present and the organisms say who is carrying it. Both are recorded
    so a card can show the second rather than only the first.
    """
    rows = [
        {
            "organism": str(o.get("organism") or ""),
            "fragments": int(o.get("fragments") or 0),
        }
        for o in (entry.get("organisms") or ())
        if o.get("organism")
    ]
    total = sum(r["fragments"] for r in rows) or 1
    for row in rows:
        row["share"] = round(row["fragments"] / total, 4)
    rows.sort(key=lambda r: (-r["fragments"], r["organism"]))
    return tuple(rows[:8])


def _entry_index(results: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Every panel entry in this sample, keyed `panel:ENTRY`."""
    out: dict[str, dict[str, Any]] = {}
    for panel in results.get("panels") or []:
        if not isinstance(panel, dict):
            continue
        name = str(panel.get("name"))
        for kind in ("genes", "decoys"):
            for entry in panel.get(kind) or []:
                if isinstance(entry, dict) and entry.get("entry_id"):
                    out[f"{name}:{entry['entry_id']}"] = entry
    return out


def _by_gene(results: Mapping[str, Any]) -> dict[str, list[str]]:
    """Gene symbol to the entry keys that measure it."""
    out: dict[str, list[str]] = {}
    for key, entry in _entry_index(results).items():
        gene = str(entry.get("gene") or "").strip()
        if gene:
            out.setdefault(gene, []).append(key)
    return out


def score(
    reading_id: str,
    genes: Sequence[str],
    *,
    results: Mapping[str, Any],
    ranges: Mapping[str, Any],
    searched: Iterable[str] | None = None,
) -> ReadingScore:
    """Place one reading against the cohort, by its limiting gene."""
    gene_ranges = ranges.get("genes") or {}
    entries = _entry_index(results)
    by_gene = _by_gene(results)
    looked_for = set(searched) if searched is not None else None

    if not genes:
        return ReadingScore(
            reading_id, None, None, RULE, (), None,
            "this reading has no required genes declared, so there is nothing to place",
        )
    if looked_for is not None and not (set(genes) & looked_for):
        return ReadingScore(
            reading_id, None, None, RULE, (), None,
            "none of this reading's genes was searched for in this run, which is "
            "different from searching and finding nothing",
        )

    contributors: list[GenePosition] = []
    for gene in genes:
        for key in by_gene.get(gene, []):
            spec = gene_ranges.get(key)
            entry = entries.get(key) or {}
            value = entry.get("copies_per_100_genomes")
            pct = (
                percentile_of(spec["percentiles"], float(value))
                if spec and isinstance(value, int | float) else None
            )
            contributors.append(GenePosition(
                key=key, gene=gene,
                label=str(entry.get("label") or gene),
                organisms=_organisms_of(entry),
                value=float(value) if isinstance(value, int | float) else None,
                percentile=None if pct is None else round(pct, 1),
                detected_in=(spec or {}).get("detected_in"),
                cohort_n=(spec or {}).get("cohort_n"),
            ))

    measured = [c for c in contributors if c.value is not None]
    if not measured:
        return ReadingScore(
            reading_id, None, None, RULE, tuple(contributors), None,
            "none of this reading's genes was measured in this sample, so it is "
            "reported without a position rather than with an invented one",
        )

    # The reading's own value: the sum of its genes, as its panel would do.
    total = float(sum(c.value or 0.0 for c in measured))

    # And its own distribution, built over the same cohort by summing the
    # same genes for every sample in it. That is what makes this a position
    # rather than an assertion.
    spec = (ranges.get("readings") or {}).get(reading_id)
    ladder = (spec or {}).get("percentiles")
    if not ladder and len(measured) == 1:
        # A sum over one gene is that gene, so its distribution is this
        # reading's distribution exactly rather than as an approximation.
        ladder = (gene_ranges.get(measured[0].key) or {}).get("percentiles")
    if not ladder:
        return ReadingScore(
            reading_id, None, total, RULE, tuple(contributors), None, NO_COMPOSITE_NOTE,
        )
    position = percentile_of(ladder, total)
    # The gene that contributes most, named so a reader can see what is
    # carrying the total rather than being told only the total.
    largest = max(measured, key=lambda c: (c.value or 0.0, c.key))
    return ReadingScore(
        reading_id=reading_id,
        percentile=round(position, 1),
        value=round(total, 4),
        rule=RULE,
        contributors=tuple(contributors),
        limiting=largest.gene,
    )


def score_all(
    readings: Sequence[Mapping[str, Any]],
    *,
    results: Mapping[str, Any],
    ranges: Mapping[str, Any],
    searched: Iterable[str] | None = None,
) -> dict[str, ReadingScore]:
    """Score every reading that declares its genes.

    ``readings`` are ``{"reading_id": ..., "genes": [...]}`` mappings, which
    is what each module already knows about its own rows.
    """
    looked = list(searched) if searched is not None else None
    out: dict[str, ReadingScore] = {}
    for reading in readings:
        rid = str(reading.get("reading_id") or "")
        if not rid:
            continue
        out[rid] = score(
            rid, [str(g) for g in reading.get("genes") or ()],
            results=results, ranges=ranges, searched=looked,
        )
    return out
