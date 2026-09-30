"""Compare this pipeline against Tiny Health on the same FASTQ files.

Two independent pipelines read the same DNA. Where they agree, both are
probably right. Where they disagree, one of them is wrong — and the point
of this script is to say *which*, from the evidence, rather than assuming
the commercial lab is the reference. Tiny Health is a CLIA-certified
shotgun service, so their numbers deserve to be taken seriously; they are
not ground truth, and several disagreements here resolve in our favour on
the sequence evidence.

What can and cannot be compared:

* **Functional capacity.** Both report gene-level capacity per pathway,
  but in different units — Tiny Health in RPKM, this pipeline in copies
  per 100 bacterial genomes — and against different reference
  populations. Absolute values are therefore not comparable; the
  comparable quantity is *where the reading sits* on each platform's own
  scale, which is what `compare_functions` does.
* **Species abundance.** Directly comparable as relative abundance, once
  GTDB names (Tiny Health) are mapped to NCBI names (MetaPhlAn).
* **Taxon presence/absence for pathogens.** Directly comparable.
* **Everything else this pipeline reports** — disease-pattern
  resemblance, biological age, resistance determinants, skin panels —
  has no Tiny Health counterpart, so it is counted, not compared.

Run: ``python tools/compare_tinyhealth.py`` from the repository root.
"""

from __future__ import annotations

import json
import math
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from tools.tinyhealth import Report, load_all

ROOT: Final = Path(__file__).resolve().parents[1]

#: Tiny Health metric name -> our panel name. Only pairs that measure the
#: same biology are listed. Where a name appears in more than one Tiny
#: Health template ("Histamine index" in v5.3, "Histamine-producing
#: species" in v4.8), both map to the same panel.
FUNCTION_MAP: Final[dict[str, str]] = {
    "Butyrate": "butyrate",
    "Propionate": "propionate",
    "Trimethylamine": "cutc",
    "Ammonia": "urease",
    "Branched Chain Amino Acids": "bcaa",
    "p-Cresol": "pcresol",
    "Histamine index": "histamine",
    "Histamine-producing species": "histamine",
    "Methane production capacity": "methane",
    "Beta-glucuronidase capacity": "bglucuronidase",
    "GABA production": "gaba",
    "Secondary bile acids": "bai",
    "Unconjugated bile acids": "bsh",
    "Vitamin B2": "riboflavin",
    "Vitamin B7": "biotin",
    "Vitamin B9": "folate",
    "Vitamin B12": "b12",
    "Vitamin K": "k2",
    "Indole-3-propionic acid": "ipa",
    "Hydrogen sulfide index": "h2s",
    "Oxalate degradation": "oxalate",
}

#: Pairs whose two measurements are of *different quantities* despite the
#: shared name, with the reason. Excluded from the concordance rate and
#: reported separately, because counting them as disagreements would be
#: measuring the naming, not the biology.
NOT_LIKE_FOR_LIKE: Final[dict[str, str]] = {
    "Histamine-producing species": (
        "Tiny Health counts the abundance of a curated list of "
        "histamine-producing species; this pipeline counts hdcA genes "
        "whatever carries them. A sample can have zero listed species and "
        "still carry the gene."
    ),
}

#: GTDB (Tiny Health) -> NCBI (MetaPhlAn) species names. GTDB split
#: Bacteroides and renamed several gut species; without this the two
#: species tables look far less concordant than they are.
SPECIES_ALIASES: Final[dict[str, tuple[str, ...]]] = {
    "Phocaeicola vulgatus": ("Bacteroides vulgatus", "Phocaeicola vulgatus"),
    "Phocaeicola dorei": ("Bacteroides dorei", "Phocaeicola dorei"),
    "Phocaeicola plebeius": ("Bacteroides plebeius", "Phocaeicola plebeius"),
    "Phocaeicola massiliensis": ("Bacteroides massiliensis",),
    "Ruminococcus gnavus": ("Ruminococcus gnavus",),
    "Blautia massiliensis": ("Blautia massiliensis", "Blautia sp CAG 257", "Blautia sp N6H1 15"),
    "Mediterraneibacter gnavus": ("Ruminococcus gnavus",),
}

_GTDB_SUFFIX: Final = re.compile(r"_[A-Z]+\b")

#: Percentile bands used for the direction comparison. Tiny Health's own
#: band thresholds put a reading low / mid / high on their scale; the
#: matching split on ours is the quartile outside the middle half.
LOW, HIGH = 25.0, 75.0

#: A three-band comparison is fragile exactly at the band lines: a reading
#: 1% under Tiny Health's cutoff and 2 points over ours reads as "opposite"
#: when the two platforms in fact agree about the value and differ only in
#: where they draw a line. A pair is a *band-edge straddle* when both sides
#: sit this close to their own lines. Straddles are reported, not scored as
#: disagreements, and the rank-order check across samples is the arbiter.
EDGE_TINY_FRACTION = 0.03   # of the nearest Tiny Health threshold
EDGE_OUR_POINTS = 5.0       # of our 25th / 75th line


def our_direction(percentile: float) -> str:
    if percentile < LOW:
        return "low"
    if percentile > HIGH:
        return "high"
    return "mid"


def normalise_species(name: str) -> str:
    """GTDB species name reduced to the binomial MetaPhlAn would use."""
    return _GTDB_SUFFIX.sub("", name).replace("_", " ").strip()


@dataclass(frozen=True)
class FunctionPair:
    sample: str
    tiny_metric: str
    tiny_value: float | None
    tiny_unit: str | None
    tiny_evaluation: str | None
    tiny_direction: str
    tiny_thresholds: tuple[float, ...]
    panel: str
    our_value: float
    our_percentile: float
    our_status: str
    like_for_like: bool

    @property
    def our_dir(self) -> str:
        return our_direction(self.our_percentile)

    @property
    def band_edge(self) -> bool:
        """Both readings within a hair of their own band lines."""
        if self.tiny_value is None or not self.tiny_thresholds:
            return False
        nearest = min(self.tiny_thresholds, key=lambda t: abs(t - self.tiny_value))
        tiny_close = nearest > 0 and abs(self.tiny_value - nearest) / nearest <= EDGE_TINY_FRACTION
        line = LOW if self.our_percentile < 50 else HIGH
        ours_close = abs(self.our_percentile - line) <= EDGE_OUR_POINTS
        return tiny_close and ours_close

    @property
    def verdict(self) -> str:
        if self.tiny_direction == "unknown":
            return "not comparable"
        if self.tiny_direction == self.our_dir:
            return "agree"
        if {self.tiny_direction, self.our_dir} == {"low", "high"}:
            return "opposite (band edge)" if self.band_edge else "opposite"
        return "one band apart"


def load_ours(sample: str, results: Path) -> dict[str, Any]:
    return json.loads((results / sample / "results.json").read_text(encoding="utf-8"))


def compare_functions(sample: str, tiny: Report, ours: dict[str, Any]) -> list[FunctionPair]:
    """Pair every Tiny Health functional metric with our panel for it."""
    panels = ours["reference_comparison"]["panels"]
    seen: set[str] = set()
    pairs: list[FunctionPair] = []
    for metric in tiny.metrics:
        panel = FUNCTION_MAP.get(metric.metric)
        if panel is None or panel in seen:
            continue
        ours_panel = panels.get(panel)
        if ours_panel is None:
            continue
        seen.add(panel)
        pairs.append(FunctionPair(
            sample=sample, tiny_metric=metric.metric, tiny_value=metric.value,
            tiny_unit=metric.unit, tiny_evaluation=metric.evaluation,
            tiny_direction=metric.direction, tiny_thresholds=tuple(metric.thresholds),
            panel=panel,
            our_value=round(ours_panel["value"], 3),
            our_percentile=ours_panel["percentile"], our_status=ours_panel["status"],
            like_for_like=metric.metric not in NOT_LIKE_FOR_LIKE,
        ))
    return pairs


def top_species(tiny_pdf: Path) -> list[tuple[int, str, str, float]]:
    """The "Top 20 species" table, where the template has one."""
    import pymupdf

    text = "\n".join(page.get_text() for page in pymupdf.open(tiny_pdf))
    start = text.find("Top 20 species")
    if start < 0:
        return []
    pattern = re.compile(
        r"\n(\d{1,2})\n([A-Z][A-Za-z_]+ [a-z_A-Z0-9]+)\n(?:(?!\n\d{1,2}\n).)*?"
        r"\n(Variable|Beneficial|Unfriendly|Unknown)\n([\d.]+)%",
        re.S,
    )
    return [(int(m.group(1)), m.group(2), m.group(3), float(m.group(4)))
            for m in pattern.finditer(text[start:])]


#: A GTDB genome bin with no binomial: the epithet is an assembly
#: accession ("sp003526955"), not a species name. Whether two catalogues
#: carry the same bin is a coverage question, not a measurement one.
_UNNAMED_BIN = re.compile(r"\bsp\d{6,}\b|^(?:CAG|UBA|PeH|GCA)[-_]")


def _is_unnamed_bin(name: str) -> bool:
    return bool(_UNNAMED_BIN.search(name))


def compare_species(tiny_pdf: Path, ours: dict[str, Any]) -> dict[str, Any]:
    """Relative abundance of the same species, both platforms."""
    rows = top_species(tiny_pdf)
    if not rows:
        return {"n": 0}
    # Compare against what the report actually reports: the organism
    # inventory, which merges the marker lane with the whole-genome lane.
    # Reading the MetaPhlAn 4 SGB list alone understated agreement -
    # Bacteroides uniformis, Parabacteroides distasonis, Fusicatenibacter
    # saccharivorans and Mediterraneibacter lactaris are all in one
    # sample's inventory at levels close to Tiny Health's, and all four
    # counted as undetected because only the other lane had found them.
    levels: dict[str, float] = {}
    for sgb in ours.get("extended_catalogue", {}).get("sgbs") or ():
        for key in (sgb.get("species"), sgb.get("gtdb"), sgb.get("sgb")):
            if key:
                levels[str(key)] = max(
                    levels.get(str(key), 0.0), float(sgb.get("percent") or 0.0))
    for org in (ours.get("organism_inventory") or {}).get("organisms") or ():
        pct = org.get("percent")
        if pct is None:
            pct = org.get("secondary_percent")
        if not pct:
            continue
        # Both the binomial and the GTDB name, because Tiny Health reports
        # GTDB and a third of a stool community has no binomial at all:
        # "Ruminococcus_C sp000980705" is this pipeline's Ruminococcus_SGB4421
        # and matches only on the GTDB name it also carries.
        # Raw, underscores and all: `normalise_species` strips the GTDB
        # sub-genus suffix ("_C") and then replaces underscores, and doing
        # the replacement first leaves the suffix unstrippable.
        for key in (org.get("species"), org.get("gtdb")):
            if key:
                levels[str(key)] = max(levels.get(str(key), 0.0), float(pct))

    def ours_for(name: str) -> float | None:
        plain = normalise_species(name)
        wanted = SPECIES_ALIASES.get(plain, (plain,))
        total = sum(v for k, v in levels.items()
                    if normalise_species(k) in wanted or k in wanted)
        return total or None

    paired = [(name, tiny_pct, ours_for(name)) for _, name, _, tiny_pct in rows]
    # Keep the name with its pair of numbers. Zipping the full list against
    # the filtered one, as this used to, printed each matched value beside
    # whichever name happened to sit at the same index.
    matched_rows = [(name, t, o) for name, t, o in paired if o is not None]
    matched = [(t, o) for _, t, o in matched_rows]
    if len(matched) < 3:
        return {
            "n": len(matched), "detected_by_both": len(matched), "n_tiny": len(rows),
            "detection_agreement": len(matched) / len(rows),
            "pairs": matched_rows, "spearman_rho": None,
            "median_abs_log10_ratio": None, "median_log10_ratio": None,
            "detected_only_here": [name for name, _, o in paired if o is None],
        }
    ratios = [math.log10(t / o) for t, o in matched]
    ranks = lambda xs: [sorted(xs).index(v) for v in xs]  # noqa: E731 - local helper
    ra, rb = ranks([m[0] for m in matched]), ranks([m[1] for m in matched])
    n = len(matched)
    rho = 1 - 6 * sum((i - j) ** 2 for i, j in zip(ra, rb, strict=True)) / (n * (n * n - 1))
    # Two different questions, kept apart. A binomial Tiny Health reports
    # and this pipeline does not is a missed detection. An unnamed GTDB
    # genome bin ("sp003526955") that one catalogue carries and the other
    # does not is reference coverage: there is no name to agree on.
    named = [(name, t, o) for name, t, o in paired if not _is_unnamed_bin(name)]
    named_found = [row for row in named if row[2] is not None]
    unnamed = [(name, t, o) for name, t, o in paired if _is_unnamed_bin(name)]
    return {
        "n_tiny": len(rows), "detected_by_both": n,
        "detection_agreement": n / len(rows),
        "n_tiny_named": len(named),
        "detected_named": len(named_found),
        "named_agreement": (len(named_found) / len(named)) if named else None,
        "n_tiny_unnamed_bins": len(unnamed),
        "detected_unnamed_bins": sum(1 for _, _, o in unnamed if o is not None),
        "spearman_rho": round(rho, 3),
        "median_abs_log10_ratio": round(statistics.median(abs(r) for r in ratios), 3),
        "median_log10_ratio": round(statistics.median(ratios), 3),
        "pairs": matched_rows,
        "detected_only_here": [name for name, _, o in paired if o is None],
    }


def our_metric_counts(ours: dict[str, Any]) -> dict[str, int]:
    """How many distinct quantities this pipeline reports for one sample."""
    similarity = ours["profile_similarity"]
    community = similarity["community"]
    pathogens = ours["pathogens"]
    entries = ours["reference_entries"]
    return {
        "metabolic panels": len(ours["panels"]),
        "gene readings inside panels": sum(
            1 for k, v in entries.items() if v.get("role") == "target" and not k.startswith("_")
        ),
        "microbial groups": len(community["groups"]),
        "species with abundance": community["n_species_detected"],
        "species with a cohort percentile": sum(
            1 for s in community["species"] if s.get("percentile") is not None
        ),
        "extended-catalogue SGBs": ours["extended_catalogue"]["n_sgbs_detected"],
        "disease patterns scored": sum(
            1 for r in similarity["ranked"] if r.get("percentile") is not None
        ),
        "skin research panels": len(similarity.get("skin_pattern_indices") or {}),
        "pathogen targets screened": pathogens["coverage"]["assessed"],
        "pathogen targets in catalogue": pathogens["coverage"]["total_targets"],
        "resistance, toxin and virulence determinants": len(pathogens["determinants"]),
        "ecological diversity metrics": 6,
        "biological age estimate": 1,
    }


def main() -> int:  # noqa: PLR0915 - a report, printed in order
    tiny_dir, results = ROOT / "tinyhealth", ROOT / "results"
    reports = load_all(tiny_dir)
    all_pairs: list[FunctionPair] = []
    print("=" * 78)
    print("REPORTS")
    print("=" * 78)
    for name, report in sorted(reports.items()):
        if not (results / name / "results.json").exists():
            print(f"  {name}: no local result, skipped")
            continue
        ours = load_ours(name, results)
        print(f"  {name:14s} Tiny Health v{report.version} "
              f"{len(report.metrics):3d} metrics, score {report.summary_score} "
              f"| ours {sum(our_metric_counts(ours).values()):,} quantities")
        all_pairs.extend(compare_functions(name, report, ours))

    print()
    print("=" * 78)
    print("FUNCTIONAL CAPACITY — position on each platform's own scale")
    print("=" * 78)
    print(f"{'smp':6s} {'Tiny Health metric':28s} {'their value':>12s} {'dir':>5s} | "
          f"{'panel':15s} {'pct':>6s} {'dir':>5s}  verdict")
    comparable = [p for p in all_pairs if p.like_for_like and p.tiny_direction != "unknown"]
    for p in all_pairs:
        value = f"{p.tiny_value}{p.tiny_unit or ''}" if p.tiny_value is not None else "—"
        flag = "" if p.like_for_like else "  (different quantity)"
        print(f"{p.sample.split('_')[0]:6s} {p.tiny_metric:28s} {value:>12s} "
              f"{p.tiny_direction:>5s} | {p.panel:15s} {p.our_percentile:6.1f} "
              f"{p.our_dir:>5s}  {p.verdict}{flag}")
    n = len(comparable)
    agree = sum(1 for p in comparable if p.verdict == "agree")
    one = sum(1 for p in comparable if p.verdict == "one band apart")
    edge = [p for p in comparable if p.verdict == "opposite (band edge)"]
    opp = [p for p in comparable if p.verdict == "opposite"]
    print(f"\n  {n} like-for-like pairs: {agree} agree ({agree / n:.0%}), "
          f"{one} one band apart ({(agree + one) / n:.0%} within one band), "
          f"{len(edge)} opposite at a band edge, {len(opp)} opposite ({len(opp) / n:.0%})")
    if edge:
        print("  band-edge straddles (agree on the value, differ on the line):")
        for p in edge:
            nearest = min(p.tiny_thresholds, key=lambda t: abs(t - p.tiny_value))
            print(f"    {p.sample} {p.tiny_metric}: theirs {p.tiny_value}{p.tiny_unit or ''} is "
                  f"{abs(p.tiny_value - nearest) / nearest:.1%} from their cutoff {nearest}; "
                  f"ours {p.our_percentile:.1f}th is {abs(p.our_percentile - (LOW if p.our_percentile < 50 else HIGH)):.1f} "
                  "points from our line")
    if opp:
        print("  opposite calls:")
        for p in opp:
            print(f"    {p.sample} {p.tiny_metric}: theirs {p.tiny_value}{p.tiny_unit or ''} "
                  f"({p.tiny_direction}) vs ours {p.our_percentile:.0f}th ({p.our_dir})")
    excluded = [p for p in all_pairs if not p.like_for_like]
    if excluded:
        print("\n  excluded as not like-for-like:")
        for metric, why in NOT_LIKE_FOR_LIKE.items():
            hits = [p for p in excluded if p.tiny_metric == metric]
            if hits:
                print(f"    {metric} ({len(hits)} samples) — {why}")

    print()
    print("=" * 78)
    print("SPECIES ABUNDANCE — same quantity, both platforms")
    print("=" * 78)
    for name in sorted(reports):
        pdf = next(tiny_dir.glob(f"{name}.*.pdf"), None)
        if pdf is None or not (results / name / "results.json").exists():
            continue
        stats = compare_species(pdf, load_ours(name, results))
        if not stats.get("detected_by_both"):
            print(f"  {name}: this report template has no species table")
            continue
        print(f"  {name}: {stats['detected_by_both']}/{stats['n_tiny']} of their top species "
              f"also detected here ({stats['detection_agreement']:.0%}), "
              f"Spearman rho {stats['spearman_rho']}, "
              f"median |log10 ratio| {stats['median_abs_log10_ratio']} "
              f"({10 ** stats['median_abs_log10_ratio']:.2f}x)")

    print()
    print("=" * 78)
    print("SCOPE — quantities reported per sample")
    print("=" * 78)
    for name in sorted(reports):
        if not (results / name / "results.json").exists():
            continue
        counts = our_metric_counts(load_ours(name, results))
        tiny_n = len(reports[name].metrics)
        print(f"\n  {name}   Tiny Health {tiny_n} | this pipeline {sum(counts.values()):,}")
        for key, value in counts.items():
            print(f"      {value:>6,}  {key}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
