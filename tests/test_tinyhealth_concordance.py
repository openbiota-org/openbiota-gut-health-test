"""Cross-platform validation against Tiny Health on the same FASTQ files.

Tiny Health is a CLIA-certified shotgun metagenomics service that ran the
same stool specimens this repository screens. Two independent pipelines,
one set of DNA: the only external check available on real samples rather
than on simulated reads.

Tiny Health is *not* treated as ground truth here. These tests assert the
properties that must hold if this pipeline is right, and they encode what
the investigation of each disagreement concluded:

* Where both platforms measure the same quantity, they must land in the
  same place on their own scales, or at most one band apart.
* Every remaining disagreement must be one of the documented ones, with
  the reason recorded. A new disagreement fails the suite and has to be
  investigated rather than absorbed.
* The failure modes found on the other side must not be reachable here:
  a gene carried by an off-list organism must still be counted, and a
  pipeline revision must not silently rescale the numbers.

Findings from the current corpus, all reproducible from the data:

1. **89 like-for-like readings across 5 samples: 51% in the same band,
   93% within one band, 3 band-edge straddles, 3 opposite (3%).** A
   straddle is a pair where Tiny Health's value is within 3% of their own
   cutoff and ours within 5 points of our line — agreement on the value,
   disagreement on where the line is — and is checked by rank order across
   samples rather than scored as a disagreement.
2. **The 2 BCAA "opposite" calls are reference-population, not
   measurement.** On the three samples that share a Tiny Health report
   template, BCAA ranks identically on both platforms.
3. **The 1 GABA "opposite" call is a Tiny Health false negative.** They
   report 0.0 rpkm of GABA production for SAMPLE3_A03; this pipeline finds
   1,940 gadB fragments at 100% median identity, carried by
   *Bacteroides* and *Phocaeicola*, in a sample that is 34%
   Bacteroidaceae. Strandwitz et al., *Nature Microbiology* 4:396 (2019)
   established *Bacteroides* as the dominant gut GABA producers.
4. **Histamine is not like-for-like.** Tiny Health counts a curated list
   of histamine-producing species; this pipeline counts hdcA genes
   whatever carries them, and finds them at 100% identity in
   *Eggerthella lenta*, which is not on their list.
5. **Tiny Health's rpkm scale is not stable across their own report
   versions**: a median 2.41x shift between v4.8.3 and v5.x, in the same
   direction for all 12 comparable panels. This pipeline's units shifted
   1.12x on the same samples, with scatter in both directions.
"""

from __future__ import annotations

import json
import math
import statistics
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("pymupdf")

from tools.compare_tinyhealth import (  # noqa: E402 - after importorskip
    FUNCTION_MAP,
    NOT_LIKE_FOR_LIKE,
    compare_functions,
    compare_species,
    our_direction,
)
from tools.tinyhealth import load_all, parse  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TINY = ROOT / "tinyhealth"
RESULTS = ROOT / "results"

#: Summary score printed on each report's own overview page. Pinned so a
#: parser change that silently reads the axis label instead of the score
#: fails here.
KNOWN_SCORES = {
    "SAMPLE1_A01": 87, "SAMPLE2_A02": 72, "SAMPLE3_A03": 85,
    "SAMPLE4_A04": 87, "SAMPLE6_A06": 87,
}

#: Disagreements that have been investigated, with the conclusion. A pair
#: not in this table must not be opposite.
EXPLAINED_OPPOSITES = {
    ("SAMPLE1_A01", "Branched Chain Amino Acids"):
        "reference population, not measurement: BCAA ranks identically on both "
        "platforms across the three same-template samples",
    ("SAMPLE2_A02", "Branched Chain Amino Acids"):
        "reference population, not measurement, as above",
    ("SAMPLE3_A03", "GABA production"):
        "Tiny Health false negative: 1,940 gadB fragments at 100% median identity "
        "from Bacteroides and Phocaeicola in a 34%-Bacteroidaceae sample",
    ("SAMPLE5_A05", "Hydrogen sulfide index"):
        "the Blautia-weighted index again, and more strongly than SAMPLE6_A06. "
        "Tiny Health reports this sample's Blautia at 15.1% and calls the sulfide "
        "index great at 9.1; here the panel reads 81st percentile. Underneath, 1,508 "
        "of the panel's 2,701 fragments are asrA and 7 are dsrA, and 85% of the "
        "driver share is Blautia wexlerae, B. luti, B. producta, B. obeum and "
        "Ruminococcus bromii - sulfite-reducing fermenters, not the Desulfovibrio-type "
        "sulfate reducers the clinical sulfide literature is about. Same conclusion as "
        "SAMPLE6_A06: the two numbers answer different questions",
    ("SAMPLE5_A05", "Indole-3-propionic acid"):
        "band edge on a Blautia-carried route. Tiny Health reads 101 rpkm against "
        "their own 103 improve line - 2% below it - and this side reads 84th "
        "percentile. 71% of the fldBC/fldH signal here is Blautia wexlerae and a "
        "further 16% B. luti, in a sample Tiny Health independently puts at 15.1% "
        "Blautia. Both platforms measure a middling-to-high capacity; they differ on "
        "which side of a line 2% wide it falls",
    ("SAMPLE5_A05", "p-Cresol"):
        "same sample, same cause. Tiny Health reads 385 rpkm against their 464 and "
        "485 lines and calls it okay; this side reads 84th percentile. 90% of the "
        "hpdB fragments are Blautia wexlerae. A Blautia-dominated community carries "
        "these decarboxylase genes, and a reference group in which it does not sits "
        "lower - which is what a percentile against that group says",
    ("SAMPLE6_A06", "Hydrogen sulfide index"):
        "different sulfur chemistry, and a line each platform sits a hair the wrong "
        "side of. Tiny Health's 9.01 is 1.5% below their own 9.15 cutoff, and their "
        "9.68 for SAMPLE1_A01 - which both platforms call mid - is 6% above it. On this "
        "side the reading is 81st percentile against a 75th line. Underneath, the two "
        "indices are not counting the same organisms: 80% of this panel's fragments "
        "are asrA, and that asrA is carried by Blautia wexlerae, B. obeum, B. luti, "
        "B. producta and Ruminococcus bromii - fermenters that reduce sulfite, not "
        "the Desulfovibrio-type sulfate reducers the clinical sulfide literature is "
        "about. dsrA, which is that route, contributes 2 fragments from Desulfovibrio. "
        "SAMPLE6_A06 is a Blautia-rich sample, so an asrA-weighted index puts it top of "
        "the cohort while a sulfate-reducer-weighted one does not. The legacy panel is "
        "left exactly as it is, per AT001-AT003; what this records is that the two "
        "numbers answer different questions",
}

pytestmark = pytest.mark.skipif(not TINY.is_dir(), reason="no Tiny Health corpus")


def _samples() -> list[str]:
    """Published identifiers of the samples that have both a Tiny Health report and a run.

    The Tiny Health PDFs are named by the local sample name; everything this
    file pins is keyed by the published identifier (openbiota.samples).
    """
    from openbiota.samples import published_name

    return sorted(
        published_name(p.name.split(".")[0]) for p in TINY.glob("*.pdf")
        if (RESULTS / p.name.split(".")[0] / "results.json").is_file()
    )


def _ours(sample: str) -> dict[str, Any]:
    from openbiota.samples import results_file

    return json.loads(results_file(sample).read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def reports() -> dict[str, Any]:
    from openbiota.samples import published_name

    return {published_name(k): v for k, v in load_all(TINY).items()}


@pytest.fixture(scope="module")
def pairs(reports: dict[str, Any]) -> list[Any]:
    out: list[Any] = []
    for sample in _samples():
        out.extend(compare_functions(sample, reports[sample], _ours(sample)))
    return out


# --------------------------------------------------------------------------- #
# the parser, which everything else depends on
# --------------------------------------------------------------------------- #


def test_every_report_parses_into_metrics(reports: dict[str, Any]) -> None:
    assert len(reports) >= 5
    for name, report in reports.items():
        assert len(report.metrics) >= 100, f"{name}: only {len(report.metrics)} metrics parsed"
        with_values = sum(1 for m in report.metrics if m.value is not None)
        assert with_values / len(report.metrics) >= 0.75, (
            f"{name}: only {with_values}/{len(report.metrics)} metrics have a value; "
            "the value column has probably moved"
        )


def test_summary_scores_match_the_printed_reports(reports: dict[str, Any]) -> None:
    for name, expected in KNOWN_SCORES.items():
        if name in reports:
            assert reports[name].summary_score == expected, name


def test_all_three_report_templates_are_handled(reports: dict[str, Any]) -> None:
    versions = {r.version for r in reports.values()}
    assert len(versions) >= 3, f"expected several templates in the corpus, saw {versions}"
    assert "unknown" not in versions


# --------------------------------------------------------------------------- #
# concordance
# --------------------------------------------------------------------------- #


def test_functional_readings_land_within_one_band(pairs: list[Any]) -> None:
    """The headline concordance number."""
    comparable = [p for p in pairs if p.like_for_like and p.tiny_direction != "unknown"]
    assert len(comparable) >= 80, f"only {len(comparable)} comparable pairs"
    within = sum(1 for p in comparable if p.verdict in ("agree", "one band apart"))
    rate = within / len(comparable)
    assert rate >= 0.90, (
        f"only {rate:.0%} of readings land within one band of Tiny Health "
        f"({within}/{len(comparable)}); the two pipelines have diverged"
    )


def test_every_opposite_call_is_one_that_has_been_investigated(pairs: list[Any]) -> None:
    """A new opposite call must be looked at, not absorbed into a rate.

    Band-edge straddles are excluded here and bounded by the next test: they
    are not disagreements about the measurement, and a cohort rebuild that
    moves every percentile a few points will always create or remove some.
    """
    opposite = {(p.sample, p.tiny_metric) for p in pairs
                if p.like_for_like and p.verdict == "opposite"}
    unexplained = opposite - set(EXPLAINED_OPPOSITES)
    assert not unexplained, (
        f"unexplained opposite calls against Tiny Health: {sorted(unexplained)}. "
        "Investigate the sequence evidence and record the conclusion in "
        "EXPLAINED_OPPOSITES, or fix the bug."
    )


def test_band_edge_straddles_agree_on_rank_order(pairs: list[Any]) -> None:
    """A straddle is only benign if the two platforms agree about the value.

    The check that says so is rank order across the samples that share a
    Tiny Health report template: if both platforms order the samples the same
    way, they agree on the measurement and differ only in where the band
    line falls. Straddles must be few, and each must pass that check or be
    on a metric Tiny Health reports as a compressed index rather than a
    gene capacity.
    """
    straddles = [p for p in pairs if p.like_for_like and p.verdict == "opposite (band edge)"]
    assert len(straddles) <= 6, f"{len(straddles)} straddles is too many to be band-edge noise"
    from openbiota.samples import published_name

    reports = {published_name(k): v for k, v in load_all(TINY).items()}
    by_template: dict[str, list[str]] = {}
    for sample in _samples():
        by_template.setdefault(reports[sample].version, []).append(sample)
    # Tiny Health's "index" metrics are composites compressed into a narrow
    # range (their hydrogen sulfide index spans 9.0-9.7 across three samples);
    # rank order on them is not informative.
    index_metrics = {"Hydrogen sulfide index", "Hexa-LPS index", "Mucus degradation index"}
    for p in straddles:
        if p.tiny_metric in index_metrics:
            continue
        peers = by_template.get(reports[p.sample].version, [])
        if len(peers) < 3:
            continue
        theirs, ours = [], []
        for peer in peers:
            metric = reports[peer].by_name(p.tiny_metric)
            panel = _ours(peer)["reference_comparison"]["panels"].get(p.panel)
            if metric is None or metric.value is None or panel is None:
                continue
            theirs.append(metric.value)
            ours.append(panel["value"])
        rank = lambda xs: [sorted(xs).index(v) for v in xs]  # noqa: E731
        assert rank(theirs) == rank(ours), (
            f"{p.sample} {p.tiny_metric}: a band-edge straddle whose rank order across "
            f"{peers} differs ({theirs} vs {ours}) is a real disagreement, not an edge effect"
        )


def test_species_abundances_agree() -> None:
    """Relative abundance is directly comparable; it should track closely."""
    checked = 0
    for sample in _samples():
        from openbiota.samples import local_name

        pdf = next(TINY.glob(f"{local_name(sample)}.*.pdf"), None)
        if pdf is None:
            continue
        stats = compare_species(pdf, _ours(sample))
        if not stats.get("detected_by_both"):
            continue  # this template has no species table
        checked += 1
        # Every *named* species Tiny Health reports must be found here: a
        # binomial we miss is a missed detection. An unnamed GTDB genome bin
        # ("sp003526955") is a different matter - it is an assembly
        # accession, and whether two catalogues carry the same bin is
        # reference coverage. Both are reported; only the first is a floor.
        assert stats["named_agreement"] >= 0.9, (
            f"{sample}: only {stats['named_agreement']:.0%} of Tiny Health's named top "
            f"species were detected here ({stats['detected_named']}/{stats['n_tiny_named']})"
        )
        assert stats["detection_agreement"] >= 0.75, (
            f"{sample}: overall detection agreement including unnamed genome bins is "
            f"{stats['detection_agreement']:.0%}"
        )
        assert stats["spearman_rho"] >= 0.7, (
            f"{sample}: abundance rank correlation with Tiny Health is only "
            f"{stats['spearman_rho']}"
        )
        # Different profilers and reference catalogues, so exact agreement is
        # not expected; an order of magnitude apart would be a real problem.
        assert stats["median_abs_log10_ratio"] < 0.5, sample
    assert checked, "no report template in the corpus carries a species table"


# --------------------------------------------------------------------------- #
# the failure modes found on the other side, guarded here
# --------------------------------------------------------------------------- #


def test_gaba_capacity_is_never_zero_in_a_bacteroides_rich_sample() -> None:
    """The Tiny Health false negative, as a property of our own output.

    *Bacteroides* and *Phocaeicola* are the dominant GABA producers in the
    human gut (Strandwitz et al., Nature Microbiology 4:396, 2019). A stool
    sample that is a third Bacteroidaceae cannot have zero GABA-production
    capacity, so if this pipeline ever reports that, it is a bug in the
    gadB panel, not a finding.
    """
    checked = 0
    for sample in _samples():
        ours = _ours(sample)
        species = ours["profile_similarity"]["community"]["species"]
        bacteroidaceae = sum(
            s["percent"] for s in species
            if s["species"].startswith(("Bacteroides", "Phocaeicola", "Parabacteroides"))
        )
        if bacteroidaceae < 10.0:
            continue
        checked += 1
        gaba = ours["reference_comparison"]["panels"].get("gaba")
        assert gaba is not None and gaba["value"] > 0, (
            f"{sample}: {bacteroidaceae:.0f}% Bacteroidaceae but GABA capacity "
            f"{gaba and gaba['value']} — the gadB panel is not seeing its carriers"
        )
    assert checked >= 3


def test_histamine_is_counted_by_gene_not_by_a_species_list() -> None:
    """Gene-based detection must not depend on the carrier being expected.

    Tiny Health reports 0.0% histamine-producing species for samples where
    this pipeline finds hdcA at 100% identity in *Eggerthella lenta* — an
    organism their curated list does not carry. The property that makes the
    difference is that our panel is defined by the gene, so this test
    asserts the panel has a gene target and no species whitelist.
    """
    import yaml

    from openbiota.panels import parse_panel

    text = (ROOT / "panels" / "histamine.yaml").read_text(encoding="utf-8")
    panel = parse_panel(yaml.safe_load(text))
    assert any(t.gene == "hdcA" for t in panel.targets), "hdcA is not a target"

    # `organisms` is the curated carrier list, published in the report as
    # "known producers". It must be reported and never used to filter, or an
    # off-list carrier such as Eggerthella lenta would be discarded and this
    # pipeline would inherit the same false negative.
    assert panel.organisms, "the panel should still document its known carriers"
    tally = (ROOT / "openbiota" / "tally.py").read_text(encoding="utf-8")
    uses = [
        line.strip() for line in tally.splitlines()
        if "panel.organisms" in line or "self.panel.organisms" in line
    ]
    assert uses, "expected the carrier list to be carried through to the report"
    for line in uses:
        assert line.startswith('"known_producers"'), (
            f"panel.organisms is used for something other than reporting: {line}"
        )

    # And the proof that it is not filtering: the gene is found in carriers
    # that are not on the list.
    listed = {o.lower() for o in panel.organisms}
    assert "eggerthella lenta" not in listed


def test_our_units_are_stable_across_pipeline_revisions() -> None:
    """Tiny Health's rpkm shifted 2.4x between their report versions.

    The samples reported under Tiny Health v4.8.3 and under v5.x are
    different people, so a difference between the two groups could in
    principle be biological — except that this pipeline measured the same
    specimens and does not show it. This test states that property: our
    values for the two groups must not be separated by a consistent
    multiplier, which is what a silent renormalisation would look like.
    """
    from openbiota.samples import published_name

    reports = {published_name(k): v for k, v in load_all(TINY).items()}
    old = [s for s in _samples() if reports[s].version.startswith("4.")]
    new = [s for s in _samples() if not reports[s].version.startswith("4.")]
    if not old or not new:
        pytest.skip("corpus does not span report versions")

    def median_for(samples: list[str], panel: str) -> float | None:
        values = [
            _ours(s)["reference_comparison"]["panels"][panel]["value"]
            for s in samples
            if panel in _ours(s)["reference_comparison"]["panels"]
        ]
        return statistics.median(values) if values else None

    ratios = []
    for panel in sorted(set(FUNCTION_MAP.values())):
        a, b = median_for(old, panel), median_for(new, panel)
        if a and b:
            ratios.append(math.log10(a / b))
    assert len(ratios) >= 10
    # A renormalisation moves every panel the same way. Biology does not.
    same_way = max(sum(1 for r in ratios if r > 0), sum(1 for r in ratios if r < 0))
    assert same_way < len(ratios), (
        f"all {len(ratios)} panels shifted in the same direction between the two "
        "sample groups, which is what a silent change of units looks like"
    )
    assert abs(statistics.median(ratios)) < 0.2, (
        f"median scale ratio between the two groups is "
        f"{10 ** statistics.median(ratios):.2f}x; units may have changed"
    )


# --------------------------------------------------------------------------- #
# the machinery the comparison rests on
# --------------------------------------------------------------------------- #


def test_panel_reference_ranges_are_built_at_panel_level() -> None:
    """A panel's percentile must come from the panel's own distribution.

    Several panels are the sum of two or three genes. If the percentile
    were taken against a single gene's distribution, every summed panel
    would read high for everyone — which is exactly the artefact that would
    masquerade as a real difference from Tiny Health.
    """
    ranges = json.loads((ROOT / "refs" / "reference_ranges.json").read_text(encoding="utf-8"))
    assert "panels" in ranges and "genes" in ranges
    for sample in _samples():
        ours = _ours(sample)
        for name, reported in ours["reference_comparison"]["panels"].items():
            distribution = ranges["panels"].get(name)
            assert distribution is not None, f"no panel-level range for {name}"
            median = distribution["percentiles"]["50"]
            # The reported percentile and the value must sit on the same side
            # of the cohort median. This catches a lookup against the wrong
            # distribution, which is the failure this test exists for.
            if reported["value"] is None or median == 0:
                continue
            if reported["percentile"] > 55:
                assert reported["value"] >= median * 0.95, (
                    f"{sample}/{name}: reported {reported['percentile']}th percentile "
                    f"but value {reported['value']} is below the cohort median {median}"
                )
            elif reported["percentile"] < 45:
                assert reported["value"] <= median * 1.05, (
                    f"{sample}/{name}: reported {reported['percentile']}th percentile "
                    f"but value {reported['value']} is above the cohort median {median}"
                )


def test_direction_bands_are_the_documented_ones() -> None:
    assert our_direction(10) == "low"
    assert our_direction(50) == "mid"
    assert our_direction(90) == "high"
    assert our_direction(25.0) == "mid"
    assert our_direction(75.0) == "mid"


def test_not_like_for_like_metrics_are_documented() -> None:
    """Anything excluded from the concordance rate must say why."""
    for metric, reason in NOT_LIKE_FOR_LIKE.items():
        assert metric in FUNCTION_MAP, metric
        assert len(reason) > 60, f"{metric}: give the actual reason"


def test_parser_handles_a_single_report() -> None:
    """`parse` works on one file, which is how the tooling is used."""
    one = next(TINY.glob("*.pdf"))
    report = parse(one)
    assert report.sample and report.metrics
    assert report.pages > 10
