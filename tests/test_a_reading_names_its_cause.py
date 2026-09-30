"""A reading must say what is behind it, and must not overstate a trace.

The colorectal virulence panel printed ABOVE AVERAGE at the 76th
percentile, and beneath the heading "What is driving this" it named
*Escherichia coli* as the whole of its signal and marked it "not found".
Four fragments. Neither of the other two genes present. None of the three
organisms the panel exists to detect in the sample at all.

Both halves of that were wrong and both are tested here: the position,
which came from a cohort where three quarters carry none, and the
attribution, which named an organism nobody found.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.drivers import Driver, Drivers
from openbiota.pdfreport import status_for
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401

REPO = Path(__file__).resolve().parent.parent
SAMPLES = ("SAMPLE2_A02", "SAMPLE1_A01", "SAMPLE3_A03", "SAMPLE4_A04", "SAMPLE6_A06")


def _ranges_path():
    from openbiota.samples import reference_ranges

    path = reference_ranges()
    if not path.is_file():
        pytest.skip("refs/reference_ranges.json is not built in this checkout (make cohort)")
    return path


def _results(sample: str) -> dict[str, Any]:
    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


# --- a position against a mostly-empty cohort -----------------------------


def test_a_trace_in_a_mostly_absent_cohort_is_not_called_high() -> None:
    """Where the median is zero, any detection outranks it by arithmetic."""
    status = status_for(percentile=76.3, higher_means="adverse", cohort_mostly_absent=True)
    assert status.label == "present"
    assert status.label != "above average"
    assert not status.assessed, "a position earned by zero-inflation is not a verdict"
    assert "carries none of this" in status.note


def test_the_same_position_in_a_real_spread_still_reads_normally() -> None:
    """The gate must not flatten readings that have a distribution."""
    status = status_for(percentile=76.3, higher_means="adverse")
    assert status.label == "above average"
    assert status.assessed


def test_nothing_detected_outranks_the_gate() -> None:
    """A measured zero is still 'not detected', not 'present'."""
    status = status_for(percentile=0.0, detected=False, cohort_mostly_absent=True)
    assert status.label == "not detected"


@pytest.mark.parametrize("sample", SAMPLES)
def test_no_reading_claims_a_verdict_off_a_zero_median_cohort(sample: str) -> None:
    """The report-level guarantee, checked against the shipped ranges."""
    ranges = json.loads(
        _ranges_path().read_text(encoding="utf-8")
    )
    verdicts = {"notably low", "low", "below average", "above average", "high", "notably high"}
    comparison = _results(sample)["reference_comparison"]["panels"]
    offenders = []
    for name, row in comparison.items():
        spec = (ranges["panels"].get(name) or {}).get("percentiles") or {}
        if spec.get("50") not in (0, 0.0):
            continue
        if str(row.get("status")) in verdicts:
            offenders.append(f"{name} reads {row['status']!r} against an all-zero median")
    assert not offenders, "; ".join(offenders)


# --- attribution ----------------------------------------------------------


def _drivers(*rows: tuple[str, int, bool, float | None]) -> Drivers:
    top = tuple(
        Driver(organism=n, fragments=f, share=f / max(sum(r[1] for r in rows), 1),
               in_sample=ins, sample_percent=pct, sample_class=None,
               genus_in_sample=False)
        for n, f, ins, pct in rows
    )
    return Drivers(panel="p", entry_ids=(), total_fragments=sum(r[1] for r in rows), top=top)


def _verdict_text(drv: Drivers) -> str:
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota.pdfreport import _drivers_verdict, _styles

    story: list[Any] = []
    _drivers_verdict(story, _styles(), drv)
    if not story:
        return ""
    SimpleDocTemplate("/tmp/_cause.pdf", pagesize=A4).build(story)
    return " ".join(
        "".join(p.extract_text() for p in pypdf.PdfReader("/tmp/_cause.pdf").pages).split()
    )


def test_a_reading_with_no_organism_present_says_so_first() -> None:
    """The exact failure. Before the table, in plain words."""
    text = _verdict_text(_drivers(("Escherichia coli", 4, False, None)))
    assert "No organism behind this reading was found in your sample" in text
    assert "trace signal with no identified source" in text


def test_a_reading_with_organisms_present_names_them_and_their_share() -> None:
    text = _verdict_text(_drivers(
        ("Collinsella aerofaciens", 60, True, 0.85),
        ("Escherichia coli", 40, False, None),
    ))
    assert "Collinsella aerofaciens" in text
    assert "0.85% of your community" in text
    assert "60%" in text, "the share carried by organisms actually present"


def test_an_absent_organism_is_never_described_as_a_cause() -> None:
    """A resemblance is not a presence, and the wording has to keep them apart."""
    text = _verdict_text(_drivers(
        ("Present bug", 50, True, 1.0), ("Absent bug", 50, False, None),
    ))
    assert "resemblance rather than a presence" in text


def test_the_organism_evidence_block_lists_only_organisms_in_the_sample() -> None:
    """Advice about a species somebody does not carry is advice for somebody else."""
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota.pdfreport import _organism_evidence_block, _styles

    story: list[Any] = []
    _organism_evidence_block(story, _styles(), _drivers(
        ("Collinsella aerofaciens", 60, True, 0.85),
        ("Escherichia coli", 40, False, None),
    ))
    if not story:
        pytest.skip("no curated record for the fixture organisms in this checkout")
    SimpleDocTemplate("/tmp/_orgev.pdf", pagesize=A4).build(story)
    text = " ".join(
        "".join(p.extract_text() for p in pypdf.PdfReader("/tmp/_orgev.pdf").pages).split()
    )
    assert "Collinsella aerofaciens" in text
    assert "Escherichia coli" not in text, "an organism not in the sample was described"


def test_association_only_evidence_is_labelled_as_such() -> None:
    """So a reader cannot mistake a correlation for a lever."""
    import pypdf
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate

    from openbiota.pdfreport import _organism_evidence_block, _styles

    story: list[Any] = []
    _organism_evidence_block(story, _styles(),
                             _drivers(("Collinsella aerofaciens", 10, True, 0.85)))
    if not story:
        pytest.skip("no curated record for Collinsella in this checkout")
    SimpleDocTemplate("/tmp/_assoc.pdf", pagesize=A4).build(story)
    text = " ".join(
        "".join(p.extract_text() for p in pypdf.PdfReader("/tmp/_assoc.pdf").pages).split()
    )
    assert "Association only" in text
