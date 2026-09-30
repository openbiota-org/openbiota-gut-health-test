"""The merged organism inventory: what may be pooled, and what may not.

One rule runs through this module and every test here exists to hold it:
detection pools every catalogue, ranking does not. "Is this organism here"
is answered by the sample, so the widest search wins. "How does this compare
with other people" is answered against a reference population, and that
population was measured on one catalogue, so a reading from a different one
cannot be ranked against it without inventing the number.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from openbiota import inventory
from openbiota.engines.metaphlan4 import RENAMED

RESULTS = sorted(Path("results").glob("*/results.json"))


def _lane(name: str, species: dict[str, float], *, rankable: bool) -> inventory.Lane:
    return inventory.Lane(name=name, species=species, rankable=rankable)


class _Row:
    """Stand-in for a reference-scored row."""

    def __init__(self, species: str, percent: float, percentile: float) -> None:
        self.species = species
        self.percent = percent
        self.percentile = percentile
        self.cohort_prevalence = 0.5
        self.groups: tuple[str, ...] = ()
        self.trace = False
        self.rare_in_cohort = False


def test_detection_pools_every_catalogue() -> None:
    """An organism either catalogue names is in the inventory."""
    inv = inventory.build([
        _lane("a", {"Alpha_one": 1.0, "Shared_two": 2.0}, rankable=True),
        inventory.Lane(name="b", species={"Shared_two": 2.5, "Beta_three": 0.5},
                       rankable=False, primary=True),
    ])
    assert {o.species for o in inv} == {"Alpha_one", "Shared_two", "Beta_three"}


def test_abundances_come_from_one_lane_and_are_never_merged() -> None:
    """Two lanes are two compositions with two denominators.

    Taking the larger reading per organism let a real inventory total 130%.
    The primary lane's value is the share; the other lane's is kept beside
    it as a secondary reading and never enters the total.
    """
    inv = inventory.build([
        _lane("a", {"Alpha_one": 1.0, "Shared_two": 2.0}, rankable=True),
        inventory.Lane(name="b", species={"Shared_two": 2.5, "Beta_three": 0.5},
                       rankable=False, primary=True),
    ])
    shared = inv.get("Shared_two")
    assert shared.in_primary and shared.percent == pytest.approx(2.5)
    assert shared.secondary_percent == pytest.approx(2.0)
    alpha = inv.get("Alpha_one")
    assert not alpha.in_primary and alpha.percent == 0.0
    assert alpha.secondary_percent == pytest.approx(1.0)
    # The composition is the primary lane's total and nothing else.
    assert inv.composition_total == pytest.approx(3.0)
    assert inv.primary_lane == "b"


def test_two_primary_lanes_is_a_configuration_error() -> None:
    with pytest.raises(ValueError, match="exactly one lane may be primary"):
        inventory.build([
            inventory.Lane(name="a", species={"X_y": 1.0}, rankable=False, primary=True),
            inventory.Lane(name="b", species={"X_y": 1.0}, rankable=False, primary=True),
        ])


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_no_real_composition_exceeds_one_hundred() -> None:
    for path in RESULTS:
        inv = inventory.from_json(json.loads(path.read_text()).get("organism_inventory"))
        if inv is None:
            continue
        # estimated shares are drawn from the unclassified band, so the whole
        # is composition + what remains unplaced, never composition + the
        # full band
        remainder = inv.unplaced_percent if inv.unplaced_percent is not None else (inv.unclassified_percent or 0.0)
        total = inv.composition_total + remainder
        assert inv.composition_total <= 100.0 + 1e-6, (path, inv.composition_total)
        assert total <= 100.5, (path, total)


def test_a_rename_is_one_organism_not_two() -> None:
    """The same organism under two names must not be counted twice."""
    old, new = next(iter((k, v) for k, v in RENAMED.items() if k != v))
    inv = inventory.build([
        _lane("a", {old: 1.0}, rankable=True),
        _lane("b", {new: 1.2}, rankable=False),
    ])
    assert len(inv) == 1, [o.species for o in inv]
    only = inv.organisms[0]
    assert only.species == new
    assert only.formerly == old
    # Findable under either name, because a reader may know the older one.
    assert inv.get(old) is only
    assert inv.get(new) is only


def test_no_organism_claims_a_former_name_identical_to_its_own() -> None:
    """The source map records checked-but-unchanged names; those are not renames."""
    for old, new in RENAMED.items():
        if old == new:
            assert inventory.former_name(new) is None, new


def test_an_unrankable_detection_never_acquires_a_percentile() -> None:
    """The line that must not be crossed.

    An organism measured only on a catalogue with no reference population
    gets no percentile, even when a row of the same name happens to exist.
    Ranking it against a cohort measured another way would produce a number
    that looks authoritative and means nothing.
    """
    inv = inventory.build(
        [_lane("unrankable", {"Ghost_species": 3.0}, rankable=False)],
        ranked_rows=[_Row("Ghost_species", 3.0, 88.0)],
    )
    only = inv.organisms[0]
    assert only.percentile is None
    assert only.rankable is False


def test_a_rankable_detection_keeps_the_abundance_its_rank_came_from() -> None:
    """A row must not mix one catalogue's amount with another's percentile."""
    inv = inventory.build(
        [
            _lane("scoring", {"Real_species": 2.0}, rankable=True),
            _lane("wider", {"Real_species": 9.0}, rankable=False),
        ],
        ranked_rows=[_Row("Real_species", 2.0, 70.0)],
    )
    only = inv.organisms[0]
    assert only.percentile == 70.0
    assert only.percent == pytest.approx(2.0), (
        "the displayed amount must be the one the percentile was computed "
        "from, or the row contradicts itself"
    )


def test_strain_fingerprints_join_by_genome_bin() -> None:
    inv = inventory.build(
        [inventory.Lane(
            name="wider", species={}, rankable=False,
            rows=[{"species": "Some_organism", "percent": 1.0,
                   "genus": "Some", "sgb": "SGB1"}],
        )],
        strain={"organisms": [{"sgb": "SGB1", "markers_resolved": 42}]},
    )
    assert inv.organisms[0].resolved_to_strain
    assert inv.organisms[0].strain["markers_resolved"] == 42


def test_the_inventory_survives_the_results_file() -> None:
    inv = inventory.build([
        _lane("a", {"Alpha_one": 1.0}, rankable=True),
        _lane("b", {"Beta_two": 0.5}, rankable=False),
    ])
    back = inventory.from_json(inv.to_json())
    assert back is not None
    assert back.counts() == inv.counts()


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_real_samples_see_more_than_one_catalogue_alone() -> None:
    """The point of the merge, checked against real output."""
    for path in RESULTS:
        results = json.loads(path.read_text())
        blob = results.get("organism_inventory")
        if not blob or blob.get("status") == "stage_error":
            continue
        inv = inventory.from_json(blob)
        assert inv is not None
        # Every rankable organism has a percentile and vice versa.
        for o in inv:
            assert (o.percentile is not None) == o.rankable
        # Pooling beat the single catalogue it is compared against.
        single = len((results.get("profile_similarity") or {})
                     .get("community", {}).get("species") or [])
        if single:
            assert len(inv) >= single, (
                f"{path}: merged inventory ({len(inv)}) is smaller than one "
                f"catalogue alone ({single}), which cannot happen if the "
                "merge is a union"
            )


# --------------------------------------------------------------------------- #
# Prose hygiene
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_no_internal_token_or_shell_command_reaches_the_reader() -> None:
    """Nothing in a rendered report may read as a defect.

    Two classes of leak have reached the page before: a bare status
    identifier dropped into the middle of a sentence
    (``(unmatched_context_only)``), and an instruction to run a command
    (``run `openbiota calibrate-detection` ``). A person reading their own
    results can act on neither. Both are caught here rather than by eye.
    """
    import re

    import pymupdf

    banned = re.compile(
        r"unmatched_context_only|comparator_unavailable|uncalibrated prior"
        r"|openbiota [a-z-]+`|stage_error|not_computable|no_validated_panel"
        r"|Traceback|NoneType|\bTODO\b|\bFIXME\b",
    )
    for path in RESULTS:
        pdf = next(path.parent.glob("*_report.pdf"), None)
        if pdf is None:
            continue
        with pymupdf.open(pdf) as doc:
            text = "".join(page.get_text() for page in doc)
        hits = sorted({m.group(0) for m in banned.finditer(text)})
        assert not hits, f"{pdf.name} shows internal text to the reader: {hits}"


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_every_headline_organism_count_is_the_same_number() -> None:
    """One number for "how many organisms", everywhere it is printed.

    Page one said 251 and the community page it linked to said 92, because
    the community page still counted one catalogue. Every place that prints
    a detection total must print the merged inventory's total. The
    reference-comparable richness is a different quantity and is labelled as
    such, so it is allowed to differ - and is checked separately below.
    """
    import re

    import pymupdf

    for path in RESULTS:
        pdf = next(path.parent.glob("*_report.pdf"), None)
        inv = json.loads(path.read_text()).get("organism_inventory") or {}
        n = inv.get("n_organisms")
        if pdf is None or not n:
            continue
        with pymupdf.open(pdf) as doc:
            p1 = doc[0].get_text()
            texts = [pg.get_text() for pg in doc]
            # Locate sections by the outline, not by matching their title in
            # the page text: page 3's contents table lists every section by
            # number and name, and reads exactly like a heading.
            toc = {
                int(title.split()[0]): page
                for _lvl, title, page in doc.get_toc()
                if title.split() and title.split()[0].isdigit()
            }
        # Page one: the donut centre and the tile.
        assert re.search(rf"\b{n}\b\s*\n?\s*ORGANISMS", p1), f"{pdf.name}: page-one donut is not {n}"
        assert re.search(rf"\b{n}\b\s*\n?\s*ALL ORGANISMS DETECTED", p1), f"{pdf.name}: page-one tile is not {n}"
        # The community section: the tile it links to. Located by the
        # bookmark, which is what the section heading itself registers, and
        # numbered from SECTIONS rather than written in - the numbers shift
        # when a section is added or, as here, unnumbered.
        from openbiota.pdfreport import SECTIONS

        community = SECTIONS["community"]
        assert community in toc, f"{pdf.name}: the community section is not in the outline"
        sec4 = texts[toc[community] - 1]
        assert "Your gut community" in sec4, (
            f"{pdf.name}: the section {community} bookmark lands elsewhere"
        )
        m = re.search(r"\b(\d+)\s*\n\s*ORGANISMS DETECTED", sec4)
        assert m, f"{pdf.name}: the community section has no 'Organisms detected' tile"
        assert int(m.group(1)) == n, (
            f"{pdf.name}: the community tile says {m.group(1)}, inventory is {n}"
        )
        # The classified-organisms section: the heading count.
        sec7 = next((t for t in texts if "Your organisms, classified" in t), "")
        assert re.search(rf"\b{n} organisms\b", sec7), f"{pdf.name}: the classified section does not show {n}"
        # Nowhere may the old single-lane phrase survive as a headline.
        for i, t in enumerate(texts):
            assert not re.search(r"\b\d+\s+species the profiler named", t), (pdf.name, i + 1)
        # The most abundant organism on the community page is the most
        # abundant organism in the composition, under the same name and
        # number. It once read "Bacteroides vulgatus 22.3%" there and
        # "Phocaeicola vulgatus 12.01%" four pages on.
        top = next((o for o in (inv.get("organisms") or []) if o.get("in_primary")), None)
        if top:
            name = (top.get("gtdb") if top.get("unnamed") and top.get("gtdb") else top["species"]).replace("_", " ")
            m4 = re.search(r"Most abundant organisms\s*\n(.*?)\n(\d+\.\d)%", sec4, re.S)
            assert m4, f"{pdf.name}: no 'Most abundant organisms' list on the community section"
            assert name.split()[0] in m4.group(1), (
                f"{pdf.name}: the community section leads with {m4.group(1).strip()[:40]!r}, composition leads with {name!r}")
            assert abs(float(m4.group(2)) - float(top["percent"])) < 0.15, (
                f"{pdf.name}: the community section says {m4.group(2)}% for {name}, composition says {top['percent']}%")
        # One Firmicutes:Bacteroidetes ratio per page.
        ratios = {m.group(1) for m in re.finditer(r"(\d+\.\d\d)\s*\n\s*FIRMICUTES : BACTEROIDETES", sec4)}
        ratios |= {m.group(1) for m in re.finditer(r"Firmicutes : Bacteroidetes\s*\n\s*(\d+\.\d\d)", sec4)}
        assert len(ratios) <= 1, f"{pdf.name}: section 4 prints two different F:B ratios: {sorted(ratios)}"


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_reference_comparable_richness_is_labelled_not_disguised() -> None:
    """The smaller count is legitimate - it is what the percentile is valid
    on - but it must say so, never masquerade as the total."""
    import pymupdf

    for path in RESULTS:
        pdf = next(path.parent.glob("*_report.pdf"), None)
        if pdf is None:
            continue
        with pymupdf.open(pdf) as doc:
            text = "".join(pg.get_text() for pg in doc)
        if "Richness, reference-comparable" in text:
            assert "It is not the total found" in text, pdf.name


@pytest.mark.skipif(not RESULTS, reason="no analysed samples on disk")
def test_page_two_carries_only_findings_with_graphics() -> None:
    """Page 2 is a grid of findings, each beside the graphic it came from.

    A caveat has no graphic. Placed on this page it sits beside an empty
    cell and reads as a rendering defect - which is exactly what happened
    when the census's own qualifications were appended as rows. Those
    sentences belong in the fine print of the section they qualify, and
    this pins that they never come back to page 2.
    """
    import pymupdf

    caveats = (
        "No age, sex or country",
        "How reliably a low-abundance organism",
        "conservative general curve",
        "The reference group is adults and this sample",
    )
    for path in RESULTS:
        pdf = next(path.parent.glob("*_report.pdf"), None)
        if pdf is None:
            continue
        with pymupdf.open(pdf) as doc:
            if doc.page_count < 2:
                continue
            page_two = doc[1].get_text()
        leaked = [c for c in caveats if c in page_two]
        assert not leaked, (
            f"{pdf.name} page 2 renders a caveat as a findings row: {leaked}"
        )
