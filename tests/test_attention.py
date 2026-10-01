"""The organisms-that-need-attention list: one list, two sections, no gaps.

Section 6 (the glance) and section 16 (the cards) were built from two
different systems and disagreed by twenty-two organisms on one sample. They
are now drawn from one list built by :func:`openbiota.pdfattention.build`.
These tests hold that invariant, and the one that bit first: an organism
known to the census under its catalogue name and to the inventory under its
current name must be listed once, not twice.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pytest

from openbiota import inventory, modulators, organisms, pdfattention

RESULTS = sorted(Path("results").glob("*/results.json"))
REPORTS = sorted(Path("results").glob("*/*_report.pdf"))

pytestmark = pytest.mark.skipif(not RESULTS, reason="no results/ to check against")


@dataclass
class _Taxon:
    species: str
    display_name: str
    percent: float
    percentile: float | None
    prevalence: float | None
    detected: bool
    detection_power: float | None
    lod95_percent: float | None
    bucket: int | None


@dataclass
class _Ref:
    n_samples: int
    n_studies: int


class _Findings:
    """Just enough of FindingsResult for build(): bucket() and reference."""

    def __init__(self, blob: dict) -> None:
        self.taxa = []
        for t in blob["findings"]["taxa"]:
            det = t.get("detection") or {}
            pos = t.get("reference_position") or {}
            prev = t.get("reference_prevalence") or {}
            self.taxa.append(_Taxon(
                species=t["species"], display_name=t.get("display_name") or t["species"].replace("_", " "),
                percent=float(t.get("percent") or 0.0), percentile=pos.get("percentile_among_carriers"),
                prevalence=prev.get("study_balanced"), detected=det.get("status") == "detected",
                detection_power=det.get("power_if_present"), lod95_percent=det.get("lod95_percent"),
                bucket=t.get("bucket"),
            ))
        ref = blob["findings"].get("reference") or {}
        self.reference = _Ref(int(ref.get("n_samples") or 0), int(ref.get("n_studies") or 0))

    def bucket(self, n: int) -> list[_Taxon]:
        return [t for t in self.taxa if t.bucket == n]


def _load(path: Path):
    blob = json.loads(path.read_text())
    inv = inventory.from_json(blob.get("organism_inventory"))
    findings = _Findings(blob["findings_and_evidence"]) if blob.get("findings_and_evidence") else None
    return inv, findings


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_no_organism_is_listed_twice_under_two_names(path: Path) -> None:
    """Ruminococcus_gnavus (census) and Mediterraneibacter_gnavus (inventory)
    are one organism and get one row."""
    inv, findings = _load(path)
    items = pdfattention.build(inv, findings)
    names = [a.species for a in items]
    assert len(names) == len(set(names)), "same species twice"
    # Collapse every alias to the inventory's current name and check again.
    canon = []
    for a in items:
        o = inv.get(a.species) if inv is not None else None
        canon.append(o.species if o is not None else a.species)
    dupes = {c for c in canon if canon.count(c) > 1}
    assert not dupes, f"one organism under two names: {sorted(dupes)}"


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_every_flagged_verdict_and_every_missing_taxon_is_in_the_list(path: Path) -> None:
    """The page is exactly the image of the one rule: every organism the
    verdict flags is listed, nothing listed as present is unflagged, and the
    census adds only what no method found."""
    inv, findings = _load(path)
    items = pdfattention.build(inv, findings)
    listed = {a.species for a in items}
    flagged = {v.organism.species for v in organisms.verdicts(list(inv.organisms)) if v.flagged}
    assert flagged <= listed, f"flagged but not listed: {sorted(flagged - listed)}"
    present_listed = {a.species for a in items if a.detected}
    assert present_listed == flagged, (
        f"listed as present without a flag: {sorted(present_listed - flagged)}; "
        f"flagged but missing: {sorted(flagged - present_listed)}")
    assert len(listed) == len(items), "an organism is listed twice"
    if findings is not None:
        for t in findings.bucket(4):
            o = inv.get(t.species)
            # an organism any method found with support is not missing, whatever the reference catalogue saw
            found = o is not None and ((o.in_primary and o.percent > 0)
                                       or (o.status == "supported" and (o.best_percent or 0.0) > 0))
            assert t.species in listed or (o is not None and o.species in listed) or found, \
                f"expected-but-absent taxon not listed: {t.species}"
        for a in items:
            if a.group != "missing":
                continue
            o = inv.get(a.species)
            assert o is None or not (o.in_primary and o.percent > 0), f"listed missing but found: {a.species}"
            assert o is None or o.status != "supported" or not (o.best_percent or 0.0), \
                f"listed missing but supported: {a.species}"


@pytest.mark.parametrize("path", RESULTS, ids=lambda p: p.parent.name)
def test_every_listed_organism_has_a_role_and_something_that_moves_it(path: Path) -> None:
    """The user's complaint was cards that said nothing about what the
    organism does or what changes it. Every card must have both."""
    inv, findings = _load(path)
    items = pdfattention.build(inv, findings)
    no_role = [a.species for a in items if a.modulators is None or not a.modulators.role]
    no_lever = [a.species for a in items
                if a.want is not None
                and (a.modulators is None or not a.modulators.for_direction(a.want))]
    assert not no_role, f"no description: {no_role}"
    # Not every organism has directional evidence - Dielma has none - but
    # the ones the reader will act on must.
    actionable = [s for s in no_lever if any(a.species == s and a.group in ("overgrown", "missing")
                                               for a in items)]
    allowed = {"Dielma_fastidiosa"}
    assert set(actionable) <= allowed, f"overgrown/missing with nothing to do: {actionable}"


def test_every_modulator_item_states_its_evidence_level_and_transit() -> None:
    table = modulators._table()
    assert table, "registry did not load"
    for m in set(table.values()):
        for ag in m.decrease + m.increase:
            assert ag.evidence in modulators.EVIDENCE, f"{m.organism}: unknown evidence {ag.evidence!r}"
            assert ag.kind in modulators.KIND, f"{m.organism}: unknown kind {ag.kind!r}"
            assert ag.source, f"{m.organism}: {ag.agent!r} has no source"
            assert ag.survives_transit, f"{m.organism}: {ag.agent!r} does not reach the gut; do not list it"


def test_genus_fallback_is_marked_as_genus_level() -> None:
    m = modulators.lookup("Blautia_wexlerae")
    assert m is not None and m.basis == "genus"
    s = modulators.lookup("Mediterraneibacter_gnavus")
    assert s is not None and s.basis == "species"
    # Old name resolves to the same species entry.
    assert modulators.lookup("Ruminococcus_gnavus") is not None
    assert modulators.lookup("Ruminococcus_gnavus").organism == s.organism
    assert modulators.lookup("Nothing_at_all") is None


@pytest.mark.parametrize("path", REPORTS, ids=lambda p: p.parent.name)
def test_rendered_cards_match_the_list(path: Path) -> None:
    """In the rendered PDF there is exactly one card per organism in the list
    - so the glance, which is drawn from the same list, cannot disagree with
    the cards, and no card is truncated away."""
    pymupdf = pytest.importorskip("pymupdf")
    inv, findings = _load(path.parent / "results.json")
    expected = len(pdfattention.build(inv, findings))
    doc = pymupdf.open(path)
    text = "".join(page.get_text() for page in doc)
    cards = text.count("What it is and what it does")
    assert cards == expected, f"{cards} cards rendered for {expected} organisms needing attention"
    assert text.count("Why it is flagged in your sample") == expected
