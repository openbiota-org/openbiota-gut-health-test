"""Every problem this report found, and everything known about changing it.

The actions section used to be five blocks that each knew about part of the
report: a lever table that only saw flagged organisms, an evidence overview
that only saw registry mappings, a planner, a model-simulation block and a
list of capabilities the run could not reach. A reader had to assemble the
answer themselves, and most findings were not in any of the blocks at all.

This module assembles it once. It takes every finding the report makes -
organisms that are overgrown, depleted or missing, metabolite pathways that
sit outside the usual range, disease-pattern resemblances, detected
pathogens, biofilm and barrier readings - and turns each one into a
:class:`Problem` with its own options attached. A problem with no option
says so in those words; it does not quietly disappear.

Options come from three registries and are treated identically once here:

* ``taxa/classification/modulators.yaml`` - what moves a named organism;
* ``interventions/*.yaml`` - the exact-protocol evidence cards, with their
  safety rules and counter-evidence;
* ``interventions/atlas_options.yaml`` - the research atlas, imported by
  ``scripts/import_atlas_evidence.py``.

No option is dropped for the design of its study. A result from a dish and a
randomised trial are both real findings about the organism; what separates
them is a word and a small marker, decided by :mod:`openbiota.delivery`, and
the fine print that goes at the foot of the table.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota import delivery as DEL

ATLAS_FILE: Final = Path("interventions/atlas_options.yaml")

#: What a problem wants done about it.
WANT_DOWN: Final = "decrease"
WANT_UP: Final = "increase"

#: Reader-facing words for each material kind, and the order they are shown
#: in: things you eat first, then things you buy, then things a prescriber
#: decides. Nothing is hidden behind a category - this is only ordering.
KIND_WORDS: Final[dict[str, tuple[str, int]]] = {
    "food": ("food", 0),
    "diet": ("way of eating", 1),
    "fibre": ("fibre", 2),
    "prebiotic": ("prebiotic fibre", 3),
    "probiotic": ("probiotic", 4),
    "spore_probiotic": ("spore probiotic", 5),
    "postbiotic": ("postbiotic", 6),
    "polyphenol": ("plant compound", 7),
    "herb": ("herb or botanical extract", 8),
    "essential_oil": ("essential oil", 9),
    "enzyme": ("enzyme", 10),
    "mineral": ("mineral or protein supplement", 11),
    "compound": ("supplement or compound", 12),
    "lifestyle": ("lifestyle", 13),
    "drug": ("medication - a prescriber's decision", 14),
    "procedure": ("procedure", 15),
    "clinical_route": ("clinical route", 16),
}

#: Kinds a reader cannot act on alone. They are shown, and marked.
PRESCRIBER = frozenset({"drug", "procedure", "clinical_route"})


# --------------------------------------------------------------------------- #
# names
# --------------------------------------------------------------------------- #

#: A GTDB rank suffix: `Blautia_A luti`, `Ruminococcus_B gnavus`. It has to
#: come off *before* underscores become spaces, or there is no underscore
#: left to recognise it by and the suffix survives as a word of its own.
_SUFFIX = re.compile(r"_[A-Z]{1,2}(?=[\s_]|$)")
_PUNCT = re.compile(r"[^a-z0-9 ]")
_WS = re.compile(r"\s+")


def norm(name: str | None) -> str:
    """One spelling for an organism, whatever registry it came from.

    GTDB writes ``Blautia_A luti`` and the literature writes ``Blautia
    luti``; the report itself writes ``Mediterraneibacter_gnavus`` where the
    research writes ``Ruminococcus gnavus``. Underscores, suffixes and case
    all go, so the three meet.
    """
    s = _SUFFIX.sub("", str(name or "")).replace("_", " ")
    return _WS.sub(" ", _PUNCT.sub("", s.lower())).strip()


def organism_keys(row: Mapping[str, Any]) -> set[str]:
    """Every name one inventory row could be matched by."""
    keys: set[str] = set()
    for value in (row.get("species"), row.get("gtdb"), row.get("formerly"),
                  (row.get("native_ids") or {}).get("scoring")):
        if value:
            keys.add(norm(value))
    for alias in (row.get("aliases") or ()):
        keys.add(norm(alias))
    for alias in (row.get("formerly_listed_as") or ()):
        if isinstance(alias, str):
            keys.add(norm(alias))
    return {k for k in keys if k}


def genus_of(name: str) -> str:
    return norm(name).split(" ")[0] if name else ""


# --------------------------------------------------------------------------- #
# options
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Option:
    """One thing that has been shown to move something, with its setting."""

    label: str
    kind: str
    direction: str
    setting: str
    delivery: str
    effect: str
    boundary: str
    source: str
    option_id: str
    #: "species" when the study named this organism, "genus" when it named
    #: only its genus, "panel"/"concern" for a pathway or a general finding.
    binding: str

    @property
    def reach(self) -> DEL.Reach:
        return DEL.assess(self.setting, self.delivery)

    @property
    def kind_words(self) -> str:
        return KIND_WORDS.get(self.kind, (self.kind.replace("_", " "), 99))[0]

    @property
    def kind_rank(self) -> int:
        return KIND_WORDS.get(self.kind, ("", 99))[1]

    @property
    def needs_prescriber(self) -> bool:
        return self.kind in PRESCRIBER

    @property
    def sort_key(self) -> tuple[Any, ...]:
        """Strongest evidence first, then the kind, then stably by name."""
        binding_rank = {"strain": -1, "species": 0, "panel": 1, "concern": 2, "genus": 3}
        return (
            binding_rank.get(self.binding, 4),
            *DEL.sort_key(self.setting, self.delivery),
            self.kind_rank,
            self.label.lower(),
        )


@lru_cache(maxsize=1)
def _atlas() -> tuple[dict[str, Any], ...]:
    """The imported atlas options, or nothing when the file is not built."""
    path = ATLAS_FILE if ATLAS_FILE.is_file() else Path(__file__).resolve().parents[1] / ATLAS_FILE
    if not path.is_file():
        return ()
    try:
        loaded = yaml.safe_load(path.read_text()) or {}
    except (OSError, yaml.YAMLError):
        return ()
    return tuple(loaded.get("atlas_options") or ())


@lru_cache(maxsize=1)
def _index() -> dict[str, dict[str, list[dict[str, Any]]]]:
    """The atlas, indexed by what each option is about."""
    by: dict[str, dict[str, list[dict[str, Any]]]] = {
        "species": {}, "genus": {}, "panel": {}, "concern": {},
    }
    for option in _atlas():
        for taxon in option.get("taxa") or ():
            by["species"].setdefault(norm(taxon), []).append(option)
            by["genus"].setdefault(genus_of(taxon), []).append(option)
        for genus in option.get("genera") or ():
            by["genus"].setdefault(norm(genus), []).append(option)
        for panel in option.get("panels") or ():
            by["panel"].setdefault(str(panel), []).append(option)
        for concern in option.get("concerns") or ():
            by["concern"].setdefault(str(concern), []).append(option)
    return by


def _as_option(raw: Mapping[str, Any], *, binding: str, direction: str | None = None) -> Option:
    source = str(raw.get("url") or "")
    if not source and raw.get("pmid"):
        source = f"https://pubmed.ncbi.nlm.nih.gov/{raw['pmid']}/"
    if not source and raw.get("doi"):
        source = f"https://doi.org/{raw['doi']}"
    return Option(
        label=str(raw.get("material") or "").strip(),
        kind=str(raw.get("kind") or "compound"),
        direction=str(direction or raw.get("direction") or "unclear"),
        setting=str(raw.get("setting") or "in_vitro"),
        delivery=str(raw.get("delivery") or "unknown"),
        effect=str(raw.get("effect") or "").strip(),
        boundary=str(raw.get("boundary") or "").strip(),
        source=source,
        option_id=str(raw.get("option_id") or ""),
        binding=binding,
    )


def atlas_options(
    *,
    species: Iterable[str] = (),
    genera: Iterable[str] = (),
    panels: Iterable[str] = (),
    concerns: Iterable[str] = (),
    want: str | None = None,
) -> list[Option]:
    """Everything the atlas has about these targets, in the wanted direction.

    ``want`` keeps only the options that move the finding the way the reader
    needs. This is strict on purpose. An earlier version also kept rows whose
    direction the source never stated, and the result was that fibres shown
    to *feed* an overgrown organism were printed under what would reduce it.
    A row that does not say which way it went is not a way to fix anything,
    so it is not offered as one.
    """
    index = _index()
    found: dict[str, Option] = {}
    for binding, keys in (("species", species), ("panel", panels),
                          ("concern", concerns), ("genus", genera)):
        for key in keys:
            lookup = norm(key) if binding in {"species", "genus"} else str(key)
            for raw in index[binding].get(lookup, ()):
                # A row that names several organisms carries a direction for
                # each of them; use the one for the organism actually asked
                # about rather than the row's overall direction.
                direction = option_direction(raw, key if binding == "species" else None)
                option = _as_option(raw, binding=binding, direction=direction)
                if not option.label or not option.reach.reportable:
                    continue
                if want and option.direction != want:
                    continue
                # One row can be reached both as "this organism" and as "its
                # genus". The species reading is the one that carries the
                # direction for the organism actually asked about, so it wins
                # - otherwise a fibre shown to feed an organism arrives under
                # both what lowers it and what raises it.
                rank = {"species": 0, "panel": 1, "concern": 2, "genus": 3}
                seen = found.get(option.option_id)
                if seen is None or rank[binding] < rank[seen.binding]:
                    found[option.option_id] = option
    return sorted(found.values(), key=lambda o: o.sort_key)


def option_direction(raw: Mapping[str, Any], taxon: str | None) -> str:
    """The direction this row reports *for this organism*, if it says."""
    if taxon:
        wanted = norm(taxon)
        for name, direction in (raw.get("taxon_direction") or {}).items():
            if norm(name) == wanted:
                return str(direction)
    return str(raw.get("direction") or "unclear")


def _opposite(want: str) -> str:
    return WANT_UP if want == WANT_DOWN else WANT_DOWN


# --------------------------------------------------------------------------- #
# problems
# --------------------------------------------------------------------------- #


def _dest(kind: str, name: str) -> str:
    """Where this report explains this particular reading.

    "See section 16" is not an answer when section 16 is twenty pages of
    organisms: the reader has to land on the organism. Every reading already
    has its own destination for exactly this reason, so the actions section
    uses it rather than the top of the section.
    """
    from openbiota import pdflinks  # noqa: PLC0415 - pdflinks imports this module's siblings

    if not name:
        return ""
    if kind == "organism":
        return pdflinks.organism_dest(name, detail=True)
    if kind == "panel":
        return pdflinks.function_dest(name, detail=True)
    if kind == "pathogen":
        return pdflinks.pathogen_dest(name, detail=True)
    if kind == "pattern":
        return pdflinks.profile_dest(name, detail=True)
    return ""


def problem_key(kind: str, name: str) -> str:
    """The identity of one finding, the same from either end.

    A detail page and the actions section both build this from the organism
    or pathway they are about, so the link between them needs nothing to be
    handed across.
    """
    return f"{kind}:{norm(name)}"


@dataclass
class Problem:
    """One finding of this report, and what is known about changing it."""

    key: str
    kind: str
    title: str
    detail: str
    want: str
    severity: float
    section: int
    #: The destination of this reading's own page in the report, so "where"
    #: lands on the finding rather than at the top of a section the reader
    #: then has to scroll.
    detail_dest: str = ""
    options: list[Option] = field(default_factory=list)
    #: Things shown to push this finding further the wrong way. Knowing what
    #: feeds an overgrowth is as useful as knowing what starves it, and it
    #: keeps those rows from being mistaken for solutions.
    worsens: list[Option] = field(default_factory=list)
    #: Set when nothing is on record, so the page can say so in words.
    nothing_known: bool = False
    #: Why there is nothing, when the registry has an explanation for it.
    notes: list[str] = field(default_factory=list)

    @property
    def markers(self) -> set[str]:
        return {m for o in (*self.options, *self.worsens) for m in o.reach.markers}

    def top(self, n: int) -> list[Option]:
        return self.options[:n]


def _detail_for(item: Any) -> str:
    """The number that makes an organism finding concrete."""
    bits: list[str] = []
    percent = getattr(item, "percent", None)
    if isinstance(percent, (int, float)) and percent > 0:
        bits.append(f"{percent:.2f}% of your community" if percent < 1 else f"{percent:.1f}% of your community")
    deviation = getattr(item, "deviation", None)
    if isinstance(deviation, (int, float)):
        if deviation >= 100:
            bits.append(f"{1 + deviation / 100:.0f}\u00d7 the typical carrier")
        elif deviation > 0:
            bits.append(f"{deviation:.0f}% above the typical carrier")
        elif deviation < 0:
            bits.append(f"{abs(deviation):.0f}% below the typical carrier")
    percentile = getattr(item, "percentile", None)
    if isinstance(percentile, (int, float)):
        bits.append(f"{percentile:.0f}th percentile")
    return "; ".join(bits)


#: How a flagged organism's group maps onto what the reader wants to happen.
_GROUP_WANT: Final[dict[str, str]] = {
    "overgrown": WANT_DOWN,
    "missing": WANT_UP,
    "depleted": WANT_UP,
    "watch": WANT_DOWN,
}

_GROUP_TITLE: Final[dict[str, str]] = {
    "overgrown": "is overgrown",
    "missing": "was not detected",
    "depleted": "is well below the usual level",
    "watch": "needs qualification",
}


def alias_map(inventory: Mapping[str, Any] | None) -> dict[str, set[str]]:
    """Every name each detected organism also answers to.

    This is what lets the research reach the report at all. The report speaks
    GTDB r232 and the literature speaks NCBI, so the organism this report
    calls *Mediterraneibacter gnavus* is the one forty years of papers call
    *Ruminococcus gnavus*, and without the crosswalk not one of those papers
    would be found. The inventory already carries the old name in
    ``formerly`` and the per-method names in ``native_ids``; this turns them
    into a lookup.
    """
    out: dict[str, set[str]] = {}
    for row in ((inventory or {}).get("organisms") or ()):
        keys = organism_keys(row)
        for key in keys:
            out.setdefault(key, set()).update(keys)
    return out


def organism_problems(
    items: Sequence[Any], *, section: int, inventory: Mapping[str, Any] | None = None
) -> list[Problem]:
    """A problem for every organism the report flagged."""
    aliases = alias_map(inventory)
    out: list[Problem] = []
    for item in items:
        group = str(getattr(item, "group", "") or "watch")
        want = _GROUP_WANT.get(group, WANT_DOWN)
        display = str(getattr(item, "display", "") or getattr(item, "species", ""))
        name = norm(getattr(item, "species", "") or display)
        names = set(aliases.get(name) or set()) | {name, norm(display)}
        species = str(getattr(item, "species", "") or display)
        problem = Problem(
            key=problem_key("organism", species),
            kind="organism",
            title=f"{display} {_GROUP_TITLE.get(group, 'needs attention')}",
            detail=_detail_for(item),
            want=want,
            severity=float(getattr(item, "weight", 0.0) or 0.0),
            section=section,
            detail_dest=_dest("organism", species),
        )
        genera = sorted({genus_of(n) for n in names if n})
        problem.options = atlas_options(species=sorted(names), genera=genera, want=want)
        problem.options += _modulator_options(getattr(item, "modulators", None), want=want)
        problem.options = _dedupe(problem.options)
        problem.worsens = _dedupe(
            atlas_options(species=sorted(names), genera=genera, want=_opposite(want))
            + _modulator_options(getattr(item, "modulators", None), want=_opposite(want))
        )
        problem.nothing_known = not problem.options
        out.append(problem)
    return out


def _modulator_options(modulators: Any, *, want: str) -> list[Option]:
    """The organism registry's own agents, in the same shape as the atlas."""
    if modulators is None:
        return []
    side = getattr(modulators, "decrease" if want == WANT_DOWN else "increase", ()) or ()
    out: list[Option] = []
    for agent in side:
        kind = str(getattr(agent, "kind", "") or "compound")
        setting = DEL.LEGACY_SETTING.get(str(getattr(agent, "evidence", "")), "in_vitro")
        stated = "oral_tested" if getattr(agent, "survives_transit", True) else "enteric_required"
        option = Option(
            label=str(getattr(agent, "agent", "")),
            kind=kind,
            direction=want,
            setting=setting,
            delivery=DEL.delivery_for_kind(kind, stated=stated),
            effect=str(getattr(agent, "note", "") or ""),
            boundary="",
            source=str(getattr(agent, "source", "") or ""),
            option_id=f"REG-{norm(getattr(agent, 'agent', ''))}",
            binding="species",
        )
        if option.label and option.reach.reportable:
            out.append(option)
    return out


def _dedupe(options: Sequence[Option]) -> list[Option]:
    """One row per material, keeping the strongest evidence for it."""
    best: dict[str, Option] = {}
    for option in sorted(options, key=lambda o: o.sort_key):
        key = norm(option.label)[:60]
        best.setdefault(key, option)
    return sorted(best.values(), key=lambda o: o.sort_key)


#: The panels a reader would want to move, and which way. A panel that is
#: only ever informational is not a problem and gets no row.
PANEL_WANT: Final[dict[str, str]] = {
    "butyrate": WANT_UP, "propionate": WANT_UP, "fermentation": WANT_UP,
    "b12": WANT_UP, "folate": WANT_UP, "riboflavin": WANT_UP, "biotin": WANT_UP,
    "k2": WANT_UP, "thiamine": WANT_UP, "niacin": WANT_UP, "pantothenate": WANT_UP,
    "b6": WANT_UP, "indole": WANT_UP, "ipa": WANT_UP, "bsh": WANT_UP,
    "glutathione": WANT_UP, "spermidine": WANT_UP, "urolithin": WANT_UP, "urda": WANT_UP,
    "h2s": WANT_DOWN, "bai": WANT_DOWN, "pcresol": WANT_DOWN, "cutc": WANT_DOWN,
    "bglucuronidase": WANT_DOWN, "histamine": WANT_DOWN, "ethanol": WANT_DOWN,
    "putrescine": WANT_DOWN, "crc_virulence": WANT_DOWN, "mucinlps": WANT_DOWN,
    "methane": WANT_DOWN, "urease": WANT_DOWN, "oxalate": WANT_DOWN, "ammonia": WANT_DOWN,
    "bcaa": WANT_DOWN, "tmao": WANT_DOWN,
}

#: Band words that mean a panel is outside the usual range and in which
#: direction it sits. Both spellings of the report's own seven-band scale
#: are accepted, since the bands are written with spaces on the page and
#: with underscores in the registries.
_LOW = {"notably_low", "low", "somewhat_low", "below_average"}
_HIGH = {"notably_high", "high", "somewhat_high", "above_average"}

#: The report's seven-band scale, as percentile bounds and band names. The
#: panel rows record a percentile, not a band, so the band is derived here
#: with the same bounds the rest of the report draws.
_BOUNDS: Final = ((5.0, "notably_low"), (15.0, "low"), (25.0, "below_average"),
                  (75.0, "typical"), (85.0, "above_average"), (95.0, "high"))


def band_of(percentile: float | None) -> str:
    """Which of the seven bands a percentile falls in."""
    if percentile is None:
        return "unknown"
    for bound, name in _BOUNDS:
        if percentile <= bound:
            return name
    return "notably_high"


def panel_problems(panels: Sequence[Mapping[str, Any]], *, section: int) -> list[Problem]:
    """A problem for every metabolite pathway sitting outside the usual range.

    A pathway is a problem when it sits away from the middle *in the
    direction that matters for it*: low butyrate is a finding, high butyrate
    is not; high hydrogen sulfide is a finding, low is not.
    """
    out: list[Problem] = []
    for panel in panels:
        key = str(panel.get("id") or panel.get("panel") or "")
        band = str(panel.get("band") or panel.get("level") or "").replace(" ", "_").lower()
        if not band or band == "unknown":
            band = band_of(panel.get("percentile"))
        want = PANEL_WANT.get(key)
        if not key or not want:
            continue
        sits_low, sits_high = band in _LOW, band in _HIGH
        if not (sits_low or sits_high):
            continue
        if (want == WANT_UP and not sits_low) or (want == WANT_DOWN and not sits_high):
            continue
        label = str(panel.get("label") or panel.get("metabolite") or key.replace("_", " "))
        words = band.replace("_", " ")
        percentile = panel.get("percentile")
        problem = Problem(
            key=problem_key("panel", key),
            kind="panel",
            title=f"{label} capacity is {words}",
            detail=(f"{percentile:.0f}th percentile against the reference cohort"
                    if isinstance(percentile, (int, float)) else str(panel.get("detail") or "")),
            want=want,
            severity=3.0 if band.startswith("notably") else 2.0,
            section=section,
            detail_dest=_dest("panel", key),
        )
        problem.options = _dedupe(atlas_options(panels=[key], want=want))
        problem.worsens = _dedupe(atlas_options(panels=[key], want=_opposite(want)))
        problem.nothing_known = not problem.options
        out.append(problem)
    return out


def pathogen_problems(records: Sequence[Mapping[str, Any]], *, section: int) -> list[Problem]:
    """A problem for every organism the pathogen screen actually called."""
    out: list[Problem] = []
    for record in records:
        name = str(record.get("species") or record.get("name") or record.get("display") or "")
        if not name:
            continue
        problem = Problem(
            key=problem_key("pathogen", name),
            kind="pathogen",
            title=f"{name} was detected by the pathogen screen",
            detail=str(record.get("amount_words") or record.get("tier") or ""),
            want=WANT_DOWN,
            severity=10.0,
            section=section,
            detail_dest=_dest("pathogen", str(record.get("target_id") or name)),
        )
        problem.options = _dedupe(atlas_options(
            species=[name], genera=[genus_of(name)], concerns=["pathogen_burden"], want=WANT_DOWN,
        ))
        problem.worsens = _dedupe(atlas_options(species=[name], want=WANT_UP))
        problem.nothing_known = not problem.options
        out.append(problem)
    return out


def concern_problems(flags: Mapping[str, Any], *, section: int) -> list[Problem]:
    """Problems that are about the gut as a whole rather than one organism."""
    out: list[Problem] = []
    for key, entry in (flags or {}).items():
        if not entry:
            continue
        title = str((entry or {}).get("title") or key.replace("_", " ").capitalize())
        problem = Problem(
            key=problem_key("concern", key),
            kind="concern",
            title=title,
            detail=str((entry or {}).get("detail") or ""),
            want=str((entry or {}).get("want") or WANT_DOWN),
            severity=float((entry or {}).get("severity") or 1.0),
            section=section,
        )
        problem.options = _dedupe(atlas_options(concerns=[key], want=problem.want))
        problem.nothing_known = not problem.options
        out.append(problem)
    return out


#: The evidence lanes of the curated registry, as settings here. The lanes
#: already record what kind of study stands behind an exact protocol.
_LANE_SETTING: Final[dict[str, str]] = {
    "A": "human_trial", "B": "human_trial", "C": "human_trial",
    "D": "human_observational", "E": "in_vitro", "X": "human_trial",
}

#: And the registry's own intervention classes, as material kinds.
_CLASS_KIND: Final[dict[str, str]] = {
    "probiotic": "probiotic", "live_biotherapeutic": "probiotic",
    "multi_ingredient_product": "compound", "botanical_supplement": "herb",
    "vitamin_supplement": "compound", "prescription_drug": "drug",
    "otc_drug": "drug", "clinical_route": "clinical_route",
    "fmt_microbiota_product": "procedure", "diet": "diet", "lifestyle": "lifestyle",
    "prebiotic": "prebiotic", "food": "food",
}


#: Registry cards that are not things to do. `no_supported_targeted_intervention`
#: says there is nothing targeted; `inform_only` explains a reading without
#: proposing anything. Both carry a context note worth reading, and neither is
#: an option: counting them as one put the sentence "no supported targeted
#: intervention" in the column of things that might help, and inflated the
#: count of studied options on page one.
NOT_AN_OPTION: Final[frozenset[str]] = frozenset({
    "no_supported_targeted_intervention", "inform_only", "not_actionable",
})


def registry_notes(results: Mapping[str, Any] | None) -> dict[str, list[str]]:
    """What the registry says about findings it has no action for, by trigger.

    The reasoning on those cards is the best explanation the report has of
    *why* there is nothing, so it is kept and shown beside "nothing on
    record" rather than discarded with the card.
    """
    out: dict[str, list[str]] = {}
    findings = (results or {}).get("findings_and_evidence") or {}
    for summary in (findings.get("evidence_summaries") or ()):
        if str(summary.get("render")) not in NOT_AN_OPTION:
            continue
        note = str(summary.get("context_note") or "").strip()
        if not note:
            continue
        for trigger_id in (summary.get("trigger_ids") or ()):
            out.setdefault(str(trigger_id), []).append(note)
    return out


def registry_options(results: Mapping[str, Any] | None) -> dict[str, list[Option]]:
    """The curated exact-protocol cards, keyed by the trigger that fired them.

    These are the registry's own evidence summaries: an exact product, the
    protocol that was actually studied, the population it was studied in and
    the safety rules that apply. They are the strongest material the report
    has, and they reach the actions section through the same shape as
    everything else rather than through a block of their own.

    Only cards that actually propose something come through here; see
    :data:`NOT_AN_OPTION`.
    """
    out: dict[str, list[Option]] = {}
    findings = (results or {}).get("findings_and_evidence") or {}
    for summary in (findings.get("evidence_summaries") or ()):
        if str(summary.get("render")) in NOT_AN_OPTION:
            continue
        label = str(summary.get("intervention_display") or "").strip()
        if not label:
            # A card with no intervention is not an option, whatever its
            # title reads like.
            continue
        kind = _CLASS_KIND.get(str(summary.get("intervention_class") or ""), "compound")
        protocol = str(summary.get("studied_protocol_text") or "")
        # In participant mode the protocol is redacted to a placeholder; it
        # carries no information and only crowds out the population, which
        # does.
        if "not shown in participant mode" in protocol:
            protocol = ""
        population = str(summary.get("population_text") or "")
        option = Option(
            label=label if len(label) <= 150 else label[:149] + "\u2026",
            kind=kind,
            direction="registry",
            setting=_LANE_SETTING.get(str(summary.get("supporting_evidence_lane") or ""), "human_trial"),
            delivery="oral_tested",
            effect=" ".join(x for x in (protocol, f"Studied in {population}." if population else "") if x)[:480],
            # The context note is where a card's limits live - who it was
            # studied in, what it is not, what it has been shown to harm.
            # It is not trimmed to fit a column.
            boundary=str(summary.get("context_note") or "")[:600],
            source=str((summary.get("sources") or [{}])[0].get("url") or "")
            if isinstance(summary.get("sources"), list) and summary.get("sources") else "",
            option_id=str(summary.get("summary_id") or summary.get("action_mapping_id") or ""),
            binding="strain",
        )
        for trigger_id in (summary.get("trigger_ids") or ()):
            out.setdefault(str(trigger_id), []).append(option)
    return out


def pattern_problems(
    results: Mapping[str, Any] | None, *, section: int
) -> list[Problem]:
    """A problem for every disease pattern this sample resembles.

    The resemblance is a research observation about a community, not a
    diagnosis, and the only evidence that bears on it is the curated
    registry: exact products studied in people who actually have the
    condition. A pattern with nothing behind it says so.
    """
    findings = (results or {}).get("findings_and_evidence") or {}
    by_trigger = registry_options(results)
    out: list[Problem] = []
    for trigger in (findings.get("triggers") or ()):
        if trigger.get("trigger_type") != "research_profile_resemblance":
            continue
        level = str(trigger.get("level") or "")
        if level not in {"high", "notably_high"}:
            continue
        display = str(trigger.get("display") or trigger.get("target") or "")
        problem = Problem(
            key=problem_key("pattern", str(trigger.get("target") or display)),
            kind="pattern",
            title=display,
            detail="A community-level resemblance to a published pattern, not a diagnosis",
            want="registry",
            severity=2.5 if level == "notably_high" else 1.5,
            section=section,
            detail_dest=_dest("pattern", str(trigger.get("target") or "")),
        )
        problem.options = _dedupe(by_trigger.get(str(trigger.get("trigger_id")), []))
        problem.nothing_known = not problem.options
        out.append(problem)
    return out


def attach_registry(problems: Sequence[Problem], results: Mapping[str, Any] | None) -> None:
    """Add the curated cards to the organism and pathway problems they match.

    A registry card fired by "butyrate gene capacity: low" belongs on the
    butyrate problem, at the top, because an exact studied protocol is the
    strongest thing the report can offer. Cards that propose nothing bring
    their reasoning instead, which is what a problem with no options has to
    say for itself.
    """
    findings = (results or {}).get("findings_and_evidence") or {}
    by_trigger = registry_options(results)
    notes_by_trigger = registry_notes(results)
    _KIND = {
        "functional_capacity_finding": "panel",
        "ecological_finding": "panel",
        "exact_microbiome_feature": "organism",
    }
    by_target: dict[tuple[str, str], list[Option]] = {}
    notes_by_target: dict[tuple[str, str], list[str]] = {}
    for trigger in (findings.get("triggers") or ()):
        kind = _KIND.get(str(trigger.get("trigger_type")))
        if not kind:
            continue
        key = (kind, norm(trigger.get("target")))
        options = by_trigger.get(str(trigger.get("trigger_id")))
        if options:
            by_target.setdefault(key, []).extend(options)
        for note in notes_by_trigger.get(str(trigger.get("trigger_id")), ()):
            notes_by_target.setdefault(key, []).append(note)
    for problem in problems:
        _, _, target = problem.key.partition(":")
        key = (problem.kind, norm(target))
        extra = by_target.get(key)
        if extra:
            problem.options = _dedupe([*extra, *problem.options])
        for note in notes_by_target.get(key, ()):
            if note not in problem.notes:
                problem.notes.append(note)
        problem.nothing_known = not problem.options


def ordered(problems: Iterable[Problem]) -> list[Problem]:
    """Worst first, and inside that the ones with something to do about them.

    A pathogen outranks an overgrowth, an overgrowth outranks a pathway
    reading. Within a tier, a problem with options comes before one without,
    so a reader meets the actionable findings first and the open ones
    together at the end.
    """
    tier = {"pathogen": 0, "organism": 1, "concern": 2, "panel": 3, "pattern": 4}
    return sorted(
        problems,
        key=lambda p: (tier.get(p.kind, 5), p.nothing_known, -p.severity, p.title.lower()),
    )


def footnotes(problems: Iterable[Problem]) -> list[str]:
    """The fine print the shown options actually earned."""
    used: set[str] = set()
    for problem in problems:
        used |= problem.markers
    return DEL.footnotes_for(used)


__all__ = [
    "Option",
    "attach_registry",
    "Problem",
    "atlas_options",
    "concern_problems",
    "footnotes",
    "norm",
    "ordered",
    "organism_keys",
    "organism_problems",
    "panel_problems",
    "pathogen_problems",
    "pattern_problems",
    "registry_options",
]
