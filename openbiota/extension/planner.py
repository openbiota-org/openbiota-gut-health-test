"""Consolidated action planner — A12, BUILD_SPEC_v0.8.3 §8.

The report already produces an evidence card for every finding: what has
been studied about changing it, in whom, with the evidence for and against
side by side. What it has not had is a view *across* those cards — the one
a reader actually wants, which is "of all of that, what would I do first?"

This builds that view without replacing anything underneath it. Every
option here points back at the cards it came from, and nothing is invented
that is not already in the registry.

Four rules shape it.

**Ranking is prioritisation, not prediction.** The sort is deterministic
and every reason is printed. A higher position means more of this person's
findings are addressed by better-matched evidence — it is not a probability
that the option will work.

**A missing trial is not a missing option.** Laboratory and animal evidence
appears in the normal list with its actual setting named. Hiding it behind
a human-only default would turn "nobody has run the trial" into "nothing
can be done".

**A contraindication removes an option from *this person's* start-here
list, and from nothing else.** The evidence record stays, visible, with the
reason it does not apply here. An unknown fact is not a contraindication.

**A study's exposure is not a dose.** The protocol is shown as what was
done in the study; no cell-culture concentration or animal mg/kg is
converted into something to take.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

SCHEMA_VERSION: Final = "openbiota.action-planner/1.0"

#: How many options the "start here" list may hold. §8.2 says one to three:
#: a list of ten is a catalogue, and a reader does ten things by doing none.
START_HERE_MAX: Final = 3

#: The categories §8.2 requires, in the order they are offered.
CATEGORIES: Final[tuple[tuple[str, str], ...]] = (
    ("food", "Foods"),
    ("prebiotic", "Prebiotics and fibres"),
    ("probiotic", "Probiotics"),
    ("supplement", "Supplements and compounds"),
    ("herb", "Herbs and botanicals"),
    ("lifestyle", "Lifestyle"),
    ("medication", "Medication routes (clinician review)"),
)
CATEGORY_LABELS: Final[Mapping[str, str]] = dict(CATEGORIES)

#: Evidence lanes, best-matched first. The letters are the registry's own.
#: A lane is a statement about what was studied, never about how well it
#: worked.
LANE_ORDER: Final[Mapping[str, int]] = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "X": 5}

LANE_WORDS: Final[Mapping[str, str]] = {
    "A": "human trial in this condition",
    "B": "human trial in a related group",
    "C": "human study, other design",
    "D": "animal or ex-vivo study",
    "E": "laboratory study",
    "X": "mechanistic or indirect",
}

#: Why an option sits where it sits. Printed, per §8.2.
RANK_REASONS: Final[tuple[str, ...]] = (
    "addresses_more_findings",
    "better_matched_evidence",
    "matches_a_flagged_finding",
    "studied_in_a_comparable_population",
    "exact_formulation_studied",
    "no_contraindication_recorded",
    "feasible_as_a_food_or_supplement",
)

#: Food-checklist groupings, §8.3.
FOOD_GROUPS: Final[tuple[tuple[str, str, tuple[str, ...]], ...]] = (
    ("plant_foods", "Vegetables and fruit",
     ("vegetable", "fruit", "leafy", "berry", "apple", "banana", "onion", "garlic",
      "artichoke", "chicory", "leek", "asparagus", "broccoli", "kiwi")),
    ("whole_grains", "Whole grains",
     ("grain", "oat", "barley", "rye", "wheat", "bran", "wholegrain", "whole-grain",
      "sorghum", "millet", "buckwheat")),
    ("legumes", "Legumes and pulses",
     ("legume", "bean", "lentil", "chickpea", "pea", "soy", "pulse")),
    ("nuts_seeds", "Nuts and seeds",
     ("nut", "seed", "almond", "walnut", "flax", "chia", "hemp", "psyllium")),
    ("fermented", "Fermented foods",
     ("fermented", "yoghurt", "yogurt", "kefir", "kimchi", "sauerkraut", "kombucha",
      "miso", "tempeh", "cheese")),
    ("fibre_preparations", "Named fibre preparations",
     ("inulin", "gos", "fos", "resistant starch", "arabinoxylan", "beta-glucan",
      "pectin", "polydextrose", "partially hydrolysed guar", "maltodextrin",
      "xylooligosaccharide", "lactulose", "hmo", "2'fl", "lnnt")),
    ("polyphenol", "Polyphenol-rich foods and extracts",
     ("polyphenol", "cocoa", "cacao", "green tea", "pomegranate", "grape", "olive",
      "turmeric", "curcumin", "resveratrol", "urolithin", "equol")),
    ("compounds", "Supplements and compounds",
     ("vitamin", "riboflavin", "folate", "thiamine", "niacin", "biotin", "cobalamin",
      "butyrate", "glutamine", "zinc", "magnesium", "postbiotic", "tributyrin")),
)

#: Alias merges for the checklist. Materially different products stay apart:
#: cocoa powder and cocoa husk are not one item however similar the word.
FOOD_ALIASES: Final[Mapping[str, str]] = {
    "cacao": "cocoa",
    "cocoa powder": "cocoa",
    "beans": "legumes",
    "pulses": "legumes",
    "yogurt": "yoghurt",
    "wholegrain": "whole grain",
    "whole-grain": "whole grain",
}

#: Products that share a word but not a thing.
NEVER_MERGE: Final[tuple[tuple[str, str], ...]] = (
    ("cocoa husk", "cocoa"),
    ("cocoa butter", "cocoa"),
    ("wheat bran", "wheat"),
    ("oat bran", "oat"),
)

_WORD = re.compile(r"[^a-z0-9]+")


class PlannerError(ValueError):
    """An option that would misstate its own evidence or applicability."""


def _slug(text: str) -> str:
    return _WORD.sub("_", str(text).lower()).strip("_")


def normalise_food(name: str) -> str:
    """One name per food, without merging things that only sound alike."""
    text = str(name).strip().lower()
    for specific, _general in NEVER_MERGE:
        if specific in text:
            return specific
    return FOOD_ALIASES.get(text, text)


def food_group(name: str) -> tuple[str, str]:
    """Which checklist group a food belongs in."""
    text = normalise_food(name)
    for group_id, label, words in FOOD_GROUPS:
        if any(word in text for word in words):
            return group_id, label
    return "other", "Other"


# --------------------------------------------------------------------------- #
# one option
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Tradeoff:
    """Something the option might also do, that a reader should weigh.

    Kept separate from the evidence for it. §8.2: studied positive effects,
    possible undesired effects and uncertainties are shown together, and a
    count of favourable edges is not a net benefit.
    """

    kind: str          # undesired_effect | uncertainty | interaction | feeds_flagged_organism
    detail: str
    source: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {"kind": self.kind, "detail": self.detail, "source": self.source}


@dataclass(frozen=True)
class Exclusion:
    """Why an option is not on *this person's* start-here list.

    It stays in the catalogue. §8.2 is explicit that a contraindication
    affects the new planner only and deletes no evidence record.
    """

    reason: str
    detail: str
    from_supplied_fact: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "reason": self.reason, "detail": self.detail,
            "from_supplied_fact": self.from_supplied_fact,
            "note": (
                "This option stays in the catalogue with its evidence intact; it is only "
                "excluded from the start-here list for this person."
            ),
        }


@dataclass
class Option:
    """One consolidated action, across every finding it addresses."""

    option_id: str
    label: str
    category: str
    #: Findings this option has evidence about, by the report's own names.
    addresses: tuple[str, ...] = ()
    #: Findings this report flagged that this option addresses.
    addresses_flagged: tuple[str, ...] = ()
    best_lane: str = "X"
    lanes: tuple[str, ...] = ()
    n_supporting: int = 0
    n_against: int = 0
    formulation: str | None = None
    studied_exposure: str | None = None
    population: str | None = None
    monitoring_question: str | None = None
    tradeoffs: tuple[Tradeoff, ...] = ()
    exclusions: tuple[Exclusion, ...] = ()
    card_ids: tuple[str, ...] = ()
    #: Things a person could buy, from the registry's components.
    foods: tuple[str, ...] = ()
    #: A named dietary pattern, which is a way of eating rather than an item.
    pattern: str | None = None
    rank_reasons: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.category not in CATEGORY_LABELS:
            raise PlannerError(f"{self.option_id}: unknown category {self.category!r}")
        if self.best_lane not in LANE_ORDER:
            raise PlannerError(f"{self.option_id}: unknown evidence lane {self.best_lane!r}")

    @property
    def eligible_for_start_here(self) -> bool:
        return not self.exclusions

    @property
    def has_human_evidence(self) -> bool:
        return self.best_lane in {"A", "B", "C"}

    @property
    def sort_key(self) -> tuple[Any, ...]:
        """Deterministic, and every component is a printed reason.

        Order: flagged findings addressed, then findings addressed, then
        how well-matched the best evidence is, then supporting count, then
        the identifier so two equal options never swap between runs.
        """
        return (
            -len(self.addresses_flagged),
            -len(self.addresses),
            LANE_ORDER[self.best_lane],
            -self.n_supporting,
            self.option_id,
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "option_id": self.option_id,
            "label": self.label,
            "category": self.category,
            "category_label": CATEGORY_LABELS[self.category],
            "addresses": list(self.addresses),
            "addresses_flagged": list(self.addresses_flagged),
            "n_findings_addressed": len(self.addresses),
            "best_evidence_lane": self.best_lane,
            "best_evidence_words": LANE_WORDS[self.best_lane],
            "evidence_lanes": list(self.lanes),
            "n_supporting": self.n_supporting,
            "n_against": self.n_against,
            "counts_note": (
                "Counts of favourable and unfavourable findings. They are not a net benefit "
                "and they do not add up to a verdict."
            ),
            "formulation": self.formulation,
            "studied_exposure": self.studied_exposure,
            "exposure_note": (
                "What was done in the study, not a dose to take. A laboratory concentration "
                "or an animal exposure is not converted into a human regimen."
            ),
            "population": self.population,
            "monitoring_question": self.monitoring_question,
            "tradeoffs": [t.to_json() for t in self.tradeoffs],
            "exclusions": [e.to_json() for e in self.exclusions],
            "eligible_for_start_here": self.eligible_for_start_here,
            "evidence_cards": list(self.card_ids),
            "foods": list(self.foods),
            "dietary_pattern": self.pattern,
            "rank_reasons": list(self.rank_reasons),
        }


def rank_reasons_for(option: Option, *, most_findings: int) -> tuple[str, ...]:
    """The reasons this option sits where it does, in the sort's own order."""
    out: list[str] = []
    if option.addresses_flagged:
        out.append("matches_a_flagged_finding")
    if len(option.addresses) >= max(2, most_findings):
        out.append("addresses_more_findings")
    if option.best_lane in {"A", "B"}:
        out.append("better_matched_evidence")
    if option.population:
        out.append("studied_in_a_comparable_population")
    if option.formulation:
        out.append("exact_formulation_studied")
    if not option.exclusions:
        out.append("no_contraindication_recorded")
    if option.category in {"food", "prebiotic", "probiotic"}:
        out.append("feasible_as_a_food_or_supplement")
    return tuple(r for r in out if r in RANK_REASONS)


# --------------------------------------------------------------------------- #
# monitoring
# --------------------------------------------------------------------------- #

#: §8.2(8): every action carries a question that could actually be answered.
#: None of them promises colonisation or a score change.
MONITORING_BY_CATEGORY: Final[Mapping[str, str]] = {
    "food": (
        "Did your symptoms and tolerance change over four weeks of eating it regularly? "
        "A repeat sample after eight weeks would show whether the community followed."
    ),
    "prebiotic": (
        "Did tolerance hold as you built up the amount, and did the fibre-using and "
        "butyrate readings move on a repeat sample after 4 to 8 weeks?"
    ),
    "probiotic": (
        "Is the strain detectable while you are taking it, and does anything it was taken "
        "for change? A sample four weeks after stopping distinguishes passing through from "
        "settling in."
    ),
    "supplement": (
        "Did the endpoint the study measured change, and did you tolerate it? Some of these "
        "have a directly measurable marker; where they do, that is the thing to measure."
    ),
    "herb": (
        "Did you tolerate it, and did the specific endpoint the study measured change? "
        "Botanicals interact with medicines more often than their availability suggests, "
        "so check that first."
    ),
    "lifestyle": (
        "Did the change hold for eight weeks at all? That comes before asking whether "
        "anything measurable followed it."
    ),
    "medication": (
        "This is a clinician's decision and their monitoring plan, not a self-directed "
        "action, and the endpoint is a clinical one rather than a microbiome reading."
    ),
}


def monitoring_question(category: str) -> str:
    return MONITORING_BY_CATEGORY.get(
        category,
        "What would you expect to change, and how would you know? An action with no "
        "answerable question is an action with no way to tell whether it worked.",
    )


# --------------------------------------------------------------------------- #
# the plan
# --------------------------------------------------------------------------- #


#: Component entity types that name something a person could actually buy.
#: A dietary *pattern* is not a shopping item and is listed separately.
SHOPPABLE: Final[frozenset[str]] = frozenset({
    "prebiotic", "probiotic_strain", "postbiotic_strain", "compound", "botanical",
    "vitamin", "multi_ingredient",
})


def shopping_items(identity: Mapping[str, Any] | None) -> tuple[tuple[str, str], ...]:
    """(item, entity type) for everything in an intervention you could buy.

    Read from the registry's own components. The display name is a study
    title - "High-fibre versus high-fermented-food diet (Wastyk 2021)" is
    not something to put on a shopping list, and printing it as one would
    make the checklist unusable.
    """
    out: list[tuple[str, str]] = []
    for component in (identity or {}).get("components") or []:
        kind = str(component.get("entity_type") or "")
        if kind not in SHOPPABLE:
            continue
        name = (
            component.get("molecule") or component.get("name")
            or component.get("strain") or component.get("species")
            or component.get("ingredient")
        )
        if name:
            out.append((str(name), kind))
    return tuple(dict.fromkeys(out))


def dietary_pattern(identity: Mapping[str, Any] | None) -> str | None:
    """The named dietary pattern, if this intervention is one."""
    for component in (identity or {}).get("components") or []:
        if str(component.get("entity_type")) == "diet" and component.get("name"):
            return str(component["name"])
    return None


@dataclass
class ChecklistItem:
    item: str
    group_id: str
    group_label: str
    for_options: tuple[str, ...]
    addresses: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "item": self.item, "group": self.group_id, "group_label": self.group_label,
            "from_options": list(self.for_options), "addresses": list(self.addresses),
        }


@dataclass
class Plan:
    """The consolidated view: what to do first, and everything else."""

    sample_id: str
    options: list[Option] = field(default_factory=list)
    goals: tuple[str, ...] = ()

    @property
    def ordered(self) -> list[Option]:
        return sorted(self.options, key=lambda o: o.sort_key)

    @property
    def start_here(self) -> list[Option]:
        """One to three options, each eligible and each doing different work.

        Two probiotics that address the same finding are one idea, not two,
        so the list prefers breadth across categories before depth.
        """
        picked: list[Option] = []
        seen_categories: set[str] = set()
        covered: set[str] = set()
        for option in self.ordered:
            if not option.eligible_for_start_here:
                continue
            fresh = set(option.addresses) - covered
            if picked and option.category in seen_categories and not fresh:
                continue
            picked.append(option)
            seen_categories.add(option.category)
            covered |= set(option.addresses)
            if len(picked) >= START_HERE_MAX:
                break
        return picked

    @property
    def excluded(self) -> list[Option]:
        return [o for o in self.ordered if o.exclusions]

    def by_category(self) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        for category, label in CATEGORIES:
            members = [o for o in self.ordered if o.category == category]
            if members:
                out.append({
                    "category": category, "category_label": label,
                    "n_options": len(members),
                    "options": [o.to_json() for o in members],
                })
        return out

    def checklist(self) -> list[dict[str, Any]]:
        """A deduplicated food list, grouped, each item naming what it is for."""
        merged: dict[str, ChecklistItem] = {}
        for option in self.ordered:
            if option.exclusions:
                continue
            for food in option.foods:
                name = normalise_food(food)
                group_id, group_label = food_group(name)
                existing = merged.get(name)
                if existing is None:
                    merged[name] = ChecklistItem(
                        item=name, group_id=group_id, group_label=group_label,
                        for_options=(option.option_id,), addresses=tuple(option.addresses),
                    )
                else:
                    merged[name] = ChecklistItem(
                        item=name, group_id=existing.group_id, group_label=existing.group_label,
                        for_options=(*existing.for_options, option.option_id),
                        addresses=tuple(dict.fromkeys((*existing.addresses, *option.addresses))),
                    )
        groups: dict[str, list[ChecklistItem]] = {}
        for item in merged.values():
            groups.setdefault(item.group_id, []).append(item)
        out: list[dict[str, Any]] = []
        for group_id, label, _words in (*FOOD_GROUPS, ("other", "Other", ())):
            members = sorted(groups.get(group_id, []), key=lambda i: i.item)
            if members:
                out.append({
                    "group": group_id, "group_label": label,
                    "items": [m.to_json() for m in members],
                })
        return out

    def patterns(self) -> list[dict[str, Any]]:
        """Named ways of eating, kept apart from items you can buy."""
        out: list[dict[str, Any]] = []
        for option in self.ordered:
            if option.pattern and not option.exclusions:
                out.append({
                    "pattern": option.pattern,
                    "from_option": option.option_id,
                    "label": option.label,
                    "addresses": list(option.addresses),
                    "studied_protocol": option.studied_exposure,
                })
        return out

    def coverage(self) -> dict[str, Any]:
        addressed = {f for o in self.options for f in o.addresses}
        flagged = {f for o in self.options for f in o.addresses_flagged}
        return {
            "n_options": len(self.options),
            "n_eligible": sum(1 for o in self.options if o.eligible_for_start_here),
            "n_excluded_for_this_person": len(self.excluded),
            "n_findings_addressed": len(addressed),
            "n_flagged_findings_addressed": len(flagged),
            "n_with_human_evidence": sum(1 for o in self.options if o.has_human_evidence),
            "n_without_human_evidence": sum(
                1 for o in self.options if not o.has_human_evidence
            ),
            "coverage_note": (
                "An option without a human trial is still an option; its evidence setting is "
                "named on the card. Absence of a trial is not absence of anything to try."
            ),
        }

    def to_json(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "feature_id": "A12",
            "sample_id": self.sample_id,
            "goals": list(self.goals),
            "start_here": [o.to_json() for o in self.start_here],
            "start_here_note": (
                f"At most {START_HERE_MAX} options, chosen to address different findings "
                "rather than to repeat one idea. This is a prioritisation with its reasons "
                "printed, not a prediction that any of them will work for you."
            ),
            "catalogue": self.by_category(),
            "excluded_for_this_person": [o.to_json() for o in self.excluded],
            "food_checklist": self.checklist(),
            "dietary_patterns": self.patterns(),
            "checklist_note": (
                "Items are read from each intervention's recorded components, so the list "
                "holds things you could buy rather than the titles of the studies they came "
                "from. A named way of eating is listed separately, because it is a pattern "
                "and not an item."
            ),
            "coverage": self.coverage(),
            "ranking_policy": {
                "order": [
                    "how many of your flagged findings it addresses",
                    "how many of your findings it addresses at all",
                    "how well the best available evidence matches the question",
                    "how much supporting evidence there is",
                    "its identifier, so two equal options never swap between runs",
                ],
                "reasons_vocabulary": list(RANK_REASONS),
                "note": (
                    "Prioritisation, not a predicted probability of benefit. Every reason is "
                    "printed on the option it applies to."
                ),
            },
            "limitations": [
                "An option is specific to the compound, preparation, strain, substrate, "
                "model and endpoint actually studied.",
                "A study exposure is not a dose. No laboratory concentration or animal "
                "exposure is converted into something to take.",
                "A combined result is not automatically an effect of either component alone.",
                "Conflicting and null evidence stays beside the positive evidence.",
            ],
        }


# --------------------------------------------------------------------------- #
# building from the existing evidence cards
# --------------------------------------------------------------------------- #

#: Registry intervention classes mapped onto the planner's categories.
CLASS_TO_CATEGORY: Final[Mapping[str, str]] = {
    "food": "food", "diet": "food", "dietary_pattern": "food", "whole_food": "food",
    "fiber": "prebiotic", "fibre": "prebiotic", "prebiotic": "prebiotic",
    "probiotic": "probiotic", "live_biotherapeutic": "probiotic", "synbiotic": "probiotic",
    "postbiotic": "supplement", "supplement": "supplement", "compound": "supplement",
    "micronutrient": "supplement", "vitamin": "supplement",
    "vitamin_supplement": "supplement", "mineral_supplement": "supplement",
    "botanical": "herb", "herb": "herb", "essential_oil": "herb",
    "botanical_supplement": "herb",
    "behaviour": "lifestyle", "behavior": "lifestyle", "exercise": "lifestyle",
    "sleep": "lifestyle", "lifestyle": "lifestyle",
    "drug": "medication", "medication": "medication", "antibiotic": "medication",
    "procedure": "medication", "fmt": "medication",
    "clinical_route": "medication", "prescription_drug": "medication",
    "fmt_microbiota_product": "medication",
}


def category_for(intervention_class: str | None, label: str = "") -> str:
    """The planner category for a registry class, defaulting conservatively."""
    if intervention_class:
        mapped = CLASS_TO_CATEGORY.get(_slug(intervention_class))
        if mapped:
            return mapped
    text = f"{intervention_class or ''} {label}".lower()
    for word, category in (
        ("lactobacill", "probiotic"), ("bifidobacter", "probiotic"),
        ("inulin", "prebiotic"), ("fibre", "prebiotic"), ("fiber", "prebiotic"),
        ("extract", "herb"), ("oil", "herb"),
    ):
        if word in text:
            return category
    return "supplement"


def addresses_flagged(displays: Sequence[str], flagged: Sequence[str]) -> tuple[str, ...]:
    """Which of these findings concern an organism the report flagged.

    A trigger reads "Faecalibacterium prausnitzii: low"; the flag is on
    "Faecalibacterium prausnitzii". Comparing the two whole strings finds
    nothing, which silently emptied the one ranking signal that matters
    most.
    """
    keys = {_slug(f) for f in flagged if f}
    out: list[str] = []
    for display in displays:
        name = _slug(str(display).split(":")[0])
        if name and name in keys:
            out.append(display)
    return tuple(out)


def build_plan(
    summaries: Iterable[Mapping[str, Any]],
    *,
    sample_id: str,
    flagged_findings: Sequence[str] = (),
    excluded: Mapping[str, Exclusion] | None = None,
    goals: Sequence[str] = (),
    identities: Mapping[str, Mapping[str, Any]] | None = None,
) -> Plan:
    """Consolidate the per-finding evidence cards into one ranked plan.

    Cards are grouped by the intervention they concern, so an option that
    appears under four findings becomes one option addressing four things
    rather than four options.
    """
    excluded = dict(excluded or {})
    identities = dict(identities or {})
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for card in summaries:
        key = str(card.get("intervention_id") or card.get("intervention_display") or "")
        if not key:
            continue
        grouped.setdefault(key, []).append(card)

    plan = Plan(sample_id=sample_id, goals=tuple(goals))
    for intervention_id, cards in grouped.items():
        first = cards[0]
        label = str(first.get("intervention_display") or intervention_id)
        category = category_for(first.get("intervention_class"), label)
        addresses = tuple(dict.fromkeys(
            display for card in cards for display in (card.get("trigger_displays") or [])
        ))
        lanes = tuple(dict.fromkeys(
            str(card.get("supporting_evidence_lane") or "X") for card in cards
        ))
        known = [lane for lane in lanes if lane in LANE_ORDER]
        best = min(known, key=lambda lane: LANE_ORDER[lane]) if known else "X"
        option = Option(
            option_id=_slug(intervention_id),
            label=label,
            category=category,
            addresses=addresses,
            addresses_flagged=addresses_flagged(addresses, flagged_findings),
            best_lane=best,
            lanes=lanes,
            n_supporting=sum(len(card.get("supporting") or []) for card in cards),
            n_against=sum(len(card.get("against") or []) for card in cards),
            formulation=first.get("studied_protocol_text"),
            studied_exposure=first.get("studied_protocol_text"),
            population=first.get("population_text"),
            monitoring_question=monitoring_question(category),
            tradeoffs=tuple(
                Tradeoff("undesired_effect", str(reason), card.get("summary_id"))
                for card in cards for reason in (card.get("against_evidence") or [])
            ),
            exclusions=(excluded[intervention_id],) if intervention_id in excluded else (),
            card_ids=tuple(str(card.get("summary_id")) for card in cards),
            foods=tuple(
                normalise_food(item) for item, kind in
                shopping_items((identities.get(intervention_id) or {}).get("identity"))
                if kind != "probiotic_strain"
            ),
            pattern=dietary_pattern((identities.get(intervention_id) or {}).get("identity")),
        )
        plan.options.append(option)

    most = max((len(o.addresses) for o in plan.options), default=0)
    for option in plan.options:
        object.__setattr__(option, "rank_reasons", rank_reasons_for(option, most_findings=most))
    return plan


__all__ = [
    "CATEGORIES",
    "CATEGORY_LABELS",
    "CLASS_TO_CATEGORY",
    "FOOD_ALIASES",
    "FOOD_GROUPS",
    "LANE_ORDER",
    "LANE_WORDS",
    "MONITORING_BY_CATEGORY",
    "NEVER_MERGE",
    "RANK_REASONS",
    "SCHEMA_VERSION",
    "START_HERE_MAX",
    "ChecklistItem",
    "Exclusion",
    "Option",
    "Plan",
    "PlannerError",
    "Tradeoff",
    "SHOPPABLE",
    "addresses_flagged",
    "build_plan",
    "category_for",
    "dietary_pattern",
    "food_group",
    "monitoring_question",
    "normalise_food",
    "rank_reasons_for",
    "shopping_items",
]
