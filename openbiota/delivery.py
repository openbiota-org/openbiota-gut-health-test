"""Which evidence counts as a real option, and what has to be said about it.

The rule here replaces the one the report used to apply by accident, which
was to lean on human trials and let everything else fall off the page. That
threw away most of what is actually known. A great deal of the useful work on
gut organisms - especially on botanicals, essential-oil constituents, enzymes
and antimicrobial food compounds - has only ever been done in culture or in
mice, and a reader who is told nothing is better off knowing that rosemary
extract inhibited this organism in a dish than being told there is nothing.

So the design of the study never excludes a finding. **Every setting counts:
a randomised human trial, a mouse study, a dish.** What matters instead is
whether the material can reach the place it would have to act:

* a study that gave the material to a living animal or person has already
  answered that question, whatever else it did or did not show;
* a result from a dish is a real result about the organism, and becomes a
  real option whenever the material can arrive in the gut intact - either
  because its chemistry survives the stomach unaided, or because it is sold
  in an enteric-protected form that gets it there;
* only a material that cannot be swallowed at all, or that is destroyed and
  absorbed long before the colon with no protected preparation available, is
  held back - and even then the finding stays on record with its reason.

What the reader sees is the finding, with the setting named in a word and a
small marker. The fine print sits once at the foot of the table. It never
goes in the introduction, and it never swallows the finding.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

# --------------------------------------------------------------------------- #
# the two axes
# --------------------------------------------------------------------------- #

#: Where the observation was made. Ordered strongest-first for sorting only;
#: no rank here removes anything from the report.
SETTINGS: Final[tuple[str, ...]] = (
    "human_trial",
    "human_observational",
    "animal",
    "ex_vivo",
    "in_vitro",
    # A documented food rather than a study: what it is made from, how it is
    # fermented and whether it still contains live organisms when you eat it.
    # That is a real thing to know and a real thing to choose, and it is not
    # a trial, so it says so.
    "food_identity",
    "mechanistic",
)

#: Reader-facing words, and the marker that carries the caveat to the
#: footnote. A human finding needs no marker.
SETTING_WORDS: Final[dict[str, tuple[str, str]]] = {
    "human_trial": ("human trial", ""),
    "human_observational": ("human cohorts", ""),
    "animal": ("animal study", "\u2020"),
    "ex_vivo": ("human tissue or stool, outside the body", "\u2021"),
    "in_vitro": ("laboratory culture", "\u2021"),
    "food_identity": ("documented fermented food", "\u00a4"),
    "mechanistic": ("mechanism only", "\u2021"),
}

SETTING_RANK: Final[dict[str, int]] = {s: i for i, s in enumerate(SETTINGS)}

#: How the material gets to where it would have to act.
#:
#: ``oral_tested``       the study itself gave it by mouth to a living body.
#: ``inherently_stable`` its chemistry survives the stomach unaided: fibres and
#:                       resistant starches, polyphenol glycosides, bacterial
#:                       spores, mineral salts, whole foods, most herbs eaten
#:                       as food.
#: ``enteric_required``  it arrives intact only in an enteric-coated or
#:                       delayed-release preparation. Those preparations are
#:                       ordinary and sold; this is a note, not a refusal.
#: ``absorbed_upstream`` it is absorbed in the small intestine, so a colonic
#:                       effect needs the protected form as well.
#: ``not_ingestible``    topical, injected, inhaled or otherwise not a thing
#:                       that is swallowed. This is the only exclusion.
#: ``unknown``           the source does not say. Carried, and said.
DELIVERIES: Final[tuple[str, ...]] = (
    "oral_tested",
    "inherently_stable",
    "enteric_required",
    "absorbed_upstream",
    "not_ingestible",
    "unknown",
)

#: The marker and the footnote each delivery earns. Only the routes that
#: change what a reader would have to buy carry a marker.
DELIVERY_WORDS: Final[dict[str, tuple[str, str]]] = {
    "oral_tested": ("taken by mouth in the study", ""),
    "inherently_stable": ("reaches the gut unaided", ""),
    "enteric_required": ("needs an enteric-protected form", "\u00a7"),
    "absorbed_upstream": ("needs an enteric-protected form", "\u00a7"),
    "not_ingestible": ("not something that is swallowed", ""),
    "unknown": ("route not stated", "\u00b6"),
}

#: The footnotes, in the order their markers first appear. One line each: the
#: reader should be able to take them in without leaving the table.
FOOTNOTES: Final[dict[str, str]] = {
    "\u2020": "Shown in animals, not yet in people.",
    "\u2021": "Shown in laboratory culture rather than in a living body: it tells you the "
              "organism is susceptible to the substance, not that eating it will do this in you.",
    "\u00a7": "Would need an enteric-coated or delayed-release preparation to reach the colon "
              "intact; the plain form is absorbed or destroyed before it arrives.",
    "\u00a4": "A documented fermented food: what it is made from, how it is fermented and "
              "whether it still contains live organisms when eaten. Not a trial showing it "
              "changes this reading.",
    "\u00b6": "The source does not say how it was delivered.",
}


@dataclass(frozen=True, slots=True)
class Reach:
    """Whether a finding is an option here, and what has to be said about it."""

    reportable: bool
    reason: str
    markers: tuple[str, ...]
    setting_words: str
    delivery_words: str

    @property
    def marker_text(self) -> str:
        """The superscript to put after the option's name, or an empty string."""
        return "".join(self.markers)


def assess(setting: str, delivery: str) -> Reach:
    """Does this finding belong in front of a reader, and with what caveat?

    A study in a living body has already shown the material can get where it
    acts. A study in a dish has not, so the material's own route decides -
    and the only thing that disqualifies it is not being swallowable at all.
    """
    setting = setting if setting in SETTING_WORDS else "mechanistic"
    delivery = delivery if delivery in DELIVERY_WORDS else "unknown"
    in_a_body = setting in {"human_trial", "human_observational", "animal"}

    s_words, s_mark = SETTING_WORDS[setting]
    d_words, d_mark = DELIVERY_WORDS[delivery]

    if delivery == "not_ingestible" and not in_a_body:
        return Reach(
            False,
            "the material is not swallowed, so a culture result about it is not an option here",
            (), s_words, d_words,
        )
    markers = tuple(m for m in (s_mark, d_mark if not in_a_body else "") if m)
    return Reach(True, "", markers, s_words, d_words)


def footnotes_for(markers: object) -> list[str]:
    """The fine print the markers used on a page actually need, in order.

    Only the notes that were earned: a page with no animal evidence on it
    does not carry a line about animal evidence.
    """
    used = set(markers if isinstance(markers, (set, list, tuple)) else ())
    return [f"{m} {FOOTNOTES[m]}" for m in FOOTNOTES if m in used]


def sort_key(setting: str, delivery: str) -> tuple[int, int]:
    """Strongest setting first, then the route that needs least of the reader.

    An ordering, not a filter: a laboratory finding sorts below a trial and
    is still on the page.
    """
    return (
        SETTING_RANK.get(setting, len(SETTINGS)),
        DELIVERIES.index(delivery) if delivery in DELIVERIES else len(DELIVERIES),
    )


#: The legacy words the organism registry uses, mapped onto the settings here.
LEGACY_SETTING: Final[dict[str, str]] = {
    "human_trial": "human_trial",
    "human_observational": "human_observational",
    "animal": "animal",
    "in_vitro": "in_vitro",
    "ex_vivo": "ex_vivo",
    "mechanistic": "mechanistic",
}

#: Material kinds, and whether that kind of thing reaches the colon unaided.
#: Used when a source names a material but never states a route.
KIND_DELIVERY: Final[dict[str, str]] = {
    "prebiotic": "inherently_stable",
    "fibre": "inherently_stable",
    "food": "inherently_stable",
    "diet": "oral_tested",
    "lifestyle": "oral_tested",
    "probiotic": "inherently_stable",
    "postbiotic": "inherently_stable",
    "spore_probiotic": "inherently_stable",
    "mineral": "inherently_stable",
    "polyphenol": "inherently_stable",
    "herb": "inherently_stable",
    "enzyme": "enteric_required",
    "essential_oil": "enteric_required",
    "compound": "unknown",
    "drug": "oral_tested",
    "procedure": "oral_tested",
    "clinical_route": "oral_tested",
}


def delivery_for_kind(kind: str, *, stated: str | None = None) -> str:
    """The delivery route to use for a material whose source did not state one."""
    if stated in DELIVERIES and stated != "unknown":
        return stated
    return KIND_DELIVERY.get(kind, "unknown")


__all__ = [
    "DELIVERIES",
    "DELIVERY_WORDS",
    "FOOTNOTES",
    "SETTINGS",
    "SETTING_WORDS",
    "Reach",
    "assess",
    "delivery_for_kind",
    "footnotes_for",
    "sort_key",
]
