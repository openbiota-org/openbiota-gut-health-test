"""Single-line labels for the page-1 dashboard.

Page 1 is a grid, and a grid only reads as one if every row is the same
height. The full descriptive labels are written for the detail sections,
where a reading gets a whole card and a two-line name costs nothing; in a
half-width summary column the same name wraps to two or three lines, and the
page becomes as ragged as whatever happened to score highly in that sample.

So each label has a short form here, and the rule is that the short form has
to fit on one line at summary width. Everywhere else in the report keeps the
long name, because that is where the precision belongs.

THE RULE IS MEASURED, NOT GUESSED. A character limit let six names through
that wrapped on the page ("Dietary fibre and carbohydrate breakdown" is 40
characters and 164 points; the column is 99). `fits_page_one` measures the
string in the font and size page 1 actually uses, against the width the
column actually has, and `tests/test_page_one_names_fit.py` runs that
measurement over every panel in `panels/` and every profile in
`profiles/` - not a sample's top ten, all of them, because any of them can
be the one that scores highly next time. Add a panel or a profile and the
test tells you whether it needs an entry here before a report is built.

Underneath the table, `pdfsummary._scale_rows` refuses to wrap a name at
all: a name that somehow reaches it over budget is drawn smaller rather
than on two lines. That is a backstop for a reader, not a licence; the
test is what keeps the table honest.

Two mechanisms, in order:

1. An explicit table, for names where only a human can decide what to drop.
   "Branched-chain amino acids (valine, leucine, isoleucine)" becomes
   "BCAAs" because that is what a reader recognises, not because any rule
   could derive it.
2. A fallback that trims the patterns the long labels share — the
   "— gut pattern" suffix every disease profile carries, and trailing
   parentheticals — so a profile added later is still short by default
   rather than silently ragged.
"""

from __future__ import annotations

import re
from typing import Final

__all__ = [
    "PAGE_ONE_FONT", "PAGE_ONE_MARGIN", "PAGE_ONE_SIZE",
    "fits_page_one", "short_metabolite", "short_profile_label", "text_width",
]

#: Disease patterns, by profile name. Keyed on the stable identifier rather
#: than the label so rewording a label cannot silently drop the short form.
_PROFILE_SHORT: Final[dict[str, str]] = {
    "acvd": "Atherosclerosis",
    "ad_clinical": "Alzheimer's disease",
    "ad_mci": "Mild cognitive impairment",
    "ad_preclinical_amyloid": "Preclinical amyloid",
    "adenoma": "Colorectal adenoma",
    "alopecia_areata": "Alopecia areata",
    "androgenetic_alopecia": "Androgenetic alopecia",
    "ankylosing_spondylitis": "Ankylosing spondylitis",
    "celiac": "Coeliac disease",
    "cirrhosis": "Liver cirrhosis",
    "ckd": "Chronic kidney disease",
    "crc": "Colorectal cancer",
    "crohns": "Crohn's disease",
    "csu": "Chronic hives",
    "hypertension": "Hypertension",
    "ibd": "IBD (composite)",
    "ibs": "Irritable bowel syndrome",
    "ibs_c": "IBS with constipation",
    "ibs_d": "IBS with diarrhoea",
    "ibs_m": "IBS, mixed type",
    "longcovid": "Long COVID",
    "masld": "Fatty liver (MASLD)",
    "mdd": "Major depression",
    "mecfs": "ME/CFS",
    "ms": "Multiple sclerosis",
    "obesity": "Obesity",
    "parkinsons": "Parkinson's disease",
    "ra": "Rheumatoid arthritis",
    "sle": "Lupus, non-specific",
    "symptomatic_dermographism": "Skin-stroking hives",
    "t1d": "Type 1 diabetes",
    "t2d": "Type 2 diabetes",
    "uc": "Ulcerative colitis",
}

#: Metabolic readings, by the display name the panel carries.
_METABOLITE_SHORT: Final[dict[str, str]] = {
    "Vitamin B12 (cobalamin)": "Vitamin B12",
    "Secondary bile acids": "Secondary bile acids",
    "Branched-chain amino acids (valine, leucine, isoleucine)": "BCAAs",
    "Deconjugated glucuronides (oestrogens, drug metabolites, bilirubin)":
        "Glucuronidase activity",
    "Biotin (vitamin B7)": "Biotin",
    "Unconjugated bile acids": "Unconjugated bile acids",
    "Butyrate (SCFA)": "Butyrate",
    "Colibactin and bacterial toxins (virulence factors, not metabolites)":
        "Colibactin and toxins",
    "Trimethylamine / TMAO": "TMAO",
    "Dopamine (from L-DOPA) and tyramine (from tyrosine)": "Dopamine and tyramine",
    "Folate (vitamin B9)": "Folate",
    "gamma-aminobutyric acid (GABA)": "GABA",
    "Hydrogen sulfide (H2S)": "Hydrogen sulfide",
    "Histamine": "Histamine",
    "Indole / indoxyl sulfate": "Indole",
    "Indole-3-propionate (IPA)": "Indole-3-propionate",
    "Vitamin K2 (menaquinone, MK-n)": "Vitamin K2",
    "Methane (CH4)": "Methane",
    "Oxalate (degraded, not produced)": "Oxalate breakdown",
    "p-Cresol / p-cresyl sulfate": "p-Cresol",
    "Propionate (SCFA)": "Propionate",
    "Riboflavin (vitamin B2)": "Riboflavin",
    "Tryptamine": "Tryptamine",
    "Imidazole propionate (ImP)": "Imidazole propionate",
    "Ammonia (from urea)": "Ammonia",
    # The six that wrapped. Measured, not guessed: each is under the page-1
    # budget with the direction arrow beside it (see `fits_page_one`).
    "Dietary fibre and carbohydrate breakdown": "Fibre breakdown",
    "Microbial ethanol formation": "Ethanol formation",
    "Fermentation routes and cross-feeding": "Fermentation routes",
    "Glucosinolate to isothiocyanate conversion": "Glucosinolate use",
    "Mucin glycans and lipid A": "Mucin glycans & lipid A",
    "Protein and nitrogen metabolism": "Nitrogen metabolism",
}

#: Suffixes every disease label shares. Stripped before anything else so the
#: fallback starts from the disease name rather than the boilerplate.
_SUFFIXES: Final[tuple[str, ...]] = (
    " — gut pattern",
    " — ecological pattern",
    " — stage-specific gut pattern",
    " — translated 16S pattern",
    " gut-microbiome similarity",
    " chronic-state community similarity",
    "-associated community pattern",
    " community similarity",
)

_TRAILING_PAREN: Final = re.compile(r"\s*\([^)]*\)\s*$")


def _trim(label: str) -> str:
    """Strip the shared boilerplate, then any trailing parenthetical."""
    text = label.strip()
    for suffix in _SUFFIXES:
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    # An em-dash clause after the disease name is always elaboration.
    text = text.split(" — ")[0]
    text = _TRAILING_PAREN.sub("", text)
    return text.strip(" —-") or label


def short_profile_label(name: str, label: str) -> str:
    """One-line name for a disease pattern on page 1."""
    explicit = _PROFILE_SHORT.get(name)
    return explicit if explicit else _trim(label)


def short_metabolite(name: str) -> str:
    """One-line name for a metabolic reading on page 1."""
    explicit = _METABOLITE_SHORT.get(name)
    return explicit if explicit else _trim(name)


#: What page 1 draws a row name in: the "cell" style, bold.
PAGE_ONE_FONT: Final = "Helvetica-Bold"
PAGE_ONE_SIZE: Final = 7.8
#: Headroom under the column width. Font metrics are deterministic, so any
#: positive margin is safe; this one exists so a name is never decided by
#: a fraction of a point.
PAGE_ONE_MARGIN: Final = 1.0


def text_width(text: str, *, font: str = PAGE_ONE_FONT, size: float = PAGE_ONE_SIZE) -> float:
    """The width the text takes on the page, from the font's own metrics."""
    from reportlab.pdfbase.pdfmetrics import stringWidth  # noqa: PLC0415

    return float(stringWidth(text, font, size))


def fits_page_one(
    text: str, budget: float, *, font: str = PAGE_ONE_FONT, size: float = PAGE_ONE_SIZE,
) -> bool:
    """Does this name sit on one line in a page-1 column `budget` points wide?"""
    return text_width(text, font=font, size=size) + PAGE_ONE_MARGIN <= budget


def self_test() -> int:
    """Every short form must fit page 1, measured, and the fallback must do real work."""
    from openbiota.pdfsummary import page_one_name_budget  # noqa: PLC0415 - avoids a cycle

    checks = 0
    with_arrow = page_one_name_budget(arrow=True)
    without = page_one_name_budget(arrow=False)
    for key, value in _METABOLITE_SHORT.items():
        assert fits_page_one(value, with_arrow), (
            f"{key}: {value!r} is {text_width(value):.1f}pt; page 1 allows {with_arrow:.1f}"
        )
        assert value == value.strip()
        checks += 1
    for key, value in _PROFILE_SHORT.items():
        assert fits_page_one(value, without), (
            f"{key}: {value!r} is {text_width(value):.1f}pt; page 1 allows {without:.1f}"
        )
        assert value == value.strip()
        checks += 1

    # The fallback strips what it claims to strip.
    assert _trim("Crohn's disease — gut pattern") == "Crohn's disease"
    assert _trim("Obesity — ecological pattern") == "Obesity"
    assert _trim("Type 2 diabetes — gut pattern (metformin-gated)") == "Type 2 diabetes"
    assert _trim("Colorectal cancer-associated community pattern") == "Colorectal cancer"
    assert _trim("Biotin (vitamin B7)") == "Biotin"
    checks += 5

    # A profile with no entry still gets something short.
    assert short_profile_label("unknown", "Some disease — gut pattern") == "Some disease"
    # And a label that trims to nothing falls back rather than vanishing.
    assert short_profile_label("x", "— gut pattern")
    checks += 2

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} short-name checks passed")
