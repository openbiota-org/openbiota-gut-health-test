"""AT005 — where each section sits, and that nothing moved by accident.

Section numbers are printed, cross-referenced in prose and used as link
destinations, so a section that shifts silently breaks references that
still look correct. The contents page is deliberately unnumbered, so the
summary is 1 and what stood out is 2; the input register sits before the
closing four sections.
"""

from __future__ import annotations

import pytest

from openbiota.pdfreport import SECTIONS

#: The order the reader meets the report in. Overview sections first, then
#: the detail pages that expand them, then the material about the test
#: itself. Written out rather than derived, so a reordering has to be made
#: here deliberately.
EXPECTED_ORDER: tuple[tuple[str, int], ...] = (
    # The contents page is unnumbered (0): it is the map, not a stop.
    ("summary", 1), ("guide", 0), ("stood_out", 2),
    ("community", 3), ("groups", 4), ("organisms", 5), ("catalogue", 6),
    ("pathogens", 7), ("age", 8), ("functions", 9), ("patterns", 10),
    ("biofilm", 11), ("mycobiome", 12), ("skin", 13), ("actions", 14),
    ("groups_detail", 15), ("organisms_detail", 16), ("strains", 17),
    ("pathogens_detail", 18), ("functions_detail", 19), ("patterns_detail", 20),
    ("biofilm_detail", 21), ("mycobiome_detail", 22), ("skin_detail", 23),
    ("context", 24),
    ("limits", 25), ("sequencing", 26), ("accuracy", 27), ("technical", 28),
)


def test_every_section_is_where_the_contract_puts_it() -> None:
    for key, number in EXPECTED_ORDER:
        assert SECTIONS.get(key) == number, (
            f"{key} is section {SECTIONS.get(key)} and the contract says {number}"
        )


def test_the_map_holds_nothing_the_contract_does_not_name() -> None:
    extra = sorted(set(SECTIONS) - {k for k, _ in EXPECTED_ORDER})
    assert not extra, f"sections nobody declared: {extra}"


def test_the_numbers_are_contiguous_and_unique() -> None:
    """A gap or a repeat makes a printed cross-reference ambiguous."""
    numbers = sorted(n for n in SECTIONS.values() if n)  # 0 is the unnumbered contents
    assert numbers == list(range(1, len(numbers) + 1))
    assert list(SECTIONS.values()).count(0) == 1


def test_the_input_register_precedes_the_closing_four() -> None:
    """§10.1's authorised insertion, and the reason the tail moved."""
    for key in ("limits", "sequencing", "accuracy", "technical"):
        assert SECTIONS[key] > SECTIONS["context"], f"{key} should follow the input register"


def test_the_contents_are_on_page_two() -> None:
    """The contents come before the findings.

    Page two used to be What stood out, which assumes the reader already
    knows where things are; the page that tells them was behind it.
    """
    # Unnumbered, and printed between the summary and what stood out: the
    # summary is section 1, what stood out is section 2.
    assert SECTIONS["guide"] == 0
    assert SECTIONS["summary"] == 1
    assert SECTIONS["stood_out"] == 2


def test_the_detail_pages_follow_every_overview() -> None:
    """A detail page before its overview would be read before its context."""
    for overview, detail in (
        ("groups", "groups_detail"), ("organisms", "organisms_detail"),
        ("pathogens", "pathogens_detail"), ("functions", "functions_detail"),
        ("patterns", "patterns_detail"), ("biofilm", "biofilm_detail"),
        ("mycobiome", "mycobiome_detail"), ("skin", "skin_detail"),
    ):
        assert SECTIONS[overview] < SECTIONS[detail], f"{detail} precedes {overview}"


def test_the_full_catalogue_follows_the_flagged_organisms() -> None:
    """Having read which organisms are flagged, "what else is in there" is
    the next question, and it should not need a jump to the back."""
    assert SECTIONS["catalogue"] == SECTIONS["organisms"] + 1


@pytest.mark.parametrize("key", [k for k, _ in EXPECTED_ORDER])
def test_every_section_has_a_link_destination(key: str) -> None:
    from openbiota.pdflinks import section_dest

    assert section_dest(key) == f"sec-{key}"


# --- the authorised title change ------------------------------------------


def test_the_functions_sections_carry_their_authorised_titles() -> None:
    """The rename was to these words, and the overview and detail pages
    have to agree on them or the cross-reference reads as a different
    section.

    "& Metabolism" came off because the heading no longer fitted on one
    line of page one, and a first page whose title wraps is a first page
    that looks different in every report.
    """
    from pathlib import Path

    source = (Path(__file__).resolve().parent.parent / "openbiota" / "pdfreport.py").read_text(
        encoding="utf-8"
    )
    assert "Your Microbial Functions at a Glance" in source
    assert "Your Microbial Functions in Detail" in source
    assert "Metabolism at a Glance" not in source, "the long title is back"
