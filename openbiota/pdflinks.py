"""Internal navigation for the report: destinations, links and the outline.

One rule, applied everywhere: a section-level element links to its
at-a-glance section, an item-level element links straight to its detail
card. So on page 1 the whole age block goes to the age section and the
whole diversity block to the community section, while each metabolite or
disease row goes to that reading's own card; in the at-a-glance sections
each row goes to its card; each card links back to its row; every
"section N" in the text goes to that section; and the sections populate
the reader's bookmarks panel.

The mechanism is the PDF specification's own — GoTo link annotations to
named destinations, in the format since version 1.1 — which every reader
follows. Destinations are named rather than numbered because section start
pages differ between samples: a name resolves to wherever its target was
drawn. ReportLab refuses to save a document containing a link to a
destination that was never defined, so a dead link is a build failure,
not something a reader discovers.

Links carry no border and no underline. Link text is set in the accent
colour, which is the report's one cue that something is clickable; table
rows that link are clickable across their whole width.

Destination names:

    sec-<key>            a section, keyed as in `pdfreport.SECTIONS`
    sec-<key>-<n>        a numbered sub-section (13.1 -> sec-groups_detail-1)
    fn-<panel>           a metabolic function's at-a-glance row
    fn-<panel>-detail    its detail card
    pr-<profile>         a disease pattern's at-a-glance row
    pr-<profile>-detail  its detail page
    grp-<group>[-detail]      a microbial group's row / card
    org-<species>[-detail]    an organism's row / card
    pth-<target>[-detail]     a pathogen finding's row / card
    skn-<panel>[-detail]      a skin panel's row / page
"""

from __future__ import annotations

import re
from typing import Any, Final

from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import Flowable, Table
from reportlab.platypus import Paragraph as _Paragraph

__all__ = [
    "Anchor",
    "LinkedTable",
    "Paragraph",
    "RowLink",
    "back_link",
    "back_link_line",
    "function_dest",
    "group_dest",
    "link_sections",
    "linked_block",
    "linked_cell",
    "organism_dest",
    "pathogen_dest",
    "profile_dest",
    "section_dest",
    "section_heading",
    "skin_dest",
]

#: Headroom above a destination, so the target is not flush with the top of
#: the reader's window when it lands.
HEADROOM: Final = 5 * mm

_SLUG: Final = re.compile(r"[^a-z0-9]+")


def _slug(text: str) -> str:
    return _SLUG.sub("-", text.lower()).strip("-")


def section_dest(key: str, sub: int | None = None) -> str:
    return f"sec-{key}" if sub is None else f"sec-{key}-{sub}"


def function_dest(panel: str, *, detail: bool = False) -> str:
    return f"fn-{_slug(panel)}" + ("-detail" if detail else "")


def profile_dest(name: str, *, detail: bool = False) -> str:
    return f"pr-{_slug(name)}" + ("-detail" if detail else "")


def group_dest(name: str, *, detail: bool = False) -> str:
    return f"grp-{_slug(name)}" + ("-detail" if detail else "")


def organism_dest(species: str, *, detail: bool = False) -> str:
    return f"org-{_slug(species)}" + ("-detail" if detail else "")


def pathogen_dest(target_id: str, *, detail: bool = False) -> str:
    return f"pth-{_slug(target_id)}" + ("-detail" if detail else "")


def skin_dest(panel_id: str, *, detail: bool = False) -> str:
    return f"skn-{_slug(panel_id)}" + ("-detail" if detail else "")


def _sections_and_colour() -> tuple[dict[int, str], str]:
    # Imported here, not at module level: pdfreport imports this module for
    # its Paragraph, so a top-level import would be circular.
    from openbiota.pdfreport import ACCENT, SECTIONS

    colour = f"#{int(ACCENT.red * 255):02X}{int(ACCENT.green * 255):02X}{int(ACCENT.blue * 255):02X}"
    return {number: key for key, number in SECTIONS.items()}, colour


# --------------------------------------------------------------------------- #
# "section N" in running text
# --------------------------------------------------------------------------- #

#: "section 8", "Section 7", "sections 4 to 6, 13, 14". The number must not
#: continue as "13.1" (a sub-section reference) or as more digits.
_SECTION_RE: Final = re.compile(
    r"\b([Ss]ections?)((?:\s|&nbsp;)+)(\d+(?:\s+to\s+\d+)?(?:,\s*\d+)*)(?!\.\d)(?!\d)"
)
_NUMBER_RE: Final = re.compile(r"\d+")


def _inside_link(text: str, pos: int) -> bool:
    """True if `pos` falls between an opening <a and its </a>."""
    return text.rfind("<a ", 0, pos) > text.rfind("</a>", 0, pos)


def link_sections(text: str) -> str:
    """Turn every "section N" mention into a link to that section.

    A single number links as a phrase ("section 8"); a list links each
    number on its own ("sections 4 to 6, 13, 14"). Numbers that are not
    sections are left alone, as is anything already inside a link.
    """
    if "ection" not in text:
        return text
    by_number, colour = _sections_and_colour()

    def anchor(number: str, label: str) -> str:
        return f'<a href="#{section_dest(by_number[int(number)])}" color="{colour}">{label}</a>'

    def replace(match: re.Match[str]) -> str:
        if _inside_link(text, match.start()):
            return match.group(0)
        word, gap, numbers = match.groups()
        found = _NUMBER_RE.findall(numbers)
        if not all(int(n) in by_number for n in found):
            return match.group(0)
        if len(found) == 1:
            return anchor(found[0], f"{word}{gap}{found[0]}")
        return word + gap + _NUMBER_RE.sub(lambda m: anchor(m.group(0), m.group(0)), numbers)

    return _SECTION_RE.sub(replace, text)


class Paragraph(_Paragraph):
    """A Paragraph whose "section N" mentions are live links to those sections.

    Every report module builds its text with this class, so a section
    reference is a link wherever it appears, with no call site having to
    remember to make it one.
    """

    def __init__(
        self,
        text: str,
        style: ParagraphStyle | None = None,
        bulletText: Any = None,  # noqa: N803 - ReportLab's own parameter names
        frags: Any = None,
        caseSensitive: int = 1,  # noqa: N803
        encoding: str = "utf8",
    ) -> None:
        if isinstance(text, str):
            text = link_sections(text)
        super().__init__(text, style, bulletText, frags, caseSensitive, encoding)


# --------------------------------------------------------------------------- #
# Destinations
# --------------------------------------------------------------------------- #


class Anchor(Flowable):
    """A zero-size flowable that names the place it sits.

    Put it in front of the thing a link should land on. With `outline` it
    also adds an entry to the reader's bookmarks panel; the entry is added
    once even if the flowable is drawn more than once.

    `_ZEROSIZE` is ReportLab's own flag for flowables that must be drawn
    despite having no height; containers such as `KeepInFrame` skip
    zero-height content without it, and the destination would silently
    never be defined.
    """

    _ZEROSIZE = True

    def __init__(
        self,
        name: str,
        *,
        above: float = HEADROOM,
        outline: tuple[str, int] | None = None,
    ) -> None:
        super().__init__()
        self.name = name
        self.above = above
        self.outline = outline
        self.width = self.height = 0

    def wrap(self, *_: float) -> tuple[float, float]:
        return 0.0, 0.0

    def draw(self) -> None:
        canv = self.canv
        canv.bookmarkHorizontal(self.name, 0, self.above)
        if self.outline is not None:
            seen = getattr(canv, "_openbiota_outlined", None)
            if seen is None:
                seen = set()
                canv._openbiota_outlined = seen  # noqa: SLF001 - our own marker on the canvas
            if self.name not in seen:
                seen.add(self.name)
                title, level = self.outline
                canv.addOutlineEntry(title, self.name, level=level)


_TAGS: Final = re.compile(r"<[^>]+>")


def _plain(markup: str) -> str:
    return _TAGS.sub("", markup).replace("&nbsp;", " ").replace("&mdash;", "—").replace("&amp;", "&")


def section_heading(
    section: int,
    title: str,
    style: ParagraphStyle,
    *,
    sub: int | None = None,
    opens_section: bool = False,
    trailing: str = "",
) -> list[Flowable]:
    """A section heading with its destination and outline entry.

    Returns the anchor(s) and the heading, to be `extend`ed onto the story.
    The heading text is `"{section} &nbsp; {title}"`, as it has always
    been; sub-sections read `"{section}.{sub}"` and sit one level down in
    the outline. A section that begins straight away with its first
    sub-section passes `opens_section=True`, and the sub-heading defines
    the section's own destination and outline entry as well. `trailing` is
    markup that goes on the heading line but not into the outline — a link
    back to where the reader came from, say.
    """
    by_number, _ = _sections_and_colour()
    key = by_number[section]
    # Section 0 is the unnumbered one (the contents page): title alone, in
    # the heading and in the outline.
    label = "" if section == 0 else (str(section) if sub is None else f"{section}.{sub}")
    lead = f"{label} &nbsp; " if label else ""
    outline_lead = f"{label}  " if label else ""
    out: list[Flowable] = []
    if sub is not None and opens_section:
        out.append(Anchor(section_dest(key), above=6 * mm, outline=(f"{section}  {_plain(title)}", 0)))
    out.append(
        Anchor(
            section_dest(key, sub),
            above=6 * mm,
            outline=(f"{outline_lead}{_plain(title)}", 0 if sub is None else 1),
        )
    )
    out.append(Paragraph(f"{lead}{title}" + (f" &nbsp;&nbsp;{trailing}" if trailing else ""), style))
    return out


def back_link(dest: str, label: str, size: float = 6.4) -> str:
    """Markup for a small link back to where the reader came from."""
    _, colour = _sections_and_colour()
    return f'<a href="#{dest}" color="{colour}"><font size="{size}">\u2039 {label}</font></a>'


def back_link_line(dest: str, label: str, *, size: float = 7.0) -> Flowable:
    """A back-link on its own line, left-aligned, above the title it serves.

    Set inline after a title the link is easy to miss, and it competes with
    the title for the eye. On its own line above, left-aligned to the same
    margin as everything else, it reads as navigation: the reader sees where
    they came from before they read where they are.
    """
    _, colour = _sections_and_colour()
    style = ParagraphStyle(
        "backlink",
        fontName="Helvetica-Bold",
        fontSize=size,
        leading=size + 2.6,
        textColor=HexColor(colour),
        alignment=TA_LEFT,
        spaceBefore=0,
        spaceAfter=1.2,
    )
    return Paragraph(
        f'<a href="#{dest}" color="{colour}">\u2039&nbsp;{label}</a>', style
    )


# --------------------------------------------------------------------------- #
# Table rows as links and destinations
# --------------------------------------------------------------------------- #


class RowLink(Flowable):
    """A zero-size marker in a row's first cell.

    `LinkedTable` reads it when drawing: the row (and the `rows - 1` rows
    under it) becomes a link to `href`, and `anchor` names the row so other
    links can land on it. Kept in the cell rather than on the table so it
    travels with its row when the table splits across pages.
    """

    _ZEROSIZE = True

    def __init__(self, *, href: str | None = None, anchor: str | None = None, rows: int = 1) -> None:
        super().__init__()
        self.href = href
        self.anchor = anchor
        self.rows = rows
        self.width = self.height = 0

    def wrap(self, *_: float) -> tuple[float, float]:
        return 0.0, 0.0

    def draw(self) -> None:
        return


def linked_cell(marker: RowLink, *content: Flowable) -> list[Flowable]:
    """First-cell content for a linked row."""
    return [marker, *content]


def linked_block(href: str, content: list[Flowable], width: float) -> Table:
    """Make a whole stack of flowables one clickable area.

    For section-level elements on page 1 — the age block, the diversity
    block, the gut-health dial, the pathogen card — where the reader should
    be able to click anywhere on the element. A single-cell table with no
    padding, so the block measures exactly as it did before.
    """
    table = LinkedTable([[linked_cell(RowLink(href=href), *content)]], colWidths=[width], hAlign="LEFT")
    table.setStyle(
        [
            ("LEFTPADDING", (0, 0), (-1, -1), 0),
            ("RIGHTPADDING", (0, 0), (-1, -1), 0),
            ("TOPPADDING", (0, 0), (-1, -1), 0),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]
    )
    return table


def _marker(row: Any) -> RowLink | None:
    if not row:
        return None
    first = row[0]
    if isinstance(first, RowLink):
        return first
    if isinstance(first, list | tuple):
        for item in first:
            if isinstance(item, RowLink):
                return item
    return None


class LinkedTable(Table):
    """A Table whose rows can be links and destinations, via `RowLink` markers.

    Drawn exactly as a Table; afterwards each marked row gets a link
    rectangle across the table's full width, and a destination just above
    its top edge.
    """

    def draw(self) -> None:
        super().draw()
        heights = list(self._rowHeights)
        y = self._height
        for index, row in enumerate(self._cellvalues):
            marker = _marker(row)
            if marker is not None:
                span = max(1, min(marker.rows, len(heights) - index))
                bottom = y - sum(heights[index:index + span])
                if marker.href:
                    self.canv.linkRect("", marker.href, (0, bottom, self._width, y), relative=1, thickness=0)
                if marker.anchor:
                    self.canv.bookmarkHorizontal(marker.anchor, 0, y + HEADROOM)
            y -= heights[index]


def self_test() -> int:
    """The text linker: what it links, and what it leaves alone."""
    checks = 0
    # Derived from SECTIONS rather than written as literals, so inserting a
    # section renumbers the fixture instead of breaking it.
    from openbiota.pdfreport import SECTIONS as _S

    age = _S["age"]
    out = link_sections(f"see section {age} for the working")
    assert 'href="#sec-age"' in out and f">section {age}</a>" in out, out
    a, b, c, d = (
        _S["community"], _S["organisms"], _S["actions"], _S["organisms_detail"],
    )
    out = link_sections(f"sections {a} to {b}, {c}, {d}")
    assert out.count("<a ") == 4, out
    assert 'href="#sec-community"' in out, out
    assert 'href="#sec-organisms_detail"' in out, out
    # A sub-section reference is not a section.
    sub = f"section {_S['actions']}.1"
    assert link_sections(f"{sub} lists them") == f"{sub} lists them"
    # Numbers that are not sections are left alone.
    assert link_sections("section 99 does not exist") == "section 99 does not exist"
    # Already linked text is not linked twice.
    path = _S["pathogens"]
    once = link_sections(f"Section {path} gives the amount")
    assert link_sections(once) == once
    # Bold and &nbsp; do not get in the way.
    out = link_sections(f"<b>Section&nbsp;{path}</b> gives")
    assert 'href="#sec-pathogens"' in out, out
    checks += 6
    assert function_dest("Butyrate (SCFA)") == "fn-butyrate-scfa"
    assert function_dest("x", detail=True).endswith("-detail")
    assert profile_dest("ibs_d") == "pr-ibs-d"
    checks += 3
    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} link checks passed")
