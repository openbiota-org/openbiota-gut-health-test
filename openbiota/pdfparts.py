"""Shared rendering primitives for the v0.8.3 report sections.

Three modules render extension sections and all three need the same four
things: a colour as hex for inline markup, a number that becomes an em dash
rather than a zero when it is absent, a share that distinguishes "none" from
"too small to round", and the thin-ruled table the sections use throughout.

They live here rather than being copied, so that a change to how an absent
value is printed changes it everywhere at once.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Final

from reportlab.lib import colors
from reportlab.platypus import Table, TableStyle

from openbiota.pdfreport import RULE

#: An absent value. Never a zero: they mean different things.
DASH: Final = "\u2014"


def hex_of(colour: colors.Color) -> str:
    """A reportlab colour as `#RRGGBB`, for inline `<font color=...>`."""
    return f"#{int(colour.red * 255):02X}{int(colour.green * 255):02X}{int(colour.blue * 255):02X}"


def num(value: float | None, digits: int = 2, *, suffix: str = "") -> str:
    """A number, or an em dash. Never a zero standing in for no answer."""
    if value is None:
        return DASH
    return f"{value:,.{digits}f}{suffix}"


def pct(value: float | None) -> str:
    """A share, distinguishing none, too-small-to-round, and a percentage."""
    if value is None:
        return DASH
    if value <= 1e-9:
        return "none"
    if value < 0.01:
        return "<1%"
    return f"{value * 100:.0f}%"


def table(
    rows: Sequence[Sequence[Any]], widths: Sequence[float], *, header: bool = True,
    linked: bool = False,
) -> Table:
    """The thin-ruled table the extension sections use.

    With ``linked`` the table is a `LinkedTable`, which turns `RowLink`
    markers in a row's first cell into a link rectangle and a destination.
    Built here rather than by the caller so a linked table and a plain one
    cannot drift apart in padding or rule weight.
    """
    cls: type[Table] = Table
    if linked:
        from openbiota.pdflinks import LinkedTable  # noqa: PLC0415 - avoids a cycle

        cls = LinkedTable
    out = cls([list(r) for r in rows], colWidths=list(widths), hAlign="LEFT",
              repeatRows=1 if header else 0)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.2),
        ("TOPPADDING", (0, 0), (-1, -1), 2.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.25, RULE),
    ]
    if header:
        style.append(("LINEBELOW", (0, 0), (-1, 0), 0.6, RULE))
    out.setStyle(TableStyle(style))
    return out


__all__ = ["DASH", "hex_of", "num", "pct", "table"]
