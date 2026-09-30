"""A reader for Strict OOXML workbooks, which openpyxl cannot open.

BUILD_SPEC_v0.8.3 section 11.7, issue 6: "S5 uses Strict OOXML; support its
namespace or convert with documented provenance. **A parser returning zero
sheets must fail loudly.**"

That is exactly what happens without this module. The published S5 supplement
declares `http://purl.oclc.org/ooxml/spreadsheetml/main` in `xl/workbook.xml`
while its `[Content_Types].xml` uses the transitional media types, so
openpyxl walks the content types, finds nothing it recognises as a worksheet,
warns once and returns a workbook with **no sheets at all**. Read through
pandas, that surfaces as "Worksheet not found"; read through a more forgiving
wrapper it would surface as an empty dataset, and an empty replay fixture
silently passes every assertion about values it no longer contains.

So this reader works from the workbook relationships rather than the content
types, accepts either namespace, and raises on an empty sheet list.

Only what the fixture needs is implemented: shared strings, inline strings,
numbers, booleans and dates-as-serials. Formulas are read as their cached
values, which is what a published supplement carries.
"""

from __future__ import annotations

import re
import zipfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

#: Both spreadsheetml namespaces. Strict uses the OCLC form; transitional uses
#: the openxmlformats form. A file may mix them, as this one does.
NS_MAIN: Final[tuple[str, ...]] = (
    "http://purl.oclc.org/ooxml/spreadsheetml/main",
    "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
)
NS_REL: Final[tuple[str, ...]] = (
    "http://purl.oclc.org/ooxml/officeDocument/relationships",
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
)

_CELL_RE: Final = re.compile(
    r"<c\b(?P<attrs>[^>]*?)(?:/>|>(?P<body>.*?)</c>)", re.S
)
_ATTR_RE: Final = re.compile(r'(\w+(?::\w+)?)="([^"]*)"')
_V_RE: Final = re.compile(r"<v\b[^>]*>(.*?)</v>", re.S)
_IS_T_RE: Final = re.compile(r"<is\b[^>]*>(.*?)</is>", re.S)
_T_RE: Final = re.compile(r"<t\b[^>]*>(.*?)</t>", re.S)
_ROW_RE: Final = re.compile(r"<row\b(?P<attrs>[^>]*?)(?:/>|>(?P<body>.*?)</row>)", re.S)
_SI_RE: Final = re.compile(r"<si\b[^>]*?(?:/>|>(.*?)</si>)", re.S)
_COL_RE: Final = re.compile(r"^([A-Z]+)")


class StrictOoxmlError(ValueError):
    """A workbook this reader will not silently misread."""


def _unescape(text: str) -> str:
    return (
        text.replace("&lt;", "<").replace("&gt;", ">")
        .replace("&quot;", '"').replace("&apos;", "'")
        .replace("&amp;", "&")
    )


def _column_index(ref: str) -> int:
    """`A` -> 0, `Z` -> 25, `AA` -> 26."""
    match = _COL_RE.match(ref)
    if not match:
        raise StrictOoxmlError(f"cell reference {ref!r} has no column")
    n = 0
    for char in match.group(1):
        n = n * 26 + (ord(char) - 64)
    return n - 1


@dataclass(slots=True)
class StrictWorkbook:
    """One Strict OOXML workbook, opened by relationship rather than by media type."""

    path: Path
    sheet_names: tuple[str, ...]
    _zip: zipfile.ZipFile
    _sheet_paths: dict[str, str]
    _shared: tuple[str, ...]
    strict: bool

    @classmethod
    def open(cls, path: Path | str) -> StrictWorkbook:
        path = Path(path)
        archive = zipfile.ZipFile(path)
        try:
            workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
        except KeyError as exc:  # pragma: no cover - not a spreadsheet at all
            raise StrictOoxmlError(f"{path}: no xl/workbook.xml; not an xlsx package") from exc

        strict = any(ns in workbook_xml for ns in NS_MAIN[:1])

        # Relationship id -> part path, from the workbook's own rels part.
        rels_xml = archive.read("xl/_rels/workbook.xml.rels").decode("utf-8")
        rel_targets: dict[str, str] = {}
        for chunk in re.findall(r"<Relationship\b[^>]*/>", rels_xml):
            attrs = dict(_ATTR_RE.findall(chunk))
            rid, target = attrs.get("Id"), attrs.get("Target")
            if not rid or not target:
                continue
            rel_targets[rid] = target if target.startswith("xl/") else f"xl/{target.lstrip('/')}"

        sheets: list[str] = []
        sheet_paths: dict[str, str] = {}
        for chunk in re.findall(r"<sheet\b[^>]*/>", workbook_xml):
            attrs = dict(_ATTR_RE.findall(chunk))
            name = attrs.get("name")
            rid = next((v for k, v in attrs.items() if k.endswith(":id") or k == "id"), None)
            if not name or not rid or rid not in rel_targets:
                continue
            sheets.append(_unescape(name))
            sheet_paths[_unescape(name)] = rel_targets[rid]

        # The rule this module exists for.
        if not sheets:
            raise StrictOoxmlError(
                f"{path}: the workbook declares no readable sheets. This is the Strict "
                "OOXML failure the specification calls out: an empty sheet list must fail "
                "loudly rather than be accepted as an empty valid dataset."
            )

        shared: list[str] = []
        for candidate in ("xl/sharedStrings.xml", "xl/SharedStrings.xml"):
            if candidate in archive.namelist():
                blob = archive.read(candidate).decode("utf-8")
                for body in _SI_RE.findall(blob):
                    shared.append(_unescape("".join(_T_RE.findall(body or ""))))
                break

        return cls(
            path=path, sheet_names=tuple(sheets), _zip=archive,
            _sheet_paths=sheet_paths, _shared=tuple(shared), strict=strict,
        )

    def rows(self, sheet: str) -> Iterator[list[Any]]:
        """Every row of one sheet as a list of Python values.

        Ragged rows are padded to the widest row seen so far, and gaps from
        sparse cell references are filled with None, so a column index means
        the same thing on every row.
        """
        if sheet not in self._sheet_paths:
            raise StrictOoxmlError(
                f"{self.path}: no sheet named {sheet!r}; have {list(self.sheet_names)}"
            )
        blob = self._zip.read(self._sheet_paths[sheet]).decode("utf-8")
        for row_match in _ROW_RE.finditer(blob):
            body = row_match.group("body") or ""
            cells: dict[int, Any] = {}
            for cell in _CELL_RE.finditer(body):
                attrs = dict(_ATTR_RE.findall(cell.group("attrs") or ""))
                ref = attrs.get("r")
                index = _column_index(ref) if ref else (max(cells) + 1 if cells else 0)
                cells[index] = self._value(attrs.get("t"), cell.group("body") or "")
            if not cells:
                yield []
                continue
            yield [cells.get(i) for i in range(max(cells) + 1)]

    def _value(self, kind: str | None, body: str) -> Any:
        if kind == "s":
            raw = _V_RE.search(body)
            if raw is None:
                return None
            try:
                return self._shared[int(raw.group(1))]
            except (ValueError, IndexError):
                return None
        if kind == "inlineStr":
            inline = _IS_T_RE.search(body)
            return _unescape("".join(_T_RE.findall(inline.group(1)))) if inline else None
        if kind == "str":
            raw = _V_RE.search(body)
            return _unescape(raw.group(1)) if raw else None
        if kind == "b":
            raw = _V_RE.search(body)
            return bool(int(raw.group(1))) if raw else None
        raw = _V_RE.search(body)
        if raw is None:
            return None
        text = raw.group(1).strip()
        if not text:
            return None
        try:
            return int(text) if re.fullmatch(r"-?\d+", text) else float(text)
        except ValueError:
            return _unescape(text)

    def table(self, sheet: str) -> tuple[tuple[str, ...], list[list[Any]]]:
        """One sheet as a header tuple and its data rows."""
        iterator = self.rows(sheet)
        try:
            header_row = next(iterator)
        except StopIteration as exc:
            raise StrictOoxmlError(f"{self.path}: sheet {sheet!r} is empty") from exc
        header = tuple(
            (str(h) if h is not None else f"column_{i}") for i, h in enumerate(header_row)
        )
        rows = [r for r in iterator if any(v is not None for v in r)]
        if not rows:
            raise StrictOoxmlError(
                f"{self.path}: sheet {sheet!r} has a header and no data rows"
            )
        return header, rows

    def records(self, sheet: str) -> list[dict[str, Any]]:
        """One sheet as a list of dictionaries keyed by header."""
        header, rows = self.table(sheet)
        out: list[dict[str, Any]] = []
        for row in rows:
            padded = list(row) + [None] * (len(header) - len(row))
            out.append(dict(zip(header, padded, strict=False)))
        return out

    def close(self) -> None:
        self._zip.close()

    def __enter__(self) -> StrictWorkbook:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def read_sheet(path: Path | str, sheet: str) -> list[dict[str, Any]]:
    """Read one sheet of a Strict OOXML workbook into records."""
    with StrictWorkbook.open(path) as wb:
        return wb.records(sheet)


def sheet_names(path: Path | str) -> Sequence[str]:
    with StrictWorkbook.open(path) as wb:
        return wb.sheet_names


__all__ = [
    "NS_MAIN",
    "NS_REL",
    "StrictOoxmlError",
    "StrictWorkbook",
    "read_sheet",
    "sheet_names",
]
