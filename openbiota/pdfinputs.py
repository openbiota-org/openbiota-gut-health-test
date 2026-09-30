"""The input register section — Your Information & Report Context (A14, spec §10.1).

One section near the end of the report, holding everything that went into it:
what you supplied, what the report assumed because you did not, what it
recorded but does not use, and which optional measurements exist. The
domain sections show what each input *changed*; this is the register and
the provenance view.

It opens with a compact summary, then the table the specification names:
Information · Value used · Supplied or assumed · What it affects. Absent
optional inputs are grouped into categories rather than given a page of
empty charts each.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final

from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import CondPageBreak, Spacer

from openbiota.extension.inputs import CATEGORIES
from openbiota.pdflinks import Paragraph, section_heading
from openbiota.pdfparts import DASH
from openbiota.pdfparts import hex_of as _hex
from openbiota.pdfparts import table as _table
from openbiota.pdfreport import AMBER, CONTENT_WIDTH, GREEN, INK_FAINT, INK_SOFT, SLATE

#: How each provenance state is shown. Colour never carries the meaning on
#: its own: the words say it too.
STATE_TONE: Final[Mapping[str, Any]] = {
    "supplied": GREEN,
    "derived": GREEN,
    "inherited_default": AMBER,
    "assumed_default": AMBER,
    "not_supplied": INK_FAINT,
    "not_applicable": INK_FAINT,
    "conflicting": AMBER,
}


def _value_cell(record: Mapping[str, Any], st: Mapping[str, ParagraphStyle]) -> Paragraph:
    value = record.get("effective_value")
    if value is None:
        return Paragraph(f"<font color='{_hex(INK_FAINT)}'>{DASH}</font>", st["cell"])
    text = ("yes" if value else "no") if isinstance(value, bool) else str(value)
    unit = record.get("unit")
    return Paragraph(f"<b>{text}</b>" + (f" {unit}" if unit else ""), st["cell"])


def _affects_cell(record: Mapping[str, Any], st: Mapping[str, ParagraphStyle]) -> Paragraph:
    if not record.get("used_in_this_report"):
        note = record.get("usage_note") or "Recorded; not used in this report"
        limits = record.get("limitations") or []
        extra = f" {limits[0]}" if limits else ""
        return Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>{note}.{extra}</font>", st["cell"])
    parts = []
    for influence in record.get("affects") or []:
        if influence.get("kind") == "none":
            continue
        verb = {
            "changes_calculation": "changes a number",
            "changes_interpretation": "changes how a reading is read",
            "adds_context": "adds context",
        }.get(str(influence.get("kind")), "is used")
        parts.append(f"<b>{verb}</b> \u2014 {influence.get('explanation')}")
    return Paragraph(f"<font size='6.4'>{' '.join(parts)}</font>", st["cell"])


def inputs_section(
    story: list[Any],
    st: Mapping[str, ParagraphStyle],
    register: Mapping[str, Any] | None,
    *,
    section: int | None = None,
) -> None:
    """Render the whole register. Its display number comes from `SECTIONS`."""
    if not register:
        return
    summary = register.get("summary") or {}
    heading = "Your Information & Report Context"
    if section is None:
        story.append(CondPageBreak(70 * mm))
        story.append(Spacer(1, 3 * mm))
        story.append(Paragraph(heading, st["h2"]))
    else:
        # The same helper every other section uses, so this one appears in
        # the bookmarks panel and can be linked to from the contents.
        story.extend(section_heading(section, heading, st["h1"]))
    story.append(Paragraph(register.get("headline") or "", st["body"]))

    # -- compact summary --------------------------------------------------- #
    counts = [
        ("Supplied", summary.get("n_supplied", 0), GREEN),
        ("Assumed", summary.get("n_assumed", 0), AMBER),
        ("Not supplied", summary.get("n_not_supplied", 0), SLATE),
        ("Laboratory results", summary.get("n_laboratory_results", 0), SLATE),
    ]
    if summary.get("n_conflicting"):
        counts.append(("Conflicting", summary["n_conflicting"], AMBER))
    story.append(Paragraph(
        " &nbsp;&nbsp;".join(
            f"<font color='{_hex(tone)}'><b>{value}</b></font> "
            f"<font size='6.6' color='{_hex(INK_SOFT)}'>{label.lower()}</font>"
            for label, value, tone in counts
        ), st["small"],
    ))

    _register_table(story, st, register)
    _missing_block(story, st, register)
    _laboratory_block(story, st, register)

    story.append(Spacer(1, 1.5 * mm))
    for line in register.get("contract") or []:
        story.append(Paragraph(
            f"<font size='6.2' color='{_hex(INK_FAINT)}'>\u00b7 {line}</font>", st["small"]))


def _register_table(
    story: list[Any], st: Mapping[str, ParagraphStyle], register: Mapping[str, Any]
) -> None:
    """Information · Value used · Supplied or assumed · What it affects."""
    groups = [
        g for g in register.get("by_category") or []
        if any(r.get("input_state") != "not_supplied" for r in g.get("records") or [])
    ]
    if not groups:
        return
    story.append(Spacer(1, 2.5 * mm))
    story.append(Paragraph("<b>Everything this report used</b>", st["h3"]))
    head = ("INFORMATION", "VALUE USED", "SUPPLIED OR ASSUMED", "WHAT IT AFFECTS")
    rows: list[list[Any]] = [[Paragraph(h, st["label"]) for h in head]]
    order = {key: i for i, (key, _) in enumerate(CATEGORIES)}
    for group in sorted(groups, key=lambda g: order.get(str(g.get("category")), 99)):
        present = [r for r in group.get("records") or []
                   if r.get("input_state") != "not_supplied"]
        if not present:
            continue
        rows.append([
            Paragraph(f"<font size='6.4' color='{_hex(INK_SOFT)}'><b>"
                      f"{str(group.get('category_label')).upper()}</b></font>", st["cell"]),
            "", "", "",
        ])
        # Defaults that all say the same thing are collapsed into one row: a
        # reader needs to know thirteen safety questions were assumed "no",
        # not to read thirteen identical lines.
        for record in _collapse(present):
            state = str(record.get("input_state"))
            rows.append([
                Paragraph(record.get("label") or record.get("input_id"), st["cell"]),
                _value_cell(record, st),
                Paragraph(
                    f"<font color='{_hex(STATE_TONE.get(state, SLATE))}'>"
                    f"{_provenance_words(state)}</font>"
                    + (f"<br/><font size='5.8' color='{_hex(INK_SOFT)}'>"
                       f"{record.get('default_rationale')}</font>"
                       if record.get("default_rationale") else ""),
                    st["cell"],
                ),
                _affects_cell(record, st),
            ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.14, 0.26, 0.38)]))


def _provenance_words(state: str) -> str:
    return {
        "supplied": "Supplied",
        "derived": "Worked out from what you supplied",
        "inherited_default": "Assumed \u2014 existing default",
        "assumed_default": "Assumed for this run",
        "not_supplied": "Not supplied",
        "not_applicable": "Does not apply",
        "conflicting": "Conflicting records",
    }.get(state, state)


def _collapse(records: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Fold a run of identical safety defaults into one honest row."""
    safety = [
        r for r in records
        if r.get("input_state") == "inherited_default"
        and r.get("effective_value") is False
        and str(r.get("category")) == "clinical_history"
    ]
    if len(safety) < 4:
        return list(records)
    rest = [r for r in records if r not in safety]
    names = ", ".join(str(r.get("label")) for r in safety)
    folded = {
        "input_id": "context.safety_flags",
        "label": f"{len(safety)} safety questions",
        "effective_value": "assumed no",
        "input_state": "inherited_default",
        "category": "clinical_history",
        "default_rationale": (
            "None were supplied. Each is assumed absent so the action cards can be produced; "
            "any one of them being true changes them."
        ),
        "used_in_this_report": True,
        "affects": [{
            "kind": "changes_interpretation",
            "explanation": (
                f"These are the safety rules on the intervention cards ({names}). If any is "
                "true, the affected probiotic, supplement, diet or drug card moves to "
                "clinician review or is withdrawn."
            ),
        }],
    }
    return [*rest, folded]


def _missing_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], register: Mapping[str, Any]
) -> None:
    """Optional inputs nobody supplied, grouped, with what each would add."""
    missing = [
        m for m in register.get("missing_categories") or []
        if m.get("category") != "laboratory"
    ]
    if not missing:
        return
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph("<b>Not supplied, and what each would add</b> &nbsp;"
                           f"<font size='7' color='{_hex(INK_FAINT)}'>every reading above was "
                           "produced without them</font>", st["h3"]))
    rows: list[list[Any]] = [[Paragraph(h, st["label"])
                              for h in ("CATEGORY", "NOT SUPPLIED", "WHAT IT WOULD ADD")]]
    for entry in missing:
        examples = ", ".join(entry.get("examples") or [])
        rows.append([
            Paragraph(f"<b>{entry.get('category_label')}</b>", st["cell"]),
            Paragraph(f"<font size='6.4'>{examples}</font>", st["cell"]),
            Paragraph(f"<font size='6.4'>{entry.get('what_it_would_add') or DASH}</font>",
                      st["cell"]),
        ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.22, 0.34, 0.44)]))


def _laboratory_block(
    story: list[Any], st: Mapping[str, ParagraphStyle], register: Mapping[str, Any]
) -> None:
    """The optional laboratory panels: what was supplied, and where it goes."""
    panels = register.get("laboratory_panels") or []
    if not panels:
        return
    supplied = [p for p in panels if p.get("supplied")]
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        "<b>External laboratory results</b> &nbsp;"
        f"<font size='7' color='{_hex(INK_FAINT)}'>"
        + ("optional; this report is complete without them" if not supplied
           else f"{sum(p.get('n_results', 0) for p in supplied)} imported")
        + "</font>", st["h3"],
    ))
    if not supplied:
        names = ", ".join(str(p.get("label")).lower() for p in panels)
        story.append(Paragraph(
            f"<font size='6.8'><b>Not supplied.</b> Every reading in this report is produced "
            f"from your sequencing data. If you have results for {names}, they can be imported "
            "and will be shown beside the finding each one relates to \u2014 as a separate, "
            "differently-labelled measurement, never merged into a genetic reading or used to "
            "overrule it. Nothing here is missing without them.</font>", st["small"],
        ))
        return
    rows: list[list[Any]] = [[Paragraph(h, st["label"])
                              for h in ("RESULT", "VALUE", "METHOD & SPECIMEN",
                                        "HOW IT RELATES TO THIS REPORT")]]
    for panel in supplied:
        for result in panel.get("results") or []:
            interval = result.get("reference_interval") or {}
            within = result.get("within_laboratory_interval")
            band = ""
            if interval.get("low") is not None or interval.get("high") is not None:
                band = (f"<br/><font size='5.8' color='{_hex(INK_SOFT)}'>laboratory interval "
                        f"{interval.get('low', DASH)}\u2013{interval.get('high', DASH)}"
                        + ("" if within is None else
                           (", within" if within else ", outside")) + "</font>")
            rows.append([
                Paragraph(f"<b>{result.get('original_name')}</b>", st["cell"]),
                Paragraph(f"<b>{result.get('display_value')}</b>{band}", st["cell"]),
                Paragraph(f"<font size='6.4'>{result.get('method') or DASH}; "
                          f"{result.get('specimen')}"
                          + (f"; {result.get('collected_at')}" if result.get("collected_at") else "")
                          + "</font>", st["cell"]),
                Paragraph(f"<font size='6.4'>{result.get('relation_to_this_report')}</font>",
                          st["cell"]),
            ])
    story.append(_table(rows, [CONTENT_WIDTH * w for w in (0.24, 0.18, 0.22, 0.36)]))
    absent = [str(p.get("label")).lower() for p in panels if not p.get("supplied")]
    if absent:
        story.append(Paragraph(
            f"<font size='6.4' color='{_hex(INK_SOFT)}'>Not supplied: {', '.join(absent)}. "
            "No value is assumed for an assay nobody ran, and no normal or abnormal "
            "classification is produced for one.</font>", st["small"],
        ))


__all__ = ["inputs_section"]
