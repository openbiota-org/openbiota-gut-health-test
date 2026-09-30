"""Read a Tiny Health gut report PDF into structured metrics.

Tiny Health is a CLIA-lab consumer microbiome test that ran shotgun
sequencing on the same stool specimens this repository screens, so their
reports are the one external check available on real samples: two
independent pipelines, one set of FASTQs. They are not ground truth —
this module exists so that where the two disagree, the disagreement can
be *located* rather than argued about.

Three report templates are in the corpus and all three are parsed:

    v4.8.3  metric name at x≈52, value at x≈220, evaluation encoded in the
            value's *colour* (green great, grey okay, amber improve, red
            support), thresholds in 6.8 pt grey at the right
    v5.1.0   as above but the evaluation is a word at x≈530
    v5.3.1   numbered rows, name at x≈80, evaluation word at x≈318, value
            in a white-on-colour pill at x≥375 or as a plain 9.8 pt figure

The parsers work from text position and colour rather than reading order,
because reading order interleaves the threshold ticks of one row with the
label of the next. Each returns the same `Metric` records, so the
comparison code never has to know which template it came from.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

#: Evaluation words, normalised. Tiny Health renamed these between
#: templates ("Needs support" -> "Support") without changing the meaning.
EVAL_WORDS: Final[dict[str, str]] = {
    "Support": "support", "Needs support": "support",
    "Improve": "improve", "Needs improving": "improve",
    "Okay": "okay", "Great": "great",
}

#: v4.8.3 encodes the evaluation only in the colour of the value.
COLOUR_EVAL: Final[dict[int, str | None]] = {
    0x009985: "great", 0x707269: "okay", 0xEDA831: "improve",
    0xDA5D5D: "support", 0x3B3941: None,
}

_RANK: Final = re.compile(r"\s*\((p|f|g|sp)\)\s*$")
_NUM: Final = re.compile(r"^([\d.,]+)\s*(%|rpkm|species|years|\w+)?$")
_VERSION: Final = re.compile(r"(?:System version:|v)\s*(\d+\.\d+\.\d+)")
#: The headline score, written "72 /100" (v5.x) or "87 Microbiome summary
#: score" (v4.8.3). The scale's own axis labels read "0 ... 100 Microbiome
#: summary score", so the "/100" form is tried first and the bare form is
#: only accepted when it is not the axis maximum.
_SCORE_SLASH: Final = re.compile(r"(\d+)\s*/\s*100\b")
_SCORE_WORDS: Final = re.compile(r"(?<!0 )(\d+) Microbiome summary score")


@dataclass(frozen=True)
class Metric:
    """One reading from a Tiny Health report."""

    metric: str
    category: str | None
    subcategory: str | None
    evaluation: str | None
    value: float | None
    unit: str | None
    rank: str | None
    page: int
    thresholds: tuple[float, ...] = ()

    @property
    def direction(self) -> str:
        """low / mid / high against Tiny Health's own band thresholds.

        Their bands are drawn as tick marks on a scale, so the thresholds
        are recoverable but their *meaning* (which side is good) is not:
        that lives in the evaluation word. This is the position only.
        """
        if self.value is None or not self.thresholds:
            return "unknown"
        lo, hi = min(self.thresholds), max(self.thresholds)
        if self.value < lo:
            return "low"
        if self.value > hi:
            return "high"
        return "mid"


@dataclass(frozen=True)
class Report:
    sample: str
    version: str
    pages: int
    summary_score: int | None
    metrics: tuple[Metric, ...] = field(default_factory=tuple)

    def by_name(self, name: str) -> Metric | None:
        return next((m for m in self.metrics if m.metric == name), None)


def _lines(page: Any) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for block in page.get_text("dict")["blocks"]:
        if block["type"] != 0:
            continue
        for line in block["lines"]:
            text = "".join(s["text"] for s in line["spans"]).strip().replace("\xa0", " ")
            if not text:
                continue
            span = line["spans"][0]
            out.append({
                "y": line["bbox"][1], "x": line["bbox"][0],
                "size": round(span["size"], 1), "colour": span.get("color", 0), "text": text,
            })
    out.sort(key=lambda r: (r["y"], r["x"]))
    return out


def _number(text: str) -> tuple[float | None, str | None]:
    match = _NUM.match(text.replace(",", ""))
    if not match:
        return None, text or None
    return float(match.group(1)), (match.group(2) or "").strip() or None


def _strip_rank(name: str) -> tuple[str, str | None]:
    match = _RANK.search(name)
    return (_RANK.sub("", name), match.group(1)) if match else (name, None)


def _thresholds(rows: list[dict[str, Any]]) -> tuple[float, ...]:
    return tuple(v for v in (_number(r["text"])[0] for r in rows) if v is not None)


def _parse_v53(doc: Any) -> list[Metric]:
    out: list[Metric] = []
    category = subcategory = None
    for page_no, page in enumerate(doc, start=1):
        lines = _lines(page)
        for r in lines:
            if (r["size"] == 13.5 and r["x"] < 40 and not r["text"].endswith("metrics")
                    and not r["text"][:1].isdigit()):
                category = r["text"]
        numbered = [r for r in lines
                    if r["size"] == 7.9 and 42 <= r["x"] <= 56 and r["text"].isdigit()]
        for row in numbered:
            y = row["y"]
            name = next((r for r in lines if r["size"] == 7.9 and 78 <= r["x"] <= 82
                         and abs(r["y"] - y) < 1.5), None)
            if name is None:
                continue
            subs = [r for r in lines if r["size"] == 7.9 and 32 <= r["x"] <= 37 and r["y"] < y]
            if subs:
                subcategory = subs[-1]["text"]
            ev = next((r for r in lines if r["size"] == 7.0 and 300 <= r["x"] <= 345
                       and abs(r["y"] - y) < 6 and r["text"] in EVAL_WORDS), None)
            # White-on-colour pill, or a plain figure where there is no scale.
            val = next((r for r in lines if r["size"] in (7.5, 7.0) and r["x"] >= 375
                        and abs(r["y"] - y) < 8 and r["colour"] == 0xFFFFFF
                        and _number(r["text"])[0] is not None), None)
            if val is None:
                val = next((r for r in lines if r["size"] == 9.8 and r["x"] >= 375
                            and abs(r["y"] - y) < 4 and _number(r["text"])[0] is not None), None)
            value, unit = _number(val["text"]) if val else (None, None)
            ticks = [r for r in lines if r["size"] == 6.0 and r["x"] >= 380 and 18 <= (r["y"] - y) <= 30]
            out.append(Metric(
                metric=name["text"], category=category, subcategory=subcategory,
                evaluation=EVAL_WORDS.get(ev["text"]) if ev else None,
                value=value, unit=unit, rank=None, page=page_no,
                thresholds=_thresholds(ticks),
            ))
    return out


def _parse_v48(doc: Any) -> list[Metric]:
    out: list[Metric] = []
    category = None
    for page_no, page in enumerate(doc, start=1):
        lines = _lines(page)
        for r in lines:
            if r["size"] == 13.5 and r["x"] < 45:
                category = r["text"]
        for name in [r for r in lines if r["size"] == 9.0 and 50 <= r["x"] <= 54]:
            val = next((r for r in lines if r["size"] == 9.0 and 215 <= r["x"] <= 225
                        and abs(r["y"] - name["y"]) < 1.5), None)
            if val is None:
                continue
            subs = [r for r in lines if r["size"] == 9.0 and 38 <= r["x"] <= 41 and r["y"] < name["y"]]
            ticks = [r for r in lines if r["size"] == 6.8 and r["x"] >= 460
                     and -14 <= (r["y"] - name["y"]) <= 14]
            label, rank = _strip_rank(name["text"])
            value, unit = _number(val["text"])
            out.append(Metric(
                metric=label, category=category,
                subcategory=subs[-1]["text"] if subs else None,
                evaluation=COLOUR_EVAL.get(val["colour"], "unknown"),
                value=value, unit=unit, rank=rank, page=page_no,
                thresholds=_thresholds(ticks),
            ))
    return out


def _parse_v51(doc: Any) -> list[Metric]:
    out: list[Metric] = []
    category = None
    for page_no, page in enumerate(doc, start=1):
        lines = _lines(page)
        for r in lines:
            if r["size"] == 13.5 and r["x"] < 40:
                category = r["text"]
        for ev in [r for r in lines if r["size"] == 9.0 and r["x"] >= 530 and r["text"] in EVAL_WORDS]:
            name = next((r for r in lines if r["size"] == 9.0 and r["x"] <= 32
                         and r["colour"] == 0x3B3942 and -2 <= (ev["y"] - r["y"]) <= 10), None)
            if name is None:
                continue
            subs = [r for r in lines if r["size"] == 10.5 and 35 <= r["x"] <= 40 and r["y"] < name["y"]]
            val = next((r for r in lines if r["size"] == 9.0 and 295 <= r["x"] <= 430
                        and -14 <= (r["y"] - name["y"]) <= 2 and _number(r["text"])[0] is not None), None)
            ticks = [r for r in lines if r["size"] == 6.8 and 320 <= r["x"] <= 545
                     and -3 <= (r["y"] - name["y"]) <= 22]
            label, rank = _strip_rank(name["text"])
            value, unit = _number(val["text"]) if val else (None, None)
            out.append(Metric(
                metric=label, category=category,
                subcategory=subs[-1]["text"] if subs else None,
                evaluation=EVAL_WORDS[ev["text"]], value=value, unit=unit,
                rank=rank, page=page_no, thresholds=_thresholds(ticks),
            ))
    return out


def parse(path: Path) -> Report:
    """Read one Tiny Health report."""
    import pymupdf

    doc = pymupdf.open(path)
    text = "\n".join(page.get_text() for page in doc)
    version = (_VERSION.search(text).group(1) if _VERSION.search(text) else "unknown")
    flat = " ".join(text.split())
    match = _SCORE_SLASH.search(flat) or _SCORE_WORDS.search(flat)
    score = int(match.group(1)) if match else None
    # 5.5 keeps 5.3's table layout, so it reads with the same parser: 126
    # metrics against 5.3's 139. An unrecognised 5.x fell through to the 4.8
    # reader and produced nothing at all, which is worse than a version guess
    # — the concordance tests then compare against an empty report.
    if version.startswith(("5.3", "5.4", "5.5", "5.6")):
        metrics = _parse_v53(doc)
    elif version.startswith("5.1"):
        metrics = _parse_v51(doc)
    else:
        metrics = _parse_v48(doc)
    return Report(
        sample=path.name.split(".")[0], version=version, pages=len(doc),
        summary_score=score, metrics=tuple(metrics),
    )


def load_all(directory: Path) -> dict[str, Report]:
    """Every Tiny Health report in a directory, keyed by sample name."""
    return {p.name.split(".")[0]: parse(p) for p in sorted(directory.glob("*.pdf"))}


if __name__ == "__main__":  # pragma: no cover
    import sys
    from collections import Counter

    root = Path(sys.argv[1] if len(sys.argv) > 1 else "tinyhealth")
    for name, report in load_all(root).items():
        counts = Counter(m.evaluation for m in report.metrics)
        print(f"{name}: v{report.version} score={report.summary_score} "
              f"{len(report.metrics)} metrics, "
              f"{sum(1 for m in report.metrics if m.value is not None)} with values, "
              f"{dict(counts)}")
