#!/usr/bin/env python
"""Write the v0.8.3 functional panels from `extension/new_panels.yaml`.

The discrimination rules for these readings were written first, in
`openbiota.extension.{substrates,fermentation,vitamins,nitrogen,biotransform}`.
Rules with nothing to grade report `not_assayed`, which is honest and
useless. These panels give them something to grade.

The content lives in YAML rather than in this file so that a reviewer reads
the biology without reading the generator, and so that adding a panel is an
edit to data rather than to code.

    .venv/bin/python scripts/write_v083_panels.py
    .venv/bin/openbiota build-db
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import yaml

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "extension" / "new_panels.yaml"
PANELS = REPO / "panels"

#: Identity floor for a target. Below this a match is a fold-level
#: resemblance rather than the same enzyme.
MIN_IDENTITY = 55.0
DECOY_IDENTITY = 50.0
LENGTH_TOLERANCE = 0.30


def wrap(text: str, width: int) -> list[str]:
    """Fold to a width, on word boundaries."""
    words = " ".join(str(text).split()).split(" ")
    out: list[str] = []
    line = ""
    for word in words:
        if line and len(line) + len(word) + 1 > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    return out


def scalar(text: str) -> str:
    """Quote a YAML scalar safely whichever quotes it contains.

    A protein name can carry an apostrophe - "pyridoxal 5'-phosphate" - which
    ends a single-quoted scalar in the middle of the query.
    """
    if "'" in text:
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return f"'{text}'"


def query_for(name: str) -> str:
    return f'(protein_name:"{name}") AND (taxonomy_id:2)'


def block(key: str, text: str, indent: int, width: int = 72) -> list[str]:
    pad = " " * indent
    return [f"{pad}{key}: >-", *[f"{pad}  {line}" for line in wrap(text, width)]]


def render(panel: dict[str, Any]) -> str:
    out: list[str] = [f"name: {panel['name']}"]
    out += block("description", panel["description"], 0)
    out.append(f"metabolite: {panel['metabolite']}")
    out += block("pathway", panel["pathway"], 0)
    out.append(f"category: {panel['category']}")
    out += ["", "extension: true", "", "targets:"]
    for spec in panel["targets"]:
        count = int(spec.get("n_uniprot") or 1000)
        # A strain-specific reaction is measured against the exact enzymes
        # described for it, so its query, identity floor and set size are
        # declared on the target rather than taken from the panel default.
        # A protein with no UniProt entry is held in this repository and
        # addressed by file and key; the taxonomy filter would be meaningless
        # on a named sequence, so it is not applied.
        literal = spec.get("query_literal")
        raw = spec.get("query_raw")
        if literal:
            query = literal
        elif raw:
            query = f"{raw} AND (taxonomy_id:2)"
        else:
            query = query_for(spec["query"])
        identity = float(spec.get("min_identity") or MIN_IDENTITY)
        cap = int(spec.get("max_sequences") or min(3000, max(200, count)))
        out += [
            f"  - id: {spec['id']}",
            f"    label: {spec['label']}",
            f"    gene: {spec['gene']}",
            f"    source: {spec.get('source') or 'uniprot'}",
            f"    query: {scalar(query)}",
            f"    min_identity: {identity}",
            f"    max_sequences: {cap}",
            f"    length_tolerance: {LENGTH_TOLERANCE}",
            "",
        ]
    if panel.get("decoys"):
        out.append("decoys:")
        for decoy in panel["decoys"]:
            raw = decoy.get("query_raw")
            joined = raw or " OR ".join(
                f'(protein_name:"{n}")' for n in decoy["query_or"])
            out += [
                f"  - id: {decoy['id']}",
                f"    label: {decoy['label']}",
                "    source: uniprot",
                f"    query: {scalar(f'({joined}) AND (taxonomy_id:2)')}",
                f"    min_identity: {DECOY_IDENTITY}",
                f"    max_sequences: {int(decoy.get('max_sequences') or 2500)}",
                f"    length_tolerance: {LENGTH_TOLERANCE}",
                *block("note", decoy["note"], 4, 66),
                "",
            ]
    out.append("aggregate: sum")
    out.append(f"aggregate_from: [{', '.join(panel['aggregate_from'])}]")
    out += block("aggregate_reason", panel["aggregate_reason"], 0)
    out += ["", "min_fragments_for_stability: 20", "min_alignment_aa: 25"]
    if panel.get("expected_copies_per_genome"):
        out.append(f"expected_copies_per_genome: {panel['expected_copies_per_genome']}")
    out.append("")
    out += block("citation", panel["citation"], 0)
    out += ["", "interpretation:"]
    for key in ("what_it_is", "made_from", "summary", "evidence_detail"):
        out += block(key, panel[key], 2, 70)
    out.append(f"  higher_means: {panel['higher_means']}")
    out.append(f"  evidence_strength: {panel.get('evidence_strength', 'limited')}")
    out.append("  implications:")
    for key in ("higher", "typical", "lower"):
        out += block(key, panel["implications"][key], 4, 66)
    out += block("citation", panel["citation"], 2, 70)
    return "\n".join(out) + "\n"


def main() -> int:
    source = yaml.safe_load(SOURCE.read_text(encoding="utf-8"))
    panels = source.get("panels") or []
    for panel in panels:
        path = PANELS / f"{panel['name']}.yaml"
        path.write_text(render(panel), encoding="utf-8")
        n_targets = len(panel["targets"])
        n_decoys = len(panel.get("decoys") or [])
        print(f"  {path.name}: {n_targets} target(s), {n_decoys} decoy set(s)")
    print(f"{len(panels)} panel(s) written from {SOURCE.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
