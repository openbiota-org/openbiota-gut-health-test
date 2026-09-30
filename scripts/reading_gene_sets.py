"""Write each reading's gene-entry keys, after the panel-curation step.

Run between the harvest and the ladder build so the ladders are built over
exactly the gene sets the scorer will use at report time. Building them
over a different set is how the urolithin reading came to sit at the 46th
percentile against a ladder for three genes while its value counted one.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.extension.engine import _defer_to_panel_curation  # noqa: E402,PLC2701
from openbiota.samples import results_file  # noqa: E402


def main() -> int:
    # Any recent run carries the reading requirements and the panel layout.
    source = results_file("SAMPLE2_A02")
    results = json.loads(source.read_text())
    stored = (results.get("extension") or {}).get("views", {}).get("reading_scores") or {}
    if not stored:
        print("no reading scores to take gene sets from")
        return 1

    requirements = _defer_to_panel_curation([
        {"reading_id": rid, "genes": [c["gene"] for c in (v.get("contributors") or [])]}
        for rid, v in stored.items()
    ], results)

    keymap: dict[str, list[str]] = {}
    for panel in results.get("panels") or []:
        for kind in ("genes", "decoys"):
            for entry in panel.get(kind) or []:
                gene = str(entry.get("gene") or "").strip()
                if gene:
                    keymap.setdefault(gene, []).append(
                        f"{panel['name']}:{entry['entry_id']}")

    out: dict[str, list[str]] = {}
    for requirement in requirements:
        keys = sorted({k for g in requirement["genes"] for k in keymap.get(str(g), [])})
        if keys:
            out[str(requirement["reading_id"])] = keys
    Path("/tmp/reading_genes.json").write_text(json.dumps(out, indent=1))
    print(f"{len(out)} reading gene sets written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
