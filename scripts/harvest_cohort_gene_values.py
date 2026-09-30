"""Recover per-sample gene values from the cohort alignments already on disk.

The reference file stores a percentile ladder per panel and per gene, and
nothing for a combination of genes. That is why a reading over several
genes had no honest position: there was no distribution to place it in.

There does not need to be a new cohort run to build one. Every alignment
from last night is still here, 200 files over 100 samples, all against the
database in use. Re-tallying them costs no download and no DIAMOND - the
search step finds its own cache - and it yields exactly the per-sample
numbers the ladders were reduced from.

Writes /tmp/cohort_gene_values.json: {gene_key: [value per sample]} plus
the sample order, so any reading's gene set can be summed per sample and
turned into a ladder of its own.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.cli import (  # noqa: E402
    Reporter,
    _entries_for_db,
    build_reference_database,
    iter_lines,
    load_panel_set,
    output_fields,
    tally,
)

COHORT = REPO / "refs" / "cohort"
def current_db() -> str:
    """The fingerprint of the database in use, read rather than assumed.

    Hardcoding it meant the harvest found nothing the moment the database
    was rebuilt, which is exactly when it needs to run."""
    from openbiota.cli import _entries_for_db, build_reference_database, load_panel_set

    panels = load_panel_set(REPO / "panels")
    db = build_reference_database(
        panels, _entries_for_db(panels, None), refs_dir=REPO / "refs",
        reporter=Reporter(verbose=False), diamond="diamond", threads=8,
    )
    return str(db.fingerprint)


def usable_samples(fingerprint: str) -> list[str]:
    """Cohort samples whose alignments were made against the current database."""
    out = []
    for marker in sorted(COHORT.glob("results/*/alignments/hits_all_R1.tsv.done.json")):
        try:
            if json.loads(marker.read_text()).get("db_fingerprint") == fingerprint:
                out.append(marker.parents[1].name)
        except (OSError, ValueError):
            continue
    return out


def main() -> int:
    reporter = Reporter(verbose=False)
    panel_set = load_panel_set(REPO / "panels")
    database = build_reference_database(
        panel_set, _entries_for_db(panel_set, None),
        refs_dir=REPO / "refs", reporter=reporter, diamond="diamond", threads=8,
    )
    fields = output_fields(with_residues=bool(database.anchors))

    samples = usable_samples(current_db())
    print(f"{len(samples)} cohort samples aligned against the database in use")

    gene_values: dict[str, list[float]] = {}
    order: list[str] = []
    started = time.monotonic()
    for i, srr in enumerate(samples, start=1):
        align = COHORT / "results" / srr / "alignments"
        try:
            # The search step is already done and its output is on disk, so
            # the tally is all that is left: no FASTQ, no DIAMOND, no network.
            streams = [
                (mate, iter_lines(align / f"hits_all_{mate}.tsv"))
                for mate in ("R1", "R2")
                if (align / f"hits_all_{mate}.tsv").is_file()
            ]
            if not streams:
                continue
            result = tally(streams=streams, database=database,
                           panel_set=panel_set, fields=fields)
        except Exception as exc:  # noqa: BLE001 - one sample cannot stop the harvest
            print(f"  [{i}/{len(samples)}] {srr}: skipped ({type(exc).__name__}: {exc})"[:150])
            continue
        if result.normalizer.fragments < 300:
            print(f"  [{i}/{len(samples)}] {srr}: below the rpoB floor, skipped")
            continue
        order.append(srr)
        for panel in result.panels:
            for gene in panel.targets:
                if gene.copies_per_100_genomes is not None:
                    gene_values.setdefault(gene.key, []).append(
                        float(gene.copies_per_100_genomes))
        if i % 10 == 0 or i == len(samples):
            print(f"  [{i}/{len(samples)}] {len(order)} usable, "
                  f"{time.monotonic() - started:.0f}s")

    payload = {"samples": order, "n": len(order), "genes": gene_values}
    Path("/tmp/cohort_gene_values.json").write_text(json.dumps(payload))
    print(f"\n{len(order)} samples, {len(gene_values)} genes -> /tmp/cohort_gene_values.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
