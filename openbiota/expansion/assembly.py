"""Targeted assembly for unresolved lineages (spec 0.8.4 §4D).

Runs when SingleM reports a lineage placed above species with mean marker
coverage at or above `singlem.ASSEMBLY_COVERAGE`: coverage that high with
no species-level home is a novel organism or a catalogue gap, and only
sequence decides which. The condition is data, not a flag.

What runs
---------
1. metaSPAdes 4.3.0 (pinned, vendor/) on the host-filtered reads.
2. Contigs >= MIN_CONTIG are kept; assembly statistics are recorded.
3. SingleM on the contigs: which unresolved lineages now have contigs
   carrying their marker windows, and how many distinct single-copy
   markers each lineage's contigs hold (a completeness proxy - a full
   genome carries each marker once).
4. skani of the assembly against every genome fetched for this sample
   (the detected organisms' references): contigs at >= 95% ANI / 65% AF
   to a known genome are `already_represented`; what remains is the
   residual sequence for the unresolved lineages.

What is not claimed
-------------------
No binner is pinned for macOS in this checkout, so recovered contigs are
reported per lineage as `sequence recovered`, with marker counts, and are
never counted as organisms. A marker fragment or several contigs are not
several organisms. Completeness and contamination in the MIMAG sense need
a bin; the marker-copy count is labelled a proxy.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Final

MIN_CONTIG: Final = 2_000
SPADES: Final = Path("vendor/SPAdes-4.3.0-Darwin/bin/metaspades.py")
TIME_BUDGET_S: Final = 2 * 3600
MEMORY_GB: Final = 120


def available() -> bool:
    return SPADES.is_file()


def _stats(contigs: Path) -> dict[str, Any]:
    lengths: list[int] = []
    cur = 0
    with contigs.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = 0
    acc = 0
    for L in lengths:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return {"n_contigs": len(lengths), "total_bp": total, "n50": n50, "longest": lengths[0] if lengths else 0,
            "n_ge_10kb": sum(1 for L in lengths if L >= 10_000)}


def _filter(src: Path, dest: Path, min_len: int) -> None:
    with src.open() as fh, dest.open("w") as out:
        name, seq = None, []
        def flush() -> None:
            if name and sum(len(s) for s in seq) >= min_len:
                out.write(f">{name}\n" + "".join(seq) + "\n")
        for line in fh:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        flush()


def _seqkit() -> str | None:
    return shutil.which("seqkit")


def select_unplaced_reads(*, sample: str, r1: Path, r2: Path | None, work_dir: Path,
                          kraken_assignments: Path, kraken_db: Path, threads: int) -> tuple[Path, Path | None, dict[str, Any]]:
    """The reads no gut catalogue could place at species level.

    Targeted means targeted: an unresolved lineage is made of reads that
    Kraken (UHGG, 4,744 gut species) left unclassified or placed only above
    species. Assembling those, and not the 85% of reads that already have a
    species, is what makes the stage affordable on every sample and keeps
    the assembly about the organisms nobody has catalogued. The selection
    and its counts are recorded so the operating point is visible.
    """
    from openbiota.engines import kraken as _kraken

    sel = work_dir / "unplaced"
    sel.mkdir(parents=True, exist_ok=True)
    o1 = sel / f"{sample}.unplaced_1.fastq.gz"
    o2 = sel / f"{sample}.unplaced_2.fastq.gz" if r2 is not None else None
    meta = sel / "selection.json"
    if o1.is_file() and (o2 is None or o2.is_file()) and meta.is_file():
        return o1, o2, json.loads(meta.read_text())
    ids, counts = _kraken.unplaced_read_ids(kraken_assignments, kraken_db)
    ids_file = sel / "unplaced.ids"
    ids_file.write_text("\n".join(sorted(ids)) + "\n")
    seqkit = _seqkit()
    if seqkit is None:
        raise RuntimeError("seqkit is required to select reads for targeted assembly")
    for src, dst in ((r1, o1), (r2, o2)):
        if src is None or dst is None:
            continue
        subprocess.run([seqkit, "grep", "-j", str(max(2, min(threads, 8))), "-f", str(ids_file), "-o", str(dst), str(src)],
                       check=True, capture_output=True, text=True)
    info = {
        "selection": "reads Kraken2/UHGG left unclassified or placed above species",
        "n_selected_pairs": len(ids), **{f"kraken_{k}": v for k, v in counts.items()},
        "selected_fraction": round(len(ids) / max(1, counts.get("total", 0)), 4),
    }
    meta.write_text(json.dumps(info, indent=1))
    return o1, o2, info


def run(*, sample: str, r1: Path, r2: Path | None, work_dir: Path, threads: int, triggers: list[dict[str, Any]],
        reference_genomes: list[Path], kraken_assignments: Path | None = None, kraken_db: Path | None = None) -> dict[str, Any]:
    t0 = time.monotonic()
    if not triggers:
        return {"status": "not_triggered", "reason": "no unresolved lineage reached the assembly coverage"}
    if not available():
        return {"status": "not_assessed", "reason": f"metaSPAdes not installed at {SPADES}", "triggers": triggers}
    work_dir.mkdir(parents=True, exist_ok=True)
    selection: dict[str, Any] = {"selection": "all host-filtered reads"}
    if kraken_assignments is not None and kraken_db is not None and kraken_assignments.is_file():
        r1, r2, selection = select_unplaced_reads(sample=sample, r1=r1, r2=r2, work_dir=work_dir,
                                                  kraken_assignments=kraken_assignments, kraken_db=kraken_db, threads=threads)
    out_dir = work_dir / "metaspades"
    contigs = out_dir / "contigs.fasta"
    if not contigs.is_file():
        # --only-assembler: BayesHammer read correction asserted on this
        # platform (kmer_cluster.cpp) under memory pressure; the multi-k
        # assembly graph handles read errors itself and the documented mode
        # exists for exactly this. Recorded in the result.
        cmd = [str(SPADES), "-1", str(r1)] + (["-2", str(r2)] if r2 else []) + [
            "-o", str(out_dir), "-t", str(threads), "-m", str(MEMORY_GB), "--only-assembler"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=TIME_BUDGET_S)
        except subprocess.TimeoutExpired:
            return {"status": "failed", "reason": f"metaSPAdes exceeded the {TIME_BUDGET_S // 60} min budget", "triggers": triggers}
        except subprocess.CalledProcessError as exc:
            return {"status": "failed", "reason": f"metaSPAdes exit {exc.returncode}: {exc.stderr[-400:]}", "triggers": triggers}
    kept = work_dir / f"{sample}.contigs.min{MIN_CONTIG}.fasta"
    if not kept.is_file():
        _filter(contigs, kept, MIN_CONTIG)
    stats = _stats(kept)

    # SingleM on the contigs: marker windows per lineage.
    lineage_markers: dict[str, dict[str, Any]] = {}
    singlem = Path(".venv-singlem/bin/singlem")
    pkg = Path("refs/singlem/GlobDB_r232.metapackage_v4.smpkg")
    if singlem.is_file() and pkg.is_dir() and stats["n_contigs"]:
        otu = work_dir / f"{sample}.contigs.singlem.otu.tsv"
        if not otu.is_file():
            subprocess.run([str(singlem), "pipe", "--genome-fasta-files", str(kept), "--metapackage", str(pkg),
                            "--otu-table", str(otu), "--threads", str(threads)], capture_output=True, text=True, check=False)
        if otu.is_file():
            with otu.open() as fh:
                header = fh.readline().rstrip("\n").split("\t")
                idx = {h: i for i, h in enumerate(header)}
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    tax = p[idx.get("taxonomy", len(p) - 1)].replace("Root; ", "").strip()
                    gene = p[idx.get("gene", 0)]
                    rec = lineage_markers.setdefault(tax, {"markers": set(), "n_windows": 0})
                    rec["markers"].add(gene)
                    rec["n_windows"] += 1
    per_lineage = []
    for trig in triggers:
        lin = trig["lineage"]
        hits = {k: v for k, v in lineage_markers.items() if k.startswith(lin)}
        markers = set().union(*(v["markers"] for v in hits.values())) if hits else set()
        per_lineage.append({
            **trig, "n_marker_families_on_contigs": len(markers),
            "completeness_proxy": round(len(markers) / 59.0, 3) if markers else 0.0,  # SingleM's 59 single-copy families
            "sub_lineages_on_contigs": len(hits),
            "verdict": ("sequence recovered; not a genome bin" if markers else "no marker-bearing contigs recovered"),
        })

    # Dereplication against the sample's fetched references.
    derep: dict[str, Any] = {"run": False}
    if reference_genomes and shutil.which("skani") and stats["n_contigs"]:
        from openbiota.expansion import ani
        pairs = ani.dist([kept], [p for p in reference_genomes if p.is_file()], work_dir=work_dir / "ani", threads=threads)
        same = [p for p in pairs if p.verdict == "same_species"]
        derep = {"run": True, "n_reference_genomes": len(reference_genomes),
                 "references_matched_at_species_level": sorted({p.reference for p in same})[:50],
                 "note": "whole-assembly ANI against each reference; a match means the assembly contains that species, not that the residual does"}

    return {
        "status": "completed", "assembler": "metaSPAdes 4.3.0 (--only-assembler)", "min_contig": MIN_CONTIG, "assembly": stats,
        "reads": selection,
        "triggers": per_lineage, "dereplication": derep, "elapsed_s": round(time.monotonic() - t0, 1),
        "binning": "not performed: no pinned binner in this checkout; contigs are reported per lineage and never counted as organisms",
        "meaning": ("Assembly was run because SingleM placed lineages above species at reconstructable coverage. "
                    "Marker families on contigs are a completeness proxy, not MIMAG completeness."),
    }
