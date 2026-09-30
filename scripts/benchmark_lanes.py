#!/usr/bin/env python3
"""Detection benchmark: simulated communities, every lane, the merged inventory.

Spec 0.8.4 §7 asks for engineering-target gates on simulated communities
(species precision >= 99%, recall >= 95%/98% at 8M/30M fragments, <= 5
false supported species per community, incremental gain over the installed
union with a confidence interval above zero) and says plainly that these
are targets, not claimed performance. This harness measures; it does not
tune after looking.

Design
------
* Truth genomes come from the HRGM2 species representatives already on disk
  (4,824 gut species, reconciled to GTDB R232). Two tiers are drawn:
    exact_reference   the HRGM2 representative is itself a GTDB R232 genome
    held_out_strain   the HRGM2 representative is a member of a known R232
                      species but not its representative (real divergence)
  Truth-to-canonical mapping (HRGM2 cluster -> R232 species) is frozen from
  the reconciliation table *before* any reads are simulated.
* Communities: `--species` organisms, half `even` (equal mass) and half
  `long_tailed` (log-normal mass, floor 0.02% of microbial DNA), each with
  a distinct seed. Reads: InSilicoSeq 2.0.1, HiSeq error model, paired,
  `--pairs` fragments.
* Every lane runs on each community exactly as the pipeline runs it, then
  the inventory is built two ways: baseline lanes only (MP3 scoring absent
  here; Jun23 + GTDB232) and all lanes. Scoring is at species level on the
  canonical R232 name; a parent/complex call is not a species true positive.
* Output: results/benchmarks/<name>/summary.json and summary.md with
  per-community precision, recall, false supported species, per-lane
  ablation, and the paired recall gain (all lanes - baseline) with a
  bootstrap CI over communities. Communities are the independent unit;
  seeds are technical replicates.

Compute is the constraint, not the design: one 8M-pair community costs
about an hour across the lanes. Run as many as the night allows and say
how many were run.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import math
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from openbiota import inventory  # noqa: E402
from openbiota.engines import (  # noqa: E402
    kraken,
    metaphlan_jan26,
    motus,
    singlem,
    sylph,
    sylph_globdb,
)
from openbiota.engines import metaphlan4 as mpa4  # noqa: E402
from openbiota.expansion import confirm, gtdb  # noqa: E402

RECON = REPO / "refs" / "expanded" / "reconciliation" / "hrgm2.tsv.gz"
HRGM2 = REPO / "refs" / "supplements" / "hrgm2" / "HRGMv2_Rep_Genome"
OUT = REPO / "results" / "benchmarks"


# --------------------------------------------------------------------------- #
# truth
# --------------------------------------------------------------------------- #

def truth_pool() -> tuple[list[dict], list[dict]]:
    """(exact_reference, held_out_strain) candidate genomes with frozen canonical species."""
    reps = {c.representative: sp for sp, c in gtdb.species_clusters("r232").items()}
    exact: list[dict] = []
    held: list[dict] = []
    with gzip.open(RECON, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["category"] != "already_represented" or not row["canonical_species"] or not row["local_path"]:
                continue
            path = Path(row["local_path"])
            if not path.is_file():
                continue
            rec = {"id": row["source_id"], "canonical": row["canonical_species"], "path": str(path),
                   "method": row["method"]}
            # An HRGM2 representative that carries a GTDB accession and is that species' representative
            # is exact; every other member is a held-out strain of a known species.
            (exact if rec["id"] in reps else held).append(rec)
    return exact, held


def design(species: int, kind: str, seed: int, exact: list[dict], held: list[dict]) -> list[dict]:
    """Half the species by their GTDB R232 representative genome (exact tier),
    half by the HRGM2 MAG of the same species (held-out strain tier)."""
    from openbiota.expansion import genomes as _genomes

    rng = random.Random(seed)
    pool = list(held)
    rng.shuffle(pool)
    chosen = pool[:species]
    reps = gtdb.representative_accessions("r232")
    exact_ids: set[str] = set()
    for i, g in enumerate(chosen):
        if i % 2 == 0:
            acc = reps.get(g["canonical"])
            if not acc:
                continue
            try:
                path = _genomes.fetch(acc)
            except Exception:  # noqa: BLE001 - fall back to the MAG for this species
                continue
            g["path"], g["id"] = str(path), acc
            exact_ids.add(acc)
    exact = [g for g in chosen if g["id"] in exact_ids]
    if kind == "even":
        masses = [1.0] * len(chosen)
    else:
        masses = [math.exp(rng.gauss(0, 1.2)) for _ in chosen]
    total = sum(masses)
    floor = 0.0002 * total  # 0.02% of microbial DNA mass
    masses = [max(m, floor) for m in masses]
    total = sum(masses)
    out = []
    for g, m in zip(chosen, masses, strict=True):
        out.append({**g, "mass_fraction": m / total, "tier": "exact_reference" if g in exact else "held_out_strain"})
    return out


# --------------------------------------------------------------------------- #
# simulation
# --------------------------------------------------------------------------- #

def _concat(genomes: list[dict], dest: Path) -> dict[str, str]:
    """One FASTA with contig names retagged to the genome id; returns contig -> id."""
    contig_of: dict[str, str] = {}
    with dest.open("w") as out:
        for g in genomes:
            opener = gzip.open if g["path"].endswith(".gz") else open
            with opener(g["path"], "rt") as fh:  # type: ignore[arg-type]
                for line in fh:
                    if line.startswith(">"):
                        name = f"{g['id']}__{line[1:].split()[0]}"
                        contig_of[name] = g["id"]
                        out.write(f">{name}\n")
                    else:
                        out.write(line)
    return contig_of


def simulate(community: list[dict], pairs: int, seed: int, work: Path, cpus: int) -> tuple[Path, Path]:
    work.mkdir(parents=True, exist_ok=True)
    r1, r2 = work / "sim_R1.fastq.gz", work / "sim_R2.fastq.gz"
    if r1.is_file() and r2.is_file():
        return r1, r2
    fasta = work / "genomes.fna"
    _concat(community, fasta)
    # iss abundance file is per *genome* (all contigs of a genome share its mass)... it wants one
    # abundance per FASTA record, so distribute each genome's mass over its contigs by length.
    lengths: dict[str, int] = {}
    with fasta.open() as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
    by_genome: dict[str, int] = {}
    for contig, n in lengths.items():
        by_genome[contig.split("__", 1)[0]] = by_genome.get(contig.split("__", 1)[0], 0) + n
    mass = {g["id"]: g["mass_fraction"] for g in community}
    ab = work / "abundance.txt"
    with ab.open("w") as out:
        for contig, n in lengths.items():
            gid = contig.split("__", 1)[0]
            out.write(f"{contig}\t{mass[gid] * n / by_genome[gid]:.10f}\n")
    iss = REPO / ".venv" / "bin" / "iss"
    cmd = [str(iss), "generate", "--genomes", str(fasta), "--abundance_file", str(ab), "--model", "hiseq",
           "--n_reads", str(pairs * 2), "--seed", str(seed), "--cpus", str(cpus), "--compress",
           "--output", str(work / "sim")]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    fasta.unlink(missing_ok=True)
    return r1, r2


# --------------------------------------------------------------------------- #
# negative controls: organism-free host and defined-food backgrounds
# --------------------------------------------------------------------------- #

NEGATIVE_BACKGROUNDS = (
    # (id, source FASTA, records to take, mass fraction)
    ("human_GRCh38", REPO / "refs/host/grch38.fna.gz", ("chr21", "chr22"), 0.70),
    ("bos_taurus", REPO / "refs/decoys/bos_taurus/genome.fna", None, 0.15),
    ("gallus_gallus", REPO / "refs/decoys/gallus_gallus/genome.fna", None, 0.15),
)


def _subset_fasta(src: Path, dest: Path, names: tuple[str, ...] | None, max_bases: int = 120_000_000) -> None:
    """Named records (or the first ones up to max_bases) of a large genome."""
    opener = gzip.open if src.suffix == ".gz" else open
    taken = 0
    keep = False
    with opener(src, "rt") as fh, dest.open("w") as out:  # type: ignore[arg-type]
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = (name in names) if names is not None else (taken < max_bases)
                if keep:
                    out.write(line)
            elif keep:
                out.write(line)
                taken += len(line) - 1
                if names is None and taken >= max_bases:
                    keep = False
    if taken == 0:
        raise ValueError(f"no sequence taken from {src}")


def negative_community(seed: int, work: Path) -> list[dict]:  # noqa: ARG001 - the background is fixed; seed is the community index in its name
    """A community with no microbial truth: human, cattle and chicken DNA.

    What a clean host/food background looks like to the lanes. The gate
    (spec §7) is zero false *supported* microbial species; provisional
    calls are reported as well because they are what a reader would see."""
    work.mkdir(parents=True, exist_ok=True)
    out = []
    for gid, src, names, mass in NEGATIVE_BACKGROUNDS:
        sub = work / f"{gid}.fna"
        if not sub.is_file():
            _subset_fasta(src, sub, names)
        out.append({"id": gid, "canonical": gid, "path": str(sub), "mass_fraction": mass, "tier": "background"})
    return out


def score_negative(results: dict, confirmation: dict | None) -> dict:
    inv = inventory.from_results(results)
    blob = inv.to_json()
    if confirmation:
        blob = confirm.apply_verdicts(blob, confirmation)
        inv = inventory.from_json(blob) or inv
    supported = [o for o in inv.organisms if o.status == "supported"]
    provisional = [o for o in inv.organisms if o.status == "provisional"]
    return {
        "n_supported": len(supported), "n_provisional": len(provisional),
        "supported": [{"organism": o.gtdb or o.species, "lanes": list(o.lanes), "percent": o.best_percent} for o in supported],
        "provisional": [{"organism": o.gtdb or o.species, "lanes": list(o.lanes), "percent": o.best_percent} for o in provisional][:50],
        "gate_zero_false_supported": len(supported) == 0,
    }


# --------------------------------------------------------------------------- #
# lanes and scoring
# --------------------------------------------------------------------------- #

def run_lanes(sample: str, r1: Path, r2: Path, work: Path, threads: int) -> dict:
    results: dict = {"sample": sample}
    t = {}
    from openbiota.engines.metaphlan4 import locate_tools4
    t0 = time.monotonic()
    try:
        ext = mpa4.run_metaphlan4(locate_tools4(), r1=r1, r2=r2, db_dir=REPO / "refs/metaphlan4_db", work_dir=work / "taxonomy",
                                  reporter=_Quiet(), threads=threads)
        results["extended_catalogue"] = ext.to_json()
    except Exception as exc:  # noqa: BLE001
        results["extended_catalogue"] = {"error": str(exc)[:200]}
    t["jun23"] = time.monotonic() - t0
    t0 = time.monotonic()
    gp = sylph.run_sylph(sample=sample, r1=r1, r2=r2, db_dir=REPO / "refs/sylph", work_dir=work / "genome", threads=threads)
    results["genome_profile"] = gp.to_json() if gp else {}
    t["gtdb_sylph"] = time.monotonic() - t0
    lanes: dict = {}
    for lane_id, fn in (
        ("metaphlan_jan26", lambda: metaphlan_jan26.run_metaphlan_jan26(sample=sample, r1=r1, r2=r2, db_dir=REPO / "refs/metaphlan4_db",
                                                                         work_dir=work / "taxonomy", strain_dir=work / "strain", threads=threads)),
        ("sylph_globdb", lambda: sylph_globdb.run_sylph_globdb(sample=sample, r1=r1, r2=r2, refs_dir=REPO / "refs", work_dir=work / "globdb", threads=threads)),
        ("motus4", lambda: motus.run_motus(sample=sample, r1=r1, r2=r2, refs_dir=REPO / "refs", work_dir=work / "motus", threads=threads)),
        ("kraken_uhgg", lambda: kraken.run_kraken(sample=sample, r1=r1, r2=r2, refs_dir=REPO / "refs", work_dir=work / "kraken", threads=threads)),
        ("singlem_globdb", lambda: singlem.run_singlem(sample=sample, r1=r1, r2=r2, refs_dir=REPO / "refs", work_dir=work / "singlem", threads=threads)),
    ):
        t0 = time.monotonic()
        try:
            r = fn()
            if r is not None:
                lanes[lane_id] = r.to_json()
        except Exception as exc:  # noqa: BLE001
            lanes[lane_id] = {"error": str(exc)[:300], "observations": []}
        t[lane_id] = time.monotonic() - t0
    results["detection"] = {"lanes": lanes}
    results["timings_s"] = {k: round(v, 1) for k, v in t.items()}
    return results


class _Quiet:
    def info(self, *_: object) -> None: ...
    def ok(self, *_: object) -> None: ...
    def record(self, *_: object) -> None: ...
    def warn(self, *_: object) -> None: ...


def _species_of(o: inventory.Organism) -> str:
    """The canonical R232 species an organism was called as.

    '' for complexes, genus-level placements and placeholder displays such
    as `Coprococcus sp. (SGB5119)`: those are not exact-species calls and
    are scored separately as higher-rank findings.
    """
    if o.count_category in ("unresolved_complex", "higher_rank"):
        return ""
    name = (o.gtdb or "").strip()
    if not name or " sp. (" in name or name.endswith(")"):
        name = o.species.replace("_", " ").strip() if not o.unnamed else ""
    return name


def _genus_of_call(o: inventory.Organism) -> str:
    g = (o.gtdb_genus or o.genus or "").strip()
    return g.split("_")[0] if g else ""


def score(results: dict, truth: list[dict], *, lanes_subset: set[str] | None, confirmation: dict | None = None) -> dict:
    """Species precision/recall for one inventory built from a subset of lanes."""
    r = dict(results)
    if lanes_subset is not None:
        r["detection"] = {"lanes": {k: v for k, v in results["detection"]["lanes"].items() if k in lanes_subset}}
        if "extended" not in lanes_subset:
            r["extended_catalogue"] = {}
        if "genome" not in lanes_subset:
            r["genome_profile"] = {}
    inv = inventory.from_results(r)
    if confirmation is not None:
        blob = confirm.apply_verdicts(inv.to_json(), confirmation)
        inv = inventory.from_json(blob) or inv
    truth_species = {g["canonical"] for g in truth}
    truth_genera = {g["canonical"].split(" ")[0].split("_")[0] for g in truth}
    called = {_species_of(o) for o in inv.organisms if o.status == "supported" and _species_of(o)}
    called_any = {_species_of(o) for o in inv.organisms if _species_of(o)}
    higher = [o for o in inv.organisms if o.status == "supported" and not _species_of(o)]
    # A GlobDB or UHGG unit is a different clustering of the same genomes. A
    # called unit whose genome is the same species (ANI/AF) as a truth genome
    # of the same genus is that organism under another catalogue's name, not
    # a false positive; it is credited to the truth species it matches.
    equivalents = _unit_equivalents(called - truth_species, truth, inv)
    called_eq = {equivalents.get(c, c) for c in called}
    tp = len(called_eq & truth_species)
    fp = len(called_eq - truth_species)
    fn = len(truth_species - called_eq)
    # A false positive whose GTDB sibling stem matches a truth species
    # ("Collinsella aerofaciens" called, "Collinsella aerofaciens_F" in
    # truth) is a call at the right species complex under the wrong split
    # label. It stays a false positive under the spec's exact-species rule
    # and is also counted on its own, so the two failure modes - sibling
    # split versus an unrelated call - are visible separately.
    truth_stems = {inventory._sibling_base(t) for t in truth_species}
    sibling_fp = sorted(f for f in (called_eq - truth_species) if inventory._sibling_base(f) in truth_stems)
    return {
        "n_truth": len(truth_species), "n_called_supported": len(called), "n_called_any": len(called_any),
        "tp": tp, "fp": fp, "fn": fn,
        "false_positive_species": sorted(called_eq - truth_species),
        "false_positive_sibling_split": sibling_fp,
        "fp_sibling_split": len(sibling_fp), "fp_other": fp - len(sibling_fp),
        "precision_complex_level": (tp + len(sibling_fp)) / len(called) if called else None,
        "missed_species": sorted(truth_species - called_eq),
        "catalogue_unit_equivalents": equivalents,
        "higher_rank_calls": len(higher),
        "higher_rank_genus_correct": sum(1 for o in higher if _genus_of_call(o) in truth_genera),
        "precision": tp / len(called) if called else None, "recall": tp / len(truth_species) if truth_species else None,
        "false_supported_species": fp,
        "recall_any_status": len(called_any & truth_species) / len(truth_species) if truth_species else None,
        "by_tier": {
            tier: (sum(1 for g in truth if g["tier"] == tier and g["canonical"] in called) /
                   max(1, sum(1 for g in truth if g["tier"] == tier)))
            for tier in ("exact_reference", "held_out_strain")
        },
    }


def _unit_equivalents(candidates: set[str], truth: list[dict], inv: inventory.Inventory) -> dict[str, str]:
    """Called catalogue units -> the truth species their genome is the same species as.

    Only non-GTDB units (GlobDB/UHGG ids) are considered, only against truth
    genomes of the same genus, by skani ANI >= 95 and AF(shorter) >= 65.
    Units whose genome is not on disk are left as they are."""
    from openbiota.expansion import ani, genomes
    if not candidates or not ani.available():
        return {}
    by_display = {(o.gtdb or o.species.replace("_", " ")): o for o in inv.organisms}
    out: dict[str, str] = {}
    queries: list[tuple[str, Path]] = []
    for name in candidates:
        o = by_display.get(name)
        if o is None or not o.unnamed:
            continue
        gid = next((v.split(":", 1)[-1] for k, v in (o.native_ids or {}).items() if k in ("globdb", "kraken")), "")
        gid = gid.split(";")[0]
        if not gid or gid.startswith(("GCA_", "GCF_", "uhgg-taxid")):
            continue
        try:
            queries.append((name, genomes.fetch(gid)))
        except Exception:  # noqa: BLE001 - pending archive or unfetchable: stays as called
            continue
    if not queries:
        return {}
    refs = {g["canonical"]: Path(g["path"]) for g in truth if Path(g["path"]).is_file()}
    pairs = ani.dist([q for _, q in queries], list(refs.values()), work_dir=Path("/tmp/ob_bench_ani"), threads=8)
    best = ani.best_by_query(pairs)
    ref_by_name = {v.name: k for k, v in refs.items()}
    for name, q in queries:
        b = best.get(q.name)
        if b is not None and b.verdict == "same_species":
            truth_sp = ref_by_name.get(b.reference)
            if truth_sp and truth_sp.split(" ")[0].split("_")[0] == name.split(" ")[0].split("_")[0]:
                out[name] = truth_sp
    return out


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def bootstrap_mean_ci(xs: list[float], n: int = 2000, seed: int = 7) -> tuple[float, float]:
    if not xs:
        return (0.0, 0.0)
    rng = random.Random(seed)
    means = sorted(statistics.fmean(rng.choices(xs, k=len(xs))) for _ in range(n))
    return (means[int(0.025 * n)], means[int(0.975 * n) - 1])


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def write_summary(args: Any, out: Path, per_community: list[dict], negatives: list[dict]) -> int:
    """summary.json + summary.md for a benchmark directory."""
    if negatives:
        (out / "negatives.json").write_text(json.dumps(negatives, indent=1))
    if not per_community:
        summary = {"name": args.name, "n_communities": 0, "negatives": negatives,
                   "negative_gate_zero_false_supported": all(n["gate_zero_false_supported"] for n in negatives) if negatives else None,
                   "real_blanks": "none available; disclosed"}
        (out / "summary.json").write_text(json.dumps(summary, indent=1))
        print(json.dumps(summary, indent=1))
        return 0
    gains = [c["recall_gain"] for c in per_community]
    tp = sum(c["all_lanes"]["tp"] for c in per_community)
    called = sum(c["all_lanes"]["n_called_supported"] for c in per_community)
    summary = {
        "name": args.name, "run_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_communities": len(per_community), "species_per_community": args.species, "pairs": args.pairs,
        "mean_recall_baseline": statistics.fmean(c["baseline"]["recall"] for c in per_community),
        "mean_recall_all_lanes": statistics.fmean(c["all_lanes"]["recall"] for c in per_community),
        "pooled_precision_all_lanes": tp / called if called else None,
        "pooled_precision_ci95": wilson(tp, called),
        "pooled_precision_complex_level": (tp + sum(c["all_lanes"]["fp_sibling_split"] for c in per_community)) / called if called else None,
        "mean_fp_sibling_split_per_community": statistics.fmean(c["all_lanes"]["fp_sibling_split"] for c in per_community),
        "mean_fp_other_per_community": statistics.fmean(c["all_lanes"]["fp_other"] for c in per_community),
        "mean_false_supported_per_community": statistics.fmean(c["all_lanes"]["false_supported_species"] for c in per_community),
        "paired_recall_gain_mean": statistics.fmean(gains), "paired_recall_gain_ci95_bootstrap": bootstrap_mean_ci(gains),
        "gates": {
            "precision_ge_0.99": (tp / called if called else 0) >= 0.99,
            "recall_ge_0.95_at_8M": statistics.fmean(c["all_lanes"]["recall"] for c in per_community) >= 0.95,
            "false_supported_le_5": all(c["all_lanes"]["false_supported_species"] <= 5 for c in per_community),
            "gain_ci_above_zero": bootstrap_mean_ci(gains)[0] > 0 if len(gains) >= 2 else None,
        },
        "communities": per_community,
        "negatives": negatives,
        "negative_gate_zero_false_supported": all(n["gate_zero_false_supported"] for n in negatives) if negatives else None,
        "real_blanks": "none available; disclosed (spec §7 rare/open-set gate)",
        "caveats": [
            f"{len(per_community)} communities is far below the 20 the spec asks for; the CI is correspondingly wide",
            "the MetaPhlAn 3 scoring lane is absent from the baseline here (its cohort context is irrelevant to detection)",
            "competitive confirmation ran inside the benchmark; candidates whose reference genome is not on disk stay provisional (pending)",
            "truth genomes are HRGM2 representatives; a species whose HRGM2 genome is a divergent member of its R232 species is the held-out tier",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    md = [f"# Detection benchmark {args.name}", "",
          f"{len(per_community)} communities x {args.species} species x {args.pairs:,} pairs", "",
          "| community | baseline recall | all-lanes recall | precision | false supported | gain |", "|---|---|---|---|---|---|"]
    for c in per_community:
        md.append(f"| {c['community']} | {c['baseline']['recall']:.3f} | {c['all_lanes']['recall']:.3f} | "
                  f"{(c['all_lanes']['precision'] or 0):.3f} | {c['all_lanes']['false_supported_species']} | {c['recall_gain']:+.3f} |")
    md += ["", f"pooled precision {summary['pooled_precision_all_lanes']:.4f} (95% CI {summary['pooled_precision_ci95'][0]:.3f}-{summary['pooled_precision_ci95'][1]:.3f}); "
           f"at species-complex level (sibling-split calls credited) {summary['pooled_precision_complex_level']:.4f}; "
           f"false supported per community: {summary['mean_fp_sibling_split_per_community']:.1f} sibling-split, {summary['mean_fp_other_per_community']:.1f} other",
           f"paired recall gain {summary['paired_recall_gain_mean']:+.3f} (bootstrap 95% CI {summary['paired_recall_gain_ci95_bootstrap'][0]:+.3f} to {summary['paired_recall_gain_ci95_bootstrap'][1]:+.3f})",
           "", "gates: " + json.dumps(summary["gates"])]
    if negatives:
        md += ["", f"organism-free backgrounds: {len(negatives)}; false supported microbial species: "
               + ", ".join(str(n["n_supported"]) for n in negatives)
               + "; provisional: " + ", ".join(str(n["n_provisional"]) for n in negatives)
               + f" (gate {'MET' if summary['negative_gate_zero_false_supported'] else 'FAILED'}); real blanks: none available"]
    md += ["", "caveats:", *[f"- {c}" for c in summary["caveats"]]]
    (out / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0


def rescore(out: Path) -> int:
    """Re-score every community in a finished benchmark directory (no reruns)."""
    per_community: list[dict] = []
    for cdir in sorted(p for p in out.iterdir() if p.is_dir() and p.name.startswith("c")):
        results = json.loads((cdir / "results.json").read_text())
        truth = json.loads((cdir / "truth.json").read_text())
        old = json.loads((cdir / "score.json").read_text()) if (cdir / "score.json").is_file() else {}
        confs = list((cdir / "confirmation").glob("*.confirmation.json"))
        conf = json.loads(confs[0].read_text()) if confs else None
        baseline = score(results, truth, lanes_subset={"extended", "genome"})
        all_lanes = score(results, truth, lanes_subset=None, confirmation=conf)
        rec = {**old, "baseline": baseline, "all_lanes": all_lanes,
               "all_lanes_before_confirmation": score(results, truth, lanes_subset=None),
               "recall_gain": (all_lanes["recall"] or 0) - (baseline["recall"] or 0), "rescored": True}
        rec.setdefault("community", cdir.name)
        rec.setdefault("kind", "even" if "even" in cdir.name else "long_tailed")
        per_community.append(rec)
        (cdir / "score.json").write_text(json.dumps(rec, indent=1))
        print(f"[{cdir.name}] recall {all_lanes['recall']:.3f} precision {all_lanes['precision']:.3f} "
              f"(complex-level {all_lanes['precision_complex_level']:.3f}) FP {all_lanes['fp']} = "
              f"{all_lanes['fp_sibling_split']} sibling-split + {all_lanes['fp_other']} other")
    negatives = json.loads((out / "negatives.json").read_text()) if (out / "negatives.json").is_file() else []
    old_summary = json.loads((out / "summary.json").read_text()) if (out / "summary.json").is_file() else {}
    class _A:  # the fields write_summary needs
        name = out.name
        species = old_summary.get("species_per_community", per_community[0]["n_species"] if per_community else 0)
        pairs = old_summary.get("pairs", per_community[0].get("pairs", 0) if per_community else 0)
    write_summary(_A(), out, per_community, negatives)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=dt.datetime.now().strftime("bench_%Y%m%d_%H%M"))
    ap.add_argument("--communities", type=int, default=2)
    ap.add_argument("--species", type=int, default=600)
    ap.add_argument("--pairs", type=int, default=8_000_000)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--keep-reads", action="store_true")
    ap.add_argument("--negative", type=int, default=0,
                    help="also run this many organism-free host/food background communities (spec §7 open-set gate)")
    ap.add_argument("--negative-pairs", type=int, default=2_000_000)
    ap.add_argument("--rescore", action="store_true",
                    help="recompute score.json and the summary for --name from its saved results, truth and confirmation")
    args = ap.parse_args()
    if args.rescore:
        return rescore(OUT / args.name)

    out = OUT / args.name
    out.mkdir(parents=True, exist_ok=True)
    negatives: list[dict] = []
    for j in range(args.negative):
        seed = args.seed * 1000 + 900 + j
        cname = f"n{j:02d}_background_s{seed}"
        cdir = out / cname
        cdir.mkdir(exist_ok=True)
        community = negative_community(seed, cdir / "backgrounds")
        (cdir / "truth.json").write_text(json.dumps([], indent=1))
        print(f"[{cname}] organism-free background; simulating {args.negative_pairs:,} pairs ...", flush=True)
        t0 = time.monotonic()
        r1, r2 = simulate(community, args.negative_pairs, seed, cdir, args.threads)
        sim_s = time.monotonic() - t0
        results = run_lanes(cname, r1, r2, cdir, args.threads)
        (cdir / "results.json").write_text(json.dumps(results, indent=1))
        conf = None
        try:
            inv_all = inventory.from_results(results)
            if confirm.available() and confirm.select_candidates(inv_all):
                conf = confirm.run_confirmation(sample=cname, inventory=inv_all, r1=r1, r2=r2,
                                                work_dir=cdir / "confirmation", threads=args.threads)
        except Exception as exc:  # noqa: BLE001
            print(f"[{cname}] confirmation failed: {exc}", flush=True)
        neg = score_negative(results, conf)
        rec = {"community": cname, "kind": "background", "seed": seed, "pairs": args.negative_pairs,
               "simulation_s": round(sim_s, 1), "timings_s": results["timings_s"], **neg}
        negatives.append(rec)
        (cdir / "score.json").write_text(json.dumps(rec, indent=1))
        print(f"[{cname}] supported microbial calls: {neg['n_supported']} (gate {'MET' if neg['gate_zero_false_supported'] else 'FAILED'}); "
              f"provisional {neg['n_provisional']}", flush=True)
        if not args.keep_reads:
            for p in (r1, r2):
                p.unlink(missing_ok=True)

    exact, held = truth_pool()
    print(f"truth pool: {len(exact)} exact-reference, {len(held)} held-out-strain genomes")
    if len(exact) + len(held) < args.species:
        print(f"pool too small for {args.species} species; using {len(exact) + len(held)}")
        args.species = len(exact) + len(held)

    per_community: list[dict] = []
    for i in range(args.communities):
        kind = "even" if i % 2 == 0 else "long_tailed"
        seed = args.seed * 1000 + i
        cname = f"c{i:02d}_{kind}_s{seed}"
        cdir = out / cname
        cdir.mkdir(exist_ok=True)
        truth = design(args.species, kind, seed, exact, held)
        (cdir / "truth.json").write_text(json.dumps(truth, indent=1))   # frozen before simulation
        print(f"[{cname}] {len(truth)} species; simulating {args.pairs:,} pairs ...", flush=True)
        t0 = time.monotonic()
        r1, r2 = simulate(truth, args.pairs, seed, cdir, args.threads)
        sim_s = time.monotonic() - t0
        print(f"[{cname}] simulated in {sim_s:.0f}s; running lanes ...", flush=True)
        results = run_lanes(cname, r1, r2, cdir, args.threads)
        (cdir / "results.json").write_text(json.dumps(results, indent=1))
        baseline = score(results, truth, lanes_subset={"extended", "genome"})
        # Competitive confirmation on the merged inventory, as the pipeline runs it.
        conf = None
        try:
            inv_all = inventory.from_results(results)
            if confirm.available() and confirm.select_candidates(inv_all):
                t0 = time.monotonic()
                conf = confirm.run_confirmation(sample=cname, inventory=inv_all, r1=r1, r2=r2,
                                                work_dir=cdir / "confirmation", threads=args.threads)
                conf["elapsed_s"] = round(time.monotonic() - t0, 1)
                print(f"[{cname}] confirmation: {conf.get('n_tested')} tested, "
                      f"{sum(1 for v in conf['verdicts'] if v['status'] == 'not_detected')} rejected, "
                      f"{sum(1 for v in conf['verdicts'] if v['status'] == 'pending')} pending a genome ({conf['elapsed_s']}s)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[{cname}] confirmation failed: {exc}", flush=True)
        all_lanes = score(results, truth, lanes_subset=None, confirmation=conf)
        all_lanes_unconfirmed = score(results, truth, lanes_subset=None)
        ablation = {}
        for lane in ("metaphlan_jan26", "sylph_globdb", "motus4", "kraken_uhgg", "singlem_globdb"):
            ablation[lane] = score(results, truth, lanes_subset={"extended", "genome", lane})
        rec = {"community": cname, "kind": kind, "seed": seed, "n_species": len(truth), "pairs": args.pairs,
               "simulation_s": round(sim_s, 1), "timings_s": results["timings_s"],
               "baseline": baseline, "all_lanes": all_lanes, "all_lanes_before_confirmation": all_lanes_unconfirmed,
               "confirmation": {k: v for k, v in (conf or {}).items() if k != "verdicts"}, "ablation_baseline_plus": ablation,
               "recall_gain": (all_lanes["recall"] or 0) - (baseline["recall"] or 0)}
        per_community.append(rec)
        (cdir / "score.json").write_text(json.dumps(rec, indent=1))
        print(f"[{cname}] baseline recall {baseline['recall']:.3f} precision {baseline['precision'] or 0:.3f} | "
              f"all lanes recall {all_lanes['recall']:.3f} precision {all_lanes['precision'] or 0:.3f} "
              f"FP {all_lanes['false_supported_species']}", flush=True)
        if not args.keep_reads:
            for p in (r1, r2):
                p.unlink(missing_ok=True)

    return write_summary(args, out, per_community, negatives)


if __name__ == "__main__":
    raise SystemExit(main())
