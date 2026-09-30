"""Kraken2 + Bracken gut rescue lane (spec 0.8.4 §4B).

A read-classification method that searches differently from marker
profiling and from genome sketching: every k-mer of every read against a
gut-focused panel. The panel here is the publisher-built UHGG v2.0.2
Kraken2 database (4,744 gut species, 289k genomes), extended by the
project's rescue panel when one has been built (see
`openbiota.expansion.reconcile`).

What this lane may and may not claim
------------------------------------
* Kraken's confidence score is not a probability; `CONFIDENCE` and
  `MIN_HIT_GROUPS` are operating points, recorded with every result and
  calibrated on tuning communities (`scripts/benchmark_lanes.py`).
* Bracken redistributes reads to species; it is an abundance step, not an
  independent confirmation. Its numbers are read fractions.
* UHGG's taxonomy uses GTDB-r95-era names and its own integer taxids. The
  taxid is the native identifier; the species representative (an MGYG
  accession) is what the crosswalk resolves through GlobDB to a current
  cluster. Names are never joined to another release by string.
* `--report-minimizer-data` is kept: distinct minimizers per taxon are the
  support metric that separates a real genome from shared k-mers.
"""

from __future__ import annotations

import contextlib
import gzip
import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Final

from openbiota.engines.lanebase import LaneResult, Observation, cache_key, input_sha256, read_cached

LANE_ID: Final = "kraken_uhgg"
DB_RELEASE: Final = "UHGG v2.0.2 Kraken2 (publisher build)"
DB_DIRNAME: Final = "uhgg_v2.0.2"

#: The two Kraken panels this engine runs. `uhgg` is the publisher's gut
#: catalogue database; `rescue` is the combined gut rescue panel built by
#: scripts/build_rescue_panel_kraken.py from the reconciliation targets,
#: their near neighbours, the detected alternatives and background decoys
#: (spec 0.8.4 §4B). Same tool, same operating point, one method.
PANELS: Final = {
    "uhgg": {
        "lane_id": "kraken_uhgg", "dirname": "uhgg_v2.0.2", "release": DB_RELEASE,
        "taxonomy": "UHGG v2.0.2 (GTDB r95-era names; MGYG representatives)", "id_prefix": "uhgg-taxid",
        "lock": "refs/expanded/locks/uhgg_kraken2.lock.json", "rep_key": "uhgg_species_rep",
    },
    "rescue": {
        "lane_id": "kraken_rescue", "dirname": "rescue_panel",
        "release": "gut rescue panel (local build: reconciliation targets, near neighbours, detected alternatives, decoys)",
        "taxonomy": "custom; GlobDB / GTDB R232 lineages, taxid <-> source id in panel_manifest.tsv", "id_prefix": "panel-taxid",
        "lock": None, "rep_key": "panel_genome",
    },
}
#: Operating point. Recorded in every result; calibrated, not guessed, by
#: the lane benchmark. A custom gut database with many close species needs
#: a floor above Kraken's permissive default of 0.
CONFIDENCE: Final = 0.15
MIN_HIT_GROUPS: Final = 3
#: Bracken species-level threshold: reads at species level needed before
#: redistribution counts a species.
BRACKEN_THRESHOLD: Final = 10
#: Reporting floor for the inventory: below this a species stays in the
#: lane's raw record (fragment assignments are retained) but is not an
#: organism. On one 5.7M-pair sample Kraken put reads on 3,179 species;
#: 279 cleared this floor. Shared k-mers between close species produce
#: exactly the long tail this removes.
INVENTORY_MIN_READS: Final = 1_000
INVENTORY_MIN_MINIMIZERS: Final = 20_000
#: The lane's own "supported": still only a candidate for confirmation.
SUPPORTED_MIN_READS: Final = 5_000
SUPPORTED_MIN_MINIMIZERS: Final = 100_000
CACHE_VERSION: Final = 2

_REQUIRED = ("hash.k2d", "opts.k2d", "taxo.k2d", "seqid2taxid.map", "taxonomy/names.dmp", "taxonomy/nodes.dmp")


def _tools() -> tuple[str | None, str | None]:
    return shutil.which("kraken2"), shutil.which("bracken")


def kraken_version() -> str:
    exe, _ = _tools()
    if not exe:
        return ""
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, check=False).stdout
    return out.splitlines()[0].replace("Kraken version", "").strip() if out else ""


def database_ready(refs_dir: Path, panel: str = "uhgg") -> bool:
    spec = PANELS[panel]
    db = refs_dir / "kraken2" / spec["dirname"]
    if not all((db / f).is_file() for f in _REQUIRED):
        return False
    if panel == "rescue":
        # built locally; build.json is written only after kraken2-build finished
        return (db / "build.json").is_file()
    # The fetcher writes the lock only when hash.k2d verified in full.
    lock = Path(spec["lock"])
    if lock.is_file():
        rec = json.loads(lock.read_text()).get("files", {}).get(f"kraken2/{spec['dirname']}/hash.k2d")
        return bool(rec and rec.get("verified"))
    return True


def available(refs_dir: Path, panel: str = "uhgg") -> bool:
    k, b = _tools()
    return k is not None and b is not None and database_ready(refs_dir, panel)


def _read_length(r1: Path, n: int = 2000) -> int:
    opener = gzip.open if r1.suffix == ".gz" else open
    total = count = 0
    with opener(r1, "rt") as fh:  # type: ignore[arg-type]
        for i, line in enumerate(fh):
            if i % 4 == 1:
                total += len(line.strip())
                count += 1
                if count >= n:
                    break
    return int(round(total / count)) if count else 150


def _taxonomy(db: Path) -> tuple[dict[int, str], dict[int, tuple[int, str]]]:
    names: dict[int, str] = {}
    with (db / "taxonomy" / "names.dmp").open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[int(parts[0])] = parts[1]
    nodes: dict[int, tuple[int, str]] = {}
    with (db / "taxonomy" / "nodes.dmp").open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 3:
                nodes[int(parts[0])] = (int(parts[1]), parts[2])
    return names, nodes


def _lineage(taxid: int, names: dict[int, str], nodes: dict[int, tuple[int, str]]) -> str:
    prefix = {"domain": "d", "superkingdom": "d", "phylum": "p", "class": "c", "order": "o",
              "family": "f", "genus": "g", "species": "s"}
    chain: list[tuple[str, str]] = []
    t = taxid
    seen = 0
    while t in nodes and t != 1 and seen < 40:
        parent, rank = nodes[t]
        if rank in prefix:
            chain.append((prefix[rank], names.get(t, str(t))))
        if parent == t:
            break
        t = parent
        seen += 1
    chain.reverse()
    return ";".join(f"{p}__{n}" for p, n in chain)


def _species_reps(db: Path) -> dict[int, str]:
    """taxid -> the MGYG accession whose contigs the database holds for it."""
    reps: dict[int, str] = {}
    with (db / "seqid2taxid.map").open() as fh:
        for line in fh:
            seqid, _, taxid = line.rstrip("\n").partition("\t")
            if not taxid:
                continue
            acc = seqid.split("|", 1)[0].rsplit("_", 1)[0]
            reps.setdefault(int(taxid), acc)
    return reps


def run_kraken(
    *,
    sample: str,
    r1: Path,
    r2: Path | None,
    refs_dir: Path,
    work_dir: Path,
    threads: int = 8,
    panel: str = "uhgg",
) -> LaneResult | None:
    spec = PANELS[panel]
    kraken, bracken = _tools()
    db = refs_dir / "kraken2" / spec["dirname"]
    if kraken is None or bracken is None or not database_ready(refs_dir, panel):
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    sha = input_sha256([r1, r2])
    version = kraken_version()
    build_stamp = ""
    if panel == "rescue":
        build_stamp = str(json.loads((db / "build.json").read_text()).get("built_at", ""))
    key = cache_key(str(CACHE_VERSION), version, spec["release"], build_stamp, str(CONFIDENCE), str(MIN_HIT_GROUPS),
                    str(BRACKEN_THRESHOLD), sha)
    cached = work_dir / f"kraken.{key}.json"
    hit = read_cached(cached)
    if hit is not None:
        return hit

    t0 = time.monotonic()
    report = work_dir / f"{sample}.kraken2.minimizers.report"
    assignments = work_dir / f"{sample}.kraken2.assignments.tsv.gz"
    cmd = [kraken, "--db", str(db), "--threads", str(threads), "--confidence", str(CONFIDENCE),
           "--minimum-hit-groups", str(MIN_HIT_GROUPS), "--report-minimizer-data",
           "--report", str(report), "--output", "/dev/stdout", "--gzip-compressed"]
    if r2 is not None:
        cmd += ["--paired", str(r1), str(r2)]
    else:
        cmd.append(str(r1))
    # Per-read assignments are retained (spec §4B) and compressed on the
    # way to disk: kraken2 writes to a pipe, gzip writes the file. Passing a
    # GzipFile as stdout would hand kraken2 the raw descriptor and leave an
    # uncompressed stream under a .gz name.
    with assignments.open("wb") as out_fh:
        gz = subprocess.Popen(["gzip", "-1"], stdin=subprocess.PIPE, stdout=out_fh)
        proc = subprocess.Popen(cmd, stdout=gz.stdin, stderr=subprocess.PIPE, text=True)
        _, kraken_stderr = proc.communicate()
        gz.stdin.close()
        gz.wait()
    if proc.returncode != 0 or gz.returncode != 0 or not report.is_file():
        raise RuntimeError(f"kraken2 failed (exit {proc.returncode}): {(kraken_stderr or '')[-2000:]}")

    # Bracken wants the six-column report; the minimizer columns are kept
    # in the retained report and read below for support.
    std_report = work_dir / f"{sample}.kraken2.report"
    minimizers: dict[int, tuple[int, int]] = {}
    with report.open() as src, std_report.open("w") as dst:
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 8:
                with contextlib.suppress(ValueError):
                    minimizers[int(parts[6])] = (int(parts[3]), int(parts[4]))
                dst.write("\t".join(parts[:3] + parts[5:]) + "\n")
            else:
                dst.write(line)
    read_len = _read_length(r1)
    kmer_len = 150 if read_len >= 125 else 100
    bracken_out = work_dir / f"{sample}.bracken.species.tsv"
    bcmd = [bracken, "-d", str(db), "-i", str(std_report), "-o", str(bracken_out), "-r", str(kmer_len),
            "-l", "S", "-t", str(BRACKEN_THRESHOLD)]
    bproc = subprocess.run(bcmd, capture_output=True, text=True, check=False)
    bracken_ok = bproc.returncode == 0 and bracken_out.is_file()

    names, nodes = _taxonomy(db)
    reps = _species_reps(db)
    observations, summary = _parse(std_report, bracken_out if bracken_ok else None, minimizers, names, nodes, reps,
                                   id_prefix=spec["id_prefix"], rep_key=spec["rep_key"])
    summary.update({"read_length_mean": read_len, "bracken_kmer_distribution": f"database{kmer_len}mers",
                    "kraken_stderr_tail": kraken_stderr.strip().splitlines()[-3:] if kraken_stderr else []})
    notes = [f"operating point: --confidence {CONFIDENCE} --minimum-hit-groups {MIN_HIT_GROUPS}; Kraken confidence is "
             "not a probability", "Bracken redistributes reads to species; it is an abundance step, not a confirmation",
             f"inventory floor: >= {INVENTORY_MIN_READS:,} read pairs and >= {INVENTORY_MIN_MINIMIZERS:,} distinct minimizers; "
             "species below it are retained here and are not organisms"]
    if not bracken_ok:
        notes.append(f"Bracken did not run ({bproc.stderr.strip()[-200:]}); species read fractions are Kraken's own")
    if panel == "rescue":
        notes.append("panel: reconciliation targets with their GTDB/GlobDB near neighbours, every genome any sample's "
                     "confirmation fetched, and food/PhiX decoys; a decoy hit is reported here and is never an organism")
    result = LaneResult(
        lane_id=spec["lane_id"], tool="kraken2+bracken", tool_version=f"kraken2 {version}; bracken {_bracken_version(bracken)}",
        reference_release_id=spec["release"] + (f" built {build_stamp}" if build_stamp else ""), taxonomy_release=spec["taxonomy"],
        observations=observations, elapsed_s=time.monotonic() - t0, cached=False,
        command=tuple(cmd + ["&&"] + bcmd), input_sha256=sha, raw_result_uri=str(report), notes=notes, summary=summary,
    )
    cached.write_text(json.dumps(result.to_json(), indent=1))
    return result


def unplaced_read_ids(assignments: Path, db: Path) -> tuple[set[str], dict[str, int]]:
    """Read ids Kraken left without a species-level home.

    Unclassified reads and reads placed only above species (an LCA at
    genus or higher) are the material a novel or uncatalogued organism is
    made of; they are what targeted assembly should assemble. Returns the
    ids and the counts behind them, so the operating point is on record.
    """
    _, nodes = _taxonomy(db)
    placed_cache: dict[int, bool] = {0: False}

    def at_species_or_below(taxid: int) -> bool:
        hit = placed_cache.get(taxid)
        if hit is not None:
            return hit
        t, seen = taxid, 0
        result = False
        while t in nodes and t != 1 and seen < 40:
            parent, rank = nodes[t]
            if rank == "species":
                result = True
                break
            t = parent
            seen += 1
        placed_cache[taxid] = result
        return result

    ids: set[str] = set()
    counts = {"total": 0, "unclassified": 0, "above_species": 0, "at_species": 0}
    opener = gzip.open if assignments.suffix == ".gz" else open
    with opener(assignments, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            parts = line.split("\t", 3)
            if len(parts) < 3:
                continue
            counts["total"] += 1
            if parts[0] == "U":
                counts["unclassified"] += 1
                ids.add(parts[1])
                continue
            try:
                taxid = int(parts[2])
            except ValueError:
                continue
            if at_species_or_below(taxid):
                counts["at_species"] += 1
            else:
                counts["above_species"] += 1
                ids.add(parts[1])
    return ids, counts


def _bracken_version(exe: str) -> str:
    out = subprocess.run([exe, "-v"], capture_output=True, text=True, check=False)
    text = (out.stdout or out.stderr or "").strip()
    return text.splitlines()[0][:40] if text else ""


def _parse(
    report: Path, bracken_out: Path | None, minimizers: dict[int, tuple[int, int]],
    names: dict[int, str], nodes: dict[int, tuple[int, str]], reps: dict[int, str],
    *, id_prefix: str = "uhgg-taxid", rep_key: str = "uhgg_species_rep",
) -> tuple[list[Observation], dict[str, Any]]:
    # Kraken report: pct, clade reads, taxon reads, rank, taxid, name
    clade_reads: dict[int, int] = {}
    direct_reads: dict[int, int] = {}
    total = unclassified = 0
    with report.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            taxid = int(parts[4])
            clade_reads[taxid] = int(parts[1])
            direct_reads[taxid] = int(parts[2])
            if taxid == 0:
                unclassified = int(parts[1])
            if taxid == 1:
                total = int(parts[1])
    total_reads = total + unclassified

    bracken: dict[int, int] = {}
    if bracken_out is not None:
        with bracken_out.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                try:
                    bracken[int(parts[idx["taxonomy_id"]])] = int(parts[idx["new_est_reads"]])
                except (KeyError, ValueError, IndexError):
                    continue
    bracken_total = sum(bracken.values())

    obs: list[Observation] = []
    for taxid, (_parent, rank) in nodes.items():
        if rank != "species":
            continue
        reads = clade_reads.get(taxid, 0)
        if reads <= 0:
            continue
        est = bracken.get(taxid)
        distinct, total_min = minimizers.get(taxid, (0, 0))
        # Support: a species is provisional until confirmation when it rests
        # on few reads or few distinct minimizers - shared k-mers between
        # close species look exactly like this.
        status = "supported" if (reads >= SUPPORTED_MIN_READS and distinct >= SUPPORTED_MIN_MINIMIZERS) else "provisional"
        if est is not None and est < BRACKEN_THRESHOLD:
            status = "provisional"
        below_floor = reads < INVENTORY_MIN_READS or distinct < INVENTORY_MIN_MINIMIZERS or (est is not None and est < BRACKEN_THRESHOLD)
        obs.append(Observation(
            native_id=f"{id_prefix}:{taxid}", lineage=_lineage(taxid, names, nodes), rank="species", status=status,
            abundance_value=(100.0 * est / bracken_total) if (est is not None and bracken_total) else (
                100.0 * reads / (total or 1)),
            abundance_unit="percent of Bracken species-assigned read pairs" if est is not None else
                           "percent of Kraken-classified read pairs",
            denominator="Bracken species-assigned read pairs" if est is not None else "Kraken-classified read pairs",
            support_metrics={
                "kraken_clade_reads": reads, "kraken_direct_reads": direct_reads.get(taxid, 0),
                "bracken_est_reads": est, "distinct_minimizers": distinct, "total_minimizers": total_min,
                rep_key: reps.get(taxid, ""), "below_inventory_floor": below_floor,
            },
        ))
    obs.sort(key=lambda o: -(o.abundance_value or 0.0))
    summary = {
        "total_read_pairs": total_reads, "classified_read_pairs": total, "unclassified_read_pairs": unclassified,
        "classified_fraction": round(total / total_reads, 4) if total_reads else None,
        "n_species_with_reads": len(obs), "bracken_species_reads": bracken_total or None,
        "n_species_above_inventory_floor": sum(1 for o in obs if not o.support_metrics.get("below_inventory_floor")),
    }
    return obs, summary
