"""StrainPhlAn 4 on the Jan26 marker alignment (spec 0.8.4 §5).

What the project had was consensus-marker reconstruction (`sample2markers`
on the Jun23 SAM): a fingerprint per organism, no comparison. This module
completes the comparative workflow, on the Jan26 SAM and nothing older:

    1. sample2markers   Jan26 SAM -> per-organism consensus markers
    2. eligibility      organisms with enough reconstructed markers
    3. extract_markers  the clade's marker sequences from the Jan26 database
    4. references       member genomes of the organism's GTDB R232 species
                        (accession-versioned, fetched once) - a tree of one
                        sample is not a tree
    5. strainphlan      alignment and phylogeny of the sample's consensus
                        against those references (PhyloPhlAn inside)
    6. placement        the sample's nearest reference and its distance,
                        read from the tree; mixture evidence from the
                        consensus' polymorphic-site rate

An organism whose markers are too few is `unresolved` with the number.
A nearest reference is a placement, not proof of an identical resident
strain; a consensus can hide co-resident strains, and the polymorphism
rate is reported for that reason. Strain results never enter species
richness: they hang off an organism the inventory already lists.

Complete: the comparative step runs for every eligible organism the
inventory lists, most abundant first. Each clade's result is cached on
disk, so a run that hits the time budget records what it left for the next
run rather than pretending it was done, and the next run resumes there.
Marker clades the inventory does not list as organisms are not placed: a
strain is a property of a detected organism, never a detection.
"""

from __future__ import annotations

import bz2
import json
import re
import shutil
import subprocess
import tempfile
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from openbiota.expansion import genomes, gtdb

INDEX: Final = "mpa_vJan26_CHOCOPhlAnSGB_202605"
MIN_MARKERS: Final = 20                 # markers reconstructed for the organism to be eligible
MAX_CLADES: Final = 500                 # comparative trees per sample per run (every eligible organism)
MAX_REFERENCES: Final = 6               # member genomes per clade
TIME_BUDGET_S: Final = 4 * 60 * 60      # per run; cached clades cost nothing, so runs converge
CACHE_VERSION: Final = 1


def _venv_bin(name: str) -> Path | None:
    root = Path(__file__).resolve().parents[2]
    exe = root / ".venv-mpa42" / "bin" / name
    return exe if exe.is_file() else None


def _path_env() -> dict[str, str]:
    """PATH with the project's vendor/bin first: raxmlHPC, blastn, makeblastdb live there."""
    import os
    root = Path(__file__).resolve().parents[2]
    venv_bin = root / ".venv-mpa42" / "bin"
    env = dict(os.environ)
    env["PATH"] = f"{venv_bin}:{root / 'vendor' / 'bin'}:/usr/local/bin:{env.get('PATH', '')}"
    return env


def externals_present() -> dict[str, bool]:
    path = _path_env()["PATH"]
    return {t: shutil.which(t, path=path) is not None for t in ("blastn", "makeblastdb", "mafft", "trimal", "raxmlHPC", "FastTree")}


def available() -> bool:
    ext = externals_present()
    return all(_venv_bin(n) for n in ("sample2markers.py", "extract_markers.py", "strainphlan")) and (
        ext["blastn"] and ext["makeblastdb"] and ext["mafft"] and ext["trimal"] and ext["raxmlHPC"])


def run_sample2markers(*, sample: str, sam: Path, db_dir: Path, work_dir: Path, threads: int) -> Path:  # noqa: ARG001 - `sample` names the run in logs and cache paths of callers
    """Jan26 consensus markers for the sample; cached on the SAM's size."""
    exe = _venv_bin("sample2markers.py")
    if exe is None:
        raise RuntimeError("sample2markers.py not found in .venv-mpa42")
    out_dir = work_dir / "jan26_consensus"
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = out_dir / "stamp.json"
    existing = sorted(out_dir.glob("*.json.bz2")) + sorted(out_dir.glob("*.pkl"))
    if existing and stamp.is_file() and json.loads(stamp.read_text()).get("sam_size") == sam.stat().st_size:
        return existing[0]
    pkl = db_dir / f"{INDEX}.pkl"
    cmd = [str(exe), "-i", str(sam), "-o", str(out_dir), "-d", str(pkl), "-n", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=_path_env())
    produced = sorted(out_dir.glob("*.json.bz2")) + sorted(out_dir.glob("*.pkl"))
    if proc.returncode != 0 or not produced:
        raise RuntimeError(f"sample2markers (Jan26) failed: {proc.stderr[-1500:]}")
    stamp.write_text(json.dumps({"sam_size": sam.stat().st_size, "cache_version": CACHE_VERSION, "command": cmd}))
    return produced[0]


def _consensus_markers(path: Path) -> dict[str, dict[str, Any]]:
    """clade -> {n_markers, polymorphic_rate} from a sample2markers output."""
    if path.suffix == ".bz2":
        with bz2.open(path, "rt") as fh:
            data = json.load(fh)
    else:
        import pickle
        with path.open("rb") as fh:
            data = pickle.load(fh)  # noqa: S301 - our own tool's output
    out: dict[str, dict[str, Any]] = {}
    markers = data.get("consensus_markers") if isinstance(data, dict) else data
    for m in markers or []:
        name = str(m.get("marker") or m.get("name") or "")
        # marker ids embed the SGB: e.g. "SGB4837__..." or "..._SGB4837_..."
        mt = re.search(r"(SGB\d+)", name)
        if not mt:
            continue
        clade = mt.group(1)
        seq = str(m.get("sequence") or "")
        rec = out.setdefault(clade, {"n_markers": 0, "bases": 0, "ambiguous": 0})
        rec["n_markers"] += 1
        rec["bases"] += len(seq)
        rec["ambiguous"] += sum(1 for ch in seq if ch not in "ACGTacgt-")
    for rec in out.values():
        rec["polymorphic_rate"] = round(rec["ambiguous"] / rec["bases"], 5) if rec["bases"] else None
    return out


def _reference_accessions(sgb: str) -> list[str]:
    """Member genomes of the SGB's GTDB R232 species: representative first."""
    from openbiota import gtdb as _g
    hit = _g.lookup(sgb, _g.JAN26)
    if hit is None or not hit.resolved:
        return []
    clusters = gtdb.species_clusters("r232")
    c = clusters.get(hit.species)
    if c is None:
        return []
    members = [c.representative] + [m for m in c.members if m != c.representative]
    return members[:MAX_REFERENCES]


def _nearest_from_tree(tree_path: Path, sample: str) -> tuple[str | None, float | None]:
    """Nearest reference leaf to the sample and the branch-length distance."""
    try:
        import dendropy
    except ImportError:
        return None, None
    try:
        tree = dendropy.Tree.get(path=str(tree_path), schema="newick", preserve_underscores=True)
    except Exception:  # noqa: BLE001
        return None, None
    pdm = tree.phylogenetic_distance_matrix()
    leaves = {t.label: t for t in tree.taxon_namespace}
    me = next((t for lbl, t in leaves.items() if lbl and sample in lbl), None)
    if me is None:
        return None, None
    best, dist = None, None
    for lbl, t in leaves.items():
        if t is me or not lbl:
            continue
        d = pdm.distance(me, t)
        if dist is None or d < dist:
            best, dist = lbl, d
    return best, (round(float(dist), 6) if dist is not None else None)


def _markers2clade(db_dir: Path) -> Path:
    """marker name -> clade, as a TSV derived once from the database pickle.

    The pickle holds five million marker records and takes minutes to load;
    the two columns StrainPhlAn extraction needs load in seconds.
    """
    tsv = db_dir / f"{INDEX}.markers2clade.tsv"
    if tsv.is_file() and tsv.stat().st_size > 0:
        return tsv
    import bz2
    import pickle

    pkl = db_dir / f"{INDEX}.pkl"
    with bz2.open(pkl, "rb") as fh:
        data = pickle.load(fh)  # noqa: S301 - MetaPhlAn's own database
    tmp = tsv.with_suffix(".tmp")
    with tmp.open("w") as out:
        for name, info in data["markers"].items():
            out.write(f"{name}\t{info.get('clade', '')}\n")
    tmp.replace(tsv)
    return tsv


def _extract_clade_markers(*, db_dir: Path, clades: list[str], markers_dir: Path, threads: int) -> None:
    """Write t__<SGB>.fna for every clade in one pass over the marker FASTA.

    MetaPhlAn's extract_markers.py re-dumps the whole 21 GB Bowtie2 index
    on every call and then parses the dump once per clade; a dozen clades
    cost an hour. The dump already exists beside the index (the pipeline
    keeps it), the marker->clade map is a TSV, and seqkit pulls every
    requested marker in a single streaming pass. The output is identical:
    the clade's marker sequences under their database names.
    """
    dump = db_dir / f"{INDEX}.fna"
    if not dump.is_file() or dump.stat().st_size < 1_000_000_000:
        raise FileNotFoundError(f"marker FASTA dump missing: {dump}")
    seqkit = shutil.which("seqkit")
    if seqkit is None:
        raise FileNotFoundError("seqkit not installed")
    wanted = {f"t__{c}" for c in clades}
    clade_of: dict[str, str] = {}
    with _markers2clade(db_dir).open() as fh:
        for line in fh:
            name, _, clade = line.rstrip("\n").partition("\t")
            if clade in wanted:
                clade_of[name] = clade
    if not clade_of:
        raise ValueError(f"no markers in the database for {sorted(wanted)}")
    with tempfile.TemporaryDirectory(dir=markers_dir) as td:
        ids = Path(td) / "ids.txt"
        ids.write_text("\n".join(clade_of) + "\n")
        pulled = Path(td) / "pulled.fna"
        subprocess.run([seqkit, "grep", "-j", str(max(2, min(threads, 12))), "-f", str(ids), "-o", str(pulled), str(dump)],
                       check=True, capture_output=True, text=True)
        handles: dict[str, Any] = {}
        try:
            current = None
            with pulled.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        clade = clade_of.get(name)
                        current = None
                        if clade is not None:
                            if clade not in handles:
                                handles[clade] = (markers_dir / f"{clade}.fna.tmp").open("w")
                            current = handles[clade]
                    if current is not None:
                        current.write(line)
        finally:
            for h in handles.values():
                h.close()
        for clade in wanted:
            tmp = markers_dir / f"{clade}.fna.tmp"
            if tmp.is_file() and tmp.stat().st_size > 0:
                tmp.replace(markers_dir / f"{clade}.fna")
            elif tmp.is_file():
                tmp.unlink()


def attach_to_inventory(blob: dict[str, Any], analysis: Mapping[str, Any] | None) -> int:
    """Write each placed clade onto its organism's `strain` record.

    An organism the marker-typing lane already resolved keeps that record
    and gains the placement fields beside it; one without a record gets the
    placement as its record. Only clades that were placed. Returns how many
    were attached. Never changes counts of organisms; a strain is a property
    of an organism that is already there.
    """
    if not analysis or analysis.get("status") != "completed":
        return 0
    placed = {c["sgb"]: c for c in analysis.get("clades") or [] if c.get("status") == "placed" and c.get("sgb")}
    if not placed:
        return 0
    n = 0
    for rec in blob.get("organisms") or []:
        sgb = str(rec.get("sgb") or "")
        if sgb not in placed:
            continue
        c = placed[sgb]
        existing = dict(rec.get("strain") or {})
        if existing.get("nearest_reference"):
            continue
        placement = {
            "placement_method": "StrainPhlAn 4 (Jan26 markers)", "nearest_reference": c.get("nearest_reference"),
            "distance_to_nearest": c.get("distance_to_nearest"), "n_markers": c.get("n_markers"),
            "polymorphic_rate": c.get("polymorphic_rate"), "references": list(c.get("references") or []),
            "mixture_evidence": c.get("mixture_evidence"), "note": c.get("note"),
        }
        rec["strain"] = {**existing, **placement} if existing else {"method": "StrainPhlAn 4 (Jan26 markers)", **placement}
        n += 1
    counts = blob.get("counts")
    if isinstance(counts, dict):
        counts["strain_resolved"] = sum(1 for r in blob.get("organisms") or [] if r.get("strain"))
    return n


def run(*, sample: str, sam: Path, db_dir: Path, work_dir: Path, threads: int, inventory: Any) -> dict[str, Any]:
    t0 = time.monotonic()
    if not available():
        return {"status": "not_assessed", "reason": "StrainPhlAn externals missing: "
                + ", ".join(k for k, v in externals_present().items() if not v), "clades": []}
    consensus = run_sample2markers(sample=sample, sam=sam, db_dir=db_dir, work_dir=work_dir, threads=threads)
    markers = _consensus_markers(consensus)
    # Eligibility by marker count, then by abundance in the inventory.
    abundance = {}
    for o in getattr(inventory, "organisms", []) or []:
        if o.sgb and o.sgb_db == INDEX:
            abundance[o.sgb] = o.best_percent
    eligible = sorted((c for c, rec in markers.items() if rec["n_markers"] >= MIN_MARKERS and c in abundance),
                      key=lambda c: (-(abundance.get(c, 0.0)), -markers[c]["n_markers"]))
    unresolved = [{"sgb": c, "n_markers": rec["n_markers"], "reason": f"{rec['n_markers']} markers < {MIN_MARKERS}"}
                  for c, rec in markers.items() if rec["n_markers"] < MIN_MARKERS]
    # Marker clades with reads but no organism in the inventory (the
    # profiler did not call them) are listed, not placed.
    not_in_inventory = sorted(c for c, rec in markers.items() if rec["n_markers"] >= MIN_MARKERS and c not in abundance)
    clades: list[dict[str, Any]] = []
    deferred: list[str] = []
    strain_dir = work_dir / "strainphlan_jan26"
    strain_dir.mkdir(parents=True, exist_ok=True)
    extract = _venv_bin("extract_markers.py")
    strainphlan = _venv_bin("strainphlan")
    pkl = db_dir / f"{INDEX}.pkl"
    todo = [c for c in eligible[:MAX_CLADES] if not (strain_dir / c / "result.json").is_file()]
    deferred.extend(eligible[MAX_CLADES:])
    # One extraction pass for every clade, into a cache shared by every
    # sample: extract_markers decompresses the whole marker database to a
    # temporary FASTA each time it runs (about 20 GB for Jan26), and a
    # clade's markers do not depend on the sample. Once an SGB has been
    # extracted for any sample it is never extracted again.
    markers_dir = db_dir / f"clade_markers_{INDEX}"
    markers_dir.mkdir(exist_ok=True)
    missing = [c for c in todo if not (markers_dir / f"t__{c}.fna").is_file()]
    extraction_error = ""
    if missing:
        try:
            _extract_clade_markers(db_dir=db_dir, clades=missing, markers_dir=markers_dir, threads=threads)
        except Exception as exc:  # noqa: BLE001 - fall back to MetaPhlAn's own extractor
            extraction_error = f"{type(exc).__name__}: {str(exc)[:200]}"
        still = [c for c in missing if not (markers_dir / f"t__{c}.fna").is_file()]
        if still:
            p = subprocess.run([str(extract), "-c", *[f"t__{c}" for c in still], "-d", str(pkl), "-o", str(markers_dir)],
                               capture_output=True, text=True, check=False, env=_path_env())
            extraction_error = p.stderr[-300:] if p.returncode != 0 else extraction_error
    # Reference genomes for every clade in this run are fetched together,
    # eight at a time: the downloads are network-bound and independent, and
    # fetching them one clade at a time put ten minutes of idle waiting in
    # front of every sample's trees. Failures are simply absent from the
    # cache; the per-clade loop below counts what it can use.
    wanted_accs = {acc for c in todo for acc in _reference_accessions(c)}
    if wanted_accs:
        from concurrent.futures import ThreadPoolExecutor

        def _prefetch(acc: str) -> None:
            try:
                genomes.fetch(acc)
            except Exception:  # noqa: BLE001 - recorded per clade below
                return

        with ThreadPoolExecutor(max_workers=4) as pool:  # NCBI rate-limits beyond a few connections
            list(pool.map(_prefetch, sorted(wanted_accs)))
    # The budget is for the comparative trees. Extraction and genome
    # fetching are cache fills and must not eat it.
    t_budget = time.monotonic()
    for sgb in eligible[:MAX_CLADES]:
        clade_dir = strain_dir / sgb
        clade_dir.mkdir(exist_ok=True)
        done = clade_dir / "result.json"
        if done.is_file():
            clades.append(json.loads(done.read_text()))
            continue
        if time.monotonic() - t_budget > TIME_BUDGET_S:
            deferred.append(sgb)
            continue
        rec: dict[str, Any] = {"sgb": sgb, "n_markers": markers[sgb]["n_markers"],
                               "polymorphic_rate": markers[sgb]["polymorphic_rate"], "abundance": abundance.get(sgb)}
        marker_fna = markers_dir / f"t__{sgb}.fna"
        if not marker_fna.is_file():
            rec.update({"status": "unresolved", "reason": f"marker extraction failed: {extraction_error or 'no markers written'}"})
            clades.append(rec)
            continue  # not cached: an extraction failure is retried next run
        # references
        refs: list[Path] = []
        accs = _reference_accessions(sgb)
        for acc in accs:
            try:
                refs.append(genomes.fetch(acc))
            except Exception:  # noqa: BLE001
                continue
        rec["references"] = [r.name for r in refs]
        if len(refs) < 3:
            rec.update({"status": "unresolved", "reason": f"only {len(refs)} reference genomes could be fetched "
                        f"({len(accs)} members in GTDB R232); a comparative tree needs at least three"})
            clades.append(rec)
            if len(accs) < 3:
                done.write_text(json.dumps(rec))  # the species has too few genomes: a fact, cached
            continue  # a fetch that failed is retried next run, not cached
        # strainphlan
        out_dir = clade_dir / "out"
        out_dir.mkdir(exist_ok=True)
        cmd = [str(strainphlan), "-d", str(pkl), "-s", str(consensus), "-m", str(marker_fna), "-o", str(out_dir), "-c", f"t__{sgb}",
               "-n", str(threads), "--mutation_rates", "--non_interactive", "--phylophlan_mode", "fast",
               "--marker_in_n_samples_perc", "50", "--sample_with_n_markers", "10",
               "--sample_with_n_markers_perc", "10",
               "--sample_with_n_markers_after_filt", "10", "--sample_with_n_markers_after_filt_perc", "10"]
        # `-r` takes every reference after one flag (nargs='+'); repeating
        # the flag keeps only the last genome and StrainPhlAn then stops
        # with "samples + references are less than 4".
        cmd += ["-r", *[str(r) for r in refs]]
        p = subprocess.run(cmd, capture_output=True, text=True, check=False, env=_path_env())
        tree = next(iter(out_dir.glob("*.tre")), None) or next(iter(out_dir.glob("*.nwk")), None) \
            or next(iter(out_dir.glob("RAxML_bestTree*")), None)
        if p.returncode != 0 or tree is None:
            # A tool failure is not a result: recorded for this run, never
            # cached, so the next run tries again.
            rec.update({"status": "unresolved", "reason": f"strainphlan did not produce a tree: {(p.stderr or p.stdout)[-400:]}"})
            clades.append(rec)
            continue
        nearest, dist = _nearest_from_tree(tree, sample)
        if nearest is None:
            info = next(iter(out_dir.glob("*.info")), None)
            kept = ""
            if info is not None:
                m = re.search(r"Number of samples after filtering:\s*(\d+)", info.read_text())
                kept = m.group(1) if m else ""
            rec.update({"status": "unresolved", "tree": str(tree),
                        "reason": ("the sample's consensus did not pass StrainPhlAn's marker filters "
                                   f"(samples kept after filtering: {kept or 'none'}); the references alone form the tree")})
            clades.append(rec)
            done.write_text(json.dumps(rec))
            continue
        rec.update({
            "status": "placed",
            "tree": str(tree), "nearest_reference": nearest, "distance_to_nearest": dist,
            "mixture_evidence": ("polymorphic sites at " + f"{markers[sgb]['polymorphic_rate']:.2%}" +
                                 "; a consensus can hide co-resident strains") if markers[sgb]["polymorphic_rate"] else "none measured",
            "note": "a nearest reference is a placement, not proof of an identical resident strain",
        })
        clades.append(rec)
        done.write_text(json.dumps(rec))
    return {
        "status": "completed", "database": INDEX, "consensus": str(consensus), "n_clades_with_markers": len(markers),
        "n_eligible": len(eligible), "n_compared": sum(1 for c in clades if c.get("status") in ("placed", "tree_built")),
        "n_unresolved": len(unresolved) + sum(1 for c in clades if c.get("status") == "unresolved"),
        "clades": clades, "unresolved": unresolved, "deferred": deferred, "not_in_inventory": not_in_inventory,
        "operating_points": {"min_markers": MIN_MARKERS, "max_clades_per_run": MAX_CLADES,
                             "max_references": MAX_REFERENCES, "time_budget_s": TIME_BUDGET_S},
        "elapsed_s": round(time.monotonic() - t0, 1),
        "meaning": ("Comparative placement of this sample's consensus markers against member genomes of the "
                    "organism's GTDB R232 species. placed: nearest reference and distance read from the tree; "
                    "unresolved: too few markers or references. Strain results never add to species counts."),
    }
