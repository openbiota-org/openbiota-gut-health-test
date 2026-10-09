"""Competitive confirmation of new and confusable organisms (spec 0.8.4 §4C).

A profiler's call is a claim. This stage tests the claims that need it -
every organism the installed baseline never saw, and every one resting on
a single method - by mapping the sample's reads to the candidate genome
*and* its closest detected relatives together, so a read has to choose.
What comes back per candidate:

    unique_fragments   read pairs that mapped to this genome with MAPQ >= 20
    breadth            fraction of the genome covered by at least one read
    mean_depth         mean depth over covered bases
    identity           1 - mismatches / aligned bases, over unique fragments
    evenness           observed breadth / the breadth a random spread of the
                       same reads would give (Poisson, 1 - exp(-mean depth));
                       1 = reads spread over the genome, near 0 = piled into
                       a few regions (a shared operon, a mobile element)
    competitors        the alternatives' unique fragments, so a reader sees
                       what the reads chose against

Operating points are constants here, recorded in every result, and set so
that a handful of reads on a mobile element or a conserved rRNA cannot
establish a species: a supported call needs many independent fragments,
real breadth, species-level identity and coverage spread over the genome.
They are starting points to be calibrated by `scripts/benchmark_lanes.py`,
not tuned after looking at a sample.

A candidate whose reference genome cannot be fetched yet stays provisional
with `pending_archive` as the reason. Nothing is promoted by default.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
import time
import zlib
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.expansion import genomes

MIN_UNIQUE_FRAGMENTS: Final = 100
MIN_BREADTH: Final = 0.03
MIN_IDENTITY: Final = 0.97
MIN_EVENNESS: Final = 0.25
PROVISIONAL_FRAGMENTS: Final = 20
MAPQ_UNIQUE: Final = 20
MAX_CANDIDATES: Final = 400   # every candidate the spec names is tested; the cap is a safety, not a budget
#: A same-genus call below this fraction of the genus's largest reading is
#: tested as a possible shadow of that larger organism.
SHADOW_FRACTION: Final = 0.20   # a single-method call under a fifth of its relative's unique reads, at lower identity, is the relative's DNA
IDENTITY_GAP: Final = 0.003     # identity below the relative's by this much marks reads that belong to the relative
MAX_ALTERNATIVES: Final = 6
CACHE_VERSION: Final = 5


@dataclass
class Verdict:
    organism: str
    genome_id: str
    status: str                       # supported | provisional | not_detected | pending
    reason: str
    unique_fragments: int = 0
    total_fragments: int = 0
    breadth: float = 0.0
    mean_depth: float = 0.0
    identity: float = 0.0
    evenness: float = 0.0
    genome_length: int = 0
    competitors: dict[str, int] = field(default_factory=dict)
    previous_status: str = ""
    #: For a rejected call: the relative the reads belong to (its share goes there).
    relative: str = ""

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        for k in ("breadth", "mean_depth", "identity", "evenness"):
            d[k] = round(d[k], 4)
        return d


def _evenness(breadth: float, mean_depth: float) -> float:
    """How evenly the reads are spread: observed breadth over the breadth the
    same mean depth would give if reads fell at random (1 - e^-depth).

    A genuine organism at strain distance from the reference covers the
    genome broadly and scores near 1 even where divergent regions drop out
    (0.5 or more in practice); a shadow - reads of a relative landing on the
    regions the two genomes share - reaches a fraction of the genome at a
    depth that would have covered most of it, and scores near 0. The
    previous definition (1 - CV of per-contig depth) scored every
    fragmented draft genome near 0 whatever the reads did.
    """
    if mean_depth <= 0 or breadth <= 0:
        return 0.0
    import math

    expected = 1.0 - math.exp(-mean_depth)
    return min(1.0, breadth / expected) if expected > 0 else 0.0


def _tools() -> tuple[str | None, str | None, str | None]:
    return shutil.which("bowtie2"), shutil.which("bowtie2-build"), shutil.which("samtools")


def available() -> bool:
    return all(_tools())


def select_candidates(inventory: Any) -> list[Any]:
    """Who gets tested (spec §4C): every newly added organism, every
    provisional or ambiguous call, and every confusable one - a call resting
    on a single method with a relative of the same genus also detected, so
    that sister species compete for the reads that could belong to either."""
    genera: dict[str, int] = {}
    for o in inventory.organisms:
        g = (o.gtdb_genus or o.genus or "").split("_")[0]
        if g:
            genera[g] = genera.get(g, 0) + 1
    # The largest reading per genus: a same-genus call at a small fraction
    # of it is a shadow candidate - reads of the dominant species spilling
    # onto a sibling's shared markers or k-mers look exactly like a minor
    # sibling, and only competitive mapping tells them apart.
    top_in_genus: dict[str, float] = {}
    for o in inventory.organisms:
        g = (o.gtdb_genus or o.genus or "").split("_")[0]
        if g:
            top_in_genus[g] = max(top_in_genus.get(g, 0.0), o.best_percent or 0.0)
    cands = []
    for o in inventory.organisms:
        g = (o.gtdb_genus or o.genus or "").split("_")[0]
        confusable = len(getattr(o, "methods", ()) or ()) <= 1 and genera.get(g, 0) >= 2
        shadow = (genera.get(g, 0) >= 2 and len(getattr(o, "methods", ()) or ()) <= 1
                  and (o.best_percent or 0.0) < SHADOW_FRACTION * top_in_genus.get(g, 0.0))
        # An organism the primary lane did not place needs the competition
        # whatever its support: its share of the genus is apportioned from
        # the reads that map uniquely to its genome (inventory.unify_shares).
        apportion = not o.in_primary and (o.secondary_percent or 0.0) > 0
        if o.incremental_gain == "expansion" or o.status in ("provisional", "ambiguous") or confusable or shadow or apportion:
            cands.append(o)
    cands.sort(key=lambda o: -(o.best_percent or 0.0))
    return cands[:MAX_CANDIDATES]


_REPS: dict[str, str] | None = None


def _r232_representative(species: str) -> str | None:
    global _REPS
    if _REPS is None:
        try:
            from openbiota.expansion import gtdb as _g
            _REPS = _g.representative_accessions("r232")
        except Exception:  # noqa: BLE001
            _REPS = {}
    return _REPS.get(species)


def _fetchable(gid: str) -> bool:
    if gid.startswith("file:"):
        return genomes.available(gid) == "fetchable"
    return bool(gid) and (
        gid.startswith(("GCA_", "GCF_", "MGYG", "HRGMV2_", "HRGMv2_", "HROM_Genome_")) or genomes.in_globdb(gid))


def _genome_id_for(o: Any) -> str | None:
    """The identifier whose genome represents this organism.

    The lanes' own genome ids come first - the GlobDB genome behind a
    cluster, the panel genome a Kraken hit was classified to, the mOTUs
    cluster's representative - then anything a native identifier names
    that a source can fetch. Anything resolved to a GTDB R232 species is
    represented by that species' representative genome, accession-versioned
    from NCBI - which makes every named call testable, sister species
    included.
    """
    for gid in getattr(o, "genome_ids", ()) or ():
        if _fetchable(gid):
            return gid
    for lane in ("globdb", "genome", "jan26", "extended", "motus", "rescue", "kraken", "singlem"):
        nid = (o.native_ids or {}).get(lane)
        if not nid:
            continue
        for token in str(nid).split(";"):
            gid = genomes.normalise(token)
            if gid.startswith(("uhgg-taxid", "rescue-taxid", "mOTUv4", "SGB")):
                continue
            if _fetchable(gid):
                return gid
    species = (getattr(o, "gtdb", None) or "").strip()
    if species and " sp. (" not in species and not species.endswith(")"):
        rep = _r232_representative(species)
        if rep:
            return rep
    return None


def _alternatives(o: Any, inventory: Any) -> list[tuple[str, str]]:
    """Detected relatives in the same genus with a fetchable genome: (organism, genome id)."""
    out: list[tuple[str, str]] = []
    genus = (o.gtdb_genus or o.genus or "").split("_")[0]
    if not genus:
        return out
    # the largest relatives first: they are the ones whose reads a shadow is made of
    relatives = sorted((x for x in inventory.organisms if x is not o
                        and (x.gtdb_genus or x.genus or "").split("_")[0] == genus),
                       key=lambda x: -(x.best_percent or 0.0))
    for other in relatives:
        gid = _genome_id_for(other)
        if gid and genomes.available(gid) in ("cached", "fetchable"):
            out.append((other.display, gid))
        if len(out) >= MAX_ALTERNATIVES:
            break
    return out


def _contig_map(fasta_gz: Path, tag: str) -> tuple[list[str], int, list[tuple[str, int]]]:
    """Contig names retagged with the genome id, so one index holds many genomes."""
    contigs: list[tuple[str, int]] = []
    out_lines: list[str] = []
    name = None
    length = 0
    with gzip.open(fasta_gz, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, length))
                name = f"{tag}|{line[1:].split()[0]}"
                length = 0
                out_lines.append(f">{name}\n")
            else:
                seq = line.strip()
                length += len(seq)
                out_lines.append(seq + "\n")
    if name is not None:
        contigs.append((name, length))
    return out_lines, sum(n for _, n in contigs), contigs


def run_confirmation(
    *,
    sample: str,
    inventory: Any,
    r1: Path,
    r2: Path | None,
    work_dir: Path,
    threads: int = 8,
) -> dict[str, Any]:
    bowtie2, build, samtools = _tools()
    if not (bowtie2 and build and samtools):
        return {"status": "not_assessed", "reason": "bowtie2/samtools not installed", "verdicts": []}
    work_dir.mkdir(parents=True, exist_ok=True)
    candidates = select_candidates(inventory)
    verdicts: list[Verdict] = []
    if not candidates:
        return {"status": "completed", "n_candidates": 0, "verdicts": [], "operating_points": _points()}

    # Gather genomes: candidates and their alternatives, one combined reference.
    t0 = time.monotonic()
    genome_of: dict[str, str] = {}          # tag -> organism display
    pending: list[tuple[Any, str]] = []
    fasta_lines: list[str] = []
    lengths: dict[str, int] = {}
    contig_genome: dict[str, str] = {}
    contig_len: dict[str, int] = {}
    roles: dict[str, str] = {}

    def add(gid: str, organism: str, role: str) -> bool:
        if gid in lengths:
            return True
        try:
            path = genomes.fetch(gid)
        except genomes.Pending:
            return False
        except Exception:  # noqa: BLE001 - one unfetchable genome must not stop the stage
            return False
        try:
            lines, total, contigs = _contig_map(path, gid)
        except (OSError, EOFError, ValueError, zlib.error) as exc:  # a cached copy that is not a readable gzip
            genomes.discard(gid)
            try:
                lines, total, contigs = _contig_map(genomes.fetch(gid), gid)
            except Exception:  # noqa: BLE001 - refilled once; still unreadable means the source is
                return False
            del exc
        fasta_lines.extend(lines)
        lengths[gid] = total
        genome_of[gid] = organism
        roles[gid] = role
        for cname, clen in contigs:
            contig_genome[cname] = gid
            contig_len[cname] = clen
        return True

    tested: list[tuple[Any, str]] = []
    for o in candidates:
        gid = _genome_id_for(o)
        if gid is None:
            verdicts.append(Verdict(o.display, "", "pending", "no reference genome identifier for this organism",
                                    previous_status=o.status))
            continue
        if not add(gid, o.display, "candidate"):
            verdicts.append(Verdict(o.display, gid, "pending", "reference genome pending (archive not yet acquired)",
                                    previous_status=o.status))
            pending.append((o, gid))
            continue
        tested.append((o, gid))
        for alt_name, alt_gid in _alternatives(o, inventory):
            add(alt_gid, alt_name, "alternative")

    if not tested:
        return {"status": "completed", "n_candidates": len(candidates), "n_tested": 0,
                "verdicts": [v.to_json() for v in verdicts], "operating_points": _points(),
                "elapsed_s": round(time.monotonic() - t0, 1)}

    # The competition's answer depends only on the reads, the genomes in the
    # index and the operating points. A rerun with the same three (a report
    # change, a merge that did not touch this sample's candidates) reuses the
    # previous alignment's result instead of aligning for twenty minutes.
    import hashlib

    def _stat(p: Path | None) -> str:
        # name and size, not mtime: the host-filtered read files are rewritten
        # by every run with the same content, and a fresh timestamp must not
        # look like new reads
        try:
            return f"{p.name}:{p.stat().st_size}" if p else "-"
        except OSError:
            return str(p)

    fingerprint = hashlib.sha256("\n".join([
        f"v{CACHE_VERSION}", json.dumps(_points(), sort_keys=True), _stat(r1), _stat(r2),
        *sorted(f"{gid}\t{genome_of[gid]}\t{roles[gid]}\t{lengths[gid]}" for gid in lengths),
        *sorted(f"cand\t{o.display}\t{gid}\t{o.status}\t{len(o.methods or ())}" for o, gid in tested),
    ]).encode()).hexdigest()[:16]
    cache_inputs = {
        "reads": [_stat(r1), _stat(r2)], "n_genomes": len(lengths), "n_candidates": len(tested),
        "genomes_sha": hashlib.sha256("\n".join(sorted(lengths)).encode()).hexdigest()[:12],
        "candidates_sha": hashlib.sha256("\n".join(sorted(f"{o.display}\t{gid}\t{o.status}\t{len(o.methods or ())}"
                                                             for o, gid in tested)).encode()).hexdigest()[:12],
        "roles_sha": hashlib.sha256("\n".join(sorted(f"{gid}\t{genome_of[gid]}\t{roles[gid]}" for gid in lengths)).encode()).hexdigest()[:12],
    }
    cache_file = work_dir / f"{sample}.confirmation.{fingerprint}.json"
    if not cache_file.is_file():
        # The fingerprint names the read files. A sample whose reads were
        # renamed (same bytes, same size) with the same genomes and the same
        # candidates is the same competition: adopt that cache under the new
        # name rather than mapping eight million pairs again.
        for other in sorted(work_dir.glob("*.confirmation.*.json")):
            try:
                prior = json.loads(other.read_text()).get("cache_inputs") or {}
            except (OSError, ValueError):
                continue
            same = (prior.get("genomes_sha") == cache_inputs["genomes_sha"]
                    and prior.get("candidates_sha") == cache_inputs["candidates_sha"]
                    and prior.get("roles_sha") == cache_inputs["roles_sha"]
                    and [x.rsplit(":", 1)[-1] for x in prior.get("reads", [])]
                    == [x.rsplit(":", 1)[-1] for x in cache_inputs["reads"]])
            if same:
                other.rename(cache_file)
                break
    if cache_file.is_file():
        try:
            cached = json.loads(cache_file.read_text())
        except (OSError, ValueError):
            cached = None
        if cached and cached.get("status") == "completed" and cached.get("genomes"):
            cached["verdicts"] = [v.to_json() for v in verdicts if v.status == "pending"] + [
                v for v in cached["verdicts"] if v.get("status") != "pending"]
            cached["served_from_cache"] = True
            cached["elapsed_s"] = round(time.monotonic() - t0, 1)
            (work_dir / f"{sample}.confirmation.json").write_text(json.dumps(cached, indent=1))
            return cached

    ref = work_dir / f"{sample}.confirm.fna"
    ref.write_text("".join(fasta_lines))
    index = work_dir / f"{sample}.confirm"
    subprocess.run([build, "--threads", str(threads), "-q", str(ref), str(index)], check=True, capture_output=True)

    # Align competitively; keep only mapped reads; stream through samtools.
    bam = work_dir / f"{sample}.confirm.bam"
    cmd = [bowtie2, "-p", str(threads), "--no-unal", "-x", str(index)]
    cmd += ["-1", str(r1), "-2", str(r2)] if r2 is not None else ["-U", str(r1)]
    with bam.open("wb") as out:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen([samtools, "sort", "-@", str(max(2, threads // 4)), "-o", "-"], stdin=p1.stdout, stdout=out, stderr=subprocess.PIPE)
        p1.stdout.close()  # type: ignore[union-attr]
        p2.communicate()
        p1.wait()
    subprocess.run([samtools, "index", str(bam)], check=True, capture_output=True)

    # Per-contig coverage, then per-genome aggregation.
    cov = subprocess.run([samtools, "coverage", str(bam)], capture_output=True, text=True, check=True).stdout
    per_genome_cov: dict[str, list[tuple[int, int, float, int]]] = defaultdict(list)  # (len, covbases, meandepth, reads)
    for line in cov.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        cname = parts[0]
        gid = contig_genome.get(cname)
        if gid is None:
            continue
        end, covbases, meandepth, nreads = int(parts[2]), int(parts[4]), float(parts[6]), int(parts[3])
        per_genome_cov[gid].append((end, covbases, meandepth, nreads))

    # Unique fragments and identity from the alignments themselves.
    view = subprocess.Popen([samtools, "view", "-q", str(MAPQ_UNIQUE), "-F", "0x904", str(bam)],
                            stdout=subprocess.PIPE, text=True)
    uniq: dict[str, int] = defaultdict(int)
    mism: dict[str, int] = defaultdict(int)
    alen: dict[str, int] = defaultdict(int)
    assert view.stdout is not None
    for line in view.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        gid = contig_genome.get(parts[2])
        if gid is None:
            continue
        uniq[gid] += 1
        nm = 0
        for tag in parts[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag[5:])
                break
        mism[gid] += nm
        alen[gid] += len(parts[9])
    view.wait()
    total_by_genome: dict[str, int] = {gid: sum(r for _, _, _, r in rows) for gid, rows in per_genome_cov.items()}
    identity_of: dict[str, float] = {g: (1.0 - mism[g] / alen[g]) if alen[g] else 0.0 for g in lengths}

    def _genus(name: str) -> str:
        return name.split(" ")[0].split("_")[0]

    for o, gid in tested:
        rows = per_genome_cov.get(gid, [])
        glen = lengths.get(gid, 0) or 1
        covbases = sum(c for _, c, _, _ in rows)
        breadth = covbases / glen
        mean_depth = (sum(d * L for L, _c, d, _n in rows) / glen) if rows else 0.0
        evenness = _evenness(breadth, mean_depth)
        identity = identity_of.get(gid, 0.0)
        u = uniq[gid]
        my_genus = _genus(o.gtdb_genus or o.genus or o.display)
        competitors = {genome_of[g]: uniq[g] for g in lengths if g != gid and _genus(genome_of[g]) == my_genus}
        # The relative test (spec §4C): a candidate is not the organism when a
        # relative in the same competition took more reads at higher identity
        # - the reads had both genomes to choose from and chose the other.
        beaten_by = None
        for g in lengths:
            if g == gid or _genus(genome_of[g]) != my_genus:
                continue
            if uniq[g] > u and identity_of.get(g, 0.0) > identity + 0.002 and (beaten_by is None or uniq[g] > beaten_by[1]):
                beaten_by = (genome_of[g], uniq[g], identity_of[g])
        # Rejected when the relative both out-recruited the candidate and did
        # so at higher identity - and the candidate looks like a shadow, not a
        # minor species. A shadow is reads of the dominant organism landing on
        # a sibling's shared regions: few fragments relative to the sibling, a
        # lower identity, and coverage piled into patches (low evenness). A
        # genuine minor species covers its genome evenly even at a strain
        # distance from the reference. Independent marker methods that saw
        # the organism outrank the ratio: an organism two or more methods
        # agree on is only rejected on the absolute floor with patchy
        # coverage, never on the ratio alone (Blautia producta, seen by four
        # methods at 0.07%, was wrongly rejected before this).
        single_method = len(getattr(o, "methods", ()) or ()) <= 1
        # Calibrated on the simulated communities (scripts/benchmark_lanes.py):
        # every false single-method call that had a relative in its competition
        # mapped at an identity 0.5-2.5% below that relative's and recruited a
        # fraction of its reads; every true minor species mapped within 0.3%
        # of its relative whatever the read ratio. The identity gap is the
        # discriminator; the read ratio and evenness no longer gate it.
        ratio_shadow = (beaten_by is not None and identity < beaten_by[2] - IDENTITY_GAP
                        and u < SHADOW_FRACTION * beaten_by[1])
        absolute_shadow = beaten_by is not None and u < MIN_UNIQUE_FRAGMENTS * 5 and (single_method or evenness < MIN_EVENNESS)
        relative = ""
        if beaten_by is not None and (absolute_shadow or (single_method and ratio_shadow)):
            status = "not_detected"
            relative = beaten_by[0]
            reason = (f"{u:,} unique fragments at {identity:.1%} identity, while {beaten_by[0]} took {beaten_by[1]:,} at "
                      f"{beaten_by[2]:.1%} in the same competition; the reads belong to the relative")
        elif u >= MIN_UNIQUE_FRAGMENTS and breadth >= MIN_BREADTH and identity >= MIN_IDENTITY and evenness >= MIN_EVENNESS:
            status, reason = "supported", (
                f"{u:,} uniquely mapping fragments over {breadth:.1%} of the genome at {identity:.1%} identity, "
                f"spread evenly (evenness {evenness:.2f}), with {len(competitors)} relatives competing")
        elif u >= PROVISIONAL_FRAGMENTS:
            short = []
            if u < MIN_UNIQUE_FRAGMENTS:
                short.append(f"only {u} unique fragments")
            if breadth < MIN_BREADTH:
                short.append(f"breadth {breadth:.1%}")
            if identity < MIN_IDENTITY:
                short.append(f"identity {identity:.1%}")
            if evenness < MIN_EVENNESS:
                short.append(f"uneven coverage ({evenness:.2f})")
            status, reason = "provisional", "; ".join(short) or "below the supported operating point"
        else:
            best_alt = max(competitors.items(), key=lambda kv: kv[1], default=None)
            status = "not_detected"
            relative = best_alt[0] if (best_alt and best_alt[1] > u) else ""
            reason = (f"{u} unique fragments; reads in this genus map to {best_alt[0]} ({best_alt[1]:,} fragments)"
                      if best_alt and best_alt[1] > u else f"{u} unique fragments after competition")
        verdicts.append(Verdict(
            organism=o.display, genome_id=gid, status=status, reason=reason, unique_fragments=u,
            total_fragments=total_by_genome.get(gid, 0), breadth=breadth, mean_depth=mean_depth, identity=identity,
            evenness=evenness, genome_length=glen, competitors=competitors, previous_status=o.status,
            relative=relative,
        ))
    # Clean the big intermediates; keep the BAM index-less summary.
    for p in (ref, bam, bam.with_suffix(".bam.bai")):
        p.unlink(missing_ok=True)
    for p in work_dir.glob(f"{sample}.confirm.*.bt2*"):
        p.unlink(missing_ok=True)
    result = {
        "status": "completed", "n_candidates": len(candidates), "n_tested": len(tested), "n_pending": len(pending),
        "n_genomes_in_competition": len(lengths), "verdicts": [v.to_json() for v in verdicts],
        # every genome in the competition, candidate or relative: what mapped
        # uniquely to it, so the inventory can apportion a genus total among
        # the populations the mapping told apart
        "genomes": {gid: {"organism": genome_of[gid], "role": roles.get(gid, ""), "unique_fragments": int(uniq[gid]),
                          "genome_length": int(lengths[gid]), "identity": round(identity_of.get(gid, 0.0), 4)}
                    for gid in lengths},
        "cache_inputs": cache_inputs,
        "operating_points": _points(), "elapsed_s": round(time.monotonic() - t0, 1),
        "meaning": (
            "Each candidate genome competed with its detected relatives for the same reads. supported: many "
            "independent fragments across the genome at species-level identity; provisional: some reads but short "
            "of that; not_detected: the reads that reached this genus belong to a relative. pending: the reference "
            "genome is not yet on disk, so nothing was decided."
        ),
    }
    (work_dir / f"{sample}.confirmation.json").write_text(json.dumps(result, indent=1))
    cache_file.write_text(json.dumps(result, indent=1))
    for old in work_dir.glob(f"{sample}.confirmation.*.json"):
        if old != cache_file:
            old.unlink(missing_ok=True)  # one competition per sample is kept
    return result


def _points() -> dict[str, Any]:
    return {"min_unique_fragments": MIN_UNIQUE_FRAGMENTS, "min_breadth": MIN_BREADTH, "min_identity": MIN_IDENTITY,
            "min_evenness": MIN_EVENNESS, "provisional_fragments": PROVISIONAL_FRAGMENTS, "mapq_unique": MAPQ_UNIQUE,
            "shadow_fraction": SHADOW_FRACTION, "identity_gap": IDENTITY_GAP,
            "max_candidates": MAX_CANDIDATES, "calibration": "initial operating points; calibrated by scripts/benchmark_lanes.py"}


_RELATIVE_PATTERNS: Final = (
    re.compile(r"^[\d,]+ unique fragments at [\d.]+% identity, while (?P<name>.+?) took [\d,]+ at "),
    re.compile(r"^\d+ unique fragments; reads in this genus map to (?P<name>.+?) \([\d,]+ fragments\)"),
)


def relative_of(verdict: dict[str, Any]) -> str:
    """The relative a rejected call's reads belong to.

    Written on the verdict since the field exists; read back out of the
    verdict's own reason line for the confirmations cached before it did,
    so a cached competition is not run again for a name it already gave.
    """
    if verdict.get("status") != "not_detected":
        return ""
    if verdict.get("relative"):
        return str(verdict["relative"])
    reason = str(verdict.get("reason") or "")
    for pattern in _RELATIVE_PATTERNS:
        m = pattern.match(reason)
        if m:
            return m.group("name")
    return ""


def apply_verdicts(inventory_json: dict[str, Any], confirmation: dict[str, Any]) -> dict[str, Any]:
    """Write confirmation outcomes back onto the inventory records."""
    by_name = {v["organism"]: v for v in confirmation.get("verdicts") or []}
    n_changed = 0
    kept: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = list(inventory_json.get("rejected") or [])
    from openbiota import inventory as _inventory

    for rec in inventory_json.get("organisms") or []:
        display = _inventory.display_of(rec)
        v = by_name.get(display)
        if v is None or v["status"] == "pending":
            if v is not None:
                rec["confidence_basis"] = (rec.get("confidence_basis") or "") + "; confirmation pending: " + v["reason"]
            kept.append(rec)
            continue
        metrics = {k: v[k] for k in ("unique_fragments", "breadth", "identity", "evenness", "competitors")}
        if v["status"] == "not_detected":
            # A rejected call leaves the organism list and is reported as a
            # rejected call with its reason - never silently dropped.
            n_changed += 1
            rejected.append({**rec, "status": "rejected", "confidence_basis": f"competitive confirmation: {v['reason']}",
                             "confirmation": metrics, "reads_belong_to": relative_of(v) or None})
            continue
        if rec.get("status") != v["status"]:
            n_changed += 1
        rec["status"] = v["status"]
        rec["confidence_basis"] = f"competitive confirmation: {v['reason']}"
        rec["confirmation"] = metrics
        kept.append(rec)
    inventory_json["organisms"] = kept
    inventory_json["rejected"] = rejected
    inventory_json["n_confirmation_changed"] = n_changed
    # Every derived block is recomputed from the records as they now stand:
    # counts, categories, incremental gain, the comparison lists. Patching a
    # few keys by hand left `expansion_supported` and friends stale.
    from openbiota import inventory as _inventory

    # One share of the whole per organism, now that the competition has
    # told the populations within each genus apart (spec 0.8.4 §6).
    _inventory.unify_shares(inventory_json, confirmation)
    rebuilt = _inventory.from_json(inventory_json)
    if rebuilt is not None:
        fresh = rebuilt.to_json()
        for key in ("n_organisms", "n_named", "n_unnamed", "n_genera", "n_rankable", "n_strain_resolved",
                    "n_in_composition", "n_secondary_only", "composition_total_percent", "unplaced_percent",
                    "share_unification", "counts", "comparison"):
            inventory_json[key] = fresh[key]
        inventory_json["counts"]["rejected"] = len(rejected)
    else:
        inventory_json["n_organisms"] = len(kept)
    comp = inventory_json.get("comparison") or {}
    comp["rejected"] = [{"organism": _inventory.display_of(r),
                         "reason": r.get("confidence_basis", ""), "detected_by": r.get("detected_by", [])} for r in rejected]
    inventory_json["comparison"] = comp
    return inventory_json
