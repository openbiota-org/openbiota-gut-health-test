"""The whole-genome fungal lane: candidates, competitive confirmation and the
fragment ledger (spec §6.2-6.3, §7.1).

Two stages.

**Candidates.** The sample's read sketch (already made for the bacterial
GTDB lane, same ``c``) is profiled against a sketch database of every
admitted fungal assembly. Anything with a containment ANI above the
candidate floor is a candidate, as is every fungus the marker lane
(EukDetect2) or the pathogen screen named.

**Confirmation.** Every read pair is aligned competitively - not a
fungi-enriched subset - against the candidate genomes, their congeneric
near neighbours, a fixed core panel of common gut and food fungi, and
decoys (gut bacteria, PhiX). Each pair is adjudicated once: identity from
``NM`` against the aligned length, uniqueness from MAPQ and from competing
alignments within a small score margin. A pair whose competitors are all
one species is species evidence; competitors across species of one genus
make genus-level fungal evidence; a competitor that is a decoy makes the
pair *cross-kingdom ambiguous*, which is reported and never counted as
fungal. Support for a taxon requires fragments spread across separated
genomic regions and more than one contig, so a contaminated contig or a
conserved rDNA locus cannot manufacture a species.

The starting support rule (20 fragments, 3 separated regions) is a locked
starting value from the spec's calibration grid, recorded as such in the
output; it is not a validated clinical threshold.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.mycobiome.references import Assembly, Lock

CACHE_VERSION: Final = 5
CANDIDATE_ANI: Final = 93.0          # sylph Adjusted_ANI floor to become a candidate
SPECIES_ANI: Final = 95.0            # containment ANI read as same-species evidence
MIN_KMERS: Final = 20
IDENTITY_MIN: Final = 0.97
SPECIES_IDENTITY_MIN: Final = 0.985   # mean identity below this: the reads are a near neighbour, call the genus
ALIGNED_FRACTION_MIN: Final = 0.80
MAPQ_UNIQUE: Final = 20
COMPETITION_MARGIN: Final = 6        # alignment-score points; within this, hits compete
REGION_BP: Final = 100_000
SUPPORT_FRAGMENTS: Final = 20
SUPPORT_REGIONS: Final = 3
PROVISIONAL_FRAGMENTS: Final = 5
PROVISIONAL_REGIONS: Final = 2
TRACE_FRAGMENTS: Final = 2            # listed as a trace, never interpreted
CONCENTRATION_MAX: Final = 0.70      # >70% of fragments in one 20-kb window -> locus-concentrated
MAX_NEIGHBOURS_PER_GENUS: Final = 6
MAX_ASSEMBLIES_PER_CANDIDATE: Final = 1   # one backbone per species: conspecific duplicates collapse MAPQ
MIN_CONTIG_BP: Final = 5_000

SUPPORT_RULE: Final = {
    "rule_id": "myco-support-v1",
    "supported": {"fragments": SUPPORT_FRAGMENTS, "regions": SUPPORT_REGIONS, "contigs": 2},
    "provisional": {"fragments": PROVISIONAL_FRAGMENTS, "regions": PROVISIONAL_REGIONS},
    "trace": {"fragments": TRACE_FRAGMENTS},
    "identity_min": IDENTITY_MIN, "species_identity_min": SPECIES_IDENTITY_MIN,
    "aligned_fraction_min": ALIGNED_FRACTION_MIN, "mapq_unique": MAPQ_UNIQUE,
    "competition_margin": COMPETITION_MARGIN, "region_bp": REGION_BP,
    "status": "starting_values_from_spec_calibration_grid_not_yet_qualified_on_mocks",
}

#: Common gut and food fungi always present in the competitive index, so a
#: candidate is judged against the organisms most likely to compete with it.
CORE_PANEL_SPECIES: Final = (
    "Saccharomyces cerevisiae", "Candida albicans", "Candida tropicalis", "Candida parapsilosis",
    "Nakaseomyces glabratus", "Debaryomyces hansenii", "Malassezia restricta", "Malassezia globosa",
    "Clavispora lusitaniae", "Pichia kudriavzevii", "Kluyveromyces marxianus", "Kluyveromyces lactis",
    "Geotrichum candidum", "Cyberlindnera jadinii", "Torulaspora delbrueckii", "Meyerozyma guilliermondii",
    "Agaricus bisporus", "Lentinula edodes", "Pleurotus ostreatus", "Aspergillus fumigatus", "Aspergillus niger",
    "Aspergillus oryzae", "Penicillium chrysogenum", "Penicillium camemberti", "Penicillium roqueforti",
    "Rhodotorula mucilaginosa", "Candidozyma auris", "Cladosporium sphaerospermum", "Fusarium graminearum",
    "Wickerhamomyces anomalus",
)


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class SketchHit:
    accession: str
    species_taxid: int | None
    species: str
    ani: float
    eff_cov: float
    sequence_abundance: float
    kmers: int

    def to_json(self) -> dict[str, Any]:
        return {"accession": self.accession, "species_taxid": self.species_taxid, "species": self.species,
                "adjusted_ani": self.ani, "effective_coverage": self.eff_cov,
                "sequence_abundance_within_fungal_db": self.sequence_abundance, "kmers": self.kmers}


@dataclass(slots=True)
class TaxonSupport:
    """Confirmed alignment support for one species (or one genus, when reads
    could not be placed below it)."""

    key: str                              # "sp:<taxid>" or "genus:<name>"
    rank: str                             # species | genus
    name: str
    species_taxid: int | None
    genus: str
    accessions: list[str]
    fragments: int = 0
    fragments_unique: int = 0
    fragments_short_contigs: int = 0
    discovery_fragments: int | None = None
    regions: int = 0
    contigs: int = 0
    top_window_fraction: float = 0.0
    identity_mean: float = 0.0
    mapq_mean: float = 0.0
    genome_bases: int = 0                 # backbone length used for coverage normalisation
    detection_state: str = "not_supported"
    reasons: list[str] = field(default_factory=list)
    backbone: str | None = None
    marker_lane: bool = False
    sketch_ani: float | None = None

    @property
    def expected_depth(self) -> float | None:
        return None

    def to_json(self) -> dict[str, Any]:
        return {
            "key": self.key, "rank": self.rank, "name": self.name, "species_taxid": self.species_taxid,
            "genus": self.genus, "accessions": self.accessions, "fragments": self.fragments,
            "fragments_unique_mapq": self.fragments_unique, "fragments_on_short_contigs": self.fragments_short_contigs,
            "discovery_fragments": self.discovery_fragments, "separated_regions_100kb": self.regions,
            "contigs_hit": self.contigs, "top_20kb_window_fraction": round(self.top_window_fraction, 4),
            "identity_mean": round(self.identity_mean, 5), "mapq_mean": round(self.mapq_mean, 2),
            "backbone_accession": self.backbone, "backbone_bases": self.genome_bases,
            "detection_state": self.detection_state, "reasons": self.reasons,
            "also_in_marker_lane": self.marker_lane, "sketch_adjusted_ani": self.sketch_ani,
        }


@dataclass(slots=True)
class Ledger:
    """Fragment accounting for the primary denominator (spec §7.1)."""

    eligible_fragments: int
    aligned_fragments: int = 0
    fungal_confident: int = 0             # F: species- or genus-level fungal, unique or fungal-only competition
    fungal_low_identity: int = 0
    cross_kingdom_ambiguous: int = 0
    decoy_fragments: int = 0
    by_taxon_fragments: dict[str, int] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "denominator_stage": "qc_nonhost_before_bacterial_filter",
            "eligible_fragments": self.eligible_fragments,
            "aligned_to_index_fragments": self.aligned_fragments,
            "fungal_supported_fragments": self.fungal_confident,
            "fungal_below_identity_fragments": self.fungal_low_identity,
            "cross_kingdom_ambiguous_fragments": self.cross_kingdom_ambiguous,
            "decoy_fragments": self.decoy_fragments,
            "unassigned_fragments": max(0, self.eligible_fragments - self.aligned_fragments),
        }


@dataclass(slots=True)
class GenomeLaneResult:
    status: str
    sketch_hits: list[SketchHit]
    taxa: list[TaxonSupport]
    ledger: Ledger | None
    index_accessions: list[str]
    decoy_accessions: list[str]
    bam: str | None
    cached: bool
    error: str | None = None
    timings_s: dict[str, float] = field(default_factory=dict)
    discovered: dict[str, int] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "status": self.status, "cached": self.cached, "error": self.error, "support_rule": SUPPORT_RULE,
            "sketch": {"candidate_ani_floor": CANDIDATE_ANI, "min_kmers": MIN_KMERS, "hits": [h.to_json() for h in self.sketch_hits]},
            "competitive_index": {"n_fungal_assemblies": len(self.index_accessions), "assemblies": self.index_accessions,
                                  "decoys": self.decoy_accessions},
            "discovery": {"min_fragments": DISCOVERY_MIN_FRAGMENTS, "species_fragments": self.discovered},
            "ledger": None if self.ledger is None else self.ledger.to_json(),
            "taxa": [t.to_json() for t in self.taxa], "bam": self.bam, "timings_s": self.timings_s,
        }


# --------------------------------------------------------------------------- #
# tools
# --------------------------------------------------------------------------- #


def _sylph() -> Path | None:
    p = Path(__file__).resolve().parents[2] / "vendor" / "sylph" / "bin" / "sylph"
    return p if p.is_file() else (Path(shutil.which("sylph")) if shutil.which("sylph") else None)


def _tool(name: str) -> str | None:
    return shutil.which(name)


def available(lock: Lock | None, base: Path) -> bool:
    return (lock is not None and bool(lock.assemblies) and _sylph() is not None and _tool("minimap2") is not None
            and _tool("samtools") is not None and lock.sketch_db is not None and (base / lock.sketch_db).is_file())


def _fp(*parts: str) -> str:
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def _acc_of(genome_file: str) -> str:
    m = re.search(r"(GC[AF]_\d+\.\d+)", genome_file)
    return m.group(1) if m else Path(genome_file).name


# --------------------------------------------------------------------------- #
# stage 1: candidates from the sketch database
# --------------------------------------------------------------------------- #


def sketch_candidates(*, sample: str, r1: Path, r2: Path | None, lock: Lock, base: Path, work_dir: Path,
                      read_sketch: Path | None, threads: int) -> list[SketchHit]:
    sylph = _sylph()
    assert sylph is not None and lock.sketch_db is not None
    work_dir.mkdir(parents=True, exist_ok=True)
    sketch = read_sketch
    if sketch is None or not sketch.is_file():
        sk_dir = work_dir / "sketches"
        sk_dir.mkdir(exist_ok=True)
        cmd = [str(sylph), "sketch", "-t", str(threads), "-c", "200", "-d", str(sk_dir)]
        cmd += ["-1", str(r1), "-2", str(r2)] if r2 is not None else ["-r", str(r1)]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        found = sorted(sk_dir.glob("*.sylsp"))
        if not found:
            raise RuntimeError("sylph produced no read sketch")
        sketch = found[0]
    out = work_dir / f"{sample}.fungi.sylph.tsv"
    proc = subprocess.run(
        [str(sylph), "profile", "-t", str(threads), "-M", str(MIN_KMERS), "-u", "--min-spacing", "10",
         str(base / lock.sketch_db), str(sketch), "-o", str(out)],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"sylph profile failed: {proc.stderr[-500:]}")
    by_acc = lock.by_accession()
    hits: list[SketchHit] = []
    if out.is_file():
        with out.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = _acc_of(row.get("Genome_file", ""))
                a = by_acc.get(acc)
                ani = float(row.get("Adjusted_ANI") or row.get("Naive_ANI") or 0.0)
                hits.append(SketchHit(
                    accession=acc, species_taxid=a.species_taxid if a else None, species=a.species if a else acc,
                    ani=ani, eff_cov=float(row.get("Eff_cov") or 0.0),
                    sequence_abundance=float(row.get("Sequence_abundance") or 0.0),
                    kmers=int(float(row.get("Containment_ind", "0/0").split("/")[0] or 0)) if "/" in (row.get("Containment_ind") or "") else 0,
                ))
    hits.sort(key=lambda h: -h.ani)
    return hits


# --------------------------------------------------------------------------- #
# stage 1b: discovery by alignment to one representative per species
# --------------------------------------------------------------------------- #

DISCOVERY_MIN_FRAGMENTS: Final = 2


def discover_candidates(*, sample: str, r1: Path, r2: Path | None, lock: Lock, base: Path, work_dir: Path,
                        threads: int, log=print) -> dict[int, int]:
    """{species_taxid: confident fragments} from aligning every pair to the
    species-representative index. Sketching cannot see a fungus at 0.003%
    of fragments; alignment can. Cached per sample."""
    import pysam

    mmi = base / "species_index" / "fungi_species.mmi"
    meta_p = base / "species_index" / "fungi_species.json"
    if not (mmi.is_file() and meta_p.is_file()):
        return {}
    meta = json.loads(meta_p.read_text())
    contig_acc: dict[str, str] = meta["contig_acc"]
    bam = work_dir / f"{sample}.fungal_discovery.{lock.lock_id}.bam"
    if not bam.is_file() or bam.stat().st_size == 0:
        align(fasta=mmi, r1=r1, r2=r2, out_bam=bam, threads=threads, log=log)
    by_acc = lock.by_accession()
    decoys = {d["accession"] for d in lock.decoys}
    counts: dict[int, int] = defaultdict(int)
    seen: set[str] = set()
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for aln in fh:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_read2:
                continue
            if aln.query_name in seen:
                continue
            seen.add(aln.query_name)
            acc = contig_acc.get(aln.reference_name, aln.reference_name.split("|")[0])
            if acc in decoys:
                continue
            ident, aligned = _identity(aln)
            if ident < IDENTITY_MIN or aligned < ALIGNED_FRACTION_MIN * (aln.query_length or 1):
                continue
            a = by_acc.get(acc)
            if a is not None and a.species_taxid:
                counts[a.species_taxid] += 1
    return {tid: n for tid, n in counts.items() if n >= DISCOVERY_MIN_FRAGMENTS}


# --------------------------------------------------------------------------- #
# stage 2: competitive index and alignment
# --------------------------------------------------------------------------- #


def _pick(assemblies: list[Assembly], n: int) -> list[Assembly]:
    """Reference genomes first, then the longest scaffold-level ones."""
    order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    ranked = sorted(assemblies, key=lambda a: (0 if a.source in ("refseq_reference", "both") else 1,
                                               order.get(a.level, 4), -a.bases))
    return ranked[:n]


def choose_index(lock: Lock, candidate_species: set[int], candidate_accessions: set[str]) -> list[Assembly]:
    """Candidates, their congeneric neighbours and the core panel."""
    by_species = lock.species_assemblies()
    by_acc = lock.by_accession()
    chosen: dict[str, Assembly] = {}
    genera: set[str] = set()
    name_to_taxid = {a.species: a.species_taxid for a in lock.assemblies if a.species_taxid}
    for name in CORE_PANEL_SPECIES:
        tid = name_to_taxid.get(name)
        if tid:
            candidate_species.add(tid)
    for tid in candidate_species:
        for a in _pick(by_species.get(tid, []), MAX_ASSEMBLIES_PER_CANDIDATE):
            chosen[a.accession] = a
            genera.add(a.genus)
    for acc in candidate_accessions:
        a = by_acc.get(acc)
        if a is not None and a.is_fungus:
            chosen[a.accession] = a
            genera.add(a.genus)
    # Near neighbours: other species of the same genera, one assembly each.
    for genus in sorted(g for g in genera if g):
        others = [tid for tid, asms in by_species.items() if asms and asms[0].genus == genus and tid not in candidate_species]
        # Prefer species with a reference assembly, then more assemblies (better known).
        others.sort(key=lambda tid: (0 if any(a.source in ("refseq_reference", "both") for a in by_species[tid]) else 1,
                                     -len(by_species[tid])))
        for tid in others[:MAX_NEIGHBOURS_PER_GENUS]:
            a = _pick(by_species[tid], 1)[0]
            chosen[a.accession] = a
    return sorted(chosen.values(), key=lambda a: a.accession)


def index_key(lock: Lock, assemblies: list[Assembly]) -> str:
    return _fp(str(CACHE_VERSION), lock.lock_id, *sorted(a.accession for a in assemblies), *sorted(d["accession"] for d in lock.decoys))


def build_index(lock: Lock, base: Path, assemblies: list[Assembly], cache_dir: Path, log=print) -> tuple[Path, dict[str, str], dict[str, int]]:
    """Concatenate the chosen fungal assemblies and the decoys into one FASTA
    with ``<ACC>|<contig>`` headers. Returns (fasta, contig->acc, contig->len).
    Cached by the sorted accession set."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    decoys = lock.decoys
    key = index_key(lock, assemblies)
    fasta = cache_dir / f"index.{key}.fa"
    meta = cache_dir / f"index.{key}.json"
    if fasta.is_file() and meta.is_file():
        m = json.loads(meta.read_text())
        return fasta, m["contig_acc"], {k: int(v) for k, v in m["contig_len"].items()}
    log(f"  building competitive index: {len(assemblies)} fungal assemblies + {len(decoys)} decoys")
    contig_acc: dict[str, str] = {}
    contig_len: dict[str, int] = {}

    def _copy(src: Path, acc: str, out) -> None:
        import gzip
        opener = gzip.open if src.suffix == ".gz" else open
        name = None
        n = 0
        with opener(src, "rt") as fh:  # type: ignore[operator]
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        contig_len[name] = n
                    raw = line[1:].split()[0]
                    name = f"{acc}|{raw}"
                    contig_acc[name] = acc
                    n = 0
                    out.write(f">{name}\n")
                else:
                    n += len(line.strip())
                    out.write(line)
            if name is not None:
                contig_len[name] = n

    with fasta.open("w") as out:
        for a in assemblies:
            _copy(base / a.path, a.accession, out)
        for d in decoys:
            p = (base / d["path"]).resolve()
            _copy(p, d["accession"], out)
    meta.write_text(json.dumps({"contig_acc": contig_acc, "contig_len": contig_len}))
    return fasta, contig_acc, contig_len


#: The species-representative index is ~90 GB in memory. One alignment
#: against it at a time per machine; a second one swaps the host.
_DISCOVERY_LOCK: Final = Path("/tmp/openbiota-mycobiome-discovery.lock")


def align(*, fasta: Path, r1: Path, r2: Path | None, out_bam: Path, threads: int, log=print) -> None:
    """minimap2 short-read mode over all pairs, secondaries kept for
    competition; only mapped records are retained. Alignments against the
    big species index are serialised machine-wide with a file lock."""
    import fcntl

    mm = _tool("minimap2")
    st = _tool("samtools")
    assert mm and st
    big = fasta.suffix == ".mmi"
    lock_fh = None
    if big:
        lock_fh = _DISCOVERY_LOCK.open("w")
        log("  waiting for the discovery-index lock (one 90 GB alignment at a time)")
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
    try:
        _align(mm, st, fasta=fasta, r1=r1, r2=r2, out_bam=out_bam, threads=threads, log=log)
    finally:
        if lock_fh is not None:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)
            lock_fh.close()


def _align(mm: str, st: str, *, fasta: Path, r1: Path, r2: Path | None, out_bam: Path, threads: int, log=print) -> None:
    cmd_mm = [mm, "-ax", "sr", "-t", str(threads), "--secondary=yes", "-N", "5", str(fasta), str(r1)]
    if r2 is not None:
        cmd_mm.append(str(r2))
    log(f"  aligning all read pairs competitively ({out_bam.name})")
    with out_bam.with_suffix(".log").open("w") as errlog:
        p1 = subprocess.Popen(cmd_mm, stdout=subprocess.PIPE, stderr=errlog)
        p2 = subprocess.Popen([st, "view", "-@", "4", "-b", "-F", "4", "-o", str(out_bam), "-"], stdin=p1.stdout, stderr=errlog)
        p1.stdout.close()  # type: ignore[union-attr]
        p2.communicate()
        p1.wait()
    if p1.returncode != 0 or p2.returncode != 0:
        raise RuntimeError("minimap2/samtools alignment failed; see " + str(out_bam.with_suffix(".log")))


# --------------------------------------------------------------------------- #
# stage 3: the fragment ledger
# --------------------------------------------------------------------------- #


def _identity(aln) -> tuple[float, int]:
    """(identity, aligned query length) from NM and the CIGAR."""
    stats = aln.get_cigar_stats()[0]
    matches = stats[0]           # M (may include mismatches)
    ins = stats[1]
    dele = stats[2]
    aligned = matches + ins
    nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
    denom = matches + ins + dele
    if denom <= 0:
        return 0.0, 0
    return max(0.0, 1.0 - nm / denom), aligned


def adjudicate(bam: Path, contig_acc: dict[str, str], contig_len: dict[str, int], lock: Lock, eligible: int,
               log=print) -> tuple[Ledger, dict[str, TaxonSupport]]:
    log(f"  adjudicating fragments in {bam.name}")
    import pysam

    by_acc = lock.by_accession()
    decoy_accs = {d["accession"] for d in lock.decoys}
    acc_species: dict[str, tuple[int | None, str, str]] = {}
    for acc in set(contig_acc.values()):
        a = by_acc.get(acc)
        if a is not None:
            acc_species[acc] = (a.species_taxid, a.species, a.genus)

    # Group alignments by fragment. Only mapped records are in the BAM.
    frags: dict[str, dict[int, list[Any]]] = defaultdict(lambda: {1: [], 2: []})
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for aln in fh:
            if aln.is_unmapped or aln.is_supplementary:
                continue
            mate = 2 if aln.is_read2 else 1
            frags[aln.query_name][mate].append((
                aln.reference_name, aln.reference_start, aln.mapping_quality, aln.is_secondary,
                aln.get_tag("AS") if aln.has_tag("AS") else 0, *_identity(aln), aln.query_length or 0,
            ))
    ledger = Ledger(eligible_fragments=eligible, aligned_fragments=len(frags))
    support: dict[str, TaxonSupport] = {}
    windows: dict[str, set[tuple[str, int]]] = defaultdict(set)
    contigs: dict[str, set[str]] = defaultdict(set)
    fine: dict[str, dict[tuple[str, int], int]] = defaultdict(lambda: defaultdict(int))
    ident_sum: dict[str, float] = defaultdict(float)
    mapq_sum: dict[str, float] = defaultdict(float)

    def slot(key: str, rank: str, name: str, tid: int | None, genus: str, acc: str) -> TaxonSupport:
        if key not in support:
            support[key] = TaxonSupport(key=key, rank=rank, name=name, species_taxid=tid, genus=genus, accessions=[])
        s = support[key]
        if acc not in s.accessions and acc in acc_species:
            s.accessions.append(acc)
        return s

    for _qname, mates in frags.items():
        recs = mates[1] or mates[2]
        if not recs:
            continue
        primary = [r for r in recs if not r[3]]
        if not primary:
            continue
        p = max(primary, key=lambda r: r[4])
        ref, pos, mapq, _sec, as_, ident, aligned, qlen = p
        acc = contig_acc.get(ref, ref.split("|")[0])
        competitors = {contig_acc.get(r[0], r[0].split("|")[0]) for r in recs if r[3] and r[4] >= as_ - COMPETITION_MARGIN}
        competitors.discard(acc)
        if acc in decoy_accs:
            ledger.decoy_fragments += 1
            continue
        if competitors & decoy_accs:
            ledger.cross_kingdom_ambiguous += 1
            continue
        if ident < IDENTITY_MIN or (qlen and aligned < ALIGNED_FRACTION_MIN * qlen):
            ledger.fungal_low_identity += 1
            continue
        short = contig_len.get(ref, 0) < MIN_CONTIG_BP
        tid, species, genus = acc_species.get(acc, (None, acc, ""))
        comp_species = {acc_species.get(c, (None, c, ""))[0] for c in competitors}
        comp_species.discard(tid)
        comp_genera = {acc_species.get(c, (None, c, ""))[2] for c in competitors}
        if comp_species and len(comp_species - {None}) > 0:
            # Competing species: resolve at the lowest supported common rank.
            if comp_genera <= {genus} and genus:
                key, rank, name = f"genus:{genus}", "genus", genus
                s = slot(key, rank, name, None, genus, acc)
            else:
                ledger.fungal_confident += 1
                s = slot("fungi:unresolved", "kingdom", "Fungi (placed above genus)", None, "", acc)
                s.fragments += 1
                if short:
                    s.fragments_short_contigs += 1
                continue
        else:
            key, rank, name = f"sp:{tid}" if tid else f"acc:{acc}", "species", species
            s = slot(key, rank, name, tid, genus, acc)
        ledger.fungal_confident += 1
        s.fragments += 1
        if short:
            s.fragments_short_contigs += 1
        if mapq >= MAPQ_UNIQUE:
            s.fragments_unique += 1
        windows[s.key].add((ref, pos // REGION_BP))
        contigs[s.key].add(ref)
        fine[s.key][(ref, pos // 20_000)] += 1
        ident_sum[s.key] += ident
        mapq_sum[s.key] += mapq

    for key, s in support.items():
        s.regions = len(windows[key])
        s.contigs = len(contigs[key])
        if s.fragments:
            s.identity_mean = ident_sum[key] / s.fragments
            s.mapq_mean = mapq_sum[key] / s.fragments
            s.top_window_fraction = max(fine[key].values()) / s.fragments if fine[key] else 0.0
        # Backbone: the reference-preferred assembly among those hit.
        cands = [by_acc[a] for a in s.accessions if a in by_acc]
        if cands:
            bb = _pick(cands, 1)[0]
            s.backbone, s.genome_bases = bb.accession, bb.bases
        single_contig = any(by_acc[a].n_contigs <= 1 for a in s.accessions if a in by_acc)
        reasons: list[str] = []
        if s.top_window_fraction > CONCENTRATION_MAX and s.fragments >= PROVISIONAL_FRAGMENTS:
            reasons.append("locus_concentrated_possible_rdna_mito_or_contaminant")
        if s.fragments and s.fragments_short_contigs / s.fragments > 0.5 and s.fragments >= PROVISIONAL_FRAGMENTS:
            reasons.append("mostly_on_short_contigs_repeat_or_organelle_likely")
        if s.fragments >= SUPPORT_FRAGMENTS and s.regions >= SUPPORT_REGIONS and (s.contigs >= 2 or single_contig) and not reasons:
            s.detection_state = "supported"
        elif s.fragments >= PROVISIONAL_FRAGMENTS and s.regions >= PROVISIONAL_REGIONS:
            s.detection_state = "provisional"
            if s.fragments < SUPPORT_FRAGMENTS:
                reasons.append("below_support_fragment_threshold")
            if s.regions < SUPPORT_REGIONS:
                reasons.append("too_few_separated_regions")
        elif s.fragments >= TRACE_FRAGMENTS:
            s.detection_state = "trace"
            reasons.append("trace_below_provisional_rule")
        else:
            s.detection_state = "not_supported"
            reasons.append("insufficient_fragments_or_regions")
        if s.rank == "species" and s.fragments >= PROVISIONAL_FRAGMENTS and s.identity_mean < SPECIES_IDENTITY_MIN and s.genus:
            # A consistent ~1.5%+ divergence from the nearest genome is a
            # relative of that species, not the species: report the genus.
            reasons.append(f"identity_{s.identity_mean:.3f}_below_species_confidence_nearest_{s.name.replace(' ', '_')}")
            s.rank = "genus"
            s.name = s.genus
            s.species_taxid = None
        if s.rank == "genus":
            s.detection_state = "ambiguous_complex" if s.detection_state == "supported" else s.detection_state
            if "reads_compete_across_species_of_one_genus" not in reasons and not any(r.startswith("identity_") for r in reasons):
                reasons.append("reads_compete_across_species_of_one_genus")
        s.reasons = reasons
        ledger.by_taxon_fragments[key] = s.fragments
    return ledger, support


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #


def run_genome_lane(*, sample: str, r1: Path, r2: Path | None, lock: Lock, base: Path, work_dir: Path,
                    eligible_fragments: int, extra_species_taxids: set[int], extra_accessions: set[str],
                    read_sketch: Path | None, threads: int = 16, log=print) -> GenomeLaneResult:
    import time

    work_dir.mkdir(parents=True, exist_ok=True)
    key = _fp(str(CACHE_VERSION), lock.lock_id, sample, str(eligible_fragments),
              *sorted(str(t) for t in extra_species_taxids), *sorted(extra_accessions),
              f"{r1.stat().st_size}", f"{r2.stat().st_size if r2 else 0}")
    cached = work_dir / f"genome_lane.{key}.json"
    if cached.is_file():
        d = json.loads(cached.read_text())
        return _from_json(d, cached=True)
    t: dict[str, float] = {}
    t0 = time.time()
    try:
        hits = sketch_candidates(sample=sample, r1=r1, r2=r2, lock=lock, base=base, work_dir=work_dir,
                                 read_sketch=read_sketch, threads=threads)
        t["sketch_profile"] = round(time.time() - t0, 1)
        td = time.time()
        discovered = discover_candidates(sample=sample, r1=r1, r2=r2, lock=lock, base=base, work_dir=work_dir,
                                         threads=threads, log=log)
        t["discovery_align"] = round(time.time() - td, 1)
        cand_species = {h.species_taxid for h in hits if h.ani >= CANDIDATE_ANI and h.species_taxid}
        cand_accs = {h.accession for h in hits if h.ani >= CANDIDATE_ANI}
        cand_species |= set(discovered)
        cand_species |= extra_species_taxids
        cand_accs |= extra_accessions
        assemblies = choose_index(lock, set(cand_species), set(cand_accs))
        t1 = time.time()
        fasta, contig_acc, contig_len = build_index(lock, base, assemblies, base / "index_cache", log=log)
        t["index"] = round(time.time() - t1, 1)
        t2 = time.time()
        bam = work_dir / f"{sample}.fungal_competitive.{index_key(lock, assemblies)}.bam"
        if not bam.is_file() or bam.stat().st_size == 0:
            align(fasta=fasta, r1=r1, r2=r2, out_bam=bam, threads=threads, log=log)
        t["align"] = round(time.time() - t2, 1)
        t3 = time.time()
        ledger, support = adjudicate(bam, contig_acc, contig_len, lock, eligible_fragments, log=log)
        t["adjudicate"] = round(time.time() - t3, 1)
        marker_species = extra_species_taxids
        for s in support.values():
            s.marker_lane = s.species_taxid in marker_species if s.species_taxid else False
            s.discovery_fragments = discovered.get(s.species_taxid) if s.species_taxid else None
            s.sketch_ani = max((h.ani for h in hits if h.species_taxid == s.species_taxid), default=None) if s.species_taxid else None
        taxa = sorted(support.values(), key=lambda s: -s.fragments)
        for s in taxa:
            s.reasons = list(s.reasons)
        res = GenomeLaneResult(status="resolved", sketch_hits=hits, taxa=taxa, ledger=ledger,
                               index_accessions=[a.accession for a in assemblies],
                               decoy_accessions=[d["accession"] for d in lock.decoys], bam=str(bam), cached=False, timings_s=t,
                               discovered={str(k): v for k, v in sorted(discovered.items())})
    except Exception as exc:  # noqa: BLE001 - the lane reports its own failure
        res = GenomeLaneResult(status="failed", sketch_hits=[], taxa=[], ledger=None, index_accessions=[],
                               decoy_accessions=[], bam=None, cached=False, error=f"{type(exc).__name__}: {exc}"[:1500], timings_s=t)
    cached.write_text(json.dumps(res.to_json()))
    return res


def _from_json(d: dict[str, Any], *, cached: bool) -> GenomeLaneResult:
    hits = [SketchHit(accession=h["accession"], species_taxid=h.get("species_taxid"), species=h["species"],
                      ani=float(h["adjusted_ani"]), eff_cov=float(h["effective_coverage"]),
                      sequence_abundance=float(h["sequence_abundance_within_fungal_db"]), kmers=int(h.get("kmers", 0)))
            for h in (d.get("sketch") or {}).get("hits", [])]
    taxa = []
    for x in d.get("taxa", []):
        s = TaxonSupport(key=x["key"], rank=x["rank"], name=x["name"], species_taxid=x.get("species_taxid"), genus=x.get("genus", ""),
                         accessions=list(x.get("accessions", [])), fragments=int(x["fragments"]),
                         fragments_unique=int(x.get("fragments_unique_mapq", 0)), fragments_short_contigs=int(x.get("fragments_on_short_contigs", 0)),
                         discovery_fragments=x.get("discovery_fragments"), regions=int(x.get("separated_regions_100kb", 0)),
                         contigs=int(x.get("contigs_hit", 0)), top_window_fraction=float(x.get("top_20kb_window_fraction", 0.0)),
                         identity_mean=float(x.get("identity_mean", 0.0)), mapq_mean=float(x.get("mapq_mean", 0.0)),
                         genome_bases=int(x.get("backbone_bases", 0)), detection_state=x["detection_state"],
                         reasons=list(x.get("reasons", [])), backbone=x.get("backbone_accession"),
                         marker_lane=bool(x.get("also_in_marker_lane", False)), sketch_ani=x.get("sketch_adjusted_ani"))
        taxa.append(s)
    led = d.get("ledger")
    ledger = None
    if led:
        ledger = Ledger(eligible_fragments=int(led["eligible_fragments"]), aligned_fragments=int(led["aligned_to_index_fragments"]),
                        fungal_confident=int(led["fungal_supported_fragments"]), fungal_low_identity=int(led["fungal_below_identity_fragments"]),
                        cross_kingdom_ambiguous=int(led["cross_kingdom_ambiguous_fragments"]), decoy_fragments=int(led["decoy_fragments"]),
                        by_taxon_fragments={t.key: t.fragments for t in taxa})
    ci = d.get("competitive_index") or {}
    return GenomeLaneResult(status=d["status"], sketch_hits=hits, taxa=taxa, ledger=ledger,
                            index_accessions=list(ci.get("assemblies", [])), decoy_accessions=list(ci.get("decoys", [])),
                            bam=d.get("bam"), cached=cached, error=d.get("error"), timings_s=d.get("timings_s") or {},
                            discovered=dict((d.get("discovery") or {}).get("species_fragments") or {}))


def poisson_interval(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson interval for a fragment fraction k/n (sampling uncertainty only)."""
    if n <= 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


__all__ = ["CORE_PANEL_SPECIES", "GenomeLaneResult", "Ledger", "SUPPORT_RULE", "SketchHit", "TaxonSupport",
           "adjudicate", "align", "available", "build_index", "choose_index", "poisson_interval", "run_genome_lane",
           "sketch_candidates"]
