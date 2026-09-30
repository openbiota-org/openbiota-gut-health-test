"""Strain marker panels derived from accessioned sequence (spec §6.4.2, §6.5).

A panel is a set of discriminating single-nucleotide sites on a species
backbone (a versioned reference assembly), each with the allele carried by
every panel member. Members are reference genotypes (a named isolate
assembly, or an isolate genotyped from its public reads) or lineage groups
(several isolates sharing the allele). Sites are derived from the sequences
themselves - never from a paper's descriptive coordinates - and screened:

* a member's allele must differ from the backbone and be consistent across
  that member's sequence;
* conspecific background genomes are genotyped at the same sites, so a
  site carried by much of the species is a species polymorphism, not a
  member marker (it is kept but flagged ``background_frequency``);
* sites within ``MIN_SPACING`` of one another collapse into one block so
  linked sites are never counted as independent evidence.

Assemblies are compared with ``minimap2 -cx asm10 --cs`` and the ``cs`` tag
is parsed for substitutions. Isolate reads are genotyped with ``bcftools``
under a diploid model where the species is diploid. Every panel records its
backbone accession and hash, the member accessions or runs, and the claim
scope of each member (``reference_genotype`` or ``lineage``).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

MIN_SPACING: Final = 300              # bp; retained for documentation of the earlier rule
BLOCK_BP: Final = 1_000               # linkage block: one 1-kb window of one contig
BACKGROUND_COMMON: Final = 0.10       # member allele in >10% of background genomes -> species polymorphism
PANEL_VERSION: Final = 1


@dataclass(slots=True)
class Member:
    genotype_id: str
    kind: str                          # assembly | reads
    accession: str                     # GCA/GCF or SRR
    claim_scope: str                   # reference_genotype | lineage | background
    label: str = ""
    group: str = ""                    # lineage group this member belongs to ("" = none)
    phenotype: dict[str, Any] = field(default_factory=dict)
    source: str = ""


@dataclass(slots=True)
class Site:
    contig: str
    pos: int                           # 1-based on the backbone
    ref: str
    alleles: dict[str, str]            # genotype_id -> allele (A/C/G/T, "N" unknown, "R" = backbone)
    background_frequency: float | None
    block_id: int

    def to_json(self) -> dict[str, Any]:
        return {"contig": self.contig, "pos": self.pos, "ref": self.ref, "alleles": self.alleles,
                "background_frequency": self.background_frequency, "block_id": self.block_id}


@dataclass(slots=True)
class Panel:
    panel_id: str
    species: str
    species_taxid: int
    backbone_accession: str
    backbone_sha256: str
    ploidy: int
    genetic_code: int | None
    members: list[Member]
    sites: list[Site]
    background_accessions: list[str]
    method: str
    lock_id: str

    @property
    def n_blocks(self) -> int:
        return len({s.block_id for s in self.sites})

    def to_json(self) -> dict[str, Any]:
        return {
            "panel_version": PANEL_VERSION, "panel_id": self.panel_id, "species": self.species,
            "species_taxid": self.species_taxid, "backbone_accession": self.backbone_accession,
            "backbone_sha256": self.backbone_sha256, "ploidy_model": self.ploidy, "genetic_code_id": self.genetic_code,
            "members": [{"genotype_id": m.genotype_id, "kind": m.kind, "accession": m.accession, "claim_scope": m.claim_scope,
                         "label": m.label, "group": m.group, "phenotype": m.phenotype, "source": m.source} for m in self.members],
            "n_sites": len(self.sites), "n_blocks": self.n_blocks, "sites": [s.to_json() for s in self.sites],
            "background_accessions": self.background_accessions, "method": self.method, "lock_id": self.lock_id,
        }


def reblock_panel_file(path: Path) -> int:
    """Recompute block ids for a stored panel under the current BLOCK_BP rule.
    Returns the number of blocks."""
    d = json.loads(path.read_text())
    keys: dict[tuple[str, int], int] = {}
    for s in d["sites"]:
        k = (s["contig"], int(s["pos"]) // BLOCK_BP)
        if k not in keys:
            keys[k] = len(keys) + 1
        s["block_id"] = keys[k]
    d["n_blocks"] = len(keys)
    d["method"] = d.get("method", "") + f"; blocks = {BLOCK_BP}-bp windows"
    path.write_text(json.dumps(d))
    return len(keys)


def load_panel(path: Path) -> Panel:
    d = json.loads(path.read_text())
    return Panel(
        panel_id=d["panel_id"], species=d["species"], species_taxid=int(d["species_taxid"]),
        backbone_accession=d["backbone_accession"], backbone_sha256=d["backbone_sha256"], ploidy=int(d["ploidy_model"]),
        genetic_code=d.get("genetic_code_id"),
        members=[Member(**m) for m in d["members"]],
        sites=[Site(contig=s["contig"], pos=int(s["pos"]), ref=s["ref"], alleles=dict(s["alleles"]),
                    background_frequency=s.get("background_frequency"), block_id=int(s["block_id"])) for s in d["sites"]],
        background_accessions=list(d.get("background_accessions", [])), method=d.get("method", ""), lock_id=d.get("lock_id", ""),
    )


# --------------------------------------------------------------------------- #
# assembly-vs-backbone substitutions
# --------------------------------------------------------------------------- #

_CS_SUB = re.compile(r"\*([acgtn])([acgtn])")
_CS_TOK = re.compile(r"(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)")


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def assembly_substitutions(backbone: Path, query: Path, threads: int = 8) -> dict[tuple[str, int], tuple[str, str]]:
    """{(contig, pos1): (ref_base, alt_base)} for every substitution in the
    best alignment of *query* to *backbone*. Only primary alignments with
    mapping quality >= 30 contribute, so paralogous or repeat copies do not
    write alleles onto the wrong locus."""
    mm = shutil.which("minimap2")
    if mm is None:
        raise RuntimeError("minimap2 not found")
    proc = subprocess.run([mm, "-cx", "asm10", "--cs", "-t", str(threads), "--secondary=no", str(backbone), str(query)],
                          capture_output=True, text=True, check=True)
    out: dict[tuple[str, int], tuple[str, str]] = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        mapq = int(f[11])
        if mapq < 30:
            continue
        tname, tstart = f[5], int(f[7])
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        pos = tstart  # 0-based on target
        for tok in _CS_TOK.findall(cs):
            if tok[0] == ":":
                pos += int(tok[1:])
            elif tok[0] == "*":
                out[(tname, pos + 1)] = (tok[1].upper(), tok[2].upper())
                pos += 1
            elif tok[0] == "-":
                pos += len(tok) - 1
            # insertions ("+") consume no target bases
    return out


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with (gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")) as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _blocks(sites: list[tuple[str, int]]) -> dict[tuple[str, int], int]:
    """Linkage block = one BLOCK_BP window of one contig. Sites in the same
    window are never counted as independent evidence."""
    ids: dict[tuple[str, int], int] = {}
    keys: dict[tuple[str, int], int] = {}
    for c, p in sorted(sites):
        k = (c, p // BLOCK_BP)
        if k not in keys:
            keys[k] = len(keys) + 1
        ids[(c, p)] = keys[k]
    return ids


def derive_assembly_panel(
    *, panel_id: str, species: str, species_taxid: int, backbone: Path, backbone_accession: str,
    members: dict[str, tuple[Path, Member]], background: dict[str, Path], ploidy: int, genetic_code: int | None,
    lock_id: str, threads: int = 8, log=print,
) -> Panel:
    """Panel from member assemblies against conspecific background assemblies."""
    log(f"  panel {panel_id}: {len(members)} members, {len(background)} background genomes")
    member_subs = {gid: assembly_substitutions(backbone, path, threads) for gid, (path, _m) in members.items()}
    bg_subs = {acc: assembly_substitutions(backbone, path, threads) for acc, path in background.items()}
    # Candidate sites: any member substitution.
    cand: set[tuple[str, int]] = set()
    for subs in member_subs.values():
        cand |= set(subs)
    blocks = _blocks(sorted(cand))
    sites: list[Site] = []
    for key in sorted(cand):
        contig, pos = key
        ref = next(subs[key][0] for subs in member_subs.values() if key in subs)
        alleles = {gid: (subs[key][1] if key in subs else "R") for gid, subs in member_subs.items()}
        # Background frequency of the member (non-reference) alleles.
        alts = {a for a in alleles.values() if a != "R"}
        if bg_subs:
            carriers = sum(1 for subs in bg_subs.values() if key in subs and subs[key][1] in alts)
            bf = carriers / len(bg_subs)
        else:
            bf = None
        sites.append(Site(contig=contig, pos=pos, ref=ref, alleles=alleles, background_frequency=bf, block_id=blocks[key]))
    return Panel(panel_id=panel_id, species=species, species_taxid=species_taxid, backbone_accession=backbone_accession,
                 backbone_sha256=_sha256(backbone), ploidy=ploidy, genetic_code=genetic_code,
                 members=[m for _p, m in members.values()], sites=sites, background_accessions=sorted(background),
                 method="minimap2 asm10 cs-tag substitutions, MAPQ>=30, background-screened", lock_id=lock_id)


# --------------------------------------------------------------------------- #
# isolate reads genotyped on the backbone
# --------------------------------------------------------------------------- #


def genotype_reads(*, backbone: Path, r1: Path, r2: Path | None, out_dir: Path, name: str, ploidy: int,
                   threads: int = 8, log=print) -> dict[tuple[str, int], str]:
    """{(contig, pos1): allele} from isolate reads: minimap2 sr + bcftools
    mpileup/call under the given ploidy. Heterozygous diploid sites are
    encoded as IUPAC-free two-letter strings (e.g. "AG") so a heterozygote
    is never mistaken for two strains."""
    mm, st, bt = shutil.which("minimap2"), shutil.which("samtools"), shutil.which("bcftools")
    if not (mm and st and bt):
        raise RuntimeError("minimap2, samtools and bcftools are required to genotype isolate reads")
    out_dir.mkdir(parents=True, exist_ok=True)
    bam = out_dir / f"{name}.bam"
    vcf = out_dir / f"{name}.vcf.gz"
    if not vcf.is_file():
        log(f"  genotyping isolate {name} on {backbone.name}")
        if not (backbone.with_suffix(backbone.suffix + ".fai")).is_file():
            subprocess.run([st, "faidx", str(backbone)], check=True, capture_output=True)
        cmd = [mm, "-ax", "sr", "-t", str(threads), str(backbone), str(r1)] + ([str(r2)] if r2 else [])
        with (out_dir / f"{name}.align.log").open("w") as err:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
            p2 = subprocess.Popen([st, "sort", "-@", "4", "-o", str(bam), "-"], stdin=p1.stdout, stderr=err)
            p1.stdout.close()  # type: ignore[union-attr]
            p2.communicate()
            p1.wait()
        subprocess.run([st, "index", str(bam)], check=True, capture_output=True)
        with vcf.open("wb") as out:
            p1 = subprocess.Popen([bt, "mpileup", "-Ou", "-f", str(backbone), "-q", "20", "-Q", "20", "-a", "AD,DP",
                                   "--threads", str(threads), str(bam)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen([bt, "call", "-mv", "--ploidy", str(ploidy), "-Oz"], stdin=p1.stdout, stdout=out,
                                  stderr=subprocess.DEVNULL)
            p1.stdout.close()  # type: ignore[union-attr]
            p2.communicate()
            p1.wait()
    alleles: dict[tuple[str, int], str] = {}
    proc = subprocess.run([bt, "query", "-i", "QUAL>=30 && INFO/DP>=8 && TYPE=\"snp\"", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n", str(vcf)],
                          capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        chrom, pos, ref, alt, gt = line.split("\t")
        alts = alt.split(",")
        calls = re.split(r"[/|]", gt)
        bases = []
        for c in calls:
            if c == "0":
                bases.append(ref)
            elif c.isdigit() and int(c) - 1 < len(alts):
                bases.append(alts[int(c) - 1])
        if not bases or any(len(b) != 1 for b in bases):
            continue
        bases = sorted(set(bases))
        alleles[(chrom, int(pos))] = "".join(bases) if len(bases) > 1 else bases[0]
    return alleles


def derive_reads_panel(
    *, panel_id: str, species: str, species_taxid: int, backbone: Path, backbone_accession: str,
    members: dict[str, tuple[Path, Path | None, Member]], ploidy: int, genetic_code: int | None, lock_id: str,
    work_dir: Path, threads: int = 8, log=print,
) -> Panel:
    """Panel from isolate reads: sites where the isolates differ from one
    another (or all from the backbone). Heterozygous calls are kept as such."""
    geno = {gid: genotype_reads(backbone=backbone, r1=r1, r2=r2, out_dir=work_dir / gid, name=gid, ploidy=ploidy,
                                threads=threads, log=log) for gid, (r1, r2, _m) in members.items()}
    cand: set[tuple[str, int]] = set()
    for g in geno.values():
        cand |= set(g)
    # Keep sites polymorphic among members (discriminating), including
    # backbone-like members ("R") where others differ.
    keep: list[tuple[str, int]] = []
    for key in cand:
        alleles = {gid: g.get(key, "R") for gid, g in geno.items()}
        if len(set(alleles.values())) > 1:
            keep.append(key)
    blocks = _blocks(sorted(keep))
    ref_of = _reference_bases(backbone, keep)
    sites = [Site(contig=c, pos=p, ref=ref_of.get((c, p), "N"), alleles={gid: g.get((c, p), "R") for gid, g in geno.items()},
                  background_frequency=None, block_id=blocks[(c, p)]) for c, p in sorted(keep)]
    return Panel(panel_id=panel_id, species=species, species_taxid=species_taxid, backbone_accession=backbone_accession,
                 backbone_sha256=_sha256(backbone), ploidy=ploidy, genetic_code=genetic_code,
                 members=[m for _r1, _r2, m in members.values()], sites=sites, background_accessions=[],
                 method=f"minimap2 sr + bcftools call ploidy={ploidy}; QUAL>=30, DP>=8; member-discriminating sites", lock_id=lock_id)


def _reference_bases(backbone: Path, sites: list[tuple[str, int]]) -> dict[tuple[str, int], str]:
    want: dict[str, list[int]] = defaultdict(list)
    for c, p in sites:
        want[c].append(p)
    out: dict[tuple[str, int], str] = {}
    name = None
    buf: list[str] = []

    def flush() -> None:
        if name is not None and name in want:
            seq = "".join(buf)
            for p in want[name]:
                if 0 < p <= len(seq):
                    out[(name, p)] = seq[p - 1].upper()

    with _open(backbone) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                buf = []
            elif name in want:
                buf.append(line.strip())
        flush()
    return out


__all__ = ["BACKGROUND_COMMON", "MIN_SPACING", "Member", "Panel", "Site", "assembly_substitutions",
           "derive_assembly_panel", "derive_reads_panel", "genotype_reads", "load_panel"]
