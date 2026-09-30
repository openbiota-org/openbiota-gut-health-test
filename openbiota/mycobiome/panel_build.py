"""Build the initial strain panels from the admitted references (spec §6.5).

Four panels, each derived from accessioned sequence and written to
``refs/mycobiome/panels/<species_taxid>.panel.json``:

* **Clavispora lusitaniae** - P4013B (GCA_023627835.1; experimentally
  protective isolate) against every other admitted C. lusitaniae assembly.
* **Saccharomyces cerevisiae** - the boulardii-like isolate assemblies
  (biocodex / I-745, ATCC MYA-796, unique28, EDRL and relatives) as one
  lineage group, with complete-genome S. cerevisiae backgrounds as members
  that must be *excluded* before a boulardii-like call is made.
* **Candida albicans** - the Li 2022 high- and low-damaging gut isolates
  genotyped from their public reads (diploid model) on SC5314, with the
  chromosome-level C. albicans assemblies as background members.
* **Candidozyma auris** - every admitted complete/chromosome assembly as a
  reference genotype on B8441; equivalence groups fall out of the data.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from openbiota.mycobiome.panels import Member, Panel, derive_assembly_panel, derive_reads_panel
from openbiota.mycobiome.references import Assembly, Lock, root

P4013B = "GCA_023627835.1"
S288C = "GCF_000146045.2"
SC5314 = "GCF_000182965.3"
B8441 = "GCA_002759435.3"
LUSITANIAE_TAXID, SCEREVISIAE_TAXID, ALBICANS_TAXID, AURIS_TAXID = 36911, 4932, 5476, 498019
BOULARDII_RE = re.compile(r"boulardii|I-745|CNCM|Unique28|EDRL|MYA-796|MYA-797|biocodex|kirkman|Unisankyo", re.I)
LI2022 = {  # verified run -> (isolate, phenotype)  (spec §6.5)
    "SRR13741117": ("IDB311", "high_damaging"), "SRR13741122": ("IDB101", "high_damaging"),
    "SRR13741104": ("IDC561", "low_damaging"), "SRR13741110": ("IDB891", "low_damaging"),
}
MAX_BACKGROUND = 40


def _backbone(lock: Lock, taxid: int, prefer: str | None = None) -> Assembly:
    asms = lock.species_assemblies().get(taxid, [])
    if prefer and any(a.accession == prefer for a in asms):
        return next(a for a in asms if a.accession == prefer)
    ref = [a for a in asms if a.source in ("refseq_reference", "both")]
    return (ref or sorted(asms, key=lambda a: -a.bases))[0]


def _write(panel: Panel, out_dir: Path) -> Path:
    p = out_dir / f"{panel.species_taxid}.panel.json"
    p.write_text(json.dumps(panel.to_json()))
    return p


def _bg(a: Assembly, label: str) -> Member:
    return Member(genotype_id=a.accession, kind="assembly", accession=a.accession, claim_scope="background",
                  label=f"{label} {a.strain}".strip(), group="", source="NCBI assembly")


def build_lusitaniae(lock: Lock, base: Path, out_dir: Path, *, threads: int, log=print) -> Path | None:
    lus = lock.species_assemblies().get(LUSITANIAE_TAXID, [])
    p4 = next((a for a in lus if a.accession == P4013B), None)
    if p4 is None:
        return None
    bb = _backbone(lock, LUSITANIAE_TAXID)
    members = {"P4013B": (base / p4.path, Member(
        genotype_id="P4013B", kind="assembly", accession=P4013B, claim_scope="reference_genotype",
        label="Clavispora lusitaniae P4013B (experimentally protective isolate)",
        phenotype={"endpoint": "colitis protection in mice via indole-3-ethanol/AHR", "direction": "favorable",
                   "label": "human_isolate_plus_experiment"},
        source="https://pmc.ncbi.nlm.nih.gov/articles/PMC12615640/"))}
    # Other C. lusitaniae isolates are background members: a P4013B call
    # requires them to be excluded, not merely P4013B to be nearest.
    for a in lus:
        if a.accession in (P4013B, bb.accession):
            continue
        members[a.accession] = (base / a.path, _bg(a, "C. lusitaniae"))
    panel = derive_assembly_panel(panel_id="clavispora_lusitaniae_p4013b_v1", species=bb.species,
                                  species_taxid=LUSITANIAE_TAXID, backbone=base / bb.path, backbone_accession=bb.accession,
                                  members=members, background={}, ploidy=1, genetic_code=bb.genetic_code,
                                  lock_id=lock.lock_id, threads=threads, log=log)
    return _write(panel, out_dir)


def build_saccharomyces(lock: Lock, base: Path, out_dir: Path, *, threads: int, log=print) -> Path | None:
    sc = lock.species_assemblies().get(SCEREVISIAE_TAXID, [])
    if not sc:
        return None
    bb = _backbone(lock, SCEREVISIAE_TAXID, prefer=S288C)
    members: dict = {}
    n_bg = 0
    for a in sc:
        if a.accession == bb.accession:
            continue
        if BOULARDII_RE.search(f"{a.strain} {a.organism}"):
            members[a.accession] = (base / a.path, Member(
                genotype_id=a.accession, kind="assembly", accession=a.accession, claim_scope="reference_genotype",
                label=f"S. cerevisiae boulardii-like isolate {a.strain}".strip(), group="boulardii_like",
                source="NCBI assembly; a lineage call is not a product identity"))
        elif a.level == "Complete Genome" and n_bg < MAX_BACKGROUND:
            members[a.accession] = (base / a.path, _bg(a, "S. cerevisiae"))
            n_bg += 1
    if not members:
        return None
    panel = derive_assembly_panel(panel_id="saccharomyces_cerevisiae_boulardii_lineage_v1", species=bb.species,
                                  species_taxid=SCEREVISIAE_TAXID, backbone=base / bb.path, backbone_accession=bb.accession,
                                  members=members, background={}, ploidy=2, genetic_code=bb.genetic_code,
                                  lock_id=lock.lock_id, threads=threads, log=log)
    return _write(panel, out_dir)


def build_albicans(lock: Lock, base: Path, out_dir: Path, *, threads: int, log=print) -> Path | None:
    """HD/LD isolates from reads (diploid), backgrounds from assemblies."""
    alb = lock.species_assemblies().get(ALBICANS_TAXID, [])
    reads_dir = base / "isolate_reads" / "li2022"
    man = reads_dir / "manifest.tsv"
    if not alb or not man.is_file():
        return None
    bb = _backbone(lock, ALBICANS_TAXID, prefer=SC5314)
    runs: dict[str, list[Path]] = {}
    for row in man.read_text().splitlines()[1:]:
        run, _iso, _ph, fname, *_ = row.split("\t")
        runs.setdefault(run, []).append(reads_dir / fname)
    members: dict = {}
    for run, (iso, ph) in LI2022.items():
        files = sorted(runs.get(run, []))
        if len(files) < 2 or not all(f.is_file() for f in files):
            log(f"  Li 2022 run {run} ({iso}) not present; skipped")
            continue
        members[iso] = (files[0], files[1], Member(
            genotype_id=iso, kind="reads", accession=run, claim_scope="reference_genotype",
            label=f"C. albicans gut isolate {iso} ({ph.replace('_', '-')})", group=ph,
            phenotype={"assay": "immune-cell damage (Li 2022)", "class": ph,
                       "note": "study assay label, not measured in any participant"},
            source="https://pmc.ncbi.nlm.nih.gov/articles/PMC9166917/"))
    if not members:
        return None
    panel = derive_reads_panel(panel_id="candida_albicans_li2022_hd_ld_v1", species=bb.species, species_taxid=ALBICANS_TAXID,
                               backbone=base / bb.path, backbone_accession=bb.accession, members=members, ploidy=2,
                               genetic_code=bb.genetic_code, lock_id=lock.lock_id, work_dir=base / "panel_work" / "albicans",
                               threads=threads, log=log)
    # Background assemblies: chromosome-level C. albicans, genotyped at the
    # same sites so a common polymorphism is not read as an isolate marker.
    from openbiota.mycobiome.panels import assembly_substitutions

    bgs = [a for a in alb if a.level in ("Complete Genome", "Chromosome") and a.accession != bb.accession][:MAX_BACKGROUND]
    for a in bgs:
        subs = assembly_substitutions(base / bb.path, base / a.path, threads)
        panel.members.append(_bg(a, "C. albicans"))
        for s in panel.sites:
            hit = subs.get((s.contig, s.pos))
            s.alleles[a.accession] = hit[1] if hit else "R"
    return _write(panel, out_dir)


def build_auris(lock: Lock, base: Path, out_dir: Path, *, threads: int, log=print) -> Path | None:
    aur = lock.species_assemblies().get(AURIS_TAXID, [])
    if not aur:
        return None
    bb = _backbone(lock, AURIS_TAXID, prefer=B8441)
    members = {
        a.accession: (base / a.path, Member(genotype_id=a.accession, kind="assembly", accession=a.accession,
                                             claim_scope="reference_genotype", label=f"C. auris {a.strain}".strip(),
                                             group="", source="NCBI assembly; clade not asserted without a verified source"))
        for a in aur if a.accession != bb.accession and a.level in ("Complete Genome", "Chromosome")
    }
    if not members:
        return None
    panel = derive_assembly_panel(panel_id="candidozyma_auris_reference_genotypes_v1", species=bb.species, species_taxid=AURIS_TAXID,
                                  backbone=base / bb.path, backbone_accession=bb.accession, members=members, background={},
                                  ploidy=1, genetic_code=bb.genetic_code, lock_id=lock.lock_id, threads=threads, log=log)
    return _write(panel, out_dir)


def build_all(lock: Lock, base: Path | None = None, *, threads: int = 16, log=print) -> list[Path]:
    base = base or root()
    out_dir = base / "panels"
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for fn in (build_lusitaniae, build_saccharomyces, build_auris, build_albicans):
        try:
            p = fn(lock, base, out_dir, threads=threads, log=log)
            if p:
                written.append(p)
                log(f"  wrote {p.name}")
        except Exception as exc:  # noqa: BLE001 - one panel failing must not lose the others
            log(f"  panel {fn.__name__} failed: {type(exc).__name__}: {exc}")
    return written


__all__ = ["build_all", "build_albicans", "build_auris", "build_lusitaniae", "build_saccharomyces"]
