"""Fungal strain resolution, run for every supported fungus (spec §6.4).

For each supported species the reads already aligned competitively are
taken to the species backbone: callable nuclear bases and depth are
measured, the discriminating sites of the species' panel (where one exists)
are read out of the pileup, and the observed alleles are compared with
every panel member - reference genotypes, lineage groups - plus the
explicit "none of the panel" hypothesis. The output is the strongest
*supported* claim and the reason finer resolution was not reached:

    unresolved                 nothing discriminating was callable
    species_only               no panel exists for this species (a reference gap)
    lineage                    a lineage group's alleles are supported, no single genotype
    reference_equivalence_group  several deposited genotypes remain indistinguishable
    reference_genotype         one panel member is supported against the others
    sample_genotype            callable sites contradict every member: a novel genotype

At the depth a stool sample gives a low-abundance fungus (often well under
0.1×) the honest result is usually ``unresolved`` with the number of loci
searched and callable stated. That is a complete analysis, not a skipped
one; the difference is recorded in ``analysis_status``.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.mycobiome.genome_lane import TaxonSupport
from openbiota.mycobiome.panels import Panel, load_panel
from openbiota.mycobiome.references import Lock
from openbiota.mycobiome.strain_calls import pileup_sites, resolve

MIN_DEPTH_CALLABLE: Final = 3
MIN_BASEQ: Final = 20
MIN_MAPQ: Final = 20
#: Supporting blocks needed before a single reference genotype is claimed.
MIN_SUPPORT_BLOCKS: Final = 3
#: Contradicting blocks tolerated for a hypothesis to remain compatible.
CONTRADICTION_TOLERANCE: Final = 0
#: Minor-allele fraction across callable sites that suggests a mixture.
MIXTURE_MINOR_FRACTION: Final = 0.20
MIXTURE_MIN_SITES: Final = 3
CACHE_VERSION: Final = 1

REASONS: Final = frozenset({
    "insufficient_discriminatory_coverage", "reference_equivalence", "novel_or_unrepresented_genotype",
    "mixed_population_unresolved", "reference_panel_missing", "failed_qc", "execution_error",
})


@dataclass(slots=True)
class LocusCall:
    contig: str
    pos: int
    depth: int
    alleles: dict[str, int]
    block_id: int

    def to_json(self) -> dict[str, Any]:
        return {"contig": self.contig, "pos": self.pos, "depth": self.depth, "alleles": self.alleles, "block_id": self.block_id}


@dataclass(slots=True)
class StrainRecord:
    taxon_key: str
    species: str
    analysis_status: str = "not_assessed"        # not_assessed | complete | partial | failed
    resolution: str = "not_assessed"
    reason_codes: list[str] = field(default_factory=list)
    panel_lock_id: str | None = None
    panel_id: str | None = None
    backbone_accession: str | None = None
    reference_accessions: list[str] = field(default_factory=list)
    reference_genotypes: list[str] = field(default_factory=list)
    reference_equivalence_group: list[str] = field(default_factory=list)
    compatible_genotypes: list[str] = field(default_factory=list)
    lineage_id: str | None = None
    genotype_fingerprint_id: str | None = None
    ploidy_model: int | None = None
    genetic_code_id: int | None = None
    callable_nuclear_bases: int | None = None
    backbone_bases: int | None = None
    breadth_ge1: float | None = None
    mean_depth: float | None = None
    expected_depth: float | None = None
    eligible_discriminatory_loci: int | None = None
    callable_discriminatory_loci: int | None = None
    callable_blocks: int | None = None
    supporting_loci: list[LocusCall] = field(default_factory=list)
    contradictory_loci: list[LocusCall] = field(default_factory=list)
    hypothesis_scores: dict[str, dict[str, int]] = field(default_factory=dict)
    mixture_status: str = "not_assessed"
    mixture_components: list[dict[str, Any]] = field(default_factory=list)
    determinants: list[dict[str, Any]] = field(default_factory=list)
    phenotype_claims: list[dict[str, Any]] = field(default_factory=list)
    recommended_analytical_followup: list[str] = field(default_factory=list)
    error: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "taxon_key": self.taxon_key, "species": self.species, "analysis_status": self.analysis_status,
            "resolution": self.resolution, "reason_codes": self.reason_codes, "panel_lock_id": self.panel_lock_id,
            "panel_id": self.panel_id, "backbone_accession": self.backbone_accession,
            "reference_accessions": self.reference_accessions, "reference_genotypes": self.reference_genotypes,
            "reference_equivalence_group": self.reference_equivalence_group, "compatible_genotypes": self.compatible_genotypes,
            "lineage_id": self.lineage_id, "genotype_fingerprint_id": self.genotype_fingerprint_id,
            "ploidy_model": self.ploidy_model, "genetic_code_id": self.genetic_code_id,
            "callable_nuclear_bases": self.callable_nuclear_bases, "backbone_bases": self.backbone_bases,
            "breadth_ge1": self.breadth_ge1, "mean_depth": self.mean_depth, "expected_depth": self.expected_depth,
            "eligible_discriminatory_loci": self.eligible_discriminatory_loci,
            "callable_discriminatory_loci": self.callable_discriminatory_loci, "callable_blocks": self.callable_blocks,
            "supporting_loci": [x.to_json() for x in self.supporting_loci[:200]],
            "contradictory_loci": [x.to_json() for x in self.contradictory_loci[:200]],
            "hypothesis_scores": self.hypothesis_scores, "mixture_status": self.mixture_status,
            "mixture_components": self.mixture_components, "phase_blocks": [], "determinants": self.determinants,
            "phenotype_claims": self.phenotype_claims, "evidence_artifact_ids": [], "validation_release_id": None,
            "recommended_analytical_followup": self.recommended_analytical_followup, "error": self.error,
        }


def panels_dir(base: Path) -> Path:
    return base / "panels"


def find_panel(base: Path, species_taxid: int | None) -> Panel | None:
    if species_taxid is None:
        return None
    p = panels_dir(base) / f"{species_taxid}.panel.json"
    return load_panel(p) if p.is_file() else None


def _sorted_indexed(bam: Path, work: Path, st: str) -> Path:
    sbam = work / (bam.stem + ".sorted.bam")
    if not sbam.is_file() or not (sbam.with_suffix(".bam.bai")).is_file():
        subprocess.run([st, "sort", "-@", "8", "-o", str(sbam), str(bam)], check=True, capture_output=True)
        subprocess.run([st, "index", str(sbam)], check=True, capture_output=True)
    return sbam


def _coverage(sbam: Path, contigs: list[str], st: str, sites: set[tuple[str, int]] | None = None,
              ) -> tuple[int, int, float, set[tuple[str, int]]]:
    """(callable bases at depth>=MIN_DEPTH_CALLABLE, bases covered >=1, summed depth,
    the panel sites that are callable) over the contigs, in one streaming pass.

    ``samtools depth`` emits only covered positions, so at the depth a stool
    sample gives a fungus this is a short stream; reading the panel sites
    out of it means only callable sites are ever piled up."""
    if not contigs:
        return 0, 0, 0.0, set()
    proc = subprocess.Popen([st, "depth", "-Q", str(MIN_MAPQ), "-q", str(MIN_BASEQ), str(sbam)],
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    want = set(contigs)
    callable_ = covered = 0
    total = 0
    hit: set[tuple[str, int]] = set()
    assert proc.stdout is not None
    for line in proc.stdout:
        c, pos, d = line.rstrip("\n").split("\t")[:3]
        if c not in want:
            continue
        dd = int(d)
        if dd >= 1:
            covered += 1
            total += dd
        if dd >= MIN_DEPTH_CALLABLE:
            callable_ += 1
            if sites is not None:
                key = (c.split("|", 1)[1] if "|" in c else c, int(pos))
                if key in sites:
                    hit.add(key)
    proc.wait()
    return callable_, covered, total, hit


def assess(*, taxon: TaxonSupport, bam: Path, lock: Lock, base: Path, work_dir: Path, eligible_bases: float,
           fungal_fraction: float | None, log=print) -> StrainRecord:
    """The mandatory strain attempt for one supported taxon."""
    rec = StrainRecord(taxon_key=taxon.key, species=taxon.name)
    log(f"  strain attempt: {taxon.name}")
    st = shutil.which("samtools")
    if st is None:
        rec.analysis_status, rec.error = "failed", "samtools not found"
        rec.reason_codes.append("execution_error")
        return rec
    try:
        import pysam

        by_acc = lock.by_accession()
        bb = by_acc.get(taxon.backbone or "")
        if bb is None:
            rec.analysis_status = "failed"
            rec.reason_codes.append("reference_panel_missing")
            rec.error = "no backbone assembly for this taxon"
            return rec
        rec.backbone_accession = bb.accession
        rec.backbone_bases = bb.bases
        rec.genetic_code_id = bb.genetic_code
        lock.taxonomy.get(bb.taxid, {})
        rec.ploidy_model = _ploidy(bb.species)
        work_dir.mkdir(parents=True, exist_ok=True)
        sbam = _sorted_indexed(bam, work_dir, st)
        with pysam.AlignmentFile(str(sbam), "rb") as fh:
            contigs = [r for r in fh.references if r.startswith(bb.accession + "|")]
        panel = find_panel(base, taxon.species_taxid)
        site_keys = {(x.contig, x.pos) for x in panel.sites} if panel is not None else None
        callable_, covered, total, callable_sites = _coverage(sbam, contigs, st, site_keys)
        rec.callable_nuclear_bases = callable_
        rec.breadth_ge1 = (covered / bb.bases) if bb.bases else None
        rec.mean_depth = (total / bb.bases) if bb.bases else None
        if fungal_fraction is not None and bb.bases:
            rec.expected_depth = eligible_bases * fungal_fraction / bb.bases

        if panel is None:
            rec.analysis_status = "complete"
            rec.resolution = "species_only"
            rec.reason_codes.append("reference_panel_missing")
            rec.recommended_analytical_followup.append(
                "No discriminating marker panel exists for this species yet; intraspecies references from CGF/NCBI "
                "would need to be assembled into one before a strain claim is possible.")
            rec.mixture_status = "not_assessed"
            return rec
        rec.panel_id, rec.panel_lock_id = panel.panel_id, panel.lock_id
        if panel.backbone_accession != bb.accession:
            # The panel is on another backbone of the same species: use it.
            with pysam.AlignmentFile(str(sbam), "rb") as fh:
                contigs = [r for r in fh.references if r.startswith(panel.backbone_accession + "|")]
            if not contigs:
                rec.analysis_status = "partial"
                rec.resolution = "species_only"
                rec.reason_codes.append("reference_panel_missing")
                rec.error = "panel backbone was not in the competitive index"
                return rec
            rec.backbone_accession = panel.backbone_accession
        rec.eligible_discriminatory_loci = len(panel.sites)
        if panel.backbone_accession != bb.accession:
            _c, _v, _t, callable_sites = _coverage(sbam, contigs, st, site_keys)
        calls = pileup_sites(sbam, panel, rec.backbone_accession, only=callable_sites)
        return resolve(rec, panel, calls)
    except Exception as exc:  # noqa: BLE001 - a failed attempt is reported, never hidden
        rec.analysis_status = "failed"
        rec.reason_codes.append("execution_error")
        rec.error = f"{type(exc).__name__}: {exc}"[:800]
        return rec


def _ploidy(species: str) -> int:
    """Biological ploidy assumption per species; assembly representation is
    recorded separately. Diploid: Candida albicans, C. tropicalis,
    C. parapsilosis complex, Saccharomyces (often), Malassezia? no - haploid."""
    diploid = ("Candida albicans", "Candida tropicalis", "Candida parapsilosis", "Candida orthopsilosis",
               "Candida metapsilosis", "Candida dubliniensis", "Saccharomyces cerevisiae", "Candida auris", "Candidozyma auris")
    return 2 if species in diploid else 1


__all__ = ["CONTRADICTION_TOLERANCE", "MIN_DEPTH_CALLABLE", "MIN_SUPPORT_BLOCKS", "REASONS", "LocusCall", "StrainRecord",
           "assess", "find_panel", "panels_dir"]
