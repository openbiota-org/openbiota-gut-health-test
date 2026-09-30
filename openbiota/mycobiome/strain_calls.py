"""Reading panel sites out of a pileup and comparing genotype hypotheses.

Split from :mod:`strains` so each half stays readable. ``pileup_sites``
returns the callable discriminating loci with allele counts; ``resolve``
scores every panel member against them by *linkage block* (sites within
one block count once), keeps the members with no contradicting block as
the compatible set, and names the resolution the evidence actually
supports. Heterozygous calls at a diploid site are not treated as two
strains; a minor allele across several blocks is reported as mixture
evidence, ambiguous under a diploid model.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from openbiota.mycobiome.panels import Panel

if TYPE_CHECKING:
    from openbiota.mycobiome.strains import LocusCall, StrainRecord

MIN_DEPTH_CALLABLE = 3
MIN_BASEQ = 20
MIN_MAPQ = 20
MIN_SUPPORT_BLOCKS = 3
CONTRADICTION_TOLERANCE = 0
MIXTURE_MINOR_FRACTION = 0.20
MIXTURE_MIN_SITES = 3


def pileup_sites(sbam: Path, panel: Panel, backbone_accession: str,
                 only: set[tuple[str, int]] | None = None) -> list[LocusCall]:
    """Allele counts at the panel sites that have coverage. ``only`` is the
    set of (contig, pos) already known to be callable from the depth pass;
    without it every site is visited, which is slow for a 500k-site panel."""
    import pysam

    from openbiota.mycobiome.strains import LocusCall

    calls: list[LocusCall] = []
    sites = panel.sites if only is None else [s for s in panel.sites if (s.contig, s.pos) in only]
    with pysam.AlignmentFile(str(sbam), "rb") as fh:
        refs = set(fh.references)
        for s in sites:
            ref = f"{backbone_accession}|{s.contig}"
            if ref not in refs:
                continue
            counts: dict[str, int] = defaultdict(int)
            for col in fh.pileup(ref, s.pos - 1, s.pos, truncate=True, min_base_quality=MIN_BASEQ,
                                 min_mapping_quality=MIN_MAPQ, stepper="samtools"):
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None:
                        continue
                    base = pr.alignment.query_sequence[pr.query_position].upper()
                    if base in "ACGT":
                        counts[base] += 1
            depth = sum(counts.values())
            if depth >= MIN_DEPTH_CALLABLE:
                calls.append(LocusCall(contig=s.contig, pos=s.pos, depth=depth, alleles=dict(counts), block_id=s.block_id))
    return calls


def _member_allele(site_allele: str, ref: str) -> set[str]:
    """The bases a member is expected to show at a site: "R" means backbone;
    a two-letter string is a heterozygous diploid genotype."""
    if site_allele == "R":
        return {ref}
    return set(site_allele)


def resolve(rec: StrainRecord, panel: Panel, calls: list[LocusCall]) -> StrainRecord:
    rec.callable_discriminatory_loci = len(calls)
    rec.callable_blocks = len({c.block_id for c in calls})
    rec.analysis_status = "complete"
    rec.reference_accessions = [m.accession for m in panel.members if m.kind == "assembly"]
    rec.reference_genotypes = [m.genotype_id for m in panel.members]
    members = {m.genotype_id: m for m in panel.members}

    if not calls:
        rec.resolution = "unresolved"
        rec.reason_codes.append("insufficient_discriminatory_coverage")
        rec.compatible_genotypes = list(members)
        rec.reference_equivalence_group = list(members)
        rec.mixture_status = "not_assessed"
        rec.recommended_analytical_followup.append(followup_text(rec, panel))
        return rec

    site_by_pos = {(x.contig, x.pos): x for x in panel.sites}
    support: dict[str, set[int]] = {g: set() for g in members}
    contradict: dict[str, set[int]] = {g: set() for g in members}
    minor_sites = 0
    for c in calls:
        site = site_by_pos[(c.contig, c.pos)]
        total = sum(c.alleles.values())
        observed = {a for a, n in c.alleles.items() if n / total >= MIXTURE_MINOR_FRACTION}
        if len(observed) > 1:
            minor_sites += 1
        for gid in members:
            expected = _member_allele(site.alleles.get(gid, "R"), site.ref)
            # Supported when the observed alleles are within the expected set;
            # contradicted when an observed allele is outside it.
            if observed <= expected:
                support[gid].add(c.block_id)
            else:
                contradict[gid].add(c.block_id)
    rec.hypothesis_scores = {g: {"supporting_blocks": len(support[g]), "contradicting_blocks": len(contradict[g])} for g in members}
    compatible = [g for g in members if len(contradict[g]) <= CONTRADICTION_TOLERANCE]
    rec.compatible_genotypes = compatible

    # Mixture evidence, read under the species' ploidy model.
    if minor_sites >= MIXTURE_MIN_SITES:
        rec.mixture_status = "ambiguous" if (rec.ploidy_model or 1) >= 2 else "mixture_supported"
        rec.mixture_components = [{"note": f"minor alleles at {minor_sites} callable sites",
                                   "interpretation": ("heterozygosity or a conspecific mixture; not separable here"
                                                      if rec.mixture_status == "ambiguous" else "two or more genotypes")}]
        if rec.mixture_status == "mixture_supported":
            rec.reason_codes.append("mixed_population_unresolved")
    else:
        rec.mixture_status = "no_mixture_evidence_at_achieved_sensitivity"

    supporting = [c for c in calls if any(c.block_id in support[g] for g in compatible)] if compatible else []
    rec.supporting_loci = supporting
    rec.contradictory_loci = [c for c in calls if all(c.block_id in contradict[g] for g in members)]

    # Name the resolution the evidence supports.
    if not compatible:
        rec.resolution = "sample_genotype"
        rec.reason_codes.append("novel_or_unrepresented_genotype")
        rec.genotype_fingerprint_id = _fingerprint(calls)
    else:
        # Blocks that actually discriminate the compatible members from the rest.
        others = [g for g in members if g not in compatible]
        disc = {c.block_id for c in calls if any(c.block_id in contradict[g] for g in others)} if others else set()
        scopes = {members[g].claim_scope for g in compatible}
        groups = {members[g].group for g in compatible}
        if compatible and scopes == {"background"} and len(disc) >= MIN_SUPPORT_BLOCKS:
            # Only background genomes remain: the sample is not any studied member.
            rec.resolution = "lineage"
            rec.lineage_id = "background_not_studied_members"
            rec.reference_equivalence_group = compatible
            rec.reason_codes.append("reference_equivalence")
        elif len(compatible) == 1 and len(disc) >= MIN_SUPPORT_BLOCKS:
            rec.resolution = "reference_genotype" if scopes == {"reference_genotype"} else "lineage"
            if rec.resolution == "lineage":
                rec.lineage_id = compatible[0]
        elif len(groups) == 1 and next(iter(groups)) and "background" not in scopes and len(disc) >= MIN_SUPPORT_BLOCKS:
            # Every compatible member is one lineage group and the rest were excluded.
            rec.resolution = "lineage"
            rec.lineage_id = next(iter(groups))
            rec.reference_equivalence_group = compatible
            rec.reason_codes.append("reference_equivalence")
        elif len(compatible) == len(members) or len(disc) < MIN_SUPPORT_BLOCKS:
            rec.resolution = "unresolved"
            rec.reason_codes.append("insufficient_discriminatory_coverage")
            rec.reference_equivalence_group = compatible
        else:
            rec.resolution = "reference_equivalence_group" if "reference_genotype" in scopes else "lineage"
            rec.reference_equivalence_group = compatible
            rec.reason_codes.append("reference_equivalence")
            if rec.resolution == "lineage":
                rec.lineage_id = "+".join(compatible)
    rec.recommended_analytical_followup.append(followup_text(rec, panel))
    return rec


def _fingerprint(calls: list[LocusCall]) -> str:
    import hashlib

    blob = ";".join(f"{c.contig}:{c.pos}:{''.join(sorted(a for a, n in c.alleles.items() if n))}" for c in calls)
    return "fp-" + hashlib.sha256(blob.encode()).hexdigest()[:12]


def followup_text(rec: StrainRecord, panel: Panel) -> str:
    n = rec.eligible_discriminatory_loci or 0
    k = rec.callable_discriminatory_loci or 0
    depth = rec.mean_depth or 0.0
    if k == 0:
        return (f"{n} discriminating sites in {panel.n_blocks} linkage blocks were searched and none was callable at "
                f"{depth:.3f}x mean depth. Deeper sequencing of this sample at the same fungal fraction, a targeted "
                "fungal capture, or culture and isolate sequencing would be needed to place this strain.")
    if rec.resolution in ("reference_genotype", "lineage"):
        return f"{k} of {n} discriminating sites callable in {rec.callable_blocks} blocks; the call rests on those blocks."
    return (f"{k} of {n} discriminating sites callable in {rec.callable_blocks} blocks - enough to exclude some "
            f"genotypes, not to separate {len(rec.reference_equivalence_group)} that remain compatible.")


__all__ = ["followup_text", "pileup_sites", "resolve"]
