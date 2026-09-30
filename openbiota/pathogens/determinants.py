"""Toxin, virulence and resistance determinant screening (spec §5.6, §8.2-8.3).

Two questions, kept apart at every stage:

1. **Is the gene sequence there?** A nucleotide question, answered by
   competitive alignment against NDARO reference families exactly as
   organisms are.
2. **Which organism carries it?** A *linkage* question, answered only by
   physical co-location — the same read pair or the same assembled contig
   as a target's sequence. Nothing else counts. A stool sample contains
   thousands of genomes; a resistance gene found in it belongs to one of
   them, and co-occurrence in the same tube is not evidence of which.

Collapsing the two is the classic failure of resistome reporting: `blaKPC`
present plus *Klebsiella pneumoniae* present becomes "carbapenem-resistant
Klebsiella", which the data does not support and which changes treatment.
This module therefore reports `unlinked` as a first-class outcome with its
own wording, never as a weaker version of a linked call.

Three further rules the spec is explicit about:

* Gene presence is not a phenotype. `clinical_phenotype_inference` stays
  `not_inferred` — resistance genes can be unexpressed, truncated or on a
  plasmid that is lost, and susceptibility testing is a different assay.
* Allele-level claims need allele-level evidence. A family whose members
  differ in clinical meaning (`blaOXA` carbapenemases versus the rest) is
  reported at family level with `requires_allele_resolution` recorded, not
  promoted to the clinically alarming member.
* Locus completeness is reported, not assumed. Finding `ybtP` is not
  finding the yersiniabactin locus, and an operon missing its B subunit is
  reported as partial.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.pathogens.align import TargetEvidence, iter_sam
from openbiota.pathogens.catalog import DeterminantSeed
from openbiota.pathogens.schema import DeterminantResult, PathogenError

__all__ = [
    "DeterminantPolicy",
    "NdaroMap",
    "adjudicate_determinant",
    "load_ndaro_map",
    "screen_determinants",
    "self_test",
]

SCHEMA_VERSION: Final = "5.0"

#: Linkage states that actually name a carrier organism. `unlinked` and
#: `ambiguous` do not, and must never be rendered as though they did.
LINKED_STATES: Final[frozenset[str]] = frozenset(
    {"contig_supported", "read_pair_supported", "genome_supported"}
)

#: Determinant statuses that mean "this gene's sequence is present".
#: `ambiguous_allele` is included: the family is there, but which member is
#: not resolved, so it is a finding that cannot carry an allele-level claim.
PRESENT_STATUSES: Final[frozenset[str]] = frozenset(
    {
        "supported_intact_sequence",
        "supported_partial_sequence",
        "ambiguous_allele",
    }
)


@dataclass(frozen=True)
class DeterminantPolicy:
    """Support thresholds for gene-level calls.

    Deliberately stricter per base than the organism policy: a resistance
    gene is short, often high-identity across families, and frequently
    present on mobile elements shared between unrelated organisms, so a
    handful of reads over a 300 bp window is not a call.
    """

    min_distinct_fragments: int = 2
    min_covered_fraction: float = 0.60
    min_identity: float = 0.95
    #: Families whose members differ in clinical meaning need this much of
    #: the reference before even a family-level statement is made.
    min_covered_fraction_heterogeneous: float = 0.80


@dataclass(frozen=True)
class NdaroFamilyMapping:
    """One determinant's mapping onto installed NDARO gene families."""

    determinant_id: str
    gene_families: tuple[str, ...]
    match_mode: str
    #: Prefix mode names a *stem* rather than exact families, because NDARO
    #: splits some determinants across many numbered families (ipaH has ten,
    #: tet has seventy-one). The stem is resolved against the installed
    #: family list at reference-build time.
    family_prefixes: tuple[str, ...] = ()
    note: str = ""
    requires_allele_resolution: bool = False
    partial_locus: bool = False
    gap_reason: str | None = None

    @property
    def available(self) -> bool:
        return self.match_mode != "absent" and bool(
            self.gene_families or self.family_prefixes
        )

    def matches(self, family: str) -> bool:
        """Whether an installed NDARO family belongs to this determinant."""
        if family in self.gene_families:
            return True
        return self.match_mode == "prefix" and any(
            family.startswith(p) for p in self.family_prefixes
        )


@dataclass(frozen=True)
class NdaroMap:
    """The frozen determinant-to-reference mapping.

    Explicit by design. Fuzzy name matching against NDARO produces
    false positives that are hard to see afterwards — `sta` is a
    streptothricin acetyltransferase, not the ETEC heat-stable enterotoxin,
    and a prefix search for `sul` picks up sulfonamide-degrading enzymes
    rather than acquired dihydropteroate synthases. Every family here was
    verified against the catalogue by name.
    """

    catalog_version: str
    by_determinant: Mapping[str, NdaroFamilyMapping]
    installed_family_count: int = 0

    def for_determinant(self, determinant_id: str) -> NdaroFamilyMapping | None:
        return self.by_determinant.get(determinant_id)


def load_ndaro_map(path: Path) -> NdaroMap:
    """Load and validate `ndaro_map.yaml`."""
    path = Path(path)
    if not path.is_file():
        raise PathogenError(f"NDARO mapping not found: {path}")
    doc = yaml.safe_load(path.read_text()) or {}
    raw = doc.get("mappings") or []
    if not raw:
        raise PathogenError(f"{path}: no mappings")

    out: dict[str, NdaroFamilyMapping] = {}
    for entry in raw:
        det_id = str(entry.get("determinant_id") or "").strip()
        if not det_id:
            raise PathogenError(f"{path}: mapping without determinant_id")
        if det_id in out:
            raise PathogenError(f"{path}: duplicate mapping for {det_id}")
        mode = str(entry.get("match_mode") or "").strip()
        if mode not in {"exact", "prefix", "absent"}:
            raise PathogenError(f"{det_id}: unknown match_mode {mode!r}")
        families = tuple(
            str(f) for f in (entry.get("ndaro_gene_families") or ()) if str(f).strip()
        )
        prefixes = tuple(
            str(f) for f in (entry.get("ndaro_family_prefixes") or ()) if str(f).strip()
        )
        if mode != "absent" and not (families or prefixes):
            raise PathogenError(f"{det_id}: {mode} mapping names no gene family")
        if mode == "absent" and (families or prefixes):
            raise PathogenError(f"{det_id}: absent mapping must name no family")
        if mode == "exact" and prefixes:
            raise PathogenError(
                f"{det_id}: exact mapping must not carry family prefixes"
            )
        out[det_id] = NdaroFamilyMapping(
            determinant_id=det_id,
            gene_families=families,
            match_mode=mode,
            family_prefixes=prefixes,
            note=str(entry.get("note") or "").strip(),
            requires_allele_resolution=bool(entry.get("requires_allele_resolution")),
            partial_locus=bool(entry.get("partial_locus")),
            gap_reason=(
                # The map file writes `reference_gap_reason`; accept the
                # shorter spelling too rather than silently defaulting a
                # declared reason to the generic one.
                str(
                    entry.get("reference_gap_reason")
                    or entry.get("gap_reason")
                    or "not_in_ndaro"
                )
                if mode == "absent"
                else None
            ),
        )
    return NdaroMap(
        catalog_version=str(doc.get("ndaro_catalog_version") or "unknown"),
        by_determinant=out,
        installed_family_count=int(doc.get("installed_family_count") or 0),
    )


# --------------------------------------------------------------------------- #
# Reference bundle
# --------------------------------------------------------------------------- #

#: `AMR_CDS.fa` headers are pipe-delimited:
#: `protein_accession|nucleotide_accession|?|?|allele|allele_or_symbol|product`.
#: Field 5 is the *allele*, not the gene family — `stxA2b` rather than
#: `stxA2` — so grouping on it would split every family into its variants
#: and make `stxA1` look absent when four of its alleles are installed.
#: The authoritative family lives in `ReferenceGeneCatalog.txt`, and the two
#: are joined on protein accession.
_PROTEIN_FIELD: Final = 0
_ALLELE_FIELD: Final = 4


def _family_by_protein(catalog_path: Path) -> dict[str, str]:
    """Map every reference protein accession to its NDARO gene family."""
    out: dict[str, str] = {}
    with Path(catalog_path).open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            i_family = header.index("gene_family")
            i_refseq = header.index("refseq_protein_accession")
            i_genbank = header.index("genbank_protein_accession")
        except ValueError as exc:  # pragma: no cover - malformed download
            raise PathogenError(f"{catalog_path}: unexpected columns ({exc})") from exc
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(i_family, i_refseq, i_genbank):
                continue
            family = cols[i_family].strip()
            if not family:
                continue
            for i in (i_refseq, i_genbank):
                accession = cols[i].strip()
                if accession:
                    out[accession] = family
    if not out:
        raise PathogenError(f"{catalog_path}: no gene families")
    return out


def parse_amr_cds(
    path: Path, catalog_path: Path
) -> dict[str, list[tuple[str, str, bytes]]]:
    """Group NDARO reference CDS by true gene family.

    Returns `family -> [(allele, accession, sequence)]`. Grouping by family
    rather than allele is what lets a determinant NDARO splits across
    numbered members (`ipaH` has ten, `stxA2` twenty) be screened as one
    thing while keeping each member's identity for subunit accounting.
    """
    from openbiota.pathogens.kmers import iter_fasta

    family_of = _family_by_protein(Path(catalog_path))
    families: dict[str, list[tuple[str, str, bytes]]] = {}
    unmatched = 0
    for header, seq in iter_fasta(Path(path)):
        fields = header.split()[0].split("|") if header else []
        if len(fields) <= _ALLELE_FIELD or not seq:
            continue
        accession = fields[_PROTEIN_FIELD].strip()
        family = family_of.get(accession)
        if family is None:
            # An unjoinable record is skipped rather than filed under a
            # guessed family: a misfiled reference produces a wrong gene
            # name on a report, which is worse than one fewer allele.
            unmatched += 1
            continue
        families.setdefault(family, []).append(
            (fields[_ALLELE_FIELD].strip() or family, accession, seq)
        )
    if not families:
        raise PathogenError(f"{path}: no usable NDARO CDS records")
    return families


@dataclass(frozen=True)
class DeterminantBundle:
    """Installed determinant references, resolved against the NDARO map."""

    bundle_id: str
    catalog_version: str
    fasta: Path
    index_prefix: Path | None
    #: `determinant_id -> installed gene families`, after prefix expansion.
    families_by_determinant: Mapping[str, tuple[str, ...]]
    #: `gene_family -> reference length in bases`, for coverage fractions.
    family_bases: Mapping[str, int]
    #: `gene_family -> determinant_id`, for attributing alignments back.
    determinant_of_family: Mapping[str, str]
    #: Determinants the map declared but NDARO could not supply.
    gaps: Mapping[str, str] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "bundle_id": self.bundle_id,
            "catalog_version": self.catalog_version,
            "fasta": str(self.fasta),
            "index_prefix": str(self.index_prefix) if self.index_prefix else None,
            "families_by_determinant": {
                k: list(v) for k, v in self.families_by_determinant.items()
            },
            "family_bases": dict(self.family_bases),
            "determinant_of_family": dict(self.determinant_of_family),
            "gaps": dict(self.gaps),
        }


def build_determinant_bundle(
    seeds: Sequence[DeterminantSeed],
    ndaro: NdaroMap,
    amr_cds: Path,
    gene_catalog: Path,
    out_dir: Path,
    *,
    threads: int = 8,
    bowtie2_build: str = "bowtie2-build",
    progress: Callable[[str], None] | None = None,
) -> DeterminantBundle:
    """Resolve the NDARO map against the installed CDS and index the result.

    A determinant whose families are not in the installed catalogue becomes a
    recorded gap rather than an empty search: the whole point of the map's
    `absent` mode is that "we cannot look for this" must survive into the
    report as distinct from "this is not here".
    """
    import hashlib
    import shutil
    import subprocess

    say = progress or (lambda _m: None)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    families = parse_amr_cds(Path(amr_cds), Path(gene_catalog))
    installed = set(families)
    say(f"{len(installed):,} NDARO gene families in the installed catalogue")

    families_by_determinant: dict[str, tuple[str, ...]] = {}
    determinant_of_family: dict[str, str] = {}
    gaps: dict[str, str] = {}
    for seed in seeds:
        mapping = ndaro.for_determinant(seed.determinant_id)
        if mapping is None:
            gaps[seed.determinant_id] = "no_reference_mapping"
            continue
        if not mapping.available:
            gaps[seed.determinant_id] = mapping.gap_reason or "not_in_ndaro"
            continue
        resolved = tuple(sorted(f for f in installed if mapping.matches(f)))
        if not resolved:
            # The map named families the installed catalogue does not have.
            # Recorded, not silently dropped.
            gaps[seed.determinant_id] = "mapped_families_not_installed"
            continue
        families_by_determinant[seed.determinant_id] = resolved
        for family in resolved:
            # First determinant wins; a family shared by two determinants is
            # attributed once so no fragment is double-counted.
            determinant_of_family.setdefault(family, seed.determinant_id)

    wanted = sorted(determinant_of_family)
    if not wanted:
        raise PathogenError("no determinant resolved to an installed gene family")

    fasta = out_dir / "determinants.fna"
    family_bases: dict[str, int] = {}
    with fasta.open("wb") as out:
        for family in wanted:
            total = 0
            for i, (symbol, accession, seq) in enumerate(families[family]):
                out.write(f">{family}|{symbol}|{accession}|{i}\n".encode())
                for j in range(0, len(seq), 80):
                    out.write(seq[j : j + 80] + b"\n")
                total = max(total, len(seq))
            # The *longest* member, not the sum: coverage is "how much of one
            # copy of this gene did we span", and summing alleles of the same
            # gene would make every result look like a fragment of a fragment.
            family_bases[family] = total
    say(
        f"{len(wanted):,} families for {len(families_by_determinant):,} determinants "
        f"({fasta.stat().st_size / 1e6:.1f} MB), {len(gaps):,} gaps"
    )

    index_prefix: Path | None = None
    exe = shutil.which(bowtie2_build)
    if exe is None:
        say("bowtie2-build missing: determinant screening will be unavailable")
    else:
        index_prefix = out_dir / "determinants"
        if not (
            (out_dir / "determinants.1.bt2").is_file()
            or (out_dir / "determinants.1.bt2l").is_file()
        ):
            proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
                [exe, "--threads", str(threads), "-f", str(fasta), str(index_prefix)],
                capture_output=True, text=True, check=False,
            )
            if proc.returncode != 0:
                tail = (proc.stderr or "").strip().splitlines()[-1:]
                say(f"bowtie2-build failed: {' '.join(tail)}")
                index_prefix = None
            else:
                say("determinant index built")
        else:
            say("reusing determinant index")

    fingerprint = hashlib.sha256(
        ("|".join(wanted) + ndaro.catalog_version).encode()
    ).hexdigest()[:16]
    bundle = DeterminantBundle(
        bundle_id=f"determinants-{fingerprint}",
        catalog_version=ndaro.catalog_version,
        fasta=fasta,
        index_prefix=index_prefix,
        families_by_determinant=families_by_determinant,
        family_bases=family_bases,
        determinant_of_family=determinant_of_family,
        gaps=gaps,
    )
    (out_dir / "manifest.json").write_text(
        __import__("json").dumps(bundle.to_json(), indent=1)
    )
    return bundle


# --------------------------------------------------------------------------- #
# Locus completeness
# --------------------------------------------------------------------------- #


def locus_completeness(
    seed: DeterminantSeed,
    present_subunits: Sequence[str],
    *,
    partial_reference: bool = False,
) -> str:
    """Whether every declared part of a multi-part determinant was found.

    Returns `complete`, `partial`, `single_component` or `not_assessed`.
    An operon reported as present on the strength of one subunit is a
    misleading result — LT-A without LT-B is not a functional heat-labile
    enterotoxin — so the parts are counted rather than assumed.
    """
    declared = tuple(seed.subunits or ())
    if partial_reference:
        return "partial_reference_only"
    if not declared:
        return "single_component" if present_subunits else "not_assessed"
    found = {s for s in present_subunits if s in declared}
    if not found:
        return "not_assessed"
    return "complete" if len(found) == len(declared) else "partial"


# --------------------------------------------------------------------------- #
# Host linkage
# --------------------------------------------------------------------------- #


def resolve_linkage(
    *,
    contig_shared_targets: Sequence[str] = (),
    read_pair_shared_targets: Sequence[str] = (),
    associated_target_ids: Sequence[str] = (),
) -> tuple[str, tuple[str, ...]]:
    """Decide the carrier organism from physical co-location only.

    `associated_target_ids` is the *catalog's* claim about which organisms
    normally carry this gene. It is used to order and label the answer,
    never to produce one: an association in the literature is not evidence
    of linkage in this sample. When nothing physically links the gene to a
    genome the result is `unlinked` with an empty carrier list, which the
    report renders as "carrier not established".

    Linkage is always computed when the evidence exists. The catalog's
    `requires_host_linkage` flag governs how loudly an unlinked result is
    caveated, not whether the question is asked — physical co-location is
    useful whether or not the determinant's meaning depends on it.
    """
    associated = tuple(associated_target_ids)

    def ordered(found: Sequence[str]) -> tuple[str, ...]:
        # Associated organisms first, then anything else, both sorted, so
        # the output is deterministic and the expected carrier reads first.
        known = sorted({t for t in found if t in associated})
        other = sorted({t for t in found if t not in associated})
        return tuple(known + other)

    if contig_shared_targets:
        return "contig_supported", ordered(contig_shared_targets)
    if read_pair_shared_targets:
        return "read_pair_supported", ordered(read_pair_shared_targets)
    return "unlinked", ()


# --------------------------------------------------------------------------- #
# Adjudication
# --------------------------------------------------------------------------- #


def _statement(
    seed: DeterminantSeed,
    status: str,
    linkage: str,
    carriers: Sequence[str],
    *,
    completeness: str,
    mapping: NdaroFamilyMapping | None,
) -> str:
    """Plain language that says exactly what was and was not established."""
    name = seed.display_name or seed.gene_family
    kind_word = {
        "amr": "antimicrobial-resistance gene",
        "toxin": "toxin gene",
        "pathotype_marker": "pathotype marker gene",
        "virulence_locus": "virulence locus",
    }.get(seed.kind, "gene")

    if status == "not_assessed":
        reason = (
            mapping.gap_reason.replace("_", " ")
            if mapping is not None and mapping.gap_reason
            else "no usable reference"
        )
        return (
            f"{name} was not assessed in this sample ({reason}). This is a gap "
            f"in what could be looked for, not evidence that {name} is absent."
        )

    if status == "not_detected":
        return (
            f"No sequence matching the {name} {kind_word} was found. Given the "
            f"depth of this library that argues against its presence, but a gene "
            f"carried by a rare community member can fall below detection."
        )

    if status == "candidate_homolog":
        return (
            f"Sequence resembling {name} was present but fell short of the "
            f"support needed to call the {kind_word} present. Treated as a "
            f"candidate signal, not a finding."
        )

    # Supported.
    bits = [f"Sequence for the {name} {kind_word} is present in this sample."]

    if completeness == "partial":
        bits.append(
            "Only some of its declared components were found, so the complete "
            "functional locus is not established."
        )
    elif completeness == "partial_reference_only":
        bits.append(
            "The installed reference covers only part of this locus, so "
            "completeness could not be assessed."
        )
    elif completeness == "complete" and seed.subunits:
        bits.append("All of its declared components were found.")

    if linkage in LINKED_STATES and carriers:
        where = (
            "the same assembled DNA fragment"
            if linkage == "contig_supported"
            else "the same read pair"
        )
        bits.append(
            f"It was found on {where} as {', '.join(carriers)}, which links the "
            f"gene to that organism in this sample."
        )
    elif linkage == "unlinked":
        bits.append(
            "Which organism carries it is <b>not established</b>. Stool contains "
            "many genomes, and unless a read or contig physically joins this "
            "gene to a named one, the carrier is unknown — including when a "
            "likely carrier was also detected."
        )
        if seed.linkage_absent_statement:
            bits.append(seed.linkage_absent_statement)

    if mapping is not None and mapping.requires_allele_resolution:
        bits.append(
            "This reference family contains members that differ in clinical "
            "meaning, and this result is at family level only, so the specific "
            "variant is not identified."
        )
    if seed.not_all_alleles_are:
        bits.append(seed.not_all_alleles_are)

    if seed.kind == "amr":
        bits.append(
            "Carrying a resistance gene is not the same as an infection that "
            "will fail treatment: the gene may not be expressed, and "
            "susceptibility testing on an isolate is the assay that answers "
            "that question."
        )
    elif seed.kind in {"toxin", "pathotype_marker", "virulence_locus"}:
        bits.append(
            "A toxin or virulence gene in stool is not itself disease; it "
            "describes what an organism in the community could do."
        )
    return " ".join(bits)


def adjudicate_determinant(
    seed: DeterminantSeed,
    *,
    sample_id: str,
    mapping: NdaroFamilyMapping | None,
    evidence: TargetEvidence | None = None,
    reference_bases: int | None = None,
    present_subunits: Sequence[str] = (),
    contig_shared_targets: Sequence[str] = (),
    read_pair_shared_targets: Sequence[str] = (),
    policy: DeterminantPolicy | None = None,
    assay_eligible: bool = True,
    qc_pass: bool = True,
    analysis_status: str = "completed",
) -> DeterminantResult:
    """Turn evidence for one determinant into a result record.

    Ordering matches the organism kernel: capability first, then evidence.
    A determinant with no installed reference family is `not_assessed`
    whatever the reads look like, because there was nothing to align to.
    """
    policy = policy or DeterminantPolicy()
    reasons: list[str] = []
    covered: float | None = None

    reference_status = (
        "genome_supported" if mapping is not None and mapping.available
        else "no_usable_reference"
    )

    unusable = (
        mapping is None
        or not mapping.available
        or not assay_eligible
        or not qc_pass
        or analysis_status != "completed"
    )
    if unusable:
        if mapping is None:
            reasons.append("no_reference_mapping")
        elif not mapping.available:
            reasons.append(mapping.gap_reason or "not_in_ndaro")
        if not assay_eligible:
            reasons.append("assay_ineligible")
        if not qc_pass:
            reasons.append("qc_failed")
        if analysis_status != "completed":
            reasons.append("analysis_incomplete")
        status, linkage, carriers = "not_assessed", "unlinked", ()
        completeness = "not_assessed"
    else:
        frags = evidence.qualifying_fragments if evidence else 0
        identity = (
            evidence.median_identity
            if evidence and evidence.median_identity is not None
            else 0.0
        )
        # Coverage is the fraction of the reference gene actually spanned.
        # Without a reference length it cannot be computed, and a gene call
        # made on read count alone is exactly the kind of result that turns
        # a conserved domain into a spurious resistance finding — so an
        # unknown length caps the outcome at `candidate`.
        if evidence is not None and reference_bases:
            covered = min(
                1.0, evidence.informative_bases_covered / float(reference_bases)
            )
        needed = (
            policy.min_covered_fraction_heterogeneous
            if mapping.requires_allele_resolution
            else policy.min_covered_fraction
        )
        present = (
            frags >= policy.min_distinct_fragments
            and covered is not None
            and covered >= needed
            and identity >= policy.min_identity
        )
        if frags == 0:
            status = "not_detected"
        elif not present:
            status = "candidate_homolog"
            if frags < policy.min_distinct_fragments:
                reasons.append("too_few_distinct_fragments")
            if covered is None:
                reasons.append("reference_length_unknown")
            elif covered < needed:
                reasons.append("insufficient_reference_coverage")
            if identity < policy.min_identity:
                reasons.append("identity_below_threshold")
        else:
            status = "supported_intact_sequence"

        completeness = (
            locus_completeness(
                seed, present_subunits, partial_reference=mapping.partial_locus
            )
            if present
            else "not_assessed"
        )
        if present:
            # The status carries the two caveats that change what may be
            # claimed: an incomplete locus, and a family whose members do
            # not mean the same thing clinically.
            if mapping.requires_allele_resolution:
                status = "ambiguous_allele"
            elif completeness in {"partial", "partial_reference_only"}:
                status = "supported_partial_sequence"
            linkage, carriers = resolve_linkage(
                contig_shared_targets=contig_shared_targets,
                read_pair_shared_targets=read_pair_shared_targets,
                associated_target_ids=seed.associated_seed_ids or (),
            )
        else:
            linkage, carriers = "unlinked", ()
        if mapping.requires_allele_resolution:
            reasons.append("family_level_only_allele_unresolved")
        if mapping.partial_locus:
            reasons.append("reference_covers_locus_partially")

    return DeterminantResult(
        schema_version=SCHEMA_VERSION,
        sample_id=sample_id,
        determinant_id=seed.determinant_id,
        display_name=seed.display_name or seed.gene_family,
        gene_family=seed.gene_family,
        kind=seed.kind,
        determinant_status=status,
        host_linkage=linkage,
        assay_eligibility="eligible" if assay_eligible else "ineligible",
        analysis_status=analysis_status,
        reference_status=reference_status,
        amr_class=seed.amr_class,
        aligned_identity=(
            evidence.median_identity
            if evidence and status in PRESENT_STATUSES
            else None
        ),
        covered_reference_fraction=(
            covered if status in PRESENT_STATUSES else None
        ),
        locus_completeness=completeness,
        required_discriminating_positions_assessed=(
            bool(mapping and not mapping.requires_allele_resolution)
            if status in PRESENT_STATUSES
            else False
        ),
        # Never inferred from gene presence. Different assay, different question.
        clinical_phenotype_inference="not_inferred",
        supporting_fragments=evidence.qualifying_fragments if evidence else 0,
        associated_target_ids=tuple(seed.associated_seed_ids or ()),
        linked_target_ids=tuple(carriers),
        reason_codes=tuple(dict.fromkeys(reasons)),
        plain_statement=_statement(
            seed, status, linkage, carriers,
            completeness=completeness, mapping=mapping,
        ),
    )


def screen_determinants(
    seeds: Sequence[DeterminantSeed],
    *,
    sample_id: str,
    ndaro: NdaroMap | None,
    evidence_by_family: Mapping[str, TargetEvidence] | None = None,
    subunits_by_determinant: Mapping[str, Sequence[str]] | None = None,
    family_bases: Mapping[str, int] | None = None,
    families_by_determinant: Mapping[str, Sequence[str]] | None = None,
    extra_gaps: Mapping[str, str] | None = None,
    contig_links: Mapping[str, Sequence[str]] | None = None,
    read_pair_links: Mapping[str, Sequence[str]] | None = None,
    policy: DeterminantPolicy | None = None,
    qc_pass: bool = True,
    analysis_status: str = "completed",
    progress: Callable[[str], None] | None = None,
) -> tuple[DeterminantResult, ...]:
    """Adjudicate every determinant seed for one sample.

    Every seed produces a record. A seed with no reference mapping yields
    `not_assessed`, which is what keeps an unscreenable gene family out of
    the "none found" column.
    """
    say = progress or (lambda _m: None)
    evidence_by_family = evidence_by_family or {}
    subunits_by_determinant = subunits_by_determinant or {}
    contig_links = contig_links or {}
    read_pair_links = read_pair_links or {}

    family_bases = family_bases or {}
    families_by_determinant = families_by_determinant or {}
    extra_gaps = extra_gaps or {}

    out: list[DeterminantResult] = []
    for seed in seeds:
        mapping = ndaro.for_determinant(seed.determinant_id) if ndaro else None
        gap = extra_gaps.get(seed.determinant_id)
        if mapping is not None and gap and not families_by_determinant.get(
            seed.determinant_id
        ):
            # The bundle build found the map's families were not installable.
            # Downgrade the mapping to a declared gap so the record says why
            # rather than reporting a search that never happened.
            mapping = NdaroFamilyMapping(
                determinant_id=mapping.determinant_id,
                gene_families=(),
                match_mode="absent",
                note=mapping.note,
                gap_reason=gap,
            )

        # A determinant maps to one or more families; the strongest family
        # evidence represents it, and the subunit list carries the detail.
        # The installed family list is preferred over the map's declaration,
        # because only the bundle knows what a prefix actually expanded to.
        installed = tuple(
            families_by_determinant.get(seed.determinant_id) or ()
        )
        if not installed and mapping is not None and mapping.available:
            installed = tuple(
                dict.fromkeys(
                    list(mapping.gene_families)
                    + [f for f in evidence_by_family if mapping.matches(f)]
                )
            )
        # Pick the family that is best *supported*, ranked by how much of its
        # reference is covered and only then by read count. Ranking by reads
        # alone lets a mosaic family that piles hundreds of fragments onto one
        # conserved block outrank a relative covered end to end, which would
        # downgrade a genuine group-level finding to a candidate.
        best: TargetEvidence | None = None
        best_family: str | None = None
        best_rank: tuple[float, int] = (-1.0, -1)
        for family in installed:
            ev = evidence_by_family.get(family)
            if ev is None:
                continue
            length = family_bases.get(family) or 0
            fraction = (
                min(1.0, ev.informative_bases_covered / float(length))
                if length
                else 0.0
            )
            rank = (fraction, ev.qualifying_fragments)
            if rank > best_rank:
                best, best_family, best_rank = ev, family, rank

        out.append(
            adjudicate_determinant(
                seed,
                sample_id=sample_id,
                mapping=mapping,
                evidence=best,
                reference_bases=(
                    family_bases.get(best_family) if best_family else None
                ),
                present_subunits=subunits_by_determinant.get(seed.determinant_id, ()),
                contig_shared_targets=contig_links.get(seed.determinant_id, ()),
                read_pair_shared_targets=read_pair_links.get(seed.determinant_id, ()),
                policy=policy,
                qc_pass=qc_pass,
                analysis_status=analysis_status,
            )
        )
    supported = sum(1 for d in out if d.determinant_status == "supported")
    unassessed = sum(1 for d in out if d.determinant_status == "not_assessed")
    say(
        f"determinants: {supported} supported, {unassessed} not assessed, "
        f"{len(out)} screened"
    )
    return tuple(out)


# --------------------------------------------------------------------------- #
# Read-level screening
# --------------------------------------------------------------------------- #


def _merged_length(spans: Sequence[tuple[int, int]]) -> int:
    """Total bases covered by a set of possibly overlapping intervals.

    Overlapping alignments must not each add their full length: twenty reads
    stacked on one conserved 150 bp window cover 150 bases, not 3,000, and
    counting them as 3,000 is how a single conserved domain gets promoted to
    a full-length gene call.
    """
    if not spans:
        return 0
    total = 0
    current_start, current_end = spans[0]
    for start, end in sorted(spans)[1:]:
        if start > current_end:
            total += current_end - current_start
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    return total + (current_end - current_start)


def screen_reads_for_determinants(
    r1: Path,
    r2: Path | None,
    bundle: DeterminantBundle,
    *,
    threads: int = 8,
    min_identity: float = 0.95,
    min_aligned_fraction: float = 0.70,
    binary: str = "bowtie2",
    progress: Callable[[str], None] | None = None,
) -> dict[str, TargetEvidence]:
    """Align a whole library against the determinant references, streaming.

    Direct alignment rather than k-mer nomination: the reference is a few
    megabases, so a single streamed bowtie2 pass is cheaper than building
    and querying a second k-mer index, and it needs no informative-region
    machinery — these are short genes, and the competition that matters
    (which allele within a family) is handled at adjudication.

    The SAM is consumed line by line and never held in memory. A full
    library produces tens of millions of records, almost all unaligned, so
    buffering them would cost more than the alignment itself.
    """
    import shutil
    import subprocess

    say = progress or (lambda _m: None)
    if bundle.index_prefix is None:
        raise PathogenError("determinant bundle has no aligner index")
    exe = shutil.which(binary)
    if exe is None:
        raise PathogenError(f"{binary} is required for determinant screening")

    argv = [
        exe,
        "-x", str(bundle.index_prefix),
        "--threads", str(threads),
        "--no-unal",          # unaligned records are the overwhelming majority
        "--local",            # genes are shorter than fragments; allow soft clips
        "-k", "5",            # a few allele hits per read, not all 579 families
        "--mm",
    ]
    if r2 is not None:
        argv += ["-1", str(r1), "-2", str(r2)]
    else:
        argv += ["-U", str(r1)]

    fragments: dict[str, set[str]] = {}
    spans: dict[str, list[tuple[int, int]]] = {}
    identities: dict[str, list[float]] = {}
    best_score: dict[str, int] = {}
    considered = 0

    say(f"aligning library against {len(bundle.family_bases):,} gene families")
    proc = subprocess.Popen(  # noqa: S603 - fixed argv, no shell
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    assert proc.stdout is not None
    try:
        for record in iter_sam(proc.stdout):
            considered += 1
            if record.query_length <= 0 or record.aligned_bases <= 0:
                continue
            identity = 1.0 - (
                record.edit_distance / float(max(1, record.aligned_bases))
            )
            aligned_fraction = record.aligned_bases / float(record.query_length)
            if identity < min_identity or aligned_fraction < min_aligned_fraction:
                continue
            family = record.target_id
            fragments.setdefault(family, set()).add(record.query)
            spans.setdefault(family, []).append(
                (record.pos, record.pos + record.aligned_bases)
            )
            identities.setdefault(family, []).append(identity)
            best_score[family] = max(best_score.get(family, 0), record.score)
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        stderr = proc.stderr.read().decode("utf-8", "replace") if proc.stderr else ""
        if proc.stderr is not None:
            proc.stderr.close()
        code = proc.wait()
    if code != 0:
        tail = stderr.strip().splitlines()[-1:] or ["no message"]
        raise PathogenError(f"bowtie2 failed (rc={code}): {tail[0]}")

    evidence: dict[str, TargetEvidence] = {}
    for family, names in fragments.items():
        values = sorted(identities[family])
        median = values[len(values) // 2]
        evidence[family] = TargetEvidence(
            target_id=family,
            qualifying_fragments=len(names),
            raw_fragments=len(identities[family]),
            informative_bases_covered=_merged_length(spans[family]),
            median_identity=median,
            max_score=best_score.get(family, 0),
            resolution_specific_support=True,
        )
    say(
        f"determinant alignment: {considered:,} records considered, "
        f"{len(evidence):,} gene families with qualifying support"
    )
    return evidence


def screen_determinants_from_reads(
    seeds: Sequence[DeterminantSeed],
    *,
    sample_id: str,
    r1: Path | None,
    r2: Path | None,
    ndaro: NdaroMap | None,
    bundle: DeterminantBundle | None,
    policy: DeterminantPolicy | None = None,
    qc_pass: bool = True,
    threads: int = 8,
    progress: Callable[[str], None] | None = None,
) -> tuple[DeterminantResult, ...]:
    """Full determinant screen for one sample, reads to adjudicated records.

    Every failure mode yields `not_assessed` records rather than an empty
    or negative panel: no reads, no bundle, no aligner, or an alignment that
    did not complete.
    """
    say = progress or (lambda _m: None)
    if not seeds:
        return ()
    if bundle is None or ndaro is None or r1 is None or not Path(r1).is_file():
        say("determinant screening not run: reads or reference unavailable")
        return screen_determinants(
            seeds, sample_id=sample_id, ndaro=ndaro,
            analysis_status="not_run", progress=say,
        )

    try:
        evidence = screen_reads_for_determinants(
            Path(r1), Path(r2) if r2 else None, bundle,
            threads=threads, progress=say,
        )
    except PathogenError as exc:
        say(f"determinant screening failed: {exc}")
        return screen_determinants(
            seeds, sample_id=sample_id, ndaro=ndaro,
            analysis_status="failed", progress=say,
        )

    # The subunits actually seen, so operon completeness is counted rather
    # than assumed from the family having any support at all.
    subunits: dict[str, list[str]] = {}
    for det_id, families in bundle.families_by_determinant.items():
        seen = [f for f in families if f in evidence]
        if seen:
            subunits[det_id] = seen

    return screen_determinants(
        seeds,
        sample_id=sample_id,
        ndaro=ndaro,
        evidence_by_family=evidence,
        subunits_by_determinant=subunits,
        family_bases=bundle.family_bases,
        families_by_determinant=bundle.families_by_determinant,
        extra_gaps=bundle.gaps,
        policy=policy,
        qc_pass=qc_pass,
        progress=say,
    )


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #


def _seed(**kw: Any) -> DeterminantSeed:
    base: dict[str, Any] = {
        "determinant_id": "amr.test",
        "source_section": "8.3",
        "display_name": "testGene",
        "gene_family": "testGene",
        "kind": "amr",
        "candidate_sources": ("NDARO",),
        "note": "",
    }
    base.update(kw)
    return DeterminantSeed(**base)


#: Reference length used throughout the self-test, so coverage fractions in
#: the cases below are exact rather than approximate.
_GENE_BASES: Final = 1_000


def _ev(frags: int, covered: float, identity: float) -> TargetEvidence:
    return TargetEvidence(
        target_id="t",
        qualifying_fragments=frags,
        informative_bases_covered=int(round(covered * _GENE_BASES)),
        median_identity=identity,
    )


def self_test() -> int:  # noqa: PLR0915 - a contract suite reads best flat
    """Invariants that protect the gene/carrier separation."""
    checks = 0
    mapping = NdaroFamilyMapping("amr.test", ("testGene",), "exact")
    strong = _ev(6, 0.95, 0.99)

    # A gene with no reference mapping is not assessed, never absent.
    r = adjudicate_determinant(_seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=None, evidence=strong)
    assert r.determinant_status == "not_assessed"
    assert "not evidence that" in r.plain_statement
    checks += 1

    # An explicit NDARO gap behaves the same way.
    absent = NdaroFamilyMapping("amr.test", (), "absent", gap_reason="not_in_ndaro")
    r = adjudicate_determinant(
        _seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=absent, evidence=strong
    )
    assert r.determinant_status == "not_assessed"
    assert "not_in_ndaro" in r.reason_codes
    checks += 1

    # No reads at all is a negative, and it is worded as one.
    r = adjudicate_determinant(_seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=None)
    assert r.determinant_status == "not_detected"
    assert r.host_linkage == "unlinked"
    checks += 1

    # Support requires fragments, coverage *and* identity.
    for weak in (_ev(1, 0.95, 0.99), _ev(6, 0.20, 0.99), _ev(6, 0.95, 0.80)):
        r = adjudicate_determinant(
            _seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=weak
        )
        assert r.determinant_status == "candidate_homolog", weak
        checks += 1

    # Supported but unlinked: the carrier must stay unnamed, and the
    # statement must say so in words.
    r = adjudicate_determinant(
        _seed(associated_seed_ids=("bacteria.kpneumoniae",)),
        sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=strong,
    )
    assert r.determinant_status in PRESENT_STATUSES
    assert r.host_linkage == "unlinked"
    assert r.linked_target_ids == ()
    assert "not established" in r.plain_statement
    # The association is preserved as metadata but is not a carrier claim.
    assert r.associated_target_ids == ("bacteria.kpneumoniae",)
    checks += 1

    # Contig linkage outranks read-pair linkage and names the carrier.
    r = adjudicate_determinant(
        _seed(associated_seed_ids=("bacteria.kpneumoniae",)),
        sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=strong,
        contig_shared_targets=("bacteria.kpneumoniae",),
        read_pair_shared_targets=("bacteria.other",),
    )
    assert r.host_linkage == "contig_supported"
    assert r.linked_target_ids == ("bacteria.kpneumoniae",)
    checks += 1

    # Phenotype is never inferred from presence.
    assert r.clinical_phenotype_inference == "not_inferred"
    assert "not the same as an infection" in r.plain_statement
    checks += 1

    # Operon completeness: one subunit of two is partial, not present.
    two = _seed(determinant_id="toxin.lt", kind="toxin", subunits=("ltA", "ltB"))
    r = adjudicate_determinant(
        two, sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=strong,
        present_subunits=("ltA",),
    )
    assert r.locus_completeness == "partial"
    assert "not established" in r.plain_statement
    r = adjudicate_determinant(
        two, sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=strong,
        present_subunits=("ltA", "ltB"),
    )
    assert r.locus_completeness == "complete"
    checks += 2

    # A heterogeneous family needs more coverage and says it is family-level.
    hetero = NdaroFamilyMapping(
        "amr.test", ("blaOXA",), "prefix", requires_allele_resolution=True
    )
    r = adjudicate_determinant(
        _seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=hetero, evidence=_ev(6, 0.70, 0.99)
    )
    assert r.determinant_status == "candidate_homolog"
    r = adjudicate_determinant(
        _seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=hetero, evidence=_ev(6, 0.95, 0.99)
    )
    assert r.determinant_status in PRESENT_STATUSES
    assert not r.required_discriminating_positions_assessed
    assert "family level only" in r.plain_statement
    checks += 2

    # QC failure and ineligibility are gaps, not negatives.
    for kw in ({"qc_pass": False}, {"assay_eligible": False},
               {"analysis_status": "failed"}):
        r = adjudicate_determinant(
            _seed(), sample_id="S", reference_bases=_GENE_BASES, mapping=mapping, evidence=strong, **kw
        )
        assert r.determinant_status == "not_assessed"
        checks += 1

    # Every seed yields exactly one record, including unmappable ones.
    seeds = (_seed(), _seed(determinant_id="amr.other", gene_family="other"))
    rows = screen_determinants(seeds, sample_id="S", ndaro=None)
    assert len(rows) == len(seeds)
    assert {r.determinant_status for r in rows} == {"not_assessed"}
    checks += 1

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"determinants self-test: {self_test()} checks passed")
