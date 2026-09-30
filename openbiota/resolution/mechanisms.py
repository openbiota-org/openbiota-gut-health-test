"""Measured candidate mechanisms, replacing the disease-name causal bridge.

The rule this module retires worked like this: if a donor sat at or above the
75th percentile of some disease profile, at least 20 points above the
recipient, and any animal-transfer paper existed for that disease, the donor
was flagged as introducing a transferable causal pattern.

That rule is wrong, and spec §8.4 requires its removal. A percentile is a
resemblance to a case-control community pattern. It is not a feature, so it
cannot be transferred, and a paper about *the disease* says nothing about
whether *this donor* carries the thing the paper studied. In the supplied
samples the rule flagged colorectal cancer for SAMPLE1 on a resemblance score while
SAMPLE1's actual measured colibactin evidence was a single `clbB` fragment — which
is not an intact biosynthetic island by any standard.

What replaces it is narrower and harder: **a candidate mechanism counts only
when the sequence evidence for it was measured in that donor.** The percentile
stays in the report as descriptive resemblance, clearly labelled, and never
blocks a donor by itself (§10).

The rules below are all versions of one principle — a gene-family fragment
count is the weakest possible evidence, and the report must say so:

* One or two `clbB` fragments are not a functional pks island. The island is
  ~55 kb across clbA-S; a single NRPS-PKS fragment is a candidate signal and
  nothing more (acceptance 31).
* A `bft` fragment count is not enterotoxigenic *B. fragilis*. The subtype
  (bft1/2/3), the gene's integrity and the carrier all remain unresolved
  (acceptance 30).
* Nothing here is confirmed, and a zero count with no calibrated limit is
  `not_detected_limit_unvalidated`, never a clean negative.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final

from .schema import (
    AMBIGUOUS_HOMOLOGY,
    ANIMAL_MECHANISTIC,
    CANDIDATE_SEQUENCE,
    CHARACTERIZED_FUNCTION,
    COMPLETED,
    HUMAN_ASSOCIATION,
    NOT_DETECTED_UNVALIDATED,
    BiologicalEvidence,
    Coverage,
    Linkage,
    Provenance,
    ResolutionCall,
    Validation,
)


@dataclass(frozen=True, slots=True)
class MechanismGene:
    """One measured gene family and what it can and cannot establish."""

    target_id: str
    gene: str
    label: str
    #: The complete locus this gene belongs to, when it belongs to one.
    locus: str
    #: Roughly how many genes the functional locus needs. A single-gene hit
    #: against a multi-gene island is the clearest case of over-reading.
    locus_gene_count: int | None
    required_for_call: str
    evidence: BiologicalEvidence
    #: Subtypes that would have to be distinguished for a real call.
    subtypes: tuple[str, ...] = ()


PKS_EVIDENCE = BiologicalEvidence(
    kind=CHARACTERIZED_FUNCTION,
    statement=(
        "The colibactin pks island produces a genotoxin that crosslinks DNA and leaves a "
        "distinctive mutational signature in human colorectal tumours."
    ),
    studied_entity="55,140-bp pks island of E. coli IHE3034 (AM229678.1); MIBiG BGC0000972.5",
    endpoint="DNA damage and mutational signature",
    sources=(
        "AM229678.1 (IHE3034 pks island)",
        "MIBiG BGC0000972.5",
    ),
    boundary=(
        "Genotoxic capacity requires the intact multi-gene island. A hit on one "
        "NRPS-PKS gene establishes neither island architecture nor a functional pathway, "
        "and the carrier may be E. coli, Klebsiella, Citrobacter or another organism."
    ),
)

BFT_EVIDENCE = BiologicalEvidence(
    kind=ANIMAL_MECHANISTIC,
    statement=(
        "B. fragilis toxin (fragilysin) cleaves E-cadherin and drives IL-17-associated "
        "mucosal inflammation in experimental colitis models."
    ),
    studied_entity="bft1 (AB026625.1), bft2 (AB026626.1), bft3 (AB026624.1) alleles",
    endpoint="E-cadherin cleavage and IL-17 colitis in animal models",
    sources=(
        "AB026625.1 / AB026626.1 / AB026624.1",
        "Allele study PMID 10612750",
    ),
    boundary=(
        "The species alone cannot distinguish enterotoxigenic from nontoxigenic "
        "populations, and related metalloprotease homologs are decoys. A fragment count "
        "does not establish the subtype, the gene's integrity or which organism carries it."
    ),
)

FADA_EVIDENCE = BiologicalEvidence(
    kind=HUMAN_ASSOCIATION,
    statement=(
        "FadA is a Fusobacterium adhesin that binds E-cadherin and activates "
        "beta-catenin signalling."
    ),
    studied_entity="F. nucleatum ATCC25586 FN0264",
    endpoint="adhesion and beta-catenin activation",
    sources=("FN0264, F. nucleatum ATCC25586",),
    boundary=(
        "Close paralogs FN1529 and FNV2159 are decoys that must be distinguished. "
        "Adhesin DNA is not a cancer diagnosis."
    ),
)

#: The mechanisms the existing functional panel actually measures. Each is
#: registered as its own target so the report can never again present a
#: combined panel title as though the whole panel were positive.
GENES: Final[tuple[MechanismGene, ...]] = (
    MechanismGene(
        target_id="mech.pks_island.clbB",
        gene="clbB",
        label="Colibactin pks island (clbB fragment evidence)",
        locus="clbA-S pks island",
        locus_gene_count=19,
        required_for_call=(
            "full-locus gene architecture across clbA-S with integrity and carrier evidence"
        ),
        evidence=PKS_EVIDENCE,
    ),
    MechanismGene(
        target_id="mech.etbf.bft",
        gene="bft",
        label="B. fragilis toxin (bft fragment evidence)",
        locus="bft within the BfPAI",
        locus_gene_count=None,
        required_for_call=(
            "subtype-discriminating CDS mapping (bft1/2/3), intact-versus-disrupted state "
            "and BfPAI context with carrier linkage"
        ),
        evidence=BFT_EVIDENCE,
        subtypes=("bft1", "bft2", "bft3"),
    ),
    MechanismGene(
        target_id="mech.fusobacterium.fadA",
        gene="fadA",
        label="Fusobacterium FadA adhesin (fadA fragment evidence)",
        locus="fadA",
        locus_gene_count=1,
        required_for_call=(
            "accession-versioned CDS match distinguished from paralogs FN1529 and FNV2159"
        ),
        evidence=FADA_EVIDENCE,
    ),
)

BY_GENE: Final[Mapping[str, MechanismGene]] = {g.gene: g for g in GENES}

#: Below this many independent fragments, a multi-gene locus cannot even be
#: argued to be partially present. Not a clinical threshold: a reporting floor.
SINGLE_GENE_FRAGMENT_FLOOR: Final = 3

#: Which nucleotide locus check corroborates or contradicts each gene panel.
#: The panel searches translated protein; the locus stage maps nucleotides
#: against exact accession-versioned references. When a protein hit at high
#: amino-acid identity has no nucleotide counterpart, the hit is to a homolog
#: rather than the target -- which is the decoy problem the panel cannot see.
CORROBORATING_LOCUS: Final[Mapping[str, str]] = {
    "mech.etbf.bft": "locus.bft.subtype",
    "mech.pks_island.clbB": "locus.pks_island.architecture",
}


def _locus_verdict(
    results: Mapping[str, Any], locus_target_id: str
) -> Mapping[str, Any] | None:
    """The nucleotide check for a gene panel target, if it ran."""
    for row in ((results.get("locus_resolution") or {}).get("targets") or []):
        if row.get("target_id") == locus_target_id:
            return row
    return None


def _panel_genes(results: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    """The crc_virulence panel's gene rows, by gene name."""
    for panel in results.get("panels") or []:
        if isinstance(panel, Mapping) and panel.get("name") == "crc_virulence":
            return {
                str(row.get("gene")): row
                for row in panel.get("genes") or []
                if isinstance(row, Mapping)
            }
    return {}


def calls_for_sample(
    sample_id: str, results: Mapping[str, Any], *, database_release: str | None = None
) -> list[ResolutionCall]:
    """Typed mechanism calls for one sample, from measured panel evidence."""
    rows = _panel_genes(results)
    out: list[ResolutionCall] = []
    for gene in GENES:
        row = rows.get(gene.gene)
        if row is None:
            # The panel did not run or does not contain this gene. Not a
            # negative: the target is simply unassessed.
            out.append(
                ResolutionCall(
                    sample_id=sample_id,
                    target_id=gene.target_id,
                    identity_kind="gene_family",
                    assay_status="not_requested",
                    reason_codes=("panel_row_absent",),
                    biological_evidence=(gene.evidence,),
                    plain=(
                        f"{gene.label}: the functional panel did not report this gene, so it "
                        "was not assessed. That is not an absence."
                    ),
                )
            )
            continue

        fragments = int(row.get("fragments") or 0)
        identity = row.get("median_identity")
        residue = row.get("residue_check")
        coverage = Coverage(
            independent_fragments=fragments,
            unique_fragments=fragments,
            # No calibrated limit of detection exists for these panels, so the
            # fields stay null rather than implying one.
            limit_of_detection=None,
            limit_scope=None,
        )
        provenance = Provenance(
            database_release=database_release,
            reference_accession_versions=tuple(gene.evidence.sources),
            parameters={
                "engine": "diamond",
                "panel": "crc_virulence",
                "median_identity": identity,
                "residue_check": residue,
                "reference_count": row.get("reference_count"),
            },
            threshold_version="crc_virulence/1",
        )

        if fragments <= 0:
            out.append(
                ResolutionCall(
                    sample_id=sample_id,
                    target_id=gene.target_id,
                    identity_kind="gene_family",
                    assay_status=COMPLETED,
                    call_id=f"{sample_id}:{gene.target_id}",
                    analytical_call=NOT_DETECTED_UNVALIDATED,
                    reason_codes=("no_fragments", "detection_limit_not_calibrated"),
                    biological_evidence=(gene.evidence,),
                    coverage=coverage,
                    validation=Validation(
                        analytical_status="not_validated",
                        reference_function_status="experimentally_characterized_reference",
                        phenotype_measured_in_sample=False,
                        clinical_predictive_status="not_established",
                    ),
                    provenance=provenance,
                    plain=(
                        f"No {gene.gene} fragments were found. This panel has no calibrated "
                        "detection limit, so that is a non-detection rather than a "
                        "demonstrated absence."
                    ),
                )
            )
            continue

        # Present as fragments. The question is what that can support.
        multi_gene = (gene.locus_gene_count or 1) > 1
        thin = fragments < SINGLE_GENE_FRAGMENT_FLOOR
        reasons = ["gene_family_fragments_only", "carrier_unresolved"]

        # Reconcile against the nucleotide check rather than reporting two
        # unrelated numbers.
        locus = _locus_verdict(results, CORROBORATING_LOCUS.get(gene.target_id, ""))
        contradicted = bool(
            locus
            and locus.get("assay_status") == COMPLETED
            and str(locus.get("analytical_call", "")).startswith("not_detected")
        )
        corroborated = bool(
            locus
            and locus.get("assay_status") == COMPLETED
            and str(locus.get("analytical_call", "")) in {
                CANDIDATE_SEQUENCE, "supported_detection"
            }
        )
        if contradicted:
            reasons.append("nucleotide_check_found_no_match")
        elif corroborated:
            reasons.append("nucleotide_check_corroborates")
        if multi_gene:
            reasons.append("locus_architecture_not_assessed")
        if gene.subtypes:
            reasons.append("subtype_not_discriminated")
        if thin:
            reasons.append("fragment_count_below_reporting_floor")

        if multi_gene:
            plain = (
                f"{fragments} {gene.gene} fragment(s) matched. The functional unit is the "
                f"{gene.locus} (~{gene.locus_gene_count} genes); one gene family's fragments "
                "cannot establish that the island is present or intact, and the carrier "
                "organism is unresolved. This is a candidate signal, not a functional "
                "capacity."
            )
        else:
            plain = (
                f"{fragments} {gene.gene} fragment(s) matched"
                + (
                    f", subtype not distinguished among {', '.join(gene.subtypes)}"
                    if gene.subtypes else ""
                )
                + ". Gene integrity and the carrier organism are unresolved, so this is a "
                "candidate signal rather than a confirmed pathotype."
            )

        if contradicted:
            plain += (
                " <b>The nucleotide check disagrees.</b> Reads were mapped against the exact "
                "reference sequence(s) for this target and none aligned, so the protein-level "
                "match is to a homolog in the same family rather than to this target itself. "
                "A high amino-acid identity to a gene-family reference does not survive that "
                "test, and this is the decoy case the translated search cannot distinguish."
            )
        elif corroborated:
            plain += (
                " The nucleotide check agrees: reads also map to the exact reference sequence, "
                "which the protein search alone could not establish."
            )

        out.append(
            ResolutionCall(
                sample_id=sample_id,
                target_id=gene.target_id,
                identity_kind="gene_family",
                assay_status=COMPLETED,
                call_id=f"{sample_id}:{gene.target_id}",
                # A contradicted hit is an ambiguous homology result, not a
                # candidate detection of this target.
                analytical_call=(
                    AMBIGUOUS_HOMOLOGY if contradicted else CANDIDATE_SEQUENCE
                ),
                reason_codes=tuple(reasons),
                biological_evidence=(gene.evidence,),
                coverage=coverage,
                # Sample co-detection is the weakest rung of the linkage
                # hierarchy and the only one a translated protein search can
                # reach (spec §6.1).
                linkage=Linkage(state="sample_co_detection"),
                validation=Validation(
                    analytical_status="not_validated",
                    reference_function_status="experimentally_characterized_reference",
                    phenotype_measured_in_sample=False,
                    clinical_predictive_status="not_established",
                ),
                provenance=provenance,
                plain=plain,
            )
        )
    return out


@dataclass(frozen=True, slots=True)
class MechanismFinding:
    """A donor's candidate mechanism, for the matcher to display per donor."""

    person_id: str
    target_id: str
    label: str
    gene: str
    fragments: int
    analytical_call: str
    required_for_call: str
    reason_codes: tuple[str, ...]
    plain: str
    sources: tuple[str, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "person_id": self.person_id,
            "target_id": self.target_id,
            "label": self.label,
            "gene": self.gene,
            "fragments": self.fragments,
            "analytical_call": self.analytical_call,
            "resolved_as_homolog": self.analytical_call == AMBIGUOUS_HOMOLOGY,
            "required_for_a_confirmed_call": self.required_for_call,
            "reason_codes": list(self.reason_codes),
            "plain": self.plain,
            "sources": list(self.sources),
            # None of these is ever a confirmed hazard, so none of them may
            # produce a probability or a numeric risk (spec §9.1).
            "confirmed": False,
            "disease_transmission_probability": None,
        }


def findings_for(
    person_id: str, results: Mapping[str, Any]
) -> list[MechanismFinding]:
    """Candidate mechanisms measured in one person's material.

    Reported per person and never pooled: one donor's clean panel does not
    offset another's candidate signal (spec §10, acceptance 61).
    """
    out: list[MechanismFinding] = []
    for call in calls_for_sample(person_id, results):
        # Ambiguous-homology rows stay in the ledger. "We mapped the exact
        # reference and nothing aligned, so the protein hit is a homolog" is a
        # resolved question and useful to a reader; dropping it would leave the
        # panel's raw fragment count as the last word.
        if not call.is_positive and call.analytical_call != AMBIGUOUS_HOMOLOGY:
            continue
        gene = next(g for g in GENES if g.target_id == call.target_id)
        out.append(
            MechanismFinding(
                person_id=person_id,
                target_id=call.target_id,
                label=gene.label,
                gene=gene.gene,
                fragments=call.coverage.independent_fragments or 0,
                analytical_call=call.analytical_call,
                required_for_call=gene.required_for_call,
                reason_codes=call.reason_codes,
                plain=call.plain,
                sources=gene.evidence.sources,
            )
        )
    return sorted(out, key=lambda f: (-f.fragments, f.target_id))


def self_test() -> int:
    """Check the over-reading guards. Returns a failure count."""
    failures = 0

    def panel(**counts: int) -> dict[str, Any]:
        return {
            "panels": [
                {
                    "name": "crc_virulence",
                    "genes": [
                        {"gene": g, "fragments": n, "median_identity": 95.0,
                         "residue_check": "pass", "reference_count": 10}
                        for g, n in counts.items()
                    ],
                }
            ]
        }

    # Acceptance 31: one or two clbB fragments are not an intact pks island.
    for n in (1, 2):
        calls = {c.target_id: c for c in calls_for_sample("S", panel(clbB=n))}
        pks = calls["mech.pks_island.clbB"]
        if pks.analytical_call != CANDIDATE_SEQUENCE:
            failures += 1
        if "locus_architecture_not_assessed" not in pks.reason_codes:
            failures += 1
        # The wording must deny the island, not assert it.
        if "cannot establish" not in pks.plain:
            failures += 1
        if "candidate signal" not in pks.plain:
            failures += 1

    # Acceptance 30: a bft count stays target-specific and unconfirmed.
    calls = {c.target_id: c for c in calls_for_sample("SAMPLE4", panel(bft=11))}
    bft = calls["mech.etbf.bft"]
    if bft.analytical_call != CANDIDATE_SEQUENCE:
        failures += 1
    if "subtype_not_discriminated" not in bft.reason_codes:
        failures += 1
    if bft.linkage.state != "sample_co_detection":
        failures += 1

    # A zero count is never a calibrated negative.
    calls = {c.target_id: c for c in calls_for_sample("S", panel(fadA=0))}
    fada = calls["mech.fusobacterium.fadA"]
    if fada.analytical_call != NOT_DETECTED_UNVALIDATED:
        failures += 1
    if fada.coverage.limit_of_detection is not None:
        failures += 1

    # A missing panel row is unassessed, not absent.
    calls = {c.target_id: c for c in calls_for_sample("S", {"panels": []})}
    if calls["mech.etbf.bft"].assay_status != "not_requested":
        failures += 1
    if calls["mech.etbf.bft"].is_negative:
        failures += 1

    # A contradicted panel hit becomes ambiguous homology, stays visible, and
    # never reads as a candidate mechanism.
    contradicted = dict(
        panel(bft=11),
        locus_resolution={"targets": [{
            "target_id": "locus.bft.subtype",
            "assay_status": COMPLETED,
            "analytical_call": NOT_DETECTED_UNVALIDATED,
        }]},
    )
    calls = {c.target_id: c for c in calls_for_sample("SAMPLE4", contradicted)}
    bft_call = calls["mech.etbf.bft"]
    if bft_call.analytical_call != AMBIGUOUS_HOMOLOGY:
        failures += 1
    if "nucleotide_check_found_no_match" not in bft_call.reason_codes:
        failures += 1
    if "match is to a homolog" not in bft_call.plain:
        failures += 1
    visible = [f for f in findings_for("SAMPLE4", contradicted) if f.gene == "bft"]
    if not visible or visible[0].to_json()["resolved_as_homolog"] is not True:
        failures += 1

    # Findings never carry a probability.
    for finding in findings_for("SAMPLE4", panel(bft=11, clbB=1)):
        blob = finding.to_json()
        if blob["confirmed"] is not False:
            failures += 1
        if blob["disease_transmission_probability"] is not None:
            failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
