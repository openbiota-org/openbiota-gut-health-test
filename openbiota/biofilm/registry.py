"""The 24 biological modules (BF-M01 .. BF-M24) and their scoring policy.

Each module names a mechanism, the features that would detect it, the axis
it may contribute to, and - critically - the conditions under which it may
contribute at all. A module that is merely interesting does not score.

Three separate status fields exist because they fail independently:

``analytical_status``
    Can we identify the sequence reliably at the stated resolution? This is
    a molecular question. It can be ``supported_at_stated_resolution``
    while the biology is still unclear.

``biological_evidence_status``
    What kind of study supports the claimed direction? Inherited from the
    module's sources and capped by them.

``clinical_validation_status``
    Is there a validated relationship to an outcome? For every module here
    the answer is ``not_established``, and the spec is explicit that this
    must not suppress an analytically supported research measurement.

A module contributes to an axis only when :meth:`Module.eligible_for_axis`
holds: the analytics must be supported, the direction must come from a
reviewed claim rather than a correlation, and the module must not be
``context_dependent``. Generic adhesion, capsule, quorum sensing and
"a beneficial species is present" are all explicitly refused a direction.

Dependence groups prevent double counting. BF-M03 and BF-M04 both measure
PNAG and therefore share ``BF-G-PNAG``: a long operon, many homologs, or
several papers about one mechanism cannot outvote a different mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from . import sources as src

#: Census states. Every module is in exactly one, and the report prints the
#: distribution so nothing is silently omitted.
CENSUS_STATES: Final = (
    "implemented",
    "candidate_not_validated",
    "context_only",
    "source_unavailable",
    "not_measurable_from_current_assay",
)

ANALYTICAL_STATES: Final = (
    "pending_validation",
    "supported_at_stated_resolution",
    "failed_validation",
)

CLINICAL_STATES: Final = ("not_established", "established")

#: Which axis a module could inform. ``context_dependent`` means the
#: mechanism has no agreed direction and must never be added to either axis.
AXES: Final = ("harmful_associated", "protective_associated", "context_dependent")

#: What kind of thing is being measured.
MEASUREMENT_KINDS: Final = (
    "gene_family",
    "locus_completeness",
    "community_gene_pattern",
    "taxonomic_pattern",
    "external_assay",
    "cross_reference",
)


@dataclass(frozen=True, slots=True)
class Module:
    """One biological module in the census."""

    id: str
    label: str
    measurement_kind: str
    mechanism_group: str
    primary_source_ids: tuple[str, ...]
    feature_symbols: tuple[str, ...]
    census_state: str
    analytical_status: str
    candidate_axis: str
    #: Dependence group. Modules sharing one cannot double-count.
    aggregation_group: str
    interpretation: str
    #: What this module explicitly does not establish.
    not_established: str
    clinical_validation_status: str = "not_established"
    #: True when naming an organism requires resolved carrier linkage.
    carrier_required_for_organism_claim: bool = True
    molecular_resolution: str = "gene_family"
    #: An explicit, source-backed adverse or protective context claim.
    #: Without one a module cannot take a direction, however many genes hit.
    directional_claim: str = ""
    report_even_if_unscored: bool = True
    notes: str = ""

    def __post_init__(self) -> None:
        if self.census_state not in CENSUS_STATES:
            raise ValueError(f"{self.id}: bad census_state {self.census_state!r}")
        if self.analytical_status not in ANALYTICAL_STATES:
            raise ValueError(f"{self.id}: bad analytical_status {self.analytical_status!r}")
        if self.candidate_axis not in AXES:
            raise ValueError(f"{self.id}: bad candidate_axis {self.candidate_axis!r}")
        if self.measurement_kind not in MEASUREMENT_KINDS:
            raise ValueError(f"{self.id}: bad measurement_kind {self.measurement_kind!r}")
        if self.clinical_validation_status not in CLINICAL_STATES:
            raise ValueError(f"{self.id}: bad clinical status")
        # A module may not assert something its own sources forbid.
        for sid in self.primary_source_ids:
            src.get(sid)
        # A directional module must say why, in words traceable to a source.
        if self.candidate_axis != "context_dependent" and not self.directional_claim:
            raise ValueError(
                f"{self.id}: axis {self.candidate_axis} requires a directional_claim; "
                "generic adhesion, capsule or quorum positivity cannot take a "
                "direction without an explicit source-backed context claim"
            )

    @property
    def biological_evidence_status(self) -> str:
        """Inherited from sources, never declared independently."""
        return src.best_evidence_class(self.primary_source_ids)

    @property
    def eligible_for_axis(self) -> bool:
        """Whether this module may contribute a percentile to an axis.

        Requires supported analytics, a real direction, and an explicit
        source-backed directional claim. Clinical validation is deliberately
        NOT required - section 8.4 is explicit that a supported gene-family
        assay plus a source-specific biological claim may be activated for
        experimental reference ranking without a human outcome trial.
        """
        return (
            self.analytical_status == "supported_at_stated_resolution"
            and self.candidate_axis in {"harmful_associated", "protective_associated"}
            and bool(self.directional_claim)
            and self.census_state == "implemented"
        )

    @property
    def ineligibility_reason(self) -> str:
        """Why this module is not contributing, recorded per section 8.3."""
        if self.eligible_for_axis:
            return ""
        if self.candidate_axis == "context_dependent":
            return "context_dependent_no_agreed_direction"
        if self.census_state != "implemented":
            return self.census_state
        if self.analytical_status != "supported_at_stated_resolution":
            return self.analytical_status
        return "no_directional_claim"

    def manifest(self) -> dict[str, object]:
        """The module manifest contract from section 6."""
        return {
            "schema_version": "openbiota.biofilm-module/1.1",
            "module_id": self.id,
            "label": self.label,
            "measurement_kind": self.measurement_kind,
            "mechanism_group": self.mechanism_group,
            "primary_source_ids": list(self.primary_source_ids),
            "feature_symbols": list(self.feature_symbols),
            "sequence_mapping_status": (
                "resolved"
                if self.analytical_status == "supported_at_stated_resolution"
                else "requires_source_identifier_resolution"
            ),
            "analytical_status": self.analytical_status,
            "biological_evidence_status": self.biological_evidence_status,
            "clinical_validation_status": self.clinical_validation_status,
            "molecular_resolution": self.molecular_resolution,
            "carrier_required_for_organism_claim": self.carrier_required_for_organism_claim,
            "census_state": self.census_state,
            "interpretation": {
                "default_valence": self.candidate_axis,
                "directional_claim": self.directional_claim,
                "phenotype_directly_measured": False,
                "summary": self.interpretation,
                "not_established": self.not_established,
            },
            "score": {
                "candidate_axis": self.candidate_axis,
                "eligible_for_axis": self.eligible_for_axis,
                "ineligibility_reason": self.ineligibility_reason,
                "aggregation_group": self.aggregation_group,
            },
            "report_even_if_unscored": self.report_even_if_unscored,
            "notes": self.notes,
        }


_M: Final = (
    Module(
        id="BF-M01",
        label="Curli-family bacterial amyloid genes",
        measurement_kind="locus_completeness",
        mechanism_group="amyloid",
        primary_source_ids=("BF-S06", "BF-S07"),
        feature_symbols=("csgA", "csgB", "csgD", "csgE", "csgF", "csgG"),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-AMYLOID",
        interpretation=(
            "Bacterial amyloid production potential. A complete supported locus "
            "means more than an isolated regulator, and healthy carriage is "
            "common: curli genes are widespread in commensal Enterobacteriaceae."
        ),
        not_established=(
            "Host fibrin, brain amyloid, Parkinson's disease, or an active "
            "biofilm. Harmful-axis eligibility would need an explicit "
            "experimentally supported adverse-context claim, which csgA "
            "detection alone does not provide."
        ),
        notes=(
            "Deliberately context_dependent. The spec is explicit that csgA "
            "detection is not by itself an adverse finding."
        ),
    ),
    Module(
        id="BF-M02",
        label="Biofilm-associated protein (BAP) family targets",
        measurement_kind="gene_family",
        mechanism_group="amyloid",
        primary_source_ids=("BF-S06",),
        feature_symbols=("bap", "bapA", "esp", "sasG"),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-AMYLOID",
        interpretation=(
            "Amyloid-associated adhesive potential from the BF-S06 Table S2 "
            "reference identifiers. Validated domains, untested homologs, "
            "organism and study model are kept apart."
        ),
        not_established=(
            "That an untested homolog behaves like the experimentally evaluated "
            "reference domain. Repetitive surface proteins are especially prone "
            "to ambiguous domain hits."
        ),
        notes=(
            "Must come from Table S2 identifiers, not from every annotation "
            "containing the string 'bap'."
        ),
    ),
    Module(
        id="BF-M03",
        label="PNAG/PGA matrix synthesis and export locus",
        measurement_kind="locus_completeness",
        mechanism_group="pnag_and_regulation",
        primary_source_ids=("BF-S05", "BF-S07"),
        feature_symbols=("pgaA", "pgaB", "pgaC", "pgaD", "icaA", "icaB", "icaC", "icaD"),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PNAG",
        interpretation=(
            "Matrix-production potential. Context-dependent on purpose: PNAG "
            "synthesis is widespread and its health direction is not agreed."
        ),
        not_established=(
            "Universal harm, or that pga and ica systems are interchangeable."
        ),
    ),
    Module(
        id="BF-M04",
        label="Bile-acid-diarrhoea-associated biofilm gene pattern",
        measurement_kind="community_gene_pattern",
        mechanism_group="pnag_and_regulation",
        primary_source_ids=("BF-S05",),
        feature_symbols=("bssS", "pgaA", "pgaB", "pgaC", "pgaD"),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PNAG",
        interpretation=(
            "One exploratory human disease-associated pattern, scored as a single "
            "grouped signal rather than six independent proofs. Shares "
            "BF-G-PNAG with BF-M03 so the same PNAG signal cannot count twice."
        ),
        not_established=(
            "A bile-acid-diarrhoea diagnosis, attribution of these community "
            "genes to R. gnavus without linkage, or that BssS is a "
            "one-direction 'more biofilm' switch. BssS is a regulator."
        ),
        notes=(
            "The abstract emphasises bssS, pgaA and pgaB; the results discuss all "
            "five. Both statements are stored."
        ),
    ),
    Module(
        id="BF-M05",
        label="Klebsiella type-1 and type-3 adhesion loci",
        measurement_kind="gene_family",
        mechanism_group="adhesion",
        primary_source_ids=("BF-S08", "BF-S09"),
        feature_symbols=("mrkA", "mrkD", "fimH", "ecpD"),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-ADHESION",
        interpretation=(
            "Persistence and adhesion potential with clinical-isolate evidence. "
            "Annotation anchors are MrkD WP_004149659.1, FimH BAH65076.1 and "
            "EcpD WP_002890060.1."
        ),
        not_established=(
            "Active gut biofilm, universal virulence, or a strain's phenotype: "
            "BF-S08 found substantial phenotype variation across 100 isolates "
            "with similar gene content. Those anchors are annotation anchors, "
            "not a sufficient diagnostic panel."
        ),
    ),
    Module(
        id="BF-M06",
        label="Extracellular polysaccharide and matrix-associated proteins",
        measurement_kind="gene_family",
        mechanism_group="matrix",
        primary_source_ids=("BF-S07",),
        feature_symbols=("eps", "bcsA", "bcsB", "wza", "wzb", "wzc"),
        census_state="context_only",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-MATRIX",
        interpretation=(
            "A broad discovery inventory. Capsule, cellulose and other matrix "
            "capabilities stay separately typed."
        ),
        not_established="That a capsule hit is an attached biofilm.",
    ),
    Module(
        id="BF-M07",
        label="eDNA, quorum sensing and cyclic-di-GMP regulatory context",
        measurement_kind="gene_family",
        mechanism_group="regulation",
        primary_source_ids=("BF-S07", "BF-S18"),
        feature_symbols=("luxS", "pfs", "sdiA", "ihfA", "ihfB", "hupA"),
        census_state="context_only",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-REGULATION",
        interpretation=(
            "Broadly conserved and pleiotropic genes, reported as annotation "
            "context only. BF-S18 showed the same regulator can act in opposite "
            "directions at different sites."
        ),
        not_established=(
            "Any harmful-score contribution from housekeeping counts. Total "
            "extracted DNA also cannot distinguish intracellular from matrix DNA."
        ),
    ),
    Module(
        id="BF-M08",
        label="C. difficile adhesion, matrix and spore context",
        measurement_kind="cross_reference",
        mechanism_group="organism_specific",
        primary_source_ids=("BF-S03",),
        feature_symbols=("cwp84", "slpA", "tcdA", "tcdB"),
        census_state="context_only",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-ORGANISM",
        interpretation=(
            "Reuses the existing validated organism and toxin assays. Biofilm "
            "potential, viable spores and toxin potential stay separate."
        ),
        not_established="Spore viability or a biofilm reservoir, neither of which stool DNA measures.",
    ),
    Module(
        id="BF-M09",
        label="Candida adhesion, matrix and regulatory context",
        measurement_kind="gene_family",
        mechanism_group="fungal",
        primary_source_ids=("BF-S03",),
        feature_symbols=("ALS1", "ALS3", "HWP1", "BCR1", "EFG1"),
        census_state="not_measurable_from_current_assay",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-FUNGAL",
        interpretation=(
            "Research fungal biofilm potential. Eukaryotic depth in stool shotgun "
            "is low and the ALS family is highly paralogous, so unresolved is the "
            "expected state rather than a failure."
        ),
        not_established=(
            "Candidiasis or invasive disease from presence. Fungal targets do not "
            "inherit bacterial detection limits."
        ),
    ),
    Module(
        id="BF-M10",
        label="L. reuteri GtfW protective mechanism",
        measurement_kind="gene_family",
        mechanism_group="protective_strain",
        primary_source_ids=("BF-S15",),
        feature_symbols=("gtfW",),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="protective_associated",
        aggregation_group="BF-G-PROTECTIVE-STRAIN",
        directional_claim=(
            "BF-S15 observed improved probiotic properties in the biofilm state "
            "of L. reuteri ATCC 23272, with GtfW implicated in the mechanism."
        ),
        molecular_resolution="strain",
        interpretation=(
            "A candidate protective mechanism conditioned on strain identity. "
            "Species-level abundance alone does not activate it."
        ),
        not_established=(
            "That the published microsphere formulation can be inferred from "
            "stool DNA, or that any L. reuteri carries the studied phenotype."
        ),
    ),
    Module(
        id="BF-M11",
        label="F. prausnitzii host-glycan protective mechanism",
        measurement_kind="gene_family",
        mechanism_group="protective_strain",
        primary_source_ids=("BF-S14",),
        feature_symbols=(),
        census_state="source_unavailable",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PROTECTIVE-STRAIN",
        interpretation=(
            "A protective evidence card and optional multimodal context. No "
            "validated specific DNA panel was established by BF-S14, so there "
            "are no feature symbols to measure."
        ),
        not_established=(
            "A protective score from generic F. prausnitzii or cobalamin "
            "abundance. Those are not substituted for the missing panel."
        ),
        notes=(
            "source_unavailable is the honest state: the mechanism is real but "
            "the sequence panel that would detect it does not exist yet."
        ),
    ),
    Module(
        id="BF-M12",
        label="Fibre-associated mixed-community mechanisms",
        measurement_kind="taxonomic_pattern",
        mechanism_group="protective_ecology",
        primary_source_ids=("BF-S16",),
        feature_symbols=(),
        census_state="context_only",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PROTECTIVE-ECOLOGY",
        interpretation=(
            "Supportive ecological evidence that fibre-associated biofilm "
            "formation can be beneficial."
        ),
        not_established=(
            "An exclusive DNA signature. Fibre-degradation genes are not proof "
            "of protective biofilm."
        ),
    ),
    Module(
        id="BF-M13",
        label="Other verified probiotic-strain adhesion and EPS mechanisms",
        measurement_kind="gene_family",
        mechanism_group="protective_strain",
        primary_source_ids=("BF-S26", "BF-S27", "BF-S28"),
        feature_symbols=(),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PROTECTIVE-STRAIN",
        interpretation=(
            "Expandable only with exact strain, gene and phenotype evidence plus "
            "a biological claim audit."
        ),
        not_established=(
            "Anything resembling a rule that Lactobacillus presence is good. "
            "That rule is explicitly forbidden."
        ),
    ),
    Module(
        id="BF-M14",
        label="S. aureus host fibrin-interaction capability",
        measurement_kind="gene_family",
        mechanism_group="fibrin",
        primary_source_ids=("BF-S17",),
        feature_symbols=("coa", "vwb", "clfA", "clfB"),
        census_state="context_only",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-FIBRIN",
        interpretation=(
            "An extraintestinal and context panel. Fibrin is a host substrate; "
            "these genes describe a capability to interact with it."
        ),
        not_established=(
            "Measured host fibrin, a cardiac or blood diagnosis, or an "
            "intestinal fibrin biofilm. This is never an automatic harmful-axis "
            "input."
        ),
    ),
    Module(
        id="BF-M15",
        label="Co-occurring pathogen, toxin and AMR calls",
        measurement_kind="cross_reference",
        mechanism_group="cross_reference",
        primary_source_ids=("BF-S20", "BF-S21"),
        feature_symbols=(),
        census_state="implemented",
        analytical_status="supported_at_stated_resolution",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-CROSSREF",
        interpretation=(
            "Cross-references the existing validated pathogen, toxin and AMR "
            "calls as independent clinical and mechanistic context."
        ),
        not_established=(
            "Same-cell linkage or a measured resistant biofilm. Co-occurrence in "
            "a community is not co-location in a cell."
        ),
    ),
    Module(
        id="BF-M16",
        label="Measured host thrombin, protein activity and tissue architecture",
        measurement_kind="external_assay",
        mechanism_group="host_assay",
        primary_source_ids=("BF-S12", "BF-S13"),
        feature_symbols=(),
        census_state="not_measurable_from_current_assay",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-HOST",
        interpretation="An optional externally measured assay layer.",
        not_established=(
            "Anything inferred from F2 DNA, taxon abundance or generic fliC "
            "expression. Host protein activity is measured or it is absent."
        ),
    ),
    Module(
        id="BF-M17",
        label="Longitudinal shared strain with biofilm-associated loci",
        measurement_kind="cross_reference",
        mechanism_group="longitudinal",
        primary_source_ids=("BF-S10",),
        feature_symbols=(),
        census_state="implemented",
        analytical_status="supported_at_stated_resolution",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-LONGITUDINAL",
        interpretation=(
            "Tracks potential and strain persistence using the existing strain "
            "resolution lane."
        ),
        not_established=(
            "That an attached biofilm or a disease was transferred. Shared "
            "strain detection is shared strain detection."
        ),
    ),
    Module(
        id="BF-M18",
        label="Two-feature biofilm-associated community proxy (H-C)",
        measurement_kind="taxonomic_pattern",
        mechanism_group="community_proxy",
        primary_source_ids=("BF-S01",),
        feature_symbols=("ecoli_complex", "m_gnavus"),
        census_state="implemented",
        analytical_status="supported_at_stated_resolution",
        candidate_axis="harmful_associated",
        aggregation_group="BF-G-COMMUNITY-PROXY",
        directional_claim=(
            "BF-S01 found E. coli / Escherichia-Shigella and R. gnavus relatively "
            "enriched in biofilm-positive participants. The direction is taken "
            "from that reported association and nothing else."
        ),
        molecular_resolution="species_complex",
        interpretation=(
            "A deliberately narrower stool-WGS translation of a tissue-derived "
            "16S association. Reported as an experimental proxy under its own "
            "heading, never as measured harmful biofilm."
        ),
        not_established=(
            "A stool biofilm probability. The transport from biopsy 16S "
            "Escherichia/Shigella to stool WGS E. coli species complex is "
            "recorded as a measurement change, not waved away."
        ),
        notes="Scored by BF-PROXY-COMMUNITY-1.0, never by the H-M mechanism axis.",
    ),
    Module(
        id="BF-M19",
        label="Four-genus biofilm-negative-associated ecological proxy (P-E)",
        measurement_kind="taxonomic_pattern",
        mechanism_group="ecology_proxy",
        primary_source_ids=("BF-S01",),
        feature_symbols=("faecalibacterium", "coprococcus", "subdoligranulum", "blautia"),
        census_state="implemented",
        analytical_status="supported_at_stated_resolution",
        candidate_axis="protective_associated",
        aggregation_group="BF-G-ECOLOGY-PROXY",
        directional_claim=(
            "BF-S01 reported these four genera relatively depleted in "
            "biofilm-positive biopsies. The direction is that reported "
            "depletion and nothing more."
        ),
        molecular_resolution="genus",
        interpretation=(
            "One correlated ecological group, not four independent protective "
            "mechanisms. Low values mean less of this narrow community pattern."
        ),
        not_established=(
            "A missing beneficial biofilm, a species deficiency, a need for "
            "transplant, or an instruction to supplement these organisms. Their "
            "stool-WGS use is a hypothesis, not validated transport."
        ),
        notes="Scored by BF-PROXY-ECOLOGY-1.0. A positive P-E never activates P-M.",
    ),
    Module(
        id="BF-M20",
        label="Linked matrix/adhesion plus independently supported damaging capability",
        measurement_kind="cross_reference",
        mechanism_group="linked_concern",
        primary_source_ids=("BF-S20",),
        feature_symbols=(),
        census_state="implemented",
        analytical_status="supported_at_stated_resolution",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-LINKED",
        interpretation=(
            "A concern-context overlay that joins a matrix or adhesion call to "
            "an independently supported damaging capability, but only when "
            "carrier linkage is direct_supported."
        ),
        not_established=(
            "Anything from co-occurrence on different unknown carriers. This is "
            "an overlay, not an extra copy of the same feature in the axis."
        ),
        notes=(
            "BF-S20 is the anchor because it found colibactin-marker carriage "
            "associated with dysplasia while generic biofilm presence was not."
        ),
    ),
    Module(
        id="BF-M21",
        label="Faeces-derived community behaviour",
        measurement_kind="external_assay",
        mechanism_group="behaviour",
        primary_source_ids=("BF-S19", "BF-S23"),
        feature_symbols=(),
        census_state="not_measurable_from_current_assay",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-BEHAVIOUR",
        interpretation=(
            "Optional externally measured biomass, dispersal, viable burden, "
            "epithelial effect and intervention response. Lab-condition scope is "
            "preserved: this is a faeces-derived community ex vivo."
        ),
        not_established=(
            "Any of these values reconstructed from stool DNA. BF-S19 found "
            "greater adverse behaviour despite lower biomass, so biomass is not "
            "a stand-in for the rest of the vector."
        ),
    ),
    Module(
        id="BF-M22",
        label="Tissue architecture and host-interface observations",
        measurement_kind="external_assay",
        mechanism_group="architecture",
        primary_source_ids=("BF-S01", "BF-S32", "BF-S04"),
        feature_symbols=(),
        census_state="not_measurable_from_current_assay",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-ARCHITECTURE",
        interpretation=(
            "Optional image or pathology-derived coverage, thickness, epithelial "
            "proximity and mucus response, per sampled site."
        ),
        not_established=(
            "A universal whole-gut amount from per-site sampling. A negative site "
            "is a negative site."
        ),
    ),
    Module(
        id="BF-M23",
        label="Matched expression and activity support",
        measurement_kind="external_assay",
        mechanism_group="expression",
        primary_source_ids=("BF-S22", "BF-S13"),
        feature_symbols=(),
        census_state="not_measurable_from_current_assay",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-EXPRESSION",
        interpretation=(
            "Optional RNA or protein support for a named DNA mechanism, requiring "
            "matched sampling and context before any joint claim."
        ),
        not_established=(
            "Physical attachment from generic expression. BF-S22's two arms were "
            "both biofilms, so it cannot supply a biofilm-versus-planktonic model."
        ),
    ),
    Module(
        id="BF-M24",
        label="Additional source-specific protective-antagonism evidence",
        measurement_kind="cross_reference",
        mechanism_group="protective_strain",
        primary_source_ids=("BF-S25", "BF-S29", "BF-S30", "BF-S31"),
        feature_symbols=(),
        census_state="candidate_not_validated",
        analytical_status="pending_validation",
        candidate_axis="context_dependent",
        aggregation_group="BF-G-PROTECTIVE-STRAIN",
        interpretation=(
            "Expands the protective evidence registry from BF-S25 to BF-S31, "
            "keeping species or family context where exact strain resolution is "
            "not supported by public references."
        ),
        not_established=(
            "That a named probiotic species implies the tested strain or "
            "formulation."
        ),
    ),
)

MODULES: Final[dict[str, Module]] = {m.id: m for m in _M}

#: Fixed dependence groups for each axis. Section 8.3 requires these to be
#: frozen before any sample is viewed, and requires BF-M03 and BF-M04 to
#: share one group because they measure the same PNAG signal.
HARMFUL_AXIS_GROUPS: Final = (("BF-M18",),)
PROTECTIVE_AXIS_GROUPS: Final = (("BF-M19",),)


def get(module_id: str) -> Module:
    try:
        return MODULES[module_id]
    except KeyError:
        raise KeyError(
            f"unknown module {module_id!r}; known IDs are BF-M01..BF-M24"
        ) from None


def census() -> dict[str, object]:
    """Completeness census across all 24 modules.

    Printed in the report so an unimplemented module is visibly
    unimplemented rather than absent. The largest bucket is
    ``candidate_not_validated``, which is the correct state for a mechanism
    with real evidence and no validated sequence panel.
    """
    by_state: dict[str, list[str]] = {state: [] for state in CENSUS_STATES}
    for module in MODULES.values():
        by_state[module.census_state].append(module.id)
    eligible = [m.id for m in MODULES.values() if m.eligible_for_axis]
    groups: dict[str, list[str]] = {}
    for module in MODULES.values():
        groups.setdefault(module.aggregation_group, []).append(module.id)
    return {
        "n_modules": len(MODULES),
        "by_census_state": dict(by_state),
        "n_eligible_for_axis": len(eligible),
        "eligible_module_ids": eligible,
        "dependence_groups": groups,
        "shared_groups": {
            name: ids for name, ids in groups.items() if len(ids) > 1
        },
        "note": (
            "Modules sharing a dependence group measure one signal and cannot "
            "count twice. BF-M03 and BF-M04 share BF-G-PNAG for exactly that "
            "reason."
        ),
    }
