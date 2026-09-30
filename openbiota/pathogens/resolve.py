"""Species resolution, pathotype gating and carriage context.

This pass runs after the organism screen and the determinant screen, with
both results in hand, and answers three questions the alignment stage cannot:

1. **Can this species actually be named?** Some organisms are not separable by
   sequence identity from a close relative that lives in most healthy guts.
   *Shigella* is the hard case: the Genome Taxonomy Database reclassified every
   *Shigella* species as a later heterotypic synonym of *Escherichia coli*
   (R07-RS207), having found that its earlier compromise renamed ~60% of
   *E. coli* genomes — including K-12 — to "*E. flexneri*". 16S cannot separate
   them (>99% identity) and neither can whole-genome ANI. A species-level
   *Shigella* call therefore requires the invasion markers, not breadth; in a
   stool sample without them, ordinary commensal *E. coli* is by far the more
   likely source of the sequence. Unresolved evidence is preserved on the
   catalogue's own combined *Shigella*/EIEC group row, which exists for exactly
   this purpose.

2. **Is the disease-causing trait present?** For an organism whose pathogenicity
   *is* a toxin or pathotype — *C. difficile*, ETBF, STEC, EPEC, ETEC, EAEC —
   the species is not the finding. Non-toxigenic *C. difficile* is not a
   pathogen; it is a common coloniser and is associated with *protection* from
   *C. difficile* infection. So the organism row reports its DNA honestly and
   the pathogen classification waits for the determinant.

3. **Is this normal?** For organisms with published healthy-carriage rates, the
   row carries that rate and its source, because "is this finding normal?" is
   the first question a reader asks and a report that cannot answer it invites
   the reader to assume the worst.

Nothing here deletes evidence. A demoted row keeps every count, gains a reason
code, and says what would settle it.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any, Final

#: Determinant statuses that count as the trait being present. A homolog or an
#: ambiguous allele is not the gene: it is a reason to test.
MARKER_SUPPORTED: Final[frozenset[str]] = frozenset(
    {"supported_intact_sequence", "supported_partial_sequence"}
)


@dataclass(frozen=True)
class Complex:
    """A group of targets that sequence identity cannot separate."""

    complex_id: str
    label: str
    #: Catalogue `family` values whose species/pathotype rows are gated.
    families: tuple[str, ...]
    #: Where unresolved evidence is preserved.
    fallback_target_id: str
    #: The near neighbour that makes the group unresolvable, for the wording.
    confounder: str
    why: str
    source: str


ESCHERICHIA_SHIGELLA: Final = Complex(
    complex_id="escherichia_shigella",
    label="Escherichia coli / Shigella complex",
    families=("Shigella/EIEC", "Diarrheagenic E. coli", "Other Escherichia"),
    fallback_target_id="bacteria.shigella_eiec_combined_group",
    confounder="commensal Escherichia coli",
    why=(
        "Shigella species and the diarrhoeagenic Escherichia coli pathotypes are the same "
        "genomic species as ordinary gut E. coli, which most healthy adults carry. Alignment "
        "identity and genome breadth cannot separate them, so a species or pathotype name "
        "needs its defining virulence marker"
    ),
    source=(
        "GTDB R07-RS207 reclassified Shigella as heterotypic synonyms of E. coli "
        "(doi:10.1101/2021.09.22.461432); 16S identity between them exceeds 99% "
        "(PMC5711669)"
    ),
)

COMPLEXES: Final[tuple[Complex, ...]] = (ESCHERICHIA_SHIGELLA,)

#: The marker evidence that licenses a species or pathotype name inside a
#: complex. Any one of the listed determinants being supported is enough;
#: an empty tuple means no validated single marker exists in this panel, so
#: the name can never be claimed from sequence alone.
SPECIES_MARKERS: Final[Mapping[str, tuple[str, ...]]] = {
    "bacteria.shigella_boydii": ("pathotype_marker.ipah",),
    "bacteria.shigella_sonnei": ("pathotype_marker.ipah",),
    "bacteria.shigella_flexneri": ("pathotype_marker.ipah",),
    "bacteria.shigella_dysenteriae": ("pathotype_marker.ipah",),
    "bacteria.eiec": ("pathotype_marker.ipah", "pathotype_marker.virf"),
    "bacteria.shigella_eiec_combined_group": ("pathotype_marker.ipah",),
    "bacteria.stec": ("toxin.stx1", "toxin.stx2"),
    "bacteria.etec": ("toxin.elta", "toxin.eltb", "toxin.esta_sta"),
    "bacteria.epec": ("pathotype_marker.eae", "pathotype_marker.bfpa"),
    "bacteria.eaec": ("pathotype_marker.aggr",),
    "bacteria.daec": (),
    "bacteria.hybrid_diarrheagenic_escherichia_coli": (),
    "bacteria.escherichia_albertii": (),
    "bacteria.escherichia_fergusonii": (),
}

#: Determinant kinds that can define pathogenicity. A resistance gene cannot:
#: S. aureus plus an unlinked mecA is not MRSA, and the catalogue says so.
DISEASE_DEFINING_KINDS: Final[frozenset[str]] = frozenset(
    {"toxin", "virulence_locus", "pathotype_marker"}
)


def pathotype_markers_for(
    target_id: str, determinants: Sequence[Any]
) -> tuple[str, ...]:
    """The disease-defining determinants the catalogue links to one target.

    Driven by the determinant screen's own `associated_target_ids` rather than
    a hand-kept table, so a target can never quietly escape the gate: adding
    *Staphylococcus aureus* to the catalogue without listing its enterotoxins
    here is exactly how it came to be counted as a disease-causing organism on
    five DNA fragments with no toxin gene in sight.
    """
    out: list[str] = []
    for d in determinants:
        if str(getattr(d, "kind", "") or "") not in DISEASE_DEFINING_KINDS:
            continue
        if target_id in tuple(getattr(d, "associated_target_ids", ()) or ()):
            out.append(str(getattr(d, "determinant_id", "")))
    return tuple(sorted(out))


@dataclass(frozen=True)
class Carriage:
    """A published healthy-carriage rate, so a row can answer "is this normal?"."""

    statement: str
    source: str


#: Only organisms with a defensible published rate appear here. Absence from
#: this table means "we do not have a carriage figure", never "abnormal".
CARRIAGE: Final[Mapping[str, Carriage]] = {
    "bacteria.clostridioides_difficile": Carriage(
        "Carried without symptoms by 0–17.5% of healthy adults; the toxin-producing form by "
        "about 1.5% (95% CI 0.7–2.6%). Toxin-negative strains are associated with protection "
        "from C. difficile infection, not with disease. Infant carriage reaches 75%.",
        "Meta-analysis of 51 studies, 39,447 participants (doi:10.1186/s13099-024-00674-0); "
        "non-toxigenic protection (doi:10.3389/fmicb.2018.01700)",
    ),
    "bacteria.shigella_eiec_combined_group": Carriage(
        "Escherichia coli, which this group cannot be separated from by sequence, is a normal "
        "member of most adult guts. Shigella itself is a notifiable infection and is not "
        "carried asymptomatically at population scale.",
        "GTDB R07-RS207 (doi:10.1101/2021.09.22.461432)",
    ),
    "protozoa.blastocystis_spp": Carriage(
        "Found in about 16% of healthy people worldwide and 20–30% of European adults, where it "
        "tracks with healthier diets rather than disease.",
        "56,989 metagenomes across 32 countries (doi:10.1016/j.cell.2024.06.018); Flemish Gut "
        "Flora Project, 30% (doi:10.1136/gutjnl-2018-316106)",
    ),
    "fungi.candida_albicans": Carriage(
        "Detected in the stool of 82.9% of 695 healthy adults by quantitative PCR.",
        "Milieu Intérieur cohort (PMC10732203)",
    ),
    "bacteria.klebsiella_pneumoniae": Carriage(
        "Gut carriage in 16.3% of a general adult population (95% CI 15.0–17.7).",
        "Tromsø 7 population study, n=2,975 (PMC8244762)",
    ),
    "bacteria.klebsiella_pneumoniae_complex": Carriage(
        "Gut carriage in 16.3% of a general adult population (95% CI 15.0–17.7).",
        "Tromsø 7 population study, n=2,975 (PMC8244762)",
    ),
    "bacteria.enterococcus_faecalis": Carriage(
        "A lifelong member of the normal gut community in the great majority of healthy adults.",
        "OpenBiota pathogen catalogue, resident lane",
    ),
    "bacteria.enterococcus_faecium": Carriage(
        "A common member of the normal gut community in healthy adults.",
        "OpenBiota pathogen catalogue, resident lane",
    ),
}


def _supported_markers(determinants: Sequence[Any]) -> set[str]:
    out: set[str] = set()
    for d in determinants:
        status = str(getattr(d, "determinant_status", "") or "")
        if status in MARKER_SUPPORTED:
            out.add(str(getattr(d, "determinant_id", "")))
    return out


def _assessed_markers(determinants: Sequence[Any]) -> set[str]:
    out: set[str] = set()
    for d in determinants:
        status = str(getattr(d, "determinant_status", "") or "")
        if status != "not_assessed":
            out.add(str(getattr(d, "determinant_id", "")))
    return out


def complex_for(family: str) -> Complex | None:
    for cx in COMPLEXES:
        if family in cx.families:
            return cx
    return None


#: The one organism the strict mode applies to.
CDIFF_TARGET: Final = "bacteria.clostridioides_difficile"

#: The statement the strict mode attaches when the species is found and the
#: toxin genes are not. Written as a positive finding that calls for action,
#: which is what the mode is for: a reader who has asked for it wants a
#: species-level detection surfaced as a pathogen call, not filed as carriage.
STRICT_CDIFF_STATEMENT: Final = (
    "This test scanned the DNA fragments in your sample and found a definite positive "
    "match to multiple DNA sequences of <i>Clostridioides difficile</i>. The genes that "
    "produce its toxins (tcdA, tcdB) were not found among these fragments. That does not "
    "rule out <i>C. difficile</i> disease: toxin genes are a small target that a stool "
    "shotgun test can miss at ordinary depth, and finding the organism at all raises the "
    "possibility of disease relative to a sample in which it was not found. It cannot "
    "confirm it. This may warrant further action."
)


def resolve(
    results: Mapping[str, Any],
    determinants: Sequence[Any],
    *,
    families: Mapping[str, str],
    requires: Mapping[str, bool] | None = None,
    strict_cdiff: bool = False,
) -> dict[str, Any]:
    """Return the results with species resolution and pathotype gating applied.

    `families` maps target_id → catalogue family and `requires` maps
    target_id → `requires_determinants`; the caller reads both from the
    compiled targets.

    ``strict_cdiff`` changes one gate. Ordinarily a *C. difficile* species
    detection without tcdA/tcdB is filed as carriage, because the
    non-toxigenic organism is a common coloniser and the literature
    associates it with protection from infection, not disease. In strict
    mode the species detection itself is surfaced as a positive pathogen
    call, stated as such, with the toxin-gene result reported alongside as
    a limitation of the assay rather than as a clearance. The mode is
    recorded in the output so the reader knows which interpretation they
    are looking at.
    """
    requires = dict(requires or {})
    supported = _supported_markers(determinants)
    assessed = _assessed_markers(determinants)
    out: dict[str, Any] = dict(results)

    positive = {"supported_sequence", "marker_signal", "candidate_signal"}
    absorbed: dict[str, list[str]] = {}

    for target_id, row in results.items():
        family = str(families.get(target_id, "") or getattr(row, "family", "") or "")
        cx = complex_for(family)
        markers = SPECIES_MARKERS.get(target_id)
        status = str(getattr(row, "sequence_status", "") or "")

        # ---- 1. species resolution inside an indistinguishable complex ---- #
        # A row with no organism-specific fragment has nothing to reassign: the
        # pathogen lane already calls it ambiguous because every fragment is
        # shared. Demoting it would print "DNA present, species cannot be
        # named" over zero evidence, which is worse than saying nothing.
        specific = int(getattr(row, "unique_supporting_fragments", 0) or 0)
        if cx is not None and target_id != cx.fallback_target_id and markers is not None:
            has_marker = bool(set(markers) & supported)
            if not has_marker and status in positive and specific > 0:
                why_no_marker = (
                    "no validated single marker distinguishes it in this panel"
                    if not markers
                    else (
                        f"{', '.join(sorted(markers))} not detected"
                        if set(markers) & assessed
                        else f"{', '.join(sorted(markers))} not assessed"
                    )
                )
                out[target_id] = _demote(
                    row,
                    complex_=cx,
                    why_no_marker=why_no_marker,
                    markers=markers,
                )
                absorbed.setdefault(cx.fallback_target_id, []).append(
                    str(getattr(row, "display_name", target_id))
                )
                continue
            if has_marker:
                out[target_id] = _annotate(
                    row,
                    species_resolution="resolved_by_marker",
                    pathotype_evidence="supported",
                    pathotype_markers=tuple(sorted(set(markers) & supported)),
                    counts_as_pathogen=status in {"supported_sequence", "marker_signal"},
                    complex_id=cx.complex_id,
                    complex_label=cx.label,
                )
                continue

        # ---- 2. pathotype-dependent organisms ----------------------------- #
        # Gated when the catalogue says the disease is the trait, not the
        # species (`requires_determinants`), using the determinants actually
        # linked to this target.
        needed = (
            pathotype_markers_for(target_id, determinants)
            if requires.get(target_id)
            else ()
        )
        # As with species resolution: a row with no organism-specific fragment
        # has nothing to gate. Calling it "carriage without the toxin genes"
        # asserts the organism is present, which zero specific fragments does
        # not support.
        if needed and status in positive and specific > 0:
            present = sorted(set(needed) & supported)
            if present:
                out[target_id] = _annotate(
                    row,
                    pathotype_evidence="supported",
                    pathotype_markers=tuple(present),
                    counts_as_pathogen=status in {"supported_sequence", "marker_signal"},
                )
            elif strict_cdiff and target_id == CDIFF_TARGET:
                # Strict mode: the species is the finding. The toxin result
                # stays on the row, described as a limit of what shotgun
                # sequencing can see rather than as evidence of safety.
                looked = bool(set(needed) & assessed)
                out[target_id] = _annotate(
                    row,
                    pathotype_evidence="not_detected" if looked else "not_assessed",
                    pathotype_markers=(),
                    counts_as_pathogen=True,
                    report_tier="pathogen",
                    strict_mode="strict_cdiff",
                    extra_reason="strict_cdiff_species_positive",
                    statement_override=STRICT_CDIFF_STATEMENT,
                )
            else:
                looked = bool(set(needed) & assessed)
                out[target_id] = _annotate(
                    row,
                    pathotype_evidence="not_detected" if looked else "not_assessed",
                    pathotype_markers=(),
                    counts_as_pathogen=False,
                    report_tier="pathotype_negative",
                    extra_reason=(
                        "defining_determinant_not_detected" if looked
                        else "defining_determinant_not_assessed"
                    ),
                    statement_suffix=(
                        f"The genes that make this organism capable of causing disease "
                        f"({', '.join(_gene_names(needed))}) were "
                        + ("looked for and not found" if looked else "not assessed")
                        + " in this sample, so this is carriage evidence rather than a "
                        "disease finding."
                    ),
                )
            continue

        # ---- 3. everything else keeps its own classification -------------- #
        out[target_id] = _annotate(row)

    # ---- the combined-group row inherits what could not be resolved ------- #
    for fallback_id, names in absorbed.items():
        row = out.get(fallback_id)
        if row is None:
            continue
        cx = next((c for c in COMPLEXES if c.fallback_target_id == fallback_id), None)
        if cx is None:
            continue
        markers = SPECIES_MARKERS.get(fallback_id) or ()
        has_marker = bool(set(markers) & supported)
        # The group's own reference is the confounder's genome, so it carries no
        # distinguishing k-mer of its own and its status is usually
        # `not_assessed`. The alignment evidence is nonetheless real: it simply
        # cannot be resolved to a species. Aggregate the strongest member
        # evidence here so the reader sees one honest row instead of a row that
        # says "not assessed" under a paragraph describing what was found.
        members = [
            results[tid] for tid in results
            if str(getattr(results[tid], "display_name", "")) in set(names)
            and int(getattr(results[tid], "unique_supporting_fragments", 0) or 0) > 0
        ]
        aggregate: dict[str, Any] = {}
        if members and str(getattr(row, "sequence_status", "")) in (
            "not_assessed", "not_detected"
        ):
            aggregate = {
                "sequence_status": "ambiguous_signal",
                "display_status": "ambiguous_signal",
                "display_qualifier": (
                    "evidence aggregated from the complex; species not resolvable"
                ),
                "unique_supporting_fragments": max(
                    int(getattr(m, "unique_supporting_fragments", 0) or 0) for m in members
                ),
                "informative_regions_supported": max(
                    int(getattr(m, "informative_regions_supported", 0) or 0) for m in members
                ),
                "informative_bases_covered": max(
                    int(getattr(m, "informative_bases_covered", 0) or 0) for m in members
                ),
                "reference_breadth_fraction": max(
                    float(getattr(m, "reference_breadth_fraction", 0.0) or 0.0) for m in members
                ),
            }
        out[fallback_id] = _annotate(
            row,
            **aggregate,
            species_resolution="group_level_only",
            complex_id=cx.complex_id,
            complex_label=cx.label,
            absorbed_from=tuple(sorted(set(names))),
            pathotype_evidence="supported" if has_marker else (
                "not_detected" if set(markers) & assessed else "not_assessed"
            ),
            pathotype_markers=tuple(sorted(set(markers) & supported)),
            counts_as_pathogen=has_marker,
            report_tier=None if has_marker else "unresolved_complex",
            statement_override=(
                (
                    f"Species-level names could not be resolved within the {cx.label}; the "
                    f"evidence reported for {', '.join(sorted(set(names)))} is carried here "
                    "instead. "
                )
                + (
                    "An invasion marker was found, so this group is reported as a finding."
                    if has_marker
                    else (
                        "No invasion marker was found. In stool, sequence like this most often "
                        f"comes from {cx.confounder}, which the alignment cannot tell apart. "
                        "A clinical stool culture or NAAT is what distinguishes them."
                    )
                )
            ),
        )
    return out


def _gene_names(determinant_ids: Sequence[str]) -> list[str]:
    """`toxin.tcda` → `tcdA`-ish, for a sentence a reader can follow."""
    out = []
    for did in determinant_ids:
        name = did.split(".", 1)[-1]
        if name.startswith("tcd") or name.startswith("cdt") or name.startswith("stx"):
            name = name[:3] + name[3:].upper()
        elif name.startswith("bft"):
            name = "bft-" + name.split("_")[-1]
        elif name == "esta_sta":
            name = "estA/sta"
        out.append(name)
    return out


def _annotate(row: Any, **fields: Any) -> Any:
    """Attach resolution metadata to a frozen result row."""
    extra_reason = fields.pop("extra_reason", None)
    suffix = fields.pop("statement_suffix", None)
    override = fields.pop("statement_override", None)
    updates: dict[str, Any] = {}
    for key, value in fields.items():
        updates[key] = value
    if extra_reason:
        updates["reason_codes"] = tuple(
            dict.fromkeys([*getattr(row, "reason_codes", ()), extra_reason])
        )
    if suffix:
        base = str(getattr(row, "plain_statement", "") or "").rstrip()
        updates["plain_statement"] = f"{base} {suffix}".strip()
    if override:
        updates["plain_statement"] = str(override).strip()
    carriage = CARRIAGE.get(str(getattr(row, "target_id", "")))
    if carriage is not None:
        updates["carriage_statement"] = carriage.statement
        updates["carriage_source"] = carriage.source
    if not updates:
        return row
    return replace(row, **updates)


def _demote(row: Any, *, complex_: Complex, why_no_marker: str, markers: Sequence[str]) -> Any:
    """Turn a species call the data cannot support into an honest group signal."""
    name = str(getattr(row, "display_name", "") or getattr(row, "target_id", ""))
    statement = (
        f"DNA consistent with {name} was aligned, but this species cannot be named from "
        f"sequence alone: {complex_.why}. Here, {why_no_marker}. The evidence is preserved on "
        f"the {complex_.label} row. {complex_.source}."
    )
    return replace(
        row,
        sequence_status="ambiguous_signal",
        display_status="ambiguous_signal",
        display_qualifier=f"species not resolvable within the {complex_.label}",
        plain_statement=statement,
        reason_codes=tuple(
            dict.fromkeys(
                [*getattr(row, "reason_codes", ()), "species_not_resolvable_within_complex"]
            )
        ),
        species_resolution="not_resolvable_within_complex",
        complex_id=complex_.complex_id,
        complex_label=complex_.label,
        pathotype_evidence="not_detected" if markers else "not_applicable",
        pathotype_markers=(),
        counts_as_pathogen=False,
        report_tier="unresolved_complex",
        confirmation_options=tuple(
            dict.fromkeys(
                [*getattr(row, "confirmation_options", ()), "clinical_stool_culture",
                 "clinical_stool_NAAT"]
            )
        ),
    )


def self_test() -> int:  # pragma: no cover - exercised by tests/test_pathogens.py
    """Contract checks for the two gates, with no I/O."""
    from openbiota.pathogens.schema import DeterminantResult, PathogenResult

    def organism(target_id: str, *, name: str, status: str = "supported_sequence") -> PathogenResult:
        return PathogenResult(
            schema_version="pathogens.result.v1",
            sample_id="S",
            target_id=target_id,
            group="bacteria",
            display_name=name,
            interpretation_class="established_enteric",
            assay_eligibility="eligible",
            analysis_status="completed",
            reference_status="genome_supported",
            sequence_status=status,
            contamination_status="controls_unavailable",
            resolution="species",
            clinical_interpretation="organism_sequence_not_infection_diagnosis",
            validation_scope="computational_research_rule",
            calling_profile_id="research_dna_v1",
            reference_bundle_id="b",
            unique_supporting_fragments=564,
            informative_regions_supported=20,
            informative_bases_covered=5000,
            display_status=status,
            plain_statement="base.",
        )

    def marker(did: str, status: str, target: str = "") -> DeterminantResult:
        return DeterminantResult(
            schema_version="5.0",
            sample_id="S",
            determinant_id=did,
            display_name=did,
            gene_family=did.split(".")[-1],
            kind=did.split(".")[0],
            associated_target_ids=(target,) if target else (),
            determinant_status=status,
            host_linkage="unlinked",
            assay_eligibility="eligible",
            analysis_status="completed",
            reference_status="genome_supported",
        )

    families = {
        "bacteria.shigella_sonnei": "Shigella/EIEC",
        "bacteria.shigella_eiec_combined_group": "Shigella/EIEC",
        "bacteria.clostridioides_difficile": "Toxin-dependent anaerobes",
    }
    rows = {
        "bacteria.shigella_sonnei": organism("bacteria.shigella_sonnei", name="Shigella sonnei"),
        "bacteria.shigella_eiec_combined_group": organism(
            "bacteria.shigella_eiec_combined_group", name="Shigella/EIEC (combined group)",
            status="candidate_signal",
        ),
        "bacteria.clostridioides_difficile": organism(
            "bacteria.clostridioides_difficile", name="Clostridioides difficile"
        ),
    }

    # No markers: the species call is demoted, C. difficile is carriage.
    requires = {"bacteria.clostridioides_difficile": True}
    cd_id = "bacteria.clostridioides_difficile"
    out = resolve(
        rows,
        [marker("pathotype_marker.ipah", "not_detected", "bacteria.shigella_sonnei"),
         marker("toxin.tcda", "not_detected", cd_id),
         marker("toxin.tcdb", "not_detected", cd_id)],
        families=families, requires=requires,
    )
    shig = out["bacteria.shigella_sonnei"]
    assert shig.sequence_status == "ambiguous_signal", shig.sequence_status
    assert shig.counts_as_pathogen is False
    assert shig.species_resolution == "not_resolvable_within_complex"
    assert shig.unique_supporting_fragments == 564, "evidence must be preserved"
    group = out["bacteria.shigella_eiec_combined_group"]
    assert group.absorbed_from == ("Shigella sonnei",)
    assert group.counts_as_pathogen is False
    assert "commensal Escherichia coli" in group.plain_statement
    cd = out["bacteria.clostridioides_difficile"]
    assert cd.pathotype_evidence == "not_detected"
    assert cd.counts_as_pathogen is False
    assert cd.report_tier == "pathotype_negative"
    assert cd.carriage_statement and "1.5%" in cd.carriage_statement
    assert "tcdA" in cd.plain_statement or "tcdB" in cd.plain_statement

    # With the invasion marker, the species call stands and counts.
    out2 = resolve(
        rows,
        [marker("pathotype_marker.ipah", "supported_intact_sequence", "bacteria.shigella_sonnei")],
        families=families, requires=requires,
    )
    assert out2["bacteria.shigella_sonnei"].counts_as_pathogen is True
    assert out2["bacteria.shigella_sonnei"].species_resolution == "resolved_by_marker"

    # With toxin genes, C. difficile is a pathogen finding again.
    out3 = resolve(
        rows,
        [marker("toxin.tcdb", "supported_intact_sequence", cd_id)],
        families=families, requires=requires,
    )
    assert out3["bacteria.clostridioides_difficile"].counts_as_pathogen is True
    assert out3["bacteria.clostridioides_difficile"].pathotype_evidence == "supported"

    # A resistance gene is not a disease-defining determinant: S. aureus plus
    # an unlinked mecA is not a finding.
    assert pathotype_markers_for(
        "bacteria.staphylococcus_aureus",
        [marker("amr.meca", "supported_intact_sequence", "bacteria.staphylococcus_aureus")],
    ) == ()

    # ---- strict C. difficile rule ---------------------------------------- #
    # Default: species without tcdA/tcdB is carriage (checked above). Strict:
    # the species is the finding, stated as a positive, toxin result kept.
    cd_dets = [marker("toxin.tcda", "not_detected", cd_id),
               marker("toxin.tcdb", "not_detected", cd_id)]
    strict = resolve(rows, cd_dets, families=families, requires=requires,
                     strict_cdiff=True)[cd_id]
    assert strict.report_tier == "pathogen" and strict.counts_as_pathogen is True
    assert strict.strict_mode == "strict_cdiff"
    assert "definite positive match" in strict.plain_statement
    assert "does not rule out" in strict.plain_statement
    assert "may warrant further action" in strict.plain_statement
    assert strict.pathotype_evidence == "not_detected", "the toxin result stays on the row"
    # It promotes a detection; it does not invent one. A row with no
    # organism-specific fragments is left untouched by the gate - the report
    # keeps it out of the headline on sequence status, which is a separate
    # filter from this one.
    absent = replace(rows[cd_id], sequence_status="not_detected", display_status="not_detected",
                     unique_supporting_fragments=0)
    out = resolve({cd_id: absent}, cd_dets, families=families, requires=requires,
                  strict_cdiff=True)[cd_id]
    assert out.strict_mode is None, "the strict rule must not fire without a detection"
    assert out.report_tier != "pathogen"
    # And it leaves every other toxin-gated organism alone.
    stec = "bacteria.stec"
    out = resolve({stec: organism(stec, name="STEC")},
                  [marker("toxin.stx1", "not_detected", stec), marker("toxin.stx2", "not_detected", stec)],
                  families={stec: "x"}, requires={stec: True}, strict_cdiff=True)[stec]
    assert out.report_tier == "pathotype_negative" and out.counts_as_pathogen is False
    return 0
